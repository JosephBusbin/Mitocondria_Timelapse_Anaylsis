#!/usr/bin/env python3
"""
collect_4d_sld_timelapses.py
Phase 1 & 2 Utility: SlideBook (.sld / .sldy) 4D Time-Lapse Discovery, Health Audit, Compression & Staging

Overview:
  Scans source directories of 3i SlideBook (.sld / .sldy) microscopy containers in strictly
  read-only mode, inspects internal series metadata via Bio-Formats ('showinf -nopix -omexml-only')
  without loading pixel buffers into memory, and filters for true 4D (X, Y, Z, T) time-lapse datasets.
  Every candidate file (readable or not) receives one row in 'sld_health_audit.csv' so multi-year
  lab archives can be triaged without re-running the scan.

Supported Operational Modes (--mode):
  1. 'scan-only'       : Generates 'sld_4d_inventory.csv', 'sld_health_audit.csv', 'movie_list.txt',
                         and 'extraction_provenance.json' reporting every series and calculating storage savings.
  2. 'extract-ometiff' : Extracts ONLY the qualifying 4D series into losslessly compressed (LZW or Deflate)
                         16-bit '.ome.tif' hyperstacks ready for Fast4DReg, ilastik, Cellpose, TrackMate 7,
                         and MitoLandscape. Automatically uses FIJI's bundled Java runtime if present.
  3. 'copy-sld'        : Copies (or hardlinks) only the raw '.sld' files (and companion '_Files' folders)
                         that contain at least one qualifying 4D series to a staging folder.

Zero-Dependency Design:
  Uses ONLY Python 3 standard library modules (argparse, concurrent.futures, csv, dataclasses,
  datetime, hashlib, json, os, pathlib, platform, re, shutil, subprocess, sys, xml.etree.ElementTree)
  so it runs out-of-the-box on a lab microscopy Windows PC or on W&M SciClone HPC via
  'module load python/3.12.7'.

SLURM Manifest Paths (movie_list.txt):
  Entries are written relative to the raw OME-TIFF root (e.g. '2025/<movie>.ome.tif') unless
  --hpc-raw-root or the MITO_RAW_ROOT environment variable supplies an absolute prefix. Relative
  entries are resolved against $MITO_RAW_ROOT at job run time, so the same manifest works no
  matter which storage location (personal or shared lab) the data lands in.

Sharded runs (SLURM job arrays, see scripts/submit_full_run.bash):
  --shard K/N      processes only every N-th candidate file (sorted list, round-robin, K = 1..N), so
                   N array tasks together cover every file exactly once.
  --file-list F    processes only the files named in F (one path per line, absolute or relative to
                   source_dir) instead of searching source_dir recursively.
  --report-dir D   where the CSV/JSON reports go (default: output_dir).
  With --shard, every report gets a '_shardK' suffix (e.g. sld_4d_inventory_shard3.csv) so
  concurrent tasks never overwrite each other; scripts/merge_collect_shards.py combines them.
  With --shard or --file-list, 'extract-ometiff' exits with code 4 if any extraction failed, so
  Slurm marks the task FAILED and it can be retried.
"""

import argparse
import concurrent.futures
import csv
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Optional, Tuple, Dict


@dataclass
class SeriesMetadata:
    """Stores dimensional, calibration, and provenance metadata for a single series inside an .sld file."""
    cohort_year: str
    sld_path: str
    sld_filename: str
    sld_file_size_gb: float
    series_index: int
    series_name: str
    size_x: int
    size_y: int
    size_z: int
    size_c: int
    size_t: int
    pixel_type: str
    bits_per_pixel: int
    phys_size_x_um: Optional[float]
    phys_size_y_um: Optional[float]
    phys_size_z_um: Optional[float]
    time_interval_s: Optional[float]
    estimated_uncompressed_gb: float
    is_4d_timelapse: bool
    calibration_complete: bool
    extracted_ometiff_relpath: str
    # Mito-construct filename classification (construct_classifier.py; empty when --mito-filter off)
    mito_label: str = ""
    mito_fp: str = ""
    mito_color: str = ""
    red_mito_flag: bool = False


@dataclass
class FileHealthRecord:
    """Audit record written to sld_health_audit.csv for every candidate .sld/.sldy file (pass or fail)."""
    cohort_year: str
    rel_path: str
    sld_filename: str
    file_size_mb: float
    modified_utc: str
    has_companion_dirs: bool
    magic_hex: str
    is_null_header: bool
    status: str
    total_series: int
    qualifying_4d_series: int
    error_summary: str


class ShowinfError(RuntimeError):
    """
    Raised when 'showinf' cannot parse a container. Carries the raw Java stderr separately from the
    human-readable hint so failures can be classified on the original exception text.
    """

    def __init__(self, hint: str, raw_stderr: str):
        super().__init__(hint)
        self.raw_stderr = raw_stderr


# OME unit conversion to the pipeline's canonical um / s (mirrors tiff_io._LENGTH_TO_UM/_TIME_TO_S;
# duplicated here because this script must stay standard-library only for the lab PC).
LENGTH_TO_UM = {
    "um": 1.0, "µm": 1.0, "μm": 1.0, "micron": 1.0, "microns": 1.0,
    "nm": 1e-3, "Å": 1e-4, "A": 1e-4, "mm": 1e3, "cm": 1e4, "m": 1e6, "pm": 1e-6,
}
TIME_TO_S = {"s": 1.0, "sec": 1.0, "ms": 1e-3, "µs": 1e-6, "us": 1e-6, "min": 60.0, "h": 3600.0}


def convert_ome_unit(value: Optional[str], unit: Optional[str], table: Dict[str, float]) -> Optional[float]:
    """Converts an OME numeric attribute to um or s. Missing unit = OME default (um / s); unknown unit = None."""
    if value is None:
        return None
    try:
        val = float(value)
    except ValueError:
        return None
    if not unit:
        return val
    factor = table.get(unit.strip())
    return val * factor if factor is not None else None


# Status values that mean showinf parsed the file successfully.
HEALTHY_STATUSES = {"OK_4D", "OK_NON_4D"}

# 3i's native SBReadFile library appends a timing trace to the hard-coded RELATIVE name
# "C:\\Java_SBReadFileJNI_openFile.txt" (a literal filename on Linux, megabytes per run) in the
# current directory. Bio-Formats subprocesses run from this throwaway directory so the trace
# never lands in the repo or the user's data folders. Set by main(); None = inherit cwd.
NATIVE_LOG_SINK_DIR: Optional[str] = None


def sanitize_filename(name: str) -> str:
    """Converts a SlideBook series name into a clean, POSIX-safe filename stem."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return re.sub(r"_+", "_", cleaned).strip("_") or "series"


def infer_cohort_year(sld_path: Path) -> str:
    """
    Infers the experimental cohort year (e.g., '2023', '2024') from parent folder names
    or filename date prefixes (e.g., '20240514_Tadpole.sld'). Defaults to 'unclassified'.
    """
    for part in reversed(sld_path.parts):
        match = re.search(r"\b(201\d|202\d)\b", part)
        if match:
            return match.group(1)
    match_stem = re.search(r"(201\d|202\d)", sld_path.stem)
    if match_stem:
        return match_stem.group(1)
    return "unclassified"


def estimate_series_gb(x: int, y: int, z: int, c: int, t: int, bits: int) -> float:
    """Estimates uncompressed byte size of a single 5D image series in gigabytes (GB)."""
    bytes_per_pixel = max(1, (bits + 7) // 8)
    total_bytes = int(x) * int(y) * int(z) * int(c) * int(t) * bytes_per_pixel
    return round(total_bytes / (1024 ** 3), 4)


def compute_sha256(file_path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Computes SHA-256 checksum of an extracted file for transfer integrity auditing."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def inspect_file_header(sld_path: Path) -> Tuple[float, str, str, bool]:
    """
    Fast read-only byte probe executed before invoking Java/showinf.
    Returns (file_size_mb, modified_utc_iso, first_16_bytes_hex, is_empty_or_null_header).
    A null header (first 4 KB all 0x00) indicates a pre-allocated file from a crashed acquisition.
    """
    try:
        st = sld_path.stat()
        size_mb = round(st.st_size / (1024 ** 2), 3)
        mtime = datetime.datetime.fromtimestamp(st.st_mtime, tz=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if st.st_size == 0:
            return 0.0, mtime, "", True
        with open(sld_path, "rb") as f:
            head = f.read(4096)
        magic_hex = head[:16].hex()
        is_null = (len(head) == 0) or all(b == 0 for b in head)
        return size_mb, mtime, magic_hex, is_null
    except OSError:
        return 0.0, "", "", True


def classify_showinf_error(stderr_text: str) -> Tuple[str, str]:
    """
    Categorizes raw showinf/Java stderr into actionable buckets for sld_health_audit.csv.
    Must be given the ORIGINAL Java stderr (ShowinfError.raw_stderr), not a rewritten hint.
    """
    cleaned = " | ".join(line.strip() for line in stderr_text.splitlines() if line.strip())
    short_err = cleaned[:300]
    low = cleaned.lower()
    if "not a valid slidebook document" in low or "problem reading the image records" in low:
        return "SLIDEBOOK_IMAGE_RECORDS_ERR", short_err
    # The open-source SlidebookReader (no '6') crashes on SlideBook 6+ block tables.
    if "arrayindexoutofboundsexception" in low or "slidebookreader" in low:
        return "MISSING_SB6_READER", short_err
    if "unsatisfiedlinkerror" in low or "noclassdeffounderror" in low or "outofmemoryerror" in low:
        return "JNI_OR_JAVA_ERR", short_err
    if "file already open" in low or "file not open" in low:
        return "SLIDEBOOK_SESSION_CONCURRENCY_ERR", short_err
    return "OTHER_ERR", short_err


def _tools_dirs() -> List[Path]:
    """
    Returns user-space pipeline tool roots: $MITO_TOOLS_DIR (HPC/custom installs) and the
    Windows lab default %LOCALAPPDATA%\\mito_pipeline_tools created by setup_lab_windows.ps1.
    """
    roots: List[Path] = []
    if os.environ.get("MITO_TOOLS_DIR"):
        roots.append(Path(os.environ["MITO_TOOLS_DIR"]))
    if os.environ.get("LOCALAPPDATA"):
        roots.append(Path(os.environ["LOCALAPPDATA"]) / "mito_pipeline_tools")
    return roots


def discover_fiji_java_home(explicit_fiji_dir: Optional[str] = None) -> Optional[str]:
    """
    Locates the OpenJDK Java runtime bundled inside FIJI (Fiji.app/java/...) on Windows,
    macOS, or Linux so the user does not need a separate system-wide Java installation.
    """
    if os.environ.get("BF_JAVA_HOME") and Path(os.environ["BF_JAVA_HOME"]).exists():
        return os.environ["BF_JAVA_HOME"]
    if os.environ.get("JAVA_HOME") and Path(os.environ["JAVA_HOME"]).exists():
        return os.environ["JAVA_HOME"]

    candidate_roots: List[Path] = []
    if explicit_fiji_dir:
        candidate_roots.append(Path(explicit_fiji_dir))

    home = Path.home()
    candidate_roots.extend([
        Path("C:/Fiji.app"),
        home / "Fiji.app",
        home / "Desktop" / "Fiji.app",
        home / "Downloads" / "Fiji.app",
        Path("/opt/Fiji.app"),
        Path("/Applications/Fiji.app"),
    ])

    for tools_root in _tools_dirs():
        candidate_roots.append(tools_root / "jre")

    is_win = sys.platform.startswith("win")
    java_exe_name = "java.exe" if is_win else "java"

    for root in candidate_roots:
        if not root.is_dir():
            continue
        direct_bin = root / "bin" / java_exe_name
        if direct_bin.is_file():
            return str(root.resolve())
        for java_bin in root.rglob(java_exe_name):
            if java_bin.is_file() and java_bin.parent.name == "bin":
                return str(java_bin.parent.parent.resolve())
    return None


def discover_slidebook_reader_jar(
    explicit_bftools_dir: Optional[str] = None,
    showinf_bin: Optional[str] = None,
    explicit_fiji_dir: Optional[str] = None,
) -> Optional[str]:
    """
    Locates 3i's official SlideBook6Reader.jar plugin containing native libraries
    (SBReadFile.dll on Windows, libSBReadFile.so on Linux) required to read modern
    SlideBook 6 & 7 (.sld / .sldy) files without open-source Bio-Formats parser errors.
    """
    # 1. Respect explicit or already set BF_CP environment variable
    if os.environ.get("BF_CP"):
        for entry in os.environ["BF_CP"].split(os.pathsep):
            p = Path(entry.strip('"'))
            if p.is_file() and "SlideBook" in p.name:
                return str(p.resolve())

    candidate_dirs: List[Path] = []
    if explicit_bftools_dir:
        candidate_dirs.append(Path(explicit_bftools_dir))
        # Sibling 'readers' folder: the jar lives outside bftools so version upgrades cannot delete it.
        candidate_dirs.append(Path(explicit_bftools_dir).parent / "readers")
    if showinf_bin and Path(showinf_bin).exists():
        candidate_dirs.append(Path(showinf_bin).resolve().parent)
        candidate_dirs.append(Path(showinf_bin).resolve().parent.parent / "readers")

    # User-space tool installs: <tools>/readers (current layout) and <tools>/bftools (legacy layout)
    for tools_root in _tools_dirs():
        candidate_dirs.append(tools_root / "readers")
        candidate_dirs.append(tools_root / "bftools")

    # Current working directory and scripts directory
    cwd = Path.cwd()
    candidate_dirs.extend([
        cwd,
        cwd / "bftools",
        Path(__file__).resolve().parent / "bftools",
        Path(__file__).resolve().parent / "scripts",
    ])

    # Search common FIJI paths for SlideBook plugin
    candidate_fiji: List[Path] = []
    if explicit_fiji_dir:
        candidate_fiji.append(Path(explicit_fiji_dir))
    home = Path.home()
    candidate_fiji.extend([
        Path("C:/Fiji.app"),
        home / "Fiji.app",
        home / "Desktop" / "Fiji.app",
        home / "Downloads" / "Fiji.app",
        Path("/opt/Fiji.app"),
        Path("/Applications/Fiji.app"),
    ])
    for froot in candidate_fiji:
        candidate_dirs.append(froot / "jars" / "bio-formats")
        candidate_dirs.append(froot / "plugins")

    for cdir in candidate_dirs:
        if cdir.is_dir():
            for jar in cdir.glob("*SlideBook*.jar"):
                if jar.is_file():
                    return str(jar.resolve())
    return None


def parse_ome_xml_metadata(
    xml_string: str,
    sld_path: Path,
    source_root: Path,
    min_z: int,
    min_t: int,
    min_c: int,
) -> List[SeriesMetadata]:
    """
    Parses OME-XML metadata returned by Bio-Formats 'showinf -nopix -omexml-only' to extract
    per-series dimensions (X, Y, Z, C, T) and physical calibrations (dx, dy, dz, dt).
    """
    results: List[SeriesMetadata] = []
    file_size_gb = round(sld_path.stat().st_size / (1024 ** 3), 4) if sld_path.exists() else 0.0
    cohort_year = infer_cohort_year(sld_path)

    try:
        rel_parent = sld_path.parent.relative_to(source_root)
    except ValueError:
        rel_parent = Path(cohort_year)
    if rel_parent == Path("."):
        # Files at the top of source_root get a cohort-year folder so every staged movie has a
        # cohort directory (the SLURM job names scratch/output folders after it).
        rel_parent = Path(cohort_year)

    # Strip opening OME namespace declarations and any prefixed attributes (e.g. xsi:schemaLocation)
    # so standard ElementTree queries work reliably without unbound namespace errors across OME versions.
    xml_clean = re.sub(r'<OME[^>]*>', '<OME>', xml_string, count=1)
    xml_clean = re.sub(r'\s\w+:\w+="[^"]*"', '', xml_clean)
    root = ET.fromstring(xml_clean)

    images = root.findall(".//Image")
    for idx, img in enumerate(images):
        series_name = img.attrib.get("Name", f"Series_{idx}")
        pixels = img.find("Pixels")
        if pixels is None:
            continue

        size_x = int(pixels.attrib.get("SizeX", 1))
        size_y = int(pixels.attrib.get("SizeY", 1))
        size_z = int(pixels.attrib.get("SizeZ", 1))
        size_c = int(pixels.attrib.get("SizeC", 1))
        size_t = int(pixels.attrib.get("SizeT", 1))
        pixel_type = pixels.attrib.get("Type", "uint16")
        bits = int(pixels.attrib.get("SignificantBits", 16 if "16" in pixel_type else 8))

        def _get_float(attr_name: str, table: Dict[str, float]) -> Optional[float]:
            # Values are converted with their OME '<attr>Unit' (default um / s) so the inventory
            # columns really are in um and seconds (SlideBook DeltaT is usually in ms).
            val = convert_ome_unit(pixels.attrib.get(attr_name), pixels.attrib.get(f"{attr_name}Unit"), table)
            return round(val, 4) if val is not None else None

        phys_x = _get_float("PhysicalSizeX", LENGTH_TO_UM)
        phys_y = _get_float("PhysicalSizeY", LENGTH_TO_UM)
        phys_z = _get_float("PhysicalSizeZ", LENGTH_TO_UM)
        time_inc = _get_float("TimeIncrement", TIME_TO_S)

        # Fallback: compute mean delta_t from <Plane> elements if TimeIncrement attribute is omitted
        if time_inc is None and size_t > 1:
            planes = pixels.findall("Plane")
            t_deltas: Dict[int, float] = {}
            for p in planes:
                t_idx = p.attrib.get("TheT")
                delta_t = convert_ome_unit(p.attrib.get("DeltaT"), p.attrib.get("DeltaTUnit"), TIME_TO_S)
                if t_idx is not None and delta_t is not None:
                    try:
                        t_i = int(t_idx)
                        if t_i not in t_deltas:
                            t_deltas[t_i] = delta_t
                    except ValueError:
                        pass
            if len(t_deltas) >= 2:
                sorted_t = [t_deltas[k] for k in sorted(t_deltas.keys())]
                diffs = [sorted_t[i + 1] - sorted_t[i] for i in range(len(sorted_t) - 1)]
                if diffs:
                    time_inc = round(sum(diffs) / len(diffs), 3)

        est_gb = estimate_series_gb(size_x, size_y, size_z, size_c, size_t, bits)
        qualifies = (size_z >= min_z) and (size_t >= min_t) and (size_c >= min_c)
        calib_ok = all(v is not None and v > 0 for v in (phys_x, phys_y, phys_z, time_inc))

        safe_series = sanitize_filename(series_name)
        out_name = f"{sld_path.stem}__S{idx:02d}_{safe_series}.ome.tif"
        rel_out_path = str((rel_parent / out_name).as_posix()) if qualifies else ""

        results.append(
            SeriesMetadata(
                cohort_year=cohort_year,
                sld_path=str(sld_path),
                sld_filename=sld_path.name,
                sld_file_size_gb=file_size_gb,
                series_index=idx,
                series_name=series_name,
                size_x=size_x,
                size_y=size_y,
                size_z=size_z,
                size_c=size_c,
                size_t=size_t,
                pixel_type=pixel_type,
                bits_per_pixel=bits,
                phys_size_x_um=phys_x,
                phys_size_y_um=phys_y,
                phys_size_z_um=phys_z,
                time_interval_s=time_inc,
                estimated_uncompressed_gb=est_gb,
                is_4d_timelapse=qualifies,
                calibration_complete=calib_ok,
                extracted_ometiff_relpath=rel_out_path,
            )
        )
    return results


def inspect_sld_via_showinf(
    sld_path: Path,
    source_root: Path,
    showinf_bin: str,
    env_vars: Dict[str, str],
    min_z: int,
    min_t: int,
    min_c: int,
) -> List[SeriesMetadata]:
    """
    Executes Bio-Formats 'showinf -nopix -omexml-only' in a read-only subprocess to parse
    series headers without touching or modifying the source .sld file.
    """
    # -nopix: read header metadata only, skipping pixel blocks; -omexml-only: output structured XML
    cmd = [showinf_bin, "-nopix", "-omexml-only", str(sld_path)]
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env_vars,
        cwd=NATIVE_LOG_SINK_DIR, check=False,
    )
    # Decode bytes explicitly: Windows default text decoding crashes on the micro sign (0xb5 in 'µm').
    stdout_str = (proc.stdout or b"").decode("utf-8", errors="replace")
    stderr_str = (proc.stderr or b"").decode("utf-8", errors="replace")
    if proc.returncode != 0 or "<OME" not in stdout_str:
        err_msg = stderr_str.strip() or stdout_str.strip()
        if "Problem reading the Image Records" in err_msg or "The file is not a valid SlideBook document" in err_msg:
            has_dir = len(find_companion_sld_dirs(sld_path)) > 0
            dir_hint = " (companion data folder found)" if has_dir else " (likely aborted capture or companion folder missing)"
            raise ShowinfError(f"Corrupted Image Records / aborted acquisition{dir_hint}", err_msg)
        if "ArrayIndexOutOfBoundsException" in err_msg or "SlidebookReader" in err_msg:
            raise ShowinfError(
                f"Bio-Formats legacy parser failed on {sld_path.name}: 3i's SlideBook6Reader was not used. "
                "Either SlideBook6Reader.jar is not on BF_CP (expected in <tools>/readers), or the jar's "
                "native library failed to load on this OS (e.g. the 2026-09 Linux build is missing a "
                "symbol; setup_hpc_env.bash pins a working build). Run showinf -debug to see which.",
                err_msg,
            )
        raise ShowinfError(f"showinf failed (exit {proc.returncode}): {err_msg[:200]}", err_msg)
    xml_start = stdout_str.find("<OME")
    xml_str = stdout_str[xml_start:]
    return parse_ome_xml_metadata(xml_str, sld_path, source_root, min_z, min_t, min_c)


def find_companion_sld_dirs(sld_path: Path) -> List[Path]:
    """
    Finds any companion data directories associated with a SlideBook .sld or .sldy container.
    SlideBook 6+ frequently stores raw pixel streams in an adjacent folder sharing the stem name.
    """
    companions: List[Path] = []
    stem = sld_path.stem
    parent = sld_path.parent
    for candidate in [
        parent / f"{stem}_data",
        parent / f"{stem}.sld_Files",
        parent / f"{stem}_Files",
        parent / f"{stem}.dir",
    ]:
        if candidate.is_dir():
            companions.append(candidate)
    return companions


def make_health_record(
    sld_path: Path,
    source_root: Path,
    status: str,
    error_summary: str = "",
    total_series: int = 0,
    qualifying_4d_series: int = 0,
    probe_header: bool = True,
) -> FileHealthRecord:
    """
    Builds one sld_health_audit.csv row. Set probe_header=False for lock/temp files so we never
    open a file that SlideBook may hold an exclusive lock on.
    """
    try:
        rel_path = str(sld_path.relative_to(source_root).as_posix())
    except ValueError:
        rel_path = sld_path.name

    if probe_header:
        size_mb, mtime_utc, magic_hex, is_null = inspect_file_header(sld_path)
    else:
        try:
            st = sld_path.stat()
            size_mb = round(st.st_size / (1024 ** 2), 3)
            mtime_utc = datetime.datetime.fromtimestamp(st.st_mtime, tz=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        except OSError:
            size_mb, mtime_utc = 0.0, ""
        magic_hex, is_null = "", False

    return FileHealthRecord(
        cohort_year=infer_cohort_year(sld_path),
        rel_path=rel_path,
        sld_filename=sld_path.name,
        file_size_mb=size_mb,
        modified_utc=mtime_utc,
        has_companion_dirs=len(find_companion_sld_dirs(sld_path)) > 0,
        magic_hex=magic_hex,
        is_null_header=is_null,
        status=status,
        total_series=total_series,
        qualifying_4d_series=qualifying_4d_series,
        error_summary=error_summary,
    )


def audit_and_inspect_sld(
    sld_path: Path,
    source_root: Path,
    showinf_bin: str,
    env_vars: Dict[str, str],
    min_z: int,
    min_t: int,
    min_c: int,
) -> Tuple[List[SeriesMetadata], FileHealthRecord, str]:
    """
    Runs the fast byte-header check followed by showinf inspection. Never raises: returns
    (parsed series, health record, human-readable console hint). The hint is "" on success.
    """
    _, _, _, is_null = inspect_file_header(sld_path)
    if is_null:
        msg = "File is 0 bytes or first 4 KB are all null (0x00) bytes."
        return [], make_health_record(sld_path, source_root, "EMPTY_OR_NULL_FILE", msg), msg

    try:
        series_list = inspect_sld_via_showinf(
            sld_path, source_root, showinf_bin, env_vars, min_z, min_t, min_c
        )
    except ShowinfError as exc:
        status, summary = classify_showinf_error(exc.raw_stderr)
        return [], make_health_record(sld_path, source_root, status, summary), str(exc)
    except (ET.ParseError, ValueError, OSError) as exc:
        summary = f"{type(exc).__name__}: {exc}"[:300]
        return [], make_health_record(sld_path, source_root, "OME_XML_PARSE_ERR", summary), summary

    n_total = len(series_list)
    n_qual = sum(1 for s in series_list if s.is_4d_timelapse)
    if n_total == 0:
        status = "EMPTY_NO_IMAGES"
    elif n_qual > 0:
        status = "OK_4D"
    else:
        status = "OK_NON_4D"
    rec = make_health_record(sld_path, source_root, status, "", n_total, n_qual)
    return series_list, rec, ""


def extract_4d_series_to_ometiff(
    meta: SeriesMetadata,
    output_dir: Path,
    bfconvert_bin: str,
    env_vars: Dict[str, str],
    compression: str = "LZW",
    compute_checksum: bool = False,
    overwrite: bool = False,
    dry_run: bool = False,
) -> Tuple[bool, str, float, Optional[str]]:
    """
    Extracts a single qualifying 4D series from an .sld file into a losslessly compressed .ome.tif.
    Never modifies the source .sld file.
    """
    out_path = output_dir / meta.extracted_ometiff_relpath
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        return True, f"[DRY-RUN] Would extract series {meta.series_index} -> {out_path}", 0.0, None

    if out_path.exists() and not overwrite:
        size_gb = round(out_path.stat().st_size / (1024 ** 3), 4)
        sha = compute_sha256(out_path) if compute_checksum else None
        return True, f"Skipped (already exists): {out_path.name}", size_gb, sha

    # -series: isolate specific 4D experiment; -compression: lossless LZW/Deflate; -nooverwrite: safety guard
    cmd = [
        bfconvert_bin,
        "-series", str(meta.series_index),
        "-compression", compression,
        "-overwrite" if overwrite else "-nooverwrite",
        str(meta.sld_path),
        str(out_path),
    ]
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env_vars,
        cwd=NATIVE_LOG_SINK_DIR, check=False,
    )
    stderr_str = (proc.stderr or b"").decode("utf-8", errors="replace")
    if proc.returncode != 0 or not out_path.exists():
        return (
            False,
            f"Failed extracting series {meta.series_index} from {meta.sld_filename}: {stderr_str.strip()}",
            0.0,
            None,
        )

    size_gb = round(out_path.stat().st_size / (1024 ** 3), 4)
    sha = compute_sha256(out_path) if compute_checksum else None
    return True, f"Extracted: {out_path.name} ({size_gb:.3f} GB)", size_gb, sha


def copy_qualifying_sld(
    sld_path: Path,
    output_dir: Path,
    source_root: Path,
    use_hardlinks: bool = False,
    dry_run: bool = False,
) -> Tuple[bool, str, float]:
    """Copies (or hardlinks) a qualifying .sld file and any SlideBook companion folders to output_dir."""
    try:
        rel_parent = sld_path.parent.relative_to(source_root)
    except ValueError:
        rel_parent = Path(infer_cohort_year(sld_path))
    if rel_parent == Path("."):
        rel_parent = Path(infer_cohort_year(sld_path))

    target_folder = output_dir / rel_parent
    dest_file = target_folder / sld_path.name

    if dry_run:
        return True, f"[DRY-RUN] Would stage {sld_path.name} -> {dest_file}", 0.0

    target_folder.mkdir(parents=True, exist_ok=True)
    total_bytes = 0
    if not dest_file.exists():
        if use_hardlinks:
            os.link(sld_path, dest_file)
        else:
            shutil.copy2(sld_path, dest_file)
    total_bytes += dest_file.stat().st_size

    for comp_dir in find_companion_sld_dirs(sld_path):
        dest_comp = target_folder / comp_dir.name
        if not dest_comp.exists():
            shutil.copytree(
                comp_dir, dest_comp, copy_function=os.link if use_hardlinks else shutil.copy2
            )
        for root_d, _, files in os.walk(dest_comp):
            for f in files:
                total_bytes += (Path(root_d) / f).stat().st_size

    size_gb = round(total_bytes / (1024 ** 3), 4)
    return True, f"Staged .sld: {dest_file.name} ({size_gb:.3f} GB)", size_gb


def resolve_bftools_binaries(bftools_dir: Optional[str]) -> Tuple[str, str]:
    """Locates 'showinf' and 'bfconvert' executables across Windows, macOS, and Linux."""
    is_win = sys.platform.startswith("win")
    showinf_name = "showinf.bat" if is_win else "showinf"
    bfconvert_name = "bfconvert.bat" if is_win else "bfconvert"

    candidate_dirs: List[Path] = []
    if bftools_dir:
        candidate_dirs.append(Path(bftools_dir))

    # Auto-check user-space tool installs ($MITO_TOOLS_DIR on HPC, %LOCALAPPDATA% on the lab PC)
    for tools_root in _tools_dirs():
        candidate_dirs.append(tools_root / "bftools")

    # Auto-check relative bftools directories
    candidate_dirs.extend([
        Path.cwd() / "bftools",
        Path(__file__).resolve().parent / "bftools",
    ])

    for bdir in candidate_dirs:
        s_bin = bdir / showinf_name
        b_bin = bdir / bfconvert_name
        if s_bin.exists() and b_bin.exists():
            return str(s_bin), str(b_bin)

    s_found = shutil.which(showinf_name) or shutil.which("showinf")
    b_found = shutil.which(bfconvert_name) or shutil.which("bfconvert")
    return s_found or showinf_name, b_found or bfconvert_name


def get_tool_version(bin_path: str, env_vars: Dict[str, str]) -> str:
    """Queries 'showinf -version' for Checkpoint 2 provenance logging."""
    try:
        proc = subprocess.run(
            [bin_path, "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env_vars,
            timeout=15,  # JVM cold start on a busy lab PC can exceed 10 s
            cwd=NATIVE_LOG_SINK_DIR,
            check=False,
        )
        stdout_str = (proc.stdout or b"").decode("utf-8", errors="replace")
        for line in stdout_str.splitlines():
            if "Version" in line or "Bio-Formats" in line:
                return line.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown"


def summarize_health_by_year(records: List[FileHealthRecord]) -> Dict[str, Dict[str, int]]:
    """Aggregates file health status counts by cohort year for console and JSON reporting."""
    by_year: Dict[str, Dict[str, int]] = {}
    for rec in records:
        yr_map = by_year.setdefault(rec.cohort_year, {})
        yr_map["total"] = yr_map.get("total", 0) + 1
        yr_map[rec.status] = yr_map.get(rec.status, 0) + 1
    return dict(sorted(by_year.items()))


def manifest_entry(hpc_root: str, relpath: str) -> str:
    """Builds one movie_list.txt line: absolute if an HPC root is known, else relative to it."""
    return f"{hpc_root}/{relpath}" if hpc_root else relpath


SLD_EXTENSIONS = {".sld", ".sldy"}


def parse_shard_spec(spec: Optional[str]) -> Optional[Tuple[int, int]]:
    """
    Parses '--shard K/N' into (K, N) with 1 <= K <= N. Returns None when spec is empty.
    Raises ValueError for malformed specs so a bad array index never silently scans everything.
    """
    if not spec:
        return None
    match = re.fullmatch(r"\s*(\d+)\s*/\s*(\d+)\s*", spec)
    if not match:
        raise ValueError(f"--shard must look like K/N (e.g. 3/20), got '{spec}'")
    k, n = int(match.group(1)), int(match.group(2))
    if n < 1 or not 1 <= k <= n:
        raise ValueError(f"--shard K/N needs 1 <= K <= N, got '{spec}'")
    return k, n


def select_shard(paths: List[Path], shard_index: int, shard_count: int) -> List[Path]:
    """
    Deterministic round-robin split: shard K (1-based) of N takes items K-1, K-1+N, K-1+2N, ...
    of the given (already sorted) list. The N shards are disjoint and together cover every item
    exactly once; round-robin (not contiguous blocks) spreads large cohort folders across shards.
    """
    return [p for i, p in enumerate(paths) if i % shard_count == shard_index - 1]


def read_file_list(list_path: Path, source_dir: Path) -> List[Path]:
    """
    Reads a --file-list: one .sld/.sldy path per line, absolute or relative to source_dir.
    Blank lines and '#' comments are ignored; Windows CRLF is tolerated. Duplicates are dropped.
    Every path must lie inside source_dir (relative paths and cohort folders derive from it).
    """
    seen = set()
    files: List[Path] = []
    with open(list_path, "r", encoding="utf-8-sig") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            p = Path(line)
            if not p.is_absolute():
                p = source_dir / p
            p = p.resolve()
            try:
                p.relative_to(source_dir)
            except ValueError:
                raise ValueError(f"--file-list entry is outside source_dir {source_dir}: {line}")
            if p not in seen:
                seen.add(p)
                files.append(p)
    return files


def shard_report_name(filename: str, suffix: str) -> str:
    """Inserts the shard suffix before the extension: 'sld_4d_inventory.csv' -> 'sld_4d_inventory_shard3.csv'."""
    if not suffix:
        return filename
    stem, dot, ext = filename.rpartition(".")
    return f"{stem}{suffix}.{ext}" if dot else f"{filename}{suffix}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Filter 3i SlideBook (.sld) files for 4D time-lapse series, audit file health, and stage/extract for W&M SciClone HPC."
    )
    parser.add_argument("source_dir", type=Path, help="Read-only root directory containing raw .sld files")
    parser.add_argument("output_dir", type=Path, help="Destination staging directory (must differ from source_dir)")
    parser.add_argument(
        "--mode",
        choices=["scan-only", "extract-ometiff", "copy-sld"],
        default="scan-only",
        help="Action mode: 'scan-only' (write CSV + movie_list.txt), 'extract-ometiff' (extract 4D series as compressed OME-TIFF), or 'copy-sld'.",
    )
    parser.add_argument("--min-z", type=int, default=15, help="Minimum Z slices to count as 3D stack (default: 15)")
    parser.add_argument("--min-t", type=int, default=10, help="Minimum T timepoints to count as time-lapse (default: 10)")
    parser.add_argument("--min-c", type=int, default=2, help="Minimum fluorescence channels required (default: 2 for mito-EGFP + tagRFP-T)")
    parser.add_argument(
        "--compression",
        choices=["LZW", "Deflate", "Uncompressed"],
        default="LZW",
        help="Lossless compression for extracted .ome.tif files (default: LZW)",
    )
    parser.add_argument("--bftools-dir", type=str, default=None, help="Path to Bio-Formats bftools directory if not in PATH")
    parser.add_argument("--fiji-dir", type=str, default=None, help="Path to Fiji.app directory to auto-use FIJI's bundled Java")
    parser.add_argument(
        "--hpc-raw-root",
        type=str,
        default=os.environ.get("MITO_RAW_ROOT", ""),
        help=(
            "Absolute raw OME-TIFF root on SciClone to prefix movie_list.txt entries with "
            "(default: $MITO_RAW_ROOT if set; otherwise entries are written relative and resolved "
            "against $MITO_RAW_ROOT by the SLURM job)"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Parallel worker threads for scanning/extraction (use 1 to rule out SlideBook JNI contention; default: 4)",
    )
    parser.add_argument("--checksum", action="store_true", help="Compute SHA-256 checksums on extracted .ome.tif files")
    parser.add_argument("--hardlink", action="store_true", help="Use hardlinks instead of copying in 'copy-sld' mode")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing extracted files")
    parser.add_argument("--dry-run", action="store_true", help="Simulate extraction/staging without writing image files")
    parser.add_argument(
        "--plasmid-csv",
        type=Path,
        default=Path(os.environ["MITO_PLASMID_CSV"]) if os.environ.get("MITO_PLASMID_CSV") else None,
        help="Lab plasmid table CSV for the mito-construct filename filter (default: $MITO_PLASMID_CSV). "
        "See manuals/construct_filter_guide.md.",
    )
    parser.add_argument("--construct-lexicon", type=Path, default=None,
                        help="Lexicon JSON for the filter (default: config/construct_lexicon.json)")
    parser.add_argument("--mito-overrides", type=Path, default=None,
                        help="Manual overrides CSV (pattern,label,note); first match wins over heuristics")
    parser.add_argument(
        "--mito-filter",
        choices=["off", "report", "require"],
        default=None,
        help="off: no classification; report: classify and report only; require: only files labelled MITO "
        "are inspected/extracted, others get health status SKIPPED_NOT_MITO "
        "(default: require when --plasmid-csv is given, else off)",
    )
    parser.add_argument(
        "--shard",
        type=str,
        default=None,
        help="Process only shard K of N (format K/N, 1-based; deterministic round-robin over the sorted "
        "candidate list). Reports get a '_shardK' suffix. Used by Slurm job arrays.",
    )
    parser.add_argument(
        "--file-list",
        type=Path,
        default=None,
        help="Process only the .sld/.sldy files listed in this text file (one per line, absolute or "
        "relative to source_dir) instead of searching source_dir recursively.",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=None,
        help="Directory for the CSV/JSON/movie_list reports (default: output_dir). Extracted images "
        "always go to output_dir.",
    )

    args = parser.parse_args()
    if args.mito_filter is None:
        args.mito_filter = "require" if args.plasmid_csv else "off"
    try:
        shard = parse_shard_spec(args.shard)
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    shard_suffix = f"_shard{shard[0]}" if shard else ""
    orchestrated = shard is not None or args.file_list is not None

    if not args.source_dir.exists():
        print(f"[ERROR] Source directory does not exist: {args.source_dir}", file=sys.stderr)
        return 1

    # Strict Safety Guard: Prevent writing into the raw source directory
    if args.source_dir.resolve() == args.output_dir.resolve():
        print("[ERROR] Safety violation: output_dir must be distinct from read-only source_dir!", file=sys.stderr)
        return 1

    # Absolute paths: Bio-Formats subprocesses run from a temp directory (see NATIVE_LOG_SINK_DIR).
    args.source_dir = args.source_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_dir = (args.report_dir or args.output_dir).resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    showinf_bin, bfconvert_bin = resolve_bftools_binaries(args.bftools_dir)
    showinf_bin = str(Path(showinf_bin).resolve()) if Path(showinf_bin).exists() else showinf_bin
    bfconvert_bin = str(Path(bfconvert_bin).resolve()) if Path(bfconvert_bin).exists() else bfconvert_bin

    # Throwaway working directory for Bio-Formats subprocesses (absorbs 3i's native trace file).
    # TemporaryDirectory removes it when main() returns.
    global NATIVE_LOG_SINK_DIR
    native_log_sink = tempfile.TemporaryDirectory(prefix="mito_bf_cwd_")
    NATIVE_LOG_SINK_DIR = native_log_sink.name

    # Configure subprocess environment (auto-detecting FIJI Java and SlideBook 6/7 native reader if needed)
    sub_env = os.environ.copy()
    java_home = discover_fiji_java_home(args.fiji_dir)
    if java_home:
        sub_env["JAVA_HOME"] = java_home
        sub_env["BF_JAVA_HOME"] = java_home
        java_bin = str(Path(java_home) / "bin")
        if java_bin not in sub_env.get("PATH", ""):
            sub_env["PATH"] = f"{java_bin}{os.pathsep}{sub_env.get('PATH', '')}"

    bftools_bin_dir = str(Path(showinf_bin).resolve().parent)
    if bftools_bin_dir not in sub_env.get("PATH", ""):
        sub_env["PATH"] = f"{bftools_bin_dir}{os.pathsep}{sub_env.get('PATH', '')}"

    sb_jar = discover_slidebook_reader_jar(args.bftools_dir, showinf_bin, args.fiji_dir)
    if sb_jar:
        curr_cp = sub_env.get("BF_CP", "")
        if sb_jar not in curr_cp:
            sub_env["BF_CP"] = f"{sb_jar}{os.pathsep}{curr_cp}".rstrip(os.pathsep)
        print(f"Using 3i SlideBook native reader JAR: {sb_jar}")

    if args.file_list is not None:
        try:
            listed = read_file_list(args.file_list, args.source_dir)
        except (OSError, ValueError) as exc:
            print(f"[ERROR] Cannot read --file-list {args.file_list}: {exc}", file=sys.stderr)
            return 1
        # Listed files that vanished still get a health row (UNREADABLE_FILE) below via stat().
        raw_candidates = sorted(p for p in listed if p.suffix.lower() in SLD_EXTENSIONS)
    else:
        raw_candidates = sorted(
            [p for p in args.source_dir.rglob("*") if p.is_file() and p.suffix.lower() in SLD_EXTENSIONS]
        )
    n_candidates_total = len(raw_candidates)
    if shard:
        raw_candidates = select_shard(raw_candidates, shard[0], shard[1])
        print(f"[INFO] Shard {shard[0]}/{shard[1]}: {len(raw_candidates)} of {n_candidates_total} candidate file(s).")
    if not raw_candidates and not orchestrated:
        print(f"[WARNING] No .sld or .sldy files found under {args.source_dir}")
        return 0
    # Orchestrated (sharded / file-list) runs continue with an empty list so every shard still writes
    # its (empty) reports; merge_collect_shards.py then knows the shard finished.

    sld_files: List[Path] = []
    health_records: List[FileHealthRecord] = []

    for p in raw_candidates:
        # Pre-filter: skip lock and temporary hidden files (never opened, SlideBook may hold a lock)
        if p.name.startswith(("~", ".~", "._", ".#", "temp")):
            health_records.append(make_health_record(
                p, args.source_dir, "SKIPPED_TEMP_OR_LOCK",
                "Temporary lock, backup, or hidden file prefix (~, .~, ._, .#, temp)",
                probe_header=False,
            ))
            continue
        try:
            sz = p.stat().st_size
        except OSError as e:
            health_records.append(make_health_record(
                p, args.source_dir, "UNREADABLE_FILE", f"OS file access error: {e}", probe_header=False,
            ))
            continue

        # Pre-filter: 0-byte or incomplete stubs < 64 KB (valid SlideBook containers are multi-megabyte/gigabyte)
        if sz < 64 * 1024:
            health_records.append(make_health_record(
                p, args.source_dir, "SKIPPED_STUB",
                f"Incomplete or zero-byte stub ({sz} bytes < 64 KB threshold)",
            ))
            continue

        sld_files.append(p)

    if health_records:
        print(f"[INFO] Pre-filtered {len(health_records)} non-data/stub file(s) (< 64 KB, temp locks, or unreadable).")

    # Mito-construct filename filter (names only: runs before any Bio-Formats call, so non-mito files
    # cost nothing under --mito-filter require).
    mito_calls: Dict[str, object] = {}
    if args.mito_filter != "off":
        try:
            import construct_classifier as cc
            lexicon = cc.load_lexicon(args.construct_lexicon)
            table = cc.load_construct_table(args.plasmid_csv, lexicon) if args.plasmid_csv else None
            overrides = cc.load_overrides(args.mito_overrides)
        except Exception as exc:  # bad CSV/lexicon/overrides must stop the run, not silently disable filtering
            print(f"[ERROR] Mito-construct filter setup failed: {exc}", file=sys.stderr)
            return 2
        kept: List[Path] = []
        for p in sld_files:
            rel = p.relative_to(args.source_dir).as_posix()
            call = cc.classify_name(rel, table, lexicon, overrides)
            mito_calls[str(p)] = call
            if args.mito_filter == "require" and call.label != "MITO":
                health_records.append(make_health_record(
                    p, args.source_dir, "SKIPPED_NOT_MITO",
                    f"{call.label}: {call.reasons}", probe_header=False,
                ))
            else:
                kept.append(p)
        if shard_suffix:
            # write_reports uses fixed file names: write into a private temp folder, then move the
            # files to their '_shardK' names so concurrent shards never clobber each other.
            tmp_report = Path(tempfile.mkdtemp(prefix=f".mito_reports{shard_suffix}_", dir=str(report_dir)))
            try:
                for written in cc.write_reports(list(mito_calls.values()), tmp_report):
                    os.replace(written, report_dir / shard_report_name(Path(written).name, shard_suffix))
            finally:
                shutil.rmtree(tmp_report, ignore_errors=True)
        else:
            cc.write_reports(list(mito_calls.values()), report_dir)
        label_counts: Dict[str, int] = {}
        for call in mito_calls.values():
            label_counts[call.label] = label_counts.get(call.label, 0) + 1
        print(f"[INFO] Mito-construct filter ({args.mito_filter}): {label_counts}; "
              f"see mito_classification{shard_suffix}.csv / mito_review{shard_suffix}.csv in {report_dir}")
        if args.mito_filter == "require":
            sld_files = kept

    total_raw_gb = sum(p.stat().st_size for p in sld_files) / (1024 ** 3)
    print(f"Found {len(sld_files)} valid SlideBook container(s) ({total_raw_gb:.2f} GB total raw container size).")
    print(f"Filtering criteria: SizeZ >= {args.min_z}, SizeT >= {args.min_t}, SizeC >= {args.min_c}")
    if java_home:
        print(f"Using Java runtime: {java_home}")

    all_series: List[SeriesMetadata] = []
    qualifying_series: List[SeriesMetadata] = []
    qualifying_sld_paths = set()

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        future_map = {
            pool.submit(
                audit_and_inspect_sld,
                sld,
                args.source_dir,
                showinf_bin,
                sub_env,
                args.min_z,
                args.min_t,
                args.min_c,
            ): sld
            for sld in sld_files
        }
        for future in concurrent.futures.as_completed(future_map):
            sld = future_map[future]
            series_list, health_rec, hint = future.result()
            health_records.append(health_rec)
            if health_rec.status not in HEALTHY_STATUSES:
                print(
                    f"[WARN] [{health_rec.status}] Skipping {sld.name} ({health_rec.file_size_mb:.1f} MB): {hint}",
                    file=sys.stderr,
                )
            call = mito_calls.get(str(sld))
            if call is not None:
                for s in series_list:
                    s.mito_label, s.mito_fp = call.label, call.mito_fp
                    s.mito_color, s.red_mito_flag = call.mito_color, call.red_mito_flag
            all_series.extend(series_list)
            for s in series_list:
                if s.is_4d_timelapse:
                    qualifying_series.append(s)
                    qualifying_sld_paths.add(Path(s.sld_path))

    all_series.sort(key=lambda s: (s.sld_path, s.series_index))
    qualifying_series.sort(key=lambda s: (s.sld_path, s.series_index))
    health_records.sort(key=lambda r: (r.cohort_year, r.rel_path))
    n_unhealthy = sum(1 for r in health_records if r.status not in HEALTHY_STATUSES)

    # 1. Write Master CSV Inventory
    manifest_csv = report_dir / shard_report_name("sld_4d_inventory.csv", shard_suffix)
    if all_series or orchestrated:
        # Orchestrated shards always write the header, even when empty, so the merge sees them.
        fieldnames = list(SeriesMetadata.__dataclass_fields__.keys())
        with open(manifest_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for s in all_series:
                writer.writerow(asdict(s))

    # 1b. Write file-level health audit (one row per candidate file: pass, skip, or failure)
    health_csv = report_dir / shard_report_name("sld_health_audit.csv", shard_suffix)
    with open(health_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(FileHealthRecord.__dataclass_fields__.keys()))
        writer.writeheader()
        for rec in health_records:
            writer.writerow(asdict(rec))

    # 2. Write SLURM Job Array Manifests (movie_list.txt + per-year lists)
    hpc_root = (args.hpc_raw_root or "").rstrip("/")
    master_list_path = report_dir / shard_report_name("movie_list.txt", shard_suffix)
    by_year: Dict[str, List[str]] = {}
    with open(master_list_path, "w", encoding="utf-8") as f:
        for s in qualifying_series:
            hpc_path = manifest_entry(hpc_root, s.extracted_ometiff_relpath)
            f.write(f"{hpc_path}\n")
            by_year.setdefault(s.cohort_year, []).append(hpc_path)

    for yr, paths in by_year.items():
        yr_list_path = report_dir / shard_report_name(f"movie_list_{yr}.txt", shard_suffix)
        with open(yr_list_path, "w", encoding="utf-8") as f:
            for p in paths:
                f.write(f"{p}\n")

    health_by_year = summarize_health_by_year(health_records)
    est_4d_uncompressed_gb = sum(s.estimated_uncompressed_gb for s in qualifying_series)

    # 3. Write Checkpoint 2 Provenance JSON
    provenance_log = {
        "timestamp_iso": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "host_os": platform.platform(),
        "python_version": platform.python_version(),
        "bioformats_version": get_tool_version(showinf_bin, sub_env),
        "java_home_used": java_home or "system_default",
        "source_dir": str(args.source_dir.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "filter_criteria": {
            "min_z": args.min_z,
            "min_t": args.min_t,
            "min_c": args.min_c,
            "compression": args.compression,
            "mode": args.mode,
            "workers": args.workers,
        },
        "slidebook_reader_jar": sb_jar or "not_found",
        "report_dir": str(report_dir),
        "shard": {"index": shard[0], "count": shard[1]} if shard else None,
        "file_list": str(args.file_list.resolve()) if args.file_list else None,
        "movie_list_root": hpc_root or "relative (resolve against $MITO_RAW_ROOT)",
        "summary_metrics": {
            "total_sld_candidates_found": len(raw_candidates),
            "total_sld_candidates_before_sharding": n_candidates_total,
            "total_sld_files_scanned": len(sld_files),
            "total_unparseable_or_skipped": n_unhealthy,
            "total_raw_sld_gb": round(total_raw_gb, 3),
            "total_series_inspected": len(all_series),
            "qualifying_4d_series_count": len(qualifying_series),
            "qualifying_sld_files_count": len(qualifying_sld_paths),
            "estimated_uncompressed_4d_gb": round(est_4d_uncompressed_gb, 3),
            "cohorts_discovered": {yr: len(paths) for yr, paths in by_year.items()},
        },
        "health_audit_by_year": health_by_year,
        "mito_construct_filter": {
            "mode": args.mito_filter,
            "plasmid_csv": str(args.plasmid_csv) if args.plasmid_csv else None,
            "overrides_csv": str(args.mito_overrides) if args.mito_overrides else None,
            "label_counts": {
                lab: sum(1 for c in mito_calls.values() if c.label == lab)
                for lab in sorted({c.label for c in mito_calls.values()})
            },
            "red_mito_files": sorted(Path(k).name for k, c in mito_calls.items() if c.red_mito_flag),
        },
    }
    prov_json_path = report_dir / shard_report_name("extraction_provenance.json", shard_suffix)
    with open(prov_json_path, "w", encoding="utf-8") as f:
        json.dump(provenance_log, f, indent=2)

    print("\n--- Scan & Manifest Summary ---")
    print(f"Total .sld files scanned          : {len(sld_files)} ({total_raw_gb:.2f} GB)")
    print(f"Unparseable or skipped files      : {n_unhealthy} of {len(raw_candidates)} candidate(s)")
    print(f"Total internal series inspected   : {len(all_series)}")
    print(f"Qualifying 4D time-lapse series   : {len(qualifying_series)} across {len(qualifying_sld_paths)} .sld file(s)")
    print(f"Estimated uncompressed 4D data    : {est_4d_uncompressed_gb:.2f} GB (~{est_4d_uncompressed_gb * 0.55:.2f} GB with lossless LZW)")
    print(f"Inventory CSV saved to            : {manifest_csv}")
    print(f"File health audit saved to        : {health_csv}")
    print(f"SLURM array list saved to         : {master_list_path}"
          + ("" if hpc_root else " (relative paths; resolved against $MITO_RAW_ROOT)"))
    print(f"Checkpoint 2 provenance saved to  : {prov_json_path}")
    if mito_calls:
        n_red = sum(1 for c in mito_calls.values() if c.red_mito_flag)
        print(f"Mito-construct filter ({args.mito_filter:7s}): "
              + ", ".join(f"{lab}={n}" for lab, n in provenance_log["mito_construct_filter"]["label_counts"].items())
              + (f"; {n_red} red-mito file(s): check channel_mito/channel_cell order" if n_red else ""))

    print("\n--- File Health by Cohort Year ---")
    for yr, counts in health_by_year.items():
        breakdown = ", ".join(f"{k}={v}" for k, v in counts.items() if k != "total")
        print(f"  {yr}: total={counts.get('total', 0)} ({breakdown})")

    if args.mode == "scan-only":
        return 0

    if args.mode == "extract-ometiff":
        print(f"\nExtracting {len(qualifying_series)} qualifying 4D series to {args.compression}-compressed .ome.tif...")
        total_out_gb = 0.0
        n_extract_failed = 0
        checksums: Dict[str, str] = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(
                    extract_4d_series_to_ometiff,
                    meta,
                    args.output_dir,
                    bfconvert_bin,
                    sub_env,
                    args.compression,
                    args.checksum,
                    args.overwrite,
                    args.dry_run,
                ): meta
                for meta in qualifying_series
            }
            for fut in concurrent.futures.as_completed(futures):
                meta = futures[fut]
                ok, msg, sz, sha = fut.result()
                print(f"  [{'OK' if ok else 'ERR'}] {msg}")
                n_extract_failed += 0 if ok else 1
                total_out_gb += sz
                if sha:
                    checksums[meta.extracted_ometiff_relpath] = sha

        if checksums:
            sha_path = report_dir / shard_report_name("extracted_checksums.sha256", shard_suffix)
            with open(sha_path, "w", encoding="utf-8") as f:
                for rel_p, digest in sorted(checksums.items()):
                    f.write(f"{digest}  {rel_p}\n")
        print(f"\nCompleted extraction: {total_out_gb:.2f} GB written to {args.output_dir}")
        if n_extract_failed and orchestrated:
            # Slurm array tasks must report failure so 'submit_full_run.bash retry' can resubmit them.
            print(f"[ERROR] {n_extract_failed} series failed to extract (see [ERR] lines above).", file=sys.stderr)
            return 4

    elif args.mode == "copy-sld":
        print(f"\nStaging {len(qualifying_sld_paths)} qualifying .sld file(s) to {args.output_dir}...")
        total_staged_gb = 0.0
        for sld_path in sorted(qualifying_sld_paths):
            ok, msg, sz = copy_qualifying_sld(
                sld_path, args.output_dir, args.source_dir, args.hardlink, args.dry_run
            )
            print(f"  [{'OK' if ok else 'ERR'}] {msg}")
            total_staged_gb += sz
        print(f"\nCompleted staging: {total_staged_gb:.2f} GB written to {args.output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())