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


# Status values that mean showinf parsed the file successfully.
HEALTHY_STATUSES = {"OK_4D", "OK_NON_4D"}


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

        def _get_float(attr_name: str) -> Optional[float]:
            val = pixels.attrib.get(attr_name)
            try:
                return round(float(val), 4) if val is not None else None
            except ValueError:
                return None

        phys_x = _get_float("PhysicalSizeX")
        phys_y = _get_float("PhysicalSizeY")
        phys_z = _get_float("PhysicalSizeZ")
        time_inc = _get_float("TimeIncrement")

        # Fallback: compute median delta_t from <Plane> elements if TimeIncrement attribute is omitted
        if time_inc is None and size_t > 1:
            planes = pixels.findall("Plane")
            t_deltas: Dict[int, float] = {}
            for p in planes:
                t_idx = p.attrib.get("TheT")
                delta_t = p.attrib.get("DeltaT")
                if t_idx is not None and delta_t is not None:
                    try:
                        t_i = int(t_idx)
                        if t_i not in t_deltas:
                            t_deltas[t_i] = float(delta_t)
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
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env_vars, check=False
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
                f"Bio-Formats legacy parser failed on {sld_path.name}. "
                "This SlideBook 6+ file requires 3i's native reader plugin (SlideBook6Reader.jar). "
                "Place SlideBook6Reader.jar in <tools>/readers or set BF_CP "
                "(scripts/setup_lab_windows.ps1 installs it automatically).",
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
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env_vars, check=False
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

    args = parser.parse_args()

    if not args.source_dir.exists():
        print(f"[ERROR] Source directory does not exist: {args.source_dir}", file=sys.stderr)
        return 1

    # Strict Safety Guard: Prevent writing into the raw source directory
    if args.source_dir.resolve() == args.output_dir.resolve():
        print("[ERROR] Safety violation: output_dir must be distinct from read-only source_dir!", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    showinf_bin, bfconvert_bin = resolve_bftools_binaries(args.bftools_dir)

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

    raw_candidates = sorted(
        [p for p in args.source_dir.rglob("*") if p.is_file() and p.suffix.lower() in {".sld", ".sldy"}]
    )
    if not raw_candidates:
        print(f"[WARNING] No .sld or .sldy files found under {args.source_dir}")
        return 0

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
    manifest_csv = args.output_dir / "sld_4d_inventory.csv"
    if all_series:
        fieldnames = list(asdict(all_series[0]).keys())
        with open(manifest_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for s in all_series:
                writer.writerow(asdict(s))

    # 1b. Write file-level health audit (one row per candidate file: pass, skip, or failure)
    health_csv = args.output_dir / "sld_health_audit.csv"
    with open(health_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(FileHealthRecord.__dataclass_fields__.keys()))
        writer.writeheader()
        for rec in health_records:
            writer.writerow(asdict(rec))

    # 2. Write SLURM Job Array Manifests (movie_list.txt + per-year lists)
    hpc_root = (args.hpc_raw_root or "").rstrip("/")
    master_list_path = args.output_dir / "movie_list.txt"
    by_year: Dict[str, List[str]] = {}
    with open(master_list_path, "w", encoding="utf-8") as f:
        for s in qualifying_series:
            hpc_path = manifest_entry(hpc_root, s.extracted_ometiff_relpath)
            f.write(f"{hpc_path}\n")
            by_year.setdefault(s.cohort_year, []).append(hpc_path)

    for yr, paths in by_year.items():
        yr_list_path = args.output_dir / f"movie_list_{yr}.txt"
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
        "movie_list_root": hpc_root or "relative (resolve against $MITO_RAW_ROOT)",
        "summary_metrics": {
            "total_sld_candidates_found": len(raw_candidates),
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
    }
    prov_json_path = args.output_dir / "extraction_provenance.json"
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

    print("\n--- File Health by Cohort Year ---")
    for yr, counts in health_by_year.items():
        breakdown = ", ".join(f"{k}={v}" for k, v in counts.items() if k != "total")
        print(f"  {yr}: total={counts.get('total', 0)} ({breakdown})")

    if args.mode == "scan-only":
        return 0

    if args.mode == "extract-ometiff":
        print(f"\nExtracting {len(qualifying_series)} qualifying 4D series to {args.compression}-compressed .ome.tif...")
        total_out_gb = 0.0
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
                total_out_gb += sz
                if sha:
                    checksums[meta.extracted_ometiff_relpath] = sha

        if checksums:
            sha_path = args.output_dir / "extracted_checksums.sha256"
            with open(sha_path, "w", encoding="utf-8") as f:
                for rel_p, digest in sorted(checksums.items()):
                    f.write(f"{digest}  {rel_p}\n")
        print(f"\nCompleted extraction: {total_out_gb:.2f} GB written to {args.output_dir}")

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