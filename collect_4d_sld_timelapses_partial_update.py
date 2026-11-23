#!/usr/bin/env python3
"""
collect_4d_sld_timelapses.py
Phase 1 & 2 Utility: SlideBook (.sld / .sldy) 4D Time-Lapse Discovery, Health Audit, Compression & Staging
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
    """Audit record written to sld_health_audit.csv for every scanned .sld/.sldy file (pass or fail)."""
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


def sanitize_filename(name: str) -> str:
    """Converts a SlideBook series name into a clean, POSIX-safe filename stem."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return re.sub(r"_+", "_", cleaned).strip("_") or "series"


def infer_cohort_year(sld_path: Path) -> str:
    """Infers the experimental cohort year (e.g., '2023', '2024') from folder or filename."""
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
    """Categorizes showinf failures into actionable buckets for sld_health_audit.csv."""
    cleaned = " | ".join(line.strip() for line in stderr_text.splitlines() if line.strip())
    short_err = cleaned[:300]
    low = cleaned.lower()
    if "not a valid slidebook document" in low or "problem reading the image records" in low:
        return "SLIDEBOOK_IMAGE_RECORDS_ERR", short_err
    if "unsatisfiedlinkerror" in low or "noclassdeffounderror" in low or "outofmemoryerror" in low:
        return "JNI_OR_JAVA_ERR", short_err
    if "file already open" in low or "file not open" in low:
        return "SLIDEBOOK_SESSION_CONCURRENCY_ERR", short_err
    return "OTHER_ERR", short_err


def discover_fiji_java_home(explicit_fiji_dir: Optional[str] = None) -> Optional[str]:
    """Locates the OpenJDK Java runtime bundled inside FIJI (Fiji.app/java/...)."""
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

    is_win = sys.platform.startswith("win")
    java_exe_name = "java.exe" if is_win else "java"

    for root in candidate_roots:
        java_dir = root / "java"
        if java_dir.is_dir():
            for java_bin in java_dir.rglob(java_exe_name):
                if java_bin.is_file() and java_bin.parent.name == "bin":
                    return str(java_bin.parent.parent)
    return None


def parse_ome_xml_metadata(
    xml_string: str,
    sld_path: Path,
    source_root: Path,
    min_z: int,
    min_t: int,
    min_c: int,
) -> List[SeriesMetadata]:
    """Parses OME-XML metadata returned by Bio-Formats 'showinf -nopix -omexml-only'."""
    results: List[SeriesMetadata] = []
    file_size_gb = round(sld_path.stat().st_size / (1024 ** 3), 4) if sld_path.exists() else 0.0
    cohort_year = infer_cohort_year(sld_path)

    try:
        rel_parent = sld_path.parent.relative_to(source_root)
    except ValueError:
        rel_parent = Path(cohort_year)

    xml_clean = re.sub(r'\sxmlns(:\w+)?="[^"]+"', "", xml_string)
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
    """Executes Bio-Formats 'showinf -nopix -omexml-only' in a read-only subprocess."""
    cmd = [showinf_bin, "-nopix", "-omexml-only", str(sld_path)]
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env_vars, check=False
    )
    if proc.returncode != 0 or "<OME" not in proc.stdout:
        raise RuntimeError(
            f"showinf failed on {sld_path.name} (exit {proc.returncode}): {proc.stderr.strip()}"
        )
    xml_start = proc.stdout.find("<OME")
    xml_str = proc.stdout[xml_start:]
    return parse_ome_xml_metadata(xml_str, sld_path, source_root, min_z, min_t, min_c)


def find_companion_sld_dirs(sld_path: Path) -> List[Path]:
    """Finds any companion data directories associated with a SlideBook .sld or .sldy container."""
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


def audit_and_inspect_sld(
    sld_path: Path,
    source_root: Path,
    showinf_bin: str,
    env_vars: Dict[str, str],
    min_z: int,
    min_t: int,
    min_c: int,
) -> Tuple[List[SeriesMetadata], FileHealthRecord]:
    """
    Runs the fast byte-header check followed by showinf inspection, returning both
    parsed SeriesMetadata (if valid) and a structured FileHealthRecord for every file.
    """
    cohort_year = infer_cohort_year(sld_path)
    try:
        rel_path = str(sld_path.relative_to(source_root).as_posix())
    except ValueError:
        rel_path = sld_path.name

    size_mb, mtime_utc, magic_hex, is_null = inspect_file_header(sld_path)
    has_companions = len(find_companion_sld_dirs(sld_path)) > 0

    if is_null:
        rec = FileHealthRecord(
            cohort_year=cohort_year,
            rel_path=rel_path,
            sld_filename=sld_path.name,
            file_size_mb=size_mb,
            modified_utc=mtime_utc,
            has_companion_dirs=has_companions,
            magic_hex=magic_hex,
            is_null_header=True,
            status="EMPTY_OR_NULL_FILE",
            total_series=0,
            qualifying_4d_series=0,
            error_summary="File is 0 bytes or first 4 KB are all null (0x00) bytes.",
        )
        return [], rec

    try:
        series_list = inspect_sld_via_showinf(
            sld_path, source_root, showinf_bin, env_vars, min_z, min_t, min_c
        )
        n_total = len(series_list)
        n_qual = sum(1 for s in series_list if s.is_4d_timelapse)
        if n_total == 0:
            status = "EMPTY_NO_IMAGES"
        elif n_qual > 0:
            status = "OK_4D"
        else:
            status = "OK_NON_4D"
        rec = FileHealthRecord(
            cohort_year=cohort_year,
            rel_path=rel_path,
            sld_filename=sld_path.name,
            file_size_mb=size_mb,
            modified_utc=mtime_utc,
            has_companion_dirs=has_companions,
            magic_hex=magic_hex,
            is_null_header=False,
            status=status,
            total_series=n_total,
            qualifying_4d_series=n_qual,
            error_summary="",
        )
        return series_list, rec
    except Exception as exc:
        err_status, err_summary = classify_showinf_error(str(exc))
        rec = FileHealthRecord(
            cohort_year=cohort_year,
            rel_path=rel_path,
            sld_filename=sld_path.name,
            file_size_mb=size_mb,
            modified_utc=mtime_utc,
            has_companion_dirs=has_companions,
            magic_hex=magic_hex,
            is_null_header=False,
            status=err_status,
            total_series=0,
            qualifying_4d_series=0,
            error_summary=err_summary,
        )
        return [], rec


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
    """Extracts a single qualifying 4D series from an .sld file into a losslessly compressed .ome.tif."""
    out_path = output_dir / meta.extracted_ometiff_relpath
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        return True, f"[DRY-RUN] Would extract series {meta.series_index} -> {out_path}", 0.0, None

    if out_path.exists() and not overwrite:
        size_gb = round(out_path.stat().st_size / (1024 ** 3), 4)
        sha = compute_sha256(out_path) if compute_checksum else None
        return True, f"Skipped (already exists): {out_path.name}", size_gb, sha

    cmd = [
        bfconvert_bin,
        "-series", str(meta.series_index),
        "-compression", compression,
        "-overwrite" if overwrite else "-nooverwrite",
        str(meta.sld_path),
        str(out_path),
    ]
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env_vars, check=False
    )
    if proc.returncode != 0 or not out_path.exists():
        return (
            False,
            f"Failed extracting series {meta.series_index} from {meta.sld_filename}: {proc.stderr.strip()}",
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

    if bftools_dir:
        bdir = Path(bftools_dir)
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
            text=True,
            env=env_vars,
            timeout=15,
            check=False,
        )
        for line in proc.stdout.splitlines():
            if "Version" in line or "Bio-Formats" in line:
                return line.strip()
    except Exception:
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
        help="Action mode: 'scan-only', 'extract-ometiff', or 'copy-sld'.",
    )
    parser.add_argument("--min-z", type=int, default=15, help="Minimum Z slices to count as 3D stack (default: 15)")
    parser.add_argument("--min-t", type=int, default=10, help="Minimum T timepoints to count as time-lapse (default: 10)")
    parser.add_argument("--min-c", type=int, default=2, help="Minimum fluorescence channels required (default: 2)")
    parser.add_argument(
        "--compression",
        choices=["LZW", "Deflate", "Uncompressed"],
        default="LZW",
        help="Lossless compression for extracted .ome.tif files (default: LZW)",
    )
    parser.add_argument("--bftools-dir", type=str, default=None, help="Path to Bio-Formats bftools directory")
    parser.add_argument("--fiji-dir", type=str, default=None, help="Path to Fiji.app directory")
    parser.add_argument(
        "--hpc-raw-root",
        type=str,
        default="/sciclone/data10/bestmanlab/raw_ometiff",
        help="Target root path on SciClone HPC used when generating movie_list.txt",
    )
    parser.add_argument("--workers", type=int, default=4, help="Parallel worker threads (use 1 to rule out JNI contention; default: 4)")
    parser.add_argument("--checksum", action="store_true", help="Compute SHA-256 checksums on extracted files")
    parser.add_argument("--hardlink", action="store_true", help="Use hardlinks in 'copy-sld' mode")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing extracted files")
    parser.add_argument("--dry-run", action="store_true", help="Simulate extraction/staging without writing image files")

    args = parser.parse_args()

    if not args.source_dir.exists():
        print(f"[ERROR] Source directory does not exist: {args.source_dir}", file=sys.stderr)
        return 1

    if args.source_dir.resolve() == args.output_dir.resolve():
        print("[ERROR] Safety violation: output_dir must be distinct from read-only source_dir!", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    showinf_bin, bfconvert_bin = resolve_bftools_binaries(args.bftools_dir)

    sub_env = os.environ.copy()
    java_home = discover_fiji_java_home(args.fiji_dir)
    if java_home:
        sub_env["JAVA_HOME"] = java_home
        sub_env["BF_JAVA_HOME"] = java_home

    sld_files = sorted(
        [p for p in args.source_dir.rglob("*") if p.is_file() and p.suffix.lower() in {".sld", ".sldy"}]
    )
    if not sld_files:
        print(f"[WARNING] No .sld or .sldy files found under {args.source_dir}")
        return 0

    total_raw_gb = sum(p.stat().st_size for p in sld_files) / (1024 ** 3)
    all_series: List[SeriesMetadata] = []
    qualifying_series: List[SeriesMetadata] = []
    qualifying_sld_paths = set()
    health_records: List[FileHealthRecord] = []

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
            series_list, health_rec = future.result()
            health_records.append(health_rec)
            if health_rec.status not in {"OK_4D", "OK_NON_4D"}:
                print(
                    f"[WARN] [{health_rec.status}] {sld.name} ({health_rec.file_size_mb} MB): {health_rec.error_summary}",
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

    # Write complete file-level health audit CSV (captures every pass AND failure).
    health_csv = args.output_dir / "sld_health_audit.csv"
    if health_records:
        health_fields = list(asdict(health_records[0]).keys())
        with open(health_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=health_fields)
            writer.writeheader()
            for rec in health_records:
                writer.writerow(asdict(rec))

    manifest_csv = args.output_dir / "sld_4d_inventory.csv"
    if all_series:
        fieldnames = list(asdict(all_series[0]).keys())
        with open(manifest_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for s in all_series:
                writer.writerow(asdict(s))

    hpc_root = args.hpc_raw_root.rstrip("/")
    master_list_path = args.output_dir / "movie_list.txt"
    by_year: Dict[str, List[str]] = {}
    with open(master_list_path, "w", encoding="utf-8") as f:
        for s in qualifying_series:
            hpc_path = f"{hpc_root}/{s.extracted_ometiff_relpath}"
            f.write(f"{hpc_path}\n")
            by_year.setdefault(s.cohort_year, []).append(hpc_path)

    for yr, paths in by_year.items():
        yr_list_path = args.output_dir / f"movie_list_{yr}.txt"
        with open(yr_list_path, "w", encoding="utf-8") as f:
            for p in paths:
                f.write(f"{p}\n")

    health_by_year = summarize_health_by_year(health_records)
    est_4d_uncompressed_gb = sum(s.estimated_uncompressed_gb for s in qualifying_series)
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
        "summary_metrics": {
            "total_sld_files_scanned": len(sld_files),
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

    print("\n=== SlideBook (.sld) Health Audit Summary by Cohort Year ===")
    for yr, counts in health_by_year.items():
        breakdown = ", ".join(f"{k}={v}" for k, v in counts.items() if k != "total")
        print(f"  {yr}: total={counts.get('total', 0)} ({breakdown})")
    print(f"Full file-level audit saved to: {health_csv}")

    if args.mode == "scan-only":
        return 0

    if args.mode == "extract-ometiff":
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

    elif args.mode == "copy-sld":
        for sld_path in sorted(qualifying_sld_paths):
            ok, msg, _ = copy_qualifying_sld(
                sld_path, args.output_dir, args.source_dir, args.hardlink, args.dry_run
            )
            print(f"  [{'OK' if ok else 'ERR'}] {msg}")

    return 0


if __name__ == "__main__":
    sys.exit(main())