#!/usr/bin/env python3
"""
merge_collect_shards.py
Combines the per-shard reports of a sharded collect_4d_sld_timelapses.py run into one set.

A Slurm array runs the collector as N tasks ('--shard K/N'), each writing reports with a '_shardK'
suffix into one shard directory. This script rebuilds exactly what a single unsharded run writes:

  sld_4d_inventory.csv         all series, sorted by (sld_path, series_index)
  sld_health_audit.csv         one row per candidate file, sorted by (cohort_year, rel_path)
  movie_list.txt               qualifying 4D series, absolute or relative (same rule as the collector)
  movie_list_<year>.txt        per-cohort-year lists
  mito_classification.csv      mito-construct filename calls (when the filter was on)
  mito_review.csv              REVIEW calls, high priority first
  extracted_checksums.sha256   SHA-256 of extracted OME-TIFFs (extract phase with --checksum)
  extraction_provenance.json   merged provenance: shard list, missing shards, summed metrics

Safety rules (the merge REFUSES, exit 2, instead of guessing):
  - the same series (sld_path + series_index) or the same output OME-TIFF appears in two shards;
  - the same file appears twice in the health audit;
  - one OME-TIFF has two different checksums;
  - shards disagree on the movie_list root.
Missing shards (a task failed or never ran) exit 3 unless --allow-missing is given, so an
'afterok' extraction never starts from an incomplete inventory.

Standard library only (runs on dev, on compute nodes, or on the Windows lab PC).

Usage:
  python3 merge_collect_shards.py SHARD_DIR OUT_DIR [--expected-shards N] [--allow-missing]
                                  [--hpc-raw-root ROOT]
"""

import argparse
import csv
import datetime
import json
import re
import sys
from pathlib import Path, PurePosixPath
from typing import Dict, List, Optional, Tuple

CHECKSUM_LINE_RE = re.compile(r"^([0-9a-fA-F]{64})\s+\*?(.+)$")
HEALTHY_STATUSES = {"OK_4D", "OK_NON_4D"}
RELATIVE_ROOT_PREFIX = "relative"


class MergeConflict(RuntimeError):
    """Raised when two shards report the same item (the shard split was not disjoint)."""


def shard_files(shard_dir: Path, base: str, ext: str) -> List[Tuple[int, Path]]:
    """Returns [(K, path)] for every '<base>_shardK.<ext>' in shard_dir, sorted by K."""
    found = []
    for p in shard_dir.glob(f"{base}_shard*.{ext}"):
        m = re.fullmatch(rf"{re.escape(base)}_shard(\d+)\.{re.escape(ext)}", p.name)
        if m:
            found.append((int(m.group(1)), p))
    return sorted(found)


def read_csv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Reads a CSV into (fieldnames, rows); an empty file gives ([], [])."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def write_csv_rows(path: Path, fieldnames: List[str], rows: List[Dict[str, str]], lineterminator: str = "\r\n") -> None:
    """Writes rows with the same CSV dialect as the collector (csv default) unless told otherwise."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator=lineterminator)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def merge_csv(
    shard_dir: Path, base: str, key_fn, sort_fn, what: str, lineterminator: str = "\r\n"
) -> Tuple[List[str], List[Dict[str, str]], List[int]]:
    """
    Concatenates '<base>_shardK.csv' files, refusing duplicate keys across shards.
    Returns (fieldnames, sorted rows, shard indices seen).
    """
    fieldnames: List[str] = []
    rows: List[Dict[str, str]] = []
    owner: Dict[object, int] = {}
    seen_shards = []
    for k, path in shard_files(shard_dir, base, "csv"):
        seen_shards.append(k)
        names, shard_rows = read_csv_rows(path)
        if names and not fieldnames:
            fieldnames = names
        elif names and names != fieldnames:
            raise MergeConflict(f"{path.name}: columns differ from the other shards")
        for row in shard_rows:
            key = key_fn(row)
            if key is not None:
                if key in owner:
                    raise MergeConflict(f"duplicate {what} {key} in shards {owner[key]} and {k}")
                owner[key] = k
            rows.append(row)
    rows.sort(key=sort_fn)
    return fieldnames, rows, seen_shards


def posix_parts(rel_path: str) -> Tuple[str, ...]:
    """Sort key matching sorted(Path) on POSIX: compare path components, not the raw string."""
    return tuple(PurePosixPath(rel_path).parts)


def is_true(value: str) -> bool:
    """CSV booleans are written by the collector as 'True'/'False'."""
    return str(value).strip().lower() in ("true", "1", "yes")


def manifest_entry(root: str, relpath: str) -> str:
    """Same rule as collect_4d_sld_timelapses.manifest_entry: absolute when a root is known."""
    return f"{root}/{relpath}" if root else relpath


def resolve_movie_list_root(provenances: List[dict], override: Optional[str]) -> str:
    """Uses --hpc-raw-root if given, else the (single) root the shards used; relative -> ''."""
    if override is not None:
        return override.rstrip("/")
    roots = {str(p.get("movie_list_root", "")) for p in provenances}
    roots = {"" if r.startswith(RELATIVE_ROOT_PREFIX) or not r else r for r in roots}
    if len(roots) > 1:
        raise MergeConflict(f"shards used different movie_list roots: {sorted(roots)}; pass --hpc-raw-root")
    return roots.pop().rstrip("/") if roots else ""


def merge_checksums(shard_dir: Path) -> Dict[str, str]:
    """Merges 'extracted_checksums_shardK.sha256'; one path with two digests is a conflict."""
    merged: Dict[str, str] = {}
    for k, path in shard_files(shard_dir, "extracted_checksums", "sha256"):
        with open(path, encoding="utf-8") as f:
            for line in f:
                m = CHECKSUM_LINE_RE.match(line.rstrip("\r\n"))
                if not m:
                    continue
                digest, rel = m.group(1).lower(), m.group(2)
                if rel in merged and merged[rel] != digest:
                    raise MergeConflict(f"checksum conflict for {rel} (shard {k})")
                merged[rel] = digest
    return merged


def summarize_health_by_year(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, int]]:
    """Same tally as collect_4d_sld_timelapses.summarize_health_by_year, from CSV rows."""
    by_year: Dict[str, Dict[str, int]] = {}
    for r in rows:
        yr_map = by_year.setdefault(r["cohort_year"], {})
        yr_map["total"] = yr_map.get("total", 0) + 1
        yr_map[r["status"]] = yr_map.get(r["status"], 0) + 1
    return dict(sorted(by_year.items()))


def merge_provenance(
    provenances: List[Tuple[int, dict]],
    inventory: List[Dict[str, str]],
    health: List[Dict[str, str]],
    by_year_lists: Dict[str, List[str]],
    root: str,
    expected: Optional[int],
    missing: List[int],
    out_dir: Path,
) -> dict:
    """Builds one provenance record: static fields from the first shard, metrics recomputed."""
    first = dict(provenances[0][1]) if provenances else {}
    qualifying = [r for r in inventory if is_true(r.get("is_4d_timelapse", ""))]
    summed = {}
    for key in ("total_sld_candidates_found", "total_sld_files_scanned", "total_raw_sld_gb"):
        summed[key] = sum((p.get("summary_metrics", {}) or {}).get(key, 0) or 0 for _, p in provenances)
    label_counts: Dict[str, int] = {}
    red_files = set()
    for _, p in provenances:
        mf = p.get("mito_construct_filter", {}) or {}
        for lab, n in (mf.get("label_counts", {}) or {}).items():
            label_counts[lab] = label_counts.get(lab, 0) + int(n)
        red_files.update(mf.get("red_mito_files", []) or [])
    merged = first
    merged["timestamp_iso"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    merged["output_dir"] = str(out_dir)
    merged["movie_list_root"] = root or "relative (resolve against $MITO_RAW_ROOT)"
    merged["shard"] = None
    merged["merged_from_shards"] = {
        "shard_indices": [k for k, _ in provenances],
        "expected_shard_count": expected,
        "missing_shards": missing,
        "shard_timestamps": {str(k): p.get("timestamp_iso") for k, p in provenances},
    }
    merged["summary_metrics"] = {
        "total_sld_candidates_found": int(summed["total_sld_candidates_found"]),
        "total_sld_files_scanned": int(summed["total_sld_files_scanned"]),
        "total_unparseable_or_skipped": sum(1 for r in health if r["status"] not in HEALTHY_STATUSES),
        "total_raw_sld_gb": round(float(summed["total_raw_sld_gb"]), 3),
        "total_series_inspected": len(inventory),
        "qualifying_4d_series_count": len(qualifying),
        "qualifying_sld_files_count": len({r["sld_path"] for r in qualifying}),
        "estimated_uncompressed_4d_gb": round(sum(float(r["estimated_uncompressed_gb"] or 0) for r in qualifying), 3),
        "cohorts_discovered": {yr: len(paths) for yr, paths in by_year_lists.items()},
    }
    merged["health_audit_by_year"] = summarize_health_by_year(health)
    if "mito_construct_filter" in merged:
        merged["mito_construct_filter"] = dict(merged["mito_construct_filter"])
        merged["mito_construct_filter"]["label_counts"] = dict(sorted(label_counts.items()))
        merged["mito_construct_filter"]["red_mito_files"] = sorted(red_files)
    return merged


def run_merge(
    shard_dir: Path,
    out_dir: Path,
    expected_shards: Optional[int] = None,
    allow_missing: bool = False,
    hpc_raw_root: Optional[str] = None,
) -> Tuple[int, dict]:
    """
    Performs the merge. Returns (exit_code, summary). Exit codes: 0 ok, 2 conflict, 3 missing shards.
    Outputs are written even when shards are missing (so they can be inspected), but exit 3 then.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    provenances: List[Tuple[int, dict]] = []
    for k, path in shard_files(shard_dir, "extraction_provenance", "json"):
        with open(path, encoding="utf-8") as f:
            provenances.append((k, json.load(f)))
    if expected_shards is None:
        counts = {((p.get("shard") or {}).get("count")) for _, p in provenances} - {None}
        expected_shards = max(counts) if counts else None
    present = sorted(k for k, _ in provenances)
    missing = sorted(set(range(1, expected_shards + 1)) - set(present)) if expected_shards else []

    inv_fields, inventory, _ = merge_csv(
        shard_dir, "sld_4d_inventory",
        key_fn=lambda r: (r["sld_path"], int(r["series_index"])),
        sort_fn=lambda r: (r["sld_path"], int(r["series_index"])),
        what="series",
    )
    out_paths = {}
    for r in inventory:
        rel = r.get("extracted_ometiff_relpath", "")
        if rel:
            if rel in out_paths:
                raise MergeConflict(f"two series map to the same OME-TIFF {rel}")
            out_paths[rel] = r
    health_fields, health, _ = merge_csv(
        shard_dir, "sld_health_audit",
        key_fn=lambda r: r["rel_path"],
        sort_fn=lambda r: (r["cohort_year"], r["rel_path"]),
        what="file",
    )
    if inv_fields:
        write_csv_rows(out_dir / "sld_4d_inventory.csv", inv_fields, inventory)
    if health_fields:
        write_csv_rows(out_dir / "sld_health_audit.csv", health_fields, health)

    root = resolve_movie_list_root([p for _, p in provenances], hpc_raw_root)
    by_year: Dict[str, List[str]] = {}
    with open(out_dir / "movie_list.txt", "w", encoding="utf-8") as f:
        for r in inventory:
            if is_true(r.get("is_4d_timelapse", "")):
                entry = manifest_entry(root, r["extracted_ometiff_relpath"])
                f.write(f"{entry}\n")
                by_year.setdefault(r["cohort_year"], []).append(entry)
    for yr, entries in by_year.items():
        with open(out_dir / f"movie_list_{yr}.txt", "w", encoding="utf-8") as f:
            for e in entries:
                f.write(f"{e}\n")

    cls_fields, cls_rows, _ = merge_csv(
        shard_dir, "mito_classification",
        key_fn=lambda r: r["rel_path"], sort_fn=lambda r: posix_parts(r["rel_path"]),
        what="classified file",
    )
    if cls_fields:
        write_csv_rows(out_dir / "mito_classification.csv", cls_fields, cls_rows, lineterminator="\n")
        review = [r for r in cls_rows if r.get("label") == "REVIEW"]
        review.sort(key=lambda r: (r.get("review_priority") != "high", r["rel_path"].lower()))
        write_csv_rows(out_dir / "mito_review.csv", cls_fields, review, lineterminator="\n")

    checksums = merge_checksums(shard_dir)
    if checksums:
        with open(out_dir / "extracted_checksums.sha256", "w", encoding="utf-8") as f:
            for rel, digest in sorted(checksums.items()):
                f.write(f"{digest}  {rel}\n")

    prov = merge_provenance(provenances, inventory, health, by_year, root, expected_shards, missing, out_dir)
    with open(out_dir / "extraction_provenance.json", "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2)

    summary = {
        "shards_present": present,
        "shards_missing": missing,
        "series": len(inventory),
        "qualifying_series": sum(len(v) for v in by_year.values()),
        "files_audited": len(health),
        "checksums": len(checksums),
    }
    if missing and not allow_missing:
        return 3, summary
    return 0, summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge sharded collect_4d_sld_timelapses.py reports into one set.")
    parser.add_argument("shard_dir", type=Path, help="Directory holding the *_shardK.* reports")
    parser.add_argument("out_dir", type=Path, help="Where the merged reports are written (may equal shard_dir)")
    parser.add_argument("--expected-shards", type=int, default=None,
                        help="Number of shards that should exist (default: read from the shards' provenance)")
    parser.add_argument("--allow-missing", action="store_true",
                        help="Exit 0 even when some shards are missing (their files are then absent from the merge)")
    parser.add_argument("--hpc-raw-root", type=str, default=None,
                        help="Prefix for movie_list.txt entries (default: the root the shards used; '' = relative)")
    args = parser.parse_args()

    if not args.shard_dir.is_dir():
        print(f"[ERROR] Shard directory not found: {args.shard_dir}", file=sys.stderr)
        return 1
    try:
        rc, summary = run_merge(
            args.shard_dir.resolve(), args.out_dir.resolve(), args.expected_shards,
            args.allow_missing, args.hpc_raw_root,
        )
    except MergeConflict as exc:
        print(f"[ERROR] Merge refused: {exc}", file=sys.stderr)
        return 2
    print(f"Merged shards {summary['shards_present']} -> {args.out_dir}")
    print(f"  series: {summary['series']} (qualifying 4D: {summary['qualifying_series']}), "
          f"files audited: {summary['files_audited']}, checksums: {summary['checksums']}")
    if summary["shards_missing"]:
        level = "WARN" if args.allow_missing else "ERROR"
        print(f"[{level}] Missing shard(s): {summary['shards_missing']} (failed or never ran). "
              "Re-run them with 'submit_full_run.bash retry', or pass --allow-missing.", file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
