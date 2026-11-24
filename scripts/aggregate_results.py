#!/usr/bin/env python3
"""
aggregate_results.py
Phase G of a full-archive run: combine every movie's deliverables into cohort-level tables.

Inputs
  --run-dir   the run directory of submit_full_run.bash (its run_manifest.json lists every movie
              list that was submitted, and the sacct job ids)
  Per-movie deliverables in $MITO_OUTPUT_ROOT/<cohort>/<movie_id>/ (copied by run_pipeline_slurm.bash)
  and the stage-done markers in $MITO_SCRATCH_ROOT/<cohort>/<movie_id>/.

Outputs (in --out-dir, default <run-dir>/results; also copied to $MITO_OUTPUT_ROOT/_runs/<run_id>/)
  all_spots.csv, all_spots_profiled.csv, all_tracks.csv, all_edges.csv, all_remodeling_events.csv,
  all_process_profiles.csv
      the per-movie CSVs concatenated, each row prefixed with movie_id, cohort, protocol, mito_label
      (protocol = acquisition calibration key 'dt<s>_dz<um>_dxy<um>', so movies recorded with the
      same settings can be grouped).
  qc_dashboard.csv
      one row per movie: Stage 3 QC status, raw mito contrast (sigma), drift, SNR, objects per
      frame, track / fission / fusion counts, stage-done marker times, sacct runtimes and peak memory.
  failures.csv
      every movie that did not finish Stage 5, with the reason (QC fail, task state, missing output).
  summary.json
      counts per status, per cohort and per protocol, plus the manifest provenance.

Standard library only.
"""

import argparse
import csv
import datetime
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import plan_full_run as plan  # noqa: E402  (stdlib-only helpers: manifest, sacct parsing, lists)

# (output table, per-movie file suffix)
TABLES = [
    ("all_spots.csv", "_spots.csv"),
    ("all_spots_profiled.csv", "_snt_mitolandscape_tracks.csv"),
    ("all_tracks.csv", "_tracks.csv"),
    ("all_edges.csv", "_edges.csv"),
    ("all_remodeling_events.csv", "_remodeling_events.csv"),
    ("all_process_profiles.csv", "_process_0.5um_profile.csv"),
]
PREFIX_COLUMNS = ["movie_id", "cohort", "protocol", "mito_label"]
PIPELINE_PHASES = ("canary", "stage3", "stage45", "overflow")


def mito_env() -> Dict[str, str]:
    """MITO_* storage paths: environment first, then the documented defaults (pipeline_config_io)."""
    user = os.environ.get("USER", "")
    data_root = os.environ.get("MITO_DATA_ROOT") or f"/sciclone/data10/{user}/mito_pipeline"
    return {
        "output": os.environ.get("MITO_OUTPUT_ROOT") or f"{data_root}/outputs",
        "scratch": os.environ.get("MITO_SCRATCH_ROOT") or f"/sciclone/scr30/{user}/mito_work",
        "raw": os.environ.get("MITO_RAW_ROOT") or f"{data_root}/raw_ometiff",
    }


def protocol_key(calib: Dict[str, object]) -> str:
    """'dt120.0_dz1.00_dxy0.188' from a calibration dict; '' when unknown."""
    try:
        dt = calib.get("dt_s")
        dz = calib.get("dz_um")
        dx = calib.get("dx_um")
        parts = []
        if dt not in (None, ""):
            parts.append(f"dt{float(dt):.1f}")
        if dz not in (None, ""):
            parts.append(f"dz{float(dz):.2f}")
        if dx not in (None, ""):
            parts.append(f"dxy{float(dx):.3f}")
        return "_".join(parts)
    except (TypeError, ValueError):
        return ""


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def marker_time(path: Path) -> str:
    """'completed=' value of a .done marker, or ''."""
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith("completed="):
                return line.split("=", 1)[1]
    except OSError:
        pass
    return ""


def movie_lists(man: plan.Manifest) -> List[Tuple[str, Dict[str, object]]]:
    """(job name, job record) for every pipeline array that has a movie list."""
    return [(n, j) for n, j in sorted(man.data["jobs"].items())
            if j.get("phase") in PIPELINE_PHASES and j.get("list")]


def sacct_by_movie(man: plan.Manifest, sacct_text: Optional[str]) -> Dict[Tuple[str, str], Dict[str, object]]:
    """
    {(movie_path, 'stage3'|'stage45'): latest sacct record} by mapping array index -> list line.
    One sacct call for the whole run (or a provided text, for tests).
    """
    jobs = movie_lists(man)
    ids = []
    for _, j in jobs:
        ids += [j.get("job_id", "")] + [r.get("job_id", "") for r in j.get("retries", [])]
    ids = [i for i in ids if i and not i.startswith("DRY_")]
    if sacct_text is None:
        sacct_text = ""
        if ids:
            try:
                proc = subprocess.run(
                    ["sacct", "-P", "-n", "--units=M", "-o", "JobID,State,ExitCode,Elapsed,MaxRSS", "-j", ",".join(ids)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120, check=False)
                sacct_text = proc.stdout if proc.returncode == 0 else ""
            except (OSError, subprocess.SubprocessError):
                sacct_text = ""
    out: Dict[Tuple[str, str], Dict[str, object]] = {}
    for _, j in jobs:
        try:
            entries = plan.read_list(Path(str(j["list"])))
        except OSError:
            continue
        kind = "stage45" if j.get("phase") == "stage45" or "4" in str(j.get("stages", "")) else "stage3"
        states = plan.latest_task_states(j, sacct_text)
        for idx, rec in states.items():
            if 1 <= idx <= len(entries):
                out[(entries[idx - 1], kind)] = rec
    return out


def collect_movies(man: plan.Manifest, raw_root: str) -> List[str]:
    """Unique movie paths across every submitted list (canary movies also appear in Stage 3 lists)."""
    seen: Dict[str, None] = {}
    for _, j in movie_lists(man):
        try:
            for e in plan.read_list(Path(str(j["list"]))):
                path = e if os.path.isabs(e) else os.path.join(raw_root, e)
                seen.setdefault(path, None)
        except OSError:
            continue
    return list(seen)


def movie_id_of(path: str) -> str:
    name = Path(path).name
    for suffix in (".ome.tif", ".ome.tiff", ".tif", ".tiff"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def concat_tables(movies: List[Dict[str, str]], out_dir: Path) -> Dict[str, int]:
    """Streams each per-movie CSV into one table (column union, prefix columns first)."""
    counts = {}
    for table, suffix in TABLES:
        sources = []
        columns: List[str] = []
        for mv in movies:
            path = Path(mv["final_dir"]) / f"{mv['movie_id']}{suffix}"
            if not path.exists():
                continue
            with open(path, newline="", encoding="utf-8") as f:
                header = next(csv.reader(f), [])
            for c in header:
                if c not in columns and c not in PREFIX_COLUMNS:
                    columns.append(c)
            sources.append((mv, path))
        n = 0
        with open(out_dir / table, "w", newline="", encoding="utf-8") as fo:
            writer = csv.DictWriter(fo, fieldnames=PREFIX_COLUMNS + columns, extrasaction="ignore")
            writer.writeheader()
            for mv, path in sources:
                with open(path, newline="", encoding="utf-8") as f:
                    for row in csv.DictReader(f):
                        row.update({k: mv[k] for k in PREFIX_COLUMNS})
                        writer.writerow(row)
                        n += 1
        counts[table] = n
    return counts


def build_dashboard_row(path: str, env: Dict[str, str], inv_by_rel: Dict[str, Dict[str, str]],
                        sacct: Dict[Tuple[str, str], Dict[str, object]]) -> Dict[str, str]:
    """One qc_dashboard.csv row for one movie."""
    mid = movie_id_of(path)
    cohort = Path(path).parent.name
    final_dir = Path(env["output"]) / cohort / mid
    work_dir = Path(env["scratch"]) / cohort / mid
    prov = read_json(final_dir / f"{mid}_provenance.json") or read_json(work_dir / f"{mid}_provenance.json")
    inv = next((r for rel, r in inv_by_rel.items() if path.endswith("/" + rel) or path == rel), {})
    calib = prov.get("calibration") or {
        "dt_s": inv.get("time_interval_s"), "dz_um": inv.get("phys_size_z_um"), "dx_um": inv.get("phys_size_x_um")}
    s3 = prov.get("stage3_preprocessing_metrics", {}) or {}
    s4 = prov.get("stage4_segmentation_metrics", {}) or {}
    s5 = prov.get("stage5_tracking_summary", {}) or {}
    contrast = (s3.get("raw_mito_contrast") or {}).get("contrast_sigma", "")
    row = {
        "movie_id": mid, "cohort": cohort, "protocol": protocol_key(calib),
        "mito_label": inv.get("mito_label", ""), "movie_path": path,
        "estimated_uncompressed_gb": inv.get("estimated_uncompressed_gb", ""),
        "size_t": inv.get("size_t", ""),
        "stage3_qc_status": s3.get("qc_status", ""),
        "raw_mito_contrast_sigma": contrast,
        "max_cumulative_drift_um": s3.get("max_cumulative_drift_um", ""),
        "estimated_snr": s3.get("estimated_snr", ""),
        "stage4a_probability_source": s4.get("stage4a_probability_source", ""),
        "instance_engine": s4.get("instance_engine", ""),
        "mean_objects_per_frame": s4.get("mean_instances_per_frame", ""),
        "total_tracks": s5.get("total_tracks_filtered", ""),
        "fission_events": s5.get("total_fission_events", ""),
        "fusion_events": s5.get("total_fusion_events", ""),
        "stage3_done": marker_time(work_dir / f"{mid}_stage3.done"),
        "stage4_done": marker_time(work_dir / f"{mid}_stage4.done"),
        "stage5_done": marker_time(work_dir / f"{mid}_stage5.done"),
        "final_dir": str(final_dir),
    }
    for kind in ("stage3", "stage45"):
        rec = sacct.get((path, kind)) or {}
        row[f"{kind}_state"] = str(rec.get("state", ""))
        row[f"{kind}_exit"] = str(rec.get("exit", "")) if rec else ""
        row[f"{kind}_elapsed_s"] = f"{float(rec['elapsed_s']):.0f}" if rec.get("elapsed_s") else ""
        row[f"{kind}_maxrss_gb"] = f"{float(rec['maxrss_mb']) / 1024:.2f}" if rec.get("maxrss_mb") else ""
    return row


def failure_reason(row: Dict[str, str]) -> Optional[str]:
    """None when Stage 5 finished; otherwise the most specific reason available."""
    if row["stage5_done"] and Path(row["final_dir"], f"{row['movie_id']}_tracks.csv").exists():
        return None
    if row["stage3_qc_status"].startswith("FAIL") or row["stage3_exit"] == "3" or row["stage45_exit"] == "3":
        return f"QC_FAIL ({row['stage3_qc_status'] or 'exit 3'}; contrast {row['raw_mito_contrast_sigma'] or '?'} sigma)"
    if row["stage45_exit"] == "5":
        return "STAGE4_IMPLAUSIBLE_FOREGROUND (foreground covers most of the cell; retrain ilastik, exit 5)"
    for kind in ("stage45", "stage3"):
        st = row[f"{kind}_state"]
        if st and st != "COMPLETED":
            return f"{kind.upper()}_{st} (exit {row[f'{kind}_exit']})"
    if not row["stage3_done"]:
        return "STAGE3_NOT_DONE"
    if not row["stage4_done"]:
        return "STAGE4_NOT_DONE"
    return "STAGE5_NOT_DONE_OR_OUTPUT_MISSING"


def run_aggregate(run_dir: Path, out_dir: Path, sacct_text: Optional[str] = None,
                  persist: bool = True) -> Dict[str, object]:
    man = plan.Manifest(run_dir / "run_manifest.json")
    if not man.data["jobs"]:
        raise SystemExit(f"[ERROR] No executed jobs in {run_dir}/run_manifest.json")
    env = mito_env()
    raw_root = str(man.data.get("raw_root") or env["raw"])
    inventory = plan.read_inventory(run_dir / "scan" / "merged" / "sld_4d_inventory.csv")
    inv_by_rel = {r["extracted_ometiff_relpath"]: r for r in inventory if r.get("extracted_ometiff_relpath")}
    sacct = sacct_by_movie(man, sacct_text)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = [build_dashboard_row(p, env, inv_by_rel, sacct) for p in collect_movies(man, raw_root)]
    with open(out_dir / "qc_dashboard.csv", "w", newline="", encoding="utf-8") as f:
        fields = list(rows[0].keys()) if rows else ["movie_id"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    failures = []
    for r in rows:
        reason = failure_reason(r)
        if reason:
            failures.append({"movie_id": r["movie_id"], "cohort": r["cohort"], "reason": reason,
                             "stage3_qc_status": r["stage3_qc_status"], "movie_path": r["movie_path"],
                             "logs": str(man.data.get("log_dir", ""))})
    with open(out_dir / "failures.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["movie_id", "cohort", "reason", "stage3_qc_status", "movie_path", "logs"])
        writer.writeheader()
        writer.writerows(failures)

    complete = [r for r in rows if failure_reason(r) is None]
    counts = concat_tables([{"movie_id": r["movie_id"], "cohort": r["cohort"], "protocol": r["protocol"],
                             "mito_label": r["mito_label"], "final_dir": r["final_dir"]} for r in complete], out_dir)

    def tally(key: str, subset: List[Dict[str, str]]) -> Dict[str, int]:
        t: Dict[str, int] = {}
        for r in subset:
            t[r[key] or "unknown"] = t.get(r[key] or "unknown", 0) + 1
        return dict(sorted(t.items()))

    summary = {
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "run_id": man.data.get("run_id"),
        "git_sha": man.data.get("git_sha"),
        "config_sha256": man.data.get("config_sha256"),
        "movies_total": len(rows),
        "movies_complete": len(complete),
        "movies_failed": len(failures),
        "failure_reasons": tally("reason", failures) if failures else {},
        "complete_by_cohort": tally("cohort", complete),
        "complete_by_protocol": tally("protocol", complete),
        "rows_per_table": counts,
        "model_refit": man.data.get("model_refit"),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    shutil.copy2(man.path, out_dir / "run_manifest.json")

    if persist:
        # scr30 may be purged after ~30 days: keep the small cohort tables on project storage too.
        dest = Path(env["output"]) / "_runs" / str(man.data.get("run_id") or run_dir.name)
        try:
            dest.mkdir(parents=True, exist_ok=True)
            for p in out_dir.iterdir():
                if p.is_file():
                    shutil.copy2(p, dest / p.name)
            summary["persisted_to"] = str(dest)
        except OSError as exc:
            print(f"[WARN] could not copy results to {dest}: {exc}", file=sys.stderr)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate per-movie pipeline outputs of a full run.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", type=Path, default=None, help="default: <run-dir>/results")
    parser.add_argument("--sacct-file", type=Path, default=None, help="use saved sacct output instead of calling sacct")
    parser.add_argument("--no-persist", action="store_true", help="do not copy results to $MITO_OUTPUT_ROOT/_runs/")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    out_dir = (args.out_dir or run_dir / "results").resolve()
    sacct_text = args.sacct_file.read_text() if args.sacct_file else None
    s = run_aggregate(run_dir, out_dir, sacct_text, persist=not args.no_persist)
    print(f"Aggregated {s['movies_total']} movie(s): {s['movies_complete']} complete, {s['movies_failed']} failed "
          f"-> {out_dir}")
    for reason, n in (s.get("failure_reasons") or {}).items():
        print(f"  {n:4d} x {reason}")
    if s.get("persisted_to"):
        print(f"Copied to project storage: {s['persisted_to']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
