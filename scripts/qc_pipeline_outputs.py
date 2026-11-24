#!/usr/bin/env python3
"""
qc_pipeline_outputs.py
Per-stage quality-control report for one movie's Stage 3-5 outputs (read-only).

Prints, for each stage, the numbers needed to judge whether the stage did what it should:
  Stage 3  : per-frame mito_norm signal inside/outside the cell mask, drift.
  Stage 4A : per-frame ilastik foreground fraction inside/outside the cell mask.
  Stage 4B : per-frame instance counts, object volume / Z-extent distributions, how many
             objects sit on zero-intensity (masked-out) voxels, foreground coverage.
  Stage 5A : spot volume/intensity, link displacement vs the configured limits, track lengths,
             split/merge counts per frame.
  Stage 5B : geodesic distance and compartment distributions, morphology classes.

Usage (on a compute node; it loads ~1.5 GB of stacks):
  python3 scripts/qc_pipeline_outputs.py --work-dir <scratch movie dir> --movie-id <id> [--stages 3,4,5]
"""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import ndimage

try:
    from tiff_io import read_tiff
except ImportError:
    from .tiff_io import read_tiff


def pct(values, qs=(5, 50, 95, 99)):
    """Compact percentile summary string for a 1D array."""
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return "n=0"
    parts = [f"p{q}={np.percentile(values, q):.3g}" for q in qs]
    return f"n={values.size} min={values.min():.3g} " + " ".join(parts) + f" max={values.max():.3g}"


def read_csv_rows(path: Path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def qc_stage3(work: Path, mid: str, mask: np.ndarray):
    print("\n=== STAGE 3: preprocessing ===")
    norm, _ = read_tiff(work / f"{mid}_mito_norm.tif", axes="TZYX")
    print(f"mito_norm shape (T,Z,Y,X) = {norm.shape}, dtype={norm.dtype}")
    print(f"cell mask: {mask.mean() * 100:.1f}% of voxels inside")
    print("frame | nonzero% in-mask | nonzero% outside | p50 / p99 in-mask (nonzero)")
    for t in range(norm.shape[0]):
        inside = norm[t][mask]
        outside = norm[t][~mask]
        nz = inside[inside > 0]
        p50, p99 = (np.percentile(nz, [50, 99]) if nz.size else (0, 0))
        print(
            f"  {t:3d} | {100 * (inside > 0).mean():6.2f} | {100 * (outside > 0).mean():6.2f} "
            f"| {p50:8.0f} / {p99:8.0f}"
        )
    drift_csv = work / f"{mid}_drift_shifts.csv"
    if drift_csv.exists():
        tot = [float(r["total_drift_um"]) for r in read_csv_rows(drift_csv)]
        print(f"drift total_drift_um per frame: {[round(v, 2) for v in tot]}")
    return norm


def qc_stage4(work: Path, mid: str, mask: np.ndarray, norm: np.ndarray, calib: dict):
    print("\n=== STAGE 4A: ilastik probabilities ===")
    prob, _ = read_tiff(work / f"{mid}_Probabilities.tif", axes="TZYX")
    print(f"prob shape = {prob.shape}, range [{float(prob.min()):.3f}, {float(prob.max()):.3f}]")
    print("frame | P>0.5 % in-mask | P>0.5 % outside | P>0.5 % where mito_norm==0")
    for t in range(prob.shape[0]):
        fg = prob[t] > 0.5
        zero = norm[t] == 0
        print(
            f"  {t:3d} | {100 * fg[mask].mean():6.2f} | {100 * fg[~mask].mean():6.2f} "
            f"| {100 * fg[zero].mean():6.2f}"
        )

    print("\n=== STAGE 4B: instance labels ===")
    labels, _ = read_tiff(work / f"{mid}_mito_labels.tif", axes="TZYX")
    vox_um3 = calib["dx_um"] * calib["dy_um"] * calib["dz_um"]
    print(f"voxel volume = {vox_um3:.4f} um^3")
    print("frame | objects | median vol um3 | p95 vol | max vol | 1-slice % | zero-intensity objs % | fg coverage %")
    all_vols, all_zext = [], []
    for t in range(labels.shape[0]):
        lab = labels[t]
        n = int(lab.max())
        if n == 0:
            print(f"  {t:3d} | 0")
            continue
        counts = np.bincount(lab.ravel(), minlength=n + 1)[1:]
        present = counts > 0
        vols = counts[present] * vox_um3
        sums = ndimage.sum_labels(norm[t].astype(np.float64), lab, index=np.arange(1, n + 1))[present]
        objs = ndimage.find_objects(lab)
        zext = np.array([sl[0].stop - sl[0].start for sl in objs if sl is not None])
        fg = prob[t] > 0.5
        coverage = 100 * (lab[fg] > 0).mean() if fg.any() else 0.0
        all_vols.append(vols)
        all_zext.append(zext)
        print(
            f"  {t:3d} | {present.sum():6d} | {np.median(vols):8.3f} | {np.percentile(vols, 95):8.2f} "
            f"| {vols.max():8.1f} | {100 * (zext == 1).mean():5.1f} | {100 * (sums == 0).mean():5.1f} "
            f"| {coverage:5.1f}"
        )
    if all_vols:
        v = np.concatenate(all_vols)
        z = np.concatenate(all_zext)
        print(f"all-frame volume um3: {pct(v)}")
        print(f"all-frame Z extent (slices): {pct(z)}")
        print(f"objects > 5 um3: {(v > 5).sum()} ({100 * (v > 5).mean():.1f}%)")
    return labels


def qc_stage5(work: Path, mid: str, cfg_track: dict, calib: dict):
    print("\n=== STAGE 5A: tracking ===")
    spots = read_csv_rows(work / f"{mid}_spots.csv")
    edges = read_csv_rows(work / f"{mid}_edges.csv")
    tracks = read_csv_rows(work / f"{mid}_tracks.csv")
    vol = np.array([float(s["VOLUME"]) for s in spots])
    ch2 = np.array([float(s["MEAN_INTENSITY_CH2"]) for s in spots])
    print(f"spots: {len(spots)}; VOLUME um3: {pct(vol)}")
    print(f"MEAN_INTENSITY_CH2: {pct(ch2)}; zero-intensity spots: {100 * (ch2 == 0).mean():.1f}%")
    in_track = sum(1 for s in spots if s.get("TRACK_ID", "") not in ("", "-1", "None"))
    print(f"spots assigned to a kept track: {in_track} ({100 * in_track / max(1, len(spots)):.1f}%)")

    by_type = defaultdict(list)
    for e in edges:
        by_type[e.get("EVENT_TYPE", "?")].append(float(e["DISPLACEMENT"]))
    limits = {
        "FRAME_LINK": cfg_track.get("linking_max_distance_um"),
        "GAP_CLOSE": cfg_track.get("gap_closing_max_distance_um"),
        "SPLIT_FISSION": cfg_track.get("splitting_max_distance_um"),
        "MERGE_FUSION": cfg_track.get("merging_max_distance_um"),
    }
    for etype, d in sorted(by_type.items()):
        print(f"edges {etype:14s}: {pct(d)} (config limit {limits.get(etype)} um)")
    spot_frame = {s["ID"]: int(s["FRAME"]) for s in spots}
    ev_per_frame = defaultdict(Counter)
    for e in edges:
        if e.get("EVENT_TYPE") not in ("FRAME_LINK", None):
            ev_per_frame[spot_frame.get(e["SPOT_SOURCE_ID"], -1)][e["EVENT_TYPE"]] += 1
    print("non-link events by source frame: " + json.dumps({k: dict(v) for k, v in sorted(ev_per_frame.items())}))
    nsp = np.array([int(t["NUMBER_SPOTS"]) for t in tracks])
    print(f"tracks: {len(tracks)}; NUMBER_SPOTS: {pct(nsp, (50, 95))}")
    print(f"tracks with splits: {sum(int(t['NUMBER_SPLITS']) > 0 for t in tracks)}, "
          f"with merges: {sum(int(t['NUMBER_MERGES']) > 0 for t in tracks)}")
    speed = np.array([float(t["TRACK_MEAN_SPEED"]) for t in tracks])
    print(f"TRACK_MEAN_SPEED um/s: {pct(speed)}")

    print("\n=== STAGE 5B: arbor / MitoLandscape ===")
    ml_path = work / f"{mid}_snt_mitolandscape_tracks.csv"
    if not ml_path.exists():
        print("  (no snt_mitolandscape_tracks.csv)")
        return
    ml = read_csv_rows(ml_path)
    geo = np.array([float(r["geodesic_soma_dist_um"]) for r in ml])
    print(f"geodesic_soma_dist_um: {pct(geo)}; exactly 0: {100 * (geo == 0).mean():.1f}%")
    comp_geo = defaultdict(list)
    for r, g in zip(ml, geo):
        comp_geo[r["compartment_name"]].append(g)
    for c, g in sorted(comp_geo.items()):
        g = np.asarray(g)
        print(f"  {c:15s}: n={g.size:6d}, geo==0 {100 * (g == 0).mean():5.1f}%, median geo {np.median(g):.2f} um")
    print("morphology_class: " + json.dumps(Counter(r["morphology_class"] for r in ml)))
    motility_col = "transport_state" if "transport_state" in ml[0] else None
    if motility_col:
        print(f"motility ({motility_col}): " + json.dumps(Counter(r[motility_col] for r in ml)))


def main() -> int:
    ap = argparse.ArgumentParser(description="Per-stage QC report for one movie (read-only)")
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--movie-id", required=True)
    ap.add_argument("--config", type=Path, default=Path(__file__).resolve().parent.parent / "pipeline_config.yaml")
    ap.add_argument("--stages", default="3,4,5")
    args = ap.parse_args()
    stages = {s.strip() for s in args.stages.split(",")}
    work, mid = args.work_dir, args.movie_id

    with open(work / f"{mid}_provenance.json", encoding="utf-8") as f:
        prov = json.load(f)
    calib = prov["calibration"]
    print(f"Movie: {mid}\ncalibration: {calib}")
    print(f"stage4a source: {prov.get('stage4_segmentation_metrics', {}).get('stage4a_probability_source')}")

    from pipeline_config_io import load_config
    cfg = load_config(args.config)
    cfg_track = cfg.get("stage5_tracking_and_mitolandscape", {}).get("sparse_lap_tracker", {})

    mask_raw, _ = read_tiff(work / f"{mid}_cell_mask.tif")
    mask = np.asarray(mask_raw).reshape(mask_raw.shape[-3:]) > 0  # static 3D cell mask (Z, Y, X)

    norm = None
    if "3" in stages or "4" in stages:
        norm = qc_stage3(work, mid, mask)
    if "4" in stages:
        qc_stage4(work, mid, mask, norm, calib)
    if "5" in stages:
        qc_stage5(work, mid, cfg_track, calib)
    return 0


if __name__ == "__main__":
    sys.exit(main())
