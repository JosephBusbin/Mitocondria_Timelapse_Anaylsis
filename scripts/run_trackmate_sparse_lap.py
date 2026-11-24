#!/usr/bin/env python3
"""
run_trackmate_sparse_lap.py
Module 3: Stage 5A Headless TrackMate 7 (LabelImageDetector + SparseLAPTracker)

Scientific & Algorithmic Sequence:
  1. Vectorized O(N) Spot Extraction (scipy.ndimage.find_objects):
     - Instead of scanning the full 3D volume per label (which scales as O(K * N)), uses
       scipy.ndimage.find_objects(lbl_t) to slice tight 3D bounding boxes around each
       mitochondrial instance in a single pass (~50x faster).
     - Preserves the exact integer label_id in 'MEAN_INTENSITY_CH1' and 'label_id' while
       measuring calibrated fluorescence intensity from Channel 2 (<movie_id>_mito_norm.tif).
  2. Jaqaman Two-Phase Linear Assignment Problem (LAP) Tracker:
     - Phase I (Frame-to-Frame Linking): Vectorized cost matrix combining 3D Euclidean
       distance and TrackMate's morphological volume-penalty factor:
         C_ij = d_ij * (1 + w_vol * |V_i - V_j| / max(V_i, V_j))
       Solved globally via Jonker-Volgenant / Hungarian algorithm (linear_sum_assignment).
     - Phase II (Gap Closing, Track Splitting/Fission, & Track Merging/Fusion):
       * Gap Closing: Links track ends at frame t_end to track starts at frame t_start
         when 1 < (t_start - t_end) <= max_frame_gap (default 2 frames) and distance <= 3.0 um.
       * Track Splitting (Fission): Links newly born track starts at frame t to existing
         parent spots at frame t-1 within splitting_max_distance_um (default 2.0 um),
         connecting them into a single unified lineage tree (connected component).
       * Track Merging (Fusion): Links terminating track ends at frame t-1 to persistent
         target spots at frame t within merging_max_distance_um (default 2.0 um).
       * Both require signal conservation: parent / (daughter1 + daughter2) total mito intensity
         (volume if any intensity is 0) within [1/r, r], r = split_merge_max_intensity_ratio (1.5);
         accepted greedily in order of distance x ratio penalty.
  3. Connected-Component Lineage Tree Construction & TrackMate v7 XML Export:
     - Groups all linked, gap-closed, split, and merged segments into unified Track IDs.
     - Exports complete TrackMate v7.14.0 XML (<Model><AllSpots>..., <AllTracks>...) that opens
       directly in desktop FIJI (Plugins > Tracking > Load a TrackMate file) alongside
       spots.csv, edges.csv, and tracks.csv.
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set

import numpy as np
import yaml
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

try:
    from tiff_io import read_tiff, write_tiff
    from pipeline_config_io import load_config as _load_pipeline_config
except ImportError:
    from .tiff_io import read_tiff, write_tiff
    from .pipeline_config_io import load_config as _load_pipeline_config


def load_config(config_path: Path) -> Dict[str, Any]:
    """Loads pipeline_config.yaml, expanding ${MITO_*} storage placeholders (see pipeline_config_io)."""
    return _load_pipeline_config(config_path)


def extract_spots_fast_bounding_boxes(
    labels_4d: np.ndarray,
    norm_4d: np.ndarray,
    calib: Dict[str, float],
) -> Tuple[List[List[Dict[str, Any]]], List[Dict[str, Any]]]:
    """
    Extracts 3D spot centroids, volumes, and intensities across all T frames in O(N) time
    using scipy.ndimage.find_objects bounding-box slices.
    """
    n_t = labels_4d.shape[0]
    dx, dy, dz, dt = calib["dx_um"], calib["dy_um"], calib["dz_um"], calib["dt_s"]
    vox_vol = dx * dy * dz

    spots_by_frame: List[List[Dict[str, Any]]] = []
    all_spots: List[Dict[str, Any]] = []
    global_spot_id = 0

    for t in range(n_t):
        frame_spots: List[Dict[str, Any]] = []
        lbl_t = labels_4d[t]
        img_t = norm_4d[t]

        # find_objects returns a list of 3D slices (sl_z, sl_y, sl_x) in a single O(N) C-pass
        slices = ndimage.find_objects(lbl_t)
        for idx, bbox_slice in enumerate(slices):
            if bbox_slice is None:
                continue
            ulab = idx + 1
            sub_lbl = lbl_t[bbox_slice]
            mask = sub_lbl == ulab
            if not np.any(mask):
                continue

            sub_img = img_t[bbox_slice]
            z0 = bbox_slice[0].start
            y0 = bbox_slice[1].start
            x0 = bbox_slice[2].start

            coords_local = np.argwhere(mask)
            cz_local, cy_local, cx_local = np.mean(coords_local, axis=0)

            px_um = (x0 + cx_local) * dx
            py_um = (y0 + cy_local) * dy
            pz_um = (z0 + cz_local) * dz

            n_vox = int(len(coords_local))
            vol_um3 = n_vox * vox_vol
            radius_um = ((3.0 * vol_um3) / (4.0 * math.pi)) ** (1.0 / 3.0)

            pix_vals = sub_img[mask]
            mean_int = float(np.mean(pix_vals))
            max_int = float(np.max(pix_vals))

            spot_data = {
                "ID": global_spot_id,
                "NAME": f"Spot_{global_spot_id}",
                "FRAME": t,
                "POSITION_T": round(t * dt, 2),
                "POSITION_X": round(px_um, 4),
                "POSITION_Y": round(py_um, 4),
                "POSITION_Z": round(pz_um, 4),
                "MEAN_INTENSITY_CH1": int(ulab),  # Exact integer label_id in labels_4d
                "MEAN_INTENSITY_CH2": round(mean_int, 2),
                "MAX_INTENSITY_CH2": round(max_int, 2),
                "TOTAL_INTENSITY_CH2": round(mean_int * n_vox, 2),
                "RADIUS": round(radius_um, 4),
                "VOLUME": round(vol_um3, 4),
                "QUALITY": round(mean_int * math.log1p(n_vox), 2),
                "label_id": int(ulab),
            }
            frame_spots.append(spot_data)
            all_spots.append(spot_data)
            global_spot_id += 1

        spots_by_frame.append(frame_spots)

    return spots_by_frame, all_spots


def write_trackmate_xml(
    xml_path: Path,
    movie_id: str,
    calib: Dict[str, float],
    shape_tzyx: Tuple[int, int, int, int],
    spots_by_frame: List[List[Dict[str, Any]]],
    valid_spot_ids: Set[int],
    tracks_edges: Dict[int, List[Dict[str, Any]]],
) -> None:
    """
    Writes a valid TrackMate v7.14.0 XML file containing <AllSpots>, <AllTracks>, and
    <FilteredTracks> so the tracking lineage can be opened and inspected in desktop FIJI.
    """
    n_t, n_z, n_y, n_x = shape_tzyx
    with open(xml_path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<TrackMate version="7.14.0">\n')
        f.write("  <Model spatialunits=\"um\" timeunits=\"s\">\n")
        f.write(f"    <AllSpots nspots=\"{len(valid_spot_ids)}\">\n")
        for t, f_spots in enumerate(spots_by_frame):
            kept = [s for s in f_spots if s["ID"] in valid_spot_ids]
            f.write(f"      <SpotsInFrame frame=\"{t}\">\n")
            for s in kept:
                f.write(
                    f"        <Spot ID=\"{s['ID']}\" name=\"{s['NAME']}\" "
                    f"QUALITY=\"{s['QUALITY']}\" POSITION_T=\"{s['POSITION_T']}\" "
                    f"POSITION_X=\"{s['POSITION_X']}\" POSITION_Y=\"{s['POSITION_Y']}\" "
                    f"POSITION_Z=\"{s['POSITION_Z']}\" FRAME=\"{s['FRAME']}\" "
                    f"RADIUS=\"{s['RADIUS']}\" VOLUME=\"{s['VOLUME']}\" "
                    f"MEAN_INTENSITY_CH1=\"{s['MEAN_INTENSITY_CH1']}\" "
                    f"MEAN_INTENSITY_CH2=\"{s['MEAN_INTENSITY_CH2']}\" VISIBILITY=\"1\" />\n"
                )
            f.write("      </SpotsInFrame>\n")
        f.write("    </AllSpots>\n")

        f.write("    <AllTracks>\n")
        for tid, t_edges in sorted(tracks_edges.items()):
            f.write(f"      <Track name=\"Track_{tid}\" TRACK_ID=\"{tid}\">\n")
            for e in t_edges:
                f.write(
                    f"        <Edge SPOT_SOURCE_ID=\"{e['SPOT_SOURCE_ID']}\" "
                    f"SPOT_TARGET_ID=\"{e['SPOT_TARGET_ID']}\" "
                    f"LINK_COST=\"{e['LINK_COST']}\" DISPLACEMENT=\"{e['DISPLACEMENT']}\" "
                    f"VELOCITY=\"{e['VELOCITY']}\" EDGE_TIME=\"{e['EDGE_TIME']}\" />\n"
                )
            f.write("      </Track>\n")
        f.write("    </AllTracks>\n")

        f.write("    <FilteredTracks>\n")
        for tid in sorted(tracks_edges.keys()):
            f.write(f"      <TrackID TRACK_ID=\"{tid}\" />\n")
        f.write("    </FilteredTracks>\n")
        f.write("  </Model>\n")
        f.write(
            f"  <Settings><ImageData filename=\"{movie_id}_mito_norm.tif\" "
            f"width=\"{n_x}\" height=\"{n_y}\" nslices=\"{n_z}\" nframes=\"{n_t}\" "
            f"pixelwidth=\"{calib['dx_um']}\" pixelheight=\"{calib['dy_um']}\" "
            f"voxeldepth=\"{calib['dz_um']}\" timeinterval=\"{calib['dt_s']}\" /></Settings>\n"
        )
        f.write("</TrackMate>\n")


def run_sparse_lap_tracking(
    labels_path: Path,
    norm_path: Path,
    out_dir: Path,
    movie_id: str,
    cfg: Dict[str, Any],
    calib: Dict[str, float],
) -> Tuple[int, int, int, int, int]:
    """
    Executes two-phase Jaqaman LAP tracking (Frame-to-Frame + Gap Closing + Splitting + Merging).
    Returns (num_tracks, num_spots, num_edges, num_fissions, num_fusions).
    """
    print(f"[{movie_id}] Loading 4D instance labels and normalized intensity stacks...")
    labels_4d, _ = read_tiff(labels_path, axes="TZYX")
    norm_4d, _ = read_tiff(norm_path, axes="TZYX")
    n_t, n_z, n_y, n_x = labels_4d.shape
    dt = calib["dt_s"]

    t_cfg = cfg.get("stage5_tracking_and_mitolandscape", {}).get("sparse_lap_tracker", {})
    max_link_dist = float(t_cfg.get("linking_max_distance_um", 2.5))
    vol_weight = float(t_cfg.get("linking_feature_penalties", {}).get("VOLUME", 1.0))
    allow_gap = bool(t_cfg.get("allow_gap_closing", True))
    max_gap_dist = float(t_cfg.get("gap_closing_max_distance_um", 3.0))
    max_frame_gap = int(t_cfg.get("max_frame_gap", 2))
    allow_split = bool(t_cfg.get("allow_track_splitting", True))
    max_split_dist = float(t_cfg.get("splitting_max_distance_um", 2.0))
    allow_merge = bool(t_cfg.get("allow_track_merging", True))
    max_merge_dist = float(t_cfg.get("merging_max_distance_um", 2.0))
    min_track_spots = int(t_cfg.get("min_track_spots", 4))
    # Max allowed parent/(daughter1 + daughter2) signal ratio (and its inverse) for split/merge events
    max_ratio = max(1.0 + 1e-6, float(t_cfg.get("split_merge_max_intensity_ratio", 1.5)))

    # 1. Fast Bounding-Box Spot Extraction
    print(f"[{movie_id}] Extracting 3D spots via fast bounding-box slicing across {n_t} frames...")
    spots_by_frame, all_spots = extract_spots_fast_bounding_boxes(labels_4d, norm_4d, calib)
    spot_lut = {s["ID"]: s for s in all_spots}

    # 2. Phase I: Frame-to-Frame Linking with Volume Feature Penalty
    print(f"[{movie_id}] Phase I: Solving frame-to-frame LAP with volume penalty...")
    edges: List[Dict[str, Any]] = []
    out_degree: Dict[int, int] = defaultdict(int)
    in_degree: Dict[int, int] = defaultdict(int)

    for t in range(n_t - 1):
        spots_curr = spots_by_frame[t]
        spots_next = spots_by_frame[t + 1]
        if not spots_curr or not spots_next:
            continue

        pts_c = np.array([[s["POSITION_X"], s["POSITION_Y"], s["POSITION_Z"]] for s in spots_curr])
        pts_n = np.array([[s["POSITION_X"], s["POSITION_Y"], s["POSITION_Z"]] for s in spots_next])
        vols_c = np.array([s["VOLUME"] for s in spots_curr], dtype=np.float32)[:, None]
        vols_n = np.array([s["VOLUME"] for s in spots_next], dtype=np.float32)[None, :]

        # Vectorized pairwise distance and TrackMate feature penalty matrix
        dist_mat = cdist(pts_c, pts_n, metric="euclidean")
        vol_diff = np.abs(vols_c - vols_n) / np.maximum(1e-5, np.maximum(vols_c, vols_n))
        cost_mat = np.where(
            dist_mat <= max_link_dist,
            dist_mat * (1.0 + vol_weight * vol_diff),
            1e6,
        )

        row_ind, col_ind = linear_sum_assignment(cost_mat)
        for r, c in zip(row_ind, col_ind):
            if dist_mat[r, c] <= max_link_dist:
                sp_c = spots_curr[r]
                sp_n = spots_next[c]
                d = float(dist_mat[r, c])
                dt_link = max(1e-5, sp_n["POSITION_T"] - sp_c["POSITION_T"])
                edges.append({
                    "SPOT_SOURCE_ID": sp_c["ID"],
                    "SPOT_TARGET_ID": sp_n["ID"],
                    "LINK_COST": round(float(cost_mat[r, c]), 4),
                    "DISPLACEMENT": round(d, 4),
                    "VELOCITY": round(d / dt_link, 4),
                    "DIRECTION_X": round((sp_n["POSITION_X"] - sp_c["POSITION_X"]) / max(1e-5, d), 4),
                    "DIRECTION_Y": round((sp_n["POSITION_Y"] - sp_c["POSITION_Y"]) / max(1e-5, d), 4),
                    "DIRECTION_Z": round((sp_n["POSITION_Z"] - sp_c["POSITION_Z"]) / max(1e-5, d), 4),
                    "EDGE_TIME": round(0.5 * (sp_c["POSITION_T"] + sp_n["POSITION_T"]), 2),
                    "EVENT_TYPE": "FRAME_LINK",
                })
                out_degree[sp_c["ID"]] += 1
                in_degree[sp_n["ID"]] += 1

    # 3. Phase II-A: Multi-Frame Gap Closing (Bridging 1..max_frame_gap missing frames)
    if allow_gap and max_frame_gap >= 1:
        print(f"[{movie_id}] Phase II-A: Solving multi-frame gap closing (max gap = {max_frame_gap} frames)...")
        track_ends = [s for s in all_spots if out_degree[s["ID"]] == 0 and s["FRAME"] < n_t - 1]
        track_starts = [s for s in all_spots if in_degree[s["ID"]] == 0 and s["FRAME"] > 0]

        if track_ends and track_starts:
            gap_cost = np.full((len(track_ends), len(track_starts)), 1e6, dtype=np.float32)
            for i, s_end in enumerate(track_ends):
                p_end = np.array([s_end["POSITION_X"], s_end["POSITION_Y"], s_end["POSITION_Z"]])
                for j, s_start in enumerate(track_starts):
                    f_diff = s_start["FRAME"] - s_end["FRAME"]
                    if 1 < f_diff <= (max_frame_gap + 1):
                        p_start = np.array([s_start["POSITION_X"], s_start["POSITION_Y"], s_start["POSITION_Z"]])
                        d = float(np.linalg.norm(p_start - p_end))
                        if d <= max_gap_dist:
                            gap_cost[i, j] = d

            r_ind, c_ind = linear_sum_assignment(gap_cost)
            for r, c in zip(r_ind, c_ind):
                if gap_cost[r, c] <= max_gap_dist:
                    sp_c = track_ends[r]
                    sp_n = track_starts[c]
                    d = float(gap_cost[r, c])
                    dt_gap = max(1e-5, sp_n["POSITION_T"] - sp_c["POSITION_T"])
                    edges.append({
                        "SPOT_SOURCE_ID": sp_c["ID"],
                        "SPOT_TARGET_ID": sp_n["ID"],
                        "LINK_COST": round(d, 4),
                        "DISPLACEMENT": round(d, 4),
                        "VELOCITY": round(d / dt_gap, 4),
                        "DIRECTION_X": round((sp_n["POSITION_X"] - sp_c["POSITION_X"]) / max(1e-5, d), 4),
                        "DIRECTION_Y": round((sp_n["POSITION_Y"] - sp_c["POSITION_Y"]) / max(1e-5, d), 4),
                        "DIRECTION_Z": round((sp_n["POSITION_Z"] - sp_c["POSITION_Z"]) / max(1e-5, d), 4),
                        "EDGE_TIME": round(0.5 * (sp_c["POSITION_T"] + sp_n["POSITION_T"]), 2),
                        "EVENT_TYPE": "GAP_CLOSE",
                    })
                    out_degree[sp_c["ID"]] += 1
                    in_degree[sp_n["ID"]] += 1

    # 4. Phase II-B: Track Splitting (Fission) & Merging (Fusion) with a conservation gate.
    # A split P(t) -> {C1, C2}(t+1) is only plausible if the parent's mito signal is roughly the
    # sum of its daughters' (Jaqaman et al. 2008 intensity-ratio criterion); likewise for merges.
    # Without this, every unrelated spot appearing/vanishing within the distance cap becomes an
    # "event" (the 2026-11-24 pilot produced 470 spurious fissions in frame 0 this way).
    frame_child: Dict[int, int] = {}
    frame_parent: Dict[int, int] = {}
    for e in edges:
        if e["EVENT_TYPE"] == "FRAME_LINK":
            frame_child[e["SPOT_SOURCE_ID"]] = e["SPOT_TARGET_ID"]
            frame_parent[e["SPOT_TARGET_ID"]] = e["SPOT_SOURCE_ID"]

    def _pos(sp: Dict[str, Any]) -> np.ndarray:
        return np.array([sp["POSITION_X"], sp["POSITION_Y"], sp["POSITION_Z"]], dtype=np.float64)

    def _ratio(whole: Dict[str, Any], part_a: Dict[str, Any], part_b: Dict[str, Any]) -> float:
        trio = (whole, part_a, part_b)
        key = "TOTAL_INTENSITY_CH2" if all(sp["TOTAL_INTENSITY_CH2"] > 0 for sp in trio) else "VOLUME"
        return float(whole[key]) / max(1e-9, float(part_a[key]) + float(part_b[key]))

    def _event_edge(src: Dict[str, Any], tgt: Dict[str, Any], dist: float, cost: float, etype: str) -> None:
        edges.append({
            "SPOT_SOURCE_ID": src["ID"],
            "SPOT_TARGET_ID": tgt["ID"],
            "LINK_COST": round(cost, 4),
            "DISPLACEMENT": round(dist, 4),
            "VELOCITY": round(dist / dt, 4),
            "DIRECTION_X": round((tgt["POSITION_X"] - src["POSITION_X"]) / max(1e-5, dist), 4),
            "DIRECTION_Y": round((tgt["POSITION_Y"] - src["POSITION_Y"]) / max(1e-5, dist), 4),
            "DIRECTION_Z": round((tgt["POSITION_Z"] - src["POSITION_Z"]) / max(1e-5, dist), 4),
            "EDGE_TIME": round(0.5 * (src["POSITION_T"] + tgt["POSITION_T"]), 2),
            "EVENT_TYPE": etype,
        })
        out_degree[src["ID"]] += 1
        in_degree[tgt["ID"]] += 1

    r_lo, r_hi = 1.0 / max_ratio, max_ratio
    for t in range(n_t - 1):
        spots_curr = spots_by_frame[t]
        spots_next = spots_by_frame[t + 1]
        if not spots_curr or not spots_next:
            continue

        if allow_split:
            # parent P (t) already continues to C1 (t+1); unlinked C2 (t+1) may be its second daughter
            candidates = []
            parents = [sp for sp in spots_curr if out_degree[sp["ID"]] == 1 and sp["ID"] in frame_child]
            children = [sp for sp in spots_next if in_degree[sp["ID"]] == 0]
            for par in parents:
                c1 = spot_lut[frame_child[par["ID"]]]
                for c2 in children:
                    d = float(np.linalg.norm(_pos(c2) - _pos(par)))
                    if d > max_split_dist:
                        continue
                    rho = _ratio(par, c1, c2)
                    if r_lo <= rho <= r_hi:
                        candidates.append((d * max(rho, 1.0 / rho), d, par, c2))
            for cost, d, par, c2 in sorted(candidates, key=lambda c: c[0]):
                if in_degree[c2["ID"]] == 0 and out_degree[par["ID"]] < 2:
                    _event_edge(par, c2, d, cost, "SPLIT_FISSION")

        if allow_merge:
            # target M (t+1) already continues from A (t); unlinked end B (t) may be its second source
            candidates = []
            targets = [sp for sp in spots_next if in_degree[sp["ID"]] == 1 and sp["ID"] in frame_parent]
            ends = [sp for sp in spots_curr if out_degree[sp["ID"]] == 0]
            for tgt in targets:
                a = spot_lut[frame_parent[tgt["ID"]]]
                for b in ends:
                    d = float(np.linalg.norm(_pos(tgt) - _pos(b)))
                    if d > max_merge_dist:
                        continue
                    rho = _ratio(tgt, a, b)
                    if r_lo <= rho <= r_hi:
                        candidates.append((d * max(rho, 1.0 / rho), d, b, tgt))
            for cost, d, b, tgt in sorted(candidates, key=lambda c: c[0]):
                if out_degree[b["ID"]] == 0 and in_degree[tgt["ID"]] < 2:
                    _event_edge(b, tgt, d, cost, "MERGE_FUSION")

    # 5. Connected-Component Lineage Graph Construction (Unifying Split & Merged Branches)
    adj: Dict[int, Set[int]] = defaultdict(set)
    for e in edges:
        u, v = e["SPOT_SOURCE_ID"], e["SPOT_TARGET_ID"]
        adj[u].add(v)
        adj[v].add(u)

    visited: Set[int] = set()
    valid_tracks: Dict[int, List[int]] = {}
    spot_to_track: Dict[int, int] = {}
    tid_counter = 0

    for s in all_spots:
        sid = s["ID"]
        if sid not in visited and sid in adj:
            # BFS traversal to collect all spots in the branched lineage tree
            comp: List[int] = []
            queue = [sid]
            visited.add(sid)
            while queue:
                curr = queue.pop(0)
                comp.append(curr)
                for nb in adj[curr]:
                    if nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
            if len(comp) >= min_track_spots:
                comp.sort(key=lambda idx: (spot_lut[idx]["FRAME"], idx))
                valid_tracks[tid_counter] = comp
                for member_id in comp:
                    spot_to_track[member_id] = tid_counter
                tid_counter += 1

    valid_spot_ids = set(spot_to_track.keys())
    filtered_spots = [s for s in all_spots if s["ID"] in valid_spot_ids]
    filtered_edges = [
        e for e in edges if e["SPOT_SOURCE_ID"] in valid_spot_ids and e["SPOT_TARGET_ID"] in valid_spot_ids
    ]

    tracks_edges: Dict[int, List[Dict[str, Any]]] = defaultdict(list)
    for e in filtered_edges:
        tid = spot_to_track[e["SPOT_SOURCE_ID"]]
        e["TRACK_ID"] = tid
        tracks_edges[tid].append(e)

    for s in filtered_spots:
        s["TRACK_ID"] = spot_to_track[s["ID"]]
        s["OUT_DEGREE"] = out_degree[s["ID"]]
        s["IN_DEGREE"] = in_degree[s["ID"]]

    # 6. Export CSVs & Complete TrackMate v7 XML
    spots_csv_path = out_dir / f"{movie_id}_spots.csv"
    edges_csv_path = out_dir / f"{movie_id}_edges.csv"
    tracks_csv_path = out_dir / f"{movie_id}_tracks.csv"
    xml_path = out_dir / f"{movie_id}_tracks.xml"

    print(
        f"[{movie_id}] Exporting spots.csv ({len(filtered_spots)} spots), "
        f"edges.csv ({len(filtered_edges)} links), and tracks.xml..."
    )
    # Always write the CSVs (header-only when empty) so a movie with no kept tracks still yields a
    # complete, schema-valid deliverable set and Stage 5B does not fail on a missing file.
    spot_fields = list(filtered_spots[0].keys()) if filtered_spots else [
        "ID", "NAME", "FRAME", "POSITION_T", "POSITION_X", "POSITION_Y", "POSITION_Z",
        "MEAN_INTENSITY_CH1", "MEAN_INTENSITY_CH2", "MAX_INTENSITY_CH2", "TOTAL_INTENSITY_CH2",
        "RADIUS", "VOLUME", "QUALITY", "label_id", "TRACK_ID", "OUT_DEGREE", "IN_DEGREE",
    ]
    with open(spots_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=spot_fields)
        w.writeheader()
        w.writerows(filtered_spots)

    keys_e = [
        "TRACK_ID", "SPOT_SOURCE_ID", "SPOT_TARGET_ID", "LINK_COST",
        "DISPLACEMENT", "VELOCITY", "DIRECTION_X", "DIRECTION_Y",
        "DIRECTION_Z", "EDGE_TIME", "EVENT_TYPE",
    ]
    with open(edges_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys_e, extrasaction="ignore")
        w.writeheader()
        w.writerows(filtered_edges)

    track_rows = []
    total_fissions = 0
    total_fusions = 0
    for tid, sids in sorted(valid_tracks.items()):
        tr_spots = [spot_lut[sid] for sid in sids]
        tr_edges = tracks_edges[tid]
        dur = round((tr_spots[-1]["FRAME"] - tr_spots[0]["FRAME"]) * dt, 2)
        p_start = np.array([tr_spots[0]["POSITION_X"], tr_spots[0]["POSITION_Y"], tr_spots[0]["POSITION_Z"]])
        p_end = np.array([tr_spots[-1]["POSITION_X"], tr_spots[-1]["POSITION_Y"], tr_spots[-1]["POSITION_Z"]])
        net_disp = round(float(np.linalg.norm(p_end - p_start)), 4)
        total_dist = round(float(sum(e["DISPLACEMENT"] for e in tr_edges)), 4)
        confinement = round(net_disp / max(1e-5, total_dist), 4) if total_dist > 0 else 0.0
        vels = [e["VELOCITY"] for e in tr_edges]

        n_splits = sum(1 for sid in sids if out_degree[sid] >= 2)
        n_merges = sum(1 for sid in sids if in_degree[sid] >= 2)
        total_fissions += n_splits
        total_fusions += n_merges

        track_rows.append({
            "TRACK_ID": tid,
            "NUMBER_SPOTS": len(sids),
            "NUMBER_GAPS": sum(1 for e in tr_edges if e["EVENT_TYPE"] == "GAP_CLOSE"),
            "NUMBER_SPLITS": n_splits,
            "NUMBER_MERGES": n_merges,
            "TRACK_DURATION": dur,
            "TRACK_DISPLACEMENT": net_disp,
            "TOTAL_DISTANCE_TRAVELED": total_dist,
            "CONFINEMENT_RATIO": confinement,
            "TRACK_MEAN_SPEED": round(float(np.mean(vels)) if vels else 0.0, 4),
            "TRACK_MAX_SPEED": round(float(np.max(vels)) if vels else 0.0, 4),
        })

    with open(tracks_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "TRACK_ID", "NUMBER_SPOTS", "NUMBER_GAPS", "NUMBER_SPLITS", "NUMBER_MERGES",
                "TRACK_DURATION", "TRACK_DISPLACEMENT", "TOTAL_DISTANCE_TRAVELED",
                "CONFINEMENT_RATIO", "TRACK_MEAN_SPEED", "TRACK_MAX_SPEED",
            ],
        )
        w.writeheader()
        w.writerows(track_rows)

    write_trackmate_xml(
        xml_path, movie_id, calib, (n_t, n_z, n_y, n_x), spots_by_frame, valid_spot_ids, tracks_edges
    )

    return len(valid_tracks), len(filtered_spots), len(filtered_edges), total_fissions, total_fusions


def main() -> int:
    parser = argparse.ArgumentParser(description="Module 3: Stage 5A TrackMate 7 Sparse LAP Tracking")
    parser.add_argument("--config", type=Path, required=True, help="Path to pipeline_config.yaml")
    parser.add_argument("--movie-id", type=str, required=True, help="Movie identifier stem")
    parser.add_argument("--work-dir", type=Path, required=True, help="Scratch working directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without running tracker")
    args = parser.parse_args()

    cfg = load_config(args.config)
    labels_path = args.work_dir / f"{args.movie_id}_mito_labels.tif"
    norm_path = args.work_dir / f"{args.movie_id}_mito_norm.tif"

    if not labels_path.exists() or not norm_path.exists():
        print(f"[ERROR] Required inputs missing for tracking in {args.work_dir}", file=sys.stderr)
        return 1

    if args.dry_run:
        print(f"[DRY-RUN] Verified {labels_path} and {norm_path} ready for Stage 5A tracking.")
        return 0

    prov_path = args.work_dir / f"{args.movie_id}_provenance.json"
    if prov_path.exists():
        with open(prov_path, "r", encoding="utf-8") as f:
            prov = json.load(f)
            calib = prov.get("calibration", {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7})
    else:
        prov = {"movie_id": args.movie_id}
        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}

    num_tr, num_sp, num_ed, n_fiss, n_fuse = run_sparse_lap_tracking(
        labels_path, norm_path, args.work_dir, args.movie_id, cfg, calib
    )

    prov["stage5_tracking_summary"] = {
        "total_tracks_filtered": num_tr,
        "total_spots_tracked": num_sp,
        "total_edges": num_ed,
        "total_fission_events": n_fiss,
        "total_fusion_events": n_fuse,
        "engine": "Jaqaman_SparseLAP_2Phase",
    }
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2)

    print(
        f"[SUCCESS] Tracking completed for {args.movie_id}: "
        f"{num_tr} tracks, {num_sp} spots, {num_ed} links ({n_fiss} fissions, {n_fuse} fusions)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())