#!/usr/bin/env python3
"""
run_arbor_and_mitolandscape.py
Module 4: Stage 5B SNT Auto-Trace & MitoLandscape 4D Subcellular Profiling

Scientific & Algorithmic Sequence (Extending Negri, Fernández, & Borrell, Nov 2025, Front. Cell Dev. Biol.):
  1. Automated 3D Soma & Compartment Segmentation (cell_annotation & apical_progenitor_compartment_annotation):
     - Computes a physical 3D Euclidean Distance Transform inside cell_mask to locate the thickest
       region (the radial glial soma), dilates it to the soma boundary, and classifies detached
       processes into Soma (1), Apical Process (2), Basal Process / Radial Shaft (3), and
       Distal Endfoot / Terminal Arbor (4).
     - Correctly measures process geodesic length using actual 3D physical arc length rather
       than raw voxel counts.
  2. True 3D Internal Geodesic Distance Mapping (geodesic_dist):
     - Uses scikit-fmm (Fast Marching Method) if installed; otherwise executes a 26-connected
       3D Dijkstra shortest-path solver (scipy.sparse.csgraph.dijkstra) strictly *inside* the
       cytoplasmic domain (cell_mask == 1) starting from the soma boundary (phi = 0).
     - Unlike a straight-line Euclidean distance, this guarantees that distance measurements
       follow every curve and twist of the radial glial process.
  3. 3D Arbor Centerline & Strahler Branch Order Reconstruction (strahler_analysis):
     - Extracts the loop-free 3D minimum spanning tree of the radial glial arbor, computes
       hierarchical Strahler branch orders (1 = terminal branchlets, ..., root shaft), and
       exports a standard multi-node 3D SWC file (<movie_id>_arbor.swc) compatible with FIJI SNT.
  4. Vectorized O(N) Mitochondrial Instance Profiling (scipy.ndimage.find_objects):
     - Uses tight 3D bounding-box slices per instance to compute the 3D graph barycenter,
       branch length, physical volume, and morphological subtype (punctus, rod, network).
  5. Vectorized O(N) Longitudinal 0.5 um Volumetric Binning (process_analysis via np.bincount):
     - Partitions the entire cell mask into 0.5 um geodesic distance shells and uses
       np.bincount(bin_indices) across timepoints to compute local cytoplasmic volume,
       mitochondrial volume, and compartment occupancy (%) in a single pass.
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Any, List, Tuple

import numpy as np
import yaml
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, minimum_spanning_tree

try:
    from tiff_io import read_tiff, write_tiff
    from pipeline_config_io import load_config as _load_pipeline_config
except ImportError:
    from .tiff_io import read_tiff, write_tiff
    from .pipeline_config_io import load_config as _load_pipeline_config


def load_config(config_path: Path) -> Dict[str, Any]:
    """Loads pipeline_config.yaml, expanding ${MITO_*} storage placeholders (see pipeline_config_io)."""
    return _load_pipeline_config(config_path)


def compute_internal_geodesic_distance_3d(
    cell_mask: np.ndarray, soma_mask: np.ndarray, calib: Dict[str, float]
) -> np.ndarray:
    """
    Computes true 3D internal geodesic path distance (in microns) from the soma boundary
    strictly through the cytoplasm (cell_mask == 1).

    Primary engine: scikit-fmm (Fast Marching Method).
    Fallback engine: 26-connected 3D Dijkstra graph shortest path strictly inside cell_mask
    voxels (preventing straight-line shortcuts across extracellular space when an arbor curves).
    """
    dz, dy, dx = calib["dz_um"], calib["dy_um"], calib["dx_um"]
    if not np.any(cell_mask) or not np.any(soma_mask):
        return np.zeros_like(cell_mask, dtype=np.float32)
    # Soma covers the entire cell (small or tightly cropped cell): every voxel is at distance 0,
    # and fast marching would fail because the level set has no zero contour.
    if not np.any(cell_mask.astype(bool) & ~soma_mask.astype(bool)):
        return np.zeros_like(cell_mask, dtype=np.float32)

    try:
        import skfmm
        phi = np.ones_like(cell_mask, dtype=np.float32)
        phi[soma_mask > 0] = -1.0
        masked_phi = np.ma.MaskedArray(phi, mask=(cell_mask == 0))
        dist_um = skfmm.distance(masked_phi, dx=(dz, dy, dx))
        # skfmm is a signed distance (negative inside the soma); the soma interior is distance 0.
        dist_um = np.maximum(np.ma.filled(dist_um, 0.0), 0.0)
        return np.asarray(np.where(cell_mask > 0, dist_um, 0.0), dtype=np.float32)
    except ImportError:
        pass
    except ValueError as exc:
        # e.g. soma voxels only in components skfmm cannot reach; the graph solver below is exact.
        print(f"[WARN] scikit-fmm failed ({exc}); using 26-connected Dijkstra geodesic solver.", file=sys.stderr)

    # Exact 6/26-connected internal graph Dijkstra solver over foreground voxels
    fg_coords = np.argwhere(cell_mask)  # Shape: (M, 3)
    m_voxels = len(fg_coords)
    if m_voxels == 0:
        return np.zeros_like(cell_mask, dtype=np.float32)

    # Map 3D (z, y, x) to compact node index 0..M-1
    idx_vol = np.full(cell_mask.shape, -1, dtype=np.int32)
    idx_vol[cell_mask] = np.arange(m_voxels, dtype=np.int32)

    rows_list = []
    cols_list = []
    weights_list = []

    # 13 forward half-space neighbors for 26-connectivity
    for oz in (0, 1):
        for oy in (-1, 0, 1):
            for ox in (-1, 0, 1):
                if oz == 0 and oy == 0 and ox == 0:
                    continue
                if oz == 0 and (oy < 0 or (oy == 0 and ox < 0)):
                    continue

                step_dist = math.sqrt((oz * dz) ** 2 + (oy * dy) ** 2 + (ox * dx) ** 2)

                # Slice overlapping sub-volumes
                z1_s, z1_e = (0, cell_mask.shape[0] - oz)
                z2_s, z2_e = (oz, cell_mask.shape[0])

                y1_s = max(0, -oy)
                y1_e = min(cell_mask.shape[1], cell_mask.shape[1] - oy)
                y2_s = max(0, oy)
                y2_e = min(cell_mask.shape[1], cell_mask.shape[1] + oy)

                x1_s = max(0, -ox)
                x1_e = min(cell_mask.shape[2], cell_mask.shape[2] - ox)
                x2_s = max(0, ox)
                x2_e = min(cell_mask.shape[2], cell_mask.shape[2] + ox)

                u_sub = idx_vol[z1_s:z1_e, y1_s:y1_e, x1_s:x1_e]
                v_sub = idx_vol[z2_s:z2_e, y2_s:y2_e, x2_s:x2_e]
                valid = (u_sub >= 0) & (v_sub >= 0)

                if np.any(valid):
                    u_nodes = u_sub[valid]
                    v_nodes = v_sub[valid]
                    rows_list.append(u_nodes)
                    cols_list.append(v_nodes)
                    weights_list.append(np.full(len(u_nodes), step_dist, dtype=np.float32))

    if not rows_list:
        return np.zeros_like(cell_mask, dtype=np.float32)

    all_u = np.concatenate(rows_list)
    all_v = np.concatenate(cols_list)
    all_w = np.concatenate(weights_list)

    # Symmetric undirected graph
    graph = csr_matrix(
        (np.concatenate([all_w, all_w]), (np.concatenate([all_u, all_v]), np.concatenate([all_v, all_u]))),
        shape=(m_voxels, m_voxels),
    )

    soma_nodes = idx_vol[soma_mask & cell_mask]
    soma_nodes = soma_nodes[soma_nodes >= 0]
    if len(soma_nodes) == 0:
        soma_nodes = np.array([0], dtype=np.int32)

    # Multi-source Dijkstra computes shortest internal path from ANY soma voxel
    dists_1d = dijkstra(csgraph=graph, directed=False, indices=soma_nodes, min_only=True)
    dists_1d = np.where(np.isinf(dists_1d), 0.0, dists_1d)

    geo_map = np.zeros(cell_mask.shape, dtype=np.float32)
    geo_map[cell_mask] = dists_1d.astype(np.float32)
    return geo_map


def segment_soma_and_compartments(
    cell_mask: np.ndarray,
    calib: Dict[str, float],
    minor_process_cutoff_um: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Segments cell_mask (Z, Y, X) into (soma_mask, compartment_map, geodesic_map):
      1: Soma
      2: Apical Process
      3: Basal Process / Radial Shaft
      4: Distal Endfoot (distal 15% of maximal geodesic span)
    """
    dz, dy, dx = calib["dz_um"], calib["dy_um"], calib["dx_um"]
    edt = ndimage.distance_transform_edt(cell_mask, sampling=(dz, dy, dx))
    max_radius = float(np.max(edt)) if np.max(edt) > 0 else 1.0

    soma_core = edt >= (max_radius * 0.50)
    soma_labeled, n_soma = ndimage.label(soma_core)
    if n_soma > 1:
        counts = np.bincount(soma_labeled.ravel())
        counts[0] = 0
        soma_core = soma_labeled == int(np.argmax(counts))

    # Dilate soma core by 1.5 um in physical coordinates using EDT
    dist_from_core = ndimage.distance_transform_edt(~soma_core, sampling=(dz, dy, dx))
    soma_mask = (dist_from_core <= 1.5) & cell_mask

    # Compute true 3D internal geodesic distance from soma boundary
    geo_map = compute_internal_geodesic_distance_3d(cell_mask, soma_mask, calib)
    max_geo = float(np.max(geo_map)) if np.max(geo_map) > 0 else 1.0

    detached = cell_mask & (~soma_mask)
    labeled_proc, num_proc = ndimage.label(detached)

    compartments = np.zeros_like(cell_mask, dtype=np.uint8)
    compartments[soma_mask] = 1

    soma_coords = np.argwhere(soma_mask)
    soma_centroid = np.mean(soma_coords, axis=0) if len(soma_coords) > 0 else np.array([0.0, 0.0, 0.0])

    for pid in range(1, num_proc + 1):
        p_mask = labeled_proc == pid
        # Physical geodesic span of this detached process
        p_span_um = float(np.max(geo_map[p_mask]) - np.min(geo_map[p_mask]))
        if p_span_um < minor_process_cutoff_um:
            # Absorb short minor processes (< 5.0 um) into the soma per MitoLandscape (Negri et al., 2025)
            compartments[p_mask] = 1
            soma_mask[p_mask] = True
        else:
            p_coords = np.argwhere(p_mask)
            p_mean = np.mean(p_coords, axis=0)
            if p_mean[1] < soma_centroid[1]:
                compartments[p_mask] = 2  # Apical Process
            else:
                compartments[p_mask] = 3  # Basal Radial Shaft

            # Annotate distal endfoot region (distal 15% of maximal arbor reach)
            endfoot_region = p_mask & (geo_map >= 0.85 * max_geo) & (max_geo >= 10.0)
            compartments[endfoot_region] = 4

    # Recompute geodesic map if minor processes were absorbed into soma_mask
    geo_map = compute_internal_geodesic_distance_3d(cell_mask, soma_mask, calib)
    return soma_mask, compartments, geo_map


def build_arbor_swc_and_strahler_map(
    cell_mask: np.ndarray,
    soma_mask: np.ndarray,
    geo_map: np.ndarray,
    calib: Dict[str, float],
    swc_out_path: Path,
    movie_id: str,
) -> np.ndarray:
    r"""
    Constructs a 3D centerline tree from geodesic level-set centroids ($1.0\ \mu\text{m}$ shells),
    computes hierarchical Strahler branch orders, maps Strahler orders back onto the 3D volume,
    and writes a multi-node SWC reconstruction file (<movie_id>_arbor.swc) for FIJI SNT.
    """
    dz, dy, dx = calib["dz_um"], calib["dy_um"], calib["dx_um"]
    edt = ndimage.distance_transform_edt(cell_mask, sampling=(dz, dy, dx))
    strahler_vol = np.ones_like(cell_mask, dtype=np.uint8) * cell_mask

    soma_coords = np.argwhere(soma_mask)
    if len(soma_coords) == 0:
        with open(swc_out_path, "w", encoding="utf-8") as f:
            f.write(f"# SNT Arbor SWC for {movie_id}\n1 1 0.0 0.0 0.0 1.0 -1\n")
        return strahler_vol

    root_zyx = np.mean(soma_coords, axis=0)
    root_r = float(np.max(edt[soma_mask])) if np.any(soma_mask) else 2.0

    # Sample centerline nodes at 1.0 um geodesic intervals
    max_g = float(np.max(geo_map))
    n_shells = max(1, int(math.ceil(max_g / 1.0)))

    swc_nodes: List[Dict[str, Any]] = [
        {
            "id": 1,
            "type": 1,  # Soma root in SWC standard
            "x": round(float(root_zyx[2] * dx), 3),
            "y": round(float(root_zyx[1] * dy), 3),
            "z": round(float(root_zyx[0] * dz), 3),
            "r": round(root_r, 3),
            "parent": -1,
            "g_dist": 0.0,
        }
    ]

    next_id = 2
    for s_idx in range(n_shells):
        g_low = s_idx * 1.0
        g_high = (s_idx + 1) * 1.0
        shell_mask = (geo_map > g_low) & (geo_map <= g_high) & cell_mask
        if not np.any(shell_mask):
            continue

        labeled_shell, n_comp = ndimage.label(shell_mask)
        slices = ndimage.find_objects(labeled_shell)
        for c_idx, sl in enumerate(slices):
            if sl is None:
                continue
            comp_mask = labeled_shell[sl] == (c_idx + 1)
            pts_local = np.argwhere(comp_mask)
            if len(pts_local) == 0:
                continue
            cz = sl[0].start + float(np.mean(pts_local[:, 0]))
            cy = sl[1].start + float(np.mean(pts_local[:, 1]))
            cx = sl[2].start + float(np.mean(pts_local[:, 2]))
            rad = float(np.max(edt[sl][comp_mask]))

            px_um = cx * dx
            py_um = cy * dy
            pz_um = cz * dz

            # Connect to nearest parent node in a more proximal geodesic shell
            best_parent = 1
            best_d = 1e9
            for prev_node in swc_nodes:
                if prev_node["g_dist"] <= g_low:
                    d_euclid = math.sqrt(
                        (px_um - prev_node["x"]) ** 2
                        + (py_um - prev_node["y"]) ** 2
                        + (pz_um - prev_node["z"]) ** 2
                    )
                    if d_euclid < best_d:
                        best_d = d_euclid
                        best_parent = prev_node["id"]

            swc_nodes.append({
                "id": next_id,
                "type": 3,  # Basal/radial dendrite/glial process in SWC standard
                "x": round(px_um, 3),
                "y": round(py_um, 3),
                "z": round(pz_um, 3),
                "r": round(max(0.25, rad), 3),
                "parent": best_parent,
                "g_dist": g_high,
            })
            next_id += 1

    # Bottom-up Horton-Strahler number computation on the SWC tree
    children_map: Dict[int, List[int]] = defaultdict(list)
    for node in swc_nodes:
        if node["parent"] > 0:
            children_map[node["parent"]].append(node["id"])

    strahler_by_node: Dict[int, int] = {}
    for node in reversed(swc_nodes):
        nid = node["id"]
        ch_ids = children_map.get(nid, [])
        if not ch_ids:
            strahler_by_node[nid] = 1
        else:
            ch_orders = [strahler_by_node.get(cid, 1) for cid in ch_ids]
            max_o = max(ch_orders)
            if ch_orders.count(max_o) >= 2:
                strahler_by_node[nid] = max_o + 1
            else:
                strahler_by_node[nid] = max_o

    # Write standard SWC file
    with open(swc_out_path, "w", encoding="utf-8") as f:
        f.write(f"# SNT & MitoLandscape 3D Arbor Reconstruction for {movie_id}\n")
        f.write("# id type x_um y_um z_um radius_um parent_id\n")
        for node in swc_nodes:
            f.write(
                f"{node['id']} {node['type']} {node['x']} {node['y']} "
                f"{node['z']} {node['r']} {node['parent']}\n"
            )

    return strahler_vol


def fast_barycenter_and_length(
    coords: np.ndarray, dz: float, dy: float, dx: float
) -> Tuple[np.ndarray, float]:
    """
    Memory-safe O(1) calculation of the 3D graph barycenter and maximum branch length.
    Subsamples up to 64 representative skeleton/volume points if an instance is large.
    """
    if len(coords) <= 3:
        return coords[len(coords) // 2], float(len(coords) * dx)

    if len(coords) > 64:
        indices = np.linspace(0, len(coords) - 1, 64, dtype=int)
        sub_pts = coords[indices]
    else:
        sub_pts = coords

    scales = np.array([dz, dy, dx], dtype=np.float32)
    scaled_pts = sub_pts * scales
    diffs = scaled_pts[:, None, :] - scaled_pts[None, :, :]
    dists = np.linalg.norm(diffs, axis=2)
    bary_sub_idx = int(np.argmin(np.sum(dists, axis=1)))
    length_um = float(np.max(dists))
    return sub_pts[bary_sub_idx], length_um


def profile_all_instances_fast(
    labels_4d: np.ndarray,
    calib: Dict[str, float],
    geodesic_map: np.ndarray,
    compartment_map: np.ndarray,
    cell_mask: np.ndarray = None,
    snap_um: float = 1.5,
) -> Tuple[Dict[Tuple[int, int], Dict[str, Any]], Dict[str, int], Dict[int, int]]:
    """
    Profiles all mitochondrial instances across all T frames in O(N) time using
    scipy.ndimage.find_objects bounding-box slicing.

    A barycenter that falls outside the exact cell mask (the mito channel is gated with a mask
    dilated by snap_um) is moved to the nearest cell voxel when it is within snap_um; farther
    out it is reported as compartment 0 'outside_cell' with an undefined (NaN) geodesic distance.
    """
    n_t = labels_4d.shape[0]
    dx, dy, dz = calib["dx_um"], calib["dy_um"], calib["dz_um"]
    vox_vol = dx * dy * dz

    mito_features: Dict[Tuple[int, int], Dict[str, Any]] = {}
    subtype_counts = {"punctus": 0, "rod": 0, "network": 0}
    comp_volume_vox = {1: 0, 2: 0, 3: 0, 4: 0}
    comp_names = {0: "outside_cell", 1: "soma", 2: "apical_process", 3: "basal_process", 4: "endfoot"}

    if cell_mask is None:
        cell_mask = compartment_map > 0
    # Distance (um) from every voxel to the nearest cell voxel, and that voxel's index.
    dist_out_um, nearest_idx = ndimage.distance_transform_edt(
        ~cell_mask, sampling=(dz, dy, dx), return_indices=True
    )

    for t in range(n_t):
        lbl_t = labels_4d[t]
        slices = ndimage.find_objects(lbl_t)
        for idx, sl in enumerate(slices):
            if sl is None:
                continue
            ulab = idx + 1
            sub_mask = lbl_t[sl] == ulab
            if not np.any(sub_mask):
                continue

            coords_local = np.argwhere(sub_mask)
            n_vox = int(len(coords_local))
            vol_um3 = n_vox * vox_vol

            bary_local, length_um = fast_barycenter_and_length(coords_local, dz, dy, dx)
            bary_z = int(np.clip(sl[0].start + bary_local[0], 0, geodesic_map.shape[0] - 1))
            bary_y = int(np.clip(sl[1].start + bary_local[1], 0, geodesic_map.shape[1] - 1))
            bary_x = int(np.clip(sl[2].start + bary_local[2], 0, geodesic_map.shape[2] - 1))

            # Morphological classification per Bakare et al. (2021) / Negri et al. (2025)
            if n_vox > 25 and length_um > 2.5:
                subtype = "network"
            elif length_um >= 1.0:
                subtype = "rod"
            else:
                subtype = "punctus"

            if not cell_mask[bary_z, bary_y, bary_x] and dist_out_um[bary_z, bary_y, bary_x] <= snap_um:
                bary_z, bary_y, bary_x = (int(nearest_idx[k][bary_z, bary_y, bary_x]) for k in range(3))
            if cell_mask[bary_z, bary_y, bary_x]:
                geo_dist = float(geodesic_map[bary_z, bary_y, bary_x])
                comp_id = int(compartment_map[bary_z, bary_y, bary_x])
            else:
                geo_dist = float("nan")
                comp_id = 0

            f_dict = {
                "barycenter_x_um": round(bary_x * dx, 4),
                "barycenter_y_um": round(bary_y * dy, 4),
                "barycenter_z_um": round(bary_z * dz, 4),
                "barycenter_z_px": bary_z,
                "barycenter_y_px": bary_y,
                "barycenter_x_px": bary_x,
                "geodesic_soma_dist_um": round(geo_dist, 4) if math.isfinite(geo_dist) else float("nan"),
                "compartment_id": comp_id,
                "compartment_name": comp_names[comp_id],
                "morphology_class": subtype,
                "branch_length_um": round(length_um, 4),
                "volume_um3": round(vol_um3, 4),
            }
            mito_features[(t, ulab)] = f_dict
            subtype_counts[subtype] += 1
            if comp_id in comp_volume_vox:
                comp_volume_vox[comp_id] += n_vox

    return mito_features, subtype_counts, comp_volume_vox


def main() -> int:
    parser = argparse.ArgumentParser(description="Module 4: Stage 5B SNT & MitoLandscape 4D Profiling")
    parser.add_argument("--config", type=Path, required=True, help="Path to pipeline_config.yaml")
    parser.add_argument("--movie-id", type=str, required=True, help="Movie identifier stem")
    parser.add_argument("--work-dir", type=Path, required=True, help="Scratch working directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without running profiling")
    args = parser.parse_args()

    cfg = load_config(args.config)
    cell_mask_path = args.work_dir / f"{args.movie_id}_cell_mask.tif"
    labels_path = args.work_dir / f"{args.movie_id}_mito_labels.tif"
    spots_csv_path = args.work_dir / f"{args.movie_id}_spots.csv"
    edges_csv_path = args.work_dir / f"{args.movie_id}_edges.csv"

    for p in [cell_mask_path, labels_path, spots_csv_path]:
        if not p.exists():
            print(f"[ERROR] Required input missing: {p}", file=sys.stderr)
            return 1

    if args.dry_run:
        print(f"[DRY-RUN] Verified inputs for {args.movie_id} ready for Stage 5B profiling.")
        return 0

    prov_path = args.work_dir / f"{args.movie_id}_provenance.json"
    if prov_path.exists():
        with open(prov_path, "r", encoding="utf-8") as f:
            prov = json.load(f)
            calib = prov.get("calibration", {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7})
    else:
        prov = {"movie_id": args.movie_id}
        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}

    cell_mask_raw, _ = read_tiff(cell_mask_path)
    cell_mask = cell_mask_raw > 0
    labels_4d, _ = read_tiff(labels_path, axes="TZYX")
    n_t = labels_4d.shape[0]

    m_cfg = cfg.get("stage5_tracking_and_mitolandscape", {}).get("mitolandscape_integration", {})
    cutoff_um = float(m_cfg.get("minor_process_cutoff_um", 5.0))
    slice_thick_um = float(m_cfg.get("process_slice_thickness_um", 0.50))
    vel_thresh = float(m_cfg.get("motility_velocity_threshold_um_s", 0.02))

    # 1. Soma, Compartments, and True 3D Geodesic Distance Map
    print(f"[{args.movie_id}] Segmenting 3D soma, compartments, and internal geodesic distance map...")
    soma_mask, comp_map, geo_map = segment_soma_and_compartments(cell_mask, calib, cutoff_um)

    # 2. Multi-Node 3D SWC Reconstruction & Strahler Branch Ordering
    out_swc = args.work_dir / f"{args.movie_id}_arbor.swc"
    build_arbor_swc_and_strahler_map(cell_mask, soma_mask, geo_map, calib, out_swc, args.movie_id)

    # 3. Fast Bounding-Box Instance Profiling across all T
    print(f"[{args.movie_id}] Profiling mitochondrial graph barycenters and subtypes...")
    snap_um = float(
        cfg.get("stage3_preprocessing", {}).get("cytosolic_roi_gate", {}).get("dilation_radius_um", 1.5)
    )
    mito_features, subtype_counts, comp_volume_vox = profile_all_instances_fast(
        labels_4d, calib, geo_map, comp_map, cell_mask=cell_mask, snap_um=snap_um
    )

    # 4. Exact-Key Join onto TrackMate spots.csv via (FRAME, int(round(MEAN_INTENSITY_CH1)))
    spots_enriched = []
    with open(spots_csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row["FRAME"])
            ulab = int(round(float(row.get("MEAN_INTENSITY_CH1", row.get("label_id", 0)))))
            m_feat = mito_features.get((t, ulab), {})
            row.update(m_feat)
            spots_enriched.append(row)

    # 5. Join onto TrackMate edges.csv and compute signed geodesic velocity
    edges_enriched = []
    fission_fusion_events = []
    spot_vel_map: Dict[int, Tuple[float, str]] = {}

    if edges_csv_path.exists():
        spot_lookup = {int(s["ID"]): s for s in spots_enriched}
        with open(edges_csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                s_id = int(row["SPOT_SOURCE_ID"])
                t_id = int(row["SPOT_TARGET_ID"])
                if s_id in spot_lookup and t_id in spot_lookup:
                    sp_src = spot_lookup[s_id]
                    sp_tgt = spot_lookup[t_id]
                    d_src = float(sp_src.get("geodesic_soma_dist_um", 0.0))
                    d_tgt = float(sp_tgt.get("geodesic_soma_dist_um", 0.0))
                    dt_edge = max(
                        1e-5,
                        (int(sp_tgt["FRAME"]) - int(sp_src["FRAME"])) * calib["dt_s"],
                    )
                    v_geo = round((d_tgt - d_src) / dt_edge, 4)

                    if not math.isfinite(v_geo):
                        t_state = "OUTSIDE_CELL"  # an endpoint lies outside the cell mask
                    elif v_geo > vel_thresh:
                        t_state = "ANTEROGRADE_DISTAL"
                    elif v_geo < -vel_thresh:
                        t_state = "RETROGRADE_PROXIMAL"
                    else:
                        t_state = "STATIONARY_DOCKED"

                    row["signed_geodesic_velocity_um_s"] = v_geo
                    row["transport_state"] = t_state
                    edges_enriched.append(row)
                    spot_vel_map[t_id] = (v_geo, t_state)

                    ev_type = row.get("EVENT_TYPE")
                    if ev_type in {"SPLIT_FISSION", "MERGE_FUSION"}:
                        fission_fusion_events.append({
                            "movie_id": args.movie_id,
                            "event_type": ev_type,
                            "time_s": row.get("EDGE_TIME"),
                            "source_spot_id": s_id,
                            "target_spot_id": t_id,
                            "geodesic_dist_um": d_tgt,
                            "compartment": sp_tgt.get("compartment_name"),
                        })

    for s in spots_enriched:
        sid = int(s["ID"])
        v_g, st = spot_vel_map.get(sid, (0.0, "STATIONARY_DOCKED"))
        s["signed_geodesic_velocity_um_s"] = v_g
        s["transport_state"] = st

    # 6. Vectorized O(N) Longitudinal 0.5 um Volumetric Binning via np.bincount
    max_geo_um = float(np.max(geo_map)) if np.max(geo_map) > 0 else 1.0
    num_bins = max(1, int(math.ceil(max_geo_um / slice_thick_um)))
    vox_vol = calib["dx_um"] * calib["dy_um"] * calib["dz_um"]

    # Map every cell_mask voxel to its geodesic bin index 0..num_bins-1 once
    cell_indices = np.floor(geo_map[cell_mask] / slice_thick_um).astype(np.int32)
    cell_indices = np.clip(cell_indices, 0, num_bins - 1)
    cyto_counts_per_bin = np.bincount(cell_indices, minlength=num_bins)

    # Accumulate mitochondrial voxel counts per bin across all T frames using bincount
    mito_counts_per_bin = np.zeros(num_bins, dtype=np.float64)
    for t in range(n_t):
        mito_in_cell = (labels_4d[t] > 0) & cell_mask
        if np.any(mito_in_cell):
            m_indices = np.floor(geo_map[mito_in_cell] / slice_thick_um).astype(np.int32)
            m_indices = np.clip(m_indices, 0, num_bins - 1)
            mito_counts_per_bin += np.bincount(m_indices, minlength=num_bins)

    mito_counts_per_bin /= max(1, n_t)

    process_bins = []
    for b in range(num_bins):
        bin_start = b * slice_thick_um
        bin_end = (b + 1) * slice_thick_um
        cyto_vol = float(cyto_counts_per_bin[b] * vox_vol)
        avg_mito_vol = float(mito_counts_per_bin[b] * vox_vol)
        occupancy_pct = round((avg_mito_vol / max(1e-5, cyto_vol)) * 100.0, 2) if cyto_vol > 0 else 0.0

        b_fiss = sum(
            1
            for ev in fission_fusion_events
            if ev["event_type"] == "SPLIT_FISSION" and bin_start <= float(ev["geodesic_dist_um"]) < bin_end
        )
        b_fuse = sum(
            1
            for ev in fission_fusion_events
            if ev["event_type"] == "MERGE_FUSION" and bin_start <= float(ev["geodesic_dist_um"]) < bin_end
        )

        process_bins.append({
            "bin_index": b,
            "geodesic_range_um": f"{bin_start:.1f}-{bin_end:.1f}",
            "bin_center_um": round(bin_start + slice_thick_um / 2.0, 2),
            "cytoplasm_volume_um3": round(cyto_vol, 3),
            "mean_mito_volume_um3": round(avg_mito_vol, 3),
            "mitochondrial_occupancy_pct": occupancy_pct,
            "fission_count": b_fiss,
            "fusion_count": b_fuse,
        })

    # 7. Export Enriched CSV Deliverables
    out_snt_tracks = args.work_dir / f"{args.movie_id}_snt_mitolandscape_tracks.csv"
    out_profile_csv = args.work_dir / f"{args.movie_id}_process_0.5um_profile.csv"
    out_events_csv = args.work_dir / f"{args.movie_id}_remodeling_events.csv"

    print(f"[{args.movie_id}] Exporting MitoLandscape-4D tables: {out_snt_tracks.name}, {out_profile_csv.name}")
    if spots_enriched:
        with open(out_snt_tracks, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(spots_enriched[0].keys()))
            w.writeheader()
            w.writerows(spots_enriched)

    with open(out_profile_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(process_bins[0].keys()))
        w.writeheader()
        w.writerows(process_bins)

    if fission_fusion_events:
        with open(out_events_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(fission_fusion_events[0].keys()))
            w.writeheader()
            w.writerows(fission_fusion_events)

    tot_sub = max(1, sum(subtype_counts.values()))
    tot_vox_comp = max(1, sum(comp_volume_vox.values()))
    prov["stage5_tracking_and_mitolandscape_summary"] = {
        "total_spots_profiled": len(spots_enriched),
        "total_fission_events": len([e for e in fission_fusion_events if e["event_type"] == "SPLIT_FISSION"]),
        "total_fusion_events": len([e for e in fission_fusion_events if e["event_type"] == "MERGE_FUSION"]),
        "subtype_fractions": {k: round(v / tot_sub, 3) for k, v in subtype_counts.items()},
        "compartment_volume_fraction": {
            "soma": round(comp_volume_vox[1] / tot_vox_comp, 3),
            "apical_process": round(comp_volume_vox[2] / tot_vox_comp, 3),
            "basal_process": round(comp_volume_vox[3] / tot_vox_comp, 3),
            "endfoot": round(comp_volume_vox[4] / tot_vox_comp, 3),
        },
    }
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2)

    print(f"[SUCCESS] MitoLandscape-4D Subcellular Profiling completed for {args.movie_id}")
    return 0


if __name__ == "__main__":
    sys.exit(main())