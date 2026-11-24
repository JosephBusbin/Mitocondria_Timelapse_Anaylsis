#!/usr/bin/env python3
"""
run_hybrid_segmentation.py
Module 2: Stage 4 Two-Stage Hybrid Segmentation (ilastik + seeded 3D watershed or Cellpose)

Scientific & Algorithmic Sequence:
  1. Stage 4A (ilastik Semantic Denoising):
     - Executes headless ilastik on <movie_id>_mito_norm.tif using explicitly bounded CPU
       threads (--threads=${SLURM_CPUS_PER_TASK}) and memory (--ram=24000) so it never
       triggers OOM kills on shared SLURM compute nodes.
     - Produces a clean 32-bit foreground probability map (<movie_id>_Probabilities.tif)
       scoring each voxel from 0.0 (background/haze) to 1.0 (mitochondrion).
  2. Stage 4B (Instance Splitting: separating TOUCHING mitochondria):
     - ilastik gives a foreground map; mitochondria packed end-to-end in a thin process form one
       connected region, so instances must be split explicitly or counts/volumes/fission-fusion
       calls are wrong.
     - Default engine "watershed": seeded 3D watershed on distance-to-edge x intensity
       (watershed_split_3d), seeds spaced in microns; CPU, seconds per frame.
     - Option engine "cellpose": Cellpose 3 (industry standard; GPU used automatically when
       available), 2D flows per slice + IoU Z-stitching, restricted to the ilastik foreground box.
       Stock 'cyto3' is a whole-cell model and under-segmented mitochondria on practice data
       (see README 2026-11-24); kept for fine-tuned models and the real-data comparison.
     - Either way, instances with mean ilastik P below min_object_mean_probability are dropped.
  3. O(N) Vectorized Physical Volume Filtering & Relabeling:
     - Uses np.bincount and lookup table indexing (lut[masks_t]) to filter out noise
       fragments < 0.05 um^3 (MitoLandscape threshold) and relabel instances 1..N in a
       single O(N) pass (~100x faster than looping over unique IDs).
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

# Bound OpenMP/BLAS threads to SLURM allocation
_slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK", "4")
os.environ.setdefault("OMP_NUM_THREADS", _slurm_cpus)
os.environ.setdefault("OPENBLAS_NUM_THREADS", _slurm_cpus)
os.environ.setdefault("MKL_NUM_THREADS", _slurm_cpus)

import numpy as np
import yaml
from scipy import ndimage

try:
    from tiff_io import read_tiff, write_tiff
    from pipeline_config_io import load_config as _load_pipeline_config
except ImportError:
    from .tiff_io import read_tiff, write_tiff
    from .pipeline_config_io import load_config as _load_pipeline_config


def load_config(config_path: Path) -> Dict[str, Any]:
    """Loads pipeline_config.yaml, expanding ${MITO_*} storage placeholders (see pipeline_config_io)."""
    return _load_pipeline_config(config_path)


HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


def is_valid_ilastik_project(project_file: Path) -> bool:
    """An .ilp is an HDF5 file; a 0-byte or truncated upload must not count as a trained model."""
    try:
        with open(project_file, "rb") as f:
            return f.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE
    except OSError:
        return False


def build_ilastik_command(
    ilastik_bin: str,
    project_file: Path,
    mito_norm_path: Path,
    raw_out_path: Path,
    num_threads: int,
    ram_mb: int,
) -> list:
    """
    Headless ilastik command line. ilastik's 'tif' export is 2D-only and 'multipage tiff' is
    3D-only, so a (T, Z, Y, X, C) probability volume must use an N-D format; 'numpy' (.npy) needs
    no extra Python dependency to read. The axis order is forced so the loader never guesses.
    """
    # --threads and --ram prevent ilastik's LazyFlow backend from spawning 64+ threads on a shared node
    return [
        ilastik_bin,
        "--headless",
        f"--project={project_file}",
        f"--threads={num_threads}",
        f"--ram={ram_mb}",
        "--export_source=Probabilities",
        "--output_format=numpy",
        "--output_axis_order=tzyxc",
        "--export_dtype=float32",
        f"--output_filename_format={raw_out_path}",
        str(mito_norm_path),
    ]


def load_ilastik_probabilities(raw_out_path: Path, fg_idx: int, expected_tzyx: Tuple[int, ...]) -> np.ndarray:
    """Loads ilastik's (T, Z, Y, X, C) .npy export and returns the foreground channel as float32 TZYX."""
    arr = np.load(raw_out_path, mmap_mode="r")
    if arr.ndim == 4:  # single-timepoint input exported as (Z, Y, X, C)
        arr = arr[np.newaxis]
    if arr.ndim != 5:
        raise ValueError(f"ilastik export has {arr.ndim} dims (shape {arr.shape}); expected tzyxc")
    if not 0 <= fg_idx < arr.shape[-1]:
        raise ValueError(
            f"foreground_class_index={fg_idx} but the ilastik project exports {arr.shape[-1]} classes"
        )
    prob = np.asarray(arr[..., fg_idx], dtype=np.float32)
    if prob.shape != tuple(expected_tzyx):
        raise ValueError(f"ilastik probability shape {prob.shape} != mito_norm shape {tuple(expected_tzyx)}")
    return prob


def run_stage4a_ilastik(
    mito_norm_path: Path,
    raw_out_path: Path,
    ilastik_bin: str,
    project_file: Path,
    num_threads: int = 4,
    ram_mb: int = 24000,
) -> bool:
    """
    Executes headless ilastik with explicit thread and RAM limits to respect SLURM cgroups.
    Writes the raw (T, Z, Y, X, C) probabilities to raw_out_path (.npy).
    """
    cmd = build_ilastik_command(ilastik_bin, project_file, mito_norm_path, raw_out_path, num_threads, ram_mb)
    print(f"[STAGE 4A] Running ilastik headless ({num_threads} threads, {ram_mb} MB RAM)...")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if proc.returncode != 0 or not raw_out_path.exists():
        print(f"[WARN] ilastik execution failed (exit {proc.returncode}): {proc.stderr.strip()}", file=sys.stderr)
        return False
    return True


def fallback_foreground_probability(mito_norm: np.ndarray, target_max: float = 50000.0) -> np.ndarray:
    """
    Computes a multi-scale difference-of-Gaussians + smoothed intensity foreground probability
    map [0.0, 1.0] per 2D slice when an external .ilp model is not yet present.

    The DoG response is scaled against a FIXED fraction (25%) of Stage 3's normalized range, not
    each slice's own percentile: a per-slice percentile always maps the brightest pixels of a
    slice to P = 1, so slices containing only residual noise produced confident "mitochondria".
    """
    scale = 0.25 * float(target_max)
    prob = np.zeros_like(mito_norm, dtype=np.float32)
    for t in range(mito_norm.shape[0]):
        for z in range(mito_norm.shape[1]):
            sl = mito_norm[t, z].astype(np.float32)
            if np.max(sl) > 0:
                # Combine fine (sigma=1.0) and tubular (sigma=2.0) Gaussian responses
                g1 = ndimage.gaussian_filter(sl, sigma=1.0)
                g_bg = ndimage.gaussian_filter(sl, sigma=6.0)
                dog = np.maximum(0.0, g1 - 0.5 * g_bg)
                prob[t, z] = np.clip(dog / scale, 0.0, 1.0)
    return prob


def watershed_split_3d(
    binary_3d: np.ndarray,
    intensity_3d: np.ndarray,
    calib: Dict[str, float],
    seed_sep_xy_um: float = 0.4,
    seed_sep_z_um: float = 1.0,
    seed_prominence_um: float = 0.05,
) -> np.ndarray:
    """
    Splits TOUCHING mitochondria inside a foreground mask with a seeded 3D watershed.

    Why this matters: ilastik (Stage 4A) only says which voxels are mitochondria. Two mitochondria
    lying end-to-end or side-by-side in a 1-2 um radial-glia process share one continuous
    foreground region, so plain thresholding/connected components would count them as ONE object
    and corrupt counts, volumes, and fission/fusion calls. The watershed separates them:
      1. elevation = physical distance to the mask edge (um, anisotropic) weighted by smoothed
         intensity -> one "hill" per organelle, highest at its bright core;
      2. seeds = h-maxima: peaks that rise at least seed_prominence_um above the saddle that
         separates them from a neighbour (a flat-topped organelle gives ONE seed, not several);
         seeds closer than seed_sep_xy_um / seed_sep_z_um are merged (both in microns, so they
         mean the same thing at any pixel size);
      3. each region is flooded from its seed until it meets a neighbour; the meeting surface,
         where foreground pinches or dims between organelles, becomes the boundary.
    Any foreground component that received no seed keeps its own label (nothing is dropped).
    """
    from skimage.morphology import h_maxima
    from skimage.segmentation import watershed

    binary_3d = np.asarray(binary_3d, dtype=bool)
    if not np.any(binary_3d):
        return np.zeros(binary_3d.shape, dtype=np.int32)
    dz, dy, dx = calib["dz_um"], calib["dy_um"], calib["dx_um"]

    edt = ndimage.distance_transform_edt(binary_3d, sampling=(dz, dy, dx))
    inten = ndimage.gaussian_filter(intensity_3d.astype(np.float32), sigma=(0.5, 1.0, 1.0))
    inten_norm = inten / max(1e-6, float(inten[binary_3d].max()))
    elevation = (edt * (0.5 + 0.5 * inten_norm)).astype(np.float32)

    peaks = h_maxima(elevation, seed_prominence_um).astype(bool) & binary_3d
    # Merge seeds closer than the minimum separation: dilating each by HALF the separation joins
    # exactly those whose centres are within the full separation of each other.
    rz = max(0, int(round(0.5 * seed_sep_z_um / dz)))
    ry = max(1, int(round(0.5 * seed_sep_xy_um / dy)))
    rx = max(1, int(round(0.5 * seed_sep_xy_um / dx)))
    zz, yy, xx = np.ogrid[-rz : rz + 1, -ry : ry + 1, -rx : rx + 1]
    footprint = (zz / max(rz, 0.5)) ** 2 + (yy / ry) ** 2 + (xx / rx) ** 2 <= 1.0
    markers, _ = ndimage.label(ndimage.binary_dilation(peaks, structure=footprint) & binary_3d)
    markers = markers.astype(np.int32)
    labels = watershed(-elevation, markers, mask=binary_3d).astype(np.int32)

    # Components without a seed (tiny or perfectly flat) keep their own label instead of vanishing
    rest = binary_3d & (labels == 0)
    if rest.any():
        rest_lab, _ = ndimage.label(rest)
        labels[rest] = rest_lab[rest] + int(labels.max())
    return labels


def watershed_split_3d_temporal(
    binary_3d: np.ndarray,
    intensity_3d: np.ndarray,
    calib: Dict[str, float],
    prev_labels: Optional[np.ndarray],
    seed_sep_xy_um: float = 0.4,
    seed_sep_z_um: float = 1.0,
    seed_prominence_um: float = 0.05,
    split_hysteresis: float = 2.0,
    merge_boundary_ratio: float = 0.95,
    min_overlap_voxels: int = 3,
    merge_boundary_percentile: float = 90.0,
) -> np.ndarray:
    """
    Seeded 3D watershed that uses the PREVIOUS frame's instances (already registered to the same
    coordinates by Stage 3) to suppress frame-to-frame segmentation flicker between touching
    mitochondria, which otherwise shows up downstream as fake fission/fusion pairs.

    Per foreground component C of this frame:
      - C overlaps >= 2 previous objects but this frame's own seeds would give fewer: the previous
        objects become the markers (keep them apart). Afterwards each adjacent pair is merged back
        only if there is no intensity dip between them (90th percentile of the boundary intensity >=
        merge_boundary_ratio x the dimmer of the two peaks) -> real fusions are still allowed.
      - C overlaps exactly 1 previous object but this frame's seeds would split it: re-seed C with
        split_hysteresis x seed_prominence_um (a split must be clearly stronger than the noise of
        the previous frame's single object) -> real fissions with a clear gap are still allowed.
      - otherwise (no overlap, e.g. fast movers or new objects): independent seeds, as in
        watershed_split_3d.
    With prev_labels None this is exactly watershed_split_3d.
    """
    from skimage.morphology import h_maxima
    from skimage.segmentation import watershed

    binary_3d = np.asarray(binary_3d, dtype=bool)
    if prev_labels is None or not np.any(prev_labels) or not np.any(binary_3d):
        return watershed_split_3d(binary_3d, intensity_3d, calib, seed_sep_xy_um, seed_sep_z_um, seed_prominence_um)
    dz, dy, dx = calib["dz_um"], calib["dy_um"], calib["dx_um"]
    edt = ndimage.distance_transform_edt(binary_3d, sampling=(dz, dy, dx))
    inten = ndimage.gaussian_filter(intensity_3d.astype(np.float32), sigma=(0.5, 1.0, 1.0))
    inten_norm = inten / max(1e-6, float(inten[binary_3d].max()))
    elevation = (edt * (0.5 + 0.5 * inten_norm)).astype(np.float32)

    rz = max(0, int(round(0.5 * seed_sep_z_um / dz)))
    ry = max(1, int(round(0.5 * seed_sep_xy_um / dy)))
    rx = max(1, int(round(0.5 * seed_sep_xy_um / dx)))
    zz, yy, xx = np.ogrid[-rz : rz + 1, -ry : ry + 1, -rx : rx + 1]
    footprint = (zz / max(rz, 0.5)) ** 2 + (yy / ry) ** 2 + (xx / rx) ** 2 <= 1.0

    def seeds(h: float, region: np.ndarray) -> np.ndarray:
        pk = h_maxima(np.where(region, elevation, 0.0).astype(np.float32), h).astype(bool) & region
        lab, _ = ndimage.label(ndimage.binary_dilation(pk, structure=footprint) & region)
        return lab.astype(np.int32)

    base = seeds(seed_prominence_um, binary_3d)
    comps, n_comp = ndimage.label(binary_3d)
    markers = np.zeros_like(base)
    next_id = 1
    keep_apart = []  # (component id, list of marker ids seeded from the previous frame)
    objs = ndimage.find_objects(comps)
    for ci, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        cmask = comps[sl] == ci
        prev_sub = prev_labels[sl][cmask]
        ids, cnt = np.unique(prev_sub[prev_sub > 0], return_counts=True)
        prev_ids = ids[cnt >= min_overlap_voxels]
        own = np.unique(base[sl][cmask])
        own = own[own > 0]
        msub = markers[sl]
        if len(prev_ids) >= 2 and len(own) < len(prev_ids):
            ids_here = []
            for pid in prev_ids:
                m = cmask & (prev_labels[sl] == pid)
                if m.any():
                    msub[m] = next_id
                    ids_here.append(next_id)
                    next_id += 1
            keep_apart.append((ci, ids_here))
            continue
        if len(prev_ids) == 1 and len(own) > 1:
            region = np.zeros_like(binary_3d)
            region[sl] = cmask
            strict = seeds(split_hysteresis * seed_prominence_um, region)[sl]
            src = strict
        else:
            src = base[sl]
        for oid in np.unique(src[cmask]):
            if oid == 0:
                continue
            msub[cmask & (src == oid)] = next_id
            next_id += 1
    labels = watershed(-elevation, markers, mask=binary_3d).astype(np.int32)

    # Allow real fusions: merge previous-frame-seeded neighbours with no intensity dip between them.
    struct = ndimage.generate_binary_structure(3, 1)
    for ci, ids_here in keep_apart:
        sl = objs[ci - 1]
        sub = labels[sl]
        isub = inten[sl]
        changed = True
        while changed and len(ids_here) > 1:
            changed = False
            for i in range(len(ids_here)):
                for j in range(i + 1, len(ids_here)):
                    a, b = ids_here[i], ids_here[j]
                    ma, mb = sub == a, sub == b
                    if not ma.any() or not mb.any():
                        continue
                    border = (ndimage.binary_dilation(ma, structure=struct) & mb) | (ndimage.binary_dilation(mb, structure=struct) & ma)
                    if not border.any():
                        continue
                    # A real dip means even the brightest part of the shared boundary is dimmer than both
                    # peaks; the boundary's mean would be dragged down by the objects' dim outer edges.
                    if float(np.percentile(isub[border], merge_boundary_percentile)) >= merge_boundary_ratio * min(float(isub[ma].max()), float(isub[mb].max())):
                        sub[mb] = a
                        ids_here.pop(j)
                        changed = True
                        break
                if changed:
                    break

    rest = binary_3d & (labels == 0)
    if rest.any():
        rest_lab, _ = ndimage.label(rest)
        labels[rest] = rest_lab[rest] + int(labels.max())
    return labels


def filter_and_relabel_fast(masks_t: np.ndarray, min_voxels: int) -> Tuple[np.ndarray, int]:
    """
    Vectorized O(N) volume filtering and sequential relabeling (1..K) using np.bincount
    and lookup table indexing, avoiding slow O(K * N) per-object boolean mask scans.
    """
    flat = masks_t.ravel()
    if flat.max() == 0:
        return np.zeros_like(masks_t, dtype=np.uint16), 0

    counts = np.bincount(flat)
    counts[0] = 0  # Background is always 0

    # Keep only labels meeting the physical volume floor (>= min_voxels)
    valid_mask = counts >= min_voxels
    num_valid = int(np.sum(valid_mask))
    if num_valid == 0:
        return np.zeros_like(masks_t, dtype=np.uint16), 0

    # Build lookup table mapping old_label -> new_sequential_label (1..num_valid)
    lut = np.zeros(len(counts), dtype=np.uint16)
    lut[valid_mask] = np.arange(1, num_valid + 1, dtype=np.uint16)
    return lut[masks_t], num_valid


EXIT_IMPLAUSIBLE_FOREGROUND = 5  # run_pipeline_slurm/plan_full_run: not retried (retrain ilastik instead)


def foreground_fraction_of_cell(prob_map: np.ndarray, cell_mask: np.ndarray, threshold: float) -> float:
    """Median over frames of the fraction of cell-mask voxels with P > threshold.

    Mitochondria fill a minority of the cytoplasm. A classifier that calls most of the cell
    "mitochondria" (e.g. trained with Background strokes only outside the cell) would make Stage 4B
    tile the whole cytoplasm into fake instances (seen with bestman_mito_2d_v2.ilp: 99% of the cell).
    """
    mask = np.asarray(cell_mask, dtype=bool)
    if not mask.any():
        return 0.0
    return float(np.median([(prob_map[t] > threshold)[mask].mean() for t in range(prob_map.shape[0])]))


def gate_instances_by_probability(masks: np.ndarray, prob: np.ndarray, min_mean_p: float) -> np.ndarray:
    """Zeroes every instance whose mean ilastik foreground probability is below min_mean_p."""
    n = int(masks.max())
    if n == 0:
        return masks
    mean_p = ndimage.mean(prob, labels=masks, index=np.arange(1, n + 1))
    keep = np.concatenate([[False], np.nan_to_num(np.asarray(mean_p), nan=0.0) >= min_mean_p])
    return np.where(keep[masks], masks, 0).astype(masks.dtype)


def _watershed_frame_job(job: Tuple[Any, ...]) -> Tuple[np.ndarray, int]:
    """
    One frame of the watershed engine: split -> probability gate -> volume filter/relabel.
    Runs on the frame's foreground bounding box only (all instances lie inside it, so gating and
    relabelling the box give the same labels as doing it on the full frame). Top-level function so
    a ProcessPoolExecutor can pickle it. Returns (relabelled uint16 box, instance count).
    """
    fg_box, inten_box, prob_box, calib, seed_xy, seed_z, seed_h, min_mean_p, min_voxels = job
    masks_box = watershed_split_3d(
        fg_box, inten_box, calib,
        seed_sep_xy_um=seed_xy, seed_sep_z_um=seed_z, seed_prominence_um=seed_h,
    )
    masks_box = gate_instances_by_probability(masks_box, prob_box, min_mean_p)
    return filter_and_relabel_fast(masks_box, min_voxels)


def watershed_worker_count(n_jobs: int, n_workers: Optional[int] = None) -> int:
    """
    Worker processes for the per-frame watershed: n_workers if given, else $SLURM_CPUS_PER_TASK
    (1 outside Slurm), never more than the number of frames with foreground.
    """
    if n_workers is None:
        try:
            n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
        except ValueError:
            n_workers = 1
    return max(1, min(int(n_workers), n_jobs))


def segment_instances(
    mito_norm: np.ndarray,
    prob_map: np.ndarray,
    cfg: Dict[str, Any],
    calib: Dict[str, float],
    fg_channel_idx: int = 0,
    n_workers: Optional[int] = None,
) -> Tuple[np.ndarray, str]:
    """
    Stage 4B instance segmentation: turns the ilastik foreground into individual mitochondria,
    splitting TOUCHING organelles, then gates by probability and filters by physical volume.

    engine (pipeline_config.yaml stage4b_instance_segmentation.engine):
      "watershed" (default) - seeded 3D watershed on the probability foreground (watershed_split_3d);
                              CPU, seconds per frame.
      "cellpose"            - Cellpose 3 (industry-standard deep-learning segmenter), 2D per slice
                              + Z-stitching; uses a GPU automatically when torch sees one. The
                              stock 'cyto3' model is trained on whole cells, and on practice
                              Series #39 it missed most mitochondria (18% foreground coverage vs
                              100% for the watershed, 433 s vs 5 s per frame on CPU); a fine-tuned
                              mitochondria model may change that (see manuals and README).
    If Cellpose is requested but cannot be imported, the watershed is used and a warning printed.
    Watershed frames are independent, so they run in parallel worker processes (n_workers, default
    $SLURM_CPUS_PER_TASK; 1 = serial). Each frame is relabelled 1..N on its own, so the result is
    identical to the serial loop. The Cellpose path is always serial (one model, maybe one GPU).
    Returns (labels_4d, engine_used).
    """
    s_cfg = cfg.get("stage4_hybrid_segmentation", {}).get("stage4b_instance_segmentation", {})
    engine = str(s_cfg.get("engine", "watershed")).strip().lower()
    ws_cfg = s_cfg.get("watershed", {})
    seed_xy = float(ws_cfg.get("seed_min_separation_xy_um", 0.4))
    seed_z = float(ws_cfg.get("seed_min_separation_z_um", 1.0))
    seed_h = float(ws_cfg.get("seed_min_prominence_um", 0.05))
    diameter = float(s_cfg.get("diameter_px", 8.0))
    flow_thresh = float(s_cfg.get("flow_threshold", 0.40))
    cellprob_thresh = float(s_cfg.get("cellprob_threshold", 0.0))
    stitch_thresh = float(s_cfg.get("stitch_threshold", 0.35))
    min_vol_um3 = float(s_cfg.get("min_volume_um3", 0.05))
    voxel_vol_um3 = calib["dx_um"] * calib["dy_um"] * calib["dz_um"]
    min_voxels = max(2, int(round(min_vol_um3 / voxel_vol_um3)))

    n_t, n_z, n_y, n_x = mito_norm.shape
    labels_4d = np.zeros((n_t, n_z, n_y, n_x), dtype=np.uint16)

    has_cellpose = False
    if engine in ("cellpose", "omnipose"):
        try:
            from cellpose import models
            import torch
            use_gpu = torch.cuda.is_available()
            model_type = s_cfg.get("model_type", "cyto3")
            print(f"[STAGE 4B] Engine: Cellpose (model: {model_type}, GPU: {use_gpu})...")
            model = models.CellposeModel(gpu=use_gpu, model_type=model_type)
            has_cellpose = True
        except Exception as exc:
            print(f"[WARN] Cellpose requested but unavailable ({exc}); using the 3D watershed.", file=sys.stderr)
    elif engine != "watershed":
        raise ValueError(f"Unknown stage4b engine '{engine}' (expected 'watershed' or 'cellpose')")
    if not has_cellpose:
        print(f"[STAGE 4B] Engine: seeded 3D watershed (seed spacing {seed_xy} um XY, {seed_z} um Z)...")
    engine_used = "cellpose" if has_cellpose else "watershed"

    fg_thresh = float(s_cfg.get("foreground_probability_threshold", 0.5))
    min_mean_p = float(s_cfg.get("min_object_mean_probability", 0.5))
    margin_px = int(round(2 * diameter))

    if not has_cellpose:
        # Watershed engine: build one job per frame with foreground (bounding box + margin, as below),
        # then run them serially or in a process pool; results are written back in frame order.
        jobs, boxes, frames = [], [], []
        for t in range(n_t):
            fg_t = prob_map[t] > fg_thresh
            if not fg_t.any():
                print(f"  [T={t+1}/{n_t}] No ilastik foreground (P > {fg_thresh}); 0 instances.")
                continue
            zs = np.flatnonzero(fg_t.any(axis=(1, 2)))
            ys = np.flatnonzero(fg_t.any(axis=(0, 2)))
            xs = np.flatnonzero(fg_t.any(axis=(0, 1)))
            box = (
                slice(max(0, int(zs[0]) - 1), min(n_z, int(zs[-1]) + 2)),
                slice(max(0, int(ys[0]) - margin_px), min(n_y, int(ys[-1]) + 1 + margin_px)),
                slice(max(0, int(xs[0]) - margin_px), min(n_x, int(xs[-1]) + 1 + margin_px)),
            )
            jobs.append((fg_t[box], mito_norm[t][box], prob_map[t][box], calib,
                         seed_xy, seed_z, seed_h, min_mean_p, min_voxels))
            boxes.append(box)
            frames.append(t)
        temporal = bool(ws_cfg.get("temporal_seeding", True))
        if temporal:
            # Frames depend on the previous frame's instances -> sequential (still seconds per frame).
            print("[STAGE 4B] Temporal seeding ON: previous-frame instances stabilise splits/merges.")
            prev_full = None
            done = {}
            for (t, box, job) in zip(frames, boxes, jobs):
                fg_box, inten_box, prob_box = job[0], job[1], job[2]
                prev_box = prev_full[box] if prev_full is not None else None
                masks_box = watershed_split_3d_temporal(
                    fg_box, inten_box, calib, prev_box,
                    seed_sep_xy_um=seed_xy, seed_sep_z_um=seed_z, seed_prominence_um=seed_h,
                    split_hysteresis=float(ws_cfg.get("temporal_split_hysteresis", 2.0)),
                    merge_boundary_ratio=float(ws_cfg.get("temporal_merge_boundary_ratio", 0.95)),
                    merge_boundary_percentile=float(ws_cfg.get("temporal_merge_boundary_percentile", 90.0)),
                )
                masks_box = gate_instances_by_probability(masks_box, prob_box, min_mean_p)
                relabeled_box, count_t = filter_and_relabel_fast(masks_box, min_voxels)
                labels_4d[t][box] = relabeled_box
                prev_full = labels_4d[t]
                done[t] = count_t
            for t in frames:
                print(f"  [T={t+1}/{n_t}] Segmented {done[t]} mitochondrial instances (>= {min_vol_um3} um^3).")
            return labels_4d, engine_used
        workers = watershed_worker_count(len(jobs), n_workers)
        if workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            print(f"[STAGE 4B] Watershed on {len(jobs)} frame(s) with {workers} parallel worker processes...")
            with ProcessPoolExecutor(max_workers=workers) as pool:
                results = list(pool.map(_watershed_frame_job, jobs))  # map() keeps frame order
        else:
            results = [_watershed_frame_job(job) for job in jobs]
        for t, box, (relabeled_box, count_t) in zip(frames, boxes, results):
            labels_4d[t][box] = relabeled_box
            print(f"  [T={t+1}/{n_t}] Segmented {count_t} mitochondrial instances (>= {min_vol_um3} um^3).")
        return labels_4d, engine_used

    for t in range(n_t):
        fg_t = prob_map[t] > fg_thresh
        if not fg_t.any():
            print(f"  [T={t+1}/{n_t}] No ilastik foreground (P > {fg_thresh}); 0 instances.")
            continue
        # Run Cellpose only on the Z-range / XY box that contains ilastik foreground (plus a margin):
        # empty slices cost the same CPU time as full ones and only contribute noise objects.
        zs = np.flatnonzero(fg_t.any(axis=(1, 2)))
        ys = np.flatnonzero(fg_t.any(axis=(0, 2)))
        xs = np.flatnonzero(fg_t.any(axis=(0, 1)))
        box = (
            # +-1 slice of Z context: helps Z-stitching and avoids single-slice inputs
            slice(max(0, int(zs[0]) - 1), min(n_z, int(zs[-1]) + 2)),
            slice(max(0, int(ys[0]) - margin_px), min(n_y, int(ys[-1]) + 1 + margin_px)),
            slice(max(0, int(xs[0]) - margin_px), min(n_x, int(xs[-1]) + 1 + margin_px)),
        )
        masks_t = np.zeros((n_z, n_y, n_x), dtype=np.int32)

        if has_cellpose:
            norm_ch0 = mito_norm[t][box].astype(np.float32) / 50000.0
            norm_ch1 = prob_map[t][box].copy()
            norm_ch1[norm_ch1 < 0.25] = 0.0
            img_t = np.stack([norm_ch0, norm_ch1], axis=1)  # (Z, C=2, Y, X)

            # Cellpose 3.x API: 'do_3D' (capital D). z_axis/channel_axis are explicit so Cellpose
            # never guesses which of the (Z, C, Y, X) axes holds the two input channels.
            masks_box, _, _ = model.eval(
                img_t,
                diameter=diameter,
                channels=[1, 2],
                channel_axis=1,
                z_axis=0,
                flow_threshold=flow_thresh,
                cellprob_threshold=cellprob_thresh,
                do_3D=False,
                stitch_threshold=stitch_thresh,
                min_size=min_voxels,
            )
            masks_t[box] = np.asarray(masks_box).reshape(masks_t[box].shape)
        else:
            masks_t[box] = watershed_split_3d(
                fg_t[box], mito_norm[t][box], calib,
                seed_sep_xy_um=seed_xy, seed_sep_z_um=seed_z, seed_prominence_um=seed_h,
            )

        # Keep only instances that ilastik agrees are mitochondria (mean P >= min_mean_p).
        masks_t = gate_instances_by_probability(masks_t, prob_map[t], min_mean_p)

        # Vectorized O(N) volume filter and sequential relabeling
        relabeled_t, count_t = filter_and_relabel_fast(masks_t, min_voxels)
        labels_4d[t] = relabeled_t
        print(f"  [T={t+1}/{n_t}] Segmented {count_t} mitochondrial instances (>= {min_vol_um3} um^3).")

    return labels_4d, engine_used


def segment_instances_cellpose(
    mito_norm: np.ndarray,
    prob_map: np.ndarray,
    cfg: Dict[str, Any],
    calib: Dict[str, float],
    fg_channel_idx: int = 0,
) -> np.ndarray:
    """Backward-compatible wrapper: runs segment_instances with the configured engine, returns labels only."""
    return segment_instances(mito_norm, prob_map, cfg, calib, fg_channel_idx)[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Module 2: Stage 4 Two-Stage Hybrid Segmentation")
    parser.add_argument("--config", type=Path, required=True, help="Path to pipeline_config.yaml")
    parser.add_argument("--movie-id", type=str, required=True, help="Movie identifier stem")
    parser.add_argument("--work-dir", type=Path, required=True, help="Scratch working directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without running inference")
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths_cfg = cfg.get("paths", {})
    models_dir = Path(paths_cfg.get("models_dir", args.work_dir))
    ilastik_bin = paths_cfg.get("ilastik_bin", "/opt/ilastik/run_ilastik.sh")

    mito_norm_path = args.work_dir / f"{args.movie_id}_mito_norm.tif"
    if not mito_norm_path.exists():
        print(f"[ERROR] Normalized stack not found: {mito_norm_path}", file=sys.stderr)
        return 1

    if args.dry_run:
        print(f"[DRY-RUN] Verified {mito_norm_path} ready for Stage 4 segmentation.")
        return 0

    prov_path = args.work_dir / f"{args.movie_id}_provenance.json"
    if prov_path.exists():
        with open(prov_path, "r", encoding="utf-8") as f:
            prov = json.load(f)
            calib = prov.get("calibration", {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7})
    else:
        prov = {"movie_id": args.movie_id}
        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}

    mito_norm, _ = read_tiff(mito_norm_path, axes="TZYX")
    prob_out_path = args.work_dir / f"{args.movie_id}_Probabilities.tif"

    ilastik_cfg = cfg.get("stage4_hybrid_segmentation", {}).get("stage4a_ilastik", {})
    project_file = models_dir / ilastik_cfg.get("project_file", "bestman_mito_2d_v1.ilp")
    fg_idx = int(ilastik_cfg.get("foreground_class_index", 0))
    n_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    # ilastik's cache must stay inside the job's cgroup: default to 60% of --mem (SLURM_MEM_PER_NODE
    # is in MB), leaving the rest for this process's copies of mito_norm and the probability map.
    ram_mb = int(
        os.environ.get("MITO_ILASTIK_RAM_MB")
        or max(2000, int(int(os.environ.get("SLURM_MEM_PER_NODE", "8192")) * 0.6))
    )

    allow_fallback = os.environ.get("MITO_ALLOW_ILASTIK_FALLBACK", "0") == "1"
    raw_npy_path = args.work_dir / f"{args.movie_id}_ilastik_tzyxc.npy"
    prob_map = None
    if is_valid_ilastik_project(project_file):
        if run_stage4a_ilastik(
            mito_norm_path, raw_npy_path, ilastik_bin, project_file, num_threads=n_threads, ram_mb=ram_mb
        ):
            try:
                prob_map = load_ilastik_probabilities(raw_npy_path, fg_idx, mito_norm.shape)
            except ValueError as exc:
                print(f"[WARN] {exc}", file=sys.stderr)
            finally:
                raw_npy_path.unlink(missing_ok=True)  # ~2x the size of mito_norm; TIFF copy kept below
        if prob_map is None and not allow_fallback:
            print(
                "[ERROR] A trained ilastik project is present but Stage 4A failed (see warning above). "
                "Set MITO_ALLOW_ILASTIK_FALLBACK=1 to continue with the built-in estimator.",
                file=sys.stderr,
            )
            return 1
    else:
        print(
            f"[WARN] No valid ilastik project at {project_file} (missing, empty, or not HDF5). "
            "Using adaptive multi-scale Hessian/Gaussian foreground estimator.",
            file=sys.stderr,
        )

    stage4a_source = "ilastik" if prob_map is not None else "fallback_estimator"
    if prob_map is None:
        target_max = float(
            cfg.get("stage3_preprocessing", {}).get("intensity_normalization", {}).get("target_max_uint16", 50000)
        )
        prob_map = fallback_foreground_probability(mito_norm, target_max)
    write_tiff(prob_out_path, prob_map.astype(np.float32), metadata=calib)

    # Plausibility guard: stop before segmenting if the foreground covers most of the cell.
    s4_cfg = cfg.get("stage4_hybrid_segmentation", {})
    fg_thr = float(s4_cfg.get("stage4b_instance_segmentation", {}).get("foreground_probability_threshold", 0.5))
    max_frac = float(s4_cfg.get("stage4a_ilastik", {}).get("max_foreground_fraction_of_cell", 0.5))
    mask_path = args.work_dir / f"{args.movie_id}_cell_mask.tif"
    if mask_path.exists():
        cell_mask_raw, _ = read_tiff(mask_path)
        cell_mask = np.asarray(cell_mask_raw).reshape(mito_norm.shape[1:]) > 0
        frac = foreground_fraction_of_cell(prob_map, cell_mask, fg_thr)
        prov["stage4a_foreground_fraction_of_cell"] = round(frac, 4)
        print(f"[STAGE 4A] Foreground (P > {fg_thr}) covers {100 * frac:.1f}% of the cell mask "
              f"(source: {stage4a_source}; limit {100 * max_frac:.0f}%).")
        if frac > max_frac and os.environ.get("MITO_ALLOW_IMPLAUSIBLE_FOREGROUND", "0") != "1":
            prov["stage4_status"] = "FAIL_IMPLAUSIBLE_FOREGROUND"
            with open(prov_path, "w", encoding="utf-8") as f:
                json.dump(prov, f, indent=2)
            print(
                f"[ERROR] {args.movie_id}: {100 * frac:.1f}% of the cell is called mitochondria (> {100 * max_frac:.0f}%). "
                "The foreground model is implausible (for ilastik: retrain with Background strokes on in-cell "
                "texture and between mitochondria). Probabilities kept for inspection: "
                f"{prob_out_path.name}. Override: MITO_ALLOW_IMPLAUSIBLE_FOREGROUND=1.",
                file=sys.stderr,
            )
            return EXIT_IMPLAUSIBLE_FOREGROUND

    labels_4d, engine_used = segment_instances(mito_norm, prob_map, cfg, calib, fg_channel_idx=fg_idx)

    labels_out_path = args.work_dir / f"{args.movie_id}_mito_labels.tif"
    print(f"[{args.movie_id}] Exporting 4D instance labels: {labels_out_path.name}")
    write_tiff(labels_out_path, labels_4d, metadata=calib)

    inst_counts = [int(np.max(labels_4d[t])) for t in range(labels_4d.shape[0])]
    prov["stage4_segmentation_metrics"] = {
        "ilastik_model": str(project_file.name),
        "stage4a_probability_source": stage4a_source,
        "instance_engine": engine_used,
        "instance_model": (
            cfg.get("stage4_hybrid_segmentation", {}).get("stage4b_instance_segmentation", {}).get("model_type", "cyto3")
            if engine_used == "cellpose" else "seeded_3d_watershed"
        ),
        "min_volume_filter_um3": float(
            cfg.get("stage4_hybrid_segmentation", {}).get("stage4b_instance_segmentation", {}).get("min_volume_um3", 0.05)
        ),
        "mean_instances_per_frame": round(float(np.mean(inst_counts)), 1),
        "total_instances_across_time": int(np.sum(inst_counts)),
    }
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2)

    print(
        f"[SUCCESS] Stage 4 Hybrid Segmentation completed for {args.movie_id} "
        f"(Mean instances/frame: {prov['stage4_segmentation_metrics']['mean_instances_per_frame']})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())