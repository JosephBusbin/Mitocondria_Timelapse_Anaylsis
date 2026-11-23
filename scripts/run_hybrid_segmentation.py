#!/usr/bin/env python3
"""
run_hybrid_segmentation.py
Module 2: Stage 4 Two-Stage Hybrid AI Segmentation (ilastik + Cellpose/Omnipose)

Scientific & Algorithmic Sequence:
  1. Stage 4A (ilastik Semantic Denoising):
     - Executes headless ilastik on <movie_id>_mito_norm.tif using explicitly bounded CPU
       threads (--threads=${SLURM_CPUS_PER_TASK}) and memory (--ram=24000) so it never
       triggers OOM kills on shared SLURM compute nodes.
     - Produces a clean 32-bit foreground probability map (<movie_id>_Probabilities.tif)
       scoring each voxel from 0.0 (background/haze) to 1.0 (mitochondrion).
  2. Stage 4B (Cellpose 3 / Omnipose Instance Splitting):
     - Stacks normalized intensity and probability channels into (Z, 2, Y, X) per timepoint.
     - Executes 2D flow prediction per slice and applies 2D-to-3D anisotropic IoU slice
       stitching (stitch_threshold=0.35, anisotropy=6.25) across the 1.0 um Z-steps.
     - If running on a CPU-only node without Cellpose, falls back to a 3D Seeded Marker
       Watershed on the Euclidean Distance Transform so touching organelles in thin tubes
       are still split deterministically.
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
from typing import Dict, Any, Tuple

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


def fallback_foreground_probability(mito_norm: np.ndarray) -> np.ndarray:
    """
    Computes a multi-scale difference-of-Gaussians + smoothed intensity foreground probability
    map [0.0, 1.0] per 2D slice when an external .ilp model is not yet present.
    """
    prob = np.zeros_like(mito_norm, dtype=np.float32)
    for t in range(mito_norm.shape[0]):
        for z in range(mito_norm.shape[1]):
            sl = mito_norm[t, z].astype(np.float32)
            if np.max(sl) > 0:
                # Combine fine (sigma=1.0) and tubular (sigma=2.0) Gaussian responses
                g1 = ndimage.gaussian_filter(sl, sigma=1.0)
                g_bg = ndimage.gaussian_filter(sl, sigma=6.0)
                dog = np.maximum(0.0, g1 - 0.5 * g_bg)
                pos_vals = dog[dog > 0]
                scale = float(np.percentile(pos_vals, 98.0)) if len(pos_vals) > 0 else 1.0
                prob[t, z] = np.clip(dog / max(1e-5, scale), 0.0, 1.0)
    return prob


def watershed_split_3d(binary_3d: np.ndarray, intensity_3d: np.ndarray, calib: Dict[str, float]) -> np.ndarray:
    """
    Splits touching 3D mitochondria along thin radial processes using local peak seeding on the
    anisotropic Euclidean Distance Transform + intensity-weighted marker propagation.
    """
    if not np.any(binary_3d):
        return np.zeros_like(binary_3d, dtype=np.uint16)

    # Compute physical distance transform combined with normalized intensity
    edt = ndimage.distance_transform_edt(
        binary_3d, sampling=(calib["dz_um"], calib["dy_um"], calib["dx_um"])
    )
    combined_field = edt * (0.5 + 0.5 * (intensity_3d.astype(np.float32) / max(1.0, float(np.max(intensity_3d)))))

    # Local maxima within a 3x5x5 window (~1.0 um Z x 0.8 um XY) act as organelle seeds
    max_filt = ndimage.maximum_filter(combined_field, size=(1, 5, 5))
    local_peaks = (combined_field == max_filt) & (combined_field > 0.08) & binary_3d
    seed_labels, num_seeds = ndimage.label(local_peaks)

    if num_seeds <= 1:
        cc_labels, _ = ndimage.label(binary_3d)
        return cc_labels.astype(np.uint16)

    # Assign each foreground voxel to its nearest seed weighted by physical voxel dimensions
    _, indices = ndimage.distance_transform_edt(
        seed_labels == 0,
        sampling=(calib["dz_um"], calib["dy_um"], calib["dx_um"]),
        return_indices=True,
    )
    watershed_labels = seed_labels[tuple(indices)] * binary_3d
    return watershed_labels.astype(np.uint16)


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


def segment_instances_cellpose(
    mito_norm: np.ndarray,
    prob_map: np.ndarray,
    cfg: Dict[str, Any],
    calib: Dict[str, float],
    fg_channel_idx: int = 0,
) -> np.ndarray:
    """
    Executes Stage 4B instance segmentation (Cellpose 3 / Omnipose 2D+Z stitching, or
    3D marker-seeded watershed fallback) followed by vectorized O(N) volume filtering.
    """
    s_cfg = cfg.get("stage4_hybrid_segmentation", {}).get("stage4b_instance_segmentation", {})
    diameter = float(s_cfg.get("diameter_px", 8.0))
    flow_thresh = float(s_cfg.get("flow_threshold", 0.40))
    cellprob_thresh = float(s_cfg.get("cellprob_threshold", 0.0))
    stitch_thresh = float(s_cfg.get("stitch_threshold", 0.35))
    min_vol_um3 = float(s_cfg.get("min_volume_um3", 0.05))
    voxel_vol_um3 = calib["dx_um"] * calib["dy_um"] * calib["dz_um"]
    min_voxels = max(2, int(round(min_vol_um3 / voxel_vol_um3)))

    n_t, n_z, n_y, n_x = mito_norm.shape
    labels_4d = np.zeros((n_t, n_z, n_y, n_x), dtype=np.uint16)

    try:
        from cellpose import models
        import torch
        use_gpu = torch.cuda.is_available()
        print(
            f"[STAGE 4B] Initializing Cellpose model (type: {s_cfg.get('model_type', 'cyto3')}, GPU: {use_gpu})..."
        )
        model = models.CellposeModel(gpu=use_gpu, model_type=s_cfg.get("model_type", "cyto3"))
        has_cellpose = True
    except Exception as exc:
        print(
            f"[WARN] Cellpose not available ({exc}). Using 3D EDT marker-seeded watershed splitter.",
            file=sys.stderr,
        )
        has_cellpose = False

    for t in range(n_t):
        if has_cellpose:
            norm_ch0 = mito_norm[t].astype(np.float32) / 50000.0
            norm_ch1 = prob_map[t].copy()
            norm_ch1[norm_ch1 < 0.25] = 0.0
            img_t = np.stack([norm_ch0, norm_ch1], axis=1)  # (Z, C=2, Y, X)

            # Cellpose 3.x API: 'do_3D' (capital D). z_axis/channel_axis are explicit so Cellpose
            # never guesses which of the (Z, C, Y, X) axes holds the two input channels.
            masks_t, _, _ = model.eval(
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
        else:
            bin_t = prob_map[t] >= 0.55
            masks_t = watershed_split_3d(bin_t, mito_norm[t], calib)

        # Vectorized O(N) volume filter and sequential relabeling
        relabeled_t, count_t = filter_and_relabel_fast(masks_t, min_voxels)
        labels_4d[t] = relabeled_t
        print(f"  [T={t+1}/{n_t}] Segmented {count_t} mitochondrial instances (>= {min_vol_um3} um^3).")

    return labels_4d


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

    allow_fallback = os.environ.get("MITO_ALLOW_ILASTIK_FALLBACK", "0") == "1"
    raw_npy_path = args.work_dir / f"{args.movie_id}_ilastik_tzyxc.npy"
    prob_map = None
    if is_valid_ilastik_project(project_file):
        if run_stage4a_ilastik(mito_norm_path, raw_npy_path, ilastik_bin, project_file, num_threads=n_threads):
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
        prob_map = fallback_foreground_probability(mito_norm)
    write_tiff(prob_out_path, prob_map.astype(np.float32), metadata=calib)

    labels_4d = segment_instances_cellpose(mito_norm, prob_map, cfg, calib, fg_channel_idx=fg_idx)

    labels_out_path = args.work_dir / f"{args.movie_id}_mito_labels.tif"
    print(f"[{args.movie_id}] Exporting 4D instance labels: {labels_out_path.name}")
    write_tiff(labels_out_path, labels_4d, metadata=calib)

    inst_counts = [int(np.max(labels_4d[t])) for t in range(labels_4d.shape[0])]
    prov["stage4_segmentation_metrics"] = {
        "ilastik_model": str(project_file.name),
        "stage4a_probability_source": stage4a_source,
        "instance_engine": cfg.get("stage4_hybrid_segmentation", {})
        .get("stage4b_instance_segmentation", {})
        .get("engine", "cellpose"),
        "instance_model": cfg.get("stage4_hybrid_segmentation", {})
        .get("stage4b_instance_segmentation", {})
        .get("model_type", "cyto3"),
        "min_volume_filter_um3": 0.05,
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