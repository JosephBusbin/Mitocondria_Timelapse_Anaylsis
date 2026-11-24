#!/usr/bin/env python3
"""
preprocess_4d_stack.py
Module 1: Stage 3 Adaptive 4D Image Preprocessing

Scientific & Algorithmic Sequence:
  1. Metadata Extraction: Reads physical voxel calibrations (PhysicalSizeX, Y, Z in um) and
     frame interval (TimeIncrement in s) from OME-XML headers, falling back to config defaults.
  2. Fast4DReg 3D Drift Registration (Pylvänäinen et al., 2023, J. Cell Sci.):
     - Estimates lateral (X,Y) and axial (Z) translation shifts strictly on Channel 2
       (tagRFP-T cytosolic cell arbor) via sub-pixel Fourier phase correlation on 2D maximum
       intensity projections.
     - Applies identical shift vectors to Channel 1 (mito-EGFP) in-place, removing sample/stage
       drift without canceling out true intracellular mitochondrial transport.
  3. 3D Auto-Cropping: Computes a padded 3D bounding box around the registered radial glial cell
     arbor, reducing downstream RAM and compute time by 60%–80%.
  4. Cytosolic 3D ROI Gating (Baldwin et al., 2025, JCB; Negri et al., 2025, Front. Cell Dev. Biol.):
     - Generates an exact 3D binary cell mask (*_cell_mask.tif) via Gaussian smoothing,
       Triangle thresholding, largest 3D connected-component filtering, and hole filling.
     - Dilates the cell mask by 1.5 um using an exact 3D Euclidean Distance Transform (O(N)
       complexity) and gates Channel 1 in-place to zero out 100% of extracellular brain
       autofluorescence (lipid droplets, skin pigment, adjacent processes).
  5. Bleach & Depth Attenuation Compensation (Miura, 2020, F1000Research):
     - Matches cumulative intensity histograms of each frame t > 0 to frame t = 0 using fast
       16-bit integer bincount lookup tables (O(N) complexity).
     - Optionally compensates for axial Z-depth scattering attenuation if deep-slice median
       intensity drops by more than the configured trigger ratio (default 30%).
  6. Shot-Noise Despeckling & Local Cytoplasmic Haze Removal (Chaudhry et al., 2020):
     - Applies a 3x3 2D Median filter to remove single-pixel sCMOS/EM-CCD photon shot noise.
     - Applies OpenCV C++-accelerated morphological white top-hat filtering (cv2.MORPH_TOPHAT)
       with an elliptical structuring element (radius = 15 px, ~2.4 um) to strip diffuse
       cytoplasmic GFP haze (~15x faster than scipy.ndimage.white_tophat).
     - Runs on the UNGATED crop, before step 4's gating (see remove_haze_then_gate_inplace).
  6b. Raw-signal QC: measure_raw_mito_contrast() checks the unprocessed mito channel stands
     above background noise; FAIL_NO_SIGNAL exits with code 3 so the Slurm job stops before Stage 4.
  7. Cross-Day Percentile Intensity Normalization:
     - Linearly rescales non-zero foreground intensities from [P_1.0, P_99.9] into a standardized
       16-bit dynamic range [0, 50000] so a single trained ilastik/Cellpose model generalizes
       across imaging sessions taken months or years apart.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Any, Tuple

# Respect SLURM CPU allocation for underlying C/Fortran OpenMP/BLAS threads
_slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK", "4")
os.environ.setdefault("OMP_NUM_THREADS", _slurm_cpus)
os.environ.setdefault("OPENBLAS_NUM_THREADS", _slurm_cpus)
os.environ.setdefault("MKL_NUM_THREADS", _slurm_cpus)

import cv2
import numpy as np
import yaml
from scipy import ndimage

try:
    from tiff_io import read_tiff, write_tiff
    from pipeline_config_io import load_config as _load_pipeline_config
except ImportError:
    from .tiff_io import read_tiff, write_tiff
    from .pipeline_config_io import load_config as _load_pipeline_config


EXIT_QC_FAIL = 3  # run_pipeline_slurm.bash stops before Stage 4 on this code (override: MITO_IGNORE_QC=1)


def load_config(config_path: Path) -> Dict[str, Any]:
    """Loads pipeline_config.yaml, expanding ${MITO_*} storage placeholders (see pipeline_config_io)."""
    return _load_pipeline_config(config_path)


def parse_metadata(meta_dict: Dict[str, Any], defaults: Dict[str, Any]) -> Dict[str, float]:
    """
    Extracts physical pixel dimensions (dx, dy, dz in um) and frame interval (dt in s)
    from parsed OME-XML attributes, falling back to pipeline_config.yaml defaults if absent.
    """
    calib = {
        "dx_um": float(defaults.get("expected_xy_um", 0.16)),
        "dy_um": float(defaults.get("expected_xy_um", 0.16)),
        "dz_um": float(defaults.get("expected_z_um", 1.00)),
        "dt_s": float(defaults.get("expected_dt_s", 11.7)),
    }
    if "PhysicalSizeX" in meta_dict:
        calib["dx_um"] = float(meta_dict["PhysicalSizeX"])
    if "PhysicalSizeY" in meta_dict:
        calib["dy_um"] = float(meta_dict["PhysicalSizeY"])
    if "PhysicalSizeZ" in meta_dict:
        calib["dz_um"] = float(meta_dict["PhysicalSizeZ"])
    if "TimeIncrement" in meta_dict:
        calib["dt_s"] = float(meta_dict["TimeIncrement"])
    return calib


def _structure_only(img: np.ndarray, smooth_sigma, bg_sigma) -> np.ndarray:
    """Light smoothing minus a broad background: keeps cell structure, drops offset/fixed pattern."""
    img = img.astype(np.float32)
    return np.maximum(ndimage.gaussian_filter(img, smooth_sigma) - ndimage.gaussian_filter(img, bg_sigma), 0.0)


def _register_series_2d(frames, window, max_step: float):
    """
    Corrective (x, y) shifts aligning each 2D frame to frame 0 by phase correlation.
    max_step bounds the change in drift between CONSECUTIVE frames (not the total drift, which can
    exceed it many times over in long movies). When the direct frame-0 estimate implies a larger
    step (a likely false lock), the frame-to-frame chained estimate is used for that frame.
    """
    n = len(frames)
    sx = np.zeros(n, dtype=np.float32)
    sy = np.zeros(n, dtype=np.float32)
    for t in range(1, n):
        # cv2.phaseCorrelate(a, b) returns the displacement from a to b; negate it to register b onto a.
        (dx0, dy0), _ = cv2.phaseCorrelate(frames[0], frames[t], window=window)
        (dx1, dy1), _ = cv2.phaseCorrelate(frames[t - 1], frames[t], window=window)
        direct = (-dx0, -dy0)
        chained = (sx[t - 1] - dx1, sy[t - 1] - dy1)
        if max(abs(direct[0] - sx[t - 1]), abs(direct[1] - sy[t - 1])) <= max_step:
            sx[t], sy[t] = direct
        else:
            step_x = float(np.clip(chained[0] - sx[t - 1], -max_step, max_step))
            step_y = float(np.clip(chained[1] - sy[t - 1], -max_step, max_step))
            sx[t], sy[t] = sx[t - 1] + step_x, sy[t - 1] + step_y
    return sx, sy


def estimate_fast4dreg_drift(
    ch_cell: np.ndarray,
    max_shift_xy: int = 25,
    max_shift_z: int = 5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Implements Fast4DReg (Pylvänäinen et al., 2023) projection phase correlation on Channel 2.

    Why projection-based registration:
      Because Bestman Lab spinning-disk stacks have anisotropic voxels (0.16 um XY vs 1.0 um Z),
      solving 2D lateral (X,Y) translation on Z max-projections first and axial (Z) translation
      on Y max-projections second prevents coarse Z steps from corrupting sub-pixel XY alignment,
      while running 4x-70x faster than full 3D FFT cross-correlation.

    Parameters:
      ch_cell      : 4D array of shape (T, Z, Y, X) representing the cytosolic reference channel.
      max_shift_xy : Maximum lateral drift CHANGE between consecutive frames (pixels).
      max_shift_z  : Maximum axial drift CHANGE between consecutive frames (slices).

    Returns:
      Tuple of 1D float32 arrays (shifts_x, shifts_y, shifts_z) representing the corrective
      translation shifts (in pixels) needed to align each frame t to frame t=0.
    """
    n_t, n_z, n_y, n_x = ch_cell.shape
    shifts_x = np.zeros(n_t, dtype=np.float32)
    shifts_y = np.zeros(n_t, dtype=np.float32)
    shifts_z = np.zeros(n_t, dtype=np.float32)

    # Step 1: Lateral (X, Y) drift on background-subtracted Z max-projections.
    # Raw projections sit on a ~100-count camera offset with fixed-pattern noise; phase correlation
    # whitens the spectrum, so that static pattern produces a dominant zero-shift peak and real
    # drift is missed (practice Series #39: raw r = 0.15 vs 0.92 after background subtraction).
    mips_xy = [_structure_only(np.max(ch_cell[t], axis=0), (2.0, 2.0), (25.0, 25.0)) for t in range(n_t)]
    hann_xy = cv2.createHanningWindow((n_x, n_y), cv2.CV_32F)
    sx, sy = _register_series_2d(mips_xy, hann_xy, max_shift_xy)
    shifts_x[:], shifts_y[:] = sx, sy

    # Step 2: Axial (Z) drift on background-subtracted (Z, X) max-projections, after removing the
    # lateral X shift (avoids a 3D volume shift).
    mips_zx = []
    for t in range(n_t):
        zx = np.max(ch_cell[t], axis=1).astype(np.float32)
        if t > 0:
            zx = ndimage.shift(zx, shift=(0.0, float(shifts_x[t])), order=1, mode="nearest")
        # filter along X only: Z stacks are short, and smoothing along Z would blur the axial signal
        mips_zx.append(_structure_only(zx, (0.0, 2.0), (0.0, 25.0)))
    hann_zx = cv2.createHanningWindow((n_x, n_z), cv2.CV_32F)
    _, sz = _register_series_2d(mips_zx, hann_zx, max_shift_z)
    shifts_z[:] = sz

    return shifts_x, shifts_y, shifts_z


def apply_drift_correction_inplace(
    vol_4d: np.ndarray, sx: np.ndarray, sy: np.ndarray, sz: np.ndarray
) -> None:
    """
    Applies corrective translation shifts (sz, sy, sx) in-place across all frames t > 0
    of a 4D stack (T, Z, Y, X), avoiding duplicate 4D array allocation in RAM.
    """
    for t in range(1, vol_4d.shape[0]):
        if abs(sx[t]) > 1e-3 or abs(sy[t]) > 1e-3 or abs(sz[t]) > 1e-3:
            vol_4d[t] = ndimage.shift(
                vol_4d[t],
                shift=(sz[t], sy[t], sx[t]),
                order=1,
                mode="nearest",
            )


def compute_cell_mask_and_crop(
    ch_cell_reg: np.ndarray,
    dx_um: float,
    dy_um: float,
    dz_um: float,
    pad_xy_um: float = 5.0,
    pad_z_um: float = 2.0,
    bg_sigma_xy_um: float = 5.0,
    bg_sigma_z_um: float = 3.0,
) -> Tuple[np.ndarray, Tuple[int, int, int, int, int, int]]:
    """
    Computes an exact 3D binary cell mask from the time-averaged cytosolic channel and
    returns (cropped_mask_3d, bbox_tuple) where bbox is (z_min, z_max, y_min, y_max, x_min, x_max).
    """
    # Time-average suppresses transient noise and highlights the stationary radial glial arbor
    time_avg = np.mean(ch_cell_reg, axis=0)  # Shape: (Z, Y, X)
    smoothed = ndimage.gaussian_filter(time_avg, sigma=(1.0, 2.0, 2.0))
    # Subtract the broad out-of-focus glow before thresholding. Without this, Triangle on the raw
    # average lands just above background: the mask swallowed the soma's halo (practice #39) and
    # merged a NEIGHBORING cell's soma into the target cell (practice #20).
    if bg_sigma_xy_um > 0:
        broad = ndimage.gaussian_filter(
            time_avg, sigma=(bg_sigma_z_um / dz_um, bg_sigma_xy_um / dy_um, bg_sigma_xy_um / dx_um)
        )
        smoothed = np.maximum(smoothed - broad, 0.0)

    # Rescale to 8-bit solely for Triangle threshold estimation (captures dim radial processes)
    s_min, s_max = float(np.min(smoothed)), float(np.max(smoothed))
    s_norm = ((smoothed - s_min) / max(1e-5, (s_max - s_min)) * 255.0).astype(np.uint8)
    th_val, _ = cv2.threshold(
        s_norm.reshape((-1, s_norm.shape[-1])), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_TRIANGLE
    )
    binary = s_norm > max(1, int(th_val))

    # Retain only the largest 3D connected component (the target electroporated radial glial cell)
    labeled, num_cc = ndimage.label(binary)
    if num_cc > 1:
        counts = np.bincount(labeled.ravel())
        counts[0] = 0  # Ignore background label 0
        largest_label = int(np.argmax(counts))
        binary = labeled == largest_label

    # Fill internal vacuoles/nuclear holes inside the soma
    binary = ndimage.binary_fill_holes(binary)

    coords = np.argwhere(binary)
    if len(coords) == 0:
        z_min, z_max = 0, ch_cell_reg.shape[1]
        y_min, y_max = 0, ch_cell_reg.shape[2]
        x_min, x_max = 0, ch_cell_reg.shape[3]
    else:
        pad_z = max(1, int(round(pad_z_um / dz_um)))
        pad_y = max(2, int(round(pad_xy_um / dy_um)))
        pad_x = max(2, int(round(pad_xy_um / dx_um)))

        z_min = max(0, int(np.min(coords[:, 0])) - pad_z)
        z_max = min(ch_cell_reg.shape[1], int(np.max(coords[:, 0])) + pad_z + 1)
        y_min = max(0, int(np.min(coords[:, 1])) - pad_y)
        y_max = min(ch_cell_reg.shape[2], int(np.max(coords[:, 1])) + pad_y + 1)
        x_min = max(0, int(np.min(coords[:, 2])) - pad_x)
        x_max = min(ch_cell_reg.shape[3], int(np.max(coords[:, 2])) + pad_x + 1)

    bbox = (z_min, z_max, y_min, y_max, x_min, x_max)
    cropped_mask = binary[z_min:z_max, y_min:y_max, x_min:x_max]
    return cropped_mask, bbox


def match_histogram_fast_16bit(image: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """
    Fast O(N) cumulative histogram matching for 16-bit fluorescence volumes (Miura, 2020).
    Uses integer bincount across [0, 65535] instead of O(N log N) np.unique sorting.
    """
    img_u16 = np.clip(image, 0, 65535).astype(np.uint16)
    ref_u16 = np.clip(reference, 0, 65535).astype(np.uint16)

    # Compute empirical CDFs via fast 65536-bin integer histograms
    hist_img = np.bincount(img_u16.ravel(), minlength=65536)
    cdf_img = np.cumsum(hist_img).astype(np.float64)
    cdf_img /= max(1.0, cdf_img[-1])

    hist_ref = np.bincount(ref_u16.ravel(), minlength=65536)
    cdf_ref = np.cumsum(hist_ref).astype(np.float64)
    cdf_ref /= max(1.0, cdf_ref[-1])

    # Build 65536-entry lookup table mapping input intensity -> reference intensity.
    # Interpolate only over intensities PRESENT in the reference: the full CDF is flat (=1.0) from
    # the reference maximum up to 65535, and np.interp on that plateau returns 65535, which used to
    # saturate each image's brightest voxels (and whole near-empty Z-slices, see synthetic seeds 2/3).
    present = np.flatnonzero(hist_ref)
    if present.size == 0:
        return np.zeros_like(img_u16, dtype=np.float32)
    lut = np.interp(cdf_img, cdf_ref[present], present.astype(np.float32))
    lut[0] = 0.0  # Keep zero-masked background strictly at 0
    return lut[img_u16]


def compensate_z_attenuation_inplace(
    vol_4d: np.ndarray,
    trigger_drop_ratio: float = 0.30,
    cell_4d: np.ndarray = None,
    cell_mask: np.ndarray = None,
    max_gain: float = 3.0,
) -> Dict[str, Any]:
    """
    Depth (Z) attenuation correction driven by the CYTOSOLIC cell channel, applied as a smooth gain.

    Why not histogram-match each slice to a bright one (the previous method): with sparse organelles a
    deep slice can be dim simply because it contains no mitochondria. Matching it to a bright slice
    stretched its noise to mitochondria-level brightness (practice Series #39: ~1,200 noise objects
    per frame). The cytosolic fill (channel 2) fills the cell at every depth, so its in-cell signal
    above background measures attenuation independently of where mitochondria are.

    Method: per Z-slice with enough cell (>= 0.2% of the slice), signal = median(cell channel inside the
    mask, time-averaged) - median(outside the mask). Fit log(signal) = a + b*z over valid slices. If
    the fitted drop across the cell's Z-range exceeds trigger_drop_ratio, multiply every mito slice by
    gain(z) = exp(-b*(z - z_ref)) relative to the brightest end, clipped to [1, max_gain].
    Returns a dict for provenance (applied, fitted_drop, gains range). Without cell_4d/cell_mask
    nothing is applied (the organelle channel alone cannot distinguish attenuation from content).
    """
    info: Dict[str, Any] = {"applied": False, "fitted_drop": None, "max_gain_applied": 1.0}
    if cell_4d is None or cell_mask is None:
        info["reason"] = "no cell channel/mask supplied"
        return info
    n_t, n_z = vol_4d.shape[:2]
    if n_z < 6:
        info["reason"] = "fewer than 6 slices"
        return info
    cell_avg = cell_4d.mean(axis=0)
    slice_vox = cell_mask.shape[1] * cell_mask.shape[2]
    min_vox = max(50, int(0.002 * slice_vox))
    zs, sig = [], []
    for z in range(n_z):
        inside = cell_mask[z]
        if inside.sum() < min_vox or (~inside).sum() < min_vox:
            continue
        s_val = float(np.median(cell_avg[z][inside]) - np.median(cell_avg[z][~inside]))
        if s_val > 0:
            zs.append(z)
            sig.append(s_val)
    if len(zs) < 6:
        info["reason"] = "fewer than 6 slices with enough cell"
        return info
    zs_a = np.asarray(zs, dtype=np.float64)
    slope, intercept = np.polyfit(zs_a, np.log(np.asarray(sig)), 1)
    fit = np.exp(intercept + slope * zs_a)
    drop = float(1.0 - fit.min() / fit.max())
    info["fitted_drop"] = round(drop, 3)
    if drop < trigger_drop_ratio:
        info["reason"] = f"fitted drop {drop:.2f} < trigger {trigger_drop_ratio}"
        return info
    z_ref = zs_a[np.argmax(fit)]
    gains = np.clip(np.exp(-slope * (np.arange(n_z) - z_ref)), 1.0, max_gain).astype(np.float32)
    vol_4d *= gains[np.newaxis, :, np.newaxis, np.newaxis]
    info.update(applied=True, max_gain_applied=round(float(gains.max()), 3))
    return info


def measure_raw_mito_contrast(
    ch_mito: np.ndarray,
    cell_mask: np.ndarray,
    calib: Dict[str, float],
    ring_um: Tuple[float, float] = (3.0, 8.0),
) -> Dict[str, float]:
    """
    Raw-signal QC on the UNPROCESSED mito channel (T, Z, Y, X), measured before any haze removal
    or normalization (both of which stretch pure noise into apparent structure):
      signal = P99.9 of the time-averaged, lightly smoothed stack inside the cell
               minus its median in a background ring 3-8 um outside the cell;
      noise  = median per-pixel temporal SD in that ring (spatial SD if T == 1).
    contrast_sigma = signal / noise is how far mitochondria stand above single-frame noise.
    """
    sampling = (calib["dz_um"], calib["dy_um"], calib["dx_um"])
    dist_out = ndimage.distance_transform_edt(~cell_mask, sampling=sampling)
    ring = (dist_out > ring_um[0]) & (dist_out < ring_um[1])
    if not ring.any():
        ring = ~cell_mask
    if not ring.any() or not cell_mask.any():
        return {"signal_counts": 0.0, "noise_sd": 0.0, "contrast_sigma": 0.0}
    tavg = ndimage.gaussian_filter(ch_mito.mean(axis=0), sigma=(0.5, 1.0, 1.0))
    signal = float(np.percentile(tavg[cell_mask], 99.9) - np.median(tavg[ring]))
    if ch_mito.shape[0] > 1:
        noise = float(np.median(ch_mito[:, ring].std(axis=0)))
    else:
        noise = float(ch_mito[0][ring].std())
    return {
        "signal_counts": round(signal, 2),
        "noise_sd": round(noise, 3),
        "contrast_sigma": round(signal / max(noise, 1e-6), 2),
    }


def remove_haze_then_gate_inplace(ch_mito: np.ndarray, gate_mask: np.ndarray, radius_px: int) -> None:
    """
    3x3 median despeckle + white top-hat on the UNGATED stack, then zero everything outside the
    dilated cell mask. Order matters: gating first creates a hard 0-edge, and a top-hat subtracts
    nothing from any bright region narrower than its kernel that is surrounded by zeros, so thin
    processes would keep 100% of their haze and the mask rim would turn into bright artifacts.
    """
    ksize = 2 * max(3, int(radius_px)) + 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (ksize, ksize))
    for t in range(ch_mito.shape[0]):
        for z in range(ch_mito.shape[1]):
            despeckled = cv2.medianBlur(np.ascontiguousarray(ch_mito[t, z], dtype=np.float32), 3)
            ch_mito[t, z] = cv2.morphologyEx(despeckled, cv2.MORPH_TOPHAT, kernel)
    ch_mito *= gate_mask[np.newaxis, :, :, :]


def preprocess_movie(
    input_path: Path, work_dir: Path, cfg: Dict[str, Any], dry_run: bool = False
) -> Dict[str, Any]:
    """Executes the complete Stage 3 adaptive preprocessing sequence on a single 4D movie."""
    meta_defaults = cfg.get("metadata_defaults", {})
    p_cfg = cfg.get("stage3_preprocessing", {})

    work_dir.mkdir(parents=True, exist_ok=True)
    movie_id = input_path.stem.replace(".ome", "")

    print(f"[{movie_id}] Reading 4D OME-TIFF: {input_path}")
    # Normalize to (T, C, Z, Y, X) using the file's own axis labels (OME DimensionOrder / ImageJ axes)
    # rather than guessing from array shape, which breaks for XYZTC exports or Z=2 stacks.
    raw_data, meta_tags = read_tiff(input_path, axes="TCZYX")
    calib = parse_metadata(meta_tags, meta_defaults)
    calib_warnings = list(meta_tags.get("_warnings", []))
    for msg in calib_warnings:
        print(f"[{movie_id}] [CALIBRATION WARNING] {msg}", file=sys.stderr)

    ch_mito_idx = int(meta_defaults.get("channel_mito", 1)) - 1
    ch_cell_idx = int(meta_defaults.get("channel_cell", 2)) - 1
    n_channels = raw_data.shape[1]
    if max(ch_mito_idx, ch_cell_idx) >= n_channels:
        raise ValueError(
            f"{movie_id}: needs channels {ch_mito_idx + 1} (mito) and {ch_cell_idx + 1} (cell) "
            f"but the file has C={n_channels} (shape TCZYX={raw_data.shape})."
        )

    if dry_run:
        print(f"[DRY-RUN] Validated {movie_id} shape {raw_data.shape} and calibration {calib}")
        return {"movie_id": movie_id, "stage3_preprocessing_metrics": {"qc_status": "DRY_RUN"}}

    ch_mito = raw_data[:, ch_mito_idx].astype(np.float32, copy=True)  # Shape: (T, Z, Y, X)
    ch_cell = raw_data[:, ch_cell_idx].astype(np.float32, copy=True)  # Shape: (T, Z, Y, X)
    del raw_data  # Free raw buffer immediately to save RAM on compute node
    n_t, n_z, n_y, n_x = ch_mito.shape

    # 1. Fast4DReg 3D Drift Registration
    drift_cfg = p_cfg.get("drift_correction", {})
    print(f"[{movie_id}] Running Fast4DReg 3D drift registration on Channel 2...")
    sx, sy, sz = estimate_fast4dreg_drift(
        ch_cell,
        max_shift_xy=int(drift_cfg.get("max_shift_xy_px", 25)),
        max_shift_z=int(drift_cfg.get("max_shift_z_px", 5)),
    )

    drift_csv = work_dir / f"{movie_id}_drift_shifts.csv"
    with open(drift_csv, "w", encoding="utf-8") as f:
        f.write("timepoint,dx_px,dy_px,dz_px,dx_um,dy_um,dz_um,total_drift_um\n")
        for t in range(n_t):
            tot_um = float(
                np.sqrt(
                    (sx[t] * calib["dx_um"]) ** 2
                    + (sy[t] * calib["dy_um"]) ** 2
                    + (sz[t] * calib["dz_um"]) ** 2
                )
            )
            f.write(
                f"{t},{sx[t]:.3f},{sy[t]:.3f},{sz[t]:.3f},"
                f"{sx[t]*calib['dx_um']:.3f},{sy[t]*calib['dy_um']:.3f},"
                f"{sz[t]*calib['dz_um']:.3f},{tot_um:.3f}\n"
            )

    apply_drift_correction_inplace(ch_cell, sx, sy, sz)
    apply_drift_correction_inplace(ch_mito, sx, sy, sz)
    max_drift_um = float(
        np.max(
            np.sqrt(
                (sx * calib["dx_um"]) ** 2
                + (sy * calib["dy_um"]) ** 2
                + (sz * calib["dz_um"]) ** 2
            )
        )
    )

    # 2. 3D Auto-Cropping & Exact Cell Mask Generation
    crop_cfg = p_cfg.get("auto_crop_3d", {})
    print(f"[{movie_id}] Computing 3D cell mask and bounding box...")
    cell_mask_exact, bbox = compute_cell_mask_and_crop(
        ch_cell,
        dx_um=calib["dx_um"],
        dy_um=calib["dy_um"],
        dz_um=calib["dz_um"],
        pad_xy_um=float(crop_cfg.get("padding_xy_um", 5.0)),
        pad_z_um=float(crop_cfg.get("padding_z_um", 2.0)),
        bg_sigma_xy_um=float(p_cfg.get("cytosolic_roi_gate", {}).get("background_sigma_xy_um", 5.0)),
        bg_sigma_z_um=float(p_cfg.get("cytosolic_roi_gate", {}).get("background_sigma_z_um", 3.0)),
    )
    z_min, z_max, y_min, y_max, x_min, x_max = bbox

    ch_cell_cropped = np.ascontiguousarray(ch_cell[:, z_min:z_max, y_min:y_max, x_min:x_max])
    ch_mito_cropped = np.ascontiguousarray(ch_mito[:, z_min:z_max, y_min:y_max, x_min:x_max])
    del ch_cell, ch_mito

    cell_mask_path = work_dir / f"{movie_id}_cell_mask.tif"
    write_tiff(cell_mask_path, (cell_mask_exact.astype(np.uint8) * 255))

    # 3. Cytosolic 3D ROI Gating via Fast O(N) Physical Euclidean Distance Transform
    gate_cfg = p_cfg.get("cytosolic_roi_gate", {})
    dil_rad_um = float(gate_cfg.get("dilation_radius_um", 1.5))
    # Compute exact 3D distance in microns outside cell_mask_exact
    dist_outside_um = ndimage.distance_transform_edt(
        ~cell_mask_exact, sampling=(calib["dz_um"], calib["dy_um"], calib["dx_um"])
    )
    dilated_mask = dist_outside_um <= dil_rad_um

    # Raw-signal QC before any processing (haze removal + normalization make noise look like signal)
    raw_contrast = measure_raw_mito_contrast(ch_mito_cropped, cell_mask_exact, calib)
    print(
        f"[{movie_id}] Raw mito-channel contrast: {raw_contrast['signal_counts']} counts over "
        f"noise SD {raw_contrast['noise_sd']} = {raw_contrast['contrast_sigma']} sigma"
    )

    # 3b. Despeckle + top-hat haze removal on the UNGATED stack, THEN zero extracellular voxels
    print(f"[{movie_id}] Applying 2D median filter and OpenCV morphological top-hat haze removal...")
    denoise_cfg = p_cfg.get("denoising", {})
    remove_haze_then_gate_inplace(
        ch_mito_cropped, dilated_mask, int(round(float(denoise_cfg.get("rolling_ball_radius_px", 15.0))))
    )

    # 4. Axial Z-Attenuation & Temporal Bleach Correction (Fast 16-bit LUT Histogram Matching)
    z_cfg = p_cfg.get("z_attenuation_correction", {})
    z_triggered = False
    z_info = None
    if z_cfg.get("enabled", True):
        z_info = compensate_z_attenuation_inplace(
            ch_mito_cropped,
            float(z_cfg.get("trigger_drop_ratio", 0.30)),
            cell_4d=ch_cell_cropped,
            cell_mask=cell_mask_exact,
            max_gain=float(z_cfg.get("max_gain", 3.0)),
        )
        z_triggered = bool(z_info["applied"])
        print(f"[{movie_id}] Z-attenuation (cell-channel fit): {z_info}")

    print(f"[{movie_id}] Applying fast 16-bit LUT Histogram Matching bleach correction across T...")
    ref_t0 = ch_mito_cropped[0].copy()
    for t in range(1, n_t):
        ch_mito_cropped[t] = match_histogram_fast_16bit(ch_mito_cropped[t], ref_t0)

    # 5. (Despeckle + top-hat now run in step 3b, before gating.)

    # 6. Cross-Day Percentile Normalization
    norm_cfg = p_cfg.get("intensity_normalization", {})
    non_zero = ch_mito_cropped[ch_mito_cropped > 0]
    if len(non_zero) > 0:
        p_low, p_median, p_high = np.percentile(
            non_zero,
            [
                float(norm_cfg.get("low_percentile", 1.0)),
                50.0,
                float(norm_cfg.get("high_percentile", 99.9)),
            ],
        )
    else:
        p_low, p_median, p_high = 0.0, 1.0, 1.0

    snr = round(float(p_high) / max(1.0, float(p_median)), 2)
    target_max = float(norm_cfg.get("target_max_uint16", 50000))

    # In-place normalization to conserve memory
    np.subtract(ch_mito_cropped, float(p_low), out=ch_mito_cropped)
    np.divide(ch_mito_cropped, max(1e-5, float(p_high - p_low)), out=ch_mito_cropped)
    np.clip(ch_mito_cropped, 0.0, 1.0, out=ch_mito_cropped)
    np.multiply(ch_mito_cropped, target_max, out=ch_mito_cropped)

    # Re-apply dilated mask so background stays strictly 0
    ch_mito_cropped *= dilated_mask[np.newaxis, :, :, :]

    ch_mito_norm = ch_mito_cropped.astype(np.uint16)
    ch_cell_norm = np.clip(ch_cell_cropped, 0, 65535).astype(np.uint16)

    mito_norm_path = work_dir / f"{movie_id}_mito_norm.tif"
    cell_reg_path = work_dir / f"{movie_id}_cell_RFP_reg.tif"

    print(f"[{movie_id}] Exporting normalized stacks: {mito_norm_path.name}, {cell_reg_path.name}")
    write_tiff(mito_norm_path, ch_mito_norm, metadata=calib)
    write_tiff(cell_reg_path, ch_cell_norm, metadata=calib)

    qc_pass = (snr >= float(norm_cfg.get("min_snr_gate", 2.5))) and (
        max_drift_um <= float(drift_cfg.get("qc_max_total_drift_um", 15.0))
    )
    min_contrast = float(p_cfg.get("raw_signal_qc", {}).get("min_contrast_sigma", 3.0))
    if raw_contrast["contrast_sigma"] < min_contrast:
        qc_status = "FAIL_NO_SIGNAL"
    else:
        qc_status = "PASS" if qc_pass else "QC_WARNING"

    metrics = {
        "movie_id": movie_id,
        "input_path": str(input_path),
        "calibration": calib,
        # Non-empty means some calibration came from pipeline_config.yaml defaults, not the file.
        "calibration_warnings": calib_warnings,
        "stage3_preprocessing_metrics": {
            "drift_engine": "Fast4DReg",
            "crop_origin_xyz_px": [int(x_min), int(y_min), int(z_min)],
            "cropped_dimensions_xyz_px": [
                int(x_max - x_min),
                int(y_max - y_min),
                int(z_max - z_min),
            ],
            "max_cumulative_drift_um": round(max_drift_um, 3),
            "z_attenuation_compensated": z_triggered,
            "z_attenuation_fitted_drop": (z_info or {}).get("fitted_drop") if z_cfg.get("enabled", True) else None,
            "pre_norm_p1_0": round(float(p_low), 1),
            "pre_norm_p99_9": round(float(p_high), 1),
            "estimated_snr": snr,
            "raw_mito_contrast": raw_contrast,
            "raw_mito_min_contrast_sigma": min_contrast,
            "qc_status": qc_status,
        },
    }
    provenance_json = work_dir / f"{movie_id}_provenance.json"
    with open(provenance_json, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Module 1: Stage 3 Adaptive 4D Image Preprocessing")
    parser.add_argument("--config", type=Path, required=True, help="Path to pipeline_config.yaml")
    parser.add_argument("--input", type=Path, required=True, help="Path to input 4D .ome.tif")
    parser.add_argument("--work-dir", type=Path, required=True, help="Scratch working directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without writing stacks")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"[ERROR] Input file does not exist: {args.input}", file=sys.stderr)
        return 1

    cfg = load_config(args.config)
    res = preprocess_movie(args.input, args.work_dir, cfg, dry_run=args.dry_run)
    status = res["stage3_preprocessing_metrics"]["qc_status"]
    if status == "FAIL_NO_SIGNAL":
        # Outputs are written (for inspection), but segmenting this movie would only segment noise.
        print(
            f"[QC FAIL] {res['movie_id']}: no detectable mito-channel signal "
            f"({res['stage3_preprocessing_metrics']['raw_mito_contrast']}). Stopping before Stage 4.",
            file=sys.stderr,
        )
        return EXIT_QC_FAIL
    print(f"[SUCCESS] Preprocessing completed for {res['movie_id']} (QC: {status})")
    return 0


if __name__ == "__main__":
    sys.exit(main())