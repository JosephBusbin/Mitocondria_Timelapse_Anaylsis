# Agentic Coding & Parameter Specification: 4D Mitochondrial HPC Pipeline

## 1. Purpose & Agentic Execution Contract

This specification is structured as a deterministic, machine-actionable blueprint for an AI coding agent (or software engineer) to generate, test, and deploy the complete HPC image-processing, segmentation, tracking, and arbor-mapping codebase without ambiguity.

* **Target Compute Environment**: William & Mary SciClone HPC Cluster (RedHat/Rocky Linux, SLURM scheduler, Apptainer/Singularity containers, NVIDIA GPU nodes for Stage 4B, multi-core CPU nodes for Stages 3, 4A, and 5).
* **Strict Data Immutability**: Scripts must NEVER modify, overwrite, or delete input `.sld` or raw `.ome.tif` files in `/sciclone/data10/`. All intermediate tensors must be written to `/sciclone/scr30/$USER/work/` and logged in a per-movie JSON manifest (`<movie_id>_provenance.json`).
* **Fail-Safe & Idempotency**: Every script must support `--dry-run`, `--overwrite`, and exit codes (`0` = success, `1` = input error, `2` = calibration missing, `3` = quality gate failure).

## 2. Master Configuration Schema (`pipeline_config.yaml`)

All generated scripts must read their numerical parameters from a single YAML/JSON configuration file matching the exact schema and default values below:

```yaml
pipeline_version: "2.2.0"

paths:
  raw_ometiff_root: "/sciclone/data10/bestmanlab/raw_ometiff"
  scratch_work_root: "/sciclone/scr30/${USER}/mito_work"
  models_dir: "/sciclone/data10/bestmanlab/models"
  final_output_root: "/sciclone/data10/bestmanlab/outputs"
  fiji_bin: "/opt/Fiji.app/ImageJ-linux64"
  ilastik_bin: "/opt/ilastik/run_ilastik.sh"
  mitolandscape_repo: "/opt/MitoLandscape" # Cloned from https://github.com/enricoenne/MitoLandscape

metadata_defaults:
  channel_mito: 1               # 1-indexed in ImageJ (mito-EGFP)
  channel_cell: 2               # 1-indexed in ImageJ (tagRFP-T cytosol)
  expected_xy_um: 0.16          # Fallback only if OME-XML PhysicalSizeX is missing
  expected_z_um: 1.00           # Fallback only if OME-XML PhysicalSizeZ is missing
  expected_dt_s: 11.7           # Fallback only if OME-XML TimeIncrement is missing

stage3_preprocessing:
  drift_correction:
    enabled: true
    engine: "Fast4DReg"         # Pylvänäinen et al. (2023) J. Cell Sci. (4x-70x faster than Correct 3D drift)
    reference_channel: 2        # Estimate drift strictly on tagRFP-T cell arbor, then apply to Ch1 (mito-EGFP)
    xy_projection_type: "Max_Intensity"
    z_projection_type: "Max_Intensity"
    time_averaging: 1           # Frames to average during cross-correlation
    max_shift_xy_px: 25         # ~4.0 um max frame-to-frame shift in X,Y
    max_shift_z_px: 5           # 5.0 um max frame-to-frame shift in Z
    qc_max_total_drift_um: 15.0 # Flag warning in JSON if cumulative drift > 15 um

  auto_crop_3d:
    enabled: true
    reference_channel: 2
    padding_xy_um: 5.0          # Pad 5 um around 3D cell bounding box in X,Y
    padding_z_um: 2.0           # Pad 2 slices above/below cell in Z

  cytosolic_roi_gate:
    enabled: true               # Isolates single-cell glial territory (Baldwin et al., 2025; Negri et al., 2025)
    gaussian_sigma_xy_px: 2.0
    gaussian_sigma_z_px: 1.0
    threshold_method: "Triangle" # Triangle thresholding captures dim radial shafts
    keep_largest_cc: true       # MorphoLibJ largest 3D connected component -> saved as *_cell_mask.tif
    dilation_radius_um: 1.5     # Dilate by 1.5 um to gate Ch1 (mito-EGFP) extracellular noise

  z_attenuation_correction:
    enabled: true
    trigger_drop_ratio: 0.30    # Apply Z-slice histogram matching if deep Z drops >30%

  bleach_correction:
    enabled: true
    method: "Histogram Matching" # Miura (2020) F1000Research; applied across T

  denoising:
    median_xy_radius_px: 1.0    # Removes 1-px sCMOS/EMCCD shot-noise spikes (Chaudhry et al., 2020)
    rolling_ball_radius_px: 15.0 # ~2.4 um radius; strips diffuse cytoplasmic GFP haze
    sliding_paraboloid: true

  intensity_normalization:
    enabled: true
    low_percentile: 1.0         # Mapped to 0 in uint16
    high_percentile: 99.9       # Mapped to 50000 in uint16 (leaves headroom to 65535)
    target_max_uint16: 50000
    min_snr_gate: 2.5           # Flag low-SNR warning if (P_99.9 / P_50) < 2.5

stage4_hybrid_segmentation:
  stage4a_ilastik:
    enabled: true
    project_file: "bestman_mito_2d_v1.ilp"
    export_source: "Probabilities"
    foreground_class_index: 0   # Channel 0 of exported probabilities = Mitochondria
    feature_mode: "2D_slices"   # Prevents 1.0 um Z-step blur across adjacent planes

  stage4b_instance_segmentation:
    engine: "cellpose"          # Options: "cellpose" (Cellpose 3, Stringer 2025) or "omnipose" (Cutler 2022)
    model_type: "cyto3"         # Or path to fine-tuned "bestman_mito_tube.pth"
    input_mode: "dual_channel"  # Stack [Ch1: mito_norm, Ch2: ilastik_prob]
    diameter_px: 8.0            # ~1.28 um characteristic diameter for thin tubules
    flow_threshold: 0.40        # Higher = more permissive instance splitting
    cellprob_threshold: 0.0     # Logit threshold (0.0 corresponds to P = 0.50)
    do_3d: false                # Disable isotropic 3D; use 2D + Z-stitching for 1.0 um Z
    stitch_threshold: 0.35      # 2D-to-3D mask IoU stitching threshold across Z slices
    anisotropy: 6.25            # PhysicalSizeZ / PhysicalSizeX = 1.0 / 0.16 = 6.25
    min_volume_um3: 0.05        # Hardcoded physical volume floor from MitoLandscape / Mitochondria Analyzer (~2 voxels)

stage5_tracking_and_mitolandscape:
  trackmate_detector:
    type: "LabelImageDetector"  # Reads integer instance masks from Stage 4B
    simplify_contours: true

  sparse_lap_tracker:
    linking_max_distance_um: 2.5
    linking_feature_penalties:
      VOLUME: 1.0               # Penalize linking objects of drastically different volume
    allow_gap_closing: true
    gap_closing_max_distance_um: 3.0
    max_frame_gap: 2
    allow_track_splitting: true # Mitochondrial Fission
    splitting_max_distance_um: 2.0
    allow_track_merging: true   # Mitochondrial Fusion
    merging_max_distance_um: 2.0
    min_track_spots: 4          # Keep trajectories present in >= 4 timepoints

  snt_autotrace:
    enabled: true
    tubeness_scales_um: [0.4, 0.8, 1.2, 1.6]
    soma_detection: "EDT_Threshold"
    max_mito_to_arbor_dist_um: 3.5 # Max radial distance to map a spot onto SNT path

  mitolandscape_integration:
    enabled: true               # Integrates Negri, Fernández, & Borrell (2025) Front. Cell Dev. Biol.
    connectivity_3d: 26         # 26-voxel neighborhood for 3D mitochondrial graph skeletonization
    minor_process_cutoff_um: 5.0 # Processes shorter than 5.0 um merged into soma compartment (id=1)
    process_slice_thickness_um: 0.5 # Geodesic distance bin width along twisted radial process (default: 0.5 um)
    apical_basal_orientation_axis: "Y" # Axis pointing toward ventricular surface for apical (id=2) vs basal (id=3)
    motility_velocity_threshold_um_s: 0.02 # |v_axial| > 0.02 um/s classified as active anterograde/retrograde
```

## 3. Module-by-Module Implementation Specifications

### Module 1: `preprocess_4d_stack.py` (Stage 3 Adaptive Preprocessing)

* **Runtime**: Headless FIJI Jython or Python (`tifffile`, `scipy.ndimage`, `scikit-image`, `SimpleITK` / `pyimagej`).
* **Input**: Single 4D, 2-channel `.ome.tif` from `collect_4d_sld_timelapses.py`.
* **Exact Algorithmic Sequence**:
  1. **Read Metadata**: Extract `PhysicalSizeX`, `PhysicalSizeY`, `PhysicalSizeZ`, and `TimeIncrement` from OME-XML tag. Abort with exit code `2` if missing and fallback is disabled.
  2. **3D Drift Registration via `Fast4DReg` ([Pylvänäinen et al., 2023](https://journals.biologists.com/jcs/article/136/4/jcs260728/287682/Fast4DReg-fast-registration-of-4D-microscopy))**:
     * Split channels into `imp_mito` (Ch1) and `imp_cell` (Ch2).
     * Run `Fast4DReg` (`Estimate drift`) on `imp_cell` (Ch2 `tagRFP-T`) with `XY_registration=true`, `Z_registration=true`, `projection_type=[Max Intensity]`, exporting `<movie_id>_drift_table.csv` and `<movie_id>_fast4dreg_settings.csv`.
     * Run `Fast4DReg` (`Apply drift`) on both `imp_cell` and `imp_mito` using `<movie_id>_fast4dreg_settings.csv`.
     * Compute maximum cumulative displacement $\sqrt{\Delta x^2 + \Delta y^2 + \Delta z^2}$ from `<movie_id>_drift_table.csv` and log to `<movie_id>_provenance.json`.
  3. **3D Auto-Crop Around Registered Cell**:
     * Compute maximum projection across $T$ of registered `imp_cell`, smooth with Gaussian ($\sigma_{xy}=2, \sigma_z=1$), apply `Triangle` threshold, and find 3D bounding box $[x_{\min}:x_{\max}, y_{\min}:y_{\max}, z_{\min}:z_{\max}]$.
     * Pad bounding box by `padding_xy_um` and `padding_z_um`, crop both `imp_mito` and `imp_cell`, and store `crop_origin_px` and `crop_origin_um` in `<movie_id>_provenance.json`.
  4. **Cytosolic 3D Masking ([Baldwin et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11527035/); [Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full))**:
     * From the time-averaged projection of registered `imp_cell` (or per-frame if process remodels), apply 3D Gaussian filter ($\sigma_{xy}=2.0, \sigma_z=1.0$), binarize via `Triangle` threshold, and retain the largest 3D connected component (`scipy.ndimage.label` / `MorphoLibJ`) as `mask_cell_exact` (saved as `<movie_id>_cell_mask.tif` for Module 4).
     * Dilate `mask_cell_exact` by $r = \lceil 1.5 / \Delta xy \rceil$ pixels in $X,Y$ and $1$ slice in $Z$ to create `mask_cell_dilated`.
     * Multiply `imp_mito[t]` by `mask_cell_dilated` so extracellular autofluorescence is set to `0`.
  5. **Bleach & Axial Attenuation Correction ([Miura, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7871415/))**:
     * Run FIJI `Bleach Correction` with `method=[Histogram Matching]` across $T$ on `imp_mito`.
  6. **Shot-Noise & Haze Removal ([Chaudhry et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7052579/))**:
     * Run 2D Median filter (`radius=1 stack`) followed by `Subtract Background...` (`rolling=15 sliding disable stack`).
  7. **Cross-Day Percentile Normalization**:
     * Compute $P_{1.0}$ and $P_{99.9}$ across non-zero masked voxels of `imp_mito`.
     * Apply linear stretch:
       $$I_{\text{norm}}(x,y,z,t) = \text{clip}\left(\frac{I(x,y,z,t) - P_{1.0}}{P_{99.9} - P_{1.0}}, 0, 1\right) \times 50000$$
     * Cast to `uint16` and save as `<movie_id>_mito_norm.tif` and `<movie_id>_cell_RFP_reg.tif`.

### Module 2: `run_hybrid_segmentation.py` (Stage 4A ilastik + Stage 4B Cellpose/Omnipose)

* **Runtime**: Python 3.10+ on HPC GPU node (`cellpose>=3.0`, `tifffile`, `numpy`, plus `ilastik` binary).
* **Exact Algorithmic Sequence**:
  1. **Stage 4A (Headless ilastik)**:
     * Execute subprocess:
       ```bash
       /opt/ilastik/run_ilastik.sh --headless \
           --project="${MODELS_DIR}/bestman_mito_2d_v1.ilp" \
           --export_source="Probabilities" \
           --output_format=tif \
           --output_filename_format="${WORK_DIR}/{nickname}_Probabilities.tif" \
           "${WORK_DIR}/${MOVIE_ID}_mito_norm.tif"
       ```
     * Extract Channel 0 (Mitochondrial foreground probability $P \in [0.0, 1.0]$).
  2. **Stage 4B (Cellpose 3 / Omnipose 3D Instance Splitting)**:
     * Construct a 2-channel array for each timepoint $t$ of shape `(Z, 2, Y, X)`:
       * Channel 0: `mito_norm[t]` (rescaled $0.0–1.0$)
       * Channel 1: `ilastik_prob[t]` ($0.0–1.0$, zeroing out voxels where $P < 0.25$)
     * Run Cellpose 3 model per timepoint $t$ with **2D-to-3D anisotropic stitching**:
       ```python
       from cellpose import models
       model = models.CellposeModel(gpu=True, model_type=cfg["model_type"])
       masks_t, flows, styles = model.eval(
           img_t,
           diameter=8.0,
           channels=[1, 2],
           flow_threshold=0.40,
           cellprob_threshold=0.0,
           do_3d=False,
           stitch_threshold=0.35, # Links 2D masks across 1.0 um Z slices via IoU
       )
       ```
     * Filter out 3D connected instances whose physical volume $V = N_{\text{voxels}} \times \Delta x \times \Delta y \times \Delta z < 0.05\ \mu\text{m}^3$ ([matching the volume floor validated in Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)).
     * Re-stack all timepoints $t=0 \dots T-1$ into a calibrated 4D `uint16` OME-TIFF `<movie_id>_mito_labels.tif` preserving physical voxel metadata $(\Delta x, \Delta y, \Delta z, \Delta t)$.

### Module 3: `run_trackmate_sparse_lap.py` (Stage 5A Headless TrackMate 7)

* **Runtime**: Headless FIJI (`ImageJ-linux64 --headless --console --run run_trackmate_sparse_lap.py`).
* **Input**: Composite 2-channel 4D hyperstack where **Channel 1** = `<movie_id>_mito_labels.tif` (instance IDs) and **Channel 2** = `<movie_id>_mito_norm.tif` (calibrated intensities).
* **Detector & Tracker Contract**:
  * **Detector**: `LabelImageDetectorFactory` targeting `TARGET_CHANNEL = 1`, `SIMPLIFY_CONTOURS = True`.
  * **Spot Analyzers**: Register `SpotIntensityMultiCAnalyzerFactory` and `SpotMorphologyAnalyzerFactory`. Because Channel 1 contains the integer label IDs ($1 \dots N_t$), `MEAN_INTENSITY_CH1` in `spots.csv` equals the exact integer `label_id` of that mitochondrion in `<movie_id>_mito_labels.tif`, while `MEAN_INTENSITY_CH2` records its normalized `mito-EGFP` fluorescence!
  * **Tracker**: `SparseLAPTrackerFactory` parameterized strictly from `pipeline_config.yaml` (`LINKING_MAX_DISTANCE = 2.5`, `GAP_CLOSING_MAX_DISTANCE = 3.0`, `MAX_FRAME_GAP = 2`, `ALLOW_TRACK_SPLITTING = True`, `SPLITTING_MAX_DISTANCE = 2.0`, `ALLOW_TRACK_MERGING = True`, `MERGING_MAX_DISTANCE = 2.0`).
  * **Outputs**: Write `<movie_id>_tracks.xml`, `<movie_id>_spots.csv`, `<movie_id>_edges.csv`, and `<movie_id>_tracks.csv`.

### Module 4: `run_arbor_and_mitolandscape.py` (Stage 5B SNT + MitoLandscape 4D Profiling)

* **Runtime**: Python 3.10+ (`numpy`, `pandas`, `scipy`, `scikit-image`, `scikit-fmm`, `networkx`, `tifffile`) importing/adapting functions from [`enricoenne/MitoLandscape`](https://github.com/enricoenne/MitoLandscape) ([Negri, Fernández, & Borrell, 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)).
* **Inputs**:
  * `<movie_id>_cell_mask.tif` (3D binary mask of radial glial cell from Module 1)
  * `<movie_id>_mito_labels.tif` (4D integer instance masks of mitochondria from Module 2)
  * `<movie_id>_spots.csv` and `<movie_id>_edges.csv` (TrackMate 4D trajectories from Module 3)
  * `<movie_id>_arbor.swc` (3D centerline + soma root from SNT Auto-Trace)
* **Exact Algorithmic & Mathematical Sequence**:

  1. **Automated Soma & Compartment Annotation (Adapting MitoLandscape `cell_annotation` & `apical_progenitor_compartment_annotation`)**:
     * Identify the 3D soma mask (`soma_mask`) using Euclidean Distance Transform (`scipy.ndimage.distance_transform_edt` weighted by `dz, dy, dx`) inside `cell_mask`, thresholding at the soma neck radius (or reading the SNT Soma ROI).
     * Subtract `soma_mask` from `cell_mask` to isolate detached process components (`scipy.ndimage.label`).
     * Reassign minor processes with geodesic length $< 5.0\ \mu\text{m}$ (`minor_process_cutoff_um`) back to `soma_mask` (`compartment_id = 1`).
     * Classify remaining primary processes based on their starting attachment vector relative to the soma and ventricular axis (`apical_progenitor_compartment_annotation`):
       * `compartment_id = 1`: **Soma**
       * `compartment_id = 2`: **Apical Process**
       * `compartment_id = 3`: **Basal Process / Radial Shaft**
       * `compartment_id = 4`: **Distal Endfeet / Terminal Arbor** (Strahler order $= 1$ distal to main shaft)

  2. **3D Internal Geodesic Distance Map (MitoLandscape `geodesic_dist` via `scikit-fmm`)**:
     * Construct a speed/phi array where `cell_mask == 1` is the valid domain and `soma_mask == 1` is the zero-contour origin ($\phi = 0$).
     * Execute Fast Marching Method:
       ```python
       import skfmm
       import numpy as np
       phi = np.ones_like(cell_mask, dtype=np.float32)
       phi[soma_mask > 0] = -1.0
       masked_phi = np.ma.MaskedArray(phi, mask=(cell_mask == 0))
       geodesic_map_um = skfmm.distance(masked_phi, dx=(dz_um, dy_um, dx_um))
       ```
     * Every cytoplasmic voxel $(z, y, x)$ now holds its exact 3D geodesic path distance $d_{\text{geo}}(z, y, x)$ in $\mu\text{m}$ from the soma boundary, faithfully following twists and turns of the radial glial process.

  3. **Strahler Branch Order & Loop-Free Cell Graph (MitoLandscape `strahler_analysis` via `networkx`)**:
     * Skeletonize `cell_mask` in 3D (`skimage.morphology.skeletonize_3d`) and build a 26-connected `networkx.Graph` with edge weights equal to physical Euclidean distance.
     * Break any spurious segmentation loops using `nx.minimum_spanning_tree`, weighting edges by $d_{\text{geo}}$ so loops are preferentially opened at the distal tips rather than near the soma.
     * Set the node with minimum $d_{\text{geo}}$ as the root and compute **Strahler numbers** (`strahler_order = 1` for terminal branches, incrementing where two branches of equal order meet). Map `branch_id` and `strahler_order` back to all voxels of `cell_mask` (`annotate_from_network`).

  4. **Per-Instance 3D Mitochondrial Graph Barycenter & Subtype Classification (MitoLandscape `mito_analysis`)**:
     * For each timepoint $t$ and each mitochondrial instance ID $m \in \{1 \dots N_t\}$ in `<movie_id>_mito_labels.tif`:
       * Extract its 3D binary sub-volume, skeletonize in 3D, and build a 26-connected `networkx.Graph` $G_m$.
       * **Graph Barycenter**: Compute all-pairs shortest path distances weighted by $(\Delta z, \Delta y, \Delta x)$ and select the skeleton voxel $\mathbf{b}_{m,t} = (z_b, y_b, x_b)$ that minimizes total distance to all other nodes in $G_m$. (Unlike a Euclidean center-of-mass, the graph barycenter always lies *inside* curved or branched mitochondria.)
       * Count skeleton endpoints ($N_{\text{end}}$, degree $= 1$) and junctions ($N_{\text{junc}}$, degree $\ge 3$), and classify morphology per [Bakare et al. (2021) / Negri et al. (2025)](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full):
         * **`punctus`**: Single node or aspect ratio $< 1.5$ and length $< 1.0\ \mu\text{m}$.
         * **`rod`**: Unbranched tubule ($N_{\text{junc}} = 0$, $N_{\text{end}} = 2$, length $\ge 1.0\ \mu\text{m}$).
         * **`network`**: Branched structure ($N_{\text{junc}} \ge 1$).
       * Lookup `geodesic_dist_um` $= d_{\text{geo}}(\mathbf{b}_{m,t})$, `compartment_id`, `branch_id`, and `strahler_order` at the barycenter coordinate $\mathbf{b}_{m,t}$.

  5. **4D Trajectory Coupling & $0.5\ \mu\text{m}$ Process Slicing (MitoLandscape `process_analysis` + TrackMate Edges)**:
     * Join the per-instance `mito_analysis` features onto TrackMate's `spots.csv` and `edges.csv` via `(FRAME, int(round(MEAN_INTENSITY_CH1))) == (t, label_id)`.
     * For every linked edge $(t \to t+1)$, compute signed geodesic velocity along the process:
       $$v_{\text{geo}} = \frac{d_{\text{geo}}(\mathbf{b}_{m,t+1}) - d_{\text{geo}}(\mathbf{b}_{m,t})}{\Delta t}$$
       * Classify transport direction: **Anterograde / Distal** ($v_{\text{geo}} > +0.02\ \mu\text{m/s}$), **Retrograde / Proximal** ($v_{\text{geo}} < -0.02\ \mu\text{m/s}$), or **Stationary / Docked** ($|v_{\text{geo}}| \le 0.02\ \mu\text{m/s}$).
     * **Longitudinal $0.5\ \mu\text{m}$ Binning (`process_analysis`)**: Partition the radial process into discrete geodesic distance shells $[k \cdot 0.5\ \mu\text{m}, (k+1) \cdot 0.5\ \mu\text{m})$. For each shell $k$ and timepoint $t$, compute:
       * `cytoplasm_volume_um3`: Voxel count in shell $k$ $\times (\Delta x \Delta y \Delta z)$.
       * `mito_volume_um3`: Mitochondrial voxel count in shell $k$ $\times (\Delta x \Delta y \Delta z)$.
       * `mito_occupancy_fraction`: $\text{mito\_volume\_um3} / \text{cytoplasm\_volume\_um3}$ ([compartment occupancy metric per Virga et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10923867/) and [Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)).
       * `fission_count` (out-degree $= 2$) and `fusion_count` (in-degree $= 2$) occurring inside shell $k$.
     * **Outputs**: Write `<movie_id>_snt_mitolandscape_tracks.csv`, `<movie_id>_process_0.5um_profile.csv`, and `<movie_id>_remodeling_events.csv`.

## 4. Per-Movie JSON Provenance Schema (`<movie_id>_provenance.json`)

Every SLURM array task must emit a companion JSON file alongside its CSV outputs with the following fields:

```json
{
  "movie_id": "20240514_Tadpole03_S02_NPC_Timelapse",
  "source_sld_file": "20240514_Tadpole03.sld",
  "source_series_index": 2,
  "calibration": {
    "dx_um": 0.16,
    "dy_um": 0.16,
    "dz_um": 1.0,
    "dt_s": 11.7,
    "dimensions_xyztc": [512, 512, 45, 30, 2]
  },
  "stage3_preprocessing_metrics": {
    "drift_engine": "Fast4DReg",
    "crop_origin_xyz_px": [112, 84, 3],
    "cropped_dimensions_xyz_px": [180, 240, 38],
    "max_cumulative_drift_um": 3.42,
    "pre_norm_p1_0": 268.0,
    "pre_norm_p99_9": 4120.5,
    "estimated_snr": 14.8,
    "qc_status": "PASS"
  },
  "stage4_segmentation_metrics": {
    "ilastik_model": "bestman_mito_2d_v1.ilp",
    "instance_engine": "cellpose_3.0.1",
    "instance_model": "cyto3",
    "min_volume_filter_um3": 0.05,
    "mean_instances_per_frame": 64.2
  },
  "stage5_tracking_and_mitolandscape_summary": {
    "total_tracks_filtered": 58,
    "total_fission_events": 7,
    "total_fusion_events": 5,
    "mean_geodesic_speed_um_s": 0.128,
    "subtype_fractions": {
      "punctus": 0.42,
      "rod": 0.49,
      "network": 0.09
    },
    "compartment_occupancy_pct": {
      "soma": 24.5,
      "apical_process": 18.2,
      "basal_process": 14.7,
      "endfoot": 31.0
    }
  }
}
```