# Reproducible Protocol: ilastik Workflow Evaluation & Pixel Classifier Training

## Part 1: Evaluation of ilastik's Built-In Workflows vs. Our Hybrid Architecture

When opening the [ilastik startup screen](https://www.ilastik.org/documentation/), eight workflows are presented. Below is the technical evaluation of why **Pixel Classification (paired with Stage 4B Cellpose/Omnipose and Stage 5A TrackMate)** is superior to ilastik's standalone object/tracking workflows for Professor Bestman's 4D radial glia data.

| ilastik Workflow | What It Does ([Manual & Docs](https://www.ilastik.org/documentation/)) | Why We Use or Reject It for Bestman Lab Radial Glia Data |
| :--- | :--- | :--- |
| **1. [Pixel Classification](https://www.ilastik.org/documentation/pixelclassification/pixelclassification)** | Trains a multi-scale Random Forest classifier on intensity, edge, and texture features to output a continuous probability map ($P \in [0.0, 1.0]$) for each class. | **SELECTED (Stage 4A Mainstay)** — Excels at separating true mitochondrial fluorescence from deep-tissue *Xenopus* brain scattering and cytoplasmic GFP haze using just 5–10 minutes of sparse brush strokes. |
| **2. [Autocontext](https://www.ilastik.org/documentation/autocontext/autocontext)** | Runs a two-stage cascaded Random Forest where Stage 2 trains on both the raw image features and the probability map from Stage 1. | **OPTIONAL UPGRADE** — Can sharpen organelle boundaries in extremely noisy deep-tissue stacks, but doubles headless CPU runtime and is largely redundant because **Stage 4B (Cellpose 3 / Omnipose)** already takes `[mito_norm, ilastik_prob]` as a 2-channel input. |
| **3. Automatic Tracking (Chaingraph / Conservation Tracking)** | Tracks objects in 2D+t or 3D+t using a conservation-of-mass factor graph solved via Integer Linear Programming (Manual Section 3.4, pp. 34–37). | **REJECTED** — **(1) Requires proprietary IBM CPLEX C++ static libraries** manually compiled and linked into `ilastik/lib` (Manual Section 1.3.1.1, pp. 6–9), which is fragile on SLURM clusters; **(2) Models cell mitosis (1 $\to$ 2 division)** rather than reversible mitochondrial **fusion (2 $\to$ 1) + fission (1 $\to$ 2)**; **(3)** TrackMate 7's Sparse LAP tracker handles both fission and fusion natively in open-source Java/Python. |
| **4. Object Classification** | Extracts connected components via simple or hysteresis thresholding on the probability map and classifies objects by shape/size (Manual Section 3.4.2.1, p. 36). | **REJECTED** — Relies on standard intensity thresholding to define objects. Inside a 1–2 $\mu\text{m}$ radial glial shaft, touching mitochondria share a continuous probability bridge ($P > 0.6$) and merge into a single object. **Cellpose 3 / Omnipose** gradient flows split touching tubules cleanly. |
| **5. Density Counting** | Estimates a continuous Gaussian density integral to count overlapping blobs in 2D boxes without segmenting them (Manual Section 2.4, pp. 20–24). | **REJECTED** — Explicitly **restricted to 2D images only** (Manual p. 20) and outputs regional floating-point counts rather than 3D $(x,y,z)$ coordinates, volumes, or trajectories. |
| **6. Boundary-Based Segmentation with Multicut** | Uses boundary probability maps and graph partitioning to segment closed membrane-stained cells (primarily for 3D Electron Microscopy). | **REJECTED** — Requires clear dark/bright boundary rings around every object; matrix-labeled `mito-EGFP` mitochondria at $0.16\ \mu\text{m}$ pixel size are solid diffraction-limited rods/puncta (~3–6 px wide) without hollow interiors. |
| **7. Carving & Neural Network Classification** | Interactive 3D seeded watershed for single objects (`Carving`) or BioImage.IO PyTorch runner (`Neural Network`). | **REJECTED** — `Carving` is manual per object, and `Stage 4B` already runs PyTorch Cellpose/Omnipose natively with 2D-to-3D anisotropic stitching. |

---

## Part 2: Step-by-Step Reproducible Protocol for Training ilastik (`bestman_mito_2d_v1.ilp`)

### Crucial Rule Before You Start: Train on Stage 3 Normalized Stacks, Not Raw `.sld` Exports!
Because `preprocess_4d_stack.py` (Module 1) strips out extracellular autofluorescence, removes rolling-ball haze, corrects photobleaching, and standardizes every movie's non-zero intensity range to `[0, 50000]` (`uint16`), **you must train ilastik on 2–3 `<movie_id>_mito_norm.tif` files output by Module 1**. Training on standardized stacks guarantees your `.ilp` classifier generalizes across experiments taken years apart.

---

### Step 1: Select & Prepare 3 Representative Training Stacks
1. Run `preprocess_4d_stack.py` on **3 diverse movies** spanning your dataset's conditions:
   * **Stack A**: A high-SNR representative movie from an early cohort (e.g., 2023).
   * **Stack B**: A deeper, higher-scattering movie (e.g., $>100\ \mu\text{m}$ depth or dimmer expression).
   * **Stack C**: A movie with a densely packed radial glial shaft where mitochondria pass closely past one another.
2. *(Recommended for GUI Responsiveness)*: Instead of loading all 30 timepoints ($30 \times 45 = 1,350$ planes per movie) into the interactive GUI, extract 2 representative 3D volumes per movie (e.g., Frame $T=0$ and Frame $T=15$) or load the cropped `<movie_id>_mito_norm.tif` directly.

---

### Step 2: Create the Project & Verify Axis Order (`Data Selection` Applet)
1. Launch the desktop ilastik application, select **Pixel Classification** under *Create New Project*, and save the project as **`bestman_mito_2d_v1.ilp`**.
2. In the **1. Input Data** applet, click **Add New... $\to$ Add separate image(s)...** and select your 3 `<movie_id>_mito_norm.tif` files.
3. **CRITICAL REPRODUCIBILITY CHECK — Verify Axes column** (Manual Section 3.4.1, p. 35, Fig. 3.2):
   * Inspect the **Axes** column in the Raw Data table.
   * For a 4D single-channel time-lapse $(T, Z, Y, X)$, ilastik must display **`tzyx`** (or **`zyx`** if you loaded a single 3D timepoint, with channel `c` implicit or `tzyxc`).
   * If ilastik misinterprets the axes (for example, showing `xyzc` instead of `tzyx`), right-click the row $\to$ **Edit properties...** and change **Axes** to **`tzyx`** (and set **Display Mode** to **Grayscale** so the Window Leveling contrast tool `t` is available).

---

### Step 3: Configure Multi-Scale Spatial Features (`Feature Selection` Applet)
1. Click **2. Feature Selection** on the left panel, then click **Select Features...** (Manual p. 21, Fig. 2.11; p. 50, Fig. 4.2).
2. **Set Computation Mode to `2D` (Checkbox at bottom of dialog: *Compute features in 2D*)**:
   * **Justification**: As noted in the [official ilastik Pixel Classification guide](https://www.ilastik.org/documentation/pixelclassification/pixelclassification), when $Z$-slices are thick ($1.0\ \mu\text{m}$ $z$-step vs. $0.16\ \mu\text{m}$ $xy$ pixel size = $6.25\times$ anisotropy), computing isotropic 3D filters blurs across adjacent $Z$-planes where mitochondria do not even exist. Computing features in **2D per slice** keeps tubule boundaries razor-sharp and matches `pipeline_config.yaml` (`feature_mode: "2D_slices"`).
3. **Check the Following Exact Feature Matrix ($\sigma = 0.3\text{ to }5.0\text{ px}$)**:
   At $0.16\ \mu\text{m/pixel}$, a mitochondrial tubule is ~2–4 pixels wide ($0.32–0.64\ \mu\text{m}$) and 6–25 pixels long ($1–4\ \mu\text{m}$):

| Feature Category | Specific Filter | $\sigma = 0.3$ | $\sigma = 0.7$ | $\sigma = 1.0$ | $\sigma = 1.6$ | $\sigma = 3.5$ | $\sigma = 5.0$ | $\sigma = 10.0$ | Biophysical Justification for Radial Glia Mitochondria |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Color / Intensity** | **Gaussian Smoothing** | `[x]` | `[x]` | `[x]` | `[x]` | `[x]` | `[ ]` | `[ ]` | Captures local normalized brightness of puncta ($\sigma \le 1.0$) and tubules ($\sigma \le 3.5$). |
| **Edge** | **Laplacian of Gaussian (LoG)** | `[ ]` | `[x]` | `[x]` | `[x]` | `[x]` | `[x]` | `[ ]` | **Blob detector** (Manual p. 20): peaks sharply at the center of punctate mitochondria and cross-sections of rods. |
| **Edge** | **Gaussian Gradient Magnitude** | `[ ]` | `[x]` | `[x]` | `[x]` | `[x]` | `[ ]` | `[ ]` | **Boundary detector**: identifies steep intensity drop-offs at the outer mitochondrial membrane. |
| **Edge** | **Difference of Gaussians (DoG)** | `[ ]` | `[x]` | `[x]` | `[x]` | `[x]` | `[x]` | `[ ]` | Band-pass filter that isolates organelle-sized structures ($0.1–0.8\ \mu\text{m}$) from residual background fluctuations. |
| **Texture** | **Structure Tensor Eigenvalues** | `[ ]` | `[x]` | `[x]` | `[x]` | `[x]` | `[x]` | `[ ]` | Quantifies local directional coherence and anisotropy along elongated mitochondrial rods. |
| **Texture** | **Hessian of Gaussian Eigenvalues** | `[ ]` | `[x]` | `[x]` | `[x]` | `[x]` | `[x]` | `[ ]` | **Tubular line/ridge detector** (Manual p. 20): essential for recognizing thin $0.35\ \mu\text{m}$ mitochondrial tubules and keeping narrow gaps between passing mitochondria open. |

*(Leave $\sigma = 10.0\text{ px}$ unchecked: at $10 \times 0.16 = 1.6\ \mu\text{m}$, $\sigma=10.0$ is wider than the entire radial glial shaft, which slows down CPU calculation and blurs adjacent mitochondria together.)*

---

### Step 4: Interactive Brush Annotation & Uncertainty-Guided Training (`Training` Applet)
1. Switch to **3. Training**. Rename the two default labels by double-clicking their names:
   * **Label 1 (Yellow/Red — Index 0)**: Rename to **`Mitochondria`** *(Must be Label 1 so Channel 0 of exported probabilities is foreground, matching `foreground_class_index: 0` in `pipeline_config.yaml`!)*.
   * **Label 2 (Blue/Green — Index 1)**: Rename to **`Background`**.
2. **Use Window Leveling (`t` key) for Display Contrast**:
   * Press **`t`** (or click the Window Leveling icon) and drag the mouse to adjust visual brightness so dim tubules are easy to see. *(Per the ilastik manual, this only adjusts screen display and never alters pixel values.)*
3. **Follow Strict "Sparse Center & Border" Brush Rules (`Brush Size = 1`)**:
   * **Rule A — Never paint across an entire mitochondrion**: Use **`Brush Size = 1`** and draw a thin 1-pixel line strictly down the **bright core/spine** of 5–8 representative mitochondria (include both bright somatic mitochondria and dimmer shaft rods). If you paint fuzzy outer halos as `Mitochondria`, the classifier will artificially bloat organelle volumes and merge neighbors.
   * **Rule B — Paint the "Clefts" Between Closely Spaced Mitochondria**: Zoom in ($400\%$) on the thin radial glial process where two mitochondria sit end-to-end or side-by-side. Select `Background` (`Brush Size = 1`) and place a 1–2 pixel dot/stroke directly in the **dim gap/trough between the two mitochondria**, plus a stroke along the immediate outer edge of the tubule and a broad stroke (`Brush Size = 3`) in the empty black background.
4. **Activate `Live Update` & Use the `Uncertainty` Layer**:
   * Click **Live Update** (Manual p. 22, Fig. 2.12).
   * In the bottom-left **Group Visibility** panel, toggle the eye icon next to **Uncertainty** (cyan/magenta overlay) and **Segmentation** (`I` key toggles overlays on/off).
   * **Active Learning Loop**: Navigate through $Z$-slices and across your 3 training stacks. Wherever the **Uncertainty** layer glows brightly (indicating the Random Forest trees are split ~50/50), turn **Live Update OFF**, add a 1-pixel corrective stroke, and turn **Live Update ON**.
   * **Stop Early (Prevent Overfitting)**: Stop as soon as 15–30 total strokes across the 3 stacks produce clean, separated probability islands. Over-labeling hundreds of pixels slows down tree evaluation without improving accuracy.

---

### Step 5: Configure `Prediction Export` Settings for Headless SciClone HPC Execution
1. Click **4. Prediction Export** on the left panel.
2. Ensure **Source** is set to **`Probabilities`** (the default).
3. Click **Choose Export Image Settings...** (Manual p. 25, Fig. 2.16) and verify:
   * **Convert to Data Type**: `floating 32-bit` (`float32`, range $0.0–1.0$) — *Do not renormalize to 8-bit integer, as Stage 4B expects floating-point probabilities.*
   * **Transpose to Axis Order**: Keep default (`tzyxc` / `zyxc`).
   * **Output File Info $\to$ Format**: Select **`tif`** (or `multipage tiff`).
4. Click **OK**, then save the project (**`Project $\to$ Save Project`** as `bestman_mito_2d_v1.ilp`) and transfer it to `/sciclone/data10/${USER}/bestmanlab/models/bestman_mito_2d_v1.ilp`.

---

## Part 3: Reproducibility Checklist for Your Lab Notebook (Checkpoint 3 Record)

Copy and fill out this block in your project's `README_provenance.md` whenever you train or update an `.ilp` model:

```yaml
ilastik_classifier_provenance:
  project_filename: "bestman_mito_2d_v1.ilp"
  date_trained: "YYYY-MM-DD"
  ilastik_version: "1.4.x" # Must match HPC container/binary version
  training_stacks_used:
    - movie_id: "<movie_1_mito_norm.tif>"
      frames_annotated: [0, 15]
      rationale: "High-SNR baseline radial glial arbor"
    - movie_id: "<movie_2_mito_norm.tif>"
      frames_annotated: [0, 10]
      rationale: "Deep-tissue (>100 um) high-scattering stack"
    - movie_id: "<movie_3_mito_norm.tif>"
      frames_annotated: [5, 20]
      rationale: "Crowded thin radial shaft with passing mitochondria"
  feature_configuration:
    computation_mode: "2D (per Z-slice)"
    gaussian_smoothing_sigmas: [0.3, 0.7, 1.0, 1.6, 3.5]
    laplacian_of_gaussian_sigmas: [0.7, 1.0, 1.6, 3.5, 5.0]
    gaussian_gradient_magnitude_sigmas: [0.7, 1.0, 1.6, 3.5]
    difference_of_gaussians_sigmas: [0.7, 1.0, 1.6, 3.5, 5.0]
    structure_tensor_eigenvalues_sigmas: [0.7, 1.0, 1.6, 3.5, 5.0]
    hessian_of_gaussian_eigenvalues_sigmas: [0.7, 1.0, 1.6, 3.5, 5.0]
  class_definitions:
    label_1_channel_0: "Mitochondria (1-px core spine strokes + dim tubule tips)"
    label_2_channel_1: "Background (1-px inter-mitochondrial clefts + halo edges + empty space)"
  export_settings:
    source: "Probabilities"
    dtype: "float32"
    format: "tif"
```
