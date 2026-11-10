# Master Pipeline & Data Lifecycle Plan: 4D Mitochondrial Tracking & Arbor Tracing

## Overview & Guiding Principles

This plan outlines the end-to-end data flow from the Bestman Lab microscopy acquisition computer (Windows OS) to William & Mary's **SciClone HPC cluster** and back, paired with a companion machine-readable **Agentic Coding & Parameter Specification** (`pipeline_agent_spec.md`) so every script can be generated deterministically.

* **Zero Risk to Raw Lab Data**: The original `.sld` acquisition directories on the lab PC are treated as **strictly read-only**. All script outputs, extracted series, and staging folders are written to a separate drive or dedicated staging directory.
* **Two-Stage Hybrid AI Segmentation (Core Mainstay)**: Because radial glial progenitor cells have narrow, crowded cytoplasmic processes (~1–2 $\mu\text{m}$ diameter) where individual mitochondria stack end-to-end or pass across each other, simple pixel thresholding merges touching organelles into single blobs ([a critical bottleneck demonstrated in radial glia by Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)). This pipeline pairs **Stage 4A: [ilastik](https://www.ilastik.org/documentation/basics/startup) Semantic Denoising** (to strip out deep-tissue brain scattering) with **Stage 4B: [TrackMate-Cellpose](https://imagej.net/plugins/trackmate/detectors/trackmate-cellpose) / [Omnipose](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636021/) Gradient-Flow Instance Segmentation** as a core mainstay to split touching mitochondria inside thin radial tubes with minimal or zero manual intervention.
* **4D Extension of [MitoLandscape](https://github.com/enricoenne/MitoLandscape) Subcellular Compartmentalization**: Rather than relying solely on 1D centerline projections, Phase 5 integrates the open-source Python geodesic distance (`scikit-fmm`), Strahler branch ordering (`networkx`), and $0.5\ \mu\text{m}$ volumetric compartment-slicing algorithms from **[MitoLandscape (Negri, Fernández, & Borrell, 2025)](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)** alongside **[SNT Auto-trace](https://imagej.net/plugins/snt/auto-tracing)**, extending them from static 3D images into full **4D ($x, y, z, t$) developmental time-lapses**.
* **Reproducibility & Auditability**: Every stage generates a versioned manifest (CSV/JSON) or parameter log so any measurement can be traced back to its exact parent `.sld` file, series index, preprocessing settings, model weights, and tracking parameters.

## Literature Foundation & Precedents Driving This Architecture

| Reference Study & Codebase | Relevance to Bestman Lab Radial Glia Data | Specific Methods & Modules Integrated into Our Plan |
| :--- | :--- | :--- |
| **[MitoLandscape — Negri, Fernández, & Borrell (Nov 2025), *Frontiers in Cell & Dev. Bio.*](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)** ([GitHub](https://github.com/enricoenne/MitoLandscape)) | Built specifically for **apical Radial Glia Cells (RGCs)** electroporated in intact embryonic brains with cytoplasmic + mitochondrial reporters. | **1.** Validates **ilastik 2-class pixel classification** as superior to local thresholding in RGCs.<br>**2.** Provides open-source Python functions (`geodesic_dist` via `scikit-fmm`, `cell_annotation`, `strahler_analysis`, `process_analysis`) to map 3D geodesic distance from the soma, label apical/basal processes and endfeet, and slice twisted processes into $0.5\ \mu\text{m}$ volumetric bins. |
| **[Fast4DReg — Pylvänäinen, Laine, & Jacquemet et al. (Feb 2023), *Journal of Cell Science*](https://journals.biologists.com/jcs/article/136/4/jcs260728/287682/Fast4DReg-fast-registration-of-4D-microscopy)** & **[Review (Oct 2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7618379/)** | Fast, projection-based 3D+t drift registration and live-cell preprocessing benchmark. | Validates **Step 3A (`Fast4DReg`)**: estimates lateral $(X,Y)$ then axial $(Z)$ drift on the structural channel (`tagRFP-T`) 4×–70× faster than `Correct 3D drift` and applies the exact transformation table to the mitochondrial channel (`mito-EGFP`). |
| **[Baldwin, Eroglu et al. (Sept 2025), *Journal of Cell Biology*](https://pmc.ncbi.nlm.nih.gov/articles/PMC11527035/)** ([JCB Article](https://rupress.org/jcb/article/224/10/e202410130/278301/Mitochondrial-fission-controls-astrocyte)) | *In vivo* developing cortical radial glia/astrocytes co-electroporated with cytosolic `mCherry` and `mito-EGFP` to quantify Drp1-driven mitochondrial fission. | Validates **Step 3C (Cytosolic 3D ROI Gating)**: building a 3D cell territory mask from the cytosolic channel first and masking the mitochondrial channel inside it to eliminate non-target brain tissue signal. |
| **[Virga, Polleux, & Lewis et al. (Mar 2024), *Nature Communications*](https://pmc.ncbi.nlm.nih.gov/articles/PMC10923867/)** | *In vivo* electroporated neurons co-expressing cytoplasmic and mitochondrial markers to measure compartment-specific fission/fusion balance and segment occupancy. | Defines key neurobiological output metrics—**segment mitochondrial occupancy (%)**, **intra-segment length variability**, and **compartment-resolved fission/fusion rates**—which they had to trace manually in 2D, and which our Phase 5 automates in 3D/4D. |
| **[MitoSkel — Zaghbani et al. (Aug 2025)](https://www.sciencedirect.com/science/article/pii/S1746809425002733)**, **[MoDL — Ding et al. (Jan 2025)](https://pubmed.ncbi.nlm.nih.gov/39820041/)**, **[Cellpose 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903308/)**, & **[Omnipose](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636021/)** | Deep-learning instance and tubular segmentation of live-cell mitochondria. | Resolves the core binarization limitation noted in [MitoLandscape (Negri et al., 2025)](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full), where permissive ilastik thresholds merge adjacent mitochondria in thin tubes, by using spatial gradient flows (**Stage 4B**) to split touching organelles. |

## Stage-by-Stage Data Lifecycle Summary

| Phase | Location | Input File Type & Status | Action / Transformation | Output File Type & Status |
| :--- | :--- | :--- | :--- | :--- |
| **1. Safety & Backup** | Lab PC + External Drive | `.sld` + companion `_Files` folders *(Raw, untouched)* | Read-only dry run on copied test files; full cold backup of raw drives | Verified backup `.sld` archive + test `sld_4d_inventory.csv` |
| **2. Staging & Transfer** | Lab PC $\to$ SciClone HPC (`/sciclone/data10/` or `scr30/`) | `.sld` *(Read-only source)* | Filter for 4D ($Z \ge 15, T \ge 10$); extract & LZW-compress 4D series; transfer via Globus | `*.ome.tif` *(16-bit, 2-channel $C,Z,T$, lossless LZW)* + `sld_4d_inventory.csv` |
| **3. Adaptive 4D Preprocessing** | SciClone HPC (Headless FIJI / Python) | `*.ome.tif` *(Dual-channel: Ch1 `mito-EGFP`, Ch2 `tagRFP-T`)* | `Fast4DReg` 3D drift correction on Ch2 applied to Ch1; Ch2 3D cell mask & cytosolic ROI gating; Ch1 bleach correction, rolling-ball + median denoising, & cross-day percentile normalization | `*_mito_norm.tif` *(16-bit standardized)* + `*_cell_RFP_reg.tif` + `*_cell_mask.tif` |
| **4. Two-Stage Hybrid Segmentation** | Lab PC (one-time setup) $\to$ SciClone HPC (batch GPU/CPU) | `*_mito_norm.tif` *(16-bit standardized)* | **Stage 4A**: Headless ilastik removes tissue haze (`*_Probabilities.tif`); **Stage 4B**: Cellpose 3 / Omnipose 3D splits touching mitochondria in thin tubes | `*_Probabilities.tif` *(32-bit float)* $\to$ `*_mito_labels.tif` *(16-bit 3D instance masks)* |
| **5. 4D Tracking & MitoLandscape Arbor Profiling** | SciClone HPC (SLURM job arrays by year/cohort) | `*_mito_labels.tif` + `*_cell_RFP_reg.tif` + `*_cell_mask.tif` | **5A**: TrackMate 7 (`LabelImageDetector` + `SparseLAPTracker`); **5B**: SNT Auto-Trace + `MitoLandscape` geodesic distance (`scikit-fmm`), Strahler branch order, & $0.5\ \mu\text{m}$ compartment slicing | `tracks.xml`, `spots.csv`, `edges.csv`, `tracks.csv` + `.swc`, `*_process_0.5um_profile.csv`, & `*_snt_mitolandscape_tracks.csv` |
| **6. Return & Archival** | SciClone HPC $\to$ Lab PC | Lightweight analytical outputs & summary movies | Transfer results back via Globus; retain parameter logs and scripts in Git | Local analysis-ready CSVs, XMLs, SWCs, and provenance logs |

## Phase 1: Safety, Local Verification & Cold Backup

1. **Isolated Local Script Testing**
   * **Action**: Copy 2–3 representative `.sld` files (plus any adjacent SlideBook companion folders) into an isolated `C:\Test_Sandbox\` folder. Run `collect_4d_sld_timelapses.py --mode scan-only` and `--mode extract-ometiff` on the sandbox copies only.
   * **Verification**: Open the extracted `.ome.tif` in FIJI (`Bio-Formats Importer`) to confirm that both channels (`mito-EGFP` and `tagRFP-T`), all $Z$-slices, all $T$-timepoints, voxel dimensions ($\mu\text{m}$), and frame intervals ($\Delta t$) match the original `.sld` metadata.
2. **Full Read-Only Backup of Raw Lab Files**
   * **Action**: Before running batch operations across the full multi-terabyte dataset, mirror the raw `.sld` directories onto an external backup drive (or institutional NAS) and set the raw acquisition folders to **Read-Only** permissions in Windows so no script or user action can accidentally modify, move, or overwrite them.

## Phase 2: Environment Setup, Filtering & Globus Transfer

3. **Lab Computer Dependencies (Assuming Base Windows OS + FIJI Only)**
   Because the lab PC already has FIJI installed, it **already contains a bundled Java Runtime (OpenJDK)** inside `Fiji.app\java\`. You only need three lightweight, non-intrusive additions:
   * **Python 3 (3.10+)**: Standard Windows installer from `python.org` (uses only built-in standard libraries—`csv`, `xml`, `subprocess`, `pathlib`, `concurrent.futures`—so **zero `pip` packages are required** for the filtering script).
   * **Bio-Formats Command Line Tools (`bftools`)**: Download `bftools.zip` from [Open Microscopy Environment](https://www.openmicroscopy.org/bio-formats/downloads/) and unzip it to `C:\bftools\`. Point its `BF_JAVA_HOME` environment variable to FIJI's existing Java folder.
   * **Globus Connect Personal (Windows)**: Lightweight background client to link the lab PC staging folder to W&M's Globus system.
4. **HPC Group & Storage Provisioning at William & Mary**
   * **Action**: Coordinate with W&M Research Computing (`hpc-help@wm.edu`) to establish a shared UNIX group (e.g., `bestmanlab`) and shared directory on `/sciclone/data10/` (persistent project storage) and `/sciclone/scr30/` (high-speed Lustre parallel scratch for active SLURM jobs).
   * **Directory Structure**: Organize by year/experiment cohort (`raw_ometiff/2023/`, `raw_ometiff/2024/`, `models/`, `outputs/`, `logs/`).
5. **Globus Endpoint Configuration**
   * **Action**: Connect **[Globus Connect Personal](https://www.wm.edu/offices/it/services/researchcomputing/using/filesandfilesystems/xfers/globus/personal/)** on the lab PC to the **`W&M RC/HPC/IO Collection`** (`io.sciclone.wm.edu`), targeting the lab's shared `/sciclone/data10/` directory with automatic checksum verification enabled.
6. **Batch Filtering, Lossless Compression & Transfer**
   * **Action**: Run `collect_4d_sld_timelapses.py --mode extract-ometiff --compression LZW` pointing from the read-only raw `.sld` drive to a dedicated staging folder (or external SSD).
   * **Data Transition**:
     * *Start*: Multi-TB `.sld` containers (mixed 2D snapshots, static 3D stacks, and 4D movies).
     * *End*: Filtered, single-series **16-bit LZW-compressed `.ome.tif` hyperstacks** (containing only 4D time-lapses) + `sld_4d_inventory.csv` master manifest transferred to SciClone via Globus.

## Phase 3: Adaptive 4D Image Preprocessing on HPC

Because movies are acquired across different days, tadpole developmental stages, imaging depths (up to $150\ \mu\text{m}$), and laser conditions, every stack passes through an automated 6-step preprocessing script (`preprocess_4d_stack.py`) on the HPC before segmentation. Every step adapts mathematically to each movie's internal statistics:

7. **Step 3A — 3D Sample & Stage Drift Registration via `Fast4DReg` (Locked to Cytosolic Channel)**
   * **Problem**: Tadpole heartbeat, gill movement, or thermal stage drift shifts the entire brain volume over the 5.5-minute acquisition. TrackMate cannot separate stage drift from true organelle motility.
   * **Solution & Precedent**: Run FIJI's **[Fast4DReg](https://journals.biologists.com/jcs/article/136/4/jcs260728/287682/Fast4DReg-fast-registration-of-4D-microscopy)** ([Pylvänäinen et al., 2023, *Journal of Cell Science*](https://journals.biologists.com/jcs/article/136/4/jcs260728/287682/Fast4DReg-fast-registration-of-4D-microscopy)). First, run `Fast4DReg` (`Estimate drift`) strictly on **Channel 2 (`tagRFP-T` cell outline)** using 2D maximum-intensity $XY$ projections to solve lateral $(X,Y)$ shifts followed by $XZ/YZ$ projections to solve axial $(Z)$ shifts—an approach specifically suited for anisotropic $1.0\ \mu\text{m}$ $Z$-stacks and 4×–70× faster than `Correct 3D drift`. Second, run `Fast4DReg` (`Apply drift`) using the exported `<movie_id>_drift_table.csv` on **Channel 1 (`mito-EGFP`)**, stabilizing the cell in 3D space without canceling out internal mitochondrial transport.
8. **Step 3B — Channel Splitting & Automatic 3D Bounding-Box Cropping**
   * **Empty-Volume Cropping**: A radial glial cell typically occupies a narrow diagonal corridor (~15%–25% of the $512 \times 512 \times 45$ voxel field of view).
   * **Solution**: Split the registered stack into `Ch1_mito` and `Ch2_cell`. Compute a coarse 3D bounding box around the `tagRFP-T` signal (padded by $5\ \mu\text{m}$ in $x,y$ and $2\ \mu\text{m}$ in $z$) and crop both channels while logging the $(x_0, y_0, z_0)$ origin offset. This cuts downstream ilastik and Cellpose compute time and RAM usage by **60%–80%**.
9. **Step 3C — Cytosolic 3D ROI Gating & Automated Cell Mask Generation**
   * **Out-of-Cell Autofluorescence Masking**: In live *Xenopus* optic tectum, autofluorescent lipid droplets, skin pigment granules, or neighboring electroporated axons outside the target radial glial cell can mimic small mitochondria ([a critical tissue isolation step emphasized in both Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full) and [Baldwin et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11527035/)).
   * **Solution**: Generate a closed 3D cell mask (`*_cell_mask.tif`) from the time-averaged projection of registered Channel 2 (`tagRFP-T`) using 3D Gaussian smoothing ($\sigma_{xy}=2.0, \sigma_z=1.0$) + Triangle/ilastik thresholding + `MorphoLibJ` largest-connected-component filtering (replacing manual 3D Slicer seed painting in [MitoLandscape](https://github.com/enricoenne/MitoLandscape)). Multiplying Channel 1 by a dilated version ($r = 1.5\ \mu\text{m}$) of this mask zeroes out **100% of non-cellular background**, while the un-dilated `*_cell_mask.tif` is saved directly for Phase 5B `MitoLandscape` compartment and geodesic distance analysis.
10. **Step 3D — Depth Attenuation ($Z$-Drop) & Photobleaching Compensation**
    * **Axial ($Z$) & Temporal ($T$) Attenuation**: Scattering across a $44\ \mu\text{m}$-thick brain stack causes deeper $z$-slices to appear systematically dimmer than superficial slices ([eLife 2026](https://elifesciences.org/reviewed-preprints/107154v2)), while 100% laser exposure causes temporal photobleaching from $t=0$ to $t=30$ ([Feng & Bestman, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8357028/)).
    * **Solution & Precedent**: Apply FIJI's **[Bleach Correction](https://pmc.ncbi.nlm.nih.gov/articles/PMC7871415/)** ([Miura, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7871415/)) using the **Histogram Matching** algorithm across time ($T$), and across depth ($Z$) if axial drop-off exceeds 30% ([validated in *C. elegans* live 4D mitochondrial imaging, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11541191/)).
11. **Step 3E — Shot-Noise Despeckling & Local Haze Subtraction**
    * **Solution & Precedent**: Following **[Mitochondria Analyzer](https://pmc.ncbi.nlm.nih.gov/articles/PMC7052579/)** ([Chaudhry et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7052579/)), apply a 3D $XY$-Median filter (`radius = 1 px`) to remove single-pixel sCMOS/EM-CCD photon shot noise while keeping tubule edges sharp, followed by **Rolling Ball Background Subtraction** (`radius = 15 px`, $\approx 2.4\ \mu\text{m}$) to strip out diffuse cytosolic GFP haze inside the soma and shaft.
12. **Step 3F — Cross-Day Percentile Intensity Normalization**
    * **Solution**: Rescale each movie's intensity range between its 1.0th percentile ($P_{1.0} \to 0$) and 99.9th percentile ($P_{99.9} \to 50,000$) in 16-bit unsigned integer (`uint16`) space so every movie enters ilastik and Cellpose with an identical dynamic range.

## Phase 4: Two-Stage Hybrid Segmentation (ilastik + Cellpose/Omnipose Mainstay)

In Benchmark Figure 10 and Supplementary Figure S2 of **[MitoLandscape (Negri et al., 2025)](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)**, the authors proved that segmentation stringency is the single largest source of variance in radial glia mitochondrial analysis: high thresholds lose dim tubules, whereas lower thresholds cause false-positive merging of adjacent mitochondria inside thin processes. Our two-stage hybrid architecture directly solves this trade-off:

13. **Stage 4A — ilastik Semantic Foreground Probability Mapping**
    * **Local One-Time Setup**: Train a 2-class ilastik Pixel Classifier (`Class 1: Mitochondria`, `Class 2: Background/Haze`) on 3 normalized training stacks from Phase 3 using **2D multi-scale features** ($\sigma = 0.7\text{ to }5.0\text{ px}$ in $X,Y$ only, avoiding $1.0\ \mu\text{m}$ $Z$-blur artifacts).
    * **HPC Batch Execution**: Run `run_ilastik.sh --headless` to convert `*_mito_norm.tif` into clean 32-bit probability maps (`*_Probabilities.tif`, values $0.0–1.0$).
14. **Stage 4B — Deep-Learning Instance Segmentation in Thin Tubes ([Cellpose 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903308/) / [Omnipose](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636021/))**
    * **Why This Mainstay Is Necessary**: While ilastik cleanly isolates mitochondrial foreground from brain tissue noise, **[Cellpose 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903308/)** ([Stringer & Pachitariu, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903308/)) and **[Omnipose](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636021/)** ([Cutler et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636021/)) predict orthogonal spatial gradient flows to split touching, elongated tubular organelles that share a continuous probability cloud in thin radial shafts.
    * **Minimal-Work / Automated Configuration for Anisotropic 3D Data**:
      * *Zero-Training Option*: Run Cellpose 3's `cyto3` model (or Omnipose's `bact_fluor_omni` tubular model) on the stacked `[mito_norm, Probabilities]` 2-channel input with `diameter = 8.0` pixels ($\approx 1.3\ \mu\text{m}$).
      * *Minimal Human-in-the-Loop Option (30 mins total)*: Open 5 cropped $Z$-slices in the Cellpose GUI, click "Run", correct 10–15 touching mitochondria in the thin radial shaft, and click "Train" to save `bestman_mito_tube.pth`.
      * *Anisotropic 3D Stitching (`stitch_threshold = 0.35`, `anisotropy = 6.25`)*: Runs 2D flow prediction per slice and links 3D instances across $Z$ using IoU overlap (`stitch_threshold = 0.35`), outputting a 16-bit 3D instance label stack (`*_mito_labels.tif`) where every mitochondrion has a unique integer ID ($1 \dots N$).

## Phase 5: Automated 4D Tracking & MitoLandscape Subcellular Profiling on HPC

15. **Step 5A — Batch TrackMate 7 Execution (`LabelImageDetector` + `SparseLAPTracker`)**
    * **Action**: Run SLURM job arrays executing headless FIJI on `*_mito_labels.tif` using TrackMate's **`LabelImageDetector`** paired with the **`SparseLAPTracker`** (fission/splitting and fusion/merging enabled), with `*_mito_norm.tif` mounted as Channel 2 to record fluorescence intensity alongside 3D volume and shape.
    * **Coordinate & ID Preservation**: `spots.csv` exports both the TrackMate `SPOT_ID` and the underlying integer `LABEL_ID` (via intensity readout on the label channel or centroid voxel lookup in `*_mito_labels.tif`) so Module 5B links every tracked spot directly to its 3D voxel mask.
    * **Outputs**: `tracks.xml`, `spots.csv`, `edges.csv`, `tracks.csv`, and `run_config.json`.
16. **Step 5B — Integrated SNT Auto-Trace & [MitoLandscape](https://github.com/enricoenne/MitoLandscape) 4D Subcellular Profiling**
    * **Action**: Run `run_arbor_and_mitolandscape.py` on the HPC, combining **[SNT Auto-trace](https://imagej.net/plugins/snt/auto-tracing)** with adapted modules from the open-source **[MitoLandscape Python repository (`enricoenne/MitoLandscape`)](https://github.com/enricoenne/MitoLandscape)** ([Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)):
      1. **Automated Soma & Compartment Segmentation (`cell_annotation` & `apical_progenitor_compartment_annotation`)**: Uses the 3D cell mask (`*_cell_mask.tif`) and SNT soma seed (`EDT_Threshold`) to automatically partition the radial glial cell into **Soma (`1`)**, **Apical Process (`2`)**, **Basal Process (`3`)**, and **Endfoot / Distal Terminal Branches (`4`)**.
      2. **3D Geodesic Distance Map (`geodesic_dist` via `scikit-fmm`)**: Computes the true internal 3D geodesic distance (in $\mu\text{m}$) from the soma boundary through the twisted 3D cytoplasm of the radial glial process for every voxel inside the cell mask.
      3. **Strahler Branch Order & Skeleton Graph (`strahler_analysis` via `networkx`)**: Skeletonizes the 3D cell mask, removes false loops via `minimum_spanning_tree` weighted by soma geodesic distance, and assigns **Strahler branch order** ($1 = \text{terminal branchlets}$, $2 = \text{parent branch}$, $\dots$, $\text{main radial shaft}$) and unique branch IDs to every segment of the arbor.
      4. **Mitochondrial Graph Barycenter & Subtype Classification (`mito_analysis`)**: Skeletonizes each 3D mitochondrial instance from Stage 4B in 26-connectivity (`networkx`), locates its 3D **graph barycenter** (true topological center along curved tubules), measures tubule length, volume, surface area (`MorphoLibJ` Crofton method), branch points, and classifies each organelle into **`punctus`**, **`rod`**, or **`network`** ([Bakare et al., 2021; Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)).
      5. **$0.5\ \mu\text{m}$ Process Volumetric Slicing & Occupancy (`process_analysis`)**: Slices the radial glial process into **$0.5\ \mu\text{m}$ geodesic distance bins** from the soma, calculating local cytoplasm volume, local mitochondrial volume, **mitochondrial volume density / segment occupancy (%)** ([matching Virga et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10923867/)), and local fission/fusion event density per bin over time.

## Phase 6: Returning Results & Documentation

17. **Packaging & Returning Deliverables to the Lab Computer**
    * **Action**: Transfer the lightweight `outputs/` directory back to the lab PC via Globus:
      * Master inventory and provenance logs (`sld_4d_inventory.csv`, `<movie_id>_provenance.json`, `run_config.json`).
      * TrackMate 4D tables (`spots.csv`, `edges.csv`, `tracks.csv`, `tracks.xml`) and SNT arbor reconstructions (`.swc`).
      * **New Enriched MitoLandscape-4D Tables**:
        * `<movie_id>_snt_mitolandscape_tracks.csv`: Every tracked mitochondrion at every timepoint $t$ annotated with its 3D graph barycenter, geodesic soma distance ($d_{\text{geo}}$), compartment label (soma, apical, basal, endfoot), Strahler branch order, morphological class (`punctus`/`rod`/`network`), signed axial velocity ($v_{\text{geo}}$), and transport state (anterograde, retrograde, docked).
        * `<movie_id>_process_0.5um_profile.csv`: Longitudinal $0.5\ \mu\text{m}$ bin profile along the radial glial arbor reporting cytoplasm cross-sectional volume, mitochondrial volume occupancy (%), subtype ratios, and localized fission/fusion event counts.

## Future Expansion Options (Out of Current Baseline Scope)

* **Self-Supervised Isotropic $Z$-Restoration ([CARE 3D](https://pure.mpg.de/rest/items/item_3188430_2/component/file_3188432/content))**: If sub-fragment 3D skeleton tracking (such as [MitoGraph](https://pmc.ncbi.nlm.nih.gov/articles/PMC6322684/) + [MitoTNT](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011060)) is ever desired in future work, CARE's self-supervised isotropic restoration mode can be trained directly on the high-resolution $0.16\ \mu\text{m}$ $XY$ slices of the existing `.ome.tif` archive to upsample the $1.0\ \mu\text{m}$ $Z$-axis to isotropic $0.16\ \mu\text{m}$ voxels. It is excluded from the baseline pipeline because Stage 4B's 2D-to-3D IoU slice stitching (`stitch_threshold = 0.35`) paired with `MitoLandscape` already resolves 3D instances and geodesic topology natively without synthetic $Z$-interpolation.

## Reproducibility & Documentation Checkpoints

Use these checkpoints as a live lab notebook checklist as you complete each phase. Recording these details in a `README_provenance.md` file inside your project repository guarantees that another student or reviewer can reproduce every figure and table years later.

### Checkpoint 1: Pre-Flight Safety & Backup Verification (After Phase 1)
* [ ] **Raw Data Inventory**: Record the exact source drive letter/path, total number of `.sld` files, total byte size, and date range of experiments covered.
* [ ] **Backup Location & Verification**: Log the physical label/serial number or NAS path of the backup drive, the date of backup completion, and confirmation that Windows Read-Only permissions are active on the original folders.
* [ ] **Sandbox Validation Notes**: Record the exact filenames of the 2–3 test `.sld` files used in `C:\Test_Sandbox\` and verify that metadata ($\Delta x, \Delta y, \Delta z, \Delta t$) in FIJI matched SlideBook's native view.

### Checkpoint 2: Environment, Extraction & Globus Transfer Log (After Phase 2)
* [ ] **Lab PC Software Versions**:
  * Windows OS build/version.
  * Python version (`python --version`) and Git commit hash of `collect_4d_sld_timelapses.py`.
  * Bio-Formats `bftools` version (`showinf -version`) and bundled Java version (`java -version`).
  * Globus Connect Personal version.
* [ ] **4D Filtering Criteria**: Document the exact command-line flags used (`--min-z`, `--min-t`, `--min-c`, `--compression LZW`), plus any manual exclusions.
* [ ] **Transfer Provenance**: Save `sld_4d_inventory.csv`, the Globus Transfer Task IDs, and the destination path on SciClone (`/sciclone/data10/...`).

### Checkpoint 3: Adaptive 4D Preprocessing Log (After Phase 3)
* [ ] **3D Drift Registration (`Fast4DReg`)**: Record `Fast4DReg` version, reference channel (`Ch2 tagRFP-T`), projection type (`Max Intensity`), and inspect the exported drift table (`_drift_table.csv`) for any movie exceeding $15\ \mu\text{m}$ excursion.
* [ ] **3D Cropping & Cytosolic Masking**: Record the bounding-box padding ($\mu\text{m}$), the thresholding/MorphoLibJ settings used to generate `*_cell_mask.tif`, the dilation radius ($r = 1.5\ \mu\text{m}$) used for extracellular noise gating, and the $(x_0, y_0, z_0)$ crop offset saved in the movie's metadata JSON.
* [ ] **Bleach & Axial Attenuation Correction**: Confirm `Histogram Matching` mode was used across $T$ (and note whether $Z$-attenuation correction was triggered).
* [ ] **Denoising & Normalization Parameters**: Record the rolling-ball radius (`15 px`), median filter radius (`1 px`), and the exact percentile bounds (`P_1.0` and `P_99.9` mapped to `[0, 50000]`) logged per movie in `<movie_id>_provenance.json`.

### Checkpoint 4: Two-Stage Hybrid Segmentation Log (ilastik + Cellpose/Omnipose) (After Phase 4)
* [ ] **Stage 4A (ilastik) Provenance**:
  * Exact ilastik version on PC and HPC, filenames/timepoints used for training, selected 2D feature scales ($\sigma$), and brush-stroke decision rules.
* [ ] **Stage 4B (Cellpose / Omnipose) Provenance**:
  * Exact package versions (`cellpose` / `omnipose`, `torch`, CUDA version).
  * Model identifier (`cyto3`, `bact_fluor_omni`, or custom `bestman_mito_tube.pth` + list of the 5–10 training slices if fine-tuned).
  * Exact inference parameters: `diameter` (px), `flow_threshold` (default `0.4`), `cellprob_threshold` (default `0.0`), `anisotropy` ($\Delta z / \Delta xy = 6.25$), and `stitch_threshold` (`0.35` for 2D-to-3D slice stitching).

### Checkpoint 5: TrackMate 7, MitoLandscape Profiling & Archival Sign-Off (After Phases 5–6)
* [ ] **TrackMate & Sparse LAP Config**: Record headless FIJI and TrackMate versions, `LabelImageDetector` settings, `LINKING_MAX_DISTANCE` ($2.5\ \mu\text{m}$), `GAP_CLOSING_MAX_DISTANCE` ($3.0\ \mu\text{m}$, `MAX_FRAME_GAP = 2`), `SPLITTING_MAX_DISTANCE` ($2.0\ \mu\text{m}$), `MERGING_MAX_DISTANCE` ($2.0\ \mu\text{m}$), and minimum track length (`NUMBER_SPOTS >= 4`).
* [ ] **SNT & MitoLandscape Config**: Record SNT version, `scikit-fmm` and `networkx` versions, soma detection settings (`EDT_Threshold`), minor-process length cutoff (default $5.0\ \mu\text{m}$), geodesic slice bin thickness (default $0.5\ \mu\text{m}$), and minimum mitochondrial volume filter ($0.05\ \mu\text{m}^3$ per [Negri et al., 2025](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2025.1668779/full)).
* [ ] **Visual Spot-Check & Scratch Purge**: Confirm visual overlay checks on 2–3 returned `tracks.xml` files per cohort in desktop FIJI, tag the Git repository commit with all `run_config.json` files, and record the date `/sciclone/scr30/` intermediates were purged.