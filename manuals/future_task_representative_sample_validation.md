# Future Task: Representative Randomized-Sample Validation of the Pipeline

**Status:** planned, not started (written 2026-11-24). **Start when:** the full lab archive is reachable from SciClone (Bestman lab storage) and the prerequisites below are done.

## Why this task exists
Every real-data decision so far was made on **one practice `.sld` file** (`05feb2025obf01avc01-B.sld`). That file is not representative:
- It is a **photoactivation** experiment. Plasmid OBF001 is mito-**paGFP**, which is dark until activated. Of its three two-channel 4D series:
  - #20 had no mito signal (0.9 σ);
  - #18 was borderline (2.6 σ, below the 3 σ QC gate);
  - only #39 had usable signal (35.8 σ).
- Its frame interval is 120 s, not the 11.7 s the lab believes is typical, and its pixel size is 0.19 µm, not 0.16 µm.
- It drifted up to 16–20 µm per movie.

The pipeline must therefore be re-validated on a **random, stratified sample of the real archive** before batch results are trusted. That includes the Stage 4B choice made on #39 alone (seeded 3D watershed as default, Cellpose as option).

## Prerequisites
1. Full archive scanned with the mito-construct filter: `collect_4d_sld_timelapses.py --mode scan-only --plasmid-csv … --mito-filter report`. The `mito_review.csv` decisions are recorded in an overrides CSV.
2. The calibration grouping decision is made (11.7 s vs ~120 s acquisitions; pixel size). The inventory's `time_interval_s` and `phys_size_x_um` columns are unit-correct.
3. The ilastik classifier is **retrained** on ≥ 3 diverse QC-PASS movies that are **not** in the validation sample (held out). The current `.ilp` was trained on noise and must not be used.
4. GPU time arranged for the Cellpose arm (user has access on other nodes). Record cluster, partition, GPU model and CUDA version.

## Sampling design
- **Frame:** all two-channel 4D series labelled `MITO` (not `MITO_PHOTOACTIVATABLE`) that pass the Stage 1 size criteria. Record the inventory file and its sha256.
- **Strata:**
  - cohort year;
  - acquisition protocol group (frame interval ≈ 11.7 s / ≈ 120 s / other);
  - mito colour (green / red, from `mito_color`).
  - Merge strata with fewer than 3 movies into their nearest neighbour and note it.
- **Size:** 24–30 movies in total, at least 3 per stratum, drawn **without replacement with a recorded random seed**. To be implemented then as `scripts/draw_validation_sample.py`: inventory CSV + seed → `validation_sample.csv` + a copy of the seed and inventory hash in provenance.
- **Annotation subset:** from 6–8 of the sampled movies (spread over strata), cut 10–15 crops (about 64 × 64 px × the Z-range of one structure, 1–2 frames each). Choose crops that **contain touching mitochondria** and crops that contain single isolated ones. A biologist hand-labels instances (Cellpose GUI or napari). These crops are the ground truth for Stage 4B and are never used for training.

## Per-movie checks (each stage, in order; use `scripts/qc_pipeline_outputs.py` + rendered overlays)
| Stage | What is checked | Pass/flag rule |
|---|---|---|
| 3 | `raw_mito_contrast.contrast_sigma`, drift curve, cell mask overlay (target cell only, no neighbour, no out-of-focus halo), registration (MAX-over-T sharp) | QC gate ≥ 3 σ; mask and registration judged by eye on a PNG |
| 4A | Foreground fraction in/out of mask; overlay of P > 0.5 on `mito_norm` | Foreground on mitochondria, not on haze/rim |
| 3→4 | **Late-frame noise amplification**: bleach correction matches every frame's histogram to frame 0, which stretches noise most in the most-bleached (last) frames. On practice #39 (built-in estimator) frame 5 had 145 objects vs 24–48, mostly soma specks. Compare objects/frame and median object volume vs frame index. | No trend of rising object count / falling median volume with frame index; if present, test a bleach model that does not amplify noise (e.g. exponential-fit ratio correction) |
| 4B | See the engine comparison below | — |
| 5A | Link displacement distribution vs `linking_max_distance_um` (per protocol group), track lengths, fission/fusion counts per minute, events by frame (no edge-frame pile-up) | Displacements not piled at the cap; event rates plausible and similar across frames |
| 5B | Fraction `outside_cell`, compartment fractions, geodesic distances ≥ 0, soma location on overlay | `outside_cell` < 5%; soma on the real soma |

## Stage 4B engine comparison: seeded 3D watershed vs Cellpose
**The distinction that matters: touching mitochondria.** ilastik only says which voxels are mitochondria. In a 1–2 µm radial-glia process, mitochondria lying end-to-end or side-by-side form **one continuous foreground region**. Stage 4B exists to split such regions into individual organelles. A method that merges touching mitochondria under-counts them, inflates volumes and creates false fusion/fission events; a method that splits single elongated mitochondria does the reverse. The comparison must measure **both errors**, not just overall coverage.

Arms (same ilastik probabilities, same min-volume and probability gates):
1. **Watershed** (current default): `engine: watershed`, default seed settings (0.4 µm XY / 1.0 µm Z separation, 0.05 µm prominence), plus one looser and one stricter seed setting.
2. **Cellpose, stock `cyto3`** on GPU: `engine: cellpose`, `diameter_px` 8 and 15.
3. **Cellpose, fine-tuned** on mitochondria from training-only movies (not the sample, not the annotation crops), on GPU. Record the training slices and epochs.

Metrics:
- **Touching-object accuracy** (primary), on annotated crops:
  - **merge errors**: one predicted object covering ≥ 2 true mitochondria;
  - **split errors**: one true mitochondrion covered by ≥ 2 predicted objects;
  - average precision at IoU 0.5 (AP50).
- **Foreground coverage**: % of P > 0.5 voxels assigned to an object.
- **Spill**: % of object voxels outside the foreground.
- **Morphology and resources**: objects per frame, volume and Z-extent distributions, runtime per frame (CPU vs GPU) and peak memory.
- **Downstream effect**: fission/fusion counts from Stage 5A under each arm.

Decision rule (user, 2026-11-24): keep Cellpose prominently as an option. **If Cellpose (stock and fine-tuned) fails the touching-object comparison on the real sample, it may be removed.** Cellpose "wins" a stratum only if it has lower merge + split error than the watershed at comparable coverage. Report per stratum, because the answer may differ between acquisition protocols.

**Hypothesis to test (why Cellpose is suspected to fail):** stock `cyto3` is trained on whole cells, which are about 30 px, textured and roundish. Mitochondria here are 2–6 px puncta and short rods, so the model operates far outside its training distribution. It must upsample every slice by about 30/diameter (slow on CPU), and on #39 it missed every punctum along the process (18% coverage, 45% spill). A fine-tuned model or Omnipose may close the gap; GPU access removes the speed objection but not the accuracy one. **Update this section and the README with the real result.**

## Stage 5A fission/fusion validation (added 2026-11-24 after the synthetic benchmark)
**Why:** on the realistic synthetic benchmark (`scripts/synthetic_benchmark.py`, 30 frames, 54 mitochondria, 3 true fissions + 3 true fusions):
- every true event was found, but **41 fissions and 44 fusions were predicted**;
- with the event persistence filter (`min_event_persistence_frames: 3`) there were still 23 fissions and 32 fusions.

The false events come from segmentation flicker among **crowded** mitochondria (< 2 µm apart: only ~61% resolved, versus ~98% of isolated ones). **Fission/fusion rates are therefore the least trustworthy pipeline output until validated on real data.**

**To do on the sample:**
- A biologist marks true events in 6–8 movies (frame, position), viewing raw `mito_norm` in Fiji; ideally blinded to the pipeline output.
- Score predicted events: precision and recall, matched within ±1 frame and 2 µm.
- Stratify by local crowding (nearest-neighbour distance) and by protocol (11.7 s vs 120 s: flicker persistence in frames depends on the frame interval).
- Re-tune `min_event_persistence_frames` and `split_merge_max_intensity_ratio` per protocol group.
- Consider reporting events only for isolated objects, or with a per-event confidence.
- **Until then, interpret event counts as upper bounds and compare conditions only within the same protocol group.**

## Outputs
- `validation_sample.csv` (seed, strata, movie list).
- A per-movie QC table.
- An overlay PNG set.
- The Stage 4B comparison table (per stratum).
- A short report with recommended defaults, linked from the README.
- An updated `pipeline_config.yaml` default if the evidence changes it.

## Resources (estimate; refine after the first 3 movies)
- **CPU (bora batch):** Stage 3 is about 1.5 min and 4 GB per movie. Stage 4 with the watershed is seconds per frame. Stage 5 is minutes.
- **GPU (other cluster):** the Cellpose arms. Stock `cyto3` on CPU took 7–17 min per frame, so the Cellpose arms run on GPU only.
