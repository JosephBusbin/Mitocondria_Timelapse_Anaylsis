# Mitocondria_Timelapse_Anaylsis
Image analysis pipeline for in vivo rapid 3D timelapse imaging of albino Xenopus laevis tadpoles neural progenitor cells. In collaboration with the Bestman Lab at The College of William and Mary.

---

## Pipeline at a glance

| Stage | Where | What it does | Script |
|---|---|---|---|
| 0. Construct filter | lab PC / SciClone | Classifies files by name (`mls`, lab plasmid IDs, common mito markers) → `MITO` / `MITO_PHOTOACTIVATABLE` / `NOT_MITO` / `REVIEW`, with mito colour | `scripts/construct_classifier.py` |
| 1–2. Scan + extract | SciClone (sharded arrays) | Health-audits `.sld` files, finds 4D series, extracts LZW OME-TIFF with checksums | `scripts/collect_4d_sld_timelapses.py` |
| 3. Preprocess | SciClone, 1 core | Drift registration (background-subtracted phase correlation), cell mask (glow-subtracted Triangle), **raw-signal QC (exit 3 = no mito signal)**, haze removal then gating, depth gain from the cell channel, bleach matching, normalisation | `scripts/preprocess_4d_stack.py` |
| 4A. Foreground | SciClone, 4 cores | ilastik Pixel Classification, headless → probability map; **exit 5 if it covers > 50% of the cell** (bad classifier) | `scripts/run_hybrid_segmentation.py` |
| 4B. Instances | same job | Splits **touching mitochondria**: seeded 3D watershed (default) or Cellpose (`engine: cellpose`, GPU on hima) | same |
| 5A. Tracking | same job | LAP linking, gap closing, fission/fusion with signal conservation + event persistence filter | `scripts/run_trackmate_sparse_lap.py` |
| 5B. Profiling | same job | Soma/process/endfoot compartments, internal geodesic distance, morphology, 0.5 µm profiles | `scripts/run_arbor_and_mitolandscape.py` |
| Aggregate | 1 job | Cohort tables, QC dashboard, failure list | `scripts/aggregate_results.py` |

## Running on SciClone

Run everything from `dev.sciclone.wm.edu` (never a login node). Heavy work goes through Slurm; run `module load slurm/bora` in every shell. The Python environment only works on compute nodes.

```bash
module load slurm/bora
source scripts/mito_env.bash                          # storage paths (change MITO_DATA_ROOT for lab storage)
srun -p batch -N1 -c4 --mem=16G -t 01:00:00 bash scripts/setup_hpc_env.bash      # one-time install
srun -p hima --gres=gpu:1 -c4 --mem=32G -t 01:00:00 bash scripts/setup_gpu_env.bash   # optional: GPU Cellpose env

# Whole archive: dry-run every phase first, then add --execute (manual: manuals/full_run_playbook.md)
bash scripts/submit_full_run.bash scan   --run-id RUN --source-dir <raw .sld dir> --plasmid-csv <plasmid table> --chain-extract
bash scripts/submit_full_run.bash canary --run-id RUN
bash scripts/submit_full_run.bash stage3 --run-id RUN
#   human gate: retrain ilastik on canary QC-PASS movies; check it with scripts/check_ilastik_project.py
bash scripts/submit_full_run.bash stage45   --run-id RUN
bash scripts/submit_full_run.bash aggregate --run-id RUN
bash scripts/submit_full_run.bash status    --run-id RUN --execute
```

Exit codes worth knowing:
- **3**: Stage 3 found no mito signal (a result, not retried).
- **5**: the Stage 4 foreground is implausible, so retrain ilastik (not retried).
- **4**: extraction failure in shard mode.

Every stage writes a `.done` marker keyed to its config, code and `.ilp`, so re-submitting only redoes what changed.

## Manuals and tools
- `manuals/full_run_playbook.md`: data transfer (Globus), phases, checks, retry, etiquette (≤ 20% of bora).
- `manuals/windows_lab_extraction_guide.md`: lab PC setup and extraction options (including the construct filter).
- `manuals/construct_filter_guide.md`: the filename filter, review CSV and overrides.
- `manuals/ilastik_training_protocol.md`: training rules (in-cell Background!) and the pre-flight check (Step 7).
- `manuals/future_task_representative_sample_validation.md`: the planned validation on a random sample of real data, including watershed vs Cellpose and fission/fusion validation.
- QC and testing tools:
  - `scripts/qc_pipeline_outputs.py`: per-stage QC report.
  - `scripts/synthetic_benchmark.py`: ground-truth benchmark.
  - `scripts/check_ilastik_project.py`: `.ilp` pre-flight check.
  - `scripts/test_pipeline_suite.py`: 43 tests; run on a compute node.
- Known limits (2026-11-24):
  - Fission/fusion counts are over-called in crowded regions (synthetic: about 15 predicted vs 3 true per type), so treat them as upper bounds until validated.
  - The only real test movie with signal so far is practice Series #39.

---

## Project Change & Activity Log

### User Modifications
- **2026-11-10 09:50**: Initial repository creation and project description setup.
- **2026-11-10 10:26–10:27**: 
  - Added Master Pipeline & Data Lifecycle Plan (`mitochondria_hpc_pipeline_plan.md`) covering 4D mitochondrial tracking, hybrid ilastik/Cellpose segmentation, and MitoLandscape integration.
  - Added Agentic Coding & Parameter Specification (`pipeline_agent_spec.md`).
  - Added default configuration file (`pipeline_config.yaml`).
  - Added pipeline scripts under `scripts/`:
    - `preprocess_4d_stack.py`: 4D registration and cytosolic gating.
    - `run_hybrid_segmentation.py`: ilastik + TrackMate-Cellpose segmentation.
    - `run_trackmate_sparse_lap.py`: 3D LAP tracking across frames.
    - `run_arbor_and_mitolandscape.py`: SNT auto-trace and 4D MitoLandscape compartment profiling.
    - `run_pipeline_slurm.bash`: Slurm batch execution script for the SciClone cluster.
- **2026-11-23 (commit `1b9bcfd`)**: Recovered and committed the previously missing `scripts/tiff_io.py`, plus two divergent drafts: `collect_4d_sld_timelapses_partial_update.py` and `setup_lab_windows_partial_update.ps1` (repo root).

### Agent Modifications
- **2026-11-10 10:13–10:20**:
  - Linked workspace rules to global `AGENTS.md`.
  - Initialized `PROMPT_LOG.md` for local agent memory and multi-agent context continuity.
  - Initialized change tracking in `README.md` to document both agent actions and user modifications.
- **2026-11-10 10:35**:
  - Cleaned up redundant alias symlink (`Mitocondria_image_anaylsis`), keeping the canonical `Mitocondria_Timelapse_Anaylsis` folder.
  - Registered project path in `~/.gemini/projects.json`.
  - Updated `PROMPT_LOG.md` and `README.md` with latest user additions and architecture files.
- **2026-11-23 (Claude Opus 5.5, branch `fix/reconcile-io-scripts`)**: Reconciled the recovered I/O scripts (plan Phase A).
  - Rationale: both `_partial_update` drafts were forked from older bases. Each added features, but each also reverted fixes already in `scripts/`. `tiff_io.py` had been missing entirely, so Stage 3–5 tests were silently skipped.
  - `scripts/collect_4d_sld_timelapses.py`:
    - Kept the SlideBook6Reader.jar support, the µm-safe decoding and the OME namespace fix.
    - Ported the draft's per-file health audit.
    - `sld_health_audit.csv` (one row per candidate) replaces `unparseable_sld_files.csv`.
    - Failures are now classified from the raw Java stderr. A new `MISSING_SB6_READER` status covers the crash from Bio-Formats' built-in SlideBook reader.
    - `movie_list.txt` is relative by default (the absolute prefix comes from `--hpc-raw-root` / `$MITO_RAW_ROOT`).
    - Jar and bftools discovery now also checks `<tools>/readers` and `$MITO_TOOLS_DIR`.
  - `scripts/setup_lab_windows.ps1`:
    - Kept the D:-drive guard, the hardened activation and the SlideBook6Reader.jar install.
    - Added `-BfToolsVersion` (default 8.5.0, in-place upgrade; downloads before replacing) and `-SbReaderUrl`.
    - The jar moved to `readers\` so upgrades cannot delete it.
    - Verification now checks the exact bftools version and the jar's size.
  - `scripts/tiff_io.py`:
    - Fixed OME-XML parsing (the old regex returned `{}` on bfconvert output, so calibration silently fell back to config defaults).
    - Added unit conversion to µm and seconds.
    - Fixed bogus 1.0 µm calibration read from uncalibrated ImageJ files.
    - Forced `photometric=minisblack` (Z=3/4 had been stored as RGB).
    - `.ome.tif` targets now write real OME-TIFF; 5D `.tif` writes an ImageJ TZCYX hyperstack.
    - New `read_tiff(path, axes="TCZYX")` normalization.
    - Missing calibration is reported in `meta["_warnings"]`.
  - `scripts/test_pipeline_suite.py`:
    - Dependency gates check third-party packages only, so a broken pipeline module now errors instead of skipping.
    - Fixed the config path (repo root).
    - Added 8 tests (01i–01k, 02b–02e). `MITO_PRACTICE_SLD` overrides the practice-file location.
  - `manuals/windows_lab_extraction_guide.md`: bftools 8.5.0, the `readers\` jar, relative movie lists, and the new health-audit status table.
  - Validation (dev server, light):
    - `py_compile` passes.
    - Suite under `module load python/3.12.7` (stdlib only): 22 tests, OK, 14 skipped.
    - Suite under a throwaway venv `/local/jbusbin/venvs/tiffio_check` (numpy 2.5.5, tifffile 2026.11.16, pyyaml): OK, 9 skipped (practice-file and Stage 3–5 tests), DeprecationWarnings raised as errors.
    - PowerShell is not available on SciClone, so `setup_lab_windows.ps1` still needs `-DryRun` and a full run on the lab PC.
- **2026-11-23 (Claude Opus 5.5, branch `fix/reconcile-io-scripts`, after user commit `85db0f5`)**: Phase B, HPC hardening and the first real-data run.
  - User actions:
    - Deleted the two `_partial_update` files.
    - Committed Phase A as `85db0f5`.
    - Approved Slurm test jobs on bora `batch` under the default account (no kuro).
  - Storage:
    - New `scripts/mito_env.bash` (bash) and `scripts/pipeline_config_io.py` (Python) define one settable root, `MITO_DATA_ROOT`, plus derived `MITO_RAW_ROOT`, `MITO_MODELS_DIR`, `MITO_OUTPUT_ROOT`, `MITO_TOOLS_DIR`, `MITO_SCRATCH_ROOT` and `MITO_ILASTIK_BIN`.
    - `pipeline_config.yaml` and `pipeline_agent_spec.md` use `${MITO_*}` placeholders. All four stage scripts load config through the shared loader, which errors on unexpanded variables.
    - Removed the unused `fiji_bin`.
  - `scripts/setup_hpc_env.bash` (rewritten):
    - Installs under `$MITO_TOOLS_DIR`: uv venv (module Python 3.12.7), CPU-only torch, `cellpose>=3.1,<4`, bftools 8.5.0, SlideBook6Reader.jar in `readers/`, ilastik 1.4.2, and MitoLandscape pinned at `45f3230`.
    - Writes `install_manifest.txt` and `env_freeze.txt`.
    - Refuses to run on dev.
    - Verifies the SlideBook native library loads.
  - `scripts/run_pipeline_slurm.bash` (rewritten):
    - Finds its repo under sbatch (`scontrol` Command=).
    - Resolves relative movie lists against `$MITO_RAW_ROOT` and strips CRLF.
    - Adds `--stages` and `--keep-intermediates`.
    - Fails fast when Stage 4 has no `.ilp`.
    - Deletes intermediates only after a complete 3–5 run.
  - Bugs found by running on the cluster:
    1. **3i's current Linux reader (SlideBook6Reader.jar-20260903101721) is broken upstream.** `libSBReadFile.so` reports `undefined symbol: _ZTI15CTileAnnotation`, so Bio-Formats silently used its built-in reader, which crashes. The HPC pins build `20250730183505` (verified to read all 49 series). Windows keeps the newer build (its DLL works).
    2. Cellpose 3 keyword is `do_3D`, not `do_3d`. `z_axis` and `channel_axis` are now explicit.
    3. Stage 5B fast marching crashed when the soma covered the whole cell (no zero contour). Added a zero-distance shortcut and a Dijkstra fallback.
    4. Activation now falls back to the node's own Java (dev has 1.8, compute nodes 11). Cellpose weights go to `$MITO_TOOLS_DIR/cellpose_models`, not `$HOME`.
    5. Fixed an invalid-escape `SyntaxWarning` in the arbor docstring.
    6. 3i's native library appends a multi-MB timing trace to a hard-coded relative file named `C:\Java_SBReadFileJNI_openFile.txt` in the current directory. On Linux this is a literal filename, and two copies landed in the repo during testing. The collector now runs Bio-Formats from a temporary directory (inputs and outputs resolved to absolute paths), and `.gitignore` excludes the name. Re-verified on bora: no new trace in the repo.
  - Behaviour changes:
    - The collector stages top-level `.sld` files under a `<year>/` folder, so each movie has a cohort directory.
    - Stage 3 reads with `axes="TCZYX"`, takes channel indices from the config, and records `calibration_warnings` in provenance.
    - Stages 4 and 5 read 4D stacks with `axes="TZYX"`.
  - Tests:
    - Added 01l (bash/Python defaults parity) and 01m (placeholder expansion).
    - On a bora compute node with `MITO_PRACTICE_SLD` set: all 24 pass (01–07 in the full run, 07–08 re-run after the 5B fix).
    - On dev: stdlib run OK (16 skipped); tifffile venv OK (9 skipped).
  - Real-data pilot (practice file `05feb2025obf01avc01-B.sld`):
    - Extracted on bora: Series #20 → `raw_ometiff/2025/…_t2r_-_4.ome.tif`, 0.456 GB, sha256 verified.
    - Stage 3 via `sbatch --stages 3` (job 703358): PASS, 83 s, MaxRSS 3.0 GB of 8 GB.
    - Cell crop 520×547×49. Max drift 0.58 µm. SNR 29.6.
    - Training stack: `$MITO_DATA_ROOT/ilastik_training/…_mito_norm_t0-4-9.tif`.
  - **Scientific flag:** the real file's calibration is XY 0.1877 µm and **dt 120 s**, against config defaults of 0.16 µm and 11.7 s. These were read from the file (no warnings). Pixel-unit parameters (`diameter_px`, `rolling_ball_radius_px`, `max_shift_*_px`) and per-frame tracking distances were tuned for the old assumptions and need review.
- **2026-11-23 (Claude Opus 5.5, cont.)**: Windows check, calibration investigation, cleanup.
  - User ran the suite on their own Windows PC (bftools 8.5.0 + the 2026-09 reader jar): 22 tests, OK, 12 skipped. The live-`.sld` tests 01e/01g/01h passed, so the Windows DLL works.
  - Manual and `setup_lab_windows.ps1`:
    - Step 1 no longer contains a user-specific path.
    - Step 4 explains that the leading `. ` (dot-source) is required.
    - The setup's final hint prints the portable `$env:LOCALAPPDATA` form.
  - Deleted the two stray `C:\Java_SBReadFileJNI_openFile.txt` trace files (user approved).
  - **Calibration verified from raw timestamps.**
    - Practice Series #20: per-plane DeltaT runs every 120 000 ms ± 2 ms, so the frame interval really is 120 s. XY is 0.1877 µm (OME has no TimeIncrement; tiff_io's DeltaT ms→s path is correct).
    - User-supplied max projections in `bestman_lab_data/`:
      - `MAXmlsEGFP_T7LFeb26_45TL_16bit.tif`: 30 frames, 0.1838 µm, finterval **11.71 s**.
      - `max_d318dec21jeb09jr19t1L_0-49.tif`: 50 frames, 0.1838 µm, no interval.
    - Both are single-channel TYX, so they can't enter the 4D pipeline.
    - Conclusion: 11.7 s is real for some acquisitions, but the practice movie uses a slower 120 s protocol. No source supports the 0.16 µm default.
  - **Fixed:** the collector's inventory ignored OME units, so `time_interval_s` held milliseconds (Series #20 showed 120000.444). The collector now converts `PhysicalSize*Unit` and `DeltaTUnit`/`TimeIncrementUnit` to µm and s (stdlib-only copy of the tiff_io tables). New test 01n. Inventory in `$MITO_RAW_ROOT` regenerated (S20 dt = 120.0 s).
- **2026-11-23 (Claude Opus 5.5, cont.)**: ilastik export/transfer fix after the user's first training attempt.
  - User trained a test classifier in ilastik 1.4.2 on the t0-4-9 stack. The training went "alright": the image is messy, with a mitochondrion merging into an amorphous mass.
    - GUI export to `tif` failed ("4D result into a 2D file format"), and so did `multipage tiff`.
    - Uploading the `.ilp` produced a **0-byte** `models/bestman_mito_2d_v1.ilp`.
  - **Root cause** (confirmed in ilastik 1.4.2's `opExportSlot.py`): `tif` export is 2D-only and `multipage tiff` is 3D-only. Probabilities for a `tzyx` stack are `tzyxc`.
  - **Pipeline bug found by the same cause:** `run_hybrid_segmentation.py` called headless ilastik with `--output_format=tif`. On real 4D data every Stage 4A run would have failed and quietly fallen back to the built-in estimator.
  - Fixes in `scripts/run_hybrid_segmentation.py`:
    - Headless export is now `--output_format=numpy --output_axis_order=tzyxc --export_dtype=float32` (no new dependency).
    - The raw `.npy` is shape-checked against `mito_norm`, reduced to the foreground class, written as `<movie>_Probabilities.tif`, and then deleted.
    - New `is_valid_ilastik_project()` (HDF5 signature check), so a 0-byte or non-HDF5 `.ilp` counts as missing.
    - If a valid `.ilp` exists but ilastik fails, Stage 4 now **errors** unless `MITO_ALLOW_ILASTIK_FALLBACK=1`.
    - Provenance records `stage4a_probability_source` (`ilastik` / `fallback_estimator`).
  - `scripts/run_pipeline_slurm.bash`: the Stage 4 fail-fast requires a non-empty HDF5 `.ilp`, not merely an existing file.
  - New test 05b (fake ilastik binary: export flags, `.npy` loading, shape/class checks, empty-`.ilp` rejection). Full suite on a bora compute node: 26 tests, OK (3 live-`.sld` tests skipped because `MITO_PRACTICE_SLD` was not set).
  - `manuals/ilastik_training_protocol.md`:
    - Step 2 now sets Storage to *Copy into project file*.
    - Step 4 adds "leave ambiguous regions unlabeled".
    - Step 5 is rewritten: no GUI export needed; why tif/multipage fail.
    - New Step 6: save, close, record SHA256, transfer with `scp` or VS Code Upload, verify with `sha256sum` and the HDF signature.
    - The provenance block is updated.
  - User decision recorded: the 11.7 s vs 120 s question is deferred until real data arrives. The working assumption stays 11.7 s, and heterogeneous data may be split into groups (an upcoming filtering task).
- **2026-11-23/24 (Claude Opus 5.5)**: First full Stage 3→5 run with the user's trained `.ilp`, plus per-stage QC.
  - Job 703461 ran Stages 3–5 successfully (1 h 56 min, 14.1 GB peak). Stage 4A used ilastik, as recorded in provenance. The ilastik `--ram` now follows the Slurm allocation.
  - New `scripts/qc_pipeline_outputs.py`: a read-only per-stage QC report.
  - QC finding: the practice movie's mito channel (ch1) has no detectable signal (in-cell minus background = −0.13 counts). Stage 4/5 outputs on it are segmentations of noise and are **not scientifically valid**.
  - Bugs found, fixes pending approval:
    - Stage 3 gates before the top-hat, which saturates thin processes and makes mask-edge artifacts.
    - The Stage 3 SNR QC passes noise.
    - Stage 5B mislabels outside-mask spots as basal with distance 0, and soma distances are negative.
  - Design weaknesses noted: Cellpose is not constrained by ilastik foreground, and the tracker calls every nearby appearance a fission.
- **2026-11-24 (Claude Opus 5.5)**: Pipeline fixes from stage-by-stage QC, plus a mito-construct filename filter.
  - **User modifications**:
    - Added `bestman_lab_data/Plasmid_Primer_List_Bestman_Lab.csv` (204 plasmids, read-only input).
    - Supplied two TIFF projections.
    - Trained and uploaded a test `.ilp`. It was trained on a series with no mito signal, so it must be retrained on a QC-PASS movie.
  - **Stage 3** (`preprocess_4d_stack.py`):
    - Raw-signal QC (`FAIL_NO_SIGNAL`, exit 3; the Slurm job stops before Stage 4 unless `MITO_IGNORE_QC=1`).
    - Top-hat now runs before gating.
    - Drift uses background-subtracted phase correlation with a per-frame step limit. Practice Series #39 drifted 79 px, which the old code reported as ≤ 5 px.
    - The cell mask subtracts out-of-focus glow before Triangle. The old mask merged a neighbouring soma in Series #20.
  - **Stage 4**: Cellpose runs only inside the ilastik foreground box, and objects are gated by mean probability. The fallback estimator uses a fixed scale.
  - **Stage 5A**: fission/fusion require signal conservation (`split_merge_max_intensity_ratio` 1.5), and CSVs are always written.
  - **Stage 5B**: outside-mask spots are snapped (≤ 1.5 µm) or labelled `outside_cell` instead of silently becoming `basal_process` at 0 µm. Soma distance ≥ 0.
  - New `scripts/qc_pipeline_outputs.py` (per-stage QC report).
  - **Mito-construct filter** (subagent-built, verified in the main session):
    - New files: `scripts/construct_classifier.py`, `config/construct_lexicon.json`, `scripts/test_construct_classifier.py` (37 tests) and `manuals/construct_filter_guide.md`.
    - File names are classified by the `mls` tag, plasmid IDs looked up in the lab table (`obf01` = OBF001, zero-padding tolerant) and common mito markers.
    - Labels: `MITO` / `MITO_PHOTOACTIVATABLE` (excluded by default) / `NOT_MITO` / `REVIEW`, plus `mito_color` and `red_mito_flag`.
    - Wired into `collect_4d_sld_timelapses.py` with `--plasmid-csv`, `--mito-filter {off,report,require}`, `--mito-overrides` and `--construct-lexicon`. Non-MITO files are skipped before Bio-Formats (`SKIPPED_NOT_MITO`).
    - New test 01o.
    - On the real folder, the practice `.sld` is classified `MITO_PHOTOACTIVATABLE` (mito-paGFP) and excluded.
- **2026-11-24 (Claude Opus 5.5, cont.)**: **Stage 4B default changed to a seeded 3D watershed; Cellpose kept as a prominent option.**
  - **Why Stage 4B exists — touching objects.** ilastik only labels which voxels are mitochondria. In a 1–2 µm radial-glia process, mitochondria lying end-to-end or side-by-side form ONE continuous foreground region. Stage 4B must split them:
    - merge errors under-count organelles, inflate volumes and fake fusion/fission;
    - split errors do the reverse.
    - Both engines are touching-object splitters.
  - **Evidence** (practice Series #39, frame 0, CPU, same probability map): watershed 100% coverage, 0% spill, 5 s; Cellpose `cyto3` d=8 7 objects, 18% coverage, 45% spill, 433 s, missing every punctum along the process. d=15 was no better.
  - **Why Cellpose is suspected to fail** (to be re-checked on real data): `cyto3` is trained on whole cells (about 30 px, textured), while mitochondria here are 2–6 px puncta and rods.
  - **Implementation** (`run_hybrid_segmentation.py`):
    - `segment_instances()` dispatches on `engine`.
    - `watershed_split_3d()` is now a real flooding watershed (`skimage`): h-maxima seeds on distance-to-edge × intensity, seed spacing and prominence set in µm (`watershed:` block in `pipeline_config.yaml`), and unseeded components are kept.
    - `engine: cellpose` uses a GPU automatically when torch sees one; if Cellpose cannot be imported, the watershed is used with a warning.
    - Provenance records `instance_engine`.
    - `min_volume_um3` raised from 0.05 to 0.10 (about 3 voxels) to drop 1–2-voxel noise specks.
    - The unused `anisotropy` key was removed.
  - Test 05d checks that a touching pair splits into 2, a 3-slice rod stays 1, two touching 3 × 3 puncta split into 2, and unknown engines are rejected. **33/33 pass on bora.**
  - With h-maxima seeds, #39 frame 0 gives 37 objects (previously 103 with the old nearest-seed partition). Endfoot overlay: `/sciclone/scr30/jbusbin/mito_qc_png/endfoot_ws_S39_t0.png`.
  - **Future task:** `manuals/future_task_representative_sample_validation.md`, a stratified random sample of the real archive. It includes a watershed vs Cellpose (stock and fine-tuned, GPU) comparison on hand-annotated touching mitochondria. The user approved removing Cellpose if it fails there.
  - Also updated `mitochondria_hpc_pipeline_plan.md` (Phase 4 status + Checkpoint 4), `pipeline_agent_spec.md` and `manuals/ilastik_training_protocol.md`.
- **2026-11-24 (Claude Opus 5.5, evening)**: full-run orchestration, dress rehearsal, ilastik v2 diagnosis, ground-truth benchmark, Stage 4/5 guards.
  - **User actions:** committed and pushed (c0981b7, 18e5d0a); confirmed the 9 filter rules (the PI will review `mito_review.csv`); trained and uploaded `bestman_mito_2d_v2.ilp` on S39 frames 0/2/5.
  - **Parallel full-run orchestration** (subagent-built, verified in the main session):
    - new `scripts/submit_full_run.bash` (phases `scan`/`merge`/`extract`/`canary`/`stage3`/`stage45`/`aggregate`/`retry`/`status`/`overflow`; dry-run by default, `--execute` to submit);
    - new `scripts/plan_full_run.py` (size classes, resource model with canary refit, throttles at 20% of the live partition, stage hashes);
    - new `scripts/merge_collect_shards.py` and `scripts/aggregate_results.py`; new `manuals/full_run_playbook.md`;
    - collector `--shard K/N`, `--file-list`, `--report-dir`;
    - `run_pipeline_slurm.bash`: per-stage `.done` markers, `TMPDIR` on scratch (bora `/local` is the 10 GB root fs), 1-CPU default;
    - watershed frames run in parallel.
    - Main-session change: stage hashes no longer include the git HEAD (any commit would have re-run every movie); the script sha256 still covers code changes.
  - **Dress rehearsal** on the practice `.sld` (run `rehearsal1`, `/sciclone/scr30/jbusbin/mito_rehearsal/`):
    - scan → merge → extract (3 series, 12 min) → canary (S18/S20 QC exit 3, S39 PASS) → stage3 (S39 skipped by its marker in 1 s; model refit memory ×1.47, time ×2.2) → stage45 (`aftercorr`; S18/S20 cancelled; S39 4 min 11 s, 12.3 GB) → aggregate;
    - outputs: `failures.csv` lists both QC failures with contrast; tables of 321 spots, 60 tracks, 16 events; protocol tag `dt120.0_dz1.00_dxy0.188`;
    - re-running stage45 was a no-op (2 s per task).
  - **ilastik v2 is not usable:** 15 Mitochondria / 16 Background label voxels, **0% of Background inside the cell**, 3D features. It calls **99.4% of the cell** foreground, so Stage 4 tiled the cytoplasm into ~95 fake objects per frame (built-in estimator: 2.1%). Fixes:
    - new **Stage 4 guard**: exit **5** `FAIL_IMPLAUSIBLE_FOREGROUND` when P>0.5 covers > 50% of the cell (`max_foreground_fraction_of_cell`); never retried; reported in status, `failures.csv` and the playbook;
    - new **`scripts/check_ilastik_project.py`** pre-flight checker (label counts, in-cell Background via `--cell-mask`, 2D features, storage, label order); manual Step 7 and the v2 lesson added.
    - **The user should train v3.**
  - **New `scripts/synthetic_benchmark.py`:** a ground-truth generator and evaluator for Stages 3–5 (radial-glia cell, drift, bleach, noise; stationary/motile/touching/fission/fusion mitochondria; registered truth). Results:
    - Small preset, 3 seeds: Stage 4 F1 0.88–0.96; **isolated objects 145/145**; crowded ~75%.
    - Realistic preset (30 × 49 × 664 × 952, 3.7 GB): **3.3 min per movie** (Stage 3 153 s, Stage 4 31 s, Stage 5 13 s), 11 GB peak; isolated 98%, crowded 61%, touching pairs 52/90; links R/P 0.90/0.85, displacement error 0.17 µm; compartments 0.86.
    - **Events: 3/3 fissions and 3/3 fusions found but 41 + 44 predicted.**
  - **Stage 5A:** new event persistence filter (`filter_transient_events`, `min_event_persistence_frames: 3`) removes transient and flicker split-merge pairs (realistic benchmark: 85 → 55 predicted events). **Crowded-region event calls remain unreliable**; validation was added to the future-task document.
  - Tests: 05e (foreground guard), 06c (flicker filter), test_06 extended to 6 frames, `TestFullRunOrchestration` 09a–09g. **42/42 pass on bora.**
- **2026-11-24 (Claude Opus 5.5, late)**: two Stage 3 bugs found by the synthetic benchmark (seeds 2/3 collapsed to F1 0.12–0.15) and fixed.
  1. `match_histogram_fast_16bit` interpolated over the reference CDF's plateau above its maximum, so each image's brightest voxels mapped to **65535**. This affects bleach matching in every movie. It now interpolates only over intensities present in the reference.
  2. **Depth (Z) attenuation correction rewritten.** The old version histogram-matched every slice to a bright one, which fired when slices were outside the cell or simply had no mitochondria, and stretched noise to mitochondria brightness (#39 briefly gave ~1,200 objects/frame). It now fits the **cytosolic channel's** in-cell signal against depth (exponential) and applies a smooth gain capped at `max_gain` 3.0 only if the fitted drop is > 30%. #39: drop 0.825 → gain 3.0; 20–37 objects/frame (70 at t=5), foreground 2.3% of cell.
  - After the fixes, the realistic synthetic seeds 1–3 give: F1 0.84 / 0.89 / 0.88; isolated objects 524/524, 580/581, 642/650; links R ≈ 0.96, P 0.89–0.93; compartments 0.97–0.99; ~3.4 min per 3.7 GB movie.
  - Test 04b (no saturation; Z correction ignores organelle content and follows the cell channel). **43/43 pass.**
  - README quick reference replaced by "Pipeline at a glance / Running / Manuals and tools".
- **2026-11-24 (Claude Opus 5.5, night)**: ilastik v3 and temporal seeding.
  - **User**: trained `bestman_mito_2d_v3.ilp` (2D features, 44% in-cell Background; 58/262 label voxels). On S39 it calls 0.4% of the cell foreground and gives 20–41 objects/frame, steady. Renders show the process chain + endfoot cluster with no soma texture. The saturated endfoot core is unlabelled, so v4 should add strokes inside it. Config → v3. The manual gained the user's Z-scroll observations (Mitochondria = distinct spots brighter than the texture; Background = the texture tracing the cell outline).
  - **Temporal seeding** (Stage 4B, default on: `watershed.temporal_seeding`), new `watershed_split_3d_temporal`:
    - frames are segmented in order; objects that were separate in the previous frame are kept apart unless no intensity dip separates them (90th percentile of the boundary ≥ 0.95 × the dimmer peak);
    - an object that was single needs a 2× more prominent peak to split.
    - Tuned on the realistic synthetic benchmark, 3 seeds, off → on:
      - F1 0.84–0.89 → 0.86–0.90;
      - crowded recall 65–74% → 74–82%;
      - touching pairs split 56–59 → 70–75 (/90);
      - link precision 0.89–0.93 → 0.93–0.95;
      - false-event load ~31 → 18–24 predicted (6 true);
      - true events still found (2–3/3 per type).
    - Stage 4 becomes sequential (~2–3 min instead of 0.5 min per realistic movie).
    - Tests 05f (keep apart / real fusion / no weak split) and 09f (pinned to independent mode). **44/44 pass.**
  - **Open**:
    - event persistence counted in frames is protocol-dependent (3 frames = 6 min at 120 s), so consider seconds;
    - consider an analysis frame window (the user prefers the cleaner early frames);
    - event calls still over-predict about 3–4× on synthetic data.
