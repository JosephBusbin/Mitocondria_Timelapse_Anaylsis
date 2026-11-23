# Mitocondria_Timelapse_Anaylsis
Image analysis pipeline for in vivo rapid 3D timelapse imaging of albino Xenopus laevis tadpoles neural progenitor cells. In collaboration with the Bestman Lab at The College of William and Mary.

---

## Running on SciClone (quick reference)

Run everything from `dev.sciclone.wm.edu`. Heavy work goes to compute nodes through Slurm. Load the cluster module in every new shell.

```bash
module load slurm/bora                      # makes sinfo/squeue/sbatch/srun target bora
source scripts/mito_env.bash                # storage paths; override with: export MITO_DATA_ROOT=...
# one-time install (≈15 min, ~4.2 GB under $MITO_TOOLS_DIR); refuses to run on dev itself
srun -p batch -N1 -c4 --mem=16G -t 01:00:00 bash scripts/setup_hpc_env.bash
# Stage 3 only (ilastik training data; keeps *_mito_norm.tif in $MITO_SCRATCH_ROOT)
sbatch --array=1-N%15 scripts/run_pipeline_slurm.bash pipeline_config.yaml "$MITO_RAW_ROOT/movie_list.txt" --stages 3
# full pipeline (needs $MITO_MODELS_DIR/bestman_mito_2d_v1.ilp)
sbatch --array=1-N%15 scripts/run_pipeline_slurm.bash pipeline_config.yaml "$MITO_RAW_ROOT/movie_list.txt"
squeue -u $USER; sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS
```

- Storage defaults to `/sciclone/data10/$USER/mito_pipeline`. To move to shared lab storage, change only `MITO_DATA_ROOT`.
- The Python environment works only on compute nodes, because its interpreter path exists there but not on dev. Use `srun` for interactive work: `srun --pty -c4 --mem=16G -t 1:00:00 bash`, then `source $MITO_TOOLS_DIR/activate_mito4d.sh`.

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
