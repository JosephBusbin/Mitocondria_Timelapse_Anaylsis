# Mitocondria_Timelapse_Anaylsis
Image analysis pipeline for in vivo rapid 3D timelapse imaging of albino Xenopus laevis tadpoles neural progenitor cells. In collaboration with the Bestman Lab at The College of William and Mary.

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
