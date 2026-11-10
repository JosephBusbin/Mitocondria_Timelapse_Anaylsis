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

### Agent Modifications
- **2026-11-10 10:13–10:20**:
  - Linked workspace rules to global `AGENTS.md`.
  - Initialized `PROMPT_LOG.md` for local agent memory and multi-agent context continuity.
  - Initialized change tracking in `README.md` to document both agent actions and user modifications.
- **2026-11-10 10:35**:
  - Cleaned up redundant alias symlink (`Mitocondria_image_anaylsis`), keeping the canonical `Mitocondria_Timelapse_Anaylsis` folder.
  - Registered project path in `~/.gemini/projects.json`.
  - Updated `PROMPT_LOG.md` and `README.md` with latest user additions and architecture files.
