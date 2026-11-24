# Full-Archive Run Playbook (lab PC → Globus → SciClone → results)

This guide takes the whole lab archive (about 200–1000 movies, up to ~6 TB of raw `.sld`) through the pipeline on W&M SciClone. Every SciClone command is run from **`dev.sciclone.wm.edu`** (never a login node such as `bora`), from the repository folder.

**Every submit command is a dry run until you add `--execute`.** A dry run prints the exact `sbatch` commands, with comments that explain the resources, throttles and dependencies, and submits nothing. Read the dry run, then repeat the same command with `--execute`.

The tools behind this guide:

| Tool | What it does |
|---|---|
| `scripts/submit_full_run.bash` | One subcommand per phase. Plans, prints and (with `--execute`) submits. |
| `scripts/plan_full_run.py` | Size classes, the memory/time model, throttles, run manifest, retries, status. |
| `scripts/collect_4d_sld_timelapses.py --shard K/N` | Each array task scans or extracts its share of the `.sld` files. |
| `scripts/merge_collect_shards.py` | Joins the per-task reports into one inventory and `movie_list.txt`. |
| `scripts/run_pipeline_slurm.bash` | Stages 3–5 for one movie per array task. A stage it has already finished is skipped (see "Resume"). |
| `scripts/aggregate_results.py` | Cohort tables, QC dashboard and failure list. |

All throttles, size classes and model numbers live in the `hpc_run:` block of `pipeline_config.yaml`.

---

## 0. One-time setup for a run

```bash
cd ~/Mitocondria_Timelapse_Anaylsis
export MITO_RUN_ID=archive_2026_12        # any name: letters, digits, . _ -
source scripts/mito_env.bash              # MITO_RAW_ROOT, MITO_SCRATCH_ROOT, ...
```

Where things go:

| What | Where |
|---|---|
| Run state: manifest, shard reports, movie lists, results | `/sciclone/scr30/$USER/mito_runs/$MITO_RUN_ID/` |
| Dry-run state (never mixed with real runs) | `/sciclone/scr30/$USER/mito_runs/$MITO_RUN_ID/dryrun/` |
| Slurm logs (`<job name>_<job id>_<task>.out/.err`) | `/sciclone/scr30/$USER/mito_logs/$MITO_RUN_ID/` |
| Extracted OME-TIFFs | `$MITO_RAW_ROOT` (data10) |
| Per-movie deliverables | `$MITO_OUTPUT_ROOT/<cohort>/<movie>/` (data10) |
| Per-movie intermediates and `.done` markers | `$MITO_SCRATCH_ROOT/<cohort>/<movie>/` (scr30) |
| Final cohort tables | `<run>/results/`, also copied to `$MITO_OUTPUT_ROOT/_runs/$MITO_RUN_ID/` |

---

## 1. Get the data onto SciClone (lab PC side)

Order: **name pre-filter → Globus straight from the backup drive → only MITO + REVIEW files.**

### 1a. Pre-filter by file name on the lab PC (minutes)
The mito-construct filter reads only file names (no Bio-Formats), so it runs in seconds even over thousands of files. Use the PowerShell environment from `manuals/windows_lab_extraction_guide.md`:

```powershell
python scripts\construct_classifier.py `
    --plasmid-csv "D:\Plasmid_Primer_List_Bestman_Lab.csv" `
    --scan-dir "E:\SlideBook_backup" `
    --out "D:\mito_filter_results"
```

- Open `mito_classification.csv` and `mito_review.csv` (details in `manuals/construct_filter_guide.md`).
- **Transfer only files labelled `MITO` or `REVIEW`** (`REVIEW` = the name alone was not conclusive). This usually cuts the volume well below 6 TB.
- Write your REVIEW decisions into an overrides CSV (`pattern,label,note`) so the SciClone scan applies them.

### 1b. Globus directly from the backup drive (no copy to the lab PC first)
1. Install **Globus Connect Personal** on the lab PC. It can usually be installed per user, without admin rights. If it cannot, ask W&M IT, or use a machine where you have rights.
2. In Globus Connect Personal → *Options → Access*, add the **external backup drive** (for example `E:\`) as an accessible folder. The lab PC then only relays the data; there is no drive → PC copy step.
3. In the Globus web app, transfer from the lab-PC collection to the W&M SciClone collection, into e.g. `/sciclone/data10/$USER/bestman_lab_data/`. Keep the original folder structure: cohort years are read from folder and file names. Tick "verify file integrity".
4. You can transfer just the MITO + REVIEW list: use the Globus batch option, or copy those files into a staging folder on the backup drive first.

### Time estimates (planning numbers; the canary phase measures the real ones)

| Step | Assumption | Time |
|---|---|---|
| Name pre-filter on the lab PC | names only | minutes |
| Backup drive → Globus (via the lab PC) | USB 3 HDD ~100–150 MB/s; lab PC 1 GbE ~90–110 MB/s, so the network is the limit | **~15–20 h** for 6 TB (USB 2 drive: ~2.5 days). Less after pre-filtering |
| Globus → SciClone write | data10 NFS, parallel streams | not the bottleneck |
| Phase A scan (`showinf` per `.sld`, sharded array) | ~1–3 min per file, 20 at once | **< 1 h** for ~500 files |
| Phase B extraction (`bfconvert`, I/O-bound) | ~9 MB/s per stream × 20 streams ≈ 180 MB/s | **~3 h per TB of 4D data** (one lab PC: ~30 h/TB) |
| Phase D Stage 3, all movies | ~2–5 min each on 1 core, up to 180 at once | **~15–30 min** for 600 movies |
| Human gate: ilastik training (3 QC-PASS movies) | starts after the canary, not after all of Stage 3 | ~1–2 h hands-on |
| Phase F Stages 4–5, all movies | ilastik dominates (~10–15 min per 30-frame movie at 4 threads); 45 tasks × 4 cores = 180 cores | **~2–4 h** for 600 movies |
| Phase G aggregate + QC report | one job | minutes |

Transfer, not computing, is the slow part.

---

## 2. Phases on SciClone

Run each command once without `--execute` (dry run), read it, then run it again with `--execute`.

### Phase A: scan (+ merge) [+ extraction]
```bash
bash scripts/submit_full_run.bash scan \
    --source-dir /sciclone/data10/$USER/bestman_lab_data \
    --plasmid-csv /sciclone/data10/$USER/bestman_lab_data/Plasmid_Primer_List_Bestman_Lab.csv \
    --mito-overrides /sciclone/data10/$USER/bestman_lab_data/mito_overrides.csv \
    --execute
```
- Leave out `--mito-overrides` if you have no overrides file.
- This submits a **scan array** (about 5 `.sld` files per task, 1 core, 4 GB, at most **20 at once**) and a **merge job** (`afterany`: it runs even if a task failed, then reports the missing shards and exits with an error).
- Add `--chain-extract` to also queue the extraction now (`afterok` on the merge: it starts only after a clean merge). Without it, you check the inventory first (recommended for the first run).

**Check before extraction:**
- `<run>/scan/merged/sld_health_audit.csv`: how many files are unreadable, and why (status column).
- `<run>/scan/merged/sld_4d_inventory.csv`: the `time_interval_s`, `phys_size_x_um` and `mito_label` columns. Are the acquisition protocols the ones you expect?
- `<run>/scan/merged/mito_review.csv`: are there still REVIEW files? Then add overrides and re-run `scan`. Only `MITO` is extracted by default (`hpc_run.collector.extract_labels`); `--labels MITO,REVIEW` widens it.
- `<run>/extract/extract_files.txt`: the `.sld` files that will be extracted.
- If the merge reported missing shards, run `retry --phase scan`, then `merge` (Section 3).

### Phase B: extraction
```bash
bash scripts/submit_full_run.bash extract --execute
```
- One task per qualifying `.sld` file: 2 cores, 6 GB, time sized from the largest file, at most **20 at once** (I/O throttle).
- Files that already exist are skipped and SHA-256 checksums are written. The checksums are merged into `<run>/extract/merged/extracted_checksums.sha256`.
- **Check:** `cd $MITO_RAW_ROOT && sha256sum -c <run>/extract/merged/extracted_checksums.sha256` (run it inside a short `srun` job if the list is long), and that `<run>/scan/merged/movie_list.txt` points to files that now exist.

### Phase C: canary (~10 movies spanning the size classes)
```bash
bash scripts/submit_full_run.bash canary --execute
```
- Stage 3 only, on about 10 movies picked across the size classes **S** (< 2 GB uncompressed), **M** (2–8 GB) and **L** (≥ 8 GB, run on hima CPU nodes, never on GPU nodes). The largest movie of each class is always included.
- It measures real memory and time. The `stage3` phase reads the measurements from `sacct` and rescales the model.
- **Check:**
  - `python3 scripts/qc_pipeline_outputs.py` on the canary movies. Which QC status did each get: PASS, QC_WARNING, or FAIL_NO_SIGNAL (exit 3, expected for dark constructs)?
  - `bash scripts/submit_full_run.bash status --execute`.

### Human gate (between the canary and Stages 4–5)
1. **Retrain ilastik** on 3 diverse **QC-PASS** canary movies (`$MITO_SCRATCH_ROOT/<cohort>/<movie>/<movie>_mito_norm.tif`). Follow `manuals/ilastik_training_protocol.md` and copy the `.ilp` to `$MITO_MODELS_DIR`. The Stage 4 hash includes the `.ilp` contents, so a new model re-runs Stage 4–5 for every movie automatically.
   - **Before Phase F, run the pre-flight check** (`scripts/check_ilastik_project.py`, manual Step 7) with `--cell-mask` from one training movie. Only continue on **READY**.
   - Then **trial the model on one canary movie**: `bash scripts/run_pipeline_slurm.bash` via `srun` with `--stages 4`. Check that the log line `[STAGE 4A] Foreground (P > 0.5) covers X% of the cell mask` shows a minority (typically a few % to ~25%), not ~100%.
   - If Stage 4 exits with **code 5** (`FAIL_IMPLAUSIBLE_FOREGROUND`), the model calls most of the cell mitochondria. Retrain; do not override.
2. **Validation sample**: plan the randomized sample as described in `manuals/future_task_representative_sample_validation.md`. Movies used for training must not be in the sample. The Cellpose comparison runs on the hima GPUs.

### Phase D: Stage 3 for every movie
```bash
bash scripts/submit_full_run.bash stage3 --execute
```
- This prints the model refit from the canary (`memory xA, time xB`).
- It creates **one array per size class**, 1 core each. Memory = the largest prediction in the class × 1.5, time = × 2, with a minimum of 10 min.
- The **20% cap** (180 of bora `batch`'s 900 cores, read live with `sinfo`) is shared between the arrays. Lists longer than 1000 movies are split into chunks.
- Canary movies finish in seconds: their `.done` markers make Stage 3 a no-op for them.

### Phase F: Stages 4–5 (after the human gate)
```bash
bash scripts/submit_full_run.bash stage45 --execute
```
- One array per Stage 3 array: same movie list, so the same task number means the same movie. Each task uses 4 cores (ilastik threads + parallel watershed frames).
- `aftercorr` starts task *i* as soon as Stage 3 task *i* has succeeded, so Stage 4 runs while Stage 3 is still finishing others.
- `--kill-on-invalid-dep=yes` cancels the tasks whose Stage 3 failed (for example QC exit 3) instead of leaving them pending forever.
- The command refuses to run until a valid `.ilp` exists (`--allow-fallback` = test run with the built-in estimator).
- While Stage 3 arrays are still queued, Stage 4–5 gets only the part of the 20% cap that is left. The dry run prints `scontrol update ... ArrayTaskThrottle=N` lines; use them to raise the throttle once Stage 3 is finished.

### Phase G: aggregate
```bash
bash scripts/submit_full_run.bash aggregate --execute
```
- Runs `afterany` on every Stage 4–5 array. It writes the following to `<run>/results/` and copies them to `$MITO_OUTPUT_ROOT/_runs/$MITO_RUN_ID/`:
  - `all_spots.csv`, `all_spots_profiled.csv`, `all_tracks.csv`, `all_edges.csv`, `all_remodeling_events.csv`, `all_process_profiles.csv`. Each row is tagged with `movie_id`, `cohort`, `protocol` (= `dt…_dz…_dxy…`) and `mito_label`.
  - `qc_dashboard.csv`, one row per movie:
    - Stage 3 QC status, raw contrast (σ), drift and SNR;
    - objects per frame, tracks, fission and fusion counts;
    - `.done` times;
    - `sacct` runtime and peak memory per phase.
  - `failures.csv` (every unfinished movie with its reason) and `summary.json`.

---

## 3. Resume, retry, status

- **Status:** `bash scripts/submit_full_run.bash status`. This is a single `sacct` call that counts tasks by state for each array (`QC_FAIL(exit3)` is shown separately). Do not poll it in a loop: the dependencies sequence the phases, and checking every 5 minutes or less often is plenty.
- **Retry:** `bash scripts/submit_full_run.bash retry --phase stage3 --execute` (or `scan`, `extract`, `canary`, `stage45`).
  - It resubmits only tasks that ended `FAILED`, `OUT_OF_MEMORY`, `TIMEOUT` or `NODE_FAIL`, keeping their original task numbers.
  - Each attempt gets memory and time × 1.5, up to 2 retries per array; the manifest keeps count.
  - QC failures (exit 3) are never retried: they are a result, not an error.
  - Implausible-foreground failures (exit 5, Stage 4) are never retried either: they mean the ilastik model must be retrained. `status` shows them as `FG_IMPLAUSIBLE(exit5)`, and `failures.csv` names them `STAGE4_IMPLAUSIBLE_FOREGROUND`.
  - For `stage45`, it also resubmits tasks that were cancelled before starting because their Stage 3 failed, once that Stage 3 has succeeded in a retry.
  - After a scan retry, re-merge with `bash scripts/submit_full_run.bash merge --execute`.
- **Resume / re-run safely:** each finished stage writes `<movie>_stage<N>.done`.
  - The marker contains a hash of that stage's config section, the stage scripts, the git commit, the input movie (and, for Stage 4, the `.ilp`).
  - Re-submitting anything skips the stages that are already done with the same settings. Changing a Stage 4 setting re-runs Stages 4–5 only. Editing a comment re-runs nothing.
  - `--force` (passed to `run_pipeline_slurm.bash`) ignores the markers. Use it, for example, with `--stages 3` to regenerate a training stack after its intermediates were cleaned up.
- **Logs:** `/sciclone/scr30/$USER/mito_logs/$MITO_RUN_ID/<job name>_<array job id>_<task>.out` and `.err`.

---

## 4. Etiquette and rules (please keep)

- **20% footprint cap.** No phase holds more than 20% of a partition's cores at once. The cap is computed live from `sinfo`, with config fallbacks. The scan and extraction phases have an I/O throttle of 20 tasks at once.
- **bora first.** The submitter always runs `module load slurm/bora`. L-class movies run on **hima CPU nodes**; the GPU nodes (hi04, hi05, hi07) are excluded automatically.
- **Overflow is opt-in and only when bora is saturated:** `bash scripts/submit_full_run.bash overflow` (report only) / `... overflow --overflow --execute`.
  - Candidates come only from `hpc_run.overflow.clusters` (kuro `serial`, astral `batch`).
  - Each candidate partition's `PreemptMode` is checked live (`scontrol show partition`), and any partition that is not `OFF` is refused (gulf is `CANCEL`).
  - Clusters on `deny_clusters` are refused: the VIMS cluster (chesapeake) and gulf.
  - A per-cluster environment smoke test (`overflow --overflow --smoke --execute`, one short `srun`) must pass before any job is moved there. If it fails, build a per-cluster environment with `bash scripts/setup_hpc_env.bash --env-name mito4d_env_<cluster>` on that cluster.
  - Only whole chunks with no started task are moved. Kuro is the Zahn lab's priority cluster, so ask before using it even though it is on the list.
- **No `/tmp`, no `/local`.** Temporary files go to a per-task folder under `$MITO_SCRATCH_ROOT/tmp/` that is deleted when the task ends (bora `/local` sits on a 10 GB root disk; hima `/tmp` is a shared 10 GB partition).
- **dev is for planning and submitting only.** Never run the pipeline itself on dev or on a login node.
- **Purge.** scr30 files untouched for ~30 days are likely purged. Raw OME-TIFFs and final outputs live on data10 (later the Bestman lab area), and `results/` is copied there automatically. Everything left on scr30 can be re-created. Decide long-term retention with the lab.
