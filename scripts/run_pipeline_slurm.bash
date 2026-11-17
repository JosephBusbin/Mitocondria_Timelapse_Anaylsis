#!/bin/bash
# ==============================================================================
# run_pipeline_slurm.sh
# SLURM Job Array Orchestrator: 4D Mitochondrial Tracking & MitoLandscape Profiling
# Target System: William & Mary SciClone HPC Cluster (default: bora cluster)
# ==============================================================================
# HPC Resource Consciousness & Small-Job Throttling:
#   - Designed to run as many small, independent single-movie tasks (4 CPUs,
#     8 GB RAM, 45 min walltime) rather than a single monolithic multi-node job.
#     Small per-task footprints allow the SLURM backfill scheduler on 'bora' to
#     slot tasks into idle CPU cores immediately with near-zero queue wait time.
#   - Always include '%15' (concurrency throttle) in the --array directive so at
#     most 15 tasks execute simultaneously, preventing I/O contention on shared
#     '/sciclone/data10/' and '/sciclone/scr30/' filesystems:
#       module load slurm/bora
#       sbatch --array=1-100%15 run_pipeline_slurm.sh [pipeline_config.yaml] [movie_list.txt]
#   - For GPU-accelerated Cellpose 3 runs on a GPU sub-cluster, override at
#     submission time without editing this file:
#       sbatch --gres=gpu:1 --array=1-100%8 run_pipeline_slurm.sh
# ==============================================================================
#SBATCH --job-name=Mito4D_Task
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --array=1-1%15
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Enable strict bash error handling:
#   -e: exit immediately if any command exits with a non-zero status
#   -u: treat unset variables as an error
#   -o pipefail: return exit code of the rightmost failed command in a pipeline
set -euo pipefail

# Explicitly set working directory to the directory from which sbatch was invoked
cd "${SLURM_SUBMIT_DIR:-.}"

CURRENT_USER="${USER:-$(whoami)}"
N_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Bound C/Fortran numerical library threads strictly to the allocated SLURM CPUs
export OMP_NUM_THREADS="${N_THREADS}"
export OPENBLAS_NUM_THREADS="${N_THREADS}"
export MKL_NUM_THREADS="${N_THREADS}"
export NUMEXPR_NUM_THREADS="${N_THREADS}"
export VECLIB_MAXIMUM_THREADS="${N_THREADS}"

# ------------------------------------------------------------------------------
# 1. Load W&M Lmod Environment Modules & Activate Conda Environment
# ------------------------------------------------------------------------------
if command -v module >/dev/null 2>&1; then
    module load python/3.12.7 2>/dev/null || true
    module load miniforge3/24.9.2-0 2>/dev/null || true
    module load uv/0.13.1 2>/dev/null || true
fi

PROJECT_ROOT="${PROJECT_DATA_ROOT:-/sciclone/data10/${CURRENT_USER}/bestmanlab}"
TOOLS_DIR="${PROJECT_TOOLS_DIR:-/sciclone/data10/${CURRENT_USER}/mito_pipeline_tools}"
ENV_PREFIX="${TOOLS_DIR}/envs/mito4d_env"

if command -v conda >/dev/null 2>&1 && [[ -d "${ENV_PREFIX}" ]]; then
    eval "$(conda shell.bash hook)"
    conda activate "${ENV_PREFIX}"
elif [[ -f "${ENV_PREFIX}/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "${ENV_PREFIX}/bin/activate"
fi

# Configure node-local cache and scratch directories (/local/$USER)
LOCAL_NODE_SCRATCH="/local/${CURRENT_USER}/slurm_job_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-1}"
mkdir -p "${LOCAL_NODE_SCRATCH}" 2>/dev/null || LOCAL_NODE_SCRATCH="/tmp/${CURRENT_USER}_job_$$"
mkdir -p "${LOCAL_NODE_SCRATCH}"
export TMPDIR="${LOCAL_NODE_SCRATCH}"
export UV_CACHE_DIR="${LOCAL_NODE_SCRATCH}/uv_cache"
export PIP_CACHE_DIR="${LOCAL_NODE_SCRATCH}/pip_cache"

# ------------------------------------------------------------------------------
# 2. Parse Arguments & Resolve Target Movie from Array Index
# ------------------------------------------------------------------------------
CONFIG_PATH="${1:-pipeline_config.yaml}"
MOVIE_LIST="${2:-movie_list.txt}"
ARRAY_IDX="${SLURM_ARRAY_TASK_ID:-1}"

if [[ ! -f "${CONFIG_PATH}" ]]; then
    echo "[ERROR] Configuration file not found: ${CONFIG_PATH}" >&2
    exit 1
fi

if [[ ! -f "${MOVIE_LIST}" ]]; then
    echo "[ERROR] Movie manifest list not found: ${MOVIE_LIST}" >&2
    exit 1
fi

# sed -n "${ARRAY_IDX}p": extract only the line matching the current SLURM array task ID
MOVIE_PATH=$(sed -n "${ARRAY_IDX}p" "${MOVIE_LIST}")
if [[ -z "${MOVIE_PATH}" ]] || [[ ! -f "${MOVIE_PATH}" ]]; then
    echo "[ERROR] Invalid or missing movie path at line ${ARRAY_IDX} of ${MOVIE_LIST}: '${MOVIE_PATH}'" >&2
    exit 1
fi

MOVIE_NAME=$(basename "${MOVIE_PATH}")
MOVIE_ID="${MOVIE_NAME%.ome.tif}"
MOVIE_ID="${MOVIE_ID%.tif}"

# Derive experimental cohort/year from parent folder name
COHORT_DIR=$(basename "$(dirname "${MOVIE_PATH}")")
SCRATCH_ROOT="${SCRATCH:-/sciclone/scr30/${CURRENT_USER}}"
WORK_DIR="${SCRATCH_ROOT}/mito_work/${COHORT_DIR}/${MOVIE_ID}"
FINAL_OUT_DIR="${PROJECT_ROOT}/outputs/${COHORT_DIR}/${MOVIE_ID}"

mkdir -p "${WORK_DIR}"
mkdir -p "${FINAL_OUT_DIR}"

# Register cleanup trap on EXIT to remove node-local scratch and bulky intermediate arrays
cleanup_scratch() {
    rm -rf "${LOCAL_NODE_SCRATCH}"
    rm -f "${WORK_DIR}/${MOVIE_ID}_Probabilities.tif" \
          "${WORK_DIR}/${MOVIE_ID}_mito_norm.tif" \
          "${WORK_DIR}/${MOVIE_ID}_cell_RFP_reg.tif"
}
trap cleanup_scratch EXIT

echo "======================================================================"
echo "W&M SciClone 4D Mitochondrial Pipeline — Array Task ${ARRAY_IDX}"
echo "Compute Node : $(hostname) (${N_THREADS} CPUs allocated)"
echo "Movie ID     : ${MOVIE_ID} (Cohort: ${COHORT_DIR})"
echo "Source File  : ${MOVIE_PATH}"
echo "Scratch Work : ${WORK_DIR}"
echo "Final Output : ${FINAL_OUT_DIR}"
echo "Start Time   : $(date --iso-8601=seconds)"
echo "======================================================================"

# ------------------------------------------------------------------------------
# 3. Execute Pipeline Stages 3 through 5B
# ------------------------------------------------------------------------------
# Stage 3: Adaptive 4D Image Preprocessing (Fast4DReg, 3D Crop, Cell Mask, Bleach, Norm)
echo "[STEP 1/4] Executing Module 1: Stage 3 Adaptive 4D Preprocessing..."
python3 preprocess_4d_stack.py \
    --config "${CONFIG_PATH}" \
    --input "${MOVIE_PATH}" \
    --work-dir "${WORK_DIR}"

# Stage 4: Two-Stage Hybrid AI Segmentation (ilastik + Cellpose/Omnipose 2D+Z Stitch)
echo "[STEP 2/4] Executing Module 2: Stage 4 Two-Stage Hybrid Segmentation..."
python3 run_hybrid_segmentation.py \
    --config "${CONFIG_PATH}" \
    --movie-id "${MOVIE_ID}" \
    --work-dir "${WORK_DIR}"

# Stage 5A: TrackMate 7 Sparse LAP 4D Tracking (Fission & Fusion)
echo "[STEP 3/4] Executing Module 3: Stage 5A Sparse LAP Tracking..."
python3 run_trackmate_sparse_lap.py \
    --config "${CONFIG_PATH}" \
    --movie-id "${MOVIE_ID}" \
    --work-dir "${WORK_DIR}"

# Stage 5B: SNT & MitoLandscape 4D Subcellular Profiling
echo "[STEP 4/4] Executing Module 4: Stage 5B MitoLandscape 4D Profiling..."
python3 run_arbor_and_mitolandscape.py \
    --config "${CONFIG_PATH}" \
    --movie-id "${MOVIE_ID}" \
    --work-dir "${WORK_DIR}"

# ------------------------------------------------------------------------------
# 4. Copy Final Analytical Deliverables to Persistent Project Storage
# ------------------------------------------------------------------------------
echo "[STAGE 6] Archiving analytical CSVs, XMLs, SWCs, and provenance JSON..."
for artifact in \
    "${MOVIE_ID}_spots.csv" \
    "${MOVIE_ID}_edges.csv" \
    "${MOVIE_ID}_tracks.csv" \
    "${MOVIE_ID}_tracks.xml" \
    "${MOVIE_ID}_drift_shifts.csv" \
    "${MOVIE_ID}_snt_mitolandscape_tracks.csv" \
    "${MOVIE_ID}_process_0.5um_profile.csv" \
    "${MOVIE_ID}_remodeling_events.csv" \
    "${MOVIE_ID}_provenance.json" \
    "${MOVIE_ID}_arbor.swc" \
    "${MOVIE_ID}_cell_mask.tif"; do
    if [[ -f "${WORK_DIR}/${artifact}" ]]; then
        cp -p "${WORK_DIR}/${artifact}" "${FINAL_OUT_DIR}/"
    fi
done

echo "======================================================================"
echo "[SUCCESS] Completed 4D Pipeline for ${MOVIE_ID}"
echo "End Time     : $(date --iso-8601=seconds)"
echo "Deliverables : ${FINAL_OUT_DIR}"
echo "======================================================================"