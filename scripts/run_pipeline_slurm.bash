#!/bin/bash
# ==============================================================================
# run_pipeline_slurm.bash
# SLURM Job Array Orchestrator: 4D Mitochondrial Tracking & MitoLandscape Profiling
# Target System: William & Mary SciClone HPC Cluster (default: bora cluster)
# ==============================================================================
# Submitting (from dev.sciclone.wm.edu or any directory; the script finds its own repo):
#   module load slurm/bora
#   sbatch --array=1-100%15 scripts/run_pipeline_slurm.bash [config.yaml] [movie_list.txt] [options]
#
# Options:
#   --stages LIST          Comma list of stages to run (default: 3,4,5 or $MITO_STAGES).
#                          '--stages 3' = ilastik-training pilot: runs Stage 3 only and KEEPS
#                          <movie>_mito_norm.tif in the scratch work dir.
#   --keep-intermediates   Never delete the large scratch intermediates.
#
# Environment:
#   MITO_DATA_ROOT (and the other MITO_* paths) -> storage locations, see scripts/mito_env.bash
#   MITO_REPO_DIR                -> repo location if auto-detection fails
#   MITO_ALLOW_ILASTIK_FALLBACK=1 -> allow Stage 4 without a trained .ilp (uses the built-in
#                                   Difference-of-Gaussians foreground estimate; testing only)
#
# Inputs:
#   movie_list.txt lines are absolute paths or paths relative to $MITO_RAW_ROOT (the default
#   written by collect_4d_sld_timelapses.py). Windows CRLF line endings are tolerated.
#
# HPC Resource Consciousness & Small-Job Throttling:
#   - Many small independent single-movie tasks (4 CPUs, 8 GB RAM, 45 min) let the SLURM
#     backfill scheduler on 'bora' slot tasks into idle cores with near-zero queue wait.
#   - Always include '%15' (concurrency throttle) in --array so at most 15 tasks run at once,
#     preventing I/O contention on shared /sciclone/data10 and /sciclone/scr30.
# ==============================================================================
#SBATCH --job-name=Mito4D_Task
#SBATCH --partition=batch
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

# Relative paths given at submission are relative to the submission directory
cd "${SLURM_SUBMIT_DIR:-.}"

# ------------------------------------------------------------------------------
# 1. Locate the repository (sbatch runs a spooled COPY of this script, so BASH_SOURCE
#    points into slurmd's spool directory; ask Slurm for the original path instead)
# ------------------------------------------------------------------------------
if [[ -n "${MITO_REPO_DIR:-}" ]]; then
    SCRIPT_DIR="${MITO_REPO_DIR}/scripts"
else
    SCRIPT_PATH=""
    if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
        # grep -o prints only the match; \K drops 'Command=' from the printed match
        SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null | grep -oP '^\s*Command=\K\S+' | head -1 || true)"
    fi
    if [[ -z "${SCRIPT_PATH}" ]] || [[ ! -f "${SCRIPT_PATH}" ]]; then
        SCRIPT_PATH="${BASH_SOURCE[0]}"
    fi
    SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_PATH}")" && pwd)"
fi
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [[ ! -f "${SCRIPT_DIR}/preprocess_4d_stack.py" ]]; then
    echo "[ERROR] Pipeline scripts not found in ${SCRIPT_DIR}. Set MITO_REPO_DIR=/path/to/repo." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 2. Parse arguments
# ------------------------------------------------------------------------------
STAGES="${MITO_STAGES:-3,4,5}"
KEEP_INTERMEDIATES=0
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stages)             STAGES="$2"; shift 2 ;;
        --keep-intermediates) KEEP_INTERMEDIATES=1; shift ;;
        *)                    POSITIONAL+=("$1"); shift ;;
    esac
done
CONFIG_PATH="${POSITIONAL[0]:-${REPO_DIR}/pipeline_config.yaml}"
MOVIE_LIST="${POSITIONAL[1]:-movie_list.txt}"
ARRAY_IDX="${SLURM_ARRAY_TASK_ID:-1}"

if [[ ! "${STAGES}" =~ ^[345](,[345])*$ ]]; then
    echo "[ERROR] --stages must be a comma list of 3, 4, 5 (got '${STAGES}')" >&2
    exit 1
fi
has_stage() { [[ ",${STAGES}," == *",$1,"* ]]; }

[[ -f "${CONFIG_PATH}" ]] || { echo "[ERROR] Configuration file not found: ${CONFIG_PATH}" >&2; exit 1; }
[[ -f "${MOVIE_LIST}" ]]  || { echo "[ERROR] Movie manifest list not found: ${MOVIE_LIST}" >&2; exit 1; }
CONFIG_PATH="$(readlink -f "${CONFIG_PATH}")"

# ------------------------------------------------------------------------------
# 3. Environment: storage paths, tools, thread limits, node-local scratch
# ------------------------------------------------------------------------------
# shellcheck source=mito_env.bash
source "${SCRIPT_DIR}/mito_env.bash"
ACTIVATE_HELPER="${MITO_TOOLS_DIR}/activate_mito4d.sh"
if [[ ! -f "${ACTIVATE_HELPER}" ]]; then
    echo "[ERROR] ${ACTIVATE_HELPER} not found. Run scripts/setup_hpc_env.bash first." >&2
    exit 1
fi
set +u   # module/venv activation scripts reference unset variables
# shellcheck disable=SC1090
source "${ACTIVATE_HELPER}"
set -u

CURRENT_USER="${USER:-$(whoami)}"
N_THREADS="${SLURM_CPUS_PER_TASK:-4}"
# Bound C/Fortran numerical library threads strictly to the allocated SLURM CPUs
export OMP_NUM_THREADS="${N_THREADS}"
export OPENBLAS_NUM_THREADS="${N_THREADS}"
export MKL_NUM_THREADS="${N_THREADS}"
export NUMEXPR_NUM_THREADS="${N_THREADS}"
export VECLIB_MAXIMUM_THREADS="${N_THREADS}"

LOCAL_NODE_SCRATCH="/local/${CURRENT_USER}/slurm_job_${SLURM_JOB_ID:-$$}_${ARRAY_IDX}"
mkdir -p "${LOCAL_NODE_SCRATCH}" 2>/dev/null || LOCAL_NODE_SCRATCH="${TMPDIR:-/tmp}/${CURRENT_USER}_job_$$"
mkdir -p "${LOCAL_NODE_SCRATCH}"
export TMPDIR="${LOCAL_NODE_SCRATCH}"

# ------------------------------------------------------------------------------
# 4. Resolve the target movie from the array index
# ------------------------------------------------------------------------------
# sed -n "Np": print only line N; tr -d '\r': tolerate Windows CRLF manifests
MOVIE_ENTRY="$(sed -n "${ARRAY_IDX}p" "${MOVIE_LIST}" | tr -d '\r')"
if [[ -z "${MOVIE_ENTRY}" ]]; then
    echo "[ERROR] No movie at line ${ARRAY_IDX} of ${MOVIE_LIST}" >&2
    exit 1
fi
if [[ "${MOVIE_ENTRY}" == /* ]]; then
    MOVIE_PATH="${MOVIE_ENTRY}"
else
    MOVIE_PATH="${MITO_RAW_ROOT}/${MOVIE_ENTRY}"
fi
if [[ ! -f "${MOVIE_PATH}" ]]; then
    echo "[ERROR] Movie file not found: ${MOVIE_PATH} (line ${ARRAY_IDX}: '${MOVIE_ENTRY}')" >&2
    exit 1
fi

MOVIE_NAME="$(basename "${MOVIE_PATH}")"
MOVIE_ID="${MOVIE_NAME%.ome.tif}"
MOVIE_ID="${MOVIE_ID%.tif}"
COHORT_DIR="$(basename "$(dirname "${MOVIE_PATH}")")"
WORK_DIR="${MITO_SCRATCH_ROOT}/${COHORT_DIR}/${MOVIE_ID}"
FINAL_OUT_DIR="${MITO_OUTPUT_ROOT}/${COHORT_DIR}/${MOVIE_ID}"
mkdir -p "${WORK_DIR}" "${FINAL_OUT_DIR}"

# Fail fast: Stage 4 without a trained classifier would silently use the fallback estimator.
if has_stage 4; then
    ILP_PATH="$(python3 - "${SCRIPT_DIR}" "${CONFIG_PATH}" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from pipeline_config_io import load_config
cfg = load_config(sys.argv[2])
ilastik_cfg = cfg.get("stage4_hybrid_segmentation", {}).get("stage4a_ilastik", {})
print(f"{cfg['paths']['models_dir']}/{ilastik_cfg.get('project_file', 'bestman_mito_2d_v1.ilp')}")
PY
)"
    # An .ilp is HDF5: require its signature, so a 0-byte/failed upload is not mistaken for a model.
    # (head -c 4 | tail -c 3 keeps bytes 2-4 of the file, which read "HDF" in every HDF5 file.)
    ILP_OK=0
    if [[ -s "${ILP_PATH}" ]] && [[ "$(head -c 4 "${ILP_PATH}" | tail -c 3)" == "HDF" ]]; then
        ILP_OK=1
    fi
    if [[ "${ILP_OK}" -ne 1 ]] && [[ "${MITO_ALLOW_ILASTIK_FALLBACK:-0}" != "1" ]]; then
        echo "[ERROR] Stage 4 requested but no valid ilastik project (missing, empty, or not HDF5): ${ILP_PATH}" >&2
        echo "        Train it first (run with --stages 3 for training data), or set" >&2
        echo "        MITO_ALLOW_ILASTIK_FALLBACK=1 to test with the built-in estimator." >&2
        exit 1
    fi
fi

# Delete bulky intermediates only after a COMPLETE run through Stage 5 (never on failure,
# never for partial --stages runs, so Stage 3 outputs survive for ilastik training).
PIPELINE_COMPLETE=0
cleanup_scratch() {
    rm -rf "${LOCAL_NODE_SCRATCH}"
    if [[ "${PIPELINE_COMPLETE}" -eq 1 ]] && [[ "${KEEP_INTERMEDIATES}" -eq 0 ]]; then
        rm -f "${WORK_DIR}/${MOVIE_ID}_Probabilities.tif" \
              "${WORK_DIR}/${MOVIE_ID}_mito_norm.tif" \
              "${WORK_DIR}/${MOVIE_ID}_cell_RFP_reg.tif"
    fi
}
trap cleanup_scratch EXIT

echo "======================================================================"
echo "W&M SciClone 4D Mitochondrial Pipeline - Array Task ${ARRAY_IDX}"
echo "Compute Node : $(hostname) (${N_THREADS} CPUs allocated)"
echo "Repository   : ${REPO_DIR} ($(git -C "${REPO_DIR}" rev-parse --short HEAD 2>/dev/null || echo 'no git'))"
echo "Config       : ${CONFIG_PATH}"
echo "Stages       : ${STAGES}"
echo "Movie ID     : ${MOVIE_ID} (Cohort: ${COHORT_DIR})"
echo "Source File  : ${MOVIE_PATH}"
echo "Scratch Work : ${WORK_DIR}"
echo "Final Output : ${FINAL_OUT_DIR}"
echo "Start Time   : $(date --iso-8601=seconds)"
echo "======================================================================"

# ------------------------------------------------------------------------------
# 5. Execute requested pipeline stages
# ------------------------------------------------------------------------------
if has_stage 3; then
    echo "[STAGE 3] Adaptive 4D Preprocessing (drift, crop, cell mask, bleach, normalization)..."
    python3 "${SCRIPT_DIR}/preprocess_4d_stack.py" \
        --config "${CONFIG_PATH}" --input "${MOVIE_PATH}" --work-dir "${WORK_DIR}"
fi
if has_stage 4; then
    echo "[STAGE 4] Two-Stage Hybrid Segmentation (ilastik + Cellpose 2D/Z-stitch)..."
    python3 "${SCRIPT_DIR}/run_hybrid_segmentation.py" \
        --config "${CONFIG_PATH}" --movie-id "${MOVIE_ID}" --work-dir "${WORK_DIR}"
fi
if has_stage 5; then
    echo "[STAGE 5A] Sparse LAP Tracking (gap closing, fission, fusion)..."
    python3 "${SCRIPT_DIR}/run_trackmate_sparse_lap.py" \
        --config "${CONFIG_PATH}" --movie-id "${MOVIE_ID}" --work-dir "${WORK_DIR}"
    echo "[STAGE 5B] SNT & MitoLandscape 4D Subcellular Profiling..."
    python3 "${SCRIPT_DIR}/run_arbor_and_mitolandscape.py" \
        --config "${CONFIG_PATH}" --movie-id "${MOVIE_ID}" --work-dir "${WORK_DIR}"
fi

# ------------------------------------------------------------------------------
# 6. Copy deliverables to persistent project storage
# ------------------------------------------------------------------------------
echo "[ARCHIVE] Copying analytical CSVs, XMLs, SWCs, masks, and provenance JSON..."
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

if has_stage 3 && has_stage 4 && has_stage 5; then
    PIPELINE_COMPLETE=1
fi

echo "======================================================================"
echo "[SUCCESS] Completed stages ${STAGES} for ${MOVIE_ID}"
echo "End Time     : $(date --iso-8601=seconds)"
echo "Deliverables : ${FINAL_OUT_DIR}"
if ! has_stage 4; then
    echo "ilastik training input kept at: ${WORK_DIR}/${MOVIE_ID}_mito_norm.tif"
fi
echo "======================================================================"
