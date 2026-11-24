#!/bin/bash
# ==============================================================================
# setup_gpu_env.bash
# Optional GPU Python environment for Stage 4B 'engine: cellpose' (SciClone 'hima' GPUs)
# ==============================================================================
# Why a separate environment:
#   The main pipeline env (setup_hpc_env.bash -> envs/mito4d_env) uses CPU-only PyTorch: bora has
#   no GPUs and CUDA wheels add several GB. This script builds envs/mito4d_gpu_env with the SAME
#   package versions as the CPU env (read from $MITO_TOOLS_DIR/env_freeze.txt) except PyTorch,
#   which comes from a CUDA wheel index. The CPU pipeline is never modified.
#
# GPU architecture check:
#   hima has Tesla P100 (compute capability 6.0, nodes hi04/hi05) and V100 (7.0, hi07). Newer
#   PyTorch CUDA builds drop old architectures, so candidate CUDA indexes are tried newest-first
#   and the first build whose kernel list covers every --require-cc value AND that runs a real
#   kernel on this node's GPU is kept.
#
# Policy (AGENTS.md): user-space only (uv venv), no conda init, caches in /local/$USER, run INSIDE
# a GPU job, never on dev:
#   module load slurm/bora
#   srun -p hima --gres=gpu:1 -c 4 --mem=32G -t 01:00:00 bash scripts/setup_gpu_env.bash
#
# Usage:
#   bash setup_gpu_env.bash [--tools-dir DIR] [--env-name NAME] [--cuda-indexes "cu130 cu128 ..."]
#                           [--require-cc "6.0 7.0"] [--dry-run]
# Outputs (under $MITO_TOOLS_DIR):
#   envs/<env-name>/            GPU venv
#   activate_mito4d_gpu.sh      Source in GPU jobs (loads the normal helper, then this venv)
#   gpu_install_manifest.txt    torch/CUDA/arch list/GPU used; gpu_env_freeze.txt exact versions
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CURRENT_USER="${USER:-$(whoami)}"
ENV_NAME="mito4d_gpu_env"
# Candidates tried in order: "<index>" = newest torch on that index, "<index>@<ver>" = pinned torch version.
# Recent builds dropped Pascal/Volta (cu130 torch 2.14 starts at sm_75), hence the pinned fallbacks.
CUDA_INDEXES="cu128 cu126 cu124 cu126@2.7.1 cu118@2.7.1 cu130"
REQUIRE_CC="6.0 7.0"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tools-dir)     export MITO_TOOLS_DIR="$2"; shift 2 ;;
        --env-name)      ENV_NAME="$2"; shift 2 ;;
        --cuda-indexes)  CUDA_INDEXES="$2"; shift 2 ;;
        --require-cc)    REQUIRE_CC="$2"; shift 2 ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       sed -n '2,32p' "$0"; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 2 ;;
    esac
done

# shellcheck disable=SC1091
source "${SCRIPT_DIR}/mito_env.bash"
TOOLS_DIR="${MITO_TOOLS_DIR}"
ENV_PREFIX="${TOOLS_DIR}/envs/${ENV_NAME}"
CPU_FREEZE="${TOOLS_DIR}/env_freeze.txt"
GPU_HELPER="${TOOLS_DIR}/activate_mito4d_gpu.sh"

# --- Guards: must run on a GPU compute node, after the CPU env exists ----------------------------
if [[ "$(hostname -s)" == dev* ]]; then
    echo "[ERROR] Refusing to run on the dev server. Use: srun -p hima --gres=gpu:1 ... bash $0" >&2
    exit 1
fi
if [[ ! -f "${CPU_FREEZE}" ]]; then
    echo "[ERROR] ${CPU_FREEZE} not found: run setup_hpc_env.bash first (it defines the pinned versions)." >&2
    exit 1
fi
if ! command -v nvidia-smi >/dev/null 2>&1 || ! nvidia-smi -L >/dev/null 2>&1; then
    echo "[ERROR] No visible GPU. Request one: srun -p hima --gres=gpu:1 ..." >&2
    exit 1
fi
GPU_LINE="$(nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader | head -1)"

# Node-local /local/$USER when it exists; otherwise scratch. NEVER /tmp: on hima GPU nodes /tmp is a
# 10 GB shared partition and CUDA wheels filled it to 100% (2026-11-24).
LOCAL_BASE="/local/${CURRENT_USER}"
if [[ ! -d "${LOCAL_BASE}" ]] || [[ ! -w "${LOCAL_BASE}" ]]; then
    LOCAL_BASE="/sciclone/scr30/${CURRENT_USER}/node_cache"
fi
mkdir -p "${LOCAL_BASE}"
export UV_CACHE_DIR="${LOCAL_BASE}/uv_cache"
export PIP_CACHE_DIR="${LOCAL_BASE}/pip_cache"
WORK_TMP="$(mktemp -d "${LOCAL_BASE}/gpu_setup_XXXXXX")"
trap 'rm -rf "${WORK_TMP}"' EXIT

echo "======================================================================"
echo "GPU environment setup (Stage 4B engine: cellpose)"
echo "Host          : $(hostname)   GPU: ${GPU_LINE}"
echo "Env prefix    : ${ENV_PREFIX}"
echo "Pinned from   : ${CPU_FREEZE}"
echo "CUDA indexes  : ${CUDA_INDEXES}   (required compute capabilities: ${REQUIRE_CC})"
echo "======================================================================"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN] Would create ${ENV_PREFIX}, try torch from each CUDA index, install the pinned"
    echo "          non-torch packages, verify a GPU kernel, and write ${GPU_HELPER} + manifest."
    exit 0
fi

set +u
if ! command -v module >/dev/null 2>&1 && [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/Modules/init/bash
fi
module load python/3.12.7
module load uv/0.13.1
set -u
PYTHON_BIN="$(command -v python3)"

VENV_PY="${ENV_PREFIX}/bin/python"

# --- 1. CUDA PyTorch: first candidate index that covers the required GPUs and runs a kernel -------
# Prints "OK <torch> <cuda> <archs>" only if every required compute capability is in the kernel
# list and a real tensor op succeeds on this node's GPU ("no kernel image" errors fail here).
CHECK_PY='
import sys, torch
need = [c.replace(".", "") for c in sys.argv[1].split()]
archs = torch.cuda.get_arch_list() if torch.cuda.is_available() else []
missing = [c for c in need if f"sm_{c}" not in archs]
if not torch.cuda.is_available() or missing:
    print("FAIL", torch.__version__, "cuda_available=%s" % torch.cuda.is_available(), "missing_sm=%s" % missing, archs)
    sys.exit(1)
x = torch.arange(1_000_000, device="cuda", dtype=torch.float32)
assert abs(float((x * 2).sum()) - 999999000000.0) / 999999000000.0 < 1e-3
print("OK", torch.__version__, torch.version.cuda, " ".join(archs))
'
TORCH_INDEX=""
for cand in ${CUDA_INDEXES}; do
    idx="${cand%@*}"
    spec="torch"
    [[ "${cand}" == *@* ]] && spec="torch==${cand#*@}"
    echo "[INFO] Trying ${spec} from https://download.pytorch.org/whl/${idx} (fresh venv) ..."
    # Fresh venv per candidate: leftover CUDA libraries from a rejected build break resolution.
    uv venv --clear -q "${ENV_PREFIX}" --python "${PYTHON_BIN}"
    if ! uv pip install --python "${VENV_PY}" --index-url "https://download.pytorch.org/whl/${idx}" \
            "${spec}" >"${WORK_TMP}/uv_${idx}.log" 2>&1; then
        echo "       not installable: $(grep -m1 -i -E 'error|no solution|because' "${WORK_TMP}/uv_${idx}.log" || tail -1 "${WORK_TMP}/uv_${idx}.log")"
        continue
    fi
    if RESULT="$("${VENV_PY}" -c "${CHECK_PY}" "${REQUIRE_CC}" 2>&1)"; then
        echo "       ${RESULT}"
        TORCH_INDEX="${idx} (${spec})"
        break
    fi
    echo "       rejected: $(echo "${RESULT}" | tail -1)"
done
if [[ -z "${TORCH_INDEX}" ]]; then
    echo "[ERROR] No candidate CUDA PyTorch build supports compute capabilities '${REQUIRE_CC}' on this node." >&2
    echo "        Retry with --require-cc \"7.0\" (V100 only) or older --cuda-indexes." >&2
    exit 1
fi

# --- 2. Every other package at the CPU env's exact versions (torch excluded) -----------------------
grep -v -i -E '^(torch|nvidia-|triton)' "${CPU_FREEZE}" > "${WORK_TMP}/pinned_no_torch.txt"
uv pip install --python "${VENV_PY}" -r "${WORK_TMP}/pinned_no_torch.txt"
uv pip freeze --python "${VENV_PY}" > "${TOOLS_DIR}/gpu_env_freeze.txt"

# --- 3. End-to-end check: Cellpose sees the GPU and segments a synthetic image on it ---------------
export CELLPOSE_LOCAL_MODELS_PATH="${TOOLS_DIR}/cellpose_models"
"${VENV_PY}" - <<'PY'
import time, numpy as np, torch
from cellpose import models
img = np.zeros((256, 256), np.float32)
yy, xx = np.mgrid[:256, :256]
for cy, cx in ((60, 60), (60, 120), (150, 90), (200, 200)):
    img[(yy - cy) ** 2 + (xx - cx) ** 2 < 15 ** 2] = 1.0
m = models.CellposeModel(gpu=True, model_type="cyto3")
t0 = time.time(); masks, _, _ = m.eval(img, diameter=30, channels=[0, 0]); dt = time.time() - t0
assert m.gpu and torch.cuda.is_available(), "Cellpose did not select the GPU"
print(f"[CHECK] Cellpose on {torch.cuda.get_device_name(0)}: {int(masks.max())} objects (expected 4) in {dt:.2f} s")
assert int(masks.max()) == 4
PY

# --- 4. Activation helper + manifest ---------------------------------------------------------------
cat <<EOF > "${GPU_HELPER}"
#!/bin/bash
# Generated by setup_gpu_env.bash on $(date --iso-8601=seconds). Source inside GPU jobs:
#   source "${GPU_HELPER}"
# Loads the normal pipeline helper (Java, bftools, ilastik, MITO_* paths), then switches Python
# to the CUDA-enabled venv. Stage 4B picks the GPU automatically when engine: cellpose.
source "${TOOLS_DIR}/activate_mito4d.sh"
source "${ENV_PREFIX}/bin/activate"
EOF
chmod +x "${GPU_HELPER}"

{
    echo "installed_utc=$(date -u --iso-8601=seconds)"
    echo "host=$(hostname)"
    echo "gpu=${GPU_LINE}"
    echo "python=$(${VENV_PY} --version 2>&1)"
    echo "torch_candidate=${TORCH_INDEX} from https://download.pytorch.org/whl/"
    echo "torch=$(${VENV_PY} -c 'import torch; print(torch.__version__, "cuda", torch.version.cuda)')"
    echo "torch_arch_list=$(${VENV_PY} -c 'import torch; print(" ".join(torch.cuda.get_arch_list()))')"
    echo "required_compute_capabilities=${REQUIRE_CC}"
    echo "cellpose=$(${VENV_PY} -c 'import importlib.metadata as m; print(m.version("cellpose"))')"
    echo "pinned_from=${CPU_FREEZE}"
} > "${TOOLS_DIR}/gpu_install_manifest.txt"
echo "[DONE] GPU env ready. Manifest:"
cat "${TOOLS_DIR}/gpu_install_manifest.txt"
