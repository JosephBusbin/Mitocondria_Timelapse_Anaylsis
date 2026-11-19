#!/bin/bash
# ==============================================================================
# setup_hpc_env.sh
# Standardized User-Space Environment & Tool Installation Script
# Target System: William & Mary SciClone HPC (run via salloc or on dev.sciclone.wm.edu)
# ==============================================================================
# Policy Compliance (AGENTS.md):
#   - Uses W&M Lmod modules: python/3.12.7, miniforge3/24.9.2-0, uv/0.13.1
#   - NEVER runs system-level installations (sudo / apt / yum)
#   - NEVER runs 'conda init' (uses 'eval "$(conda shell.bash hook)"')
#   - Places package/pip/uv/conda caches in local node storage (/local/$USER)
#     rather than shared NFS $HOME to avoid quota and I/O bottlenecks
#   - Saves with strict Unix LF line endings
#
# Usage:
#   bash setup_hpc_env.sh [--install-dir "/sciclone/data10/$USER/tools"] [--dry-run]
# ==============================================================================

set -euo pipefail

CURRENT_USER="${USER:-$(whoami)}"
INSTALL_DIR="${PROJECT_TOOLS_DIR:-/sciclone/data10/${CURRENT_USER}/mito_pipeline_tools}"
ENV_NAME="mito4d_env"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-dir)
            INSTALL_DIR="$2"
            shift 2
            ;;
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            echo "Usage: bash setup_hpc_env.sh [--install-dir <DIR>] [--env-name <NAME>] [--dry-run]"
            echo "  --install-dir : Target directory on project storage (default: /sciclone/data10/\$USER/mito_pipeline_tools)"
            echo "  --env-name    : Conda/uv environment name (default: mito4d_env)"
            echo "  --dry-run     : Print planned actions without executing installations"
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------------------
# 2. Configure Local Node Cache & Temporary Directories (/local/$USER)
# ------------------------------------------------------------------------------
# Per W&M RC policy, place caches and temporary downloads in /local/$USER
# rather than shared NFS $HOME. Clean up temporary download folder on exit.
LOCAL_BASE="/local/${CURRENT_USER}"
if [[ ! -d "${LOCAL_BASE}" ]] || [[ ! -w "${LOCAL_BASE}" ]]; then
    LOCAL_BASE="${TMPDIR:-/tmp}/${CURRENT_USER}_local"
fi

TMP_DL_DIR="${LOCAL_BASE}/tmp_install_$$"
export UV_CACHE_DIR="${LOCAL_BASE}/uv_cache"
export PIP_CACHE_DIR="${LOCAL_BASE}/pip_cache"
export CONDA_PKGS_DIRS="${LOCAL_BASE}/conda_pkgs"

# Ensure temporary download artifacts are cleaned up when script exits
trap 'rm -rf "${TMP_DL_DIR}"' EXIT

echo "======================================================================"
echo "W&M SciClone HPC Environment Setup for 4D Mitochondrial Pipeline"
echo "Host             : $(hostname)"
echo "Target Tools Dir : ${INSTALL_DIR}"
echo "Environment Path : ${INSTALL_DIR}/envs/${ENV_NAME}"
echo "Local Cache Dir  : ${LOCAL_BASE}"
echo "Dry Run Mode     : ${DRY_RUN}"
echo "======================================================================"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN] Would load Lmod modules: python/3.12.7, miniforge3/24.9.2-0, uv/0.13.1"
    echo "[DRY-RUN] Would create environment at: ${INSTALL_DIR}/envs/${ENV_NAME}"
    echo "[DRY-RUN] Would install Python packages via uv: numpy, scipy, pandas, pyyaml, tifffile, scikit-image, scikit-fmm, networkx, opencv-python-headless, cellpose"
    echo "[DRY-RUN] Would clone MitoLandscape repository into: ${INSTALL_DIR}/MitoLandscape"
    exit 0
fi

mkdir -p "${TMP_DL_DIR}" "${UV_CACHE_DIR}" "${PIP_CACHE_DIR}" "${CONDA_PKGS_DIRS}"
mkdir -p "${INSTALL_DIR}/envs"

# ------------------------------------------------------------------------------
# 3. Load W&M Lmod Software Modules
# ------------------------------------------------------------------------------
if command -v module >/dev/null 2>&1; then
    echo "[INFO] Loading W&M Lmod modules..."
    module load python/3.12.7 || echo "[WARN] module python/3.12.7 not found, continuing..."
    module load miniforge3/24.9.2-0 || echo "[WARN] module miniforge3/24.9.2-0 not found, continuing..."
    module load uv/0.13.1 || echo "[WARN] module uv/0.13.1 not found, continuing..."
else
    echo "[WARN] 'module' command not detected in current shell; using PATH binaries."
fi

# ------------------------------------------------------------------------------
# 4. Initialize Conda Subshell Hook (NEVER run 'conda init')
# ------------------------------------------------------------------------------
ENV_PREFIX="${INSTALL_DIR}/envs/${ENV_NAME}"

if command -v conda >/dev/null 2>&1; then
    echo "[INFO] Hooking conda into current subshell..."
    eval "$(conda shell.bash hook)"

    if [[ ! -d "${ENV_PREFIX}" ]]; then
        echo "[INFO] Creating conda environment at ${ENV_PREFIX} (Python 3.12)..."
        # -y: automatic confirmation; --prefix: install in project directory instead of $HOME/.conda
        conda create --prefix "${ENV_PREFIX}" python=3.12 -y
    else
        echo "[INFO] Existing conda environment found at ${ENV_PREFIX}"
    fi
    conda activate "${ENV_PREFIX}"
elif command -v uv >/dev/null 2>&1; then
    echo "[INFO] Creating virtual environment via uv at ${ENV_PREFIX}..."
    uv venv "${ENV_PREFIX}" --python 3.12
    # shellcheck disable=SC1091
    source "${ENV_PREFIX}/bin/activate"
else
    echo "[INFO] Creating standard Python venv at ${ENV_PREFIX}..."
    python3 -m venv "${ENV_PREFIX}"
    # shellcheck disable=SC1091
    source "${ENV_PREFIX}/bin/activate"
fi

# ------------------------------------------------------------------------------
# 5. Install Python Dependencies via uv (Fast Resolver) or pip
# ------------------------------------------------------------------------------
PYTHON_PKGS=(
    "numpy>=1.26.0"
    "scipy>=1.12.0"
    "pandas>=2.2.0"
    "pyyaml>=6.0.1"
    "tifffile>=2024.2.12"
    "imagecodecs>=2024.1.1"
    "pillow>=10.2.0"
    "opencv-python-headless>=4.9.0"
    "scikit-image>=0.22.0"
    "scikit-fmm>=2023.4.2"
    "networkx>=3.2.1"
    "cellpose>=3.0.0"
)

if command -v uv >/dev/null 2>&1; then
    echo "[INFO] Installing Python scientific & segmentation packages via uv..."
    uv pip install "${PYTHON_PKGS[@]}"
else
    echo "[INFO] Installing Python scientific & segmentation packages via pip..."
    python3 -m pip install --upgrade pip
    python3 -m pip install "${PYTHON_PKGS[@]}"
fi

# ------------------------------------------------------------------------------
# 6. Clone or Update MitoLandscape Repository (Negri et al., 2025)
# ------------------------------------------------------------------------------
MITOLANDSCAPE_DIR="${INSTALL_DIR}/MitoLandscape"
if [[ ! -d "${MITOLANDSCAPE_DIR}" ]]; then
    if command -v git >/dev/null 2>&1; then
        echo "[INFO] Cloning MitoLandscape repository into ${MITOLANDSCAPE_DIR}..."
        git clone --depth 1 https://github.com/enricoenne/MitoLandscape.git "${MITOLANDSCAPE_DIR}"
    fi
else
    echo "[INFO] MitoLandscape repository already present at ${MITOLANDSCAPE_DIR}"
fi

# ------------------------------------------------------------------------------
# 7. Download Bio-Formats Command Line Tools (bftools) if Absent
# ------------------------------------------------------------------------------
BFTOOLS_DIR="${INSTALL_DIR}/bftools"
if [[ ! -d "${BFTOOLS_DIR}" ]]; then
    if command -v curl >/dev/null 2>&1 && command -v unzip >/dev/null 2>&1; then
        echo "[INFO] Downloading Bio-Formats CLI tools (bftools) to ${BFTOOLS_DIR}..."
        # -L: follow redirects; -sS: silent but show errors; -o: output archive path
        curl -LsS "https://downloads.openmicroscopy.org/bio-formats/7.3.0/artifacts/bftools.zip" \
            -o "${TMP_DL_DIR}/bftools.zip"
        unzip -q "${TMP_DL_DIR}/bftools.zip" -d "${INSTALL_DIR}"
        chmod +x "${BFTOOLS_DIR}/showinf" "${BFTOOLS_DIR}/bfconvert"
    fi
else
    echo "[INFO] Bio-Formats CLI tools already present at ${BFTOOLS_DIR}"
fi

# ------------------------------------------------------------------------------
# 8. Write Environment Activation Helper for Compute Scripts
# ------------------------------------------------------------------------------
ACTIVATE_HELPER="${INSTALL_DIR}/activate_mito4d.sh"
cat <<EOF > "${ACTIVATE_HELPER}"
#!/bin/bash
# Source this file inside SLURM jobs or interactive allocations:
#   source "${ACTIVATE_HELPER}"
if command -v module >/dev/null 2>&1; then
    module load python/3.12.7 2>/dev/null || true
    module load miniforge3/24.9.2-0 2>/dev/null || true
    module load uv/0.13.1 2>/dev/null || true
fi
if command -v conda >/dev/null 2>&1; then
    eval "\$(conda shell.bash hook)"
    conda activate "${ENV_PREFIX}"
elif [[ -f "${ENV_PREFIX}/bin/activate" ]]; then
    source "${ENV_PREFIX}/bin/activate"
fi
export PATH="${BFTOOLS_DIR}:\${PATH}"
export PYTHONPATH="${MITOLANDSCAPE_DIR}:\${PYTHONPATH:-}"
export UV_CACHE_DIR="/local/\${USER:-\$(whoami)}/uv_cache"
export PIP_CACHE_DIR="/local/\${USER:-\$(whoami)}/pip_cache"
EOF
chmod +x "${ACTIVATE_HELPER}"

echo "======================================================================"
echo "[SUCCESS] Environment setup complete!"
echo "To activate in an interactive compute allocation (salloc):"
echo "  source ${ACTIVATE_HELPER}"
echo "======================================================================"