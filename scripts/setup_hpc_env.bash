#!/bin/bash
# ==============================================================================
# setup_hpc_env.bash
# Standardized User-Space Environment & Tool Installation Script
# Target System: William & Mary SciClone HPC (run inside a compute allocation)
# ==============================================================================
# Policy Compliance (AGENTS.md):
#   - Uses W&M Lmod modules: python/3.12.7, uv/0.13.1
#   - NEVER runs system-level installations (sudo / apt / yum) or 'conda init'
#   - Downloads and unpacking (~1 GB ilastik, ~1 GB PyTorch CPU) are heavy, so the script
#     refuses to run on dev.sciclone.wm.edu unless --allow-dev is given. Run it with:
#       module load slurm/bora
#       srun -N1 -c4 --mem=16G -t 01:00:00 bash scripts/setup_hpc_env.bash
#   - Caches and temporary downloads go to node-local storage (/local/$USER or $TMPDIR),
#     not shared NFS $HOME
#   - Saves with strict Unix LF line endings
#
# Where things are installed:
#   Everything goes under $MITO_TOOLS_DIR (default: $MITO_DATA_ROOT/tools, see mito_env.bash):
#     envs/<env-name>/                Python 3.12.7 virtual environment (uv)
#     bftools/                        Bio-Formats CLI tools (showinf / bfconvert)
#     readers/SlideBook6Reader.jar    3i SlideBook 6 reader (kept outside bftools/ so upgrades keep it)
#     ilastik-<ver>-Linux/ + ilastik  ilastik headless binary (symlink 'ilastik' -> versioned dir)
#     MitoLandscape/                  Pinned reference clone (Negri et al., 2025)
#     jre/                            Portable Java 17, only if no usable system Java exists
#     activate_mito4d.sh              Source this in jobs to get everything above on PATH
#     install_manifest.txt, env_freeze.txt   Exact versions installed (reproducibility record)
#
# Usage:
#   bash setup_hpc_env.bash [--tools-dir DIR] [--env-name NAME] [--bftools-version VER]
#                           [--ilastik-version VER] [--mitolandscape-ref SHA]
#                           [--sb-reader-url URL] [--skip-ilastik] [--allow-dev] [--dry-run]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CURRENT_USER="${USER:-$(whoami)}"

ENV_NAME="mito4d_env"
BFTOOLS_VERSION="8.5.0"
ILASTIK_VERSION="1.4.2"          # .ilp projects should be trained in the same ilastik version
MITOLANDSCAPE_REF="45f32305265c302ccba4f65ca4f65668b993027d"
# Pinned to the 2025-07-30 build: the newer 2026-09-03 build's Linux libSBReadFile.so fails to load
# ("undefined symbol: _ZTI15CTileAnnotation"), so Bio-Formats silently falls back to its built-in
# SlideBook reader, which crashes on SlideBook 6 files. Verified on bora 2026-11-23 (49 series read).
# The Windows lab PC (setup_lab_windows.ps1) uses SBReadFile.dll from the newer build, which works.
SB_READER_URL="https://sites.imagej.net/SlideBook/jars/bio-formats/SlideBook6Reader.jar-20250730183505"
TEMURIN_JRE_URL="https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jre/hotspot/normal/eclipse"
SKIP_ILASTIK=0
ALLOW_DEV=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tools-dir)          export MITO_TOOLS_DIR="$2"; shift 2 ;;
        --install-dir)        export MITO_TOOLS_DIR="$2"; shift 2 ;;   # backwards-compatible alias
        --env-name)           ENV_NAME="$2"; shift 2 ;;
        --bftools-version)    BFTOOLS_VERSION="$2"; shift 2 ;;
        --ilastik-version)    ILASTIK_VERSION="$2"; shift 2 ;;
        --mitolandscape-ref)  MITOLANDSCAPE_REF="$2"; shift 2 ;;
        --sb-reader-url)      SB_READER_URL="$2"; shift 2 ;;
        --skip-ilastik)       SKIP_ILASTIK=1; shift ;;
        --allow-dev)          ALLOW_DEV=1; shift ;;
        --dry-run)            DRY_RUN=1; shift ;;
        -h|--help)
            sed -n '2,35p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------------------
# 1. Resolve storage paths (single source of truth: mito_env.bash)
# ------------------------------------------------------------------------------
# shellcheck source=mito_env.bash
source "${SCRIPT_DIR}/mito_env.bash"
TOOLS_DIR="${MITO_TOOLS_DIR}"
ENV_PREFIX="${TOOLS_DIR}/envs/${ENV_NAME}"
BFTOOLS_DIR="${TOOLS_DIR}/bftools"
READERS_DIR="${TOOLS_DIR}/readers"
SB_JAR="${READERS_DIR}/SlideBook6Reader.jar"
ILASTIK_DIR="${TOOLS_DIR}/ilastik-${ILASTIK_VERSION}-Linux"
ILASTIK_LINK="${TOOLS_DIR}/ilastik"
MITOLANDSCAPE_DIR="${TOOLS_DIR}/MitoLandscape"
JRE_DIR="${TOOLS_DIR}/jre"
ACTIVATE_HELPER="${TOOLS_DIR}/activate_mito4d.sh"

# ------------------------------------------------------------------------------
# 2. Node-local cache & temporary directories (cleaned up on exit)
# ------------------------------------------------------------------------------
LOCAL_BASE="/local/${CURRENT_USER}"
if [[ ! -d "${LOCAL_BASE}" ]] || [[ ! -w "${LOCAL_BASE}" ]]; then
    LOCAL_BASE="${TMPDIR:-/tmp}/${CURRENT_USER}_local"
fi
TMP_DL_DIR="${LOCAL_BASE}/tmp_install_$$"
export UV_CACHE_DIR="${LOCAL_BASE}/uv_cache"
export PIP_CACHE_DIR="${LOCAL_BASE}/pip_cache"
trap 'rm -rf "${TMP_DL_DIR}"' EXIT

echo "======================================================================"
echo "W&M SciClone HPC Environment Setup for 4D Mitochondrial Pipeline"
echo "Host             : $(hostname)"
echo "MITO_DATA_ROOT   : ${MITO_DATA_ROOT}"
echo "Tools Dir        : ${TOOLS_DIR}"
echo "Python Env       : ${ENV_PREFIX}"
echo "Bio-Formats      : ${BFTOOLS_VERSION} (+ SlideBook6Reader.jar in ${READERS_DIR})"
echo "ilastik          : $([[ ${SKIP_ILASTIK} -eq 1 ]] && echo skipped || echo "${ILASTIK_VERSION}")"
echo "MitoLandscape    : ${MITOLANDSCAPE_REF}"
echo "Local Cache Dir  : ${LOCAL_BASE}"
echo "Dry Run Mode     : ${DRY_RUN}"
echo "======================================================================"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY-RUN] Would load modules python/3.12.7 + uv/0.13.1 and create a uv venv at ${ENV_PREFIX}"
    echo "[DRY-RUN] Would install CPU-only PyTorch, then numpy/scipy/pandas/pyyaml/tifffile/imagecodecs/pillow/"
    echo "          opencv-python-headless/scikit-image/scikit-fmm/networkx and cellpose>=3.1,<4 (cyto3 model)"
    echo "[DRY-RUN] Would install bftools ${BFTOOLS_VERSION} to ${BFTOOLS_DIR} and SlideBook6Reader.jar to ${SB_JAR}"
    [[ ${SKIP_ILASTIK} -eq 0 ]] && echo "[DRY-RUN] Would install ilastik ${ILASTIK_VERSION} to ${ILASTIK_DIR} (symlink ${ILASTIK_LINK})"
    echo "[DRY-RUN] Would clone MitoLandscape @ ${MITOLANDSCAPE_REF} into ${MITOLANDSCAPE_DIR}"
    echo "[DRY-RUN] Would write ${ACTIVATE_HELPER}, ${TOOLS_DIR}/install_manifest.txt and env_freeze.txt"
    exit 0
fi

# Refuse heavy installs on the shared dev server (W&M RC policy: dev is not a compute node).
if [[ "$(hostname -s)" == dev* ]] && [[ "${ALLOW_DEV}" -eq 0 ]]; then
    echo "[ERROR] Refusing to run heavy installs on $(hostname -s). Run inside a compute allocation:" >&2
    echo "        module load slurm/bora && srun -N1 -c4 --mem=16G -t 01:00:00 bash ${BASH_SOURCE[0]}" >&2
    exit 1
fi

mkdir -p "${TMP_DL_DIR}" "${UV_CACHE_DIR}" "${PIP_CACHE_DIR}" "${TOOLS_DIR}/envs" "${READERS_DIR}" \
         "${TOOLS_DIR}/cellpose_models"
# Cache and venv sit on different filesystems, so uv cannot hardlink; copy silently instead.
export UV_LINK_MODE=copy

# ------------------------------------------------------------------------------
# 3. Load W&M Lmod modules
# ------------------------------------------------------------------------------
if ! command -v module >/dev/null 2>&1 && [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/Modules/init/bash
fi
# Module init scripts reference unset variables, so relax 'set -u' just around them.
set +u
module load python/3.12.7
module load uv/0.13.1
set -u
PYTHON_BIN="$(command -v python3)"
echo "[INFO] Python interpreter: ${PYTHON_BIN} ($(${PYTHON_BIN} --version 2>&1))"

# ------------------------------------------------------------------------------
# 4. Java runtime for Bio-Formats (system Java >= 8, else portable Temurin 17)
# ------------------------------------------------------------------------------
JAVA_HOME_USED=""
java_major() {
    # 'java -version' prints e.g. 'openjdk version "1.8.0_504"' or '"17.0.13"'; 1.x means major x.
    local v
    v="$("$1" -version 2>&1 | awk -F'"' '/version/ {print $2; exit}')"
    [[ -z "${v}" ]] && { echo 0; return; }
    if [[ "${v}" == 1.* ]]; then echo "${v}" | cut -d. -f2; else echo "${v%%.*}"; fi
}
if [[ -x "${JRE_DIR}/bin/java" ]]; then
    JAVA_HOME_USED="${JRE_DIR}"
elif command -v java >/dev/null 2>&1 && [[ "$(java_major java)" -ge 8 ]]; then
    JAVA_HOME_USED="$(dirname "$(dirname "$(readlink -f "$(command -v java)")")")"
else
    echo "[INFO] No usable Java >= 8 found; downloading portable Temurin JRE 17..."
    # -L follow redirects; -f fail on HTTP errors; -sS quiet but show errors
    curl -LfsS "${TEMURIN_JRE_URL}" -o "${TMP_DL_DIR}/jre17.tar.gz"
    mkdir -p "${JRE_DIR}"
    # --strip-components=1 drops the versioned top-level folder so java lives at jre/bin/java
    tar -xzf "${TMP_DL_DIR}/jre17.tar.gz" -C "${JRE_DIR}" --strip-components=1
    JAVA_HOME_USED="${JRE_DIR}"
fi
export JAVA_HOME="${JAVA_HOME_USED}"
export PATH="${JAVA_HOME}/bin:${PATH}"
JAVA_VERSION_LINE="$("${JAVA_HOME}/bin/java" -version 2>&1 | head -1)"
echo "[INFO] Java: ${JAVA_HOME} (${JAVA_VERSION_LINE})"

# ------------------------------------------------------------------------------
# 5. Bio-Formats CLI tools (version-checked, upgraded in place) + 3i SlideBook 6 reader
# ------------------------------------------------------------------------------
# SlideBook6Reader.source records which URL the installed jar came from, so changing the pin
# (--sb-reader-url) replaces the jar instead of silently keeping an old one.
SB_SOURCE_FILE="${READERS_DIR}/SlideBook6Reader.source"
if [[ ! -f "${SB_JAR}" ]] || [[ "$(cat "${SB_SOURCE_FILE}" 2>/dev/null)" != "${SB_READER_URL}" ]]; then
    echo "[INFO] Downloading 3i SlideBook6Reader.jar (${SB_READER_URL##*-}) -> ${SB_JAR}"
    curl -LfsS "${SB_READER_URL}" -o "${TMP_DL_DIR}/SlideBook6Reader.jar"
    mv "${TMP_DL_DIR}/SlideBook6Reader.jar" "${SB_JAR}"
    echo "${SB_READER_URL}" > "${SB_SOURCE_FILE}"
fi
export BF_CP="${SB_JAR}"

NEEDS_BF=1
if [[ -x "${BFTOOLS_DIR}/showinf" ]] && "${BFTOOLS_DIR}/showinf" -version 2>&1 | grep -qF "${BFTOOLS_VERSION}"; then
    NEEDS_BF=0
    echo "[INFO] bftools ${BFTOOLS_VERSION} already present at ${BFTOOLS_DIR}"
fi
if [[ "${NEEDS_BF}" -eq 1 ]]; then
    echo "[INFO] Installing bftools ${BFTOOLS_VERSION} -> ${BFTOOLS_DIR}"
    # Download first, then replace, so a failed download never removes a working install.
    curl -LfsS "https://downloads.openmicroscopy.org/bio-formats/${BFTOOLS_VERSION}/artifacts/bftools.zip" \
        -o "${TMP_DL_DIR}/bftools.zip"
    unzip -q "${TMP_DL_DIR}/bftools.zip" -d "${TMP_DL_DIR}/bf_unpack"
    rm -rf "${BFTOOLS_DIR}"
    mv "${TMP_DL_DIR}/bf_unpack/bftools" "${BFTOOLS_DIR}"
fi
BF_VERSION_LINE="$("${BFTOOLS_DIR}/showinf" -version 2>&1 | grep -m1 -E 'Version' || true)"
if ! grep -qF "${BFTOOLS_VERSION}" <<< "${BF_VERSION_LINE}"; then
    echo "[ERROR] bftools verification failed (expected ${BFTOOLS_VERSION}): ${BF_VERSION_LINE}" >&2
    exit 1
fi
echo "[INFO] Bio-Formats: ${BF_VERSION_LINE}"

# Verify 3i's native Linux library actually loads here; if it does not, Bio-Formats silently falls
# back to its built-in SlideBook reader (which crashes on SlideBook 6 files). dlopen via ctypes
# resolves the same symbols Java's System.load does, so a broken build fails this check.
unzip -q -o "${SB_JAR}" 'META-INF/lib/linux_64/libSBReadFile.so' -d "${TMP_DL_DIR}/sb6_native"
if ! SB_LOAD_ERR="$("${PYTHON_BIN}" -c 'import ctypes, sys; ctypes.CDLL(sys.argv[1])' \
        "${TMP_DL_DIR}/sb6_native/META-INF/lib/linux_64/libSBReadFile.so" 2>&1)"; then
    echo "[ERROR] SlideBook6Reader native library failed to load on $(hostname -s):" >&2
    echo "        ${SB_LOAD_ERR##*: }" >&2
    echo "        Pin a different build with --sb-reader-url (https://sites.imagej.net/SlideBook/jars/bio-formats/)." >&2
    exit 1
fi
echo "[INFO] SlideBook6Reader native library loads: ${SB_READER_URL##*/}"

# ------------------------------------------------------------------------------
# 6. ilastik headless (pinned release; .ilp files must be trained in the same version)
# ------------------------------------------------------------------------------
if [[ "${SKIP_ILASTIK}" -eq 0 ]]; then
    if [[ ! -x "${ILASTIK_DIR}/run_ilastik.sh" ]]; then
        echo "[INFO] Downloading ilastik ${ILASTIK_VERSION} (~820 MB) to node-local storage..."
        curl -LfsS "https://files.ilastik.org/ilastik-${ILASTIK_VERSION}-Linux.tar.bz2" \
            -o "${TMP_DL_DIR}/ilastik.tar.bz2"
        echo "[INFO] Unpacking ilastik into ${TOOLS_DIR} ..."
        tar -xjf "${TMP_DL_DIR}/ilastik.tar.bz2" -C "${TOOLS_DIR}"
    else
        echo "[INFO] ilastik ${ILASTIK_VERSION} already present at ${ILASTIK_DIR}"
    fi
    # -s symbolic, -f replace existing, -n treat an existing symlink-to-dir as a file
    ln -sfn "${ILASTIK_DIR}" "${ILASTIK_LINK}"
    [[ -x "${MITO_ILASTIK_BIN}" ]] || { echo "[ERROR] ilastik binary missing: ${MITO_ILASTIK_BIN}" >&2; exit 1; }
fi

# ------------------------------------------------------------------------------
# 7. Python 3.12.7 virtual environment (uv) with CPU-only PyTorch
# ------------------------------------------------------------------------------
if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
    echo "[INFO] Creating uv venv at ${ENV_PREFIX} from ${PYTHON_BIN}"
    uv venv "${ENV_PREFIX}" --python "${PYTHON_BIN}"
fi
VENV_PY="${ENV_PREFIX}/bin/python"

# bora has no GPUs: the CPU wheel avoids ~3 GB of CUDA libraries from the default PyPI torch.
echo "[INFO] Installing CPU-only PyTorch..."
uv pip install --python "${VENV_PY}" --index-url https://download.pytorch.org/whl/cpu torch

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
    "cellpose>=3.1,<4"      # Cellpose 4 (Cellpose-SAM) dropped the 'cyto3' model used in pipeline_config.yaml
)
echo "[INFO] Installing scientific & segmentation packages..."
uv pip install --python "${VENV_PY}" "${PYTHON_PKGS[@]}"
uv pip freeze --python "${VENV_PY}" > "${TOOLS_DIR}/env_freeze.txt"

# ------------------------------------------------------------------------------
# 8. MitoLandscape reference clone, pinned to a commit
# ------------------------------------------------------------------------------
if [[ ! -d "${MITOLANDSCAPE_DIR}/.git" ]]; then
    git clone --quiet https://github.com/enricoenne/MitoLandscape.git "${MITOLANDSCAPE_DIR}"
fi
git -C "${MITOLANDSCAPE_DIR}" fetch --quiet origin
git -C "${MITOLANDSCAPE_DIR}" -c advice.detachedHead=false checkout --quiet "${MITOLANDSCAPE_REF}"

# ------------------------------------------------------------------------------
# 9. Activation helper + reproducibility manifest
# ------------------------------------------------------------------------------
cat <<EOF > "${ACTIVATE_HELPER}"
#!/bin/bash
# Generated by setup_hpc_env.bash on $(date --iso-8601=seconds). Source inside jobs or srun shells:
#   source "${ACTIVATE_HELPER}"
if ! command -v module >/dev/null 2>&1 && [[ -f /usr/share/Modules/init/bash ]]; then
    source /usr/share/Modules/init/bash
fi
module load python/3.12.7 2>/dev/null || true
source "${ENV_PREFIX}/bin/activate"
export MITO_TOOLS_DIR="\${MITO_TOOLS_DIR:-${TOOLS_DIR}}"
source "${SCRIPT_DIR}/mito_env.bash"
# Java recorded at install time; system Java differs between dev (1.8) and compute nodes (11),
# so fall back to whatever 'java' is on PATH if that directory is absent on this node.
if [[ -x "${JAVA_HOME_USED}/bin/java" ]]; then
    export JAVA_HOME="${JAVA_HOME_USED}"
    export PATH="${JAVA_HOME_USED}/bin:\${PATH}"
elif command -v java >/dev/null 2>&1; then
    export JAVA_HOME="\$(dirname "\$(dirname "\$(readlink -f "\$(command -v java)")")")"
fi
export BF_JAVA_HOME="\${JAVA_HOME:-}"
export BF_CP="${SB_JAR}"
export PATH="${BFTOOLS_DIR}:\${PATH}"
# Cellpose downloads model weights on first use; keep them off the shared \$HOME quota.
export CELLPOSE_LOCAL_MODELS_PATH="${TOOLS_DIR}/cellpose_models"
export UV_CACHE_DIR="/local/\${USER:-\$(whoami)}/uv_cache"
export PIP_CACHE_DIR="/local/\${USER:-\$(whoami)}/pip_cache"
EOF
chmod +x "${ACTIVATE_HELPER}"

{
    echo "installed_utc=$(date -u --iso-8601=seconds)"
    echo "host=$(hostname)"
    echo "python=$(${VENV_PY} --version 2>&1)"
    echo "java_home=${JAVA_HOME_USED}"
    echo "java=${JAVA_VERSION_LINE}"
    echo "bftools=${BF_VERSION_LINE}"
    echo "slidebook6reader_url=${SB_READER_URL}"
    echo "slidebook6reader_sha256=$(sha256sum "${SB_JAR}" | cut -d' ' -f1)"
    echo "ilastik=$([[ ${SKIP_ILASTIK} -eq 1 ]] && echo skipped || echo "${ILASTIK_VERSION}")"
    echo "mitolandscape_commit=$(git -C "${MITOLANDSCAPE_DIR}" rev-parse HEAD)"
    echo "cellpose=$(${VENV_PY} -c 'import cellpose, importlib.metadata as m; print(m.version("cellpose"))')"
    echo "torch=$(${VENV_PY} -c 'import torch; print(torch.__version__)')"
} > "${TOOLS_DIR}/install_manifest.txt"

echo "======================================================================"
echo "[SUCCESS] Environment setup complete!"
cat "${TOOLS_DIR}/install_manifest.txt"
echo "Activate in jobs or srun shells:"
echo "  source ${ACTIVATE_HELPER}"
echo "======================================================================"
