#!/bin/bash
# ==============================================================================
# submit_full_run.bash
# Full-archive orchestration for the 4D mitochondrial pipeline on W&M SciClone (bora by default).
# ==============================================================================
# DRY RUN BY DEFAULT: every phase prints the exact sbatch commands (with comments explaining
# resources, throttles and dependencies) and submits NOTHING. Add --execute to submit.
# Run it from dev.sciclone.wm.edu (light work only: it plans and submits, it computes nothing).
#
# Usage:
#   bash scripts/submit_full_run.bash <phase> --run-id ID [options] [--execute]
#
# Phases (in order; see manuals/full_run_playbook.md):
#   scan       A : sharded .sld metadata scan (array, 1 core, %20) + merge (afterany)
#                  [--chain-extract also submits 'extract' with afterok on the merge]
#   merge      A': re-submit only the merge of the scan shards (e.g. after retrying a shard)
#   extract    B : OME-TIFF extraction, one task per qualifying .sld (2 cores, %20, afterok merge)
#   canary     C : Stage 3 on ~10 movies spanning the size classes (measures memory/time)
#   stage3     D : Stage 3 for all movies, one array per size class (model refit from the canary)
#   stage45    F : Stages 4-5, aftercorr on the matching Stage 3 array (needs the retrained .ilp)
#   aggregate  G : cohort tables + QC dashboard + failure list (afterany on Stages 4-5)
#   retry      H : resubmit FAILED / OUT_OF_MEMORY / TIMEOUT tasks with memory/time x1.5 (max 2)
#   status       : one sacct summary of every job in the run (reads the executed manifest)
#   overflow     : report other clusters (deny list, live PreemptMode, smoke test); with
#                  --overflow move not-yet-started chunks there when bora is saturated
#
# Common options:
#   --run-id ID          run name (default $MITO_RUN_ID); state in <run_root>/<ID>/, logs in
#                        /sciclone/scr30/$USER/mito_logs/<ID>/
#   --config FILE        pipeline_config.yaml (default: repo copy); hpc_run: block = all throttles/model
#   --execute            really submit (otherwise dry run; dry-run state goes to <run>/dryrun/)
#   --raw-root DIR       raw OME-TIFF root (default $MITO_RAW_ROOT)
#   --no-module          do not 'module load slurm/bora' (tests, or Slurm already configured)
# Phase options:
#   scan:     --source-dir DIR [--file-list F] [--plasmid-csv CSV] [--mito-overrides CSV] [--chain-extract]
#   extract:  [--source-dir DIR] [--labels MITO,REVIEW]   (--labels also applies to scan --chain-extract,
#             merge, canary and stage3; default hpc_run.collector.extract_labels)
#   canary/stage3: [--movie-list F] [--inventory CSV] [--no-refit] [--sacct-file F]
#   stage45:  [--allow-fallback]   (test without a trained .ilp; built-in estimator)
#   retry:    --phase scan|extract|canary|stage3|stage45 [--sacct-file F]
#   overflow: [--overflow] [--smoke]  (--smoke --execute runs the per-cluster env smoke test)
#
# Internal (job bodies; sbatch runs this same script): __scan __merge __extract __aggregate __smoke
# ==============================================================================
set -euo pipefail

die() { echo "[ERROR] $*" >&2; exit 1; }

# ------------------------------------------------------------------------------
# Locate the repository. Jobs get MITO_REPO_DIR through --export (sbatch runs a spooled copy).
# ------------------------------------------------------------------------------
if [[ -n "${MITO_REPO_DIR:-}" ]]; then
    SCRIPT_DIR="${MITO_REPO_DIR}/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
SELF="${SCRIPT_DIR}/submit_full_run.bash"
PLANNER="${SCRIPT_DIR}/plan_full_run.py"
[[ -f "${PLANNER}" ]] || die "plan_full_run.py not found in ${SCRIPT_DIR} (set MITO_REPO_DIR)"

[[ $# -ge 1 ]] || { sed -n '2,40p' "${SELF}"; exit 1; }
PHASE="$1"; shift

# ==============================================================================
# Part 1: job bodies (run on compute nodes)
# ==============================================================================
activate_env() {
    # Pipeline venv + Java/bftools + MITO_* paths (the helper written by setup_hpc_env.bash).
    # shellcheck source=mito_env.bash
    source "${SCRIPT_DIR}/mito_env.bash"
    local helper="${MITO_TOOLS_DIR}/activate_mito4d.sh"
    [[ -f "${helper}" ]] || die "${helper} not found; run scripts/setup_hpc_env.bash first"
    set +u   # module/venv activation scripts reference unset variables
    # shellcheck disable=SC1090
    source "${helper}"
    set -u
    # Temporary files go to a per-task folder on scr30, never /tmp or /local (both 10 GB).
    TASK_TMPDIR="${MITO_SCRATCH_ROOT}/tmp/${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-pid$$}}_${SLURM_ARRAY_TASK_ID:-0}"
    mkdir -p "${TASK_TMPDIR}"
    export TMPDIR="${TASK_TMPDIR}"
    trap 'rm -rf "${TASK_TMPDIR}"' EXIT
}

parse_task_args() {
    T_RUN_DIR=""; T_SOURCE=""; T_LIST=""; T_SHARDS=""; T_RAW=""; T_CONFIG=""; T_PLASMID=""; T_PHASE=""; T_OVERRIDES=""; T_LABELS=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --run-dir)     T_RUN_DIR="$2"; shift 2 ;;
            --source-dir)  T_SOURCE="$2"; shift 2 ;;
            --file-list)   T_LIST="$2"; shift 2 ;;
            --shards)      T_SHARDS="$2"; shift 2 ;;
            --raw-root)    T_RAW="$2"; shift 2 ;;
            --config)      T_CONFIG="$2"; shift 2 ;;
            --plasmid-csv) T_PLASMID="$2"; shift 2 ;;
            --mito-overrides) T_OVERRIDES="$2"; shift 2 ;;
            --labels)      T_LABELS="$2"; shift 2 ;;
            --phase)       T_PHASE="$2"; shift 2 ;;
            *) die "unknown task option $1" ;;
        esac
    done
}

case "${PHASE}" in
    __scan)
        parse_task_args "$@"
        activate_env
        K="${SLURM_ARRAY_TASK_ID:?__scan must run as an array task}"
        OUT="${T_RUN_DIR}/scan/shards"
        mkdir -p "${OUT}"
        read -r -a COLLECTOR_ARGS <<< "$(python3 "${PLANNER}" collector-args --config "${T_CONFIG}")"
        FILTER_ARGS=(--mito-filter off)
        # 'report' labels every file (MITO/REVIEW/...) but still inspects all of them, so the
        # inventory is complete and the label filter is applied at extraction time.
        [[ -n "${T_PLASMID}" ]] && FILTER_ARGS=(--plasmid-csv "${T_PLASMID}" --mito-filter report)
        [[ -n "${T_OVERRIDES}" ]] && FILTER_ARGS+=(--mito-overrides "${T_OVERRIDES}")
        # --workers 1: one showinf at a time per task (the array provides the parallelism and
        # avoids contention in SlideBook's native JNI reader).
        python3 "${SCRIPT_DIR}/collect_4d_sld_timelapses.py" "${T_SOURCE}" "${OUT}" \
            --mode scan-only --file-list "${T_LIST}" --shard "${K}/${T_SHARDS}" --report-dir "${OUT}" \
            --workers 1 --hpc-raw-root "${T_RAW}" "${COLLECTOR_ARGS[@]}" "${FILTER_ARGS[@]}"
        exit 0 ;;
    __extract)
        parse_task_args "$@"
        activate_env
        K="${SLURM_ARRAY_TASK_ID:?__extract must run as an array task}"
        [[ -f "${T_LIST}" ]] || die "extract list ${T_LIST} missing (did the merge succeed?)"
        N_FILES="$(grep -c . "${T_LIST}" || true)"
        if [[ "${K}" -gt "${N_FILES}" ]]; then
            echo "Task ${K}: the extract list has only ${N_FILES} file(s); nothing to do."
            exit 0
        fi
        OUT="${T_RUN_DIR}/extract/shards"
        mkdir -p "${OUT}" "${T_RAW}"
        read -r -a COLLECTOR_ARGS <<< "$(python3 "${PLANNER}" collector-args --config "${T_CONFIG}")"
        # Skip-existing is the collector default (no --overwrite); --checksum writes SHA-256 per file;
        # the shard is taken over the extract list, so task K extracts file K.
        python3 "${SCRIPT_DIR}/collect_4d_sld_timelapses.py" "${T_SOURCE}" "${T_RAW}" \
            --mode extract-ometiff --file-list "${T_LIST}" --shard "${K}/${T_SHARDS}" --report-dir "${OUT}" \
            --workers 1 --checksum --hpc-raw-root "${T_RAW}" --mito-filter off "${COLLECTOR_ARGS[@]}"
        exit 0 ;;
    __merge)
        parse_task_args "$@"
        activate_env
        RC=0
        if [[ "${T_PHASE}" == "scan" ]]; then
            python3 "${SCRIPT_DIR}/merge_collect_shards.py" "${T_RUN_DIR}/scan/shards" "${T_RUN_DIR}/scan/merged" \
                --expected-shards "${T_SHARDS}" --hpc-raw-root "${T_RAW}" || RC=$?
            # The extract list is written even after a partial merge, for inspection; the exit code
            # (3 = missing shards) still stops an 'afterok' extraction.
            python3 "${PLANNER}" extract-list --config "${T_CONFIG}" --run-dir "${T_RUN_DIR}" \
                ${T_LABELS:+--labels "${T_LABELS}"}
        else
            # Tasks beyond the real extract list (chained submission) write no reports: expect one
            # shard per listed file, not per task.
            N_FILES="$(grep -c . "${T_RUN_DIR}/extract/extract_files.txt" || true)"
            python3 "${SCRIPT_DIR}/merge_collect_shards.py" "${T_RUN_DIR}/extract/shards" "${T_RUN_DIR}/extract/merged" \
                --expected-shards "${N_FILES:-${T_SHARDS}}" --hpc-raw-root "${T_RAW}" || RC=$?
        fi
        exit "${RC}" ;;
    __aggregate)
        parse_task_args "$@"
        activate_env
        python3 "${SCRIPT_DIR}/aggregate_results.py" --run-dir "${T_RUN_DIR}" --out-dir "${T_RUN_DIR}/results"
        exit 0 ;;
    __smoke)
        # Per-cluster environment check before any overflow use: imports + ilastik start-up.
        activate_env
        echo "host $(hostname), $(uname -m), $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
        python3 -c "import numpy, scipy, skimage, tifffile, yaml; print('python', __import__('sys').version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'skimage', skimage.__version__)"
        timeout 600 "${MITO_ILASTIK_BIN}" --headless --version
        echo "SMOKE_OK"
        exit 0 ;;
esac

# ==============================================================================
# Part 2: orchestration (runs on dev; plans with plan_full_run.py, then prints or submits)
# ==============================================================================
RUN_ID="${MITO_RUN_ID:-}"
CONFIG="${REPO_DIR}/pipeline_config.yaml"
EXECUTE=0
LOAD_MODULE=1
SMOKE=0
RETRY_PHASE=""
PLAN_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-id)       RUN_ID="$2"; shift 2 ;;
        --config)       CONFIG="$2"; shift 2 ;;
        --execute)      EXECUTE=1; shift ;;
        --no-module)    LOAD_MODULE=0; shift ;;
        --smoke)        SMOKE=1; shift ;;
        --phase)        RETRY_PHASE="$2"; shift 2 ;;
        --raw-root|--source-dir|--file-list|--plasmid-csv|--mito-overrides|--labels|--movie-list|--inventory|--sacct-file)
                        PLAN_ARGS+=("$1" "$2"); shift 2 ;;
        --chain-extract|--no-refit|--allow-fallback|--overflow)
                        PLAN_ARGS+=("$1"); shift ;;
        -h|--help)      sed -n '2,40p' "${SELF}"; exit 0 ;;
        *)              die "unknown option '$1' (see --help)" ;;
    esac
done
[[ -n "${RUN_ID}" ]] || die "--run-id is required (or export MITO_RUN_ID)"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]] || die "--run-id may contain only letters, digits, '.', '_' and '-'"
[[ -f "${CONFIG}" ]] || die "config not found: ${CONFIG}"
CONFIG="$(readlink -f "${CONFIG}")"

# Always talk to the bora cluster (dev's Slurm commands are cluster-agnostic until a module is loaded).
if [[ "${LOAD_MODULE}" -eq 1 ]]; then
    if ! type module >/dev/null 2>&1 && [[ -f /usr/share/Modules/init/bash ]]; then
        source /usr/share/Modules/init/bash
    fi
    module load slurm/bora >/dev/null 2>&1 || die "'module load slurm/bora' failed"
fi
# shellcheck source=mito_env.bash
source "${SCRIPT_DIR}/mito_env.bash"

# Run directory: hpc_run.run_root, else $MITO_SCRATCH_ROOT/../mito_runs (= /sciclone/scr30/$USER/mito_runs)
RUN_ROOT="$(python3 -c 'import sys; sys.path.insert(0, sys.argv[1]); import plan_full_run as p; print(p.load_hpc_run(sys.argv[2]).get("run_root") or "")' "${SCRIPT_DIR}" "${CONFIG}")"
RUN_ROOT="${RUN_ROOT:-$(dirname "${MITO_SCRATCH_ROOT}")/mito_runs}"
RUN_DIR="${RUN_ROOT}/${RUN_ID}"
if [[ "${EXECUTE}" -eq 1 ]]; then
    MANIFEST="${RUN_DIR}/run_manifest.json"
    DRY_FLAG=()
else
    MANIFEST="${RUN_DIR}/dryrun/run_manifest.json"
    DRY_FLAG=(--dry-run)
fi

if [[ "${PHASE}" == "status" ]]; then
    # Status reads the real manifest; it falls back to the dry-run one only if nothing was executed.
    if [[ -f "${RUN_DIR}/run_manifest.json" ]]; then
        python3 "${PLANNER}" status --run-dir "${RUN_DIR}"
    else
        python3 "${PLANNER}" status --run-dir "${RUN_DIR}" --dry-run
    fi
    exit 0
fi

case "${PHASE}" in
    scan|merge|extract|canary|stage3|stage45|aggregate|overflow) PLAN_PHASE="${PHASE}" ;;
    retry)  PLAN_PHASE="retry"
            [[ -n "${RETRY_PHASE}" ]] || die "retry needs --phase scan|extract|canary|stage3|stage45"
            PLAN_ARGS+=(--retry-phase "${RETRY_PHASE}") ;;
    *)      die "unknown phase '${PHASE}' (see --help)" ;;
esac

# --- helpers the generated program calls ----------------------------------------------------------
declare -A JOB_IDS=()
N_SUBMITTED=0

note() { echo "# $*"; }

ensure_dir() {
    if [[ "${EXECUTE}" -eq 1 ]]; then mkdir -p "$@"; else echo "mkdir -p $*"; fi
}

# Prints a command ready to copy-paste: one option per line, an option's separate value on the same
# line. Arguments made only of shell-safe characters print as-is; others are quoted with printf %q.
print_cmd() {
    local first=1 prev_opt=0 a q
    for a in "$@"; do
        if [[ "${a}" =~ ^[A-Za-z0-9_./=:,%@+-]+$ ]]; then q="${a}"; else q="$(printf '%q' "${a}")"; fi
        if [[ "${first}" -eq 1 ]]; then
            printf '%s' "${q}"; first=0
        elif [[ "${prev_opt}" -eq 1 ]] && [[ "${a}" != -* ]]; then
            printf ' %s' "${q}"            # value of the preceding '--option'
        else
            printf ' \\\n    %s' "${q}"
        fi
        # an option without '=' may take the next argument as its value
        if [[ "${a}" == -* ]] && [[ "${a}" != *=* ]]; then prev_opt=1; else prev_opt=0; fi
    done
    printf '\n'
}

# submit [--cluster C] NAME sbatch-args... : '@@other@@' in an argument is replaced by the job id of
# the job NAMEd 'other' submitted earlier in this program (dry run: a DRY_<name> placeholder).
submit() {
    local cluster="" name a ref id
    if [[ "$1" == "--cluster" ]]; then cluster="$2"; shift 2; fi
    name="$1"; shift
    local args=()
    for a in "$@"; do
        while [[ "${a}" =~ @@([A-Za-z0-9_.-]+)@@ ]]; do
            ref="${BASH_REMATCH[1]}"
            id="${JOB_IDS[${ref}]:-}"
            [[ -n "${id}" ]] || die "dependency on ${ref}, which was not submitted"
            a="${a//@@${ref}@@/${id}}"
        done
        args+=("${a}")
    done
    if [[ -n "${cluster}" ]]; then
        echo "# (on cluster ${cluster}: module unload slurm; module load slurm/${cluster})"
    fi
    if [[ "${EXECUTE}" -eq 1 ]]; then
        # --parsable prints 'jobid' or 'jobid;cluster'
        if [[ -n "${cluster}" ]]; then
            id="$( module unload slurm >/dev/null 2>&1; module load "slurm/${cluster}" >/dev/null 2>&1 && sbatch --parsable "${args[@]}" )"
        else
            id="$(sbatch --parsable "${args[@]}")"
        fi
        id="${id%%;*}"
        [[ "${id}" =~ ^[0-9]+$ ]] || die "sbatch for ${name} returned '${id}'"
        print_cmd sbatch "${args[@]}"
        echo "#   -> submitted ${name} = ${id}"
    else
        id="DRY_${name}"
        print_cmd sbatch "${args[@]}"
    fi
    JOB_IDS["${name}"]="${id}"
    N_SUBMITTED=$((N_SUBMITTED + 1))
    python3 "${PLANNER}" record --manifest "${MANIFEST}" --name "${name}" --job-id "${id}"
    echo
}

# cancel CLUSTER ID : used only by 'overflow --overflow' to move a fully pending chunk.
cancel() {
    local cluster="$1" id="$2"
    if [[ "${EXECUTE}" -eq 1 ]]; then
        if [[ -n "${cluster}" ]]; then ( module unload slurm >/dev/null 2>&1; module load "slurm/${cluster}" >/dev/null 2>&1 && scancel "${id}" ); else scancel "${id}"; fi
        echo "# cancelled ${id}"
    else
        echo "scancel ${id}"
    fi
}

# smoke CLUSTER PARTITION : per-cluster environment test (only with --smoke --execute).
smoke() {
    local cluster="$1" part="$2" result="FAIL"
    # -N1 -n1 -c1: one task on one core; --mem/-t: small and short (a check, not a workload)
    local cmd=(srun -p "${part}" -N1 -n1 -c1 --mem=2G -t 00:15:00 --export=ALL,MITO_REPO_DIR="${REPO_DIR}" bash "${SELF}" __smoke)
    if [[ "${EXECUTE}" -eq 1 ]] && [[ "${SMOKE}" -eq 1 ]]; then
        echo "# smoke test on ${cluster}/${part}"
        if ( module unload slurm >/dev/null 2>&1; module load "slurm/${cluster}" >/dev/null 2>&1 && "${cmd[@]}" ) | tee /dev/stderr | grep -q SMOKE_OK; then
            result="PASS"
        fi
        python3 "${PLANNER}" smoke-record --manifest "${MANIFEST}" --cluster "${cluster}" --partition "${part}" --result "${result}"
        echo "# smoke test ${cluster}: ${result}"
    else
        echo "# smoke test (runs with --smoke --execute):"
        echo "module unload slurm; module load slurm/${cluster}"
        print_cmd "${cmd[@]}"
    fi
}

# --- plan and run ---------------------------------------------------------------------------------
if [[ "${EXECUTE}" -eq 1 ]]; then
    echo "# EXECUTE mode: jobs WILL be submitted (run ${RUN_ID})"
else
    echo "# DRY RUN (nothing is submitted; add --execute to submit). Run ${RUN_ID}"
fi
echo "# run dir: ${RUN_DIR}"

PROGRAM="$(python3 "${PLANNER}" plan "${PLAN_PHASE}" --config "${CONFIG}" --run-dir "${RUN_DIR}" "${DRY_FLAG[@]}" "${PLAN_ARGS[@]}")" \
    || die "planning failed (see the message above)"
# The program consists only of note/submit/cancel/smoke/ensure_dir lines with every argument
# shell-quoted by Python's shlex.quote, so eval runs exactly those helper functions.
eval "${PROGRAM}"

if [[ "${EXECUTE}" -eq 0 ]]; then
    echo "# DRY RUN complete: ${N_SUBMITTED} sbatch command(s) printed, none submitted."
else
    echo "# Submitted ${N_SUBMITTED} job(s). Check with: bash scripts/submit_full_run.bash status --run-id ${RUN_ID} --execute"
fi
