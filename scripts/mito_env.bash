#!/bin/bash
# ==============================================================================
# mito_env.bash
# Single place that defines where the 4D Mitochondrial Pipeline keeps its data.
# ==============================================================================
# Source this file (do not execute it):
#   source scripts/mito_env.bash
#
# Every path derives from ONE root, MITO_DATA_ROOT. To move the whole pipeline
# (e.g. to shared Bestman lab storage), export a new root before sourcing:
#   export MITO_DATA_ROOT=/sciclone/data10/<lab_group>/mito_pipeline
# Any individual path below can still be overridden by exporting it first.
#
# The same defaults are mirrored in scripts/pipeline_config_io.py (Python side);
# test_pipeline_suite.py test_01l fails if the two ever drift apart.
# ==============================================================================

_MITO_USER="${USER:-$(whoami)}"

# Persistent project storage (backed up): raw OME-TIFFs, models, deliverables, tools
export MITO_DATA_ROOT="${MITO_DATA_ROOT:-/sciclone/data10/${_MITO_USER}/mito_pipeline}"
export MITO_RAW_ROOT="${MITO_RAW_ROOT:-${MITO_DATA_ROOT}/raw_ometiff}"
export MITO_MODELS_DIR="${MITO_MODELS_DIR:-${MITO_DATA_ROOT}/models}"
export MITO_OUTPUT_ROOT="${MITO_OUTPUT_ROOT:-${MITO_DATA_ROOT}/outputs}"
export MITO_TOOLS_DIR="${MITO_TOOLS_DIR:-${MITO_DATA_ROOT}/tools}"

# High-speed scratch (NOT backed up, may be purged): per-movie intermediates
export MITO_SCRATCH_ROOT="${MITO_SCRATCH_ROOT:-/sciclone/scr30/${_MITO_USER}/mito_work}"

# Tool locations inside MITO_TOOLS_DIR (written by setup_hpc_env.bash)
export MITO_ILASTIK_BIN="${MITO_ILASTIK_BIN:-${MITO_TOOLS_DIR}/ilastik/run_ilastik.sh}"

unset _MITO_USER
