#!/usr/bin/env python3
"""
pipeline_config_io.py
Shared loader for pipeline_config.yaml with storage-root environment variables.

pipeline_config.yaml refers to storage through ${MITO_*} placeholders so the same config works
for any user or for shared lab storage. This module:
  1. Fills in default MITO_* variables (mirroring scripts/mito_env.bash) when they are unset,
     so interactive runs and tests behave the same as SLURM jobs.
  2. Expands ${VAR} placeholders in the YAML text.
  3. Raises a clear error if any placeholder is still unexpanded, instead of letting a literal
     '${MITO_RAW_ROOT}/...' path cause a confusing 'file not found' much later.
"""

from __future__ import annotations

import getpass
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

# Order matters: later defaults reference earlier variables.
# Keep in sync with scripts/mito_env.bash (checked by test_pipeline_suite.py test_01l).
MITO_ENV_DEFAULTS: List[Tuple[str, str]] = [
    ("MITO_DATA_ROOT", "/sciclone/data10/{user}/mito_pipeline"),
    ("MITO_RAW_ROOT", "{MITO_DATA_ROOT}/raw_ometiff"),
    ("MITO_MODELS_DIR", "{MITO_DATA_ROOT}/models"),
    ("MITO_OUTPUT_ROOT", "{MITO_DATA_ROOT}/outputs"),
    ("MITO_TOOLS_DIR", "{MITO_DATA_ROOT}/tools"),
    ("MITO_SCRATCH_ROOT", "/sciclone/scr30/{user}/mito_work"),
    ("MITO_ILASTIK_BIN", "{MITO_TOOLS_DIR}/ilastik/run_ilastik.sh"),
]

_UNEXPANDED = re.compile(r"\$\{?[A-Za-z_][A-Za-z0-9_]*\}?")


def _current_user() -> str:
    return os.environ.get("USER") or getpass.getuser()


def resolve_mito_env(environ: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Returns the effective MITO_* values: existing environment first, documented defaults otherwise."""
    env = os.environ if environ is None else environ
    resolved: Dict[str, str] = {}
    for key, template in MITO_ENV_DEFAULTS:
        if env.get(key):
            resolved[key] = env[key]
        else:
            resolved[key] = template.format(user=_current_user(), **resolved)
    return resolved


def apply_mito_env_defaults() -> Dict[str, str]:
    """Exports default MITO_* values into os.environ for any that are unset; returns all values."""
    resolved = resolve_mito_env()
    for key, val in resolved.items():
        os.environ.setdefault(key, val)
    return resolved


def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """Loads pipeline_config.yaml, expanding ${MITO_*} and other environment placeholders."""
    apply_mito_env_defaults()
    with open(config_path, "r", encoding="utf-8") as f:
        raw_text = f.read()
    expanded = os.path.expandvars(raw_text)

    # Ignore comments when checking for leftovers so documentation lines can mention ${VAR}.
    leftovers = sorted({
        m.group(0)
        for line in expanded.splitlines()
        for m in _UNEXPANDED.finditer(line.split("#", 1)[0])
    })
    if leftovers:
        raise ValueError(
            f"Unexpanded environment variable(s) in {config_path}: {', '.join(leftovers)}. "
            "Export them (or `source scripts/mito_env.bash`) before running the pipeline."
        )
    cfg = yaml.safe_load(expanded) or {}
    return cfg
