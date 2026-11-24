#!/usr/bin/env python3
"""
plan_full_run.py
Planner behind scripts/submit_full_run.bash: decides WHAT to submit for each phase of a full-archive
run and prints it as a short bash program of 'note' / 'submit' / 'cancel' / 'smoke' lines, which
submit_full_run.bash either prints (dry run, the default) or executes (--execute).

Standard library only, so it runs on dev (module load python/3.12.7), inside jobs, or anywhere.

What it owns
  * Size classes and the memory/time model: every movie gets a class (S/M/L by
    estimated_uncompressed_gb) and a predicted peak memory / run time per stage; each class becomes
    its own job array whose request = the largest prediction in the class x safety factors.
  * Canary selection (~10 movies spanning the classes) and the refit of the model from the
    canary's sacct MaxRSS / Elapsed.
  * Throttles: 20% of a partition's cores (queried live with sinfo; config fallbacks), I/O throttles
    for scan/extract, lists chunked at <= 1000 tasks per array.
  * Dependencies: merge afterany:scan, extract afterok:merge, stage45 aftercorr:stage3 (same list,
    same index), aggregate afterany:stage45. Dependencies on jobs that already left the queue are
    dropped (Slurm rejects unknown job ids).
  * Stage-done hashes for run_pipeline_slurm.bash (.done markers), the run manifest, retries,
    status, and the opt-in overflow check (deny list, live PreemptMode, per-cluster smoke test).

Commands (all take --config, --run-dir; plan commands also --dry-run):
  plan {scan,merge,extract,canary,stage3,stage45,aggregate,retry,overflow}   print the submission program
  extract-list   write <run>/extract/extract_files.txt from the merged inventory (merge task)
  collector-args print the collector filter flags (--min-z/-t/-c, --compression) from the config
  record         store a submitted job id in the manifest
  smoke-record   store a per-cluster environment smoke-test result
  status         one sacct summary of every job in the manifest
  stage-hash     print the Stage 3/4/5 hashes used by the .done markers
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
DEFAULT_CONFIG = REPO_DIR / "pipeline_config.yaml"

# Slurm states that a retry should resubmit (CANCELLED is handled separately for Stage 4-5).
RETRY_STATES = {"FAILED", "OUT_OF_MEMORY", "TIMEOUT", "NODE_FAIL", "BOOT_FAIL", "DEADLINE"}
QC_FAIL_EXIT = 3  # run_pipeline_slurm.bash: Stage 3 FAIL_NO_SIGNAL -> not retried
FG_FAIL_EXIT = 5  # run_hybrid_segmentation.py: implausible foreground (bad ilastik model) -> not retried
PIPELINE_PHASES = ("canary", "stage3", "stage45")

# ================================================================================================
# 1. Minimal YAML reader for the hpc_run block (stdlib only)
# ================================================================================================


def strip_yaml_comment(line: str) -> str:
    """Removes a '#' comment unless the '#' sits inside single or double quotes."""
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
        elif ch == "#" and (i == 0 or line[i - 1] in " \t"):
            return line[:i].rstrip()
    return line.rstrip()


def parse_yaml_scalar(text: str) -> Any:
    """Parses one YAML scalar or inline [list]: quoted strings, true/false, null, int, float, bare text."""
    text = text.strip()
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        return [parse_yaml_scalar(item) for item in inner.split(",")] if inner else []
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        return text[1:-1]
    low = text.lower()
    if low in ("true", "yes", "on"):
        return True
    if low in ("false", "no", "off"):
        return False
    if low in ("null", "~", ""):
        return None
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def parse_simple_yaml(lines: Iterable[str]) -> Dict[str, Any]:
    """
    Parses indentation-nested 'key: value' mappings with inline scalars/lists (the subset used by
    the hpc_run block). Raises ValueError on anything else (e.g. '- item' block lists).
    """
    root: Dict[str, Any] = {}
    stack: List[Tuple[int, Dict[str, Any]]] = [(-1, root)]
    for lineno, raw in enumerate(lines, 1):
        line = strip_yaml_comment(raw.rstrip("\r\n"))
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        content = line.strip()
        if content.startswith("- "):
            raise ValueError(f"line {lineno}: block lists are not supported here; use inline [a, b]")
        key, sep, value = content.partition(":")
        if not sep:
            raise ValueError(f"line {lineno}: expected 'key: value', got '{content}'")
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        key = key.strip().strip("'\"")
        if value.strip() == "":
            child: Dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = parse_yaml_scalar(value)
    return root


def top_level_blocks(text: str) -> Dict[str, str]:
    """Splits YAML text into {top-level key: its full text block (key line + indented body)}."""
    blocks: Dict[str, List[str]] = {}
    current = None
    for line in text.splitlines():
        if line and not line[0].isspace() and not line.startswith("#"):
            key = line.split(":", 1)[0].strip()
            current = key
            blocks[current] = [line]
        elif current is not None:
            blocks[current].append(line)
    return {k: "\n".join(v) for k, v in blocks.items()}


def load_hpc_run(config_path: Path) -> Dict[str, Any]:
    """Returns the hpc_run block of pipeline_config.yaml as nested dicts (DEFAULT_HPC_RUN merged in)."""
    text = Path(config_path).read_text(encoding="utf-8")
    block = top_level_blocks(text).get("hpc_run", "")
    parsed = parse_simple_yaml(block.splitlines()).get("hpc_run", {}) if block else {}
    return deep_merge(DEFAULT_HPC_RUN, parsed or {})


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively overlays override onto a copy of base."""
    out = json.loads(json.dumps(base))
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


# Built-in defaults (identical to the shipped pipeline_config.yaml hpc_run block), so a config
# without the block still plans sensibly.
DEFAULT_HPC_RUN: Dict[str, Any] = {
    "primary_cluster": "bora",
    "primary_partition": "batch",
    "large_partition": "hima",
    "large_exclude_gpu_nodes": True,
    "large_gpu_nodes_fallback": ["hi04", "hi05", "hi07"],
    "footprint_fraction": 0.20,
    "partition_cores_fallback": {"batch": 900, "hima": 128},
    "partition_max_mem_gb": {"batch": 120, "hima": 1000},
    "scan_throttle": 20,
    "extract_throttle": 20,
    "max_array_size": 1000,
    "run_root": "",
    "log_root": "",
    "job_name_prefix": "mito",
    "collector": {"min_z": 15, "min_t": 10, "min_c": 2, "compression": "LZW",
                  "extract_labels": ["MITO"], "files_per_scan_shard": 5},
    "scan": {"cpus": 1, "mem_gb": 4, "time": "00:30:00"},
    "merge": {"cpus": 1, "mem_gb": 2, "time": "00:15:00"},
    "extract": {"cpus": 2, "mem_gb": 6, "sec_per_gb": 115, "sec_base": 300, "time_unknown": "04:00:00"},
    "aggregate": {"cpus": 1, "mem_gb": 4, "time": "00:30:00"},
    "size_classes": {"names": ["S", "M", "L"], "edges_gb": [2.0, 8.0], "large_class": "L"},
    "model": {
        "stage3_mem_gb_per_gb": 2.75, "stage3_mem_gb_base": 0.5, "stage3_sec_per_gb": 65, "stage3_cpus": 1,
        "stage45_cpus": 4, "stage45_mem_gb_per_gb": 2.75, "stage45_mem_gb_base": 1.5,
        "stage45_mem_per_worker_frame": 8.0, "ilastik_sec_per_channel_gb": 290, "ilastik_thread_saturation": 4,
        "watershed_sec_per_channel_gb": 74, "stage5_sec": 60, "default_size_t": 10, "default_size_c": 2,
        "time_factor": 2.0, "min_time_min": 10, "mem_factor": 1.5, "min_mem_gb": 2,
        "refit_min_scale": 0.5, "refit_max_scale": 4.0,
    },
    "canary_size": 10,
    "retry_factor": 1.5,
    "max_retries": 2,
    "overflow": {"clusters": [], "partitions": {}, "deny_clusters": [], "require_preempt_mode": "OFF",
                 "pending_threshold_min": 30},
}

# ================================================================================================
# 2. Stage-done hashes (.done markers written by run_pipeline_slurm.bash)
# ================================================================================================

STAGE_SECTIONS = {
    3: ["pipeline_version", "metadata_defaults", "stage3_preprocessing"],
    4: ["stage4_hybrid_segmentation"],
    5: ["stage5_tracking_and_mitolandscape"],
}
STAGE_SCRIPTS = {
    3: ["preprocess_4d_stack.py", "tiff_io.py", "pipeline_config_io.py"],
    4: ["run_hybrid_segmentation.py", "tiff_io.py", "pipeline_config_io.py"],
    5: ["run_trackmate_sparse_lap.py", "run_arbor_and_mitolandscape.py", "tiff_io.py", "pipeline_config_io.py"],
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def normalized_section_text(block: str) -> str:
    """A config block with comments and trailing spaces removed, so comment edits do not force re-runs."""
    lines = [strip_yaml_comment(line) for line in block.splitlines()]
    return "\n".join(line for line in lines if line.strip())


def git_head(repo_dir: Path) -> str:
    """Current git HEAD SHA of the repo, or 'nogit'."""
    out = run_cmd(["git", "-C", str(repo_dir), "rev-parse", "HEAD"])
    return out.strip() if out else "nogit"


def file_fingerprint(path: Optional[Path]) -> str:
    """Cheap identity of a (possibly multi-GB) input file: absolute path + size + mtime."""
    if not path:
        return "none"
    try:
        st = Path(path).stat()
    except OSError:
        return f"missing:{path}"
    return f"{Path(path).resolve()}|{st.st_size}|{int(st.st_mtime)}"


def ilp_fingerprint(path: Optional[Path]) -> str:
    """SHA-256 of a valid ilastik project (HDF5); retraining the .ilp must re-run Stage 4."""
    if not path:
        return "no-ilp"
    p = Path(path)
    try:
        with open(p, "rb") as f:
            if f.read(8) != b"\x89HDF\r\n\x1a\n":
                return "no-ilp"
    except OSError:
        return "no-ilp"
    return sha256_file(p)


def compute_stage_hashes(
    config_path: Path,
    repo_dir: Path,
    input_path: Optional[Path],
    ilp_path: Optional[Path] = None,
    head: Optional[str] = None,
) -> Dict[int, str]:
    """
    Stage hashes chained downstream: h3 = H(config sections 3, stage-3 scripts, input file);
    h4 = H(h3, sections 4, scripts 4, .ilp); h5 = H(h4, sections 5, scripts 5). Changing anything a
    stage depends on changes its hash AND every later stage's hash, so later stages re-run too.
    The git HEAD is NOT included by default: the sha256 of every stage script already captures the
    code, and including HEAD would re-run every finished movie after any commit (even docs-only).
    Pass head=... explicitly to pin a hash to a commit. The run manifest still records git_sha.
    """
    text = Path(config_path).read_text(encoding="utf-8")
    blocks = top_level_blocks(text)
    head = head if head is not None else "not-hashed"
    scripts_dir = Path(repo_dir) / "scripts"
    hashes: Dict[int, str] = {}
    previous = ""
    for stage in (3, 4, 5):
        h = hashlib.sha256()
        h.update(f"stage={stage}\nupstream={previous}\nhead={head}\n".encode())
        for sec in STAGE_SECTIONS[stage]:
            h.update(f"[{sec}]\n{normalized_section_text(blocks.get(sec, ''))}\n".encode())
        for name in STAGE_SCRIPTS[stage]:
            script = scripts_dir / name
            digest = sha256_file(script) if script.exists() else "missing"
            h.update(f"script:{name}={digest}\n".encode())
        if stage == 3:
            h.update(f"input={file_fingerprint(input_path)}\n".encode())
        if stage == 4:
            h.update(f"ilp={ilp_fingerprint(ilp_path)}\n".encode())
        hashes[stage] = h.hexdigest()
        previous = hashes[stage]
    return hashes


def marker_matches(marker_path: Path, expected_hash: str) -> bool:
    """True when a .done marker exists and its first line equals the expected stage hash."""
    try:
        with open(marker_path, encoding="utf-8") as f:
            return f.readline().strip() == expected_hash
    except OSError:
        return False


def write_marker(marker_path: Path, stage_hash: str, stage: int, extra: Optional[Dict[str, str]] = None) -> None:
    """Writes a .done marker: the hash on line 1, then key=value metadata lines (atomic rename)."""
    lines = [stage_hash, f"stage={stage}",
             f"completed={datetime.datetime.now(datetime.timezone.utc).isoformat()}"]
    for k, v in (extra or {}).items():
        lines.append(f"{k}={v}")
    tmp = marker_path.with_name(marker_path.name + ".tmp")
    tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(tmp, marker_path)


# ================================================================================================
# 3. Movies, size classes, and the memory/time model
# ================================================================================================


@dataclass
class Movie:
    entry: str                 # line as it appears in a movie list (absolute path preferred)
    path: str                  # absolute path to the OME-TIFF
    gb: float                  # estimated uncompressed size (GB)
    size_t: int
    size_c: int
    size_source: str           # "inventory" or "file_size"
    cohort: str = ""
    mito_label: str = ""
    sld_path: str = ""
    size_class: str = ""
    extra: Dict[str, str] = field(default_factory=dict)

    @property
    def movie_id(self) -> str:
        name = Path(self.path).name
        for suffix in (".ome.tif", ".ome.tiff", ".tif", ".tiff"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return name


def size_class(gb: float, edges: Sequence[float], names: Sequence[str]) -> str:
    """Class by uncompressed GB: names[i] for edges[i-1] <= gb < edges[i] (len(names) = len(edges)+1)."""
    if len(names) != len(edges) + 1:
        raise ValueError("size_classes: need exactly one more name than edges")
    for i, edge in enumerate(edges):
        if gb < float(edge):
            return names[i]
    return names[-1]


def predict_stage3(m: Movie, model: Dict[str, Any]) -> Tuple[float, float]:
    """Predicted (peak memory GB, seconds) for Stage 3 on one core."""
    mem = float(model["stage3_mem_gb_per_gb"]) * m.gb + float(model["stage3_mem_gb_base"])
    sec = float(model["stage3_sec_per_gb"]) * m.gb
    return mem, sec


def predict_stage45(m: Movie, model: Dict[str, Any], cpus: Optional[int] = None) -> Tuple[float, float]:
    """
    Predicted (peak memory GB, seconds) for Stages 4-5 with `cpus` cores:
      ilastik   = k_il * channel_GB * saturation / min(cpus, saturation)   (threads help up to 4)
      watershed = k_ws * channel_GB * ceil(T / workers) / T                 (frames in parallel)
      Stage 5   = constant
      memory    = a * GB + b + workers * c * frame_GB  (each watershed worker holds a few frame copies)
    channel_GB = movie GB / SizeC (ilastik and the watershed see only the mito channel).
    """
    cpus = int(cpus or model["stage45_cpus"])
    t = max(1, int(m.size_t or model["default_size_t"]))
    c = max(1, int(m.size_c or model["default_size_c"]))
    chan_gb = m.gb / c
    frame_gb = chan_gb / t
    sat = max(1, int(model["ilastik_thread_saturation"]))
    workers = max(1, min(cpus, t))
    ilastik = float(model["ilastik_sec_per_channel_gb"]) * chan_gb * sat / min(cpus, sat)
    watershed = float(model["watershed_sec_per_channel_gb"]) * chan_gb * math.ceil(t / workers) / t
    sec = ilastik + watershed + float(model["stage5_sec"])
    mem = (float(model["stage45_mem_gb_per_gb"]) * m.gb + float(model["stage45_mem_gb_base"])
           + workers * float(model["stage45_mem_per_worker_frame"]) * frame_gb)
    return mem, sec


def to_request(mem_gb: float, sec: float, model: Dict[str, Any],
               mem_scale: float = 1.0, time_scale: float = 1.0) -> Tuple[int, int]:
    """Turns a prediction into a Slurm request: (memory GB, minutes) with factors and floors."""
    mem_req = max(float(model["min_mem_gb"]), mem_gb * mem_scale * float(model["mem_factor"]))
    min_req = max(float(model["min_time_min"]), sec * time_scale * float(model["time_factor"]) / 60.0)
    return int(math.ceil(mem_req)), int(math.ceil(min_req))


def fmt_time(minutes: float) -> str:
    """Minutes -> Slurm D-HH:MM:SS / HH:MM:SS."""
    total = int(math.ceil(minutes)) * 60
    days, rem = divmod(total, 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return f"{days}-{h:02d}:{m:02d}:{s:02d}" if days else f"{h:02d}:{m:02d}:{s:02d}"


def parse_time_to_min(text: str) -> float:
    """Slurm time (D-HH:MM:SS, HH:MM:SS, MM:SS) -> minutes."""
    text = str(text).strip()
    days = 0
    if "-" in text:
        d, text = text.split("-", 1)
        days = int(d)
    parts = [int(float(p)) for p in text.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    h, m, s = parts[-3:]
    return days * 1440 + h * 60 + m + s / 60.0


def chunk(items: Sequence[Any], size: int) -> List[List[Any]]:
    """Splits items into consecutive chunks of at most `size` (Slurm MaxArraySize)."""
    size = max(1, int(size))
    return [list(items[i: i + size]) for i in range(0, len(items), size)]


def split_throttle(cap_tasks: int, sizes: Sequence[int]) -> List[int]:
    """
    Shares a concurrency cap between arrays in proportion to their task counts (at least 1 each,
    never more than the array has tasks).
    """
    total = sum(sizes)
    if total == 0:
        return [0 for _ in sizes]
    out = []
    for n in sizes:
        share = int(math.floor(cap_tasks * n / total))
        out.append(max(1, min(n, share)))
    return out


def select_canary(movies: Sequence[Movie], n: int) -> List[Movie]:
    """
    Picks ~n movies spanning the size classes: every non-empty class gets at least one slot
    (proportional otherwise); inside a class, movies are evenly spaced over the size-sorted list and
    the LARGEST is always included (it decides the class's memory request).
    """
    by_class: Dict[str, List[Movie]] = {}
    for m in movies:
        by_class.setdefault(m.size_class, []).append(m)
    if not by_class:
        return []
    total = len(movies)
    n = max(n, len(by_class))
    alloc = {c: max(1, int(round(n * len(ms) / total))) for c, ms in by_class.items()}
    # Trim the largest allocations until the total fits n (each class keeps >= 1).
    while sum(alloc.values()) > n:
        c = max(alloc, key=lambda k: (alloc[k], k))
        if alloc[c] <= 1:
            break
        alloc[c] -= 1
    picked: List[Movie] = []
    for c in sorted(by_class):
        ms = sorted(by_class[c], key=lambda m: (m.gb, m.path))
        k = min(alloc[c], len(ms))
        if k == 1:
            idx = [len(ms) - 1]
        else:
            idx = sorted({int(round(i * (len(ms) - 1) / (k - 1))) for i in range(k)})
        picked.extend(ms[i] for i in idx)
    return picked


# ================================================================================================
# 4. Inventory / movie lists
# ================================================================================================


def read_inventory(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_list(path: Path) -> List[str]:
    """Movie/file list: one entry per line; blank lines and '#' comments ignored; CRLF tolerated."""
    out = []
    with open(path, encoding="utf-8-sig") as f:
        for raw in f:
            line = raw.strip()
            if line and not line.startswith("#"):
                out.append(line)
    return out


def write_list(path: Path, entries: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for e in entries:
            f.write(f"{e}\n")


def label_allowed(label: str, allowed: Sequence[str]) -> bool:
    """Empty label = filter was off -> allowed; otherwise must be in the allowed set."""
    return (not label) or label in set(allowed)


def build_movies(
    entries: Sequence[str],
    inventory: Sequence[Dict[str, str]],
    raw_root: str,
    model: Dict[str, Any],
    classes: Dict[str, Any],
    allowed_labels: Optional[Sequence[str]] = None,
) -> Tuple[List[Movie], List[str], List[str]]:
    """
    Turns movie-list entries into Movie records with sizes from the inventory (matched by
    extracted_ometiff_relpath) or, failing that, from the file size / 0.55 (LZW ratio).
    Returns (movies, missing_files, label_excluded).
    """
    by_rel = {r.get("extracted_ometiff_relpath", ""): r for r in inventory if r.get("extracted_ometiff_relpath")}
    movies, missing, excluded = [], [], []
    for entry in entries:
        path = entry if os.path.isabs(entry) else os.path.join(raw_root or "", entry)
        row = by_rel.get(entry)
        if row is None:
            for rel, r in by_rel.items():
                if path.endswith("/" + rel):
                    row = r
                    break
        if row is not None and allowed_labels is not None and not label_allowed(row.get("mito_label", ""), allowed_labels):
            excluded.append(entry)
            continue
        if row is not None:
            gb = float(row.get("estimated_uncompressed_gb") or 0.0)
            m = Movie(entry=path, path=path, gb=gb, size_t=int(row.get("size_t") or 0),
                      size_c=int(row.get("size_c") or 0), size_source="inventory",
                      cohort=row.get("cohort_year", ""), mito_label=row.get("mito_label", ""),
                      sld_path=row.get("sld_path", ""))
        else:
            try:
                gb = os.path.getsize(path) / (1024 ** 3) / 0.55
            except OSError:
                gb = 0.0
            m = Movie(entry=path, path=path, gb=round(gb, 4), size_t=0, size_c=0, size_source="file_size",
                      cohort=Path(path).parent.name)
        if not os.path.exists(path):
            missing.append(path)
            continue
        m.size_class = size_class(m.gb, classes["edges_gb"], classes["names"])
        movies.append(m)
    return movies, missing, excluded


# ================================================================================================
# 5. Slurm queries (live, with fallbacks)
# ================================================================================================


def run_cmd(argv: Sequence[str], timeout: int = 60, env: Optional[Dict[str, str]] = None) -> Optional[str]:
    """Runs a command and returns stdout, or None if it is missing, fails, or times out."""
    try:
        proc = subprocess.run(list(argv), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, timeout=timeout, check=False, env=env)
    except (OSError, subprocess.SubprocessError):
        return None
    return proc.stdout if proc.returncode == 0 else None


def run_on_cluster(cluster: str, argv: Sequence[str], primary: str) -> Optional[str]:
    """
    Runs a Slurm client command against another SciClone cluster by loading its module in a
    throwaway bash (module files set SLURM_CONF). The primary cluster uses the current environment.
    """
    if cluster == primary:
        return run_cmd(argv)
    script = ('type module >/dev/null 2>&1 || source /usr/share/Modules/init/bash 2>/dev/null; '
              'module unload slurm >/dev/null 2>&1; '   # slurm modules conflict with each other
              'module load "slurm/$1" >/dev/null 2>&1 || exit 97; shift; exec "$@"')
    return run_cmd(["bash", "-c", script, "_", cluster, *argv])


def partition_info(partition: str, exclude_gpu: bool, cluster: str, primary: str,
                   hpc: Dict[str, Any]) -> Dict[str, Any]:
    """
    Live partition size: total cores, largest node memory, and GPU nodes (to exclude) from
    'sinfo -N'. Falls back to the config values when sinfo is unavailable.
    """
    out = run_on_cluster(cluster, ["sinfo", "-h", "-N", "-p", partition, "-o", "%N|%c|%m|%f|%G"], primary)
    nodes: Dict[str, Tuple[int, int, bool]] = {}
    if out:
        for line in out.splitlines():
            parts = line.strip().split("|")
            if len(parts) < 5:
                continue
            name, cpus, mem, feats, gres = parts[:5]
            try:
                is_gpu = "gpu" in feats.lower().split(",") or gres.lower().startswith("gpu")
                nodes[name] = (int(cpus), int(mem), is_gpu)
            except ValueError:
                continue
    if nodes:
        gpu_nodes = sorted(n for n, (_, _, g) in nodes.items() if g)
        usable = {n: v for n, v in nodes.items() if not (exclude_gpu and v[2])}
        return {
            "source": "sinfo",
            "cores": sum(v[0] for v in usable.values()),
            "max_mem_gb": max(v[1] for v in usable.values()) / 1024.0 if usable else 0.0,
            "exclude": gpu_nodes if exclude_gpu else [],
            "nodes": len(usable),
        }
    fallback_cores = int((hpc.get("partition_cores_fallback") or {}).get(partition, 0) or 0)
    return {
        "source": "config_fallback",
        "cores": fallback_cores,
        "max_mem_gb": float((hpc.get("partition_max_mem_gb") or {}).get(partition, 0) or 0),
        "exclude": list(hpc.get("large_gpu_nodes_fallback", [])) if exclude_gpu else [],
        "nodes": 0,
    }


def cap_tasks(cores: int, fraction: float, cpus_per_task: int) -> int:
    """Concurrent tasks allowed by the footprint cap: floor(cores * fraction / cpus), at least 1."""
    return max(1, int(math.floor(cores * float(fraction) / max(1, cpus_per_task))))


def active_job_ids(job_ids: Sequence[str]) -> List[str]:
    """
    Job ids still known to the scheduler (pending or running). Placeholder ids from a dry run
    (DRY_*) count as active. Unknown/finished ids are dropped: Slurm rejects dependencies on them.
    """
    real = [j for j in job_ids if j and not j.startswith("DRY_")]
    active = [j for j in job_ids if j and j.startswith("DRY_")]
    if real:
        out = run_cmd(["squeue", "-h", "-j", ",".join(real), "-o", "%F"])
        if out is not None:
            alive = {line.strip() for line in out.splitlines() if line.strip()}
            active.extend(j for j in real if j in alive)
    return active


def parse_sacct(text: str) -> Dict[Tuple[str, int], Dict[str, Any]]:
    """
    Parses 'sacct -P -n --units=M -o JobID,State,ExitCode,Elapsed,MaxRSS' into
    {(array_job_id, task_index): {state, exit, elapsed_s, maxrss_mb}}; MaxRSS comes from the
    .batch step, state/exit/elapsed from the task line.
    """
    tasks: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) < 4:
            continue
        jobid, state, exitcode, elapsed = parts[:4]
        rss = parts[4] if len(parts) > 4 else ""
        m = re.fullmatch(r"(\d+)_(\d+)(?:\.(\w+))?", jobid)
        if not m:
            continue
        key = (m.group(1), int(m.group(2)))
        rec = tasks.setdefault(key, {"state": "", "exit": 0, "elapsed_s": 0.0, "maxrss_mb": 0.0})
        step = m.group(3)
        if step is None:
            rec["state"] = state.split()[0] if state else ""
            try:
                rec["exit"] = int(exitcode.split(":")[0])
            except ValueError:
                rec["exit"] = 0
            rec["elapsed_s"] = parse_time_to_min(elapsed) * 60.0 if elapsed else 0.0
        mr = re.fullmatch(r"([\d.]+)([KMGT]?)", rss.strip())
        if mr:
            scale = {"K": 1 / 1024.0, "M": 1.0, "G": 1024.0, "T": 1024.0 ** 2, "": 1 / (1024.0 ** 2)}[mr.group(2)]
            rec["maxrss_mb"] = max(rec["maxrss_mb"], float(mr.group(1)) * scale)
    return tasks


# ================================================================================================
# 6. Manifest
# ================================================================================================


class Manifest:
    """run_manifest.json: provenance + every array submitted (ids, lists, resources, retries)."""

    def __init__(self, path: Path):
        self.path = path
        self.data: Dict[str, Any] = {}
        if path.exists():
            self.data = json.loads(path.read_text(encoding="utf-8"))
        self.data.setdefault("jobs", {})
        self.data.setdefault("phases", {})

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(self.path.name + ".tmp")
        tmp.write_text(json.dumps(self.data, indent=2, sort_keys=True), encoding="utf-8")
        os.replace(tmp, self.path)

    def jobs(self, phase: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        return {n: j for n, j in self.data["jobs"].items() if phase is None or j.get("phase") == phase}

    def add_job(self, name: str, **info: Any) -> None:
        info.setdefault("job_id", "")
        info.setdefault("retries", [])
        self.data["jobs"][name] = info


def run_dirs(run_dir: Path, dry_run: bool) -> Tuple[Path, Path]:
    """(state dir holding manifest + generated lists, manifest path). Dry runs use <run>/dryrun."""
    state = run_dir / "dryrun" if dry_run else run_dir
    return state, state / "run_manifest.json"


def init_manifest_provenance(man: Manifest, run_dir: Path, config: Path, repo: Path, dry_run: bool) -> None:
    """Fills/refreshes run-level provenance (git SHA, config hashes)."""
    text = config.read_text(encoding="utf-8")
    d = man.data
    d.setdefault("created", datetime.datetime.now(datetime.timezone.utc).isoformat())
    d["run_id"] = run_dir.name
    d["run_dir"] = str(run_dir)
    d["dry_run"] = dry_run
    d["repo"] = str(repo)
    d["git_sha"] = git_head(repo)
    d["config_path"] = str(config)
    d["config_sha256"] = sha256_bytes(text.encode())
    d["hpc_run_sha256"] = sha256_bytes(top_level_blocks(text).get("hpc_run", "").encode())
    d["updated"] = datetime.datetime.now(datetime.timezone.utc).isoformat()


# ================================================================================================
# 7. Program emission (consumed by submit_full_run.bash)
# ================================================================================================


class Program:
    """Collects the bash lines that submit_full_run.bash evaluates."""

    def __init__(self) -> None:
        self.lines: List[str] = []

    def note(self, text: str = "") -> None:
        self.lines.append("note " + shlex.quote(text))

    def submit(self, name: str, sbatch_args: Sequence[str], cluster: Optional[str] = None) -> None:
        head = ["submit"] + (["--cluster", cluster] if cluster else []) + [name]
        self.lines.append(" ".join(shlex.quote(a) for a in head + list(sbatch_args)))

    def cancel(self, job_id: str, cluster: Optional[str] = None) -> None:
        self.lines.append(" ".join(shlex.quote(a) for a in ["cancel", cluster or "", job_id]))

    def smoke(self, cluster: str, partition: str) -> None:
        self.lines.append(" ".join(shlex.quote(a) for a in ["smoke", cluster, partition]))

    def mkdirs(self, *dirs: Path) -> None:
        self.lines.append(" ".join(shlex.quote(a) for a in ["ensure_dir", *map(str, dirs)]))

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"


@dataclass
class Ctx:
    """Everything a phase planner needs."""
    config: Path
    repo: Path
    run_dir: Path
    log_dir: Path
    state_dir: Path
    man: Manifest
    hpc: Dict[str, Any]
    dry_run: bool
    raw_root: str
    submit_script: Path

    @property
    def prefix(self) -> str:
        return f"{self.hpc.get('job_name_prefix', 'mito')}_{self.run_dir.name}"

    def job_name(self, *parts: str) -> str:
        return "_".join([self.prefix, *parts])

    def common_args(self, name: str, partition: str, cpus: int, mem_gb: int, time_min: float,
                    exclude: Sequence[str] = (), array: Optional[str] = None) -> List[str]:
        """sbatch options shared by every job: name, partition, resources, per-run log files."""
        args = [f"--job-name={name}", f"--partition={partition}", "--nodes=1", "--ntasks=1",
                f"--cpus-per-task={cpus}", f"--mem={mem_gb}G", f"--time={fmt_time(time_min)}"]
        if array:
            args.append(f"--array={array}")
            pattern = "%x_%A_%a"
        else:
            pattern = "%x_%j"
        args += [f"--output={self.log_dir}/{pattern}.out", f"--error={self.log_dir}/{pattern}.err"]
        if exclude:
            args.append(f"--exclude={','.join(exclude)}")
        args.append(f"--export=ALL,MITO_REPO_DIR={self.repo},MITO_RUN_DIR={self.run_dir}")
        return args


def dependency_arg(kind: str, ids: Sequence[str], prog: Program) -> List[str]:
    """
    '--dependency=kind:id1:id2' for the ids still in the queue. 'ids' may contain '@@name@@'
    placeholders for jobs submitted earlier in the SAME program (always kept).
    """
    placeholders = [i for i in ids if i.startswith("@@")]
    live = active_job_ids([i for i in ids if not i.startswith("@@")])
    dropped = sorted(set(i for i in ids if not i.startswith("@@")) - set(live))
    if dropped:
        prog.note(f"dependency {kind} on {','.join(dropped)} dropped: no longer in the queue (finished)")
    keep = placeholders + live
    return [f"--dependency={kind}:{':'.join(keep)}"] if keep else []


def ph(name: str) -> str:
    """Placeholder for the id of a job submitted earlier in the same program."""
    return f"@@{name}@@"


# ---------------------------------------------------------------- scan / merge / extract


def plan_scan(ctx: Ctx, prog: Program, source_dir: Path, file_list: Optional[Path],
              plasmid_csv: Optional[Path], chain_extract: bool, overrides: Optional[Path] = None) -> None:
    """Phase A (sharded scan) + A' (merge, afterany) [+ B (extract, afterok) when chained]."""
    col = ctx.hpc["collector"]
    if file_list:
        files = sorted(set(read_list(file_list)))
    else:
        files = sorted(str(p) for p in source_dir.rglob("*") if p.is_file() and p.suffix.lower() in (".sld", ".sldy"))
    scan_dir = ctx.state_dir / "scan"
    scan_list = scan_dir / "scan_files.txt"
    write_list(scan_list, files)
    n_shards = max(1, min(int(ctx.hpc["max_array_size"]),
                          math.ceil(len(files) / max(1, int(col["files_per_scan_shard"])))))
    sc = ctx.hpc["scan"]
    throttle = min(int(ctx.hpc["scan_throttle"]), n_shards)
    shard_dir = ctx.run_dir / "scan" / "shards"
    merged_dir = ctx.run_dir / "scan" / "merged"
    ctx.man.data["phases"]["scan"] = {
        "source_dir": str(source_dir), "file_list": str(scan_list), "n_files": len(files),
        "n_shards": n_shards, "shard_dir": str(shard_dir), "merged_dir": str(merged_dir),
        "plasmid_csv": str(plasmid_csv) if plasmid_csv else None, "raw_root": ctx.raw_root,
        "mito_overrides": str(overrides) if overrides else None,
    }
    prog.note(f"Phase A  scan: {len(files)} .sld/.sldy file(s) -> {n_shards} shard(s) "
              f"(~{col['files_per_scan_shard']} files each), at most {throttle} at once (NFS read throttle)")
    if not files:
        prog.note("No .sld/.sldy files found: nothing to submit.")
        return
    prog.mkdirs(ctx.log_dir, shard_dir, merged_dir)
    name = ctx.job_name("scan")
    task = [str(ctx.submit_script), "__scan", "--run-dir", str(ctx.run_dir), "--source-dir", str(source_dir),
            "--file-list", str(scan_list), "--shards", str(n_shards), "--raw-root", ctx.raw_root,
            "--config", str(ctx.config)]
    if plasmid_csv:
        task += ["--plasmid-csv", str(plasmid_csv)]
    if overrides:
        task += ["--mito-overrides", str(overrides)]
    args = ctx.common_args(name, ctx.hpc["primary_partition"], int(sc["cpus"]), int(sc["mem_gb"]),
                           parse_time_to_min(sc["time"]), array=f"1-{n_shards}%{throttle}")
    prog.submit(name, args + task)
    ctx.man.add_job(name, phase="scan", n_tasks=n_shards, list=str(scan_list), cpus=int(sc["cpus"]),
                    mem_gb=int(sc["mem_gb"]), time_min=parse_time_to_min(sc["time"]),
                    partition=ctx.hpc["primary_partition"], array=f"1-{n_shards}%{throttle}",
                    cluster=ctx.hpc["primary_cluster"], task=task)

    mname = ctx.job_name("merge")
    mg = ctx.hpc["merge"]
    prog.note("Phase A' merge shard reports -> one inventory + movie_list + extract_files.txt "
              "(afterany: runs even if a shard failed, and then reports the missing shards)")
    margs = ctx.common_args(mname, ctx.hpc["primary_partition"], int(mg["cpus"]), int(mg["mem_gb"]),
                            parse_time_to_min(mg["time"]))
    mtask = [str(ctx.submit_script), "__merge", "--run-dir", str(ctx.run_dir), "--phase", "scan",
             "--shards", str(n_shards), "--raw-root", ctx.raw_root, "--config", str(ctx.config),
             "--labels", ",".join(col["extract_labels"])]
    prog.submit(mname, margs + [f"--dependency=afterany:{ph(name)}"] + mtask)
    ctx.man.add_job(mname, phase="merge", n_tasks=1, cluster=ctx.hpc["primary_cluster"], task=mtask)

    if chain_extract:
        plan_extract(ctx, prog, source_dir, after_merge=mname, upper_bound=len(files))


def plan_merge(ctx: Ctx, prog: Program) -> None:
    """
    Re-submits the scan merge alone (e.g. after 'retry --phase scan' fixed a failed shard); waits
    (afterany) for any scan job or retry still in the queue.
    """
    scan = ctx.man.data["phases"].get("scan")
    if not scan:
        raise SystemExit("[ERROR] No scan recorded in the manifest: run 'scan' first.")
    ids = []
    for j in ctx.man.jobs("scan").values():
        ids += [j.get("job_id", "")] + [r.get("job_id", "") for r in j.get("retries", [])]
    attempt = sum(1 for n in ctx.man.jobs("merge") if n.startswith(ctx.job_name("merge_rerun")))
    name = ctx.job_name("merge_rerun", str(attempt + 1))
    mg = ctx.hpc["merge"]
    prog.note(f"Phase A' merge (re-run): {scan['n_shards']} shard(s) -> {scan['merged_dir']}")
    prog.mkdirs(ctx.log_dir)
    args = ctx.common_args(name, ctx.hpc["primary_partition"], int(mg["cpus"]), int(mg["mem_gb"]),
                           parse_time_to_min(mg["time"]))
    task = [str(ctx.submit_script), "__merge", "--run-dir", str(ctx.run_dir), "--phase", "scan",
            "--shards", str(scan["n_shards"]), "--raw-root", ctx.raw_root, "--config", str(ctx.config),
            "--labels", ",".join(ctx.hpc["collector"]["extract_labels"])]
    prog.submit(name, args + dependency_arg("afterany", [i for i in ids if i], prog) + task)
    ctx.man.add_job(name, phase="merge", n_tasks=1, cluster=ctx.hpc["primary_cluster"], task=task)


def plan_extract(ctx: Ctx, prog: Program, source_dir: Path, after_merge: Optional[str] = None,
                 upper_bound: Optional[int] = None) -> None:
    """
    Phase B: one extraction task per qualifying .sld (skip-existing, checksums), %20 throttle,
    afterok on the merge. When the merge has not run yet (chained submission), the array is sized
    to the number of scanned files and tasks beyond the real list exit at once.
    """
    ex = ctx.hpc["extract"]
    col = ctx.hpc["collector"]
    extract_list = ctx.run_dir / "extract" / "extract_files.txt"
    shard_dir = ctx.run_dir / "extract" / "shards"
    deps: List[str] = []
    if after_merge:
        n_tasks = max(1, int(upper_bound or 1))
        time_min = parse_time_to_min(ex["time_unknown"])
        prog.note(f"Phase B  extract: list not known until the merge finishes -> {n_tasks} task(s) "
                  f"(one per scanned file; extra tasks exit at once), {fmt_time(time_min)} each")
        deps = [ph(after_merge)]
    else:
        inv_path = ctx.run_dir / "scan" / "merged" / "sld_4d_inventory.csv"
        inventory = read_inventory(inv_path)
        if not inventory:
            raise SystemExit(f"[ERROR] No merged inventory at {inv_path}: run 'scan' (and its merge) first.")
        files, per_file_gb = extract_selection(inventory, col["extract_labels"])
        if not ctx.dry_run:
            write_list(extract_list, files)  # (re)written so a --labels override takes effect
        n_tasks = len(files)
        if n_tasks == 0:
            prog.note("Phase B  extract: no qualifying files with an allowed mito label: nothing to submit.")
            return
        worst = max(per_file_gb.values())
        time_min = max(float(ctx.hpc["model"]["min_time_min"]),
                       (float(ex["sec_base"]) + float(ex["sec_per_gb"]) * worst)
                       * float(ctx.hpc["model"]["time_factor"]) / 60.0)
        prog.note(f"Phase B  extract: {n_tasks} file(s), {sum(per_file_gb.values()):.1f} GB of 4D series "
                  f"(largest file {worst:.2f} GB -> {fmt_time(time_min)} per task)")
        merge_jobs = [j.get("job_id", "") for n, j in ctx.man.jobs("merge").items()
                      if j.get("job_id") and "_extract_merge" not in n]  # scan merges only
        deps = merge_jobs
    if ctx.dry_run and not after_merge:
        prog.note(f"(dry run: {extract_list} is written by --execute)")
    throttle = min(int(ctx.hpc["extract_throttle"]), n_tasks)
    prog.mkdirs(ctx.log_dir, shard_dir, Path(ctx.raw_root))
    name = ctx.job_name("extract")
    task = [str(ctx.submit_script), "__extract", "--run-dir", str(ctx.run_dir), "--source-dir", str(source_dir),
            "--file-list", str(extract_list), "--shards", str(n_tasks), "--raw-root", ctx.raw_root,
            "--config", str(ctx.config)]
    args = ctx.common_args(name, ctx.hpc["primary_partition"], int(ex["cpus"]), int(ex["mem_gb"]), time_min,
                           array=f"1-{n_tasks}%{throttle}")
    prog.note(f"         compression {col['compression']}, at most {throttle} at once (I/O throttle), afterok on the merge")
    prog.submit(name, args + dependency_arg("afterok", deps, prog) + task)
    ctx.man.add_job(name, phase="extract", n_tasks=n_tasks, list=str(extract_list), cpus=int(ex["cpus"]),
                    mem_gb=int(ex["mem_gb"]), time_min=time_min, partition=ctx.hpc["primary_partition"],
                    array=f"1-{n_tasks}%{throttle}", cluster=ctx.hpc["primary_cluster"], task=task)

    mname = ctx.job_name("extract", "merge")
    mg = ctx.hpc["merge"]
    margs = ctx.common_args(mname, ctx.hpc["primary_partition"], int(mg["cpus"]), int(mg["mem_gb"]),
                            parse_time_to_min(mg["time"]))
    mtask = [str(ctx.submit_script), "__merge", "--run-dir", str(ctx.run_dir), "--phase", "extract",
             "--shards", str(n_tasks), "--raw-root", ctx.raw_root, "--config", str(ctx.config)]
    prog.note("         then merge the extraction checksums (afterany)")
    prog.submit(mname, margs + [f"--dependency=afterany:{ph(name)}"] + mtask)
    ctx.man.add_job(mname, phase="merge", n_tasks=1, cluster=ctx.hpc["primary_cluster"], task=mtask)


def extract_selection(inventory: Sequence[Dict[str, str]], labels: Sequence[str]) -> Tuple[List[str], Dict[str, float]]:
    """Unique .sld paths holding >= 1 qualifying 4D series with an allowed label, and their 4D GB."""
    per_file: Dict[str, float] = {}
    for r in inventory:
        if str(r.get("is_4d_timelapse", "")).strip().lower() != "true":
            continue
        if not label_allowed(r.get("mito_label", ""), labels):
            continue
        per_file[r["sld_path"]] = per_file.get(r["sld_path"], 0.0) + float(r.get("estimated_uncompressed_gb") or 0)
    return sorted(per_file), per_file


# ---------------------------------------------------------------- canary / stage3 / stage45


def resolve_movie_inputs(ctx: Ctx, movie_list: Optional[Path], inventory_path: Optional[Path]
                         ) -> Tuple[Path, List[Dict[str, str]]]:
    """Default movie list / inventory: the merged scan outputs; else the inventory next to the list."""
    merged = ctx.run_dir / "scan" / "merged"
    ml = movie_list or merged / "movie_list.txt"
    if not ml.exists():
        raise SystemExit(f"[ERROR] Movie list not found: {ml} (pass --movie-list)")
    inv = inventory_path
    if inv is None:
        for cand in (merged / "sld_4d_inventory.csv", ml.parent / "sld_4d_inventory.csv"):
            if cand.exists():
                inv = cand
                break
    return ml, read_inventory(inv) if inv else []


def class_partition(ctx: Ctx, cls: str) -> str:
    large = ctx.hpc["size_classes"].get("large_class", "L")
    return ctx.hpc["large_partition"] if cls == large else ctx.hpc["primary_partition"]


def model_with_refit(ctx: Ctx) -> Tuple[Dict[str, Any], float, float]:
    """The model plus the canary refit scales stored in the manifest (1.0 when not refitted)."""
    refit = ctx.man.data.get("model_refit") or {}
    return ctx.hpc["model"], float(refit.get("mem_scale", 1.0)), float(refit.get("time_scale", 1.0))


def plan_pipeline_arrays(ctx: Ctx, prog: Program, phase: str, movies: Sequence[Movie]) -> List[str]:
    """
    Writes per-class, per-chunk movie lists and submits one Stage 3 array per list (canary or full).
    Requests per class = the largest prediction in the class; throttles share 20% of each
    partition's cores between the arrays on that partition.
    """
    model, mem_scale, time_scale = model_with_refit(ctx)
    cpus = int(model["stage3_cpus"])
    max_arr = int(ctx.hpc["max_array_size"])
    by_class: Dict[str, List[Movie]] = {}
    for m in movies:
        by_class.setdefault(m.size_class, []).append(m)
    arrays = []  # (name, cls, chunk_idx, movies, partition)
    for cls in ctx.hpc["size_classes"]["names"]:
        ms = sorted(by_class.get(cls, []), key=lambda m: m.path)
        for ci, part in enumerate(chunk(ms, max_arr), 1):
            arrays.append((ctx.job_name(phase, cls, f"c{ci}"), cls, ci, part, class_partition(ctx, cls)))
    if not arrays:
        prog.note(f"{phase}: no movies to run.")
        return []
    pinfo = {}
    for p in sorted({a[4] for a in arrays}):
        exclude_gpu = p == ctx.hpc["large_partition"] and bool(ctx.hpc["large_exclude_gpu_nodes"])
        pinfo[p] = partition_info(p, exclude_gpu, ctx.hpc["primary_cluster"], ctx.hpc["primary_cluster"], ctx.hpc)
    throttles: Dict[str, int] = {}
    for p, info in pinfo.items():
        on_p = [a for a in arrays if a[4] == p]
        cap = cap_tasks(info["cores"], ctx.hpc["footprint_fraction"], cpus)
        for a, t in zip(on_p, split_throttle(cap, [len(a[3]) for a in on_p])):
            throttles[a[0]] = t
        prog.note(f"partition {p}: {info['cores']} cores ({info['source']}), cap "
                  f"{float(ctx.hpc['footprint_fraction']):.0%} -> {cap} concurrent {cpus}-core task(s)"
                  + (f", excluding GPU nodes {','.join(info['exclude'])}" if info["exclude"] else ""))
    deps: List[str] = []
    if phase == "stage3":
        canary_ids = [j["job_id"] for j in ctx.man.jobs("canary").values() if j.get("job_id")]
        deps = canary_ids
        still_running = [i for i in active_job_ids(canary_ids) if not i.startswith("DRY_")]
        if still_running:
            prog.note(f"[WARN] canary {','.join(still_running)} still running: the model refit used partial "
                      "numbers; Stage 3 waits for it (afterany). Re-plan after it finishes for exact requests.")
    names = []
    for name, cls, ci, part, partition in arrays:
        preds = [predict_stage3(m, model) for m in part]
        reqs = [to_request(mem, sec, model, mem_scale, time_scale) for mem, sec in preds]
        mem_req = max(r[0] for r in reqs)
        time_req = max(r[1] for r in reqs)
        info = pinfo[partition]
        if info["max_mem_gb"] and mem_req > info["max_mem_gb"]:
            prog.note(f"[WARN] {name}: {mem_req} GB exceeds the largest {partition} node "
                      f"({info['max_mem_gb']:.0f} GB); raise size_classes edges or move these movies")
        lst = ctx.state_dir / "lists" / f"{phase}_{cls}_c{ci}.txt"
        write_list(lst, [m.entry for m in part])
        thr = throttles[name]
        gbs = [m.gb for m in part]
        prog.note(f"{phase} class {cls} chunk {ci}: {len(part)} movie(s), {min(gbs):.2f}-{max(gbs):.2f} GB "
                  f"-> {cpus} cpu, {mem_req} GB, {fmt_time(time_req)}, %{thr}")
        args = ctx.common_args(name, partition, cpus, mem_req, time_req, exclude=info["exclude"],
                               array=f"1-{len(part)}%{thr}")
        task = [str(ctx.repo / "scripts" / "run_pipeline_slurm.bash"), str(ctx.config), str(lst), "--stages", "3"]
        prog.submit(name, args + dependency_arg("afterany", deps, prog) + task)
        ctx.man.add_job(name, phase=phase, size_class=cls, chunk=ci, n_tasks=len(part), list=str(lst),
                        cpus=cpus, mem_gb=mem_req, time_min=time_req, partition=partition,
                        exclude=info["exclude"], array=f"1-{len(part)}%{thr}", throttle=thr,
                        cluster=ctx.hpc["primary_cluster"], stages="3", task=task)
        names.append(name)
    return names


def plan_canary(ctx: Ctx, prog: Program, movie_list: Optional[Path], inventory_path: Optional[Path]) -> None:
    """Phase C: Stage 3 on ~10 movies spanning the size classes (measures the real model)."""
    ml, inventory = resolve_movie_inputs(ctx, movie_list, inventory_path)
    movies, missing, excluded = build_movies(read_list(ml), inventory, ctx.raw_root, ctx.hpc["model"],
                                             ctx.hpc["size_classes"], ctx.hpc["collector"]["extract_labels"])
    report_inputs(prog, ml, movies, missing, excluded)
    picked = select_canary(movies, int(ctx.hpc["canary_size"]))
    ctx.man.data["phases"]["canary"] = {"movie_list": str(ml), "movies": [m.path for m in picked]}
    prog.note(f"Phase C  canary: {len(picked)} movie(s) spanning classes "
              + ", ".join(f"{c}={sum(1 for m in picked if m.size_class == c)}" for c in ctx.hpc["size_classes"]["names"]))
    plan_pipeline_arrays(ctx, prog, "canary", picked)
    prog.note("After the canary: check QC (scripts/qc_pipeline_outputs.py), retrain ilastik on QC-PASS "
              "outputs, then run 'stage3' (it refits the model from the canary's sacct numbers).")


def refit_from_canary(ctx: Ctx, prog: Program, sacct_text: Optional[str]) -> None:
    """
    Compares the canary's measured MaxRSS / Elapsed with the unscaled model and stores
    mem_scale / time_scale (max observed/predicted ratio, clamped) in the manifest.
    """
    canary = ctx.man.jobs("canary")
    ids = [j["job_id"] for j in canary.values() if j.get("job_id") and not j["job_id"].startswith("DRY_")]
    if sacct_text is None:
        if not ids:
            prog.note("model: no completed canary in the manifest -> using the seeded model unchanged")
            return
        sacct_text = run_cmd(["sacct", "-P", "-n", "--units=M", "-o", "JobID,State,ExitCode,Elapsed,MaxRSS",
                              "-j", ",".join(ids)]) or ""
    tasks = parse_sacct(sacct_text)
    model = ctx.hpc["model"]
    mem_ratios, time_ratios = [], []
    inventory = read_inventory(ctx.run_dir / "scan" / "merged" / "sld_4d_inventory.csv")
    for name, job in canary.items():
        entries = read_list(Path(job["list"])) if Path(job.get("list", "")).exists() else []
        movies, _, _ = build_movies(entries, inventory, ctx.raw_root, model, ctx.hpc["size_classes"])
        by_path = {m.path: m for m in movies}
        for (jid, idx), rec in tasks.items():
            if jid != job.get("job_id") or idx < 1 or idx > len(entries):
                continue
            # Completed tasks and QC failures (exit 3: Stage 3 ran to the end) are valid measurements.
            if not (rec["state"] == "COMPLETED" or (rec["state"] == "FAILED" and rec["exit"] == QC_FAIL_EXIT)):
                continue
            m = by_path.get(entries[idx - 1])
            if m is None or m.gb <= 0:
                continue
            mem, sec = predict_stage3(m, model)
            if rec["maxrss_mb"] > 0:
                mem_ratios.append(rec["maxrss_mb"] / 1024.0 / mem)
            if rec["elapsed_s"] > 0 and sec > 0:
                time_ratios.append(rec["elapsed_s"] / sec)
    lo, hi = float(model["refit_min_scale"]), float(model["refit_max_scale"])
    clamp = lambda v: min(hi, max(lo, v))  # noqa: E731
    if not mem_ratios and not time_ratios:
        prog.note("model: canary has no usable measurements yet -> seeded model unchanged")
        return
    refit = {
        "mem_scale": round(clamp(max(mem_ratios)), 3) if mem_ratios else 1.0,
        "time_scale": round(clamp(max(time_ratios)), 3) if time_ratios else 1.0,
        "n_tasks_used": max(len(mem_ratios), len(time_ratios)),
        "raw_mem_ratio_max": round(max(mem_ratios), 3) if mem_ratios else None,
        "raw_time_ratio_max": round(max(time_ratios), 3) if time_ratios else None,
        "fitted": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    ctx.man.data["model_refit"] = refit
    prog.note(f"model refit from {refit['n_tasks_used']} canary task(s): memory x{refit['mem_scale']}, "
              f"time x{refit['time_scale']} (observed/predicted max, clamped to [{lo}, {hi}])")


def report_inputs(prog: Program, ml: Path, movies: Sequence[Movie], missing: Sequence[str], excluded: Sequence[str]) -> None:
    prog.note(f"movie list {ml}: {len(movies)} movie(s) found"
              + (f", {len(missing)} not on disk (not extracted?) skipped" if missing else "")
              + (f", {len(excluded)} excluded by mito label" if excluded else ""))
    for p in list(missing)[:5]:
        prog.note(f"  missing: {p}")


def plan_stage3(ctx: Ctx, prog: Program, movie_list: Optional[Path], inventory_path: Optional[Path],
                sacct_text: Optional[str], no_refit: bool) -> None:
    """Phase D: Stage 3 for every movie, arrays per size class (canary movies skip via .done markers)."""
    if not no_refit:
        refit_from_canary(ctx, prog, sacct_text)
    ml, inventory = resolve_movie_inputs(ctx, movie_list, inventory_path)
    movies, missing, excluded = build_movies(read_list(ml), inventory, ctx.raw_root, ctx.hpc["model"],
                                             ctx.hpc["size_classes"], ctx.hpc["collector"]["extract_labels"])
    report_inputs(prog, ml, movies, missing, excluded)
    ctx.man.data["phases"]["stage3"] = {"movie_list": str(ml), "n_movies": len(movies)}
    prog.note(f"Phase D  Stage 3: {len(movies)} movie(s); canary movies finish in seconds (.done markers)")
    plan_pipeline_arrays(ctx, prog, "stage3", movies)
    prog.note("Stage 4-5 arrays: run 'stage45' once the ilastik project is retrained (human gate).")


def ilp_path_from_config(ctx: Ctx) -> Path:
    """models_dir/<project_file> as run_pipeline_slurm.bash resolves it (MITO_MODELS_DIR)."""
    models = os.environ.get("MITO_MODELS_DIR") or os.path.join(
        os.environ.get("MITO_DATA_ROOT") or f"/sciclone/data10/{os.environ.get('USER', '')}/mito_pipeline", "models")
    text = ctx.config.read_text(encoding="utf-8")
    m = re.search(r'^\s*project_file:\s*"?([^"#\s]+)"?', text, re.M)
    return Path(models) / (m.group(1) if m else "bestman_mito_2d_v1.ilp")


def plan_stage45(ctx: Ctx, prog: Program, allow_fallback: bool) -> None:
    """
    Phase F: Stages 4-5 per movie, one array per Stage 3 array (same list => same index), with
    aftercorr so task i starts as soon as Stage 3 task i succeeded; --kill-on-invalid-dep cancels
    the tasks whose Stage 3 failed (e.g. QC exit 3) instead of leaving them pending forever.
    """
    s3 = ctx.man.jobs("stage3")
    if not s3:
        raise SystemExit("[ERROR] No stage3 arrays in the manifest: run 'stage3' first.")
    ilp = ilp_path_from_config(ctx)
    if ilp_fingerprint(ilp) == "no-ilp" and not allow_fallback:
        msg = f"no valid ilastik project at {ilp} (human gate: retrain on canary QC-PASS outputs first)"
        if ctx.dry_run:
            prog.note(f"[WARN] {msg}; --execute will refuse (or pass --allow-fallback)")
        else:
            raise SystemExit(f"[ERROR] {msg}; or pass --allow-fallback for a test run")
    model, mem_scale, time_scale = model_with_refit(ctx)
    cpus = int(model["stage45_cpus"])
    inventory = read_inventory(ctx.run_dir / "scan" / "merged" / "sld_4d_inventory.csv")
    pinfo = {}
    for p in sorted({j["partition"] for j in s3.values()}):
        exclude_gpu = p == ctx.hpc["large_partition"] and bool(ctx.hpc["large_exclude_gpu_nodes"])
        pinfo[p] = partition_info(p, exclude_gpu, ctx.hpc["primary_cluster"], ctx.hpc["primary_cluster"], ctx.hpc)
    # Stage 3 arrays still running hold part of the cap: Stage 4-5 gets the rest (min 1 task each).
    active3 = set(active_job_ids([j["job_id"] for j in s3.values() if j.get("job_id")]))
    throttles: Dict[str, int] = {}
    full_throttles: Dict[str, int] = {}
    for p, info in pinfo.items():
        cap_cores = int(math.floor(info["cores"] * float(ctx.hpc["footprint_fraction"])))
        busy = sum(int(j.get("throttle", 0)) * int(j.get("cpus", 1)) for j in s3.values()
                   if j.get("partition") == p and j.get("job_id") in active3)
        cap = max(1, (cap_cores - busy) // cpus)
        on_p = [n for n, j in s3.items() if j["partition"] == p]
        sizes = [int(s3[x]["n_tasks"]) for x in on_p]
        for n, t, full in zip(on_p, split_throttle(cap, sizes), split_throttle(max(1, cap_cores // cpus), sizes)):
            throttles[n] = t
            full_throttles[n] = full
        prog.note(f"partition {p}: {info['cores']} cores ({info['source']}), cap {cap_cores} cores"
                  + (f" minus {busy} held by Stage 3 arrays still queued/running" if busy else "")
                  + f" -> {cap} concurrent {cpus}-core task(s)")
    prog.note(f"Phase F  Stages 4-5: {len(s3)} array(s), aftercorr on the matching Stage 3 array")
    for s3name, job in sorted(s3.items()):
        entries = read_list(Path(job["list"]))
        movies, _, _ = build_movies(entries, inventory, ctx.raw_root, model, ctx.hpc["size_classes"])
        reqs = [to_request(*predict_stage45(m, model, cpus), model, mem_scale, time_scale) for m in movies] or [(int(model["min_mem_gb"]), int(model["min_time_min"]))]
        mem_req = max(r[0] for r in reqs)
        time_req = max(r[1] for r in reqs)
        name = s3name.replace("_stage3_", "_stage45_", 1)
        thr = min(throttles[s3name], int(job["n_tasks"]))
        prog.note(f"stage45 class {job['size_class']} chunk {job['chunk']}: {job['n_tasks']} movie(s) -> "
                  f"{cpus} cpu, {mem_req} GB, {fmt_time(time_req)}, %{thr}")
        args = ctx.common_args(name, job["partition"], cpus, mem_req, time_req, exclude=job.get("exclude", []),
                               array=f"1-{job['n_tasks']}%{thr}")
        dep = dependency_arg("aftercorr", [job["job_id"]] if job.get("job_id") else [], prog)
        if dep:
            dep.append("--kill-on-invalid-dep=yes")
        if allow_fallback:
            args[-1] += ",MITO_ALLOW_ILASTIK_FALLBACK=1"
        task = [str(ctx.repo / "scripts" / "run_pipeline_slurm.bash"), str(ctx.config), job["list"], "--stages", "4,5"]
        prog.submit(name, args + dep + task)
        ctx.man.add_job(name, phase="stage45", size_class=job["size_class"], chunk=job["chunk"],
                        n_tasks=int(job["n_tasks"]), list=job["list"], cpus=cpus, mem_gb=mem_req,
                        time_min=time_req, partition=job["partition"], exclude=job.get("exclude", []),
                        array=f"1-{job['n_tasks']}%{thr}", throttle=thr, cluster=job.get("cluster"),
                        stages="4,5", pair=s3name, task=task, allow_fallback=allow_fallback)
        full = min(full_throttles[s3name], int(job["n_tasks"]))
        if full > thr:
            # Stage 3 still holds part of the 20% cap; once it has finished, Stage 4-5 may use all of it.
            prog.note(f"  when Stage 3 has finished, raise it: scontrol update JobId=$(squeue -h -u $USER "
                      f"-n {name} -o %F | head -1) ArrayTaskThrottle={full}")


def plan_aggregate(ctx: Ctx, prog: Program) -> None:
    """Phase G: one aggregation job afterany on every Stage 4-5 array (and its retries)."""
    ids = []
    for j in ctx.man.jobs("stage45").values():
        if j.get("job_id"):
            ids.append(j["job_id"])
        ids.extend(r["job_id"] for r in j.get("retries", []) if r.get("job_id"))
    ag = ctx.hpc["aggregate"]
    name = ctx.job_name("aggregate")
    out_dir = ctx.run_dir / "results"
    prog.note(f"Phase G  aggregate -> {out_dir} (afterany on {len(ids)} Stage 4-5 job(s))")
    prog.mkdirs(ctx.log_dir)
    args = ctx.common_args(name, ctx.hpc["primary_partition"], int(ag["cpus"]), int(ag["mem_gb"]),
                           parse_time_to_min(ag["time"]))
    task = [str(ctx.submit_script), "__aggregate", "--run-dir", str(ctx.run_dir), "--config", str(ctx.config)]
    prog.submit(name, args + dependency_arg("afterany", ids, prog) + task)
    ctx.man.add_job(name, phase="aggregate", n_tasks=1, cluster=ctx.hpc["primary_cluster"], task=task)


# ---------------------------------------------------------------- retry


def latest_task_states(job: Dict[str, Any], sacct_text: str) -> Dict[int, Dict[str, Any]]:
    """Per array index, the record from the most recent attempt (original job, then retries in order)."""
    tasks = parse_sacct(sacct_text)
    attempts = [job.get("job_id", "")] + [r.get("job_id", "") for r in job.get("retries", [])]
    states: Dict[int, Dict[str, Any]] = {}
    for jid in attempts:
        for (tj, idx), rec in tasks.items():
            if tj == jid and rec["state"]:
                states[idx] = rec
    return states


def retry_indices(phase: str, states: Dict[int, Dict[str, Any]],
                  pair_states: Optional[Dict[int, Dict[str, Any]]] = None) -> List[int]:
    """
    Indices to resubmit: FAILED / OUT_OF_MEMORY / TIMEOUT / NODE_FAIL... but never a pipeline QC
    failure (exit 3). For Stage 4-5 also tasks CANCELLED before starting (killed because their
    Stage 3 task failed) whose Stage 3 has since COMPLETED (e.g. after a Stage 3 retry).
    """
    out = []
    for idx, rec in sorted(states.items()):
        if rec["state"] in RETRY_STATES:
            if phase in PIPELINE_PHASES and rec["exit"] in (QC_FAIL_EXIT, FG_FAIL_EXIT) and rec["state"] == "FAILED":
                continue
            out.append(idx)
        elif (phase == "stage45" and rec["state"] == "CANCELLED" and rec["elapsed_s"] == 0
              and pair_states is not None and (pair_states.get(idx) or {}).get("state") == "COMPLETED"):
            out.append(idx)
    return out


def compress_indices(idx: Sequence[int]) -> str:
    """[1,2,3,7,9,10] -> '1-3,7,9-10' (Slurm --array syntax)."""
    out, start, prev = [], None, None
    for i in sorted(set(idx)):
        if start is None:
            start = prev = i
        elif i == prev + 1:
            prev = i
        else:
            out.append(f"{start}-{prev}" if prev != start else f"{start}")
            start = prev = i
    if start is not None:
        out.append(f"{start}-{prev}" if prev != start else f"{start}")
    return ",".join(out)


def plan_retry(ctx: Ctx, prog: Program, phase: str, sacct_text: Optional[str]) -> None:
    """Phase H: resubmit only failed indices with memory/time x retry_factor per attempt (max_retries)."""
    jobs = ctx.man.jobs(phase)
    if not jobs:
        raise SystemExit(f"[ERROR] No '{phase}' jobs in the manifest.")
    all_ids = []
    for j in ctx.man.data["jobs"].values():
        all_ids += [j.get("job_id", "")] + [r.get("job_id", "") for r in j.get("retries", [])]
    all_ids = [i for i in all_ids if i and not i.startswith("DRY_")]
    if sacct_text is None:
        sacct_text = run_cmd(["sacct", "-P", "-n", "--units=M", "-o", "JobID,State,ExitCode,Elapsed,MaxRSS",
                              "-j", ",".join(all_ids)]) if all_ids else ""
        sacct_text = sacct_text or ""
    factor = float(ctx.hpc["retry_factor"])
    max_retries = int(ctx.hpc["max_retries"])
    any_submitted = False
    for name, job in sorted(jobs.items()):
        states = latest_task_states(job, sacct_text)
        pair_states = None
        if phase == "stage45" and job.get("pair") in ctx.man.data["jobs"]:
            pair_states = latest_task_states(ctx.man.data["jobs"][job["pair"]], sacct_text)
        idx = retry_indices(phase, states, pair_states)
        if not idx:
            prog.note(f"{name}: nothing to retry")
            continue
        attempt = len(job.get("retries", [])) + 1
        if attempt > max_retries:
            prog.note(f"{name}: {len(idx)} failed task(s) {compress_indices(idx)} but max_retries={max_retries} "
                      "reached; inspect the logs")
            continue
        mem = int(math.ceil(float(job["mem_gb"]) * factor ** attempt))
        tmin = float(job["time_min"]) * factor ** attempt
        array = compress_indices(idx)
        thr = min(len(idx), int(job.get("throttle") or len(idx)))
        rname = f"{name}_r{attempt}"
        reasons = {}
        for i in idx:
            reasons[states[i]["state"]] = reasons.get(states[i]["state"], 0) + 1
        prog.note(f"{name}: retry {attempt}/{max_retries} of {len(idx)} task(s) [{array}] "
                  f"({', '.join(f'{k}={v}' for k, v in sorted(reasons.items()))}) -> {mem} GB, {fmt_time(tmin)}")
        args = ctx.common_args(rname, job["partition"], int(job["cpus"]), mem, tmin,
                               exclude=job.get("exclude", []), array=f"{array}%{thr}")
        deps: List[str] = []
        if phase == "stage45" and job.get("pair") in ctx.man.data["jobs"]:
            pair = ctx.man.data["jobs"][job["pair"]]
            deps = [r["job_id"] for r in pair.get("retries", []) if r.get("job_id")]
        if job.get("allow_fallback"):
            args[-1] += ",MITO_ALLOW_ILASTIK_FALLBACK=1"
        prog.mkdirs(ctx.log_dir)
        prog.submit(rname, args + dependency_arg("afterany", deps, prog) + list(job["task"]),
                    cluster=job.get("cluster") if job.get("cluster") != ctx.hpc["primary_cluster"] else None)
        job.setdefault("retries", []).append({"name": rname, "job_id": "", "indices": array, "mem_gb": mem,
                                              "time_min": tmin, "attempt": attempt})
        any_submitted = True
    if not any_submitted:
        prog.note(f"retry {phase}: no task needs a retry.")


# ---------------------------------------------------------------- overflow


def plan_overflow(ctx: Ctx, prog: Program, enabled: bool) -> None:
    """
    Opt-in overflow: report every candidate cluster (deny list, live PreemptMode, idle cores, env
    smoke status); when --overflow is given and bora is saturated (our tasks pending longer than the
    threshold), move whole not-yet-started Stage 3 (+ paired 4-5) chunks to an eligible cluster.
    """
    ov = ctx.hpc.get("overflow") or {}
    clusters = list(ov.get("clusters") or [])
    deny = set(ov.get("deny_clusters") or [])
    primary = ctx.hpc["primary_cluster"]
    if not clusters:
        prog.note("overflow: hpc_run.overflow.clusters is empty -> overflow disabled")
        return
    need_mode = str(ov.get("require_preempt_mode", "OFF")).upper()
    eligible = []
    smoke = ctx.man.data.get("overflow_smoke", {})
    for c in clusters:
        part = (ov.get("partitions") or {}).get(c, "batch")
        if c in deny:
            prog.note(f"overflow {c}: REFUSED (on deny list)")
            continue
        if c == primary:
            continue
        out = run_on_cluster(c, ["scontrol", "show", "partition", part], primary)
        mode = (re.search(r"PreemptMode=(\S+)", out or "") or [None, "UNKNOWN"])[1]
        sinfo = run_on_cluster(c, ["sinfo", "-h", "-p", part, "-o", "%C"], primary) or ""
        alloc_idle = sinfo.strip().splitlines()[0].split("/") if sinfo.strip() else []
        idle = int(alloc_idle[1]) if len(alloc_idle) >= 4 else 0
        total = int(alloc_idle[3]) if len(alloc_idle) >= 4 else 0
        status = smoke.get(c, {}).get("result", "NOT_RUN")
        if mode.upper() != need_mode:
            prog.note(f"overflow {c}/{part}: REFUSED (PreemptMode={mode}, required {need_mode})")
            continue
        prog.note(f"overflow {c}/{part}: PreemptMode={mode}, {idle}/{total} cores idle, env smoke test: {status}")
        if status != "PASS":
            prog.note(f"  run the smoke test first: submit_full_run.bash overflow --overflow --smoke --execute "
                      f"(if it fails: bash scripts/setup_hpc_env.bash --env-name mito4d_env_{c} on a {c} node)")
            prog.smoke(c, part)
            continue
        eligible.append((c, part, idle, total))
    if not enabled:
        prog.note("overflow: report only (pass --overflow to plan moves)")
        return
    if not eligible:
        prog.note("overflow: no eligible cluster (all refused or smoke test not passed)")
        return
    # Saturation: our pending tasks on the primary cluster and how long they have waited.
    out = run_cmd(["squeue", "-h", "-r", "-u", os.environ.get("USER", ""), "-t", "PD,R,CG", "-o", "%F|%i|%T|%V"]) or ""
    now = datetime.datetime.now()
    oldest_wait = 0.0
    pending: Dict[str, int] = {}
    started: Dict[str, int] = {}
    for line in out.splitlines():
        parts = line.strip().split("|")
        if len(parts) < 4:
            continue
        arr, _, state, submit = parts[:4]
        if state == "PENDING":
            pending[arr] = pending.get(arr, 0) + 1
            try:
                wait = (now - datetime.datetime.fromisoformat(submit)).total_seconds() / 60.0
                oldest_wait = max(oldest_wait, wait)
            except ValueError:
                pass
        else:
            started[arr] = started.get(arr, 0) + 1
    threshold = float(ov.get("pending_threshold_min", 30))
    if oldest_wait < threshold:
        prog.note(f"overflow: {primary} not saturated (oldest pending task waited {oldest_wait:.0f} min "
                  f"< {threshold:.0f} min) -> nothing moved")
        return
    model, mem_scale, time_scale = model_with_refit(ctx)
    ilp_ok = ilp_fingerprint(ilp_path_from_config(ctx)) != "no-ilp"
    moves = 0
    for name, job in sorted(ctx.man.jobs("stage3").items()):
        jid = job.get("job_id", "")
        if not jid or job.get("moved_to") or started.get(jid) or pending.get(jid, 0) < int(job["n_tasks"]):
            continue  # some task already ran/started -> leave the chunk where it is
        partner = next((n for n, j in ctx.man.jobs("stage45").items() if j.get("pair") == name), None)
        pj = ctx.man.data["jobs"].get(partner) if partner else None
        if pj and (started.get(pj.get("job_id", "")) or pending.get(pj.get("job_id", ""), 0) < int(pj["n_tasks"])):
            continue
        c, part, idle, total = eligible[moves % len(eligible)]
        with_45 = bool(pj) and (ilp_ok or pj.get("allow_fallback"))
        stages = "3,4,5" if with_45 else "3"
        cpus = int(model["stage45_cpus"]) if with_45 else int(model["stage3_cpus"])
        mem = max(int(job["mem_gb"]), int(pj["mem_gb"])) if with_45 else int(job["mem_gb"])
        tmin = float(job["time_min"]) + (float(pj["time_min"]) if with_45 else 0.0)
        thr = min(int(job["n_tasks"]), cap_tasks(total, ctx.hpc["footprint_fraction"], cpus))
        oname = f"{name}_ov_{c}"
        prog.note(f"move {name} ({job['n_tasks']} pending) -> {c}/{part}: stages {stages}, {cpus} cpu, "
                  f"{mem} GB, {fmt_time(tmin)}, %{thr}")
        prog.cancel(jid)
        if pj:
            prog.cancel(pj["job_id"])
        args = ctx.common_args(oname, part, cpus, mem, tmin, array=f"1-{job['n_tasks']}%{thr}")
        if pj and pj.get("allow_fallback"):
            args[-1] += ",MITO_ALLOW_ILASTIK_FALLBACK=1"
        task = [str(ctx.repo / "scripts" / "run_pipeline_slurm.bash"), str(ctx.config), job["list"], "--stages", stages]
        prog.submit(oname, args + task, cluster=c)
        ctx.man.add_job(oname, phase="stage3" if not with_45 else "overflow", size_class=job.get("size_class"),
                        chunk=job.get("chunk"), n_tasks=int(job["n_tasks"]), list=job["list"], cpus=cpus,
                        mem_gb=mem, time_min=tmin, partition=part, exclude=[], throttle=thr, cluster=c,
                        stages=stages, task=task, moved_from=name)
        job["moved_to"] = oname
        if pj:
            pj["moved_to"] = oname
        moves += 1
    if not moves:
        prog.note("overflow: no fully-pending chunk to move")


# ---------------------------------------------------------------- status


def status_report(man: Manifest) -> str:
    """One sacct call over every job id in the manifest, summarised per array."""
    ids = {}
    for name, j in man.data["jobs"].items():
        for jid, label in [(j.get("job_id", ""), name)] + [(r.get("job_id", ""), r.get("name", name))
                                                            for r in j.get("retries", [])]:
            if jid and not jid.startswith("DRY_"):
                ids[jid] = label
    if not ids:
        return "No submitted jobs in the manifest (dry run only?).\n"
    out = run_cmd(["sacct", "-X", "-P", "-n", "-o", "JobID,State,ExitCode", "-j", ",".join(ids)])
    if out is None:
        return "sacct is unavailable (module load slurm/bora?).\n"
    counts: Dict[str, Dict[str, int]] = {}
    for line in out.splitlines():
        parts = line.strip().split("|")
        if len(parts) < 3:
            continue
        base = re.split(r"[_.]", parts[0])[0]
        label = ids.get(base, base)
        m = re.match(r"\d+_\[(.+)\]", parts[0])
        n = 1
        if m:  # a pending range such as 123_[5-40%20]
            n = 0
            for rng in m.group(1).split("%")[0].split(","):
                a, _, b = rng.partition("-")
                n += (int(b) - int(a) + 1) if b else 1
        state = parts[1].split()[0]
        if state == "FAILED" and parts[2].startswith(f"{QC_FAIL_EXIT}:"):
            state = "QC_FAIL(exit3)"
        elif state == "FAILED" and parts[2].startswith(f"{FG_FAIL_EXIT}:"):
            state = "FG_IMPLAUSIBLE(exit5)"
        counts.setdefault(label, {})
        counts[label][state] = counts[label].get(state, 0) + n
    lines = [f"{'array':<48} states"]
    for label in sorted(counts):
        lines.append(f"{label:<48} " + ", ".join(f"{k}={v}" for k, v in sorted(counts[label].items())))
    return "\n".join(lines) + "\n"


# ================================================================================================
# 8. CLI
# ================================================================================================


def make_ctx(args: argparse.Namespace) -> Ctx:
    config = Path(args.config).resolve()
    run_dir = Path(args.run_dir).resolve()
    dry = bool(getattr(args, "dry_run", False))
    state_dir, man_path = run_dirs(run_dir, dry)
    hpc = load_hpc_run(config)
    user = os.environ.get("USER", "")
    log_root = hpc.get("log_root") or f"/sciclone/scr30/{user}/mito_logs"
    raw_root = getattr(args, "raw_root", None) or os.environ.get("MITO_RAW_ROOT", "")
    man = Manifest(man_path)
    init_manifest_provenance(man, run_dir, config, REPO_DIR, dry)
    man.data["log_dir"] = str(Path(log_root) / run_dir.name)
    if raw_root:
        man.data.setdefault("raw_root", raw_root)
    return Ctx(config=config, repo=REPO_DIR, run_dir=run_dir, log_dir=Path(log_root) / run_dir.name,
               state_dir=state_dir, man=man, hpc=hpc, dry_run=dry,
               raw_root=raw_root or man.data.get("raw_root", ""),
               submit_script=SCRIPT_DIR / "submit_full_run.bash")


def cmd_plan(args: argparse.Namespace) -> int:
    ctx = make_ctx(args)
    prog = Program()
    sacct_text = Path(args.sacct_file).read_text() if getattr(args, "sacct_file", None) else None
    phase = args.phase
    if args.labels:  # mito labels to extract (scan --chain-extract, merge, extract)
        ctx.hpc["collector"]["extract_labels"] = [s for s in args.labels.split(",") if s]
    if phase == "scan":
        if not args.source_dir:
            raise SystemExit("[ERROR] scan needs --source-dir")
        plan_scan(ctx, prog, Path(args.source_dir).resolve(), Path(args.file_list) if args.file_list else None,
                  Path(args.plasmid_csv).resolve() if args.plasmid_csv else None, args.chain_extract,
                  Path(args.mito_overrides).resolve() if args.mito_overrides else None)
    elif phase == "merge":
        plan_merge(ctx, prog)
    elif phase == "extract":
        src = args.source_dir or (ctx.man.data["phases"].get("scan") or {}).get("source_dir") \
            or (Manifest(run_dirs(ctx.run_dir, False)[1]).data["phases"].get("scan") or {}).get("source_dir")
        if not src:
            raise SystemExit("[ERROR] extract needs --source-dir (or a scan recorded in the manifest)")
        plan_extract(ctx, prog, Path(src).resolve())
    elif phase == "canary":
        plan_canary(ctx, prog, Path(args.movie_list) if args.movie_list else None,
                    Path(args.inventory) if args.inventory else None)
    elif phase == "stage3":
        plan_stage3(ctx, prog, Path(args.movie_list) if args.movie_list else None,
                    Path(args.inventory) if args.inventory else None, sacct_text, args.no_refit)
    elif phase == "stage45":
        plan_stage45(ctx, prog, args.allow_fallback)
    elif phase == "aggregate":
        plan_aggregate(ctx, prog)
    elif phase == "retry":
        plan_retry(ctx, prog, args.retry_phase, sacct_text)
    elif phase == "overflow":
        plan_overflow(ctx, prog, args.overflow)
    ctx.man.save()
    prog.note(f"manifest: {ctx.man.path}")
    sys.stdout.write(prog.text())
    return 0


def cmd_record(args: argparse.Namespace) -> int:
    man = Manifest(Path(args.manifest))
    job_id = args.job_id.split(";")[0].strip()
    for name, j in man.data["jobs"].items():
        if name == args.name:
            j["job_id"] = job_id
            j["submitted"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            break
        for r in j.get("retries", []):
            if r.get("name") == args.name:
                r["job_id"] = job_id
                r["submitted"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
                break
        else:
            continue
        break
    else:
        print(f"[WARN] {args.name} not in manifest; recorded under 'unplanned'", file=sys.stderr)
        man.data.setdefault("unplanned", {})[args.name] = job_id
    man.save()
    return 0


def cmd_smoke_record(args: argparse.Namespace) -> int:
    man = Manifest(Path(args.manifest))
    man.data.setdefault("overflow_smoke", {})[args.cluster] = {
        "result": args.result, "partition": args.partition,
        "when": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    man.save()
    return 0


def cmd_extract_list(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).resolve()
    hpc = load_hpc_run(Path(args.config))
    labels = [s for s in args.labels.split(",") if s] if args.labels else hpc["collector"]["extract_labels"]
    inventory = read_inventory(run_dir / "scan" / "merged" / "sld_4d_inventory.csv")
    files, per_file = extract_selection(inventory, labels)
    out = run_dir / "extract" / "extract_files.txt"
    write_list(out, files)
    print(f"extract list: {len(files)} file(s), {sum(per_file.values()):.1f} GB of 4D series -> {out}")
    return 0


def cmd_collector_args(args: argparse.Namespace) -> int:
    """Prints the collector filter flags from hpc_run.collector (identical for scan and extract)."""
    col = load_hpc_run(Path(args.config))["collector"]
    print(" ".join(shlex.quote(a) for a in [
        "--min-z", str(col["min_z"]), "--min-t", str(col["min_t"]), "--min-c", str(col["min_c"]),
        "--compression", str(col["compression"])]))
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    man = Manifest(run_dirs(Path(args.run_dir).resolve(), args.dry_run)[1])
    sys.stdout.write(status_report(man))
    return 0


def cmd_stage_hash(args: argparse.Namespace) -> int:
    hashes = compute_stage_hashes(Path(args.config), Path(args.repo), Path(args.input) if args.input else None,
                                  Path(args.ilp) if args.ilp else None)
    print(hashes[3], hashes[4], hashes[5])
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plan / record / report full-archive Slurm runs (see submit_full_run.bash).")
    sub = p.add_subparsers(dest="command", required=True)

    pl = sub.add_parser("plan", help="print the submission program for one phase")
    pl.add_argument("phase", choices=["scan", "merge", "extract", "canary", "stage3", "stage45", "aggregate", "retry", "overflow"])
    pl.add_argument("--config", default=str(DEFAULT_CONFIG))
    pl.add_argument("--run-dir", required=True)
    pl.add_argument("--dry-run", action="store_true", help="plan into <run-dir>/dryrun (no real submissions)")
    pl.add_argument("--raw-root", default=None, help="raw OME-TIFF root (default $MITO_RAW_ROOT)")
    pl.add_argument("--source-dir", default=None, help="scan/extract: directory holding the .sld files")
    pl.add_argument("--file-list", default=None, help="scan: explicit .sld list instead of searching source-dir")
    pl.add_argument("--plasmid-csv", default=None, help="scan: lab plasmid table for the mito-construct labels")
    pl.add_argument("--mito-overrides", default=None, help="scan: overrides CSV (pattern,label,note) for the labels")
    pl.add_argument("--chain-extract", action="store_true", help="scan: also submit extraction (afterok on merge)")
    pl.add_argument("--labels", default=None,
                    help="scan/merge/extract/canary/stage3: comma list of mito labels to extract and run (default from config)")
    pl.add_argument("--movie-list", default=None, help="canary/stage3: movie list (default: merged scan list)")
    pl.add_argument("--inventory", default=None, help="canary/stage3: sld_4d_inventory.csv for sizes")
    pl.add_argument("--no-refit", action="store_true", help="stage3: keep the seeded model")
    pl.add_argument("--sacct-file", default=None, help="stage3/retry: read sacct output from a file (offline/tests)")
    pl.add_argument("--allow-fallback", action="store_true", help="stage45: run without a trained .ilp (testing)")
    pl.add_argument("--retry-phase", default="stage3", choices=["scan", "extract", "canary", "stage3", "stage45"])
    pl.add_argument("--overflow", action="store_true", help="overflow: plan moves (otherwise report only)")
    pl.set_defaults(func=cmd_plan)

    rc = sub.add_parser("record", help="store a submitted job id in the manifest")
    rc.add_argument("--manifest", required=True)
    rc.add_argument("--name", required=True)
    rc.add_argument("--job-id", required=True)
    rc.set_defaults(func=cmd_record)

    sr = sub.add_parser("smoke-record", help="store a per-cluster smoke test result")
    sr.add_argument("--manifest", required=True)
    sr.add_argument("--cluster", required=True)
    sr.add_argument("--partition", required=True)
    sr.add_argument("--result", required=True, choices=["PASS", "FAIL"])
    sr.set_defaults(func=cmd_smoke_record)

    el = sub.add_parser("extract-list", help="write extract_files.txt from the merged inventory")
    el.add_argument("--config", default=str(DEFAULT_CONFIG))
    el.add_argument("--run-dir", required=True)
    el.add_argument("--labels", default=None, help="comma list of mito labels (default: config extract_labels)")
    el.set_defaults(func=cmd_extract_list)

    ca = sub.add_parser("collector-args", help="print the collector filter flags from the config")
    ca.add_argument("--config", default=str(DEFAULT_CONFIG))
    ca.set_defaults(func=cmd_collector_args)

    st = sub.add_parser("status", help="one sacct summary of the run")
    st.add_argument("--run-dir", required=True)
    st.add_argument("--dry-run", action="store_true")
    st.set_defaults(func=cmd_status)

    sh = sub.add_parser("stage-hash", help="print the Stage 3, 4, 5 hashes for .done markers")
    sh.add_argument("--config", required=True)
    sh.add_argument("--repo", default=str(REPO_DIR))
    sh.add_argument("--input", default=None)
    sh.add_argument("--ilp", default=None)
    sh.set_defaults(func=cmd_stage_hash)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
