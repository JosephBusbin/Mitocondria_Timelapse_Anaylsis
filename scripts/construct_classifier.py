#!/usr/bin/env python3
"""Classify microscopy time-lapse files as mitochondrial fluorescence from NAMES only.

The Bestman lab's 4D time-lapses (.sld SlideBook files, .tif exports, ...) are
named inconsistently. Before the expensive extraction pipeline runs, this tool
decides which files most likely image MITOCHONDRIAL fluorescence, using only
file/folder names (it never opens a data file):

* mito keywords in the name (``mls`` = Mitochondrial Localization Signal,
  ``mito``, ``cox8``, ``tomm20``);
* lab plasmid IDs embedded in the name (``obf01`` = OBF001, ``jr19`` = JR19,
  any zero padding, any case) looked up in the lab plasmid table (CSV);
* common mito dyes/reporters/proteins (TMRM, MitoTracker, Drp1, ...) which are
  never trusted automatically and are routed to human REVIEW (high priority).

Labels
------
MITO                   constitutive (or sensor) mito-targeted fluorescence.
MITO_PHOTOACTIVATABLE  only photoactivatable/photoconvertible mito FPs
                       (paGFP, dendra2, mEos2...): dark or color-changing.
NOT_MITO               plasmid IDs recognized, none mito, no mito keyword.
REVIEW                 nothing recognized, unknown IDs, conflicts, or only a
                       common mito marker (review_priority = high).

Outputs (in --out DIR): ``mito_classification.csv`` (every file) and
``mito_review.csv`` (REVIEW rows only, high priority first).

Standard library only; runs on Linux (HPC) and Windows (Python 3.10+).

Usage examples
--------------
Linux / HPC::

    module load python/3.12.7
    python3 scripts/construct_classifier.py \\
        --plasmid-csv /path/to/Plasmid_Primer_List_Bestman_Lab.csv \\
        --scan-dir /path/to/raw_movies --out results/mito_filter

    # classify a plain list of names (one per line)
    python3 scripts/construct_classifier.py --plasmid-csv plasmids.csv \\
        --names-file names.txt --out results/mito_filter

    # build a synthetic test tree (empty files + expected_labels.csv)
    python3 scripts/construct_classifier.py --plasmid-csv plasmids.csv \\
        --make-synthetic /local/$USER/synthetic_names --seed 7

Windows PowerShell (lab PC)::

    . "$env:LOCALAPPDATA\\mito_pipeline_tools\\activate_mito4d.ps1"
    python scripts\\construct_classifier.py `
        --plasmid-csv "D:\\lab\\Plasmid_Primer_List_Bestman_Lab.csv" `
        --scan-dir "D:\\SlideBook_data" --out "D:\\mito_filter_results" `
        --overrides "D:\\mito_filter_results\\overrides.csv"

The plasmid CSV can also be given with the environment variable
``MITO_PLASMID_CSV``. Exit code 0 on success, 2 on bad input.
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import json
import os
import random
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path, PurePath
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LEXICON_PATH = REPO_ROOT / "config" / "construct_lexicon.json"

LABELS = ("MITO", "MITO_PHOTOACTIVATABLE", "NOT_MITO", "REVIEW")
REPORT_COLUMNS = [
    "rel_path", "label", "review_priority", "mito_fp", "fp_class", "mito_color",
    "red_mito_flag", "matched_ids", "evidence", "reasons",
]
CLASSIFICATION_CSV = "mito_classification.csv"
REVIEW_CSV = "mito_review.csv"
EXPECTED_LABELS_CSV = "expected_labels.csv"

# Characters treated as separators between a marker and an adjacent FP name.
_SEP_CHARS = "-_: \t.:"


class ClassifierInputError(ValueError):
    """Raised for unusable user input (missing columns, bad override label...)."""


# --------------------------------------------------------------------------- #
# Lexicon
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class FPInfo:
    """A fluorescent protein: canonical name, class and display channel color."""

    canonical: str
    cls: str = "unknown"
    color: str = "unknown"
    note: str = ""


@dataclass
class KeywordEntry:
    """A direct keyword or common marker with its compiled alias regexes."""

    token: str
    evidence: str
    category: str
    rationale: str
    regexes: List["re.Pattern[str]"]


@dataclass
class Lexicon:
    """Parsed lexicon JSON plus compiled regular expressions."""

    raw: dict
    path: Optional[Path]
    columns: Dict[str, str]
    localization_synonyms: Dict[str, str]
    mito_localization: str
    infer_mito_from_blank: bool
    direct_keywords: List[KeywordEntry]
    common_markers: List[KeywordEntry]
    construct_markers: List[Tuple[str, "re.Pattern[str]"]]
    fp_aliases: List[Tuple["re.Pattern[str]", FPInfo]]
    constitutive_classes: Set[str]
    photo_classes: Set[str]
    numeric_min: int
    numeric_max: int
    numeric_prefix: str
    sep_prefix_min_len: int
    unknown_min_prefix_len: int
    unknown_ignore_prefixes: Set[str]
    file_extensions: List[str]
    skip_name_patterns: List[str]


def _alias_pattern(alias: str) -> str:
    """Regex source for an alias; '-', '_' or ' ' inside it match an optional separator."""
    parts = [re.escape(p) for p in re.split(r"[-_ ]+", alias.lower()) if p]
    return r"[-_ ]?".join(parts)


def _compile_alias(alias: str, boundary: str) -> "re.Pattern[str]":
    """Compile an alias with the requested letter-boundary rule (none|left|both)."""
    left = r"(?<![a-z])" if boundary in ("left", "both") else ""
    right = r"(?![a-z])" if boundary in ("right", "both") else ""
    return re.compile(left + _alias_pattern(alias) + right)


def _keyword_entries(items: Sequence[dict], default_prefix: str) -> List[KeywordEntry]:
    """Build KeywordEntry objects from the lexicon list form."""
    entries = []
    for item in items:
        token = str(item["token"]).lower()
        aliases = item.get("aliases") or [token]
        boundary = item.get("boundary", "none")
        if boundary not in ("none", "left", "right", "both"):
            raise ClassifierInputError(f"lexicon: bad boundary {boundary!r} for token {token!r}")
        entries.append(KeywordEntry(
            token=token,
            evidence=item.get("evidence", f"{default_prefix}:{token}"),
            category=item.get("category", "keyword"),
            rationale=item.get("rationale", ""),
            regexes=[_compile_alias(a, boundary) for a in aliases],
        ))
    return entries


def load_lexicon(path: Optional[os.PathLike] = None) -> Lexicon:
    """Load and compile the lexicon JSON (default: repo config/construct_lexicon.json)."""
    lex_path = Path(path) if path else DEFAULT_LEXICON_PATH
    try:
        with open(lex_path, "r", encoding="utf-8-sig") as fh:
            raw = json.load(fh)
    except FileNotFoundError as exc:
        raise ClassifierInputError(f"lexicon not found: {lex_path}") from exc
    except json.JSONDecodeError as exc:
        raise ClassifierInputError(f"lexicon is not valid JSON ({lex_path}): {exc}") from exc

    required = ["columns", "localization_synonyms", "direct_mito_keywords", "common_markers",
                "fp_vocabulary", "numeric_id", "file_extensions"]
    missing = [k for k in required if k not in raw]
    if missing:
        raise ClassifierInputError(f"lexicon {lex_path} is missing keys: {missing}")

    fp_aliases: List[Tuple[re.Pattern, FPInfo]] = []
    for item in raw["fp_vocabulary"]:
        info = FPInfo(item["canonical"], item.get("class", "unknown"),
                      item.get("color", "unknown"), item.get("note", ""))
        for alias in item["aliases"]:
            fp_aliases.append((re.compile(_alias_pattern(alias)), info))

    markers = [(m.lower(), re.compile(re.escape(m.lower())))
               for m in raw.get("construct_name_mito_markers", ["mls", "mito", "cox8", "tomm20"])]
    idm = raw.get("id_matching", {})
    num = raw["numeric_id"]
    return Lexicon(
        raw=raw,
        path=lex_path,
        columns={k: str(v) for k, v in raw["columns"].items()},
        localization_synonyms={str(k).strip().lower(): str(v)
                               for k, v in raw["localization_synonyms"].items()},
        mito_localization=raw.get("mito_localization", "mito"),
        infer_mito_from_blank=bool(raw.get("infer_mito_from_name_when_localization_blank", False)),
        direct_keywords=_keyword_entries(raw["direct_mito_keywords"], "keyword"),
        common_markers=_keyword_entries(raw["common_markers"], "common_marker"),
        construct_markers=markers,
        fp_aliases=fp_aliases,
        constitutive_classes=set(raw.get("mito_constitutive_classes", ["constitutive", "sensor"])),
        photo_classes=set(raw.get("mito_photo_classes", ["photoactivatable", "photoconvertible"])),
        numeric_min=int(num.get("min_digits", 5)),
        numeric_max=int(num.get("max_digits", 6)),
        numeric_prefix=str(num.get("evidence_prefix", "addgene")),
        sep_prefix_min_len=int(idm.get("allow_separator_for_prefix_min_len", 3)),
        unknown_min_prefix_len=int(idm.get("unknown_id_min_prefix_len", 2)),
        unknown_ignore_prefixes={p.lower() for p in idm.get("unknown_id_ignore_prefixes", ["p"])},
        file_extensions=[e.lower() for e in raw["file_extensions"]],
        skip_name_patterns=list(raw.get("skip_name_patterns", ["~$*", ".*"])),
    )


# --------------------------------------------------------------------------- #
# Fluorescent-protein lookup helpers
# --------------------------------------------------------------------------- #
def fp_at(text: str, pos: int, lexicon: Lexicon) -> Optional[Tuple[FPInfo, int]]:
    """Return the longest FP alias that starts exactly at ``pos`` in lowercase text."""
    best: Optional[Tuple[FPInfo, int]] = None
    for rx, info in lexicon.fp_aliases:
        m = rx.match(text, pos)
        if m and (best is None or m.end() > best[1]):
            best = (info, m.end())
    return best


def fp_ending_at(text: str, end: int, lexicon: Lexicon) -> Optional[FPInfo]:
    """Return the longest FP alias that ends exactly at ``end`` in lowercase text."""
    for start in range(max(0, end - 24), end):  # earliest start = longest alias
        for rx, info in lexicon.fp_aliases:
            if rx.fullmatch(text, start, end):
                return info
    return None


def lookup_fp_cell(value: str, lexicon: Lexicon) -> Optional[FPInfo]:
    """Interpret a table FP cell such as "tagRFPt (3')" or a multi-line "mTagBFP2\\n  "."""
    if not value or not value.strip():
        return None
    first = value.strip().splitlines()[0].strip()
    cleaned = re.sub(r"\(.*?\)", "", first).lower()
    cleaned = re.sub(r"\s+", "", cleaned)
    if not cleaned:
        return None
    hit = fp_at(cleaned, 0, lexicon)
    if hit:
        return hit[0]
    return FPInfo(canonical=first, cls="unknown", color="unknown")


def fp_adjacent_to(text: str, start: int, end: int, lexicon: Lexicon) -> Optional[FPInfo]:
    """FP right after a marker span (skipping separators/digits), else right before it."""
    pos = end
    while pos < len(text) and text[pos] in _SEP_CHARS:
        pos += 1
    hit = fp_at(text, pos, lexicon)
    if hit:
        return hit[0]
    # 'mito7EGFP' style: skip a short digit run after the marker.
    m = re.match(r"\d{1,2}[-_: .]*", text[pos:])
    if m:
        hit = fp_at(text, pos + m.end(), lexicon)
        if hit:
            return hit[0]
    back = start
    while back > 0 and text[back - 1] in _SEP_CHARS:
        back -= 1
    info = fp_ending_at(text, back, lexicon)
    if info:
        return info
    # 'CEPIA4mt' style: skip trailing digits before the marker.
    digits = re.search(r"\d{1,2}$", text[:back])
    if digits:
        return fp_ending_at(text, digits.start(), lexicon)
    return None


def mito_fp_from_construct_name(name: str, lexicon: Lexicon) -> Optional[Tuple[FPInfo, str]]:
    """Find the FP fused to a mito marker in a construct name (markers in priority order)."""
    low = name.lower()
    for token, rx in lexicon.construct_markers:
        for m in rx.finditer(low):
            info = fp_adjacent_to(low, m.start(), m.end(), lexicon)
            if info:
                return info, token
    return None


# --------------------------------------------------------------------------- #
# Plasmid table
# --------------------------------------------------------------------------- #
IdKey = Tuple[str, int, str]  # (letter prefix lowercase or '' for numeric, integer, suffix)

_ID_RX = re.compile(r"([a-z](?:[a-z0-9]*[a-z])?)[-_ ]?(\d+)(?:[-_](\d+(?:\.\d+)?))?")


def normalize_id(raw: str) -> Optional[IdKey]:
    """Normalize a plasmid ID to (prefix, int, suffix); OBF001 == obf01 == obf0001.

    Pure numbers (Addgene IDs) return ('', n, ''). Returns None when the text is
    not ID-shaped. Trailing '?' / '*' / whitespace (e.g. 'p116?') are ignored.
    """
    s = (raw or "").strip().lower()
    s = re.sub(r"[?*\s]+$", "", s)
    if not s:
        return None
    if s.isdigit():
        return ("", int(s), "")
    m = _ID_RX.fullmatch(s)
    if not m:
        return None
    return (m.group(1), int(m.group(2)), m.group(3) or "")


def format_key(key: IdKey) -> str:
    """Human-readable form of a normalized key (used for unknown IDs)."""
    prefix, num, suffix = key
    return f"{prefix}{num}" + (f"-{suffix}" if suffix else "")


@dataclass
class PlasmidRecord:
    """One plasmid (duplicate table IDs are merged into one record)."""

    key: IdKey
    table_ids: List[str] = field(default_factory=list)
    names: List[str] = field(default_factory=list)
    localizations_raw: List[str] = field(default_factory=list)
    localizations: Set[str] = field(default_factory=set)
    fps: List[FPInfo] = field(default_factory=list)
    is_mito: bool = False
    mito_fp: Optional[FPInfo] = None
    mito_fp_source: str = ""

    @property
    def display_id(self) -> str:
        return self.table_ids[0] if self.table_ids else format_key(self.key)

    @property
    def display_name(self) -> str:
        return " | ".join(dict.fromkeys(n for n in self.names if n))

    @property
    def localization_label(self) -> str:
        raws = [r if r else "blank" for r in dict.fromkeys(self.localizations_raw)]
        return "/".join(raws) if raws else "blank"


@dataclass
class ConstructTable:
    """Plasmid lookup structures built from the lab CSV."""

    records: Dict[IdKey, PlasmidRecord] = field(default_factory=dict)
    prefixes: Dict[str, Set[int]] = field(default_factory=dict)
    suffixed: Dict[Tuple[str, int], Set[str]] = field(default_factory=dict)
    numeric: Set[int] = field(default_factory=set)
    source: str = ""
    quirks: List[str] = field(default_factory=list)
    _prefix_rx: Dict[str, "re.Pattern[str]"] = field(default_factory=dict, repr=False)


def _read_csv_rows(csv_path: Path) -> List[List[str]]:
    """Read CSV rows, tolerating a UTF-8 BOM or Windows-1252 encoded exports."""
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            with open(csv_path, "r", encoding=enc, newline="") as fh:
                return list(csv.reader(fh))
        except UnicodeDecodeError:
            continue
    raise ClassifierInputError(f"cannot decode {csv_path}")


def _finalize_record(rec: PlasmidRecord, lexicon: Lexicon) -> None:
    """Decide whether a merged record is mito and which FP is the mito-targeted one."""
    rec.is_mito = lexicon.mito_localization in rec.localizations
    if not rec.is_mito and lexicon.infer_mito_from_blank and rec.localizations <= {"none"} \
            and all(not r for r in rec.localizations_raw):
        strong = {"mls", "mito", "mit0", "cox8", "tomm20", "tom20"}
        if any(t in n.lower() for n in rec.names for t in strong):
            rec.is_mito = True
    if not rec.is_mito:
        return
    for name in rec.names:
        found = mito_fp_from_construct_name(name, lexicon)
        if found:
            rec.mito_fp, marker = found
            rec.mito_fp_source = f"construct name ({marker})"
            return
    if rec.fps:
        rec.mito_fp = rec.fps[0]
        rec.mito_fp_source = "Fluorescent Protein 1"


def load_construct_table(csv_path: os.PathLike, lexicon: Lexicon) -> ConstructTable:
    """Read the lab plasmid CSV (read-only) and build the ID lookup table."""
    path = Path(csv_path)
    if not path.is_file():
        raise ClassifierInputError(f"plasmid CSV not found: {path}")
    rows = _read_csv_rows(path)
    if not rows:
        raise ClassifierInputError(f"plasmid CSV is empty: {path}")
    header = [h.strip().lower() for h in rows[0]]
    col_idx: Dict[str, int] = {}
    for logical, wanted in lexicon.columns.items():
        try:
            col_idx[logical] = header.index(wanted.strip().lower())
        except ValueError:
            if logical in ("id", "localization"):
                raise ClassifierInputError(
                    f"plasmid CSV {path} has no column {wanted!r} (lexicon columns.{logical}); "
                    f"header is {rows[0]}") from None
    table = ConstructTable(source=str(path))

    def cell(row: List[str], logical: str) -> str:
        i = col_idx.get(logical)
        return row[i].strip() if i is not None and i < len(row) else ""

    seen_raw: Dict[IdKey, int] = {}
    for line_no, row in enumerate(rows[1:], start=2):
        raw_id = cell(row, "id")
        if not raw_id:
            continue
        key = normalize_id(raw_id)
        if key is None:
            table.quirks.append(f"row {line_no}: unparseable plasmid ID {raw_id!r} skipped")
            continue
        if raw_id.strip() != raw_id.strip().rstrip("?*"):
            table.quirks.append(f"row {line_no}: ID {raw_id!r} has a trailing '?'; matched as {format_key(key)}")
        rec = table.records.setdefault(key, PlasmidRecord(key=key))
        seen_raw[key] = seen_raw.get(key, 0) + 1
        if seen_raw[key] == 2:
            table.quirks.append(f"duplicate plasmid ID {raw_id!r}: rows merged (union of evidence)")
        rec.table_ids.append(raw_id.strip())
        rec.names.append(" ".join(cell(row, "name").split()))
        loc_raw = " ".join(cell(row, "localization").lower().split())
        rec.localizations_raw.append(loc_raw)
        rec.localizations.add(lexicon.localization_synonyms.get(loc_raw, "other" if loc_raw else "none"))
        for logical in ("fp1", "fp2"):
            info = lookup_fp_cell(cell(row, logical), lexicon)
            if info and info not in rec.fps:
                rec.fps.append(info)

    for key, rec in table.records.items():
        _finalize_record(rec, lexicon)
        prefix, num, suffix = key
        if prefix == "":
            table.numeric.add(num)
        else:
            table.prefixes.setdefault(prefix, set()).add(num)
            if suffix:
                table.suffixed.setdefault((prefix, num), set()).add(suffix)
    return table


def _prefix_regex(prefix: str, table: ConstructTable, lexicon: Lexicon) -> "re.Pattern[str]":
    """Regex for '<prefix><digits>[-_.<suffix>]' not preceded by a letter (cached)."""
    rx = table._prefix_rx.get(prefix)
    if rx is None:
        sep = r"[-_ ]?" if len(prefix) >= lexicon.sep_prefix_min_len else ""
        rx = re.compile(r"(?<![a-z])" + re.escape(prefix) + sep +
                        r"(\d+)(?!\d)(?:[-_.](\d+(?:\.\d+)?)(?![\d]))?")
        table._prefix_rx[prefix] = rx
    return rx


# --------------------------------------------------------------------------- #
# Name scanning
# --------------------------------------------------------------------------- #
@dataclass
class IdHit:
    """A plasmid ID found in a name."""

    key: IdKey
    text: str
    span: Tuple[int, int]


def find_ids_in_name(name: str, table: ConstructTable, lexicon: Lexicon
                     ) -> Tuple[List[IdHit], List[str]]:
    """Find known plasmid IDs and unknown ID-like tokens in a (lowercased) name.

    Returns (hits, unknown_tokens). Prefixed IDs: known prefix not preceded by a
    letter, then digits with any zero padding, not followed by a digit. Numeric
    IDs: isolated digit runs of numeric_min..numeric_max digits present in table.
    """
    low = name.lower()
    hits: List[IdHit] = []
    unknown: List[str] = []
    spans: List[Tuple[int, int]] = []
    for prefix in sorted(table.prefixes, key=len, reverse=True):
        for m in _prefix_regex(prefix, table, lexicon).finditer(low):
            if any(m.start() < e and s < m.end() for s, e in spans):
                continue
            num = int(m.group(1))
            suffix = m.group(2) or ""
            variants = table.suffixed.get((prefix, num), set())
            key: Optional[IdKey] = None
            end = m.end() if suffix else m.end(1)
            if suffix and suffix in variants:
                key = (prefix, num, suffix)
            elif (prefix, num, "") in table.records:
                key = (prefix, num, "")
                end = m.end(1)
            if key is not None:
                hits.append(IdHit(key, low[m.start():end], (m.start(), end)))
                spans.append((m.start(), end))
            elif variants:
                unknown.append(low[m.start():m.end()] + "(ambiguous:" +
                               "/".join(sorted(variants)) + ")")
                spans.append((m.start(), m.end(1)))
            elif len(prefix) >= lexicon.unknown_min_prefix_len and \
                    prefix not in lexicon.unknown_ignore_prefixes:
                unknown.append(low[m.start():m.end(1)])
                spans.append((m.start(), m.end(1)))
    if table.numeric:
        num_rx = re.compile(r"(?<!\d)(\d{%d,%d})(?!\d)" % (lexicon.numeric_min, lexicon.numeric_max))
        for m in num_rx.finditer(low):
            if any(m.start() < e and s < m.end() for s, e in spans):
                continue
            n = int(m.group(1))
            if n in table.numeric:
                hits.append(IdHit(("", n, ""), m.group(1), m.span()))
    unique: Dict[IdKey, IdHit] = {}
    for h in hits:
        unique.setdefault(h.key, h)
    return list(unique.values()), list(dict.fromkeys(unknown))


@dataclass
class KeywordHit:
    """A keyword or common-marker match in a name."""

    entry: KeywordEntry
    text: str
    span: Tuple[int, int]
    source: str  # the (lowercased) string the span refers to


def find_keywords(name: str, lexicon: Lexicon) -> Tuple[List[KeywordHit], List[KeywordHit]]:
    """Return (direct_mito_keyword_hits, common_marker_hits) for one name.

    A direct keyword inside a common-marker match (e.g. 'mito' inside
    'mitotracker' or 'mito-gcamp') is suppressed: the marker rule wins.
    """
    low = name.lower()
    common: List[KeywordHit] = []
    for entry in lexicon.common_markers:
        for rx in entry.regexes:
            for m in rx.finditer(low):
                common.append(KeywordHit(entry, m.group(0), m.span(), low))
    direct: List[KeywordHit] = []
    for entry in lexicon.direct_keywords:
        for rx in entry.regexes:
            for m in rx.finditer(low):
                if any(m.start() < c.span[1] and c.span[0] < m.end() for c in common):
                    continue
                if any(m.start() < d.span[1] and d.span[0] < m.end() for d in direct):
                    continue
                direct.append(KeywordHit(entry, m.group(0), m.span(), low))
    # 'mt' inside a longer marker such as 'mt-keima' is redundant; keep the longer one.
    kept: List[KeywordHit] = []
    for c in sorted(common, key=lambda h: h.span[0] - h.span[1]):
        if not any(k.span[0] <= c.span[0] and c.span[1] <= k.span[1] for k in kept):
            kept.append(c)
    kept.sort(key=lambda h: h.span)
    return direct, kept


# --------------------------------------------------------------------------- #
# Overrides
# --------------------------------------------------------------------------- #
@dataclass
class Override:
    """One manual override row: fnmatch pattern -> label."""

    pattern: str
    label: str
    note: str = ""


def load_overrides(path: Optional[os.PathLike]) -> List[Override]:
    """Read an overrides CSV with columns pattern,label,note ('#' lines ignored)."""
    if not path:
        return []
    p = Path(path)
    if not p.is_file():
        raise ClassifierInputError(f"overrides file not found: {p}")
    with open(p, "r", encoding="utf-8-sig", newline="") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    if not reader.fieldnames or "pattern" not in [f.strip().lower() for f in reader.fieldnames]:
        raise ClassifierInputError(f"overrides file {p} needs a header 'pattern,label,note'")
    out: List[Override] = []
    for i, row in enumerate(reader, start=2):
        row = {(k or "").strip().lower(): (v or "").strip() for k, v in row.items()}
        pattern, label = row.get("pattern", ""), row.get("label", "").upper()
        if not pattern:
            continue
        if label not in LABELS:
            raise ClassifierInputError(
                f"overrides file {p}, data row {i - 1}: label {label!r} is not one of {LABELS}")
        out.append(Override(pattern, label, row.get("note", "")))
    return out


def match_override(rel_path: str, overrides: Sequence[Override]) -> Optional[Override]:
    """First override whose pattern matches the relative path or the bare file name."""
    rp = rel_path.replace("\\", "/").lower()
    base = rp.rsplit("/", 1)[-1]
    for ov in overrides:
        pat = ov.pattern.replace("\\", "/").lower()
        if fnmatch.fnmatchcase(rp, pat) or fnmatch.fnmatchcase(base, pat):
            return ov
    return None


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
@dataclass
class ConstructCall:
    """Classification result for one file."""

    rel_path: str
    label: str
    review_priority: str = ""
    mito_fp: str = ""
    fp_class: str = "unknown"
    mito_color: str = "unknown"
    red_mito_flag: bool = False
    matched_ids: str = ""
    evidence: str = ""
    reasons: str = ""

    def as_row(self) -> Dict[str, str]:
        """Dict for csv.DictWriter."""
        return {
            "rel_path": self.rel_path, "label": self.label,
            "review_priority": self.review_priority, "mito_fp": self.mito_fp,
            "fp_class": self.fp_class, "mito_color": self.mito_color,
            "red_mito_flag": "TRUE" if self.red_mito_flag else "FALSE",
            "matched_ids": self.matched_ids, "evidence": self.evidence, "reasons": self.reasons,
        }


@dataclass
class _MitoEvidence:
    source: str          # 'plasmid' or 'keyword'
    label: str           # e.g. JR19 or mls
    fp: Optional[FPInfo]
    is_mls: bool = False


def _strip_extension(rel_path: str, lexicon: Lexicon) -> str:
    """Drop a known data extension so '.tif' / '.nd2' never feed the scanners."""
    low = rel_path.lower()
    for ext in sorted(lexicon.file_extensions, key=len, reverse=True):
        if low.endswith(ext):
            return rel_path[: len(rel_path) - len(ext)]
    return rel_path


def classify_name(rel_path: str, table: Optional[ConstructTable], lexicon: Lexicon,
                  overrides: Optional[Sequence[Override]] = None,
                  extra_names: Iterable[str] = ()) -> ConstructCall:
    """Classify one file from its relative path (folders included) and extra names.

    ``extra_names`` are additional strings to scan, e.g. SlideBook series names.
    """
    rel_norm = str(rel_path).replace("\\", "/")
    table = table or ConstructTable()
    texts = [_strip_extension(rel_norm, lexicon)] + [str(n) for n in extra_names if n]

    id_hits: Dict[IdKey, IdHit] = {}
    unknown_ids: List[str] = []
    direct: List[KeywordHit] = []
    common: List[KeywordHit] = []
    for text in texts:
        hits, unk = find_ids_in_name(text, table, lexicon)
        for h in hits:
            id_hits.setdefault(h.key, h)
        unknown_ids.extend(u for u in unk if u not in unknown_ids)
        d, c = find_keywords(text, lexicon)
        direct.extend(d)
        common.extend(c)

    records = [table.records[k] for k in id_hits]
    evidence: List[str] = []
    reasons: List[str] = []
    mito_ev: List[_MitoEvidence] = []
    nofp_mito: List[PlasmidRecord] = []
    nonmito: List[PlasmidRecord] = []

    # Direct keywords in the name, with the FP glued to them if any.
    seen_kw: Set[str] = set()
    for hit in direct:
        fp = fp_adjacent_to(hit.source, hit.span[0], hit.span[1], lexicon)
        mito_ev.append(_MitoEvidence("keyword", hit.entry.token, fp, hit.entry.token == "mls"))
        if hit.entry.evidence not in seen_kw:
            seen_kw.add(hit.entry.evidence)
            evidence.append(hit.entry.evidence)

    # Plasmid IDs.
    matched: List[str] = []
    for rec in records:
        tag = lexicon.numeric_prefix if rec.key[0] == "" else "plasmid"
        matched.append(f"{rec.display_id}={rec.display_name}")
        if rec.is_mito and rec.mito_fp is not None:
            evidence.append(f"{tag}:{rec.display_id}(mito,{rec.mito_fp.canonical})")
            mito_ev.append(_MitoEvidence("plasmid", rec.display_id, rec.mito_fp))
        elif rec.is_mito:
            evidence.append(f"{tag}:{rec.display_id}(mito,no_FP)")
            nofp_mito.append(rec)
        else:
            evidence.append(f"{tag}:{rec.display_id}({rec.localization_label})")
            nonmito.append(rec)

    # Common markers: 'explained' when a matched plasmid's own name contains it.
    unexplained: List[str] = []
    for hit in common:
        tok = hit.entry.token
        explained = any(rx.search(rec.display_name.lower()) for rec in records
                        for rx in hit.entry.regexes)
        ev = f"common_marker:{tok}" + ("(explained_by_plasmid)" if explained else "")
        if ev not in evidence:
            evidence.append(ev)
        if not explained and tok not in unexplained:
            unexplained.append(tok)
    for u in unknown_ids:
        evidence.append(f"unknown_id:{u}")

    def is_const(e: _MitoEvidence) -> bool:
        return e.fp is not None and e.fp.cls in lexicon.constitutive_classes

    def is_photo(e: _MitoEvidence) -> bool:
        return e.fp is not None and e.fp.cls in lexicon.photo_classes

    const_ev = [e for e in mito_ev if is_const(e)]
    photo_ev = [e for e in mito_ev if is_photo(e)]
    unknown_ev = [e for e in mito_ev if not is_const(e) and not is_photo(e)]

    call = ConstructCall(rel_path=rel_norm, label="REVIEW", review_priority="normal")
    primary: Optional[_MitoEvidence] = None

    def describe(e: _MitoEvidence) -> str:
        fp = e.fp.canonical if e.fp else "unknown FP"
        return f"{'plasmid ' if e.source == 'plasmid' else 'keyword '}{e.label} ({fp})"

    if const_ev:
        call.label, call.review_priority = "MITO", ""
        primary = sorted(const_ev, key=lambda e: e.source != "plasmid")[0]
        reasons.append("Mito-targeted fluorescence: " + ", ".join(describe(e) for e in const_ev) + ".")
        if photo_ev:
            reasons.append("Also photoactivatable/convertible mito evidence ("
                           + ", ".join(describe(e) for e in photo_ev)
                           + "); constitutive evidence wins.")
    elif photo_ev:
        call.label, call.review_priority = "MITO_PHOTOACTIVATABLE", ""
        primary = sorted(photo_ev, key=lambda e: e.source != "plasmid")[0]
        reasons.append("Only photoactivatable/photoconvertible mito FP: "
                       + ", ".join(describe(e) for e in photo_ev)
                       + ". Mito signal may be dark or change color; excluded from the default list.")
    elif unknown_ev:
        weak = [e for e in unknown_ev if e.source == "keyword" and not e.is_mls and e.fp is None]
        if weak and len(weak) == len(unknown_ev) and nonmito:
            call.label, call.review_priority = "REVIEW", "normal"
            reasons.append("Conflict: mito keyword " + ", ".join(e.label for e in weak)
                           + " without an attached FP, but recognized plasmids are not mito ("
                           + ", ".join(r.display_id for r in nonmito) + ").")
        else:
            call.label, call.review_priority = "MITO", ""
            primary = unknown_ev[0]
            reasons.append("Mito evidence without a recognized FP: "
                           + ", ".join(describe(e) for e in unknown_ev) + ".")
    elif unexplained:
        call.review_priority = "high"
        reasons.append("Only common mito marker(s) " + ", ".join(unexplained)
                       + " (dye/reporter/protein not verified by the plasmid table); check by eye.")
    elif nofp_mito:
        call.review_priority = "high"
        reasons.append("Mito-localized construct(s) without a fluorescent protein ("
                       + ", ".join(r.display_id for r in nofp_mito)
                       + "); the mito signal must come from something else.")
    elif unknown_ids:
        reasons.append("ID-like token(s) not in the plasmid table: " + ", ".join(unknown_ids) + ".")
    elif nonmito:
        call.label, call.review_priority = "NOT_MITO", ""
        reasons.append("Recognized plasmid(s) are not mito-localized: "
                       + ", ".join(f"{r.display_id} ({r.localization_label})" for r in nonmito) + ".")
    else:
        reasons.append("No mito keyword, plasmid ID or marker recognized in the name.")

    # Secondary notes for positive calls.
    if call.label.startswith("MITO"):
        if unexplained:
            reasons.append("Also mentions common marker(s): " + ", ".join(unexplained) + ".")
        if unknown_ids:
            reasons.append("Unrecognized ID-like token(s): " + ", ".join(unknown_ids) + ".")
        if nonmito:
            reasons.append("Non-mito plasmid(s) also present: "
                           + ", ".join(r.display_id for r in nonmito) + ".")
        if any(e.source == "keyword" and e.fp is None for e in mito_ev) and primary \
                and primary.fp is not None and call.label == "MITO_PHOTOACTIVATABLE":
            reasons.append("Keyword evidence without FP is consistent with the plasmid.")

    if primary is not None and primary.fp is not None:
        call.mito_fp = primary.fp.canonical
        call.fp_class = primary.fp.cls
        call.mito_color = primary.fp.color
        if primary.fp.note:
            reasons.append(f"{primary.fp.canonical}: {primary.fp.note}")
        colors = {e.fp.color for e in (const_ev if call.label == "MITO" else photo_ev) if e.fp}
        if len(colors) > 1:
            reasons.append("Mito FPs in different channels (" + ", ".join(sorted(colors)) + ").")
    call.red_mito_flag = call.mito_color == "red"
    if call.red_mito_flag:
        reasons.append("RED mito channel: the pipeline assumes channel 1 = mito; check channel order.")

    # Manual override (first match wins) replaces the label but keeps the FP info.
    ov = match_override(rel_norm, overrides or [])
    if ov is not None:
        heuristic = call.label
        call.label = ov.label
        call.review_priority = "normal" if ov.label == "REVIEW" else ""
        evidence.insert(0, f"override:{ov.pattern}")
        reasons.insert(0, f"Manual override '{ov.pattern}' -> {ov.label}"
                       + (f" ({ov.note})" if ov.note else "")
                       + f"; heuristic would say {heuristic}.")

    call.matched_ids = "; ".join(matched)
    call.evidence = ";".join(evidence)
    call.reasons = " ".join(reasons)
    return call


# --------------------------------------------------------------------------- #
# Directory scanning and reports
# --------------------------------------------------------------------------- #
def _is_candidate(name: str, lexicon: Lexicon) -> bool:
    low = name.lower()
    if any(fnmatch.fnmatchcase(low, p.lower()) for p in lexicon.skip_name_patterns):
        return False
    return any(low.endswith(ext) for ext in lexicon.file_extensions)


def scan_directory(root: os.PathLike, lexicon: Lexicon) -> Iterator[str]:
    """Yield relative (forward-slash) paths of candidate movie files under root.

    Only names are inspected; files are never opened. Lock/temp/hidden names
    are skipped. A directory with a data extension (SlideBook 7 '.sldy') is
    yielded as one item and not descended into.
    """
    root_path = Path(root)
    if not root_path.is_dir():
        raise ClassifierInputError(f"scan directory not found: {root_path}")
    for dirpath, dirnames, filenames in os.walk(root_path):
        dirnames.sort()
        here = Path(dirpath)
        keep = []
        for d in dirnames:
            if _is_candidate(d, lexicon):
                yield (here / d).relative_to(root_path).as_posix()
            elif not d.startswith("."):
                keep.append(d)
        dirnames[:] = keep
        for f in sorted(filenames):
            if _is_candidate(f, lexicon):
                yield (here / f).relative_to(root_path).as_posix()


def write_reports(calls: Sequence[ConstructCall], out_dir: os.PathLike) -> Tuple[Path, Path]:
    """Write mito_classification.csv (all) and mito_review.csv (REVIEW, high first)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    all_path = out / CLASSIFICATION_CSV
    review_path = out / REVIEW_CSV
    with open(all_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, lineterminator="\n")
        w.writeheader()
        for c in calls:
            w.writerow(c.as_row())
    review = [c for c in calls if c.label == "REVIEW"]
    review.sort(key=lambda c: (c.review_priority != "high", c.rel_path.lower()))
    with open(review_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, lineterminator="\n")
        w.writeheader()
        for c in review:
            w.writerow(c.as_row())
    return all_path, review_path


# --------------------------------------------------------------------------- #
# Synthetic test tree
# --------------------------------------------------------------------------- #
_FILLERS = ["05feb2025", "20250205", "2024-11-03", "dec21", "mar2025", "Feb26", "cell3",
            "fish2", "embryo1", "MAX", "proj", "45TL", "16bit", "t1L", "0-49", "z30",
            "series12", "small", "walls", "acta2", "smtx", "mtsat", "rgc", "stack",
            "timelapse", "dpf3", "2x5ms"]
_FOLDERS = ["2025-02-05", "exp A", "fish 1", "day3", "raw data", "sessionB", "2024_dec"]
_SYN_EXTS = [".sld", ".tif", ".tiff", ".ome.tif"]
_SEPARATORS = ["_", "-", " "]


def _is_neutral(token: str, table: ConstructTable, lexicon: Lexicon) -> bool:
    """True when a filler token alone triggers no ID, keyword or marker."""
    hits, unknown = find_ids_in_name(token, table, lexicon)
    direct, common = find_keywords(token, lexicon)
    return not (hits or unknown or direct or common)


def _rand_case(rng: random.Random, s: str) -> str:
    return rng.choice([s.lower(), s.upper(), s, s[:1].upper() + s[1:].lower()])


def _format_id(rng: random.Random, rec: PlasmidRecord) -> str:
    """Render a table ID with random zero padding and case (suffix kept)."""
    prefix, num, suffix = rec.key
    if prefix == "":
        return str(num)
    width = rng.randint(len(str(num)), max(len(str(num)), 4))
    text = prefix + str(num).zfill(width)
    if suffix:
        sep = "-" if "-" in rec.display_id else "_"
        text += sep + suffix
    return _rand_case(rng, text)


def make_synthetic_tree(root: os.PathLike, seed: int = 0, table: Optional[ConstructTable] = None,
                        lexicon: Optional[Lexicon] = None, per_recipe: int = 6
                        ) -> List[Tuple[str, str]]:
    """Create a nested tree of EMPTY movie files whose names mix tags, IDs and decoys.

    Returns [(rel_path, expected_label), ...] and writes expected_labels.csv in
    root. A few non-candidate files (lock file, hidden file, notes.txt) are also
    written; they must be skipped by scan_directory.
    """
    if table is None or not table.records:
        raise ClassifierInputError("make_synthetic_tree needs a loaded plasmid table")
    lexicon = lexicon or load_lexicon()
    rng = random.Random(seed)
    root_path = Path(root)
    root_path.mkdir(parents=True, exist_ok=True)

    recs = [r for r in table.records.values() if r.key[0] in table.prefixes or r.key[0] == ""]
    const_ids = [r for r in recs if r.is_mito and r.mito_fp and r.mito_fp.cls in lexicon.constitutive_classes]
    red_ids = [r for r in const_ids if r.mito_fp.color == "red"]
    photo_ids = [r for r in recs if r.is_mito and r.mito_fp and r.mito_fp.cls in lexicon.photo_classes]
    nonmito_ids = [r for r in recs if not r.is_mito]
    fillers = [f for f in _FILLERS if _is_neutral(f, table, lexicon)]
    folders = [f for f in _FOLDERS if _is_neutral(f, table, lexicon)]
    # Common-marker aliases (from the raw lexicon) used verbatim as name tokens.
    markers = [a for item in lexicon.raw["common_markers"]
               for a in item.get("aliases", [item["token"]]) if " " not in a]
    used_prefixes = set(table.prefixes)
    uid_letter = next((c for c in "qvwzyr" if c not in used_prefixes), "q")

    # Unknown IDs: a known multi-letter prefix with a number absent from the table.
    unknown_ids: List[str] = []
    for prefix in sorted(table.prefixes):
        if len(prefix) < lexicon.unknown_min_prefix_len or prefix in lexicon.unknown_ignore_prefixes:
            continue
        for n in range(90, 1000):
            if n not in table.prefixes[prefix] and (prefix, n) not in table.suffixed:
                unknown_ids.append(prefix + str(n))
                break

    def pick(pool: List[PlasmidRecord]) -> Optional[PlasmidRecord]:
        return rng.choice(pool) if pool else None

    def mls_token() -> str:
        fp = rng.choice(["EGFP", "mCherry", "GFP", ""])
        tok = rng.choice(["mls", "MLS", "Mls"]) + rng.choice(["", "-", "_"]) * bool(fp) + fp
        return rng.choice(["", "MAX"]) + tok

    # Each recipe returns (list of (kind, token), expected label); kind 'id' allows glueing.
    recipes = {
        "mls_tag": lambda: ([("kw", mls_token())], "MITO"),
        "mito_constitutive_id": lambda: ([("id", _format_id(rng, pick(const_ids)))], "MITO"),
        "red_mito_id": lambda: ([("id", _format_id(rng, pick(red_ids)))], "MITO"),
        "photo_id": lambda: ([("id", _format_id(rng, pick(photo_ids)))], "MITO_PHOTOACTIVATABLE"),
        "photo_plus_nonmito": lambda: ([("id", _format_id(rng, pick(photo_ids))),
                                        ("id", _format_id(rng, pick(nonmito_ids)))],
                                       "MITO_PHOTOACTIVATABLE"),
        "photo_plus_constitutive": lambda: ([("id", _format_id(rng, pick(photo_ids))),
                                             ("id", _format_id(rng, pick(const_ids)))], "MITO"),
        "nonmito_plus_constitutive": lambda: ([("id", _format_id(rng, pick(nonmito_ids))),
                                               ("id", _format_id(rng, pick(const_ids)))], "MITO"),
        "nonmito_ids": lambda: ([("id", _format_id(rng, r)) for r in
                                 rng.sample(nonmito_ids, min(len(nonmito_ids), rng.randint(1, 2)))],
                                "NOT_MITO"),
        "common_marker_only": lambda: ([("kw", _rand_case(rng, rng.choice(markers)))], "REVIEW"),
        "unknown_id": lambda: ([("id", _rand_case(rng, rng.choice(unknown_ids)))], "REVIEW"),
        "nothing": lambda: ([], "REVIEW"),
    }
    required = {"mito_constitutive_id": const_ids, "red_mito_id": red_ids, "photo_id": photo_ids,
                "photo_plus_nonmito": photo_ids and nonmito_ids,
                "photo_plus_constitutive": photo_ids and const_ids,
                "nonmito_plus_constitutive": nonmito_ids and const_ids,
                "nonmito_ids": nonmito_ids, "common_marker_only": markers, "unknown_id": unknown_ids}

    results: List[Tuple[str, str, str]] = []
    counter = 0
    for recipe, make in recipes.items():
        if recipe in required and not required[recipe]:
            continue
        for _ in range(per_recipe):
            counter += 1
            parts, expected = make()
            for _ in range(rng.randint(1, 3)):
                parts.insert(rng.randint(0, len(parts)), ("filler", _rand_case(rng, rng.choice(fillers))))
            parts.append(("filler", f"{uid_letter}{counter:03d}"))
            name = parts[0][1]
            for (pkind, ptext), (kind, text) in zip(parts, parts[1:]):
                glue_ok = pkind == "id" and ptext[-1:].isdigit() and text[:1].isalpha() and kind == "id"
                sep = rng.choice(_SEPARATORS + [""] * 2) if glue_ok else rng.choice(_SEPARATORS)
                name += sep + text
            depth = rng.randint(0, 2)
            sub = [rng.choice(folders) for _ in range(depth)]
            rel = "/".join(sub + [name + rng.choice(_SYN_EXTS)])
            target = root_path.joinpath(*rel.split("/"))
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b"")
            results.append((rel, expected, recipe))

    # Non-candidates that scan_directory must ignore.
    for junk in ("~$lockfile.sld", ".hidden_mls.tif", "notes.txt"):
        (root_path / junk).write_bytes(b"")
    with open(root_path / EXPECTED_LABELS_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["rel_path", "expected_label", "recipe"])
        w.writerows(results)
    return [(r, e) for r, e, _ in results]


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _read_names_file(path: Path) -> List[str]:
    if not path.is_file():
        raise ClassifierInputError(f"names file not found: {path}")
    with open(path, "r", encoding="utf-8-sig") as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]


def _summary(calls: Sequence[ConstructCall]) -> str:
    counts = {lab: 0 for lab in LABELS}
    for c in calls:
        counts[c.label] = counts.get(c.label, 0) + 1
    high = sum(1 for c in calls if c.label == "REVIEW" and c.review_priority == "high")
    red = sum(1 for c in calls if c.red_mito_flag)
    lines = [f"Classified {len(calls)} file(s):"]
    lines += [f"  {lab:<22}{counts[lab]:>6}" + (f"   ({high} high priority)" if lab == "REVIEW" else "")
              for lab in LABELS]
    lines.append(f"  red-mito (check channel order): {red}")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Filename-only classifier: which time-lapse files image mitochondrial fluorescence.",
        epilog="Labels: MITO, MITO_PHOTOACTIVATABLE, NOT_MITO, REVIEW. Exit 0 ok, 2 bad input.")
    p.add_argument("--plasmid-csv", default=os.environ.get("MITO_PLASMID_CSV"),
                   help="lab plasmid table CSV (default: $MITO_PLASMID_CSV)")
    p.add_argument("--lexicon", default=None,
                   help=f"lexicon JSON (default: {DEFAULT_LEXICON_PATH})")
    p.add_argument("--overrides", default=None, help="overrides CSV with columns pattern,label,note")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--scan-dir", metavar="DIR", help="folder to scan recursively (names only)")
    src.add_argument("--names-file", metavar="FILE", help="text file, one file name per line")
    src.add_argument("--make-synthetic", metavar="DIR",
                     help="write a synthetic tree of empty files + expected_labels.csv")
    p.add_argument("--out", metavar="DIR", help="output folder for the CSV reports")
    p.add_argument("--seed", type=int, default=0, help="random seed for --make-synthetic")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    args = build_parser().parse_args(argv)
    try:
        lexicon = load_lexicon(args.lexicon)
        table: Optional[ConstructTable] = None
        if args.plasmid_csv:
            table = load_construct_table(args.plasmid_csv, lexicon)
            print(f"Plasmid table: {table.source} ({len(table.records)} unique IDs)")
            for q in table.quirks:
                print(f"  note: {q}")
        if args.make_synthetic:
            if table is None:
                raise ClassifierInputError("--make-synthetic needs --plasmid-csv (or $MITO_PLASMID_CSV)")
            made = make_synthetic_tree(args.make_synthetic, seed=args.seed, table=table, lexicon=lexicon)
            print(f"Wrote {len(made)} synthetic files + {EXPECTED_LABELS_CSV} under {args.make_synthetic}")
            return 0
        if not args.out:
            raise ClassifierInputError("--out is required with --scan-dir / --names-file")
        if table is None:
            print("WARNING: no plasmid CSV given; only keywords and markers are used.", file=sys.stderr)
        overrides = load_overrides(args.overrides)
        names = list(scan_directory(args.scan_dir, lexicon)) if args.scan_dir \
            else _read_names_file(Path(args.names_file))
        calls = [classify_name(n, table, lexicon, overrides) for n in names]
        all_path, review_path = write_reports(calls, args.out)
    except ClassifierInputError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(_summary(calls))
    print(f"Wrote {all_path}\nWrote {review_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
