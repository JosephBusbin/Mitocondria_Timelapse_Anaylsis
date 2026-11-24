#!/usr/bin/env python3
"""Unit tests for scripts/construct_classifier.py (standard library only).

Run (Linux):    module load python/3.12.7 && python3 scripts/test_construct_classifier.py -v
Run (Windows):  python scripts\\test_construct_classifier.py -v

Tests use an inline mini plasmid table written to a temporary folder, and
additionally the real lab table when it exists (path from $MITO_PLASMID_CSV or
the default SciClone location). Temporary files live in TemporaryDirectory.
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import construct_classifier as cc  # noqa: E402

CLI = SCRIPT_DIR / "construct_classifier.py"
REAL_CSV = Path(os.environ.get(
    "MITO_PLASMID_CSV",
    "/sciclone/data10/jbusbin/bestman_lab_data/Plasmid_Primer_List_Bestman_Lab.csv"))

HEADER = ["s", "sequenced? (date)", "does this plasmid make FP positive cells?", "glycerol stock?",
          "-20 DNA duplicate?", "addgene submission?", "plasmid ID", "name (Prom_GOI_3')",
          "gateway compnent?", "localization sequence?", "GOI? tag?", "Fluorescent Protein 1 (ME)",
          "Fluorescent Protein 2 (3')", "p2a", "uas, MCS, promoter?", "antibiotic",
          "plasmid description/elements", "full map?", "long term use", "Forward Primer",
          "Reverse Primer", "product size "]

# (plasmid ID, name, localization, FP1, FP2) -- a slice of the real table plus test rows.
MINI_ROWS = [
    ("OBF001", "psox2bd::mitopaGFP_p2a-tagRFPt", "mito", "paGFP", "tagRFPT"),
    ("avc01", "p4XUAS_tagRFPt_p2a-mRuby2", "-", "tagRFPt (middle)", "ruby2"),
    ("JEB009", "pSox2bd_tagRFPt_p2a-tagRFPt", "no", "tagRFPT (middle)", "tagRFPT (3')"),
    ("JR19", "pUAS-mlsEGFPpA", "mito", "EGFP", ""),
    ("MF002", "pCMV_EGFP_p2a-mitoTAGRFPT", "mito", "EGFP", "tagRFPt"),
    ("MF008", "p4XUAS_EGFP_p2a-mit0TagRFPT", "mito", "EGFP", "tagRFPt"),
    ("54282", "mEmeraldTOMM20N10", "mito", "mEmerald", ""),
    ("55796", "mitodendra2", "mito", "dendra2", ""),
    ("45161", "pcDNA3_Drp1K38A", "mito", "", ""),
    ("jeb026", "pME_drp1k38a", "no", "", ""),
    ("MON22-3", "pME_tagBFP2-Farn_(stop)_p2a", "mem", "tagBFP2", ""),
    ("MON22-4", "pSox2bd_tagBFP2-farn-p2a_polyA", "mem", "tagBFP2", ""),
    ("MON23-3.2", "pSoxbd_gng4_p2a-eGFP", "", "EGFP", ""),
    ("p117", "psox2_nlsturboRFP (good one is \"4\")", "nuc", "turboRFP", ""),
    ("p117", "sox2bd_gal4uas_nlsturboRFP", "nuc", "turboRFP", ""),
    ("55295", "mTagBFP2-Farnesyl-5", "mem", "mTagBFP2\n                        \n   ", ""),
    ("250205", "test row: 6-digit numeric ID that hides inside the date 20250205", "mito", "EGFP", ""),
    ("t2k233", "pME_nlsmCherry", "nucleus", "mCherry", ""),
]


def write_mini_csv(path: Path, bom: bool = False) -> Path:
    """Write the mini table with the real 22-column header (optionally with a BOM)."""
    with open(path, "w", encoding="utf-8-sig" if bom else "utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        for i, (pid, name, loc, fp1, fp2) in enumerate(MINI_ROWS):
            row = [""] * len(HEADER)
            row[0], row[6], row[7], row[9], row[11], row[12] = str(i), pid, name, loc, fp1, fp2
            row[16] = "multi-line\ndescription cell"
            w.writerow(row)
    return path


def run_cli(*args: str) -> subprocess.CompletedProcess:
    """Run the classifier CLI in a subprocess with the current interpreter."""
    env = dict(os.environ)
    env.pop("MITO_PLASMID_CSV", None)
    return subprocess.run([sys.executable, str(CLI), *args], capture_output=True, text=True, env=env)


def read_csv(path: Path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def confusion_text(pairs):
    """Expected-vs-observed table for failure messages."""
    counts = Counter(pairs)
    labels = list(cc.LABELS)
    lines = ["expected \\ observed".ljust(24) + "".join(l[:12].rjust(14) for l in labels)]
    for e in labels:
        lines.append(e.ljust(24) + "".join(str(counts.get((e, o), 0)).rjust(14) for o in labels))
    return "\n".join(lines)


class MiniTableTests(unittest.TestCase):
    """Behaviour on the inline mini table."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.tmp_path = Path(cls.tmp.name)
        cls.csv_path = write_mini_csv(cls.tmp_path / "mini.csv")
        cls.lex = cc.load_lexicon()
        cls.table = cc.load_construct_table(cls.csv_path, cls.lex)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def classify(self, name, overrides=None, extra=()):
        return cc.classify_name(name, self.table, self.lex, overrides, extra)

    # --- required real-name cases -------------------------------------------
    def test_obf01_photoactivatable(self):
        c = self.classify("05feb2025obf01avc01-B.sld")
        self.assertEqual(c.label, "MITO_PHOTOACTIVATABLE")
        self.assertEqual((c.mito_fp, c.fp_class, c.mito_color), ("paGFP", "photoactivatable", "green"))
        self.assertIn("plasmid:OBF001(mito,paGFP)", c.evidence)
        self.assertIn("plasmid:avc01(-)", c.evidence)

    def test_jeb09_jr19_mito_egfp(self):
        c = self.classify("max_d318dec21jeb09jr19t1L_0-49.tif")
        self.assertEqual(c.label, "MITO")
        self.assertEqual((c.mito_fp, c.mito_color), ("EGFP", "green"))
        self.assertIn("plasmid:JR19(mito,EGFP)", c.evidence)
        self.assertIn("plasmid:JEB009(no)", c.evidence)
        self.assertIn("JR19=pUAS-mlsEGFPpA", c.matched_ids)

    def test_mls_tag(self):
        c = self.classify("MAXmlsEGFP_T7LFeb26_45TL_16bit.tif")
        self.assertEqual(c.label, "MITO")
        self.assertEqual(c.mito_fp, "EGFP")
        self.assertIn("mls_tag", c.evidence)

    def test_red_mito_mf002(self):
        for name in ("cell_mf02_x.sld", "MF002.tif", "mf0002 2025.tif"):
            c = self.classify(name)
            self.assertEqual(c.label, "MITO", name)
            self.assertEqual((c.mito_fp, c.mito_color), ("tagRFPT", "red"), name)
            self.assertTrue(c.red_mito_flag, name)

    def test_misspelled_mit0(self):
        c = self.classify("mf8_run.sld")
        self.assertEqual((c.label, c.mito_fp), ("MITO", "tagRFPT"))

    def test_fp_before_marker(self):
        c = self.classify("54282_cell.sld")
        self.assertEqual((c.label, c.mito_fp), ("MITO", "mEmerald"))

    def test_jeb09_only_not_mito(self):
        self.assertEqual(self.classify("jeb09.sld").label, "NOT_MITO")

    def test_nothing_review_normal(self):
        c = self.classify("random_cells_2024.sld")
        self.assertEqual((c.label, c.review_priority), ("REVIEW", "normal"))

    def test_common_marker_review_high(self):
        c = self.classify("cell3_TMRM_2024.sld")
        self.assertEqual((c.label, c.review_priority), ("REVIEW", "high"))
        self.assertIn("common_marker:tmrm", c.evidence)

    def test_drp1_never_mito(self):
        c = self.classify("drp1k38a_only.sld")
        self.assertEqual((c.label, c.review_priority), ("REVIEW", "high"))
        # Explained by a recognized non-mito Drp1 plasmid -> NOT_MITO.
        c2 = self.classify("jeb26_drp1k38a.sld")
        self.assertEqual(c2.label, "NOT_MITO")
        # Mito-localized row with no FP (45161 Drp1K38A) is not fluorescence evidence.
        c3 = self.classify("45161_fish.sld")
        self.assertEqual((c3.label, c3.review_priority), ("REVIEW", "high"))

    def test_photoconvertible(self):
        c = self.classify("55796_conv.tif")
        self.assertEqual((c.label, c.fp_class), ("MITO_PHOTOACTIVATABLE", "photoconvertible"))
        self.assertIn("red after", c.reasons)

    def test_constitutive_beats_photo(self):
        c = self.classify("obf01_jr19.sld")
        self.assertEqual(c.label, "MITO")
        self.assertIn("constitutive evidence wins", c.reasons)

    def test_mito_keyword_photo_fp(self):
        self.assertEqual(self.classify("mito-paGFP_act.sld").label, "MITO_PHOTOACTIVATABLE")
        self.assertEqual(self.classify("mito_only.sld").label, "MITO")

    def test_reporter_with_mito_word_goes_to_review(self):
        for name in ("mitotracker_red.sld", "mito-GCaMP6.tif", "cell_mitoTimer.tif"):
            c = self.classify(name)
            self.assertEqual((c.label, c.review_priority), ("REVIEW", "high"), name)

    def test_conflict_keyword_vs_nonmito_plasmid(self):
        c = self.classify("mito_jeb09.sld")
        self.assertEqual(c.label, "REVIEW")
        self.assertIn("Conflict", c.reasons)

    def test_override_beats_heuristics(self):
        ov_path = self.tmp_path / "ov.csv"
        ov_path.write_text("pattern,label,note\n# comment line\n*JEB09*,MITO,lab says mito dye\n"
                           "*.tif,NOT_MITO,catch-all\n", encoding="utf-8")
        ovs = cc.load_overrides(ov_path)
        c = self.classify("sub/dir/jeb09.sld", ovs)
        self.assertEqual(c.label, "MITO")
        self.assertTrue(c.evidence.startswith("override:*JEB09*"))
        self.assertEqual(self.classify("MAXmlsEGFP.tif", ovs).label, "NOT_MITO")
        # Matches the relative path too.
        ovs2 = [cc.Override("exp1/*", "REVIEW")]
        self.assertEqual(self.classify("exp1/mlsEGFP.sld", ovs2).label, "REVIEW")

    def test_override_bad_label(self):
        bad = self.tmp_path / "bad.csv"
        bad.write_text("pattern,label,note\n*x*,MAYBE,\n", encoding="utf-8")
        with self.assertRaises(cc.ClassifierInputError):
            cc.load_overrides(bad)

    # --- ID normalisation ------------------------------------------------------
    def test_normalize_id_padding(self):
        self.assertEqual(cc.normalize_id("obf01"), cc.normalize_id("OBF001"))
        self.assertEqual(cc.normalize_id("obf0001"), cc.normalize_id("OBF001"))
        self.assertEqual(cc.normalize_id("p116?"), ("p", 116, ""))
        self.assertEqual(cc.normalize_id("MON23-3.2"), ("mon", 23, "3.2"))
        self.assertEqual(cc.normalize_id("MON23_6.3"), ("mon", 23, "6.3"))
        self.assertEqual(cc.normalize_id("t2k218"), ("t2k", 218, ""))
        self.assertEqual(cc.normalize_id("23348"), ("", 23348, ""))
        for name in ("x_obf01_y.sld", "x_OBF001_y.sld", "x_obf0001_y.sld"):
            hits, _ = cc.find_ids_in_name(name, self.table, self.lex)
            self.assertEqual([h.key for h in hits], [("obf", 1, "")], name)

    def test_dates_never_ids(self):
        for name in ("feb2025.sld", "dec21.sld", "mar2025.sld", "05feb2025_dec21_mar2025.sld"):
            hits, unknown = cc.find_ids_in_name(name, self.table, self.lex)
            self.assertEqual((hits, unknown), ([], []), name)

    def test_eight_digit_date_not_addgene(self):
        hits, _ = cc.find_ids_in_name("20250205_cell.sld", self.table, self.lex)
        self.assertEqual(hits, [])
        hits, _ = cc.find_ids_in_name("x250205_cell.sld", self.table, self.lex)
        self.assertEqual([h.key for h in hits], [("", 250205, "")])

    def test_suffixed_ids(self):
        self.assertIn("plasmid:MON22-3(mem)", self.classify("mon22-3_a.sld").evidence)
        self.assertIn("MON23-3.2", self.classify("MON23-3.2.tif").evidence)
        c = self.classify("mon22.sld")  # base absent, several suffix variants
        self.assertEqual(c.label, "REVIEW")
        self.assertIn("unknown_id:mon22", c.evidence)

    def test_unknown_id(self):
        c = self.classify("jeb99_cell.sld")
        self.assertEqual((c.label, c.review_priority), ("REVIEW", "normal"))
        self.assertIn("unknown_id:jeb99", c.evidence)

    def test_decoys(self):
        for name in ("acta2_cell.sld", "smtx.sld", "walls.sld", "small.sld", "mtsat.sld",
                     "p2a_construct.sld"):
            c = self.classify(name)
            self.assertEqual(c.evidence, "", name)
            self.assertEqual(c.label, "REVIEW", name)

    def test_letter_boundary_short_markers(self):
        self.assertIn("common_marker:mt", self.classify("mt-x_2025.sld").evidence)
        self.assertIn("common_marker:mcu", self.classify("5mcu3.sld").evidence)
        self.assertNotIn("mcu", self.classify("smcurve.sld").evidence)

    def test_extra_names(self):
        c = self.classify("Series_20.sld", extra=["mlsEGFP photactivate"])
        self.assertEqual(c.label, "MITO")

    # --- table loading -------------------------------------------------------
    def test_duplicate_ids_merged(self):
        rec = self.table.records[("p", 117, "")]
        self.assertEqual(len(rec.table_ids), 2)
        self.assertTrue(any("duplicate" in q for q in self.table.quirks))

    def test_multiline_cell_and_bom(self):
        bom_csv = write_mini_csv(self.tmp_path / "mini_bom.csv", bom=True)
        t = cc.load_construct_table(bom_csv, self.lex)
        self.assertEqual(len(t.records), len(self.table.records))
        self.assertEqual(t.records[("", 55295, "")].fps[0].canonical, "BFP")

    def test_missing_column_error(self):
        bad = self.tmp_path / "bad_header.csv"
        bad.write_text("a,b\n1,2\n", encoding="utf-8")
        with self.assertRaises(cc.ClassifierInputError):
            cc.load_construct_table(bad, self.lex)

    # --- scanning, reports and CLI ---------------------------------------------
    def test_scan_directory_skips_junk(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "a").mkdir()
            for n in ("a/mls1.sld", "a/x.TIF", "~$lock.sld", ".hid.tif", "notes.txt", "b.ome.tif"):
                (root / n).write_bytes(b"")
            (root / "movie.sldy").mkdir()
            (root / "movie.sldy" / "inner.tif").write_bytes(b"")
            got = sorted(cc.scan_directory(root, self.lex))
            self.assertEqual(got, ["a/mls1.sld", "a/x.TIF", "b.ome.tif", "movie.sldy"])

    def test_names_file_cli_and_columns(self):
        with tempfile.TemporaryDirectory() as d:
            names = Path(d) / "names.txt"
            names.write_text("05feb2025obf01avc01-B.sld\n# comment\n\ncell3_TMRM_2024.sld\n"
                             "random_cells_2024.sld\njeb09.sld\n", encoding="utf-8")
            out = Path(d) / "out"
            r = run_cli("--plasmid-csv", str(self.csv_path), "--names-file", str(names), "--out", str(out))
            self.assertEqual(r.returncode, 0, r.stderr)
            rows = read_csv(out / cc.CLASSIFICATION_CSV)
            self.assertEqual(list(rows[0].keys()), cc.REPORT_COLUMNS)
            self.assertEqual([x["label"] for x in rows],
                             ["MITO_PHOTOACTIVATABLE", "REVIEW", "REVIEW", "NOT_MITO"])
            review = read_csv(out / cc.REVIEW_CSV)
            self.assertEqual(list(review[0].keys()), cc.REPORT_COLUMNS)
            self.assertEqual([x["review_priority"] for x in review], ["high", "normal"])

    def test_cli_bad_input_exit_2(self):
        with tempfile.TemporaryDirectory() as d:
            r = run_cli("--plasmid-csv", str(Path(d) / "missing.csv"), "--names-file",
                        str(Path(d) / "n.txt"), "--out", d)
            self.assertEqual(r.returncode, 2)
            r = run_cli("--make-synthetic", d)  # needs a plasmid CSV
            self.assertEqual(r.returncode, 2)

    def _synthetic_roundtrip(self, csv_path: Path, seed: int):
        table = cc.load_construct_table(csv_path, self.lex)
        with tempfile.TemporaryDirectory() as d:
            root, out = Path(d) / "tree", Path(d) / "out"
            expected = dict(cc.make_synthetic_tree(root, seed=seed, table=table, lexicon=self.lex))
            self.assertGreater(len(expected), 50)
            self.assertTrue(any("/" in k for k in expected), "tree should be nested")
            r = run_cli("--plasmid-csv", str(csv_path), "--scan-dir", str(root), "--out", str(out))
            self.assertEqual(r.returncode, 0, r.stderr)
            observed = {row["rel_path"]: row["label"] for row in read_csv(out / cc.CLASSIFICATION_CSV)}
            self.assertEqual(set(observed), set(expected), "scan found a different file set")
            pairs = [(expected[k], observed[k]) for k in expected]
            wrong = [f"  {k}: expected {expected[k]}, observed {observed[k]}"
                     for k in sorted(expected) if expected[k] != observed[k]]
            if wrong:
                self.fail("synthetic mismatches:\n" + "\n".join(wrong) + "\n" + confusion_text(pairs))

    def test_synthetic_tree_cli_scan(self):
        for seed in (1, 2, 3):
            with self.subTest(seed=seed):
                self._synthetic_roundtrip(self.csv_path, seed)

    def test_synthetic_cli_writes_expected(self):
        with tempfile.TemporaryDirectory() as d:
            r = run_cli("--plasmid-csv", str(self.csv_path), "--make-synthetic", d, "--seed", "5")
            self.assertEqual(r.returncode, 0, r.stderr)
            rows = read_csv(Path(d) / cc.EXPECTED_LABELS_CSV)
            self.assertEqual(set(rows[0].keys()), {"rel_path", "expected_label", "recipe"})


@unittest.skipUnless(REAL_CSV.is_file(), f"real plasmid CSV not found at {REAL_CSV}")
class RealTableTests(MiniTableTests.__base__):
    """Same key cases against the real lab table (read-only)."""

    @classmethod
    def setUpClass(cls):
        cls.lex = cc.load_lexicon()
        cls.table = cc.load_construct_table(REAL_CSV, cls.lex)

    def test_real_names(self):
        expect = {
            "05feb2025obf01avc01-B.sld": ("MITO_PHOTOACTIVATABLE", "paGFP", "green"),
            "max_d318dec21jeb09jr19t1L_0-49.tif": ("MITO", "EGFP", "green"),
            "MAXmlsEGFP_T7LFeb26_45TL_16bit.tif": ("MITO", "EGFP", "green"),
            "mf02_cell.sld": ("MITO", "tagRFPT", "red"),
        }
        for name, (label, fp, color) in expect.items():
            c = cc.classify_name(name, self.table, self.lex)
            self.assertEqual((c.label, c.mito_fp, c.mito_color), (label, fp, color), name)
        self.assertEqual(cc.classify_name("jeb09.sld", self.table, self.lex).label, "NOT_MITO")
        self.assertEqual(cc.classify_name("drp1k38a_only.sld", self.table, self.lex).label, "REVIEW")

    def test_real_table_shape(self):
        self.assertGreater(len(self.table.records), 150)
        mito = [r for r in self.table.records.values() if r.is_mito]
        self.assertGreaterEqual(len(mito), 45)
        for pid in ("MF002", "MF005", "MF010"):
            rec = self.table.records[cc.normalize_id(pid)]
            self.assertEqual(rec.mito_fp.canonical, "tagRFPT", pid)

    def test_real_dates_never_ids(self):
        for name in ("05feb2025.sld", "dec21.sld", "mar2025.sld", "20250205.sld", "walls_small.sld"):
            hits, unknown = cc.find_ids_in_name(name, self.table, self.lex)
            self.assertEqual((hits, unknown), ([], []), name)

    def test_real_synthetic_tree_cli_scan(self):
        for seed in (11, 12):
            with self.subTest(seed=seed):
                MiniTableTests._synthetic_roundtrip(self, REAL_CSV, seed)


if __name__ == "__main__":
    unittest.main()
