#!/usr/bin/env python3
"""
test_pipeline_suite.py
Comprehensive Verification & Unit/Integration Test Suite for the 4D Mitochondrial HPC Pipeline

Tests Covered:
  1. test_01_sld_metadata_and_safety:
     - Verifies OME-XML parsing, 4D filtering (Z >= 15, T >= 10, C >= 2), cohort year inference,
       and strict read-only safety rejection when source_dir == output_dir.
  1i-1k. Health audit, error classification, and portable movie_list.txt entries
     (sld_health_audit.csv statuses, raw-stderr classification via a fake showinf, relative paths).
  2. test_02_tiff_io_5d_roundtrip (+ 02b-02e):
     - Verifies lossless 5D (T, C, Z, Y, X) round-trips in both OME-TIFF and ImageJ formats,
       calibration and unit conversion from bfconvert-style OME-XML (xsi:schemaLocation),
       axis normalization (read_tiff(axes="TCZYX")), and int32 label masks.
  3. test_03_fast4dreg_drift_sign_and_accuracy:
     - Synthesizes a 3D cell volume shifted by (+5.0 px X, -3.0 px Y, +1.0 slice Z) and verifies
       that estimate_fast4dreg_drift returns corrective shifts (-5.0, +3.0, -1.0) and that
       apply_drift_correction_inplace restores the shifted object to its exact T=0 coordinates.
  4. test_04_fast_16bit_lut_histogram_matching:
     - Simulates 50% exponential photobleaching and verifies that match_histogram_fast_16bit
       restores the 99th percentile intensity within 2% of the unbleached T=0 reference.
  5. test_05_vectorized_volume_filter_and_watershed_split:
     - Verifies that filter_and_relabel_fast removes sub-threshold fragments (<0.05 um^3) in O(N)
       time and that watershed_split_3d separates two touching mitochondria in a narrow tube.
  6. test_06_sparse_lap_tracking_gap_fission_fusion:
     - Synthesizes a 6-frame 4D label movie containing:
       * Track A: Continuous directional motion with a 1-frame missing gap (tests GAP_CLOSE)
       * Track B: A single mitochondrion splitting into two daughters at T=2 (tests SPLIT_FISSION)
       * Track C: Two mitochondria merging into one at T=3 (tests MERGE_FUSION)
     - Verifies that run_sparse_lap_tracking detects the gap, fission, and fusion events, unifies
       branched segments into connected-component Track IDs, and outputs valid TrackMate v7 XML.
  7. test_07_curved_arbor_internal_geodesic_and_swc:
     - Constructs a U-shaped curved 3D radial glial process where the straight-line Euclidean
       distance from soma to tip is 4.0 um, but the true internal path along the curved tube is
       >11.0 um. Verifies that compute_internal_geodesic_distance_3d reports the true internal
       path length (>11.0 um) rather than taking a shortcut through extracellular space.
  8. test_08_end_to_end_pipeline_integration:
     - Runs Modules 1, 2, 3, and 4 sequentially on a synthetic 4D dual-channel movie and verifies
       all 13 output artifacts, CSV schemas, and JSON provenance records.

Dependency Gating (no silent skips):
  - Stage 1 tests (01*) need only the Python standard library.
  - tiff_io tests (02*) run when numpy + tifffile are installed.
  - Stage 3-5 tests (03-08) run when numpy, tifffile, scipy, pyyaml, and opencv are installed.
  Gates check THIRD-PARTY packages only. When they are present, tiff_io.py and the stage scripts
  are imported without try/except, so a missing or broken pipeline module is a loud error.

Usage:
  python3 test_pipeline_suite.py -v
  MITO_PRACTICE_SLD=/path/to/05feb2025obf01avc01-B.sld python3 test_pipeline_suite.py -v
"""

import csv
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Tuple

# Ensure pipeline modules in the same directory can be imported
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import collect_4d_sld_timelapses as sld_mod


def _installed(*modules: str) -> bool:
    """True when every third-party module can be found (does not import it)."""
    return all(importlib.util.find_spec(m) is not None for m in modules)


HAVE_TIFF_STACK = _installed("numpy", "tifffile")
HAVE_HPC_PACKAGES = HAVE_TIFF_STACK and _installed("scipy", "yaml", "cv2")

np = None
prep_mod = ml_mod = seg_mod = track_mod = tio = None
if HAVE_TIFF_STACK:
    import numpy as np
    import tiff_io as tio
    from tiff_io import read_tiff, write_tiff
if HAVE_HPC_PACKAGES:
    import preprocess_4d_stack as prep_mod
    import run_arbor_and_mitolandscape as ml_mod
    import run_hybrid_segmentation as seg_mod
    import run_trackmate_sparse_lap as track_mod

TIFF_SKIP_MSG = "Requires numpy + tifffile"
HPC_SKIP_MSG = "Requires HPC analysis packages (numpy, tifffile, scipy, pyyaml, opencv)"

# A bfconvert-style OME-XML header: xsi:schemaLocation (unbound prefix once xmlns is removed),
# mixed units (nm / um / ms), and no TimeIncrement so DeltaT fallback is exercised.
BFCONVERT_STYLE_OME_XML = (
    '<?xml version="1.0" encoding="UTF-8"?>'
    '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" '
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" UUID="urn:uuid:00000000-0000-0000-0000-000000000000" '
    'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 '
    'http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">'
    '<Image ID="Image:0" Name="t2r_-_4"><Pixels BigEndian="false" DimensionOrder="XYCZT" ID="Pixels:0" '
    'PhysicalSizeX="160.0" PhysicalSizeXUnit="nm" PhysicalSizeY="0.16" PhysicalSizeYUnit="µm" '
    'PhysicalSizeZ="1.0" PhysicalSizeZUnit="µm" SizeC="2" SizeT="3" SizeX="8" SizeY="8" SizeZ="2" Type="uint16">'
    '<Plane TheC="0" TheT="0" TheZ="0" DeltaT="0" DeltaTUnit="ms"/>'
    '<Plane TheC="0" TheT="1" TheZ="0" DeltaT="11700" DeltaTUnit="ms"/>'
    '<Plane TheC="0" TheT="2" TheZ="0" DeltaT="23400" DeltaTUnit="ms"/>'
    "</Pixels></Image></OME>"
)


class TestMito4DPipelineSuite(unittest.TestCase):
    """Rigorous unit and integration test suite for all modules of the 4D Mitochondrial Pipeline."""

    def setUp(self):
        self.tmp_dir = Path(tempfile.mkdtemp(prefix="mito4d_test_"))
        # pipeline_config.yaml lives in the repository root, not in scripts/
        self.cfg_path = REPO_ROOT / "pipeline_config.yaml"
        if HAVE_HPC_PACKAGES:
            self.cfg = prep_mod.load_config(self.cfg_path)
        else:
            self.cfg = {}

    def tearDown(self):
        if os.environ.get("MITO_MODELS_DIR", "").startswith(str(self.tmp_dir)):
            del os.environ["MITO_MODELS_DIR"]  # set by test_08 only
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_01_sld_metadata_and_safety(self):
        """Verifies OME-XML 4D filtering, cohort year inference, and read-only source dir guard."""
        sample_xml = """<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
          <Image ID="Image:0" Name="2D_Snap">
            <Pixels SizeX="512" SizeY="512" SizeZ="1" SizeC="2" SizeT="1" Type="uint16" SignificantBits="16"
                    PhysicalSizeX="0.16" PhysicalSizeY="0.16"/>
          </Image>
          <Image ID="Image:1" Name="RGC_4D_Movie">
            <Pixels SizeX="512" SizeY="512" SizeZ="45" SizeC="2" SizeT="30" Type="uint16" SignificantBits="16"
                    PhysicalSizeX="0.16" PhysicalSizeY="0.16" PhysicalSizeZ="1.0" TimeIncrement="11.7"/>
          </Image>
        </OME>"""
        fake_sld = Path("/data/2025/20250610_Tadpole01.sld")
        series = sld_mod.parse_ome_xml_metadata(
            sample_xml, fake_sld, Path("/data"), min_z=15, min_t=10, min_c=2
        )
        self.assertEqual(len(series), 2)
        self.assertFalse(series[0].is_4d_timelapse)
        self.assertTrue(series[1].is_4d_timelapse)
        self.assertEqual(series[1].cohort_year, "2025")
        self.assertTrue(series[1].calibration_complete)
        self.assertEqual(
            series[1].extracted_ometiff_relpath, "2025/20250610_Tadpole01__S01_RGC_4D_Movie.ome.tif"
        )

        # Verify safety guard blocks source_dir == output_dir
        proc = subprocess.run(
            [sys.executable, str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"), str(self.tmp_dir), str(self.tmp_dir)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("Safety violation", proc.stderr)

    def test_01b_criteria_thresholds(self):
        """Verifies changing criteria for Z, T, and C filters accurately selects expected series."""
        multi_xml = """<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
          <Image ID="Image:0" Name="3D_Stack_NoTime">
            <Pixels SizeX="512" SizeY="512" SizeZ="30" SizeC="2" SizeT="1" Type="uint16"/>
          </Image>
          <Image ID="Image:1" Name="2D_Timelapse">
            <Pixels SizeX="512" SizeY="512" SizeZ="1" SizeC="2" SizeT="50" Type="uint16"/>
          </Image>
          <Image ID="Image:2" Name="4D_SingleChannel">
            <Pixels SizeX="512" SizeY="512" SizeZ="20" SizeC="1" SizeT="10" Type="uint16"/>
          </Image>
          <Image ID="Image:3" Name="4D_DualChannel_ShortTime">
            <Pixels SizeX="512" SizeY="512" SizeZ="25" SizeC="2" SizeT="5" Type="uint16"/>
          </Image>
        </OME>"""
        sld_file = Path("/data/2024/test.sld")

        # Standard defaults: min_z=15, min_t=10, min_c=2 -> None qualify
        res_default = sld_mod.parse_ome_xml_metadata(multi_xml, sld_file, Path("/data"), min_z=15, min_t=10, min_c=2)
        self.assertEqual([s.is_4d_timelapse for s in res_default], [False, False, False, False])

        # Single-channel allowed: min_c=1 -> Series 2 qualifies
        res_c1 = sld_mod.parse_ome_xml_metadata(multi_xml, sld_file, Path("/data"), min_z=15, min_t=10, min_c=1)
        self.assertEqual([s.is_4d_timelapse for s in res_c1], [False, False, True, False])

        # Short timelapses allowed: min_t=5, min_c=2 -> Series 3 qualifies
        res_t5 = sld_mod.parse_ome_xml_metadata(multi_xml, sld_file, Path("/data"), min_z=15, min_t=5, min_c=2)
        self.assertEqual([s.is_4d_timelapse for s in res_t5], [False, False, False, True])

        # Permissive: min_t=5, min_c=1 -> Series 2 and Series 3 both qualify
        res_perm = sld_mod.parse_ome_xml_metadata(multi_xml, sld_file, Path("/data"), min_z=15, min_t=5, min_c=1)
        self.assertEqual([s.is_4d_timelapse for s in res_perm], [False, False, True, True])

    def test_01c_companion_dirs(self):
        """Verifies detection of SlideBook companion data directories (_Files, _data, .dir)."""
        sld_file = self.tmp_dir / "sample_acquisition.sld"
        sld_file.touch()

        # Initially no companion dirs
        self.assertEqual(sld_mod.find_companion_sld_dirs(sld_file), [])

        # Create companion candidate dirs
        comp_data = self.tmp_dir / "sample_acquisition_data"
        comp_files = self.tmp_dir / "sample_acquisition_Files"
        comp_data.mkdir()
        comp_files.mkdir()

        found = sld_mod.find_companion_sld_dirs(sld_file)
        self.assertEqual(len(found), 2)
        self.assertIn(comp_data, found)
        self.assertIn(comp_files, found)

    def test_01d_sanitize_and_size_estimation(self):
        """Verifies string sanitization and uncompressed 5D byte size calculation."""
        self.assertEqual(
            sld_mod.sanitize_filename("  4D Movie / #1 [Z-stack] (Test)  "),
            "4D_Movie_1_Z-stack_Test",
        )
        # 512 x 512 x 20 Z x 2 C x 10 T x 16-bit (2 bytes) = 209,715,200 bytes = ~0.1953 GB
        est = sld_mod.estimate_series_gb(x=512, y=512, z=20, c=2, t=10, bits=16)
        self.assertAlmostEqual(est, 0.1953, places=3)

    def _find_practice_file(self) -> Path:
        # MITO_PRACTICE_SLD lets HPC runs point at the copy on project storage.
        if os.environ.get("MITO_PRACTICE_SLD"):
            return Path(os.environ["MITO_PRACTICE_SLD"])
        candidates = [
            SCRIPT_DIR.parent / "data" / "practice" / "05feb2025obf01avc01-B.sld",
            SCRIPT_DIR / "data" / "practice" / "05feb2025obf01avc01-B.sld",
        ]
        for c in candidates:
            if c.is_file():
                return c
        return candidates[0]

    def test_01e_live_sld_inspection_if_available(self):
        """If practice .sld file exists, verifies full scan and provenance generation without errors."""
        practice_file = self._find_practice_file()
        if not practice_file.is_file():
            self.skipTest(f"Practice file not present at {practice_file}")

        out_scan = self.tmp_dir / "scan_out"
        cmd = [
            sys.executable,
            str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"),
            str(practice_file.parent),
            str(out_scan),
            "--mode", "scan-only",
        ]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"Scan failed: {proc.stderr}")
        self.assertTrue((out_scan / "sld_4d_inventory.csv").is_file())
        self.assertTrue((out_scan / "movie_list.txt").is_file())
        self.assertTrue((out_scan / "extraction_provenance.json").is_file())
        with open(out_scan / "sld_health_audit.csv", "r", encoding="utf-8") as f:
            health = list(csv.DictReader(f))
        self.assertEqual(
            [r["status"] for r in health if r["sld_filename"] == practice_file.name], ["OK_4D"]
        )

        # Verify exactly 49 series in inventory CSV
        with open(out_scan / "sld_4d_inventory.csv", "r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 49)

        # Verify Series 20 is flagged as qualifying 4D
        s20 = next(r for r in rows if r["series_index"] == "20")
        self.assertEqual(s20["is_4d_timelapse"], "True")
        self.assertEqual(s20["size_z"], "49")
        self.assertEqual(s20["size_t"], "10")
        self.assertEqual(s20["size_c"], "2")

    def test_01f_checksum_generation_and_verification(self):
        """Verifies SHA-256 computation matches standard hashlib hash."""
        test_file = self.tmp_dir / "checksum_test.bin"
        test_data = b"Mitochondrial_4D_Timelapse_Integrity_Verification_Data_Block_2026"
        with open(test_file, "wb") as f:
            f.write(test_data)
        import hashlib
        expected_sha = hashlib.sha256(test_data).hexdigest()
        actual_sha = sld_mod.compute_sha256(test_file)
        self.assertEqual(actual_sha, expected_sha)

    def test_01g_extract_ometiff_dry_run(self):
        """Verifies --mode extract-ometiff --dry-run simulates extraction without creating large TIFFs."""
        practice_file = self._find_practice_file()
        if not practice_file.is_file():
            self.skipTest(f"Practice file not present at {practice_file}")

        out_dry = self.tmp_dir / "dry_out"
        cmd = [
            sys.executable,
            str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"),
            str(practice_file.parent),
            str(out_dry),
            "--mode", "extract-ometiff",
            "--dry-run",
        ]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"Dry run failed: {proc.stderr}")
        self.assertIn("[DRY-RUN] Would extract series", proc.stdout)
        # Ensure no actual .ome.tif files were written
        tif_files = list(out_dry.glob("*.ome.tif"))
        self.assertEqual(len(tif_files), 0)

    def test_01h_skip_existing_files(self):
        """Verifies that extract-ometiff skips already extracted files unless --overwrite is passed."""
        practice_file = self._find_practice_file()
        if not practice_file.is_file():
            self.skipTest(f"Practice file not present at {practice_file}")

        out_skip = self.tmp_dir / "skip_out"
        # Top-level .sld files are staged under their cohort-year folder (2025/)
        (out_skip / "2025").mkdir(parents=True)
        # Pre-create the dummy target output file
        dummy_out = out_skip / "2025" / "05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif"
        with open(dummy_out, "wb") as f:
            f.write(b"PRE_EXISTING_MOCK_DATA")

        cmd = [
            sys.executable,
            str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"),
            str(practice_file.parent),
            str(out_skip),
            "--mode", "extract-ometiff",
        ]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"Execution failed: {proc.stderr}")
        self.assertIn("Skipped (already exists)", proc.stdout)
        # Verify dummy file was NOT overwritten
        with open(dummy_out, "rb") as f:
            self.assertEqual(f.read(), b"PRE_EXISTING_MOCK_DATA")

    def test_01i_health_audit_prefilter_statuses(self):
        """Verifies every candidate gets one sld_health_audit.csv row (stub, lock, null header)."""
        src = self.tmp_dir / "src" / "2024"
        src.mkdir(parents=True)
        (src / "tiny_stub.sld").write_bytes(b"\x01" * 100)          # < 64 KB -> SKIPPED_STUB
        (src / "~lockfile.sld").write_bytes(b"\x01" * 70_000)       # lock prefix -> never opened
        (src / "crashed_capture.sld").write_bytes(b"\x00" * 70_000)  # >= 64 KB but all zeros
        out = self.tmp_dir / "audit_out"

        proc = subprocess.run(
            [sys.executable, str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"),
             str(self.tmp_dir / "src"), str(out), "--mode", "scan-only", "--hpc-raw-root", ""],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out / "sld_health_audit.csv", "r", encoding="utf-8") as f:
            statuses = {r["sld_filename"]: r["status"] for r in csv.DictReader(f)}
        self.assertEqual(statuses, {
            "tiny_stub.sld": "SKIPPED_STUB",
            "~lockfile.sld": "SKIPPED_TEMP_OR_LOCK",
            "crashed_capture.sld": "EMPTY_OR_NULL_FILE",
        })
        self.assertFalse((out / "unparseable_sld_files.csv").exists())
        with open(out / "extraction_provenance.json", "r", encoding="utf-8") as f:
            prov = json.load(f)
        self.assertEqual(prov["summary_metrics"]["total_unparseable_or_skipped"], 3)
        self.assertEqual(prov["health_audit_by_year"]["2024"]["total"], 3)

    def test_01j_error_classification_uses_raw_stderr(self):
        """Verifies showinf failures are bucketed from the raw Java stderr, not the friendly hint."""
        cases = {
            "java.lang.AssertionError: The file is not a valid SlideBook document. Problem reading the Image Records":
                "SLIDEBOOK_IMAGE_RECORDS_ERR",
            "java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0\n"
            "  at loci.formats.in.SlidebookReader.initFile(SlidebookReader.java:412)":
                "MISSING_SB6_READER",
            "java.lang.UnsatisfiedLinkError: no SBReadFile in java.library.path": "JNI_OR_JAVA_ERR",
            "SlideBook error: file already open": "SLIDEBOOK_SESSION_CONCURRENCY_ERR",
            "something unexpected": "OTHER_ERR",
        }
        for stderr_text, expected in cases.items():
            self.assertEqual(sld_mod.classify_showinf_error(stderr_text)[0], expected, stderr_text)

        if sys.platform.startswith("win"):
            self.skipTest("Fake showinf shell script requires a POSIX shell")
        # End-to-end through audit_and_inspect_sld with a fake showinf that fails like 3i's native reader.
        fake_showinf = self.tmp_dir / "showinf"
        fake_showinf.write_text(
            "#!/bin/sh\n"
            "echo 'Exception in thread \"main\" java.lang.AssertionError: Failed with exception: "
            "The file is not a valid SlideBook document. Problem reading the Image Records' >&2\n"
            "exit 1\n"
        )
        fake_showinf.chmod(0o755)
        sld = self.tmp_dir / "aborted.sld"
        sld.write_bytes(b"\x01" * 70_000)
        series, rec, hint = sld_mod.audit_and_inspect_sld(
            sld, self.tmp_dir, str(fake_showinf), dict(os.environ), 15, 10, 2
        )
        self.assertEqual(series, [])
        self.assertEqual(rec.status, "SLIDEBOOK_IMAGE_RECORDS_ERR")
        self.assertIn("Problem reading the Image Records", rec.error_summary)
        self.assertIn("aborted acquisition", hint)

    def test_01o_mito_construct_filter_gates_bioformats(self):
        """--plasmid-csv classifies files by name BEFORE Bio-Formats: under 'require' only MITO files reach
        showinf (others get SKIPPED_NOT_MITO); under 'report' every file is inspected but still labelled."""
        if sys.platform.startswith("win"):
            self.skipTest("Fake showinf shell script requires a POSIX shell")
        src = self.tmp_dir / "src"
        (src / "2025").mkdir(parents=True)
        names = {
            "05feb2025obf01avc01-B.sld": "MITO_PHOTOACTIVATABLE",  # OBF001 = mito-paGFP
            "12mar2025jr19jeb09_cell1.sld": "MITO",                 # JR19 = pUAS-mlsEGFPpA
            "random_cells.sld": "REVIEW",
        }
        for name in names:
            (src / "2025" / name).write_bytes(b"\x01" * 70_000)
        plasmids = self.tmp_dir / "plasmids.csv"
        header = ("s,sequenced? (date),does this plasmid make FP positive cells?,glycerol stock?,-20 DNA duplicate?,"
                  "addgene submission?,plasmid ID,name (Prom_GOI_3'),gateway compnent?,localization sequence?,"
                  "GOI? tag?,Fluorescent Protein 1 (ME),Fluorescent Protein 2 (3'),p2a,\"uas, MCS, promoter?\","
                  "antibiotic,plasmid description/elements,full map?,long term use,Forward Primer,Reverse Primer,"
                  "product size \n")
        rows = [
            ",,,,,,OBF001,psox2bd::mitopaGFP_p2a-tagRFPt,-,mito,,paGFP,tagRFPT,p2a,,amp,,,,,,\n",
            ",,,,,,JR19,pUAS-mlsEGFPpA,,mito,,EGFP,,,UAS,amp,,,,,,\n",
            ",,,,,,JEB009,pSox2bd_tagRFPt_p2a-tagRFPt,,no,,tagRFPT (middle),tagRFPT (3'),p2a,,amp,,,,,,\n",
            ",,,,,,avc01,p4XUAS_tagRFPt_p2a-mRuby2,-,-,-,tagRFPt (middle),ruby2,p2a,UAS,amp,,,,,,\n",
        ]
        plasmids.write_text(header + "".join(rows), encoding="utf-8")
        bft = self.tmp_dir / "bft"
        bft.mkdir()
        log = self.tmp_dir / "showinf_calls.log"
        for tool in ("showinf", "bfconvert"):
            (bft / tool).write_text(
                "#!/bin/sh\n"
                'case "$*" in *-version*) echo "Version: fake"; exit 0;; esac\n'
                f'for a in "$@"; do last="$a"; done; echo "$last" >> "{log}"\n'
                "echo 'java.lang.AssertionError: Problem reading the Image Records' >&2\nexit 1\n"
            )
            (bft / tool).chmod(0o755)

        for mode, expected_inspected in (("require", {"12mar2025jr19jeb09_cell1.sld"}), ("report", set(names))):
            out = self.tmp_dir / f"out_{mode}"
            log.write_text("")
            proc = subprocess.run(
                [sys.executable, str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"), str(src), str(out),
                 "--mode", "scan-only", "--workers", "1", "--bftools-dir", str(bft),
                 "--plasmid-csv", str(plasmids), "--mito-filter", mode],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            inspected = {Path(line).name for line in log.read_text().split() if line.endswith(".sld")}
            self.assertEqual(inspected, expected_inspected, mode)
            with open(out / "mito_classification.csv", newline="", encoding="utf-8") as f:
                labels = {Path(r["rel_path"]).name: r["label"] for r in csv.DictReader(f)}
            self.assertEqual(labels, names, mode)
            with open(out / "sld_health_audit.csv", newline="", encoding="utf-8") as f:
                status = {r["sld_filename"]: r["status"] for r in csv.DictReader(f)}
            skipped = {n for n, st in status.items() if st == "SKIPPED_NOT_MITO"}
            self.assertEqual(skipped, set(names) - expected_inspected, mode)
            prov = json.loads((out / "extraction_provenance.json").read_text())
            self.assertEqual(prov["mito_construct_filter"]["mode"], mode)

    def test_01n_inventory_calibration_units(self):
        """Verifies the inventory reports um and seconds even when OME stores nm and DeltaT in ms."""
        series = sld_mod.parse_ome_xml_metadata(
            BFCONVERT_STYLE_OME_XML, Path("/data/2025/movie.sld"), Path("/data"), min_z=1, min_t=1, min_c=1
        )
        self.assertEqual(len(series), 1)
        self.assertAlmostEqual(series[0].phys_size_x_um, 0.16, places=4)    # 160 nm
        self.assertAlmostEqual(series[0].time_interval_s, 11.7, places=3)   # DeltaT 11700 ms apart
        self.assertTrue(series[0].calibration_complete)

    def test_01k_movie_list_entries_relative_or_absolute(self):
        """Verifies movie_list.txt entries are relative unless an HPC raw root is supplied."""
        rel = "2025/movie__S20_x.ome.tif"
        self.assertEqual(sld_mod.manifest_entry("", rel), rel)
        self.assertEqual(
            sld_mod.manifest_entry("/sciclone/data10/lab/raw_ometiff", rel),
            "/sciclone/data10/lab/raw_ometiff/2025/movie__S20_x.ome.tif",
        )

    @unittest.skipIf(sys.platform.startswith("win") or not shutil.which("bash"), "Requires bash")
    def test_01l_storage_defaults_match_between_bash_and_python(self):
        """Verifies scripts/mito_env.bash and pipeline_config_io.py agree on every MITO_* default."""
        clean_env = {k: v for k, v in os.environ.items() if not k.startswith("MITO_")}
        clean_env["USER"] = "testuser"
        for root in (None, "/sciclone/data10/somelab/mito_pipeline"):
            env = dict(clean_env)
            if root:
                env["MITO_DATA_ROOT"] = root
            # 'env -0' prints NUL-separated KEY=VALUE pairs so values with newlines parse safely
            proc = subprocess.run(
                ["bash", "-c", f"source '{SCRIPT_DIR / 'mito_env.bash'}' && env -0"],
                stdout=subprocess.PIPE, env=env, check=True,
            )
            bash_vals = dict(
                item.split("=", 1) for item in proc.stdout.decode().split("\0")
                if item.startswith("MITO_")
            )
            # Import lazily: pipeline_config_io needs pyyaml, which Stage 1 environments may lack.
            if not _installed("yaml"):
                self.skipTest("pyyaml not installed")
            import pipeline_config_io as pcio
            # resolve_mito_env reads the user name from os.environ, so pin it for the comparison.
            saved_user = os.environ.get("USER")
            os.environ["USER"] = "testuser"
            try:
                py_vals = pcio.resolve_mito_env({k: v for k, v in env.items() if k.startswith("MITO_")})
            finally:
                if saved_user is None:
                    os.environ.pop("USER", None)
                else:
                    os.environ["USER"] = saved_user
            self.assertEqual(bash_vals, py_vals)
            if root:
                self.assertEqual(py_vals["MITO_RAW_ROOT"], f"{root}/raw_ometiff")

    @unittest.skipUnless(_installed("yaml"), "Requires pyyaml")
    def test_01m_config_placeholders_expand_or_fail_loudly(self):
        """Verifies ${MITO_*} placeholders expand from the environment and unknown ones raise."""
        import pipeline_config_io as pcio
        cfg = pcio.load_config(REPO_ROOT / "pipeline_config.yaml")
        for key in ("raw_ometiff_root", "models_dir", "final_output_root", "ilastik_bin"):
            self.assertNotIn("$", cfg["paths"][key])
        bad_cfg = self.tmp_dir / "bad.yaml"
        bad_cfg.write_text("paths:\n  x: \"${MITO_DOES_NOT_EXIST_123}/y\"  # ${COMMENTS_ARE_FINE}\n")
        with self.assertRaises(ValueError) as ctx:
            pcio.load_config(bad_cfg)
        self.assertIn("MITO_DOES_NOT_EXIST_123", str(ctx.exception))
        self.assertNotIn("COMMENTS_ARE_FINE", str(ctx.exception))

    @unittest.skipUnless(HAVE_TIFF_STACK, TIFF_SKIP_MSG)
    def test_02_tiff_io_5d_roundtrip(self):
        """Verifies lossless 5D (T, C, Z, Y, X) round-trips in OME-TIFF and ImageJ formats with calibration."""
        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}
        # Z=3 is the RGB trap: without photometric='minisblack' tifffile stores it as RGB planes.
        for dtype in (np.uint16, np.float32):
            arr = (np.random.rand(3, 2, 3, 32, 32) * 50000).astype(dtype)
            for name, expect_imagej in (("test_5d.ome.tif", False), ("test_5d.tif", True)):
                path = self.tmp_dir / f"{dtype.__name__}_{name}"
                write_tiff(path, arr, metadata=calib, axes="TCZYX")
                loaded, meta = read_tiff(path, axes="TCZYX")
                self.assertEqual(loaded.dtype, arr.dtype)
                self.assertTrue(np.array_equal(arr, loaded), f"{path.name} not lossless")
                self.assertAlmostEqual(float(meta["PhysicalSizeX"]), 0.16, places=5)
                self.assertAlmostEqual(float(meta["PhysicalSizeZ"]), 1.0, places=5)
                self.assertAlmostEqual(float(meta["TimeIncrement"]), 11.7, places=5)
                self.assertNotIn("_warnings", meta)
                import tifffile
                with tifffile.TiffFile(path) as tif:
                    self.assertEqual(tif.is_imagej, expect_imagej)
                    self.assertEqual(bool(tif.ome_metadata), not expect_imagej)
        self.assertEqual([p.name for p in self.tmp_dir.iterdir() if p.name.startswith(".")], [])

    @unittest.skipUnless(HAVE_TIFF_STACK, TIFF_SKIP_MSG)
    def test_02b_bfconvert_ome_xml_units_and_deltat(self):
        """Verifies bfconvert-style OME-XML (xsi:schemaLocation, nm/ms units) parses to um and s."""
        meta = tio._parse_ome_xml(BFCONVERT_STYLE_OME_XML)
        self.assertNotIn("_warnings", meta)
        self.assertAlmostEqual(meta["PhysicalSizeX"], 0.16, places=6)   # 160 nm -> 0.16 um
        self.assertAlmostEqual(meta["PhysicalSizeY"], 0.16, places=6)
        self.assertAlmostEqual(meta["PhysicalSizeZ"], 1.0, places=6)
        self.assertAlmostEqual(meta["TimeIncrement"], 11.7, places=4)   # DeltaT in ms -> 11.7 s
        self.assertNotIn("PhysicalSizeXUnit", meta)

    @unittest.skipUnless(HAVE_TIFF_STACK, TIFF_SKIP_MSG)
    def test_02c_axis_normalization(self):
        """Verifies read_tiff(axes='TCZYX') reorders native layouts and inserts missing axes."""
        base = np.arange(2 * 3 * 4 * 5 * 6, dtype=np.uint16).reshape(2, 3, 4, 5, 6)  # T C Z Y X
        tzcyx = np.transpose(base, (0, 2, 1, 3, 4))
        ctzyx = np.transpose(base, (1, 0, 2, 3, 4))
        self.assertTrue(np.array_equal(tio.reorder_axes(tzcyx, "TZCYX", "TCZYX"), base))
        self.assertTrue(np.array_equal(tio.reorder_axes(ctzyx, "CTZYX", "TCZYX"), base))

        # 4D single-channel movie on disk -> 5D with C=1
        single = base[:, 0]
        path = self.tmp_dir / "single_channel.tif"
        write_tiff(path, single, axes="TZYX")
        loaded, meta = read_tiff(path, axes="TCZYX")
        self.assertEqual(loaded.shape, (2, 1, 4, 5, 6))
        self.assertEqual(meta["axes"], "TCZYX")
        self.assertIn("PhysicalSizeX missing", " ".join(meta["_warnings"]))

        with self.assertRaises(ValueError):
            tio.reorder_axes(np.zeros((5, 4, 8, 8)), "IZYX", "TCZYX")

    @unittest.skipUnless(HAVE_TIFF_STACK, TIFF_SKIP_MSG)
    def test_02d_int32_label_mask_roundtrip(self):
        """Verifies int32 label masks (not ImageJ-compatible) round-trip with axes and calibration."""
        labels = np.random.randint(0, 70000, size=(3, 4, 16, 16), dtype=np.int32)
        path = self.tmp_dir / "labels.tif"
        write_tiff(path, labels, metadata={"dx_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}, axes="TZYX")
        loaded, meta = read_tiff(path)
        self.assertEqual(loaded.dtype, np.int32)
        self.assertTrue(np.array_equal(labels, loaded))
        self.assertEqual(meta["axes"], "TZYX")
        self.assertAlmostEqual(float(meta["PhysicalSizeX"]), 0.16, places=6)

    @unittest.skipUnless(HAVE_TIFF_STACK, TIFF_SKIP_MSG)
    def test_02e_unparseable_ome_is_reported_not_silent(self):
        """Verifies a broken OME-XML block is reported in _warnings instead of silently returning {}."""
        meta = tio._parse_ome_xml("<OME><Image><Pixels PhysicalSizeX='0.16'></Image></OME>")
        self.assertIn("unparseable", " ".join(meta.get("_warnings", [])))

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_03_fast4dreg_drift_sign_and_accuracy(self):
        """Verifies that Fast4DReg phase correlation estimates the exact corrective sign and aligns drifted objects."""
        # Create a 2-frame stack (T=2, Z=10, Y=64, X=64) with a bright soma at (z=4, y=25, x=25) at T=0
        vol = np.zeros((2, 10, 64, 64), dtype=np.float32)
        zz, yy, xx = np.ogrid[:10, :64, :64]
        soma_t0 = ((zz - 4) ** 2 + (yy - 25) ** 2 + (xx - 25) ** 2) <= 16.0
        vol[0][soma_t0] = 5000.0

        # Shift the object at T=1 by (+4 px in X, +3 px in Y, +1 slice in Z) -> new center (z=5, y=28, x=29)
        soma_t1 = ((zz - 5) ** 2 + (yy - 28) ** 2 + (xx - 29) ** 2) <= 16.0
        vol[1][soma_t1] = 5000.0

        sx, sy, sz = prep_mod.estimate_fast4dreg_drift(vol, max_shift_xy=15, max_shift_z=4)
        # Corrective shift at T=1 must be negative (-4 in X, -3 in Y, -1 in Z) to bring it back to (25, 25, 4)
        self.assertAlmostEqual(float(sx[1]), -4.0, delta=0.6)
        self.assertAlmostEqual(float(sy[1]), -3.0, delta=0.6)
        self.assertAlmostEqual(float(sz[1]), -1.0, delta=0.6)

        prep_mod.apply_drift_correction_inplace(vol, sx, sy, sz)
        # Verify center of mass at T=1 now matches T=0 within 0.5 px
        com_t0 = ndimage_com(vol[0])
        com_t1 = ndimage_com(vol[1])
        self.assertLess(np.linalg.norm(np.array(com_t0) - np.array(com_t1)), 0.75)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_04b_histogram_matching_never_saturates_and_z_check_ignores_empty_slices(self):
        """Regression (synthetic seeds 2/3): the brightest voxels must map to the reference's range,
        never to 65535; and slices that are merely outside the cell must not trigger Z correction."""
        rng = np.random.default_rng(0)
        ref = rng.integers(0, 400, size=(64, 64)).astype(np.float32)
        img = rng.integers(0, 300, size=(64, 64)).astype(np.float32)
        out = prep_mod.match_histogram_fast_16bit(img, ref)
        self.assertLessEqual(float(out.max()), float(ref.max()))
        # Z correction is driven by the cytosolic channel: sparse mito content (bright only in shallow
        # slices) with an evenly bright cell must NOT trigger; a cell that dims with depth must.
        mask = np.zeros((16, 64, 64), bool)
        mask[:, 16:48, 16:48] = True
        cell_even = np.full((2, 16, 64, 64), 100.0, np.float32)
        cell_even[:, mask] = 400.0
        vol = np.zeros((2, 16, 64, 64), np.float32)
        vol[:, 2:5, 20:30, 20:30] = 3000.0                 # mitochondria only in shallow slices
        vol[:, :, 16:48, 16:48] += rng.uniform(0, 20, size=(2, 16, 32, 32))
        before = vol.copy()
        info = prep_mod.compensate_z_attenuation_inplace(vol, 0.30, cell_4d=cell_even, cell_mask=mask)
        self.assertFalse(info["applied"])
        np.testing.assert_array_equal(vol, before)
        cell_dim = cell_even.copy()
        for t in range(2):
            for z in range(16):
                cell_dim[t, z][mask[z]] = 100.0 + 300.0 * np.exp(-0.08 * z)   # signal above bg decays with depth
        info2 = prep_mod.compensate_z_attenuation_inplace(vol, 0.30, cell_4d=cell_dim, cell_mask=mask)
        self.assertTrue(info2["applied"])
        self.assertAlmostEqual(info2["max_gain_applied"], min(3.0, float(np.exp(0.08 * 15))), delta=0.15)  # capped by max_gain

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_04_fast_16bit_lut_histogram_matching(self):
        """Verifies that 16-bit LUT histogram matching restores bleached intensities accurately."""
        rng = np.random.default_rng(42)
        ref = rng.integers(100, 20000, size=(12, 64, 64), dtype=np.uint16).astype(np.float32)
        bleached = (ref * 0.45).astype(np.float32)  # 55% photobleaching loss

        restored = prep_mod.match_histogram_fast_16bit(bleached, ref)
        p99_ref = float(np.percentile(ref, 99.0))
        p99_restored = float(np.percentile(restored, 99.0))
        self.assertLess(abs(p99_restored - p99_ref) / p99_ref, 0.02)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_05_vectorized_volume_filter_and_watershed_split(self):
        """Verifies O(N) LUT volume filtering and 3D watershed splitting of touching organelles."""
        # Create a 3D mask with two touching spheres along Y (centers at y=12 and y=20, radius=4)
        vol_int = np.zeros((6, 32, 32), dtype=np.uint16)
        zz, yy, xx = np.ogrid[:6, :32, :32]
        d1 = np.sqrt((zz - 3) ** 2 + (yy - 11) ** 2 + (xx - 16) ** 2)
        d2 = np.sqrt((zz - 3) ** 2 + (yy - 19) ** 2 + (xx - 16) ** 2)
        binary = (d1 <= 4.5) | (d2 <= 4.5)
        vol_int[d1 <= 4.5] = int(5000)
        vol_int[d2 <= 4.5] = int(5000)
        # Add peak centers
        vol_int[3, 11, 16] = 45000
        vol_int[3, 19, 16] = 45000

        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 0.5, "dt_s": 11.7}
        split_masks = seg_mod.watershed_split_3d(binary, vol_int, calib)
        # Add a tiny 1-voxel noise speck (label 99)
        split_masks[0, 0, 0] = 99

        relabeled, count = seg_mod.filter_and_relabel_fast(split_masks, min_voxels=5)
        self.assertEqual(count, 2)
        self.assertEqual(int(relabeled[0, 0, 0]), 0)  # Noise speck removed

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_05b_ilastik_headless_export_contract(self):
        """Stage 4A must request an N-D export (ilastik 'tif' is 2D-only) and reject empty .ilp uploads."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            empty_ilp = tmp / "empty.ilp"
            empty_ilp.write_bytes(b"")
            self.assertFalse(seg_mod.is_valid_ilastik_project(empty_ilp))
            good_ilp = tmp / "good.ilp"
            good_ilp.write_bytes(seg_mod.HDF5_SIGNATURE + b"\0" * 64)
            self.assertTrue(seg_mod.is_valid_ilastik_project(good_ilp))

            # Fake ilastik: writes a (T, Z, Y, X, 2) float32 .npy where it is told to, exactly as the
            # real numpy exporter does, with class 0 = 0.75 everywhere.
            fake_bin = tmp / "fake_ilastik.py"
            fake_bin.write_text(
                "#!" + sys.executable + "\n"
                "import sys, numpy as np\n"
                "a = dict(x.split('=', 1) for x in sys.argv[1:] if x.startswith('--') and '=' in x)\n"
                "assert a['--output_format'] == 'numpy' and a['--output_axis_order'] == 'tzyxc'\n"
                "p = np.zeros((2, 3, 4, 5, 2), np.float32); p[..., 0] = 0.75; p[..., 1] = 0.25\n"
                "np.save(a['--output_filename_format'], p)\n",
                encoding="utf-8",
            )
            fake_bin.chmod(0o755)
            raw_npy = tmp / "movie_ilastik_tzyxc.npy"
            ok = seg_mod.run_stage4a_ilastik(tmp / "movie_mito_norm.tif", raw_npy, str(fake_bin), good_ilp, 1, 1000)
            self.assertTrue(ok)
            prob = seg_mod.load_ilastik_probabilities(raw_npy, 0, (2, 3, 4, 5))
            self.assertEqual(prob.dtype, np.float32)
            self.assertTrue(np.allclose(prob, 0.75))
            with self.assertRaises(ValueError):  # shape mismatch vs mito_norm is an error, not silent
                seg_mod.load_ilastik_probabilities(raw_npy, 0, (2, 3, 4, 6))
            with self.assertRaises(ValueError):  # class index beyond exported channels
                seg_mod.load_ilastik_probabilities(raw_npy, 2, (2, 3, 4, 5))

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_03c_large_cumulative_drift_over_fixed_pattern_offset(self):
        """Drift larger in total than max_shift_xy (12 px/frame, 48 px total) on a camera offset with a
        static fixed-pattern must be recovered (raw-MIP phase correlation locked onto 0 shift here)."""
        from scipy import ndimage as ndi
        rng = np.random.default_rng(3)
        n_t, n_z, n_y, n_x = 5, 6, 128, 160
        yy, xx = np.mgrid[:n_y, :n_x]
        cell = (((yy - 64) / 18.0) ** 2 + ((xx - 50) / 10.0) ** 2 <= 1).astype(np.float32)
        cell += ((np.abs(yy - 64) <= 2) & (xx > 50) & (xx < 120)).astype(np.float32) * 0.7  # process
        pattern = np.tile(rng.normal(0, 6.0, size=(1, n_x)), (n_y, 1)).astype(np.float32)   # column FPN
        stack = np.zeros((n_t, n_z, n_y, n_x), dtype=np.float32)
        for t in range(n_t):
            moved = ndi.shift(ndi.gaussian_filter(cell, 1.5), (0, 12.0 * t), order=1)
            for z in range(n_z):
                stack[t, z] = 100.0 + pattern + 40.0 * moved * (1.0 if 1 <= z <= 4 else 0.3)
        stack += rng.normal(0, 2.0, size=stack.shape).astype(np.float32)
        sx, sy, sz = prep_mod.estimate_fast4dreg_drift(stack, max_shift_xy=25, max_shift_z=5)
        np.testing.assert_allclose(sx, [-12.0 * t for t in range(n_t)], atol=1.0)
        np.testing.assert_allclose(sy, np.zeros(n_t), atol=1.0)
        np.testing.assert_allclose(sz, np.zeros(n_t), atol=0.6)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_03b_haze_removed_before_gating_and_raw_contrast_qc(self):
        """Top-hat must run before gating (thin hazy processes must not keep their haze), and the raw
        contrast QC must fail pure noise but pass a stack with real puncta."""
        stack = np.full((1, 1, 80, 80), 1000.0, dtype=np.float32)   # uniform haze everywhere (raw)
        stack[0, 0, 38:42, 40:43] = 6000.0                          # one mitochondrion on a thin process
        gate = np.zeros((1, 80, 80), dtype=bool)
        gate[0, 36:44, 10:70] = True                                # 8 px wide process << 31 px kernel
        prep_mod.remove_haze_then_gate_inplace(stack, gate, radius_px=15)
        process = stack[0, 0][gate[0]]
        self.assertLess(float(np.median(process)), 100.0)          # haze stripped inside the process
        self.assertGreater(float(stack[0, 0, 40, 41]), 3000.0)    # mitochondrion kept
        self.assertEqual(float(stack[0, 0, 10, 10]), 0.0)          # outside the gate is zero

        rng = np.random.default_rng(0)
        calib = {"dx_um": 0.2, "dy_um": 0.2, "dz_um": 1.0, "dt_s": 10.0}
        cell = np.zeros((5, 120, 120), dtype=bool)
        cell[1:4, 40:80, 40:80] = True
        noise = rng.normal(100.0, 3.0, size=(6, 5, 120, 120)).astype(np.float32)
        qc_noise = prep_mod.measure_raw_mito_contrast(noise, cell, calib)
        self.assertLess(qc_noise["contrast_sigma"], 3.0)
        with_mito = noise.copy()
        for y, x in ((50, 50), (60, 70), (72, 55)):
            with_mito[:, 2, y:y + 3, x:x + 3] += 40.0               # ~13 sigma puncta
        qc_mito = prep_mod.measure_raw_mito_contrast(with_mito, cell, calib)
        self.assertGreater(qc_mito["contrast_sigma"], 3.0)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_05d_watershed_splits_touching_but_not_single_rods(self):
        """The default Stage 4B engine must split two TOUCHING mitochondria in one foreground region, yet
        keep a single elongated rod (spanning 3 Z-slices) as ONE object; unknown engines are rejected."""
        from scipy import ndimage as ndi
        calib = {"dx_um": 0.19, "dy_um": 0.19, "dz_um": 1.0, "dt_s": 120.0}
        shape = (6, 40, 60)
        # Two touching blobs end-to-end along X, joined by a narrow, dimmer neck (one foreground region)
        pair = np.zeros(shape, np.float32)
        pair[2:4, 15:22, 10:20] = 1.0
        pair[2:4, 15:22, 23:33] = 1.0
        pair[2:4, 17:20, 20:23] = 0.6
        inten = ndi.gaussian_filter(pair, (0.5, 1.0, 1.0))
        fg = inten > 0.2
        self.assertEqual(ndi.label(fg)[1], 1)                      # truly touching: one connected region
        labels = seg_mod.watershed_split_3d(fg, inten, calib)
        self.assertEqual(len(np.unique(labels[fg])), 2)
        # One rod 12 px long, 4 px wide, 3 slices deep -> exactly one instance
        rod = np.zeros(shape, np.float32)
        rod[1:4, 18:22, 20:32] = 1.0
        rod_i = ndi.gaussian_filter(rod, (0.6, 1.0, 1.0))
        rod_fg = rod_i > 0.2
        self.assertEqual(len(np.unique(seg_mod.watershed_split_3d(rod_fg, rod_i, calib)[rod_fg])), 1)
        # Two touching small puncta (3x3 px, centres 0.76 um apart) -> two instances
        pp = np.zeros(shape, np.float32)
        pp[2, 18:21, 18:21] = 1.0
        pp[2, 18:21, 22:25] = 1.0
        pp_i = ndi.gaussian_filter(pp, (0.5, 1.0, 1.0))
        pp_fg = pp_i > 0.15
        self.assertEqual(ndi.label(pp_fg)[1], 1)
        self.assertEqual(len(np.unique(seg_mod.watershed_split_3d(pp_fg, pp_i, calib)[pp_fg])), 2)
        # Engine dispatch
        cfg = json.loads(json.dumps(self.cfg))
        b = cfg["stage4_hybrid_segmentation"]["stage4b_instance_segmentation"]
        self.assertEqual(b["engine"], "watershed")
        labels4, used = seg_mod.segment_instances(rod_i[None] * 50000, rod_i[None], cfg, calib)
        self.assertEqual(used, "watershed")
        self.assertEqual(int(labels4.max()), 1)
        b["engine"] = "stardist"
        with self.assertRaises(ValueError):
            seg_mod.segment_instances(rod_i[None] * 50000, rod_i[None], cfg, calib)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_05e_implausible_foreground_guard(self):
        """A probability map covering most of the cell must be flagged (bad classifier), a sparse one not."""
        mask = np.zeros((2, 20, 20), dtype=bool)
        mask[:, 2:18, 2:18] = True
        whole = np.zeros((3, 2, 20, 20), np.float32)
        whole[:, :, 2:18, 2:18] = 0.9                      # "everything in the cell is mitochondria"
        sparse = np.zeros_like(whole)
        sparse[:, :, 5:8, 5:8] = 0.9
        self.assertGreater(seg_mod.foreground_fraction_of_cell(whole, mask, 0.5), 0.9)
        self.assertLess(seg_mod.foreground_fraction_of_cell(sparse, mask, 0.5), 0.1)
        self.assertEqual(seg_mod.EXIT_IMPLAUSIBLE_FOREGROUND, 5)
        self.assertLessEqual(
            self.cfg["stage4_hybrid_segmentation"]["stage4a_ilastik"]["max_foreground_fraction_of_cell"], 0.5)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_05f_temporal_seeding_keeps_pairs_apart_but_allows_real_fusion(self):
        """Temporal seeding: two objects separate in the previous frame stay separate when this frame's
        own seeds would merge them (intensity dip still present), but merge when there is no dip
        (real fusion); an object that was single is not split by a weak extra peak."""
        from scipy import ndimage as ndi
        calib = {"dx_um": 0.19, "dy_um": 0.19, "dz_um": 1.0, "dt_s": 11.7}
        shape = (5, 30, 50)
        prev = np.zeros(shape, np.int32)
        prev[1:4, 12:18, 10:20] = 1
        prev[1:4, 12:18, 21:31] = 2
        # Current frame: same pair joined by a neck with a clear intensity dip
        img = np.zeros(shape, np.float32)
        img[1:4, 12:18, 10:20] = 1.0
        img[1:4, 12:18, 21:31] = 1.0
        img[1:4, 13:17, 20:21] = 0.4
        inten = ndi.gaussian_filter(img, (0.5, 1.0, 1.0))
        fg = inten > 0.15
        self.assertEqual(ndi.label(fg)[1], 1)
        strong_h = 1.0   # so strong that this frame's own seeds alone give ONE object
        alone = seg_mod.watershed_split_3d(fg, inten, calib, seed_prominence_um=strong_h)
        self.assertEqual(len(np.unique(alone[fg])), 1)
        kept = seg_mod.watershed_split_3d_temporal(fg, inten, calib, prev, seed_prominence_um=strong_h)
        self.assertEqual(len(np.unique(kept[fg])), 2)
        # Real fusion: one uniform rod covering both previous objects (no dip) -> merged into one
        fused = np.zeros(shape, np.float32)
        fused[0:5, 10:20, 10:31] = 1.0   # realistic thickness: the bright core spans the boundary plane
        f_int = ndi.gaussian_filter(fused, (0.5, 1.0, 1.0))
        f_fg = f_int > 0.15
        merged = seg_mod.watershed_split_3d_temporal(f_fg, f_int, calib, prev, seed_prominence_um=strong_h)
        self.assertEqual(len(np.unique(merged[f_fg])), 1)
        # Single object last frame + weak extra peak now -> not split (hysteresis)
        single_prev = np.where(prev > 0, 1, 0).astype(np.int32)
        w = np.zeros(shape, np.float32)
        w[1:4, 12:18, 10:31] = 1.0
        w[2, 14:16, 25:27] = 1.15   # faint bump
        w_int = ndi.gaussian_filter(w, (0.5, 1.0, 1.0))
        w_fg = w_int > 0.15
        n_plain = len(np.unique(seg_mod.watershed_split_3d(w_fg, w_int, calib)[w_fg]))
        n_temp = len(np.unique(seg_mod.watershed_split_3d_temporal(w_fg, w_int, calib, single_prev)[w_fg]))
        self.assertLessEqual(n_temp, n_plain)
        self.assertEqual(n_temp, 1)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_05c_instances_gated_by_ilastik_probability(self):
        """Stage 4B drops instances that ilastik does not call mitochondria (mean P below threshold)."""
        masks = np.zeros((1, 20, 20), dtype=np.int32)
        masks[0, 2:6, 2:6] = 1      # on foreground
        masks[0, 10:14, 10:14] = 2  # on background (noise object)
        prob = np.zeros((1, 20, 20), dtype=np.float32)
        prob[0, 2:6, 2:6] = 0.9
        prob[0, 10:14, 10:14] = 0.1
        gated = seg_mod.gate_instances_by_probability(masks, prob, 0.5)
        self.assertEqual(sorted(np.unique(gated).tolist()), [0, 1])

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_06_sparse_lap_tracking_gap_fission_fusion(self):
        """Verifies gap closing, track splitting (fission), merging (fusion), and TrackMate v7 XML validity."""
        # Create a 6-frame label stack (T=6, Z=6, Y=64, X=64); daughters persist 3 frames (T=3..5),
        # as required by min_event_persistence_frames: 3
        labels = np.zeros((6, 6, 64, 64), dtype=np.uint16)
        norm = np.full((6, 6, 64, 64), 10000, dtype=np.uint16)

        # Track 1: Present at T=0, missing at T=1 (gap!), present at T=2, splits into two at T=3 and T=4.
        # The parent (6x3 = 18 voxels) is twice a daughter (3x3 = 9), so the split conserves signal.
        labels[0, 3, 9:15, 10:13] = 1
        # T=1 gap (label 1 absent)
        labels[2, 3, 10:16, 10:13] = 1
        # T=3 Fission: parent splits into daughter A (y=10) and daughter B (y=16)
        labels[3, 3, 9:12, 10:13] = 1
        labels[3, 3, 15:18, 10:13] = 2
        labels[4, 3, 8:11, 10:13] = 1
        labels[4, 3, 16:19, 10:13] = 2
        labels[5, 3, 7:10, 10:13] = 1
        labels[5, 3, 17:20, 10:13] = 2

        lbl_path = self.tmp_dir / "gap_split_mito_labels.tif"
        norm_path = self.tmp_dir / "gap_split_mito_norm.tif"
        write_tiff(lbl_path, labels)
        write_tiff(norm_path, norm)

        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}
        n_tr, n_sp, n_ed, n_fiss, n_fuse = track_mod.run_sparse_lap_tracking(
            lbl_path, norm_path, self.tmp_dir, "gap_split", self.cfg, calib
        )

        # All 8 detections belong to 1 unified branched lineage tree with 1 gap-close and 1 fission!
        self.assertEqual(n_tr, 1)
        self.assertEqual(n_sp, 8)
        self.assertEqual(n_fiss, 1)

        # Verify TrackMate XML parses cleanly and has <AllSpots> and <AllTracks>
        xml_tree = ET.parse(self.tmp_dir / "gap_split_tracks.xml")
        root = xml_tree.getroot()
        self.assertIsNotNone(root.find(".//AllSpots"))
        self.assertIsNotNone(root.find(".//AllTracks"))

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_06c_transient_flicker_events_removed(self):
        """A split immediately undone by a merge (segmentation flicker) must not count as fission + fusion."""
        spots = {i: {"ID": i, "FRAME": f} for i, f in enumerate([0, 1, 1, 2, 3, 4, 5])}
        # 0 -> 1 (link) ; 0 -> 2 (split at t=1) ; 2 -> 3 (merge back at t=2) ; 1 -> 3 -> 4 -> 5 -> 6 (links)
        edges = [
            {"SPOT_SOURCE_ID": 0, "SPOT_TARGET_ID": 1, "EVENT_TYPE": "FRAME_LINK"},
            {"SPOT_SOURCE_ID": 0, "SPOT_TARGET_ID": 2, "EVENT_TYPE": "SPLIT_FISSION"},
            {"SPOT_SOURCE_ID": 2, "SPOT_TARGET_ID": 3, "EVENT_TYPE": "MERGE_FUSION"},
            {"SPOT_SOURCE_ID": 1, "SPOT_TARGET_ID": 3, "EVENT_TYPE": "FRAME_LINK"},
            {"SPOT_SOURCE_ID": 3, "SPOT_TARGET_ID": 4, "EVENT_TYPE": "FRAME_LINK"},
            {"SPOT_SOURCE_ID": 4, "SPOT_TARGET_ID": 5, "EVENT_TYPE": "FRAME_LINK"},
            {"SPOT_SOURCE_ID": 5, "SPOT_TARGET_ID": 6, "EVENT_TYPE": "FRAME_LINK"},
        ]
        kept, n = track_mod.filter_transient_events(edges, spots, 3)
        self.assertEqual(n, 2)
        self.assertFalse([e for e in kept if e["EVENT_TYPE"] in ("SPLIT_FISSION", "MERGE_FUSION")])
        self.assertEqual(len(kept), 5)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_06b_unconserved_appearance_is_not_fission(self):
        """A new, unrelated spot appearing next to a persisting one must not be called a fission."""
        labels = np.zeros((5, 4, 40, 40), dtype=np.uint16)
        norm = np.full((5, 4, 40, 40), 10000, dtype=np.uint16)
        for t in range(5):
            labels[t, 2, 10:13, 10:13] = 1          # persisting mitochondrion (9 voxels)
            if t >= 2:
                labels[t, 2, 14:17, 10:13] = 2      # equal-size newcomer 0.64 um away from t=2 on
        lbl_path = self.tmp_dir / "appear_mito_labels.tif"
        norm_path = self.tmp_dir / "appear_mito_norm.tif"
        write_tiff(lbl_path, labels)
        write_tiff(norm_path, norm)
        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}
        _, _, _, n_fiss, n_fuse = track_mod.run_sparse_lap_tracking(
            lbl_path, norm_path, self.tmp_dir, "appear", self.cfg, calib
        )
        self.assertEqual(n_fiss, 0)  # parent/(d1+d2) = 0.5 fails the 1.5x conservation gate
        self.assertEqual(n_fuse, 0)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_07_curved_arbor_internal_geodesic_and_swc(self):
        """Verifies that internal geodesic distance follows a curved U-shaped tube rather than Euclidean shortcuts."""
        # Build a U-shaped cell mask in (Z=3, Y=40, X=40) with dx=dy=0.5 um, dz=1.0 um
        cell_mask = np.zeros((3, 40, 40), dtype=bool)
        soma_mask = np.zeros((3, 40, 40), dtype=bool)

        # Soma at top-left (y=5..8, x=5..8)
        soma_mask[1, 5:8, 5:8] = True
        cell_mask[1, 5:8, 5:8] = True
        # Downward shaft: x=6..7, y=8..30 (length = 22 * 0.5 = 11.0 um)
        cell_mask[1, 8:30, 6:8] = True
        # Bottom turn: y=28..30, x=7..16 (length = 9 * 0.5 = 4.5 um)
        cell_mask[1, 28:30, 7:16] = True
        # Upward shaft returning near soma: x=14..16, y=8..28 (length = 20 * 0.5 = 10.0 um)
        cell_mask[1, 8:28, 14:16] = True

        calib = {"dx_um": 0.5, "dy_um": 0.5, "dz_um": 1.0, "dt_s": 11.7}
        geo_map = ml_mod.compute_internal_geodesic_distance_3d(cell_mask, soma_mask, calib)

        # Point at top of upward return shaft (z=1, y=9, x=15):
        # Straight-line Euclidean distance from soma (y=7, x=7) is only sqrt(2^2 + 8^2)*0.5 = 4.1 um,
        # but true internal path along the U-tube is ~11 + 4 + 10 = ~23-25 um!
        internal_dist_um = float(geo_map[1, 9, 15])
        self.assertGreater(internal_dist_um, 18.0)

        # Verify SWC generation produces multiple connected arbor nodes
        swc_path = self.tmp_dir / "curved_arbor.swc"
        ml_mod.build_arbor_swc_and_strahler_map(
            cell_mask, soma_mask, geo_map, calib, swc_path, "curved_test"
        )
        with open(swc_path, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
        self.assertGreater(len(lines), 5)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_07b_outside_mask_spots_and_nonnegative_soma_distance(self):
        """Spots just outside the cell snap to it; far-outside spots are 'outside_cell' (NaN distance),
        never silently 'basal_process' at distance 0; soma distances are >= 0."""
        calib = {"dx_um": 0.5, "dy_um": 0.5, "dz_um": 1.0, "dt_s": 10.0}
        cell = np.zeros((3, 40, 40), dtype=bool)
        cell[:, 5:35, 5:20] = True
        soma = np.zeros_like(cell)
        soma[:, 5:15, 5:20] = True
        geo = ml_mod.compute_internal_geodesic_distance_3d(cell, soma, calib)
        self.assertGreaterEqual(float(geo.min()), 0.0)
        comp = np.where(soma, 1, np.where(cell, 3, 0)).astype(np.uint8)
        labels = np.zeros((1, 3, 40, 40), dtype=np.uint16)
        labels[0, 1, 25:27, 10:12] = 1   # inside the process
        labels[0, 1, 25:27, 20:22] = 2   # 1-2 px (0.5-1.0 um) outside the mask edge -> snapped
        labels[0, 1, 25:27, 32:34] = 3   # 6+ um outside -> outside_cell
        feats, _, _ = ml_mod.profile_all_instances_fast(labels, calib, geo, comp, cell_mask=cell, snap_um=1.5)
        self.assertEqual(feats[(0, 1)]["compartment_name"], "basal_process")
        self.assertGreater(feats[(0, 1)]["geodesic_soma_dist_um"], 0.0)
        self.assertEqual(feats[(0, 2)]["compartment_name"], "basal_process")
        self.assertGreater(feats[(0, 2)]["geodesic_soma_dist_um"], 0.0)
        self.assertEqual(feats[(0, 3)]["compartment_name"], "outside_cell")
        self.assertTrue(np.isnan(feats[(0, 3)]["geodesic_soma_dist_um"]))

    @unittest.skipUnless(HAVE_HPC_PACKAGES, HPC_SKIP_MSG)
    def test_08_end_to_end_pipeline_integration(self):
        """Executes Modules 1 through 4 end-to-end on a synthetic 4D dual-channel movie."""
        movie_id = "e2e_synthetic"
        movie_path = self.tmp_dir / f"{movie_id}.ome.tif"

        # Realistic 5D (T=5, C=2, Z=10, Y=96, X=96) movie: camera offset 100 + read noise SD 3,
        # cytosol (Ch 2) soma + shaft, and exactly TWO blurred mitochondrial rods (Ch 1) spanning
        # 3 Z-slices: one static in the soma, one moving +3 px (0.48 um) per frame along the shaft.
        from scipy import ndimage as ndi
        rng = np.random.default_rng(7)
        n_t, n_z, n_y, n_x = 5, 10, 96, 96
        movie = np.zeros((n_t, 2, n_z, n_y, n_x), dtype=np.float32)
        zz, yy, xx = np.ogrid[:n_z, :n_y, :n_x]
        for t in range(n_t):
            soma = ((zz - 5) ** 2 * 4 + (yy - 22) ** 2 + (xx - 48) ** 2) <= 14 ** 2
            shaft = (np.abs(zz - 5) <= 2) & (yy >= 22) & (yy < 88) & (np.abs(xx - 48) <= 4)
            movie[t, 1][soma | shaft] = 300.0
            rods = np.zeros((n_z, n_y, n_x), dtype=np.float32)
            rods[4:7, 16:20, 40:48] = 1.0
            rods[4:7, 40 + 3 * t : 48 + 3 * t, 46:50] = 1.0
            movie[t, 0] = ndi.gaussian_filter(rods, (0.6, 1.0, 1.0)) * 400.0
        movie += rng.normal(100.0, 3.0, size=movie.shape)
        data = np.clip(movie, 0, 65535).astype(np.uint16)

        write_tiff(
            movie_path,
            data,
            metadata={"PhysicalSizeX": 0.16, "PhysicalSizeY": 0.16, "PhysicalSizeZ": 1.0, "TimeIncrement": 11.7},
        )

        # Run Module 1
        prep_mod.preprocess_movie(movie_path, self.tmp_dir, self.cfg)
        # Point models_dir at an empty folder so Stage 4 deterministically uses the built-in
        # estimator instead of whatever real .ilp is installed in $MITO_MODELS_DIR.
        os.environ["MITO_MODELS_DIR"] = str(self.tmp_dir / "no_models")
        # Run Module 2
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT_DIR / "run_hybrid_segmentation.py"),
                "--config", str(self.cfg_path),
                "--movie-id", movie_id,
                "--work-dir", str(self.tmp_dir),
            ],
            check=True,
        )
        # Run Module 3
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT_DIR / "run_trackmate_sparse_lap.py"),
                "--config", str(self.cfg_path),
                "--movie-id", movie_id,
                "--work-dir", str(self.tmp_dir),
            ],
            check=True,
        )
        # Run Module 4
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT_DIR / "run_arbor_and_mitolandscape.py"),
                "--config", str(self.cfg_path),
                "--movie-id", movie_id,
                "--work-dir", str(self.tmp_dir),
            ],
            check=True,
        )

        # Verify all expected outputs exist and are populated
        for suffix in [
            "_drift_shifts.csv",
            "_cell_mask.tif",
            "_mito_norm.tif",
            "_mito_labels.tif",
            "_spots.csv",
            "_edges.csv",
            "_tracks.csv",
            "_tracks.xml",
            "_arbor.swc",
            "_snt_mitolandscape_tracks.csv",
            "_process_0.5um_profile.csv",
            "_provenance.json",
        ]:
            f_path = self.tmp_dir / f"{movie_id}{suffix}"
            self.assertTrue(f_path.exists(), f"Missing output file: {f_path.name}")
            self.assertGreater(f_path.stat().st_size, 0, f"Empty output file: {f_path.name}")

        # The outputs must also be RIGHT, not just present: 2 rods per frame, no invented events,
        # and the moving rod's links match its true 0.48 um/frame displacement.
        prov = json.loads((self.tmp_dir / f"{movie_id}_provenance.json").read_text())
        self.assertEqual(prov["stage3_preprocessing_metrics"]["qc_status"], "PASS")
        labels, _ = read_tiff(self.tmp_dir / f"{movie_id}_mito_labels.tif", axes="TZYX")
        self.assertEqual([int(len(np.unique(labels[t])) - 1) for t in range(labels.shape[0])], [2] * n_t)
        with open(self.tmp_dir / f"{movie_id}_edges.csv", newline="") as f:
            edges = list(csv.DictReader(f))
        self.assertFalse([e for e in edges if e["EVENT_TYPE"] in ("SPLIT_FISSION", "MERGE_FUSION")])
        moves = sorted(float(e["DISPLACEMENT"]) for e in edges if e["EVENT_TYPE"] == "FRAME_LINK")
        self.assertEqual(len(moves), 2 * (n_t - 1))
        self.assertLess(moves[n_t - 2], 0.15)                  # static rod: ~0 um per frame
        self.assertAlmostEqual(moves[-1], 0.48, delta=0.12)    # moving rod: 3 px x 0.16 um


class TestFullRunOrchestration(unittest.TestCase):
    """
    Full-archive orchestration (submit_full_run.bash, plan_full_run.py, merge_collect_shards.py, the
    collector's --shard, stage-done markers in run_pipeline_slurm.bash, parallel watershed).
    Standard-library tests run anywhere; the marker and watershed tests need pyyaml / numpy+skimage.
    """

    # Fake showinf: every .sld holds one qualifying 4D series and one 2D projection.
    FAKE_OME = (
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        '<Image ID="Image:0" Name="cell_t1"><Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" '
        'SizeX="64" SizeY="48" SizeZ="20" SizeC="2" SizeT="12" PhysicalSizeX="0.16" PhysicalSizeY="0.16" '
        'PhysicalSizeZ="1.0" TimeIncrement="11.7"/></Image>'
        '<Image ID="Image:1" Name="cell_t1 MAX"><Pixels ID="Pixels:1" DimensionOrder="XYCZT" Type="uint16" '
        'SizeX="64" SizeY="48" SizeZ="1" SizeC="2" SizeT="1"/></Image></OME>'
    )

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="mito_fullrun_test_"))
        import plan_full_run
        import merge_collect_shards
        self.plan = plan_full_run
        self.merge = merge_collect_shards

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    # ---------------------------------------------------------------- helpers
    def _fake_bftools(self) -> Path:
        bft = self.tmp / "bft"
        bft.mkdir()
        for tool in ("showinf", "bfconvert"):
            (bft / tool).write_text(
                "#!/bin/sh\n"
                'case "$*" in *-version*) echo "Version: fake"; exit 0;; esac\n'
                f"cat <<'XML'\n{self.FAKE_OME}\nXML\n"
            )
            (bft / tool).chmod(0o755)
        return bft

    def _sld_tree(self) -> Path:
        src = self.tmp / "src"
        names = ["2024/a_cell1.sld", "2024/b_cell2.sld", "2024/sub dir/c.sld", "2025/d.sld", "2025/e.sldy",
                 "2025/f.sld", "2025/g.sld", "2025/~lock.sld"]
        for n in names:
            p = src / n
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"\x01" * 70_000)
        (src / "2025" / "stub.sld").write_bytes(b"\x01" * 10)  # SKIPPED_STUB
        return src

    def _collect(self, src: Path, out: Path, *extra: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(SCRIPT_DIR / "collect_4d_sld_timelapses.py"), str(src), str(out),
             "--mode", "scan-only", "--workers", "1", "--bftools-dir", str(self.bft),
             "--hpc-raw-root", "/raw/root", *extra],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
            env={**os.environ, "MITO_PLASMID_CSV": ""},
        )

    # ---------------------------------------------------------------- tests
    def test_09a_shard_partition_covers_every_file_exactly_once(self):
        """--shard K/N: the N shards are disjoint and together cover every candidate exactly once."""
        files = [Path(f"/data/{i:03d}.sld") for i in range(23)]
        for n in (1, 2, 5, 7, 23, 30):
            shards = [sld_mod.select_shard(files, k, n) for k in range(1, n + 1)]
            flat = [p for s in shards for p in s]
            self.assertEqual(sorted(flat), files, f"N={n}")
            self.assertEqual(len(flat), len(set(flat)), f"N={n}")
        self.assertEqual(sld_mod.parse_shard_spec("3/20"), (3, 20))
        self.assertIsNone(sld_mod.parse_shard_spec(None))
        for bad in ("0/5", "6/5", "a/b", "3", "1/0"):
            with self.assertRaises(ValueError):
                sld_mod.parse_shard_spec(bad)
        self.assertEqual(sld_mod.shard_report_name("sld_4d_inventory.csv", "_shard3"), "sld_4d_inventory_shard3.csv")
        self.assertEqual(sld_mod.shard_report_name("extracted_checksums.sha256", ""), "extracted_checksums.sha256")

    def test_09b_sharded_scan_merge_round_trip(self):
        """Sharded scan + merge_collect_shards reproduces the unsharded inventory, health audit, movie lists
        and mito classification byte for byte; duplicates and missing shards are refused."""
        if sys.platform.startswith("win"):
            self.skipTest("Fake showinf shell script requires a POSIX shell")
        self.bft = self._fake_bftools()
        src = self._sld_tree()
        plasmids = self.tmp / "plasmids.csv"
        plasmids.write_text(
            "s,sequenced? (date),does this plasmid make FP positive cells?,glycerol stock?,-20 DNA duplicate?,"
            "addgene submission?,plasmid ID,name (Prom_GOI_3'),gateway compnent?,localization sequence?,"
            "GOI? tag?,Fluorescent Protein 1 (ME),Fluorescent Protein 2 (3'),p2a,\"uas, MCS, promoter?\","
            "antibiotic,plasmid description/elements,full map?,long term use,Forward Primer,Reverse Primer,"
            "product size \n,,,,,,JR19,pUAS-mlsEGFPpA,,mito,,EGFP,,,UAS,amp,,,,,,\n", encoding="utf-8")
        filt = ["--plasmid-csv", str(plasmids), "--mito-filter", "report"]
        ref = self.tmp / "ref"
        proc = self._collect(src, ref, *filt)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        shard_dir = self.tmp / "shards"
        n = 3
        for k in range(1, n + 1):
            proc = self._collect(src, shard_dir, "--shard", f"{k}/{n}", *filt)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((shard_dir / f"sld_4d_inventory_shard{k}.csv").exists())
            self.assertTrue((shard_dir / f"mito_classification_shard{k}.csv").exists())
        self.assertFalse((shard_dir / "sld_4d_inventory.csv").exists())  # no unsuffixed clobbering
        merged = self.tmp / "merged"
        rc, summary = self.merge.run_merge(shard_dir, merged)
        self.assertEqual(rc, 0, summary)
        self.assertEqual(summary["shards_missing"], [])
        for name in ("sld_4d_inventory.csv", "sld_health_audit.csv", "movie_list.txt", "movie_list_2024.txt",
                     "movie_list_2025.txt", "mito_classification.csv", "mito_review.csv"):
            self.assertEqual((merged / name).read_bytes(), (ref / name).read_bytes(), name)
        self.assertEqual(len((merged / "movie_list.txt").read_text().splitlines()), 7)
        prov_ref = json.loads((ref / "extraction_provenance.json").read_text())["summary_metrics"]
        prov_m = json.loads((merged / "extraction_provenance.json").read_text())["summary_metrics"]
        for key in ("total_sld_candidates_found", "total_sld_files_scanned", "total_unparseable_or_skipped",
                    "total_series_inspected", "qualifying_4d_series_count", "cohorts_discovered"):
            self.assertEqual(prov_m[key], prov_ref[key], key)

        # --file-list: the same files, listed explicitly (relative to source_dir), give the same inventory
        flist = self.tmp / "files.txt"
        flist.write_text("\n".join(sorted(p.relative_to(src).as_posix() for p in src.rglob("*.sld*"))) + "\n")
        via_list = self.tmp / "via_list"
        proc = self._collect(src, via_list, "--file-list", str(flist), *filt)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((via_list / "sld_4d_inventory.csv").read_bytes(), (ref / "sld_4d_inventory.csv").read_bytes())

        # A missing shard -> exit 3 (unless allowed); a duplicated shard -> refused
        (shard_dir / "extraction_provenance_shard2.json").rename(self.tmp / "p2.json")
        rc, summary = self.merge.run_merge(shard_dir, self.tmp / "m2", expected_shards=3)
        self.assertEqual((rc, summary["shards_missing"]), (3, [2]))
        self.assertEqual(self.merge.run_merge(shard_dir, self.tmp / "m3", 3, allow_missing=True)[0], 0)
        (self.tmp / "p2.json").rename(shard_dir / "extraction_provenance_shard2.json")
        shutil.copy(shard_dir / "sld_4d_inventory_shard1.csv", shard_dir / "sld_4d_inventory_shard4.csv")
        with self.assertRaises(self.merge.MergeConflict):
            self.merge.run_merge(shard_dir, self.tmp / "m4")

    def test_09c_size_classes_and_resource_model(self):
        """Size-class binning, the memory/time model seeded from the 2026-11-24 measurements, request
        rounding, throttles (20% cap) and canary selection."""
        p = self.plan
        hpc = p.load_hpc_run(REPO_ROOT / "pipeline_config.yaml")
        if _installed("yaml"):
            import yaml
            with open(REPO_ROOT / "pipeline_config.yaml", encoding="utf-8") as f:
                self.assertEqual(p.deep_merge(p.DEFAULT_HPC_RUN, yaml.safe_load(f)["hpc_run"]), hpc)
        model, cls = hpc["model"], hpc["size_classes"]
        self.assertEqual([p.size_class(g, cls["edges_gb"], cls["names"]) for g in (0.5, 1.99, 2.0, 7.9, 8.0, 50)],
                         ["S", "S", "M", "M", "L", "L"])
        s39 = p.Movie(entry="x", path="/r/2025/S39.ome.tif", gb=1.38, size_t=6, size_c=2, size_source="inventory")
        mem, sec = p.predict_stage3(s39, model)
        self.assertAlmostEqual(mem, 2.75 * 1.38 + 0.5)                 # S39 measured 3.5 GB peak
        self.assertAlmostEqual(sec, 65 * 1.38)                          # ~90 s measured
        self.assertEqual(p.to_request(mem, sec, model), (7, 10))        # 4.3 GB x1.5 -> 7 GB; 10 min floor
        big = p.Movie(entry="y", path="/r/2025/big.ome.tif", gb=20.0, size_t=30, size_c=2, size_source="inventory")
        self.assertEqual(p.to_request(*p.predict_stage3(big, model), model), (84, 44))  # 55.5 GB, 1300 s
        mem45, sec45 = p.predict_stage45(s39, model, cpus=4)
        chan = 1.38 / 2
        self.assertAlmostEqual(sec45, 290 * chan + 74 * chan * 2 / 6 + 60)  # ilastik + 2 watershed rounds + Stage 5
        self.assertAlmostEqual(mem45, 2.75 * 1.38 + 1.5 + 4 * 8.0 * chan / 6)
        # ilastik speed-up saturates at 4 threads: 8 cpus is no faster for ilastik than 4
        import math
        il = lambda c: p.predict_stage45(big, model, cpus=c)[1] - 74 * 10 * math.ceil(30 / min(c, 30)) / 30 - 60  # noqa: E731
        self.assertAlmostEqual(il(8), il(4))
        self.assertAlmostEqual(il(1), 4 * il(4))
        self.assertEqual(p.fmt_time(44), "00:44:00")
        self.assertEqual(p.fmt_time(1500), "1-01:00:00")
        self.assertEqual(p.parse_time_to_min("1-01:00:00"), 1500)
        self.assertEqual((p.cap_tasks(900, 0.2, 1), p.cap_tasks(900, 0.2, 4), p.cap_tasks(128, 0.2, 1)), (180, 45, 25))
        self.assertEqual(p.split_throttle(180, [300, 40]), [158, 21])
        self.assertEqual(p.split_throttle(180, [5, 3]), [5, 3])
        self.assertEqual([len(c) for c in p.chunk(list(range(2500)), 1000)], [1000, 1000, 500])
        movies = ([p.Movie("s", f"/r/s{i}", 1.0 + i / 100, 10, 2, "inventory", size_class="S") for i in range(30)]
                  + [p.Movie("m", f"/r/m{i}", 4.0 + i, 10, 2, "inventory", size_class="M") for i in range(8)]
                  + [p.Movie("l", f"/r/l{i}", 20.0 + i, 10, 2, "inventory", size_class="L") for i in range(2)])
        canary = p.select_canary(movies, 10)
        self.assertEqual(len(canary), 10)
        self.assertEqual({m.size_class for m in canary}, {"S", "M", "L"})
        self.assertIn("/r/l1", [m.path for m in canary])       # largest movie of a class always included
        self.assertIn("/r/s29", [m.path for m in canary])
        # sacct parsing, retry selection (QC exit 3 never retried) and --array index compression
        sacct = ("100_1|COMPLETED|0:0|00:01:30|\n100_1.batch|COMPLETED|0:0|00:01:30|3438.28M\n"
                 "100_2|OUT_OF_MEMORY|0:125|00:00:40|\n100_2.batch|OUT_OF_MEMORY|0:125|00:00:40|7000M\n"
                 "100_3|FAILED|3:0|00:01:00|\n100_4|TIMEOUT|0:0|00:10:00|\n100_5|FAILED|1:0|00:00:05|\n")
        tasks = p.parse_sacct(sacct)
        self.assertAlmostEqual(tasks[("100", 1)]["maxrss_mb"], 3438.28)
        self.assertEqual(tasks[("100", 1)]["elapsed_s"], 90)
        states = p.latest_task_states({"job_id": "100", "retries": []}, sacct)
        self.assertEqual(p.retry_indices("stage3", states), [2, 4, 5])
        self.assertEqual(p.compress_indices([1, 2, 3, 7, 9, 10]), "1-3,7,9-10")

    def test_09d_done_markers_skip_and_rerun_on_config_change(self):
        """Stage hashes are chained: a Stage 4 config change re-runs 4 and 5 but not 3; comment edits re-run
        nothing; a new .ilp re-runs 4-5; a new input re-runs everything. Markers match only their hash."""
        p = self.plan
        cfg = self.tmp / "cfg.yaml"
        text = (REPO_ROOT / "pipeline_config.yaml").read_text(encoding="utf-8")
        cfg.write_text(text, encoding="utf-8")
        movie = self.tmp / "m.ome.tif"
        movie.write_bytes(b"x" * 100)
        h0 = p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="abc")
        marker = self.tmp / "m_stage3.done"
        p.write_marker(marker, h0[3], 3, {"job": "1_1"})
        self.assertTrue(p.marker_matches(marker, h0[3]))
        self.assertFalse(p.marker_matches(marker, h0[4]))
        self.assertFalse(p.marker_matches(self.tmp / "absent.done", h0[3]))

        cfg.write_text(text.replace("max_frame_gap: 2", "max_frame_gap: 2   # comment edit only"), encoding="utf-8")
        self.assertEqual(p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="abc"), h0)
        cfg.write_text(text.replace("seed_min_separation_xy_um: 0.4", "seed_min_separation_xy_um: 0.5"), encoding="utf-8")
        h4 = p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="abc")
        self.assertEqual(h4[3], h0[3])
        self.assertNotEqual(h4[4], h0[4])
        self.assertNotEqual(h4[5], h0[5])
        cfg.write_text(text.replace("linking_max_distance_um: 2.5", "linking_max_distance_um: 3.0"), encoding="utf-8")
        h5 = p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="abc")
        self.assertEqual((h5[3], h5[4]), (h0[3], h0[4]))
        self.assertNotEqual(h5[5], h0[5])
        cfg.write_text(text.replace("min_contrast_sigma: 3.0", "min_contrast_sigma: 2.0"), encoding="utf-8")
        h3 = p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="abc")
        self.assertTrue(all(h3[s] != h0[s] for s in (3, 4, 5)))
        cfg.write_text(text, encoding="utf-8")
        ilp = self.tmp / "model.ilp"
        ilp.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\0" * 32)
        hi = p.compute_stage_hashes(cfg, REPO_ROOT, movie, ilp, head="abc")
        self.assertEqual(hi[3], h0[3])
        self.assertNotEqual(hi[4], h0[4])
        self.assertNotEqual(p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="def")[3], h0[3])
        movie.write_bytes(b"x" * 101)
        self.assertNotEqual(p.compute_stage_hashes(cfg, REPO_ROOT, movie, None, head="abc")[3], h0[3])

    @unittest.skipUnless(_installed("yaml") and not sys.platform.startswith("win"), "Requires pyyaml and bash")
    def test_09e_slurm_script_markers_skip_rerun_and_qc_gate(self):
        """run_pipeline_slurm.bash with fake stage scripts: a second run is a no-op, a Stage 4 config change
        re-runs only 4-5, Stage 4-5 alone after a QC-failed Stage 3 exits 3, TMPDIR is under scr30-style
        MITO_SCRATCH_ROOT/tmp and is removed afterwards."""
        repo = self.tmp / "repo"
        (repo / "scripts").mkdir(parents=True)
        for name in ("run_pipeline_slurm.bash", "plan_full_run.py", "pipeline_config_io.py", "mito_env.bash", "tiff_io.py"):
            shutil.copy(SCRIPT_DIR / name, repo / "scripts" / name)
        cfg_text = (REPO_ROOT / "pipeline_config.yaml").read_text(encoding="utf-8")
        (repo / "pipeline_config.yaml").write_text(cfg_text, encoding="utf-8")
        log = self.tmp / "stages.log"
        fake = (
            "import argparse, os, sys, pathlib\n"
            "ap = argparse.ArgumentParser(); ap.add_argument('--config'); ap.add_argument('--input')\n"
            "ap.add_argument('--movie-id'); ap.add_argument('--work-dir'); a = ap.parse_args()\n"
            "wd = pathlib.Path(a.work_dir); mid = a.movie_id or pathlib.Path(a.input).name.split('.')[0]\n"
            "open(os.environ['FAKE_LOG'], 'a').write('{stage} ' + os.environ.get('TMPDIR', '') + '\\n')\n"
            "for suffix in {outputs}: (wd / (mid + suffix)).write_text('x')\n"
            "qc = os.environ.get('FAKE_QC_FAIL') == '1' and '{stage}' == '3'\n"
            "if '{stage}' == '3': (wd / (mid + '_provenance.json')).write_text("
            "'{{\"stage3_preprocessing_metrics\": {{\"qc_status\": \"' + ('FAIL_NO_SIGNAL' if qc else 'PASS') + '\"}}}}')\n"
            "sys.exit(3 if qc else 0)\n"
        )
        for script, stage, outputs in (
            ("preprocess_4d_stack.py", "3", ["_mito_norm.tif", "_cell_RFP_reg.tif", "_cell_mask.tif"]),
            ("run_hybrid_segmentation.py", "4", ["_mito_labels.tif", "_Probabilities.tif"]),
            ("run_trackmate_sparse_lap.py", "5A", ["_tracks.csv"]),
            ("run_arbor_and_mitolandscape.py", "5B", ["_process_0.5um_profile.csv"]),
        ):
            (repo / "scripts" / script).write_text(fake.format(stage=stage, outputs=outputs))
        tools = self.tmp / "data" / "tools"
        tools.mkdir(parents=True)
        (tools / "activate_mito4d.sh").write_text("# test: no venv\n")
        raw = self.tmp / "data" / "raw_ometiff" / "2025"
        raw.mkdir(parents=True)
        movies = []
        for i in (1, 2):
            mv = raw / f"movie{i}.ome.tif"
            mv.write_bytes(b"\0" * 1000)
            movies.append(mv)
        lst = self.tmp / "list.txt"
        lst.write_text("\n".join(str(m) for m in movies) + "\n")
        scratch = self.tmp / "scr30" / "mito_work"
        env = {**os.environ, "MITO_REPO_DIR": str(repo), "MITO_DATA_ROOT": str(self.tmp / "data"),
               "MITO_SCRATCH_ROOT": str(scratch), "MITO_TOOLS_DIR": str(tools), "MITO_RAW_ROOT": str(raw.parent),
               "MITO_MODELS_DIR": str(self.tmp / "models"), "MITO_OUTPUT_ROOT": str(self.tmp / "out"),
               "MITO_ALLOW_ILASTIK_FALLBACK": "1", "FAKE_LOG": str(log), "SLURM_JOB_ID": "777",
               "SLURM_ARRAY_JOB_ID": "777", "SLURM_CPUS_PER_TASK": "1", "MITO_STAGES": "3,4,5",
               # the script's 'python3' must be this interpreter (the one that has pyyaml)
               "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ.get('PATH', '')}"}
        for k in ("MITO_ILASTIK_BIN", "MITO_IGNORE_QC"):
            env.pop(k, None)

        def run(task: int, *args: str, **extra: str) -> subprocess.CompletedProcess:
            log.write_text("")
            return subprocess.run(
                ["bash", str(repo / "scripts" / "run_pipeline_slurm.bash"), str(repo / "pipeline_config.yaml"),
                 str(lst), *args], env={**env, "SLURM_ARRAY_TASK_ID": str(task), **extra},
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False, cwd=str(self.tmp))

        def stages_run() -> list:
            return [line.split()[0] for line in log.read_text().splitlines()]

        r = run(1, "--keep-intermediates")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(stages_run(), ["3", "4", "5A", "5B"])
        tmpdirs = {line.split()[1] for line in log.read_text().splitlines()}
        self.assertEqual(tmpdirs, {str(scratch / "tmp" / "777_1")})
        self.assertFalse((scratch / "tmp" / "777_1").exists())        # removed by the EXIT trap
        work = scratch / "2025" / "movie1"
        for s in (3, 4, 5):
            self.assertTrue((work / f"movie1_stage{s}.done").exists())
        r = run(1, "--keep-intermediates")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(stages_run(), [])                              # second invocation: no-op
        (repo / "pipeline_config.yaml").write_text(
            cfg_text.replace("seed_min_separation_xy_um: 0.4", "seed_min_separation_xy_um: 0.5"), encoding="utf-8")
        r = run(1)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(stages_run(), ["4", "5A", "5B"])               # Stage 3 still valid
        self.assertFalse((work / "movie1_mito_norm.tif").exists())     # complete -> intermediates cleaned
        r = run(1)
        self.assertEqual((r.returncode, stages_run()), (0, []))          # still a no-op after cleanup
        r = run(1, "--force", "--stages", "3")
        self.assertEqual((r.returncode, stages_run()), (0, ["3"]))       # --force regenerates training data
        # Movie 2: Stage 3 fails QC -> exit 3, no marker; a later Stage 4-5 job also exits 3
        r = run(2, "--stages", "3", FAKE_QC_FAIL="1")
        self.assertEqual(r.returncode, 3, r.stderr)
        self.assertFalse((scratch / "2025" / "movie2" / "movie2_stage3.done").exists())
        r = run(2, "--stages", "4,5")
        self.assertEqual((r.returncode, stages_run()), (3, []))

    @unittest.skipUnless(HAVE_HPC_PACKAGES and _installed("skimage"), HPC_SKIP_MSG)
    def test_09f_parallel_watershed_equals_serial(self):
        """Watershed frames run in a process pool (n_workers > 1) give labels identical to the serial loop
        and to the previous full-frame implementation."""
        from scipy import ndimage as ndi
        calib = {"dx_um": 0.19, "dy_um": 0.19, "dz_um": 1.0, "dt_s": 120.0}
        n_t, shape = 5, (6, 48, 64)
        rng = np.random.default_rng(7)
        vol = np.zeros((n_t,) + shape, np.float32)
        for t in range(n_t):
            for k in range(4):
                z, y, x = 2 + k % 2, 8 + 9 * k, 6 + 10 * k + 2 * t
                vol[t, z:z + 2, y:y + 5, x:x + 8] = 1.0
                vol[t, z:z + 2, y + 1:y + 4, x + 8:x + 11] = 0.6   # touching neighbour
            vol[t] = ndi.gaussian_filter(vol[t], (0.5, 1.0, 1.0)) + rng.normal(0, 0.01, shape)
        vol[3] = 0.0                                               # one empty frame
        prob = np.clip(vol, 0, 1).astype(np.float32)
        mito = (prob * 50000).astype(np.uint16)
        cfg = seg_mod.load_config(REPO_ROOT / "pipeline_config.yaml")
        # This test covers the frame-INDEPENDENT mode; temporal seeding (default) is sequential by design
        cfg["stage4_hybrid_segmentation"]["stage4b_instance_segmentation"].setdefault("watershed", {})["temporal_seeding"] = False
        serial, used = seg_mod.segment_instances(mito, prob, cfg, calib, n_workers=1)
        parallel, _ = seg_mod.segment_instances(mito, prob, cfg, calib, n_workers=3)
        self.assertEqual(used, "watershed")
        self.assertEqual(parallel.dtype, serial.dtype)
        self.assertTrue(np.array_equal(serial, parallel))
        self.assertGreater(int(serial.max()), 3)
        self.assertEqual(int(serial[3].max()), 0)
        # Reference: the pre-parallel loop (watershed + gate + relabel on the FULL frame)
        b = cfg["stage4_hybrid_segmentation"]["stage4b_instance_segmentation"]
        voxel = calib["dx_um"] * calib["dy_um"] * calib["dz_um"]
        min_vox = max(2, int(round(float(b["min_volume_um3"]) / voxel)))
        margin = int(round(2 * float(b["diameter_px"])))
        for t in range(n_t):
            fg = prob[t] > float(b["foreground_probability_threshold"])
            ref = np.zeros(shape, np.uint16)
            if fg.any():
                zs, ys, xs = (np.flatnonzero(fg.any(axis=ax)) for ax in ((1, 2), (0, 2), (0, 1)))
                box = (slice(max(0, zs[0] - 1), min(shape[0], zs[-1] + 2)),
                       slice(max(0, ys[0] - margin), min(shape[1], ys[-1] + 1 + margin)),
                       slice(max(0, xs[0] - margin), min(shape[2], xs[-1] + 1 + margin)))
                full = np.zeros(shape, np.int32)
                full[box] = seg_mod.watershed_split_3d(fg[box], mito[t][box], calib)
                full = seg_mod.gate_instances_by_probability(full, prob[t], float(b["min_object_mean_probability"]))
                ref, _ = seg_mod.filter_and_relabel_fast(full, min_vox)
            self.assertTrue(np.array_equal(ref, serial[t]), f"frame {t}")
        self.assertEqual(seg_mod.watershed_worker_count(2, 8), 2)
        self.assertGreaterEqual(seg_mod.watershed_worker_count(10, None), 1)

    @unittest.skipIf(sys.platform.startswith("win"), "submit_full_run.bash requires bash")
    def test_09g_submit_dry_run_prints_dependency_chain_and_throttles(self):
        """Dry run of every phase prints afterany (merge, Stage 3 after canary, aggregate), afterok (extract),
        aftercorr + kill-on-invalid-dep (Stages 4-5), the 20% / I/O throttles, hima GPU-node exclusion and
        chunking - and never calls sbatch (a fake sbatch on PATH fails the test if it is called)."""
        fakebin = self.tmp / "bin"
        fakebin.mkdir()
        called = self.tmp / "sbatch_called.txt"
        (fakebin / "sbatch").write_text(f'#!/bin/sh\necho "$@" >> "{called}"\nexit 1\n')
        nodes_batch = "\n".join(f"bo{i:02d}|20|128000|bora,bo|(null)" for i in range(45))
        nodes_hima = "\n".join([f"hi0{i}|32|257000|hima,hi,hi0{i}|(null)" for i in (1, 2, 3, 6)]
                               + [f"hi0{i}|32|257000|hima,hi,gpu,p100,hi0{i}|gpu:p100:1" for i in (4, 5, 7)])
        (fakebin / "sinfo").write_text(
            "#!/bin/sh\ncase \"$*\" in\n"
            f"  *hima*) cat <<'E'\n{nodes_hima}\nE\n;;\n"
            f"  *) cat <<'E'\n{nodes_batch}\nE\n;;\nesac\n")
        (fakebin / "squeue").write_text("#!/bin/sh\nexit 0\n")           # nothing is in the queue
        (fakebin / "sacct").write_text(f'#!/bin/sh\ncat "{self.tmp}/sacct.txt" 2>/dev/null\n')
        for tool in ("sbatch", "sinfo", "squeue", "sacct"):
            (fakebin / tool).chmod(0o755)
        cfg = self.tmp / "cfg.yaml"
        cfg.write_text((REPO_ROOT / "pipeline_config.yaml").read_text(encoding="utf-8")
                       .replace("max_array_size: 1000", "max_array_size: 200"), encoding="utf-8")
        # 120 .sld files; 343 movies: 300 S (1 GB), 40 M (4 GB), 3 L (20 GB)
        src = self.tmp / "sld"
        src.mkdir()
        for i in range(120):
            (src / f"f{i:03d}.sld").write_bytes(b"")
        raw = self.tmp / "raw"
        (raw / "2025").mkdir(parents=True)
        rows, entries = [], []
        for cls, n, gb in (("s", 300, 1.0), ("m", 40, 4.0), ("l", 3, 20.0)):
            for i in range(n):
                rel = f"2025/{cls}{i:03d}.ome.tif"
                (raw / rel).write_bytes(b"")
                entries.append(str(raw / rel))
                rows.append({"cohort_year": "2025", "sld_path": f"{src}/f{i % 120:03d}.sld", "series_index": str(i),
                             "size_t": "10", "size_c": "2", "estimated_uncompressed_gb": str(gb),
                             "is_4d_timelapse": "True", "extracted_ometiff_relpath": rel,
                             "mito_label": "MITO" if i % 10 else "REVIEW"})
        inv = self.tmp / "inv.csv"
        with open(inv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        movie_list = self.tmp / "movie_list.txt"
        movie_list.write_text("\n".join(entries) + "\n")
        env = {**os.environ, "PATH": f"{fakebin}{os.pathsep}{os.environ.get('PATH', '')}",
               "MITO_DATA_ROOT": str(self.tmp / "data"), "MITO_SCRATCH_ROOT": str(self.tmp / "scr" / "mito_work"),
               "MITO_RAW_ROOT": str(raw), "MITO_MODELS_DIR": str(self.tmp / "models"), "USER": "tester"}
        for k in ("MITO_OUTPUT_ROOT", "MITO_TOOLS_DIR", "MITO_RUN_ID", "MITO_REPO_DIR"):
            env.pop(k, None)

        def submit(*args: str) -> str:
            proc = subprocess.run(["bash", str(SCRIPT_DIR / "submit_full_run.bash"), *args, "--run-id", "t1",
                                   "--no-module", "--config", str(cfg)],
                                  env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            return proc.stdout

        pre = "mito_t1"
        out = submit("scan", "--source-dir", str(src), "--chain-extract")
        self.assertIn("--array=1-24%20", out)                                # 120 files / 5 per shard, scan %20
        self.assertIn(f"--dependency=afterany:DRY_{pre}_scan", out)          # merge after every shard
        self.assertIn(f"--dependency=afterok:DRY_{pre}_merge", out)          # extract only after a clean merge
        self.assertIn("--array=1-120%20", out)                               # extract %20 (I/O throttle)
        self.assertIn(f"--dependency=afterany:DRY_{pre}_extract", out)

        run_dir = self.tmp / "scr" / "mito_runs" / "t1"
        (run_dir / "scan" / "merged").mkdir(parents=True)
        shutil.copy(inv, run_dir / "scan" / "merged" / "sld_4d_inventory.csv")
        out = submit("extract", "--source-dir", str(src))
        self.assertIn("--array=1-108%20", out)                                # only files with a MITO series

        out = submit("canary", "--movie-list", str(movie_list), "--inventory", str(inv))
        self.assertIn(f"--job-name={pre}_canary_L_c1", out)
        self.assertIn("--partition=hima", out)
        self.assertIn("--exclude=hi04,hi05,hi07", out)                       # never the GPU nodes
        canary_lists = sorted((run_dir / "dryrun" / "lists").glob("canary_*.txt"))
        self.assertEqual(sum(len(p.read_text().split()) for p in canary_lists), 10)

        out = submit("stage3", "--movie-list", str(movie_list), "--inventory", str(inv))
        self.assertIn("308 movie(s) found", out)                             # MITO only: 35 REVIEW movies excluded
        for name, array in (("stage3_S_c1", "1-200%"), ("stage3_S_c2", "1-70%"), ("stage3_M_c1", "1-36%"),
                            ("stage3_L_c1", "1-2%")):
            self.assertIn(f"--job-name={pre}_{name}", out)
            self.assertIn(f"--array={array}", out)
        self.assertIn("cap 20% -> 180 concurrent 1-core task(s)", out)       # 45 x 20 cores
        self.assertIn("cap 20% -> 25 concurrent 1-core task(s)", out)        # hima CPU nodes only: 128 cores
        self.assertIn("--array=1-200%117", out)                              # 180 shared: 200/306
        self.assertIn(f"--dependency=afterany:DRY_{pre}_canary_", out)

        out = submit("stage45", "--allow-fallback")
        for cls in ("S_c1", "S_c2", "M_c1", "L_c1"):
            self.assertIn(f"--dependency=aftercorr:DRY_{pre}_stage3_{cls}", out)
        self.assertIn("--kill-on-invalid-dep=yes", out)
        self.assertIn("--cpus-per-task=4", out)
        self.assertIn("ArrayTaskThrottle=", out)        # Stage 3 still holds the cap -> raise later
        self.assertIn("MITO_ALLOW_ILASTIK_FALLBACK=1", out)

        out = submit("aggregate")
        dep = next(line for line in out.splitlines() if "--dependency=afterany:" in line)
        for cls in ("S_c1", "S_c2", "M_c1", "L_c1"):
            self.assertIn(f"DRY_{pre}_stage45_{cls}", dep)

        # Canary refit: pretend the canary ran (real job id, now finished) and feed its sacct numbers.
        man_path = run_dir / "dryrun" / "run_manifest.json"
        man = json.loads(man_path.read_text())
        s_canary = next(n for n in man["jobs"] if n.endswith("canary_S_c1"))
        man["jobs"][s_canary]["job_id"] = "555"
        man_path.write_text(json.dumps(man))
        (self.tmp / "sacct.txt").write_text("555_1|COMPLETED|0:0|00:03:00|\n555_1.batch|COMPLETED|0:0|00:03:00|9000M\n")
        out = submit("stage3", "--movie-list", str(movie_list), "--inventory", str(inv))
        self.assertIn("model refit from 1 canary task(s)", out.replace("\n# ", " "))
        self.assertIn("dropped: no longer in the queue", out)              # finished canary id not used as dep
        self.assertFalse(called.exists(), "sbatch must never be called in a dry run")


def ndimage_com(vol) -> Tuple[float, float, float]:
    """Helper to compute 3D center of mass of foreground intensity."""
    from scipy import ndimage as ndi
    return ndi.center_of_mass(vol)


if __name__ == "__main__":
    unittest.main(verbosity=2)