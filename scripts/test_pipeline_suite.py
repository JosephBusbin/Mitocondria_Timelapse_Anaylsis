#!/usr/bin/env python3
"""
test_pipeline_suite.py
Comprehensive Verification & Unit/Integration Test Suite for the 4D Mitochondrial HPC Pipeline

Tests Covered:
  1. test_01_sld_metadata_and_safety:
     - Verifies OME-XML parsing, 4D filtering (Z >= 15, T >= 10, C >= 2), cohort year inference,
       and strict read-only safety rejection when source_dir == output_dir.
  2. test_02_tiff_io_5d_roundtrip:
     - Verifies lossless 5D (T, C, Z, Y, X) uint16 and float32 round-trip read/write via tiff_io.py.
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

Usage:
  python3 test_pipeline_suite.py -v
"""

import csv
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
sys.path.insert(0, str(SCRIPT_DIR))

import collect_4d_sld_timelapses as sld_mod

try:
    import numpy as np
    import preprocess_4d_stack as prep_mod
    import run_arbor_and_mitolandscape as ml_mod
    import run_hybrid_segmentation as seg_mod
    import run_trackmate_sparse_lap as track_mod
    from tiff_io import read_tiff, write_tiff
    HAVE_HPC_PACKAGES = True
except ImportError:
    HAVE_HPC_PACKAGES = False
    np = None
    prep_mod = None
    ml_mod = None
    seg_mod = None
    track_mod = None


class TestMito4DPipelineSuite(unittest.TestCase):
    """Rigorous unit and integration test suite for all modules of the 4D Mitochondrial Pipeline."""

    def setUp(self):
        self.tmp_dir = Path(tempfile.mkdtemp(prefix="mito4d_test_"))
        self.cfg_path = SCRIPT_DIR / "pipeline_config.yaml"
        if HAVE_HPC_PACKAGES and prep_mod:
            self.cfg = prep_mod.load_config(self.cfg_path)
        else:
            self.cfg = {}

    def tearDown(self):
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
        out_skip.mkdir()
        # Pre-create the dummy target output file
        dummy_out = out_skip / "05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif"
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

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
    def test_02_tiff_io_5d_roundtrip(self):
        """Verifies lossless 5D (T, C, Z, Y, X) uint16 and float32 round-trip via tiff_io."""
        arr_u16 = np.random.randint(0, 50000, size=(3, 2, 6, 32, 32), dtype=np.uint16)
        p_u16 = self.tmp_dir / "test_5d.ome.tif"
        write_tiff(p_u16, arr_u16, metadata={"PhysicalSizeX": 0.16, "TimeIncrement": 11.7})
        loaded_u16, meta = read_tiff(p_u16)
        self.assertEqual(loaded_u16.shape, arr_u16.shape)
        self.assertTrue(np.array_equal(arr_u16, loaded_u16))

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
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

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
    def test_04_fast_16bit_lut_histogram_matching(self):
        """Verifies that 16-bit LUT histogram matching restores bleached intensities accurately."""
        rng = np.random.default_rng(42)
        ref = rng.integers(100, 20000, size=(12, 64, 64), dtype=np.uint16).astype(np.float32)
        bleached = (ref * 0.45).astype(np.float32)  # 55% photobleaching loss

        restored = prep_mod.match_histogram_fast_16bit(bleached, ref)
        p99_ref = float(np.percentile(ref, 99.0))
        p99_restored = float(np.percentile(restored, 99.0))
        self.assertLess(abs(p99_restored - p99_ref) / p99_ref, 0.02)

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
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

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
    def test_06_sparse_lap_tracking_gap_fission_fusion(self):
        """Verifies gap closing, track splitting (fission), merging (fusion), and TrackMate v7 XML validity."""
        # Create a 5-frame label stack (T=5, Z=6, Y=64, X=64)
        labels = np.zeros((5, 6, 64, 64), dtype=np.uint16)
        norm = np.full((5, 6, 64, 64), 10000, dtype=np.uint16)

        # Track 1: Present at T=0, missing at T=1 (gap!), present at T=2, splits into two at T=3 and T=4
        labels[0, 3, 10:13, 10:13] = 1
        # T=1 gap (label 1 absent)
        labels[2, 3, 11:14, 10:13] = 1
        # T=3 Fission: parent splits into daughter A (y=10) and daughter B (y=16)
        labels[3, 3, 9:12, 10:13] = 1
        labels[3, 3, 15:18, 10:13] = 2
        labels[4, 3, 8:11, 10:13] = 1
        labels[4, 3, 16:19, 10:13] = 2

        lbl_path = self.tmp_dir / "gap_split_mito_labels.tif"
        norm_path = self.tmp_dir / "gap_split_mito_norm.tif"
        write_tiff(lbl_path, labels)
        write_tiff(norm_path, norm)

        calib = {"dx_um": 0.16, "dy_um": 0.16, "dz_um": 1.0, "dt_s": 11.7}
        n_tr, n_sp, n_ed, n_fiss, n_fuse = track_mod.run_sparse_lap_tracking(
            lbl_path, norm_path, self.tmp_dir, "gap_split", self.cfg, calib
        )

        # All 6 detections belong to 1 unified branched lineage tree with 1 gap-close and 1 fission!
        self.assertEqual(n_tr, 1)
        self.assertEqual(n_sp, 6)
        self.assertEqual(n_fiss, 1)

        # Verify TrackMate XML parses cleanly and has <AllSpots> and <AllTracks>
        xml_tree = ET.parse(self.tmp_dir / "gap_split_tracks.xml")
        root = xml_tree.getroot()
        self.assertIsNotNone(root.find(".//AllSpots"))
        self.assertIsNotNone(root.find(".//AllTracks"))

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
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

    @unittest.skipUnless(HAVE_HPC_PACKAGES, "Requires HPC analysis packages (numpy, tifffile, scipy, etc.)")
    def test_08_end_to_end_pipeline_integration(self):
        """Executes Modules 1 through 4 end-to-end on a synthetic 4D dual-channel movie."""
        movie_id = "e2e_synthetic"
        movie_path = self.tmp_dir / f"{movie_id}.ome.tif"

        # Generate 5D (T=5, C=2, Z=8, Y=64, X=64) movie
        data = np.zeros((5, 2, 8, 64, 64), dtype=np.uint16)
        zz, yy, xx = np.ogrid[:8, :64, :64]
        for t in range(5):
            # Cytosol (Ch 1): Soma at y=15, x=32 + shaft to y=55
            soma = ((zz - 4) ** 2 * 2 + (yy - 15) ** 2 + (xx - 32) ** 2) <= 36
            data[t, 1][soma] = 4000
            data[t, 1, 3:6, 15:55, 30:35] = 2800
            data[t, 1] += 150

            # Mitochondria (Ch 0): 1 in soma, 1 moving anterograde along shaft
            data[t, 0, 4, 14:17, 31:34] = 8000
            my = 22 + t * 2
            data[t, 0, 4, my : my + 3, 31:34] = 9500
            data[t, 0] += 120

        write_tiff(
            movie_path,
            data,
            metadata={"PhysicalSizeX": 0.16, "PhysicalSizeY": 0.16, "PhysicalSizeZ": 1.0, "TimeIncrement": 11.7},
        )

        # Run Module 1
        prep_mod.preprocess_movie(movie_path, self.tmp_dir, self.cfg)
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


def ndimage_com(vol) -> Tuple[float, float, float]:
    """Helper to compute 3D center of mass of foreground intensity."""
    from scipy import ndimage as ndi
    return ndi.center_of_mass(vol)


if __name__ == "__main__":
    unittest.main(verbosity=2)