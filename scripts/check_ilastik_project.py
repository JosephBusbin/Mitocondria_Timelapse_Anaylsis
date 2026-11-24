#!/usr/bin/env python3
"""
check_ilastik_project.py
Pre-flight check of a trained ilastik Pixel Classification project (.ilp) BEFORE it is used on the
cluster. Catches the failure seen with bestman_mito_2d_v2.ilp (2026-11-24): only 15 Mitochondria and
16 Background voxels were labelled, so the classifier called 99% of the cell "mitochondria" and
Stage 4 tiled the whole cytoplasm into fake instances.

Checks (each PASS / WARN / FAIL):
  - ilastik version matches the HPC binary (1.4.2) and the workflow is Pixel Classification
  - label order: label 1 = "Mitochondria", label 2 = "Background" (foreground_class_index 0)
  - a trained classifier is stored in the project (Live Update was run before saving)
  - raw data stored inside the project ("Copy into project file"), so it is self-contained
  - features computed in 2D (manual Step 3; 1 um Z-steps blur 3D filters)
  - enough labels: >= 150 Mitochondria and >= 300 Background voxels
  - Background labelled INSIDE the cell (needs --cell-mask, the Stage 3 *_cell_mask.tif of the
    training movie, cropped identically): >= 30% of Background voxels must lie inside the cell.
    Without in-cell Background the classifier learns "inside the cell = mitochondria".

Needs h5py: run it with ilastik's own Python on SciClone (no install needed), e.g.
  $MITO_TOOLS_DIR/ilastik-1.4.2-Linux/bin/python scripts/check_ilastik_project.py \
      $MITO_MODELS_DIR/bestman_mito_2d_v3.ilp [--cell-mask <movie>_cell_mask.tif]
(ilastik's bundled Python also has numpy and tifffile for --cell-mask.)
Exit code: 0 if no FAIL, 1 otherwise.
"""

import argparse
import sys

import numpy as np

MIN_MITO_VOXELS = 150
MIN_BG_VOXELS = 300
MIN_BG_INSIDE_FRACTION = 0.30


def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def main() -> int:
    ap = argparse.ArgumentParser(description="Pre-flight check of an ilastik Pixel Classification .ilp")
    ap.add_argument("ilp")
    ap.add_argument("--cell-mask", default=None, help="Stage 3 *_cell_mask.tif (ZYX) of the training movie")
    ap.add_argument("--expect-version", default="1.4.2")
    args = ap.parse_args()
    try:
        import h5py
    except ImportError:
        print("h5py missing: run with ilastik's bundled Python ($MITO_TOOLS_DIR/ilastik-<ver>-Linux/bin/python)")
        return 2

    results = []

    def report(level, msg):
        results.append(level)
        print(f"[{level:4s}] {msg}")

    with h5py.File(args.ilp, "r") as f:
        ver = _s(f["ilastikVersion"][()]) if "ilastikVersion" in f else "?"
        wf = _s(f["workflowName"][()]) if "workflowName" in f else "?"
        report("PASS" if ver == args.expect_version else "FAIL", f"ilastik version {ver} (HPC runs {args.expect_version})")
        report("PASS" if wf == "Pixel Classification" else "FAIL", f"workflow: {wf}")
        pc = f["PixelClassification"]
        names = [_s(x) for x in pc["LabelNames"][()]]
        ok = len(names) >= 2 and names[0].lower().startswith("mito") and names[1].lower().startswith("back")
        report("PASS" if ok else "FAIL", f"label order {names} (need ['Mitochondria', 'Background'])")
        report("PASS" if "ClassifierForests" in pc else "FAIL",
               "trained classifier stored" if "ClassifierForests" in pc else
               "no trained classifier stored: turn Live Update on, wait for predictions, then save")
        if "FeatureSelections/ComputeIn2d" in f:
            in2d = np.asarray(f["FeatureSelections/ComputeIn2d"][()], dtype=bool)
            report("PASS" if in2d.all() else "WARN",
                   "features computed in 2D" if in2d.all() else
                   f"features computed in 3D for {int((~in2d).sum())}/{in2d.size} scales (manual Step 3 asks for 2D)")
        infos = f.get("Input Data/infos")
        shapes = []
        if infos is not None:
            for lane in infos:
                g = infos[lane].get("Raw Data")
                if g is None:
                    continue
                loc = _s(g["location"][()]) if "location" in g else "?"
                report("PASS" if loc == "ProjectInternal" else "WARN",
                       f"lane {lane}: storage {loc}" + ("" if loc == "ProjectInternal" else
                                                         " (use 'Copy into project file')"))
                if "shape" in g:
                    shapes.append(tuple(int(v) for v in g["shape"][()]))

        counts = {}
        bg_inside = bg_total = 0
        mask = None
        if args.cell_mask:
            import tifffile
            mask = tifffile.imread(args.cell_mask) > 0
        for lane in pc["LabelSets"]:
            for blk_name in pc["LabelSets"][lane]:
                blk = pc["LabelSets"][lane][blk_name]
                data = blk[()]
                vals, cnt = np.unique(data, return_counts=True)
                for v, c in zip(vals, cnt):
                    counts[int(v)] = counts.get(int(v), 0) + int(c)
                if mask is not None and "blockSlice" in blk.attrs:
                    # blockSlice like "[0:3,10:20,...]" in the lane's axis order (t,z,y,x[,c])
                    sl = _s(blk.attrs["blockSlice"]).strip("[]").split(",")
                    starts = [int(s.split(":")[0]) for s in sl]
                    nz = np.argwhere(np.squeeze(data, axis=-1) == 2 if data.ndim > len(mask.shape) + 1 else data == 2)
                    for idx in nz:
                        zyx = (idx + np.asarray(starts[: idx.size]))[-3:]
                        bg_total += 1
                        if all(0 <= zyx[k] < mask.shape[k] for k in range(3)) and mask[tuple(zyx)]:
                            bg_inside += 1
        n_mito, n_bg = counts.get(1, 0), counts.get(2, 0)
        report("PASS" if n_mito >= MIN_MITO_VOXELS else "FAIL",
               f"Mitochondria label voxels: {n_mito} (need >= {MIN_MITO_VOXELS}; label dim AND bright ones)")
        report("PASS" if n_bg >= MIN_BG_VOXELS else "FAIL",
               f"Background label voxels: {n_bg} (need >= {MIN_BG_VOXELS}; include in-cell texture and gaps)")
        if mask is not None and bg_total:
            frac = bg_inside / bg_total
            report("PASS" if frac >= MIN_BG_INSIDE_FRACTION else "FAIL",
                   f"Background voxels inside the cell: {100 * frac:.0f}% (need >= {100 * MIN_BG_INSIDE_FRACTION:.0f}%)")
        elif mask is None:
            report("WARN", "in-cell Background not checked (pass --cell-mask <training movie>_cell_mask.tif)")

    n_fail = results.count("FAIL")
    print(f"\nVerdict: {'READY' if n_fail == 0 else 'NOT READY'} ({n_fail} FAIL, {results.count('WARN')} WARN)")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
