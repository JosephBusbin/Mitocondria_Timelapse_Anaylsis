#!/usr/bin/env python3
"""
synthetic_benchmark.py
Ground-truth benchmark for Stages 3-5: generates a realistic synthetic 4D two-channel radial-glia movie
with KNOWN mitochondria (positions, identities, motion, fission/fusion events, compartments), then
scores the pipeline's outputs against that truth.

Why: the practice .sld has one usable series (#39) and no ground truth, so Stage 4B/5 correctness
(touching-object splitting, linking, event calls, geodesic compartments) could only be judged by eye.

What is simulated (all parameters in microns/seconds, converted with the chosen calibration):
  - Cell (channel 2, cytosolic fill): ellipsoidal soma -> curved radial process (tube along a
    quadratic Bezier path) -> endfoot blob; out-of-focus glow; whole-sample drift (random walk in
    XY + occasional Z steps) applied to BOTH channels, so Stage 3 registration is exercised.
  - Mitochondria (channel 1): blurred rods/puncta inside the cell. Classes:
      stationary (soma/endfoot), motile (move along the process, +/- direction, 0.02-0.15 um/s),
      fission pairs (one rod splits into two halves that drift apart; intensity conserved),
      fusion pairs (two halves approach and merge into one rod; intensity conserved),
      touching pairs (two rods end-to-end with a dim neck: must be split by Stage 4B).
  - Dim cytosolic haze in channel 1, exponential photobleaching, camera offset 100 + Gaussian read
    noise + Poisson-like shot noise.

Ground truth (JSON next to the movie): per frame, every mitochondrion's id, centroid in REGISTERED
frame-0 coordinates (um), volume (um^3), compartment, and parent/child ids for events.

Evaluation (after running Stages 3-5 on the movie):
  detection  : per frame, Hungarian matching of predicted spot centroids (spots.csv, crop -> frame-0
               coordinates via provenance crop_origin) to truth within 0.6 um XY / 1.5 um Z;
               precision, recall, F1; touching-pair split rate.
  linking    : truth links (same id in t, t+1) vs predicted FRAME_LINK/GAP_CLOSE edges between
               matched spots -> link precision/recall.
  events     : truth fission/fusion vs predicted SPLIT_FISSION/MERGE_FUSION (+/-1 frame, 2 um).
  motion     : predicted vs true displacement per matched link (median abs error, um).
  compartment: predicted compartment of matched spots vs truth (soma / process / endfoot).

Usage:
  python3 synthetic_benchmark.py generate --out DIR [--preset small|realistic] [--seed 1]
  python3 synthetic_benchmark.py evaluate --truth DIR/<id>_truth.json --work-dir WORK [--json OUT]
Presets: small (T=8, Z=14, 160x200 px; seconds per stage, used by the unit test) and realistic
(T=30, Z=49, 664x952 px, dt 11.7 s - the lab's typical protocol; tests speed/memory at real size).
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tiff_io import write_tiff  # noqa: E402

PRESETS = {
    "small": dict(T=8, Z=14, Y=160, X=200, dx=0.1877, dz=1.0, dt=11.7, n_motile=4, n_stationary=4,
                  n_fission=1, n_fusion=1, n_touching=1),
    "realistic": dict(T=30, Z=49, Y=664, X=952, dx=0.1877, dz=1.0, dt=11.7, n_motile=14, n_stationary=16,
                      n_fission=3, n_fusion=3, n_touching=3),
}


# ----------------------------------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------------------------------
def bezier(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Quadratic Bezier points for parameters s in [0, 1] (rows = points)."""
    s = s[:, None]
    return (1 - s) ** 2 * p0 + 2 * (1 - s) * s * p1 + s ** 2 * p2


class Cell:
    """Radial-glia geometry in voxel coordinates (z, y, x) of the UNdrifted frame-0 field."""

    def __init__(self, shape: Tuple[int, int, int], dx: float, dz: float, rng: np.random.Generator):
        Z, Y, X = shape
        self.dx, self.dz = dx, dz
        zc = Z / 2.0
        self.soma_c = np.array([zc, Y * 0.78, X * 0.72])
        self.soma_r_um = np.array([4.0, 6.0, 7.5])                      # semi-axes (z, y, x) um
        self.p0 = self.soma_c.copy()
        self.p2 = np.array([zc + rng.uniform(-2, 2), Y * 0.16, X * 0.22])  # endfoot centre
        self.p1 = np.array([zc, Y * 0.30 + rng.uniform(-0.05, 0.05) * Y, X * 0.75])  # curvature
        self.proc_r_um = 1.1
        self.endfoot_r_um = np.array([2.5, 4.5, 4.5])
        s = np.linspace(0, 1, 400)
        self.path = bezier(self.p0, self.p1, self.p2, s)                 # voxel coords
        seg = np.diff(self.path * np.array([dz, dx, dx]), axis=0)
        self.arclen_um = np.concatenate([[0.0], np.cumsum(np.linalg.norm(seg, axis=1))])
        self.length_um = float(self.arclen_um[-1])

    def point_at_um(self, d_um: float) -> np.ndarray:
        """Voxel coordinate on the process centreline at arc length d_um from the soma centre."""
        d = float(np.clip(d_um, 0.0, self.length_um))
        return np.array([np.interp(d, self.arclen_um, self.path[:, k]) for k in range(3)])

    def tangent_at_um(self, d_um: float) -> np.ndarray:
        a, b = self.point_at_um(max(0.0, d_um - 0.5)), self.point_at_um(min(self.length_um, d_um + 0.5))
        v = (b - a) * np.array([self.dz, self.dx, self.dx])
        n = np.linalg.norm(v)
        return v / n if n > 0 else np.array([0.0, 1.0, 0.0])

    def soma_exit_um(self) -> float:
        """Arc length where the process leaves the soma ellipsoid (+0.5 um margin)."""
        for d in np.linspace(0, self.length_um, 400):
            p = (self.point_at_um(d) - self.soma_c) * np.array([self.dz, self.dx, self.dx])
            if np.sum((p / self.soma_r_um) ** 2) > 1.0:
                return float(d) + 0.5
        return 0.0

    def mask(self, shape: Tuple[int, int, int]) -> np.ndarray:
        """Binary cell volume (soma + tube + endfoot) on the given grid."""
        Z, Y, X = shape
        zz, yy, xx = np.ogrid[:Z, :Y, :X]
        um = lambda a, c, s: ((a - c) * s)  # noqa: E731
        soma = ((um(zz, self.soma_c[0], self.dz) / self.soma_r_um[0]) ** 2
                + (um(yy, self.soma_c[1], self.dx) / self.soma_r_um[1]) ** 2
                + (um(xx, self.soma_c[2], self.dx) / self.soma_r_um[2]) ** 2) <= 1.0
        endf = ((um(zz, self.p2[0], self.dz) / self.endfoot_r_um[0]) ** 2
                + (um(yy, self.p2[1], self.dx) / self.endfoot_r_um[1]) ** 2
                + (um(xx, self.p2[2], self.dx) / self.endfoot_r_um[2]) ** 2) <= 1.0
        tube = np.zeros(shape, dtype=bool)
        r_vox_xy = int(math.ceil(self.proc_r_um / self.dx)) + 1
        r_vox_z = int(math.ceil(self.proc_r_um / self.dz)) + 1
        for p in self.path[::2]:
            z0, y0, x0 = np.round(p).astype(int)
            zs = slice(max(0, z0 - r_vox_z), min(Z, z0 + r_vox_z + 1))
            ys = slice(max(0, y0 - r_vox_xy), min(Y, y0 + r_vox_xy + 1))
            xs = slice(max(0, x0 - r_vox_xy), min(X, x0 + r_vox_xy + 1))
            lz, ly, lx = np.ogrid[zs, ys, xs]
            d2 = (((lz - p[0]) * self.dz) ** 2 + ((ly - p[1]) * self.dx) ** 2 + ((lx - p[2]) * self.dx) ** 2)
            tube[zs, ys, xs] |= d2 <= self.proc_r_um ** 2
        return soma | tube | endf

    def compartment_of(self, d_um: float, where: str) -> str:
        return where if where in ("soma", "endfoot") else "process"


# ----------------------------------------------------------------------------------------------
# Generation
# ----------------------------------------------------------------------------------------------
def rod_kernel(length_um: float, width_um: float, direction: np.ndarray, dx: float, dz: float) -> np.ndarray:
    """Binary rod (capsule) voxels centred at origin, as an (N, 3) array of voxel offsets."""
    half = length_um / 2.0 + width_um
    rz, rxy = int(math.ceil(half / dz)) + 1, int(math.ceil(half / dx)) + 1
    zz, yy, xx = np.mgrid[-rz:rz + 1, -rxy:rxy + 1, -rxy:rxy + 1]
    pts = np.stack([zz * dz, yy * dx, xx * dx], axis=-1).reshape(-1, 3)
    u = direction / max(np.linalg.norm(direction), 1e-9)
    t = np.clip(pts @ u, -length_um / 2.0, length_um / 2.0)
    d = np.linalg.norm(pts - t[:, None] * u[None, :], axis=1)
    keep = d <= width_um / 2.0
    return np.stack([zz.ravel()[keep], yy.ravel()[keep], xx.ravel()[keep]], axis=1)


def generate(out_dir: Path, preset: str = "small", seed: int = 1, movie_id: str = None) -> Path:
    """Writes <id>.ome.tif (TCZYX uint16) and <id>_truth.json; returns the truth path."""
    P = dict(PRESETS[preset])
    rng = np.random.default_rng(seed)
    T, Z, Y, X, dx, dz, dt = P["T"], P["Z"], P["Y"], P["X"], P["dx"], P["dz"], P["dt"]
    movie_id = movie_id or f"synthetic_{preset}_seed{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    cell = Cell((Z, Y, X), dx, dz, rng)
    cell_mask = cell.mask((Z, Y, X))
    d_exit = cell.soma_exit_um()
    d_endfoot = cell.length_um - 1.5

    # --- mitochondria definitions (arc-length parametrised along the process, or fixed points) ---
    objs: List[Dict] = []
    next_id = [1]

    def new(kind, where, d0, v_um_s=0.0, length=1.2, width=0.45, amp=1.0, t0=0, t1=T - 1, offset=None):
        o = dict(id=next_id[0], kind=kind, where=where, d0=d0, v=v_um_s, length=length, width=width,
                 amp=amp, t0=t0, t1=t1, offset=np.zeros(3) if offset is None else np.asarray(offset, float),
                 parent=None, children=[])
        next_id[0] += 1
        objs.append(o)
        return o

    def soma_point():
        while True:
            p = cell.soma_c + rng.uniform(-1, 1, 3) * cell.soma_r_um * 0.7 / np.array([dz, dx, dx])
            if cell_mask[tuple(np.clip(np.round(p).astype(int), 0, [Z - 1, Y - 1, X - 1]))]:
                return p

    def endfoot_point():
        while True:
            p = cell.p2 + rng.uniform(-1, 1, 3) * cell.endfoot_r_um * 0.6 / np.array([dz, dx, dx])
            if cell_mask[tuple(np.clip(np.round(p).astype(int), 0, [Z - 1, Y - 1, X - 1]))]:
                return p

    proc_lo, proc_hi = d_exit + 1.0, d_endfoot - 3.0
    total_time = dt * max(1, T - 1)
    # Each process object (or pair) gets its own window along the process: window length = the
    # object's extent + 1 um gap, so the TRUTH never overlaps unless an object is deliberately fast.
    # Windows are laid out in random order; objects that do not fit are dropped (with a note).
    plan = (["motile_slow"] * (P["n_motile"] - P["n_motile"] // 2) + ["motile_fast"] * (P["n_motile"] // 2)
            + ["touching"] * P["n_touching"] + ["fission"] * P["n_fission"] + ["fusion"] * P["n_fusion"])
    rng.shuffle(plan)
    widths = {"motile_slow": 3.0, "motile_fast": 3.0, "touching": 3.6, "fission": 4.2, "fusion": 4.2}
    cursor = proc_lo
    placed = []
    for kind in plan:
        w = widths[kind]
        if cursor + w > proc_hi:
            print(f"[GENERATE] note: process too short for another '{kind}' window; dropped")
            continue
        placed.append((kind, cursor + w / 2.0, w))
        cursor += w

    for i in range(P["n_stationary"]):
        where = "soma" if i % 2 == 0 else "endfoot"
        p = soma_point() if where == "soma" else endfoot_point()
        new("stationary", where, None, length=rng.uniform(0.4, 1.4), offset=p)
    touching_pairs = []
    t_ev = max(2, T // 2)
    for kind, d, w in placed:
        if kind == "motile_slow":   # stays inside its window for the whole movie
            length = rng.uniform(0.8, 1.6)
            v = rng.choice([-1, 1]) * rng.uniform(0.3, 0.9) * (w - length - 1.0) / 2.0 / total_time
            new("motile", "process", d - v * total_time / 2.0, v_um_s=v, length=length)
        elif kind == "motile_fast":  # realistic transport speed; may pass neighbours ("crowded" cases)
            v = rng.choice([-1, 1]) * rng.uniform(0.02, 0.15)
            new("motile_fast", "process", d, v_um_s=v, length=rng.uniform(0.8, 1.6))
        elif kind == "touching":
            a = new("touching", "process", d - 0.75, length=1.2, amp=1.0)
            b = new("touching", "process", d + 0.75, length=1.2, amp=0.9)
            touching_pairs.append([a["id"], b["id"]])
    for kind, d, w in placed:
        if kind != "fission":
            continue
        parent = new("fission_parent", "process", d, length=2.0, t1=t_ev - 1)
        for sgn in (-1, 1):
            v_out = min(0.04, 0.9 / max(dt * (T - 1 - t_ev), 1e-6))   # stays within its window
            ch = new("fission_child", "process", d + sgn * 0.6, v_um_s=sgn * v_out, length=1.0, t0=t_ev)
            ch["parent"] = parent["id"]
            parent["children"].append(ch["id"])
    for kind, d, w in placed:
        if kind != "fusion":
            continue
        srcs = []
        for sgn in (-1, 1):
            v_in = min(0.04, 0.9 / max(dt * (t_ev - 1), 1e-6))   # start <= 1.5 um from centre
            s = new("fusion_source", "process", d + sgn * (0.6 + v_in * dt * (t_ev - 1)), v_um_s=-sgn * v_in,
                    length=1.0, t1=t_ev - 1)
            srcs.append(s)
        child = new("fusion_product", "process", d, length=2.0, t0=t_ev)
        for s in srcs:
            s["children"].append(child["id"])
        child["parent"] = [s["id"] for s in srcs]

    # --- drift (applied to both channels): random walk XY ~1.5 px/frame, a Z step mid-movie ---
    drift = np.zeros((T, 3))
    for t in range(1, T):
        drift[t] = drift[t - 1] + np.array([0.0, *rng.normal(0, 1.5, 2)])
        if t == T // 3:
            drift[t, 0] += 1.0

    # --- render ---
    movie = np.zeros((T, 2, Z, Y, X), dtype=np.float32)
    cell_blur = ndimage.gaussian_filter(cell_mask.astype(np.float32), (0.8, 1.5, 1.5))
    glow = ndimage.gaussian_filter(cell_mask.astype(np.float32), (4, 25, 25))
    truth_frames = []
    bleach = 0.985
    for t in range(T):
        mito = np.zeros((Z, Y, X), dtype=np.float32)
        frame_truth = []
        for o in objs:
            if not (o["t0"] <= t <= o["t1"]):
                continue
            if o["d0"] is None:
                c = o["offset"]
                direction = np.array([0.0, 1.0, 1.0])
            else:
                d = o["d0"] + o["v"] * dt * (t - o["t0"] if o["kind"] != "fusion_source" else t)
                d = float(np.clip(d, d_exit, d_endfoot))
                c = cell.point_at_um(d)
                direction = cell.tangent_at_um(d)
            vox = rod_kernel(o["length"], o["width"], direction, dx, dz) + np.round(c).astype(int)
            ok = np.all((vox >= 0) & (vox < np.array([Z, Y, X])), axis=1)
            vox = vox[ok]
            mito[vox[:, 0], vox[:, 1], vox[:, 2]] = np.maximum(mito[vox[:, 0], vox[:, 1], vox[:, 2]], o["amp"])
            comp = o["where"]
            frame_truth.append(dict(id=o["id"], kind=o["kind"], compartment=comp,
                                    z_um=float(c[0] * dz), y_um=float(c[1] * dx), x_um=float(c[2] * dx),
                                    volume_um3=float(len(vox) * dx * dx * dz)))
        truth_frames.append(frame_truth)
        mito = ndimage.gaussian_filter(mito, (0.6, 1.2, 1.2))
        if mito.max() > 0:
            mito /= mito.max()
        b = bleach ** t
        ch1 = (55.0 * mito + 4.0 * cell_blur + 3.0 * glow) * b
        ch2 = (60.0 * cell_blur + 25.0 * glow) * b
        # drift: shift the rendered frame (sample moved); truth stays in frame-0 coordinates
        if t > 0:
            ch1 = ndimage.shift(ch1, drift[t], order=1, mode="constant")
            ch2 = ndimage.shift(ch2, drift[t], order=1, mode="constant")
        for c_idx, ch in enumerate((ch1, ch2)):
            shot = rng.normal(0.0, 1.0, ch.shape).astype(np.float32) * np.sqrt(np.maximum(ch, 0) * 0.5)
            movie[t, c_idx] = 100.0 + ch + shot + rng.normal(0, 2.5, ch.shape).astype(np.float32)
    raw = np.clip(movie, 0, 65535).astype(np.uint16)
    tif = out_dir / f"{movie_id}.ome.tif"
    write_tiff(tif, raw, metadata={"PhysicalSizeX": dx, "PhysicalSizeY": dx, "PhysicalSizeZ": dz,
                                   "TimeIncrement": dt}, axes="TCZYX")

    events = []
    for o in objs:
        if o["kind"] == "fission_parent":
            events.append(dict(type="fission", t_from=o["t1"], parent=o["id"], children=o["children"]))
        if o["kind"] == "fusion_product":
            events.append(dict(type="fusion", t_to=o["t0"], sources=o["parent"], child=o["id"]))
    truth = dict(movie_id=movie_id, preset=preset, seed=seed, shape_TCZYX=list(raw.shape),
                 calibration=dict(dx_um=dx, dy_um=dx, dz_um=dz, dt_s=dt),
                 drift_px=drift.tolist(), process_length_um=cell.length_um, frames=truth_frames,
                 events=events, touching_pair_ids=touching_pairs,
                 movie_path=str(tif))
    tp = out_dir / f"{movie_id}_truth.json"
    tp.write_text(json.dumps(truth, indent=1))
    print(f"[GENERATE] {tif} {raw.shape} ({raw.nbytes / 1e9:.2f} GB) objects={len(objs)} events={len(events)}")
    return tp


# ----------------------------------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------------------------------
def _read_csv(p: Path) -> List[Dict[str, str]]:
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def score_labels(truth: Dict, labels: np.ndarray, off: np.ndarray, tol_xy_um: float = 0.6,
                 tol_z_um: float = 1.5) -> Dict:
    """Stage-4 instance detection vs truth (label centroids, Hungarian matching), with per-kind recall
    and the touching-pair split rate. Used by evaluate() and by parameter sweeps."""
    cal = truth["calibration"]
    T = len(truth["frames"])
    vz, vy, vx = cal["dz_um"], cal["dy_um"], cal["dx_um"]
    s4 = dict(tp=0, fp=0, fn=0)
    kind_hit: Dict[str, List[int]] = {}
    pair_split = pair_total = 0
    pairs = truth.get("touching_pair_ids", [])
    for t in range(min(T, labels.shape[0])):
        tr = truth["frames"][t]
        lab = labels[t]
        n = int(lab.max())
        cents = np.array(ndimage.center_of_mass(lab > 0, lab, range(1, n + 1))) if n else np.zeros((0, 3))
        if n:
            cents = cents * np.array([vz, vy, vx]) + off
        hit: Dict[int, int] = {}
        if tr and n:
            A = np.array([[o["z_um"], o["y_um"], o["x_um"]] for o in tr])
            dxy = np.linalg.norm(A[:, None, 1:] - cents[None, :, 1:], axis=2)
            dzz = np.abs(A[:, None, 0] - cents[None, :, 0])
            cost = np.where((dxy <= tol_xy_um) & (dzz <= tol_z_um), dxy, 1e6)
            r, c = linear_sum_assignment(cost)
            ok = cost[r, c] < 1e6
            hit = {int(i): int(j) for i, j in zip(r[ok], c[ok])}
            s4["tp"] += len(hit)
            s4["fp"] += n - len(hit)
        else:
            s4["fp"] += n
        s4["fn"] += len(tr) - len(hit)
        idx = {o["id"]: i for i, o in enumerate(tr)}
        if tr:
            A = np.array([[o["z_um"], o["y_um"], o["x_um"]] for o in tr])
            dmat = np.linalg.norm(A[:, None, :] - A[None, :, :], axis=2) + np.eye(len(tr)) * 1e9
            nearest = dmat.min(axis=1)
        for i, o in enumerate(tr):
            kind_hit.setdefault(o["kind"], [0, 0])
            kind_hit[o["kind"]][0] += int(i in hit)
            kind_hit[o["kind"]][1] += 1
            grp = "isolated(>=2um)" if nearest[i] >= 2.0 else "crowded(<2um)"
            kind_hit.setdefault(grp, [0, 0])
            kind_hit[grp][0] += int(i in hit)
            kind_hit[grp][1] += 1
        for a_id, b_id in pairs:
            if a_id in idx and b_id in idx:
                pair_total += 1
                pair_split += int(idx[a_id] in hit and idx[b_id] in hit)
    p4 = s4["tp"] / max(1, s4["tp"] + s4["fp"])
    r4 = s4["tp"] / max(1, s4["tp"] + s4["fn"])
    s4.update(precision=round(p4, 3), recall=round(r4, 3), f1=round(2 * p4 * r4 / max(1e-9, p4 + r4), 3),
              touching_pairs_split=f"{pair_split}/{pair_total}",
              recall_by_kind={k: f"{a}/{b}" for k, (a, b) in sorted(kind_hit.items())})
    return s4


def evaluate(truth_path: Path, work_dir: Path, tol_xy_um: float = 0.6, tol_z_um: float = 1.5) -> Dict:
    truth = json.loads(Path(truth_path).read_text())
    mid = truth["movie_id"]
    prov = json.loads((work_dir / f"{mid}_provenance.json").read_text())
    cx, cy, cz = prov["stage3_preprocessing_metrics"]["crop_origin_xyz_px"]
    cal = truth["calibration"]
    off = np.array([cz * cal["dz_um"], cy * cal["dy_um"], cx * cal["dx_um"]])  # crop -> frame-0 um
    ml_path = work_dir / f"{mid}_snt_mitolandscape_tracks.csv"
    spots = _read_csv(ml_path if ml_path.exists() else work_dir / f"{mid}_spots.csv")
    edges = _read_csv(work_dir / f"{mid}_edges.csv")
    T = len(truth["frames"])
    by_frame_pred: Dict[int, List[Dict]] = {t: [] for t in range(T)}
    for s in spots:
        by_frame_pred[int(s["FRAME"])].append(s)

    match: Dict[Tuple[int, int], str] = {}   # (t, truth_id) -> predicted spot ID
    tp = fp = fn = 0
    comp_ok = comp_n = 0
    for t in range(T):
        tr = truth["frames"][t]
        pr = by_frame_pred[t]
        if not tr or not pr:
            fn += len(tr)
            fp += len(pr)
            continue
        A = np.array([[o["z_um"], o["y_um"], o["x_um"]] for o in tr])
        B = np.array([[float(s["POSITION_Z"]), float(s["POSITION_Y"]), float(s["POSITION_X"])] for s in pr]) + off
        dxy = np.linalg.norm(A[:, None, 1:] - B[None, :, 1:], axis=2)
        dzz = np.abs(A[:, None, 0] - B[None, :, 0])
        cost = np.where((dxy <= tol_xy_um) & (dzz <= tol_z_um), dxy, 1e6)
        r, c = linear_sum_assignment(cost)
        ok = cost[r, c] < 1e6
        tp += int(ok.sum())
        fn += len(tr) - int(ok.sum())
        fp += len(pr) - int(ok.sum())
        for i, j in zip(r[ok], c[ok]):
            match[(t, tr[i]["id"])] = pr[j]["ID"]
            name = pr[j].get("compartment_name")
            if name:
                comp_n += 1
                pred = "process" if name in ("apical_process", "basal_process") else name
                comp_ok += int(pred == tr[i]["compartment"])

    # linking
    pred_links = {(e["SPOT_SOURCE_ID"], e["SPOT_TARGET_ID"]) for e in edges
                  if e.get("EVENT_TYPE") in ("FRAME_LINK", "GAP_CLOSE")}
    l_tp = l_n = 0
    disp_err = []
    pos_pred = {s["ID"]: np.array([float(s["POSITION_Z"]), float(s["POSITION_Y"]), float(s["POSITION_X"])])
                for s in spots}
    for t in range(T - 1):
        ids_t = {o["id"]: o for o in truth["frames"][t]}
        for o in truth["frames"][t + 1]:
            if o["id"] in ids_t and (t, o["id"]) in match and (t + 1, o["id"]) in match:
                l_n += 1
                a, b = match[(t, o["id"])], match[(t + 1, o["id"])]
                if (a, b) in pred_links:
                    l_tp += 1
                    true_d = np.linalg.norm(np.array([o["z_um"], o["y_um"], o["x_um"]])
                                            - np.array([ids_t[o["id"]][k] for k in ("z_um", "y_um", "x_um")]))
                    disp_err.append(abs(float(np.linalg.norm(pos_pred[b] - pos_pred[a])) - true_d))
    matched_ids = set(match.values())
    pred_links_on_matched = {(a, b) for a, b in pred_links if a in matched_ids and b in matched_ids}

    # events
    def near(ev_type, t_ref, ids):
        pts = [np.array([o["z_um"], o["y_um"], o["x_um"]]) for tt in (t_ref, t_ref + 1) if 0 <= tt < T
               for o in truth["frames"][tt] if o["id"] in ids]
        if not pts:
            return False
        c = np.mean(pts, axis=0)
        for e in edges:
            if e.get("EVENT_TYPE") != ev_type:
                continue
            for sid in (e["SPOT_SOURCE_ID"], e["SPOT_TARGET_ID"]):
                if sid in pos_pred and np.linalg.norm(pos_pred[sid] + off - c) <= 2.0:
                    return True
        return False

    fis = [ev for ev in truth["events"] if ev["type"] == "fission"]
    fus = [ev for ev in truth["events"] if ev["type"] == "fusion"]
    fis_found = sum(near("SPLIT_FISSION", ev["t_from"], [ev["parent"], *ev["children"]]) for ev in fis)
    fus_found = sum(near("MERGE_FUSION", ev["t_to"] - 1, [*ev["sources"], ev["child"]]) for ev in fus)
    n_pred_fis = sum(1 for e in edges if e.get("EVENT_TYPE") == "SPLIT_FISSION")
    n_pred_fus = sum(1 for e in edges if e.get("EVENT_TYPE") == "MERGE_FUSION")

    # touching pairs: both members detected as SEPARATE spots
    tp_split = tp_total = 0
    for a, b in truth.get("touching_pair_ids", []):
        for t in range(T):
            ids = {o["id"] for o in truth["frames"][t]}
            if a in ids and b in ids:
                tp_total += 1
                tp_split += int((t, a) in match and (t, b) in match and match[(t, a)] != match[(t, b)])

    prec = tp / max(1, tp + fp)
    rec = tp / max(1, tp + fn)

    # Stage 4 detection from the label stack itself (independent of Stage 5's track-length filter)
    from tiff_io import read_tiff
    labels, _ = read_tiff(work_dir / f"{mid}_mito_labels.tif", axes="TZYX")
    s4 = score_labels(truth, labels, off, tol_xy_um, tol_z_um)

    res = dict(
        stage4_detection=s4,
        movie_id=mid,
        tracked_detection=dict(tp=tp, fp=fp, fn=fn, precision=round(prec, 3), recall=round(rec, 3),
                       f1=round(2 * prec * rec / max(1e-9, prec + rec), 3)),
        touching_pairs_split=f"{tp_split}/{tp_total}",
        linking=dict(true_links_evaluable=l_n, recovered=l_tp, recall=round(l_tp / max(1, l_n), 3),
                     predicted_links_between_matched=len(pred_links_on_matched),
                     precision=round(l_tp / max(1, len(pred_links_on_matched)), 3),
                     median_abs_displacement_error_um=round(float(np.median(disp_err)), 3) if disp_err else None),
        events=dict(fission_true=len(fis), fission_found=fis_found, fission_predicted_total=n_pred_fis,
                    fusion_true=len(fus), fusion_found=fus_found, fusion_predicted_total=n_pred_fus),
        compartment_accuracy=round(comp_ok / comp_n, 3) if comp_n else None,
    )
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description="Synthetic ground-truth benchmark for Stages 3-5")
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--out", type=Path, required=True)
    g.add_argument("--preset", choices=sorted(PRESETS), default="small")
    g.add_argument("--seed", type=int, default=1)
    g.add_argument("--movie-id", default=None)
    e = sub.add_parser("evaluate")
    e.add_argument("--truth", type=Path, required=True)
    e.add_argument("--work-dir", type=Path, required=True)
    e.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    if args.cmd == "generate":
        generate(args.out, args.preset, args.seed, args.movie_id)
        return 0
    res = evaluate(args.truth, args.work_dir)
    print(json.dumps(res, indent=1))
    if args.json:
        args.json.write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
