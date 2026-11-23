#!/usr/bin/env python3
"""
tiff_io.py
Universal, Atomic Multi-Dimensional TIFF & OME-TIFF Reader/Writer

Designed for the Bestman Lab 4D Mitochondrial Pipeline (Windows Lab PC + W&M SciClone HPC).

Key Guarantees:
  1. Consistent API across all pipeline modules:
       arr, meta = read_tiff(path)                      # native axis order (meta["axes"])
       arr, meta = read_tiff(path, axes="TCZYX")        # normalized axis order
       write_tiff(path, arr, metadata=None, axes=None, compression=None)
  2. Dual-backend transparency:
       - Primary: 'tifffile' (reads native OME-XML, ImageJ hyperstack metadata, BigTIFF,
         and LZW/Deflate compressed stacks via 'imagecodecs').
       - Fallback: 'PIL.Image' + 'PIL.ImageSequence' with a '.shape.json' sidecar for
         minimal environments lacking 'tifffile' (calibration is NOT embedded in the TIFF).
  3. Calibration preservation (always returned in micrometers / seconds):
       - Extracts PhysicalSizeX/Y/Z and TimeIncrement from OME-XML <Pixels>, converting
         PhysicalSize*Unit / TimeIncrementUnit (nm, mm, ms, min, ...) to um and s, with
         automatic fallback to <Plane TheT="..." DeltaT="..."> intervals (common in 3i
         SlideBook exports) and ImageJ TIFF tags (XResolution/YResolution/spacing/finterval).
       - Calibration that could not be found is listed in meta["_warnings"] so callers can
         record that config defaults were used instead of silently assuming them.
  4. Output format chosen from the FINAL filename:
       - '*.ome.tif' / '*.ome.tiff' -> OME-TIFF (any dtype, BigTIFF when large); opens as a
         calibrated hyperstack in FIJI (Bio-Formats) and ilastik.
       - other '*.tif' -> ImageJ hyperstack (axes reordered to ImageJ's TZCYX) for
         uint8/uint16/float32 under 3.5 GB; otherwise a tifffile 'shaped' TIFF (e.g. int32
         label masks) with axes + calibration stored as JSON in ImageDescription.
       - Grayscale is always forced (photometric='minisblack') so a Z or C size of 3/4 is
         never misinterpreted as RGB(A) samples.
  5. Atomic writes:
       - Writes to a hidden temporary file in the target directory and replaces the target
         atomically via os.replace(), preventing corrupt partial TIFFs on SLURM timeouts.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import tifffile
    HAS_TIFFFILE = True
except ImportError:
    HAS_TIFFFILE = False

try:
    import PIL.Image
    import PIL.ImageSequence
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

# ImageJ hyperstacks strictly support uint8, uint16, and float32, and cannot use BigTIFF (>4 GB).
_IMAGEJ_DTYPES = {np.dtype("uint8"), np.dtype("uint16"), np.dtype("float32")}
_IMAGEJ_AXIS_ORDER = "TZCYX"
_BIGTIFF_THRESHOLD_BYTES = int(3.5 * (1024 ** 3))

# Conversion factors to the pipeline's canonical units (micrometers, seconds).
_LENGTH_TO_UM = {
    "um": 1.0, "µm": 1.0, "μm": 1.0, "micron": 1.0, "microns": 1.0,
    "nm": 1e-3, "Å": 1e-4, "A": 1e-4, "mm": 1e3, "cm": 1e4, "m": 1e6, "pm": 1e-6,
}
_TIME_TO_S = {"s": 1.0, "sec": 1.0, "ms": 1e-3, "µs": 1e-6, "us": 1e-6, "min": 60.0, "h": 3600.0}

# Calibration keys callers rely on; missing ones are reported in meta["_warnings"].
_CALIBRATION_KEYS = ("PhysicalSizeX", "PhysicalSizeY", "PhysicalSizeZ", "TimeIncrement")

# OME <Pixels> attributes tifffile's OME writer understands; anything else is dropped for OME output.
_OME_WRITE_KEYS = (
    "PhysicalSizeX", "PhysicalSizeXUnit", "PhysicalSizeY", "PhysicalSizeYUnit",
    "PhysicalSizeZ", "PhysicalSizeZUnit", "TimeIncrement", "TimeIncrementUnit", "Name",
)


def _infer_axes(ndim: int) -> str:
    """Infers canonical microscopy axis labels matching the pipeline's (T, C, Z, Y, X) convention."""
    mapping = {
        2: "YX",
        3: "ZYX",
        4: "TZYX",
        5: "TCZYX",
    }
    return mapping.get(ndim, "Q" * max(0, ndim - 2) + "YX")


def _is_ome_path(path: Path) -> bool:
    """True when the target filename asks for OME-TIFF (e.g. 'movie.ome.tif')."""
    name = path.name.lower()
    return name.endswith(".ome.tif") or name.endswith(".ome.tiff")


def clean_ome_xml(ome_xml: str) -> str:
    """
    Strips namespace declarations and prefixed attributes (e.g. xsi:schemaLocation) from OME-XML.
    Removing only 'xmlns' leaves 'xsi:' attributes with an unbound prefix, which ElementTree
    rejects (bfconvert output always carries xsi:schemaLocation). Mirrors the cleaning used in
    collect_4d_sld_timelapses.parse_ome_xml_metadata.
    """
    xml_clean = re.sub(r"<OME[^>]*>", "<OME>", ome_xml, count=1)
    return re.sub(r'\s\w+:\w+="[^"]*"', "", xml_clean)


def _to_canonical(value: Any, unit: Optional[str], table: Dict[str, float]) -> Optional[float]:
    """Converts a numeric attribute to um or s. Returns None for unknown units or bad values."""
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if not unit:
        return val
    factor = table.get(unit.strip())
    return val * factor if factor is not None else None


def reorder_axes(arr: np.ndarray, src_axes: str, dst_axes: str) -> np.ndarray:
    """
    Reorders an array from src_axes to dst_axes (e.g. 'TZCYX' -> 'TCZYX').
    Axes in dst but not in src are inserted as singletons; axes in src but not in dst are
    dropped only if they have size 1. Raises ValueError for ambiguous non-singleton axes.
    """
    src_axes = src_axes.upper()
    dst_axes = dst_axes.upper()
    if len(src_axes) != arr.ndim:
        raise ValueError(f"Axes '{src_axes}' do not match array ndim {arr.ndim} (shape {arr.shape}).")
    if len(set(dst_axes)) != len(dst_axes):
        raise ValueError(f"Target axes '{dst_axes}' contain duplicates.")

    out = arr
    axes = src_axes
    for ax in src_axes:
        if ax not in dst_axes:
            idx = axes.index(ax)
            if out.shape[idx] != 1:
                raise ValueError(
                    f"Cannot map axis '{ax}' (size {out.shape[idx]}) from '{src_axes}' onto '{dst_axes}'. "
                    "The file's dimension labels are ambiguous; re-export it with explicit T/C/Z axes."
                )
            out = np.squeeze(out, axis=idx)
            axes = axes[:idx] + axes[idx + 1:]
    for ax in dst_axes:
        if ax not in axes:
            out = out[np.newaxis, ...]
            axes = ax + axes
    return np.transpose(out, [axes.index(ax) for ax in dst_axes])


def _parse_ome_xml(ome_xml: str) -> Dict[str, Any]:
    """
    Extracts pixel calibration and time interval metadata from an OME-XML string.
    Physical sizes are returned in um and TimeIncrement in s (unit attributes are consumed).
    """
    meta: Dict[str, Any] = {}
    if not ome_xml:
        return meta

    try:
        root = ET.fromstring(clean_ome_xml(ome_xml))
    except ET.ParseError as exc:
        meta["_warnings"] = [f"OME-XML present but unparseable ({exc}); calibration not read from file."]
        return meta

    pixels = root.find(".//Pixels")
    if pixels is None:
        return meta

    attrs = dict(pixels.attrib)
    for key, val in attrs.items():
        if not key.endswith("Unit"):
            meta[key] = val

    for key in ("PhysicalSizeX", "PhysicalSizeY", "PhysicalSizeZ"):
        if key in attrs:
            conv = _to_canonical(attrs[key], attrs.get(f"{key}Unit"), _LENGTH_TO_UM)
            if conv is None:
                meta.pop(key, None)
                meta.setdefault("_warnings", []).append(
                    f"{key}={attrs[key]} with unsupported unit '{attrs.get(f'{key}Unit')}' ignored."
                )
            else:
                meta[key] = conv

    if "TimeIncrement" in attrs:
        conv = _to_canonical(attrs["TimeIncrement"], attrs.get("TimeIncrementUnit"), _TIME_TO_S)
        if conv is None:
            meta.pop("TimeIncrement", None)
        else:
            meta["TimeIncrement"] = conv

    # SlideBook .sld -> OME-TIFF exports frequently omit TimeIncrement on <Pixels>
    # and instead store per-plane timestamps in <Plane TheT="..." DeltaT="...">.
    if "TimeIncrement" not in meta:
        t_deltas: Dict[int, float] = {}
        for plane in pixels.findall("Plane"):
            t_str = plane.attrib.get("TheT")
            dt_val = _to_canonical(plane.attrib.get("DeltaT"), plane.attrib.get("DeltaTUnit"), _TIME_TO_S)
            if t_str is None or dt_val is None:
                continue
            try:
                t_idx = int(t_str)
            except ValueError:
                continue
            t_deltas.setdefault(t_idx, dt_val)
        if len(t_deltas) >= 2:
            ordered = [t_deltas[k] for k in sorted(t_deltas.keys())]
            diffs = [ordered[i + 1] - ordered[i] for i in range(len(ordered) - 1)]
            meta["TimeIncrement"] = round(sum(diffs) / len(diffs), 4)

    return meta


def _parse_imagej_metadata(ij_meta: Optional[Dict[str, Any]], tif: Any) -> Dict[str, Any]:
    """
    Extracts physical calibration tags from an ImageJ hyperstack header if OME-XML is absent.
    ImageJ only records a length 'unit' when the image is calibrated; without it, XResolution
    (tifffile default 1/1) and 'spacing' are pixel units and must NOT be reported as micrometers.
    """
    meta: Dict[str, Any] = {}
    if not ij_meta:
        return meta
    for key, val in ij_meta.items():
        meta[key] = val

    unit = str(ij_meta.get("unit", "")).strip()
    length_factor = _LENGTH_TO_UM.get(unit)
    if "finterval" in ij_meta and "TimeIncrement" not in meta:
        # ImageJ frame intervals are seconds unless 'tunit' says otherwise.
        dt = _to_canonical(ij_meta["finterval"], ij_meta.get("tunit") or "s", _TIME_TO_S)
        if dt is not None and dt > 0:
            meta["TimeIncrement"] = dt
    if length_factor is None:
        return meta

    if "spacing" in ij_meta and "PhysicalSizeZ" not in meta:
        dz = _to_canonical(ij_meta["spacing"], unit, _LENGTH_TO_UM)
        if dz is not None and dz > 0:
            meta["PhysicalSizeZ"] = dz

    # XY pixel size from TIFF resolution tags (pixels per 'unit' in ImageJ files).
    page0 = tif.pages[0]
    for tag_name, key in (("XResolution", "PhysicalSizeX"), ("YResolution", "PhysicalSizeY")):
        tag = page0.tags.get(tag_name)
        if tag is None:
            continue
        try:
            num, den = tag.value
            if num and den:
                meta.setdefault(key, round(float(den) / float(num) * length_factor, 6))
        except (TypeError, ValueError, ZeroDivisionError):
            continue

    return meta


def _build_imagej_write_kwargs(
    axes: str,
    metadata: Optional[Dict[str, Any]],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Prepares ImageJ-compatible metadata and resolution tags without losing custom keys."""
    ij_meta: Dict[str, Any] = {"axes": axes}
    extra_kwargs: Dict[str, Any] = {}
    if not metadata:
        return ij_meta, extra_kwargs

    for k, v in metadata.items():
        if k != "axes" and isinstance(v, (str, int, float, bool)):
            ij_meta[k] = v

    dx, dy, dz, dt = _calibration_from(metadata)
    if dx is not None and dy is not None:
        extra_kwargs["resolution"] = (1.0 / dx, 1.0 / dy)
        ij_meta["unit"] = "um"
    if dz is not None:
        ij_meta["spacing"] = dz
        ij_meta.setdefault("unit", "um")
    if dt is not None:
        ij_meta["finterval"] = dt

    return ij_meta, extra_kwargs


def _calibration_from(metadata: Dict[str, Any]) -> Tuple[Optional[float], ...]:
    """
    Returns (dx_um, dy_um, dz_um, dt_s) from either pipeline keys (dx_um/dy_um/dz_um/dt_s) or
    OME keys (PhysicalSizeX/Y/Z, TimeIncrement). Non-positive or non-numeric values become None.
    """
    def _pos(*keys: str) -> Optional[float]:
        for key in keys:
            if metadata.get(key) is None:
                continue
            try:
                val = float(metadata[key])
            except (TypeError, ValueError):
                continue
            if val > 0:
                return val
        return None

    dx = _pos("dx_um", "PhysicalSizeX")
    dy = _pos("dy_um", "PhysicalSizeY") or dx
    return dx, dy, _pos("dz_um", "PhysicalSizeZ"), _pos("dt_s", "TimeIncrement")


def _build_ome_metadata(axes: str, metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Builds a tifffile OME metadata dict (only keys its OME writer understands)."""
    ome_meta: Dict[str, Any] = {"axes": axes}
    for key in _OME_WRITE_KEYS:
        if key in metadata and metadata[key] is not None:
            ome_meta[key] = metadata[key]
    dx, dy, dz, dt = _calibration_from(metadata)
    for key, val, unit_key, unit in (
        ("PhysicalSizeX", dx, "PhysicalSizeXUnit", "µm"),
        ("PhysicalSizeY", dy, "PhysicalSizeYUnit", "µm"),
        ("PhysicalSizeZ", dz, "PhysicalSizeZUnit", "µm"),
        ("TimeIncrement", dt, "TimeIncrementUnit", "s"),
    ):
        if val is not None:
            ome_meta[key] = val
            ome_meta[unit_key] = unit
    return ome_meta


def _json_safe(metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Keeps only JSON-serializable values for the 'shaped' ImageDescription."""
    safe: Dict[str, Any] = {}
    for k, v in metadata.items():
        try:
            json.dumps(v)
        except (TypeError, ValueError):
            continue
        safe[k] = v
    return safe


def read_tiff(
    path: Union[str, Path],
    axes: Optional[str] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Reads a 2D-5D TIFF or OME-TIFF file and returns (numpy_array, metadata_dict).

    Parameters
    ----------
    path : str or Path
        Path to the input .tif / .tiff / .ome.tif file.
    axes : str, optional
        Target axis order (e.g. 'TCZYX'). When given, the array is reordered using the file's
        own dimension labels, missing axes are inserted as size-1, and meta['axes'] is updated.
        When omitted, the array is returned in the file's native order (meta['axes']).

    Returns
    -------
    arr : np.ndarray
        Image data preserving original dtype.
    meta : dict
        Combined metadata dictionary (OME-XML <Pixels>, ImageJ tags, axis labels,
        and any '.shape.json' sidecar fields). Calibration is in um / s. 'meta["_warnings"]'
        lists any calibration keys that were absent from the file.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"TIFF file not found: {path}")

    meta: Dict[str, Any] = {}
    warn_list: List[str] = []
    arr: np.ndarray

    if HAS_TIFFFILE:
        with tifffile.TiffFile(path) as tif:
            arr = tif.asarray()
            if tif.series:
                meta["axes"] = tif.series[0].axes
            if tif.ome_metadata:
                ome_meta = _parse_ome_xml(tif.ome_metadata)
                warn_list.extend(ome_meta.pop("_warnings", []))
                meta.update(ome_meta)
            if tif.imagej_metadata:
                ij_parsed = _parse_imagej_metadata(tif.imagej_metadata, tif)
                for k, v in ij_parsed.items():
                    meta.setdefault(k, v)
            # JSON-encoded metadata placed in ImageDescription by tifffile 'shaped' writes.
            desc = tif.pages[0].description
            if desc and desc.lstrip().startswith("{"):
                try:
                    parsed_desc = json.loads(desc)
                except ValueError:
                    parsed_desc = None
                if isinstance(parsed_desc, dict):
                    for k, v in parsed_desc.items():
                        if k != "shape":
                            meta.setdefault(k, v)
                    dx, dy, dz, dt = _calibration_from(parsed_desc)
                    for key, val in (("PhysicalSizeX", dx), ("PhysicalSizeY", dy),
                                     ("PhysicalSizeZ", dz), ("TimeIncrement", dt)):
                        if val is not None:
                            meta.setdefault(key, val)
    elif HAS_PIL:
        with PIL.Image.open(path) as img:
            pages = [np.array(page) for page in PIL.ImageSequence.Iterator(img)]
            arr = np.stack(pages, axis=0) if len(pages) > 1 else pages[0]
    else:
        raise ImportError("Neither 'tifffile' nor 'Pillow' is installed; cannot read TIFF files.")

    # Restore exact shape, dtype, and metadata if a PIL fallback sidecar exists.
    shape_file = path.with_suffix(path.suffix + ".shape.json")
    if shape_file.is_file():
        try:
            with open(shape_file, "r", encoding="utf-8") as f:
                saved = json.load(f)
        except (OSError, ValueError) as exc:
            warn_list.append(f"Ignoring unreadable sidecar {shape_file.name}: {exc}")
            saved = {}
        orig_shape = saved.get("shape")
        if orig_shape and int(np.prod(orig_shape)) == arr.size:
            arr = arr.reshape(orig_shape)
        orig_dtype = saved.get("dtype")
        if orig_dtype and str(arr.dtype) != orig_dtype:
            arr = arr.astype(np.dtype(orig_dtype))
        if saved.get("axes"):
            meta["axes"] = saved["axes"]
        for k, v in saved.get("metadata", {}).items():
            meta.setdefault(k, v)

    meta.setdefault("axes", _infer_axes(arr.ndim))
    if len(meta["axes"]) != arr.ndim:
        meta["axes"] = _infer_axes(arr.ndim)

    if axes is not None:
        arr = reorder_axes(arr, meta["axes"], axes)
        meta["axes"] = axes.upper()

    for key in _CALIBRATION_KEYS:
        if key not in meta:
            warn_list.append(f"{key} missing from {path.name}; caller must fall back to config defaults.")
    if warn_list:
        meta["_warnings"] = warn_list
    return arr, meta


def write_tiff(
    path: Union[str, Path],
    arr: np.ndarray,
    metadata: Optional[Dict[str, Any]] = None,
    axes: Optional[str] = None,
    compression: Optional[str] = None,
) -> None:
    """
    Atomically writes a 2D-5D NumPy array to a TIFF file.

    The format is chosen from the final filename: '*.ome.tif' -> OME-TIFF; otherwise an
    ImageJ hyperstack (uint8/uint16/float32, <3.5 GB, axes reordered to TZCYX) or a 'shaped'
    TIFF for other dtypes (e.g. int32 label masks) or very large stacks.

    Parameters
    ----------
    path : str or Path
        Destination .tif / .ome.tif path (parent directories created automatically).
    arr : np.ndarray
        2D-5D image or label array.
    metadata : dict, optional
        Calibration and provenance dictionary (e.g., dx_um, dz_um, dt_s, PhysicalSizeX).
    axes : str, optional
        Explicit axis order string of `arr` (e.g., 'TCZYX', 'TZYX', 'ZYX'). Inferred if omitted.
    compression : str, optional
        Optional tifffile compression codec (e.g., 'zlib', 'lzw'). Uncompressed by default
        for maximum speed and native ImageJ compatibility.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.ascontiguousarray(arr)
    axes_str = (axes or _infer_axes(arr.ndim)).upper()
    if len(axes_str) != arr.ndim:
        raise ValueError(f"Axes '{axes_str}' do not match array ndim {arr.ndim} (shape {arr.shape}).")
    meta_out = {k: v for k, v in (metadata or {}).items() if k != "_warnings"}

    # Write to a temporary file in the same directory so os.replace() is atomic on POSIX/Windows.
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.stem}_", suffix=".tmp.tif", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)

    try:
        if HAS_TIFFFILE:
            use_bigtiff = arr.nbytes >= _BIGTIFF_THRESHOLD_BYTES
            # photometric='minisblack' stops tifffile treating a size-3/4 Z or C axis as RGB(A).
            common = {"compression": compression, "photometric": "minisblack"}

            if _is_ome_path(path):
                # The temp name is not '*.ome.tif', so OME output must be requested explicitly.
                tifffile.imwrite(
                    tmp_path, arr, ome=True, bigtiff=use_bigtiff,
                    metadata=_build_ome_metadata(axes_str, meta_out), **common,
                )
            elif (
                not use_bigtiff
                and arr.dtype in _IMAGEJ_DTYPES
                and set(axes_str) <= set(_IMAGEJ_AXIS_ORDER)
            ):
                ij_axes = "".join(ax for ax in _IMAGEJ_AXIS_ORDER if ax in axes_str)
                ij_arr = reorder_axes(arr, axes_str, ij_axes)
                ij_meta, extra_kwargs = _build_imagej_write_kwargs(ij_axes, meta_out)
                tifffile.imwrite(
                    tmp_path, ij_arr, imagej=True, metadata=ij_meta, **extra_kwargs, **common,
                )
            else:
                tifffile.imwrite(
                    tmp_path, arr, ome=False, bigtiff=use_bigtiff,
                    metadata={**_json_safe(meta_out), "axes": axes_str}, **common,
                )
        elif HAS_PIL:
            warnings.warn(
                f"tifffile unavailable: writing {path.name} with Pillow; calibration and axes are "
                "stored only in the '.shape.json' sidecar, not inside the TIFF.",
                RuntimeWarning,
                stacklevel=2,
            )
            flat = arr.reshape((-1, arr.shape[-2], arr.shape[-1])) if arr.ndim > 2 else arr[np.newaxis, ...]
            pages = []
            try:
                for i in range(flat.shape[0]):
                    sl = flat[i]
                    if arr.dtype == np.uint16:
                        pages.append(PIL.Image.fromarray(sl, mode="I;16"))
                    elif arr.dtype == np.uint8:
                        pages.append(PIL.Image.fromarray(sl, mode="L"))
                    elif arr.dtype == np.float32:
                        pages.append(PIL.Image.fromarray(sl, mode="F"))
                    else:
                        pages.append(PIL.Image.fromarray(sl.astype(np.float32), mode="F"))
                pages[0].save(tmp_path, save_all=True, append_images=pages[1:], compression=None)
            finally:
                for im in pages:
                    im.close()

            # PIL does not store multi-dimensional shapes or custom calibration dicts natively.
            shape_file = path.with_suffix(path.suffix + ".shape.json")
            with open(shape_file, "w", encoding="utf-8") as f:
                json.dump(
                    {"shape": list(arr.shape), "dtype": str(arr.dtype), "axes": axes_str,
                     "metadata": _json_safe(meta_out)},
                    f,
                    indent=2,
                )
        else:
            raise ImportError("Neither 'tifffile' nor 'Pillow' is installed; cannot write TIFF files.")

        os.replace(tmp_path, path)
    except BaseException:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise
