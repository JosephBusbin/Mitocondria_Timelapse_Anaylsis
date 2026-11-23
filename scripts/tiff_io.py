#!/usr/bin/env python3
"""
tiff_io.py
Universal, Atomic Multi-Dimensional TIFF & OME-TIFF Reader/Writer

Designed for the Bestman Lab 4D Mitochondrial Pipeline (Windows Lab PC + W&M SciClone HPC).

Key Guarantees:
  1. Consistent API across all pipeline modules:
       arr, meta = read_tiff(path)
       write_tiff(path, arr, metadata=None, axes=None, compression=None)
  2. Dual-backend transparency:
       - Primary: 'tifffile' (reads native OME-XML, ImageJ hyperstack metadata, BigTIFF,
         and LZW/Deflate compressed stacks via 'imagecodecs').
       - Fallback: 'PIL.Image' + 'PIL.ImageSequence' with a '.shape.json' sidecar for
         minimal environments lacking 'tifffile'.
  3. Calibration preservation:
       - Extracts PhysicalSizeX/Y/Z and TimeIncrement from OME-XML <Pixels>, with automatic
         fallback to <Plane TheT="..." DeltaT="..."> intervals (common in 3i SlideBook
         exports) and ImageJ TIFF tags (XResolution/YResolution/spacing/finterval).
  4. Safe ImageJ vs. BigTIFF routing:
       - Writes FIJI/ilastik-friendly ImageJ hyperstacks for uint8/uint16/float32 arrays
         under 3.5 GB, and automatically switches to standard/BigTIFF for int32/uint32
         label masks or >3.5 GB volumes so writes never crash.
  5. Atomic writes:
       - Writes to a hidden temporary file in the target directory and replaces the target
         atomically via os.replace(), preventing corrupt partial TIFFs on SLURM timeouts.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

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
_BIGTIFF_THRESHOLD_BYTES = int(3.5 * (1024 ** 3))


def _infer_axes(ndim: int) -> str:
    """Infers canonical microscopy axis labels matching the pipeline's (T, C, Z, Y, X) convention."""
    mapping = {
        2: "YX",
        3: "ZYX",
        4: "TZYX",
        5: "TCZYX",
    }
    return mapping.get(ndim, "Q" * max(0, ndim - 2) + "YX")


def _parse_ome_xml(ome_xml: str) -> Dict[str, Any]:
    """Extracts pixel calibration and time interval metadata from an OME-XML string."""
    meta: Dict[str, Any] = {}
    if not ome_xml:
        return meta

    xml_clean = re.sub(r'\sxmlns(:\w+)?="[^"]+"', "", ome_xml)
    try:
        root = ET.fromstring(xml_clean)
    except ET.ParseError:
        return meta

    pixels = root.find(".//Pixels")
    if pixels is None:
        return meta

    for key, val in pixels.attrib.items():
        meta[key] = val

    # SlideBook .sld -> OME-TIFF exports frequently omit TimeIncrement on <Pixels>
    # and instead store per-plane timestamps in <Plane TheT="..." DeltaT="...">.
    if "TimeIncrement" not in meta:
        t_deltas: Dict[int, float] = {}
        for plane in pixels.findall("Plane"):
            t_str = plane.attrib.get("TheT")
            dt_str = plane.attrib.get("DeltaT")
            if t_str is not None and dt_str is not None:
                try:
                    t_idx = int(t_str)
                    if t_idx not in t_deltas:
                        t_deltas[t_idx] = float(dt_str)
                except ValueError:
                    continue
        if len(t_deltas) >= 2:
            ordered = [t_deltas[k] for k in sorted(t_deltas.keys())]
            diffs = [ordered[i + 1] - ordered[i] for i in range(len(ordered) - 1)]
            if diffs:
                meta["TimeIncrement"] = round(sum(diffs) / len(diffs), 4)

    return meta


def _parse_imagej_metadata(ij_meta: Optional[Dict[str, Any]], tif: Any) -> Dict[str, Any]:
    """Extracts physical calibration tags from an ImageJ hyperstack header if OME-XML is absent."""
    meta: Dict[str, Any] = {}
    if ij_meta:
        for key, val in ij_meta.items():
            meta[key] = val
        if "spacing" in ij_meta and "PhysicalSizeZ" not in meta:
            meta["PhysicalSizeZ"] = ij_meta["spacing"]
        if "finterval" in ij_meta and "TimeIncrement" not in meta:
            meta["TimeIncrement"] = ij_meta["finterval"]

    # Extract XY pixel size from standard TIFF resolution tags (pixels per micron in ImageJ).
    try:
        page0 = tif.pages[0]
        if "XResolution" in page0.tags:
            num, den = page0.tags["XResolution"].value
            if num and den:
                meta.setdefault("PhysicalSizeX", round(float(den) / float(num), 6))
        if "YResolution" in page0.tags:
            num, den = page0.tags["YResolution"].value
            if num and den:
                meta.setdefault("PhysicalSizeY", round(float(den) / float(num), 6))
    except Exception:
        pass

    return meta


def _build_imagej_write_kwargs(
    arr: np.ndarray,
    axes: str,
    metadata: Optional[Dict[str, Any]],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Prepares ImageJ-compatible metadata and resolution tags without losing custom keys."""
    ij_meta: Dict[str, Any] = {"axes": axes}
    extra_kwargs: Dict[str, Any] = {}
    if not metadata:
        return ij_meta, extra_kwargs

    for k, v in metadata.items():
        if isinstance(v, (str, int, float, bool)):
            ij_meta[k] = v

    # Map pipeline calibration keys (dx_um/dz_um/dt_s or PhysicalSizeX/Z/TimeIncrement) to ImageJ tags.
    dx = metadata.get("dx_um") or metadata.get("PhysicalSizeX")
    dy = metadata.get("dy_um") or metadata.get("PhysicalSizeY") or dx
    dz = metadata.get("dz_um") or metadata.get("PhysicalSizeZ")
    dt = metadata.get("dt_s") or metadata.get("TimeIncrement")

    try:
        if dx is not None and float(dx) > 0 and dy is not None and float(dy) > 0:
            extra_kwargs["resolution"] = (1.0 / float(dx), 1.0 / float(dy))
            ij_meta["unit"] = "um"
        if dz is not None and float(dz) > 0:
            ij_meta["spacing"] = float(dz)
            ij_meta.setdefault("unit", "um")
        if dt is not None and float(dt) > 0:
            ij_meta["finterval"] = float(dt)
    except (TypeError, ValueError):
        pass

    return ij_meta, extra_kwargs


def read_tiff(path: Union[str, Path]) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Reads a 2D-5D TIFF or OME-TIFF file and returns (numpy_array, metadata_dict).

    Parameters
    ----------
    path : str or Path
        Path to the input .tif / .tiff / .ome.tif file.

    Returns
    -------
    arr : np.ndarray
        Image data preserving original dtype and dimensionality.
    meta : dict
        Combined metadata dictionary (OME-XML <Pixels>, ImageJ tags, axis labels,
        and any '.shape.json' sidecar fields).
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"TIFF file not found: {path}")

    meta: Dict[str, Any] = {}
    arr: np.ndarray

    if HAS_TIFFFILE:
        with tifffile.TiffFile(path) as tif:
            arr = tif.asarray()
            if tif.series:
                meta["axes"] = tif.series[0].axes
            if tif.ome_metadata:
                meta.update(_parse_ome_xml(tif.ome_metadata))
            if tif.imagej_metadata:
                ij_parsed = _parse_imagej_metadata(tif.imagej_metadata, tif)
                for k, v in ij_parsed.items():
                    meta.setdefault(k, v)
            # Check JSON-encoded metadata placed in ImageDescription for non-ImageJ writes.
            try:
                desc = tif.pages[0].description
                if desc and desc.lstrip().startswith("{"):
                    parsed_desc = json.loads(desc)
                    if isinstance(parsed_desc, dict):
                        for k, v in parsed_desc.items():
                            meta.setdefault(k, v)
            except Exception:
                pass
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
            orig_shape = saved.get("shape")
            if orig_shape and int(np.prod(orig_shape)) == arr.size:
                arr = arr.reshape(orig_shape)
            orig_dtype = saved.get("dtype")
            if orig_dtype and str(arr.dtype) != orig_dtype:
                arr = arr.astype(np.dtype(orig_dtype))
            for k, v in saved.get("metadata", {}).items():
                meta.setdefault(k, v)
        except Exception:
            pass

    meta.setdefault("axes", _infer_axes(arr.ndim))
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

    Uses ImageJ hyperstack format when compatible (uint8/uint16/float32 and <3.5 GB) so
    files open directly in FIJI and ilastik with calibrated axes, and automatically falls
    back to standard/BigTIFF for 32-bit integer label masks or >3.5 GB stacks.

    Parameters
    ----------
    path : str or Path
        Destination .tif / .ome.tif path (parent directories created automatically).
    arr : np.ndarray
        2D-5D image or label array.
    metadata : dict, optional
        Calibration and provenance dictionary (e.g., dx_um, dz_um, dt_s, PhysicalSizeX).
    axes : str, optional
        Explicit axis order string (e.g., 'TCZYX', 'TZYX', 'ZYX'). Inferred if omitted.
    compression : str, optional
        Optional tifffile compression codec (e.g., 'zlib', 'lzw'). Note that ImageJ
        hyperstacks are written uncompressed if compression is None for maximum speed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.ascontiguousarray(arr)
    axes_str = axes or _infer_axes(arr.ndim)
    meta_out = dict(metadata) if metadata else {}

    # Write to a temporary file in the same directory so os.replace() is atomic on POSIX/Windows.
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.stem}_", suffix=".tmp.tif", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)

    try:
        if HAS_TIFFFILE:
            use_bigtiff = arr.nbytes >= _BIGTIFF_THRESHOLD_BYTES
            use_imagej = (
                (not use_bigtiff)
                and (arr.dtype in _IMAGEJ_DTYPES)
                and (arr.ndim <= 5)
            )

            if use_imagej:
                ij_meta, extra_kwargs = _build_imagej_write_kwargs(arr, axes_str, meta_out)
                # ImageJ hyperstack expects TZCYX order; transpose TCZYX -> TZCYX for strict FIJI compatibility
                # only when writing 5D ImageJ stacks, otherwise keep native shape via non-ImageJ TIFF.
                if arr.ndim == 5 and axes_str == "TCZYX":
                    tifffile.imwrite(
                        tmp_path,
                        arr,
                        bigtiff=False,
                        compression=compression,
                        metadata={"axes": axes_str, **meta_out},
                    )
                else:
                    tifffile.imwrite(
                        tmp_path,
                        arr,
                        imagej=True,
                        compression=compression,
                        metadata=ij_meta,
                        **extra_kwargs,
                    )
            else:
                tifffile.imwrite(
                    tmp_path,
                    arr,
                    bigtiff=use_bigtiff,
                    compression=compression,
                    metadata={"axes": axes_str, **meta_out},
                )
        elif HAS_PIL:
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
                    {"shape": list(arr.shape), "dtype": str(arr.dtype), "axes": axes_str, "metadata": meta_out},
                    f,
                    indent=2,
                )
        else:
            raise ImportError("Neither 'tifffile' nor 'Pillow' is installed; cannot write TIFF files.")

        os.replace(tmp_path, path)
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise