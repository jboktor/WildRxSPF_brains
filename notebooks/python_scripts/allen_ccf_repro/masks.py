"""Tissue, damage, and artifact masks in cropped-and-padded fixed-image space.

Three separate ideas are kept distinct on purpose:

* blurred-but-present tissue must stay in the metric, otherwise the atlas is
  pulled off good anatomy (roi1_run1 dorsal cortex);
* genuinely missing tissue must leave the metric, otherwise the atlas is
  compressed into empty space (the dorsal cuts in roi2_run4, roi3_run1,
  roi3_run3, roi4_run1, roi4_run2, roi4_run3);
* bright non-tissue artifacts must leave the metric in both directions.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.ndimage
import skimage.draw
import skimage.filters
import skimage.measure
import skimage.morphology
import yaml

from . import cell_io, config


def load_overrides(sample_id: str) -> dict:
    path = config.MASK_OVERRIDE_ROOT / f"{sample_id}.yaml"
    if not path.exists():
        return {}
    with open(path) as handle:
        data = yaml.safe_load(handle) or {}
    return data


def _rasterise(shape: tuple, rectangles=(), polygons=()) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    for x0, y0, x1, y1 in rectangles or ():
        rows = slice(max(int(y0), 0), min(int(y1), shape[0]))
        columns = slice(max(int(x0), 0), min(int(x1), shape[1]))
        mask[rows, columns] = True
    for polygon in polygons or ():
        points = np.asarray(polygon, dtype=float)
        rows, columns = skimage.draw.polygon(points[:, 1], points[:, 0], shape)
        mask[rows, columns] = True
    return mask


def cell_supported_tissue_mask(
    fixed: np.ndarray,
    cells,
    dapi_quantile: float = 0.55,
    density_threshold: float = 0.012,
    min_object_px: int = 400,
) -> dict:
    """Tissue mask from DAPI signal unioned with cell-supported territory.

    Cell positions come from the independent segmentation, so regions that are
    optically blurred but biologically present still enter the mask.
    """
    image = np.asarray(fixed, dtype=float)
    positive = image[image > 0]
    if positive.size:
        try:
            level = float(skimage.filters.threshold_otsu(positive))
        except ValueError:
            level = float(np.quantile(positive, dapi_quantile))
    else:
        level = 0.0
    blurred = skimage.filters.gaussian(image, 2.0, preserve_range=True)
    dapi_mask = blurred > (0.5 * level)

    points = cells[["fixed_x", "fixed_y"]].dropna().to_numpy(dtype=float)
    density = cell_io.cell_density_image(points, image.shape, sigma=6.0)
    scale = float(np.quantile(density[density > 0], 0.60)) if np.any(density > 0) else 0.0
    cell_mask = density > max(density_threshold * max(scale, 1e-9), 1e-9)

    tissue = np.logical_or(dapi_mask, cell_mask)
    tissue = skimage.morphology.binary_closing(tissue, skimage.morphology.disk(5))
    tissue = scipy.ndimage.binary_fill_holes(tissue)
    tissue = skimage.morphology.remove_small_objects(tissue, min_object_px)
    return {
        "tissue": tissue,
        "dapi_only": dapi_mask,
        "cell_supported": np.logical_and(cell_mask, np.logical_not(dapi_mask)),
        "otsu_level": level,
    }


def artifact_mask(
    fixed: np.ndarray,
    tissue: np.ndarray,
    cells,
    overrides: dict,
    bright_quantile: float = 0.995,
    min_area: int = 250,
) -> np.ndarray:
    """Bright blocks with no cells are imaging artifacts, not tissue."""
    image = np.asarray(fixed, dtype=float)
    inside = image[tissue]
    if inside.size:
        level = float(np.quantile(inside, bright_quantile))
    else:
        level = float(image.max())
    points = cells[["fixed_x", "fixed_y"]].dropna().to_numpy(dtype=float)
    density = cell_io.cell_density_image(points, image.shape, sigma=6.0)
    empty = density <= (0.02 * float(density.max()) if density.max() > 0 else 0.0)

    candidate = np.logical_and(image > level, empty)
    candidate = skimage.morphology.binary_closing(candidate, skimage.morphology.disk(3))
    candidate = skimage.morphology.remove_small_objects(candidate, min_area)

    detected = np.zeros(image.shape, dtype=bool)
    for region in skimage.measure.regionprops(skimage.measure.label(candidate)):
        if region.area < min_area:
            continue
        if region.solidity < 0.55:
            continue
        detected[tuple(region.coords.T)] = True

    manual = _rasterise(
        image.shape,
        overrides.get("artifact_rects", ()),
        overrides.get("artifact_polygons", ()),
    )
    return np.logical_or(detected, manual)


def missing_tissue_exclusion(
    tissue: np.ndarray,
    atlas_tissue_in_fixed: np.ndarray,
    cells,
    overrides: dict,
    halo_px: int = 6,
    min_area: int = 500,
) -> dict:
    """Atlas territory with neither DAPI signal nor cells is missing tissue.

    Both conditions are required: blurred tissue keeps its cells and therefore
    never enters the exclusion.
    """
    shape = tissue.shape
    points = cells[["fixed_x", "fixed_y"]].dropna().to_numpy(dtype=float)
    density = cell_io.cell_density_image(points, shape, sigma=6.0)
    has_cells = density > (0.02 * float(density.max()) if density.max() > 0 else 0.0)

    covered = skimage.morphology.binary_dilation(
        np.logical_or(tissue, has_cells), skimage.morphology.disk(halo_px)
    )
    candidate = np.logical_and(atlas_tissue_in_fixed, np.logical_not(covered))
    candidate = skimage.morphology.remove_small_objects(candidate, min_area)

    auto = np.zeros(shape, dtype=bool)
    components = []
    for region in skimage.measure.regionprops(skimage.measure.label(candidate)):
        if region.area < min_area:
            continue
        auto[tuple(region.coords.T)] = True
        components.append(
            {
                "area_px": int(region.area),
                "bbox": [int(value) for value in region.bbox],
                "centroid_xy": [float(region.centroid[1]), float(region.centroid[0])],
            }
        )

    manual = _rasterise(
        shape,
        overrides.get("exclusion_rects", ()),
        overrides.get("exclusion_polygons", ()),
    )
    keep = _rasterise(
        shape,
        overrides.get("keep_rects", ()),
        overrides.get("keep_polygons", ()),
    )
    exclusion = np.logical_and(np.logical_or(auto, manual), np.logical_not(keep))
    exclusion = skimage.morphology.binary_dilation(exclusion, skimage.morphology.disk(halo_px))
    exclusion = np.logical_and(exclusion, np.logical_not(tissue))
    return {
        "exclusion": exclusion,
        "components": components,
        "manual_used": bool(manual.any()),
        "area_px": int(exclusion.sum()),
    }


def registration_mask(
    tissue: np.ndarray,
    exclusion: np.ndarray,
    artifact: np.ndarray,
    grow_px: int = 0,
) -> np.ndarray:
    """Fixed mask for elastix: the whole field of view minus the dead zones.

    Everything outside the dead zones is kept, including the background next to
    intact tissue: that background is what pulls the atlas surface onto the real
    tissue surface, and removing it measurably degrades the fit. Excluded zones
    remove that pull exactly where tissue is missing, so the atlas is not
    compressed into the cut. When a slice has no dead zones, `None` is returned
    and elastix runs unmasked, which is the historical behaviour.
    """
    dead = np.logical_or(exclusion, artifact)
    if grow_px > 0:
        dead = skimage.morphology.binary_dilation(dead, skimage.morphology.disk(grow_px))
    elif grow_px < 0:
        dead = skimage.morphology.binary_erosion(dead, skimage.morphology.disk(-grow_px))
    if not dead.any():
        return None
    return np.logical_not(dead).astype(np.uint8)


def write_override_template(sample_id: str, notes: str = "") -> Path:
    path = config.MASK_OVERRIDE_ROOT / f"{sample_id}.yaml"
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "sample_id": sample_id,
        "notes": notes,
        "artifact_rects": [],
        "exclusion_rects": [],
        "keep_rects": [],
    }
    with open(path, "w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
    return path
