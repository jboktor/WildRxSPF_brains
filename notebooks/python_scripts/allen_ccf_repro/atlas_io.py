"""Atlas loading with explicit, per-slice medial boundary handling.

The historical pipeline hard-coded `[:, :, :236]`, i.e. the atlas midline (228)
plus 8 voxels (200 um) of contralateral tissue. Sections that carry more than
200 um of contralateral tissue had nowhere to put those cells, so medial
structures (medial habenula, third ventricle, dentate medial blade) piled up
against the crop edge. Here the medial margin is a parameter, it is chosen per
slice from the measured contralateral overhang, and the pile-up is measured.
"""

from __future__ import annotations

from dataclasses import dataclass

import nrrd
import numpy as np
import skimage.io

from . import config


@dataclass
class Atlas:
    """Full-resolution 25 um atlas volumes plus the medial crop in use."""

    annotation: np.ndarray
    nissl: np.ndarray
    borders: np.ndarray
    midline: int = config.ATLAS_MIDLINE

    @property
    def shape(self) -> tuple:
        return self.annotation.shape

    def medial_limit(self, medial_extra: int) -> int:
        return int(min(self.annotation.shape[2], self.midline + int(medial_extra)))

    def coronal(self, slice_index: int, medial_extra: int) -> dict:
        """Return one coronal atlas slice cropped to the requested medial margin."""
        limit = self.medial_limit(medial_extra)
        index = int(np.clip(slice_index, 0, self.annotation.shape[0] - 1))
        return {
            "annotation": self.annotation[index, :, :limit],
            "nissl": self.nissl[index, :, :limit],
            "borders": self.borders[index, :, :limit],
            "slice_index": index,
            "medial_limit": limit,
            "column_slice_index": np.full(limit, index, dtype=int),
        }

    def oblique(self, slice_index: int, gradient: float, medial_extra: int) -> dict:
        """Return an oblique coronal slice whose AP level varies with the ML column.

        `gradient` is in atlas slices per 25 um column, positive meaning the
        section runs more posterior towards the midline.
        """
        limit = self.medial_limit(medial_extra)
        columns = np.arange(limit)
        centre = self.midline / 2.0
        levels = np.rint(slice_index + gradient * (columns - centre)).astype(int)
        levels = np.clip(levels, 0, self.annotation.shape[0] - 1)
        rows = np.arange(self.annotation.shape[1])
        row_grid = rows[:, None]
        column_grid = columns[None, :]
        level_grid = levels[None, :]
        return {
            "annotation": self.annotation[level_grid, row_grid, column_grid],
            "nissl": self.nissl[level_grid, row_grid, column_grid],
            "borders": self.borders[level_grid, row_grid, column_grid],
            "slice_index": int(slice_index),
            "medial_limit": limit,
            "column_slice_index": levels,
        }

    def plane(self, slice_index: int, gradient: float, medial_extra: int) -> dict:
        if abs(float(gradient)) < 1e-9:
            return self.coronal(slice_index, medial_extra)
        return self.oblique(slice_index, gradient, medial_extra)


_ATLAS: Atlas | None = None


def load_atlas() -> Atlas:
    """Load the 25 um volumes once per process."""
    global _ATLAS
    if _ATLAS is not None:
        return _ATLAS
    annotation, _ = nrrd.read(str(config.ANNOTATION_25), index_order="F")
    nissl, _ = nrrd.read(str(config.NISSL_25), index_order="F")
    borders = skimage.io.imread(str(config.BORDERS_25))
    _ATLAS = Atlas(
        annotation=np.ascontiguousarray(annotation),
        nissl=np.ascontiguousarray(nissl),
        borders=np.ascontiguousarray(borders),
    )
    return _ATLAS


def load_annotation_10() -> np.ndarray:
    annotation, _ = nrrd.read(str(config.ANNOTATION_10), index_order="F")
    return annotation


def region_mask(annotation_slice: np.ndarray, ids) -> np.ndarray:
    if not len(ids):
        return np.zeros(annotation_slice.shape, dtype=bool)
    return np.isin(annotation_slice, np.asarray(list(ids), dtype=annotation_slice.dtype))


def medial_extra_from_overhang(cell_columns_25um: np.ndarray, midline: int = config.ATLAS_MIDLINE) -> int:
    """Pick the smallest medial margin that still contains the contralateral tissue."""
    if not len(cell_columns_25um):
        return config.MIN_MEDIAL_EXTRA
    overhang = float(np.quantile(cell_columns_25um, 0.995)) - midline
    extra = int(np.ceil(max(overhang, 0.0))) + config.MEDIAL_EXTRA_MARGIN
    return int(np.clip(extra, config.MIN_MEDIAL_EXTRA, config.MAX_MEDIAL_EXTRA))


def boundary_pileup_pct(cell_columns_25um: np.ndarray, medial_limit: int, tolerance: int = 0) -> float:
    """Fraction of cells sitting in the last atlas column before the crop edge."""
    if not len(cell_columns_25um):
        return float("nan")
    return float(100 * np.mean(np.asarray(cell_columns_25um) >= (medial_limit - 1 - tolerance)))


CLAMP_MIN_CELLS = 10


def boundary_clamp_ratio(cell_columns_25um: np.ndarray, medial_limit: int, window: int = 10) -> float:
    """Density in the edge column relative to the columns just lateral to it.

    Medial structures legitimately reach the midline, so the fraction of cells
    near the edge says more about where a structure lives than about the fit.
    Clamping is a different signature: cells that belong past the crop edge are
    stacked into the final column, so its density spikes above its neighbours.
    A ratio near one means the cells simply thin out normally at the edge.
    """
    columns = np.asarray(cell_columns_25um, dtype=float)
    if not len(columns):
        return float("nan")
    edge = medial_limit - 1
    at_edge = float(np.sum(columns >= edge))
    neighbours = np.array(
        [np.sum((columns >= edge - step - 1) & (columns < edge - step)) for step in range(window)],
        dtype=float,
    )
    reference = float(np.median(neighbours))
    if at_edge < CLAMP_MIN_CELLS:
        return 0.0
    if reference <= 0:
        return float("inf")
    return at_edge / reference
