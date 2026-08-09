"""Cell table access and density-based inlier filtering."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.ndimage

from . import config


def read_sample(sample_id: str, path=None) -> pd.DataFrame:
    """Stream the large baseline table and keep only one sample."""
    path = config.BASELINE_CCF if path is None else path
    chunks = []
    for chunk in pd.read_csv(path, index_col=0, chunksize=250_000):
        selected = chunk[chunk["sample_id"].eq(sample_id)]
        if len(selected):
            chunks.append(selected.copy())
    if not chunks:
        raise RuntimeError(f"No cells found for sample {sample_id}")
    return pd.concat(chunks, axis=0)


def read_slice_parameters(sample_id: str) -> tuple:
    """Return the parameter workbook and the row index for one sample."""
    frame = pd.read_excel(config.SOURCE_WORKBOOK, sheet_name="slice_parameters")
    frame = frame.reset_index(drop=True)
    matches = frame.index[frame["cell_metadata"].eq(sample_id)]
    if not len(matches):
        raise RuntimeError(f"Sample {sample_id} absent from {config.SOURCE_WORKBOOK}")
    return frame, int(matches[0])


def class_points(cells: pd.DataFrame, classes, columns=("fixed_x", "fixed_y")) -> np.ndarray:
    """Points for a set of cell classes in fixed-image space."""
    if not len(classes):
        return np.zeros((0, 2))
    subset = cells[cells["subclass_label_transfer"].isin(list(classes))]
    subset = subset[list(columns)].dropna()
    return subset.to_numpy(dtype=float)


def density_inliers(
    points: np.ndarray,
    shape: tuple,
    sigma: float = 6.0,
    keep_quantile: float = 0.25,
    min_component_fraction: float = 0.08,
) -> np.ndarray:
    """Drop label-transfer outliers by keeping cells inside dense components.

    Cells of a regional class are expected to form one or a few compact clouds.
    Isolated cells scattered through cortex are transfer noise, and they used to
    drag landmark centroids away from the true structure.
    """
    if len(points) < 10:
        return np.ones(len(points), dtype=bool)
    grid = np.zeros(shape, dtype=float)
    rows = np.clip(np.rint(points[:, 1]).astype(int), 0, shape[0] - 1)
    columns = np.clip(np.rint(points[:, 0]).astype(int), 0, shape[1] - 1)
    np.add.at(grid, (rows, columns), 1.0)
    density = scipy.ndimage.gaussian_filter(grid, sigma)
    values = density[rows, columns]
    threshold = np.quantile(values, keep_quantile)
    dense = density >= max(threshold, 1e-12)
    labels, count = scipy.ndimage.label(dense)
    if count == 0:
        return np.ones(len(points), dtype=bool)
    membership = labels[rows, columns]
    sizes = np.bincount(membership, minlength=count + 1).astype(float)
    sizes[0] = 0
    keep_labels = np.where(sizes >= min_component_fraction * len(points))[0]
    if not len(keep_labels):
        keep_labels = np.array([int(np.argmax(sizes))])
    return np.isin(membership, keep_labels)


def cell_density_image(points: np.ndarray, shape: tuple, sigma: float = 3.0) -> np.ndarray:
    grid = np.zeros(shape, dtype=float)
    if not len(points):
        return grid
    rows = np.clip(np.rint(points[:, 1]).astype(int), 0, shape[0] - 1)
    columns = np.clip(np.rint(points[:, 0]).astype(int), 0, shape[1] - 1)
    np.add.at(grid, (rows, columns), 1.0)
    return scipy.ndimage.gaussian_filter(grid, sigma)
