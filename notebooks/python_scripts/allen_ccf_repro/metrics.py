"""Registration quality metrics computed in the atlas plane of each slice."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.ndimage

from . import atlas_io, config


def plane_points(cells: pd.DataFrame, columns=("atlas_x", "atlas_y")) -> np.ndarray:
    subset = cells[list(columns)].dropna()
    return subset.to_numpy(dtype=float)


def _inside_metrics(points: np.ndarray, mask: np.ndarray) -> dict:
    if not len(points) or not mask.any():
        return {
            "inside_pct": float("nan"),
            "median_distance_25um": float("nan"),
            "p90_distance_25um": float("nan"),
        }
    indices = np.rint(points).astype(int)
    valid = (
        (indices[:, 1] >= 0)
        & (indices[:, 1] < mask.shape[0])
        & (indices[:, 0] >= 0)
        & (indices[:, 0] < mask.shape[1])
    )
    indices = indices[valid]
    if not len(indices):
        return {
            "inside_pct": float("nan"),
            "median_distance_25um": float("nan"),
            "p90_distance_25um": float("nan"),
        }
    distance = scipy.ndimage.distance_transform_edt(~mask)
    values = distance[indices[:, 1], indices[:, 0]]
    return {
        "inside_pct": float(100 * np.mean(values == 0)),
        "median_distance_25um": float(np.median(values)),
        "p90_distance_25um": float(np.quantile(values, 0.90)),
    }


def tissue_containment(cells: pd.DataFrame, annotation_plane: np.ndarray) -> float:
    points = plane_points(cells)
    return _inside_metrics(points, annotation_plane > 0)["inside_pct"]


def region_metrics(
    cells: pd.DataFrame,
    annotation_plane: np.ndarray,
    registry: dict,
    allen_name_to_annots: dict,
) -> dict:
    """Per-region containment and distance for every registry region with cells."""
    results = {}
    for name, entry in registry["regions"].items():
        classes = entry.get("cell_classes", [])
        if not classes:
            continue
        ids = config.resolve_region_ids(entry, allen_name_to_annots)
        mask = atlas_io.region_mask(annotation_plane, ids)
        subset = cells[cells["subclass_label_transfer"].isin(classes)]
        points = plane_points(subset)
        summary = _inside_metrics(points, mask)
        summary.update({"n_cells": int(len(points)), "atlas_pixels": int(mask.sum()), "role": entry.get("role")})
        results[name] = summary
    return results


def _paired_displacement(baseline: pd.DataFrame, candidate: pd.DataFrame) -> pd.Series:
    """Per-cell displacement in 25 um voxels for cells present in both fits."""
    joined = (
        baseline[["ccfy", "ccfz"]]
        .join(candidate[["ccfy", "ccfz"]], how="inner", lsuffix="_base", rsuffix="_new")
        .dropna()
    )
    if not len(joined):
        return pd.Series(dtype=float)
    delta = (
        joined[["ccfy_new", "ccfz_new"]].to_numpy() - joined[["ccfy_base", "ccfz_base"]].to_numpy()
    ) / config.CCF_PIXEL_SIZE
    return pd.Series(np.sqrt(np.sum(delta**2, axis=1)), index=joined.index)


def displacement_metrics(baseline: pd.DataFrame, candidate: pd.DataFrame) -> dict:
    """Displacement of shared cells between two registrations, in 25 um voxels."""
    displacement = _paired_displacement(baseline, candidate)
    if not len(displacement):
        return {
            "median_displacement_25um": float("nan"),
            "p90_displacement_25um": float("nan"),
        }
    return {
        "median_displacement_25um": float(displacement.median()),
        "p90_displacement_25um": float(displacement.quantile(0.90)),
    }


def anchor_cell_index(
    cells: pd.DataFrame,
    annotation_plane: np.ndarray,
    registry: dict,
    allen_name_to_annots: dict,
) -> pd.Index:
    """Cells that already sit inside the atlas structure their class belongs to.

    These are the anchors: a later stage is free to move the cells that are still
    outside their structure, since that is the error being corrected, but moving
    the ones that are already right is collateral damage.
    """
    index = pd.Index([])
    for entry in registry["regions"].values():
        classes = entry.get("cell_classes", [])
        if not classes:
            continue
        mask = atlas_io.region_mask(
            annotation_plane, config.resolve_region_ids(entry, allen_name_to_annots)
        )
        if not mask.any():
            continue
        subset = cells[cells["subclass_label_transfer"].isin(classes)][["atlas_x", "atlas_y"]].dropna()
        if not len(subset):
            continue
        columns = np.rint(subset["atlas_x"].to_numpy()).astype(int)
        rows = np.rint(subset["atlas_y"].to_numpy()).astype(int)
        valid = (
            (rows >= 0) & (rows < mask.shape[0]) & (columns >= 0) & (columns < mask.shape[1])
        )
        inside = np.zeros(len(subset), dtype=bool)
        inside[valid] = mask[rows[valid], columns[valid]]
        index = index.union(subset.index[inside])
    return index


def anchor_displacement(
    baseline: pd.DataFrame,
    candidate: pd.DataFrame,
    anchors: pd.Index,
) -> dict:
    """How far the anchor cells moved between two fits of the same atlas plane."""
    displacement = _paired_displacement(baseline, candidate)
    shared = displacement.index.intersection(anchors)
    if not len(shared):
        return {
            "anchor_median_displacement_25um": float("nan"),
            "anchor_p90_displacement_25um": float("nan"),
            "n_anchor_cells": 0,
        }
    subset = displacement.loc[shared]
    return {
        "anchor_median_displacement_25um": float(subset.median()),
        "anchor_p90_displacement_25um": float(subset.quantile(0.90)),
        "n_anchor_cells": int(len(subset)),
    }


def boundary_pileup(cells: pd.DataFrame, medial_limit: int, protected_classes=()) -> dict:
    """Medial boundary diagnostics: clamping, overhang, and protected regions."""
    columns = cells["atlas_x"].dropna().to_numpy(dtype=float)
    result = {
        "medial_limit": int(medial_limit),
        "boundary_pileup_pct": atlas_io.boundary_pileup_pct(columns, medial_limit),
        "boundary_clamp_ratio": atlas_io.boundary_clamp_ratio(columns, medial_limit),
        "max_atlas_x": float(np.max(columns)) if len(columns) else float("nan"),
        "p995_atlas_x": float(np.quantile(columns, 0.995)) if len(columns) else float("nan"),
    }
    result["medial_overhang_25um"] = result["p995_atlas_x"] - config.ATLAS_MIDLINE
    if len(protected_classes):
        protected = cells[cells["subclass_label_transfer"].isin(list(protected_classes))]
        protected_columns = protected["atlas_x"].dropna().to_numpy(dtype=float)
        result["protected_boundary_pileup_pct"] = atlas_io.boundary_pileup_pct(
            protected_columns, medial_limit
        )
        result["protected_boundary_clamp_ratio"] = atlas_io.boundary_clamp_ratio(
            protected_columns, medial_limit
        )
    return result


def structure_thickness(warped_annotation: np.ndarray, ids, column_mask: np.ndarray = None) -> float:
    """Mean per-column pixel count of a structure in warped atlas space."""
    mask = atlas_io.region_mask(warped_annotation, ids)
    counts = mask.sum(axis=0).astype(float)
    if column_mask is not None:
        counts = counts[column_mask]
    counts = counts[counts > 0]
    return float(counts.mean()) if counts.size else float("nan")


def exclusion_columns(exclusion: np.ndarray, pad: int = 12) -> np.ndarray:
    """Columns that sit under or beside an excluded region."""
    columns = exclusion.any(axis=0)
    if not columns.any():
        return columns
    return scipy.ndimage.binary_dilation(columns, np.ones(2 * pad + 1, dtype=bool))


def summarise(prefix: str, values: dict) -> dict:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def flatten_region_metrics(region_results: dict) -> dict:
    flat = {}
    for name, values in region_results.items():
        for key, value in values.items():
            if key == "role":
                continue
            flat[f"{name}_{key}"] = value
    return flat
