"""Shared per-slice machinery: fixed image, atlas plane, and one registration fit."""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pandas as pd

from . import config, elastix_stage, nb_functions


def parse_list(value) -> list:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    text = str(value)
    if text.strip().lower() in ("nan", "none", ""):
        return []
    parsed = ast.literal_eval(text)
    return [] if parsed is None else list(parsed)


def prepare_fixed(frame: pd.DataFrame, row_index: int, cells_raw: pd.DataFrame) -> dict:
    """Enhanced DAPI image plus cell positions in cropped-and-padded fixed space."""
    row = frame.iloc[row_index]
    cell_types = parse_list(row["cell_types"])
    fixed_full = nb_functions.modify_dapi(
        frame,
        row_index,
        cells_raw,
        cell_types=cell_types or None,
        factor=float(row["dapi_enhance_factor"]),
    )
    bbox, fixed, _ = nb_functions.crop_and_pad_image(
        fixed_full, pad_width=config.PAD_WIDTH, area_thresh=float(row["area_thresh"])
    )
    cells = nb_functions.get_cell_metadata_for_slice_index(
        frame,
        row_index,
        cells_raw,
        ccf_pixel_size=config.CCF_PIXEL_SIZE,
        bbox=bbox,
    )
    return {
        "image": fixed,
        "bbox": bbox,
        "cells": cells,
        "shape": fixed.shape,
        "cell_types": cell_types,
        "row": row,
    }


def prepare_moving(atlas, frame: pd.DataFrame, row_index: int, slice_index: int, gradient: float, medial_extra: int) -> dict:
    """Enhanced Allen Nissl plane plus the matching annotation and borders."""
    row = frame.iloc[row_index]
    plane = atlas.plane(slice_index, gradient, medial_extra)
    annotations = parse_list(row["annots_to_amplify"])
    moving_full = nb_functions.modify_nissl(
        plane["nissl"],
        plane["annotation"],
        factor=float(row["nissl_enhance_factor"]),
        annot_dict=annotations or None,
    )
    area_thresh = float(row["area_thresh"])
    _, _, borders = nb_functions.crop_and_pad_image(
        moving_full, plane["borders"], pad_width=config.PAD_WIDTH, area_thresh=area_thresh
    )
    bbox, moving, annotation = nb_functions.crop_and_pad_image(
        moving_full, plane["annotation"], pad_width=config.PAD_WIDTH, area_thresh=area_thresh
    )
    return {
        "image": moving,
        "annotation": annotation,
        "borders": borders,
        "bbox": bbox,
        "plane": plane,
        "slice_index": int(slice_index),
        "gradient": float(gradient),
        "medial_extra": int(medial_extra),
    }


def slice_parameters(row, cfg) -> dict:
    """Per-slice elastix settings, preferring the values tuned in the sheet.

    `spline_grid_size`, `iterations`, and `histogram_bins` were tuned per slice
    for the historical fits and remain the best available prior; the config
    values are only fallbacks for slices that never got one.
    """

    def value(name, fallback):
        if name not in row.index:
            return fallback
        entry = row[name]
        return fallback if entry is None or (isinstance(entry, float) and np.isnan(entry)) else entry

    return {
        "grid_spacing_voxels": float(value("spline_grid_size", cfg.spline_grid_voxels)),
        "spline_iterations": int(value("iterations", cfg.spline_iterations)),
        "histogram_bins": int(value("histogram_bins", 32)),
        "affine_iterations": cfg.affine_iterations,
        "bending_weight": cfg.bending_energy_weight,
    }


def plane_coordinates(mapped_points: np.ndarray, moving_bbox) -> np.ndarray:
    """Convert transformix output into atlas-plane pixel coordinates."""
    points = np.asarray(mapped_points, dtype=float).copy()
    points -= config.PAD_WIDTH
    points += np.array([moving_bbox[1], moving_bbox[0]], dtype=float)
    return points


def cells_in_atlas(cells: pd.DataFrame, plane_points: np.ndarray, moving_ctx: dict) -> pd.DataFrame:
    """Attach atlas-plane and CCF coordinates to each cell."""
    plane = moving_ctx["plane"]
    annotation = plane["annotation"]
    columns = plane["column_slice_index"]
    frame = cells.copy()
    frame["atlas_x"] = plane_points[:, 0]
    frame["atlas_y"] = plane_points[:, 1]

    x_index = np.rint(frame["atlas_x"].to_numpy()).astype(int)
    y_index = np.rint(frame["atlas_y"].to_numpy()).astype(int)
    inside = (
        (x_index >= 0)
        & (x_index < annotation.shape[1])
        & (y_index >= 0)
        & (y_index < annotation.shape[0])
    )
    slice_index = np.full(len(frame), np.nan)
    slice_index[inside] = columns[x_index[inside]]
    annotation_value = np.full(len(frame), np.nan)
    annotation_value[inside] = annotation[y_index[inside], x_index[inside]]

    frame["atlas_slice_index"] = slice_index
    frame["annotation"] = annotation_value
    frame["ccfx"] = config.CCF_PIXEL_SIZE * slice_index
    frame["ccfy"] = np.where(inside, config.CCF_PIXEL_SIZE * frame["atlas_y"], np.nan)
    frame["ccfz"] = np.where(inside, config.CCF_PIXEL_SIZE * frame["atlas_x"], np.nan)
    frame.loc[~inside, ["atlas_x", "atlas_y"]] = np.nan
    return frame


def fit(
    fixed_ctx: dict,
    moving_ctx: dict,
    work_dir: Path,
    stage: str = "full",
    fixed_mask: np.ndarray = None,
    moving_mask: np.ndarray = None,
    landmark_points: dict = None,
    landmark_weight: float = 0.0,
    bending_weight: float = 0.0,
    grid_spacing_voxels: float = 20,
    affine_iterations: int = 500,
    spline_iterations: int = 800,
    histogram_bins: int = 32,
    use_rigid: bool = True,
) -> dict:
    """Run one registration and map every cell into atlas-plane coordinates."""
    fixed_points = None
    moving_points = None
    if landmark_points is not None and len(landmark_points.get("fixed", [])):
        fixed_points = landmark_points["fixed"]
        moving_points = landmark_points["moving"]
    else:
        landmark_weight = 0.0

    parameter_maps = elastix_stage.build_parameter_maps(
        stage=stage,
        iterations_affine=affine_iterations,
        iterations_spline=spline_iterations,
        grid_spacing_voxels=grid_spacing_voxels,
        histogram_bins=histogram_bins,
        landmark_weight=landmark_weight,
        bending_weight=bending_weight if stage == "full" else 0.0,
        use_rigid=use_rigid,
    )
    transform, warped_moving = elastix_stage.register(
        fixed_ctx["image"],
        moving_ctx["image"],
        parameter_maps,
        work_dir=work_dir,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
        fixed_points=fixed_points,
        moving_points=moving_points,
    )
    mapped = elastix_stage.transform_points(
        fixed_ctx["cells"][["fixed_x", "fixed_y"]].to_numpy(dtype=float),
        transform,
        moving_ctx["image"],
        work_dir,
    )
    plane_points = plane_coordinates(mapped, moving_ctx["bbox"])
    cells_atlas = cells_in_atlas(fixed_ctx["cells"], plane_points, moving_ctx)
    return {
        "transform": transform,
        "warped_moving": warped_moving,
        "cells": cells_atlas,
        "landmark_weight": landmark_weight,
        "n_landmarks": 0 if fixed_points is None else int(len(fixed_points)),
    }


def warped_annotation_in_fixed(moving_ctx: dict, transform) -> np.ndarray:
    return elastix_stage.transform_labels(moving_ctx["annotation"].astype(np.int64), transform)


def warped_borders_in_fixed(moving_ctx: dict, transform) -> np.ndarray:
    return elastix_stage.transform_image(moving_ctx["borders"], transform)


def update_baseline_table(cells_raw: pd.DataFrame, cells_atlas: pd.DataFrame) -> pd.DataFrame:
    """Write the new coordinates back onto the full sample table."""
    updated = cells_raw.copy()
    for column in ("ccfx", "ccfy", "ccfz", "annotation"):
        updated.loc[:, column] = np.nan
    shared = cells_atlas.index.intersection(updated.index)
    for column in ("ccfx", "ccfy", "ccfz", "annotation"):
        updated.loc[shared, column] = cells_atlas.loc[shared, column]
    return updated
