#!/usr/bin/env python
"""Run guarded regional-landmark registration candidates for all CCF slices.

All mutable products are isolated under ``Allen_CCF_regional_tests`` and
``figures/Allen_CCF_alignment_optimized/iteration_11|12``.  Canonical
parameter workbooks and CCF coordinate tables are read-only inputs.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import os
import pickle
import shutil
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nrrd
import numpy as np
import pandas as pd
import scipy.ndimage
import skimage.io
import skimage.measure


BASE = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
REGIONAL_SCRIPT = SCRIPT_DIR / "allen_ccf_regional_parameter_sweep.py"
TEST_ROOT = BASE / "data/interim/registration/Allen_CCF_regional_tests"
WORKBOOK = TEST_ROOT / "slice_positions_25um_all_slices_regional_landmarks.xlsx"
SOURCE_WORKBOOK = (
    BASE
    / "data/interim/registration/Allen_CCF_optimization/slice_positions_25um_final.xlsx"
)
BASELINE_CCF = BASE / "data/interim/registration/Allen_CCF_optimization/ccf3d.csv"
REFERENCE_ROOT = BASE / "data/input/allen_registration_ref"
FIGURE_ROOT = BASE / "figures/Allen_CCF_alignment_optimized"
METRICS_CSV = TEST_ROOT / "all_slices_regional_landmark_metrics.csv"
DECISIONS_CSV = TEST_ROOT / "all_slices_regional_landmark_decisions.csv"
ELIGIBILITY_CSV = TEST_ROOT / "all_slices_regional_landmark_eligibility.csv"
BASELINE_METRICS_CSV = TEST_ROOT / "all_slices_regional_baseline_metrics.csv"

TARGET_CLASSES = {
    "DG": "037 DG Glut",
    "MH": "145 MH Tac2 Glut",
    "LH": "146 LH Pou4f1 Sox1 Glut",
    "Ependymal": "323 Ependymal NN",
}
TARGET_IDS = {
    "DG": [632, 758, 790, 823],
    "MH": [483],
    "LH": [186],
    "Ependymal": [73, 81, 89, 98, 108, 116, 124, 129, 140, 145, 153, 164],
}
MIN_CELLS = {"DG": 80, "MH": 20, "LH": 20}
MIN_ATLAS_PIXELS = {"DG": 100, "MH": 1, "LH": 1}
ITERATIONS = {0.03: 11, 0.08: 12}


def load_regional_module():
    """Import the prior sweep and load its notebook-defined helpers."""
    spec = importlib.util.spec_from_file_location("allen_ccf_regional_sweep", REGIONAL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {REGIONAL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.load_notebook_functions()
    with open(REFERENCE_ROOT / "allen_name_to_annots.pkl", "rb") as handle:
        module.allen_name_to_annots = pickle.load(handle)
    return module


def parse_list(value: Any) -> list:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple, np.ndarray)):
        return list(value)
    parsed = ast.literal_eval(str(value))
    return [] if parsed is None else list(parsed)


def load_baseline_cells(samples: set[str] | None = None) -> pd.DataFrame:
    """Stream the baseline once, retaining the 13 workbook-matched samples."""
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(BASELINE_CCF, index_col=0, chunksize=125_000):
        if samples is not None:
            chunk = chunk[chunk["sample_id"].isin(samples)]
        if len(chunk):
            chunks.append(chunk)
    if not chunks:
        raise RuntimeError("No baseline CCF cells matched the requested samples")
    return pd.concat(chunks, axis=0)


def safe_target_metrics(cells: pd.DataFrame, annotation_25: np.ndarray) -> dict[str, Any]:
    """Measure regional fit without failing on absent classes or atlas regions."""
    valid = cells.dropna(subset=["ccfx", "ccfy", "ccfz"])
    result: dict[str, Any] = {}
    if valid.empty:
        result["atlas_slice"] = np.nan
        result["all_tissue_containment_pct"] = np.nan
        for label in TARGET_CLASSES:
            result.update(
                {
                    f"{label}_n": 0,
                    f"{label}_atlas_pixels": 0,
                    f"{label}_median_distance_25um": np.nan,
                    f"{label}_p90_distance_25um": np.nan,
                    f"{label}_inside_pct": np.nan,
                }
            )
        return result

    atlas_z = int(np.clip(np.rint(np.median(valid["ccfx"]) / 25), 0, annotation_25.shape[0] - 1))
    atlas_slice = annotation_25[atlas_z]
    result["atlas_slice"] = atlas_z
    all_points = np.rint(valid[["ccfy", "ccfz"]].to_numpy() / 25).astype(int)
    in_bounds = (
        (all_points[:, 0] >= 0)
        & (all_points[:, 0] < atlas_slice.shape[0])
        & (all_points[:, 1] >= 0)
        & (all_points[:, 1] < atlas_slice.shape[1])
    )
    all_points = all_points[in_bounds]
    result["all_tissue_containment_pct"] = (
        float(100 * np.mean(atlas_slice[all_points[:, 0], all_points[:, 1]] > 0))
        if len(all_points)
        else np.nan
    )

    for label, cell_class in TARGET_CLASSES.items():
        group = valid[valid["subclass_label_transfer"].eq(cell_class)]
        points = np.rint(group[["ccfy", "ccfz"]].to_numpy() / 25).astype(int)
        point_valid = (
            (points[:, 0] >= 0)
            & (points[:, 0] < atlas_slice.shape[0])
            & (points[:, 1] >= 0)
            & (points[:, 1] < atlas_slice.shape[1])
        )
        points = points[point_valid]
        mask = np.isin(atlas_slice, TARGET_IDS[label])
        result[f"{label}_n"] = int(len(points))
        result[f"{label}_atlas_pixels"] = int(mask.sum())
        if not len(points) or not mask.any():
            result[f"{label}_median_distance_25um"] = np.nan
            result[f"{label}_p90_distance_25um"] = np.nan
            result[f"{label}_inside_pct"] = np.nan
            continue
        distance = scipy.ndimage.distance_transform_edt(~mask)
        distances = distance[points[:, 0], points[:, 1]]
        result[f"{label}_median_distance_25um"] = float(np.median(distances))
        result[f"{label}_p90_distance_25um"] = float(np.quantile(distances, 0.90))
        result[f"{label}_inside_pct"] = float(100 * np.mean(distances == 0))
    return result


def landmark_eligibility(metrics: dict[str, Any]) -> dict[str, Any]:
    """Determine supported landmarks and whether regional correction is needed."""
    eligible: dict[str, bool] = {}
    reasons: list[str] = []
    for label in ("DG", "MH", "LH"):
        cells = int(metrics.get(f"{label}_n", 0))
        pixels = int(metrics.get(f"{label}_atlas_pixels", 0))
        eligible[label] = cells >= MIN_CELLS[label] and pixels >= MIN_ATLAS_PIXELS[label]
        if not eligible[label]:
            reasons.append(f"{label}:cells={cells},atlas_pixels={pixels}")

    needs = False
    if eligible["DG"]:
        needs |= (
            float(metrics["DG_p90_distance_25um"]) > 3.0
            or float(metrics["DG_inside_pct"]) < 68.0
        )
    for label in ("MH", "LH"):
        if eligible[label] and pd.notna(metrics[f"{label}_inside_pct"]):
            needs |= float(metrics[f"{label}_inside_pct"]) < 60.0
    return {
        "DG_eligible": eligible["DG"],
        "MH_eligible": eligible["MH"],
        "LH_eligible": eligible["LH"],
        "eligible_classes": ",".join(label for label, ok in eligible.items() if ok),
        "needs_regional_optimization": bool(needs),
        "ineligibility_reasons": "; ".join(reasons),
    }


def build_workbook(samples: list[str], baseline_metrics: pd.DataFrame, eligibility: pd.DataFrame) -> None:
    """Clone source sheets and add isolated all-slice candidate parameter sheets."""
    source_sheets = pd.read_excel(SOURCE_WORKBOOK, sheet_name=None)
    params = source_sheets["slice_parameters"]
    params = params[params["cell_metadata"].isin(samples)].copy().reset_index(drop=True)
    if set(params["cell_metadata"]) != set(samples):
        missing = sorted(set(samples) - set(params["cell_metadata"]))
        raise RuntimeError(f"Samples absent from optimization workbook: {missing}")

    candidate_sheets: dict[str, pd.DataFrame] = {}
    eligibility_by_sample = eligibility.set_index("sample_id")
    for weight, sheet_name in ((0.03, "landmarks003"), (0.08, "landmarks008")):
        candidate = params.copy()
        candidate["cell_types"] = candidate["cell_types"].astype(object)
        candidate["annots_to_amplify"] = candidate["annots_to_amplify"].astype(object)
        candidate["landmark_weight"] = weight
        for idx, row in candidate.iterrows():
            sample_id = row["cell_metadata"]
            eligible_classes = eligibility_by_sample.loc[sample_id, "eligible_classes"].split(",")
            eligible_classes = [item for item in eligible_classes if item]
            cell_types = parse_list(row["cell_types"])
            annotations = parse_list(row["annots_to_amplify"])
            for label in eligible_classes:
                cell_class = TARGET_CLASSES[label]
                if cell_class not in cell_types:
                    cell_types.append(cell_class)
                for annotation_id in TARGET_IDS[label]:
                    if annotation_id not in annotations:
                        annotations.append(annotation_id)
            candidate.at[idx, "cell_types"] = repr(cell_types)
            candidate.at[idx, "annots_to_amplify"] = repr(annotations)
            candidate.at[idx, "dapi_enhance_factor"] = 3
            candidate.at[idx, "nissl_enhance_factor"] = 3
            candidate.at[idx, "spline_grid_size"] = 13
        candidate_sheets[sheet_name] = candidate

    TEST_ROOT.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(WORKBOOK, engine="openpyxl", mode="w") as writer:
        for name, frame in source_sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)
        params.to_excel(writer, sheet_name="baseline_copy", index=False)
        for name, frame in candidate_sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)
        baseline_metrics.to_excel(writer, sheet_name="baseline_metrics", index=False)
        eligibility.to_excel(writer, sheet_name="eligibility", index=False)


def _quantile_landmarks(
    fixed_cloud: np.ndarray, moving_cloud: np.ndarray
) -> tuple[list[list[float]], list[list[float]]]:
    """Pair blade quantiles, adapting bin count to sparse or split DG clouds."""
    if len(fixed_cloud) < 20 or len(moving_cloud) < 20:
        return [], []
    for bin_count in (6, 5, 4):
        fixed_points: list[list[float]] = []
        moving_points: list[list[float]] = []
        fixed_edges = np.quantile(fixed_cloud[:, 0], np.linspace(0.08, 0.92, bin_count + 1))
        moving_edges = np.quantile(moving_cloud[:, 0], np.linspace(0.08, 0.92, bin_count + 1))
        for index in range(bin_count):
            fixed_bin = fixed_cloud[
                (fixed_cloud[:, 0] >= fixed_edges[index])
                & (fixed_cloud[:, 0] <= fixed_edges[index + 1])
            ]
            moving_bin = moving_cloud[
                (moving_cloud[:, 0] >= moving_edges[index])
                & (moving_cloud[:, 0] <= moving_edges[index + 1])
            ]
            if len(fixed_bin) < 10 or len(moving_bin) < 10:
                continue
            for y_quantile in (0.20, 0.80):
                fixed_points.append(
                    [float(np.median(fixed_bin[:, 0])), float(np.quantile(fixed_bin[:, 1], y_quantile))]
                )
                moving_points.append(
                    [float(np.median(moving_bin[:, 0])), float(np.quantile(moving_bin[:, 1], y_quantile))]
                )
        if len(fixed_points) >= 6:
            return fixed_points, moving_points
    return [], []


def make_robust_landmarks(
    fixed_cells: pd.DataFrame,
    moving_annotation: np.ndarray,
    eligible: dict[str, bool],
    min_cells: dict[str, int] | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Create reflected-safe and posterior-safe DG/MH/LH correspondences."""
    thresholds = {**MIN_CELLS, **(min_cells or {})}
    fixed_points: list[list[float]] = []
    moving_points: list[list[float]] = []
    classes: list[str] = []
    class_counts: dict[str, int] = {"DG": 0, "MH": 0, "LH": 0}

    if eligible.get("DG", False):
        dg_fixed = fixed_cells[
            fixed_cells["subclass_label_transfer"].eq(TARGET_CLASSES["DG"])
        ][["fixed_x", "fixed_y"]].dropna().to_numpy()
        moving_y, moving_x = np.where(np.isin(moving_annotation, TARGET_IDS["DG"]))
        dg_moving = np.column_stack([moving_x, moving_y])
        if len(dg_fixed) >= thresholds["DG"] and len(dg_moving) >= MIN_ATLAS_PIXELS["DG"]:
            labels, component_count = scipy.ndimage.label(
                np.isin(moving_annotation, TARGET_IDS["DG"])
            )
            component_sizes = np.bincount(labels.ravel())[1:]
            split_posterior = int(np.sum(component_sizes >= 25)) >= 2
            if split_posterior:
                # Only posterior sections have two substantive DG components.
                # Split both spaces independently to prevent cross-component pairing.
                fixed_center_x = float(np.median(fixed_cells["fixed_x"]))
                moving_center_x = moving_annotation.shape[1] / 2
                side_pairs = [
                    (
                        dg_fixed[dg_fixed[:, 0] < fixed_center_x],
                        dg_moving[dg_moving[:, 0] < moving_center_x],
                    ),
                    (
                        dg_fixed[dg_fixed[:, 0] >= fixed_center_x],
                        dg_moving[dg_moving[:, 0] >= moving_center_x],
                    ),
                ]
            else:
                side_pairs = [(dg_fixed, dg_moving)]
            for fixed_side, moving_side in side_pairs:
                fp, mp = _quantile_landmarks(fixed_side, moving_side)
                fixed_points.extend(fp)
                moving_points.extend(mp)
            class_counts["DG"] = len(fixed_points)
            if class_counts["DG"]:
                classes.append("DG")

    for label in ("MH", "LH"):
        if not eligible.get(label, False):
            continue
        group = fixed_cells[
            fixed_cells["subclass_label_transfer"].eq(TARGET_CLASSES[label])
        ][["fixed_x", "fixed_y"]].dropna()
        moving_y, moving_x = np.where(np.isin(moving_annotation, TARGET_IDS[label]))
        if len(group) >= thresholds[label] and len(moving_x) >= MIN_ATLAS_PIXELS[label]:
            fixed_points.append([float(group["fixed_x"].median()), float(group["fixed_y"].median())])
            moving_points.append([float(np.median(moving_x)), float(np.median(moving_y))])
            class_counts[label] = 1
            classes.append(label)

    fixed_array = np.asarray(fixed_points, dtype=float).reshape(-1, 2)
    moving_array = np.asarray(moving_points, dtype=float).reshape(-1, 2)
    diagnostics = {
        "landmark_count": int(len(fixed_array)),
        "landmark_classes": ",".join(classes),
        "DG_landmark_count": class_counts["DG"],
        "MH_landmark_count": class_counts["MH"],
        "LH_landmark_count": class_counts["LH"],
    }
    return fixed_array, moving_array, diagnostics


def plot_landmark_qc(
    sample_id: str,
    fixed: np.ndarray,
    moving_annotation: np.ndarray,
    fixed_points: np.ndarray,
    moving_points: np.ndarray,
    diagnostics: dict[str, Any],
) -> None:
    output_dir = TEST_ROOT / "landmarks_qc"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].imshow(fixed, cmap="gray")
    axes[0].scatter(fixed_points[:, 0], fixed_points[:, 1], c=np.arange(len(fixed_points)), cmap="turbo", s=24)
    axes[0].set_title("Fixed cells + landmarks")
    axes[1].imshow(np.isin(moving_annotation, sum((TARGET_IDS[x] for x in ("DG", "MH", "LH")), [])), cmap="gray")
    axes[1].scatter(moving_points[:, 0], moving_points[:, 1], c=np.arange(len(moving_points)), cmap="turbo", s=24)
    axes[1].set_title("Atlas targets + corresponding landmarks")
    for axis in axes:
        axis.axis("off")
    fig.suptitle(f"{sample_id}: {diagnostics['landmark_classes']} ({diagnostics['landmark_count']} points)")
    fig.tight_layout()
    fig.savefig(output_dir / f"{sample_id}_landmark_qc.png", dpi=220)
    plt.close(fig)


def displacement_metrics(baseline: pd.DataFrame, candidate: pd.DataFrame) -> dict[str, float]:
    joined = baseline[["ccfy", "ccfz"]].join(
        candidate[["ccfy", "ccfz"]], how="inner", lsuffix="_base", rsuffix="_candidate"
    ).dropna()
    delta = (
        joined[["ccfy_candidate", "ccfz_candidate"]].to_numpy()
        - joined[["ccfy_base", "ccfz_base"]].to_numpy()
    ) / 25
    displacement = np.sqrt(np.sum(delta**2, axis=1))
    base_points = joined[["ccfy_base", "ccfz_base"]].to_numpy() / 25
    anchor = (base_points[:, 1] < np.quantile(base_points[:, 1], 0.30)) | (
        base_points[:, 0] < np.quantile(base_points[:, 0], 0.20)
    )
    return {
        "median_displacement_25um": float(np.median(displacement)),
        "p90_displacement_25um": float(np.quantile(displacement, 0.90)),
        "left_top_anchor_median_displacement_25um": float(np.median(displacement[anchor])),
    }


def plot_outputs(
    sample_id: str,
    iteration: int,
    cells: pd.DataFrame,
    annotation_25: np.ndarray,
    annotation_10: np.ndarray,
    fixed: np.ndarray,
    moving_spline: np.ndarray,
    transformed_borders: np.ndarray,
    regional,
) -> None:
    output_dir = FIGURE_ROOT / f"iteration_{iteration:02d}" / sample_id
    output_dir.mkdir(parents=True, exist_ok=True)
    atlas_z = int(np.rint(np.nanmedian(cells["ccfx"]) / 25))
    atlas_slice = annotation_25[atlas_z]
    title = f"{sample_id}: iteration {iteration:02d}"

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(cells["ccfz"] / 25, cells["ccfy"] / 25, s=0.08, c=cells["subclass_color"])
    for contour in regional.make_atlas_outlines(atlas_slice):
        ax.plot(contour[:, 1], contour[:, 0], "k", linewidth=0.35)
    ax.set(xlim=(0, atlas_slice.shape[1]), ylim=(atlas_slice.shape[0], 0), title=title)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_dir / f"{sample_id} 2d reg.png", dpi=300)
    plt.close(fig)

    sampled = regional.sample_true_atlas_surface(cells, annotation_10)
    fig, ax = plt.subplots(figsize=(7, 10))
    ax.scatter(cells["ccfz"] / 10, cells["ccfy"] / 10, s=0.5, c=cells["subclass_color"])
    for contour in regional.make_atlas_outlines(sampled):
        ax.plot(contour[:, 1], contour[:, 0], "k", linewidth=1)
    ax.set(xlim=(0, 1140), ylim=(800, 0), title=title)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_dir / f"{sample_id} 3d reg.png", dpi=300)
    plt.close(fig)

    colors = {"DG": "#ff7f0e", "MH": "#1f77b4", "LH": "#d62728", "Ependymal": "#2ca02c"}
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(cells["ccfz"] / 25, cells["ccfy"] / 25, s=0.08, c="0.75")
    for label, cell_class in TARGET_CLASSES.items():
        group = cells[cells["subclass_label_transfer"].eq(cell_class)]
        ax.scatter(group["ccfz"] / 25, group["ccfy"] / 25, s=2, c=colors[label], label=f"{label} cells")
        for contour_index, contour in enumerate(
            skimage.measure.find_contours(np.isin(atlas_slice, TARGET_IDS[label]).astype(float), 0.5)
        ):
            ax.plot(contour[:, 1], contour[:, 0], c=colors[label], linewidth=1.5,
                    label=f"{label} atlas" if contour_index == 0 else None)
    ax.set(xlim=(0, atlas_slice.shape[1]), ylim=(atlas_slice.shape[0], 0), title=title)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", markerscale=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / f"{sample_id} regional targets.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    axes[0].imshow(fixed, cmap="gray")
    axes[0].set_title("Target-cell enhanced fixed image")
    axes[1].imshow(moving_spline, cmap="gray")
    axes[1].set_title("Warped Allen Nissl")
    axes[2].imshow(fixed, cmap="gray")
    axes[2].imshow(np.ma.masked_less_equal(transformed_borders, 0), cmap="magma", alpha=0.8)
    axes[2].set_title("Fixed image + warped Allen borders")
    for axis in axes:
        axis.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_dir / f"{sample_id} regional registration QC.png", dpi=250)
    plt.close(fig)


def run_candidate(
    sample_id: str,
    weight: float,
    baseline_cells: pd.DataFrame,
    eligibility_row: pd.Series,
    annotation_25: np.ndarray,
    annotation_10: np.ndarray,
    nissl: np.ndarray,
    borders: np.ndarray,
    regional,
    iteration: int | None = None,
    sheet_name: str | None = None,
    min_cells: dict[str, int] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame]:
    iteration = ITERATIONS[weight] if iteration is None else iteration
    run_dir = TEST_ROOT / f"iteration_{iteration:02d}" / sample_id
    if run_dir.exists():
        raise FileExistsError(
            f"Refusing to reuse non-unique run directory {run_dir}; remove it or use --resume"
        )
    run_dir.mkdir(parents=True)
    prior_cwd = Path.cwd()
    os.chdir(run_dir)
    try:
        if sheet_name is None:
            sheet_name = "landmarks003" if weight == 0.03 else "landmarks008"
        df = pd.read_excel(WORKBOOK, sheet_name=sheet_name).reset_index(drop=True)
        num = int(df.index[df["cell_metadata"].eq(sample_id)][0])
        row = df.iloc[num]
        sample_cells = baseline_cells[baseline_cells["sample_id"].eq(sample_id)].copy()

        fixed = regional.modify_dapi(
            df, num, baseline_cells, cell_types=parse_list(row["cell_types"]),
            factor=float(row["dapi_enhance_factor"])
        )
        moving_annotation = annotation_25[int(row["allen_slice_num"])]
        moving_borders = borders[int(row["allen_slice_num"])]
        moving = regional.modify_nissl(
            nissl[int(row["allen_slice_num"])], moving_annotation,
            factor=float(row["nissl_enhance_factor"]),
            annot_dict=parse_list(row["annots_to_amplify"]),
        )
        params_affine, params_spline = regional.params_from_df(df, num)
        pad_width = int(row["pad_width"]) if pd.notna(row["pad_width"]) else 20
        fixed_bbox, fixed, _ = regional.crop_and_pad_image(
            fixed, pad_width=pad_width, area_thresh=float(row["area_thresh"])
        )
        _, _, moving_borders = regional.crop_and_pad_image(
            moving, moving_borders, pad_width=pad_width, area_thresh=float(row["area_thresh"])
        )
        moving_bbox, moving, moving_annotation = regional.crop_and_pad_image(
            moving, moving_annotation, pad_width=pad_width, area_thresh=float(row["area_thresh"])
        )
        cells = regional.get_cell_metadata_for_slice_index(
            df, num, baseline_cells, ccf_pixel_size=25, bbox=fixed_bbox, pad_width=pad_width
        )
        eligible = {label: bool(eligibility_row[f"{label}_eligible"]) for label in ("DG", "MH", "LH")}
        fixed_landmarks, moving_landmarks, diagnostics = make_robust_landmarks(
            cells, moving_annotation, eligible, min_cells=min_cells
        )
        if not len(fixed_landmarks):
            raise RuntimeError(f"{sample_id}: no robust landmarks generated")
        regional.write_pts_file(fixed_landmarks, name="fix.pts")
        regional.write_pts_file(moving_landmarks, name="mov.pts")
        plot_landmark_qc(
            sample_id, fixed, moving_annotation, fixed_landmarks, moving_landmarks, diagnostics
        )
        for parameter_map in (params_affine, params_spline):
            parameter_map["Registration"] = regional.L2P(["MultiMetricMultiResolutionRegistration"])
            parameter_map["Metric"] = regional.L2P(
                ["AdvancedMattesMutualInformation", "CorrespondingPointsEuclideanDistanceMetric"]
            )
            parameter_map["Metric0Weight"] = regional.L2P([1 - weight])
            parameter_map["Metric1Weight"] = regional.L2P([weight])

        transform, moving_spline = regional.register_images(
            fixed, moving, params_affine, params_spline
        )
        transformed_borders = regional.transform_image(moving_borders, transform)
        regional.write_pts_file(cells[["fixed_x", "fixed_y"]].to_numpy())
        transformix = regional.sitk.TransformixImageFilter()
        transformix.SetTransformParameterMap(transform)
        transformix.SetMovingImage(regional.sitk.GetImageFromArray(moving))
        transformix.SetFixedPointSetFileName("points.pts")
        transformix.Execute()
        moving_points = regional.read_outputpoints_file()[:, 3]
        moving_points -= pad_width
        moving_points += np.array([moving_bbox[1], moving_bbox[0]])

        updated = sample_cells.copy()
        for column in ("ccfx", "ccfy", "ccfz", "annotation"):
            updated.loc[:, column] = np.nan
        z = np.full(len(moving_points), int(row["allen_slice_num"]), dtype=float)
        y, x = moving_points[:, 1], moving_points[:, 0]
        zi, yi, xi = z.astype(int), np.rint(y).astype(int), np.rint(x).astype(int)
        valid = (
            (yi >= 0) & (yi < annotation_25.shape[1])
            & (xi >= 0) & (xi < annotation_25.shape[2])
        )
        indices = cells.index[valid]
        updated.loc[indices, "ccfx"] = 25 * z[valid]
        updated.loc[indices, "ccfy"] = 25 * y[valid]
        updated.loc[indices, "ccfz"] = 25 * x[valid]
        updated.loc[indices, "annotation"] = annotation_25[zi[valid], yi[valid], xi[valid]]
        updated.to_csv(run_dir / f"{sample_id}_ccf2d.csv", index=True)

        metrics = safe_target_metrics(updated, annotation_25)
        metrics.update(displacement_metrics(sample_cells, updated))
        metrics.update(
            {
                "sample_id": sample_id,
                "iteration": iteration,
                "weight": weight,
                **diagnostics,
            }
        )
        plot_outputs(
            sample_id, iteration, updated, annotation_25, annotation_10,
            fixed, moving_spline, transformed_borders, regional
        )
        (run_dir / "run_metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=True))
        return metrics, updated
    finally:
        os.chdir(prior_cwd)


def evaluate_candidate(
    baseline: pd.Series, candidate: pd.Series, eligible: pd.Series
) -> tuple[str, str, bool]:
    """Apply the plan's hard-reject and target-specific acceptance gates."""
    hard_rejects: list[str] = []
    anchor = float(candidate["left_top_anchor_median_displacement_25um"])
    if anchor > 8.0:
        hard_rejects.append(f"anchor {anchor:.2f}>8")
    if (
        pd.notna(baseline["all_tissue_containment_pct"])
        and pd.notna(candidate["all_tissue_containment_pct"])
        and candidate["all_tissue_containment_pct"] < baseline["all_tissue_containment_pct"] - 1.0
    ):
        hard_rejects.append("tissue containment dropped >1 pp")
    for label in ("MH", "LH"):
        if bool(eligible[f"{label}_eligible"]):
            b, c = baseline[f"{label}_inside_pct"], candidate[f"{label}_inside_pct"]
            if pd.notna(b) and pd.notna(c) and c < b - 10:
                hard_rejects.append(f"{label} inside dropped >10 pp")
    if hard_rejects:
        decision = "reject_overwarp" if anchor > 8 else "reject_hard_gate"
        return decision, "; ".join(hard_rejects), False

    target_gates: dict[str, bool] = {}
    improvements: list[float] = []
    if bool(eligible["DG_eligible"]) and "DG" in str(candidate["landmark_classes"]):
        b_p90, c_p90 = baseline["DG_p90_distance_25um"], candidate["DG_p90_distance_25um"]
        b_in, c_in = baseline["DG_inside_pct"], candidate["DG_inside_pct"]
        target_gates["DG"] = bool(
            c_p90 <= max(3.0, 0.35 * b_p90)
            and (c_in >= b_in + 8.0 or c_in >= 68.0)
        )
        improvements.append(float((b_p90 - c_p90) + max(0, c_in - b_in) / 10))
    for label, delta in (("LH", 30.0), ("MH", 20.0)):
        if bool(eligible[f"{label}_eligible"]) and label in str(candidate["landmark_classes"]):
            b_in, c_in = baseline[f"{label}_inside_pct"], candidate[f"{label}_inside_pct"]
            target_gates[label] = bool(c_in >= 60.0 or c_in >= b_in + delta)
            improvements.append(float(c_in - b_in))

    if anchor <= 6.5 and target_gates and all(target_gates.values()):
        return "accept", "all applicable regional and anchor gates passed", False
    under_corrected = (
        anchor <= 6.5
        and bool(target_gates)
        and not all(target_gates.values())
        and any(value > 0 for value in improvements)
    )
    failed = ",".join(label for label, passed in target_gates.items() if not passed)
    return "under_corrected" if under_corrected else "reject_no_benefit", f"failed target gates: {failed}", under_corrected


def append_workbook_results(
    metrics: pd.DataFrame, decisions: pd.DataFrame
) -> None:
    with pd.ExcelWriter(WORKBOOK, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        metrics.to_excel(writer, sheet_name="candidate_metrics", index=False)
        decisions.to_excel(writer, sheet_name="decisions", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--triage-only", action="store_true", help="Build workbook, metrics, eligibility, and QC inputs without Elastix")
    parser.add_argument("--samples", nargs="+", help="Optional subset of baseline sample IDs")
    parser.add_argument("--resume", action="store_true", help="Reuse completed candidate CSV/JSON products instead of rerunning")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    regional = load_regional_module()
    annotation_25, _ = nrrd.read(str(REFERENCE_ROOT / "annotation_25.nrrd"), index_order="F")
    annotation_25 = annotation_25[:, :, :236]

    source_params = pd.read_excel(SOURCE_WORKBOOK, sheet_name="slice_parameters")
    baseline_sample_ids = set(
        pd.read_csv(BASELINE_CCF, usecols=["sample_id"])["sample_id"].dropna().unique()
    )
    samples = [
        sample for sample in source_params["cell_metadata"].tolist()
        if sample in baseline_sample_ids
    ]
    if args.samples:
        requested = set(args.samples)
        unknown = sorted(requested - set(samples))
        if unknown:
            raise ValueError(f"Unknown requested samples: {unknown}")
        samples = [sample for sample in samples if sample in requested]
    baseline_cells = load_baseline_cells(set(samples))

    baseline_records = []
    eligibility_records = []
    for sample_id in samples:
        sample_cells = baseline_cells[baseline_cells["sample_id"].eq(sample_id)]
        metrics = {"sample_id": sample_id, **safe_target_metrics(sample_cells, annotation_25)}
        baseline_records.append(metrics)
        eligibility_records.append({"sample_id": sample_id, **landmark_eligibility(metrics)})
    baseline_metrics = pd.DataFrame(baseline_records)
    eligibility = pd.DataFrame(eligibility_records)
    baseline_metrics.to_csv(BASELINE_METRICS_CSV, index=False)
    eligibility.to_csv(ELIGIBILITY_CSV, index=False)
    build_workbook(samples, baseline_metrics, eligibility)
    print(f"Parameterized {len(samples)} samples in {WORKBOOK}", flush=True)
    print(eligibility.to_string(index=False), flush=True)
    if args.triage_only:
        return

    nissl, _ = nrrd.read(str(REFERENCE_ROOT / "ara_nissl_25.nrrd"), index_order="F")
    nissl = nissl[:, :, :236]
    borders = skimage.io.imread(str(REFERENCE_ROOT / "annotations_25_border_only.tif"))[:, :, :236]
    annotation_10, _ = nrrd.read(str(REFERENCE_ROOT / "annotation_10.nrrd"), index_order="F")
    baseline_by_sample = baseline_metrics.set_index("sample_id")
    eligibility_by_sample = eligibility.set_index("sample_id")
    candidate_records: list[dict[str, Any]] = []
    decision_records: list[dict[str, Any]] = []

    for sample_id in samples:
        e = eligibility_by_sample.loc[sample_id]
        baseline = baseline_by_sample.loc[sample_id]
        if not bool(e["needs_regional_optimization"]):
            decision_records.append(
                {"sample_id": sample_id, "decision": "retain_baseline", "accepted_iteration": np.nan,
                 "reason": "baseline passes all applicable absolute gates"}
            )
            continue
        if not e["eligible_classes"]:
            decision_records.append(
                {"sample_id": sample_id, "decision": "retain_baseline", "accepted_iteration": np.nan,
                 "reason": "no eligible regional landmarks"}
            )
            continue

        accepted = False
        for weight in (0.03, 0.08):
            iteration = ITERATIONS[weight]
            run_dir = TEST_ROOT / f"iteration_{iteration:02d}" / sample_id
            if weight == 0.08 and not decision_records[-1].get("try_rescue", False):
                break
            if args.resume and (run_dir / "run_metrics.json").exists() and (run_dir / f"{sample_id}_ccf2d.csv").exists():
                metrics = json.loads((run_dir / "run_metrics.json").read_text())
                updated = pd.read_csv(run_dir / f"{sample_id}_ccf2d.csv", index_col=0)
            else:
                metrics, updated = run_candidate(
                    sample_id, weight, baseline_cells, e,
                    annotation_25, annotation_10, nissl, borders, regional
                )
            candidate_records.append(metrics)
            candidate = pd.Series(metrics)
            decision, reason, try_rescue = evaluate_candidate(baseline, candidate, e)
            decision_row = {
                "sample_id": sample_id,
                "iteration": iteration,
                "weight": weight,
                "decision": decision,
                "accepted_iteration": iteration if decision == "accept" else np.nan,
                "try_rescue": bool(weight == 0.03 and try_rescue),
                "reason": reason,
            }
            decision_records.append(decision_row)
            if decision == "accept":
                accepted_path = run_dir / f"{sample_id}_accepted_ccf2d.csv"
                shutil.copy2(run_dir / f"{sample_id}_ccf2d.csv", accepted_path)
                accepted = True
                break
            if weight == 0.03 and not try_rescue:
                break
        if not accepted and decision_records[-1]["decision"] == "under_corrected":
            decision_records[-1]["decision"] = "retain_baseline_after_under_correction"

    candidate_metrics = pd.DataFrame(candidate_records)
    decisions = pd.DataFrame(decision_records)
    candidate_metrics.to_csv(METRICS_CSV, index=False)
    decisions.to_csv(DECISIONS_CSV, index=False)
    append_workbook_results(candidate_metrics, decisions)
    print(candidate_metrics.to_string(index=False), flush=True)
    print(decisions.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
