#!/usr/bin/env python
"""Targeted roi1_run2 registration sweep for DG, ventricular, and habenular landmarks.

This script reads only the isolated experimental workbook and writes only to
Allen_CCF_regional_tests and Allen_CCF_alignment_optimized test iterations.
It never modifies the canonical or optimized registration workbooks/results.
"""

from __future__ import annotations

import ast
import json
import os
import pickle
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nrrd
import numpy as np
import pandas as pd
import scipy.ndimage
import skimage.io
import skimage.measure
import skimage.morphology
import skimage.segmentation


BASE = Path(__file__).resolve().parents[2]
NOTEBOOK = BASE / "notebooks/0X_Allen_CCF_Registration_optimized.ipynb"
TEST_ROOT = BASE / "data/interim/registration/Allen_CCF_regional_tests"
WORKBOOK = TEST_ROOT / "slice_positions_25um_roi1_run2_regional_tests.xlsx"
FIGURE_ROOT = BASE / "figures/Allen_CCF_alignment_optimized"
REFERENCE_ROOT = BASE / "data/input/allen_registration_ref"
BASELINE_CCF = BASE / "data/interim/registration/Allen_CCF_optimization/ccf3d.csv"
SAMPLE_ID = "roi1_run2"
VENTRICLE_IDS = [73, 81, 89, 98, 108, 116, 124, 129, 140, 145, 153, 164]
TARGET_CLASSES = {
    "MH": "145 MH Tac2 Glut",
    "LH": "146 LH Pou4f1 Sox1 Glut",
    "Ependymal": "323 Ependymal NN",
    "DG": "037 DG Glut",
}
DG_GRANULE_IDS = [632, 758, 790, 823]
RUNS = [
    ("loop_04_target5_grid13", 4),
    ("loop_05_target10_grid13", 5),
    ("loop_06_target8_grid9", 6),
    ("loop_07_target3_grid13", 7),
    ("loop_08_DGgranule_target3", 8),
    ("loop_09_landmarks003", 9),
    ("loop_10_landmarks008", 10),
]


def load_notebook_functions() -> None:
    """Load the workflow's original registration helpers without executing it."""
    notebook = json.loads(NOTEBOOK.read_text())
    namespace = globals()
    for cell_index in (0, 2, 4, 35):
        source = "".join(notebook["cells"][cell_index]["source"])
        source = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("%")
        )
        exec(compile(source, f"{NOTEBOOK}:cell_{cell_index}", "exec"), namespace)


def read_sample(path: Path) -> pd.DataFrame:
    """Read one sample from the large cell table without retaining other samples."""
    chunks = []
    for chunk in pd.read_csv(path, index_col=0, chunksize=125_000):
        selected = chunk[chunk["sample_id"].eq(SAMPLE_ID)]
        if len(selected):
            chunks.append(selected.copy())
    return pd.concat(chunks, axis=0)


def parse_list(value) -> list:
    if pd.isna(value):
        return []
    if isinstance(value, list):
        return value
    return list(ast.literal_eval(str(value)))


def make_atlas_outlines(annotation_slice: np.ndarray) -> list[np.ndarray]:
    outlines = []
    for value in np.unique(annotation_slice):
        if value == 0:
            continue
        outlines.extend(
            skimage.measure.find_contours((annotation_slice == value).astype(float), 0.99)
        )
    return outlines


def sample_true_atlas_surface(
    cells: pd.DataFrame, annotation_10: np.ndarray, suffix: str = ""
) -> np.ndarray:
    columns = [f"ccfx{suffix}", f"ccfy{suffix}", f"ccfz{suffix}"]
    valid = cells.dropna(subset=columns)
    image_shape = annotation_10[0].shape
    ap_surface = np.zeros(image_shape, dtype=np.uint16)
    y = np.rint(valid[columns[1]].to_numpy() / 10).astype(int)
    x = np.rint(valid[columns[2]].to_numpy() / 10).astype(int)
    z = np.rint(valid[columns[0]].to_numpy() / 10).astype(int)
    in_bounds = (
        (y >= 0)
        & (y < image_shape[0])
        & (x >= 0)
        & (x < image_shape[1])
        & (z >= 0)
        & (z < annotation_10.shape[0])
    )
    ap_surface[y[in_bounds], x[in_bounds]] = z[in_bounds]
    ap_surface = skimage.segmentation.expand_labels(ap_surface, distance=30)
    yy, xx = np.mgrid[0 : image_shape[0], 0 : image_shape[1]]
    return annotation_10[ap_surface, yy, xx]


def target_distance_metrics(
    cells: pd.DataFrame, annotation_25: np.ndarray, suffix: str = ""
) -> dict[str, float]:
    xcol, ycol, zcol = f"ccfx{suffix}", f"ccfy{suffix}", f"ccfz{suffix}"
    valid = cells.dropna(subset=[xcol, ycol, zcol])
    atlas_z = int(np.rint(np.median(valid[xcol]) / 25))
    atlas_slice = annotation_25[atlas_z]
    result: dict[str, float] = {"atlas_slice": atlas_z}

    points_all = np.rint(valid[[ycol, zcol]].to_numpy() / 25).astype(int)
    in_bounds = (
        (points_all[:, 0] >= 0)
        & (points_all[:, 0] < atlas_slice.shape[0])
        & (points_all[:, 1] >= 0)
        & (points_all[:, 1] < atlas_slice.shape[1])
    )
    points_all = points_all[in_bounds]
    result["all_tissue_containment_pct"] = float(
        100 * np.mean(atlas_slice[points_all[:, 0], points_all[:, 1]] > 0)
    )

    region_ids = {
        "MH": [483],
        "LH": [186],
        "Ependymal": VENTRICLE_IDS,
        "DG": DG_GRANULE_IDS,
    }
    for label, cell_class in TARGET_CLASSES.items():
        group = valid[valid["subclass_label_transfer"].eq(cell_class)]
        points = np.rint(group[[ycol, zcol]].to_numpy() / 25).astype(int)
        point_valid = (
            (points[:, 0] >= 0)
            & (points[:, 0] < atlas_slice.shape[0])
            & (points[:, 1] >= 0)
            & (points[:, 1] < atlas_slice.shape[1])
        )
        points = points[point_valid]
        mask = np.isin(atlas_slice, region_ids[label])
        distance = scipy.ndimage.distance_transform_edt(~mask)
        distances = distance[points[:, 0], points[:, 1]]
        result[f"{label}_n"] = int(len(points))
        result[f"{label}_median_distance_25um"] = float(np.median(distances))
        result[f"{label}_p90_distance_25um"] = float(np.quantile(distances, 0.9))
        result[f"{label}_inside_pct"] = float(100 * np.mean(distances == 0))
    return result


def displacement_metrics(
    baseline: pd.DataFrame, candidate: pd.DataFrame
) -> dict[str, float]:
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
        "p90_displacement_25um": float(np.quantile(displacement, 0.9)),
        "left_top_anchor_median_displacement_25um": float(
            np.median(displacement[anchor])
        ),
    }


def make_corresponding_landmarks(
    fixed_cells: pd.DataFrame, moving_annotation: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Pair DG blades plus MH/LH centroids in cropped registration space."""
    fixed_points: list[list[float]] = []
    moving_points: list[list[float]] = []

    dg_fixed = fixed_cells[
        fixed_cells["subclass_label_transfer"].eq(TARGET_CLASSES["DG"])
    ][["fixed_x", "fixed_y"]].dropna().to_numpy()
    dg_y, dg_x = np.where(np.isin(moving_annotation, DG_GRANULE_IDS))
    dg_moving = np.column_stack([dg_x, dg_y])
    for fixed_cloud, moving_cloud in ((dg_fixed, dg_moving),):
        fixed_edges = np.quantile(fixed_cloud[:, 0], np.linspace(0.08, 0.92, 7))
        moving_edges = np.quantile(moving_cloud[:, 0], np.linspace(0.08, 0.92, 7))
        for bin_index in range(len(fixed_edges) - 1):
            fixed_bin = fixed_cloud[
                (fixed_cloud[:, 0] >= fixed_edges[bin_index])
                & (fixed_cloud[:, 0] <= fixed_edges[bin_index + 1])
            ]
            moving_bin = moving_cloud[
                (moving_cloud[:, 0] >= moving_edges[bin_index])
                & (moving_cloud[:, 0] <= moving_edges[bin_index + 1])
            ]
            if len(fixed_bin) < 10 or len(moving_bin) < 10:
                continue
            for y_quantile in (0.2, 0.8):
                fixed_points.append(
                    [
                        float(np.median(fixed_bin[:, 0])),
                        float(np.quantile(fixed_bin[:, 1], y_quantile)),
                    ]
                )
                moving_points.append(
                    [
                        float(np.median(moving_bin[:, 0])),
                        float(np.quantile(moving_bin[:, 1], y_quantile)),
                    ]
                )

    for label, annotation_ids in (("MH", [483]), ("LH", [186])):
        group = fixed_cells[
            fixed_cells["subclass_label_transfer"].eq(TARGET_CLASSES[label])
        ][["fixed_x", "fixed_y"]].dropna()
        moving_y, moving_x = np.where(np.isin(moving_annotation, annotation_ids))
        if len(group) and len(moving_x):
            fixed_points.append(group.median().to_list())
            moving_points.append(
                [float(np.median(moving_x)), float(np.median(moving_y))]
            )

    return np.asarray(fixed_points), np.asarray(moving_points)


def plot_registration(
    output_dir: Path,
    cells: pd.DataFrame,
    annotation_10: np.ndarray,
    title: str,
) -> None:
    sampled = sample_true_atlas_surface(cells, annotation_10)
    outlines = make_atlas_outlines(sampled)
    fig = plt.figure(figsize=(7, 10))
    plt.scatter(
        cells["ccfz"] / 10,
        cells["ccfy"] / 10,
        s=0.5,
        c=cells["subclass_color"],
    )
    for outline in outlines:
        plt.plot(*np.flip(outline.T), "k", linewidth=1)
    ax = plt.gca()
    ax.set_title(title)
    ax.set_ylim(0, 800)
    ax.set_xlim(0, 1140)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    plt.axis("off")
    plt.tight_layout()
    fig.savefig(output_dir / f"{SAMPLE_ID} 3d reg.png", dpi=300)
    plt.close(fig)


def plot_targets(
    output_dir: Path,
    cells: pd.DataFrame,
    annotation_25: np.ndarray,
    title: str,
) -> None:
    atlas_z = int(np.rint(np.nanmedian(cells["ccfx"]) / 25))
    atlas_slice = annotation_25[atlas_z]
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(cells["ccfz"] / 25, cells["ccfy"] / 25, s=0.08, c="0.75")
    colors = {
        "MH": "#1f77b4",
        "LH": "#d62728",
        "Ependymal": "#2ca02c",
        "DG": "#ff7f0e",
    }
    ids = {
        "MH": [483],
        "LH": [186],
        "Ependymal": VENTRICLE_IDS,
        "DG": DG_GRANULE_IDS,
    }
    for label, cell_class in TARGET_CLASSES.items():
        group = cells[cells["subclass_label_transfer"].eq(cell_class)]
        ax.scatter(
            group["ccfz"] / 25,
            group["ccfy"] / 25,
            s=2,
            c=colors[label],
            label=f"{label} cells",
        )
        mask = np.isin(atlas_slice, ids[label])
        contours = skimage.measure.find_contours(mask.astype(float), 0.5)
        for contour_index, contour in enumerate(contours):
            ax.plot(
                contour[:, 1],
                contour[:, 0],
                color=colors[label],
                linewidth=1.5,
                label=f"{label} atlas" if contour_index == 0 else None,
            )
    ax.set_title(title)
    ax.set_xlim(0, atlas_slice.shape[1])
    ax.set_ylim(atlas_slice.shape[0], 0)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", markerscale=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / f"{SAMPLE_ID} regional targets.png", dpi=300)
    plt.close(fig)


def run_candidate(
    sheet_name: str,
    iteration: int,
    baseline_cells: pd.DataFrame,
    annotation_25: np.ndarray,
    annotation_10: np.ndarray,
    nissl: np.ndarray,
    borders: np.ndarray,
) -> dict:
    output_dir = FIGURE_ROOT / f"iteration_{iteration:02d}"
    run_dir = TEST_ROOT / f"iteration_{iteration:02d}"
    output_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(run_dir)

    df = pd.read_excel(WORKBOOK, sheet_name=sheet_name).reset_index(drop=True)
    num = int(df.index[df["cell_metadata"].eq(SAMPLE_ID)][0])
    row = df.iloc[num]
    cell_types = parse_list(row["cell_types"])
    annotations = parse_list(row["annots_to_amplify"])

    fixed = modify_dapi(
        df,
        num,
        baseline_cells,
        cell_types=cell_types,
        factor=float(row["dapi_enhance_factor"]),
    )
    moving_annotation = annotation_25[int(row["allen_slice_num"])]
    moving_borders = borders[int(row["allen_slice_num"])]
    moving = modify_nissl(
        nissl[int(row["allen_slice_num"])],
        moving_annotation,
        factor=float(row["nissl_enhance_factor"]),
        annot_dict=annotations,
    )
    params_affine, params_spline = params_from_df(df, num)
    pad_width = 20
    fixed_bbox, fixed, _ = crop_and_pad_image(
        fixed, pad_width=pad_width, area_thresh=float(row["area_thresh"])
    )
    _, _, moving_borders = crop_and_pad_image(
        moving,
        moving_borders,
        pad_width=pad_width,
        area_thresh=float(row["area_thresh"]),
    )
    moving_bbox, moving, moving_annotation = crop_and_pad_image(
        moving,
        moving_annotation,
        pad_width=pad_width,
        area_thresh=float(row["area_thresh"]),
    )

    cells = get_cell_metadata_for_slice_index(
        df,
        num,
        baseline_cells,
        ccf_pixel_size=25,
        bbox=fixed_bbox,
    )
    landmark_weight = (
        float(row["landmark_weight"])
        if "landmark_weight" in row.index and not pd.isna(row["landmark_weight"])
        else 0.0
    )
    landmark_count = 0
    if landmark_weight > 0:
        fixed_landmarks, moving_landmarks = make_corresponding_landmarks(
            cells, moving_annotation
        )
        landmark_count = len(fixed_landmarks)
        write_pts_file(fixed_landmarks, name="fix.pts")
        write_pts_file(moving_landmarks, name="mov.pts")
        for parameter_map in (params_affine, params_spline):
            parameter_map["Registration"] = L2P(
                ["MultiMetricMultiResolutionRegistration"]
            )
            parameter_map["Metric"] = L2P(
                [
                    "AdvancedMattesMutualInformation",
                    "CorrespondingPointsEuclideanDistanceMetric",
                ]
            )
            parameter_map["Metric0Weight"] = L2P([1 - landmark_weight])
            parameter_map["Metric1Weight"] = L2P([landmark_weight])

    transform, moving_spline = register_images(
        fixed, moving, params_affine, params_spline
    )
    transformed_borders = transform_image(moving_borders, transform)

    fixed_points = cells[["fixed_x", "fixed_y"]].to_numpy()
    write_pts_file(fixed_points)
    transformix = sitk.TransformixImageFilter()
    transformix.SetTransformParameterMap(transform)
    transformix.SetMovingImage(sitk.GetImageFromArray(moving))
    transformix.SetFixedPointSetFileName("points.pts")
    transformix.Execute()
    output_points = read_outputpoints_file()
    moving_points = output_points[:, 3]
    moving_points -= pad_width
    moving_points += np.array([moving_bbox[1], moving_bbox[0]])

    updated = baseline_cells.copy()
    updated.loc[:, "ccfx"] = np.nan
    updated.loc[:, "ccfy"] = np.nan
    updated.loc[:, "ccfz"] = np.nan
    updated.loc[:, "annotation"] = np.nan
    z = np.full(len(moving_points), int(row["allen_slice_num"]), dtype=float)
    y = moving_points[:, 1]
    x = moving_points[:, 0]
    zi, yi, xi = z.astype(int), y.astype(int), x.astype(int)
    valid = (
        (zi >= 0)
        & (zi < annotation_25.shape[0])
        & (yi >= 0)
        & (yi < annotation_25.shape[1])
        & (xi >= 0)
        & (xi < annotation_25.shape[2])
    )
    indices = cells.index[valid]
    updated.loc[indices, "ccfx"] = 25 * z[valid]
    updated.loc[indices, "ccfy"] = 25 * y[valid]
    updated.loc[indices, "ccfz"] = 25 * x[valid]
    updated.loc[indices, "annotation"] = annotation_25[
        zi[valid], yi[valid], xi[valid]
    ]
    updated.to_csv(run_dir / f"{SAMPLE_ID}_ccf2d.csv", index=True)

    metrics = target_distance_metrics(updated, annotation_25)
    metrics.update(displacement_metrics(baseline_cells, updated))
    metrics.update(
        {
            "iteration": iteration,
            "sheet": sheet_name,
            "dapi_enhance_factor": float(row["dapi_enhance_factor"]),
            "nissl_enhance_factor": float(row["nissl_enhance_factor"]),
            "spline_grid_size": float(row["spline_grid_size"]),
            "landmark_weight": landmark_weight,
            "landmark_count": landmark_count,
        }
    )
    title = (
        f"{SAMPLE_ID}: iteration {iteration:02d}; target factor "
        f"{row['dapi_enhance_factor']}; spline grid {row['spline_grid_size']}"
    )
    plot_registration(output_dir, updated, annotation_10, title)
    plot_targets(output_dir, updated, annotation_25, title)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    axes[0].imshow(fixed, cmap="gray")
    axes[0].set_title("Target-cell enhanced fixed image")
    axes[1].imshow(moving_spline, cmap="gray")
    axes[1].set_title("Warped Allen Nissl")
    axes[2].imshow(fixed, cmap="gray")
    axes[2].imshow(
        np.ma.masked_less_equal(transformed_borders, 0),
        cmap="magma",
        alpha=0.8,
    )
    axes[2].set_title("Fixed image + warped Allen borders")
    for axis in axes:
        axis.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_dir / f"{SAMPLE_ID} regional registration QC.png", dpi=250)
    plt.close(fig)
    return metrics


def main() -> None:
    load_notebook_functions()
    global allen_name_to_annots
    with open(REFERENCE_ROOT / "allen_name_to_annots.pkl", "rb") as handle:
        allen_name_to_annots = pickle.load(handle)
    TEST_ROOT.mkdir(parents=True, exist_ok=True)
    baseline_cells = read_sample(BASELINE_CCF)
    annotation_25, _ = nrrd.read(
        str(REFERENCE_ROOT / "annotation_25.nrrd"), index_order="F"
    )
    annotation_25 = annotation_25[:, :, :236]
    nissl, _ = nrrd.read(str(REFERENCE_ROOT / "ara_nissl_25.nrrd"), index_order="F")
    nissl = nissl[:, :, :236]
    borders = skimage.io.imread(
        str(REFERENCE_ROOT / "annotations_25_border_only.tif")
    )[:, :, :236]
    annotation_10, _ = nrrd.read(
        str(REFERENCE_ROOT / "annotation_10.nrrd"), index_order="F"
    )

    baseline_metrics = target_distance_metrics(baseline_cells, annotation_25)
    baseline_metrics.update(
        {
            "iteration": 3,
            "sheet": "existing_optimized_baseline",
            "dapi_enhance_factor": np.nan,
            "nissl_enhance_factor": np.nan,
            "spline_grid_size": 13,
            "median_displacement_25um": 0.0,
            "p90_displacement_25um": 0.0,
            "left_top_anchor_median_displacement_25um": 0.0,
        }
    )
    results = [baseline_metrics]
    for sheet_name, iteration in RUNS:
        print(f"Running {sheet_name} -> iteration_{iteration:02d}", flush=True)
        results.append(
            run_candidate(
                sheet_name,
                iteration,
                baseline_cells,
                annotation_25,
                annotation_10,
                nissl,
                borders,
            )
        )
    metrics = pd.DataFrame(results).sort_values("iteration")
    metrics.to_csv(TEST_ROOT / "roi1_run2_regional_sweep_metrics.csv", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
