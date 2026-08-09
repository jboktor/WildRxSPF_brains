"""Masked rigid -> affine -> B-spline registration with anti-compression control."""

from __future__ import annotations

import contextlib
import os
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from . import nb_functions


def _L2P(values):
    return tuple(str(value) for value in values)


@contextlib.contextmanager
def working_directory(path: Path):
    previous = os.getcwd()
    Path(path).mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield Path(path)
    finally:
        os.chdir(previous)


def _add_point_metric(parameter_map, landmark_weight: float, bending_weight: float) -> None:
    metrics = ["AdvancedMattesMutualInformation"]
    weights = [1.0 - landmark_weight]
    if bending_weight > 0:
        metrics.append("TransformBendingEnergyPenalty")
        weights.append(bending_weight)
    if landmark_weight > 0:
        metrics.append("CorrespondingPointsEuclideanDistanceMetric")
        weights.append(landmark_weight)
    if len(metrics) == 1:
        return
    parameter_map["Registration"] = _L2P(["MultiMetricMultiResolutionRegistration"])
    parameter_map["Metric"] = _L2P(metrics)
    for index, weight in enumerate(weights):
        parameter_map[f"Metric{index}Weight"] = _L2P([weight])


def rigid_map(iterations: int, histogram_bins: int = 32):
    parameter_map = sitk.GetDefaultParameterMap("rigid")
    parameter_map["NumberOfResolutions"] = _L2P([4])
    parameter_map["MaximumNumberOfIterations"] = _L2P([iterations])
    parameter_map["NumberOfHistogramBins"] = _L2P([histogram_bins])
    parameter_map["AutomaticTransformInitialization"] = _L2P(["true"])
    parameter_map["AutomaticTransformInitializationMethod"] = _L2P(["CenterOfGravity"])
    parameter_map["ErodeMask"] = _L2P(["false"])
    return parameter_map


def similarity_map(iterations: int, histogram_bins: int = 32):
    """Rotation, translation, and one uniform scale: no shear, no aspect change.

    Used to detect missing tissue. An affine fit absorbs a missing dorsal chunk
    by squeezing the atlas along that axis, which is the defect being looked
    for, so the detector has to use a transform that cannot do it.
    """
    parameter_map = sitk.GetDefaultParameterMap("rigid")
    parameter_map["Transform"] = _L2P(["SimilarityTransform"])
    parameter_map["NumberOfResolutions"] = _L2P([4])
    parameter_map["MaximumNumberOfIterations"] = _L2P([iterations])
    parameter_map["NumberOfHistogramBins"] = _L2P([histogram_bins])
    parameter_map["AutomaticTransformInitialization"] = _L2P(["true"])
    parameter_map["AutomaticTransformInitializationMethod"] = _L2P(["CenterOfGravity"])
    parameter_map["AutomaticScalesEstimation"] = _L2P(["true"])
    parameter_map["ErodeMask"] = _L2P(["false"])
    return parameter_map


def affine_map(iterations: int, histogram_bins: int = 32, landmark_weight: float = 0.0):
    parameter_map = sitk.GetDefaultParameterMap("affine")
    parameter_map["NumberOfResolutions"] = _L2P([4])
    parameter_map["MaximumNumberOfIterations"] = _L2P([iterations])
    parameter_map["NumberOfHistogramBins"] = _L2P([histogram_bins])
    parameter_map["ErodeMask"] = _L2P(["false"])
    _add_point_metric(parameter_map, landmark_weight, 0.0)
    return parameter_map


def spline_map(
    iterations: int,
    grid_spacing_voxels: float,
    histogram_bins: int = 32,
    landmark_weight: float = 0.0,
    bending_weight: float = 0.35,
):
    parameter_map = sitk.GetDefaultParameterMap("bspline")
    parameter_map["NumberOfResolutions"] = _L2P([4])
    parameter_map["GridSpacingSchedule"] = _L2P([20, 10, 5, 2])
    parameter_map["NumberOfHistogramBins"] = _L2P([histogram_bins])
    parameter_map["MaximumNumberOfIterations"] = _L2P([iterations])
    parameter_map["FinalGridSpacingInPhysicalUnits"] = _L2P([])
    parameter_map["FinalGridSpacingInVoxels"] = _L2P([grid_spacing_voxels])
    parameter_map["ErodeMask"] = _L2P(["false"])
    _add_point_metric(parameter_map, landmark_weight, bending_weight)
    return parameter_map


def build_parameter_maps(
    stage: str,
    iterations_affine: int,
    iterations_spline: int,
    grid_spacing_voxels: float,
    histogram_bins: int,
    landmark_weight: float,
    bending_weight: float,
    use_rigid: bool = True,
):
    """`stage` is 'similarity' for damage detection, 'affine' for search and
    priming, or 'full' for the final fit."""
    maps = [rigid_map(min(iterations_affine, 500), histogram_bins)] if use_rigid else []
    if stage == "similarity":
        maps.append(similarity_map(iterations_affine, histogram_bins))
        return maps
    if not maps:
        # Elastix needs an initialised starting point that the rigid stage would
        # otherwise provide.
        affine_first = affine_map(iterations_affine, histogram_bins, landmark_weight)
        affine_first["AutomaticTransformInitialization"] = _L2P(["true"])
        affine_first["AutomaticTransformInitializationMethod"] = _L2P(["CenterOfGravity"])
        maps.append(affine_first)
        if stage != "full":
            return maps
        maps.append(
            spline_map(
                iterations_spline,
                grid_spacing_voxels,
                histogram_bins,
                landmark_weight,
                bending_weight,
            )
        )
        return maps
    maps.append(affine_map(iterations_affine, histogram_bins, landmark_weight))
    if stage == "full":
        maps.append(
            spline_map(
                iterations_spline,
                grid_spacing_voxels,
                histogram_bins,
                landmark_weight,
                bending_weight,
            )
        )
    return maps


def register(
    fixed: np.ndarray,
    moving: np.ndarray,
    parameter_maps,
    work_dir: Path,
    fixed_mask: np.ndarray = None,
    moving_mask: np.ndarray = None,
    fixed_points: np.ndarray = None,
    moving_points: np.ndarray = None,
    log_to_file: bool = True,
):
    """Run elastix and return (transform parameter maps, resampled moving image)."""
    with working_directory(work_dir):
        if fixed_points is not None and len(fixed_points):
            nb_functions.write_pts_file(np.asarray(fixed_points), name="fix.pts")
            nb_functions.write_pts_file(np.asarray(moving_points), name="mov.pts")

        elastix = sitk.ElastixImageFilter()
        elastix.LogToConsoleOff()
        if log_to_file:
            elastix.LogToFileOn()
            elastix.SetOutputDirectory(str(work_dir))
        elastix.SetParameterMap(parameter_maps[0])
        for parameter_map in parameter_maps[1:]:
            elastix.AddParameterMap(parameter_map)
        elastix.SetFixedImage(sitk.GetImageFromArray(np.asarray(fixed, dtype=np.float32)))
        elastix.SetMovingImage(sitk.GetImageFromArray(np.asarray(moving, dtype=np.float32)))
        if fixed_points is not None and len(fixed_points):
            elastix.SetFixedPointSetFileName("fix.pts")
            elastix.SetMovingPointSetFileName("mov.pts")
        if fixed_mask is not None:
            elastix.SetFixedMask(sitk.Cast(sitk.GetImageFromArray(fixed_mask.astype(np.uint8)), sitk.sitkUInt8))
        if moving_mask is not None:
            elastix.SetMovingMask(sitk.Cast(sitk.GetImageFromArray(moving_mask.astype(np.uint8)), sitk.sitkUInt8))
        elastix.Execute()
        transform = _with_rotation_centre(elastix.GetTransformParameterMap(), fixed.shape)
        result = sitk.GetArrayFromImage(elastix.GetResultImage())
    return transform, result


def _with_rotation_centre(transform_maps, fixed_shape):
    """Restore the rotation centre elastix omits for SimilarityTransform.

    Without it transformix refuses the parameter file, so points and label maps
    cannot be pushed through a similarity stage.
    """
    centre = _L2P([(size - 1) / 2 for size in reversed(fixed_shape)])
    for parameter_map in transform_maps:
        if "CenterOfRotationPoint" not in parameter_map:
            parameter_map["CenterOfRotationPoint"] = centre
    return transform_maps


def transform_points(points_xy: np.ndarray, transform, moving: np.ndarray, work_dir: Path) -> np.ndarray:
    """Map fixed-image points into moving-image (atlas) coordinates."""
    points_xy = np.asarray(points_xy, dtype=float)
    if not len(points_xy):
        return points_xy
    with working_directory(work_dir):
        nb_functions.write_pts_file(points_xy, name="points.pts")
        transformix = sitk.TransformixImageFilter()
        transformix.LogToConsoleOff()
        transformix.SetTransformParameterMap(transform)
        transformix.SetMovingImage(sitk.GetImageFromArray(np.asarray(moving, dtype=np.float32)))
        transformix.SetFixedPointSetFileName("points.pts")
        transformix.Execute()
        output = nb_functions.read_outputpoints_file()
    return output[:, 3]


def transform_image(image: np.ndarray, transform, nearest: bool = False) -> np.ndarray:
    maps = transform
    if nearest:
        maps = list(transform)
        for parameter_map in maps:
            parameter_map["ResampleInterpolator"] = _L2P(["FinalNearestNeighborInterpolator"])
    transformix = sitk.TransformixImageFilter()
    transformix.LogToConsoleOff()
    transformix.SetTransformParameterMap(maps)
    transformix.SetMovingImage(sitk.GetImageFromArray(np.asarray(image, dtype=np.float32)))
    transformix.Execute()
    return sitk.GetArrayFromImage(transformix.GetResultImage())


def transform_labels(labels: np.ndarray, transform) -> np.ndarray:
    """Warp an integer label image without inventing intermediate labels."""
    values, inverse = np.unique(labels, return_inverse=True)
    compact = inverse.reshape(labels.shape).astype(np.float32)
    warped = transform_image(compact, transform, nearest=True)
    indices = np.clip(np.rint(warped).astype(int), 0, len(values) - 1)
    return values[indices]


def spatial_jacobian(
    transform,
    shape: tuple,
    moving: np.ndarray,
    work_dir: Path,
    step: int = 6,
) -> dict:
    """Determinant of the fixed-to-moving Jacobian on a coarse grid.

    Values above one mean atlas area is being squeezed into the fixed image,
    which is the signature of the atlas collapsing into a region where tissue is
    missing.
    """
    rows = np.arange(step, shape[0] - step, step)
    columns = np.arange(step, shape[1] - step, step)
    if len(rows) < 3 or len(columns) < 3:
        return {"determinant": np.zeros((0, 0)), "rows": rows, "columns": columns}
    grid_columns, grid_rows = np.meshgrid(columns, rows)
    points = np.column_stack([grid_columns.ravel(), grid_rows.ravel()]).astype(float)
    mapped = transform_points(points, transform, moving, work_dir)
    mapped_x = mapped[:, 0].reshape(grid_rows.shape)
    mapped_y = mapped[:, 1].reshape(grid_rows.shape)
    dxdc, dxdr = np.gradient(mapped_x, step, step, axis=(1, 0))
    dydc, dydr = np.gradient(mapped_y, step, step, axis=(1, 0))
    determinant = dxdc * dydr - dxdr * dydc
    return {
        "determinant": determinant,
        "rows": rows,
        "columns": columns,
    }


def jacobian_summary(jacobian: dict, region_mask: np.ndarray = None) -> dict:
    determinant = jacobian.get("determinant")
    if determinant is None or determinant.size == 0:
        return {"jacobian_min": float("nan"), "jacobian_p99": float("nan"), "jacobian_mean": float("nan")}
    summary = {
        "jacobian_min": float(np.nanmin(determinant)),
        "jacobian_p99": float(np.nanquantile(determinant, 0.99)),
        "jacobian_mean": float(np.nanmean(determinant)),
        "jacobian_negative_fraction": float(np.mean(determinant <= 0)),
    }
    if region_mask is not None and region_mask.any():
        rows = jacobian["rows"]
        columns = jacobian["columns"]
        sampled = region_mask[np.ix_(rows, columns)]
        if sampled.any():
            summary["jacobian_mean_in_region"] = float(np.nanmean(determinant[sampled]))
            summary["jacobian_p95_in_region"] = float(np.nanquantile(determinant[sampled], 0.95))
        else:
            summary["jacobian_mean_in_region"] = float("nan")
            summary["jacobian_p95_in_region"] = float("nan")
    return summary
