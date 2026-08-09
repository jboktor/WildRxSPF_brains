"""Robust multi-region corresponding landmarks.

Every point pair is built from a density-filtered cell cloud on the fixed side
and the matching Allen structure on the moving side. Regions contribute a
bounded number of points so that no single structure dominates the deformation.
"""

from __future__ import annotations

import numpy as np
import scipy.ndimage
import skimage.measure

from . import atlas_io, cell_io, config

MAX_POINTS_PER_REGION = 12


def _robust_centre(points: np.ndarray) -> list:
    return [float(np.median(points[:, 0])), float(np.median(points[:, 1]))]


def _blade_edges(cloud: np.ndarray, x: float, y_quantiles=(0.15, 0.85)) -> list:
    """Upper and lower blade points at one mediolateral station."""
    return [
        [float(x), float(np.quantile(cloud[:, 1], q))]
        for q in y_quantiles
    ]


def _medial_tip(fixed_cloud: np.ndarray, moving_cloud: np.ndarray, fraction: float = 0.18) -> tuple:
    """Anchor the medial crest and the two blade tips of the dentate.

    The standard ladder stops short of the extreme medial quantiles, so the
    crest where the blades meet was free to drift even when the body of the
    granule ribbon was well constrained. Tip stations are paired in a
    shape-normalised frame so a long thin blade cannot drag the habenula
    when absolute mediolateral extents differ.
    """
    fixed_points: list = []
    moving_points: list = []
    if len(fixed_cloud) < 30 or len(moving_cloud) < 30:
        return fixed_points, moving_points
    fixed_cut = float(np.quantile(fixed_cloud[:, 0], 1.0 - fraction))
    moving_cut = float(np.quantile(moving_cloud[:, 0], 1.0 - fraction))
    fixed_tip = fixed_cloud[fixed_cloud[:, 0] >= fixed_cut]
    moving_tip = moving_cloud[moving_cloud[:, 0] >= moving_cut]
    if len(fixed_tip) < 12 or len(moving_tip) < 12:
        return fixed_points, moving_points

    def stations(cloud: np.ndarray) -> np.ndarray:
        centre = np.median(cloud, axis=0)
        scale = np.maximum(np.std(cloud, axis=0), 1e-6)
        normalised = (cloud - centre) / scale
        # Crest = most medial; blades = high/low DV at the medial half of the tip.
        crest_idx = int(np.argmax(normalised[:, 0]))
        medial_half = cloud[normalised[:, 0] >= np.median(normalised[:, 0])]
        if len(medial_half) < 8:
            medial_half = cloud
        return np.array(
            [
                cloud[crest_idx],
                [np.median(medial_half[:, 0]), np.quantile(medial_half[:, 1], 0.2)],
                [np.median(medial_half[:, 0]), np.quantile(medial_half[:, 1], 0.8)],
            ],
            dtype=float,
        )

    fixed_stations = stations(fixed_tip)
    moving_stations = stations(moving_tip)
    # Emit crest twice so the tip has more pull than a single body-ladder station.
    for point, pair in list(zip(fixed_stations, moving_stations)) + [
        (fixed_stations[0], moving_stations[0])
    ]:
        fixed_points.append([float(point[0]), float(point[1])])
        moving_points.append([float(pair[0]), float(pair[1])])
    return fixed_points, moving_points


def _ladder(fixed_cloud: np.ndarray, moving_cloud: np.ndarray, bins: int = 6, with_tip: bool = True) -> tuple:
    """Paired upper/lower edges along the mediolateral axis (dentate blades).

    Quantile coverage runs near the medial and lateral extremes. When
    `with_tip` is set, dedicated crest anchors are appended so the medial tip
    cannot free-float. Tip anchors are skipped on rescue slices whose
    geometry-only fit is still broken, because they over-pull the habenula.
    """
    fixed_points: list = []
    moving_points: list = []
    bins = int(np.clip(min(len(fixed_cloud), len(moving_cloud)) // 40, 3, bins))
    edges = np.linspace(0.06, 0.94, bins + 1)
    fixed_edges = np.quantile(fixed_cloud[:, 0], edges)
    moving_edges = np.quantile(moving_cloud[:, 0], edges)
    for index in range(bins):
        fixed_bin = fixed_cloud[
            (fixed_cloud[:, 0] >= fixed_edges[index]) & (fixed_cloud[:, 0] <= fixed_edges[index + 1])
        ]
        moving_bin = moving_cloud[
            (moving_cloud[:, 0] >= moving_edges[index]) & (moving_cloud[:, 0] <= moving_edges[index + 1])
        ]
        if len(fixed_bin) < 10 or len(moving_bin) < 10:
            continue
        for quantile in (0.2, 0.8):
            fixed_points.append(
                [float(np.median(fixed_bin[:, 0])), float(np.quantile(fixed_bin[:, 1], quantile))]
            )
            moving_points.append(
                [float(np.median(moving_bin[:, 0])), float(np.quantile(moving_bin[:, 1], quantile))]
            )
    if with_tip:
        tip_fixed, tip_moving = _medial_tip(fixed_cloud, moving_cloud)
        fixed_points.extend(tip_fixed)
        moving_points.extend(tip_moving)
    return fixed_points, moving_points


def _ribbon(fixed_cloud: np.ndarray, moving_cloud: np.ndarray, bins: int = 6) -> tuple:
    """Paired medians of a thin cell ribbon along the mediolateral axis."""
    fixed_points: list = []
    moving_points: list = []
    bins = int(np.clip(min(len(fixed_cloud), len(moving_cloud)) // 30, 3, bins))
    edges = np.linspace(0.10, 0.90, bins + 1)
    fixed_edges = np.quantile(fixed_cloud[:, 0], edges)
    moving_edges = np.quantile(moving_cloud[:, 0], edges)
    for index in range(bins):
        fixed_bin = fixed_cloud[
            (fixed_cloud[:, 0] >= fixed_edges[index]) & (fixed_cloud[:, 0] <= fixed_edges[index + 1])
        ]
        moving_bin = moving_cloud[
            (moving_cloud[:, 0] >= moving_edges[index]) & (moving_cloud[:, 0] <= moving_edges[index + 1])
        ]
        if len(fixed_bin) < 8 or len(moving_bin) < 8:
            continue
        fixed_points.append(_robust_centre(fixed_bin))
        moving_points.append(_robust_centre(moving_bin))
    return fixed_points, moving_points


def _lattice(fixed_cloud: np.ndarray, moving_cloud: np.ndarray, nx: int = 3, ny: int = 3) -> tuple:
    """Dorsal/ventral edge ladder plus a few interior stations for caudoputamen.

    Four ML medians left the CP boundary free to drift. This keeps the pull
    local to the striatum without a dense 2D grid that overpowers neighbouring
    protected structures on hard slices.
    """
    fixed_points: list = []
    moving_points: list = []
    if len(fixed_cloud) < 40 or len(moving_cloud) < 40:
        return fixed_points, moving_points
    # Boundary ladder: upper / mid / lower CP stations along ML.
    bins = int(np.clip(nx + 1, 3, 5))
    edges = np.linspace(0.12, 0.88, bins + 1)
    fixed_x = np.quantile(fixed_cloud[:, 0], edges)
    moving_x = np.quantile(moving_cloud[:, 0], edges)
    for index in range(bins):
        fixed_bin = fixed_cloud[
            (fixed_cloud[:, 0] >= fixed_x[index]) & (fixed_cloud[:, 0] <= fixed_x[index + 1])
        ]
        moving_bin = moving_cloud[
            (moving_cloud[:, 0] >= moving_x[index]) & (moving_cloud[:, 0] <= moving_x[index + 1])
        ]
        if len(fixed_bin) < 12 or len(moving_bin) < 12:
            continue
        for q in (0.12, 0.88)[:ny] if ny <= 2 else (0.15, 0.5, 0.85)[:ny]:
            fixed_points.append(
                [float(np.median(fixed_bin[:, 0])), float(np.quantile(fixed_bin[:, 1], q))]
            )
            moving_points.append(
                [float(np.median(moving_bin[:, 0])), float(np.quantile(moving_bin[:, 1], q))]
            )
    return fixed_points, moving_points


def _normalised(points: np.ndarray) -> np.ndarray:
    centre = np.median(points, axis=0)
    scale = np.maximum(np.std(points, axis=0), 1e-6)
    return (points - centre) / scale


def _components(
    fixed_cloud: np.ndarray,
    moving_mask: np.ndarray,
    min_pixels: int,
    max_pair_distance: float,
) -> tuple:
    """Pair each atlas component with its own cell cluster.

    Aggregate centroids of ventricular labels are anatomically meaningless, so
    the lateral ventricle, third ventricle, and choroid components are matched
    separately using shape-normalised positions.
    """
    fixed_points: list = []
    moving_points: list = []
    labels, count = scipy.ndimage.label(moving_mask)
    regions = [
        region
        for region in skimage.measure.regionprops(labels)
        if region.area >= min_pixels
    ]
    if not regions or len(fixed_cloud) < 20:
        return fixed_points, moving_points

    moving_centroids = np.array([[region.centroid[1], region.centroid[0]] for region in regions])
    cluster_count = min(len(regions), 4)
    fixed_clusters = _cluster_points(fixed_cloud, cluster_count)
    if not fixed_clusters:
        return fixed_points, moving_points
    fixed_centroids = np.array([np.median(cluster, axis=0) for cluster in fixed_clusters])

    moving_normalised = _normalised(moving_centroids)
    fixed_normalised = _normalised(fixed_centroids)
    used_moving = set()
    for fixed_index, fixed_point in enumerate(fixed_normalised):
        distances = np.linalg.norm(moving_normalised - fixed_point, axis=1)
        order = np.argsort(distances)
        for moving_index in order:
            if moving_index in used_moving:
                continue
            if distances[moving_index] > max_pair_distance:
                break
            used_moving.add(int(moving_index))
            fixed_points.append([float(fixed_centroids[fixed_index][0]), float(fixed_centroids[fixed_index][1])])
            moving_points.append(
                [float(moving_centroids[moving_index][0]), float(moving_centroids[moving_index][1])]
            )
            break
    return fixed_points, moving_points


def _cluster_points(points: np.ndarray, count: int) -> list:
    """Split a cloud into `count` spatial clusters without extra dependencies."""
    if count <= 1 or len(points) < 40:
        return [points]
    centres = points[np.linspace(0, len(points) - 1, count).astype(int)]
    for _ in range(25):
        distances = np.linalg.norm(points[:, None, :] - centres[None, :, :], axis=2)
        assignment = np.argmin(distances, axis=1)
        updated = np.array(
            [
                points[assignment == index].mean(axis=0) if np.any(assignment == index) else centres[index]
                for index in range(count)
            ]
        )
        if np.allclose(updated, centres):
            break
        centres = updated
    clusters = [points[assignment == index] for index in range(count)]
    return [cluster for cluster in clusters if len(cluster) >= 15]


def region_clouds(
    cells,
    annotation_slice: np.ndarray,
    entry: dict,
    allen_name_to_annots: dict,
    image_shape: tuple,
) -> dict:
    ids = config.resolve_region_ids(entry, allen_name_to_annots)
    moving_mask = atlas_io.region_mask(annotation_slice, ids)
    raw = cell_io.class_points(cells, entry.get("cell_classes", []))
    if len(raw):
        keep = cell_io.density_inliers(raw, image_shape)
        fixed_cloud = raw[keep]
    else:
        fixed_cloud = raw
    moving_rows, moving_columns = np.where(moving_mask)
    moving_cloud = np.column_stack([moving_columns, moving_rows]).astype(float)
    return {
        "ids": ids,
        "moving_mask": moving_mask,
        "fixed_cloud": fixed_cloud,
        "moving_cloud": moving_cloud,
        "n_cells_raw": int(len(raw)),
        "n_cells_kept": int(len(fixed_cloud)),
        "atlas_pixels": int(moving_mask.sum()),
    }


def eligibility(entry: dict, clouds: dict, allen_slice_num: int) -> dict:
    """Decide whether a region can contribute landmarks, and record why not."""
    reasons = []
    if clouds["n_cells_kept"] < int(entry.get("min_cells", 40)):
        reasons.append(f"cells={clouds['n_cells_kept']}<{entry.get('min_cells', 40)}")
    if clouds["atlas_pixels"] < int(entry.get("min_atlas_pixels", 25)):
        reasons.append(f"atlas_pixels={clouds['atlas_pixels']}<{entry.get('min_atlas_pixels', 25)}")
    if entry.get("max_allen_slice") is not None and allen_slice_num > int(entry["max_allen_slice"]):
        reasons.append(f"allen_slice={allen_slice_num}>{entry['max_allen_slice']}")
    if entry.get("role") in ("mask_only", "unsafe") or entry.get("landmark") == "none":
        reasons.append(f"role={entry.get('role')}")
    return {"eligible": not reasons, "reasons": "; ".join(reasons)}


def build(
    cells,
    annotation_slice: np.ndarray,
    registry: dict,
    allen_name_to_annots: dict,
    image_shape: tuple,
    allen_slice_num: int,
    holdout: str = None,
    skip_regions=(),
    dg_tip: bool = True,
) -> dict:
    """Create all corresponding point pairs for one slice.

    `dg_tip` enables medial-crest anchors on the dentate ladder. Disable it on
    slices whose geometry-only fit is still broken; tip pull then folds the
    habenula before the body of the dentate is seated.
    """
    fixed_points: list = []
    moving_points: list = []
    report = {}
    for name, entry in registry["regions"].items():
        clouds = region_clouds(cells, annotation_slice, entry, allen_name_to_annots, image_shape)
        status = eligibility(entry, clouds, allen_slice_num)
        record = {
            "role": entry.get("role"),
            "landmark": entry.get("landmark"),
            "n_cells_raw": clouds["n_cells_raw"],
            "n_cells_kept": clouds["n_cells_kept"],
            "atlas_pixels": clouds["atlas_pixels"],
            "eligible": status["eligible"],
            "reason": status["reasons"],
            "n_points": 0,
            "used_for_fitting": False,
        }
        if not status["eligible"]:
            report[name] = record
            continue
        if name == holdout:
            record["reason"] = "held out for validation"
            report[name] = record
            continue
        if name in skip_regions:
            record["reason"] = "skipped by caller"
            report[name] = record
            continue

        strategy = entry.get("landmark")
        if strategy == "ladder":
            pairs = _ladder(
                clouds["fixed_cloud"],
                clouds["moving_cloud"],
                with_tip=bool(dg_tip and name == "DG_granule"),
            )
        elif strategy == "ribbon":
            pairs = _ribbon(clouds["fixed_cloud"], clouds["moving_cloud"])
        elif strategy == "bins":
            pairs = _ribbon(clouds["fixed_cloud"], clouds["moving_cloud"], bins=4)
        elif strategy == "lattice":
            pairs = _lattice(
                clouds["fixed_cloud"],
                clouds["moving_cloud"],
                nx=int(entry.get("lattice_nx", 4)),
                ny=int(entry.get("lattice_ny", 3)),
            )
        elif strategy == "components":
            pairs = _components(
                clouds["fixed_cloud"],
                clouds["moving_mask"],
                int(entry.get("component_min_pixels", 25)),
                float(entry.get("component_max_pair_distance_25um", 18)) / 10.0,
            )
            if not pairs[0]:
                # A fragmented structure still constrains its overall position.
                pairs = (
                    [_robust_centre(clouds["fixed_cloud"])],
                    [_robust_centre(clouds["moving_cloud"])],
                )
        elif strategy == "centroid":
            pairs = (
                [_robust_centre(clouds["fixed_cloud"])],
                [_robust_centre(clouds["moving_cloud"])],
            )
        else:
            pairs = ([], [])

        region_fixed, region_moving = pairs
        max_points = int(entry.get("max_points", MAX_POINTS_PER_REGION))
        if len(region_fixed) > max_points:
            keep = np.linspace(0, len(region_fixed) - 1, max_points).astype(int)
            region_fixed = [region_fixed[index] for index in keep]
            region_moving = [region_moving[index] for index in keep]
        fixed_points.extend(region_fixed)
        moving_points.extend(region_moving)
        record["n_points"] = len(region_fixed)
        record["used_for_fitting"] = len(region_fixed) > 0
        if name == "DG_granule":
            record["dg_tip"] = bool(dg_tip)
        if not record["n_points"]:
            record["reason"] = "no stable point pairs from cloud geometry"
        report[name] = record

    return {
        "fixed": np.asarray(fixed_points, dtype=float).reshape(-1, 2),
        "moving": np.asarray(moving_points, dtype=float).reshape(-1, 2),
        "report": report,
        "dg_tip": bool(dg_tip),
    }


def choose_holdout(registry: dict, report: dict) -> str:
    """Pick a validation region that is eligible and not needed for the fit."""
    candidates = registry.get("holdout", {}).get("candidates", [])
    minimum = int(registry.get("holdout", {}).get("min_fitting_regions", 3))
    eligible = [name for name, record in report.items() if record.get("eligible") and record.get("n_points")]
    if len(eligible) <= minimum:
        return None
    for name in candidates:
        if name in eligible:
            return name
    return None
