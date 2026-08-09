"""AP level search and obliquity handling with structure-specific scoring.

The historical workbook fixed one `allen_slice_num` per slice by eye. Where a
section is cut obliquely, no single level fits: medial structures prefer one
level and lateral structures another.

Levels and gradients are scored by re-reading the atlas at the candidate
geometry while holding one masked affine alignment fixed. That keeps the
comparison free of elastix run-to-run variation, which is larger than the
between-level differences being measured. Only the medial margin, which changes
the moving image itself, is scored with its own registration.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import accept, config, metrics, pipeline

ROLE_WEIGHT = {"hard": 1.0, "soft": 0.5}
SCORING_REGIONS = (
    "DG_granule",
    "CA1",
    "CA3",
    "MH",
    "LH",
    "ventricle_walls",
    "caudoputamen",
    "VMH",
)


def score_fit(
    cells_atlas: pd.DataFrame,
    annotation_plane: np.ndarray,
    registry: dict,
    allen_name_to_annots: dict,
    subset: np.ndarray = None,
) -> dict:
    """Combine tissue containment with structure-specific containment."""
    frame = cells_atlas if subset is None else cells_atlas[subset]
    containment = metrics.tissue_containment(frame, annotation_plane)
    region_results = metrics.region_metrics(frame, annotation_plane, registry, allen_name_to_annots)

    total_weight = 0.0
    total = 0.0
    detail = {}
    for name in SCORING_REGIONS:
        values = region_results.get(name)
        if values is None:
            continue
        if not values["n_cells"] or not values["atlas_pixels"] or np.isnan(values["inside_pct"]):
            continue
        weight = ROLE_WEIGHT.get(values["role"], 0.5)
        penalty = min(float(values["p90_distance_25um"]), 20.0)
        total += weight * (float(values["inside_pct"]) - 2.0 * penalty)
        total_weight += weight
        detail[name] = float(values["inside_pct"])
    region_score = total / total_weight if total_weight else float("nan")
    combined = 0.4 * containment + 0.6 * (region_score if not np.isnan(region_score) else containment)
    return {
        "score": float(combined),
        "tissue_containment_pct": float(containment),
        "region_score": float(region_score) if not np.isnan(region_score) else float("nan"),
        "region_detail": detail,
    }


def _register(
    atlas,
    frame: pd.DataFrame,
    row_index: int,
    fixed_ctx: dict,
    fixed_mask: np.ndarray,
    slice_index: int,
    gradient: float,
    medial_extra: int,
    work_dir: Path,
    cfg: config.RunConfig,
) -> tuple:
    moving_ctx = pipeline.prepare_moving(atlas, frame, row_index, slice_index, gradient, medial_extra)
    result = pipeline.fit(
        fixed_ctx,
        moving_ctx,
        work_dir=work_dir,
        stage="affine",
        fixed_mask=fixed_mask,
        moving_mask=None,
        affine_iterations=cfg.affine_iterations,
    )
    return moving_ctx, result


def search(
    atlas,
    frame: pd.DataFrame,
    row_index: int,
    fixed_ctx: dict,
    fixed_mask: np.ndarray,
    registry: dict,
    allen_name_to_annots: dict,
    work_root: Path,
    cfg: config.RunConfig,
) -> dict:
    """Choose AP level, AP gradient, and medial margin for one slice."""
    baseline_level = int(frame.iloc[row_index]["allen_slice_num"])
    base_extra = config.MIN_MEDIAL_EXTRA

    moving_ctx, priming = _register(
        atlas,
        frame,
        row_index,
        fixed_ctx,
        fixed_mask,
        baseline_level,
        0.0,
        base_extra,
        work_root / "level_probe",
        cfg,
    )
    cells_atlas = priming["cells"]
    limit = moving_ctx["plane"]["medial_limit"]

    # Split on atlas columns rather than image columns so the medial half is
    # the medial half regardless of how the section was reflected.
    columns = cells_atlas["atlas_x"].to_numpy(dtype=float)
    split = np.nanmedian(columns)
    medial = columns >= split
    lateral = columns < split

    level_records = []
    for offset in range(-8, 9):
        level = int(np.clip(baseline_level + offset, 0, atlas.annotation.shape[0] - 1))
        annotation = atlas.annotation[level, :, :limit]
        overall = score_fit(cells_atlas, annotation, registry, allen_name_to_annots)
        level_records.append(
            {
                "slice_index": level,
                "gradient": 0.0,
                "medial_extra": base_extra,
                "score": overall["score"],
                "tissue_containment_pct": overall["tissue_containment_pct"],
                "region_score": overall["region_score"],
                "medial_score": score_fit(cells_atlas, annotation, registry, allen_name_to_annots, medial)["score"],
                "lateral_score": score_fit(cells_atlas, annotation, registry, allen_name_to_annots, lateral)["score"],
                "source": "atlas lookup",
            }
        )
    best_level = max(level_records, key=lambda item: item["score"])
    medial_best = max(level_records, key=lambda item: item["medial_score"])
    lateral_best = max(level_records, key=lambda item: item["lateral_score"])
    level_disagreement = int(medial_best["slice_index"] - lateral_best["slice_index"])

    medial_centre = float(np.nanmedian(columns[medial])) if medial.any() else np.nan
    lateral_centre = float(np.nanmedian(columns[lateral])) if lateral.any() else np.nan

    gradient_records = []
    if abs(level_disagreement) >= 2 and np.isfinite(medial_centre) and np.isfinite(lateral_centre):
        span = medial_centre - lateral_centre
        if abs(span) > 20:
            estimate = level_disagreement / span
            candidates = {round(value, 3) for value in (estimate, estimate / 2)}
            candidates.update(cfg.ap_gradients)
            for gradient in sorted(candidates):
                if abs(gradient) < 0.005 or abs(gradient) > 0.08:
                    continue
                plane = atlas.oblique(best_level["slice_index"], gradient, base_extra)
                overall = score_fit(cells_atlas, plane["annotation"], registry, allen_name_to_annots)
                gradient_records.append(
                    {
                        "slice_index": best_level["slice_index"],
                        "gradient": float(gradient),
                        "medial_extra": base_extra,
                        "score": overall["score"],
                        "tissue_containment_pct": overall["tissue_containment_pct"],
                        "region_score": overall["region_score"],
                        "source": "atlas lookup",
                    }
                )

    gradient = 0.0
    if gradient_records:
        best_gradient = max(gradient_records, key=lambda item: item["score"])
        if best_gradient["score"] > best_level["score"] + cfg.ap_gradient_min_gain:
            gradient = float(best_gradient["gradient"])

    # The medial margin changes the moving image, so each candidate needs its
    # own registration.
    extra_records = []
    for extra in cfg.medial_extras:
        candidate_ctx, candidate = _register(
            atlas,
            frame,
            row_index,
            fixed_ctx,
            fixed_mask,
            best_level["slice_index"],
            gradient,
            extra,
            work_root / f"medial_{extra:02d}",
            cfg,
        )
        overall = score_fit(
            candidate["cells"], candidate_ctx["plane"]["annotation"], registry, allen_name_to_annots
        )
        pileup = metrics.boundary_pileup(candidate["cells"], candidate_ctx["plane"]["medial_limit"])
        annotation = candidate_ctx["plane"]["annotation"]
        dg_mask = annotation == 632
        if dg_mask.any():
            dg_max = int(np.where(dg_mask)[1].max())
            dg_tip_gap = int(candidate_ctx["plane"]["medial_limit"] - 1 - dg_max)
        else:
            dg_tip_gap = 99
        extra_records.append(
            {
                "slice_index": best_level["slice_index"],
                "gradient": gradient,
                "medial_extra": int(extra),
                "score": overall["score"],
                "tissue_containment_pct": overall["tissue_containment_pct"],
                "region_score": overall["region_score"],
                "boundary_pileup_pct": pileup["boundary_pileup_pct"],
                "boundary_clamp_ratio": pileup["boundary_clamp_ratio"],
                "max_atlas_x": pileup["max_atlas_x"],
                "dg_tip_gap": dg_tip_gap,
                "region_detail": overall["region_detail"],
                "source": "registration",
            }
        )
    # A margin that stacks cells against the crop edge is not a usable margin,
    # however well it scores otherwise.
    usable = [
        item for item in extra_records if item["boundary_clamp_ratio"] <= accept.BOUNDARY_CLAMP_LIMIT
    ] or extra_records
    peak = max(item["score"] for item in usable)
    # Prefer a few voxels of room past the DG crest when scores are essentially
    # tied, so the medial tip is not flush against the crop edge in the figures.
    near_peak = [item for item in usable if item["score"] >= peak - 1.5]
    best_extra = max(
        near_peak,
        key=lambda item: (min(item["dg_tip_gap"], 8), item["score"]),
    )

    return {
        "baseline_slice_index": baseline_level,
        "slice_index": int(best_level["slice_index"]),
        "gradient": float(gradient),
        "medial_extra": int(best_extra["medial_extra"]),
        "score": float(best_extra["score"]),
        "planar_best_slice_index": int(best_level["slice_index"]),
        "planar_best_score": float(best_level["score"]),
        "medial_best_slice_index": int(medial_best["slice_index"]),
        "lateral_best_slice_index": int(lateral_best["slice_index"]),
        "level_disagreement": level_disagreement,
        "obliquity_applied": bool(gradient != 0.0),
        "medial_extra_scores": {
            str(item["medial_extra"]): item["score"] for item in extra_records
        },
        "medial_extra_clamp": {
            str(item["medial_extra"]): item["boundary_clamp_ratio"] for item in extra_records
        },
        "records": level_records + gradient_records + extra_records,
    }
