"""One command that rebuilds every Allen CCF registration with full provenance.

    python -m allen_ccf_repro.cli --all

Stages per slice: masks, AP/obliquity search, masked affine priming, damage and
artifact masks, robust landmarks with a holdout region, masked B-spline fit,
acceptance gates, mask sensitivity check, and figure rendering.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "allen_ccf_repro"

from . import (  # noqa: E402
    accept,
    ap_search,
    atlas_io,
    cell_io,
    config,
    elastix_stage,
    landmarks as landmarks_module,
    masks as masks_module,
    metrics,
    nb_functions,
    pipeline,
    provenance,
    render,
)

RENDER_COLUMNS = [
    "subclass_label_transfer",
    "subclass_color",
    "fixed_x",
    "fixed_y",
    "atlas_x",
    "atlas_y",
    "ccfx",
    "ccfy",
    "ccfz",
    "annotation",
]


def baseline_reference(sample_id: str, cells_raw: pd.DataFrame, atlas, registry: dict, names: dict, frame, row_index: int):
    """Metrics for the historical fit, evaluated in its own atlas plane.

    Returns the metrics and, separately, the cells and atlas plane behind them,
    which are needed when the historical fit has to serve as a reference.
    """
    baseline_level = int(frame.iloc[row_index]["allen_slice_num"])
    plane = atlas.coronal(baseline_level, config.MIN_MEDIAL_EXTRA)
    annotation = plane["annotation"]
    frame_cells = cells_raw.dropna(subset=["ccfy", "ccfz"]).copy()
    frame_cells["atlas_x"] = frame_cells["ccfz"] / config.CCF_PIXEL_SIZE
    frame_cells["atlas_y"] = frame_cells["ccfy"] / config.CCF_PIXEL_SIZE
    record = {
        "tissue_containment_pct": metrics.tissue_containment(frame_cells, annotation),
        "boundary_pileup_pct": atlas_io.boundary_pileup_pct(
            frame_cells["atlas_x"].to_numpy(), plane["medial_limit"]
        ),
        "boundary_clamp_ratio": atlas_io.boundary_clamp_ratio(
            frame_cells["atlas_x"].to_numpy(), plane["medial_limit"]
        ),
        "allen_slice_num": baseline_level,
    }
    record.update(
        metrics.flatten_region_metrics(
            metrics.region_metrics(frame_cells, annotation, registry, names)
        )
    )
    return record, {"cells": frame_cells, "annotation": annotation}


def candidate_metrics(
    result: dict,
    moving_ctx: dict,
    cells_raw: pd.DataFrame,
    registry: dict,
    names: dict,
    holdout: str,
    exclusion: np.ndarray,
    reference_annotation: np.ndarray,
    work_dir: Path,
    displacement_reference: dict = None,
) -> dict:
    """All numbers the acceptance gates need for one candidate fit."""
    cells_atlas = result["cells"]
    annotation = moving_ctx["plane"]["annotation"]
    record = {
        "tissue_containment_pct": metrics.tissue_containment(cells_atlas, annotation),
        "n_landmarks": result["n_landmarks"],
        "landmark_weight": result["landmark_weight"],
    }
    region_results = metrics.region_metrics(cells_atlas, annotation, registry, names)
    record.update(metrics.flatten_region_metrics(region_results))
    protected_classes = [
        cell_class
        for entry in registry["regions"].values()
        if entry.get("protected")
        for cell_class in entry.get("cell_classes", [])
    ]
    record.update(
        metrics.boundary_pileup(cells_atlas, moving_ctx["plane"]["medial_limit"], protected_classes)
    )

    # Displacement only carries information when both fits share one atlas
    # plane: it then measures how far the landmark term dragged cells away from
    # the image-driven fit. Against the historical result it would mostly report
    # the deliberate change of AP level and medial margin, so that comparison is
    # recorded as a plain shift and never gated.
    updated = pipeline.update_baseline_table(cells_raw, cells_atlas)
    shift = metrics.displacement_metrics(cells_raw, updated)
    record["geometry_shift_25um"] = shift["median_displacement_25um"]
    if displacement_reference is not None:
        same_plane = metrics.displacement_metrics(displacement_reference["cells"], updated)
        record["landmark_pull_median_25um"] = same_plane["median_displacement_25um"]
        record["landmark_pull_p90_25um"] = same_plane["p90_displacement_25um"]
        record.update(
            metrics.anchor_displacement(
                displacement_reference["cells"], updated, displacement_reference["anchors"]
            )
        )

    if holdout and holdout in region_results:
        record["holdout_region"] = holdout
        record["holdout_inside_pct"] = region_results[holdout]["inside_pct"]
        record["holdout_median_distance_25um"] = region_results[holdout]["median_distance_25um"]
    else:
        record["holdout_region"] = holdout

    warped = pipeline.warped_annotation_in_fixed(moving_ctx, result["transform"])
    columns = metrics.exclusion_columns(exclusion) if exclusion.any() else None
    if columns is not None and columns.any():
        for label, region in (("cortex", "isocortex"), ("callosum", "corpus_callosum")):
            ids = config.resolve_region_ids(registry["regions"][region], names)
            new = metrics.structure_thickness(warped, ids, columns)
            old = metrics.structure_thickness(reference_annotation, ids, columns)
            ratio = float(new / old) if old and np.isfinite(old) and old > 0 else float("nan")
            record[f"{label}_thickness_ratio_near_mask"] = ratio

    # The Jacobian grid lives in fixed-image space, which is where the exclusion
    # zones are defined, so the two can be compared directly.
    jacobian = elastix_stage.spatial_jacobian(
        result["transform"], exclusion.shape, moving_ctx["image"], work_dir
    )
    record.update(elastix_stage.jacobian_summary(jacobian, exclusion if exclusion.any() else None))
    return record


def process_sample(sample_id: str, cfg: config.RunConfig) -> dict:
    """Run the whole pipeline for one slice and write its artefacts."""
    registry = config.load_registry()
    names = nb_functions.allen_name_to_annots()
    atlas = atlas_io.load_atlas()
    frame, row_index = cell_io.read_slice_parameters(sample_id)
    cells_raw = cell_io.read_sample(sample_id)

    sample_root = config.WORK_ROOT / sample_id
    sample_root.mkdir(parents=True, exist_ok=True)
    record: dict = {"sample_id": sample_id}

    fixed_ctx = pipeline.prepare_fixed(frame, row_index, cells_raw)
    overrides = masks_module.load_overrides(sample_id)
    tissue = masks_module.cell_supported_tissue_mask(fixed_ctx["image"], fixed_ctx["cells"])
    artifact = masks_module.artifact_mask(
        fixed_ctx["image"], tissue["tissue"], fixed_ctx["cells"], overrides
    )
    priming_mask = masks_module.registration_mask(
        tissue["tissue"], np.zeros(tissue["tissue"].shape, dtype=bool), artifact
    )

    search = ap_search.search(
        atlas,
        frame,
        row_index,
        fixed_ctx,
        priming_mask,
        registry,
        names,
        sample_root / "ap_search",
        cfg,
    )
    baseline, baseline_ctx = baseline_reference(
        sample_id, cells_raw, atlas, registry, names, frame, row_index
    )
    baseline_level = int(frame.iloc[row_index]["allen_slice_num"])
    fit_parameters = pipeline.slice_parameters(frame.iloc[row_index], cfg)

    def build_geometry(slice_index: int, gradient: float, medial_extra: int, tag: str) -> dict:
        """Masks, landmarks, and the geometry-only fit for one atlas plane."""
        moving_ctx = pipeline.prepare_moving(atlas, frame, row_index, slice_index, gradient, medial_extra)
        # Undistorted placement of the atlas: it overhangs wherever tissue is
        # missing instead of shrinking to fit, which is what makes the cut
        # detectable, and it is also the reference the thickness ratios need.
        priming = pipeline.fit(
            fixed_ctx,
            moving_ctx,
            work_dir=sample_root / tag / "priming",
            stage="similarity",
            fixed_mask=priming_mask,
            affine_iterations=cfg.affine_iterations,
        )
        priming_annotation = pipeline.warped_annotation_in_fixed(moving_ctx, priming["transform"])
        damage = masks_module.missing_tissue_exclusion(
            tissue["tissue"],
            priming_annotation > 0,
            fixed_ctx["cells"],
            overrides,
            halo_px=cfg.exclusion_halo_px,
        )
        fixed_mask = masks_module.registration_mask(tissue["tissue"], damage["exclusion"], artifact)
        report = landmarks_module.build(
            fixed_ctx["cells"],
            moving_ctx["annotation"],
            registry,
            names,
            fixed_ctx["image"].shape,
            slice_index,
        )
        holdout = landmarks_module.choose_holdout(registry, report["report"])
        landmark_fit = landmarks_module.build(
            fixed_ctx["cells"],
            moving_ctx["annotation"],
            registry,
            names,
            fixed_ctx["image"].shape,
            slice_index,
            holdout=holdout,
        )
        return {
            "tag": tag,
            "slice_index": int(slice_index),
            "gradient": float(gradient),
            "medial_extra": int(medial_extra),
            "moving_ctx": moving_ctx,
            "priming_annotation": priming_annotation,
            "damage": damage,
            "fixed_mask": fixed_mask,
            "landmark_report": report["report"],
            "landmark_fit": landmark_fit,
            "holdout": holdout,
        }

    def run_candidate(
        geometry: dict,
        weight: float,
        label: str,
        displacement_reference=None,
        stiffening: tuple = None,
    ) -> dict:
        work_dir = sample_root / geometry["tag"] / label
        parameters = dict(fit_parameters)
        if stiffening is not None:
            bending, grid_scale = stiffening
            parameters["bending_weight"] = bending
            parameters["grid_spacing_voxels"] *= grid_scale
        result = pipeline.fit(
            fixed_ctx,
            geometry["moving_ctx"],
            work_dir=work_dir,
            stage="full",
            fixed_mask=geometry["fixed_mask"],
            landmark_points=geometry["landmark_fit"] if weight > 0 else None,
            landmark_weight=weight,
            **parameters,
        )
        measurements = candidate_metrics(
            result,
            geometry["moving_ctx"],
            cells_raw,
            registry,
            names,
            geometry["holdout"],
            geometry["damage"]["exclusion"],
            geometry["priming_annotation"],
            work_dir,
            displacement_reference=displacement_reference,
        )
        return {
            "label": label,
            "weight": weight,
            "bending_weight": parameters["bending_weight"],
            "grid_spacing_voxels": parameters["grid_spacing_voxels"],
            "result": result,
            "metrics": measurements,
        }

    # Geometry ladder: the searched plane first, then the historical plane. A
    # searched AP level that loses a gated structure is not an improvement, and
    # falling back keeps the rebuild from trading one region for another.
    ladder = [(search["slice_index"], search["gradient"], search["medial_extra"], "searched")]
    if search["slice_index"] != baseline_level or search["gradient"] != 0.0:
        ladder.append((baseline_level, 0.0, search["medial_extra"], "baseline_level"))
    if search["medial_extra"] != config.MIN_MEDIAL_EXTRA:
        ladder.append((baseline_level, 0.0, config.MIN_MEDIAL_EXTRA, "historical_geometry"))

    # Only structural failures disqualify a geometry, because no landmark
    # weight recovers from folding, lost tissue containment, or clamping at the
    # atlas edge. Regional containment losses are exactly what the landmark
    # stage exists to repair, so a geometry carrying them still goes forward.
    geometry_trials = []
    for slice_index, gradient, medial_extra, tag in ladder:
        trial = build_geometry(slice_index, gradient, medial_extra, tag)
        candidate = run_candidate(trial, 0.0, "geometry_only")
        candidate["decision"] = accept.evaluate(candidate["metrics"], baseline, registry)
        candidate["gains"] = accept.improvement(candidate["metrics"], baseline, registry)
        trial["candidate"] = candidate
        geometry_trials.append(trial)
    geometry = max(
        geometry_trials,
        key=lambda item: (
            not item["candidate"]["decision"]["structural_failures"],
            -len(item["candidate"]["decision"]["regional_failures"]),
            np.nan_to_num(item["candidate"]["gains"]["mean_gated_gain"], nan=-99),
        ),
    )

    moving_ctx = geometry["moving_ctx"]
    damage = geometry["damage"]
    fixed_mask = geometry["fixed_mask"]
    holdout = geometry["holdout"]
    landmark_fit = geometry["landmark_fit"]
    landmark_all = {"report": geometry["landmark_report"]}
    priming_annotation = geometry["priming_annotation"]
    search["chosen_geometry"] = geometry["tag"]
    search["slice_index"] = geometry["slice_index"]
    search["gradient"] = geometry["gradient"]
    search["medial_extra"] = geometry["medial_extra"]
    search["obliquity_applied"] = bool(geometry["gradient"] != 0.0)
    search["geometry_trials"] = [
        {
            "tag": trial["tag"],
            "slice_index": trial["slice_index"],
            "gradient": trial["gradient"],
            "medial_extra": trial["medial_extra"],
            "accepted": trial["candidate"]["decision"]["accepted"],
            "reason": trial["candidate"]["decision"]["reason"],
            "tissue_containment_pct": trial["candidate"]["metrics"]["tissue_containment_pct"],
            "mean_gated_gain": trial["candidate"]["gains"]["mean_gated_gain"],
        }
        for trial in geometry_trials
    ]

    # The anchors come from the geometry-only fit, but only if that fit is sound.
    # When it collapsed there are no correctly placed cells to protect, and
    # measuring against it would punish the landmark stage for repairing it, so
    # the historical result becomes the reference instead.
    geometry_candidate = geometry["candidate"]
    geometry_sound = geometry_candidate["metrics"]["tissue_containment_pct"] >= (
        baseline["tissue_containment_pct"] - accept.CONTAINMENT_LOSS_LIMIT
    )
    if geometry_sound:
        reference_cells = geometry_candidate["result"]["cells"]
        reference_annotation_plane = moving_ctx["plane"]["annotation"]
    else:
        reference_cells = baseline_ctx["cells"]
        reference_annotation_plane = baseline_ctx["annotation"]
    geometry_reference = {
        "cells": pipeline.update_baseline_table(cells_raw, reference_cells),
        "anchors": metrics.anchor_cell_index(
            reference_cells, reference_annotation_plane, registry, names
        ),
        "sound_geometry": bool(geometry_sound),
    }
    # Tip anchors refine an already-usable dentate body. On rescue slices whose
    # geometry-only fit collapsed, the same tip pull folds the habenula, so the
    # landmark stage uses the body ladder alone until the slice is seated.
    use_dg_tip = bool(geometry_sound)
    if not use_dg_tip:
        landmark_fit = landmarks_module.build(
            fixed_ctx["cells"],
            moving_ctx["annotation"],
            registry,
            names,
            fixed_ctx["image"].shape,
            search["slice_index"],
            holdout=holdout,
            dg_tip=False,
        )
        landmark_all = {"report": landmark_fit["report"]}
    geometry["landmark_fit"] = landmark_fit
    geometry["landmark_report"] = landmark_all["report"]
    geometry["dg_tip"] = use_dg_tip

    # Stiffening is off by default because it measurably drags whole structures
    # off target, so it is spent only where it buys something: a fit that folds
    # is stiffened until the fold is gone.
    candidates = [geometry_candidate]
    for weight in cfg.landmark_weight_ladder:
        for stiffening in (None,) + tuple(cfg.stiffening_ladder):
            suffix = "" if stiffening is None else f"_stiff{stiffening[0]:.2f}x{stiffening[1]:g}"
            candidate = run_candidate(
                geometry, weight, f"landmarks_{weight:.3f}{suffix}", geometry_reference, stiffening
            )
            candidate["decision"] = accept.evaluate(candidate["metrics"], baseline, registry)
            candidate["gains"] = accept.improvement(candidate["metrics"], baseline, registry)
            candidates.append(candidate)
            if candidate["metrics"].get("jacobian_min", 1.0) > 0:
                break

    accepted = [item for item in candidates if item["decision"]["accepted"]]
    if accepted:
        # Among fits that are essentially tied on gated gain, prefer the one that
        # actually seats the dentate tip and caudoputamen, then protected regions.
        best_gain = max(
            np.nan_to_num(item["gains"]["mean_gated_gain"], nan=-99) for item in accepted
        )
        near_best = [
            item
            for item in accepted
            if np.nan_to_num(item["gains"]["mean_gated_gain"], nan=-99) >= best_gain - 3.0
        ]
        chosen = max(
            near_best,
            key=lambda item: (
                -np.nan_to_num(item["metrics"].get("DG_granule_p90_distance_25um"), nan=99),
                np.nan_to_num(item["metrics"].get("DG_granule_inside_pct"), nan=-1),
                np.nan_to_num(item["metrics"].get("caudoputamen_inside_pct"), nan=-1),
                np.nan_to_num(item["gains"]["mean_protected_gain"], nan=-99),
                np.nan_to_num(item["gains"]["mean_gated_gain"], nan=-99),
                item["metrics"]["tissue_containment_pct"],
            ),
        )
        status = (
            "accepted geometry-only fit"
            if chosen["weight"] == 0
            else f"accepted landmark fit (weight {chosen['weight']:.3f})"
        )
    else:
        chosen = candidates[0]
        status = "rejected by gates; baseline retained"

    # A decision that flips when the damage mask grows or shrinks by a few
    # voxels is not a decision, so slices with masks are re-fitted at both
    # margins and must agree.
    sensitivity = []
    if accepted and damage["exclusion"].any():
        for delta in cfg.mask_sensitivity_px:
            probe_mask = masks_module.registration_mask(
                tissue["tissue"], damage["exclusion"], artifact, grow_px=delta
            )
            probe_dir = sample_root / f"sensitivity_{delta:+d}"
            probe = pipeline.fit(
                fixed_ctx,
                moving_ctx,
                work_dir=probe_dir,
                stage="full",
                fixed_mask=probe_mask,
                landmark_points=landmark_fit if chosen["weight"] > 0 else None,
                landmark_weight=chosen["weight"],
                **fit_parameters,
            )
            probe_metrics = candidate_metrics(
                probe,
                moving_ctx,
                cells_raw,
                registry,
                names,
                holdout,
                damage["exclusion"],
                priming_annotation,
                probe_dir,
                displacement_reference=geometry_reference if chosen["weight"] > 0 else None,
            )
            probe_decision = accept.evaluate(probe_metrics, baseline, registry)
            sensitivity.append(
                {
                    "delta_px": delta,
                    "accepted": probe_decision["accepted"],
                    "reason": probe_decision["reason"],
                    "tissue_containment_pct": probe_metrics["tissue_containment_pct"],
                }
            )
        if not accept.sensitivity_consistent([item["accepted"] for item in sensitivity] + [True]):
            status = "rejected: decision depends on mask margin; baseline retained"

    final_accepted = status.startswith("accepted")
    warped = render.warped_images(moving_ctx, chosen["result"]["transform"])
    render_inputs = config.WORK_ROOT / sample_id / "render_inputs.npz"
    np.savez_compressed(
        render_inputs,
        fixed=fixed_ctx["image"].astype(np.float32),
        moving=moving_ctx["image"].astype(np.float32),
        moving_affine=warped["affine"].astype(np.float32),
        moving_spline=warped["spline"].astype(np.float32),
        borders=moving_ctx["borders"].astype(np.float32),
        borders_spline=warped["borders_spline"].astype(np.float32),
        annotation_plane=(
            moving_ctx["plane"]["annotation"] if final_accepted else baseline_ctx["annotation"]
        ).astype(np.int32),
        tissue=tissue["tissue"],
        cell_supported=tissue["cell_supported"],
        exclusion=damage["exclusion"],
        artifact=artifact,
        fixed_mask=(
            np.ones(fixed_ctx["image"].shape, dtype=np.uint8) if fixed_mask is None else fixed_mask
        ),
        landmarks_fixed=landmark_fit["fixed"],
        landmarks_moving=landmark_fit["moving"],
    )
    # Retaining the baseline has to mean retaining its coordinates: writing the
    # rejected candidate here would put a fit the gates refused into the very
    # table downstream analysis reads.
    if final_accepted:
        cells_out = chosen["result"]["cells"]
    else:
        # The historical table has no fixed-image coordinates, which the QC
        # panels draw, so they come from this run's fixed context.
        cells_out = baseline_ctx["cells"].join(
            fixed_ctx["cells"][["fixed_x", "fixed_y"]], how="left"
        )
    keep = [column for column in RENDER_COLUMNS if column in cells_out.columns]
    cells_out[keep].to_csv(config.WORK_ROOT / sample_id / "cells_final.csv", index=True)
    pipeline.update_baseline_table(cells_raw, cells_out).to_csv(
        config.WORK_ROOT / sample_id / "ccf2d_rebuild.csv", index=True
    )

    record.update(
        {
            "status": status,
            "accepted": final_accepted,
            "slice_number": int(frame.iloc[row_index]["Slice"]),
            "baseline_allen_slice": int(frame.iloc[row_index]["allen_slice_num"]),
            # The geometry of the result that was written out, which is the
            # historical one whenever the candidates were refused.
            "allen_slice_num": search["slice_index"] if final_accepted else baseline_level,
            "ap_gradient": search["gradient"] if final_accepted else 0.0,
            "medial_extra": search["medial_extra"] if final_accepted else config.MIN_MEDIAL_EXTRA,
            "medial_limit": (
                moving_ctx["plane"]["medial_limit"]
                if final_accepted
                else config.ATLAS_MIDLINE + config.MIN_MEDIAL_EXTRA
            ),
            "candidate_allen_slice_num": search["slice_index"],
            "candidate_medial_extra": search["medial_extra"],
            "chosen_geometry": search["chosen_geometry"],
            "obliquity_applied": search["obliquity_applied"] if final_accepted else False,
            "ap_level_disagreement": search["level_disagreement"],
            "ap_score": search["score"],
            "holdout_region": holdout,
            "landmark_weight": chosen["weight"] if final_accepted else 0.0,
            "n_landmarks": int(len(landmark_fit["fixed"])),
            "exclusion_area_px": damage["area_px"],
            "exclusion_components": len(damage["components"]),
            "artifact_area_px": int(artifact.sum()),
            "cell_supported_area_px": int(tissue["cell_supported"].sum()),
            "manual_mask_override": damage["manual_used"] or bool(overrides.get("artifact_rects")),
            "decision_reason": chosen["decision"]["reason"],
            "decision_warnings": "; ".join(chosen["decision"]["warnings"]),
        }
    )
    written = chosen["metrics"] if final_accepted else baseline
    record.update({f"final_{key}": value for key, value in written.items()})
    record.update({f"baseline_{key}": value for key, value in baseline.items()})
    record.update({f"candidate_{key}": value for key, value in chosen["metrics"].items()})

    detail = {
        "sample_id": sample_id,
        "record": record,
        "ap_search": search,
        "landmark_report": landmark_all["report"],
        "holdout_region": holdout,
        "mask": {
            "exclusion_components": damage["components"],
            "artifact_area_px": int(artifact.sum()),
            "overrides_used": overrides,
        },
        "candidates": [
            {
                "label": item["label"],
                "weight": item["weight"],
                "accepted": item["decision"]["accepted"],
                "reason": item["decision"]["reason"],
                "warnings": item["decision"]["warnings"],
                "gains": item["gains"],
                "metrics": item["metrics"],
            }
            for item in candidates
        ],
        "sensitivity": sensitivity,
    }
    with open(config.WORK_ROOT / sample_id / "detail.json", "w") as handle:
        json.dump(config.as_json_safe(detail), handle, indent=2)
    return record


def _worker(sample_id: str, cfg_dict: dict) -> dict:
    cfg = config.RunConfig(**cfg_dict)
    try:
        return process_sample(sample_id, cfg)
    except Exception as error:  # noqa: BLE001 - one bad slice must not kill the run
        return {
            "sample_id": sample_id,
            "status": f"failed: {error}",
            "accepted": False,
            "traceback": traceback.format_exc(),
        }


def render_all(records: list, output_root: Path) -> list:
    """Render every slice once, with the 10 um annotation loaded a single time."""
    registry = config.load_registry()
    names = nb_functions.allen_name_to_annots()
    annotation_10 = atlas_io.load_annotation_10()
    output_root.mkdir(parents=True, exist_ok=True)
    library_records = []
    for record in sorted(records, key=lambda item: item["sample_id"]):
        sample_id = record["sample_id"]
        inputs_path = config.WORK_ROOT / sample_id / "render_inputs.npz"
        cells_path = config.WORK_ROOT / sample_id / "cells_final.csv"
        if not inputs_path.exists() or not cells_path.exists():
            continue
        data = np.load(inputs_path, allow_pickle=False)
        cells = pd.read_csv(cells_path, index_col=0)
        # Figure names carry the Allen level, so a changed level leaves the
        # previous level's files behind and both look current.
        for stale in output_root.glob(f"{sample_id}*"):
            stale.unlink()
        title = (
            f"{sample_id} — {record['status']}\n"
            f"Allen {record['allen_slice_num']} (baseline {record['baseline_allen_slice']}), "
            f"gradient {record['ap_gradient']:+.3f}, medial margin {record['medial_extra']} voxels"
        )
        render.coordinate_figures(
            output_root,
            sample_id,
            title,
            cells,
            data["annotation_plane"],
            annotation_10,
            registry,
            names,
        )
        registration_figure, cells_figure = render.registration_qc_figures(
            output_root,
            sample_id,
            int(record["slice_number"]),
            int(record["allen_slice_num"]),
            record["status"],
            data["fixed"],
            data["moving"],
            data["moving_affine"],
            data["moving_spline"],
            data["borders"],
            data["borders_spline"],
            cells,
        )
        render.mask_qc_figure(
            output_root,
            sample_id,
            data["fixed"],
            data["tissue"],
            data["cell_supported"],
            data["exclusion"],
            data["artifact"],
            data["fixed_mask"],
            title,
        )
        render.landmark_qc_figure(
            output_root,
            sample_id,
            data["fixed"],
            data["moving"],
            {"fixed": data["landmarks_fixed"], "moving": data["landmarks_moving"]},
            title,
        )
        library_records.append(
            {
                "sample_id": sample_id,
                "status": record["status"],
                "registration_figure": registration_figure,
                "cells_figure": cells_figure,
                "surface_figure": f"{sample_id} 3d reg.png",
            }
        )
    render.library(output_root, library_records)
    return library_records


def assemble_final_coordinate_tables(sample_ids: list | None = None) -> dict:
    """Concatenate per-slice rebuild tables into usable ccf2d.csv / ccf3d.csv.

    Historical naming: ccf3d carries primary CCF columns plus optional `_2`
    duplicates; ccf2d is the same cells with only the primary CCF columns.
    Both files are written under RUN_ROOT so downstream analysis has a single
    path rather than thirteen per-slice CSVs.
    """
    samples = list(sample_ids) if sample_ids else list(config.SAMPLES)
    frames = []
    for sample_id in samples:
        path = config.WORK_ROOT / sample_id / "ccf2d_rebuild.csv"
        if not path.exists():
            raise FileNotFoundError(f"missing per-slice rebuild table: {path}")
        frame = pd.read_csv(path, index_col=0)
        if "sample_id" not in frame.columns:
            frame["sample_id"] = sample_id
        frames.append(frame)
    combined = pd.concat(frames, axis=0, ignore_index=False)
    # Primary columns are the rebuild result; keep `_2` in sync so a reader of
    # either set of columns gets the accepted (or retained-baseline) fit.
    for column in ("ccfx", "ccfy", "ccfz", "annotation"):
        duplicate = f"{column}_2"
        if column in combined.columns:
            combined[duplicate] = combined[column]

    ccf3d_path = config.RUN_ROOT / "ccf3d.csv"
    ccf2d_path = config.RUN_ROOT / "ccf2d.csv"
    combined.to_csv(ccf3d_path, index=True)
    drop_secondary = [column for column in ("ccfx_2", "ccfy_2", "ccfz_2", "annotation_2") if column in combined.columns]
    combined.drop(columns=drop_secondary).to_csv(ccf2d_path, index=True)
    return {
        "ccf3d": ccf3d_path,
        "ccf2d": ccf2d_path,
        "n_cells": int(len(combined)),
        "n_samples": int(combined["sample_id"].nunique()) if "sample_id" in combined.columns else len(samples),
    }


def promote(records: list, source: Path, destination: Path) -> list:
    """Replace the current finals with the rebuild figures.

    Refused slices are promoted too: their figures now show the retained
    historical fit and say so, and leaving the previous renderer's output in
    place beside the new ones is how stale finals get mistaken for current.
    """
    destination.mkdir(parents=True, exist_ok=True)
    promoted = []
    for record in records:
        sample_id = record["sample_id"]
        for existing in destination.glob(f"{sample_id}*"):
            existing.unlink()
        for path in sorted(source.glob(f"{sample_id}*")):
            shutil.copy2(path, destination / path.name)
            promoted.append(path.name)
    for name in ("FINAL_FITS_LIBRARY.jpg", "FINAL_FIGURE_MANIFEST.csv"):
        candidate = source / name
        if candidate.exists():
            shutil.copy2(candidate, destination / name)
    legacy = destination / "render_summary.json"
    if legacy.exists():
        legacy.unlink()
    _write_final_readme(records, destination)
    return promoted


def _write_final_readme(records: list, destination: Path) -> None:
    accepted = [r["sample_id"] for r in records if r.get("accepted")]
    retained = [r["sample_id"] for r in records if not r.get("accepted")]
    lines = [
        "# Allen CCF registration — current final figures",
        "",
        f"Rebuilt {datetime.now():%Y-%m-%d} by `python -m allen_ccf_repro.cli --all --promote`.",
        "Every figure here comes from that single run; nothing is carried over from",
        "the earlier iteration directories.",
        "",
        "Seven figures per slice, plus a three-column `FINAL_FITS_LIBRARY.jpg` that",
        "stacks registration diagnostics, cells+borders, and the 3d CCF surface:",
        "",
        "1. `2d reg.png` — cells with fixed-AP Allen contours",
        "2. `3d reg.png` — cells with curved-surface Allen contours",
        "3. `regional targets.png` — DG, MH, LH, and ependymal diagnostics",
        "4. `mask QC.png` — tissue, cell-supported, missing-tissue, and artifact masks",
        "5. `landmark QC.png` — the correspondence points driving the spline",
        "6. `reg_5.jpg` — fixed/affine/spline registration diagnostic",
        "7. `reg_5_cells_all.jpg` — Allen borders alone and over all cells",
        "",
        f"Rebuilt fit accepted ({len(accepted)}): " + (", ".join(sorted(accepted)) or "none"),
        f"Historical baseline retained ({len(retained)}): "
        + (", ".join(sorted(retained)) or "none"),
        "",
        "Retained slices are rendered from their baseline fit and their titles say so.",
        "Per-slice geometry, gates, and metrics are in `FINAL_FIGURE_MANIFEST.csv` and",
        "the run manifest under the rebuild work root.",
    ]
    (destination / "README.md").write_text("\n".join(lines) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Rebuild all Allen CCF registrations.")
    parser.add_argument("--all", action="store_true", help="Process every sample.")
    parser.add_argument("--samples", nargs="*", default=None, help="Subset of sample ids.")
    parser.add_argument("--jobs", type=int, default=None, help="Parallel slices.")
    parser.add_argument("--skip-fits", action="store_true", help="Reuse existing fits and only render.")
    parser.add_argument("--no-render", action="store_true", help="Skip figure rendering.")
    parser.add_argument("--promote", action="store_true", help="Copy accepted figures over the current finals.")
    parser.add_argument("--spline-iterations", type=int, default=None)
    parser.add_argument("--search-iterations", type=int, default=None)
    args = parser.parse_args(argv)

    cfg = config.RunConfig()
    if args.jobs:
        cfg.jobs = args.jobs
    if args.spline_iterations:
        cfg.spline_iterations = args.spline_iterations
    if args.search_iterations:
        cfg.search_iterations = args.search_iterations
    if args.samples:
        cfg.samples = list(args.samples)
    elif not args.all:
        parser.error("pass --all or --samples")

    config.ensure_directories()
    registry = config.load_registry()

    records: list = []
    if args.skip_fits:
        for sample_id in cfg.samples:
            path = config.WORK_ROOT / sample_id / "detail.json"
            if path.exists():
                with open(path) as handle:
                    records.append(json.load(handle)["record"])
    else:
        cfg_dict = cfg.to_dict()
        if cfg.jobs <= 1:
            for sample_id in cfg.samples:
                print(f"[{sample_id}] start", flush=True)
                records.append(_worker(sample_id, cfg_dict))
                print(f"[{sample_id}] {records[-1]['status']}", flush=True)
        else:
            with ProcessPoolExecutor(max_workers=cfg.jobs) as pool:
                futures = {
                    pool.submit(_worker, sample_id, cfg_dict): sample_id for sample_id in cfg.samples
                }
                for future in as_completed(futures):
                    record = future.result()
                    records.append(record)
                    print(f"[{record['sample_id']}] {record['status']}", flush=True)

    table = pd.DataFrame(records).sort_values("sample_id")
    table.to_csv(config.RUN_ROOT / "rebuild_decisions.csv", index=False)

    # Usable project-level coordinate tables (not only per-slice rebuild CSVs).
    assembled = assemble_final_coordinate_tables([record["sample_id"] for record in records])
    print(
        f"wrote {assembled['n_cells']} cells across {assembled['n_samples']} samples to "
        f"{assembled['ccf2d'].name} and {assembled['ccf3d'].name}",
        flush=True,
    )

    if not args.no_render:
        library_records = render_all(records, config.REBUILD_FIGURES)
        pd.DataFrame(library_records).to_csv(
            config.REBUILD_FIGURES / "FINAL_FIGURE_MANIFEST.csv", index=False
        )

    if args.promote:
        promoted = promote(records, config.REBUILD_FIGURES, config.PROMOTED_FIGURES)
        print(f"promoted {len(promoted)} figures", flush=True)

    provenance.write_manifest(
        config.RUN_ROOT / "run_manifest.json",
        cfg,
        registry,
        [
            {
                key: value
                for key, value in record.items()
                if key
                in (
                    "sample_id",
                    "status",
                    "accepted",
                    "allen_slice_num",
                    "baseline_allen_slice",
                    "ap_gradient",
                    "medial_extra",
                    "holdout_region",
                    "landmark_weight",
                    "n_landmarks",
                    "exclusion_area_px",
                    "decision_reason",
                )
            }
            for record in sorted(records, key=lambda item: item["sample_id"])
        ],
    )
    accepted = int(table["accepted"].sum()) if "accepted" in table else 0
    print(f"{accepted}/{len(table)} slices accepted", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
