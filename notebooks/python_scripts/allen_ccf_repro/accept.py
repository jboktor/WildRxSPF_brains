"""Per-region acceptance gates for a candidate registration."""

from __future__ import annotations

import numpy as np

ANCHOR_DISPLACEMENT_LIMIT = 8.0
ANCHOR_MIN_CELLS = 200
CONTAINMENT_LOSS_LIMIT = 1.0
PROTECTED_LOSS_LIMIT = 2.0
HOLDOUT_RELATIVE_LIMIT = 0.20
THICKNESS_RATIO_LIMIT = 0.85
COMPRESSION_RATIO_LIMIT = 1.25
BOUNDARY_CLAMP_LIMIT = 3.0
GATED_REGION_LOSS_LIMIT = 3.0
# Distances come from a voxel grid, so a median that shifts by less than one
# voxel is quantisation rather than movement. Protected regions are held to a
# stricter containment loss instead of a finer distance than can be measured.
REGION_DISTANCE_LIMIT = 1.0
PROTECTED_DISTANCE_LIMIT = 1.0
GROSS_REGION_LOSS_LIMIT = 20.0
# Per-region gates forgive losses that the distance check calls quantisation,
# which lets a candidate accumulate many forgiven losses and still be worse
# than the fit it would replace. These bound the total.
NET_GATED_GAIN_LIMIT = -1.0
NET_PROTECTED_GAIN_LIMIT = -5.0


def _value(mapping: dict, key: str, default=np.nan) -> float:
    value = mapping.get(key, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def evaluate(candidate: dict, baseline: dict, registry: dict) -> dict:
    """Return the accept/reject decision with an explicit reason for each gate.

    Failures are split into structural ones, which no later stage can repair,
    and regional ones, which the landmark stage can still fix.
    """
    structural: list = []
    regional: list = []
    warnings: list = []

    # With only a handful of correctly placed cells to start from, their median
    # displacement says more about which few cells happened to land right than
    # about the candidate, so the gate stands down and reports instead.
    anchor = _value(candidate, "anchor_median_displacement_25um")
    anchor_count = _value(candidate, "n_anchor_cells", 0)
    if np.isfinite(anchor) and anchor > ANCHOR_DISPLACEMENT_LIMIT:
        message = f"already-correct cells moved {anchor:.2f} voxels > {ANCHOR_DISPLACEMENT_LIMIT}"
        if anchor_count >= ANCHOR_MIN_CELLS:
            structural.append(message)
        else:
            warnings.append(f"{message} (only {int(anchor_count)} such cells)")

    containment_new = _value(candidate, "tissue_containment_pct")
    containment_old = _value(baseline, "tissue_containment_pct")
    if np.isfinite(containment_new) and np.isfinite(containment_old):
        loss = containment_old - containment_new
        if loss > CONTAINMENT_LOSS_LIMIT:
            structural.append(
                f"tissue containment fell {loss:.2f} pp ({containment_old:.2f} to {containment_new:.2f})"
            )

    # Containment is a binary in/out test against a structure boundary, so a
    # sub-voxel shift can flip a large fraction of cells that sit packed against
    # that boundary. Median distance measures the physical error instead, and a
    # containment loss only counts as a failure when the cells actually moved
    # away from the structure, or when the loss is so large that the region
    # clearly slipped off its target.
    for name, entry in registry["regions"].items():
        new = _value(candidate, f"{name}_inside_pct")
        old = _value(baseline, f"{name}_inside_pct")
        if not (np.isfinite(new) and np.isfinite(old)):
            continue
        loss = old - new
        distance_new = _value(candidate, f"{name}_median_distance_25um")
        distance_old = _value(baseline, f"{name}_median_distance_25um")
        drift = (
            distance_new - distance_old
            if np.isfinite(distance_new) and np.isfinite(distance_old)
            else float("nan")
        )
        moved = (np.isfinite(drift) and drift > REGION_DISTANCE_LIMIT) or loss > GROSS_REGION_LOSS_LIMIT
        protected_moved = (
            np.isfinite(drift) and drift > PROTECTED_DISTANCE_LIMIT
        ) or loss > GROSS_REGION_LOSS_LIMIT
        detail = f"lost {loss:.2f} pp containment"
        if np.isfinite(drift):
            detail += f", median distance {distance_old:.2f} to {distance_new:.2f}"
        if entry.get("protected") and loss > PROTECTED_LOSS_LIMIT and protected_moved:
            regional.append(f"protected region {name} {detail}")
        elif entry.get("gate") and loss > GATED_REGION_LOSS_LIMIT and moved:
            regional.append(f"gated region {name} {detail}")
        elif loss > GATED_REGION_LOSS_LIMIT:
            warnings.append(f"{name} {detail}")

    net = improvement(candidate, baseline, registry)
    gated_gain = net["mean_gated_gain"]
    if np.isfinite(gated_gain) and gated_gain < NET_GATED_GAIN_LIMIT:
        regional.append(
            f"net regression across gated regions ({gated_gain:+.2f} pp mean containment)"
        )
    protected_gain = net["mean_protected_gain"]
    if np.isfinite(protected_gain) and protected_gain < NET_PROTECTED_GAIN_LIMIT:
        regional.append(
            f"net regression across protected regions ({protected_gain:+.2f} pp mean containment)"
        )

    holdout = candidate.get("holdout_region")
    if holdout:
        new_distance = _value(candidate, "holdout_median_distance_25um")
        old_distance = _value(baseline, f"{holdout}_median_distance_25um")
        if np.isfinite(new_distance) and np.isfinite(old_distance) and old_distance > 0:
            change = (new_distance - old_distance) / old_distance
            if change > HOLDOUT_RELATIVE_LIMIT:
                regional.append(
                    f"holdout {holdout} median distance worsened {100 * change:.0f}%"
                )
        new_inside = _value(candidate, "holdout_inside_pct")
        old_inside = _value(baseline, f"{holdout}_inside_pct")
        if np.isfinite(new_inside) and np.isfinite(old_inside) and old_inside > 0:
            relative = (old_inside - new_inside) / old_inside
            if relative > HOLDOUT_RELATIVE_LIMIT:
                regional.append(f"holdout {holdout} containment fell {100 * relative:.0f}%")

    for label in ("cortex_thickness_ratio_near_mask", "callosum_thickness_ratio_near_mask"):
        ratio = _value(candidate, label)
        if np.isfinite(ratio) and ratio < THICKNESS_RATIO_LIMIT:
            structural.append(f"{label} {ratio:.2f} < {THICKNESS_RATIO_LIMIT}")

    jacobian_min = _value(candidate, "jacobian_min")
    if np.isfinite(jacobian_min) and jacobian_min <= 0:
        structural.append("deformation folds (non-positive Jacobian)")
    global_mean = _value(candidate, "jacobian_mean")
    in_region = _value(candidate, "jacobian_p95_in_region")
    if np.isfinite(global_mean) and np.isfinite(in_region) and global_mean > 0:
        ratio = in_region / global_mean
        if ratio > COMPRESSION_RATIO_LIMIT:
            structural.append(f"atlas compressed into excluded tissue (ratio {ratio:.2f})")

    clamp = _value(candidate, "boundary_clamp_ratio")
    if np.isfinite(clamp) and clamp > BOUNDARY_CLAMP_LIMIT:
        structural.append(f"cells stacked against the medial atlas edge (ratio {clamp:.2f})")
    protected_clamp = _value(candidate, "protected_boundary_clamp_ratio")
    if np.isfinite(protected_clamp) and protected_clamp > BOUNDARY_CLAMP_LIMIT:
        structural.append(
            f"habenular cells stacked against the medial atlas edge (ratio {protected_clamp:.2f})"
        )
    overhang = _value(candidate, "medial_overhang_25um")
    if np.isfinite(overhang) and overhang > 10:
        warnings.append(f"tissue extends {overhang:.1f} voxels past the atlas midline")

    failures = structural + regional
    return {
        "accepted": not failures,
        "failures": failures,
        "structural_failures": structural,
        "regional_failures": regional,
        "warnings": warnings,
        "reason": "; ".join(failures) if failures else "all gates passed",
    }


def improvement(candidate: dict, baseline: dict, registry: dict) -> dict:
    """Summarise where the candidate actually gained, for the decision record."""
    gains = {}
    for name, entry in registry["regions"].items():
        key = f"{name}_inside_pct"
        new = _value(candidate, key)
        old = _value(baseline, key)
        if np.isfinite(new) and np.isfinite(old):
            gains[name] = new - old
    gains["tissue_containment_pct"] = _value(candidate, "tissue_containment_pct") - _value(
        baseline, "tissue_containment_pct"
    )
    protected = [
        name for name, entry in registry["regions"].items() if entry.get("protected") and name in gains
    ]
    gated = [name for name, entry in registry["regions"].items() if entry.get("gate") and name in gains]
    return {
        "gains": gains,
        "mean_gated_gain": float(np.mean([gains[name] for name in gated])) if gated else float("nan"),
        "mean_protected_gain": float(np.mean([gains[name] for name in protected]))
        if protected
        else float("nan"),
    }


def sensitivity_consistent(decisions) -> bool:
    """Mask sensitivity check: the decision must not depend on the mask margin."""
    values = [bool(item) for item in decisions if item is not None]
    return len(set(values)) <= 1
