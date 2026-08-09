#!/usr/bin/env python
"""Continue the isolated all-slice landmark sweep through adaptive iterations.

This runner consumes iterations 11/12, runs only the plan-authorized adaptive
branches (13--17), re-evaluates every candidate with protected-region loss
gates, and writes one final decision per sample. Canonical coordinates,
workbooks, and root-level final images are read-only.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
from pathlib import Path
from typing import Any

import nrrd
import numpy as np
import pandas as pd
import skimage.io


BASE = Path(__file__).resolve().parents[2]
ALL_SLICE_SCRIPT = Path(__file__).with_name("allen_ccf_all_slice_landmark_sweep.py")
TEST_ROOT = BASE / "data/interim/registration/Allen_CCF_regional_tests"
FIGURE_ROOT = BASE / "figures/Allen_CCF_alignment_optimized"
WORKBOOK = TEST_ROOT / "slice_positions_25um_all_slices_regional_landmarks.xlsx"
METRICS_CSV = TEST_ROOT / "all_slices_regional_landmark_metrics.csv"
DECISIONS_CSV = TEST_ROOT / "all_slices_regional_landmark_decisions.csv"
FINAL_DECISIONS_CSV = TEST_ROOT / "all_slices_regional_landmark_final_decisions.csv"
REFERENCE_ROOT = BASE / "data/input/allen_registration_ref"

ADAPTIVE_SHEETS = {
    13: ("landmarks001", 0.01),
    14: ("landmarks005", 0.05),
    15: ("landmarks002", 0.02),
    16: ("landmarks004", 0.04),
    17: ("roi4_run1_lh10_003", 0.03),
}


def load_all_slice_module():
    spec = importlib.util.spec_from_file_location("allen_ccf_all_slice", ALL_SLICE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {ALL_SLICE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def create_adaptive_sheets(sweep) -> None:
    """Append candidate parameter sheets without replacing prior history."""
    source = pd.read_excel(WORKBOOK, sheet_name="landmarks003")
    sheets: dict[str, pd.DataFrame] = {}
    for iteration, (sheet_name, weight) in ADAPTIVE_SHEETS.items():
        candidate = source.copy()
        candidate["landmark_weight"] = weight
        if iteration == 17:
            row_index = candidate.index[candidate["cell_metadata"].eq("roi4_run1")]
            if len(row_index) != 1:
                raise RuntimeError("roi4_run1 is missing or duplicated in the isolated workbook")
            index = int(row_index[0])
            cell_types = sweep.parse_list(candidate.at[index, "cell_types"])
            annotations = sweep.parse_list(candidate.at[index, "annots_to_amplify"])
            if sweep.TARGET_CLASSES["LH"] not in cell_types:
                cell_types.append(sweep.TARGET_CLASSES["LH"])
            for annotation_id in sweep.TARGET_IDS["LH"]:
                if annotation_id not in annotations:
                    annotations.append(annotation_id)
            candidate.at[index, "cell_types"] = repr(cell_types)
            candidate.at[index, "annots_to_amplify"] = repr(annotations)
        sheets[sheet_name] = candidate
    with pd.ExcelWriter(WORKBOOK, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False)


def protected_loss_reasons(baseline: pd.Series, candidate: pd.Series) -> list[str]:
    """Protect well-contained MH/LH populations, including low-count classes."""
    reasons: list[str] = []
    for label in ("MH", "LH"):
        baseline_inside = baseline[f"{label}_inside_pct"]
        candidate_inside = candidate[f"{label}_inside_pct"]
        cell_count = baseline[f"{label}_n"]
        if (
            pd.notna(baseline_inside)
            and pd.notna(candidate_inside)
            and cell_count >= 10
            and baseline_inside >= 60
            and candidate_inside < baseline_inside - 30
        ):
            reasons.append(
                f"{label} protected loss {baseline_inside:.1f}->{candidate_inside:.1f} pp"
            )
    return reasons


def evaluate_candidate(
    baseline: pd.Series,
    candidate: pd.Series,
    eligibility: pd.Series,
) -> tuple[str, str]:
    """Apply hard safety gates and applicable regional correction gates."""
    anchor = float(candidate["left_top_anchor_median_displacement_25um"])
    hard_rejects: list[str] = []
    if anchor > 8.0:
        hard_rejects.append(f"anchor {anchor:.2f}>8")
    if (
        pd.notna(baseline["all_tissue_containment_pct"])
        and pd.notna(candidate["all_tissue_containment_pct"])
        and candidate["all_tissue_containment_pct"]
        < baseline["all_tissue_containment_pct"] - 1.0
    ):
        hard_rejects.append("tissue containment dropped >1 pp")

    hard_rejects.extend(protected_loss_reasons(baseline, candidate))
    for label in ("MH", "LH"):
        if bool(eligibility[f"{label}_eligible"]):
            baseline_inside = baseline[f"{label}_inside_pct"]
            candidate_inside = candidate[f"{label}_inside_pct"]
            if (
                pd.notna(baseline_inside)
                and pd.notna(candidate_inside)
                and candidate_inside < baseline_inside - 10
            ):
                hard_rejects.append(f"{label} inside dropped >10 pp")

    clear_win = False
    if pd.notna(baseline["DG_p90_distance_25um"]) and pd.notna(candidate["DG_p90_distance_25um"]):
        clear_win |= candidate["DG_p90_distance_25um"] <= 0.8 * baseline["DG_p90_distance_25um"]
    for label in ("MH", "LH"):
        if pd.notna(baseline[f"{label}_inside_pct"]) and pd.notna(candidate[f"{label}_inside_pct"]):
            clear_win |= candidate[f"{label}_inside_pct"] >= baseline[f"{label}_inside_pct"] + 20
    if candidate["median_displacement_25um"] > 12 and not clear_win:
        hard_rejects.append("median displacement >12 without clear regional win")

    if hard_rejects:
        decision = "reject_overwarp" if anchor > 8 else "reject_hard_gate"
        return decision, "; ".join(dict.fromkeys(hard_rejects))

    target_gates: dict[str, bool] = {}
    classes = set(str(candidate["landmark_classes"]).split(","))
    if bool(eligibility["DG_eligible"]) and "DG" in classes:
        baseline_p90 = baseline["DG_p90_distance_25um"]
        candidate_p90 = candidate["DG_p90_distance_25um"]
        baseline_inside = baseline["DG_inside_pct"]
        candidate_inside = candidate["DG_inside_pct"]
        target_gates["DG"] = bool(
            candidate_p90 <= max(3.0, 0.35 * baseline_p90)
            and (candidate_inside >= baseline_inside + 8 or candidate_inside >= 68)
        )
    for label, improvement in (("LH", 30.0), ("MH", 20.0)):
        if bool(eligibility[f"{label}_eligible"]) and label in classes:
            baseline_inside = baseline[f"{label}_inside_pct"]
            candidate_inside = candidate[f"{label}_inside_pct"]
            target_gates[label] = bool(
                candidate_inside >= 60 or candidate_inside >= baseline_inside + improvement
            )

    if anchor <= 6.5 and target_gates and all(target_gates.values()):
        return "accept", "all applicable regional, protected-region, and anchor gates passed"
    if anchor <= 8.0 and target_gates and not all(target_gates.values()) and clear_win:
        failed = ",".join(label for label, passed in target_gates.items() if not passed)
        return "under_corrected", f"safe improvement but failed target gates: {failed}"
    failed = ",".join(label for label, passed in target_gates.items() if not passed)
    return "reject_no_benefit", f"failed target gates: {failed or 'none'}"


def candidate_score(baseline: pd.Series, candidate: pd.Series, eligibility: pd.Series) -> float:
    """Rank multiple accepted candidates reproducibly against their baseline."""
    components: list[float] = []
    if bool(eligibility["DG_eligible"]) and pd.notna(baseline["DG_p90_distance_25um"]):
        denominator = max(float(baseline["DG_p90_distance_25um"]), 1.0)
        components.append(
            (baseline["DG_p90_distance_25um"] - candidate["DG_p90_distance_25um"])
            / denominator
            + (candidate["DG_inside_pct"] - baseline["DG_inside_pct"]) / 100
        )
    for label in ("MH", "LH"):
        if bool(eligibility[f"{label}_eligible"]):
            components.append(
                (candidate[f"{label}_inside_pct"] - baseline[f"{label}_inside_pct"]) / 100
            )
    regional_score = float(np.mean(components)) if components else 0.0
    return regional_score - float(candidate["left_top_anchor_median_displacement_25um"]) / 20


def run_or_resume(
    sweep,
    sample_id: str,
    iteration: int,
    weight: float,
    sheet_name: str,
    baseline_cells: pd.DataFrame,
    eligibility: pd.Series,
    atlas: dict[str, Any],
    resume: bool,
    min_cells: dict[str, int] | None = None,
) -> dict[str, Any]:
    run_dir = TEST_ROOT / f"iteration_{iteration:02d}" / sample_id
    metrics_path = run_dir / "run_metrics.json"
    coordinate_path = run_dir / f"{sample_id}_ccf2d.csv"
    if resume and metrics_path.exists() and coordinate_path.exists():
        return json.loads(metrics_path.read_text())
    metrics, _ = sweep.run_candidate(
        sample_id,
        weight,
        baseline_cells,
        eligibility,
        atlas["annotation_25"],
        atlas["annotation_10"],
        atlas["nissl"],
        atlas["borders"],
        atlas["regional"],
        iteration=iteration,
        sheet_name=sheet_name,
        min_cells=min_cells,
    )
    return metrics


def eligibility_for_iteration(
    base_eligibility: pd.Series, sample_id: str, iteration: int
) -> pd.Series:
    eligibility = base_eligibility.copy()
    if sample_id == "roi4_run1" and iteration == 17:
        eligibility["LH_eligible"] = True
        classes = set(str(eligibility["eligible_classes"]).split(","))
        classes.add("LH")
        eligibility["eligible_classes"] = ",".join(sorted(item for item in classes if item))
    return eligibility


def build_history_and_final(
    all_metrics: pd.DataFrame,
    baseline_metrics: pd.DataFrame,
    eligibility: pd.DataFrame,
    prior_decisions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline_by_sample = baseline_metrics.set_index("sample_id")
    eligibility_by_sample = eligibility.set_index("sample_id")
    prior_lookup = {
        (row.sample_id, int(row.iteration)): row.decision
        for row in prior_decisions.dropna(subset=["iteration"]).itertuples()
    }
    history_records: list[dict[str, Any]] = []
    for candidate in all_metrics.sort_values(["sample_id", "iteration"]).itertuples(index=False):
        candidate_series = pd.Series(candidate._asdict())
        sample_id = candidate_series["sample_id"]
        iteration = int(candidate_series["iteration"])
        current_eligibility = eligibility_for_iteration(
            eligibility_by_sample.loc[sample_id], sample_id, iteration
        )
        decision, reason = evaluate_candidate(
            baseline_by_sample.loc[sample_id], candidate_series, current_eligibility
        )
        history_records.append(
            {
                "sample_id": sample_id,
                "iteration": iteration,
                "weight": candidate_series["weight"],
                "prior_decision": prior_lookup.get((sample_id, iteration), np.nan),
                "decision": decision,
                "reason": reason,
                "candidate_score": (
                    candidate_score(
                        baseline_by_sample.loc[sample_id],
                        candidate_series,
                        current_eligibility,
                    )
                    if decision == "accept"
                    else np.nan
                ),
            }
        )
    history = pd.DataFrame(history_records)

    final_records: list[dict[str, Any]] = []
    for sample_id in baseline_metrics["sample_id"]:
        sample_history = history[history["sample_id"].eq(sample_id)]
        accepted = sample_history[sample_history["decision"].eq("accept")]
        if len(accepted):
            selected = accepted.sort_values(
                ["candidate_score", "iteration"], ascending=[False, True]
            ).iloc[0]
            selected_metrics = all_metrics[
                all_metrics["sample_id"].eq(sample_id)
                & all_metrics["iteration"].eq(selected["iteration"])
            ].iloc[0]
            final_records.append(
                {
                    "sample_id": sample_id,
                    "final_decision": "accept_candidate",
                    "accepted_iteration": int(selected["iteration"]),
                    "accepted_weight": selected["weight"],
                    "candidate_score": selected["candidate_score"],
                    "anchor_displacement_25um": selected_metrics[
                        "left_top_anchor_median_displacement_25um"
                    ],
                    "DG_p90_distance_25um": selected_metrics["DG_p90_distance_25um"],
                    "DG_inside_pct": selected_metrics["DG_inside_pct"],
                    "MH_inside_pct": selected_metrics["MH_inside_pct"],
                    "LH_inside_pct": selected_metrics["LH_inside_pct"],
                    "reason": selected["reason"],
                }
            )
        else:
            rejection_summary = " | ".join(
                f"i{int(row.iteration)}:{row.decision} ({row.reason})"
                for row in sample_history.itertuples()
            )
            baseline = baseline_by_sample.loc[sample_id]
            final_records.append(
                {
                    "sample_id": sample_id,
                    "final_decision": "retain_baseline",
                    "accepted_iteration": np.nan,
                    "accepted_weight": np.nan,
                    "candidate_score": np.nan,
                    "anchor_displacement_25um": 0.0,
                    "DG_p90_distance_25um": baseline["DG_p90_distance_25um"],
                    "DG_inside_pct": baseline["DG_inside_pct"],
                    "MH_inside_pct": baseline["MH_inside_pct"],
                    "LH_inside_pct": baseline["LH_inside_pct"],
                    "reason": rejection_summary,
                }
            )
    return history, pd.DataFrame(final_records)


def ensure_accepted_artifacts(final_decisions: pd.DataFrame) -> None:
    """Mark selected coordinates and verify all nested selected figures."""
    for row in final_decisions[final_decisions["final_decision"].eq("accept_candidate")].itertuples():
        iteration = int(row.accepted_iteration)
        run_dir = TEST_ROOT / f"iteration_{iteration:02d}" / row.sample_id
        source = run_dir / f"{row.sample_id}_ccf2d.csv"
        accepted = run_dir / f"{row.sample_id}_accepted_ccf2d.csv"
        if not accepted.exists():
            shutil.copy2(source, accepted)
        figure_dir = FIGURE_ROOT / f"iteration_{iteration:02d}" / row.sample_id
        expected = [
            f"{row.sample_id} 2d reg.png",
            f"{row.sample_id} 3d reg.png",
            f"{row.sample_id} regional targets.png",
            f"{row.sample_id} regional registration QC.png",
        ]
        missing = [name for name in expected if not (figure_dir / name).exists()]
        if missing:
            raise RuntimeError(f"Missing accepted figures for {row.sample_id}: {missing}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resume", action="store_true", help="Reuse completed iteration 13--17 run products"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sweep = load_all_slice_module()
    create_adaptive_sheets(sweep)

    baseline_metrics = pd.read_csv(sweep.BASELINE_METRICS_CSV)
    eligibility = pd.read_csv(sweep.ELIGIBILITY_CSV)
    prior_metrics = pd.read_csv(METRICS_CSV)
    prior_decisions = pd.read_csv(DECISIONS_CSV)
    samples = baseline_metrics["sample_id"].tolist()
    baseline_cells = sweep.load_baseline_cells(set(samples))
    baseline_by_sample = baseline_metrics.set_index("sample_id")
    eligibility_by_sample = eligibility.set_index("sample_id")

    regional = sweep.load_regional_module()
    annotation_25, _ = nrrd.read(
        str(REFERENCE_ROOT / "annotation_25.nrrd"), index_order="F"
    )
    annotation_25 = annotation_25[:, :, :236]
    annotation_10, _ = nrrd.read(
        str(REFERENCE_ROOT / "annotation_10.nrrd"), index_order="F"
    )
    nissl, _ = nrrd.read(str(REFERENCE_ROOT / "ara_nissl_25.nrrd"), index_order="F")
    nissl = nissl[:, :, :236]
    borders = skimage.io.imread(
        str(REFERENCE_ROOT / "annotations_25_border_only.tif")
    )[:, :, :236]
    atlas = {
        "regional": regional,
        "annotation_25": annotation_25,
        "annotation_10": annotation_10,
        "nissl": nissl,
        "borders": borders,
    }

    iteration11_decisions = prior_decisions[prior_decisions["iteration"].eq(11)]
    overwarp_samples = iteration11_decisions[
        iteration11_decisions["decision"].eq("reject_overwarp")
    ]["sample_id"].tolist()
    under_corrected_samples = iteration11_decisions[
        iteration11_decisions["decision"].eq("under_corrected")
    ]["sample_id"].tolist()
    adaptive_records: list[dict[str, Any]] = []

    for sample_id in overwarp_samples:
        print(f"{sample_id}: iteration 13 weight 0.01", flush=True)
        adaptive_records.append(
            run_or_resume(
                sweep,
                sample_id,
                13,
                0.01,
                "landmarks001",
                baseline_cells,
                eligibility_by_sample.loc[sample_id],
                atlas,
                args.resume,
            )
        )

    for sample_id in under_corrected_samples:
        current_eligibility = eligibility_by_sample.loc[sample_id]
        print(f"{sample_id}: iteration 14 weight 0.05", flush=True)
        metrics14 = run_or_resume(
            sweep,
            sample_id,
            14,
            0.05,
            "landmarks005",
            baseline_cells,
            current_eligibility,
            atlas,
            args.resume,
        )
        adaptive_records.append(metrics14)
        decision14, _ = evaluate_candidate(
            baseline_by_sample.loc[sample_id], pd.Series(metrics14), current_eligibility
        )
        if decision14 == "under_corrected":
            print(f"{sample_id}: iteration 15 weight 0.02", flush=True)
            adaptive_records.append(
                run_or_resume(
                    sweep,
                    sample_id,
                    15,
                    0.02,
                    "landmarks002",
                    baseline_cells,
                    current_eligibility,
                    atlas,
                    args.resume,
                )
            )
        elif decision14 == "reject_overwarp":
            print(f"{sample_id}: iteration 16 weight 0.04", flush=True)
            adaptive_records.append(
                run_or_resume(
                    sweep,
                    sample_id,
                    16,
                    0.04,
                    "landmarks004",
                    baseline_cells,
                    current_eligibility,
                    atlas,
                    args.resume,
                )
            )

    roi4_eligibility = eligibility_for_iteration(
        eligibility_by_sample.loc["roi4_run1"], "roi4_run1", 17
    )
    print("roi4_run1: iteration 17 weight 0.03 with LH>=10", flush=True)
    adaptive_records.append(
        run_or_resume(
            sweep,
            "roi4_run1",
            17,
            0.03,
            "roi4_run1_lh10_003",
            baseline_cells,
            roi4_eligibility,
            atlas,
            args.resume,
            min_cells={"LH": 10},
        )
    )

    adaptive_metrics = pd.DataFrame(adaptive_records)
    adaptive_metrics.to_csv(
        TEST_ROOT / "all_slices_regional_landmark_adaptive_metrics.csv", index=False
    )
    all_metrics = (
        pd.concat([prior_metrics, adaptive_metrics], ignore_index=True)
        .drop_duplicates(["sample_id", "iteration"], keep="last")
        .sort_values(["sample_id", "iteration"])
        .reset_index(drop=True)
    )
    history, final_decisions = build_history_and_final(
        all_metrics, baseline_metrics, eligibility, prior_decisions
    )
    all_metrics.to_csv(METRICS_CSV, index=False)
    history.to_csv(DECISIONS_CSV, index=False)
    final_decisions.to_csv(FINAL_DECISIONS_CSV, index=False)
    ensure_accepted_artifacts(final_decisions)

    with pd.ExcelWriter(
        WORKBOOK, engine="openpyxl", mode="a", if_sheet_exists="replace"
    ) as writer:
        all_metrics.to_excel(writer, sheet_name="candidate_metrics", index=False)
        history.to_excel(writer, sheet_name="decision_history", index=False)
        final_decisions.to_excel(writer, sheet_name="final_decisions", index=False)

    print("\nFINAL DECISIONS", flush=True)
    print(final_decisions.to_string(index=False), flush=True)
    print("\nCANDIDATE HISTORY", flush=True)
    print(history.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
