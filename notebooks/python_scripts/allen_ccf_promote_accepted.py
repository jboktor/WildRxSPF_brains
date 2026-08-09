#!/usr/bin/env python
"""Publish accepted per-slice registrations into the optimized figure root."""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops, ImageDraw


BASE = Path(__file__).resolve().parents[2]
FIGURE_ROOT = BASE / "figures/Allen_CCF_alignment_optimized"
CANONICAL_ROOT = BASE / "figures/Allen_CCF_alignment"
TEST_ROOT = BASE / "data/interim/registration/Allen_CCF_regional_tests"
DEFAULT_DECISIONS = TEST_ROOT / "all_slices_regional_landmark_decisions.csv"
MANIFEST_PATH = FIGURE_ROOT / "accepted_manifest.json"
PROMOTION_LOG = FIGURE_ROOT / "accepted_promotions.csv"


def _guard_optimized_path(path: Path) -> None:
    resolved = path.resolve()
    if resolved == CANONICAL_ROOT.resolve() or CANONICAL_ROOT.resolve() in resolved.parents:
        raise ValueError(f"Refusing to write canonical output: {resolved}")
    if resolved != FIGURE_ROOT.resolve() and FIGURE_ROOT.resolve() not in resolved.parents:
        raise ValueError(f"Promotion target is outside optimized figures: {resolved}")


def load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"updated_at": None, "samples": {}}


def save_manifest(manifest: dict) -> None:
    _guard_optimized_path(MANIFEST_PATH)
    manifest["updated_at"] = datetime.now(timezone.utc).isoformat()
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def _trim_white(image: Image.Image) -> Image.Image:
    rgb = image.convert("RGB")
    background = Image.new("RGB", rgb.size, "white")
    box = ImageChops.difference(rgb, background).getbbox()
    return rgb.crop(box) if box else rgb


def promote_sample(
    sample_id: str,
    iteration: int,
    decision: str,
    note: str,
    manifest: dict,
) -> None:
    iteration_dir = FIGURE_ROOT / f"iteration_{iteration:02d}"
    nested_dir = iteration_dir / sample_id
    source_dir = nested_dir if nested_dir.is_dir() else iteration_dir
    if not source_dir.is_dir():
        raise FileNotFoundError(f"Missing accepted iteration directory: {source_dir}")

    sources = sorted(source_dir.glob(f"{sample_id} *"))
    expected = {
        f"{sample_id} 2d reg.png",
        f"{sample_id} 3d reg.png",
    }
    found = {path.name for path in sources}
    missing = expected - found
    if missing:
        raise FileNotFoundError(
            f"{source_dir} lacks final images required for promotion: {sorted(missing)}"
        )

    # Remove optional diagnostics from a previously accepted candidate before
    # publishing this decision. Also remove baseline registration JPGs when a
    # regional candidate is accepted, because they depict stale coordinates.
    for stale in FIGURE_ROOT.glob(f"{sample_id} regional *.png"):
        _guard_optimized_path(stale)
        stale.unlink()
    if iteration != 3:
        for stale in FIGURE_ROOT.glob(f"{sample_id}_slice_*.jpg"):
            _guard_optimized_path(stale)
            stale.unlink()

    copied = []
    for source in sources:
        if source.suffix.lower() not in {".png", ".jpg", ".jpeg"}:
            continue
        target = FIGURE_ROOT / source.name
        _guard_optimized_path(target)
        shutil.copy2(source, target)
        copied.append(target.name)

    timestamp = datetime.now(timezone.utc).isoformat()
    manifest["samples"][sample_id] = {
        "accepted_iteration": int(iteration),
        "decision": decision,
        "note": note,
        "source": str(source_dir.relative_to(BASE)),
        "files": copied,
        "promoted_at": timestamp,
    }
    save_manifest(manifest)

    row = pd.DataFrame(
        [
            {
                "timestamp": timestamp,
                "sample_id": sample_id,
                "accepted_iteration": int(iteration),
                "decision": decision,
                "note": note,
                "files": "|".join(copied),
            }
        ]
    )
    if PROMOTION_LOG.exists():
        row.to_csv(PROMOTION_LOG, mode="a", header=False, index=False)
    else:
        row.to_csv(PROMOTION_LOG, index=False)


def rebuild_library(manifest: dict) -> None:
    sample_ids = sorted(manifest["samples"])
    if not sample_ids:
        return

    cell_width, cell_height = 760, 500
    title_height = 60
    canvas = Image.new(
        "RGB",
        (2 * cell_width, title_height + len(sample_ids) * cell_height),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    draw.text((20, 20), "Iteration 03 baseline", fill="black")
    draw.text((cell_width + 20, 20), "Current accepted final", fill="black")

    for row_index, sample_id in enumerate(sample_ids):
        current = manifest["samples"][sample_id]
        iteration = int(current["accepted_iteration"])
        paths = [
            FIGURE_ROOT / "iteration_03" / f"{sample_id} 3d reg.png",
            FIGURE_ROOT / f"{sample_id} 3d reg.png",
        ]
        y0 = title_height + row_index * cell_height
        for column_index, path in enumerate(paths):
            image = _trim_white(Image.open(path))
            image.thumbnail((cell_width - 40, cell_height - 55))
            x = column_index * cell_width + (cell_width - image.width) // 2
            y = y0 + 40 + (cell_height - 45 - image.height) // 2
            canvas.paste(image, (x, y))
        draw.text(
            (10, y0 + 8),
            f"{sample_id} (accepted iteration {iteration:02d})",
            fill="black",
        )

    library = FIGURE_ROOT / "alignment_iteration_library.png"
    comparison = FIGURE_ROOT / "accepted_vs_baseline_comparison.png"
    _guard_optimized_path(library)
    canvas.save(library)
    canvas.save(comparison)


def rebuild_iteration_log() -> None:
    history_path = TEST_ROOT / "all_slices_regional_landmark_decisions.csv"
    output_path = FIGURE_ROOT / "iteration_log.csv"
    if not history_path.exists():
        return
    history = pd.read_csv(history_path)
    prior = pd.read_csv(output_path) if output_path.exists() else pd.DataFrame()
    if len(prior) and "iteration" in prior:
        prior = prior[prior["iteration"] < 11]
    rows = []
    for iteration, group in history.groupby("iteration", sort=True):
        decision_counts = group["decision"].value_counts().to_dict()
        accepted = sorted(group.loc[group["decision"].eq("accept"), "sample_id"])
        rows.append(
            {
                "iteration": int(iteration),
                "scope": "all-slice regional landmark refinement",
                "changes": (
                    f"Evaluated {len(group)} slice candidate(s) at landmark weight "
                    f"{group['weight'].iloc[0]:g}; preserved slice-specific baseline geometry."
                ),
                "impact": "; ".join(
                    f"{decision}={count}"
                    for decision, count in sorted(decision_counts.items())
                ),
                "decision": (
                    f"Accepted: {', '.join(accepted)}"
                    if accepted
                    else "No candidates accepted; retained prior accepted/baseline images."
                ),
            }
        )
    rebuilt = pd.concat([prior, pd.DataFrame(rows)], ignore_index=True)
    rebuilt.to_csv(output_path, index=False)


def _selected_iteration(row: pd.Series) -> int:
    decision = str(row.get("decision", "")).lower()
    if decision.startswith("retain"):
        return 3
    for column in ("accepted_iteration", "selected_iteration", "iteration"):
        value = row.get(column)
        if pd.notna(value):
            return int(value)
    raise ValueError(f"No accepted iteration in decision row for {row.get('sample_id')}")


def consolidate_decisions(decisions: pd.DataFrame) -> pd.DataFrame:
    """Return exactly one accepted or retained-baseline row per sample."""
    final_rows: list[dict] = []
    for sample_id, group in decisions.groupby("sample_id", sort=True):
        accepted = group[group["decision"].eq("accept")]
        if len(accepted):
            row = accepted.iloc[-1].to_dict()
            row["note"] = str(row.get("reason", "accepted guarded regional candidate"))
        else:
            row = group.iloc[-1].to_dict()
            history = "; ".join(
                f"iteration {int(item.iteration):02d}: {item.decision} ({item.reason})"
                for item in group.itertuples()
            )
            row.update(
                {
                    "iteration": 3,
                    "weight": 0.0,
                    "decision": "retain_baseline",
                    "accepted_iteration": 3,
                    "try_rescue": False,
                    "reason": history,
                    "note": history,
                }
            )
        final_rows.append(row)
    return pd.DataFrame(final_rows).sort_values("sample_id").reset_index(drop=True)


def promote_from_decisions(decisions_path: Path, only_sample: str | None) -> None:
    raw_decisions = pd.read_csv(decisions_path)
    if "decision" not in raw_decisions and "final_decision" in raw_decisions:
        raw_decisions["decision"] = np.where(
            raw_decisions["final_decision"].astype(str).str.startswith("accept"),
            "accept",
            "retain_baseline",
        )
        raw_decisions["iteration"] = raw_decisions["accepted_iteration"].fillna(3)
        raw_decisions["weight"] = raw_decisions.get("accepted_weight", 0).fillna(0)
        raw_decisions["try_rescue"] = False
    all_decisions = consolidate_decisions(raw_decisions)
    decisions = all_decisions
    if only_sample:
        decisions = decisions[decisions["sample_id"].eq(only_sample)]
    if decisions.empty:
        raise ValueError("No decision rows selected for promotion")

    manifest = load_manifest()
    for _, row in decisions.sort_values("sample_id").iterrows():
        decision = str(row.get("decision", ""))
        if decision.startswith("reject") or decision == "needs_manual":
            continue
        note = str(row.get("note", row.get("reason", "")))
        promote_sample(
            str(row["sample_id"]),
            _selected_iteration(row),
            decision,
            note,
            manifest,
        )
        rebuild_library(manifest)

    all_decisions.to_csv(FIGURE_ROOT / "final_slice_assessment.csv", index=False)
    promotion_summary = pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                **entry,
            }
            for sample_id, entry in sorted(manifest["samples"].items())
        ]
    )
    promotion_summary.to_csv(FIGURE_ROOT / "current_final_manifest.csv", index=False)
    rebuild_iteration_log()
    rebuild_library(manifest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    parser.add_argument("--sample", help="Promote only one finalized sample")
    args = parser.parse_args()
    promote_from_decisions(args.decisions, args.sample)


if __name__ == "__main__":
    main()
