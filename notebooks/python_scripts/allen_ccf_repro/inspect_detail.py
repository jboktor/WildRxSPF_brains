"""Print the candidate ladder and per-region metrics from a slice's detail.json.

    python -m allen_ccf_repro.inspect_detail roi1_run2
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "allen_ccf_repro"

from . import config  # noqa: E402

NAN = float("nan")


def show(sample_id: str) -> None:
    detail = json.loads((config.WORK_ROOT / sample_id / "detail.json").read_text())
    record = detail["record"]
    candidates = detail["candidates"]

    print(f"=== {sample_id}: {record['status']}")
    print(
        f"    slice {record['slice_number']} baseline plane {record['baseline_allen_slice']} "
        f"-> plane {record['allen_slice_num']} gradient {record['ap_gradient']:+.3f} "
        f"medial_extra {record['medial_extra']} ({record.get('chosen_geometry')})"
    )
    for trial in detail["ap_search"].get("geometry_trials", []):
        print(
            f"    trial {trial['tag']:18s} z={trial['slice_index']} g={trial['gradient']:+.3f} "
            f"m={trial['medial_extra']:2d} gain={trial['mean_gated_gain']:6.2f} "
            f":: {trial['reason'][:110]}"
        )

    regions = sorted(
        key[len("baseline_") : -len("_inside_pct")]
        for key in record
        if key.startswith("baseline_") and key.endswith("_inside_pct")
    )
    print(f"{'region':20s} {'base':>6s}" + "".join(f"{item['label'][-9:]:>10s}" for item in candidates))
    for region in regions:
        row = f"{region:20s} {record.get('baseline_' + region + '_inside_pct', NAN):6.1f}"
        for item in candidates:
            row += f"{item['metrics'].get(region + '_inside_pct', NAN):10.1f}"
        print(row)

    print()
    for item in candidates:
        measures = item["metrics"]
        print(
            f"    {item['label']:16s} accepted={item['accepted']!s:5s} "
            f"tissue={measures['tissue_containment_pct']:.2f} "
            f"gain={item['gains']['mean_gated_gain']:6.2f} "
            f"shift={measures.get('geometry_shift_25um', NAN):5.2f} "
            f"anchor={measures.get('anchor_median_displacement_25um', NAN):5.2f} "
            f"jmin={measures.get('jacobian_min', NAN):.2f}"
        )
        print(f"        {item['reason'][:220]}")
    for probe in detail.get("sensitivity", []):
        print(
            f"    mask {probe['delta_px']:+d}px accepted={probe['accepted']} :: {probe['reason'][:110]}"
        )


def main() -> None:
    names = sys.argv[1:] or [
        path.parent.name for path in sorted(config.WORK_ROOT.glob("*/detail.json"))
    ]
    for name in names:
        if (config.WORK_ROOT / name / "detail.json").exists():
            show(name)
            print()


if __name__ == "__main__":
    main()
