#!/usr/bin/env python3
"""Warp sibling 2DG maps through accepted Template→CCF registrations.

Reuses per-slice pixel tables and transforms from the scalar Template run
(2026-11-11). Does not re-fit Elastix. Sibling RGB overlays are decoded with
the validated 254-color LUT when present.

Run with spatialomics Python:

  python warp_2dg_sibling_maps_to_ccf.py all
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
SOURCE_ANALYSIS_DATE = "2026-11-23c"
ANALYSIS_DATE = "2026-11-23c"
T_SIGNIFICANCE = 1.725

INPUT_ROOT = WKDIR / "data/input/2DG_data"
SOURCE_ROOT = WKDIR / "data/interim/2DG_CCF_Integration" / SOURCE_ANALYSIS_DATE
OUTPUT_ROOT = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
TABLE_ROOT = WKDIR / "tables/2DG_CCF_Integration"
FIGURE_ROOT = WKDIR / "figures/2DG_CCF_Integration"
LUT_PATH = TABLE_ROOT / "2DG_tiff_rgb_to_t_lut_validated_2026-11-07.csv"
ANNOTATION_NAMES = TABLE_ROOT / f"Allen_CCF_annotation_names_hierarchy_{SOURCE_ANALYSIS_DATE}.csv"

REGISTRATION_SUMMARY = (
    SOURCE_ROOT / f"2dg_registration_summary_strict_{SOURCE_ANALYSIS_DATE}.csv"
)
INTRUDER_REGION_SUMMARY = (
    TABLE_ROOT / f"2DG_CCF_region_summary_strict_{SOURCE_ANALYSIS_DATE}.csv"
)

MAP_FAMILIES = {
    "baseline_wildr_vs_spf": {
        "root": INPUT_ROOT / "WilderR vs SPF _Baseline",
        "kind": "legacy_rgb",
        "pattern": "template42{source:03d}.tif",
        "t_status": "decoded_from_complete_254_color_palette_and_email_formula",
    },
    "wildr_intruder_vs_baseline": {
        "root": INPUT_ROOT / "WildR_IntrudervBaseline",
        "kind": "rgb",
        "pattern": "WildRIntruderVBaseline{source:03d}.tif",
        "t_status": "decoded_from_complete_254_color_palette_and_email_formula",
    },
    "spf_intruder_vs_baseline": {
        "root": INPUT_ROOT / "SPF_IntrudervBaseline",
        "kind": "rgb",
        "pattern": "SPF_IntrudervBaseline{source:03d}.tif",
        "t_status": "decoded_from_complete_254_color_palette_and_email_formula",
    },
}

ROI_SEEDS = [
    {
        "roi": "L_LS",
        "source_number": 48,
        "overlay_path": INPUT_ROOT / "L_LS_WildR_Intruder" / "L_LS_WildR_Intruder048.tif",
        "target_name_regex": r"(?i)lateral septal|septal nucleus",
    },
    {
        "roi": "L_NAc",
        "source_number": 59,
        "overlay_path": INPUT_ROOT / "L_NAc_WildR_Intruder" / "L_NAc_WildR_Intruder059.tif",
        "target_name_regex": r"(?i)nucleus accumbens|accumbens",
    },
    {
        "roi": "R_BMA",
        "source_number": 44,
        "overlay_path": INPUT_ROOT / "R_BMA_WildR_Intruder" / "R_BMA_WildR_Intruder044.tif",
        "target_name_regex": r"(?i)basomedial|basomedial amygdalar",
    },
    {
        "roi": "R_Cg",
        "source_number": 50,
        "overlay_path": INPUT_ROOT / "R_Cg_WildR_Intruder" / "R_Cg_WildR_Intruder050.tif",
        "target_name_regex": r"(?i)cingulate|anterior cingulate",
    },
    {
        "roi": "R_PVN",
        "source_number": 45,
        "overlay_path": INPUT_ROOT / "R_PVN_WildR_Intruder" / "R_PVN_WildR_Intruder045.tif",
        "target_name_regex": r"(?i)paraventricular.*(hypothal|nucleus of the hypothal)|periventricular",
    },
    {
        "roi": "R_vHPC",
        "source_number": 26,
        "overlay_path": INPUT_ROOT / "R_vHPC_WildR_Intruder" / "R_vHPC_WildR_Intruder026.tif",
        "target_name_regex": r"(?i)hippocamp|subiculum|CA[1-3]|dentate",
    },
]


def log(message: str) -> None:
    print(f"[2DG-warp] {message}", flush=True)


def write_csv_new(path: Path, frame: pd.DataFrame, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    frame.to_csv(path, index=False, compression=compression)


def write_json_new(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def load_lut() -> pd.DataFrame:
    lut = pd.read_csv(LUT_PATH)
    required = {"overlay_r", "overlay_g", "overlay_b", "tiff_y", "t_value"}
    missing = required - set(lut.columns)
    if missing:
        raise ValueError(f"LUT missing columns: {sorted(missing)}")
    return lut


def load_rgb(path: Path) -> np.ndarray:
    array = np.asarray(Image.open(path).convert("RGB"), dtype=np.uint8)
    if array.ndim != 3 or array.shape[2] != 3:
        raise ValueError(f"Expected RGB TIFF at {path}, got {array.shape}")
    return array


def grayscale_mask(rgb: np.ndarray) -> np.ndarray:
    return (rgb[..., 0] == rgb[..., 1]) & (rgb[..., 1] == rgb[..., 2])


def accepted_slices() -> pd.DataFrame:
    summary = pd.read_csv(REGISTRATION_SUMMARY)
    accepted = summary.loc[summary.accepted_final.astype(bool)].copy()
    if accepted.empty:
        raise RuntimeError("No accepted slices in registration summary.")
    return accepted.sort_values("slice_position")


def pixel_table_path(slice_position: int) -> Path:
    return SOURCE_ROOT / "work" / f"slice_{slice_position:03d}" / "pixels_ccf.csv.gz"


def sibling_path(family: str, source_number: int) -> Path | None:
    spec = MAP_FAMILIES[family]
    path = spec["root"] / spec["pattern"].format(source=source_number)
    return path if path.exists() else None


def decode_overlay_at_pixels(
    rgb: np.ndarray, pixels: pd.DataFrame, lut: pd.DataFrame, t_status: str
) -> pd.DataFrame:
    ys = pixels.source_y.to_numpy(dtype=int)
    xs = pixels.source_x.to_numpy(dtype=int)
    if (
        ys.min() < 0
        or xs.min() < 0
        or ys.max() >= rgb.shape[0]
        or xs.max() >= rgb.shape[1]
    ):
        raise ValueError("Source coordinates fall outside sibling image bounds.")
    sample = rgb[ys, xs]
    frame = pixels.drop(
        columns=["overlay_r", "overlay_g", "overlay_b", "tiff_y", "t_value", "t_value_status"],
        errors="ignore",
    ).copy()
    frame["overlay_r"] = sample[:, 0]
    frame["overlay_g"] = sample[:, 1]
    frame["overlay_b"] = sample[:, 2]
    gray = (
        (frame.overlay_r == frame.overlay_g) & (frame.overlay_g == frame.overlay_b)
    )
    decoded = frame.merge(
        lut[
            ["overlay_r", "overlay_g", "overlay_b", "tiff_y", "t_value"]
        ],
        on=["overlay_r", "overlay_g", "overlay_b"],
        how="left",
        validate="many_to_one",
    )
    decoded["overlay_present"] = (~gray) & decoded.t_value.notna()
    decoded.loc[~decoded.overlay_present, "t_value"] = np.nan
    decoded.loc[~decoded.overlay_present, "tiff_y"] = np.nan
    decoded["t_value_status"] = np.where(
        decoded.overlay_present, t_status, "not_an_overlay_pixel"
    )
    decoded["significant_overlay"] = decoded.overlay_present & (
        decoded.t_value.abs() >= T_SIGNIFICANCE
    )
    return decoded


def warp_family(family: str, lut: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    accepted = accepted_slices()
    overlay_frames = []
    region_frames = []
    missing = []
    for row in accepted.itertuples(index=False):
        path = sibling_path(family, int(row.source_number))
        if path is None:
            missing.append(int(row.source_number))
            continue
        pixels = pd.read_csv(pixel_table_path(int(row.slice_position)))
        inside = pixels.loc[
            pixels.inside_allen_tissue.astype(bool) & (pixels.annotation != 0)
        ].copy()
        if inside.empty:
            continue
        rgb = load_rgb(path)
        decoded = decode_overlay_at_pixels(
            rgb, inside, lut, MAP_FAMILIES[family]["t_status"]
        )
        decoded["map_family"] = family
        significant = decoded.loc[decoded.significant_overlay].copy()
        if not significant.empty:
            overlay_frames.append(
                significant[
                    [
                        "map_family",
                        "slice_position",
                        "source_number",
                        "ccfx",
                        "ccfy",
                        "ccfz",
                        "annotation",
                        "annotation_name",
                        "overlay_r",
                        "overlay_g",
                        "overlay_b",
                        "tiff_y",
                        "t_value",
                        "t_value_status",
                    ]
                ]
            )
        region = (
            decoded.groupby(["annotation", "annotation_name"], as_index=False)
            .agg(
                n_template_pixels=("annotation", "size"),
                n_overlay_pixels=("significant_overlay", "sum"),
                mean_t_value=("t_value", "mean"),
                median_t_value=("t_value", "median"),
                mean_positive_t=(
                    "t_value",
                    lambda s: float(np.nanmean(s[s > 0])) if np.any(s > 0) else np.nan,
                ),
                mean_absolute_t=("t_value", lambda s: float(np.nanmean(np.abs(s)))),
            )
        )
        region["map_family"] = family
        region["slice_position"] = int(row.slice_position)
        region["source_number"] = int(row.source_number)
        region["selected_allen_index"] = int(row.selected_allen_index)
        region_frames.append(region)
    if not region_frames:
        raise RuntimeError(f"No warped regions produced for {family}")
    regions_by_slice = pd.concat(region_frames, ignore_index=True)
    region_summary = (
        regions_by_slice.groupby(["map_family", "annotation", "annotation_name"], as_index=False)
        .agg(
            n_slices=("slice_position", "nunique"),
            n_template_pixels=("n_template_pixels", "sum"),
            n_overlay_pixels=("n_overlay_pixels", "sum"),
            mean_t_value=("mean_t_value", "mean"),
            median_t_value=("median_t_value", "median"),
            mean_positive_t=("mean_positive_t", "mean"),
            mean_absolute_t=("mean_absolute_t", "mean"),
        )
    )
    overlays = (
        pd.concat(overlay_frames, ignore_index=True)
        if overlay_frames
        else pd.DataFrame(
            columns=[
                "map_family",
                "slice_position",
                "source_number",
                "ccfx",
                "ccfy",
                "ccfz",
                "annotation",
                "annotation_name",
                "overlay_r",
                "overlay_g",
                "overlay_b",
                "tiff_y",
                "t_value",
                "t_value_status",
            ]
        )
    )
    log(
        f"{family}: slices_used={regions_by_slice.slice_position.nunique()} "
        f"missing_sources={missing} overlay_points={len(overlays)} regions={len(region_summary)}"
    )
    return overlays, region_summary


def warp_roi_seeds() -> pd.DataFrame:
    accepted = accepted_slices().set_index("source_number")
    names = pd.read_csv(ANNOTATION_NAMES)
    rows = []
    for spec in ROI_SEEDS:
        source_number = int(spec["source_number"])
        if source_number not in accepted.index:
            rows.append(
                {
                    "roi": spec["roi"],
                    "source_number": source_number,
                    "slice_accepted": False,
                    "n_overlay_pixels": 0,
                    "n_inside_allen": 0,
                    "n_target_pixels": 0,
                    "iou_target": np.nan,
                    "decision_note": "seed_slice_not_in_accepted_registration",
                }
            )
            continue
        meta = accepted.loc[source_number]
        slice_position = int(meta.slice_position)
        pixels = pd.read_csv(pixel_table_path(slice_position))
        inside = pixels.loc[
            pixels.inside_allen_tissue.astype(bool) & (pixels.annotation != 0)
        ].copy()
        rgb = load_rgb(Path(spec["overlay_path"]))
        ys = inside.source_y.to_numpy(dtype=int)
        xs = inside.source_x.to_numpy(dtype=int)
        sample = rgb[ys, xs]
        gray = (sample[:, 0] == sample[:, 1]) & (sample[:, 1] == sample[:, 2])
        footprint = inside.loc[~gray].copy()
        target = names.loc[
            names.annotation_name.str.contains(spec["target_name_regex"], regex=True, na=False),
            "annotation",
        ].astype(int)
        target_set = set(target.tolist())
        n_target = int(footprint.annotation.isin(target_set).sum()) if len(footprint) else 0
        n_inside = int(len(footprint))
        # IoU vs all inside pixels labeled as target on this slice.
        slice_target = inside.loc[inside.annotation.isin(target_set)]
        union = len(
            set(zip(footprint.source_x, footprint.source_y))
            | set(zip(slice_target.source_x, slice_target.source_y))
        )
        inter = n_target
        iou = float(inter / union) if union else np.nan
        rows.append(
            {
                "roi": spec["roi"],
                "source_number": source_number,
                "slice_position": slice_position,
                "selected_allen_index": int(meta.selected_allen_index),
                "slice_accepted": True,
                "n_overlay_pixels": n_inside,
                "n_inside_allen": n_inside,
                "n_target_pixels": n_target,
                "n_slice_target_pixels": int(len(slice_target)),
                "iou_target": iou,
                "target_name_regex": spec["target_name_regex"],
                "n_target_annotations": int(len(target_set)),
                "decision_note": "ok",
            }
        )
    return pd.DataFrame(rows)


def copy_intruder_reference() -> None:
    """Point hypothesis code at Intruder tables with a dated sidecar index."""
    path = OUTPUT_ROOT / f"2dg_map_family_index_{ANALYSIS_DATE}.json"
    if path.exists():
        return
    payload = {
        "analysis_date": ANALYSIS_DATE,
        "source_registration": SOURCE_ANALYSIS_DATE,
        "intruder_region_summary": str(INTRUDER_REGION_SUMMARY),
        "intruder_overlay_points": str(
            SOURCE_ROOT
            / f"2dg_overlay_points_ccf_scalar-t_{SOURCE_ANALYSIS_DATE}.csv.gz"
        ),
        "registration_summary": str(REGISTRATION_SUMMARY),
        "map_families": list(MAP_FAMILIES) + ["intruder_wildr_vs_spf"],
    }
    write_json_new(path, payload)


def run_all() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    lut = load_lut()
    copy_intruder_reference()

    all_overlays = []
    all_regions = []
    for family in MAP_FAMILIES:
        overlays, regions = warp_family(family, lut)
        all_overlays.append(overlays)
        all_regions.append(regions)
        write_csv_new(
            OUTPUT_ROOT / f"2dg_overlay_points_ccf_{family}_{ANALYSIS_DATE}.csv.gz",
            overlays,
            compression="gzip",
        )
        write_csv_new(
            TABLE_ROOT / f"2DG_CCF_region_summary_{family}_{ANALYSIS_DATE}.csv",
            regions,
        )

    combined_overlays = pd.concat(all_overlays, ignore_index=True)
    combined_regions = pd.concat(all_regions, ignore_index=True)
    write_csv_new(
        OUTPUT_ROOT / f"2dg_overlay_points_ccf_siblings_{ANALYSIS_DATE}.csv.gz",
        combined_overlays,
        compression="gzip",
    )
    write_csv_new(
        TABLE_ROOT / f"2DG_CCF_region_summary_siblings_{ANALYSIS_DATE}.csv",
        combined_regions,
    )

    roi = warp_roi_seeds()
    write_csv_new(TABLE_ROOT / f"2DG_ROI_seed_ccf_overlap_{ANALYSIS_DATE}.csv", roi)
    write_roi_figure(roi)
    write_json_new(
        OUTPUT_ROOT / f"2dg_sibling_warp_provenance_{ANALYSIS_DATE}.json",
        {
            "analysis_date": ANALYSIS_DATE,
            "source_registration": SOURCE_ANALYSIS_DATE,
            "lut_path": str(LUT_PATH),
            "n_overlay_points": int(len(combined_overlays)),
            "n_regions": int(len(combined_regions)),
            "families": {
                family: int((combined_regions.map_family == family).sum())
                for family in MAP_FAMILIES
            },
            "roi_seeds": roi.to_dict(orient="records"),
        },
    )
    log(f"Wrote sibling warp products under {OUTPUT_ROOT} and {TABLE_ROOT}")


def write_roi_figure(roi: pd.DataFrame) -> None:
    path = FIGURE_ROOT / f"2DG_ROI_seed_ccf_overlap_{ANALYSIS_DATE}.png"
    if path.exists():
        return
    plot = roi.copy()
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    colors = np.where(plot.slice_accepted, "#FF6C0C", "#999999")
    ax.barh(plot.roi, plot.iou_target.fillna(0), color=colors)
    ax.set_xlabel("IoU of ROI overlay footprint vs target Allen parcels (seed slice)")
    ax.set_title("H5 spatial QC: ROI seed footprints in expected Allen parcels")
    ax.set_xlim(0, max(0.2, float(np.nanmax(plot.iou_target.values)) * 1.25))
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["all", "roi-only"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "all":
        run_all()
    elif args.command == "roi-only":
        OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
        roi = warp_roi_seeds()
        write_csv_new(TABLE_ROOT / f"2DG_ROI_seed_ccf_overlap_{ANALYSIS_DATE}.csv", roi)
        write_roi_figure(roi)


if __name__ == "__main__":
    main()
