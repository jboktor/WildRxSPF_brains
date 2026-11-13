#!/usr/bin/env python3
"""Register the November 2026 raw-scalar 2DG export to Allen CCF.

Inputs:
  data/input/2DG_data/Template/TemplateNNN.tif
  data/input/2DG_data/WildRvSPF Intruder/WildRvSPF_IntruderNNN.tif

The first pair is an explicit sentinel (Template001 is zero; Intruder001 is
uniform y=142, approximately t=0) and is excluded. Pairs 002..066 are aligned
after removing the one-pixel border from the 277x377 Intruder images. Template
indices are converted with 2DG=y/1.494; Intruder indices are converted directly
with t=(y/18.409)-7.762.

This is a versioned branch. It imports the tested registration machinery from
register_2dg_to_allen_ccf.py but writes only to raw_scalar_2026-11-11 paths.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

import register_2dg_to_allen_ccf as core


RAW_ANALYSIS_DATE = "2026-11-11"
WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
TEMPLATE_ROOT = WKDIR / "data/input/2DG_data/Template"
INTRUDER_ROOT = WKDIR / "data/input/2DG_data/WildRvSPF Intruder"
DEFAULT_OUTPUT_ROOT = (
    WKDIR
    / "data/interim/2DG_CCF_Integration"
    / f"raw_scalar_{RAW_ANALYSIS_DATE}"
)
DEFAULT_FIGURE_ROOT = (
    WKDIR / "figures/2DG_CCF_Integration" / f"raw_scalar_{RAW_ANALYSIS_DATE}"
)
DEFAULT_TABLE_ROOT = (
    WKDIR / "tables/2DG_CCF_Integration" / f"raw_scalar_{RAW_ANALYSIS_DATE}"
)

FIRST_SOURCE_NUMBER = 2
LAST_SOURCE_NUMBER = 66
FIRST_BREGMA_MM = 2.16
SPACING_MM = 0.14
TEMPLATE_SCALE = 1.494
T_SCALE = 18.409
T_OFFSET = 7.762

# Ensure imported helper filenames use the raw branch date.
core.ANALYSIS_DATE = RAW_ANALYSIS_DATE


def log(message: str) -> None:
    print(f"[2DG-RAW-CCF] {message}", flush=True)


def load_index_tiff(path: Path) -> np.ndarray:
    image = Image.open(path)
    array = np.asarray(image, dtype=np.uint8)
    if array.ndim != 2:
        raise ValueError(f"Expected one indexed channel, got {array.shape}: {path}")
    return array


def source_paths(source_number: int) -> tuple[Path, Path]:
    return (
        TEMPLATE_ROOT / f"Template{source_number:03d}.tif",
        INTRUDER_ROOT / f"WildRvSPF_Intruder{source_number:03d}.tif",
    )


def tissue_mask(template_y: np.ndarray) -> np.ndarray:
    candidate = template_y > 0
    candidate = ndimage.binary_closing(candidate, iterations=2)
    candidate = ndimage.binary_fill_holes(candidate)
    labels, count = ndimage.label(candidate)
    if not count:
        return candidate
    sizes = np.bincount(labels.ravel())
    keep = np.argsort(sizes[1:])[-4:] + 1
    result = np.isin(labels, keep)
    result = ndimage.binary_closing(result, iterations=2)
    return ndimage.binary_fill_holes(result)


def load_pair(source_number: int) -> dict:
    template_path, intruder_path = source_paths(source_number)
    template_y = load_index_tiff(template_path)
    intruder_full = load_index_tiff(intruder_path)
    if intruder_full.shape != (
        template_y.shape[0] + 2,
        template_y.shape[1] + 2,
    ):
        raise ValueError(
            f"Expected a one-pixel Intruder border for source {source_number}: "
            f"template={template_y.shape}, intruder={intruder_full.shape}"
        )
    intruder_y = intruder_full[1:-1, 1:-1]
    if intruder_y.shape != template_y.shape:
        raise ValueError(f"Core grids do not match for source {source_number}.")
    tissue = tissue_mask(template_y)
    t_value = intruder_y.astype(np.float32) / T_SCALE - T_OFFSET
    overlay = tissue.copy()
    branch = np.zeros(template_y.shape, dtype=np.int8)
    branch[tissue & (t_value < 0)] = -1
    branch[tissue & (t_value > 0)] = 1
    display_magnitude = np.full(template_y.shape, np.nan, dtype=np.float32)
    display_magnitude[tissue] = np.abs(t_value[tissue]) / max(
        abs(-T_OFFSET),
        abs(255 / T_SCALE - T_OFFSET),
    )
    intruder_rgb = np.repeat(intruder_y[:, :, None], 3, axis=2)
    return {
        "template_path": template_path,
        "intruder_path": intruder_path,
        "template_y": template_y,
        "template_2dg": template_y.astype(np.float32) / TEMPLATE_SCALE,
        "tissue": tissue.astype(np.uint8),
        "intruder_rgb": intruder_rgb,
        "overlay": overlay.astype(np.uint8),
        "overlay_branch": branch,
        "display_magnitude": display_magnitude,
        "tiff_y": intruder_y,
        "t_value": t_value,
    }


def build_manifest() -> pd.DataFrame:
    rows = []
    for position, source_number in enumerate(
        range(FIRST_SOURCE_NUMBER, LAST_SOURCE_NUMBER + 1)
    ):
        pair = load_pair(source_number)
        bregma = FIRST_BREGMA_MM - position * SPACING_MM
        rows.append(
            {
                "slice_position": position,
                "source_number": source_number,
                "slice_id": f"2dg_raw_{source_number:03d}",
                "intruder_path": str(pair["intruder_path"]),
                "baseline_path": str(pair["template_path"]),
                "template_path": str(pair["template_path"]),
                "height_px": pair["template_y"].shape[0],
                "width_px": pair["template_y"].shape[1],
                "nominal_bregma_mm": bregma,
                "nominal_allen_index": core.nominal_allen_index(bregma),
                "intruder_gray_fraction": 1.0,
                "baseline_gray_fraction": 1.0,
                "hidden_by_both_fraction": 0.0,
                "shared_gray_exact": False,
                "email_expected_slice_count": core.EXPECTED_EMAIL_SLICE_COUNT,
                "observed_slice_count": LAST_SOURCE_NUMBER - FIRST_SOURCE_NUMBER + 1,
                "slice_count_matches_email": False,
                "source_image_left_is_anatomical_right": True,
                "source_flipped_for_ccf": core.FLIP_SOURCE_LEFT_RIGHT,
                "t_value_status": "direct_scalar_index_and_email_formula",
                "input_export_version": "raw_scalar_2026-11-11",
                "intruder_border_crop": "remove one pixel on every edge",
                "excluded_source_001": True,
                "n_tissue_pixels": int(pair["tissue"].sum()),
                "template_y_min": int(pair["template_y"].min()),
                "template_y_max": int(pair["template_y"].max()),
                "intruder_y_min": int(pair["tiff_y"].min()),
                "intruder_y_max": int(pair["tiff_y"].max()),
                "t_min": float(pair["t_value"][pair["tissue"].astype(bool)].min()),
                "t_max": float(pair["t_value"][pair["tissue"].astype(bool)].max()),
            }
        )
    return pd.DataFrame(rows)


def prepare(
    output_root: Path,
    figure_root: Path,
) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "prepared").mkdir(parents=True, exist_ok=True)
    (output_root / "work").mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)
    manifest_file = core.manifest_path(output_root)
    if manifest_file.exists():
        manifest = pd.read_csv(manifest_file)
        log(f"Using existing manifest without modification: {manifest_file}")
    else:
        manifest = build_manifest()
        core.write_csv_new(manifest_file, manifest)
        log(f"Wrote {len(manifest)}-slice raw manifest: {manifest_file}")

    for row in manifest.itertuples(index=False):
        destination = core.prepared_path(output_root, int(row.slice_position))
        if destination.exists():
            continue
        pair = load_pair(int(row.source_number))
        arrays = {
            key: pair[key]
            for key in (
                "template_y",
                "template_2dg",
                "tissue",
                "intruder_rgb",
                "overlay",
                "overlay_branch",
                "display_magnitude",
                "tiff_y",
                "t_value",
            )
        }
        if core.FLIP_SOURCE_LEFT_RIGHT:
            arrays = {
                name: np.flip(value, axis=1).copy()
                for name, value in arrays.items()
            }
        np.savez_compressed(destination, **arrays)

    provenance_path = (
        output_root / f"2dg_raw_preparation_qc_{RAW_ANALYSIS_DATE}.json"
    )
    if not provenance_path.exists():
        sentinel_template, sentinel_intruder = source_paths(1)
        sentinel_template_y = load_index_tiff(sentinel_template)
        sentinel_intruder_y = load_index_tiff(sentinel_intruder)[1:-1, 1:-1]
        core.write_json_new(
            provenance_path,
            {
                "analysis_date": RAW_ANALYSIS_DATE,
                "input_export_version": "raw_scalar_2026-11-11",
                "template_root": TEMPLATE_ROOT,
                "intruder_root": INTRUDER_ROOT,
                "source_pairs_used": "002..066",
                "source_001_exclusion": {
                    "template_unique_values": np.unique(
                        sentinel_template_y
                    ).tolist(),
                    "intruder_unique_values": np.unique(
                        sentinel_intruder_y
                    ).tolist(),
                    "reason": "blank template and uniform near-zero t sentinel",
                },
                "intruder_crop": "[1:-1, 1:-1] to remove one-pixel border",
                "prepared_shape": [275, 375],
                "nominal_first_bregma_mm": FIRST_BREGMA_MM,
                "spacing_mm": SPACING_MM,
                "nominal_last_bregma_mm": float(
                    manifest.nominal_bregma_mm.iloc[-1]
                ),
                "template_formula": "2DG = y / 1.494",
                "t_formula": "t = (y / 18.409) - 7.762",
                "source_orientation": "image left is anatomical right (email)",
                "ccf_orientation_action": "horizontal flip before registration",
            },
        )
    core.write_preparation_figure(output_root, figure_root, manifest)
    log(f"Prepared {len(manifest)} raw-scalar slices.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--figure-root", type=Path, default=DEFAULT_FIGURE_ROOT)
    parser.add_argument("--table-root", type=Path, default=DEFAULT_TABLE_ROOT)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("prepare")
    register_parser = subparsers.add_parser("register-slice")
    register_parser.add_argument("--slice-position", type=int, required=True)
    subparsers.add_parser("aggregate")
    subparsers.add_parser("status")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare(args.output_root, args.figure_root)
    elif args.command == "register-slice":
        core.register_slice(args.output_root, args.slice_position)
    elif args.command == "aggregate":
        core.aggregate(args.output_root, args.figure_root, args.table_root)
    elif args.command == "status":
        core.status(args.output_root)


if __name__ == "__main__":
    main()
