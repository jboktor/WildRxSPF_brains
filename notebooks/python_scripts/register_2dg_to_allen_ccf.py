#!/usr/bin/env python3
"""Register the scalar 2DG Template + WildRvSPF Intruder stacks to Allen CCFv3.

Primary inputs (preferred over the older RGB overlay exports):

* ``Template/Template###.tif`` — identity-palette grayscale anatomy (y).
  Background is 0; tissue is y > 0. Convert with ``2DG = y / 1.494``.
* ``WildRvSPF Intruder/WildRvSPF_Intruder###.tif`` — identity-palette
  grayscale statistical map (y), one pixel larger on each side with border
  value 142. Crop to the Template FOV, then ``t = (y / 18.409) - 7.762``.
* ``2DG_Template_Allen_manual_matches_2026-11-23b_final.csv`` — manual
  bregma anchors per Template index. AP search is restricted to
  ``AP_SEARCH_OFFSETS`` around each anchor (±0.1 mm).

The older RGB Baseline/Intruder overlays remain on disk as secondary archives
and are not used by the primary pipeline. Existing dated artifacts are never
overwritten.

Run with the spatialomics Python binary.  Typical sequence:

  python register_2dg_to_allen_ccf.py prepare
  python register_2dg_to_allen_ccf.py register-slice --slice-position 0
  python register_2dg_to_allen_ccf.py aggregate

The independent slice registrations are intended for the accompanying SLURM
array script.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nrrd
import numpy as np
import pandas as pd
import SimpleITK as sitk
from PIL import Image
from scipy import ndimage
from skimage.segmentation import find_boundaries


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
ANALYSIS_DATE = "2026-11-23c"
EFFECT_ANALYSIS_DATE = "2026-11-23c"
INPUT_ROOT = WKDIR / "data/input/2DG_data"
TEMPLATE_ROOT = INPUT_ROOT / "Template"
INTRUDER_SCALAR_ROOT = INPUT_ROOT / "WildRvSPF Intruder"
MANUAL_MATCHES_PATH = (
    INPUT_ROOT / "2DG_Template_Allen_manual_matches_2026-11-23b_final.csv"
)
# Legacy RGB exports retained for inventory/cross-checks only.
LEGACY_INTRUDER_RGB_ROOT = INPUT_ROOT / "WildR vs SPF _Intruder"
LEGACY_BASELINE_RGB_ROOT = INPUT_ROOT / "WilderR vs SPF _Baseline"
ATLAS_ROOT = WKDIR / "data/input/allen_registration_ref"
DEFAULT_OUTPUT_ROOT = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
FIGURE_ROOT = WKDIR / "figures/2DG_CCF_Integration"
TABLE_ROOT = WKDIR / "tables/2DG_CCF_Integration"

NISSL_PATH = ATLAS_ROOT / "ara_nissl_25.nrrd"
ANNOTATION_PATH = ATLAS_ROOT / "annotation_25.nrrd"
ANNOTATION_NAMES_PATH = ATLAS_ROOT / "allen_name_to_annots.pkl"

CCF_VOXEL_UM = 25.0
CCF_BREGMA_ZERO_INDEX = 216
NOMINAL_FIRST_BREGMA_MM = 2.16
NOMINAL_SECTION_SPACING_MM = 0.14
# Email quoted ~60 sections; Template001 is blank, Template002 aligns with the
# historical first tissue section (old template42002).
FIRST_TISSUE_SOURCE_NUMBER = 2
EXPECTED_EMAIL_SLICE_COUNT = 60
TEMPLATE_BACKGROUND = 0
INTRUDER_BACKGROUND = 142
TEMPLATE_2DG_SCALE = 1.494
T_Y_SCALE = 18.409
T_Y_OFFSET = 7.762
T_SIGNIFICANCE = 1.725
T_VALUE_STATUS = "scalar_tiff_y_and_email_formula"

# The email states that the image's left side is the brain's right side.  Allen
# CCF z=0 is the anatomical left side, so source arrays are mirrored before fit.
FLIP_SOURCE_LEFT_RIGHT = True

MIN_DICE = 0.45
MIN_INSIDE_FRACTION = 0.80
MIN_JACOBIAN = 0.15
MAX_JACOBIAN = 6.0
# Manual bregma anchors are trusted; allow only ±0.1 mm (±4 voxels at 25 µm).
AP_SEARCH_OFFSETS = (-4, -2, 0, 2, 4)


@dataclass
class Atlas:
    nissl: np.ndarray
    annotation: np.ndarray


_ATLAS: Atlas | None = None


def log(message: str) -> None:
    print(f"[2DG-CCF] {message}", flush=True)


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def write_json_new(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    with path.open("w") as handle:
        json.dump(json_safe(payload), handle, indent=2, sort_keys=True)


def write_csv_new(path: Path, frame: pd.DataFrame, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    frame.to_csv(path, index=False, compression=compression)


def source_files() -> list[tuple[Path, Path]]:
    templates = sorted(TEMPLATE_ROOT.glob("Template*.tif"))
    intruders = sorted(INTRUDER_SCALAR_ROOT.glob("WildRvSPF_Intruder*.tif"))
    if not templates or not intruders:
        raise FileNotFoundError(
            "Scalar Template / WildRvSPF Intruder TIFF stacks were not found."
        )
    by_number: dict[int, Path] = {}
    for path in intruders:
        number = int(path.stem.replace("WildRvSPF_Intruder", ""))
        by_number[number] = path
    pairs = []
    for path in templates:
        number = int(path.stem.replace("Template", ""))
        if number not in by_number:
            raise ValueError(f"Missing Intruder pair for {path.name}")
        pairs.append((path, by_number[number]))
    if len(pairs) != len(intruders):
        raise ValueError("Template and Intruder TIFF counts do not match one-to-one.")
    return pairs


def load_grayscale_y(path: Path) -> np.ndarray:
    image = Image.open(path)
    array = np.asarray(image)
    if array.ndim != 2:
        raise ValueError(f"Expected grayscale TIFF, found shape {array.shape}: {path}")
    if image.mode not in {"P", "L", "I;16", "I"} and array.dtype != np.uint8:
        raise ValueError(f"Unexpected TIFF mode {image.mode} at {path}")
    return array.astype(np.uint8, copy=False)


def crop_intruder_to_template(intruder: np.ndarray, template_shape: tuple[int, int]) -> np.ndarray:
    height, width = template_shape
    if intruder.shape == template_shape:
        return intruder
    if intruder.shape == (height + 2, width + 2):
        # Intruder stack is Template FOV plus a 1-pixel border of background 142.
        return intruder[1 : height + 1, 1 : width + 1].copy()
    raise ValueError(
        f"Cannot align Intruder shape {intruder.shape} to Template shape {template_shape}"
    )


def y_to_t(y: np.ndarray) -> np.ndarray:
    return y.astype(np.float32) / T_Y_SCALE - T_Y_OFFSET


def t_to_display_rgb(template_y: np.ndarray, t_value: np.ndarray, overlay: np.ndarray) -> np.ndarray:
    """False-color QC image: grayscale anatomy with signed-t highlights."""
    rgb = np.stack([template_y, template_y, template_y], axis=-1).astype(np.uint8)
    if not overlay.any():
        return rgb
    magnitude = np.clip(np.abs(t_value) / max(T_SIGNIFICANCE * 2.0, 1e-6), 0.0, 1.0)
    positive = overlay & (t_value > 0)
    negative = overlay & (t_value < 0)
    for mask, color in (
        (negative, np.array([220, 40, 40], dtype=np.float32)),
        (positive, np.array([40, 90, 220], dtype=np.float32)),
    ):
        if not mask.any():
            continue
        alpha = magnitude[mask][:, None]
        base = rgb[mask].astype(np.float32)
        rgb[mask] = np.clip((1.0 - alpha) * base + alpha * color, 0, 255).astype(np.uint8)
    return rgb


def tissue_mask_from_template(template_y: np.ndarray) -> np.ndarray:
    candidate = template_y > TEMPLATE_BACKGROUND
    if not candidate.any():
        return candidate
    candidate = ndimage.binary_closing(candidate, iterations=2)
    candidate = ndimage.binary_fill_holes(candidate)
    labels, count = ndimage.label(candidate)
    if not count:
        return candidate
    sizes = np.bincount(labels.ravel())
    keep_n = min(4, count)
    keep_labels = np.argsort(sizes[1:])[-keep_n:] + 1
    tissue = np.isin(labels, keep_labels)
    tissue = ndimage.binary_closing(tissue, iterations=1)
    tissue = ndimage.binary_fill_holes(tissue)
    return tissue


def load_template_intruder_pair(template_path: Path, intruder_path: Path) -> dict:
    template_y = load_grayscale_y(template_path)
    intruder_raw = load_grayscale_y(intruder_path)
    intruder_y = crop_intruder_to_template(intruder_raw, template_y.shape)
    tissue = tissue_mask_from_template(template_y)
    t_value = y_to_t(intruder_y)
    # Significant Intruder voxels inside anatomy; background 142 → t ≈ 0.
    overlay = tissue & (intruder_y != INTRUDER_BACKGROUND) & (
        np.abs(t_value) >= T_SIGNIFICANCE
    )
    branch = np.zeros(template_y.shape, dtype=np.int8)
    branch[overlay & (t_value < 0)] = -1
    branch[overlay & (t_value > 0)] = 1
    display_magnitude = np.full(template_y.shape, np.nan, dtype=np.float32)
    display_magnitude[overlay] = np.abs(t_value[overlay]) / max(T_SIGNIFICANCE, 1e-6)
    intruder_rgb = t_to_display_rgb(template_y, t_value, overlay)
    return {
        "template_y": template_y,
        "template_2dg": template_y.astype(np.float32) / TEMPLATE_2DG_SCALE,
        "tissue": tissue.astype(bool),
        "intruder_y": intruder_y,
        "t_value": t_value.astype(np.float32),
        "tiff_y": intruder_y.astype(np.float32),
        "intruder_rgb": intruder_rgb,
        "overlay": overlay,
        "overlay_branch": branch,
        "display_magnitude": display_magnitude,
        "template_tissue_fraction": float(tissue.mean()),
        "intruder_signal_fraction": float((intruder_y != INTRUDER_BACKGROUND).mean()),
        "significant_overlay_fraction": float(overlay.mean()),
        "intruder_raw_shape": tuple(intruder_raw.shape),
        "blank_template": not bool(tissue.any()),
    }


def nominal_allen_index(bregma_mm: float) -> int:
    return int(round(CCF_BREGMA_ZERO_INDEX - bregma_mm / (CCF_VOXEL_UM / 1000.0)))


def load_manual_bregma_anchors() -> pd.DataFrame:
    if not MANUAL_MATCHES_PATH.exists():
        raise FileNotFoundError(
            f"Manual Template↔Allen matches CSV not found: {MANUAL_MATCHES_PATH}"
        )
    anchors = pd.read_csv(MANUAL_MATCHES_PATH)
    required = {"template_index", "bregma"}
    missing = required - set(anchors.columns)
    if missing:
        raise ValueError(
            f"Manual matches CSV missing columns {sorted(missing)}: {MANUAL_MATCHES_PATH}"
        )
    anchors = anchors.rename(
        columns={"template_index": "source_number", "bregma": "manual_bregma_mm"}
    )
    anchors["source_number"] = anchors["source_number"].astype(int)
    anchors["manual_bregma_mm"] = pd.to_numeric(
        anchors["manual_bregma_mm"], errors="coerce"
    )
    if anchors["source_number"].duplicated().any():
        raise ValueError("Manual matches CSV has duplicated template_index values")
    return anchors[["source_number", "manual_bregma_mm"]]


def build_manifest() -> pd.DataFrame:
    rows = []
    pairs = source_files()
    anchors = load_manual_bregma_anchors()
    tissue_source_numbers = []
    for template_path, _ in pairs:
        source_number = int(template_path.stem.replace("Template", ""))
        preview = load_grayscale_y(template_path)
        if (preview > TEMPLATE_BACKGROUND).any():
            tissue_source_numbers.append(source_number)
    first_tissue = (
        min(tissue_source_numbers)
        if tissue_source_numbers
        else FIRST_TISSUE_SOURCE_NUMBER
    )
    for position, (template_path, intruder_path) in enumerate(pairs):
        rebuilt = load_template_intruder_pair(template_path, intruder_path)
        source_number = int(template_path.stem.replace("Template", ""))
        anchor_rows = anchors.loc[anchors.source_number == source_number]
        if len(anchor_rows):
            manual_bregma = anchor_rows.manual_bregma_mm.iloc[0]
        else:
            manual_bregma = np.nan
        has_manual_anchor = bool(pd.notna(manual_bregma))
        if has_manual_anchor:
            bregma = float(manual_bregma)
            allen_index = nominal_allen_index(bregma)
            anchor_source = "manual_csv"
        else:
            bregma = np.nan
            allen_index = -1
            anchor_source = "missing"
        skip_registration = bool(rebuilt["blank_template"] or not has_manual_anchor)
        rows.append(
            {
                "slice_position": position,
                "source_number": source_number,
                "slice_id": f"2dg_{source_number:03d}",
                "template_path": str(template_path),
                "intruder_path": str(intruder_path),
                "height_px": int(rebuilt["template_y"].shape[0]),
                "width_px": int(rebuilt["template_y"].shape[1]),
                "intruder_raw_height_px": int(rebuilt["intruder_raw_shape"][0]),
                "intruder_raw_width_px": int(rebuilt["intruder_raw_shape"][1]),
                "manual_bregma_mm": (
                    float(manual_bregma) if has_manual_anchor else np.nan
                ),
                "anchor_source": anchor_source,
                "nominal_bregma_mm": bregma,
                "nominal_allen_index": allen_index,
                "template_tissue_fraction": rebuilt["template_tissue_fraction"],
                "intruder_signal_fraction": rebuilt["intruder_signal_fraction"],
                "significant_overlay_fraction": rebuilt["significant_overlay_fraction"],
                "blank_template": rebuilt["blank_template"],
                "skip_registration": skip_registration,
                "email_expected_slice_count": EXPECTED_EMAIL_SLICE_COUNT,
                "observed_slice_count": len(pairs),
                "n_tissue_slices": len(tissue_source_numbers),
                "first_tissue_source_number": first_tissue,
                "slice_count_matches_email": len(tissue_source_numbers)
                == EXPECTED_EMAIL_SLICE_COUNT,
                "source_image_left_is_anatomical_right": True,
                "source_flipped_for_ccf": FLIP_SOURCE_LEFT_RIGHT,
                "t_value_status": T_VALUE_STATUS,
                "t_formula": "t = (y / 18.409) - 7.762",
                "template_2dg_formula": "2DG = y / 1.494",
                "t_significance_threshold": T_SIGNIFICANCE,
                "manual_matches_path": str(MANUAL_MATCHES_PATH),
                "ap_search_offsets": ",".join(str(x) for x in AP_SEARCH_OFFSETS),
            }
        )
    return pd.DataFrame(rows)


def prepared_path(output_root: Path, slice_position: int) -> Path:
    return output_root / "prepared" / f"slice_{slice_position:03d}.npz"


def manifest_path(output_root: Path) -> Path:
    return output_root / f"2dg_slice_manifest_{ANALYSIS_DATE}.csv"


def prepare(output_root: Path, figure_root: Path) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "prepared").mkdir(parents=True, exist_ok=True)
    (output_root / "work").mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)
    manifest_file = manifest_path(output_root)
    if manifest_file.exists():
        manifest = pd.read_csv(manifest_file)
        log(f"Using existing manifest without modification: {manifest_file}")
    else:
        manifest = build_manifest()
        write_csv_new(manifest_file, manifest)
        log(f"Wrote {len(manifest)}-slice manifest: {manifest_file}")

    for row in manifest.itertuples(index=False):
        destination = prepared_path(output_root, int(row.slice_position))
        if destination.exists():
            continue
        rebuilt = load_template_intruder_pair(
            Path(row.template_path), Path(row.intruder_path)
        )
        arrays = {
            "template_y": rebuilt["template_y"],
            "template_2dg": rebuilt["template_2dg"],
            "tissue": rebuilt["tissue"].astype(np.uint8),
            "intruder_y": rebuilt["intruder_y"],
            "tiff_y": rebuilt["tiff_y"],
            "t_value": rebuilt["t_value"],
            "intruder_rgb": rebuilt["intruder_rgb"],
            "overlay": rebuilt["overlay"].astype(np.uint8),
            "overlay_branch": rebuilt["overlay_branch"],
            "display_magnitude": rebuilt["display_magnitude"],
        }
        if FLIP_SOURCE_LEFT_RIGHT:
            arrays = {name: np.flip(value, axis=1).copy() for name, value in arrays.items()}
        np.savez_compressed(destination, **arrays)
    write_preparation_qc(output_root, manifest)
    write_preparation_figure(output_root, figure_root, manifest)
    log(f"Prepared {len(manifest)} slices under {output_root / 'prepared'}")


def write_preparation_qc(output_root: Path, manifest: pd.DataFrame) -> None:
    path = output_root / f"2dg_preparation_qc_{ANALYSIS_DATE}.json"
    if path.exists():
        return
    tissue = manifest.loc[~manifest.blank_template.astype(bool)]
    payload = {
        "analysis_date": ANALYSIS_DATE,
        "primary_inputs": {
            "template_root": str(TEMPLATE_ROOT),
            "intruder_root": str(INTRUDER_SCALAR_ROOT),
        },
        "legacy_rgb_archives": {
            "intruder": str(LEGACY_INTRUDER_RGB_ROOT),
            "baseline": str(LEGACY_BASELINE_RGB_ROOT),
        },
        "observed_slice_count": int(len(manifest)),
        "n_tissue_slices": int((~manifest.blank_template.astype(bool)).sum()),
        "email_expected_slice_count": EXPECTED_EMAIL_SLICE_COUNT,
        "slice_count_matches_email": bool(
            int((~manifest.blank_template.astype(bool)).sum())
            == EXPECTED_EMAIL_SLICE_COUNT
        ),
        "first_tissue_source_number": int(manifest.first_tissue_source_number.iloc[0]),
        "manual_matches_path": str(MANUAL_MATCHES_PATH),
        "n_manual_anchors": int(
            manifest.anchor_source.eq("manual_csv").sum()
            if "anchor_source" in manifest.columns
            else 0
        ),
        "ap_search_offsets": list(AP_SEARCH_OFFSETS),
        "first_nominal_bregma_mm": float(
            tissue.loc[tissue.nominal_bregma_mm.notna(), "nominal_bregma_mm"].iloc[0]
        )
        if len(tissue.loc[tissue.nominal_bregma_mm.notna()])
        else None,
        "last_nominal_bregma_mm": float(
            tissue.loc[tissue.nominal_bregma_mm.notna(), "nominal_bregma_mm"].iloc[-1]
        )
        if len(tissue.loc[tissue.nominal_bregma_mm.notna()])
        else None,
        "spacing_mm": NOMINAL_SECTION_SPACING_MM,
        "source_orientation": "image left is anatomical right (email)",
        "ccf_orientation_action": "horizontal flip before registration",
        "t_value_status": T_VALUE_STATUS,
        "t_formula": "t = (y / 18.409) - 7.762",
        "template_2dg_formula": "2DG = y / 1.494",
        "intruder_crop": "center-crop 1px border (277x377 -> 275x375)",
        "significant_overlay_rule": (
            f"tissue & intruder_y != {INTRUDER_BACKGROUND} & |t| >= {T_SIGNIFICANCE}"
        ),
        "mean_significant_overlay_fraction": float(
            manifest.significant_overlay_fraction.mean()
        ),
    }
    write_json_new(path, payload)


def write_preparation_figure(
    output_root: Path, figure_root: Path, manifest: pd.DataFrame
) -> None:
    path = figure_root / f"2DG_template_intruder_qc_{ANALYSIS_DATE}.png"
    if path.exists():
        return
    tissue_positions = manifest.loc[
        ~manifest.blank_template.astype(bool), "slice_position"
    ].astype(int).tolist()
    if not tissue_positions:
        raise RuntimeError("No tissue slices available for preparation QC figure.")
    positions = sorted(
        {
            tissue_positions[0],
            tissue_positions[len(tissue_positions) // 2],
            tissue_positions[-1],
        }
    )
    fig, axes = plt.subplots(len(positions), 3, figsize=(12, 3.6 * len(positions)))
    if len(positions) == 1:
        axes = np.asarray([axes])
    for row_index, position in enumerate(positions):
        prepared = np.load(prepared_path(output_root, position))
        template = prepared["template_y"]
        display = prepared["intruder_rgb"]
        t_map = prepared["t_value"]
        tissue = prepared["tissue"].astype(bool)
        axes[row_index, 0].imshow(template, cmap="gray", vmin=0, vmax=255)
        axes[row_index, 0].set_title(f"Slice {position}: Template anatomy")
        axes[row_index, 1].imshow(
            np.ma.masked_where(~tissue, t_map),
            cmap="coolwarm",
            vmin=-T_SIGNIFICANCE * 2,
            vmax=T_SIGNIFICANCE * 2,
        )
        axes[row_index, 1].set_title("Scalar Intruder t map")
        axes[row_index, 2].imshow(display)
        axes[row_index, 2].set_title(f"|t| ≥ {T_SIGNIFICANCE} overlay")
        for axis in axes[row_index]:
            axis.axis("off")
    fig.suptitle(
        "Primary 2DG inputs: Template anatomy + WildRvSPF Intruder scalar t",
        y=1.0,
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def load_atlas() -> Atlas:
    global _ATLAS
    if _ATLAS is None:
        log("Loading Allen 25 um Nissl and annotation volumes.")
        nissl, _ = nrrd.read(str(NISSL_PATH), index_order="F")
        annotation, _ = nrrd.read(str(ANNOTATION_PATH), index_order="F")
        _ATLAS = Atlas(
            nissl=np.asarray(nissl, dtype=np.float32),
            annotation=np.asarray(annotation),
        )
    return _ATLAS


def normalize_image(image: np.ndarray, mask: np.ndarray) -> np.ndarray:
    result = np.zeros(image.shape, dtype=np.float32)
    values = np.asarray(image, dtype=np.float32)[mask]
    if not len(values):
        return result
    low, high = np.quantile(values, [0.01, 0.99])
    if high <= low:
        high = low + 1.0
    result = np.clip((image.astype(np.float32) - low) / (high - low), 0, 1)
    result[~mask] = 0
    return result


def parameter_map(kind: str, iterations: int):
    if kind == "similarity":
        result = sitk.GetDefaultParameterMap("rigid")
        result["Transform"] = ("SimilarityTransform",)
    elif kind == "affine":
        result = sitk.GetDefaultParameterMap("affine")
    else:
        raise ValueError(kind)
    result["NumberOfResolutions"] = ("4",)
    result["MaximumNumberOfIterations"] = (str(iterations),)
    result["NumberOfHistogramBins"] = ("32",)
    result["Metric"] = ("AdvancedMeanSquares",)
    result["AutomaticTransformInitialization"] = ("true",)
    result["AutomaticTransformInitializationMethod"] = ("CenterOfGravity",)
    result["AutomaticScalesEstimation"] = ("true",)
    result["ErodeMask"] = ("false",)
    result["DefaultPixelValue"] = ("0",)
    return result


def register_images(
    fixed: np.ndarray,
    moving: np.ndarray,
    fixed_mask: np.ndarray | None,
    moving_mask: np.ndarray | None,
    maps: list,
    work_dir: Path,
) -> tuple[list, np.ndarray]:
    work_dir.mkdir(parents=True, exist_ok=True)
    elastix = sitk.ElastixImageFilter()
    elastix.LogToConsoleOff()
    elastix.LogToFileOn()
    elastix.SetOutputDirectory(str(work_dir))
    elastix.SetFixedImage(sitk.GetImageFromArray(fixed.astype(np.float32)))
    elastix.SetMovingImage(sitk.GetImageFromArray(moving.astype(np.float32)))
    if fixed_mask is not None:
        elastix.SetFixedMask(
            sitk.Cast(sitk.GetImageFromArray(fixed_mask.astype(np.uint8)), sitk.sitkUInt8)
        )
    if moving_mask is not None:
        elastix.SetMovingMask(
            sitk.Cast(sitk.GetImageFromArray(moving_mask.astype(np.uint8)), sitk.sitkUInt8)
        )
    elastix.SetParameterMap(maps[0])
    for item in maps[1:]:
        elastix.AddParameterMap(item)
    elastix.Execute()
    transform = list(elastix.GetTransformParameterMap())
    centre = tuple(str((size - 1) / 2) for size in reversed(fixed.shape))
    for item in transform:
        if "CenterOfRotationPoint" not in item:
            item["CenterOfRotationPoint"] = centre
    return transform, sitk.GetArrayFromImage(elastix.GetResultImage())


def transform_image(image: np.ndarray, transform: list, nearest: bool = False) -> np.ndarray:
    maps = []
    for parameter_map_item in transform:
        copied = sitk.ParameterMap()
        for key, value in parameter_map_item.items():
            copied[key] = tuple(value)
        if nearest:
            copied["ResampleInterpolator"] = ("FinalNearestNeighborInterpolator",)
        maps.append(copied)
    transformix = sitk.TransformixImageFilter()
    transformix.LogToConsoleOff()
    transformix.SetTransformParameterMap(maps)
    transformix.SetMovingImage(sitk.GetImageFromArray(np.asarray(image, dtype=np.float32)))
    transformix.Execute()
    return sitk.GetArrayFromImage(transformix.GetResultImage())


def write_points_file(path: Path, points_xy: np.ndarray) -> None:
    with path.open("w") as handle:
        handle.write("point\n")
        handle.write(f"{len(points_xy)}\n")
        for x, y in points_xy:
            handle.write(f"{float(x):.6f} {float(y):.6f}\n")


def transform_points(
    points_xy: np.ndarray, transform: list, moving: np.ndarray, work_dir: Path
) -> np.ndarray:
    work_dir.mkdir(parents=True, exist_ok=True)
    point_file = work_dir / "source_points.pts"
    write_points_file(point_file, points_xy)
    transformix = sitk.TransformixImageFilter()
    transformix.LogToConsoleOff()
    transformix.SetTransformParameterMap(transform)
    transformix.SetMovingImage(sitk.GetImageFromArray(moving.astype(np.float32)))
    transformix.SetFixedPointSetFileName(str(point_file))
    transformix.SetOutputDirectory(str(work_dir))
    transformix.Execute()
    output_file = work_dir / "outputpoints.txt"
    mapped = []
    with output_file.open() as handle:
        for line in handle:
            token = line.split("OutputPoint = [", 1)[1].split("]", 1)[0]
            mapped.append([float(value) for value in token.split()])
    return np.asarray(mapped, dtype=float)


def dice_score(left: np.ndarray, right: np.ndarray) -> float:
    left = left.astype(bool)
    right = right.astype(bool)
    denominator = left.sum() + right.sum()
    return float(2 * np.logical_and(left, right).sum() / denominator) if denominator else 0.0


def normalized_mutual_information(left: np.ndarray, right: np.ndarray, mask: np.ndarray) -> float:
    x = left[mask].astype(float)
    y = right[mask].astype(float)
    if len(x) < 100:
        return float("nan")
    histogram, _, _ = np.histogram2d(x, y, bins=32)
    probability = histogram / histogram.sum()
    px = probability.sum(axis=1)
    py = probability.sum(axis=0)
    nz = probability > 0
    mutual_information = np.sum(
        probability[nz]
        * np.log(probability[nz] / (px[:, None] * py[None, :])[nz])
    )
    hx = -np.sum(px[px > 0] * np.log(px[px > 0]))
    hy = -np.sum(py[py > 0] * np.log(py[py > 0]))
    return float(2 * mutual_information / (hx + hy)) if (hx + hy) else 0.0


def affine_jacobian_summary(
    transform: list, shape: tuple[int, int], moving: np.ndarray, work_dir: Path
) -> dict:
    step = 24
    rows = np.arange(step, shape[0] - step, step)
    columns = np.arange(step, shape[1] - step, step)
    xx, yy = np.meshgrid(columns, rows)
    points = np.column_stack([xx.ravel(), yy.ravel()])
    mapped = transform_points(points, transform, moving, work_dir)
    mx = mapped[:, 0].reshape(xx.shape)
    my = mapped[:, 1].reshape(yy.shape)
    dxdc, dxdr = np.gradient(mx, step, step, axis=(1, 0))
    dydc, dydr = np.gradient(my, step, step, axis=(1, 0))
    determinant = dxdc * dydr - dxdr * dydc
    return {
        "jacobian_min": float(np.nanmin(determinant)),
        "jacobian_max": float(np.nanmax(determinant)),
        "jacobian_mean": float(np.nanmean(determinant)),
        "jacobian_negative_fraction": float(np.mean(determinant <= 0)),
    }


def slice_work_dir(output_root: Path, position: int) -> Path:
    return output_root / "work" / f"slice_{position:03d}"


def fit_candidate(
    fixed: np.ndarray,
    fixed_mask: np.ndarray,
    atlas: Atlas,
    atlas_index: int,
    work_dir: Path,
    iterations: int,
    final: bool,
) -> dict:
    annotation = atlas.annotation[atlas_index]
    moving_mask = annotation > 0
    moving_nissl = normalize_image(np.log1p(atlas.nissl[atlas_index]), moving_mask)
    # Cross-modality intensity registration allowed the much larger Allen brain
    # to surround the 2DG section without matching its outline. Register smooth
    # within-tissue distance fields instead; Nissl remains an independent QC
    # image and is never used to manufacture correspondence.
    moving = ndimage.gaussian_filter(moving_mask.astype(np.float32), sigma=3)
    maps = [parameter_map("similarity", iterations)]
    if final:
        affine = parameter_map("affine", max(250, iterations))
        affine["AutomaticTransformInitialization"] = ("false",)
        maps.append(affine)
    transform, warped = register_images(
        fixed, moving, None, None, maps, work_dir
    )
    warped_mask = transform_image(moving_mask.astype(np.float32), transform, nearest=True) > 0.5
    warped_nissl = transform_image(moving_nissl, transform)
    overlap = fixed_mask | warped_mask
    return {
        "atlas_index": atlas_index,
        "transform": transform,
        "warped": warped,
        "warped_nissl": warped_nissl,
        "warped_mask": warped_mask,
        "dice": dice_score(fixed_mask, warped_mask),
        "nmi": normalized_mutual_information(fixed, warped, overlap),
        "moving": moving,
        "annotation": annotation,
    }


def annotation_name_map() -> dict[int, str]:
    with ANNOTATION_NAMES_PATH.open("rb") as handle:
        hierarchy = pickle.load(handle)
    direct_names: dict[int, str] = {}
    # Prefer direct/leaf names (a singleton annotation list), then the shortest
    # available label for IDs that only occur in a broader hierarchy entry.
    for name, values in hierarchy.items():
        if len(values) == 1:
            direct_names[int(values[0])] = str(name)
    names = dict(direct_names)
    specificity = {value: 1 for value in direct_names}
    for name, values in hierarchy.items():
        for value in values:
            value = int(value)
            if value not in direct_names and (
                value not in specificity or len(values) < specificity[value]
            ):
                names[value] = str(name)
                specificity[value] = len(values)
    names[0] = "outside atlas"
    return names


def write_annotation_map(table_root: Path) -> None:
    path = table_root / f"Allen_CCF_annotation_names_hierarchy_{ANALYSIS_DATE}.csv"
    if path.exists():
        log(f"Using existing Allen annotation name map without modification: {path}")
        return
    names = annotation_name_map()
    frame = pd.DataFrame(
        {
            "annotation": list(names.keys()),
            "annotation_name": list(names.values()),
        }
    ).sort_values("annotation")
    write_csv_new(path, frame)
    log(f"Wrote Allen annotation name map: {path}")


def export_scalar_t_overlay(output_root: Path, table_root: Path) -> None:
    """Copy strict overlay points into the scalar-t artifact expected by R."""
    source = output_root / f"2dg_overlay_points_ccf_strict_{ANALYSIS_DATE}.csv.gz"
    destination = (
        output_root
        / f"2dg_overlay_points_ccf_scalar-t_{EFFECT_ANALYSIS_DATE}.csv.gz"
    )
    provenance_path = (
        output_root / f"2dg_scalar_t_export_provenance_{EFFECT_ANALYSIS_DATE}.json"
    )
    for path in (destination, provenance_path):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    if not source.exists():
        raise FileNotFoundError(f"Strict CCF overlay table is missing: {source}")
    overlay = pd.read_csv(source)
    if overlay.empty:
        raise ValueError("Strict overlay table is empty; cannot export scalar t points.")
    if not np.isfinite(overlay["t_value"]).all():
        raise ValueError("Strict overlay table contains non-finite t_value entries.")
    overlay = overlay.copy()
    overlay["t_value_status"] = T_VALUE_STATUS
    write_csv_new(destination, overlay, compression="gzip")
    write_json_new(
        provenance_path,
        {
            "analysis_date": EFFECT_ANALYSIS_DATE,
            "source_overlay": source,
            "scalar_overlay": destination,
            "t_formula": "t = (y / 18.409) - 7.762",
            "t_value_status": T_VALUE_STATUS,
            "n_overlay_pixels": int(len(overlay)),
            "t_min": float(overlay.t_value.min()),
            "t_max": float(overlay.t_value.max()),
            "note": (
                "Primary WildRvSPF Intruder stack already stores scalar y; "
                "no RGB palette decode is required."
            ),
        },
    )
    log(f"Exported {len(overlay)} scalar-t CCF overlay pixels: {destination}")


def register_slice(output_root: Path, slice_position: int) -> None:
    manifest = pd.read_csv(manifest_path(output_root))
    if slice_position < 0 or slice_position >= len(manifest):
        raise IndexError(f"slice-position must be between 0 and {len(manifest) - 1}")
    row = manifest.iloc[slice_position]
    work_dir = slice_work_dir(output_root, slice_position)
    metrics_path = work_dir / "registration_metrics.csv"
    pixel_path = work_dir / "pixels_ccf.csv.gz"
    qc_path = work_dir / "registration_qc.png"
    if metrics_path.exists() and pixel_path.exists():
        log(
            f"Slice {slice_position} already complete under {work_dir}; "
            "leaving existing artifacts untouched."
        )
        return
    if metrics_path.exists() or pixel_path.exists() or qc_path.exists():
        raise FileExistsError(
            f"Incomplete existing slice result under {work_dir}; refusing to overwrite"
        )
    work_dir.mkdir(parents=True, exist_ok=True)

    prepared = np.load(prepared_path(output_root, slice_position))
    tissue = prepared["tissue"].astype(bool)
    blank_or_empty = bool(getattr(row, "blank_template", False)) or not tissue.any()
    missing_anchor = bool(getattr(row, "skip_registration", False)) and not blank_or_empty
    if blank_or_empty or missing_anchor:
        nominal_bregma = (
            float(row.nominal_bregma_mm)
            if pd.notna(row.nominal_bregma_mm)
            else float("nan")
        )
        nominal_index = (
            int(row.nominal_allen_index)
            if pd.notna(row.nominal_allen_index) and int(row.nominal_allen_index) >= 0
            else -1
        )
        metrics = {
            "slice_position": slice_position,
            "source_number": int(row.source_number),
            "manual_bregma_mm": nominal_bregma,
            "anchor_source": getattr(row, "anchor_source", "missing"),
            "nominal_bregma_mm": nominal_bregma,
            "nominal_allen_index": nominal_index,
            "selected_allen_index": nominal_index,
            "selected_bregma_mm": nominal_bregma,
            "ap_offset_voxels": 0,
            "dice": 0.0,
            "nmi": float("nan"),
            "inside_allen_tissue_fraction": 0.0,
            "n_tissue_pixels": 0,
            "n_overlay_pixels": 0,
            "accepted_slice_fit": False,
            "rejection_reason": (
                "blank_or_empty_template" if blank_or_empty else "no_manual_anchor"
            ),
            "jacobian_min": float("nan"),
            "jacobian_max": float("nan"),
            "jacobian_mean": float("nan"),
            "jacobian_negative_fraction": float("nan"),
        }
        write_csv_new(metrics_path, pd.DataFrame([metrics]))
        write_csv_new(
            pixel_path,
            pd.DataFrame(
                columns=[
                    "slice_position",
                    "source_number",
                    "nominal_bregma_mm",
                    "selected_allen_index",
                    "ccfx",
                    "source_x",
                    "source_y",
                    "registration_x",
                    "registration_y",
                    "atlas_x",
                    "atlas_y",
                    "ccfy",
                    "ccfz",
                    "annotation",
                    "annotation_name",
                    "template_y",
                    "template_2dg",
                    "overlay_present",
                    "overlay_branch",
                    "overlay_signed_relative",
                    "overlay_r",
                    "overlay_g",
                    "overlay_b",
                    "tiff_y",
                    "t_value",
                    "t_value_status",
                    "inside_allen_tissue",
                ]
            ),
            compression="gzip",
        )
        write_json_new(
            work_dir / "provenance.json",
            {
                "metrics": metrics,
                "input_template": row.template_path,
                "input_intruder": row.intruder_path,
                "skipped": True,
            },
        )
        log(
            f"Slice {slice_position} skipped "
            f"({metrics['rejection_reason']}); metrics={metrics_path}"
        )
        return
    fixed = ndimage.gaussian_filter(tissue.astype(np.float32), sigma=3)
    atlas = load_atlas()
    nominal = int(row.nominal_allen_index)
    if nominal < 0:
        raise ValueError(
            f"Slice {slice_position} lacks a usable manual Allen AP anchor"
        )
    candidate_indices = sorted(
        {
            int(np.clip(nominal + offset, 0, atlas.annotation.shape[0] - 1))
            for offset in AP_SEARCH_OFFSETS
        }
    )
    candidates = []
    for atlas_index in candidate_indices:
        log(f"Slice {slice_position}: testing Allen AP index {atlas_index}")
        candidate = fit_candidate(
            fixed,
            tissue,
            atlas,
            atlas_index,
            work_dir / "search" / f"allen_{atlas_index:03d}",
            iterations=100,
            final=False,
        )
        candidates.append(candidate)
    best_search = max(
        candidates,
        key=lambda item: item["dice"] + 0.25 * np.nan_to_num(item["nmi"]),
    )
    selected_index = int(best_search["atlas_index"])
    log(f"Slice {slice_position}: final fit at Allen AP index {selected_index}")
    fit = fit_candidate(
        fixed,
        tissue,
        atlas,
        selected_index,
        work_dir / "final_fit",
        iterations=500,
        final=True,
    )
    jacobian = affine_jacobian_summary(
        fit["transform"], fixed.shape, fit["moving"], work_dir / "jacobian"
    )

    yy, xx = np.where(tissue)
    source_points = np.column_stack([xx, yy]).astype(float)
    mapped = transform_points(
        source_points, fit["transform"], fit["moving"], work_dir / "map_pixels"
    )
    atlas_x = mapped[:, 0]
    atlas_y = mapped[:, 1]
    xi = np.rint(atlas_x).astype(int)
    yi = np.rint(atlas_y).astype(int)
    inside_bounds = (
        (xi >= 0)
        & (xi < fit["annotation"].shape[1])
        & (yi >= 0)
        & (yi < fit["annotation"].shape[0])
    )
    annotation = np.zeros(len(xi), dtype=np.int64)
    annotation[inside_bounds] = fit["annotation"][yi[inside_bounds], xi[inside_bounds]]
    inside_tissue = inside_bounds & (annotation > 0)
    inside_fraction = float(inside_tissue.mean())

    names = annotation_name_map()
    rgb = prepared["intruder_rgb"][yy, xx]
    branch = prepared["overlay_branch"][yy, xx].astype(int)
    magnitude = prepared["display_magnitude"][yy, xx].astype(float)
    if "tiff_y" in prepared.files:
        tiff_y = prepared["tiff_y"][yy, xx].astype(float)
    else:
        tiff_y = np.full(len(yy), np.nan)
    if "t_value" in prepared.files:
        t_value = prepared["t_value"][yy, xx].astype(float)
    else:
        t_value = np.full(len(yy), np.nan)
    t_value_status = getattr(row, "t_value_status", T_VALUE_STATUS)
    original_x = prepared["template_y"].shape[1] - 1 - xx if FLIP_SOURCE_LEFT_RIGHT else xx
    pixels = pd.DataFrame(
        {
            "slice_position": slice_position,
            "source_number": int(row.source_number),
            "nominal_bregma_mm": float(row.nominal_bregma_mm),
            "selected_allen_index": selected_index,
            "ccfx": selected_index * CCF_VOXEL_UM,
            "source_x": original_x,
            "source_y": yy,
            "registration_x": xx,
            "registration_y": yy,
            "atlas_x": atlas_x,
            "atlas_y": atlas_y,
            "ccfy": atlas_y * CCF_VOXEL_UM,
            "ccfz": atlas_x * CCF_VOXEL_UM,
            "annotation": annotation,
            "annotation_name": [names.get(int(value), f"Allen ID {value}") for value in annotation],
            "template_y": prepared["template_y"][yy, xx],
            "template_2dg": prepared["template_2dg"][yy, xx],
            "overlay_present": branch != 0,
            "overlay_branch": np.select(
                [branch == -1, branch == 1, branch == 2],
                ["red", "blue", "other"],
                default="none",
            ),
            "overlay_signed_relative": np.where(
                branch == 0, np.nan, np.where(branch == -1, -magnitude, magnitude)
            ),
            "overlay_r": rgb[:, 0],
            "overlay_g": rgb[:, 1],
            "overlay_b": rgb[:, 2],
            "tiff_y": tiff_y,
            "t_value": t_value,
            "t_value_status": t_value_status,
            "inside_allen_tissue": inside_tissue,
        }
    )
    write_csv_new(pixel_path, pixels, compression="gzip")

    accepted = bool(
        fit["dice"] >= MIN_DICE
        and inside_fraction >= MIN_INSIDE_FRACTION
        and jacobian["jacobian_min"] >= MIN_JACOBIAN
        and jacobian["jacobian_max"] <= MAX_JACOBIAN
        and jacobian["jacobian_negative_fraction"] == 0
    )
    metrics = {
        "slice_position": slice_position,
        "source_number": int(row.source_number),
        "manual_bregma_mm": float(row.nominal_bregma_mm),
        "anchor_source": getattr(row, "anchor_source", "manual_csv"),
        "nominal_bregma_mm": float(row.nominal_bregma_mm),
        "nominal_allen_index": nominal,
        "selected_allen_index": selected_index,
        "selected_bregma_mm": (CCF_BREGMA_ZERO_INDEX - selected_index) * 0.025,
        "ap_offset_voxels": selected_index - nominal,
        "dice": fit["dice"],
        "nmi": fit["nmi"],
        "inside_allen_tissue_fraction": inside_fraction,
        "n_tissue_pixels": len(pixels),
        "n_overlay_pixels": int(pixels.overlay_present.sum()),
        "accepted_slice_fit": accepted,
        "rejection_reason": "" if accepted else "one or more prespecified QC gates failed",
        **jacobian,
    }
    write_csv_new(metrics_path, pd.DataFrame([metrics]))
    write_transform_maps(work_dir / "transform", fit["transform"])
    write_slice_qc(qc_path, prepared, fit, metrics)
    write_json_new(
        work_dir / "provenance.json",
        {
            "metrics": metrics,
            "search_candidates": [
                {
                    "atlas_index": item["atlas_index"],
                    "dice": item["dice"],
                    "nmi": item["nmi"],
                }
                for item in candidates
            ],
            "input_template": row.template_path,
            "input_intruder": row.intruder_path,
            "source_flipped_left_right": FLIP_SOURCE_LEFT_RIGHT,
            "value_conversions": {
                "template_2dg": "template_y / 1.494",
                "t_value": "t = (y / 18.409) - 7.762",
            },
        },
    )
    log(f"Slice {slice_position} complete; accepted={accepted}; metrics={metrics_path}")


def write_transform_maps(root: Path, transform: list) -> None:
    root.mkdir(parents=True, exist_ok=True)
    for index, parameter_map_item in enumerate(transform):
        path = root / f"TransformParameters.{index}.txt"
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite transform: {path}")
        sitk.WriteParameterFile(parameter_map_item, str(path))


def write_slice_qc(path: Path, prepared, fit: dict, metrics: dict) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite QC figure: {path}")
    boundaries = find_boundaries(fit["warped_mask"], mode="outer")
    overlay_rgb = prepared["intruder_rgb"].copy()
    template = prepared["template_y"]
    t_map = prepared["t_value"] if "t_value" in prepared.files else None
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    axes[0].imshow(template, cmap="gray", vmin=0, vmax=255)
    axes[0].set_title("Template anatomy")
    if t_map is not None:
        axes[1].imshow(
            np.ma.masked_where(~prepared["tissue"].astype(bool), t_map),
            cmap="coolwarm",
            vmin=-T_SIGNIFICANCE * 2,
            vmax=T_SIGNIFICANCE * 2,
        )
        axes[1].set_title("Scalar Intruder t")
    else:
        axes[1].imshow(overlay_rgb)
        axes[1].set_title("Intruder display")
    axes[2].imshow(template, cmap="gray", vmin=0, vmax=255)
    axes[2].contour(boundaries, levels=[0.5], colors="tab:orange", linewidths=0.7)
    axes[2].set_title("Allen tissue boundary")
    axes[3].imshow(fit["warped_nissl"], cmap="gray")
    axes[3].set_title(
        f"Allen Nissl warped\nDice={metrics['dice']:.2f}; inside={metrics['inside_allen_tissue_fraction']:.2f}"
    )
    for axis in axes:
        axis.axis("off")
    fig.suptitle(
        f"2DG slice {metrics['source_number']} → Allen {metrics['selected_allen_index']} "
        f"({'accepted' if metrics['accepted_slice_fit'] else 'rejected'})"
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def read_completed_metrics(output_root: Path) -> pd.DataFrame:
    paths = sorted((output_root / "work").glob("slice_*/registration_metrics.csv"))
    if not paths:
        raise FileNotFoundError("No completed per-slice registration metrics were found.")
    return pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)


def add_continuity_qc(metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = metrics.sort_values("slice_position").reset_index(drop=True)
    delta = metrics.selected_allen_index.diff()
    metrics["selected_ap_step_voxels"] = delta
    expected_step = NOMINAL_SECTION_SPACING_MM / (CCF_VOXEL_UM / 1000.0)
    max_offset = max(abs(x) for x in AP_SEARCH_OFFSETS)

    # Infer stack AP travel from accepted anchors/fits. Manual anchors run
    # posterior→anterior (Allen AP index decreases as slice_position increases).
    fit_ok = metrics.loc[
        metrics.accepted_slice_fit.astype(bool)
        & (metrics.selected_allen_index >= 0)
    ]
    if len(fit_ok) >= 2:
        ap_direction = int(
            np.sign(
                float(fit_ok.selected_allen_index.iloc[-1])
                - float(fit_ok.selected_allen_index.iloc[0])
            )
        )
    else:
        ap_direction = 1
    if ap_direction == 0:
        ap_direction = 1
    metrics["ap_travel_direction"] = ap_direction
    metrics["ap_monotonic"] = delta.isna() | (
        (delta * ap_direction) >= 0
    )
    metrics["ap_step_plausible"] = delta.isna() | (
        ((delta * ap_direction) >= 0)
        & ((delta.abs()) <= expected_step + 2 * max_offset)
    )
    metrics["accepted_final"] = False
    metrics["strict_ap_delta_from_last_accepted"] = np.nan
    metrics["continuity_note"] = ""
    last_index = None
    last_position = None
    for row_index in metrics.index:
        if not bool(metrics.at[row_index, "accepted_slice_fit"]):
            metrics.at[row_index, "continuity_note"] = "per-slice registration gate failed"
            continue
        current_index = int(metrics.at[row_index, "selected_allen_index"])
        if current_index < 0:
            metrics.at[row_index, "continuity_note"] = "missing Allen AP selection"
            continue
        current_position = int(metrics.at[row_index, "slice_position"])
        if last_index is None:
            metrics.at[row_index, "accepted_final"] = True
            last_index = current_index
            last_position = current_position
            continue
        strict_delta = current_index - last_index
        metrics.at[row_index, "strict_ap_delta_from_last_accepted"] = strict_delta
        position_delta = current_position - last_position
        maximum_step = position_delta * expected_step + 2 * max_offset
        directed_step = strict_delta * ap_direction
        if directed_step < 0 or abs(strict_delta) > maximum_step:
            metrics.at[row_index, "continuity_note"] = (
                "strict monotonic AP continuity gate failed"
            )
            continue
        metrics.at[row_index, "accepted_final"] = True
        last_index = current_index
        last_position = current_position
    return metrics


def aggregate(output_root: Path, figure_root: Path, table_root: Path) -> None:
    metrics = add_continuity_qc(read_completed_metrics(output_root))
    summary_path = output_root / f"2dg_registration_summary_strict_{ANALYSIS_DATE}.csv"
    write_csv_new(summary_path, metrics)
    accepted_positions = set(
        metrics.loc[metrics.accepted_final, "slice_position"].astype(int).tolist()
    )
    if not accepted_positions:
        raise RuntimeError("No slices passed final QC; aggregation stopped.")

    region_frames = []
    overlay_frames = []
    for position in sorted(accepted_positions):
        path = slice_work_dir(output_root, position) / "pixels_ccf.csv.gz"
        frame = pd.read_csv(path)
        inside = frame[frame.inside_allen_tissue.astype(bool) & (frame.annotation != 0)].copy()
        region = (
            inside.groupby(["annotation", "annotation_name"], as_index=False)
            .agg(
                n_template_pixels=("template_2dg", "size"),
                mean_template_2dg=("template_2dg", "mean"),
                median_template_2dg=("template_2dg", "median"),
                n_overlay_pixels=("overlay_present", "sum"),
                mean_overlay_signed_relative=("overlay_signed_relative", "mean"),
                median_overlay_signed_relative=("overlay_signed_relative", "median"),
                mean_t_value=("t_value", "mean"),
                median_t_value=("t_value", "median"),
            )
        )
        region["slice_position"] = position
        region["source_number"] = int(frame.source_number.iloc[0])
        region["selected_allen_index"] = int(frame.selected_allen_index.iloc[0])
        region_frames.append(region)
        overlay_frames.append(
            inside.loc[
                inside.overlay_present.astype(bool),
                [
                    "slice_position",
                    "source_number",
                    "ccfx",
                    "ccfy",
                    "ccfz",
                    "annotation",
                    "annotation_name",
                    "overlay_branch",
                    "overlay_signed_relative",
                    "overlay_r",
                    "overlay_g",
                    "overlay_b",
                    "tiff_y",
                    "t_value",
                    "t_value_status",
                ],
            ]
        )
    regions_by_slice = pd.concat(region_frames, ignore_index=True)
    overlay_points = pd.concat(overlay_frames, ignore_index=True)
    region_summary = (
        regions_by_slice.groupby(["annotation", "annotation_name"], as_index=False)
        .agg(
            n_slices=("slice_position", "nunique"),
            n_template_pixels=("n_template_pixels", "sum"),
            mean_template_2dg=("mean_template_2dg", "mean"),
            n_overlay_pixels=("n_overlay_pixels", "sum"),
            mean_overlay_signed_relative=("mean_overlay_signed_relative", "mean"),
            median_overlay_signed_relative=("median_overlay_signed_relative", "median"),
            mean_t_value=("mean_t_value", "mean"),
            median_t_value=("median_t_value", "median"),
        )
    )

    table_root.mkdir(parents=True, exist_ok=True)
    figure_root.mkdir(parents=True, exist_ok=True)
    write_csv_new(
        table_root / f"2DG_CCF_region_summary_strict_{ANALYSIS_DATE}.csv",
        region_summary,
    )
    write_csv_new(
        table_root / f"2DG_CCF_region_by_slice_strict_{ANALYSIS_DATE}.csv.gz",
        regions_by_slice,
        compression="gzip",
    )
    write_csv_new(
        output_root / f"2dg_overlay_points_ccf_strict_{ANALYSIS_DATE}.csv.gz",
        overlay_points,
        compression="gzip",
    )
    write_aggregate_figures(metrics, overlay_points, figure_root)
    write_json_new(
        output_root / f"2dg_aggregation_provenance_strict_{ANALYSIS_DATE}.json",
        {
            "n_completed_slices": int(len(metrics)),
            "n_accepted_slices": int(metrics.accepted_final.sum()),
            "accepted_slice_positions": sorted(accepted_positions),
            "n_overlay_points": int(len(overlay_points)),
            "n_allen_regions": int(len(region_summary)),
            "t_value_status": sorted(
                overlay_points.t_value_status.dropna().astype(str).unique().tolist()
            ),
            "registration_status": "provisional: TIFF spacing/origin assumed from email",
        },
    )
    log(f"Aggregated {len(accepted_positions)} accepted slices; summary={summary_path}")


def write_aggregate_figures(
    metrics: pd.DataFrame, overlay_points: pd.DataFrame, figure_root: Path
) -> None:
    qc_path = figure_root / f"2DG_CCF_registration_qc_strict_{ANALYSIS_DATE}.png"
    if qc_path.exists():
        raise FileExistsError(f"Refusing to overwrite figure: {qc_path}")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    colors = np.where(metrics.accepted_final, "tab:blue", "tab:red")
    axes[0].scatter(metrics.source_number, metrics.dice, c=colors, s=24)
    axes[0].axhline(MIN_DICE, color="0.4", linestyle="--", linewidth=1)
    axes[0].set_ylabel("Tissue Dice")
    axes[0].set_title("2DG Template → Allen CCF registration QC")
    axes[1].plot(
        metrics.source_number,
        metrics.nominal_allen_index,
        color="0.6",
        label="nominal",
    )
    axes[1].scatter(
        metrics.source_number,
        metrics.selected_allen_index,
        c=colors,
        s=24,
        label="selected",
    )
    axes[1].set_xlabel("Source slice number")
    axes[1].set_ylabel("Allen AP index")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(qc_path, dpi=200, facecolor="white")
    plt.close(fig)

    coverage_path = figure_root / f"2DG_CCF_overlay_3d_strict_{ANALYSIS_DATE}.png"
    if coverage_path.exists():
        raise FileExistsError(f"Refusing to overwrite figure: {coverage_path}")
    sample = overlay_points
    if len(sample) > 50000:
        sample = sample.sample(50000, random_state=42)
    fig = plt.figure(figsize=(10, 7))
    axis = fig.add_subplot(111, projection="3d")
    if "t_value" in sample.columns and np.isfinite(sample.t_value).any():
        color = np.where(sample.t_value >= 0, "tab:blue", "tab:red")
    else:
        color = np.where(sample.overlay_signed_relative >= 0, "tab:blue", "tab:red")
    axis.scatter(
        sample.ccfx,
        sample.ccfz,
        sample.ccfy,
        c=color,
        s=1.2,
        alpha=0.35,
        linewidth=0,
    )
    axis.set_xlabel("CCF AP (µm)")
    axis.set_ylabel("CCF ML (µm)")
    axis.set_zlabel("CCF DV (µm)")
    axis.invert_zaxis()
    axis.set_title("WildRvSPF Intruder scalar t overlay in Allen CCF")
    fig.tight_layout()
    fig.savefig(coverage_path, dpi=200, facecolor="white")
    plt.close(fig)


def status(output_root: Path) -> None:
    manifest = pd.read_csv(manifest_path(output_root))
    complete = sorted(
        int(path.parent.name.replace("slice_", ""))
        for path in (output_root / "work").glob("slice_*/registration_metrics.csv")
    )
    print(
        json.dumps(
            {
                "output_root": str(output_root),
                "manifest_slices": int(len(manifest)),
                "completed_slices": len(complete),
                "complete_positions": complete,
                "pending_positions": sorted(set(manifest.slice_position) - set(complete)),
            },
            indent=2,
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Dated interim output directory.",
    )
    parser.add_argument(
        "--figure-root",
        type=Path,
        default=FIGURE_ROOT,
        help="Final figure directory.",
    )
    parser.add_argument(
        "--table-root",
        type=Path,
        default=TABLE_ROOT,
        help="Final table directory.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("inventory")
    subparsers.add_parser("prepare")
    register_parser = subparsers.add_parser("register-slice")
    register_parser.add_argument("--slice-position", type=int, required=True)
    subparsers.add_parser("aggregate")
    subparsers.add_parser("annotation-map")
    subparsers.add_parser("export-scalar-t")
    # Alias retained so older notebook chunks still resolve.
    subparsers.add_parser("decode-t-overlay")
    subparsers.add_parser("status")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "inventory":
        args.output_root.mkdir(parents=True, exist_ok=True)
        write_csv_new(manifest_path(args.output_root), build_manifest())
    elif args.command == "prepare":
        prepare(args.output_root, args.figure_root)
    elif args.command == "register-slice":
        register_slice(args.output_root, args.slice_position)
    elif args.command == "aggregate":
        aggregate(args.output_root, args.figure_root, args.table_root)
    elif args.command == "annotation-map":
        write_annotation_map(args.table_root)
    elif args.command in {"export-scalar-t", "decode-t-overlay"}:
        export_scalar_t_overlay(args.output_root, args.table_root)
    elif args.command == "status":
        status(args.output_root)


if __name__ == "__main__":
    main()
