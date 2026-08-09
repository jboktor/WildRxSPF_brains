#!/usr/bin/env python
"""Regenerate a clean, complete figure suite for the final accepted CCF fits."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nrrd
import numpy as np
import pandas as pd
import skimage.measure
import skimage.io
from PIL import Image, ImageChops, ImageDraw


BASE = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
ALL_SLICE_SCRIPT = SCRIPT_DIR / "allen_ccf_all_slice_landmark_sweep.py"
TEST_ROOT = BASE / "data/interim/registration/Allen_CCF_regional_tests"
WORKBOOK = TEST_ROOT / "slice_positions_25um_all_slices_regional_landmarks.xlsx"
FINAL_DECISIONS = TEST_ROOT / "all_slices_regional_landmark_final_decisions.csv"
BASELINE_CCF = BASE / "data/interim/registration/Allen_CCF_optimization/ccf3d.csv"
REFERENCE_ROOT = BASE / "data/input/allen_registration_ref"
FIGURE_ROOT = BASE / "figures/Allen_CCF_alignment_optimized"
FINAL_ROOT = FIGURE_ROOT / "final"
WORK_ROOT = TEST_ROOT / "final_render_work"


def load_all_slice_module():
    spec = importlib.util.spec_from_file_location("all_slice_sweep", ALL_SLICE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {ALL_SLICE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def trim_white(image: Image.Image) -> Image.Image:
    rgb = image.convert("RGB")
    background = Image.new("RGB", rgb.size, "white")
    box = ImageChops.difference(rgb, background).getbbox()
    return rgb.crop(box) if box else rgb


def final_coordinates(
    sample_id: str, accepted_iteration: int, baseline_cells: pd.DataFrame
) -> pd.DataFrame:
    if accepted_iteration == 3:
        return baseline_cells[baseline_cells["sample_id"].eq(sample_id)].copy()
    run_dir = TEST_ROOT / f"iteration_{accepted_iteration:02d}" / sample_id
    accepted = run_dir / f"{sample_id}_accepted_ccf2d.csv"
    if not accepted.exists():
        accepted = run_dir / f"{sample_id}_ccf2d.csv"
    return pd.read_csv(accepted, index_col=0)


def plot_final_coordinate_figures(
    sample_id: str,
    status: str,
    cells: pd.DataFrame,
    annotation_25: np.ndarray,
    annotation_10: np.ndarray,
    all_slice,
    regional,
) -> None:
    atlas_z = int(np.rint(np.nanmedian(cells["ccfx"]) / 25))
    atlas_slice = annotation_25[atlas_z]
    title = f"{sample_id} — FINAL {status}"

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(
        cells["ccfz"] / 25,
        cells["ccfy"] / 25,
        s=0.08,
        c=cells["subclass_color"],
    )
    for contour in regional.make_atlas_outlines(atlas_slice):
        ax.plot(contour[:, 1], contour[:, 0], "k", linewidth=0.35)
    ax.set(xlim=(0, atlas_slice.shape[1]), ylim=(atlas_slice.shape[0], 0), title=title)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FINAL_ROOT / f"{sample_id} 2d reg.png", dpi=300)
    plt.close(fig)

    sampled = regional.sample_true_atlas_surface(cells, annotation_10)
    fig, ax = plt.subplots(figsize=(7, 10))
    ax.scatter(
        cells["ccfz"] / 10,
        cells["ccfy"] / 10,
        s=0.5,
        c=cells["subclass_color"],
    )
    for contour in regional.make_atlas_outlines(sampled):
        ax.plot(contour[:, 1], contour[:, 0], "k", linewidth=1)
    ax.set(xlim=(0, 1140), ylim=(800, 0), title=title)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FINAL_ROOT / f"{sample_id} 3d reg.png", dpi=300)
    plt.close(fig)

    colors = {
        "DG": "#ff7f0e",
        "MH": "#1f77b4",
        "LH": "#d62728",
        "Ependymal": "#2ca02c",
    }
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(cells["ccfz"] / 25, cells["ccfy"] / 25, s=0.08, c="0.75")
    for label, cell_class in all_slice.TARGET_CLASSES.items():
        group = cells[cells["subclass_label_transfer"].eq(cell_class)]
        ax.scatter(
            group["ccfz"] / 25,
            group["ccfy"] / 25,
            s=2,
            c=colors[label],
            label=f"{label} cells",
        )
        contours = skimage.measure.find_contours(
            np.isin(atlas_slice, all_slice.TARGET_IDS[label]).astype(float), 0.5
        )
        for contour_index, contour in enumerate(contours):
            ax.plot(
                contour[:, 1],
                contour[:, 0],
                c=colors[label],
                linewidth=1.5,
                label=f"{label} atlas" if contour_index == 0 else None,
            )
    ax.set(xlim=(0, atlas_slice.shape[1]), ylim=(atlas_slice.shape[0], 0), title=title)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", markerscale=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FINAL_ROOT / f"{sample_id} regional targets.png", dpi=300)
    plt.close(fig)


def plot_requested_qc_figures(
    sample_id: str,
    slice_number: int,
    allen_slice: int,
    status: str,
    fixed: np.ndarray,
    moving: np.ndarray,
    moving_rigid: np.ndarray,
    moving_spline: np.ndarray,
    moving_borders: np.ndarray,
    moving_borders_spline: np.ndarray,
    cells: pd.DataFrame,
    regional,
) -> tuple[str, str]:
    percent_high = 0.99
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(15, 10))
    fig.suptitle(
        f"{sample_id} — FINAL {status}\nslice {slice_number} to Allen {allen_slice}",
        fontsize=16,
    )
    panels = [
        ("fixed slice", fixed, "gray"),
        (
            "fixed + rigid",
            regional.imageoverlay(
                regional.imagescPercent(regional.scale_result(fixed), 0, percent_high),
                regional.imagescPercent(
                    regional.scale_result(moving_rigid), 0, percent_high
                ),
            ),
            None,
        ),
        (
            "fixed + spline",
            regional.imageoverlay(
                regional.imagescPercent(regional.scale_result(fixed), 0, percent_high),
                regional.imagescPercent(
                    regional.scale_result(moving_spline), 0, percent_high
                ),
            ),
            None,
        ),
    ]
    for axis, (panel_title, image, cmap) in zip(axes[0, :3], panels):
        axis.set_title(panel_title)
        axis.imshow(image, cmap=cmap)
    axes[0, 3].set_title("fixed + borders spline")
    axes[0, 3].imshow(regional.imagescPercent(fixed, 0, percent_high), cmap="gray")
    axes[0, 3].imshow(
        regional.border_transparency(moving_borders_spline, RGBval=[1, 0, 1])
    )

    for axis, panel_title, image in (
        (axes[1, 0], f"moving Allen {allen_slice}", moving),
        (axes[1, 1], "moving rigid", moving_rigid),
        (axes[1, 2], "moving spline", moving_spline),
    ):
        axis.set_title(panel_title)
        axis.imshow(regional.imagescPercent(image, 0, percent_high), cmap="gray")
    axes[1, 3].axis("off")

    axes[2, 0].set_title("borders")
    axes[2, 0].set_facecolor("black")
    axes[2, 0].imshow(regional.border_transparency(moving_borders, RGBval=[1, 0, 1]))
    axes[2, 1].set_title("moving + borders")
    axes[2, 1].imshow(regional.imagescPercent(moving, 0, percent_high), cmap="gray")
    axes[2, 1].imshow(regional.border_transparency(moving_borders, RGBval=[1, 0, 1]))
    axes[2, 2].set_title("moving + borders spline")
    axes[2, 2].imshow(
        regional.imagescPercent(moving_spline, 0, percent_high), cmap="gray"
    )
    axes[2, 2].imshow(
        regional.border_transparency(moving_borders_spline, RGBval=[1, 0, 1])
    )
    axes[2, 3].axis("off")
    for axis in axes.flat:
        axis.set_xticks([])
        axis.set_yticks([])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    reg_name = (
        f"{sample_id}_slice_{slice_number:03d}_allen_{allen_slice:03d}_reg_5.jpg"
    )
    fig.savefig(FINAL_ROOT / reg_name, dpi=160, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15, 10))
    axes[0].set_title(f"FINAL Allen borders — {status}")
    axes[0].imshow(
        regional.border_transparency(moving_borders_spline, RGBval=[0, 0, 0])
    )
    axes[1].set_title(f"FINAL cells + Allen borders — {status}")
    axes[1].set_facecolor("white")
    axes[1].scatter(
        cells["fixed_x"],
        cells["fixed_y"],
        s=0.3,
        c=cells["subclass_color"],
    )
    axes[1].imshow(
        regional.border_transparency(moving_borders_spline, RGBval=[0, 0, 0])
    )
    for axis in axes:
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    fig.suptitle(
        f"{sample_id} — FINAL {status} — slice {slice_number} to Allen {allen_slice}"
    )
    fig.tight_layout()
    cells_name = (
        f"{sample_id}_slice_{slice_number:03d}_allen_{allen_slice:03d}"
        "_reg_5_cells_all.jpg"
    )
    fig.savefig(FINAL_ROOT / cells_name, dpi=300, facecolor="white")
    plt.close(fig)
    return reg_name, cells_name


def make_library(records: list[dict]) -> None:
    width, height = 850, 560
    title_height = 70
    canvas = Image.new(
        "RGB", (2 * width, title_height + len(records) * height), "white"
    )
    draw = ImageDraw.Draw(canvas)
    draw.text((20, 22), "FINAL registration diagnostics", fill="black")
    draw.text((width + 20, 22), "FINAL cells + Allen borders", fill="black")
    for row_index, record in enumerate(records):
        y0 = title_height + row_index * height
        for column_index, key in enumerate(("registration_figure", "cells_figure")):
            image = trim_white(Image.open(FINAL_ROOT / record[key]))
            image.thumbnail((width - 30, height - 55))
            x = column_index * width + (width - image.width) // 2
            y = y0 + 42 + (height - 48 - image.height) // 2
            canvas.paste(image, (x, y))
        draw.text(
            (10, y0 + 8),
            f"{record['sample_id']} — {record['status']}",
            fill="black",
        )
    canvas.save(FINAL_ROOT / "FINAL_FITS_LIBRARY.jpg", quality=92)


def remove_root_legacy_images() -> list[str]:
    removed = []
    for extension in ("*.png", "*.jpg", "*.jpeg"):
        for path in FIGURE_ROOT.glob(extension):
            path.unlink()
            removed.append(path.name)
    return sorted(removed)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cleanup-root",
        action="store_true",
        help="Remove loose legacy images from the optimized root after verification.",
    )
    args = parser.parse_args()

    all_slice = load_all_slice_module()
    regional = all_slice.load_regional_module()
    decisions = pd.read_csv(FINAL_DECISIONS).sort_values("sample_id")
    baseline_cells = all_slice.load_baseline_cells(set(decisions["sample_id"]))
    eligibility = pd.read_csv(all_slice.ELIGIBILITY_CSV).set_index("sample_id")

    annotation_25, _ = nrrd.read(
        str(REFERENCE_ROOT / "annotation_25.nrrd"), index_order="F"
    )
    annotation_25 = annotation_25[:, :, :236]
    annotation_10, _ = nrrd.read(
        str(REFERENCE_ROOT / "annotation_10.nrrd"), index_order="F"
    )
    nissl, _ = nrrd.read(
        str(REFERENCE_ROOT / "ara_nissl_25.nrrd"), index_order="F"
    )
    nissl = nissl[:, :, :236]
    borders = skimage.io.imread(
        str(REFERENCE_ROOT / "annotations_25_border_only.tif")
    )[:, :, :236]

    if FINAL_ROOT.exists():
        shutil.rmtree(FINAL_ROOT)
    FINAL_ROOT.mkdir(parents=True)
    WORK_ROOT.mkdir(parents=True, exist_ok=True)

    records = []
    for decision in decisions.itertuples():
        sample_id = decision.sample_id
        accepted = str(decision.final_decision).startswith("accept")
        iteration = int(decision.accepted_iteration) if accepted else 3
        weight = float(decision.accepted_weight) if accepted else 0.0
        sheet_name = "landmarks003" if accepted else "baseline_copy"
        status = (
            f"accepted landmark fit (iteration {iteration:02d})"
            if accepted
            else "retained baseline fit (iteration 03)"
        )
        df = pd.read_excel(WORKBOOK, sheet_name=sheet_name).reset_index(drop=True)
        num = int(df.index[df["cell_metadata"].eq(sample_id)][0])
        row = df.iloc[num]
        final_cells = final_coordinates(sample_id, iteration, baseline_cells)
        plot_final_coordinate_figures(
            sample_id,
            status,
            final_cells,
            annotation_25,
            annotation_10,
            all_slice,
            regional,
        )

        work_dir = WORK_ROOT / sample_id
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True)
        prior_cwd = Path.cwd()
        os.chdir(work_dir)
        try:
            fixed = regional.modify_dapi(
                df,
                num,
                baseline_cells,
                cell_types=all_slice.parse_list(row["cell_types"]),
                factor=float(row["dapi_enhance_factor"]),
            )
            moving_annotation = annotation_25[int(row["allen_slice_num"])]
            moving_borders = borders[int(row["allen_slice_num"])]
            moving = regional.modify_nissl(
                nissl[int(row["allen_slice_num"])],
                moving_annotation,
                factor=float(row["nissl_enhance_factor"]),
                annot_dict=all_slice.parse_list(row["annots_to_amplify"]),
            )
            params_affine, params_spline = regional.params_from_df(df, num)
            pad_width = (
                int(row["pad_width"]) if pd.notna(row["pad_width"]) else 20
            )
            fixed_bbox, fixed, _ = regional.crop_and_pad_image(
                fixed,
                pad_width=pad_width,
                area_thresh=float(row["area_thresh"]),
            )
            _, _, moving_borders = regional.crop_and_pad_image(
                moving,
                moving_borders,
                pad_width=pad_width,
                area_thresh=float(row["area_thresh"]),
            )
            _, moving, moving_annotation = regional.crop_and_pad_image(
                moving,
                moving_annotation,
                pad_width=pad_width,
                area_thresh=float(row["area_thresh"]),
            )
            fixed_cells = regional.get_cell_metadata_for_slice_index(
                df,
                num,
                baseline_cells,
                ccf_pixel_size=25,
                bbox=fixed_bbox,
                pad_width=pad_width,
            )
            if accepted:
                eligibility_row = eligibility.loc[sample_id]
                enabled = {
                    label: bool(eligibility_row[f"{label}_eligible"])
                    for label in ("DG", "MH", "LH")
                }
                fixed_points, moving_points, diagnostics = (
                    all_slice.make_robust_landmarks(
                        fixed_cells, moving_annotation, enabled
                    )
                )
                if not len(fixed_points):
                    raise RuntimeError(f"No final landmarks generated for {sample_id}")
                regional.write_pts_file(fixed_points, name="fix.pts")
                regional.write_pts_file(moving_points, name="mov.pts")
                for parameter_map in (params_affine, params_spline):
                    parameter_map["Registration"] = regional.L2P(
                        ["MultiMetricMultiResolutionRegistration"]
                    )
                    parameter_map["Metric"] = regional.L2P(
                        [
                            "AdvancedMattesMutualInformation",
                            "CorrespondingPointsEuclideanDistanceMetric",
                        ]
                    )
                    parameter_map["Metric0Weight"] = regional.L2P([1 - weight])
                    parameter_map["Metric1Weight"] = regional.L2P([weight])
            transform, moving_spline = regional.register_images(
                fixed, moving, params_affine, params_spline
            )
            moving_rigid = regional.transform_image(moving, transform[0])
            moving_borders_spline = regional.transform_image(
                moving_borders, transform, interpolation=True
            )
        finally:
            os.chdir(prior_cwd)

        reg_name, cells_name = plot_requested_qc_figures(
            sample_id,
            int(row["Slice"]),
            int(row["allen_slice_num"]),
            status,
            fixed,
            moving,
            moving_rigid,
            moving_spline,
            moving_borders,
            moving_borders_spline,
            fixed_cells,
            regional,
        )
        records.append(
            {
                "sample_id": sample_id,
                "slice": int(row["Slice"]),
                "allen_slice_num": int(row["allen_slice_num"]),
                "status": status,
                "source_iteration": iteration,
                "registration_figure": reg_name,
                "cells_figure": cells_name,
            }
        )
        print(f"Rendered {sample_id}: {status}", flush=True)

    records_frame = pd.DataFrame(records)
    records_frame.to_csv(FINAL_ROOT / "FINAL_FIGURE_MANIFEST.csv", index=False)
    make_library(records)
    readme_lines = [
        "# Allen CCF registration — current final figures",
        "",
        "This folder contains only the final selected fit for each of the 13 slices.",
        "Every slice has five freshly rendered figures:",
        "",
        "1. `2d reg.png` — cells with fixed-AP Allen contours",
        "2. `3d reg.png` — cells with curved-surface Allen contours",
        "3. `regional targets.png` — DG, MH, LH, and ependymal diagnostics",
        "4. `reg_5.jpg` — full fixed/rigid/spline registration diagnostic",
        "5. `reg_5_cells_all.jpg` — Allen borders alone and over all cells",
        "",
        "Accepted regional fits: roi1_run2, roi2_run2, roi4_run3 (iteration 11).",
        "All other slices use the retained safe iteration-03 baseline.",
        "",
        "Historical iteration directories outside this folder are audit records, not finals.",
    ]
    (FINAL_ROOT / "README.md").write_text("\n".join(readme_lines) + "\n")

    expected = len(records) * 5 + 3
    actual = len([path for path in FINAL_ROOT.iterdir() if path.is_file()])
    if actual != expected:
        raise RuntimeError(f"Final folder verification failed: {actual} != {expected}")

    removed = remove_root_legacy_images() if args.cleanup_root else []
    (FIGURE_ROOT / "CURRENT_FINAL_FIGURES_ARE_IN_final.txt").write_text(
        "Use figures/Allen_CCF_alignment_optimized/final/ for current final fits.\n"
        "iteration_* directories are historical audit records.\n"
    )
    (FINAL_ROOT / "render_summary.json").write_text(
        json.dumps(
            {
                "slice_count": len(records),
                "figures_per_slice": 5,
                "removed_legacy_root_images": removed,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
