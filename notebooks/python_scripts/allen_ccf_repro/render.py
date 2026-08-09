"""Final figures, mask and landmark QC overlays, and the library montage."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skimage.measure
import skimage.segmentation
from PIL import Image, ImageChops, ImageDraw

from . import atlas_io, config, elastix_stage, nb_functions

REGIONAL_COLOURS = {
    "DG_granule": "#ff7f0e",
    "MH": "#1f77b4",
    "LH": "#d62728",
    "ventricle_walls": "#2ca02c",
    "CA1_pyramidal": "#9467bd",
    "CA3_pyramidal": "#8c564b",
    "caudoputamen": "#17becf",
}


def atlas_outlines(annotation_slice: np.ndarray) -> list:
    outlines = []
    for value in np.unique(annotation_slice):
        if value == 0:
            continue
        outlines.extend(
            skimage.measure.find_contours((annotation_slice == value).astype(float), 0.99)
        )
    return outlines


def sample_atlas_surface(cells: pd.DataFrame, annotation_10: np.ndarray) -> np.ndarray:
    """Build the 10 um annotation surface that each cell's AP level points at."""
    valid = cells.dropna(subset=["ccfx", "ccfy", "ccfz"])
    shape = annotation_10[0].shape
    surface = np.zeros(shape, dtype=np.uint16)
    rows = np.rint(valid["ccfy"].to_numpy() / 10).astype(int)
    columns = np.rint(valid["ccfz"].to_numpy() / 10).astype(int)
    levels = np.rint(valid["ccfx"].to_numpy() / 10).astype(int)
    inside = (
        (rows >= 0)
        & (rows < shape[0])
        & (columns >= 0)
        & (columns < shape[1])
        & (levels >= 0)
        & (levels < annotation_10.shape[0])
    )
    surface[rows[inside], columns[inside]] = levels[inside]
    surface = skimage.segmentation.expand_labels(surface, distance=30)
    row_grid, column_grid = np.mgrid[0 : shape[0], 0 : shape[1]]
    return annotation_10[surface, row_grid, column_grid]


def _frame(x: np.ndarray, y: np.ndarray, mask: np.ndarray, margin: int = 15):
    """Axis limits and a matching figure size around the drawn content.

    The atlas plane keeps a wide empty surround, so framing on the content
    rather than the array bounds is what stops the panels from being mostly
    background.
    """
    rows, columns = np.nonzero(mask)
    x_min = min(np.nanmin(x), columns.min()) - margin if columns.size else np.nanmin(x) - margin
    x_max = max(np.nanmax(x), columns.max()) + margin if columns.size else np.nanmax(x) + margin
    y_min = min(np.nanmin(y), rows.min()) - margin if rows.size else np.nanmin(y) - margin
    y_max = max(np.nanmax(y), rows.max()) + margin if rows.size else np.nanmax(y) + margin
    width, height = x_max - x_min, y_max - y_min
    scale = 9.0 / max(width, height)
    return (x_min, x_max), (y_max, y_min), (max(width * scale, 3.0), max(height * scale, 3.0) + 0.6)


def coordinate_figures(
    output_root: Path,
    sample_id: str,
    title: str,
    cells: pd.DataFrame,
    annotation_plane: np.ndarray,
    annotation_10: np.ndarray,
    registry: dict,
    allen_name_to_annots: dict,
) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    valid = cells.dropna(subset=["ccfy", "ccfz"])

    xlim, ylim, size = _frame(valid["ccfz"] / 25, valid["ccfy"] / 25, annotation_plane > 0)
    figure, axis = plt.subplots(figsize=size)
    axis.scatter(valid["ccfz"] / 25, valid["ccfy"] / 25, s=0.08, c=valid["subclass_color"])
    for contour in atlas_outlines(annotation_plane):
        axis.plot(contour[:, 1], contour[:, 0], "k", linewidth=0.35)
    axis.set(xlim=xlim, ylim=ylim, title=title)
    axis.set_aspect("equal")
    axis.axis("off")
    figure.tight_layout()
    figure.savefig(output_root / f"{sample_id} 2d reg.png", dpi=300)
    plt.close(figure)

    sampled = sample_atlas_surface(valid, annotation_10)
    xlim, ylim, size = _frame(valid["ccfz"] / 10, valid["ccfy"] / 10, sampled > 0, margin=40)
    figure, axis = plt.subplots(figsize=size)
    axis.scatter(valid["ccfz"] / 10, valid["ccfy"] / 10, s=0.5, c=valid["subclass_color"])
    for contour in atlas_outlines(sampled):
        axis.plot(contour[:, 1], contour[:, 0], "k", linewidth=1)
    axis.set(xlim=xlim, ylim=ylim, title=title)
    axis.set_aspect("equal")
    axis.axis("off")
    figure.tight_layout()
    figure.savefig(output_root / f"{sample_id} 3d reg.png", dpi=300)
    plt.close(figure)

    xlim, ylim, size = _frame(valid["ccfz"] / 25, valid["ccfy"] / 25, annotation_plane > 0)
    figure, axis = plt.subplots(figsize=size)
    axis.scatter(valid["ccfz"] / 25, valid["ccfy"] / 25, s=0.08, c="0.78")
    for name, colour in REGIONAL_COLOURS.items():
        entry = registry["regions"].get(name)
        if entry is None:
            continue
        group = valid[valid["subclass_label_transfer"].isin(entry.get("cell_classes", []))]
        if len(group):
            axis.scatter(group["ccfz"] / 25, group["ccfy"] / 25, s=2, c=colour, label=f"{name} cells")
        ids = config.resolve_region_ids(entry, allen_name_to_annots)
        mask = atlas_io.region_mask(annotation_plane, ids)
        if not mask.any():
            continue
        for index, contour in enumerate(skimage.measure.find_contours(mask.astype(float), 0.5)):
            axis.plot(
                contour[:, 1],
                contour[:, 0],
                c=colour,
                linewidth=1.4,
                label=f"{name} atlas" if index == 0 else None,
            )
    axis.set(xlim=xlim, ylim=ylim, title=title)
    axis.set_aspect("equal")
    axis.axis("off")
    axis.legend(loc="lower right", markerscale=4, fontsize=7)
    figure.tight_layout()
    figure.savefig(output_root / f"{sample_id} regional targets.png", dpi=300)
    plt.close(figure)


def registration_qc_figures(
    output_root: Path,
    sample_id: str,
    slice_number: int,
    allen_slice: int,
    status: str,
    fixed: np.ndarray,
    moving: np.ndarray,
    moving_affine: np.ndarray,
    moving_spline: np.ndarray,
    borders: np.ndarray,
    borders_spline: np.ndarray,
    cells: pd.DataFrame,
) -> tuple:
    percent_high = 0.99
    figure, axes = plt.subplots(nrows=3, ncols=4, figsize=(15, 10))
    figure.suptitle(f"{sample_id} — {status}\nslice {slice_number} to Allen {allen_slice}", fontsize=15)

    axes[0, 0].set_title("fixed slice")
    axes[0, 0].imshow(fixed, cmap="gray")
    axes[0, 1].set_title("fixed + affine")
    axes[0, 1].imshow(
        nb_functions.imageoverlay(
            nb_functions.imagescPercent(nb_functions.scale_result(fixed), 0, percent_high),
            nb_functions.imagescPercent(nb_functions.scale_result(moving_affine), 0, percent_high),
        )
    )
    axes[0, 2].set_title("fixed + spline")
    axes[0, 2].imshow(
        nb_functions.imageoverlay(
            nb_functions.imagescPercent(nb_functions.scale_result(fixed), 0, percent_high),
            nb_functions.imagescPercent(nb_functions.scale_result(moving_spline), 0, percent_high),
        )
    )
    axes[0, 3].set_title("fixed + borders spline")
    axes[0, 3].imshow(nb_functions.imagescPercent(fixed, 0, percent_high), cmap="gray")
    axes[0, 3].imshow(nb_functions.border_transparency(borders_spline, RGBval=[1, 0, 1]))

    for axis, panel_title, image in (
        (axes[1, 0], f"moving Allen {allen_slice}", moving),
        (axes[1, 1], "moving affine", moving_affine),
        (axes[1, 2], "moving spline", moving_spline),
    ):
        axis.set_title(panel_title)
        axis.imshow(nb_functions.imagescPercent(image, 0, percent_high), cmap="gray")
    axes[1, 3].axis("off")

    axes[2, 0].set_title("borders")
    axes[2, 0].set_facecolor("black")
    axes[2, 0].imshow(nb_functions.border_transparency(borders, RGBval=[1, 0, 1]))
    axes[2, 1].set_title("moving + borders")
    axes[2, 1].imshow(nb_functions.imagescPercent(moving, 0, percent_high), cmap="gray")
    axes[2, 1].imshow(nb_functions.border_transparency(borders, RGBval=[1, 0, 1]))
    axes[2, 2].set_title("moving + borders spline")
    axes[2, 2].imshow(nb_functions.imagescPercent(moving_spline, 0, percent_high), cmap="gray")
    axes[2, 2].imshow(nb_functions.border_transparency(borders_spline, RGBval=[1, 0, 1]))
    axes[2, 3].axis("off")
    for axis in axes.flat:
        axis.set_xticks([])
        axis.set_yticks([])
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    registration_name = f"{sample_id}_slice_{slice_number:03d}_allen_{allen_slice:03d}_reg_5.jpg"
    figure.savefig(output_root / registration_name, dpi=160, facecolor="white")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(15, 10))
    axes[0].set_title(f"Allen borders — {status}")
    axes[0].imshow(nb_functions.border_transparency(borders_spline, RGBval=[0, 0, 0]))
    axes[1].set_title(f"cells + Allen borders — {status}")
    axes[1].set_facecolor("white")
    axes[1].scatter(cells["fixed_x"], cells["fixed_y"], s=0.3, c=cells["subclass_color"])
    axes[1].imshow(nb_functions.border_transparency(borders_spline, RGBval=[0, 0, 0]))
    for axis in axes:
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    figure.suptitle(f"{sample_id} — {status} — slice {slice_number} to Allen {allen_slice}")
    figure.tight_layout()
    cells_name = f"{sample_id}_slice_{slice_number:03d}_allen_{allen_slice:03d}_reg_5_cells_all.jpg"
    figure.savefig(output_root / cells_name, dpi=300, facecolor="white")
    plt.close(figure)
    return registration_name, cells_name


def mask_qc_figure(
    output_root: Path,
    sample_id: str,
    fixed: np.ndarray,
    tissue: np.ndarray,
    cell_supported: np.ndarray,
    exclusion: np.ndarray,
    artifact: np.ndarray,
    fixed_mask: np.ndarray,
    title: str,
) -> str:
    figure, axes = plt.subplots(1, 3, figsize=(16, 6))
    axes[0].set_title("fixed image")
    axes[0].imshow(nb_functions.imagescPercent(fixed, 0, 0.99), cmap="gray")

    overlay = np.zeros(fixed.shape + (4,))
    overlay[tissue] = [0.2, 0.6, 1.0, 0.28]
    overlay[cell_supported] = [0.1, 0.9, 0.3, 0.55]
    overlay[artifact] = [1.0, 0.85, 0.1, 0.65]
    overlay[exclusion] = [1.0, 0.1, 0.1, 0.55]
    axes[1].set_title("tissue (blue), cell-supported (green),\nartifact (yellow), missing tissue (red)")
    axes[1].imshow(nb_functions.imagescPercent(fixed, 0, 0.99), cmap="gray")
    axes[1].imshow(overlay)

    axes[2].set_title("elastix fixed mask")
    axes[2].imshow(fixed_mask, cmap="gray")
    for axis in axes:
        axis.set_xticks([])
        axis.set_yticks([])
    figure.suptitle(title)
    figure.tight_layout()
    name = f"{sample_id} mask QC.png"
    figure.savefig(output_root / name, dpi=170, facecolor="white")
    plt.close(figure)
    return name


def landmark_qc_figure(
    output_root: Path,
    sample_id: str,
    fixed: np.ndarray,
    moving: np.ndarray,
    landmarks: dict,
    title: str,
) -> str:
    figure, axes = plt.subplots(1, 2, figsize=(14, 7))
    axes[0].set_title("fixed landmarks")
    axes[0].imshow(nb_functions.imagescPercent(fixed, 0, 0.99), cmap="gray")
    axes[1].set_title("moving (Allen) landmarks")
    axes[1].imshow(nb_functions.imagescPercent(moving, 0, 0.99), cmap="gray")
    fixed_points = landmarks["fixed"]
    moving_points = landmarks["moving"]
    if len(fixed_points):
        axes[0].scatter(fixed_points[:, 0], fixed_points[:, 1], s=26, c="#ff2d55", marker="x")
        axes[1].scatter(moving_points[:, 0], moving_points[:, 1], s=26, c="#ff2d55", marker="x")
        for index, (point, pair) in enumerate(zip(fixed_points, moving_points)):
            axes[0].annotate(str(index), point, color="#ffd60a", fontsize=6)
            axes[1].annotate(str(index), pair, color="#ffd60a", fontsize=6)
    for axis in axes:
        axis.set_xticks([])
        axis.set_yticks([])
    figure.suptitle(title)
    figure.tight_layout()
    name = f"{sample_id} landmark QC.png"
    figure.savefig(output_root / name, dpi=170, facecolor="white")
    plt.close(figure)
    return name


def trim_white(image: Image.Image) -> Image.Image:
    background = Image.new(image.mode, image.size, (255, 255, 255))
    difference = ImageChops.difference(image, background)
    box = difference.getbbox()
    return image.crop(box) if box else image


def library(output_root: Path, records: list, name: str = "FINAL_FITS_LIBRARY.jpg") -> Path:
    width, height = 780, 560
    columns = 3
    title_height = 70
    canvas = Image.new("RGB", (columns * width, title_height + len(records) * height), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((20, 22), "Registration diagnostics", fill="black")
    draw.text((width + 20, 22), "Cells + Allen borders", fill="black")
    draw.text((2 * width + 20, 22), "3D CCF registration", fill="black")
    for row_index, record in enumerate(records):
        top = title_height + row_index * height
        panel_keys = (
            ("registration_figure", None),
            ("cells_figure", None),
            ("surface_figure", f"{record['sample_id']} 3d reg.png"),
        )
        for column_index, (key, fallback) in enumerate(panel_keys):
            path_name = record.get(key) or fallback
            if not path_name:
                continue
            path = output_root / path_name
            if not path.exists():
                continue
            image = trim_white(Image.open(path).convert("RGB"))
            image.thumbnail((width - 30, height - 55))
            left = column_index * width + (width - image.width) // 2
            canvas.paste(image, (left, top + 42 + (height - 48 - image.height) // 2))
        draw.text((10, top + 8), f"{record['sample_id']} — {record['status']}", fill="black")
    path = output_root / name
    canvas.save(path, quality=92)
    return path


def warped_images(moving_ctx: dict, transform) -> dict:
    affine_only = tuple(transform)[: max(len(transform) - 1, 1)]
    return {
        "affine": elastix_stage.transform_image(moving_ctx["image"], affine_only),
        "spline": elastix_stage.transform_image(moving_ctx["image"], transform),
        "borders_spline": elastix_stage.transform_image(moving_ctx["borders"], transform),
    }
