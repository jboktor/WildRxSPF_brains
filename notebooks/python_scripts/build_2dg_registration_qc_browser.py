#!/usr/bin/env python3
"""Build a per-slice 2DG → Allen CCF registration QC browser.

The browser is a visual audit of every delivered coronal plane, including
continuity-rejected and blank sections. Each frame shows:

1. original Template anatomy;
2. source-space Intruder or Baseline WildR−SPF effect;
3. selected Allen Nissl plane with annotation borders;
4. warped Template anatomy with source/atlas tissue edges;
5. the selected effect projected onto the Allen plane; and
6. sagittal AP context plus the full front-to-back registration trajectory.

This is display-only. It reuses the 2026-11-11 registration products and the
2026-11-12 Baseline decoding LUT; it does not refit Elastix.

Run with the spatialomics Python binary:

  /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
    notebooks/python_scripts/build_2dg_registration_qc_browser.py
"""

from __future__ import annotations

import base64
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nrrd
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from PIL import Image
from scipy import ndimage
from skimage.segmentation import find_boundaries


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
SOURCE_DATE = "2026-11-11"
CONTROL_DATE = "2026-11-12"
ANALYSIS_DATE = "2026-11-21c"

SOURCE_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / SOURCE_DATE
OUT_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
FIGURE_DIR = WKDIR / "figures/2DG_CCF_Integration"
INPUT_ROOT = WKDIR / "data/input/2DG_data"
ATLAS_ROOT = WKDIR / "data/input/allen_registration_ref"
TABLE_DIR = WKDIR / "tables/2DG_CCF_Integration"

REG_SUMMARY = SOURCE_DIR / f"2dg_registration_summary_strict_{SOURCE_DATE}.csv"
NISSL_PATH = ATLAS_ROOT / "ara_nissl_25.nrrd"
ANNOTATION_PATH = ATLAS_ROOT / "annotation_25.nrrd"
BASELINE_ROOT = INPUT_ROOT / "WilderR vs SPF _Baseline"
LUT_PATH = TABLE_DIR / "2DG_tiff_rgb_to_t_lut_validated_2026-11-07.csv"

HTML_OUT = FIGURE_DIR / f"2DG_registration_qc_browser_{ANALYSIS_DATE}.html"
MANIFEST_OUT = OUT_DIR / f"2DG_registration_qc_browser_manifest_{ANALYSIS_DATE}.json"
THUMBNAIL_OUT = FIGURE_DIR / f"2DG_registration_qc_browser_preview_{ANALYSIS_DATE}.png"

T_THRESHOLD = 1.725
T_CLIP = 3.5
CCF_VOXEL_UM = 25.0


def log(message: str) -> None:
    print(f"[2DG-registration-QC] {message}", flush=True)


def refuse_overwrite(*paths: Path) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing artifacts:\n" + "\n".join(existing))


def load_atlas() -> tuple[np.ndarray, np.ndarray]:
    nissl, _ = nrrd.read(str(NISSL_PATH))
    annotation, _ = nrrd.read(str(ANNOTATION_PATH))
    return nissl.astype(np.float32), annotation


def grayscale(array: np.ndarray, *, positive_only: bool = False) -> np.ndarray:
    values = np.asarray(array, dtype=np.float32)
    finite = np.isfinite(values)
    if positive_only:
        finite &= values > 0
    if not finite.any():
        return np.zeros(values.shape, dtype=np.float32)
    lo, hi = np.nanpercentile(values[finite], [1.5, 98.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(values[finite])), float(np.nanmax(values[finite]))
    if hi <= lo:
        return np.zeros(values.shape, dtype=np.float32)
    return np.clip((values - lo) / (hi - lo), 0, 1)


def decode_rgb_t(path: Path, lut_map: dict[int, float], shape: tuple[int, int]) -> np.ndarray:
    if not path.exists():
        return np.full(shape, np.nan, dtype=np.float32)
    rgb = np.asarray(Image.open(path).convert("RGB"), dtype=np.uint8)
    if rgb.shape[:2] != shape:
        raise ValueError(f"Baseline shape {rgb.shape[:2]} does not match Template {shape}: {path}")
    packed = (
        (rgb[..., 0].astype(np.int32) << 16)
        | (rgb[..., 1].astype(np.int32) << 8)
        | rgb[..., 2].astype(np.int32)
    )
    out = np.full(shape, np.nan, dtype=np.float32)
    for key, value in lut_map.items():
        out[packed == key] = value
    gray = (rgb[..., 0] == rgb[..., 1]) & (rgb[..., 1] == rgb[..., 2])
    out[gray] = np.nan
    return out


def load_lut() -> dict[int, float]:
    lut = pd.read_csv(LUT_PATH)
    red = lut["overlay_r"].to_numpy(dtype=np.int32)
    green = lut["overlay_g"].to_numpy(dtype=np.int32)
    blue = lut["overlay_b"].to_numpy(dtype=np.int32)
    packed = (
        (red << 16)
        | (green << 8)
        | blue
    )
    return dict(zip(packed, lut["t_value"].to_numpy(dtype=float)))


def rasterize_mean(
    values: np.ndarray,
    ys: np.ndarray,
    xs: np.ndarray,
    shape: tuple[int, int],
) -> np.ndarray:
    out = np.full(shape, np.nan, dtype=np.float32)
    sums = np.zeros(shape, dtype=np.float64)
    counts = np.zeros(shape, dtype=np.int32)
    yi = ys.astype(int)
    xi = xs.astype(int)
    values = np.asarray(values, dtype=float)
    ok = (
        np.isfinite(values)
        & (yi >= 0)
        & (yi < shape[0])
        & (xi >= 0)
        & (xi < shape[1])
    )
    np.add.at(sums, (yi[ok], xi[ok]), values[ok])
    np.add.at(counts, (yi[ok], xi[ok]), 1)
    occupied = counts > 0
    out[occupied] = (sums[occupied] / counts[occupied]).astype(np.float32)
    return out


def effect_overlay(base: np.ndarray, effect: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(effect) & (np.abs(effect) >= T_THRESHOLD)
    return base, np.ma.masked_where(~mask, effect)


def registration_edge_rgb(
    warped_template: np.ndarray,
    allen_tissue: np.ndarray,
) -> np.ndarray:
    present = np.isfinite(warped_template)
    weights = ndimage.gaussian_filter(present.astype(np.float32), sigma=0.8)
    smoothed = ndimage.gaussian_filter(
        np.nan_to_num(warped_template, nan=0.0),
        sigma=0.8,
    )
    reconstructed = np.full(warped_template.shape, np.nan, dtype=np.float32)
    usable = weights > 0.04
    reconstructed[usable] = smoothed[usable] / weights[usable]
    base = grayscale(reconstructed)
    rgb = np.stack([base, base, base], axis=-1)
    source_mask = ndimage.binary_dilation(present, iterations=2)
    source_mask = ndimage.binary_closing(source_mask, iterations=3)
    source_mask = ndimage.binary_fill_holes(source_mask)
    source_edge = find_boundaries(source_mask, mode="outer")
    atlas_edge = find_boundaries(allen_tissue.astype(bool), mode="outer")
    # Cyan = warped source tissue; orange = Allen tissue.
    rgb[source_edge] = np.array([0.00, 0.72, 0.78])
    rgb[atlas_edge] = np.array([1.00, 0.42, 0.05])
    overlap = source_edge & atlas_edge
    rgb[overlap] = np.array([0.10, 0.72, 0.25])
    return rgb


def sagittal_context(nissl: np.ndarray, ap_index: int) -> np.ndarray:
    mid_ml = nissl.shape[2] // 2
    sagittal = np.log1p(np.maximum(nissl[:, :, mid_ml], 0)).T
    rgb = np.stack([grayscale(sagittal, positive_only=True)] * 3, axis=-1)
    ap = int(np.clip(ap_index, 0, rgb.shape[1] - 1))
    rgb[:, max(0, ap - 1) : ap + 2] = np.array([1.0, 0.42, 0.05])
    return rgb


def qc_status(row: pd.Series) -> tuple[str, str]:
    if bool(row.get("accepted_final", False)):
        return "accepted", "Accepted"
    reason = str(row.get("rejection_reason", "") or "").strip()
    continuity = str(row.get("continuity_note", "") or "").strip()
    if reason == "blank_or_empty_template":
        return "blank", "Blank input"
    if bool(row.get("accepted_slice_fit", False)):
        return "continuity", "Continuity reject"
    detail = reason.replace("_", " ") if reason else continuity.replace("_", " ")
    return "fit-failed", detail.title() if detail else "Fit rejected"


def build_frame(
    row: pd.Series,
    map_key: str,
    nissl: np.ndarray,
    annotation: np.ndarray,
    reg: pd.DataFrame,
    lut_map: dict[int, float],
) -> bytes:
    slice_position = int(row["slice_position"])
    source_number = int(row["source_number"])
    ap_index = int(row["selected_allen_index"])
    prepared_path = SOURCE_DIR / "prepared" / f"slice_{slice_position:03d}.npz"
    pixels_path = SOURCE_DIR / "work" / f"slice_{slice_position:03d}" / "pixels_ccf.csv.gz"

    prepared = np.load(prepared_path)
    template = prepared["template_y"].astype(np.float32)
    template_2dg = prepared["template_2dg"].astype(np.float32)
    tissue = prepared["tissue"].astype(bool)
    intruder_t = prepared["t_value"].astype(np.float32)
    baseline_path = BASELINE_ROOT / f"template42{source_number:03d}.tif"
    baseline_available = baseline_path.exists()
    baseline_t = decode_rgb_t(baseline_path, lut_map, template.shape)
    source_effect = intruder_t if map_key == "intruder" else baseline_t

    allen_plane = nissl[ap_index].astype(np.float32)
    allen_tissue = annotation[ap_index] > 0
    atlas_shape = allen_plane.shape
    warped_template = np.full(atlas_shape, np.nan, dtype=np.float32)
    warped_effect = np.full(atlas_shape, np.nan, dtype=np.float32)

    if pixels_path.exists():
        pixels = pd.read_csv(
            pixels_path,
            usecols=[
                "source_x",
                "source_y",
                "atlas_x",
                "atlas_y",
                "template_2dg",
                "t_value",
            ],
        )
        atlas_y = pixels["atlas_y"].to_numpy()
        atlas_x = pixels["atlas_x"].to_numpy()
        warped_template = rasterize_mean(
            pixels["template_2dg"].to_numpy(),
            atlas_y,
            atlas_x,
            atlas_shape,
        )
        if map_key == "intruder":
            values = pixels["t_value"].to_numpy()
        else:
            sy = pixels["source_y"].to_numpy(dtype=int)
            sx = pixels["source_x"].to_numpy(dtype=int)
            values = baseline_t[sy, sx]
        warped_effect = rasterize_mean(values, atlas_y, atlas_x, atlas_shape)

    template_gray = grayscale(template, positive_only=True)
    allen_gray = grayscale(np.log1p(np.maximum(allen_plane, 0)), positive_only=True)
    allen_with_border = np.stack([allen_gray] * 3, axis=-1)
    annotation_edge = find_boundaries(allen_tissue, mode="outer")
    allen_with_border[annotation_edge] = np.array([0.18, 0.18, 0.18])
    edge_overlay = registration_edge_rgb(warped_template, allen_tissue)
    sag = sagittal_context(nissl, ap_index)

    status_key, status_label = qc_status(row)
    color = {
        "accepted": "#20966d",
        "continuity": "#e36f14",
        "blank": "#777777",
        "fit-failed": "#b33b3b",
    }[status_key]
    map_label = "Intruder WildR−SPF" if map_key == "intruder" else "Baseline WildR−SPF"

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 6.9), facecolor="white")
    for ax in axes.flat:
        ax.set_axis_off()
        ax.set_facecolor("white")

    axes[0, 0].imshow(template_gray, cmap="gray", vmin=0, vmax=1)
    axes[0, 0].contour(tissue, levels=[0.5], colors=["#00a7b5"], linewidths=0.7)
    axes[0, 0].set_title("1 · Original Template anatomy", fontsize=10.5, weight="bold")

    source_base, source_masked = effect_overlay(template_gray, source_effect)
    axes[0, 1].imshow(source_base, cmap="gray", vmin=0, vmax=1)
    axes[0, 1].imshow(
        source_masked,
        cmap="RdBu_r",
        norm=Normalize(-T_CLIP, T_CLIP),
        alpha=0.82,
    )
    axes[0, 1].set_title(f"2 · Source {map_label} t", fontsize=10.5, weight="bold")

    axes[0, 2].imshow(allen_with_border)
    axes[0, 2].set_title(
        f"3 · Selected Allen plane · AP {ap_index} · {float(row['selected_bregma_mm']):+.2f} mm",
        fontsize=10.5,
        weight="bold",
    )

    axes[1, 0].imshow(edge_overlay)
    axes[1, 0].set_title(
        "4 · Registration edge agreement\ncyan=warped Template · orange=Allen · green=agreement",
        fontsize=9.4,
        weight="bold",
    )

    registered_base, registered_masked = effect_overlay(allen_gray, warped_effect)
    axes[1, 1].imshow(registered_base, cmap="gray", vmin=0, vmax=1)
    axes[1, 1].imshow(
        registered_masked,
        cmap="RdBu_r",
        norm=Normalize(-T_CLIP, T_CLIP),
        alpha=0.82,
    )
    axes[1, 1].contour(allen_tissue, levels=[0.5], colors=["#333333"], linewidths=0.35)
    axes[1, 1].set_title(
        f"5 · {map_label} projected into CCF",
        fontsize=10.5,
        weight="bold",
    )
    if map_key == "baseline" and not baseline_available:
        for ax in (axes[0, 1], axes[1, 1]):
            ax.text(
                0.5,
                0.5,
                "Baseline TIFF unavailable\n(no effect data shown)",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=12,
                weight="bold",
                color="#8c2f1b",
                bbox={
                    "boxstyle": "round,pad=0.55",
                    "facecolor": "#fff6f2",
                    "edgecolor": "#d98b73",
                    "alpha": 0.96,
                },
            )

    axes[1, 2].imshow(sag)
    axes[1, 2].set_title("6 · AP context and front→back trajectory", fontsize=10.5, weight="bold")
    trajectory = axes[1, 2].inset_axes([0.08, 0.03, 0.84, 0.30])
    accepted = reg["accepted_final"].astype(bool).to_numpy()
    trajectory.plot(
        reg["slice_position"],
        reg["selected_allen_index"],
        color="#b8bcc3",
        linewidth=1.1,
        zorder=1,
    )
    trajectory.scatter(
        reg.loc[accepted, "slice_position"],
        reg.loc[accepted, "selected_allen_index"],
        s=9,
        color="#20966d",
        label="accepted",
        zorder=2,
    )
    trajectory.scatter(
        reg.loc[~accepted, "slice_position"],
        reg.loc[~accepted, "selected_allen_index"],
        s=11,
        color="#e36f14",
        label="rejected",
        zorder=2,
    )
    trajectory.scatter(
        [slice_position],
        [ap_index],
        s=46,
        facecolors="white",
        edgecolors="#111111",
        linewidths=1.2,
        zorder=3,
    )
    trajectory.set_facecolor((1, 1, 1, 0.90))
    trajectory.tick_params(labelsize=6, length=2)
    trajectory.set_xlabel("source slice", fontsize=6, labelpad=1)
    trajectory.set_ylabel("Allen AP", fontsize=6, labelpad=1)
    trajectory.spines[["top", "right"]].set_visible(False)

    dice = float(row["dice"])
    inside = float(row["inside_allen_tissue_fraction"])
    fig.suptitle(
        f"Template {source_number:03d} · {map_label} · {status_label}",
        x=0.02,
        y=0.985,
        ha="left",
        fontsize=15,
        weight="bold",
        color="#111111",
    )
    fig.text(
        0.98,
        0.976,
        (
            f"Dice {dice:.3f}  ·  inside Allen tissue {inside:.1%}  ·  "
            f"Jacobian {float(row['jacobian_min']):.2f}–{float(row['jacobian_max']):.2f}"
            if np.isfinite(float(row["jacobian_min"]))
            else f"Dice {dice:.3f}  ·  {status_label}"
        ),
        ha="right",
        va="top",
        fontsize=9,
        color=color,
        weight="bold",
    )
    fig.text(
        0.02,
        0.015,
        (
            "Registration metrics describe anatomical alignment, not WildR−SPF evidence. "
            f"Effect overlay shows |t| ≥ {T_THRESHOLD}; blue = WildR↑, red = SPF↑."
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#4f5560",
    )
    fig.subplots_adjust(left=0.02, right=0.985, top=0.91, bottom=0.06, wspace=0.08, hspace=0.19)

    buf = io.BytesIO()
    fig.savefig(buf, format="webp", dpi=115, facecolor="white", pil_kwargs={"quality": 84})
    plt.close(fig)
    return buf.getvalue()


def encoded_webp(data: bytes) -> str:
    return "data:image/webp;base64," + base64.b64encode(data).decode("ascii")


def build_frames(
    nissl: np.ndarray,
    annotation: np.ndarray,
) -> tuple[list[dict], dict, bytes]:
    reg = pd.read_csv(REG_SUMMARY).sort_values("slice_position").reset_index(drop=True)
    lut_map = load_lut()
    frames: list[dict] = []
    preview: bytes | None = None

    for index, row in reg.iterrows():
        status_key, status_label = qc_status(row)
        source_number = int(row["source_number"])
        baseline_available = (
            BASELINE_ROOT / f"template42{source_number:03d}.tif"
        ).exists()
        images = {}
        for map_key in ("intruder", "baseline"):
            image = build_frame(row, map_key, nissl, annotation, reg, lut_map)
            images[map_key] = encoded_webp(image)
            if preview is None and bool(row["accepted_final"]) and map_key == "intruder":
                preview = image
        frames.append(
            {
                "index": int(index),
                "slice_position": int(row["slice_position"]),
                "source_number": source_number,
                "selected_allen_index": int(row["selected_allen_index"]),
                "selected_bregma_mm": float(row["selected_bregma_mm"]),
                "dice": float(row["dice"]),
                "inside_fraction": float(row["inside_allen_tissue_fraction"]),
                "accepted_final": bool(row["accepted_final"]),
                "status_key": status_key,
                "status_label": status_label,
                "rejection_reason": str(row.get("rejection_reason", "") or ""),
                "continuity_note": str(row.get("continuity_note", "") or ""),
                "baseline_available": baseline_available,
                "images": images,
            }
        )
        if (index + 1) % 10 == 0:
            log(f"Rendered {index + 1}/{len(reg)} sections")

    manifest = {
        "analysis_date": ANALYSIS_DATE,
        "source_registration": SOURCE_DATE,
        "control_maps": CONTROL_DATE,
        "n_slices": len(frames),
        "n_accepted": int(reg["accepted_final"].sum()),
        "n_rejected": int((~reg["accepted_final"].astype(bool)).sum()),
        "map_keys": ["intruder", "baseline"],
        "map_labels": {
            "intruder": "Intruder WildR−SPF",
            "baseline": "Baseline WildR−SPF",
        },
        "n_baseline_available": sum(
            frame["baseline_available"] for frame in frames
        ),
        "baseline_unavailable_templates": [
            frame["source_number"]
            for frame in frames
            if not frame["baseline_available"]
        ],
    }
    if preview is None:
        raise RuntimeError("No accepted section available for preview")
    return frames, manifest, preview


def write_html(frames: list[dict], manifest: dict, path: Path) -> None:
    payload = json.dumps({"manifest": manifest, "frames": frames}, separators=(",", ":"))
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>2DG → Allen CCF per-slice registration QC</title>
<style>
  :root {{
    --ink:#17191c; --muted:#606771; --line:#dfe3e8; --soft:#f5f6f8;
    --accent:#ff6c0c; --accepted:#16845d; --rejected:#d56412;
  }}
  * {{ box-sizing:border-box; }}
  html,body {{ margin:0; padding:0; background:#fff; color:var(--ink);
    font-family:Inter,Helvetica,Arial,sans-serif; }}
  .wrap {{ max-width:1540px; margin:0 auto; padding:14px 18px 18px; }}
  .header {{ display:flex; gap:18px; align-items:flex-start; justify-content:space-between; }}
  h1 {{ margin:0 0 5px; font-size:21px; letter-spacing:-0.015em; }}
  .sub {{ margin:0; max-width:980px; color:var(--muted); font-size:12.5px; line-height:1.45; }}
  .counts {{ display:flex; gap:6px; flex-wrap:wrap; justify-content:flex-end; }}
  .pill {{ border:1px solid var(--line); border-radius:999px; padding:5px 9px;
    font-size:11px; white-space:nowrap; background:#fff; }}
  .pill.accepted {{ color:var(--accepted); border-color:#b7dfd0; background:#f0faf6; }}
  .pill.rejected {{ color:var(--rejected); border-color:#f3d0b8; background:#fff7f1; }}
  .controls {{ position:sticky; top:0; z-index:10; margin:12px 0 9px;
    padding:9px 10px; background:rgba(255,255,255,.97); border:1px solid var(--line);
    border-radius:10px; box-shadow:0 3px 12px rgba(30,40,55,.06); }}
  .control-row {{ display:flex; align-items:center; gap:7px; flex-wrap:wrap; }}
  .control-row + .control-row {{ margin-top:7px; }}
  button,select {{ border:1px solid var(--line); border-radius:7px; background:#fff;
    color:var(--ink); padding:6px 10px; font-size:11.5px; cursor:pointer; }}
  button:hover,button.active {{ border-color:var(--accent); background:#fff5ee; }}
  button.active {{ color:#9d3b00; font-weight:700; }}
  #prev,#next {{ font-weight:700; }}
  input[type=range] {{ flex:1 1 360px; min-width:220px; accent-color:var(--accent); }}
  .meta {{ min-width:360px; color:var(--muted); font-size:11.5px; }}
  .stage {{ border:1px solid var(--line); border-radius:12px; background:#fff;
    overflow:hidden; box-shadow:0 3px 16px rgba(30,40,55,.06); }}
  .stage img {{ width:100%; height:auto; max-height:calc(100vh - 205px);
    min-height:420px; object-fit:contain; display:block; background:#fff; }}
  .hint {{ margin-top:7px; color:var(--muted); font-size:11px; }}
  @media (max-width:800px) {{
    .header {{ display:block; }}
    .counts {{ justify-content:flex-start; margin-top:8px; }}
    .meta {{ min-width:100%; }}
    .stage img {{ min-height:0; }}
  }}
</style>
</head>
<body>
<main class="wrap">
  <div class="header">
    <div>
      <h1>2DG → Allen CCF registration QC, section by section</h1>
      <p class="sub">Move front→back through every delivered coronal plane. Compare the original Template,
      source contrast, selected Allen plane, edge agreement, registered contrast, and AP trajectory in one view.
      Registration metrics assess anatomy only; they are not genotype-effect statistics.</p>
    </div>
    <div class="counts">
      <span class="pill">{manifest['n_slices']} delivered</span>
      <span class="pill accepted">{manifest['n_accepted']} accepted</span>
      <span class="pill rejected">{manifest['n_rejected']} rejected / blank</span>
    </div>
  </div>
  <section class="controls">
    <div class="control-row">
      <button id="prev" type="button">◀ Previous</button>
      <input id="slider" type="range" min="0" max="{len(frames)-1}" value="0"/>
      <button id="next" type="button">Next ▶</button>
      <select id="sectionSelect" aria-label="Jump to Template section"></select>
    </div>
    <div class="control-row">
      <strong style="font-size:11.5px">Map</strong>
      <button type="button" class="map-btn active" data-map="intruder">Intruder WildR−SPF</button>
      <button type="button" class="map-btn" data-map="baseline">Baseline WildR−SPF</button>
      <strong style="font-size:11.5px;margin-left:8px">Show</strong>
      <button type="button" class="filter-btn active" data-filter="all">All sections</button>
      <button type="button" class="filter-btn" data-filter="accepted">Accepted only</button>
      <button type="button" class="filter-btn" data-filter="rejected">Rejected / blank</button>
      <span class="meta" id="meta"></span>
    </div>
  </section>
  <section class="stage"><img id="view" alt="Six-panel 2DG registration QC view"/></section>
  <div class="hint">Keyboard: ←/→ moves one section · I = Intruder · B = Baseline · A = all sections.
  Orange trajectory points were withheld from downstream regional analyses.</div>
</main>
<script>
const DATA={payload};
let mapKey="intruder";
let filterKey="all";
let visible=DATA.frames.map((_,i)=>i);
let position=0;
const view=document.getElementById("view");
const slider=document.getElementById("slider");
const select=document.getElementById("sectionSelect");
const meta=document.getElementById("meta");

function rebuildVisible(keepCurrent=true) {{
  const current=visible[position];
  visible=DATA.frames.map((fr,i)=>[fr,i]).filter(([fr]) =>
    filterKey==="all" || (filterKey==="accepted" ? fr.accepted_final : !fr.accepted_final)
  ).map(([,i])=>i);
  position=Math.max(0, keepCurrent ? visible.indexOf(current) : 0);
  if(position<0) position=0;
  slider.max=String(Math.max(visible.length-1,0));
  select.innerHTML=visible.map((i,j)=>{{
    const fr=DATA.frames[i];
    return `<option value="${{j}}">Template ${{String(fr.source_number).padStart(3,"0")}} · ${{fr.status_label}}</option>`;
  }}).join("");
  render();
}}
function render() {{
  if(!visible.length) return;
  const fr=DATA.frames[visible[position]];
  view.src=fr.images[mapKey];
  slider.value=String(position);
  select.value=String(position);
  const baselineNote=(mapKey==="baseline"&&!fr.baseline_available)?" · Baseline TIFF unavailable":"";
  meta.textContent=`Template ${{String(fr.source_number).padStart(3,"0")}} · ${{fr.status_label}} · Allen AP ${{fr.selected_allen_index}} · bregma ${{fr.selected_bregma_mm.toFixed(2)}} mm · ${{position+1}}/${{visible.length}}${{baselineNote}}`;
  document.querySelectorAll(".map-btn").forEach(b=>b.classList.toggle("active",b.dataset.map===mapKey));
  document.querySelectorAll(".filter-btn").forEach(b=>b.classList.toggle("active",b.dataset.filter===filterKey));
}}
function move(delta) {{
  position=Math.max(0,Math.min(visible.length-1,position+delta));
  render();
}}
document.getElementById("prev").onclick=()=>move(-1);
document.getElementById("next").onclick=()=>move(1);
slider.oninput=e=>{{position=Number(e.target.value);render();}};
select.onchange=e=>{{position=Number(e.target.value);render();}};
document.querySelectorAll(".map-btn").forEach(b=>b.onclick=()=>{{mapKey=b.dataset.map;render();}});
document.querySelectorAll(".filter-btn").forEach(b=>b.onclick=()=>{{filterKey=b.dataset.filter;rebuildVisible(false);}});
window.addEventListener("keydown",e=>{{
  if(e.key==="ArrowLeft") move(-1);
  if(e.key==="ArrowRight") move(1);
  if(e.key.toLowerCase()==="i") {{mapKey="intruder";render();}}
  if(e.key.toLowerCase()==="b") {{mapKey="baseline";render();}}
  if(e.key.toLowerCase()==="a") {{filterKey="all";rebuildVisible(false);}}
}});
rebuildVisible(false);
</script>
</body>
</html>
"""
    path.write_text(html)


def main() -> None:
    refuse_overwrite(HTML_OUT, MANIFEST_OUT, THUMBNAIL_OUT)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    log("Loading Allen Nissl and annotation volumes")
    nissl, annotation = load_atlas()
    log(f"Atlas shape: {nissl.shape}")
    frames, manifest, preview = build_frames(nissl, annotation)
    slim = {
        **manifest,
        "html": str(HTML_OUT),
        "preview": str(THUMBNAIL_OUT),
        "frame_index": [
            {
                key: frame[key]
                for key in (
                    "slice_position",
                    "source_number",
                    "selected_allen_index",
                    "selected_bregma_mm",
                    "dice",
                    "inside_fraction",
                    "accepted_final",
                    "status_key",
                    "status_label",
                    "rejection_reason",
                    "continuity_note",
                    "baseline_available",
                )
            }
            for frame in frames
        ],
    }
    MANIFEST_OUT.write_text(json.dumps(slim, indent=2))
    write_html(frames, manifest, HTML_OUT)
    Image.open(io.BytesIO(preview)).save(THUMBNAIL_OUT)
    log(f"Wrote browser: {HTML_OUT}")
    log(f"Wrote manifest: {MANIFEST_OUT}")
    log(f"Wrote preview: {THUMBNAIL_OUT}")


if __name__ == "__main__":
    main()
