#!/usr/bin/env python
"""Build an interactive Allen AP z-stack registration browser.

Layout per frame:
  - Left: input Template anatomy with selected metric painted (source space)
  - Right: Allen Nissl coronal plane with the same metric warped into CCF
  - Upper-right inset: mid-sagittal Nissl overview with current AP cut

Selectable metrics (WildR−SPF polarity: positive/blue = WildR higher):
  Intruder t, Baseline WildR−SPF, WildR Intruder−Baseline, SPF Intruder−Baseline,
  Template 2DG, registration tissue mask

Usage:
  /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \\
    notebooks/python_scripts/build_2dg_ccf_zstack_browser.py
"""

from __future__ import annotations

import base64
import io
import json
from pathlib import Path

import nrrd
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image

WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
SOURCE_DATE = "2026-11-11"
CONTROL_DATE = "2026-11-12"
ANALYSIS_DATE = "2026-11-14e"

SOURCE_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / SOURCE_DATE
CONTROL_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / CONTROL_DATE
OUT_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
FIGURE_DIR = WKDIR / "figures/2DG_CCF_Integration"
ATLAS_ROOT = WKDIR / "data/input/allen_registration_ref"

NISSL_PATH = ATLAS_ROOT / "ara_nissl_25.nrrd"
ANNOTATION_PATH = ATLAS_ROOT / "annotation_25.nrrd"
REG_SUMMARY = SOURCE_DIR / f"2dg_registration_summary_strict_{SOURCE_DATE}.csv"

METRIC_FILES = {
    "intruder_t": SOURCE_DIR / f"2dg_overlay_points_ccf_scalar-t_{SOURCE_DATE}.csv.gz",
    "baseline_wildr_vs_spf": CONTROL_DIR
    / f"2dg_overlay_points_ccf_baseline_wildr_vs_spf_{CONTROL_DATE}.csv.gz",
    "wildr_intruder_vs_baseline": CONTROL_DIR
    / f"2dg_overlay_points_ccf_wildr_intruder_vs_baseline_{CONTROL_DATE}.csv.gz",
    "spf_intruder_vs_baseline": CONTROL_DIR
    / f"2dg_overlay_points_ccf_spf_intruder_vs_baseline_{CONTROL_DATE}.csv.gz",
}

METRIC_LABELS = {
    "intruder_t": "Intruder WildR−SPF t",
    "baseline_wildr_vs_spf": "Baseline WildR−SPF t",
    "wildr_intruder_vs_baseline": "WildR Intruder−Baseline",
    "spf_intruder_vs_baseline": "SPF Intruder−Baseline",
    "template_2dg": "Template 2DG (anatomy)",
    "tissue_mask": "Tissue mask (registration)",
}

HTML_OUT = FIGURE_DIR / f"2DG_CCF_zstack_browser_{ANALYSIS_DATE}.html"
MANIFEST_OUT = OUT_DIR / f"2DG_CCF_zstack_browser_manifest_{ANALYSIS_DATE}.json"
FRAME_DIR = OUT_DIR / "zstack_frames"

T_CLIP = 3.5
CCF_VOXEL_UM = 25.0


def log(msg: str) -> None:
    print(f"[zstack] {msg}", flush=True)


def refuse_overwrite(*paths: Path) -> None:
    existing = [str(p) for p in paths if p.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite:\n" + "\n".join(existing))


def load_atlas():
    nissl, _ = nrrd.read(str(NISSL_PATH))
    annotation, _ = nrrd.read(str(ANNOTATION_PATH))
    return nissl.astype(np.float32), annotation


def to_uint8_img(rgb: np.ndarray) -> Image.Image:
    rgb = np.clip(rgb, 0, 255).astype(np.uint8)
    return Image.fromarray(rgb)


def encode_png(rgb: np.ndarray, max_width: int = 1100) -> str:
    img = to_uint8_img(rgb)
    if img.width > max_width:
        h = int(img.height * max_width / img.width)
        img = img.resize((max_width, h), Image.Resampling.BILINEAR)
    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def gray_to_rgb(gray: np.ndarray, vmin=None, vmax=None) -> np.ndarray:
    g = gray.astype(np.float32)
    if vmin is None:
        vmin = np.nanpercentile(g[np.isfinite(g) & (g > 0)], 2) if np.any(g > 0) else 0
    if vmax is None:
        vmax = np.nanpercentile(g[np.isfinite(g) & (g > 0)], 98) if np.any(g > 0) else 1
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = 0.0, 1.0
    norm = (g - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0, 1)
    norm = np.nan_to_num(norm, nan=0.0)
    return (np.stack([norm, norm, norm], axis=-1) * 255.0).astype(np.float32)


def diverging_overlay(base_rgb: np.ndarray, values: np.ndarray, mask: np.ndarray, alpha: float = 0.72) -> np.ndarray:
    """Paint signed values with RdBu_r-like mapping onto RGB base."""
    out = base_rgb.astype(np.float32).copy()
    if not np.any(mask):
        return out
    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=-T_CLIP, vmax=T_CLIP)
    rgba = cmap(norm(np.clip(values, -T_CLIP, T_CLIP)))
    for c in range(3):
        channel = out[..., c]
        painted = rgba[..., c] * 255.0
        channel[mask] = (1 - alpha) * channel[mask] + alpha * painted[mask]
        out[..., c] = channel
    return out


def rasterize_source(values: np.ndarray, ys: np.ndarray, xs: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    grid = np.full(shape, np.nan, dtype=np.float32)
    counts = np.zeros(shape, dtype=np.float32)
    yi = np.clip(ys.astype(int), 0, shape[0] - 1)
    xi = np.clip(xs.astype(int), 0, shape[1] - 1)
    ok = np.isfinite(values)
    yi, xi, values = yi[ok], xi[ok], values[ok]
    for y, x, v in zip(yi, xi, values):
        if np.isnan(grid[y, x]):
            grid[y, x] = v
            counts[y, x] = 1
        else:
            counts[y, x] += 1
            grid[y, x] += (v - grid[y, x]) / counts[y, x]
    return grid


def rasterize_atlas(values: np.ndarray, atlas_y: np.ndarray, atlas_x: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    return rasterize_source(values, atlas_y, atlas_x, shape)


def sagittal_overview(nissl: np.ndarray, ap_index: int, width: int = 220) -> np.ndarray:
    # nissl: (AP, DV, ML) -> mid-sagittal ML plane
    mid_ml = nissl.shape[2] // 2
    sag = nissl[:, :, mid_ml].T  # (DV, AP)
    rgb = gray_to_rgb(np.log1p(np.maximum(sag, 0)))
    # draw AP cut
    ap = int(np.clip(ap_index, 0, rgb.shape[1] - 1))
    rgb[:, max(0, ap - 1) : ap + 2, :] = np.array([255, 80, 0], dtype=np.float32)
    img = to_uint8_img(rgb)
    h = int(img.height * width / img.width)
    return np.asarray(img.resize((width, h), Image.Resampling.BILINEAR), dtype=np.float32)


def colorbar_strip(height: int = 280, width: int = 70) -> np.ndarray:
    fig, ax = plt.subplots(figsize=(1.1, 3.4))
    sm = ScalarMappable(norm=Normalize(vmin=-T_CLIP, vmax=T_CLIP), cmap="RdBu_r")
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax)
    cb.set_ticks([-T_CLIP, 0, T_CLIP])
    cb.set_ticklabels([f"SPF↑\n{-T_CLIP:g}", "0", f"WildR↑\n{T_CLIP:g}"])
    cb.set_label("signed effect\n(WildR − SPF)\nL← mid →R", fontsize=7)
    fig.tight_layout(pad=0.2)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    img = Image.open(buf).convert("RGB")
    img = img.resize((width, height), Image.Resampling.BILINEAR)
    return np.asarray(img, dtype=np.float32)


def compose_frame(
    input_rgb: np.ndarray,
    atlas_rgb: np.ndarray,
    overview_rgb: np.ndarray,
    cbar_rgb: np.ndarray,
    title: str,
    subtitle: str,
) -> np.ndarray:
    # Normalize panel heights
    target_h = 420
    def fit_h(arr, h=target_h):
        img = to_uint8_img(arr)
        w = int(img.width * h / img.height)
        return np.asarray(img.resize((w, h), Image.Resampling.BILINEAR), dtype=np.float32)

    left = fit_h(input_rgb)
    right = fit_h(atlas_rgb)
    gap = np.full((target_h, 16, 3), 255.0, dtype=np.float32)
    mid = np.concatenate([left, gap, right], axis=1)

    # Canvas with header + overview inset
    pad_top = 78
    pad_right = overview_rgb.shape[1] + cbar_rgb.shape[1] + 36
    canvas = np.full((mid.shape[0] + pad_top, mid.shape[1] + pad_right, 3), 255.0, dtype=np.float32)
    canvas[pad_top:, : mid.shape[1], :] = mid

    # Place overview and colorbar on the right
    ov = overview_rgb
    cb = cbar_rgb
    x0 = mid.shape[1] + 12
    canvas[pad_top : pad_top + ov.shape[0], x0 : x0 + ov.shape[1], :] = ov
    cb_y = pad_top + ov.shape[0] + 12
    cb_h = min(cb.shape[0], canvas.shape[0] - cb_y - 8)
    if cb_h > 20:
        cb_img = to_uint8_img(cb).resize((cb.shape[1], cb_h), Image.Resampling.BILINEAR)
        canvas[cb_y : cb_y + cb_h, x0 : x0 + cb.shape[1], :] = np.asarray(cb_img, dtype=np.float32)

    # Titles via matplotlib text render onto strip
    fig = plt.figure(figsize=(canvas.shape[1] / 100, pad_top / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.01, 0.70, title, fontsize=13, fontweight="bold", va="center")
    ax.text(0.01, 0.28, subtitle, fontsize=9, va="center", color="#333333")
    ax.text(0.72, 0.70, "Input slice", fontsize=10, ha="center", fontweight="bold")
    ax.text(0.90, 0.70, "Allen CCF plane", fontsize=10, ha="center", fontweight="bold")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor="white")
    plt.close(fig)
    header = np.asarray(Image.open(buf).convert("RGB").resize((canvas.shape[1], pad_top)), dtype=np.float32)
    canvas[:pad_top, :, :] = header

    # Panel labels on image
    return canvas


def load_metric_points(metric_key: str) -> pd.DataFrame | None:
    path = METRIC_FILES.get(metric_key)
    if path is None or not path.exists():
        return None
    cols = ["source_number", "ccfy", "ccfz", "t_value"]
    df = pd.read_csv(path, usecols=lambda c: c in cols + ["t_value_status"])
    return df


def build_frames(nissl: np.ndarray, annotation: np.ndarray) -> tuple[list[dict], dict]:
    reg = pd.read_csv(REG_SUMMARY)
    accepted = reg[reg["accepted_final"]].sort_values("slice_position").reset_index(drop=True)
    log(f"Accepted slices: {len(accepted)}")

    metric_points = {k: load_metric_points(k) for k in METRIC_FILES}
    cbar = colorbar_strip()

    metric_keys = [
        "intruder_t",
        "baseline_wildr_vs_spf",
        "wildr_intruder_vs_baseline",
        "spf_intruder_vs_baseline",
        "template_2dg",
        "tissue_mask",
    ]

    frames_meta = []
    FRAME_DIR.mkdir(parents=True, exist_ok=True)

    for i, row in accepted.iterrows():
        slice_pos = int(row["slice_position"])
        source_number = int(row["source_number"])
        ap_index = int(row["selected_allen_index"])
        dice = float(row["dice"])
        bregma = float(row["selected_bregma_mm"])
        prepared_path = SOURCE_DIR / "prepared" / f"slice_{slice_pos:03d}.npz"
        pixels_path = SOURCE_DIR / "work" / f"slice_{slice_pos:03d}" / "pixels_ccf.csv.gz"
        if not prepared_path.exists() or not pixels_path.exists():
            log(f"Skipping missing slice {slice_pos}")
            continue

        prepared = np.load(prepared_path)
        template = prepared["template_y"]
        tissue = prepared["tissue"].astype(bool)
        template_2dg = prepared["template_2dg"]
        t_full = prepared["t_value"]
        pixels = pd.read_csv(
            pixels_path,
            usecols=[
                "source_x",
                "source_y",
                "atlas_x",
                "atlas_y",
                "t_value",
                "template_2dg",
                "overlay_present",
                "inside_allen_tissue",
            ],
        )

        allen_plane = nissl[ap_index]
        allen_rgb = gray_to_rgb(np.log1p(np.maximum(allen_plane, 0)))
        overview = sagittal_overview(nissl, ap_index)
        base_input = gray_to_rgb(template, vmin=0, vmax=255)

        # Dense Intruder / Template fields from prepared + pixels for atlas
        atlas_shape = allen_plane.shape
        atlas_intruder = rasterize_atlas(
            pixels["t_value"].to_numpy(),
            pixels["atlas_y"].to_numpy(),
            pixels["atlas_x"].to_numpy(),
            atlas_shape,
        )

        per_metric_images = {}
        for metric in metric_keys:
            if metric == "template_2dg":
                input_rgb = gray_to_rgb(template_2dg)
                atlas_vals = rasterize_atlas(
                    pixels["template_2dg"].to_numpy(),
                    pixels["atlas_y"].to_numpy(),
                    pixels["atlas_x"].to_numpy(),
                    atlas_shape,
                )
                atlas_rgb = gray_to_rgb(atlas_vals)
            elif metric == "tissue_mask":
                input_rgb = base_input.copy()
                input_rgb[~tissue] = input_rgb[~tissue] * 0.25
                atlas_rgb = allen_rgb.copy()
                inside = pixels["inside_allen_tissue"].astype(bool).to_numpy()
                if inside.any():
                    ay = np.clip(pixels.loc[inside, "atlas_y"].to_numpy().astype(int), 0, atlas_shape[0] - 1)
                    ax = np.clip(pixels.loc[inside, "atlas_x"].to_numpy().astype(int), 0, atlas_shape[1] - 1)
                    atlas_rgb[ay, ax, 0] = 255
                    atlas_rgb[ay, ax, 1] = 140
                    atlas_rgb[ay, ax, 2] = 0
            elif metric == "intruder_t":
                mask = tissue & np.isfinite(t_full) & (np.abs(t_full) >= 1.725)
                input_rgb = diverging_overlay(base_input, t_full, mask)
                amask = np.isfinite(atlas_intruder) & (np.abs(atlas_intruder) >= 1.725)
                atlas_rgb = diverging_overlay(allen_rgb, atlas_intruder, amask)
            else:
                pts = metric_points.get(metric)
                input_rgb = base_input.copy()
                atlas_rgb = allen_rgb.copy()
                if pts is not None:
                    sub = pts[pts["source_number"] == source_number]
                    if len(sub):
                        # Approximate source paint via atlas→source inverse is unavailable;
                        # show atlas paint from CCF coords and faint source silhouette.
                        ay = (sub["ccfy"].to_numpy() / CCF_VOXEL_UM).astype(int)
                        ax = (sub["ccfz"].to_numpy() / CCF_VOXEL_UM).astype(int)
                        vals = sub["t_value"].to_numpy()
                        grid = np.full(atlas_shape, np.nan, dtype=np.float32)
                        ok = (
                            (ay >= 0)
                            & (ay < atlas_shape[0])
                            & (ax >= 0)
                            & (ax < atlas_shape[1])
                            & np.isfinite(vals)
                        )
                        for y, x, v in zip(ay[ok], ax[ok], vals[ok]):
                            grid[y, x] = v if np.isnan(grid[y, x]) else 0.5 * (grid[y, x] + v)
                        amask = np.isfinite(grid)
                        atlas_rgb = diverging_overlay(allen_rgb, np.nan_to_num(grid), amask, alpha=0.85)
                        # Sibling maps lack dense source-space warps: left = Template only.
                        input_rgb = base_input.copy()

            title = (
                f"Slice {source_number}  →  Allen AP index {ap_index}  "
                f"(bregma {bregma:+.2f} mm)  |  Dice {dice:.3f}"
            )
            subtitle = (
                f"Metric: {METRIC_LABELS[metric]}   ·   "
                "blue = WildR↑, red = SPF↑   ·   "
                "Inset: mid-sagittal Nissl + AP cut   ·   "
                "ML: left ← midline(5700µm) → right"
            )
            if metric not in {"intruder_t", "template_2dg", "tissue_mask"}:
                subtitle += "   ·   left=Template anatomy only; metric painted on Allen panel"

            frame_rgb = compose_frame(input_rgb, atlas_rgb, overview, cbar, title, subtitle)
            out_png = FRAME_DIR / f"slice_{slice_pos:03d}__{metric}.png"
            to_uint8_img(frame_rgb).save(out_png, optimize=True)
            per_metric_images[metric] = encode_png(frame_rgb)

        frames_meta.append(
            {
                "index": int(i),
                "slice_position": slice_pos,
                "source_number": source_number,
                "selected_allen_index": ap_index,
                "selected_bregma_mm": bregma,
                "dice": dice,
                "images": per_metric_images,
            }
        )
        if (i + 1) % 10 == 0:
            log(f"Rendered {i + 1}/{len(accepted)} slices")

    manifest = {
        "analysis_date": ANALYSIS_DATE,
        "n_slices": len(frames_meta),
        "metrics": metric_keys,
        "metric_labels": METRIC_LABELS,
        "polarity": {
            "signed_maps": "WildR minus SPF",
            "positive_or_blue": "WildR higher",
            "negative_or_red": "SPF higher",
        },
    }
    return frames_meta, manifest


def write_html(frames: list[dict], manifest: dict, path: Path) -> None:
    # Compact JSON payload (base64 images already in frames)
    payload = {"manifest": manifest, "frames": frames}
    payload_json = json.dumps(payload)
    metric_buttons = "\n".join(
        f'<button type="button" data-metric="{k}" class="metric-btn">{METRIC_LABELS[k]}</button>'
        for k in manifest["metrics"]
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>2DG × Allen CCF z-stack browser</title>
<style>
  :root {{
    --bg: #0f1115;
    --panel: #1a1f27;
    --text: #f2f2f2;
    --muted: #a8b0bd;
    --accent: #ff6c0c;
    --line: #2c3440;
  }}
  html, body {{ margin:0; padding:0; background:var(--bg); color:var(--text);
    font-family: Helvetica, Arial, sans-serif; }}
  .wrap {{ max-width: 1500px; margin: 0 auto; padding: 16px 18px 28px; }}
  h1 {{ font-size: 20px; margin: 0 0 6px; }}
  .sub {{ color: var(--muted); font-size: 13px; margin-bottom: 14px; line-height: 1.4; }}
  .bar {{ display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-bottom:10px; }}
  button {{ background:#243041; color:var(--text); border:1px solid var(--line);
    border-radius:8px; padding:8px 12px; cursor:pointer; font-size:12px; }}
  button.active, button:hover {{ background:var(--accent); color:#111; border-color:var(--accent); font-weight:700; }}
  input[type=range] {{ width: min(720px, 70vw); }}
  .meta {{ color:var(--muted); font-size:12px; min-width: 280px; }}
  .stage {{ background: var(--panel); border:1px solid var(--line); border-radius:12px;
    padding:10px; overflow:auto; }}
  .stage img {{ width:100%; height:auto; display:block; border-radius:8px; background:#fff; }}
  .hint {{ margin-top:10px; color:var(--muted); font-size:12px; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>Allen CCF z-stack registration browser</h1>
  <div class="sub">
    Scroll accepted Template→Allen slices. Left = input section; right = Allen coronal plane with warped metric;
    upper-right inset marks the AP cut on a mid-sagittal Nissl overview.
    Signed maps: <b>blue = WildR higher</b>, <b>red = SPF higher</b>.
  </div>
  <div class="bar">
    <button type="button" id="prev">◀ Prev</button>
    <input id="slider" type="range" min="0" max="{max(len(frames)-1, 0)}" value="0"/>
    <button type="button" id="next">Next ▶</button>
    <span class="meta" id="meta"></span>
  </div>
  <div class="bar" id="metrics">
    {metric_buttons}
  </div>
  <div class="stage"><img id="view" alt="z-stack frame"/></div>
  <div class="hint">Keyboard: ←/→ change slice · 1–6 change metric. {manifest['n_slices']} accepted slices.</div>
</div>
<script>
const PAYLOAD = {payload_json};
let idx = 0;
let metric = PAYLOAD.manifest.metrics[0];
const view = document.getElementById('view');
const slider = document.getElementById('slider');
const meta = document.getElementById('meta');

function setMetricButtons() {{
  document.querySelectorAll('.metric-btn').forEach(btn => {{
    btn.classList.toggle('active', btn.dataset.metric === metric);
  }});
}}

function render() {{
  const fr = PAYLOAD.frames[idx];
  if (!fr) return;
  view.src = fr.images[metric];
  meta.textContent = `slice ${{fr.source_number}} · Allen AP ${{fr.selected_allen_index}} · bregma ${{fr.selected_bregma_mm.toFixed(2)}} mm · Dice ${{fr.dice.toFixed(3)}} · ${{idx+1}}/${{PAYLOAD.frames.length}}`;
  slider.value = String(idx);
  setMetricButtons();
}}

document.getElementById('prev').onclick = () => {{ idx = Math.max(0, idx-1); render(); }};
document.getElementById('next').onclick = () => {{ idx = Math.min(PAYLOAD.frames.length-1, idx+1); render(); }};
slider.oninput = (e) => {{ idx = parseInt(e.target.value, 10); render(); }};
document.querySelectorAll('.metric-btn').forEach(btn => {{
  btn.onclick = () => {{ metric = btn.dataset.metric; render(); }};
}});
window.addEventListener('keydown', (e) => {{
  if (e.key === 'ArrowLeft') {{ idx = Math.max(0, idx-1); render(); }}
  if (e.key === 'ArrowRight') {{ idx = Math.min(PAYLOAD.frames.length-1, idx+1); render(); }}
  const num = parseInt(e.key, 10);
  if (num >= 1 && num <= PAYLOAD.manifest.metrics.length) {{
    metric = PAYLOAD.manifest.metrics[num-1];
    render();
  }}
}});
render();
</script>
</body>
</html>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html)


def main() -> None:
    refuse_overwrite(HTML_OUT, MANIFEST_OUT)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    log("Loading Allen atlas")
    nissl, annotation = load_atlas()
    log(f"Atlas shape {nissl.shape}")
    frames, manifest = build_frames(nissl, annotation)
    if not frames:
        raise RuntimeError("No frames rendered")
    # Manifest without huge base64
    slim = {
        **manifest,
        "frame_index": [
            {
                "slice_position": f["slice_position"],
                "source_number": f["source_number"],
                "selected_allen_index": f["selected_allen_index"],
                "dice": f["dice"],
            }
            for f in frames
        ],
        "html": str(HTML_OUT),
        "frame_dir": str(FRAME_DIR),
    }
    if MANIFEST_OUT.exists():
        raise FileExistsError(MANIFEST_OUT)
    MANIFEST_OUT.write_text(json.dumps(slim, indent=2))
    log("Writing HTML browser")
    write_html(frames, manifest, HTML_OUT)
    log(f"Done: {HTML_OUT}")


if __name__ == "__main__":
    main()
