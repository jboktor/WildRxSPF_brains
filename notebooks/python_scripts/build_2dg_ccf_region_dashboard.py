#!/usr/bin/env python3
"""Build the linked 2DG Allen-domain ranking and anatomy dashboard.

The dashboard is descriptive. It ranks the curated Allen CCF domains and
subdomains by the voxel-weighted mean signed WildR-minus-SPF t value among
significant voxels from accepted registrations. It does not manufacture
parcel-level p values, FDR values, confidence intervals, or animal-level
replication from group statistical maps.

Outputs are date-stamped and never overwrite existing artifacts:

  tables/2DG_CCF_Integration/
    2DG_CCF_domain_subdomain_scores_2026-11-22e.csv
  data/interim/2DG_CCF_Integration/2026-11-22e/
    2DG_CCF_region_dashboard_manifest_2026-11-22e.json
  figures/2DG_CCF_Integration/
    2DG_CCF_region_dashboard_2026-11-22e.html

Run with:

  /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
    notebooks/python_scripts/build_2dg_ccf_region_dashboard.py
"""

from __future__ import annotations

import base64
import io
import json
import pickle
import sys
from pathlib import Path

import nrrd
import numpy as np
import pandas as pd
from PIL import Image
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_2dg_ccf_region_meshes import TAXONOMY, safe_key


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
SOURCE_DATE = "2026-11-23c"
CONTROL_DATE = "2026-11-23c"
MESH_DATE = "2026-11-21e"
ANALYSIS_DATE = "2026-11-23k"

SOURCE_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / SOURCE_DATE
CONTROL_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / CONTROL_DATE
MESH_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / MESH_DATE
OUT_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
TABLE_DIR = WKDIR / "tables/2DG_CCF_Integration"
FIGURE_DIR = WKDIR / "figures/2DG_CCF_Integration"
ATLAS_DIR = WKDIR / "data/input/allen_registration_ref"

REGISTRATION_PATH = SOURCE_DIR / f"2dg_registration_summary_strict_{SOURCE_DATE}.csv"
SLICE_MANIFEST_PATH = SOURCE_DIR / f"2dg_slice_manifest_{SOURCE_DATE}.csv"
INTRUDER_PATH = SOURCE_DIR / f"2dg_overlay_points_ccf_scalar-t_{SOURCE_DATE}.csv.gz"
BASELINE_PATH = (
    CONTROL_DIR
    / f"2dg_overlay_points_ccf_baseline_wildr_vs_spf_{CONTROL_DATE}.csv.gz"
)
ANNOTATION_PATH = ATLAS_DIR / "annotation_25.nrrd"
NISSL_PATH = ATLAS_DIR / "ara_nissl_25.nrrd"
ALLEN_LUT_PATH = ATLAS_DIR / "allen_lut.csv"
HIERARCHY_PATH = ATLAS_DIR / "allen_name_to_annots.pkl"
BASELINE_ROOT = WKDIR / "data/input/2DG_data/WilderR vs SPF _Baseline"
MESH_MANIFEST = (
    MESH_DIR / f"allen_ccf_region_mesh_manifest_{MESH_DATE}.json"
)
MESH_VERTICES = (
    MESH_DIR / f"allen_ccf_region_mesh_vertices_{MESH_DATE}.csv.gz"
)
MESH_FACES = MESH_DIR / f"allen_ccf_region_mesh_faces_{MESH_DATE}.csv.gz"
PLOTLY_JS = Path(
    "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/"
    "spatialomics/lib/R/library/plotly/htmlwidgets/lib/plotlyjs/"
    "plotly-latest.min.js"
)

SCORES_OUT = (
    TABLE_DIR / f"2DG_CCF_domain_subdomain_scores_{ANALYSIS_DATE}.csv"
)
MANIFEST_OUT = (
    OUT_DIR / f"2DG_CCF_region_dashboard_manifest_{ANALYSIS_DATE}.json"
)
HTML_OUT = FIGURE_DIR / f"2DG_CCF_region_dashboard_{ANALYSIS_DATE}.html"

VOXEL_UM = 25.0
T_THRESHOLD = 1.725
T_CLIP = 4.0
ROOT_FACE_CAP = 6_000
DOMAIN_FACE_CAP = 3_500
SUBDOMAIN_FACE_CAP = 1_800


def log(message: str) -> None:
    print(f"[2DG-region-dashboard] {message}", flush=True)


def refuse_overwrite(*paths: Path) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing artifacts:\n" + "\n".join(existing)
        )


def clean_json(value):
    if isinstance(value, dict):
        return {str(key): clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else round(float(value), 4)
    return value


def image_data_url(array: np.ndarray, fmt: str, **save_kwargs) -> str:
    buffer = io.BytesIO()
    Image.fromarray(array).save(buffer, format=fmt, **save_kwargs)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    mime = "image/webp" if fmt.upper() == "WEBP" else "image/png"
    return f"data:{mime};base64,{encoded}"


def grayscale(array: np.ndarray) -> np.ndarray:
    values = np.log1p(np.maximum(np.asarray(array, dtype=np.float32), 0))
    finite = np.isfinite(values) & (values > 0)
    if not finite.any():
        return np.zeros(values.shape, dtype=np.uint8)
    lo, hi = np.percentile(values[finite], [1.5, 98.5])
    scaled = np.clip((values - lo) / max(hi - lo, 1e-8), 0, 1)
    return np.round(scaled * 255).astype(np.uint8)


def build_regions(hierarchy: dict[str, list[int]]) -> list[dict]:
    regions: list[dict] = []
    for domain_order, (domain_label, spec) in enumerate(TAXONOMY.items()):
        subdomains = []
        for sub_order, (acronym, (source_name, color)) in enumerate(
            spec["subdomains"].items()
        ):
            if source_name not in hierarchy:
                raise KeyError(f"Missing Allen hierarchy entry: {source_name}")
            subdomains.append(
                {
                    "key": f"region__{acronym}",
                    "label": f"{acronym} · {source_name}",
                    "short_label": acronym,
                    "group": domain_label,
                    "level": "subdomain",
                    "color": color,
                    "annotation_ids": sorted(
                        set(map(int, hierarchy[source_name]))
                    ),
                    "domain_order": domain_order,
                    "subdomain_order": sub_order,
                }
            )
        source_name = spec["source"]
        if source_name is None:
            domain_ids = sorted(
                set().union(
                    *(set(item["annotation_ids"]) for item in subdomains)
                )
            )
        else:
            if source_name not in hierarchy:
                raise KeyError(f"Missing Allen hierarchy entry: {source_name}")
            domain_ids = sorted(set(map(int, hierarchy[source_name])))
        regions.append(
            {
                "key": f"domain__{safe_key(domain_label)}",
                "label": domain_label,
                "short_label": domain_label,
                "group": domain_label,
                "level": "domain",
                "color": spec["color"],
                "annotation_ids": domain_ids,
                "domain_order": domain_order,
                "subdomain_order": -1,
            }
        )
        regions.extend(subdomains)
    if len(regions) != 89:
        raise ValueError(f"Expected 89 domain/subdomain regions, found {len(regions)}")
    return regions


def load_overlay(path: Path) -> pd.DataFrame:
    columns = [
        "slice_position",
        "source_number",
        "ccfx",
        "ccfy",
        "ccfz",
        "annotation",
        "t_value",
    ]
    frame = pd.read_csv(path, usecols=columns)
    finite = np.isfinite(
        frame[["ccfx", "ccfy", "ccfz", "t_value"]].to_numpy(dtype=float)
    ).all(axis=1)
    frame = frame.loc[finite].copy()
    frame["annotation"] = frame["annotation"].astype(int)
    if frame.empty:
        raise ValueError(f"No finite overlay points: {path}")
    if (frame["t_value"].abs() < T_THRESHOLD - 1e-6).any():
        raise ValueError(f"Overlay includes sub-threshold values: {path}")
    return frame


def aggregate_scores(
    overlays: dict[str, pd.DataFrame], regions: list[dict]
) -> pd.DataFrame:
    rows: list[dict] = []
    for contrast, frame in overlays.items():
        for region in regions:
            subset = frame[frame["annotation"].isin(region["annotation_ids"])]
            score = float(subset["t_value"].mean()) if len(subset) else np.nan
            rows.append(
                {
                    "contrast": contrast,
                    "key": region["key"],
                    "label": region["label"],
                    "short_label": region["short_label"],
                    "group": region["group"],
                    "level": region["level"],
                    "color": region["color"],
                    "domain_order": region["domain_order"],
                    "subdomain_order": region["subdomain_order"],
                    "mean_t": score,
                    "mean_abs_t": (
                        float(subset["t_value"].abs().mean())
                        if len(subset)
                        else np.nan
                    ),
                    "n_voxels": int(len(subset)),
                    "n_slices": int(subset["source_number"].nunique()),
                    "frac_wildr_positive": (
                        float((subset["t_value"] > 0).mean())
                        if len(subset)
                        else np.nan
                    ),
                }
            )
    scores = pd.DataFrame(rows)
    if scores.duplicated(["contrast", "key"]).any():
        raise ValueError("Duplicate score rows")
    return scores


def compact_points(frame: pd.DataFrame) -> dict[str, list]:
    return {
        "s": frame["source_number"].astype(int).tolist(),
        "a": frame["annotation"].astype(int).tolist(),
        "x": frame["ccfx"].round(1).tolist(),
        "y": frame["ccfz"].round(1).tolist(),
        "z": frame["ccfy"].round(1).tolist(),
        "t": frame["t_value"].clip(-T_CLIP, T_CLIP).round(3).tolist(),
    }


def compact_meshes() -> dict[str, dict]:
    manifest = json.loads(MESH_MANIFEST.read_text())
    vertices = pd.read_csv(MESH_VERTICES)
    faces = pd.read_csv(MESH_FACES)
    info = {item["key"]: item for item in manifest["regions"]}
    output: dict[str, dict] = {}
    for key, region_info in info.items():
        region_vertices = vertices.loc[
            vertices["region"] == key, ["vertex_id", "ap", "dv", "ml"]
        ]
        region_faces = faces.loc[faces["region"] == key, ["i", "j", "k"]]
        if key == "root":
            cap = ROOT_FACE_CAP
        elif region_info["level"] == "domain":
            cap = DOMAIN_FACE_CAP
        else:
            cap = SUBDOMAIN_FACE_CAP
        if len(region_faces) > cap:
            take = np.linspace(0, len(region_faces) - 1, cap).round().astype(int)
            region_faces = region_faces.iloc[np.unique(take)]
        used = np.unique(region_faces[["i", "j", "k"]].to_numpy().ravel())
        selected = region_vertices.set_index("vertex_id").loc[used]
        remap = {int(old): new for new, old in enumerate(used)}
        remapped = region_faces[["i", "j", "k"]].replace(remap).astype(int)
        output[key] = {
            "label": region_info["label"],
            "group": region_info["group"],
            "level": region_info["level"],
            "color": region_info["color"],
            "x": selected["ap"].round(1).tolist(),
            "y": selected["ml"].round(1).tolist(),
            "z": selected["dv"].round(1).tolist(),
            "i": remapped["i"].tolist(),
            "j": remapped["j"].tolist(),
            "k": remapped["k"].tolist(),
        }
    expected = {"root"} | {item["key"] for item in info.values() if item["key"] != "root"}
    if set(output) != expected:
        raise ValueError("Compact mesh set differs from mesh manifest")
    return output


def status_for_row(row: pd.Series) -> tuple[str, str]:
    if bool(row["accepted_final"]):
        return "accepted", "Accepted"
    reason = str(row.get("rejection_reason", "") or "")
    if reason == "blank_or_empty_template":
        return "blank", "Blank input"
    if bool(row.get("accepted_slice_fit", False)):
        return "continuity", "Continuity reject"
    return "fit_failed", "Fit rejected"


def atlas_frames(
    registration: pd.DataFrame,
    slice_manifest: pd.DataFrame,
    annotation: np.ndarray,
    nissl: np.ndarray,
    lut: np.ndarray,
    regions: list[dict],
    overlays: dict[str, pd.DataFrame],
) -> list[dict]:
    # Base anatomy matches registration QC: grayscale Nissl + dark annotation
    # borders only. No Allen palette fill and no selected-region color wash —
    # signed-t dots carry the contrast on a neutral grey brain.
    del lut  # retained in signature for call-site compatibility
    region_ids = {
        region["key"]: set(region["annotation_ids"]) for region in regions
    }
    source_inputs = slice_manifest.set_index("source_number")
    frames: list[dict] = []
    for _, row in registration.sort_values("slice_position").iterrows():
        ap_index = int(row["selected_allen_index"])
        labels = annotation[ap_index].astype(np.int32)
        gray = grayscale(nissl[ap_index])
        brain = labels > 0
        blended = np.stack([gray, gray, gray], axis=-1).astype(np.float32)
        blended[~brain] = 255
        # Soften outside-tissue a touch so the brain reads clearly.
        blended[brain] = 0.92 * blended[brain] + 0.08 * 255.0
        boundary = find_boundaries(labels, mode="inner")
        tissue_edge = find_boundaries(brain, mode="outer")
        blended[boundary] = np.array([42, 46, 52], dtype=np.float32)
        blended[tissue_edge] = np.array([24, 26, 30], dtype=np.float32)
        base = image_data_url(
            np.clip(blended, 0, 255).astype(np.uint8),
            "WEBP",
            quality=86,
            method=6,
        )
        encoded_labels = np.zeros((*labels.shape, 3), dtype=np.uint8)
        encoded_labels[..., 0] = labels & 255
        encoded_labels[..., 1] = (labels >> 8) & 255
        encoded_labels[..., 2] = (labels >> 16) & 255
        label_map = image_data_url(
            encoded_labels,
            "PNG",
            optimize=True,
            compress_level=9,
        )
        status_key, status_label = status_for_row(row)
        source_number = int(row["source_number"])
        source_input = source_inputs.loc[source_number]
        intruder_available = (
            Path(str(source_input["intruder_path"])).exists()
            and not bool(source_input["blank_template"])
        )
        baseline_available = (
            BASELINE_ROOT / f"template42{source_number:03d}.tif"
        ).exists()
        present_ids = set(map(int, np.unique(labels)))
        intersections = [
            key
            for key, ids in region_ids.items()
            if ids.intersection(present_ids)
        ]
        frames.append(
            {
                "slice_position": int(row["slice_position"]),
                "source_number": source_number,
                "ap_index": ap_index,
                "bregma_mm": float(row["selected_bregma_mm"]),
                "status_key": status_key,
                "status_label": status_label,
                "accepted_final": bool(row["accepted_final"]),
                "intruder_available": intruder_available,
                "baseline_available": baseline_available,
                "intersections": intersections,
                "base": base,
                "labels": label_map,
            }
        )
    if len(frames) != 66:
        raise ValueError(f"Expected 66 delivered sections, found {len(frames)}")
    return frames


def build_manifest(
    regions: list[dict],
    scores: pd.DataFrame,
    registration: pd.DataFrame,
    frames: list[dict],
) -> dict:
    return {
        "analysis_date": ANALYSIS_DATE,
        "score_definition": (
            "Voxel-weighted mean signed t among significant voxels from "
            "accepted registrations; positive = WildR higher, negative = SPF higher"
        ),
        "coronal2d_style": (
            "grayscale Allen Nissl + annotation borders only; no region fill or "
            "selected-domain tint — signed-t voxels carry the color"
        ),
        "inferential_warning": (
            "Descriptive group-map summary only; no parcel-level p, FDR, SE, "
            "confidence interval, or animal-level model is available"
        ),
        "t_threshold": T_THRESHOLD,
        "n_regions": len(regions),
        "n_domains": sum(item["level"] == "domain" for item in regions),
        "n_subdomains": sum(item["level"] == "subdomain" for item in regions),
        "n_delivered_slices": len(registration),
        "n_accepted_slices": int(registration["accepted_final"].sum()),
        "baseline_unavailable_sources": [
            frame["source_number"]
            for frame in frames
            if not frame["baseline_available"]
        ],
        "domain_semantics": (
            "Allen hierarchy composites; domains can overlap and are not additive"
        ),
        "inputs": {
            "intruder_overlay": str(INTRUDER_PATH),
            "baseline_overlay": str(BASELINE_PATH),
            "registration": str(REGISTRATION_PATH),
            "slice_manifest": str(SLICE_MANIFEST_PATH),
            "annotation": str(ANNOTATION_PATH),
            "nissl": str(NISSL_PATH),
            "hierarchy": str(HIERARCHY_PATH),
            "mesh_manifest": str(MESH_MANIFEST),
        },
        "score_rows": int(len(scores)),
        "outputs": {
            "scores": str(SCORES_OUT),
            "html": str(HTML_OUT),
        },
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>2DG Allen regional signature explorer</title>
<script>__PLOTLY_JS__</script>
<style>
:root{--ink:#181a1e;--muted:#626a75;--line:#dfe3e8;--soft:#f5f6f8;--accent:#ff6c0c;--wildr:#2166ac;--spf:#b2182b}
*{box-sizing:border-box}html,body{margin:0;height:100%;overflow:hidden;background:#fff;color:var(--ink);font-family:Inter,Helvetica,Arial,sans-serif}
button,select,input{font:inherit}.top{height:58px;display:flex;align-items:center;justify-content:space-between;padding:8px 14px;border-bottom:1px solid var(--line);background:#fff}
.top h1{font-size:18px;margin:0;letter-spacing:-.02em}.top p{font-size:11px;color:var(--muted);margin:2px 0 0}.badge{font-size:10px;padding:5px 8px;border:1px solid #f0c7b1;background:#fff6f1;border-radius:99px;color:#8c3b14}
.controls{height:72px;display:flex;gap:12px;align-items:end;padding:8px 14px;border-bottom:1px solid var(--line);background:rgba(255,255,255,.98);overflow-x:auto}
label{display:block;font-size:10px;color:var(--muted);white-space:nowrap}select,input[type=search]{display:block;margin-top:3px;border:1px solid #cfd5dc;border-radius:7px;background:#fff;padding:6px 8px;color:var(--ink);min-width:130px}
.layout{height:calc(100% - 130px);display:grid;grid-template-columns:minmax(420px,46%) 1fr;overflow:hidden}
.rank-panel{border-right:1px solid var(--line);display:flex;flex-direction:column;min-width:0;min-height:0}.rank-head{padding:9px 12px 5px}.rank-head h2{font-size:14px;margin:0}.rank-head p{font-size:10.5px;color:var(--muted);margin:3px 0}.rank-scroll{flex:1;min-height:0;overflow-y:auto;overflow-x:hidden}
#rankPlot{width:100%}.anatomy-panel{position:relative;display:flex;flex-direction:column;min-width:0;min-height:0;overflow:hidden}
.anatomy-head{height:74px;padding:8px 12px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;gap:10px;background:#fff}
.selection h2{font-size:15px;margin:0}.selection p{font-size:10.5px;color:var(--muted);margin:3px 0}.mode-buttons{display:flex;gap:5px;align-items:start}
.mode-buttons button,.slice-buttons button{border:1px solid #cfd5dc;border-radius:7px;background:#fff;padding:6px 8px;font-size:11px;cursor:pointer}.mode-buttons button.active,.slice-buttons button.active{border-color:var(--accent);background:#fff4ec;color:#8b3500;font-weight:700}
.view{position:relative;flex:1;overflow:hidden}.view.hidden{display:none!important}#brain3d{width:100%;height:100%}
#sliceView{display:grid;grid-template-rows:1fr auto;height:100%;background:#f6f7f8}.canvas-wrap{position:relative;min-height:0;display:grid;place-items:center;overflow:hidden}.canvas-wrap canvas{max-width:100%;max-height:100%;image-rendering:auto;background:#fff;box-shadow:0 3px 16px rgba(25,35,45,.12)}
.slice-ui{padding:8px 12px 10px;background:#fff;border-top:1px solid var(--line)}.slice-row{display:flex;align-items:center;gap:8px}.slice-row input{flex:1;accent-color:var(--accent)}.slice-meta{font-size:11px;color:var(--muted);margin-top:4px}.slice-buttons{display:flex;gap:5px}
.legend{position:absolute;right:12px;bottom:12px;z-index:4;width:220px;padding:8px 10px;border:1px solid var(--line);border-radius:8px;background:rgba(255,255,255,.95);box-shadow:0 3px 14px rgba(25,35,45,.08);font-size:10px}.gradient{height:8px;border-radius:99px;background:linear-gradient(90deg,var(--spf),#f7f7f7 50%,var(--wildr));border:1px solid #d8dce2}.legend-labels{display:flex;justify-content:space-between;color:var(--muted);margin-top:2px}
.warning{position:absolute;left:12px;bottom:12px;z-index:4;max-width:330px;padding:7px 9px;border-radius:7px;background:rgba(255,248,243,.95);border:1px solid #efcdb9;color:#80411e;font-size:10px}
.coronal-mode .warning{bottom:70px}
@media(max-width:900px){html,body{overflow:auto}.top{height:auto}.controls{height:auto;align-items:start;flex-wrap:wrap}.layout{height:auto;display:block}.rank-panel{height:650px;border-right:0;border-bottom:1px solid var(--line)}.anatomy-panel{height:720px}.viewer-note{display:none}}
</style>
</head>
<body>
<header class="top"><div><h1>2DG regional signature explorer</h1><p>Linked Allen domain ranking · accepted-map evidence · delivered-section atlas context</p></div><span class="badge">Descriptive group map · not animal-level inference</span></header>
<section class="controls" aria-label="Dashboard controls">
  <label>Contrast<select id="contrast"><option value="intruder">Intruder WildR−SPF</option><option value="baseline">Baseline WildR−SPF</option></select></label>
  <label>Rank by<select id="rankMode"><option value="wildr">WildR enriched</option><option value="spf">SPF enriched</option><option value="magnitude">Absolute difference</option><option value="atlas">Atlas order</option></select></label>
  <label>Metric<select id="metric"><option value="mean_t">Mean signed t</option><option value="mean_abs_t">Mean |t|</option><option value="n_voxels">Significant voxels</option><option value="n_slices">Accepted slices</option></select></label>
  <label>Subdomains/group<select id="topN"><option value="5">Top 5</option><option value="8" selected>Top 8</option><option value="all">All enabled</option></select></label>
  <label>Minimum voxels<select id="minVoxels"><option value="0">No minimum</option><option value="10">10</option><option value="25" selected>25</option><option value="100">100</option></select></label>
  <label>Find region<input id="search" type="search" placeholder="PVH, amygdala…"></label>
</section>
<main class="layout">
  <section class="rank-panel">
    <div class="rank-head"><h2>Domain and subdomain ranking</h2><p>Domain rows summarize their Allen descendant IDs; indented rows are the currently enabled subdomains. Click any bar to inspect anatomy.</p></div>
    <div class="rank-scroll"><div id="rankPlot"></div></div>
  </section>
  <section class="anatomy-panel">
    <div class="anatomy-head">
      <div class="selection"><h2 id="selectionTitle">Hypothalamus</h2><p id="selectionMeta"></p></div>
      <div class="mode-buttons" role="group" aria-label="Anatomy view">
        <button data-mode="slices3d" class="active">3D slices</button>
        <button data-mode="average3d">3D average</button>
        <button data-mode="coronal2d">2D coronal</button>
      </div>
    </div>
    <div class="view" id="threeView"><div id="brain3d"></div></div>
    <div class="view hidden" id="sliceView">
      <div class="canvas-wrap"><canvas id="sliceCanvas" width="456" height="320"></canvas></div>
      <div class="slice-ui">
        <div class="slice-row"><div class="slice-buttons"><button id="prevSlice">◀</button><button id="nextSlice">▶</button><button id="jumpIntersect">Next intersecting</button></div><input id="sliceSlider" type="range" min="0" max="65" value="0"></div>
        <div class="slice-meta" id="sliceMeta"></div>
      </div>
    </div>
    <div class="legend"><div><b>Signed t · WildR − SPF</b></div><div class="gradient"></div><div class="legend-labels"><span>SPF↑ · t&lt;0</span><span>0</span><span>WildR↑ · t&gt;0</span></div><div style="margin-top:4px;color:#626a75;line-height:1.35">Same scale for ranking bars, 3D mesh tint, and slice voxels. 2D coronal keeps a grey brain + region outlines so difference colors stand out. Positive t = WildR higher; negative t = SPF higher.</div></div>
    <div class="warning" id="viewWarning">Scores use accepted registrations only. The coronal browser includes all 66 delivered sections for QC context.</div>
  </section>
</main>
<script>
const DATA=__PAYLOAD__;
const $=id=>document.getElementById(id);
const state={contrast:"intruder",rankMode:"wildr",metric:"mean_t",topN:"8",minVoxels:25,search:"",selected:"domain__hypothalamus",mode:"slices3d",slice:0};
const scoresBy={};
DATA.scores.forEach(r=>scoresBy[`${r.contrast}|${r.key}`]=r);
const regionsBy=Object.fromEntries(DATA.regions.map(r=>[r.key,r]));
const imageCache=new Map();

function scoreRow(key){return scoresBy[`${state.contrast}|${key}`]||null}
function numeric(v){return Number.isFinite(Number(v))?Number(v):0}
// Canonical diverging map for this dashboard:
//   t > 0 → WildR↑ → blue (#2166ac)
//   t < 0 → SPF↑   → red  (#b2182b)
// Ranking bars, region mesh tint, 3D voxels, and 2D coronal dots all use this.
function signedColor(v){const x=Math.max(-4,Math.min(4,numeric(v)))/4;if(x<0){const q=1+x;return `rgb(${Math.round(178+(247-178)*q)},${Math.round(24+(247-24)*q)},${Math.round(43+(247-43)*q)})`}return `rgb(${Math.round(247+(33-247)*x)},${Math.round(247+(102-247)*x)},${Math.round(247+(172-247)*x)})`}
function signedLabel(t){const v=numeric(t);if(!Number.isFinite(v)||v===0)return "t ~ 0";return v>0?`t ${v>=0?"+":""}${v.toFixed(2)} · WildR↑`:`t ${v.toFixed(2)} · SPF↑`}
// Plotly's built-in RdBu maps high→red; reverse so high/WildR→blue matches signedColor/legend.
const SIGNED_T_COLORSCALE=[[0,"#b2182b"],[0.5,"#f7f7f7"],[1,"#2166ac"]];
function metricValue(row){if(!row)return 0;return numeric(row[state.metric])}
function sortValue(row){const v=numeric(row&&row.mean_t);if(state.rankMode==="spf")return v;if(state.rankMode==="magnitude")return -Math.abs(v);if(state.rankMode==="atlas")return numeric(row.domain_order)*100+numeric(row.subdomain_order);return -v}

function rankedRows(){
  const min=Number(state.minVoxels);const q=state.search.trim().toLowerCase();let groups=[];
  DATA.domains.forEach(domain=>{
    const drow=scoreRow(domain.key);let subs=DATA.regions.filter(r=>r.level==="subdomain"&&r.group===domain.label).map(r=>scoreRow(r.key)).filter(Boolean);
    subs=subs.filter(r=>numeric(r.n_voxels)>=min);
    if(q)subs=subs.filter(r=>(r.label+" "+r.group).toLowerCase().includes(q));
    subs.sort((a,b)=>sortValue(a)-sortValue(b));
    if(state.topN!=="all")subs=subs.slice(0,Number(state.topN));
    if(!q||domain.label.toLowerCase().includes(q)||subs.length)groups.push({domain,row:drow,subs});
  });
  groups.sort((a,b)=>state.rankMode==="atlas"?numeric(a.row.domain_order)-numeric(b.row.domain_order):sortValue(a.row)-sortValue(b.row));
  return groups.flatMap(g=>[g.row,...g.subs]);
}

function drawRanks(){
  const rows=rankedRows();const metric=state.metric;
  const labels=rows.map(r=>r.level==="domain"?`<b>${r.label}</b>`:`   ${r.short_label} · ${r.label.split(" · ").slice(1).join(" · ")}`);
  const values=rows.map(r=>metricValue(r));
  const colors=rows.map(r=>metric==="mean_t"?signedColor(r.mean_t):(r.level==="domain"?r.color:"#8f99a5"));
  const widths=rows.map(r=>r.level==="domain"?1:0);
  const text=rows.map(r=>`${r.label}<br>${state.contrast==="intruder"?"Intruder":"Baseline"} WildR−SPF<br>mean ${r.mean_t==null?"NA":signedLabel(r.mean_t)}<br>${r.n_voxels} significant voxels · ${r.n_slices} accepted slices`);
  const height=Math.max(460,rows.length*25+80);
  Plotly.react("rankPlot",[{type:"bar",orientation:"h",x:values,y:labels,customdata:rows.map(r=>r.key),text:text,hoverinfo:"text",marker:{color:colors,line:{color:"#ffffff",width:widths}},selectedpoints:rows.map((r,i)=>r.key===state.selected?i:null).filter(x=>x!=null)}],{height,margin:{l:210,r:24,t:12,b:45},paper_bgcolor:"#fff",plot_bgcolor:"#fff",font:{family:"Inter,Helvetica,Arial",size:10,color:"#222"},xaxis:{title:metric==="mean_t"?"Mean signed t (WildR − SPF)":metric.replaceAll("_"," "),zeroline:true,zerolinecolor:"#555",gridcolor:"#eceef1"},yaxis:{autorange:"reversed",automargin:false},showlegend:false,dragmode:false},{displayModeBar:false,responsive:true});
}

function selectedPoints(){
  const ids=new Set(regionsBy[state.selected].annotation_ids);const p=DATA.points[state.contrast];const out={x:[],y:[],z:[],t:[],s:[]};
  for(let n=0;n<p.a.length;n++)if(ids.has(p.a[n])){out.x.push(p.x[n]);out.y.push(p.y[n]);out.z.push(p.z[n]);out.t.push(p.t[n]);out.s.push(p.s[n])}
  return out;
}

function draw3d(){
  const region=regionsBy[state.selected],mesh=DATA.meshes[state.selected],root=DATA.meshes.root,row=scoreRow(state.selected),points=selectedPoints();
  const rootTrace={type:"mesh3d",x:root.x,y:root.y,z:root.z,i:root.i,j:root.j,k:root.k,color:"#c8c8cf",opacity:.07,hoverinfo:"skip",name:"Whole brain",showscale:false};
  const meanT=row&&row.mean_t!=null?Number(row.mean_t):null;
  const regionColor=signedColor(meanT);
  const selectedTrace={type:"mesh3d",x:mesh.x,y:mesh.y,z:mesh.z,i:mesh.i,j:mesh.j,k:mesh.k,color:regionColor,opacity:state.mode==="average3d"?.78:.22,hovertemplate:`${region.label}<br>mean ${meanT==null?"NA":signedLabel(meanT)}<extra></extra>`,name:region.label,showscale:false};
  const traces=[rootTrace,selectedTrace];
  if(state.mode==="slices3d"){
    const unique=[...new Set(points.x)].sort((a,b)=>a-b);const lx=[],ly=[],lz=[];
    unique.forEach(ap=>{lx.push(ap,ap,ap,ap,ap,null);ly.push(0,11400,11400,0,0,null);lz.push(0,0,8000,8000,0,null)});
    traces.push({type:"scatter3d",mode:"lines",x:lx,y:ly,z:lz,line:{color:"#9aa0a8",width:1},opacity:.16,hoverinfo:"skip",showlegend:false});
    traces.push({type:"scatter3d",mode:"markers",x:points.x,y:points.y,z:points.z,marker:{size:2.5,opacity:.78,color:points.t,colorscale:SIGNED_T_COLORSCALE,cmin:-4,cmax:4,showscale:false},text:points.t.map((t,i)=>`Template ${String(points.s[i]).padStart(3,"0")}<br>${signedLabel(t)}`),hoverinfo:"text",showlegend:false});
  }
  Plotly.react("brain3d",traces,{margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"#fff",scene:{xaxis:{title:"CCF AP (µm)",backgroundcolor:"#fbfbfc",gridcolor:"#eceef1"},yaxis:{title:"CCF ML (µm)",backgroundcolor:"#fbfbfc",gridcolor:"#eceef1"},zaxis:{title:"CCF DV (µm)",autorange:"reversed",backgroundcolor:"#fbfbfc",gridcolor:"#eceef1"},aspectmode:"data",camera:{eye:{x:1.55,y:1.45,z:.85},up:{x:0,y:0,z:1}}},showlegend:false},{displaylogo:false,responsive:true,scrollZoom:true,modeBarButtonsToRemove:["toImage","select3d","lasso3d"]});
  $("viewWarning").textContent=state.mode==="slices3d"
    ?"Slice voxels use the same signed-t colors as the ranking bars (blue=WildR↑ / red=SPF↑). The translucent mesh is tinted by the region's mean signed t."
    :"Selected region mesh is shaded by mean signed t (blue=WildR↑ / red=SPF↑), matching the ranking bars; descriptive only, not parcel-level inference.";
}

function loadImage(url){if(imageCache.has(url))return imageCache.get(url);const promise=new Promise((resolve,reject)=>{const im=new Image();im.onload=()=>resolve(im);im.onerror=reject;im.src=url});imageCache.set(url,promise);return promise}
function pointColor(t,alpha){const x=Math.max(-4,Math.min(4,t))/4;if(x<0){const q=1+x;return `rgba(${Math.round(178+(247-178)*q)},${Math.round(24+(247-24)*q)},${Math.round(43+(247-43)*q)},${alpha})`}return `rgba(${Math.round(247+(33-247)*x)},${Math.round(247+(102-247)*x)},${Math.round(247+(172-247)*x)},${alpha})`}

async function drawSlice(){
  const frame=DATA.frames[state.slice],canvas=$("sliceCanvas"),ctx=canvas.getContext("2d"),region=regionsBy[state.selected],ids=new Set(region.annotation_ids);
  const [base,labelImage]=await Promise.all([loadImage(frame.base),loadImage(frame.labels)]);
  canvas.width=base.width;canvas.height=base.height;ctx.clearRect(0,0,canvas.width,canvas.height);ctx.drawImage(base,0,0);
  // Label mask is used only for selection outline + voxel filtering — no region fill.
  const off=document.createElement("canvas");off.width=labelImage.width;off.height=labelImage.height;const ox=off.getContext("2d",{willReadFrequently:true});ox.drawImage(labelImage,0,0);const pixels=ox.getImageData(0,0,off.width,off.height).data;
  const mask=new Uint8Array(off.width*off.height);for(let p=0,n=0;p<pixels.length;p+=4,n++){const id=pixels[p]+pixels[p+1]*256+pixels[p+2]*65536;if(ids.has(id))mask[n]=1}
  const image=ctx.getImageData(0,0,canvas.width,canvas.height);
  // Emphasize selected-region outline only (no colored fill) so signed-t dots stand out.
  for(let y=1;y<canvas.height-1;y++)for(let x=1;x<canvas.width-1;x++){
    const n=y*canvas.width+x;if(!mask[n])continue;
    if(!mask[n-1]||!mask[n+1]||!mask[n-canvas.width]||!mask[n+canvas.width]){
      const p=n*4;image.data[p]=12;image.data[p+1]=14;image.data[p+2]=18;
    }
  }
  ctx.putImageData(image,0,0);
  if(frame.accepted_final){
    const p=DATA.points[state.contrast];
    // Pass 1: soft halo for selected voxels (QC-style readability on grayscale).
    for(let n=0;n<p.s.length;n++){
      if(p.s[n]!==frame.source_number||!ids.has(p.a[n]))continue;
      const x=p.y[n]/25,y=p.z[n]/25;
      ctx.beginPath();ctx.arc(x,y,4.2,0,Math.PI*2);ctx.fillStyle=pointColor(p.t[n],0.22);ctx.fill();
    }
    // Pass 2: crisp cores — selected larger/opaque; other-plane voxels faint context.
    for(let n=0;n<p.s.length;n++){
      if(p.s[n]!==frame.source_number)continue;
      const selected=ids.has(p.a[n]);
      const x=p.y[n]/25,y=p.z[n]/25;
      ctx.beginPath();ctx.arc(x,y,selected?2.8:1.35,0,Math.PI*2);ctx.fillStyle=pointColor(p.t[n],selected?0.95:0.2);ctx.fill();
    }
  }
  $("sliceSlider").value=state.slice;
  const available=state.contrast==="intruder"?frame.intruder_available:frame.baseline_available;
  const intersects=frame.intersections.includes(state.selected);
  $("sliceMeta").textContent=`Template ${String(frame.source_number).padStart(3,"0")} · ${frame.status_label} · Allen AP ${frame.ap_index} · bregma ${Number(frame.bregma_mm).toFixed(2)} mm · ${intersects?"selection intersects this plane":"selection absent from this plane"}${available?"":" · chosen contrast unavailable"}`;
  $("viewWarning").textContent=frame.accepted_final?(available?"Grey Allen anatomy with region outlines only (no regional fill); voxels use blue=WildR↑ / red=SPF↑, same scale as the ranking bars.":"This contrast TIFF is unavailable for this delivered section."):"This delivered section was withheld from scoring and 3D evidence; atlas anatomy is shown for QC context only.";
}

function updateSelection(){
  const region=regionsBy[state.selected],row=scoreRow(state.selected);
  $("selectionTitle").textContent=region.label;
  $("selectionMeta").textContent=row?`${state.contrast==="intruder"?"Intruder":"Baseline"} · mean ${row.mean_t==null?"NA":signedLabel(row.mean_t)} · ${row.n_voxels} significant voxels · ${row.n_slices} accepted slices`:"No score";
  drawRanks();
  if(state.mode==="coronal2d")drawSlice();else draw3d();
}

function setMode(mode){
  state.mode=mode;document.body.classList.toggle("coronal-mode",mode==="coronal2d");document.querySelectorAll(".mode-buttons button").forEach(b=>b.classList.toggle("active",b.dataset.mode===mode));$("threeView").classList.toggle("hidden",mode==="coronal2d");$("sliceView").classList.toggle("hidden",mode!=="coronal2d");
  if(mode==="coronal2d")drawSlice();else draw3d();
}

$("contrast").onchange=e=>{state.contrast=e.target.value;updateSelection()};
$("rankMode").onchange=e=>{state.rankMode=e.target.value;drawRanks()};
$("metric").onchange=e=>{state.metric=e.target.value;drawRanks()};
$("topN").onchange=e=>{state.topN=e.target.value;drawRanks()};
$("minVoxels").onchange=e=>{state.minVoxels=Number(e.target.value);drawRanks()};
$("search").oninput=e=>{state.search=e.target.value;drawRanks()};
document.querySelectorAll(".mode-buttons button").forEach(b=>b.onclick=()=>setMode(b.dataset.mode));
$("sliceSlider").oninput=e=>{state.slice=Number(e.target.value);drawSlice()};
$("prevSlice").onclick=()=>{state.slice=Math.max(0,state.slice-1);drawSlice()};
$("nextSlice").onclick=()=>{state.slice=Math.min(DATA.frames.length-1,state.slice+1);drawSlice()};
$("jumpIntersect").onclick=()=>{for(let d=1;d<=DATA.frames.length;d++){const i=(state.slice+d)%DATA.frames.length;if(DATA.frames[i].intersections.includes(state.selected)){state.slice=i;drawSlice();break}}};
window.addEventListener("keydown",e=>{if(state.mode!=="coronal2d")return;if(e.key==="ArrowLeft")$("prevSlice").click();if(e.key==="ArrowRight")$("nextSlice").click()});
window.addEventListener("resize",()=>{Plotly.Plots.resize("rankPlot");Plotly.Plots.resize("brain3d")});

drawRanks();
document.getElementById("rankPlot").on("plotly_click",event=>{const key=event.points[0].customdata;if(regionsBy[key]){state.selected=key;updateSelection()}});
updateSelection();
</script>
</body>
</html>
"""


def main() -> None:
    refuse_overwrite(SCORES_OUT, MANIFEST_OUT, HTML_OUT)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    log("Loading Allen hierarchy and curated taxonomy")
    with HIERARCHY_PATH.open("rb") as handle:
        hierarchy: dict[str, list[int]] = pickle.load(handle)
    regions = build_regions(hierarchy)

    log("Loading accepted Intruder and Baseline overlay points")
    overlays = {
        "intruder": load_overlay(INTRUDER_PATH),
        "baseline": load_overlay(BASELINE_PATH),
    }
    registration = pd.read_csv(REGISTRATION_PATH)
    slice_manifest = pd.read_csv(SLICE_MANIFEST_PATH)
    accepted_sources = set(
        registration.loc[registration["accepted_final"], "source_number"].astype(int)
    )
    for contrast, frame in overlays.items():
        unexpected = set(frame["source_number"].astype(int)) - accepted_sources
        if unexpected:
            raise ValueError(
                f"{contrast} overlay contains rejected sources: {sorted(unexpected)}"
            )

    log("Computing voxel-weighted domain and subdomain scores")
    scores = aggregate_scores(overlays, regions)

    log("Loading atlas volumes and rendering 66 delivered-section contexts")
    annotation, _ = nrrd.read(str(ANNOTATION_PATH))
    nissl, _ = nrrd.read(str(NISSL_PATH))
    lut = pd.read_csv(ALLEN_LUT_PATH, header=None).to_numpy()
    frames = atlas_frames(
        registration,
        slice_manifest,
        annotation,
        nissl,
        lut,
        regions,
        overlays,
    )

    log("Compacting Allen meshes for linked 3D selection")
    meshes = compact_meshes()
    expected_meshes = {"root"} | {region["key"] for region in regions}
    if set(meshes) != expected_meshes:
        raise ValueError("Dashboard mesh keys differ from taxonomy keys")

    manifest = build_manifest(regions, scores, registration, frames)
    payload = clean_json(
        {
            "analysis_date": ANALYSIS_DATE,
            "regions": regions,
            "domains": [r for r in regions if r["level"] == "domain"],
            "scores": scores.to_dict(orient="records"),
            "frames": frames,
            "points": {
                contrast: compact_points(frame)
                for contrast, frame in overlays.items()
            },
            "meshes": meshes,
            "manifest": manifest,
        }
    )
    payload_text = json.dumps(payload, separators=(",", ":"))
    html = (
        HTML_TEMPLATE.replace("__PLOTLY_JS__", PLOTLY_JS.read_text())
        .replace("__PAYLOAD__", payload_text)
    )
    if "__PAYLOAD__" in html or "__PLOTLY_JS__" in html:
        raise ValueError("Dashboard template replacement failed")

    log("Writing new dated dashboard artifacts")
    scores_temp = SCORES_OUT.with_name(f".tmp-{SCORES_OUT.name}")
    manifest_temp = MANIFEST_OUT.with_name(f".tmp-{MANIFEST_OUT.name}")
    html_temp = HTML_OUT.with_name(f".tmp-{HTML_OUT.name}")
    scores.to_csv(scores_temp, index=False)
    manifest_temp.write_text(json.dumps(clean_json(manifest), indent=2))
    html_temp.write_text(html)
    for temp, final in (
        (scores_temp, SCORES_OUT),
        (manifest_temp, MANIFEST_OUT),
        (html_temp, HTML_OUT),
    ):
        if final.exists():
            raise FileExistsError(f"Refusing to overwrite: {final}")
        temp.rename(final)
    log(f"Wrote {SCORES_OUT}")
    log(f"Wrote {MANIFEST_OUT}")
    log(f"Wrote {HTML_OUT} ({HTML_OUT.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
