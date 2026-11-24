#!/usr/bin/env python3
"""Build a self-contained two-pane Template-to-Allen manual matching browser.

The left pane reproduces panel 1 from the registration QC browser: original
Template anatomy with its detected tissue edge. The right pane reproduces
panel 3: Allen Nissl anatomy with the atlas tissue edge. Unlike the QC browser,
the Allen pane includes every coronal plane in ``ara_nissl_25.nrrd`` so a user
can manually inspect and record the best AP/bregma match for each delivered
Template section.

Mappings start blank. Template and Allen panes scroll independently unless the
user enables linked review after saving matches.

Output:
  figures/2DG_CCF_Integration/
    2DG_Template_Allen_manual_matcher_2026-11-23b.html

Run with:
  /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
    notebooks/python_scripts/build_2dg_template_allen_manual_matcher.py
"""

from __future__ import annotations

import base64
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from skimage.segmentation import find_boundaries

from build_2dg_registration_qc_browser import grayscale, load_atlas, qc_status


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
SOURCE_DATE = "2026-11-11"
ANALYSIS_DATE = "2026-11-23b"
PRIOR_HTML = (
    WKDIR
    / "figures/2DG_CCF_Integration"
    / "2DG_Template_Allen_manual_matcher_2026-11-23.html"
)
SOURCE_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / SOURCE_DATE
FIGURE_DIR = WKDIR / "figures/2DG_CCF_Integration"
REGISTRATION_PATH = (
    SOURCE_DIR / f"2dg_registration_summary_strict_{SOURCE_DATE}.csv"
)
HTML_OUT = (
    FIGURE_DIR
    / f"2DG_Template_Allen_manual_matcher_{ANALYSIS_DATE}.html"
)

BREGMA_INTERCEPT_MM = 5.4
BREGMA_STEP_MM = 0.025


def log(message: str) -> None:
    print(f"[2DG-manual-matcher] {message}", flush=True)


def refuse_overwrite(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")


def image_data_url(rgb: np.ndarray, quality: int) -> str:
    array = np.clip(np.round(rgb * 255), 0, 255).astype(np.uint8)
    buffer = io.BytesIO()
    Image.fromarray(array, mode="RGB").save(
        buffer,
        format="WEBP",
        quality=quality,
        method=6,
    )
    return "data:image/webp;base64," + base64.b64encode(
        buffer.getvalue()
    ).decode("ascii")


def template_image(slice_position: int) -> str:
    prepared_path = (
        SOURCE_DIR / "prepared" / f"slice_{slice_position:03d}.npz"
    )
    prepared = np.load(prepared_path)
    template = prepared["template_y"].astype(np.float32)
    tissue = prepared["tissue"].astype(bool)
    gray = grayscale(template, positive_only=True)
    rgb = np.stack([gray, gray, gray], axis=-1)
    edge = find_boundaries(tissue, mode="outer")
    rgb[edge] = np.array([0.0, 0.66, 0.71], dtype=np.float32)
    return image_data_url(rgb, quality=92)


def allen_image(
    nissl: np.ndarray,
    annotation: np.ndarray,
    ap_index: int,
) -> str:
    plane = nissl[ap_index].astype(np.float32)
    tissue = annotation[ap_index] > 0
    gray = grayscale(
        np.log1p(np.maximum(plane, 0)),
        positive_only=True,
    )
    rgb = np.stack([gray, gray, gray], axis=-1)
    edge = find_boundaries(tissue, mode="outer")
    rgb[edge] = np.array([0.18, 0.18, 0.18], dtype=np.float32)
    return image_data_url(rgb, quality=86)


def bregma_for_ap(ap_index: int) -> float:
    return BREGMA_INTERCEPT_MM - BREGMA_STEP_MM * ap_index


def load_prior_payload(path: Path) -> dict | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8")
    match = re.search(r"const DATA=(\{.*?\});\s*let templateIndex", text, re.S)
    if not match:
        return None
    payload = json.loads(match.group(1))
    if len(payload.get("templates", [])) != 66:
        return None
    if len(payload.get("allen", [])) != 528:
        return None
    return payload


def build_payload() -> dict:
    prior = load_prior_payload(PRIOR_HTML)
    registration = (
        pd.read_csv(REGISTRATION_PATH)
        .sort_values("slice_position")
        .reset_index(drop=True)
    )
    if len(registration) != 66:
        raise ValueError(
            f"Expected 66 delivered Template sections, found {len(registration)}"
        )

    if prior is not None:
        log(f"Reusing encoded images from {PRIOR_HTML.name}")
        templates = []
        for index, row in registration.iterrows():
            prior_row = prior["templates"][int(index)]
            if int(prior_row["source_number"]) != int(row["source_number"]):
                raise ValueError("Prior Template order does not match registration")
            status_key, status_label = qc_status(row)
            templates.append(
                {
                    "index": int(index),
                    "template_index": int(row["source_number"]),
                    "slice_position": int(row["slice_position"]),
                    "source_number": int(row["source_number"]),
                    "status_key": status_key,
                    "status_label": status_label,
                    "accepted_final": bool(row["accepted_final"]),
                    "image": prior_row["image"],
                }
            )
        allen = [
            {
                "ap_index": int(item["ap_index"]),
                "bregma_mm": round(float(item["bregma_mm"]), 3),
                "image": item["image"],
            }
            for item in prior["allen"]
        ]
    else:
        nissl, annotation = load_atlas()
        if nissl.shape != annotation.shape:
            raise ValueError(
                f"Nissl shape {nissl.shape} != annotation shape {annotation.shape}"
            )
        log("Encoding 66 original Template anatomy images")
        templates = []
        for index, row in registration.iterrows():
            status_key, status_label = qc_status(row)
            templates.append(
                {
                    "index": int(index),
                    "template_index": int(row["source_number"]),
                    "slice_position": int(row["slice_position"]),
                    "source_number": int(row["source_number"]),
                    "status_key": status_key,
                    "status_label": status_label,
                    "accepted_final": bool(row["accepted_final"]),
                    "image": template_image(int(row["slice_position"])),
                }
            )
        log(f"Encoding all {nissl.shape[0]} Allen coronal planes")
        allen = []
        for ap_index in range(nissl.shape[0]):
            allen.append(
                {
                    "ap_index": ap_index,
                    "bregma_mm": round(bregma_for_ap(ap_index), 3),
                    "image": allen_image(nissl, annotation, ap_index),
                }
            )
            if (ap_index + 1) % 50 == 0:
                log(f"Encoded {ap_index + 1}/{nissl.shape[0]} Allen planes")

    return {
        "manifest": {
            "analysis_date": ANALYSIS_DATE,
            "source_registration": SOURCE_DATE,
            "n_templates": len(templates),
            "n_allen_planes": len(allen),
            "allen_ap_min": 0,
            "allen_ap_max": len(allen) - 1,
            "bregma_formula": "bregma_mm = 5.400 - 0.025 * allen_ap_index",
            "bregma_anterior_mm": bregma_for_ap(0),
            "bregma_posterior_mm": bregma_for_ap(len(allen) - 1),
            "mappings_start_blank": True,
        },
        "templates": templates,
        "allen": allen,
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>2DG Template ↔ Allen full-coronal manual matcher</title>
<style>
:root{--ink:#181a1e;--muted:#626a75;--line:#dfe3e8;--soft:#f5f6f8;--accent:#ff6c0c;--teal:#00a7b5;--green:#16845d}
*{box-sizing:border-box}html,body{margin:0;height:100%;overflow:hidden;background:#f5f6f8;color:var(--ink);font-family:Inter,Helvetica,Arial,sans-serif}
button,select,input{font:inherit}
.app{height:100%;display:grid;grid-template-rows:auto auto 1fr;gap:8px;padding:10px}
.top{display:flex;justify-content:space-between;gap:18px;align-items:start;background:#fff;border:1px solid var(--line);border-radius:11px;padding:10px 13px}
h1{font-size:18px;margin:0 0 3px;letter-spacing:-.02em}.sub{font-size:11px;line-height:1.4;color:var(--muted);margin:0;max-width:980px}
.range-note{font-size:11px;text-align:right;color:var(--muted);white-space:nowrap}.range-note strong{display:block;color:var(--ink);font-size:13px}
.toolbar{display:grid;grid-template-columns:1.2fr auto 1fr;gap:10px;align-items:center;background:#fff;border:1px solid var(--line);border-radius:11px;padding:8px 12px}
.summary{font-size:11px;color:var(--muted)}.summary strong{display:block;color:var(--ink);font-size:14px;margin-bottom:2px}
.link-box{display:flex;flex-direction:column;gap:6px;min-width:280px}
.toggle{display:flex;align-items:center;gap:8px;font-size:12px;font-weight:700}
.toggle input{width:16px;height:16px;accent-color:var(--accent)}
.pair-row{display:flex;align-items:center;gap:6px}.pair-row input[type=range]{flex:1;accent-color:var(--accent)}
.pair-row[disabled]{opacity:.45;pointer-events:none}
.save-tools{display:flex;align-items:center;justify-content:flex-end;gap:6px;flex-wrap:wrap}
button,select,input[type=number]{border:1px solid #cfd5dc;border-radius:7px;background:#fff;color:var(--ink);padding:6px 8px;font-size:11px}
button{cursor:pointer}button:hover{border-color:var(--accent);background:#fff5ee}.primary{background:#fff4ec;border-color:#f0a776;color:#8b3500;font-weight:700}
.progress{font-size:10px;color:var(--green);font-weight:700;min-width:74px;text-align:right}
.workspace{display:grid;grid-template-columns:minmax(220px,260px) 1fr 1fr;gap:8px;min-height:0}
.sidebar{display:grid;grid-template-rows:auto 1fr auto;background:#fff;border:1px solid var(--line);border-radius:12px;overflow:hidden;min-height:0}
.sidebar-head{padding:10px 11px;border-bottom:1px solid var(--line)}.sidebar-head h2{font-size:14px;margin:0 0 2px}.sidebar-head p{font-size:10.5px;color:var(--muted);margin:0}
.table-wrap{overflow:auto;min-height:0}.map-table{width:100%;border-collapse:collapse;font-size:11px}
.map-table th,.map-table td{border-bottom:1px solid var(--line);padding:5px 7px;text-align:left;white-space:nowrap}
.map-table th{position:sticky;top:0;background:#fafbfc;z-index:1;font-size:10px;color:var(--muted);font-weight:700}
.map-table tr{cursor:pointer}.map-table tr:hover{background:#fff7f1}.map-table tr.active{background:#fff0e6}
.map-table input{width:100%;min-width:78px;padding:4px 5px}
.sidebar-foot{padding:8px 10px;border-top:1px solid var(--line);font-size:10px;color:var(--muted)}
.pane{display:grid;grid-template-rows:auto auto 1fr auto;background:#fff;border:1px solid var(--line);border-radius:12px;overflow:hidden;min-width:0;min-height:0}
.pane-head{display:flex;justify-content:space-between;align-items:start;gap:10px;padding:9px 11px;border-bottom:1px solid var(--line)}
.pane-head h2{font-size:14px;margin:0 0 2px}.pane-head p{font-size:10.5px;color:var(--muted);margin:0}.status{font-size:10px;border:1px solid var(--line);border-radius:999px;padding:4px 7px;white-space:nowrap}
.status.accepted{color:#16845d;background:#f0faf6;border-color:#b7dfd0}.status.rejected{color:#9b4612;background:#fff7f1;border-color:#f3d0b8}.status.blank{color:#626a75;background:#f5f6f8;border-color:#dfe3e8}
.big-label{text-align:center;padding:7px 10px;border-bottom:1px solid var(--line);background:#fafbfc}
.big-label .value{font-size:25px;font-weight:800;letter-spacing:-.025em}.big-label .detail{font-size:11px;color:var(--muted);margin-top:1px}
.image-wrap{min-height:0;display:grid;place-items:center;padding:8px;background:#090a0c;overflow:hidden;overscroll-behavior:contain}
.image-wrap img{display:block;max-width:100%;max-height:100%;width:auto;height:auto;object-fit:contain}
.controls{padding:8px 10px;border-top:1px solid var(--line);background:#fff}.control-row{display:flex;align-items:center;gap:6px}
.control-row+ .control-row{margin-top:6px}.control-row input[type=range]{flex:1;accent-color:var(--accent);min-width:90px}.control-row select{flex:1;min-width:0}
.axis-labels{display:flex;justify-content:space-between;color:var(--muted);font-size:9.5px;margin-top:3px}
.hint{font-size:10px;color:var(--muted)}.manual-saved{color:var(--green);font-weight:700}
@media(max-width:1100px){html,body{overflow:auto}.app{height:auto}.workspace{grid-template-columns:1fr}.sidebar,.pane{height:560px}.toolbar{grid-template-columns:1fr}.save-tools{justify-content:flex-start}.range-note{text-align:left}}
</style>
</head>
<body>
<main class="app">
  <header class="top">
    <div>
      <h1>Template ↔ full Allen coronal atlas matcher</h1>
      <p class="sub">Start from a blank mapping table. Scroll Template and Allen panes independently, then press <b>Link current bregma → template</b> to save a match. Enable <b>Linked review</b> only after you have matches, to flip both panes together with the third scrollbar.</p>
    </div>
    <div class="range-note"><strong>528 Allen planes · 25 µm spacing</strong>Anterior +5.400 mm → posterior −7.775 mm</div>
  </header>
  <section class="toolbar">
    <div class="summary"><strong id="summaryTemplate">Template 001</strong><span id="summaryStatus">No manual match yet</span></div>
    <div class="link-box">
      <label class="toggle"><input id="linkToggle" type="checkbox"> Linked review of saved mappings</label>
      <div class="pair-row" id="pairRow" disabled>
        <button id="pairPrev" type="button">◀</button>
        <input id="pairSlider" type="range" min="0" max="65" value="0" aria-label="Linked mapping review">
        <button id="pairNext" type="button">▶</button>
      </div>
      <div class="hint" id="pairHint">Independent scrolling is on. Linked review stays off until you enable it.</div>
    </div>
    <div class="save-tools">
      <span class="progress" id="progress">0 / 66 saved</span>
      <button id="saveMatch" type="button" class="primary">Link current bregma → template</button>
      <button id="clearMatch" type="button">Clear this template</button>
      <button id="downloadCsv" type="button">Download CSV</button>
    </div>
  </section>
  <section class="workspace">
    <aside class="sidebar">
      <div class="sidebar-head">
        <h2>Manual map table</h2>
        <p>Blank start. Edit bregma directly, or save from the Allen pane.</p>
      </div>
      <div class="table-wrap">
        <table class="map-table">
          <thead><tr><th>template_index</th><th>bregma</th></tr></thead>
          <tbody id="mapBody"></tbody>
        </table>
      </div>
      <div class="sidebar-foot">CSV export columns: <code>template_index,bregma</code></div>
    </aside>
    <article class="pane" id="templatePane">
      <div class="pane-head">
        <div><h2>1 · Original Template anatomy</h2><p>Wheel here or use A / D. Does not move Allen unless Linked review is on.</p></div>
        <span class="status" id="templateStatus"></span>
      </div>
      <div class="big-label">
        <div class="value" id="templateValue">Template 001</div>
        <div class="detail" id="templateDetail"></div>
      </div>
      <div class="image-wrap" id="templateWheel"><img id="templateImage" alt="Original Template anatomy"></div>
      <div class="controls">
        <div class="control-row">
          <button id="templatePrev" type="button">◀</button>
          <input id="templateSlider" type="range" min="0" max="65" value="0" aria-label="Template section">
          <button id="templateNext" type="button">▶</button>
          <select id="templateSelect" aria-label="Jump to Template section"></select>
        </div>
        <div class="axis-labels"><span>Template 001 · anterior</span><span>Template 066 · posterior</span></div>
      </div>
    </article>
    <article class="pane" id="allenPane">
      <div class="pane-head">
        <div><h2>3 · Allen reference plane</h2><p>Wheel here or use ← / → through every atlas plane.</p></div>
        <span class="status accepted">Full atlas</span>
      </div>
      <div class="big-label">
        <div class="value" id="allenValue">Bregma +5.400 mm</div>
        <div class="detail" id="allenDetail">Allen AP index 0 / 527</div>
      </div>
      <div class="image-wrap" id="allenWheel"><img id="allenImage" alt="Allen Nissl reference plane"></div>
      <div class="controls">
        <div class="control-row">
          <button id="allenPrev" type="button">◀</button>
          <input id="allenSlider" type="range" min="0" max="527" value="0" aria-label="Allen AP plane">
          <button id="allenNext" type="button">▶</button>
          <label class="hint" for="bregmaInput">Bregma mm</label>
          <input id="bregmaInput" type="number" min="-7.775" max="5.400" step="0.025" value="5.400">
        </div>
        <div class="axis-labels"><span>AP 0 · anterior · +5.400 mm</span><span>AP 527 · posterior · −7.775 mm</span></div>
      </div>
    </article>
  </section>
</main>
<script>
const DATA=__PAYLOAD__;
let templateIndex=0;
let allenIndex=0;
let linkedReview=false;
let matches={};
const storageKey=`2dg-template-allen-manual-matches-blank-${DATA.manifest.analysis_date}`;
try{matches=JSON.parse(localStorage.getItem(storageKey)||"{}")}catch(_){matches={}}
const $=id=>document.getElementById(id);

function clamp(value,min,max){return Math.max(min,Math.min(max,value))}
function signed(value,digits=3){return `${value>=0?"+":""}${Number(value).toFixed(digits)}`}
function currentTemplate(){return DATA.templates[templateIndex]}
function currentAllen(){return DATA.allen[allenIndex]}
function savedFor(template){return matches[String(template.template_index)]||null}
function savedForCurrent(){return savedFor(currentTemplate())}
function persist(){try{localStorage.setItem(storageKey,JSON.stringify(matches))}catch(_){}}
function apFromBregma(bregma){return clamp(Math.round((5.4-Number(bregma))/0.025),0,DATA.allen.length-1)}
function statusClass(template){return template.accepted_final?"accepted":(template.status_key==="blank"?"blank":"rejected")}

function setMatch(template,bregma){
  const ap=apFromBregma(bregma);
  const allen=DATA.allen[ap];
  matches[String(template.template_index)]={
    template_index:template.template_index,
    bregma:Number(allen.bregma_mm.toFixed(3)),
    allen_ap_index:allen.ap_index
  };
  persist();
}

function clearMatch(template){
  delete matches[String(template.template_index)];
  persist();
}

function renderTable(){
  $("mapBody").innerHTML=DATA.templates.map((t,i)=>{
    const saved=savedFor(t);
    const value=saved?Number(saved.bregma).toFixed(3):"";
    return `<tr data-index="${i}" class="${i===templateIndex?"active":""}"><td>${t.template_index}</td><td><input data-index="${i}" type="number" step="0.025" min="-7.775" max="5.400" value="${value}" placeholder="—"></td></tr>`;
  }).join("");
  $("mapBody").querySelectorAll("tr").forEach(row=>{
    row.onclick=e=>{if(e.target.tagName==="INPUT")return;templateIndex=Number(row.dataset.index);renderTemplate();};
  });
  $("mapBody").querySelectorAll("input").forEach(input=>{
    input.onchange=()=>{
      const idx=Number(input.dataset.index);
      const template=DATA.templates[idx];
      if(input.value===""){clearMatch(template)}
      else{setMatch(template,input.value)}
      if(idx===templateIndex){
        const saved=savedFor(template);
        if(saved)allenIndex=saved.allen_ap_index;
        renderAll();
      }else{renderProgress();renderTable()}
    };
  });
}

function renderProgress(){
  const n=Object.keys(matches).length;
  $("progress").textContent=`${n} / ${DATA.templates.length} saved`;
  $("pairHint").textContent=linkedReview
    ?"Linked review on: the third scrollbar flips Template and its saved Allen plane together."
    :"Independent scrolling is on. Linked review stays off until you enable it.";
}

function renderTemplate(){
  const t=currentTemplate();
  const saved=savedForCurrent();
  $("templateImage").src=t.image;
  $("templateSlider").value=String(templateIndex);
  $("templateSelect").value=String(templateIndex);
  $("pairSlider").value=String(templateIndex);
  $("templateValue").textContent=`Template ${String(t.template_index).padStart(3,"0")}`;
  $("templateDetail").innerHTML=`Delivered ${templateIndex+1}/${DATA.templates.length} · ${t.status_label}${saved?` · <span class="manual-saved">manual bregma ${signed(saved.bregma)}</span>`:" · no manual match yet"}`;
  $("templateStatus").textContent=t.status_label;
  $("templateStatus").className=`status ${statusClass(t)}`;
  $("summaryTemplate").textContent=`Template ${String(t.template_index).padStart(3,"0")}`;
  $("summaryStatus").textContent=saved?`Manual match saved at bregma ${signed(saved.bregma)} mm`:"No manual match yet";
  if(linkedReview && saved)allenIndex=saved.allen_ap_index;
  renderAllen();
  renderProgress();
  renderTable();
}

function renderAllen(){
  allenIndex=clamp(allenIndex,0,DATA.allen.length-1);
  const a=currentAllen();
  $("allenImage").src=a.image;
  $("allenSlider").value=String(allenIndex);
  $("allenValue").textContent=`Bregma ${signed(a.bregma_mm)} mm`;
  $("allenDetail").textContent=`Allen AP index ${a.ap_index} / ${DATA.allen.length-1} · 25 µm coronal spacing`;
  $("bregmaInput").value=Number(a.bregma_mm).toFixed(3);
}

function renderAll(){renderTemplate()}
function moveTemplate(delta){templateIndex=clamp(templateIndex+delta,0,DATA.templates.length-1);renderTemplate()}
function moveAllen(delta){allenIndex=clamp(allenIndex+delta,0,DATA.allen.length-1);renderAllen()}
function movePair(delta){templateIndex=clamp(templateIndex+delta,0,DATA.templates.length-1);renderTemplate()}
function wheelStep(event,move){event.preventDefault();move(event.deltaY>0?1:-1)}
function syncLinkControls(){
  linkedReview=$("linkToggle").checked;
  $("pairRow").toggleAttribute("disabled",!linkedReview);
  renderProgress();
  if(linkedReview)renderTemplate();
}

$("templatePrev").onclick=()=>moveTemplate(-1);
$("templateNext").onclick=()=>moveTemplate(1);
$("templateSlider").oninput=e=>{templateIndex=Number(e.target.value);renderTemplate()};
$("templateSelect").onchange=e=>{templateIndex=Number(e.target.value);renderTemplate()};
$("allenPrev").onclick=()=>moveAllen(-1);
$("allenNext").onclick=()=>moveAllen(1);
$("allenSlider").oninput=e=>{allenIndex=Number(e.target.value);renderAllen()};
$("bregmaInput").onchange=e=>{allenIndex=apFromBregma(e.target.value);renderAllen()};
$("templateWheel").addEventListener("wheel",e=>wheelStep(e,moveTemplate),{passive:false});
$("allenWheel").addEventListener("wheel",e=>wheelStep(e,moveAllen),{passive:false});
$("linkToggle").onchange=syncLinkControls;
$("pairPrev").onclick=()=>movePair(-1);
$("pairNext").onclick=()=>movePair(1);
$("pairSlider").oninput=e=>{templateIndex=Number(e.target.value);renderTemplate()};
$("saveMatch").onclick=()=>{setMatch(currentTemplate(),currentAllen().bregma_mm);renderAll()};
$("clearMatch").onclick=()=>{clearMatch(currentTemplate());renderAll()};
$("downloadCsv").onclick=()=>{
  const rows=DATA.templates.map(t=>{
    const saved=savedFor(t);
    return {template_index:t.template_index,bregma:saved?Number(saved.bregma).toFixed(3):""};
  });
  const csv=["template_index,bregma",...rows.map(r=>`${r.template_index},${r.bregma}`)].join("\n");
  const url=URL.createObjectURL(new Blob([csv],{type:"text/csv"}));
  const link=document.createElement("a");link.href=url;
  link.download=`2DG_Template_Allen_manual_matches_${DATA.manifest.analysis_date}.csv`;
  link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
};
window.addEventListener("keydown",event=>{
  if(["INPUT","SELECT","TEXTAREA"].includes(document.activeElement.tagName))return;
  if(event.key==="ArrowLeft"){event.preventDefault();moveAllen(-1)}
  if(event.key==="ArrowRight"){event.preventDefault();moveAllen(1)}
  if(event.key.toLowerCase()==="a")moveTemplate(-1);
  if(event.key.toLowerCase()==="d")moveTemplate(1);
  if(event.key==="Enter"){event.preventDefault();$("saveMatch").click()}
});
$("templateSelect").innerHTML=DATA.templates.map((t,i)=>`<option value="${i}">Template ${String(t.template_index).padStart(3,"0")} · ${t.status_label}</option>`).join("");
syncLinkControls();
renderAll();
</script>
</body>
</html>
"""


def write_html(payload: dict, path: Path) -> None:
    compact = json.dumps(payload, separators=(",", ":"))
    path.write_text(
        HTML_TEMPLATE.replace("__PAYLOAD__", compact),
        encoding="utf-8",
    )


def main() -> None:
    refuse_overwrite(HTML_OUT)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    write_html(payload, HTML_OUT)
    size_mb = HTML_OUT.stat().st_size / 1_000_000
    log(
        f"Wrote {HTML_OUT} ({size_mb:.1f} MB; "
        f"{payload['manifest']['n_templates']} Templates; "
        f"{payload['manifest']['n_allen_planes']} Allen planes; "
        "blank mapping table)"
    )


if __name__ == "__main__":
    main()
