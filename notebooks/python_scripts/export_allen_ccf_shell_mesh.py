#!/usr/bin/env python
"""Export a plotly-sized translucent Allen CCF whole-brain shell mesh.

Source: local Allen `annotation_25.nrrd` (non-zero tissue mask → marching cubes).
Also archives the official Allen root structure OBJ (id 997) for provenance.

Outputs under data/interim/2DG_CCF_Integration/<ANALYSIS_DATE>/:
  allen_ccf_shell_mesh_plotly.npz
  allen_ccf_shell_vertices.csv.gz
  allen_ccf_shell_faces.csv.gz
  allen_ccf_shell_edges.csv.gz
  allen_ccf_shell_mesh_provenance.json

Usage:
  /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \\
    notebooks/python_scripts/export_allen_ccf_shell_mesh.py
"""

from __future__ import annotations

import json
from pathlib import Path

import nrrd
import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import measure

WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
ANALYSIS_DATE = "2026-11-14d"
OUT_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
ANNOTATION_PATH = WKDIR / "data/input/allen_registration_ref/annotation_25.nrrd"
ALLEN_ROOT_OBJ = OUT_DIR / "allen_ccf_root_997.obj"

# 200 µm voxels → ~14k verts / ~28k faces — interactive in Plotly
DOWNSAMPLE_STEP = 8


def log(msg: str) -> None:
    print(f"[ccf-shell] {msg}", flush=True)


def refuse(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npz_path = OUT_DIR / "allen_ccf_shell_mesh_plotly.npz"
    vert_path = OUT_DIR / "allen_ccf_shell_vertices.csv.gz"
    face_path = OUT_DIR / "allen_ccf_shell_faces.csv.gz"
    edge_path = OUT_DIR / "allen_ccf_shell_edges.csv.gz"
    prov_path = OUT_DIR / "allen_ccf_shell_mesh_provenance.json"
    for p in (npz_path, vert_path, face_path, edge_path, prov_path):
        refuse(p)

    log(f"Loading {ANNOTATION_PATH}")
    ann, header = nrrd.read(str(ANNOTATION_PATH))
    mask = (ann > 0)[::DOWNSAMPLE_STEP, ::DOWNSAMPLE_STEP, ::DOWNSAMPLE_STEP]
    mask = ndimage.binary_closing(mask, iterations=1)
    verts, faces, normals, _ = measure.marching_cubes(
        mask.astype(np.float32), level=0.5, step_size=1, allow_degenerate=False
    )
    scale = 25.0 * DOWNSAMPLE_STEP
    V = (verts * scale).astype(np.float32)  # columns: AP, DV, ML µm
    F = faces.astype(np.int32)
    edges = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0).astype(np.int32)

    np.savez_compressed(
        npz_path,
        vertices_ap_dv_ml_um=V,
        faces=F,
        edges=edges,
        normals=normals.astype(np.float32),
        voxel_um=scale,
        downsample_step=DOWNSAMPLE_STEP,
        source="annotation_25.nrrd non-zero mask",
        plotly_map="x=AP, y=ML, z=DV",
    )
    pd.DataFrame({"ap": V[:, 0], "dv": V[:, 1], "ml": V[:, 2]}).to_csv(
        vert_path, index=False
    )
    pd.DataFrame({"i": F[:, 0], "j": F[:, 1], "k": F[:, 2]}).to_csv(
        face_path, index=False
    )
    pd.DataFrame({"a": edges[:, 0], "b": edges[:, 1]}).to_csv(edge_path, index=False)

    prov = {
        "analysis_date": ANALYSIS_DATE,
        "annotation_path": str(ANNOTATION_PATH),
        "allen_root_obj_archived": str(ALLEN_ROOT_OBJ) if ALLEN_ROOT_OBJ.exists() else None,
        "allen_root_structure_id": 997,
        "downsample_step": DOWNSAMPLE_STEP,
        "voxel_um": scale,
        "n_vertices": int(len(V)),
        "n_faces": int(len(F)),
        "n_edges": int(len(edges)),
        "plotly_coordinates": {
            "x": "AP µm",
            "y": "ML µm (annotation axis 2)",
            "z": "DV µm (annotation axis 1, autorange reversed in Plotly)",
        },
        "aesthetic": (
            "Translucent whole-brain shell with optional wire edges; "
            "tunable opacity to match Allen Brain Explorer–style outline views."
        ),
        "nrrd_header_keys": sorted(header.keys()),
    }
    prov_path.write_text(json.dumps(prov, indent=2))
    log(f"Wrote {npz_path} ({len(V)} verts, {len(F)} faces)")


if __name__ == "__main__":
    main()
