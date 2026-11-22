#!/usr/bin/env python3
"""Export compact Allen CCF domain/subdomain meshes for the modular 2DG viewer.

This adapts the supplied BrainGlobe/brainrender reference aesthetic to the
project's local ``annotation_25.nrrd`` so no additional atlas packages are
required. Parent structures include all descendant annotation IDs from the
local Allen hierarchy.

Outputs are date-stamped and never overwrite existing artifacts:

  data/interim/2DG_CCF_Integration/2026-11-21e/
    allen_ccf_region_mesh_manifest_2026-11-21e.json
    allen_ccf_region_mesh_vertices_2026-11-21e.csv.gz
    allen_ccf_region_mesh_faces_2026-11-21e.csv.gz
    allen_ccf_region_mesh_edges_2026-11-21e.csv.gz

Run with:

  /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
    notebooks/python_scripts/export_2dg_ccf_region_meshes.py
"""

from __future__ import annotations

import json
import os
import pickle
from pathlib import Path

import nrrd
import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import measure


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
ANALYSIS_DATE = "2026-11-21e"
ATLAS_DIR = WKDIR / "data/input/allen_registration_ref"
OUT_DIR = WKDIR / "data/interim/2DG_CCF_Integration" / ANALYSIS_DATE
ANNOTATION_PATH = ATLAS_DIR / "annotation_25.nrrd"
HIERARCHY_PATH = ATLAS_DIR / "allen_name_to_annots.pkl"

MANIFEST_OUT = OUT_DIR / f"allen_ccf_region_mesh_manifest_{ANALYSIS_DATE}.json"
VERTICES_OUT = OUT_DIR / f"allen_ccf_region_mesh_vertices_{ANALYSIS_DATE}.csv.gz"
FACES_OUT = OUT_DIR / f"allen_ccf_region_mesh_faces_{ANALYSIS_DATE}.csv.gz"
EDGES_OUT = OUT_DIR / f"allen_ccf_region_mesh_edges_{ANALYSIS_DATE}.csv.gz"

VOXEL_UM = 25.0
MAX_ROOT_EDGES = 8_000
MAX_DOMAIN_EDGES = 4_000
MAX_SUBDOMAIN_EDGES = 2_000
SEED = 42


# Label, Allen hierarchy name. Domain masks use the named parent or, for
# Amygdala, the union of its curated subdomains because Allen splits amygdala
# across cortical-subplate and striatal branches.
TAXONOMY: dict[str, dict] = {
    "Isocortex": {
        "source": "Isocortex",
        "color": "#7bc8b5",
        "subdomains": {
            "MOp": ("Primary motor area", "#bce4db"),
            "MOs": ("Secondary motor area", "#a7dacf"),
            "SSp": ("Primary somatosensory area", "#8dcec0"),
            "ACA": ("Anterior cingulate area", "#69bbaa"),
            "PL": ("Prelimbic area", "#54aa9a"),
            "ILA": ("Infralimbic area", "#429989"),
            "ORB": ("Orbital area", "#2f8778"),
            "AI": ("Agranular insular area", "#237467"),
            "RSP": ("Retrosplenial area", "#176256"),
            "VISp": ("Primary visual area", "#0f5046"),
            "AUDp": ("Primary auditory area", "#093f37"),
            "ENTl": ("Entorhinal area, lateral part", "#06302a"),
        },
    },
    "Olfactory areas": {
        "source": "Olfactory areas",
        "color": "#76c9d6",
        "subdomains": {
            "MOB": ("Main olfactory bulb", "#b7e4ea"),
            "AOB": ("Accessory olfactory bulb", "#8ed4de"),
            "AON": ("Anterior olfactory nucleus", "#64c1cf"),
            "PIR": ("Piriform area", "#3cabbc"),
            "OT": ("Olfactory tubercle", "#238e9f"),
        },
    },
    "Hippocampal formation": {
        "source": "Hippocampal formation",
        "color": "#aa87d8",
        "subdomains": {
            "HIP": ("Hippocampal region", "#d7c4eb"),
            "CA1": ("Field CA1", "#c3a7e3"),
            "CA3": ("Field CA3", "#ad8ad8"),
            "DG": ("Dentate gyrus", "#906ac5"),
            "SUB": ("Subiculum", "#714ba8"),
        },
    },
    "Amygdala": {
        "source": None,
        "color": "#d1495b",
        "subdomains": {
            "LA": ("Lateral amygdalar nucleus", "#f4a261"),
            "BLA": ("Basolateral amygdalar nucleus", "#e76f51"),
            "BMA": ("Basomedial amygdalar nucleus", "#d1495b"),
            "CEA": ("Central amygdalar nucleus", "#9b2226"),
            "MEA": ("Medial amygdalar nucleus", "#7f1525"),
            "COA": ("Cortical amygdalar area", "#c9184a"),
            "PA": ("Posterior amygdalar nucleus", "#721329"),
        },
    },
    "Striatum": {
        "source": "Striatum",
        "color": "#d7bd58",
        "subdomains": {
            "CP": ("Caudoputamen", "#f3df8d"),
            "ACB": ("Nucleus accumbens", "#e7ce6c"),
            "FS": ("Fundus of striatum", "#d8b849"),
            "LSX": ("Lateral septal complex", "#c5a12e"),
            "STRa": ("Striatum-like amygdalar nuclei", "#a98520"),
        },
    },
    "Pallidum": {
        "source": "Pallidum",
        "color": "#bfad55",
        "subdomains": {
            "GPe": ("Globus pallidus, external segment", "#e3d789"),
            "GPi": ("Globus pallidus, internal segment", "#cabe68"),
            "BST": ("Bed nuclei of the stria terminalis", "#ad9844"),
            "SI": ("Substantia innominata", "#8e7b32"),
            "MS": ("Medial septal complex", "#716126"),
        },
    },
    "Thalamus": {
        "source": "Thalamus",
        "color": "#66b77c",
        "subdomains": {
            "LGd": ("Dorsal part of the lateral geniculate complex", "#a8d7ad"),
            "MG": ("Medial geniculate complex", "#8bcb98"),
            "VPM": ("Ventral posteromedial nucleus of the thalamus", "#6fbd83"),
            "VPL": ("Ventral posterolateral nucleus of the thalamus", "#55ad6e"),
            "MD": ("Mediodorsal nucleus of thalamus", "#3d995a"),
            "PVT": ("Paraventricular nucleus of the thalamus", "#2d8349"),
            "RT": ("Reticular nucleus of the thalamus", "#206e3c"),
            "LH": ("Lateral habenula", "#155a31"),
            "MH": ("Medial habenula", "#0b4727"),
        },
    },
    "Hypothalamus": {
        "source": "Hypothalamus",
        "color": "#3a86ff",
        "subdomains": {
            "LHA": ("Lateral hypothalamic area", "#90e0ef"),
            "ZI": ("Zona incerta", "#72d2ec"),
            "ARH": ("Arcuate hypothalamic nucleus", "#48cae4"),
            "MPO": ("Medial preoptic area", "#28b9dc"),
            "SO": ("Supraoptic nucleus", "#11a4cf"),
            "VMH": ("Ventromedial hypothalamic nucleus", "#0096c7"),
            "SCH": ("Suprachiasmatic nucleus", "#087eb3"),
            "PMv": ("Ventral premammillary nucleus", "#12679d"),
            "DMH": ("Dorsomedial nucleus of the hypothalamus", "#185188"),
            "MM": ("Medial mammillary nucleus", "#123d70"),
            "PVH": ("Paraventricular hypothalamic nucleus", "#023e8a"),
        },
    },
    "Midbrain": {
        "source": "Midbrain",
        "color": "#cc6ca8",
        "subdomains": {
            "VTA": ("Ventral tegmental area", "#efb0d1"),
            "SNc": ("Substantia nigra, compact part", "#e794c2"),
            "SNr": ("Substantia nigra, reticular part", "#da78b2"),
            "PAG": ("Periaqueductal gray", "#c75da0"),
            "SCs": ("Superior colliculus, sensory related", "#b14991"),
            "SCm": ("Superior colliculus, motor related", "#a23e83"),
            "IC": ("Inferior colliculus", "#923476"),
            "DR": ("Dorsal nucleus raphe", "#752861"),
            "RN": ("Red nucleus", "#5a1e4c"),
        },
    },
    "Pons": {
        "source": "Pons",
        "color": "#b79b87",
        "subdomains": {
            "LC": ("Locus ceruleus", "#dac8bb"),
            "PB": ("Parabrachial nucleus", "#c7ad9a"),
            "PG": ("Pontine gray", "#a98b76"),
        },
    },
    "Medulla": {
        "source": "Medulla",
        "color": "#907768",
        "subdomains": {
            "NTS": ("Nucleus of the solitary tract", "#c3ada0"),
            "DMX": ("Dorsal motor nucleus of the vagus nerve", "#a58a7a"),
            "AP": ("Area postrema", "#866c5e"),
            "IO": ("Inferior olivary complex", "#6e574c"),
        },
    },
    "Cerebellum": {
        "source": "Cerebellum",
        "color": "#9fc26b",
        "subdomains": {
            "CBX": ("Cerebellar cortex", "#c9dfa5"),
            "CBN": ("Cerebellar nuclei", "#91b858"),
        },
    },
}


def log(message: str) -> None:
    print(f"[ccf-region-mesh] {message}", flush=True)


def refuse_overwrite(*paths: Path) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing artifacts:\n" + "\n".join(existing))


def safe_key(value: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in value).strip("_")


def mesh_from_ids(
    annotation: np.ndarray,
    annotation_ids: set[int],
    *,
    step: int,
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isin(annotation, np.fromiter(annotation_ids, dtype=annotation.dtype))
    occupied = [
        np.flatnonzero(mask.any(axis=axes))
        for axes in ((1, 2), (0, 2), (0, 1))
    ]
    if any(not len(indices) for indices in occupied):
        raise ValueError("mask is empty")
    lo = np.maximum(np.array([indices[0] for indices in occupied]) - 2 * step, 0)
    hi = np.minimum(
        np.array([indices[-1] for indices in occupied]) + 2 * step + 1,
        annotation.shape,
    )
    crop = mask[lo[0] : hi[0] : step, lo[1] : hi[1] : step, lo[2] : hi[2] : step]
    crop = ndimage.binary_closing(crop, iterations=1)
    if min(crop.shape) < 3 or crop.all() or not crop.any():
        if step > 1:
            return mesh_from_ids(annotation, annotation_ids, step=max(1, step // 2))
        raise ValueError(f"mask cannot be meshed at crop shape {crop.shape}")
    verts, faces, _, _ = measure.marching_cubes(
        crop.astype(np.float32),
        level=0.5,
        step_size=1,
        allow_degenerate=False,
    )
    verts = (verts * (VOXEL_UM * step)) + (lo * VOXEL_UM)
    return verts.astype(np.float32), faces.astype(np.int32)


def sampled_edges(faces: np.ndarray, maximum: int, rng: np.random.Generator) -> np.ndarray:
    edges = np.vstack((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    if len(edges) > maximum:
        edges = edges[rng.choice(len(edges), size=maximum, replace=False)]
    return edges.astype(np.int32)


def main() -> None:
    refuse_overwrite(MANIFEST_OUT, VERTICES_OUT, FACES_OUT, EDGES_OUT)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log("Loading local Allen annotation and hierarchy")
    annotation, _ = nrrd.read(str(ANNOTATION_PATH))
    with HIERARCHY_PATH.open("rb") as handle:
        hierarchy: dict[str, list[int]] = pickle.load(handle)
    required_names = {
        source_name
        for spec in TAXONOMY.values()
        for source_name in (
            [spec["source"]] if spec["source"] is not None else []
        )
    }
    required_names.update(
        source_name
        for spec in TAXONOMY.values()
        for source_name, _ in spec["subdomains"].values()
    )
    missing_names = sorted(required_names - set(hierarchy))
    if missing_names:
        raise KeyError(
            "Allen hierarchy is missing required taxonomy entries:\n"
            + "\n".join(missing_names)
        )

    regions: list[dict] = []
    vertices_rows: list[pd.DataFrame] = []
    faces_rows: list[pd.DataFrame] = []
    edges_rows: list[pd.DataFrame] = []
    rng = np.random.default_rng(SEED)

    def add_region(
        key: str,
        label: str,
        group: str,
        level: str,
        color: str,
        ids: set[int],
        step: int,
        max_edges: int,
    ) -> bool:
        try:
            verts, faces = mesh_from_ids(annotation, ids, step=step)
        except ValueError as exc:
            log(f"SKIP {key}: {exc}")
            return False
        edges = sampled_edges(faces, max_edges, rng)
        vertices_rows.append(
            pd.DataFrame(
                {
                    "region": key,
                    "vertex_id": np.arange(len(verts), dtype=np.int32),
                    "ap": verts[:, 0],
                    "dv": verts[:, 1],
                    "ml": verts[:, 2],
                }
            )
        )
        faces_rows.append(
            pd.DataFrame(
                {
                    "region": key,
                    "i": faces[:, 0],
                    "j": faces[:, 1],
                    "k": faces[:, 2],
                }
            )
        )
        edge_verts = verts[edges]
        edges_rows.append(
            pd.DataFrame(
                {
                    "region": key,
                    "ap1": edge_verts[:, 0, 0],
                    "dv1": edge_verts[:, 0, 1],
                    "ml1": edge_verts[:, 0, 2],
                    "ap2": edge_verts[:, 1, 0],
                    "dv2": edge_verts[:, 1, 1],
                    "ml2": edge_verts[:, 1, 2],
                }
            )
        )
        regions.append(
            {
                "key": key,
                "label": label,
                "group": group,
                "level": level,
                "color": color,
                "n_annotation_ids": len(ids),
                "n_vertices": len(verts),
                "n_faces": len(faces),
                "n_wire_edges": len(edges),
            }
        )
        log(f"{key}: {len(verts):,} vertices, {len(faces):,} faces")
        return True

    root_ids = set(np.unique(annotation)) - {0}
    add_region(
        "root",
        "Whole brain",
        "Whole brain",
        "root",
        "#c8c8cf",
        root_ids,
        step=8,
        max_edges=MAX_ROOT_EDGES,
    )

    for domain_label, spec in TAXONOMY.items():
        valid_subdomains: list[tuple[str, str, str, set[int]]] = []
        for acronym, (source_name, color) in spec["subdomains"].items():
            valid_subdomains.append(
                (acronym, source_name, color, set(map(int, hierarchy[source_name])))
            )
        source_name = spec["source"]
        if source_name is None:
            domain_ids = set().union(*(item[3] for item in valid_subdomains))
        elif source_name in hierarchy:
            domain_ids = set(map(int, hierarchy[source_name]))
        domain_key = f"domain__{safe_key(domain_label)}"
        domain_ok = add_region(
            domain_key,
            domain_label,
            domain_label,
            "domain",
            spec["color"],
            domain_ids,
            step=6,
            max_edges=MAX_DOMAIN_EDGES,
        )
        if not domain_ok:
            continue
        for acronym, source_name, color, ids in valid_subdomains:
            add_region(
                f"region__{acronym}",
                f"{acronym} · {source_name}",
                domain_label,
                "subdomain",
                color,
                ids,
                step=5,
                max_edges=MAX_SUBDOMAIN_EDGES,
            )

    manifest = {
        "analysis_date": ANALYSIS_DATE,
        "annotation_path": str(ANNOTATION_PATH),
        "hierarchy_path": str(HIERARCHY_PATH),
        "voxel_um": VOXEL_UM,
        "viewer_coordinates": {"x": "AP", "y": "ML", "z": "DV"},
        "palette": "Curated anatomical hue families adapted from supplied OKLCH palette",
        "domain_semantics": (
            "Allen hierarchy composites; domains can overlap and are not a disjoint partition"
        ),
        "regions": regions,
    }
    outputs = {
        VERTICES_OUT: pd.concat(vertices_rows, ignore_index=True),
        FACES_OUT: pd.concat(faces_rows, ignore_index=True),
        EDGES_OUT: pd.concat(edges_rows, ignore_index=True),
    }
    staged: dict[Path, Path] = {}
    for final_path, frame in outputs.items():
        temp_path = final_path.with_name(
            f".tmp-{os.getpid()}-{final_path.name}"
        )
        frame.to_csv(temp_path, index=False, compression="gzip")
        staged[final_path] = temp_path
    manifest_temp = MANIFEST_OUT.with_name(
        f".tmp-{os.getpid()}-{MANIFEST_OUT.name}"
    )
    manifest_temp.write_text(json.dumps(manifest, indent=2))
    for final_path, temp_path in staged.items():
        if final_path.exists():
            raise FileExistsError(f"Refusing to overwrite: {final_path}")
        temp_path.rename(final_path)
    if MANIFEST_OUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {MANIFEST_OUT}")
    manifest_temp.rename(MANIFEST_OUT)
    log(f"Wrote {len(regions)} region meshes to {OUT_DIR}")


if __name__ == "__main__":
    main()
