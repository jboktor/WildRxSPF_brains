"""Paths, run settings, and the frozen region registry for the CCF rebuild."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

BASE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_ROOT = PACKAGE_ROOT / "configs"
REGISTRY_FILE = CONFIG_ROOT / "region_registry.yaml"
MASK_OVERRIDE_ROOT = CONFIG_ROOT / "mask_overrides"

NOTEBOOK = BASE / "notebooks/0X_Allen_CCF_Registration_optimized.ipynb"
REFERENCE_ROOT = BASE / "data/input/allen_registration_ref"
SOURCE_WORKBOOK = BASE / "data/interim/registration/Allen_CCF_optimization/slice_positions_25um_final.xlsx"
BASELINE_CCF = BASE / "data/interim/registration/Allen_CCF_optimization/ccf3d.csv"

RUN_ROOT = BASE / "data/interim/registration/Allen_CCF_rebuild"
WORK_ROOT = RUN_ROOT / "work"
FIGURE_ROOT = BASE / "figures/Allen_CCF_alignment_optimized"
REBUILD_FIGURES = FIGURE_ROOT / "final_rebuild"
PROMOTED_FIGURES = FIGURE_ROOT / "final"
CANONICAL_FIGURES = BASE / "figures/Allen_CCF_alignment"

ANNOTATION_25 = REFERENCE_ROOT / "annotation_25.nrrd"
ANNOTATION_10 = REFERENCE_ROOT / "annotation_10.nrrd"
NISSL_25 = REFERENCE_ROOT / "ara_nissl_25.nrrd"
BORDERS_25 = REFERENCE_ROOT / "annotations_25_border_only.tif"
ALLEN_NAME_MAP = REFERENCE_ROOT / "allen_name_to_annots.pkl"

CCF_PIXEL_SIZE = 25
ATLAS_MIDLINE = 228
MAX_MEDIAL_EXTRA = 40
MIN_MEDIAL_EXTRA = 8
MEDIAL_EXTRA_MARGIN = 6
PAD_WIDTH = 20

SAMPLES = [
    "roi1_run1",
    "roi1_run2",
    "roi2_run1",
    "roi2_run2",
    "roi2_run4",
    "roi3_run1",
    "roi3_run2",
    "roi3_run3",
    "roi3_run4",
    "roi4_run1",
    "roi4_run2",
    "roi4_run3",
    "roi4_run4",
]


@dataclass
class RunConfig:
    """Everything that changes the numerical result of a pipeline run."""

    samples: list = field(default_factory=lambda: list(SAMPLES))
    ap_coarse_offsets: tuple = (-6, -3, 0, 3, 6)
    ap_fine_radius: int = 2
    medial_extras: tuple = (0, 4, 8, 14, 22)
    ap_gradients: tuple = (-0.05, -0.03, 0.0, 0.03, 0.05)
    ap_gradient_min_gain: float = 0.4
    landmark_weight: float = 0.03
    landmark_weight_ladder: tuple = (0.03, 0.015, 0.01, 0.05)
    # A bending penalty was measured to stiffen the spline enough to leave whole
    # structures behind: on roi2_run1 caudoputamen containment fell from 80% to
    # 4% as the weight went from 0 to 0.35, while the Jacobian stayed positive
    # throughout. Folding is therefore controlled by the Jacobian gate instead.
    bending_energy_weight: float = 0.0
    # Applied only to a fit that folds, as (bending weight, control grid scale).
    # Widening the grid removes degrees of freedom outright, which stops folding
    # more reliably than penalising curvature does.
    stiffening_ladder: tuple = ((0.15, 1.5), (0.25, 2.0), (0.35, 2.5))
    spline_grid_voxels: int = 13
    affine_iterations: int = 500
    spline_iterations: int = 2000
    search_iterations: int = 250
    tissue_margin_px: int = 14
    exclusion_halo_px: int = 6
    mask_sensitivity_px: tuple = (-3, 3)
    jobs: int = 4
    seed: int = 0

    def to_dict(self) -> dict:
        return {key: list(value) if isinstance(value, tuple) else value for key, value in self.__dict__.items()}


def load_registry(path: Path = REGISTRY_FILE) -> dict:
    """Read the frozen region registry."""
    with open(path) as handle:
        registry = yaml.safe_load(handle)
    for name, entry in registry["regions"].items():
        entry.setdefault("role", "soft")
        entry.setdefault("cell_classes", [])
        entry.setdefault("allen_ids", [])
        entry.setdefault("allen_names", [])
        entry.setdefault("min_cells", 40)
        entry.setdefault("min_atlas_pixels", 25)
        entry.setdefault("landmark", "none")
        entry["name"] = name
    return registry


def resolve_region_ids(entry: dict, allen_name_to_annots: dict) -> list:
    """Expand a registry entry into a sorted list of Allen annotation ids."""
    ids = set(int(value) for value in entry.get("allen_ids", []))
    for name in entry.get("allen_names", []):
        if name not in allen_name_to_annots:
            raise KeyError(f"Allen structure name absent from reference map: {name}")
        ids.update(int(value) for value in allen_name_to_annots[name])
    scope = entry.get("allen_name_scope")
    pattern = entry.get("allen_name_regex")
    if pattern:
        import re

        compiled = re.compile(pattern)
        scope_ids = None
        if scope:
            scope_ids = set(int(value) for value in allen_name_to_annots[scope])
        for name, values in allen_name_to_annots.items():
            if not compiled.search(str(name)):
                continue
            selected = set(int(value) for value in values)
            if scope_ids is not None:
                selected &= scope_ids
            ids.update(selected)
    return sorted(ids)


def ensure_directories() -> None:
    for path in (RUN_ROOT, WORK_ROOT, REBUILD_FIGURES):
        path.mkdir(parents=True, exist_ok=True)


def sample_work_dir(sample_id: str, stage: str) -> Path:
    path = WORK_ROOT / sample_id / stage
    path.mkdir(parents=True, exist_ok=True)
    return path


def env_python() -> str:
    return os.environ.get(
        "ALLEN_CCF_PYTHON",
        "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python",
    )


def as_json_safe(value: Any) -> Any:
    """Convert numpy scalars and paths so json.dump never fails on a manifest."""
    import numpy as np

    if isinstance(value, dict):
        return {str(key): as_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [as_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.ndarray):
        return as_json_safe(value.tolist())
    return value
