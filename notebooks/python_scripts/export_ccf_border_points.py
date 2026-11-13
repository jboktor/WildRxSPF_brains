#!/usr/bin/env python3
"""Export a subsampled Allen CCF border point cloud for interactive Plotly overlays."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
ANALYSIS_DATE = "2026-11-12"
BORDER_PATH = (
    WKDIR / "data/input/allen_registration_ref/annotations_25_border_only.tif"
)
OUT_PATH = (
    WKDIR
    / "data/interim/2DG_CCF_Integration"
    / ANALYSIS_DATE
    / f"allen_ccf_border_points_subsampled_{ANALYSIS_DATE}.csv.gz"
)
CCF_VOXEL_UM = 25.0
TARGET_POINTS = 18000
STRIDE_AP = 4
RNG = np.random.default_rng(42)


def main() -> None:
    if OUT_PATH.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUT_PATH}")
    image = Image.open(BORDER_PATH)
    n_frames = image.n_frames
    points = []
    for ap in range(0, n_frames, STRIDE_AP):
        image.seek(ap)
        frame = np.asarray(image)
        ys, xs = np.where(frame > 0)
        if len(ys) == 0:
            continue
        # Downsample within plane.
        keep = RNG.choice(len(ys), size=min(80, len(ys)), replace=False)
        for y, x in zip(ys[keep], xs[keep]):
            points.append(
                (
                    ap * CCF_VOXEL_UM,
                    float(y) * CCF_VOXEL_UM,
                    float(x) * CCF_VOXEL_UM,
                )
            )
    frame = pd.DataFrame(points, columns=["ccfx", "ccfy", "ccfz"])
    if len(frame) > TARGET_POINTS:
        frame = frame.sample(TARGET_POINTS, random_state=42)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT_PATH, index=False, compression="gzip")
    print(f"[ccf-border] wrote {len(frame)} points -> {OUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
