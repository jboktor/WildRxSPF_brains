#!/usr/bin/env python
"""Rename legacy tissue__{slug(supertype_name)} dirs → tissue__{supertype_label}.

Uses cell_meta_keep_cell.csv for the name↔label map. Safe no-op if already migrated.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import pandas as pd

WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
DSPIN = WKDIR / "data/interim/dspin"
META = DSPIN / "cell_meta_keep_cell.csv"


def slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")


def main():
    md = pd.read_csv(META, usecols=["tissue", "supertype_label", "supertype_name"])
    map_df = (
        md.drop_duplicates(["tissue", "supertype_name"])
        .assign(old_stem=lambda d: d["tissue"] + "__" + d["supertype_name"].map(slug))
        .assign(new_stem=lambda d: d["tissue"] + "__" + d["supertype_label"].astype(str))
    )
    n_move = n_skip = n_miss = n_conflict = 0
    for _, row in map_df.iterrows():
        old = DSPIN / row["old_stem"]
        new = DSPIN / row["new_stem"]
        if old.resolve() == new.resolve():
            n_skip += 1
            continue
        if not old.is_dir():
            # already at label path, or never exported
            if new.is_dir():
                n_skip += 1
            else:
                n_miss += 1
            continue
        if new.exists():
            print(f"CONFLICT keep both: {old.name} and {new.name}")
            n_conflict += 1
            continue
        print(f"mv {old.name} → {new.name}")
        old.rename(new)
        n_move += 1
    print(f"done move={n_move} skip={n_skip} miss={n_miss} conflict={n_conflict}")


if __name__ == "__main__":
    main()
