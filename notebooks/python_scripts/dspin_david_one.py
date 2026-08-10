#!/usr/bin/env python
"""Run DSPIN DAVID annotation for one full_manifest unit (SLURM array task).

Reads gene_programs_total_*_prior_0.csv from the unit's dspin_run/, cleans
mouse gene symbols, calls dspin.annotate.automatic_david_annotation, and
writes a tidy program_annotations_david.csv for the Analysis notebook.

Auth: DAVID_TOKEN env var, or first line of .secrets/david_token
(email registered at https://davidbioinformatics.nih.gov/).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import ssl
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Cluster Python often lacks system CA bundle; pin certifi before mygene/httpx/suds.
try:
    import certifi

    _ca = certifi.where()
    os.environ.setdefault("SSL_CERT_FILE", _ca)
    os.environ.setdefault("REQUESTS_CA_BUNDLE", _ca)
    os.environ.setdefault("CURL_CA_BUNDLE", _ca)
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=_ca)
except Exception:
    pass

import pandas as pd

WKDIR = Path("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains")
DSPIN_DIR = WKDIR / "data/interim/dspin"
SECRET = WKDIR / ".secrets" / "david_token"


def resolve_token() -> str:
    tok = os.environ.get("DAVID_TOKEN") or os.environ.get("DAVID_EMAIL")
    if tok and tok.strip():
        return tok.strip()
    if SECRET.exists():
        line = SECRET.read_text().strip().splitlines()[0].strip()
        if line and not line.startswith("#"):
            return line
    raise SystemExit(
        "Missing DAVID token. Set DAVID_TOKEN to your DAVID-registered email, "
        f"or write it to {SECRET}"
    )


def unit_dir(row) -> Path:
    for key in ("unit_path", "slice_path"):
        if key in row and pd.notna(row[key]):
            return Path(row[key])
    return DSPIN_DIR / f"{row['tissue']}__{row['supertype_label']}"


def find_program_csv(sdir: Path) -> Path:
    run = sdir / "dspin_run"
    cands = sorted(run.glob("gene_programs_total_*_prior_0.csv"))
    if not cands:
        raise FileNotFoundError(f"No gene_programs_total_*_prior_0.csv under {run}")
    # Prefer exact total_20 if present
    for c in cands:
        if "total_20" in c.name:
            return c
    return cands[0]


def clean_symbol(g: str) -> str:
    if g is None or (isinstance(g, float) and pd.isna(g)):
        return ""
    s = str(g).strip()
    if not s or s.lower() == "nan":
        return ""
    # Strip Ensembl / GRCm39 suffixes from Seurat gene names
    s = re.sub(r"-GRCm39(-\d+)?$", "", s)
    s = re.sub(r"^ENSMUSG\d+$", "", s)
    return s


def prepare_clean_csv(raw_csv: Path, out_csv: Path) -> Path:
    df = pd.read_csv(raw_csv)
    cleaned = df.apply(lambda col: col.map(clean_symbol))
    # Drop empty cells → NaN so DSPIN dropna works
    cleaned = cleaned.replace("", pd.NA)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_csv(out_csv, index=False)
    return out_csv


def tidy_from_david_results(results_dir: Path, program_csv: Path) -> pd.DataFrame:
    """Parse list_*.csv DAVID cluster dumps → top term per program."""
    prog_cols = list(pd.read_csv(program_csv, nrows=0).columns)
    rows = []
    for col in prog_cols:
        f = results_dir / f"list_{col}.csv"
        if not f.exists():
            rows.append(
                {
                    "program_index": col,
                    "program_id": f"P{col}",
                    "top_term": None,
                    "top_category": None,
                    "top_pvalue": None,
                    "enrichment_score_cluster1": None,
                    "n_terms_p05": 0,
                    "status": "missing_list",
                }
            )
            continue
        # Files mix cluster headers + TSV rows; parse loosely
        text = f.read_text(errors="replace").splitlines()
        cluster1_score = None
        terms = []
        i = 0
        while i < len(text):
            line = text[i]
            if line.startswith("Annotation Cluster"):
                m = re.search(r"EnrichmentScore:([0-9.eE+-]+)", line)
                score = float(m.group(1)) if m else None
                if cluster1_score is None and score is not None:
                    cluster1_score = score
                i += 1
                if i < len(text) and line.startswith("Annotation Cluster"):
                    # header row follows
                    if i < len(text) and text[i].startswith("Category"):
                        i += 1
                continue
            if line.startswith("Category\t"):
                i += 1
                continue
            parts = line.split("\t")
            if len(parts) >= 5 and parts[0] and not parts[0].startswith("Annotation"):
                try:
                    pval = float(parts[4])
                except ValueError:
                    i += 1
                    continue
                terms.append(
                    {
                        "category": parts[0],
                        "term": parts[1].split("~")[-1] if parts[1] else parts[1],
                        "pvalue": pval,
                    }
                )
            i += 1
        sig = [t for t in terms if t["pvalue"] <= 0.05]
        sig_sorted = sorted(sig or terms, key=lambda t: t["pvalue"])
        top = sig_sorted[0] if sig_sorted else None
        rows.append(
            {
                "program_index": col,
                "program_id": f"P{col}",
                "top_term": None if top is None else top["term"],
                "top_category": None if top is None else top["category"],
                "top_pvalue": None if top is None else top["pvalue"],
                "enrichment_score_cluster1": cluster1_score,
                "n_terms_p05": len(sig),
                "status": "ok" if top is not None else "no_terms",
            }
        )
    return pd.DataFrame(rows)


def run_one(row, token: str, force: bool = False, max_attempts: int = 4) -> dict:
    from dspin import annotate

    sdir = unit_dir(row)
    status_path = sdir / "david_status.json"
    tidy_path = sdir / "program_annotations_david.csv"

    if (
        not force
        and status_path.exists()
        and tidy_path.exists()
        and json.loads(status_path.read_text()).get("status") == "ok"
    ):
        print("SKIP already ok:", sdir.name, flush=True)
        return json.loads(status_path.read_text())

    raw = find_program_csv(sdir)
    clean = sdir / "dspin_run" / "gene_programs_clean_for_david.csv"
    prepare_clean_csv(raw, clean)

    t0 = time.time()
    status = {
        "tissue": row.get("tissue"),
        "supertype_label": row.get("supertype_label"),
        "supertype_name": row.get("supertype_name"),
        "unit_id": sdir.name,
        "program_csv": str(raw),
        "clean_csv": str(clean),
        "started": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "error": None,
    }
    print(f"DAVID annotate {sdir.name} from {clean.name}", flush=True)

    data_folder = str(clean.parent) + "/"
    file_name = clean.name
    last_err = None  # type: Exception | None
    for attempt in range(1, max_attempts + 1):
        try:
            annotate.automatic_david_annotation(
                data_folder=data_folder,
                file_name=file_name,
                token_name=token,
                species="mouse",
            )
            last_err = None
            break
        except Exception as e:
            last_err = e
            wait = min(120, 15 * attempt)
            print(
                f"attempt {attempt}/{max_attempts} failed: {type(e).__name__}: {e}; "
                f"sleep {wait}s",
                flush=True,
            )
            if attempt < max_attempts:
                time.sleep(wait)
    if last_err is not None:
        raise last_err

    results_dir = clean.parent / file_name.replace(".csv", "_david_results")
    tidy = tidy_from_david_results(results_dir, clean)
    tidy.insert(0, "tissue", row.get("tissue"))
    tidy.insert(1, "supertype_label", row.get("supertype_label"))
    tidy.insert(2, "supertype_name", row.get("supertype_name"))
    tidy.to_csv(tidy_path, index=False)

    status["status"] = "ok"
    status["n_programs"] = int(len(tidy))
    status["n_with_term"] = int(tidy["top_term"].notna().sum())
    status["elapsed_s"] = round(time.time() - t0, 1)
    status["finished"] = datetime.now(timezone.utc).isoformat()
    status["tidy"] = str(tidy_path)
    status_path.write_text(json.dumps(status, indent=2))
    print("DONE", status, flush=True)
    return status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--max-attempts", type=int, default=4)
    args = ap.parse_args()

    man = pd.read_csv(args.manifest)
    if args.index < 0 or args.index >= len(man):
        raise SystemExit(f"index {args.index} out of range 0..{len(man)-1}")
    row = man.iloc[args.index]
    token = resolve_token()
    sdir = unit_dir(row)
    try:
        run_one(row, token=token, force=args.force, max_attempts=args.max_attempts)
    except Exception as e:
        sdir.mkdir(parents=True, exist_ok=True)
        err = {
            "status": "error",
            "tissue": row.get("tissue"),
            "supertype_label": None
            if pd.isna(row.get("supertype_label"))
            else str(row.get("supertype_label")),
            "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
        }
        (sdir / "david_status.json").write_text(json.dumps(err, indent=2))
        print("ERROR", err["error"], flush=True)
        raise


if __name__ == "__main__":
    main()
