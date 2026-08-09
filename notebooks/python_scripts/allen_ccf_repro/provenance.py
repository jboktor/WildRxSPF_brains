"""Run provenance: input hashes, git state, and the reproducibility manifest."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from . import config


def file_digest(path: Path, chunk: int = 1 << 20) -> dict:
    path = Path(path)
    if not path.exists():
        return {"path": str(path), "exists": False}
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
            size += len(block)
    return {"path": str(path), "exists": True, "sha256": digest.hexdigest(), "bytes": size}


def git_state() -> dict:
    def run(args):
        try:
            return subprocess.run(
                args, cwd=config.BASE, capture_output=True, text=True, check=True
            ).stdout.strip()
        except Exception as error:  # noqa: BLE001 - provenance must never break a run
            return f"unavailable: {error}"

    return {
        "commit": run(["git", "rev-parse", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty": bool(run(["git", "status", "--porcelain"])),
    }


def input_digests() -> dict:
    return {
        "notebook": file_digest(config.NOTEBOOK),
        "workbook": file_digest(config.SOURCE_WORKBOOK),
        "baseline_ccf": file_digest(config.BASELINE_CCF),
        "annotation_25": file_digest(config.ANNOTATION_25),
        "nissl_25": file_digest(config.NISSL_25),
        "borders_25": file_digest(config.BORDERS_25),
        "region_registry": file_digest(config.REGISTRY_FILE),
    }


def package_digests() -> dict:
    return {
        path.name: file_digest(path)
        for path in sorted(config.PACKAGE_ROOT.glob("*.py"))
    }


def mask_override_digests() -> dict:
    return {
        path.name: file_digest(path)
        for path in sorted(config.MASK_OVERRIDE_ROOT.glob("*.yaml"))
    }


def write_manifest(path: Path, run_config, registry: dict, samples: list, extra: dict = None) -> Path:
    payload = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "host": platform.node(),
        "command": "python -m allen_ccf_repro.cli " + " ".join(sys.argv[1:]),
        # A --skip-fits invocation rewrites this file from cached fits, so say so
        # rather than let the recorded command imply it produced them.
        "fits_computed_this_run": "--skip-fits" not in sys.argv[1:],
        "refit_command": "python -m allen_ccf_repro.cli --all --jobs 5 --promote",
        "git": git_state(),
        "run_config": run_config.to_dict(),
        "region_registry_version": registry.get("version"),
        "inputs": input_digests(),
        "pipeline_modules": package_digests(),
        "mask_overrides": mask_override_digests(),
        "samples": samples,
    }
    if extra:
        payload.update(extra)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(config.as_json_safe(payload), handle, indent=2, sort_keys=False)
    return path
