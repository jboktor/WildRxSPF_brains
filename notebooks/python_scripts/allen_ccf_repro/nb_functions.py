"""Load the original notebook's registration helpers without executing the notebook.

The rebuild deliberately reuses the geometry helpers that produced the canonical
results (`image_from_df`, `get_cell_metadata_for_slice_index`, `modify_dapi`,
`modify_nissl`, `crop_and_pad_image`, `write_pts_file`, `read_outputpoints_file`)
so that fixed-image space is identical to the historical pipeline.
"""

from __future__ import annotations

import json
import pickle
import threading
from types import ModuleType

from . import config

_LOCK = threading.Lock()
_NAMESPACE: ModuleType | None = None
_HELPER_CELLS = (0, 2, 4, 35)


def notebook_namespace() -> ModuleType:
    """Return a module-like namespace holding the notebook helper functions."""
    global _NAMESPACE
    with _LOCK:
        if _NAMESPACE is not None:
            return _NAMESPACE
        module = ModuleType("allen_ccf_notebook_helpers")
        notebook = json.loads(config.NOTEBOOK.read_text())
        for index in _HELPER_CELLS:
            source = "".join(notebook["cells"][index]["source"])
            source = "\n".join(
                line for line in source.splitlines() if not line.lstrip().startswith("%")
            )
            exec(compile(source, f"{config.NOTEBOOK}:cell_{index}", "exec"), module.__dict__)
        with open(config.ALLEN_NAME_MAP, "rb") as handle:
            module.allen_name_to_annots = pickle.load(handle)
        _NAMESPACE = module
        return module


def allen_name_to_annots() -> dict:
    return notebook_namespace().allen_name_to_annots


def __getattr__(name: str):
    namespace = notebook_namespace()
    if hasattr(namespace, name):
        return getattr(namespace, name)
    raise AttributeError(name)
