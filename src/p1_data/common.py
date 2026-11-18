# Shared helpers used by every Person 1 file.
# Loads the experiment config, resolves repo paths, and loads/verifies the source dataset.

# %%
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml


# Find the repo root (the folder that contains config/experiment.yaml).
def find_repo_root() -> Path:
    try:
        start = Path(__file__).resolve().parent
    except NameError:
        start = Path.cwd()
    for p in [start, *start.parents]:
        if (p / "config" / "experiment.yaml").exists():
            return p
    raise FileNotFoundError("could not find config/experiment.yaml above " + str(start))


REPO_ROOT = find_repo_root()


# %%
# Read the shared experiment settings from config/experiment.yaml.
def load_config(path: str | Path | None = None) -> dict:
    path = path or REPO_ROOT / "config" / "experiment.yaml"
    with open(path) as fh:
        return yaml.safe_load(fh)


# Turn a repo-relative path from the config into an absolute path.
def repo_path(rel: str | Path) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else REPO_ROOT / p


# Compute a file's SHA-256 hash (used to confirm we have the right dataset).
def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# Load the shipped human-limb dataset, downloading it first if it is missing.
def load_source(cfg: dict, check_hash: bool = True):
    import anndata as ad

    path = repo_path(cfg["data"]["source_file"])
    if not path.exists():
        import urllib.request

        path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading source data to {path} ...")
        urllib.request.urlretrieve(cfg["data"]["source_url"], path)
    if check_hash:
        got = sha256(path)
        want = cfg["data"].get("expected_sha256")
        if want and got != want:
            raise ValueError(f"sha256 mismatch for {path}: got {got}, expected {want}")
    return ad.read_h5ad(path)


# Convert a sparse or dense matrix to a dense numpy array.
def dense(M) -> np.ndarray:
    return M.toarray() if sp.issparse(M) else np.asarray(M)


# Build the standard run name, e.g. "targeted_r0".
def condition_id(condition: str, repeat: int) -> str:
    return f"{condition}_r{repeat}"
