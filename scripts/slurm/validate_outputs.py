"""Check completed fit artifacts before dependent work; no score thresholds."""

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from p2_fit import load_config, repo_path, run_identity, sha256


def imported_full0_authorization(record):
    """Authorize one laptop fit only through a separately pinned transfer receipt."""
    if os.environ.get("P2_REUSE_LOCAL_FULL0", "0") != "1":
        return None
    receipt_path = os.environ.get("P2_IMPORT_COMPLETED")
    expected_hash = os.environ.get("P2_IMPORT_COMPLETED_SHA256")
    if not receipt_path or not expected_hash:
        raise ValueError("Local-fit reuse requires a receipt path and trusted receipt SHA-256")
    body = Path(receipt_path).read_bytes()
    digest = hashlib.sha256(body).hexdigest()
    if digest != expected_hash.lower():
        raise ValueError("Imported-fit receipt does not match the trusted SHA-256")
    receipt = json.loads(body)
    if receipt.get("schema_version") != 1 or receipt.get("run") != "X/full_r0":
        raise ValueError("Import receipt must authorize only X/full_r0")
    gpu = record.get("gpu", "")
    if receipt.get("source_gpu") != gpu or "RTX 4090 Laptop GPU" not in gpu:
        raise ValueError("Imported full-data fit must retain its recorded RTX 4090 Laptop GPU provenance")
    required = {f"results/p2/X/full_r0/{name}" for name in
                ("run.json", "training.csv", "fitted.h5ad", "model/model.pt")}
    files = receipt.get("files", [])
    if len(files) != len(required) or {entry["path"] for entry in files} != required:
        raise ValueError("Import receipt must hash exactly the four full_r0 artifacts")
    for entry in files:
        path = repo_path(entry["path"])
        if not path.is_file() or path.is_symlink() or sha256(path) != entry["sha256"].lower():
            raise ValueError(f"Imported artifact hash mismatch: {entry['path']}")
    return {"receipt_sha256": digest, "source_gpu": gpu, "run": "X/full_r0"}


def planned_runs(cfg):
    return [(pca, condition, repeat)
            for repeat in range(cfg["sampling"]["n_repeats"])
            for condition in cfg["sampling"]["conditions"]
            for pca in ("X", "shipped_pca")]


def validate_run(cfg, pca, condition, repeat, expected_gpu=""):
    folder = repo_path("results/p2") / pca / f"{condition}_r{repeat}"
    record = json.loads((folder / "run.json").read_text(encoding="utf-8"))
    if record.get("status") != "complete":
        raise ValueError(f"Fit is not complete: {folder}")
    if record.get("identity") != run_identity(cfg, condition, repeat, pca):
        raise ValueError(f"Fit inputs/settings differ: {folder}")
    imported = imported_full0_authorization(record) if (pca, condition, repeat) == ("X", "full", 0) else None
    if not imported and expected_gpu and expected_gpu.lower() not in record.get("gpu", "").lower():
        raise ValueError(f"Unexpected GPU at {folder}: {record.get('gpu')}")
    if not math.isfinite(float(record["stage_spearman"])):
        raise ValueError(f"Nonfinite stage correlation: {folder}")
    for name in ("model/model.pt", "fitted.h5ad", "training.csv"):
        artifact = folder / name
        if not artifact.is_file() or artifact.stat().st_size == 0:
            raise ValueError(f"Missing or empty artifact: {artifact}")
    history = pd.read_csv(folder / "training.csv")
    if history.empty:
        raise ValueError(f"Empty training history: {folder}")
    for name in ("elbo_train", "elbo_validation"):
        values = pd.to_numeric(history[name], errors="raise").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite {name}: {folder}")
    fitted = ad.read_h5ad(folder / "fitted.h5ad", backed="r")
    try:
        if fitted.n_obs != record["n_cells"] or fitted.n_vars != record["n_genes"]:
            raise ValueError(f"Fitted dimensions differ from run record: {folder}")
        if not fitted.obs_names.is_unique:
            raise ValueError(f"Duplicate fitted cell IDs: {folder}")
        times = fitted.obs["latent_time"].to_numpy(dtype=float)
        stages = fitted.obs[cfg["data"]["stage_key"]].astype(str).str.removeprefix("Pcw").astype(float)
        if not np.isfinite(times).all() or not np.isfinite(stages).all():
            raise ValueError(f"Nonfinite fitted cell times/stages: {folder}")
    finally:
        fitted.file.close()
    return {"run": f"{pca}/{condition}_r{repeat}", "gpu": record.get("gpu"),
            "epochs": record["epochs"], "stage_spearman": record["stage_spearman"],
            "imported_local": imported}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope", choices=("pilots", "all", "local-full0"), default="pilots")
    parser.add_argument("--expected-gpu", default="")
    args = parser.parse_args()
    cfg = load_config(args.config)
    runs = [(pca, "full", 0) for pca in ("X", "shipped_pca")]
    if args.scope == "all":
        runs = planned_runs(cfg)
    elif args.scope == "local-full0":
        if os.environ.get("P2_REUSE_LOCAL_FULL0", "0") != "1":
            raise ValueError("Local-fit validation requires P2_REUSE_LOCAL_FULL0=1")
        runs = [("X", "full", 0)]
    checked = [validate_run(cfg, *run, args.expected_gpu) for run in runs]
    print(json.dumps({"validated": checked}, indent=2), flush=True)


if __name__ == "__main__":
    main()
