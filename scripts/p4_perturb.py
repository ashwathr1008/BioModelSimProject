"""Compute one independent RegVelo perturbation and its cell-level fate effects."""
import argparse
import hashlib
import json
import sys
import time
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from filelock import FileLock
from p4_inference import run_perturbation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from p3_analysis.fate_inference import ENDPOINTS, infer_fate_probabilities, load_terminal_anchors


def run_name(condition, repeat):
    if condition == "tutorial":
        if repeat != 0:
            raise ValueError("The archived tutorial fit has repeat 0.")
        return "tutorial"
    return f"{condition}_r{repeat}"


def run(condition, repeat, regulator, batch_size):
    started = time.monotonic()
    name = run_name(condition, repeat)
    tutorial = condition == "tutorial"
    folder = ROOT / "results/p2/tutorial" if tutorial else ROOT / "results/p2/X" / name
    p3_dir = ROOT / "results/p3" / ("tutorial" if tutorial else "fate")
    baseline_file = p3_dir / f"{name}_fate.csv"
    anchor_file = ROOT / "results/p3/terminal_anchors/terminal_anchors.json"
    baseline = pd.read_csv(baseline_file, index_col="cell_id")
    if tutorial:
        baseline = baseline.assign(condition=condition, repeat=repeat)
    anchors = load_terminal_anchors(anchor_file)
    adata = ad.read_h5ad(folder / "fitted.h5ad")
    record = json.loads((folder / "run.json").read_text())
    if (record["condition"], record["repeat"]) != (condition, repeat):
        raise ValueError(f"Incorrect checkpoint metadata for {name}")
    if not (baseline.index.is_unique and adata.obs_names.is_unique and
            set(baseline.index) == set(adata.obs_names) and
            (baseline["condition"] == condition).all() and (baseline["repeat"] == repeat).all()):
        raise ValueError(f"Baseline identity mismatch for {name}")
    if not adata.var_names.is_unique or regulator not in adata.var_names:
        raise ValueError(f"Cannot identify {regulator} uniquely")
    baseline = baseline.loc[adata.obs_names]
    if not np.array_equal(baseline["cell_type"], adata.obs["leiden_R_celltype"]):
        raise ValueError("Cell types are misaligned")
    draws = record["settings"]["posterior_samples"]
    seed = 4100 + repeat
    print(f"{name} {regulator}: computing perturbed velocity, {draws} samples", flush=True)
    pert, model, checks = run_perturbation(folder, adata, regulator, draws, batch_size, seed)
    if not pert.obs_names.equals(adata.obs_names) or not pert.var_names.equals(adata.var_names):
        raise ValueError("Unexpected change to fitted cell or gene ordering")
    weights = model.module.v_encoder.fc1.weight.detach().cpu()
    original = torch.load(folder / "model/model.pt", map_location="cpu", weights_only=False)["model_state_dict"]["v_encoder.fc1.weight"]
    position = adata.var_names.get_loc(regulator)
    edges_removed = torch.count_nonzero(original[:, position]).item()
    original[:, position] = 0
    if not edges_removed or not torch.equal(weights, original):
        raise ValueError(f"Unexpected GRN edit for {regulator}")
    difference = np.asarray(pert.layers["velocity"]) - np.asarray(adata.layers["velocity"])
    if not np.isfinite(difference).all() or not np.any(difference):
        raise ValueError("Missing or unchanged perturbed velocity")
    print(f"{name} {regulator}: running Person 3 fate inference", flush=True)
    fate = infer_fate_probabilities(pert, anchors).loc[baseline.index]
    tables = []
    for endpoint in ENDPOINTS:
        table = baseline[["condition", "repeat", "cell_type"]].copy()
        table["regulator"] = regulator
        table["endpoint"] = endpoint
        table["p_unperturbed"] = baseline[endpoint]
        table["p_perturbed"] = fate[endpoint]
        table["delta_p"] = table["p_perturbed"] - table["p_unperturbed"]
        tables.append(table)
    out = ROOT / "results/p4" / name
    out.mkdir(parents=True, exist_ok=True)
    pending = out / f"{regulator}.pending.csv.gz"
    pd.concat(tables).to_csv(pending, index_label="cell_id")
    pending.replace(out / f"{regulator}.csv.gz")
    info = dict(condition=condition, repeat=repeat, regulator=regulator, seed=seed,
                posterior_samples=draws, batch_size=batch_size, cutoff=0, effects=0,
                outgoing_weights_removed=edges_removed,
                baseline_sha256=hashlib.sha256(baseline_file.read_bytes()).hexdigest(),
                anchors_sha256=hashlib.sha256(anchor_file.read_bytes()).hexdigest(),
                n_cells=len(fate), device=str(model.device), seconds=time.monotonic()-started, **checks)
    pending_info = out / f"{regulator}.pending.json"
    pending_info.write_text(json.dumps(info, indent=2) + "\n")
    pending_info.replace(out / f"{regulator}.json")
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("condition", choices=["full", "targeted", "random", "tutorial"])
    parser.add_argument("repeat", type=int, choices=range(3))
    parser.add_argument("regulator", choices=["MSC", "MYOG"])
    parser.add_argument("--batch-size", type=int, default=256)
    args = parser.parse_args()
    out = ROOT / "results/p4" / run_name(args.condition, args.repeat)
    key = hashlib.sha256(str(out / args.regulator).encode()).hexdigest()
    with FileLock(Path(tempfile.gettempdir()) / f"regvelo-p4-{key}.lock"):
        metadata = out / f"{args.regulator}.json"
        if metadata.exists() and (out / f"{args.regulator}.csv.gz").is_file():
            saved = json.loads(metadata.read_text())
            if saved.get("batch_size") == args.batch_size and saved.get("reference_exact") is True:
                print(f"Already complete: {out.name} {args.regulator}")
                sys.exit(0)
        run(args.condition, args.repeat, args.regulator, args.batch_size)
