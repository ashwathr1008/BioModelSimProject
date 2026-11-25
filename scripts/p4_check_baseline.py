"""Check the full_r0 unperturbed probabilities against Person 3's published table."""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from p3_analysis.fate_inference import ENDPOINTS, infer_fate_probabilities, load_terminal_anchors


def main():
    adata = ad.read_h5ad(ROOT / "results/p2/X/full_r0/fitted.h5ad")
    anchors = load_terminal_anchors(ROOT / "results/p3/terminal_anchors/terminal_anchors.json")
    result = infer_fate_probabilities(adata, anchors)
    reference = pd.read_csv(ROOT / "results/p3/fate/full_r0_fate.csv", index_col="cell_id")
    reference = reference.loc[result.index, list(ENDPOINTS)]
    differences = np.abs(result.to_numpy() - reference.to_numpy())
    record = {
        "model": "full_r0",
        "cells": len(result),
        "anchors": sum(map(len, anchors.values())),
        "max_abs_difference": float(differences.max()),
        "mean_abs_difference": float(differences.mean()),
        "per_endpoint_max": dict(zip(ENDPOINTS, differences.max(axis=0).tolist())),
    }
    print(json.dumps(record, indent=2))
    out = ROOT / "results/p4/qc"
    out.mkdir(exist_ok=True, parents=True)
    (out / "baseline_reproduction.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
