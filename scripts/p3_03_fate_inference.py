
from pathlib import Path
from multiprocessing import freeze_support
import sys

import anndata as ad


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from p3_analysis.fate_inference import (
    load_terminal_anchors,
    infer_fate_probabilities,
)


def main():

    condition = "full"
    repeat = 0

    data_path = (
        PROJECT_ROOT
        / "results/p2/X"
        / f"{condition}_r{repeat}"
        / "fitted.h5ad"
    )

    anchor_path = (
        PROJECT_ROOT
        / "results/p3/terminal_anchors"
        / "terminal_anchors.json"
    )

    output_dir = PROJECT_ROOT / "results/p3/fate"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the fitted dataset and fixed terminal anchors.
    adata = ad.read_h5ad(data_path)
    anchors = load_terminal_anchors(anchor_path)

    print("\n===== Dataset =====")
    print(adata.shape)

    print("\n===== Terminal Anchors =====")
    for endpoint, cells in anchors.items():
        print(endpoint, len(cells))

    # Compute fate probabilities.
    print("\n===== Fate Inference =====")

    probabilities = infer_fate_probabilities(
        adata,
        anchors,
        n_jobs=1,
    )

    print("\n===== Fate Probability Shape =====")
    print(probabilities.shape)

    print("\n===== First Five Cells =====")
    print(probabilities.head())

    print("\n===== Probability Summary =====")
    print(probabilities.describe())

    print("\n===== Row Sum Validation =====")
    print(probabilities.sum(axis=1).describe())

    # Add cell-level metadata.
    results = probabilities.copy()

    results["cell_type"] = adata.obs[
        "leiden_R_celltype"
    ]

    results["condition"] = condition
    results["repeat"] = repeat

    results["latent_time"] = adata.obs["latent_time"]
    results["adj_stage_num"] = adata.obs["adj_stage_num"]

    output_path = (
        output_dir
        / f"{condition}_r{repeat}_fate.csv"
    )

    results.to_csv(output_path)

    print("\n===== Fate Inference Complete =====")
    print("Saved to:", output_path)


if __name__ == "__main__":
    freeze_support()
    main()
