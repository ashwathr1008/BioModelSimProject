
from pathlib import Path
from multiprocessing import freeze_support
import sys

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from p3_analysis.fate_inference import (
    load_terminal_anchors,
    build_kernel,
    infer_fate_probabilities,
)

from p3_analysis.trajectory import (
    compute_latent_time_score,
    compute_cbc_scores,
)


def main():

    tutorial_dir = PROJECT_ROOT / "results/p2/tutorial"

    output_dir = PROJECT_ROOT / "results/p3/tutorial"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Locate the fitted tutorial dataset.
    paths = list(tutorial_dir.rglob("fitted.h5ad"))

    if len(paths) != 1:
        raise ValueError(
            f"Expected one tutorial dataset, found {len(paths)}"
        )

    adata = ad.read_h5ad(paths[0])

    print("\n===== Tutorial Dataset =====")
    print("Path:", paths[0])
    print(adata)

    # Load the same fixed anchors used in the main experiment.
    anchors = load_terminal_anchors(
        PROJECT_ROOT
        / "results/p3/terminal_anchors/terminal_anchors.json"
    )

    print("\n===== Anchor Validation =====")

    for endpoint, cell_ids in anchors.items():

        missing = set(cell_ids) - set(adata.obs_names)

        print(
            f"{endpoint}: "
            f"{len(cell_ids)} anchors, "
            f"{len(missing)} missing"
        )

        if missing:
            raise ValueError(
                f"Missing tutorial anchors: {endpoint}"
            )

    # Evaluate developmental timing.
    print("\n===== Latent Time =====")

    timing = compute_latent_time_score(adata)
    print(timing)

    pd.DataFrame([timing]).to_csv(
        output_dir / "tutorial_timing.csv",
        index=False,
    )

    # Construct the transition kernel.
    print("\n===== Building Kernel =====")

    kernel = build_kernel(adata, n_jobs=1)

    # Compute fate probabilities.
    print("\n===== Fate Inference =====")

    fate = infer_fate_probabilities(
        adata,
        anchors,
        kernel=kernel,
        n_jobs=1,
    )

    print(fate.describe())

    fate["cell_type"] = adata.obs["leiden_R_celltype"]
    fate["latent_time"] = adata.obs["latent_time"]
    fate["adj_stage_num"] = adata.obs["adj_stage_num"]

    fate.to_csv(output_dir / "tutorial_fate.csv")

    # Calculate transition-specific CBC.
    print("\n===== CBC Analysis =====")

    cbc = compute_cbc_scores(kernel, adata)

    print(cbc.to_string(index=False))

    cbc.to_csv(
        output_dir / "tutorial_cbc.csv",
        index=False,
    )

    print("\n===== Tutorial Validation Complete =====")
    print("Saved to:", output_dir)


if __name__ == "__main__":
    freeze_support()
    main()
