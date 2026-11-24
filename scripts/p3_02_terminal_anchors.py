
from pathlib import Path
from multiprocessing import freeze_support

import anndata as ad
import cellrank as cr
import pandas as pd
import json

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = PROJECT_ROOT / "results/p2/X/full_r0/fitted.h5ad"


def main():
    # Load the fitted full-data model output.
    adata = ad.read_h5ad(DATA_PATH)

    # Construct the velocity kernel.
    print("\n===== Velocity Kernel =====")

    vk = cr.kernels.VelocityKernel(
        adata,
        xkey="Ms",
        vkey="velocity",
    )

    vk.compute_transition_matrix(n_jobs=1)

    # Construct the connectivity kernel.
    print("\n===== Connectivity Kernel =====")

    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    # Combine velocity and connectivity kernels.
    print("\n===== Combined Kernel =====")

    kernel = 0.8 * vk + 0.2 * ck
    kernel.compute_transition_matrix()

    print("Transition matrix shape:")
    print(kernel.transition_matrix.shape)

    # Initialize the GPCCA estimator.
    estimator = cr.estimators.GPCCA(kernel)

    print("\n===== GPCCA Macrostates =====")

    estimator.compute_macrostates(
        n_states=9,
        cluster_key="leiden_R_celltype",
    )

    print("\n===== Identified Macrostates =====")
    print(estimator.macrostates.value_counts())
    # Define the mapping from macrostates to terminal lineages.
    terminal_mapping = {
        "PAX7+MyoProg1": [
            "PAX7+MyoProg1_1",
            "PAX7+MyoProg1_2",
        ],
        "MYH3+MyoC": ["MYH3+MyoC"],
        "MYL3+MyoC": ["MYL3+MyoC"],
    }

    output_dir = PROJECT_ROOT / "results/p3/terminal_anchors"
    output_dir.mkdir(parents=True, exist_ok=True)

    anchors = {}
    records = []

    macrostates = estimator.macrostates

    print("\n===== Terminal Anchor Selection =====")

    for endpoint, states in terminal_mapping.items():
        mask = macrostates.isin(states)

        candidate_ids = macrostates.index[mask]

        # Keep only cells with the matching biological annotation.
        cell_types = adata.obs.loc[
            candidate_ids, "leiden_R_celltype"
        ]

        cell_ids = cell_types.index[
            cell_types == endpoint
        ].tolist()

        if not cell_ids:
            raise ValueError(f"No valid anchors for {endpoint}")

        anchors[endpoint] = cell_ids

        print(
            f"{endpoint}: {len(cell_ids)} valid anchors "
            f"out of {len(candidate_ids)} candidates"
        )

        for cell_id in cell_ids:
            records.append({
                "cell_id": cell_id,
                "endpoint": endpoint,
                "macrostate": str(macrostates.loc[cell_id]),
            })

    # Validate anchor availability across all experimental runs.
    print("\n===== Cross-condition Validation =====")

    all_anchor_ids = {
        cell_id
        for ids in anchors.values()
        for cell_id in ids
    }

    for condition in ["full", "targeted", "random"]:
        for repeat in range(3):
            run = f"{condition}_r{repeat}"

            path = (
                PROJECT_ROOT
                / "results/p2/X"
                / run
                / "fitted.h5ad"
            )

            subset = ad.read_h5ad(path, backed="r")

            try:
                missing = all_anchor_ids - set(subset.obs_names)

                print(f"{run}: {len(missing)} missing anchors")

                if missing:
                    raise ValueError(
                        f"{run} is missing terminal anchors"
                    )
            finally:
                subset.file.close()

    # Save the validated anchor manifest.
    with open(
        output_dir / "terminal_anchors.json",
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(anchors, f, indent=2)

    pd.DataFrame(records).to_csv(
        output_dir / "terminal_anchors.csv",
        index=False,
    )

    print("\n===== Anchor Selection Complete =====")
    print("Total anchors:", len(all_anchor_ids))
    print("Saved to:", output_dir)

if __name__ == "__main__":
    freeze_support()
    main()
