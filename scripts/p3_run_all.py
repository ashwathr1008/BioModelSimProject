
from pathlib import Path
from multiprocessing import freeze_support
import sys
import gc

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


CONDITIONS = ("full", "targeted", "random")
REPEATS = (0, 1, 2)


def main():

    input_dir = PROJECT_ROOT / "results/p2/X"
    output_dir = PROJECT_ROOT / "results/p3"

    fate_dir = output_dir / "fate"
    trajectory_dir = output_dir / "trajectory"
    summary_dir = output_dir / "summary"

    for directory in (
        fate_dir,
        trajectory_dir,
        summary_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    anchors = load_terminal_anchors(
        output_dir
        / "terminal_anchors"
        / "terminal_anchors.json"
    )

    all_cbc = []
    all_timing = []

    for repeat in REPEATS:
        for condition in CONDITIONS:

            run = f"{condition}_r{repeat}"

            print("\n" + "=" * 60)
            print(f"Processing: {run}")
            print("=" * 60)

            data_path = input_dir / run / "fitted.h5ad"

            adata = ad.read_h5ad(data_path)

            print("Dataset shape:", adata.shape)

            # Build the transition kernel once per run.
            print("\n===== Building Kernel =====")

            kernel = build_kernel(
                adata,
                n_jobs=1,
            )

            # Compute cell fate probabilities.
            print("\n===== Fate Inference =====")

            fate = infer_fate_probabilities(
                adata,
                anchors,
                n_jobs=1,
                kernel=kernel,
            )

            fate["cell_type"] = adata.obs[
                "leiden_R_celltype"
            ]

            fate["condition"] = condition
            fate["repeat"] = repeat

            fate["latent_time"] = adata.obs[
                "latent_time"
            ]

            fate["adj_stage_num"] = adata.obs[
                "adj_stage_num"
            ]

            fate.to_csv(
                fate_dir / f"{run}_fate.csv"
            )

            # Compute developmental timing.
            print("\n===== Latent Time =====")

            timing = compute_latent_time_score(adata)

            timing["condition"] = condition
            timing["repeat"] = repeat

            timing_df = pd.DataFrame([timing])

            timing_df.to_csv(
                trajectory_dir / f"{run}_timing.csv",
                index=False,
            )

            all_timing.append(timing_df)

            # Compute cross-boundary correctness.
            print("\n===== CBC Analysis =====")

            cbc = compute_cbc_scores(
                kernel,
                adata,
            )

            cbc["condition"] = condition
            cbc["repeat"] = repeat

            cbc.to_csv(
                trajectory_dir / f"{run}_cbc.csv",
                index=False,
            )

            all_cbc.append(cbc)

            print(cbc.to_string(index=False))

            print(f"\nCompleted: {run}")

            # Release memory before loading the next run.
            del adata, kernel, fate, cbc
            gc.collect()

    # Combine trajectory results across all runs.
    combined_cbc = pd.concat(
        all_cbc,
        ignore_index=True,
    )

    combined_timing = pd.concat(
        all_timing,
        ignore_index=True,
    )

    combined_cbc.to_csv(
        summary_dir / "all_cbc.csv",
        index=False,
    )

    combined_timing.to_csv(
        summary_dir / "all_timing.csv",
        index=False,
    )

    print("\n===== All CBC Results =====")
    print(combined_cbc.to_string(index=False))

    print("\n===== All Timing Results =====")
    print(combined_timing.to_string(index=False))

    print("\n===== All Runs Complete =====")
    print("Results saved to:", output_dir)


if __name__ == "__main__":
    freeze_support()
    main()
