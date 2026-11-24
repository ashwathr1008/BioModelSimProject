
from pathlib import Path
from multiprocessing import freeze_support
import sys
import gc

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from p3_analysis.fate_inference import build_kernel

from p3_analysis.trajectory import compute_cell_level_cbc


CONDITIONS = ("full", "targeted", "random")
REPEATS = (0, 1, 2)


def main():

    input_dir = PROJECT_ROOT / "results/p2/X"

    output_dir = (
        PROJECT_ROOT / "results/p3/cell_level_cbc"
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []

    for repeat in REPEATS:
        for condition in CONDITIONS:

            run = f"{condition}_r{repeat}"

            print("\n" + "=" * 60)
            print(f"Processing: {run}")
            print("=" * 60)

            # Load the matching fitted dataset.
            adata = ad.read_h5ad(
                input_dir / run / "fitted.h5ad"
            )

            # Build the same kernel used in the previous analysis.
            kernel = build_kernel(adata, n_jobs=1)

            # Compute individual boundary-cell CBC scores.
            results = compute_cell_level_cbc(
                kernel,
                adata,
            )

            results["condition"] = condition
            results["repeat"] = repeat

            # Compare against previously saved CBC summaries.
            reference = pd.read_csv(
                PROJECT_ROOT
                / "results/p3/trajectory"
                / f"{run}_cbc.csv"
            )

            for _, row in reference.iterrows():

                mask = (
                    (results["source"] == row["source"])
                    & (results["target"] == row["target"])
                    & results["valid"]
                )

                scores = results.loc[mask, "cbc"]

                if len(scores) != row["n_valid"]:
                    raise ValueError(
                        f"Boundary count mismatch: {run}"
                    )

                if len(scores) and not np.isclose(
                    scores.mean(),
                    row["mean_cbc"],
                    atol=1e-6,
                    rtol=0,
                ):
                    raise ValueError(
                        f"CBC mean mismatch: {run}"
                    )

            print("CBC summary validation: PASSED")

            results.to_csv(
                output_dir / f"{run}_cell_cbc.csv",
                index=False,
            )

            all_results.append(results)

            print("Boundary records:", len(results))
            print(f"Completed: {run}")

            del adata, kernel, results
            gc.collect()

    # Combine all cell-level results.
    combined = pd.concat(
        all_results,
        ignore_index=True,
    )

    keys = [
        "condition",
        "repeat",
        "cell_id",
        "source",
        "target",
    ]

    if combined.duplicated(keys).any():
        raise ValueError("Duplicated CBC records.")

    combined.to_csv(
        output_dir / "all_cell_cbc.csv",
        index=False,
    )

    print("\n===== Cell-level CBC Complete =====")
    print("Total records:", len(combined))
    print("Saved to:", output_dir)


if __name__ == "__main__":
    freeze_support()
    main()
