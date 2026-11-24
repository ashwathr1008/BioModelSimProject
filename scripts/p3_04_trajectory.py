
from pathlib import Path
from multiprocessing import freeze_support
import sys

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from p3_analysis.fate_inference import build_kernel

from p3_analysis.trajectory import (
    compute_latent_time_score,
    compute_cbc_scores,
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

    output_dir = PROJECT_ROOT / "results/p3/trajectory"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the fitted dataset.
    adata = ad.read_h5ad(data_path)

    print("\n===== Dataset =====")
    print(adata.shape)

    # Evaluate latent time.
    print("\n===== Latent Time =====")

    timing = compute_latent_time_score(adata)

    print(timing)

    timing_df = pd.DataFrame([{
        "condition": condition,
        "repeat": repeat,
        **timing,
    }])

    # Construct the shared transition kernel.
    print("\n===== Building Kernel =====")

    kernel = build_kernel(adata, n_jobs=1)

    # Calculate per-transition CBC scores.
    print("\n===== CBC Analysis =====")

    cbc = compute_cbc_scores(kernel, adata)

    cbc["condition"] = condition
    cbc["repeat"] = repeat

    print(cbc.to_string(index=False))

    # Save trajectory analysis results.
    timing_df.to_csv(
        output_dir / f"{condition}_r{repeat}_timing.csv",
        index=False,
    )

    cbc.to_csv(
        output_dir / f"{condition}_r{repeat}_cbc.csv",
        index=False,
    )

    print("\n===== Trajectory Analysis Complete =====")
    print("Saved to:", output_dir)


if __name__ == "__main__":
    freeze_support()
    main()

