
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]

CONDITIONS = ("full", "targeted", "random")
REPEATS = (0, 1, 2)

ENDPOINTS = (
    "PAX7+MyoProg1",
    "MYH3+MyoC",
    "MYL3+MyoC",
)


def main():
    output_dir = PROJECT_ROOT / "results/p3/figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    for repeat in REPEATS:

        print(f"Processing repeat {repeat}")

        datasets = {}
        fate_tables = {}

        # Load datasets and align fate probabilities by cell ID.
        for condition in CONDITIONS:

            run = f"{condition}_r{repeat}"

            adata = ad.read_h5ad(
                PROJECT_ROOT
                / "results/p2/X"
                / run
                / "fitted.h5ad"
            )

            fate = pd.read_csv(
                PROJECT_ROOT
                / "results/p3/fate"
                / f"{run}_fate.csv",
                index_col="cell_id",
            )

            if not fate.index.is_unique:
                raise ValueError(f"Duplicated cell IDs: {run}")

            if set(fate.index) != set(adata.obs_names):
                raise ValueError(f"Cell ID mismatch: {run}")

            datasets[condition] = adata
            fate_tables[condition] = fate.loc[adata.obs_names]

        # Use identical UMAP limits across conditions.
        reference = datasets["full"].obsm["X_umap_ref"]

        xlim = (
            reference[:, 0].min() - 1,
            reference[:, 0].max() + 1,
        )

        ylim = (
            reference[:, 1].min() - 1,
            reference[:, 1].max() + 1,
        )

        fig, axes = plt.subplots(
            3,
            3,
            figsize=(15, 15),
            constrained_layout=True,
        )

        for i, condition in enumerate(CONDITIONS):

            adata = datasets[condition]
            fate = fate_tables[condition]

            umap = adata.obsm["X_umap_ref"]

            for j, endpoint in enumerate(ENDPOINTS):

                ax = axes[i, j]

                probabilities = fate[endpoint].to_numpy()
                order = np.argsort(probabilities)

                scatter = ax.scatter(
                    umap[order, 0],
                    umap[order, 1],
                    c=probabilities[order],
                    cmap="viridis",
                    vmin=0,
                    vmax=1,
                    s=4,
                )

                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.set_aspect("equal")

                if i == 0:
                    ax.set_title(endpoint)

                if j == 0:
                    ax.set_ylabel(
                        f"{condition.capitalize()}\nUMAP 2"
                    )

                if i == 2:
                    ax.set_xlabel("UMAP 1")

        fig.colorbar(
            scatter,
            ax=axes.ravel().tolist(),
            label="Fate probability",
            shrink=0.7,
        )

        fig.suptitle(
            f"Cell Fate Comparison - Repeat {repeat}",
            fontsize=18,
        )

        output_path = (
            output_dir
            / f"fate_comparison_r{repeat}.png"
        )

        fig.savefig(
            output_path,
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(fig)

        print("Saved:", output_path)


if __name__ == "__main__":
    main()
