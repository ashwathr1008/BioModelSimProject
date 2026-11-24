
from pathlib import Path

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]

ENDPOINTS = (
    "PAX7+MyoProg1",
    "MYH3+MyoC",
    "MYL3+MyoC",
)


def main():
    run = "full_r0"

    data_path = (
        PROJECT_ROOT / "results/p2/X"
        / run / "fitted.h5ad"
    )

    fate_path = (
        PROJECT_ROOT / "results/p3/fate"
        / f"{run}_fate.csv"
    )

    output_dir = PROJECT_ROOT / "results/p3/figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the fitted dataset and fate probabilities.
    adata = ad.read_h5ad(data_path)
    fate = pd.read_csv(fate_path, index_col="cell_id")

    # Verify cell IDs before aligning the results.
    if not fate.index.is_unique:
        raise ValueError("Duplicated fate cell IDs.")

    if set(fate.index) != set(adata.obs_names):
        raise ValueError("Cell ID mismatch.")

    fate = fate.loc[adata.obs_names]

    # Use the reference UMAP for consistent visualization.
    umap = adata.obsm["X_umap_ref"]

    # Create the cell-type reference plot.
    fig, ax = plt.subplots(figsize=(10, 8))

    cell_types = adata.obs["leiden_R_celltype"]

    for cell_type in cell_types.unique():
        mask = (cell_types == cell_type).to_numpy()

        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=5,
            alpha=0.7,
            label=cell_type,
        )

    ax.set_title("Human Limb Development - Cell Types")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        markerscale=2,
        frameon=False,
    )

    ax.set_aspect("equal")
    fig.tight_layout()

    fig.savefig(
        output_dir / f"{run}_cell_types.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    # Plot the probability of reaching each endpoint.
    fig, axes = plt.subplots(
        1, 3,
        figsize=(18, 6),
        constrained_layout=True,
    )

    for ax, endpoint in zip(axes, ENDPOINTS):

        probabilities = fate[endpoint].to_numpy()

        # Draw high-probability cells last for visibility.
        order = probabilities.argsort()

        scatter = ax.scatter(
            umap[order, 0],
            umap[order, 1],
            c=probabilities[order],
            cmap="viridis",
            vmin=0,
            vmax=1,
            s=5,
        )

        ax.set_title(endpoint)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal")

    fig.colorbar(
        scatter,
        ax=axes.tolist(),
        label="Fate probability",
        shrink=0.8,
    )

    fig.suptitle("Cell Fate Probabilities - Full r0")

    fig.savefig(
        output_dir / f"{run}_fate_umap.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print("Visualization complete.")
    print("Saved to:", output_dir)


if __name__ == "__main__":
    main()
