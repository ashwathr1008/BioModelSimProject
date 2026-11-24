
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


TRANSITIONS = (
    ("PAX3+MyoProg", "MyoB1", "targeted_pathway"),
    ("MyoB1", "MyoC1", "targeted_pathway"),
    ("PAX3+PAX7+MyoProg", "MyoB2", "control_pathway"),
    ("MyoB2", "MyoC2", "control_pathway"),
)


def compute_latent_time_score(
    adata: ad.AnnData,
) -> dict:
    """Evaluate latent time against developmental stage."""

    rho, p_value = spearmanr(
        adata.obs["latent_time"],
        adata.obs["adj_stage_num"],
        nan_policy="omit",
    )

    return {
        "spearman_rho": float(rho),
        "p_value": float(p_value),
        "n_cells": adata.n_obs,
    }


def compute_cbc_scores(
    kernel,
    adata: ad.AnnData,
) -> pd.DataFrame:
    """Calculate CBC separately for each transition."""

    records = []

    for source, target, pathway in TRANSITIONS:

        print(f"Computing CBC: {source} -> {target}")

        scores = kernel.cbc(
            source=source,
            target=target,
            cluster_key="leiden_R_celltype",
            rep="X_pca",
            graph_key="distances",
        )

        scores = np.asarray(scores).ravel()

        # Exclude non-finite scores from the summary.
        valid = scores[np.isfinite(scores)]

        records.append({
            "source": source,
            "target": target,
            "pathway": pathway,
            "mean_cbc": (
                float(np.mean(valid))
                if len(valid) else np.nan
            ),
            "median_cbc": (
                float(np.median(valid))
                if len(valid) else np.nan
            ),
            "std_cbc": (
                float(np.std(valid))
                if len(valid) else np.nan
            ),
            "n_valid": len(valid),
            "n_invalid": len(scores) - len(valid),
            "n_source": int(
                (adata.obs["leiden_R_celltype"] == source).sum()
            ),
            "n_target": int(
                (adata.obs["leiden_R_celltype"] == target).sum()
            ),
        })

    return pd.DataFrame(records)

def compute_cell_level_cbc(
    kernel,
    adata: ad.AnnData,
) -> pd.DataFrame:
    """Compute CBC scores with their corresponding cell IDs."""

    records = []

    for source, target, pathway in TRANSITIONS:

        print(f"Cell-level CBC: {source} -> {target}")

        # Identify boundary cells in CellRank's original order.
        boundary_ids = kernel._get_boundary(
            source=source,
            target=target,
            cluster_key="leiden_R_celltype",
            graph_key="distances",
        )

        # Compute CBC using the original CellRank implementation.
        scores = kernel.cbc(
            source=source,
            target=target,
            cluster_key="leiden_R_celltype",
            rep="X_pca",
            graph_key="distances",
        )

        scores = np.asarray(scores).ravel()

        # Verify alignment between cell IDs and CBC scores.
        if len(boundary_ids) != len(scores):
            raise ValueError(
                f"Boundary-score mismatch: {source} -> {target}"
            )

        cell_ids = adata.obs_names[
            np.asarray(boundary_ids, dtype=int)
        ]

        for cell_id, score in zip(cell_ids, scores):

            records.append({
                "cell_id": cell_id,
                "source": source,
                "target": target,
                "pathway": pathway,
                "cbc": float(score),
                "valid": bool(np.isfinite(score)),
            })

    return pd.DataFrame(records)