
from pathlib import Path
import json

import anndata as ad
import cellrank as cr
import numpy as np
import pandas as pd


ENDPOINTS = (
    "PAX7+MyoProg1",
    "MYH3+MyoC",
    "MYL3+MyoC",
)


def load_terminal_anchors(path: str | Path) -> dict:
    """Load fixed terminal cell IDs."""
    with open(path, "r", encoding="utf-8") as f:
        anchors = json.load(f)

    if set(anchors) != set(ENDPOINTS):
        raise ValueError("Unexpected terminal endpoints.")

    return {endpoint: anchors[endpoint] for endpoint in ENDPOINTS}


def build_kernel(
    adata: ad.AnnData,
    n_jobs: int = 1,
):
    """Construct the combined CellRank transition kernel."""

    vk = cr.kernels.VelocityKernel(
        adata,
        xkey="Ms",
        vkey="velocity",
    )

    vk.compute_transition_matrix(n_jobs=n_jobs)

    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    kernel = 0.8 * vk + 0.2 * ck
    kernel.compute_transition_matrix()

    return kernel


def infer_fate_probabilities(
    adata: ad.AnnData,
    anchors: dict,
    n_jobs: int = 1,
    kernel=None,
) -> pd.DataFrame:
    """Compute fate probabilities using fixed terminal anchors."""

    if not adata.obs_names.is_unique:
        raise ValueError("Cell IDs must be unique.")

    available_ids = set(adata.obs_names)

    for endpoint, cell_ids in anchors.items():
        missing = set(cell_ids) - available_ids

        if missing:
            raise ValueError(
                f"{endpoint}: {len(missing)} missing anchors."
            )

    # Build the transition matrix.
    if kernel is None:
        kernel = build_kernel(adata, n_jobs=n_jobs)

    # Set terminal states without recomputing macrostates.
    estimator = cr.estimators.GPCCA(kernel)

    estimator.set_terminal_states(anchors)

    # Compute absorption probabilities using a sparse solver.
    estimator.compute_fate_probabilities(
        solver="gmres",
        use_petsc=False,
        n_jobs=1,
        tol=1e-6,
    )

    fate = estimator.fate_probabilities

    probabilities = pd.DataFrame(
        np.asarray(fate),
        index=adata.obs_names,
        columns=fate.names,
    )

    # Preserve the agreed endpoint order.
    probabilities = probabilities.loc[:, list(ENDPOINTS)]
    probabilities.index.name = "cell_id"

    # Validate numerical results.
    values = probabilities.to_numpy()

    if not np.isfinite(values).all():
        raise ValueError("Non-finite fate probabilities.")

    if values.min() < -1e-6 or values.max() > 1 + 1e-6:
        raise ValueError("Fate probabilities outside valid bounds.")

    if not np.allclose(
        values.sum(axis=1),
        1.0,
        atol=1e-3,
        rtol=0,
    ):
        raise ValueError("Fate probabilities do not sum to one.")

    # Validate the absorbing states.
    for endpoint, cell_ids in anchors.items():
        if not np.allclose(
            probabilities.loc[cell_ids, endpoint],
            1.0,
            atol=1e-3,
            rtol=0,
        ):
            raise ValueError(
                f"Invalid absorption probabilities: {endpoint}"
            )

    return probabilities
