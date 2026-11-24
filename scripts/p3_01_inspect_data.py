
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_PATH = PROJECT_ROOT / "results/p2/X/full_r0/fitted.h5ad"

adata = ad.read_h5ad(DATA_PATH)

# Inspect the dataset dimensions.
print("\n===== Dataset Information =====")
print("Number of cells:", adata.n_obs)
print("Number of genes:", adata.n_vars)

# Inspect the cell-type distribution.
print("\n===== Cell Types =====")
print(adata.obs["leiden_R_celltype"].value_counts())

# Inspect the developmental stages.
print("\n===== Developmental Stages =====")
print(adata.obs["adj_stage"].value_counts())

# Inspect terminal populations.
print("\n===== Terminal Populations =====")

terminal_states = [
    "PAX7+MyoProg1",
    "MYH3+MyoC",
    "MYL3+MyoC",
]

for state in terminal_states:
    mask = adata.obs["leiden_R_celltype"] == state
    print(f"{state}: {mask.sum()} cells")

# Inspect the RNA velocity matrix.
print("\n===== RNA Velocity =====")

velocity = adata.layers["velocity"]

print("Shape:", velocity.shape)
print("All values finite:", np.isfinite(velocity).all())
print("Minimum:", velocity.min())
print("Maximum:", velocity.max())

# Inspect the cell connectivity graph.
print("\n===== Connectivity Graph =====")
print("Shape:", adata.obsp["connectivities"].shape)

# Evaluate latent time against developmental stage.
print("\n===== Latent Time Evaluation =====")

rho, p_value = spearmanr(
    adata.obs["latent_time"],
    adata.obs["adj_stage_num"],
    nan_policy="omit",
)

print(f"Spearman correlation: {rho:.4f}")
print(f"P-value: {p_value:.4e}")

# Inspect cell identifiers.
print("\n===== Cell IDs =====")
print(adata.obs_names[:5].tolist())
print("Unique cell IDs:", adata.obs_names.is_unique)
