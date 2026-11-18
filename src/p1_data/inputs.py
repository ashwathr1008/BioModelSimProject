# Checks on the shipped dataset: what each RNA layer really is, how the gene network is
# oriented, cell counts by type/stage/sample, and how many cells sit on each transition boundary.

# %%
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    from .common import dense
except ImportError:  # running this file directly in VS Code interactive mode
    from common import dense


# %%
# Undo scVelo's per-cell normalization to get raw counts back, and report how many are integers.
def recover_raw_counts(adata, layer: str, size_key: str) -> tuple[np.ndarray, float, float]:
    X = dense(adata.layers[layer])
    init = adata.obs[size_key].to_numpy(dtype=float)
    const = float(np.median(init))
    raw = X * init[:, None] / const
    nz = raw[raw > 0]
    frac_int = float(np.mean(np.isclose(nz, np.round(nz), atol=1e-2))) if nz.size else float("nan")
    return np.round(raw), const, frac_int


# Average each cell's values over its graph neighbours (same as scVelo pp.moments).
def moments_from_graph(connectivities, X: np.ndarray) -> np.ndarray:
    C = sp.csr_matrix(connectivities, copy=True)
    C.data[:] = 1.0
    C = C + sp.identity(C.shape[0], format="csr")
    C = sp.diags(1.0 / np.asarray(C.sum(1)).ravel()) @ C
    return np.asarray(C @ X)


# Scale each gene (column) to the range 0-1.
def minmax(Z: np.ndarray) -> np.ndarray:
    lo, hi = Z.min(0), Z.max(0)
    rng = np.where(hi - lo == 0, 1.0, hi - lo)
    return (Z - lo) / rng


# %%
# Verify what the spliced/unspliced/matrix/X/Ms/Mu layers contain and how the shipped graph was built.
def provenance_checks(adata) -> dict:
    out = {}
    rawS, cS, fS = recover_raw_counts(adata, "spliced", "initial_size_spliced")
    rawU, cU, fU = recover_raw_counts(adata, "unspliced", "initial_size_unspliced")
    out["spliced_is_per_cell_normalized_raw"] = fS
    out["spliced_norm_constant"] = cS
    out["unspliced_is_per_cell_normalized_raw"] = fU
    out["unspliced_norm_constant"] = cU
    if "matrix" in adata.layers and "ambiguous" in adata.layers:
        M = dense(adata.layers["matrix"])
        A = dense(adata.layers["ambiguous"])
        out["matrix_equals_S_plus_U_plus_ambiguous"] = float(np.mean(np.isclose(M, rawS + rawU + A)))
    x = dense(adata.X)
    out["X_min"] = float(x.min())
    out["X_max"] = float(x.max())
    out["X_equals_log1p_spliced_maxabs"] = float(np.abs(x - np.log1p(dense(adata.layers["spliced"]))).max())
    for lay, src in [("Ms", "spliced"), ("Mu", "unspliced")]:
        if lay in adata.layers and "connectivities" in adata.obsp:
            rec = minmax(moments_from_graph(adata.obsp["connectivities"], dense(adata.layers[src])))
            out[f"{lay}_reproduced_from_shipped_graph_maxabs"] = float(np.abs(rec - dense(adata.layers[lay])).max())
    if "X_pca" in adata.obsm:
        xc = x - x.mean(0)
        U, s, _ = np.linalg.svd(xc, full_matrices=False)
        pcs = U[:, :5] * s[:5]
        out["shipped_pca_vs_606gene_pca_abs_corr_pc1to5"] = [
            round(float(abs(np.corrcoef(pcs[:, k], adata.obsm["X_pca"][:, k])[0, 1])), 3) for k in range(5)
        ]
    out["neighbors_params"] = dict(adata.uns.get("neighbors", {}).get("params", {}))
    return out


# Check the prior gene network: shape, orientation, edge counts, and MSC/MYOG connections.
def grn_checks(adata, endpoints: list[str]) -> dict:
    sk = adata.uns["skeleton"]
    out = {"skeleton_type": type(sk).__name__}
    if isinstance(sk, pd.DataFrame):
        out["skeleton_index_matches_var"] = list(sk.index) == list(adata.var_names)
        out["skeleton_columns_match_var"] = list(sk.columns) == list(adata.var_names)
        v = sk.to_numpy()
    else:
        v = np.asarray(sk)
    tfs = set(adata.var_names[adata.var["TF"].astype(bool)])
    rows_nz = set(adata.var_names[v.sum(1) > 0])
    out.update(
        shape=list(v.shape),
        values=[float(u) for u in np.unique(v)],
        n_edges=int((v != 0).sum()),
        self_loops=int(np.trace(v != 0)),
        n_TF_annotated=len(tfs),
        n_rows_with_outgoing_edges=len(rows_nz),
        orientation_rows_are_regulators=rows_nz <= tfs,
        genes_without_incoming_edges=int((v.sum(0) == 0).sum()),
    )
    for g in ["MSC", "MYOG"]:
        if g in adata.var_names:
            i = list(adata.var_names).index(g)
            out[f"{g}_is_TF"] = g in tfs
            out[f"{g}_out_degree"] = int((v[i] != 0).sum())
            out[f"{g}_in_degree"] = int((v[:, i] != 0).sum())
    return out


# %%
# Count cells by cell type, stage, and sample.
def count_tables(adata, ct: str, stage: str, sample: str) -> dict[str, pd.DataFrame]:
    o = adata.obs
    return {
        "celltype_counts": o[ct].value_counts().rename_axis(ct).to_frame("n_cells"),
        "celltype_by_stage": pd.crosstab(o[ct], o[stage]),
        "celltype_by_sample": pd.crosstab(o[ct], o[sample]),
        "sample_by_stage": pd.crosstab(o[sample], o[stage]),
    }


# For each transition, count source cells with at least one neighbour in the target type (cells CBC can score).
def boundary_support(connectivities, labels: pd.Series, transitions) -> pd.DataFrame:
    C = sp.csr_matrix(connectivities)
    lab = labels.astype(str).to_numpy()
    rows = []
    for src, tgt in transitions:
        s_idx = np.where(lab == src)[0]
        tmask = (lab == tgt).astype(float)
        n_tgt_nbrs = np.asarray(C[s_idx].astype(bool).astype(float) @ tmask).ravel() if s_idx.size else np.array([])
        rows.append(
            dict(
                source=src,
                target=tgt,
                n_source=int(s_idx.size),
                n_target=int((lab == tgt).sum()),
                n_source_with_target_neighbor=int((n_tgt_nbrs > 0).sum()),
                frac_source_with_target_neighbor=float((n_tgt_nbrs > 0).mean()) if s_idx.size else float("nan"),
            )
        )
    return pd.DataFrame(rows)
