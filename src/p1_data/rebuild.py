# Rebuilds the model inputs for one sampling condition using only the retained cells:
# new PCA, new neighbour graph, new smoothed Ms/Mu layers, then QC checks on the result.

# %%
from __future__ import annotations

import json
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

try:
    from .common import dense
    from .inputs import boundary_support, minmax, moments_from_graph
except ImportError:  # running this file directly in VS Code interactive mode
    from common import dense
    from inputs import boundary_support, minmax, moments_from_graph

UNS_KEEP = ["skeleton", "regulators", "targets", "leiden_R_celltype_colors", "adj_stage_colors"]
REF_OBSM = {"X_pca": "X_pca_ref", "X_umap": "X_umap_ref", "X_draw_graph_fa": "X_draw_graph_fa_ref"}


# %%
# Choose the expression matrix that PCA is computed on (shipped X or log1p of spliced).
def _pca_input(sub, how: str):
    if how == "X":
        return sub.X.copy()
    if how == "log1p_spliced":
        S = sub.layers["spliced"]
        return S.log1p() if sp.issparse(S) else np.log1p(S)
    raise ValueError(f"unknown pca_input {how!r}")


# Build a clean dataset from retained cells and recompute PCA, graph, and Ms/Mu from scratch.
def build_condition(src: ad.AnnData, manifest: pd.DataFrame, cfg: dict) -> ad.AnnData:
    import scanpy as sc

    p = cfg["preprocessing"]
    keep = manifest.loc[manifest.status == "retained", "cell_id"].to_numpy()
    view = src[keep]

    sub = ad.AnnData(X=view.X.copy(), obs=view.obs.copy(), var=src.var.copy())
    for lay in ["spliced", "unspliced"]:
        sub.layers[lay] = view.layers[lay].copy()
    for k in UNS_KEEP:
        if k in src.uns:
            sub.uns[k] = src.uns[k].copy() if hasattr(src.uns[k], "copy") else src.uns[k]
    for k, new in REF_OBSM.items():
        if k in view.obsm:
            sub.obsm[new] = np.asarray(view.obsm[k]).copy()

    if p["pca_input"] == "shipped_pca":
        tmp = ad.AnnData(X=np.zeros((sub.n_obs, 1), dtype=np.float32))
        tmp.obs_names = sub.obs_names
        tmp.obsm["X_pca"] = sub.obsm["X_pca_ref"].copy()
        tmp.uns["pca"] = {}
    else:
        tmp = ad.AnnData(X=_pca_input(sub, p["pca_input"]))
        tmp.obs_names = sub.obs_names
        sc.pp.pca(tmp, n_comps=p["n_comps"], svd_solver="arpack", random_state=p["preprocessing_seed"])
    sc.pp.neighbors(
        tmp,
        n_neighbors=p["n_neighbors"],
        n_pcs=p["n_pcs"],
        use_rep="X_pca",
        method=p["neighbors_method"],
        metric=p["metric"],
        transformer=p.get("knn_transformer"),
        random_state=p["preprocessing_seed"],
    )
    sub.obsm["X_pca"] = tmp.obsm["X_pca"]
    sub.uns["pca"] = {k: v for k, v in tmp.uns.get("pca", {}).items() if k != "params"}
    sub.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {k: v for k, v in tmp.uns["neighbors"]["params"].items() if v is not None},
    }
    sub.obsp["connectivities"] = tmp.obsp["connectivities"]
    sub.obsp["distances"] = tmp.obsp["distances"]

    conn = sub.obsp["connectivities"]
    Ms = moments_from_graph(conn, dense(sub.layers["spliced"]))
    Mu = moments_from_graph(conn, dense(sub.layers["unspliced"]))
    if p["minmax_scale"]:
        Ms, Mu = minmax(Ms), minmax(Mu)
    sub.layers["Ms"] = Ms.astype(np.float32)
    sub.layers["Mu"] = Mu.astype(np.float32)

    m0 = manifest.iloc[0]
    for col in ["condition", "repeat", "sampling_seed", "training_seed"]:
        sub.obs[col] = m0[col]
    meta = dict(
        condition=str(m0["condition"]),
        repeat=int(m0["repeat"]),
        sampling_seed=int(m0["sampling_seed"]),
        training_seed=int(m0["training_seed"]),
        n_retained=int(len(keep)),
        n_removed=int((manifest.status == "removed").sum()),
        preprocessing=p,
        sampling=cfg["sampling"],
        endpoints=cfg["endpoints"],
        packages={k: version(k) for k in ["anndata", "scanpy", "numpy", "scipy", "pandas"]},
    )
    sub.uns["p1_prep_json"] = json.dumps(meta)
    return sub


# %%
# How similar the rebuilt neighbour graph is to the shipped one (mean Jaccard overlap per cell).
def knn_overlap_with_reference(sub: ad.AnnData, src: ad.AnnData) -> float:
    pos = pd.Series(np.arange(src.n_obs), index=src.obs_names)
    keep_pos = pos[sub.obs_names].to_numpy()
    ref = sp.csr_matrix(src.obsp["distances"])[keep_pos][:, keep_pos].tocsr()
    new = sp.csr_matrix(sub.obsp["distances"]).tocsr()
    jac = []
    for i in range(sub.n_obs):
        a = set(new.indices[new.indptr[i] : new.indptr[i + 1]])
        b = set(ref.indices[ref.indptr[i] : ref.indptr[i + 1]])
        if a or b:
            jac.append(len(a & b) / len(a | b))
    return float(np.mean(jac))


# Run QC on a rebuilt dataset: finite values, gene/network alignment, graph connectivity, counts, boundaries.
def check_condition(sub: ad.AnnData, src: ad.AnnData, cfg: dict) -> tuple[dict, pd.DataFrame]:
    ct = cfg["data"]["celltype_key"]
    out = {}
    for lay in ["spliced", "unspliced", "Ms", "Mu"]:
        out[f"{lay}_all_finite"] = bool(np.isfinite(dense(sub.layers[lay])).all())
    for lay in ["Ms", "Mu"]:
        Z = dense(sub.layers[lay])
        out[f"{lay}_n_constant_genes"] = int((Z.max(0) == Z.min(0)).sum())
        out[f"{lay}_range"] = [float(Z.min()), float(Z.max())]
    out["gene_order_matches_source"] = list(sub.var_names) == list(src.var_names)
    sk_s, sk_n = src.uns["skeleton"], sub.uns["skeleton"]
    out["skeleton_identical_to_source"] = bool(np.array_equal(np.asarray(sk_s), np.asarray(sk_n)))
    if isinstance(sk_n, pd.DataFrame):
        out["skeleton_aligned_to_var"] = list(sk_n.index) == list(sub.var_names) == list(sk_n.columns)
    ncomp, lab = connected_components(sub.obsp["connectivities"], directed=False)
    out["graph_n_components"] = int(ncomp)
    out["graph_largest_component_frac"] = float(np.bincount(lab).max() / sub.n_obs)
    out["endpoint_counts"] = {e: int((sub.obs[ct] == e).sum()) for e in cfg["endpoints"]}
    out["celltype_counts"] = {k: int(v) for k, v in sub.obs[ct].value_counts().items()}
    out["knn_jaccard_vs_shipped_graph"] = knn_overlap_with_reference(sub, src)
    bs = boundary_support(sub.obsp["connectivities"], sub.obs[ct], cfg["transitions"])
    return out, bs
