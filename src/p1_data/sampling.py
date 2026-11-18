# Builds the three sampling conditions (full, targeted, random) as explicit lists of
# retained and removed cell IDs, plus summary tables of what each condition removed.

# %%
from __future__ import annotations

import numpy as np
import pandas as pd

CONDITION_CODES = {"full": 0, "targeted": 1, "random": 2, "random_stagematched": 3}


# %%
# Make a seeded random number generator unique to one (repeat seed, condition) pair.
def _rng(sampling_seed: int, condition: str) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([sampling_seed, CONDITION_CODES[condition]]))


# Number of cells to remove: removal_fraction x size of the target population.
def n_to_remove(obs: pd.DataFrame, cfg: dict) -> int:
    ct = cfg["data"]["celltype_key"]
    n_target = int((obs[ct] == cfg["sampling"]["target_population"]).sum())
    return int(round(cfg["sampling"]["removal_fraction"] * n_target))


# Pick which cell IDs to remove for one condition and repeat (endpoint cells are never removed).
def removed_ids(obs: pd.DataFrame, condition: str, repeat: int, cfg: dict) -> np.ndarray:
    ct, st = cfg["data"]["celltype_key"], cfg["data"]["stage_key"]
    s = cfg["sampling"]
    seed = s["sampling_seed_base"] + repeat
    endpoints = set(cfg["endpoints"])
    target = s["target_population"]
    if target in endpoints:
        raise ValueError("target population must not be a protected endpoint")
    labels = obs[ct].astype(str)
    nonterminal = ~labels.isin(endpoints)
    n = n_to_remove(obs, cfg)

    if condition == "full":
        return np.array([], dtype=object)
    if condition == "targeted":
        pool = obs.index[labels == target].to_numpy()
        return _rng(seed, condition).choice(pool, size=n, replace=False)
    if condition == "random":
        pool = obs.index[nonterminal].to_numpy()
        return _rng(seed, condition).choice(pool, size=n, replace=False)
    if condition == "random_stagematched":
        tgt_removed = removed_ids(obs, "targeted", repeat, cfg)
        per_stage = obs.loc[tgt_removed, st].astype(str).value_counts()
        rng = _rng(seed, condition)
        picks = []
        for stage, k in sorted(per_stage.items()):
            pool = obs.index[nonterminal & (obs[st].astype(str) == stage)].to_numpy()
            if pool.size < k:
                raise ValueError(f"stage {stage}: need {k} nonterminal cells, only {pool.size}")
            picks.append(rng.choice(pool, size=k, replace=False))
        return np.concatenate(picks)
    raise ValueError(f"unknown condition {condition!r}")


# Build the manifest table: one row per cell with its status (retained/removed), labels, and seeds.
def make_manifest(obs: pd.DataFrame, condition: str, repeat: int, cfg: dict) -> pd.DataFrame:
    d = cfg["data"]
    removed = set(removed_ids(obs, condition, repeat, cfg))
    m = pd.DataFrame(
        {
            "cell_id": obs.index.to_numpy(),
            "condition": condition,
            "repeat": repeat,
            "sampling_seed": cfg["sampling"]["sampling_seed_base"] + repeat,
            "training_seed": cfg["model"]["training_seeds"][repeat],
            "celltype": obs[d["celltype_key"]].astype(str).to_numpy(),
            "stage": obs[d["stage_key"]].astype(str).to_numpy(),
            "sample": obs[d["sample_key"]].astype(str).to_numpy(),
        }
    )
    m["is_endpoint"] = m["celltype"].isin(cfg["endpoints"])
    m["status"] = np.where(m["cell_id"].isin(removed), "removed", "retained")
    assert not (m["is_endpoint"] & (m["status"] == "removed")).any(), "endpoint cell removed"
    return m


# %%
# One summary row per manifest: retained/removed counts for all cells, the target, and endpoints.
def summarize(manifests: list[pd.DataFrame], cfg: dict) -> pd.DataFrame:
    rows = []
    target = cfg["sampling"]["target_population"]
    for m in manifests:
        r = m[m.status == "retained"]
        rows.append(
            dict(
                condition=m.condition.iloc[0],
                repeat=int(m.repeat.iloc[0]),
                sampling_seed=int(m.sampling_seed.iloc[0]),
                training_seed=int(m.training_seed.iloc[0]),
                n_retained=len(r),
                n_removed=int((m.status == "removed").sum()),
                n_target_retained=int((r.celltype == target).sum()),
                n_target_removed=int(((m.status == "removed") & (m.celltype == target)).sum()),
                n_endpoint_retained=int(r.is_endpoint.sum()),
            )
        )
    return pd.DataFrame(rows)


# Count removed cells by cell type, stage, or sample for every condition and repeat.
def removed_composition(manifests: list[pd.DataFrame], by: str) -> pd.DataFrame:
    parts = []
    for m in manifests:
        rem = m[m.status == "removed"]
        if rem.empty:
            continue
        c = rem[by].value_counts().rename("n_removed").to_frame()
        c["condition"], c["repeat"] = m.condition.iloc[0], int(m.repeat.iloc[0])
        parts.append(c.reset_index().rename(columns={"index": by}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
