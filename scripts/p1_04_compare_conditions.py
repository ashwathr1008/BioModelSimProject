# Person 1, step 4: compare the prepared datasets before and after removing cells.
# Plots (1) how many cells remain on each transition boundary, (2) how much the smoothed model
# inputs (Ms) change for cells kept in every condition, and (3) MYOG/MSC phase portraits.
# Needs steps 2 and 3 first. Figures go to results/p1/figures/.
# Run: python scripts/p1_04_compare_conditions.py   (or run the cells in VS Code interactive mode)

# %%
import argparse
import sys
from pathlib import Path

import matplotlib

INTERACTIVE = "ipykernel" in sys.modules
if not INTERACTIVE:
    matplotlib.use("Agg")
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    _here = Path(__file__).resolve().parent
except NameError:
    _here = Path.cwd()
REPO = next(p for p in [_here, *_here.parents] if (p / "config" / "experiment.yaml").exists())
sys.path.insert(0, str(REPO / "src"))

from p1_data.common import condition_id, dense, load_config, repo_path

CONFIG = None
GENES = ["MYOG", "MSC"]   # genes for the phase portraits

ap = argparse.ArgumentParser(allow_abbrev=False)
ap.add_argument("--config", default=CONFIG)
ap.add_argument("--genes", nargs="*", default=GENES)
args, _ = ap.parse_known_args()

# %%
# Load config and set paths.
cfg = load_config(args.config)
ct = cfg["data"]["celltype_key"]
res = repo_path(cfg["output"]["results_dir"])
prep = repo_path(cfg["output"]["prepared_dir"])
figs = res / "figures"
figs.mkdir(parents=True, exist_ok=True)
repeats = range(cfg["sampling"]["n_repeats"])
reduced = [c for c in cfg["sampling"]["conditions"] if c != "full"]
COND_COLORS = {"full": "#4c72b0", "random": "#dd8452", "targeted": "#c44e52"}


# Save a figure and show it in interactive mode.
def finish(fig, name):
    fig.tight_layout()
    fig.savefig(figs / name, dpi=150)
    plt.show() if INTERACTIVE else plt.close(fig)
    print("saved", figs / name)


# Load one prepared dataset (e.g. "targeted_r0").
def load_prepared(cid):
    return ad.read_h5ad(prep / f"{cid}.h5ad")


# %%
# Figure 1: boundary support per transition (cells CBC can score), full vs random vs targeted.
bs = []
for c in ["full"] + reduced:
    for r in repeats:
        t = pd.read_csv(res / "qc" / f"{condition_id(c, r)}_boundary_support.csv")
        t["condition"], t["repeat"] = c, r
        bs.append(t)
bs = pd.concat(bs)
bs["transition"] = bs.source + " → " + bs.target
agg = bs.groupby(["transition", "condition"]).n_source_with_target_neighbor.agg(["mean", "min", "max"])
transitions = [f"{s} → {t}" for s, t in cfg["transitions"]]

fig, ax = plt.subplots(figsize=(10, 4.5))
width = 0.8 / (1 + len(reduced))
for i, c in enumerate(["full"] + reduced):
    a = agg.xs(c, level="condition").loc[transitions]
    x = np.arange(len(transitions)) + i * width
    ax.bar(x, a["mean"], width, color=COND_COLORS.get(c), label=c,
           yerr=[a["mean"] - a["min"], a["max"] - a["mean"]], capsize=3)
ax.set_xticks(np.arange(len(transitions)) + width * len(reduced) / 2)
ax.set_xticklabels(transitions, rotation=15, ha="right")
ax.set_ylabel("source cells with a neighbour\nin the target type")
ax.set_title("Boundary support per transition (bars: mean over repeats, whiskers: range)")
ax.legend(frameon=False)
finish(fig, "boundary_support_by_condition.png")

# %%
# Figure 2: change in smoothed spliced input (Ms) vs full data, per cell type, for cells kept in all conditions.
full = load_prepared(condition_id("full", 0))
full_ms = pd.DataFrame(dense(full.layers["Ms"]), index=full.obs_names)
rows = []
for r in repeats:
    subs = {c: load_prepared(condition_id(c, r)) for c in reduced}
    shared = full.obs_names
    for c in reduced:
        shared = shared.intersection(subs[c].obs_names)
    for c in reduced:
        ms = pd.DataFrame(dense(subs[c].layers["Ms"]), index=subs[c].obs_names).loc[shared]
        delta = np.abs(ms.to_numpy() - full_ms.loc[shared].to_numpy()).mean(1)
        d = pd.DataFrame({"celltype": full.obs.loc[shared, ct].astype(str).to_numpy(), "delta": delta})
        for cell_type, v in d.groupby("celltype").delta.mean().items():
            rows.append(dict(condition=c, repeat=r, celltype=cell_type, mean_abs_delta_Ms=v))
    del subs
delta_tab = pd.DataFrame(rows)
delta_tab.to_csv(res / "delta_Ms_vs_full_by_celltype.csv", index=False)

order = list(full.obs[ct].cat.categories)
g = delta_tab.groupby(["celltype", "condition"]).mean_abs_delta_Ms.agg(["mean", "min", "max"])
fig, ax = plt.subplots(figsize=(11, 4.5))
width = 0.8 / len(reduced)
for i, c in enumerate(reduced):
    a = g.xs(c, level="condition").reindex(order)
    x = np.arange(len(order)) + i * width
    ax.bar(x, a["mean"], width, color=COND_COLORS.get(c), label=c,
           yerr=[a["mean"] - a["min"], a["max"] - a["mean"]], capsize=3)
ax.set_xticks(np.arange(len(order)) + width * (len(reduced) - 1) / 2)
ax.set_xticklabels(order, rotation=30, ha="right")
ax.set_ylabel("mean |ΔMs| vs full data\n(per cell, averaged over genes)")
ax.set_title("How much the model input changes for cells kept in every condition")
ax.legend(frameon=False)
finish(fig, "delta_Ms_vs_full_by_celltype.png")

# %%
# Figure 3: phase portraits (Ms vs Mu) for selected genes, full vs each reduced condition (repeat 0).
datasets = {"full": full, **{c: load_prepared(condition_id(c, 0)) for c in reduced}}
cats = list(full.obs[ct].cat.categories)
colors = dict(zip(cats, full.uns[f"{ct}_colors"]))
genes = [gname for gname in args.genes if gname in full.var_names]
fig, axes = plt.subplots(len(genes), len(datasets), figsize=(5 * len(datasets), 4.2 * len(genes)), squeeze=False)
for i, gname in enumerate(genes):
    for j, (name, a) in enumerate(datasets.items()):
        ax = axes[i, j]
        k = list(a.var_names).index(gname)
        ms, mu = dense(a.layers["Ms"])[:, k], dense(a.layers["Mu"])[:, k]
        order = np.random.default_rng(0).permutation(a.n_obs)   # random draw order so no cell type hides another
        point_colors = a.obs[ct].astype(str).map(colors).to_numpy()
        ax.scatter(ms[order], mu[order], s=2, c=point_colors[order], rasterized=True)
        ax.set_title(f"{gname}: {name} ({a.n_obs} cells)")
        ax.set_xlabel("Ms (smoothed spliced)")
        ax.set_ylabel("Mu (smoothed unspliced)")
handles = [plt.Line2D([], [], marker="o", ls="", color=colors[c], label=c) for c in cats]
fig.legend(handles=handles, loc="center right", fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0, 0.88, 1))
fig.savefig(figs / "phase_portraits_r0.png", dpi=150)
plt.show() if INTERACTIVE else plt.close(fig)
print("saved", figs / "phase_portraits_r0.png")

# %%
