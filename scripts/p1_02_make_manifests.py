# Person 1, step 2: define every experimental condition as an exact list of cell IDs.
# Writes one manifest per condition x repeat (results/p1/manifests/*.csv.gz), summary tables,
# and sampling QC figures (results/p1/figures/). Commit the manifests: they define the experiment.
# Run: python scripts/p1_02_make_manifests.py   (or run the cells in VS Code interactive mode)

# %%
import argparse
import sys
from pathlib import Path

import matplotlib

INTERACTIVE = "ipykernel" in sys.modules
if not INTERACTIVE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

try:
    _here = Path(__file__).resolve().parent
except NameError:
    _here = Path.cwd()
REPO = next(p for p in [_here, *_here.parents] if (p / "config" / "experiment.yaml").exists())
sys.path.insert(0, str(REPO / "src"))

from p1_data.common import condition_id, load_config, load_source, repo_path
from p1_data.sampling import make_manifest, n_to_remove, removed_composition, summarize

CONFIG = None

ap = argparse.ArgumentParser(allow_abbrev=False)
ap.add_argument("--config", default=CONFIG)
args, _ = ap.parse_known_args()

# %%
# Load config and cell metadata.
cfg = load_config(args.config)
s = cfg["sampling"]
res = repo_path(cfg["output"]["results_dir"])
(res / "manifests").mkdir(parents=True, exist_ok=True)
(res / "figures").mkdir(parents=True, exist_ok=True)
src = load_source(cfg)
obs = src.obs

# %%
# Feasibility check: enough target cells must remain after removal.
target = s["target_population"]
n_target = int((obs[cfg["data"]["celltype_key"]] == target).sum())
n_rm = n_to_remove(obs, cfg)
retained_target = n_target - n_rm
print(f"Target {target}: {n_target} cells; removing {n_rm} ({s['removal_fraction']:.0%}); {retained_target} retained")
if retained_target < s["min_retained_target_cells"]:
    raise RuntimeError(f"FEASIBILITY FAIL: {retained_target} < min_retained_target_cells={s['min_retained_target_cells']}")

# %%
# Build and save one manifest per condition x repeat.
conditions = list(s["conditions"]) + list(s.get("optional_conditions") or [])
manifests = []
for r in range(s["n_repeats"]):
    for c in conditions:
        m = make_manifest(obs, c, r, cfg)
        m.to_csv(res / "manifests" / f"{condition_id(c, r)}.csv.gz", index=False)
        manifests.append(m)

# %%
# Summary tables: counts per manifest, and what was removed by cell type / stage / sample.
summ = summarize(manifests, cfg)
summ.to_csv(res / "sampling_summary.csv", index=False)
print(summ.to_string(index=False))
for by in ["celltype", "stage", "sample"]:
    removed_composition(manifests, by).to_csv(res / f"removed_by_{by}.csv", index=False)

# %%
# How much the removed sets overlap between repeats.
rows = []
for c in conditions:
    if c == "full":
        continue
    sets = {m.repeat.iloc[0]: set(m.cell_id[m.status == "removed"]) for m in manifests if m.condition.iloc[0] == c}
    for a in sets:
        for b in sets:
            if a < b:
                rows.append(dict(condition=c, repeat_a=a, repeat_b=b,
                                 jaccard=len(sets[a] & sets[b]) / len(sets[a] | sets[b])))
overlap = pd.DataFrame(rows)
overlap.to_csv(res / "removed_set_overlap_between_repeats.csv", index=False)
print(overlap)

# %%
# Figures for repeat 0: removed cells by stage and cell type, and retained cells per condition.
m0 = [m for m in manifests if m.repeat.iloc[0] == 0 and m.condition.iloc[0] != "full"]
for by in ["stage", "celltype"]:
    tab = pd.concat(
        {m.condition.iloc[0]: m[m.status == "removed"][by].value_counts() for m in m0}, axis=1
    ).fillna(0).sort_index()
    ax = tab.plot.bar(figsize=(9, 4))
    ax.set_ylabel("cells removed (repeat 0)")
    ax.set_title(f"Removed cells by {by}: targeted vs matched random")
    plt.tight_layout()
    plt.savefig(res / "figures" / f"removed_by_{by}_r0.png", dpi=150)
    plt.show() if INTERACTIVE else plt.close()

tab = pd.concat(
    {m.condition.iloc[0]: m[m.status == "retained"].celltype.value_counts()
     for m in manifests if m.repeat.iloc[0] == 0}, axis=1
).fillna(0)
ax = tab.plot.bar(figsize=(10, 4))
ax.set_ylabel("cells retained (repeat 0)")
ax.set_title("Retained cells by cell type per condition")
plt.tight_layout()
plt.savefig(res / "figures" / "retained_by_celltype_r0.png", dpi=150)
plt.show() if INTERACTIVE else plt.close()
print(f"Wrote manifests, summaries and figures to {res}")

# %%
# Figure: where the removed cells sit on the shipped UMAP (full data vs targeted vs random, repeat 0).
ct = cfg["data"]["celltype_key"]
umap = src.obsm["X_umap"]
cats = list(obs[ct].cat.categories)
colors = dict(zip(cats, src.uns[f"{ct}_colors"]))
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for c in cats:
    sel = (obs[ct] == c).to_numpy()
    axes[0].scatter(umap[sel, 0], umap[sel, 1], s=2, color=colors[c], label=c, rasterized=True)
axes[0].set_title(f"Full data ({obs.shape[0]} cells)")
axes[0].legend(markerscale=5, fontsize=7, loc="best", frameon=False)
for ax, cond in zip(axes[1:], ["targeted", "random"]):
    m = next(m for m in manifests if m.condition.iloc[0] == cond and m.repeat.iloc[0] == 0)
    gone = (m.status == "removed").to_numpy()
    ax.scatter(umap[~gone, 0], umap[~gone, 1], s=2, color="lightgrey", rasterized=True, label="retained")
    ax.scatter(umap[gone, 0], umap[gone, 1], s=4, color="crimson", rasterized=True, label="removed")
    ax.set_title(f"{cond}: {gone.sum()} cells removed (repeat 0)")
    ax.legend(markerscale=5, fontsize=8, frameon=False)
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
plt.tight_layout()
plt.savefig(res / "figures" / "umap_removed_cells_r0.png", dpi=150)
plt.show() if INTERACTIVE else plt.close()
