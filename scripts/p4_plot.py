"""Plot the paired divergence of TF effects, with one point per repeat."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
directory = ROOT / "results/p4"
data = pd.read_csv(directory / "shared_cell_effects_by_repeat.csv")
data = data.loc[data.source_group == "all_shared_nonterminal"]
endpoints = ["PAX7+MyoProg1", "MYH3+MyoC", "MYL3+MyoC"]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
for ax, tf in zip(axes, ["MSC", "MYOG"]):
    for repeat, color in enumerate(["#2a678c", "#ae5929", "#28735a"]):
        values = data.loc[(data.regulator == tf) & (data["repeat"] == repeat)].set_index("endpoint")
        ax.scatter([x + 0.11 * (repeat-1) for x in range(3)],
                   100 * values.loc[endpoints, "targeted_minus_random_mae"],
                   color=color, label=f"Repeat {repeat}", s=38)
    ax.axhline(0, color="0.4", lw=1)
    ax.set_xticks(range(3), endpoints, rotation=15)
    ax.set_title(tf)
    ax.margins(x=0.15)
axes[0].set_ylabel("Targeted − random mean absolute\nchange from full (percentage points)")
axes[1].legend(frameon=False)
fig.tight_layout()
fig.savefig(directory / "perturbation_stability.png", dpi=180)
plt.close(fig)

fig, axes = plt.subplots(2, 3, figsize=(10, 5.6), sharex=True, sharey=True)
for row, tf in enumerate(["MSC", "MYOG"]):
    for col, endpoint in enumerate(endpoints):
        ax = axes[row, col]
        subset = data.loc[(data.regulator == tf) & (data.endpoint == endpoint)].sort_values("repeat")
        for (_, run), color in zip(subset.iterrows(), ["#2a678c", "#ae5929", "#28735a"]):
            shifts = 100 * run[["mean_full", "mean_targeted", "mean_random"]].to_numpy(dtype=float)
            ax.plot(range(3), shifts, "o-", color=color, lw=1, label=f"Repeat {run['repeat']}")
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xticks(range(3), ["Full", "Targeted", "Random"])
        if row == 0:
            ax.set_title(endpoint)
        if col == 0:
            ax.set_ylabel(f"{tf}\nMean fate change (pp)")
axes[0, -1].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(directory / "signed_fate_changes.png", dpi=180)
plt.close(fig)
