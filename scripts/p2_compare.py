"""Summarize the tutorial and configured PCA experiments, including incomplete runs."""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from p1_data.common import condition_id, load_config, repo_path
from p2_run_all import configured_pcas


PCA_LABELS = {"X": "Recomputed PCA", "shipped_pca": "Shipped PCA"}
COLORS = {"X": "#2166ac", "shipped_pca": "#d6604d", "tutorial": "#555555"}


def planned_runs(cfg, root):
    yield "tutorial", "tutorial", 0, root / "tutorial"
    for repeat in range(cfg["sampling"]["n_repeats"]):
        for pca in configured_pcas(cfg):
            for condition in cfg["sampling"]["conditions"]:
                yield pca, condition, repeat, root / pca / condition_id(condition, repeat)


def correlation(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 3 or not np.isfinite(a).all() or not np.isfinite(b).all():
        return np.nan
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def read_cells(path, cfg):
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        obs = adata.obs
        stage = obs[cfg["data"]["stage_key"]].astype(str).str.removeprefix("Pcw").astype(float)
        cells = pd.DataFrame({"latent_time": obs["latent_time"].astype(float), "stage": stage})
    finally:
        adata.file.close()
    if not cells.index.is_unique or not np.isfinite(cells.to_numpy()).all():
        raise ValueError("Cell IDs must be unique and latent time and stage must be finite")
    return cells


def collect_runs(cfg, root):
    rows, cells = [], {}
    fields = ("training_seed", "n_cells", "n_genes", "epochs", "elapsed_seconds",
              "peak_cpu_memory_mib", "gpu_peak_allocated_mib", "stage_spearman", "device",
              "elbo_train_last", "elbo_train_minimum", "elbo_validation_last", "elbo_validation_minimum")
    for pca, condition, repeat, folder in planned_runs(cfg, root):
        run_id = str(folder.relative_to(root)).replace("\\", "/")
        row = dict(run=run_id, pca_input=pca, condition=condition, repeat=repeat,
                   status="planned", error="", **{field: np.nan for field in fields})
        record_path = folder / "run.json"
        try:
            if record_path.exists():
                record = json.loads(record_path.read_text(encoding="utf-8"))
                row.update({field: record.get(field, np.nan) for field in fields})
                row.update(status=record["status"], error=record.get("error", ""))
            if row["status"] == "complete":
                missing = [name for name in ("fitted.h5ad", "model/model.pt", "training.csv")
                           if not (folder / name).is_file()]
                if missing:
                    row.update(status="missing outputs", error=", ".join(missing))
                else:
                    frame = read_cells(folder / "fitted.h5ad", cfg)
                    row["n_cells"] = len(frame)
                    row["stage_spearman"] = correlation(frame.latent_time, frame.stage)
                    cells[(pca, condition, repeat)] = frame
            if (folder / "training.csv").is_file():
                history = pd.read_csv(folder / "training.csv")
                for column in ("elbo_train", "elbo_validation"):
                    values = pd.to_numeric(history[column], errors="coerce").dropna()
                    row[f"{column}_last"] = values.iloc[-1] if len(values) else np.nan
                    row[f"{column}_minimum"] = values.min() if len(values) else np.nan
        except (OSError, ValueError, KeyError) as error:
            row.update(status="invalid outputs", error=f"{type(error).__name__}: {error}")
            cells.pop((pca, condition, repeat), None)
        rows.append(row)
    return pd.DataFrame(rows), cells


def align_cells(frames):
    """Use cell identities, never row positions, and require consistent stage labels."""
    common = frames[0].index
    for frame in frames[1:]:
        common = common.intersection(frame.index, sort=False)
    if len(common) < 3:
        raise ValueError("Fewer than three shared cells")
    aligned = [frame.loc[common] for frame in frames]
    for frame in aligned[1:]:
        if not np.array_equal(aligned[0].stage.to_numpy(), frame.stage.to_numpy()):
            raise ValueError("Developmental stages disagree for shared cell IDs")
    return aligned


def compare_cells(cfg, cells):
    stages, pairs, availability = [], [], []
    pca_variants = configured_pcas(cfg)
    groups = [("baseline_pilots", 0, [("tutorial", "tutorial", 0)]
               + [(pca, "full", 0) for pca in pca_variants])]
    for repeat in range(cfg["sampling"]["n_repeats"]):
        keys = [(pca, condition, repeat) for pca in pca_variants
                for condition in cfg["sampling"]["conditions"]]
        groups.append(("shared_repeat", repeat, keys))
    for scope, repeat, keys in groups:
        missing = [f"{pca}/{condition}_r{r}" for pca, condition, r in keys if (pca, condition, r) not in cells]
        if missing:
            availability.append(dict(scope=scope, repeat=repeat, status="incomplete", n_shared=0,
                                     detail="Missing complete outputs: " + ", ".join(missing)))
            continue
        try:
            aligned = align_cells([cells[key] for key in keys])
        except ValueError as error:
            availability.append(dict(scope=scope, repeat=repeat, status="invalid", n_shared=0, detail=str(error)))
            continue
        n_shared = len(aligned[0])
        availability.append(dict(scope=scope, repeat=repeat, status="complete", n_shared=n_shared, detail=""))
        metrics = {}
        for key, frame in zip(keys, aligned):
            metrics[key] = correlation(frame.latent_time, frame.stage)
            stages.append(dict(scope=scope, repeat=repeat, pca_input=key[0], condition=key[1],
                               n_shared=n_shared, stage_spearman=metrics[key]))
        for i, (left, a) in enumerate(zip(keys, aligned)):
            for right, b in zip(keys[i + 1:], aligned[i + 1:]):
                if scope == "shared_repeat" and left[0] != right[0] and left[1] != right[1]:
                    continue
                pairs.append(dict(scope=scope, repeat=repeat, left_pca=left[0], left_condition=left[1],
                                  right_pca=right[0], right_condition=right[1], n_shared=n_shared,
                                  latent_time_spearman=correlation(a.latent_time, b.latent_time),
                                  stage_spearman_right_minus_left=metrics[right] - metrics[left]))
    stage_cols = ["scope", "repeat", "pca_input", "condition", "n_shared", "stage_spearman"]
    pair_cols = ["scope", "repeat", "left_pca", "left_condition", "right_pca", "right_condition",
                 "n_shared", "latent_time_spearman", "stage_spearman_right_minus_left"]
    return (pd.DataFrame(stages, columns=stage_cols), pd.DataFrame(pairs, columns=pair_cols),
            pd.DataFrame(availability))


def save_plots(runs, stages, output, pca_variants=("X", "shipped_pca")):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    completed = runs.loc[runs.status == "complete"]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    baseline = stages.loc[stages.scope == "baseline_pilots"]
    if len(baseline):
        for x, (_, row) in enumerate(baseline.iterrows()):
            ax.scatter(x, row.stage_spearman, s=90, color=COLORS[row.pca_input])
        labels = {"tutorial": "Tutorial", **PCA_LABELS}
        ax.set_xticks(range(len(baseline)), [labels[pca] for pca in baseline.pca_input])
        ax.set_title(f"Seed-0 baseline fits: {int(baseline.iloc[0].n_shared):,} shared cells")
    else:
        ax.text(.5, .6, "Baseline comparison pending: selected completed fits required",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_title("Seed-0 baseline comparison")
    ax.set_ylim(-1.05, 1.05)
    ax.axhline(0, color="#bbbbbb", linewidth=.8)
    ax.set_ylabel("Signed stage Spearman correlation")
    fig.tight_layout()
    fig.savefig(output / "baseline_stage_correlation.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    shared = stages.loc[stages.scope == "shared_repeat"]
    for ax, condition in zip(axes, ("full", "targeted", "random")):
        subset = shared.loc[shared.condition == condition]
        for repeat, block in subset.groupby("repeat"):
            values = block.set_index("pca_input").reindex(pca_variants)
            ax.plot(range(len(pca_variants)), values.stage_spearman, marker="o", label=f"Repeat {repeat}")
        ax.set_xticks(range(len(pca_variants)), [PCA_LABELS[pca] for pca in pca_variants])
        ax.set_xlim(-.3, max(.3, len(pca_variants) - .7))
        ax.set_title(condition.capitalize())
        ax.axhline(0, color="#bbbbbb", linewidth=.8)
        ax.set_ylim(-1.05, 1.05)
        if subset.empty:
            ax.text(.5, .6, "Awaiting complete repeat blocks", ha="center", va="center", transform=ax.transAxes)
    axes[0].set_ylabel("Signed stage Spearman correlation")
    if len(shared):
        axes[-1].legend(frameon=False)
    fig.suptitle("Stage association on shared cells within each repeat")
    fig.tight_layout()
    fig.savefig(output / "shared_cell_stage_correlation.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, max(3, .32 * len(completed) + 1)))
    if len(completed):
        ax.barh(completed.run, completed.elapsed_seconds / 60,
                color=[COLORS[pca] for pca in completed.pca_input])
        ax.invert_yaxis()
    else:
        ax.text(.5, .5, "No completed fits", ha="center", va="center", transform=ax.transAxes)
        ax.set_yticks([])
    ax.set_xlabel("Elapsed minutes (loading, training, sampling, and saving)")
    ax.set_title("Completed-fit computational cost")
    fig.tight_layout()
    fig.savefig(output / "fit_runtime.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=repo_path("config/p2_x.yaml"))
    parser.add_argument("--results-dir", type=Path, default=repo_path("results/p2"))
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    cfg = load_config(args.config)
    pca_variants = configured_pcas(cfg)
    output = args.output_dir or args.results_dir / "comparison"
    output.mkdir(parents=True, exist_ok=True)
    runs, cells = collect_runs(cfg, args.results_dir)
    stages, pairs, availability = compare_cells(cfg, cells)
    runs.to_csv(output / "run_summary.csv", index=False)
    stages.to_csv(output / "shared_cell_stage.csv", index=False)
    pairs.to_csv(output / "paired_latent_time.csv", index=False)
    availability.to_csv(output / "comparison_status.csv", index=False)
    save_plots(runs, stages, output, pca_variants)
    notes = {
        "completed_runs": int((runs.status == "complete").sum()),
        "planned_runs": len(runs),
        "pca_variants": pca_variants,
        "evaluation": f"Baseline pilots use their shared cells. Each main repeat uses the intersection of all {len(pca_variants) * len(cfg['sampling']['conditions'])} selected fits ({', '.join(pca_variants)}; {', '.join(cfg['sampling']['conditions'])}), including endpoint cells.",
        "interpretation": "Signed stage correlation is a secondary biological check, not predictive accuracy. Latent-time Spearman measures rank stability. These are computational repeats, not independent biological replicates. Loss magnitudes across differently preprocessed inputs are not a model-selection score.",
        "pca_guidance": "Recomputed PCA follows the removal experiment without retaining full-data PCA axes. Shipped PCA is an optional sensitivity analysis closer to the tutorial and retains axes learned from all cells. Only configured variants are included. Do not select the variant by the size or significance of the thinning effect.",
        "pending_downstream": "Fate probabilities, transition CBC, and MSC/MYOG perturbation comparisons require the Person 3/4 outputs and are not assessed here.",
    }
    (output / "comparison_notes.json").write_text(json.dumps(notes, indent=2) + "\n", encoding="utf-8")
    print(runs[["run", "status", "stage_spearman", "elapsed_seconds"]].to_string(index=False))
    print(f"Saved comparisons to {output}")


if __name__ == "__main__":
    main()
