"""Summarize TF effects on the same retained nonterminal cells within each repeat."""
from pathlib import Path
import json
import hashlib

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
P4 = ROOT / "results/p4"
ENDPOINTS = ["PAX7+MyoProg1", "MYH3+MyoC", "MYL3+MyoC"]

def depletion_auc(frame):
    n = len(frame)
    return mannwhitneyu(frame.p_unperturbed, frame.p_perturbed).statistic / (n*n)


def summarize_tutorial():
    folder = P4 / "tutorial"
    # Published values in the human limb tutorial, ordered by ENDPOINTS.
    published = {"MSC": [0.671613, 0.327251, 0.439194], "MYOG": [0.218008, 0.784527, 0.747119]}
    rows = []
    for tf, scores in published.items():
        table = pd.read_csv(folder / f"{tf}.csv.gz")
        for endpoint, score in zip(ENDPOINTS, scores):
            cells = table.loc[table.endpoint == endpoint]
            rows.append(dict(regulator=tf, endpoint=endpoint, n_cells=len(cells), fixed_anchor_auc=depletion_auc(cells),
                             published_tutorial_auc=score, mean_delta_p=cells.delta_p.mean()))
    pd.DataFrame(rows).to_csv(folder / "depletion_check.csv", index=False)


def summarize():
    rows = []
    scores = []
    manifest = []
    for repeat in range(3):
        for tf in ["MSC", "MYOG"]:
            runs = {}
            for cond in ["full", "targeted", "random"]:
                folder = P4 / f"{cond}_r{repeat}"
                table = pd.read_csv(folder / f"{tf}.csv.gz")
                info = json.loads((folder / f"{tf}.json").read_text())
                if not ((table.condition == cond) & (table["repeat"] == repeat) & (table.regulator == tf)).all():
                    raise ValueError(f"Incorrect labels in {folder.name}/{tf}")
                if set(table.endpoint) != set(ENDPOINTS) or len(table) != 3 * info["n_cells"]:
                    raise ValueError(f"Incomplete endpoints or cells in {folder.name}/{tf}")
                values = table[["p_unperturbed", "p_perturbed", "delta_p"]].to_numpy()
                if not np.isfinite(values).all() or not np.allclose(values[:, 1] - values[:, 0], values[:, 2], atol=1e-12):
                    raise ValueError(f"Invalid cell-level changes in {folder.name}/{tf}")
                if not (table.groupby("cell_id").size() == len(ENDPOINTS)).all():
                    raise ValueError(f"Missing cell/endpoint combinations in {folder.name}/{tf}")
                path = folder / f"{tf}.csv.gz"
                manifest.append(dict(condition=cond, repeat=repeat, regulator=tf, cells=info["n_cells"],
                                     rows=len(table), seconds=info["seconds"], path=str(path.relative_to(ROOT)),
                                     sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
                runs[cond] = table.set_index(["cell_id", "endpoint"], verify_integrity=True)
                for endpoint, subset in table.groupby("endpoint"):
                    n = len(subset)
                    auc = depletion_auc(subset)
                    scores.append(dict(condition=cond, repeat=repeat, regulator=tf, endpoint=endpoint,
                                       n_cells=n, depletion_auc=auc))
            shared = runs["full"].index.intersection(runs["targeted"].index).intersection(runs["random"].index)
            frame = runs["full"].loc[shared, ["cell_type", "delta_p"]].rename(columns={"delta_p": "full"})
            for cond in ["targeted", "random"]:
                reduced = runs[cond].loc[shared]
                if not frame.cell_type.equals(reduced.cell_type):
                    raise ValueError("Different cell labels across conditions")
                frame[cond] = reduced.delta_p
            frame = frame.loc[~frame.cell_type.isin(ENDPOINTS)].reset_index()
            groups = [("all_shared_nonterminal", frame), *frame.groupby("cell_type")]
            for group, subset in groups:
                for endpoint, part in subset.groupby("endpoint"):
                    targeted = np.abs(part.targeted-part.full).mean()
                    random = np.abs(part.random-part.full).mean()
                    rows.append(dict(repeat=repeat, regulator=tf, endpoint=endpoint, source_group=group,
                                     n_cells=len(part), mean_full=part.full.mean(), mean_targeted=part.targeted.mean(),
                                     mean_random=part.random.mean(), median_full=part.full.median(),
                                     median_targeted=part.targeted.median(), median_random=part.random.median(),
                                     mae_targeted=targeted, mae_random=random,
                                     targeted_minus_random_mae=targeted-random))
    effects = pd.DataFrame(rows)
    effects.to_csv(P4 / "shared_cell_effects_by_repeat.csv", index=False)
    pooled = effects.loc[effects.source_group == "all_shared_nonterminal"]
    contrasts = pooled.groupby(["regulator", "endpoint"])["targeted_minus_random_mae"]
    intervals = contrasts.agg(["median", "min", "max", "count"])
    intervals["repeats_targeted_larger"] = contrasts.apply(lambda values: int((values > 0).sum()))
    intervals.to_csv(P4 / "paired_contrast_across_repeats.csv")
    pd.DataFrame(scores).to_csv(P4 / "depletion_by_run.csv", index=False)
    pd.DataFrame(manifest).to_csv(P4 / "run_manifest.csv", index=False)


if __name__ == "__main__":
    summarize()
    if all((P4 / "tutorial" / f"{tf}.json").exists() for tf in ("MSC", "MYOG")):
        summarize_tutorial()
