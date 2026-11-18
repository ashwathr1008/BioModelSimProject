# Person 1, step 3: rebuild the model inputs for every manifest and run QC on each.
# Writes data/prepared/{condition}_r{repeat}.h5ad (gitignored; hand these to Person 2) and
# QC files under results/p1/qc/. First confirms our smoothing recipe reproduces the shipped Ms/Mu.
# Run: python scripts/p1_03_build_inputs.py [--only full_r0 targeted_r0] [--pca-input shipped_pca]
#      (or run the cells in VS Code interactive mode)

# %%
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    _here = Path(__file__).resolve().parent
except NameError:
    _here = Path.cwd()
REPO = next(p for p in [_here, *_here.parents] if (p / "config" / "experiment.yaml").exists())
sys.path.insert(0, str(REPO / "src"))

from p1_data.common import dense, load_config, load_source, repo_path
from p1_data.inputs import minmax, moments_from_graph
from p1_data.rebuild import build_condition, check_condition

CONFIG = None
ONLY = None        # e.g. ["full_r0", "targeted_r0"]; None builds every manifest
PCA_INPUT = None   # None uses the config; or "X", "log1p_spliced", "shipped_pca"

ap = argparse.ArgumentParser(allow_abbrev=False)
ap.add_argument("--config", default=CONFIG)
ap.add_argument("--only", nargs="*", default=ONLY)
ap.add_argument("--pca-input", default=PCA_INPUT, choices=["X", "log1p_spliced", "shipped_pca"])
args, _ = ap.parse_known_args()

# %%
# Load config, set output folders (a non-default PCA input writes to separate folders).
cfg = load_config(args.config)
res = repo_path(cfg["output"]["results_dir"])
prep = repo_path(cfg["output"]["prepared_dir"])
qc = res / "qc"
if args.pca_input and args.pca_input != cfg["preprocessing"]["pca_input"]:
    cfg["preprocessing"]["pca_input"] = args.pca_input
    prep = prep.with_name(f"{prep.name}_{args.pca_input}")
    qc = res / f"qc_{args.pca_input}"
qc.mkdir(parents=True, exist_ok=True)
prep.mkdir(parents=True, exist_ok=True)
src = load_source(cfg)

# %%
# Recipe check: the shipped graph + our smoothing must reproduce the shipped Ms/Mu.
for lay, raw in [("Ms", "spliced"), ("Mu", "unspliced")]:
    rec = minmax(moments_from_graph(src.obsp["connectivities"], dense(src.layers[raw])))
    err = float(np.abs(rec - dense(src.layers[lay])).max())
    print(f"recipe check {lay}: max abs error vs shipped = {err:.2e}")
    if err > 1e-4:
        raise RuntimeError("moments recipe does not reproduce shipped layers; stop and investigate")

# %%
# Find the manifests to build.
files = sorted((res / "manifests").glob("*.csv.gz"))
if args.only:
    files = [f for f in files if f.name.replace(".csv.gz", "") in set(args.only)]
if not files:
    raise RuntimeError("no manifests found; run p1_02_make_manifests.py first")
print([f.name for f in files])

# %%
# Build, check, and save each condition.
rows = []
for f in files:
    cid = f.name.replace(".csv.gz", "")
    man = pd.read_csv(f)
    print(f"\n== {cid}: {(man.status == 'retained').sum()} retained ==")
    sub = build_condition(src, man, cfg)
    checks, bs = check_condition(sub, src, cfg)
    sub.write_h5ad(prep / f"{cid}.h5ad", compression="gzip")
    with open(qc / f"{cid}_checks.json", "w") as fh:
        json.dump(checks, fh, indent=2)
    bs.to_csv(qc / f"{cid}_boundary_support.csv", index=False)
    print(json.dumps({k: v for k, v in checks.items() if k != "celltype_counts"}, indent=1))
    print(bs.to_string(index=False))
    flat = {k: v for k, v in checks.items() if not isinstance(v, (dict, list))}
    flat["cid"] = cid
    for _, b in bs.iterrows():
        flat[f"boundary_{b.source}->{b.target}"] = b.n_source_with_target_neighbor
    rows.append(flat)

# %%
# Save the combined QC table.
qc_summary = pd.DataFrame(rows).set_index("cid")
qc_summary.to_csv(qc / "qc_summary.csv")
print(qc_summary)
print(f"Wrote prepared inputs to {prep} and QC to {qc}")

# %%
