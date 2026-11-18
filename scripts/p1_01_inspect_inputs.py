# Person 1, step 1: inspects the shipped human-limb dataset before using it.
# Checks what each RNA layer contains, how the gene network is oriented, cell counts by
# type/stage/sample, and boundary support on the shipped graph. Writes to results/p1/inputs/.

# %%
import argparse
import json
import sys
from pathlib import Path

try:
    _here = Path(__file__).resolve().parent
except NameError:
    _here = Path.cwd()
REPO = next(p for p in [_here, *_here.parents] if (p / "config" / "experiment.yaml").exists())
sys.path.insert(0, str(REPO / "src"))

from p1_data.common import load_config, load_source, repo_path
from p1_data.inputs import boundary_support, count_tables, grn_checks, provenance_checks

CONFIG = None

ap = argparse.ArgumentParser(allow_abbrev=False)
ap.add_argument("--config", default=CONFIG)
args, _ = ap.parse_known_args()

# %%
# Load config and data.
cfg = load_config(args.config)
d = cfg["data"]
out = repo_path(cfg["output"]["results_dir"]) / "inputs"
out.mkdir(parents=True, exist_ok=True)
adata = load_source(cfg)
print(adata)

# %%
# Layer provenance and gene-network checks.
report = {
    "shape": list(adata.shape),
    "layers": [str(k) for k in adata.layers.keys()],
    "obs_columns": list(adata.obs.columns),
    "provenance": provenance_checks(adata),
    "grn": grn_checks(adata, cfg["endpoints"]),
}
print(json.dumps(report, indent=2, default=str))

# %%
# Count tables and boundary support on the shipped graph.
for name, tab in count_tables(adata, d["celltype_key"], d["stage_key"], d["sample_key"]).items():
    tab.to_csv(out / f"{name}.csv")
bs = boundary_support(adata.obsp["connectivities"], adata.obs[d["celltype_key"]], cfg["transitions"])
bs.to_csv(out / "boundary_support_shipped_graph.csv", index=False)
print(bs.to_string(index=False))

# %%
# Save the report.
with open(out / "inputs_report.json", "w") as fh:
    json.dump(report, fh, indent=2, default=str)
print(f"Wrote {out}")

# %%
