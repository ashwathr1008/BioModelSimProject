# BioModelSimProject: RegVelo under incomplete developmental sampling

02-512/02-712 final project. Plan: `regvelo_project_plan.md` (audited version:
`regvelo_project_plan_audited.md`). Roles: `regvelo_team_work_allocation.md`.

## Layout

```
config/experiment.yaml   shared frozen settings (sampling, preprocessing, seeds, endpoints)
src/p1_data/             Person 1: input checks, sampling manifests, per-condition rebuild
scripts/p1_0*.py         Person 1 entry points (run in order)
docs/p1_data_notes.md    Person 1 findings, decisions and the hand-off contract
results/p1/              manifests, QC tables and figures (committed; small)
data/                    raw + prepared AnnData (gitignored; regenerate with the scripts)
```

## Quick start (Person 1 pipeline, about 2 minutes on a CPU)

```bash
pip install -r requirements-p1.txt
python scripts/p1_01_inspect_inputs.py
python scripts/p1_02_make_manifests.py
python scripts/p1_03_build_inputs.py
python scripts/p1_04_compare_conditions.py   # before/after comparison plots
```

Every `.py` file is split into `# %%` cells, so you can also open a script in VS Code and use
**Run Cell** (Interactive Window). Options such as `ONLY` or `PCA_INPUT` sit at the top of each script.
