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
