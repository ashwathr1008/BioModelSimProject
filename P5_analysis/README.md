# P5 analysis — final, report/presentation-ready package

**Start here: [PROJECT_FINAL_GUIDE.md](PROJECT_FINAL_GUIDE.md).** This is the team's final scientific handoff: executed P1→P5 workflow, quantitative results, methodological limits, source-file links, PPT slide plan and report-ready English text.

This directory was updated from the user's **locally completed** P5 outputs on 2026-11-25. The analyzed results are **frozen**, not just theoretical code or GitHub-only preliminary summaries.

## Authoritative files

| Purpose | Location |
|---|---|
| **Whole-project scientific guide** | [`PROJECT_FINAL_GUIDE.md`](PROJECT_FINAL_GUIDE.md) |
| Final local P5 results: 18 tables | [`outputs/tables/`](outputs/tables/) |
| Final local P5 graphs: 3 main figures | [`outputs/figures/`](outputs/figures/) |
| Local run summary | [`outputs/reports/findings.md`](outputs/reports/findings.md) |
| Complete evidence/analysis status | [`STATUS.md`](STATUS.md) |
| Frozen SHA-256 evidence snapshot | [`FINAL_FREEZE_MANIFEST.json`](FINAL_FREEZE_MANIFEST.json) |
| Freeze verifier | [`verify_final_package.py`](verify_final_package.py) |
| Core reproducible analysis | [`run_analysis.py`](run_analysis.py) |
| Supplemental earlier independent source checks | [`outputs/verified_from_github/`](outputs/verified_from_github/) |

**Source-of-truth rule:** report numeric P5 results from `outputs/tables/` and the original P1–P4 source paths cited in `PROJECT_FINAL_GUIDE.md`. Some historical `.md` files under `outputs/verified_from_github/` and older P1/P2 handoff docs still describe a *pre-completion* state; do **not** use their outdated status statements.

## Using this after unzipping into the team repo

Unzip with the `P5_analysis/` folder at repository root (alongside `config/`, `results/`, and `scripts/`). The ZIP does **not** add, rewrite or delete P1–P4 source files.

```bash
cd BioModelSimProject-main
python3 P5_analysis/verify_final_package.py
python3 -m pytest -q P5_analysis/tests
# Optional: rerun the locally completed P5 metrics if the P3/P4 inputs are available
python3 P5_analysis/run_analysis.py
```

Running `run_analysis.py` **recomputes** P5 outputs; if you want to preserve the exact frozen version for Git, verify and commit it without rerunning. The checksums will flag any changed file. Large P2 fitted models remain in the external [P2 release](https://github.com/ashwathr1008/BioModelSimProject/releases/tag/p2-x-completed-2026-11-22). The nine 10×10 transition matrices were **not run** and are **not part of the final Results**; the optional `--matrix` code remains for future work.

**Reproducibility caveat:** Only three computational fit blocks were run. No new stage-matched models, biological replicates, experimental TF knockouts, or GRN parameter-recovery experiments are implied by this package.
