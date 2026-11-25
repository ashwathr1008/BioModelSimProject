# FINAL P5 execution and evidence status — 2026-11-25

**Read [PROJECT_FINAL_GUIDE.md](PROJECT_FINAL_GUIDE.md) for the team's report/PPT reference.** This status supersedes historical pre-local-run status notes contained in early versions of this package. It does not rewrite the team's P1–P4 files.

| Deliverable | Final state | Canonical evidence |
|---|---|---|
| P1 sampled inputs and graph QC | Completed by P1 | `../docs/p1_data_notes.md`, `../results/p1/qc/` |
| P2 nine experimental X fits + tutorial | Completed by P2; large binaries separately released | `../docs/p2_model_fitting.md`, `../results/p2/comparison/run_summary.csv` |
| P3 fate inference, 106 anchors, CBC and timing | Completed by P3 | `../results/p3/fate/`, `../results/p3/cell_level_cbc/`, `../results/p3/terminal_anchors/` |
| P4 separate MSC/MYOG blockade, 18 main analyses | Completed by P4 | `../docs/p4_results.md`, `../results/p4/qc/main_outputs.json`, `../results/p4/run_manifest.csv` |
| P5 baseline fate TVD and all shared-cell comparisons | **COMPUTED on user's local run** | `outputs/tables/baseline_fate_by_repeat.csv`, `baseline_fate_per_cell.csv.gz` |
| P5 per-cell CBC / coverage, celltype, stage and non-MyoB1 checks | **COMPUTED** | `outputs/tables/cbc_matched_boundary_cells.csv`, `baseline_fate_by_cell_type.csv`, `baseline_fate_by_stage.csv`, `fate_composition_robustness.csv` |
| P5 Full fit seed sensitivity | **COMPUTED** | `outputs/tables/full_training_seed_variability.csv` |
| P5 P4 effect comparison and cross-layer marginal error | **COMPUTED** | `outputs/tables/perturbation_effects_existing_p4.csv`, `cross_layer_baseline_vs_perturbation.csv` |
| P5 **direct per-cell** P3/P4 cross-layer association | **COMPUTED from 18 P4 .csv.gz inputs in user's local run** | `outputs/tables/cross_layer_per_cell_association.csv`, `cross_layer_by_celltype.csv` |
| P5 short result report and primary figures | **COMPUTED** | `outputs/reports/findings.md`, `outputs/figures/` |
| 10×10 transition matrices across nine fits | **NOT EXECUTED**; P2 fitted H5AD/CellRank are required | Optional code in `run_analysis.py`, `transition_aggregation.py` |
| Stage-matched random controls | **NOT FITTED**; no new training budget | P1 stage imbalance remains a design limitation |
| Full learned-GRN weight stability | **NOT EXECUTED** | Requires P2 model checkpoints; optional future work |

**Final findings:** larger mean baseline fate TVD after Targeted than Random in all 3 repeats, reduction in MyoB1-associated absolute CBC coverage, larger paired CBC shifts along both MyoB1-related transitions, persistent major MSC/MYOG signed mean effects but no consistent Targeted excess TF-response error. Meaningful training-seed variation and developmental-stage imbalance constrain causal attribution. Numbers, tables, graphical assets, and the exact allowed vs unsupported statements are indexed in [the guide](PROJECT_FINAL_GUIDE.md).

**Immutable snapshot:** `FINAL_FREEZE_MANIFEST.json` checksums the bundled scientific data/plots/code/docs. `python3 P5_analysis/verify_final_package.py` verifies those bytes. This does *not* certify missing external model binaries or the biological correctness of the inference.
