# Person 3 — RegVelo Analysis

## 1. Scope and evidence

Person 3 analyzes the fitted RegVelo models for the human limb development project: terminal anchor definition, cell fate inference, developmental timing, cross-boundary correctness (CBC), and visualization. Person 4 consumes the shared fate-inference interface for perturbation work. Person 5 performs aligned cross-condition comparisons, statistical analysis, and integrated figures.

This document is based on the retrievable text of [介紹Project與Person3工作](chatgpt-conversation://6b050308-72b0-83ea-8988-733a8ec2c3e1), read on 2026-11-24. The retrieval returned five recent exchanges and no older-page cursor. It includes complete proposed plotting scripts, user-pasted tutorial execution output, and assistant summaries of earlier analysis. Earlier core scripts, the nine-run numerical tables, and the contents of attached images were not available in that text.

Evidence is distinguished throughout:

- **Recorded output:** numerical output pasted by the user in the source conversation.
- **Conversation report:** a statement or interpretation in an earlier assistant response, not independently verified here.
- **Script behavior:** behavior visible in the supplied plotting code; this does not establish that the code was saved or executed.
- **Pending verification:** information absent from the retrieved source.

Only this Markdown document was created for this handoff. No analysis was rerun, no existing project artifacts were verified on disk, and nothing was pushed to GitHub. Paths below describe the source conversation's project layout, not verified files in the current workspace.

## 2. Experimental design and inputs

The described comparison comprises three conditions and three computational repeats, for nine fitted models:

| Condition | Description in the conversation | Run names |
|---|---|---|
| Full | Complete experimental dataset | `full_r0`, `full_r1`, `full_r2` |
| Targeted | Removal of 75% of MyoB1 cells | `targeted_r0`, `targeted_r1`, `targeted_r2` |
| Random | Random-removal comparison | `random_r0`, `random_r1`, `random_r2` |

The retrieved text does not specify the random-removal count, sampling eligibility, seeds, or whether repeat seeds pair conditions. These must be recovered from the preprocessing/model-fitting records before describing the control as exactly size matched.

The plotting scripts read each fitted dataset from:

```text
results/p2/X/{condition}_r{repeat}/fitted.h5ad
```

Relevant fields visible in the scripts are `obs_names` (cell identifiers), `obs['leiden_R_celltype']` (cell types), and `obsm['X_umap_ref']` (reference visualization coordinates).

## 3. Fixed terminal anchors

The conversation specifies 106 fixed terminal anchors shared across analyses:

| Terminal lineage | Anchors | Missing in tutorial validation |
|---|---:|---:|
| `PAX7+MyoProg1` | 53 | 0 |
| `MYH3+MyoC` | 23 | 0 |
| `MYL3+MyoC` | 30 | 0 |
| **Total** | **106** | **0** |

Counts and tutorial presence are supported by recorded output. The conversation describes the anchors as designated absorbing states and names `terminal_anchors.json` as a handoff artifact. Its directory, JSON schema, cell identifiers, and exact selection algorithm were not included in the retrieved text. Selection thresholds, stage filters, tie breaking, and seed settings must therefore not be inferred.

Fixed anchors provide a common endpoint definition, but unequal counts (53/23/30) remain a methodological limitation. Presence in the tutorial dataset does not prove presence, uniqueness, correct labeling, or exclusion from removal in all nine experimental datasets. Those checks still require the actual anchor file and datasets.

## 4. Cell fate inference

The conversation reports CellRank fate inference using the weighted transition kernel:

```text
K = 0.8 × VelocityKernel + 0.2 × ConnectivityKernel
```

The resulting fate probabilities describe predicted absorption into the three specified terminal lineages. They are conditional on the fitted model, graph, kernel, and selected endpoints; they are not directly observed lineage outcomes.

`fate_inference.py` is named as the shared implementation for Persons 3 and 4. Its actual source was absent from the retrieved exchanges. The callable name, argument signature, return types, estimator class, solver, convergence settings, velocity representation, and package versions cannot be documented reliably from this record. No runnable API example is invented here.

The plotting code establishes the following CSV consumer contract:

```text
results/p3/fate/{run}_fate.csv
```

- A `cell_id` column is read as the index.
- Required probability columns are `PAX7+MyoProg1`, `MYH3+MyoC`, and `MYL3+MyoC`.
- Fate-table cell IDs must be unique.
- Their set must equal the corresponding fitted dataset's `obs_names` set.
- Rows are explicitly reordered with `fate.loc[adata.obs_names]` before plotting.

These scripts do not explicitly check finite probabilities, probability bounds, row sums, or uniqueness of `adata.obs_names`. Such validation should be confirmed in the core inference/output validation code before handoff.

## 5. Trajectory measurements

### Developmental timing

The analysis reports Spearman correlation for latent time and developmental timing. The tutorial dataset contains `obs['latent_time']` and `obs['adj_stage_num']`, but the correlation implementation is not in the retrieved record. The exact time vector, stage mapping, missing-value handling, and aggregation rules require confirmation from the original script.

The tutorial output gives:

```text
spearman_rho = 0.6518723654863241
p_value      = 0.0
n_cells      = 12207
```

The printed p-value should be treated as the numerical output of the computation, not proof of an exactly zero probability. Cell-level correlation also does not make individual cells independent biological replicates.

### Cross-boundary correctness (CBC)

The four reported transitions are:

| Pathway | Source | Target |
|---|---|---|
| Targeted | `PAX3+MyoProg` | `MyoB1` |
| Targeted | `MyoB1` | `MyoC1` |
| Control | `PAX3+PAX7+MyoProg` | `MyoB2` |
| Control | `MyoB2` | `MyoC2` |

CBC is used to assess local directional consistency across these boundaries. The available output contains mean, median, standard deviation, and valid/invalid/source/target counts. The exact formula, neighbor selection, coordinate and velocity projection, weighting, and invalid-score criteria are absent. In particular, `n_valid` should not be treated as the total source population: the reported counts differ substantially.

Cell-level CBC was reported complete and `all_cell_cbc.csv` was named for handoff. Its schema and values were not supplied, so paired-cell analyses cannot be reconstructed from this document alone.

## 6. Tutorial baseline validation: recorded results

The user-pasted log identifies the historical input as:

```text
C:\Users\willy\Documents\master3\StatiscalMethod\project\BioModelSimProject\results\p2\tutorial\fitted.h5ad
```

This is a historical log path, not a claim that the file exists locally now.

The recorded AnnData shape is **12,207 cells × 606 genes**. Available fields include `latent_time`, `adj_stage_num`, `leiden_R_celltype`, `X_pca`, `X_umap`, `velocity`, `Ms`, `Mu`, `fit_t`, `latent_time_velovi`, `connectivities`, and `distances`. The tutorial field listing does not include `X_umap_ref`; the experimental plotting scripts should therefore not be assumed to work on this tutorial file unchanged.

All 106 anchors were found. Kernel building and fate inference completed in the pasted log, which ended with `Tutorial Validation Complete` and a reported save directory of `results/p3/tutorial`.

### Fate probabilities across tutorial cells

| Statistic | PAX7+MyoProg1 | MYH3+MyoC | MYL3+MyoC |
|---|---:|---:|---:|
| Count | 12207 | 12207 | 12207 |
| Mean | 0.708070 | 0.264881 | 0.027049 |
| Standard deviation | 0.365430 | 0.333574 | 0.098620 |
| Minimum | 0.000000 | 0.000000 | 0.000000 |
| 25th percentile | 0.457765 | 0.003435 | 0.000161 |
| Median | 0.926481 | 0.067792 | 0.003001 |
| 75th percentile | 0.996361 | 0.466813 | 0.020222 |
| Maximum | 1.000000 | 1.000000 | 1.000000 |

These marginal summaries do not establish per-cell probability normalization or agreement with another model.

### Tutorial CBC

| Source → target | Mean | Median | SD | Valid | Invalid | Source cells | Target cells |
|---|---:|---:|---:|---:|---:|---:|---:|
| PAX3+MyoProg → MyoB1 | 0.773615 | 0.790816 | 0.146638 | 121 | 0 | 583 | 759 |
| MyoB1 → MyoC1 | 0.815511 | 0.840730 | 0.107139 | 77 | 0 | 759 | 412 |
| PAX3+PAX7+MyoProg → MyoB2 | 0.676316 | 0.680150 | 0.154201 | 399 | 0 | 1790 | 1527 |
| MyoB2 → MyoC2 | 0.767185 | 0.805685 | 0.136250 | 283 | 0 | 1527 | 639 |

The earlier assistant described a tutorial reference Spearman range of 0.63–0.68. The observed 0.651872 lies within that stated range, but no original tutorial citation or external validation of the range was available. This supports agreement with the conversation's stated benchmark, not a separately verified reproduction of every tutorial method. The conversation explicitly notes that the project used its shared anchors rather than the original tutorial anchors.

## 7. Experimental findings and their evidence limits

The conversation reports that fate inference, timing, CBC, and cell-level CBC were completed for all nine models. Complete nine-run numerical results were not returned by the thread reader. No nine-run results table or cross-repeat effect estimate can therefore be supplied faithfully here.

One earlier assistant response provides these rounded whole-dataset means for `full_r0`:

| Lineage | Tutorial mean, recorded output | Full r0 mean, conversation report |
|---|---:|---:|
| PAX7+MyoProg1 | 70.8070% | 73.93% |
| MYH3+MyoC | 26.4881% | 24.05% |
| MYL3+MyoC | 2.7049% | 2.02% |

The Full r0 figures are secondary reports, not values re-read from a CSV. Similar population means do not demonstrate agreement for matched individual cells.

The earlier assistant summarized the experimental comparisons as follows:

- Global developmental timing did not consistently decline.
- Directional consistency into MyoB1 decreased slightly, while the MyoB2 pathway was relatively stable.
- In the Repeat 0 comparison image, PAX7 fate patterns looked broadly similar; MYH3 fate showed more noticeable local differences and apparently lower probabilities in parts of the targeted condition; the MYL3 high-probability region remained spatially localized with apparent probability changes.

These are conversation-level qualitative interpretations. The images were not independently inspected for this document, and numerical support for the timing/CBC comparisons was not present in the retrieved record. They must not be presented as confirmed effect sizes or statistically significant results.

The proposed interpretation—that local trajectory measures may be relatively robust while long-term fate predictions change—is a hypothesis for Person 5 to test. Repeat 0 images alone, with different retained cell sets across conditions, cannot establish it.

## 8. Visualization scripts and reproducible commands

Two complete plotting scripts were supplied in the source conversation. Their behavior can be documented, but their existence and execution in a repository were not verified here.

### `scripts/p3_07_visualize_fate.py`

The supplied script reads `full_r0` fitted data and its fate CSV, checks fate IDs for uniqueness and set equality, aligns rows by cell ID, and plots using `X_umap_ref`.

```powershell
# Run from the original project root after confirming the script and inputs exist.
python scripts/p3_07_visualize_fate.py
```

Declared outputs:

```text
results/p3/figures/full_r0_cell_types.png
results/p3/figures/full_r0_fate_umap.png
```

The cell-type plot uses `leiden_R_celltype`. Fate panels use `viridis`, a fixed 0–1 color scale, equal aspect ratio, and high-probability cells drawn last. Images are saved at 300 dpi. The conversation subsequently contains image attachments and assistant interpretations of both plots, but their files were not checked here.

### `scripts/p3_08_compare_fate.py`

The supplied script loops over repeats 0, 1, and 2 and conditions `full`, `targeted`, and `random`. It performs the same fate-ID checks and row alignment for each dataset.

```powershell
# Run from the original project root after confirming all nine input pairs exist.
python scripts/p3_08_compare_fate.py
```

Each figure has condition rows and terminal-lineage columns. All panels use the Full dataset's UMAP axis limits for that repeat, expanded by one unit, and a 0–1 color scale. The script retains each condition's entire available cell set and saves:

```text
results/p3/figures/fate_comparison_r0.png
results/p3/figures/fate_comparison_r1.png
results/p3/figures/fate_comparison_r2.png
```

The conversation discusses the Repeat 0 image. Generation of Repeat 1 and Repeat 2 images was explicitly conditional in the latest summary and remains unconfirmed.

Both plotting scripts only read existing fitted datasets and fate tables; they do not rebuild kernels or recompute fate probabilities. Reference UMAP is a visualization coordinate system, not evidence that inferred trajectories are unchanged. The scripts use the supplied reference coordinates but do not verify that common cells have identical coordinates across datasets.

### Reproduction boundary

The source supports the two commands above. It does not expose executable commands or filenames for anchor selection, core inference, nine-run aggregation, cell-level CBC, or tutorial validation. An end-to-end reproducible workflow requires those original scripts, the anchor file, input data, package versions, seed settings, and actual function signatures. These are missing provenance items, not instructions to substitute guessed implementations.

## 9. Handoff interfaces

### Person 4: perturbation analysis

The conversation names these deliverables:

- `fate_inference.py`: shared fate-inference implementation; path and callable contract require inspection.
- `terminal_anchors.json`: fixed endpoint IDs; path and schema require inspection.
- Nine unperturbed fate CSVs, consumed using the cell-ID and lineage-column contract in Section 4.

Before perturbation comparisons, Person 4 should confirm the shared implementation, preserve the agreed endpoint definitions and 0.8/0.2 kernel settings, and verify anchor availability in each input. The interface must be read from the actual module rather than reconstructed from this document. Perturbed and unperturbed predictions must be aligned by stable cell IDs and explicit lineage names.

### Person 5: comparisons and statistics

The conversation names `all_cell_cbc.csv`, `all_cbc.csv`, `all_timing.csv`, and nine fate CSVs. The aggregate CSV directories and exact column schemas were not shown. Only the fate CSV location and consumer schema are supported by the plotting code.

Person 5 should:

1. Recover run/condition/repeat identifiers and verify completeness against the nine expected runs.
2. Align fate comparisons on common retained cell IDs and use the agreed non-terminal population. Define whether “non-terminal” excludes only anchors or all endpoint-labeled populations; the retrieved text does not resolve this distinction.
3. Compare cell-level CBC only after checking boundary-cell eligibility and membership across conditions; report changes in valid cell counts alongside scores.
4. Preserve the distinction between computational repeats and independent biological replication when choosing statistical procedures.
5. Quantify effects across repeats before interpreting the qualitative Repeat 0 observations as robust findings.

### Artifact status ledger

| Artifact | Evidence available here |
|---|---|
| `terminal_anchors.json` | Named in handoff; tutorial log verifies counts/presence, not JSON contents |
| `fate_inference.py` | Named shared module; source/signature unavailable |
| Nine fitted experimental datasets | Paths referenced by plotting code; not locally verified |
| Nine fate CSVs | Reported complete; paths/schema visible in plotting code; values unavailable |
| Aggregate timing/CBC and cell-level CBC CSVs | Named and reported complete; files/schema unavailable |
| Tutorial validation | User-pasted execution output with numerical results and reported save directory |
| Single-model plots and Repeat 0 comparison | Attachments discussed in source conversation; not inspected here |
| Repeat 1 and Repeat 2 comparison plots | Expected outputs; completion unconfirmed |
| `p3_analysis.md` | Created in the current workspace for this request |

## 10. Limitations and outstanding validation

- **Unequal anchor counts:** the 53/23/30 allocation may affect endpoint comparisons; no anchor-count sensitivity analysis is documented in the retrieved record.
- **Preprocessing differences:** the conversation mentions PCA, neighbor graph, moments, and fitting differences between tutorial and reconstructed experimental inputs. A timing difference cannot be assigned solely to PCA; exact settings and causal attribution are unresolved.
- **Changing cell composition:** full-condition and removal-condition distributions reflect different cells. Whole-dataset fate means and plots are not paired comparisons.
- **Changing boundary composition:** CBC changes may reflect altered boundary-cell eligibility and neighborhood structure as well as altered velocities.
- **Computational repeats:** three repeats do not establish independent biological replication or justify treating all cells as independent experimental units.
- **Imposed endpoints:** high fate probabilities near absorbing anchors are expected by construction and do not independently validate biological correctness.
- **Tutorial comparability:** shared project anchors differ from the original tutorial anchors; the stated reference range lacks its original citation in the retrieved material.
- **Incomplete numerical provenance:** the nine-run results, core implementations, software versions, and seeds are unavailable in the retrieved source. Stronger quantitative conclusions and complete reproduction instructions require those records.
- **Unverified local artifacts:** no historical data, source code, aggregate CSV, or figure file was opened or validated during this documentation task. Reported save messages establish what the earlier run printed, not present-day file availability.

Before declaring the analysis handoff fully validated, inspect the original core scripts and schemas, verify all nine outputs and anchor coverage, check fate finiteness/bounds/normalization and cell-ID integrity, and confirm the remaining figures. These are outstanding verification tasks; this document does not claim they have been completed.

## 11. Source traceability

All project-specific claims above come from the linked conversation, with the following evidence anchors:

| Source exchange | Material used |
|---|---|
| User message beginning `===== Tutorial Dataset =====` | Dataset dimensions/fields, anchor coverage, exact timing result, fate summary, CBC table, reported completion |
| Assistant response to tutorial output | `p3_07_visualize_fate.py`, tutorial-range claim, rounded Full r0 means, preprocessing caveats |
| Assistant response to Full r0 fate image | Qualitative single-model interpretation and reported work status |
| Assistant response to cell-type image | Complete `p3_08_compare_fate.py`, expected comparison outputs, common-cell caution |
| Assistant response to Repeat 0 comparison image | Qualitative condition comparisons, 106-anchor and 0.8/0.2 settings, Person 4/5 handoff inventory |
| Final documentation-planning response | Intended English document and historical destination `BioModelSimProject/docs/p3_analysis.md`; explicit lack of confirmed document generation/GitHub upload at that point |

No external literature validation or new model computation was performed for this documentation-only task.
