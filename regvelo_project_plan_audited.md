# RegVelo under incomplete developmental sampling: audited project plan

Course: 02-512 / 02-712 Computational Methods for Biological Modeling and Simulation

Audit date: November 15, 2026. Scope: a two-week computational project with AI assistance.

This is a corrected companion to `regvelo_project_plan_Opus.md`; the original is unchanged. The audit checks the plan against the supplied final Cell paper, official tutorials, current documented APIs, and available source code. No dataset was downloaded, model trained, or end-to-end implementation tested for this audit. API details must be checked against the environment actually installed.

## 1. Verdict and changes

Keep the central experiment. It has a useful question, a documented dataset, and a manageable scope: one intermediate population, three conditions, three computational repeats, and two selected TF perturbations. The experiment needs tighter controls and more careful interpretation, not a larger model or more datasets.

| Item in the Opus plan | Audit | Revision |
|---|---|---|
| Human limb tutorial, supplied genes/prior, MSC and MYOG | Good starting point | Keep. This measures the stability of two known-regulator predictions, not genome-wide driver discovery. |
| Three sampling conditions and repeated fits | Appropriate scope | Keep nine main fits and eighteen separate perturbation analyses. |
| Fixed model variant and regularization | Good control | Keep `soft_regularized`, `lam2=1.0`; fix other settings too. |
| Rebuild moments after thinning | Essential, but sketch is incomplete | Explicitly rebuild neighbors. Updating PCA does not guarantee that `moments` replaces an inherited graph. |
| Raw layers inferred from their names | Not established | Verify that RNA layers are unsmoothed per-cell measurements and document their transformations. Normalized unsmoothed layers may be usable. |
| Full-data and reduced-data metrics computed on different populations | Can confound changed predictions with changed composition | Use the same surviving evaluation cells within each repeat for primary comparisons. |
| Random removal from the whole dataset | Also changes endpoint availability | In the revised primary design, protect all three endpoint populations and randomly remove eligible nonterminal cells. State this narrower comparator explicitly. |
| Fixed terminal names | Names do not fix absorbing cell identities | Select and save terminal anchor cell IDs once; reuse them for every condition and perturbation. |
| Depletion likelihood means lineage shrinkage | Incorrect | It is a rank-based summary of modeled fate-probability shifts, not predicted cell loss. |
| Lost statistical significance means a weakened regulator effect | Incorrect as a standalone criterion | Compare effect size, direction, and consistency across repeats. |
| MyoB2 is an unaffected internal control | Too strong | It is a comparison branch sharing the fitted GRN and graph; indirect changes are possible. |
| Match tutorial numbers exactly | Unnecessarily strict | Check qualitative behavior and investigate large discrepancies; do not tune toward exact values. |
| Eight-to-nine macrostate warning shows instability | Unsupported | This reflects a real-Schur complex-conjugate-block constraint, not evidence of biological instability. |
| Three seed pairs separate sampling and training noise | They do not | Treat them as three computational blocks with combined sources of variability. |
| No difference proves robustness | Too strong | Limit conclusions to the tested thinning scheme, fixed inputs, and observed variation. |

## 2. Question and estimand

**With the published gene panel and regulatory prior held fixed, does selectively reducing one intermediate myoblast population change RegVelo's trajectory and MSC/MYOG perturbation predictions more than removing an equal number of randomly chosen nonterminal cells?**

This is a controlled sensitivity study of the full analysis pipeline. Each reduced dataset gets its own preprocessing, fitted model, and transition graph. Differences can therefore arise through altered smoothing, scaling, neighborhood structure, or fitted dynamics. Do not attribute every change exclusively to an unidentified regulatory edge.

The complete-data fit is a reference for stability, not biological ground truth. Known differentiation transitions and developmental-stage annotations provide additional, imperfect biological checks. This is not a held-out-cell forecasting study.

## 3. Data and fixed choices

Use `rgv.datasets.human_limb()`. The documented object contains 12,207 cells and 606 retained genes, myogenic cell labels in `leiden_R_celltype`, stage labels in `adj_stage`, sample labels in `adj_sample`, TF flags, and a supplied regulatory skeleton. The tutorial covers post-conception weeks 5-9. Do not assume that `adj_sample` necessarily means an independent donor without checking metadata.

Fix before looking at targeted-versus-random outcomes:

- The ordered 606-gene panel, TF definitions, and supplied GRN, including verified matrix orientation and gene alignment.
- The selected model configuration, other regularization settings, optimizer, stopping rule, maximum training budget, and output sampling settings.
- The normalization policy, PCA gene mask and dimensionality, neighbor count, and preprocessing random seed.
- Fate inference using the tutorial's `0.8 * VelocityKernel + 0.2 * ConnectivityKernel` mixture, with the same construction in unperturbed and perturbed analyses.
- The three reference endpoints: `PAX7+MyoProg1`, `MYH3+MyoC`, and `MYL3+MyoC`.
- The target population, thinning fraction, sampling protocol, evaluation cell definitions, and metric aggregation.

`PAX7+MyoProg1` is a progenitor endpoint/reservoir in this analysis; calling it a terminal state does not mean that it is mature muscle or an irreversible biological endpoint.

The fixed genes and prior are intentional background information. The conclusion is conditional on their availability. Rebuilding the graph removes deleted cells' direct contribution through neighborhood averaging; it does not erase all information contributed by the original dataset to the published gene selection or prior.

## 4. Pilot and baseline acceptance

### A. Tutorial smoke test

Run one fit using the published selected configuration, rather than first reproducing the entire three-mode by three-repeat model-selection grid. Confirm that the data load, training finishes, fate probabilities are valid, and separate MSC/MYOG perturbations run. Measure the total time for training, output generation, CellRank, and both perturbations.

The tutorial reports stage correlations around 0.63-0.68 for its regularized fits and depletion likelihoods around 0.67 for MSC toward the progenitor endpoint and 0.78/0.75 for MYOG toward the two mature endpoints. Treat these as references, not acceptance thresholds. Check broadly consistent directions and plausible trajectories; investigate major discrepancies rather than tuning until numbers match.

### B. Rebuilt experimental baseline

Process the full dataset through the same clean rebuilding pipeline used for both reduced conditions. This full-data version is the baseline for the actual experiment. It can differ from the shipped object, for example if the original graph used more genes before the final 606-gene selection. Keep the shipped-data smoke test and rebuilt-data comparison conceptually separate.

Document and freeze the working package versions. Current docs and older tutorial cells do not always use identical interfaces. A notebook that imports successfully is not yet a completed end-to-end smoke test.

## 5. Sampling conditions and repeats

| Condition | Removal rule |
|---|---|
| Full | Retain all cells. |
| Targeted | Remove 75% of MyoB1 cells, subject to the count/connectivity check below. |
| Matched random nonterminal | Remove exactly the same number of cells randomly from all eligible nonterminal cells, including MyoB1. Protect all cells bearing the three endpoint labels. |

Protect the same endpoint populations in both reduced conditions. This is an explicit change from global random removal: it isolates sensitivity to sampling the nonterminal part of the trajectory while retaining endpoint availability.

MyoB1 is justified by the annotated route `PAX3+MyoProg -> MyoB1 -> MyoC1`. The comparison route is `PAX3+PAX7+MyoProg -> MyoB2 -> MyoC2`. MyoB2 is not guaranteed unaffected: both routes share the model and regulatory network.

Before fixing the removal rule, tabulate cell type by stage/sample. Check target-population size, retained counts, and boundary support. Choose MyoB1 unless the count/connectivity pilot shows the proposed thinning makes its transitions unassessable; then reduce the fraction or select MyoB2 under a recorded, count-based rule. Do not choose the population or fraction because it produces a larger regulatory effect. Freeze this choice before the main comparison.

If feasible, sample within observed sample/stage strata so targeted removal does not accidentally eliminate a specimen. A random control matched to the targeted deletions' sample/stage counts is an optional refinement without extra model fits. If using simple random sampling, report the resulting stage/sample composition and avoid attributing the result solely to cell identity when specimen coverage also changes.

Use three repeat blocks. In block i, create independent targeted and random deletion sets using recorded sampling seeds; use the same training-seed label for full, targeted, and random fits. Keep sampling and training seeds in separate recorded fields. This is blocking, not exact cancellation of stochastic training: different cell sets create different minibatches. With one fit per subset, sampling and optimization variation remain combined.

Save retained and removed cell IDs, cell counts by type/stage/sample, seeds, and settings. The main budget is nine trained models and eighteen TF perturbation analyses, plus any nonreusable pilot work.

## 6. Correct rebuilding contract

The original subsetting sketch should not be used unchanged. A sliced AnnData retains graph and other derived state. Recomputing PCA alone does not guarantee replacement of the old neighbor graph when `scv.pp.moments` runs.

Use this sequence:

1. Verify provenance and transformations of `X`, `spliced`, and `unspliced`. Layer names alone do not prove raw counts. Avoid repeated normalization or log transformation. RNA layers must be valid unsmoothed per-cell inputs.
2. Choose retained cell IDs and create a clean AnnData carrying only the intended per-cell expression inputs, annotations, ordered genes, TF flags, and aligned GRN. Do not inherit fitted-model registration, old moments, velocities, transition matrices, or fate probabilities.
3. Apply the declared normalization policy consistently. Never silently reuse old smoothing.
4. Recompute PCA and explicitly call the neighbors routine on the retained cells.
5. Compute `Ms`/`Mu` from the newly built graph and verified RNA layers.
6. Apply RegVelo's declared scaling without reselecting the genes. Its documented preprocessing helper supports `filter_on_r2=False` when preserving the fixed panel. Do not apply scaling twice.
7. Check exact gene order, prior dimensions/orientation, TF presence, finite values, constant genes, graph connectivity, and retained-cell counts before fitting.

Source-informed illustration of the essential graph step, not an execution-tested replacement notebook:

```python
# Preconditions: clean sub object, verified normalization, fixed PCA gene mask.
sc.pp.pca(sub, n_comps=30, svd_solver="arpack", random_state=0)
sc.pp.neighbors(
    sub, n_pcs=30, n_neighbors=chosen_k,
    use_rep="X_pca", random_state=0,
)
scv.pp.moments(sub, n_pcs=None, n_neighbors=None)
sub = rgv.pp.preprocess_data(sub, filter_on_r2=False)
```

Check these signatures against the installed versions and set the PCA gene mask explicitly. Choose `chosen_k` once during the pilot rather than assume 30 reproduces the original data preparation. If only smoothed RNA values are available, obtain valid unsmoothed inputs or stop this thinning analysis; documenting leakage does not make its conclusions valid.

## 7. Terminal anchors and separate TF simulations

Once, using the unperturbed rebuilt full-data reference, identify suitable terminal anchor cells belonging to the three endpoint labels. Save a dictionary from endpoint name to exact `obs_names`. A tutorial-style macrostate procedure can be used for this initial selection. Review label agreement and keep the anchor selection fixed thereafter; never choose anchors using perturbed outcomes.

Because endpoint populations are protected, these anchors survive in every condition. For each fit:

1. Build the prescribed mixed kernel and call `GPCCA.set_terminal_states(anchor_dict)` before computing unperturbed fate probabilities.
2. Simulate MSC and MYOG independently, each starting from the saved unperturbed model. Do not apply MYOG to an already MSC-perturbed model.
3. Recompute perturbed velocity kernels and fate probabilities using the same condition-specific connectivity kernel, the same blend, and the same anchor IDs.
4. Align fate-probability columns by endpoint name and rows by cell ID.

The current `cellfate_perturbation` API does accept `terminal_indices`; that part of the original plan is valid. However, do not rely on a summary-function argument to correct fate matrices already computed with different boundaries. Control anchors explicitly upstream and verify any helper's behavior in the installed version.

The tutorial's eight-to-nine macrostate adjustment accommodates complex conjugate eigenvalue blocks. It is not, by itself, evidence of unstable endpoint detection. After anchors are fixed, endpoint discovery is not an outcome of this experiment; do not claim that specified endpoints were independently recovered.

## 8. Evaluation that compares like with like

For each repeat i, define `E_i` as the cell IDs shared by full, targeted, and random conditions. Use the same ordered `E_i` to score all three fits. For fate and TF-effect summaries, use shared nonterminal cells, and summarize prespecified source groups separately where possible. This avoids diluting the analysis with absorbing anchors and prevents changing cell mixtures from masquerading as changing predictions. These cells participated in fitting; this is a stability comparison, not out-of-sample validation.

### Primary readout: changes in predicted TF effects

For each cell c, endpoint k, TF, and sampling condition j, calculate:

`delta_p[c,k,TF,j] = p_perturbed[c,k,TF,j] - p_unperturbed[c,k,j]`.

Report mean/median signed effects within the same shared source-cell groups. Also measure the mean absolute per-cell difference from the corresponding full-data TF effect. Compare this error for targeted versus random thinning within each repeat. Keep effects for all three endpoints, not only whichever produces the largest difference.

MSC affecting the progenitor endpoint and MYOG affecting maturation endpoints are useful tutorial checks. Greater MYOG sensitivity is a hypothesis, not a required result: other branches may preserve its signal, and MSC can change indirectly.

Retain tutorial depletion likelihood as a secondary summary. It is a normalized rank/AUC comparison: values above 0.5 indicate a depletion-like shift in modeled fate probabilities. It is not a fraction of cells lost, knockout success probability, or terminal abundance forecast. A score of 0.67 does not mean 67% depletion. Loss of a significance flag alone is not evidence that the biological mechanism weakened.

### Supporting readout: unperturbed fate stability

On the same shared nonterminal cells, summarize the absolute difference between each reduced fit's unperturbed fate probabilities and the full fit's probabilities. This distinguishes altered baseline fate mapping from altered TF-response predictions.

### Supporting readout: local transition agreement

Compute CBC separately for the four prespecified transitions. `kernel.cbc(source, target, cluster_key, rep=...)` is a valid CellRank API and returns values for source boundary cells, not a single overall score.

Save a fixed reference PCA representation for evaluation only and attach its rows by cell ID as `X_eval_pca` when scoring. Continue rebuilding each condition's training PCA and graph. The fixed scoring representation is a shared yardstick; it is not a training input and does not turn this into independent validation.

CBC still depends on the condition's transition graph and available boundary cells. Record boundary-cell IDs, counts, finite-value counts, and the fraction of source cells with boundary support. Compare the intersection of evaluable source cells when possible, alongside the loss of boundary coverage. If no valid boundary remains, report unassessable/disconnected rather than zero or silently averaging an empty array. Interpret CBC as a graph-dependent pipeline diagnostic, not an isolated measure of regulatory-network accuracy.

Use a fixed display embedding for plots as well, or clearly label independently refitted embeddings. Visual rearrangement alone is not changed developmental direction.

### Secondary readout: stage association

Compute cell-level latent time using the documented aggregation consistently; the paper averages gene-specific latent times. Compare its Spearman association with PCW stage on the same shared cells, optionally within prespecified broad cell groups.

Specimen age is useful side information, not exact differentiation time for every cell. Do not calibrate latent time in hours, optimize the model using this metric and then call it an independent test, or treat an arbitrary 0.05 correlation change as a failure threshold.

## 9. Analysis and claims

- Plot every computational repeat and paired targeted-minus-random contrasts. Report effect sizes and the median/range across repeats.
- Three repeats do not justify strong population-level significance claims or treating thousands of correlated cells as independent biological replicates. Formal p-values are optional; within-fit tutorial p-values are descriptive, not the main evidence for the sampling comparison.
- The main comparison is whether targeted thinning consistently changes the selected readouts more than matched random nonterminal thinning. Quantify the size, not just the direction, of that difference.
- The fixed prior, selected genes, endpoint protection, sampling fraction, and two TFs limit generalization. Do not claim genome-wide driver-ranking accuracy or general superiority over other methods.
- If little changes, conclude: "We observed limited sensitivity under this particular thinning scheme and configuration across three computational repeats." Do not infer universal robustness or absence of smaller effects.

The original paper already includes extensive sampling, temporal-coverage, prior-ablation, and repeat-fit analyses. Describe this as a focused course extension examining structured intermediate-state loss and matched-cell TF effects, not as the first robustness analysis of RegVelo.

## 10. Two-week execution and stopping rules

| Time | Work and decision |
|---|---|
| Days 1-2 | Load and inspect data; pin a working environment; run one selected-configuration pilot including both TF analyses; measure cost. Confirm valid unsmoothed inputs. |
| Day 3 | Validate the clean rebuild on full data; choose target/fraction using counts and graph support; freeze anchors, settings, sampling protocol, and evaluation definitions. |
| Days 4-8 | Complete the nine main fits and eighteen separate perturbation analyses; save checkpoints, manifests, and metrics after each completed run. |
| Days 9-10 | Resolve failed runs and finalize comparison plots. Only add work if all required outputs already exist. |
| Days 11-14 | Interpret results, write methods/results/limitations, and prepare the presentation. |

Do not automatically rerun the tutorial's nine-fit model-selection grid. Reuse a pilot as a main full-data fit only if its preprocessing and all frozen settings exactly match the experimental protocol; otherwise label it separate pilot work.

If setup or valid preprocessing is unresolved by the end of day 2, or one full fit plus downstream analyses makes the nine-fit budget unrealistic, simplify before spending a week. Prefer a clearly labeled smaller pilot or a reproduction-focused comparison over a rushed claim of robustness. Missing unsmoothed inputs is a blocker for the thinning experiment, not a caveat that can be ignored. Colab or another available GPU environment is a possible execution option, not a guarantee that dependency problems disappear.

If the main experiment finishes early, prioritize more repeat blocks before adding a 90% thinning condition. More repetitions address the core variability question; an additional severity level expands the scope.

## 11. Deliverables and course connection

Minimum outputs:

1. A reproducible notebook/script and saved environment/version information.
2. A sampling table by cell type, stage, and sample, including exact retained counts.
3. Matched-cell MSC/MYOG effect plots for the three endpoint states.
4. Fate-stability and per-transition CBC plots with boundary-support counts.
5. A secondary latent-time/stage plot and a short report with explicit limitations.

The course connection is the inference of a coupled ODE system, constraints from regulatory priors, identifiability under incomplete observations, and Markov-chain fate inference. The analogy with conditioning in linear inverse problems is useful intuition, not a demonstrated diagnosis of this nonlinear model. Do not assign lecture numbers unless verified against the actual lecture materials.

## 12. Sources

- Supplied final paper: Wang et al. (2026), Cell 189, 3773-3800, especially Figure 3, STAR Methods on perturbation scores, and the supplementary robustness analyses. https://doi.org/10.1016/j.cell.2026.04.022
- Human limb tutorial and documented dataset/results: https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html
- Official workflow, including explicit neighbors before moments: https://regvelo.readthedocs.io/en/stable/getting_started.html
- Preprocessing API: https://regvelo.readthedocs.io/en/stable/api/genapi/regvelo.preprocessing.preprocess_data.html
- Preprocessing implementation: https://github.com/theislab/regvelo/blob/main/src/regvelo/preprocessing/_preprocess_data.py
- Current perturbation-summary API: https://regvelo.readthedocs.io/en/stable/api/genapi/regvelo.metrics.cellfate_perturbation.html
- CellRank velocity/CBC API: https://cellrank.readthedocs.io/en/stable/api/_autosummary/kernels/cellrank.kernels.VelocityKernel.html
- CBC implementation: https://cellrank.readthedocs.io/en/stable/_modules/cellrank/kernels/_base_kernel.html
- CellRank endpoint API: https://cellrank.readthedocs.io/en/stable/api/_autosummary/estimators/cellrank.estimators.GPCCA.html
- scVelo moments documentation: https://scvelo.readthedocs.io/en/stable/scvelo.pp.moments.html
- Example of the neighbor-reuse conditional in scVelo versions: https://github.com/theislab/scvelo/issues/1212

Documentation was checked on November 15, 2026. Pin actual installed versions and retain any implementation deviations with the results.
