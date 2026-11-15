# Robustness of RegVelo to incomplete sampling during human muscle development

**Course:** 02-512 / 02-712 Computational Methods for Biological Modeling and Simulation  
**Duration:** Two weeks  
**Approach:** Reproduce an established analysis and perform a controlled sensitivity experiment

## 1. Overview and motivation

This project tests whether under-sampling an intermediate cell population changes RegVelo's predictions of muscle development more than removing the same number of cells randomly. Using the published human limb myogenesis dataset, we will examine developmental transitions, cell-fate probabilities, and the predicted effects of two transcription factors, MSC and MYOG.

Single-cell experiments provide snapshots of a developmental process. Transient intermediate populations can be rare or unevenly sampled, leaving gaps between progenitor and differentiated states. These gaps may affect both the inferred trajectory and the regulatory mechanisms inferred to drive it.

RegVelo connects RNA-processing dynamics with a gene regulatory network. For each gene, it models unspliced RNA production and splicing, together with spliced RNA production and degradation. Transcription depends on the expression of regulatory genes, coupling the genes into a system of ordinary differential equations. The model is fitted to the selected dataset. CellRank then uses its inferred dynamics to construct a Markov model of cell transitions and estimate probabilities of reaching specified endpoint states.

The project workflow is to:

1. Reproduce the selected model from the human muscle-development tutorial.
2. Establish a full-data baseline using a consistent preprocessing pipeline.
3. Partially remove one intermediate myoblast population.
4. Compare with removing the same number of randomly selected nonterminal cells.
5. Refit RegVelo and compare predictions on cells retained across conditions.
6. Examine separate simulated perturbations of MSC and MYOG.

The main experiment comprises **three conditions, three computational repeats per condition, nine model fits, and eighteen TF perturbation analyses**.

## 2. Research question and hypotheses

**With the published gene panel and regulatory prior held fixed, does selectively reducing one intermediate myoblast population change RegVelo's trajectory and MSC/MYOG perturbation predictions more than removing an equal number of randomly chosen nonterminal cells?**

The working hypothesis is that targeted thinning will disrupt some predictions more than matched random thinning because it removes information from a specific developmental transition. We will also examine whether trajectory predictions remain relatively stable while TF perturbation predictions change.

The complete-data model provides a reference for stability, rather than biological ground truth. Annotated transitions and developmental stages provide additional biological checks. This is a focused course extension of the paper's existing robustness analyses.

## 3. Dataset and biological system

Use the processed human lower-limb myogenesis dataset supplied through `rgv.datasets.human_limb()` and its accompanying tutorial.

| Item | Specification |
|---|---|
| System | Human embryonic/fetal limb muscle development |
| Developmental coverage | Post-conception weeks 5-9 |
| Documented size | 12,207 cells and 606 retained genes |
| Cell types | `adata.obs['leiden_R_celltype']` |
| Developmental stage | `adata.obs['adj_stage']` |
| Sample labels | `adata.obs['adj_sample']` |
| RNA measurements | Verify usable, unsmoothed per-cell `spliced` and `unspliced` layers |
| Regulatory inputs | Supplied `adata.uns['skeleton']` and TF annotations |
| Model | `soft_regularized`, with `lam2=1.0` |

The candidate target is **MyoB1**, on the annotated route `PAX3+MyoProg -> MyoB1 -> MyoC1`. The route `PAX3+PAX7+MyoProg -> MyoB2 -> MyoC2` provides a within-dataset comparison. Both routes share the fitted regulatory network, so MyoB2 is not an independent or necessarily unaffected control.

The endpoint/reference states for fate inference are:

- `PAX7+MyoProg1`: a progenitor endpoint or reservoir in this analysis.
- `MYH3+MyoC`: a differentiated myocyte endpoint.
- `MYL3+MyoC`: a differentiated myocyte endpoint.

An endpoint in the modeled dataset is not necessarily an irreversible biological terminal state. Developmental stage records specimen age, which is informative but does not specify the exact differentiation age of every cell. Sample labels will be interpreted using their metadata rather than automatically treated as independent donor identities.

## 4. Experimental design

### Conditions

| Condition | Definition |
|---|---|
| Full data | Retain all cells and run the complete preprocessing and fitting workflow. |
| Targeted thinning | Remove 75% of MyoB1 cells, subject to the feasibility check below. |
| Matched random nonterminal thinning | Remove exactly the same number of randomly selected nonterminal cells, including eligible MyoB1 cells. |

Retain all cells belonging to the three endpoint populations in both reduced conditions. This holds endpoint availability constant while testing sampling of the nonterminal part of development.

Before fixing the target and fraction, tabulate cell counts by type, stage, and sample and inspect support for the annotated transitions. Use MyoB1 at 75% removal if enough retained cells and boundary connections remain for evaluation. If the count/connectivity pilot shows that the condition is unassessable, reduce the fraction or use MyoB2 according to a documented count-based rule. Freeze this choice before examining comparative TF effects.

The primary protocol is uniform random selection within the target population or eligible nonterminal control pool. Inspect and report the resulting stage/sample composition. If specimen imbalance requires stratification, define and freeze the stratified sampling protocol during the pilot, before the main experiment.

### Repeats and fixed settings

Use three computational repeat blocks. Each block contains a full-data fit, a targeted-thinning fit, and a matched random-thinning fit. Record separate sampling and training seeds, and use the same training-seed label across the three conditions within a block.

Hold constant:

- The ordered gene panel, TF annotations, and aligned prior GRN.
- Model variant, regularization, optimizer, stopping rule, and training budget.
- Normalization policy, PCA settings and gene mask, neighbor count, and preprocessing seed.
- Velocity-output sampling settings.
- CellRank kernel construction, terminal anchor identities, and metric definitions.

Sampling and training variability are combined within these blocks. They are computational repeats, not independent biological replicates. Save exact retained/removed cell IDs and counts for every condition.

## 5. Implementation workflow

### Step 1: Environment and tutorial pilot

Establish a compatible environment containing RegVelo, scvi-tools, PyTorch, Scanpy, scVelo, and CellRank. Record package versions and hardware and preserve the working environment specification.

Run one fit using the tutorial's selected regularized configuration. Confirm successful training, plausible trajectories, valid fate probabilities, and separate MSC/MYOG perturbations. Measure runtime and memory for the entire workflow, including downstream analyses.

The tutorial's stage correlations around 0.63-0.68 and MSC/MYOG effect directions are reference checks. Numerical agreement is approximate because preprocessing, software, and stochastic fitting can differ. Investigate substantial discrepancies before proceeding. The training budget uses the selected configuration directly rather than a new model-selection grid.

### Step 2: Verify inputs and prepare subsets

Verify the transformations applied to `X`, `spliced`, and `unspliced`. Layer names do not establish that values are raw counts; normalized but unsmoothed per-cell layers may be suitable if their provenance is understood. Avoid repeated normalization or log transformation.

Check gene order, TF presence, GRN orientation, stage values, and cell-type/sample counts. Create sampling manifests containing the retained and removed cell IDs for each condition and repeat.

### Step 3: Rebuild preprocessing for every condition

Process full and reduced datasets through the same sequence:

1. Construct a clean AnnData object from retained cells' verified expression inputs, annotations, fixed genes, and prior network.
2. Apply the declared normalization policy.
3. Recompute PCA and explicitly rebuild neighbors on retained cells.
4. Compute the smoothed RNA layers `Ms` and `Mu` from the new graph.
5. Scale the model inputs without condition-specific gene selection. When using RegVelo's supported preprocessing helper, `filter_on_r2=False` preserves the fixed panel.
6. Check finite values, constant genes, graph connectivity, and gene/GRN alignment.

Do not carry inherited moments, graphs, velocity outputs, transition matrices, fate probabilities, or model-registration state into the new fits. Recomputing PCA alone does not guarantee that a moments call replaces an existing graph.

The rebuilt full dataset is the experimental baseline. It may differ from the shipped tutorial object, for example if the published graph was constructed before final gene filtering. Validate it using the same biological and numerical checks applied to reduced conditions.

### Step 4: Fit RegVelo and infer unperturbed fates

Fit and save each model using the frozen configuration. Generate velocities and latent times consistently.

Use the tutorial's CellRank mixture: **0.8 x VelocityKernel + 0.2 x ConnectivityKernel**.

Select terminal anchor cells once from the unperturbed rebuilt full-data reference, verifying agreement with the three endpoint labels. Save a dictionary mapping endpoint names to exact cell IDs. Reuse it in every condition through `GPCCA.set_terminal_states(anchor_dict)` before calculating fate probabilities.

Protecting the endpoint populations ensures these anchors survive thinning. Because the endpoints are specified, discovering endpoint identities is not an outcome of this experiment.

### Step 5: Simulate MSC and MYOG separately

For every fitted model, simulate MSC and MYOG independently, each starting from the saved unperturbed model. These simulations remove the selected TF's outgoing regulatory effects and generate model-based predictions.

Recompute each perturbed velocity kernel and fate probabilities with the same terminal anchors, mixture weights, and condition-specific connectivity kernel as its unperturbed analysis. Align all outputs by cell ID and endpoint name.

The tutorial associates MSC perturbation with lower probability of the progenitor endpoint and MYOG perturbation with lower probabilities of differentiated endpoints. Greater MYOG sensitivity to myoblast thinning is a hypothesis to examine, not a required outcome.

## 6. Evaluation

### Shared evaluation cells

Within each repeat, compare the same cell IDs retained in all three conditions. For fate and TF-effect summaries, use shared nonterminal cells and prespecified source-cell groups where sample sizes allow. This keeps the composition of evaluated cells constant across the compared fits.

These cells participated in model fitting. The evaluation measures stability rather than out-of-sample prediction accuracy.

### Primary outcome: TF perturbation-effect stability

For each cell, endpoint, TF, and sampling condition, calculate:

`delta_p = perturbed fate probability - unperturbed fate probability`.

Negative values indicate a lower modeled probability of reaching the endpoint after perturbation; positive values indicate a higher probability.

Report mean and median signed effects on the shared evaluation cells, including source-group summaries. For each reduced condition, calculate its mean absolute per-cell deviation from the corresponding full-data TF response:

`D_condition = mean(abs(delta_p_condition - delta_p_full))`.

Calculate the paired contrast `D_targeted - D_random` within each repeat. A positive contrast indicates greater departure from the full-data response after targeted thinning. It does not establish which model is biologically correct.

Report both TFs across all three endpoints. Retain the tutorial's depletion likelihood as a secondary summary: values above 0.5 indicate a depletion-like shift in modeled fate probabilities, not a fraction of cells lost or a direct forecast of lineage abundance.

### Supporting outcome: unperturbed fate stability

On the shared nonterminal cells, calculate absolute differences between each reduced fit's unperturbed fate probabilities and the full fit's probabilities. This separates changes in baseline fate mapping from changes in TF-response predictions.

### Supporting outcome: local transition agreement

Calculate cross-boundary correctness (CBC) separately for:

1. `PAX3+MyoProg -> MyoB1`
2. `MyoB1 -> MyoC1`
3. `PAX3+PAX7+MyoProg -> MyoB2`
4. `MyoB2 -> MyoC2`

Use a fixed reference PCA representation for evaluation, indexed by cell ID, while independently rebuilding each condition's training PCA and graph. The reference representation is used only as a common scoring coordinate system. Use a fixed display embedding for visual comparisons as well.

CBC provides scores for source cells with evaluable connections across a boundary. Record their identities, counts, fraction of source cells, and finite-value counts. Compare shared evaluable boundary cells when possible and report boundary coverage separately. Mark a transition with no valid boundary as unassessable or disconnected rather than assigning zero.

CBC reflects the condition's graph and boundary support as well as fitted dynamics, so it is a diagnostic of the complete analysis pipeline.

### Secondary outcome: developmental-stage association

Obtain cell-level latent time by consistently aggregating gene-specific latent times, following the documented mean-across-genes definition. Compute Spearman correlation with numeric post-conception stage on the shared cells.

Stage association provides biological context. Latent time is not calibrated in hours, and this metric will not be used to tune individual conditions or define an arbitrary failure threshold.

## 7. Analysis and interpretation

Plot every repeat and paired targeted-minus-random contrast. Summarize effect magnitudes with the median and range across repeats. Report cell counts and stage/sample composition for each condition.

The primary evidence is the magnitude, direction, and consistency of targeted-versus-random differences. A lost within-fit significance flag alone is insufficient to establish a weaker regulator effect because significance also depends on sample size and variability. Three computational repeats support an exploratory analysis, and individual smoothed cells are not independent biological replicates.

Possible results include larger disruptions after targeted thinning, stable trajectories with altered TF responses, similar sensitivity in both reduced conditions, or local transitions becoming unassessable. Each is informative for the research question.

Conclusions are conditional on the selected population, removal fraction, fixed genes/prior, model configuration, endpoints, and two TFs. Fixing the published genes and prior supplies background information that would not necessarily be available in a completely new experiment. Effects include preprocessing and graph changes as well as changes in fitted dynamics. The study therefore assesses pipeline sensitivity under these conditions rather than universal robustness or genome-wide driver accuracy.

## 8. Two-week schedule

| Days | Tasks | Deliverable |
|---|---|---|
| 1-2 | Establish environment; verify inputs; run one selected-configuration pilot including both TF analyses; measure cost | Working workflow and feasibility estimate |
| 3 | Validate rebuilt full data; freeze sampling, settings, anchors, and evaluation rules | Configuration and sampling manifests |
| 4-8 | Complete nine main fits and eighteen perturbation analyses; save outputs | Models, fate predictions, and metric tables |
| 9-10 | Resolve failed analyses; finalize comparisons and figures | Final figures and result summaries |
| 11-14 | Interpret results, write the report, prepare slides | Technical report, presentation, and reproducible analysis |

A pilot counts as a main full-data fit only if its preprocessing and settings match the frozen protocol. Otherwise, it is separate feasibility work. If time remains after the core experiment, prioritize additional repeat blocks before another thinning severity.

## 9. Risks and contingencies

| Risk | Response |
|---|---|
| Environment problems | Establish and pin a working environment during days 1-2; use available local, cluster, or notebook GPU resources as appropriate. |
| Unclear RNA-layer provenance | Trace preparation and transformations before rebuilding. |
| Only smoothed RNA inputs available | Obtain valid unsmoothed inputs; if unavailable within the schedule, use a reproduction-focused project instead of the thinning experiment. |
| Too few retained target cells or boundary connections | Apply the count/connectivity-based target or fraction rule before the main comparison. |
| Excessive runtime | Simplify early and describe any smaller study as exploratory; preserve time for analysis and writing. |
| Fate calculations fail | Check anchor IDs, graph reachability, probability validity, and solver behavior before interpreting results. |
| Small or inconsistent differences | Report effect sizes and computational variability without treating an inconclusive comparison as proof of general robustness. |

By the end of day 2, establish a successful end-to-end pilot or decide how to simplify. Runtime projections must include downstream fate inference and both TF perturbations.

## 10. Reproducibility and deliverables

Record versions, hardware, configuration, sampling/training seeds, retained/removed IDs, gene order, prior-network identity, terminal anchor IDs, and evaluation-cell definitions. Save checkpoints and outputs after each completed run. Validate row/column alignment, finite probabilities, and graph support before aggregation.

Deliverables:

1. A reproducible notebook or script workflow and environment specification.
2. A retained-cell table by condition, type, stage, and sample.
3. MSC/MYOG fate-probability effect plots and paired stability contrasts.
4. Unperturbed fate-stability and transition-specific CBC plots with boundary support.
5. A secondary latent-time/stage comparison.
6. A technical report and presentation explaining methods, findings, and limitations.

AI assistance can support implementation, debugging, documentation, and figure preparation. Data alignment, preprocessing, biological interpretation, and the distinction between simulated predictions and experimental evidence remain explicit verification steps.

## 11. Connection to course material

The project combines coupled differential equations, parameter inference, regulatory priors, incomplete observations, and Markov-chain fate models. It tests whether reduced information about part of a biological process destabilizes fitted dynamics or their mechanistic interpretation, connecting directly to identifiability and uncertainty in biological modeling and simulation.

## 12. References and resources

- Wang, W., Hu, Z., Weiler, P., et al. (2026). *RegVelo: Gene-regulatory-informed dynamics of single cells.* Cell 189, 3773-3800. [Paper](https://doi.org/10.1016/j.cell.2026.04.022).
- [Human limb tutorial](https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html).
- [RegVelo getting-started workflow](https://regvelo.readthedocs.io/en/stable/getting_started.html).
- [RegVelo preprocessing API](https://regvelo.readthedocs.io/en/stable/api/genapi/regvelo.preprocessing.preprocess_data.html).
- [RegVelo perturbation-summary API](https://regvelo.readthedocs.io/en/stable/api/genapi/regvelo.metrics.cellfate_perturbation.html).
- [CellRank velocity kernel and CBC API](https://cellrank.readthedocs.io/en/stable/api/_autosummary/kernels/cellrank.kernels.VelocityKernel.html).
- [CellRank GPCCA and terminal-state API](https://cellrank.readthedocs.io/en/stable/api/_autosummary/estimators/cellrank.estimators.GPCCA.html).
- [scVelo moments documentation](https://scvelo.readthedocs.io/en/stable/scvelo.pp.moments.html).
- [RegVelo reproducibility notebooks](https://theislab.github.io/regvelo_reproducibility).
- [RegVelo source repository](https://github.com/theislab/regvelo).
