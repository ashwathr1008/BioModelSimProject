# P5 analysis: generated results

**Estimand:** within-repeat common nonterminal cells; Full is a reference, not ground truth.

## Primary finding: unperturbed fate probability stability

| Repeat | Targeted mean TVD | Random mean TVD | T − R |
|---|---:|---:|---:|
| 0 | 0.06001 | 0.04217 | +0.01785 |
| 1 | 0.04774 | 0.02537 | +0.02237 |
| 2 | 0.02899 | 0.02140 | +0.00758 |

Across all three computational blocks, Targeted induced a larger mean fate-distribution shift on matched nonterminal cells. This does **not** establish a biological effect or isolate cell identity from developmental stage.

## Additional sensitivity checks

- R0, excluding MyoB1 as source: Targeted mean TVD=0.0568, Random=0.0404, n=8,055 identical cells.
- R1, excluding MyoB1 as source: Targeted mean TVD=0.0466, Random=0.0248, n=8,056 identical cells.
- R2, excluding MyoB1 as source: Targeted mean TVD=0.0279, Random=0.0205, n=8,066 identical cells.

The above reweighting/exclusion addresses whether the result is *only* present in remaining MyoB1 source cells, not the confounding of removed cell identities with developmental stage.

### Full-data repeat variation (unchanged cell input)

These comparisons change training seed (and may include finite inference-sampling effects). Their dispersion is descriptive, not a formal sampling null distribution.

| Full comparison | Mean TVD, evaluation cohort R0 |
|---|---:|
| full_r0_vs_full_r1 | 0.01871 |
| full_r0_vs_full_r2 | 0.04628 |
| full_r1_vs_full_r2 | 0.06012 |

The Full-vs-Full variation overlaps the Targeted-vs-Full variability in magnitude; avoid interpreting the contrast as a direct, isolated effect of deleting MyoB1.

## Existing P3 direction and coverage checks

P3 CBC means must be evaluated alongside source boundary counts; the evaluable MyoB1→MyoC1 boundary declines sharply under Targeted thinning. See the generated paired CBC tables.

## Existing P4 perturbation effects

- MSC → PAX7+MyoProg1: targeted-minus-random discrepancy (percentage points): R0=-0.60, R1=+0.74, R2=+0.90.
- MYOG → MYH3+MyoC: targeted-minus-random discrepancy (percentage points): R0=-4.61, R1=+2.13, R2=-0.44.
No uniformly larger Targeted discrepancy in P4 across the three computational blocks.

## Optional whole-network transition analysis

Not executed. P2 fitted.h5ad files and CellRank/RegVelo were not available; the ordinary GitHub ZIP excludes these large assets.

## Direct per-cell cross-layer comparison

Computed from all 18 original P4 compressed tables; see tables/cross_layer_per_cell_association.csv.

## Critical constraints

- Three computational repeats are not biological replicates. No cell-level p-values are interpreted as biological inference.
- The fixed 606-gene panel and prior GRN include information extracted from the original full study.
- Deleting MyoB1 also disproportionately removes earlier PCW stages; the current random control is not stage matched.
- P3 CellRank fate probabilities are conditional on fixed terminal anchors; P4 blockades manipulate learned TF regulatory weights, not empirical knockouts.
- Paired cells participated in model fitting; these results measure within-dataset prediction sensitivity, not held-out predictive accuracy.
