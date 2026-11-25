# P5 scientific synthesis — source-derived results

**Date:** 2026-11-25. **Data:** connected public team GitHub repository (`ashwathr1008/BioModelSimProject`) and its P3 fate/CBC results plus P4 summary tables. **Unit of analysis:** matched source cells within a computational repeat. The three blocks are not independent biological replicates.

## Main question

With fixed GRN prior and 606-gene panel, does deleting 75% of intermediate MyoB1 cells change RegVelo's inferred fate, trajectory, and in-silico TF-response predictions more than deleting the same number of randomly selected nonterminal cells?

## Core matched-cell results

| Repeat | Shared nonterminal cells | Targeted mean baseline fate TVD | Random mean baseline fate TVD | T−R |
|---|---:|---:|---:|---:|
| R0 | 8,234 | 0.06001 | 0.04217 | +0.01785 |
| R1 | 8,235 | 0.04774 | 0.02537 | +0.02237 |
| R2 | 8,243 | 0.02899 | 0.02140 | +0.00758 |

Fate TVD = half the L1 distance between the **three-endpoint probability vectors** for the *same* cell under a reduced fit and its Full fit. Positive T−R in all three computational repeats means targeted removal caused a larger **within-dataset prediction shift** in these runs, not that the targeted model is less biologically correct.

### Is the result driven solely by remaining MyoB1 cells?

No. Recomputing the cell-weighted TVD among common nonterminal cells **excluding all remaining MyoB1 sources** gives:

| Repeat | Common non-MyoB1 sources | Targeted mean TVD | Random mean TVD | T−R |
|---|---:|---:|---:|---:|
| R0 | 8,055 | 0.05676 | 0.04045 | +0.01632 |
| R1 | 8,056 | 0.04664 | 0.02479 | +0.02185 |
| R2 | 8,066 | 0.02792 | 0.02046 | +0.00745 |

Thus, the difference affects predictions for other source cell populations; it is **not** evidence of a specific intercellular signaling mechanism. Both graph/moment preprocessing and model fitting changed after thinning. The seven-cell-type **equal-weighted** sensitivity check also retains a positive T−R in all three repeats; see `fate_composition_robustness_verified.csv`.

### Essential qualification: Full model repeat variability

On the **same 8,234 nonterminal cells** used in the R0 three-condition contrast, evaluating different Full fits yields:

| Full-model comparison | Mean fate TVD |
|---|---:|
| Full r0 vs Full r1 | 0.01871 |
| Full r0 vs Full r2 | 0.04628 |
| Full r1 vs Full r2 | 0.06012 |

These three Full models are fit on the same cells with different training seeds; finite posterior inference sampling may contribute to variation. In other words, between-fit variation with no cell deletion is sometimes as large as or larger than the apparent targeted-versus-Full discrepancy. **This constrains claims of specificity and robustness.** The within-repeat Targeted/Random comparison shares a training-seed label with its Full fit, so cross-seed Full-vs-Full is a diagnostic reference, not a matched statistical null for that contrast.

For MyoB1 source cells, Full repeat variability is particularly large: on the R0 evaluation cohort, mean TVD ranges from 0.0345 to 0.1267 across Full seed pairs. See `full_training_seed_variability_verified.csv` and `figures/verified_full_seed_variability.png`.

## CBC and trajectory

Matched P3 CBC boundary-cell results show targeted removal produces larger absolute CBC changes than random removal for both **MyoB1-linked paths** in each computational block, notwithstanding small shifts in the unpaired mean CBC. For MyoB1→MyoC1, matched boundary cells number 31, 30 and 22 across repeats, and Targeted versus Random mean absolute CBC changes are 0.099/0.045, 0.084/0.033 and 0.094/0.035. These are graph-dependent model diagnostics; changing boundary eligibility matters.

The previously computed common-cell latent-time rank correlations are high, but global rank consistency cannot establish cell-fate, velocity, GRN-parameter or perturbation-response stability.

## TF response (P4, imported—not rederived from compressed per-cell files)

P4's headline perturbation-effect discrepancy contrasts (Targeted MAE − Random MAE, in probability points) switch sign across repeats:

| Repeat | MSC→PAX7+MyoProg1 | MYOG→MYH3+MyoC |
|---|---:|---:|
| R0 | −0.00602 | −0.04608 |
| R1 | +0.00740 | +0.02128 |
| R2 | +0.00899 | −0.00444 |

The main *mean* directions of MSC and MYOG blockades remain consistent across fits, but there is no consistent targeted excess in the magnitude of the **TF-effect error**. The baseline-fate and perturbation-effect MAEs have different definitions; P4 effect error subtracts two simulated perturbation effects, each of which itself subtracts a baseline prediction. Therefore, per-cell correlations between them require care: mathematical coupling can arise from shared baseline terms.

## What this supports

Structured loss of a developmental intermediate cell population causes more baseline fate divergence than equal-sized random thinning in all three computational blocks, and it changes the evaluation coverage and paired CBC of MyoB1-linked paths. However, the regulator-blockade excess error is not consistently targeted-dominant. Crucially, there is substantial Full-model repeat variability even with unchanged cell inputs.

## What this does **not** establish

- The learned GRN edges are causally identifiable or biologically correct.
- A true MSC/MYOG knockout will reproduce the model's simulated fate changes.
- Selectively missing MyoB1 cell **identity**, rather than over-represented missing **early developmental time**, uniquely caused the contrast.
- Thousands of correlated cells or three computational seeds amount to independent biological replication.
- The complete 10×10 CellRank transition structure was measured here; it needs the nine P2 fitted AnnData inputs.
- Direct *per-cell* association between P3 baseline fate error and P4 perturbation effect error has already been computed; the 18 P4 compressed original tables were not available in this sandbox.

## Recommended report framing

**Constrained pipeline sensitivity, not a general robustness claim.** Present the consistent Targeted baseline shift together with the substantial Full-fit variation. Show P3 CBC coverage and matched-boundary errors, then the mixed-sign P4 effect contrast. Conclude that fate prediction and simulated TF response are not interchangeable robustness readouts. Use the P1 developmental-stage composition imbalance, fixed prior, non-independent repeats and differing preprocessing as explicit caveats.
