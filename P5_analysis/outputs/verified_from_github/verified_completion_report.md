# P5: verified analytical update

Data provenance: `results/p3/fate/*.csv`, `results/p3/cell_level_cbc/all_cell_cbc.csv`, and `results/p4/shared_cell_effects_by_repeat.csv` fetched from the connected GitHub repository. Aggregated input tables in this folder were derived from complete text blobs (not synthetic examples).

## Main findings

**1. Unperturbed fate stability.** In all three computational blocks, targeted removal produced larger mean TVD than equal-sized random removal on the same nonterminal cells.

| Repeat | Targeted TVD | Random TVD | T−R |
|---|---:|---:|---:|
| R0 | 0.0600 | 0.0422 | +0.0178 |
| R1 | 0.0477 | 0.0254 | +0.0224 |
| R2 | 0.0290 | 0.0214 | +0.0076 |

**2. Local trajectory CBC.** When CBC is compared on a common evaluable boundary-cell subset, both MyoB1-linked transitions have larger mean absolute changes under Targeted than Random in each of three repeats. This can be obscured by nearly unchanged unpaired CBC mean values. Boundary eligibility and size change substantially.

| Transition | R0 T/R | R1 T/R | R2 T/R |
|---|---:|---:|---:|
| PAX3+MyoProg → MyoB1 | 0.141 / 0.060 | 0.108 / 0.034 | 0.096 / 0.041 |
| MyoB1 → MyoC1 | 0.099 / 0.045 | 0.084 / 0.033 | 0.094 / 0.035 |
| PAX3+PAX7+MyoProg → MyoB2 | 0.047 / 0.048 | 0.049 / 0.037 | 0.050 / 0.038 |
| MyoB2 → MyoC2 | 0.037 / 0.036 | 0.045 / 0.030 | 0.040 / 0.037 |

**3. Stage stratification.** Within each repeat, targeted thinning causes larger mean baseline TVD than random thinning among early cells (PCW ≤6.1). This is descriptive; the stage effect is confounded with MyoB1 identity in the removal design.

**4. TF-effect stability.** P4 contrasts for the two headline axes have mixed signs over repeats. Cross-layer mean absolute discrepancies are combined in `cross_layer_marginal_errors_verified.csv` on the same shared nonterminal cohort. These marginal means cannot establish a *per-cell association* between baseline error and TF-effect error. For the per-cell joint relationship, the original 18 compressed P4 result tables are needed.

**5. 10×10 transition matrices.** Not computed because none of the nine `results/p2/X/*/fitted.h5ad` files exist in this runtime and the RegVelo/CellRank environment is unavailable. The necessary model archive is 1,429,028,483 bytes in the P2 GitHub release. The tested Python implementation in `run_analysis.py --matrix` will compute these after installation and extraction. Do not substitute CBC matrices or invented transition probabilities.

## Constraints

- All three repeats come from one biological dataset; they are not independent biological replicates. No population-level significance claims.
- Full is a computational reference, not ground truth.
- GRN prior and gene panel were fixed using original data information.
- Stage and cell identity are confounded under the present targeted vs random removal design.
- P4 manipulates TF outgoing regulatory weights in silico; it does not reproduce experimental knockout effects.
