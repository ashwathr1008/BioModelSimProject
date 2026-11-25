# Do missing MyoB1 cells change RegVelo’s perturbation predictions?

Removing 75% of MyoB1 cells did not consistently change the predicted effects of MSC or MYOG more than removing the same number of random nonterminal cells. The predicted directions for MSC toward the PAX7+ progenitor fate and MYOG toward the MYH3+ mature fate held across all nine main fits. The magnitude varied across fits, especially for MYOG.

## Analysis

We used three fits each for the full dataset, targeted removal of 569 MyoB1 cells, and matched random removal of 569 nonterminal cells. For each fitted model, we separately silenced the outgoing regulatory weights of MSC and MYOG, drew 30 posterior samples, and recomputed fate probabilities. We used the same 106 terminal anchors and the same 80% velocity and 20% connectivity kernel mixture as the unperturbed analysis. For each cell and endpoint, the effect was the perturbed probability minus that model’s unperturbed probability. Within each repeat, we compared the same nonterminal cells across all three conditions. This gave 8,234, 8,235, and 8,243 cells in repeats 0, 1, and 2.

## Changes in predicted fate

On the shared cells, silencing MSC reduced the mean probability of reaching PAX7+MyoProg1 by 2.67 to 6.44 percentage points across the nine fits. Silencing MYOG reduced the mean probability of reaching MYH3+MyoC by 3.00 to 12.95 points. The corresponding full-data effects had medians of −5.60 and −8.05 points across the three repeats. The mean MYL3+MyoC change stayed within 0.07 points of zero for both regulators in every main fit. A near-zero signed mean can include opposing cell-level changes. [Signed effects by condition and repeat](../results/p4/signed_fate_changes.png)

## Targeted versus random removal

For each reduced fit we calculated the mean absolute difference from the full-data effect on matched cells. The table gives the targeted discrepancy minus the random discrepancy, in percentage points. A positive entry means the targeted removal changed the prediction more.

| Regulator | Endpoint | Repeat 0 | Repeat 1 | Repeat 2 | Median |
|---|---|---:|---:|---:|---:|
| MSC | PAX7+MyoProg1 | −0.60 | +0.74 | +0.90 | +0.74 |
| MSC | MYH3+MyoC | −0.57 | +0.73 | +0.88 | +0.73 |
| MSC | MYL3+MyoC | +0.056 | +0.012 | −0.020 | +0.012 |
| MYOG | PAX7+MyoProg1 | −4.61 | +2.12 | −0.41 | −0.41 |
| MYOG | MYH3+MyoC | −4.61 | +2.13 | −0.44 | −0.44 |
| MYOG | MYL3+MyoC | +0.104 | +0.010 | −0.006 | +0.010 |

All six contrasts changed sign across the three computational repeats. The largest swing concerned MYOG toward MYH3+MyoC. Random removal changed the effect more in repeat 0, while targeted removal changed it more in repeat 1. Within the surviving MyoB1 cells, the same contrast was −10.90, +6.04, and −6.27 points. These results support a consistent direction for the two main regulator–fate predictions, with uncertainty about how strongly thinning shifts their magnitudes. The full fit serves as a reference for stability. The three computational repeats describe variation for this dataset, gene panel, and thinning design. [Repeat-level comparison](../results/p4/perturbation_stability.png) · [Source-group results](../results/p4/shared_cell_effects_by_repeat.csv)

## Comparison with the human limb tutorial

Using the separately archived tutorial fit and Person 3’s fixed anchors, the depletion AUC was 0.568 for MSC toward PAX7+MyoProg1 and 0.662 for MYOG toward MYH3+MyoC. The published values were 0.672 and 0.785. MYOG toward MYL3+MyoC had an AUC of 0.525, close to the neutral score of 0.5 and well below the published 0.747. For MSC toward MYL3+MyoC the score crossed 0.5 in the opposite direction, from 0.439 in the tutorial to 0.564 here. The course pipeline uses its own fitted checkpoint, fixed anchors and mixed kernel, so the source of each difference cannot be isolated from this comparison. [Tutorial comparison](../results/p4/tutorial/depletion_check.csv) · [Published tutorial](https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html)

## Files for Charlie

The main run manifest identifies 18 compressed tables with 638,694 cell–regulator–endpoint rows in total. Join conditions using cell ID, repeat, regulator, and endpoint. Each table includes the unperturbed probability, perturbed probability, signed change, and source cell type. Probabilities in the CSV files are on the 0–1 scale. Multiplying a difference by 100 gives percentage points. [Run manifest](../results/p4/run_manifest.csv) · [Methods and reproduction commands](p4_perturbations.md)
