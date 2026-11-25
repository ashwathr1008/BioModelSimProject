# P5 interim results — verified from complete GitHub CSVs

Date: 2026-11-25. This is **not** a claim of having executed inside the user's Mac Desktop folder. Input files were read from the connected GitHub repository's exact blob content. Numeric tables have been rounded to 8 decimal places; the full local script regenerates them with native precision.

## A. Baseline, unperturbed fate probability sensitivity

A common set of **8,234 / 8,235 / 8,243 nonterminal cell IDs** was selected within R0 / R1 / R2 across Full, Targeted, Random. Fate probabilities towards the three P3 endpoints were compared with the Full fit on the **same cells**.

\[
\mathrm{TVD}_i^{c}=0.5\sum_{k\in\{\mathrm{PAX7,MYH3,MYL3}\}}|p_{i,k}^{c}-p_{i,k}^{Full}|.
\]

| Repeat | Targeted mean TVD | Random mean TVD | T − R | Targeted flip rate | Random flip rate |
|---|---:|---:|---:|---:|---:|
| R0 | 0.06001 | 0.04217 | **+0.01785** | 3.49% | 3.17% |
| R1 | 0.04774 | 0.02537 | **+0.02237** | 2.89% | 2.70% |
| R2 | 0.02899 | 0.02140 | **+0.00758** | 1.84% | 0.90% |

Targeted makes unperturbed fate distributions more different from Full than Random in all three computational blocks. This is a **descriptive within-dataset sensitivity result**, not independent biological evidence.

**MyoB1 subset** (cells still present in all conditions):

| Repeat | n MyoB1 | Targeted mean TVD | Random mean TVD | Difference |
|---|---:|---:|---:|---:|
| R0 | 179 | 0.20628 | 0.11962 | +0.08666 |
| R1 | 179 | 0.09718 | 0.05133 | +0.04585 |
| R2 | 177 | 0.07765 | 0.06409 | +0.01356 |

This local effect is stronger than the overall average, especially in R0. But not **every** cell type shows targeted excess in every repeat, so stratified results are needed. See `baseline_fate_by_celltype.csv`.

**Directional warning:** Although Total Variation excess is consistently positive, signed fate shifts are not consistent. For Targeted vs Full, the mean MYH3+MyoC probability shifts by **−5.93 pp in R0**, **−4.54 pp in R1**, but **+2.58 pp in R2**. Do not claim a consistently lower mature fate probability after MyoB1 thinning.

## B. Link to P4 perturbation sensitivity

P4 already computed *mean absolute discrepancy of TF effects* on the same shared nonterminal cells, and compared Targeted minus Random.

| Repeat | MSC → PAX7+MyoProg1 (pp) | MYOG → MYH3+MyoC (pp) |
|---|---:|---:|
| R0 | −0.60 | −4.61 |
| R1 | +0.74 | +2.13 |
| R2 | +0.90 | −0.44 |

Thus, a **consistent Targeted excess in unperturbed fate instability** coexists with **non-consistent Targeted excess in simulated TF-effect instability**. These are different quantities, and the size of errors is not directly attributable to distinct biological mechanisms. To establish a direct cross-layer comparison, run the full package against the raw P4 files, producing matched per-cell baseline-versus-TF-effect error tables.

## C. Limitations and next steps

1. The 3 repeats share the same original biological dataset. They do not give valid biological-level p-values.
2. MyoB1 is an early developmental population. Targeted removal also changes PCW sampling, so stage effects may confound cell identity.
3. Fixed GRN prior and gene panel retain upstream knowledge from complete data.
4. P3 fate probabilities depend on designated absorbing anchors; P4 regulatory blockade is a *model intervention*.
5. Full model is a reference, not an experimentally observed ground truth.
6. Execute `run_analysis.py` against the extracted project to regenerate full-resolution, cell-level and stage-stratified analyses; extract P2 fitted H5AD release assets for optional 10×10 matrix analysis.
