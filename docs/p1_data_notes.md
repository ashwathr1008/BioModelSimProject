# Person 1: data, sampling and preprocessing notes

Status as of 2026-11-18. All numbers below come from the scripts in `scripts/p1_0*.py`
run on the shipped object (sha256 `27ee1272…f7f2fb8c`).

## How to run

```bash
pip install -r requirements-p1.txt
python scripts/p1_01_inspect_inputs.py   # provenance, GRN and count checks -> results/p1/inputs/
python scripts/p1_02_make_manifests.py   # retained/removed cell IDs        -> results/p1/manifests/
python scripts/p1_03_build_inputs.py     # rebuilt model inputs + QC        -> data/prepared/, results/p1/qc/
python scripts/p1_04_compare_conditions.py  # before/after comparison plots -> results/p1/figures/
```
Step 1 downloads the source file to `data/raw/` if it is missing. The whole pipeline
takes about 2 minutes on a laptop CPU. `data/` is gitignored; manifests and QC are small and committed.

## 1. What the input layers are (verified, not assumed from names)

| Layer | What it is | Evidence |
|---|---|---|
| `spliced`, `unspliced` | Raw counts divided by each cell's total (`initial_size_spliced/unspliced`, computed on the full gene set) times a global constant (6601 / 1569) | Multiplying back recovers integers for 100% of nonzero entries |
| `matrix` | Raw total counts | Equals recovered spliced + unspliced + `ambiguous` exactly |
| `Ms`, `Mu` | kNN moments of `spliced`/`unspliced` over the shipped graph, then per-gene min-max scaled | Reproduced to 4e-7 by our recipe (binarized connectivities + self-loop, row-normalized) |
| `X` | Per-cell normalized log expression (min 0, max 6.56); not log1p(spliced) | Exact provenance unknown; no negative values, so no regression/scaling |

**Decision:** use `spliced`/`unspliced` as shipped. Their normalization is per cell (plus one
global constant that min-max scaling cancels), so removing other cells does not change them.
No re-normalization or extra log transform.

## 2. Gene panel and regulatory prior

- 606 genes, 77 annotated TFs. `uns['skeleton']`: 606×606 binary, index and columns equal `var_names`.
- Orientation: **rows = regulators, columns = targets** (only TF rows have edges). RegVelo's
  `ModelComparison` passes `W.T` (targets × regulators), which is consistent.
- 3,834 edges, no self-loops, every gene has at least one incoming edge.
- MSC: TF, 283 targets. MYOG: TF, 297 targets.
- Gene order and skeleton are identical in every prepared file (checked automatically).

## 3. Counts, stages and samples

| Cell type | n | | Cell type | n |
|---|---|---|---|---|
| PAX7+MyoProg2 | 3627 | | MyoC2 | 639 |
| PAX7+MyoProg1 (endpoint) | 2317 | | PAX3+MyoProg | 583 |
| PAX3+PAX7+MyoProg | 1790 | | MyoC1 | 412 |
| MyoB2 | 1527 | | MYL3+MyoC (endpoint) | 374 |
| **MyoB1 (target)** | **759** | | MYH3+MyoC (endpoint) | 179 |

- 10 stages (PCW 5.1–9.3), 25 samples. **Each sample belongs to exactly one stage.** Early stages
  are split by limb region (Proximal/Middle/Distal, Upper/Lower), so samples are tissue
  pieces or replicates, not independent donors.
- Full tables: `results/p1/inputs/celltype_by_stage.csv`, `celltype_by_sample.csv`, `sample_by_stage.csv`.

## 4. Sampling conditions (frozen defaults in `config/experiment.yaml`)

- Targeted: remove 569 of 759 MyoB1 cells (75%), so **190 retained** (feasibility minimum 100: pass).
- Random: remove 569 cells uniformly from all non-endpoint cells (MyoB1 included). This removes
  46–57 MyoB1 cells per repeat.
- Endpoints (`PAX7+MyoProg1`, `MYH3+MyoC`, `MYL3+MyoC`) are never removed (asserted).
- Sampling seed for repeat r = 1000 + r; training seed = r. Each condition draws from its own
  seeded stream (`SeedSequence([seed, condition_code])`).
- Removed sets differ across repeats: Jaccard about 0.6 for targeted (it is 75% of one
  population) and about 0.03 for random.
- `full_r0/r1/r2` contain identical cells; they differ only in the training seed.

### Stage confound: needs a team decision

MyoB1 is an early population: 449 of the 569 removed targeted cells come from PCW 5.4–6.1,
whereas random removal mostly takes late cells (PCW 8.4–9.3 dominate). See
`results/p1/figures/removed_by_stage_r0.png`. So "targeted vs random" also means "early vs
late cells removed". This matters most for the latent-time vs stage outcome.

Option: add `random_stagematched` (already implemented). It removes the same number of random
non-endpoint cells with the same stage distribution as the targeted removal. Enable it by adding it to
`sampling.optional_conditions`. That is 3 extra fits, so Person 2 must agree on the compute.
**Decide before the main runs.**

## 5. Rebuild recipe and the baseline caveat

For every condition (full data included), `src/p1_data/rebuild.py` builds a clean AnnData from
retained cells and runs these steps:
1. PCA (50 comps, arpack, seed 0) on `X`
2. Exact kNN graph (30 neighbours, 30 PCs, UMAP connectivities, sklearn transformer)
3. `Ms`/`Mu` from the new graph
4. Per-gene min-max scaling

Nothing graph-derived is inherited. The shipped PCA and embeddings are kept only as
`X_pca_ref`, `X_umap_ref` and `X_draw_graph_fa_ref`, for common-coordinate CBC scoring and
display (Persons 3 and 5).

**Finding:** the shipped graph was not built from a PCA of these 606 genes. Only PC1 matches
(|r| = 0.95); PCs 2–5 do not. So it was most likely built before final gene filtering.
Consequences:
- The rebuilt full-data graph shares only about 12% of kNN neighbours with the shipped graph.
- The rebuilt `Ms`/`Mu` correlate 0.83 / 0.78 with the shipped layers.
- **The rebuilt full baseline is therefore not the tutorial object.**

A sensitivity variant reuses the authors' PCA coordinates but still rebuilds the graph and
moments on retained cells only:
`python scripts/p1_03_build_inputs.py --pca-input shipped_pca` (writes to `data/prepared_shipped_pca/`).
- Its full-data graph overlaps the shipped one at 93%, so it is close to the tutorial.
- Caveat: the PCA axes were fit with all cells, so removed cells still shape the coordinates
  (but not the graph or the smoothing).
- `full_r0`, `targeted_r0` and `random_r0` are already built for this variant.

**Ask for Person 2's pilot:** fit `full_r0` from both `data/prepared/` and
`data/prepared_shipped_pca/` with the tutorial configuration. Compare the latent-time vs stage
correlation (tutorial: 0.63–0.68) and the MSC/MYOG directions. Freeze one variant **before**
looking at any targeted/random contrast.

## 6. QC on the prepared inputs (`results/p1/qc/qc_summary.csv`)

All 9 files: finite values, no constant genes in `Ms`/`Mu`, one connected graph component, gene
order and skeleton identical to source, endpoints intact.

Boundary support is the number of source cells with at least one neighbour in the target
cell type, i.e. the cells CBC can score:

| Transition | full | random (r0–r2) | targeted (r0–r2) |
|---|---|---|---|
| PAX3+MyoProg → MyoB1 | 222 | 201–218 | 148–165 |
| MyoB1 → MyoC1 | 163 | 151–159 | **52–63** |
| PAX3+PAX7+MyoProg → MyoB2 | 764 | 719–739 | 786–792 |
| MyoB2 → MyoC2 | 469 | 441–456 | 479–481 |

MyoB1 → MyoC1 stays evaluable (more than 50 boundary cells) under targeted thinning. That is the
transition where an effect is most expected.

## 7. Hand-off contract to Person 2: `data/prepared/{condition}_r{repeat}.h5ad`

- Layers: `spliced`, `unspliced` (as shipped); `Ms`, `Mu` (rebuilt, min-max scaled; train on these).
- `var`: same 606 genes in the same order, `var['TF']` (bool).
- `uns['skeleton']`: regulators × targets DataFrame (pass `W.T` to RegVelo, as `ModelComparison` does).
- `obs`: all original columns + `condition`, `repeat`, `sampling_seed`, `training_seed`.
- `obsp`/`obsm['X_pca']`: the rebuilt graph and PCA. Reference coordinates are in `*_ref` keys.
- `uns['p1_prep_json']`: all preprocessing and sampling settings, plus package versions.
- Exact cell IDs: `results/p1/manifests/{cid}.csv.gz` (`cell_id`, `status`, `celltype`, `stage`, `sample`, seeds).

Not included (by design): velocities, latent time, transition matrices, fate probabilities.
These are recomputed downstream for each condition.

## 8. Before/after comparison plots (`results/p1/figures/`)

- `umap_removed_cells_r0.png`: targeted removal is one compact region (MyoB1); random removal is spread across the map.
- `boundary_support_by_condition.png`: targeted removal cuts MyoB1 → MyoC1 support from 163 to about 57 cells and PAX3+MyoProg → MyoB1 from 222 to about 156. Random removal barely changes them. The MyoB2 route is unaffected by either.
- `delta_Ms_vs_full_by_celltype.png`: for cells kept in every condition, the smoothed input Ms of the remaining MyoB1 cells changes about 3× more under targeted removal (mean |ΔMs| ≈ 0.024) than under random removal (≈ 0.009). Other cell types change by similar small amounts in both. Note that Ms is min-max scaled per condition, so part of every change is rescaling.
- `phase_portraits_r0.png`: MYOG and MSC Ms-vs-Mu plots. The MyoB1 cells (light blue) on the MYOG induction branch thin out under targeted removal.
