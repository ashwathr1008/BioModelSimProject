# Person 2: completed RegVelo fits and team handoff

All **nine main fits plus the tutorial baseline** completed successfully on
November 22, 2026, at 5:17 p.m. Eastern. The main experiment uses recomputed PCA
(`X`), three conditions (`full`, `targeted`, `random`), and three repeats (`0`, `1`,
`2`). The final eight-fit batch took 20 hours 20 minutes on an NVIDIA RTX 4090 Laptop
GPU. The tutorial and first full-data fit had completed earlier.

The final audit verified all 40 output files: matching input/settings identities,
finite training histories with the recorded epoch counts, readable fitted AnnData
files with the expected dimensions and layers, and finite cell-level latent times
and stages. The fitting and comparison queue has stopped.

Fate inference, transition CBC, and MSC/MYOG perturbation analyses are **not yet
included**. Those are the next handoffs to Persons 3 and 4, followed by Person 5's
combined analysis.

## Selected experiment

[`config/p2_x.yaml`](../config/p2_x.yaml) records the executed Person 2 settings and
is the default for the queue and comparison scripts. It retains Aario's conditions,
sampling manifests, protected endpoints, gene panel, prior GRN, and preprocessing
settings. Training seeds are 0, 1, and 2; sampling seeds are 1000, 1001, and 1002.
The same training seed is used across conditions within each repeat.

PCA, neighbors, and moments are rebuilt using each condition's retained cells.
Recomputed PCA was chosen to match the removal experiment. This is a design choice,
not a claim that it is more accurate than the supplied representation. The authors'
PCA coordinates were learned using all cells, including cells later removed.

The earlier `config/p2_both_pca.yaml` preserves an expanded, superseded plan.
Prepared shipped-PCA inputs exist, but their rebuilt full-data model and the
shipped-PCA experimental grid were not fitted. The tutorial baseline uses the
original supplied object and moments; it is not a substitute for that missing
rebuilt shipped-PCA fit. Delta jobs were cancelled before running, and all completed
fits in this handoff were trained locally.

## Files for the team

| Path | Purpose |
|---|---|
| `config/p2_x.yaml` | Executed settings for the X-only experiment. |
| `requirements-p2.txt` | Compatible RegVelo, PyTorch, scvi-tools, and analysis dependencies. |
| `scripts/p2_fit.py` | Fit one condition/repeat; validate inputs and saved results before reusing a completed fit. |
| `scripts/p2_memory.py` | Reduce ODE autograd memory while retaining the full-batch training objective. |
| `scripts/p2_run_all.py` | Run the selected fits sequentially, skip matching completed fits, and stop on errors. |
| `scripts/p2_compare.py` | Summarize fit status, runtime, stage correlations, and paired latent-time agreement. |
| `results/p2/X/` | Nine experimental output folders. |
| `results/p2/tutorial/` | Original tutorial baseline outputs. |
| `results/p2/preprocessing/` | Exact sampling manifests and input QC. |
| `results/p2/comparison/` | Final comparison tables and figures. |

Each experimental folder is named
`results/p2/X/{condition}_r{repeat}/`, for example `results/p2/X/targeted_r1/`.

| File in each fit folder | Contents |
|---|---|
| `model/model.pt` | Saved RegVelo checkpoint. |
| `fitted.h5ad` | Matching AnnData with the condition's cells, graph, moments, fitted velocities and latent times. |
| `training.csv` | Training and validation ELBO by epoch. |
| `run.json` | Status, settings, seeds, input identities, software, hardware, runtime, memory, and stage correlation. |

**Keep each checkpoint paired with its own fitted AnnData.** The large `.pt` and
`.h5ad` files are excluded from Git history. Download the archives from the
[completed X fits release](https://github.com/ashwathr1008/BioModelSimProject/releases/tag/p2-x-completed-2026-11-22):

- `regvelo-p2-x-fits-2026-11-22.zip`: all four artifacts for each of the ten fits.
- `regvelo-p2-x-inputs-2026-11-22.zip`: the raw source, nine exact prepared X inputs,
  and nine sampling manifests, for reproducing or validating training.
- `SHA256SUMS`: archive checksums to verify after downloading.

An ordinary `git clone` provides the code and small results, not these large files.
Extract archives into the repository root so the `results/p2/` and `data/` paths
are restored. Persons 3 and 4 can start with the fit archive; the input archive is
needed to reuse completed fits through the training script's input-hash checks.

### Po-Wei: fate and trajectory analysis

Read the matching `fitted.h5ad` for every condition/repeat. It contains `velocity`,
gene-level `fit_t` and `latent_time_velovi`, cell-level `obs["latent_time"]`, numeric
`obs["adj_stage_num"]`, and the condition-specific graph. Preserve cell IDs and use
the team's fixed terminal anchors. Fate probabilities and CBC still need to be
computed; completion of model fitting does not establish their validity.

### Kevin: MSC and MYOG perturbations

Reload each model with its matching AnnData in the pinned fitting environment:

```python
from pathlib import Path
import anndata as ad
from regvelo import REGVELOVI

folder = Path("results/p2/X/full_r0")
adata = ad.read_h5ad(folder / "fitted.h5ad")
model = REGVELOVI.load(str(folder / "model"), adata=adata)
```

Start from the unperturbed checkpoint separately for MSC and MYOG. Coordinate the
fate calculation and terminal anchors with Po-Wei. Return condition/repeat/cell-ID
tables for Charlie; expected tutorial directions are checks to investigate, not
outcomes to enforce.

### Charlie: comparisons and figures

The existing Person 2 outputs are:

| File in `results/p2/comparison/` | Contents |
|---|---|
| `run_summary.csv` | Status, seeds, dimensions, epochs, runtime, memory, and per-fit stage correlation for all ten fits. |
| `shared_cell_stage.csv` | Stage correlations on the common cells within each comparison group. |
| `paired_latent_time.csv` | Rank agreement between fitted latent times on common cells. |
| `comparison_status.csv` | Availability and common-cell counts for each comparison. |
| `comparison_notes.json` | Scope and interpretation limits. |
| `baseline_stage_correlation.png` | Tutorial and rebuilt X full-data baseline correlations. |
| `shared_cell_stage_correlation.png` | Three-condition stage correlations for each repeat. |
| `fit_runtime.png` | Recorded fit runtimes. |

The main repeat comparisons contain 11,104, 11,105, and 11,113 common cells for
repeats 0, 1, and 2. Person 2's current stage/latent-time summaries include endpoints;
apply the agreed common retained **nonterminal-cell** population for the primary
downstream effect comparison. Stage correlation is a secondary diagnostic, not
predictive accuracy or evidence of a regulator effect. Loss values across differently
preprocessed inputs are not a PCA-selection score. These computational repeats are
not independent biological replicates.

## Environment and reproduction

The completed fits used Python 3.12.6, PyTorch 2.5.1+cu124, RegVelo 1.0.0,
scvi-tools 1.2.0, and NumPy 1.26.4. Use `requirements-p2.txt` and the
[Windows installation guide](install_regvelo_windows.md). Package versions are also
recorded in every `run.json` and in `results/p2/environment-lock.txt`.

The model uses the tutorial's soft regularization (`lam2=1.0`), AdamW with learning
rate 0.01, a 90% training split, full-batch training, early stopping on validation
ELBO with patience 45, a cap of 1,500 epochs, and 30 posterior samples.

The executed inputs are `data/prepared_p2/{condition}_r{repeat}.h5ad`; these files
and the raw dataset are outside Git history. Their manifests are in
`results/p2/preprocessing/manifests/`. The source dataset is checked against its
configured SHA-256. Aario's original prepared files and QC remain separate from
the inputs rebuilt in the Person 2 environment.

Run commands from the repository root with the pinned environment active. Inspect
the selected scope and refresh comparisons from downloaded outputs without fitting:

```shell
python scripts/p2_run_all.py --dry-run
python scripts/p2_fit.py --config config/p2_x.yaml status --pca-input X
python scripts/p2_compare.py --config config/p2_x.yaml
```

To reproduce training, first obtain or rebuild the exact prepared inputs. Rebuilding
with the committed manifests uses:

```shell
python scripts/p1_03_build_inputs.py --config config/p2_x.yaml --pca-input X
```

Then run a single fit or the whole selected queue:

```shell
python scripts/p2_fit.py --config config/p2_x.yaml targeted --repeat 0 --pca-input X
python scripts/p2_run_all.py --config config/p2_x.yaml
```

The queue validates input hashes and scientific settings before reusing a completed
fit. A regenerated H5AD can differ byte-for-byte even if its contents are similar;
do not replace the original prepared inputs during a queue. A checkpoint and its
matching fitted AnnData suffice for downstream work without retraining.

### GPU memory method

`scripts/p2_memory.py` checkpoints the ODE calculation in 2,048-cell chunks.
Encoder/decoder normalization, across-cell losses, global sorting, and optimizer
updates still use the full batch. Completed runs recorded about 6 GiB of peak
PyTorch allocations; this is not total device memory usage.

[`results/p2/memory_validation.json`](../results/p2/memory_validation.json) records
the numerical comparison of outputs, loss, and parameter gradients against the
unmodified implementation within the stated tolerances. Floating-point ordering
can differ; bitwise equivalence is not claimed. The helper's SHA-256 is part of each
fit's saved identity.

## References

- [Aario's preprocessing and pilot notes](p1_data_notes.md)
- [RegVelo human-limb tutorial](https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html)
- [RegVelo 1.0.0 model comparison implementation](https://github.com/theislab/regvelo/blob/v1.0.0/src/regvelo/ModelComparison.py)
