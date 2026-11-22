# RegVelo project chat handoff

> Historical handoff from before the completed GPU runs. The environment and run
> status below are superseded: all nine X fits and the tutorial completed on
> November 22. See [the current fitting guide and results handoff](docs/p2_model_fitting.md).
> Preserve this file as background, not as the current execution plan.

**November 21, 2026 · Ashwath Raghuram · CMU 02-712 Computational Methods for Biological Modeling and Simulation**

**Where the chat ended.** The human-limb dataset and all nine primary training inputs are ready locally. Ashwath’s fitting workflow has been written and tested with a stand-in model. The project has a Python 3.10 environment. Package downloads failed in the chat session, so the next task is to install RegVelo and verify it on the GPU. **Completed scientific fits: 0.**

## 1. The original request and Ashwath’s role

Ashwath asked Codex to sync the shared repository, review Aario’s contributions, and implement Person 2’s work for review in VS Code. His coding preferences were simple, efficient Python, descriptive names, straightforward control flow, and short comments explaining non-obvious decisions. The fitting script reuses Aario’s project helpers.

Repository: [ashwathr1008/BioModelSimProject](https://github.com/ashwathr1008/BioModelSimProject)

Local checkout: `C:\Users\ashwa\Documents\02-712 Computational Methods for Biological Modeling and Simulation\BioModelSimProject`. All relative paths below start here.

| Person | Owner | Responsibility |
|---|---|---|
| 1 | Aario (listed as Arrio in the assignment) | Inspect the data, define sampling, rebuild model inputs, and check preprocessing. |
| 2 | Ashwath | Set up the environment, reproduce the published tutorial baseline, fit the experimental models, record computational cost, and hand off each completed fit. |
| 3 | Po-Wei | Define shared terminal anchors and analyze cell fate and trajectories. |
| 4 | Kevin | Simulate MSC and MYOG perturbations separately for each fit using Po-Wei’s fate procedure. |
| 5 | Charlie | Align shared cells, compare conditions within repeats, and produce statistics and figures. |

The project asks whether removing cells from one intermediate myoblast population changes RegVelo’s predictions more than removing the same number of randomly chosen nonterminal cells. The core design is **three conditions × three repeats = nine fits**, followed by two perturbations per fit. Sampling seeds are 1000, 1001, and 1002. Training seeds are 0, 1, and 2, paired across conditions within each repeat. The rebuilt full-data model is the reference for stability.

## 2. Repository sync and Aario’s work

The local checkout started at `04c6b7d`. Aario’s data pipeline was committed in `6841586` and merged in `59ec64a`. Direct Git access from the shell failed DNS resolution, so Codex used verified Git objects to fast-forward the checkout. Local `main` and its `origin/main` tracking ref reached `23a5f9f`. A pre-existing local deletion of `regvelo_project_plan_audited.md` was preserved. The Person 2 changes from this chat remain uncommitted and unpushed.

Aario pushed the input-inspection and preprocessing scripts, nine manifests identifying retained and removed cells, QC tables, and comparison figures. His code rebuilds neighbors and the smoothed spliced and unspliced layers (`Ms` and `Mu`) independently on each retained dataset. Details are in `docs/p1_data_notes.md` and `results/p1/`.

| Condition | Cells retained | MyoB1 retained | Protected endpoint cells |
|---|---:|---|---:|
| Full | 12,207 | All 759 | 2,870 |
| Targeted | 11,638 | 190 after removing 569 MyoB1 cells | 2,870 |
| Random nonterminal | 11,638 | 713, 712, and 702 in repeats 0, 1, and 2 | 2,870 |

**The baseline decision.** Rebuilding PCA on the 606 project genes yields 11.5% average neighbor overlap with the shipped graph. Reusing the authors’ PCA coordinates and rebuilding neighbors yields 93.3%. The latter retains information from the full dataset in the PCA axes. Aario asked Ashwath to compare the two rebuilt full-data variants with the original tutorial at seed 0 before choosing the main experimental variant.

**The stage confound.** MyoB1 cells concentrate at earlier developmental stages, so targeted removal changes the stage distribution. Aario implemented an optional stage-matched random control. It remains disabled. Adding it would require three further fits and a team decision before the main comparison.

## 3. What Codex implemented and tested for Person 2

| Project file | Purpose |
|---|---|
| `scripts/p2_fit.py` | One entry point for the tutorial, full, targeted, and random fits. Checks the manifest, cell order, seeds, and prior alignment and records status for each run. |
| `config/experiment.yaml` | Holds training seeds, `lam2=1.0`, the 1,500-epoch cap, and 30 posterior samples. |
| `requirements-p2.txt` | Specifies the fitting stack, including RegVelo 1.0.0, scvi-tools 1.2.0, PyTorch 2.5.1, and NumPy 1.26.4. |
| `docs/p2_model_fitting.md` | Explains the pilot comparison, experimental run order, outputs, and checkpoint reloading. |
| `docs/install_regvelo_windows.md` | Gives VS Code Codex the Windows install commands and a real-data GPU initialization check. An identical copy was also present at the repo root when this handoff was prepared. |
| `README.md` and `.gitignore` | Link the guides and keep the local virtual environment and fitting locks out of Git. |

The script uses the tutorial’s soft-regularized model. It registers `Ms` and `Mu` and transposes the stored prior so RegVelo gets targets in rows and regulators in columns. The tutorial configuration uses full-batch training, AdamW, learning rate 0.01, a 90% training split, and early stopping on validation ELBO with patience 45. Cell-level latent time is the mean of the gene-scaled times in `fit_t`, rescaled across cells. The script records signed Spearman correlation with developmental stage. The tutorial’s published **absolute** correlations are 0.629, 0.679, and 0.672 for seeds 0, 1, and 2. These are reference values for the pilots. [Human-limb tutorial](https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html)

Integration checks used 72 real manifest cell IDs with synthetic expression, Aario’s preprocessing, real PyTorch tensors and AnnData files, RegVelo’s published output routine, and a stand-in for model training. They exercised all three conditions, the tutorial path, alternate PCA, prior orientation, preservation of cell and gene order, H5AD round trips, and training-history export. Tests also covered failed-run retries, completed-run skips, missing outputs, and simultaneous attempts to run the same fit. Scientific validation still requires an actual fit.

## 4. Why the data were missing and how that was resolved

Git ignores `data/` and `*.h5ad`. Syncing therefore brought down the scripts, manifests, and QC. The project loader can download and checksum the expression dataset. This session could not resolve Google Drive’s hostname, so Ashwath downloaded the file locally. [Source dataset](https://drive.google.com/file/d/17ZaqWx91w--iTtwWMR-SCh6wZO3bzMtM/view)

Codex confirmed that the download, `adata_processed_velo.h5ad`, matched Aario’s SHA-256, then renamed it to `data/raw/adata_humanlimb_preprocessed.h5ad`. The checksum is `27ee1272d31a13c57926f71c3456acfd9e81c8e2ef897dcd643c4879f7f2fb8c`. The file contains 12,207 cells, 606 genes, spliced and unspliced RNA, cell-type/stage/sample annotations, TF flags, and the regulatory network. It is the only source-data download required.

| Local directory | Verified contents |
|---|---|
| `data/raw/` | One original human-limb H5AD for the tutorial. |
| `data/prepared/` | Nine H5ADs. Each of `full`, `targeted`, and `random` has `r0`, `r1`, and `r2`. |
| `data/prepared_shipped_pca/` | One alternate `full_r0.h5ad` for the baseline pilot. |

Codex generated the ten prepared files with Aario’s builder and committed manifests. The local preprocessing environment used Scanpy 1.12.4 and NumPy 2.5.4. Aario used Scanpy 1.11.5 and NumPy 2.2.6. All ten files matched the expected cell IDs and order, protected endpoint counts, and all four transition-boundary counts. The default full-data neighbor-overlap statistic differed from his by about 0.000002. His committed QC tables were preserved.

## 5. Environment installation and the VS Code handoff

The machine has an **NVIDIA GeForce RTX 4090 Laptop GPU** with 16,376 MiB of memory. Codex found an existing 64-bit Python 3.10.11 installation and created a project-local virtual environment. The interpreter at `.venv\Scripts\python.exe` and pip 23.0.1 were verified. A separate local Python environment supplied the packages used for preprocessing.

The GPU PyTorch install failed when the shell could not resolve `download.pytorch.org`. A RegVelo package check independently failed to resolve `pypi.org`. The local caches lacked the required fitting packages. At this handoff, RegVelo and CUDA-enabled PyTorch still need installation in the project environment. GPU tensor execution and model initialization on the real baseline remain pending.

**Prompt to paste into VS Code Codex after opening the local checkout**

> Read `docs/install_regvelo_windows.md` and carry it out using the existing `.venv`. Install the dependencies, verify CUDA on my RTX 4090, and confirm RegVelo can initialize on `full_r0`. Stop before starting the full tutorial fit.

The Windows guide installs PyTorch from its CUDA 12.4 wheel index, installs the remaining packages from `requirements-p2.txt`, runs `pip check`, executes a CUDA tensor operation, and constructs the RegVelo model on the GPU using the real dataset.

## 6. Scientific work after installation

Run the three seed-0 pilots from the repository root in PowerShell.

```powershell
& .\.venv\Scripts\python.exe scripts/p2_fit.py tutorial
& .\.venv\Scripts\python.exe scripts/p2_fit.py full --repeat 0 --pca-input X
& .\.venv\Scripts\python.exe scripts/p2_fit.py full --repeat 0 --pca-input shipped_pca
```

1. Compare the shipped fit to the published tutorial. Evaluate the two rebuilt baselines with Po-Wei and Kevin using stage correlation, valid fate probabilities, and MSC/MYOG effect directions.
2. Record the PCA choice before examining removal effects. The selected rebuilt `full_r0` counts as the first of the nine experimental fits. Completing the three pilots and the nine-fit design therefore takes **eleven distinct fits** overall.
3. Finish repeat 0 with its targeted and random fits. If shipped PCA is selected, generate the remaining prepared inputs for that variant. Validate the first block through the downstream comparisons before starting repeats 1 and 2.
4. Hand off each completed fit as it becomes available. The checkpoint `model/model.pt` goes with its `fitted.h5ad`, which contains velocities and latent times. `training.csv` records training and validation ELBO. `run.json` records status, settings, seeds, versions, hardware, elapsed time, memory, warnings, and stage correlation.

The shipped fit will go under `results/p2/tutorial/`. Experimental fits will go under `results/p2/{pca_input}/{condition}_r{repeat}/`. Po-Wei uses the fitted velocities, Kevin reloads the matching checkpoint independently for each perturbation, and Charlie compares retained cells shared across the three conditions within each repeat. Runtime estimates for the complete workflow will also need their downstream timings.
