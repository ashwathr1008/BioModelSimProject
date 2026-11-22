# BioModelSimProject: RegVelo under incomplete developmental sampling

- Final project for 02-512/02-712 Computational Methods for Biological Modeling and Simulation (CMU, Fall 2026).
- Question: if one intermediate cell population is under-sampled, do RegVelo's predictions change more than when the same number of random cells is removed?
- Model: RegVelo (Wang et al., *Cell* 2026), which couples RNA splicing dynamics with a gene regulatory network; CellRank turns its output into cell-fate probabilities.
- Data: human embryonic limb muscle development, 12,207 cells and 606 genes from post-conception weeks 5–9.
- Design: three conditions, each run three times (nine model fits):
  - full data;
  - 75% of MyoB1 myoblasts removed;
  - the same number of random cells removed.
- Endpoint cell types (PAX7+MyoProg1, MYH3+MyoC, MYL3+MyoC) are never removed.
- Readouts: predicted effects of switching off MSC and MYOG, cell-fate probabilities, transition scores, and latent time versus developmental stage.

See the [data preparation notes](docs/p1_data_notes.md) and [model fitting guide](docs/p2_model_fitting.md) for inputs, commands, and outputs.

Person 2 completed all nine main fits using recomputed PCA (`X`) and the separate
tutorial baseline on November 22. The selected settings are in
[`config/p2_x.yaml`](config/p2_x.yaml). Training records, QC, and comparison figures
are in [`results/p2/`](results/p2/).

Download the fitted datasets/checkpoints and reproducible inputs from the
[completed X-run release](https://github.com/ashwathr1008/BioModelSimProject/releases/tag/p2-x-completed-2026-11-22).
Extract the archives into the repository root. The large model and H5AD files are
release assets; they are not included in a normal clone. The fitting guide explains
verification and downstream use. Fate inference, transition scoring, and MSC/MYOG
perturbation analysis remain the next steps; the current comparison figures
summarize fitting and latent time.
