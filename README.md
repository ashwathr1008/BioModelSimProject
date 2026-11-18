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
