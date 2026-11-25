# Person 4 handoff — MSC and MYOG perturbations

Hi everyone, I’ve finished the MSC and MYOG perturbation analysis. I ran each regulator separately on all nine of Ashwath’s fitted models (full, targeted, and random, with three repeats each), giving 18 main analyses. I also ran both perturbations on the original tutorial fit. Each run starts from its own unperturbed checkpoint and uses Po-Wei’s fate function, the same 106 terminal anchors, and the same 80% velocity / 20% connectivity kernel.

Across all nine fits, blocking MSC lowered the mean predicted PAX7+MyoProg1 fate probability by 2.67–6.44 percentage points on matched nonterminal cells. Blocking MYOG lowered the mean MYH3+MyoC probability by 3.00–12.95 points. Whether targeted MyoB1 removal changed those effects more than random removal varied across repeats. One difference worth mentioning in the presentation is that our tutorial check gave a much weaker MYOG effect toward MYL3+MyoC (depletion AUC 0.525 versus 0.747 in the published tutorial).

The main files are listed below.

- [scripts/p4_perturb.py](../scripts/p4_perturb.py) runs one condition, repeat, and regulator. [scripts/p4_inference.py](../scripts/p4_inference.py) checks the faster inference path against RegVelo’s tutorial function before each run.
- [results/p4/run_manifest.csv](../results/p4/run_manifest.csv) lists all 18 main cell-level tables and their checksums. Each table includes cell ID, cell type, condition, repeat, regulator, endpoint, and the probabilities before and after perturbation. A negative `delta_p` means the perturbation lowered the predicted probability of that fate.
- [results/p4/shared_cell_effects_by_repeat.csv](../results/p4/shared_cell_effects_by_repeat.csv) contains the matched comparisons for all shared nonterminal cells and for each source cell type. [scripts/p4_summarize.py](../scripts/p4_summarize.py) and [scripts/p4_plot.py](../scripts/p4_plot.py) regenerate the summaries and both figures.
- [docs/p4_results.md](p4_results.md) contains the numbers and figures for the report. [docs/p4_perturbations.md](p4_perturbations.md) has the methods and reproduction commands. The [tutorial comparison](../results/p4/tutorial/depletion_check.csv) is reported separately.

Charlie, start with the run manifest if you want to merge the cell-level results with Po-Wei’s outputs. Within each repeat, keep the nonterminal cell IDs present in all three conditions, then compare `delta_p` for the same cell, regulator, and endpoint. The matched sets contain 8,234, 8,235, and 8,243 cells for repeats 0, 1, and 2. To use the finished comparisons directly, filter the repeat-level table to `source_group == "all_shared_nonterminal"`. A positive `targeted_minus_random_mae` means targeted removal changed the predicted effect more than random removal. Multiply by 100 to report percentage points.

All result tables and figures are in Git. Run `git pull origin main` to get them. The large Person 2 model downloads are only needed to rerun inference.

Kevin
