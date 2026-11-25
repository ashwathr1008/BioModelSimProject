# MSC and MYOG perturbations

For each of the nine models, simulate MSC and MYOG separately with RegVelo 1.0.0. Each call reloads the checkpoint from the X fitting run with its matching fitted AnnData. The human limb tutorial specifies a zero cutoff. Set all nonzero outgoing weights of the selected regulator to zero, use 30 posterior draws as in the fitting configuration, and compute a new velocity field. These are model interventions on regulatory weights.

Use Person 3's shared fate inference with the 106 predefined terminal anchors. For each run, recompute the velocity kernel with the perturbed velocity, retaining the run's own graph and the 0.8 velocity / 0.2 connectivity weights. Align by cell ID and endpoint name. Subtract the corresponding published Person 3 baseline probability. A difference of −0.10 means a ten-percentage-point decrease in the model's probability of reaching that endpoint. The two TF runs use the same posterior seed within a repeat. The stored baseline was generated separately, so finite posterior sampling also contributes to small differences.

The tutorial associates MSC with progenitor maintenance and MYOG with terminal differentiation. A prior expectation is a decrease in the PAX7+ progenitor fate probability after MSC perturbation, and a decrease in the MYH3+ or MYL3+ mature fate probability after MYOG perturbation. Evaluate all three directions for each regulator. [Human limb tutorial](https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html)

The local inference shortcut omits RegVelo's unused expression likelihood during posterior decoding. The velocity routine still calls the original ODE solver with the decoder's sampled times, followed by the library's clipping and scaling. Each job first compares the shortcut to the complete public tutorial function using the same seed and actual cells from that model. The code requires identical latent times, scaled times and velocities before starting all 30 draws.

## Run

Start from team repository commit `f46fb71` or a compatible checkout containing Person 3's baselines and anchors. Extract the handoff ZIP into that checkout. Install `requirements-p2.txt` in a Python 3.12 environment and extract the fits archive into the repository root. Use each checkpoint with its own fitted AnnData. [Team repository](https://github.com/ashwathr1008/BioModelSimProject) · [Person 2 fits release](https://github.com/ashwathr1008/BioModelSimProject/releases/tag/p2-x-completed-2026-11-22)

From the repository root, run the two perturbations for the first full-data fit.

```sh
python scripts/p4_perturb.py full 0 MSC
python scripts/p4_perturb.py full 0 MYOG
```

Repeat for full, targeted and random with repeat numbers 0, 1 and 2. A new process for each TF starts from the unperturbed checkpoint. To check the platform against the published unperturbed probabilities, run `python scripts/p4_check_baseline.py`. After all 18 jobs finish, run the summaries and figures.

```sh
python scripts/p4_summarize.py
python scripts/p4_plot.py
python scripts/p4_validate.py
```

## Handoff

For a separate tutorial check, restore `results/p2/tutorial/model/model.pt` and the corresponding `fitted.h5ad` from the same release. Run `python scripts/p4_perturb.py tutorial 0 MSC` and `python scripts/p4_perturb.py tutorial 0 MYOG`. These use Person 3's tutorial fate table and the shared fixed-anchor definition. The matched sampling contrasts use the nine rebuilt X fits exclusively.

Each run writes one gzipped CSV per regulator under `results/p4/<condition>_r<repeat>/`. Each row represents one cell, regulator and endpoint. Columns include cell ID, annotation, condition, repeat, regulator, endpoint, baseline probability, perturbed probability and signed change (delta_p). Run settings and duration are recorded alongside the table.

For each repeat, compare only nonterminal cells retained in all three fits. The primary discrepancy is the mean absolute difference between a reduced model's delta_p and the corresponding full model's delta_p, calculated for the same regulator, endpoint and cells. Subtract the random discrepancy from the targeted discrepancy: positive values indicate a larger change after targeted removal. Retain the three repeat-level results. The full fit is a reference for stability. Cell-type summaries are descriptive.

The secondary rank-based depletion AUC compares unperturbed and perturbed fate probabilities for all cells within each fit, assigning tied pairs half weight. Values above 0.5 indicate a shift toward lower predicted probabilities after perturbation. It is a rank score. No cell-level significance test is used to infer biological replication.

[Human limb tutorial](https://regvelo.readthedocs.io/en/stable/tutorials/humanlimb/Limb_ModelComp.html)
