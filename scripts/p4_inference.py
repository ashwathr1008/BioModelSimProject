"""Accelerate the posterior-only output path and check it against RegVelo."""
from importlib.metadata import version
from math import ceil
from time import monotonic
from unittest.mock import patch

import numpy as np
import regvelo
import scvi
from regvelo._module import VELOVAE


def posterior_decoder(module, z, gamma, beta, alpha_1, latent_dim=None):
    # Both posterior output routines read only px_rho from generative().
    # get_velocity() then calls the original ODE solver with those times.
    return {"px_rho": module.decoder(z, latent_dim=latent_dim)}


def run_perturbation(folder, adata, regulator, draws, batch_size, seed):
    if version("regvelo") != "1.0.0":
        raise RuntimeError("The posterior decoder shortcut requires RegVelo 1.0.0.")
    options = dict(model=str(folder / "model"), TF=regulator, effects=0, cutoff=0, batch_size=batch_size)
    # Test the entire public perturbation call on real cells from this fit.
    idx = np.linspace(0, adata.n_obs - 1, min(batch_size, adata.n_obs), dtype=int)
    sample = adata[idx].copy()
    scvi.settings.seed = seed
    reference, _ = regvelo.tl.in_silico_block_simulation(adata=sample.copy(), n_samples=2, **options)
    with patch.object(VELOVAE, "generative", posterior_decoder):
        scvi.settings.seed = seed
        t0 = monotonic()
        check, _ = regvelo.tl.in_silico_block_simulation(adata=sample.copy(), n_samples=2, **options)
        elapsed = monotonic() - t0
        for layer in ["velocity", "latent_time_regvelo", "fit_t"]:
            np.testing.assert_array_equal(reference.layers[layer], check.layers[layer], err_msg=layer)
        print(f"Reference reproduced for {len(idx)} cells and 2 draws ({elapsed:.2f}s on the fast path).", flush=True)
        scvi.settings.seed = seed
        solve = VELOVAE._get_induction_unspliced_spliced
        total = ceil(adata.n_obs / batch_size) * draws
        completed = 0
        velocity_start = None

        def timed_solve(module, times):
            nonlocal completed, velocity_start
            if velocity_start is None:
                velocity_start = monotonic()
            values = solve(module, times)
            completed += 1
            if completed == draws or completed % (5 * draws) == 0 or completed == total:
                seconds = monotonic() - velocity_start
                remaining = seconds / completed * (total - completed)
                print(f"Velocity draws {completed}/{total}, estimated {remaining / 60:.1f} min left", flush=True)
            return values

        with patch.object(VELOVAE, "_get_induction_unspliced_spliced", timed_solve):
            perturbed, model = regvelo.tl.in_silico_block_simulation(adata=adata, n_samples=draws, **options)
        if completed != total:
            raise RuntimeError(f"Expected {total} ODE evaluations, observed {completed}")
    return perturbed, model, {"reference_cells": len(idx), "reference_draws": 2, "reference_exact": True}
