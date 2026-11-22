"""Bound RegVelo's ODE autograd memory while retaining full-batch training."""

from importlib.metadata import version
from types import MethodType

import torch
from torch.utils.checkpoint import checkpoint


def enable_ode_checkpointing(model, chunk_size=256):
    """Checkpoint cell chunks inside the ODE only; leave the full-batch loss intact.

    Apply after model construction, before training. All genes of each cell stay
    together because their ODEs are coupled. Encoder/decoder normalization and
    across-cell loss terms still see the complete training batch. This changes
    memory use and floating-point summation details, not the intended objective.
    """
    if version("regvelo") != "1.0.0" or version("torchode") != "1.0.1":
        raise RuntimeError("ODE checkpointing was validated only for RegVelo 1.0.0 / torchode 1.0.1")
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or chunk_size < 2:
        raise ValueError("ODE chunk_size must be an integer of at least two cells")
    module = model.module
    if type(module).__module__ != "regvelo._module" or type(module).__name__ != "VELOVAE":
        raise TypeError("Expected the RegVelo VELOVAE module")
    if hasattr(module, "_p2_ode_chunk_size"):
        raise RuntimeError("ODE checkpointing is already enabled on this model")
    original = module._get_induction_unspliced_spliced

    def chunked_induction(self, t, eps=1e-6):
        if t.ndim != 2 or t.shape[0] < 2 or t.shape[1] != self.n_targets:
            raise ValueError("ODE checkpointing requires at least two cells and the complete ordered target-gene panel")
        if self.v_batch.num_g != self.n_targets or self.x0.shape != (self.n_targets, 2):
            raise ValueError("Unsupported RegVelo ODE state shape")

        def solve(times):
            unspliced, spliced, _ = original(times, eps=eps)
            return unspliced, spliced

        unspliced_parts, spliced_parts = [], []
        start = 0
        while start < t.shape[0]:
            stop = min(start + chunk_size, t.shape[0])
            # The upstream one-cell branch has different shape semantics.
            if t.shape[0] - stop == 1:
                stop += 1
            times = t[start:stop]
            if torch.is_grad_enabled():
                unspliced, spliced = checkpoint(solve, times, use_reentrant=False)
            else:
                unspliced, spliced = solve(times)
            unspliced_parts.append(unspliced)
            spliced_parts.append(spliced)
            start = stop

        # The correlation term consumes global cell indices, never chunk indices.
        index = torch.sort(t, dim=0).indices.T
        return torch.cat(unspliced_parts, dim=0), torch.cat(spliced_parts, dim=0), index

    module._get_induction_unspliced_spliced = MethodType(chunked_induction, module)
    module._p2_ode_chunk_size = chunk_size
    return dict(method="checkpoint_ode_by_cell", chunk_size=chunk_size,
                use_reentrant=False, training_batch="full", checkpoint_scope="ODE only")
