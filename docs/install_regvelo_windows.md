# Install RegVelo for the human-limb project

This is the Windows environment used for the completed Person 2 fits. See the
[fitting guide](p2_model_fitting.md) for the completed results and release downloads.
Installing or verifying the environment does not start training.

## Install

Use PowerShell in the `BioModelSimProject` repository root. The tested interpreter
was 64-bit Python 3.12.6, with an NVIDIA GeForce RTX 4090 Laptop GPU. Use a project
environment at `.venv`; the separate class environment is not modified. On a fresh
checkout, install Python 3.12 and create the environment if it does not exist.
Download and extract the release's inputs archive before the data checks below.

```powershell
if (-not (Test-Path '.\.venv\Scripts\python.exe')) { py -3.12 -m venv .venv }
$python = (Resolve-Path '.\.venv\Scripts\python.exe').Path
& $python -c "import sys; print(sys.executable); assert sys.version_info[:2] == (3, 12) and sys.prefix != sys.base_prefix"
& $python -m pip install --upgrade pip
& $python -m pip install "torch==2.5.1" --index-url https://download.pytorch.org/whl/cu124
& $python -m pip install -r requirements-p2.txt
& $python -m pip check
```

Run each command after the previous one succeeds. The CUDA wheel index is for the PyTorch command only. Use the normal PyPI index for the project requirements. If a download fails with `getaddrinfo failed`, retry from a terminal with network access. For a resolver or import error, use the first failing dependency to diagnose it, apply the smallest compatible fix in this environment, and repeat the checks.

## Verify GPU and data access

Run the following from the same repository root. It uses the fitting scriptâ€™s input loader to check the rebuilt baseline against Aarioâ€™s cell manifest. It also executes a CUDA operation and constructs the projectâ€™s RegVelo model on the GPU.

```powershell
@'
import sys
from importlib.metadata import version

import cellrank
import scvelo
import torch
from regvelo import REGVELOVI

sys.path.insert(0, "scripts")
from p2_fit import load_config, load_input

print(sys.executable)
for package in ("regvelo", "scvi-tools", "torch", "scvelo", "cellrank"):
    print(package, version(package))
assert torch.version.cuda == "12.4" and torch.cuda.is_available()
assert torch.arange(4, device="cuda").square().sum().item() == 14
print(torch.cuda.get_device_name(0))

cfg = load_config("config/p2_x.yaml")
adata = load_input(cfg, "full", 0, "X")
assert adata.shape == (12207, 606)
REGVELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
prior = torch.tensor(adata.uns["skeleton"].to_numpy().T, dtype=torch.int32)
model = REGVELOVI(adata, W=prior, regulators=adata.var_names[adata.var["TF"]],
                  soft_constraint=True, lam2=cfg["model"]["lam2"])
model.to_device("cuda:0")
assert model.device.type == "cuda"
print("Loaded full_r0:", adata.shape, "Model device:", model.device)
'@ | & $python -
```

In VS Code, select `.venv\Scripts\python.exe`. All planned fits are already complete;
downstream users can load the released outputs without retraining. The
[model fitting guide](p2_model_fitting.md) covers reproduction commands and output
handoffs. `python scripts/p2_run_all.py --dry-run` displays the selected scope
without starting a fit.
