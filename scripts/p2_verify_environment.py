"""Verify the GPU environment and every prepared input before starting the fit queue."""

import argparse
import json
import platform
import subprocess
import sys
from importlib.metadata import version

from p2_fit import load_config, load_input, repo_path
from p2_run_all import configured_pcas
from p1_data.rebuild import check_condition


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/p2_x.yaml")
    parser.add_argument("--inputs-only", action="store_true", help="Validate inputs without opening another GPU context")
    args = parser.parse_args()
    cfg = load_config(args.config)
    pca_variants = configured_pcas(cfg)
    subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
    source = load_input(cfg, "tutorial", 0, "shipped")
    assert source.shape == (12207, 606)
    endpoints = cfg["endpoints"]
    celltype_key = cfg["data"]["celltype_key"]
    expected_endpoints = {name: int(source.obs[celltype_key].eq(name).sum()) for name in endpoints}
    checks_by_run = {}
    for pca in pca_variants:
        for repeat in range(cfg["sampling"]["n_repeats"]):
            for condition in cfg["sampling"]["conditions"]:
                adata = load_input(cfg, condition, repeat, pca)
                checks, _ = check_condition(adata, source, cfg)
                assert all(checks[f"{layer}_all_finite"] for layer in ("Ms", "Mu", "spliced", "unspliced"))
                assert checks["Ms_n_constant_genes"] == checks["Mu_n_constant_genes"] == 0
                assert checks["gene_order_matches_source"] and checks["skeleton_identical_to_source"]
                assert checks["skeleton_aligned_to_var"] and checks["graph_n_components"] == 1
                assert checks["endpoint_counts"] == expected_endpoints
                expected_shape = (12207 if condition == "full" else 11638, 606)
                assert adata.shape == expected_shape
                name = f"{pca}/{condition}_r{repeat}"
                checks_by_run[name] = {"shape": list(adata.shape), **checks}
                print(f"Verified {name}: {adata.shape}", flush=True)
    folder = repo_path("results/p2")
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "input_validation.json").write_text(json.dumps(checks_by_run, indent=2) + "\n", encoding="utf-8")
    if args.inputs_only:
        print(f"All {len(checks_by_run)} selected prepared inputs verified; GPU check omitted by request.")
        return
    import cellrank
    import scvelo
    import torch
    from regvelo import REGVELOVI

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable")
    assert torch.arange(4, device="cuda").square().sum().item() == 14
    adata = load_input(cfg, "full", 0, pca_variants[0])
    REGVELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    prior = torch.tensor(adata.uns["skeleton"].to_numpy().T, dtype=torch.int32)
    model = REGVELOVI(
        adata, W=prior, regulators=adata.var_names[adata.var["TF"]],
        soft_constraint=True, lam2=cfg["model"]["lam2"],
    )
    model.to_device("cuda:0")
    assert model.device.type == "cuda"
    packages = ("regvelo", "scvi-tools", "torch", "torchode", "lightning", "scvelo", "cellrank",
                "anndata", "scanpy", "numpy", "scipy", "pandas", "scikit-learn", "jax", "jaxlib")
    record = {
        "status": "verified", "python": platform.python_version(), "interpreter": sys.executable,
        "gpu": torch.cuda.get_device_name(0), "cuda": torch.version.cuda,
        "packages": {name: version(name) for name in packages}, "prepared_inputs": checks_by_run,
        "pca_variants": pca_variants,
    }
    (folder / "environment.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    freeze = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (folder / "environment-lock.txt").write_text(freeze, encoding="utf-8")
    print(f"GPU model initialized on {record['gpu']}; all {len(checks_by_run)} selected prepared inputs verified.")


if __name__ == "__main__":
    main()
