"""Fit the selected human-limb RegVelo model to one dataset."""

import argparse
import json
import platform
import sys
import time
import warnings
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from p1_data.common import condition_id, dense, load_config, load_source, repo_path, sha256


# Explicit RegVelo 1.0.0 defaults keep both PCA variants on the same training protocol.
TRAINING_OPTIONS = dict(
    optimizer="AdamW", lr=0.01, weight_decay=1e-5, eps=1e-16,
    train_size=0.9, validation_size=None, early_stopping=True,
    gradient_clip_val=10, early_stopping_patience=45,
    early_stopping_monitor="elbo_validation", accelerator="gpu", devices=[0],
)


def input_path(cfg, condition, repeat, pca):
    if condition == "tutorial":
        return repo_path(cfg["data"]["source_file"])
    folder = repo_path(cfg["output"]["prepared_dir"])
    if pca != cfg["preprocessing"]["pca_input"]:
        folder = folder.with_name(f"{folder.name}_{pca}")
    return folder / f"{condition_id(condition, repeat)}.h5ad"


def manifest_path(cfg, condition, repeat):
    return repo_path(cfg["output"]["results_dir"]) / "manifests" / f"{condition_id(condition, repeat)}.csv.gz"


def run_identity(cfg, condition, repeat, pca):
    """A completed fit is reusable only for the same inputs and scientific settings."""
    tutorial = condition == "tutorial"
    return dict(
        condition=condition, repeat=repeat, pca_input=pca,
        input_sha256=sha256(input_path(cfg, condition, repeat, pca)),
        manifest_sha256=None if tutorial else sha256(manifest_path(cfg, condition, repeat)),
        model_settings=cfg["model"], training_options=TRAINING_OPTIONS,
        preprocessing=None if tutorial else {**cfg["preprocessing"], "pca_input": pca},
        sampling=None if tutorial else cfg["sampling"],
        endpoints=cfg["endpoints"], stage_key=cfg["data"]["stage_key"],
        ode_implementation_sha256=(sha256(Path(__file__).with_name("p2_memory.py"))
                                   if cfg["model"].get("ode_chunk_size") else None),
    )


def load_input(cfg, condition, repeat, pca):
    """Load the tutorial data or check a rebuilt dataset against its sampling manifest."""
    if condition == "tutorial":
        return load_source(cfg)

    import anndata as ad

    name = condition_id(condition, repeat)
    adata = ad.read_h5ad(input_path(cfg, condition, repeat, pca))
    prep = json.loads(adata.uns["p1_prep_json"])
    expected = {
        "condition": condition, "repeat": repeat,
        "sampling_seed": cfg["sampling"]["sampling_seed_base"] + repeat,
        "training_seed": cfg["model"]["training_seeds"][repeat],
    }
    for key, value in expected.items():
        if prep[key] != value or not adata.obs[key].eq(value).all():
            raise ValueError(f"{name}: inconsistent {key}, expected {value!r}")
    expected_preprocessing = {**cfg["preprocessing"], "pca_input": pca}
    if prep["preprocessing"] != expected_preprocessing:
        raise ValueError(f"{name}: preprocessing differs from the requested configuration")
    for key in ("sampling", "endpoints"):
        if prep[key] != cfg[key]:
            raise ValueError(f"{name}: {key} differs from the requested configuration")

    manifest = manifest_path(cfg, condition, repeat)
    cells = pd.read_csv(manifest, usecols=["cell_id", "status"])
    retained = cells.loc[cells["status"] == "retained", "cell_id"]
    if not adata.obs_names.equals(pd.Index(retained)):
        raise ValueError(f"{name}: cell IDs or order differ from {manifest.name}")
    return adata


def fit(adata, cfg, record, folder):
    """Train RegVelo and save the checkpoint and cell-level outputs together."""
    import regvelo as rgv
    import scvi
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for these GPU fits; check the project interpreter and PyTorch installation")

    genes = adata.var_names
    prior = adata.uns["skeleton"]
    if not (genes.is_unique and genes.equals(prior.index) and genes.equals(prior.columns)):
        raise ValueError("Regulatory prior rows and columns must match the ordered gene panel")
    regulators = genes[adata.var["TF"]]
    if not {"MSC", "MYOG"}.issubset(regulators):
        raise ValueError("Expected both MSC and MYOG in the transcription factor panel")
    for layer in ("Ms", "Mu"):
        values = np.asarray(dense(adata.layers[layer]), dtype=np.float32, order="C")
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite values in {layer}")
        adata.layers[layer] = values

    settings = cfg["model"]
    scvi.settings.seed = record["training_seed"]
    rgv.REGVELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    # The stored prior has regulators in rows. RegVelo expects them in columns.
    weights = torch.tensor(prior.to_numpy().T, dtype=torch.int32)
    model = rgv.REGVELOVI(
        adata, W=weights, regulators=regulators, soft_constraint=True, lam2=settings["lam2"]
    )
    model.to_device("cuda:0")
    if settings.get("ode_chunk_size"):
        from p2_memory import enable_ode_checkpointing

        record["memory_optimization"] = enable_ode_checkpointing(model, settings["ode_chunk_size"])
    torch.cuda.reset_peak_memory_stats(model.device)
    record.update(n_cells=adata.n_obs, n_genes=adata.n_vars, batch_size=adata.n_obs)
    model.train(max_epochs=settings["max_epochs"], batch_size=adata.n_obs, **TRAINING_OPTIONS)
    record["device"] = str(model.device)
    if model.device.type != "cuda":
        raise RuntimeError(f"Expected CUDA training, got {model.device}")
    record["gpu"] = torch.cuda.get_device_name(model.device)
    train_loss = model.history["elbo_train"]
    validation_loss = model.history["elbo_validation"]
    history = pd.concat([train_loss, validation_loss], axis=1)
    history.to_csv(folder / "training.csv", index_label="epoch")
    record["epochs"] = len(train_loss)
    model.save(str(folder / "model"), overwrite=True)

    rgv.tl.set_output(adata, model, n_samples=settings["posterior_samples"], batch_size=adata.n_obs)
    for layer in ("velocity", "fit_t", "latent_time_velovi"):
        if not np.isfinite(adata.layers[layer]).all():
            raise ValueError(f"Nonfinite model output in {layer}")

    # The tutorial averages the gene-scaled times before rescaling across cells.
    cell_time = np.asarray(adata.layers["fit_t"]).mean(axis=1)
    stage = adata.obs[cfg["data"]["stage_key"]].astype(str).str.removeprefix("Pcw").astype(float)
    if np.ptp(cell_time) == 0 or not np.isfinite(stage).all() or stage.nunique() < 2:
        raise ValueError("Stage correlation needs finite, varying times and stages")
    adata.obs["latent_time"] = (cell_time - cell_time.min()) / np.ptp(cell_time)
    adata.obs["adj_stage_num"] = stage
    record["stage_spearman"] = float(spearmanr(adata.obs["latent_time"], stage).statistic)
    if model.device.type == "cuda":
        record["gpu_peak_allocated_mib"] = round(torch.cuda.max_memory_allocated(model.device) / 2**20)
    adata.write_h5ad(folder / "fitted.h5ad", compression="gzip")


def write_record(folder, record):
    path = folder / "run.json.tmp"
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    path.replace(folder / "run.json")


def run_status(folder):
    path = folder / "run.json"
    if not path.exists():
        return "planned"
    status = json.loads(path.read_text(encoding="utf-8"))["status"]
    if status == "complete":
        for artifact in ("model/model.pt", "fitted.h5ad", "training.csv"):
            if not (folder / artifact).is_file():
                return "missing local outputs"
    return status


def peak_memory_mib():
    if sys.platform == "win32":
        try:
            import psutil
        except ImportError:
            return None
        return round(psutil.Process().memory_info().peak_wset / 2**20)
    if sys.platform.startswith("linux"):
        import resource
        return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    return None


def run(cfg, condition, repeat, pca, folder):
    from filelock import FileLock

    folder.mkdir(parents=True, exist_ok=True)
    with FileLock(folder / ".fit.lock", timeout=0):
        if run_status(folder) == "complete":
            previous = json.loads((folder / "run.json").read_text(encoding="utf-8"))
            if previous.get("identity") != run_identity(cfg, condition, repeat, pca):
                raise ValueError(f"Completed fit has different or unverified inputs/settings: {folder}")
            print(f"Already complete: {folder}")
            return
        started = time.perf_counter()
        record = dict(
            status="running", condition=condition, repeat=repeat, pca_input=pca,
            training_seed=cfg["model"]["training_seeds"][repeat],
            sampling_seed=None if condition == "tutorial" else cfg["sampling"]["sampling_seed_base"] + repeat,
            model="soft_regularized", settings=cfg["model"],
            training_options=TRAINING_OPTIONS, interpreter=sys.executable,
            python=platform.python_version(), platform=platform.platform(),
            cpu=platform.processor() or platform.machine(),
        )
        write_record(folder, record)
        packages = ("regvelo", "scvi-tools", "torch", "lightning", "torchode", "anndata",
                    "scanpy", "scvelo", "cellrank", "numpy", "scipy", "pandas")
        with warnings.catch_warnings(record=True) as caught:
            try:
                record["packages"] = {name: version(name) for name in packages}
                adata = load_input(cfg, condition, repeat, pca)
                record["identity"] = run_identity(cfg, condition, repeat, pca)
                fit(adata, cfg, record, folder)
                record["status"] = "complete"
            except BaseException as error:
                record.update(status="failed", error=f"{type(error).__name__}: {error}")
                raise
            finally:
                record["elapsed_seconds"] = round(time.perf_counter() - started, 1)
                record["warnings"] = sorted({str(w.message) for w in caught})
                record["peak_cpu_memory_mib"] = peak_memory_mib()
                write_record(folder, record)
    print(f"Saved {folder} (stage Spearman = {record['stage_spearman']:.3f})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("condition", help="tutorial, status, or a condition from the configuration")
    parser.add_argument("--config", help="Path to an experiment YAML; defaults to config/experiment.yaml")
    parser.add_argument("--repeat", type=int, default=0)
    parser.add_argument("--pca-input", choices=["X", "shipped_pca"])
    args = parser.parse_args()
    cfg = load_config(args.config)
    conditions = cfg["sampling"]["conditions"]
    if args.condition not in ["tutorial", *conditions, "status"]:
        parser.error(f"unknown condition {args.condition!r}; configured conditions: {conditions}")
    if args.repeat not in range(len(cfg["model"]["training_seeds"])):
        parser.error("repeat must index the configured training seeds")
    pca = args.pca_input or cfg["preprocessing"]["pca_input"]
    root = repo_path("results/p2")
    if args.condition == "status":
        print("run\tstatus\toutput")
        folders = [root / "tutorial"]
        for repeat in range(cfg["sampling"]["n_repeats"]):
            folders.extend(root / pca / condition_id(c, repeat) for c in conditions)
        for folder in folders:
            print(f"{folder.name}\t{run_status(folder)}\t{folder.relative_to(repo_path('.'))}")
    elif args.condition == "tutorial":
        if args.repeat != 0 or args.pca_input:
            parser.error("the shipped tutorial pilot uses repeat 0 and its original PCA")
        run(cfg, "tutorial", 0, "shipped", root / "tutorial")
    else:
        folder = root / pca / condition_id(args.condition, args.repeat)
        run(cfg, args.condition, args.repeat, pca, folder)


if __name__ == "__main__":
    main()
