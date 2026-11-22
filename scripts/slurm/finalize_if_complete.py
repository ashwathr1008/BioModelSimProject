"""Let the last successful array task write the shared comparison exactly once."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from filelock import FileLock

from validate_outputs import load_config, planned_runs, repo_path, validate_run


def finalize(config, expected_sha256, expected_gpu=""):
    config = Path(config).resolve()
    actual_sha256 = hashlib.sha256(config.read_bytes()).hexdigest()
    if actual_sha256 != expected_sha256:
        raise ValueError("Configuration changed before comparison")
    cfg = load_config(config)
    receipt_dir = repo_path("results/p2/slurm")
    receipt_dir.mkdir(parents=True, exist_ok=True)
    marker = receipt_dir / "comparison_complete.json"
    with FileLock(receipt_dir / ".comparison.lock"):
        runs = planned_runs(cfg)
        for pca, condition, repeat in runs:
            folder = repo_path("results/p2") / pca / f"{condition}_r{repeat}"
            record_path = folder / "run.json"
            if not record_path.is_file():
                print("Other fits remain; comparison deferred", flush=True)
                return False
            record = json.loads(record_path.read_text(encoding="utf-8"))
            if record.get("status") != "complete":
                print("Other fits remain; comparison deferred", flush=True)
                return False
            for name in ("model/model.pt", "fitted.h5ad", "training.csv"):
                if not (folder / name).is_file():
                    raise ValueError(f"Completed fit is missing an artifact: {folder / name}")
        if marker.is_file():
            previous = json.loads(marker.read_text(encoding="utf-8"))
            if previous.get("config_sha256") != expected_sha256:
                raise ValueError("Comparison marker belongs to a different configuration")
            import_hash = (os.environ.get("P2_IMPORT_COMPLETED_SHA256", "").lower()
                           if os.environ.get("P2_REUSE_LOCAL_FULL0", "0") == "1" else None)
            if previous.get("import_receipt_sha256") != import_hash:
                raise ValueError("Comparison marker belongs to different imported-fit provenance")
            print(f"Comparison already complete: {marker}", flush=True)
            return True
        checked = [validate_run(cfg, *run, expected_gpu) for run in runs]
        subprocess.run(
            [sys.executable, str(repo_path("scripts/p2_compare.py")), "--config", str(config)],
            cwd=repo_path("."), check=True,
        )
        record = {"status": "complete", "config_sha256": expected_sha256,
                  "completed_utc": datetime.now(timezone.utc).isoformat(),
                  "validated_runs": checked, "output_dir": "results/p2/comparison",
                  "mixed_hardware": len({run["gpu"] for run in checked}) > 1,
                  "import_receipt_sha256": next((run["imported_local"]["receipt_sha256"]
                                                 for run in checked if run.get("imported_local")), None)}
        temporary = marker.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        temporary.replace(marker)
        print(f"Comparison complete: {marker}", flush=True)
        return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--expected-gpu", default="")
    args = parser.parse_args()
    finalize(args.config, args.config_sha256, args.expected_gpu)


if __name__ == "__main__":
    main()
