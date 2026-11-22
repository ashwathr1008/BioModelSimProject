"""Run the tutorial and configured PCA experiments sequentially on the project GPU."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results" / "p2"


def timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def configured_pcas(cfg):
    """Old configurations retain both variants; the selected X plan specifies one."""
    variants = cfg.get("execution", {}).get("pca_variants", ["X", "shipped_pca"])
    if (not isinstance(variants, list) or not variants
            or any(pca not in ("X", "shipped_pca") for pca in variants)
            or len(set(variants)) != len(variants)):
        raise ValueError("execution.pca_variants must list unique values from X and shipped_pca")
    return variants


def make_jobs(n_repeats, pca_variants=("X", "shipped_pca")):
    """Put selected seed-zero baselines first, then finish each repeat block."""
    jobs = [("tutorial", "tutorial", 0, None)]
    jobs.extend((f"{pca}/full_r0", "full", 0, pca) for pca in pca_variants)
    for repeat in range(n_repeats):
        for condition in ("full", "targeted", "random"):
            if repeat == 0 and condition == "full":
                continue
            for pca in pca_variants:
                jobs.append((f"{pca}/{condition}_r{repeat}", condition, repeat, pca))
    return jobs


def write_queue(record):
    record["updated_at"] = timestamp()
    temporary = RESULTS / "queue.json.tmp"
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(RESULTS / "queue.json")


def run_command(command, log_path):
    """Append output and terminate our own child if the queue is interrupted."""
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{timestamp()}] {subprocess.list2cmdline(command)}\n")
        log.flush()
        process = subprocess.Popen(
            command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
        )
        try:
            return process.wait()
        except BaseException:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise


def confirm_outputs(job_id):
    folder = RESULTS / job_id
    record = json.loads((folder / "run.json").read_text(encoding="utf-8"))
    if record.get("status") != "complete":
        raise RuntimeError(f"{job_id}: fitting process exited without a complete run record")
    for artifact in ("model/model.pt", "fitted.h5ad", "training.csv"):
        if not (folder / artifact).is_file():
            raise FileNotFoundError(f"{job_id}: missing {artifact}")


def confirm_config(config, expected_sha256):
    if hashlib.sha256(config.read_bytes()).hexdigest() != expected_sha256:
        raise RuntimeError(f"Configuration changed after the queue started: {config}; stop to avoid mixed settings")


def run_queue(config, n_repeats, pca_variants=("X", "shipped_pca")):
    jobs = make_jobs(n_repeats, pca_variants)
    logs = RESULTS / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    record = {
        "status": "running", "started_at": timestamp(), "config": str(config),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "pca_variants": list(pca_variants),
        "interpreter": sys.executable, "pid": os.getpid(), "current_job": None,
        "phase": "fitting", "completed_jobs": [], "failed_jobs": [],
        "jobs": {job[0]: {"status": "planned"} for job in jobs},
    }
    write_queue(record)
    try:
        for job_id, condition, repeat, pca in jobs:
            log_path = logs / f"{job_id.replace('/', '_')}.log"
            record.update(current_job=job_id, phase="fitting")
            job_record = record["jobs"][job_id]
            job_record.update(status="running", started_at=timestamp(), log=str(log_path))
            write_queue(record)
            print(f"[{timestamp()}] Fitting {job_id}; log: {log_path}", flush=True)
            command = [sys.executable, "-u", str(REPO / "scripts/p2_fit.py"), condition,
                       "--config", str(config)]
            if pca is not None:
                command.extend(["--repeat", str(repeat), "--pca-input", pca])
            confirm_config(config, record["config_sha256"])
            code = run_command(command, log_path)
            job_record["returncode"] = code
            if code:
                raise RuntimeError(f"{job_id}: fitting exited with code {code}; see {log_path}")
            confirm_outputs(job_id)
            job_record.update(status="complete", finished_at=timestamp())
            record["completed_jobs"].append(job_id)
            record["phase"] = "comparison"
            write_queue(record)
            confirm_config(config, record["config_sha256"])
            code = run_command(
                [sys.executable, "-u", str(REPO / "scripts/p2_compare.py"), "--config", str(config)],
                logs / "comparison.log",
            )
            if code:
                raise RuntimeError(f"Comparison exited with code {code}; see {logs / 'comparison.log'}")
            print(f"[{timestamp()}] Complete: {job_id} ({len(record['completed_jobs'])}/{len(jobs)})", flush=True)
        record.update(status="complete", current_job=None, phase="complete", finished_at=timestamp())
    except BaseException as error:
        interrupted = isinstance(error, KeyboardInterrupt)
        record.update(status="interrupted" if interrupted else "failed",
                      error=f"{type(error).__name__}: {error}", finished_at=timestamp())
        job_id = record["current_job"]
        if job_id and record["jobs"][job_id]["status"] == "running":
            record["jobs"][job_id].update(status=record["status"], finished_at=timestamp())
            record["failed_jobs"].append(job_id)
        raise
    finally:
        write_queue(record)


def main():
    import yaml
    from filelock import FileLock

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=REPO / "config/p2_x.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the selected jobs and saved status without starting fits or changing files")
    args = parser.parse_args()
    config = args.config.resolve()
    cfg = yaml.safe_load(config.read_text(encoding="utf-8"))
    try:
        pca_variants = configured_pcas(cfg)
    except ValueError as error:
        parser.error(str(error))
    if cfg["sampling"]["conditions"] != ["full", "targeted", "random"]:
        parser.error("this queue expects the full, targeted, random experiment")
    n_repeats = cfg["sampling"]["n_repeats"]
    if n_repeats < 1 or len(cfg["model"]["training_seeds"]) != n_repeats:
        parser.error("each repeat must have exactly one configured training seed")
    for filename in ("p2_fit.py", "p2_compare.py"):
        if not (REPO / "scripts" / filename).is_file():
            parser.error(f"missing scripts/{filename}")
    if args.dry_run:
        jobs = make_jobs(n_repeats, pca_variants)
        completed = 0
        for job_id, _, _, _ in jobs:
            record_path = RESULTS / job_id / "run.json"
            status = "planned"
            if record_path.is_file():
                status = json.loads(record_path.read_text(encoding="utf-8")).get("status", "unknown")
            completed += status == "complete"
            print(f"{job_id}\t{status}")
        print(f"{len(jobs)} selected runs; {completed} saved complete; {len(jobs) - completed} not complete.")
        print("Completed runs will be validated against inputs and settings before reuse. No fits started.")
        return
    RESULTS.mkdir(parents=True, exist_ok=True)
    with FileLock(RESULTS / ".queue.lock", timeout=0):
        run_queue(config, n_repeats, pca_variants)


if __name__ == "__main__":
    main()
