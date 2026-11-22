"""Stream progress over the existing authenticated connection; never manage jobs."""

import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/p2"
EXPECTED = [f"{pca}/{condition}_r{repeat}"
            for repeat in range(3)
            for condition in ("full", "targeted", "random")
            for pca in ("X", "shipped_pca")]


def snapshot():
    job_ids = {path.read_text().strip().split(";")[0]
               for path in (RESULTS / "slurm").glob("campaign.*/*.job_id")}
    if not job_ids or not all(job.isdigit() for job in job_ids):
        raise RuntimeError("Missing or invalid campaign job receipts")
    runs = {}
    for name in EXPECTED:
        path = RESULTS / name / "run.json"
        if not path.exists():
            runs[name] = {"status": "planned"}
            continue
        try:
            record = json.loads(path.read_text())
            runs[name] = {key: record[key] for key in
                          ("status", "epochs", "elapsed_seconds", "gpu", "stage_spearman", "error")
                          if key in record}
        except (OSError, ValueError) as error:
            runs[name] = {"status": "unreadable", "error": str(error)}
    queue = subprocess.run(
        ["squeue", "-u", "araghuram", "--noheader", "--format=%i|%j|%T|%M|%R"],
        capture_output=True, text=True, check=False,
    )
    own_queue = [line for line in queue.stdout.splitlines()
                 if line.split("|", 1)[0].split("_", 1)[0].split(".", 1)[0] in job_ids]
    marker = RESULTS / "slurm/comparison_complete.json"
    comparison = json.loads(marker.read_text()) if marker.exists() else {"status": "pending"}
    return {"updated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "project_dir": str(ROOT),
            "completed": sum(run["status"] == "complete" for run in runs.values()),
            "failed": sum(run["status"] == "failed" for run in runs.values()),
            "runs": runs, "job_ids": sorted(job_ids), "comparison": comparison,
            "queue_returncode": queue.returncode,
            "queue": own_queue, "queue_error": queue.stderr.strip()}


def main():
    status_path = RESULTS / "slurm" / "live_status.json"
    while True:
        status = snapshot()
        temporary = status_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(status, indent=2) + "\n")
        temporary.replace(status_path)
        print("P2_STATUS " + json.dumps(status, separators=(",", ":")), flush=True)
        if status["completed"] == len(EXPECTED) and status["comparison"].get("status") == "complete":
            print("All 18 fits and the comparison are complete.", flush=True)
            return
        if status["queue_returncode"] == 0 and not status["queue"]:
            print("No jobs remain in your queue. Inspect run records and Slurm logs for unfinished fits.", flush=True)
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
