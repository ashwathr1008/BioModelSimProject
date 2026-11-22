"""Finish the active fit, then retire its suspended local queue supervisor."""

import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import psutil
from filelock import FileLock


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "p2"


def timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_json(path, record):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def same_process(pid, created):
    try:
        process = psutil.Process(pid)
        return process if process.create_time() == created and process.is_running() else None
    except psutil.NoSuchProcess:
        return None


def verify_completed_fit(folder):
    path = folder / "run.json"
    if not path.is_file():
        return False, "Fit process ended without run.json"
    run = json.loads(path.read_text(encoding="utf-8"))
    if run.get("status") != "complete":
        return False, run.get("error") or f"Fit process ended with status {run.get('status')!r}"
    missing = [name for name in ("model/model.pt", "fitted.h5ad", "training.csv")
               if not (folder / name).is_file() or (folder / name).stat().st_size == 0]
    return (False, "Missing completed-fit outputs: " + ", ".join(missing)) if missing else (True, None)


def main():
    status_path = RESULTS / "finish_current.json"
    queue_path = RESULTS / "queue.json"
    with FileLock(RESULTS / ".finish_current.lock", timeout=0):
        queue = json.loads(queue_path.read_text(encoding="utf-8"))
        if not (queue.get("stop_after_current") and queue.get("supervisor_suspended")):
            raise RuntimeError("The local queue must already be suspended for finish-current mode")
        if queue.get("current_job") != "tutorial":
            raise RuntimeError("This handoff watcher expects the active tutorial fit")
        parent = psutil.Process(queue["pid"])
        if not any(Path(argument).name == "p2_run_all.py" for argument in parent.cmdline()):
            raise RuntimeError("Supervisor command does not match the local queue")
        if parent.status() != psutil.STATUS_STOPPED:
            raise RuntimeError("Supervisor is not suspended; refusing to interfere with running jobs")
        parent_created = parent.create_time()
        identities = queue["active_fit_processes"]
        if not identities:
            raise RuntimeError("No recorded active fit processes")
        for item in identities:
            process = same_process(item["pid"], item["created"])
            if process is not None:
                command = process.cmdline()
                if not (any(Path(argument).name == "p2_fit.py" for argument in command)
                        and "tutorial" in command):
                    raise RuntimeError("A recorded fit process has an unexpected command")
        record = dict(status="waiting_for_current_fit", watcher_pid=os.getpid(),
                      supervisor_pid=parent.pid, supervisor_created=parent_created,
                      active_fit_processes=identities, current_job="tutorial", started_at=timestamp())
        write_json(status_path, record)
        print("Watching tutorial; its training process remains running and the supervisor stays suspended.", flush=True)
        try:
            while any(same_process(item["pid"], item["created"]) for item in identities):
                supervisor = same_process(parent.pid, parent_created)
                if supervisor is not None and supervisor.status() != psutil.STATUS_STOPPED:
                    raise RuntimeError("Queue supervisor was resumed unexpectedly")
                record["updated_at"] = timestamp()
                write_json(status_path, record)
                time.sleep(10)

            complete, error = verify_completed_fit(RESULTS / "tutorial")
            supervisor = same_process(parent.pid, parent_created)
            if supervisor is not None:
                if supervisor.status() != psutil.STATUS_STOPPED:
                    raise RuntimeError("Refusing to terminate a resumed queue supervisor")
                supervisor.terminate()
                supervisor.wait(timeout=15)

            queue = json.loads(queue_path.read_text(encoding="utf-8"))
            if queue["pid"] != parent.pid or queue.get("current_job") != "tutorial":
                raise RuntimeError("Queue identity changed during the finish-current handoff")
            queue.update(status="stopped_after_current" if complete else "failed",
                         phase="stopped", current_job=None, supervisor_suspended=False,
                         finished_at=timestamp(), updated_at=timestamp(),
                         remaining_destination="NCSA Delta")
            queue["jobs"]["tutorial"].update(status="complete" if complete else "failed", finished_at=timestamp())
            collection = queue["completed_jobs"] if complete else queue["failed_jobs"]
            if "tutorial" not in collection:
                collection.append("tutorial")
            if error:
                queue["error"] = error
            write_json(queue_path, queue)
            with (RESULTS / "logs" / "comparison.log").open("a", encoding="utf-8") as log:
                log.write(f"\n[{timestamp()}] Refresh after finishing local tutorial\n")
                log.flush()
                result = subprocess.run(
                    [sys.executable, "-u", str(ROOT / "scripts/p2_compare.py"),
                     "--config", queue["config"]], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
                )
            record.update(status="finished", tutorial_complete=complete,
                          comparison_returncode=result.returncode, finished_at=timestamp(), error=error)
            write_json(status_path, record)
            print(f"Tutorial complete={complete}; local queue stopped. Remaining fits will not start locally.", flush=True)
        except BaseException as error:
            record.update(status="error", error=f"{type(error).__name__}: {error}", updated_at=timestamp())
            write_json(status_path, record)
            raise


if __name__ == "__main__":
    main()
