"""Run exactly X/full_r0 after the current tutorial handoff finishes, then stop."""

import ctypes
import hashlib
import json
import os
import subprocess
import sys
import time
from ctypes import wintypes
from datetime import datetime, timezone
from pathlib import Path

import psutil
from filelock import FileLock

from p2_finish_current import verify_completed_fit


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "p2"
CONFIG = ROOT / "config" / "p2_both_pca.yaml"
JOB = "X/full_r0"


def timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_status(record):
    record["updated_at"] = timestamp()
    path = RESULTS / "local_followup.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def check_config(expected_hash):
    if hashlib.sha256(CONFIG.read_bytes()).hexdigest() != expected_hash:
        raise RuntimeError("Frozen configuration changed; refusing to mix protocols")


def execute(command, log_path, record):
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{timestamp()}] {subprocess.list2cmdline(command)}\n")
        log.flush()
        process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        record["child_pid"] = process.pid
        write_status(record)
        try:
            while process.poll() is None:
                time.sleep(10)
                write_status(record)
            record["child_returncode"] = process.returncode
            record["child_pid"] = None
            if process.returncode:
                raise RuntimeError(f"Command failed with code {process.returncode}; see {log_path}")
        except BaseException:
            if process.poll() is None:
                parent = psutil.Process(process.pid)
                children = parent.children(recursive=True)
                for child in reversed(children):
                    try:
                        child.terminate()
                    except psutil.NoSuchProcess:
                        pass
                process.terminate()
                process.wait(timeout=15)
            raise


def main():
    if os.name != "nt":
        raise RuntimeError("This local follow-up runner is for the Windows laptop")
    with FileLock(RESULTS / ".local_followup.lock", timeout=0):
        config_hash = hashlib.sha256(CONFIG.read_bytes()).hexdigest()
        api = ctypes.WinDLL("kernel32", use_last_error=True).SetThreadExecutionState
        api.argtypes = [wintypes.DWORD]
        api.restype = wintypes.DWORD
        previous_state = api(0x80000001)
        if not previous_state:
            raise ctypes.WinError(ctypes.get_last_error())
        record = dict(status="waiting_for_tutorial", pid=os.getpid(), job=JOB,
                      started_at=timestamp(), config=str(CONFIG), config_sha256=config_hash,
                      maximum_additional_fits=1, keep_awake_active=True,
                      keep_display_on=False, execution_state_api_return=int(previous_state))
        try:
            write_status(record)
            print("Waiting for the tutorial handoff; then exactly X/full_r0 will run locally.", flush=True)
            while True:
                check_config(config_hash)
                handoff = json.loads((RESULTS / "finish_current.json").read_text(encoding="utf-8"))
                if handoff["status"] == "finished":
                    if not handoff.get("tutorial_complete"):
                        raise RuntimeError("Tutorial did not complete successfully; follow-up will not start")
                    if handoff.get("comparison_returncode") != 0:
                        raise RuntimeError("Tutorial comparison refresh failed; follow-up will not start")
                    break
                if handoff["status"] == "error":
                    raise RuntimeError(f"Tutorial handoff failed: {handoff.get('error')}")
                try:
                    watcher = psutil.Process(handoff["watcher_pid"])
                    if not any(Path(argument).name == "p2_finish_current.py" for argument in watcher.cmdline()):
                        raise RuntimeError("Tutorial watcher process identity changed")
                except psutil.NoSuchProcess:
                    # The watcher may publish completion and exit between these reads.
                    latest = json.loads((RESULTS / "finish_current.json").read_text(encoding="utf-8"))
                    if latest["status"] == "finished":
                        continue
                    raise RuntimeError("Tutorial watcher exited before completing the handoff") from None
                time.sleep(10)
                write_status(record)

            complete, error = verify_completed_fit(RESULTS / "tutorial")
            if not complete:
                raise RuntimeError(error)
            queue = json.loads((RESULTS / "queue.json").read_text(encoding="utf-8"))
            if queue["status"] != "stopped_after_current":
                raise RuntimeError("Original local queue has not stopped safely")
            check_config(config_hash)
            record.update(status="running", phase="fitting", fit_started_at=timestamp())
            write_status(record)
            execute([sys.executable, "-u", str(ROOT / "scripts/p2_fit.py"), "full",
                     "--repeat", "0", "--pca-input", "X", "--config", str(CONFIG)],
                    RESULTS / "logs" / "X_full_r0.log", record)
            complete, error = verify_completed_fit(RESULTS / JOB)
            if not complete:
                raise RuntimeError(error)
            record.update(phase="comparison", fit_complete=True)
            write_status(record)
            check_config(config_hash)
            execute([sys.executable, "-u", str(ROOT / "scripts/p2_compare.py"),
                     "--config", str(CONFIG)], RESULTS / "logs" / "comparison.log", record)
            record.update(status="complete", phase="stopped", finished_at=timestamp())
            print("X/full_r0 is complete. No further local fits will be launched.", flush=True)
        except BaseException as error:
            record.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                          error=f"{type(error).__name__}: {error}", finished_at=timestamp())
            raise
        finally:
            api(0x80000000)
            record["keep_awake_active"] = False
            write_status(record)


if __name__ == "__main__":
    main()
