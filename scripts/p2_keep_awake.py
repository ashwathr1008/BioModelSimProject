"""Prevent idle sleep only while the existing RegVelo run queue is alive."""

import ctypes
import json
import os
import time
from ctypes import wintypes
from datetime import datetime, timezone
from pathlib import Path

import psutil
from filelock import FileLock


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "p2"
ES_CONTINUOUS = 0x80000000
ES_SYSTEM_REQUIRED = 0x00000001


def write_status(record):
    record["updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    path = RESULTS / "keep_awake.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main():
    if os.name != "nt":
        raise RuntimeError("This sleep-prevention helper requires Windows")
    with FileLock(RESULTS / ".keep_awake.lock", timeout=0):
        queue = json.loads((RESULTS / "queue.json").read_text(encoding="utf-8"))
        if queue["status"] != "running":
            raise RuntimeError("The training queue is not running")
        process = psutil.Process(queue["pid"])
        command = process.cmdline()
        if not any(Path(argument).name == "p2_run_all.py" for argument in command):
            raise RuntimeError("Queue PID does not belong to the expected training script")
        started = process.create_time()
        api = ctypes.WinDLL("kernel32", use_last_error=True).SetThreadExecutionState
        api.argtypes = [wintypes.DWORD]
        api.restype = wintypes.DWORD
        if not api(ES_CONTINUOUS | ES_SYSTEM_REQUIRED):
            raise ctypes.WinError(ctypes.get_last_error())
        record = dict(status="active", helper_pid=os.getpid(), queue_pid=process.pid,
                      queue_process_created=started, prevent_idle_sleep=True,
                      keep_display_on=False, persistent_settings_changed=False)
        try:
            write_status(record)
            print(f"Preventing idle sleep while training queue PID {process.pid} runs.", flush=True)
            while process.is_running() and process.create_time() == started:
                state = json.loads((RESULTS / "queue.json").read_text(encoding="utf-8"))
                if state.get("pid") != process.pid or state.get("status") != "running":
                    break
                time.sleep(10)
                write_status(record)
        finally:
            api(ES_CONTINUOUS)
            record.update(status="released", prevent_idle_sleep=False)
            write_status(record)
            print("Training queue ended; released the idle-sleep request.", flush=True)


if __name__ == "__main__":
    main()
