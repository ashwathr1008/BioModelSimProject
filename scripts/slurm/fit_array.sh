#!/usr/bin/env bash
# Scheduler account, partition, CPUs, memory, GPU and time are supplied by sbatch.
# Exactly one independent model uses one GPU; no MPI or multi-GPU training.
set -euo pipefail

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${SLURM_ARRAY_TASK_ID:?An array index from 0 through 17 is required}"
: "${PROJECT_DIR:?Export the absolute project path}"
: "${P2_ENV:?Export the prepared Linux virtual environment path}"
: "${P2_CONFIG:?Export the frozen configuration snapshot path}"
: "${P2_CONFIG_SHA256:?Export the snapshot SHA-256}"
[[ "$SLURM_ARRAY_TASK_ID" =~ ^([0-9]|1[0-7])$ ]] || {
    echo 'Array task must be between 0 and 17' >&2; exit 2;
}
cd -- "$PROJECT_DIR"
[[ -f "$P2_ENV/.p2_environment_ready.json" ]] || {
    echo 'Run setup_environment.sh once before submitting fits' >&2; exit 2;
}
# Keep setup from mutating the shared environment while any fit is using it.
exec 9>"${P2_ENV}.setup.lock"
flock --shared 9
PYTHON="$P2_ENV/bin/python"
export PYTHONUNBUFFERED=1 PYTHONUTF8=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS" OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
"$PYTHON" - "$PROJECT_DIR" "$P2_ENV" "$P2_CONFIG" "$P2_CONFIG_SHA256" <<'PY'
import hashlib, json, math, os, shutil, sys
from pathlib import Path
import torch
project, environment, config = map(Path, sys.argv[1:4])
minimum_free_gib = float(os.environ.get("P2_MIN_FREE_GIB", "2"))
if not math.isfinite(minimum_free_gib) or minimum_free_gib <= 0:
    raise ValueError("P2_MIN_FREE_GIB must be a positive finite number")
free_bytes = shutil.disk_usage(project).free
if free_bytes < minimum_free_gib * 2**30:
    raise RuntimeError(f"Project filesystem has only {free_bytes / 2**30:.2f} GiB free; require {minimum_free_gib:g} GiB")
print(f"Project filesystem free: {free_bytes / 2**30:.2f} GiB", flush=True)
record = json.loads((environment / ".p2_environment_ready.json").read_text())
assert sys.version_info[:2] == (3, 12), sys.version
for name, expected in record["requirements_sha256"].items():
    assert hashlib.sha256((project / name).read_bytes()).hexdigest() == expected, name
assert hashlib.sha256(config.read_bytes()).hexdigest() == sys.argv[4], "Configuration changed"
assert torch.cuda.is_available() and torch.cuda.device_count() == 1, "Require one Slurm-assigned GPU"
assert torch.arange(4, device="cuda").square().sum().item() == 14
gpu = torch.cuda.get_device_name(0)
expected_gpu = os.environ.get("P2_EXPECT_GPU", "")
assert not expected_gpu or expected_gpu.lower() in gpu.lower(), (expected_gpu, gpu)
print("Allocated GPU:", gpu, flush=True)
PY
if [[ "${P2_REQUIRE_PILOTS:-0}" = 1 ]]; then
    "$PYTHON" scripts/slurm/validate_outputs.py --config "$P2_CONFIG" \
        --scope pilots --expected-gpu "${P2_EXPECT_GPU:-}"
fi

# Within each repeat: X/full, shipped/full, X/targeted, shipped/targeted,
# X/random, shipped/random. Tutorial is intentionally excluded.
task="$SLURM_ARRAY_TASK_ID"
repeat=$((task / 6))
condition_index=$(((task % 6) / 2))
pca_index=$((task % 2))
conditions=(full targeted random)
pcas=(X shipped_pca)
condition="${conditions[$condition_index]}"
pca="${pcas[$pca_index]}"
printf 'Job %s task %s: %s/%s_r%s on %s\n' "$SLURM_JOB_ID" "$task" "$pca" "$condition" "$repeat" "$(hostname)"
job_work_dir="$(mktemp -d "${TMPDIR:-/tmp}/regvelo-p2.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}.XXXXXXXX")"
printf 'Training working directory: %s\n' "$job_work_dir"
srun --ntasks=1 --chdir="$job_work_dir" "$PYTHON" -u "$PROJECT_DIR/scripts/p2_fit.py" "$condition" \
    --repeat "$repeat" --pca-input "$pca" --config "$P2_CONFIG"
if [[ "${P2_COMPARE_WHEN_COMPLETE:-0}" = 1 ]]; then
    "$PYTHON" scripts/slurm/finalize_if_complete.py --config "$P2_CONFIG" \
        --config-sha256 "$P2_CONFIG_SHA256" --expected-gpu "${P2_EXPECT_GPU:-}"
fi
