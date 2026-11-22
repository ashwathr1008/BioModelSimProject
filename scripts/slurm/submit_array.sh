#!/usr/bin/env bash
# Run from the remote project root after checking allocation and partition access.
set -euo pipefail

: "${P2_ACCOUNT:?Set the verified Delta GPU allocation account, not your username}"
: "${P2_PARTITION:?Set a GPU partition verified for your allocation}"
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
PROJECT_DIR="$(cd -- "$PROJECT_DIR" && pwd -P)"
P2_ENV="${P2_ENV:-$PROJECT_DIR/.venv-delta}"
source_config="${P2_SOURCE_CONFIG:-$PROJECT_DIR/config/p2_both_pca.yaml}"
[[ -f "$P2_ENV/.p2_environment_ready.json" ]] || {
    echo 'The Linux environment has not passed setup checks' >&2; exit 2;
}
[[ -f "$source_config" && -f "$PROJECT_DIR/scripts/p2_fit.py" ]] || {
    echo 'Project source/configuration missing' >&2; exit 2;
}
mkdir -p -- "$PROJECT_DIR/results/p2/slurm"
submission_dir="$(mktemp -d "$PROJECT_DIR/results/p2/slurm/submission.XXXXXXXX")"
P2_CONFIG="$submission_dir/config.yaml"
cp -- "$source_config" "$P2_CONFIG"
P2_CONFIG_SHA256="$(sha256sum -- "$P2_CONFIG" | cut -d ' ' -f 1)"
export PROJECT_DIR P2_ENV P2_CONFIG P2_CONFIG_SHA256

# First submit P2_ARRAY=0-1%2 for the two full-data pilots. After those pass,
# P2_ARRAY=2-17%N submits the rest, with N chosen from live allocation limits.
command=(sbatch --parsable --job-name=regvelo-p2
    --account="$P2_ACCOUNT" --partition="$P2_PARTITION"
    --nodes=1 --ntasks=1 --gpus-per-node=1
    --cpus-per-task="${P2_CPUS:-8}" --mem="${P2_MEMORY:-32G}"
    --time="${P2_TIME:-04:00:00}" --array="${P2_ARRAY:-0-1%2}"
    --chdir="$PROJECT_DIR" --export=ALL
    --output="$PROJECT_DIR/results/p2/slurm/%A_%a.out"
    --error="$PROJECT_DIR/results/p2/slurm/%A_%a.err")
if [[ -n "${P2_DEPENDENCY:-}" ]]; then
    command+=(--dependency="$P2_DEPENDENCY" --kill-on-invalid-dep=yes)
fi
command+=("$PROJECT_DIR/scripts/slurm/fit_array.sh")
printf '%q ' "${command[@]}" > "$submission_dir/command.txt"
printf '\n' >> "$submission_dir/command.txt"
cat -- "$submission_dir/command.txt"
if [[ "${P2_PRINT_ONLY:-0}" = 1 ]]; then
    echo "Prepared submission only: $submission_dir"
    exit 0
fi
job_id="$("${command[@]}")"
[[ "$job_id" =~ ^[0-9]+(\;[[:alnum:]_-]+)?$ ]] || {
    echo "Unexpected sbatch response: $job_id" >&2; exit 2;
}
printf '%s\n' "$job_id" > "$submission_dir/job_id.txt"
if [[ -n "${P2_JOB_ID_FILE:-}" ]]; then
    (set -o noclobber; printf '%s\n' "$job_id" > "$P2_JOB_ID_FILE")
fi
printf 'Submitted %s; snapshot and receipt: %s\n' "$job_id" "$submission_dir"
