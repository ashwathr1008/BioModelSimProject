#!/usr/bin/env bash
# --check is read-only; --auto compares controller ETAs; --switch forces A100 selection.
set -euo pipefail
mode="${1:---check}"
[[ "$mode" = --check || "$mode" = --switch || "$mode" = --auto ]] || { echo 'Use --check, --auto, or --switch' >&2; exit 2; }
PROJECT_DIR="$(cd -- "${2:?Supply the deployed project directory}" && pwd -P)"
old_campaign="$(cd -- "${3:?Supply the original campaign receipt directory}" && pwd -P)"
case "$old_campaign" in
    "$PROJECT_DIR"/results/p2/slurm/campaign.*) ;;
    *) echo 'Campaign receipts must belong to this project' >&2; exit 2 ;;
esac
old_ids=()
for name in pilots remainder; do
    job="$(< "$old_campaign/$name.job_id")"; job="${job%%;*}"
    [[ "$job" =~ ^[0-9]+$ ]] || { echo 'Invalid original job receipt' >&2; exit 2; }
    old_ids+=("$job")
done
[[ "${old_ids[0]}" != "${old_ids[1]}" ]] || { echo 'Original job IDs must differ' >&2; exit 2; }
account="${P2_ACCOUNT:-bjvz-delta-gpu}"
test_partition() {
    sbatch --test-only --account="$account" --partition="$1" --nodes=1 --ntasks=1 \
        --gpus-per-node=1 --cpus-per-task=8 --mem=32G --time=04:00:00 --wrap=true
}
printf 'Login-node UTC: '; date -u +%FT%TZ
printf '\nOriginal controller job timestamps (SubmitTime, EligibleTime, StartTime):\n'
old_job="$(scontrol show job -o "${old_ids[0]}")"
printf '%s\n' "$old_job"
squeue --start -r -j "${old_ids[0]},${old_ids[1]}" -o '%i|%P|%T|%S|%R'
declare -A estimates
for partition in gpuH200x8 gpuA100x4 gpuA100x8 gpuA40x4; do
    printf '\nStart estimate for %s:\n' "$partition"
    if estimate="$(test_partition "$partition" 2>&1)"; then
        estimates[$partition]="$estimate"
    else
        estimates[$partition]=''
    fi
    printf '%s\n' "$estimate"
done
echo 'Compare controller StartTime with controller SubmitTime/EligibleTime; login-node and controller clocks may differ.'
[[ "$mode" != --check ]] || exit 0
selected_partition='gpuA100x4,gpuA100x8'
expected_gpu=A100
decision_json='{"mode":"explicit","partition":"gpuA100x4,gpuA100x8"}'
if [[ "$mode" = --auto ]]; then
    for job in "${old_ids[@]}"; do
        states="$(squeue -r -h -j "$job" -o '%T')"
        if [[ -z "$states" ]] || [[ "${states//PENDING/}" =~ [^[:space:]] ]]; then
            echo "AUTO: original job $job is no longer wholly pending; keeping the current campaign unchanged."
            exit 0
        fi
    done
    P2_ENV="$(< "$PROJECT_DIR/results/p2/slurm/environment_path.txt")"
    selection="$("$P2_ENV/bin/python" - "$old_job" \
        "${estimates[gpuH200x8]}" "${estimates[gpuA100x4]}" \
        "${estimates[gpuA100x8]}" "${estimates[gpuA40x4]}" <<'PY'
import datetime as dt
import json
import re
import sys

stamp = r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})"
def parse(text, pattern):
    match = re.search(pattern + r"\s*" + stamp, text, re.I)
    try:
        return dt.datetime.fromisoformat(match.group(1)) if match else None
    except ValueError:
        return None

old = sys.argv[1]
submitted = parse(old, r"SubmitTime=")
original_start = parse(old, r"StartTime=")
partitions = ("gpuH200x8", "gpuA100x4", "gpuA100x8", "gpuA40x4")
estimates = {p: parse(text, r"(?:to start at|Estimated start(?: time)?[:=]?|StartTime=)")
             for p, text in zip(partitions, sys.argv[2:])}
baseline = original_start or estimates["gpuH200x8"]
baseline_source = "existing_job" if original_start is not None else "test_only_estimate"
estimates["gpuH200x8"] = baseline
invalid = [p for p, value in estimates.items() if value is not None and submitted and value < submitted]
known = {p: value for p, value in estimates.items() if value is not None}
selected, gpu, reason = "UNKNOWN", "", "No trustworthy current H200/controller ETA; keeping jobs unchanged"
if invalid:
    reason = "Controller ETAs precede the original controller SubmitTime; clock uncertainty, keeping jobs unchanged"
elif baseline is not None:
    earliest = min(known.values())
    slack = dt.timedelta(minutes=30)
    if baseline <= earliest + slack:
        selected, gpu, reason = "KEEP", "H200", "H200 starts within 30 minutes of the earliest parseable estimate"
    else:
        a100 = [(value, p) for p, value in known.items() if p.startswith("gpuA100")]
        if a100 and min(a100)[0] <= earliest + slack:
            _, selected = min(a100)
            gpu, reason = "A100", "A100 starts within 30 minutes of the earliest estimate; H200 does not"
        elif estimates["gpuA40x4"] is not None:
            selected, gpu, reason = "gpuA40x4", "A40", "A40 is more than 30 minutes earlier than the faster GPU tiers with known ETAs"
chosen_start = baseline if selected == "KEEP" else estimates.get(selected)
record = {"mode": "auto", "selection": selected, "expected_gpu": gpu, "reason": reason,
          "h200_baseline_source": baseline_source,
          "controller_submit_time": submitted.isoformat() if submitted else None,
          "controller_start_estimates": {p: value.isoformat() if value else None for p, value in estimates.items()},
          "unknown_partitions": [p for p, value in estimates.items() if value is None],
          "minutes_earlier_than_original": (baseline - chosen_start).total_seconds() / 60 if baseline and chosen_start else None,
          "illustrative_finish_at_three_4h_waves": (chosen_start + dt.timedelta(hours=12)).isoformat() if chosen_start else None,
          "projection_caveat": "Controller clock only; assumes a 4h pilot then two 4h waves with immediate resources. This is not a completion prediction."}
print(selected + "\t" + (gpu or "UNKNOWN") + "\t" + json.dumps(record, separators=(",", ":")))
PY
)"
    IFS=$'\t' read -r selected_partition expected_gpu decision_json <<< "$selection"
    printf 'P2_SWITCH_DECISION %s\n' "$decision_json"
    [[ "$selected_partition" != KEEP && "$selected_partition" != UNKNOWN ]] || exit 0
fi
test_partition "$selected_partition"
cd -- "$PROJECT_DIR"
exec 9>results/p2/slurm/.switch.lock
flock --exclusive 9
history="$(mktemp -d "$PROJECT_DIR/results/p2/slurm/switch.XXXXXXXX")"
printf 'Switch receipts: %s\n' "$history"
printf '%s\n' "$decision_json" > "$history/decision.json"
printf '%s\n' "${old_ids[@]}" > "$history/original_job_ids.txt"
export P2_SWITCH_NEW_IDS="$history/replacement_job_ids.txt"
: > "$P2_SWITCH_NEW_IDS"
phase=checking
pending_only() {
    local rows state
    rows="$(squeue -r -h -j "$1" -o '%T')"
    [[ -n "$rows" ]] || { echo "Job $1 is absent from the queue" >&2; return 1; }
    while IFS= read -r state; do
        [[ "$state" = PENDING ]] || { echo "Job $1 has a $state task; refusing to switch" >&2; return 1; }
    done <<< "$rows"
}
rollback() {
    local status=$? job unsafe=0
    trap - EXIT
    if (( status != 0 )) && [[ "$phase" = submitting ]]; then
        while IFS= read -r job; do
            job="${job%%;*}"
            [[ "$job" =~ ^[0-9]+$ ]] || continue
            if ! pending_only "$job"; then unsafe=1; fi
            scancel --state=PENDING "$job" || true
        done < "$P2_SWITCH_NEW_IDS"
        if (( unsafe == 0 )); then
            for job in "${old_ids[@]}"; do scontrol release "$job" || true; done
            echo 'Replacement submission failed; original pending jobs were released.' >&2
        else
            echo 'A replacement is no longer wholly pending; originals remain held for inspection.' >&2
        fi
    elif (( status != 0 )) && [[ "$phase" = committing ]]; then
        echo "Replacement IDs are recorded in $history; inspect their held/released states before retrying." >&2
    fi
    exit "$status"
}
trap rollback EXIT
trap 'exit 130' HUP INT TERM
for job in "${old_ids[@]}"; do pending_only "$job"; done
phase=submitting
for job in "${old_ids[@]}"; do scontrol hold "$job"; done
for job in "${old_ids[@]}"; do pending_only "$job"; done
export PROJECT_DIR
P2_ENV="$(< "$PROJECT_DIR/results/p2/slurm/environment_path.txt")"
export P2_ENV P2_ACCOUNT="$account" P2_PARTITION="$selected_partition" P2_EXPECT_GPU="$expected_gpu"
export P2_CPUS=8 P2_MEMORY=32G P2_TIME=04:00:00 P2_CAMPAIGN_CONCURRENCY=8
export P2_SOURCE_CONFIG="$old_campaign/config.yaml" P2_REUSE_LOCAL_FULL0=1
export P2_IMPORT_COMPLETED="$old_campaign/imported_full0.json"
export P2_IMPORT_COMPLETED_SHA256=63923d1800c81cc0c5d90b77b73754526554b493c13becb474d4e8c367196ff2
(
    # Hold replacements until both IDs are recorded, even if one could start immediately.
    sbatch() {
        local reply
        reply="$(command sbatch --hold "$@")" || return $?
        [[ "$reply" =~ ^[0-9]+(\;[[:alnum:]_-]+)?$ ]] || { echo "Unexpected sbatch reply: $reply" >&2; return 2; }
        printf '%s\n' "$reply" >> "$P2_SWITCH_NEW_IDS"
        printf '%s\n' "$reply"
    }
    export -f sbatch
    bash "$PROJECT_DIR/scripts/slurm/submit_campaign.sh"
) | tee "$history/submission.log"
mapfile -t replacement_ids < "$P2_SWITCH_NEW_IDS"
[[ "${#replacement_ids[@]}" = 2 ]] || { echo 'Expected exactly two replacement job receipts' >&2; exit 2; }
for index in 0 1; do replacement_ids[$index]="${replacement_ids[$index]%%;*}"; done
[[ "${replacement_ids[0]}" != "${replacement_ids[1]}" ]] || { echo 'Replacement job IDs must differ' >&2; exit 2; }
for job in "${old_ids[@]}" "${replacement_ids[@]}"; do pending_only "$job"; done
phase=committing
# State filtering prevents cancellation of a task that raced into RUNNING.
scancel --state=PENDING "${old_ids[@]}"
for job in "${replacement_ids[@]}"; do scontrol release "$job"; done
phase=complete
printf 'pilots\t%s\nremainder\t%s\n' "${replacement_ids[0]}" "${replacement_ids[1]}" > "$history/replacement_jobs.tsv"
printf 'Switched pending campaign to %s: %s, %s. Receipts: %s\n' "$expected_gpu" "${replacement_ids[0]}" "${replacement_ids[1]}" "$history"
