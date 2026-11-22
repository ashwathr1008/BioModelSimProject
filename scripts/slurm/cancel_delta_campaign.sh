#!/usr/bin/env bash
# Cancel only this project's recorded Delta jobs; preserve data and results.
set -euo pipefail
exec 2>&1
project=${1:?Supply the deployed project directory}
[[ $(id -un) == araghuram ]] || { echo 'Unexpected SSH user'; exit 2; }
[[ "$project" == /projects/bjvz/araghuram/BioModelSimProject-20261121T222509Z_ea0e56dc ]] || {
    echo 'Unexpected project directory'; exit 2;
}
declare -A expected=([24333181]=1 [24333183]=1 [24333584]=1 [24333585]=1)
for campaign in campaign.7WX4WmQO campaign.maqqyj3K; do
    for role in pilots remainder; do
        receipt="$project/results/p2/slurm/$campaign/$role.job_id"
        job=$(< "$receipt")
        job=${job%%;*}
        [[ "$job" =~ ^[0-9]+$ && -n ${expected[$job]:-} ]] || {
            echo "Unexpected job receipt: $receipt"; exit 2;
        }
    done
done
queue=$(squeue --noheader --user=araghuram --format='%i|%j|%a|%T')
printf 'Queue before cancellation:\n%s\n' "$queue"
declare -A active=()
while IFS='|' read -r task name account state; do
    job=${task%%_*}
    [[ -n "$job" && -n ${expected[$job]:-} ]] || continue
    [[ "$name" == regvelo-p2 && "$account" == bjvz-delta-gpu ]] || {
        echo "Job identity mismatch: $task $name $account"; exit 2;
    }
    active[$job]=1
done <<< "$queue"
if (( ${#active[@]} )); then
    # Cancel the dependent remainder before its pilot, avoiding an orphaned array.
    for job in 24333585 24333584 24333183 24333181; do
        if [[ -n ${active[$job]:-} ]]; then
            scancel --user=araghuram --account=bjvz-delta-gpu --name=regvelo-p2 "$job" || true
        fi
    done
fi
for attempt in {1..10}; do
    queue=$(squeue --noheader --user=araghuram --format='%i|%j|%a|%T')
    remaining=0
    while IFS='|' read -r task name account state; do
        job=${task%%_*}
        if [[ -n "$job" && -n ${expected[$job]:-} ]]; then
            printf 'Still in queue: %s|%s|%s|%s\n' "$task" "$name" "$account" "$state"
            remaining=$((remaining + 1))
        fi
    done <<< "$queue"
    if (( remaining == 0 )); then
        sacct -X --jobs=24333181,24333183,24333584,24333585 --noheader --parsable2 \
            --format=JobIDRaw,JobName,State,Elapsed || true
        printf '%s\n' 'P2_CANCEL_RESULT {"status":"verified_no_active_jobs","job_ids":["24333181","24333183","24333584","24333585"],"remaining":0}'
        exit 0
    fi
    sleep 2
done
echo 'Cancellation requested, but some jobs remain in the queue; verification incomplete.'
exit 3
