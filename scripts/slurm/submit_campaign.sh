#!/usr/bin/env bash
# Submit all 18 fits in one session; Slurm releases the rest after both pilots.
set -euo pipefail
: "${P2_ACCOUNT:?Set the verified GPU allocation account}"
: "${P2_PARTITION:?Set the verified GPU partition}"
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
PROJECT_DIR="$(cd -- "$PROJECT_DIR" && pwd -P)"
P2_ENV="${P2_ENV:-$PROJECT_DIR/.venv-delta}"
source_config="${P2_SOURCE_CONFIG:-$PROJECT_DIR/config/p2_both_pca.yaml}"
[[ -f "$P2_ENV/.p2_environment_ready.json" && -f "$source_config" ]] || {
    echo 'Prepared Linux environment and configuration are required' >&2; exit 2;
}
mkdir -p -- "$PROJECT_DIR/results/p2/slurm"
campaign_dir="$(mktemp -d "$PROJECT_DIR/results/p2/slurm/campaign.XXXXXXXX")"
printf 'Campaign receipts: %s\n' "$campaign_dir"
trap 'printf "Submission stopped; inspect receipts in %s\n" "$campaign_dir" >&2' ERR
P2_SOURCE_CONFIG="$campaign_dir/config.yaml"
cp -- "$source_config" "$P2_SOURCE_CONFIG"
P2_EXPECT_GPU="${P2_EXPECT_GPU:-H200}"
P2_REUSE_LOCAL_FULL0="${P2_REUSE_LOCAL_FULL0:-0}"
export PROJECT_DIR P2_ENV P2_SOURCE_CONFIG P2_EXPECT_GPU P2_REUSE_LOCAL_FULL0
pilot_array='0-1%2'
case "$P2_REUSE_LOCAL_FULL0" in
    0) ;;
    1)
        : "${P2_IMPORT_COMPLETED:?Supply the authorized X/full_r0 transfer receipt}"
        : "${P2_IMPORT_COMPLETED_SHA256:?Supply its independently verified SHA-256}"
        cp -- "$P2_IMPORT_COMPLETED" "$campaign_dir/imported_full0.json"
        P2_IMPORT_COMPLETED="$campaign_dir/imported_full0.json"
        export P2_IMPORT_COMPLETED P2_IMPORT_COMPLETED_SHA256
        "$P2_ENV/bin/python" "$PROJECT_DIR/scripts/slurm/validate_outputs.py" \
            --config "$P2_SOURCE_CONFIG" --scope local-full0 --expected-gpu "$P2_EXPECT_GPU" \
            | tee "$campaign_dir/import_validation.json"
        pilot_array='1%1'
        echo 'Reusing the validated laptop X/full_r0; submitting 17 H200 fits.'
        ;;
    *) echo 'P2_REUSE_LOCAL_FULL0 must be 0 or 1' >&2; exit 2 ;;
esac

P2_ARRAY="$pilot_array" P2_REQUIRE_PILOTS=0 P2_COMPARE_WHEN_COMPLETE=0 \
    P2_JOB_ID_FILE="$campaign_dir/pilots.job_id" \
    bash "$PROJECT_DIR/scripts/slurm/submit_array.sh"
if [[ "${P2_PRINT_ONLY:-0}" = 1 ]]; then
    pilot_id=1  # Numeric placeholder only for the printed dependency plan.
else
    pilot_id="$(< "$campaign_dir/pilots.job_id")"
    pilot_id="${pilot_id%%;*}"
fi
P2_ARRAY="2-17%${P2_CAMPAIGN_CONCURRENCY:-4}" P2_REQUIRE_PILOTS=1 \
    P2_COMPARE_WHEN_COMPLETE=1 P2_DEPENDENCY="afterok:$pilot_id" \
    P2_JOB_ID_FILE="$campaign_dir/remainder.job_id" \
    bash "$PROJECT_DIR/scripts/slurm/submit_array.sh"
if [[ "${P2_PRINT_ONLY:-0}" = 1 ]]; then
    printf 'Plan only; dependency ID 1 is a placeholder. Receipts: %s\n' "$campaign_dir"
    exit 0
fi
remainder_id="$(< "$campaign_dir/remainder.job_id")"
remainder_id="${remainder_id%%;*}"
printf 'pilots\t%s\nremainder\t%s\n' "$pilot_id" "$remainder_id" > "$campaign_dir/jobs.tsv"
cat -- "$campaign_dir/jobs.tsv"
printf 'Comparison completion marker: %s/results/p2/slurm/comparison_complete.json\n' "$PROJECT_DIR"
printf 'Submitted campaign; receipts: %s\n' "$campaign_dir"
