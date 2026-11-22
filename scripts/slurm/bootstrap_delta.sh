#!/usr/bin/env bash
# Called after the authenticated upload has been checksum-verified and extracted.
set -eo pipefail
PROJECT_DIR="$(cd -- "${1:?Supply the extracted project directory}" && pwd -P)"
export PROJECT_DIR
[[ "$(id -un)" = araghuram ]] || { echo 'Unexpected Delta user' >&2; exit 2; }
allowed_location=0
case "$PROJECT_DIR" in
    /work/hdd/bjvz/araghuram/BioModelSimProject-*) allowed_location=1 ;;
esac
# Resolve the verified project alias rather than guessing its /taiga mount path.
if [[ -d /projects/bjvz/araghuram ]]; then
    projects_root="$(cd -- /projects/bjvz/araghuram && pwd -P)"
    case "$PROJECT_DIR" in
        "$projects_root"/*) allowed_location=1 ;;
    esac
fi
[[ "$allowed_location" = 1 ]] || { echo 'Unexpected project location' >&2; exit 2; }
cd -- "$PROJECT_DIR"
mkdir -p results/p2/slurm
exec > >(tee -a results/p2/slurm/deployment.log) 2>&1
echo "Delta deployment started: $(date -u +%FT%TZ)"
environment_record="$PROJECT_DIR/results/p2/slurm/environment_path.txt"
if [[ -s "$environment_record" ]]; then
    P2_ENV="$(< "$environment_record")"
else
    P2_ENV="${P2_ENV:-/u/araghuram/.local/share/regvelo/p2-py312-cu124}"
fi
[[ "$P2_ENV" = /* ]] || { echo 'The environment path must be absolute' >&2; exit 2; }
export P2_ENV
shopt -s nullglob
existing_pilots=(results/p2/slurm/campaign.*/pilots.job_id)
if (( ${#existing_pilots[@]} )); then
    echo 'An existing campaign receipt was found; monitoring without submitting duplicates.'
    [[ -x "$P2_ENV/bin/python" ]] || { echo "Recorded environment is unavailable: $P2_ENV" >&2; exit 2; }
    "$P2_ENV/bin/python" scripts/slurm/watch_campaign.py
    exit 0
fi
hostname
accounts
quota
scontrol show partition gpuH200x8
sacctmgr -P show qos where name=bjvz-delta-gpu \
    format=Name,Flags,MaxJobsPU,MaxSubmitJobsPU,MaxTRESPU,MaxTRESPA,GrpJobs,GrpSubmitJobs,GrpTRES,GrpTRESMins,GrpTRESRunMins,MaxWall || true
sacctmgr -P show assoc where account=bjvz-delta-gpu \
    format=Account,User,Partition,QOS,DefaultQOS,GrpJobs,GrpSubmitJobs,GrpTRES,MaxJobs,MaxSubmitJobs,MaxWall || true

if ! type module >/dev/null 2>&1; then
    source /etc/profile
fi
module spider miniforge3-python
module load miniforge3-python
module list
command -v conda
set -u
umask 077
export CONDA_PKGS_DIRS=/u/araghuram/.cache/regvelo/conda
export PIP_CACHE_DIR=/u/araghuram/.cache/regvelo/pip
mkdir -p -- "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR"
python_prefix=/u/araghuram/.local/share/regvelo/python-3.12.6
mkdir -p -- "$(dirname -- "$python_prefix")"
exec 8>"${python_prefix}.setup.lock"
flock --exclusive 8
if [[ ! -x "$python_prefix/bin/python" ]]; then
    conda create --yes --prefix "$python_prefix" --override-channels -c conda-forge --no-default-packages python=3.12.6 pip
fi
"$python_prefix/bin/python" -c 'import sys; assert sys.version_info[:3] == (3, 12, 6), sys.version'
flock --unlock 8
bash scripts/slurm/setup_environment.sh "$PROJECT_DIR" "$python_prefix/bin/python" "$P2_ENV"
printf '%s\n' "$P2_ENV" > "$environment_record"
export P2_ACCOUNT=bjvz-delta-gpu P2_PARTITION=gpuH200x8
export P2_CPUS=8 P2_MEMORY=32G P2_TIME=04:00:00
sbatch --test-only --account="$P2_ACCOUNT" --partition="$P2_PARTITION" \
    --nodes=1 --ntasks=1 --gpus-per-node=1 --cpus-per-task="$P2_CPUS" \
    --mem="$P2_MEMORY" --time="$P2_TIME" --wrap=true
export P2_CAMPAIGN_CONCURRENCY="${P2_CAMPAIGN_CONCURRENCY:-8}"
bash scripts/slurm/submit_campaign.sh
echo "Delta submission finished: $(date -u +%FT%TZ)"
echo 'Slurm jobs continue even if this SSH window is closed.'
"$P2_ENV/bin/python" scripts/slurm/watch_campaign.py
