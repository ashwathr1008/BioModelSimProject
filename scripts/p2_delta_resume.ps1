# Reuse the uploaded inputs and completed local fits in the verified allocation.
param([switch]$SkipUpload)
$ErrorActionPreference = 'Stop'
$projectDir = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectDir
$transferDir = Join-Path $projectDir 'results/p2/delta_transfer'
$plan = Get-Content (Join-Path $transferDir 'resume_plan.json') -Raw | ConvertFrom-Json
if ($plan.account -ne 'bjvz-delta-gpu' -or $plan.partition -ne 'gpuH200x8' -or $plan.remaining_fits -ne 17) {
    throw 'Unexpected allocation or run count.'
}
if ($plan.remote_project -notmatch '^/projects/bjvz/araghuram/BioModelSimProject-[A-Za-z0-9_-]+$' -or
    $plan.base_remote_archive -notmatch '^/u/araghuram/regvelo-p2-[A-Za-z0-9_-]+\.tar$' -or
    $plan.overlay_remote_archive -notmatch '^/u/araghuram/regvelo-overlay-[A-Za-z0-9_-]+\.tar$') {
    throw 'Unexpected remote path.'
}
foreach ($field in @('base_sha256', 'overlay_sha256', 'import_receipt_sha256')) {
    if ($plan.$field -notmatch '^[a-f0-9]{64}$') { throw "Invalid checksum: $field" }
}
$fingerprint = & ssh-keygen -lf 'config/delta_known_hosts'
if ($LASTEXITCODE -ne 0 -or $fingerprint -notmatch 'SHA256:mIv7OqNSuuWnu1tlY//wG4gWAn4z0mbE11KkBXDAY/g') {
    throw 'Delta host key verification failed.'
}
$archive = (Resolve-Path -LiteralPath $plan.overlay_local_archive).Path
if ((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $plan.overlay_sha256) {
    throw 'The overlay archive has changed.'
}
$statePath = Join-Path $transferDir 'resume_status.json'
$livePath = Join-Path $transferDir 'remote_live_status.json'
$logPath = Join-Path $transferDir ('resume-session-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '.log')
$state = [ordered]@{status='ready'; remote_project=$plan.remote_project; remaining_fits=17;
    account=$plan.account; partition=$plan.partition; log=$logPath; started_at=[DateTime]::UtcNow.ToString('o')}
function Save-ResumeState {
    $state.updated_at = [DateTime]::UtcNow.ToString('o')
    $state | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $statePath -Encoding UTF8
}
$sshOptions = @('-F', 'NUL', '-o', 'StrictHostKeyChecking=yes',
    '-o', 'UserKnownHostsFile=config/delta_known_hosts',
    '-o', 'PreferredAuthentications=keyboard-interactive,password',
    '-o', 'ConnectTimeout=20', '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=6')
try {
    if (-not $SkipUpload) {
        $state.status = 'uploading_completed_fits_and_code'
        Save-ResumeState
        Write-Host 'Upload completed fits and updated code: enter your NCSA password here, then complete Duo.'
        & scp @sshOptions $archive ('araghuram@login.delta.ncsa.illinois.edu:' + $plan.overlay_remote_archive)
        if ($LASTEXITCODE -ne 0) { throw "Overlay upload failed with code $LASTEXITCODE." }
        $state.upload_finished_at = [DateTime]::UtcNow.ToString('o')
    }
    $remoteScript = @'
set -euo pipefail
exec 2>&1
umask 077
[[ $(id -un) = araghuram ]] || exit 2
base=__BASE__
overlay=__OVERLAY__
project=__PROJECT__
printf '%s  %s\n' __BASE_HASH__ "$base" | sha256sum --check -
printf '%s  %s\n' __OVERLAY_HASH__ "$overlay" | sha256sum --check -
echo 'Global and storage-pool quota diagnostics:'
lfs quota -h -p 23080 /work/hdd/bjvz || true
lfs quota -h -p 23080 --pool ddn_hdd /work/hdd/bjvz || true
lfs quota -h -p 23080 --pool ddn_ssd /work/nvme/bjvz || true
quota_report=$(quota)
printf '%s\n' "$quota_report"
df -h /projects/bjvz /u/araghuram
if [[ ! -e "$project" ]]; then
    quota_free_kib=$(printf '%s\n' "$quota_report" | awk -F'|' '
        function kib(value, suffix, number) {
            gsub(/[[:space:]]/, "", value)
            suffix = substr(value, length(value), 1)
            number = substr(value, 1, length(value)-1) + 0
            if (suffix == "K") return number
            if (suffix == "M") return number * 1024
            if (suffix == "G") return number * 1024 * 1024
            if (suffix == "T") return number * 1024 * 1024 * 1024
            return -1
        }
        $2 ~ /^[[:space:]]*\/projects\/bjvz[[:space:]]*$/ {
            used = kib($3); limit = kib($4)
            if (used >= 0 && limit > 0) printf "%.0f\n", limit - used - 128 * 1024
        }')
    [[ "$quota_free_kib" =~ ^[0-9]+$ ]] || { echo 'Cannot determine project quota headroom'; exit 2; }
    (( quota_free_kib >= 6 * 1024 * 1024 )) || { echo 'Need at least 6 GiB project quota headroom before extraction'; exit 2; }
    available_kib=$(df -Pk /projects/bjvz | awk 'NR == 2 {print $4}')
    [[ "$available_kib" =~ ^[0-9]+$ ]] || { echo 'Cannot determine project free space'; exit 2; }
    (( available_kib >= 6 * 1024 * 1024 )) || { echo 'Need at least 6 GiB project headroom before extraction'; exit 2; }
    mkdir -p -- "$project"
    tar --no-same-owner --no-same-permissions --no-overwrite-dir -xf "$base" -C "$project"
    (cd -- "$project" && sha256sum --check SHA256SUMS)
    printf '%s\n' __BASE_HASH__ > "$project/.base_verified_sha256"
fi
[[ -f "$project/.base_verified_sha256" ]] || { echo 'Incomplete prior extraction; preserve it and use a new destination'; exit 2; }
[[ $(cat "$project/.base_verified_sha256") = __BASE_HASH__ ]] || { echo 'Base receipt differs'; exit 2; }
if [[ ! -f "$project/.overlay_verified_sha256" ]]; then
    tar --no-same-owner --no-same-permissions --no-overwrite-dir -xf "$overlay" -C "$project"
    (cd -- "$project" && sha256sum --check OVERLAY_SHA256SUMS)
    printf '%s\n' __OVERLAY_HASH__ > "$project/.overlay_verified_sha256"
fi
[[ $(cat "$project/.overlay_verified_sha256") = __OVERLAY_HASH__ ]] || { echo 'Overlay receipt differs'; exit 2; }
cd -- "$project"
sha256sum --check OVERLAY_SHA256SUMS
export P2_REUSE_LOCAL_FULL0=1
export P2_IMPORT_COMPLETED="$project/results/p2/slurm/imported_local_full0.json"
export P2_IMPORT_COMPLETED_SHA256=__IMPORT_HASH__
export P2_MIN_FREE_GIB=2
bash scripts/slurm/bootstrap_delta.sh "$project"
'@
    $remoteScript = $remoteScript.Replace('__BASE__', $plan.base_remote_archive).Replace('__OVERLAY__', $plan.overlay_remote_archive).Replace('__PROJECT__', $plan.remote_project).Replace('__BASE_HASH__', $plan.base_sha256).Replace('__OVERLAY_HASH__', $plan.overlay_sha256).Replace('__IMPORT_HASH__', $plan.import_receipt_sha256).Replace("`r", '')
    $encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($remoteScript))
    $remoteCommand = 'printf %s ' + $encoded + ' | base64 --decode | bash -l'
    $state.status = 'connecting_for_storage_check_setup_and_submission'
    Save-ResumeState
    Write-Host 'Storage check, setup and submission: enter your NCSA password again and complete Duo.'
    Write-Host 'Keep this window open during setup; it will stream Slurm progress afterward.'
    & ssh @sshOptions 'araghuram@login.delta.ncsa.illinois.edu' $remoteCommand |
        Tee-Object -FilePath $logPath | ForEach-Object {
            $line = $_.ToString()
            if ($line.StartsWith('P2_STATUS ')) {
                $line.Substring(10) | ConvertFrom-Json | ConvertTo-Json -Depth 12 |
                    Set-Content -LiteralPath $livePath -Encoding UTF8
                $state.status = 'monitoring_slurm_jobs'
                Save-ResumeState
            }
            Write-Host $line
        }
    $sshCode = $LASTEXITCODE
    if ($sshCode -ne 0) { throw "Remote setup/session exited with code $sshCode; see $logPath" }
    $state.status = 'remote_session_finished'
    Save-ResumeState
} catch {
    $state.status = 'error'
    $state.error = $_.Exception.Message
    Save-ResumeState
    throw
}
