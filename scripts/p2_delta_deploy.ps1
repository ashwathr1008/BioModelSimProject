# Authenticate only in this terminal. No password or Duo code is stored by this script.
param([switch]$SkipUpload)
$ErrorActionPreference = 'Stop'
$projectDir = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectDir
$planPath = Join-Path $projectDir 'results/p2/delta_transfer/deployment.json'
$plan = Get-Content -LiteralPath $planPath -Raw | ConvertFrom-Json
if ($plan.account -ne 'bjvz-delta-gpu' -or $plan.partition -ne 'gpuH200x8') {
    throw 'Deployment differs from the verified Delta allocation.'
}
if ($plan.remote_project -notmatch '^/work/hdd/bjvz/araghuram/BioModelSimProject-[A-Za-z0-9_-]+$' -or
    $plan.remote_archive -notmatch '^/u/araghuram/regvelo-p2-[A-Za-z0-9_-]+\.tar$' -or
    $plan.archive_sha256 -notmatch '^[a-f0-9]{64}$') {
    throw 'Invalid deployment path or checksum.'
}
$fingerprint = & ssh-keygen -lf 'config/delta_known_hosts'
if ($LASTEXITCODE -ne 0 -or $fingerprint -notmatch 'SHA256:mIv7OqNSuuWnu1tlY//wG4gWAn4z0mbE11KkBXDAY/g') {
    throw 'Delta host key does not match the verified NCSA fingerprint.'
}
$archive = (Resolve-Path -LiteralPath $plan.archive_path).Path
if ((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $plan.archive_sha256) {
    throw 'The local upload archive has changed.'
}
$transferDir = Split-Path -Parent $planPath
$statePath = Join-Path $transferDir 'deployment_status.json'
$livePath = Join-Path $transferDir 'remote_live_status.json'
$logPath = Join-Path $transferDir ('session-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '.log')
$state = [ordered]@{ status = 'ready'; remote_project = $plan.remote_project; account = $plan.account;
    partition = $plan.partition; log = $logPath; started_at = [DateTime]::UtcNow.ToString('o') }
function Save-DeploymentState {
    $state.updated_at = [DateTime]::UtcNow.ToString('o')
    $state | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $statePath -Encoding UTF8
}
$sshOptions = @('-F', 'NUL', '-o', 'StrictHostKeyChecking=yes',
    '-o', 'UserKnownHostsFile=config/delta_known_hosts',
    '-o', 'PreferredAuthentications=keyboard-interactive,password',
    '-o', 'ConnectTimeout=20', '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=6')
try {
    if (-not $SkipUpload) {
        $state.status = 'uploading'
        Save-DeploymentState
        Write-Host 'Upload: enter your NCSA password here and complete Duo. Never paste credentials into chat.'
        $destination = 'araghuram@login.delta.ncsa.illinois.edu:' + $plan.remote_archive
        & scp @sshOptions $archive $destination
        if ($LASTEXITCODE -ne 0) { throw "Upload failed with exit code $LASTEXITCODE." }
        $state.upload_finished_at = [DateTime]::UtcNow.ToString('o')
    }
    $remoteScript = @'
set -euo pipefail
exec 2>&1
umask 077
upload=__ARCHIVE__
project=__PROJECT__
printf '%s  %s\n' __SHA256__ "$upload" | sha256sum --check -
if [[ ! -e "$project" ]]; then
    mkdir -p -- "$project"
    tar --no-same-owner -xf "$upload" -C "$project"
fi
cd -- "$project"
sha256sum --check SHA256SUMS
bash scripts/slurm/bootstrap_delta.sh "$project"
'@
    $remoteScript = $remoteScript.Replace('__ARCHIVE__', $plan.remote_archive).Replace('__PROJECT__', $plan.remote_project).Replace('__SHA256__', $plan.archive_sha256).Replace("`r", '')
    $encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($remoteScript))
    $remoteCommand = 'printf %s ' + $encoded + ' | base64 --decode | bash -l'
    $state.status = 'connecting_for_setup_and_submission'
    Save-DeploymentState
    Write-Host 'Setup and submission: enter your password here again and complete Duo.'
    Write-Host 'Keep this window open during installation. After submission it will stream job status.'
    & ssh @sshOptions 'araghuram@login.delta.ncsa.illinois.edu' $remoteCommand |
        Tee-Object -FilePath $logPath | ForEach-Object {
            $line = $_.ToString()
            if ($line.StartsWith('P2_STATUS ')) {
                $live = $line.Substring(10) | ConvertFrom-Json
                $live | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $livePath -Encoding UTF8
                $state.status = 'monitoring_slurm_jobs'
                Save-DeploymentState
            }
            Write-Host $line
        }
    $sshCode = $LASTEXITCODE
    if ($sshCode -ne 0) { throw "SSH setup/session exited with code $sshCode. See $logPath" }
    $state.status = 'remote_session_finished'
    Save-DeploymentState
    Write-Host "Session finished. Local status and log: $transferDir"
} catch {
    $state.status = 'error'
    $state.error = $_.Exception.Message
    Save-DeploymentState
    throw
}
