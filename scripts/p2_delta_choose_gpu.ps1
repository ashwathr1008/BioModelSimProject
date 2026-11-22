# Compare queue estimates and use the earliest practical supported GPU queue.
param([switch]$CheckOnly)
$ErrorActionPreference = 'Stop'
$projectDir = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectDir
$transferDir = Join-Path $projectDir 'results/p2/delta_transfer'
$plan = Get-Content (Join-Path $transferDir 'resume_plan.json') -Raw | ConvertFrom-Json
if ($plan.remote_project -notmatch '^/projects/bjvz/araghuram/BioModelSimProject-[A-Za-z0-9_-]+$') {
    throw 'Unexpected remote project path.'
}
$fingerprint = & ssh-keygen -lf 'config/delta_known_hosts'
if ($LASTEXITCODE -ne 0 -or $fingerprint -notmatch 'SHA256:mIv7OqNSuuWnu1tlY//wG4gWAn4z0mbE11KkBXDAY/g') {
    throw 'Delta host key verification failed.'
}
$source = Get-Content 'scripts/slurm/switch_pending_campaign.sh' -Raw
$encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($source.Replace("`r", '')))
$mode = if ($CheckOnly) { '--check' } else { '--auto' }
$oldCampaign = $plan.remote_project + '/results/p2/slurm/campaign.7WX4WmQO'
$remoteCommand = 'printf %s ' + $encoded + ' | base64 --decode | bash -l -s -- ' + $mode + ' ' + $plan.remote_project + ' ' + $oldCampaign
$logPath = Join-Path $transferDir ('gpu-queue-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '.log')
$statePath = Join-Path $transferDir 'gpu_queue_status.json'
$state = [ordered]@{status='connecting'; mode=$mode; log=$logPath;
    remote_project=$plan.remote_project; started_at=[DateTime]::UtcNow.ToString('o')}
$state | ConvertTo-Json | Set-Content -LiteralPath $statePath -Encoding UTF8
Write-Host 'Enter your NCSA password here and complete Duo for the queue check.'
& ssh -F NUL -o StrictHostKeyChecking=yes -o UserKnownHostsFile=config/delta_known_hosts -o PreferredAuthentications=keyboard-interactive,password -o ConnectTimeout=20 -o ServerAliveInterval=30 'araghuram@login.delta.ncsa.illinois.edu' $remoteCommand | Tee-Object -FilePath $logPath
$sshCode = $LASTEXITCODE
$state.status = if ($sshCode -eq 0) { 'finished' } else { 'error' }
$state.ssh_exit_code = $sshCode
$state.finished_at = [DateTime]::UtcNow.ToString('o')
$state | ConvertTo-Json | Set-Content -LiteralPath $statePath -Encoding UTF8
Write-Host "Queue check exited with code $sshCode. Saved to $logPath"
if ($sshCode -ne 0) { throw "Queue check failed; inspect $logPath" }
