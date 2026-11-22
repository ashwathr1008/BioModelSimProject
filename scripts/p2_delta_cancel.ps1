# Cancel this project's Delta jobs and save verification locally. Does not start fits.
$ErrorActionPreference = 'Stop'
$projectDir = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectDir
$transferDir = Join-Path $projectDir 'results/p2/delta_transfer'
$plan = Get-Content (Join-Path $transferDir 'resume_plan.json') -Raw | ConvertFrom-Json
if ($plan.remote_project -ne '/projects/bjvz/araghuram/BioModelSimProject-20261121T222509Z_ea0e56dc') {
    throw 'Unexpected remote project path.'
}
$fingerprint = & ssh-keygen -lf 'config/delta_known_hosts'
if ($LASTEXITCODE -ne 0 -or $fingerprint -notmatch 'SHA256:mIv7OqNSuuWnu1tlY//wG4gWAn4z0mbE11KkBXDAY/g') {
    throw 'Delta host key verification failed.'
}
$source = Get-Content 'scripts/slurm/cancel_delta_campaign.sh' -Raw
$encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($source.Replace("`r", '')))
$remoteCommand = 'printf %s ' + $encoded + ' | base64 --decode | bash -l -s -- ' + $plan.remote_project
$logPath = Join-Path $transferDir ('cancel-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '.log')
$statePath = Join-Path $transferDir 'cancel_status.json'
$state = [ordered]@{status='connecting'; log=$logPath; remote_project=$plan.remote_project;
    started_at=[DateTime]::UtcNow.ToString('o')}
$state | ConvertTo-Json | Set-Content -LiteralPath $statePath -Encoding UTF8
Write-Host 'Enter your NCSA password here and complete Duo to cancel the RegVelo jobs.'
& ssh -F NUL -o StrictHostKeyChecking=yes -o UserKnownHostsFile=config/delta_known_hosts -o PreferredAuthentications=keyboard-interactive,password -o ConnectTimeout=20 -o ServerAliveInterval=30 'araghuram@login.delta.ncsa.illinois.edu' $remoteCommand | Tee-Object -FilePath $logPath
$sshCode = $LASTEXITCODE
$resultLine = @(Get-Content -LiteralPath $logPath | Where-Object { $_.StartsWith('P2_CANCEL_RESULT ') })
$verified = $false
if ($sshCode -eq 0 -and $resultLine.Count -eq 1) {
    $result = $resultLine[0].Substring(17) | ConvertFrom-Json
    $verified = ($result.status -eq 'verified_no_active_jobs' -and $result.remaining -eq 0)
    $state.result = $result
}
$state.status = if ($verified) { 'verified_no_active_jobs' } else { 'verification_failed' }
$state.ssh_exit_code = $sshCode
$state.finished_at = [DateTime]::UtcNow.ToString('o')
$state | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statePath -Encoding UTF8
if (-not $verified) { throw "Cancellation not verified; inspect $logPath" }
Write-Host 'Verified: no RegVelo campaign jobs remain in the Delta queue. Local training remains stopped.'
