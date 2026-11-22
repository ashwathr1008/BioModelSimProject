# Interactive, read-only Delta access check. Password and Duo stay in this terminal.
$ErrorActionPreference = 'Stop'
$projectDir = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectDir
$knownHosts = Join-Path $projectDir 'config/delta_known_hosts'
if (-not (Test-Path -LiteralPath $knownHosts)) {
    throw 'Missing verified Delta host key file.'
}
$fingerprint = & ssh-keygen -lf $knownHosts
if ($fingerprint -notmatch 'SHA256:mIv7OqNSuuWnu1tlY//wG4gWAn4z0mbE11KkBXDAY/g') {
    throw 'Delta host key does not match the published NCSA ED25519 fingerprint.'
}
$commands = @'
hostname
id -un
accounts
quota
sinfo -p gpuH200x8,gpuA100x4,gpuA100x8,gpuA40x4
scontrol show partition gpuH200x8
squeue -u araghuram
sacctmgr -nP show assoc where user=araghuram format=Account,Partition,QOS,DefaultQOS
module spider python 2>&1
module spider cuda 2>&1
'@
$commands = $commands.Replace("`r", '')
if ($commands.Contains("'")) { throw 'Unexpected quoting in remote inventory command' }
$remoteCommand = 'bash -lc ' + "'" + $commands + "'"
$output = Join-Path $projectDir 'results/p2/delta_access.txt'
Write-Host 'Enter your NCSA password and complete Duo below. Do not paste credentials into chat.'
& ssh -F NUL -o StrictHostKeyChecking=yes -o UserKnownHostsFile=config/delta_known_hosts -o PreferredAuthentications=keyboard-interactive,password -o ConnectTimeout=20 araghuram@login.delta.ncsa.illinois.edu $remoteCommand | Tee-Object -FilePath $output
$sshCode = $LASTEXITCODE
Write-Host "SSH exited with code $sshCode. Inventory saved to $output"
Read-Host 'Press Enter to close this access-check window'
