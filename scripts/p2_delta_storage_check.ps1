# Read-only diagnosis of the failed Delta extraction.
$ErrorActionPreference = 'Stop'
$projectDir = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectDir
$plan = Get-Content 'results/p2/delta_transfer/deployment.json' -Raw | ConvertFrom-Json
if ($plan.remote_project -notmatch '^/work/hdd/bjvz/araghuram/BioModelSimProject-[A-Za-z0-9_-]+$') {
    throw 'Unexpected remote project path.'
}
$fingerprint = & ssh-keygen -lf 'config/delta_known_hosts'
if ($LASTEXITCODE -ne 0 -or $fingerprint -notmatch 'SHA256:mIv7OqNSuuWnu1tlY//wG4gWAn4z0mbE11KkBXDAY/g') {
    throw 'Delta host key does not match the verified NCSA fingerprint.'
}
$commands = @'
set +e
exec 2>&1
failed=__FAILED__
printf '\n--- Identity and quota ---\n'
hostname
id
quota
printf '\n--- Failed extraction ownership ---\n'
stat -c '%A %a %U %G %u:%g %s %n' "$failed" "$failed/data" "$failed/data/prepared_p2/full_r0.h5ad"
find "$failed" -maxdepth 4 -name '*.h5ad' -printf '%u %g %s %p\n'
for base in /work/hdd/bjvz /work/nvme/bjvz; do
    printf '\n--- Storage area: %s ---\n' "$base"
    stat -c '%A %a %U %G %u:%g %s %n' "$base" "$base/araghuram"
    df -h "$base"
    df -i "$base"
    if command -v lfs >/dev/null 2>&1; then
        lfs quota -h -u araghuram "$base"
        group_id=$(stat -c '%g' "$base")
        lfs quota -h -g "$group_id" "$base"
        lfs project -d "$base" "$base/araghuram"
        lfs getstripe -d "$base"
        lfs pool_list "$base"
    fi
done
if command -v lfs >/dev/null 2>&1; then
    printf '\n--- Failed path Lustre metadata ---\n'
    lfs project -d "$failed" "$failed/data/prepared_p2/full_r0.h5ad"
    lfs getstripe -d "$failed"
    lfs getstripe "$failed/data/prepared_p2/full_r0.h5ad"
    file_group=$(stat -c '%g' "$failed/data/prepared_p2/full_r0.h5ad")
    lfs quota -h -g "$file_group" "$failed"
fi
printf '\n--- Current Slurm jobs ---\n'
squeue -u araghuram
printf '\nStorage diagnostics finished.\n'
'@
$commands = $commands.Replace('__FAILED__', $plan.remote_project).Replace("`r", '')
$encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($commands))
$remoteCommand = 'printf %s ' + $encoded + ' | base64 --decode | bash -l'
$outputPath = Join-Path $projectDir 'results/p2/delta_transfer/storage_check.log'
Write-Host 'Enter your NCSA password here and complete Duo for this storage check.'
& ssh -F NUL -o StrictHostKeyChecking=yes -o UserKnownHostsFile=config/delta_known_hosts -o PreferredAuthentications=keyboard-interactive,password -o ConnectTimeout=20 'araghuram@login.delta.ncsa.illinois.edu' $remoteCommand | Tee-Object -FilePath $outputPath
$sshCode = $LASTEXITCODE
Write-Host "Storage check exited with code $sshCode. Saved to $outputPath"
if ($sshCode -ne 0) { throw "SSH storage check failed with code $sshCode." }
