$ErrorActionPreference = 'Stop'
$taskProject = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$taskDownloads = Join-Path $taskProject '.tools\linux-runtime'
$taskDistro = 'PBMC-Ubuntu'
$taskArchive = Join-Path $taskDownloads 'ubuntu-24.04.5-wsl-amd64.wsl'

if (Test-Path 'HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Component Based Servicing\RebootPending') {
    throw 'Restart Windows before completing the Linux installation. This script does not restart the computer.'
}

$taskInstalled = ((& wsl.exe --list --quiet) -join "`n").Replace("`0", '')
if ($taskInstalled -notmatch '(?m)^PBMC-Ubuntu\s*$') {
    if (!(Test-Path -LiteralPath $taskArchive)) { throw 'The prepared Ubuntu installer is missing.' }
    & wsl.exe --install --from-file $taskArchive --name $taskDistro --location (Join-Path $env:LOCALAPPDATA 'PBMC-Ubuntu') --no-launch
    if ($LASTEXITCODE -ne 0) { throw 'Ubuntu registration failed.' }
}

$taskMambaBinary = Join-Path $taskDownloads 'micromamba'
if (!(Test-Path -LiteralPath $taskMambaBinary)) {
    & (Join-Path $taskProject '.venv\Scripts\python.exe') -c 'import pathlib,sys,tarfile; p=pathlib.Path(sys.argv[1]); t=tarfile.open(p / "micromamba-linux-64.tar.bz2"); (p / "micromamba").write_bytes(t.extractfile("bin/micromamba").read()); t.close()' $taskDownloads
    if ($LASTEXITCODE -ne 0) { throw 'Micromamba extraction failed.' }
}

$taskLinuxProject = '/mnt/' + $taskProject.Substring(0,1).ToLowerInvariant() + $taskProject.Substring(2).Replace('\','/')
& wsl.exe --distribution $taskDistro --user root --exec bash "$taskLinuxProject/environment/salmon/setup.sh" $taskLinuxProject
if ($LASTEXITCODE -ne 0) { throw 'Salmon environment setup failed.' }
