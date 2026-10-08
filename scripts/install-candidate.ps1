param(
    [Parameter(Mandatory=$true)][string]$CandidateDll,
    [Parameter(Mandatory=$true)][string]$KspDirectory,
    [string]$BackupDirectory = (Join-Path (Split-Path -Parent $PSScriptRoot) 'build\install-backups'),
    [string]$SourceCommit
)

$ErrorActionPreference = 'Stop'
function Assert-KspClosed {
    $running = @(Get-Process -Name KSP_x64 -ErrorAction SilentlyContinue)
    if ($running.Count -gt 0) { throw "Close KSP before installation. Running PID(s): $($running.Id -join ', ')." }
}

Assert-KspClosed
$source = (Get-Item -LiteralPath $CandidateDll).FullName
$target = Join-Path (Get-Item -LiteralPath $KspDirectory).FullName 'GameData\PatchedConicFixes\Plugins\PatchedConicFixes.dll'
if (!(Test-Path -LiteralPath $target)) { throw "Existing PatchedConicFixes DLL was not found: $target" }
if ([string]::Equals($source, $target, [StringComparison]::OrdinalIgnoreCase)) { throw 'Candidate must be a separate build artifact.' }
$assembly = [Reflection.AssemblyName]::GetAssemblyName($source)
if ($assembly.Name -ne 'PatchedConicFixes') { throw 'Candidate is not a PatchedConicFixes assembly.' }
$originalHash = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
$candidateHash = (Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
$originalVersion = [Reflection.AssemblyName]::GetAssemblyName($target).Version.ToString()
New-Item -ItemType Directory -Force -Path $BackupDirectory | Out-Null
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$backup = Join-Path $BackupDirectory "PatchedConicFixes-before-$stamp.dll"
Copy-Item -LiteralPath $target -Destination $backup
Assert-KspClosed
Copy-Item -LiteralPath $source -Destination $target -Force
$installedHash = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
if ($installedHash -ne $candidateHash) {
    Copy-Item -LiteralPath $backup -Destination $target -Force
    throw "Installed hash did not match; restored the previous DLL from $backup."
}
$manifest = [ordered]@{
    installedAtUtc = [DateTime]::UtcNow.ToString('o')
    sourceCommit = $SourceCommit
    candidate = $source
    installed = $target
    candidateVersion = $assembly.Version.ToString()
    candidateSha256 = $candidateHash
    backup = $backup
    previousVersion = $originalVersion
    previousSha256 = $originalHash
}
$manifestPath = Join-Path $BackupDirectory "installation-$stamp.json"
$manifest | ConvertTo-Json | Set-Content -LiteralPath $manifestPath -Encoding UTF8
$manifest | ConvertTo-Json
