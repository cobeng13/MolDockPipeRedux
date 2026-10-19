param(
    [string]$EnvironmentName = 'moldockpipe-clean',
    [switch]$CheckOnly
)
$ErrorActionPreference = 'Stop'
$repositoryRoot = Split-Path -Parent $PSScriptRoot
$env:PYTHONPATH = Join-Path $repositoryRoot 'src'
$condaExecutable = $null
if ($env:CONDA_EXE -and (Test-Path -LiteralPath $env:CONDA_EXE)) {
    $condaExecutable = $env:CONDA_EXE
}
if (-not $condaExecutable) {
    $command = Get-Command conda.exe -ErrorAction SilentlyContinue
    if ($command) { $condaExecutable = $command.Source }
}
if (-not $condaExecutable) {
    foreach ($base in @($env:USERPROFILE, $env:LOCALAPPDATA, $env:ProgramData)) {
        if (-not $base) { continue }
        foreach ($distribution in @('miniconda3', 'anaconda3', 'miniforge3')) {
            $candidate = Join-Path $base "$distribution\Scripts\conda.exe"
            if (Test-Path -LiteralPath $candidate) {
                $condaExecutable = $candidate
                break
            }
        }
        if ($condaExecutable) { break }
    }
}
if (-not $condaExecutable) {
    throw 'Conda was not found. Install the course Conda environment first, or run this launcher from a Conda terminal.'
}
Push-Location -LiteralPath $repositoryRoot
try {
    if ($CheckOnly) {
        & $condaExecutable run --no-capture-output -n $EnvironmentName python -c "from pathlib import Path; import moldockpipe; from moldockpipe.receptors.ad4zn import check_ad4zn_environment; e=check_ad4zn_environment(Path.cwd()); print('Application:', moldockpipe.__file__); print('AD4Zn tools:', e.components); print('Setup:', '; '.join(e.problems) or 'Files found; runtime checks run before docking')"
    } else {
        & $condaExecutable run --no-capture-output -n $EnvironmentName python -m moldockpipe ui
    }
    $launchExitCode = $LASTEXITCODE
    if ($launchExitCode -ne 0) {
        Write-Host "Launch/check failed. Verify the '$EnvironmentName' Conda environment using README.md."
    }
} finally {
    Pop-Location
}
exit $launchExitCode
