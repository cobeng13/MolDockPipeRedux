$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

try {
    $condaHook = & conda shell.powershell hook
    if ($LASTEXITCODE -ne 0) {
        throw 'Conda PowerShell initialization failed.'
    }
    $condaHook | Out-String | Invoke-Expression

    $global:LASTEXITCODE = 0
    conda activate moldockpipe-clean
    if ($LASTEXITCODE -ne 0 -or $env:CONDA_DEFAULT_ENV -ne 'moldockpipe-clean') {
        throw 'Conda environment activation failed.'
    }
}
catch {
    Write-Warning "Conda could not activate the moldockpipe-clean environment. $($_.Exception.Message)"
    Write-Host 'Ensure Conda is installed and available, and that the environment exists.'
    Write-Host 'Try running this launcher from an Anaconda PowerShell Prompt.'
    Read-Host 'Press Enter to close'
    exit 1
}

try {
    $global:LASTEXITCODE = 0
    & moldockpipe
    $launchExit = $LASTEXITCODE
    if ($launchExit -ne 0) {
        Write-Warning "MolDockPipe exited with error code $launchExit."
        Read-Host 'Press Enter to close'
    }
    exit $launchExit
}
catch {
    Write-Warning "MolDockPipe could not run. $($_.Exception.Message)"
    Read-Host 'Press Enter to close'
    exit 1
}
