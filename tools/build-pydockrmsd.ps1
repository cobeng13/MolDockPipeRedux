param([string] $Python = "python")

$ErrorActionPreference = "Stop"
$revision = "6642effa194abfdb5f216ee80c66d65e23c6eef5"
$root = Split-Path -Parent $PSScriptRoot
$source = Join-Path $root "build/pydockrmsd-$revision"
$wheels = Join-Path $root "build/wheels"

& $Python -c "import platform, sys; sys.exit(not (sys.platform == 'win32' and sys.version_info[:2] == (3, 13) and platform.machine().upper() == 'AMD64'))"
if ($LASTEXITCODE -ne 0) { throw "Use Windows x64 CPython 3.13." }

New-Item -ItemType Directory -Path (Split-Path -Parent $source), $wheels -Force | Out-Null
if (-not (Test-Path -LiteralPath $source)) {
    git clone --quiet https://github.com/neudinger/pyDockRMSD.git $source
    if ($LASTEXITCODE -ne 0) { throw "Could not clone PyDockRMSD." }
}
git -C $source checkout --quiet --detach $revision
if ($LASTEXITCODE -ne 0) { throw "Could not check out the pinned PyDockRMSD revision." }
if (git -C $source status --porcelain --untracked-files=no) {
    throw "The pinned PyDockRMSD source has local changes."
}

$previousVersion = $env:version
try {
    $env:version = "1.0.0"
    & $Python -m pip wheel --no-deps --wheel-dir $wheels $source
    if ($LASTEXITCODE -ne 0) { throw "PyDockRMSD wheel build failed." }
} finally {
    $env:version = $previousVersion
}
