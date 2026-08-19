param([switch]$Force)
$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "16"
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$exe = Join-Path $project "dreamplace_cpp.exe"
$benchRoot = (Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005")).Path
$outRoot = Join-Path $project "output\epsilon_continuation_20260819\repro"
$runs = @(
    @{ name="adaptec1"; bins=512; iterations=1400; active=125; refine=50 },
    @{ name="adaptec2"; bins=1024; iterations=1600; active=165; refine=65 }
)
foreach ($run in $runs) {
    $out = Join-Path $outRoot $run.name
    if ((Test-Path (Join-Path $out "summary.txt")) -and -not $Force) {
        Write-Host "[skip] $($run.name)"
        continue
    }
    New-Item -ItemType Directory -Force $out | Out-Null
    $benchmark = Join-Path (Join-Path $benchRoot $run.name) $run.name
    $args = @($benchmark, "--exact-hpwl", "--optimizer", "adam", "--refine-optimizer", "adam",
        "--iterations", $run.iterations, "--bins", $run.bins,
        "--active-set-radius", $run.active, "--refine-active-set-radius", $run.refine,
        "--active-set-power", 4, "--step-fraction", "0.003", "--feasible-refinement",
        "--epsilon-continuation-iterations", 200,
        "--epsilon-continuation-span-ratio", "0.15",
        "--epsilon-continuation-optimizer", "adam",
        "--epsilon-continuation-lr-scale", "0.75",
        "--exact-subgradient-iterations", 100,
        "--exact-subgradient-optimizer", "amsgrad",
        "--exact-subgradient-lr-scale", "0.1",
        "--exact-subgradient-filter-backtracks", 10,
        "--legal-checkpoints", "--legal-checkpoint-interval", 100,
        "--log-every", 100, "--output", $out)
    Write-Host "[run] $($run.name)"
    $sw = [Diagnostics.Stopwatch]::StartNew()
    & $exe @args 2>&1 | Tee-Object -FilePath (Join-Path $out "run.log")
    $exit = $LASTEXITCODE
    $sw.Stop()
    "wall_seconds=$($sw.Elapsed.TotalSeconds)" | Set-Content (Join-Path $out "timing.txt")
    if ($exit -ne 0) { throw "$($run.name) failed with exit code $exit" }
}
