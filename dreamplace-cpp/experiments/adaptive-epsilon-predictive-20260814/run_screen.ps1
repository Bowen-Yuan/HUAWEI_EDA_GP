param([string[]]$Variants = @("temporal", "span050", "span035"))

$ErrorActionPreference = "Stop"
$env:Path = "D:\MingGW\ucrt64\bin;" + $env:Path
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$exe = Join-Path $project "dreamplace_cpp.exe"
$benchmarkDir = (Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005\adaptec1")).Path
$benchmark = Join-Path $benchmarkDir "adaptec1"
$root = Join-Path $project "output\adaptive_epsilon_predictive_20260814\screen"
$spanCaps = @{ temporal = 0.0; span050 = 0.50; span035 = 0.35 }

foreach ($variant in $Variants) {
    if (-not $spanCaps.ContainsKey($variant)) { throw "Unknown variant: $variant" }
    $output = Join-Path $root $variant
    New-Item -ItemType Directory -Force $output | Out-Null
    $args = @(
        $benchmark, "--exact-hpwl", "--optimizer", "adam", "--refine-optimizer", "adam",
        "--iterations", 800, "--bins", 512, "--active-set-radius", 125,
        "--refine-active-set-radius", 50, "--active-set-power", 4,
        "--step-fraction", 0.003, "--feasible-refinement",
        "--adaptive-active-predictive", "--adaptive-active-refinement",
        "--adaptive-active-min-scale", 0.60, "--adaptive-active-max-scale", 1.35,
        "--adaptive-active-window", 50, "--adaptive-active-interval", 25,
        "--adaptive-active-gain", 0.55, "--adaptive-active-deadband", 0.0025,
        "--adaptive-active-max-step", 0.025, "--adaptive-active-refine-max-step", 0.006,
        "--adaptive-active-span-cap", $spanCaps[$variant],
        "--dreamplace-detailed", "--dp-passes", 5, "--legal-refine-rounds", 5,
        "--cell-insertion-passes", 2, "--legal-projected-passes", 4,
        "--row-relegalization-passes", 2, "--independent-set-size", 8,
        "--hungarian-matching", "--log-every", 100, "--output", $output
    )
    $wall = [Diagnostics.Stopwatch]::StartNew()
    & $exe @args 2>&1 | Tee-Object -FilePath (Join-Path $output "run.log")
    if ($LASTEXITCODE -ne 0) { throw "$variant failed: $LASTEXITCODE" }
    $wall.Stop()
    "wall_seconds=$($wall.Elapsed.TotalSeconds)" | Set-Content (Join-Path $output "timing.txt")
}
