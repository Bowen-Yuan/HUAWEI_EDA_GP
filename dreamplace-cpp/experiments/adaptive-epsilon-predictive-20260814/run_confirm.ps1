param([string[]]$Datasets = @("adaptec1", "adaptec2"))

$ErrorActionPreference = "Stop"
$env:Path = "D:\MingGW\ucrt64\bin;" + $env:Path
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$exe = Join-Path $project "dreamplace_cpp.exe"
$root = Join-Path $project "output\adaptive_epsilon_predictive_20260814\confirm"
$specs = @{
    adaptec1 = @{ bins = 512; iterations = 1400; active = 125; refine = 50 }
    adaptec2 = @{ bins = 1024; iterations = 1600; active = 165; refine = 65 }
}
foreach ($dataset in $Datasets) {
    if (-not $specs.ContainsKey($dataset)) { throw "Unknown dataset: $dataset" }
    $s = $specs[$dataset]
    $benchmark = Join-Path ((Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005\$dataset")).Path) $dataset
    $output = Join-Path $root $dataset
    New-Item -ItemType Directory -Force $output | Out-Null
    $args = @(
        $benchmark, "--exact-hpwl", "--optimizer", "adam", "--refine-optimizer", "adam",
        "--iterations", $s.iterations, "--bins", $s.bins,
        "--active-set-radius", $s.active, "--refine-active-set-radius", $s.refine,
        "--active-set-power", 4, "--step-fraction", 0.003, "--feasible-refinement",
        "--adaptive-active-predictive", "--adaptive-active-refinement",
        "--adaptive-active-min-scale", 0.75, "--adaptive-active-max-scale", 1.20,
        "--adaptive-active-window", 50, "--adaptive-active-interval", 25,
        "--adaptive-active-gain", 0.35, "--adaptive-active-deadband", 0.0025,
        "--adaptive-active-max-step", 0.015, "--adaptive-active-refine-max-step", 0.004,
        "--dreamplace-detailed", "--dp-passes", 5, "--legal-refine-rounds", 5,
        "--cell-insertion-passes", 2, "--legal-projected-passes", 4,
        "--row-relegalization-passes", 2, "--independent-set-size", 8,
        "--hungarian-matching", "--log-every", 100, "--output", $output
    )
    $sw = [Diagnostics.Stopwatch]::StartNew()
    & $exe @args 2>&1 | Tee-Object -FilePath (Join-Path $output "run.log")
    if ($LASTEXITCODE -ne 0) { throw "$dataset failed: $LASTEXITCODE" }
    $sw.Stop()
    "wall_seconds=$($sw.Elapsed.TotalSeconds)" | Set-Content (Join-Path $output "timing.txt")
}
