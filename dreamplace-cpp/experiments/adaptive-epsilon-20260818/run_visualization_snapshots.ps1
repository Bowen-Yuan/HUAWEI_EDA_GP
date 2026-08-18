param(
    [string[]]$Datasets = @(
        "adaptec1", "adaptec2", "adaptec3", "adaptec4",
        "bigblue1", "bigblue2"
    ),
    [int]$SnapshotEvery = 100
)

$ErrorActionPreference = "Stop"
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$executable = Join-Path $project "dreamplace_cpp.exe"
$benchmarkRoot = (Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005")).Path
$outputRoot = Join-Path $project "output\adaptive_epsilon_20260818_visualization"
$specs = @{
    adaptec1 = @{ bins = 512; iterations = 1400; active = 125; refine = 50 }
    adaptec2 = @{ bins = 1024; iterations = 1600; active = 165; refine = 65 }
    adaptec3 = @{ bins = 1024; iterations = 1600; active = 265; refine = 105 }
    adaptec4 = @{ bins = 1024; iterations = 1600; active = 265; refine = 105 }
    bigblue1 = @{ bins = 512; iterations = 1600; active = 125; refine = 50 }
    bigblue2 = @{ bins = 1024; iterations = 1600; active = 212; refine = 85 }
}

if ($SnapshotEvery -le 0) { throw "SnapshotEvery must be positive" }
New-Item -ItemType Directory -Force $outputRoot | Out-Null
foreach ($dataset in $Datasets) {
    if (-not $specs.ContainsKey($dataset)) { throw "Unknown dataset: $dataset" }
    $spec = $specs[$dataset]
    $benchmark = Join-Path (Join-Path $benchmarkRoot $dataset) $dataset
    $output = Join-Path $outputRoot $dataset
    $snapshots = Join-Path $output "snapshots"
    New-Item -ItemType Directory -Force $output,$snapshots | Out-Null
    $arguments = @(
        $benchmark, "--exact-hpwl", "--optimizer", "adam", "--refine-optimizer", "adam",
        "--iterations", $spec.iterations, "--bins", $spec.bins,
        "--active-set-radius", $spec.active, "--refine-active-set-radius", $spec.refine,
        "--active-set-power", 4, "--step-fraction", 0.003, "--feasible-refinement",
        "--legal-checkpoints", "--legal-checkpoint-interval", 50,
        "--adaptive-active-smart", "--adaptive-active-min-scale", 0.70,
        "--adaptive-active-max-scale", 1.10, "--adaptive-active-window", 50,
        "--adaptive-active-interval", 25, "--adaptive-active-gain", 0.50,
        "--adaptive-active-deadband", 0.0025, "--adaptive-active-max-step", 0.03,
        "--snapshot-every", $SnapshotEvery, "--snapshot-dir", $snapshots,
        "--dreamplace-detailed", "--dp-passes", 5, "--legal-refine-rounds", 5,
        "--cell-insertion-passes", 2, "--legal-projected-passes", 4,
        "--row-relegalization-passes", 2, "--independent-set-size", 8,
        "--hungarian-matching", "--log-every", 100, "--output", $output
    )
    $sw = [Diagnostics.Stopwatch]::StartNew()
    & $executable @arguments *> (Join-Path $output "run.log")
    $exitCode = $LASTEXITCODE
    $sw.Stop()
    if ($exitCode -ne 0) { throw "$dataset failed with exit code $exitCode" }
    "visualization_runtime_seconds=$($sw.Elapsed.TotalSeconds)" | Set-Content (Join-Path $output "runtime.txt")
    Write-Host ("[complete] {0} snapshot_every={1} runtime_seconds={2:N3}" -f $dataset, $SnapshotEvery, $sw.Elapsed.TotalSeconds)
}
