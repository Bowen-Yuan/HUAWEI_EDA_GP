param(
    [string[]]$Datasets = @(
        "adaptec1", "adaptec2", "adaptec3", "adaptec4",
        "bigblue1", "bigblue2"
    )
)

$ErrorActionPreference = "Stop"
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$executable = Join-Path $project "dreamplace_cpp.exe"
$benchmarkRoot = (Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005")).Path
$outputRoot = Join-Path $project "output\adaptive_epsilon_20260818_performance"

$specs = @{
    adaptec1 = @{ bins = 512; iterations = 1400; active = 125; refine = 50 }
    adaptec2 = @{ bins = 1024; iterations = 1600; active = 165; refine = 65 }
    adaptec3 = @{ bins = 1024; iterations = 1600; active = 265; refine = 105 }
    adaptec4 = @{ bins = 1024; iterations = 1600; active = 265; refine = 105 }
    bigblue1 = @{ bins = 512; iterations = 1600; active = 125; refine = 50 }
    bigblue2 = @{ bins = 1024; iterations = 1600; active = 212; refine = 85 }
}

New-Item -ItemType Directory -Force $outputRoot | Out-Null
$records = @()
foreach ($dataset in $Datasets) {
    if (-not $specs.ContainsKey($dataset)) { throw "Unknown dataset: $dataset" }
    $spec = $specs[$dataset]
    $benchmark = Join-Path (Join-Path $benchmarkRoot $dataset) $dataset
    $output = Join-Path $outputRoot $dataset
    New-Item -ItemType Directory -Force $output | Out-Null
    $arguments = @(
        $benchmark,
        "--exact-hpwl",
        "--optimizer", "adam",
        "--refine-optimizer", "adam",
        "--iterations", $spec.iterations,
        "--bins", $spec.bins,
        "--active-set-radius", $spec.active,
        "--refine-active-set-radius", $spec.refine,
        "--active-set-power", 4,
        "--step-fraction", 0.003,
        "--feasible-refinement",
        "--legal-checkpoints",
        "--legal-checkpoint-interval", 50,
        "--adaptive-active-smart",
        "--adaptive-active-min-scale", 0.70,
        "--adaptive-active-max-scale", 1.10,
        "--adaptive-active-window", 50,
        "--adaptive-active-interval", 25,
        "--adaptive-active-gain", 0.50,
        "--adaptive-active-deadband", 0.0025,
        "--adaptive-active-max-step", 0.03,
        "--dreamplace-detailed",
        "--dp-passes", 5,
        "--legal-refine-rounds", 5,
        "--cell-insertion-passes", 2,
        "--legal-projected-passes", 4,
        "--row-relegalization-passes", 2,
        "--independent-set-size", 8,
        "--hungarian-matching",
        "--log-every", 50,
        "--output", $output
    )
    $log = Join-Path $output "run.log"
    $sw = [Diagnostics.Stopwatch]::StartNew()
    & $executable @arguments *> $log
    $exitCode = $LASTEXITCODE
    $sw.Stop()
    if ($exitCode -ne 0) { throw "$dataset failed with exit code $exitCode" }
    $summary = Join-Path $output "summary.txt"
    $iteration = $spec.iterations
    if (Test-Path $summary) {
        $summaryMap = @{}
        foreach ($line in Get-Content $summary) {
            if ($line -match '^([^=]+)=(.*)$') { $summaryMap[$matches[1]] = $matches[2] }
        }
        if ($summaryMap.ContainsKey("iterations")) { $iteration = [int]$summaryMap["iterations"] }
    }
    $runtime = $sw.Elapsed.TotalSeconds
    "runtime_seconds=$runtime" | Set-Content (Join-Path $output "runtime.txt")
    $records += [pscustomobject]@{
        dataset = $dataset
        iterations = $iteration
        runtime_seconds = $runtime
        output = $output
    }
    Write-Host ("[complete] {0} iterations={1} runtime_seconds={2:N3}" -f $dataset, $iteration, $runtime)
}
$records | Export-Csv (Join-Path $outputRoot "runtime_summary.csv") -NoTypeInformation
