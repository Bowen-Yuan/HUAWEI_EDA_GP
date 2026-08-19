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
$outputRoot = Join-Path $project "output\state_triggered_epsilon_20260819\frozen15"
$perl = "D:\latex\101\texlive\2023\tlpkg\tlperl\bin\perl.exe"
$legalChecker = Join-Path $benchmarkRoot "legal2.pl\legal2.pl"

$specs = @{
    adaptec1 = @{ bins = 512; iterations = 1400; fallback = 1100; active = 125; refine = 50 }
    adaptec2 = @{ bins = 1024; iterations = 1600; fallback = 1300; active = 165; refine = 65 }
    adaptec3 = @{ bins = 1024; iterations = 1600; fallback = 1300; active = 265; refine = 105 }
    adaptec4 = @{ bins = 1024; iterations = 1600; fallback = 1300; active = 265; refine = 105 }
    bigblue1 = @{ bins = 512; iterations = 1600; fallback = 1300; active = 125; refine = 50 }
    bigblue2 = @{ bins = 1024; iterations = 1600; fallback = 1300; active = 212; refine = 85 }
}

function Read-Summary([string]$path) {
    $values = @{}
    foreach ($line in Get-Content $path) {
        if ($line -match '^([^=]+)=(.*)$') { $values[$matches[1]] = $matches[2] }
    }
    return $values
}

function Read-LegalErrors([string]$path) {
    $errors = @(0, 0, 0, 0)
    foreach ($line in Get-Content $path) {
        if ($line -match '^\s*([0-3])\s+(\d+)\s*$') {
            $errors[[int]$matches[1]] = [int]$matches[2]
        }
    }
    return $errors
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
        "--epsilon-continuation-start-iteration", $spec.fallback,
        "--epsilon-continuation-iterations", 300,
        "--epsilon-continuation-span-ratio", 0.15,
        "--epsilon-continuation-lr-scale", 0.75,
        "--epsilon-continuation-optimizer", "adam",
        "--epsilon-continuation-state-trigger",
        "--epsilon-trigger-min-refine", 250,
        "--epsilon-trigger-window", 50,
        "--epsilon-trigger-overflow-upper", 0.085,
        "--epsilon-trigger-overflow-range", 0.004,
        "--epsilon-trigger-min-hpwl-drop", 0.0003,
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
    $stopwatch = [Diagnostics.Stopwatch]::StartNew()
    & $executable @arguments *> $log
    $exitCode = $LASTEXITCODE
    $stopwatch.Stop()
    if ($exitCode -ne 0) { throw "$dataset failed: $exitCode" }
    $runtime = $stopwatch.Elapsed.TotalSeconds
    "runtime_seconds=$runtime" | Set-Content (Join-Path $output "runtime.txt")

    $legalLog = Join-Path $output "legal2.log"
    Push-Location $output
    try {
        & $perl $legalChecker "$benchmark.nodes" "$benchmark.pl" `
            (Join-Path $output "final.pl") "$benchmark.scl" 20 *> $legalLog
        if ($LASTEXITCODE -ne 0) { throw "$dataset legal2.pl failed: $LASTEXITCODE" }
    } finally {
        Pop-Location
    }
    $summary = Read-Summary (Join-Path $output "summary.txt")
    $legalErrors = @(Read-LegalErrors $legalLog)
    $records += [pscustomobject]@{
        dataset = $dataset
        iterations = $spec.iterations
        actual_transition = [int]$summary["epsilon_continuation_actual_start"]
        selected_iteration = [int]$summary["selected_legal_iteration"]
        runtime_seconds = $runtime
        gp_seconds = [double]$summary["gp_seconds"]
        gp_hpwl = [double]$summary["gp_hpwl"]
        gp_overflow = [double]$summary["gp_overflow"]
        final_legal_hpwl = [double]$summary["detailed_hpwl"]
        internal_legal = $summary["internal_legal"]
        legal2_type0 = $legalErrors[0]
        legal2_type1 = $legalErrors[1]
        legal2_type2 = $legalErrors[2]
        legal2_type3 = $legalErrors[3]
        output = $output
    }
    $records | Export-Csv (Join-Path $outputRoot "results_partial.csv") -NoTypeInformation
    Write-Host ("[complete] {0} legal={1:N0} runtime={2:N2}s checker={3}/{4}/{5}/{6}" -f `
        $dataset, [double]$summary["detailed_hpwl"], $runtime,
        $legalErrors[0], $legalErrors[1], $legalErrors[2], $legalErrors[3])
}
