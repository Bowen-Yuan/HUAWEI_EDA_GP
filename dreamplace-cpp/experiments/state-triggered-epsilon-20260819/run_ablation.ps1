param(
    [string[]]$Datasets = @("adaptec1", "adaptec2")
)

$ErrorActionPreference = "Stop"
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$executable = Join-Path $project "dreamplace_cpp.exe"
$benchmarkRoot = (Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005")).Path
$outputRoot = Join-Path $project "output\state_triggered_epsilon_20260819\ablation"

$specs = @{
    adaptec1 = @{ bins = 512; iterations = 1400; fallback = 1100; active = 125; refine = 50 }
    adaptec2 = @{ bins = 1024; iterations = 1600; fallback = 1300; active = 165; refine = 65 }
}
$variants = @(
    @{ name = "state15"; ratio = 0.15; continuation = 300; exact = 0 },
    @{ name = "state10"; ratio = 0.10; continuation = 300; exact = 0 },
    @{ name = "state15_exact50"; ratio = 0.15; continuation = 250; exact = 50 }
)

function Read-Summary([string]$path) {
    $values = @{}
    foreach ($line in Get-Content $path) {
        if ($line -match '^([^=]+)=(.*)$') { $values[$matches[1]] = $matches[2] }
    }
    return $values
}

New-Item -ItemType Directory -Force $outputRoot | Out-Null
$records = @()
foreach ($dataset in $Datasets) {
    if (-not $specs.ContainsKey($dataset)) { throw "Unknown dataset: $dataset" }
    $spec = $specs[$dataset]
    foreach ($variant in $variants) {
        $benchmark = Join-Path (Join-Path $benchmarkRoot $dataset) $dataset
        $output = Join-Path $outputRoot (Join-Path $variant.name $dataset)
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
            "--epsilon-continuation-iterations", $variant.continuation,
            "--epsilon-continuation-span-ratio", $variant.ratio,
            "--epsilon-continuation-lr-scale", 0.75,
            "--epsilon-continuation-optimizer", "adam",
            "--epsilon-continuation-state-trigger",
            "--epsilon-trigger-min-refine", 250,
            "--epsilon-trigger-window", 50,
            "--epsilon-trigger-overflow-upper", 0.085,
            "--epsilon-trigger-overflow-range", 0.004,
            "--epsilon-trigger-min-hpwl-drop", 0.0003,
            "--exact-subgradient-iterations", $variant.exact,
            "--exact-subgradient-lr-scale", 0.10,
            "--exact-subgradient-optimizer", "amsgrad",
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
        if ($exitCode -ne 0) { throw "$dataset/$($variant.name) failed: $exitCode" }
        $summary = Read-Summary (Join-Path $output "summary.txt")
        $runtime = $stopwatch.Elapsed.TotalSeconds
        "runtime_seconds=$runtime" | Set-Content (Join-Path $output "runtime.txt")
        $records += [pscustomobject]@{
            dataset = $dataset
            variant = $variant.name
            total_iterations = $spec.iterations
            actual_transition = [int]$summary["epsilon_continuation_actual_start"]
            selected_iteration = [int]$summary["selected_legal_iteration"]
            runtime_seconds = $runtime
            gp_seconds = [double]$summary["gp_seconds"]
            gp_hpwl = [double]$summary["gp_hpwl"]
            gp_overflow = [double]$summary["gp_overflow"]
            final_legal_hpwl = [double]$summary["detailed_hpwl"]
            internal_legal = $summary["internal_legal"]
            boundary_errors = [int]$summary["legal_boundary_errors"]
            alignment_errors = [int]$summary["legal_alignment_errors"]
            overlap_errors = [int]$summary["legal_overlap_errors"]
        }
        $records | Export-Csv (Join-Path $outputRoot "results_summary.csv") -NoTypeInformation
        Write-Host ("[complete] {0}/{1} legal={2:N0} runtime={3:N2}s" -f `
            $dataset, $variant.name, [double]$summary["detailed_hpwl"], $runtime)
    }
}
