param(
    [string[]]$Datasets = @("adaptec1", "adaptec2"),
    [string[]]$Cases = @("rel10-p4-f1", "rel15-p4-f1")
)

$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "16"
$project = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$executable = Join-Path $project "dreamplace_cpp.exe"
$benchmarkRoot = (Resolve-Path (Join-Path $project "..\alg-electronic\ispd2005")).Path
$outputRoot = Join-Path $project "output\relative_epsilon_20260819"

$datasetsConfig = @{
    adaptec1 = @{ bins = 512; iterations = 1400; active = 125; refine = 50 }
    adaptec2 = @{ bins = 1024; iterations = 1600; active = 165; refine = 65 }
}
$casesConfig = @{
    "rel10-p4-f1" = @{ ratio = 0.10; power = 4; floor = 1.0; small = 10.0; density = 8.0e-5 }
    "rel15-p4-f1" = @{ ratio = 0.15; power = 4; floor = 1.0; small = 10.0; density = 8.0e-5 }
    "rel15-p2-f1" = @{ ratio = 0.15; power = 2; floor = 1.0; small = 10.0; density = 8.0e-5 }
    "rel15-p1-f1" = @{ ratio = 0.15; power = 1; floor = 1.0; small = 10.0; density = 8.0e-5 }
    "rel15-p2-lambda-match" = @{ ratio = 0.15; power = 2; floor = 1.0; small = 10.0; density = 2.25e-5 }
    "rel15-p2-f5-small50" = @{ ratio = 0.15; power = 2; floor = 5.0; small = 50.0; density = 8.0e-5 }
    "rel15-p2-f10-small100" = @{ ratio = 0.15; power = 2; floor = 10.0; small = 100.0; density = 8.0e-5 }
    "rel15-p2-f10-to-f1" = @{ ratio = 0.15; power = 2; floor = 10.0; small = 100.0; refineFloor = 1.0; refineSmall = 10.0; density = 8.0e-5 }
    "rel15-p2-f15-to-strict" = @{ ratio = 0.15; power = 2; floor = 15.0; small = 100.0; refineFloor = 0.0; refineSmall = 0.0; density = 8.0e-5 }
    "rel15-p2-f20-to-strict" = @{ ratio = 0.15; power = 2; floor = 20.0; small = 133.3333333333; refineFloor = 0.0; refineSmall = 0.0; density = 8.0e-5 }
}

New-Item -ItemType Directory -Force $outputRoot | Out-Null
foreach ($caseName in $Cases) {
    if (-not $casesConfig.ContainsKey($caseName)) { throw "Unknown case: $caseName" }
    $case = $casesConfig[$caseName]
    foreach ($dataset in $Datasets) {
        if (-not $datasetsConfig.ContainsKey($dataset)) { throw "Unknown dataset: $dataset" }
        $spec = $datasetsConfig[$dataset]
        $benchmark = Join-Path (Join-Path $benchmarkRoot $dataset) $dataset
        $output = Join-Path (Join-Path $outputRoot $caseName) $dataset
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
            "--active-set-power", $case.power,
            "--active-set-span-ratio", $case.ratio,
            "--active-set-min-radius", $case.floor,
            "--active-set-small-span-threshold", $case.small,
            "--density-weight-scale", $case.density,
            "--step-fraction", 0.003,
            "--feasible-refinement",
            "--legal-checkpoints",
            "--legal-checkpoint-interval", 50,
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
        if ($case.ContainsKey("refineFloor")) {
            $arguments += @(
                "--refine-active-set-min-radius", $case.refineFloor,
                "--refine-active-set-small-span-threshold", $case.refineSmall
            )
        }
        $log = Join-Path $output "run.log"
        $timer = [Diagnostics.Stopwatch]::StartNew()
        & $executable @arguments *> $log
        $exitCode = $LASTEXITCODE
        $timer.Stop()
        if ($exitCode -ne 0) { throw "$caseName/$dataset failed: exit $exitCode" }
        "runtime_seconds=$($timer.Elapsed.TotalSeconds)" |
            Set-Content (Join-Path $output "runtime.txt")
        Write-Host ("[complete] {0}/{1} runtime={2:N3}s" -f
            $caseName, $dataset, $timer.Elapsed.TotalSeconds)
    }
}
