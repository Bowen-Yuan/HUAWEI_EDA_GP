param(
    [int]$Threads = 16,
    [int]$Iterations = 1000,
    [string]$Optimizer = "adam",
    [string]$Tag = "baseline"
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
$Workspace = Split-Path -Parent $Root
$Executable = Join-Path $Root "build\epsilon_active.exe"
$Benchmarks = @{
    adaptec1 = @{ Base = "dataset\adaptec1\adaptec1"; Bins = 512; Epsilon = 125 }
    adaptec2 = @{ Base = "alg\ispd2005\adaptec2\adaptec2"; Bins = 1024; Epsilon = 165 }
    adaptec3 = @{ Base = "dataset\adaptec3\adaptec3"; Bins = 1024; Epsilon = 265 }
    adaptec4 = @{ Base = "alg\ispd2005\adaptec4\adaptec4"; Bins = 1024; Epsilon = 265 }
}

if ($Threads -lt 1 -or $Threads -gt 40) {
    throw "Threads must be between 1 and 40"
}
foreach ($Name in "adaptec1", "adaptec2", "adaptec3", "adaptec4") {
    $Spec = $Benchmarks[$Name]
    $Benchmark = Join-Path $Workspace $Spec.Base
    $Output = Join-Path $Root "output\$Tag\$Name"
    & $Executable --benchmark $Benchmark --output $Output `
        --bins $Spec.Bins --iterations $Iterations --threads $Threads `
        --optimizer $Optimizer --lambda-policy dreamplace `
        --hpwl-epsilon $Spec.Epsilon --active-power 4 --step-fraction 0.003
    if ($LASTEXITCODE -notin 0, 2) {
        throw "$Name failed with exit code $LASTEXITCODE"
    }
}

