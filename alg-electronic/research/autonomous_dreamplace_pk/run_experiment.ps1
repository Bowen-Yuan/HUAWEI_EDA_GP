param(
    [Parameter(Mandatory = $true)]
    [string]$RunName,

    [Parameter(Mandatory = $true)]
    [ValidateSet("adaptec1", "adaptec2")]
    [string]$Benchmark,

    [string]$BenchmarkStem = "",

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$SolverArgs
)

$ErrorActionPreference = "Stop"
$studyRoot = $PSScriptRoot
$placerRoot = (Resolve-Path (Join-Path $studyRoot "..\..")).Path
$executable = Join-Path $placerRoot "nsp_placer.exe"
$benchmarkStemNative = if ($BenchmarkStem) {
    [System.IO.Path]::GetFullPath($BenchmarkStem)
} else {
    Join-Path $placerRoot "ispd2005\$Benchmark\$Benchmark"
}
$benchmarkStem = $benchmarkStemNative.Replace("\", "/")
$runDir = Join-Path $studyRoot "experiments\runs\$RunName"

if (-not (Test-Path -LiteralPath $executable)) {
    throw "Missing solver executable: $executable"
}
if (-not (Test-Path -LiteralPath "$benchmarkStemNative.aux")) {
    throw "Missing benchmark: $benchmarkStemNative.aux"
}
if (Test-Path -LiteralPath $runDir) {
    throw "Run directory already exists: $runDir"
}

New-Item -ItemType Directory -Path $runDir | Out-Null
$commandLine = @($executable, $benchmarkStem) + $SolverArgs
$commandLine -join " " | Set-Content -Encoding ascii (Join-Path $runDir "command.txt")

$started = Get-Date
Push-Location $runDir
try {
    & $executable $benchmarkStem @SolverArgs *> "stdout.log"
    if ($LASTEXITCODE -ne 0) {
        throw "Solver exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
$finished = Get-Date

$generatedPlacement = "$benchmarkStemNative.nsp.pl"
if (-not (Test-Path -LiteralPath $generatedPlacement)) {
    throw "Solver did not produce placement: $generatedPlacement"
}
$runPlacement = Join-Path $runDir "placement.nsp.pl"
Copy-Item -LiteralPath $generatedPlacement -Destination $runPlacement

# legal2.pl is the official ISPD 2005 legality checker. The C++ solver performs
# legalization; this checker is the final gate for fixed-terminal movement,
# row-boundary violations, row alignment, and overlap.
$legalChecker = Join-Path $placerRoot "ispd2005\legal2.pl\legal2.pl"
if (-not (Test-Path -LiteralPath $legalChecker)) {
    throw "Missing ISPD legality checker: $legalChecker"
}
$perlCommand = Get-Command perl -ErrorAction SilentlyContinue
$perlExecutable = if ($perlCommand) {
    $perlCommand.Source
} elseif (Test-Path -LiteralPath "D:\MingGW\usr\bin\perl.exe") {
    "D:\MingGW\usr\bin\perl.exe"
} else {
    throw "Perl is required to run $legalChecker"
}

$legalStdout = Join-Path $runDir "legal2.stdout.log"
Push-Location $runDir
try {
    & $perlExecutable $legalChecker `
        "$benchmarkStemNative.nodes" `
        "$benchmarkStemNative.pl" `
        $runPlacement `
        "$benchmarkStemNative.scl" `
        20 *> $legalStdout
    if ($LASTEXITCODE -ne 0) {
        throw "legal2.pl exited with code $LASTEXITCODE; see $legalStdout"
    }
}
finally {
    Pop-Location
}

$legalityCounts = @{}
Get-Content -LiteralPath $legalStdout | ForEach-Object {
    if ($_ -match '^\s*([0-3])\s+(\d+)\s*$') {
        $legalityCounts[[int]$Matches[1]] = [int64]$Matches[2]
    }
}
if ($legalityCounts.Count -ne 4) {
    throw "Could not parse all four legal2.pl error counts; see $legalStdout"
}
$legalityFailures = 0..3 | Where-Object { $legalityCounts[$_] -ne 0 }
if ($legalityFailures.Count -ne 0) {
    $summary = ($legalityFailures | ForEach-Object {
        "Type$_=$($legalityCounts[$_])"
    }) -join ", "
    throw "Placement is not legal: $summary; see $legalStdout"
}

$metadata = [ordered]@{
    run_name = $RunName
    benchmark = $Benchmark
    started = $started.ToString("o")
    finished = $finished.ToString("o")
    solver_and_legalization_wall_time_sec = [math]::Round(($finished - $started).TotalSeconds, 3)
    exit_code = 0
    legal2_type0 = $legalityCounts[0]
    legal2_type1 = $legalityCounts[1]
    legal2_type2 = $legalityCounts[2]
    legal2_type3 = $legalityCounts[3]
}
$verified = Get-Date
$metadata.verification_finished = $verified.ToString("o")
$metadata.total_wall_time_sec = [math]::Round(($verified - $started).TotalSeconds, 3)
$legalHpwlLine = Select-String -LiteralPath (Join-Path $runDir "stdout.log") `
    -Pattern '\[Legalizer\] Exact HPWL before=([0-9.eE+-]+) after=([0-9.eE+-]+)' |
    Select-Object -Last 1
if ($legalHpwlLine -and $legalHpwlLine.Line -match `
    '\[Legalizer\] Exact HPWL before=([0-9.eE+-]+) after=([0-9.eE+-]+)') {
    $metadata.pre_legal_hpwl = [double]$Matches[1]
    $metadata.post_legal_hpwl = [double]$Matches[2]
}
$metadata | ConvertTo-Json | Set-Content -Encoding ascii (Join-Path $runDir "metadata.json")

Select-String -Path (Join-Path $runDir "stdout.log") `
    -Pattern "Restored best feasible|Done:|HPWL final=" | ForEach-Object Line
