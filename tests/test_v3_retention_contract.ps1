param(
    [string]$ExperimentDirectory,
    [string]$ExpectedSavedStage = ''
)
$ErrorActionPreference = 'Stop'
if (-not $ExperimentDirectory) { throw 'pass -ExperimentDirectory <metrics-only experiment directory>' }
$names = @(Get-ChildItem -LiteralPath $ExperimentDirectory -File | ForEach-Object Name | Sort-Object)
$expected = @('experiment.md', 'params.json', 'trajectory.csv')
if (Compare-Object $names $expected) { throw "retention contract violated: $($names -join ', ')" }
$trajectory = @(Import-Csv (Join-Path $ExperimentDirectory 'trajectory.csv'))
$hasPercentage = $trajectory.Count -gt 0 -and (
    $null -ne $trajectory[0].overflow_percent -or
    $null -ne $trajectory[0].overflow_percent_after
)
if (-not $hasPercentage) { throw 'trajectory lacks percentage overflow metrics' }

$directories = @(Get-ChildItem -LiteralPath $ExperimentDirectory -Directory)
if ($ExpectedSavedStage) {
    if ($directories.Count -ne 1 -or $directories[0].Name -ne 'saved') {
        throw 'explicit save must create only the saved directory'
    }
    $saved = @(Get-ChildItem -LiteralPath $directories[0].FullName -File)
    if ($saved.Count -ne 1 -or $saved[0].Name -ne "$ExpectedSavedStage.pl") {
        throw 'explicit save retained an unexpected placement set'
    }
} elseif ($directories.Count -ne 0) {
    throw "metrics-only experiment contains directories: $($directories.Name -join ', ')"
}

$params = Get-Content (Join-Path $ExperimentDirectory 'params.json') -Raw | ConvertFrom-Json
$checkpointPath = $params.input_checkpoint
$checkpointHash = $params.input_checkpoint_sha256
if ($checkpointPath -and $checkpointHash) {
    $actual = (Get-FileHash -LiteralPath $checkpointPath -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne $checkpointHash.ToLowerInvariant()) { throw 'external checkpoint hash changed' }
}
$workspace = Join-Path (Join-Path ([IO.Path]::GetTempPath()) 'nonsmooth-gp') $params.run_id
if (Test-Path -LiteralPath $workspace) { throw "temporary workspace was not cleaned: $workspace" }
Write-Output 'V3 retention contract passed'
