param([string]$ExperimentDirectory)
$ErrorActionPreference = 'Stop'
if (-not $ExperimentDirectory) { throw 'pass -ExperimentDirectory <metrics-only experiment directory>' }
$names = @(Get-ChildItem -LiteralPath $ExperimentDirectory -File | ForEach-Object Name | Sort-Object)
$expected = @('experiment.md', 'params.json', 'trajectory.csv')
if (Compare-Object $names $expected) { throw "retention contract violated: $($names -join ', ')" }
$trajectory = Import-Csv (Join-Path $ExperimentDirectory 'trajectory.csv')
if ($trajectory.Count -lt 1 -or -not $trajectory[0].overflow_percent) { throw 'trajectory lacks percentage overflow metrics' }
Write-Output 'V3 retention contract passed'
