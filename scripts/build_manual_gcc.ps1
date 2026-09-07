# Manual MinGW build for machines without CMake.  This script mirrors the
# CMakeLists.txt targets (nsgp, nsgp_tests, nsgp_adaptive_tests) one to one;
# CMakeLists.txt remains the authoritative build description.  Run from the
# repository root:
#   powershell -ExecutionPolicy Bypass -File scripts/build_manual_gcc.ps1
param(
    [string]$Gpp = 'D:/MingGW/ucrt64/bin/g++.exe',
    [string]$GitCommit = 'unknown',
    [string]$GitBranch = 'unknown',
    [string]$GitRemote = 'unknown',
    [string]$GitDirty = 'not_checked'
)
$ErrorActionPreference = 'Stop'

# PowerShell mangles embedded quotes in -D macros, so the Git metadata is
# injected through a generated forced-include header instead.  Every source
# file guards these macros with #ifndef, matching the CMake definitions.
$metadataHeader = 'build/nsgp_git_metadata.h'
$metadataLines = @(
  '#pragma once',
  ('#define NSGP_GIT_COMMIT "' + $GitCommit + '"'),
  ('#define NSGP_GIT_BRANCH "' + $GitBranch + '"'),
  ('#define NSGP_GIT_DIRTY "' + $GitDirty + '"'),
  ('#define NSGP_GIT_REMOTE "' + $GitRemote + '"'))
New-Item -ItemType Directory -Force -Path 'build' | Out-Null
Set-Content -Path $metadataHeader -Value $metadataLines -Encoding ascii

$flags = @('-std=c++17','-O3','-Wall','-Wextra','-Wpedantic','-fopenmp',
  '-Iframework/kernel/include','-Iframework/code','-Ithird_party',
  '-Imodules/global_view_gp/code','-Imodules/adaptive_lambda_gp/code',
  '-include', $metadataHeader)

function Convert-ToObjectPath([string]$Source) {
  $dotted = $Source.Replace('/', '_').Replace('.cpp', '.o')
  return (Join-Path 'build/obj' $dotted)
}

$kernelSources = @(
  'framework/kernel/src/bookshelf.cpp',
  'framework/kernel/src/density.cpp',
  'framework/kernel/src/hpwl.cpp',
  'framework/kernel/src/optimizer.cpp',
  'framework/kernel/src/step_policy.cpp',
  'modules/exact_joint_gp/code/batch_acceptance.cpp',
  'modules/exact_joint_gp/code/lambda_controller.cpp',
  'modules/exact_joint_gp/code/placer.cpp',
  'modules/exact_recovery/code/compact_recovery.cpp',
  'modules/exact_recovery/code/recovery.cpp',
  'modules/surplus_bisection/code/bisection.cpp',
  'modules/equal_shape_swap/code/swap_recovery.cpp',
  'modules/global_capacity_transport/code/coarse_flow.cpp',
  'modules/global_capacity_transport/code/transport.cpp',
  'modules/density_coordinate/code/density_coordinate.cpp')

$nsgpSources = @(
  'framework/code/main.cpp',
  'framework/code/microkernel.cpp',
  'modules/global_view_gp/code/global_view_lab.cpp',
  'modules/layout_init/code/module.cpp',
  'modules/hpwl_adam/code/module.cpp',
  'modules/exact_joint_gp/code/module.cpp',
  'modules/exact_recovery/code/module.cpp',
  'modules/surplus_bisection/code/module.cpp',
  'modules/equal_shape_swap/code/module.cpp',
  'modules/global_capacity_transport/code/module.cpp',
  'modules/density_coordinate/code/module.cpp',
  'modules/adaptive_lambda_gp/code/adaptive_lambda_gp.cpp',
  'modules/adaptive_lambda_gp/code/module.cpp')

New-Item -ItemType Directory -Force -Path 'build/obj' | Out-Null

$allSources = $kernelSources + $nsgpSources
foreach ($source in $allSources) {
  $objectPath = Convert-ToObjectPath $source
  & $Gpp @flags -c $source -o $objectPath
  if ($LASTEXITCODE -ne 0) { throw "compile failed: $source" }
}
Write-Output 'objects compiled'

$allObjects = $allSources | ForEach-Object { Convert-ToObjectPath $_ }
$kernelObjects = $kernelSources | ForEach-Object { Convert-ToObjectPath $_ }

& $Gpp -fopenmp @allObjects -o build/nsgp.exe -ladvapi32
if ($LASTEXITCODE -ne 0) { throw 'link failed: nsgp' }
& $Gpp @flags tests/test_numeric_contracts.cpp @kernelObjects -o build/nsgp_tests.exe
if ($LASTEXITCODE -ne 0) { throw 'link failed: nsgp_tests' }
$adaptiveObject = Convert-ToObjectPath 'modules/adaptive_lambda_gp/code/adaptive_lambda_gp.cpp'
& $Gpp @flags tests/test_adaptive_lambda_contracts.cpp $adaptiveObject @kernelObjects -o build/nsgp_adaptive_tests.exe
if ($LASTEXITCODE -ne 0) { throw 'link failed: nsgp_adaptive_tests' }
Write-Output 'ALL TARGETS BUILT'
