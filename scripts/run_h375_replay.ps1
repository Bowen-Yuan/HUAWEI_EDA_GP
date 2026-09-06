param(
    [string]$DatasetRoot = 'D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005',
    [int]$Threads = 16,
    [string]$ExperimentName = 'h375_full_replay'
)

$ErrorActionPreference = 'Stop'
$project = Split-Path -Parent $PSScriptRoot
$core = Join-Path $project 'bin\nsgp_legacy_stage.exe'
$homotopy = Join-Path $project 'bin\nsgp_historical_homotopy.exe'
$benchmark = Join-Path (Join-Path $DatasetRoot 'adaptec1') 'adaptec1'
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$root = Join-Path $project "framework\experiment_logs\${stamp}_${ExperimentName}"
New-Item -ItemType Directory -Force -Path $root | Out-Null
if (!(Test-Path "$benchmark.aux")) { throw "Missing benchmark: $benchmark" }
if (!(Test-Path $core) -or !(Test-Path $homotopy)) { throw 'Build the local stage and homotopy executables first.' }

function Run-Stage([string]$Name, [string]$Exe, [string[]]$StageArgs) {
    $stage = Join-Path $root $Name
    New-Item -ItemType Directory -Force -Path $stage | Out-Null
    $log = Join-Path $stage 'stdout_stderr.log'
    @{ name=$Name; executable=$Exe; arguments=$StageArgs; started_at=(Get-Date).ToString('o') } |
        ConvertTo-Json -Depth 4 | Set-Content -Encoding utf8 (Join-Path $stage 'invocation.json')
    & $Exe @StageArgs *>> $log
    if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne 2) { throw "$Name failed ($LASTEXITCODE): $log" }
    return $stage
}

# The explicit stage list is the historical H375 chain. historical DCT/Poisson
# stages are deliberately marked in the root manifest rather than being hidden.
$s219 = Run-Stage 'h219_hpwl_seed' $core @('--benchmark',$benchmark,'--output',(Join-Path $root 'h219_hpwl_seed'),'--bins','512','--iterations','200','--threads',$Threads,'--optimizer','adam','--lambda-policy','dreamplace','--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--step-fraction','0.003','--seed','219','--sigma-ratio','0.001','--log-every','1')
$c219 = Run-Stage 'h219_dct_poisson_128' $homotopy @('--benchmark',$benchmark,'--initial-placement',(Join-Path $s219 'global.pl'),'--output',(Join-Path $root 'h219_dct_poisson_128'),'--bins','128','--stages','24','--stage-iterations','50','--snapshot-every','25','--threads',$Threads,'--mu-start','5.36','--mu-decay','0.99','--lambda-overlap','1','--hpwl-epsilon','125','--hpwl-active-power','1','--adaptive-weights','--gradient-balanced-weights','--adaptive-warmup-iterations','100','--adaptive-update-interval','10','--stability-patience','1000','--adaptive-mu-max','5.36','--adaptive-mu-handoff-overflow','0.6','--adaptive-handoff-on-overflow','--adaptive-lambda-step','15','--adaptive-lambda-max','150','--adaptive-final-lambda-band','--select-feasible-hpwl','--anchor-weight','0.0001','--anchor-mode','centroid','--step-fraction','0.002','--seed','219')
$f221 = Run-Stage 'h221_dct_poisson_512' $homotopy @('--benchmark',$benchmark,'--initial-placement',(Join-Path $c219 'snapshots\iter_1100.pl'),'--output',(Join-Path $root 'h221_dct_poisson_512'),'--bins','512','--stages','8','--stage-iterations','50','--snapshot-every','25','--threads',$Threads,'--mu-start','964','--mu-decay','0.99','--lambda-overlap','1','--hpwl-epsilon','125','--hpwl-active-power','1','--adaptive-weights','--gradient-balanced-weights','--adaptive-warmup-iterations','0','--adaptive-update-interval','10','--stability-patience','1000','--adaptive-mu-max','964','--adaptive-mu-handoff-overflow','0.35','--adaptive-handoff-on-overflow','--adaptive-lambda-step','520','--adaptive-lambda-max','1950','--adaptive-final-lambda-band','--select-feasible-hpwl','--anchor-weight','0.00001','--anchor-mode','initial','--step-fraction','0.002','--seed','221')
$h252 = Run-Stage 'h252_cap15_recovery' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $f221 'best.pl'),'--output',(Join-Path $root 'h252_cap15_recovery'),'--bins','512','--iterations','1','--threads',$Threads,'--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--stop-overflow','0.15','--recovery-sweeps','2','--recovery-line-search','4','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-net-blocks','--recovery-net-block-density-direction','--recovery-net-block-degree-limit','32','--recovery-net-block-max-nodes','16','--recovery-net-block-max-blocks','10000','--recovery-overflow-cap','0.15')
$h253 = Run-Stage 'h253_bridge' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $h252 'global.pl'),'--output',(Join-Path $root 'h253_bridge'),'--bins','512','--iterations','81','--threads',$Threads,'--density-weight-scale','0.5','--net-batch-weight','1','--net-batch-degree-limit','32','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--step-fraction','0.002','--stop-overflow','0.07')
$h254 = Run-Stage 'h254_retighten' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $h253 'global.pl'),'--output',(Join-Path $root 'h254_retighten'),'--bins','512','--iterations','120','--threads',$Threads,'--density-weight-scale','4.75','--net-batch-weight','1','--net-batch-degree-limit','32','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--step-fraction','0.001','--stop-overflow','0.07')
$h255 = Run-Stage 'h255_polish' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $h254 'global.pl'),'--output',(Join-Path $root 'h255_polish'),'--bins','512','--iterations','1','--threads',$Threads,'--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recovery-sweeps','4','--recovery-line-search','6','--recovery-degree-limit','100','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-overflow-cap','0.07')
$h257 = Run-Stage 'h257_polish' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $h255 'global.pl'),'--output',(Join-Path $root 'h257_polish'),'--bins','512','--iterations','1','--threads',$Threads,'--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recovery-sweeps','1','--recovery-line-search','4','--recovery-degree-limit','100','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-overflow-cap','0.07')
$b372 = Run-Stage 'h372_surplus_bisection' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $h257 'global.pl'),'--output',(Join-Path $root 'h372_surplus_bisection'),'--bins','512','--iterations','1','--threads','40','--density-weight-scale','0','--hpwl-epsilon','110','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recursive-bisection','--bisection-leaf-bins','8','--bisection-degree-limit','32','--bisection-fm-passes','1','--bisection-position-seeded','--bisection-surplus-only','--bisection-leaf-nearest-capacity','--bisection-leaf-hpwl-guided')
$r372 = Run-Stage 'h372_recovery10' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $b372 'global.pl'),'--output',(Join-Path $root 'h372_recovery10'),'--bins','512','--iterations','1','--threads','40','--density-weight-scale','0','--hpwl-epsilon','85','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recovery-sweeps','10','--recovery-line-search','8','--recovery-degree-limit','16','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-net-blocks','--recovery-net-block-density-direction','--recovery-net-block-contraction','--recovery-net-block-degree-limit','16','--recovery-net-block-max-nodes','8','--recovery-net-block-max-blocks','10000','--recovery-overflow-cap','0.07')
$h375 = Run-Stage 'h375_equal_shape_swap5' $core @('--benchmark',$benchmark,'--initial-placement',(Join-Path $r372 'global.pl'),'--output',(Join-Path $root 'h375_equal_shape_swap5'),'--bins','512','--iterations','1','--threads','40','--density-weight-scale','0','--hpwl-epsilon','85','--active-power','4','--density-epsilon','0','--fixed-epsilon','--swap-sweeps','5','--swap-radius-bins','64','--swap-candidates','32','--swap-exact-shortlist','16','--swap-degree-limit','1000','--swap-net-aware','--stop-overflow','0.07')

$summary = Get-Content (Join-Path $h375 'summary.txt')
$summary | Set-Content -Encoding utf8 (Join-Path $root 'final_summary.txt')
$values = @{}
foreach ($line in $summary) { if ($line -match '^([^=]+)=(.+)$') { $values[$Matches[1]] = $Matches[2] } }
$finalHpwl = [double]$values['gp_hpwl']
$finalOverflowPercent = 100.0 * [double]$values['gp_overflow']
@{ experiment_name=$ExperimentName; dataset='adaptec1'; dataset_initialization='raw Bookshelf .pl'; chain=@('h219_hpwl_seed','h219_dct_poisson_128','h221_dct_poisson_512','h252_cap15_recovery','h253_bridge','h254_retighten','h255_polish','h257_polish','h372_surplus_bisection','h372_recovery10','h375_equal_shape_swap5'); uses_historical_smooth_surrogate=$true; overflow_display='percent (raw ratio x 100)'; final_stage='h375_equal_shape_swap5'; final_hpwl=$finalHpwl; final_overflow_percent=$finalOverflowPercent; finished_at=(Get-Date).ToString('o') } | ConvertTo-Json -Depth 5 | Set-Content -Encoding utf8 (Join-Path $root 'experiment_manifest.json')
@"
# Experiment Summary

Experiment: $ExperimentName
Dataset initialization: raw Bookshelf placement
Final HPWL: $finalHpwl
Final overflow: $([Math]::Round($finalOverflowPercent, 8))%
Historical smooth surrogate: yes (H219/H221 DCT/Poisson only)

The full stage invocation, input checkpoint chain, stdout/stderr, and selected checkpoints are in sibling stage folders.
"@ | Set-Content -Encoding utf8 (Join-Path $root 'run_summary.md')
Write-Host "H375 replay complete: $root"
