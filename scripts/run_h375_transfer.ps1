param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('adaptec1', 'adaptec2', 'adaptec3', 'adaptec4')]
    [string]$Benchmark,
    [string]$DatasetRoot = 'D:\codex_project\HUAWEI_EDA\alg\ispd2005',
    [string]$OutputRoot = 'D:\codex_project\HUAWEI_EDA\epsilon-active\output\h375_transfer_20260830'
)

$ErrorActionPreference = 'Stop'
$project = Split-Path -Parent $PSScriptRoot
$core = Join-Path $project 'build\epsilon_active.exe'
$homotopy = Join-Path $project 'experiments\h72_homotopy_electrostatic\homotopy.exe'
$benchmarkBase = Join-Path (Join-Path $DatasetRoot $Benchmark) $Benchmark
$root = Join-Path $OutputRoot $Benchmark

if (!(Test-Path -LiteralPath $core)) { throw "Missing core executable: $core" }
if (!(Test-Path -LiteralPath $homotopy)) { throw "Missing archived homotopy executable: $homotopy" }
if (!(Test-Path -LiteralPath "$benchmarkBase.aux")) { throw "Missing benchmark: $benchmarkBase" }

New-Item -ItemType Directory -Force -Path $root | Out-Null
function Run([string]$label, [string]$exe, [string[]]$arguments) {
    $log = Join-Path $root "$label.log"
    "COMMAND: $exe $($arguments -join ' ')" | Set-Content -Encoding utf8 $log
    & $exe @arguments *>> $log
    # The core executable returns 2 when the selected checkpoint is infeasible;
    # its global.pl/summary are still valid inputs to the scheduled next stage.
    if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne 2) {
        throw "$label failed with exit code $LASTEXITCODE; see $log"
    }
}

# H219: reproducible fresh HPWL seed.
$s219 = Join-Path $root 'h219_seed'
Run 'h219_seed' $core @('--benchmark',$benchmarkBase,'--output',$s219,'--bins','512','--iterations','200','--threads','16','--optimizer','adam','--lambda-policy','dreamplace','--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--step-fraction','0.003','--seed','219','--sigma-ratio','0.001','--log-every','1')

# H219 coarse electrostatic homotopy.  The selected continuation is iter_1100.
$c219 = Join-Path $root 'h219_coarse'
Run 'h219_coarse' $homotopy @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $s219 'global.pl'),'--output',$c219,'--bins','128','--stages','24','--stage-iterations','50','--snapshot-every','25','--threads','16','--mu-start','5.36','--mu-decay','0.99','--lambda-overlap','1','--hpwl-epsilon','125','--hpwl-active-power','1','--adaptive-weights','--gradient-balanced-weights','--adaptive-warmup-iterations','100','--adaptive-update-interval','10','--stability-patience','1000','--adaptive-mu-max','5.36','--adaptive-mu-handoff-overflow','0.6','--adaptive-handoff-on-overflow','--adaptive-lambda-step','15','--adaptive-lambda-max','150','--adaptive-final-lambda-band','--select-feasible-hpwl','--anchor-weight','0.0001','--anchor-mode','centroid','--step-fraction','0.002','--seed','219')

$coarseCheckpoint = Join-Path $c219 'snapshots\iter_1100.pl'
if (!(Test-Path -LiteralPath $coarseCheckpoint)) { throw "Missing H219 continuation checkpoint: $coarseCheckpoint" }

# H221 fine electrostatic homotopy.
$f221 = Join-Path $root 'h221_fine'
Run 'h221_fine' $homotopy @('--benchmark',$benchmarkBase,'--initial-placement',$coarseCheckpoint,'--output',$f221,'--bins','512','--stages','8','--stage-iterations','50','--snapshot-every','25','--threads','16','--mu-start','964','--mu-decay','0.99','--lambda-overlap','1','--hpwl-epsilon','125','--hpwl-active-power','1','--adaptive-weights','--gradient-balanced-weights','--adaptive-warmup-iterations','0','--adaptive-update-interval','10','--stability-patience','1000','--adaptive-mu-max','964','--adaptive-mu-handoff-overflow','0.35','--adaptive-handoff-on-overflow','--adaptive-lambda-step','520','--adaptive-lambda-max','1950','--adaptive-final-lambda-band','--select-feasible-hpwl','--anchor-weight','0.00001','--anchor-mode','initial','--step-fraction','0.002','--seed','221')

# H252--H257: exact recovery, bridge, retighten, and 7% polish.
$h252 = Join-Path $root 'h252'; Run 'h252' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $f221 'best.pl'),'--output',$h252,'--bins','512','--iterations','1','--threads','16','--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--stop-overflow','0.15','--recovery-sweeps','2','--recovery-line-search','4','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-net-blocks','--recovery-net-block-density-direction','--recovery-net-block-degree-limit','32','--recovery-net-block-max-nodes','16','--recovery-net-block-max-blocks','10000','--recovery-overflow-cap','0.15')
$h253 = Join-Path $root 'h253'; Run 'h253' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $h252 'global.pl'),'--output',$h253,'--bins','512','--iterations','81','--threads','16','--density-weight-scale','0.5','--net-batch-weight','1','--net-batch-degree-limit','32','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--step-fraction','0.002','--stop-overflow','0.07')
$h254 = Join-Path $root 'h254'; Run 'h254' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $h253 'global.pl'),'--output',$h254,'--bins','512','--iterations','120','--threads','16','--density-weight-scale','4.75','--net-batch-weight','1','--net-batch-degree-limit','32','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--step-fraction','0.001','--stop-overflow','0.07')
$h255 = Join-Path $root 'h255'; Run 'h255' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $h254 'global.pl'),'--output',$h255,'--bins','512','--iterations','1','--threads','16','--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recovery-sweeps','4','--recovery-line-search','6','--recovery-degree-limit','100','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-overflow-cap','0.07')
$h257 = Join-Path $root 'h257'; Run 'h257' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $h255 'global.pl'),'--output',$h257,'--bins','512','--iterations','1','--threads','16','--density-weight-scale','0','--hpwl-epsilon','125','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recovery-sweeps','1','--recovery-line-search','4','--recovery-degree-limit','100','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-overflow-cap','0.07')

# H372: surplus-only capacity assignment, then ten-sweep exact recovery.
$b372 = Join-Path $root 'h372_bisect'; Run 'h372_bisect' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $h257 'global.pl'),'--output',$b372,'--bins','512','--iterations','1','--threads','40','--density-weight-scale','0','--hpwl-epsilon','110','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recursive-bisection','--bisection-leaf-bins','8','--bisection-degree-limit','32','--bisection-fm-passes','1','--bisection-position-seeded','--bisection-surplus-only','--bisection-leaf-nearest-capacity','--bisection-leaf-hpwl-guided')
$r372 = Join-Path $root 'h372_recovery10'; Run 'h372_recovery10' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $b372 'global.pl'),'--output',$r372,'--bins','512','--iterations','1','--threads','40','--density-weight-scale','0','--hpwl-epsilon','85','--active-power','4','--density-epsilon','0','--fixed-epsilon','--recovery-sweeps','10','--recovery-line-search','8','--recovery-degree-limit','16','--recovery-axis-separated','--recovery-breakpoint-oracle','--recovery-compact-directions','--recovery-net-blocks','--recovery-net-block-density-direction','--recovery-net-block-contraction','--recovery-net-block-degree-limit','16','--recovery-net-block-max-nodes','8','--recovery-net-block-max-blocks','10000','--recovery-overflow-cap','0.07')

# H375: density-neutral exact equal-shape exchange.
$h375 = Join-Path $root 'h375_exchange'; Run 'h375_exchange' $core @('--benchmark',$benchmarkBase,'--initial-placement',(Join-Path $r372 'global.pl'),'--output',$h375,'--bins','512','--iterations','1','--threads','40','--density-weight-scale','0','--hpwl-epsilon','85','--active-power','4','--density-epsilon','0','--fixed-epsilon','--swap-sweeps','5','--swap-radius-bins','64','--swap-candidates','32','--swap-exact-shortlist','16','--swap-degree-limit','1000','--swap-net-aware','--stop-overflow','0.07')

$final = Get-Content -LiteralPath (Join-Path $h375 'summary.txt')
$final | Set-Content -Encoding utf8 (Join-Path $root 'final_summary.txt')
Write-Host "H375 transfer replay finished: $root"
