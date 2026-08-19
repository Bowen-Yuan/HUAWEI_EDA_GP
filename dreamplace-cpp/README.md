# DREAMPlace C++

Independent, CPU-oriented C++17 reconstruction of the DREAMPlace placement
flow.  It reads raw ISPD 2005 Bookshelf data, performs center-Gaussian global
placement, legalizes onto obstacle-aware row segments, and applies legal
detailed-placement refinement.

The project does not link against `alg-electronic` and does not consume ePlace
or prior NSP placement artifacts.  The only accepted placement input is the
raw `<benchmark>.pl` named by the benchmark base path.  Any path containing
`eplace`, `nsp`, or `gp.pl` is rejected.

Implemented algorithm:

1. Center initialization with 0.1% Gaussian noise (paper Section III).
2. Weighted-average wirelength with exact HPWL used for reporting.
3. Stretched-cell electrostatic density on a DCT Poisson grid.
4. Gradient-norm density-weight initialization and RePlAce/DREAMPlace update.
5. Nesterov accelerated gradient with Barzilai-Borwein step estimation.
6. Fixed-obstacle row segmentation and left-to-right Tetris-like Greedy
   legalization followed by segment-aware Abacus isotonic compaction.
7. Exact pin-offset-HPWL adjacent reorder that preserves legal intervals.
8. Optional CPU DREAMPlace-style K-Reorder, equal-size global swap, and
   independent-set slot matching; every move is accepted only after an exact
   affected-net HPWL check and a full legality check follows every stage.
9. Bookshelf output suitable for `ispd2005/legal2.pl`.

The original smooth flow remains the default.  An opt-in migration layer also
supports the exact nonsmooth objective

```text
exact HPWL (including pin offsets) + lambda * electrostatic density
```

with Heavy Ball, Adam, AMSGrad, or AdaGrad; DREAMPlace-style, trajectory, or
ratio lambda control; and optional HPWL-only bundle aggregation.  Density is
never included in the bundle because its weight changes during optimization.
The exact mode is selected explicitly, for example:

```powershell
.\dreamplace_cpp.exe ..\alg-electronic\ispd2005\adaptec1\adaptec1 `
  --exact-hpwl --optimizer adam --lambda-policy dreamplace `
  --bins 64 --iterations 800 --output output\exact_adam
```

Use `--bundle-hpwl` with `--lambda-policy trajectory` to enable the migrated
bundle-inspired direction.  Existing commands without these flags still use
weighted-average wirelength and DREAMPlace BB-Nesterov.

An opt-in two-stage experiment can reuse a placement produced by this program
without changing the default center initialization.  `--initial-pl FILE`
loads movable coordinates only, keeps fixed objects at their raw Bookshelf
locations, requires every movable node to be present, and disables the exact
solver's additional random perturbation.  For example:

```powershell
.\dreamplace_cpp.exe ..\alg-electronic\ispd2005\adaptec1\adaptec1 `
  --initial-pl output\adaptec1_papergrid_600\global.pl `
  --exact-hpwl --bundle-hpwl --serious-bundle --feasible-refinement `
  --bins 512 --iterations 200 --output output\smooth_warm_exact_bundle
```

The post-feasible controller and the extended detailed placer are also opt-in:

```powershell
.\dreamplace_cpp.exe ..\alg-electronic\ispd2005\adaptec1\adaptec1 `
  --exact-hpwl --optimizer adam --bins 512 --iterations 800 `
  --feasible-refinement --bundle-hpwl --bundle-start-overflow 0.12 `
  --dreamplace-detailed --output output\exact_512_refine_bundle
```

`--feasible-refinement` resets Adam's moments on the first state at or below
the stop-overflow threshold, switches to AMSGrad, reduces the step, and uses a
bounded lambda feedback controller over the 6.5%-7.0% band.  It never replaces
exact HPWL by a smooth surrogate.  `--dreamplace-detailed` is a CPU algorithmic
counterpart of the upstream operator sequence, not a line-for-line port of the
CUDA kernels.

Additional legalization-aware research controls are opt-in and therefore
reversible without changing the command format.  For example, the strongest
detailed-placement ablation on the validated adaptec1 checkpoint used:

```powershell
--dreamplace-detailed --legal-refine-rounds 3 `
  --cell-insertion-passes 1 --cell-insertion-window 16 `
  --legal-projected-passes 2 --legal-projected-step-sites 4 `
  --row-relegalization-passes 1 `
  --independent-set-size 8 --hungarian-matching
```

Mixed-spectral fields, multilevel density, exact gradient sampling,
serious/null bundle trials, legal checkpoint ranking, and late legal
projection are also available through `--help`.  They remain disabled by
default.  Removing the added flags restores the previous path; the full
implementation and negative as well as positive ablations are documented in
`report/legalization_aware_exact_hpwl_extensions_report.pdf`.

The epsilon-active direction uses a cached flat pin topology internally. This
is an implementation-only optimization and preserves the existing public
function and command-line behavior. An experimental adaptive radius controller
can be enabled with `--adaptive-active-set`; it is disabled by default and is
documented, together with a1/a2 measurements, in
`experiments/epsilon-active-optimization-20260813/analysis.md`.

The independent `dreamplace-cpp-relative-epsilon` experiment can constrain
the trial-direction radius separately for each net and coordinate axis:

```text
epsilon(e,d) = min(global_radius, span_ratio * span(e,d))
```

Use `--active-set-span-ratio 0.10` or `0.15` to enable it. A nonzero
`--active-set-min-radius` is applied only when the directional span is below
`--active-set-small-span-threshold`; this is the sole case in which epsilon
may exceed the selected span ratio. All three options default to zero, so
existing commands retain their original behavior.

The smooth-to-exact continuation study has two additional opt-in components.
`--bundle-groups N` partitions nets by stable net index and keeps an
independent short HPWL bundle for each group; nearby exact-subgradient samples
are blended into the trial direction, while exact HPWL and the overflow filter
remain authoritative.  `--bundle-overflow-tolerance X` is only a trial band:
a step crossing 7% from a feasible center is bisected back to the strict
boundary, so saved feasible checkpoints still obey `--stop-overflow`.

After legalization, `--legal-bundle-passes N` enables a fixed-row,
fixed-obstacle-segment, fixed-order x-only proximal bundle.  Every trial is
projected and site-snapped by the segment Abacus operator, then accepted only
after strict exact-HPWL improvement and a full legality check.  The feature is
disabled by default and does not change existing detailed-placement commands.
The strongest adaptec1 exploratory command was:

```powershell
.\dreamplace_cpp.exe ..\alg-electronic\ispd2005\adaptec1\adaptec1 `
  --initial-pl output\adaptec1_papergrid_600\global.pl `
  --exact-hpwl --optimizer amsgrad --bins 512 --iterations 1 `
  --feasible-refinement --refine-lr-scale 0.01 `
  --dreamplace-detailed --legal-refine-rounds 5 `
  --cell-insertion-passes 1 --legal-projected-passes 2 `
  --independent-set-size 8 --hungarian-matching `
  --legal-bundle-passes 16 --legal-bundle-size 12 `
  --legal-bundle-step-sites 4 `
  --output output\smooth_warm_legal_bundle_control_v4
```

It produced a legal 74.425850M placement.  A matched five-round command with
the legal bundle disabled produced 74.520409M, so the bundle's measured net
improvement was 94.559K (0.127%).  A 20-step pre-legal grouped-bundle run
reduced GP HPWL from 71.235743M to 71.219645M at 7%, but worsened the matched
three-round final legal result by 9.103K.  The grouped pre-legal stage is
therefore retained as an experimental negative result, not enabled in the
best pipeline.

## Progressive exact-HPWL legalization experiment

`--progressive-legalization` enables a separate, default-off five-stage
research path. It never replaces exact HPWL by a weighted-average or
log-sum-exp surrogate:

1. a short exact-HPWL-only stage;
2. exact HPWL plus electrostatic spreading with hybrid density-weight control;
3. a continuous ramp of fixed-macro escape, row proximity, and obstacle-free
   segment-capacity forces;
4. the existing Greedy/Abacus discrete legalizer; and
5. DREAMPlace-style legal detailed placement plus optional fixed-topology
   legal bundle refinement.

Stage 3 does not assign a row, segment, site, or cell order. Fixed-macro
overlap is penalized by the squared distance to the nearest escape face, and
row/segment terms are continuous directions only. Discrete topology changes
remain in stages 4 and 5. The hybrid lambda controller uses the original
DREAMPlace trajectory while overflow is above 11%, then treats lambda as a
band-dual variable for the configured 7%-8% target. On the first entry into
the target band, Adam moments are cleared, the optimizer switches to AMSGrad,
and the learning rate is reduced.

The strongest center-initialized exact-HPWL adaptec1 run currently gives
83.254M at 7.744% overflow before legalization and 88.895M after detailed
placement. The external ISPD `legal2.pl` checker reports Type 0/1/2/3 =
0/0/0/0. This improves the previous center-initialized 512x512 exact result
of 85.401M / 6.958% and 91.030M legal, but it does not reach the requested
75-80M pre-legal range. The result and all negative ablations are documented
in `report/progressive_legalization_exact_hpwl_report.pdf`.

Reproduction command:

```powershell
.\dreamplace_cpp.exe ..\alg-electronic\ispd2005\adaptec1\adaptec1 `
  --progressive-legalization --optimizer adam `
  --bins 512 --iterations 1000 `
  --hpwl-only 100 --progressive-density-iterations 475 `
  --progressive-obstacle-iterations 200 `
  --hpwl-step-fraction 0.004 --step-fraction 0.002 `
  --refine-lr-scale 0.05 `
  --progressive-hpwl-target 75000000 `
  --progressive-overflow-lower 0.07 `
  --progressive-overflow-upper 0.08 `
  --progressive-obstacle-scale 0.2 `
  --progressive-row-force 0.03 `
  --progressive-segment-force 0.05 `
  --progressive-congestion-gain 2.0 `
  --snapshot-every 25 `
  --dreamplace-detailed --dp-passes 1 `
  --legal-refine-rounds 3 `
  --cell-insertion-passes 1 --cell-insertion-window 16 `
  --legal-projected-passes 2 --legal-projected-step-sites 4 `
  --row-relegalization-passes 1 `
  --independent-set-size 8 --hungarian-matching `
  --legal-bundle-passes 12 --legal-bundle-size 12 `
  --legal-bundle-step-sites 4 `
  --output output\progressive_a1_512_1000_refine
```

Removing `--progressive-legalization` restores the prior solver path. The
progressive-only CSV columns are appended only when this switch is active, so
existing parsers see their original schema for old commands.

Build with MinGW:

```powershell
D:\MingGW\ucrt64\bin\mingw32-make.exe
```

Run from the project directory:

```powershell
.\dreamplace_cpp.exe ..\alg-electronic\ispd2005\adaptec1\adaptec1
```

The default output is written below `output/<benchmark>/`.  Use `--help` for
configuration flags.  The paper-reference settings are 512 bins for adaptec1,
1024 bins for adaptec2, 1000 global iterations, target density 1.0, and stop
overflow 0.07.

## Validated results

All runs below start from the paper's center Gaussian initialization with seed
1000.  The final column compares legal post-refinement HPWL against DREAMPlace
Table II, not against the lower pre-legal GP number.

| case | grid / steps | selected GP (HPWL / overflow) | legal HPWL | paper | gap | legal2 Type 0/1/2/3 |
|---|---:|---:|---:|---:|---:|---:|
| adaptec1 | 512x512 / 600 | 71.236M / 6.965% | 74.819M | 73.220M | +2.18% | 0 / 0 / 0 / 0 |
| adaptec2 | 1024x1024 / 650 | 80.889M / 6.930% | 84.427M | 82.220M | +2.68% | 0 / 0 / 0 / 0 |

Artifacts are in `output/adaptec1_papergrid_600` and
`output/adaptec2_papergrid_650`.  `global_metrics.csv` contains every GP
iteration.  `global.pl` is the selected pre-legal checkpoint and `final.pl` is
the externally validated legal result.

The implementation is not a byte-for-byte port.  In particular, the electric
field currently differentiates the DCT potential with finite differences,
while upstream DREAMPlace uses mixed sine/cosine spectral transforms.  The
detailed-placement operators are CPU counterparts and remain smaller than the
upstream GPU implementations; the paper's NTUplace3 final stage is not
included.  These remain quality and runtime opportunities.

## Nonsmooth migration status

At 64x64 and 800 iterations, exact HPWL + Adam + DREAMPlace lambda selected
90.069M at 6.952% overflow and produced a legal 114.017M placement.  The
official checker reported Type 0/1/2/3 = 0/0/0/0.  Trajectory lambda plus four
HPWL bundle cuts also reached 6.824% overflow, but its legal HPWL was worse at
118.153M.  These runs validate the migrated interfaces and legalization, not
paper-level quality; the exact mode still trails the smooth baseline and needs
optimizer/step-size refinement.

The high-resolution exact-HPWL ablation improves this substantially.  At
512x512 and 800 iterations, Adam plus the post-feasible bundle/refinement path
selected 85.401M at 6.958% overflow and produced a legal 91.030M result after
the extended detailed placer.  The official checker reported Type 0/1/2/3 =
0/0/0/0.  The corresponding artifacts are in
`output/exact_adam_512_800_refine_bundle_dp`; the full ablation and limitations
are documented in `report/high_resolution_detailed_placement_report.pdf`.
