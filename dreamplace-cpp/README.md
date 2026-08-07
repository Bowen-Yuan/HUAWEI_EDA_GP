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
