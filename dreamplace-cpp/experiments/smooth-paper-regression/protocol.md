# Smooth DREAMPlace Paper Regression Protocol

## Question

Does the current default smooth C++ path still reproduce the expected
DREAMPlace-style behavior on ISPD 2005 adaptec1 through adaptec4 after the
exact-HPWL and progressive-legalization extensions?

## Locked Configuration

- Wirelength: weighted-average (the default smooth path).
- Optimizer: DREAMPlace/RePlAce Nesterov with BB step estimation.
- Density: electrostatic DCT Poisson model with fillers.
- Initialization: center Gaussian, seed 1000, sigma ratio 0.001 in this C++
  reconstruction.
- Iterations: 1000 for every design.
- Target density: 1.0.
- Stop-overflow checkpoint: 0.07.
- Net-degree cutoff: 100.
- Grid: 512x512 for adaptec1; 1024x1024 for adaptec2, adaptec3, adaptec4.
- No exact-HPWL, progressive-legalization, bundle, gradient-sampling,
  multilevel-density, or legal-projection switches.
- Global snapshots every 25 iterations. All movable cells will be rendered.

The grids and 1000-step schedule are taken from the official DREAMPlace
`test/ispd2005/adaptec*.json` configurations. The C++ reconstruction does not
claim bitwise equivalence to upstream CUDA/PyTorch; this is a regression of
the reconstructed smooth algorithmic path.

## Paper Reference

DREAMPlace Table II final HPWL values (float64, after the paper's full
placement flow):

| design | paper final HPWL |
|---|---:|
| adaptec1 | 73.22M |
| adaptec2 | 82.22M |
| adaptec3 | 193.72M |
| adaptec4 | 174.08M |

These are final post-DP values. Pre-legal GP HPWL is reported separately and
must not be compared directly against Table II as though it were final HPWL.

## Outputs

Each run must contain `global_metrics.csv`, `global.pl`, `final.pl`,
`summary.txt`, global snapshots, and legal-stage snapshots. The aggregate
study must include per-design HPWL/overflow curves, a cross-design normalized
curve, and one all-cell convergence/legalization animation per design.
