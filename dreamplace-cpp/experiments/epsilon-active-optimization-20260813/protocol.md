# Epsilon-active optimization experiment

Date: 2026-08-13

## Question

Measure the effect of a flat-buffer/OpenMP optimization of the epsilon-active
HPWL direction and test an opt-in adaptive radius controller. Exact HPWL is
used for the objective and all final checks; no smooth HPWL surrogate or
bundle is enabled.

## Controlled variables

- center-Gaussian initialization, seed 1000, sigma ratio 0.001;
- Adam global optimizer and feasible refinement;
- adaptec1: 512 x 512, 1400 iterations, radius 125 and refinement radius 50;
- adaptec2: 1024 x 1024, 1600 iterations, radius 165 and refinement radius 65;
- no legal-aware checkpoint during GP, so the timing isolates the normal GP
  path; both configurations use the same final legalizer;
- adaptive mode updates a multiplicative radius scale every 25 iterations,
  clamped to [0.40, 1.20]. It widens the active band when overflow is above
  the target band and narrows it when feasible HPWL is rising.

## Variants

1. `fixed-r`: legacy fixed epsilon path.
2. `adaptive-r`: `--adaptive-active-set` with the controller above.

The default command-line behavior is unchanged: adaptive mode is disabled
unless the new flag is supplied.
