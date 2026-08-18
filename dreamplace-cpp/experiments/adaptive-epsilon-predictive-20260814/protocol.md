# Predictive epsilon-active experiment

## Objective

Improve the exact, nonsmooth HPWL placement path without replacing the
objective by a smooth surrogate.  Epsilon only defines a trial direction;
all metrics, feasibility tests, checkpoints, and legal-placement comparisons
use exact HPWL.

## Controller

The opt-in predictive controller combines two scales.

1. A temporal controller updates the global epsilon scale from windowed
   overflow lag, overflow pressure, and exact-HPWL trend.
2. A per-net cap limits epsilon in each coordinate to a fraction of that
   net's bounding-box span.  This prevents a fixed absolute radius from
   selecting both sides of a short net and cancelling its direction.

The controller has density, boundary, and refinement phases.  Refinement
updates use a smaller log-step cap.  Legacy fixed and adaptive modes remain
unchanged unless `--adaptive-active-predictive` is passed.

## Screening plan

- Dataset: adaptec1, seed 1000, 512 x 512 bins, 800 iterations.
- Common optimizer: Adam, exact HPWL, electrostatic density, feasible-band
  refinement, no bundle, no legal-checkpoint sampling.
- Variants: temporal-only, span cap 0.50, span cap 0.35.
- Promote the best feasible candidate by final legal HPWL, then GP HPWL.

## Confirmation plan

- Run the promoted candidate for 1400 iterations on adaptec1 and 1600 on
  adaptec2.
- Compare against controlled fixed epsilon, legacy adaptive epsilon, and the
  first-generation smart controller.
- Freeze parameters before running the eight ISPD 2005 benchmarks.

