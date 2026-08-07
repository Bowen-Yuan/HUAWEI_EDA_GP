# H0 End-To-End Pipeline Protocol

Status: locked before execution

## Hypothesis

A faithful center initialization, pin-offset-aware WA objective, electrostatic
density, and DREAMPlace-style LG/DP pipeline should produce legal adaptec1 and
adaptec2 placements from raw Bookshelf input without any warm start.

## Sanity Gates

1. Log the exact `.pl` path and reject forbidden artifact names.
2. Verify the initial movable-cell centroid is at the layout center.
3. Unit-check WA and DCT gradients on synthetic cases.
4. Require finite metrics and decreasing overflow during GP.
5. Require internal legality after Greedy, Abacus, and every DP pass.
6. Require Type 0-3 counts of zero from `legal2.pl`.

## Metrics

- exact pin-offset HPWL before legalization
- electrostatic overflow
- exact HPWL after Greedy and Abacus
- exact HPWL after detailed placement
- mean/max legalization displacement
- wall-clock time

