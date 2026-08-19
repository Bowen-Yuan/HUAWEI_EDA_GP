# State-triggered per-net epsilon experiment (2026-08-19)

## Question

Can a state-triggered transition from the frozen epsilon-active baseline to a
per-net, per-axis span cap improve final legal HPWL without increasing the
original iteration budget?

## Locked evaluation

The primary metric is final exact HPWL after the full legalizer. A retained
placement must pass the internal legality check and ISPD 2005 `legal2.pl` with
Type 0/1/2/3 errors equal to 0/0/0/0. Secondary metrics are selected GP HPWL,
overflow, selected iteration, transition iteration, objective iteration count,
and measured wall time. Global-placement curves are recorded as CSV; no
snapshots or animation are generated during timed runs.

The frozen 2026-08-18 smart adaptive epsilon results are the comparison
baseline. The grid, seed, initialization, fixed radii, base iteration budgets,
legal checkpoint interval, and detailed-placement settings remain unchanged.

## State trigger

Stage 2 may start after at least 250 refinement iterations when a trailing
50-iteration window satisfies all of the following:

- maximum overflow is at most 8.5%;
- overflow range is at most 0.4 percentage points;
- exact HPWL decreases by at least 0.03% across the window.

The configured start iteration remains a fallback deadline. If the state
trigger fires early, continuation is stretched to the same fixed phase end, so
the total iteration budget does not change.

## a1/a2 ablation

The following confirmatory variants are selected before inspecting results:

1. state-triggered continuation to a 15% directional net-span cap;
2. state-triggered continuation to a 10% directional net-span cap;
3. state-triggered continuation to 15%, followed by 50 low-step AMSGrad exact
   subgradient iterations with the existing exact HPWL/overflow filter and
   legal-aware rollback.

For adaptec1, the fixed fallback is iteration 1100; for adaptec2 it is 1300.
Variants 1 and 2 continue to the original 1400/1600 endpoint. Variant 3 ends
continuation 50 iterations earlier and uses those 50 iterations for the exact
subgradient trial, preserving the same total budget.

## Frozen six-design run

One variant is selected using only adaptec1/adaptec2 final legal HPWL, with
legality as a hard constraint and runtime as a tie-breaker. It is then frozen
and run once on adaptec1--4 and bigblue1--2. No parameters are changed after
the six-design run starts.

The success criterion is a lower geometric-mean legal HPWL ratio than the
frozen epsilon baseline, with no legality regression and no increase in the
base objective iteration budget.
