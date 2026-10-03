# Method and evaluation notes

## Objective protocol

`rho_target` is a density capacity parameter, not the overflow definition.
For a positive-capacity bin, the linear density penalty is
`max(occupancy - rho_target * available_area, 0)`.  The legacy protocol sums
only those bins.  The strict protocol additionally charges every movable area
inside zero-capacity bins.  The controller, archive and final selection all
use strict overflow; legacy is exported side-by-side.

## Controller

The solver starts with a state machine.  Catastrophic normalized HPWL invokes
HPWL rescue; high strict overflow invokes density rescue; when both occur, the
larger normalized violation chooses the temporary priority.  Rescue acceptance
is guarded: HPWL rescue has an overflow budget and density rescue has both
block and total HPWL damage budgets.  Only after hysteresis conditions are met
does the controller use `HPWL + lambda * D_strict`, where `lambda` is based on
EMA gradient-norm balancing and strict-overflow pressure.

The local density direction is full 2-D by default.  A single short-axis
direction is used only if its two overlap components have a configurable
dominance ratio.  Plateau escape moves cells through the nearest feasible
cardinal direction and accepts the batch only after an exact strict-oracle
check.

## Experiment interpretation

Short smoke runs validate the code path; they do not establish placement
quality.  In particular, a capacity-random adaptec3 start can reach low strict
overflow while remaining in HPWL rescue.  Any algorithmic conclusion should
compare matched wall-clock budgets, multiple seeds and strict-Pareto fronts.

## Scope

This folder owns adaptive search, artifacts and visualization.  `code/gpplacer`
owns Bookshelf parsing, placement projection and exact oracle kernels.  The
strict zero-capacity extension is intentionally in the shared evaluator so all
experiments can be checked by the same definition.
