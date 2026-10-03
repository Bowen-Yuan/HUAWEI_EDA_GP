# Architecture and Mathematical Contract

## Scope

The executable performs global placement only. Its output is a continuous
Bookshelf placement named `global.pl`. There is no row snapping, site
alignment, overlap removal, legalization, detailed placement, or legality
gate in the executable or core library.

The implementation is independent of `dreamplace-cpp`. That project was read
to understand its data conventions, epsilon-active HPWL direction, adaptive
epsilon controller, and DREAMPlace-style lambda update. No source module is
included, linked, or copied as a framework dependency.

## Module Boundaries

| Module | Responsibility |
|---|---|
| `bookshelf` | Parse raw ISPD 2005 Bookshelf files and write continuous GP placements |
| `hpwl` | Evaluate exact pin-offset HPWL and an epsilon-active direction |
| `density` | Evaluate exact rectangle/bin overlap, overflow, and active-piece directions |
| `lambda_controller` | Initialize and adapt the density multiplier |
| `optimizer` | Stateful optimizer interface and Adam/AMSGrad/AdaGrad/HeavyBall/SGD implementations |
| `placer` | Compose exact oracles, checkpoint feasible GP states, log metrics, and emit snapshots |

## Exact Objective

For each net `e`, pin coordinates include their Bookshelf offsets. The
wirelength term is

```text
HPWL_e = max_i x_i - min_i x_i + max_i y_i - min_i y_i.
```

For a grid bin `b`, occupancy is the sum of exact rectangle intersection areas
of movable instances and physical fixed macros. `terminal_NI` objects are
ports, so they do not consume physical bin capacity. With bin area `A_b`,

```text
rho_b = occupancy_b / A_b
D = sum_b 0.5 * A_b * max(rho_b - target_density, 0)^2.
```

Reported overflow is

```text
sum_b A_b * max(rho_b - target_density, 0) / total_movable_area.
```

Fixed macros therefore reduce the remaining capacity of every bin by their
exact intersection area. They are immutable, and only movable instances
receive density directions.

## Epsilon-Active Direction

Epsilon never changes `HPWL_e`, `rho_b`, `D`, or reported overflow. It only
selects and weights affine pieces whose value lies within epsilon of the
current active `min`/`max` piece.

For HPWL, pins within epsilon of a net extremum receive normalized triangular
weights raised to `active_power`. For rectangle/bin overlap, both affine
branches of each `min` and `max` endpoint are admitted when their values are
within epsilon. Density uses an independent `density_active_power` (linear by
default), because fourth-power net-extremum selectivity suppresses useful bin
boundary pieces. This gives a useful direction when a small instance lies in a
locally flat overlap region, while all acceptance, checkpointing, and metrics
remain exact.

## Parallelism

OpenMP parallelizes net evaluation, node gradient gathering, movable/bin
intersection work, reductions, and optimizer coordinates. The CLI accepts
1-40 threads and rejects larger values. A value of zero selects the smaller of
the runtime maximum and 40.

## Extensibility

Optimizers implement a single `compute_delta` interface. Adding a solver does
not require changing either oracle. Lambda policies are similarly isolated
from objective evaluation. This separation is essential because future
experiments may replace Adam or the scalar lambda policy without changing the
nonsmooth model.
