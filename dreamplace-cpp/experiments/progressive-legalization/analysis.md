# Progressive legalization analysis

## Question

Can a center-initialized exact, nonsmooth HPWL flow first spread density,
then introduce legality continuously, and thereby reduce the final discrete
legalization shock without using a smooth wirelength surrogate?

## Best observed result

The strongest run is `output/progressive_a1_512_1000_refine`:

| metric | value |
|---|---:|
| selected global HPWL | 83.254M |
| selected global overflow | 7.744% |
| Greedy HPWL | 89.971M |
| Abacus HPWL | 90.181M |
| final detailed HPWL | 88.895M |
| external legal2 Type 0/1/2/3 | 0 / 0 / 0 / 0 |

The requested 75-80M global range was not reached. Relative to the previous
center-initialized 512x512 exact baseline (85.401M/6.958%, 91.030M legal), the
new run improves global HPWL by 1.747M and final legal HPWL by 2.135M.

## Mechanism

The accepted controller is hybrid:

- above 11% overflow: DREAMPlace/RePlAce multiplicative lambda update;
- below 11% overflow: PI-like logarithmic band-dual update for 7%-8%;
- first entry into 8%: reset Adam moments, switch to AMSGrad, reduce the step
  by 20x;
- stage 3: smoothstep ramp of continuous fixed-macro, row, and segment forces;
- every trial: exact-HPWL/overflow/obstacle merit filter with backtracking.

The fixed-macro term is the squared distance to the complement of a fixed
rectangle expanded by the movable cell half-size. It always chooses the
nearest escape axis. No row, segment, site, or order is assigned until the
Greedy/Abacus phase.

## Evidence and rejected alternatives

| configuration | global HPWL | overflow | final legal HPWL | conclusion |
|---|---:|---:|---:|---|
| old exact 512, 800 | 85.401M | 6.958% | 91.030M | control |
| pure BandDual, 800 | 101.277M | 7.266% | 106.415M | density pressure too early |
| hybrid, 800 | 83.701M | 7.530% | 89.237M | large improvement |
| hybrid+refine, 1000 | 83.254M | 7.744% | 88.895M | best |
| sampling+bundle, 1000 | 83.577M | 7.741% | 89.186M | rejected |

Hard HPWL caps were also rejected. An 80M cap froze the trial filter around
81.9% overflow; a dynamic cap froze around 69.3%. At such high overlap, a
density-reducing move necessarily increases HPWL, so a global cap prevents
the layout from ever reaching a useful density state.

## Interpretation

The five-stage hypothesis is partially supported. Progressive obstacle and
segment forces substantially reduce the last legalization displacement, but
the short HPWL-only stage creates an almost fully overlapped center cluster.
The later density stage must reconstruct spatial order using sparse extrema
subgradients, which costs roughly 40M HPWL before the target band is reached.
The remaining gap is therefore primarily a global-topology problem, not a
legalizer-only problem.

The next credible experiment is a constrained density-tangent exact-HPWL
direction: project a sampled/bundle HPWL direction onto the local linearized
overflow descent half-space, and solve a small trust-region master. This can
use non-extreme pins through sampled active sets while keeping exact HPWL as
the acceptance oracle. Merely extending the current 1000-step trajectory is
not supported by the data because all trials after iteration 799 are rejected.
