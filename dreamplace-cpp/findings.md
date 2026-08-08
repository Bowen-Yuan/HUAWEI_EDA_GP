# Findings

## Current Understanding

- DREAMPlace does not use ePlace bound-to-bound initialization in the reported
  flow.  Movable cells start at the layout center with Gaussian standard
  deviation equal to 0.1% of each layout dimension.
- Fixed objects enter the initial electrostatic density map.  Small movable
  cells are stretched to at least sqrt(2) bin widths/heights while preserving
  area through a density ratio.
- Legalization is a separate discrete phase: macro legalization, Tetris-like
  Greedy legalization, then Abacus.  Detailed placement is required for the
  final Table II HPWL.

## Lessons And Constraints

- Comparing a pre-legal HPWL to paper Table II is invalid.
- A legalizer that ignores pin offsets can select the wrong local moves.
- A capacity-only density metric is insufficient evidence of geometric
  legality; the official `legal2.pl` result is authoritative.
- Overflow must be evaluated without fillers.  Fillers participate in the
  electrostatic objective, but DREAMPlace passes zero filler nodes to its
  overflow operator.
- BB curvature must compare gradients evaluated with the same density weight
  and wirelength gamma.  Mixing gradients from consecutive outer objectives
  stalls spreading.
- Greedy legalization must preserve the GP left-to-right order.  A global
  width-first ordering changed adaptec1 post-legal HPWL from 74.819M to
  123.992M at the same GP checkpoint.

## Validated Outcome

- adaptec1: 71.235743M at 6.9654% overflow before legalization; 74.819421M
  after legalization/refinement; `legal2.pl` Type 0/1/2/3 all zero.
- adaptec2: 80.889298M at 6.9302% overflow before legalization; 84.426599M
  after legalization/refinement; `legal2.pl` Type 0/1/2/3 all zero.
- Relative to paper Table II, the final legal gaps are +2.18% and +2.68%.

## Exploratory Smooth-To-Exact Continuation

- A smooth global placement is a strong warm start for exact HPWL: the exact
  bundle phase begins around 71.24M instead of the roughly 85M reached by the
  center-initialized exact path.
- On adaptec1, the current serious/null bundle lowered the warm-start GP HPWL
  by about 15K while staying at 7% overflow, but a matched enhanced-legalizer
  control was 1.16K better after legalization.  Lower pre-legal HPWL is not by
  itself evidence of a better final placement.
- Most later trials reduced exact HPWL but exceeded the hard 7% guard by only
  a small numerical amount and were classified as null steps.  A constrained
  filter with a narrow overflow hysteresis band is more promising than simply
  extending this trajectory.

## Grouped And Legal Bundle Study

- Four-group HPWL bundle histories plus nearby exact-subgradient sampling
  lowered the smooth warm-start GP result by 16.098K at a strict 7% boundary.
  However, the matched final legal HPWL was 9.103K worse.  Better continuous
  HPWL still does not imply a better legal topology.
- Bisection of a feasible-to-infeasible trial successfully preserves strict
  7% overflow, but once the center lies on the boundary the accepted fraction
  tends to zero.  A density-tangent master is needed for further pre-legal
  progress.
- The post-legal fixed-row, fixed-segment, fixed-order bundle is effective.
  Against a matched five-round no-bundle result of 74.520409M, it reached
  74.425850M, a 94.559K (0.127%) net improvement, with full internal legality.
- Default-off regression reproduced the previous 74.548626M result exactly.

## Progressive Legalization Findings

- A 100-step exact-HPWL-only phase is numerically effective but structurally
  destructive: HPWL falls to 45.714M while overflow remains 99.43%. The
  density stage must then rebuild almost the entire spatial topology.
- The key lambda result is a two-regime controller. DREAMPlace-style
  multiplicative updates spread the extreme-overflow state more effectively;
  a PI-like band-dual update is useful only after overflow drops below 11%.
  Pure band-dual control from the start overweights density and gives 101.277M
  at the same 512x512 resolution.
- Continuous fixed-macro and segment forces reduce the final discrete
  legalization jump to 5.641M: 83.254M GP becomes 88.895M legal. This is much
  smaller than the old 64x64 jump, but it does not solve the global-topology
  loss created during density spreading.
- The current best center-initialized, exact-only path improves the previous
  512x512 result by 1.747M before legalization and 2.135M after legalization.
  It still misses the requested 75-80M GP band by 3.254M at its upper edge.
- More iterations with the current filter do not help. After iteration 799,
  all trials are rejected while lambda changes, so the state is frozen.
  Further work needs density-tangent directions, a constrained master problem,
  or an active-set method rather than another scalar lambda schedule.
- Gradient sampling plus the current global serious/null bundle is not enough:
  it is 0.323M worse pre-legal and 0.291M worse legal than the Adam/AMSGrad
  control. Bundle remains useful after legalization, where row, segment, and
  order are fixed and the feasible polytope is much better defined.
