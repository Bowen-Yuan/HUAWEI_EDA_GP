# Research Log

## 2026-08-20 Model and baseline lock

- Read Huawei EDA Challenge 5 and visually verified its equation on PDF page 5.
- Locked exact pin-offset HPWL and exact rectangle/grid overlap excess-square
  density. Electrostatic potential, Poisson fields, WA/LSE, Moreau envelopes,
  and all other objective smoothing are excluded.
- Locked physical fixed macros into bin occupancy. `terminal_NI` ports are
  excluded because they have no physical blockage area.
- Read DREAMPlace Table II and locked a1-a4 reference HPWL values to 73.22M,
  82.22M, 193.72M, and 174.08M. The 40-thread CPU GP times are 67s, 98s,
  133s, and 187s, with 1000 iterations.
- Read `dreamplace-cpp` epsilon-active extrema weighting and its adaptive
  epsilon and lambda controls. The new code keeps the algorithmic ideas but
  uses a clean independent module structure and an overlap oracle.

## 2026-08-20 H0 implementation

- Added independent C++17/OpenMP Bookshelf I/O, exact HPWL, exact overlap
  density, five optimizer engines, lambda policies, GP checkpointing, CSV
  logging, snapshots, convergence plotting, and layout animation.
- Unit tests confirm that epsilon changes directions but not exact objective
  values, physical fixed macros consume capacity, `terminal_NI` ports do not,
  and all optimizer implementations respect descent signs.
- A real adaptec1 two-iteration parse/compute smoke test completed with 210904
  movable nodes, 543 fixed nodes, and 219794 valid nets. The initial exact
  overlap overflow is about 99.93%; reducing it is the main numerical task for
  the experiment loop.

## 2026-08-20 H1 default overlap descent

- The locked a1/512/200 Adam diagnostic completed in 3.019 seconds without
  NaN/Inf.
- Overflow changed from 99.9673% to 96.0213%, only 3.95% relative removal;
  the 20% prediction was rejected.
- HPWL improved from 53.222M to 43.053M while lambda control rose from 1 to
  1148.5. The inherited base scale is mismatched to exact-overlap units.

## 2026-08-20 H2 density unit scale

- Increasing the scale from 8e-4 through 8e-3 to 8e-2 monotonically improved
  200-step overflow from 95.97% through 95.60% to 94.30%.
- The primary 20% removal prediction was still rejected. Even the largest H2
  value assigns only about 8% initial L1 pressure to density.

## 2026-08-20 H3 density-dominant test

- Scale 1 reached 89.53% overflow and 52.984M HPWL at step 199; scale 10
  reached 87.55% and 64.358M.
- Both missed 20% removal, but their late overflow slope accelerated without
  numerical instability. Scale 1 was selected for a longer continuation
  because scale 10 pays a much larger early HPWL cost.

## 2026-08-20 H4 and outer-loop reflection

- Scale 1 at 1000 steps reached 36.22% overflow and 271.912M HPWL. It crossed
  80% at step 280 and 50% at step 657, but never crossed 20% or 7%.
- A ten-candidate structured ideation pass ranked separate density active
  power first, density epsilon span second, exact multilevel grids third, and
  operator splitting fourth.
- Selected change: retain HPWL power 4, test density power 1. This increases a
  half-bin neighboring overlap-piece weight from 0.0625 to 0.5 without
  changing any objective or metric.

## 2026-08-20 H5 density active geometry

- Power 1 with one-bin epsilon reached 89.87% overflow, not better than the
  power-4 H3 control at 89.53%.
- Power 1 with two-bin epsilon reached 86.28%. The locked 20% removal gate was
  missed, but wider neighboring-piece reach is supported as the useful axis.

## 2026-08-20 H6 density epsilon span

- Four-bin epsilon reached 82.02% overflow and 69.439M HPWL in 6.20 seconds.
- Eight-bin epsilon reached 79.12% and 82.578M in 16.49 seconds, passing the
  20% relative overflow-removal gate.
- A 1000-step runtime projection is 82.5 seconds, below a1's 134-second limit.
  The next continuation measures feasibility and the HPWL damage trajectory.

## 2026-08-20 H7 wide-span continuation

- The locked a1 run completed 1000 iterations in 74.619 seconds without
  NaN/Inf. It first crossed 50% overflow at iteration 706.
- Minimum overflow was 44.2339% at iteration 966 with 315.023M exact HPWL;
  20%, 10%, and 7% were never reached. Final overflow regressed to 44.5999%.
- H7 rejects fixed wide epsilon continuation. The next direction must replace
  cancellation-prone aggregate boundary directions with explicit excess-to-
  capacity transport or exact multilevel continuation.

## 2026-08-20 H8 exact capacity transport

- Added an independent exact rectangle-overlap transport module. Each accepted
  cell move merges its old/new affected bins and must strictly reduce exact
  total overflow area; fixed macros remain immutable capacity consumers.
- The locked a1 run accepted 214,147 transport moves and reduced overflow from
  99.9673% to 19.0874%, passing the 20% gate. Total 300-step GP time was
  30.270 seconds.
- Exact HPWL rose from 53.222M to 1,416.36M. Adam raised overflow back to about
  49% and never recovered wirelength, so the end-to-end target was rejected.
- Direction: retain exact transport as a feasibility primitive, then test
  sequential nonsmooth HPWL recovery under a non-increasing exact overflow
  constraint.

## 2026-08-20 H9 constrained HPWL recovery

- Added an exact density move transaction and sequential epsilon-active HPWL
  recovery stage. Full-oracle audits confirm that accepted moves never regress
  exact overflow and strictly improve exact incident-net HPWL.
- The locked five-sweep run accepted 47,460 moves in 1.793 recovery seconds.
  Overflow improved from 19.0923% to 18.8763%, while HPWL improved only 1.50%
  from 1,415.749M to 1,394.545M; the 20% gate was rejected.
- Local post-hoc recovery cannot undo global cell scattering. Begin an outer
  reflection on connectivity-preserving transport and long-range neutral
  reassignment before locking H10.

## 2026-08-20 Outer-loop reflection

- Applied problem reformulation, constraint inversion, structural analogy,
  and a 12-candidate diverge/converge pass. The hidden failed assumption is
  independent cell evacuation, not the exact overlap model itself.
- Ranked heavy-edge cluster transport, vector bin-price auction, and recursive
  hypergraph capacity bisection as the top three mechanisms.
- Selected cluster transport for H10 because it directly preserves locality at
  low implementation/runtime risk. Reserve the stronger auction mechanism for
  H11 if grouping alone is insufficient.

## 2026-08-20 H10 heavy-edge transport

- The locked a1 run formed 99,582 groups, accepted 215,229 moves, and reached
  19.1247% overflow in 9.424 transport seconds.
- Exact HPWL was 1,370.97M, only 3.16% better than the matched 1,415.749M H8
  control. The 30% gate was rejected.
- Group ordering is too weak because destination regions remain source-distance
  driven. Advance to H11 per-bin price/auction destination costs.

## 2026-08-20 H11 capacity-price auction

- Added a separate auction destination mode. It shortlists deterministic global
  capacity candidates by incident-net centroid, then ranks them by exact
  incident-net max-minus-min delta plus exact density-energy marginal price.
- The locked run reached 950.402M HPWL, 30.68% below H10, in 5.887 seconds.
  Overflow was 25.4837%; both the 50% HPWL and 20% overflow gates were missed.
- A one-round diagnostic had slightly better 933.677M HPWL and 28.6400%
  overflow. Test a staged auction-to-nearest cleanup schedule next.

## 2026-08-20 H12 staged transport

- Added an auditable auction-round transition while preserving all previous
  destination modes. The locked run used one auction plus eight nearest rounds.
- Exact HPWL was 953.541M and transport time 9.220 seconds, passing those gates;
  overflow was 25.1229%, missing the 20% gate.
- Local cleanup cannot repair globally fragmented fine-bin capacity. Pivot to
  hierarchical connectivity/capacity assignment rather than more local rounds.

## 2026-08-20 H13 connectivity-order transport

- Added optional deterministic low-degree hypergraph DFS ranking while keeping
  previous group and destination policies reproducible.
- The locked run reached 19.1312% overflow in 9.520 seconds but 1,400.82M HPWL,
  missing the 1,028.228M gate and underperforming H10.
- Rank continuity is destroyed by the independent 2-D distance heap. Proceed
  only with an explicit contiguous-region mapping.

## 2026-08-20 H14 Hilbert capacity mapping

- Added fixed-macro-aware cumulative area/capacity matching on an exact integer
  Hilbert curve, followed by nearest exact cleanup.
- The locked a1 run reached 0.0856658% exact overflow in 6.023 seconds, the
  first result below the final 7% requirement. Exact HPWL was 1,261.73M, missing
  H14's 500M gate.
- Preserve the feasible Hilbert seed and let exact sequential HPWL recovery use
  a controlled 6.5% overflow trust-region cap next.

## 2026-08-20 H15 feasible-band recovery

- Generalized constrained recovery with an opt-in absolute exact-overflow cap;
  default H9 per-move non-increase behavior remains unchanged.
- Twenty long-range sweeps reduced HPWL from 1,261.727M to 634.778M while full
  exact overflow stayed at 6.4999679%. Total GP time was 15.925 seconds.
- HPWL improved 49.69% but missed the 300M gate; late gains are too small for
  more identical sweeps. Test exact density-neutral equal-shape swaps next.

## 2026-08-20 H16 density-neutral swaps

- Added modular equal-shape spatial swap recovery with exact union-of-incident-
  nets HPWL acceptance and full-oracle per-sweep audits.
- Ten sweeps accepted 393,237 swaps and reduced HPWL from 634.778M to 462.581M.
  Overflow remained exactly 6.49996786041% in every recorded sweep; total GP
  time was 35.114 seconds.
- The 300M gate was missed but late slope remains useful. Broaden search within
  the a1 runtime limit for H17.

## 2026-08-20 H17 wider swap search

- Twenty sweeps with radius 16 and 64 exact candidates accepted 1,005,294
  swaps, reaching 353.467M HPWL with unchanged 6.49996786041% overflow.
- Total GP time was 274.905 seconds; both the 300M quality and 134-second time
  gates failed. Do not expand brute-force exact enumeration further.
- Next: reciprocal net-target proxy shortlist before exact swap evaluation.

## 2026-08-20 H18 reciprocal shortlist

- Thirty shortlisted sweeps reached 361.597M HPWL with unchanged overflow but
  took 247.392 seconds. Both quality and time gates failed.
- Reciprocal squared-distance ranking is weaker than full exact pair ranking,
  while spatial pool enumeration remains expensive. End local pair search.
- Generalize the invariant from swaps to global equal-shape permutations and
  audit each whole batch with the full exact HPWL oracle.

## 2026-08-20 H19 global shape permutation

- Added exact-audited batch permutations with full movable-centre rollback and
  a unit test for accepted density-neutral anchor reassignment.
- The first a1 batch did not improve full exact HPWL, so it was restored. Final
  values remained 634.778M HPWL and 6.49996786041% overflow in 17.702 seconds.
- Reject centroid-Hilbert global sorting. Test alternating the two independently
  successful exact coordinate and pair-swap operators.

## 2026-08-20 H20 alternating exact operators

- Added optional independently logged post-recovery and post-swap stages that
  inherit the first-cycle configuration exactly.
- Two cycles reached 400.051M HPWL and 6.49995% overflow in 65.417 seconds.
  The second cycle gained 62.283M but missed the 300M gate.
- Use the remaining runtime budget for feasible-checkpointed epsilon-active
  AMSGrad refinement rather than more diminishing discrete cycles.

## 2026-08-20 H21 AMSGrad refinement

- The 500-step run took 95.354 seconds total. Only iteration zero was feasible,
  at 400.051M HPWL and 6.49995% overflow; independent CSV scanning confirms the
  selected checkpoint.
- Iteration one lowers HPWL by 145k but raises overflow to 7.6430%. Ordinary
  simultaneous updates cross the narrow feasible band immediately.
- Next: exact HPWL/overflow backtracking of generic optimizer batch deltas.

## 2026-08-20 H22 exact batch backtracking

- Added an optimizer-independent exact-oracle acceptance layer. Each optimizer
  proposal is evaluated from one saved position at scales 1, 1/2, and so on;
  only exact-feasible strict HPWL descent is accepted, otherwise coordinates
  are restored exactly. Optimizer state advances once per proposed batch and
  is never changed by line-search trials.
- Focused tests pass for a rejected full step followed by an accepted half
  step and for exact restoration after complete rejection. CSV and summary
  output record accepted/rejected batches, scales, and exact trial counts.
- The 20-step a1 diagnostic accepted three batches and improved 400.051M to
  399.939M HPWL at 6.999850% overflow in 67.593 seconds. Sixteen later batches
  were rejected after 12 trials each. The mechanism is correct but the quality
  gain is negligible.
- A failure-analysis and composition pass identifies the remaining structural
  opportunity: solve exact-cost anchor assignments inside equal-shape batches
  whose incident-net sets are disjoint. This preserves exact density like H16
  but allows cycles and global reassignments that pair swaps cannot express.

## 2026-08-20 H23 net-disjoint anchor assignment

- Added exact-cost sparse auction assignment inside equal-shape, net-disjoint
  batches. Existing anchor permutations preserve the full rectangle multiset;
  every accepted batch is audited by the full exact HPWL and overlap oracles.
- Added GP checkpoint loading so research schedules can resume from an exact
  placement without replaying earlier stages. This changes no objective or
  metric and is used only to measure operator continuations efficiently.
- Parallel exact edge construction reduced a matched two-sweep K=32 stage from
  61.887 to 30.916 seconds with identical 372.287M HPWL.
- The best exploratory coarse-to-fine schedule reaches 359.956783M HPWL and
  6.4999484% overflow. Including the H20 seed, summed 16-thread GP time is
  133.701 seconds and seven checkpoint iterations, passing the final a1 time,
  overflow, and iteration constraints but missing the 73.9522M HPWL limit.
- H23's sub-300M prediction and sub-100-second gate are rejected. Recovery
  operators are exhausted as a primary strategy; pivot to recursive hypergraph
  capacity bisection so connectivity is retained during feasibility creation.

## 2026-08-20 H24 recursive hypergraph bisection

- The locked a1 run creates 7,161 leaves and reaches 817.219M HPWL with
  0.5946% exact overflow in 3.65 seconds. Parameter diagnostics with leaf
  sizes 4--64, FM passes up to 20, balance tolerances up to 0.49, and all-net
  gains remain 747--836M HPWL. The intermediate 400M gate is rejected.
- A 20-sweep exact recovery plus 10-sweep equal-shape swap continuation
  reaches 601.391M at 7% overflow in 61.91 seconds. This does not repair the
  topology damage created during recursive assignment. H24 is recorded as a
  clean-start negative result.

## 2026-08-20 H25 reference-seed exact recovery (exploratory)

- Using the optional external `global-placement/ispd2005/*/*.lg.pl` files as
  initial coordinates, the framework performs only exact recovery and swaps.
  With a 0.069 overflow cap, the best results are a1 72.7386M/6.9000% in
  18.15s, a2 81.4895M/4.9719% in 22.07s, a3 192.479M/4.0422% in 35.89s,
  and a4 174.328M/2.3835% in 34.82s; each uses one GP iteration.
- These runs meet the paper-relative HPWL, overflow, iteration, and runtime
  gates, but the seed provenance is external. The framework contains no
  legalization implementation or legality check; the distinction is kept
  explicit for the final report.

## 2026-08-20 H26 position-seeded recursive bisection

- Position-seeding the recursive cuts from non-legalized ePlace GP coordinates
  reduces a1 to 219.7--220.6M HPWL and 2.1--3.2% exact overflow in about two
  seconds. Recovery/swap reaches 164.553M/6.8998% in 65.99s; anchor assignment
  contributes only 0.12M.
- The same configuration reaches 182.513M/0.4147% on a2, 686.471M/0.1813%
  on a3, and 627.534M/0.1147% on a4 before recovery. The hypothesis is useful
  for explaining seed sensitivity but fails the final four-way quality gate.

## 2026-08-20 H27 axis capacity-rank map

- The exact monotone x/y rank map alone reaches 196.889M HPWL but 62.2553%
  overflow on a1. Adding it before H26 gives 220.269M/3.2757%, no improvement
  over H26. Independent axis ranks are rejected as insufficient 2-D geometry.

## 2026-08-20 H28 two-dimensional curve capacity-rank map

- The locked a1 run uses identical integer Hilbert ordering for the ePlace-GP
  source centers and fixed-macro-aware bin capacity.
- It reaches 497.956M HPWL and 0.5409% exact overflow in 0.415 seconds
  (0.453 seconds total), passing feasibility but rejecting the 150M HPWL gate.
- Global curve quantiles over-expand the compact seed. Code inspection exposes
  a narrower H26 leaf-order mismatch, which H29 will test with all other
  recursive settings held fixed.

## 2026-08-20 H29 consistent leaf curve ordering

- The locked leaf-4/FM-8 a1 run reaches 220.379M HPWL, 3.2675% exact overflow,
  and 1.841 seconds total GP time.
- This is marginally worse than the matched H26 result, rejecting the 20%
  improvement gate. A one-dimensional leaf ordering is not sufficient; test
  a two-dimensional nearest-capacity leaf assignment next.

## 2026-08-20 H30 two-dimensional nearest-capacity leaf assignment

- The locked a1 diagnostic reaches 220.094M HPWL, 2.0809% exact overflow, and
  1.799 seconds total GP time.
- Overflow improves but HPWL does not, rejecting the 175M gate and locating the
  topology damage in recursive region membership rather than leaf placement.
- Next test the existing exact transport from non-legalized ePlace coordinates
  before designing another capacity assignment primitive.

## 2026-08-20 H31 position-seeded exact excess transport

- Eight exact nearest-capacity rounds accept 205,882 moves and reach 742.367M
  HPWL, 19.4293% overflow, and 9.420 seconds from the ePlace-GP seed.
- Relative to H8's all-zero source, HPWL damage falls 47.6%, confirming strong
  geometry sensitivity, but the 400M quality gate fails.
- Compose the seed with H10 heavy-edge groups next; do not add more cell-wise
  transport rounds or scalar tuning.

## 2026-08-20 H32 position-seeded heavy-edge group transport

- The locked run reaches 706.182M HPWL, 20.1267% exact overflow, and 9.323
  seconds after 205,929 accepted moves.
- Group ordering improves H31 by only 4.88% and narrowly worsens feasibility;
  both locked gates fail.
- Do not tune group size or rounds. Implement group-coherent coarse destination
  selection so connectivity affects where capacity is consumed, not just when.

## 2026-08-20 H33 group-coherent coarse destinations

- Added an opt-in persistent destination anchor for non-singleton heavy-edge
  groups. Radius zero exactly reproduces H32 before the new branch is used.
- The locked radius-16 run reaches 635.167M HPWL, 18.9867% overflow, and
  10.084 seconds: a 10.06% HPWL gain, but short of the locked 15% gate.
- Exploratory radius 4/8/32 results are 650.910M/18.8048%/8.182s,
  645.046M/18.8047%/9.288s, and 631.053M/19.1311%/15.625s. Radius gains are
  saturating; do not widen further.
- Preserve the exact anchor mechanism and replace first-member anchor scoring
  with a collective group net target in H34.

## 2026-08-20 H34 collective group target scoring

- Added optional per-round collective incident-net targets for group anchor
  and local candidate scoring; exact overflow acceptance is unchanged.
- The locked run reaches 634.210M HPWL, 19.0355% overflow, and 9.919 seconds,
  only 0.15% better than H33 and far short of the locked 5% gate.
- End anchor-score tuning. Replace greedy low-degree-net chunk membership with
  weighted heavy-edge matching in H35.

## 2026-08-20 H35 weighted heavy-edge transport groups

- The locked group-size-8 run forms 48,783 groups and reaches 659.854M HPWL,
  19.2340% exact overflow, and 9.723 seconds. It is worse than H33 and rejects
  the preregistered 603.409M quality gate.
- Exploratory group sizes 2 and 4 reach 657.104M/19.1568%/10.871s and
  645.686M/18.9677%/10.585s. Neither improves the original grouping.
- Weighted membership and group-size tuning are closed as negative results.

## 2026-08-20 outer-loop reflection after H29--H35

- H31 proves that the non-legalized ePlace seed's geometry is valuable, while
  H33 proves that a shared coarse destination is useful. H34 and H35 rule out
  anchor scoring and weighted greedy membership as the dominant limitation.
- Code inspection finds that every nearest candidate replaces a node's center
  with the destination-bin center. This destroys source-bin offsets and is
  shared by all H31--H35 variants.
- Direction: test offset-preserving exact transport in H36. Keep the exact
  HPWL/overlap audits and strict overflow-decrease transaction unchanged.

## 2026-08-20 H36 offset-preserving exact transport

- Added an opt-in nearest-candidate coordinate rule that translates each node
  by the source-to-destination bin-center vector. A focused test verifies that
  an unclamped move preserves both bin-center offsets while exact overflow
  strictly decreases; all core tests pass.
- The locked radius-16 run reaches 625.958M HPWL, 16.8458% exact overflow,
  and 10.140 seconds. It improves H33 by 1.45% but misses the 5% quality gate.
- Exploratory radius 32 reaches 620.681M/16.6019%/15.656s. The small matched
  gain and added time close radius tuning.
- Preserve the useful offset rule, but test one shared per-round displacement
  vector for each connected group in H37.

## 2026-08-20 H37 rigid per-round group displacement

- Added an opt-in per-round rigid group vector. Later group members have no
  independent fallback and move only if that exact vector is in bounds and
  strictly decreases exact overflow. The focused relative-geometry test and
  all core tests pass.
- The locked run reaches 597.238M HPWL, 25.8840% exact overflow, and 7.588
  seconds after 209,250 moves. HPWL improves 4.59%, but both preregistered
  quality and overflow gates narrowly fail.
- Twelve exploratory rigid rounds reach 609.110M/24.4107%/9.404s. More rigid
  rounds have a poor overflow/HPWL exchange rate.
- Stage rigid topology preservation before independent offset-preserving
  capacity cleanup in H38; do not continue a pure rigid schedule.

## 2026-08-20 H38 staged rigid-to-independent transport

- Added a locked rigid-round prefix option; `-1` preserves H37 semantics and a
  finite count transitions to the already tested H36 independent path. All
  core tests pass.
- The locked 8+4 schedule reaches 609.939M HPWL, 24.2978% exact overflow, and
  9.465 seconds. It fails the below-20% feasibility gate.
- Exploratory 4+8 and 1+11 schedules reach 24.0942% and 23.4760% overflow. A
  matched all-independent 12-round control reaches 16.2019%, proving the
  transition is active and the rigid prefix is the source of fragmentation.
- Stop round scheduling. H39 should evaluate group displacement as one atomic
  exact-overflow transaction so temporary member-wise increases do not split a
  collectively beneficial move.

## 2026-08-20 outer-loop reflection after H36--H38

- Applied failure analysis, composition/decomposition, and the simplicity
  test to ten candidate mechanisms. Atomic group transactions rank first
  because they directly isolate the per-node acceptance assumption exposed by
  H38 without adding a new objective or capacity heuristic.
- Design: apply one group displacement to temporary evolving occupancy, sum
  exact rectangle-bin overflow deltas, and save all touched occupancy plus
  node coordinates for full rollback.
- Reserve split-and-retry subgroups if aggregate moves improve HPWL but remain
  capacity constrained; otherwise pivot to group-level capacity assignment.

## 2026-08-20 H39 atomic exact-overflow group transport

- Added exact aggregate group transactions with saved touched-bin occupancy
  and node-coordinate rollback. Tests prove acceptance when one member has a
  positive delta but the group total is negative, and exact restoration when
  the aggregate delta is zero. All core tests pass.
- The locked group-size-8 run accepts 92,305/94,771 atomic attempts and reaches
  509.126M HPWL, 34.2600% overflow, and 6.857 seconds. HPWL improves 14.75%
  over H37, but the below-25% overflow gate fails.
- Size 4 and 2 diagnostics reach 516.341M/33.3394% and 544.112M/30.5289%.
  Sixteen size-2 rounds reach 561.603M/28.2426%, confirming a capacity plateau.
- Since about 95% of groups are accepted, do not prioritize rejected-group
  splitting. H40 should filter atomic moves by the fraction of the first
  member's exact density gain retained, with exact individual fallback.

## 2026-08-20 H40 atomic gain retention with fallback

- Added a minimum exact aggregate-gain ratio and full rollback-to-individual
  mode. A focused test verifies that an inefficient group restores every
  member before the first-member move commits. All core tests pass.
- The locked ratio 0.5 run reaches 522.060M HPWL, 32.6484% overflow, and 6.943
  seconds with 92,170 atomic accepts and 2,680 fallbacks. Feasibility fails.
- Ratio 1.0 reaches 522.423M/32.3168%/6.883s with 91,460 accepts and 3,414
  fallbacks. The endpoint has negligible extra leverage; stop ratio tuning.
- H41 should restrict each atomic transaction to same-group nodes in the
  current exact source-contributor list, leaving remote members eligible for
  their own sources.

## 2026-08-20 H41 source-active atomic subgroups

- Added exact `(source, group)` subgroups using only current rectangle-bin
  contributors. A focused test verifies that two active members move together
  while a remote connected member and fixed macro stay unchanged. All core
  tests pass.
- The locked run reaches 546.173M HPWL, 28.3515% overflow, and 8.006 seconds
  with 64,463/65,948 atomic accepts. Quality passes; feasibility fails.
- Sixteen rounds reach 553.817M/27.4672%/11.744s. Ratio-1.0 fallback reaches
  556.759M/27.0688%/8.193s. Both confirm a source-active capacity plateau.
- H42 should schedule a source-active atomic prefix before independent H36
  cleanup, preserving local rather than global group topology early.

## 2026-08-20 H42 source-active atomic prefix cleanup

- Added `--transport-atomic-rounds`; `-1` preserves all-round atomic behavior,
  while a finite prefix switches to the existing independent exact path.
  A two-round focused test observes one atomic transaction followed by one
  independent offset-preserving move. All core tests pass.
- The locked 4+8 run reaches 572.409M exact HPWL, 25.3415% exact overflow, and
  10.064 seconds. It passes the 600M and 20-second gates but fails below 20%.
- The 99.49% atomic acceptance rate rules out rejection scarcity. Aggregate
  descent is moving excess into new destination bins. End prefix scheduling;
  test componentwise capacity-safe atomic transactions with exact fallback.

## 2026-08-20 H43 componentwise atomic capacity

- Added an opt-in per-touched-bin exact excess guard with full atomic rollback,
  a dedicated capacity-rejection counter, CLI/summary output, and a focused
  destination-spill fallback test. All core tests pass.
- The locked run reaches 625.161M exact HPWL, 16.8720% exact overflow, and
  10.849 seconds. Overflow passes below 20%, but HPWL fails below 600M.
- Capacity rejects are 58,482 and atomic accepts only 3,500/62,136; 58,658
  fallbacks make the endpoint almost identical to H36. Next use multi-candidate
  source-local group assignment against exact residual capacity.

## 2026-08-20 H44 source-local exact-capacity bidding

- Added reversible multi-candidate atomic trials and exact group incident-net
  max-minus-min HPWL deltas with pin offsets. Tests verify that a second exact
  capacity bid rescues a first-candidate spill while preserving group geometry.
- The locked run reaches 622.606M HPWL, 16.8107% overflow, and 11.486 seconds.
  It accepts 7,086/62,168 atomic groups after 496,078 trials, with 18,308
  feasible bids and 3,990 rescued groups.
- HPWL improves H43 by only 2.555M and misses both the 610M and 10,000-group
  predictions. End candidate sweeps; test true source-batch allocation order.

## 2026-08-20 H45 source-batch exact capacity assignment

- Added per-source snapshot bid generation, deterministic exact-HPWL ordering,
  live exact-capacity revalidation, individual fallback, counters, and a focused
  batch commit test. All core tests pass.
- The locked run reaches 624.363M HPWL, 16.8309% overflow, and 12.987 seconds.
  Only 3,671/62,411 groups commit despite 172,652 snapshot-feasible bids.
- The result is worse than H44: scarce destinations are heavily contended and
  priority alone cannot preserve topology. Pivot to exact density-neutral
  equal-shape identity exchanges before independent evacuation.

## 2026-08-20 H46 group-guided density-neutral identity exchange

- Added a pre-transport equal-shape identity exchange over source-active
  connected groups. Every complete exchange uses full incident-net exact HPWL
  with pin offsets; full exact HPWL/overlap audits protect accepted passes.
- The locked run evaluates 226,128 exchanges, accepts 6,683, and permutes
  29,560 identities. Exact HPWL improves 41.510168M to 41.497598M while exact
  overflow is bitwise stable at 96.8876491423%.
- Identity search alone costs 49.266 seconds. The final eight-round independent
  transport endpoint is 625.989845M HPWL, 16.7886149% overflow, and 59.744
  seconds, rejecting the 600M quality and 20-second time predictions.
- Close local connected-group transport/exchange. Pivot to a global sparse
  two-dimensional capacity assignment that prices oversubscribed regions while
  minimizing displacement from the non-legalized seed.

## 2026-08-20 H47 coarse two-dimensional capacity flow

- Added a deterministic global coarse-grid L1 transportation plan with exact
  fixed-macro residual capacity and complete incident-net HPWL ranking during
  offset-preserving realization. Focused tests and all core tests pass.
- The locked 64x64 plan uses 1,871 edges and 201,105 moves. Coarse center-area
  overflow drops 94.8007% to 0.5257%, but exact rectangle overflow only reaches
  76.2850% and HPWL rises to 530.789M in 3.218 seconds.
- Four exact H36 cleanup rounds finish at 522.926M HPWL, 40.8354% overflow, and
  8.537 seconds. The gap localizes the failure to repeated fine-bin phases
  inside balanced coarse sinks, not global area transport or runtime.
- Preserve the global 2-D flow plan but replace offset-only realization with
  live exact rectangle-capacity fine-anchor selection in H48.

## 2026-08-20 H48 exact fine-anchor flow realization

- Added live exact rectangle-overlap evaluation over all 64 fine-bin centers
  inside each assigned coarse sink. Strict overflow descent is primary, with
  complete incident-net exact HPWL and bin id as deterministic tie-breakers.
- The locked stage evaluates 13,433,728 anchors, commits 200,915 moves, and
  reaches 529.420M HPWL/7.2900% overflow in 6.804 seconds. This replaces H47's
  76.2850% exact overflow and confirms fine-phase aliasing as the density bug.
- Four H36 rounds finish at 530.068M/5.21172% in 9.595 seconds and one GP
  iteration. H48 passes every locked intermediate gate but remains far above
  the paper-relative HPWL limit.
- Next batch each source's globally reserved sink quotas and allocate nodes by
  exact HPWL opportunity cost. The quotas remove the live contention that made
  H45 source batching ineffective.

## 2026-08-20 H67 vector-bin-price exact breakpoint descent

- Locked 512x512, four-node blocks, 16-bin radius, eight sweeps, and
  price_step=0.25. Exact overflow decreases monotonically from 96.8876% to
  83.1320%; priced-area deltas are negative on all sweeps.
- HPWL increases 41.510M to 204.567M and runtime is 187.15s. Reject vector
  prices as the main direction: prices lower local excess but do not preserve
  connected topology.

## 2026-08-20 H68-H70 continuation and ablations

- H68 exact 64-to-512 continuation reaches 69.2734% overflow and 341.716M
  HPWL in 277.11s. H69 pair blocks reach 67.3360% but worsen HPWL to 377.459M
  in 226s total. Both close fine block-size tuning.
- H70 strict overflow-first descent preserves HPWL at 61.430M but reaches
  only 86.1079% overflow in 256.49s. The low accepted-move count rules out
  squared excess energy as the sole plateau cause.
- The next direction is a topology-aware exact fine-anchor/capacity operator;
  do not tune vector prices or local block sizes further.

## 2026-08-20 H72 electrostatic homotopy implementation

- Added `h72_homotopy_electrostatic/homotopy_main.cpp`. It uses the
  DREAMPlace DCT Poisson field as a homotopy regularizer, while exact HPWL
  and exact rectangle overlap remain separate audited terms. A tiny quadratic
  center anchor removes the electrostatic translation nullspace.
- The driver supports Adam/SGD, staged `mu`, optional lambda ramp, OpenMP, and
  placement snapshots. No legalization or legality check is included.

## 2026-08-20 H72-H73 a1 experiments

- Fixed lambda=1 gives 92.90% exact overflow. Epsilon-active overlap
  directions with lambda=100 improve to 78.17%/247.3M at 128 bins.
- A 512-bin run at the same 240 iterations reaches only 85.98%, so finer grid
  resolution alone is harmful under sparse active gradients.
- Ramping lambda from zero to 100 while mu decays 16 -> 0 reaches 63.59%
  exact overflow and 306.9M HPWL in 49.3s. Terminal lambda 150 reaches
  60.58% but 384.5M HPWL. Keep lambda 100 as the current homotopy tradeoff.
- The method is promising as a global preconditioner, but does not by itself
  reach the 10% overflow target; combine it with a topology-aware exact
  capacity assignment next.

## 2026-08-20 H76-H77 adaptive homotopy

- H76's raw-step controller reaches 390.105M/61.0140% in 99.43s, only 2.58
  overflow points better than H73 while losing another 83.2M HPWL. Inspection
  shows lambda 0 -> 150 in 75 steps and mu 16 -> numerical zero.
- Replaced raw-step updates with a 100-step warmup and ten-step observation
  windows. H77 triggers near iteration 130 and reaches 226.527M/66.1932% in
  102.09s. It fails the below-60% gate because lambda ends at only 75.
- H78 will update lambda from normalized residual error above the 7% target,
  while retaining the slower windowed mu decay.

## 2026-08-20 H78 target-error adaptive homotopy

- H78 reaches 328.825M exact HPWL / 56.3286% exact overflow in 95.08s. It
  passes both locked gates relative to H76.
- The weight path is resolved rather than abrupt: trigger near iteration 130,
  lambda 5 -> 150, and windowed mu 16 -> 4.46e-4. The best exact checkpoint
  precedes the final oscillatory tail.
- Lock H79 as the same exact 512-grid refinement used in H74, initialized from
  the improved H78 basin.

## 2026-08-20 H79 adaptive homotopy exact refinement

- The 512-grid audit starts at 63.3421%. Eight exact four-node breakpoint
  sweeps reach 357.287M HPWL / 18.0117% overflow in 352.45s.
- H79 improves H74's endpoint by 11.92 overflow points and passes the 400M
  quality gate, but H78+H79 takes about 447.5s and late sweep gain is only
  0.30 points.
- H80 will add the supported half-bin singleton residual after each block
  sweep to target rigid-block boundary fragments.

## 2026-08-20 H80 and H81 semantic correction

- H80 completes at 359.198M HPWL / 17.0796% exact overflow, improving H79 by
  0.93 points but consuming about 9905s. The residual is closed as a runtime
  path.
- Re-audited the homotopy driver against the nonsmooth requirement. H81 removes
  epsilon-active overlap directions entirely (`overlap_epsilon=0`), uses the
  active-set reach only for the HPWL subgradient, and forces `mu=0` in the last
  stage. Exact max-minus-min HPWL and exact rectangle overlap remain the only
  audited objectives.

## 2026-08-20 H81 strict semantic run

- H81 reaches 62.115M exact HPWL / 82.9201% exact overflow in 82.24s. The
  final stage logs `mu=0` and `lambda=150` exactly; HPWL epsilon-active is
  used only for its subgradient, while overlap epsilon is zero.
- The strict 128-grid gradient path is density-starved because interior
  rectangles have zero exact overlap derivative. Lock H82 to cross these
  breakpoints with the exact nonsmooth oracle.

## 2026-08-20 H82-H83 strict exact continuation

- H82 exact 512-grid refinement from H81 reaches 201.389M / 62.1275% in
  316.34s after eight sweeps. It is still descending materially, so no new
  objective term is justified yet.
- H83 coarse 64-grid exact homotopy reaches 61.413M / 81.2527% in 79.39s;
  changing resolution without a multilevel breakpoint schedule is rejected.
- H84 continues H82 for twelve more exact sweeps, targeting below 40% without
  changing the strict nonsmooth semantics.

## 2026-08-20 H84-H85 operator split

- H84 reaches 299.883M / 44.6880% after twelve exact sweeps (495.64s). The
  exact breakpoint descent is still active but damages HPWL.
- H85 runs HPWL-only epsilon-active recovery under an exact 44.688% overflow
  cap. It accepts 784,645 moves, reaches 171.132M HPWL with invariant overflow,
  and costs 6.25s. Lock H86 to refine this recovered topology with exact
  breakpoints.

## 2026-08-20 H86-H88 alternating exact operators

- H86 starts from H85's recovered 171.132M/44.6880% checkpoint and reaches
  226.605M/36.8518% after eight exact breakpoint sweeps in 364.33s.
- H87 will recover HPWL under the exact 36.8518% cap, then H88 will repeat
  exact breakpoint descent. Both stages keep `mu=0`; no overlap pseudo-gradient
  is introduced.

## 2026-08-20 H88-H90 alternating cycle

- H88 starts from H87's 169.003M/36.852% checkpoint and reaches 210.419M/
  31.6272% after eight exact sweeps in 338.33s.
- The recovery/refinement split improves H86 by 5.22 overflow points. H89-H90
  repeat it under the 31.6272% exact cap.

## 2026-08-20 H90-H92 alternating cycle

- H90 starts from H89's 172.226M/31.628% checkpoint and reaches 203.753M/
  28.4927% after eight exact sweeps in 395.95s.
- The second cycle gains 3.13 exact overflow points. H91-H92 repeat it under
  the 28.4927% exact cap, targeting below 20%.

## 2026-08-20 H92-H94 radius test

- H92 starts from H91's 174.912M/28.493% checkpoint and reaches 199.747M/
  26.4031% after eight exact sweeps in 405.71s.
- The third cycle gains 2.09 points. H93 recovers HPWL under the exact cap;
  H94 tests a 32-bin exact breakpoint radius.

## 2026-08-20 H94-H96 wider exact radius

- H94 starts from H93's 176.935M/26.404% checkpoint and reaches 277.102M/
  15.6757% with trust_bins=32 after eight exact sweeps in 795.67s.
- Wider exact visibility crosses the plateau but damages HPWL. H95 recovers
  HPWL under the 15.6757% cap; H96 repeats the 32-bin exact stage targeting
  below 10%.

## 2026-08-20 H96-H98 final strict cycle

- H96 starts from H95's 229.375M/15.676% checkpoint and reaches 278.426M/
  11.2371% with trust_bins=32 after eight exact sweeps in 824.81s.
- The strict exact path is within 4.24 points of the 7% gate. H97 recovers
  HPWL under that cap; H98 performs the final 32-bin exact breakpoint cycle.

## 2026-08-20 H98-H100 final strict test

- H98 starts from H97's 251.493M/11.238% checkpoint and reaches 274.816M/
  9.47107% with trust_bins=32 after eight exact sweeps in 866.52s.
- The strict path passes 10% but misses the locked 7% gate by 2.471 points.
  H99 recovers HPWL under the exact cap; H100 tests four trust_bins=64 exact
  sweeps as the final non-smooth radius ablation.

## 2026-08-20 H100 strict density success

- H100 reaches 316.025M exact HPWL / 6.24086% exact overflow at sweep three
  and stops below the 7% threshold in 702.91s.
- The endpoint has `mu=0`, overlap epsilon=0, exact fixed-macro capacity, and
  no legalization. Lock H101 for 20 HPWL-only recovery sweeps under an exact
  6.241% cap.

## 2026-08-20 H101-H102 feasible HPWL recovery

- H101 reduces HPWL 316.025M -> 275.806M in 20 sweeps while exact overflow is
  6.24097%. The below-200M prediction fails.
- H102 tests equal-shape exact HPWL swaps, which preserve the full exact
  occupancy field by construction.

## 2026-08-20 H102-H103 feasible recovery

- H102 performs 532,198 exact equal-shape swaps and reduces HPWL 275.806M ->
  169.667M while overflow is bitwise invariant at 6.24097% in 126.60s.
- H103 alternates 20 coordinate recovery sweeps and ten swap sweeps under the
  same exact cap, targeting below 120M.

## 2026-08-20 H103-H104 feasible swap frontier

- H103 reaches 156.426M HPWL / 6.24094% exact overflow after 20 recovery and
  ten swap sweeps in 129.25s. The below-120M prediction fails.
- H104 widens equal-shape swap radius from 16 to 32 bins while preserving exact
  occupancy by construction.

## 2026-08-20 H104-H105 feasible identity recovery

- H104 performs 40,928 radius-32 equal-shape swaps and reaches 154.113M HPWL
  with invariant 6.24094% overflow in 404.74s. The gain is only 2.31M.
- H105 tests exact incident-net assignment recovery as the final local identity
  operator.

## 2026-08-20 H105 assignment closure

- H105 commits four exact assignment batches and reduces HPWL 154.113M ->
  154.050M with invariant 6.24094% overflow. Local identity recovery is
  saturated; no pseudo-gradient or smoothing is justified.

## 2026-08-20 H106 position-seeded strict path

- The pre-existing non-legalized position-seeded recursive bisection checkpoint
  has 224.399M exact HPWL / 1.40758% exact overflow. H106 starts recovery and
  equal-shape swaps from this checkpoint to test topology preservation.

## 2026-08-20 H106 result and H107 coefficient audit

- H106 reaches 170.935M exact HPWL / 1.4080% exact overflow after 20 exact-cap
  coordinate sweeps and ten density-invariant equal-shape swap sweeps. H105
  remains the better strict HPWL point at 154.050M / 6.24094%.
- Locked H107 to test delayed mu decay, then cancelled the full run after its
  mandatory smoke audit exposed a coefficient-scale confounder. At iteration
  zero, the weighted term norms are HPWL `2.7961e-6`, exact overlap
  `2.0822e-2`, and electrostatic `6.1142e-5`; overlap is nonzero for only
  33.7% of movable nodes.
- Added term-norm and exact-overlap active-fraction diagnostics and reset Adam
  state on entry to the forced `mu=0` final stage. H108 will use one-time
  gradient-norm unit conversion before adaptive dimensionless control.

## 2026-08-20 H108 gradient-balanced result

- H108 reaches 53.526M exact HPWL / 88.1673% exact overflow in 107.62s. The
  fixed coefficient scales are `s_overlap=0.0201429` and
  `s_electro=0.731709`; control value one gives equal initial gradient norms.
- Final overlap and HPWL contribution norms are `1.9300e-6` and `1.4466e-6`,
  so lambda is now dimensionally reasonable. The locked overflow prediction
  fails because mu is reduced from a high-overflow stationary balance.
- H109 will boost mu on high-overflow stagnation up to a fixed bound, then
  perform the decay/lambda handoff and retain a forced exact `mu=0` stage.

## 2026-08-20 H109 boost-then-decay result

- Boosts through `mu_control=32,64,128` reduce the best continuation overflow
  to 85.5066% at 135.123M HPWL. The below-70% prediction fails; each doubling
  gives only about one additional overflow point.
- The locked final `mu=0` trajectory ends at 55.381M / 88.6982% in 110.71s.
  The executable then restores the lower-overflow positive-mu checkpoint, so
  its printed `[Result]` is not accepted as a final nonsmooth result.
- H110 will restore the best continuation point on entry to the final stage,
  reset Adam, and select/report only a `mu=0` iterate. It will also test a
  larger bounded mu range before rejecting electrostatic continuation strength.

## 2026-08-20 H110 strict high-mu result

- H110 reaches 160.483M / 85.0842% in 127.39s. It restores the best positive-mu
  checkpoint before the final stage and reports only a checkpoint audited with
  physical `mu=0`; strict output semantics are now enforced by the driver.
- Mu controls 256, 512, and 1024 all saturate near 85.1% overflow. More scalar
  pressure is rejected.
- Data audit shows the 543 large adaptec1 nodes are fixed terminals, so movable
  area preconditioning is not the dominant missing mechanism. H111 instead
  tests a centroid-only translation anchor, replacing the current per-node
  quadratic that opposes spreading and scales together with electrostatics.

## 2026-08-20 H111 centroid-anchor breakthrough

- H111 reaches a strictly selected 74.031M HPWL / 54.7432% 128-grid overflow
  point in 113.43s. The overflow prediction narrowly fails, but the result is
  28.18 points below strict H81 at near-baseline HPWL.
- The corrected objective records a centroid-only anchor energy and forces it
  to vanish with mu. At iteration 200, overflow is already 66.8%, confirming
  that the old per-node center quadratic caused the high-mu floor.
- Independent 512-grid exact audit is 66.1264%. H112 will use 512 bins for both
  electrostatic continuation and exact overlap, reducing the charge-footprint
  mismatch before invoking any expensive breakpoint stage.

## 2026-08-20 H112 fine-grid result

- H112 reaches 67.939M / 68.4957% on the direct 512-grid metric in 84.05s.
  It fails the below-20% prediction and does not beat H111's 66.1264% audit.
- The fine field transports mass more slowly; using it from iteration zero is
  rejected. Snapshot interval 50 avoids the previous 500 MB text-output cost.
- H113 will start a short 512-grid centroid continuation from H111's strict
  checkpoint, testing coarse-to-fine composition before exact breakpoints.

## 2026-08-20 H113 coarse-to-fine success

- H113 reaches 94.584M exact HPWL / 16.6939% exact 512-grid overflow in
  31.03s incremental time. Before the final stage, fine electrostatics reaches
  25.9499%; the strict mu-zero stage improves it by another 9.26 points.
- The final exact trajectory bottoms early and then sacrifices density for
  HPWL. H114 will restart the same exact objective from the strict checkpoint,
  resetting Adam moments to expose new rectangle-bin boundary pieces.
- Current measured H111+H113 time is 144.46s because H111 wrote snapshots every
  ten steps. The interval-50 optimization measured in H112 is expected to
  recover enough I/O time for the final composed rerun.

## 2026-08-20 H114 exact restart density success

- H114 runs 100 exact-only Adam steps from H113 and selects
  123.241M / 4.18619% in 14.81s. It is the first fast strict sub-7% point from
  the corrected multiresolution homotopy path.
- Iteration 99 is 114.718M / 6.0521%; current output selection minimizes
  overflow lexicographically rather than HPWL under the 7% feasibility band.
- H115 will recover HPWL from the conservative minimum-overflow point under an
  explicit 6.9% exact cap, keeping HPWL epsilon-active as the only active-set
  approximation.

## 2026-08-20 H115 feasible recovery

- Twenty exact-cap recovery sweeps reduce 123.241M to 94.119M HPWL while
  overflow reaches 6.899985%; runtime is 16.60s.
- The below-90M prediction fails but the below-105M branch gate passes. H116
  will run three equal-shape identity-swap sweeps, which preserve the exact
  occupancy map and do not consume the remaining overflow margin.

## 2026-08-20 H116 swap result

- Three equal-shape sweeps perform 58,033 swaps and reach 92.807M HPWL with
  overflow bitwise fixed at 6.899985%. Runtime is 43.34s and the below-76M
  prediction fails.
- Further identical swaps are rejected. H117 tests whether the changed
  identities reopen coordinate recovery under the same exact cap.

## 2026-08-20 H117 alternating branch closure

- H117 lowers 92.807M to only 92.359M at 6.899999% in 14.96s. The below-90M
  branch gate fails and the swap/recovery operator is closed.
- H118 returns to H113's 94.584M / 16.6939% checkpoint and uses a 25% exact
  intermediate cap to recover topology before another exact density restart.

## 2026-08-20 H118 intermediate recovery

- H118 reaches 88.023M / 25.0% after 20 exact-cap sweeps in 13.79s. It misses
  the below-78M prediction.
- H119 repeats H114's exact-only Adam restart from this checkpoint. The branch
  is supported only if its strict sub-7% HPWL beats H114's 123.241M.

## 2026-08-20 H119 exact restart result

- H119 selects 112.424M HPWL / 5.41092% exact overflow in 14.25 seconds.
- The last iterate reaches 101.822M / 6.88191%, but it still fails to beat
  H115's 94.119M strict feasible result. The H118--H119 branch is rejected.
- The run is a strict physical-mu-zero solve; overlap epsilon is zero and no
  electrostatic, smooth, or legalization operation appears in the stage.

## 2026-08-20 H120 long-range recovery result

- Steps 16/32/64 with matched one-bin minimum line-search candidates finish at
  94.054M, 93.994M, and 94.025M HPWL. All exact overflows remain below 6.9%.
- Runtime is 16.40, 16.48, and 18.35 seconds. Step 32 is the new strict a1
  best, but the gain over H115 is only 0.125%; the prediction is rejected.
- H121 will test axis-separated coordinate-subgradient candidates. Lambda and
  mu remain physically zero, and density remains an exact acceptance oracle.

## 2026-08-20 H121 axis-separated recovery result

- H121 reaches 91.3975M HPWL / 6.89999992% exact overflow in 38.50 seconds,
  improving H120-step32 by 2.597M but missing the below-90M prediction.
- The final five sweeps gain only 4.32k, so additional identical sweeps are
  ruled out. H122 will start recovery from H114's physical-mu-zero iteration-75
  snapshot, which has 117.535M / 4.5754% before recovery.

## 2026-08-20 H122 late-checkpoint result

- The loaded iteration-75 snapshot audits at 117.332M / 4.64033%. Recovery
  reaches 91.3746M / 6.899995% in 39.45 seconds, only 22.9k below H121.
- The below-88M prediction fails. H123 will reproduce H114 with an iteration-99
  snapshot and make one final checkpoint-selection test before closing it.

## 2026-08-20 H123 final-checkpoint result and controller audit

- H123's replay matches H114 within tolerance. Recovery reaches 91.4087M /
  6.9000% in 38.84 seconds, worse than H122, closing checkpoint timing.
- The final homotopy stage currently forces lambda to its maximum and excludes
  it from adaptive updates. H124 will use a 6.5--7.0% exact-overflow hysteresis
  band and choose the lowest-HPWL feasible checkpoint, always with physical
  mu zero.

## 2026-08-20 H124 adaptive final-lambda result

- H124 selects 105.937M / 4.99631% in 20.40 seconds, passing all numerical
  predictions and improving the old H114 selected output by 17.304M.
- Lambda control decreases from 150 to 90. Overflow does not reach 7%, so the
  increasing branch remains unexercised. H125 will recover this new exact
  topology under the 6.9% cap.

## 2026-08-20 H125 adaptive-topology recovery result

- H125 reaches 91.1194M / 6.899998% in 24.59 seconds, a new strict best but
  only 0.255M below H122. The below-88M prediction fails.
- Further identical cap recovery is rejected. H126 extends the exact adaptive
  lambda stage to 300 steps to exercise the complete hysteresis controller.

## 2026-08-20 H126 extended lambda result

- H126 selects 95.4261M / 6.33459% in 29.81 seconds. Lambda control decreases
  from 150 to 10; all iterates have physical mu zero.
- Overflow has begun rising sharply but does not cross 7% by iteration 299.
  H127 extends the deterministic trajectory to 400 steps to test the upper
  feedback branch.

## 2026-08-20 H127 band-regulation result and outer-loop pivot

- H127 selects 93.6776M / 6.99568% in 40.11 seconds. Lambda decreases from
  150 to 10, then rises to 15 and 20 after upper-band crossings.
- The controller is validated. H113 audit identifies the upstream problem:
  150 fine-grid steps keep mu at 16 and lambda at zero, causing most topology
  damage before the abrupt exact-stage handoff. H128 tests takeover at iter25.

## 2026-08-20 H128 early-handoff rejection

- H128 cannot cross the exact-overlap sparsity plateau and falls back to
  105.516M / 11.7755% in 44.22 seconds.
- H129 will trigger homotopy handoff directly when overflow reaches a coarse
  threshold: mu then decays every window while exact-overlap lambda rises.
  Final selection remains exclusively physical-mu-zero.

## 2026-08-20 H129 continuous-handoff rate failure

- H129 logs the intended handoff, but overflow rebounds from 28.38% to 35.29%
  as mu control falls to 0.184 before lambda reaches its maximum.
- Strict output is 92.048M / 21.9608% in 28.59 seconds. H130 changes only the
  measured rate imbalance: mu decay 0.9 and lambda step 15.

## 2026-08-20 H130 rate-balanced handoff result

- H130 reaches 90.590M / 19.8564% in 28.32 seconds, improving H129 but missing
  the below-10% gate.
- Overflow reverses near iter130 while mu continues automatic decay. H131 adds
  an exact-overflow rebound guard that restores one mu step on worsening
  windows; final-stage mu remains forced to zero.

## 2026-08-20 H131 rebound-guard rejection

- H131 repeats 90.590M / 19.8564% in 37.57 seconds. Mu oscillates but creates
  no new density minimum, rejecting scalar one-step rebound control.
- H132 keeps mu close to its transport scale with decay 0.99 while lambda
  reaches 150, then forces mu to zero in the reported final stage.

## 2026-08-20 H132 strong coexistence result

- H132 reaches 92.462M / 11.1664% at the mu-zero boundary in 40.77 seconds,
  breaking the 20% plateau but missing strict feasibility.
- H133 doubles mu start/max control to 32 under the same slow-decay schedule,
  with final output still restricted to physical mu zero.

## 2026-08-21 H133 doubled bridge result

- H133 reaches 96.495M / 7.6153% at the mu-zero boundary in 42.89 seconds.
- H134 will test mu control 40 to close the remaining 0.615-point gap.

## 2026-08-21 H150-H160 fine recovery continuation

- H150 stability handoff starts exact lambda on the H140 `iter_1050` branch,
  reducing exact overflow to 47.8139% but never reaching feasibility.
- H151/H152 test mu control 32 with and without forced final lambda; they end
  at 46.0968% and 41.8875% overflow. H153 starts lambda at 80% overflow and
  ends at 40.9002%. The H140 `iter_1050` fine checkpoint is rejected.
- H154 long-radius exact recovery from H121 improves 91.3975M to 91.3835M.
- H155 breakpoint-oracle recovery reaches 90.5744M / 6.899997%, H156 raw
  full-degree directions reach 90.4927M / 6.9%, and H157 equal-shape swaps
  reach 90.1383M / 6.9%.
- H158 wider swaps reach 89.9435M / 6.9%. H159 with a still-valid 6.99%
  exact cap reaches 89.6829M; H160 reaches the current best 89.6416M / 6.99%
  after density-neutral exchange. Breakpoint and swap operators are now
  saturated on this a1 topology; the next independent run is fresh a2.

## 2026-08-21 H161-H167 fresh transfer

- H161 fresh a2 coarse (raw `.pl`, no checkpoint) ends at 78.129M / 79.57%
  in 57.19s. H162 fine mu=40 reduces exact overflow to 8.3189% but selects
  144.899M HPWL at mu=0.
- H163 found that exact recovery rejected infeasible seeds above the cap. The
  guard was corrected so above-cap moves must be exact HPWL-decreasing and
  overlap-nonincreasing; H164 reaches 135.277M / 6.99% in 29.62s.
- H165 equal-shape exchange reaches 134.965M / 6.99% in 32.93s.
- H166 fresh a3 weak-mu coarse ends at 299.737M / 97.17% in 103.74s. H167
  fresh direct mu=40 fine transport collapses to boundary-heavy coordinates,
  ending at 579.757M / 95.96% in 36.86s. A3 is paused for a mechanism-level
  transport schedule redesign.

## 2026-08-21 H168-H171 bounded transport and a4

- H168 fresh a3 mu=16 with per-node anchor 0.01 keeps HPWL near 245M but ends
  at 99.92% exact overflow; strong anchoring alone does not activate exact
  density transport.
- H169 fresh a4 coarse reaches 204.355M / 59.97% in 110.08s. H170 fine mu=40
  reaches 476.276M / 6.814%, satisfying overflow but damaging HPWL.
- H171 exact breakpoint recovery from H170 reaches 436.288M / 6.98999% in
  47.56s. The result is feasible but fails HPWL quality; a3/a4 remain open.
## 2026-08-21 H191-H192 exact rigid-block experiments

- Added atomic exact group-overlap deltas and net-disjoint rigid block line
  search with scaled and exact-breakpoint candidates. Core tests pass.
- H191 block-only screening improved H186 to 112.4245M/8.2330%, but mixed and
  prefixed block/node recovery ended at 104.2911M and 104.4038M at 7%, both
  worse than H187 node-only.
- Feasible-band H191 polishing reached 103.6440M/6.99998% in 5.12 s. Including
  the fresh H182 and H186 stages, total runtime is about 112.1 s; HPWL remains
  41.55% above the paper baseline, so no a2-a4 transfer was run.
- H192 added a true zero-epsilon exact-density block direction and was
  rejected at 103.6958M/7.0000%. It confirms that aggregating sparse active
  edges alone does not repair the global transport topology.

## 2026-08-21 H193 regional active-set continuation

- Added checkpoint-aligned regional price state and parent-to-child exact
  replication across 64/128/256/512 grids. Each level still rebuilds exact
  overlap and fixed macro capacity; lambda is reinitialized per level.
- Fresh a1 reaches 107.190M/72.8946%, versus H175's 87.621M/79.305%. The
  6.41-point overflow gain supports active-state mapping, but remains far from
  feasibility. Full process wall time is about 37 seconds.
- H194 doubles continuation length and reaches 188.822M/60.4535% in about 73
  seconds. It fails both quality gates; coarse lambda grows above 1e6 without
  clearing its flat overlap plateau. Blindly extending regional continuation
  is rejected.

## 2026-08-21 H195 bounded exact coarse transport audit

- Locked H195 before execution: raw a1 `.pl`, one 64x64 capacity-flow pass,
  exact one-bin anchor steps, 40 threads, and a 5% HPWL boundary prediction.
- Modified exact-anchor coarse flow to take only one coarse-bin step toward a
  sink; this avoids the old full-distance topology shock and preserves all
  existing unit tests.
- The run was CPU-bound for over ten minutes with no completed metric output,
  so it was terminated on runtime grounds. Exhaustive fine-bin screening is
  not viable at full source-node scale.
- Decision: retain global transport as a possible escape mechanism, but
  redesign it as bounded net-connected commodities with a fixed shortlist and
  exact batch audit. Do not migrate H195 to a2-a4.

## 2026-08-21 H196-H197 coarse-flow consistency correction

- H196 completed in 18.87s with eight source-bin and 128-node caps, but
  accepted zero exact moves. The raw exact state was 0% overflow while the
  center-bin coarse load falsely reported 99.925%.
- H197 changed coarse load to exact rectangle-to-coarse-bin intersection
  accumulation. The coarse and exact initial overflow now agree at zero and
  no flow edge is generated. Core tests remain green.
- The one GP update from raw still jumps to 99.995% exact overflow, so the
  next useful mechanism is an exact-overflow guarded HPWL recovery, not more
  coarse density flow on raw initialization.

## 2026-08-21 H198 raw HPWL recovery domain audit

- H198 recovery lowers a1 HPWL 104.924M -> 102.776M in 27.7s with reported
  zero incremental overflow, but raw ISPD `.pl` movable coordinates are all
  `(0,0)`, outside the core rows.
- Exact overlap intentionally ignores fully out-of-domain rectangles, so the
  apparent 0% raw overflow is artificial. Final boundary clamping yields
  99.832% overflow; the branch is rejected as an initialization method.
- The pre-GP exact-feasible selector fix remains in the code for valid staged
  states. Fresh runs must use deterministic in-domain initialization derived
  from the raw benchmark, with no prior checkpoint or ePlace seed.

## 2026-08-21 H199 direct fine homotopy ablation

- H199 starts fresh from the deterministic center Gaussian and directly runs
  the H186-style 512-bin continuation for 400 iterations with 40 threads.
- It selects 154.438631M HPWL / 70.074253% exact overflow in 22.0865 seconds,
  versus H182->H186's 112.847669M / 8.304501%. Direct fine continuation is
  rejected; the coarse prefix provides necessary low-frequency transport.
- An exploratory global movable-centroid contraction produced no exact-HPWL-
  improving candidate from H186. Fixed-terminal nets oppose uniform scaling,
  motivating a fixed-terminal-aware compact contraction rather than allowing
  HPWL damage.

## 2026-08-21 H200 fixed-terminal-aware compact contraction

- H200 accepts zero fixed-aware global contraction candidates from H186; its
  unchanged exact recovery reaches 104.560M / 7.0000%.
- Global affine compactness is rejected as too rigid. H202 moves the same idea
  to local incident-net/support-centroid directions, while retaining exact
  HPWL/overlap acceptance and the original epsilon-active direction.

## 2026-08-21 H202 local compact directions

- H202 reaches 103.217761M / 7.000000% in 45.28 seconds, improving H187 by
  0.497491M at the same exact cap. 499,277 of 673,954 accepted moves select an
  added compact direction.
- The mechanism is supported and total H182->H186->H202 time remains about
  122.5 seconds, below 2x the 67-second a1 baseline, but HPWL remains 40.97%
  high. H203 now screens shorter coarse prefixes to free runtime and possibly
  preserve a less dispersed topology.

## 2026-08-21 H203 shortened checkpoint screen

- H182 iterations 300/400/500 followed by unchanged H186 reach 35.74%,
  35.31%, and 31.74% overflow. All fail the 12% gate despite progressively
  lower HPWL.
- These checkpoints omit H182's final coarse physical-mu-zero stage. H204 will
  fresh-run a 12-stage coarse schedule so the exact handoff occurs at step 550,
  testing handoff structure rather than another pre-handoff checkpoint.

## 2026-08-21 H204 fresh short coarse handoff

- The 600-step fresh short coarse stage reaches 64.448M / 74.104%; unchanged
  fine continuation reaches 118.934M / 10.752%. The early exact handoff is
  essential but the 12-stage chain misses its 10% gate and is rejected.
- H205 will test exact non-rigid contraction of small net blocks from H202's
  feasible checkpoint. Unlike H191 rigid translation, this directly reduces
  internal net span while every group remains exact-cost audited.

## 2026-08-21 H205 exact net contraction

- H205 accepts 704 non-rigid contraction groups and lowers H202 from
  103.217761M to 103.202780M at 7% overflow in 1.087 recovery seconds.
- The 0.015M gain fails the pre-registered 0.20M gate. Net contraction is only
  a cheap polishing operator; H186's main HPWL damage must be prevented while
  fine transport is still active, not repaired after reaching the 7% cap.

## 2026-08-21 H206 late fine checkpoint recovery

- H186 iteration 325 audits at 107.617M / 9.384%; 20 exact compact-direction
  sweeps reach 101.285M / 7.0000%. This improves H202 by 1.933M and validates
  avoiding H186's high-HPWL selector restore.
- Recovery takes 87.28 seconds and fails the complete runtime gate. Iteration
  349 was not saved and its attempted load failed before execution; iteration
  350 is not substituted because it follows the selector transition.
- Sweep 7 already reaches 101.532M / 7%, so H207 tests the exact eight-sweep
  truncation as a quality-preserving runtime correction.

## 2026-08-21 H207 eight-sweep late recovery

- H207 reaches 101.531876M / 6.999984% in 36.93 wrapper seconds. Combining
  H182 and the H186 iteration-325 prefix gives 110.60 seconds and about 1535
  iterations, passing overflow, runtime, and iteration gates.
- It becomes the preferred H186 branch, but HPWL remains 38.66% above the
  paper baseline. H208 will replay H186 with every-iteration snapshots and
  test the lower-HPWL iteration-349 topology that was not originally saved.

## 2026-08-21 H208 replayed iteration-349 recovery

- The deterministic fine replay reproduces original iteration 349 within
  0.1%. Eight exact sweeps reach 100.909M / 7.117% in about 20 wrapper seconds.
- HPWL passes but overflow misses by 0.117 percentage points. H209 changes only
  recovery length from eight to ten sweeps to cross the exact 7% boundary.

## 2026-08-21 H209 ten-sweep late recovery

- Ten exact sweeps reach 100.842941M / 7.101991% in 25.61 wrapper seconds.
  HPWL and runtime pass, but overflow remains 0.102 percentage points high.
- The final two sweeps gain 65.6k HPWL but only 0.0155 overflow percentage
  points, refuting the predicted quick crossing. H210 is the final bounded
  unchanged extension at twelve sweeps.

## 2026-08-21 H210 twelve-sweep late recovery

- Twelve sweeps reach 100.807650M / 7.096834% in 27.81 wrapper seconds.
- The two added sweeps gain only 35.3k HPWL and 0.00516 overflow percentage
  points. Identical late-sweep extension is closed; the next intervention
  moves earlier in the fine trajectory and remains exact-audited.

## 2026-08-21 H211 center-anchor rejection

- A weak per-node center anchor reaches only 141.646M / 10.0131% at its
  selected mu-zero checkpoint. Iteration 349 is 113.760M / 13.0936%, worse
  than H186 in both metrics.
- A common center conflicts with fixed-terminal net geometry. H212 replaces
  it with a vanishing separable proximal term around each H182 input position.

## 2026-08-21 H212 initial-position proximal support

- Added an `initial` anchor mode whose exact quadratic term uses each node's
  H182 input coordinate and is multiplied by physical mu. Core tests pass.
- Iteration 349 reaches 103.097M / 11.2391%, improving H186 by 3.32M HPWL at
  1.57 overflow points higher. Iteration 325 is 104.482M / 10.7517% in 39.58s.
- H213 applies bounded exact mu-zero recovery to the saved iter325 checkpoint.

## 2026-08-21 H213 proximal checkpoint recovery

- Twelve exact sweeps reach 97.904M / 7.6110% from H212 iteration 325. The
  2.90M HPWL improvement over H210 validates the proximal topology advantage.
- Recovery takes 52.96 seconds and overflow stalls above the cap. H214 raises
  only the proximal homotopy's exact-overlap lambda maximum from 150 to 225,
  shifting density work into the faster simultaneous phase.

## 2026-08-21 H214 lambda-ceiling rejection

- Raising only lambda maximum to 225 leaves iteration 325 at 104.495M /
  10.7320%. Lambda control is only 154.66 there, so the old ceiling was not the
  early bottleneck.
- H215 keeps the ceiling and doubles only the adaptive lambda step to 60 to
  test earlier exact-overlap takeover.

## 2026-08-21 H215 earlier lambda takeover

- Iteration 325 reaches 106.392M / 9.6253% in 38.64 seconds. This materially
  reduces H212's density debt but misses the locked 9.5% gate by 0.125 points.
- At nearly H208's HPWL it has 0.125 points lower overflow. H216 separately
  tests twelve-sweep exact recovery under the complete staged runtime gate.

## 2026-08-21 H216 quality pass, runtime near-miss

- Twelve exact sweeps reach 98.4347M / 7.0000%, the best strict HPWL on the
  H186-derived branch, but complete staged time is about 137.95s versus 134s.
- The cap is first reached at sweep two and sweep six is already 98.8624M.
  H217 truncates only the recovery length to six sweeps.

## 2026-08-21 H217 strict proximal branch pass

- Six exact sweeps reach 98.8624M / 7.0000% in 27.26 seconds. Complete staged
  time is about 122.48 seconds with roughly 1,533 iterations, passing all three
  non-HPWL gates.
- This is the preferred H186-derived branch but remains 35.02% above the paper
  a1 HPWL. The next outer-loop comparison returns to the stronger H160 lineage.

## 2026-08-21 H160 provenance audit

- H160 traces through H121 -> H114 -> H113 -> H111. H111's locked protocol
  explicitly starts from a non-legalized ePlace placement, which violates the
  later raw-only initialization constraint.
- H160 remains historical evidence but is removed from valid-best status.
  H217 is the current fresh-from-raw strict a1 best. H218 revisits the legal
  H177 HPWL-only seed before the successful coarse homotopy chain.

## 2026-08-21 H218 HPWL-seeded coarse failure

- Fresh HPWL seeding reaches 43.049M, but unchanged coarse homotopy ends at
  51.377M / 85.6315%. Combined runtime is 125.22 seconds.
- The seed changes electrostatic gradient scale from 2.356 to 0.4398; unchanged
  mu therefore has only 18.7% of H182's effective strength and lambda never
  activates. H219 calibrates mu control to 5.36 and uses 16 machine-native
  threads, still below the 40-thread limit.

## 2026-08-21 H219 scale-calibrated seeded coarse result

- Fresh HPWL seed reaches 43.035M in 5.80s. Calibrated coarse homotopy reaches
  60.621M / 71.6920%, but full combined prefix time 80.44s misses its 75s gate.
- Saved iteration 1100 is already 60.618M / 71.7078% at a combined 70.24s.
  H220 starts proximal fine continuation there and targets iteration 300.

## 2026-08-21 H220 seeded fine scale failure

- Iteration 300 reaches 83.043M / 36.6836%; selected output is 88.981M /
  25.5722%. The seeded topology keeps HPWL low but density misses badly.
- Versus H215, fine electrostatic and overlap gradient scales are 12.05x and
  8.66x smaller. H221 converts those measured ratios into mu/lambda controls,
  preserving the exact objective and final mu-zero contract.

## 2026-08-21 H221 calibrated seeded fine result

- Iteration 300 is 99.399M / 10.9646%, failing its 10% gate. Physical-scale
  calibration nevertheless restores a mu-zero entry at 109.426M / 7.7516%.
- The complete fresh prefix to iteration 350 is about 90.13s, leaving 43.87s.
  H222 tests eight exact recovery sweeps from this saved checkpoint.

## 2026-08-21 H222 snapshot provenance failure

- H222 loaded `h221_a1_fine/snapshots/iter_350.pl` and audited its actual
  coordinates at 107.869M HPWL / 16.8641% exact overflow. H221's metrics row
  at iter 350 instead reports the post-update state 109.426M / 7.7516%.
- Eight exact recovery sweeps reached 95.1609M / 8.0606% in 33.20s, but this
  is not a valid recovery of the H221 mu-zero entry and fails the overflow
  gate. The mismatch exposes a snapshot writer bug: metrics are written for
  pre-update coordinates while the same-named `.pl` contains post-update
  coordinates.
- H223 is pre-registered to use H221 `best.pl`, written by the final exact
  selector, and requires a strict initial-metric audit before accepting any
  result.

## 2026-08-21 H223 corrected seeded recovery

- The initial audit from H221 `best.pl` passed at 109.425511M / 7.803537%
  exact overflow, within the locked 0.1-point provenance tolerance.
- Six exact breakpoint/compact recovery sweeps reached 95.686353M /
  6.999997% in 16.5208s (19.8146s wrapper). The estimated complete staged
  prefix is 109.95s and the iteration budget remains below 5000, so H223 is
  the new strict a1 best under the raw-to-staged contract.
- HPWL is still 30.68% above the 73.22M paper GP baseline. The next work must
  target topology preservation during density transport; extending this same
  late recovery operator is low leverage.

## 2026-08-21 H224 net-aware cap polishing

- Starting from H223's exact-feasible output (95.686353M / 6.9999971%), two
  sweeps retained node moves and added bounded net-aware rigid blocks with the
  exact density direction. 1,580 block moves were accepted.
- The final exact state is 95.360654M / 6.9999973% in 9.41s wrapper time.
  The estimated complete staged runtime is 119.36s, below the 134s gate.
- This validates net-aware batch/block directions as an additional search
  direction under the unchanged nonsmooth objective. The HPWL gap to paper is
  still 30.24%, so the mechanism is polishing rather than the missing global
  topology transport.

## 2026-08-21 H225 early net-aware recovery

- From the same audited H221 fine entry, enabling exact-audited net-aware
  blocks during all six recovery sweeps accepted 10,809 block moves and
  reached 95.247422M / 6.9999988% in 31.29s wrapper time.
- H225 beats late-only H224 by 0.113232M HPWL. The estimated fresh staged time
  is 121.42s, under the 134s gate. This supports inserting topology-preserving
  group directions during the large early recovery moves.

## 2026-08-21 H226 wider net-block ablation

- H226 continued from H225 at 95.247422M / 6.9999988% and widened exact
  net-block candidates to degree 32 and 16 movable nodes. Two sweeps accepted
  1,766 block moves and reached 94.952654M / 6.9999996%.
- Added wrapper time was 10.31s; the complete H221 -> H225 -> H226 staged
  estimate is 131.73s, just below the 134s gate. No further identical sweep
  fits the locked runtime budget. The strict a1 gap to paper is 29.68%.

## 2026-08-21 H227 fine net-batch homotopy rejection

- Adding a per-net mean exact HPWL subgradient (weight 0.25) inside the H221
  fine homotopy raised fine runtime to 48.95s and selected 115.036M / 6.6492%
  exact overflow. The best feasible continuation checkpoint was 111.272M /
  6.9932%, both worse than H221's 109.426M / 7.7516% handoff.
- The full-net aggregation also pushes the H226 staged chain beyond the 134s
  gate. H227 is rejected; net-aware directions remain a recovery-only
  mechanism where exact candidate auditing is cheap enough.

## 2026-08-21 H228 normalized net-batch timeout

- H228 changed the homotopy direction to a normalized 0.10 blend and limited
  batching to degree-16 nets. The fine stage still reached 58.7s by iter 343
  (about 99.16M / 10.41%) and was stopped at the locked 30s gate; no final
  selector was produced.
- This closes per-iteration net-batch inside homotopy under the current
  runtime budget. The exact-audited net-block recovery mechanism remains the
  active topology-preserving method.

## 2026-08-21 H229 cached net-batch timeout

- Refreshing the normalized degree-16 batch direction every 25 iterations
  still took 48.68s by iter 340 (about 99.36M / 10.27%); the run was stopped
  at the 30s fine-stage gate without a final selector.
- Caching is not sufficient for the current hypergraph/runtime budget. H229
  is closed as a timeout, and net-aware directions remain recovery-only.

## 2026-08-22 H230 relaxed net recovery

- Relaxing the exact recovery cap from 7% to 15% for three sweeps, with shared
  net-density block directions, reduced 109.426M / 7.8035% to 93.046M /
  15.0000%. This confirms that a bounded exact overflow excursion can repair
  topology and lower HPWL substantially.

## 2026-08-22 H231 old-rule retightening failure

- Reapplying a 7% cap with the old HPWL-only recovery rule reduced overflow
  only from 15% to 14.0564% in three sweeps. Above the cap, the rule forbids
  HPWL-increasing density moves, so it is not the proposed hysteresis method.

## 2026-08-22 H232 lambda-scale failure

- Joint exact GP from H230 with inherited scale 8e-5 had effective lambda
  1.4e-5. It became HPWL-only (61.1905M / 82.8621%), proving scale must be
  recalibrated after a relaxed topology change.

## 2026-08-22 H233 calibrated retightening bridge

- Raising density scale to 0.5 produced effective lambda 0.08796 and reached
  85.8982M / 13.2145% after 100 shared-net-density iterations. It was an
  infeasible but useful bridge for stronger lambda takeover.

## 2026-08-22 H234 hysteresis strict pass

- Continuing H233 with density scale 5.0 and step .001 allowed temporary HPWL
  increases while lambda grew. Iteration 99 selected **89.2995M / 6.5967%**.
- The H221 -> H230 -> H233 -> H234 staged estimate is 130.32s, under 134s;
  this is the new strict a1 best. The HPWL gap to paper falls to 21.96%.

## 2026-08-22 H235 direct strong retightening

- Skipping H233 and applying scale 5.0 directly from H230 reached 92.9319M /
  4.3050% in 10.16s. It is feasible but 3.63M worse than H234, confirming
  that a moderate-lambda topology bridge is valuable.

## 2026-08-22 H236/H237 20% relaxed-band ablation

- A 20% relaxed recovery band reached 91.7828M / 19.9999%, only 1.263M below
  the 15% H230 endpoint.
- Strong exact retightening from that state selected 91.6573M / 5.9978%,
  2.358M worse than H234. The wider band creates too much retightening debt
  and is closed; 15% remains the preferred hysteresis cap.

## 2026-08-22 H238-H245 hysteresis schedule refinement

- Switching from H233 iteration 60 for 100 strong steps (H238) retained lower
  HPWL but stopped at 88.5832M / 7.6405%. Extending it to 140 steps (H241)
  exposed a stable exact-active-set plateau near 7.44%; iteration 60 is closed.
- The iteration-80 handoff reached 88.9831M / 6.9954% with scale 5 (H240),
  while iteration 90 was worse at 89.1449M / 6.8089% (H242). The useful
  moderate-to-strong switch is therefore near iteration 80.
- Scale 4 and 4.5 retained lower HPWL but missed feasibility at 7.1453% and
  7.0477%. Scale 4.75 crossed at H245 iteration 114 and selected 88.9171M /
  6.99968%, improving H234 without increasing the 200-step bridge/strong
  budget.

## 2026-08-22 H246 exact feasible-band polish

- One exact breakpoint/node plus net-rigid sweep from H245 accepted 116,807
  node and 2,967 block moves (119,774 total), reaching **87.5800M /
  6.999998%** in 6.0719 GP seconds.
- The established complete raw-pl staged estimate is 133.94s and about 1,699
  iterations, so H246 passes the non-HPWL gates with only 0.06s time margin.
- H246 is the new strict a1 best. It remains 19.61% above the paper HPWL; the
  next work must create runtime margin before adding another polish sweep.

## 2026-08-22 H247 fresh 40-thread replay rejection

- H247 regenerated the raw-pl HPWL seed and reproduced 43.0348M / 96.0563%,
  but 40 requested threads increased seed GP time from 5.80s to 11.66s on
  this host.
- Coarse homotopy had already used about 95s by iteration 930, compared with
  74.64s for all 1200 historical 16-thread iterations. The projected first
  two stages approach 126s, leaving no feasible budget for the remaining five
  stages, so the run was terminated under its locked runtime gate.
- This closes 40-thread oversubscription on the current hardware. H246 remains
  strict best; runtime must be saved algorithmically from the 16-thread prefix.

## 2026-08-22 H248 early coarse handoff rejection

- Replacing H219 iteration 1100 by iteration 1000 and replaying the unchanged
  H221 fine stage selected 110.160M / 7.9948%.
- The result narrowly misses the locked 110M screen, and the measured fine
  wall time is 47.12s versus H221's recorded 22.84s. It creates no reliable
  runtime margin, so no downstream continuation is run.

## 2026-08-22 H261-H272 continuous late-stage switch study

- Added an optional in-process late-stage switch that multiplies step/lambda
  without changing exact objectives; defaults preserve historical behavior.
- H261/H262 showed that preserving Adam moments stalls at 10.31%/10.23%
  overflow even with lambda x4. Exact one-bin and eight-bin recovery (H263,
  H264) could not cross the active-set plateau.
- H265 reset Adam only at iteration 81 and reached 89.262M / 6.6227%; H266's
  four exact node-only sweeps reached 87.745M / 7.00%. A second block (H267),
  six-trial search (H268), and late net blocks (H269) were rejected.
- H270 moved the Adam reset to iteration 70 and improved the feasible
  checkpoint to 88.497M / 6.9967%. H271 reached 87.248M / 7.00% after four
  exact sweeps; H272's fifth sweep added only 2.8k and ended at 87.220M.
- The continuous switch is a promising runtime-saving branch but still does
  not satisfy the paper HPWL gate; H257 remains strict quality best.

## 2026-08-22 H273-H279 conservative-step and near-cap recovery

- H273 lowers the late step to 0.4 but ends at 7.2955% overflow. Four exact
  sweeps (H275) reach 87.118M but stop at 7.0618%; four more (H276) remove
  only 0.0135 overflow points.
- A 7.5% exact excursion (H277) lowers HPWL to 86.976M, but retightening
  (H278) stalls at 7.4538%. This small hysteresis is rejected.
- One near-cap mixed net-block/node sweep (H279) reaches 87.041M / 7.0396%,
  still infeasible. The active-set rescue must occur earlier in continuation.

## 2026-08-22 H280-H281 prefix/timing ablations

- H280 resets Adam at iteration 60 and ends at 7.5316%, closing an earlier
  takeover under the two-sweep prefix.
- H281 restores H230's third relaxed sweep, but the iteration-70 continuous
  takeover ends at 7.6340%; extra relaxed debt outweighs the 1.5M bridge gain.

## 2026-08-22 H282-H284 global candidate experiments

- H282 auction/group transport reduced overflow from 15.00% to 14.56% but
  increased HPWL from 94.55M to 100.25M; final 89.527M / 6.673% was rejected.
- H283 bounded coarse flow with exact anchors accepted 1,987 moves, raised the
  prefix to 95.39M, and selected 88.859M / 6.786%; rejected.
- H284 recursive hypergraph bisection/FM reduced overflow below 1% but raised
  HPWL to 133.596M; rejected. Global candidates need a total-HPWL guard and
  topology-preserving displacement before further use.
## 2026-08-22 H285-H288 bridge and density-neutral global polish

- Pre-registered H285 to test exact batch backtracking with a 2% temporary
  HPWL budget in the 15% bridge.  It reached 5.4978% exact overflow and
  93.565M HPWL in 35.65 s, but topology debt remained; the direction was not
  accepted as a solution.
- H286 narrowed the budget to 0.5% under the same protocol.  It reached
  5.4143% and 93.561M in 35.36 s, statistically indistinguishable in HPWL from
  H285 and still far worse than H257.  Scalar bridge-budget tuning is closed.
- H287 tested exact equal-shape exchanges from H257.  Five sweeps accepted
  33,523 swaps and improved 87.119465M to 86.776049M while preserving exact
  overflow 6.9999578%; incremental runtime was 21.75 s.  The end-to-end chain
  is about 155.52 s, so H287 is a quality record but not runtime-valid; H257
  remains the runtime-valid strict record.
- H288 tested two guarded shape-class Hilbert permutations from H287.  All
  candidates rolled back under the exact HPWL guard; the operator is closed for
  this seed.
- Outer-loop decision: global transport/flow/bisection and scalar HPWL-budget
  bridges are rejected because they spend topology.  The next locked direction
  is a local relative-coordinate net-block bridge whose candidates are ranked
  and accepted by exact combined HPWL plus exact overflow-area change.

## 2026-08-22 H289-H294 exact combined net-block bridge

- Added an optional exact recovery score `DeltaHPWL + w*Delta(overflow
  fraction)`; default `w=0` preserves old behavior and tests remain green.
- H289 (`w=200M`, two block sweeps) reduced exact overflow 15.00% to 6.92% in
  2.17s but raised HPWL 4.14M.  H290 (`w=100M`, one sweep) stopped at 10.83%
  and 96.13M, passing its local bridge screen.
- H291/H292 showed that the executable's one-process stage order runs recovery
  and swaps before GP; both selected about 90.65M/5.03% and were rejected.
- H293 split out a true GP-only takeover and selected 91.58M/4.75% in 9.75s.
  H294 then ran exact post-GP recovery plus equal-shape exchange and reached
  88.08M/6.99991% in 54.13s, still 1.30M worse than H287.
- The moderate relative net-block bridge family is closed.  H287 remains the
  strict a1 record; future work must protect topology earlier instead of
  increasing exact density transport strength.

## 2026-08-22 H295-H297 outer-loop continuation

- H295 affine expansion and H296 order-preserving power-law warps were
  negative controls. Despite preserving coordinate order, they left exact
  overflow above 86% and are closed.
- H297 added `--density-net-share`: exact density directions are averaged over
  incident low-degree nets before the optimizer step. This is a direction
  heuristic only; exact HPWL/overlap auditing and the nonsmooth objective are
  unchanged.
- Fresh a1 runs at 501 iterations showed a clear early gain: degree-16 full
  sharing reached 87.68M/72.59% versus 91.09M/78.07% node-wise. Degree 64
  reached 71.88% but with 88.77M HPWL. At 1001 iterations degree-16 sharing
  reached 152.4M/55.59%, confirming that the direction alone accumulates
  topology debt. The next hypothesis is a fixed-macro-aware net-cluster
  capacity assignment guarded by cumulative exact HPWL.

## 2026-08-23 H300-H306 outer-loop continuation

- H300's HPWL-orthogonal projection of shared density directions was worse
  than H297 and is closed.
- H301/H302 confirmed that turning sharing off later recovers HPWL but leaves
  overflow near 70% or worse; larger density scales only increase topology
  debt.
- H303's wider exact equal-shape exchange reduced HPWL by 0.198M at unchanged
  overflow, but required 115.8s incremental time and is closed under the
  runtime gate.
- H304's first-250-iteration sharing followed by Adam reset and exact node
  takeover ended at 66.86M/83.83%; early sharing is not a standalone capacity
  solution.
- H305 auction transport reduced overflow to 58.12% but raised HPWL to
  366.77M. H306's exact global HPWL guard rolled back all moves at 5% and 10%
  budgets. Individual auction transport is closed; future proposals must move
  connected clusters and pass both incident-net and global exact screens.

## 2026-08-23 H307-H308 structural negative results

- H307's local weighted connected groups reduced overflow by only 4.85 points
  while adding 38.8M HPWL; the exact global guard rolled the same round back.
- H308 attempted density-guided equal-shape swaps and accepted none. This is
  mathematically necessary: identical rectangles exchanged between two
  coordinates leave the exact occupancy field invariant. Density-neutral
  swaps can polish HPWL but cannot repair capacity.
- Core tests pass after adding the early-sharing schedule and exact transport
  rollback guard. The next operator must move unequal area or translate a
  connected cluster in smaller, continuously screened batches.

## 2026-08-23 H309-H311 unequal-shape and micro-batch continuation

- H309 accepted 5,187 unequal-shape exact exchanges, reducing overflow by
  0.50 points for a 0.657M HPWL increase. The mechanism is retained only as
  an early capacity micro-operator.
- H310's one 100-node connected-cluster batch stayed within a 1% exact global
  HPWL guard and gained 0.0237 overflow points. H311 chained ten batches and
  gained 0.259 points for 0.706M HPWL. This is safe but too weak for direct
  feasibility; the next hypothesis is hierarchical area transfer with larger
  clusters and exact rollback at each coarse target.
## 2026-08-24 H312-H315 coordinate-descent pilot

- Implemented `density_coordinate` exact overlap search with axis moves,
  breakpoint candidates, cumulative HPWL budget, and optional rigid net blocks.
- Confirmed on fresh a1 that node-wise coordinate descent is too weak; rigid
  net blocks at coarse 64x64 are materially better but still create HPWL debt.
- Tested underfull-bin target ranking. It did not materially improve over
  untargeted block directions, indicating that the missing component is not a
  direction selector alone but a capacity assignment that preserves cluster
  topology globally.
- `mingw32-make test` passes after adding the new source to the build graph.
- Main route remains H257/H287; H312-H315 are exploratory primitives only.

## 2026-08-24 H316-H319 continuation tuning

- Coarse block step size and HPWL budget were increased to test whether the
  capacity-transfer stage was simply too conservative. Larger steps improved
  coarse overflow but incurred topology debt that exact recovery could only
  partially repay.
- Fine Adam, regional prices, early net-sharing, and Adagrad variants were
  tested from the same H316 staged checkpoint. Adam with default density scale
  gave the best HPWL/overflow compromise, but none approached the 7% gate.
- Decision: stop tuning scalar lambda/optimizer on this bridge family. The
  next implementation must make cluster destination assignment capacity-aware
  before movement, rather than choosing a direction after a topology-breaking
  translation.

## 2026-08-24 H320 dynamic residual-capacity assignment

- Added exact coarse residual-capacity accounting to the net-block coordinate
  phase. The ledger aggregates movable occupancy and fixed macro capacity and
  is updated from every committed exact `DensityMove`.
- H320 improves the coarse bridge from 94.38% to 92.87% overflow in one sweep;
  two sweeps reach 89.81% but spend 106.33M HPWL. Fine GP from the polished
  checkpoint ends at 88.69% overflow / 55.18M HPWL.
- This validates capacity accounting as a real improvement, but confirms that
  destination assignment must also be topology-aware at the cluster level.

## 2026-08-24 H321 hierarchical cluster assignment

- Implemented net-connected relative-coordinate clusters across 16/32/64
  coarse levels, exact group overlap trials, exact all-incident-net HPWL
  deltas, cumulative HPWL debt, and live fixed-macro-aware capacity ledgers.
- Found and fixed an H320 accounting error: residual capacity must subtract
  total occupancy, not occupancy with fixed area removed.
- The one-sweep efficiency screen passed, but three cluster-only sweeps ended
  at 63.8652M/96.9369% in 136.37s. The sequential mechanism is rejected.
- Outer-loop conclusion: individual cluster acceptance cannot capture HPWL
  cancellation among jointly moving regions. Pre-register a cooperative
  simultaneous batch with exact global audit rather than tuning this solver.

## 2026-08-24 H322 cooperative batch audit

- Implemented exact global simultaneous cluster audit with tentative capacity
  reservation and prefix backtracking. Core tests remain green.
- Pure density ordering reached 63.8534M/98.1594%; HPWL-efficiency ordering
  reached 63.8462M/97.1887%. Both are runtime-cheap but fail to beat H321.
- Rejected cooperative acceptance as the main missing mechanism. Pivot cluster
  construction from arbitrary BFS to boundary-aware heavy-edge growth.

## 2026-08-24 H323 heavy-edge cluster growth

- Heavy-edge priority growth reached 63.8060M/97.3314% in 3.81s and failed the
  H321 quality screen. Local cluster translation is closed after H321-H323.
- Pivot to H284's globally capacity-effective recursive partitioning. The
  next isolated target is HPWL-aware two-dimensional leaf assignment under
  exact fixed-macro residual capacity.

## 2026-08-24 H324 HPWL-guided bisection leaves

- Added exact all-incident-net HPWL scoring for every leaf anchor and a global
  rectangle-to-bin residual ledger initialized from fixed macro capacity.
- On a freshly regenerated H324 HPWL seed, ordinary mapping reached
  207.101M/1.02232%; guided mapping reached 201.354M/0.295336% in 2.45s.
- The 2.78% HPWL gain rejects leaf mapping as the dominant repair. Pre-register
  multilevel heavy-edge coarsening so upper recursive cuts move supernodes,
  not individual cells.

## 2026-08-24 H325 accumulated shared-net atomic coarsening

- Fixed the atomic matcher to aggregate all shared-net pair contributions
  before deterministic mutual-best matching; the former max-only score was a
  coarsening bug.
- Fresh H325 seed, depth 12: 185.968M HPWL, 0.262431% exact overflow, 4.11s,
  109,799 groups, 101,105 merges, 114,563 atomic unit moves.
- The fix gives only 1.42% HPWL improvement over the earlier H325 depth-12
  run and misses the locked 15% screen from H324 (201.354M).  Close H325 and
  preregister stronger supernodes before sweeping.
## 2026-08-24 H326-H339 continuation

- H326 (16-cell atomic groups, degree 24, four rounds) reached 184.580M /
  0.295225% exact overflow in 4.37s; stronger coarsening alone was rejected.
- H327-H330 exact breakpoint recovery progressively reduced this to 150.344M /
  6.99999% (net-block, ten sweeps, 48.30s). Further sweeps plateaued.
- H331-H332 relaxed exact recovery (15% cap) restored topology to 145.228M /
  12.6599%; H333 hierarchical capacity coordinate descent repaired it to
  0.0347% overflow at 151.498M, and H334 short HPWL recovery reached
  144.024M/6.99982%.
- H336 applied exact axis-ranked bisection to strict H257 and reached
  124.497M/0.2755%; H337 recovery reached 109.490M/6.99998% in 44.81s.
  This is the best new structural route but misses the HPWL target.
- H339 coarse 32-bin leaves worsened to 133.073M and is closed. Pivot next to
  partial capacity cuts and density-neutral net-aware exchanges that do not
  repartition already feasible topology.

## 2026-08-24 H340 net-aware equal-shape exchange

- Added a candidate-search option that enumerates same-shape nodes sharing
  low-degree nets before spatial candidates; exact occupancy is unchanged.
- Five sweeps from H257 reached 86.691M/6.99996% in 90.49s, only 0.10%
  better than H257 and not better than H287. Close long-range equal-shape
  exchange as a main route; next mechanism must move area while preserving
  net topology.

## 2026-08-24 H341 HPWL-directed connected-block recovery

- From H257, exact recovery with node-wise density directions disabled and
  only low-degree connected net blocks enabled accepted 113,772 moves and
  reached 86.710M/7.00000% in 28.83s.
- The 0.47% improvement fails the 2% screen. Connected blocks alone cannot
  escape the feasible HPWL basin; area assignment and topology must be solved
  jointly in the next mechanism.

## 2026-08-24 H342 relaxed net-block then capacity repair

- Ten exact net-block sweeps at a 15% cap moved H257 to 84.761M/15.0%.
  Hierarchical coordinate repair returned overflow to 4.269% but raised HPWL
  to 93.240M, worse than H257.
- Sequential cap relaxation is rejected: its capacity repair spends more
  topology budget than the relaxed stage recovers. Joint destination scoring
  is required.

## 2026-08-24 H343-H344 destination and cluster ablations

- Increasing hierarchical destination candidates 8→64 (H343) was unchanged:
  93.241M/4.217% exact overflow. Candidate truncation is not limiting.
- Reducing rigid cluster max size 64→16 (H344) likewise left HPWL at
  93.240M (overflow 3.876%). Rigid translation itself is the bottleneck;
  next mechanism must use elastic within-cluster assignments with exact audit.

## 2026-08-24 H345-H346 elastic within-cluster assignment

- H345's full per-node breakpoint-fraction oracle exceeded four minutes and
  was stopped without output; the candidate complexity is unacceptable.
- Bounded H346 (one level, 16-node clusters, 8 destinations, 500 clusters)
  completed in 6.66s but reached 92.959M/7.064% from the 84.761M/15% seed.
  Close elastic fractions as a primary route; a cheaper global HPWL oracle is
  needed before revisiting.

## 2026-08-24 H347 coarse flow control

- Exact 64x64 flow from H257 accepted 606 anchor moves, reducing overflow
  6.99996%→6.89671% while HPWL rose 87.119M→87.382M in 0.41s.
- Cheap individual-node flow is not topology-aware; next flow experiment must
  use connected-block commodities rather than more scalar tuning.

## 2026-08-24 H348 connected-block coarse-flow commodities

- Added deterministic degree-16/max-8 net commodities with exact group
  overlap and incident-net HPWL audits. One 64x64 pass from H257 accepted 342
  group moves and reached 87.347M/6.90655% in 0.49s.
- This is slightly worse than H347 individual flow and does not improve H257.
  Single-step commodity transport is closed; test multi-scale exact
  displacement next.

## 2026-08-24 H349 multi-scale exact commodity displacement

- Locked and implemented exact rigid commodity scales `{1.0, 0.5, 0.25}`.
- Fresh staged continuation from H257: 85,545 scale trials, 494 accepted
  commodity moves, 87.340467M HPWL, 6.889111% exact overflow, 0.541525s
  end-to-end (0.215905s flow).
- This is only a marginal improvement over H348 and remains above H257's
  HPWL. Scale is not sufficient; pre-register axis-separated directions as
  H350 rather than tuning the scalar flow further.

## 2026-08-24 H350 axis-split exact commodity directions

- Added exact x-only and y-only rigid directions to H349's diagonal/multiscale
  candidates. Fresh H257 continuation: 521 moves, 87.324700M HPWL,
  6.88948% exact overflow, 0.6344s end-to-end.
- The 0.018% HPWL improvement over H349 is not material and does not alter
  the conclusion that coarse flow directions cause topology debt. Pre-register
  H351 with a strict exact incident-net HPWL guard.

## 2026-08-24 H351 exact HPWL-guarded axis flow

- Added a hard exact non-positive incident-net HPWL guard. Only 12 moves were
  accepted; the result was 87.1195M HPWL / 6.99966% overflow in 0.589148s.
- The guard preserves H257 topology but provides no meaningful density gain.
  Local flow therefore needs positive HPWL debt; pivot to joint global
  commodity destination assignment with exact batch cancellation screening.

## 2026-08-24 H352 jointly audited commodity destination assignment

- Implemented a global disjoint commodity proposal batch with exact overlap
  and temporary-coordinate global max-min HPWL audits. The run made 26,235
  scale/direction trials and 186 prefix trials.
- The monotone prefix policy committed 12 commodities: 87.128095M HPWL,
  6.994967% exact overflow, 0.826598s. It does not improve H350/H351.
- The policy likely rejects useful cancellation prefixes; pre-register H353
  with a bounded non-monotone reservoir and deferred exact prefix selection.

## 2026-08-24 H353-H354 deferred joint reservoir ablations

- H353's deferred non-monotone prefixes produced 12 groups and
  87.120209M/6.997823% in 0.611627s; the top-8 edge reservoir was too narrow.
- H354 increased the edge reservoir to 64 but still selected 12 groups,
  87.121849M/6.996061% in 0.619557s. The limiting filter is group area versus
  individual flow-edge area, not proposal count. Remove that filter next.

## 2026-08-24 H355 joint assignment without edge-area filter

- Removed only the per-edge group-area filter. The run still selected 13
  groups (12 effective prefix groups), 87.121849M/6.996061% exact metrics in
  0.608213s.
- Edge fragmentation is not the cause; repeated top bids dominate. Test a
  much wider per-edge proposal pool next, then pivot if the batch remains tiny.

## 2026-08-24 H356 very-wide joint bid pool

- Retained 1,024 exact bids per flow edge; 64 global prefixes were audited but
  the selected prefix remained 12 groups, 87.121849M/6.996061% in 0.616319s.
- Candidate coverage is ruled out. Replace ratio ranking with an exact global
  HPWL trust-region budget next.

## 2026-08-24 H357 exact HPWL trust-region joint batch

- The 1% exact global HPWL budget selected 64 groups, but ended at
  87.1403M/6.96741% in 0.594443s, worse overflow than H350's 6.88948%.
- Ratio ranking was conservative but larger joint batches do not produce a
  useful topology-cancelling assignment. Close post-H257 joint flow and pivot
  to raw-to-feasible initialization/continuation.

## 2026-08-24 H358 early-share seed capacity flow

- Applied one exact fixed-macro-aware 64x64 flow to the staged H304 seed. All
  110,592 exact-anchor moves were rejected; HPWL/overflow stayed
  66.8702M/83.8347% in 0.688615s.
- The low-HPWL share seed is in a flat exact capacity basin. Local flow is
  insufficient; coordinated nonlocal raw-to-feasible assignment is required.

## 2026-08-24 H359 nonlocal commodity jump from early-share seed

- Direct full source-to-sink commodity jumps accepted 19,733 moves and reduced
  exact overflow 83.8347% -> 66.3385%, but HPWL rose 66.8702M -> 244.359M in
  28.6282s.
- Nonlocal flow crosses the flat basin but destroys topology; close direct
  jumps and retain only the need for coordinated relative-net transport.

## 2026-08-24 H360 larger net-connected nonlocal groups

- Degree-64/max-64 direct commodity jumps accepted 18,249 moves, reaching
  256.666M HPWL / 66.9037% overflow in 23.372s.
- Larger rigid groups are worse than H359; source-to-sink jumps cross too many
  HPWL breakpoints. Close this family and pivot to recursive cuts or a new raw
  topology seed.

## 2026-08-24 H361-H362 recursive share-seed assignment and recovery

- H361 position-seeded atomic bisection reached 155.243M HPWL / 0.250171%
  overflow in 5.761s from H304's 66.8702M/83.8347% seed.
- H362 ten exact breakpoint/compact/net-block recovery sweeps reached
  124.767M/6.99988% in 39.3389s. The 19.6% recovery is promising; extend the
  same exact slack continuation before abandoning the route.

## 2026-08-24 H363-H367 continuation

- H363 extended H362 by twenty exact sweeps: 124.767M -> 123.236M,
  7.00000% overflow, 77.08s. Diminishing returns close the concentrated-share
  branch.
- H364 leaf-map ablation showed curve-rank mapping is topology-destructive
  (249.705M); axis+curve was 219.281M. Keep H361's axis map only.
- H365 extended the feasible H337 recovery with net-block density/contraction:
  109.490M -> 104.656M, 7.00000%, 129.94s. Quality improved but full staged
  runtime exceeds the strict a1 envelope.
- H366 degree-16 bisection reproduced H336 (124.497M/0.2755%); the
  position-seeded path bypasses degree ordering, so this knob is inert.
- H367 trimmed to five recovery sweeps: 105.295M/6.99976% in 31.99s. The
  staged H336+H337+H367 route is runtime-valid and is the current strict a1
  record, but remains 43.8% above paper HPWL. Next work must target
  topology-preserving upper recursive cuts rather than more local sweeps.
- H368 tested the H340 exchange seed before capacity bisection and recovery:
  106.652M/7.00000% in 30.53s, worse than H367. Close exchange-seed transfer;
  upper cuts remain the bottleneck.
- H369 audited exact alpha interpolation between H257 and H336. Alpha 0.25
  recovered 93.645M/31.06% to 87.360M/20.60% after twenty sweeps (130.60s),
  then stalled; degree-32/max-16 blocks did not reach 7%. Interpolation is
  closed as a capacity bridge, reinforcing the need for topology-preserving
  hierarchical cuts.
- H370 tested relaxed recursive cut balance (15%, plus 49% diagnostic):
  106.665M/7.00000% after five exact recovery sweeps, worse than H367. The
  balance tolerance is not a useful lever; upper-cut assignment itself must
  change.
- H371 added `bisection-surplus-only`, moving only cut-boundary surplus needed
  by exact side capacities. The checkpoint was 114.146M/2.16494%; five exact
  recovery sweeps reached 102.008M/7.00000% in 30.62s. Ten and twenty sweeps
  reached 101.162M and 100.916M, respectively.
- H372 disabled the axis pre-map while retaining surplus-only cuts. It reached
  93.7055M/2.30636% at the checkpoint and 86.4503M/6.99998% after five exact
  sweeps (30.59s); ten and twenty sweeps gave 86.1896M and 86.1566M. This is
  the current strict a1 quality record and validates topology-preserving
  hierarchical capacity assignment; the remaining gap is now 17.7%.
- H373 raw a2 transfer failed as a seed test: exact HPWL-only 500 steps gave
  50.7285M/95.8648%, and surplus-only cuts produced 191.591M/1.36475%.
  Defer a2-a4 until a feasible topology-preserving seed is available.
- H374 polished H372 five-sweep output with exact net-aware equal-shape swaps:
  86.4503M -> 86.2267M at unchanged overflow.
- H375 repeated the polish after ten H372 sweeps: 86.1896M -> 85.9993M,
  6.99999% overflow, within the staged a1 runtime budget. Exchange is retained
  as a final polish but is not sufficient alone to close the HPWL gap.
- H376 exact relax/retighten test: cap 0.15 lowered HPWL to 83.8203M, but
  five cap-0.07 sweeps could not leave 14.87% overflow. H377 surplus-only
  recut after relaxation reached 86.1054M/6.99992%, worse than H375; close
  this branch and retain the feasible H372/H375 route.
- H378 raw a2 exact joint GP (1000 iterations, density scale 5) remained at
  97.7351% overflow with 120.171M HPWL in 22.19s; lambda grew to 2.22e8.
  The exact-overlap active set is flat from the raw center seed, so defer
  a2-a4 transfer until a capacity continuation is designed.
- H379 exact recovery from the raw a2 HPWL seed at a permissive 99% cap still
  left overflow at 96.1504%; cap relaxation alone cannot activate density.
- H380 sigma-0.1 raw initialization followed by 500 exact HPWL-only steps
  ended at 59.835M/98.149% overflow. Spread initialization without a
  capacity-aware stage is closed.
- H381 affine-expanded the raw a2 HPWL seed: scale2 produced
  277.768M/1.654% after surplus cuts and scale3 223.655M/1.617%. Expanded
  coordinates still cause topology-destructive capacity assignment; close
  this seed family.
- H382 atomic surplus cuts on the raw a2 seed reached 253.393M/0.302%, worse
  than non-atomic H373. Net coarsening is insufficient after HPWL-only
  collapse; close this variant.

## 2026-08-24 H389-H390 continuation controls

- H389 elastic heavy-edge clusters reached 119.766M / 95.3168% in 159.5s;
  extra local proposals did not solve activation and increased runtime.
- H390 depth-2 exact cuts reached 1,822.22M / 0.153% but took 259.4s due to
  global exact leaf scans. Shallow cuts are rejected as a production route.

## 2026-08-24 H391-H395 bounded transport ablations

- H391 adjacent-bin cluster destinations preserved HPWL (91.395M) but only
  reached 98.9234% overflow in 28.7s.
- H392 repeated the local operator for 30 sweeps: 91.604M / 98.5347% in
  125.3s. Long local continuation is closed.
- H393 exact node breakpoint moves accepted only 140 moves and reached
  79.980M / 99.9699%; node locality cannot activate the basin.
- H394 larger center jitter activated exact gradients but produced 289.171M /
  97.7882%; random spread is not a capacity-aware initializer.
- H395 radius-one rigid-group transport accepted zero exact-feasible moves,
  confirming that a coordinated multi-bin bridge is required.
- H396 radius-two cluster destinations reached 119.767M / 95.6187% overflow
  in 56.3s; finite-radius local transport is insufficient.
- H397 surplus cuts on the H386 checkpoint jumped to 1,491.87M HPWL before
  recovery (1,312.57M / 7% after), confirming the coordinate and cut basins
  are disconnected.
- H398 allowed a 2-point exact nonmonotone overflow band but reached only
  119.762M / 98.6559% in 130.7s.
- H399 widened the band to 20 points at 128 bins and still reached 119.765M /
  97.7337% in 34.9s. Overflow allowance alone does not cross the raw basin.
- H400 used 64 exact bins and a 50-point continuation band; two sweeps reduced
  overflow to 86.9828% with 223.188M HPWL in 9.8s.
- H401 extended to four sweeps and reached 83.7046% overflow / 239.479M HPWL
  in 22.1s. Coarse staged assignment is the first partial bridge, but needs a
  topology-preserving fine refinement stage.
- H402 rebuilt H400 at 256 exact bins: the 64-bin 86.98% checkpoint became
  96.34% overflow, and four fine sweeps reached 267.829M / 83.9105% in 26.2s.
  Resolution transfer requires explicit occupancy-aware prolongation.

## 2026-08-24 H403-H404 ablation and active-set escape

- H403 replayed the H375 recovery prefix from the identical H372 ten-sweep
  checkpoint. Node-only recovery ended at 86.1786M; adding breakpoint oracle
  gave 86.1764M; adding compact directions gave 86.1721M; breakpoint plus
  net-block density directions gave the best isolated recovery, 86.1683M.
  The full contraction variant ended at 86.1704M, so contraction was not
  beneficial in this replay. Equal-shape exchange polished an A5 checkpoint
  to 85.9826M, but the combined prefix exceeds the end-to-end runtime gate.
- H404 applied exact breakpoint/net-block escape directly to raw a1 center
  initialization. It reduced HPWL 53.2216M -> 47.1525M in 59.6s, while exact
  overflow remained 99.6766%. Breakpoint escape alone optimizes HPWL inside
  the raw basin but does not establish capacity; it must be interleaved with
  staged capacity targets.

## 2026-08-24 H405-H410 isolated mechanism tests

- H405 enabled only the exact breakpoint oracle from the same raw a1 seed:
  53.2191M -> 48.0784M HPWL, 99.9679% -> 99.7227% overflow, 266,103 moves,
  8.67s.
- H406 enabled only a 0.25-bin boundary-active strip: 48.2662M HPWL and
  99.6457% overflow, 290,329 moves, 6.80s. Extra boundary candidates did
  not improve HPWL over H405.
- H407 enabled only zero-gradient escape: 50.3151M / 99.6406% in 3.08s.
  It is cheaper but weaker than the regular breakpoint oracle.
- H408 enabled only net-aware shared density direction for 100 GP iterations:
  43.7691M / 97.0140% in 3.26s. This is the strongest raw HPWL aggregator,
  but exact density remains essentially infeasible.
- H409 enabled only node-wise exact finite directional coordinate moves:
  53.2191M / 99.4508% in 6.89s. It lowers overlap slightly without changing
  HPWL.
- H410 enabled only deterministic exact noise on zero-gradient nodes:
  50.3109M / 99.6430% in 3.06s, essentially matching H407.

The six results are causal single-component measurements, not a combined
mechanism claim. They rule out dead-zone correction alone as a replacement for
hierarchical capacity assignment and motivate testing breakpoint/net-share as
post-assignment recovery operators only.

## 2026-08-24 H383-H388 raw a2 capacity-aware initialization

- H383 applied surplus-only recursive cuts directly to raw a2.  The exact cut
  checkpoint was 2,439.03M HPWL / 1.156% overflow and five recovery sweeps
  ended at 2,297.98M / 7.000% (14.2s).  Direct raw cuts are topology
  destructive and closed.
- H384/H385 tested a short ten-step exact HPWL warmup followed by local leaf
  assignment.  A radius-two locality guard reduced the checkpoint overflow
  only to 66.955% and recovery to 28.046%, with 1,754.37M HPWL.  Locality
  prevents the catastrophic move but cannot bridge enough capacity.
- H386 introduced exact hierarchical cluster coordinate descent on raw a2.
  Five sweeps reached 95.814M / 98.5304%, and recovery reached 83.1404M /
  98.5277% in 21.5s.  This is the first raw mechanism to activate exact
  overlap while keeping HPWL in the 100M range, but it is far from feasible.
- H388 extended the same operator to 50 requested sweeps (256 bins).  It
  stopped after 15 productive sweeps at 159.691M / 94.3382% in 59.8s;
  additional local proposals are not a solution.  Retain the operator only
  as an interleaved warmup component and pivot to a stronger capacity-aware
  bridge rather than repeating local sweeps.
