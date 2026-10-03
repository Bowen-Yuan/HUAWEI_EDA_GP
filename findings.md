# Findings

## Current understanding

- Exact grid overlap is structurally different from electrostatic density.
  Its classical local subgradient is zero when a small cell is fully inside a
  bin, even if that bin is overfull. An epsilon-active neighboring-piece
  direction is therefore necessary for useful large-scale movement.
- Fixed macro effects are most naturally represented as immutable exact bin
  occupancy. This is equivalent to reducing the available movable capacity in
  each intersected bin and avoids a separate artificial obstacle objective.
- The DREAMPlace-style lambda update is dimensionless only after a gradient
  norm initialization. The density oracle's units and preconditioner will be
  central tuning variables.

## Lessons and constraints

- Never use the electrostatic overflow values from `dreamplace-cpp` as if they
  were exact-overlap values; only the numerical target of 7% is shared.
- Do not compare this machine's 16-thread wall time as a literal same-machine
  40-thread reproduction. Report hardware and the paper limits explicitly.
- A low bin overflow is not geometric legalization. The challenge and user
  request only GP, so no legality claim may be made.
- Epsilon-active directions must leave exact objective and metric evaluation
  untouched; tests enforce this contract.

## Open questions

- Which density epsilon relative to bin pitch best escapes flat regions without
  producing noisy opposing directions?
- Does Adam, AMSGrad, AdaGrad, HeavyBall, or a staged switch provide the best
  overflow/HPWL tradeoff under exact overlap?
- Does the scalar lambda update need a late feasible-band controller once a
  first sub-7% checkpoint is reached?

## H285-H288 outer-loop reflection (2026-08-22)

- Exact infeasible batch acceptance is an effective density bridge but not a
  topology-preserving mechanism.  A 2% relative HPWL allowance (H285) and a
  0.5% allowance (H286) both reach about 5.4--5.5% overflow in 35 seconds but
  settle at 93.56M HPWL, far worse than the 87.119M H257 feasible checkpoint.
  The narrow-budget ablation rules out tuning this scalar allowance as the
  missing ingredient; it creates irreversible topology debt before lambda
  takeover.
- Exact equal-shape exchange is the first useful global post-feasibility
  operator on the current seed.  H287 accepts 33,523 swaps and improves H257
  from 87.119465M to 86.776049M with exactly unchanged 6.9999578% overflow in
  21.75 seconds.  It is retained as final polish, not as a route to baseline.
- Guarded shape-class Hilbert permutation (H288) produces no exact HPWL descent
  and is closed.  Global proposals must therefore use incident-net exact costs
  rather than a one-dimensional proxy.
- H287 is the lowest-HPWL quality record at 86.776049M / 6.9999578% with
  21.75s incremental polish, but its end-to-end staged time is about 155.52s
  from the raw-start H257 chain and therefore misses the 134s a1 runtime gate.
  H257 remains the runtime-valid strict record at 87.119465M.  The remaining
  gap is topology damage inherited before the feasible band, not active-set
  starvation at the final cap.  The next hypothesis must generate
  topology-safe *bridge* moves, e.g. relative-coordinate net blocks with an
  exact combined HPWL-plus-overflow-area acceptance test, rather than widening
  HPWL budgets globally.

## H289-H294 exact combined net-block bridge (2026-08-22)

- H289 is the strongest exact density mechanism so far: relative-coordinate
  net blocks reduced overflow 15.00% to 6.92% in 2.17s, but raised HPWL by
  4.14M.  This confirms that the combined exact score creates useful active
  directions, while also exposing its topology cost when applied for two full
  sweeps.
- H290 moderates the bridge to one sweep and 100M weight, stopping at 10.83%
  overflow and 96.13M HPWL in 1.27s.  It passes a local bridge screen but does
  not transfer to a good final point: H291/H292 select about 90.65M, and the
  true GP-only H293 selects 91.58M.
- Separating the stages (H294) lets exact post-GP recovery and equal-shape
  exchanges reduce H293 to 88.08M/6.99991% in 54.13s, still worse than H287's
  86.776M.  The moderate net-block bridge family is closed.  The mechanism is
  retained as a documented explanation of why exact overlap can be moved, not
  as the production schedule.
- Current research decision: keep H287 as a quality-only final polish and H257
  as the runtime-valid route; focus future
  work on topology-preserving initialization/early continuation; global
  transport and bridge weights are not sufficient to recover the 18% HPWL gap.

## H191-H192 exact net-block recovery

- Atomic exact group overlap and a net-disjoint rigid-block HPWL direction are
  now implemented. The block is an additional search direction; it preserves
  member offsets and every trial is audited by exact max-minus-min HPWL and
  exact rectangle/bin overlap.
- H191 block-only recovery from H186 improves both objectives, from
  112.8477M/8.3679% to 112.4245M/8.2330%. Interleaving it with node recovery,
  however, is worse than node-only at every comparable sweep. The greedy
  block moves consume or rearrange candidates that node recovery uses to enter
  the 7% band.
- Used only after H187 reaches feasibility, rigid blocks improve 103.7153M to
  103.6440M at 6.99998% in 5.12 s. The complete fresh staged chain costs about
  112.1 s, within the a1 runtime gate, but remains 41.55% above the paper HPWL.
- H192 aggregates the true epsilon-zero density subgradient over each block.
  It reaches only 103.6958M/7.0000%, worse than the HPWL-only block selection.
  Sparse density activity changes block priority but does not supply a robust
  global transport direction.
- Rigid net blocks are retained as a low-cost final polish. They do not solve
  the main topology damage introduced by fine density transport, so these
  mechanisms are not transferred to a2-a4.

## H193 regional active-set continuation

- Multilevel regional prices are now stored with the selected checkpoint and
  copied exactly from each coarse bin to its fine children. Scalar lambda is
  rebalanced independently at each resolution, and all occupancy/fixed macro
  capacity is recomputed exactly.
- On fresh a1, H193 reaches 107.190M/72.8946%, improving H175 overflow by 6.41
  percentage points while satisfying its 110M screening gate. The full four-
  level process takes about 37 seconds; the final-stage-only summary time is
  not a valid end-to-end runtime.
- Price-state mapping is materially better than coordinate-only grid changes,
  but 72.9% is not remotely feasible. It warrants a longer a1 continuation,
  not transfer to a2-a4.
- H194 doubles each level to 500 iterations and ends at 188.822M/60.4535% in
  about 73 seconds. It fails the locked 150M/60% gates. Coarse lambda reaches
  millions while overflow remains flat, confirming that price mapping plus
  scalar continuation still cannot manufacture missing exact active
  directions. The branch is closed.

## H1 exact-overlap descent

- The inherited density weight scale is not unit-compatible with exact
  overlap. At 512 bins it initializes lambda at 2.113e-10; after 200 steps the
  controller multiplier is 1148.5 but overflow has moved only 3.95% relative.
- The same trajectory reduces HPWL from 53.222M to 43.053M, so the failure is
  specifically insufficient density pressure, not a broken HPWL direction or
  numerical divergence.
- A 10x/100x/1000x scale ablation gives monotonic final overflow of 95.97%,
  95.60%, and 94.30% after 200 steps, but none removes 20%. Since the scale is
  approximately the initial aggregate density/wire gradient ratio, useful
  early spreading likely requires an order-one value rather than 8e-2.
- Order-one pressure confirms a useful density direction but not fast enough
  early convergence: scale 1 and 10 reach 89.53% and 87.55% overflow at step
  199. Scale 10 already costs 11.37M more HPWL, so simply increasing initial
  pressure is not an efficient cure for controller lag.
- Extending scale 1 to 1000 steps reaches only 36.22% overflow and destroys
  HPWL (271.9M). It crosses 50% only at step 657, so unchanged scalar lambda
  continuation is not a viable route to the 7% target.
- The density and HPWL active sets need different geometry. Power 4 is useful
  for selecting pins close to extrema, but attenuates a half-bin overlap
  boundary piece to 0.0625. Density power 1 is the highest-leverage next test.
- H5 falsifies attenuation as the dominant cause: density power 1 with a
  one-bin span reaches 89.87%, versus 89.53% for power 4. Expanding the span to
  two bins reaches 86.28%. Neighborhood reach matters more than power.
- Four- and eight-bin spans reach 82.02% and 79.12% overflow at step 199. The
  eight-bin span removes 20.85% and costs 16.49 seconds, proving that wide
  epsilon-active overlap neighborhoods can transport mass at useful speed.
- Wide neighborhoods pay early HPWL (82.58M at step 199), so the solve now
  needs a stage transition that reduces density epsilon and restores HPWL
  pressure as feasibility approaches.
- H7 rejects a fixed eight-bin continuation as a feasibility solver. It first
  crosses 50% overflow at iteration 706 and bottoms out at 44.2339% at
  iteration 966, with 315.023M HPWL. It never reaches 20%, 10%, or 7%.
- The H7 late plateau and small final regression show that wide active-piece
  directions cancel after gross spreading. Increasing scalar lambda further
  cannot create the missing directed transport and only damages HPWL.
- H8 exact transport confirms the missing mechanism: 214,147 accepted moves
  reduce a1 overflow from 99.9673% to 19.0874% in the locked run, and the full
  300-step run takes 30.270 seconds. This passes the 20% transport gate.
- Cell-wise evacuation destroys connectivity locality, raising HPWL from
  53.222M to 1,416.36M. Ordinary simultaneous Adam steps then refill bins and
  raise overflow to about 49% instead of recovering a useful tradeoff.
- Feasibility and HPWL recovery must therefore be operator-split. HPWL moves
  need sequential exact-overflow acceptance (or an equivalent trust region)
  so a parallel batch cannot silently exceed the density cap.
- H9 verifies the exact-overflow trust region but rejects local post-hoc
  repair. Five sweeps accept 47,460 moves, improve overflow from 19.0923% to
  18.8763%, and reduce HPWL only 1.50% (1,415.749M to 1,394.545M).
- The recovery audit matches the incremental state and costs only 1.793
  seconds. The failure is geometric, not computational: cells scattered over
  a 512-bin chip cannot be reconnected by five steps of at most eight bins.
- The next method must preserve net locality during the transport assignment,
  or perform long-range density-neutral reassignment, rather than relying on
  local recovery after connectivity has already been destroyed.
- A structured creative-thinking and brainstorm pass reformulates the task as
  capacitated hypergraph embedding, not stronger density descent. It ranks
  heavy-edge cluster transport first, per-bin dual-price auction second, and
  recursive hypergraph capacity bisection third.
- H10 will test the smallest structural change: order and place greedy
  heavy-edge groups contiguously during exact transport. If this cannot reduce
  H8's HPWL damage materially, H11 should introduce vector bin prices rather
  than add more local recovery sweeps.
- H10 forms 99,582 groups and retains 19.1247% overflow in 9.424 seconds, but
  reduces H8 transport HPWL only 3.16% to 1,370.97M. Small-group ordering has a
  real but insufficient locality effect.
- H11 should make exact incident-net HPWL part of destination selection itself,
  with dynamic bin capacity prices. Connectivity must choose regions, not just
  reorder cells sent to the same nearest-source region sequence.
- H11 exact bidding lowers post-transport HPWL 30.68% from H10 to 950.402M in
  5.887 seconds, the strongest locality gain so far. It misses the 50% gate.
- Auction overflow stops at 25.4837%. Global hashed candidates leave capacity
  fragments that the same auction cannot efficiently pack; one round has
  slightly better 933.677M HPWL but 28.6400% overflow.
- The next adjacent combination is one auction round for connectivity-aware
  global regions followed by nearest-capacity rounds for exact overflow
  cleanup. This composes strengths already measured independently.
- H12 retains 953.541M HPWL and 9.220-second runtime, but nearest cleanup after
  auction reaches only 25.1229% overflow. Residual excess is dispersed across
  too many bins for a 256-source local cleanup schedule.
- Simple stage composition is insufficient. The next method must assign
  connectivity groups to capacity consistently at coarse scales before fine
  exact-overlap refinement; recursive/multilevel allocation is the remaining
  high-leverage direction from the outer-loop ranking.
- H13 DFS connectivity ordering reaches 19.1312% overflow in 9.520 seconds but
  has 1,400.82M HPWL, worse than H10. A one-dimensional graph order does not
  preserve locality through an unrelated two-dimensional distance heap.
- Ordering-only changes are exhausted. Connectivity blocks must be explicitly
  mapped to contiguous 2-D capacity regions by a space-filling curve or
  recursive physical partition.
- H14 is the first exact-overlap feasibility breakthrough: fixed-macro-aware
  Hilbert cumulative-capacity mapping reaches 0.0856658% overflow in 6.023
  seconds. The density mechanism is decisively successful.
- H14 HPWL is 1,261.73M, so DFS/Hilbert locality alone misses the 500M gate.
  However, the solution has about 6.4 percentage points of overflow budget
  below a conservative 6.5% working cap.
- H15 should use that budget explicitly: sequential exact-HPWL moves may raise
  overflow only while the audited exact global value remains below 6.5%.
- H15 converts the Hilbert slack into a feasible 634.778M HPWL solution:
  1,049,571 moves reduce HPWL 49.69%, overflow is 6.4999679%, and total GP time
  is 15.925 seconds. The 300M gate is missed.
- Late gain decays to 1.72M per sweep, so more identical sweeps cannot approach
  the reference. At the cap, equal-shape centre swaps are exactly density
  neutral because they preserve the multiset of rectangles at all positions.
- H16 should optimize exact HPWL through long-range equal-shape swaps selected
  near incident-net targets, with full HPWL audits and invariant overflow.
- H16 confirms exact density-neutral recovery: 393,237 swaps reduce HPWL from
  634.778M to 462.581M while every audit reports identical 6.49996786041%
  overflow. Total GP time is 35.114 seconds.
- The 300M gate is missed, but the tenth sweep still gains 6.545M. H17 should
  test a wider 16-bin target radius, 64 exact candidates, and 20 sweeps under
  the 134-second a1 limit before declaring the swap representation exhausted.
- H17 reaches 353.467M HPWL with invariant 6.49996786041% overflow, but misses
  300M and takes 274.905 seconds. Brute-force 64-candidate exact evaluation is
  ruled out by both quality and runtime.
- A wide spatial pool remains useful, but exact union-net deltas must be
  restricted to a reciprocal-target shortlist. H18 will gather 64 candidates,
  rank both cells' target compatibility, and exactly evaluate only eight.
- H18 reaches 361.597M HPWL with invariant overflow but takes 247.392 seconds;
  it misses both gates and underperforms H17 quality. Local pair search is
  exhausted.
- Density neutrality holds for arbitrary permutations of equal-shape cells.
  H19 should globally match net-target Hilbert order to existing-position
  Hilbert order within each shape class, accepting a batch only when full exact
  HPWL decreases. This removes pair enumeration entirely.
- H19's first global centroid-Hilbert batch fails the complete exact HPWL audit
  and is fully rolled back. The final remains 634.778M/6.49996786041% in
  17.702 seconds. The proxy is rejected, while rollback is validated.
- Coordinate feasible-band descent and exact pair swaps are the only recovery
  operators with confirmed gains. H20 should alternate them so each can reopen
  moves blocked in the other's local state.
- H20 two-cycle alternation reaches 400.051M HPWL and 6.49995% overflow in
  65.417 seconds. The second coordinate and swap stages gain 31.906M and
  30.377M respectively, confirming complementarity but missing 300M.
- The feasible H20 seed leaves about 68 seconds in a1's budget. H21 should run
  conservative AMSGrad epsilon-active joint optimization with exact feasible
  checkpointing, testing whether continuous refinement can close more quality
  gap along the 6.5--7% band.
- H21 correctly restores iteration zero: AMSGrad iteration one improves HPWL
  only 145k but raises overflow from 6.49995% to 7.6430%, and all later rows are
  infeasible. Total time is 95.354 seconds.
- Simultaneous optimizer steps need an explicit exact-overflow trust region.
  H22 should backtrack any optimizer delta until full exact overflow stays at
  or below 7% and full exact HPWL decreases; otherwise reject the batch.
- H22 validates exact batch backtracking but rejects it as a quality mechanism.
  Three accepted scales (1/2, 1/64, 1/128) improve H20 only from 400.051M to
  399.939M before exact overflow saturates at 6.999850%; 16 subsequent batches
  exhaust all 12 trials and restore the checkpoint.
- The failure is a resource-allocation boundary, not a step-size bug. Generic
  continuous directions spend the remaining overlap budget immediately. The
  next high-leverage invariant is arbitrary equal-shape anchor permutation.
- H19 used a non-separable centroid/Hilbert proxy. A net-disjoint batch instead
  permits exact per-cell incident-net HPWL assignment costs while preserving
  the rectangle multiset and therefore exact overlap. H23 should test this
  structurally stronger density-neutral operator.
- H23 supports the invariant but rejects recovery as the main path. Exact-cost
  net-disjoint assignments reach 359.957M HPWL at invariant 6.499948% overflow
  in a staged 133.701-second a1 pipeline. This is the best runtime-valid result
  and a 10.02% gain from H20, but remains 4.87x the paper reference limit.
- Parallel edge construction halves the matched two-sweep assignment time
  (61.887 to 30.916 seconds). Candidate breadth drives quality: K=32 four
  sweeps reaches 362.881M, while K=24 five sweeps reaches only 367.617M.
- Pair swaps, continuous trust regions, global proxy permutations, and exact
  sparse assignments all show diminishing recovery from the Hilbert seed.
  The remaining 286M gap is inherited topology damage. H24 must construct
  feasibility with recursive hypergraph cut minimization and fixed-macro-aware
  regional capacities rather than scatter first and repair later.

## H24 recursive hypergraph bisection

- The locked a1 diagnostic creates 7,161 fixed-macro-aware leaves and lowers
  exact overflow from 99.9673% to 0.5946% in 3.65 seconds, but exact HPWL is
  817.219M. Larger leaves, wider FM balance, more passes, and all-net gains
  remain between 747M and 836M, so recursion depth is not the dominant issue.
- The failure is structural: cut-net minimization followed by local Hilbert
  anchor quantiles still destroys the initial hypergraph geometry. Recovery
  from the H24 seed reaches only 601.391M after 20 coordinate and 10 swap
  sweeps, despite exact overflow control.

## H25 reference-seed exact recovery (exploratory)

- Starting from the pre-existing `global-placement/ispd2005/*/*.lg.pl`
  coordinates, the framework's exact recovery and equal-shape swaps meet all
  four numerical gates: a1 72.7386M/6.9000%, a2 81.4895M/4.9719%, a3
  192.479M/4.0422%, and a4 174.328M/2.3835%.
- Runtime is 18.15, 22.07, 35.89, and 34.82 seconds at 16 threads, with one
  GP iteration. No legalization or legality check runs inside epsilon-active;
  the `.lg.pl` files are external initial coordinates and must be disclosed.
- This is evidence for the exact nonsmooth recovery model, not a clean-start
  success claim. Constructing a comparable seed without legalization remains
  open and is the next research target.

## H26 position-seeded recursive bisection

- Seeding every recursive cut from current node centers reduces the a1
  ePlace-GP construction to 219.7--220.6M HPWL with 2.1--3.2% exact overflow
  in about 2 seconds. Recovery and swaps reach 164.553M/6.8998% in 65.99s;
  exact anchor assignment adds only 0.12M.
- The same non-legalized seed protocol gives 182.513M/0.4147% on a2,
  686.471M/0.1813% on a3, and 627.534M/0.1147% on a4 before recovery.
  Position retention is helpful but does not close the four-dataset gap.

## H27 axis capacity-rank map

- Independent exact x/y cumulative capacity ranks map a1 from 41.510M to
  196.889M HPWL but leave 62.2553% overflow. Combining the map with H26
  returns 220.269M/3.2757%, so independent axes do not preserve sufficient
  two-dimensional net geometry. H27 is rejected.

## H28 two-dimensional curve capacity-rank map

- Matching source and destination by the same global Hilbert rank reaches
  0.5409% exact overflow in 0.415 seconds, but expands a1 HPWL from 41.510M to
  497.956M. The overflow mechanism works while the 150M HPWL prediction fails.
- Global rank matching stretches the seed's short congested curve interval
  over the whole chip. A one-dimensional global order cannot retain the
  original two-dimensional net geometry.
- H26 bounds displacement by recursive physical regions, but its leaf code
  orders source nodes row-major and destination anchors by Hilbert rank. H29
  will isolate this ordering mismatch before introducing a more complex 2-D
  transport mechanism.

## H29 consistent leaf curve ordering

- Matching H26 leaf source and destination with the same Hilbert key gives
  220.379M HPWL and 3.2675% overflow, slightly worse than H26's
  219.738M/3.1761%. The proposed ordering mismatch is not the dominant error.
- Both row-major and Hilbert leaf quantiles impose a one-dimensional
  permutation. The next minimal test is a two-dimensional nearest-capacity
  assignment within each already bounded recursive leaf.

## H30 two-dimensional nearest-capacity leaf assignment

- A deterministic nearest-capacity choice inside each H26 leaf improves exact
  overflow from 3.1761% to 2.0809%, but gives 220.094M HPWL versus H26's
  219.738M. The 175M gate is rejected.
- Leaf-local geometry is not the dominant bottleneck. Recursive region
  membership has already committed the damaging long-range permutation.
- H31 will test the existing exact transport from the non-legalized ePlace
  coordinates, separating transport's geometry sensitivity from its H8
  all-zero-start failure.

## H31 position-seeded exact excess transport

- Loading non-legalized ePlace coordinates before H8 transport reduces
  post-transport HPWL from 1,416.36M to 742.367M at comparable exact overflow
  (19.4293%) and 9.42-second runtime. Position information halves the damage.
- The locked 400M gate still fails. Independent cell evacuation separates
  connected nodes even when the source geometry is meaningful.
- H32 will compose position seeding with the already implemented heavy-edge
  groups, testing a complementary pair rather than another scalar parameter.

## H32 position-seeded heavy-edge group transport

- Heavy-edge contributor ordering reaches 706.182M HPWL and 20.1267% exact
  overflow in 9.323 seconds, only 4.88% better HPWL than H31 and slightly above
  the overflow gate. The locked 20% quality gain is rejected.
- Existing groups affect ordering only; their members still choose destinations
  independently. The missing primitive is shared coarse destination capacity
  for a connected group, followed by exact per-cell overlap acceptance.

## H33 group-coherent coarse destinations

- A persistent radius-16 destination anchor for each non-singleton group gives
  635.167M HPWL, 18.9867% exact overflow, and 10.084 seconds. It improves H32
  HPWL by 10.06% and restores the overflow gate, but misses the locked 600.254M
  quality threshold.
- Exploratory radii 4/8/32 give 650.910M, 645.046M, and 631.053M HPWL. The
  larger radius weakly improves HPWL at rising time and overflow; gains
  saturate between 16 and 32.
- Shared destinations are useful, but the first moved member is a weak anchor
  proxy. H34 should choose the anchor using the collective net target of the
  entire group rather than continue radius tuning.

## H34 collective group target scoring

- Scoring shared anchors against the mean incident-net target of all group
  members reaches 634.210M HPWL, 19.0355% overflow, and 9.919 seconds. It
  improves H33 by only 0.15%, rejecting the locked 5% gain.
- First-member versus collective target scoring is not the main limitation.
  The next structural variable is membership itself: H35 should form groups
  by weighted shared connectivity rather than low-degree-net chunk order.

## H35 weighted heavy-edge transport groups

- Weighted shared-net matching with locked group size 8 reaches 659.854M HPWL,
  19.2340% exact overflow, and 9.723 seconds. It is 3.89% worse than H33 and
  rejects the preregistered 603.409M threshold.
- Exploratory group sizes 2 and 4 reach 657.104M and 645.686M HPWL. All three
  weighted variants underperform the original H33 membership, so group-size
  tuning and weighted greedy membership are exhausted.
- The H29--H35 outer reflection identifies a shared hidden assumption: exact
  transport always snaps node centers to destination-bin centers. This erases
  the non-legalized seed's sub-bin relative geometry even when a connected
  group shares a coarse destination.
- H36 will preserve each node's offset from the processed source-bin center by
  translating it by the source-to-destination bin-center displacement. Every
  candidate still requires strict exact-overflow descent.

## H36 offset-preserving exact transport

- The locked radius-16 run reaches 625.958M HPWL, 16.8458% exact overflow,
  and 10.140 seconds. Offset preservation improves H33 by 1.45% in HPWL and
  2.14 overflow percentage points, but rejects the locked 603.409M gate.
- Exploratory radius 32 reaches 620.681M/16.6019%/15.656s. Its gain over the
  matched H33 radius-32 result is only 1.64%, confirming that radius expansion
  does not unlock a larger offset-preservation benefit.
- Bin-center snapping is a real secondary loss, not the dominant one. The
  remaining inconsistency is that connected group members preserve their own
  offsets but still choose different displacement vectors.
- H37 should give a group one shared per-round rigid displacement proposal,
  retaining independent strict exact-overflow acceptance for every member.

## H37 rigid per-round group displacement

- Eight rigid rounds reach 597.238M HPWL, 25.8840% exact overflow, and 7.588
  seconds. HPWL improves 4.59% over H36, narrowly missing the 5% gate, but
  exact overflow fails the below-20% gate.
- Accepted moves increase slightly to 209,250, so the feasibility loss is not
  a move-count collapse. A shared vector cannot exploit heterogeneous residual
  capacity across all group members.
- Twelve rigid rounds reach 609.110M/24.4107%/9.404s. The weak overflow gain
  and HPWL regression rule out pure rigid continuation.
- Rigid displacement and independent offset-preserving transport are
  complementary: the former protects topology and the latter evacuates
  capacity. H38 should test a preregistered rigid-to-independent schedule.

## H38 staged rigid-to-independent transport

- The locked 8+4 schedule reaches 609.939M HPWL, 24.2978% exact overflow, and
  9.465 seconds. It retains the HPWL advantage but fails the below-20% gate.
- Earlier transitions 4+8 and 1+11 still leave 24.0942% and 23.4760% overflow,
  whereas a 12-round all-independent control reaches 16.2019%. Even one rigid
  prefix round creates residual capacity fragmentation that per-node cleanup
  does not undo.
- Scheduling is not the missing composition. The stronger hidden restriction
  is per-node strict acceptance: a rigid group's aggregate exact overflow can
  decrease even if one intermediate member's delta is positive.
- End rigid-round tuning. The next test should transact a whole group
  displacement atomically against exact aggregate overflow, with full rollback
  on rejection and exact final audits.

## Outer reflection after H36--H38

- A ten-candidate failure-analysis pass ranks atomic group transactions ahead
  of capacity reservations, group auctions, min-cost flow, subgroup retries,
  and group cycles. It is the simplest direct test of H38's transaction-
  granularity failure.
- H39 will apply one candidate displacement to a whole group against evolving
  temporary occupancy, then commit only if the summed exact overflow-area
  delta is strictly negative. Original coordinates and touched-bin occupancies
  provide exact rollback.
- This changes acceptance granularity, not the density objective: positive
  member-wise intermediate deltas are permitted only when the complete group
  has strict exact overflow descent.

## H39 atomic exact-overflow group transport

- The locked group-size-8 run accepts 92,305/94,771 atomic attempts and reaches
  509.126M HPWL, 34.2600% exact overflow, and 6.857 seconds. Atomic transactions
  improve H37 HPWL by 14.75% but fail feasibility.
- Group sizes 4 and 2 give 516.341M/33.3394% and 544.112M/30.5289%. Sixteen
  size-2 rounds reach 561.603M/28.2426%; granularity and continuation form a
  smooth tradeoff but do not reach 20%.
- Atomic acceptance remains about 95%, so rejected-group splitting is not the
  main bottleneck. Accepted groups can preserve only a tiny fraction of the
  first member's exact density gain while still satisfying strict aggregate
  descent.
- H40 should require a minimum exact aggregate-gain fraction and fall back to
  the first member's already-audited individual move when the group is density
  inefficient.

## H40 atomic gain retention with fallback

- The locked 0.5 ratio accepts 92,170/94,805 groups, falls back 2,680, and
  reaches 522.060M HPWL, 32.6484% exact overflow, and 6.943 seconds. HPWL
  passes but overflow remains far above 25%.
- The maximum ratio 1.0 still accepts 91,460/94,817 groups and reaches
  522.423M/32.3168%/6.883s. A global retention ratio has almost no additional
  feasibility leverage, so ratio tuning is closed.
- Atomic transactions currently move remote group members that do not overlap
  the active overfull source, consume unrelated capacity, and mark those nodes
  unavailable for their own source later in the round.
- H41 should form the atomic subgroup from the current source's exact
  contributors only; non-active members remain eligible elsewhere.

## H41 source-active atomic subgroups

- The locked run accepts 64,463/65,948 source-active atomic attempts and
  reaches 546.173M HPWL, 28.3515% exact overflow, and 8.006 seconds. It passes
  600M and improves H39 overflow by 5.91 points, but misses 20%.
- Sixteen rounds reach 553.817M/27.4672%/11.744s, showing a plateau. Combining
  ratio-1.0 fallback reaches 556.759M/27.0688%/8.193s; most local subgroups
  already pass the exact gain filter.
- Restricting atomic membership to current exact contributors is validated,
  but it is an early topology operator rather than a full feasibility solver.
- H42 should use a fixed source-active atomic prefix followed by fully
  independent H36 rounds; H38's negative global-rigid schedule does not move
  remote members in the same way and does not rule this out.

## H42 source-active atomic prefix cleanup

- The locked four-atomic plus eight-independent run reaches 572.409M HPWL,
  25.3415% overflow, and 10.064 seconds after 219,955 moves. Quality and time
  pass their proxy gates, but feasibility misses below 20% by 5.34 points.
- The suffix improves H41 overflow by only 3.01 points and remains 9.14 points
  worse than the matched all-independent control. A focused test proves the
  transition executes, so atomic-prefix scheduling is closed.
- Atomic acceptance is 63,171/63,493 (99.49%). Aggregate descent therefore
  accepts moves that reduce source excess while creating smaller excess in
  destination bins. H43 should require componentwise non-increasing exact
  overflow over every touched bin, with individual exact fallback.

## H43 componentwise atomic capacity

- The locked run reaches 625.161M HPWL, 16.8720% overflow, and 10.849 seconds.
  Feasibility passes below 20%, but quality misses the 600M gate.
- The exact guard rejects 58,482 trials; only 3,500/62,136 atomic attempts are
  accepted and 58,658 fall back. The endpoint nearly reproduces H36 independent
  transport, proving destination spill caused the atomic density plateau while
  a one-candidate hard guard discards the connectivity advantage.
- H44 should assign source-active subgroups over multiple candidate destinations
  using explicit exact residual capacities, not tune the guard or group size.

## H44 source-local exact-capacity bidding

- Eight candidates yield 622.606M HPWL, 16.8107% overflow, and 11.486 seconds.
  Density and time pass, but HPWL misses 610M and only 7,086 atomic groups are
  accepted versus the locked 10,000 prediction.
- Of 496,078 exact bids, 18,308 are feasible and 3,990 groups are rescued after
  their first candidate fails. The mechanism works, but improves H43 HPWL by
  only 0.41%; candidate breadth is not the dominant remaining limitation.
- Stop candidate-count tuning. H45 should batch all groups in a source against
  one occupancy snapshot and allocate scarce capacity by exact HPWL value.

## H45 source-batch exact capacity assignment

- The locked run reaches 624.363M HPWL, 16.8309% overflow, and 12.987 seconds,
  missing the 610M quality prediction and underperforming H44.
- Source snapshots expose 172,652 feasible bids from 497,424 trials, yet live
  capacity permits only 3,671/62,411 group commits. Most bids compete for the
  same residual capacity; exact-HPWL greedy ordering does not solve the conflict.
- Close strictly capacity-safe assignment. H46 should use exact density-neutral
  equal-shape anchor exchanges to change which cell identities occupy congested
  source positions before independent evacuation.

## H46 group-guided density-neutral identity exchange

- The invariant is validated: 6,683 accepted complete exchanges permute 29,560
  identities, reduce exact HPWL from 41.510168M to 41.497598M, and leave exact
  overlap bitwise unchanged at 96.8876491423%.
- The pass evaluates 226,128 complete exchanges and costs 49.266 seconds for
  only 0.01257M HPWL gain. Eight subsequent H36 rounds end at
  625.990M/16.7886% in 59.744 seconds, failing the locked 600M and 20-second
  gates.
- Local connected-group identity exchange is closed as a clean-start primary
  direction. The next representation must allocate two-dimensional global
  capacity while explicitly minimizing displacement/topology damage, rather
  than tune more local group candidates.

## H47 coarse two-dimensional capacity flow

- The global L1 plan balances its center-area model efficiently: 1,871 edges
  and 201,105 moves reduce 64x64 coarse overflow from 94.8007% to 0.5257% in
  3.218 seconds.
- Exact rectangle overlap remains 76.2850% and HPWL reaches 530.789M; four fine
  cleanup rounds end at 522.926M/40.8354% in 8.537 seconds. All locked quality
  and feasibility predictions fail.
- Offset preservation clones the seed's congested fine-bin phase into each
  coarse sink. Center-area balance cannot represent exact rectangle occupancy;
  H48 must allocate live exact fine-bin anchors inside each coarse flow sink.

## H48 exact fine-anchor flow realization

- Live rectangle-anchor search resolves H47's phase aliasing: 200,915 moves
  reduce exact overlap from 96.8876% to 7.2900% in 6.804 seconds while the same
  coarse plan is retained.
- Four exact cleanup rounds reach 530.068M HPWL and 5.21172% overflow in 9.595
  seconds with one GP iteration. All H48 intermediate gates pass, but HPWL is
  still 7.17x the final a1 limit of 73.9522M.
- Exact fine capacity is solved; anonymous coarse commodity assignment now
  dominates topology damage. H49 should jointly assign each source's fixed
  sink quotas by exact HPWL opportunity cost rather than process edges greedily.

## H67-H70 exact breakpoint density-only continuation

- H67 vector-bin prices lower exact overflow monotonically from 96.8876% to
  83.1320%, with negative priced-area deltas on every sweep, but HPWL rises
  from 41.510M to 204.567M in 187.15s. Per-bin pressure does not preserve
  two-dimensional net connectivity and is rejected as the primary solver.
- H68 exact 64-to-512 continuation reaches 69.2734% final overflow and
  341.716M HPWL in 277.11s. Coarse scale crosses fine-bin plateaus, but the
  coarse stage damages seed topology too strongly for local recovery.
- H69 pair blocks at fine scale reach 67.3360% overflow but worsen HPWL to
  377.459M; fine-stage block-size tuning is closed.
- H70 strict exact overflow-first ordering preserves HPWL at 61.430M but
  reaches only 86.1079% overflow in 256.49s. Sparse acceptance rules out
  squared excess energy as the dominant plateau cause.
- The remaining limitation is a global capacity-assignment/topology conflict,
  not smoothing or a missing scalar price. Future work should use an exact
  topology-aware fine-anchor operator rather than more local tuning.

## H72-H73 electrostatic homotopy

- Added a reusable homotopy driver with exact HPWL subgradient, exact
  rectangle-bin overlap energy/overflow, DREAMPlace DCT Poisson electrostatic
  energy, and a small quadratic center anchor for strong convexity. `mu` is
  staged from 16 to zero; Adam and SGD are selectable.
- With fixed `lambda_overlap=1`, a1 improves only to 92.90% exact overflow.
  Adding epsilon-active overlap directions (4 bin pitches) and
  `lambda_overlap=100` reaches 78.17% at 247.3M HPWL after 240 iterations.
- Increasing the grid from 128 to 512 at the same 240 iterations worsens
  exact overflow to 85.98%; fine-grid active-set sparsity dominates the
  nominal resolution gain.
- The effective improvement is a two-parameter homotopy: start with
  electrostatic spreading and ramp the exact-overlap weight from zero to its
  terminal value as `mu` decreases. At 600 iterations, terminal lambda 100
  reaches 63.59% exact overflow / 306.9M HPWL in 49.3s; lambda 150 reaches
  60.58% / 384.5M. The lambda-100 point is the best practical tradeoff.
- Homotopy substantially improves the local density/HPWL frontier, but it is
  still far above the 10% overflow target. A topology-aware exact capacity
  stage remains necessary after the electrostatic continuation.

## H76-H77 adaptive homotopy

- H76 updates weights at every raw step after its trigger. Lambda reaches 150
  in 75 iterations while mu collapses from 16 to numerical zero; the best a1
  point is 390.105M HPWL / 61.0140% overflow in 99.43s. This is too abrupt to
  be a well-resolved continuation.
- H77 changes weights only on ten-iteration windows and reaches 226.527M /
  66.1932% in 102.09s. The lower HPWL confirms that slow mu decay preserves
  topology, but lambda reaches only 75 because improving windows freeze it.
- The next controller should make lambda proportional to the remaining error
  above the fixed 7% target while retaining windowed mu decay.
- H78 confirms that controller: 328.825M HPWL / 56.3286% overflow in 95.08s,
  improving H76 by 61.3M HPWL and 4.69 overflow percentage points. Lambda
  reaches 150 smoothly while mu falls to 4.46e-4 at window boundaries.
- H79 then applies eight exact 512-grid block-breakpoint sweeps and reaches
  357.287M / 18.0117% in 352.45s. This improves H74 by 11.92 overflow points,
  but late gains diminish to 0.30 points per sweep and total runtime is high.
- H80 adds the previously supported half-bin singleton residual after every
  block sweep. It reaches 359.198M / 17.0796%, but takes about 9905s, so the
  mechanism is useful for diagnosis and not acceptable as the main runtime
  path.
- A semantic audit found the H72 driver was also using epsilon-active overlap
  directions. The corrected H81 driver sets overlap epsilon to exactly zero,
  confines epsilon-active selection to the HPWL subgradient, and forces mu to
  exactly zero in the final stage.
- H81's strict run reaches 62.115M / 82.9201% in 82.24s. The final 100 steps
  are exactly `mu=0, lambda=150`; exact overlap active-set starvation, rather
  than a hidden smoothing term, is now the limiting mechanism.
- H83 tests a 64x64 exact initial grid and reaches 61.413M / 81.2527% in
  79.39s, so coarse resolution alone does not cure the starvation.
- H82's exact 512-grid breakpoint stage reaches 201.389M / 62.1275% after
  eight sweeps, with 2.2--4.5 overflow points of gain per sweep. Continue the
  same exact oracle before changing the non-smooth model.
- H84 continuation reaches 299.883M / 44.6880% after twelve additional
  sweeps; gains remain above one point per sweep but HPWL is high.
- H85 exact-feasible HPWL recovery lowers that endpoint to 171.132M while
  preserving 44.6880% overflow in 6.25s. H86 tests breakpoint refinement from
  this better topology instead of adding any pseudo-gradient.
- H86 reaches 226.605M / 36.8518% after eight exact sweeps from H85. The
  alternating recovery/refinement frontier is improving, so H87-H88 continue
  with an exact 36.85% cap.
- H88 reaches 210.419M / 31.6272% after recovery then eight exact sweeps,
  improving H86 by 5.22 overflow points. Continue one more alternating cycle.
- H90 reaches 203.753M / 28.4927% after the second cycle, another 3.13-point
  exact overflow gain. The alternating frontier remains active.
- H92 reaches 199.747M / 26.4031% after the third cycle, but the gain is only
  2.09 points. H94 tests wider exact breakpoint visibility after recovery.
- H94's 32-bin exact radius reaches 277.102M / 15.6757% in 795.67s, a much
  stronger density move but with HPWL damage. H95-H96 recover under the exact
  cap before repeating the wider radius.
- H96's second 32-bin cycle reaches 278.426M / 11.2371% in 824.81s. The strict
  exact path is now within 4.24 points of the 7% gate; H97-H98 make one final
  recovery/refinement cycle.
- H98 reaches 274.816M / 9.47107% after the final 32-bin cycle. It passes the
  10% exploratory target but misses the locked 7% gate by 2.471 points.
  H100 tests a 64-bin exact breakpoint radius after cap recovery.
- H100 reaches 316.025M / 6.24086% at the third 64-bin exact sweep and stops
  below 7%. This is the first strict nonsmooth, fixed-macro-aware feasible
  checkpoint; HPWL and runtime remain far above target.
- H101 extended coordinate recovery reaches 275.806M / 6.24097% in 16.44s;
  it preserves feasibility but misses the below-200M prediction.
- H102 equal-shape swaps reach 169.667M with exact overflow bitwise invariant
  at 6.24097%. H103 alternates coordinate and swap recovery.
- H103 reaches 156.426M / 6.24094%; the below-120M prediction fails and
  coordinate recovery is nearly stationary. H104 widens the invariant swap
  radius.
- H104 radius-32 swaps reach 154.113M / 6.24094% in 404.74s, only 2.31M
  improvement. H105 tests exact incident-net assignment recovery.
- H105 assignment recovery reaches 154.050M / 6.24094%; only 0.063M HPWL
  improvement. The strict feasible topology is locally saturated, and the
  remaining a1 HPWL gap is inherited from density construction.
- The existing non-legalized position-seeded recursive bisection checkpoint
  is 224.399M HPWL / 1.40758% exact overflow. H106 evaluates recovery from
  this better-constrained topology as a separate construction path.
- H106 reaches 170.935M HPWL / 1.4080% exact overflow in 189.62s. It is more
  feasible than H105 but has 16.88M worse HPWL, so H105 remains the strict a1
  frontier point at 154.050M / 6.24094%.
- H107's pre-run gradient audit finds that the old normalized coefficients are
  not comparable: `lambda=150` makes the sparse exact-overlap term about
  7,450x the HPWL term and 340x the `mu=16` electrostatic term on the initial
  state, while only 33.7% of nodes are overlap-active. Scalar magnitude and
  active-set coverage are separate failure modes.
- The next homotopy controller must initialize physical coefficients from
  gradient-norm ratios, then adapt dimensionless lambda/mu multipliers. The
  final stage still forces the physical electrostatic coefficient to zero.
- H108 validates the unit conversion but reaches only 53.526M / 88.1673%.
  Initially, control one gives equal `2.7961e-6` norms for all three terms;
  at the final `lambda_control=150`, overlap is only 1.33x HPWL rather than
  H107's 7,450x imbalance.
- The remaining controller error is directional: at iteration 120,
  `mu_control=16` has reached a high-overflow equilibrium where electrostatic
  and HPWL norms are both about `5.7e-6`. Decaying mu at this point raises
  overflow. High-overflow stagnation must first boost mu; decay belongs after
  a spreading milestone or bounded maximum.
- H109 confirms the boost direction but not its magnitude: control values
  16, 32, 64, and 128 move the best exact overflow from about 88.8% to
  85.5066% at 135.123M HPWL. Decay reverses the density gain, and the final
  `mu=0` endpoint is 55.381M / 88.6982%.
- A global best-overflow restore can select a positive-mu checkpoint after the
  final stage. Such a point is only a continuation seed. Strict reporting must
  restore it before the final stage and then select output exclusively among
  `mu=0` iterates.
- H110 enforces that reporting rule and reaches 160.483M / 85.0842% in
  127.39s. Increasing mu control through 256, 512, and 1024 saturates near
  85.1%, ruling out insufficient scalar range.
- All 543 nodes above area 10,000 in adaptec1 are fixed terminals; movable
  macro area dominance is not the explanation. The stronger conflict is the
  per-node center quadratic, which scales with mu and directly penalizes the
  spreading electrostatics is meant to create. A centroid-only translation
  anchor preserves relative spreading and is the correct nullspace treatment.
- H111 validates that correction: the strict 128-grid final point is
  74.031M / 54.7432% in 113.43s, versus H109's roughly 89% final overflow.
  At iteration 200 it already reaches 66.8% with 71.3M HPWL.
- Exact 512-grid audit of H111 is 66.1264%, an 11.38-point gap from the coarse
  metric. The remaining mismatch is spatial resolution: 128-bin charge
  enlargement is about 118 units while standard cells are roughly 13 x 12.
  Repeating the corrected homotopy at 512 bins is now justified despite old
  pre-correction 512-grid failures.
- H112 rejects 512 bins for the full trajectory: it reaches 67.939M / 68.4957%
  in 84.05s and is still worse than H111's 66.1264% fine-grid audit. Fine
  fields improve fidelity but lose long-range transport.
- Snapshot interval 50 cuts substantial placement-text I/O and keeps the
  512-grid 600-step run below 85s. The next supported design is multiresolution:
  coarse global spreading followed by a short fine electrostatic continuation,
  then an exact mu-zero stage.
- H113 validates multiresolution: 150 fine positive-mu steps move H111's
  66.1264% audit to 25.9499%, and the strict mu-zero stage selects
  94.584M / 16.6939% in 31.03s.
- Exact Adam reaches its density minimum early and then trades overflow for
  HPWL. Restarting from the exact checkpoint resets moments and exposes a new
  piecewise overlap active set while preserving the exact nonsmooth model.
  This is far cheaper than the block-breakpoint sweeps used in H82-H100.
- H114 validates the restart: a 14.81s exact-only run reaches
  123.241M / 4.18619%. Its iteration 99 is 114.718M / 6.0521%, already showing
  a better feasible HPWL point than the minimum-overflow selector returns.
- The remaining task on a1 is now feasible HPWL recovery, not density. An
  exact 6.9% cap provides 2.71 overflow points of slack for coordinate moves
  while preserving the final nonsmooth model.
- H115 converts that slack into a 94.119M / 6.899985% point in 16.60s, a
  23.6% HPWL reduction from H114. Coordinate capacity is now saturated.
- Equal-shape swaps are the next exact invariant: they preserve every bin's
  occupancy while changing cell identities. H102's first three sweeps were
  strongly front-loaded, suggesting a short quality test can approach the a1
  HPWL target without another density solve.
- H116 rejects that prediction on the new topology: three sweeps cost 43.34s
  and improve only 94.119M to 92.807M at bitwise-invariant overflow. The swap
  representation is not the primary remaining HPWL mechanism.
- Identity changes can nevertheless reopen exact-cap coordinate moves. One
  final alternating recovery is justified; repeated swaps are not.
- H117 gains only 0.448M after H116 and ends the swap/recovery branch at
  92.359M / 6.9%. The final feasible topology is locally saturated.
- Most HPWL damage enters during H114's direct 16.7% to 4.2% density push.
  A wider intermediate recovery band before the next exact restart is the
  next topology-preserving split; it changes scheduling, not the model.
- H118 reaches 88.023M / 25.0% in 13.79s, missing its 78M prediction but
  improving the comparable pre-restart topology by about 4.9M. Only a fresh
  exact density restart can determine whether that gain survives feasibility.
- H119 selects 112.424M / 5.41092% in 14.25 seconds. Its final iterate is
  101.822M / 6.88191%, but both are worse than H115's 94.119M / 6.899985%.
- Protecting HPWL at a 25% intermediate cap does not improve the final feasible
  topology. The H118--H119 branch is closed; H115 remains the strict a1 best.
- H120 will test whether H115's coordinate recovery saturates because its
  largest candidate is only eight bins and its four-level line search couples
  x and y. All candidates remain exact-HPWL improving and exact-overflow
  audited.
- H120's 16-, 32-, and 64-bin variants reach 94.054M, 93.994M, and 94.025M,
  all below 6.9% overflow in 16.40--18.35 seconds. Step 32 is a new strict best
  but improves H115 by only 0.125% and misses the below-90M prediction.
- Fixed-direction step range is rejected as the main limitation. H121 will
  retain the joint negative-subgradient candidate and add x-only and y-only
  coordinate projections, each subject to the same exact objective and exact
  overlap-cap audit.
- H121 confirms axis coupling as a real obstruction: coordinate projections
  improve H120-step32 by 2.597M to 91.3975M / 6.89999992% in 38.50 seconds.
  It misses the below-90M prediction but establishes a new strict a1 best.
- The last five sweeps gain only 4.32k, so more identical recovery is rejected.
  H114 iteration 75 is a stronger seed: 117.535M / 4.5754%, versus the selected
  minimum-overflow seed's 123.241M / 4.1862%. H122 tests this checkpoint choice.
- H122 audits the saved iteration-75 placement at 117.332M / 4.6403% and ends
  at 91.3746M / 6.899995% in 39.45 seconds, only 22.9k better than H121.
- Moderate checkpoint timing is therefore not important. H123 will test the
  unsaved H114 iteration-99 point at 114.718M / 6.0521%; if it returns to the
  same 91M basin, checkpoint selection is closed.
- H123 reproduces the final checkpoint and recovers to 91.4087M / 6.9000%,
  worse than H122. Three distinct final-stage checkpoints converge to the same
  basin, so checkpoint timing is closed.
- Code audit finds that final-stage lambda is forced to its maximum and barred
  from adaptive updates by `stage != final`. The output selector also minimizes
  overflow rather than HPWL within the feasible band. H124 will correct both
  while keeping mu exactly zero and overlap exact.
- H124's final mu-zero controller returns 105.937M / 4.9963% in 20.40 seconds,
  improving H114's old selected output by 17.304M. Lambda control falls from
  150 to 90 as exact overflow stays below the 6.5% lower threshold.
- Quality, feasibility, and runtime predictions pass, but the upper lambda
  branch is not reached within 150 steps. H125 will test whether the topology
  gain survives axis-separated exact-cap recovery to 6.9%.
- H125 ends at 91.1194M / 6.899998% in 24.59 seconds. It is a new strict best,
  but improves H122 by only 0.255M and misses the below-88M prediction.
- Distinct exact-stage seeds still collapse to the same cap-recovery basin.
  H126 will extend the adaptive exact trajectory to 300 steps so lambda can
  reach the upper band and exercise both controller directions.
- H126 reaches 95.4261M / 6.3346% in 29.81 seconds, with lambda control reduced
  from 150 to 10. Quality, feasibility, and runtime predictions pass.
- Overflow rises from 5.40% at iteration 280 to 6.33% at 299 but does not yet
  cross 7%, so H127 extends the same run to 400 steps to audit the controller's
  upper response and feasible-band stability.
- H127 selects 93.6776M / 6.99568% in 40.11 seconds. Lambda control follows
  150 to 10 to 15 to 20, and overflow remains bounded near 7% after crossing.
  The bidirectional exact-stage controller is validated.
- H113's preceding fine continuation is not adaptive in practice: mu control
  stays 16 and lambda stays zero for 150 steps, raising HPWL from 74.03M to
  92.88M before an abrupt handoff. H128 tests an earlier iter25 handoff to the
  validated physical-mu-zero exact controller.
- H128's direct iter25 handoff stalls near 12.7--13.5% and falls back to
  105.516M / 11.7755% after 450 exact steps. Exact overlap alone is too sparse
  from a 39% state.
- The two extremes bracket the needed mechanism: retain electrostatic transport
  past 39%, but trigger continuous mu decay and exact-lambda growth by overflow
  rather than waiting for progress to stall. H129 implements that handoff.
- H129 executes the handoff but is rate-imbalanced: mu control drops 16 to
  0.184 while lambda rises only to 93.19. HPWL reaches 84.118M, but overflow
  rebounds from 28.38% to 35.29%, and strict output is 92.048M / 21.9608%.
- H130 keeps mu longer with decay 0.9 and makes lambda take over three times
  faster with step 15. This targets the observed rebound without changing any
  objective term or final-stage semantics.
- H130 improves the strict minimum to 90.590M / 19.8564%, but misses 10%.
  Overflow bottoms near iter130 and then worsens while mu keeps decaying from
  about 6.2 to 1.95 even though lambda is already near its maximum.
- H131 will make mu decay response-aware: a worsening exact-overflow window
  restores one decay step, while an improving window continues decay. The
  reported final stage still forces mu exactly to zero.
- H131's rebound guard oscillates mu around 2.7--5.0 without improving the
  19.86% minimum; strict output repeats H130. Scalar one-step rebound control
  is rejected.

## H195 bounded-runtime audit of coarse exact transport

- H195 was pre-registered from the raw a1 `.pl` and launched with 40 threads,
  one 64x64 coarse-flow pass, and exact fine-bin anchor screening.
- The implementation was changed so exact-anchor mode takes one coarse-bin
  step instead of a full source-to-sink jump. This preserves the intended
  local transport semantics and all core tests still pass.
- The run remained CPU-bound for more than ten minutes before producing any
  metric file, so it was stopped on the runtime gate. The bottleneck is the
  product of many source nodes, flow edges, and fine-bin exact trials.
- Global optimization remains warranted as a direction, but it must be a
  bounded candidate generator: cap source regions and commodities, evaluate a
  fixed shortlist, and audit one rigid/net-connected move at a time. Full
  individual-node exact flow is rejected as a main-stage mechanism.
- H132 tests strong coexistence instead: lambda still reaches 150 quickly, but
  mu decays by only 0.99 per window and remains a long-range transport bridge
  until it is forced exactly to zero in the final stage.
- H132 reaches 92.462M / 11.1664% before/following the mu-zero boundary in
  40.77 seconds. It passes the pre-final 12% gate but remains infeasible.
- H133 doubles mu control to 32 while retaining slow decay and fast exact-
  lambda takeover, testing whether stronger transport enters 7% before the
  final term is removed.
- H133 reaches 96.495M / 7.6153% at the mu-zero boundary in 42.89 seconds,
  improving H132 by 3.55 overflow points but missing feasibility by 0.615.
- H134 will test mu control 40 as a narrow strength extrapolation; final output
  remains restricted to physical mu zero.

## H150-H160 exact fine recovery and breakpoint oracle

- H150 activates stability handoff on the H140 `iter_1050` continuation. Exact
  lambda reaches 150 and overflow improves from 81.9% to 47.8%, but no
  physical-mu-zero feasible checkpoint appears.
- H151/H152 increase the bridge to mu control 32. Without forced final lambda
  the result is 46.1%; forcing lambda=150 still stops at 41.9%. H153 starts
  exact lambda at 80% overflow while retaining the strong bridge, but ends at
  40.9%. These runs close the `iter_1050` fine checkpoint: its 82% fine-grid
  concentration cannot be repaired by stronger electrostatic transport.
- H154 tests 128-bin exact recovery from the stronger H121 physical checkpoint
  and improves only 0.014M, showing fixed radius is not the main obstruction.
- H155 adds an exact breakpoint oracle that enumerates rectangle-edge grid
  crossings along the HPWL subgradient. It lowers 91.397M to 90.574M at 6.9%
  overflow without changing the objective or adding smoothing.
- H156 removes the degree-100 exclusion and uses raw exact HPWL subgradients;
  it reaches 90.493M. H157 equal-shape exact identity exchange lowers this to
  90.138M without changing overlap. H158's wider exchange reaches 89.943M.
- H159 spends the permitted 6.99% exact-overflow slack and reaches 89.683M;
  H160 adds another density-neutral exchange stage and reaches the current
  strict a1 best, 89.642M / 6.99% in 22.86s. Exchange and breakpoint gains
  saturate, so this branch is closed as a route to the 73.22M baseline.
- All H150-H160 reported candidates use exact max-minus-min HPWL, exact
  rectangle-bin overlap with fixed macro capacity, epsilon=0 for overlap,
  physical mu=0 for final candidates, and no legalization or legality check.

## H161-H167 fresh a2/a3 transfer

- H161 is a fresh a2 raw-`.pl` coarse chain. The a1 weak-mu schedule reaches
  only 79.57% overflow (78.13M HPWL) after 57.2s, showing benchmark-specific
  transport sensitivity.
- H162 continues from that own checkpoint with mu=40 fine transport and reaches
  8.3189% overflow at the physical mu-zero boundary (144.899M HPWL).
- H163 initially exposes a recovery guard bug: an infeasible seed above the
  requested cap had no accepted moves. The exact guard is now corrected to
  allow only overlap-nonincreasing HPWL moves until the cap is crossed.
- H164 then enters the 6.99% band at 135.277M in 29.6s; H165 equal-shape
  exchange reaches 134.965M / 6.99% in 32.9s. This is the current a2 best,
  but HPWL remains 64.1% above the 82.22M paper baseline.
- H166 fresh a3 weak-mu coarse ends at 97.17% overflow; H167 fresh direct
  mu=40 fine transport ends at 95.96% and 579.8M HPWL. The former starves
  overlap transport, while the latter pushes cells to die boundaries. A3
  requires a bounded transport/anchor schedule before further recovery.

## H168-H171 bounded transport and a4 transfer

- H168 adds a strong per-node anchor to fresh a3 mu=16 transport. It prevents
  the HPWL explosion, but exact overflow remains 99.92%, so the a3 issue is
  not only boundary collapse; exact density transport is effectively inactive.
- H169 fresh a4 coarse reaches 59.97% overflow and 204.36M HPWL in 110.1s.
  H170 fine mu=40 crosses the 7% requirement (6.814%) but damages topology to
  476.28M HPWL. H171 exact breakpoint recovery spends the 6.99% cap and lowers
  this to 436.29M / 6.99% in 47.6s, still far above the 174.08M baseline.
- A3/a4 therefore need a calibrated density transport implementation or
  multi-resolution schedule; continuing scalar mu sweeps would not be
  evidence-driven.

## H196-H197 coarse-flow accounting correction

- H196 bounded the candidate generator and finished in 18.87s, but accepted
  no exact moves because the coarse load was inconsistent with the exact
  rectangle oracle. The raw `.pl` audited at 0% exact overflow while center-bin
  area assignment reported 99.925% coarse overflow.
- H197 changed coarse load to exact rectangle-to-coarse-bin intersection area.
  The audit now agrees at 0% and produces no flow edges from the exact-feasible
  raw state. This is a correctness fix, not a smoothing change.
- The first GP update still destabilizes the raw feasible state (96.669M and
  99.995% overflow), so raw initialization should be treated as an HPWL
  recovery problem with an exact-overflow guard, not as a density-flow source.

## H198 raw HPWL recovery and initialization-domain audit

- The exact recovery stage reduced a1 HPWL from 104.924M to 102.776M in
  27.7s with reported incremental overflow zero, but the raw ISPD `.pl` puts
  all movable cells at `(0,0)`, outside the core rows.
- Because exact density ignores rectangles fully outside the layout, the raw
  state falsely audits at 0% overflow. Final boundary clamping exposes the
  actual 99.832% overflow. H198 is rejected as a valid raw-coordinate solve.
- The pre-GP feasible-checkpoint selector is retained for in-domain staged
  states. Fresh experiments must use the deterministic in-domain initializer
  from raw benchmark metadata, never a prior result or ePlace seed.

## H199 direct fine ablation and compactness diagnosis

- Fresh direct H186-style 512-bin homotopy reaches only 154.439M / 70.0743%,
  while H182->H186 reaches 112.848M / 8.3045%. The coarse prefix is a necessary
  low-frequency transport stage, though its 1200-step duration can still be
  shortened in a separate ablation.
- Uniform contraction of all movable nodes is not an HPWL descent direction
  because nets tied to fixed macros and terminals oppose global scaling. A
  useful compactness direction must preserve those anchors and remain subject
  to exact HPWL and exact overlap acceptance.
- Even excluding direct fixed-terminal neighbors yields zero accepted global
  contraction moves. Compactness must be localized per node/net; a single
  affine contraction cannot represent the placement's heterogeneous anchors.
- Local incident-net and support-centroid directions are effective when used
  as exact-audited additions to node-wise recovery: H202 improves H187 from
  103.715M to 103.218M at 7% overflow. The gain is real but far smaller than
  the 30M gap to the paper baseline, so compact directions are a polishing
  mechanism rather than the missing global topology mechanism.
- H182 pre-final checkpoints at 300--500 steps cannot feed H186 below 31.7%
  overflow. The missing event is the coarse physical-mu-zero/lambda handoff,
  so prefix shortening must move that handoff earlier rather than simply cut
  the trajectory at an electrostatic checkpoint.
- A fresh early coarse exact handoff at step 550 improves the subsequent fine
  result to 10.75%, confirming the handoff mechanism, but is dominated by the
  full prefix in quality and measured runtime. Prefix shortening alone does
  not solve the HPWL gap.
- Exact non-rigid contraction of small nets is feasible and fast, but after
  H202 it gains only 0.015M. The remaining 30M paper-baseline gap is primarily
  a global topology/path issue introduced during density transport, not merely
  excess geometric support that a feasible-band postprocessor can collapse.
- H186's final-stage selector itself accounts for a meaningful part of the
  loss: recovery from its earlier iteration-325 topology reaches 101.285M,
  versus 103.218M from the selected checkpoint. Most of this gain appears in
  eight sweeps, making checkpoint choice plus early stopping a stronger lever
  than additional post-cap contraction operators.
- Eight sweeps from iteration 325 retain nearly all useful recovery and meet
  the full staged runtime/iteration gates at 101.532M / 7%. Late fine
  checkpoint selection is now a validated mechanism, not merely a diagnostic.
- The later iteration-349 topology is stronger: eight sweeps reach 100.909M
  but stop at 7.117% overflow. This supports late checkpoint selection and
  indicates that only a small exact-recovery extension is needed for feasibility.
- Extending iteration-349 recovery to ten sweeps reaches 100.843M / 7.102%,
  but the last two sweeps remove only 0.0155 overflow percentage points. The
  late local operator is entering a strongly diminishing tail rather than
  approaching the cap linearly.
- Twelve sweeps stop at 100.808M / 7.0968%; the last two remove only another
  0.00516 percentage points. More identical late recovery is rejected. The
  remaining exact-overlap violation must be addressed before the local active
  set is exhausted, ideally while the fine trajectory still has 10--15%
  overflow and more spatial slack.
- A weak per-node common-center anchor is dominated by H186: iteration 349 is
  113.760M / 13.094%. Compactness must preserve the incoming net/fixed-terminal
  topology; a universal geometric center is the wrong reference. A proximal
  anchor to each node's H182 position is the next minimal vanishing homotopy.
- Initial-position proximal continuation confirms that preserving incoming
  topology matters: at iteration 349 it saves 3.32M HPWL versus H186 while
  using 1.57 additional overflow percentage points. The gain must now survive
  exact mu-zero recovery; the auxiliary anchor itself is absent from the final
  model.
- Exact recovery from the proximal path reaches 97.904M HPWL, a further 2.90M
  gain over H210, but stalls at 7.611% and costs 52.96 seconds. Proximal
  topology is supported, while its density debt must be reduced during the
  simultaneous exact-lambda takeover rather than by longer local recovery.
- Raising proximal lambda's ceiling alone is inert through the useful path:
  control reaches only 154.66 by iteration 325 and the checkpoint remains
  104.495M / 10.732%. The next control variable is takeover rate, not maximum.
- Doubling the takeover step shifts iteration 325 to 106.392M / 9.625%, a
  clear density/HPWL trade from H212 and a slightly denser-feasible frontier
  than H208 at similar HPWL. Earlier exact lambda is effective; whether its
  topology survives strict recovery is the next discriminator.
- The early-lambda proximal checkpoint survives exact recovery: H216 reaches
  98.435M / 7% and improves the prior branch by 2.37M. Its only failure is a
  3.95-second staged runtime overrun; six of twelve sweeps already provide
  98.862M / 7%, making exact early stopping the next direct correction.
- Six-sweep early stopping passes the complete staged gate at 98.862M / 7%
  in about 122.48 seconds and 1,533 iterations. Proximal continuation is now a
  validated H186-branch mechanism, but its 35.02% paper-HPWL gap means further
  late local polishing is not the high-leverage next step.
- The numerically lower H160 result is not valid under the current initialization
  contract: its lineage reaches H111, whose protocol explicitly uses a
  non-legalized ePlace seed. It is retained as historical evidence only. The
  valid fresh-from-raw best is H217 until improved by a clean staged run.
- A fresh 43.049M HPWL seed does not work with H182's unmodified control:
  gradient balancing reduces electrostatic scale by 5.36x, overflow remains
  85.63%, and lambda never hands off. Any initialization change must recalibrate
  physical mu from measured gradient scale, not reuse its dimensionless value.
- Multiplying seeded coarse mu by the measured 5.36 ratio restores transport:
  iteration 1100 reaches 60.618M / 71.708%, about 5M lower HPWL than H182 at
  comparable coarse density. This validates scale-normalized HPWL seeding;
  early checkpointing is required to retain the complete runtime budget.
- The same scale issue recurs at the seeded fine handoff: uncalibrated H220
  reaches only 36.68% overflow at iteration 300 despite 83.04M HPWL. Both mu
  and exact lambda must be normalized to physical gradient coefficients after
  every initialization/resolution change; raw controller values do not transfer.
- Calibrating both fine coefficients recovers a useful path: H221's mu-zero
  entry is 109.426M / 7.752% with a 90.13-second complete fresh prefix. The
  iteration-300 prediction misses at 10.965%, but low residual density debt
  makes exact recovery the decisive follow-up.

## H230-H246 overflow hysteresis refinement

- A 15% exact-overflow excursion with shared net-density rigid moves repairs
  topology, while a 20% excursion creates too much retightening debt.
- The moderate bridge has an active-set role, not just a lambda-ramp role.
  Switching at H233 iteration 60 eventually starves at about 7.44% even with
  very large effective lambda, while iteration 80 crosses the exact gate.
- At the iteration-80 handoff, scale 4 and 4.5 preserve better HPWL but stop
  just above feasibility. Scale 4.75 is the measured boundary: H245 reaches
  88.9171M / 6.99968% within the same 200 bridge-plus-strong steps as H234.
- One exact feasible-band sweep is strongly complementary to simultaneous
  Adam: H246 improves H245 by 1.3371M to 87.5800M while exact overflow remains
  at 6.999998%. This is the new strict a1 best under the estimated 133.94s /
  1,699-iteration staged accounting.
- Runtime margin is now effectively exhausted. Further quality stages require
  deleting work from the prefix or demonstrating a faster 40-thread replay;
  simply appending another recovery sweep would violate the current time gate.
- A fresh 40-thread replay shows that the current host is oversubscribed at 40
  OpenMP threads: HPWL seeding doubles from 5.80s to 11.66s and coarse
  homotopy is also substantially slower. The allowed maximum is not the
  machine-optimal setting; 16 threads remains the valid local runtime path.
- Moving the coarse-to-fine handoff 100 iterations earlier does not preserve
  the strict chain efficiently: H248 selects 110.160M / 7.9948% and its fine
  replay takes 47.12s. The H221 iteration-1100 handoff remains necessary under
  the current controls.

## H261-H272 in-process late-stage continuation (2026-08-22)

- A same-process moderate-to-strong switch avoids the H253→H254 coordinate-only
  restart damage, but preserving Adam moments alone stalls at 10.31% overflow
  (H261). Increasing late lambda fourfold without changing optimizer state
  still stalls at 10.23% (H262), identifying inherited Adam history as a
  principal active-set obstruction.
- Resetting only Adam at the switch, while preserving coordinates and the
  lambda controller, crosses the exact 7% band at 89.262M (H265). This is a
  valid non-smooth continuation mechanism: it changes search state, not the
  objective or density oracle.
- Earlier takeover is better. Switching at iteration 70 reaches 88.497M /
  6.997% (H270), and four exact node-only sweeps reach 87.248M (H271). A
  fifth sweep gives only 2.8k additional gain (H272), so the late active set
  is entering a diminishing tail.
- Long node moves (H264), wider line search (H268), and one late net-block
  sweep (H269) do not beat H271. Net-aware blocks remain useful earlier in the
  route, not after this feasible band.
- H272 is about 17s cheaper than the H257 route under the same prefix estimate,
  but its 87.220M HPWL remains 18.0% above the paper's 73.22M baseline. H257
  remains strict quality best (87.119M); neither meets the paper HPWL gate.

## H273-H279 conservative-step and near-cap recovery (2026-08-22)

- Reducing the iteration-70 late step multiplier from 0.5 to 0.4 lowers the
  HPWL rebound but ends at 7.2955% overflow (H273). One and four exact
  recoveries reach 7.1332% and 7.0618%; H275 attains 87.118M HPWL but remains
  infeasible.
- Four more hard-cap sweeps (H276) remove only 0.0135 overflow points. A 7.5%
  exact hysteresis excursion (H277) lowers HPWL to 86.976M, but cap-7%
  retightening (H278) stalls at 7.4538%. Small hysteresis therefore creates
  topology debt faster than exact recovery can repay it.
- A near-cap mixed net-block/node sweep (H279) improves H275 to 87.041M /
  7.0396%, but still misses feasibility. Net blocks are marginally useful near
  the cap, not a complete active-set rescue.

## H280-H281 earlier takeover and relaxed-prefix ablations (2026-08-22)

- Moving the Adam reset to iteration 60 (H280) fails at 7.5316%, confirming
  that the useful switch boundary is later than 60 for the two-sweep prefix.
- Restoring H230's third relaxed sweep (H281) lowers the bridge HPWL but ends
  at 7.6340% under the same 201-step iteration-70 takeover. The extra relaxed
  topology carries too much exact density debt; H252 remains the preferred
  runtime/feasibility prefix.

## H282-H284 global candidate experiments (2026-08-22)

- Existing exact-audited auction/group transport (H282) reduced overflow only
  0.44 points but raised prefix HPWL by 5.7M; the selected result was 89.527M.
  Connectivity ordering alone does not make local destination moves topology
  safe.
- Bounded 64x64 coarse flow with exact anchors (H283) accepted 1,987 moves,
  but still raised prefix HPWL by 0.84M and selected 88.859M. Exact anchors
  preserve capacity accounting, not enough net topology.
- Recursive hypergraph bisection/FM (H284) overcorrected capacity: overflow
  fell below 1%, but HPWL exploded to 133.596M. Region partitioning must be
  coupled to continuous displacement/HPWL cost before it can be useful.
- All three methods were tested as candidate generators with unchanged exact
  HPWL/overlap auditing; none is accepted for the main route. The next global
  mechanism must screen total HPWL before committing a block/region move and
  operate earlier with topology-preserving relative coordinates.

## H295-H297 topology-preserving initialization and shared density directions

- Uniform affine and monotone power-law coordinate warps preserve local order
  but cannot allocate movable area around fixed macros. H295/H296 both fail
  their pre-registered overflow screens and are closed.
- H297 implements a direction-only heuristic: for each low-degree net, exact
  rectangle-overlap subgradients are averaged over movable pins and blended
  back into the per-node density direction. The exact max-minus-min HPWL and
  exact rectangle/bin overlap oracles are unchanged.
- On a fresh HPWL-only a1 seed, degree-16 full sharing improves the 501-step
  frontier from 91.0M/78.07% to 87.68M/72.59%; degree 64 reaches 71.88%
  overflow but with 88.77M HPWL. This confirms active-set starvation as a
  real mechanism, not merely a lambda-scale artifact.
- Unchecked continuation to 1001 steps reaches 152.4M/55.59%, so shared
  directions alone still create cumulative topology debt. Retain sharing
  only in an early continuation and combine it with a capacity-aware,
  fixed-macro-aware cluster assignment plus cumulative exact HPWL screening.

## H300-H306 guarded directions and transport audit

- HPWL-orthogonal projection of the shared density direction (H300) did not
  help: alpha 0.5 reached 87.08M/74.21% and alpha 1 reached 99.86M/77.55%,
  both worse than unprojected H297. The useful density direction is not
  separable from the HPWL subgradient by a simple projection.
- H304 made sharing explicitly early-only (first 250 iterations), then reset
  Adam and handed control to exact node-wise GP. It ended at 66.86M/83.83%,
  showing that takeover recovers HPWL but cannot invent missing capacity.
- H305's incident-net auction transport reduced overflow to 58.12% in one
  round but exploded HPWL to 366.77M. Local net screening misses cross-net
  accumulation.
- H306 added a per-round exact global HPWL rollback guard. At 5% and 10%
  budgets it rolled back every auction move, proving the H305 damage is not a
  small drift. Keep the guard as a safety primitive, but close individual
  auction transport. The open mechanism is a connected-cluster proposal with
  fixed-macro-aware capacity and exact global screening.

## H307-H308 local groups and identity-swap limit

- Weighted net-connected groups restricted to four-bin destinations (H307)
  still raised HPWL from 42.98M to 81.74M for only 4.85 overflow points; the
  5% exact global guard correctly rolled back the round.
- H308 establishes a structural limit: exchanging two identical rectangles
  leaves every rectangle-to-bin occupancy value unchanged exactly, so an
  equal-shape swap can never reduce exact overflow. H287 remains useful only
  as density-neutral HPWL polish. Any capacity-repair exchange must involve
  unequal shapes, empty capacity, or a translated connected cluster.

## H309-H311 unequal-shape and micro-batch capacity repair

- H309 enables unequal-shape exact exchanges into nearest underfull bins. It
  accepts 5,187 swaps and lowers overflow by 0.50 points for a 0.657M HPWL
  increase. This is topology-safe enough to retain, but not remotely enough
  for the 7% target.
- H310 limits weighted connected-cluster transport to 100 nodes under a 1%
  exact global HPWL guard; one batch gains only 0.0237 overflow points for
  0.074M HPWL.
- H311 chains ten such independently audited batches: overflow improves by
  0.259 points and HPWL rises 0.706M. The monotone safe behavior validates the
  guard, but thousands of batches would exhaust the HPWL budget. The next
  mechanism must transfer larger area hierarchically while preserving net
  coordinates inside each cluster.

## H312-H315 exact overlap coordinate descent and net-block targets (2026-08-24)

- Added an optional `density_coordinate` phase. It performs exact axis/breakpoint
  coordinate descent and accepts only exact overflow-decreasing moves under a
  cumulative HPWL budget. This is an active-set search phase; it does not alter
  the HPWL or rectangle/bin overlap objective.
- Fresh a1, 512x512, one node-coordinate sweep with a 10% HPWL budget moved 789
  nodes and changed overflow 99.9673% -> 99.8705% (HPWL 53.22M -> 53.41M).
  Three sweeps with a 20% budget reached 99.0417% but required 4,634 moves and
  54.91M HPWL. Single-node capacity transfer is too weak.
- Added rigid net-block coordinate descent. On a fresh 64x64 grid, one sweep
  with 5,000 blocks reduced overflow 99.7009% -> 97.2479% (HPWL 53.22M ->
  58.48M); three sweeps reached 95.9516% but raised HPWL to 63.87M. A
  capacity-target variant that ranks exact underfull bins was similar
  (97.1634% after one sweep), so direction targeting alone does not solve the
  topology/capacity conflict.
- A 512x512 block sweep improved only 0.21 percentage points, confirming that
  coarse resolution is necessary for area transfer. Passing the coarse result
  to fine GP caused density debt to reappear; a 500-step trajectory run ended
  at 77.87% overflow and 103.92M HPWL.
- Coordinate descent is therefore retained as a breakpoint escape primitive,
  but not as the main a1 route. The next mechanism still needs a true coarse
  capacity assignment with relative-coordinate clusters and a stronger global
  HPWL guard.

## H316-H319 coarse block bridge plus fine continuation (2026-08-24)

- Increasing the coarse block step to 16 bins reduced fresh 64x64 overflow to
  94.3770% at 79.76M HPWL; four exact HPWL-preserving recovery sweeps reduced
  this to 93.4910% and 66.16M. The bridge is useful as a topology-aware
  preconditioner but remains far from feasible.
- Fine 512x512 GP from that checkpoint with default density scale reached
  89.4969% overflow and 54.39M HPWL in 500 iterations. With density scale 5,
  it reached 78.9254% but 107.04M HPWL. Thus lambda strength trades capacity
  against topology in the expected way and cannot close the gap by tuning
  alone.
- Regional prices plus early net-sharing from the same checkpoint reached
  86.229% / 67.29M at 300 iterations, but the 500-step continuation ended at
  76.6178% / 111.13M. Early sharing is helpful only as an activator; retaining
  it or regional prices late creates the same topology debt.
- Adagrad with density scale 5 was inferior to Adam (300 steps: 97.96%
  overflow / 76.51M HPWL), so optimizer swapping alone is not sufficient.

## H320 dynamic coarse capacity assignment (2026-08-24)

- Added a dynamic coarse residual-capacity ledger. Fixed macro occupancy is
  subtracted exactly, and every accepted cluster move updates the residual
  capacity of all affected coarse bins. Destination candidates are therefore
  selected from capacity that has not already been consumed by earlier moves.
- On fresh a1 at 64x64, one 16-bin rigid-block sweep improved overflow to
  92.8733% at 79.84M HPWL, versus 94.3770% without the residual ledger. Two
  sweeps reached 89.8148% but raised HPWL to 106.33M; exact recovery reduced
  this to 88.49M/89.1231%.
- Fine default GP from the polished H320 checkpoint ended at 55.18M/88.6936%
  after 500 iterations. The coarse assignment is measurably better than the
  previous bridge but still does not preserve enough topology for the 7%
  target. The residual ledger is retained as a component of the next
  hierarchical cluster assignment, not as a complete solution.

## H321 hierarchical relative-coordinate clusters (2026-08-24)

- Corrected the coarse residual ledger to use capacity minus total occupancy;
  the former H320 expression accidentally removed fixed macro area from the
  occupancy term and could rank obstructed destinations as available.
- Added coarse-to-fine overloaded-source clusters grown through low-degree
  nets. Cluster translations preserve relative coordinates and every move is
  audited with exact overlap and all incident nets.
- One fresh a1 sweep with node fallback reached 57.3375M/97.3552% in 2.55s.
  Cluster-only isolated 158 accepted moves and reached 54.7245M/98.7623%.
- Three cluster-only sweeps reached 63.8652M/96.9369% in 136.37s, exhausting
  both the 20% HPWL budget and the a1 runtime gate. Sequential cluster commit
  is therefore closed.
- The remaining structural issue is cooperative cut cancellation: adjacent
  clusters that will move together pay temporary cut-net HPWL when screened
  separately. The next test is stage-level simultaneous cluster assignment
  with exact global batch audit and prefix backtracking.

## H322 cooperative exact cluster batches (2026-08-24)

- Added tentative destination-capacity reservation, simultaneous disjoint
  cluster placement, exact global HPWL/overlap audit, rollback, and sorted
  prefix backtracking.
- Sorting by raw exact density gain reached 63.8534M/98.1594% in 5.13s.
  Sorting by density gain per individual HPWL cost improved to
  63.8462M/97.1887% in 3.60s, but still missed H321's 96.9369% screen.
- Cooperative acceptance is closed as the primary mechanism. The dominant
  remaining defect is arbitrary BFS truncation: clusters have high hypergraph
  boundary cost before destination assignment. The next cluster generator
  must explicitly maximize heavy-edge internal connectivity.

## H323 heavy-edge cluster growth (2026-08-24)

- Replaced FIFO growth with accumulated `net_weight/(degree-1)` max-priority
  growth while retaining H322's exact cooperative audit.
- Fresh a1 reached 63.8060M/97.3314% in 3.81s, slightly worse than ordinary
  BFS. Better local internal connectivity does not create the global area
  transport required from the concentrated initialization.
- Close local rigid cluster translation. Reopen recursive global hypergraph
  capacity partitioning, whose prior H284 result proves capacity success, and
  target its main defect: HPWL-destructive leaf coordinate assignment.

## H324 exact-HPWL-guided bisection leaves (2026-08-24)

- Regenerated a fresh a1 HPWL-only checkpoint inside H324: 43.1844M/96.1513%
  in 3.86s. No historical placement file was used.
- Ordinary position-seeded bisection reached 207.101M/1.02232%. Replacing only
  leaf mapping with exact incident-net HPWL plus rectangle-level residual
  capacity reached 201.354M/0.295336% in 2.45s.
- The 2.78% HPWL improvement fails the locked 10% screen. Leaf mapping matters,
  but upper uncoarsened cell-by-cell partitioning dominates topology damage.
- The next global mechanism must add multilevel hypergraph coarsening and keep
  strongly connected movable cells atomic across upper capacity cuts.

## H325 accumulated shared-net atomic coarsening (2026-08-24)

- Corrected the H325 heavy-edge matcher to sort and reduce all shared-net pair
  contributions before mutual-best matching.  The prior implementation kept
  only the largest single net for a pair, producing almost no coarsening
  (113,037 groups from 210,904 cells).
- On H325's own fresh HPWL seed, the corrected depth-12 run reached
  185.968M HPWL, 0.262431% exact overflow, and 4.11s.  It formed 109,799
  groups with 101,105 merges and 114,563 atomic unit moves.
- This is a real but small 1.42% HPWL improvement over the previous H325
  depth-12 result and remains far below the preregistered 15% improvement
  screen over H324's 201.354M.  The remaining loss is therefore in upper
  recursive capacity cuts, not leaf anchors or release depth.
- H325 is closed.  The next mechanism will test stronger multilevel atomicity
  (larger supernodes and one extra matching round) under a new protocol before
  any parameter sweep.

## H326-H334 staged recovery and hierarchical capacity (2026-08-24)

- H326 stronger atomicity (16-cell groups, degree 24, four matching rounds)
  reached 184.580M/0.295225% in 4.37s; larger supernodes alone are not enough.
- Exact recovery is effective when allowed to spend the 7% slack: H327 from
  H326 reached 162.736M/6.99963% in 10.14s, and H328's ten sweeps reached
  153.616M/6.99999% in 35.54s.  Net-block recovery (H329-H330) reached
  150.344M/6.99999% in 48.30s, only a modest improvement.
- Relaxing to 15% (H331) restored topology to 150.891M/14.9999%; eight-step
  retightening stalled at 12.8905% and twenty more sweeps only reached
  12.6599% (145.228M).  A hierarchical exact coordinate repair (H333) then
  assigned residual capacity and reached 151.498M/0.0347% in 5.51s; a short
  post-recovery (H334) reached 144.024M/6.99982%.  This validates hierarchical
  capacity assignment but does not close the HPWL gap.
- H336 is the best new structural route: starting from strict H257, exact
  axis-ranked bisection reached 124.497M/0.2755%; ten exact net-block recovery
  sweeps (H337) reached 109.490M/6.99998% in 44.81s.  It remains above the
  paper target but is materially better than HPWL-seed bisection.
- Coarse leaf bisection (H339, 32-bin leaves) worsened HPWL to 133.073M, so
  leaf coarsening is closed.  The next mechanism should preserve H257 topology
  while using density-neutral, net-aware global exchanges or partial capacity
  cuts rather than repartitioning every cell.

## H340 net-aware equal-shape exchange (2026-08-24)

- Added an optional exact candidate neighborhood that enumerates same-shape
  nodes sharing low-degree nets before the spatial ring. Equal-shape swaps
  leave the exact occupancy field invariant, so this changes only the search.
- Five sweeps from H257 with 64 candidates and a 32-entry exact shortlist
  produced 86.691M/6.99996% in 90.49s, only 0.10% below H257 and not better
  than H287. Long-range density-neutral swaps are closed as the primary
  mechanism; area-changing topology-preserving assignment is required next.

## H341 HPWL-directed connected-block recovery (2026-08-24)

- Disabled node-wise density directions and searched only exact connected
  net-block moves (degree 16, block size 8) from H257 under the 7% cap.
- Ten sweeps accepted 113,772 block moves and reached 86.710M/7.00000% in
  28.83s, only 0.47% below H257.  Connected blocks alone do not escape the
  feasible HPWL basin; the missing mechanism must coordinate area transfer
  and topology in a joint global assignment.

## H342 relaxed net-block then capacity repair (2026-08-24)

- From H257, relaxing the exact cap to 15% for ten connected-block HPWL
  sweeps reached 84.761M/15.0%. The subsequent hierarchical capacity repair
  reduced overflow to 4.269% but raised HPWL to 93.240M.
- The two-stage relaxation therefore does not beat H257: the capacity repair
  spends more topology budget than the relaxed stage gains. Future methods
  must score area-changing destinations jointly with net HPWL, rather than
  applying a sequential cap relaxation.

## H343-H344 destination shortlist and cluster granularity (2026-08-24)

- Expanding the hierarchical destination shortlist from 8 to 64 (H343) left
  the result unchanged at 93.241M/4.217% overflow, showing that candidate
  truncation is not the bottleneck.
- Reducing rigid clusters from 64 to 16 nodes (H344) also left HPWL unchanged
  at 93.240M, although overflow reached 3.876%. Rigid translation itself,
  rather than cluster size or destination count, is the topology bottleneck.
- The next mechanism must allow elastic within-cluster assignment (different
  node displacements) while jointly auditing exact HPWL and capacity.

## H345-H346 elastic within-cluster assignment (2026-08-24)

- Implemented an exact elastic operator in which each node chooses a discrete
  breakpoint fraction independently, followed by one group overlap audit.
  The full H345 setting was computationally infeasible (>4 minutes) because
  per-node incident-net oracles multiplied the candidate cost.
- A bounded one-level profile (H346: 16-node clusters, 8 destinations, 500
  cluster budget) completed in 6.66s but raised HPWL 84.761M→92.959M and
  stopped at 7.064% overflow. The elastic fractions do not solve the topology
  tradeoff at this scale; close the mechanism unless a cheaper global oracle is
  designed.

## H347 coarse flow control (2026-08-24)

- Exact 64x64 coarse flow from H257 accepted 606 anchor moves and reduced
  exact overflow only from 6.99996% to 6.89671%, while HPWL rose to 87.382M
  (0.31%) in 0.41s. Individual-node flow is cheap but not topology-aware;
  connected-block commodities are required for any further flow design.

## H348 connected-block coarse-flow commodities (2026-08-24)

- Implemented deterministic low-degree net commodities (degree <=16, max 8
  cells) with exact group overlap and incident-net HPWL screening on each
  coarse flow edge.
- From H257, one 64x64 pass accepted 342 group moves and reached
  87.347M/6.90655% in 0.49s. This is slightly worse than the individual-flow
  control (87.382M/6.89671%) and does not improve H257. Single-step commodity
  transport is closed; the next variant will test exact multi-scale
  displacement for each commodity.

## H349 multi-scale exact commodity displacement (2026-08-24)

- H349 screened 1.0, 0.5, and 0.25 rigid fractions for each connected
  commodity using exact rectangle/bin overlap and exact incident-net HPWL.
- From H257 it committed 494 moves and reached 87.340467M HPWL / 6.889111%
  overflow in 0.541525s, versus H348's 87.347M / 6.90655% in 0.487668s.
- The extra scales provide only a 0.007% HPWL improvement over H348 and still
  leave a 0.221M HPWL debt relative to H257. Scalar displacement magnitude is
  therefore not the bottleneck; the coarse flow's diagonal direction is.
- Keep H349 as an opt-in primitive but reject it as the production route. The
  next hypothesis is to screen axis-separated `(dx,0)` and `(0,dy)` rigid
  commodity directions in addition to the diagonal.

## H350 axis-split exact commodity directions (2026-08-24)

- H350 screened diagonal, x-only, and y-only rigid commodity vectors at
  scales 1.0, 0.5, and 0.25 with exact overlap and HPWL audits.
- From H257 it reached 87.324700M HPWL / 6.88948% overflow in 0.6344s with
  521 moves. This improves H349 by only 15.8k (0.018%) and leaves the H257
  topology debt essentially unchanged.
- Axis decomposition is retained as a useful direction primitive, but scalar
  flow plus direction screening is not enough. H351 will add an exact
  non-positive incident-net HPWL guard to quantify the density gain available
  without paying topology debt.

## H351 exact HPWL-guarded axis flow (2026-08-24)

- H351 added a hard exact incident-net guard `delta_HPWL <= 0` to H350's
  diagonal/axis/multiscale commodity candidates.
- Only 12 moves survived; the fresh H257 continuation ended at 87.1195M HPWL
  and 6.99966% exact overflow in 0.589s, essentially unchanged from H257.
- Positive HPWL moves are therefore necessary for local flow to move area. The
  missing mechanism is joint global assignment whose cut-net deltas can cancel
  in a single exact batch, not another scalar or direction sweep.
- Close H351 as a solver and pivot to a jointly audited commodity assignment;
  keep H257/H287 as the current strict/quality records.

## H352 jointly audited commodity destination assignment (2026-08-24)

- H352 generated 26,235 exact direction/scale trials and audited 186 global
  prefixes of disjoint commodity proposals. It committed one 12-commodity
  batch, reaching 87.128095M HPWL / 6.994967% overflow in 0.826598s.
- The exact global audit is valid, but requiring every prefix to decrease
  overflow rejects the cooperative cancellation pattern the method was meant
  to expose. The monotone prefix rule, not the exact objective, is now the
  bottleneck.
- H352 is rejected as configured. H353 will keep a bounded non-monotone
  candidate reservoir and choose only the best final exact prefix.

## H353-H354 deferred joint reservoir ablations (2026-08-24)

- H353 allowed non-monotone intermediate prefixes, but the per-edge top-8
  reservoir still produced only 12 disjoint commodities: 87.120209M /
  6.997823% in 0.611627s.
- H354 widened the reservoir to 64 per edge and obtained the same 12-group
  batch: 87.121849M / 6.996061% in 0.619557s. Candidate ranking is not the
  bottleneck; the earlier per-edge `group area <= edge.area` filter discards
  larger connected commodities.
- Next remove that fragment filter and enforce total source-area reservation
  only at the exact global batch stage.

## H355 joint assignment without edge-fragment area filter (2026-08-24)

- Removing the per-edge group-area rejection still yielded only 13 accepted
  groups (12 effective prefix groups): 87.121849M HPWL / 6.996061% overflow in
  0.608213s.
- The unchanged result rules out edge-fragment area as the primary bottleneck;
  repeated high-score groups dominate each edge's proposal pool. A very wide
  per-edge bid pool is the next diagnostic before closing this assignment
  family.

## H356 very-wide joint bid pool (2026-08-24)

- Increasing exact bids per edge from 64 to 1,024 still selected the same
  12-group prefix: 87.121849M HPWL / 6.996061% overflow in 0.616319s.
- Proposal-count coverage is not the bottleneck. The overflow-gain/HPWL ratio
  favors a tiny prefix; the next test uses an explicit exact HPWL trust-region
  budget and maximizes overflow decrease within it.

## H357 exact HPWL trust-region joint batch (2026-08-24)

- A 1% global exact HPWL budget selected 64 disjoint commodities instead of
  H356's 12, confirming that ratio ranking was conservative.
- The larger batch reached 87.1403M HPWL / 6.96741% overflow in 0.594443s,
  worse overflow than H350's 6.88948% and still far above the paper HPWL.
- Joint coarse-flow assignment is therefore closed as a post-H257 route. The
  dominant loss is inherited before feasibility; the next mechanism must alter
  raw-to-feasible initialization/continuation while retaining exact audits.

## H358 early-share seed capacity-flow test (2026-08-24)

- The H304 staged seed (66.8702M HPWL / 83.8347% overflow) had 432 coarse
  edges, but all 110,592 exact-anchor moves were rejected; output stayed
  66.8702M / 83.8347% in 0.688615s.
- This confirms a flat exact capacity basin around the low-HPWL share seed.
  Local coarse flow cannot bridge it; raw-to-feasible continuation needs a
  coordinated nonlocal assignment before fine-grid recovery.

## H359 nonlocal commodity jump from early-share seed (2026-08-24)

- Direct source-to-destination connected commodity jumps crossed the H358 flat
  basin: 19,733 moves reduced exact overflow 83.8347% -> 66.3385%.
- HPWL exploded 66.8702M -> 244.359M and runtime was 28.6282s. Nonlocal flow
  is capacity-effective but topology-destructive; direct jumps are closed.
- The remaining raw-to-feasible problem is coordinated large-group transport
  with relative-net geometry preservation.

## H361-H362 recursive share-seed assignment and recovery (2026-08-24)

- H361 position-seeded atomic recursive bisection from H304 reduced overflow
  83.8347% -> 0.250171% in 5.761s, with HPWL 66.8702M -> 155.243M.
- H362 exact breakpoint/compact/net-block recovery spent the 7% slack and
  reached 124.767M/6.99988% in 39.339s (726,870 moves), a 19.6% HPWL recovery.
- This is materially better than direct nonlocal flow and remains within the
  a1 runtime envelope. Continue recovery to test whether the topology debt can
  be reduced below 110M before closing or transferring the route.

## H363-H367 feasible-branch continuation and mapping ablations (2026-08-24)

- H363 extended the concentrated-share branch by twenty exact recovery sweeps,
  improving 124.767M to 123.236M at 7.00000% overflow in 77.08s. The small
  gain closes that branch as a route to the paper target.
- H364 replayed recursive assignment with alternative exact leaf maps. The
  H361 axis-rank map remained best (155.243M); curve-rank and axis+curve maps
  reached 249.705M and 219.281M, respectively. Curve ordering breaks local
  topology and is rejected.
- H365 continued H337's feasible-topology recovery with net-block density
  directions and contraction for twenty sweeps: 109.490M -> 104.656M at
  7.00000% overflow, 129.94s. This is the best a1 quality but misses the
  strict runtime gate when appended to the full staged prefix.
- H366 lowered the bisection degree limit to 16, but position-seeded cuts do
  not use degree-limited connectivity ordering; it exactly reproduced H336
  (124.497M/0.2755%). This parameter is not a lever.
- H367 used five H337 recovery sweeps and reached 105.295M/6.99976% in
  31.99s. Including the H336/H337 prefix it is runtime-valid for a1 and is
  the current strict quality record. The remaining 43.8% HPWL gap is incurred
  before local recovery, so further scalar recovery tuning is unlikely to
  meet the paper target; a new topology-preserving upper-cut assignment is
  required before a4 transfer.
- H368 fed the H340 exact exchange seed through H336-style capacity cuts and
  five recovery sweeps. It reached 106.652M/7.00000% in 30.53s, worse than
  H367; the small pre-bisection topology gain is erased by upper cuts.
- H369 tested exact coordinate interpolation between H257 and H336. The best
  alpha=0.25 started at 93.645M/31.06% and recovered to 87.360M/20.60%, but
  plateaued above the 7% cap even with twenty sweeps and larger net blocks.
  Partial interpolation preserves HPWL topology yet lacks a capacity bridge;
  hierarchical assignment must move area across cuts while retaining net
  geometry.
- H370 relaxed recursive cut balance from 2% to 15% (and an exploratory 49%).
  The bisection checkpoint changed only to 124.001M/0.37% overflow, and five
  exact recovery sweeps ended at 106.665M/7.00000%, worse than H367. Cut
  imbalance tolerance is not the topology bottleneck.
- H371 implemented surplus-only hierarchical cuts: position-induced sides are
  preserved and only boundary surplus violating exact side capacities moves.
  From H257 this reduced the checkpoint to 114.146M/2.16% overflow and the
  five-sweep recovery to 102.008M/7.00000%.
- H372 removed the axis pre-map from H371. The checkpoint became 93.7055M/
  2.31%, and exact recovery reached 86.4503M in five sweeps (30.59s),
  86.1896M in ten, and 86.1566M in twenty. This is a major topology-preserving
  gain and the current strict a1 record, but remains 17.7% above paper HPWL.
- H373 attempted raw a2 transfer. A 500-step exact HPWL-only seed reached
  50.7285M/95.86% overflow; surplus-only cuts reduced overflow to 1.36% but
  raised HPWL to 191.591M. Without a feasible topology seed, the mechanism
  does not transfer, so a2-a4 claims remain deferred.
- H374/H375 exact net-aware equal-shape exchange polished H372: 86.4503M to
  86.2267M after five-sweep recovery, and 86.1896M to 85.9993M after
  ten-sweep recovery. H375 is the current strict a1 record at 7% overflow;
  exchange gains are modest and the remaining gap is 17.4%.
- H376 temporarily relaxed the exact overflow cap to 15%, lowering H375 to
  83.8203M, but direct exact retightening remained at 14.87% overflow. H377
  recutting that relaxed placement recovered to 86.1054M/6.99992%, worse than
  H375. Relaxation is useful for topology exploration but its feasible basin
  is disconnected; retain H375 as the production branch.
- H378 tested a raw a2 exact joint GP seed (1000 iterations, density scale 5)
  without electrostatics. Exact overlap gradients remained effectively flat:
  overflow was 97.7351% and HPWL 120.171M after 22.19s. A2 still requires a
  dedicated raw-to-capacity continuation before H372 can transfer.
- H379 applied exact recovery directly to the raw a2 HPWL seed with a loose
  99% cap; overflow slightly worsened to 96.1504%, confirming that the flat
  exact active set is not caused solely by the 7% constraint.
- H380 widened raw initialization (sigma 0.1), but HPWL-only descent still
  collapsed to 59.835M/98.149% overflow. Spatial spread without capacity
  prices is not a viable a2 seed.
- H381 affine-expanded the raw a2 HPWL seed (scales 2 and 3) before
  surplus-only cuts. Scale 2 reached 277.768M/1.654% and scale 3 reached
  223.655M/1.617%; expansion does not preserve enough net topology for
  capacity assignment.
- H382 enabled atomic net coarsening during a2 surplus cuts, but HPWL rose to
  253.393M/0.302%, worse than non-atomic H373. Atomic grouping cannot rescue
  the collapsed raw seed; a2 needs a fundamentally different capacity-aware
  initialization.

## H360 larger net-connected nonlocal groups (2026-08-24)

- Raising commodity limits to degree <=64/max 64 still accepted 18,249 direct
  jumps and reduced overflow only to 66.9037%.
- HPWL rose 66.8702M -> 256.666M in 23.372s, worse than H359. Larger rigid
  groups do not cure topology damage; the full source-to-sink jump is the
  defect. Close rigid nonlocal jumps and pivot to recursive capacity cuts or a
topology-aware raw initialization.

## H383-H388 raw capacity-aware continuation (2026-08-24)

The raw a2 experiments isolate the remaining failure mode.  Direct
surplus-only recursive cuts preserve exact capacity but destroy HPWL topology
(H383: 79.8M -> 2,439M).  A short HPWL warmup plus local leaf restriction
avoids the billion-scale jump but leaves 28--67% overflow (H385).  Exact
hierarchical cluster coordinate descent is qualitatively better: it crosses
the zero-gradient basin with only a 20% HPWL increase and reaches 98.53%
overflow (H386).  Extending it does not help (H388: 94.34% overflow and
159.7M HPWL).  Thus coordinate descent is a useful activation operator, not a
complete raw-to-feasible initializer.  The next design must interleave
coarse capacity assignment with these local net-aware clusters, moving
neighboring groups across adjacent bins and auditing exact HPWL at each batch.

H389 confirms elastic/heavy-edge cluster growth still leaves 95.32% overflow
and costs 159.5s. H390's shallow exact cut reaches 0.153% overflow only with
1,822M HPWL and 259.4s from global leaf scans. The next mechanism must bound
candidate destinations to adjacent bins and preserve net geometry batchwise.

H391--H395 clarify the locality trade-off. Adjacent-bin clusters preserve
HPWL but cannot move enough area (98.5--98.9% overflow); node moves and
radius-one rigid groups are effectively frozen, while random jitter destroys
HPWL. The next bridge must allow a small number of multi-bin transfers under
an exact global HPWL trust region, with intermediate capacity targets rather
than a one-shot 7% assignment.

H396 confirms radius two remains insufficient (95.62% overflow). H397 shows
that chaining the coordinate checkpoint into recursive cuts causes a 1,300M
HPWL topology jump. The required method is one interleaved multi-bin bridge
with intermediate exact capacity targets and a global exact HPWL trust region.

H398--H399 tested exact nonmonotone overflow bands of 2 and 20 percentage
points; both remained above 97% overflow. Allowing temporary objective
increase is not enough because the candidate graph itself lacks a feasible
multi-bin path. The next mechanism must construct explicit staged target
capacities and coordinated net-group assignments in one bridge.

H400--H401 provide the first partial bridge: coarse exact staged assignment
at 64 bins lowers overflow to 86.98% and 83.70%, but HPWL rises to 223--239M.
This branch is useful only as an intermediate checkpoint; exact fine
refinement must recover topology without undoing the occupancy gain.

H403 isolates the H375 recovery components from the identical H372 checkpoint.
Node moves provide a small baseline gain; breakpoint oracle and compact
directions add only marginal improvements, while net-block density directions
are the strongest isolated recovery component. Contraction is not consistently
helpful, and equal-shape exchange is a small but expensive final polish.
H404 confirms the proposed dead-zone escape can cross a flat HPWL basin from
raw initialization, but leaves overflow at 99.68%. Thus boundary-active search
is useful as an HPWL recovery operator, not as a standalone density solution;
it must be coupled to staged exact capacity targets.

## H405-H410 isolated zero-gradient mechanism ablation (2026-08-24)

To correct the previous combined test, six mechanisms were run separately from
the identical raw a1 center-Gaussian seed (seed 123, sigma 0.001), with exact
overlap and no smoothing, electrostatics, legalization, or prior placement
checkpoint:

| experiment | only enabled mechanism | HPWL (M) | overflow | time (s) |
|---|---|---:|---:|---:|
| H405 | breakpoint oracle | 48.0784 | 99.7227% | 8.67 |
| H406 | 0.25-bin boundary strip | 48.2662 | 99.6457% | 6.80 |
| H407 | zero-gradient dead-zone escape | 50.3151 | 99.6406% | 3.08 |
| H408 | net-aware shared density direction (100 GP iters) | 43.7691 | 97.0140% | 3.26 |
| H409 | node-wise exact finite directional oracle | 53.2191 | 99.4508% | 6.89 |
| H410 | controlled deterministic exact noise | 50.3109 | 99.6430% | 3.06 |

H405 is the strongest isolated HPWL recovery operator; H408 is the strongest
raw-seed HPWL aggregator but still leaves density nearly infeasible. H409 is
the only isolated mechanism whose primary effect is overlap reduction, and it
does so only by 0.52 percentage points. H407 and H410 are nearly identical,
showing that dead-zone escape and controlled exact noise are interchangeable at
this scale. The causal conclusion is that zero-gradient fixes are useful
search directions, not a capacity bridge: a topology-preserving staged exact
capacity assignment must be interleaved before these operators can approach
the 7% target.
