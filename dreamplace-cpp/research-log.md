# Research Log

## 2026-08-06 Bootstrap

- Locked the standalone implementation boundary: no source or binary dependency
  on `alg-electronic`.
- Locked raw Bookshelf input.  Existing `.eplace-ip.pl`, `.eplace.gp.pl`, and
  `.nsp.pl` files are forbidden as placement inputs.
- Paper Table II targets are 73.22M for adaptec1 and 82.22M for adaptec2 after
  the complete GP/LG/DP flow.  These are distinct from pre-legal values quoted
  by later papers.
- Official repository configurations use weighted-average wirelength,
  Nesterov, target density 1.0, stop overflow 0.07, seed 1000, fillers, and
  512x512 / 1024x1024 bins.
- H0 is an end-to-end correctness experiment before quality tuning.

## 2026-08-06 H0 Validation

- Corrected BB step estimation to reevaluate both reference gradients under
  the current lambda and gamma, matching upstream `step_bb` semantics.
- Separated objective density (movables plus fillers) from overflow density
  (movables only), matching `ElectricPotential.forward(mode="overflow")`.
- At 64x64, adaptec1 reached 6.9178% overflow by iteration 546, confirming the
  expected high-overflow plateau followed by rapid electrostatic spreading.
- At 512x512, adaptec1 selected iteration 595: GP 71.235743M / 6.9654%.
- Width-first legalization produced 123.991762M and was rejected.  Preserving
  left-to-right GP order reduced the same checkpoint to 74.819421M.
- At 1024x1024, adaptec2 selected iteration 614: GP 80.889298M / 6.9302%,
  final legal 84.426599M.
- Both final placements passed the ISPD `legal2.pl` script with zero Type
  0/1/2/3 errors.

## 2026-08-06 Nonsmooth Migration

- Added an exact pin-offset HPWL oracle with tied-extrema subgradient sharing.
- Added an opt-in solver layer for Heavy Ball, Adam, AMSGrad and AdaGrad, three
  lambda controllers, and delayed HPWL-only bundle aggregation.  Smooth WA +
  BB-Nesterov remains the default path.
- A 64x64/800-step Adam run with the initial dimensionless lambda cap of 1e6
  stalled at 13.315% overflow.  Raising the numerical cap to 1e12 allowed the
  same seed and parameters to select iteration 521 at 90.069M/6.952%.
- The selected checkpoint legalized to 114.017M and passed `legal2.pl` with
  Type 0/1/2/3 all zero.
- Trajectory lambda plus four delayed HPWL cuts selected 95.187M/6.824% and
  legalized to 118.153M.  Bundle integration is operational but did not
  improve quality in this experiment.

## 2026-08-07 Smooth Warm Start And Exact Bundle Refinement

- Added an opt-in `--initial-pl` path for placements produced by this program.
  It loads movable cells only, preserves raw fixed objects, requires a complete
  movable-node set, and disables the exact solver's additional random noise.
- Started from the validated smooth adaptec1 checkpoint at
  71.235743M/6.9654% and ran 200 exact-HPWL iterations with 12 cuts,
  serious/null acceptance, a 0.01 feasible-refinement step scale, and
  legal-aware checkpoint selection.
- The best exact GP state reached 71.220005M/6.9990%.  The legal-aware state
  selected at iteration 10 was 71.220652M/6.9976% and finalized at 74.549788M.
- A matched direct-legalization control, using the same enhanced detailed
  placement but no effective global refinement, finalized at 74.548626M.
  Therefore the current bundle reduced pre-legal HPWL by about 15.1K but did
  not improve final legal HPWL; the 1.16K difference is slightly negative.
- Only six trial steps were accepted before the feasible guard rejected trials
  near 7.0000% overflow.  This identifies filter/hysteresis design and
  legal-aware trial acceptance as the next bundle issues, rather than a need
  for more iterations with the same policy.

## 2026-08-07 Grouped And Fixed-Topology Bundle Study

- Added default-off net-group HPWL bundle histories.  Grouped exact
  subgradients sum exactly to the original oracle, and nearby sampled exact
  subgradients contribute only to trial directions.
- Added trial-band boundary bisection.  A trial that leaves a strictly
  feasible center but remains within the configured tolerance is projected
  back to the strict 7% boundary before exact serious/null evaluation.
- Added a default-off post-legal x-only bundle.  Row assignment, obstacle-free
  segment ownership, and cell order remain fixed; every trial passes through
  segment Abacus and a full legality check.
- The grouped pre-legal stage improved GP by 16.098K but worsened a matched
  final legal result by 9.103K, so it was not selected for the best pipeline.
- The post-legal bundle reached 74.425850M.  The matched five-round no-bundle
  control reached 74.520409M, giving a 94.559K (0.127%) net improvement.
- A command with all new options omitted reproduced the prior 74.548626M
  three-round result exactly, confirming rollback compatibility.

## 2026-08-08 Progressive Exact-HPWL Legalization

- Added the default-off `--progressive-legalization` path. It preserves the
  original command and CSV schemas unless explicitly enabled.
- Stage 1 accepts only non-increasing exact pin-offset HPWL trials. Stage 2
  combines the exact subgradient with the existing electrostatic grid model.
  Stage 3 ramps a continuous fixed-macro escape potential, row proximity, and
  obstacle-free segment-capacity directions; it performs no row/site snap or
  cell-order assignment.
- A hard 80M HPWL cap locked the layout at 81.9% overflow. A dynamic cap near
  110M still locked it at 69.3%. Both policies were rejected; the HPWL cap is
  now applied only near the target overflow band.
- Pure band-dual control at 512x512 reached 101.277M/7.266% and 106.415M
  legal. Switching to the DREAMPlace update above 11% overflow and band-dual
  control below 11% improved the 800-step result to 83.701M/7.530% and
  89.237M legal.
- Adding 200 nominal refinement iterations produced no further accepted step
  after iteration 799. The best state is 83.254M/7.744% and 88.895M legal.
  The target 75-80M GP range was not reached.
- A combined gradient-sampling/serious-bundle run was slightly worse at
  83.577M/7.741% and 89.186M legal. Null steps eventually reduced the trial
  scale to 0.01 and the predicted decrease to zero.
- Continuous segment guidance reduced the 64x64 segment-overflow proxy from
  38.89% to 30.89% and improved legal HPWL from 125.489M to 124.650M, despite
  a small GP regression. This supports retaining it as a legalizability term.
- The final 88.895M placement passed the external ISPD `legal2.pl` checker
  with Type 0/1/2/3 = 0/0/0/0.
