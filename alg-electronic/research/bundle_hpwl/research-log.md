# Bundle HPWL research log

## 2026-07-31 - Bootstrap

- Locked metric: lowest exact HPWL with fine-grid overflow at or below 6%.
- Historical 3000-step baseline: 79,422,944 HPWL, 0.059985 overflow.
- Chosen first method: a small proximal cutting-plane bundle for exact HPWL,
  solved in the cut-weight simplex with Frank-Wolfe line search.
- Bundle is applied only in fine refinement; density model and lambda feedback
  remain unchanged.
- Git protocol cannot be used because this workspace is not exposed as a valid
  Git worktree to the current environment. The protocol is locked in files
  before experiments instead.

## 2026-07-31 - H1 confirmatory results

- The unchanged 800-step phase schedule is not yet feasible: baseline ended at
  51.842M / 39.54% overflow. This invalidated 800 steps as a feasible-ranking
  proxy, but it remains useful as a divergence screen.
- B1 (size 4, interval 5, prox 1, current mix 0.35): 64.397M / 34.82% at 800.
  Newest-cut weight stayed near 0.95, so little true aggregation occurred.
- B2 (size 4, interval 5, prox 2, current mix 0.50): 56.825M / 31.72% at 800.
  Newest-cut weight was about 0.33 and aggregate norm about 0.41 of the current
  subgradient norm.
- At 1800 steps B2 achieved 83.009M / 5.9956%, versus the contemporaneous
  baseline 90.117M / 5.9997%.
- At 3000 steps B2 achieved 76.522M / 5.9999%, improving the historical
  79.423M baseline by 3.65%.

## Next locked ablation

- H2a: prox 2, current mix 0.35 tests stronger history while keeping the same
  bundle geometry.
- H2b: prox 3, current mix 0.50 tests a stronger minimum-norm bundle solution
  without reducing the explicit current-gradient fraction.
- Screen both at 1800 steps and promote only the better feasible result.

## 2026-07-31 - H2 bundle-strength ablation

- H2a (prox 2, current mix 0.35) reached 81.523M / 5.9991% at 1800,
  outperforming B2 at that budget. At 3000 it regressed to 77.768M / 6.0000%,
  worse than B2's 76.522M.
- H2b (prox 3, current mix 0.50) reached 82.779M / 5.9706% at 1800. It
  over-regularized the HPWL direction and retained excess density margin.
- The 1800-step ranking is not a reliable substitute for full 3000-step
  validation. It remains a useful rejection screen only.
- Next: hold prox=2/current-mix=0.50 and test bundle sizes 3 and 6. This isolates
  how many distinct extrema active sets are useful.

## 2026-07-31 - H3 bundle capacity

- Size 3: 81.821M / 5.4368% at 1800 and 76.164M / 5.9999% at 3000.
  This is the new best, 4.10% below the historical exact-HPWL baseline.
- Size 6: 83.318M / 6.1460% at 1800, with aggregate/current norm ratio 0.34.
  Too many cuts over-regularized the HPWL direction and delayed feasibility.
- Next: size 3 with intervals 3 and 7 to isolate the temporal span of extrema
  active-set memory. Other controls remain fixed.

## 2026-07-31 - H4 extrema-driven cut timing

- Size 3 / interval 3: 82.981M / 5.9803% at 1800. Cuts were too similar and
  the low aggregate norm did not improve the feasible direction.
- Size 3 / interval 7: 83.349M / 5.9972% at 1800. Older active sets were stale.
- Adaptive trigger with min age 2 measured mean cosine -0.454 between current
  exact HPWL subgradient and the newest two-step-old cut; every eligible step
  triggered, with minimum cosine -0.903.
- This reveals a strong extrema-driven two-cycle. Next test: retain only two
  cuts at interval 5 to represent the alternating active sets directly.

## 2026-07-31 - H5 two-state bundle

- Size 2 reached 80.383M / 4.6955% at 1800, the lowest mid-run HPWL but with
  substantial density slack.
- At 3000 it reached 76.120M / 5.9968%, narrowly better than size 3's 76.164M.
  The 0.06% gap is below known run-to-run variation and is not yet significant.
- The robust conclusion is that 2-3 cuts are sufficient; the exact best size
  needs replication.
- Next: repeat size 2, then test a lower lambda increase cap because size 2
  over-spreads in the middle of the run.

## 2026-07-31 - H6 lambda coupling and extrema-aware mix protocol

- Size-2 repeat improved to 75.677M / 5.9982%, confirming a stable large gain
  over the non-bundle baseline but also showing about 0.44M run variation.
- Lambda-up 1.10 was still at 11.93% overflow after 1800 steps. Lambda-up 1.13
  reached 6.873% at 1800 but only 76.332M / 5.9975% at 3000. The default 1.15
  remains better.
- Next confirmatory method: extrema-aware current mix. Use
  mix=clip(0.50+0.30*cosine, 0.20, 0.80), where cosine compares the current
  exact HPWL subgradient with the newest cut. Conflict increases bundle use;
  alignment increases current-extrema use. All other best settings are fixed.

## 2026-07-31 - H7 extrema-aware mix result and decay protocol

- Full-run extrema-aware mix reached 76.920M / 5.9998%, worse than fixed mix.
  The insertion-time cosine remained near -0.85 late in the run, holding the
  effective current mix near 0.246 and preventing final HPWL recovery.
- Next: use the same cosine response early, linearly decay its gain to zero by
  fine step 1400, then finish with the validated fixed mix 0.50.

## 2026-07-31 - H7 result and outer-loop synthesis

- Decayed extrema mix reached 76.146M / 5.9996%. It recovered the fixed-mix
  performance range but did not beat the 75.677M best repeat.
- Conclude this research cycle: the supported contribution is a 2-3 cut
  proximal aggregate bundle for exact HPWL. Direct cosine control is useful as
  a diagnostic of extrema switching but not as a persistent controller.
- Default command behavior remains unchanged. Bundle is opt-in until broader
  benchmark and multi-repeat validation are available.

## 2026-08-01 - H8 relaxed-overflow cycle bootstrap

- New locked metric: lowest exact HPWL with fine-grid overflow <= 0.065.
- The existing CLI exposed the guarded recovery target but hard-coded the
  controller normalization and best-feasible checkpoint limit to 0.060.
- Add an optional `--overflow-baseline` parameter with unchanged 0.060 default.
- H8a is confirmatory: use the validated size-2 bundle and set limit=0.065,
  guarded recovery target=0.0648, interval=25. Other controls stay fixed.
- Subsequent changes will be based on whether the late trajectory wastes
  overflow headroom or repeatedly crosses the relaxed boundary.

## 2026-08-01 - H8a result and H8b role separation

- H8a reached 75.881M / 0.064991, slightly worse than the old 75.677M best.
- Its late mean effective lambda was 4.10e-5 versus 3.43e-5 in the old best.
  Raising the shared overflow baseline changed early ratio feedback and
  strengthened the eventual density pressure, so relaxed feasibility was not
  isolated.
- H8b adds a separate optional `--overflow-limit`. Controller normalization
  remains 0.060, the best-state limit is 0.065, and the fine recovery target is
  0.0648. This preserves the validated early trajectory while exposing only
  the extra final capacity.

## 2026-08-01 - H8b result and H9a late-step protocol

- H8b reached 75.762M / 0.064799. It is 0.084M above the prior best, well
  inside known run-to-run variation.
- Overflow was stable while HPWL fell by 1.533M over the last 400 fine steps;
  the final learning rate cooled from 8.81 to 2.14.
- H9a keeps the H8b lambda target/limit and bundle fixed, changing only the
  final cooldown floor from 0.15 to 0.25. This tests whether the 3000-step
  budget ends before bundle-guided HPWL recovery is complete.

## 2026-08-01 - H9a result and H9b cooldown expansion

- H9a reached 75.599M / 0.064875, the new best by about 0.078M.
- The lowest late checkpoints recur every five iterations, synchronized with
  bundle cut insertion. The best remains at fine iteration 2596, while the
  final logged iterate rebounds, so stronger late steps can help only if the
  checkpoint mechanism captures the nonsmooth cycle minimum.
- H9b changes only the final cooldown floor to 0.35. A positive result will be
  followed by interval contraction rather than an unbounded floor increase.

## 2026-08-01 - H9b rejection and H10a boundary controller

- H9b regressed to 76.298M / 0.064790. A 0.35 cooldown floor amplifies the
  exact-HPWL bundle cycle more than it accelerates recovery.
- Return to floor 0.25. H10a raises the recovery target from 0.0648 to 0.0649
  and shortens fine feedback from 25 to 10 steps, testing a higher-bandwidth
  controller that spends more of the 0.065 feasibility budget.

## 2026-08-01 - H10a rejection and H10b target isolation

- H10a regressed sharply to 77.675M / 0.064918. Lambda rose to about 7.45,
  versus 5.16 for H9a. The controller gains are per update, so shortening the
  interval without rescaling gains multiplied the effective closed-loop gain.
- H10b restores interval 25 and floor 0.25, changing only the target to 0.0649.
  This is the clean test of whether the remaining static headroom has value.

## 2026-08-01 - H10b result and H11 confirmatory repeat

- H10b reached 75.757M / 0.064907, worse than H9a and effectively tied with
  H8b. Spending the last 0.0001 target headroom is not a stable lever.
- H11 repeats H9a unchanged to test whether its 75.599M result survives known
  OpenMP/run-to-run variation. This repeat is confirmatory.

## 2026-08-01 - H11 confirms floor 0.25; H12a interval contraction

- H11 reached 75.638M / 0.064801. Together with H9a at 75.599M / 0.064875,
  the floor-0.25 repeat range is only 0.039M and both beat the old 6% best.
- H12a tests floor 0.28 with target 0.0648. This contracts the interval between
  the confirmed 0.25 point and rejected 0.35 point.

## 2026-08-01 - H12a result and H13 current-cut test

- H12a reached 75.685M / 0.064793, worse than both floor-0.25 repeats. Retain
  floor 0.25.
- H13 changes only bundle current mix from 0.50 to 0.55 at the confirmed lambda
  and cooldown settings, testing whether relaxed density permits a more direct
  exact-HPWL direction.

## 2026-08-01 - H13 result and relaxed-overflow outer-loop synthesis

- H13 reached 75.745M / 0.064833; current mix 0.55 is worse than the replicated
  mix-0.50 range and is rejected.
- Conclude the cycle with the H9a/H11 configuration: controller baseline 0.060,
  feasibility limit 0.065, recovery target 0.0648, guard interval 25, cooldown
  floor 0.25, and the validated two-cut bundle.
- Independent results are 75.599M / 0.064875 and 75.638M / 0.064801. The best
  improves the prior bundle result by 0.104%; both runs improve it, with only a
  0.039M repeat range.
- Relaxing overflow alone is not a large lever. The meaningful mechanism is
  preserving the early 6% controller trajectory while spending relaxed capacity
  during fine recovery, then preventing the learning rate from cooling too fast.
- All nine 3000-step outcomes are stored in
  `research/bundle_hpwl/data/relaxed_overflow_results.csv`.
# 2026-08-02: ISPD2005 eight-benchmark breadth run

- Locked a one-run-per-benchmark protocol before execution.
- Verified the DREAMPlace paper does not use its Section III-F 20%/10% values
  as final GP feasibility thresholds. The official repository at commit
  `6627f3327e6cc17db7782c0b90073a498531ca3c` defines default
  `stop_overflow=0.1`.
- Used 10% only for selected-state feasibility. Kept the validated 6% lambda
  baseline and 6.48% recovery target unchanged.
- Completed all eight 3000-step runs successfully in 2506 seconds wall time.
- All selected states met 10%; HPWL quality generalized poorly on instances
  lacking ePlace-IP initialization, while density feasibility generalized.
- Generated 8 GIFs (61 frames each), 8 per-benchmark convergence figures, two
  aggregate figures, raw result CSV, Markdown index, and PDF report.
- Classified this as exploratory breadth evidence because initialization is
  inconsistent and there is one run per benchmark.

## 2026-08-03/04 - adaptec2--4 initialization, lambda, and bundle cycle

- Locked a new protocol before code changes and runs. Primary metric is lowest
  exact HPWL at this solver's overflow <=10%; DREAMPlace Table II is contextual.
- Target alignment alone improved supplied-`.pl` results to 110.524M,
  305.871M, and 304.352M. This ruled out the 6.48% target as the sole cause.
- Added opt-in DREAMPlace-style center initialization and gradient diagnostics.
  An 800-step screen showed a large positive effect and was promoted.
- Center + guarded bundle reached 89.074M, 243.180M, and 206.646M. Paired
  no-bundle runs regressed to 92.213M, 258.925M, and 223.772M.
- On adaptec2, cooldown 0.50 and current mix 0.65 regressed; trajectory PI-D
  reached feasibility earlier and improved modestly. Reallocating coarse and
  medium iterations to fine refinement gave another 0.20M.
- A final narrow bundle ablation rejected current mix 0.35 and size 3. Prox 3
  improved by 0.10M on adaptec2 and transferred small gains to adaptec3/4.
- Final best results are 88.367M/9.9496%, 241.559M/9.9547%, and
  203.894M/9.9602%. The DREAMPlace reference level was not reached.
- Electric-density, baseline-lambda, and trajectory-lambda tests passed. The
  exact HPWL implementation remains max-minus-min with no smoothing operators.

## 2026-08-04 - conventional gradient and lambda cycle

- Locked the no-bundle protocol before adding options or running experiments.
- Compared plain normalized subgradient, heavy-ball, local adaptive restart,
  Adam, guarded lambda, trajectory PI-D and projected dual ascent on adaptec2.
- Added backward-compatible Adam parameter flags and an opt-in projected dual
  controller with an explicit unit test. Old invocation behavior is unchanged.
- Tuned Adam beta1 to 0.90 and found a step-scale plateau at 3--5. A controlled
  3000-step a2 run improves heavy-ball 93.075M to Adam 87.162M.
- Projected dual rate 1 was infeasible; rates 3--8 were feasible but all worse
  than trajectory. Guarded feedback overspread to 9.36% at the screen budget.
- Promoted Adam to adaptec3/4. Final best results are 87.162M/9.9536%,
  230.189M/9.9492% and 201.272M/9.9143%, all below the prior bundle values.
- Repeated a2 scale 5 at 87.227M, within 0.065M of the first run. Recorded all
  screen and promoted results, including negative runs, under
  `experiments/gradient_lambda_optimization`.
- Generated five vector/raster figures and numerical summary CSV files. The
  cycle concludes with tuned Adam + trajectory lambda as the current best
  exact-HPWL non-bundle method; DREAMPlace reference values remain unreached.

## 2026-08-04 - seven-percent adaptec1/adaptec2 cycle

- Locked the formal feasibility limit at 0.0700 using the pre-legalization
  DREAMPlace values in Table 5 of the stochastic-subgradient paper:
  70.3M/6.92% for adaptec1 and 79.3M/6.89% for adaptec2.
- Kept 150 exact-HPWL-only steps with lambda zero, removed coarse/medium
  stages, and used a 600-step electrostatic-density trajectory.
- Screened Adam/AMSGrad, trajectory targets 0.0695/0.0698, predictive lambda
  look-ahead, and adaptive step scales. No HPWL smoothing was introduced.
- Predictive lambda reduced a1 overshoot but worsened HPWL. This rejects the
  hypothesis that the smoothest overflow trajectory is automatically best.
- Full-budget tuning was necessary because 1200- and 3000-step runs begin
  cooldown at different fine iterations. The a1 full-scale optimum was 10,
  while a2 Adam scales 5 and 7 were statistically tied.
- Hard late scale switches were negative on both instances and are not
  promoted. Retaining optimizer moments does not make an abrupt scale change
  trajectory-neutral.
- Final best results are 72.576M/6.9488% for adaptec1 and 90.425M/6.9795% for
  adaptec2, gaps of 3.24% and 14.03% to the paper reference HPWL values.
- Archived all commands, stdout, convergence, lambda and gradient diagnostics.
  All formal files contain 3000 rows, no NaN/Inf tokens were found, and all six
  relevant unit-test targets passed.
- Fixed the experiment runner to copy each generated placement into its own
  run directory. Repeated the two final configurations: a1 reached
  72.774M/6.9492% and a2 reached 90.465M/6.9805% with archived placements.
