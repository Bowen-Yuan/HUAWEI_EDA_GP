# Findings

## Current understanding

- Exact HPWL is convex and piecewise linear, so every evaluated subgradient
  defines a globally valid affine lower cut.
- Heavy-ball already filters the final combined gradient, but it does not use
  HPWL function values or linearization errors. A bundle can reject stale cuts
  whose affine model is poor at the current point.
- The density objective and lambda vary over time, so a full bundle model of the
  combined objective would not be stationary. The first implementation bundles
  only exact HPWL, then combines the selected HPWL direction with current
  electrostatic density gradient.
- H1 is supported: a four-cut proximal aggregate improved 3000-step feasible
  HPWL from 79.423M to 76.522M without changing the density or lambda models.
- The useful regime is not simple temporal averaging. With prox scale 1 the
  newest cut retained about 95% weight; prox scale 2 produced a high-entropy
  mixture, about one-third newest-cut weight, and a 0.40 norm ratio.
- At 1800 steps the bundle entered feasibility slightly later than baseline
  but with much lower HPWL and less density overshoot. Its benefit is therefore
  better movement along the HPWL-density tradeoff, not merely faster spreading.
- Stronger explicit history (current mix 0.35) looked better at 1800 but worse
  at 3000. Mid-run feasible HPWL and final recovery quality are not monotonic.
- Prox scale 3 produced nearly the same aggregate norm ratio as B2 but worse
  overflow slack. Aggregate norm alone does not explain performance; cut age
  and trajectory-dependent lambda coupling matter.
- Three cuts outperform four at 3000 steps (76.164M vs 76.522M), while six cuts
  are too conservative. A short memory of recent extrema active sets appears
  more useful than a richer approximation of older HPWL geometry.
- Exact HPWL subgradients exhibit a pronounced sign-changing cycle: at a
  two-step separation, global cosine averages about -0.45. This quantitatively
  confirms that extreme-pin switching is a central optimizer issue. A naive
  cosine trigger fires constantly and therefore provides no selectivity.
- A two-cut model performs as well as the three-cut best (76.120M vs 76.164M),
  supporting the alternating-active-set interpretation. However, this small
  difference is not statistically resolved.
- A repeat of the two-cut model reached 75.677M, so the bundle improvement over
  79.423M is robust, while sub-million differences between bundle variants are
  not.
- Reducing lambda-up to 1.13 worsened final bundle HPWL. The original continuous
  lambda controller should remain fixed while testing direction mechanisms.
- Direct cosine-adaptive mixing is too reactive because the negative-extrema
  cycle persists throughout optimization. If used, it must be annealed away so
  late iterations recover wirelength with the fixed bundle/current balance.

## Lessons and constraints

- Never smooth HPWL.
- Compare strict feasible HPWL, not raw final iterate HPWL.
- Screen at 800 steps, then promote only plausible configurations to 3000.
- The unchanged staged solver is infeasible at 800 steps; use 1800-step feasible
  screening after the initial numerical sanity check.

## Open questions

- What bundle size and cut insertion interval preserve useful extrema history?
- How much current-cut mixing is required near the 6% boundary?
- Does a bundle improve HPWL recovery or merely reduce gradient norm?

## Final synthesis for this cycle

- Two independent size-2 bundle runs achieved 76.120M and 75.677M at strict
  feasible overflow, both substantially below the 79.423M historical baseline.
- Best observed improvement is 4.72%; the conservative repeat-range conclusion
  is an improvement of 3.30M to 3.75M.
- Mechanism evidence: exact HPWL gradients separated by two steps have strongly
  negative cosine, and small bundles cancel this extrema-induced oscillation
  using valid affine cuts and their linearization errors.
- Three cuts are statistically similar; four are slightly worse; six over-
  regularize. Interval five is better than intervals three and seven.
- Keep the method opt-in pending adaptec2-4 validation and deterministic
  multi-run statistics.

## 2026-08-01 relaxed-overflow cycle

### New supported result

- Separating the controller normalization baseline from the best-state
  feasibility limit is necessary. Raising one shared baseline to 0.065 changed
  the early spread trajectory and increased late effective lambda; it did not
  convert extra overflow into lower HPWL.
- Keep the controller baseline at 0.060, set the fine feasibility limit to
  0.065, and target 0.0648 with the guarded controller. This preserves the
  validated early trajectory and uses extra capacity only during fine recovery.
- Raising the final learning-rate multiplier from 0.15 to 0.25 was the useful
  optimizer change. Two independent runs reached 75.599M/0.064875 and
  75.638M/0.064801, a tight 0.039M range.
- The best observed gain over the prior 75.677M bundle result is 0.078M
  (0.104%). This is real but modest; the stronger conclusion is the replicated
  75.599--75.638M range, not a large benefit from relaxing overflow.

### Mechanistic interpretation

- At cooldown floor 0.15, overflow and lambda were already stable while HPWL
  still fell by 1.533M over the last 400 fine steps. The 3000-step budget ended
  during ongoing wirelength recovery.
- Floor 0.25 accelerates that late recovery without destabilizing the 0.065
  boundary. Best checkpoints remain near fine step 2596.
- Best late HPWL checkpoints recur every five steps, synchronized with bundle
  cut insertion. Floor 0.35 amplifies this nonsmooth cycle and regresses.
- Lambda gains are defined per controller update. Reducing the interval from 25
  to 10 without rescaling gains drove lambda from about 5.2 to 7.45 and worsened
  HPWL to 77.675M. Feedback frequency and gain cannot be tuned independently.
- Moving the target from 0.0648 to 0.0649 did not help, and increasing current
  exact-cut mix from 0.50 to 0.55 also regressed. The extra 0.5 percentage-point
  feasibility budget is not the main remaining bottleneck.

### Updated constraints

- Recommended relaxed configuration: two cuts, interval 5, prox 2.0, current
  mix 0.50, overflow limit 0.065, recovery target 0.0648, guard interval 25,
  cooldown floor 0.25.
- Do not use floor >= 0.35 or unnormalized high-frequency lambda updates.
- HPWL remains exact max-minus-min; only the density model is smooth.

## 2026-08-02 eight-benchmark generalization cycle

- All eight ISPD2005 runs completed 3000 exact-HPWL bundle steps and produced a
  selected state below the DREAMPlace-guided 10% reference limit.
- Final selected overflow is robust across seven instances (6.410--6.520%);
  bigblue3 is higher at 7.162% and reaches 10% only at global iteration 2162.
- Density feasibility transfers better than HPWL quality. Relative to the
  different-flow DREAMPlace Table II numbers, pre-legalization HPWL ratios span
  1.03--2.33 and worsen on bigblue3/4.
- Initialization is a major unresolved confounder: only adaptec1 has an
  ePlace-IP file. The other seven fall back to supplied Bookshelf `.pl` files.
- The fixed controller is not scale invariant. Its 7e7 HPWL baseline and
  effective density-weight scaling were tuned on adaptec1; bigblue3 visibly
  saturates lambda before becoming feasible.
- A valid next comparison requires consistent initialization plus per-instance
  gradient-scale normalization, followed by paired bundle/no-bundle runs.

## 2026-08-04 adaptec2--4 optimization cycle

- The missing ePlace-IP files were the dominant confounder. A deterministic
  DREAMPlace-style center initialization reduced post-Phase-0 HPWL from
  127/347/353M to 68.8/170.0/165.4M on adaptec2/3/4.
- Best 3000-step exact-HPWL results at this solver's overflow <=10% are
  88.367M, 241.559M, and 203.894M: improvements of 23.74%, 25.98%, and 35.73%
  over the previous fixed-configuration runs.
- These values do not reach DREAMPlace Table II. Remaining gaps are 7.48%,
  24.69%, and 17.13%, and the comparison remains different-flow/contextual.
- Aligning the lambda target to 10% gives only 4--6% improvement by itself;
  initialization gives 19--32% after target alignment.
- Bundle is still useful on all three. Removing it worsens final HPWL by
  4--9%; however, both stronger history (mix 0.35) and more current direction
  (mix 0.65) regress on adaptec2.
- A 600-step trajectory PI-D controller reaches feasibility earlier and holds
  9.95% accurately. Removing the coarse/medium phases reallocates 250 steps to
  constrained fine recovery and provides a small additional gain.
- In best runs, bundle/current HPWL RMS is about 0.49 and late density/HPWL
  force RMS is only 0.067--0.127. The remaining gap is therefore not primarily
  caused by excessive late lambda; it is a constrained descent-efficiency and
  flow-equivalence problem.
- Increasing cooldown floor to 0.50, using 3 cuts, lengthening Phase 0, and
  current mixes 0.35/0.65 are negative results. Prox 3 is only marginally
  positive.

## 2026-08-04 conventional-gradient and lambda cycle

- A no-bundle Adam method with exact HPWL improves the previous bundle result
  on adaptec2/3/4. Best 3000-step results are 87.162M, 230.189M and 201.272M at
  9.91--9.95% overflow, improving bundle by 1.36%, 4.71% and 1.29%.
- A controlled 3000-step adaptec2 run with identical initialization, phase
  allocation and trajectory lambda gives 93.075M for heavy-ball versus
  87.162M for tuned Adam. The new gain is therefore not a phase/controller
  confounder.
- Adam does not smooth HPWL. It receives the current exact max-minus-min
  subgradient and applies coordinate-wise first/second-moment scaling after
  evaluation. Beta1=0.90 filters extrema chatter better than 0.80; 0.95 is too
  stale.
- The original Adam step multiplier 0.70 was much too conservative after the
  solver's global RMS normalization. Scales 3--5 form the useful region;
  adaptec2/3 prefer 5 while adaptec4 prefers 3, so the method is not yet scale
  invariant.
- Trajectory PI-D lambda beats both guarded proportional feedback and a new
  projected dual rule. On the 1200-step screen it gives 89.984M, versus
  92.860M guarded and 93.173M for the best feasible dual rate.
- Projected dual rate 1 is infeasible. Rates 3--8 become feasible, but faster
  feasibility does not yield better HPWL recovery; one gain cannot serve both
  the density ramp and boundary regulation.
- Adaptive restart is strongly harmful (146.263M at 1200 steps) because
  negative exact-HPWL gradient dot products are intrinsic active-pin switches,
  not exceptional optimizer failure.
- All promoted Adam runs select their last fine-grid state. The final 400 steps
  still reduce HPWL by 2.36M/6.82M/4.50M, while density/HPWL force ratios are
  only 0.052/0.045/0.029. Remaining gaps are budget/constrained-descent issues.
- DREAMPlace Table II remains unreached: tuned-Adam gaps are 6.01%, 18.83% and
  15.62%, under a still different placement and overflow flow.

## 2026-08-04 seven-percent adaptec1/adaptec2 cycle

- The relevant external PK row is Table 5 of the stochastic-subgradient paper:
  DREAMPlace reports pre-legalization global placement at 70.3M/6.92% for
  adaptec1 and 79.3M/6.89% for adaptec2.
- With a strict 7% feasibility limit, the best 3000-step exact-HPWL results are
  72.576M/6.9488% for adaptec1 and 90.425M/6.9795% for adaptec2. The gaps to
  the paper references are 3.24% and 14.03%.
- Final instance-specific settings are AMSGrad scale 10 with target 0.0695 for
  adaptec1, and Adam scale 5 with target 0.0698 for adaptec2. Both retain the
  150-step lambda-zero HPWL-only phase and a 600-fine-step density trajectory.
- Short-budget optimizer rankings do not transfer directly because the
  cooldown starts at fine step 250 for a 1200-step run but at 1800 for a
  3000-step run. Formal tuning must compare full-budget runs.
- Predictive overflow look-ahead reduced a1 lambda peak and density overshoot,
  but worsened feasible HPWL. Temporary over-spreading is useful for later
  exact-HPWL recovery in this solver.
- Hard late scale switches also regress on both instances. The optimizer
  moments and the lambda/overflow orbit are coupled; a discontinuous scale
  increase disturbs the constrained trajectory.
- The best a2 run still drops about 2.42M HPWL in its final 400 steps while the
  late density/HPWL force ratio is only about 0.059. The a2 gap is therefore not
  explained by excessive late lambda alone.
- Final repeats reached 72.774M/6.9492% and 90.465M/6.9805%. Their differences
  from the lowest observations are 0.27% and 0.044%; a1 has more visible
  trajectory variation, but neither repeat changes the optimizer conclusion.
