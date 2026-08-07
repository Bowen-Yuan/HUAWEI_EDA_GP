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
