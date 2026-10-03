# Autonomous DREAMPlace PK Log

## 2026-08-04 Bootstrap

- Locked metric: lowest exact HPWL state with overflow no greater than 0.0700.
- External targets: adaptec1 70.3M/6.92%; adaptec2 79.3M/6.89%.
- Exact HPWL smoothing remains forbidden.
- Git preregistration is unavailable because the repository `.git` directory is empty.
- Ranked mechanisms: unfinished recovery, inefficient boundary-normal direction,
  discontinuous late scale, then density discretization/macro normalization.

## 2026-08-04 Budget Revision

- User fixed the budget at 800 total iterations. The running 6000-step H1 jobs
  were terminated and are excluded from all comparisons.
- The official DREAMPlace `test/ispd2005/adaptec1.json` and `adaptec2.json`
  both set `target_density=1.0` and `stop_overflow=0.07`; the official bin counts
  are 512x512 and 1024x1024 respectively.
- The current solver had used its generic default target density 0.8. H0 tests
  official density-capacity alignment under the strict 800-step budget.

## 2026-08-04 H0/H2 Reflection

- H0: a1 75.853M/6.947%; a2 97.690M/6.986%. Official density alignment alone
  is insufficient. First feasible fine steps were 470 and 417.
- H2 control (a2, horizon 200, scale 9): 99.605M/6.976%.
- H2 wide projection: 110.096M/1.804%; it removes too much HPWL resistance and
  grossly over-spreads. Wide-band projection is rejected.
- A dormant project artifact, `adaptec2.eplace.gp.pl`, evaluates at
  85.029M/4.4788% with zero new iterations. H3 explores it as an explicitly
  labeled warm-start, not a from-center result.

## 2026-08-05 Strict 800-Step Completion

- Added opt-in feasible-boundary trust. It keeps a feasible checkpoint, rejects
  overflow crossings and HPWL-increasing candidates, rolls back, and contracts
  the local step. Existing invocations are unchanged.
- adaptec1 crossed the 70.3M reference in three independent 800-step runs. The
  final visualized repeat reached 68.002M/6.8662% from center initialization.
- A low-scale adaptec2 ePlace-GP warm-start reached 73.809M/6.9999% in a final
  visualized 800-step run. It is explicitly not a center-start comparison.
- The best adaptec2 center-start result was 85.481M/6.9653% using 100 medium
  preconditioning steps, degree-100 gradient filtering, and feasible trust.
- Shorter trajectory horizons, strict/delayed merit corridors, early bundle,
  and post-feasible bundle did not cross the adaptec2 center reference.
- Generated final GIF animations plus vector/raster convergence and baseline
  comparison figures under `to_human/final_800`.
