# Epsilon-active performance run, 2026-08-18

This run reproduces the previously frozen smart adaptive epsilon-active
configuration on adaptec1--4 and bigblue1--2.  No parameter was changed
during the runs.

Common settings:

- exact nonsmooth HPWL, electrostatic density, Adam before and after feasible refinement;
- center Gaussian initialization, seed 1000;
- epsilon active-set power 4;
- smart controller: window 50, update interval 25, gain 0.50,
  deadband 0.0025, scale range 0.70--1.10, maximum log step 0.03;
- feasible refinement, legal checkpoint interval 50;
- DREAMPlace-style detailed placement: 5 outer passes, 5 legal refinement
  rounds, 2 insertion passes, 4 projected passes, 2 row re-legalization
  passes, independent-set size 8, Hungarian matching;
- no bundle, no gradient sampling, no snapshots during the timed run.

| design | bins | iterations | global epsilon | refinement epsilon |
|---|---:|---:|---:|---:|
| adaptec1 | 512x512 | 1400 | 125 | 50 |
| adaptec2 | 1024x1024 | 1600 | 165 | 65 |
| adaptec3 | 1024x1024 | 1600 | 265 | 105 |
| adaptec4 | 1024x1024 | 1600 | 265 | 105 |
| bigblue1 | 512x512 | 1600 | 125 | 50 |
| bigblue2 | 1024x1024 | 1600 | 212 | 85 |

`results_summary.csv` contains runtime, GP HPWL/overflow, final detailed HPWL,
and internal legality counts.  The static convergence figures and GIF are
generated after all processes finish from `global_metrics.csv`; their render
time is not included in runtime.
