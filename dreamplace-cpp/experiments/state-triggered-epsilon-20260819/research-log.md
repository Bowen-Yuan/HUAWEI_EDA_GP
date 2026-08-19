# Research log

## 2026-08-19: protocol and implementation

- Locked the protocol before running confirmatory experiments.
- Added an opt-in state trigger; legacy calls retain the original constant
  stage boundary.
- Preserved the total objective iteration budget by stretching continuation
  from an early trigger to the original fixed endpoint.
- Build, core tests, and a short fallback-transition smoke test passed.
- No Perl interpreter is available, so official `legal2.pl` validation is
  deferred; internal boundary, site-alignment, and overlap checks remain hard
  constraints.

## 2026-08-19: correction and frozen six-design run

- Corrected the Perl environment: TeX Live provides a usable interpreter at
  `D:\\latex\\101\\texlive\\2023\\tlpkg\\tlperl\\bin\\perl.exe`.
- The user requested an immediate frozen run, so no further tuning was done.
- Froze the state-triggered 15% span-cap variant after adaptec1/adaptec2.
- Completed adaptec1--4 and bigblue1--2; all six passed official `legal2.pl`
  with Type 0/1/2/3 equal to 0/0/0/0.
- Improved the geometric-mean legal HPWL by approximately 0.151% relative to
  the frozen adaptive epsilon baseline.
