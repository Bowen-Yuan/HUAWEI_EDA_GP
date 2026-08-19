# State-triggered epsilon continuation results

The frozen production variant is state-triggered continuation to a 15%
per-net, per-axis span cap. It keeps the original iteration budgets (1400 for
adaptec1 and 1600 for the other designs), uses Adam in the fixed and
continuation phases, and does not force a pure subgradient tail. The state
trigger requires at least 250 refinement iterations, a 50-step window with
maximum overflow <= 8.5%, overflow range <= 0.4 percentage points, and a
relative exact-HPWL drop >= 0.03%; the configured stage boundary remains the
fallback.

## Final results

| design | iterations | transition | GP HPWL (M) | GP overflow | legal HPWL (M) | runtime (s) | GP (s) | legal2 Type 0/1/2/3 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| adaptec1 | 1400 | 898 | 72.445 | 6.991% | 76.666 | 239.1 | 131.7 | 0/0/0/0 |
| adaptec2 | 1600 | 1025 | 79.472 | 6.992% | 82.240 | 170.9 | 126.1 | 0/0/0/0 |
| adaptec3 | 1600 | 1199 | 188.273 | 6.998% | 199.512 | 575.1 | 422.6 | 0/0/0/0 |
| adaptec4 | 1600 | 1077 | 174.431 | 6.992% | 182.498 | 699.5 | 513.5 | 0/0/0/0 |
| bigblue1 | 1600 | 941 | 87.109 | 6.998% | 89.815 | 146.4 | 100.8 | 0/0/0/0 |
| bigblue2 | 1600 | 1070 | 135.528 | 6.999% | 142.870 | 850.0 | 390.6 | 0/0/0/0 |

Runtime includes global placement and the full detailed legalizer, but excludes
the subsequent official Perl checker; `GP (s)` is the global-placement time
reported by the solver. The checker is
the ISPD 2005 `legal2.pl` invoked with the TeX Live Perl interpreter at
`D:\\latex\\101\\texlive\\2023\\tlpkg\\tlperl\\bin\\perl.exe`.

## Comparison with the frozen adaptive epsilon baseline

The previous legal HPWL values were 76.635M, 82.372M, 199.907M, 183.252M,
89.743M and 143.239M for adaptec1--4 and bigblue1--2. The new relative changes
are respectively +0.040%, -0.160%, -0.198%, -0.412%, +0.081% and -0.257%.
The six-design geometric-mean legal HPWL improves by approximately 0.151%.

The improvement is therefore real but modest: the state trigger avoids a
fixed early switch that hurt adaptec2, while the nonzero 15% cap preserves the
stable epsilon-active direction. The remaining limitation is legal-quality
conversion: lower GP HPWL does not always imply lower legal HPWL, as seen on
adaptec1 and bigblue1.

## Curves and artifacts

`plot_convergence.py` creates `fig_state_triggered_hpwl.{png,pdf}` and
`fig_state_triggered_overflow.{png,pdf}` from the recorded global metrics.
The dashed vertical markers show the actual state-triggered transition.
