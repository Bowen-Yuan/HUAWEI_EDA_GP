# Exact-HPWL epsilon continuation results

The opt-in flow preserves the historical fixed-epsilon trajectory, then appends
200 Adam steps that linearly interpolate each net/axis radius to 15% of its
current span. A final 100 AMSGrad steps use the radius-zero exact HPWL
subgradient and accept a backtracked trial only when exact HPWL does not rise
and overflow remains at or below 7%.

Legal-aware checkpoints are evaluated every 100 iterations. Both final files
pass the internal checker and the official ISPD 2005 `legal2.pl` checker with
Type 0/1/2/3 counts equal to 0/0/0/0.

The continuation improves feasible GP HPWL substantially, but most of that
gain does not survive legalization. The selected adaptec1 checkpoint is fixed
stage iteration 1300; adaptec2 selects the restored fixed-stage solution at the
start of continuation. This is evidence that the remaining issue is
legalization topology preservation rather than failure of the exact HPWL
direction to reduce its stated objective.
