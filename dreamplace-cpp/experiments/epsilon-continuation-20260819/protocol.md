# Three-stage exact-HPWL epsilon continuation

The final trial keeps the historical fixed-epsilon command and appends two
opt-in phases:

1. The original 1400/1600 iterations use the unchanged fixed radius, Adam,
   lambda controller and learning-rate schedule.
2. 200 iterations reset optimizer moments and linearly blend the fixed radius
   to a per-net/per-axis radius capped at 15% of the current net span. These
   steps use Adam with a 0.75 learning-rate scale.
3. 100 iterations set the active radius to zero, remove the density gradient,
   and use the exact HPWL subgradient with AMSGrad at scale 0.1. A trial is
   accepted only after exact HPWL and overflow filters pass.

Feasible checkpoints are ranked by the existing quick legalizer every 100
iterations. The final full legalizer is run only on the selected checkpoint.
The original command remains unchanged when the appended counts are zero.

The 20-iteration regression test compared legacy and extended CSV prefixes and
reported zero HPWL, overflow, or learning-rate differences. Both final a1/a2
placements passed internal legality and ISPD 2005 `legal2.pl` with Type
0/1/2/3 equal to 0/0/0/0.
