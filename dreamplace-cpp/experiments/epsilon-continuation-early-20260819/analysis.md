# Early epsilon-continuation experiments

All comparisons keep the total iteration budget and the original learning-rate
clock unchanged. Entering stage 2 earlier improves adaptec1 final legal HPWL by
0.010M relative to the late-switch continuation, but worsens adaptec2 by
0.036M. The effect therefore does not generalize across the two designs.

On adaptec1, continuously reducing epsilon to zero reaches the lowest feasible
GP HPWL (71.888M), but legal-aware selection returns to iteration 1200 and the
final legal HPWL is 76.991M. Stopping at the 15% per-net target is slightly
better after legalization (76.988M). Direct degeneration to zero helps the
continuous exact objective but not the final legal objective.

The retained early-to-15% a1/a2 placements pass ISPD 2005 `legal2.pl` with
Type 0/1/2/3 equal to 0/0/0/0.
