# Results And Analysis

## Pre-Legal Grouped Bundle

The 20-step exact stage used four net groups, three cuts per group, two nearby
exact-subgradient samples every five iterations, and exact serious/null
filtering.  Boundary bisection kept the selected state strictly feasible.

- Warm start: 71.235743M / 6.9654%
- Selected exact state: 71.219645M / 7.0000%
- GP improvement: 16.098K

After reaching the boundary, most proposed directions pointed outside the
feasible set and bisection fractions collapsed to zero.  More iterations with
the same direction model do not help; a true density-tangent projection would
be required.

With three detailed outer rounds and eight legal-bundle passes, the grouped
pre-legal run finalized at 74.484991M.  The matched run without pre-legal
refinement finalized at 74.475888M.  Thus the lower GP HPWL worsened final
legal HPWL by 9.103K.  The pre-legal hypothesis is not supported.

## Post-Legal Fixed-Topology Bundle

All trials keep the current row, obstacle-free segment, and left-to-right
order.  The bundle master produces an x direction, segment Abacus projects it
onto non-overlap and site constraints, and a complete legality plus exact
HPWL check controls serious/null acceptance.

Three-round matched comparison:

| Configuration | Final legal HPWL |
|---|---:|
| No legal bundle | 74.548626M |
| 8 passes, 8 cuts | 74.475888M |
| 16 passes, 12 cuts | 74.445861M |

Five-round matched comparison:

| Configuration | Final legal HPWL |
|---|---:|
| No legal bundle | 74.520409M |
| 16 passes, 12 cuts | 74.425850M |

The five-round net improvement attributable to the legal bundle is 94.559K,
or 0.127%.  Later rounds contain genuine null steps and proximal-radius
contraction, demonstrating that the method is using bundle history rather
than merely repeating the old projected subgradient operator.

## Conclusion

The fixed-topology post-legal bundle is supported and remains opt-in.  The
grouped pre-legal bundle improves the continuous GP metric but is rejected as
the best final pipeline because that improvement does not survive
legalization.  The next mathematically relevant pre-legal extension is a
density-tangent constrained master, not additional iterations or a wider
overflow tolerance.
