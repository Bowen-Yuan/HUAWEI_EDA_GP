# Per-net relative epsilon experiment (2026-08-19)

## Question

Can an exact-HPWL epsilon-active direction retain the existing a1/a2 quality
when its support is made genuinely local?  For each net and coordinate axis,

```text
epsilon(e,d) = min(global_epsilon,
                   max(alpha * span(e,d), epsilon_floor_if_small))
```

The floor is allowed only when `span(e,d) < small_span_threshold`.

## Frozen controls

- Raw ISPD 2005 Bookshelf input, center-Gaussian initialization, seed 1000.
- Exact HPWL evaluation and Adam trial optimizer.
- Original DREAMPlace lambda policy and 6.5%-7.0% feasible refinement.
- Same grid, iteration budget, legal-aware checkpoint selection, Greedy,
  Abacus and DREAMPlace-style detailed placement as the 2026-08-18 run.
- The old global adaptive epsilon controller is disabled.
- No bundle and no smooth HPWL objective.

## First ablation

| Case | alpha | power | small-span floor | small-span threshold |
|---|---:|---:|---:|---:|
| rel10-p4-f1 | 0.10 | 4 | 1 DBU | 10 DBU |
| rel15-p4-f1 | 0.15 | 4 | 1 DBU | 10 DBU |
| rel15-p2-f1 | 0.15 | 2 | 1 DBU | 10 DBU |
| rel15-p1-f1 | 0.15 | 1 | 1 DBU | 10 DBU |
| rel15-p2-lambda-match | 0.15 | 2 | 1 DBU | 10 DBU |
| rel15-p2-f5-small50 | 0.15 | 2 | 5 DBU | 50 DBU |
| rel15-p2-f10-small100 | 0.15 | 2 | 10 DBU | 100 DBU |
| rel15-p2-f10-to-f1 | 0.15 | 2 | 10 -> 1 DBU | 100 -> 10 DBU |
| rel15-p2-f15-to-strict | 0.15 | 2 | 15 -> 0 DBU | 100 -> 0 DBU |
| rel15-p2-f20-to-strict | 0.15 | 2 | 20 -> 0 DBU | 133.33 -> 0 DBU |

The normal-net epsilon cannot exceed 10% or 15% of its current directional
span. Only spans below 10 DBU may use the 1 DBU floor.

The `p=2` and `p=1` cases are promoted only after `p=4` showed that changing
10% to 15% alone did not restore quality. They keep the same 15% support and
only redistribute weight inside that local boundary band.

The promoted `lambda-match` case changes the initial density scale from
`8.0e-5` to `2.25e-5`. This compensates for the local direction's roughly
3.5x larger initial wire-gradient L1 norm and matches the old controller's
initial effective lambda; all subsequent DREAMPlace lambda updates remain
unchanged.

The `f5-small50` case addresses the center-Gaussian initialization, where
many net spans are only a few DBU. On the measured a1 initial state the floor
is used on about 79% of directions, but complete-span coverage is limited to
about 14% rather than the old method's 79%--90%.

The promoted `f10-small100` case raises initial complete-span coverage to
about 27% on a1. The 100 DBU exception threshold is below 1% of the roughly
10.7K DBU placement width; nets outside the exception remain strictly capped
at 15% per coordinate axis.

The `f10-to-f1` case confines the larger exception to density spreading.
When feasible refinement starts, Adam state is already reset by the frozen
pipeline and the floor changes to 1 DBU for spans below 10 DBU. Normal nets
remain capped at 15% throughout both stages.

The final `f15-to-strict` candidate uses the continuous early rule
`max(15 DBU, 0.15 * span)` below 100 DBU. Initial complete-span coverage on
a1 is about 40%, still below half. The exception is removed entirely at
refinement, leaving the strict per-net 15% cap.

The bounded final screen raises the continuous early floor to 20 DBU below
133.33 DBU spans. It covers about 53%--55% of a1 directions completely at
the measured initial state, but the exception is again removed entirely at
refinement. Larger floors are excluded from this study because they approach
the rejected majority/full-support regime.

Primary metrics are exact GP HPWL at overflow no greater than 7%, final legal
exact HPWL, internal/external legality and runtime.  Existing smart-adaptive
results are reference data and are not rerun.
