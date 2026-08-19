# Findings

## Current understanding

The previous early fixed-iteration switch improved adaptec1 slightly but
worsened adaptec2. Continuous degeneration to a zero radius lowered feasible
GP HPWL without lowering final legal HPWL. The current experiment therefore
tests a state-based transition, keeps a nonzero 10% or 15% per-net cap, and
treats strict-subgradient optimization as a rollback-safe candidate phase.

## Open questions

- Does state triggering choose meaningfully different transition iterations?
- Is 10% or 15% more robust after full legalization?
- Can a short strict-subgradient phase improve legal-aware checkpoint quality?

## Final findings

- State-triggered 15% continuation improves four of six designs relative to
  the frozen adaptive epsilon baseline.
- The geometric-mean legal HPWL improvement is about 0.151%.
- A short strict-subgradient tail was worse on adaptec1 and was not retained.
- All retained placements pass the official checker with zero errors.
- GP-to-legal conversion remains the limiting mechanism on adaptec1 and
  bigblue1, where GP improvements did not reduce final legal HPWL.
