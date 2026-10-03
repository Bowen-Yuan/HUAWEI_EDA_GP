# Findings

## Current Understanding

- Under the strict 800-step budget, feasible-boundary trust changes the late
  exact-subgradient dynamics materially. Three adaptec1 runs reach 68.00--68.45M
  at overflow below 7%, versus the 70.3M reference.
- adaptec2 reaches 73.81M/6.9999% from the existing ePlace-GP placement with a
  local Adam scale of 0.25. This is a warm-start result. The best center-start
  result is 85.48M/6.9653% and does not beat the 79.3M reference.
- The trust mechanism checkpoints a feasible state, rejects boundary crossings
  or HPWL-increasing candidates, rolls back, and contracts the local step. It
  preserves exact max-minus-min HPWL and acts only as an acceptance rule.

- The trajectory controller already reaches and tracks the 7% feasibility boundary.
- Official DREAMPlace ISPD2005 configurations use target density 1.0, while the
  current solver used 0.8. This changes the feasible region and is the strongest
  identified explanation for excessive spreading and HPWL.
- adaptec2 continues to reduce exact HPWL by about 2.42M over its last 400 fine
  steps while the density/HPWL force ratio is about 0.059. The 3000-step state is
  therefore budget-limited or constrained-direction-limited, not simply over-penalized.
- Hard late step-scale switches previously regressed both adaptec1 and adaptec2.
- Predictive lambda reduced overflow overshoot but worsened HPWL, so temporary
  over-spreading should not be eliminated blindly.

## Lessons and Constraints

- Never smooth HPWL.
- Select the lowest feasible state, not the final state and not the lowest-overflow state.
- Do not reuse short-budget optimizer rankings without matching the full schedule.
- Do not use discontinuous optimizer-scale changes without state treatment.
- Do not force 150 HPWL-only steps on an already feasible mature GP warm-start;
  the phase deliberately clusters cells and destroys the useful density state.
- DREAMPlace's degree-100 gradient filter improves adaptec2 center runs, but only
  from 89.60M to 87.70M under the current 800-step flow.
- Strict and delayed exact-merit corridors stall before feasibility; their local
  electric directions cannot reduce the proposed merit at the activation point.
- Bundle cuts are harmful before 7% and remain negative when activated only
  after feasibility in this short-budget setting.

## Open Questions

- Can a non-destructive initialization reproduce the adaptec2 GP warm-start
  quality inside the same 800-step accounting?
- Can a serious/null-step bundle or a true feasible line search escape the
  85.48M center-start active-set basin without smoothing HPWL?
- How much of the remaining comparison gap comes from non-identical overflow,
  legalization, and detailed-placement operators?
