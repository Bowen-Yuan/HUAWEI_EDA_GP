# Findings

## Current Understanding

- DREAMPlace does not use ePlace bound-to-bound initialization in the reported
  flow.  Movable cells start at the layout center with Gaussian standard
  deviation equal to 0.1% of each layout dimension.
- Fixed objects enter the initial electrostatic density map.  Small movable
  cells are stretched to at least sqrt(2) bin widths/heights while preserving
  area through a density ratio.
- Legalization is a separate discrete phase: macro legalization, Tetris-like
  Greedy legalization, then Abacus.  Detailed placement is required for the
  final Table II HPWL.

## Lessons And Constraints

- Comparing a pre-legal HPWL to paper Table II is invalid.
- A legalizer that ignores pin offsets can select the wrong local moves.
- A capacity-only density metric is insufficient evidence of geometric
  legality; the official `legal2.pl` result is authoritative.
- Overflow must be evaluated without fillers.  Fillers participate in the
  electrostatic objective, but DREAMPlace passes zero filler nodes to its
  overflow operator.
- BB curvature must compare gradients evaluated with the same density weight
  and wirelength gamma.  Mixing gradients from consecutive outer objectives
  stalls spreading.
- Greedy legalization must preserve the GP left-to-right order.  A global
  width-first ordering changed adaptec1 post-legal HPWL from 74.819M to
  123.992M at the same GP checkpoint.

## Validated Outcome

- adaptec1: 71.235743M at 6.9654% overflow before legalization; 74.819421M
  after legalization/refinement; `legal2.pl` Type 0/1/2/3 all zero.
- adaptec2: 80.889298M at 6.9302% overflow before legalization; 84.426599M
  after legalization/refinement; `legal2.pl` Type 0/1/2/3 all zero.
- Relative to paper Table II, the final legal gaps are +2.18% and +2.68%.
