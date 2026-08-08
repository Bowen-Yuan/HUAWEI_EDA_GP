# Progressive Legalization Experiment Protocol

## Research Question

Can an exact, nonsmooth-HPWL global placer reduce the legalization shock by
separating wirelength contraction, electrostatic spreading, and continuous
macro avoidance before the existing discrete legalizer runs?

## Locked Constraints

- The reported and accepted wirelength objective is exact HPWL with pin
  offsets. Weighted-average or log-sum-exp wirelength is not used by this
  experimental path.
- The existing command line remains valid. The new path is disabled unless
  `--progressive-legalization` is supplied.
- Phase 3 may add continuous gradients and reject unsafe trials, but it may
  not copy, blend, snap, or project coordinates from Greedy/Abacus.
- Final quality is reported only after the full internal legality check.
- Negative and partially successful results are retained.

## Hypotheses

### H1: Short Exact-HPWL Contraction

Seventy-five to one hundred exact-HPWL-only AMSGrad steps should reduce the
initial center placement enough to give the density controller wirelength
headroom, without spending so long in the highly overlapped state that Adam's
moments become unusable for spreading.

Prediction: 100 steps will have lower HPWL than 75 steps, but may require a
moment reset at density activation. The preferred setting is selected by the
final phase-2 HPWL/overflow pair, not by phase-1 HPWL alone.

### H2: Band-Dual Density Control

A log-domain proportional-integral density multiplier with an explicit
`7%-8%` overflow band and a 75M HPWL budget should avoid the runaway density
weight observed with monotone DREAMPlace-style multiplication.

Prediction: compared with the existing trajectory controller, the band-dual
controller will enter the target band without continuing to increase lambda
after overflow falls below its midpoint.

### H3: Continuous Macro Avoidance

A gradually ramped squared distance-to-obstacle-complement penalty, combined
with a weak row-attraction term, should reduce movable/fixed overlap before
Greedy legalization without imposing row, site, segment, or order constraints.

Prediction: phase 3 will reduce both the movable area overlapping fixed macros
and the number of overlapping movable cells. The Greedy HPWL increase should
also be smaller than for a matched run that omits the continuous obstacle
penalty.

## Five Stages

1. Exact nonsmooth HPWL only for 75 or 100 iterations.
2. Exact HPWL plus electrostatic density, controlled toward 7.5% overflow.
3. Continue stage 2 while smoothly ramping macro and row penalties.
4. Existing obstacle-segment Greedy plus Abacus legalization.
5. Existing DREAMPlace-style K-Reorder, insertion, global swap,
   independent-set matching, and fixed-topology legal Bundle refinement.

## Primary Metrics

- exact HPWL at each stage boundary;
- grid overflow and electrostatic density weight;
- movable/fixed overlap area divided by movable area;
- number of movable cells overlapping fixed objects;
- mean row distance and segment overflow proxy;
- HPWL after Greedy, Abacus, and detailed placement;
- internal boundary, site-alignment, fixed-overlap, and movable-overlap errors.

## Evaluation Matrix

The initial search uses adaptec1 and a fixed seed of 1000. Short ablations use
a reduced iteration budget and retain the same density grid where practical.
The final selected configuration is compared against:

- the existing center-initialized exact path;
- the same progressive schedule without phase-3 macro/row penalties;
- the validated smooth DREAMPlace-style flow (context only, not an exact-HPWL
  ablation);
- the previous best legal fixed-topology Bundle result, 74.425850M.

The target `75M / 7%-8%` is treated as an engineering goal, not assumed to be
achieved. Any miss is reported with its stage-wise failure mechanism.
