# Analysis

The exact bundle run selected iteration 10 by quick legal HPWL:

- GP HPWL: 71.220652M
- Overflow: 6.9976%
- Final legal HPWL: 74.549788M

The matched direct-legalization control produced:

- GP HPWL: 71.235743M
- Overflow: 6.9654%
- Final legal HPWL: 74.548626M

The hypothesis is not supported by this configuration.  Exact GP HPWL improved
by 15.090K, but final legal HPWL worsened by 1.162K.  Six serious steps were
accepted early.  Later trials repeatedly lowered HPWL while landing just above
the hard 7% overflow boundary, causing the proximal scale to collapse and the
center to stop moving.  Future work should test a filter/hysteresis constraint
controller and rank serious trial points by matched full legalization rather
than by GP HPWL alone.
