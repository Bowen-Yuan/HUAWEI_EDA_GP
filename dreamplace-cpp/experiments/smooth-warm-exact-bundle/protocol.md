# Smooth Warm Start And Exact Bundle Refinement

Status: exploratory

## Hypothesis

A validated smooth 512x512 global placement provides the missing global
topology, allowing a short exact-HPWL bundle phase to reduce exact wirelength
without losing the 7% density constraint or final legal quality.

## Locked Comparison

- Benchmark: adaptec1
- Warm start: `output/adaptec1_papergrid_600/global.pl`
- Exact refinement: 200 iterations, 512x512 density grid, 12 bundle cuts,
  serious/null acceptance, 0.01 feasible-refinement step scale
- Primary metric: final legal exact HPWL
- Secondary metrics: pre-legal exact HPWL and overflow
- Control: the same warm start and enhanced detailed placer with no effective
  exact global refinement

The experiment is successful only if final legal HPWL improves over the
matched control.  A pre-legal-only reduction is not sufficient.
