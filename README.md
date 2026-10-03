# Adaptive Pareto Soft-Penalty Placement

This project is the adaptive HPWL/density experiment built on the shared
`code/gpplacer` parser and exact density oracle.  Its primary metric is the
**strict overflow protocol**: conventional positive-capacity overflow plus
movable area that occupies zero-capacity (blocked) bins.  The legacy metric is
still recorded to make earlier results comparable, but it is not used for
selection or controller decisions.

## Main capabilities

- Guarded HPWL rescue, density rescue, dual emergency handling and adaptive
  joint soft-penalty descent.
- Hybrid density directions: use a short-axis direction only for strongly
  anisotropic local overlap; otherwise retain the two-dimensional direction.
- Exact-oracle-accepted, local cardinal density evacuation rather than global
  jumps to the largest-slack bin.
- Recursive hypergraph coarsening, a 32-seed diversified candidate portfolio,
  successive halving and a small dependency-free ridge proxy.
- Bounded strict Pareto archive, crash recovery checkpoints and run comparison.
- PNG/PDF figures plus an offline, self-contained Plotly report.

## Run

From `D:\codex_project\HUAWEI_EDA`:

```powershell
$env:PYTHONPATH = 'D:\codex_project\HUAWEI_EDA\code;D:\codex_project\HUAWEI_EDA\adaptive_pareto_soft\src'
python -m adaptive_pareto.cli solve `
  --aux dataset\adaptec1\adaptec1.aux `
  --placement code-hard\runs_alternating_25\adaptec1\20260720T144806Z\solution.pl `
  --config adaptive_pareto_soft\configs\adaptec1_smoke.json
```

For a portfolio/multilevel run, omit `--placement` and use a configuration
with `portfolio.enabled=true` (the default).  Resume a standalone seeded run:

```powershell
python -m adaptive_pareto.cli solve --resume-run adaptive_pareto_soft\results\adaptec1\<timestamp>
```

Render or compare previously completed runs:

```powershell
python -m adaptive_pareto.cli visualize-run --run-dir <run-dir>
python -m adaptive_pareto.cli compare-runs --run-dir <run-a> --run-dir <run-b> --output adaptive_pareto_soft\results\comparison
```

Every run contains placements, exact evaluation traces, snapshots,
checkpoints, Pareto archive, static figures and `report.html`.
