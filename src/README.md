# Solver source

`pareto_continuation.py` is the independent adaptive-penalty solver. It is
separate from the hard-constraint implementations under `code-hard/`.

The solver needs the shared placement parser and exact oracle in `../code`:

```powershell
$env:PYTHONPATH = 'D:\codex_project\HUAWEI_EDA\code'
python adaptive_pareto_soft\src\pareto_continuation.py `
  --aux dataset\adaptec1\adaptec1.aux `
  --placement code-hard\runs_alternating_25\adaptec1\20260720T144806Z\solution.pl `
  --warm-cycles 3 --joint-steps 520 --target-overflow-percent 25 `
  --output-root adaptive_pareto_soft\results
```

See the project [README](../README.md) and
[PROJECT_INTRO.md](../PROJECT_INTRO.md) before modifying parameters.
