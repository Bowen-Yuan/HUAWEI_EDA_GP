# nonsmooth-gp

C++17 research code for exact, non-smooth global placement. Reported HPWL is weighted pin-offset max-minus-min HPWL; density uses exact rectangle/bin overlap. Smooth wirelength and density surrogates are not part of the challenge pipeline.

Default data root: `D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005`. Data is external and is never committed.

Build with `cmake -S . -B build` then `cmake --build build --config Release`. Run the small default-chain check with `build/Release/nsgp.exe run --case adaptec1 --pipeline framework/params/pipelines/smoke.json --threads 8` (or the corresponding single-config executable path). Audit a checkpoint with `nsgp audit --case adaptec1 --placement selected.pl`.

`run` defaults to one case. The complete benchmark set only runs through `nsgp batch --all-cases --pipeline ...`.  The runner is a small registry-driven microkernel: it loads data, invokes stage modules, performs a fresh exact audit at each boundary, and writes only `params.json`, `experiment.md`, and `trajectory.csv`.  Add `--save-stage MODULE` only when an audited placement is explicitly needed.

The project uses a compact static module registry.  Each optimization stage owns a `modules/<name>/code/module.cpp` adapter and JSON parameters; shared Bookshelf I/O, exact evaluators, and optimizers live in `framework/kernel`. `historical_dct_poisson` remains deliberately disabled and is rejected by `challenge_nonsmooth` pipelines.

## V3 retained global-view experiments

`global_view_gp` is a normal registry stage and can be placed anywhere in a JSON
pipeline. `nsgp lab` remains a thin, convenient single-stage V3 entry point that
uses the same search implementation. It reads the existing H375
`adaptec1` 7% checkpoint as a read-only external input, records its SHA-256,
and keeps layout hand-offs in memory.  A normal run deliberately produces only
`params.json`, `experiment.md`, and `trajectory.csv` under
`framework/results/experiments/<run-id>`; overflow in all reports is a
percentage.  For example:

```powershell
nsgp lab --case adaptec1 --iterations 50 --optimizer adam --active-ensemble --bundle
```

The equivalent composable form is `nsgp run --case adaptec1 --placement
<checkpoint.pl> --pipeline framework/params/pipelines/global_view_smoke.json`;
edit or add only the module parameter JSON and pipeline when composing a new
numerical experiment.

Use `--save-stage NAME` only when a placement artifact is explicitly needed.
The lab supports `adam`, `amsgrad`, `adagrad`, `heavy-ball`, `sgd`,
`normalized-sgd`, and `dual-averaging`; `--step-policy trust` enables a
trust-radius cap and exact backtracking.  `--capacity-transport` enables the
coarse-grid candidate / canonical-grid exact-audit transport stage.
