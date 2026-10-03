# Electrostatic density variant

This directory preserves the command-line interface of `alg` while replacing
the local exact-overlap density subgradient with a smooth electrostatic model.

The density path is:

1. Stretch small movable cells to at least `sqrt(2) * bin_size` for the density
   footprint while preserving their original deposited area.
2. Deposit movable cells with a compact separable triangular kernel and fixed
   cells with exact rectangle/bin overlap.
3. Apply an orthonormal 2D cosine transform to the density map.
4. Solve the Neumann Poisson equation spectrally, with the DC mode removed.
5. Recover both electric-field components with mixed sine/cosine transforms.
6. Integrate the field over each movable cell footprint to obtain the negative
   electric force, which is the density-energy gradient.

HPWL and density gradients are kept separate and combined as

```text
gradient = hpwl_gradient + lambda_effective * electric_density_gradient
```

The effective coefficient uses DREAMPlace-style initial gradient balancing:

```text
lambda_effective = lambda_control * ||grad(HPWL)||_1 / ||grad(density)||_1
```

`overflow` remains an independent reporting/stopping metric computed from the
smooth density map. It is not used as the differentiable density objective.

## Baseline-driven lambda schedule

Phase 0 remains HPWL-only with lambda fixed at zero. Starting from the first
density phase, the controller evaluates every 50 density steps using

```text
hpwl_ratio     = 7e7 / current_hpwl
overflow_ratio = 0.06 / current_overflow
```

A larger ratio means that metric is relatively closer to its baseline. In the
coarse and medium phases, if the HPWL ratio exceeds the overflow ratio by more
than the 5% hysteresis band, density-control lambda is multiplied by 1.15. If
the overflow ratio is larger by more than 5%, lambda is halved. Otherwise it
is held. Fine-grid updates run every 25 steps and use smooth primal-dual
feedback around 5.98% overflow. Lambda is clamped to `[1e-4, 100]`, and each
decision is recorded in `lambda_strategy.csv`.

Build and run exactly as before:

```powershell
mingw32-make -j4
.\nsp_placer.exe ispd2005/adaptec1/adaptec1
```

The current default keeps HPWL exact in every phase. Phase 0 is strictly
HPWL-only (`lambda=0`) and restores its lowest-HPWL state before spreading.
Density phases use the original max-minus-min HPWL subgradient together with
`p1_lambda=0.02`, per-grid gradient normalization, heavy-ball updates, and smooth
overflow feedback around 5.98%. No LogSumExp wirelength approximation or
gamma annealing is active in the default path. The default budget is 3000
iterations, and the best feasible state (`overflow <= 0.06`) is restored
before output.

The fine-grid feedback uses an exponential primal-dual update. It tracks an
EMA of normalized overflow error, increases lambda when density is infeasible,
and gently decreases lambda when density is feasible but HPWL remains above
the `7e7` reference. Every update is clamped by `--lambda-up` and
`--lambda-down`; lambda remains within `[1e-4, --lambda-max]`.

An optional `--trajectory-lambda` controller was also implemented. It follows
a smoothstep overflow trajectory using PI-D error, overflow descent velocity,
integral anti-windup, and a horizon capped by the remaining fine-step budget.
It produced a more stable overflow path, but its best feasible HPWL was worse
than the continuous baseline controller, so it is retained as an ablation and
is not the default.

The main solver has no HPWL-smoothing switch: every phase always calls the
exact max-minus-min HPWL subgradient. Historical smooth-HPWL CSV files remain
under `experiments/` only as archived ablation data. `--heavy-ball` and
`--no-lambda-guard` remain available for optimizer/controller ablations.

With exact HPWL, the verified default 3000-iteration result is HPWL 79,882,384
and overflow 0.059998; the best repeated run is HPWL 79,422,944 at overflow
0.059985. A separately tuned short-budget run with 150 HPWL-only iterations
followed directly by 650 fine-grid steps gives HPWL 106,027,552 at overflow
0.059873. The short run uses Adam, faster lambda growth, and no extra
late-stage cooldown; it does not reuse the 3000-step optimizer schedule.
See `experiments/EXACT_HPWL_RESULTS.md` for the complete ablation table.

The benchmark inputs are hard-linked from `../alg/ispd2005`; generated output
files are created inside this directory and do not overwrite those inputs.

Run the focused electrostatic test with:

```powershell
mingw32-make test-electric
```
