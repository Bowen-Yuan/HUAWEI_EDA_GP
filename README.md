# Epsilon-Active Exact Nonsmooth Global Placer

This directory is an independent C++17 global-placement framework for Huawei
EDA Challenge 5. It does not compile or link any source file from
`dreamplace-cpp`.

The optimized model is evaluated exactly:

```text
sum_e weight_e * HPWL_e + lambda *
sum_b 0.5 * bin_area * max(rho_b - target_density, 0)^2
```

- `HPWL_e` is the exact pin-offset `max-min` value.
- `rho_b` is computed from exact rectangle/bin intersections.
- Fixed macros consume bin capacity. Nonphysical `terminal_NI` ports do not.
- Epsilon-active weights choose a direction from nearby nonsmooth pieces; they
  never change the objective value or reported metrics.
- There is no weighted-average, log-sum-exp, Poisson, DCT, legalization, or
  legality-checking path.

## Build

```powershell
mingw32-make -j16
mingw32-make test
```

## Run

```powershell
.\build\epsilon_active.exe `
  --benchmark ..\dataset\adaptec1\adaptec1 `
  --output output\adaptec1 `
  --bins 512 --iterations 1000 --threads 16 `
  --optimizer adam --lambda-policy dreamplace `
  --hpwl-epsilon 125 --active-power 4 --density-active-power 1 `
  --snapshot-every 10
```

The only placement result is `global.pl`. The program intentionally does not
produce a legalized `final.pl`.

Optional exact-overlap active-set escape:

```powershell
--density-coordinate-sweeps 1 `
--density-coordinate-net-blocks `
--density-coordinate-capacity-assignment `
--density-coordinate-hierarchical-clusters `
--density-coordinate-cooperative-batch `
--density-coordinate-heavy-edge-growth `
--density-coordinate-no-node-moves `
--density-coordinate-assignment-bins 16 `
--density-coordinate-hierarchy-levels 3 `
--density-coordinate-cluster-max-nodes 128 `
--density-coordinate-block-degree 16 `
--density-coordinate-block-max-nodes 32 `
--density-coordinate-step-bins 4 `
--density-coordinate-hpwl-budget 0.10
```

This phase performs exact axis/breakpoint coordinate descent. It is a search
heuristic with an exact HPWL budget; it does not smooth the objective or add a
new density term. It is most useful on a coarse grid as an early bridge before
fine-grid epsilon-active GP.

The hierarchical option grows disjoint net-connected clusters in overloaded
coarse regions and translates them rigidly toward live residual capacity.
Every candidate is re-audited with exact HPWL and exact rectangle overlap
before commit; fixed macro occupancy is included in the capacity ledger.
