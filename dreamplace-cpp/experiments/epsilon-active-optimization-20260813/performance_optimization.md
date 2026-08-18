# Epsilon-active performance optimization

Date: 2026-08-13

## Benchmark

The controlled benchmark is `adaptec1` with the same exact-HPWL path used by
the existing epsilon-active experiment:

```text
--exact-hpwl --optimizer adam --iterations 1400 --bins 512
--active-set-radius 125 --refine-active-set-radius 50
--active-set-power 4 --step-fraction 0.003
--feasible-refinement --no-detailed --no-abacus --profile
```

The process used `OMP_NUM_THREADS=16` on the existing Intel Core i5-13500H
machine. No smooth HPWL surrogate, bundle mode, reduced grid, or reduced
iteration budget was enabled.

## Changes retained

1. The epsilon-active wirelength kernel uses cached flat pin arrays, cached
   net weights and CSR-style node-to-pin incidence. The exact objective and
   active-set direction are unchanged.
2. Position capture/apply uses one parallel scan instead of four separate
   OpenMP regions. Boundary clamping is fused into apply.
3. Movable/filler gradient assembly is indexed directly and parallelized.
4. Adam bias-correction denominators are computed once per iteration instead
   of recomputing `pow(beta, iteration)` for every coordinate.
5. Filler effective geometry is cached after initialization.
6. In the finite-difference electric-field path, the transformed potential
   buffer is swapped instead of copied before IDCT.
7. Overflow, density energy, L1 norms and RMS reductions use OpenMP reductions.
8. Tangent-refinement work buffers are allocated once and reused.

These changes preserve the command-line interface and the numerical model.
The adaptive epsilon controller remains opt-in and disabled by default.

## Rejected experiments

Thread-private grid deposition and atomic grid deposition were tested. On this
Windows/MinGW build they were slower or statistically indistinguishable from
the serial rectangle deposition because each rectangle touches only a small
number of bins and the extra grid merge/atomic contention dominates. They are
not retained.

OpenMP affinity settings were also tested. `libgomp` reports that affinity is
not supported in the current configuration, and the measured difference was
negligible.

## Measurements

| version | iterations | wall time | loop time | density | wirelength | other | GP HPWL | GP overflow | legal HPWL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| before optimization (recorded) | 1400 | 271.6 s | - | - | - | - | 72.879 M | 6.996% | GP-only |
| optimized, 100-step profile | 100 | 4.32 s | 4.32 s | 2.54 s | 0.87 s | 0.91 s | 43.598 M | 97.207% | diagnostic |
| optimized final run | 1400 | 136.1 s | 131.4 s | 75.35 s | 29.68 s | 26.41 s | 72.904 M | 6.995% | 77.557 M |

The end-to-end speedup is approximately `2.06x` while the final placement
quality and internal legality remain in the same range as the original run:

```text
GP HPWL       72.9041 M
GP overflow    6.99528 %
legal HPWL    77.5569 M
internal legal = yes
```

The remaining dominant cost is the 512x512 finite-difference density solve,
especially its two-dimensional DCT/IDCT and grid deposition. Replacing those
operations with a lower-resolution grid or a different density approximation
would no longer be a pure performance optimization and therefore was not
enabled in the default path.
