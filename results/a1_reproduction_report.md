# adaptec1 (a1) RBSM reproduction report

Date: 2026-07-17. Implementation: the independent `rbsm_place` package in this
directory, run with `D:\anaconda\envs\deeplearn\python.exe`, PyTorch 2.7.0+cu118,
and an RTX 3050 Laptop GPU (4GB).

## Paper target

Table 5 reports the pre-legalization RBSM result for adaptec1 as HPWL `5.05e7`,
density overflow `26.8%`, 85 epochs, and overlap `10.36%`.

The evaluator here uses pin-offset HPWL, bin excess divided by usable row capacity,
and movable/movable rectangle overlap divided by total movable area. Fixed terminals
are accounted for as density blockages. Density is exact for the stated bin grid.
When a grid bucket is pathologically dense, pair-overlap is a deterministic bounded
sample estimate rather than an exhaustive O(n^2) count; it is therefore diagnostic,
not the deciding criterion below.

## Executed experiments

| Run | Initialization / key setting | Epochs | HPWL | Density overflow | Pair overlap | Result |
|---|---|---:|---:|---:|---:|---|
| `smoke4` | random, 1 inner step | 1 | 2.214e9 | 16.71% (64 bins) | 16.26% | end-to-end CUDA smoke pass |
| `warm-full-seed2026` | existing `.pl` warm start, alpha=5 | 85 | 4.761e7 | 57.24% (128 bins) | 5.34e6% sample estimate | HPWL target met; density failed |
| `density-full-seed2026` | warm start, alpha=0, 200k local pairs/refresh | 85 | 5.587e7 | 53.1744% (128 bins) | 3.08e5% sample estimate | failed |
| `density-full-seed2027` | same configuration | 85 | 5.591e7 | 53.1743% (128 bins) | 3.10e5% sample estimate | failed |
| `density-full-seed2028` | same configuration | 85 | 5.609e7 | 53.1777% (128 bins) | 2.27e5% sample estimate | failed |

The three density-focused 85-epoch runs have median HPWL `5.591e7` and median density
overflow `53.1744%`. Thus, the current implementation can push HPWL below the paper
value in a separate alpha=5 experiment, but cannot satisfy the reported density-overflow
target. The density configuration was selected after a 15-epoch tuning run lowered
overflow only from about 57% to 53%; the three completed seeds confirm that this gap is
stable rather than a random-seed accident.

## Conclusion

**The paper's adaptec1 Table 5 indicators were not reproduced on this implementation
and hardware.** The blocking metric is density overflow, not HPWL. The paper omits
several ISPD-specific implementation details (bin configuration, sampling temperature,
initialization/capacity spreading, and exact overlap evaluator), so this is evidence
about the documented algorithm rather than a claim that the authors' implementation is
incorrect.

## Re-run commands

```powershell
Set-Location D:\codex_project\HUAWEI_EDA\dai-code
& D:\anaconda\envs\deeplearn\python.exe -m unittest discover -s tests -v
& D:\anaconda\envs\deeplearn\python.exe rbsm.py solve `
  --aux ..\dataset\adaptec1\adaptec1.aux `
  --config configs\adaptec1_density_full.json `
  --output-dir runs\density-full-seed2026
```

For a strict formula-only run, use `configs/adaptec1_paper_equations.json`. It uses
random initialization and leaves the paper's stated raw update behavior unmodified; it
is intentionally separate from the resolved GPU-scalable configuration.
