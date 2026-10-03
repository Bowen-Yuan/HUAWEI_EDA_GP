# Adaptive-penalty Pareto experiment

## Purpose

This is an independent experiment from the earlier overflow-hard-constraint
work. It starts from the previous 25%-cap solution, alternates short
single-objective blocks, and then switches to a soft-penalty joint descent.
No candidate is rejected because its overflow exceeds a threshold.

The joint objective is

\[
F(x)=H(x)+\lambda(x)D(x),\qquad
\lambda(x)=\frac{\|g_H\|}{\max(\|g_D\|,10^{-12})}
\exp\!\left[0.7\,\mathrm{clip}\!\left(\frac{O(x)-25}{5},-4,4\right)\right].
\]

Here, `H` is exact HPWL, `D` is the shared oracle's linear overflow area, and
`O=100D/A_available`. The gradient-norm ratio makes the two gradient terms
commensurate; the exponential term raises density pressure above the 25%
soft target and relaxes it below the target. Each move is restricted by a
component-wise trust radius of 4--16 row heights. A bounded archive retains
the nondominated points observed during the run.

## Reproducible run

Command:

```powershell
$env:PYTHONPATH = 'D:\codex_project\HUAWEI_EDA\code'
python adaptive_pareto_soft\src\pareto_continuation.py `
  --aux dataset\adaptec1\adaptec1.aux `
  --placement code-hard\runs_alternating_25\adaptec1\20260720T144806Z\solution.pl `
  --warm-cycles 3 --joint-steps 520 --target-overflow-percent 25 `
  --output-root adaptive_pareto_soft\results
```

Output: [20260720T152105Z](../results/adaptec1/20260720T152105Z).
Runtime was 195.17 s for 955 recorded iterations.

| State | HPWL | Overflow |
|---|---:|---:|
| Input 25%-cap solution | 1,020.438 M | 23.498% |
| After pure HPWL block 1 | 417.386 M | 29.933% |
| After pure density block 1 | 1,768.869 M | 0.006% |
| After pure HPWL block 2 | 772.937 M | 16.986% |
| After pure density block 2 | 1,847.691 M | 0.004% |
| After pure HPWL block 3 | 798.736 M | 18.133% |
| After pure density block 3 | 798.795 M | 16.708% |
| Final adaptive joint descent | **88.538 M** | **18.706%** |

The soft target range was 20--30%. The final point is slightly below it, not
because of a constraint, but because the observed joint descent found a
strictly better low-HPWL state there. The archive contains no nondominated
20--30% point with lower HPWL than this final point; it is therefore also the
reported target representative. This is an *empirical* Pareto result for the
visited trajectory, not a proof of global Pareto optimality.

## Visual evidence

![Alternating and joint convergence](../results/adaptec1/20260720T152105Z/figures/adaptive_convergence.png)

![Observed nondominated frontier](../results/adaptec1/20260720T152105Z/figures/pareto_frontier.png)

The sharp lambda spike near the second density phase occurs where the exact
density gradient is nearly zero; the denominator safeguard keeps the run
finite, and the component trust region limits the corresponding displacement.

## Files

- [pareto_continuation.py](pareto_continuation.py): independent solver.
- [solution_final.pl](../results/adaptec1/20260720T152105Z/solution_final.pl): final joint solution.
- [solution_pareto_target.pl](../results/adaptec1/20260720T152105Z/solution_pareto_target.pl): soft-target representative.
- [iterations.csv](../results/adaptec1/20260720T152105Z/iterations.csv): complete trajectory, lambda, and accept/reject records.
- [pareto_archive.csv](../results/adaptec1/20260720T152105Z/pareto_archive.csv): retained nondominated samples.

## Scope note

This implementation intentionally calls the shared `gpplacer` exact oracle,
so its density semantics and its previously identified blocked-bin caveat are
unchanged. This directory changes the search controller only; it does not
alter the common parser or objective implementation.
