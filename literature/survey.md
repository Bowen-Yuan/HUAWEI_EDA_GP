# Scoping literature review: exact/non-smooth density directions in VLSI placement

Generated: 2026-08-23  
Review type: scoping review  
Question: Which ideas from analytical global placement and nonsmooth constrained
optimization can guide an exact rectangle-overlap, non-smoothed HPWL solver?

## Search protocol

Databases: Crossref REST API and Semantic Scholar REST API. Queries were
`DREAMPlace Deep Learning Toolkit Enabled GPU Acceleration for Modern VLSI
Placement`, `RePlAce Advancing Solution Quality and Routability Validation in
Global Placement`, `ePlace Electrostatics based Placement`, `global placement
exact overlap constraint ADMM`, and `nonsmooth optimization exact penalty
active set`. Search date: 2026-08-23. Inclusion: primary technical papers with
placement, density/congestion, or nonsmooth constrained-optimization methods.
Exclusion: legalization-only methods, smooth surrogate-only claims, and papers
without a verifiable DOI or publisher record.

## Evidence summary

| Study | DOI | Relevant evidence | Implication here |
|---|---|---|---|
| Lu et al., DREAMPlace (TCAD 2021) | 10.1109/TCAD.2020.3003843 | GPU/autodiff implementation of analytical placement; electrostatic density is a smooth computational model | Useful runtime and continuation baseline, but its electrostatic term cannot remain in our final objective |
| Jiang et al., RePlAce (TCAD 2019) | 10.1109/TCAD.2018.2859220 | Routability-aware global placement and density control benefit from staged parameter schedules | Supports stage-level controller changes; does not justify smoothing exact overlap |
| Lu et al., ePlace-MS (TCAD 2015) | 10.1109/TCAD.2015.2391263 | Electrostatics provides long-range global movement and FFT acceleration | Explains why exact local overlap is sparse; only use as historical motivation, not final cost |
| Curtis & Que, nonsmooth BFGS-SQP (2017) | 10.1080/10556788.2016.1208749 | Nonsmooth constrained optimization can use safeguarded quasi-Newton/SQP directions with line-search acceptance | Motivates exact-audited direction proposals and active-set restarts |
| Distributed nonsmooth optimization with exact penalty (CDC 2023) | 10.1109/CDC49753.2023.10383876 | Exact penalties and distributed block directions can preserve nonsmooth constraints without smoothing | Supports block/group directions, but our overlap audit remains the authoritative acceptance test |

## Synthesis and gap

The placement literature overwhelmingly uses smooth density/electrostatic
surrogates, while nonsmooth optimization literature emphasizes safeguarded
directions and exact penalties rather than a rectangle-grid overlap oracle.
This leaves a gap directly relevant to the project: a net-aware, capacity-aware
candidate generator whose every accepted state is checked by exact HPWL and
exact overlap. H297 is consistent with this gap: sharing exact density
directions addresses active-set starvation, but needs a cumulative HPWL guard
and fixed-macro-aware capacity assignment to avoid topology debt.

## Search log

| Database | Query | Result |
|---|---|---:|
| Crossref | DREAMPlace title | 3 records; first DOI verified |
| Crossref | RePlAce title | 3 records; first DOI verified |
| Crossref | ePlace title | 3 records; first DOI verified |
| Crossref | exact overlap/ADMM placement | mostly unrelated constraint programming; no direct exact-overlap VLSI method |
| Semantic Scholar | nonsmooth optimization exact penalty active set | rate-limited for placement queries; nonsmooth exact-penalty records returned |
