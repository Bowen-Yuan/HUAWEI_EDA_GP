# Independent RBSM reproduction

This directory is a standalone implementation of Yue, Dai, and Yu's random-batch
splitting method (RBSM) for nonsmooth VLSI global placement. It reads ISPD2005
Bookshelf files directly and does not import the repository's older `gpplacer` code.

Run from this directory with the available CUDA environment:

```powershell
conda run -n deeplearn python rbsm.py inspect --aux ..\dataset\adaptec1\adaptec1.aux
conda run -n deeplearn python rbsm.py solve --aux ..\dataset\adaptec1\adaptec1.aux --config configs\adaptec1_smoke.json --output-dir runs\smoke
conda run -n deeplearn python rbsm.py evaluate --aux ..\dataset\adaptec1\adaptec1.aux --placement runs\smoke\solution.pl --bins 512
```

`adaptec1_paper_equations.json` retains the paper's stated parameters and random
initialization. `adaptec1_resolved.json` makes the unresolved scalability choices
explicit: connectivity warm start, stable weighted sampling, candidate-pair grid and
a maximum physical displacement per update. The resolved cap prevents the stated
`lr=0.1, gamma=1000` combination from moving a standard cell hundreds of sites.

The primary report compares pre-legalization HPWL, bin density overflow and exact
pairwise rectangle overlap with Table 5: 5.05e7, 26.8%, and 10.36%, respectively.
