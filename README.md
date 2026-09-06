# nonsmooth-gp

C++17 research code for exact, non-smooth global placement. Reported HPWL is weighted pin-offset max-minus-min HPWL; density uses exact rectangle/bin overlap. Smooth wirelength and density surrogates are not part of the challenge pipeline.

Default data root: `D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005`. Data is external and is never committed.

Build with `cmake -S . -B build` then `cmake --build build --config Release`. Run the small default-chain check with `build/Release/nsgp.exe run --case adaptec1 --pipeline framework/params/pipelines/smoke.json --threads 8` (or the corresponding single-config executable path). Audit a checkpoint with `nsgp audit --case adaptec1 --placement selected.pl`.

`run` defaults to one case. The complete benchmark set only runs through `nsgp batch --all-cases --pipeline ...`. Each module boundary writes a fresh-audited `selected.pl`, CSV summary, and manifest under `framework/results/`.

The project uses a compact static module registry: `layout_init`, `hpwl_adam`, `exact_joint_gp`, `exact_recovery`, `surplus_bisection`, and `equal_shape_swap`. `historical_dct_poisson` remains deliberately disabled and is rejected by `challenge_nonsmooth` pipelines.
