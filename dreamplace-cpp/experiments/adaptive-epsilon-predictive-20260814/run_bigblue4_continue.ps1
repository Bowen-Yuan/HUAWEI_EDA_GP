$ErrorActionPreference = "Stop"
$env:Path = "D:\MingGW\ucrt64\bin;" + $env:Path
$project = "D:\codex_project\HUAWEI_EDA\dreamplace-cpp"
$exe = Join-Path $project "dreamplace_cpp.exe"
$bench = "D:\codex_project\HUAWEI_EDA\alg-electronic\ispd2005\bigblue4\bigblue4"
$init = Join-Path $project "output\adaptive_epsilon_ispd2005_20260814\bigblue4\global.pl"
$out = Join-Path $project "output\adaptive_epsilon_predictive_20260814\bigblue4_continue"
New-Item -ItemType Directory -Force $out | Out-Null
& $exe $bench --initial-pl $init --exact-hpwl --optimizer adam --refine-optimizer adam --iterations 1600 --bins 2048 --active-set-radius 370 --refine-active-set-radius 150 --active-set-power 4 --step-fraction 0.0015 --refine-start-overflow 0.15 --feasible-refinement --adaptive-active-predictive --adaptive-active-min-scale 0.75 --adaptive-active-max-scale 1.20 --adaptive-active-window 50 --adaptive-active-interval 25 --adaptive-active-gain 0.35 --adaptive-active-deadband 0.0025 --adaptive-active-max-step 0.015 --dreamplace-detailed --dp-passes 5 --legal-refine-rounds 5 --cell-insertion-passes 2 --legal-projected-passes 4 --row-relegalization-passes 2 --independent-set-size 8 --hungarian-matching --log-every 100 --output $out 2>&1 | Tee-Object (Join-Path $out "run.log")
if ($LASTEXITCODE -ne 0) { throw "bigblue4 continuation failed: $LASTEXITCODE" }
