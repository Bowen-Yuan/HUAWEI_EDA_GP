#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$ROOT_DIR/output"

mkdir -p "$OUTPUT_DIR"
python "$ROOT_DIR/check_env.py"

cd "$ROOT_DIR/suanfa"
python -u compare.py 2>&1 | tee "$OUTPUT_DIR/train.log"

python "$ROOT_DIR/plot_metrics.py" \
    --input "$OUTPUT_DIR/convergence_metrics.csv" \
    --output "$OUTPUT_DIR/convergence_curves.png"

cp RBSM_sparse_placement.npy "$OUTPUT_DIR/RBSM_sparse_placement.npy"
if [[ -d prof_log ]]; then
    mkdir -p "$OUTPUT_DIR/prof_log"
    cp -a prof_log/. "$OUTPUT_DIR/prof_log/"
fi

python - <<'PY'
from pathlib import Path
import numpy as np

result = Path("../output/RBSM_sparse_placement.npy").resolve()
positions = np.load(result)
print(f"Saved: {result}")
print(f"Shape: {positions.shape}")
print(f"Dtype: {positions.dtype}")
print(f"Finite: {np.isfinite(positions).all()}")
PY
