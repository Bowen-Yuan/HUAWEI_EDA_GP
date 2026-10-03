# AutoDL adaptec1 RBSM GPU package

This package runs the current `compare.py` GPU reference implementation on the
ISPD2005 adaptec1 benchmark. Only filesystem paths were adapted. The numerical
model, optimizer, loss functions, and default parameters were not changed.

Current defaults in `suanfa/compare.py`:

- case: adaptec1
- outer steps: 100
- distance SGD learning rate: 5
- penalty SGD learning rate: 5
- inner steps per outer step: 25

These are the current code defaults, not an exact reconstruction of the paper's
Table 5 configuration.

## 1. Create the AutoDL instance

AutoDL does not let you freely combine CPU, system memory, Ubuntu, Python,
CUDA, and PyTorch. First select a machine with a fixed hardware configuration,
then select one complete image from the image list.

### Machine selection

- Filter GPU model by `RTX 3090 (24GB)`.
- Among the available 3090 machines, prefer one whose details show at least
  24 GiB of system memory. 32 GiB or more is safer. CPU core count is secondary.
- Do not search for separate 8-core/64-GB selectors: CPU and memory are normally
  fixed properties of each listed machine.
- The default system and data disk capacities are sufficient for this package.

### Image selection

Select an AutoDL official image whose framework is PyTorch. Use the image list
actually shown for the selected machine and apply this priority:

1. PyTorch 2.4 or newer: preferred.
2. PyTorch 2.3: supported and sufficient.
3. PyTorch 2.0-2.2: do not use with the current code because the required
   `torch.amp.GradScaler("cuda")` interface may be absent.
4. PyTorch 1.x or a plain CUDA image: do not select.

CUDA 11.8, 12.1, and 12.4 PyTorch images are all suitable. Python 3.9, 3.10,
3.11, and 3.12 are also suitable. Ubuntu does not need to be selected
separately; it is part of the image. A typical acceptable label is therefore
`PyTorch 2.3+ / Python 3.x / CUDA 11.8+`, but the exact versions do not have to
match this example.

The CUDA version printed by `nvidia-smi` is the driver capability and may be
newer than `torch.version.cuda`; this is normal. Do not reinstall PyTorch in an
official image. Run `check_env.py` to test the APIs that this code actually
uses.

Use `/root/autodl-tmp` for the package and training outputs. It is faster than
running directly from network storage.

## 2. Upload and extract

Upload `autodl_rbsm_adaptec1.zip` with the AutoDL file browser, then run:

```bash
cd /root/autodl-tmp
unzip autodl_rbsm_adaptec1.zip
cd autodl_rbsm_adaptec1
```

If the archive was uploaded elsewhere, replace `/root/autodl-tmp` with its
actual location.

## 3. Verify and install the small Python dependency

```bash
nvidia-smi
python --version
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
pip install -r requirements.txt
python check_env.py
```

`check_env.py` must report the RTX 3090, CUDA availability, all input files,
and `Environment check: OK`.

## 4. Run in tmux

Use `tmux` so an SSH disconnect does not terminate training:

```bash
cd /root/autodl-tmp/autodl_rbsm_adaptec1
tmux new -s rbsm
bash run_compare.sh
```

The run uses `RBSM_SEED=2026` and evaluates full HPWL after every inner step by
default. Both settings can be overridden without editing the file:

```bash
RBSM_SEED=2027 RBSM_EVAL_EVERY_STEP=1 bash run_compare.sh
```

Exact full-HPWL evaluation at every inner step can substantially increase the
runtime. Use `RBSM_EVAL_EVERY_STEP=5` or `25` for a faster, coarser curve.

Detach from tmux with `Ctrl-b`, then `d`. Reconnect with:

```bash
tmux attach -t rbsm
```

In another terminal, monitor the run with:

```bash
watch -n 2 nvidia-smi
tail -f /root/autodl-tmp/autodl_rbsm_adaptec1/output/train.log
```

The first preprocessing stage can consume more than 12 GiB of system memory
and may produce no log output for several minutes.

## 5. Results

After a successful run, download the `output` directory. It contains:

- `RBSM_sparse_placement.npy`: `(210904, 2)` movable-cell coordinates
- `convergence_metrics.csv`: full HPWL and approximate density overflow history
- `convergence_curves.png`: HPWL and overflow convergence curves
- `train.log`: complete console log
- `prof_log/`: PyTorch profiler trace, usually tens of megabytes

Check the final output again with:

```bash
python - <<'PY'
import numpy as np
p = np.load('output/RBSM_sparse_placement.npy')
print(p.shape, p.dtype, np.isfinite(p).all())
print('x:', p[:, 0].min(), p[:, 0].max())
print('y:', p[:, 1].min(), p[:, 1].max())
PY
```

AutoDL instances continue billing while powered on. Download the results and
stop the instance after training completes.
