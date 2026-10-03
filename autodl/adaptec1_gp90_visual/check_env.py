from pathlib import Path
import sys

import numpy as np
import torch


ROOT = Path(__file__).resolve().parent
CASE_DIR = ROOT / "ispd2005" / "adaptec1"
REQUIRED_FILES = [
    CASE_DIR / "adaptec1.nodes",
    CASE_DIR / "adaptec1.nets",
    CASE_DIR / "adaptec1.pl",
    CASE_DIR / "adaptec1.scl",
    ROOT / "suanfa" / "compare.py",
]


def system_memory_gib():
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return None
    for line in meminfo.read_text(encoding="ascii").splitlines():
        if line.startswith("MemTotal:"):
            return int(line.split()[1]) / 1024**2
    return None


missing = [str(path) for path in REQUIRED_FILES if not path.is_file()]
if missing:
    raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))

if not torch.cuda.is_available():
    raise RuntimeError("CUDA is unavailable. Select a GPU PyTorch image in AutoDL.")

if not hasattr(torch.Tensor, "scatter_reduce_"):
    raise RuntimeError("PyTorch is too old: scatter_reduce_ is required.")

try:
    # Exercise the exact CUDA and AMP interfaces used by compare.py.
    values = torch.tensor([3.0, 1.0], device="cuda")
    indices = torch.tensor([0, 0], dtype=torch.long, device="cuda")
    reduced = torch.full((1,), float("inf"), device="cuda")
    reduced.scatter_reduce_(0, indices, values, reduce="amin")
    torch.amp.GradScaler("cuda", enabled=False)
    with torch.amp.autocast("cuda", dtype=torch.float16):
        amp_result = values * values
    torch.cuda.synchronize()
    if reduced.item() != 1.0 or not torch.isfinite(amp_result).all():
        raise RuntimeError("CUDA API self-test returned an invalid result.")
except (AttributeError, TypeError, RuntimeError) as exc:
    raise RuntimeError(
        "The selected PyTorch image does not support the CUDA/AMP APIs used by "
        "compare.py. Select an official PyTorch 2.3 or newer image."
    ) from exc

device = torch.cuda.get_device_properties(0)
gpu_memory = device.total_memory / 1024**3
host_memory = system_memory_gib()

print(f"Python: {sys.version.split()[0]}")
print(f"PyTorch: {torch.__version__}")
print(f"NumPy: {np.__version__}")
print(f"CUDA runtime: {torch.version.cuda}")
print(f"GPU: {device.name}")
print(f"GPU memory: {gpu_memory:.1f} GiB")
if host_memory is not None:
    print(f"System memory: {host_memory:.1f} GiB")

if gpu_memory < 12:
    print("WARNING: less than 12 GiB GPU memory may cause severe paging or OOM.")
if host_memory is not None and host_memory < 24:
    print("WARNING: less than 24 GiB system memory may fail during dense net preprocessing.")

print("Input files: OK")
print("Environment check: OK")
