"""Interpolate movable coordinates between two UCLA .pl files for staged tests."""
from pathlib import Path
import re
import sys

if len(sys.argv) != 4:
    raise SystemExit("usage: interpolate_pl.py base.pl target.pl alpha")
base_path, target_path, alpha_text = sys.argv[1:]
alpha = float(alpha_text)
line_re = re.compile(r"^(\S+)(\s+)([-+0-9.eE]+)(\s+)([-+0-9.eE]+)(.*)$")

def read(path):
    out = {}
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines(True)
    for line in lines:
        m = line_re.match(line.rstrip("\r\n"))
        if m:
            out[m.group(1)] = (float(m.group(3)), float(m.group(5)), line)
    return lines, out

base_lines, base = read(base_path)
_, target = read(target_path)
for i, line in enumerate(base_lines):
    m = line_re.match(line.rstrip("\r\n"))
    if not m or m.group(1) not in target:
        continue
    name = m.group(1)
    bx, by, _ = base[name]
    tx, ty, _ = target[name]
    # Fixed nodes are marked with /FIXED or /FIXED_NI in the suffix.
    suffix = m.group(6)
    if "FIXED" in suffix.upper():
        continue
    x = (1.0 - alpha) * bx + alpha * tx
    y = (1.0 - alpha) * by + alpha * ty
    newline = f"{name}{m.group(2)}{x:.10f}{m.group(4)}{y:.10f}{suffix}\n"
    base_lines[i] = newline
sys.stdout.write("".join(base_lines))
