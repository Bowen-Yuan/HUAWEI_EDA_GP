"""Run the unmodified ISPD reference evaluators.

This module deliberately does *not* reimplement the contest metrics.  Results
published by this project are obtained by launching the Perl scripts supplied
with the benchmark/contest distribution and parsing their stdout.  This keeps
the reported value independent from the differentiable density surrogate used
inside the placers.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re
import shutil
import subprocess


@dataclass(frozen=True, slots=True)
class OfficialEvaluation:
    """Metrics emitted by the ISPD 2005 HPWL and ISPD 2006 density scripts."""

    hpwl: int
    density_target: float
    bins: int
    bins_x: int
    bins_y: int
    violation_bins: int
    violation_fraction: float
    average_overflow: float
    maximum_overflow: float
    overflow_per_bin: float
    total_overflow_amount: float
    scaled_overflow_per_bin: float

    def to_dict(self) -> dict[str, float | int]:
        return asdict(self)


def _bookshelf_files(aux: str | Path) -> dict[str, Path]:
    aux_path = Path(aux).resolve()
    for raw in aux_path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("UCLA") or ":" not in line:
            continue
        _, names = line.split(":", 1)
        files = {Path(name).suffix: aux_path.parent / name for name in names.split()}
        required = {".nodes", ".nets", ".pl", ".scl"}
        if required <= set(files):
            return files
    raise ValueError(f"{aux_path} does not declare a Bookshelf RowBasedPlacement input set")


def _run(perl: str, script: Path, arguments: list[Path | str | float]) -> str:
    executable = shutil.which(perl) if not Path(perl).exists() else perl
    if executable is None:
        raise RuntimeError(
            f"Perl executable {perl!r} was not found. Supply --perl pointing to a Perl runtime; "
            "the official ISPD scripts cannot be replaced by an internal evaluator."
        )
    completed = subprocess.run(
        [str(executable), str(script), *(str(item) for item in arguments)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"Official evaluator failed ({script.name}):\n{completed.stdout}")
    return completed.stdout


def _match(pattern: str, output: str, label: str) -> tuple[str, ...]:
    found = re.search(pattern, output, flags=re.MULTILINE)
    if found is None:
        raise RuntimeError(f"Could not parse {label} from official evaluator output:\n{output}")
    return found.groups()


def run_official_evaluation(
    aux: str | Path,
    placement: str | Path,
    *,
    density_target: float,
    perl: str = "perl",
    hpwl_script: str | Path | None = None,
    density_script: str | Path | None = None,
) -> OfficialEvaluation:
    """Evaluate a solution by executing the original contest scripts.

    The ISPD 2006 density checker is defined only for a target density in
    ``(0, 1)``.  The target is therefore intentionally mandatory rather than
    inheriting the project's historical ``0.8`` surrogate default.
    """
    if not 0.0 < density_target < 1.0:
        raise ValueError("density_target must lie strictly between 0 and 1 for the ISPD 2006 checker")
    root = Path(__file__).resolve().parents[2]
    hpwl = Path(hpwl_script) if hpwl_script else root / "dataset" / "hpwl.pl" / "hpwl.pl"
    density = (Path(density_script) if density_script else
               root / "dataset" / "official_ispd2006" / "check_density_target.pl")
    if not hpwl.is_file() or not density.is_file():
        raise FileNotFoundError("Official ISPD evaluator script is missing from dataset/")
    files = _bookshelf_files(aux)
    hpwl_out = _run(perl, hpwl, [files[".nodes"], files[".pl"], Path(placement), files[".nets"]])
    density_out = _run(perl, density, [files[".nodes"], Path(placement), files[".scl"], density_target])
    hpwl_value = int(_match(r"Total HPWL:\s*([0-9]+)", hpwl_out, "HPWL")[0])
    bins, bins_x, bins_y = _match(r"Total\s+(\d+)\s+\((\d+)\s+x\s+(\d+)\)\s+bins", density_out, "bin count")
    violations, fraction, average, maximum = _match(
        r"Violation num:\s*(\d+)\s+\(([0-9.eE+-]+)\)\s+Avg overflow:\s*([0-9.eE+-]+)\s+Max overflow:\s*([0-9.eE+-]+)",
        density_out, "density violations",
    )
    per_bin, total = _match(
        r"Overflow per bin:\s*([0-9.eE+-]+)\s+Total overflow amount:\s*([0-9.eE+-]+)",
        density_out, "overflow amount",
    )
    scaled = _match(r"Scaled Overflow per bin:\s*([0-9.eE+-]+)", density_out, "scaled overflow")[0]
    return OfficialEvaluation(
        hpwl=hpwl_value, density_target=density_target, bins=int(bins), bins_x=int(bins_x), bins_y=int(bins_y),
        violation_bins=int(violations), violation_fraction=float(fraction), average_overflow=float(average),
        maximum_overflow=float(maximum), overflow_per_bin=float(per_bin),
        total_overflow_amount=float(total), scaled_overflow_per_bin=float(scaled),
    )
