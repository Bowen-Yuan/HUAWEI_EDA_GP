"""Regression tests for parsing output from the unmodified ISPD scripts."""

from __future__ import annotations

from gpplacer.official import _match


def test_parse_official_hpwl_and_density_fields() -> None:
    hpwl = "Phase 3: Net file processing is done.\n         Total HPWL: 104924229\n"
    density = (
        "\tTotal 12 (3 x 4) bins. Target density: 0.600000\n"
        "\tViolation num: 2 (0.166667)\tAvg overflow: 0.250000\tMax overflow: 0.300000\n"
        "\tOverflow per bin: 4.000000\tTotal overflow amount: 48.000000\n"
        "\tScaled Overflow per bin: 0.000123\n"
    )
    assert _match(r"Total HPWL:\s*([0-9]+)", hpwl, "HPWL") == ("104924229",)
    assert _match(r"Total\s+(\d+)\s+\((\d+)\s+x\s+(\d+)\)\s+bins", density, "bins") == ("12", "3", "4")
    assert _match(r"Scaled Overflow per bin:\s*([0-9.eE+-]+)", density, "scaled") == ("0.000123",)
