"""Regression tests for the supplied ISPD2005 Bookshelf data."""

from __future__ import annotations

from pathlib import Path

from gpplacer.io.bookshelf import load_bookshelf


ROOT = Path(__file__).resolve().parents[2]


def test_adaptec1_metadata() -> None:
    db = load_bookshelf(ROOT / "dataset" / "adaptec1" / "adaptec1.aux")
    assert (db.node_count, int(db.fixed.sum()), db.net_count, db.pin_count, len(db.rows)) == (
        211447, 543, 221142, 944053, 890,
    )
    assert db.core_bounds == (459.0, 459.0, 11151.0, 11139.0)


def test_adaptec3_metadata() -> None:
    db = load_bookshelf(ROOT / "dataset" / "adaptec3" / "adaptec3.aux")
    assert (db.node_count, int(db.fixed.sum()), db.net_count, db.pin_count, len(db.rows)) == (
        451650, 723, 466758, 1875039, 1944,
    )
    assert db.core_bounds == (36.0, 58.0, 23226.0, 23386.0)
