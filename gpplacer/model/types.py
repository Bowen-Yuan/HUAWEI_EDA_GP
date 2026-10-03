"""In-memory representation of a Bookshelf row-based placement instance.

Coordinates inside the solver are cell centres.  Bookshelf ``.pl`` files use
lower-left coordinates, so conversion happens only in the parser and writer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class SiteRow:
    """One usable horizontal row segment from a Bookshelf ``.scl`` file."""

    x: float
    y: float
    width: float
    height: float
    site_width: float


@dataclass(slots=True)
class PlacementDB:
    """Compact topology, geometry, and fixed-placement state.

    ``net_start`` and ``node_net_start`` are CSR offsets.  Pin coordinates are
    ``centres[pin_node] + pin_offset``.  Net weights default to one if ``.wts``
    is empty, as in the supplied adaptec cases.
    """

    source_aux: str
    node_names: tuple[str, ...]
    net_names: tuple[str, ...]
    width: FloatArray
    height: FloatArray
    fixed: BoolArray
    initial_centres: FloatArray
    orientation: tuple[str, ...]
    pin_node: IntArray
    pin_offset: FloatArray
    net_start: IntArray
    net_weight: FloatArray
    node_net_start: IntArray
    node_nets: IntArray
    rows: tuple[SiteRow, ...]

    @property
    def node_count(self) -> int:
        return len(self.node_names)

    @property
    def net_count(self) -> int:
        return len(self.net_names)

    @property
    def pin_count(self) -> int:
        return int(self.pin_node.size)

    @property
    def movable(self) -> BoolArray:
        return ~self.fixed

    @property
    def row_height(self) -> float:
        heights = {row.height for row in self.rows}
        if len(heights) != 1:
            raise ValueError("The first implementation requires uniform row heights.")
        return next(iter(heights))

    @property
    def core_bounds(self) -> tuple[float, float, float, float]:
        """Return the bounding rectangle of all legal row segments."""
        return (
            min(row.x for row in self.rows),
            min(row.y for row in self.rows),
            max(row.x + row.width for row in self.rows),
            max(row.y + row.height for row in self.rows),
        )

    def node_index(self, name: str) -> int:
        """Resolve a Bookshelf object name; intended for diagnostics/tests."""
        try:
            return self.node_names.index(name)
        except ValueError as exc:
            raise KeyError(f"Unknown node {name!r}") from exc

    def pin_slice(self, net: int) -> slice:
        """Return the contiguous pin range belonging to ``net``."""
        return slice(int(self.net_start[net]), int(self.net_start[net + 1]))
