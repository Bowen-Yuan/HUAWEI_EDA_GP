"""Exact HPWL value and deterministic Clarke-subgradient implementation."""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def hpwl_value_gradient(
    centres: np.ndarray, pin_node: np.ndarray, pin_offset: np.ndarray,
    net_start: np.ndarray, net_weight: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return weighted HPWL, a deterministic subgradient, and net extrema.

    For tied extrema, each tied pin receives equal coefficient.  ``extrema``
    stores the first pin index at min-x, max-x, min-y, max-y for inexpensive
    active-configuration change detection.
    """
    node_count = centres.shape[0]
    net_count = net_weight.size
    gradient = np.zeros((node_count, 2), dtype=np.float64)
    extrema = np.empty((net_count, 4), dtype=np.int64)
    total = 0.0
    for net in range(net_count):
        first = net_start[net]
        last = net_start[net + 1]
        min_x = 1.0e300
        max_x = -1.0e300
        min_y = 1.0e300
        max_y = -1.0e300
        min_x_pin = first
        max_x_pin = first
        min_y_pin = first
        max_y_pin = first
        for pin in range(first, last):
            node = pin_node[pin]
            x = centres[node, 0] + pin_offset[pin, 0]
            y = centres[node, 1] + pin_offset[pin, 1]
            if x < min_x:
                min_x, min_x_pin = x, pin
            if x > max_x:
                max_x, max_x_pin = x, pin
            if y < min_y:
                min_y, min_y_pin = y, pin
            if y > max_y:
                max_y, max_y_pin = y, pin
        extrema[net, 0] = min_x_pin
        extrema[net, 1] = max_x_pin
        extrema[net, 2] = min_y_pin
        extrema[net, 3] = max_y_pin
        total += net_weight[net] * ((max_x - min_x) + (max_y - min_y))
        count_min_x = 0
        count_max_x = 0
        count_min_y = 0
        count_max_y = 0
        for pin in range(first, last):
            node = pin_node[pin]
            x = centres[node, 0] + pin_offset[pin, 0]
            y = centres[node, 1] + pin_offset[pin, 1]
            if x == min_x:
                count_min_x += 1
            if x == max_x:
                count_max_x += 1
            if y == min_y:
                count_min_y += 1
            if y == max_y:
                count_max_y += 1
        weight = net_weight[net]
        for pin in range(first, last):
            node = pin_node[pin]
            x = centres[node, 0] + pin_offset[pin, 0]
            y = centres[node, 1] + pin_offset[pin, 1]
            if x == min_x:
                gradient[node, 0] -= weight / count_min_x
            if x == max_x:
                gradient[node, 0] += weight / count_max_x
            if y == min_y:
                gradient[node, 1] -= weight / count_min_y
            if y == max_y:
                gradient[node, 1] += weight / count_max_y
    return total, gradient, extrema
