from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rbsm_place.bookshelf import Design
from rbsm_place.metrics import hpwl, overlap_pairs
from rbsm_place.solver import RBSMSolver, SolverConfig


class MetricsTest(unittest.TestCase):
    def tiny(self) -> Design:
        return Design(Path("tiny.aux"), ("a", "b"), np.array([2., 2.]), np.array([2., 2.]), np.array([False, False]),
                      np.array([[1., 1.], [4., 1.]]), np.array([0, 1]), np.zeros((2, 2)), np.array([0, 2]), np.ones(1), ((0., 0., 10., 10.),))

    def test_hpwl(self):
        d = self.tiny(); self.assertEqual(hpwl(d, d.centres), 3.0)

    def test_overlap(self):
        d = self.tiny(); c = d.centres.copy(); c[1] = (1.5, 1.)
        area, count = overlap_pairs(d, c)
        self.assertAlmostEqual(area, 3.0); self.assertEqual(count, 1)

    def test_resolved_step_is_bounded(self):
        d = self.tiny()
        solver = RBSMSolver(d, SolverConfig(epochs=1, inner_steps=1, batch_fraction=1,
                           device="cpu", initialization="bookshelf", noise_scale=0, max_displacement=2))
        solver.run(ROOT / "tmp_test")
        self.assertLessEqual(float(np.abs(solver.positions.detach().numpy() - d.centres).max()), 4.01)


if __name__ == "__main__":
    unittest.main()
