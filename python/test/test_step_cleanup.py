import os
import unittest
from pathlib import Path
from .pytest_common import pushd, create_black_oil_simulator

class TestStepCleanup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1a")

    def test_cleanup_without_step(self):
        # step_cleanup() directly after step_init(), without any step() in
        # between, must not crash
        with pushd(self.data_dir):
            sim = create_black_oil_simulator(
                filename="SPE1CASE1.DATA", args=["--output-dir=cleanup_without_step"])
            sim.step_init()
            self.assertEqual(sim.step_cleanup(), 0)
