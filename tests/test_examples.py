"""Optional graphical example checks; opt in after installing .[examples].

PAD_RUN_EXAMPLES=1 python tests/test_examples.py
Cartopy's coastline data may require network access on the first run.
"""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipUnless(os.environ.get("PAD_RUN_EXAMPLES") == "1",
                     "Set PAD_RUN_EXAMPLES=1 to execute graphical examples.")
class ExampleTests(unittest.TestCase):
    def run_example(self, filename):
        script = Path(__file__).resolve().parents[1] / "examples" / filename
        with tempfile.TemporaryDirectory() as workdir:
            environment = dict(os.environ, MPLBACKEND="Agg", MPLCONFIGDIR=workdir)
            # Running elsewhere also checks that sample-data paths are script-relative.
            result = subprocess.run(
                [sys.executable, str(script)], cwd=workdir, env=environment,
                capture_output=True, text=True, timeout=600,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("PAD distance:", result.stdout)

    def test_synthetic_example(self):
        self.run_example("example_PAD_01.py")

    def test_netcdf_example(self):
        self.run_example("example_PAD_02.py")


if __name__ == "__main__":
    unittest.main(verbosity=2)
