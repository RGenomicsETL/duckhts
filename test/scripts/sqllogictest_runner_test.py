"""Exercise the real runner's process status for valid and broken test inputs."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class RunnerStatusTest(unittest.TestCase):
    def test_process_status(self):
        root = Path(__file__).resolve().parents[2]
        runner = root / "scripts" / "run_sqllogictest.py"
        with tempfile.TemporaryDirectory(prefix="duckhts-sqllogic-gate-") as directory:
            folder = Path(directory)
            valid = folder / "valid.test"
            broken = folder / "broken.test"
            listing = folder / "tests.txt"
            valid.write_text("query I\nSELECT 42;\n----\n42\n", encoding="utf-8")
            listing.write_text(f"{valid}\n{broken}\n", encoding="utf-8")

            def run(*arguments):
                return subprocess.run(
                    [sys.executable, str(runner), "--test-dir", str(folder), *map(str, arguments)],
                    cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                )

            result = run("--file-path", valid)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("SUCCESS", result.stdout)
            for content, diagnostic in (
                ("statement error\nSELECT 42;\n", "Invalid SQLLogicTest"),
                ("query I\nSELECT 42;\n----\n43\n", "ERROR"),
            ):
                broken.write_text(content, encoding="utf-8")
                for arguments in (
                    ("--file-path", broken), ("--file-list", listing), ()
                ):
                    with self.subTest(content=content, arguments=arguments):
                        result = run(*arguments)
                        self.assertNotEqual(result.returncode, 0, result.stdout)
                        self.assertIn(diagnostic, result.stdout)


if __name__ == "__main__":
    unittest.main()
