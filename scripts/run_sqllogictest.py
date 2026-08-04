#!/usr/bin/env python3
"""Run DuckHTS SQLLogicTest files with per-test file cleanup.

The upstream duckdb-sqllogictest runner cleans temporary databases but does not
remove files produced by SQL statements such as COPY, bgzip, index builders,
mosdepth, or idxstats.  This runner adds a narrow teardown layer around that
runner.  It only removes output paths declared by the test SQL, skips tracked
files and fixture directories, and keeps /tmp cleanup limited to DuckHTS-named
files.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Iterable

import duckdb_sqllogictest.python_runner as sqllogic_runner


_WORKING_DIRECTORY = "__WORKING_DIRECTORY__/"
_OUTPUT_ASSIGNMENT = re.compile(
    r"\b(?:output_path|output|index_path|gzi_path|log_path)\s*:=\s*'([^']+)'",
    re.IGNORECASE,
)
_COPY_OUTPUT = re.compile(r"\bTO\s+'([^']+)'", re.IGNORECASE)
_MOSDEPTH_PREFIX = re.compile(r"\bduckhts_mosdepth\s*\(\s*'([^']+)'", re.IGNORECASE)
_MOSDEPTH_SUFFIXES = (
    ".mosdepth.summary.txt",
    ".per-base.bed.gz",
    ".per-base.bed.gz.csi",
    ".regions.bed.gz",
    ".regions.bed.gz.csi",
    ".region.dist.txt",
    ".thresholds.bed.gz",
    ".thresholds.bed.gz.csi",
    ".quantized.bed.gz",
    ".quantized.bed.gz.csi",
)
_SIDECAR_SUFFIXES = (".fai", ".gzi", ".tbi", ".csi", ".bai", ".crai")


def _tracked_files(root: Path) -> set[str]:
    try:
        output = subprocess.check_output(
            ["git", "ls-files", "-z"], cwd=root, stderr=subprocess.DEVNULL
        )
    except (OSError, subprocess.CalledProcessError):
        return set()
    return {item for item in output.decode().split("\0") if item}


def _replace_working_directory(value: str, root: Path) -> Path | None:
    if value.startswith(_WORKING_DIRECTORY):
        return (root / value[len(_WORKING_DIRECTORY) :]).resolve()
    if value.startswith("/tmp/duckhts_"):
        return Path(value).resolve()
    return None


def _allowed_output(path: Path, root: Path, tracked: set[str]) -> bool:
    if path.is_absolute() and path.parent == Path("/tmp"):
        return path.name.startswith("duckhts_")

    try:
        relative = path.relative_to(root)
    except ValueError:
        return False

    relative_name = relative.as_posix()
    if relative_name in tracked:
        return False
    if relative_name == "test/data" or relative_name.startswith("test/data/"):
        return False
    if relative_name == "test/duckvep" or relative_name.startswith("test/duckvep/"):
        return False
    if relative_name == "third_party" or relative_name.startswith("third_party/"):
        return False
    return True


def _declared_artifacts(test_path: Path, root: Path, tracked: set[str]) -> set[Path]:
    try:
        sql = test_path.read_text(encoding="utf-8")
    except OSError:
        return set()

    values = set(_OUTPUT_ASSIGNMENT.findall(sql))
    values.update(_COPY_OUTPUT.findall(sql))
    mosdepth_prefixes = set(_MOSDEPTH_PREFIX.findall(sql))

    artifacts: set[Path] = set()
    for raw in values:
        path = _replace_working_directory(raw, root)
        if path is None or not _allowed_output(path, root, tracked):
            continue
        artifacts.add(path)
        artifacts.update(Path(str(path) + suffix) for suffix in _SIDECAR_SUFFIXES)

    for raw in mosdepth_prefixes:
        prefix = _replace_working_directory(raw, root)
        if prefix is None or not _allowed_output(prefix, root, tracked):
            continue
        artifacts.update(Path(str(prefix) + suffix) for suffix in _MOSDEPTH_SUFFIXES)

    return artifacts


def _cleanup(paths: Iterable[Path]) -> None:
    if os.environ.get("DUCKHTS_SQLLOGIC_KEEP_ARTIFACTS") == "1":
        return
    verbose = os.environ.get("DUCKHTS_SQLLOGIC_CLEANUP_VERBOSE") == "1"
    for path in sorted(set(paths), key=str):
        try:
            if path.is_symlink() or path.is_file():
                path.unlink()
                if verbose:
                    print(f"[sqllogic cleanup] removed {path}")
            elif path.is_dir():
                shutil.rmtree(path)
                if verbose:
                    print(f"[sqllogic cleanup] removed {path}/")
        except FileNotFoundError:
            pass


class DuckHTSSQLLogicTestExecutor(sqllogic_runner.SQLLogicTestExecutor):
    def __init__(self, build_directory=None, default_test_dir=None):
        super().__init__(build_directory, default_test_dir)
        self._repo_root = Path.cwd().resolve()
        self._tracked = _tracked_files(self._repo_root)

    def execute_test(self, test):
        artifacts = _declared_artifacts(Path(test.path), self._repo_root, self._tracked)
        try:
            return super().execute_test(test)
        finally:
            _cleanup(artifacts)


def main() -> None:
    # Reuse the pinned upstream CLI and parser, replacing only its executor class.
    sqllogic_runner.SQLLogicTestExecutor = DuckHTSSQLLogicTestExecutor
    sqllogic_runner.SQLLogicPythonRunner().run()


if __name__ == "__main__":
    main()
