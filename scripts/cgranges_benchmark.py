#!/usr/bin/env python3
"""Compatibility launcher for the R-owned cgranges benchmark.

Use ``r/duckhtsbench/scripts/cgranges_benchmark.R`` directly for new automation.
This launcher preserves the historic Python path while forwarding all options to
an optparse-based R entry point.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    script = root / "r" / "duckhtsbench" / "scripts" / "cgranges_benchmark.R"
    os.execvp("Rscript", ["Rscript", str(script), *sys.argv[1:]])


if __name__ == "__main__":
    main()
