#!/usr/bin/env python3
"""Reject machine-local paths in benchmark chunks that render or execute."""

from __future__ import annotations

import re
import sys
from pathlib import Path

CHUNK_START = re.compile(r"^```\{([^}]*)}\s*$")
EVAL_FALSE = re.compile(r"\beval\s*=\s*FALSE\b", re.IGNORECASE)
EXECUTABLE_LANGUAGES = {"r", "bash", "sh", "python", "perl"}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    errors: list[str] = []

    for path in sorted((root / "benchmarks").glob("*.Rmd")):
        active = False
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            chunk = CHUNK_START.match(line)
            if chunk:
                header = chunk.group(1)
                engine = header.split(",", 1)[0].strip()
                language = engine.split(None, 1)[0].lower()
                active = language in EXECUTABLE_LANGUAGES and not EVAL_FALSE.search(header)
                continue
            if line.startswith("```"):
                active = False
                continue
            if active and "/root/" in line:
                errors.append(f"{path.relative_to(root)}:{line_number}: machine-local path in an active chunk")

    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("benchmark portability: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
