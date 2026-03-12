#!/usr/bin/env python3
"""Render package-bundled function catalog artifacts from functions.yaml.

The manifest is stored in JSON-formatted YAML so this script can use only the
Python standard library.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from collections import OrderedDict
from pathlib import Path


def die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def escape_md(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def load_manifest(path: Path) -> list[dict[str, object]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=OrderedDict)
    except json.JSONDecodeError as exc:
        die(f"Failed to parse {path}: {exc}")

    functions = payload.get("functions")
    if not isinstance(functions, list):
        die(f"{path} is missing a top-level 'functions' array")

    required = {"name", "kind", "category", "signature", "returns", "r_wrapper", "description", "examples"}
    seen = set()
    for index, entry in enumerate(functions):
        if not isinstance(entry, dict):
            die(f"functions[{index}] must be an object")
        missing = sorted(required - set(entry))
        if missing:
            die(f"functions[{index}] is missing required fields: {', '.join(missing)}")
        name = entry["name"]
        if not isinstance(name, str) or not name:
            die(f"functions[{index}].name must be a non-empty string")
        if name in seen:
            die(f"Duplicate function entry for {name}")
        seen.add(name)
        examples = entry["examples"]
        if not isinstance(examples, list) or not all(isinstance(x, str) for x in examples):
            die(f"functions[{index}].examples must be a list of strings")
    return functions


def render_markdown(functions: list[dict[str, object]]) -> str:
    by_category: OrderedDict[str, list[dict[str, object]]] = OrderedDict()
    for function in functions:
        by_category.setdefault(str(function["category"]), []).append(function)

    lines: list[str] = []
    lines.append("## Extension Function Catalog")
    lines.append("")
    lines.append("This section is generated from `functions.yaml`.")
    lines.append("")

    for category, entries in by_category.items():
        lines.append(f"### {category}")
        lines.append("")
        lines.append("| Function | Kind | Returns | R helper | Description |")
        lines.append("| --- | --- | --- | --- | --- |")
        for entry in entries:
            name = f"`{entry['name']}`"
            kind = escape_md(str(entry["kind"]))
            returns = escape_md(str(entry["returns"]))
            wrapper = str(entry["r_wrapper"]).strip()
            wrapper_md = f"`{wrapper}`" if wrapper else ""
            description = escape_md(str(entry["description"]))
            lines.append(f"| {name} | {kind} | {returns} | {wrapper_md} | {description} |")
        lines.append("")

    return "\n".join(lines)


def write_tsv(path: Path, functions: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["name", "kind", "category", "signature", "returns", "r_wrapper", "description", "examples"]
        )
        for entry in functions:
            writer.writerow(
                [
                    entry["name"],
                    entry["kind"],
                    entry["category"],
                    entry["signature"],
                    entry["returns"],
                    entry["r_wrapper"],
                    entry["description"],
                    " || ".join(entry["examples"]),
                ]
            )


def main(argv: list[str]) -> int:
    repo_root = Path(argv[1]).resolve() if len(argv) > 1 else Path(__file__).resolve().parents[1]
    manifest_path = repo_root / "functions.yaml"
    if not manifest_path.exists():
        die(f"Manifest not found: {manifest_path}")

    functions = load_manifest(manifest_path)
    catalog_dir = repo_root / "r" / "Rduckhts" / "inst" / "function_catalog"
    catalog_dir.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(manifest_path, catalog_dir / "functions.yaml")
    (catalog_dir / "functions.md").write_text(render_markdown(functions) + "\n", encoding="utf-8")
    write_tsv(catalog_dir / "functions.tsv", functions)
    print(f"Rendered {len(functions)} function entries into {catalog_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
