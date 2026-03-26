#!/usr/bin/env python3
"""Render generated documentation artifacts from functions.yaml.

The manifest is stored in JSON-formatted YAML so this script can use only the
Python standard library.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path


def die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def escape_md(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def load_root_description(path: Path) -> dict[str, object]:
    if not path.exists():
        die(f"Extension metadata not found: {path}")

    metadata: dict[str, object] = {}
    current_list_key: str | None = None
    current_list: list[str] | None = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if raw_line.startswith((" ", "\t")):
            if current_list_key is None or current_list is None:
                die(f"Unsupported indentation in {path}: {raw_line!r}")
            item = stripped
            if not item.startswith("- "):
                die(f"Unsupported list item syntax in {path}: {raw_line!r}")
            current_list.append(item[2:].strip())
            continue

        current_list_key = None
        current_list = None
        key, sep, value = raw_line.partition(":")
        if not sep:
            die(f"Invalid metadata line in {path}: {raw_line!r}")

        key = key.strip()
        value = value.strip()
        if not key:
            die(f"Invalid metadata key in {path}: {raw_line!r}")

        if value:
            metadata[key] = value
            continue

        metadata[key] = []
        current_list_key = key
        current_list = metadata[key]

    version = metadata.get("version")
    if not isinstance(version, str) or not version:
        die(f"{path} is missing a non-empty 'version' field")

    return metadata


def load_manifest(path: Path) -> OrderedDict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=OrderedDict)
    except json.JSONDecodeError as exc:
        die(f"Failed to parse {path}: {exc}")

    community_extension = payload.get("community_extension")
    if not isinstance(community_extension, dict):
        die(f"{path} is missing a top-level 'community_extension' object")

    extension = community_extension.get("extension")
    if not isinstance(extension, dict):
        die(f"{path} is missing community_extension.extension")
    required_extension = {"name", "description", "language", "build", "license", "maintainers"}
    missing_extension = sorted(required_extension - set(extension))
    if missing_extension:
        die(
            "community_extension.extension is missing required fields: "
            + ", ".join(missing_extension)
        )
    if not isinstance(extension["maintainers"], list) or not all(
        isinstance(x, str) and x for x in extension["maintainers"]
    ):
        die("community_extension.extension.maintainers must be a list of non-empty strings")

    repo = community_extension.get("repo")
    if not isinstance(repo, dict):
        die(f"{path} is missing community_extension.repo")
    required_repo = {"github"}
    missing_repo = sorted(required_repo - set(repo))
    if missing_repo:
        die("community_extension.repo is missing required fields: " + ", ".join(missing_repo))

    docs = community_extension.get("docs")
    if not isinstance(docs, dict):
        die(f"{path} is missing community_extension.docs")
    for field in ("hello_world_lines", "extended_intro", "feature_notes"):
        value = docs.get(field)
        if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
            die(f"community_extension.docs.{field} must be a list of strings")

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
    return payload


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


def quote_yaml_scalar(value: str) -> str:
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def build_extended_description(
    functions: list[dict[str, object]], docs: dict[str, object]
) -> list[str]:
    by_category: OrderedDict[str, list[dict[str, object]]] = OrderedDict()
    for function in functions:
        by_category.setdefault(str(function["category"]), []).append(function)

    lines: list[str] = []
    for paragraph in docs["extended_intro"]:
        lines.extend(str(paragraph).splitlines())
        lines.append("")

    lines.append("Functions included in this extension:")
    lines.append("")
    for category, entries in by_category.items():
        lines.append(f"### {category}")
        for entry in entries:
            lines.append(f"- `{entry['signature']}`: {entry['description']}")
        lines.append("")

    lines.append("Operational notes:")
    for note in docs["feature_notes"]:
        lines.append(f"- {note}")

    return lines


def resolve_repo_ref(repo_root: Path, repo: dict[str, object]) -> str:
    explicit_ref = repo.get("ref")
    if isinstance(explicit_ref, str) and explicit_ref:
        return explicit_ref

    ref_source = repo.get("ref_source", "git_head")
    if ref_source != "git_head":
        die(f"Unsupported community_extension.repo.ref_source: {ref_source}")

    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            cwd=repo_root,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        die(f"Failed to resolve git HEAD for community descriptor: {exc}")

    return result.stdout.strip()


def render_description_yaml(
    repo_root: Path,
    manifest: OrderedDict[str, object],
    functions: list[dict[str, object]],
    root_description: dict[str, object],
) -> str:
    community_extension = manifest["community_extension"]
    extension = community_extension["extension"]
    repo = community_extension["repo"]
    docs = community_extension["docs"]
    repo_ref = resolve_repo_ref(repo_root, repo)
    extended_description_lines = build_extended_description(functions, docs)

    lines: list[str] = []
    lines.append("extension:")
    lines.append(f"  name: {quote_yaml_scalar(str(extension['name']))}")
    lines.append(f"  description: {quote_yaml_scalar(str(extension['description']))}")
    lines.append(f"  version: {quote_yaml_scalar(str(root_description['version']))}")
    for field in ("language", "build", "license"):
        lines.append(f"  {field}: {quote_yaml_scalar(str(extension[field]))}")
    optional_fields = ("requires_toolchains", "excluded_platforms")
    for field in optional_fields:
        value = extension.get(field)
        if isinstance(value, str) and value:
            lines.append(f"  {field}: {quote_yaml_scalar(value)}")
    lines.append("  maintainers:")
    for maintainer in extension["maintainers"]:
        lines.append(f"    - {quote_yaml_scalar(str(maintainer))}")
    lines.append("")

    lines.append("repo:")
    lines.append(f"  github: {quote_yaml_scalar(str(repo['github']))}")
    lines.append(f"  ref: {quote_yaml_scalar(repo_ref)}")
    lines.append("")

    lines.append("docs:")
    lines.append("  hello_world: |")
    for line in docs["hello_world_lines"]:
        lines.append(f"    {line}" if line else "    ")
    lines.append("  extended_description: |")
    for line in extended_description_lines:
        lines.append(f"    {line}" if line else "    ")
    lines.append("")
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    repo_root = Path(argv[1]).resolve() if len(argv) > 1 else Path(__file__).resolve().parents[1]
    manifest_path = repo_root / "functions.yaml"
    if not manifest_path.exists():
        die(f"Manifest not found: {manifest_path}")
    root_description_path = repo_root / "description.yml"

    manifest = load_manifest(manifest_path)
    root_description = load_root_description(root_description_path)
    functions = manifest["functions"]
    catalog_dir = repo_root / "r" / "Rduckhts" / "inst" / "function_catalog"
    catalog_dir.mkdir(parents=True, exist_ok=True)
    description_path = repo_root / "community-extensions" / "extensions" / "duckhts" / "description.yml"
    description_path.parent.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(manifest_path, catalog_dir / "functions.yaml")
    (catalog_dir / "functions.md").write_text(render_markdown(functions) + "\n", encoding="utf-8")
    write_tsv(catalog_dir / "functions.tsv", functions)
    description_path.write_text(
        render_description_yaml(repo_root, manifest, functions, root_description),
        encoding="utf-8",
    )
    print(f"Rendered {len(functions)} function entries into {catalog_dir}")
    print(f"Rendered community extension descriptor into {description_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
