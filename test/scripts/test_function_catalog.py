#!/usr/bin/env python3
"""Network-free checks for the function catalog's generated documentation."""

from __future__ import annotations

import contextlib
import copy
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import render_function_catalog as catalog  # noqa: E402


def manifest_fixture() -> dict[str, object]:
    return {
        "community_extension": {
            "extension": {
                "name": "example",
                "description": "Example extension",
                "language": "C",
                "build": "cmake",
                "license": "MIT",
                "maintainers": ["example"],
            },
            "repo": {"github": "example/catalog", "ref": "test-revision"},
            "docs": {
                "hello_world_lines": ["SELECT sample_scan();"],
                "extended_intro": ["An example catalog."],
                "feature_notes": ["Explicit per-call settings."],
            },
        },
        "functions": [
            {
                "name": "sample_scan",
                "kind": "table",
                "category": "Readers",
                "signature": "sample_scan(path := 'a|b', marker := '```')",
                "returns": "STRUCT(value VARCHAR, markers VARCHAR /* ```` */)",
                "r_wrapper": "r_sample_scan",
                "description": "Scan records | preserve fields.",
                "details": {
                    "Input semantics": "Use `a|b` unchanged.\n\nPreserve source order.",
                    "Limits": "A full buffer returns a partial chunk; an oversized item errors.",
                },
                "examples": ["SELECT 'a|b', '```';\nSELECT 'second statement';"],
            },
            {
                "name": "sample_version",
                "kind": "scalar",
                "category": "Diagnostics",
                "signature": "sample_version()",
                "returns": "VARCHAR",
                "r_wrapper": "",
                "description": "Return the runtime version.",
                "examples": [],
            },
        ],
    }


class FunctionCatalogTests(unittest.TestCase):
    def load_fixture(self, manifest: dict[str, object]) -> dict[str, object]:
        with tempfile.TemporaryDirectory(prefix="duckhts-catalog-test-") as directory:
            path = Path(directory) / "functions.yaml"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            return catalog.load_manifest(path)

    def test_details_are_optional_and_keep_heading_order(self) -> None:
        manifest = self.load_fixture(manifest_fixture())
        self.assertEqual(list(manifest["functions"][0]["details"]), ["Input semantics", "Limits"])
        self.assertNotIn("details", manifest["functions"][1])

    def test_malformed_details_fail_before_generation(self) -> None:
        cases = [None, [], "prose", 1, {"": "prose"}, {" ": "prose"},
                 {"two\nlines": "prose"}, {"two\rlines": "prose"},
                 {"Heading": None}, {"Heading": []}, {"Heading": ""}, {"Heading": " \n"}]
        for details in cases:
            with self.subTest(details=details):
                manifest = manifest_fixture()
                manifest["functions"][0]["details"] = details
                error = io.StringIO()
                with contextlib.redirect_stderr(error), self.assertRaises(SystemExit) as raised:
                    self.load_fixture(manifest)
                self.assertEqual(raised.exception.code, 1)
                self.assertIn("functions[0].details", error.getvalue())

    def test_duplicate_headings_are_not_silently_dropped(self) -> None:
        with tempfile.TemporaryDirectory(prefix="duckhts-catalog-test-") as directory:
            path = Path(directory) / "functions.yaml"
            path.write_text('{"details": {"Input": "first", "Input": "second"}}', encoding="utf-8")
            error = io.StringIO()
            with contextlib.redirect_stderr(error), self.assertRaises(SystemExit):
                catalog.load_manifest(path)
            self.assertIn("Duplicate manifest property: Input", error.getvalue())

    def test_reference_retains_complete_fields_and_safe_code_fences(self) -> None:
        manifest = manifest_fixture()
        functions = manifest["functions"]
        before = copy.deepcopy(functions)
        reference = catalog.render_reference(functions)
        self.assertEqual(functions, before)
        self.assertEqual(reference.count("\n## "), 2)
        self.assertIn("\n## sample_scan\n", reference)
        self.assertIn("\n## sample_version\n", reference)
        for entry in functions:
            for field in ("signature", "returns", "description"):
                self.assertIn(entry[field], reference)
            for example in entry["examples"]:
                self.assertIn(example, reference)
            for heading, prose in entry.get("details", {}).items():
                self.assertIn(f"### {heading}\n\n{prose}\n", reference)
        self.assertLess(reference.index("### Input semantics"), reference.index("### Limits"))
        self.assertIn("````sql\n" + functions[0]["signature"] + "\n````\n", reference)
        self.assertIn("`````\n" + functions[0]["returns"] + "\n`````\n", reference)
        self.assertIn("````sql\n" + functions[0]["examples"][0] + "\n````\n", reference)

    def test_summary_links_omit_return_schemas_and_details(self) -> None:
        functions = manifest_fixture()["functions"]
        summary = catalog.render_markdown(functions)
        self.assertIn("| Function | Kind | R helper | Description |", summary)
        self.assertNotIn("| Returns |", summary)
        self.assertIn("[`sample_scan`](reference.md#sample_scan)", summary)
        self.assertIn("[`sample_version`](reference.md#sample_version)", summary)
        self.assertIn("Scan records \\| preserve fields.", summary)
        self.assertNotIn(functions[0]["returns"], summary)
        self.assertNotIn(functions[0]["details"]["Limits"], summary)

    def test_macro_identifiers_are_checked_for_c_generation(self) -> None:
        manifest = manifest_fixture()
        manifest["functions"][1]["name"] = 'unsafe"name'
        manifest["functions"][1]["kind"] = "scalar_macro"
        with tempfile.TemporaryDirectory(prefix="duckhts-catalog-test-") as directory:
            root = Path(directory)
            (root / "functions.yaml").write_text(json.dumps(manifest), encoding="utf-8")
            (root / "description.yml").write_text("version: 1.5.2\n", encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()) as error:
                with self.assertRaises(SystemExit):
                    catalog.main(["render_function_catalog.py", str(root)])
            self.assertIn("Invalid public macro identifier", error.getvalue())

    def test_checked_in_public_macro_header_matches_the_manifest(self) -> None:
        repo = Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory(prefix="duckhts-catalog-test-") as directory:
            root = Path(directory)
            # The header depends only on the function list; pin the descriptor's git
            # ref so rendering needs no repository.
            manifest = json.loads((repo / "functions.yaml").read_text(encoding="utf-8"))
            manifest["community_extension"]["repo"]["ref"] = "test-revision"
            (root / "functions.yaml").write_text(json.dumps(manifest), encoding="utf-8")
            (root / "description.yml").write_text(
                (repo / "description.yml").read_text(encoding="utf-8"), encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(catalog.main(["render_function_catalog.py", str(root)]), 0)
            self.assertEqual(
                (root / "src/include/duckhts_public_macros.h").read_text(encoding="utf-8"),
                (repo / "src/include/duckhts_public_macros.h").read_text(encoding="utf-8"))

    def test_main_writes_reference_and_keeps_the_tsv_schema(self) -> None:
        manifest = manifest_fixture()
        manifest["functions"][1]["kind"] = "scalar_macro"
        with tempfile.TemporaryDirectory(prefix="duckhts-catalog-test-") as directory:
            root = Path(directory)
            (root / "functions.yaml").write_text(json.dumps(manifest), encoding="utf-8")
            (root / "description.yml").write_text("version: 1.5.2\n", encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(catalog.main(["render_function_catalog.py", str(root)]), 0)
            generated = root / "r/Rduckhts/inst/function_catalog"
            macro_header = (root / "src/include/duckhts_public_macros.h").read_text()
            self.assertIn('    "sample_version",', macro_header)
            self.assertEqual(json.loads((generated / "functions.yaml").read_text()), manifest)
            self.assertEqual((generated / "reference.md").read_text(),
                             catalog.render_reference(manifest["functions"]) + "\n")
            self.assertFalse((generated / "reference.md").read_text().endswith("\n\n"))
            self.assertEqual((generated / "functions.md").read_text(),
                             catalog.render_markdown(manifest["functions"]) + "\n")
            with (generated / "functions.tsv").open(newline="") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            fields = ["name", "kind", "category", "signature", "returns", "r_wrapper",
                      "description", "examples"]
            self.assertEqual(rows[0], fields)
            for row, entry in zip(rows[1:], manifest["functions"], strict=True):
                expected = [entry[field] for field in fields[:-1]] + [" || ".join(entry["examples"])]
                self.assertEqual(row, expected)
            descriptor = (root / "community-extensions/extensions/duckhts/description.yml").read_text()
            self.assertIn("https://github.com/example/catalog/blob/test-revision/"
                          "r/Rduckhts/inst/function_catalog/reference.md", descriptor)
            self.assertIn("`sample_scan`: Scan records | preserve fields.", descriptor)
            self.assertNotIn(manifest["functions"][0]["signature"], descriptor)
            self.assertNotIn(manifest["functions"][0]["details"]["Limits"], descriptor)


if __name__ == "__main__":
    unittest.main()
