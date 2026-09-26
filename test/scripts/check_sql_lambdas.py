#!/usr/bin/env python3
"""Reject deprecated SQL lambdas in extension C string literals."""

from pathlib import Path
import re

C_STRING = re.compile(r'"(?:\\[^\n]|[^"\\\n])*"')
SQL_STRING = re.compile(r"'(?:''|[^'])*'")
# Any `name ->` or `(a, b) ->` outside SQL string literals is flagged, whatever the body
# (`x -> 1` and `x -> 'c'` are deprecated lambdas too). Extension SQL does not use the JSON
# arrow operators; write json_extract(...) instead so this guard stays exact.
OLD_LAMBDA = re.compile(
    r"(?<![\w.])(?:[A-Za-z_]\w*|\(\s*[A-Za-z_]\w*"
    r"(?:\s*,\s*[A-Za-z_]\w*)+\s*\))\s*->")


def deprecated_lambdas(source):
    for literal in C_STRING.finditer(source):
        sql = SQL_STRING.sub("''", literal.group())
        for match in OLD_LAMBDA.finditer(sql):
            yield source.count("\n", 0, literal.start()) + 1, match.group()


def main():
    assert list(deprecated_lambdas('x->member; "list_transform(xs, x -> x + 1)"')) == [
        (1, "x ->")]
    assert list(deprecated_lambdas('"list_reduce(xs, (a, b) -> a + b)"')) == [
        (1, "(a, b) ->")]
    assert [e for _, e in deprecated_lambdas('"list_transform(xs, x -> 1)" "list_transform(xs, x -> \'c\')"')] == [
        "x ->", "x ->"]
    assert not list(deprecated_lambdas('"\'a -> b\'" "lambda x: x + 1"'))

    root = Path(__file__).resolve().parents[2] / "src"
    findings = [
        f"{path}:{line}: {expression}"
        for path in sorted(root.rglob("*.c"))
        for line, expression in deprecated_lambdas(path.read_text(encoding="utf-8"))
    ]
    if findings:
        raise SystemExit("Deprecated SQL lambda in C literal:\n" + "\n".join(findings))
    print("SQL lambda syntax guard: passed")


if __name__ == "__main__":
    main()
