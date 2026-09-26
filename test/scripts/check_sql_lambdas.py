#!/usr/bin/env python3
"""Reject deprecated SQL lambdas in extension C string literals."""

from pathlib import Path
import re

C_STRING = re.compile(r'"(?:\\[^\n]|[^"\\\n])*"')
SQL_STRING = re.compile(r"'(?:''|[^'])*'")
OLD_LAMBDA = re.compile(
    r"(?<![\w.])(?:[A-Za-z_]\w*|\(\s*[A-Za-z_]\w*"
    r"(?:\s*,\s*[A-Za-z_]\w*)+\s*\))\s*->(?!\s*['\"\d])")


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
    assert not list(deprecated_lambdas('"j->\'key\'" "j->0" "\'a -> b\'" "lambda x: x + 1"'))

    root = Path(__file__).resolve().parents[2] / "src"
    findings = [
        f"{path}:{line}: {expression}"
        for path in sorted(root.rglob("*.c"))
        for line, expression in deprecated_lambdas(path.read_text())
    ]
    if findings:
        raise SystemExit("Deprecated SQL lambda in C literal:\n" + "\n".join(findings))
    print("SQL lambda syntax guard: passed")


if __name__ == "__main__":
    main()
