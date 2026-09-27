#!/usr/bin/env python3
"""Reject deprecated SQL lambdas in extension C sources and repository scripts."""

from pathlib import Path
import re
import subprocess

C_STRING = re.compile(r'"(?:\\[^\n]|[^"\\\n])*"')
SQL_STRING = re.compile(r"'(?:''|[^'])*'")
# Any `name ->` or `(a, b) ->` outside SQL string literals is flagged, whatever the body
# (`x -> 1` and `x -> 'c'` are deprecated lambdas too). Extension SQL does not use the JSON
# arrow operators; write json_extract(...) instead so this guard stays exact.
OLD_LAMBDA = re.compile(
    r"(?<![\w.])(?:[A-Za-z_]\w*|\(\s*[A-Za-z_]\w*"
    r"(?:\s*,\s*[A-Za-z_]\w*)*\s*\))\s*->")


# Adjacent C string literals (separated only by whitespace or comments) are one
# string after translation, so a lambda may be split across them.
BETWEEN_PARTS = re.compile(r"(?:\s|/\*.*?\*/|//[^\n]*)*\Z", re.S)


def sql_strings(source):
    """Yield (start offset, text) for each run of adjacent C string literals."""
    start = end = None
    parts = []
    for literal in C_STRING.finditer(source):
        if parts and BETWEEN_PARTS.match(source[end:literal.start()]):
            parts.append(literal.group()[1:-1])
        else:
            if parts:
                yield start, "".join(parts)
            start, parts = literal.start(), [literal.group()[1:-1]]
        end = literal.end()
    if parts:
        yield start, "".join(parts)


def deprecated_lambdas(source):
    for start, text in sql_strings(source):
        sql = SQL_STRING.sub("''", text)
        for match in OLD_LAMBDA.finditer(sql):
            yield source.count("\n", 0, start) + 1, match.group()


# Scripts (R, Rmd, Python, shell, JavaScript, HTML, SQL) carry SQL in string
# literals, and their prose also uses arrows ("history -> {path}", comments). There
# only string literals are scanned, and only an arrow in argument position (after
# "(" or ",") counts, which is where a lambda parameter always sits.
SCRIPT_SUFFIXES = {".R", ".Rmd", ".py", ".sql", ".sh", ".js", ".mjs", ".html"}
SCRIPT_STRING = re.compile(r'"(?:\\.|[^"\\])*"', re.S)
SQL_CHUNK = re.compile(r"^```\{sql[^}]*\}\n(.*?)^```", re.S | re.M)
ARGUMENT_LAMBDA = re.compile(
    r"[(,]\s*((?:[A-Za-z_]\w*|\(\s*[A-Za-z_]\w*(?:\s*,\s*[A-Za-z_]\w*)*\s*\))\s*->)")
EXCLUDED = ("third_party/", "r/Rduckhts/inst/duckhts_extension/", "test/scripts/check_sql_lambdas.py")


def script_lambdas(source, suffix):
    if suffix == ".sql":
        texts = [(0, source)]
    else:
        texts = [(m.start(), m.group()) for m in SCRIPT_STRING.finditer(source)]
        if suffix == ".Rmd":
            texts += [(m.start(1), m.group(1)) for m in SQL_CHUNK.finditer(source)]
    for start, text in texts:
        for match in ARGUMENT_LAMBDA.finditer(SQL_STRING.sub("''", text)):
            yield source.count("\n", 0, start) + 1, match.group(1)


def main():
    assert list(deprecated_lambdas('x->member; "list_transform(xs, x -> x + 1)"')) == [
        (1, "x ->")]
    assert list(deprecated_lambdas('"list_reduce(xs, (a, b) -> a + b)"')) == [
        (1, "(a, b) ->")]
    assert [e for _, e in deprecated_lambdas('"list_transform(xs, x -> 1)" "list_transform(xs, x -> \'c\')"')] == [
        "x ->", "x ->"]
    assert [e for _, e in deprecated_lambdas('"list_transform(xs, (x) -> x)"')] == ["(x) ->"]
    assert [e for _, e in deprecated_lambdas('"list_transform(xs, x"\n  " -> x)"')] == ["x ->"]
    assert [e for _, e in deprecated_lambdas('"list_transform(xs, x" /* part */ " -> x)"')] == ["x ->"]
    assert not list(deprecated_lambdas('"list_transform(xs, x"; y->z; " -> ignored-by-c"'))
    assert not list(deprecated_lambdas('"\'a -> b\'" "lambda x: x + 1"'))

    assert [e for _, e in script_lambdas('q <- "SELECT list_transform(xs,x->x)"', ".R")] == ["x->"]
    assert [e for _, e in script_lambdas('"list_reduce(xs, (a, b) -> a + b)"', ".py")] == ["(a, b) ->"]
    assert not list(script_lambdas('cat(glue("history -> {path}")) # empty -> NULL', ".R"))
    assert [e for _, e in script_lambdas('```{sql}\nSELECT list_filter(xs, x -> x > 1);\n```\n', ".Rmd")] == ["x ->"]
    assert [e for _, e in script_lambdas("SELECT list_filter(xs, x -> x > 1);", ".sql")] == ["x ->"]

    repo = Path(__file__).resolve().parents[2]
    findings = [
        f"{path}:{line}: {expression}"
        for path in sorted((repo / "src").rglob("*.c"))
        for line, expression in deprecated_lambdas(path.read_text(encoding="utf-8"))
    ]
    tracked = subprocess.run(["git", "-C", str(repo), "ls-files"], check=True,
                             capture_output=True, text=True).stdout.splitlines()
    for name in tracked:
        path = repo / name
        if path.suffix not in SCRIPT_SUFFIXES or name.startswith(EXCLUDED) or not path.is_file():
            continue
        findings += [f"{name}:{line}: {expression}" for line, expression in
                     script_lambdas(path.read_text(encoding="utf-8", errors="replace"), path.suffix)]
    if findings:
        raise SystemExit("Deprecated SQL lambda (use lambda x: ...):\n" + "\n".join(findings))
    print("SQL lambda syntax guard: passed")


if __name__ == "__main__":
    main()
