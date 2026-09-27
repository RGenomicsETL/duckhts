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


# Scripts and SQL tests embed SQL in every string form their language has (single,
# double and triple quotes, R raw strings, JavaScript templates, sqllogictest
# records). Rather than parse each form, scan all text outside comments and flag
# only an arrow in argument position (right after "(" or ","), which is where a
# lambda parameter always sits; prose such as "history -> {path}" never does.
SCRIPT_SUFFIXES = {".R", ".Rmd", ".py", ".sql", ".sh", ".js", ".mjs", ".html", ".test"}
HASH_COMMENTS = {".R", ".Rmd", ".py", ".sh", ".test"}
ARGUMENT_LAMBDA = re.compile(
    r"[(,]\s*((?:[A-Za-z_]\w*|\(\s*[A-Za-z_]\w*(?:\s*,\s*[A-Za-z_]\w*)*\s*\))\s*->)")
EXCLUDED = ("third_party/", "r/Rduckhts/inst/duckhts_extension/", "test/scripts/check_sql_lambdas.py")


def strip_comments(source, suffix):
    lines = []
    for line in source.split("\n"):
        if suffix in HASH_COMMENTS:
            line = "" if line.lstrip().startswith("#") else re.sub(r"\s#\s.*$", "", line)
        if suffix in (".sql", ".test"):
            line = re.sub(r"--.*$", "", line)
        if suffix in (".js", ".mjs", ".html"):
            line = re.sub(r"(^|\s)//.*$", "", line)
        lines.append(line)
    return "\n".join(lines)


def script_lambdas(source, suffix):
    text = strip_comments(source, suffix)
    for match in ARGUMENT_LAMBDA.finditer(text):
        yield text.count("\n", 0, match.start()) + 1, match.group(1)


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
    assert [e for _, e in script_lambdas("q <- 'SELECT list_transform(xs, x -> x)'", ".R")] == ["x ->"]
    assert [e for _, e in script_lambdas("q = '''SELECT list_filter(xs, (x) -> x > 1)'''", ".py")] == ["(x) ->"]
    assert [e for _, e in script_lambdas("q = `SELECT list_reduce(xs, (a, b) -> a + b)`;", ".mjs")] == ["(a, b) ->"]
    assert [e for _, e in script_lambdas('```{sql}\nSELECT list_filter(xs, x -> x > 1);\n```\n', ".Rmd")] == ["x ->"]
    assert [e for _, e in script_lambdas("SELECT list_filter(xs, x -> x > 1);", ".sql")] == ["x ->"]
    assert [e for _, e in script_lambdas("query I\nSELECT list_transform([1], x -> x);\n----\n[1]\n", ".test")] == ["x ->"]
    assert not list(script_lambdas('cat(glue("history -> {path}"))\n# empty, a -> b\nx <- 1  # A,M -> NULL', ".R"))
    assert not list(script_lambdas("-- list_transform(xs, x -> x)\n# (x) -> y", ".test"))

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
