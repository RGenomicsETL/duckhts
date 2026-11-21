#!/usr/bin/env python3
"""Reject htslib read-open calls outside the file-access module.

Every reader opens its inputs and their indexes through
src/hts_file_access.c, which is what lets `SET duckhts_file_system = 'duckdb'`
apply to all of them. A source file that calls one of the htslib functions
below reads with htslib's own file access whatever the setting says.

ALLOWED names the call sites that do so on purpose. An entry that no longer
matches a call fails too, so the list cannot outlive its sites.
"""

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]

# The module, and the htslib scheme handler it opens routed names through.
MODULE = {"src/hts_file_access.c", "src/duckdb_fs_hfile.c"}

READ_OPEN_FUNCTIONS = (
    "hts_open", "hts_open_format", "hts_hopen", "sam_open", "sam_open_format",
    "bcf_open", "vcf_open", "cram_open", "bgzf_open", "bgzf_hopen", "hopen",
    "fai_load", "fai_load3", "fai_load_format", "fai_load3_format",
    "sam_index_load", "sam_index_load2", "sam_index_load3",
    "bcf_index_load", "bcf_index_load2", "bcf_index_load3",
    "tbx_index_load", "tbx_index_load2", "tbx_index_load3",
    "hts_idx_load", "hts_idx_load2", "hts_idx_load3",
    "hts_set_fai_filename", "bcf_sr_add_reader", "bwOpen", "bbOpen",
)
# hts_set_opt(fp, CRAM_OPT_REFERENCE, path) opens the reference by name.
READ_OPEN_OPTIONS = ("CRAM_OPT_REFERENCE",)

ALLOWED = {
    ("src/bgzip.c", "bgzf_open"):
        "bgunzip writes its output with this process; bgzip reads its input with stdio",
    ("src/reference_cache.c", "fai_load3_format"):
        "FAI_CREATE builds a missing index, and the cached handle outlives the connection",
    ("src/libbigwig_hfile_io.c", "bwOpen"):
        "the libBigWig adapter: its urlOpen() opens the file through the module",
}

COMMENT_OR_LITERAL = re.compile(
    r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'', re.S)
CALL = re.compile(r"(?<![\w.>])(" + "|".join(READ_OPEN_FUNCTIONS) + r")\s*\(")
OPTION = re.compile(r"\b(" + "|".join(READ_OPEN_OPTIONS) + r")\b")


def blank(match):
    """Comments and literals become spaces; newlines stay, so line numbers hold."""
    return re.sub(r"[^\n]", " ", match.group())


def read_opens(source):
    """Yield (line number, name) for each read-open call or option in C source."""
    code = COMMENT_OR_LITERAL.sub(blank, source)
    for pattern in (CALL, OPTION):
        for match in pattern.finditer(code):
            yield code.count("\n", 0, match.start()) + 1, match.group(1)


def main():
    sources = sorted(
        path for pattern in ("src/**/*.c", "src/**/*.h") for path in ROOT.glob(pattern))
    if len(sources) < 100:
        sys.exit("file access guard: found %d sources under src/; wrong directory?" % len(sources))
    violations = []
    used = set()
    for path in sources:
        name = path.relative_to(ROOT).as_posix()
        if name in MODULE:
            continue
        for line, function in read_opens(path.read_text(encoding="utf-8", errors="replace")):
            if (name, function) in ALLOWED:
                used.add((name, function))
            else:
                violations.append("%s:%d: %s" % (name, line, function))
    stale = sorted(set(ALLOWED) - used)
    for entry in violations:
        print(entry)
    for name, function in stale:
        print("%s: allowed call to %s is no longer there; remove the entry" % (name, function))
    if violations or stale:
        sys.exit(
            "file access guard: %d call(s) outside src/hts_file_access.c, %d stale allowance(s).\n"
            "Open inputs with the duckhts_file_access_* functions (src/include/hts_file_access.h)."
            % (len(violations), len(stale)))
    print("file access guard: %d sources, %d allowed htslib-only call sites: passed"
          % (len(sources), len(ALLOWED)))


if __name__ == "__main__":
    main()
