#!/usr/bin/env python3
"""Reject htslib calls that open a file by name outside the file-access module.

Every reader opens its inputs and their indexes through
src/hts_file_access.c, which is what lets `SET duckhts_file_system = 'duckdb'`
apply to all of them. A source file that calls one of the htslib functions
below opens that file with htslib's own file access whatever the setting says.

A call that does so on purpose carries a marker comment on the line of the
call or on the line directly above it:

    /* duckhts-file-access: allow <reason> */

The guard fails on a listed call without a marker, and on a marker with no
listed call on its line or the next one, so a marker cannot outlive its call.
"""

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]

# The module, and the htslib scheme handler it opens routed names through.
MODULE = {"src/hts_file_access.c", "src/duckdb_fs_hfile.c"}

# htslib (and libBigWig) functions that open a file named by an argument.
# hts_readlist and hts_readlines open their argument when it names a file,
# regidx_init and bcf_sr_regions_init when given a file name, and
# bcf_sr_set_regions, bcf_sr_set_targets and bcf_sr_set_samples when told the
# argument is a file. bam_index_load and bam_index_build are sam.h macros over
# the sam_index functions; the guard sees the macro name in source, so it lists
# both. The list covers the functions found in the htslib and libBigWig headers
# that open a file by name; a new one needs adding here.
NAME_OPEN_FUNCTIONS = (
    "hts_open", "hts_open_format", "hts_hopen", "sam_open", "sam_open_format",
    "bcf_open", "vcf_open", "cram_open", "bgzf_open", "bgzf_hopen", "hopen",
    "fai_load", "fai_load3", "fai_load_format", "fai_load3_format",
    "fai_build", "fai_build3",
    "sam_index_load", "sam_index_load2", "sam_index_load3",
    "bcf_index_load", "bcf_index_load2", "bcf_index_load3",
    "tbx_index_load", "tbx_index_load2", "tbx_index_load3",
    "hts_idx_load", "hts_idx_load2", "hts_idx_load3",
    "sam_index_build", "sam_index_build2", "sam_index_build3", "sam_index_build4",
    "bcf_index_build", "bcf_index_build2", "bcf_index_build3",
    "tbx_index_build", "tbx_index_build2", "tbx_index_build3",
    "bgzf_index_load", "bgzf_index_dump",
    "hts_set_fai_filename", "bcf_sr_add_reader",
    "bcf_sr_set_regions", "bcf_sr_set_targets", "bcf_sr_set_samples",
    "hts_readlist", "hts_readlines", "regidx_init", "bcf_sr_regions_init",
    "bcf_hdr_set", "hts_file_type", "bgzf_is_bgzf", "knet_open",
    "bcf_sweep_init", "bam_index_load", "bam_index_build",
    "bwOpen", "bbOpen", "bwIsBigWig", "bbIsBigBed",
)
# hts_set_opt(fp, CRAM_OPT_REFERENCE, path) opens the reference by name.
NAME_OPEN_OPTIONS = ("CRAM_OPT_REFERENCE",)

MARKER = re.compile(r"/\*\s*duckhts-file-access:\s*allow\s+\S[^*]*\*/")
COMMENT_OR_LITERAL = re.compile(
    r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'', re.S)
CALL = re.compile(r"(?<![\w.>])(" + "|".join(NAME_OPEN_FUNCTIONS) + r")\s*\(")
OPTION = re.compile(r"\b(" + "|".join(NAME_OPEN_OPTIONS) + r")\b")


def blank(match):
    """Comments and literals become spaces; newlines stay, so line numbers hold."""
    return re.sub(r"[^\n]", " ", match.group())


def line_of(text, offset):
    return text.count("\n", 0, offset) + 1


def name_opens(source):
    """Yield (line number, name) for each listed call or option in C source."""
    code = COMMENT_OR_LITERAL.sub(blank, source)
    for pattern in (CALL, OPTION):
        for match in pattern.finditer(code):
            yield line_of(code, match.start()), match.group(1)


def markers(source):
    """Line numbers that carry an allow marker."""
    return {line_of(source, match.start()) for match in MARKER.finditer(source)}


def main():
    sources = sorted(
        path for pattern in ("src/**/*.c", "src/**/*.h") for path in ROOT.glob(pattern))
    if len(sources) < 100:
        sys.exit("file access guard: found %d sources under src/; wrong directory?" % len(sources))
    violations = []
    stale = []
    allowed = 0
    for path in sources:
        name = path.relative_to(ROOT).as_posix()
        if name in MODULE:
            continue
        source = path.read_text(encoding="utf-8", errors="replace")
        marked = markers(source)
        calls = {}
        for line, function in name_opens(source):
            calls.setdefault(line, []).append(function)
        for line, functions in sorted(calls.items()):
            if line in marked or line - 1 in marked:
                allowed += len(functions)
            else:
                violations.extend("%s:%d: %s" % (name, line, function) for function in functions)
        stale.extend("%s:%d" % (name, line) for line in sorted(marked)
                     if line not in calls and line + 1 not in calls)
    for entry in violations:
        print(entry)
    for entry in stale:
        print("%s: allow marker without a listed call on its line or the next; remove it" % entry)
    if violations or stale:
        sys.exit(
            "file access guard: %d call(s) outside src/hts_file_access.c without an allow marker, "
            "%d stale marker(s).\n"
            "Open inputs with the duckhts_file_access_* functions (src/include/hts_file_access.h)."
            % (len(violations), len(stale)))
    print("file access guard: %d sources, %d marked htslib-only call sites: passed"
          % (len(sources), allowed))


if __name__ == "__main__":
    main()
