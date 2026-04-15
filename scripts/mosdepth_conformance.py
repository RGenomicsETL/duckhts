#!/usr/bin/env python3
"""
Native mosdepth conformance test for DuckHTS.

Validates that the native `duckhts_mosdepth(...)` table function reproduces the
currently implemented mosdepth-compatible fast-mode output set:

  1. per-base depth BED.gz
  2. summary.txt
  3. global.dist.txt
  4. window/BED region means
  5. region.dist.txt
  6. optional thresholds.bed.gz

All DuckHTS runs go through the native table function. There is no SQL
reconstruction fallback in this script.

Usage:
    python3 scripts/mosdepth_conformance.py <alignment_path> [--chrom 11] [--fasta ref.fa] [--extension path]

Requires: mosdepth, duckdb CLI, samtools
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile


def run(cmd, *, env=None, timeout=600, check=True):
    """Run a shell command; return stdout."""
    merged = {**os.environ, **(env or {})}
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, env=merged, timeout=timeout
    )
    if check and result.returncode != 0:
        print(f"FAIL: {cmd}\nstderr: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def sql_quote_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def duckdb_run_native(
    *,
    ext,
    prefix,
    bam,
    chrom=None,
    by=None,
    quantize=None,
    fasta=None,
    thresholds=None,
    threads=4,
    precision_digits=6,
    tmpdir,
):
    sql = [f"LOAD {sql_quote_string(ext)};"]
    args = [
        sql_quote_string(prefix),
        sql_quote_string(bam),
    ]
    if chrom is not None:
        args.append(f"chrom := {sql_quote_string(chrom)}")
    if by is not None:
        args.append(f"by := {sql_quote_string(by)}")
    if fasta is not None:
        args.append(f"fasta := {sql_quote_string(fasta)}")
    if quantize is not None:
        args.append(f"quantize := {sql_quote_string(quantize)}")
    if thresholds is not None:
        args.append(f"thresholds := {sql_quote_string(thresholds)}")
    args.append(f"threads := {threads}")
    args.append(f"precision_digits := {precision_digits}")
    args.append("overwrite := TRUE")
    sql.append(f"SELECT * FROM duckhts_mosdepth({', '.join(args)});")

    sqlfile = os.path.join(tmpdir, "query.sql")
    with open(sqlfile, "w", encoding="utf-8") as handle:
        handle.write("\n".join(sql) + "\n")
    run(f"duckdb -unsigned < {sqlfile}", timeout=1200)


def get_chrom_length(bam, chrom):
    """Get chromosome length from BAM/CRAM index."""
    for line in run(f"samtools idxstats {bam}").splitlines():
        parts = line.split("\t")
        if parts[0] == chrom:
            return int(parts[1])
    sys.exit(f"Chromosome {chrom} not found in {bam}")


def read_bed(path, chrom=None):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if chrom and parts[0] != chrom:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return rows


def read_dist(path, chrom):
    dist = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            contig, depth, proportion = line.rstrip("\n").split("\t")
            if contig == chrom:
                dist[int(depth)] = proportion
    return dist


def read_summary_row(path, chrom):
    with open(path, encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == chrom:
                return {
                    "length": parts[1],
                    "bases": parts[2],
                    "mean": parts[3],
                    "min": parts[4],
                    "max": parts[5],
                }
    raise RuntimeError(f"Chromosome {chrom} not found in summary {path}")


def compare_perbase(native_bed_gz, mosdepth_bed_gz, chrom):
    """
    Compare non-zero per-base BED intervals.

    The current native fast-mode writer can still emit zero-depth runs more
    explicitly than upstream mosdepth. Summary/distribution comparisons cover
    the zero-depth accounting; this comparison checks the non-zero intervals.
    """
    native_rows = [row for row in read_bed(native_bed_gz, chrom) if row[3] != "0"]
    mosdepth_rows = [row for row in read_bed(mosdepth_bed_gz, chrom) if row[3] != "0"]
    if native_rows == mosdepth_rows:
        return True, len(native_rows), 0
    diffs = sum(1 for a, b in zip(native_rows, mosdepth_rows) if a != b)
    diffs += abs(len(native_rows) - len(mosdepth_rows))
    return False, len(native_rows), diffs


def compare_distribution(native_dist, mosdepth_dist, chrom):
    native = read_dist(native_dist, chrom)
    mosdepth = read_dist(mosdepth_dist, chrom)
    mismatches = 0
    for depth in sorted(mosdepth):
        if native.get(depth) != mosdepth[depth]:
            mismatches += 1
    return mismatches == 0, len(mosdepth), mismatches


def compare_windows(native_bed_gz, mosdepth_bed_gz, chrom):
    native_rows = [
        row for row in read_bed(native_bed_gz, chrom) if row[3] != "0.000000"
    ]
    mosdepth_rows = [
        row for row in read_bed(mosdepth_bed_gz, chrom) if row[3] != "0.000000"
    ]
    if native_rows == mosdepth_rows:
        return True, len(native_rows), 0
    diffs = sum(1 for a, b in zip(native_rows, mosdepth_rows) if a != b)
    diffs += abs(len(native_rows) - len(mosdepth_rows))
    return False, len(native_rows), diffs


def compare_tsv(native_path, mosdepth_path, chrom=None, skip_header=False):
    native_rows = read_tsv_rows(native_path, chrom, skip_header=skip_header)
    mosdepth_rows = read_tsv_rows(mosdepth_path, chrom, skip_header=skip_header)
    if native_rows == mosdepth_rows:
        return True, len(native_rows), 0
    diffs = sum(1 for a, b in zip(native_rows, mosdepth_rows) if a != b)
    diffs += abs(len(native_rows) - len(mosdepth_rows))
    return False, len(native_rows), diffs


def read_tsv_rows(path, chrom=None, skip_header=False):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt", encoding="utf-8") as handle:
        if skip_header:
            next(handle, None)
        for line in handle:
            parts = tuple(line.rstrip("\n").split("\t"))
            if chrom and parts and parts[0] != chrom:
                continue
            rows.append(parts)
    return rows


def compare_thresholds(native_bed_gz, mosdepth_bed_gz, chrom):
    native_rows = read_tsv_rows(native_bed_gz, chrom, skip_header=True)
    mosdepth_rows = read_tsv_rows(mosdepth_bed_gz, chrom, skip_header=True)
    if native_rows == mosdepth_rows:
        return True, len(native_rows), 0
    diffs = sum(1 for a, b in zip(native_rows, mosdepth_rows) if a != b)
    diffs += abs(len(native_rows) - len(mosdepth_rows))
    return False, len(native_rows), diffs


def main():
    parser = argparse.ArgumentParser(description="Native mosdepth conformance test")
    parser.add_argument("bam", help="Input BAM or CRAM file (indexed)")
    parser.add_argument("--chrom", default="11", help="Chromosome to test (default: 11)")
    parser.add_argument("--fasta", default=None, help="Reference FASTA for CRAM input when required")
    parser.add_argument("--extension", default=None, help="Path to duckhts.duckdb_extension")
    parser.add_argument(
        "--window-size",
        type=int,
        default=10000,
        help="Window size for region means (default: 10000)",
    )
    parser.add_argument(
        "--by-bed",
        default=None,
        help="BED file for mosdepth --by comparisons; overrides --window-size when set",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Decompression threads for both mosdepth and duckhts_mosdepth (default: 4)",
    )
    parser.add_argument(
        "--thresholds",
        default=None,
        help="Optional comma-separated threshold list for --by comparisons (e.g. 1,10,20)",
    )
    parser.add_argument(
        "--quantize",
        default=None,
        help="Optional mosdepth quantize spec (e.g. :1:4:)",
    )
    parser.add_argument(
        "--keep-tmp", action="store_true", help="Keep temporary files for inspection"
    )
    args = parser.parse_args()

    bam = os.path.abspath(args.bam)
    chrom = args.chrom
    fasta = os.path.abspath(args.fasta) if args.fasta else None
    if bam.endswith(".cram") and not fasta:
        sys.exit("--fasta is required for CRAM conformance runs")

    if args.extension:
        ext = os.path.abspath(args.extension)
    else:
        ext = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "build",
                "release",
                "duckhts.duckdb_extension",
            )
        )
    if not os.path.exists(ext):
        sys.exit(f"Extension not found: {ext}")

    chrom_len = get_chrom_length(bam, chrom)
    print(f"Alignment:  {bam}")
    print(f"Chromosome: {chrom} ({chrom_len:,} bp)")
    print(f"Extension:  {ext}")
    if args.by_bed:
        print(f"By BED:     {os.path.abspath(args.by_bed)}")
    else:
        print(f"Window:     {args.window_size:,} bp")
    print(f"Threads:    {args.threads}")
    if args.thresholds:
        print(f"Thresholds: {args.thresholds}")
    if args.quantize:
        print(f"Quantize:   {args.quantize}")
    if fasta:
        print(f"Reference:  {fasta}")
    print()

    tmpdir = tempfile.mkdtemp(prefix="mosdepth_conf_")
    results = {}

    print("--- Running mosdepth (ground truth) ---")
    md_env = {"MOSDEPTH_PRECISION": "6"}
    fasta_flag = f" -f {fasta}" if fasta else ""
    by_value = os.path.abspath(args.by_bed) if args.by_bed else str(args.window_size)
    quantize_flag = f" --quantize {args.quantize}" if args.quantize else ""
    md_fast = os.path.join(tmpdir, "md_fast")
    md_windows = os.path.join(tmpdir, "md_fast_windows")
    run(
        f"mosdepth --fast-mode{quantize_flag} -t {args.threads}{fasta_flag} {md_fast} {bam}",
        env=md_env,
        timeout=3600,
    )
    run(
        f"mosdepth --fast-mode --by {by_value}{quantize_flag}"
        f"{(' -T ' + args.thresholds) if args.thresholds else ''}"
        f" -t {args.threads}{fasta_flag} {md_windows} {bam}",
        env=md_env,
        timeout=3600,
    )
    print("  mosdepth runs complete.")

    print("\n--- Running native duckhts_mosdepth(...) ---")
    native_prefix = os.path.join(tmpdir, "native_fast")
    duckdb_run_native(
        ext=ext,
        prefix=native_prefix,
        bam=bam,
        chrom=chrom,
        by=by_value,
        quantize=args.quantize,
        fasta=fasta,
        thresholds=args.thresholds,
        threads=args.threads,
        precision_digits=6,
        tmpdir=tmpdir,
    )
    print("  native run complete.")

    native_perbase = f"{native_prefix}.per-base.bed.gz"
    native_summary = f"{native_prefix}.mosdepth.summary.txt"
    native_global = f"{native_prefix}.mosdepth.global.dist.txt"
    native_regions = f"{native_prefix}.regions.bed.gz"
    native_region_dist = f"{native_prefix}.mosdepth.region.dist.txt"
    native_quantized = f"{native_prefix}.quantized.bed.gz"
    native_thresholds = f"{native_prefix}.thresholds.bed.gz"

    print("\n--- Test 1: Fast-mode per-base depth ---")
    ok, n_rows, diffs = compare_perbase(native_perbase, f"{md_fast}.per-base.bed.gz", chrom)
    results["fast_perbase"] = ok
    print(f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} non-zero intervals, {diffs} differences")

    print("\n--- Test 2: Fast-mode summary ---")
    native_summary_row = read_summary_row(native_summary, chrom)
    mosdepth_summary_row = read_summary_row(f"{md_fast}.mosdepth.summary.txt", chrom)
    ok_bases = native_summary_row["bases"] == mosdepth_summary_row["bases"]
    ok_mean = native_summary_row["mean"] == mosdepth_summary_row["mean"]
    ok_max = native_summary_row["max"] == mosdepth_summary_row["max"]
    ok = ok_bases and ok_mean and ok_max
    results["fast_summary"] = ok
    print(
        f"  bases: native={native_summary_row['bases']} mosdepth={mosdepth_summary_row['bases']} "
        f"{'PASS' if ok_bases else 'FAIL'}"
    )
    print(
        f"  mean:  native={native_summary_row['mean']} mosdepth={mosdepth_summary_row['mean']} "
        f"{'PASS' if ok_mean else 'FAIL'}"
    )
    print(
        f"  max:   native={native_summary_row['max']} mosdepth={mosdepth_summary_row['max']} "
        f"{'PASS' if ok_max else 'FAIL'}"
    )

    print("\n--- Test 3: Fast-mode global distribution ---")
    ok, n_levels, mismatches = compare_distribution(
        native_global, f"{md_fast}.mosdepth.global.dist.txt", chrom
    )
    results["fast_distribution"] = ok
    print(f"  {'PASS' if ok else 'FAIL'}: {n_levels} depth levels, {mismatches} mismatches")

    if args.by_bed:
        print("\n--- Test 4: Fast-mode BED region means ---")
        ok, n_rows, diffs = compare_tsv(native_regions, f"{md_windows}.regions.bed.gz", chrom)
        region_label = "rows"
    else:
        print(f"\n--- Test 4: Fast-mode {args.window_size}-bp region means ---")
        ok, n_rows, diffs = compare_windows(
            native_regions, f"{md_windows}.regions.bed.gz", chrom
        )
        region_label = "non-zero windows"
    results["fast_windows"] = ok
    print(f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} {region_label}, {diffs} differences")

    print("\n--- Test 5: Fast-mode region distribution ---")
    ok, n_levels, mismatches = compare_distribution(
        native_region_dist, f"{md_windows}.mosdepth.region.dist.txt", chrom
    )
    results["fast_region_distribution"] = ok
    print(f"  {'PASS' if ok else 'FAIL'}: {n_levels} depth levels, {mismatches} mismatches")

    if args.thresholds:
        print("\n--- Test 6: Fast-mode thresholds ---")
        ok, n_rows, diffs = compare_thresholds(
            native_thresholds, f"{md_windows}.thresholds.bed.gz", chrom
        )
        results["fast_thresholds"] = ok
        print(f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} rows, {diffs} differences")

    if args.quantize:
        print("\n--- Test 7: Fast-mode quantized output ---")
        ok, n_rows, diffs = compare_tsv(
            native_quantized, f"{md_fast}.quantized.bed.gz", chrom
        )
        results["fast_quantize"] = ok
        print(f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} rows, {diffs} differences")

    print("\n" + "=" * 60)
    passed = sum(1 for ok in results.values() if ok)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")
    for name, ok in results.items():
        print(f"  {'PASS' if ok else 'FAIL'}: {name}")

    if args.keep_tmp:
        print(f"\nTemp files kept at: {tmpdir}")
    else:
        shutil.rmtree(tmpdir)
        print("\nTemp files cleaned up.")

    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
