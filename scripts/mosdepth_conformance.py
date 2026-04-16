#!/usr/bin/env python3
"""
Native mosdepth conformance test for DuckHTS.

Validates that the native `duckhts_mosdepth(...)` table function reproduces the
implemented mosdepth-compatible output set for the selected mode:

  1. per-base depth BED.gz
  2. summary.txt
  3. global.dist.txt
  4. window/BED region means
  5. region.dist.txt
  6. optional thresholds.bed.gz

All DuckHTS runs go through the native table function. There is no SQL
reconstruction fallback in this script.

Usage:
    python3 scripts/mosdepth_conformance.py <alignment_path> [--mode fast|default|fragment] [--chrom 11] [--fasta ref.fa] [--extension path]

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
    read_groups=None,
    flag=1796,
    include_flag=0,
    mapq=0,
    min_frag_len=None,
    max_frag_len=None,
    fast_mode=False,
    fragment_mode=False,
    use_median=False,
    threads=4,
    processing_threads=None,
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
    if read_groups is not None:
        args.append(f"read_groups := {sql_quote_string(read_groups)}")
    args.append(f"flag := {flag}")
    args.append(f"include_flag := {include_flag}")
    args.append(f"mapq := {mapq}")
    if quantize is not None:
        args.append(f"quantize := {sql_quote_string(quantize)}")
    if thresholds is not None:
        args.append(f"thresholds := {sql_quote_string(thresholds)}")
    if min_frag_len is not None:
        args.append(f"min_frag_len := {min_frag_len}")
    if max_frag_len is not None:
        args.append(f"max_frag_len := {max_frag_len}")
    args.append(f"fast_mode := {'TRUE' if fast_mode else 'FALSE'}")
    args.append(f"fragment_mode := {'TRUE' if fragment_mode else 'FALSE'}")
    args.append(f"use_median := {'TRUE' if use_median else 'FALSE'}")
    args.append(f"threads := {threads}")
    args.append(f"precision_digits := {precision_digits}")
    if processing_threads is not None and processing_threads > 0:
        args.append(f"processing_threads := {processing_threads}")
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
    """Compare per-base BED intervals (full comparison including zero-depth)."""
    native_rows = read_bed(native_bed_gz, chrom)
    mosdepth_rows = read_bed(mosdepth_bed_gz, chrom)
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
    """Compare windowed region BED intervals (full comparison including zero-depth)."""
    native_rows = read_bed(native_bed_gz, chrom)
    mosdepth_rows = read_bed(mosdepth_bed_gz, chrom)
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
    parser.add_argument(
        "--mode",
        choices=("fast", "default", "fragment"),
        default="fast",
        help="Mosdepth mode to validate (default: fast)",
    )
    parser.add_argument(
        "--chrom", default="11", help="Chromosome to test (default: 11)"
    )
    parser.add_argument(
        "--fasta", default=None, help="Reference FASTA for CRAM input when required"
    )
    parser.add_argument(
        "--extension", default=None, help="Path to duckhts.duckdb_extension"
    )
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
        "--processing-threads",
        type=int,
        default=None,
        help="Number of parallel contig processing threads for duckhts_mosdepth (default: 0 = sequential)",
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
        "--read-groups",
        default=None,
        help="Optional comma-separated RG filter",
    )
    parser.add_argument(
        "--flag",
        type=int,
        default=1796,
        help="Exclude reads with any of these SAM flag bits set (default: 1796)",
    )
    parser.add_argument(
        "--include-flag",
        type=int,
        default=0,
        help="Only include reads with any of these SAM flag bits set (default: 0)",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=0,
        help="Ignore reads with MAPQ less than this threshold (default: 0)",
    )
    parser.add_argument(
        "--min-frag-len",
        type=int,
        default=None,
        help="Optional minimum absolute template length filter",
    )
    parser.add_argument(
        "--max-frag-len",
        type=int,
        default=None,
        help="Optional maximum absolute template length filter",
    )
    parser.add_argument(
        "--use-median",
        action="store_true",
        help="Validate mosdepth -m / use_median region output",
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
    print(f"Mode:       {args.mode}")
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
    if args.read_groups:
        print(f"ReadGroups: {args.read_groups}")
    if args.flag != 1796:
        print(f"Flag:       {args.flag}")
    if args.include_flag != 0:
        print(f"IncludeFlg: {args.include_flag}")
    if args.mapq != 0:
        print(f"MAPQ:       {args.mapq}")
    if args.min_frag_len is not None:
        print(f"MinFragLen: {args.min_frag_len}")
    if args.max_frag_len is not None:
        print(f"MaxFragLen: {args.max_frag_len}")
    if args.use_median:
        print("UseMedian:  TRUE")
    if fasta:
        print(f"Reference:  {fasta}")
    print()

    tmpdir = tempfile.mkdtemp(prefix="mosdepth_conf_")
    results = {}

    print("--- Running mosdepth (ground truth) ---")
    md_env = {"MOSDEPTH_PRECISION": "6"}
    fasta_flag = f" -f {fasta}" if fasta else ""
    chrom_flag = f" --chrom {chrom}" if chrom else ""
    by_value = os.path.abspath(args.by_bed) if args.by_bed else str(args.window_size)
    quantize_flag = f" --quantize {args.quantize}" if args.quantize else ""
    if args.mode == "fast":
        mode_flag = " --fast-mode"
    elif args.mode == "fragment":
        mode_flag = " --fragment-mode"
    else:
        mode_flag = ""
    read_groups_flag = f" -R {args.read_groups}" if args.read_groups else ""
    flag_opt = f" -F {args.flag}" if args.flag != 1796 else ""
    include_flag_opt = f" -i {args.include_flag}" if args.include_flag != 0 else ""
    mapq_opt = f" -Q {args.mapq}" if args.mapq != 0 else ""
    min_frag_flag = f" -l {args.min_frag_len}" if args.min_frag_len is not None else ""
    max_frag_flag = f" -u {args.max_frag_len}" if args.max_frag_len is not None else ""
    median_flag = " -m" if args.use_median else ""
    md_main = os.path.join(tmpdir, f"md_{args.mode}")
    md_windows = os.path.join(tmpdir, f"md_{args.mode}_windows")
    run(
        f"mosdepth{mode_flag}{quantize_flag}{read_groups_flag}{flag_opt}{include_flag_opt}{mapq_opt}{min_frag_flag}{max_frag_flag}{median_flag}"
        f" -t {args.threads}{fasta_flag}{chrom_flag} {md_main} {bam}",
        env=md_env,
        timeout=3600,
    )
    run(
        f"mosdepth{mode_flag}{median_flag} --by {by_value}{quantize_flag}"
        f"{(' -T ' + args.thresholds) if args.thresholds else ''}{read_groups_flag}{flag_opt}{include_flag_opt}{mapq_opt}"
        f"{min_frag_flag}{max_frag_flag} -t {args.threads}{fasta_flag}{chrom_flag} {md_windows} {bam}",
        env=md_env,
        timeout=3600,
    )
    print("  mosdepth runs complete.")

    print("\n--- Running native duckhts_mosdepth(...) ---")
    native_prefix = os.path.join(tmpdir, f"native_{args.mode}")
    duckdb_run_native(
        ext=ext,
        prefix=native_prefix,
        bam=bam,
        chrom=chrom,
        by=by_value,
        quantize=args.quantize,
        fasta=fasta,
        thresholds=args.thresholds,
        read_groups=args.read_groups,
        flag=args.flag,
        include_flag=args.include_flag,
        mapq=args.mapq,
        min_frag_len=args.min_frag_len,
        max_frag_len=args.max_frag_len,
        fast_mode=args.mode == "fast",
        fragment_mode=args.mode == "fragment",
        use_median=args.use_median,
        threads=args.threads,
        processing_threads=args.processing_threads,
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

    if args.mode == "fast":
        mode_label = "Fast-mode"
    elif args.mode == "fragment":
        mode_label = "Fragment-mode"
    else:
        mode_label = "Default-mode"

    print(f"\n--- Test 1: {mode_label} per-base depth ---")
    ok, n_rows, diffs = compare_perbase(
        native_perbase, f"{md_main}.per-base.bed.gz", chrom
    )
    results[f"{args.mode}_perbase"] = ok
    print(f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} intervals, {diffs} differences")

    print(f"\n--- Test 2: {mode_label} summary ---")
    native_summary_row = read_summary_row(native_summary, chrom)
    mosdepth_summary_row = read_summary_row(f"{md_main}.mosdepth.summary.txt", chrom)
    ok_bases = native_summary_row["bases"] == mosdepth_summary_row["bases"]
    ok_mean = native_summary_row["mean"] == mosdepth_summary_row["mean"]
    ok_max = native_summary_row["max"] == mosdepth_summary_row["max"]
    ok = ok_bases and ok_mean and ok_max
    results[f"{args.mode}_summary"] = ok
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

    print(f"\n--- Test 3: {mode_label} global distribution ---")
    ok, n_levels, mismatches = compare_distribution(
        native_global, f"{md_main}.mosdepth.global.dist.txt", chrom
    )
    results[f"{args.mode}_distribution"] = ok
    print(
        f"  {'PASS' if ok else 'FAIL'}: {n_levels} depth levels, {mismatches} mismatches"
    )

    if args.by_bed:
        print(f"\n--- Test 4: {mode_label} BED region means ---")
        ok, n_rows, diffs = compare_tsv(
            native_regions, f"{md_windows}.regions.bed.gz", chrom
        )
        region_label = "rows"
    else:
        print(f"\n--- Test 4: {mode_label} {args.window_size}-bp region means ---")
        ok, n_rows, diffs = compare_windows(
            native_regions, f"{md_windows}.regions.bed.gz", chrom
        )
        region_label = "windows"
    results[f"{args.mode}_windows"] = ok
    print(
        f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} {region_label}, {diffs} differences"
    )

    print(f"\n--- Test 5: {mode_label} region distribution ---")
    ok, n_levels, mismatches = compare_distribution(
        native_region_dist, f"{md_windows}.mosdepth.region.dist.txt", chrom
    )
    results[f"{args.mode}_region_distribution"] = ok
    print(
        f"  {'PASS' if ok else 'FAIL'}: {n_levels} depth levels, {mismatches} mismatches"
    )

    if args.thresholds:
        print(f"\n--- Test 6: {mode_label} thresholds ---")
        ok, n_rows, diffs = compare_thresholds(
            native_thresholds, f"{md_windows}.thresholds.bed.gz", chrom
        )
        results[f"{args.mode}_thresholds"] = ok
        print(f"  {'PASS' if ok else 'FAIL'}: {n_rows:,} rows, {diffs} differences")

    if args.quantize:
        print(f"\n--- Test 7: {mode_label} quantized output ---")
        ok, n_rows, diffs = compare_tsv(
            native_quantized, f"{md_main}.quantized.bed.gz", chrom
        )
        results[f"{args.mode}_quantize"] = ok
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
