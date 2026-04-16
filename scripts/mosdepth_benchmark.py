#!/usr/bin/env python3
"""
Mosdepth vs native DuckHTS mosdepth benchmark.

Times upstream `mosdepth` against the native `duckhts_mosdepth(...)` table
function. The benchmark intentionally exercises only the native implementation;
there is no SQL reconstruction path in this script.

Usage:
    python3 scripts/mosdepth_benchmark.py <alignment_path> [options]

Examples:
    python3 scripts/mosdepth_benchmark.py NA12878.bam
    python3 scripts/mosdepth_benchmark.py my.bam --chrom 11 --verify
    python3 scripts/mosdepth_benchmark.py my.cram --fasta ref.fa --window-size 1000 --mode default
    python3 scripts/mosdepth_benchmark.py my.bam --mode fragment --chrom 11

Requires: mosdepth, duckdb CLI, samtools
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import time


def run(cmd, *, env=None, timeout=3600, check=True):
    """Run a shell command; return (stdout, elapsed_seconds, peak_rss_mb)."""
    merged = {**os.environ, **(env or {})}
    start = time.monotonic()
    # Use GNU time to capture per-process peak RSS accurately
    time_cmd = f"/usr/bin/time -v sh -c {subprocess.list2cmdline([cmd])} "
    result = subprocess.run(
        time_cmd,
        shell=True,
        capture_output=True,
        text=True,
        env=merged,
        timeout=timeout,
    )
    elapsed = time.monotonic() - start
    # Parse peak RSS from GNU time stderr output
    peak_rss_mb = 0.0
    stderr_lines = result.stderr.split("\n")
    actual_stderr = []
    for line in stderr_lines:
        if "Maximum resident set size" in line:
            # Format: "Maximum resident set size (kbytes): 12345"
            try:
                peak_rss_mb = int(line.strip().split(":")[-1].strip()) / 1024.0
            except (ValueError, IndexError):
                pass
        else:
            actual_stderr.append(line)
    if check and result.returncode != 0:
        print(f"FAIL: {cmd}\nstderr: {''.join(actual_stderr)}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip(), elapsed, peak_rss_mb


def sql_quote_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def duckdb_run_native(
    *,
    ext,
    prefix,
    bam,
    chrom=None,
    fasta=None,
    by=None,
    no_per_base=False,
    flag=1796,
    include_flag=0,
    fast_mode=False,
    fragment_mode=False,
    mapq=0,
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
    if no_per_base:
        args.append("no_per_base := TRUE")
    args.append(f"flag := {flag}")
    args.append(f"include_flag := {include_flag}")
    args.append(f"fast_mode := {'TRUE' if fast_mode else 'FALSE'}")
    args.append(f"fragment_mode := {'TRUE' if fragment_mode else 'FALSE'}")
    args.append(f"mapq := {mapq}")
    args.append(f"use_median := {'TRUE' if use_median else 'FALSE'}")
    args.append(f"threads := {threads}")
    if processing_threads is not None and processing_threads > 0:
        args.append(f"processing_threads := {processing_threads}")
    args.append(f"precision_digits := {precision_digits}")
    args.append("overwrite := TRUE")
    sql.append(f"SELECT * FROM duckhts_mosdepth({', '.join(args)});")

    sqlfile = os.path.join(tmpdir, "query.sql")
    with open(sqlfile, "w", encoding="utf-8") as handle:
        handle.write("\n".join(sql) + "\n")
    _, elapsed, peak_rss = run(f"duckdb -unsigned < {sqlfile}", timeout=3600)
    return elapsed, peak_rss


def get_chrom_lengths(bam):
    out, _, _ = run(f"samtools idxstats {bam}")
    lengths = {}
    for line in out.splitlines():
        parts = line.split("\t")
        if parts[0] != "*" and int(parts[1]) > 0:
            lengths[parts[0]] = int(parts[1])
    return lengths


def get_total_reads(bam):
    out, _, _ = run(f"samtools idxstats {bam}")
    mapped = 0
    unmapped = 0
    for line in out.splitlines():
        parts = line.split("\t")
        if parts[0] == "*":
            unmapped += int(parts[3])
            continue
        mapped += int(parts[2])
        unmapped += int(parts[3])
    return mapped + unmapped


def fmt_time(seconds):
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes = int(seconds) // 60
    remainder = seconds - (minutes * 60)
    return f"{minutes}m{remainder:.1f}s"


def fmt_ratio(mosdepth_s, native_s):
    if native_s == 0:
        return "inf"
    ratio = mosdepth_s / native_s
    if ratio >= 1:
        return f"{ratio:.2f}x faster"
    return f"{1 / ratio:.2f}x slower"


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


def compare_bed(native_bed_gz, mosdepth_bed_gz, chrom):
    native_rows = read_bed(native_bed_gz, chrom)
    mosdepth_rows = read_bed(mosdepth_bed_gz, chrom)
    if native_rows == mosdepth_rows:
        return True, len(native_rows), 0
    diffs = sum(1 for a, b in zip(native_rows, mosdepth_rows) if a != b)
    diffs += abs(len(native_rows) - len(mosdepth_rows))
    return False, len(native_rows), diffs


def main():
    parser = argparse.ArgumentParser(description="Mosdepth vs native DuckHTS benchmark")
    parser.add_argument("bam", help="Input BAM or CRAM file (indexed)")
    parser.add_argument(
        "--mode",
        choices=("fast", "default", "fragment"),
        default="fast",
        help="Mosdepth mode to benchmark (default: fast)",
    )
    parser.add_argument(
        "--fasta", default=None, help="Reference FASTA for CRAM input when required"
    )
    parser.add_argument(
        "--chrom", default=None, help="Single chromosome to test (default: whole file)"
    )
    parser.add_argument(
        "--extension", default=None, help="Path to duckhts.duckdb_extension"
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
        "--runs",
        type=int,
        default=1,
        help="Number of timing runs for each tool (default: 1)",
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
        "--window-size",
        type=int,
        default=None,
        help="If set, benchmark `--by <window-size>` region output instead of per-base output",
    )
    parser.add_argument(
        "--by-bed",
        default=None,
        help="BED file for mosdepth --by comparisons; overrides --window-size when set",
    )
    parser.add_argument(
        "--no-per-base",
        action="store_true",
        help="Pass no_per_base := TRUE to duckhts_mosdepth and -n to mosdepth",
    )
    parser.add_argument(
        "--use-median",
        action="store_true",
        help="Pass use_median := TRUE / -m when benchmarking region output",
    )
    parser.add_argument(
        "--verify", action="store_true", help="Compare first-run outputs"
    )
    parser.add_argument("--keep-tmp", action="store_true", help="Keep temporary files")
    args = parser.parse_args()

    bam = os.path.abspath(args.bam)
    fasta = os.path.abspath(args.fasta) if args.fasta else None
    by_bed = os.path.abspath(args.by_bed) if args.by_bed else None
    if bam.endswith(".cram") and not fasta:
        sys.exit("--fasta is required for CRAM benchmarks")

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

    chrom_lengths = get_chrom_lengths(bam)
    if args.chrom and args.chrom not in chrom_lengths:
        sys.exit(f"Chromosome {args.chrom} not found in {bam}")
    total_reads = get_total_reads(bam)
    scope_bp = chrom_lengths[args.chrom] if args.chrom else sum(chrom_lengths.values())

    print("=" * 70)
    print("Mosdepth vs native DuckHTS benchmark")
    print("=" * 70)
    print(f"Alignment:    {os.path.basename(bam)}")
    print(f"Mode:         {args.mode}")
    print(f"Reads:        {total_reads:,}")
    print(f"Scope:        {args.chrom or 'whole file'} ({scope_bp:,} bp)")
    print(f"Threads:      {args.threads}")
    print(f"Runs:         {args.runs}")
    print(f"Extension:    {ext}")
    if fasta:
        print(f"Reference:    {fasta}")
    if by_bed:
        print(f"By BED:       {by_bed}")
    elif args.window_size is not None:
        print(f"Window size:  {args.window_size:,}")
    if args.no_per_base:
        print("Per-base:     disabled")
    if args.flag != 1796:
        print(f"Flag:         {args.flag}")
    if args.include_flag != 0:
        print(f"Include flag: {args.include_flag}")
    if args.mapq != 0:
        print(f"MAPQ:         {args.mapq}")
    if args.use_median:
        print("Use median:   TRUE")
    print()

    tmpdir = tempfile.mkdtemp(prefix="mosdepth_bench_")
    md_env = {"MOSDEPTH_PRECISION": "6"}
    fasta_flag = f" -f {fasta}" if fasta else ""
    chrom_flag = f" --chrom {args.chrom}" if args.chrom else ""
    # Determine the --by value for both mosdepth and duckhts_mosdepth
    if by_bed:
        by_value = by_bed
    elif args.window_size is not None:
        by_value = str(args.window_size)
    else:
        by_value = None
    window_flag = f" --by {by_value}" if by_value is not None else ""
    no_per_base_flag = " -n" if args.no_per_base else ""
    flag_opt = f" -F {args.flag}" if args.flag != 1796 else ""
    include_flag_opt = f" -i {args.include_flag}" if args.include_flag != 0 else ""
    mapq_opt = f" -Q {args.mapq}" if args.mapq != 0 else ""
    if args.mode == "fast":
        mode_flag = " --fast-mode"
    elif args.mode == "fragment":
        mode_flag = " --fragment-mode"
    else:
        mode_flag = ""
    median_flag = " -m" if args.use_median else ""

    mosdepth_times = []
    native_times = []
    mosdepth_rss_values = []
    native_rss_values = []
    md_prefix_first = None
    native_prefix_first = None

    for run_idx in range(args.runs):
        label = f"run {run_idx + 1}/{args.runs}" if args.runs > 1 else "run"
        md_prefix = os.path.join(tmpdir, f"md_r{run_idx}")
        native_prefix = os.path.join(tmpdir, f"native_r{run_idx}")

        print(f"mosdepth {label}...")
        _, md_time, md_rss = run(
            f"mosdepth{mode_flag}{median_flag}{flag_opt}{include_flag_opt}{mapq_opt} -t {args.threads}{fasta_flag}{chrom_flag}{window_flag}{no_per_base_flag} {md_prefix} {bam}",
            env=md_env,
            timeout=3600,
        )
        mosdepth_times.append(md_time)
        mosdepth_rss_values.append(md_rss)
        print(f"  {fmt_time(md_time)}  (peak RSS: {md_rss:.0f} MB)")

        print(f"duckhts_mosdepth {label}...")
        native_time, native_rss = duckdb_run_native(
            ext=ext,
            prefix=native_prefix,
            bam=bam,
            chrom=args.chrom,
            fasta=fasta,
            by=by_value,
            no_per_base=args.no_per_base,
            flag=args.flag,
            include_flag=args.include_flag,
            fast_mode=args.mode == "fast",
            fragment_mode=args.mode == "fragment",
            mapq=args.mapq,
            use_median=args.use_median,
            threads=args.threads,
            processing_threads=args.processing_threads,
            precision_digits=6,
            tmpdir=tmpdir,
        )
        native_times.append(native_time)
        native_rss_values.append(native_rss)
        print(f"  {fmt_time(native_time)}  (peak RSS: {native_rss:.0f} MB)")

        if md_prefix_first is None:
            md_prefix_first = md_prefix
        if native_prefix_first is None:
            native_prefix_first = native_prefix

    mosdepth_avg = sum(mosdepth_times) / len(mosdepth_times)
    native_avg = sum(native_times) / len(native_times)

    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    print(f"{'Tool':<20} {'Average':>12} {'Best':>12} {'Peak RSS':>12}")
    print("-" * 70)
    print(
        f"{'mosdepth':<20} {fmt_time(mosdepth_avg):>12} {fmt_time(min(mosdepth_times)):>12} {max(mosdepth_rss_values):>9.0f} MB"
    )
    print(
        f"{'duckhts_mosdepth':<20} {fmt_time(native_avg):>12} {fmt_time(min(native_times)):>12} {max(native_rss_values):>9.0f} MB"
    )
    print("-" * 70)
    print(f"Ratio: DuckHTS is {fmt_ratio(mosdepth_avg, native_avg)} vs mosdepth")

    if args.verify and md_prefix_first and native_prefix_first:
        print("\nVerification:")
        verify_chrom = args.chrom or next(iter(chrom_lengths))
        if by_value is None:
            ok, n_rows, diffs = compare_bed(
                f"{native_prefix_first}.per-base.bed.gz",
                f"{md_prefix_first}.per-base.bed.gz",
                verify_chrom,
            )
            print(
                f"  per-base ({verify_chrom}): {'PASS' if ok else 'FAIL'} "
                f"({n_rows:,} intervals, {diffs} diffs)"
            )
        else:
            ok, n_rows, diffs = compare_bed(
                f"{native_prefix_first}.regions.bed.gz",
                f"{md_prefix_first}.regions.bed.gz",
                verify_chrom,
            )
            by_label = os.path.basename(by_bed) if by_bed else f"{args.window_size}bp"
            print(
                f"  regions [{by_label}] ({verify_chrom}): {'PASS' if ok else 'FAIL'} "
                f"({n_rows:,} rows, {diffs} diffs)"
            )

    if args.keep_tmp:
        print(f"\nTemp dir: {tmpdir}")
    else:
        shutil.rmtree(tmpdir)
        print("\nTemp dir cleaned up.")


if __name__ == "__main__":
    main()
