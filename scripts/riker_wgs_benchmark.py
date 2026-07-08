#!/usr/bin/env python3
"""Harmonized WGS coverage benchmark: riker wgs vs mosdepth vs duckhts_mosdepth.

Rendered report: benchmarks/benchmark_riker_wgs.Rmd. Stage the input BAM with
scripts/stage_riker_wgs_bam.sh first.

Mirrors the methodology of fulcrumgenomics/riker's benchmark-pipeline
(config/mosdepth-compare.config.yaml, workflow/rules/run_{riker,mosdepth}.smk)
at riker v0.4.0:

  * same sample (HG00188_30x, ERR3240174, transcoded CRAM -> BAM)
  * cold page cache before every run (requires root to drop caches)
  * harmonized read selection: MAPQ 20, orphans counted, dup+secondary+
    supplementary excluded, qcfail counted, mosdepth in accurate mode
  * replicates for min/max bars

Deliberate differences from upstream's harness:

  * every tool is hard-capped to the same set of *distinct physical cores*
    with taskset (--cpus). On a hybrid CPU (Intel P/E-cores, or SMT siblings)
    an uncapped thread budget would let a tool with more internal threads
    (duckhts_mosdepth spawns processing_threads workers, each with its own
    htslib decode) silently outspend the others. Pass one core id per thread
    of the largest budget; pick distinct physical cores, not SMT siblings.
  * riker is additionally run with --min-bq 0. Upstream states BQ filtering is
    "the one difference we couldn't harmonize away" because *mosdepth* has no
    such switch -- but riker does. riker-bq0 is the apples-to-apples
    computation; riker (default) reproduces their published configuration.

Example:

    sudo python3 scripts/riker_wgs_benchmark.py \\
        --bam riker-bench-data/stage/HG00188_30x/input.bam \\
        --ref riker-bench-data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \\
        --ext build/release/duckhts.duckdb_extension \\
        --riker-bin riker \\
        --outdir riker-bench-data/results --tsv riker-bench-data/results/bench.tsv \\
        --tools mosdepth,duckhts,riker,riker-bq0 --threads 1,2,4 --reps 3
"""

import argparse
import os
import re
import shutil
import subprocess
import time


def cpulist(cpus, nthreads):
    if nthreads > len(cpus):
        raise SystemExit(
            f"thread budget {nthreads} exceeds {len(cpus)} pinned cores; "
            f"pass a longer --cpus list"
        )
    return ",".join(str(c) for c in cpus[:nthreads])


def drop_caches():
    subprocess.run(["sync"], check=True)
    with open("/proc/sys/vm/drop_caches", "w") as fh:
        fh.write("3\n")


def parse_gnu_time(path):
    """Return (wall_seconds, max_rss_mb, pct_cpu) from /usr/bin/time -v output."""
    wall = rss = pct = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Elapsed (wall clock) time"):
                v = line.split(": ", 1)[1]
                wall = 0.0
                for p in v.split(":"):
                    wall = wall * 60 + float(p)
            elif line.startswith("Maximum resident set size"):
                rss = int(line.split(":")[-1].strip()) / 1024.0
            elif line.startswith("Percent of CPU this job got"):
                m = re.match(r".*: (\d+)%", line)
                if m:
                    pct = int(m.group(1))
    return wall, rss, pct


def run_one(cmd, rundir, cpus, timeout):
    os.makedirs(rundir, exist_ok=True)
    timefile = os.path.join(rundir, "time.txt")
    logfile = os.path.join(rundir, "tool.log")
    full = ["taskset", "-c", cpus, "/usr/bin/time", "-v", "-o", timefile] + cmd
    with open(os.path.join(rundir, "cmdline.txt"), "w") as fh:
        fh.write(" ".join(full) + "\n")
    drop_caches()
    t0 = time.monotonic()
    with open(logfile, "wb") as lf:
        proc = subprocess.run(full, stdout=lf, stderr=subprocess.STDOUT, timeout=timeout)
    elapsed = time.monotonic() - t0
    if proc.returncode != 0:
        tail = open(logfile).read()[-2000:]
        raise SystemExit(f"FAILED rc={proc.returncode}: {' '.join(full)}\n{tail}")
    wall, rss, pct = parse_gnu_time(timefile)
    return (wall if wall is not None else elapsed), rss, pct


def build_cmd(tool, nthreads, args, rundir):
    prefix = os.path.join(rundir, "out")
    if tool == "mosdepth":
        # -t is EXTRA decompression threads on top of the main thread.
        return ["mosdepth", "-t", str(max(nthreads - 1, 0)),
                "-Q", "20", "-F", "3332", "--no-per-base", prefix, args.bam]

    if tool == "duckhts":
        # processing_threads = contig workers; threads = htslib decode threads
        # per worker. N workers x inline decode == N threads total.
        sql = (
            "SET threads TO 1; "
            f"LOAD '{args.ext}'; "
            "SELECT count(*) FROM duckhts_mosdepth("
            f"'{prefix}', '{args.bam}', no_per_base := TRUE, mapq := 20, flag := 3332, "
            f"threads := 0, processing_threads := {nthreads}, overwrite := TRUE);"
        )
        return ["duckdb", "-unsigned", "-c", sql]

    if tool in ("riker", "riker-bq0"):
        cmd = [args.riker_bin, "wgs", "--input", args.bam, "--output", prefix,
               "--reference", args.ref, "--threads", str(nthreads),
               "--include-unpaired-reads", "--min-mapq", "20"]
        if tool == "riker-bq0":
            cmd += ["--min-bq", "0"]
        return cmd

    raise ValueError(tool)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True, help="reference FASTA (indexed) for riker + CRAM")
    ap.add_argument("--ext", required=True, help="built duckhts.duckdb_extension")
    ap.add_argument("--riker-bin", required=True, help="riker binary (name on PATH or path)")
    ap.add_argument("--outdir", required=True, help="per-run scratch (time.txt/tool.log)")
    ap.add_argument("--tsv", required=True, help="results TSV (appended)")
    ap.add_argument("--tools", default="mosdepth,duckhts,riker,riker-bq0")
    ap.add_argument("--threads", default="1,2,4", help="comma-separated thread budgets")
    ap.add_argument("--cpus", default="0,2,4,6,8,10",
                    help="distinct physical core ids to pin to (one per thread of "
                         "the largest budget; avoid SMT siblings / E-cores)")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--timeout", type=int, default=14400)
    args = ap.parse_args()

    if os.geteuid() != 0:
        raise SystemExit("must be root to drop the page cache before each run")

    tools = args.tools.split(",")
    budgets = [int(t) for t in args.threads.split(",")]
    cpus = [int(c) for c in args.cpus.split(",")]

    os.makedirs(os.path.dirname(os.path.abspath(args.tsv)), exist_ok=True)
    new = not os.path.exists(args.tsv)
    tsv = open(args.tsv, "a", buffering=1)
    if new:
        tsv.write("tool\tthreads\trep\twall_s\tmax_rss_mb\tpct_cpu\n")

    for nthreads in budgets:
        pin = cpulist(cpus, nthreads)
        for tool in tools:
            for rep in range(1, args.reps + 1):
                rundir = os.path.join(args.outdir, f"t{nthreads}", tool, f"rep{rep}")
                if os.path.exists(rundir):
                    shutil.rmtree(rundir)
                cmd = build_cmd(tool, nthreads, args, rundir)
                wall, rss, pct = run_one(cmd, rundir, pin, args.timeout)
                tsv.write(f"{tool}\t{nthreads}\t{rep}\t{wall:.2f}\t{rss:.1f}\t{pct}\n")
                print(f"[t{nthreads}] {tool:12s} rep{rep}  {wall:8.2f}s  "
                      f"{rss / 1024:6.2f} GB  {pct}% cpu", flush=True)

    tsv.close()


if __name__ == "__main__":
    main()
