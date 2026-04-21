#!/usr/bin/env python3
"""
Benchmark DuckHTS cgranges overlap probes against bedtk.

This benchmark is honest about the current public API shape:
- DuckHTS builds a session-scoped immutable cgranges index once
- queries are issued as repeated bound overlap probes
- bedtk runs whole-query-file existence filtering with `bedtk flt`

The script reports end-to-end timing and peak RSS for both tools, while also
splitting DuckHTS build/index time from repeated query-pass time.
"""

import argparse
import csv
import json
import os
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def run_timed(cmd, *, timeout=3600, check=True):
    start = time.monotonic()
    time_cmd = f"/usr/bin/time -v sh -c {subprocess.list2cmdline([cmd])} "
    result = subprocess.run(
        time_cmd,
        shell=True,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - start
    peak_rss_mb = 0.0
    stderr_lines = result.stderr.splitlines()
    actual_stderr = []
    for line in stderr_lines:
        if "Maximum resident set size" in line:
            try:
                peak_rss_mb = int(line.split(":", 1)[1].strip()) / 1024.0
            except (ValueError, IndexError):
                pass
        else:
            actual_stderr.append(line)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"command failed: {cmd}\nstdout:\n{result.stdout}\nstderr:\n" + "\n".join(actual_stderr)
        )
    return result.stdout, elapsed, peak_rss_mb


def generate_rows(subjects, queries, seed):
    rng = random.Random(seed)
    chroms = ["chr1", "chr2", "chr3", "chr4"]
    chrom_span = 5_000_000

    subject_rows = []
    for i in range(subjects):
        chrom = chroms[i % len(chroms)]
        width = 50 + (i % 200)
        start = rng.randrange(0, chrom_span - width)
        end = start + width
        subject_rows.append((chrom, start, end, f"s{i}"))

    query_rows = []
    for i in range(queries):
        chrom = chroms[(i * 3) % len(chroms)]
        width = 30 + (i % 120)
        start = rng.randrange(0, chrom_span - width)
        end = start + width
        query_rows.append((chrom, start, end, f"q{i}"))

    return subject_rows, query_rows


def write_bed(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def _worker_duckhts(args):
    import duckdb

    con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    try:
        con.execute(f"LOAD '{str(Path(args.extension).resolve()).replace("'", "''")}'")
        build_start = time.monotonic()
        subject_sql = (
            "CREATE TABLE cgr_subject AS "
            "SELECT column0 AS chrom, CAST(column1 AS BIGINT) AS start, "
            "CAST(column2 AS BIGINT) AS \"end\", column3 AS label "
            f"FROM read_csv_auto('{str(Path(args.subject_bed).resolve()).replace("'", "''")}', "
            "delim := '\\t', header := false)"
        )
        con.execute(subject_sql)
        con.execute(
            "SELECT duckhts_cgranges_from_query('bench_idx', "
            "'SELECT chrom, start, \"end\", label FROM cgr_subject', "
            "'chrom', 'start', 'end', 'label')"
        )
        con.execute("SELECT duckhts_cgranges_index('bench_idx')")
        build_elapsed = time.monotonic() - build_start

        query_rows = []
        with open(args.query_bed, "r", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for chrom, start, end, name in reader:
                query_rows.append((chrom, int(start), int(end), name))

        pass_seconds = []
        matched = None
        total_hits = None
        for _ in range(args.passes):
            t0 = time.monotonic()
            pass_matched = 0
            pass_hits = 0
            for chrom, start, end, _name in query_rows:
                rows = con.execute(
                    "SELECT interval_ordinal FROM duckhts_cgranges_overlaps(?, ?, ?, ?)",
                    ["bench_idx", chrom, start, end],
                ).fetchall()
                if rows:
                    pass_matched += 1
                    pass_hits += len(rows)
            pass_seconds.append(time.monotonic() - t0)
            if matched is None:
                matched = pass_matched
                total_hits = pass_hits
            elif matched != pass_matched or total_hits != pass_hits:
                raise RuntimeError("DuckHTS repeated passes produced inconsistent results")

        con.execute("SELECT duckhts_cgranges_destroy('bench_idx')")

        payload = {
            "tool": "duckhts_cgranges",
            "subjects": args.subjects,
            "queries": args.queries,
            "passes": args.passes,
            "build_index_sec": build_elapsed,
            "query_total_sec": sum(pass_seconds),
            "query_pass_1_sec": pass_seconds[0] if pass_seconds else 0.0,
            "matched_query_intervals": matched or 0,
            "total_hits": total_hits or 0,
        }
        Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    finally:
        con.close()


def _worker_bedtk(args):
    bedtk = str(Path(args.bedtk).resolve())
    pass_seconds = []
    matched = None
    for _ in range(args.passes):
        t0 = time.monotonic()
        proc = subprocess.run(
            [bedtk, "flt", args.subject_bed, args.query_bed],
            check=True,
            capture_output=True,
            text=True,
        )
        pass_seconds.append(time.monotonic() - t0)
        pass_matched = sum(1 for line in proc.stdout.splitlines() if line.strip())
        if matched is None:
            matched = pass_matched
        elif matched != pass_matched:
            raise RuntimeError("bedtk repeated passes produced inconsistent results")

    payload = {
        "tool": "bedtk",
        "subjects": args.subjects,
        "queries": args.queries,
        "passes": args.passes,
        "build_index_sec": 0.0,
        "query_total_sec": sum(pass_seconds),
        "query_pass_1_sec": pass_seconds[0] if pass_seconds else 0.0,
        "matched_query_intervals": matched or 0,
        "total_hits": None,
    }
    Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def bench_main(args):
    repo_root = Path(__file__).resolve().parents[1]
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    subject_rows, query_rows = generate_rows(args.subjects, args.queries, args.seed)
    subject_bed = out_dir / "subject.bed"
    query_bed = out_dir / "query.bed"
    write_bed(subject_bed, subject_rows)
    write_bed(query_bed, query_rows)

    script_path = Path(__file__).resolve()
    duck_json = out_dir / "duckhts.json"
    bedtk_json = out_dir / "bedtk.json"

    print("DuckHTS cgranges vs bedtk benchmark")
    print(f"Subjects:   {args.subjects:,}")
    print(f"Queries:    {args.queries:,}")
    print(f"Passes:     {args.passes}")
    print(f"Seed:       {args.seed}")
    print(f"Extension:  {Path(args.extension).resolve()}")
    print(f"bedtk:      {Path(args.bedtk).resolve()}")

    duck_cmd = subprocess.list2cmdline([
        sys.executable,
        str(script_path),
        "--worker-duckhts",
        "--extension",
        str(Path(args.extension).resolve()),
        "--subject-bed",
        str(subject_bed),
        "--query-bed",
        str(query_bed),
        "--passes",
        str(args.passes),
        "--subjects",
        str(args.subjects),
        "--queries",
        str(args.queries),
        "--json-out",
        str(duck_json),
    ])
    print("\n--- Running DuckHTS cgranges benchmark ---")
    _, duck_elapsed, duck_rss = run_timed(duck_cmd, timeout=args.timeout)
    duck_payload = json.loads(duck_json.read_text(encoding="utf-8"))
    duck_payload["total_elapsed_sec"] = duck_elapsed
    duck_payload["peak_rss_mb"] = duck_rss

    bedtk_cmd = subprocess.list2cmdline([
        sys.executable,
        str(script_path),
        "--worker-bedtk",
        "--bedtk",
        str(Path(args.bedtk).resolve()),
        "--subject-bed",
        str(subject_bed),
        "--query-bed",
        str(query_bed),
        "--passes",
        str(args.passes),
        "--subjects",
        str(args.subjects),
        "--queries",
        str(args.queries),
        "--json-out",
        str(bedtk_json),
    ])
    print("\n--- Running bedtk benchmark ---")
    _, bedtk_elapsed, bedtk_rss = run_timed(bedtk_cmd, timeout=args.timeout)
    bedtk_payload = json.loads(bedtk_json.read_text(encoding="utf-8"))
    bedtk_payload["total_elapsed_sec"] = bedtk_elapsed
    bedtk_payload["peak_rss_mb"] = bedtk_rss

    if duck_payload["matched_query_intervals"] != bedtk_payload["matched_query_intervals"]:
        raise RuntimeError(
            "DuckHTS/bedtk mismatch on matched query intervals: "
            f"{duck_payload['matched_query_intervals']} vs {bedtk_payload['matched_query_intervals']}"
        )

    summary_rows = [duck_payload, bedtk_payload]
    csv_path = out_dir / "summary.csv"
    fieldnames = [
        "tool",
        "subjects",
        "queries",
        "passes",
        "build_index_sec",
        "query_total_sec",
        "query_pass_1_sec",
        "total_elapsed_sec",
        "peak_rss_mb",
        "matched_query_intervals",
        "total_hits",
    ]
    with open(csv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print("\nSummary")
    for row in summary_rows:
        print(
            f"- {row['tool']}: total={row['total_elapsed_sec']:.2f}s, "
            f"peak_rss={row['peak_rss_mb']:.1f} MB, matched_queries={row['matched_query_intervals']}"
        )
    print(f"\nWrote {csv_path}")


def build_parser():
    parser = argparse.ArgumentParser(description="Benchmark DuckHTS cgranges against bedtk")
    parser.add_argument("--extension", default="build/release/duckhts.duckdb_extension")
    parser.add_argument("--bedtk", default=".sync/bedtk/bedtk")
    parser.add_argument("--subjects", type=int, default=50000)
    parser.add_argument("--queries", type=int, default=5000)
    parser.add_argument("--passes", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out-dir", default=".tmp/cgranges_benchmark")
    parser.add_argument("--timeout", type=int, default=3600)

    parser.add_argument("--worker-duckhts", action="store_true")
    parser.add_argument("--worker-bedtk", action="store_true")
    parser.add_argument("--subject-bed")
    parser.add_argument("--query-bed")
    parser.add_argument("--json-out")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.worker_duckhts:
        _worker_duckhts(args)
        return
    if args.worker_bedtk:
        _worker_bedtk(args)
        return
    bench_main(args)


if __name__ == "__main__":
    main()
