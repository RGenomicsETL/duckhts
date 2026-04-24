#!/usr/bin/env python3
"""Benchmark DuckHTS cgranges scalar probes against bedtk and bedtools.

The benchmark exercises the streaming-provider path, not the older
``overlaps_bulk(query_string)`` or one-SQL-call-per-query shapes:

* DuckHTS builds a session-scoped cgranges index once from the subject BED.
* Query BED rows are streamed through ``read_bed(...)``.
* ``duckhts_cgranges_has_overlap(...)`` filters query rows in a vectorized
  scalar predicate.
* ``duckhts_cgranges_count_overlaps(...)`` annotates query rows with overlap
  counts.
* ``bedtk flt`` is used as the overlap-existence baseline.
* ``bedtools intersect -u`` and ``bedtools intersect -c`` are included when a
  bedtools executable is available.

The default mode creates deterministic synthetic BED files so the rendered
benchmark is reproducible and quick. Pass ``--subject-bed`` and ``--query-bed``
to benchmark real BED files instead.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable


def sql_quote(path: Path | str) -> str:
    return str(path).replace("'", "''")


def command_exists(path_or_name: str | None) -> str | None:
    if not path_or_name:
        return None
    p = Path(path_or_name)
    if p.exists() and os.access(p, os.X_OK):
        return str(p.resolve())
    found = shutil.which(path_or_name)
    return found


def parse_time_verbose(stderr: str) -> float:
    for line in stderr.splitlines():
        if "Maximum resident set size" in line:
            try:
                return int(line.split(":", 1)[1].strip()) / 1024.0
            except (ValueError, IndexError):
                return 0.0
    return 0.0


def run_peak_rss(argv: list[str], timeout: int) -> tuple[str, float, float]:
    start = time.monotonic()
    result = subprocess.run(
        ["/usr/bin/time", "-v", *argv],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - start
    peak_rss_mb = parse_time_verbose(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            "command failed: " + subprocess.list2cmdline(argv) +
            f"\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout, elapsed, peak_rss_mb


def generate_rows(subjects: int, queries: int, seed: int):
    rng = random.Random(seed)
    chroms = ["chr1", "chr2", "chr3", "chr4"]
    chrom_span = 5_000_000

    subject_rows = []
    for i in range(subjects):
        chrom = chroms[i % len(chroms)]
        width = 50 + (i % 200)
        start = rng.randrange(0, chrom_span - width)
        subject_rows.append((chrom, start, start + width, f"s{i}"))

    query_rows = []
    for i in range(queries):
        chrom = chroms[(i * 3) % len(chroms)]
        width = 30 + (i % 120)
        start = rng.randrange(0, chrom_span - width)
        query_rows.append((chrom, start, start + width, f"q{i}"))

    return subject_rows, query_rows


def write_bed(path: Path, rows: Iterable[tuple[str, int, int, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def copy_bed_prefix(src: Path, dst: Path, limit: int | None) -> int:
    n = 0
    with src.open("r", encoding="utf-8") as inp, dst.open("w", encoding="utf-8") as out:
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            out.write(line)
            n += 1
            if limit and n >= limit:
                break
    return n


def count_bed_rows(path: Path) -> int:
    n = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def _worker_duckhts(args: argparse.Namespace) -> None:
    import duckdb

    con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    try:
        con.execute(f"LOAD '{sql_quote(Path(args.extension).resolve())}'")
        subject_path = sql_quote(Path(args.subject_bed).resolve())
        query_path = sql_quote(Path(args.query_bed).resolve())

        build_start = time.monotonic()
        con.execute(
            "CREATE TABLE cgr_subject AS "
            "SELECT chrom, start, \"end\" "
            f"FROM read_bed('{subject_path}')"
        )
        con.execute(
            "SELECT duckhts_cgranges_from_query("
            "  'bench_idx', "
            "  'SELECT chrom, start, \"end\" FROM cgr_subject', "
            "  'chrom', 'start', 'end'"
            ")"
        )
        # from_query already finalizes today; this remains harmless and documents
        # the intended index-ready state for readers of the benchmark.
        con.execute("SELECT duckhts_cgranges_index('bench_idx')")
        build_elapsed = time.monotonic() - build_start

        pass_seconds: list[float] = []
        matched: int | None = None
        total_hits: int | None = None
        for _ in range(args.passes):
            q0 = time.monotonic()
            if args.duckhts_mode == "filter":
                rows = con.execute(
                    "SELECT count(*)::BIGINT "
                    f"FROM read_bed('{query_path}') q "
                    "WHERE duckhts_cgranges_has_overlap('bench_idx', q.chrom, q.start, q.\"end\")"
                ).fetchone()
                pass_matched = int(rows[0])
                pass_hits = None
            elif args.duckhts_mode == "count":
                rows = con.execute(
                    "WITH counts AS ("
                    "  SELECT duckhts_cgranges_count_overlaps('bench_idx', q.chrom, q.start, q.\"end\") AS n "
                    f"  FROM read_bed('{query_path}') q"
                    ") "
                    "SELECT count(*) FILTER (WHERE n > 0)::BIGINT, COALESCE(sum(n), 0)::BIGINT FROM counts"
                ).fetchone()
                pass_matched = int(rows[0])
                pass_hits = int(rows[1])
            else:
                raise RuntimeError(f"unknown duckhts mode: {args.duckhts_mode}")
            pass_seconds.append(time.monotonic() - q0)
            if matched is None:
                matched = pass_matched
                total_hits = pass_hits
            elif matched != pass_matched or total_hits != pass_hits:
                raise RuntimeError("DuckHTS repeated passes produced inconsistent results")

        con.execute("SELECT duckhts_cgranges_destroy('bench_idx')")
        payload = {
            "tool": "duckhts",
            "variant": f"scalar_{args.duckhts_mode}",
            "subject_intervals": args.subject_intervals,
            "query_intervals": args.query_intervals,
            "passes": args.passes,
            "build_index_sec": build_elapsed,
            "query_total_sec": sum(pass_seconds),
            "query_pass_1_sec": pass_seconds[0] if pass_seconds else 0.0,
            "matched_query_intervals": matched or 0,
            "total_hits": total_hits,
        }
        Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    finally:
        con.close()


def _worker_bedtk(args: argparse.Namespace) -> None:
    bedtk = command_exists(args.bedtk)
    if not bedtk:
        raise SystemExit(f"bedtk executable not found: {args.bedtk}")

    pass_seconds = []
    matched: int | None = None
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
        "variant": "flt",
        "subject_intervals": args.subject_intervals,
        "query_intervals": args.query_intervals,
        "passes": args.passes,
        "build_index_sec": 0.0,
        "query_total_sec": sum(pass_seconds),
        "query_pass_1_sec": pass_seconds[0] if pass_seconds else 0.0,
        "matched_query_intervals": matched or 0,
        "total_hits": None,
    }
    Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _worker_bedtools(args: argparse.Namespace) -> None:
    bedtools = command_exists(args.bedtools)
    if not bedtools:
        raise SystemExit(f"bedtools executable not found: {args.bedtools}")

    pass_seconds = []
    matched: int | None = None
    total_hits: int | None = None
    for _ in range(args.passes):
        t0 = time.monotonic()
        if args.bedtools_mode == "u":
            proc = subprocess.run(
                [bedtools, "intersect", "-a", args.query_bed, "-b", args.subject_bed, "-u"],
                check=True,
                capture_output=True,
                text=True,
            )
            pass_matched = sum(1 for line in proc.stdout.splitlines() if line.strip())
            pass_hits = None
        elif args.bedtools_mode == "c":
            proc = subprocess.run(
                [bedtools, "intersect", "-a", args.query_bed, "-b", args.subject_bed, "-c"],
                check=True,
                capture_output=True,
                text=True,
            )
            pass_matched = 0
            pass_hits = 0
            for line in proc.stdout.splitlines():
                if not line.strip():
                    continue
                count = int(line.rstrip("\n").split("\t")[-1])
                if count > 0:
                    pass_matched += 1
                    pass_hits += count
        else:
            raise RuntimeError(f"unknown bedtools mode: {args.bedtools_mode}")
        pass_seconds.append(time.monotonic() - t0)
        if matched is None:
            matched = pass_matched
            total_hits = pass_hits
        elif matched != pass_matched or total_hits != pass_hits:
            raise RuntimeError("bedtools repeated passes produced inconsistent results")

    payload = {
        "tool": "bedtools",
        "variant": f"intersect_{args.bedtools_mode}",
        "subject_intervals": args.subject_intervals,
        "query_intervals": args.query_intervals,
        "passes": args.passes,
        "build_index_sec": 0.0,
        "query_total_sec": sum(pass_seconds),
        "query_pass_1_sec": pass_seconds[0] if pass_seconds else 0.0,
        "matched_query_intervals": matched or 0,
        "total_hits": total_hits,
    }
    Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def prepare_inputs(args: argparse.Namespace, out_dir: Path) -> tuple[Path, Path, int, int, str]:
    if args.subject_bed and args.query_bed:
        subject_src = Path(args.subject_bed).resolve()
        query_src = Path(args.query_bed).resolve()
        subject_bed = out_dir / "subject.bed"
        query_bed = out_dir / "query.bed"
        subject_intervals = copy_bed_prefix(subject_src, subject_bed, args.limit_subjects)
        query_intervals = copy_bed_prefix(query_src, query_bed, args.limit_queries)
        label = args.label or "provided_bed"
        return subject_bed, query_bed, subject_intervals, query_intervals, label

    subject_rows, query_rows = generate_rows(args.subjects, args.queries, args.seed)
    subject_bed = out_dir / "subject.bed"
    query_bed = out_dir / "query.bed"
    write_bed(subject_bed, subject_rows)
    write_bed(query_bed, query_rows)
    return subject_bed, query_bed, len(subject_rows), len(query_rows), args.label or "synthetic"


def append_result(results: list[dict], payload_path: Path, elapsed: float, rss: float) -> None:
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    payload["total_elapsed_sec"] = elapsed
    payload["peak_rss_mb"] = rss
    results.append(payload)


def bench_main(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    subject_bed, query_bed, subject_intervals, query_intervals, dataset_label = prepare_inputs(args, out_dir)

    print("DuckHTS cgranges scalar-provider benchmark")
    print(f"Dataset           : {dataset_label}")
    print(f"Subject BED       : {subject_bed} ({subject_intervals:,} intervals)")
    print(f"Query BED         : {query_bed} ({query_intervals:,} intervals)")
    print(f"Passes            : {args.passes}")
    print(f"Extension         : {Path(args.extension).resolve()}")
    print(f"bedtk             : {command_exists(args.bedtk) or 'not found'}")
    print(f"bedtools          : {command_exists(args.bedtools) or 'not found'}")

    script_path = Path(__file__).resolve()
    results: list[dict] = []

    jobs: list[tuple[str, Path, list[str]]] = []
    for mode in ("filter", "count"):
        json_out = out_dir / f"duckhts_scalar_{mode}.json"
        jobs.append((
            f"DuckHTS scalar {mode}",
            json_out,
            [
                sys.executable, str(script_path),
                "--worker-duckhts", "--duckhts-mode", mode,
                "--extension", str(Path(args.extension).resolve()),
                "--subject-bed", str(subject_bed),
                "--query-bed", str(query_bed),
                "--passes", str(args.passes),
                "--subject-intervals", str(subject_intervals),
                "--query-intervals", str(query_intervals),
                "--json-out", str(json_out),
            ],
        ))

    if command_exists(args.bedtk):
        json_out = out_dir / "bedtk_flt.json"
        jobs.append((
            "bedtk flt",
            json_out,
            [
                sys.executable, str(script_path),
                "--worker-bedtk",
                "--bedtk", str(command_exists(args.bedtk)),
                "--subject-bed", str(subject_bed),
                "--query-bed", str(query_bed),
                "--passes", str(args.passes),
                "--subject-intervals", str(subject_intervals),
                "--query-intervals", str(query_intervals),
                "--json-out", str(json_out),
            ],
        ))

    bedtools_path = command_exists(args.bedtools)
    if bedtools_path:
        for mode in ("u", "c"):
            json_out = out_dir / f"bedtools_intersect_{mode}.json"
            jobs.append((
                f"bedtools intersect -{mode}",
                json_out,
                [
                    sys.executable, str(script_path),
                    "--worker-bedtools", "--bedtools-mode", mode,
                    "--bedtools", bedtools_path,
                    "--subject-bed", str(subject_bed),
                    "--query-bed", str(query_bed),
                    "--passes", str(args.passes),
                    "--subject-intervals", str(subject_intervals),
                    "--query-intervals", str(query_intervals),
                    "--json-out", str(json_out),
                ],
            ))

    for label, json_out, argv in jobs:
        print(f"\n--- {label} ---")
        _, elapsed, rss = run_peak_rss(argv, args.timeout)
        append_result(results, json_out, elapsed, rss)
        row = results[-1]
        print(
            f"{row['tool']}:{row['variant']} "
            f"query={row['query_total_sec']:.3f}s total={row['total_elapsed_sec']:.3f}s "
            f"rss={row['peak_rss_mb']:.1f}MB matched={row['matched_query_intervals']} "
            f"hits={row.get('total_hits')}"
        )

    # Semantic checks across comparable outputs.
    filter_rows = [r for r in results if r["variant"] in {"scalar_filter", "flt", "intersect_u"}]
    if filter_rows:
        expected = filter_rows[0]["matched_query_intervals"]
        for row in filter_rows[1:]:
            if row["matched_query_intervals"] != expected:
                raise RuntimeError(
                    "matched-query mismatch: "
                    + ", ".join(f"{r['tool']}:{r['variant']}={r['matched_query_intervals']}" for r in filter_rows)
                )
    count_rows = [r for r in results if r["variant"] in {"scalar_count", "intersect_c"}]
    if len(count_rows) > 1:
        expected_hits = count_rows[0].get("total_hits")
        for row in count_rows[1:]:
            if row.get("total_hits") != expected_hits:
                raise RuntimeError(
                    "overlap-count mismatch: "
                    + ", ".join(f"{r['tool']}:{r['variant']}={r.get('total_hits')}" for r in count_rows)
                )

    csv_path = out_dir / "summary.csv"
    fieldnames = [
        "tool", "variant", "subject_intervals", "query_intervals", "passes",
        "build_index_sec", "query_total_sec", "query_pass_1_sec", "total_elapsed_sec",
        "peak_rss_mb", "matched_query_intervals", "total_hits",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow({k: row.get(k) for k in fieldnames})

    metadata = {
        "dataset": dataset_label,
        "subject_bed": str(subject_bed),
        "query_bed": str(query_bed),
        "subject_intervals": subject_intervals,
        "query_intervals": query_intervals,
        "passes": args.passes,
        "bedtk_available": command_exists(args.bedtk) is not None,
        "bedtools_available": bedtools_path is not None,
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("\nSummary")
    for row in results:
        print(
            f"- {row['tool']:8s} {row['variant']:16s} "
            f"build={row['build_index_sec']:7.3f}s "
            f"query={row['query_total_sec']:8.3f}s "
            f"total={row['total_elapsed_sec']:8.3f}s "
            f"rss={row['peak_rss_mb']:7.1f} MB "
            f"matched={row['matched_query_intervals']} hits={row.get('total_hits')}"
        )
    print(f"\nWrote {csv_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark DuckHTS cgranges scalar probes against bedtk/bedtools")
    parser.add_argument("--extension", default="build/release/duckhts.duckdb_extension")
    parser.add_argument("--bedtk", default=".sync/bedtk/bedtk")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--subjects", type=int, default=50000)
    parser.add_argument("--queries", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--passes", type=int, default=3)
    parser.add_argument("--subject-bed")
    parser.add_argument("--query-bed")
    parser.add_argument("--limit-subjects", type=int)
    parser.add_argument("--limit-queries", type=int)
    parser.add_argument("--label")
    parser.add_argument("--out-dir", default=".tmp/cgranges_benchmark")
    parser.add_argument("--timeout", type=int, default=3600)

    parser.add_argument("--worker-duckhts", action="store_true")
    parser.add_argument("--worker-bedtk", action="store_true")
    parser.add_argument("--worker-bedtools", action="store_true")
    parser.add_argument("--duckhts-mode", choices=["filter", "count"], default="filter")
    parser.add_argument("--bedtools-mode", choices=["u", "c"], default="u")
    parser.add_argument("--subject-intervals", type=int, default=0)
    parser.add_argument("--query-intervals", type=int, default=0)
    parser.add_argument("--json-out")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.worker_duckhts:
        _worker_duckhts(args)
        return
    if args.worker_bedtk:
        _worker_bedtk(args)
        return
    if args.worker_bedtools:
        _worker_bedtools(args)
        return
    bench_main(args)


if __name__ == "__main__":
    main()
