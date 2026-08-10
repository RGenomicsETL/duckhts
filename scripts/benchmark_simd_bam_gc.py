#!/usr/bin/env python3
"""Benchmark DuckHTS seq_gc_content() on real BAM read sequences.

This complements the synthetic seq_gc_content SIMD benchmark with a BAM-backed
workload.  It reports two measurements per backend request:

* bam_scan: read BAM records and compute seq_gc_content(SEQ) in the same query.
* materialized_seq: materialize BAM SEQ strings once, then time only the real
  BAM sequence strings flowing through seq_gc_content().

Each backend request runs in a fresh Python/DuckDB process so process-wide SIMD
selection is clean and comparable across requests.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def source_sql(bam: str, max_reads: int) -> str:
    base = f"SELECT SEQ FROM read_bam({sql_string(bam)}) WHERE SEQ IS NOT NULL"
    if max_reads > 0:
        base += f" LIMIT {max_reads}"
    return f"({base})"


def aggregate_sql(source: str) -> str:
    return (
        "SELECT count(*) AS reads, "
        "coalesce(sum(length(SEQ)), 0) AS total_bases, "
        "coalesce(sum(coalesce(seq_gc_content(SEQ), 0.0)), 0.0) AS gc_fraction_sum "
        f"FROM {source}"
    )


def as_result_tuple(row: tuple[Any, Any, Any]) -> tuple[int, int, float]:
    return (int(row[0]), int(row[1]), float(row[2]))


def benchmark_query(con: Any, query: str, iterations: int) -> tuple[tuple[int, int, float], list[float]]:
    expected = as_result_tuple(con.execute(query).fetchone())
    timings: list[float] = []
    for _ in range(iterations):
        t0 = time.perf_counter()
        observed = as_result_tuple(con.execute(query).fetchone())
        timings.append(time.perf_counter() - t0)
        if observed[:2] != expected[:2] or abs(observed[2] - expected[2]) > 1e-9:
            raise RuntimeError(f"unstable aggregate result: {observed!r} != {expected!r}")
    return expected, timings


def make_record(
    *,
    benchmark: str,
    args: argparse.Namespace,
    requested: str | None,
    selected: str | None,
    kernel_row: tuple[Any, Any, Any] | None,
    available: bool,
    skipped: bool,
    skip_reason: str | None = None,
    load_sec: float | None = None,
    result: tuple[int, int, float] | None = None,
    timings: list[float] | None = None,
) -> dict[str, Any]:
    timings = timings or []
    reads = result[0] if result else None
    total_bases = result[1] if result else None
    median_sec = statistics.median(timings) if timings else None
    min_sec = min(timings) if timings else None
    return {
        "benchmark": benchmark,
        "backend_request": args.backend,
        "requested_backend": requested,
        "selected_backend": selected,
        "kernel_backend": kernel_row[0] if kernel_row else None,
        "kernel_capability": kernel_row[1] if kernel_row else None,
        "kernel_scalar_fallback": bool(kernel_row[2]) if kernel_row else None,
        "available": available,
        "skipped": skipped,
        "skip_reason": skip_reason,
        "bam_path": str(Path(args.bam).resolve()),
        "max_reads": args.max_reads,
        "threads": args.threads,
        "reads": reads,
        "total_bases": total_bases,
        "iterations": args.iterations,
        "load_sec": load_sec,
        "gc_fraction_sum": result[2] if result else None,
        "timings_sec": timings,
        "median_sec": median_sec,
        "min_sec": min_sec,
        "reads_per_sec_median": (reads / median_sec) if reads is not None and median_sec else None,
        "mbases_per_sec_median": ((total_bases / median_sec) / 1e6) if total_bases is not None and median_sec else None,
    }


def worker(args: argparse.Namespace) -> int:
    import duckdb  # imported only in worker processes

    con = duckdb.connect(":memory:", config={"allow_unsigned_extensions": "true"})
    con.execute(f"PRAGMA threads={args.threads}")
    con.execute(f"LOAD {sql_string(str(Path(args.extension).resolve()))}")

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    available = args.backend == "auto" or bool(
        con.execute("SELECT duckhts_simd_backend_available(?)", [args.backend]).fetchone()[0]
    )
    if not available:
        records = [
            make_record(
                benchmark=mode,
                args=args,
                requested=None,
                selected=None,
                kernel_row=None,
                available=False,
                skipped=True,
                skip_reason="backend request is not available in this process",
            )
            for mode in modes
        ]
        print(json.dumps(records, sort_keys=True))
        return 0

    selected_after_set = con.execute(
        "SELECT backend FROM duckhts_simd_set_backend(?)",
        [args.backend],
    ).fetchone()[0]
    selected = con.execute("SELECT duckhts_simd_backend()").fetchone()[0]
    requested = con.execute("SELECT duckhts_simd_requested_backend()").fetchone()[0]
    kernel_row = con.execute(
        "SELECT selected_backend, selected_capability, scalar_fallback "
        "FROM duckhts_simd_kernel_info() WHERE kernel = 'seq_base_counts'"
    ).fetchone()
    if selected_after_set != selected:
        raise RuntimeError(f"set_backend returned {selected_after_set!r}, diagnostics returned {selected!r}")

    src = source_sql(args.bam, args.max_reads)
    records: list[dict[str, Any]] = []

    if "bam_scan" in modes:
        result, timings = benchmark_query(con, aggregate_sql(src), args.iterations)
        records.append(
            make_record(
                benchmark="bam_scan",
                args=args,
                requested=requested,
                selected=selected,
                kernel_row=kernel_row,
                available=True,
                skipped=False,
                result=result,
                timings=timings,
            )
        )

    if "materialized_seq" in modes:
        con.execute("DROP TABLE IF EXISTS bam_gc_seqs")
        t0 = time.perf_counter()
        con.execute(f"CREATE TEMP TABLE bam_gc_seqs AS SELECT SEQ FROM {src}")
        load_sec = time.perf_counter() - t0
        result, timings = benchmark_query(con, aggregate_sql("bam_gc_seqs"), args.iterations)
        records.append(
            make_record(
                benchmark="materialized_seq",
                args=args,
                requested=requested,
                selected=selected,
                kernel_row=kernel_row,
                available=True,
                skipped=False,
                load_sec=load_sec,
                result=result,
                timings=timings,
            )
        )

    print(json.dumps(records, sort_keys=True))
    return 0


def run_backend(args: argparse.Namespace, backend: str) -> list[dict[str, Any]]:
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--extension",
        args.extension,
        "--bam",
        args.bam,
        "--max-reads",
        str(args.max_reads),
        "--iterations",
        str(args.iterations),
        "--threads",
        str(args.threads),
        "--modes",
        args.modes,
        "--backend",
        backend,
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"backend {backend!r} failed with exit {proc.returncode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout.strip().splitlines()[-1])


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", default="build/release/duckhts.duckdb_extension")
    parser.add_argument("--bam", default="test/data/nanopore.bam")
    parser.add_argument("--max-reads", type=int, default=0, help="0 means all reads")
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--modes", default="bam_scan,materialized_seq")
    parser.add_argument("--backends", default="scalar,auto,avx2,avx512")
    parser.add_argument("--backend", default="auto", help=argparse.SUPPRESS)
    parser.add_argument("--json-out")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args(argv)

    if args.max_reads < 0:
        raise SystemExit("--max-reads must be >= 0")
    if args.iterations < 1:
        raise SystemExit("--iterations must be positive")
    if args.threads < 1:
        raise SystemExit("--threads must be positive")
    if args.worker:
        return worker(args)

    results: list[dict[str, Any]] = []
    for backend in [b.strip() for b in args.backends.split(",") if b.strip()]:
        results.extend(run_backend(args, backend))

    scalar_rates = {
        r["benchmark"]: r["mbases_per_sec_median"]
        for r in results
        if not r.get("skipped") and r.get("backend_request") == "scalar"
    }

    print("benchmark\tbackend_request\tdispatch_label\tseq_kernel\tstatus\tmedian_s\tMbases/s\tspeedup_vs_scalar")
    for r in results:
        if r.get("skipped"):
            print(f"{r['benchmark']}\t{r['backend_request']}\t-\t-\tskipped\t-\t-\t-")
            continue
        scalar_rate = scalar_rates.get(r["benchmark"])
        speedup = ""
        if scalar_rate and r["mbases_per_sec_median"]:
            speedup = f"{r['mbases_per_sec_median'] / scalar_rate:.3f}"
        print(
            f"{r['benchmark']}\t{r['backend_request']}\t{r['selected_backend']}\t{r.get('kernel_backend', '-')}\tok\t"
            f"{r['median_sec']:.6f}\t{r['mbases_per_sec_median']:.1f}\t{speedup}"
        )

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
