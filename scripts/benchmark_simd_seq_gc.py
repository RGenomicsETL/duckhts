#!/usr/bin/env python3
"""Benchmark DuckHTS seq_gc_content() across explicit SIMD backend requests.

Each requested backend is run in a fresh Python/DuckDB process for clean state,
but backend selection itself is done through SQL with duckhts_simd_set_backend().
Unsupported backends are reported as skipped instead of failing, so the same
benchmark can run on x86, ARM, and wasm-capable builds with different SIMD sets.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path


def worker(args: argparse.Namespace) -> int:
    import duckdb  # imported only in the worker process

    con = duckdb.connect(":memory:", config={"allow_unsigned_extensions": "true"})
    con.execute(f"LOAD {sql_string(str(Path(args.extension).resolve()))}")

    # duckhts_simd_backend_available() is for concrete backends.  "auto" is
    # a request that selects the best available backend, so measure it directly.
    available = args.backend == "auto" or bool(
        con.execute(
            "SELECT duckhts_simd_backend_available(?)",
            [args.backend],
        ).fetchone()[0]
    )
    if not available:
        print(json.dumps({
            "backend_request": args.backend,
            "requested_backend": None,
            "selected_backend": None,
            "available": False,
            "skipped": True,
            "skip_reason": "backend request is not available in this process",
            "rows": args.rows,
            "seq_len": args.seq_len,
            "total_bases": args.rows * args.seq_len,
            "iterations": args.iterations,
            "timings_sec": [],
        }, sort_keys=True))
        return 0

    selected_after_set = con.execute(
        "SELECT duckhts_simd_set_backend(?)",
        [args.backend],
    ).fetchone()[0]
    selected = con.execute("SELECT duckhts_simd_backend()").fetchone()[0]
    requested = con.execute("SELECT duckhts_simd_requested_backend()").fetchone()[0]
    if selected_after_set != selected:
        raise RuntimeError(f"set_backend returned {selected_after_set!r}, diagnostics returned {selected!r}")

    pattern = (args.pattern * ((args.seq_len // len(args.pattern)) + 1))[: args.seq_len]
    con.execute("CREATE TEMP TABLE seqs(seq VARCHAR)")
    # Keep data generation out of the timed loop and materialize real strings so
    # constant folding cannot hide seq_gc_content() work.
    con.executemany("INSERT INTO seqs VALUES (?)", [(pattern,) for _ in range(args.rows)])
    expected = con.execute("SELECT sum(seq_gc_content(seq)) FROM seqs").fetchone()[0]

    timings = []
    for _ in range(args.iterations):
        t0 = time.perf_counter()
        value = con.execute("SELECT sum(seq_gc_content(seq)) FROM seqs").fetchone()[0]
        timings.append(time.perf_counter() - t0)
        if value != expected:
            raise RuntimeError(f"unstable result: {value!r} != {expected!r}")

    total_bases = args.rows * args.seq_len
    median_sec = statistics.median(timings)
    result = {
        "backend_request": args.backend,
        "requested_backend": requested,
        "selected_backend": selected,
        "available": True,
        "skipped": False,
        "rows": args.rows,
        "seq_len": args.seq_len,
        "total_bases": total_bases,
        "iterations": args.iterations,
        "result": expected,
        "timings_sec": timings,
        "median_sec": median_sec,
        "min_sec": min(timings),
        "mbases_per_sec_median": (total_bases / median_sec) / 1e6,
    }
    print(json.dumps(result, sort_keys=True))
    return 0


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def run_backend(args: argparse.Namespace, backend: str) -> dict:
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--extension",
        args.extension,
        "--rows",
        str(args.rows),
        "--seq-len",
        str(args.seq_len),
        "--iterations",
        str(args.iterations),
        "--pattern",
        args.pattern,
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
    parser.add_argument("--rows", type=int, default=200_000)
    parser.add_argument("--seq-len", type=int, default=512)
    parser.add_argument("--iterations", type=int, default=7)
    parser.add_argument("--pattern", default="ACGTNNacgtnn")
    parser.add_argument("--backends", default="scalar,auto,avx2,avx512")
    parser.add_argument("--backend", default="auto", help=argparse.SUPPRESS)
    parser.add_argument("--json-out")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args(argv)

    if args.worker:
        return worker(args)

    results = []
    for backend in [b.strip() for b in args.backends.split(",") if b.strip()]:
        results.append(run_backend(args, backend))

    scalar = next((r for r in results if not r.get("skipped") and r.get("backend_request") == "scalar"), None)
    scalar_mbps = scalar["mbases_per_sec_median"] if scalar else None

    print("backend_request\tselected\tstatus\tmedian_s\tMbases/s\tspeedup_vs_scalar")
    for r in results:
        if r.get("skipped"):
            print(f"{r['backend_request']}\t-\tskipped\t-\t-\t-")
            continue
        speedup = ""
        if scalar_mbps and r["mbases_per_sec_median"]:
            speedup = f"{r['mbases_per_sec_median'] / scalar_mbps:.3f}"
        print(
            f"{r['backend_request']}\t{r['selected_backend']}\tok\t"
            f"{r['median_sec']:.6f}\t{r['mbases_per_sec_median']:.1f}\t{speedup}"
        )

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
