#!/usr/bin/env python3
"""Convenience wrapper for a real BED cgranges benchmark.

This delegates to ``scripts/cgranges_benchmark.py`` and supplies the repo's
DuckBedQC GRCh38 exon/clinical-region BED files by default. The delegated
benchmark uses streaming ``read_bed(...)`` provider rows plus the vectorized
``duckhts_cgranges_has_overlap`` / ``duckhts_cgranges_count_overlaps`` scalar
helpers and ``duckhts_cgranges_overlaps_list`` + ``UNNEST`` expansion; it does
not generate giant UNION ALL bulk-query SQL.
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

from duckhts_cache import duckhts_cache_subdir


def registry_artifact_path(repo: Path, artifact_id: str) -> Path:
    configured = os.environ.get("DUCKHTSBENCH_REGISTRY")
    registry = Path(configured) if configured else repo / "r" / "duckhtsbench" / "inst" / "benchmark_registry.tsv"
    if not registry.is_file():
        raise FileNotFoundError(f"benchmark registry not found: {registry}")
    with registry.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["id"] == artifact_id:
                return duckhts_cache_subdir(*Path(row["cache_relpath"]).parts)
    raise ValueError(f"artifact {artifact_id!r} is not registered in {registry}")


def build_parser() -> argparse.ArgumentParser:
    repo = Path(__file__).resolve().parents[1]
    duckbedqc_dir = Path(os.environ["DUCKBEDQC_DIR"]) if "DUCKBEDQC_DIR" in os.environ else registry_artifact_path(repo, "duckbedqc_118fc21")
    p = argparse.ArgumentParser(description="Run the real DuckBedQC cgranges benchmark")
    p.add_argument("--extension", default=str(repo / "build" / "release" / "duckhts.duckdb_extension"))
    p.add_argument("--bedtk", default=str(repo / ".sync" / "bedtk" / "bedtk"))
    p.add_argument("--bedtools", default="bedtools")
    p.add_argument("--subject-bed", default=str(duckbedqc_dir / "data" / "GRCh38_exons.bed"))
    p.add_argument("--query-bed", default=str(duckbedqc_dir / "data" / "GRCh38_illumina_clinical_regions_v100.39.0.bed"))
    p.add_argument("--limit-subjects", type=int, help="Use only the first N subject intervals")
    p.add_argument("--limit-queries", type=int, help="Use only the first N query intervals")
    p.add_argument("--passes", type=int, default=1)
    p.add_argument("--out-dir", default=str(duckhts_cache_subdir("benchmarks", "cgranges-real")))
    p.add_argument("--timeout", type=int, default=7200)
    return p


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    args = build_parser().parse_args()
    cmd = [
        sys.executable,
        str(repo / "scripts" / "cgranges_benchmark.py"),
        "--extension", args.extension,
        "--bedtk", args.bedtk,
        "--bedtools", args.bedtools,
        "--subject-bed", args.subject_bed,
        "--query-bed", args.query_bed,
        "--passes", str(args.passes),
        "--out-dir", args.out_dir,
        "--timeout", str(args.timeout),
        "--label", "DuckBedQC_GRCh38_exons_vs_clinical_regions",
    ]
    if args.limit_subjects is not None:
        cmd.extend(["--limit-subjects", str(args.limit_subjects)])
    if args.limit_queries is not None:
        cmd.extend(["--limit-queries", str(args.limit_queries)])
    return subprocess.call(cmd, cwd=str(repo))


if __name__ == "__main__":
    raise SystemExit(main())
