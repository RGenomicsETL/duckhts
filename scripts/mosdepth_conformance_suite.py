#!/usr/bin/env python3
"""
Run a native duckhts_mosdepth conformance matrix over local fixture files.

This wraps scripts/mosdepth_conformance.py so we can exercise a small but more
diverse fixture set, including upstream mosdepth edge cases that are now
vendored into test/data/.
"""

import argparse
import os
import subprocess
import sys


CASES = [
    {
        "name": "range_bam",
        "bam": "test/data/range.bam",
        "mode": "fast",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "thresholds": "1,2",
    },
    {
        "name": "range_cram",
        "bam": "test/data/range.cram",
        "mode": "fast",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "thresholds": "1,2",
        "fasta": "test/data/ce.fa",
    },
    {
        "name": "big_csi",
        "bam": "test/data/big.bam",
        "mode": "fast",
        "chrom": "ref",
        "window_size": "100000",
    },
    {
        "name": "empty_tids",
        "bam": "test/data/empty-tids.bam",
        "mode": "fast",
        "chrom": "HPV18",
        "window_size": "1000",
        "thresholds": "1,10,100",
    },
    {
        "name": "overlapping_pairs_fast",
        "bam": "test/data/overlapping-pairs.bam",
        "mode": "fast",
        "chrom": "1",
        "window_size": "1000",
    },
    {
        "name": "empty_tids_bed_quantize",
        "bam": "test/data/empty-tids.bam",
        "mode": "fast",
        "chrom": "HPV18",
        "by_bed": "test/data/empty-tids.bed",
        "thresholds": "1,10,100",
        "quantize": ":1:4:",
    },
    {
        "name": "range_bam_default",
        "bam": "test/data/range.bam",
        "mode": "default",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "thresholds": "1,2",
    },
    {
        "name": "range_bam_fragment",
        "bam": "test/data/range.bam",
        "mode": "fragment",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "thresholds": "1,2",
    },
    {
        "name": "range_cram_default",
        "bam": "test/data/range.cram",
        "mode": "default",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "thresholds": "1,2",
        "fasta": "test/data/ce.fa",
    },
    {
        "name": "range_cram_fragment",
        "bam": "test/data/range.cram",
        "mode": "fragment",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "thresholds": "1,2",
        "fasta": "test/data/ce.fa",
    },
    {
        "name": "overlapping_pairs_default",
        "bam": "test/data/overlapping-pairs.bam",
        "mode": "default",
        "chrom": "1",
        "window_size": "1000",
    },
    {
        "name": "range_read_groups_missing",
        "bam": "test/data/range.bam",
        "mode": "default",
        "chrom": "CHROMOSOME_I",
        "window_size": "1000",
        "read_groups": "missing",
    },
    {
        "name": "range_mapq_include_flag",
        "bam": "test/data/range.bam",
        "mode": "default",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "mapq": 60,
        "include_flag": 64,
    },
    {
        "name": "range_exclude_proper_pairs",
        "bam": "test/data/range.bam",
        "mode": "default",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "flag": 2,
    },
    {
        "name": "range_use_median",
        "bam": "test/data/range.bam",
        "mode": "default",
        "chrom": "CHROMOSOME_II",
        "window_size": "1000",
        "use_median": True,
    },
]


def run_case(repo_root, extension, threads, case):
    cmd = [
        sys.executable,
        os.path.join(repo_root, "scripts", "mosdepth_conformance.py"),
        os.path.join(repo_root, case["bam"]),
        "--mode",
        case.get("mode", "fast"),
        "--chrom",
        case["chrom"],
        "--threads",
        str(threads),
        "--extension",
        extension,
    ]
    if "window_size" in case:
        cmd.extend(["--window-size", case["window_size"]])
    if "by_bed" in case:
        cmd.extend(["--by-bed", os.path.join(repo_root, case["by_bed"])])
    if "thresholds" in case:
        cmd.extend(["--thresholds", case["thresholds"]])
    if "quantize" in case:
        cmd.extend(["--quantize", case["quantize"]])
    if "read_groups" in case:
        cmd.extend(["--read-groups", case["read_groups"]])
    if "flag" in case:
        cmd.extend(["--flag", str(case["flag"])])
    if "include_flag" in case:
        cmd.extend(["--include-flag", str(case["include_flag"])])
    if "mapq" in case:
        cmd.extend(["--mapq", str(case["mapq"])])
    if case.get("use_median"):
        cmd.append("--use-median")
    if "fasta" in case:
        cmd.extend(["--fasta", os.path.join(repo_root, case["fasta"])])

    print(f"\n=== {case['name']} ===")
    result = subprocess.run(cmd, cwd=repo_root)
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="Run native mosdepth conformance fixture suite")
    parser.add_argument(
        "--extension",
        default=os.path.join("build", "debug", "duckhts.duckdb_extension"),
        help="Path to duckhts.duckdb_extension (default: build/debug/duckhts.duckdb_extension)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Decompression threads for both tools (default: 4)",
    )
    args = parser.parse_args()

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    extension = os.path.abspath(os.path.join(repo_root, args.extension))
    if not os.path.exists(extension):
        sys.exit(f"Extension not found: {extension}")

    passed = 0
    for case in CASES:
        if run_case(repo_root, extension, args.threads, case):
            passed += 1

    total = len(CASES)
    print(f"\nSuite results: {passed}/{total} cases passed")
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
