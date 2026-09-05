#!/usr/bin/env python3
"""Fail every direct DuckHTS duckdb_malloc call on real reader query paths.

The shim is test-only and POSIX-only; no global malloc interposition. Exact
successful schemas/rows must survive recovery. The ledger covers first-party
DuckDB-owned allocations, not HTSlib/libc buffers or DuckDB output storage.
ASan/UBSan additionally check invalid access and cleanup in the full extension.
"""
import argparse
import ctypes
from pathlib import Path

import duckdb


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--probe", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]

    def source(reader, fixture, options=""):
        path = str(root / "test/data" / fixture).replace("'", "''")
        return f"{reader}('{path}'{options})"

    cases = {
        "bcf-indexed-wide": source("read_bcf", "vcf_file.bcf"),
        "bcf-sequential-wide": source("read_bcf", "bcf_cache_lifecycle.vcf", ", scan_mode:='sequential'"),
        "bcf-tidy": source("read_bcf", "bcf_cache_lifecycle.vcf", ", tidy_format:=true"),
        "bcf-fallback-gt": source("read_bcf", "default_gt.vcf"),
        "bcf-regions": source("read_bcf", "region_union.bcf", ", region:='chr1:11-11,chr1:19-19'"),
        "vcf-tabix": source("read_bcf", "region_union.vcf.gz"),
        "bam-indexed": source("read_bam", "bam_scan_mixed.bam", ", decompression_threads:=0"),
        "bam-read-groups": source("read_bam", "bam_read_groups.sam", ", scan_mode:='sequential', decompression_threads:=0"),
        "cram-indexed": source("read_bam", "bam_scan_mixed.cram", ", decompression_threads:=0"),
        "fasta-sequential": source("read_fasta", "region_names.fa", ", scan_mode:='sequential'"),
        "fasta-indexed": source("read_fasta", "region_names.fa.gz", ", region:='{chr,part}:1-5,{chr:part}:1-5'"),
        "fastq": source("read_fastq", "r1.fq"),
        "fastq-paired": source("read_fastq", "r1.fq", ", mate_path:='" + str(root / "test/data/r2.fq") + "'"),
    }
    con = duckdb.connect(config={"allow_unsigned_extensions": "true", "threads": "1"})
    extension = str(args.extension.resolve())
    con.execute("LOAD '" + extension.replace("'", "''") + "'")
    probe = ctypes.CDLL(str(args.probe.resolve()))
    probe.reader_alloc_open.argtypes = [ctypes.c_char_p]
    probe.reader_alloc_arm.argtypes = [ctypes.c_long]
    for name in ("attempts", "failures", "live"):
        getattr(probe, "reader_alloc_" + name).restype = ctypes.c_long
    assert probe.reader_alloc_open(extension.encode()) == 0
    assert probe.reader_alloc_arm(0) == 0
    assert probe.reader_alloc_helper_checks() == 0
    probe.reader_alloc_disarm()

    def read(sql):
        result = con.execute(sql)
        schema = [(col[0], str(col[1])) for col in result.description]
        rows = result.fetchall()
        # Full physical multiset: no implicit order promise for parallel readers.
        return schema, sorted(rows, key=repr)

    failures = 0
    try:
        for name, relation in cases.items():
            for projection in ("*", "count(*)"):
                # FILE_OFFSET currently misreads non-BGZF SAM/CRAM handles.
                # https://github.com/RGenomicsETL/duckhts/issues/194
                # That separate transport defect is not a recovery oracle;
                # retain every biological field, including SAMPLE_ID.
                if name in ("cram-indexed", "bam-read-groups") and projection == "*":
                    projection = "* EXCLUDE(FILE_OFFSET)"
                sql = f"SELECT {projection} FROM {relation}"
                assert probe.reader_alloc_arm(0) == 0
                expected = read(sql)
                con.execute("SELECT 4242").fetchall()  # release the prior result
                count = probe.reader_alloc_attempts()
                assert probe.reader_alloc_live() == 0, (name, "control leak")
                probe.reader_alloc_disarm()
                assert count > 0, (name, "probe did not intercept allocations")
                for nth in range(1, count + 1):
                    print(f"{name} {projection}: fail {nth}/{count}", flush=True)
                    assert probe.reader_alloc_arm(nth) == 0
                    try:
                        read(sql)
                    except duckdb.Error as error:
                        assert "out of memory" in str(error).lower(), (name, nth, error)
                    else:
                        raise AssertionError((name, nth, "allocation failure silently succeeded"))
                    assert probe.reader_alloc_failures() == 1, (name, nth, "injection missed")
                    assert con.execute("SELECT 4242").fetchone() == (4242,)
                    assert probe.reader_alloc_live() == 0, (name, nth, "first-party leak")
                    probe.reader_alloc_disarm()
                    assert read(sql) == expected, (name, nth, "recovery changed schema/rows")
                    con.execute("SELECT 4242").fetchall()
                    failures += 1
    finally:
        probe.reader_alloc_close()
        con.close()
    print(f"reader allocation failures: {failures} errors, zero tracked leaks, exact schema/row recovery: OK")


if __name__ == "__main__":
    main()
