#!/usr/bin/env python3
"""Check actual DuckDB child-list reservations for a long numeric CIGAR token."""
import argparse
import ctypes
from pathlib import Path

import duckdb


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--probe", required=True, type=Path)
    args = parser.parse_args()
    extension = str(args.extension.resolve())

    con = duckdb.connect(config={"allow_unsigned_extensions": "true", "threads": "1"})
    con.execute("LOAD '" + extension.replace("'", "''") + "'")
    probe = ctypes.CDLL(str(args.probe.resolve()))
    probe.reader_alloc_open.argtypes = [ctypes.c_char_p]
    probe.reader_list_arm.argtypes = [ctypes.c_long, ctypes.c_long]
    probe.reader_list_attempts.restype = ctypes.c_long
    probe.reader_list_max_reserve.restype = ctypes.c_uint64
    assert probe.reader_alloc_open(extension.encode()) == 0
    try:
        for cigar, label in (
            ("repeat('0', 2000000) || '1M'", "text"),
            ("[16]::UINTEGER[]", "packed"),
        ):
            assert probe.reader_list_arm(1, 0) == 0
            try:
                rows = con.execute(f"SELECT cigar_aligned_blocks({cigar}, 42)").fetchall()
                assert rows == [({"ref_start": [42], "query_start": [0], "width": [1]},)]
                assert probe.reader_list_attempts() == 3, label
                assert probe.reader_list_max_reserve() == 1, label
            finally:
                probe.reader_list_disarm()
    finally:
        probe.reader_alloc_close()
        con.close()
    print("CIGAR child-list reserve: three one-block reservations per overload: OK")


if __name__ == "__main__":
    main()
