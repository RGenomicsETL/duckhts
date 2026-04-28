#!/usr/bin/env python3
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import duckdb

REPO = Path(__file__).resolve().parents[2]
BEDTK = REPO / ".sync" / "bedtk" / "bedtk"
EXTENSION = REPO / "build" / "release" / "duckhts.duckdb_extension"

SUBJECT_ROWS = [
    ("chr1", 10, 20, "a"),
    ("chr1", 15, 25, "b"),
    ("chr1", 30, 40, "c"),
    ("chr2", 5, 10, "d"),
]
QUERY_ROWS = [
    ("chr1", 5, 12, "q1"),
    ("chr1", 18, 19, "q2"),
    ("chr1", 26, 29, "q3"),
    ("chr1", 30, 40, "q4"),
    ("chr2", 1, 7, "q5"),
]


def write_bed(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(row)


def run_bedtk(subject_bed: Path, query_bed: Path):
    proc = subprocess.run(
        [str(BEDTK), "flt", str(subject_bed), str(query_bed)],
        check=True,
        capture_output=True,
        text=True,
    )
    hits = set()
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        chrom, start, end, name = line.split("\t")[:4]
        hits.add((chrom, int(start), int(end), name))
    return hits


def run_duckhts(subject_rows, query_rows):
    con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    try:
        con.execute(f"LOAD '{str(EXTENSION).replace("'", "''")}'")
        con.execute(
            "CREATE TABLE cgr_subject(chrom VARCHAR, start BIGINT, \"end\" BIGINT, label VARCHAR)"
        )
        con.executemany("INSERT INTO cgr_subject VALUES (?, ?, ?, ?)", subject_rows)
        con.execute(
            "SELECT duckhts_cgranges_from_query('bedtk_idx', "
            "'SELECT chrom, start, \"end\", label FROM cgr_subject', "
            "'chrom', 'start', 'end', 'label')"
        )
        con.execute("SELECT duckhts_cgranges_index('bedtk_idx')")

        hits = set()
        for chrom, start, end, name in query_rows:
            rows = con.execute(
                "SELECT interval_ordinal FROM duckhts_cgranges_overlaps(?, ?, ?, ?)",
                ["bedtk_idx", chrom, start, end],
            ).fetchall()
            if rows:
                hits.add((chrom, start, end, name))
        return hits
    finally:
        con.close()


def main() -> int:
    if not BEDTK.exists():
        raise SystemExit(f"bedtk binary not found: {BEDTK}")
    if not EXTENSION.exists():
        raise SystemExit(f"extension not found: {EXTENSION}; build release first")

    with tempfile.TemporaryDirectory(prefix="duckhts_cgranges_bedtk_") as tmpdir:
        tmpdir = Path(tmpdir)
        subject_bed = tmpdir / "subject.bed"
        query_bed = tmpdir / "query.bed"
        write_bed(subject_bed, SUBJECT_ROWS)
        write_bed(query_bed, QUERY_ROWS)

        bedtk_hits = run_bedtk(subject_bed, query_bed)
        duckhts_hits = run_duckhts(SUBJECT_ROWS, QUERY_ROWS)

    if bedtk_hits != duckhts_hits:
        print("bedtk hits:", sorted(bedtk_hits), file=sys.stderr)
        print("duckhts hits:", sorted(duckhts_hits), file=sys.stderr)
        raise SystemExit("duckhts_cgranges overlap-existence parity failed")

    print("duckhts_cgranges overlap-existence parity passed")
    print(f"matched_query_intervals={len(duckhts_hits)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
