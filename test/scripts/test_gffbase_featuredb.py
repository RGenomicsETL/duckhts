#!/usr/bin/env python3
"""Exercise physical duplicate IDs and keyed feature-database query parity."""

import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import duckdb

ROOT = Path(__file__).resolve().parents[2]
DRIVER = ROOT / "scripts/gffbase_featuredb_benchmark.py"
EXTENSION = ROOT / "build/release/duckhts.duckdb_extension"


def build(engine, source, database):
    output = subprocess.run(
        [sys.executable, str(DRIVER), "--worker", engine, "--input", str(source),
         "--database", str(database), "--extension", str(EXTENSION)],
        capture_output=True, text=True, check=True
    )
    return json.loads(output.stdout)


def main():
    with tempfile.TemporaryDirectory(prefix="gffbase-featuredb-test-") as temporary:
        root = Path(temporary)
        source = root / "duplicates.gff3.gz"
        with gzip.open(source, "wt") as output:
            output.write("##gff-version 3\n")
            output.write("chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1;Alias=a;Alias=b\n")
            output.write("chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n")
            output.write("chr1\ttest\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t1\n")
            output.write("chr1\ttest\texon\t20\t30\t.\t+\t.\tID=e1;Parent=t1\n")
        duck = build("DuckHTS", source, root / "duck.duckdb")
        base = build("gffbase", source, root / "base.duckdb")
        for key in ("features", "anchors", "windows", "anchors_sha256",
                    "windows_sha256", "descendants", "descendants_sha256",
                    "region_hits", "region_sha256"):
            assert duck[key] == base[key], (key, duck[key], base[key])
        assert duck["features"] == 4
        assert duck["descendants"] == 3
        assert duck["region_hits"] == 4
        with duckdb.connect(str(root / "duck.duckdb"), read_only=True) as database:
            attributes = database.execute("SELECT attributes FROM features WHERE id = 'g1'").fetchone()[0]
        assert [(pair["value"], pair["idx"]) for pair in attributes
                if pair["key"] == "Alias"] == [("a", 0), ("b", 1)]
    print("feature database duplicate-ID parity: passed")


if __name__ == "__main__":
    main()
