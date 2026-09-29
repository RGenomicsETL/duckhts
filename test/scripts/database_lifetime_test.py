"""A closed file database must be released: reopenable, and on Linux no descriptor."""

import gc
import os
import pathlib
import sys
import tempfile

import duckdb


ROOT = pathlib.Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "build/release/duckhts.duckdb_extension"
DATA = ROOT / "test/data"
CONFIG = {"allow_unsigned_extensions": "true"}


def open_descriptors(database_file):
    """Descriptors of this process that refer to database_file (Linux only)."""
    target = os.path.realpath(database_file)
    found = []
    for entry in os.listdir("/proc/self/fd"):
        try:
            if os.path.realpath(os.readlink(f"/proc/self/fd/{entry}")) == target:
                found.append(entry)
        except OSError:
            continue
    return found


def use_registered_functions(con, work):
    """Call every function family that keeps state between LOAD and close."""
    for _, sql in con.execute(
        "SELECT name, sql FROM duckhts_macro_definitions() ORDER BY install_order"
    ).fetchall():
        con.execute(sql)

    con.execute(
        "CREATE TEMP TABLE targets AS "
        "SELECT * FROM (VALUES ('chr1', 10, 20, 'a'), ('chr1', 15, 30, 'b')) "
        "AS t(chrom, start, \"end\", label)"
    )
    indexed = con.execute(
        "SELECT * FROM duckhts_cgranges_from_table("
        "'lifetime_idx', 'targets', 'chrom', 'start', 'end', 'label')"
    ).fetchone()
    assert indexed == (True,), indexed
    hits = con.execute(
        "SELECT count(*) FROM (SELECT unnest(duckhts_cgranges_overlaps_list("
        "'lifetime_idx', 'chr1', 16, 18)))"
    ).fetchone()
    assert hits == (2,), hits
    con.execute("SELECT duckhts_cgranges_destroy('lifetime_idx')")

    panel = work / "panel.parquet"
    con.execute(
        "COPY (SELECT 'WBcel235' AS assembly, 0::UBIGINT AS site_index, "
        "'CHROMOSOME_I' AS region, 914::UBIGINT AS position, "
        "'A' AS allele_a, 'C' AS allele_b) "
        f"TO '{panel}' (FORMAT parquet)"
    )
    counts = con.execute(
        "SELECT status FROM duckhts_somalier_bam_counts("
        f"'{DATA / 'range.bam'}', NULL, 'sample-1', '{DATA / 'ce.fa'}', "
        f"panel_parquet := '{panel}', "
        f"index_path := '{DATA / 'range.bam.bai'}', "
        f"reference_index_path := '{DATA / 'ce.fa.fai'}')"
    ).fetchall()
    assert counts == [("measured",)], counts


def main():
    with tempfile.TemporaryDirectory(prefix="duckhts_lifetime_") as tmp:
        work = pathlib.Path(tmp)
        database_file = work / "lifetime.duckdb"

        con = duckdb.connect(str(database_file), config=CONFIG)
        con.execute(f"LOAD '{EXTENSION}'")
        if sys.platform == "linux":
            assert open_descriptors(database_file), "database is open before close"
        use_registered_functions(con, work)
        con.close()
        del con
        gc.collect()

        if sys.platform == "linux":
            leaked = open_descriptors(database_file)
            assert not leaked, f"closed database still open as descriptors {leaked}"

        reopened = duckdb.connect(str(database_file), config=CONFIG)
        assert reopened.execute("SELECT 42").fetchone() == (42,)
        reopened.close()
    print("database lifetime: closed database released and reopened")


if __name__ == "__main__":
    main()
