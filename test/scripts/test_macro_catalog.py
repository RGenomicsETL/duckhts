"""Multi-process DuckHTS macro installation and file-catalog checks."""

import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile

import duckdb


ROOT = pathlib.Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "build/release/duckhts.duckdb_extension"
CONFIG = {"allow_unsigned_extensions": "true"}


def connect(path, read_only=False):
    return duckdb.connect(str(path), read_only=read_only, config=CONFIG)


def load(con):
    con.execute(f"LOAD '{EXTENSION}'")


def names(con, database):
    return {row[0] for row in con.execute(
        "SELECT function_name FROM duckdb_functions() "
        "WHERE function_type IN ('macro', 'table_macro') AND database_name = ?",
        [database],
    ).fetchall()}


def definitions(con):
    rows = con.execute(
        "SELECT name, public, sql, definitions_sha256 "
        "FROM duckhts_macro_definitions() ORDER BY install_order"
    ).fetchall()
    assert len(rows) == 36
    assert len({row[3] for row in rows}) == 1
    ordered = b"".join(
        bytes([int(exposed)]) +
        sql.replace("CREATE OR REPLACE TEMP MACRO ",
                    "CREATE OR REPLACE MACRO ", 1).encode() + b"\0"
        for _, exposed, sql, _ in rows
    )
    assert rows[0][3] == hashlib.sha256(ordered).hexdigest()
    return rows


def install(con):
    for _, _, sql, _ in definitions(con):
        con.execute(sql)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_file(path):
    con = connect(path)
    load(con)
    rows = definitions(con)
    catalog = con.execute("SELECT current_database()").fetchone()[0]
    assert names(con, catalog) & {row[0] for row in rows} == set()
    assert con.execute("SELECT duckhts_htslib_version()").fetchone()[0]
    con.close()


def clean_reopened(path):
    con = connect(path, read_only=True)
    load(con)
    rows = definitions(con)
    catalog = con.execute("SELECT current_database()").fetchone()[0]
    assert names(con, catalog) & {row[0] for row in rows} == set()
    con.close()


def writable(path):
    con = connect(path)
    con.execute("CREATE MACRO duckhts_quote_ident(x) AS 'old'")
    con.close()

    first = connect(path)
    second = connect(path)
    load(first)
    load(second)
    rows = definitions(first)
    public = {name for name, exposed, _, _ in rows if exposed}
    catalog_entries = json.loads((ROOT / "functions.yaml").read_text())["functions"]
    declared = {entry["name"] for entry in catalog_entries
                if entry["kind"] in ("scalar_macro", "table_macro")}
    assert public == declared
    catalog = first.execute("SELECT current_database()").fetchone()[0]
    assert names(first, catalog) & {row[0] for row in rows} == {"duckhts_quote_ident"}
    assert names(first, "temp") & public == set()
    assert first.execute("SELECT duckhts_htslib_version()").fetchone()[0]
    install(first)
    install(first)
    assert names(first, "temp") & public == public
    assert first.execute("SELECT duckhts_quote_ident('a')").fetchone()[0] == '"a"'
    assert second.execute("SELECT duckhts_quote_ident('a')").fetchone()[0] == "old"
    assert names(second, "temp") & public == set()
    second.execute("CREATE TEMP TABLE inputs AS SELECT 'a' AS column_name")
    install(second)
    second.execute(
        "CREATE TEMP TABLE norm_input AS SELECT 'chrS' AS chrom, 2 AS pos, "
        "'T' AS ref, 'C' AS alt"
    )
    fasta = ROOT / "test/data/liftover_repeat_src.fa"
    assert second.execute(
        f"SELECT count(*) FROM duckhts_bcftools_norm('norm_input', '{fasta}')"
    ).fetchone()[0] == 1
    assert second.execute(
        "WITH columns AS (SELECT column_name FROM inputs) "
        "SELECT duckhts_quote_ident(column_name) FROM columns"
    ).fetchone()[0] == '"a"'
    first.close()
    second.close()


def reopened(path):
    con = connect(path, read_only=True)
    load(con)
    rows = definitions(con)
    catalog = con.execute("SELECT current_database()").fetchone()[0]
    assert names(con, catalog) & {row[0] for row in rows} == {"duckhts_quote_ident"}
    assert names(con, "temp") & {row[0] for row in rows} == set()
    assert con.execute("SELECT duckhts_htslib_version()").fetchone()[0]
    install(con)
    con.execute("CREATE TEMP TABLE norm_input AS "
                "SELECT 'chrS' AS chrom, 2 AS pos, 'T' AS ref, 'C' AS alt")
    fasta = ROOT / "test/data/liftover_repeat_src.fa"
    assert con.execute(
        f"SELECT count(*) FROM duckhts_bcftools_norm('norm_input', '{fasta}')"
    ).fetchone()[0] == 1
    con.close()


def attached(path):
    con = duckdb.connect(":memory:", config=CONFIG)
    load(con)
    quoted = str(path).replace("'", "''")
    con.execute(f"ATTACH '{quoted}' AS attached (READ_ONLY)")
    assert "duckhts_quote_ident" in names(con, "memory")
    assert names(con, "attached") & {row[0] for row in definitions(con)} == {
        "duckhts_quote_ident"
    }
    con.close()


def rollback(path):
    con = connect(path)
    load(con)
    public = {name for name, exposed, _, _ in definitions(con) if exposed}
    con.execute("BEGIN TRANSACTION")
    install(con)
    assert names(con, "temp") & public == public
    con.execute("ROLLBACK")
    assert names(con, "temp") & public == set()
    install(con)
    assert names(con, "temp") & public == public
    con.close()


def main():
    with tempfile.TemporaryDirectory(prefix="duckhts_macro_catalog_") as tmp:
        clean = pathlib.Path(tmp) / "clean.duckdb"
        subprocess.run([sys.executable, __file__, "clean_file", str(clean)], check=True)
        subprocess.run([sys.executable, __file__, "clean_reopened", str(clean)], check=True)
        path = pathlib.Path(tmp) / "catalog.duckdb"
        subprocess.run([sys.executable, __file__, "writable", str(path)], check=True)
        before = digest(path)
        subprocess.run([sys.executable, __file__, "reopened", str(path)], check=True)
        assert digest(path) == before
        subprocess.run([sys.executable, __file__, "attached", str(path)], check=True)
        assert digest(path) == before
        subprocess.run([sys.executable, __file__, "rollback", str(path)], check=True)
        print("macro catalog: 36 definitions, 27 public, file/readonly/rollback/two connections OK")


if __name__ == "__main__":
    actions = {"clean_file": clean_file, "clean_reopened": clean_reopened,
               "writable": writable, "reopened": reopened,
               "attached": attached, "rollback": rollback}
    if len(sys.argv) == 3 and sys.argv[1] in actions:
        actions[sys.argv[1]](pathlib.Path(sys.argv[2]))
    else:
        main()
