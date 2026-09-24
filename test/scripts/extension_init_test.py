"""Exercise extension initialization through pinned native DuckDB runtimes."""

import argparse
from pathlib import Path
import sys
import unittest

import duckdb


class ExtensionInitialization(unittest.TestCase):
    def test_runtime_admission_and_macro_dependencies(self):
        with duckdb.connect(config={"allow_unsigned_extensions": "true"}) as connection:
            load = "LOAD '" + str(args.extension).replace("'", "''") + "'"
            if args.expect == "unsupported":
                with self.assertRaisesRegex(
                    duckdb.Error, "DuckHTS requires DuckDB 1.4.0 or newer"
                ) as failure:
                    connection.execute(load)
                self.assertNotIsInstance(failure.exception, duckdb.FatalException)
                registered = connection.execute(
                    "SELECT count(*) FROM duckdb_functions() WHERE function_name IN "
                    "('read_bam', 'read_hts_header', 'duckhts_somalier_bam_counts')"
                ).fetchone()
                self.assertEqual(registered, (0,))
            else:
                connection.execute(load)
                sites = Path(__file__).resolve().parents[1] / "data" / "somalier_sites.vcf"
                self.assertEqual(
                    connection.execute(
                        "SELECT count(*) FROM duckhts_somalier_import_sites(?, 'GRCh38')",
                        [str(sites)],
                    ).fetchone(),
                    (3,),
                )
                self.assertEqual(
                    connection.execute(
                        "SELECT cigar_query_length('5M2I'), duckhts_quote_ident('two words')"
                    ).fetchone(),
                    (7, '"two words"'),
                )
            self.assertEqual(connection.execute("SELECT 42").fetchone(), (42,))

    def test_sql_registration_failure_preserves_cause_and_connection(self):
        with duckdb.connect(config={"allow_unsigned_extensions": "true"}) as connection:
            load = "LOAD '" + str(args.failure_extension).replace("'", "''") + "'"
            with self.assertRaisesRegex(
                duckdb.Error, "__duckhts_missing_init_relation"
            ) as failure:
                connection.execute(load)
            self.assertNotIsInstance(failure.exception, duckdb.FatalException)
            self.assertEqual(connection.execute("SELECT 42").fetchone(), (42,))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", type=Path, required=True)
    parser.add_argument("--failure-extension", type=Path, required=True)
    parser.add_argument("--expect", choices=("supported", "unsupported"), default="supported")
    args = parser.parse_args()
    args.extension = args.extension.resolve(strict=True)
    args.failure_extension = args.failure_extension.resolve(strict=True)
    print("DuckDB", duckdb.__version__, "expected", args.expect, flush=True)
    unittest.main(argv=[sys.argv[0]], verbosity=2)
