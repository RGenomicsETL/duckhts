"""read_bam through the DuckDB file system (SET duckhts_file_system = 'duckdb').

htslib's own file access is the oracle: every committed BAM and CRAM fixture
must give the same rows through both settings. The remaining tests cover what
only the DuckDB route does: DuckDB's access policy, its error text, index
lookup that ignores the working directory, and HTTP secrets.
"""

import argparse
import contextlib
import os
import shutil
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import duckdb

DATA = Path(__file__).resolve().parents[1] / "data"
INDEXED = DATA / "bam_scan_mixed.bam"
ROW_DIGEST = (
    "SELECT count(*), sum(hash(QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL))"
    " FROM read_bam(?{options})"
)


def connect(file_system=None):
    connection = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    connection.execute("LOAD '" + str(args.extension).replace("'", "''") + "'")
    if file_system:
        connection.execute("SET duckhts_file_system = '" + file_system + "'")
    return connection


def outcome(connection, path, options=""):
    """The row digest of one read, or the class of the error it raises."""
    try:
        return connection.execute(ROW_DIGEST.format(options=options), [str(path)]).fetchone()
    except duckdb.Error as failure:
        return type(failure).__name__


def contigs(connection, path):
    rows = connection.execute(
        "SELECT DISTINCT RNAME FROM read_bam(?) WHERE RNAME IS NOT NULL ORDER BY 1", [str(path)]
    ).fetchall()
    return [row[0] for row in rows]


class RangeHandler(BaseHTTPRequestHandler):
    """Serves DATA with byte ranges; paths under /protected/ need 'Bearer sesame'."""

    protocol_version = "HTTP/1.1"

    def log_message(self, *arguments):
        pass

    def reply(self, status, body=b"", content_range=None):
        self.send_response(status)
        self.send_header("Accept-Ranges", "bytes")
        if content_range:
            self.send_header("Content-Range", content_range)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def serve(self):
        name = self.path
        if name.startswith("/protected/"):
            if self.headers.get("Authorization") != "Bearer sesame":
                self.reply(403)
                return
            name = name[len("/protected") :]
        target = DATA / name.lstrip("/")
        if "/" in name.lstrip("/") or not target.is_file():
            self.reply(404)
            return
        content = target.read_bytes()
        requested = self.headers.get("Range")
        if not requested:
            self.reply(200, content)
            return
        first_text, last_text = requested.split("=", 1)[1].split("-", 1)
        first = int(first_text)
        last = min(int(last_text), len(content) - 1) if last_text else len(content) - 1
        if first >= len(content):
            self.reply(416, content_range=f"bytes */{len(content)}")
            return
        self.reply(206, content[first : last + 1], f"bytes {first}-{last}/{len(content)}")

    do_GET = serve
    do_HEAD = serve


class DuckdbFileSystemRoute(unittest.TestCase):
    def test_default_is_htslib(self):
        with connect() as connection:
            value = connection.execute("SELECT current_setting('duckhts_file_system')").fetchone()[0]
        self.assertEqual(value, "htslib")

    def test_fixtures_read_the_same_through_both_file_systems(self):
        fixtures = sorted(DATA.glob("*.bam")) + sorted(DATA.glob("*.cram"))
        self.assertGreaterEqual(len(fixtures), 30)
        compared_rows = 0
        compared_regions = 0
        with connect("htslib") as oracle, connect("duckdb") as routed:
            for fixture in fixtures:
                with self.subTest(fixture=fixture.name, scan="full"):
                    expected = outcome(oracle, fixture)
                    self.assertEqual(outcome(routed, fixture), expected)
                    if isinstance(expected, tuple):
                        compared_rows += expected[0]
                if not isinstance(expected, tuple):
                    continue
                for contig in contigs(oracle, fixture):
                    options = ", region := '" + contig.replace("'", "''") + "'"
                    with self.subTest(fixture=fixture.name, region=contig):
                        expected_region = outcome(oracle, fixture, options)
                        self.assertEqual(outcome(routed, fixture, options), expected_region)
                        if isinstance(expected_region, tuple):
                            compared_regions += 1
        self.assertGreater(compared_rows, 1000)
        self.assertGreater(compared_regions, 10)

    def test_invalid_setting_is_a_bind_error(self):
        with connect("tabix") as connection:
            with self.assertRaisesRegex(duckdb.Error, "duckhts_file_system must be 'htslib' or 'duckdb'"):
                connection.execute(ROW_DIGEST.format(options=""), [str(INDEXED)])

    def test_missing_file_reports_the_path_and_duckdb_reason(self):
        missing = str(INDEXED) + ".missing"
        with connect("duckdb") as connection:
            with self.assertRaises(duckdb.Error) as failure:
                connection.execute(ROW_DIGEST.format(options=""), [missing])
        message = str(failure.exception)
        self.assertIn("read_bam: failed to open SAM/BAM/CRAM file: " + missing, message)
        self.assertIn("No such file or directory", message)
        self.assertNotIn("duckdbfs:", message)

    def test_duckdb_access_policy_applies_to_routed_reads(self):
        with connect("duckdb") as connection:
            expected = outcome(connection, INDEXED)
            self.assertIsInstance(expected, tuple)
            connection.execute("SET enable_external_access = false")
            with self.assertRaisesRegex(duckdb.Error, "(?i)disabled|permission"):
                connection.execute(ROW_DIGEST.format(options=""), [str(INDEXED)])

    def test_index_lookup_ignores_the_working_directory(self):
        with connect("htslib") as oracle:
            contig = contigs(oracle, INDEXED)[0]
            options = ", region := '" + contig + "'"
            expected = outcome(oracle, INDEXED, options)
        self.assertIsInstance(expected, tuple)
        with tempfile.TemporaryDirectory() as directory, contextlib.chdir(directory):
            for suffix in (".bai", ".csi"):
                Path(INDEXED.name + suffix).write_bytes(b"not an index")
            with connect("duckdb") as routed:
                self.assertEqual(outcome(routed, INDEXED, options), expected)
            self.assertEqual(sorted(os.listdir(".")), [INDEXED.name + ".bai", INDEXED.name + ".csi"])

    def test_index_forms(self):
        with connect("htslib") as oracle:
            contig = contigs(oracle, INDEXED)[0]
            options = ", region := '" + contig + "'"
            expected = outcome(oracle, INDEXED, options)
        bai = DATA / (INDEXED.name + ".bai")
        csi = DATA / (INDEXED.name + ".csi")
        with tempfile.TemporaryDirectory() as directory, connect("duckdb") as routed:
            directory = Path(directory)
            copy = directory / "reads.bam"
            shutil.copy(INDEXED, copy)
            with self.subTest(index="none"):
                with self.assertRaisesRegex(duckdb.Error, "Region query requires an index"):
                    routed.execute(ROW_DIGEST.format(options=options), [str(copy)])
            with self.subTest(index="explicit index_path"):
                explicit = options + ", index_path := '" + str(bai) + "'"
                self.assertEqual(outcome(routed, copy, explicit), expected)
            for name, source in (("reads.bai", bai), ("reads.csi", csi), ("reads.bam.bai", bai), ("reads.bam.csi", csi)):
                with self.subTest(index=name):
                    shutil.copy(source, directory / name)
                    self.assertEqual(outcome(routed, copy, options), expected)
                    (directory / name).unlink()

    def test_remote_block_size_default_and_bounds(self):
        with connect() as connection:
            default = connection.execute("SELECT current_setting('duckhts_remote_block_size')").fetchone()[0]
            self.assertEqual(default, 0)
            for invalid in (-1, 1, 131071, 268435457):
                with self.subTest(invalid=invalid):
                    connection.execute(f"SET duckhts_remote_block_size = {invalid}")
                    with self.assertRaisesRegex(
                        duckdb.Error,
                        "read_bam: duckhts_remote_block_size must be 0 or between 131072 and 268435456 bytes",
                    ):
                        connection.execute(ROW_DIGEST.format(options=""), [str(INDEXED)])
            for valid in (0, 131072, 268435456):
                with self.subTest(valid=valid):
                    connection.execute(f"SET duckhts_remote_block_size = {valid}")
                    self.assertIsInstance(outcome(connection, INDEXED), tuple)

    def test_http_block_sizes_and_secret(self):
        with connect("duckdb") as routed:
            try:
                routed.execute("LOAD httpfs")
            except duckdb.Error as failure:
                self.skipTest("httpfs is not installed for this DuckDB: " + str(failure).splitlines()[0])
        server = ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            base = "http://127.0.0.1:" + str(server.server_address[1])
            with connect("htslib") as oracle:
                contig = contigs(oracle, INDEXED)[0]
                options = ", region := '" + contig + "'"
                expected = outcome(oracle, INDEXED), outcome(oracle, INDEXED, options)
            with connect("duckdb") as routed:
                routed.execute("LOAD httpfs")
                open_url = base + "/" + INDEXED.name
                self.assertEqual((outcome(routed, open_url), outcome(routed, open_url, options)), expected)
                # A block size changes how the bytes are fetched, never the rows.
                for file_system in ("htslib", "duckdb"):
                    routed.execute("SET duckhts_file_system = '" + file_system + "'")
                    for block_size in (131072, 4194304):
                        with self.subTest(file_system=file_system, block_size=block_size):
                            routed.execute(f"SET duckhts_remote_block_size = {block_size}")
                            self.assertEqual(
                                (outcome(routed, open_url), outcome(routed, open_url, options)), expected
                            )
                routed.execute("RESET duckhts_remote_block_size")
                protected = base + "/protected/" + INDEXED.name
                with self.assertRaisesRegex(duckdb.Error, "403"):
                    routed.execute(ROW_DIGEST.format(options=""), [protected])
                routed.execute(
                    "CREATE SECRET (TYPE http, BEARER_TOKEN 'sesame', SCOPE '" + base + "/protected')"
                )
                self.assertEqual((outcome(routed, protected), outcome(routed, protected, options)), expected)
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    args, remaining = parser.parse_known_args()
    args.extension = args.extension.resolve()
    unittest.main(argv=[__file__] + remaining)
