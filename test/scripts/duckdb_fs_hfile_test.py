"""DuckHTS readers through the DuckDB file system (SET duckhts_file_system = 'duckdb').

htslib's own file access is the oracle: every reader must give the same rows
for the committed fixtures through both settings, on full scans and on index
or region paths, and must fail on the same inputs. The remaining tests cover
what only the DuckDB route does: DuckDB's access policy, its error text, index
lookup that ignores the working directory, a read-only route that builds no
missing index, and HTTP secrets, including the reference of a CRAM. The scalar
functions that read a FASTA row by row, bcftools_norm_row and
bcftools_munge_row, are compared row for row over variants cut from the
fixture FASTAs, and their handles are shown to stay with the query that
opened them.
"""

import argparse
import collections
import contextlib
import csv
import gzip
import hashlib
import os
import re
import shutil
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "test" / "data"
BIGWIG = ROOT / "third_party" / "libBigWig" / "test" / "test.bw"
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


# DuckDB's reason for a failed open, which only the DuckDB route appends to the reader's message.
DUCKDB_REASON = re.compile(r' \((?:Cannot open file|HTTP|IO Error|Permission Error)[^\n]*\)$')


def relation_digest(connection, relation, parameters=()):
    """Row count and sum of row hashes of a relation, or the class and message of its error.

    The message is the reader's own first line: DuckDB's reason for a failed
    open is removed, so both routes must report the same failure.
    """
    try:
        return connection.execute(
            "SELECT count(*), sum(hash(t)) FROM (" + relation + ") t", list(parameters)
        ).fetchone()
    except duckdb.Error as failure:
        message = str(failure).splitlines()[0]
        return type(failure).__name__ + ": " + DUCKDB_REASON.sub("", message)


def first_value(connection, query, path):
    """The first value of a one-column query over path, or None when it fails or is empty."""
    try:
        row = connection.execute(query, [str(path)]).fetchone()
    except duckdb.Error:
        return None
    return row[0] if row and row[0] is not None else None


def quoted(text):
    return "'" + str(text).replace("'", "''") + "'"


def fixtures(*patterns, directory=DATA):
    found = []
    for pattern in patterns:
        found.extend(directory.glob(pattern))
    return sorted(set(found))


def file_digests(directory):
    return {
        entry.name: hashlib.sha256(entry.read_bytes()).hexdigest()
        for entry in sorted(Path(directory).iterdir())
        if entry.is_file()
    }


def contigs(connection, path):
    rows = connection.execute(
        "SELECT DISTINCT RNAME FROM read_bam(?) WHERE RNAME IS NOT NULL ORDER BY 1", [str(path)]
    ).fetchall()
    return [row[0] for row in rows]


# One row of each scalar function that reads a FASTA; the one ? is the FASTA path.
NORM_ONE_ROW = "SELECT bcftools_norm_row('chrF', 2, 'C', 'CT', ?, NULL::BIGINT, NULL::BIGINT, NULL, NULL)"
MUNGE_ONE_ROW = (
    "SELECT bcftools_munge_row('chrF', 2, 'A', 'C', 'rs1', 0.01, NULL, NULL, 0.12, 1000, NULL, NULL, NULL, 0.3,"
    " NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, '+-', ?, NULL, NULL, NULL, NULL, NULL)"
)


def norm_cases(where="", fasta="fasta", fai="NULL", gzi="NULL", alt="alt", cases="norm_cases"):
    """bcftools_norm_row over the cases load_reference_cases creates; fasta, fai, gzi and alt are SQL."""
    return (
        "SELECT id, kind, bcftools_norm_row(chrom, pos, ref, " + alt + ", " + fasta + ", end_pos, svlen, "
        + fai + ", " + gzi + ") AS normed FROM " + cases + where + " ORDER BY id"
    )


def munge_cases(where="", fasta="fasta", cases="munge_cases"):
    """bcftools_munge_row over the cases load_reference_cases creates; fasta is SQL."""
    return (
        "SELECT id, kind, bcftools_munge_row(chrom, pos, a1, a2, 'rs' || id, 0.01, NULL, NULL, 0.12, 1000, NULL,"
        " NULL, NULL, 0.3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, '+-', " + fasta + ", NULL, NULL, NULL,"
        " NULL, NULL) AS munged FROM " + cases + where + " ORDER BY id"
    )


REFERENCE_WINDOW_BASES = 65536  # DUCKHTS_REFERENCE_WINDOW_BASES in src/include/reference_cache.h


def fasta_sequences(path):
    """The sequences of a FASTA read by Python alone, as (name, bases) in file order."""
    opener = gzip.open if path.suffix == ".gz" else open
    sequences = []
    with opener(path, "rt") as lines:
        for line in lines:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                sequences.append([line[1:].split()[0] if line[1:].split() else "", []])
            elif sequences:
                sequences[-1][1].append(line)
    return [(name, "".join(parts)) for name, parts in sequences if name]


def sampled_positions(length):
    """1-based positions: all of a short sequence; the ends and each cache-window edge of a long one."""
    if length <= 2000:
        return list(range(1, length + 1))
    positions = set(range(1, 301)) | set(range(length - 299, length + 1))
    for edge in range(REFERENCE_WINDOW_BASES, length, REFERENCE_WINDOW_BASES):
        if edge > 4 * REFERENCE_WINDOW_BASES:
            break
        positions |= set(range(edge - 150, edge + 150))
    return sorted(positions)


def other_base(base, skip=""):
    return next(candidate for candidate in "ACGT" if candidate != base.upper() and candidate not in skip)


def reference_cases(fasta):
    """Variants and summary-statistic rows cut from the bases of one FASTA.

    The alleles come from the file as Python reads it, so they do not depend on
    either file-access route. Deletions and insertions inside a repeat
    left-align; one tandem pair is multi-allelic.
    """
    norm_rows = []
    munge_rows = []
    for name, bases in fasta_sequences(fasta):
        for position in sampled_positions(len(bases)):
            start = position - 1
            base = bases[start]
            three = bases[start : start + 3]
            other = other_base(base)
            norm_rows.append(("snp", name, position, base, other, None, None))
            norm_rows.append(("ref_mismatch", name, position, other, base.upper(), None, None))
            norm_rows.append(("symbolic_del", name, position, base, "<DEL>", position + 4, None))
            norm_rows.append(("symbolic_dup", name, position, base, "<DUP>", None, 3))
            if len(bases) - start >= 2:
                norm_rows.append(("deletion_shared_tail", name, position, bases[start : start + 2], bases[start + 1],
                                  None, None))
            if len(three) == 3:
                norm_rows.append(("deletion", name, position, three, base, None, None))
                norm_rows.append(("insertion", name, position, base, three, None, None))
                norm_rows.append(("multi_allelic", name, position, three, base + "," + three + three[1:], None, None))
            munge_rows.append(("a2_is_reference", name, position, other, base))
            munge_rows.append(("a1_is_reference", name, position, base, other))
            munge_rows.append(("neither_is_reference", name, position, other, other_base(base, skip=other)))
            if len(bases) - start >= 2:
                munge_rows.append(("both_match", name, position, bases[start : start + 2], base))
        alias = name[3:] if name.lower().startswith("chr") else "chr" + name
        if alias and bases:
            norm_rows.append(("contig_alias", alias, 1, bases[0], other_base(bases[0]), None, None))
            munge_rows.append(("contig_alias", alias, 1, other_base(bases[0]), bases[0]))
    norm_rows.append(("missing_contig", "no_such_contig", 1, "A", "C", None, None))
    munge_rows.append(("missing_contig", "no_such_contig", 1, "A", "C"))
    return norm_rows, munge_rows


def load_reference_cases(connections, directory, cases_by_fasta):
    """Creates norm_cases and munge_cases, the same rows in every connection.

    cases_by_fasta maps the FASTA path each row names to reference_cases() of
    the file its bases were cut from. Returns the number of rows of each table.
    """
    norm_file = Path(directory) / "norm_cases.csv"
    munge_file = Path(directory) / "munge_cases.csv"
    counts = []
    for target, table, width in ((norm_file, 0, 7), (munge_file, 1, 5)):
        with open(target, "w", newline="") as handle:
            writer = csv.writer(handle)
            identifier = 0
            for fasta, cases in cases_by_fasta.items():
                for row in cases[table]:
                    identifier += 1
                    assert len(row) == width
                    writer.writerow([identifier, str(fasta)] + ["" if value is None else value for value in row])
            counts.append(identifier)
    for connection in connections:
        connection.execute(
            "CREATE OR REPLACE TEMP TABLE norm_cases AS SELECT * FROM read_csv(?, header = false, quote = '\"',"
            " nullstr = '', columns = {'id': 'BIGINT', 'fasta': 'VARCHAR', 'kind': 'VARCHAR', 'chrom': 'VARCHAR',"
            " 'pos': 'BIGINT', 'ref': 'VARCHAR', 'alt': 'VARCHAR', 'end_pos': 'BIGINT', 'svlen': 'BIGINT'})",
            [str(norm_file)],
        )
        connection.execute(
            "CREATE OR REPLACE TEMP TABLE munge_cases AS SELECT * FROM read_csv(?, header = false, quote = '\"',"
            " nullstr = '', columns = {'id': 'BIGINT', 'fasta': 'VARCHAR', 'kind': 'VARCHAR', 'chrom': 'VARCHAR',"
            " 'pos': 'BIGINT', 'a1': 'VARCHAR', 'a2': 'VARCHAR'})",
            [str(munge_file)],
        )
    return tuple(counts)


class RangeHandler(BaseHTTPRequestHandler):
    """Serves DATA with byte ranges.

    Paths under /protected/ need 'Bearer sesame'. Paths under /signed/ stand for
    presigned URLs: each needs the query string '?sig=ok' and is refused with 403
    without it.
    """

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

    requested_paths = collections.Counter()

    def serve(self):
        name = self.path
        RangeHandler.requested_paths[name] += 1
        if name.startswith("/protected/"):
            if self.headers.get("Authorization") != "Bearer sesame":
                self.reply(403)
                return
            name = name[len("/protected") :]
        if name.startswith("/signed/"):
            name, _, query = name[len("/signed") :].partition("?")
            if query != "sig=ok":
                self.reply(403)
                return
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

    def test_index_of_a_url_with_a_query_string(self):
        """Through the DuckDB route a presigned URL finds its index with the suffix before the query.

        The oracle is htslib reading the local file. htslib's own route is not
        compared: it fetches the .csi with the query string and then asks for it
        again without one, so a region query of a presigned URL fails there.
        """
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
            signed = base + "/signed/" + INDEXED.name + "?sig=ok"
            with connect("htslib") as oracle:
                contig = contigs(oracle, INDEXED)[0]
                options = ", region := '" + contig + "'"
                expected = outcome(oracle, INDEXED, options)
            self.assertIsInstance(expected, tuple)
            with connect("duckdb") as connection:
                connection.execute("LOAD httpfs")
                RangeHandler.requested_paths.clear()
                self.assertEqual(outcome(connection, signed, options), expected)
                index_requests = [path for path in RangeHandler.requested_paths
                                  if path.startswith("/signed/" + INDEXED.name + ".")]
                # htslib's order puts .csi first; the fixture has one.
                self.assertEqual(index_requests, ["/signed/" + INDEXED.name + ".csi?sig=ok"])
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


class ReadersThroughBothFileSystems(unittest.TestCase):
    """Every reader gives htslib's rows, or htslib's failure, through the DuckDB file system."""

    compared = collections.Counter()

    @classmethod
    def setUpClass(cls):
        cls.oracle = connect("htslib")
        cls.routed = connect("duckdb")
        cls.scratch = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.oracle.close()
        cls.routed.close()
        cls.scratch.cleanup()
        print("\nrelations equal through both file systems (reader: rows-equal + errors-equal):")
        for reader in sorted({name for name, _ in cls.compared}):
            print(f"  {reader}: {cls.compared[(reader, 'rows')]} + {cls.compared[(reader, 'error')]}")

    def compare(self, reader, relation, parameters=(), label=None):
        """Asserts one relation is the same through both routes; returns htslib's outcome."""
        with self.subTest(reader=reader, case=label or relation):
            expected = relation_digest(self.oracle, relation, parameters)
            observed = relation_digest(self.routed, relation, parameters)
            if observed != expected and isinstance(observed, str) and isinstance(expected, str):
                # A failing parallel scan reports the worker that fails first, which
                # varies between runs of one route. One thread makes it one message.
                for connection in (self.oracle, self.routed):
                    connection.execute("SET threads = 1")
                expected = relation_digest(self.oracle, relation, parameters)
                observed = relation_digest(self.routed, relation, parameters)
                for connection in (self.oracle, self.routed):
                    connection.execute("RESET threads")
            self.assertEqual(observed, expected)
            self.compared[(reader, "rows" if isinstance(expected, tuple) else "error")] += 1
        return expected

    def compare_over(self, reader, paths, relation, region_query=None, region_option="region"):
        """relation has one ? for the path; region_query finds a region with htslib."""
        for path in paths:
            self.compare(reader, relation.format(options=""), [str(path)], path.name)
            if region_query is None:
                continue
            region = first_value(self.oracle, region_query, path)
            if region is None:
                continue
            options = ", " + region_option + " := " + quoted(region)
            self.compare(reader + " region", relation.format(options=options), [str(path)], path.name)

    def assert_readers_compared(self, minimums):
        for reader, minimum in minimums.items():
            with self.subTest(reader=reader):
                self.assertGreaterEqual(self.compared[(reader, "rows")], minimum)

    def test_alignment_readers(self):
        alignments = fixtures("*.bam", "*.cram", "*.sam")
        insert_sizes = fixtures("*.bam", "*.cram", directory=DATA / "insert_size")
        self.assertGreaterEqual(len(alignments), 40)
        first_contig = "SELECT RNAME FROM read_bam(?) WHERE RNAME IS NOT NULL LIMIT 1"
        self.compare_over("bam_bin_counts", alignments, "SELECT * FROM bam_bin_counts(?, 100{options})")
        self.compare_over("read_pileup", alignments, "SELECT * FROM read_pileup(?{options})", first_contig)
        self.compare_over(
            "duckhts_bam_insert_size_counts",
            alignments + insert_sizes,
            "SELECT * FROM duckhts_bam_insert_size_counts(?{options})",
        )
        self.compare_over(
            "duckhts_bam_insert_size_counts sequential",
            insert_sizes,
            "SELECT * FROM duckhts_bam_insert_size_counts(?, scan_mode := 'sequential'{options})",
        )
        self.compare_over(
            "duckhts_bam_mate_overlap_counts",
            alignments,
            "SELECT * FROM duckhts_bam_mate_overlap_counts(?{options})",
            first_contig,
        )
        for name in ("read_hts_header", "read_hts_index", "read_hts_index_spans"):
            self.compare_over(name, alignments, "SELECT * FROM " + name + "(?{options})")
        self.assert_readers_compared(
            {
                "bam_bin_counts": 20,
                "read_pileup region": 10,
                "duckhts_bam_insert_size_counts": 20,
                "duckhts_bam_mate_overlap_counts": 10,
                "read_hts_header": 30,
                "read_hts_index": 15,
            }
        )

    def test_alignment_readers_with_a_reference(self):
        reference = quoted(DATA / "fixture_ref.fa")
        worm = quoted(DATA / "ce.fa")
        mismatch = quoted(DATA / "bam_mismatch.fa")
        bed = quoted(DATA / "fixture_mixed_regions.bed")
        for alignment, fasta in (("fixture_mixed.cram", reference), ("range.cram", worm)):
            path = [str(DATA / alignment)]
            self.compare("read_bam reference", "SELECT * FROM read_bam(?, reference := " + fasta + ")", path)
            self.compare(
                "bam_bin_counts reference", "SELECT * FROM bam_bin_counts(?, 100, reference := " + fasta + ")", path
            )
            self.compare(
                "duckhts_bam_insert_size_counts reference",
                "SELECT * FROM duckhts_bam_insert_size_counts(?, reference := " + fasta + ")",
                path,
            )
            self.compare(
                "duckhts_bam_mate_overlap_counts reference",
                "SELECT * FROM duckhts_bam_mate_overlap_counts(?, reference := " + fasta + ")",
                path,
            )
        for alignment in ("fixture_mixed.bam", "fixture_mixed.cram"):
            self.compare(
                "duckhts_bam_bed_coverage",
                "SELECT * FROM duckhts_bam_bed_coverage(?, " + bed + ", reference := " + reference + ")",
                [str(DATA / alignment)],
            )
        for alignment in ("bam_mismatch.bam", "bam_mismatch.cram"):
            for options in ("", ", region := 'ref1'", ", region := 'ref2:1-40'"):
                self.compare(
                    "duckhts_bam_mismatch_counts",
                    "SELECT * FROM duckhts_bam_mismatch_counts(?, " + mismatch + options + ")",
                    [str(DATA / alignment)],
                    alignment + options,
                )
        for mask in fixtures("bam_mismatch*.vcf.gz", "bam_mismatch*.bcf"):
            self.compare(
                "duckhts_bam_mismatch_counts mask",
                "SELECT * FROM duckhts_bam_mismatch_counts(?, " + mismatch + ", mask := " + quoted(mask) + ")",
                [str(DATA / "bam_mismatch.bam")],
                mask.name,
            )
        panel = Path(self.scratch.name) / "somalier_panel.parquet"
        self.oracle.execute(
            "COPY (SELECT * FROM (VALUES"
            " ('WBcel235', 0, 'CHROMOSOME_I', 1, 'A', 'G'), ('WBcel235', 1, 'CHROMOSOME_I', 914, 'A', 'C'),"
            " ('WBcel235', 2, 'CHROMOSOME_I', 2, 'A', 'G'), ('WBcel235', 3, 'NOT_IN_REFERENCE', 1, 'A', 'C'))"
            " panel(assembly, site_index, region, position, allele_a, allele_b)) TO "
            + quoted(panel)
            + " (FORMAT parquet)"
        )
        for alignment in ("range.bam", "range.cram"):
            self.compare(
                "duckhts_somalier_bam_counts",
                "SELECT * FROM duckhts_somalier_bam_counts(?, NULL, 'sample-1', "
                + worm
                + ", panel_parquet := "
                + quoted(panel)
                + ")",
                [str(DATA / alignment)],
            )
        self.assert_readers_compared(
            {
                "read_bam reference": 2,
                "bam_bin_counts reference": 2,
                "duckhts_bam_bed_coverage": 2,
                "duckhts_bam_mismatch_counts": 6,
                "duckhts_somalier_bam_counts": 2,
            }
        )

    def test_alignment_writers_read_their_input_through_both_file_systems(self):
        """idxstats and mosdepth write with this process; their input follows the setting."""
        scratch = Path(self.scratch.name)
        written = {}
        for route, connection in (("htslib", self.oracle), ("duckdb", self.routed)):
            directory = scratch / ("written_" + route)
            directory.mkdir()
            rows = []
            for alignment in fixtures("range.bam", "range.cram", "bam_scan_mixed.bam", "fixture_mixed.bam"):
                rows.append(
                    connection.execute(
                        "SELECT success, used_index_fast_path FROM duckhts_samtools_idxstats(?, output := ?)",
                        [str(alignment), str(directory / (alignment.name + ".idxstats"))],
                    ).fetchall()
                )
            for alignment in fixtures("range.bam", "bam_scan_mixed.bam"):
                rows.append(
                    connection.execute(
                        "SELECT success FROM duckhts_mosdepth(?, ?, by := '500')",
                        [str(directory / (alignment.name + ".mosdepth")), str(alignment)],
                    ).fetchall()
                )
            written[route] = (rows, file_digests(directory))
        self.assertGreaterEqual(len(written["htslib"][1]), 10)
        self.assertEqual(written["duckdb"], written["htslib"])
        self.compared[("duckhts_samtools_idxstats", "rows")] += 4
        self.compared[("duckhts_mosdepth", "rows")] += 2

    def test_variant_readers(self):
        variants = fixtures("*.vcf", "*.vcf.gz", "*.bcf")
        self.assertGreaterEqual(len(variants), 100)
        first_contig = "SELECT CHROM FROM read_bcf(?) LIMIT 1"
        self.compare_over("read_bcf", variants, "SELECT * FROM read_bcf(?{options})", first_contig)
        self.compare_over("read_geno", variants, "SELECT * FROM read_geno(?{options})", first_contig)
        self.compare_over("read_bcf_samples", variants, "SELECT * FROM read_bcf_samples(?{options})")
        for name in ("read_hts_header", "read_hts_index", "read_hts_index_spans"):
            self.compare_over(name + " variants", variants, "SELECT * FROM " + name + "(?{options})")
        self.assert_readers_compared(
            {
                "read_bcf": 80,
                "read_bcf region": 30,
                "read_geno": 60,
                "read_geno region": 20,
                "read_bcf_samples": 80,
                "read_hts_index variants": 30,
            }
        )

    def test_tabix_and_interval_readers(self):
        tables = fixtures("*.tsv", "*.txt", "*.bed", "*.bed.gz", "*.tsv.gz", "*.vcf.gz", "*.gff3.gz", "*.gtf.gz")
        first_sequence = "SELECT column0 FROM read_tabix(?) LIMIT 1"
        self.compare_over("read_tabix", tables, "SELECT * FROM read_tabix(?{options})", first_sequence)
        self.compare_over(
            "read_gff",
            fixtures("*.gff3", "*.gff3.gz"),
            "SELECT * FROM read_gff(?{options})",
            "SELECT seqname FROM read_gff(?) LIMIT 1",
        )
        self.compare_over(
            "read_gtf",
            fixtures("*.gtf", "*.gtf.gz"),
            "SELECT * FROM read_gtf(?{options})",
            "SELECT seqname FROM read_gtf(?) LIMIT 1",
        )
        self.compare_over(
            "read_bed",
            fixtures("*.bed", "*.bed.gz"),
            "SELECT * FROM read_bed(?{options})",
            "SELECT chrom FROM read_bed(?) LIMIT 1",
        )
        indexed_bed = Path(self.scratch.name) / "regions.bed.gz"
        self.oracle.execute(
            "SELECT * FROM bgzip(?, output_path := ?)", [str(DATA / "fixture_mixed_regions.bed"), str(indexed_bed)]
        ).fetchall()
        self.oracle.execute("SELECT * FROM tabix_index(?, preset := 'bed')", [str(indexed_bed)]).fetchall()
        self.compare_over(
            "read_bed", [indexed_bed], "SELECT * FROM read_bed(?{options})", "SELECT chrom FROM read_bed(?) LIMIT 1"
        )
        self.compare(
            "fasta_nuc bed region",
            "SELECT * FROM fasta_nuc(?, bed_path := ?, region := '11')",
            [str(DATA / "fixture_ref.fa"), str(indexed_bed)],
        )
        self.compare_over(
            "read_hts_index tabix", fixtures("*.bed.gz", "*.gff3.gz", "*.gtf.gz"), "SELECT * FROM read_hts_index(?{options})"
        )
        self.assert_readers_compared(
            {"read_tabix": 60, "read_tabix region": 10, "read_gff": 8, "read_gff region": 1,
             "read_gtf": 2, "read_gtf region": 1, "read_bed": 10, "read_bed region": 1,
             "fasta_nuc bed region": 1}
        )

    def test_sequence_readers(self):
        fastas = fixtures("*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz")
        fastqs = fixtures("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz")
        self.assertGreaterEqual(len(fastas), 25)
        first_sequence = "SELECT NAME FROM read_fasta(?) LIMIT 1"
        self.compare_over("read_fasta", fastas, "SELECT * FROM read_fasta(?{options})", first_sequence)
        self.compare_over("read_fasta count", fastas, "SELECT count(*) FROM read_fasta(?{options})")
        self.compare_over("read_fastq", fastqs, "SELECT * FROM read_fastq(?{options})")
        self.compare_over("detect_quality_encoding", fastqs, "SELECT * FROM detect_quality_encoding(?{options})")
        self.compare_over("fasta_nuc", fastas, "SELECT * FROM fasta_nuc(?, bin_width := 50{options})", first_sequence)
        self.compare_over("read_hts_index sequences", fastas, "SELECT * FROM read_hts_index(?{options})")
        self.compare(
            "read_fastq paired",
            "SELECT * FROM read_fastq(?, mate_path := ?)",
            [str(DATA / "r1.fq"), str(DATA / "r2.fq")],
        )
        self.compare(
            "fasta_nuc bed",
            "SELECT * FROM fasta_nuc(?, bed_path := ?)",
            [str(DATA / "fixture_ref.fa"), str(DATA / "fixture_mixed_regions.bed")],
        )
        self.assert_readers_compared(
            {"read_fasta": 25, "read_fasta region": 20, "read_fasta count": 25, "read_fastq": 5,
             "detect_quality_encoding": 5, "fasta_nuc": 20, "fasta_nuc region": 20,
             "read_fastq paired": 1, "fasta_nuc bed": 1}
        )

    def test_genbank_bigwig_score_and_liftover(self):
        self.compare_over("read_genbank", fixtures("*.gb"), "SELECT * FROM read_genbank(?{options})")
        self.compare("read_bigwig", "SELECT * FROM read_bigwig(?)", [str(BIGWIG)])
        self.compare("read_bigwig region", "SELECT * FROM read_bigwig(?, region := '1')", [str(BIGWIG)])
        self.compare("read_bigwig", "SELECT * FROM read_bigwig(?)", [str(DATA / "range.bam")], "not a bigWig")
        score = "SELECT * FROM bcftools_score(?, ?, use := 'GT', columns := 'PLINK'{options})"
        inputs = [str(DATA / "score_input.vcf"), str(DATA / "score_summary.tsv")]
        self.compare("bcftools_score", score.format(options=""), inputs)
        self.compare(
            "bcftools_score",
            score.format(options=", regions_file := " + quoted(DATA / "score_regions.txt")),
            inputs,
            "regions_file",
        )
        self.compare(
            "bcftools_score",
            "SELECT * FROM bcftools_score(?, [?, ?], use := 'GT', columns := 'PLINK')",
            inputs + [str(DATA / "score_summary_na.tsv")],
            "two summaries",
        )
        self.compare(
            "bcftools_score",
            "SELECT * FROM bcftools_score(?, ?, use := 'GT')",
            [str(DATA / "score_input.vcf"), str(DATA / "score_gwas_summary.vcf")],
            "GWAS-VCF summary",
        )
        # include and exclude read a @file sample list through the file access.
        samples = Path(self.scratch.name) / "score_samples.txt"
        samples.write_text("S2\n")
        compressed = Path(self.scratch.name) / "score_samples.txt.gz"
        with gzip.open(compressed, "wt") as handle:
            handle.write("S1\n\nS2\n")
        for option, sample_list in (("include", samples), ("exclude", samples), ("include", compressed)):
            expression = "GT[@" + str(sample_list) + ']="alt"'
            self.compare(
                "bcftools_score",
                score.format(options=", " + option + " := " + quoted(expression)),
                inputs,
                option + " with a @file sample list " + sample_list.name,
            )
        self.compare(
            "bcftools_score",
            score.format(options=", include := " + quoted('GT[@' + str(samples) + '.missing]="alt"')),
            inputs,
            "a missing @file sample list",
        )
        for chain in fixtures("*.chain"):
            destination = chain.with_name(chain.stem + "_dst.fa")
            source = chain.with_name(chain.stem + "_src.fa")
            if not destination.exists():
                destination = source = chain.with_suffix(".fa")
            if not destination.exists() or not source.exists():
                continue
            self.compare(
                "bcftools_liftover",
                "SELECT bcftools_liftover(NAME, position, 'A', 'C', ?, ?, ?, 1, 250, false, NULL::BIGINT, false)"
                " FROM read_fasta(?), range(1, 9) positions(position)",
                [str(chain), str(destination), str(source), str(source)],
                chain.name,
            )
        scratch = Path(self.scratch.name)
        converted = {}
        for route, connection in (("htslib", self.oracle), ("duckdb", self.routed)):
            directory = scratch / ("fasta_" + route)
            directory.mkdir()
            rows = [
                connection.execute(
                    "SELECT success, records_written FROM genbank_to_fasta(?, output_path := ?)",
                    [str(record), str(directory / (record.name + ".fa"))],
                ).fetchall()
                for record in fixtures("phix174.gb", "genbank_circular_wrap.gb")
            ]
            converted[route] = (rows, file_digests(directory))
        self.assertEqual(len(converted["htslib"][1]), 2)
        self.assertEqual(converted["duckdb"], converted["htslib"])
        self.compared[("genbank_to_fasta", "rows")] += 2
        self.assert_readers_compared(
            {"read_genbank": 10, "read_bigwig": 1, "read_bigwig region": 1, "bcftools_score": 7, "bcftools_liftover": 5}
        )

    def test_a_missing_input_fails_through_both_file_systems_with_duckdbs_reason(self):
        missing = str(DATA / "no-such-input")
        present = {
            "bam": str(INDEXED),
            "fasta": str(DATA / "fixture_ref.fa"),
            "bed": str(DATA / "fixture_mixed_regions.bed"),
            "vcf": str(DATA / "score_input.vcf"),
            "chain": str(DATA / "liftover.chain"),
        }
        scratch = Path(self.scratch.name)
        # A real panel, so duckhts_somalier_bam_counts reaches the alignment open.
        panel = scratch / "missing_input_panel.parquet"
        self.oracle.execute(
            "COPY (SELECT * FROM (VALUES ('WBcel235', 0, 'CHROMOSOME_I', 1, 'A', 'G'))"
            " panel(assembly, site_index, region, position, allele_a, allele_b)) TO "
            + quoted(panel) + " (FORMAT parquet)"
        )
        present["panel"] = str(panel)
        # Each relation names the missing file once; the message must carry DuckDB's reason for it.
        with_reason = {
            "read_bam": ("SELECT * FROM read_bam(?)", [missing]),
            "bam_bin_counts": ("SELECT * FROM bam_bin_counts(?, 100)", [missing]),
            "read_pileup": ("SELECT * FROM read_pileup(?, region := 'chr1')", [missing]),
            "duckhts_bam_insert_size_counts": ("SELECT * FROM duckhts_bam_insert_size_counts(?)", [missing]),
            "duckhts_bam_mate_overlap_counts": ("SELECT * FROM duckhts_bam_mate_overlap_counts(?)", [missing]),
            "duckhts_bam_mismatch_counts": (
                "SELECT * FROM duckhts_bam_mismatch_counts(?, ?)", [missing, present["fasta"]]),
            "duckhts_bam_mismatch_counts mask": (
                "SELECT * FROM duckhts_bam_mismatch_counts(?, ?, mask := ?)",
                [str(DATA / "bam_mismatch.bam"), str(DATA / "bam_mismatch.fa"), missing]),
            "duckhts_bam_bed_coverage": ("SELECT * FROM duckhts_bam_bed_coverage(?, ?)", [missing, present["bed"]]),
            "duckhts_bam_bed_coverage bed": (
                "SELECT * FROM duckhts_bam_bed_coverage(?, ?)", [str(DATA / "fixture_mixed.bam"), missing]),
            "duckhts_samtools_idxstats": (
                "SELECT * FROM duckhts_samtools_idxstats(?, output := ?)", [missing, str(scratch / "none.idxstats")]),
            "duckhts_mosdepth": ("SELECT * FROM duckhts_mosdepth(?, ?)", [str(scratch / "none.mosdepth"), missing]),
            "read_hts_header": ("SELECT * FROM read_hts_header(?)", [missing]),
            "read_hts_index": ("SELECT * FROM read_hts_index(?)", [missing]),
            "read_hts_index_spans": ("SELECT * FROM read_hts_index_spans(?)", [missing]),
            "read_bcf": ("SELECT * FROM read_bcf(?)", [missing]),
            "read_geno": ("SELECT * FROM read_geno(?)", [missing]),
            "read_bcf_samples": ("SELECT * FROM read_bcf_samples(?)", [missing]),
            "read_tabix": ("SELECT * FROM read_tabix(?)", [missing]),
            "read_gff": ("SELECT * FROM read_gff(?)", [missing]),
            "read_gtf": ("SELECT * FROM read_gtf(?)", [missing]),
            "read_bed": ("SELECT * FROM read_bed(?)", [missing]),
            "read_fasta": ("SELECT * FROM read_fasta(?)", [missing]),
            "read_fastq": ("SELECT * FROM read_fastq(?)", [missing]),
            "read_fastq mate": ("SELECT * FROM read_fastq(?, mate_path := ?)", [str(DATA / "r1.fq"), missing]),
            "detect_quality_encoding": ("SELECT * FROM detect_quality_encoding(?)", [missing]),
            "fasta_nuc bed": ("SELECT * FROM fasta_nuc(?, bed_path := ?)", [present["fasta"], missing]),
            "fasta_nuc": ("SELECT * FROM fasta_nuc(?, bin_width := 50)", [missing]),
            "fasta_nuc bed missing fasta": (
                "SELECT * FROM fasta_nuc(?, bed_path := ?)", [missing, present["bed"]]),
            "read_genbank": ("SELECT * FROM read_genbank(?)", [missing]),
            "genbank_to_fasta": (
                "SELECT * FROM genbank_to_fasta(?, output_path := ?)", [missing, str(scratch / "none.fa")]),
            "read_bigwig": ("SELECT * FROM read_bigwig(?)", [missing]),
            "bcftools_score": (
                "SELECT * FROM bcftools_score(?, ?, use := 'GT', columns := 'PLINK')",
                [missing, str(DATA / "score_summary.tsv")]),
            "bcftools_score summary": (
                "SELECT * FROM bcftools_score(?, ?, use := 'GT', columns := 'PLINK')", [present["vcf"], missing]),
            "bcftools_liftover": (
                "SELECT bcftools_liftover('chrF', 2, 'C', 'T', ?, ?, ?, 1, 250, false, NULL::BIGINT, false)",
                [missing, str(DATA / "liftover_dst.fa"), str(DATA / "liftover_src.fa")]),
            "bcftools_norm_row": (NORM_ONE_ROW, [missing]),
            "bcftools_munge_row": (MUNGE_ONE_ROW, [missing]),
        }
        # These report the failure without a slot for the reason.
        class_only = {
            "duckhts_somalier_bam_counts": (
                "SELECT * FROM duckhts_somalier_bam_counts(?, NULL, 'sample-1', ?, panel_parquet := ?)",
                [missing, str(DATA / "ce.fa"), present["panel"]]),
        }
        for reader, (relation, parameters) in {**with_reason, **class_only}.items():
            with self.subTest(reader=reader):
                expected = relation_digest(self.oracle, relation, parameters)
                self.assertIsInstance(expected, str)
                with self.assertRaises(duckdb.Error) as failure:
                    self.routed.execute(relation, parameters).fetchall()
                message = str(failure.exception)
                if reader == "bcftools_score summary":
                    # This message quotes strerror(errno), and a failed open through
                    # DuckDB sets EIO whatever DuckDB's reason is.
                    expected = expected.replace("No such file or directory", "Input/output error")
                self.assertEqual(
                    type(failure.exception).__name__ + ": " + DUCKDB_REASON.sub("", message.splitlines()[0]), expected
                )
                self.assertNotIn("duckdbfs:", message)
                if reader in with_reason:
                    self.assertIn("No such file or directory", message)
                self.compared[(reader, "error")] += 1

    def test_duckdb_access_policy_applies_to_every_reader_family(self):
        relations = {
            "read_bcf": ("SELECT * FROM read_bcf(?)", [str(DATA / "score_input.vcf")]),
            "read_tabix": ("SELECT * FROM read_tabix(?)", [str(DATA / "fixture_mixed_regions.bed")]),
            "read_bed": ("SELECT * FROM read_bed(?)", [str(DATA / "fixture_mixed_regions.bed")]),
            "read_fasta": ("SELECT * FROM read_fasta(?)", [str(DATA / "fixture_ref.fa")]),
            "read_fasta region": ("SELECT * FROM read_fasta(?, region := '11')", [str(DATA / "fixture_ref.fa")]),
            "read_genbank": ("SELECT * FROM read_genbank(?)", [str(DATA / "phix174.gb")]),
            "read_bigwig": ("SELECT * FROM read_bigwig(?)", [str(BIGWIG)]),
            "read_hts_header": ("SELECT * FROM read_hts_header(?)", [str(INDEXED)]),
            "bam_bin_counts": ("SELECT * FROM bam_bin_counts(?, 100)", [str(INDEXED)]),
            "bcftools_liftover": (
                "SELECT bcftools_liftover('chrF', 2, 'C', 'T', ?, ?, ?, 1, 250, false, NULL::BIGINT, false)",
                [str(DATA / "liftover.chain"), str(DATA / "liftover_dst.fa"), str(DATA / "liftover_src.fa")]),
            "bcftools_norm_row": (NORM_ONE_ROW, [str(DATA / "liftover_src.fa")]),
            "bcftools_munge_row": (MUNGE_ONE_ROW, [str(DATA / "liftover_src.fa")]),
        }
        with connect("duckdb") as connection:
            for reader, (relation, parameters) in relations.items():
                with self.subTest(reader=reader, access="allowed"):
                    self.assertIsInstance(relation_digest(connection, relation, parameters), tuple)
            connection.execute("SET enable_external_access = false")
            for reader, (relation, parameters) in relations.items():
                with self.subTest(reader=reader, access="disabled"):
                    with self.assertRaises(duckdb.Error):
                        connection.execute(relation, parameters).fetchall()

    def test_duckdb_access_policy_applies_to_a_somalier_panel(self):
        """panel_parquet is read inside allowed_directories or allowed_paths only."""
        with tempfile.TemporaryDirectory() as root:
            inside, outside = Path(root, "in'side"), Path(root, "outside")
            inside.mkdir()
            outside.mkdir()
            for name in ("range.bam", "range.bam.bai", "ce.fa", "ce.fa.fai"):
                shutil.copy(DATA / name, inside / name)
            for directory in (inside, outside):
                self.oracle.execute(
                    "COPY (SELECT * FROM (VALUES ('WBcel235', 0, 'CHROMOSOME_I', 1, 'A', 'G'),"
                    " ('WBcel235', 1, 'CHROMOSOME_I', 914, 'A', 'C'))"
                    " panel(assembly, site_index, region, position, allele_a, allele_b)) TO "
                    + quoted(directory / "panel.parquet") + " (FORMAT parquet)"
                )
            relation = (
                "SELECT * FROM duckhts_somalier_bam_counts(?, NULL, 's', ?, panel_parquet := ?) ORDER BY ALL"
            )

            def counts(file_system, panel, restricted, allowed_paths=()):
                with connect(file_system) as connection:
                    if restricted:
                        connection.execute("SET allowed_directories = [" + quoted(str(inside) + "/") + "]")
                        if allowed_paths:
                            connection.execute(
                                "SET allowed_paths = [" + ", ".join(quoted(path) for path in allowed_paths) + "]")
                        connection.execute("SET enable_external_access = false")
                    parameters = [str(inside / "range.bam"), str(inside / "ce.fa"), str(panel)]
                    return connection.execute(relation, parameters).fetchall()

            expected = counts("htslib", inside / "panel.parquet", restricted=False)
            self.assertEqual(len(expected), 2)
            for panel in (inside / "panel.parquet", outside / "panel.parquet"):
                with self.subTest(route="htslib", panel=panel.parent.name):
                    self.assertEqual(counts("htslib", panel, restricted=True), expected)
            with self.subTest(route="duckdb", panel="inside"):
                self.assertEqual(counts("duckdb", inside / "panel.parquet", restricted=True), expected)
            with self.subTest(route="duckdb", panel="outside"):
                with self.assertRaisesRegex(
                        duckdb.Error, 'Permission Error: Cannot access file "' + re.escape(str(outside / "panel.parquet"))):
                    counts("duckdb", outside / "panel.parquet", restricted=True)
            with self.subTest(route="duckdb", panel="outside, in allowed_paths"):
                self.assertEqual(
                    counts("duckdb", outside / "panel.parquet", restricted=True,
                           allowed_paths=[outside / "panel.parquet"]),
                    expected)

    def test_duckdb_access_policy_applies_to_a_score_sample_list(self):
        """A @file sample list in include is read like the inputs: inside allowed_directories only."""
        with tempfile.TemporaryDirectory() as root:
            inside, outside = Path(root, "inside"), Path(root, "outside")
            inside.mkdir()
            outside.mkdir()
            for name in ("score_input.vcf", "score_summary.tsv"):
                shutil.copy(DATA / name, inside / name)
            for directory in (inside, outside):
                (directory / "samples.txt").write_text("S2\n")
            relation = (
                "SELECT * FROM bcftools_score(?, ?, use := 'GT', columns := 'PLINK', include := ?) ORDER BY ALL"
            )

            def score(file_system, sample_list, restricted):
                with connect(file_system) as connection:
                    if restricted:
                        connection.execute("SET allowed_directories = [" + quoted(str(inside) + "/") + "]")
                        connection.execute("SET enable_external_access = false")
                    expression = "GT[@" + str(sample_list) + ']="alt"'
                    parameters = [str(inside / "score_input.vcf"), str(inside / "score_summary.tsv"), expression]
                    return connection.execute(relation, parameters).fetchall()

            expected = score("htslib", inside / "samples.txt", restricted=False)
            self.assertEqual(len(expected), 2)
            # The sample list must change the result, or reading it would prove nothing.
            (inside / "everyone.txt").write_text("S1\nS2\n")
            self.assertNotEqual(score("htslib", inside / "everyone.txt", restricted=False), expected)
            for sample_list in (inside / "samples.txt", outside / "samples.txt"):
                with self.subTest(route="htslib", sample_list=sample_list.parent.name):
                    self.assertEqual(score("htslib", sample_list, restricted=True), expected)
            with self.subTest(route="duckdb", sample_list="inside"):
                self.assertEqual(score("duckdb", inside / "samples.txt", restricted=True), expected)
            with self.subTest(route="duckdb", sample_list="outside"):
                # The error names the list and carries DuckDB's reason for refusing it.
                refused = re.escape(str(outside / "samples.txt"))
                with self.assertRaisesRegex(duckdb.Error,
                                            "Could not read: " + refused + r' \(Cannot access file "' + refused):
                    score("duckdb", outside / "samples.txt", restricted=True)

    def test_variant_and_tabix_index_forms(self):
        """An index beside the file under each name htslib looks for, an explicit one, and none."""
        cases = (
            ("read_bcf", DATA / "formatcols.vcf.gz", "variants.vcf.gz",
             ("variants.vcf.gz.tbi", "variants.vcf.gz.csi", "variants.vcf.csi", "variants.vcf.tbi"), "CHROM"),
            ("read_bcf", DATA / "geno_sites.bcf", "variants.bcf", ("variants.bcf.csi", "variants.csi"), "CHROM"),
            ("read_gff", DATA / "gff_named_attributes.gff3.gz", "features.gff3.gz",
             ("features.gff3.gz.tbi", "features.gff3.tbi"), "seqname"),
        )
        for reader, source, name, index_names, region_column in cases:
            source_index = next(
                index for index in (Path(str(source) + ".tbi"), Path(str(source) + ".csi")) if index.exists()
            )
            relation = "SELECT * FROM " + reader + "(?{options})"
            region = first_value(self.oracle, "SELECT " + region_column + " FROM " + reader + "(?) LIMIT 1", source)
            self.assertIsNotNone(region)
            options = ", region := " + quoted(region)
            expected = relation_digest(self.oracle, relation.format(options=options), [str(source)])
            self.assertIsInstance(expected, tuple)
            with tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                copy = directory / name
                shutil.copy(source, copy)
                with self.subTest(reader=reader, file=name, index="none"):
                    without_index = relation_digest(self.oracle, relation.format(options=options), [str(copy)])
                    self.assertIsInstance(without_index, str)
                    self.assertEqual(
                        relation_digest(self.routed, relation.format(options=options), [str(copy)]), without_index
                    )
                with self.subTest(reader=reader, file=name, index="explicit index_path"):
                    explicit = options + ", index_path := " + quoted(source_index)
                    self.assertEqual(
                        relation_digest(self.routed, relation.format(options=explicit), [str(copy)]), expected
                    )
                for index_name in index_names:
                    if not index_name.endswith(source_index.suffix):
                        continue
                    with self.subTest(reader=reader, file=name, index=index_name):
                        shutil.copy(source_index, directory / index_name)
                        beside = relation_digest(self.oracle, relation.format(options=options), [str(copy)])
                        self.assertEqual(beside, expected)
                        self.assertEqual(
                            relation_digest(self.routed, relation.format(options=options), [str(copy)]), expected
                        )
                        (directory / index_name).unlink()

    def test_variant_index_lookup_ignores_the_working_directory(self):
        source = DATA / "formatcols.vcf.gz"
        relation = "SELECT * FROM read_bcf(?, region := '1')"
        expected = relation_digest(self.oracle, relation, [str(source)])
        self.assertIsInstance(expected, tuple)
        with tempfile.TemporaryDirectory() as directory, contextlib.chdir(directory):
            decoys = [source.name + suffix for suffix in (".tbi", ".csi")]
            for decoy in decoys:
                Path(decoy).write_bytes(b"not an index")
            with connect("duckdb") as routed:
                self.assertEqual(relation_digest(routed, relation, [str(source)]), expected)
            self.assertEqual(sorted(os.listdir(".")), sorted(decoys))

    def test_the_duckdb_route_builds_no_missing_fasta_index(self):
        """bcftools_liftover builds a missing .fai with htslib; the DuckDB route is read-only."""
        relation = "SELECT bcftools_liftover('chrF', 2, 'C', 'T', ?, ?, ?, 1, 250, false, NULL::BIGINT, false)"
        expected = relation_digest(
            self.oracle,
            relation,
            [str(DATA / "liftover.chain"), str(DATA / "liftover_dst.fa"), str(DATA / "liftover_src.fa")],
        )
        self.assertIsInstance(expected, tuple)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for name in ("liftover.chain", "liftover_dst.fa", "liftover_src.fa"):
                shutil.copy(DATA / name, directory / name)
            parameters = [str(directory / "liftover.chain"), str(directory / "liftover_dst.fa"),
                          str(directory / "liftover_src.fa")]
            with self.assertRaisesRegex(duckdb.Error, "failed to load destination FASTA index"):
                self.routed.execute(relation, parameters).fetchall()
            self.assertEqual(sorted(entry.name for entry in directory.iterdir()),
                             ["liftover.chain", "liftover_dst.fa", "liftover_src.fa"])
            self.assertEqual(relation_digest(self.oracle, relation, parameters), expected)
            self.assertIn("liftover_dst.fa.fai", [entry.name for entry in directory.iterdir()])
            self.assertEqual(relation_digest(self.routed, relation, parameters), expected)

    def test_http_readers_and_cram_reference_use_the_connection_secret(self):
        with connect("duckdb") as probe:
            try:
                probe.execute("LOAD httpfs")
            except duckdb.Error as failure:
                self.skipTest("httpfs is not installed for this DuckDB: " + str(failure).splitlines()[0])
        relations = {
            "read_bam reference": (
                "SELECT * FROM read_bam({0}range.cram', reference := {0}ce.fa')",
                ("/protected/range.cram", "/protected/ce.fa", "/protected/ce.fa.fai")),
            "read_bam reference region": (
                "SELECT * FROM read_bam({0}range.cram', reference := {0}ce.fa', region := 'CHROMOSOME_I')",
                ("/protected/range.cram.crai", "/protected/ce.fa")),
            "read_bcf region": (
                "SELECT * FROM read_bcf({0}formatcols.vcf.gz', region := '1')",
                ("/protected/formatcols.vcf.gz", "/protected/formatcols.vcf.gz.csi")),
            "read_gff region": (
                "SELECT * FROM read_gff({0}gff_named_attributes.gff3.gz', region := 'chr1')",
                ("/protected/gff_named_attributes.gff3.gz.tbi",)),
            "read_fasta region": (
                "SELECT * FROM read_fasta({0}ce.fa', region := 'CHROMOSOME_I:1-200')", ("/protected/ce.fa.fai",)),
            "read_bed": ("SELECT * FROM read_bed({0}fixture_mixed_regions.bed')", ("/protected/fixture_mixed_regions.bed",)),
            "read_genbank": ("SELECT * FROM read_genbank({0}phix174.gb')", ("/protected/phix174.gb",)),
            "bcftools_liftover": (
                "SELECT bcftools_liftover('chrF', 2, 'C', 'T', {0}liftover.chain', {0}liftover_dst.fa',"
                " {0}liftover_src.fa', 1, 250, false, NULL::BIGINT, false)",
                ("/protected/liftover.chain", "/protected/liftover_dst.fa.fai")),
        }
        server = ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            base = "http://127.0.0.1:" + str(server.server_address[1])
            with connect("duckdb") as routed:
                routed.execute("LOAD httpfs")
                for reader, (relation, _) in relations.items():
                    with self.subTest(reader=reader, secret="absent"):
                        with self.assertRaises(duckdb.Error):
                            routed.execute(relation.format("'" + base + "/protected/")).fetchall()
                routed.execute("CREATE SECRET (TYPE http, BEARER_TOKEN 'sesame', SCOPE '" + base + "/protected')")
                for reader, (relation, requested) in relations.items():
                    with self.subTest(reader=reader, secret="present"):
                        expected = relation_digest(self.oracle, relation.format("'" + str(DATA) + "/"))
                        self.assertIsInstance(expected, tuple)
                        RangeHandler.requested_paths.clear()
                        self.assertEqual(
                            relation_digest(routed, relation.format("'" + base + "/protected/")), expected
                        )
                        for path in requested:
                            self.assertGreater(RangeHandler.requested_paths[path], 0, path)
        finally:
            server.shutdown()
            server.server_close()
            thread.join()

    def test_invalid_settings_are_errors_in_every_reader_family(self):
        relations = (
            ("read_bcf", "SELECT * FROM read_bcf(?)", [str(DATA / "score_input.vcf")]),
            ("read_tabix", "SELECT * FROM read_tabix(?)", [str(DATA / "fixture_mixed_regions.bed")]),
            ("read_fasta", "SELECT * FROM read_fasta(?)", [str(DATA / "fixture_ref.fa")]),
            ("read_bigwig", "SELECT * FROM read_bigwig(?)", [str(BIGWIG)]),
            ("bcftools_liftover",
             "SELECT bcftools_liftover('chrF', 2, 'C', 'T', ?, ?, ?, 1, 250, false, NULL::BIGINT, false)",
             [str(DATA / "liftover.chain"), str(DATA / "liftover_dst.fa"), str(DATA / "liftover_src.fa")]),
            ("bcftools_norm_row", NORM_ONE_ROW, [str(DATA / "liftover_src.fa")]),
            ("bcftools_munge_row", MUNGE_ONE_ROW, [str(DATA / "liftover_src.fa")]),
        )
        with connect("tabix") as connection:
            for reader, relation, parameters in relations:
                with self.subTest(reader=reader, setting="duckhts_file_system"):
                    with self.assertRaisesRegex(
                        duckdb.Error, reader + ": duckhts_file_system must be 'htslib' or 'duckdb'"
                    ):
                        connection.execute(relation, parameters).fetchall()
        with connect() as connection:
            connection.execute("SET duckhts_remote_block_size = 1")
            for reader, relation, parameters in relations:
                with self.subTest(reader=reader, setting="duckhts_remote_block_size"):
                    with self.assertRaisesRegex(
                        duckdb.Error,
                        reader + ": duckhts_remote_block_size must be 0 or between 131072 and 268435456 bytes",
                    ):
                        connection.execute(relation, parameters).fetchall()


class ReferenceScalarsThroughBothFileSystems(unittest.TestCase):
    """bcftools_norm_row and bcftools_munge_row give htslib's rows through the DuckDB file system."""

    @classmethod
    def setUpClass(cls):
        cls.oracle = connect("htslib")
        cls.routed = connect("duckdb")
        cls.scratch = tempfile.TemporaryDirectory()
        # Only FASTAs whose indexes are committed, so that htslib builds nothing in test/data.
        cls.fastas = [
            fasta
            for fasta in fixtures("*.fa", "*.fa.gz")
            if Path(str(fasta) + ".fai").is_file()
            and (fasta.suffix != ".gz" or Path(str(fasta) + ".gzi").is_file())
        ]
        cls.cases = {fasta: reference_cases(fasta) for fasta in cls.fastas}
        cls.norm_count, cls.munge_count = load_reference_cases(
            (cls.oracle, cls.routed), cls.scratch.name, cls.cases
        )

    @classmethod
    def tearDownClass(cls):
        cls.oracle.close()
        cls.routed.close()
        cls.scratch.cleanup()

    def rows_through_both(self, relation, parameters=()):
        """Asserts the rows of relation are htslib's through the DuckDB route; returns them."""
        expected = self.oracle.execute(relation, list(parameters)).fetchall()
        observed = self.routed.execute(relation, list(parameters)).fetchall()
        self.assertEqual(len(observed), len(expected))
        for observed_row, expected_row in zip(observed, expected):
            self.assertEqual(observed_row, expected_row)
        return expected

    def test_fixture_set(self):
        compressed = [fasta for fasta in self.fastas if fasta.suffix == ".gz"]
        self.assertGreaterEqual(len(self.fastas), 25)
        self.assertGreaterEqual(len(compressed), 2)
        self.assertGreater(self.norm_count, 20000)
        self.assertGreater(self.munge_count, 10000)

    def test_bcftools_norm_row_rows_are_the_same(self):
        rows = self.rows_through_both(norm_cases())
        self.assertEqual(len(rows), self.norm_count)
        status = collections.Counter((kind, normed["norm_status"]) for _, kind, normed in rows)
        left_aligned = collections.Counter()
        positions = dict(self.oracle.execute("SELECT id, pos FROM norm_cases").fetchall())
        for identifier, kind, normed in rows:
            if normed["norm_status"] == "Normalized" and normed["pos_normed"] < positions[identifier]:
                left_aligned[kind] += 1
        print("\nbcftools_norm_row: %d rows equal over %d FASTAs; statuses %s; left-aligned %s" % (
            len(rows), len(self.fastas), dict(collections.Counter(s for _, s in status.elements())),
            dict(left_aligned)))
        # The oracle's rows must exercise normalisation, not only pass-through.
        for kind in ("deletion", "deletion_shared_tail", "insertion", "multi_allelic", "symbolic_del"):
            with self.subTest(kind=kind):
                self.assertGreaterEqual(left_aligned[kind], 100)
        self.assertGreaterEqual(status[("snp", "Unchanged")], 1000)
        self.assertGreaterEqual(status[("ref_mismatch", "RefMismatch")], 1000)
        self.assertGreaterEqual(status[("contig_alias", "Unchanged")], 10)
        self.assertEqual(status[("missing_contig", "MissingContig")], len(self.fastas))
        as_list = self.rows_through_both(norm_cases(alt="string_split(alt, ',')"))
        self.assertEqual(as_list, rows)

    def test_bcftools_munge_row_rows_are_the_same(self):
        rows = self.rows_through_both(munge_cases())
        self.assertEqual(len(rows), self.munge_count)
        outcome_of = collections.Counter(
            (kind, munged["alleles_swapped"], munged["filter"]) for _, kind, munged in rows
        )
        print("\nbcftools_munge_row: %d rows equal over %d FASTAs; (swapped, filter) %s" % (
            len(rows), len(self.fastas),
            dict(collections.Counter((swapped, tag) for _, swapped, tag in outcome_of.elements()))))
        self.assertGreaterEqual(outcome_of[("a2_is_reference", False, None)], 1000)
        self.assertGreaterEqual(outcome_of[("a1_is_reference", True, None)], 1000)
        self.assertGreaterEqual(outcome_of[("neither_is_reference", False, "REF_MISMATCH")], 1000)
        self.assertGreaterEqual(outcome_of[("both_match", False, "IFFY")], 1000)
        self.assertGreaterEqual(outcome_of[("contig_alias", False, None)], 10)
        self.assertEqual(outcome_of[("missing_contig", False, "MissingContig")], len(self.fastas))

    def test_more_references_than_one_worker_keeps_open(self):
        """Rows that name another FASTA each time, over more files than a worker holds handles for."""
        plain = [fasta for fasta in self.fastas if fasta.suffix != ".gz"][:12]
        self.assertGreater(len(plain), 8)  # DUCKHTS_REFERENCE_CACHE_ENTRIES
        chosen = ", ".join(quoted(fasta) for fasta in plain)
        # Three rows of each file; one thread evaluates them file after file, three times round.
        round_robin = (
            "(SELECT * FROM (SELECT *, row_number() OVER (PARTITION BY fasta ORDER BY id) AS turn FROM {table}"
            " WHERE fasta IN (" + chosen + ")) WHERE turn <= 3 ORDER BY turn, fasta) cases"
        )
        for connection in (self.oracle, self.routed):
            connection.execute("SET threads = 1")
        try:
            for relation in (norm_cases(cases=round_robin.format(table="norm_cases")),
                             munge_cases(cases=round_robin.format(table="munge_cases"))):
                self.assertEqual(len(self.rows_through_both(relation)), 3 * len(plain))
        finally:
            for connection in (self.oracle, self.routed):
                connection.execute("RESET threads")

    def test_explicit_index_paths(self):
        """An index named by the fai and gzi arguments, kept apart from the FASTA."""
        for source in (DATA / "region_names.fa.gz", DATA / "bam_mismatch.long.fa.gz", DATA / "liftover_src.fa"):
            with self.subTest(fasta=source.name), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                (directory / "index").mkdir()
                copy = directory / ("reference" + "".join(source.suffixes[-2:] if source.suffix == ".gz"
                                                         else source.suffixes[-1:]))
                shutil.copy(source, copy)
                fai = directory / "index" / "names.fai"
                shutil.copy(str(source) + ".fai", fai)
                gzi = "NULL"
                if source.suffix == ".gz":
                    shutil.copy(str(source) + ".gzi", directory / "index" / "blocks.gzi")
                    gzi = quoted(directory / "index" / "blocks.gzi")
                listing = sorted(str(entry.relative_to(directory)) for entry in directory.rglob("*"))
                beside = " WHERE fasta = " + quoted(source)
                expected = self.oracle.execute(norm_cases(beside)).fetchall()
                self.assertGreater(len(expected), 30)
                self.assertTrue(any(normed["norm_status"] in ("Normalized", "Unchanged") for _, _, normed in expected))
                # The same rows, each naming the copy instead of the fixture.
                apart = norm_cases(beside, fasta=quoted(copy), fai=quoted(fai), gzi=gzi)
                self.assertEqual(self.rows_through_both(apart), expected)
                with self.assertRaises(duckdb.Error) as failure:
                    self.routed.execute(norm_cases(beside, fasta=quoted(copy))).fetchall()
                self.assertIn("bcftools_norm_row: failed to load FASTA index: FASTA '" + str(copy) + "'",
                              str(failure.exception))
                self.assertIn("No such file or directory", str(failure.exception))
                self.assertEqual(
                    sorted(str(entry.relative_to(directory)) for entry in directory.rglob("*")), listing)

    def test_the_duckdb_route_builds_no_missing_fasta_index(self):
        """Both functions build a missing .fai with htslib; the DuckDB route is read-only."""
        source = DATA / "liftover_src.fa"
        where = " WHERE fasta = " + quoted(source)
        for function, cases in (("bcftools_norm_row", norm_cases), ("bcftools_munge_row", munge_cases)):
            with self.subTest(function=function), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                copy = directory / "reference.fa"
                absent = directory / "absent.fa"
                shutil.copy(source, copy)
                expected = self.oracle.execute(cases(where)).fetchall()
                self.assertGreater(len(expected), 30)
                for path, missing in ((copy, str(copy) + ".fai"), (absent, str(absent) + ".fai")):
                    with self.subTest(fasta=path.name):
                        with self.assertRaises(duckdb.Error) as failure:
                            self.routed.execute(cases(where, fasta=quoted(path))).fetchall()
                        message = str(failure.exception)
                        self.assertIn(function + ": failed to load FASTA index: FASTA '" + str(path) + "'", message)
                        self.assertIn('(Cannot open file "' + missing + '": No such file or directory)', message)
                        self.assertNotIn("duckdbfs:", message)
                        self.assertEqual([entry.name for entry in directory.iterdir()], ["reference.fa"])
                relation = cases(where, fasta=quoted(copy))
                self.assertEqual(self.oracle.execute(relation).fetchall(), expected)
                self.assertEqual(sorted(entry.name for entry in directory.iterdir()),
                                 ["reference.fa", "reference.fa.fai"])
                self.assertEqual(self.routed.execute(relation).fetchall(), expected)

    def test_routed_references_stay_with_their_query(self):
        """A reference opened through one connection serves no later query, of either route."""
        fasta = DATA / "ce.fa"
        where = " WHERE fasta = " + quoted(fasta) + " AND chrom = 'CHROMOSOME_I' AND pos {range}"
        near = where.format(range="<= 300")
        far = where.format(range="> 4 * %d" % REFERENCE_WINDOW_BASES)
        for function, cases in (("bcftools_norm_row", norm_cases), ("bcftools_munge_row", munge_cases)):
            with self.subTest(function=function):
                expected_near = self.oracle.execute(cases(near)).fetchall()
                expected_far = self.oracle.execute(cases(far)).fetchall()
                self.assertGreater(len(expected_near), 500)
                self.assertGreater(len(expected_far), 500)
                first = connect("duckdb")
                second = connect("htslib")
                try:
                    # One thread runs every query below, so they share one htslib thread cache.
                    for connection in (first, second):
                        connection.execute("SET threads = 1")
                        load_reference_cases((connection,), self.scratch.name, self.cases)
                    self.assertEqual(first.execute(cases(near)).fetchall(), expected_near)
                    first.execute("SET duckhts_file_system = 'htslib'")
                    self.assertEqual(first.execute(cases(near)).fetchall(), expected_near)
                    self.assertEqual(second.execute(cases(near)).fetchall(), expected_near)
                    first.execute("SET duckhts_file_system = 'duckdb'")
                    self.assertEqual(first.execute(cases(far)).fetchall(), expected_far)
                    first.close()
                    # Bases outside every window read so far need a live handle.
                    self.assertEqual(second.execute(cases(far)).fetchall(), expected_far)
                    self.assertEqual(second.execute(cases(near)).fetchall(), expected_near)
                finally:
                    first.close()
                    second.close()

    def test_routed_references_close_with_their_query(self):
        """After a query through DuckDB this process holds the FASTA open no longer; htslib's cache does."""
        descriptors = Path("/proc/self/fd")
        if not descriptors.is_dir():
            self.skipTest("no /proc/self/fd to list open files")
        source = DATA / "liftover_src.fa"
        where = " WHERE fasta = " + quoted(source)

        def open_copies(directory):
            held = []
            for descriptor in descriptors.iterdir():
                try:
                    target = os.readlink(descriptor)
                except OSError:
                    continue
                if target.startswith(str(directory)):
                    held.append(Path(target).name)
            return sorted(held)

        for function, cases in (("bcftools_norm_row", norm_cases), ("bcftools_munge_row", munge_cases)):
            with self.subTest(function=function), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory).resolve()
                copy = directory / "reference.fa"
                shutil.copy(source, copy)
                shutil.copy(str(source) + ".fai", str(copy) + ".fai")
                expected = self.oracle.execute(cases(where)).fetchall()
                relation = cases(where, fasta=quoted(copy))
                with connect("duckdb") as connection:
                    connection.execute("SET threads = 1")
                    load_reference_cases((connection,), self.scratch.name, self.cases)
                    self.assertEqual(connection.execute(relation).fetchall(), expected)
                    connection.execute("SELECT 1").fetchall()
                    self.assertEqual(open_copies(directory), [])
                    # The same listing sees the handle htslib's thread cache keeps.
                    connection.execute("SET duckhts_file_system = 'htslib'")
                    self.assertEqual(connection.execute(relation).fetchall(), expected)
                    connection.execute("SELECT 1").fetchall()
                    self.assertEqual(open_copies(directory), ["reference.fa"])

    def test_a_reference_htslib_keeps_open_does_not_serve_the_duckdb_route(self):
        """htslib's cached handle must not answer a query that DuckDB's access policy forbids."""
        fasta = DATA / "liftover_src.fa"
        where = " WHERE fasta = " + quoted(fasta)
        for function, cases in (("bcftools_norm_row", norm_cases(where)), ("bcftools_munge_row", munge_cases(where))):
            with self.subTest(function=function), connect("htslib") as connection:
                connection.execute("SET threads = 1")
                load_reference_cases((connection,), self.scratch.name, self.cases)
                expected = self.oracle.execute(cases).fetchall()
                self.assertEqual(connection.execute(cases).fetchall(), expected)
                connection.execute("SET duckhts_file_system = 'duckdb'")
                self.assertEqual(connection.execute(cases).fetchall(), expected)
                connection.execute("SET enable_external_access = false")
                with self.assertRaisesRegex(duckdb.Error, "(?i)disabled|permission"):
                    connection.execute(cases).fetchall()
                connection.execute("SET duckhts_file_system = 'htslib'")
                self.assertEqual(connection.execute(cases).fetchall(), expected)

    def test_http_reference_uses_the_connection_secret(self):
        with connect("duckdb") as probe:
            try:
                probe.execute("LOAD httpfs")
            except duckdb.Error as failure:
                self.skipTest("httpfs is not installed for this DuckDB: " + str(failure).splitlines()[0])
        fasta = DATA / "ce.fa"
        where = " WHERE fasta = " + quoted(fasta) + " AND chrom = 'CHROMOSOME_II' AND pos <= 200"
        server = ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            url = "http://127.0.0.1:" + str(server.server_address[1]) + "/protected/ce.fa"
            with connect("duckdb") as routed:
                routed.execute("LOAD httpfs")
                load_reference_cases((routed,), self.scratch.name, self.cases)
                for function, cases in (("bcftools_norm_row", norm_cases), ("bcftools_munge_row", munge_cases)):
                    relation = cases(where, fasta=quoted(url))
                    with self.subTest(function=function, secret="absent"):
                        with self.assertRaisesRegex(duckdb.Error, function + ": failed to load FASTA index.*403"):
                            routed.execute(relation).fetchall()
                routed.execute(
                    "CREATE SECRET (TYPE http, BEARER_TOKEN 'sesame', SCOPE '" + url.rsplit("/", 1)[0] + "')")
                for function, cases in (("bcftools_norm_row", norm_cases), ("bcftools_munge_row", munge_cases)):
                    relation = cases(where, fasta=quoted(url))
                    with self.subTest(function=function, secret="present"):
                        expected = self.oracle.execute(cases(where)).fetchall()
                        self.assertGreater(len(expected), 500)
                        RangeHandler.requested_paths.clear()
                        self.assertEqual(routed.execute(relation).fetchall(), expected)
                        for path in ("/protected/ce.fa", "/protected/ce.fa.fai"):
                            self.assertGreater(RangeHandler.requested_paths[path], 0, path)
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


class VolatileScalarsDoNoWorkAtPlanning(unittest.TestCase):
    """bcftools_munge_row and bcftools_liftover with constant arguments wait for the query to run."""

    @staticmethod
    def calls(fasta):
        fasta = quoted(fasta)
        return {
            "bcftools_munge_row": (
                "SELECT bcftools_munge_row('chrF', 2, 'A', 'C', 'rs1', 0.01, NULL, NULL, 0.12, 1000, NULL, NULL, NULL,"
                " 0.3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, '+-', " + fasta + ", NULL, NULL, NULL, NULL, NULL)"
            ),
            "bcftools_liftover": (
                "SELECT bcftools_liftover('chrF', 2, 'C', 'T', " + quoted(DATA / "liftover.chain") + ", " + fasta + ", "
                + fasta + ", 1, 250, false, NULL::BIGINT, false)"
            ),
        }

    def test_prepare_and_explain_open_no_reference(self):
        for file_system in ("htslib", "duckdb"):
            with tempfile.TemporaryDirectory() as directory, connect(file_system) as connection:
                absent = Path(directory) / "absent.fa"
                for function, call in self.calls(absent).items():
                    with self.subTest(file_system=file_system, function=function):
                        connection.execute("PREPARE planned AS " + call)
                        connection.execute("EXPLAIN " + call).fetchall()
                        with self.assertRaisesRegex(duckdb.Error, function):
                            connection.execute(call).fetchall()
                        with self.assertRaisesRegex(duckdb.Error, function):
                            connection.execute("EXECUTE planned").fetchall()
                        self.assertEqual(os.listdir(directory), [])

    def test_prepare_and_explain_build_no_index(self):
        for function in ("bcftools_munge_row", "bcftools_liftover"):
            with self.subTest(function=function), tempfile.TemporaryDirectory() as directory:
                copy = Path(directory) / "reference.fa"
                shutil.copy(DATA / "liftover_src.fa", copy)
                call = self.calls(copy)[function]
                with connect("htslib") as connection:
                    connection.execute("PREPARE planned AS " + call)
                    connection.execute("EXPLAIN " + call).fetchall()
                    self.assertEqual(os.listdir(directory), ["reference.fa"])
                    self.assertEqual(len(connection.execute(call).fetchall()), 1)
                    self.assertEqual(sorted(os.listdir(directory)), ["reference.fa", "reference.fa.fai"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    args, remaining = parser.parse_known_args()
    args.extension = args.extension.resolve()
    unittest.main(argv=[__file__] + remaining)
