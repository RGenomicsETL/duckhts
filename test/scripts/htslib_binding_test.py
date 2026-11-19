"""Check that the host process cannot replace the htslib that DuckHTS bundles.

Run with test/scripts/htslib_binding_host.c preloaded (LD_PRELOAD). The host
library defines hopen() and seq_nt16_str ahead of the extension in the dynamic
loader's search order; read_fastq must still open the file through the bundled
hopen() and decode bases through the bundled table.
"""

import argparse
import ctypes
from pathlib import Path
import unittest

import duckdb

HOST_TABLE = b"X" * 16


def fastq_sequences(path):
    """Sequence lines of a FASTQ file whose records are four lines each."""
    lines = path.read_text().splitlines()
    return [line.upper() for line in lines[1::4]]


class HtslibBinding(unittest.TestCase):
    def test_host_symbols_do_not_replace_bundled_htslib(self):
        process = ctypes.CDLL(None)
        host = ctypes.CDLL(str(args.host.resolve()))

        # The preload is in effect: the process-wide symbols are the host's.
        self.assertEqual(
            ctypes.cast(process.hopen, ctypes.c_void_p).value,
            ctypes.cast(host.hopen, ctypes.c_void_p).value,
        )
        self.assertEqual((ctypes.c_char * 16).in_dll(process, "seq_nt16_str").raw, HOST_TABLE)

        fixture = Path(__file__).resolve().parents[1] / "data" / "interleaved.fq"
        expected = fastq_sequences(fixture)
        self.assertGreater(len(expected), 0)
        self.assertTrue(all(set(sequence) <= set("ACGTN") for sequence in expected))

        with duckdb.connect(config={"allow_unsigned_extensions": "true"}) as connection:
            connection.execute("LOAD '" + str(args.extension).replace("'", "''") + "'")
            rows = connection.execute(
                "SELECT sequence FROM read_fastq(?)", [str(fixture)]
            ).fetchall()

        self.assertEqual([row[0] for row in rows], expected)
        self.assertEqual(host.htslib_binding_host_hopen_calls(), 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--host", required=True, type=Path)
    args, remaining = parser.parse_known_args()
    unittest.main(argv=[__file__] + remaining)
