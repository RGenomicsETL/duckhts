#!/usr/bin/env python3
"""Create an indexed BAM with a valid first contig and an invalid second record."""

import pathlib
import struct
import tempfile
import zlib

import pysam

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE = ROOT / "test/data/coverage_decode_error.sam"
TARGET = SOURCE.with_suffix(".bam")

with pysam.AlignmentFile(str(SOURCE), "r") as source:
    with pysam.AlignmentFile(str(TARGET), "wb", header=source.header) as target:
        for record in source:
            target.write(record)
pysam.index(str(TARGET))

# Preserve BGZF block starts and in-block record offsets used by the good index.
# The malformed record is in the final data block, so changing its compressed
# size cannot move a subsequent record. Its CIGAR consumes one base, not two.
encoded = TARGET.read_bytes()
blocks = []
offset = 0
changed = False
while offset < len(encoded):
    size = struct.unpack_from("<H", encoded, offset + 16)[0] + 1
    block = encoded[offset:offset + size]
    data = bytearray(zlib.decompress(block, wbits=31))
    name_offset = data.find(b"bad_read\0")
    if name_offset >= 0:
        record_offset = name_offset - 36
        assert struct.unpack_from("<i", data, record_offset + 20)[0] == 1
        struct.pack_into("<i", data, record_offset + 20, 2)
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            compressed = pathlib.Path(temporary) / "block.bgz"
            with pysam.BGZFile(str(compressed), "wb") as out:
                out.write(bytes(data))
            # BGZFile.close appends the standard 28-byte EOF block.
            block = compressed.read_bytes()[:-28]
        changed = True
    elif changed:
        assert not data, "The corrupt record must be in the final data block"
    blocks.append(block)
    offset += size
assert changed
TARGET.write_bytes(b"".join(blocks))
