#!/usr/bin/env python3
"""Maintainer-only, network-free FILE_OFFSET fixtures; run from the repo root.

Inputs are the committed bam_scan_mixed.bam (synthetic authority in
prepare_bam_scan_fixtures.R) and nanopore.bam. samtools 1.23 / HTSlib 1.23
supplies SAM decoding, bgzip supplies BGZF SAM, and Python gzip supplies
ordinary gzip. Every derived alignment file must preserve decoded SAM rows.

The offset oracle walks physical BGZF blocks and BAM length prefixes without
DuckHTS or an HTSlib tell call. It records the position *after* each record,
normalizing an exhausted block to the next compressed block's start.
"""

import csv
import gzip
import shutil
import struct
import subprocess
from pathlib import Path


def record_ends(path):
    compressed = path.read_bytes()
    blocks = []
    decoded = bytearray()
    position = 0
    while position < len(compressed):
        # These committed BGZF inputs use the canonical six-byte BC extra field.
        assert compressed[position:position + 4] == b"\x1f\x8b\x08\x04"
        assert compressed[position + 10:position + 16] == b"\x06\x00BC\x02\x00"
        size = struct.unpack_from("<H", compressed, position + 16)[0] + 1
        block = gzip.decompress(compressed[position:position + size])
        blocks.append((len(decoded), len(decoded) + len(block), position, position + size))
        decoded.extend(block)
        position += size
    assert position == len(compressed) and decoded[:4] == b"BAM\1"
    offset = 8 + struct.unpack_from("<i", decoded, 4)[0]
    references = struct.unpack_from("<i", decoded, offset)[0]
    offset += 4
    for _ in range(references):
        offset += 8 + struct.unpack_from("<i", decoded, offset)[0]
    while offset < len(decoded):
        size = struct.unpack_from("<i", decoded, offset)[0]
        assert size >= 32 and offset + 4 + size <= len(decoded)
        name_len = decoded[offset + 12]
        qname = decoded[offset + 36:offset + 36 + name_len - 1].decode("ascii")
        offset += 4 + size
        start, end, address, next_address = next(b for b in blocks if b[0] < offset <= b[1])
        virtual = next_address << 16 if offset == end else (address << 16) | (offset - start)
        yield path.name, qname, virtual
    assert offset == len(decoded)


def main():
    data = Path("test/data")
    package = Path("r/Rduckhts/inst/extdata")
    source = data / "bam_scan_mixed.bam"
    sam = subprocess.check_output(["samtools", "view", "--no-PG", "-h", str(source)])
    plain_bam = gzip.decompress(source.read_bytes())
    files = {
        "bam_offset.sam": sam,
        "bam_offset.sam.gz": gzip.compress(sam, mtime=0),
        "bam_offset.sam.bgz": subprocess.check_output(["bgzip", "-c"], input=sam),
        "bam_offset_uncompressed.bam": plain_bam,
        "bam_offset_gzip.bam": gzip.compress(plain_bam, mtime=0),
    }
    expected = subprocess.check_output(["samtools", "view", str(source)])
    for name, content in files.items():
        path = data / name
        path.write_bytes(content)
        assert subprocess.check_output(["samtools", "view", str(path)]) == expected
        shutil.copyfile(path, package / name)
    oracle = data / "bam_record_ends.tsv"
    with oracle.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "QNAME", "FILE_OFFSET"))
        for name, count in (("bam_scan_mixed.bam", 5), ("nanopore.bam", 186)):
            rows = list(record_ends(data / name))
            assert len(rows) == count
            writer.writerows(rows)
    shutil.copyfile(oracle, package / oracle.name)


if __name__ == "__main__":
    main()
