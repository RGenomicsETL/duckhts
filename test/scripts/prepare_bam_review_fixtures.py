"""Rebuild tiny BAM review fixtures with pysam (run from the repository root)."""
from pathlib import Path
import shutil
import pysam

DATA = Path("test/data")


def bam_from_sam(name, index=True):
    with pysam.AlignmentFile(str(DATA / f"{name}.sam"), "r") as source:
        with pysam.AlignmentFile(str(DATA / f"{name}.bam"), "wb", header=source.header) as output:
            for record in source:
                # SAM parsing drops RNAME for POS 0; set the BAM position directly.
                if name == "bam_pos_zero" and record.query_name.startswith("u"):
                    record.reference_start = -1
                output.write(record)
    if index:
        pysam.index(str(DATA / f"{name}.bam"))


for name in ("bam_pos_zero", "bam_wrong_tags", "bam_stream_contigs"):
    bam_from_sam(name)
shutil.copyfile(DATA / "bam_stream_contigs.bam", DATA / "bam_stream_contigs_noindex.bam")
bam_from_sam("bam_short_index")
with (DATA / "bam_stale_header.sam").open("w") as source:
    source.write("@HD\tVN:1.6\tSO:coordinate\n")
    for tid in range(1000):
        source.write(f"@SQ\tSN:chr{tid}\tLN:100\n")
bam_from_sam("bam_stale_header", index=False)
shutil.copyfile(DATA / "bam_short_index.bam.bai", DATA / "bam_stale_header.bam.bai")
