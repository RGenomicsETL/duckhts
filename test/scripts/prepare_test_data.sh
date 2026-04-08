#!/usr/bin/env bash
# Prepare indexed test data for duckhts SQL tests.
#
# Copies files from vendored htslib test suite into test/data/ and
# builds the required indexes (BAI, CSI, TBI, FAI) so that region
# queries work without stderr noise.
#
# Prerequisites: samtools, bcftools, bgzip, tabix (all from htslib/samtools).
#
# Usage:  ./test/scripts/prepare_test_data.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SRC="$REPO_ROOT/third_party/htslib/test"
DST="$REPO_ROOT/test/data"
PKG_DST="$REPO_ROOT/r/Rduckhts/inst/extdata"

mkdir -p "$DST"
mkdir -p "$PKG_DST"

echo "==> Preparing test data in $DST"

# ---- BAM (copy + index) ----
cp "$SRC/range.bam" "$DST/range.bam"
samtools index "$DST/range.bam"
echo "  range.bam + .bai"

# ---- Parallel empty-contig BAM regression fixture ----
for out_dir in "$DST" "$PKG_DST"; do
  samtools view -b -o "$out_dir/parallel_empty_contigs.bam" \
    "$SCRIPT_DIR/parallel_empty_contigs.sam"
  samtools index "$out_dir/parallel_empty_contigs.bam"
done
echo "  parallel_empty_contigs.bam + .bai"

# ---- VCF → bgzipped VCF + index ----
bcftools view "$SRC/formatcols.vcf" -Oz -o "$DST/formatcols.vcf.gz"
bcftools index "$DST/formatcols.vcf.gz"
echo "  formatcols.vcf.gz + .csi"

# ---- BCF (copy + index) ----
cp "$SRC/tabix/vcf_file.bcf" "$DST/vcf_file.bcf"
bcftools index "$DST/vcf_file.bcf"
echo "  vcf_file.bcf + .csi"

# ---- FASTA (copy + index) ----
cp "$SRC/ce.fa" "$DST/ce.fa"
samtools faidx "$DST/ce.fa"
echo "  ce.fa + .fai"

# ---- FASTQ (copy as-is; no index needed) ----
cp "$SRC/fastq/r1.fq" "$DST/r1.fq"
echo "  r1.fq"

# ---- Legacy FASTQ quality fixture (Illumina Phred+64 text) ----
printf '@legacy1\nACGT\n+\nhhhh\n' > "$DST/legacy_phred64.fq"
echo "  legacy_phred64.fq"

# ---- GFF → bgzipped + tabix ----
bgzip -c "$SRC/tabix/gff_file.gff" > "$DST/gff_file.gff.gz"
tabix -p gff "$DST/gff_file.gff.gz"
echo "  gff_file.gff.gz + .tbi"

# ---- vcfppR-generated VCF compliance fixtures ----
Rscript "$SCRIPT_DIR/vcfpp.R"
echo "  vcfppR-generated VCF spec/mapping fixtures + manifest"

# ---- Parallel empty-contig VCF regression fixture (bgzip + tabix) ----
for out_dir in "$DST" "$PKG_DST"; do
  bgzip -c "$out_dir/parallel_empty_contigs.vcf" > "$out_dir/parallel_empty_contigs.vcf.gz"
  tabix -f -p vcf "$out_dir/parallel_empty_contigs.vcf.gz"
done
echo "  parallel_empty_contigs.vcf.gz + .tbi"

echo "==> Done. $(ls "$DST" | wc -l) files in test/data/"
