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

# ---- Upstream mosdepth conformance fixtures (copy as-is) ----
MOSDEPTH_SRC="$REPO_ROOT/.sync/mosdepth/tests"
for out_dir in "$DST" "$PKG_DST"; do
  cp "$MOSDEPTH_SRC/bad.bed" "$out_dir/bad.bed"
  cp "$MOSDEPTH_SRC/big.bam" "$out_dir/big.bam"
  cp "$MOSDEPTH_SRC/big.bam.csi" "$out_dir/big.bam.csi"
  cp "$MOSDEPTH_SRC/empty-tids.bam" "$out_dir/empty-tids.bam"
  cp "$MOSDEPTH_SRC/empty-tids.bam.bai" "$out_dir/empty-tids.bam.bai"
  cp "$MOSDEPTH_SRC/empty-tids.bed" "$out_dir/empty-tids.bed"
  cp "$MOSDEPTH_SRC/full-fragment-pairs.bam" "$out_dir/full-fragment-pairs.bam"
  cp "$MOSDEPTH_SRC/full-fragment-pairs.bam.bai" "$out_dir/full-fragment-pairs.bam.bai"
  cp "$MOSDEPTH_SRC/missing.bed" "$out_dir/missing.bed"
  cp "$MOSDEPTH_SRC/nanopore.bam" "$out_dir/nanopore.bam"
  cp "$MOSDEPTH_SRC/nanopore.bam.bai" "$out_dir/nanopore.bam.bai"
  cp "$MOSDEPTH_SRC/overlapping-pairs.bam" "$out_dir/overlapping-pairs.bam"
  cp "$MOSDEPTH_SRC/overlapping-pairs.bam.bai" "$out_dir/overlapping-pairs.bam.bai"
  cp "$MOSDEPTH_SRC/ovl.bam" "$out_dir/ovl.bam"
  cp "$MOSDEPTH_SRC/ovl.bam.bai" "$out_dir/ovl.bam.bai"
  cp "$MOSDEPTH_SRC/track.bed" "$out_dir/track.bed"
  cp "$MOSDEPTH_SRC/unordered.bed" "$out_dir/unordered.bed"
done
echo "  mosdepth upstream fixtures + indexes"

# ---- WisecondorX fixed-bin counting fixtures (copy as-is) ----
RWX_SRC="$REPO_ROOT/RWisecondorX/inst/extdata"
for out_dir in "$DST" "$PKG_DST"; do
  cp "$RWX_SRC/fixture_paired.bam" "$out_dir/fixture_paired.bam"
  cp "$RWX_SRC/fixture_paired.bam.bai" "$out_dir/fixture_paired.bam.bai"
  cp "$RWX_SRC/fixture_single.bam" "$out_dir/fixture_single.bam"
  cp "$RWX_SRC/fixture_single.bam.bai" "$out_dir/fixture_single.bam.bai"
  cp "$RWX_SRC/fixture_mixed.bam" "$out_dir/fixture_mixed.bam"
  cp "$RWX_SRC/fixture_mixed.bam.bai" "$out_dir/fixture_mixed.bam.bai"
  cp "$RWX_SRC/fixture_mixed.cram" "$out_dir/fixture_mixed.cram"
  cp "$RWX_SRC/fixture_mixed.cram.crai" "$out_dir/fixture_mixed.cram.crai"
  cp "$RWX_SRC/fixture_ref.fa" "$out_dir/fixture_ref.fa"
  cp "$RWX_SRC/fixture_ref.fa.fai" "$out_dir/fixture_ref.fa.fai"
done
echo "  WisecondorX fixed-bin counting fixtures + indexes"

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

# ---- GFF3 strict-validation fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/gff_strict_valid.gff3" <<'EOF'
##gff-version 3
chr1	src	region	.	.	.	.	.	ID=r
chr1	src	exon	1	10	.	?	.	ID=e;Note=hello%20world
chr1	src	CDS	1	9	.	+	0	ID=cds
EOF
  cat > "$out_dir/gff_strict_invalid.gff3" <<'EOF'
##gff-version 3
chr1	src	exon	1	10	.	+	.	ID=ok
chr1	src	exon	0	10	.	+	.	ID=bad
EOF
  cat > "$out_dir/gff_strict_invalid_attr.gff3" <<'EOF'
##gff-version 3
chr1	src	exon	1	10	.	+	.	ID=ok;broken
EOF
  cat > "$out_dir/gff_strict_extra_field.gff3" <<'EOF'
##gff-version 3
chr1	src	exon	1	10	.	+	.	ID=ok	extra
EOF
  cat > "$out_dir/gff_strict_invalid_end.gff3" <<'EOF'
##gff-version 3
chr1	src	region	.	0	.	.	.	ID=bad_end
EOF
  cat > "$out_dir/gff_attrs.gff3" <<'EOF'
##gff-version 3
chr1	src	gene	1	10	.	+	.	ID=g;Dbxref=GeneID:1,HGNC:HGNC:1;Alias=a;Alias=b;Note=hello%20world
EOF
  cat > "$out_dir/gtf_attrs.gtf" <<'EOF'
chr1	src	exon	1	10	.	+	.	gene_id "G1"; transcript_id "T1"; note "weird; semi";
EOF
done
echo "  GFF/GTF strict-validation and attribute fixtures"

# ---- vcfppR-generated VCF fixtures (spec/mapping/regression) + manifest ----
Rscript "$SCRIPT_DIR/vcfpp.R"
echo "  vcfppR-generated VCF spec/mapping/regression fixtures + manifest"

# ---- Parallel empty-contig VCF regression fixture (bgzip + tabix) ----
for out_dir in "$DST" "$PKG_DST"; do
  bgzip -c "$out_dir/parallel_empty_contigs.vcf" > "$out_dir/parallel_empty_contigs.vcf.gz"
  tabix -f -p vcf "$out_dir/parallel_empty_contigs.vcf.gz"
done
echo "  parallel_empty_contigs.vcf.gz + .tbi"

# ---- bcftools_score multi-summary list / directory / generated-name collision fixtures ----
# List entries intentionally use repo-root-relative paths to match upstream
# bcftools +score --summaries list-file semantics (entries are not resolved
# relative to the list file location).
cp "$DST/score_summary.tsv" "$DST/score_summary_CNT.tsv"
cp "$PKG_DST/score_summary.tsv" "$PKG_DST/score_summary_CNT.tsv"
for out_dir in "$DST" "$PKG_DST"; do
  rm -rf "$out_dir/score_summary_dir"
  mkdir -p "$out_dir/score_summary_dir/nested.tsv"
  cp "$out_dir/score_summary.tsv" "$out_dir/score_summary_dir/score_summary.tsv"
  cp "$out_dir/score_summary_na.tsv" "$out_dir/score_summary_dir/score_summary_na.tsv"
  printf 'sidecar\n' > "$out_dir/score_summary_dir/score_summary.tsv.tbi"
  printf 'sidecar\n' > "$out_dir/score_summary_dir/score_summary_na.tsv.csi"
  printf 'not a summary\n' > "$out_dir/score_summary_dir/README.md"
  printf 'nested ignored\n' > "$out_dir/score_summary_dir/nested.tsv/ignored.txt"
done
cat > "$DST/score_summaries.list" <<'EOF'
test/data/score_summary.tsv
test/data/score_summary_na.tsv
EOF
cat > "$PKG_DST/score_summaries.list" <<'EOF'
score_summary.tsv
score_summary_na.tsv
EOF
echo "  score_summaries.list + score_summary_dir"

echo "==> Done. $(ls "$DST" | wc -l) files in test/data/"
