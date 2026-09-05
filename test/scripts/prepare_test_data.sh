#!/usr/bin/env bash
# Maintainer-only reconstruction of the committed DuckHTS SQL fixture closure.
# A normal fresh clone already has the small fixtures and indexes it needs; use
# `make test_release`, not this script. Full reconstruction also requires the
# local pinned upstream mirrors named below.
#
# Prerequisites: samtools, bcftools, bgzip, tabix (all from htslib/samtools).
#
# Usage:  ./test/scripts/prepare_test_data.sh [--duckvep-only]

set -euo pipefail

MODE="${1:-}"
if [[ -n "$MODE" && "$MODE" != "--duckvep-only" ]]; then
  echo "usage: $0 [--duckvep-only]" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SRC="$REPO_ROOT/third_party/htslib/test"
DST="$REPO_ROOT/test/data"
PKG_DST="$REPO_ROOT/r/Rduckhts/inst/extdata"

mkdir -p "$DST"
mkdir -p "$PKG_DST"

echo "==> Preparing test data in $DST"

# ---- Pinned Ensembl model-build acceptance fixture ----
DUCKVEP_FIXTURE_SRC="$DST/duckvep/ensembl_core"
DUCKVEP_FIXTURE_DST="$PKG_DST/duckvep/ensembl_core"
rm -rf "$DUCKVEP_FIXTURE_DST"
mkdir -p "$DUCKVEP_FIXTURE_DST"
cp -a "$DUCKVEP_FIXTURE_SRC/." "$DUCKVEP_FIXTURE_DST/"
echo "  DuckVEP Ensembl core fixture"

# ---- DuckVEP terminal partial-codon fixtures ----
mkdir -p "$DST/duckvep"
cat > "$DST/duckvep/terminal_partial.fa" <<'EOF'
>partial
AAAAAAAAAAAAAAAAAAAAATGCCCAAAGGGTCAAAAAAAAAAAAAAAAAAAA
EOF
samtools faidx "$DST/duckvep/terminal_partial.fa"
# Exact Ensembl 116 GRCh38 1:45011916-45014703 reference slice used by the
# ENST00000650713 / ClinVar 1:45013701:C:CTAG regression. The committed bases
# are the fixture authority; verify their source digest before rebuilding the
# index or copying the fixture into the R package.
CLINVAR_PARTIAL_FASTA="$DST/duckvep/clinvar_terminal_partial.fa"
CLINVAR_PARTIAL_SHA256="55af53489858b09a2d482d8953320bdf7cecc33d74c1272dfcd1e02fc4d1de0b"
observed_clinvar_partial_sha256="$(
  { tail -n +2 "$CLINVAR_PARTIAL_FASTA" || true; } |
    tr -d '\n' |
    sha256sum |
    cut -d' ' -f1
)"
if [[ "$observed_clinvar_partial_sha256" != "$CLINVAR_PARTIAL_SHA256" ]]; then
  echo "ClinVar terminal-partial reference fixture digest mismatch" >&2
  exit 1
fi
samtools faidx "$CLINVAR_PARTIAL_FASTA"
cp "$DST/duckvep/terminal_partial.fa" "$PKG_DST/duckvep_terminal_partial.fa"
cp "$DST/duckvep/terminal_partial.fa.fai" "$PKG_DST/duckvep_terminal_partial.fa.fai"
cp "$CLINVAR_PARTIAL_FASTA" "$PKG_DST/duckvep_clinvar_terminal_partial.fa"
cp "$CLINVAR_PARTIAL_FASTA.fai" "$PKG_DST/duckvep_clinvar_terminal_partial.fa.fai"
echo "  DuckVEP terminal partial-codon fixtures"

if [[ "$MODE" == "--duckvep-only" ]]; then
  echo "==> DuckVEP fixture sync complete"
  exit 0
fi

# ---- BAM (copy + index) ----
# Synthetic region-union authority is the committed VCF, including its two
# identical physical records and REF/END spans across disjoint requests.
bgzip -c "$DST/region_union.vcf" > "$DST/region_union.vcf.gz"
tabix -f -p vcf "$DST/region_union.vcf.gz"
bcftools view --no-version -Ob -o "$DST/region_union.bcf" "$DST/region_union.vcf"
bcftools index -f "$DST/region_union.bcf"
cp "$DST/region_union.vcf" "$DST/region_union.vcf.gz" "$DST/region_union.vcf.gz.tbi" \
   "$DST/region_union.bcf" "$DST/region_union.bcf.csi" "$PKG_DST/"

(cd "$REPO_ROOT" && Rscript test/scripts/prepare_bam_scan_fixtures.R)
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

# ---- Malformed VCF regression fixture ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/malformed_bad_pos.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	10	.	A	C	.	PASS	.	GT	0/1
1	bad	.	A	G	.	PASS	.	GT	0/1
EOF
done
echo "  malformed_bad_pos.vcf"

# ---- Corrupt BCF header-vs-payload type-clash regression fixtures ----
TMP_BCF_CLASH_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_BCF_CLASH_DIR"' EXIT
cat > "$TMP_BCF_CLASH_DIR/info_string.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=String,Description="depth encoded as string">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	C	.	PASS	DP=abc
EOF
bcftools view --no-version -Ob -o "$TMP_BCF_CLASH_DIR/info_string.bcf" "$TMP_BCF_CLASH_DIR/info_string.vcf"
bcftools view --no-version -h "$TMP_BCF_CLASH_DIR/info_string.bcf" | perl -pe 's/##INFO=<ID=DP,Number=1,Type=String,Description="depth encoded as string">/##INFO=<ID=DP,Number=1,Type=Integer,Description="depth claimed as integer">/' > "$TMP_BCF_CLASH_DIR/info_int.hdr"
bcftools reheader -h "$TMP_BCF_CLASH_DIR/info_int.hdr" -o "$TMP_BCF_CLASH_DIR/bcf_info_type_clash.bcf" "$TMP_BCF_CLASH_DIR/info_string.bcf"
cat > "$TMP_BCF_CLASH_DIR/format_string.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=XX,Number=1,Type=String,Description="format encoded as string">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	10	.	A	C	.	PASS	.	XX	foo	bar
EOF
bcftools view --no-version -Ob -o "$TMP_BCF_CLASH_DIR/format_string.bcf" "$TMP_BCF_CLASH_DIR/format_string.vcf"
bcftools view --no-version -h "$TMP_BCF_CLASH_DIR/format_string.bcf" | perl -pe 's/##FORMAT=<ID=XX,Number=1,Type=String,Description="format encoded as string">/##FORMAT=<ID=XX,Number=1,Type=Integer,Description="format claimed as integer">/' > "$TMP_BCF_CLASH_DIR/format_int.hdr"
bcftools reheader -h "$TMP_BCF_CLASH_DIR/format_int.hdr" -o "$TMP_BCF_CLASH_DIR/bcf_format_type_clash.bcf" "$TMP_BCF_CLASH_DIR/format_string.bcf"
# Reverse clash: numeric payload with the header reheadered to claim String.
# htslib's string decode path copies the payload without checking info->type,
# so unless the INFO string preflight runs under every policy this returns raw
# bytes instead of NULL under the default 'null' policy.
cat > "$TMP_BCF_CLASH_DIR/info_int.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##INFO=<ID=NN,Number=1,Type=Integer,Description="value encoded as integer">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	C	.	PASS	NN=42
EOF
bcftools view --no-version -Ob -o "$TMP_BCF_CLASH_DIR/info_int_payload.bcf" "$TMP_BCF_CLASH_DIR/info_int.vcf"
bcftools view --no-version -h "$TMP_BCF_CLASH_DIR/info_int_payload.bcf" | perl -pe 's/##INFO=<ID=NN,Number=1,Type=Integer,Description="value encoded as integer">/##INFO=<ID=NN,Number=1,Type=String,Description="value claimed as string">/' > "$TMP_BCF_CLASH_DIR/info_str.hdr"
bcftools reheader -h "$TMP_BCF_CLASH_DIR/info_str.hdr" -o "$TMP_BCF_CLASH_DIR/bcf_info_str_clash.bcf" "$TMP_BCF_CLASH_DIR/info_int_payload.bcf"
for out_dir in "$DST" "$PKG_DST"; do
  cp "$TMP_BCF_CLASH_DIR/bcf_info_type_clash.bcf" "$out_dir/bcf_info_type_clash.bcf"
  cp "$TMP_BCF_CLASH_DIR/bcf_format_type_clash.bcf" "$out_dir/bcf_format_type_clash.bcf"
  cp "$TMP_BCF_CLASH_DIR/bcf_info_str_clash.bcf" "$out_dir/bcf_info_str_clash.bcf"
done
rm -rf "$TMP_BCF_CLASH_DIR"
trap - EXIT
echo "  bcf_info_type_clash.bcf + bcf_format_type_clash.bcf + bcf_info_str_clash.bcf"

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

# ---- FASTQ query-name projection regression ----
for out_dir in "$DST" "$PKG_DST"; do
  cp "$SCRIPT_DIR/fastq_long_qname.fq" "$out_dir/fastq_long_qname.fq"
done
echo "  fastq_long_qname.fq"

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

# ---- Liftover repeat-run indel swap regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/liftover_repeat_src.fa" <<'EOF'
>chrS
GTTTTCT
EOF
  cat > "$out_dir/liftover_repeat_dst.fa" <<'EOF'
>chrD
GTTTTTCT
EOF
  cat > "$out_dir/liftover_repeat.chain" <<'EOF'
chain 100 chrS 7 + 0 7 chrD 8 + 0 8 1
1 0 1
6
EOF
  samtools faidx "$out_dir/liftover_repeat_src.fa"
  samtools faidx "$out_dir/liftover_repeat_dst.fa"
done
echo "  liftover repeat-run indel swap fixtures + .fai"

# ---- Liftover clip-pad negative-score regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/liftover_clip_pad_src.fa" <<'EOF'
>chrS
TCCGTCTCAAAAA
EOF
  cat > "$out_dir/liftover_clip_pad_dst.fa" <<'EOF'
>chrD
TCCGTCTTAAAAA
EOF
  cat > "$out_dir/liftover_clip_pad.chain" <<'EOF'
chain 100 chrS 13 + 0 13 chrD 13 + 0 12 1
8 1 0
4
EOF
  samtools faidx "$out_dir/liftover_clip_pad_src.fa"
  samtools faidx "$out_dir/liftover_clip_pad_dst.fa"
done
echo "  liftover clip-pad negative-score fixtures + .fai"

# ---- Liftover clip-pad allocation-limit regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  python3 - "$out_dir" <<'PY'
from pathlib import Path
import sys

out_dir = Path(sys.argv[1])
first_block = 2201
source = ("ACGT" * 1100)[:4005]
destination = source[:first_block] + source[first_block + 1:]
(out_dir / "liftover_nw_limit_src.fa").write_text(">chrS\n" + source + "\n")
(out_dir / "liftover_nw_limit_dst.fa").write_text(">chrD\n" + destination + "\n")
(out_dir / "liftover_nw_limit.chain").write_text(
    "chain 100 chrS 4005 + 0 4005 chrD 4004 + 0 4004 1\n"
    "2201 1 0\n"
    "1803\n"
)
PY
  samtools faidx "$out_dir/liftover_nw_limit_src.fa"
  samtools faidx "$out_dir/liftover_nw_limit_dst.fa"
done
echo "  liftover clip-pad allocation-limit fixtures + .fai"

# ---- Liftover chr23/X source-FASTA alias regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/liftover_chr23_alias_src.fa" <<'EOF'
>chrX
ACGTACGTAA
EOF
  cat > "$out_dir/liftover_chr23_alias_dst.fa" <<'EOF'
>chrLiftX
ACGTACGTAA
EOF
  cat > "$out_dir/liftover_chr23_alias.chain" <<'EOF'
chain 100 X 10 + 0 10 chrLiftX 10 + 0 10 1
10
EOF
  samtools faidx "$out_dir/liftover_chr23_alias_src.fa"
  samtools faidx "$out_dir/liftover_chr23_alias_dst.fa"
done
echo "  liftover chr23/X source-FASTA alias fixtures + .fai"

# ---- Liftover spanning-deletion swap / ref-add regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/liftover_star_swap_src.fa" <<'EOF'
>chrS
TCCAC
EOF
  cat > "$out_dir/liftover_star_swap_dst.fa" <<'EOF'
>chrD
TCTGC
EOF
  cat > "$out_dir/liftover_star_swap.chain" <<'EOF'
chain 1 chrS 5 + 0 5 chrD 5 + 0 5 1
5
EOF
  cat > "$out_dir/liftover_star_refadd_src.fa" <<'EOF'
>chrS
GTGCGTGGGTGGGC
EOF
  cat > "$out_dir/liftover_star_refadd_dst.fa" <<'EOF'
>chrD
GTGCGGCCGGGGGGGC
EOF
  cat > "$out_dir/liftover_star_refadd.chain" <<'EOF'
chain 1 chrS 14 + 0 14 chrD 16 + 0 16 1
5 0 2
9
EOF
  samtools faidx "$out_dir/liftover_star_swap_src.fa"
  samtools faidx "$out_dir/liftover_star_swap_dst.fa"
  samtools faidx "$out_dir/liftover_star_refadd_src.fa"
  samtools faidx "$out_dir/liftover_star_refadd_dst.fa"
done
echo "  liftover spanning-deletion swap / ref-add fixtures + .fai"

# ---- vcfppR-generated VCF fixtures (spec/mapping/regression) + manifest ----
Rscript "$SCRIPT_DIR/vcfpp.R"
echo "  vcfppR-generated VCF spec/mapping/regression fixtures + manifest"

# ---- read_bcf regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  cat > "$out_dir/genotype_ploidy_edge_cases.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BIG	DIP	HAP	TET	TRI
chr1	1	ploidy_edges	A	C,G,T,AA,AC,AG,AT,CA,CC,CG	.	PASS	.	GT	0/10	0|1	1	0/1/2/3	0/1/2
EOF

  cat > "$out_dir/test_vep_tidy.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL|DISTANCE">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	vep_tidy	A	C	.	PASS	CSQ=C|missense_variant|GENE1|12	GT	0/1	1/1
EOF
done

Rscript --vanilla - "$DST/tidy_chunk_boundary.vcf" "$PKG_DST/tidy_chunk_boundary.vcf" <<'RS'
args <- commandArgs(TRUE)
samples <- sprintf("S%04d", seq_len(2053))
lines <- c(
  "##fileformat=VCFv4.3",
  "##contig=<ID=chr1,length=1000>",
  "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL|DISTANCE\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"),
  paste(c("chr1", "10", "tidy_chunk", "A", "C", ".", "PASS", "CSQ=C|missense_variant|GENE1|12", "GT", rep("0/1", length(samples))), collapse = "\t")
)
for (path in args) writeLines(lines, path)
RS

echo "  read_bcf projection/ploidy regression fixtures"

# ---- DuckVEP / bcftools csq fixture plumbing ----
mkdir -p "$DST/duckvep"
cat > "$DST/duckvep/minimal.fa" <<'EOF'
>chrDuck
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
cat > "$DST/duckvep/minimal.gff3" <<'EOF'
##gff-version 3
chrDuck	duckvep	gene	100	300	.	+	.	ID=gene:DUCK1;Name=DUCK1;biotype=protein_coding
chrDuck	duckvep	mRNA	100	300	.	+	.	ID=transcript:DUCK1-201;Parent=gene:DUCK1;Name=DUCK1-201;biotype=protein_coding
chrDuck	duckvep	exon	100	150	.	+	.	ID=exon:DUCK1-1;Parent=transcript:DUCK1-201;rank=1
chrDuck	duckvep	CDS	120	240	.	+	0	ID=CDS:DUCK1;Parent=transcript:DUCK1-201
chrDuck	duckvep	exon	200	300	.	+	.	ID=exon:DUCK1-2;Parent=transcript:DUCK1-201;rank=2
EOF
cat > "$DST/duckvep/minimal_bcsq.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chrDuck,length=256>
##INFO=<ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq fixture. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrDuck	125	duck_missense	A	C	.	PASS	BCSQ=missense_variant|DUCK1|DUCK1-201|protein_coding|+|M1T|125A>C
chrDuck	160	duck_intron	C	T	.	PASS	BCSQ=intron_variant|DUCK1|DUCK1-201|protein_coding|+||160C>T
chrDuck	350	duck_intergenic	G	A	.	PASS	BCSQ=intergenic_variant||||||350G>A
EOF
cat > "$DST/duckvep/ensembl_release_consequences.vcf" <<'EOF'
##fileformat=VCFv4.2
##source=ensembl;version=116;url=https://e116.ensembl.org/homo_sapiens
##INFO=<ID=VE,Number=.,Type=String,Description="Variant effect of a variant overlapping a sequence feature as computed by the ensembl variant effect pipeline. Format=Consequence|Index|Feature_type|Feature_id. Index indentifies for which variant sequence the effect is described for.">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl's Variant Effect Pipeline. Format=Allele|Consequence|Feature_type|Feature|Amino_acids|SIFT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	15528196	rs1312204123	A	G	.	.	VE=missense_variant|0|mRNA|ENST00000643195,non_coding_transcript_variant|0|ncRNA|ENST00000775238,intron_variant|0|primary_transcript|ENST00000775238;CSQ=G|intron_variant|primary_transcript|ENST00000775238||,G|missense_variant|mRNA|ENST00000643195|N/S|deleterious_-_low_confidence(0.05)
22	15528198	rs1985999326	G	A,T	.	.	VE=missense_variant|0|mRNA|ENST00000643195,missense_variant|1|mRNA|ENST00000643195,non_coding_transcript_variant|0|ncRNA|ENST00000775238,intron_variant|0|primary_transcript|ENST00000775238,non_coding_transcript_variant|1|ncRNA|ENST00000775238,intron_variant|1|primary_transcript|ENST00000775238;CSQ=A|missense_variant|mRNA|ENST00000643195|V/I|tolerated_-_low_confidence(1),A|intron_variant|primary_transcript|ENST00000775238||,T|intron_variant|primary_transcript|ENST00000775238||,T|missense_variant|mRNA|ENST00000643195|V/F|tolerated_-_low_confidence(0.16)
EOF
cp "$DST/duckvep/ensembl_release_consequences.vcf" \
  "$PKG_DST/ensembl_release_consequences.vcf"
echo "  DuckVEP/csq minimal and Ensembl release consequence fixtures"

# ---- Parallel empty-contig VCF/BCF regression fixtures ----
for out_dir in "$DST" "$PKG_DST"; do
  bgzip -c "$out_dir/parallel_empty_contigs.vcf" > "$out_dir/parallel_empty_contigs.vcf.gz"
  tabix -f -p vcf "$out_dir/parallel_empty_contigs.vcf.gz"
  bcftools view --no-version -Ob -o "$out_dir/parallel_empty_contigs.bcf" \
    "$out_dir/parallel_empty_contigs.vcf"
  bcftools index -f "$out_dir/parallel_empty_contigs.bcf"
done
echo "  parallel_empty_contigs.vcf.gz + .tbi and .bcf + .csi"

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
