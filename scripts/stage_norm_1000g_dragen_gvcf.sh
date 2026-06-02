#!/usr/bin/env bash
set -euo pipefail

# Stage a local, benchmark-friendly 1000 Genomes DRAGEN gVCF slice.
# Rendering benchmarks remains offline/non-network-dependent; run this script
# explicitly when refreshing the optional public gVCF workload.

BCFTOOLS_BIN=${BCFTOOLS_BIN:-bcftools}
TABIX_BIN=${TABIX_BIN:-tabix}
NORM_GVCF_SAMPLE=${NORM_GVCF_SAMPLE:-HG00096}
NORM_GVCF_REGION=${NORM_GVCF_REGION:-chr22:20000000-30000000}
NORM_GVCF_SOURCE=${NORM_GVCF_SOURCE:-https://s3.amazonaws.com/1000genomes-dragen-v3.7.6/data/individuals/hg38-graph-based/${NORM_GVCF_SAMPLE}/${NORM_GVCF_SAMPLE}.hard-filtered.gvcf.gz}
NORM_GVCF_OUT_DIR=${NORM_GVCF_OUT_DIR:-.tmp/1000g_dragen_hg00096}
NORM_GVCF_THREADS=${NORM_GVCF_THREADS:-2}

safe_region=${NORM_GVCF_REGION//:/_}
safe_region=${safe_region//-/_}
NORM_GVCF_OUT=${NORM_GVCF_OUT:-${NORM_GVCF_OUT_DIR}/${NORM_GVCF_SAMPLE}.hard-filtered.${safe_region}.g.vcf.gz}

mkdir -p "$(dirname "$NORM_GVCF_OUT")"

"$BCFTOOLS_BIN" view \
  --threads "$NORM_GVCF_THREADS" \
  -r "$NORM_GVCF_REGION" \
  -s "$NORM_GVCF_SAMPLE" \
  -Oz \
  -o "$NORM_GVCF_OUT" \
  "$NORM_GVCF_SOURCE"

"$TABIX_BIN" -f -p vcf "$NORM_GVCF_OUT"

records=$("$BCFTOOLS_BIN" index -n "$NORM_GVCF_OUT")
samples=$("$BCFTOOLS_BIN" query -l "$NORM_GVCF_OUT" | paste -sd, -)

printf 'Staged 1000G DRAGEN gVCF slice:\n'
printf '  source=%s\n' "$NORM_GVCF_SOURCE"
printf '  sample=%s\n' "$samples"
printf '  region=%s\n' "$NORM_GVCF_REGION"
printf '  records=%s\n' "$records"
printf '  output=%s\n' "$NORM_GVCF_OUT"
printf '  index=%s.tbi\n' "$NORM_GVCF_OUT"
printf '\nUse with benchmark_norm.Rmd:\n'
printf '  NORM_GVCF_VCF=%q NORM_GVCF_FASTA=/path/to/hg38.fa make bench-norm\n' "$NORM_GVCF_OUT"
