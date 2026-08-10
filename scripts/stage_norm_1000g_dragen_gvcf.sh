#!/usr/bin/env bash
set -euo pipefail

# Stage a benchmark-friendly 1000 Genomes DRAGEN gVCF slice beneath the
# DuckHTS cache. Rendering benchmarks remains offline/non-network-dependent;
# run this explicit staging step when the cached workload is absent.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

BCFTOOLS_BIN=${BCFTOOLS_BIN:-bcftools}
TABIX_BIN=${TABIX_BIN:-tabix}
NORM_GVCF_SAMPLE=${NORM_GVCF_SAMPLE:-HG00096}
NORM_GVCF_REGION=${NORM_GVCF_REGION:-chr22:20000000-30000000}
NORM_GVCF_SOURCE=${NORM_GVCF_SOURCE:-https://s3.amazonaws.com/1000genomes-dragen-v3.7.6/data/individuals/hg38-graph-based/${NORM_GVCF_SAMPLE}/${NORM_GVCF_SAMPLE}.hard-filtered.gvcf.gz}
NORM_GVCF_OUT_DIR=${NORM_GVCF_OUT_DIR:-$(duckhts_cache_subdir "benchmarks/norm/1000g-dragen/${NORM_GVCF_SAMPLE}")}
NORM_GVCF_THREADS=${NORM_GVCF_THREADS:-2}

safe_region=${NORM_GVCF_REGION//:/_}
safe_region=${safe_region//-/_}
NORM_GVCF_OUT=${NORM_GVCF_OUT:-${NORM_GVCF_OUT_DIR}/${NORM_GVCF_SAMPLE}.hard-filtered.${safe_region}.g.vcf.gz}
NORM_GVCF_PROVENANCE=${NORM_GVCF_PROVENANCE:-${NORM_GVCF_OUT}.provenance.tsv}

mkdir -p "$(dirname "$NORM_GVCF_OUT")"

if [[ ! -s "$NORM_GVCF_OUT" ]]; then
  "$BCFTOOLS_BIN" view \
    --threads "$NORM_GVCF_THREADS" \
    -r "$NORM_GVCF_REGION" \
    -s "$NORM_GVCF_SAMPLE" \
    -Oz \
    -o "$NORM_GVCF_OUT" \
    "$NORM_GVCF_SOURCE"
else
  echo "Using cached gVCF slice: $NORM_GVCF_OUT"
fi

if [[ ! -s "${NORM_GVCF_OUT}.tbi" ]]; then
  "$TABIX_BIN" -f -p vcf "$NORM_GVCF_OUT"
fi

if [[ ! -f "$NORM_GVCF_PROVENANCE" ]]; then
  duckhts_write_provenance "$NORM_GVCF_PROVENANCE" \
    "workload=normalization_gvcf_slice" \
    "source_url=$NORM_GVCF_SOURCE" \
    "source_release=1000genomes-dragen-v3.7.6" \
    "sample=$NORM_GVCF_SAMPLE" \
    "region=$NORM_GVCF_REGION" \
    "output_vcf=$NORM_GVCF_OUT" \
    "output_index=${NORM_GVCF_OUT}.tbi" \
    "staging_command=bcftools view --threads $NORM_GVCF_THREADS -r $NORM_GVCF_REGION -s $NORM_GVCF_SAMPLE -Oz" \
    "bcftools=$($BCFTOOLS_BIN --version | head -n 1)" \
    "tabix=$($TABIX_BIN --version 2>&1 | head -n 1)"
fi

records=$("$BCFTOOLS_BIN" index -n "$NORM_GVCF_OUT")
samples=$("$BCFTOOLS_BIN" query -l "$NORM_GVCF_OUT" | paste -sd, -)

printf 'Staged 1000G DRAGEN gVCF slice:\n'
printf '  source=%s\n' "$NORM_GVCF_SOURCE"
printf '  sample=%s\n' "$samples"
printf '  region=%s\n' "$NORM_GVCF_REGION"
printf '  records=%s\n' "$records"
printf '  output=%s\n' "$NORM_GVCF_OUT"
printf '  index=%s.tbi\n' "$NORM_GVCF_OUT"
printf '  provenance=%s\n' "$NORM_GVCF_PROVENANCE"
printf '\nUse with benchmark_norm.Rmd:\n'
printf '  NORM_GVCF_VCF=%q NORM_GVCF_PROVENANCE=%q NORM_GVCF_FASTA=/path/to/hg38.fa make bench-norm\n' \
  "$NORM_GVCF_OUT" "$NORM_GVCF_PROVENANCE"
