#!/usr/bin/env bash
# Stage the riker-vs-duckhts WGS benchmark input, exactly as riker's
# benchmark-pipeline does: download the NYGC 1000G high-coverage CRAM for
# HG00188 (ERR3240174) plus the GRCh38 decoy/HLA reference it was aligned to,
# then transcode CRAM -> BAM and index:
#   samtools view -@ N -b -T <ref> -o input.bam source.cram
#   samtools index input.bam
#
# Raw downloads, derived BAM, and provenance stay below the DuckHTS cache.
# Network access is used only here (a staging script), never in the build.
#
# Usage: scripts/stage_riker_wgs_bam.sh [BASE_DIR] [THREADS]
#   BASE_DIR  where to stage (default $DUCKHTS_CACHE_DIR/benchmarks/riker-wgs)
#   THREADS   samtools decode/index threads (default 8)
#
# Rendered report: benchmarks/benchmark_riker_wgs.Rmd
# Driver:          scripts/riker_wgs_benchmark.py
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

BASE="${1:-$(duckhts_cache_subdir "benchmarks/riker-wgs")}"
THREADS="${2:-8}"

REF_URL="https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CRAM_URL="https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240174/HG00188.final.cram"

REF_DIR="$BASE/reference"
DOWNLOAD_DIR="$BASE/downloads"
STAGE="$BASE/stage/HG00188_30x"
REF="$REF_DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CRAM="$DOWNLOAD_DIR/HG00188.final.cram"
BAM="$STAGE/input.bam"
PROVENANCE="$BASE/provenance.tsv"

mkdir -p "$REF_DIR" "$DOWNLOAD_DIR" "$STAGE"

fetch() {  # url dest
  if [[ -s "$2" ]]; then echo "have $(basename "$2")"; return; fi
  echo "downloading $(basename "$2") ..."
  curl -sSfL --retry 5 --retry-delay 10 -C - -o "$2" "$1"
}

fetch "$REF_URL" "$REF"
if [[ ! -s "$REF.fai" ]]; then
  echo "indexing reference ..."
  samtools faidx "$REF"
fi

if [[ ! -s "$BAM" ]]; then
  fetch "$CRAM_URL" "$CRAM"
  samtools quickcheck -v "$CRAM" && echo "CRAM quickcheck OK"

  echo "transcoding CRAM -> BAM ..."
  samtools view -@ "$THREADS" -b -T "$REF" -o "$BAM" "$CRAM"
  samtools quickcheck -v "$BAM" && echo "BAM quickcheck OK"
fi

if [[ ! -s "$BAM.bai" ]]; then
  echo "indexing BAM ..."
  samtools index -@ "$THREADS" "$BAM"
fi

if [[ ! -f "$PROVENANCE" ]]; then
  duckhts_write_provenance "$PROVENANCE" \
    "workload=riker_wgs" \
    "sample=HG00188" \
    "run_accession=ERR3240174" \
    "dataset_release=1000G_2504_high_coverage" \
    "reference_url=$REF_URL" \
    "cram_url=$CRAM_URL" \
    "reference_fasta=$REF" \
    "source_cram=$CRAM" \
    "derived_bam=$BAM" \
    "derived_bai=${BAM}.bai" \
    "staging_command=samtools view -@ $THREADS -b -T $REF -o $BAM $CRAM; samtools index -@ $THREADS $BAM" \
    "samtools=$(samtools --version | head -n 1)"
fi

echo "STAGING COMPLETE: $BAM"
echo "PROVENANCE: $PROVENANCE"
