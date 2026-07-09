#!/usr/bin/env bash
# Stage the riker-vs-duckhts WGS benchmark input, exactly as riker's
# benchmark-pipeline does: download the NYGC 1000G high-coverage CRAM for
# HG00188 (ERR3240174) plus the GRCh38 decoy/HLA reference it was aligned to,
# then transcode CRAM -> BAM and index:
#   samtools view -@ N -b -T <ref> -o input.bam source.cram && rm cram
#   samtools index input.bam
#
# The CRAM is deleted after transcode (riker's "don't double-store" rule).
# Network access is used only here (a staging script), never in the build.
#
# Usage: scripts/stage_riker_wgs_bam.sh [BASE_DIR] [THREADS]
#   BASE_DIR  where to stage (default ./riker-bench-data)
#   THREADS   samtools decode/index threads (default 8)
#
# Rendered report: benchmarks/benchmark_riker_wgs.Rmd
# Driver:          scripts/riker_wgs_benchmark.py
set -euo pipefail

BASE="${1:-./riker-bench-data}"
THREADS="${2:-8}"

REF_URL="https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CRAM_URL="https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240174/HG00188.final.cram"
CRAM_BYTES=14348265497   # exact size; guards against a truncated download

REF_DIR="$BASE/ref"
STAGE="$BASE/stage/HG00188_30x"
REF="$REF_DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CRAM="$STAGE/HG00188.final.cram"
BAM="$STAGE/input.bam"

mkdir -p "$REF_DIR" "$STAGE"

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
  got=$(stat -c %s "$CRAM")
  if [[ "$got" != "$CRAM_BYTES" ]]; then
    echo "FATAL: CRAM size $got != expected $CRAM_BYTES (truncated download)" >&2
    exit 1
  fi
  echo "CRAM size verified: $got bytes"
  samtools quickcheck -v "$CRAM" && echo "CRAM quickcheck OK"

  echo "transcoding CRAM -> BAM ..."
  samtools view -@ "$THREADS" -b -T "$REF" -o "$BAM" "$CRAM"
  samtools quickcheck -v "$BAM" && echo "BAM quickcheck OK"
  rm -f "$CRAM"   # don't double-store
  echo "removed CRAM"
fi

if [[ ! -s "$BAM.bai" ]]; then
  echo "indexing BAM ..."
  samtools index -@ "$THREADS" "$BAM"
fi

echo "STAGING COMPLETE: $BAM"
