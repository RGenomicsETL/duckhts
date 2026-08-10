#!/usr/bin/env bash
# Stage the pinned DuckBedQC BED providers used by the real cgranges benchmark.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

DUCKBEDQC_REPOSITORY=${DUCKBEDQC_REPOSITORY:-https://github.com/sounkou-bioinfo/DuckBedQC.git}
DUCKBEDQC_REF=${DUCKBEDQC_REF:-118fc21c6cde9d680989dd4d1b613789539469f3}
DUCKBEDQC_DIR=${DUCKBEDQC_DIR:-$(duckhts_cache_subdir datasets/duckbedqc)}
PROVENANCE="$DUCKBEDQC_DIR/provenance.tsv"

command -v git >/dev/null 2>&1 || { echo "git is required" >&2; exit 1; }
if [[ ! -d "$DUCKBEDQC_DIR/.git" ]]; then
    mkdir -p "$(dirname "$DUCKBEDQC_DIR")"
    git clone "$DUCKBEDQC_REPOSITORY" "$DUCKBEDQC_DIR"
fi
git -C "$DUCKBEDQC_DIR" fetch --quiet origin "$DUCKBEDQC_REF"
git -C "$DUCKBEDQC_DIR" checkout --quiet --detach "$DUCKBEDQC_REF"
if [[ "$(git -C "$DUCKBEDQC_DIR" rev-parse HEAD)" != "$DUCKBEDQC_REF" ]]; then
    echo "DuckBedQC checkout does not match requested ref: $DUCKBEDQC_REF" >&2
    exit 1
fi
if [[ ! -f "$DUCKBEDQC_DIR/data/GRCh38_exons.bed" ||
      ! -f "$DUCKBEDQC_DIR/data/GRCh38_illumina_clinical_regions_v100.39.0.bed" ]]; then
    echo "Pinned DuckBedQC ref lacks the required GRCh38 cgranges providers" >&2
    exit 1
fi
duckhts_write_provenance "$PROVENANCE" \
        "workload=cgranges_real_benchmark" \
        "source_repository=$DUCKBEDQC_REPOSITORY" \
        "source_revision=$DUCKBEDQC_REF" \
        "checkout=$DUCKBEDQC_DIR" \
        "subject_bed=$DUCKBEDQC_DIR/data/GRCh38_exons.bed" \
        "query_bed=$DUCKBEDQC_DIR/data/GRCh38_illumina_clinical_regions_v100.39.0.bed" \
        "staging_command=git clone; git fetch; git checkout --detach"
printf 'DUCKBEDQC_DIR=%s\n' "$DUCKBEDQC_DIR"
printf 'DUCKBEDQC_PROVENANCE=%s\n' "$PROVENANCE"
