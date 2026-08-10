#!/usr/bin/env bash
# Stage the GRCh37/GRCh38 liftover references below the DuckHTS cache.
# Network access is confined to this explicit command; benchmarks and
# conformance runners only consume the staged paths it prints.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

SAMTOOLS_BIN=${SAMTOOLS_BIN:-samtools}
LIFTOVER_REFERENCE_DIR=${LIFTOVER_REFERENCE_DIR:-$(duckhts_cache_subdir "references/liftover/grch37-to-grch38")}
LIFTOVER_SRC_FASTA_URL=${LIFTOVER_SRC_FASTA_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz}
LIFTOVER_DST_FASTA_URL=${LIFTOVER_DST_FASTA_URL:-https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta}
LIFTOVER_CHAIN_URL=${LIFTOVER_CHAIN_URL:-https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz}

for tool in curl gzip "$SAMTOOLS_BIN"; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "required command not found: $tool" >&2
        exit 1
    }
done

mkdir -p "$LIFTOVER_REFERENCE_DIR/downloads"
SRC_GZ="$LIFTOVER_REFERENCE_DIR/downloads/human_g1k_v37.fasta.gz"
SRC_FASTA="$LIFTOVER_REFERENCE_DIR/human_g1k_v37.fasta"
DST_FASTA="$LIFTOVER_REFERENCE_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHAIN="$LIFTOVER_REFERENCE_DIR/GRCh37_to_GRCh38.chain.gz"
PROVENANCE="$LIFTOVER_REFERENCE_DIR/provenance.tsv"

fetch() { # URL DESTINATION
    if [[ -s "$2" ]]; then
        echo "Using cached $(basename "$2")"
        return
    fi
    echo "Downloading $1"
    curl --fail --location --retry 5 --retry-delay 5 --continue-at - --output "$2" "$1"
}

fetch "$LIFTOVER_SRC_FASTA_URL" "$SRC_GZ"
if [[ ! -s "$SRC_FASTA" ]]; then
    gzip --decompress --keep --stdout "$SRC_GZ" >"${SRC_FASTA}.partial.$$"
    mv -f "${SRC_FASTA}.partial.$$" "$SRC_FASTA"
fi
fetch "$LIFTOVER_DST_FASTA_URL" "$DST_FASTA"
fetch "$LIFTOVER_CHAIN_URL" "$CHAIN"

for fasta in "$SRC_FASTA" "$DST_FASTA"; do
    if [[ ! -s "${fasta}.fai" ]]; then
        "$SAMTOOLS_BIN" faidx "$fasta"
    fi
done

if [[ ! -f "$PROVENANCE" ]]; then
    duckhts_write_provenance "$PROVENANCE" \
        "workload=liftover_reference_bundle" \
        "source_assembly=GRCh37_human_g1k_v37" \
        "source_fasta_release=1000_Genomes_phase2_reference" \
        "source_fasta_url=$LIFTOVER_SRC_FASTA_URL" \
        "source_fasta_download=$SRC_GZ" \
        "source_fasta=$SRC_FASTA" \
        "destination_assembly=GRCh38_GCA_000001405.15_no_alt_analysis_set" \
        "destination_fasta_release=Broad_hg38_v0" \
        "destination_fasta_url=$LIFTOVER_DST_FASTA_URL" \
        "destination_fasta=$DST_FASTA" \
        "chain_source=Ensembl_GRCh37_to_GRCh38_assembly_mapping" \
        "chain_url=$LIFTOVER_CHAIN_URL" \
        "chain=$CHAIN" \
        "staging_command=gzip --decompress; samtools faidx" \
        "samtools=$($SAMTOOLS_BIN --version | head -n 1)"
fi

cat <<EOF
Liftover references staged under: $LIFTOVER_REFERENCE_DIR
Provenance: $PROVENANCE

Use with benchmark or conformance commands:
  LIFTOVER_CHAIN_PATH=$CHAIN
  LIFTOVER_SRC_FASTA=$SRC_FASTA
  LIFTOVER_DST_FASTA=$DST_FASTA
EOF
