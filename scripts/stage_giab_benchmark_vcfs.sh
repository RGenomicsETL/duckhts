#!/usr/bin/env bash
# Stage public GIAB NIST v4.2.1 GRCh37 benchmark VCFs into the DuckHTS cache.
# The source table is the authority for release URLs and sample identities.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

BCFTOOLS_BIN=${BCFTOOLS_BIN:-bcftools}
CASES_TSV=${LIFTOVER_CASES_TSV:-$SCRIPT_DIR/conformance_case_table.tsv}
GIAB_SAMPLES=${GIAB_SAMPLES:-HG001,HG002,HG006}
GIAB_DIR=${GIAB_DIR:-$(duckhts_cache_subdir datasets/giab/nist-v4.2.1/grch37)}

command -v curl >/dev/null 2>&1 || { echo "curl is required" >&2; exit 1; }
command -v "$BCFTOOLS_BIN" >/dev/null 2>&1 || {
    echo "bcftools not found: $BCFTOOLS_BIN" >&2
    exit 1
}
[[ -f "$CASES_TSV" ]] || { echo "case table not found: $CASES_TSV" >&2; exit 1; }
mkdir -p "$GIAB_DIR"

sample_requested() {
    local sample="$1" requested
    IFS=',' read -r -a requested <<<"$GIAB_SAMPLES"
    for candidate in "${requested[@]}"; do
        [[ "$sample" == "${candidate//[[:space:]]/}" ]] && return 0
    done
    return 1
}

stage_one() { # case_id sample source URL
    local case_id="$1" sample="$2" source="$3"
    local output provenance
    output="$GIAB_DIR/$(basename "${source%%\?*}")"
    provenance="${output}.provenance.tsv"

    if [[ ! -s "$output" ]]; then
        echo "Downloading $sample NIST v4.2.1 GRCh37 benchmark VCF"
        curl --fail --location --retry 5 --retry-delay 5 --continue-at - --output "$output" "$source"
    else
        echo "Using cached $sample benchmark VCF: $output"
    fi
    if [[ ! -s "${output}.tbi" ]]; then
        "$BCFTOOLS_BIN" index -f -t "$output"
    fi
    if [[ ! -f "$provenance" ]]; then
        duckhts_write_provenance "$provenance" \
            "workload=giab_nist_v4.2.1_grch37" \
            "case_id=$case_id" \
            "sample=$sample" \
            "source_url=$source" \
            "release=NISTv4.2.1" \
            "assembly=GRCh37" \
            "cached_vcf=$output" \
            "cached_index=${output}.tbi" \
            "staging_command=curl; bcftools index -t" \
            "bcftools=$($BCFTOOLS_BIN --version | head -n 1)"
    fi
}

count=0
while IFS=$'\t' read -r case_id dataset sample description default_region source; do
    [[ "$case_id" == "case_id" || -z "$case_id" ]] && continue
    [[ "$dataset" == "giab" ]] || continue
    sample_requested "$sample" || continue
    stage_one "$case_id" "$sample" "$source"
    count=$((count + 1))
done < "$CASES_TSV"

if [[ "$count" -eq 0 ]]; then
    echo "no GIAB samples selected from $CASES_TSV: $GIAB_SAMPLES" >&2
    exit 1
fi

echo "GIAB VCF cache: $GIAB_DIR"
