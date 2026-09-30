#!/usr/bin/env bash
# Prints the `bcftools roh` RG segments for every scenario in test/sql/roh.test, as SQL
# VALUES rows (sample, chrom, start, end, length, n_markers, quality). The expected rows
# in the test were produced by this script with bcftools 1.23.1-70-g6dbd8fef.
# Run from the repository root: bash test/scripts/roh_bcftools_expected.sh
set -euo pipefail

vcf=test/data/roh_fixture.vcf.gz
af_file=test/data/roh_af.tsv.gz
map_mask='test/data/roh_map_{CHROM}.txt'
tmp=$(mktemp -d)
trap 'rm -r "$tmp"' EXIT
cp test/data/roh_map_chr1.txt "$tmp/chr1_only_chr1.txt"

scenario() {
  local name=$1
  shift
  echo "-- $name: bcftools roh $*"
  bcftools roh "$@" -O r "$vcf" 2>/dev/null | awk -F'\t' '
    $1 == "RG" { printf "%s('\''%s'\'', '\''%s'\'', %s, %s, %s, %s, %s)", sep, $2, $3, $4, $5, $6, $7, $8; sep = ",\n " }
    END { print ";" }' | sed '1s/^/ /'
}

scenario pl_af_tag --AF-tag AF
scenario gt30_af_tag -G30 --AF-tag AF
scenario pl_af_file --AF-file "$af_file"
scenario rec_rate --AF-tag AF -M 1e-6
scenario genetic_map --AF-tag AF -m "$map_mask"
scenario genetic_map_rec_rate --AF-tag AF -m "$map_mask" -M 0.5
scenario genetic_map_chr1_only --AF-tag AF -m "$tmp/chr1_only_{CHROM}.txt"
scenario transitions --AF-tag AF -a 1e-5 -H 1e-6
scenario gt30_af_file_map -G30 --AF-file "$af_file" -m "$map_mask"
scenario samples --AF-tag AF -s S2,S4
