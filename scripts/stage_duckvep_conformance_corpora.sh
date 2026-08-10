#!/usr/bin/env bash
# Stage external GRCh38 consequence corpora without downloading their complete
# cohort callsets. Network access is confined to this explicit staging step.
# Every derived VCF is deterministic, sites-only, indexed with tabix, and
# checksum-locked; a changed upstream object or tool output fails closed.
#
# Usage:
#   scripts/stage_duckvep_conformance_corpora.sh [OUTPUT_DIR] [CORPUS]
#
# OUTPUT_DIR defaults to $DUCKHTS_CACHE_DIR/corpora/duckvep. CORPUS is one
# of: all (default), hprc-african4-chr22, sniffles2-chr22, dbvar-chr22.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

OUTPUT_DIR="${1:-$(duckhts_cache_subdir "corpora/duckvep")}"
CORPUS="${2:-all}"

for tool in bcftools tabix sha256sum curl; do
	if ! command -v "$tool" >/dev/null 2>&1; then
		echo "FATAL: required tool is not installed: $tool" >&2
		exit 1
	fi
done

sha256_of() {
	sha256sum "$1" | awk '{print $1}'
}

verify_output() { # path expected_sha256 expected_bytes
	local path="$1"
	local expected_sha256="$2"
	local expected_bytes="$3"
	local actual_sha256 actual_bytes

	if [[ ! -s "$path" ]]; then
		return 1
	fi
	actual_sha256="$(sha256_of "$path")"
	actual_bytes="$(stat -c %s "$path")"
	if [[ "$actual_sha256" != "$expected_sha256" ]]; then
		echo "FATAL: $path SHA-256 $actual_sha256 != $expected_sha256" >&2
		exit 1
	fi
	if [[ "$actual_bytes" != "$expected_bytes" ]]; then
		echo "FATAL: $path bytes $actual_bytes != $expected_bytes" >&2
		exit 1
	fi
	return 0
}

finish_vcf() { # partial output expected VCF sha/bytes expected TBI sha/bytes
	local partial="$1"
	local output="$2"
	local vcf_sha="$3"
	local vcf_bytes="$4"
	local tbi_sha="$5"
	local tbi_bytes="$6"

	tabix -f -p vcf "$partial"
	verify_output "$partial" "$vcf_sha" "$vcf_bytes"
	verify_output "$partial.tbi" "$tbi_sha" "$tbi_bytes"
	mv "$partial" "$output"
	mv "$partial.tbi" "$output.tbi"
}

write_receipt() { # output source source_id
	local output="$1"
	local source="$2"
	local source_id="$3"
	local receipt="$output.receipt.tsv"

	printf 'field\tvalue\n' >"$receipt"
	printf 'source_url\t%s\n' "$source" >>"$receipt"
	printf 'source_id\t%s\n' "$source_id" >>"$receipt"
	printf 'output_path\t%s\n' "$output" >>"$receipt"
	printf 'output_bytes\t%s\n' "$(stat -c %s "$output")" >>"$receipt"
	printf 'output_sha256\t%s\n' "$(sha256_of "$output")" >>"$receipt"
	printf 'index_sha256\t%s\n' "$(sha256_of "$output.tbi")" >>"$receipt"
	printf 'bcftools\t%s\n' "$(bcftools --version | head -n 1)" >>"$receipt"
	printf 'htslib\t%s\n' "$(bcftools --version | sed -n '2p')" >>"$receipt"
}

stage_hprc() {
	local dir="$OUTPUT_DIR/hprc_v2_grch38"
	local output="$dir/hprc_v2_grch38_chr22_african4_carried.vcf.gz"
	local partial="$output.partial.$$"
	local base="https://human-pangenomics.s3.us-west-2.amazonaws.com/pangenomes/scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-grch38"
	local vcf_version="LnGRAq5AlCslA2a31J72oUZMxvR_5yuI"
	local tbi_version="FuarSRln2LCGzuG.waMWXsSLQ9UoI5nQ"
	local vcf_url="$base/hprc-v2.0-mc-grch38.vcf.gz?versionId=$vcf_version"
	local tbi_url="$base/hprc-v2.0-mc-grch38.vcf.gz.tbi?versionId=$tbi_version"
	local source="${vcf_url}##idx##${tbi_url}"
	local samples="HG02055,HG02145,HG02723,HG03098"
	local vcf_sha="c7210035d99961e6243e432ffe513ec18520cd097af2b09cd945e5a2674aa309"
	local tbi_sha="8c6870f3e5bed25666812376dfea6cce9fdff8f8bc9aa6703c10bf6165ae1130"

	mkdir -p "$dir"
	if verify_output "$output" "$vcf_sha" 9321604 &&
		verify_output "$output.tbi" "$tbi_sha" 18525; then
		echo "have $output"
	else
		bcftools view --no-version -r chr22 -s "$samples" -Ou "$source" |
			bcftools view --no-version -i 'GT="alt"' -Ou |
			bcftools norm --no-version -m -any -Ou |
			bcftools view --no-version -i 'GT="alt"' -Ou |
			bcftools view --no-version -G -Ou |
			bcftools annotate --no-version \
				--rename-chrs <(printf 'chr22\t22\n') \
				-Oz -o "$partial"
		finish_vcf "$partial" "$output" "$vcf_sha" 9321604 \
			"$tbi_sha" 18525
	fi
	write_receipt "$output" "$vcf_url" \
		"vcf_version=$vcf_version;tbi_version=$tbi_version;region=chr22;samples=$samples;carried_alt_only=true;split_multiallelic=true;genotypes_removed=true;chrom_alias=chr22:22"
}

stage_sniffles2() {
	local dir="$OUTPUT_DIR/sniffles2_1kgp"
	local output="$dir/sniffles2_joint_chr22_sites.vcf.gz"
	local partial="$output.partial.$$"
	local source="https://1kgp-sv-imputation.s3.eu-west-1.amazonaws.com/sv_calls/sniffles2_joint_sv_calls.vcf.gz"
	local source_etag='"131c622662ade1c745a7fad0b3b40be7-183"'
	local vcf_sha="814ff6ca8556a61ae3140c6f47cadd73b8f9ceb131b1dc7fc7ccadae7dfd9e32"
	local tbi_sha="0aaf1f0f52bed48d0dbeaabbc5944da72f717c52f7a838f878b98d786283ae00"
	local actual_etag

	mkdir -p "$dir"
	actual_etag="$(curl -fsSIL "$source" | awk 'tolower($1) == "etag:" {gsub("\\r", "", $2); print $2}' | tail -n 1)"
	if [[ "$actual_etag" != "$source_etag" ]]; then
		echo "FATAL: Sniffles2 source ETag $actual_etag != $source_etag" >&2
		exit 1
	fi
	if verify_output "$output" "$vcf_sha" 3629446 &&
		verify_output "$output.tbi" "$tbi_sha" 11871; then
		echo "have $output"
	else
		bcftools view --no-version -r chr22 -G -Ou "$source" |
			bcftools annotate --no-version \
				--rename-chrs <(printf 'chr22\t22\n') \
				-Oz -o "$partial"
		finish_vcf "$partial" "$output" "$vcf_sha" 3629446 \
			"$tbi_sha" 11871
	fi
	write_receipt "$output" "$source" \
		"etag=$source_etag;source=Sniffles2_2.0.7;region=chr22;genotypes_removed=true;chrom_alias=chr22:22"
}

stage_dbvar() {
	local dir="$OUTPUT_DIR/dbvar/GRCh38_20260127"
	local output="$dir/GRCh38.variant_region.all.chr22.vcf.gz"
	local partial="$output.partial.$$"
	local source="https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_region.all.vcf.gz"
	local manifest="https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/@md5Sum.md5"
	local source_md5="56b1ce7d343c8982bd84076b63e0cd81"
	local vcf_sha="60e880b451e8878a047b983f3a207ffb4a5426d86ff8ab4409a88430e4062119"
	local tbi_sha="fdb77069b9bb5b6b9a37bce3765eb40a9d492d436864a3671cacc0302d5b9237"

	mkdir -p "$dir"
	if ! curl -fsSL "$manifest" |
		grep -Fq "$source_md5  GRCh38.variant_region.all.vcf.gz"; then
		echo "FATAL: dbVar source MD5 is absent from $manifest" >&2
		exit 1
	fi
	if verify_output "$output" "$vcf_sha" 2031057 &&
		verify_output "$output.tbi" "$tbi_sha" 1771; then
		echo "have $output"
	else
		bcftools view --no-version -r 22 -Oz -o "$partial" "$source"
		finish_vcf "$partial" "$output" "$vcf_sha" 2031057 \
			"$tbi_sha" 1771
	fi
	write_receipt "$output" "$source" \
		"fileDate=20260127;source_md5=$source_md5;region=22"
}

case "$CORPUS" in
	all)
		stage_hprc
		stage_sniffles2
		stage_dbvar
		;;
	hprc-african4-chr22)
		stage_hprc
		;;
	sniffles2-chr22)
		stage_sniffles2
		;;
	dbvar-chr22)
		stage_dbvar
		;;
	*)
		echo "FATAL: unknown corpus: $CORPUS" >&2
		exit 1
		;;
esac

echo "DuckVEP external corpora staged under $OUTPUT_DIR"
