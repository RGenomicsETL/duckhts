#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFORMANCE_DIR="$ROOT_DIR/third_party/conformance"
TP_DIR="$ROOT_DIR/third_party"

mkdir -p "$CONFORMANCE_DIR"

cat > "$CONFORMANCE_DIR/README.md" <<'EOF'
# DuckHTS Conformance Data

This directory holds datasets used for conformance testing. Data is staged from
vendored source trees (htslib/bcftools/samtools) to keep tests offline and
reproducible.

Expected contents (examples):
- VCF/BCF: small indexed cohort files
- BAM/CRAM: small alignment files with indexes
- Reference FASTA (if needed for CRAM)

Add data using a reproducible, checksum-verified vendoring script.
EOF

copy_if_exists() {
	local src="$1"
	local dest="$2"
	if [[ -d "$src" ]]; then
		rm -rf "$dest"
		mkdir -p "$(dirname "$dest")"
		cp -a "$src" "$dest"
		echo "Copied $src -> $dest"
	fi
}

copy_if_exists "$TP_DIR/htslib/test" "$CONFORMANCE_DIR/htslib/test"
copy_if_exists "$TP_DIR/bcftools/test" "$CONFORMANCE_DIR/bcftools/test"
copy_if_exists "$TP_DIR/samtools/test" "$CONFORMANCE_DIR/samtools/test"

# Generate a checksum manifest for reproducibility
if command -v sha256sum >/dev/null 2>&1; then
	(cd "$CONFORMANCE_DIR" && find . -type f -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256)
fi

echo "Prepared $CONFORMANCE_DIR from vendored test data."
