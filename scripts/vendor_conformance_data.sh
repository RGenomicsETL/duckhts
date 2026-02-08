#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFORMANCE_DIR="$ROOT_DIR/third_party/conformance"

mkdir -p "$CONFORMANCE_DIR"

cat > "$CONFORMANCE_DIR/README.md" <<'EOF'
# DuckHTS Conformance Data

This directory holds small, indexed datasets used for conformance testing.

Expected contents (examples):
- VCF/BCF: small indexed cohort files
- BAM/CRAM: small alignment files with indexes
- Reference FASTA (if needed for CRAM)

Add data using a reproducible, checksum-verified vendoring script.
EOF

echo "Prepared $CONFORMANCE_DIR (data vendoring not yet implemented)."
