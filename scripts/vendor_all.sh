#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/vendor_htslib.sh"
"$SCRIPT_DIR/vendor_bcftools.sh"
"$SCRIPT_DIR/vendor_samtools.sh"
"$SCRIPT_DIR/vendor_cgranges.sh"

if [[ "${VENDOR_CONFORMANCE:-0}" == "1" ]]; then
  "$SCRIPT_DIR/vendor_conformance_data.sh"
fi
