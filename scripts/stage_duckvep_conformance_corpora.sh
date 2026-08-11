#!/usr/bin/env bash
# Compatibility entry point. r/duckhtsbench owns DuckVEP corpus source
# identities, cache paths, deterministic derivations, indexes, and receipts.
#
# Usage:
#   scripts/stage_duckvep_conformance_corpora.sh [OUTPUT_DIR] [CORPUS]
#
# OUTPUT_DIR defaults to the registry paths below $DUCKHTS_CACHE_DIR. CORPUS is
# one of: all (default), hprc-african4-chr22, sniffles2-chr22, dbvar-chr22.
set -euo pipefail

if [[ $# -gt 2 ]]; then
  echo "usage: $0 [OUTPUT_DIR] [CORPUS]" >&2
  exit 2
fi
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="${1:-}"
CORPUS="${2:-all}"
ARGS=(--corpus "$CORPUS")
if [[ -n "$OUTPUT_DIR" ]]; then
  ARGS+=(--output-dir "$OUTPUT_DIR")
fi
exec Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_duckvep_conformance_corpora.R" "${ARGS[@]}"
