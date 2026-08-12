#!/usr/bin/env bash
# Compatibility entry point. r/duckhtsbench owns the DuckBedQC source pin,
# checkout validation, cache destination, and provenance receipt.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DUCKHTS_ROOT="$ROOT_DIR" exec Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_duckbedqc.R" "$@"
