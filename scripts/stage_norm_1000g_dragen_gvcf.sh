#!/usr/bin/env bash
# Compatibility entry point. r/duckhtsbench owns the 1000G DRAGEN source,
# cache output, receipt, and slice/index derivation.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DUCKHTS_ROOT="$ROOT_DIR" exec Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_norm.R" "$@"
