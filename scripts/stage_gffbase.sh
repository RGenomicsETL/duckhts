#!/usr/bin/env bash
# Compatibility entry point. r/duckhtsbench owns the GFFBase registry,
# cache destination, receipt, and staging implementation.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DUCKHTS_ROOT="$ROOT_DIR" exec Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_gffbase.R" "$@"
