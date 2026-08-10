#!/usr/bin/env bash
# Compatibility entry point. r/duckhtsbench owns the liftover registry,
# supplier-identity checks, cache destination, receipt, and stage implementation.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DUCKHTS_ROOT="$ROOT_DIR" exec Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_liftover.R" "$@"
