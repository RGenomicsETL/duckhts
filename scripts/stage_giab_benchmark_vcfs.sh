#!/usr/bin/env bash
# Compatibility entry point. r/duckhtsbench owns the GIAB source registry,
# cache destination, provenance record, and staging implementation.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_giab.R" --root "$ROOT_DIR"
