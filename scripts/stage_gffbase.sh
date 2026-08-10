#!/usr/bin/env bash
# Install the benchmark's pinned GFFBase package under the DuckHTS cache.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

PYTHON_BIN=${PYTHON_BIN:-python3}
GFFBASE_VERSION=${GFFBASE_VERSION:-0.1.0}
GFFBASE_SITE=${GFFBASE_SITE:-$(duckhts_cache_subdir benchmarks/gffbase/python-site)}
PIP_CACHE_DIR=${PIP_CACHE_DIR:-$(duckhts_cache_subdir pip)}
export PIP_CACHE_DIR
PROVENANCE="${GFFBASE_SITE}/provenance.tsv"

mkdir -p "$GFFBASE_SITE" "$PIP_CACHE_DIR"
if [[ ! -d "$GFFBASE_SITE/gffbase" ]]; then
    "$PYTHON_BIN" -m pip install --target "$GFFBASE_SITE" "gffbase==$GFFBASE_VERSION"
else
    echo "Using cached GFFBase $GFFBASE_VERSION: $GFFBASE_SITE"
fi
if [[ ! -f "$PROVENANCE" ]]; then
    duckhts_write_provenance "$PROVENANCE" \
        "workload=gffbase_conformance_benchmark" \
        "package=gffbase" \
        "package_version=$GFFBASE_VERSION" \
        "source=PyPI" \
        "site_directory=$GFFBASE_SITE" \
        "pip_cache=$PIP_CACHE_DIR" \
        "staging_command=$PYTHON_BIN -m pip install --target $GFFBASE_SITE gffbase==$GFFBASE_VERSION" \
        "python=$($PYTHON_BIN --version)"
fi
printf 'GFFBASE_SITE=%s\n' "$GFFBASE_SITE"
printf 'GFFBASE_PROVENANCE=%s\n' "$PROVENANCE"
