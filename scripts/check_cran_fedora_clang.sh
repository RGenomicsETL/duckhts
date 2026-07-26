#!/usr/bin/env bash
# Reproduce the CRAN-like Fedora 44 / Clang 22 package check.
# See scripts/check_cran_fedora.sh for the pinned image and scope.
set -euo pipefail

export DUCKHTS_FEDORA_PROFILE=clang
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/check_cran_fedora.sh"
