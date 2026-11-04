#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
# Built by CMake with the extension's htslib link libraries (Makefile test target).
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cp cmake_build/release/duckhts_fasta_nuc_alloc_test "$work/test"
"$work/test"
