#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-staging-cache-test.XXXXXX")"
FAKE_BIN="$TMP_DIR/bin"
CACHE_DIR="$TMP_DIR/cache"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$FAKE_BIN"

cat >"$FAKE_BIN/curl" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
out=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output) out="$2"; shift 2 ;;
    *) url="$1"; shift ;;
  esac
done
mkdir -p "$(dirname "$out")"
if [[ "$url" == *"human_g1k_v37.fasta.gz" ]]; then
  printf '>1\nACGT\n' | gzip -c >"$out"
else
  printf 'staged from %s\n' "$url" >"$out"
fi
EOF

cat >"$FAKE_BIN/bcftools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
  --version) printf 'bcftools mock 1.0\n' ;;
  view)
    out=""
    while [[ $# -gt 0 ]]; do
      case "$1" in -o) out="$2"; shift 2 ;; *) shift ;; esac
    done
    mkdir -p "$(dirname "$out")"
    printf '##fileformat=VCFv4.2\n' >"$out"
    ;;
  index)
    for arg in "$@"; do case "$arg" in *.vcf.gz) target="$arg" ;; esac; done
    touch "${target}.tbi"
    ;;
  query) printf 'HG00096\n' ;;
  *) exit 0 ;;
esac
EOF

cat >"$FAKE_BIN/tabix" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--version" ]]; then printf 'tabix mock 1.0\n'; exit 0; fi
for arg in "$@"; do case "$arg" in *.vcf.gz) target="$arg" ;; esac; done
touch "${target}.tbi"
EOF

cat >"$FAKE_BIN/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
  --version) printf 'samtools mock 1.0\n' ;;
  faidx) touch "${2}.fai" ;;
  quickcheck) : ;;
  view)
    out=""
    while [[ $# -gt 0 ]]; do
      case "$1" in -o) out="$2"; shift 2 ;; *) shift ;; esac
    done
    mkdir -p "$(dirname "$out")"
    printf 'BAM\n' >"$out"
    ;;
  index) touch "${@: -1}.bai" ;;
  *) exit 0 ;;
esac
EOF
chmod +x "$FAKE_BIN"/*

PATH="$FAKE_BIN:$PATH" DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  bash "$ROOT_DIR/scripts/stage_norm_1000g_dragen_gvcf.sh" >/dev/null
[[ -f "$CACHE_DIR/benchmarks/norm/1000g-dragen/HG00096/HG00096.hard-filtered.chr22_20000000_30000000.g.vcf.gz.provenance.tsv" ]]

PATH="$FAKE_BIN:$PATH" DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  bash "$ROOT_DIR/scripts/stage_liftover_references.sh" >/dev/null
[[ -f "$CACHE_DIR/references/liftover/grch37-to-grch38/provenance.tsv" ]]

PATH="$FAKE_BIN:$PATH" DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  bash "$ROOT_DIR/scripts/stage_riker_wgs_bam.sh" >/dev/null
[[ -f "$CACHE_DIR/benchmarks/riker-wgs/downloads/HG00188.final.cram" ]]
[[ -f "$CACHE_DIR/benchmarks/riker-wgs/provenance.tsv" ]]

echo "DuckHTS staging cache defaults: OK"
