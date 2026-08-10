#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-liftover-registry.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT
CACHE_DIR="$TMP_DIR/cache"
BIN_DIR="$TMP_DIR/bin"
mkdir -p "$BIN_DIR"

cat >"$BIN_DIR/curl" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output|-o) output="$2"; shift 2 ;;
    *) shift ;;
  esac
done
mkdir -p "$(dirname "$output")"
printf '##fileformat=VCFv4.2\n' >"$output"
EOF
cat >"$BIN_DIR/bcftools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  index) touch "${@: -1}.tbi" ;;
  --version) printf 'bcftools fake 1.0\n' ;;
esac
EOF
chmod +x "$BIN_DIR/curl" "$BIN_DIR/bcftools"

registry="$TMP_DIR/registry.tsv"
printf '%s\n' \
  $'id\tworkload\trole\trelease\tlocator\taccess\tcache_relpath\ttransform\tconsumer\tstage_order\tsupplier_identity' \
  $'fixture\tliftover\tbenchmark_vcf\ttest\thttps://example.test/fixture.vcf.gz\tpublic\tdatasets/giab/fixture.vcf.gz\tdirect_download\ttest\t1\t' \
  >"$registry"
cases="$TMP_DIR/cases.tsv"
printf '%s\n' \
  $'case_id\tdataset\tsample\tdescription\tdefault_region\tartifact_id' \
  $'fixture_case\tgiab\tHGTEST\tfixture\t20\tfixture' \
  >"$cases"
for file in chain source.fa destination.fa; do printf 'fixture\n' >"$TMP_DIR/$file"; done

PATH="$BIN_DIR:$PATH" \
  DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  DUCKHTSBENCH_REGISTRY="$registry" \
  LIFTOVER_CASES_TSV="$cases" \
  LIFTOVER_OUT_DIR="$TMP_DIR/out" \
  LIFTOVER_CHAIN_PATH="$TMP_DIR/chain" \
  LIFTOVER_SRC_FASTA="$TMP_DIR/source.fa" \
  LIFTOVER_DST_FASTA="$TMP_DIR/destination.fa" \
  LIFTOVER_STAGE_ONLY=1 \
  bash "$ROOT_DIR/scripts/liftover_conformance_batch.sh" fixture_case >/dev/null

output="$CACHE_DIR/datasets/giab/fixture.vcf.gz"
[[ -s "$output" ]]
[[ -f "$output.tbi" ]]
[[ -f "$output.provenance.tsv" ]]
grep -Fq $'artifact_id\tfixture' "$output.provenance.tsv"
echo "Liftover registry batch staging: OK"
