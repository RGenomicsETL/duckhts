#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-conformance-plugin.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT
SOURCE_DIR="$TMP_DIR/source"
PLUGIN_DIR="$TMP_DIR/cache/plugins"
mkdir -p "$SOURCE_DIR/plugins"

for plugin in score.so munge.so liftover.so; do
  printf 'original %s\n' "$plugin" >"$SOURCE_DIR/plugins/$plugin"
done
(
  cd "$SOURCE_DIR/plugins"
  zip -q "$SOURCE_DIR/score_1.22-20250819.zip" score.so munge.so liftover.so
)

# shellcheck source=scripts/conformance/common.sh
source "$ROOT_DIR/scripts/conformance/common.sh"
DUCKHTS_SCORE_PLUGIN_DIR="$PLUGIN_DIR"
resolved="$(conformance_resolve_plugin_dir "$SOURCE_DIR")"
[[ "$resolved" == "$PLUGIN_DIR" ]]
[[ -f "$PLUGIN_DIR/.duckhts-score-plugin.identity.tsv" ]]
grep -Fxq 'original score.so' "$PLUGIN_DIR/score.so"

# A persistent cache must not accept arbitrary replacement oracle bytes.
printf 'modified score.so\n' >"$PLUGIN_DIR/score.so"
conformance_resolve_plugin_dir "$SOURCE_DIR" >/dev/null
grep -Fxq 'original score.so' "$PLUGIN_DIR/score.so"

echo "Conformance plugin cache identity: OK"
