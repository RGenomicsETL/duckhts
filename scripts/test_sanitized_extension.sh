#!/usr/bin/env bash
set -euo pipefail

sanitizer=${1:-asan}
case "$sanitizer" in
  asan)
    compile_flags=(-O1 -g -fsanitize=address -fno-omit-frame-pointer)
    link_flags=(-fsanitize=address)
    runtime=(env "ASAN_OPTIONS=detect_leaks=0:halt_on_error=1" "LD_PRELOAD=$(${CC:-cc} -print-file-name=libasan.so)")
    ;;
  ubsan)
    compile_flags=(-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined)
    link_flags=(-fsanitize=undefined)
    runtime=(env "UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
    ;;
  *) echo "usage: $0 {asan|ubsan}" >&2; exit 2 ;;
esac

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
build_dir="cmake_build/sanitizer-$sanitizer"
out_dir="build/sanitizer-$sanitizer"
rm -rf "$build_dir" "$out_dir"
cmake -S . -B "$build_dir" \
  -DEXTENSION_NAME=duckhts \
  -DTARGET_DUCKDB_VERSION_MAJOR=1 \
  -DTARGET_DUCKDB_VERSION_MINOR=2 \
  -DTARGET_DUCKDB_VERSION_PATCH=0 \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="${compile_flags[*]}" \
  -DCMAKE_SHARED_LINKER_FLAGS="${link_flags[*]}"
cmake --build "$build_dir" --config Debug -j2
mkdir -p "$out_dir"
cp "$build_dir/libduckhts.so" "$out_dir/libduckhts.so"
./configure/venv/bin/python3 extension-ci-tools/scripts/append_extension_metadata.py \
  -l "$out_dir/libduckhts.so" \
  -o "$out_dir/duckhts.duckdb_extension" \
  -n duckhts -dv v1.2.0 \
  -evf configure/extension_version.txt -pf configure/platform.txt

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
duckdb_bin=${DUCKDB_BIN:-duckdb}
extension="$root/$out_dir/duckhts.duckdb_extension"

run_sql() {
  local expected=$1 sql=$2 rc
  printf "LOAD '%s';\n%s\nSELECT 4242;\n" "$extension" "$sql" >"$tmp/input.sql"
  set +e
  "${runtime[@]}" "$duckdb_bin" -unsigned -csv <"$tmp/input.sql" >"$tmp/stdout" 2>"$tmp/stderr"
  rc=$?
  set -e
  if grep -Eq 'AddressSanitizer|runtime error:|UndefinedBehaviorSanitizer' "$tmp/stderr"; then
    cat "$tmp/stderr" >&2
    return 1
  fi
  grep -qx '4242' "$tmp/stdout" || { cat "$tmp/stdout" >&2; cat "$tmp/stderr" >&2; return 1; }
  if [ -n "$expected" ]; then
    if [ "$rc" -eq 0 ] || ! grep -Fq "$expected" "$tmp/stderr"; then
      cat "$tmp/stderr" >&2
      return 1
    fi
  elif [ "$rc" -ne 0 ]; then
    cat "$tmp/stderr" >&2
    return 1
  fi
}

run_sql "clip-pad alignment exceeds the 4194304-cell limit" \
  "SELECT bcftools_liftover('chrS',2,substr(repeat('ACGT',551),2,2200),'C','test/data/liftover_nw_limit.chain','test/data/liftover_nw_limit_dst.fa','test/data/liftover_nw_limit_src.fa',1,250,false,NULL::BIGINT,true);"
run_sql "mode must be 'overlap' or 'contain'" \
  "CREATE TABLE cgr_san(chrom VARCHAR,start BIGINT,stop BIGINT); INSERT INTO cgr_san VALUES ('chr1',0,1); SELECT duckhts_cgranges_from_query('cgr_san_idx','SELECT * FROM cgr_san','chrom','start','stop'); SELECT * FROM duckhts_cgranges_overlaps('cgr_san_idx','chr1',0,1,mode:='bad');"
run_sql "read_bcf: failed to read or parse BCF/VCF record" \
  "SELECT count(*) FROM read_bcf('test/data/malformed_bad_pos.vcf',tidy_format:=true);"
run_sql "precision_digits must be between 0 and 18" \
  "SELECT * FROM duckhts_mosdepth('$tmp/mosdepth','test/data/range.bam',precision_digits:=-1,overwrite:=true);"
run_sql "" \
  "SELECT count(*) FROM bcftools_score('test/data/score_input.vcf','test/data/score_summary.tsv',use:='GT',columns:='PLINK');"

${CC:-cc} -std=c11 -Wall -Wextra -Werror "${compile_flags[@]}" \
  -Isrc/include -Ithird_party/htslib test/scripts/bcftools_filter_recovery_probe.c \
  -L"$out_dir" -Wl,-rpath,"$root/$out_dir" -lduckhts "${link_flags[@]}" \
  -o "$tmp/filter_recovery"
if [ "$sanitizer" = asan ]; then
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
    LD_PRELOAD="$(${CC:-cc} -print-file-name=libasan.so)" "$tmp/filter_recovery"
else
  "${runtime[@]}" "$tmp/filter_recovery"
fi
echo "$sanitizer complete-extension gates: OK"
