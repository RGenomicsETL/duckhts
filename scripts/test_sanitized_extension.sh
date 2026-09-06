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
for format in bcf full.vcf.gz; do
  run_sql "" \
    "SET threads=4; SELECT * FROM bcf_index('test/data/bcf_scan_contigs.$format',index_path:='$tmp/reload.$format.index',threads:=1); PREPARE scan AS SELECT CASE WHEN list(POS ORDER BY POS)=[10,20,20,20,40] THEN true ELSE error('index snapshot lost rows') END FROM read_bcf('test/data/bcf_scan_contigs.$format',index_path:='$tmp/reload.$format.index',decompression_threads:=0); COPY (SELECT 'broken index') TO '$tmp/reload.$format.index' (FORMAT CSV,HEADER FALSE); EXECUTE scan;"
done
for header in full partial none; do
  run_sql "" \
    "SET threads=4; CREATE TABLE scanned AS SELECT * FROM read_bcf('test/data/bcf_scan_contigs.$header.vcf.gz',index_path:='test/data/bcf_scan_contigs.$header.vcf.gz.index.tbi',decompression_threads:=0); SELECT CASE WHEN count(*)=5 AND count(DISTINCT CHROM)=2 THEN true ELSE error('lost contig records') END FROM scanned;"
done
run_sql "region list: empty item" \
  "SELECT count(*) FROM read_bam('test/data/range.bam',region:=',,');"
run_sql "region list: invalid item" \
  "SELECT count(*) FROM read_bcf('test/data/region_union.bcf',region:='chr1:1-20,chr1:20-10');"
run_sql "region list: invalid item" \
  "SELECT count(*) FROM read_bcf('test/data/region_union_no_contig.vcf.gz',index_path:='test/data/region_union_no_contig.index.tbi',region:='chr1:1-20,chr1:20-10');"
run_sql "region list: empty item" \
  "SELECT count(*) FROM read_fasta('test/data/ce.fa',region:='CHROMOSOME_I:1-10,');"
run_sql "" \
  "SELECT CASE WHEN count(*)=2 AND min(NAME)='chr,part' AND min(SEQUENCE)='ACGTA' THEN true ELSE error('quoted FASTA header name') END FROM read_fasta('test/data/region_names.fa.gz',region:='{chr,part}:1-5,{chr,part}:1-5');"
run_sql "region list: invalid item" \
  "SELECT count(*) FROM read_tabix('test/data/gff_file.gff.gz',region:='X:20-10');"
run_sql "" \
  "SELECT CASE WHEN count(*)=4 AND count(*) FILTER (WHERE ID='long_ref')=2 AND count(*) FILTER (WHERE ID='long_end')=2 THEN true ELSE error('BCF long-record region union') END FROM read_bcf('test/data/region_union.bcf',region:='chr1:19-19,chr1:11-11,chr1:19-19',tidy_format:=true); SELECT CASE WHEN count(*)=4 THEN true ELSE error('VCF stored duplicates') END FROM read_bcf('test/data/region_union.vcf.gz',region:='chr1:18-18,chr1:18-18',tidy_format:=true) WHERE ID='duplicate';"
run_sql "read_bam: failed to read SAM/BAM/CRAM record" \
  "SELECT QNAME FROM read_bam('test/data/bam_scan_malformed.sam',scan_mode:='sequential',decompression_threads:=0);"
run_sql "" \
  "SET threads=4; SELECT CASE WHEN count(QNAME)=2053 THEN true ELSE error('legacy BAM tail ownership') END FROM read_bam('test/data/bam_scan_all_unplaced.bam',index_path:='test/data/bam_scan_all_unplaced.bam.legacy.bai',decompression_threads:=0);"
run_sql "" \
  "SET threads=4; SELECT CASE WHEN count(QNAME)=5 AND count(*) FILTER (WHERE QNAME='unplaced')=2 THEN true ELSE error('BAM tail ownership') END FROM read_bam('test/data/bam_scan_mixed.bam',decompression_threads:=0); SELECT CASE WHEN count(QNAME)=2053 THEN true ELSE error('CRAM tail ownership') END FROM read_bam('test/data/bam_scan_all_unplaced.cram',decompression_threads:=0);"
for input in bam_offset.sam bam_offset.sam.gz bam_offset.sam.bgz bam_offset_uncompressed.bam bam_offset_gzip.bam bam_scan_mixed.cram; do
  for mode in auto sequential; do
    run_sql "" \
      "SET threads=4; SELECT CASE WHEN count(*)=5 AND count(FILE_OFFSET)=0 AND count(*) FILTER (WHERE QNAME='unplaced')=2 THEN true ELSE error('unsupported FILE_OFFSET transport') END FROM read_bam('test/data/$input',scan_mode:='$mode');"
  done
done
run_sql "" \
  "SET threads=4; CREATE TEMP TABLE actual AS SELECT QNAME,FILE_OFFSET FROM read_bam('test/data/nanopore.bam'); CREATE TEMP TABLE expected AS SELECT QNAME,FILE_OFFSET FROM read_csv('test/data/bam_record_ends.tsv') WHERE file='nanopore.bam'; SELECT CASE WHEN count(*)=0 AND (SELECT count(*) FROM actual)=186 THEN true ELSE error('BAM post-record offset') END FROM ((SELECT * FROM actual EXCEPT ALL SELECT * FROM expected) UNION ALL (SELECT * FROM expected EXCEPT ALL SELECT * FROM actual));"
run_sql "precision_digits must be between 0 and 18" \
  "SELECT * FROM duckhts_mosdepth('$tmp/mosdepth','test/data/range.bam',precision_digits:=-1,overwrite:=true);"
run_sql "" \
  "SELECT count(*) FROM bcftools_score('test/data/score_input.vcf','test/data/score_summary.tsv',use:='GT',columns:='PLINK');"
run_sql "" \
  "SELECT CASE WHEN count(*)=2053 AND count(*) FILTER (WHERE VEP_SYMBOL=['GENE1'] AND VEP_DISTANCE=[12])=2053 THEN true ELSE error('tidy CSQ record lifetime') END FROM read_bcf('test/data/tidy_chunk_boundary.vcf',tidy_format:=true); SELECT POS,SAMPLE_ID,VEP_SYMBOL,INFO_DP,FORMAT_GT,FORMAT_AD,FORMAT_ST FROM read_bcf('test/data/bcf_cache_lifecycle.vcf',tidy_format:=true);"

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

# Exercise the shared first-party cache without DuckDB execution, retaining
# leak checks and exact byte oracles across eviction and concurrent fetches.
${CC:-cc} -std=c11 -D_DEFAULT_SOURCE -Wall -Wextra -Werror "${compile_flags[@]}" \
  -Isrc/include -Ithird_party/htslib test/scripts/reference_cache_test.c \
  -L"$out_dir" -Wl,-rpath,"$root/$out_dir" -lduckhts -pthread "${link_flags[@]}" \
  -o "$tmp/reference_cache_test"
if [ "$sanitizer" = asan ]; then
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
    LD_PRELOAD="$(${CC:-cc} -print-file-name=libasan.so)" "$tmp/reference_cache_test" "$tmp"
else
  "${runtime[@]}" "$tmp/reference_cache_test" "$tmp"
fi
cmake --build "$build_dir" --target duckhts_region_list_test -j2
cmake --build "$build_dir" --target duckhts_bcf_scan_test -j2
if [ "$sanitizer" = asan ]; then
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
    LD_PRELOAD="$(${CC:-cc} -print-file-name=libasan.so)" "$build_dir/duckhts_region_list_test"
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
    LD_PRELOAD="$(${CC:-cc} -print-file-name=libasan.so)" "$build_dir/duckhts_bcf_scan_test" "$tmp"
else
  "${runtime[@]}" "$build_dir/duckhts_region_list_test"
  "${runtime[@]}" "$build_dir/duckhts_bcf_scan_test" "$tmp"
fi
cmake --build "$build_dir" --target duckhts_reader_alloc_probe -j2
${CC:-cc} -std=c11 -UNDEBUG -Wall -Wextra -Werror "${compile_flags[@]}" \
  -Isrc/include -Ithird_party/htslib test/scripts/bam_format_test.c \
  -L"$out_dir" -Wl,-rpath,"$root/$out_dir" -lduckhts "${link_flags[@]}" \
  -o "$tmp/bam_format_test"
"${runtime[@]}" "$tmp/bam_format_test"
python_runtime=("${runtime[@]}")
if [ "$sanitizer" = asan ]; then
  # Python does not link C++; load it at startup so ASan can resolve
  # __cxa_throw before DuckDB is imported and exercises query errors.
  python_runtime=(env "ASAN_OPTIONS=detect_leaks=0:halt_on_error=1"
    "LD_PRELOAD=$(${CC:-cc} -print-file-name=libasan.so):$(${CXX:-c++} -print-file-name=libstdc++.so)")
fi
"${python_runtime[@]}" ./configure/venv/bin/python3 test/scripts/reader_alloc_test.py \
  --extension "$extension" --probe "$build_dir/libduckhts_reader_alloc_probe.so"
# Reuse the complete typed oracle for nested output as well as the native decoder.
"${python_runtime[@]}" ./configure/venv/bin/python3 scripts/run_sqllogictest.py \
  --test-dir test/sql --file-path test/sql/geno.test --external-extension "$extension"
echo "$sanitizer complete-extension gates: OK"
