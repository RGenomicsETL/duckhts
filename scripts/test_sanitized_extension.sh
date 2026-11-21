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

# The C API version comes from the Makefile's TARGET_DUCKDB_VERSION, for
# example v1.5.6; `make test-sanitized-extension` passes it.
target=${TARGET_DUCKDB_VERSION:?set TARGET_DUCKDB_VERSION, or run make test-sanitized-extension}
if [[ ! $target =~ ^v([0-9]+)\.([0-9]+)\.([0-9]+)$ ]]; then
  echo "TARGET_DUCKDB_VERSION must look like v1.5.6, not '$target'" >&2
  exit 2
fi
target_major=${BASH_REMATCH[1]}
target_minor=${BASH_REMATCH[2]}
target_patch=${BASH_REMATCH[3]}

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
build_dir="cmake_build/sanitizer-$sanitizer"
out_dir="build/sanitizer-$sanitizer"
rm -rf "$build_dir" "$out_dir"
cmake -S . -B "$build_dir" \
  -DEXTENSION_NAME=duckhts \
  -DTARGET_DUCKDB_VERSION_MAJOR="$target_major" \
  -DTARGET_DUCKDB_VERSION_MINOR="$target_minor" \
  -DTARGET_DUCKDB_VERSION_PATCH="$target_patch" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="${compile_flags[*]}" \
  -DCMAKE_SHARED_LINKER_FLAGS="${link_flags[*]}"
cmake --build "$build_dir" --config Debug -j2
mkdir -p "$out_dir"
cp "$build_dir/libduckhts.so" "$out_dir/libduckhts.so"
./configure/venv/bin/python3 extension-ci-tools/scripts/append_extension_metadata.py \
  -l "$out_dir/libduckhts.so" \
  -o "$out_dir/duckhts.duckdb_extension" \
  -n duckhts -dv "$target" \
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
  "CREATE TABLE cgr_san(chrom VARCHAR,start BIGINT,stop BIGINT); INSERT INTO cgr_san VALUES ('chr1',0,1); SELECT * FROM duckhts_cgranges_from_table('cgr_san_idx','cgr_san','chrom','start','stop'); SELECT * FROM duckhts_cgranges_overlaps('cgr_san_idx','chr1',0,1,mode:='bad');"
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
  "CREATE TEMP TABLE roh_map AS SELECT 'chr1' AS chrom, position AS pos, \"Genetic_Map(cM)\" AS cm FROM read_csv('test/data/roh_map_chr1.txt', delim = ' ') UNION ALL SELECT 'chr2', position, \"Genetic_Map(cM)\" FROM read_csv('test/data/roh_map_chr2.txt', delim = ' '); SELECT CASE WHEN count(*) = 4 AND sum(n_markers) = 597 THEN true ELSE error('roh PL segments') END FROM duckhts_roh('test/data/roh_fixture.vcf.gz', 'AF'); SELECT CASE WHEN count(*) = 4 AND sum(n_markers) = 596 THEN true ELSE error('roh GT segments') END FROM duckhts_roh('test/data/roh_fixture.vcf.gz', 'AF', gt_error := 30); SELECT CASE WHEN count(*) = 4 AND sum(n_markers) = 569 THEN true ELSE error('roh map segments') END FROM duckhts_roh_map('test/data/roh_fixture.vcf.gz', 'AF', 'roh_map', rec_rate := 0.5); SELECT CASE WHEN count(*) = 16 THEN true ELSE error('roh transition segments') END FROM duckhts_roh('test/data/roh_fixture.vcf.gz', 'AF', hw_to_az := 1e-5, az_to_hw := 1e-6); SELECT CASE WHEN len(duckhts_roh_segments([]::BIGINT[], []::DOUBLE[], []::INTEGER[][], NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9)) = 0 AND len(duckhts_roh_segments([5000], [0.05], [[60, 30, 0]], NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9)) = 1 THEN true ELSE error('roh degenerate lists') END;"
run_sql "positions must be sorted ascending" \
  "SELECT duckhts_roh_segments([200, 100], [0.3, 0.3], [[0, 30, 60], [0, 30, 60]], NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9);"
run_sql "positions and af lists differ in length" \
  "SELECT duckhts_roh_segments([100, 200], [0.3], [[0, 30, 60], [0, 30, 60]], NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9);"
run_sql "genetic map positions must be strictly ascending" \
  "SELECT duckhts_roh_segments([100], [0.3], [[0, 30, 60]], [200, 100]::BIGINT[], [0.1, 0.2]::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9);"
run_sql "" \
  "SELECT CASE WHEN count(*)=2053 AND count(*) FILTER (WHERE VEP_SYMBOL=['GENE1'] AND VEP_DISTANCE=[12])=2053 THEN true ELSE error('tidy CSQ record lifetime') END FROM read_bcf('test/data/tidy_chunk_boundary.vcf',tidy_format:=true); SELECT POS,SAMPLE_ID,VEP_SYMBOL,INFO_DP,FORMAT_GT,FORMAT_AD,FORMAT_ST FROM read_bcf('test/data/bcf_cache_lifecycle.vcf',tidy_format:=true);"

# Each relation through the DuckDB file system of the connection must equal the
# one htslib's own file access gives in the same session; a failed open must
# carry DuckDB's reason. This runs the routed opens, index lookups, reference
# loads and per-query liftover contexts under the sanitizer.
while IFS= read -r relation; do
  run_sql "" \
    "SET threads=4; CREATE TEMP TABLE expected AS SELECT count(*) AS n, sum(hash(t)) AS row_hashes FROM ($relation) t; SET duckhts_file_system='duckdb'; SELECT CASE WHEN (SELECT n FROM expected) > 0 AND (SELECT (count(*), sum(hash(t))) FROM ($relation) t) = (SELECT (n, row_hashes) FROM expected) THEN true ELSE error('DuckDB file system route differs from htslib') END;"
done <<'RELATIONS'
SELECT * FROM read_bam('test/data/bam_scan_mixed.bam', region := 'chr1')
SELECT * FROM read_bam('test/data/range.cram', reference := 'test/data/ce.fa')
SELECT * FROM bam_bin_counts('test/data/fixture_mixed.cram', 100, reference := 'test/data/fixture_ref.fa')
SELECT * FROM duckhts_bam_mismatch_counts('test/data/bam_mismatch.cram', 'test/data/bam_mismatch.fa', mask := 'test/data/bam_mismatch.mask.vcf.gz')
SELECT * FROM read_bcf('test/data/formatcols.vcf.gz', region := '1')
SELECT * FROM read_geno('test/data/geno_sites.bcf', region := 'chr1')
SELECT * FROM read_gff('test/data/gff_file.gff.gz', region := 'X:2934816-2935190')
SELECT * FROM read_fasta('test/data/ce.fa', region := 'CHROMOSOME_I:1-200')
SELECT * FROM fasta_nuc('test/data/ce.fa', bed_path := 'test/data/targets.bed')
SELECT * FROM read_hts_index('test/data/range.bam')
SELECT * FROM read_genbank('test/data/phix174.gb')
SELECT * FROM read_bigwig('third_party/libBigWig/test/test.bw')
SELECT * FROM bcftools_score('test/data/score_input.vcf', 'test/data/score_summary.tsv', use := 'GT', columns := 'PLINK')
SELECT bcftools_liftover('chrF', position, 'C', 'T', 'test/data/liftover.chain', 'test/data/liftover_dst.fa', 'test/data/liftover_src.fa', 1, 250, false, NULL::BIGINT, false) FROM range(1, 4000) positions(position)
RELATIONS
while IFS= read -r relation; do
  run_sql "No such file or directory" "SET duckhts_file_system='duckdb'; $relation;"
done <<'RELATIONS'
SELECT * FROM read_bam('test/data/no_such_input')
SELECT * FROM read_bcf('test/data/no_such_input')
SELECT * FROM read_fasta('test/data/no_such_input')
SELECT * FROM read_bigwig('test/data/no_such_input')
SELECT bcftools_liftover('chrF', 2, 'C', 'T', 'test/data/no_such_input', 'test/data/liftover_dst.fa', 'test/data/liftover_src.fa', 1, 250, false, NULL::BIGINT, false)
RELATIONS

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
if [ "$(uname -s)" = Linux ]; then
  cmake --build "$build_dir" --target duckhts_bcf_info_oom_test -j2
  if [ "$sanitizer" = asan ]; then
    ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
      LD_PRELOAD="$(${CC:-cc} -print-file-name=libasan.so)" "$build_dir/duckhts_bcf_info_oom_test"
    cmake --build "$build_dir" --target duckhts_hts_region_ownership_test -j2
    ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
      LD_PRELOAD="$(${CC:-cc} -print-file-name=libasan.so)" \
      "$build_dir/duckhts_hts_region_ownership_test" \
        test/data/range.bam test/data/range.bam.bai \
        test/data/range.cram test/data/range.cram.crai
  else
    "${runtime[@]}" "$build_dir/duckhts_bcf_info_oom_test"
  fi
fi
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
"${python_runtime[@]}" ./configure/venv/bin/python3 scripts/run_sqllogictest.py \
  --test-dir test/sql --file-path test/sql/bcf_malformed_error.test --external-extension "$extension"
# CIGAR checks assert numeric results as well as memory safety: UBSan does not
# diagnose unsigned wraparound, and invalid rows must not leak list entries.
for cigar_test in cigar_utils cigar_aligned_blocks bam_cigar_binary; do
  "${python_runtime[@]}" ./configure/venv/bin/python3 scripts/run_sqllogictest.py \
    --test-dir test/sql --file-path "test/sql/$cigar_test.test" --external-extension "$extension"
done
echo "$sanitizer complete-extension gates: OK"
