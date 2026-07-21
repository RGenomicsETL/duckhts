DuckHTS VariantKey Benchmark and Conformance
================

<!-- benchmark_variantkey_conformance.md is generated from benchmark_variantkey_conformance.Rmd. -->

# Benchmark and conformance

This report separates four different workloads so the full conformance
harness is not mistaken for raw VariantKey encode throughput:

1.  conservative contig-key normalization
2.  pure in-memory VariantKey encoding
3.  real `read_bcf(...)` scan + VariantKey encoding
4.  full DuckHTS-vs-bcftools `%VKX` conformance

DuckHTS uses the vendored official VariantKey C API and follows the
bcftools `%VKX` / `+add-variantkey` convention of accepting 1-based VCF
`POS` at the SQL boundary while encoding the upstream 0-based field
internally.

VariantKey / RegionKey are based on Nicola Asuni’s work:
<https://doi.org/10.1101/473744>.

## Run

``` sh
make bench-variantkey
```

Useful overrides:

- `BCFTOOLS_BIN`: optional bcftools executable override
- `VARIANTKEY_REAL_CASES`: optional comma-separated local VCF/BCF paths
- `VARIANTKEY_OUT_DIR`: optional directory for emitted TSVs and
  summaries
- `DUCKHTS_EXTENSION`: optional extension path override
- `VARIANTKEY_BENCH_ROWS`: synthetic benchmark rows, default `4000000`
- `VARIANTKEY_BENCH_RUNS`: timed repeats for the benchmark sections,
  default `3`
- `VARIANTKEY_BENCH_CASE`: optional real VCF/BCF path for the
  scan+encode section
- `VARIANTKEY_DUCK_THREADS_LIST`: DuckDB thread grid for the benchmark
  sections, default `1,2,4`
- `VARIANTKEY_CONFORMANCE_DUCK_THREADS`: DuckDB threads for the full
  conformance harness, default `max(VARIANTKEY_DUCK_THREADS_LIST)`
- `VARIANTKEY_DUCK_DECOMPRESSION_THREADS`: `read_bcf(...)` worker
  threads, default `0`

If `VARIANTKEY_REAL_CASES` is unset, this report looks for the locally
downloaded whole-WGS GIAB files already used in earlier DuckHTS
conformance work under `/root/giab_norm/`.

## Benchmark settings

``` r
benchmark_settings
#>                           setting                  value
#> 1                  synthetic_rows                4000000
#> 2                  benchmark_runs                      3
#> 3             duckdb_threads_grid                  1,2,4
#> 4 full_conformance_duckdb_threads                      4
#> 5  read_bcf_decompression_threads                      0
#> 6                scan_encode_case giab_hg001_grch37_full
```

## 1. Conservative contig join key

This section forces materialization of `duckhts_contig_key(...)` over a
mixed synthetic stream and compares it with the equivalent built-in SQL
expression. Both paths remove one non-empty leading `chr`, normalize
`M`/`MT`, and uppercase `X`/`Y`; neither maps accessions or numeric sex
chromosomes. The result maximum is retained so DuckDB cannot reduce the
query to a row count that skips the key.

``` r
contig_source_sql <- sprintf(
  paste(
    "SELECT CASE i %% 8",
    "  WHEN 0 THEN 'chr1' WHEN 1 THEN '1'",
    "  WHEN 2 THEN 'chrM' WHEN 3 THEN 'MT'",
    "  WHEN 4 THEN 'chrX' WHEN 5 THEN 'Y'",
    "  WHEN 6 THEN 'chrGL000220.1' ELSE 'NC_000001.11'",
    "END AS contig FROM range(%d) t(i)"
  ),
  bench_rows
)

contig_key_results <- rbind(
  measure_thread_grid(
    "duckhts_contig_key",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          "SELECT count(*) AS n, max(duckhts_contig_key(contig)) AS max_key FROM (%s)",
          contig_source_sql
        )
      )
    }
  ),
  measure_thread_grid(
    "equivalent_builtin_sql",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          paste(
            "WITH source AS (%s), stripped AS (",
            "  SELECT CASE WHEN length(contig) > 3 AND lower(substr(contig, 1, 3)) = 'chr'",
            "              THEN substr(contig, 4) ELSE contig END AS contig",
            "  FROM source",
            ")",
            "SELECT count(*) AS n,",
            "       max(CASE WHEN lower(contig) IN ('m', 'mt') THEN 'MT'",
            "                WHEN lower(contig) = 'x' THEN 'X'",
            "                WHEN lower(contig) = 'y' THEN 'Y'",
            "                ELSE contig END) AS max_key",
            "FROM stripped"
          ),
          contig_source_sql
        )
      )
    }
  )
)

contig_key_results
#>   threads              operation  rows median_seconds rows_per_second
#> 1       1     duckhts_contig_key 4e+06          0.104        38461538
#> 2       2     duckhts_contig_key 4e+06          0.099        40404040
#> 3       4     duckhts_contig_key 4e+06          0.101        39603960
#> 4       1 equivalent_builtin_sql 4e+06          0.268        14925373
#> 5       2 equivalent_builtin_sql 4e+06          0.267        14981273
#> 6       4 equivalent_builtin_sql 4e+06          0.267        14981273
```

## 2. Pure in-memory encode kernel

This section isolates the VariantKey UDF itself from HTS decoding, file
I/O, and conformance comparison work.

``` r
pure_kernel_results <- rbind(
  measure_thread_grid(
    "numeric_variantkey",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          "SELECT count(*) AS n, max(variantkey('1', i, 'A', 'C')) AS max_vk FROM range(1, %d) t(i)",
          bench_rows + 1L
        )
      )
    }
  ),
  measure_thread_grid(
    "variantkey_hex",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          paste(
            "SELECT count(*) AS n,",
            "max(variantkey_hex(variantkey('1', i, 'A', 'C'))) AS max_vkx",
            "FROM range(1, %d) t(i)"
          ),
          bench_rows + 1L
        )
      )
    }
  )
)

pure_kernel_results
#>   threads          operation  rows median_seconds rows_per_second
#> 1       1 numeric_variantkey 4e+06          0.086        46511628
#> 2       2 numeric_variantkey 4e+06          0.087        45977011
#> 3       4 numeric_variantkey 4e+06          0.085        47058824
#> 4       1     variantkey_hex 4e+06          0.310        12903226
#> 5       2     variantkey_hex 4e+06          0.312        12820513
#> 6       4     variantkey_hex 4e+06          0.311        12861736
```

## 3. Real BCF scan + encode

This section measures the cost of scanning a real indexed GIAB BCF/VCF
through `read_bcf(...)` and computing VariantKey values, without TSV
materialization or bcftools-side conformance comparison.

``` r
scan_case_sql <- sql_string(bench_case_path)
scan_results <- rbind(
  measure_thread_grid(
    "read_bcf_count",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          "SELECT count(*) AS n FROM read_bcf(%s, decompression_threads := %d)",
          scan_case_sql,
          duck_decompression_threads
        )
      )
    }
  ),
  measure_thread_grid(
    "read_bcf_plus_numeric_variantkey",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          paste(
            "SELECT count(*) AS n,",
            "max(variantkey(CHROM, POS, REF, coalesce(ALT[1], '.'))) AS max_vk",
            "FROM read_bcf(%s, decompression_threads := %d)"
          ),
          scan_case_sql,
          duck_decompression_threads
        )
      )
    }
  ),
  measure_thread_grid(
    "read_bcf_plus_variantkey_hex",
    function() {
      DBI::dbGetQuery(
        con,
        sprintf(
          paste(
            "SELECT count(*) AS n,",
            "max(variantkey_hex(variantkey(CHROM, POS, REF, coalesce(ALT[1], '.')))) AS max_vkx",
            "FROM read_bcf(%s, decompression_threads := %d)"
          ),
          scan_case_sql,
          duck_decompression_threads
        )
      )
    }
  )
)
scan_results <- cbind(case = bench_case_label, scan_results, stringsAsFactors = FALSE)
scan_results
#>                     case threads                        operation    rows
#> 1 giab_hg001_grch37_full       1                   read_bcf_count 3891440
#> 2 giab_hg001_grch37_full       2                   read_bcf_count 3891440
#> 3 giab_hg001_grch37_full       4                   read_bcf_count 3891440
#> 4 giab_hg001_grch37_full       1 read_bcf_plus_numeric_variantkey 3891440
#> 5 giab_hg001_grch37_full       2 read_bcf_plus_numeric_variantkey 3891440
#> 6 giab_hg001_grch37_full       4 read_bcf_plus_numeric_variantkey 3891440
#> 7 giab_hg001_grch37_full       1     read_bcf_plus_variantkey_hex 3891440
#> 8 giab_hg001_grch37_full       2     read_bcf_plus_variantkey_hex 3891440
#> 9 giab_hg001_grch37_full       4     read_bcf_plus_variantkey_hex 3891440
#>   median_seconds rows_per_second
#> 1          4.671        833106.4
#> 2          2.392       1626856.2
#> 3          1.274       3054505.5
#> 4          5.128        758861.2
#> 5          2.641       1473472.2
#> 6          1.393       2793567.8
#> 7          5.359        726150.4
#> 8          2.749       1415583.8
#> 9          1.449       2685603.9
```

## 4. Full DuckHTS vs bcftools conformance harness

This final section is the end-to-end parity workflow. It includes both
sides of the comparison and is therefore much slower than the pure
encode kernel.

``` r
results <- do.call(
  rbind,
  lapply(seq_len(nrow(cases)), function(i) run_case(cases$label[i], cases$path[i]))
)
results
#>                    label
#> 1 giab_hg001_grch37_full
#> 2 giab_hg002_grch37_full
#> 3 giab_hg006_grch37_full
#>                                                        path duckdb_threads
#> 1 /root/giab_norm/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz              4
#> 2 /root/giab_norm/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz              4
#> 3 /root/giab_norm/HG006_GRCh37_1_22_v4.2.1_benchmark.vcf.gz              4
#>   bcftools_rows duckhts_rows mismatch_groups elapsed_seconds rows_per_second
#> 1       3891440      3891440               0           4.578        850030.6
#> 2       4033796      4033796               0           5.104        790320.5
#> 3       3878664      3878664               0           4.427        876138.2
```

All listed cases should have `mismatch_groups == 0`.

## Notes

- For paper-style raw VariantKey throughput comparisons, use the **pure
  in-memory encode kernel** section rather than the full conformance
  harness.
- The **real BCF scan + encode** section includes HTS decoding and
  DuckDB table function overhead, so it is the right comparison for
  DuckHTS-as-reader performance rather than bare header-only VariantKey
  microbenchmarks.
- The **full conformance harness** includes bcftools `%VKX` export,
  DuckHTS export, TSV materialization, and rowwise comparison, so it
  should be read as a parity-validation workflow, not as a bare encoding
  benchmark.
- This framework compares against bcftools `%VKX`, which is the same
  VariantKey hexadecimal output exposed by `+add-variantkey`.
- DuckHTS currently compares `ALT[1]` from `read_bcf(...)`, mirroring
  the biallelic `%VKX` field that bcftools derives from the first ALT
  allele on the current record.
- Large, ambiguous, and symbolic alleles are expected to fall back to
  the official hashed nonreversible VariantKey mode. That still
  preserves bcftools `%VKX` parity, but it does not encode `END`,
  `SVLEN`, mate breakend coordinates, or other SV metadata.
