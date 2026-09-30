DuckHTS VariantKey Join and Overlap Benchmark
================

<!-- benchmark_variantkey_join_overlap.md is generated from benchmark_variantkey_join_overlap.Rmd. -->

# Benchmark

This report measures synthetic join and interval workloads implemented
by DuckHTS. It compares exact reversible and mixed `VariantKey` joins,
exact `RegionKey` span joins, and interval overlap using
RegionKey-bounded candidate ranges, DuckDB range joins, and
`duckhts_cgranges_*`. The separate
`benchmark_variantkey_conformance.Rmd` report measures key encoding and
`%VKX` parity; this report measures join and overlap throughput.

## Run

``` sh
make bench-variantkey-join
```

Useful overrides:

- `VARIANTKEY_JOIN_ROWS`: synthetic exact-join rows, default `1000000`
- `VARIANTKEY_INTERVAL_ROWS`: synthetic interval rows, default `250000`
- `VARIANTKEY_JOIN_RUNS`: timed repeats, default `3`
- `VARIANTKEY_JOIN_THREADS_LIST`: DuckDB thread grid, default `1,2,4`
- `DUCKHTS_EXTENSION`: optional extension path override

## Benchmark settings

``` r
benchmark_settings
#>               setting   value
#> 1   variant_join_rows 1000000
#> 2       interval_rows  250000
#> 3      benchmark_runs       3
#> 4 duckdb_threads_grid   1,2,4
#> 5 mixed_hash_fraction     20%
```

``` r
benchmark_environment <- data.frame(
  field = c("source_revision", "cpu", "process_affinity", "duckdb_version"),
  value = c(
    source_revision,
    cpu_model,
    cpu_affinity,
    DBI::dbGetQuery(con, "SELECT version() AS v")$v[1L]
  ),
  stringsAsFactors = FALSE
)
benchmark_environment
#>              field                                     value
#> 1  source_revision  83554766f5181c22625a53b2f2cbb0863f428204
#> 2              cpu       13th Gen Intel(R) Core(TM) i5-13500
#> 3 process_affinity pid 3914165's current affinity list: 0-19
#> 4   duckdb_version                                    v1.5.5
```

## 1. Exact reversible VariantKey joins

This section benchmarks exact self-joins over normalized reversible
variants.

``` r
reversible_multicol <- measure_thread_grid(
  "reversible_multicolumn_join",
  join_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM vk_probe_reversible p",
        "JOIN vk_ann_reversible a",
        "  ON p.chrom = a.chrom AND p.pos = a.pos AND p.ref = a.ref AND p.alt = a.alt"
      )
    )
  }
)

reversible_stored_vk <- measure_thread_grid(
  "reversible_stored_vk_join",
  join_rows,
  function() {
    DBI::dbGetQuery(
      con,
      "SELECT count(*) AS n FROM vk_probe_reversible p JOIN vk_ann_reversible a USING (vk)"
    )
  }
)

reversible_on_the_fly_vk <- measure_thread_grid(
  "reversible_on_the_fly_vk_join",
  join_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM vk_probe_reversible p",
        "JOIN vk_ann_reversible a",
        "  ON variantkey(p.chrom, p.pos, p.ref, p.alt) = a.vk"
      )
    )
  }
)

reversible_join_results <- rbind(
  reversible_multicol,
  reversible_stored_vk,
  reversible_on_the_fly_vk
)
reversible_join_results
#>                       operation threads input_rows result_rows median_seconds
#> 1   reversible_multicolumn_join       1    1000000       1e+06          0.206
#> 2   reversible_multicolumn_join       2    1000000       1e+06          0.117
#> 3   reversible_multicolumn_join       4    1000000       1e+06          0.071
#> 4     reversible_stored_vk_join       1    1000000       1e+06          0.050
#> 5     reversible_stored_vk_join       2    1000000       1e+06          0.037
#> 6     reversible_stored_vk_join       4    1000000       1e+06          0.032
#> 7 reversible_on_the_fly_vk_join       1    1000000       1e+06          0.092
#> 8 reversible_on_the_fly_vk_join       2    1000000       1e+06          0.053
#> 9 reversible_on_the_fly_vk_join       4    1000000       1e+06          0.041
#>   input_rows_per_second result_rows_per_second
#> 1               4854369                4854369
#> 2               8547009                8547009
#> 3              14084507               14084507
#> 4              20000000               20000000
#> 5              27027027               27027027
#> 6              31250000               31250000
#> 7              10869565               10869565
#> 8              18867925               18867925
#> 9              24390244               24390244
```

## 2. Mixed joins with hashed-row refinement

This section benchmarks the safe join policy for mixed data:

- join on stored `vk`
- refine only when the probe row is hashed / nonreversible

``` r
mixed_multicol <- measure_thread_grid(
  "mixed_multicolumn_join",
  join_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM vk_probe_mixed p",
        "JOIN vk_ann_mixed a",
        "  ON p.chrom = a.chrom AND p.pos = a.pos AND p.ref = a.ref AND p.alt = a.alt"
      )
    )
  }
)

mixed_stored_vk <- measure_thread_grid(
  "mixed_stored_vk_join_only",
  join_rows,
  function() {
    DBI::dbGetQuery(
      con,
      "SELECT count(*) AS n FROM vk_probe_mixed p JOIN vk_ann_mixed a USING (vk)"
    )
  }
)

mixed_safe_refine <- measure_thread_grid(
  "mixed_stored_vk_join_plus_hash_refine",
  join_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM vk_probe_mixed p",
        "JOIN vk_ann_mixed a",
        "  ON p.vk = a.vk",
        " AND (NOT p.is_hash OR (p.chrom = a.chrom AND p.pos = a.pos AND p.ref = a.ref AND p.alt = a.alt))"
      )
    )
  }
)

mixed_join_results <- rbind(
  mixed_multicol,
  mixed_stored_vk,
  mixed_safe_refine
)
mixed_join_results
#>                               operation threads input_rows result_rows
#> 1                mixed_multicolumn_join       1    1000000       1e+06
#> 2                mixed_multicolumn_join       2    1000000       1e+06
#> 3                mixed_multicolumn_join       4    1000000       1e+06
#> 4             mixed_stored_vk_join_only       1    1000000       1e+06
#> 5             mixed_stored_vk_join_only       2    1000000       1e+06
#> 6             mixed_stored_vk_join_only       4    1000000       1e+06
#> 7 mixed_stored_vk_join_plus_hash_refine       1    1000000       1e+06
#> 8 mixed_stored_vk_join_plus_hash_refine       2    1000000       1e+06
#> 9 mixed_stored_vk_join_plus_hash_refine       4    1000000       1e+06
#>   median_seconds input_rows_per_second result_rows_per_second
#> 1          0.199               5025126                5025126
#> 2          0.114               8771930                8771930
#> 3          0.074              13513514               13513514
#> 4          0.050              20000000               20000000
#> 5          0.035              28571429               28571429
#> 6          0.029              34482759               34482759
#> 7          0.143               6993007                6993007
#> 8          0.088              11363636               11363636
#> 9          0.058              17241379               17241379
```

## 3. Exact RegionKey span joins

Exact equality is only one RegionKey operation. RegionKey also sorts by
chromosome and start coordinate, can be binary-searched by that prefix,
and provides exact overlap predicates. The next section measures those
properties without claiming that a bare 64-bit key is an augmented
interval tree.

``` r
region_multicol <- measure_thread_grid(
  "region_multicolumn_join",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM rk_probe_eq p",
        "JOIN rk_ann a",
        "  ON p.chrom = a.chrom AND p.start = a.start AND p.\"end\" = a.\"end\""
      )
    )
  }
)

region_rk <- measure_thread_grid(
  "region_stored_rk_join",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      "SELECT count(*) AS n FROM rk_probe_eq p JOIN rk_ann a USING (rk)"
    )
  }
)

region_join_results <- rbind(region_multicol, region_rk)
region_join_results
#>                 operation threads input_rows result_rows median_seconds
#> 1 region_multicolumn_join       1     250000      250000          0.024
#> 2 region_multicolumn_join       2     250000      250000          0.015
#> 3 region_multicolumn_join       4     250000      250000          0.017
#> 4   region_stored_rk_join       1     250000      250000          0.011
#> 5   region_stored_rk_join       2     250000      250000          0.009
#> 6   region_stored_rk_join       4     250000      250000          0.009
#>   input_rows_per_second result_rows_per_second
#> 1              10416667               10416667
#> 2              16666667               16666667
#> 3              14705882               14705882
#> 4              22727273               22727273
#> 5              27777778               27777778
#> 6              27777778               27777778
```

## 4. Interval overlap implementations

The subject intervals have a declared maximum span of 30 bases. That
lets the RegionKey query derive a correct lower start bound,
range-search the sortable chromosome/start prefix, and apply
`are_overlapping_regionkeys(...)` to the candidate rows. The 30-base
maximum subject span is part of the synthetic input contract; without
it, the bounded query would miss intervals that start before the query.

The same probes and subjects are also evaluated through two vanilla
DuckDB forms and the immutable cgranges index. Because this isolated
distribution has one chromosome, the direct start/end query needs only
the two inequalities and plans as DuckDB’s `IE_JOIN`. Packing chromosome
and position into the RegionKey-compatible 33-bit coordinate generalizes
that two-range plan across chromosomes. All pair-producing methods must
return the same overlap count.

``` r
duckdb_interval_overlap <- measure_thread_grid(
  "duckdb_single_contig_interval_iejoin_pairs",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM ov_probe p",
        "JOIN ov_subject s",
        "  ON p.start < s.\"end\"",
        " AND p.\"end\" > s.start"
      )
    )
  }
)

duckdb_packed_interval_overlap <- measure_thread_grid(
  "duckdb_packed_coordinate_iejoin_pairs",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM ov_probe p",
        "JOIN ov_subject s",
        "  ON p.rk_chrom_start < s.rk_chrom_end",
        " AND p.rk_chrom_end > s.rk_chrom_start"
      )
    )
  }
)

regionkey_bounded_overlap <- measure_thread_grid(
  "regionkey_bounded_start_range_pairs",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n",
        "FROM ov_probe p",
        "JOIN ov_subject s",
        "  ON s.rk_chrom_start >=",
        "     (regionkey(p.chrom, greatest(0, p.start - 30 + 1)::BIGINT, greatest(0, p.start - 30 + 1)::BIGINT) >> 31)",
        " AND s.rk_chrom_start < (regionkey(p.chrom, p.\"end\"::BIGINT, p.\"end\"::BIGINT) >> 31)",
        " AND are_overlapping_regionkeys(p.rk, s.rk)"
      )
    )
  }
)
```

``` r
cgranges_has_overlap <- measure_thread_grid(
  "cgranges_has_overlap_filter",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM ov_probe p",
          "WHERE duckhts_cgranges_has_overlap('%s', p.chrom, p.start, p.\"end\")"
        ),
        idx_name
      )
    )
  }
)

cgranges_count_overlap <- measure_thread_grid(
  "cgranges_count_overlaps_sum",
  interval_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT sum(duckhts_cgranges_count_overlaps('%s', p.chrom, p.start::BIGINT, p.\"end\"::BIGINT))::BIGINT AS n",
          "FROM ov_probe p"
        ),
        idx_name
      )
    )
  }
)

cgranges_overlap_results <- rbind(cgranges_has_overlap, cgranges_count_overlap)
cgranges_overlap_results
#>                     operation threads input_rows result_rows median_seconds
#> 1 cgranges_has_overlap_filter       1     250000      250000          0.024
#> 2 cgranges_has_overlap_filter       2     250000      250000          0.013
#> 3 cgranges_has_overlap_filter       4     250000      250000          0.014
#> 4 cgranges_count_overlaps_sum       1     250000      250000          0.029
#> 5 cgranges_count_overlaps_sum       2     250000      250000          0.016
#> 6 cgranges_count_overlaps_sum       4     250000      250000          0.015
#>   input_rows_per_second result_rows_per_second
#> 1              10416667               10416667
#> 2              19230769               19230769
#> 3              17857143               17857143
#> 4               8620690                8620690
#> 5              15625000               15625000
#> 6              16666667               16666667
```

``` r
overlap_pair_results <- rbind(
  duckdb_interval_overlap,
  duckdb_packed_interval_overlap,
  regionkey_bounded_overlap,
  cgranges_count_overlap
)

stopifnot(length(unique(overlap_pair_results$result_rows)) == 1L)
overlap_pair_results
#>                                     operation threads input_rows result_rows
#> 1  duckdb_single_contig_interval_iejoin_pairs       1     250000      250000
#> 2  duckdb_single_contig_interval_iejoin_pairs       2     250000      250000
#> 3  duckdb_single_contig_interval_iejoin_pairs       4     250000      250000
#> 4       duckdb_packed_coordinate_iejoin_pairs       1     250000      250000
#> 5       duckdb_packed_coordinate_iejoin_pairs       2     250000      250000
#> 6       duckdb_packed_coordinate_iejoin_pairs       4     250000      250000
#> 7         regionkey_bounded_start_range_pairs       1     250000      250000
#> 8         regionkey_bounded_start_range_pairs       2     250000      250000
#> 9         regionkey_bounded_start_range_pairs       4     250000      250000
#> 10                cgranges_count_overlaps_sum       1     250000      250000
#> 11                cgranges_count_overlaps_sum       2     250000      250000
#> 12                cgranges_count_overlaps_sum       4     250000      250000
#>    median_seconds input_rows_per_second result_rows_per_second
#> 1           0.115               2173913                2173913
#> 2           0.063               3968254                3968254
#> 3           0.041               6097561                6097561
#> 4           0.114               2192982                2192982
#> 5           0.062               4032258                4032258
#> 6           0.043               5813953                5813953
#> 7           0.131               1908397                1908397
#> 8           0.080               3125000                3125000
#> 9           0.055               4545455                4545455
#> 10          0.029               8620690                8620690
#> 11          0.016              15625000               15625000
#> 12          0.015              16666667               16666667
```
