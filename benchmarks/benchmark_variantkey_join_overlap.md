DuckHTS VariantKey Join and Overlap Benchmark
================

<!-- benchmark_variantkey_join_overlap.md is generated from benchmark_variantkey_join_overlap.Rmd. -->

# Benchmark

This report focuses on the workloads that matter for a future
DuckVEP-style annotation layer:

1.  exact variant joins with stored `VariantKey`
2.  safe mixed joins where hashed/nonreversible rows refine on exact
    alleles
3.  exact span joins with stored `RegionKey`
4.  interval overlap through RegionKey-ordered candidate ranges, DuckDB
    range joins, and `duckhts_cgranges_*`

This is intentionally separate from
`benchmark_variantkey_conformance.Rmd`. That report measures encoding
and `%VKX` parity. This one measures join and overlap throughput.

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
#>              field                                    value
#> 1  source_revision fab39c5a6faa52211bdfd70855b24ca2832dfd0a
#> 2              cpu      13th Gen Intel(R) Core(TM) i5-13500
#> 3 process_affinity pid 151766's current affinity list: 0-19
#> 4   duckdb_version                                   v1.5.3
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
#> 1   reversible_multicolumn_join       1    1000000       1e+06          0.190
#> 2   reversible_multicolumn_join       2    1000000       1e+06          0.109
#> 3   reversible_multicolumn_join       4    1000000       1e+06          0.067
#> 4     reversible_stored_vk_join       1    1000000       1e+06          0.041
#> 5     reversible_stored_vk_join       2    1000000       1e+06          0.034
#> 6     reversible_stored_vk_join       4    1000000       1e+06          0.026
#> 7 reversible_on_the_fly_vk_join       1    1000000       1e+06          0.080
#> 8 reversible_on_the_fly_vk_join       2    1000000       1e+06          0.050
#> 9 reversible_on_the_fly_vk_join       4    1000000       1e+06          0.043
#>   input_rows_per_second result_rows_per_second
#> 1               5263158                5263158
#> 2               9174312                9174312
#> 3              14925373               14925373
#> 4              24390244               24390244
#> 5              29411765               29411765
#> 6              38461538               38461538
#> 7              12500000               12500000
#> 8              20000000               20000000
#> 9              23255814               23255814
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
#> 1          0.195               5128205                5128205
#> 2          0.108               9259259                9259259
#> 3          0.066              15151515               15151515
#> 4          0.042              23809524               23809524
#> 5          0.030              33333333               33333333
#> 6          0.030              33333333               33333333
#> 7          0.137               7299270                7299270
#> 8          0.081              12345679               12345679
#> 9          0.053              18867925               18867925
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
#> 1 region_multicolumn_join       1     250000      250000          0.023
#> 2 region_multicolumn_join       2     250000      250000          0.016
#> 3 region_multicolumn_join       4     250000      250000          0.014
#> 4   region_stored_rk_join       1     250000      250000          0.011
#> 5   region_stored_rk_join       2     250000      250000          0.008
#> 6   region_stored_rk_join       4     250000      250000          0.010
#>   input_rows_per_second result_rows_per_second
#> 1              10869565               10869565
#> 2              15625000               15625000
#> 3              17857143               17857143
#> 4              22727273               22727273
#> 5              31250000               31250000
#> 6              25000000               25000000
```

## 4. Interval overlap implementations

The subject intervals have a declared maximum span of 30 bases. That
lets the RegionKey query derive a correct lower start bound,
range-search the sortable chromosome/start prefix, and apply
`are_overlapping_regionkeys(...)` to the candidate rows. The maximum
span is part of the index receipt: omitting it would miss long intervals
that start before the query.

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
        "     (regionkey(p.chrom, greatest(0, p.start - 30 + 1), greatest(0, p.start - 30 + 1)) >> 31)",
        " AND s.rk_chrom_start < (regionkey(p.chrom, p.\"end\", p.\"end\") >> 31)",
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
          "SELECT sum(duckhts_cgranges_count_overlaps('%s', p.chrom, p.start, p.\"end\"))::BIGINT AS n",
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
#> 1 cgranges_has_overlap_filter       1     250000      250000          0.023
#> 2 cgranges_has_overlap_filter       2     250000      250000          0.013
#> 3 cgranges_has_overlap_filter       4     250000      250000          0.012
#> 4 cgranges_count_overlaps_sum       1     250000      250000          0.028
#> 5 cgranges_count_overlaps_sum       2     250000      250000          0.014
#> 6 cgranges_count_overlaps_sum       4     250000      250000          0.014
#>   input_rows_per_second result_rows_per_second
#> 1              10869565               10869565
#> 2              19230769               19230769
#> 3              20833333               20833333
#> 4               8928571                8928571
#> 5              17857143               17857143
#> 6              17857143               17857143
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
#> 1           0.111               2252252                2252252
#> 2           0.061               4098361                4098361
#> 3           0.041               6097561                6097561
#> 4           0.108               2314815                2314815
#> 5           0.063               3968254                3968254
#> 6           0.044               5681818                5681818
#> 7           0.125               2000000                2000000
#> 8           0.077               3246753                3246753
#> 9           0.053               4716981                4716981
#> 10          0.028               8928571                8928571
#> 11          0.014              17857143               17857143
#> 12          0.014              17857143               17857143
```

## External provider composition

Real supplementary-provider joins are intentionally outside this
benchmark. The former report depended on an undeclared local mixture of
ClinVar, GIAB, Ensembl, gnomAD, REVEL, AlphaMissense, and
ClinvArbitration files. A result cannot be portable merely because paths
are environment variables. This workload returns only after one explicit
cache staging command declares every source release, access requirement,
transformation, cached output, and consumer variable.

## Notes

- The exact-join sections are the relevant comparison class for future
  `echtvar`-style annotation throughput work.
- The `mixed_stored_vk_join_plus_hash_refine` query is the safe SQL
  pattern for nonreversible VariantKey rows: equal keys still refine on
  exact `REF/ALT`.
- VariantKey does not normalize alleles. A real provider comparison must
  name one normalized biallelic representation before timing a join.
- RegionKey’s chromosome/start ordering supports bounded range
  discovery; a general nested-interval workload needs a complete overlap
  index such as cgranges or an explicit maximum-end structure.
- DuckDB’s interval join is the zero-registry SQL baseline. The
  RegionKey range query and cgranges have different preparation and
  repeated-query trade-offs. A staged real workload, not the key format
  alone, decides between them.
