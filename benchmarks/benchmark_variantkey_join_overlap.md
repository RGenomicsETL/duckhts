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
4.  true interval overlap queries through `duckhts_cgranges_*`

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
#>                       operation threads    rows median_seconds rows_per_second
#> 1   reversible_multicolumn_join       1 1000000          0.186         5376344
#> 2   reversible_multicolumn_join       2 1000000          0.106         9433962
#> 3   reversible_multicolumn_join       4 1000000          0.072        13888889
#> 4     reversible_stored_vk_join       1 1000000          0.048        20833333
#> 5     reversible_stored_vk_join       2 1000000          0.030        33333333
#> 6     reversible_stored_vk_join       4 1000000          0.024        41666667
#> 7 reversible_on_the_fly_vk_join       1 1000000          0.078        12820513
#> 8 reversible_on_the_fly_vk_join       2 1000000          0.051        19607843
#> 9 reversible_on_the_fly_vk_join       4 1000000          0.036        27777778
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
#>                               operation threads    rows median_seconds
#> 1                mixed_multicolumn_join       1 1000000          0.185
#> 2                mixed_multicolumn_join       2 1000000          0.105
#> 3                mixed_multicolumn_join       4 1000000          0.063
#> 4             mixed_stored_vk_join_only       1 1000000          0.042
#> 5             mixed_stored_vk_join_only       2 1000000          0.029
#> 6             mixed_stored_vk_join_only       4 1000000          0.025
#> 7 mixed_stored_vk_join_plus_hash_refine       1 1000000          0.133
#> 8 mixed_stored_vk_join_plus_hash_refine       2 1000000          0.080
#> 9 mixed_stored_vk_join_plus_hash_refine       4 1000000          0.049
#>   rows_per_second
#> 1         5405405
#> 2         9523810
#> 3        15873016
#> 4        23809524
#> 5        34482759
#> 6        40000000
#> 7         7518797
#> 8        12500000
#> 9        20408163
```

## 3. Exact RegionKey span joins

`RegionKey` is useful for exact stored span equality. It is not the
general interval-overlap engine; that comes next.

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
#>                 operation threads   rows median_seconds rows_per_second
#> 1 region_multicolumn_join       1 250000          0.023        10869565
#> 2 region_multicolumn_join       2 250000          0.015        16666667
#> 3 region_multicolumn_join       4 250000          0.017        14705882
#> 4   region_stored_rk_join       1 250000          0.010        25000000
#> 5   region_stored_rk_join       2 250000          0.008        31250000
#> 6   region_stored_rk_join       4 250000          0.009        27777778
```

## 4. cgranges overlap probes

This is the actual overlap path of interest for
transcript/gene/regulatory annotation candidate generation.

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
#>                     operation threads   rows median_seconds rows_per_second
#> 1 cgranges_has_overlap_filter       1 250000          0.023        10869565
#> 2 cgranges_has_overlap_filter       2 250000          0.013        19230769
#> 3 cgranges_has_overlap_filter       4 250000          0.014        17857143
#> 4 cgranges_count_overlaps_sum       1 250000          0.027         9259259
#> 5 cgranges_count_overlaps_sum       2 250000          0.014        17857143
#> 6 cgranges_count_overlaps_sum       4 250000          0.014        17857143
```

## Notes

- The exact-join sections are the relevant comparison class for future
  `echtvar`-style annotation throughput work.
- The `mixed_stored_vk_join_plus_hash_refine` query represents the safe
  SQL pattern for hashed/nonreversible VariantKey rows when exact
  `REF/ALT` verification is required.
- `RegionKey` equality is useful for exact stored spans, but real
  overlap work should go through `duckhts_cgranges_*` or later overlap
  kernels rather than a plain `rk = rk` join.
- For richer overlap comparisons against external tools such as `bedtk`
  or `bedtools`, see the existing `benchmarks/benchmark_cgranges.Rmd`
  report.
