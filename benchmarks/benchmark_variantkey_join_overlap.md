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
- `VARIANTKEY_CLINVAR_VCF`: dated ClinVar GRCh38 VCF
- `VARIANTKEY_GIAB_VCF`: GIAB HG002 GRCh38 benchmark VCF
- `VARIANTKEY_REGULATORY_PARQUET`: Ensembl regulatory relation
- `VARIANTKEY_TRANSCRIPTS_PARQUET`: Ensembl transcript relation
- `VARIANTKEY_GNOMAD_CONSTRAINT`: gnomAD gene-constraint table
- `VARIANTKEY_REAL_CACHE_DIR`: generated benchmark Parquet directory

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
#>              field                                        value
#> 1  source_revision     ccf6971dac948fb0d6d7bb88a53ff3079a223c77
#> 2              cpu          13th Gen Intel(R) Core(TM) i5-13500
#> 3 process_affinity pid 3475186's current affinity list: 0,2,4,6
#> 4   duckdb_version                                       v1.5.3

source_receipts <- data.frame(
  source = names(real_source_paths),
  path = unname(real_source_paths),
  bytes = unname(file.info(real_source_paths)$size),
  sha256 = vapply(real_source_paths, sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
source_receipts
#>                                                            source
#> clinvar_vcf                                           clinvar_vcf
#> giab_vcf                                                 giab_vcf
#> regulatory_parquet                             regulatory_parquet
#> transcripts_parquet                           transcripts_parquet
#> gnomad_constraint                               gnomad_constraint
#> revel_v13_derived_parquet               revel_v13_derived_parquet
#> alphamissense_v2_hg38                       alphamissense_v2_hg38
#> clinvarbitration_zenodo_16792026 clinvarbitration_zenodo_16792026
#>                                                                                                                        path
#> clinvar_vcf                                             /root/duckvep/data/corpora/clinvar/20260706/clinvar_20260706.vcf.gz
#> giab_vcf                                             /root/duckvep/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
#> regulatory_parquet                                          /root/duckvep/data/cache/homo_sapiens.116.38/regulatory.parquet
#> transcripts_parquet                                        /root/duckvep/data/cache/homo_sapiens.116.38/transcripts.parquet
#> gnomad_constraint                                                /root/bioconnect-sprint-py/.cache/gnomad_constraint.txt.gz
#> revel_v13_derived_parquet                                            /root/bioconnect-sprint-py/.cache/revel_grch37.parquet
#> alphamissense_v2_hg38                                  /root/duckhts/.tmp/variantkey_join_overlap/AlphaMissense_hg38.tsv.gz
#> clinvarbitration_zenodo_16792026 /root/.cache/RClinVarbitration/zenodo/16792026/clinvarbitration_data/clinvar_decisions.tsv
#>                                      bytes
#> clinvar_vcf                      192290992
#> giab_vcf                         156252944
#> regulatory_parquet                 9513375
#> transcripts_parquet               13808698
#> gnomad_constraint                  4609488
#> revel_v13_derived_parquet        171593942
#> alphamissense_v2_hg38            642961469
#> clinvarbitration_zenodo_16792026 131533304
#>                                                                                            sha256
#> clinvar_vcf                      59a83b34d425daf58cd0dd463d6f2952f0a833ddf8fe6698fd30010642e5e1e9
#> giab_vcf                         adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c
#> regulatory_parquet               86ec412b252b446235d34092a82b37c9cff91aee972e035c8f1f5f83d374aa57
#> transcripts_parquet              90510fe83d9a96895ee7405aae1e9a211154d7dda14869e8669bd346112ef521
#> gnomad_constraint                153031d34b6794e8e99eb0306bc3c50b13b18accda8b0ffef91c2623dd3affd5
#> revel_v13_derived_parquet        1e7bfb57344be97dcf644fb2fc23317263e26a104f4b9b53d4d918610199a535
#> alphamissense_v2_hg38            0516cfd71c0767ac8f9c469252d429000e94e02c008b6e3a46d4b4646fcd3475
#> clinvarbitration_zenodo_16792026 3e28bf60f3c934887803877b67dde60cfd8dd03a4617a80b106696cd1198db7c
```

AlphaMissense is the official v2 GRCh38 genomic file from [Zenodo record
8360242](https://zenodo.org/records/8360242); its published MD5 is
`9fd167735f16a1b87da6eb3e4c25fcb5`. REVEL is the locally derived,
column-pruned Parquet form of the official [REVEL v1.3
download](https://sites.google.com/site/revelgenomics/downloads). The
ClinvArbitration profile is the extracted TSV from [Zenodo record
16792026](https://zenodo.org/records/16792026). That TSV does not carry
an assembly column; the ClinvArbitration/Talos release contract declares
GRCh38, and the provider receipt below supplies that semantic input
before keying. The fact does not become inferable from the coordinate
columns.

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
#> 1   reversible_multicolumn_join       1    1000000       1e+06          0.188
#> 2   reversible_multicolumn_join       2    1000000       1e+06          0.106
#> 3   reversible_multicolumn_join       4    1000000       1e+06          0.066
#> 4     reversible_stored_vk_join       1    1000000       1e+06          0.042
#> 5     reversible_stored_vk_join       2    1000000       1e+06          0.030
#> 6     reversible_stored_vk_join       4    1000000       1e+06          0.028
#> 7 reversible_on_the_fly_vk_join       1    1000000       1e+06          0.087
#> 8 reversible_on_the_fly_vk_join       2    1000000       1e+06          0.052
#> 9 reversible_on_the_fly_vk_join       4    1000000       1e+06          0.036
#>   input_rows_per_second result_rows_per_second
#> 1               5319149                5319149
#> 2               9433962                9433962
#> 3              15151515               15151515
#> 4              23809524               23809524
#> 5              33333333               33333333
#> 6              35714286               35714286
#> 7              11494253               11494253
#> 8              19230769               19230769
#> 9              27777778               27777778
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
#> 1          0.191               5235602                5235602
#> 2          0.108               9259259                9259259
#> 3          0.062              16129032               16129032
#> 4          0.040              25000000               25000000
#> 5          0.030              33333333               33333333
#> 6          0.024              41666667               41666667
#> 7          0.135               7407407                7407407
#> 8          0.080              12500000               12500000
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
#> 2 region_multicolumn_join       2     250000      250000          0.015
#> 3 region_multicolumn_join       4     250000      250000          0.016
#> 4   region_stored_rk_join       1     250000      250000          0.011
#> 5   region_stored_rk_join       2     250000      250000          0.008
#> 6   region_stored_rk_join       4     250000      250000          0.009
#>   input_rows_per_second result_rows_per_second
#> 1              10869565               10869565
#> 2              16666667               16666667
#> 3              15625000               15625000
#> 4              22727273               22727273
#> 5              31250000               31250000
#> 6              27777778               27777778
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
#> 1 cgranges_has_overlap_filter       1     250000      250000          0.025
#> 2 cgranges_has_overlap_filter       2     250000      250000          0.013
#> 3 cgranges_has_overlap_filter       4     250000      250000          0.013
#> 4 cgranges_count_overlaps_sum       1     250000      250000          0.028
#> 5 cgranges_count_overlaps_sum       2     250000      250000          0.015
#> 6 cgranges_count_overlaps_sum       4     250000      250000          0.015
#>   input_rows_per_second result_rows_per_second
#> 1              10000000               10000000
#> 2              19230769               19230769
#> 3              19230769               19230769
#> 4               8928571                8928571
#> 5              16666667               16666667
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
#> 1           0.109               2293578                2293578
#> 2           0.061               4098361                4098361
#> 3           0.043               5813953                5813953
#> 4           0.108               2314815                2314815
#> 5           0.061               4098361                4098361
#> 6           0.043               5813953                5813953
#> 7           0.126               1984127                1984127
#> 8           0.077               3246753                3246753
#> 9           0.049               5102041                5102041
#> 10          0.028               8928571                8928571
#> 11          0.015              16666667               16666667
#> 12          0.015              16666667               16666667
```

## 5. Real supplementary-annotation composition

The synthetic sections isolate operator costs. This section exercises
the storage and join architecture on real release data:

- 2026-07-06 ClinVar GRCh38 as an exact and positional provider;
- GIAB HG002 GRCh38 v4.2.1 as a whole-genome probe relation;
- Ensembl 116 regulatory features as the interval provider; and
- gnomAD v2.1.1 constraint as a gene provider joined to the full Ensembl
  116 transcript relation;
- REVEL v1.3’s 77.97 million GRCh37 SNV scores; and
- the official AlphaMissense v2 GRCh38 genomic release’s 71.70 million
  rows.

ClinvArbitration Zenodo record 16792026 is also exercised as a GRCh38
exact provider. Its TSV has 3.65 million rows and omits an assembly
column; GRCh38 is therefore an explicit receipt fact, not a guessed
property of `chr`-prefixed coordinates.

VCF is decoded once into narrow, typed, ZSTD Parquet. VariantKey and
RegionKey are prepared columns, not recomputed for every downstream
provider join. The Parquet files are generated benchmark artifacts and
are not committed.

``` r
real_case_path <- file.path(real_cache_dir, "giab_hg002_grch38_keyed.parquet")
real_clinvar_path <- file.path(real_cache_dir, "clinvar_20260706_grch38_keyed.parquet")
real_regulatory_path <- file.path(real_cache_dir, "ensembl116_grch38_regulatory_keyed.parquet")
real_constraint_path <- file.path(real_cache_dir, "gnomad_v211_constraint_gene.parquet")
real_case_reversible_path <- file.path(real_cache_dir, "giab_hg002_reversible.parquet")
real_case_hashed_path <- file.path(real_cache_dir, "giab_hg002_hashed.parquet")
real_clinvar_reversible_path <- file.path(real_cache_dir, "clinvar_20260706_reversible.parquet")
real_clinvar_hashed_path <- file.path(real_cache_dir, "clinvar_20260706_hashed.parquet")
real_revel_path <- file.path(real_cache_dir, "revel_grch37_variantkey.parquet")
real_revel_probe_path <- file.path(real_cache_dir, "revel_probe_hash19.parquet")
real_alphamissense_path <- file.path(real_cache_dir, "alphamissense_hg38_variantkey.parquet")
real_clinvarbitration_path <- file.path(real_cache_dir, "clinvarbitration_16792026_grch38_keyed.parquet")
real_clinvarbitration_reversible_path <- file.path(
  real_cache_dir,
  "clinvarbitration_16792026_reversible.parquet"
)
real_clinvarbitration_hashed_path <- file.path(
  real_cache_dir,
  "clinvarbitration_16792026_hashed.parquet"
)
real_clinvar_reversible_rg32k_path <- file.path(
  real_cache_dir,
  "clinvar_20260706_reversible_rg32768.parquet"
)
real_clinvar_reversible_rg1m_path <- file.path(
  real_cache_dir,
  "clinvar_20260706_reversible_rg1048576.parquet"
)

canonical_contigs_sql <- paste(
  sql_string(c(as.character(1:22), "X", "Y", "MT")),
  collapse = ", "
)

giab_stage_query <- sprintf(
  paste(
    "WITH alleles AS (",
    "  SELECT duckhts_contig_key(CHROM) AS chrom,",
    "         POS::BIGINT AS pos, REF AS ref, allele.alt AS alt",
    "  FROM read_bcf(%s, scan_mode := 'sequential')",
    "  CROSS JOIN UNNEST(ALT) AS allele(alt)",
    "), keyed AS (",
    "  SELECT alleles.*,",
    "         variantkey(chrom, pos, ref, alt) AS vk,",
    "         regionkey(chrom, pos - 1, pos - 1 + length(ref)) AS rk",
    "  FROM alleles",
    "  WHERE chrom IN (%s)",
    ")",
    "SELECT chrom, pos, ref, alt, vk,",
    "       mod(extract_variantkey_refalt(vk), 2) = 1 AS is_hash,",
    "       (rk >> 31) AS rk_chrom_start,",
    "       (regionkey(chrom, pos - 1 + length(ref), pos - 1 + length(ref)) >> 31) AS rk_chrom_end,",
    "       pos - 1 AS start0, pos - 1 + length(ref) AS end0, rk",
    "FROM keyed",
    "WHERE vk IS NOT NULL AND rk IS NOT NULL",
    "ORDER BY vk"
  ),
  sql_path(giab_vcf),
  canonical_contigs_sql
)

clinvar_stage_query <- sprintf(
  paste(
    "WITH alleles AS (",
    "  SELECT duckhts_contig_key(CHROM) AS chrom,",
    "         POS::BIGINT AS pos, REF AS ref, allele.alt AS alt,",
    "         INFO_ALLELEID AS allele_id, INFO_CLNSIG AS clinical_significance,",
    "         INFO_GENEINFO AS gene_info",
    "  FROM read_bcf(%s, scan_mode := 'sequential')",
    "  CROSS JOIN UNNEST(ALT) AS allele(alt)",
    "), keyed AS (",
    "  SELECT alleles.*, variantkey(chrom, pos, ref, alt) AS vk",
    "  FROM alleles",
    "  WHERE chrom IN (%s)",
    ")",
    "SELECT chrom, pos, ref, alt, vk,",
    "       mod(extract_variantkey_refalt(vk), 2) = 1 AS is_hash,",
    "       allele_id, clinical_significance, gene_info",
    "FROM keyed",
    "WHERE vk IS NOT NULL",
    "ORDER BY vk"
  ),
  sql_path(clinvar_vcf),
  canonical_contigs_sql
)

regulatory_stage_query <- sprintf(
  paste(
    "WITH source AS (",
    "  SELECT stable_id, duckhts_contig_key(chrom) AS chrom,",
    "         start - 1 AS start0, \"end\" AS end0,",
    "         feature_type, so_term, so_accession",
    "  FROM read_parquet(%s)",
    "), keyed AS (",
    "  SELECT source.*, regionkey(chrom, start0, end0) AS rk",
    "  FROM source",
    "  WHERE chrom IN (%s)",
    ")",
    "SELECT keyed.*, rk >> 31 AS rk_chrom_start,",
    "       regionkey(chrom, end0, end0) >> 31 AS rk_chrom_end",
    "FROM keyed",
    "WHERE rk IS NOT NULL",
    "ORDER BY rk"
  ),
  sql_path(regulatory_parquet),
  canonical_contigs_sql
)

constraint_stage_query <- sprintf(
  paste(
    "SELECT gene_id, gene AS gene_symbol, transcript,",
    "       try_cast(pLI AS DOUBLE) AS pLI,",
    "       try_cast(oe_lof_upper AS DOUBLE) AS oe_lof_upper",
    "FROM read_csv_auto(%s, delim = '\\t', header = TRUE, nullstr = 'NA')",
    "WHERE gene_id IS NOT NULL"
  ),
  sql_path(gnomad_constraint)
)

revel_stage_query <- sprintf(
  paste(
    "WITH keyed AS (",
    "  SELECT variantkey(chrom, pos, ref, alt) AS vk, revel",
    "  FROM read_parquet(%s)",
    ")",
    "SELECT vk, revel FROM keyed WHERE vk IS NOT NULL"
  ),
  sql_path(revel_parquet)
)

alphamissense_columns_sql <- paste0(
  "{'chrom':'VARCHAR','pos':'BIGINT','ref':'VARCHAR','alt':'VARCHAR',",
  "'genome':'VARCHAR','uniprot_id':'VARCHAR','transcript_id':'VARCHAR',",
  "'protein_variant':'VARCHAR','am_pathogenicity':'DOUBLE','am_class':'VARCHAR'}"
)
alphamissense_stage_query <- sprintf(
  paste(
    "WITH keyed AS (",
    "  SELECT variantkey(duckhts_contig_key(chrom), pos, ref, alt) AS vk,",
    "         am_pathogenicity, am_class",
    "  FROM read_csv(%s, delim = '\\t', skip = 4, header = FALSE, columns = %s)",
    ")",
    "SELECT vk, am_pathogenicity, am_class FROM keyed WHERE vk IS NOT NULL"
  ),
  sql_path(alphamissense_gz),
  alphamissense_columns_sql
)

clinvarbitration_stage_query <- sprintf(
  paste(
    "WITH source AS (",
    "  SELECT duckhts_contig_key(contig) AS chrom,",
    "         try_cast(position AS BIGINT) AS pos,",
    "         reference AS ref, alternate AS alt,",
    "         clinical_significance, try_cast(gold_stars AS UTINYINT) AS gold_stars,",
    "         try_cast(allele_id AS BIGINT) AS allele_id",
    "  FROM read_csv_auto(%s, delim = '\\t', header = TRUE, all_varchar = TRUE)",
    "), keyed AS (",
    "  SELECT source.*, variantkey(chrom, pos, ref, alt) AS vk",
    "  FROM source",
    ")",
    "SELECT keyed.*, mod(extract_variantkey_refalt(vk), 2) = 1 AS is_hash",
    "FROM keyed WHERE vk IS NOT NULL"
  ),
  sql_path(clinvarbitration_tsv)
)

staging_results <- rbind(
  cbind(source = "GIAB HG002 VCF to keyed Parquet", copy_timed(giab_stage_query, real_case_path)),
  cbind(source = "ClinVar VCF to keyed Parquet", copy_timed(clinvar_stage_query, real_clinvar_path)),
  cbind(source = "Ensembl regulatory to keyed Parquet", copy_timed(regulatory_stage_query, real_regulatory_path)),
  cbind(source = "gnomAD constraint to gene Parquet", copy_timed(constraint_stage_query, real_constraint_path)),
  cbind(
    source = "REVEL v1.3 to VariantKey Parquet",
    copy_timed(revel_stage_query, real_revel_path, "FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880")
  ),
  cbind(
    source = "AlphaMissense v2 to VariantKey Parquet",
    copy_timed(
      alphamissense_stage_query,
      real_alphamissense_path,
      "FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880"
    )
  ),
  cbind(
    source = "ClinvArbitration 16792026 to keyed Parquet",
    copy_timed(
      clinvarbitration_stage_query,
      real_clinvarbitration_path,
      "FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880"
    )
  )
)
staging_results$rows_per_second <- staging_results$rows / staging_results$seconds
staging_results$output_mb_per_second <- staging_results$bytes / staging_results$seconds / 1e6
staging_results
#>                                       source
#> 1            GIAB HG002 VCF to keyed Parquet
#> 2               ClinVar VCF to keyed Parquet
#> 3        Ensembl regulatory to keyed Parquet
#> 4          gnomAD constraint to gene Parquet
#> 5           REVEL v1.3 to VariantKey Parquet
#> 6     AlphaMissense v2 to VariantKey Parquet
#> 7 ClinvArbitration 16792026 to keyed Parquet
#>                                                                                        path
#> 1                /root/duckhts/.tmp/variantkey_join_overlap/giab_hg002_grch38_keyed.parquet
#> 2          /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_grch38_keyed.parquet
#> 3     /root/duckhts/.tmp/variantkey_join_overlap/ensembl116_grch38_regulatory_keyed.parquet
#> 4            /root/duckhts/.tmp/variantkey_join_overlap/gnomad_v211_constraint_gene.parquet
#> 5                /root/duckhts/.tmp/variantkey_join_overlap/revel_grch37_variantkey.parquet
#> 6          /root/duckhts/.tmp/variantkey_join_overlap/alphamissense_hg38_variantkey.parquet
#> 7 /root/duckhts/.tmp/variantkey_join_overlap/clinvarbitration_16792026_grch38_keyed.parquet
#>       rows     bytes seconds rows_per_second output_mb_per_second
#> 1  4096123  90299692   6.996       585495.00           12.9073316
#> 2  4438467  39774942   7.672       578528.02            5.1844294
#> 3   643528   9958593   0.136      4731823.53           73.2249485
#> 4    19704    368158   0.474        41569.62            0.7767046
#> 5 77966138 201658023   2.272     34316081.87           88.7579327
#> 6 71697556 254370824  12.721      5636157.22           19.9961343
#> 7  3647840  30075314   0.566      6444947.00           53.1365972

clinvarbitration_profile <- data.frame(
  source = "ClinvArbitration Zenodo 16792026 TSV",
  rows = as.numeric(DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT count(*) AS n FROM read_csv_auto(%s, delim = '\\t', header = TRUE,",
        "all_varchar = TRUE, ignore_errors = TRUE)"
      ),
      sql_path(clinvarbitration_tsv)
    )
  )$n[1L]),
  bytes = file.info(clinvarbitration_tsv)$size,
  assembly_in_relation = FALSE,
  receipt_assembly = "GRCh38",
  stringsAsFactors = FALSE
)
clinvarbitration_profile
#>                                 source    rows     bytes assembly_in_relation
#> 1 ClinvArbitration Zenodo 16792026 TSV 3647840 131533304                FALSE
#>   receipt_assembly
#> 1           GRCh38

lane_staging_results <- rbind(
  cbind(
    source = "GIAB reversible lane",
    copy_timed(
      sprintf(
        "SELECT * EXCLUDE (is_hash) FROM read_parquet(%s) WHERE NOT is_hash ORDER BY vk",
        sql_path(real_case_path)
      ),
      real_case_reversible_path
    )
  ),
  cbind(
    source = "GIAB hashed lane",
    copy_timed(
      sprintf(
        "SELECT * EXCLUDE (is_hash) FROM read_parquet(%s) WHERE is_hash ORDER BY vk",
        sql_path(real_case_path)
      ),
      real_case_hashed_path
    )
  ),
  cbind(
    source = "ClinVar reversible numeric lane",
    copy_timed(
      sprintf(
        paste(
          "SELECT vk, allele_id, clinical_significance, gene_info",
          "FROM read_parquet(%s) WHERE NOT is_hash ORDER BY vk"
        ),
        sql_path(real_clinvar_path)
      ),
      real_clinvar_reversible_path
    )
  ),
  cbind(
    source = "ClinVar hashed refinement lane",
    copy_timed(
      sprintf(
        paste(
          "SELECT chrom, pos, ref, alt, vk, allele_id, clinical_significance, gene_info",
          "FROM read_parquet(%s) WHERE is_hash ORDER BY vk"
        ),
        sql_path(real_clinvar_path)
      ),
      real_clinvar_hashed_path
    )
  ),
  cbind(
    source = "ClinvArbitration reversible numeric lane",
    copy_timed(
      sprintf(
        paste(
          "SELECT vk, allele_id, clinical_significance, gold_stars",
          "FROM read_parquet(%s) WHERE NOT is_hash ORDER BY vk"
        ),
        sql_path(real_clinvarbitration_path)
      ),
      real_clinvarbitration_reversible_path
    )
  ),
  cbind(
    source = "ClinvArbitration hashed refinement lane",
    copy_timed(
      sprintf(
        paste(
          "SELECT chrom, pos, ref, alt, vk, allele_id, clinical_significance, gold_stars",
          "FROM read_parquet(%s) WHERE is_hash ORDER BY vk"
        ),
        sql_path(real_clinvarbitration_path)
      ),
      real_clinvarbitration_hashed_path
    )
  )
)
lane_staging_results$rows_per_second <- lane_staging_results$rows / lane_staging_results$seconds
lane_staging_results$output_mb_per_second <- lane_staging_results$bytes / lane_staging_results$seconds / 1e6
lane_staging_results
#>                                     source
#> 1                     GIAB reversible lane
#> 2                         GIAB hashed lane
#> 3          ClinVar reversible numeric lane
#> 4           ClinVar hashed refinement lane
#> 5 ClinvArbitration reversible numeric lane
#> 6  ClinvArbitration hashed refinement lane
#>                                                                                      path
#> 1                /root/duckhts/.tmp/variantkey_join_overlap/giab_hg002_reversible.parquet
#> 2                    /root/duckhts/.tmp/variantkey_join_overlap/giab_hg002_hashed.parquet
#> 3          /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible.parquet
#> 4              /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_hashed.parquet
#> 5 /root/duckhts/.tmp/variantkey_join_overlap/clinvarbitration_16792026_reversible.parquet
#> 6     /root/duckhts/.tmp/variantkey_join_overlap/clinvarbitration_16792026_hashed.parquet
#>      rows    bytes seconds rows_per_second output_mb_per_second
#> 1 4036258 88244417   0.531       7601239.2            166.18534
#> 2   59865  2091959   0.117        511666.7             17.87999
#> 3 4399443 27029899   0.275      15997974.5             98.29054
#> 4   39024  2825872   0.107        364710.3             26.41002
#> 5 3620269 21214375   0.204      17746416.7            103.99203
#> 6   27571   599913   0.056        492339.3             10.71273

revel_probe_staging <- copy_timed(
  sprintf(
    paste(
      "SELECT vk FROM read_parquet(%s)",
      "WHERE hash(vk) %% 19 = 0 ORDER BY vk"
    ),
    sql_path(real_revel_path)
  ),
  real_revel_probe_path
)
revel_probe_staging
#>                                                                    path    rows
#> 1 /root/duckhts/.tmp/variantkey_join_overlap/revel_probe_hash19.parquet 4103497
#>     bytes seconds
#> 1 9371710   0.384

row_group_staging_results <- rbind(
  cbind(
    row_group_size = 32768L,
    copy_timed(
      sprintf("SELECT * FROM read_parquet(%s)", sql_path(real_clinvar_reversible_path)),
      real_clinvar_reversible_rg32k_path,
      "FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 32768"
    )
  ),
  cbind(
    row_group_size = 1048576L,
    copy_timed(
      sprintf("SELECT * FROM read_parquet(%s)", sql_path(real_clinvar_reversible_path)),
      real_clinvar_reversible_rg1m_path,
      "FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1048576"
    )
  )
)

default_row_groups <- DBI::dbGetQuery(
  con,
  sprintf(
    "SELECT count(DISTINCT row_group_id)::BIGINT AS n FROM parquet_metadata(%s)",
    sql_path(real_clinvar_reversible_path)
  )
)$n[1L]
row_group_staging_results$row_groups <- vapply(
  row_group_staging_results$path,
  function(path) {
    as.numeric(DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT count(DISTINCT row_group_id)::BIGINT AS n FROM parquet_metadata(%s)",
        sql_path(path)
      )
    )$n[1L])
  },
  numeric(1L)
)
row_group_layouts <- rbind(
  data.frame(
    row_group_size = 122880L,
    path = real_clinvar_reversible_path,
    rows = lane_staging_results$rows[
      lane_staging_results$source == "ClinVar reversible numeric lane"
    ],
    bytes = file.info(real_clinvar_reversible_path)$size,
    seconds = NA_real_,
    row_groups = as.numeric(default_row_groups),
    stringsAsFactors = FALSE
  ),
  row_group_staging_results[, c(
    "row_group_size", "path", "rows", "bytes", "seconds", "row_groups"
  )]
)
row_group_layouts
#>   row_group_size
#> 1         122880
#> 2          32768
#> 3        1048576
#>                                                                                       path
#> 1           /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible.parquet
#> 2   /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible_rg32768.parquet
#> 3 /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible_rg1048576.parquet
#>      rows    bytes seconds row_groups
#> 1 4399443 27029899      NA         36
#> 2 4399443 26945479   0.275        135
#> 3 4399443 27004414   0.371          5
```

### Exact and positional ClinVar joins

The exact join uses the collision-safe policy: reversible rows join on
the 64-bit VariantKey alone; hashed/nonreversible rows additionally
refine on the normalized allele tuple. The multicolumn query is the
correctness oracle. The positional query intentionally has a different,
larger result cardinality.

``` r
real_case_rows <- staging_results$rows[staging_results$source == "GIAB HG002 VCF to keyed Parquet"]

real_exact_multicol <- measure_thread_grid(
  "real_clinvar_exact_multicolumn_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) a",
          "  ON q.chrom = a.chrom AND q.pos = a.pos",
          " AND q.ref = a.ref AND q.alt = a.alt"
        ),
        sql_path(real_case_path),
        sql_path(real_clinvar_path)
      )
    )
  }
)

real_exact_variantkey <- measure_thread_grid(
  "real_clinvar_exact_variantkey_refined_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) a",
          "  ON q.vk = a.vk",
          " AND (NOT (q.is_hash OR a.is_hash)",
          "      OR (q.chrom = a.chrom AND q.pos = a.pos AND q.ref = a.ref AND q.alt = a.alt))"
        ),
        sql_path(real_case_path),
        sql_path(real_clinvar_path)
      )
    )
  }
)

real_exact_variantkey_lanes <- measure_thread_grid(
  "real_clinvar_exact_variantkey_split_lanes_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste0(
          "SELECT sum(n)::BIGINT AS n FROM (",
          "SELECT count(*)::BIGINT AS n ",
          "FROM read_parquet(%s) q JOIN read_parquet(%s) a USING (vk) ",
          "UNION ALL ",
          "SELECT count(*)::BIGINT AS n ",
          "FROM read_parquet(%s) q JOIN read_parquet(%s) a ",
          "ON q.vk = a.vk AND q.chrom = a.chrom AND q.pos = a.pos ",
          "AND q.ref = a.ref AND q.alt = a.alt)"
        ),
        sql_path(real_case_reversible_path),
        sql_path(real_clinvar_reversible_path),
        sql_path(real_case_hashed_path),
        sql_path(real_clinvar_hashed_path)
      )
    )
  }
)

real_positional <- measure_thread_grid(
  "real_clinvar_positional_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) a",
          "  ON q.chrom = a.chrom AND q.pos = a.pos"
        ),
        sql_path(real_case_path),
        sql_path(real_clinvar_path)
      )
    )
  }
)

real_exact_results <- rbind(
  real_exact_multicol,
  real_exact_variantkey,
  real_exact_variantkey_lanes,
  real_positional
)
stopifnot(
  length(unique(rbind(
    real_exact_multicol,
    real_exact_variantkey,
    real_exact_variantkey_lanes
  )$result_rows)) == 1L
)
real_exact_results
#>                                          operation threads input_rows
#> 1             real_clinvar_exact_multicolumn_pairs       1    4096123
#> 2             real_clinvar_exact_multicolumn_pairs       2    4096123
#> 3             real_clinvar_exact_multicolumn_pairs       4    4096123
#> 4      real_clinvar_exact_variantkey_refined_pairs       1    4096123
#> 5      real_clinvar_exact_variantkey_refined_pairs       2    4096123
#> 6      real_clinvar_exact_variantkey_refined_pairs       4    4096123
#> 7  real_clinvar_exact_variantkey_split_lanes_pairs       1    4096123
#> 8  real_clinvar_exact_variantkey_split_lanes_pairs       2    4096123
#> 9  real_clinvar_exact_variantkey_split_lanes_pairs       4    4096123
#> 10                   real_clinvar_positional_pairs       1    4096123
#> 11                   real_clinvar_positional_pairs       2    4096123
#> 12                   real_clinvar_positional_pairs       4    4096123
#>    result_rows median_seconds input_rows_per_second result_rows_per_second
#> 1        44561          0.581               7050126               76697.07
#> 2        44561          0.310              13213300              143745.16
#> 3        44561          0.183              22383186              243502.73
#> 4        44561          0.606               6759279               73533.00
#> 5        44561          0.317              12921524              140570.98
#> 6        44561          0.209              19598675              213210.53
#> 7        44561          0.341              12012091              130677.42
#> 8        44561          0.196              20898587              227352.04
#> 9        44561          0.109              37579110              408816.51
#> 10       57488          0.381              10750979              150887.14
#> 11       57488          0.210              19505348              273752.38
#> 12       57488          0.128              32000961              449125.00
```

The split-lane layout processes 12.01 million GIAB alleles/s on one
thread and 37.58 million/s on four threads. It is faster than the
multicolumn oracle because only the 1.46% GIAB hashed lane and 0.88%
ClinVar hashed lane carry allele strings into the refinement join.

The same collision-safe split is applied to the GRCh38 ClinvArbitration
release. Its classification and star columns are provider payload;
DuckHTS does not recompute the ClinvArbitration policy in this
benchmark.

``` r
real_clinvarbitration_exact <- measure_thread_grid(
  "ClinvArbitration_16792026_exact_split_lanes",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste0(
          "SELECT sum(n)::BIGINT AS n FROM (",
          "SELECT count(*)::BIGINT AS n FROM read_parquet(%s) q ",
          "JOIN read_parquet(%s) a USING (vk) UNION ALL ",
          "SELECT count(*)::BIGINT AS n FROM read_parquet(%s) q ",
          "JOIN read_parquet(%s) a ON q.vk = a.vk AND q.chrom = a.chrom ",
          "AND q.pos = a.pos AND q.ref = a.ref AND q.alt = a.alt)"
        ),
        sql_path(real_case_reversible_path),
        sql_path(real_clinvarbitration_reversible_path),
        sql_path(real_case_hashed_path),
        sql_path(real_clinvarbitration_hashed_path)
      )
    )
  }
)
real_clinvarbitration_exact
#>                                     operation threads input_rows result_rows
#> 1 ClinvArbitration_16792026_exact_split_lanes       1    4096123       40318
#> 2 ClinvArbitration_16792026_exact_split_lanes       2    4096123       40318
#> 3 ClinvArbitration_16792026_exact_split_lanes       4    4096123       40318
#>   median_seconds input_rows_per_second result_rows_per_second
#> 1          0.296              13838253               136209.5
#> 2          0.151              27126642               267006.6
#> 3          0.084              48763369               479976.2
```

### Dense pathogenicity-score providers

REVEL and AlphaMissense are intentionally measured as separate exact
providers, not folded into the consequence kernel. The AlphaMissense run
uses the real GRCh38 GIAB relation. The REVEL v1.3 source is GRCh37, so
its throughput probe is a deterministic `hash(vk) % 19 = 0` sample
spanning the provider’s complete chromosome/key range. This prevents a
leading-key sample from making DuckDB dynamically prune most provider
row groups. It is a storage and execution measurement, not a biological
cross-build comparison.

``` r
alpha_probe_rows <- lane_staging_results$rows[
  lane_staging_results$source == "GIAB reversible lane"
]
revel_provider_rows <- staging_results$rows[
  staging_results$source == "REVEL v1.3 to VariantKey Parquet"
]
alpha_provider_rows <- staging_results$rows[
  staging_results$source == "AlphaMissense v2 to VariantKey Parquet"
]

real_revel_join <- measure_thread_grid(
  "REVEL_v1.3_exact_score_join",
  revel_probe_staging$rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n, sum(a.revel) AS score_sum",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) a USING (vk)"
        ),
        sql_path(real_revel_probe_path),
        sql_path(real_revel_path)
      )
    )
  }
)

real_alphamissense_join <- measure_thread_grid(
  "AlphaMissense_v2_exact_score_join",
  alpha_probe_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n, sum(a.am_pathogenicity) AS score_sum",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) a USING (vk)"
        ),
        sql_path(real_case_reversible_path),
        sql_path(real_alphamissense_path)
      )
    )
  }
)

dense_score_results <- rbind(real_revel_join, real_alphamissense_join)
dense_score_results
#>                           operation threads input_rows result_rows
#> 1       REVEL_v1.3_exact_score_join       1    4103497     4103497
#> 2       REVEL_v1.3_exact_score_join       2    4103497     4103497
#> 3       REVEL_v1.3_exact_score_join       4    4103497     4103497
#> 4 AlphaMissense_v2_exact_score_join       1    4036258        9225
#> 5 AlphaMissense_v2_exact_score_join       2    4036258        9225
#> 6 AlphaMissense_v2_exact_score_join       4    4036258        9225
#>   median_seconds input_rows_per_second result_rows_per_second
#> 1          3.239               1266902            1266902.439
#> 2          1.698               2416665            2416664.900
#> 3          0.924               4441014            4441014.069
#> 4          2.647               1524842               3485.077
#> 5          1.409               2864626               6547.197
#> 6          0.777               5194669              11872.587
```

The slower AlphaMissense path still processes 1.52 million GIAB
alleles/s on one thread and 5.19 million/s on four. The provider has
71,697,556 rows and its hot Parquet projection contains only VariantKey,
score, and class. REVEL’s hot projection has 77,966,138 rows.

### Row-group zone-map pruning for a local delta

Whole-genome joins necessarily visit the whole provider. Incremental
annotation is different: a new local batch can carry an explicit
VariantKey range, allowing Parquet row-group statistics to reject
unrelated provider pages. This real 5 Mb chromosome-17 slice compares
three row-group sizes while holding rows, ordering, compression, query,
and result cardinality fixed.

``` r
delta_chrom <- "17"
delta_pos_min <- 40000000L
delta_pos_max <- 45000000L
delta_rows <- as.numeric(DBI::dbGetQuery(
  con,
  sprintf(
    paste(
      "SELECT count(*) AS n FROM read_parquet(%s)",
      "WHERE chrom = %s AND pos BETWEEN %d AND %d"
    ),
    sql_path(real_case_reversible_path),
    sql_string(delta_chrom),
    delta_pos_min,
    delta_pos_max
  )
)$n[1L])

row_group_layouts$candidate_row_groups <- vapply(
  row_group_layouts$path,
  function(path) {
    as.numeric(DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "WITH bounds AS (SELECT variantkey_range(%s, %d, %d) AS b)",
          "SELECT count(DISTINCT row_group_id)::BIGINT AS n",
          "FROM parquet_metadata(%s), bounds",
          "WHERE path_in_schema = 'vk'",
          "  AND try_cast(stats_max_value AS UBIGINT) >= b.min",
          "  AND try_cast(stats_min_value AS UBIGINT) <= b.max"
        ),
        sql_string(delta_chrom),
        delta_pos_min,
        delta_pos_max,
        sql_path(path)
      )
    )$n[1L])
  },
  numeric(1L)
)
row_group_layouts
#>   row_group_size
#> 1         122880
#> 2          32768
#> 3        1048576
#>                                                                                       path
#> 1           /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible.parquet
#> 2   /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible_rg32768.parquet
#> 3 /root/duckhts/.tmp/variantkey_join_overlap/clinvar_20260706_reversible_rg1048576.parquet
#>      rows    bytes seconds row_groups candidate_row_groups
#> 1 4399443 27029899      NA         36                    1
#> 2 4399443 26945479   0.275        135                    2
#> 3 4399443 27004414   0.371          5                    1

row_group_query_results <- do.call(
  rbind,
  lapply(seq_len(nrow(row_group_layouts)), function(i) {
    provider_path <- row_group_layouts$path[i]
    measure_thread_grid(
      sprintf("clinvar_local_delta_rg%d", row_group_layouts$row_group_size[i]),
      delta_rows,
      function() {
        DBI::dbGetQuery(
          con,
          sprintf(
            paste(
              "SELECT count(*) AS n",
              "FROM read_parquet(%s) q",
              "JOIN (",
              "  SELECT vk FROM read_parquet(%s)",
              "  WHERE vk BETWEEN",
              "        (variantkey_range(%s, %d, %d)).min AND",
              "        (variantkey_range(%s, %d, %d)).max",
              ") a USING (vk)",
              "WHERE q.chrom = %s AND q.pos BETWEEN %d AND %d"
            ),
            sql_path(real_case_reversible_path),
            sql_path(provider_path),
            sql_string(delta_chrom),
            delta_pos_min,
            delta_pos_max,
            sql_string(delta_chrom),
            delta_pos_min,
            delta_pos_max,
            sql_string(delta_chrom),
            delta_pos_min,
            delta_pos_max
          )
        )
      },
      runs = max(bench_runs, 10L)
    )
  })
)
stopifnot(length(unique(row_group_query_results$result_rows)) == 1L)
row_group_query_results
#>                       operation threads input_rows result_rows median_seconds
#> 1  clinvar_local_delta_rg122880       1       6023         240         0.0085
#> 2  clinvar_local_delta_rg122880       2       6023         240         0.0090
#> 3  clinvar_local_delta_rg122880       4       6023         240         0.0090
#> 4   clinvar_local_delta_rg32768       1       6023         240         0.0090
#> 5   clinvar_local_delta_rg32768       2       6023         240         0.0090
#> 6   clinvar_local_delta_rg32768       4       6023         240         0.0090
#> 7 clinvar_local_delta_rg1048576       1       6023         240         0.0200
#> 8 clinvar_local_delta_rg1048576       2       6023         240         0.0200
#> 9 clinvar_local_delta_rg1048576       4       6023         240         0.0200
#>   input_rows_per_second result_rows_per_second
#> 1              708588.2               28235.29
#> 2              669222.2               26666.67
#> 3              669222.2               26666.67
#> 4              669222.2               26666.67
#> 5              669222.2               26666.67
#> 6              669222.2               26666.67
#> 7              301150.0               12000.00
#> 8              301150.0               12000.00
#> 9              301150.0               12000.00
```

The default 122,880-row layout admits 1 of 36 provider row groups and
answers the 6,023-allele local delta in 0.0085 seconds on one thread.
The 1,048,576-row layout also admits one row group, but that group
contains much more unrelated data and is slower. Smaller row groups do
not help this 5 Mb query enough to offset their metadata and scheduling
overhead.

### Regulatory interval joins

Ensembl coordinates are converted from one-based inclusive to zero-based
half-open before RegionKey construction. The maximum source span is read
from the staged relation and used as the proven RegionKey look-back
bound. A per-contig vanilla SQL query exposes typed start/end
inequalities directly to DuckDB’s interval join, while the
packed-coordinate form exposes the same problem in one numeric domain.

``` r
real_regulatory_stats <- DBI::dbGetQuery(
  con,
  sprintf(
    paste(
      "SELECT count(*) AS rows, max(end0 - start0)::BIGINT AS max_span,",
      "       sum(end0 - start0)::UBIGINT AS summed_span",
      "FROM read_parquet(%s)"
    ),
    sql_path(real_regulatory_path)
  )
)
real_regulatory_stats
#>     rows max_span summed_span
#> 1 643528    48907   316683589

real_idx_name <- "variantkey_real_regulatory_idx"
try(
  DBI::dbGetQuery(con, sprintf("SELECT duckhts_cgranges_destroy('%s') AS ok", real_idx_name)),
  silent = TRUE
)
real_index_query <- sprintf(
  "SELECT chrom, start0, end0, stable_id FROM read_parquet(%s)",
  sql_path(real_regulatory_path)
)
real_index_seconds <- system.time({
  DBI::dbGetQuery(
    con,
    paste0(
      "SELECT duckhts_cgranges_from_query(",
      sql_string(real_idx_name), ", ",
      sql_string(real_index_query), ", ",
      sql_string("chrom"), ", ", sql_string("start0"), ", ",
      sql_string("end0"), ", ", sql_string("stable_id"),
      ") AS ok"
    )
  )
  DBI::dbGetQuery(con, sprintf("SELECT duckhts_cgranges_index('%s') AS ok", real_idx_name))
})[["elapsed"]]

data.frame(
  operation = "cgranges_build_and_finalize",
  intervals = real_regulatory_stats$rows,
  seconds = unname(real_index_seconds),
  intervals_per_second = real_regulatory_stats$rows / real_index_seconds
)
#>                     operation intervals seconds intervals_per_second
#> 1 cgranges_build_and_finalize    643528   0.098              6566612
```

``` r
real_contigs <- DBI::dbGetQuery(
  con,
  sprintf("SELECT DISTINCT chrom FROM read_parquet(%s) ORDER BY chrom", sql_path(real_regulatory_path))
)$chrom

real_typed_contig_queries <- vapply(
  real_contigs,
  function(chrom) {
    sprintf(
      paste(
        "SELECT count(*)::BIGINT AS n",
        "FROM read_parquet(%s) q",
        "JOIN read_parquet(%s) s",
        "  ON q.start0 < s.end0 AND q.end0 > s.start0",
        "WHERE q.chrom = %s AND s.chrom = %s"
      ),
      sql_path(real_case_path),
      sql_path(real_regulatory_path),
      sql_string(chrom),
      sql_string(chrom)
    )
  },
  character(1L)
)
real_typed_interval_sql <- paste0(
  "SELECT sum(n)::BIGINT AS n FROM (",
  paste(real_typed_contig_queries, collapse = " UNION ALL "),
  ")"
)

real_typed_interval_plan <- DBI::dbGetQuery(
  con,
  paste("EXPLAIN", real_typed_contig_queries[[1L]])
)$explain_value
real_typed_iejoin_count <- sum(lengths(regmatches(
  real_typed_interval_plan,
  gregexpr("IE_JOIN", real_typed_interval_plan, fixed = TRUE)
)))
stopifnot(real_typed_iejoin_count == 1L)

data.frame(
  physical_operator = "IE_JOIN",
  representative_contig = real_contigs[[1L]],
  operators_in_representative_plan = real_typed_iejoin_count,
  contig_shards_in_full_query = length(real_contigs)
)
#>   physical_operator representative_contig operators_in_representative_plan
#> 1           IE_JOIN                     1                                1
#>   contig_shards_in_full_query
#> 1                          24

real_typed_interval <- measure_thread_grid(
  "real_regulatory_typed_per_contig_iejoin_pairs",
  real_case_rows,
  function() DBI::dbGetQuery(con, real_typed_interval_sql)
)

real_packed_interval <- measure_thread_grid(
  "real_regulatory_packed_coordinate_iejoin_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) s",
          "  ON q.rk_chrom_start < s.rk_chrom_end",
          " AND q.rk_chrom_end > s.rk_chrom_start"
        ),
        sql_path(real_case_path),
        sql_path(real_regulatory_path)
      )
    )
  }
)

real_max_span <- as.numeric(real_regulatory_stats$max_span[1L])
real_regionkey_interval <- measure_thread_grid(
  "real_regulatory_regionkey_bounded_start_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM read_parquet(%s) q",
          "JOIN read_parquet(%s) s",
          "  ON s.rk_chrom_start >=",
          "     (regionkey(q.chrom, greatest(0, q.start0 - %d + 1), greatest(0, q.start0 - %d + 1)) >> 31)",
          " AND s.rk_chrom_start < (regionkey(q.chrom, q.end0, q.end0) >> 31)",
          " AND are_overlapping_regionkeys(q.rk, s.rk)"
        ),
        sql_path(real_case_path),
        sql_path(real_regulatory_path),
        real_max_span,
        real_max_span
      )
    )
  }
)

real_cgranges_interval <- measure_thread_grid(
  "real_regulatory_cgranges_count_pairs",
  real_case_rows,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT sum(duckhts_cgranges_count_overlaps('%s', chrom, start0, end0))::BIGINT AS n",
          "FROM read_parquet(%s)"
        ),
        real_idx_name,
        sql_path(real_case_path)
      )
    )
  }
)

real_interval_results <- rbind(
  real_typed_interval,
  real_packed_interval,
  real_regionkey_interval,
  real_cgranges_interval
)
stopifnot(length(unique(real_interval_results$result_rows)) == 1L)
real_interval_results
#>                                         operation threads input_rows
#> 1   real_regulatory_typed_per_contig_iejoin_pairs       1    4096123
#> 2   real_regulatory_typed_per_contig_iejoin_pairs       2    4096123
#> 3   real_regulatory_typed_per_contig_iejoin_pairs       4    4096123
#> 4  real_regulatory_packed_coordinate_iejoin_pairs       1    4096123
#> 5  real_regulatory_packed_coordinate_iejoin_pairs       2    4096123
#> 6  real_regulatory_packed_coordinate_iejoin_pairs       4    4096123
#> 7   real_regulatory_regionkey_bounded_start_pairs       1    4096123
#> 8   real_regulatory_regionkey_bounded_start_pairs       2    4096123
#> 9   real_regulatory_regionkey_bounded_start_pairs       4    4096123
#> 10           real_regulatory_cgranges_count_pairs       1    4096123
#> 11           real_regulatory_cgranges_count_pairs       2    4096123
#> 12           real_regulatory_cgranges_count_pairs       4    4096123
#>    result_rows median_seconds input_rows_per_second result_rows_per_second
#> 1       407430          1.178               3477184               345865.9
#> 2       407430          0.750               5461497               543240.0
#> 3       407430          0.383              10694838              1063785.9
#> 4       407430          1.117               3667075               364753.8
#> 5       407430          0.643               6370331               633639.2
#> 6       407430          0.411               9966236               991313.9
#> 7       407430          3.580               1144168               113807.3
#> 8       407430          1.939               2112493               210123.8
#> 9       407430          1.078               3799743               377949.9
#> 10      407430          0.370              11070603              1101162.2
#> 11      407430          0.191              21445670              2133141.4
#> 12      407430          0.104              39385798              3917596.2

real_interval_speedups <- merge(
  real_interval_results[, c("operation", "threads", "input_rows_per_second")],
  real_typed_interval[, c("threads", "input_rows_per_second")],
  by = "threads",
  suffixes = c("", "_typed_baseline")
)
real_interval_speedups$speedup_vs_typed_iejoin <-
  real_interval_speedups$input_rows_per_second /
  real_interval_speedups$input_rows_per_second_typed_baseline
real_interval_speedups
#>    threads                                      operation input_rows_per_second
#> 1        1  real_regulatory_typed_per_contig_iejoin_pairs               3477184
#> 2        1           real_regulatory_cgranges_count_pairs              11070603
#> 3        1  real_regulatory_regionkey_bounded_start_pairs               1144168
#> 4        1 real_regulatory_packed_coordinate_iejoin_pairs               3667075
#> 5        2 real_regulatory_packed_coordinate_iejoin_pairs               6370331
#> 6        2  real_regulatory_typed_per_contig_iejoin_pairs               5461497
#> 7        2           real_regulatory_cgranges_count_pairs              21445670
#> 8        2  real_regulatory_regionkey_bounded_start_pairs               2112493
#> 9        4  real_regulatory_regionkey_bounded_start_pairs               3799743
#> 10       4 real_regulatory_packed_coordinate_iejoin_pairs               9966236
#> 11       4  real_regulatory_typed_per_contig_iejoin_pairs              10694838
#> 12       4           real_regulatory_cgranges_count_pairs              39385798
#>    input_rows_per_second_typed_baseline speedup_vs_typed_iejoin
#> 1                               3477184               1.0000000
#> 2                               3477184               3.1837838
#> 3                               3477184               0.3290503
#> 4                               3477184               1.0546106
#> 5                               5461497               1.1664075
#> 6                               5461497               1.0000000
#> 7                               5461497               3.9267016
#> 8                               5461497               0.3867973
#> 9                              10694838               0.3552876
#> 10                             10694838               0.9318735
#> 11                             10694838               1.0000000
#> 12                             10694838               3.6826923
```

The physical-plan assertion above matters. DuckDB’s [IEJoin
description](https://duckdb.org/2022/05/27/iejoin) explains that a
double-inequality join sorts both relations on the first condition,
constructs merged and permutation arrays for the second condition, and
uses a bitmap to identify matching pairs. DuckDB parallelizes sorted
blocks and can spool join blocks to disk. The measured 665/728 MiB is
therefore the resident cost of real parallel `IE_JOIN`
planning/execution over these Parquet relations, not a nested loop
mislabeled as an interval join. The per-contig `UNION ALL` is
deliberate: it makes contig equality a static shard filter and leaves
exactly the two coordinate inequalities for each `IE_JOIN` operator.

For this real source, cgranges is 3.18 times the typed DuckDB interval
join on one thread and 3.68 times on four threads. The globally bounded
RegionKey range is correct but loses because the single 48,907-base
source interval widens candidate discovery for every GIAB allele.

### Gene annotations at gene cardinality

Gene annotations are not a multi-million-row allele problem. They should
be joined after affected transcript rows have been reduced to distinct
genes. This comparison uses the same full Ensembl model both ways to
expose the avoidable transcript multiplicity.

``` r
gene_transcript_join <- measure_thread_grid(
  "gnomad_constraint_join_per_transcript",
  644292,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM read_parquet(%s) t",
          "JOIN read_parquet(%s) g USING (gene_id)"
        ),
        sql_path(transcripts_parquet),
        sql_path(real_constraint_path)
      )
    )
  }
)

gene_distinct_join <- measure_thread_grid(
  "gnomad_constraint_join_after_gene_dedup",
  644292,
  function() {
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT count(*) AS n",
          "FROM (SELECT DISTINCT gene_id FROM read_parquet(%s)) t",
          "JOIN read_parquet(%s) g USING (gene_id)"
        ),
        sql_path(transcripts_parquet),
        sql_path(real_constraint_path)
      )
    )
  }
)

rbind(gene_transcript_join, gene_distinct_join)
#>                                 operation threads input_rows result_rows
#> 1   gnomad_constraint_join_per_transcript       1     644292      426382
#> 2   gnomad_constraint_join_per_transcript       2     644292      426382
#> 3   gnomad_constraint_join_per_transcript       4     644292      426382
#> 4 gnomad_constraint_join_after_gene_dedup       1     644292       19082
#> 5 gnomad_constraint_join_after_gene_dedup       2     644292       19082
#> 6 gnomad_constraint_join_after_gene_dedup       4     644292       19082
#>   median_seconds input_rows_per_second result_rows_per_second
#> 1          0.016              40268250               26648875
#> 2          0.010              64429200               42638200
#> 3          0.008              80536500               53297750
#> 4          0.012              53691000                1590167
#> 5          0.008              80536500                2385250
#> 6          0.006             107382000                3180333
```

### Materialized provider output

Counting is useful for operator isolation, but the clinical pipeline
must keep the matching annotation rows. These four-thread runs write the
exact ClinVar and regulatory results to ZSTD Parquet and report both row
and byte output.

``` r
set_duck_threads(4L)
#> [1] 0
real_exact_output_path <- file.path(real_cache_dir, "giab_clinvar_exact_matches.parquet")
real_regulatory_output_path <- file.path(real_cache_dir, "giab_regulatory_overlaps.parquet")

real_exact_output_query <- sprintf(
  paste0(
    "SELECT q.chrom, q.pos, q.ref, q.alt, ",
    "a.allele_id, a.clinical_significance, a.gene_info ",
    "FROM read_parquet(%s) q JOIN read_parquet(%s) a USING (vk) ",
    "UNION ALL ",
    "SELECT q.chrom, q.pos, q.ref, q.alt, ",
    "a.allele_id, a.clinical_significance, a.gene_info ",
    "FROM read_parquet(%s) q JOIN read_parquet(%s) a ",
    "ON q.vk = a.vk AND q.chrom = a.chrom AND q.pos = a.pos ",
    "AND q.ref = a.ref AND q.alt = a.alt"
  ),
  sql_path(real_case_reversible_path),
  sql_path(real_clinvar_reversible_path),
  sql_path(real_case_hashed_path),
  sql_path(real_clinvar_hashed_path)
)

real_regulatory_output_query <- sprintf(
  paste(
    "SELECT q.chrom, q.pos, q.ref, q.alt,",
    "       s.stable_id, s.feature_type, s.so_term, s.so_accession",
    "FROM read_parquet(%s) q",
    "JOIN read_parquet(%s) s",
    "  ON q.rk_chrom_start < s.rk_chrom_end",
    " AND q.rk_chrom_end > s.rk_chrom_start"
  ),
  sql_path(real_case_path),
  sql_path(real_regulatory_path)
)

materialized_results <- rbind(
  cbind(output = "GIAB exact ClinVar annotations", copy_timed(real_exact_output_query, real_exact_output_path)),
  cbind(output = "GIAB Ensembl regulatory overlaps", copy_timed(real_regulatory_output_query, real_regulatory_output_path))
)
materialized_results$rows_per_second <- materialized_results$rows / materialized_results$seconds
materialized_results$input_rows <- real_case_rows
materialized_results$input_rows_per_second <- real_case_rows / materialized_results$seconds
materialized_results$output_mb_per_second <- materialized_results$bytes / materialized_results$seconds / 1e6
stopifnot(materialized_results$rows[1L] == real_exact_variantkey_lanes$result_rows[1L])
stopifnot(materialized_results$rows[2L] == real_packed_interval$result_rows[1L])
materialized_results
#>                             output
#> 1   GIAB exact ClinVar annotations
#> 2 GIAB Ensembl regulatory overlaps
#>                                                                            path
#> 1 /root/duckhts/.tmp/variantkey_join_overlap/giab_clinvar_exact_matches.parquet
#> 2   /root/duckhts/.tmp/variantkey_join_overlap/giab_regulatory_overlaps.parquet
#>     rows   bytes seconds rows_per_second input_rows input_rows_per_second
#> 1  44561  440154   0.238        187231.1    4096123              17210601
#> 2 407430 2230407   0.551        739437.4    4096123               7433980
#>   output_mb_per_second
#> 1             1.849387
#> 2             4.047926
```

## 6. Isolated peak memory and population-scale projection

The preceding timings run in one long-lived R/DuckDB process. Peak
memory is measured separately because cumulative allocator state would
make an in-process RSS sample meaningless. Each row below is a fresh
DuckDB CLI process wrapped by `/usr/bin/time -v`, with the same CPU
pinning and an in-memory database. Reported RSS includes the DuckDB
process, decompression buffers, operator state, cgranges construction
transients, and the loaded extension; it excludes the operating system
page cache.

``` r
rss_cpu_1 <- Sys.getenv("VARIANTKEY_RSS_CPU_1", unset = "0")
rss_cpu_4 <- Sys.getenv("VARIANTKEY_RSS_CPU_4", unset = "0,2,4,6")

exact_lanes_sql <- sprintf(
  paste0(
    "SELECT sum(n)::BIGINT AS n FROM (",
    "SELECT count(*)::BIGINT AS n FROM read_parquet(%s) q ",
    "JOIN read_parquet(%s) a USING (vk) UNION ALL ",
    "SELECT count(*)::BIGINT AS n FROM read_parquet(%s) q ",
    "JOIN read_parquet(%s) a ON q.vk = a.vk AND q.chrom = a.chrom ",
    "AND q.pos = a.pos AND q.ref = a.ref AND q.alt = a.alt)"
  ),
  sql_path(real_case_reversible_path),
  sql_path(real_clinvar_reversible_path),
  sql_path(real_case_hashed_path),
  sql_path(real_clinvar_hashed_path)
)
clinvarbitration_lanes_sql <- sprintf(
  paste0(
    "SELECT sum(n)::BIGINT AS n FROM (",
    "SELECT count(*)::BIGINT AS n FROM read_parquet(%s) q ",
    "JOIN read_parquet(%s) a USING (vk) UNION ALL ",
    "SELECT count(*)::BIGINT AS n FROM read_parquet(%s) q ",
    "JOIN read_parquet(%s) a ON q.vk = a.vk AND q.chrom = a.chrom ",
    "AND q.pos = a.pos AND q.ref = a.ref AND q.alt = a.alt)"
  ),
  sql_path(real_case_reversible_path),
  sql_path(real_clinvarbitration_reversible_path),
  sql_path(real_case_hashed_path),
  sql_path(real_clinvarbitration_hashed_path)
)

revel_join_sql <- sprintf(
  paste(
    "SELECT count(*) AS n, sum(a.revel) AS score_sum",
    "FROM read_parquet(%s) q JOIN read_parquet(%s) a USING (vk)"
  ),
  sql_path(real_revel_probe_path),
  sql_path(real_revel_path)
)
alpha_join_sql <- sprintf(
  paste(
    "SELECT count(*) AS n, sum(a.am_pathogenicity) AS score_sum",
    "FROM read_parquet(%s) q JOIN read_parquet(%s) a USING (vk)"
  ),
  sql_path(real_case_reversible_path),
  sql_path(real_alphamissense_path)
)

alpha_chr1_rows <- as.numeric(DBI::dbGetQuery(
  con,
  sprintf(
    "SELECT count(*) AS n FROM read_parquet(%s) WHERE extract_variantkey_chrom(vk) = 1",
    sql_path(real_alphamissense_path)
  )
)$n[1L])
alpha_chr1_join_sql <- sprintf(
  paste(
    "SELECT count(*) AS n, sum(a.am_pathogenicity) AS score_sum",
    "FROM (SELECT vk FROM read_parquet(%s)",
    "      WHERE extract_variantkey_chrom(vk) = 1) q",
    "JOIN (SELECT vk, am_pathogenicity FROM read_parquet(%s)",
    "      WHERE extract_variantkey_chrom(vk) = 1) a USING (vk)"
  ),
  sql_path(real_alphamissense_path),
  sql_path(real_alphamissense_path)
)

real_cgranges_sql <- paste0(
  "SELECT duckhts_cgranges_from_query('rss_reg', ",
  sql_string(real_index_query), ", 'chrom', 'start0', 'end0', 'stable_id'); ",
  "SELECT duckhts_cgranges_index('rss_reg'); ",
  sprintf(
    paste(
      "SELECT sum(duckhts_cgranges_count_overlaps('rss_reg', chrom, start0, end0))",
      "FROM read_parquet(%s)"
    ),
    sql_path(real_case_path)
  )
)

cgranges_synthetic_sql <- function(n, label_kind = "bigint") {
  label_sql <- switch(
    label_kind,
    none = "",
    bigint = ", i::BIGINT AS label",
    varchar = ", 'feature_' || i::VARCHAR AS label",
    stop("unknown cgranges label kind", call. = FALSE)
  )
  label_arg <- if (label_kind == "none") "" else ", 'label'"
  query <- sprintf(
    "SELECT '1' AS chrom, i * 10 AS start0, i * 10 + 5 AS end0%s FROM range(%d) t(i)",
    label_sql,
    as.numeric(n)
  )
  paste0(
    "SELECT duckhts_cgranges_from_query('rss_syn', ", sql_string(query),
    ", 'chrom', 'start0', 'end0'", label_arg, "); ",
    "SELECT duckhts_cgranges_index('rss_syn')"
  )
}

rss_alpha_stage_path <- file.path(real_cache_dir, "rss_alphamissense_stage.parquet")
rss_revel_stage_path <- file.path(real_cache_dir, "rss_revel_stage.parquet")
alpha_stage_sql <- sprintf(
  "COPY (%s) TO %s (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)",
  alphamissense_stage_query,
  sql_path(rss_alpha_stage_path)
)
revel_stage_sql <- sprintf(
  "COPY (%s) TO %s (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)",
  revel_stage_query,
  sql_path(rss_revel_stage_path)
)

isolated_memory <- rbind(
  isolated_duckdb_measure("baseline", "SELECT 1", 1L, rss_cpu_1),
  isolated_duckdb_measure("baseline", "SELECT 1", 4L, rss_cpu_4),
  isolated_duckdb_measure("ClinVar split exact join", exact_lanes_sql, 1L, rss_cpu_1),
  isolated_duckdb_measure("ClinVar split exact join", exact_lanes_sql, 4L, rss_cpu_4),
  isolated_duckdb_measure(
    "ClinvArbitration split exact join",
    clinvarbitration_lanes_sql,
    1L,
    rss_cpu_1
  ),
  isolated_duckdb_measure(
    "ClinvArbitration split exact join",
    clinvarbitration_lanes_sql,
    4L,
    rss_cpu_4
  ),
  isolated_duckdb_measure("REVEL 77.97M exact join", revel_join_sql, 1L, rss_cpu_1),
  isolated_duckdb_measure("REVEL 77.97M exact join", revel_join_sql, 4L, rss_cpu_4),
  isolated_duckdb_measure("AlphaMissense 71.70M exact join", alpha_join_sql, 1L, rss_cpu_1),
  isolated_duckdb_measure("AlphaMissense 71.70M exact join", alpha_join_sql, 4L, rss_cpu_4),
  isolated_duckdb_measure("AlphaMissense chromosome 1 shard self-join", alpha_chr1_join_sql, 1L, rss_cpu_1),
  isolated_duckdb_measure("AlphaMissense chromosome 1 shard self-join", alpha_chr1_join_sql, 4L, rss_cpu_4),
  isolated_duckdb_measure("regulatory cgranges build and full GIAB probe", real_cgranges_sql, 1L, rss_cpu_1),
  isolated_duckdb_measure("regulatory DuckDB per-contig IEJoin", real_typed_interval_sql, 1L, rss_cpu_1),
  isolated_duckdb_measure("regulatory DuckDB per-contig IEJoin", real_typed_interval_sql, 4L, rss_cpu_4),
  isolated_duckdb_measure("cgranges 1M no labels", cgranges_synthetic_sql(1000000, "none"), 1L, rss_cpu_1),
  isolated_duckdb_measure("cgranges 1M BIGINT labels", cgranges_synthetic_sql(1000000, "bigint"), 1L, rss_cpu_1),
  isolated_duckdb_measure("cgranges 1M VARCHAR labels", cgranges_synthetic_sql(1000000, "varchar"), 1L, rss_cpu_1),
  isolated_duckdb_measure("cgranges 10M BIGINT labels", cgranges_synthetic_sql(10000000, "bigint"), 1L, rss_cpu_1),
  isolated_duckdb_measure(
    "AlphaMissense provider preparation",
    alpha_stage_sql,
    4L,
    rss_cpu_4,
    cleanup_paths = rss_alpha_stage_path
  ),
  isolated_duckdb_measure(
    "REVEL provider preparation",
    revel_stage_sql,
    4L,
    rss_cpu_4,
    cleanup_paths = rss_revel_stage_path
  )
)
unlink(c(rss_alpha_stage_path, rss_revel_stage_path))

baseline_rss <- isolated_memory[isolated_memory$operation == "baseline", c("threads", "peak_rss_mib")]
names(baseline_rss)[2L] <- "baseline_rss_mib"
isolated_memory <- merge(isolated_memory, baseline_rss, by = "threads", all.x = TRUE, sort = FALSE)
isolated_memory$incremental_peak_rss_mib <- pmax(
  0,
  isolated_memory$peak_rss_mib - isolated_memory$baseline_rss_mib
)
isolated_memory
#>    threads                                     operation repeats peak_rss_mib
#> 1        1                                      baseline       1     36.04297
#> 2        1           regulatory DuckDB per-contig IEJoin       1    664.76562
#> 3        1                      ClinVar split exact join       1    185.25000
#> 4        1                         cgranges 1M no labels       1    129.20703
#> 5        1             ClinvArbitration split exact join       1    175.42578
#> 6        1                    cgranges 1M VARCHAR labels       1    196.70312
#> 7        1                       REVEL 77.97M exact join       1    183.21484
#> 8        1               AlphaMissense 71.70M exact join       1    181.33203
#> 9        1    AlphaMissense chromosome 1 shard self-join       1    307.30469
#> 10       1 regulatory cgranges build and full GIAB probe       1    145.64062
#> 11       1                    cgranges 10M BIGINT labels       1   1281.66016
#> 12       1                     cgranges 1M BIGINT labels       1    144.98438
#> 13       4                      ClinVar split exact join       1    191.89062
#> 14       4                                      baseline       1     35.78125
#> 15       4                       REVEL 77.97M exact join       1    189.70312
#> 16       4             ClinvArbitration split exact join       1    179.90625
#> 17       4    AlphaMissense chromosome 1 shard self-join       1    315.09375
#> 18       4               AlphaMissense 71.70M exact join       1    188.63672
#> 19       4           regulatory DuckDB per-contig IEJoin       1    711.98438
#> 20       4            AlphaMissense provider preparation       1   1128.80078
#> 21       4                    REVEL provider preparation       1   1278.75391
#>    median_peak_rss_mib median_elapsed_seconds baseline_rss_mib
#> 1             36.04297                   0.13         36.04297
#> 2            664.76562                   1.14         36.04297
#> 3            185.25000                   0.46         36.04297
#> 4            129.20703                   0.18         36.04297
#> 5            175.42578                   0.37         36.04297
#> 6            196.70312                   0.26         36.04297
#> 7            183.21484                   2.87         36.04297
#> 8            181.33203                   2.43         36.04297
#> 9            307.30469                   1.49         36.04297
#> 10           145.64062                   0.60         36.04297
#> 11          1281.66016                   1.11         36.04297
#> 12           144.98438                   0.22         36.04297
#> 13           191.89062                   0.13         35.78125
#> 14            35.78125                   0.02         35.78125
#> 15           189.70312                   0.85         35.78125
#> 16           179.90625                   0.11         35.78125
#> 17           315.09375                   0.44         35.78125
#> 18           188.63672                   0.72         35.78125
#> 19           711.98438                   0.36         35.78125
#> 20          1128.80078                  12.77         35.78125
#> 21          1278.75391                   2.31         35.78125
#>    incremental_peak_rss_mib
#> 1                   0.00000
#> 2                 628.72266
#> 3                 149.20703
#> 4                  93.16406
#> 5                 139.38281
#> 6                 160.66016
#> 7                 147.17188
#> 8                 145.28906
#> 9                 271.26172
#> 10                109.59766
#> 11               1245.61719
#> 12                108.94141
#> 13                156.10938
#> 14                  0.00000
#> 15                153.92188
#> 16                144.12500
#> 17                279.31250
#> 18                152.85547
#> 19                676.20312
#> 20               1093.01953
#> 21               1242.97266
```

The real 643,528-row regulatory cgranges index uses string labels
because an expanding annotation needs a stable ordinal or identifier.
Its peak includes index construction and the complete GIAB count probe.
The synthetic rows expose the current API’s payload cost: besides
cgranges’ 16-byte core interval, DuckHTS retains return coordinates,
label validity, and currently one owned chromosome string per interval.
A count-only compact constructor that interns chromosome names could
reduce this, but it is an optimization—not permission to put dbSNP or
TOPMed into an interval index.

``` r
cgranges_memory_density <- data.frame(
  operation = c(
    "cgranges 1M no labels",
    "cgranges 1M BIGINT labels",
    "cgranges 1M VARCHAR labels",
    "cgranges 10M BIGINT labels"
  ),
  intervals = c(1e6, 1e6, 1e6, 1e7),
  stringsAsFactors = FALSE
)
cgranges_memory_density <- merge(
  cgranges_memory_density,
  isolated_memory[, c("operation", "peak_rss_mib", "incremental_peak_rss_mib")],
  by = "operation",
  sort = FALSE
)
cgranges_memory_density$incremental_peak_bytes_per_interval <-
  cgranges_memory_density$incremental_peak_rss_mib * 1024^2 / cgranges_memory_density$intervals
cgranges_memory_density
#>                    operation intervals peak_rss_mib incremental_peak_rss_mib
#> 1      cgranges 1M no labels     1e+06     129.2070                 93.16406
#> 2  cgranges 1M BIGINT labels     1e+06     144.9844                108.94141
#> 3 cgranges 1M VARCHAR labels     1e+06     196.7031                160.66016
#> 4 cgranges 10M BIGINT labels     1e+07    1281.6602               1245.61719
#>   incremental_peak_bytes_per_interval
#> 1                             97.6896
#> 2                            114.2333
#> 3                            168.4644
#> 4                            130.6124
```

### Storage and logical time at catalog scale

The provider hot files below intentionally exclude source presentation
columns. They are serving projections, not replacements for the source
receipts. Exact providers are scanned one at a time; therefore twelve
independent annotations have approximately the sum of their elapsed
times but the maximum of their individual RSS values. Loading twelve
indexes concurrently is neither required nor the design.

``` r
provider_storage <- data.frame(
  provider = c(
    "ClinVar reversible exact lane",
    "ClinvArbitration exact lanes",
    "REVEL v1.3",
    "AlphaMissense v2"
  ),
  rows = c(
    lane_staging_results$rows[lane_staging_results$source == "ClinVar reversible numeric lane"],
    sum(lane_staging_results$rows[
      lane_staging_results$source %in% c(
        "ClinvArbitration reversible numeric lane",
        "ClinvArbitration hashed refinement lane"
      )
    ]),
    revel_provider_rows,
    alpha_provider_rows
  ),
  bytes = c(
    lane_staging_results$bytes[lane_staging_results$source == "ClinVar reversible numeric lane"],
    sum(lane_staging_results$bytes[
      lane_staging_results$source %in% c(
        "ClinvArbitration reversible numeric lane",
        "ClinvArbitration hashed refinement lane"
      )
    ]),
    staging_results$bytes[staging_results$source == "REVEL v1.3 to VariantKey Parquet"],
    staging_results$bytes[staging_results$source == "AlphaMissense v2 to VariantKey Parquet"]
  ),
  stringsAsFactors = FALSE
)
provider_storage$compressed_bytes_per_row <- provider_storage$bytes / provider_storage$rows
provider_storage
#>                        provider     rows     bytes compressed_bytes_per_row
#> 1 ClinVar reversible exact lane  4399443  27029899                 6.143937
#> 2  ClinvArbitration exact lanes  3647840  21814288                 5.980056
#> 3                    REVEL v1.3 77966138 201658023                 2.586482
#> 4              AlphaMissense v2 71697556 254370824                 3.547831

catalog_scenarios <- data.frame(
  catalog = c("TOPMed Freeze 8 BRAVO", "dbSNP Build 157 live RefSNP"),
  rows = c(705486649, 1172689405),
  stringsAsFactors = FALSE
)
catalog_scenarios$narrow_hot_gb_low <- catalog_scenarios$rows *
  min(provider_storage$compressed_bytes_per_row) / 1e9
catalog_scenarios$narrow_hot_gb_high <- catalog_scenarios$rows *
  max(provider_storage$compressed_bytes_per_row) / 1e9

dense_rates <- dense_score_results[dense_score_results$threads %in% c(1L, 4L), ]
logical_time <- do.call(
  rbind,
  lapply(seq_len(nrow(catalog_scenarios)), function(i) {
    do.call(
      rbind,
      lapply(c(1L, 4L), function(th) {
        rates <- dense_rates$input_rows_per_second[dense_rates$threads == th]
        data.frame(
          catalog = catalog_scenarios$catalog[i],
          threads = th,
          catalog_rows = catalog_scenarios$rows[i],
          twelve_dense_providers_minutes_fast = catalog_scenarios$rows[i] / max(rates) * 12 / 60,
          twelve_dense_providers_minutes_slow = catalog_scenarios$rows[i] / min(rates) * 12 / 60,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
catalog_scenarios
#>                       catalog       rows narrow_hot_gb_low narrow_hot_gb_high
#> 1       TOPMed Freeze 8 BRAVO  705486649          1.824729           4.334465
#> 2 dbSNP Build 157 live RefSNP 1172689405          3.033140           7.204929
logical_time
#>                       catalog threads catalog_rows
#> 1       TOPMed Freeze 8 BRAVO       1    705486649
#> 2       TOPMed Freeze 8 BRAVO       4    705486649
#> 3 dbSNP Build 157 live RefSNP       1   1172689405
#> 4 dbSNP Build 157 live RefSNP       4   1172689405
#>   twelve_dense_providers_minutes_fast twelve_dense_providers_minutes_slow
#> 1                            92.53240                           111.37190
#> 2                            27.16195                            31.77142
#> 3                           153.81122                           185.12703
#> 4                            45.14972                            52.81179
```

These are deliberately labelled logical-time projections: they linearly
apply the measured REVEL-to-AlphaMissense throughput envelope to twelve
equally dense providers. They do not claim that all real sources are
that large, equally compressible, or equally expensive. For a
705-million-row or 1.17-billion-row query, the execution macro must
compile chromosome or coordinate-tile tasks so one provider shard is the
hash-build side and the population relation streams. The measured
largest AlphaMissense chromosome has 7,183,937 source rows; its isolated
self-join is the concrete memory denominator for that scheduling model.
Four chromosome workers may share immutable Parquet pages, but each owns
its join state, so concurrency must be bounded by measured RSS rather
than by the number of available source files.

## Notes

- The exact-join sections are the relevant comparison class for future
  `echtvar`-style annotation throughput work.
- The `mixed_stored_vk_join_plus_hash_refine` query represents the safe
  SQL pattern for hashed/nonreversible VariantKey rows when exact
  `REF/ALT` verification is required.
- VariantKey does not normalize alleles. Exact provider joins require
  the same normalized biallelic representation on both sides; the real
  ClinVar result counts exact stored representations, not merely
  equivalent biological edits.
- The reversible exact lane joins only `UBIGINT vk`; chromosome and
  allele strings are absent from that provider file. A general non-human
  store should use a release-scoped numeric contig ordinal and retain
  the ordinal-to-name map in the dataset manifest, Parquet key/value
  metadata, or a companion relation. Official VariantKey/RegionKey
  chromosome codes remain limited to canonical human contigs.
- RegionKey is more than exact equality: its chromosome/start ordering
  supports binary search, ordered merge, and bounded start-range
  discovery. Arbitrary nested intervals still require a proven look-back
  bound, an augmented max-end index, or a complete overlap index such as
  cgranges.
- DuckDB’s native interval join is the zero-registry SQL baseline. The
  RegionKey range query and cgranges are alternatives with different
  storage, preparation, and repeated-query trade-offs; the measured
  report, not the key format alone, decides which path is preferable.
- The Ensembl regulatory source has one 48,907-base interval, so the
  globally bounded RegionKey query admits a wide candidate window for
  every probe. The result is correct but slower than both cgranges and
  DuckDB’s interval join; span-class partitioning or an augmented
  max-end index would be required to improve that RegionKey search
  without dropping long features.
- Sorted VariantKey row groups already provide useful min/max pruning.
  Physical chromosome/tile partitions become worthwhile when a
  distributed delta can be compiled into static touched-prefix scans; a
  plain hash join over sparse keys does not by itself guarantee file or
  row-group pruning.
- For richer overlap comparisons against external tools such as `bedtk`
  or `bedtools`, see the existing `benchmarks/benchmark_cgranges.Rmd`
  report.
- The real-source section is the 1.5.0 SQL-composition proof. It
  measures actual provider and case relations rather than claiming that
  the existence of scalar keys alone constitutes a
  supplementary-annotation system.
