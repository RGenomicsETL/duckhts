ROH macro scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the
runs-of-homozygosity entry points added for \#318 and the per-sample fitted
count decode. Since \#329 the sites of each sample and chromosome are held in
native buffers that DuckHTS bounds, instead of DuckDB lists; the last sections
compare this revision with the list-aggregate revision and with the previous
render of this report.

- `duckhts_roh_counts` over a relation of read counts. Its growth dimensions are
  the sites per sample and chromosome (the length of one kernel list) and the
  number of samples. Synthetic counts give the 1×/2×/4× series. A real public
  workload adds read counts from three 1000 Genomes 30× CRAMs at the 108,757
  chr20 sites of the ROH ancestry evaluation, decoded for 1, 2 and 3 samples.
- `duckhts_roh_counts_fitted` over the same inputs, with a fit relation that
  gives every sample `seq_error` 0.001 and `contamination` 0, the scalar
  defaults. It is measured on the same series as `duckhts_roh_counts`.
- `duckhts_roh_ancestry` over the registered chr20 children BCF of the ROH
  ancestry evaluation (849,143 records). It is measured along three series:
  - the number of samples decoded together (8, 16, 32) with the full reference;
  - the reference rows, with 8 samples and the chr20 reference cut at position
    quantiles to a quarter, a half and all of its sites (the BCF is unchanged,
    so only the reference stage grows);
  - both together: (8, quarter), (16, half) and (32, all).

  A configuration shared by several series is measured once and listed in each.
- `duckhts_roh_af_table` and `duckhts_roh_af_table_bcftools` over the same BCF,
  with GT evidence (`gt_error := 30`) and a frequency relation built from the
  registered matched chr20 sites (`roh_ancestry_chr20_af_sites`, 108,757 rows,
  one per position; `af` is the first INFO/AF value). Both are measured over
  8, 16 and 32 samples decoded together. The two macros differ in how they
  choose the record of a listed position: `duckhts_roh_af_table` takes the
  record whose alleles match the row, and `duckhts_roh_af_table_bcftools`
  takes the record bcftools roh `--AF-file` takes, at the cost of a second,
  sites-only scan with `FILE_OFFSET` and a grouping-sets aggregate.

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`. The
JSON profiler gives query latency and DuckDB’s peak buffer and temporary
bytes, and `/usr/bin/time` gives peak RSS. `max_temp_directory_size` is zero,
so a spill fails the render. Each cell has three repetitions at one thread and
at four threads, run in an order that alternates sizes. Each query reduces the
runs it returns to a count and a total length, so the macro’s output is
consumed.

The synthetic read counts are deterministic. They are generated at render
time into a temporary directory, which is untimed staging. Each sample has a
2,000-site homozygous stretch every 10,000 sites over Hardy-Weinberg
background, with genotypes drawn at each site’s frequency and 30 reads per site,
so every input yields runs. The real counts are registry artifacts
`roh_counts_chr20_na18507`, `roh_counts_chr20_hg00403` and
`roh_counts_chr20_hg00188`: `duckhts_somalier_bam_counts` at the chr20
evaluation sites of the registered public CRAMs, read by range requests, with
MAPQ ≥ 20, base quality ≥ 13 and mate-overlap suppression
(`benchmarks/roh_counts_stage.R`). Each is checked against its registered row
count and counts digest before anything is measured. The ancestry
inputs are registry artifacts `roh_ancestry_chr20_children_bcf`,
`roh_ancestry_chr20_reference_long` and `roh_ancestry_chr20_q`. Their samples
are the first 8, 16 or 32 sample IDs of the registered proportions, sorted.
The reference cuts are the 25% and 50% position quantiles of its distinct sites.
The frequency relation is a view over the registry artifact
`roh_ancestry_chr20_af_sites`, checked against its registered SHA-256, so
each scan of it reads the Parquet file and nothing is held between queries.

The two frequency-relation macros return the same runs only where their rules
choose the same record. Before measuring, one untimed scan of the BCF without
samples finds, at each listed position, the records in file order
(`FILE_OFFSET`). The render checks that every such record and every row is a
single-base substitution with one ALT, so the indel-class skip of bcftools’
reader cannot act and bcftools’ rule reduces to “the first record of the
position claims it, and the position is skipped unless that record’s alleles
equal the row’s”. A position where that first record differs from the
record whose alleles match the row is one where the rules disagree. For each
sample set, an untimed fresh process then runs both macros over the relation
without those positions and asserts that they return identical runs, row for
row; over the full relation it counts the runs that differ and reports the
count. Every timed observation of these two macros also returns a digest of
its runs, and the render asserts that each macro’s digest is the same at both
thread counts and in every repetition.

`ROH_SCALING_SIZE=reduced` runs a driver check instead of the report: each
cell has one repetition instead of three, and the observations go to a
temporary file instead of `benchmark_roh_scaling_runs.csv`. The inputs keep
their full sizes, because the budgets leave out costs that do not grow with
the input: with lists a hundredth as long, the count decodes pass their
ceilings by a few MiB of fixed cost, and with 2 or 4 samples
`duckhts_roh_ancestry` takes a plan that peaks higher than with 8 (about 770
and 990 MiB against 640 MiB at one thread with the full reference). Every gate
and check still runs. A reduced run is not a measurement. The default is
`full`.

## Budgets, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the extension and run no scan.

Each ceiling is overhead + 3 × the live state below, and each is a gate: the
render stops if any observation exceeds its ceiling or spills.

- **`duckhts_roh_counts`.** The live state is:
  - every sample’s packed site list, 24 bytes per row (position, frequency and
    two counts), held natively until the scan ends and then by DuckDB as one
    BLOB per sample and chromosome;
  - the kernel workspace of each thread for the longest list, 36.25 bytes per
    site (position, two emissions, the forward-backward values and two
    traceback bits).

  That is rows × 24 + threads × sites per list × 36.25 bytes.
- **`duckhts_roh_counts_fitted`.** The same live state as `duckhts_roh_counts`.
  The fit lookup holds one row per sample, and the sample set is a second scan
  of the input that keeps one row per sample, so neither grows with the sites.
- **`duckhts_roh_ancestry`.** A record without a reference site is dropped
  before the sites are collected, so the lists hold reference sites only: the
  reference rows divided by its number of groups. The live state is:
  - the packed site lists, 16 bytes per sample and reference site (position,
    frequency and genotype evidence);
  - the distinct called sites and the hash table that joins the reference to
    them, 64 bytes per BCF record each;
  - the reference join and the per-site frequency lists, 64 bytes per reference
    row;
  - the kernel workspace, 36.25 bytes per reference site and thread.

  That is samples × reference sites × 16 + records × 128 + reference rows × 64 +
  threads × reference sites × 36.25 bytes.
- **`duckhts_roh_af_table`.** A record without a matching frequency row is
  dropped by an inner join before the sites are collected, and a row matches
  at most one record here, so each sample’s list holds at most one site per
  frequency row. The genotype scan streams into the lists. The live state is:
  - the packed site lists, 16 bytes per sample and frequency row;
  - the hash table that joins the frequency relation to the calls, built on
    the frequency side, 128 bytes per frequency row (the five typed columns,
    the hash and the table’s slots);
  - the kernel workspace, 36.25 bytes per frequency row and thread.

  That is samples × frequency rows × 16 + frequency rows × 128 +
  threads × frequency rows × 36.25 bytes.
- **`duckhts_roh_af_table_bcftools`.** The live state of
  `duckhts_roh_af_table`, plus the stages that choose each listed position’s
  record. The sites-only scan streams. With L the records at a listed POS
  (counted by the untimed scan above) the stages hold:
  - the distinct listed positions and the hash table of the left join that
    marks the sites-only scan with them, 64 bytes per frequency row each;
  - the frequency relation grouped per position and the hash table that joins
    it to the marked records, 128 bytes per frequency row each;
  - the duplicate-row guard’s grouping of the relation, 64 bytes per frequency
    row;
  - the grouping-sets aggregate, one group per contig and per record at a
    listed POS, 256 bytes per group (the hash entry, three offsets, a position,
    REF and the ALT list); the contig groups are negligible on one contig;
  - the window stages over that result (the next contig’s first POS, the
    contig’s and the position’s last offsets, the claiming row per position),
    each of which materializes its input at about 128 bytes per row, with two
    consecutive stages live while one feeds the next, so 256 bytes per listed
    record;
  - the hash table of the claimed records that the genotype scan probes,
    64 bytes per listed record.

  That is the `duckhts_roh_af_table` state + frequency rows × 448 + L × 576
  bytes.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_scaling.Rmd")'
# Driver check only; writes the page outside benchmarks/:
ROH_SCALING_SIZE=reduced taskset -c 0-3 Rscript -e \
  'rmarkdown::render("benchmarks/benchmark_roh_scaling.Rmd", output_dir = tempdir())'
```

Source revision: 2355ab28d18b9ab5ccaf6dad934087cea1b08488; `src` tree 9ada6d903f188bc52198290590ce9b0599c765ca. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 0682f5dded2dda4a5a3a720e615769d250029ad1a02449cb5d192ae959e61ce5. Host load before the run: 0.71, 1.03, 1.15. CPU affinity: 8-15. Empty-run overhead: 39.3 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                         | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:------------------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_af_table          | samples      |     1 |       1 |       8 |    6793144 | 1553     |          2.589 |       2.557 |       2.621 |           73.1 |        73.5 |      130.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       1 |      16 |   13586288 | 3112     |          3.789 |       3.777 |       3.970 |           89.2 |        89.2 |      170.0 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.300 |       6.285 |       6.427 |          134.9 |       134.9 |      249.7 | TRUE          |
| duckhts_roh_af_table          | samples      |     1 |       4 |       8 |    6793144 | 1553     |          2.535 |       2.532 |       2.579 |           82.7 |        86.9 |      164.0 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       4 |      16 |   13586288 | 3112     |          3.770 |       3.738 |       3.818 |          120.2 |       120.6 |      203.9 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.117 |       6.067 |       6.246 |          166.3 |       166.9 |      283.5 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       1 |       8 |    6793144 | 1553     |          3.543 |       3.537 |       3.667 |          109.4 |       112.0 |      449.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       1 |      16 |   13586288 | 3112     |          4.488 |       4.477 |       4.549 |          136.1 |       137.3 |      489.2 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.458 |       6.410 |       6.472 |          172.5 |       174.1 |      568.9 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       4 |       8 |    6793144 | 1553     |          3.490 |       3.488 |       3.556 |          141.4 |       146.0 |      483.2 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       4 |      16 |   13586288 | 3112     |          4.514 |       4.422 |       4.518 |          167.4 |       174.5 |      523.1 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.341 |       6.308 |       6.497 |          232.1 |       244.9 |      602.7 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.746 |       3.740 |       3.756 |          243.2 |       244.1 |      470.0 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.119 |       5.118 |       5.350 |          382.6 |       382.9 |      610.2 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.244 |       8.239 |       8.278 |          686.1 |       694.9 |      951.5 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.684 |       3.654 |       3.703 |          270.8 |       277.2 |      478.7 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.088 |       5.037 |       5.098 |          478.8 |       479.0 |      627.5 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.125 |       8.031 |       8.193 |          794.4 |       837.9 |      986.0 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.746 |       3.740 |       3.756 |          243.2 |       244.1 |      470.0 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       1 |       8 |    6793144 | 641      |          3.996 |       3.976 |       3.998 |          374.5 |       376.7 |      589.9 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.549 |       4.506 |       4.637 |          617.0 |       632.3 |      829.5 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.684 |       3.654 |       3.703 |          270.8 |       277.2 |      478.7 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.984 |       3.942 |       4.003 |          474.3 |       490.8 |      607.1 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.494 |       4.454 |       4.505 |          775.5 |       810.6 |      864.0 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.549 |       4.506 |       4.637 |          617.0 |       632.3 |      829.5 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.830 |       5.774 |       5.885 |          659.7 |       661.4 |      870.1 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.244 |       8.239 |       8.278 |          686.1 |       694.9 |      951.5 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.494 |       4.454 |       4.505 |          775.5 |       810.6 |      864.0 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.702 |       5.701 |       5.718 |          757.0 |       762.3 |      904.7 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.125 |       8.031 |       8.193 |          794.4 |       837.9 |      986.0 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.019 |       0.019 |           49.9 |        50.1 |       58.0 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.034 |       0.032 |       0.034 |           52.0 |        52.1 |       65.5 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.048 |       0.047 |       0.054 |           53.5 |        53.6 |       73.0 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.031 |       0.020 |       0.034 |           50.9 |        50.9 |       91.9 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.027 |       0.021 |       0.030 |           60.4 |        60.7 |       99.3 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.034 |       0.032 |       0.037 |           63.5 |        64.8 |      106.8 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.256 |       0.256 |       0.263 |           96.2 |        96.5 |      202.5 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.509 |       0.500 |       0.510 |          143.5 |       143.7 |      339.9 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.999 |       0.998 |       1.003 |          237.6 |       237.6 |      614.5 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.139 |       0.133 |       0.140 |          148.4 |       150.8 |      280.3 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.262 |       0.262 |       0.268 |          210.5 |       232.8 |      417.6 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.489 |       0.445 |       0.489 |          334.3 |       351.3 |      692.3 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.150 |       0.149 |       0.151 |          122.3 |       122.4 |      211.6 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.311 |       0.307 |       0.315 |          187.5 |       187.5 |      384.0 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.616 |       0.614 |       0.623 |          334.2 |       334.4 |      728.8 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.148 |       0.146 |       0.149 |          129.8 |       129.8 |      522.8 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.310 |       0.303 |       0.313 |          196.2 |       200.1 |     1006.3 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.602 |       0.586 |       0.626 |          337.3 |       342.8 |     1973.3 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.023 |       0.022 |       0.023 |           56.2 |        56.3 |       58.0 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.038 |       0.037 |       0.040 |           58.9 |        58.9 |       65.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.055 |       0.054 |       0.055 |           61.9 |        61.9 |       73.0 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.025 |       0.024 |       0.032 |           56.1 |        56.1 |       91.9 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.040 |       0.038 |       0.040 |           63.2 |        63.5 |       99.3 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.039 |       0.037 |       0.040 |           73.6 |        75.4 |      106.8 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.291 |       0.290 |       0.292 |          106.1 |       106.2 |      202.5 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.568 |       0.568 |       0.579 |          163.0 |       163.1 |      339.9 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       1 |      32 |    8000000 | 800      |          1.129 |       1.124 |       1.137 |          262.9 |       263.1 |      614.5 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.199 |       0.194 |       0.215 |          144.6 |       152.7 |      280.3 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.313 |       0.306 |       0.342 |          223.4 |       226.6 |      417.6 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.597 |       0.578 |       0.616 |          372.6 |       378.4 |      692.3 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.169 |       0.165 |       0.170 |          127.6 |       127.7 |      211.6 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.351 |       0.346 |       0.352 |          192.7 |       192.7 |      384.0 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.692 |       0.692 |       0.701 |          338.9 |       339.2 |      728.8 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.157 |       0.155 |       0.160 |          137.6 |       138.8 |      522.8 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.328 |       0.317 |       0.334 |          201.2 |       205.7 |     1006.3 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.623 |       0.618 |       0.626 |          347.5 |       350.3 |     1973.3 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                         | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:------------------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry          | joint     |       1 |            0.451 |            0.687 |        1.573 |        2.821 |         3.746 |
| duckhts_roh_ancestry          | reference |       1 |            0.093 |            0.187 |        1.540 |        2.537 |         3.746 |
| duckhts_roh_af_table          | samples   |       1 |            0.550 |            0.733 |        1.220 |        1.845 |         2.589 |
| duckhts_roh_af_table_bcftools | samples   |       1 |            0.341 |            0.525 |        1.244 |        1.577 |         3.543 |
| duckhts_roh_ancestry          | samples   |       1 |            0.358 |            0.500 |        1.069 |        1.112 |         4.549 |
| duckhts_roh_counts            | samples   |       1 |            0.992 |            0.972 |        1.492 |        2.470 |         0.256 |
| duckhts_roh_counts_fitted     | samples   |       1 |            0.964 |            0.990 |        1.536 |        2.478 |         0.291 |
| duckhts_roh_counts            | sites     |       1 |            1.053 |            0.988 |        1.533 |        2.733 |         0.150 |
| duckhts_roh_counts_fitted     | sites     |       1 |            1.050 |            0.982 |        1.510 |        2.656 |         0.169 |
| duckhts_roh_ancestry          | joint     |       4 |            0.466 |            0.675 |        1.768 |        2.934 |         3.684 |
| duckhts_roh_ancestry          | reference |       4 |            0.113 |            0.174 |        1.751 |        2.864 |         3.684 |
| duckhts_roh_af_table          | samples   |       4 |            0.573 |            0.698 |        1.453 |        2.011 |         2.535 |
| duckhts_roh_af_table_bcftools | samples   |       4 |            0.371 |            0.490 |        1.184 |        1.641 |         3.490 |
| duckhts_roh_ancestry          | samples   |       4 |            0.344 |            0.511 |        0.976 |        1.024 |         4.494 |
| duckhts_roh_counts            | samples   |       4 |            0.920 |            0.899 |        1.418 |        2.253 |         0.139 |
| duckhts_roh_counts_fitted     | samples   |       4 |            0.653 |            0.932 |        1.545 |        2.577 |         0.199 |
| duckhts_roh_counts            | sites     |       4 |            1.065 |            0.957 |        1.512 |        2.599 |         0.148 |
| duckhts_roh_counts_fitted     | sites     |       4 |            1.059 |            0.928 |        1.462 |        2.525 |         0.157 |

## Findings

Every macro stays within its declared budget at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 49.9 MiB
peak RSS and three samples in 0.048 s at
53.5 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_counts_fitted` returns the same runs as `duckhts_roh_counts` in
every cell, which the render checks, and stays within the same budget. Its
peak RSS over the scalar decode at the largest inputs, one thread, median of
three: 4.7 MiB with 4 million sites of one sample and
25.3 MiB with 32 samples of 250,000 sites.

`duckhts_roh_ancestry`: memory still grows with the number of samples decoded
together, because one pass over the BCF needs every sample’s sites. Each
sample now adds 16 bytes per reference site. `samples :=` chooses the batch, and
`max_sites` and `max_site_bytes` turn an oversized one into an explicit error.

`duckhts_roh_af_table` and `duckhts_roh_af_table_bcftools`: the relation has
108,757 rows, one per position, and the BCF has
109,105 records at those positions; at
347 positions it has
more than one record, as where a multiallelic site is written as several
biallelic records.
At 4 of those positions the first record in file order is not the record whose alleles match the row, so bcftools’ rule skips a position that the matching rule uses. Over the full relation the two macros therefore return different runs: 22 runs differ (counted in both directions) with 8 samples and 78 with 32 samples. The render reports these counts and does not assert equality there. Without the positions where the rules disagree, both
macros return identical runs for every sample set, which the render asserts:

| samples | agree_runs_matching | agree_runs_bcftools | agree_runs_differing | runs_matching | runs_bcftools | runs_differing |
|--------:|--------------------:|--------------------:|---------------------:|--------------:|--------------:|---------------:|
|       8 |                1553 |                1553 |                    0 |          1553 |          1553 |             22 |
|      16 |                3112 |                3112 |                    0 |          3112 |          3112 |             40 |
|      32 |                6037 |                6037 |                    0 |          6037 |          6037 |             78 |

With 32 samples, one thread, median of
3: `duckhts_roh_af_table` takes
6.3 s at 134.9 MiB
peak RSS (ceiling 249.7 MiB), and
`duckhts_roh_af_table_bcftools` takes 6.458 s
at 172.5 MiB (ceiling 568.9 MiB).
The difference is the cost of the bcftools site rule: the sites-only scan
with `FILE_OFFSET` and the stages over the records at listed positions.

## Change from the previous revision

The previous revision, a5f3b47a, held each sample’s sites in a
DuckDB list aggregate, kept the records without a reference site in those
lists, built the ancestry reference with ordered list aggregates and gave the
planner no row estimate for `read_bcf`. Its
observations are read from git; they were made with the same driver and inputs
on the same host. One thread, median of three:

| macro                | dimension    | scale | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|:---------------------|:-------------|------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.746 |       243.215 |      0.600 |     0.137 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.119 |       382.574 |      0.505 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.244 |       686.117 |      0.453 |     0.104 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.746 |       243.215 |      0.600 |     0.137 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         3.996 |       374.543 |      0.610 |     0.174 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.549 |       616.977 |      0.627 |     0.227 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.549 |       616.977 |      0.627 |     0.227 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.830 |       659.711 |      0.537 |     0.163 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.244 |       686.117 |      0.453 |     0.104 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        49.941 |      0.863 |     0.917 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.034 |        51.992 |      0.837 |     0.796 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.048 |        53.516 |      0.802 |     0.675 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.256 |        96.223 |      0.779 |     0.377 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.509 |       143.492 |      0.771 |     0.307 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.999 |       237.574 |      0.747 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.150 |       122.316 |      0.844 |     0.753 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.311 |       187.453 |      0.919 |     0.604 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.616 |       334.234 |      0.877 |     0.476 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 686 MiB and
from 18.2 s to 8.2 s.
The segments are unchanged: the chr20 evaluation arms reproduce their 51,083
(ancestry) and 60,204 (pooled) segments exactly.

## Change from the previous render

The previous render, 7d61c878, measured the same cells with the
same driver on the same host, with `taskset -c 8-15`. The macros it measured
are unchanged since; `read_bcf` gained an opt-in `FILE_OFFSET` column that
their scans do not request. `duckhts_roh_af_table` and
`duckhts_roh_af_table_bcftools` are new here and have no earlier cells. Host
load differs between renders. One thread, median of three:

| macro                     | dimension    | scale | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|:--------------------------|:-------------|------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
| duckhts_roh_ancestry      | joint        |     1 |          3.698 |        244.484 |         3.746 |       243.215 |      1.013 |     0.995 |
| duckhts_roh_ancestry      | joint        |     2 |          5.060 |        385.195 |         5.119 |       382.574 |      1.012 |     0.993 |
| duckhts_roh_ancestry      | joint        |     4 |          8.205 |        698.496 |         8.244 |       686.117 |      1.005 |     0.982 |
| duckhts_roh_ancestry      | reference    |     1 |          3.698 |        244.484 |         3.746 |       243.215 |      1.013 |     0.995 |
| duckhts_roh_ancestry      | reference    |     2 |          3.966 |        372.727 |         3.996 |       374.543 |      1.007 |     1.005 |
| duckhts_roh_ancestry      | reference    |     4 |          4.487 |        629.414 |         4.549 |       616.977 |      1.014 |     0.980 |
| duckhts_roh_ancestry      | samples      |     1 |          4.487 |        629.414 |         4.549 |       616.977 |      1.014 |     0.980 |
| duckhts_roh_ancestry      | samples      |     2 |          5.729 |        647.402 |         5.830 |       659.711 |      1.018 |     1.019 |
| duckhts_roh_ancestry      | samples      |     4 |          8.205 |        698.496 |         8.244 |       686.117 |      1.005 |     0.982 |
| duckhts_roh_counts        | real_samples |     1 |          0.019 |         49.262 |         0.019 |        49.941 |      0.979 |     1.014 |
| duckhts_roh_counts        | real_samples |     2 |          0.033 |         51.418 |         0.034 |        51.992 |      1.020 |     1.011 |
| duckhts_roh_counts        | real_samples |     3 |          0.048 |         52.871 |         0.048 |        53.516 |      0.996 |     1.012 |
| duckhts_roh_counts        | samples      |     1 |          0.253 |         96.043 |         0.256 |        96.223 |      1.013 |     1.002 |
| duckhts_roh_counts        | samples      |     2 |          0.499 |        142.938 |         0.509 |       143.492 |      1.021 |     1.004 |
| duckhts_roh_counts        | samples      |     4 |          0.991 |        236.582 |         0.999 |       237.574 |      1.008 |     1.004 |
| duckhts_roh_counts        | sites        |     1 |          0.150 |        121.785 |         0.150 |       122.316 |      1.000 |     1.004 |
| duckhts_roh_counts        | sites        |     2 |          0.306 |        186.762 |         0.311 |       187.453 |      1.014 |     1.004 |
| duckhts_roh_counts        | sites        |     4 |          0.612 |        333.688 |         0.616 |       334.234 |      1.006 |     1.002 |
| duckhts_roh_counts_fitted | real_samples |     1 |          0.022 |         55.723 |         0.023 |        56.238 |      1.015 |     1.009 |
| duckhts_roh_counts_fitted | real_samples |     2 |          0.038 |         58.586 |         0.038 |        58.871 |      0.996 |     1.005 |
| duckhts_roh_counts_fitted | real_samples |     3 |          0.054 |         61.102 |         0.055 |        61.855 |      1.013 |     1.012 |
| duckhts_roh_counts_fitted | samples      |     1 |          0.291 |        103.973 |         0.291 |       106.109 |      1.002 |     1.021 |
| duckhts_roh_counts_fitted | samples      |     2 |          0.570 |        159.020 |         0.568 |       162.965 |      0.998 |     1.025 |
| duckhts_roh_counts_fitted | samples      |     4 |          1.126 |        266.082 |         1.129 |       262.879 |      1.002 |     0.988 |
| duckhts_roh_counts_fitted | sites        |     1 |          0.170 |        127.102 |         0.169 |       127.586 |      0.998 |     1.004 |
| duckhts_roh_counts_fitted | sites        |     2 |          0.348 |        192.742 |         0.351 |       192.664 |      1.007 |     1.000 |
| duckhts_roh_counts_fitted | sites        |     4 |          0.688 |        338.613 |         0.692 |       338.945 |      1.007 |     1.001 |
