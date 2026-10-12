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

Source revision: 68dfa0332ad9ea34eea556eddea231cd83d842f6; `src` tree 2de704e472ceeda276caa0eeef7093221378f3c4. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 66598566776aa44faffb68a8e0ba62961346c11db160962d115c30291d33586e. Host load before the run: 1.95, 1.08, 0.60. CPU affinity: 8-15. Empty-run overhead: 39.5 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                         | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:------------------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_af_table          | samples      |     1 |       1 |       8 |    6793144 | 1553     |          2.531 |       2.520 |       2.534 |           73.2 |        73.5 |      130.4 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       1 |      16 |   13586288 | 3112     |          3.794 |       3.783 |       3.804 |           89.7 |        89.8 |      170.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.248 |       6.238 |       6.351 |          133.3 |       133.4 |      249.9 | TRUE          |
| duckhts_roh_af_table          | samples      |     1 |       4 |       8 |    6793144 | 1553     |          2.529 |       2.518 |       2.530 |           82.5 |        85.2 |      164.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       4 |      16 |   13586288 | 3112     |          3.730 |       3.725 |       3.771 |          113.3 |       114.4 |      204.1 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.162 |       6.108 |       6.176 |          176.5 |       180.0 |      283.7 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       1 |       8 |    6793144 | 1553     |          3.524 |       3.509 |       3.529 |          108.1 |       109.2 |      449.6 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       1 |      16 |   13586288 | 3112     |          4.512 |       4.471 |       4.533 |          135.8 |       136.3 |      489.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.433 |       6.413 |       6.442 |          171.5 |       183.1 |      569.1 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       4 |       8 |    6793144 | 1553     |          3.489 |       3.483 |       3.499 |          144.6 |       145.6 |      483.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       4 |      16 |   13586288 | 3112     |          4.450 |       4.422 |       4.462 |          184.6 |       188.6 |      523.3 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.312 |       6.289 |       6.348 |          244.0 |       259.2 |      602.9 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.725 |       3.720 |       3.736 |          244.6 |       245.0 |      470.2 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.101 |       5.091 |       5.109 |          382.5 |       382.9 |      610.4 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.195 |       8.128 |       8.216 |          698.5 |       699.3 |      951.7 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.693 |       3.665 |       3.716 |          276.1 |       290.4 |      478.9 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.083 |       5.038 |       5.100 |          478.1 |       481.2 |      627.6 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.123 |       8.118 |       8.203 |          788.0 |       820.3 |      986.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.725 |       3.720 |       3.736 |          244.6 |       245.0 |      470.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       1 |       8 |    6793144 | 641      |          3.987 |       3.986 |       3.999 |          373.0 |       376.1 |      590.0 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.482 |       4.477 |       4.497 |          631.0 |       631.2 |      829.6 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.693 |       3.665 |       3.716 |          276.1 |       290.4 |      478.9 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.971 |       3.919 |       3.981 |          482.8 |       484.5 |      607.3 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.489 |       4.474 |       4.549 |          767.1 |       812.1 |      864.2 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.482 |       4.477 |       4.497 |          631.0 |       631.2 |      829.6 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.713 |       5.705 |       5.758 |          659.7 |       660.9 |      870.3 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.195 |       8.128 |       8.216 |          698.5 |       699.3 |      951.7 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.489 |       4.474 |       4.549 |          767.1 |       812.1 |      864.2 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.697 |       5.661 |       5.717 |          763.4 |       783.6 |      904.9 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.123 |       8.118 |       8.203 |          788.0 |       820.3 |      986.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.018 |       0.020 |           49.9 |        49.9 |       58.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.033 |       0.032 |       0.034 |           51.9 |        52.1 |       65.7 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.048 |       0.047 |       0.048 |           53.6 |        53.6 |       73.1 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.026 |       0.020 |       0.033 |           52.1 |        52.4 |       92.0 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.034 |       0.030 |       0.034 |           58.7 |        60.4 |       99.5 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.034 |       0.034 |       0.035 |           63.4 |        63.6 |      107.0 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.254 |       0.254 |       0.254 |           96.5 |        96.7 |      202.7 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.500 |       0.500 |       0.503 |          143.4 |       143.6 |      340.0 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.998 |       0.995 |       0.999 |          237.2 |       237.3 |      614.7 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.141 |       0.137 |       0.148 |          142.5 |       143.5 |      280.5 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.276 |       0.270 |       0.292 |          209.4 |       240.5 |      417.8 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.489 |       0.463 |       0.492 |          326.5 |       328.8 |      692.5 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.149 |       0.148 |       0.149 |          122.2 |       122.3 |      211.8 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.306 |       0.304 |       0.310 |          187.5 |       187.5 |      384.2 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.613 |       0.610 |       0.615 |          334.5 |       334.6 |      729.0 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.145 |       0.144 |       0.146 |          131.2 |       132.3 |      523.0 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.310 |       0.307 |       0.316 |          196.2 |       196.2 |     1006.5 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.607 |       0.597 |       0.610 |          338.7 |       346.3 |     1973.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.023 |       0.022 |       0.023 |           56.1 |        56.1 |       58.2 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.038 |       0.038 |       0.039 |           58.9 |        58.9 |       65.7 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.054 |       0.054 |       0.054 |           61.7 |        61.7 |       73.1 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.025 |       0.025 |       0.027 |           56.4 |        56.9 |       92.0 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.026 |       0.024 |       0.040 |           66.4 |        67.1 |       99.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.036 |       0.035 |       0.039 |           70.0 |        71.4 |      107.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.290 |       0.290 |       0.291 |          104.5 |       104.7 |      202.7 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.570 |       0.570 |       0.571 |          163.7 |       164.0 |      340.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       1 |      32 |    8000000 | 800      |          1.133 |       1.126 |       1.138 |          262.8 |       263.2 |      614.7 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.194 |       0.190 |       0.197 |          145.6 |       149.5 |      280.5 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.337 |       0.304 |       0.342 |          226.7 |       232.8 |      417.8 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.581 |       0.577 |       0.624 |          362.2 |       365.9 |      692.5 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.169 |       0.167 |       0.169 |          127.9 |       128.0 |      211.8 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.350 |       0.349 |       0.351 |          193.1 |       193.2 |      384.2 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.691 |       0.689 |       0.694 |          339.1 |       339.1 |      729.0 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.157 |       0.157 |       0.157 |          138.7 |       138.7 |      523.0 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.317 |       0.311 |       0.321 |          203.7 |       205.5 |     1006.5 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.630 |       0.618 |       0.631 |          347.0 |       348.3 |     1973.5 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                         | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:------------------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry          | joint     |       1 |            0.454 |            0.684 |        1.564 |        2.856 |         3.725 |
| duckhts_roh_ancestry          | reference |       1 |            0.098 |            0.169 |        1.525 |        2.580 |         3.725 |
| duckhts_roh_af_table          | samples   |       1 |            0.584 |            0.720 |        1.225 |        1.821 |         2.531 |
| duckhts_roh_af_table_bcftools | samples   |       1 |            0.357 |            0.512 |        1.256 |        1.586 |         3.524 |
| duckhts_roh_ancestry          | samples   |       1 |            0.350 |            0.520 |        1.045 |        1.107 |         4.482 |
| duckhts_roh_counts            | samples   |       1 |            0.978 |            0.996 |        1.486 |        2.458 |         0.254 |
| duckhts_roh_counts_fitted     | samples   |       1 |            0.976 |            0.989 |        1.567 |        2.515 |         0.290 |
| duckhts_roh_counts            | sites     |       1 |            1.040 |            1.001 |        1.534 |        2.737 |         0.149 |
| duckhts_roh_counts_fitted     | sites     |       1 |            1.054 |            0.981 |        1.510 |        2.651 |         0.169 |
| duckhts_roh_ancestry          | joint     |       4 |            0.461 |            0.676 |        1.732 |        2.854 |         3.693 |
| duckhts_roh_ancestry          | reference |       4 |            0.105 |            0.177 |        1.749 |        2.778 |         3.693 |
| duckhts_roh_af_table          | samples   |       4 |            0.561 |            0.724 |        1.373 |        2.139 |         2.529 |
| duckhts_roh_af_table_bcftools | samples   |       4 |            0.351 |            0.504 |        1.277 |        1.687 |         3.489 |
| duckhts_roh_ancestry          | samples   |       4 |            0.344 |            0.512 |        0.995 |        1.027 |         4.489 |
| duckhts_roh_counts            | samples   |       4 |            0.971 |            0.824 |        1.469 |        2.291 |         0.141 |
| duckhts_roh_counts_fitted     | samples   |       4 |            0.796 |            0.786 |        1.557 |        2.488 |         0.194 |
| duckhts_roh_counts            | sites     |       4 |            1.095 |            0.969 |        1.495 |        2.582 |         0.145 |
| duckhts_roh_counts_fitted     | sites     |       4 |            1.018 |            0.989 |        1.469 |        2.502 |         0.157 |

## Findings

Every macro stays within its declared budget at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 49.9 MiB
peak RSS and three samples in 0.048 s at
53.6 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_counts_fitted` returns the same runs as `duckhts_roh_counts` in
every cell, which the render checks, and stays within the same budget. Its
peak RSS over the scalar decode at the largest inputs, one thread, median of
three: 4.6 MiB with 4 million sites of one sample and
25.6 MiB with 32 samples of 250,000 sites.

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
6.248 s at 133.3 MiB
peak RSS (ceiling 249.9 MiB), and
`duckhts_roh_af_table_bcftools` takes 6.433 s
at 171.5 MiB (ceiling 569.1 MiB).
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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.725 |       244.605 |      0.596 |     0.137 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.101 |       382.531 |      0.504 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.195 |       698.543 |      0.450 |     0.106 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.725 |       244.605 |      0.596 |     0.137 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         3.987 |       372.992 |      0.609 |     0.173 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.482 |       630.973 |      0.618 |     0.232 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.482 |       630.973 |      0.618 |     0.232 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.713 |       659.734 |      0.526 |     0.163 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.195 |       698.543 |      0.450 |     0.106 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        49.863 |      0.856 |     0.916 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.033 |        51.945 |      0.813 |     0.795 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.048 |        53.637 |      0.797 |     0.677 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.254 |        96.520 |      0.773 |     0.379 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.500 |       143.387 |      0.757 |     0.307 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.998 |       237.227 |      0.746 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.149 |       122.207 |      0.839 |     0.752 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.306 |       187.488 |      0.906 |     0.604 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.613 |       334.520 |      0.872 |     0.476 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 699 MiB and
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
| duckhts_roh_ancestry      | joint        |     1 |          3.698 |        244.484 |         3.725 |       244.605 |      1.007 |     1.000 |
| duckhts_roh_ancestry      | joint        |     2 |          5.060 |        385.195 |         5.101 |       382.531 |      1.008 |     0.993 |
| duckhts_roh_ancestry      | joint        |     4 |          8.205 |        698.496 |         8.195 |       698.543 |      0.999 |     1.000 |
| duckhts_roh_ancestry      | reference    |     1 |          3.698 |        244.484 |         3.725 |       244.605 |      1.007 |     1.000 |
| duckhts_roh_ancestry      | reference    |     2 |          3.966 |        372.727 |         3.987 |       372.992 |      1.005 |     1.001 |
| duckhts_roh_ancestry      | reference    |     4 |          4.487 |        629.414 |         4.482 |       630.973 |      0.999 |     1.002 |
| duckhts_roh_ancestry      | samples      |     1 |          4.487 |        629.414 |         4.482 |       630.973 |      0.999 |     1.002 |
| duckhts_roh_ancestry      | samples      |     2 |          5.729 |        647.402 |         5.713 |       659.734 |      0.997 |     1.019 |
| duckhts_roh_ancestry      | samples      |     4 |          8.205 |        698.496 |         8.195 |       698.543 |      0.999 |     1.000 |
| duckhts_roh_counts        | real_samples |     1 |          0.019 |         49.262 |         0.019 |        49.863 |      0.971 |     1.012 |
| duckhts_roh_counts        | real_samples |     2 |          0.033 |         51.418 |         0.033 |        51.945 |      0.991 |     1.010 |
| duckhts_roh_counts        | real_samples |     3 |          0.048 |         52.871 |         0.048 |        53.637 |      0.990 |     1.014 |
| duckhts_roh_counts        | samples      |     1 |          0.253 |         96.043 |         0.254 |        96.520 |      1.004 |     1.005 |
| duckhts_roh_counts        | samples      |     2 |          0.499 |        142.938 |         0.500 |       143.387 |      1.003 |     1.003 |
| duckhts_roh_counts        | samples      |     4 |          0.991 |        236.582 |         0.998 |       237.227 |      1.007 |     1.003 |
| duckhts_roh_counts        | sites        |     1 |          0.150 |        121.785 |         0.149 |       122.207 |      0.994 |     1.003 |
| duckhts_roh_counts        | sites        |     2 |          0.306 |        186.762 |         0.306 |       187.488 |      0.999 |     1.004 |
| duckhts_roh_counts        | sites        |     4 |          0.612 |        333.688 |         0.613 |       334.520 |      1.001 |     1.002 |
| duckhts_roh_counts_fitted | real_samples |     1 |          0.022 |         55.723 |         0.023 |        56.129 |      1.011 |     1.007 |
| duckhts_roh_counts_fitted | real_samples |     2 |          0.038 |         58.586 |         0.038 |        58.879 |      0.981 |     1.005 |
| duckhts_roh_counts_fitted | real_samples |     3 |          0.054 |         61.102 |         0.054 |        61.699 |      1.007 |     1.010 |
| duckhts_roh_counts_fitted | samples      |     1 |          0.291 |        103.973 |         0.290 |       104.520 |      0.998 |     1.005 |
| duckhts_roh_counts_fitted | samples      |     2 |          0.570 |        159.020 |         0.570 |       163.703 |      1.001 |     1.029 |
| duckhts_roh_counts_fitted | samples      |     4 |          1.126 |        266.082 |         1.133 |       262.789 |      1.005 |     0.988 |
| duckhts_roh_counts_fitted | sites        |     1 |          0.170 |        127.102 |         0.169 |       127.895 |      0.993 |     1.006 |
| duckhts_roh_counts_fitted | sites        |     2 |          0.348 |        192.742 |         0.350 |       193.117 |      1.006 |     1.002 |
| duckhts_roh_counts_fitted | sites        |     4 |          0.688 |        338.613 |         0.691 |       339.062 |      1.005 |     1.001 |
