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
  gives every sample `seq_error` 0.001, `contamination` 0, `allele_balance`
  0.5 and `artefact_weight` 0, the scalar defaults, and `artefact_weight_sd`
  0 (it still decodes each list three times for `stable_fraction`). It is measured on the same
  series as `duckhts_roh_counts`.
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

Source revision: e7667ff353725f27c1dc1c64732239fdd5a9ec21; `src` tree 0340e562109701f322a733450b16be60a5daba10. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: cfcb1394601761f2a5ebd1a0d6745fea6099fc6b13232b388accb88938decf97. Host load before the run: 1.62, 1.08, 1.26. CPU affinity: 8-15. Empty-run overhead: 39.5 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                         | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:------------------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_af_table          | samples      |     1 |       1 |       8 |    6793144 | 1553     |          2.570 |       2.562 |       2.572 |           73.4 |        73.5 |      130.4 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       1 |      16 |   13586288 | 3112     |          3.791 |       3.769 |       3.795 |           91.3 |        91.4 |      170.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.261 |       6.205 |       6.311 |          128.5 |       131.7 |      249.9 | TRUE          |
| duckhts_roh_af_table          | samples      |     1 |       4 |       8 |    6793144 | 1553     |          2.566 |       2.544 |       2.579 |           84.0 |        84.5 |      164.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       4 |      16 |   13586288 | 3112     |          3.740 |       3.701 |       3.806 |          108.8 |       110.3 |      204.1 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.247 |       6.190 |       6.252 |          172.0 |       177.2 |      283.7 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       1 |       8 |    6793144 | 1553     |          3.562 |       3.551 |       3.567 |          110.3 |       111.7 |      449.6 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       1 |      16 |   13586288 | 3112     |          4.486 |       4.459 |       4.503 |          137.0 |       137.2 |      489.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.475 |       6.455 |       6.498 |          180.1 |       180.3 |      569.1 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       4 |       8 |    6793144 | 1553     |          3.554 |       3.545 |       3.559 |          144.2 |       148.6 |      483.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       4 |      16 |   13586288 | 3112     |          4.495 |       4.443 |       4.529 |          160.4 |       172.4 |      523.3 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.399 |       6.376 |       6.408 |          224.8 |       251.3 |      602.9 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.780 |       3.762 |       3.782 |          246.0 |       246.5 |      470.2 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.102 |       5.065 |       5.142 |          381.1 |       382.6 |      610.4 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.206 |       8.185 |       8.208 |          696.6 |       698.9 |      951.7 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.709 |       3.702 |       3.730 |          272.1 |       286.7 |      478.9 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.055 |       5.035 |       5.085 |          485.3 |       485.8 |      627.6 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.110 |       8.060 |       8.111 |          791.3 |       792.2 |      986.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.780 |       3.762 |       3.782 |          246.0 |       246.5 |      470.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       1 |       8 |    6793144 | 641      |          4.015 |       4.012 |       4.017 |          374.3 |       374.4 |      590.0 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.550 |       4.524 |       4.559 |          630.0 |       631.6 |      829.6 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.709 |       3.702 |       3.730 |          272.1 |       286.7 |      478.9 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.967 |       3.966 |       3.982 |          476.1 |       488.5 |      607.3 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.508 |       4.504 |       4.534 |          767.5 |       773.1 |      864.2 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.550 |       4.524 |       4.559 |          630.0 |       631.6 |      829.6 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.735 |       5.706 |       5.779 |          659.2 |       659.8 |      870.3 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.206 |       8.185 |       8.208 |          696.6 |       698.9 |      951.7 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.508 |       4.504 |       4.534 |          767.5 |       773.1 |      864.2 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.674 |       5.635 |       5.674 |          755.8 |       760.9 |      904.9 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.110 |       8.060 |       8.111 |          791.3 |       792.2 |      986.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.018 |       0.019 |           50.1 |        50.3 |       58.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.033 |       0.032 |       0.033 |           52.1 |        52.4 |       65.7 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.047 |       0.047 |       0.048 |           53.8 |        54.0 |       73.1 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.025 |       0.020 |       0.036 |           51.1 |        52.5 |       92.0 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.032 |       0.021 |       0.032 |           58.8 |        61.1 |       99.5 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.035 |       0.035 |       0.041 |           63.7 |        66.0 |      107.0 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.254 |       0.251 |       0.256 |           97.1 |        97.2 |      202.7 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.499 |       0.498 |       0.506 |          143.7 |       143.7 |      340.0 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.996 |       0.991 |       0.998 |          237.2 |       237.7 |      614.7 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.158 |       0.134 |       0.161 |          134.3 |       138.3 |      280.5 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.277 |       0.273 |       0.281 |          209.5 |       229.8 |      417.8 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.479 |       0.478 |       0.530 |          341.8 |       357.1 |      692.5 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.151 |       0.150 |       0.152 |          122.5 |       122.7 |      211.8 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.309 |       0.307 |       0.314 |          187.8 |       187.9 |      384.2 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.617 |       0.612 |       0.623 |          334.3 |       334.4 |      729.0 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.148 |       0.145 |       0.148 |          132.6 |       134.3 |      523.0 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.308 |       0.306 |       0.318 |          196.3 |       197.9 |     1006.5 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.611 |       0.610 |       0.612 |          336.0 |       336.4 |     1973.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.039 |       0.039 |       0.040 |           56.8 |        57.0 |       58.2 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.072 |       0.071 |       0.072 |           59.5 |        59.6 |       65.7 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.104 |       0.103 |       0.104 |           62.2 |        62.3 |       73.1 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.050 |       0.042 |       0.050 |           58.3 |        58.6 |       92.0 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.044 |       0.042 |       0.050 |           68.7 |        69.9 |       99.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.062 |       0.060 |       0.069 |           76.3 |        77.8 |      107.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.583 |       0.582 |       0.586 |          106.6 |       106.8 |      202.7 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       1 |      16 |    4000000 | 400      |          1.149 |       1.149 |       1.158 |          172.3 |       172.5 |      340.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       1 |      32 |    8000000 | 800      |          2.293 |       2.288 |       2.302 |          263.1 |       263.2 |      614.7 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.385 |       0.381 |       0.419 |          159.7 |       166.1 |      280.5 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.568 |       0.545 |       0.630 |          242.2 |       243.9 |      417.8 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       4 |      32 |    8000000 | 800      |          1.148 |       0.920 |       1.154 |          347.2 |       375.7 |      692.5 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.316 |       0.314 |       0.316 |          127.9 |       128.2 |      211.8 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.639 |       0.629 |       0.643 |          193.2 |       193.3 |      384.2 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       1 |       1 |    4000000 | 400      |          1.269 |       1.268 |       1.280 |          339.6 |       339.7 |      729.0 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.299 |       0.296 |       0.300 |          139.0 |       141.2 |      523.0 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.598 |       0.596 |       0.606 |          196.8 |       198.0 |     1006.5 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       4 |       1 |    4000000 | 400      |          1.190 |       1.182 |       1.197 |          341.0 |       342.7 |     1973.5 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                         | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:------------------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry          | joint     |       1 |            0.433 |            0.686 |        1.549 |        2.832 |         3.780 |
| duckhts_roh_ancestry          | reference |       1 |            0.087 |            0.181 |        1.522 |        2.561 |         3.780 |
| duckhts_roh_af_table          | samples   |       1 |            0.561 |            0.724 |        1.244 |        1.751 |         2.570 |
| duckhts_roh_af_table_bcftools | samples   |       1 |            0.333 |            0.529 |        1.242 |        1.633 |         3.562 |
| duckhts_roh_ancestry          | samples   |       1 |            0.334 |            0.517 |        1.046 |        1.106 |         4.550 |
| duckhts_roh_counts            | samples   |       1 |            0.977 |            0.997 |        1.480 |        2.443 |         0.254 |
| duckhts_roh_counts_fitted     | samples   |       1 |            0.979 |            0.996 |        1.616 |        2.468 |         0.583 |
| duckhts_roh_counts            | sites     |       1 |            1.036 |            0.997 |        1.533 |        2.729 |         0.151 |
| duckhts_roh_counts_fitted     | sites     |       1 |            1.017 |            0.989 |        1.511 |        2.655 |         0.316 |
| duckhts_roh_ancestry          | joint     |       4 |            0.447 |            0.682 |        1.784 |        2.908 |         3.709 |
| duckhts_roh_ancestry          | reference |       4 |            0.097 |            0.185 |        1.750 |        2.821 |         3.709 |
| duckhts_roh_af_table          | samples   |       4 |            0.544 |            0.740 |        1.295 |        2.048 |         2.566 |
| duckhts_roh_af_table_bcftools | samples   |       4 |            0.339 |            0.510 |        1.112 |        1.559 |         3.554 |
| duckhts_roh_ancestry          | samples   |       4 |            0.332 |            0.515 |        0.985 |        1.031 |         4.508 |
| duckhts_roh_counts            | samples   |       4 |            0.810 |            0.789 |        1.560 |        2.545 |         0.158 |
| duckhts_roh_counts_fitted     | samples   |       4 |            0.560 |            1.015 |        1.517 |        2.174 |         0.385 |
| duckhts_roh_counts            | sites     |       4 |            1.058 |            0.989 |        1.480 |        2.534 |         0.148 |
| duckhts_roh_counts_fitted     | sites     |       4 |            1.001 |            0.992 |        1.416 |        2.453 |         0.299 |

## Findings

Every macro stays within its declared budget at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 50.1 MiB
peak RSS and three samples in 0.047 s at
53.8 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_counts_fitted` returns the same runs as `duckhts_roh_counts` in
every cell, which the render checks, and stays within the same budget. Its
peak RSS over the scalar decode at the largest inputs, one thread, median of
three: 5.3 MiB with 4 million sites of one sample and
25.9 MiB with 32 samples of 250,000 sites.

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
6.261 s at 128.5 MiB
peak RSS (ceiling 249.9 MiB), and
`duckhts_roh_af_table_bcftools` takes 6.475 s
at 180.1 MiB (ceiling 569.1 MiB).
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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.780 |       246.008 |      0.605 |     0.138 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.102 |       381.094 |      0.504 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.206 |       696.629 |      0.450 |     0.106 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.780 |       246.008 |      0.605 |     0.138 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         4.015 |       374.277 |      0.613 |     0.174 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.550 |       630.016 |      0.627 |     0.232 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.550 |       630.016 |      0.627 |     0.232 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.735 |       659.188 |      0.528 |     0.163 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.206 |       696.629 |      0.450 |     0.106 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        50.145 |      0.856 |     0.921 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.033 |        52.117 |      0.812 |     0.798 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.047 |        53.848 |      0.788 |     0.679 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.254 |        97.094 |      0.772 |     0.381 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.499 |       143.680 |      0.756 |     0.308 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.996 |       237.191 |      0.745 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.151 |       122.496 |      0.849 |     0.754 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.309 |       187.844 |      0.914 |     0.605 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.617 |       334.348 |      0.878 |     0.476 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 697 MiB and
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
| duckhts_roh_ancestry      | joint        |     1 |          3.698 |        244.484 |         3.780 |       246.008 |      1.022 |     1.006 |
| duckhts_roh_ancestry      | joint        |     2 |          5.060 |        385.195 |         5.102 |       381.094 |      1.008 |     0.989 |
| duckhts_roh_ancestry      | joint        |     4 |          8.205 |        698.496 |         8.206 |       696.629 |      1.000 |     0.997 |
| duckhts_roh_ancestry      | reference    |     1 |          3.698 |        244.484 |         3.780 |       246.008 |      1.022 |     1.006 |
| duckhts_roh_ancestry      | reference    |     2 |          3.966 |        372.727 |         4.015 |       374.277 |      1.012 |     1.004 |
| duckhts_roh_ancestry      | reference    |     4 |          4.487 |        629.414 |         4.550 |       630.016 |      1.014 |     1.001 |
| duckhts_roh_ancestry      | samples      |     1 |          4.487 |        629.414 |         4.550 |       630.016 |      1.014 |     1.001 |
| duckhts_roh_ancestry      | samples      |     2 |          5.729 |        647.402 |         5.735 |       659.188 |      1.001 |     1.018 |
| duckhts_roh_ancestry      | samples      |     4 |          8.205 |        698.496 |         8.206 |       696.629 |      1.000 |     0.997 |
| duckhts_roh_counts        | real_samples |     1 |          0.019 |         49.262 |         0.019 |        50.145 |      0.971 |     1.018 |
| duckhts_roh_counts        | real_samples |     2 |          0.033 |         51.418 |         0.033 |        52.117 |      0.990 |     1.014 |
| duckhts_roh_counts        | real_samples |     3 |          0.048 |         52.871 |         0.047 |        53.848 |      0.979 |     1.018 |
| duckhts_roh_counts        | samples      |     1 |          0.253 |         96.043 |         0.254 |        97.094 |      1.003 |     1.011 |
| duckhts_roh_counts        | samples      |     2 |          0.499 |        142.938 |         0.499 |       143.680 |      1.000 |     1.005 |
| duckhts_roh_counts        | samples      |     4 |          0.991 |        236.582 |         0.996 |       237.191 |      1.005 |     1.003 |
| duckhts_roh_counts        | sites        |     1 |          0.150 |        121.785 |         0.151 |       122.496 |      1.005 |     1.006 |
| duckhts_roh_counts        | sites        |     2 |          0.306 |        186.762 |         0.309 |       187.844 |      1.008 |     1.006 |
| duckhts_roh_counts        | sites        |     4 |          0.612 |        333.688 |         0.617 |       334.348 |      1.007 |     1.002 |
| duckhts_roh_counts_fitted | real_samples |     1 |          0.022 |         55.723 |         0.039 |        56.773 |      1.767 |     1.019 |
| duckhts_roh_counts_fitted | real_samples |     2 |          0.038 |         58.586 |         0.072 |        59.453 |      1.862 |     1.015 |
| duckhts_roh_counts_fitted | real_samples |     3 |          0.054 |         61.102 |         0.104 |        62.246 |      1.920 |     1.019 |
| duckhts_roh_counts_fitted | samples      |     1 |          0.291 |        103.973 |         0.583 |       106.605 |      2.006 |     1.025 |
| duckhts_roh_counts_fitted | samples      |     2 |          0.570 |        159.020 |         1.149 |       172.320 |      2.018 |     1.084 |
| duckhts_roh_counts_fitted | samples      |     4 |          1.126 |        266.082 |         2.293 |       263.074 |      2.036 |     0.989 |
| duckhts_roh_counts_fitted | sites        |     1 |          0.170 |        127.102 |         0.316 |       127.941 |      1.862 |     1.007 |
| duckhts_roh_counts_fitted | sites        |     2 |          0.348 |        192.742 |         0.639 |       193.168 |      1.836 |     1.002 |
| duckhts_roh_counts_fitted | sites        |     4 |          0.688 |        338.613 |         1.269 |       339.562 |      1.845 |     1.003 |
