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

Source revision: 6a651f6aa1cdbc03d4d9692ee5d05b9150be8be5; `src` tree 0f829a860b92dbf98069944479775481bf9d43e8. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 2693fbb5f1b4d7c2b7f6ae7b08aec19919d3f1864dfae4051c68cb79987a2dc2. Host load before the run: 1.30, 1.19, 1.43. CPU affinity: 8-15. Empty-run overhead: 39.4 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                         | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:------------------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_af_table          | samples      |     1 |       1 |       8 |    6793144 | 1553     |          2.535 |       2.526 |       2.553 |           73.4 |        73.4 |      130.4 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       1 |      16 |   13586288 | 3112     |          3.816 |       3.777 |       3.818 |           89.9 |        90.0 |      170.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.246 |       6.198 |       6.290 |          131.9 |       133.5 |      249.8 | TRUE          |
| duckhts_roh_af_table          | samples      |     1 |       4 |       8 |    6793144 | 1553     |          2.533 |       2.527 |       2.590 |           83.2 |        84.7 |      164.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       4 |      16 |   13586288 | 3112     |          3.776 |       3.685 |       4.292 |          110.8 |       121.0 |      204.0 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.278 |       6.193 |       6.414 |          167.1 |       189.5 |      283.7 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       1 |       8 |    6793144 | 1553     |          3.507 |       3.482 |       3.528 |          117.7 |       117.8 |      449.6 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       1 |      16 |   13586288 | 3112     |          4.544 |       4.493 |       4.552 |          135.3 |       136.8 |      489.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.514 |       6.412 |       6.530 |          169.3 |       181.6 |      569.0 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       4 |       8 |    6793144 | 1553     |          3.562 |       3.544 |       3.589 |          145.1 |       145.8 |      483.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       4 |      16 |   13586288 | 3112     |          4.487 |       4.478 |       5.011 |          166.6 |       184.2 |      523.2 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.384 |       6.372 |       6.562 |          234.8 |       245.5 |      602.9 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.749 |       3.743 |       3.762 |          245.2 |       245.5 |      470.2 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.127 |       5.105 |       5.133 |          381.8 |       382.3 |      610.3 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.240 |       8.233 |       8.268 |          697.8 |       701.6 |      951.6 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.688 |       3.676 |       3.738 |          286.8 |       286.8 |      478.8 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.082 |       5.040 |       5.098 |          467.4 |       475.4 |      627.6 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.234 |       8.148 |       8.276 |          799.9 |       815.6 |      986.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.749 |       3.743 |       3.762 |          245.2 |       245.5 |      470.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       1 |       8 |    6793144 | 641      |          3.994 |       3.994 |       4.029 |          373.7 |       374.6 |      590.0 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.520 |       4.510 |       4.532 |          630.8 |       631.5 |      829.6 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.688 |       3.676 |       3.738 |          286.8 |       286.8 |      478.8 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       4 |       8 |    6793144 | 641      |          3.963 |       3.961 |       3.974 |          477.1 |       479.8 |      607.3 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.500 |       4.472 |       4.517 |          761.8 |       771.8 |      864.2 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.520 |       4.510 |       4.532 |          630.8 |       631.5 |      829.6 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.787 |       5.747 |       5.800 |          648.0 |       657.3 |      870.3 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.240 |       8.233 |       8.268 |          697.8 |       701.6 |      951.6 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.500 |       4.472 |       4.517 |          761.8 |       771.8 |      864.2 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       4 |      16 |   13586288 | 2637     |          5.720 |       5.665 |       5.728 |          754.5 |       767.8 |      904.9 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.234 |       8.148 |       8.276 |          799.9 |       815.6 |      986.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.018 |       0.019 |           50.2 |        50.3 |       58.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.033 |       0.032 |       0.034 |           52.5 |        52.5 |       65.6 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.048 |       0.048 |       0.048 |           53.7 |        53.9 |       73.1 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.026 |       0.025 |       0.030 |           51.2 |        52.2 |       92.0 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.032 |       0.025 |       0.036 |           60.7 |        61.1 |       99.5 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.035 |       0.034 |       0.041 |           64.8 |        67.3 |      107.0 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.254 |       0.251 |       0.255 |           96.9 |        97.1 |      202.7 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.503 |       0.502 |       0.505 |          143.7 |       143.9 |      340.0 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       1 |      32 |    8000000 | 800      |          0.993 |       0.992 |       0.998 |          237.6 |       237.8 |      614.7 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.166 |       0.152 |       0.183 |          143.8 |       144.6 |      280.5 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.275 |       0.261 |       0.280 |          209.0 |       210.6 |      417.8 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.462 |       0.446 |       0.474 |          335.2 |       348.9 |      692.5 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.149 |       0.148 |       0.151 |          122.3 |       122.5 |      211.8 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.309 |       0.307 |       0.309 |          187.7 |       187.7 |      384.2 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.620 |       0.620 |       0.621 |          334.6 |       334.7 |      728.9 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.148 |       0.147 |       0.148 |          129.9 |       134.5 |      522.9 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.306 |       0.305 |       0.313 |          196.3 |       196.4 |     1006.5 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.610 |       0.596 |       0.628 |          336.9 |       342.5 |     1973.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.040 |       0.040 |       0.041 |           56.6 |        56.9 |       58.2 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.073 |       0.073 |       0.073 |           59.3 |        59.4 |       65.6 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.105 |       0.103 |       0.109 |           62.2 |        62.3 |       73.1 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.042 |       0.041 |       0.044 |           58.2 |        58.3 |       92.0 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.051 |       0.051 |       0.059 |           68.3 |        70.4 |       99.5 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.064 |       0.062 |       0.064 |           76.6 |        77.5 |      107.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.586 |       0.578 |       0.588 |          104.8 |       106.4 |      202.7 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       1 |      16 |    4000000 | 400      |          1.154 |       1.146 |       1.156 |          167.8 |       168.1 |      340.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       1 |      32 |    8000000 | 800      |          2.304 |       2.296 |       2.310 |          264.4 |       264.5 |      614.7 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.395 |       0.387 |       0.424 |          155.0 |       157.1 |      280.5 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.542 |       0.533 |       0.628 |          237.0 |       241.5 |      417.8 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       4 |      32 |    8000000 | 800      |          1.143 |       1.019 |       1.196 |          368.9 |       391.9 |      692.5 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.318 |       0.316 |       0.320 |          128.0 |       128.1 |      211.8 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.644 |       0.644 |       0.650 |          193.2 |       193.2 |      384.2 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       1 |       1 |    4000000 | 400      |          1.286 |       1.274 |       1.290 |          339.5 |       339.5 |      728.9 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.299 |       0.297 |       0.300 |          139.6 |       141.2 |      522.9 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.600 |       0.597 |       0.604 |          196.9 |       197.8 |     1006.5 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       4 |       1 |    4000000 | 400      |          1.183 |       1.178 |       1.185 |          342.8 |       344.0 |     1973.5 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                         | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:------------------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry          | joint     |       1 |            0.452 |            0.685 |        1.557 |        2.846 |         3.749 |
| duckhts_roh_ancestry          | reference |       1 |            0.091 |            0.179 |        1.524 |        2.573 |         3.749 |
| duckhts_roh_af_table          | samples   |       1 |            0.590 |            0.711 |        1.225 |        1.797 |         2.535 |
| duckhts_roh_af_table_bcftools | samples   |       1 |            0.374 |            0.520 |        1.150 |        1.438 |         3.507 |
| duckhts_roh_ancestry          | samples   |       1 |            0.356 |            0.510 |        1.027 |        1.106 |         4.520 |
| duckhts_roh_counts            | samples   |       1 |            0.988 |            0.982 |        1.483 |        2.452 |         0.254 |
| duckhts_roh_counts_fitted     | samples   |       1 |            0.977 |            0.998 |        1.601 |        2.523 |         0.586 |
| duckhts_roh_counts            | sites     |       1 |            1.051 |            1.005 |        1.535 |        2.736 |         0.149 |
| duckhts_roh_counts_fitted     | sites     |       1 |            1.019 |            0.998 |        1.509 |        2.652 |         0.318 |
| duckhts_roh_ancestry          | joint     |       4 |            0.463 |            0.696 |        1.630 |        2.789 |         3.688 |
| duckhts_roh_ancestry          | reference |       4 |            0.104 |            0.183 |        1.664 |        2.656 |         3.688 |
| duckhts_roh_af_table          | samples   |       4 |            0.576 |            0.733 |        1.332 |        2.008 |         2.533 |
| duckhts_roh_af_table_bcftools | samples   |       4 |            0.333 |            0.509 |        1.148 |        1.618 |         3.562 |
| duckhts_roh_ancestry          | samples   |       4 |            0.346 |            0.525 |        0.990 |        1.050 |         4.500 |
| duckhts_roh_counts            | samples   |       4 |            0.731 |            0.751 |        1.453 |        2.331 |         0.166 |
| duckhts_roh_counts_fitted     | samples   |       4 |            0.454 |            1.078 |        1.529 |        2.380 |         0.395 |
| duckhts_roh_counts            | sites     |       4 |            1.052 |            0.995 |        1.511 |        2.594 |         0.148 |
| duckhts_roh_counts_fitted     | sites     |       4 |            1.003 |            0.980 |        1.410 |        2.456 |         0.299 |

## Findings

Every macro stays within its declared budget at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 50.2 MiB
peak RSS and three samples in 0.048 s at
53.7 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_counts_fitted` returns the same runs as `duckhts_roh_counts` in
every cell, which the render checks, and stays within the same budget. Its
peak RSS over the scalar decode at the largest inputs, one thread, median of
three: 4.9 MiB with 4 million sites of one sample and
26.8 MiB with 32 samples of 250,000 sites.

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
6.246 s at 131.9 MiB
peak RSS (ceiling 249.8 MiB), and
`duckhts_roh_af_table_bcftools` takes 6.514 s
at 169.3 MiB (ceiling 569 MiB).
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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.749 |       245.180 |      0.600 |     0.138 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.127 |       381.805 |      0.506 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.240 |       697.816 |      0.452 |     0.106 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.749 |       245.180 |      0.600 |     0.138 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         3.994 |       373.723 |      0.610 |     0.173 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.520 |       630.840 |      0.623 |     0.232 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.520 |       630.840 |      0.623 |     0.232 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.787 |       647.965 |      0.533 |     0.160 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.240 |       697.816 |      0.452 |     0.106 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        50.195 |      0.864 |     0.922 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.033 |        52.457 |      0.836 |     0.803 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.048 |        53.703 |      0.805 |     0.677 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.254 |        96.938 |      0.771 |     0.380 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.503 |       143.668 |      0.761 |     0.308 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         0.993 |       237.594 |      0.743 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.149 |       122.344 |      0.841 |     0.753 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.309 |       187.695 |      0.915 |     0.605 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.620 |       334.594 |      0.883 |     0.477 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 698 MiB and
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
| duckhts_roh_ancestry      | joint        |     1 |          3.698 |        244.484 |         3.749 |       245.180 |      1.014 |     1.003 |
| duckhts_roh_ancestry      | joint        |     2 |          5.060 |        385.195 |         5.127 |       381.805 |      1.013 |     0.991 |
| duckhts_roh_ancestry      | joint        |     4 |          8.205 |        698.496 |         8.240 |       697.816 |      1.004 |     0.999 |
| duckhts_roh_ancestry      | reference    |     1 |          3.698 |        244.484 |         3.749 |       245.180 |      1.014 |     1.003 |
| duckhts_roh_ancestry      | reference    |     2 |          3.966 |        372.727 |         3.994 |       373.723 |      1.007 |     1.003 |
| duckhts_roh_ancestry      | reference    |     4 |          4.487 |        629.414 |         4.520 |       630.840 |      1.007 |     1.002 |
| duckhts_roh_ancestry      | samples      |     1 |          4.487 |        629.414 |         4.520 |       630.840 |      1.007 |     1.002 |
| duckhts_roh_ancestry      | samples      |     2 |          5.729 |        647.402 |         5.787 |       647.965 |      1.010 |     1.001 |
| duckhts_roh_ancestry      | samples      |     4 |          8.205 |        698.496 |         8.240 |       697.816 |      1.004 |     0.999 |
| duckhts_roh_counts        | real_samples |     1 |          0.019 |         49.262 |         0.019 |        50.195 |      0.981 |     1.019 |
| duckhts_roh_counts        | real_samples |     2 |          0.033 |         51.418 |         0.033 |        52.457 |      1.020 |     1.020 |
| duckhts_roh_counts        | real_samples |     3 |          0.048 |         52.871 |         0.048 |        53.703 |      0.999 |     1.016 |
| duckhts_roh_counts        | samples      |     1 |          0.253 |         96.043 |         0.254 |        96.938 |      1.003 |     1.009 |
| duckhts_roh_counts        | samples      |     2 |          0.499 |        142.938 |         0.503 |       143.668 |      1.008 |     1.005 |
| duckhts_roh_counts        | samples      |     4 |          0.991 |        236.582 |         0.993 |       237.594 |      1.003 |     1.004 |
| duckhts_roh_counts        | sites        |     1 |          0.150 |        121.785 |         0.149 |       122.344 |      0.996 |     1.005 |
| duckhts_roh_counts        | sites        |     2 |          0.306 |        186.762 |         0.309 |       187.695 |      1.009 |     1.005 |
| duckhts_roh_counts        | sites        |     4 |          0.612 |        333.688 |         0.620 |       334.594 |      1.013 |     1.003 |
| duckhts_roh_counts_fitted | real_samples |     1 |          0.022 |         55.723 |         0.040 |        56.582 |      1.805 |     1.015 |
| duckhts_roh_counts_fitted | real_samples |     2 |          0.038 |         58.586 |         0.073 |        59.348 |      1.893 |     1.013 |
| duckhts_roh_counts_fitted | real_samples |     3 |          0.054 |         61.102 |         0.105 |        62.219 |      1.940 |     1.018 |
| duckhts_roh_counts_fitted | samples      |     1 |          0.291 |        103.973 |         0.586 |       104.770 |      2.016 |     1.008 |
| duckhts_roh_counts_fitted | samples      |     2 |          0.570 |        159.020 |         1.154 |       167.809 |      2.025 |     1.055 |
| duckhts_roh_counts_fitted | samples      |     4 |          1.126 |        266.082 |         2.304 |       264.363 |      2.045 |     0.994 |
| duckhts_roh_counts_fitted | sites        |     1 |          0.170 |        127.102 |         0.318 |       127.996 |      1.872 |     1.007 |
| duckhts_roh_counts_fitted | sites        |     2 |          0.348 |        192.742 |         0.644 |       193.160 |      1.850 |     1.002 |
| duckhts_roh_counts_fitted | sites        |     4 |          0.688 |        338.613 |         1.286 |       339.484 |      1.870 |     1.003 |
