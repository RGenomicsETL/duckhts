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

Source revision: d39b5f7606946ddcd98c4244e3e647ab2eff8b54; `src` tree d8d9a75ea2378d58e4524438db2e204272c6b76f. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: e983edbbe3c312f6e8eac9fef495a0bb060cfa711bf2f9e1c2b2d6f885de48d3. Host load before the run: 1.25, 1.20, 0.89. CPU affinity: 8-15. Empty-run overhead: 39.6 MiB. Every observation is in [`benchmark_roh_scaling_runs.csv`](benchmark_roh_scaling_runs.csv).

| macro                         | dimension    | scale | threads | samples | input_rows | runs_out | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib | within_budget |
|:------------------------------|:-------------|------:|--------:|--------:|-----------:|:---------|---------------:|------------:|------------:|---------------:|------------:|-----------:|:--------------|
| duckhts_roh_af_table          | samples      |     1 |       1 |       8 |    6793144 | 1553     |          2.635 |       2.562 |       3.001 |           74.0 |        74.1 |      130.5 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       1 |      16 |   13586288 | 3112     |          4.156 |       3.845 |       4.774 |           90.1 |        90.2 |      170.3 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.478 |       6.432 |       6.718 |          128.3 |       128.8 |      250.0 | TRUE          |
| duckhts_roh_af_table          | samples      |     1 |       4 |       8 |    6793144 | 1553     |          2.656 |       2.580 |       2.746 |           80.1 |        81.6 |      164.4 | TRUE          |
| duckhts_roh_af_table          | samples      |     2 |       4 |      16 |   13586288 | 3112     |          3.863 |       3.845 |       4.648 |          113.7 |       118.9 |      204.2 | TRUE          |
| duckhts_roh_af_table          | samples      |     4 |       4 |      32 |   27172576 | 6037     |          7.161 |       6.323 |       7.541 |          168.6 |       169.1 |      283.8 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       1 |       8 |    6793144 | 1553     |          4.205 |       3.575 |       4.737 |          111.5 |       113.0 |      449.7 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       1 |      16 |   13586288 | 3112     |          4.682 |       4.677 |       5.785 |          139.0 |       139.0 |      489.5 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       1 |      32 |   27172576 | 6037     |          6.681 |       6.669 |       7.160 |          177.6 |       182.0 |      569.2 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     1 |       4 |       8 |    6793144 | 1553     |          3.737 |       3.582 |       3.801 |          143.2 |       147.0 |      483.6 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     2 |       4 |      16 |   13586288 | 3112     |          4.634 |       4.465 |       5.721 |          177.1 |       177.1 |      523.4 | TRUE          |
| duckhts_roh_af_table_bcftools | samples      |     4 |       4 |      32 |   27172576 | 6037     |          6.636 |       6.467 |       7.246 |          220.9 |       255.6 |      603.0 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       1 |       8 |    6793144 | 356      |          3.843 |       3.770 |       4.087 |          243.5 |       244.4 |      470.4 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       1 |      16 |   13586288 | 1256     |          5.167 |       5.140 |       6.073 |          383.2 |       384.2 |      610.5 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       1 |      32 |   27172576 | 5120     |          8.650 |       8.501 |       8.668 |          696.5 |       699.8 |      951.8 | TRUE          |
| duckhts_roh_ancestry          | joint        |     1 |       4 |       8 |    6793144 | 356      |          3.750 |       3.739 |       3.789 |          275.9 |       279.2 |      479.0 | TRUE          |
| duckhts_roh_ancestry          | joint        |     2 |       4 |      16 |   13586288 | 1256     |          5.252 |       5.245 |       5.727 |          487.9 |       497.2 |      627.8 | TRUE          |
| duckhts_roh_ancestry          | joint        |     4 |       4 |      32 |   27172576 | 5120     |          8.520 |       8.500 |       8.831 |          789.1 |       793.0 |      986.3 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       1 |       8 |    6793144 | 356      |          3.843 |       3.770 |       4.087 |          243.5 |       244.4 |      470.4 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       1 |       8 |    6793144 | 641      |          4.040 |       4.010 |       4.356 |          373.9 |       374.7 |      590.2 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       1 |       8 |    6793144 | 1337     |          4.701 |       4.578 |       4.838 |          632.2 |       633.1 |      829.8 | TRUE          |
| duckhts_roh_ancestry          | reference    |     1 |       4 |       8 |    6793144 | 356      |          3.750 |       3.739 |       3.789 |          275.9 |       279.2 |      479.0 | TRUE          |
| duckhts_roh_ancestry          | reference    |     2 |       4 |       8 |    6793144 | 641      |          4.112 |       4.022 |       4.172 |          481.0 |       485.9 |      607.4 | TRUE          |
| duckhts_roh_ancestry          | reference    |     4 |       4 |       8 |    6793144 | 1337     |          4.737 |       4.672 |       5.058 |          782.3 |       820.6 |      864.3 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       1 |       8 |    6793144 | 1337     |          4.701 |       4.578 |       4.838 |          632.2 |       633.1 |      829.8 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       1 |      16 |   13586288 | 2637     |          5.796 |       5.793 |       5.873 |          659.0 |       662.3 |      870.4 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       1 |      32 |   27172576 | 5120     |          8.650 |       8.501 |       8.668 |          696.5 |       699.8 |      951.8 | TRUE          |
| duckhts_roh_ancestry          | samples      |     1 |       4 |       8 |    6793144 | 1337     |          4.737 |       4.672 |       5.058 |          782.3 |       820.6 |      864.3 | TRUE          |
| duckhts_roh_ancestry          | samples      |     2 |       4 |      16 |   13586288 | 2637     |          6.049 |       6.049 |       6.603 |          780.5 |       800.6 |      905.0 | TRUE          |
| duckhts_roh_ancestry          | samples      |     4 |       4 |      32 |   27172576 | 5120     |          8.520 |       8.500 |       8.831 |          789.1 |       793.0 |      986.3 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.019 |       0.018 |       0.019 |           49.8 |        49.9 |       58.3 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.034 |       0.033 |       0.035 |           52.0 |        52.1 |       65.8 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.051 |       0.048 |       0.052 |           53.7 |        53.9 |       73.3 | TRUE          |
| duckhts_roh_counts            | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.026 |       0.025 |       0.027 |           51.3 |        52.2 |       92.2 | TRUE          |
| duckhts_roh_counts            | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.033 |       0.028 |       0.035 |           58.8 |        62.1 |       99.6 | TRUE          |
| duckhts_roh_counts            | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.036 |       0.034 |       0.046 |           63.7 |        67.0 |      107.1 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.268 |       0.261 |       0.305 |           96.9 |        97.1 |      202.8 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       1 |      16 |    4000000 | 400      |          0.521 |       0.520 |       0.539 |          143.9 |       144.1 |      340.2 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       1 |      32 |    8000000 | 800      |          1.032 |       1.015 |       1.041 |          237.6 |       237.6 |      614.8 | TRUE          |
| duckhts_roh_counts            | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.166 |       0.157 |       0.183 |          133.5 |       141.5 |      280.6 | TRUE          |
| duckhts_roh_counts            | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.370 |       0.312 |       0.377 |          222.5 |       234.5 |      418.0 | TRUE          |
| duckhts_roh_counts            | samples      |     4 |       4 |      32 |    8000000 | 800      |          0.580 |       0.549 |       0.603 |          325.7 |       334.0 |      692.6 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.156 |       0.155 |       0.157 |          122.5 |       122.9 |      212.0 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.323 |       0.318 |       0.363 |          187.7 |       187.9 |      384.3 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       1 |       1 |    4000000 | 400      |          0.670 |       0.655 |       0.687 |          334.6 |       334.9 |      729.1 | TRUE          |
| duckhts_roh_counts            | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.158 |       0.157 |       0.164 |          129.7 |       130.0 |      523.1 | TRUE          |
| duckhts_roh_counts            | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.328 |       0.325 |       0.339 |          196.6 |       196.7 |     1006.6 | TRUE          |
| duckhts_roh_counts            | sites        |     4 |       4 |       1 |    4000000 | 400      |          0.651 |       0.649 |       0.659 |          336.6 |       336.8 |     1973.6 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       1 |       1 |     108757 | 154      |          0.040 |       0.039 |       0.042 |           56.7 |        57.0 |       58.3 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       1 |       2 |     217514 | 387      |          0.074 |       0.071 |       0.076 |           59.2 |        59.4 |       65.8 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       1 |       3 |     326271 | 602      |          0.110 |       0.103 |       0.116 |           62.3 |        62.4 |       73.3 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     1 |       4 |       1 |     108757 | 154      |          0.046 |       0.043 |       0.053 |           58.3 |        58.3 |       92.2 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     2 |       4 |       2 |     217514 | 387      |          0.060 |       0.045 |       0.068 |           68.0 |        69.4 |       99.6 | TRUE          |
| duckhts_roh_counts_fitted     | real_samples |     3 |       4 |       3 |     326271 | 602      |          0.071 |       0.066 |       0.071 |           77.1 |        78.4 |      107.1 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       1 |       8 |    2000000 | 200      |          0.596 |       0.595 |       0.607 |          105.3 |       106.6 |      202.8 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       1 |      16 |    4000000 | 400      |          1.201 |       1.168 |       1.271 |          160.5 |       165.1 |      340.2 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       1 |      32 |    8000000 | 800      |          2.369 |       2.345 |       2.477 |          262.9 |       268.4 |      614.8 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     1 |       4 |       8 |    2000000 | 200      |          0.417 |       0.378 |       0.652 |          155.5 |       161.3 |      280.6 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     2 |       4 |      16 |    4000000 | 400      |          0.789 |       0.688 |       0.941 |          237.4 |       240.0 |      418.0 | TRUE          |
| duckhts_roh_counts_fitted     | samples      |     4 |       4 |      32 |    8000000 | 800      |          1.185 |       1.052 |       1.519 |          369.4 |       383.4 |      692.6 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       1 |       1 |    1000000 | 100      |          0.328 |       0.327 |       0.330 |          128.1 |       128.2 |      212.0 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       1 |       1 |    2000000 | 200      |          0.672 |       0.657 |       0.723 |          193.3 |       193.4 |      384.3 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       1 |       1 |    4000000 | 400      |          1.333 |       1.317 |       1.483 |          339.5 |       339.6 |      729.1 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     1 |       4 |       1 |    1000000 | 100      |          0.318 |       0.313 |       0.321 |          136.8 |       141.0 |      523.1 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     2 |       4 |       1 |    2000000 | 200      |          0.635 |       0.612 |       0.639 |          196.4 |       197.4 |     1006.6 | TRUE          |
| duckhts_roh_counts_fitted     | sites        |     4 |       4 |       1 |    4000000 | 400      |          1.228 |       1.227 |       1.267 |          340.3 |       344.1 |     1973.6 | TRUE          |

Time exponents per doubling (log2 of the median ratio) and peak-RSS ratios to 1×:

| macro                         | dimension | threads | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | one_x_seconds |
|:------------------------------|:----------|--------:|-----------------:|-----------------:|-------------:|-------------:|--------------:|
| duckhts_roh_ancestry          | joint     |       1 |            0.427 |            0.743 |        1.574 |        2.860 |         3.843 |
| duckhts_roh_ancestry          | reference |       1 |            0.072 |            0.219 |        1.536 |        2.596 |         3.843 |
| duckhts_roh_af_table          | samples   |       1 |            0.657 |            0.640 |        1.218 |        1.734 |         2.635 |
| duckhts_roh_af_table_bcftools | samples   |       1 |            0.155 |            0.513 |        1.247 |        1.593 |         4.205 |
| duckhts_roh_ancestry          | samples   |       1 |            0.302 |            0.578 |        1.042 |        1.102 |         4.701 |
| duckhts_roh_counts            | samples   |       1 |            0.957 |            0.987 |        1.485 |        2.452 |         0.268 |
| duckhts_roh_counts_fitted     | samples   |       1 |            1.011 |            0.980 |        1.524 |        2.497 |         0.596 |
| duckhts_roh_counts            | sites     |       1 |            1.047 |            1.053 |        1.532 |        2.731 |         0.156 |
| duckhts_roh_counts_fitted     | sites     |       1 |            1.037 |            0.988 |        1.509 |        2.650 |         0.328 |
| duckhts_roh_ancestry          | joint     |       4 |            0.486 |            0.698 |        1.768 |        2.860 |         3.750 |
| duckhts_roh_ancestry          | reference |       4 |            0.133 |            0.204 |        1.743 |        2.835 |         3.750 |
| duckhts_roh_af_table          | samples   |       4 |            0.540 |            0.891 |        1.419 |        2.105 |         2.656 |
| duckhts_roh_af_table_bcftools | samples   |       4 |            0.310 |            0.518 |        1.237 |        1.543 |         3.737 |
| duckhts_roh_ancestry          | samples   |       4 |            0.353 |            0.494 |        0.998 |        1.009 |         4.737 |
| duckhts_roh_counts            | samples   |       4 |            1.161 |            0.647 |        1.667 |        2.440 |         0.166 |
| duckhts_roh_counts_fitted     | samples   |       4 |            0.920 |            0.586 |        1.527 |        2.376 |         0.417 |
| duckhts_roh_counts            | sites     |       4 |            1.054 |            0.990 |        1.516 |        2.595 |         0.158 |
| duckhts_roh_counts_fitted     | sites     |       4 |            1.000 |            0.951 |        1.436 |        2.488 |         0.318 |

## Findings

Every macro stays within its declared budget at every size and thread count,
synthetic and real, and the render enforces that.

`duckhts_roh_counts`: its one-thread 1× runs are under the 5-second floor, so
it gets a memory verdict, not a timing verdict. On the real workload, one
thread decodes one sample’s 108,757 chr20 sites in
0.019 s at 49.8 MiB
peak RSS and three samples in 0.051 s at
53.7 MiB. The operating scale this report supports
is what it measures: per-chromosome lists of up to 4 million synthetic sites,
up to 32 synthetic samples, and three real 30× samples at 108,757 sites;
larger panels and cohorts are not measured here.

`duckhts_roh_counts_fitted` returns the same runs as `duckhts_roh_counts` in
every cell, which the render checks, and stays within the same budget. Its
peak RSS over the scalar decode at the largest inputs, one thread, median of
three: 4.9 MiB with 4 million sites of one sample and
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
6.478 s at 128.3 MiB
peak RSS (ceiling 250 MiB), and
`duckhts_roh_af_table_bcftools` takes 6.681 s
at 177.6 MiB (ceiling 569.2 MiB).
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
| duckhts_roh_ancestry | joint        |     1 |          6.247 |       1780.504 |         3.843 |       243.484 |      0.615 |     0.137 |
| duckhts_roh_ancestry | joint        |     2 |         10.128 |       3440.414 |         5.167 |       383.191 |      0.510 |     0.111 |
| duckhts_roh_ancestry | joint        |     4 |         18.219 |       6601.320 |         8.650 |       696.512 |      0.475 |     0.106 |
| duckhts_roh_ancestry | reference    |     1 |          6.247 |       1780.504 |         3.843 |       243.484 |      0.615 |     0.137 |
| duckhts_roh_ancestry | reference    |     2 |          6.551 |       2156.871 |         4.040 |       373.934 |      0.617 |     0.173 |
| duckhts_roh_ancestry | reference    |     4 |          7.254 |       2717.016 |         4.701 |       632.223 |      0.648 |     0.233 |
| duckhts_roh_ancestry | samples      |     1 |          7.254 |       2717.016 |         4.701 |       632.223 |      0.648 |     0.233 |
| duckhts_roh_ancestry | samples      |     2 |         10.858 |       4040.809 |         5.796 |       659.043 |      0.534 |     0.163 |
| duckhts_roh_ancestry | samples      |     4 |         18.219 |       6601.320 |         8.650 |       696.512 |      0.475 |     0.106 |
| duckhts_roh_counts   | real_samples |     1 |          0.022 |         54.445 |         0.019 |        49.840 |      0.853 |     0.915 |
| duckhts_roh_counts   | real_samples |     2 |          0.040 |         65.332 |         0.034 |        52.047 |      0.857 |     0.797 |
| duckhts_roh_counts   | real_samples |     3 |          0.060 |         79.270 |         0.051 |        53.676 |      0.847 |     0.677 |
| duckhts_roh_counts   | samples      |     1 |          0.329 |        254.898 |         0.268 |        96.863 |      0.816 |     0.380 |
| duckhts_roh_counts   | samples      |     2 |          0.661 |        466.645 |         0.521 |       143.922 |      0.788 |     0.308 |
| duckhts_roh_counts   | samples      |     4 |          1.337 |        903.488 |         1.032 |       237.562 |      0.772 |     0.263 |
| duckhts_roh_counts   | sites        |     1 |          0.177 |        162.422 |         0.156 |       122.488 |      0.881 |     0.754 |
| duckhts_roh_counts   | sites        |     2 |          0.338 |        310.410 |         0.323 |       187.707 |      0.956 |     0.605 |
| duckhts_roh_counts   | sites        |     4 |          0.702 |        702.176 |         0.670 |       334.625 |      0.954 |     0.477 |

With 32 samples and the full reference, `duckhts_roh_ancestry` went from
6601 MiB to 697 MiB and
from 18.2 s to 8.6 s.
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
| duckhts_roh_ancestry      | joint        |     1 |          3.698 |        244.484 |         3.843 |       243.484 |      1.039 |     0.996 |
| duckhts_roh_ancestry      | joint        |     2 |          5.060 |        385.195 |         5.167 |       383.191 |      1.021 |     0.995 |
| duckhts_roh_ancestry      | joint        |     4 |          8.205 |        698.496 |         8.650 |       696.512 |      1.054 |     0.997 |
| duckhts_roh_ancestry      | reference    |     1 |          3.698 |        244.484 |         3.843 |       243.484 |      1.039 |     0.996 |
| duckhts_roh_ancestry      | reference    |     2 |          3.966 |        372.727 |         4.040 |       373.934 |      1.019 |     1.003 |
| duckhts_roh_ancestry      | reference    |     4 |          4.487 |        629.414 |         4.701 |       632.223 |      1.048 |     1.004 |
| duckhts_roh_ancestry      | samples      |     1 |          4.487 |        629.414 |         4.701 |       632.223 |      1.048 |     1.004 |
| duckhts_roh_ancestry      | samples      |     2 |          5.729 |        647.402 |         5.796 |       659.043 |      1.012 |     1.018 |
| duckhts_roh_ancestry      | samples      |     4 |          8.205 |        698.496 |         8.650 |       696.512 |      1.054 |     0.997 |
| duckhts_roh_counts        | real_samples |     1 |          0.019 |         49.262 |         0.019 |        49.840 |      0.967 |     1.012 |
| duckhts_roh_counts        | real_samples |     2 |          0.033 |         51.418 |         0.034 |        52.047 |      1.045 |     1.012 |
| duckhts_roh_counts        | real_samples |     3 |          0.048 |         52.871 |         0.051 |        53.676 |      1.052 |     1.015 |
| duckhts_roh_counts        | samples      |     1 |          0.253 |         96.043 |         0.268 |        96.863 |      1.061 |     1.009 |
| duckhts_roh_counts        | samples      |     2 |          0.499 |        142.938 |         0.521 |       143.922 |      1.043 |     1.007 |
| duckhts_roh_counts        | samples      |     4 |          0.991 |        236.582 |         1.032 |       237.562 |      1.042 |     1.004 |
| duckhts_roh_counts        | sites        |     1 |          0.150 |        121.785 |         0.156 |       122.488 |      1.043 |     1.006 |
| duckhts_roh_counts        | sites        |     2 |          0.306 |        186.762 |         0.323 |       187.707 |      1.053 |     1.005 |
| duckhts_roh_counts        | sites        |     4 |          0.612 |        333.688 |         0.670 |       334.625 |      1.094 |     1.003 |
| duckhts_roh_counts_fitted | real_samples |     1 |          0.022 |         55.723 |         0.040 |        56.730 |      1.779 |     1.018 |
| duckhts_roh_counts_fitted | real_samples |     2 |          0.038 |         58.586 |         0.074 |        59.227 |      1.936 |     1.011 |
| duckhts_roh_counts_fitted | real_samples |     3 |          0.054 |         61.102 |         0.110 |        62.332 |      2.040 |     1.020 |
| duckhts_roh_counts_fitted | samples      |     1 |          0.291 |        103.973 |         0.596 |       105.316 |      2.049 |     1.013 |
| duckhts_roh_counts_fitted | samples      |     2 |          0.570 |        159.020 |         1.201 |       160.492 |      2.108 |     1.009 |
| duckhts_roh_counts_fitted | samples      |     4 |          1.126 |        266.082 |         2.369 |       262.926 |      2.103 |     0.988 |
| duckhts_roh_counts_fitted | sites        |     1 |          0.170 |        127.102 |         0.328 |       128.074 |      1.931 |     1.008 |
| duckhts_roh_counts_fitted | sites        |     2 |          0.348 |        192.742 |         0.672 |       193.297 |      1.931 |     1.003 |
| duckhts_roh_counts_fitted | sites        |     4 |          0.688 |        338.613 |         1.333 |       339.453 |      1.939 |     1.002 |
