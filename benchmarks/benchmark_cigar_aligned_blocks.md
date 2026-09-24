CIGAR aligned blocks benchmark
================

This report compares `cigar_aligned_blocks(CIGAR, POS)` with a plain binary-CIGAR
scan and a SQL block derivation using `list_slice` prefix sums inside lambdas.
The plain scan is the denominator; the SQL derivation is quadratic in the op
count of each alignment record. All three are timed on the same host in one render,
through one in-memory DuckDB connection, over two real alignments staged
through the duckhtsbench registry:

- the `ont-ecoli-k12` workload: ENA run `ERR14686255`, 25,950 MinION reads of
  an *Escherichia coli* K-12 MG1655 derivative, aligned by `minimap2 -x map-ont`
  to the NCBI RefSeq assembly `GCF_000005845.2` and sorted by `samtools` at
  staging time. The genome is one contig, so `read_bam` scans it on one thread.
- the `riker-wgs` workload: the 1000 Genomes high-coverage sample HG00188
  (`ERR3240174`) transcoded from CRAM to BAM at staging time, read through the
  `chr22` region of its index. A region scan is one `read_bam` thread by
  construction, and it keeps the short-read condition to a few million reads.

A measured render requires exact, NULL-aware equality of the three block lists
for every physical record before timing, including duplicate read identities
and missing CIGARs. Timed passes check exact counts and two fingerprints keyed
by physical ordinal, including explicit NULL markers. These checks complement
the fixture oracle in `test/sql/cigar_aligned_blocks.test`.

The recorded measurements at `1139f23ddde1` have aggregate-XOR validation only,
not per-record conformance evidence. XOR can conceal duplicate errors and does
not distinguish NULL lists from empty lists in this DuckDB runtime.
[The recorded data](https://github.com/ryandward/duckhts/blob/54995e6f3a211786b3421a629edfd87d8929e371/benchmarks/data/cigar_aligned_blocks.csv)
retain their measured timings and denominators. The recorded revision below
identifies the measurement shown; rendering stored data does not validate it.

The focused current-source comparison is
[`benchmark_cigar_validation.md`](benchmark_cigar_validation.md). The nearest report
on the short-read input is [`benchmark_riker_wgs.md`](benchmark_riker_wgs.md),
which times whole-file `read_bam` scans against Riker; it shares the registered
BAM but not the region or the projection, so it is context rather than a
comparison.

## Reproduction

Stage both inputs once, with network access only in this step. The long-read
staging needs `minimap2` and `samtools`; the short-read staging needs
`samtools` and downloads a 14 GB CRAM plus the GRCh38 reference:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_ont_ecoli()'
Rscript -e 'duckhtsbench::duckhts_bench_stage_riker()'
```

`DUCKHTS_EXTENSION` selects the extension to measure (default:
`build/release/duckhts.duckdb_extension`). Render from a tracked-clean revision;
rendering resolves the staged paths through the registry and does not download:

``` sh
taskset -c 8-11 Rscript -e 'rmarkdown::render("benchmarks/benchmark_cigar_aligned_blocks.Rmd")'
```

Both inputs run on one `read_bam` thread, so the benchmark pins each pass to
the highest logical CPU in the affinity mask of the rendering process; wrap the
render in `taskset -c` to choose it. Input staging, extension loading,
connection setup, warm-up and result verification are excluded from timing.
To render the stored measurements without claiming a new run, use
`rmarkdown::render("benchmarks/benchmark_cigar_aligned_blocks.Rmd", params = list(measure = FALSE))`.
This mode reads the checked-in results and metadata; it does not run validation.

| Property                                     | Recorded value                                                                                                             |
|:---------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------|
| revision                                     | 1139f23ddde1                                                                                                               |
| src tree                                     | f51e800bb6dc98f1513e5117f0db2c015ce35e00                                                                                   |
| extension MD5                                | 2201f56c1dbd4da138407022d61cd6a8                                                                                           |
| run date                                     | 2026-09-22                                                                                                                 |
| long-read input                              | duckhtsbench registry ont-ecoli-k12: ENA ERR14686255 reads on NCBI RefSeq GCF_000005845.2, minimap2 map-ont, samtools sort |
| long-read reads (ERR14686255)                | 27,377                                                                                                                     |
| long-read ops                                | 2,817,619                                                                                                                  |
| long-read blocks                             | 1,408,917                                                                                                                  |
| long-read aligner                            | minimap2 2.28-r1209                                                                                                        |
| long-read sorter                             | samtools 1.21                                                                                                              |
| short-read input                             | duckhtsbench registry riker-wgs: 1000 Genomes HG00188 ERR3240174 CRAM transcoded to BAM, region chr22                      |
| short-read reads (chr22)                     | 9,569,553                                                                                                                  |
| short-read ops                               | 10,964,914                                                                                                                 |
| short-read blocks                            | 9,971,254                                                                                                                  |
| block fingerprint, long-read (SQL = scalar)  | EC0D17B9470FF589                                                                                                           |
| block fingerprint, short-read (SQL = scalar) | 944928B6155FAD95                                                                                                           |
| DuckDB version                               | v1.5.5                                                                                                                     |
| htslib version                               | 1.24                                                                                                                       |
| CPU                                          | Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz                                                                                   |

## Results

Each workload scans the same rows and returns one row to R, so the timing
covers reading, decoding and projection rather than transporting all blocks.
The Reads column counts physical alignment records, including duplicate names.
Measured runs first compare every physical record against the SQL oracle,
then check exact record/op/block totals and record-keyed XOR and sum fingerprints
on every timed pass. NULL markers distinguish missing CIGARs from valid CIGARs
with no aligned blocks. Recorded-data rendering does not execute those checks;
the recorded `1139f23ddde1` data have the aggregate-XOR limits described above.

| Input                               | Workload             | Threads | CPU affinity | Runs |     Reads |  CIGAR ops |    Blocks | Minimum seconds | Median seconds | Maximum seconds | Median reads/s |
|:------------------------------------|:---------------------|--------:|-------------:|-----:|----------:|-----------:|----------:|----------------:|---------------:|----------------:|---------------:|
| ONT E. coli K-12 (ERR14686255)      | read_bam only        |       1 |           11 |    9 |    27,377 |  2,817,619 | 1,408,917 |           0.278 |          0.283 |           0.302 |         96,739 |
| ONT E. coli K-12 (ERR14686255)      | blocks in SQL        |       1 |           11 |    9 |    27,377 |  2,817,619 | 1,408,917 |          29.959 |         30.415 |          32.293 |            900 |
| ONT E. coli K-12 (ERR14686255)      | cigar_aligned_blocks |       1 |           11 |    9 |    27,377 |  2,817,619 | 1,408,917 |           0.284 |          0.295 |           0.298 |         92,803 |
| Illumina HG00188 chr22 (ERR3240174) | read_bam only        |       1 |           11 |    5 | 9,569,553 | 10,964,914 | 9,971,254 |           4.518 |          4.550 |           4.632 |      2,103,198 |
| Illumina HG00188 chr22 (ERR3240174) | blocks in SQL        |       1 |           11 |    5 | 9,569,553 | 10,964,914 | 9,971,254 |           8.573 |          8.591 |           8.719 |      1,113,904 |
| Illumina HG00188 chr22 (ERR3240174) | cigar_aligned_blocks |       1 |           11 |    5 | 9,569,553 | 10,964,914 | 9,971,254 |           4.981 |          5.001 |           5.914 |      1,913,528 |

| Input                               | Plain scan | Blocks in SQL | cigar_aligned_blocks | SQL cost over plain | Scalar cost over plain |
|:------------------------------------|-----------:|--------------:|---------------------:|--------------------:|-----------------------:|
| ONT E. coli K-12 (ERR14686255)      |      0.283 |        30.415 |                0.295 |              30.132 |                  0.012 |
| Illumina HG00188 chr22 (ERR3240174) |      4.550 |         8.591 |                5.001 |               4.041 |                  0.451 |

The cost of each derivation is its median minus the plain scan’s median on the
same input, which is the time the projection itself adds. The SQL derivation
grows with the square of each read’s op count, so its cost is set by the
long-read input’s op distribution; the scalar’s cost is linear in ops and is
bounded by the writes of three lists per read.

The measurement does not include staging, extension loading, or any consumer
of the block lists, and it does not measure memory. Both inputs run on one
`read_bam` thread, the long-read input because its genome is one contig and the
short-read input because a region scan is single-threaded, so the numbers say
nothing about scaling across contigs. Small timings vary with scheduling; the
differences are for this host and these inputs only.
