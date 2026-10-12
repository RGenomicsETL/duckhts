Parquet conversion scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for the
COPY statements that the `*_convert_parquet_sql` builders return. Issue
<https://github.com/RGenomicsETL/duckhts/issues/298> lists these executed COPY
workloads for BAM, VCF and GFF as not measured. A builder only returns SQL;
the caller executes it, as the R wrappers do through DBI. Here each
observation builds the statement, untimed, and executes it in a fresh DuckDB
CLI process. No earlier report measures these workloads; this is their first
baseline.

## Workloads

Every input is a registry artifact of `r/duckhtsbench`. Before anything is
measured, each is checked against the identity the registry records: its
checksum and size, or its counted content. The ONT alignment records neither,
so the report gives its SHA-256. The builders run with their defaults, plus
the region or column list that cuts the input and, for the CRAM, the
reference.

- **BAM, short reads** (`duckhts_bam_convert_parquet_sql`, which reads SAM,
  BAM and CRAM): the chr1 CRAM slice of the 1000 Genomes 30× sample HG00403
  (`bam_mismatch_hg00403_chr1_cram`, 150-base reads), with the registered
  GRCh38 FASTA as the reference (`liftover_grch38_fasta`). The 1×, 2× and 4×
  regions start at the beginning of chr1 and end at the 12.5%, 25% and 50%
  quantiles of the alignment start positions, so each holds twice the
  alignments of the one before.
- **BAM, long reads**: the minimap2 alignment of an ONT E. coli K-12 run
  (`ont_ecoli_k12_bam`). Its regions end at the 25% and 50% quantiles of the
  mapped start positions and at the end of the chromosome. The whole file has
  fewer reads than one row group holds, so this series measures wide rows
  that never fill a row group. Its runs take well under 5 s, so it gets a
  memory verdict only.
- **VCF, records** (`duckhts_bcf_convert_parquet_sql`): the GIAB HG002 v4.2.1
  phased benchmark VCF on GRCh38 (`geno_giab_phased_source`, indexed by
  `geno_giab_phased_source_tbi`), one sample on chr1 to chr22. The 1× and 2×
  regions list the chromosomes in header order up to the record at the 25%
  and 50% rank; the 4× region lists all of them.
- **VCF, samples**: the 1000 Genomes phase 3 chr22 cut of the SQL scaling
  audit (`sql_scaling_phase3_1x`, biallelic SNVs of 2,504 samples). Its records
  are fixed. `columns :=` keeps the record columns and the GT columns of the
  first 626, 1,252 and 2,504 samples in header order. In the builders’ layout
  every sample is one column, so the samples set the width of a row.
- **GFF** (`duckhts_gff_convert_parquet_sql`): the GENCODE v49 basic
  annotation (`gencode_v49_basic_gff3`). The 1× input is the file itself; the
  2× and 4× inputs are two and four copies of it, concatenated at render time.
  Each copy is a whole gzip member, and the reader decompresses the members in
  sequence.

The cuts, the concatenated copies and the input denominators (rows and
decoded bytes) are computed at render time, untimed. Those passes read every
input in full, so every measured run reads from a warm page cache.

## Observations

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`, which
gives the peak RSS. The JSON profiler gives the latency of the COPY statement,
its CPU time, the time of the reader and of the COPY operator, and DuckDB’s
peak buffer and temporary bytes. `max_temp_directory_size` is zero, so a spill
fails the render; the memory limit is DuckDB’s default, which the report
records. Each run writes its output into a temporary directory and deletes it
once its size, rows, row groups and DuckHTS metadata are read. The render
stops if less than 10 GiB is free on that file system before a run. Each cell
has three repetitions at one thread and at four threads, run in an order that
alternates sizes. The render command pins the processes to four cores and
runs them at scheduling priority −15. The builder query that returns the
statement is run before each observation and is not timed; it reads only the
header.

## Budget, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the extension and run no COPY.

The ceiling of each observation is overhead + 3 × the live state below, and it
is a gate: the render stops if any observation exceeds its ceiling or spills.
The live state follows the reader contracts and the Parquet writer of DuckDB
v1.5.1 (`7dbb2e646f`, the CLI used here):

- **One writer.** The generated COPY keeps insertion order, and DuckHTS’s
  table functions, registered through the C API, report no batch index:
  `src/execution/operator/scan/physical_table_scan.cpp:255` needs
  `get_partition_data`, which `src/main/capi/table_function-c.cpp` never sets.
  `ParquetWriteExecutionMode` then picks `REGULAR_COPY_TO_FILE`
  (`extension/parquet/parquet_extension.cpp:661-668`), whose sink is not
  parallel
  (`src/include/duckdb/execution/operator/persistent/physical_copy_to_file.hpp:88-89`).
  One thread scans and buffers, whatever the `threads` setting.
- **The row group.** The sink buffers decoded rows until it holds
  `row_group_size` rows, 100,000 by default in the builders, and checks this
  after each chunk of 2,048 rows (`parquet_extension.cpp:371-390`). A byte limit
  applies only when `ROW_GROUP_SIZE_BYTES` is given, and DuckDB refuses that
  option while insertion order is kept (`parquet_extension.cpp:333-337`). So a
  row group holds R = min(rows, 102,047) rows of w bytes. w is the decoded
  width of a row in DuckDB’s vector layout: 16 bytes for each string plus its
  bytes beyond 12, 16 bytes for each list plus its elements, and the fixed
  width of every other value. Validity bits, under 1%, are left out. w is
  counted from each input before measuring.
- **Its encoded pages.** Each page is compressed into a buffer of ZSTD’s bound
  of its plain-encoded size (`extension/parquet/column_writer.cpp:156-163`), and
  the writer keeps all of them until the row group is written
  (`extension/parquet/parquet_writer.cpp:512-599` and `637-690`). Plain encoding
  is never larger than the decoded layout for these column types, so this is
  at most R × w. The uncompressed buffers of the pages being encoded are freed
  page by page; the factor of 3 covers them.
- **Column dictionaries.** The state of each leaf column holds a dictionary of
  NextPowerOfTwo(2 × ⌊R/5⌋) slots of at most 24 bytes and a target buffer of at
  most 1 MiB until the row group is written
  (`extension/parquet/include/writer/templated_column_writer.hpp:66-75`,
  `src/include/duckdb/common/primitive_dictionary.hpp:37-62`): at most 2.5 MiB
  per column when R is 102,047.
- **The reader.** One output chunk of 2,048 rows, 2,048 × w bytes, and the
  state of htslib 1.24 (`third_party/htslib`):
  - for a CRAM read with the builders’ default of two decompression threads,
    htslib shares one copy of the reference between the threads and loads the
    whole contig (`cram/cram_io.c:3253` and `5573-5585`): 248,956,422 bytes for
    chr1. It decodes slices ahead into a queue of twice the threads and reads
    one slice past it (`cram/cram_decode.c:3517-3523`), so up to 7 slices of at
    most 10,000 alignments (`cram/cram_structs.h:87`) are held. For each
    decoded alignment, with its name, bases, qualities and tags, 1 KiB is
    allowed;
  - for BAM, VCF and GFF, 8 MiB is allowed for the file buffers (64 KiB BGZF or
    gzip blocks, and at most twice the threads in queued decompression jobs:
    `bgzf.c:1756-1758`), the header and one record.

That is 2 × R × w + leaf columns × the dictionary + 2,048 × w + the htslib
state.

The render also stops if any output’s row count differs from its input’s or
lacks the DuckHTS reader metadata. For every series whose 1× run takes at
least 5 s on one thread, which is the least that `STYLE.md` accepts for a
timing verdict, it stops if a time exponent is above 1.25 at either doubling.
A series under 5 s gets the memory verdict alone.

``` sh
nice -n -15 taskset -c 16-19 Rscript -e 'rmarkdown::render("benchmarks/benchmark_convert_parquet_scaling.Rmd")'
```

`CONVERT_PARQUET_BENCH_SIZE=reduced` runs the same code on small inputs, to
check the report end to end: CRAM regions of 0.2% to 0.8% of the alignments,
VCF regions of 1/64 to 1/16 of the records, 8 to 32 samples, and copies of
the MANE annotation (`mane_v15_ensembl_gff3`) instead of GENCODE. That run
writes its observations to a temporary directory, so pass `output_dir` to keep
its report out of `benchmarks/` too; its numbers are not evidence.

``` sh
CONVERT_PARQUET_BENCH_SIZE=reduced taskset -c 0-3 Rscript -e 'rmarkdown::render("benchmarks/benchmark_convert_parquet_scaling.Rmd", output_dir = "/tmp")'
```

Source revision: 68dfa0332ad9ea34eea556eddea231cd83d842f6; `src` tree 2de704e472ceeda276caa0eeef7093221378f3c4. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 66598566776aa44faffb68a8e0ba62961346c11db160962d115c30291d33586e. Host load at start: 2.31, 2.24, 1.89. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0. Every observation is in [`benchmark_convert_parquet_scaling_runs.csv`](benchmark_convert_parquet_scaling_runs.csv).

Overhead: 39 MiB.

## Inputs

The CRAM slice has 62,660,978 alignments in 6,270 slices; the ONT file has 27,377 reads, 25,689 of them mapped, and SHA-256 8bda421bd0a523b2b76ab3d9e984b9ba665c5ec799482378ef02965317c019ef. `rows` are the rows each builder’s SELECT returns, counted untimed; `decoded_mib` and `bytes_per_row` are their decoded size; `row_group_rows` is R, `leaf_columns` the columns of the output, and `live_state_mib` and `budget_mib` follow from the model above.

| series          | scale | input                                                                                                                      |       rows | leaf_columns | decoded_mib | bytes_per_row | row_group_rows | live_state_mib | budget_mib |
|:----------------|------:|:---------------------------------------------------------------------------------------------------------------------------|-----------:|-------------:|------------:|--------------:|---------------:|---------------:|-----------:|
| bam_short_reads |     1 | chr1:1-29661267                                                                                                            |  7,832,623 |           14 |     4,059.9 |         543.5 |        102,047 |          447.6 |    1,381.9 |
| bam_short_reads |     2 | chr1:1-61075319                                                                                                            | 15,665,245 |           14 |     8,119.6 |         543.5 |        102,047 |          447.6 |    1,381.9 |
| bam_short_reads |     4 | chr1:1-123524417                                                                                                           | 31,330,489 |           14 |    16,239.1 |         543.5 |        102,047 |          447.6 |    1,381.9 |
| bam_long_reads  |     1 | NC_000913.3:1-1175904                                                                                                      |      6,423 |           14 |        30.5 |       4,974.3 |          6,423 |           94.0 |      320.9 |
| bam_long_reads  |     2 | NC_000913.3:1-2502986                                                                                                      |     12,845 |           14 |        61.3 |       5,005.2 |         12,845 |          157.0 |      510.1 |
| bam_long_reads  |     4 | NC_000913.3:1-4641652                                                                                                      |     25,689 |           14 |       121.1 |       4,941.7 |         25,689 |          279.0 |      876.1 |
| vcf_records     |     1 | chr1,chr2,chr3,chr4:1-47831439                                                                                             |  1,012,086 |           31 |     1,262.3 |       1,307.8 |        102,047 |          342.6 |    1,066.8 |
| vcf_records     |     2 | chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8:1-10303467                                                                         |  2,024,171 |           31 |     2,500.9 |       1,295.5 |        102,047 |          340.2 |    1,059.6 |
| vcf_records     |     4 | chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 |  4,048,342 |           31 |     5,019.9 |       1,300.2 |        102,047 |          341.1 |    1,062.4 |
| vcf_samples     |     1 | 626 samples                                                                                                                |    124,658 |          660 |     1,247.2 |      10,491.0 |        102,047 |        3,720.5 |   11,200.4 |
| vcf_samples     |     2 | 1,252 samples                                                                                                              |    124,658 |        1,286 |     2,437.9 |      20,507.0 |        102,047 |        7,254.5 |   21,802.6 |
| vcf_samples     |     4 | 2,504 samples                                                                                                              |    124,658 |        2,538 |     4,819.4 |      40,539.0 |        102,047 |       14,322.7 |   43,007.0 |
| gff_records     |     1 | 1 copy                                                                                                                     |  5,866,158 |            9 |     2,876.4 |         514.2 |        102,047 |          131.6 |      433.8 |
| gff_records     |     2 | 2 copies                                                                                                                   | 11,732,316 |            9 |     5,752.7 |         514.2 |        102,047 |          131.6 |      433.8 |
| gff_records     |     4 | 4 copies                                                                                                                   | 23,464,632 |            9 |    11,505.5 |         514.2 |        102,047 |          131.6 |      433.8 |

## Measurements

`seconds` are the latency of the COPY statement over three repetitions. `scan_share` is the share of it spent in the reader, the rest being the Parquet writer. `cpu_per_second` is DuckDB’s CPU time per second of latency; htslib’s decompression threads are not in it. `row_groups` and `output_mib` describe the Parquet file.

| series          | scale | threads |       rows | row_groups | output_mib | median_seconds | min_seconds | max_seconds | scan_share | cpu_per_second | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:----------------|------:|--------:|-----------:|-----------:|-----------:|---------------:|------------:|------------:|-----------:|---------------:|---------------:|------------:|---------------:|-----------:|
| bam_short_reads |     1 |       1 |  7,832,623 |         79 |      245.6 |         16.876 |      16.135 |      16.887 |       0.40 |           0.99 |          604.9 |       610.5 |          142.8 |    1,381.9 |
| bam_short_reads |     2 |       1 | 15,665,245 |        157 |      482.7 |         32.008 |      31.941 |      33.104 |       0.42 |           1.00 |          651.3 |       654.0 |          142.8 |    1,381.9 |
| bam_short_reads |     4 |       1 | 31,330,489 |        313 |      951.5 |         64.796 |      64.655 |      65.677 |       0.42 |           1.00 |          729.7 |       731.0 |          142.8 |    1,381.9 |
| bam_short_reads |     1 |       4 |  7,832,623 |         79 |      245.6 |         16.127 |      16.088 |      16.134 |       0.42 |           0.99 |          610.6 |       615.4 |          142.9 |    1,381.9 |
| bam_short_reads |     2 |       4 | 15,665,245 |        157 |      482.7 |         32.318 |      32.005 |      32.355 |       0.41 |           1.00 |          654.2 |       659.2 |          142.9 |    1,381.9 |
| bam_short_reads |     4 |       4 | 31,330,489 |        313 |      951.5 |         65.379 |      64.196 |      65.667 |       0.41 |           1.00 |          732.2 |       743.3 |          142.8 |    1,381.9 |
| bam_long_reads  |     1 |       1 |      6,423 |          1 |       11.5 |          0.252 |       0.250 |       0.253 |       0.16 |           0.20 |          112.4 |       112.4 |           64.5 |      320.9 |
| bam_long_reads  |     2 |       1 |     12,845 |          1 |       23.1 |          0.500 |       0.497 |       0.503 |       0.15 |           0.20 |          178.8 |       179.0 |          133.6 |      510.1 |
| bam_long_reads  |     4 |       1 |     25,689 |          1 |       45.5 |          0.950 |       0.940 |       0.982 |       0.15 |           0.20 |          293.5 |       293.8 |          273.2 |      876.1 |
| bam_long_reads  |     1 |       4 |      6,423 |          1 |       11.5 |          0.248 |       0.248 |       0.249 |       0.16 |           0.21 |          115.3 |       115.4 |           64.5 |      320.9 |
| bam_long_reads  |     2 |       4 |     12,845 |          1 |       23.1 |          0.487 |       0.485 |       0.488 |       0.15 |           0.20 |          183.5 |       183.7 |          133.6 |      510.1 |
| bam_long_reads  |     4 |       4 |     25,689 |          1 |       45.5 |          0.963 |       0.954 |       0.990 |       0.14 |           0.20 |          299.5 |       299.7 |          273.3 |      876.1 |
| vcf_records     |     1 |       1 |  1,012,086 |         11 |       18.1 |          8.784 |       8.779 |       8.789 |       0.68 |           0.99 |          316.8 |       316.9 |          296.0 |    1,066.8 |
| vcf_records     |     2 |       1 |  2,024,171 |         21 |       35.7 |         17.228 |      17.214 |      17.236 |       0.69 |           0.99 |          314.2 |       314.3 |          296.0 |    1,059.6 |
| vcf_records     |     4 |       1 |  4,048,342 |         41 |       72.0 |         34.460 |      34.382 |      34.505 |       0.69 |           0.99 |          323.6 |       324.0 |          296.0 |    1,062.4 |
| vcf_records     |     1 |       4 |  1,012,086 |         11 |       18.1 |          8.817 |       8.802 |       8.878 |       0.68 |           0.99 |          322.4 |       323.7 |          296.2 |    1,066.8 |
| vcf_records     |     2 |       4 |  2,024,171 |         21 |       35.7 |         17.299 |      17.068 |      17.304 |       0.69 |           0.99 |          320.8 |       328.2 |          296.2 |    1,059.6 |
| vcf_records     |     4 |       4 |  4,048,342 |         41 |       72.0 |         34.520 |      34.459 |      34.600 |       0.69 |           0.99 |          349.2 |       349.7 |          296.4 |    1,062.4 |
| vcf_samples     |     1 |       1 |    124,658 |          2 |        7.1 |         15.903 |      15.903 |      15.943 |       0.72 |           0.95 |        2,265.4 |     2,268.5 |        1,859.7 |   11,200.4 |
| vcf_samples     |     2 |       1 |    124,658 |          2 |       13.5 |         24.751 |      24.543 |      24.986 |       0.60 |           0.93 |        4,301.1 |     4,303.3 |        3,630.5 |   21,802.6 |
| vcf_samples     |     4 |       1 |    124,658 |          2 |       25.8 |         46.109 |      45.962 |      47.204 |       0.48 |           0.91 |        8,817.5 |     8,817.5 |        7,174.0 |   43,007.0 |
| vcf_samples     |     1 |       4 |    124,658 |          2 |        7.1 |         15.953 |      15.909 |      16.095 |       0.71 |           0.95 |        2,481.9 |     2,481.9 |        1,859.7 |   11,200.4 |
| vcf_samples     |     2 |       4 |    124,658 |          2 |       13.5 |         24.661 |      24.651 |      24.718 |       0.59 |           0.93 |        4,537.2 |     4,537.3 |        3,630.7 |   21,802.6 |
| vcf_samples     |     4 |       4 |    124,658 |          2 |       25.8 |         46.425 |      46.238 |      46.929 |       0.47 |           0.91 |        9,618.7 |     9,636.2 |        7,174.2 |   43,007.0 |
| gff_records     |     1 |       1 |  5,866,158 |         59 |       84.9 |          9.862 |       9.835 |       9.864 |       0.49 |           1.00 |          176.6 |       176.7 |           90.9 |      433.8 |
| gff_records     |     2 |       1 | 11,732,316 |        117 |      169.8 |         19.785 |      19.762 |      19.851 |       0.49 |           0.99 |          182.4 |       182.4 |           90.9 |      433.8 |
| gff_records     |     4 |       1 | 23,464,632 |        234 |      339.6 |         39.368 |      39.327 |      39.371 |       0.48 |           1.00 |          186.6 |       186.8 |           90.9 |      433.8 |
| gff_records     |     1 |       4 |  5,866,158 |         59 |       84.9 |          9.846 |       9.835 |       9.928 |       0.48 |           1.00 |          180.0 |       183.6 |           96.8 |      433.8 |
| gff_records     |     2 |       4 | 11,732,316 |        117 |      169.8 |         19.673 |      19.671 |      19.767 |       0.48 |           0.99 |          181.7 |       183.0 |           96.8 |      433.8 |
| gff_records     |     4 |       4 | 23,464,632 |        234 |      339.6 |         39.390 |      39.341 |      40.256 |       0.48 |           1.00 |          187.6 |       190.3 |           96.8 |      433.8 |

Growth of the median time against the input rows, or against the samples of the `vcf_samples` series (an exponent of 1 is linear), the size ratio at each doubling, and peak-RSS ratios to 1×:

| series          | threads | one_x_seconds | ratio_2x | ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|:----------------|--------:|--------------:|---------:|---------:|-----------------:|-----------------:|-------------:|-------------:|
| bam_short_reads |       1 |        16.876 |        2 |        2 |            0.923 |            1.017 |        1.077 |        1.206 |
| bam_short_reads |       4 |        16.127 |        2 |        2 |            1.003 |            1.016 |        1.071 |        1.199 |
| bam_long_reads  |       1 |         0.252 |        2 |        2 |            0.991 |            0.926 |        1.591 |        2.611 |
| bam_long_reads  |       4 |         0.248 |        2 |        2 |            0.970 |            0.985 |        1.592 |        2.598 |
| vcf_records     |       1 |         8.784 |        2 |        2 |            0.972 |            1.000 |        0.992 |        1.021 |
| vcf_records     |       4 |         8.817 |        2 |        2 |            0.972 |            0.997 |        0.995 |        1.083 |
| vcf_samples     |       1 |        15.903 |        2 |        2 |            0.638 |            0.898 |        1.899 |        3.892 |
| vcf_samples     |       4 |        15.953 |        2 |        2 |            0.628 |            0.913 |        1.828 |        3.876 |
| gff_records     |       1 |         9.862 |        2 |        2 |            1.004 |            0.993 |        1.033 |        1.057 |
| gff_records     |       4 |         9.846 |        2 |        2 |            0.999 |            1.002 |        1.009 |        1.042 |

## Findings

Every observation is within its ceiling, none spills, and every output has its input’s rows and the DuckHTS reader metadata; the render enforces all of this. Every COPY ran as one `COPY_TO_FILE` operator, the writer that keeps insertion order without batch indexes. DuckDB’s CPU time is at most 1.00 s per second of latency in any run, and four threads change the median time of a cell by a factor of 0.96 to 1.01, so one DuckDB thread does the work whatever the setting.

On one thread the peak RSS above the overhead is between 0.59 (`vcf_samples` at 2×) and 1.54 (`bam_short_reads` at 4×) times the declared live state, against the factor of 3 that the ceiling allows.

In the three records series the writer flushes a row group at 100,000 rows, so the model holds no state that grows with the input, and DuckDB’s peak buffer differs by at most 0.0 MiB between the sizes of a series. From 1× to 4× on one thread the peak RSS changes by +124.8 MiB (a ratio of 1.21) for the CRAM, +6.8 MiB (1.02) for the VCF and +10.0 MiB (1.06) for the GFF; the model does not account for that change, and it is not memory that DuckDB’s buffer manager holds.

The time exponents of the series with a timing verdict (`bam_short_reads`, `vcf_records`, `vcf_samples`, `gff_records`) are at most 1.02. `bam_long_reads` takes less than 5 s at 1× and gets no timing verdict.

On one thread, 31,330,489 CRAM alignments convert in 64.8 s (483,523 alignments per second, 42% of the time in the reader) at 729.7 MiB peak RSS; the model attributes 237 MiB of that to the chr1 reference, which htslib loads whole for its decompression threads. 4,048,342 VCF records take 34.5 s (69% in the reader) at 323.6 MiB, and 23,464,632 GFF features 39.4 s (48% in the reader) at 186.6 MiB.

Memory grows with the width of a row. In the builders’ layout each sample of a VCF is one column, and a row group holds the genotypes of every sample for its records: on one thread 626 samples take 2265.4 MiB peak RSS and 2,504 samples 8817.5 MiB, about 3.5 MiB per sample, for an input file of 21 MiB and an output of 25.8 MiB. Long reads make wide rows too: at 4× the ONT reads make one row group of 25,689 rows, and the peak RSS is 293.5 MiB, 2.61 times the 1× peak.

The operating scale this report supports is what it measures: 31,330,489 alignments of one 30× chromosome, a whole-genome single-sample VCF of 4,048,342 records, four copies of a full human annotation (23,464,632 features) and 2,504 samples of 124,658 records, all with the default row group of 100,000 rows. More samples, wider rows or a larger row group need memory in proportion; `row_group_size` is the builders’ parameter for it, and its effect is not measured here.
