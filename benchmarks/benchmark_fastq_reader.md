DuckHTS FASTQ reader and fused QC benchmark
================

<!-- benchmark_fastq_reader.md is generated from benchmark_fastq_reader.Rmd. Do not edit the .md by hand. -->

Status: current measurement report for the direct FASTQ parser and
`duckhts_fastq_qc(...)`. The measurements are tied to the source
revisions, artifacts, inputs, and one-core protocol below; they are not
claims about other machines, compression libraries, or thread counts.

## Questions

This report answers three separate questions:

1.  Does parsing projected FASTQ fields directly improve
    `read_fastq(...)` over the previous temporary-`bam1_t` path?
2.  Does a bounded C aggregate remove the cost of expressing per-cycle
    QC as one relational row per base?
3.  After that correction, where does the one-core DuckHTS path stand
    against a current fastp build doing QC without output, filtering,
    trimming, or duplication analysis?

The workloads are intentionally kept separate. A `COUNT(*)`, a
materialized sequence/quality scan, a complete per-cycle QC reduction,
and fastp are not the same computation.

## Host and software

| item                                | value                                                                                                     |
|:------------------------------------|:----------------------------------------------------------------------------------------------------------|
| host                                | 13th Gen Intel Core i5-13500; Linux x86-64                                                                |
| reader/SQL-negative-control pinning | CPU 2 (P-core 1), one thread, normal priority; page cache warmed before timing                            |
| paired-QC pinning                   | CPU 2 (P-core 1), one thread, nice -10 for both DuckHTS and fastp; page cache warmed before timing        |
| DuckDB                              | 1.5.1 (7dbb2e646f)                                                                                        |
| htslib                              | 1.24, vendored                                                                                            |
| previous DuckHTS                    | 247e5c8                                                                                                   |
| direct-parser DuckHTS               | a88fe105056c                                                                                              |
| fused-QC artifact                   | worktree-final; sha256 198b6f28cf31576022d2fbe42d1c49298b5eda6b6d573a45eab469a608546221                   |
| fastp                               | 1.3.6-4-gd517536; commit d517536b021bca0916cf33cb456f4e4b8aa24456; -O3; Highway, ISA-L 2.31.0, libdeflate |

The fastp command used `-w 1`, `--disable_adapter_trimming`,
`--disable_quality_filtering`, `--disable_length_filtering`,
`--disable_trim_poly_g`, and `--dont_eval_duplication`, with no read
output. DuckHTS used `PRAGMA threads=1`. Both processes were pinned with
`taskset`. The exact-artifact paired QC runs used the same elevated
scheduler priority because unrelated unpinned host workers otherwise
entered CPU 2; a diagnostic first-file DuckHTS run changed from 32.33
seconds under contention to 25.53 seconds at the recorded priority.
Reader before/after and the SQL negative control retain their original
normal-priority measurements and are not mixed into the paired fastp
ratio.

## Reader-only fixtures

The uncompressed fixture contains 2,000,000 ordinary four-line, 150-base
reads (300,000,000 bases; 646 MB). Each DuckHTS cell is the median of
seven timed runs. The compressed fixture is the first file in the exome
manifest below (20,203,002 reads; 2,040,503,202 bases; 1,914,722,761
compressed bytes), with three timed runs per cell.

| input                  | workload | implementation | median_seconds | speedup_vs_previous |
|:-----------------------|:---------|:---------------|---------------:|--------------------:|
| uncompressed synthetic | count    | previous       |          0.501 |               1.000 |
| uncompressed synthetic | count    | direct         |          0.324 |               1.546 |
| uncompressed synthetic | name     | previous       |          0.517 |               1.000 |
| uncompressed synthetic | name     | direct         |          0.335 |               1.543 |
| uncompressed synthetic | strings  | previous       |          1.363 |               1.000 |
| uncompressed synthetic | strings  | direct         |          0.991 |               1.375 |
| uncompressed synthetic | packed   | previous       |          1.229 |               1.000 |
| uncompressed synthetic | packed   | direct         |          0.960 |               1.280 |
| gzip HG001 exome R1    | strings  | previous       |         27.003 |               1.000 |
| gzip HG001 exome R1    | strings  | direct         |         25.686 |               1.051 |
| gzip HG001 exome R1    | packed   | previous       |         26.670 |               1.000 |
| gzip HG001 exome R1    | packed   | direct         |         27.122 |               0.983 |

The direct parser clearly improves the local uncompressed path: 1.55x
for count/name projection, 1.38x for all text fields, and 1.28x for
packed sequence and quality. On the real gzip file, text projection
improves by 4.9%, while the packed projection is 1.7% slower and must be
treated as statistically flat or a small regression rather than a
universal parser win. Decompression and line-oriented transport dominate
that compressed workload.

## Full HG001 exome QC

The retained manifest is `benchmarks/data/fastq_exome_inputs.csv`: eight
Garvan HG001 exome FASTQ files, four R1/R2 lane pairs, checksum-pinned
to the GIAB/NCBI downloads. Together they contain 168,453,484 reads,
17,013,801,884 bases, and 16,161,610,259 compressed bytes. Each
full-corpus cell is one warmed-cache run; reruns are required before
using small differences as optimization evidence.

The DuckHTS formulations are:

- `per-base SQL`:
  `read_fastq(..., sequence_encoding := 'nt16',   quality_representation := 'phred')`,
  followed by per-base list expansion and SQL grouping;
- `fused aggregate`: projected sequence/quality strings flow directly
  into `duckhts_fastq_qc(...)`, and only its 808 cycle rows are
  expanded;
- `fastp QC`: eight one-thread fastp processes run sequentially, one per
  file; the table reports their summed wall time and largest per-process
  RSS.

| engine  | workload        |  seconds | peak_rss_kib |     reads |       bases | threads | cpu | nice | million_reads_per_second | million_bases_per_second | relative_to_fastp |
|:--------|:----------------|---------:|-------------:|----------:|------------:|--------:|----:|-----:|-------------------------:|-------------------------:|------------------:|
| DuckHTS | per-base SQL    | 1236.185 |       166360 | 168453484 | 17013801884 |       1 |   2 |    0 |                    0.136 |                   13.763 |             0.146 |
| DuckHTS | fused aggregate |  215.149 |       169668 | 168453484 | 17013801884 |       1 |   2 |  -10 |                    0.783 |                   79.079 |             0.839 |
| fastp   | QC, no output   |  180.580 |        45192 | 168453484 | 17013801884 |       1 |   2 |  -10 |                    0.933 |                   94.218 |             1.000 |

The final fused aggregate reaches 782,962 reads/s and 79.08 million
bases/s while returning all global and per-cycle sufficient statistics.
It is **19.1% slower than fastp** on the exact-artifact, equal-priority
one-core gzip comparison. The earlier normal-priority series measured
1,236.185 seconds for the per-base SQL plan and 216.609 seconds for the
aggregate (5.71x); the final aggregate’s 215.149 seconds independently
preserves that conclusion, but the table labels the scheduler difference
instead of presenting the 5.75x quotient as a strictly paired
measurement. The result validates the aggregate architecture, not a
claim that DuckHTS has already won the parser race.

The original reader-only gzip text scan took a 25.686-second median for
the first file; the same source tree plus fused QC took 25.816 seconds
in its smoke run. The final exact-artifact aggregate probe took 25.53
seconds under the paired-QC scheduler protocol. Together these
measurements indicate that gzip decode, line parsing, normalization, and
DuckDB string materialization dominate here. The eleven exact per-cycle
updates are not the observed 19.1% gap to fastp. The next reader
experiment should target a bounded ordinary-four-line block path with
multiline fallback, while retaining htslib transport and the current
error contract.

## Correctness

All eight fused results matched fastp exactly for:

- read and base totals;
- Q20, Q30, and Q40 base totals;
- maximum cycle and mean read length.

Across 808 `(file, read_end, cycle)` rows, DuckHTS integer sufficient
statistics reproduced fastp’s rendered curves within fastp’s output
rounding:

| field_family                        | maximum_absolute_difference |
|:------------------------------------|----------------------------:|
| global integer totals               |                     0.0e+00 |
| global GC fraction                  |                     4.9e-07 |
| per-cycle mean-quality curves       |                     5.0e-05 |
| per-cycle nucleotide-content curves |                     5.0e-07 |

The aggregate stores integer counts and quality sums. Floating means and
fractions are late SQL/report projections, so the integer result is not
coupled to fastp’s JSON precision.

## SIMD scope

The `fastq_qc` logical kernel has scalar, AVX2, ARM NEON, and wasm
SIMD128 registrations. Scalar remains the oracle. Seeded standalone
contracts compare every runnable backend on deterministic offsets/tails
and randomized reads; the NEON and wasm SIMD128 translation units are
also cross-target compiled on this x86 host. Only AVX2 executed in this
report. No ARM or browser performance claim follows from successful
cross-compilation.

SIMD reduces block classification and threshold counting, but exact
per-cycle statistics still update a small scalar state array. More
importantly, the end-to-end result above is input/decompression
dominated. A microkernel multiplier would not justify claiming the same
multiplier for `read_fastq()`.

## Reproduction

`benchmarks/benchmark_fastq_reader.R` owns fixture validation, CPU
affinity, DuckDB settings, fastp flags, checksums, raw timing CSVs, and
correctness CSVs. The full inputs are deliberately not committed.
Recreate them from the URLs and checksums in
`benchmarks/data/fastq_exome_inputs.csv`, then run the script’s
`synthetic`, `gzip`, `fastp-exome`, `duckhts-exome`, and
`duckhts-qc-exome` modes. Keep before/after extensions as separate
artifacts; do not rebuild an old revision in place and call it the
baseline.
