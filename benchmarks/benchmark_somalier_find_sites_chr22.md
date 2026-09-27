Somalier find-sites: 1000 Genomes 30x phased chr22 and chrX
================

This is a whole-chromosome differential of DuckHTS and pinned Somalier
[v0.3.4](https://github.com/brentp/somalier/tree/ff58fdade8f4f8293d904f10e0a4a13f1fac808d)
(binary SHA-256
`18717c205a9c4b65d479f1d2cf069a30b047a5378edf544d403d4081f06a3a78`). The
unmodified [1000 Genomes high-coverage 30x phased chr22
VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz)
is staged through `duckhtsbench` (SHA-256
`70ff24cfd01eab24c6f4c6e825f55a2b3475ff595f6e12e179c4ef06b9cc4df6`,
519,930,289 bytes, GRCh38). Its declared INFO/AF and INFO/AN are used
without recomputing from genotypes. Somalier skips sample genotypes;
DuckHTS reads the physical VCF in sequential mode with `samples := ''`.
The cache VCF is not committed.

## Selection rules

| Rule          | Somalier v0.3.4 and DuckHTS default                                                                                           | DuckHTS named alternative                                              |
|:--------------|:------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------|
| Human contigs | Exclude autosomal REF=C; X zero-based positions 2,781,479..154,931,044                                                        | `compatibility_mode = "generic"` does not classify human sex contigs   |
| AF rank ties  | Retain input order on equal AF scores; a reversed-position two-site probe selects the first input record with pinned Somalier | `tie_order = "lexical"` orders by position and allele                  |
| Spacing       | Greedy autosomal distance 10,000 by default; X/Y selections do not update spacing state                                       | `sex_spacing = "enforced"` applies X distance 1,000 and Y distance 200 |
| Caps          | Up to 65,535 autosomal, 10,001 X, and 5,001 Y selections after ranking                                                        | Named per-sex maximum arguments                                        |

- Source revision: 99d2f4438e52fae5d47c9e3e0b523408a9f67a5b. Workload:
  `--min-AN 6000`, `--snp-dist 10000`, minimum AF 0.15, AF target 0.48,
  no interval or gnotate exclusion. One run per tool, warm filesystem
  cache, GNU `/usr/bin/time -v`.
- Denominators: 1,070,401 input records; 40,521 Somalier candidates and
  40,521 DuckHTS candidates before nearby-variant exclusion; 38,058
  chr22 spacing inputs. The native spacing kernel’s per-list limit is
  1,000,000 candidates; this chromosome’s largest list is below it.
- Selected: 2402 Somalier sites and 2402 DuckHTS sites; 2402 shared
  four-field keys `(region, position, source_ref, source_alt)`.
  Input-order and lexical tie policies differed on 0 selected sites on
  this chromosome.
- Threads: Somalier uses two reader/decompression threads and one writer
  thread (pinned source); DuckHTS uses one DuckDB thread and a
  sequential reader. The timed DuckHTS process also verifies its cached
  input and starts R/DuckDB; Somalier timing covers the CLI only. These
  single-process measurements are not matched-concurrency speedups.

| gate            | records | somalier_records |
|:----------------|--------:|-----------------:|
| src             | 1070401 |          1070401 |
| eligible        |  742856 |               NA |
| pass_snv        |  625286 |               NA |
| af_an           |   44354 |               NA |
| annotation_gate |   44354 |               NA |
| qc_gate         |   40521 |               NA |
| interval_gate   |   40521 |               NA |
| gated           |   40521 |            40521 |
| indels          |   31422 |               NA |
| snps            |  133277 |               NA |
| indel_clear     |   39521 |               NA |
| neighbors       |   38058 |               NA |
| kept            |    2402 |               NA |
| selected        |    2402 |             2402 |

Whole-chromosome stage counts. Indels and snps are exclusion-neighbor
inventories, not sequential pass counts. NA denotes stages the pinned
Somalier executable does not expose; its own log reports the gated
candidate count.

| tool            | wall_seconds | peak_rss_kib |
|:----------------|-------------:|-------------:|
| Somalier v0.3.4 |         6.07 |       108912 |
| DuckHTS         |        11.08 |       452576 |

Wall seconds and peak RSS KiB for each timed process, including its own
initialization and I/O. Single runs are not a variance estimate.

| class         | records |
|:--------------|--------:|
| DuckHTS only  |       0 |
| Somalier only |       0 |

Keyed disagreement classes (all records counted).

The retained `somalier_find_sites_chr22_disagreements.tsv` contains the
first 100 disagreement keys and origins; 0 rows are present. An empty
table means no disagreements, not an absent comparison. The
`somalier_find_sites_chr22_observations.tsv` and gate table retain the
full comparison denominators. The previous
[`benchmark_somalier_find_sites.md`](benchmark_somalier_find_sites.md)
uses two small fixture slices and is not a comparable whole-chromosome
baseline.

## Chromosome X, including physical PAR records

The unmodified [1000 Genomes 30x phased chrX
VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz)
(SHA-256
`f42c786cba67ff74c30e9bd37254b5a605a5957cfe3a5dc54e7dd41d9c9d8edc`,
2,800,398,832 bytes) contains 118,327 records in PAR1 (before POS
2,781,480) and 19,227 records in PAR2 (after POS 154,931,045). The
physical first and last positions are 13,189 and 156,030,223. Somalier’s
pinned X-range rule admits only zero-based positions
2,781,479..154,931,044. The 137,554 PAR records account for the drop
from source to eligible records. Somalier’s log says “0 autosomal
variants” for this X-only input; the actual compressed output VCF, not
that autosomal-only log counter, is used for the selected-site
denominator and keyed comparison.

- Source revision: 99d2f4438e52fae5d47c9e3e0b523408a9f67a5b.
- Denominators: 2,858,184 input records; 250,641 Somalier candidates and
  250,641 DuckHTS gated candidates; 234,367 spacing inputs (below the
  1,000,000 limit).
- Selected: 10,001 Somalier and 10,001 DuckHTS sites; 10,001 matching
  four-field keys. Input and lexical tie policies differed on 0 selected
  sites. 0 DuckHTS-only and 0 Somalier-only keys.

| gate            | records | somalier_records |
|:----------------|--------:|-----------------:|
| src             | 2858184 |          2858184 |
| eligible        | 2720630 |               NA |
| pass_snv        | 2286442 |               NA |
| af_an           |  250641 |               NA |
| annotation_gate |  250641 |               NA |
| qc_gate         |  250641 |               NA |
| interval_gate   |  250641 |               NA |
| gated           |  250641 |           250641 |
| indels          |  110112 |               NA |
| snps            |  441703 |               NA |
| indel_clear     |  243716 |               NA |
| neighbors       |  234367 |               NA |
| kept            |  234367 |               NA |
| selected        |   10001 |            10001 |

chrX stage counts; NA denotes stages not reported by the pinned
executable.

| tool            | wall_seconds | peak_rss_kib |
|:----------------|-------------:|-------------:|
| Somalier v0.3.4 |        15.51 |       436696 |
| DuckHTS         |        23.25 |      1226680 |

chrX whole-process wall seconds and peak RSS KiB; thread and measurement
caveats above also apply.

| class         | records |
|:--------------|--------:|
| DuckHTS only  |       0 |
| Somalier only |       0 |

chrX keyed disagreement classes. First 100 rows in
somalier_find_sites_chrx_disagreements.tsv; retained: 0 rows.

## Operator profile and matched-thread selection

The before/after chrX profile uses source revisions `396f9599` and
`99d2f443` on the same 2,858,184 physical VCF records, with one DuckDB
thread, warm filesystem cache, no sample genotype decode, and the same
250,641 gated candidates and 234,367 spacing inputs. Each complete query
was profiled with DuckDB JSON profiling (query latency and operator
timing); `EXPLAIN ANALYZE` also checked the ordered-stream plan. The 14
named gate counts were collected by one query with shared materialized
CTEs. Operator times are pipelined and **not additive**. No gate count
or keyed denominator was dropped to obtain these timings.

| component                               | baseline_seconds | ordered_seconds | measurement                                              |
|:----------------------------------------|-----------------:|----------------:|:---------------------------------------------------------|
| Full selection query                    |         219.6400 |         21.2760 | DuckDB JSON profile latency                              |
| All named gate counts                   |         220.3170 |         21.3400 | DuckDB JSON profile latency                              |
| READ_BCF scan and projected INFO decode |          19.2370 |         19.3300 | READ_BCF operator in selection                           |
| Eligible position and allele filters    |           0.0970 |          0.0930 | Two FILTER operators in gate-count query                 |
| PASS SNV gate                           |           0.2080 |          0.2060 | FILTER operator in gate-count query                      |
| AF and AN gate                          |           0.0170 |          0.0170 | FILTER operator in gate-count query                      |
| Annotation gate                         |           0.0050 |          0.0050 | FILTER operator in gate-count query                      |
| QC gate                                 |           0.0010 |          0.0010 | FILTER operator in gate-count query                      |
| Optional interval and gnotate gates     |           0.0012 |          0.0012 | No intervals or gnotate input; remaining FILTER operator |
| Indel inventory filter                  |           0.0236 |          0.0236 | FILTER operator in gate-count query                      |
| SNP inventory filter                    |           0.0340 |          0.0365 | FILTER operator in gate-count query                      |
| Nearby-indel exclusion                  |          40.9320 |          0.0740 | Range HASH_JOIN versus ordered WINDOW in selection       |
| Nearby-SNP counting                     |         158.5790 |          0.7920 | Range HASH_JOIN versus grouped RANGE WINDOW in selection |
| AF ranking window                       |           0.1240 |          0.1250 | WINDOW operator in selection                             |
| List aggregation                        |           0.0110 |          0.0130 | HASH_GROUP_BY with list in selection                     |
| Spacing kernel projection               |           0.0064 |          0.0003 | Single-row mask PROJECTION in selection                  |
| Selection cap filter                    |           0.0001 |          0.0004 | FILTER operator in gate-count query                      |

chrX seconds. Named gate rows are FILTER operators; join/window rows
isolate nearby-variant work; the mask projection bounds the spacing
kernel and its projection rather than timing the kernel alone.

A separate `count(*)` over `read_bcf` took 16.700 s; decoding and
summing only `INFO_AF[1]` and `INFO_AN` took 18.265 s. The selector’s
`READ_BCF` operator took 19.237 s before and 19.330 s after. Its
projected columns are `CHROM`, `POS`, `REF`, `ALT`, `FILTER`, `INFO_AF`,
`INFO_AN`, `INFO_BaseQRankSum`, `INFO_ClippingRankSum`, `INFO_FS`,
`INFO_MQ`, `INFO_MQRankSum`, `INFO_QD`, and `INFO_ReadPosRankSum`; no
genotype columns or other INFO fields are materialized. These are
separate full scans, not additive decode components. The positional and
QC predicates run after the reader; this table function does not push
them into VCF decoding. Sampling with `perf` was unavailable
(`perf_event_paranoid=4`).

For the matched-thread runs below, DuckHTS uses one SQL worker and
either zero or two htslib BGZF decompression workers; Somalier uses its
pinned two readers and one writer. DuckDB is set to one thread in both
DuckHTS runs. A separate warm-cache process profiles one selection per
row; `process_seconds` adds R startup, the connection, and result
transfer to `selection_seconds`, but excludes the input cache checksum.
Somalier’s CLI process includes its own startup and VCF writer; these
distinct scopes are reported, not treated as a controlled throughput
ratio. Values are single runs, not variance estimates. The official
one-thread `select` workflow *includes* cache verification and took
11.08 s on chr22 and 23.25 s on chrX. Isolated warm-cache
`duckhts_bench_fetch()` calls took 0.285 s and 1.440 s, respectively; R
and library startup added about 0.3 s and opening a DuckDB connection
about 0.04 s. All 2,402 chr22 and 10,001 chrX keys also agree between
one-worker and three-worker DuckHTS selections.

| chromosome | tool            | active_threads          | selection_seconds | process_seconds | peak_rss_kib | input_records | selected_sites |
|:-----------|:----------------|:------------------------|------------------:|----------------:|-------------:|--------------:|---------------:|
| chr22      | Somalier v0.3.4 | 2 readers + 1 writer    |                NA |            6.07 |       108912 |       1070401 |           2402 |
| chr22      | DuckHTS         | 1 SQL + 0 decompression |            10.404 |           10.81 |       480948 |       1070401 |           2402 |
| chr22      | DuckHTS         | 1 SQL + 2 decompression |             5.001 |            5.43 |       481220 |       1070401 |           2402 |
| chrX       | Somalier v0.3.4 | 2 readers + 1 writer    |                NA |           15.51 |       436696 |       2858184 |          10001 |
| chrX       | DuckHTS         | 1 SQL + 0 decompression |            21.183 |           21.67 |      1254476 |       2858184 |          10001 |
| chrX       | DuckHTS         | 1 SQL + 2 decompression |            13.379 |           13.81 |      1256376 |       2858184 |          10001 |

Warm-cache whole-chromosome selection: input/output denominators, wall
seconds and peak process RSS KiB. NA: Somalier has no separate
in-process SQL selection timing.

The ordered position passes reduce the chrX selection query from 219.640
s to 21.276 s at one thread. At three active workers, DuckHTS takes
13.81 s versus Somalier’s 15.51 s; peak RSS remains about 1.25 GB rather
than 0.44 GB. `READ_BCF` and its necessary projected INFO fields now
dominate the single-thread wall time. Reducing this memory or decode
cost would require reader-level filtering or a different candidate
representation, not another change to the nearby-variant join.

Reproduce from the repository root with installed `Rduckhts` and
`duckhtsbench`, the pinned binary staged in `duckhtsbench`, and GNU
time:

``` bash
export DUCKHTSBENCH_REGISTRY="$PWD/r/duckhtsbench/inst/benchmark_registry.tsv"
Rscript test/scripts/somalier_find_sites_chromosome.R upstream "$PWD" chr22
time_file="$(Rscript -e 'cat(file.path(dirname(duckhtsbench::duckhts_bench_artifact_path("somalier_find_sites_chr22_source")), "chr22.duckhts.time"))')"
/usr/bin/time -v -o "$time_file" \
  Rscript test/scripts/somalier_find_sites_chromosome.R select "$PWD" chr22
Rscript test/scripts/somalier_find_sites_chromosome.R compare "$PWD" chr22
Rscript test/scripts/somalier_find_sites_chromosome.R upstream "$PWD" chrX
x_time_file="$(Rscript -e 'cat(file.path(dirname(duckhtsbench::duckhts_bench_artifact_path("somalier_find_sites_chrx_source")), "chrX.duckhts.time"))')"
/usr/bin/time -v -o "$x_time_file" \
  Rscript test/scripts/somalier_find_sites_chromosome.R select "$PWD" chrX
Rscript test/scripts/somalier_find_sites_chromosome.R compare "$PWD" chrX
# Selection-only process timings exclude the cache checksum:
for chrom in chr22 chrX; do
  for decompression_workers in 0 2; do
    /usr/bin/time -v Rscript scripts/benchmark_somalier_find_sites_profile.R \
      "$chrom" "$decompression_workers"
  done
done
cd benchmarks && Rscript -e "rmarkdown::render('benchmark_somalier_find_sites_chr22.Rmd')"
```
