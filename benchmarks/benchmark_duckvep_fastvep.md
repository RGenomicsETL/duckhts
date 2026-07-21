DuckVEP and FastVEP end-to-end comparison
================

<!-- benchmark_duckvep_fastvep.md is generated from benchmark_duckvep_fastvep.Rmd. Do not edit the .md by hand. -->

Status: current measured comparison of DuckVEP and FastVEP on the
complete GIAB HG002 GRCh38 small-variant benchmark VCF. It includes VCF
decoding and real uncompressed tabular output. DuckVEP additionally pays
for an explicit global coordinate sort. The consequence speed lane
excludes HGVS, regulatory/motif features, and supplementary annotation
for both engines.

## Result

On one pinned physical performance core, DuckVEP processed the full file
and wrote 47.63 million annotation rows in a median **64.54 seconds**
over three runs. Native FastVEP wrote 47.83 million rows in a median
**164.38 seconds** over three runs. The observed FastVEP/DuckVEP median
wall-time ratio is **2.55**.

That ratio describes these two native operational pipelines; it is not a
controlled engine-only speedup. DuckVEP enforces input order inside its
measured query and spills that sort under a 4 GB DuckDB memory limit,
while FastVEP trusts the already sorted VCF. FastVEP also computes
values for native tab fields that the current DuckVEP projection leaves
as placeholders. The exact times, row counts, and output widths are the
comparison; the report does not subtract either tool’s required work to
manufacture an artificial common kernel number.

| engine  | runs | median seconds | minimum seconds | maximum seconds | source records/s | source ALT alleles/s | output rows | output rows/s | median CPU % | maximum RSS GiB |
|:--------|-----:|---------------:|----------------:|----------------:|-----------------:|---------------------:|:------------|--------------:|-------------:|----------------:|
| DuckVEP |    3 |          64.54 |           63.75 |           66.44 |            62726 |                63466 | 47,629,345  |        737982 |           99 |            5.47 |
| FastVEP |    3 |         164.38 |          164.23 |          164.75 |            24628 |                24919 | 47,828,122  |        290961 |           99 |            3.17 |

The comparable rates use all 4,048,342 source records and all 4,096,123
source ALT alleles for both engines. DuckVEP annotates 4,095,611 literal
A/C/G/T/N alleles at 63,458 accepted alleles/s; 512 symbolic or
otherwise non-literal ALT alleles remain outside this particular
small-allele projection rather than disappearing from the shared source
denominator.

## Four-core comparison

The four-core comparison pins both processes to physical performance
cores 2, 4, 6, and 8. DuckVEP uses four DuckDB task threads while
retaining `decompression_threads := 0`; FastVEP uses four Rayon workers.
DuckVEP has three observations, while the FastVEP scaling point has one
and therefore no dispersion estimate.

| engine  | runs | median seconds | minimum seconds | maximum seconds | source ALT alleles/s | output rows/s | median CPU % | maximum RSS GiB |
|:--------|-----:|---------------:|----------------:|----------------:|---------------------:|--------------:|-------------:|----------------:|
| DuckVEP |    3 |          32.06 |           31.99 |           32.39 |               127764 |       1485631 |          251 |            5.57 |
| FastVEP |    1 |          68.09 |           68.09 |           68.09 |                60157 |        702425 |          239 |            3.17 |

DuckVEP’s four-thread median is **32.06 seconds**, versus **68.09
seconds** for FastVEP. The observed same-core-count wall-time ratio is
therefore **2.12** in DuckVEP’s favor. DuckVEP delivers a 2.01-fold
speedup over its one-thread median, or 50.3% elapsed-time parallel
efficiency, and reaches only 251% median aggregate CPU use. The measured
pipeline therefore remains limited by serial or partially parallel
reader, sort, coordination, or writer work rather than saturating four
cores with consequence computation.

The implementation contract shares the immutable resident model across
DuckDB workers while giving each worker its own mutable annotation
workspace. Measured maximum process RSS was 5.57 GiB at four threads,
versus 5.47 GiB at one. That process-level observation rules out a
fourfold increase in total RSS; it does not attribute individual
allocations to the model or worker workspaces.

### Four-core operator profile

An additional diagnostic run used the same four-core command and real
output contract with DuckDB JSON profiling enabled. It recorded 29.68
seconds of query latency and 64.27 aggregate CPU-seconds, or 2.17 cores
averaged across the query. The six largest local operator CPU totals
were:

| stage                              | aggregate CPU seconds | output rows |
|:-----------------------------------|----------------------:|:------------|
| TSV COPY                           |                 16.59 | 1           |
| transcript-label hash join         |                 14.56 | 47,629,345  |
| annotation-list expansion / UNNEST |                  8.66 | 47,629,345  |
| 17-field text projection           |                  7.78 | 47,629,345  |
| sequential read_bcf scan           |                  4.34 | 4,048,342   |
| FILTER                             |                  2.86 | 4,095,611   |

These are aggregate operator CPU times, not additive wall-clock phases,
and do not by themselves prove which operator serialized a pipeline.
They do show that this query has no chromosome-per-worker scheduler:
DuckDB schedules vector-sized work around a global sort. The largest
measured costs are the 47.63-million-row transcript-label join, TSV
writer, annotation-list expansion, and final string projection. The
one-handle sequential VCF scan is visible but smaller. The observed 251%
median process utilization is therefore not a largest-chromosome tail,
and the end-to-end scaling ceiling is not attributable to the
consequence kernel alone.

If explicit genomic parallelism is introduced later, dynamically
scheduled coordinate tiles are preferable to one static task per
chromosome. Static chromosome tasks would create exactly the
chromosome-size and transcript-density tail that this current query does
not have.

## FastVEP scaling

FastVEP was compiled with `RUSTFLAGS='-C target-cpu=native'`; its
release profile also enables optimization level 3, fat LTO, and one
code-generation unit. `RAYON_NUM_THREADS` and `taskset` independently
control its worker count and CPU affinity. Only the six physical
performance cores were used; no SMT sibling or efficiency core was mixed
into this series.

| threads | seconds | ALT alleles/s | output rows/s | CPU % |  speedup | parallel efficiency % | maximum RSS GiB |
|--------:|--------:|--------------:|--------------:|------:|---------:|----------------------:|----------------:|
|       1 |  164.75 |         24863 |        290307 |    99 | 1.000000 |                 100.0 |            3.17 |
|       2 |  100.53 |         40745 |        475760 |   161 | 1.638814 |                  81.9 |            3.17 |
|       4 |   68.09 |         60157 |        702425 |   239 | 2.419592 |                  60.5 |            3.17 |
|       6 |   58.65 |         69840 |        815484 |   289 | 2.809037 |                  46.8 |            3.17 |

Six Rayon workers reduced wall time to 58.65 seconds, but consumed only
289% aggregate CPU and delivered 2.81 times the one-thread throughput.
The measured reader, ordered result coordinator, and writer therefore
limit scaling before six cores are saturated. This report does not
extrapolate that curve beyond the six homogeneous performance cores.
Each scaling point is one observation from the same sequence of runs;
unlike the one-core comparison, the scaling table has no repeated-run
dispersion estimate.

The six-core FastVEP wall time is 5.89 seconds shorter than the one-core
DuckVEP wall time. That is useful capacity information, but it is not a
same-core-count comparison; this report has no end-to-end six-core
DuckVEP row.

DuckDB was fixed at `PRAGMA threads=1` and the complete R worker was
pinned to CPU 2 for the one-core comparison. The four-core DuckVEP row
changes only `PRAGMA threads` and the affinity set. Neither row uses
htslib decompression threads, so the report does not relabel
decompression work as annotation parallelism.

## What was timed

Both lanes begin with the same compressed, coordinate-sorted GIAB VCF on
the same local filesystem and end with a real uncompressed tab file. The
full file is intentionally read with `scan_mode := 'sequential'`;
indexed region subsetting is a separate `read_bcf(...)` transport path
and is not part of this whole-file measurement.

DuckVEP includes:

1.  R, DuckDB, extension, and resident-model initialization;
2.  `read_bcf(..., scan_mode := 'sequential', decompression_threads := 0)`
    VCF decoding with no htslib worker pool;
3.  ALT expansion and chromosome-alias mapping;
4.  an explicit global order by sequence-region ordinal, position,
    record, and ALT ordinal;
5.  transcript candidate discovery and consequence annotation at a
    5,000-base transcript distance; and
6.  projection and `COPY` of 17 tab-separated fields to a real file.

FastVEP includes its process startup, cache load, VCF reader, transcript
candidate discovery, consequence annotation, ordered result handling,
and native 17-field tab writer. It consumes the already sorted VCF and
does not perform an explicit sort.

| tool    | output_contract       | row_count  | bytes         | lines      | observed_sha256                                                  | observed_runs | byte_order_stable |
|:--------|:----------------------|:-----------|:--------------|:-----------|:-----------------------------------------------------------------|--------------:|:------------------|
| duckvep | matched_tab_17_fields | 47,629,345 | 6,174,109,710 | 47,629,346 | cab3099378b4ccb0de4c2ecabc4cf6cac58c32aeb456015c8a6657782c24e8e3 |             6 | FALSE             |
| fastvep | native_tab_17_fields  | 47,828,122 | 6,376,454,054 | 47,828,124 | aec1d995ecd519e86f449372f4c9b6975f72dc9284afc16c646778e9c14d201f |             6 | TRUE              |

FastVEP’s tab schema is its native fixed 17-field projection. DuckVEP
matches the same physical column count and writes comparable text
volume, but it does not fabricate FastVEP-only presentation values such
as codon strings, flags, or existing-variation labels. The SHA-256 value
is an observed file receipt, not an assertion that the two engines write
identical text. FastVEP produced the same bytes in all six checked
thread/repeat observations. DuckVEP produced multiple byte orders across
6 observations because the final relational join has no output
`ORDER BY`; a same-parser order-independent fingerprint over all
47,629,345 rows was identical in all 6 runs checked for multiset
stability. Its checked fingerprint is XOR 15135830789387217416,
low-32-bit sum 102267960154848972, and high-32-bit sum
102341672868090308.

The FastVEP cache contains 645,457 transcripts; the exact VEP-filtered
DuckVEP model contains 644,427. The 1,030-transcript difference remains
visible in the output cardinalities and the VEP differential below. An
intersection-only performance or concordance denominator would conceal
it.

## Concordance against Ensembl VEP 116

Speed alone is not the compatibility contract. Both engines were
compared with the same Ensembl VEP 116 executable on a deterministic
held-out ClinVar chromosome-21 sample: 1,864 variants sampled by allele
shape with seed 113. The exact 113 KB sample is retained as
`benchmarks/data/duckvep_fastvep/clinvar_chr21_seed113.vcf`;
reproduction does not depend on regenerating it from an unrecorded
mutable source. VEP ran from its indexed cache with HGVS enabled.
FastVEP ran its current native binary with `--hgvs --output-format vcf`.
DuckVEP’s oracle Parquet came from the same VEP run and the same
resident model.

| engine  | metric                              | result          | percent |
|:--------|:------------------------------------|:----------------|--------:|
| duckvep | union_keys                          | 56,998 / 56,998 | 100.000 |
| duckvep | shared_keys                         | 56,998 / 56,998 | 100.000 |
| duckvep | missing_keys                        | 0 / 56,998      |   0.000 |
| duckvep | extra_keys                          | 0 / 56,998      |   0.000 |
| duckvep | consequence_exact                   | 56,998 / 56,998 | 100.000 |
| duckvep | consequence_discordant              | 0 / 56,998      |   0.000 |
| duckvep | hgvsc_suffix_exact_including_absent | 56,998 / 56,998 | 100.000 |
| duckvep | hgvsc_suffix_discordant             | 0 / 56,998      |   0.000 |
| duckvep | hgvsp_suffix_exact_including_absent | 56,998 / 56,998 | 100.000 |
| duckvep | hgvsp_suffix_discordant             | 0 / 56,998      |   0.000 |
| fastvep | union_keys                          | 57,317 / 57,317 | 100.000 |
| fastvep | shared_keys                         | 56,998 / 57,317 |  99.443 |
| fastvep | missing_keys                        | 0 / 57,317      |   0.000 |
| fastvep | extra_keys                          | 319 / 57,317    |   0.557 |
| fastvep | consequence_exact                   | 54,045 / 57,317 |  94.291 |
| fastvep | consequence_discordant              | 2,953 / 56,998  |   5.181 |
| fastvep | hgvsc_suffix_exact_including_absent | 43,149 / 56,998 |  75.703 |
| fastvep | hgvsc_suffix_discordant             | 13,849 / 56,998 |  24.297 |
| fastvep | hgvsp_suffix_exact_including_absent | 46,847 / 56,998 |  82.191 |
| fastvep | hgvsp_suffix_discordant             | 10,151 / 56,998 |  17.809 |

DuckVEP reproduced all 56,998 VEP transcript pairs, consequence sets,
HGVSc suffixes, and HGVSp suffixes, including mutually absent HGVS
values. FastVEP emitted 319 extra transcript keys and no missing VEP
keys; 54,045 consequence sets were exact on the 57,317-key union. Its
HGVS suffix comparison is materially less exact, which is why the
FastVEP native tab speed lane must not be described as VEP HGVS parity.
Accession and version prefixes are removed before this comparison and
are not covered by these HGVS metrics.

This held-out differential is a compatibility witness for this
comparison, not a replacement for DuckVEP’s full multi-corpus
conformance ledger or statistical state exploration. The first 100
discordant FastVEP pairs are retained in
`benchmarks/data/duckvep_fastvep/conformance_examples.csv`.

Both speed lanes use the shared 5,000-base transcript-distance default.
That choice makes the tool comparison like-for-like; it does not
establish distance-invariant performance. DuckVEP’s separate
[dense-region report](duckvep_throughput.md) retains zero-, 10,000-, and
50,000-base measurements and row fingerprints so this end-to-end
workload does not become the only performance authority.

## Revisions and input receipts

| item                             | receipt                                                          |
|:---------------------------------|:-----------------------------------------------------------------|
| DuckHTS measured checkout        | 2a1c37cf2938e8226a078f19ca429ffeff84de73                         |
| DuckHTS extension binary         | 53af1444f11bd22092001fe361e36f89bd397f7fcb52af927fefb69378c5281b |
| DuckVEP benchmark worker         | 4af9f71500e4ffd081c6853ee96c7e9a298e9ed3bd52d32e3e9dfdcc9f0497e2 |
| FastVEP checkout                 | 7038e7c17708e7d2226149e78e0bb297bcc6d1d6                         |
| FastVEP native binary            | b4cb538537646a4eaa494e0ab29978e8ead73009f643e369b4f8ee447e392d5a |
| FastVEP rebuilt transcript cache | 00a3357ea30325c9d93f53ce0dabc81cb6542a0fd6d8741e895331935f89f962 |
| DuckVEP Ensembl 116 GRCh38 model | 4c2077c83958f7a300119225f91ebf288baba49497d6efe691a2a9898eb4710a |
| GIAB HG002 input VCF             | adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c |
| ClinVar 2026-07-06 source VCF    | 59a83b34d425daf58cd0dd463d6f2952f0a833ddf8fe6698fd30010642e5e1e9 |
| held-out ClinVar input VCF       | 7ecec9a7507166c2f6d4db240ba27b6f07ed9da3b0d95936882888941bcf3bf4 |
| Ensembl 116 primary FASTA        | 1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b |
| Ensembl 116 primary FASTA index  | 0998f61682f4041b11f0d156e1db6dae3e4c743e26643a3f45ea7faea70cb604 |
| Ensembl 116 GFF3                 | 08e881d96ab6385a2c31f063a018be4b2c36860b323f2724be07022deeef21ce |

The extension binary was clean-built at `c977cba`; there are no
extension, public-catalog, or build-wiring changes from that revision
through the measured `2a1c37c` checkout. The separate worker SHA-256
identifies the exact benchmark harness used for every retained DuckVEP
timing; the report aborts if its checked worker differs. FastVEP was
updated to the latest upstream checkout before the run and its
transcript cache was rebuilt from the receipt-pinned Ensembl 116 GFF3
and FASTA rather than reusing an older cache.

The host is an Intel Core i5-13500 with six physical performance cores.
CPU 2 is the single-core lane; FastVEP uses CPU sets `2`, `2,4`,
`2,4,6,8`, and `0,2,4,6,8,10`. GNU `/usr/bin/time -v` supplies wall
time, CPU use, and process-level maximum RSS. All measurements used a
warm local filesystem page cache and do not claim cold-storage or
`fsync` latency.

## Reproduction

The exact measured command lines are retained verbatim in the checked-in
GNU time files. The DuckVEP worker is
`benchmarks/benchmark_duckvep_fastvep_worker.R`; it owns the SQL and
exposes `--threads`, `--distance`, `--memory-limit`, and optional
`--profile-json`. FastVEP was built with:

``` bash
env RUSTFLAGS='-C target-cpu=native' cargo build \
  --manifest-path .sync/fastVEP/Cargo.toml --release --locked
```

FastVEP creates its binary transcript cache as a side effect of
annotation when the requested cache does not exist. The measured cache
was rebuilt without source-label or synonym overrides as follows; the
retained held-out VCF is only the trigger input for this
cache-construction pass:

``` bash
cache=/tmp/homo_sapiens_grch38_116_7038e7c.cache
test ! -e "$cache"
.sync/fastVEP/target/release/fastvep annotate \
  --input benchmarks/data/duckvep_fastvep/clinvar_chr21_seed113.vcf \
  --output /dev/null --output-format tab --no-progress \
  --gff3 /root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.116.gff3.gz \
  --fasta /root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --transcript-cache "$cache"
test "$(sha256sum "$cache" | cut -d' ' -f1)" = \
  00a3357ea30325c9d93f53ce0dabc81cb6542a0fd6d8741e895331935f89f962
```

The single-core invocations were:

``` bash
/usr/bin/time -v \
  -o /tmp/duckvep-fastvep-benchmark/duckvep_giab_1_final.time \
  taskset -c 2 Rscript benchmarks/benchmark_duckvep_fastvep_worker.R \
  --extension build/release/duckhts.duckdb_extension \
  --model /root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb \
  --input /root/duckvep/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --output /tmp/duckvep-fastvep-benchmark/duckvep_giab_1_final.tab \
  --threads 1 --distance 5000 --memory-limit 4GB

/usr/bin/time -v \
  -o /tmp/duckvep-fastvep-benchmark/fastvep_giab_1.time \
  taskset -c 2 env RAYON_NUM_THREADS=1 \
  .sync/fastVEP/target/release/fastvep annotate \
  --input /root/duckvep/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --output /tmp/duckvep-fastvep-benchmark/fastvep_giab_1.tab \
  --output-format tab --distance 5000 --no-progress \
  --transcript-cache /root/duckvep/data/cache/fastvep/homo_sapiens_grch38_116_7038e7c.cache

Rscript benchmarks/benchmark_duckvep_fastvep_receipt.R \
  --input /tmp/duckvep-fastvep-benchmark/duckvep_giab_1_final.tab \
  --tool duckvep --output-contract matched_tab_17_fields \
  --threads 1 --run 1 \
  --timing-file benchmarks/data/duckvep_fastvep/duckvep_giab_threads1.time \
  --output /tmp/duckvep-fastvep-benchmark/duckvep_receipt_run1.csv

Rscript benchmarks/benchmark_duckvep_fastvep_receipt.R \
  --input /tmp/duckvep-fastvep-benchmark/fastvep_giab_1.tab \
  --tool fastvep --output-contract native_tab_17_fields \
  --threads 1 --run 1 --skip-lines 1 \
  --timing-file benchmarks/data/duckvep_fastvep/fastvep_giab_threads1.time \
  --output /tmp/duckvep-fastvep-benchmark/fastvep_receipt_run1.csv
```

The four-core DuckVEP measurements change the affinity set and DuckDB
task count without enabling htslib decompression workers. Run suffixes
and timing destinations distinguish the three repetitions:

``` bash
/usr/bin/time -v \
  -o /tmp/duckvep-fastvep-benchmark/duckvep_giab_4_run1.time \
  taskset -c 2,4,6,8 \
  Rscript benchmarks/benchmark_duckvep_fastvep_worker.R \
  --extension build/release/duckhts.duckdb_extension \
  --model /root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb \
  --input /root/duckvep/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --output /tmp/duckvep-fastvep-benchmark/duckvep_giab_4.tab \
  --threads 4 --distance 5000 --memory-limit 4GB

taskset -c 2,4,6,8 \
  Rscript benchmarks/benchmark_duckvep_fastvep_worker.R \
  --extension build/release/duckhts.duckdb_extension \
  --model /root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb \
  --input /root/duckvep/data/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --output /tmp/duckvep-fastvep-benchmark/duckvep_giab_4_profile.tab \
  --threads 4 --distance 5000 --memory-limit 4GB \
  --profile-json /tmp/duckvep-fastvep-benchmark/duckvep_giab_4_profile.json

cp /tmp/duckvep-fastvep-benchmark/duckvep_giab_4_profile.json \
  benchmarks/data/duckvep_fastvep/duckvep_giab_threads4_profile.json
```

The same receipt command is run before an output is replaced. Its
per-run rows are retained in
`benchmarks/data/duckvep_fastvep/output_observations.csv`; the render
derives observed-run counts, byte-order stability, and multiset
stability from those rows rather than from an aggregate assertion.

`benchmarks/benchmark_duckvep_fastvep_conformance.R` collects the
VEP-oracle comparison without dropping missing or extra transcript keys.
The held-out oracle, FastVEP HGVS output, and collector were produced
with:

``` bash
out=/tmp/duckvep-fastvep-benchmark/conformance200.CCg7sP
Rscript test/duckvep/conformance/corpus_differential.R \
  --corpus clinvar_chr21_seed113 \
  --vcf /root/duckvep/data/corpora/clinvar/20260706/clinvar_20260706.vcf.gz \
  --source-version 20260706 \
  --source-checksum sha256:59a83b34d425daf58cd0dd463d6f2952f0a833ddf8fe6698fd30010642e5e1e9 \
  --cache-dir /root/duckvep/data/vep_cache \
  --assembly GRCh38 --species homo_sapiens \
  --fasta /root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --database /root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb \
  --model-sql= --model-name fastvep_conformance \
  --sample-per-shape 200 --max-allele-length 50 --seed 113 --chrom 21 \
  --distance 5000 --hgvs --fork 8 --vep-buffer-size 5000 \
  --sample-vcf "$out/clinvar_chr21_seed113.vcf" --keep-sample-vcf \
  --annotations-out "$out/duckvep_annotations.parquet" \
  --hgvs-out "$out/duckvep_hgvs_summary.csv" \
  --hgvs-pairs-out "$out/duckvep_hgvs_pairs.parquet" \
  --extension build/release/duckhts.duckdb_extension

test "$(sha256sum "$out/clinvar_chr21_seed113.vcf" | cut -d' ' -f1)" = \
  7ecec9a7507166c2f6d4db240ba27b6f07ed9da3b0d95936882888941bcf3bf4

taskset -c 2 env RAYON_NUM_THREADS=1 \
  .sync/fastVEP/target/release/fastvep annotate \
  --input "$out/clinvar_chr21_seed113.vcf" --output "$out/fastvep.vcf" \
  --output-format vcf --hgvs --distance 5000 --no-progress \
  --fasta /root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --transcript-cache /root/duckvep/data/cache/fastvep/homo_sapiens_grch38_116_7038e7c.cache

Rscript benchmarks/benchmark_duckvep_fastvep_conformance.R \
  --extension build/release/duckhts.duckdb_extension \
  --fastvep-vcf "$out/fastvep.vcf" \
  --input-vcf "$out/clinvar_chr21_seed113.vcf" \
  --oracle-parquet "$out/duckvep_annotations.parquet" \
  --output "$out/conformance.csv" \
  --examples "$out/conformance_examples.csv"
```

Render this report only from the R Markdown source:

``` bash
Rscript -e "rmarkdown::render('benchmarks/benchmark_duckvep_fastvep.Rmd', output_format='github_document')"
```
