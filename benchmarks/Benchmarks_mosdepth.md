DuckHTS Mosdepth Benchmark
================

<!-- Benchmarks_mosdepth.md is generated from Benchmarks_mosdepth.Rmd. -->

# Benchmark

This benchmark:

- times native `duckhts_mosdepth(...)` against upstream `mosdepth`
- uses the native benchmark harness in `scripts/mosdepth_benchmark.py`
- verifies output equality on the first run of each case
- exercises whole-genome `fast`, `default`, `fragment`, and filtered
  default mode

# Run

`make bench-mosdepth`

Override defaults with:

- `MOSDEPTH_BENCH_BAM`
- `MOSDEPTH_BENCH_RUNS`
- `MOSDEPTH_BENCH_THREADS`
- `MOSDEPTH_BENCH_WINDOW`
- `MOSDEPTH_BENCH_EXTENSION`

# Inputs

    #>                                                                       bam
    #> 1 /root/duckhts/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
    #>                                              extension runs threads window_size
    #> 1 /root/duckhts/build/release/duckhts.duckdb_extension    1       4     1000000

# Whole-Genome Cases

    #>                                    case     mode                        filter
    #> 1                     whole_genome_fast     fast               default filters
    #> 2                  whole_genome_default  default               default filters
    #> 3                 whole_genome_fragment fragment               default filters
    #> 4 whole_genome_default_mapq60_include64  default MAPQ >= 60, include_flag = 64
    #>   mosdepth_sec duckhts_sec                  speedup
    #> 1        32.94       27.66 1.19x faster vs mosdepth
    #> 2        37.25       28.76 1.30x faster vs mosdepth
    #> 3        36.15       27.95 1.29x faster vs mosdepth
    #> 4        37.59       28.71 1.31x faster vs mosdepth
    #>                                     verification
    #> 1 regions (1): PASS (230 non-zero rows, 0 diffs)
    #> 2 regions (1): PASS (230 non-zero rows, 0 diffs)
    #> 3 regions (1): PASS (230 non-zero rows, 0 diffs)
    #> 4 regions (1): PASS (230 non-zero rows, 0 diffs)

# Commands

    #>                                    case
    #> 1                     whole_genome_fast
    #> 2                  whole_genome_default
    #> 3                 whole_genome_fragment
    #> 4 whole_genome_default_mapq60_include64
    #>                                                                                                                                                                                                                                                                                                                           command
    #> 1                                        python3 'scripts/mosdepth_benchmark.py' '/root/duckhts/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam' '--mode' 'fast' '--window-size' '1000000' '--no-per-base' '--threads' '4' '--runs' '1' '--extension' '/root/duckhts/build/release/duckhts.duckdb_extension' '--verify'
    #> 2                                     python3 'scripts/mosdepth_benchmark.py' '/root/duckhts/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam' '--mode' 'default' '--window-size' '1000000' '--no-per-base' '--threads' '4' '--runs' '1' '--extension' '/root/duckhts/build/release/duckhts.duckdb_extension' '--verify'
    #> 3                                    python3 'scripts/mosdepth_benchmark.py' '/root/duckhts/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam' '--mode' 'fragment' '--window-size' '1000000' '--no-per-base' '--threads' '4' '--runs' '1' '--extension' '/root/duckhts/build/release/duckhts.duckdb_extension' '--verify'
    #> 4 python3 'scripts/mosdepth_benchmark.py' '/root/duckhts/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam' '--mode' 'default' '--window-size' '1000000' '--no-per-base' '--threads' '4' '--runs' '1' '--extension' '/root/duckhts/build/release/duckhts.duckdb_extension' '--verify' '--mapq' '60' '--include-flag' '64'

# Raw Output

## whole_genome_fast

``` text
======================================================================
Mosdepth vs native DuckHTS benchmark
======================================================================
Alignment:    NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
Mode:         fast
Reads:        172,724,240
Scope:        whole file (3,137,454,505 bp)
Threads:      4
Runs:         1
Extension:    /root/duckhts/build/release/duckhts.duckdb_extension
Window size:  1,000,000
Per-base:     disabled

mosdepth run...
  32.94s
duckhts_mosdepth run...
  27.66s

======================================================================
BENCHMARK SUMMARY
======================================================================
Tool                      Average         Best
----------------------------------------------------------------------
mosdepth                   32.94s       32.94s
duckhts_mosdepth           27.66s       27.66s
----------------------------------------------------------------------
Ratio: DuckHTS is 1.19x faster vs mosdepth

Verification:
  regions (1): PASS (230 non-zero rows, 0 diffs)

Temp dir cleaned up.
```

## whole_genome_default

``` text
======================================================================
Mosdepth vs native DuckHTS benchmark
======================================================================
Alignment:    NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
Mode:         default
Reads:        172,724,240
Scope:        whole file (3,137,454,505 bp)
Threads:      4
Runs:         1
Extension:    /root/duckhts/build/release/duckhts.duckdb_extension
Window size:  1,000,000
Per-base:     disabled

mosdepth run...
  37.25s
duckhts_mosdepth run...
  28.76s

======================================================================
BENCHMARK SUMMARY
======================================================================
Tool                      Average         Best
----------------------------------------------------------------------
mosdepth                   37.25s       37.25s
duckhts_mosdepth           28.76s       28.76s
----------------------------------------------------------------------
Ratio: DuckHTS is 1.30x faster vs mosdepth

Verification:
  regions (1): PASS (230 non-zero rows, 0 diffs)

Temp dir cleaned up.
```

## whole_genome_fragment

``` text
======================================================================
Mosdepth vs native DuckHTS benchmark
======================================================================
Alignment:    NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
Mode:         fragment
Reads:        172,724,240
Scope:        whole file (3,137,454,505 bp)
Threads:      4
Runs:         1
Extension:    /root/duckhts/build/release/duckhts.duckdb_extension
Window size:  1,000,000
Per-base:     disabled

mosdepth run...
  36.15s
duckhts_mosdepth run...
  27.95s

======================================================================
BENCHMARK SUMMARY
======================================================================
Tool                      Average         Best
----------------------------------------------------------------------
mosdepth                   36.15s       36.15s
duckhts_mosdepth           27.95s       27.95s
----------------------------------------------------------------------
Ratio: DuckHTS is 1.29x faster vs mosdepth

Verification:
  regions (1): PASS (230 non-zero rows, 0 diffs)

Temp dir cleaned up.
```

## whole_genome_default_mapq60_include64

``` text
======================================================================
Mosdepth vs native DuckHTS benchmark
======================================================================
Alignment:    NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
Mode:         default
Reads:        172,724,240
Scope:        whole file (3,137,454,505 bp)
Threads:      4
Runs:         1
Extension:    /root/duckhts/build/release/duckhts.duckdb_extension
Window size:  1,000,000
Per-base:     disabled
Include flag: 64
MAPQ:         60

mosdepth run...
  37.59s
duckhts_mosdepth run...
  28.71s

======================================================================
BENCHMARK SUMMARY
======================================================================
Tool                      Average         Best
----------------------------------------------------------------------
mosdepth                   37.59s       37.59s
duckhts_mosdepth           28.71s       28.71s
----------------------------------------------------------------------
Ratio: DuckHTS is 1.31x faster vs mosdepth

Verification:
  regions (1): PASS (230 non-zero rows, 0 diffs)

Temp dir cleaned up.
```
