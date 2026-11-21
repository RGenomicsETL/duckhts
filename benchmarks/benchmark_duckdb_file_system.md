read_bam: htslib file access and the DuckDB file system over HTTP
================

`read_bam` reads a file in one of two ways. By default htslib opens it
with its own file access, which for an `http://` URL is htslib’s libcurl
backend. After `SET duckhts_file_system = 'duckdb'` the same htslib
decoder reads through the DuckDB file system of the connection, which
for an `http://` URL is the httpfs extension (`src/duckdb_fs_hfile.c`).
This report compares the two routes on the same file, the same server
and the same queries: record counts, wall time, peak memory, and the
HTTP requests each route sends.

## What is held equal

Both routes run the same BAM decoder with the same htslib tuning.
DuckHTS chooses that tuning from the path the user gave, not from the
route (`src/include/hts_io_tuning.h`):

| Scan                            | hFILE block (MiB) | BGZF block cache (MiB) |
|:--------------------------------|------------------:|-----------------------:|
| indexed region or parallel scan |                 1 |                     64 |
| single-handle streaming scan    |                 4 |                    256 |

The hFILE block is the size of each read htslib asks its file backend
for. `SET duckhts_remote_block_size` replaces it; every section but the
one that varies it leaves that setting at 0, the automatic size in the
table. The BGZF cache size has no setting. htslib’s libcurl backend has
further controls that this benchmark leaves at their defaults: the retry
and low-speed environment variables (`HTS_RETRY_MAX`, `HTS_RETRY_DELAY`,
`HTS_RETRY_MAX_DELAY`, `HTS_LOW_SPEED_LIMIT`, `HTS_LOW_SPEED_TIME`) and
its rule to read on, without a new request, across a forward seek
shorter than 1 MB (`MIN_SEEK_FORWARD` in
`third_party/htslib/hfile_libcurl.c`). httpfs runs with its default
settings, except in the section that varies them.

Both routes load httpfs before the timed query, so neither pays for
loading it.

## Inputs and method

| item                        | value                                                                                   |
|:----------------------------|:----------------------------------------------------------------------------------------|
| source revision             | 35dface142aac9fc18c0f93ddfeee767ec5a0563                                                |
| source tree                 | clean                                                                                   |
| extension md5               | 0ed3a10111594faee8ef4e7603bb29bf                                                        |
| htslib                      | 1.24                                                                                    |
| R                           | R version 4.6.2 (2026-10-29)                                                            |
| DuckDB R package            | 2.0.0 (DuckDB v2.0.0)                                                                   |
| httpfs                      | 77392a7207, default settings except where a section varies them                         |
| CPU                         | 13th Gen Intel(R) Core(TM) i5-13500                                                     |
| logical CPUs                | 20                                                                                      |
| DuckDB memory limit         | 2 GiB                                                                                   |
| DuckDB temp directory limit | 0 B (no spill allowed)                                                                  |
| server                      | scripts/benchmark_range_server.py on 127.0.0.1, one thread per connection               |
| repetitions                 | 5 per region workload, 3 per scan workload; one fresh R process each, routes alternated |

The inputs are the chromosome 20 alignments of 1000 Genomes sample
HG00403 and their doubled and quadrupled copies, staged by
`r/duckhtsbench`. `samtools view -c` gives the record count that every
measurement must return:

| scale | artifact                         |         bytes | index_bytes |
|:------|:---------------------------------|--------------:|------------:|
| 1x    | insert_size_hg00403_chr20_bam    |   907,071,774 |     220,800 |
| 2x    | insert_size_hg00403_chr20_bam_2x | 1,151,649,969 |     222,736 |
| 4x    | insert_size_hg00403_chr20_bam_4x | 1,494,806,083 |     226,328 |

| Scale | Region                  | Records (samtools) |
|:------|:------------------------|-------------------:|
| 1x    | chr20:10000000-10100000 |             23,544 |
| 1x    | chr20:10000000-12500000 |            582,233 |
| 1x    | chr20:10000000-15000000 |          1,169,601 |
| 1x    | chr20:10000000-20000000 |          2,367,216 |
| 1x    | whole file              |         16,816,430 |
| 2x    | whole file              |         33,632,860 |
| 4x    | whole file              |         67,265,720 |

Each query is
`SELECT count(*), sum(hash(QNAME, FLAG, POS, CIGAR)) FROM read_bam(url [, region := ...])`.
It decodes four fields of every record and returns one row. A run passes
only when its count equals the samtools count, its hash equals the hash
of every other run on that input, DuckDB spilled nothing, and its peak
resident memory stayed under the ceiling set before measuring.

The server adds a fixed delay before it answers each request: 0 ms, and
20 ms to stand for a network round trip. It does not limit bandwidth.
Client and server share one host, so these timings show what each route
costs per request and per byte on a loopback connection; they are not
timings for an object store or a wide-area network.

## Memory budget

The idle process (extension and httpfs loaded, nothing read) peaks at
153 MiB. The ceiling for a run is that overhead plus three times the I/O
state of its open handles: one hFILE block, one BGZF block cache and one
32 MiB read-ahead buffer per handle. A region query and a single-thread
scan hold one handle; a scan on four threads holds four. A run with a
block size set budgets that size for the hFILE block. The query itself
keeps no rows: it aggregates as it reads.

## Region queries

One thread, the 1× file, regions that start at chr20:10,000,000.

| region                  | scale | delay ms | threads | route  |   records | median s | range s     | peak RSS MiB | ceiling MiB | GETs | open-ended | closed ranges | largest range MiB | body MiB written |
|:------------------------|:------|---------:|--------:|:-------|----------:|---------:|:------------|-------------:|------------:|-----:|-----------:|--------------:|------------------:|-----------------:|
| chr20:10000000-10100000 | 1x    |        0 |       1 | htslib |    23,544 |    0.049 | 0.048-0.071 |      154.430 |      444.48 |   22 |          5 |             0 |                NA |           27.389 |
| chr20:10000000-10100000 | 1x    |        0 |       1 | duckdb |    23,544 |    0.332 | 0.321-0.335 |      157.059 |      444.48 |   20 |          0 |            20 |                 1 |            4.257 |
| chr20:10000000-12500000 | 1x    |        0 |       1 | htslib |   582,233 |    0.108 | 0.097-0.116 |      154.336 |      444.48 |   22 |          5 |             0 |                NA |           54.889 |
| chr20:10000000-12500000 | 1x    |        0 |       1 | duckdb |   582,233 |    0.412 | 0.408-0.421 |      157.180 |      444.48 |   45 |          0 |            45 |                 1 |           29.257 |
| chr20:10000000-15000000 | 1x    |        0 |       1 | htslib | 1,169,601 |    0.179 | 0.167-0.190 |      154.320 |      444.48 |   22 |          5 |             0 |                NA |           80.389 |
| chr20:10000000-15000000 | 1x    |        0 |       1 | duckdb | 1,169,601 |    0.486 | 0.483-0.498 |      157.359 |      444.48 |   72 |          0 |            72 |                 1 |           56.257 |
| chr20:10000000-20000000 | 1x    |        0 |       1 | htslib | 2,367,216 |    0.316 | 0.307-0.329 |      155.785 |      444.48 |   22 |          5 |             0 |                NA |          137.451 |
| chr20:10000000-20000000 | 1x    |        0 |       1 | duckdb | 2,367,216 |    0.646 | 0.639-0.662 |      157.145 |      444.48 |  128 |          0 |           128 |                 1 |          112.257 |
| chr20:10000000-10100000 | 1x    |       20 |       1 | htslib |    23,544 |    0.529 | 0.522-0.534 |      155.406 |      444.48 |   22 |          5 |             0 |                NA |           40.826 |
| chr20:10000000-10100000 | 1x    |       20 |       1 | duckdb |    23,544 |    0.898 | 0.889-0.901 |      157.203 |      444.48 |   20 |          0 |            20 |                 1 |            4.257 |
| chr20:10000000-12500000 | 1x    |       20 |       1 | htslib |   582,233 |    0.593 | 0.588-0.603 |      156.965 |      444.48 |   22 |          5 |             0 |                NA |           69.139 |
| chr20:10000000-12500000 | 1x    |       20 |       1 | duckdb |   582,233 |    1.775 | 1.444-1.781 |      156.828 |      444.48 |   45 |          0 |            45 |                 1 |           29.257 |
| chr20:10000000-15000000 | 1x    |       20 |       1 | htslib | 1,169,601 |    0.661 | 0.629-0.664 |      155.965 |      444.48 |   22 |          5 |             0 |                NA |           95.451 |
| chr20:10000000-15000000 | 1x    |       20 |       1 | duckdb | 1,169,601 |    2.687 | 2.269-2.741 |      159.141 |      444.48 |   72 |          0 |            72 |                 1 |           56.257 |
| chr20:10000000-20000000 | 1x    |       20 |       1 | htslib | 2,367,216 |    0.787 | 0.782-0.816 |      155.602 |      444.48 |   22 |          5 |             0 |                NA |          151.326 |
| chr20:10000000-20000000 | 1x    |       20 |       1 | duckdb | 2,367,216 |    4.668 | 4.637-4.697 |      156.949 |      444.48 |  128 |          0 |           128 |                 1 |          112.257 |

“Open-ended” counts `Range: bytes=N-` requests, which ask for everything
from an offset to the end of the file. “Body MiB written” is what the
server handed to the socket. For an open-ended range that the client
abandons, it includes bytes the client never read, so on the htslib
route it is an upper bound on the bytes received.

## Whole-file scans

Served with no added delay at 1×, 2× and 4×, on one and four DuckDB
threads. The last two rows read the 1× file from the local disk through
each route.

| region     | scale | delay ms | threads | route  |    records | median s | range s     | peak RSS MiB | ceiling MiB |  GETs | open-ended | closed ranges | largest range MiB | body MiB written |
|:-----------|:------|---------:|--------:|:-------|-----------:|---------:|:------------|-------------:|------------:|------:|-----------:|--------------:|------------------:|-----------------:|
| whole file | 1x    |        0 |       1 | htslib | 16,816,430 |    2.053 | 2.049-2.066 |      154.609 |    1,029.48 |    33 |         16 |             0 |                NA |          994.859 |
| whole file | 1x    |        0 |       1 | duckdb | 16,816,430 |    2.711 | 2.675-2.727 |      155.508 |    1,029.48 |   899 |          0 |           899 |                 1 |          882.228 |
| whole file | 2x    |        0 |       1 | htslib | 33,632,860 |    3.858 | 3.838-3.907 |      154.195 |    1,029.48 |    42 |         25 |             0 |                NA |        1,321.769 |
| whole file | 2x    |        0 |       1 | duckdb | 33,632,860 |    4.671 | 4.658-4.713 |      155.551 |    1,029.48 | 1,152 |          0 |         1,152 |                 1 |        1,134.382 |
| whole file | 4x    |        0 |       1 | htslib | 67,265,720 |    7.473 | 7.419-7.510 |      154.520 |    1,029.48 |    44 |         27 |             0 |                NA |        1,672.708 |
| whole file | 4x    |        0 |       1 | duckdb | 67,265,720 |    8.492 | 8.450-8.506 |      155.945 |    1,029.48 | 1,481 |          0 |         1,481 |                 1 |        1,463.623 |
| whole file | 1x    |        0 |       4 | htslib | 16,816,430 |    0.851 | 0.841-0.859 |      169.078 |    1,317.48 |    96 |         64 |             0 |                NA |        1,405.910 |
| whole file | 1x    |        0 |       4 | duckdb | 16,816,430 |    1.306 | 1.297-1.327 |      168.070 |    1,317.48 |   930 |          0 |           930 |                 1 |          898.801 |
| whole file | 2x    |        0 |       4 | htslib | 33,632,860 |    1.502 | 1.421-1.623 |      172.676 |    1,317.48 |   112 |         80 |             0 |                NA |        1,816.820 |
| whole file | 2x    |        0 |       4 | duckdb | 33,632,860 |    1.936 | 1.898-1.942 |      168.559 |    1,317.48 | 1,203 |          0 |         1,203 |                 1 |        1,172.183 |
| whole file | 4x    |        0 |       4 | htslib | 67,265,720 |    2.949 | 2.939-3.028 |      169.590 |    1,317.48 |   129 |         97 |             0 |                NA |        2,313.083 |
| whole file | 4x    |        0 |       4 | duckdb | 67,265,720 |    3.219 | 3.215-3.256 |      171.293 |    1,317.48 | 1,536 |          0 |         1,536 |                 1 |        1,505.026 |
| whole file | 1x    |       NA |       1 | htslib | 16,816,430 |    1.942 | 1.918-2.005 |      146.785 |    1,029.48 |     0 |          0 |             0 |                NA |            0.000 |
| whole file | 1x    |       NA |       1 | duckdb | 16,816,430 |    1.923 | 1.921-1.948 |      147.180 |    1,029.48 |     0 |          0 |             0 |                NA |            0.000 |

## httpfs settings

The 100 kb and 10 Mb regions of the 1× file on one thread, on the DuckDB
route under httpfs settings that each differ from the defaults in one
respect, 3 runs each. The `htslib` rows are the htslib route at its
defaults, for reference.

| region                  | delay ms | setting                    | SET                                                                                               | median s | range s     | peak RSS MiB | GETs | largest range MiB | body MiB written |
|:------------------------|---------:|:---------------------------|:--------------------------------------------------------------------------------------------------|---------:|:------------|-------------:|-----:|------------------:|-----------------:|
| chr20:10000000-10100000 |        0 | htslib                     |                                                                                                   |    0.046 | 0.035-0.051 |      154.352 |   22 |                NA |           26.451 |
| chr20:10000000-10100000 |        0 | defaults                   |                                                                                                   |    0.328 | 0.324-0.334 |      157.359 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | force_download             | force_download = true                                                                             |    0.739 | 0.738-0.751 |    1,023.211 |    8 |                NA |        1,730.523 |
| chr20:10000000-10100000 |        0 | no_external_file_cache     | enable_external_file_cache = false                                                                |    0.335 | 0.334-0.336 |      157.258 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    0.329 | 0.322-0.337 |      157.090 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    0.330 | 0.323-0.333 |      156.957 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    0.336 | 0.328-0.338 |      157.586 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | no_keep_alive              | http_keep_alive = false                                                                           |    0.072 | 0.046-0.081 |      156.320 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    0.069 | 0.066-0.102 |      157.906 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    0.334 | 0.331-0.339 |      157.543 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | htslib                     |                                                                                                   |    0.528 | 0.524-0.531 |      155.504 |   22 |                NA |           41.576 |
| chr20:10000000-10100000 |       20 | defaults                   |                                                                                                   |    0.902 | 0.901-0.903 |      156.891 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | force_download             | force_download = true                                                                             |    0.938 | 0.926-0.943 |    1,023.270 |    8 |                NA |        1,730.523 |
| chr20:10000000-10100000 |       20 | no_external_file_cache     | enable_external_file_cache = false                                                                |    0.898 | 0.891-0.904 |      157.836 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    0.901 | 0.896-0.903 |      157.125 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    0.897 | 0.893-0.900 |      159.211 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    0.902 | 0.893-0.902 |      157.352 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | no_keep_alive              | http_keep_alive = false                                                                           |    0.690 | 0.687-0.697 |      157.273 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    0.691 | 0.688-0.693 |      156.930 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    0.903 | 0.895-0.905 |      157.633 |   20 |                 1 |            4.257 |
| chr20:10000000-20000000 |        0 | htslib                     |                                                                                                   |    0.324 | 0.317-0.334 |      155.371 |   22 |                NA |          138.201 |
| chr20:10000000-20000000 |        0 | defaults                   |                                                                                                   |    0.646 | 0.641-0.648 |      156.953 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | force_download             | force_download = true                                                                             |    1.003 | 1.001-1.004 |    1,023.449 |    8 |                NA |        1,730.523 |
| chr20:10000000-20000000 |        0 | no_external_file_cache     | enable_external_file_cache = false                                                                |    0.652 | 0.633-0.668 |      157.422 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    0.648 | 0.645-0.650 |      157.145 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    0.651 | 0.633-0.657 |      157.312 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    0.641 | 0.634-0.657 |      157.680 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | no_keep_alive              | http_keep_alive = false                                                                           |    0.430 | 0.406-0.434 |      157.453 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    0.440 | 0.431-0.445 |      157.402 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    0.639 | 0.630-0.663 |      160.555 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | htslib                     |                                                                                                   |    0.795 | 0.787-0.800 |      157.391 |   22 |                NA |          151.326 |
| chr20:10000000-20000000 |       20 | defaults                   |                                                                                                   |    4.642 | 4.608-4.669 |      156.941 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | force_download             | force_download = true                                                                             |    1.191 | 1.174-1.205 |    1,022.680 |    8 |                NA |        1,730.523 |
| chr20:10000000-20000000 |       20 | no_external_file_cache     | enable_external_file_cache = false                                                                |    4.684 | 4.649-4.687 |      157.445 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    4.679 | 4.660-4.691 |      157.402 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    4.642 | 4.631-4.709 |      157.309 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    4.682 | 4.651-4.685 |      157.836 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | no_keep_alive              | http_keep_alive = false                                                                           |    4.587 | 4.581-4.590 |      157.320 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    4.600 | 4.555-4.607 |      158.957 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    4.665 | 4.658-4.688 |      157.340 |  128 |                 1 |          112.257 |

## Block size

The same two regions on one thread, on both routes, with
`duckhts_remote_block_size` at each value below, 3 runs each. `auto` is
0, which for a region query is 1 MiB. The environment variable
`DUCKHTS_BENCH_EXTRA_R_LIBS` names further R libraries that each hold
another version of the duckdb package; the same extension binary is
loaded into each.

| DuckDB | httpfs     |
|:-------|:-----------|
| v2.0.0 | 77392a7207 |
| v1.5.6 | 4bc690d    |

Median seconds:

| DuckDB | region                  | delay ms | route  |  auto | 0.125 MiB | 1 MiB | 4 MiB | 16 MiB | 64 MiB |
|:-------|:------------------------|---------:|:-------|------:|----------:|------:|------:|-------:|-------:|
| v2.0.0 | chr20:10000000-10100000 |        0 | htslib | 0.055 |     0.043 | 0.050 | 0.059 |  0.066 |  0.104 |
| v2.0.0 | chr20:10000000-10100000 |        0 | duckdb | 0.332 |     0.333 | 0.330 | 0.337 |  0.341 |  0.370 |
| v2.0.0 | chr20:10000000-20000000 |        0 | htslib | 0.307 |     0.297 | 0.307 | 0.339 |  0.343 |  0.374 |
| v2.0.0 | chr20:10000000-20000000 |        0 | duckdb | 0.655 |     0.658 | 0.672 | 0.630 |  0.629 |  0.656 |
| v2.0.0 | chr20:10000000-10100000 |       20 | htslib | 0.528 |     0.522 | 0.530 | 0.534 |  0.544 |  0.606 |
| v2.0.0 | chr20:10000000-10100000 |       20 | duckdb | 0.899 |     1.069 | 0.895 | 0.883 |  0.888 |  0.913 |
| v2.0.0 | chr20:10000000-20000000 |       20 | htslib | 0.795 |     0.777 | 0.790 | 0.809 |  0.837 |  0.869 |
| v2.0.0 | chr20:10000000-20000000 |       20 | duckdb | 4.683 |    20.689 | 4.671 | 2.110 |  1.356 |  1.220 |
| v1.5.6 | chr20:10000000-10100000 |        0 | htslib | 0.028 |     0.023 | 0.028 | 0.033 |  0.054 |  0.099 |
| v1.5.6 | chr20:10000000-10100000 |        0 | duckdb | 0.162 |     0.167 | 0.171 | 0.160 |  0.168 |  0.280 |
| v1.5.6 | chr20:10000000-20000000 |        0 | htslib | 0.306 |     0.293 | 0.307 | 0.319 |  0.336 |  0.375 |
| v1.5.6 | chr20:10000000-20000000 |        0 | duckdb | 0.564 |     0.573 | 0.566 | 0.485 |  0.483 |  0.588 |
| v1.5.6 | chr20:10000000-10100000 |       20 | htslib | 0.358 |     0.355 | 0.360 | 0.366 |  0.381 |  0.444 |
| v1.5.6 | chr20:10000000-10100000 |       20 | duckdb | 0.458 |     0.453 | 0.452 | 0.439 |  0.452 |  0.579 |
| v1.5.6 | chr20:10000000-20000000 |       20 | htslib | 0.632 |     0.625 | 0.630 | 0.652 |  0.674 |  0.717 |
| v1.5.6 | chr20:10000000-20000000 |       20 | duckdb | 1.060 |     1.050 | 1.047 | 1.749 |  0.952 |  0.890 |

GET requests (median):

| DuckDB | region                  | delay ms | route  | auto | 0.125 MiB | 1 MiB | 4 MiB | 16 MiB | 64 MiB |
|:-------|:------------------------|---------:|:-------|-----:|----------:|------:|------:|-------:|-------:|
| v2.0.0 | chr20:10000000-10100000 |        0 | htslib |   22 |        22 |    22 |    22 |     22 |     22 |
| v2.0.0 | chr20:10000000-10100000 |        0 | duckdb |   20 |        29 |    20 |    19 |     19 |     19 |
| v2.0.0 | chr20:10000000-20000000 |        0 | htslib |   22 |        22 |    22 |    22 |     22 |     22 |
| v2.0.0 | chr20:10000000-20000000 |        0 | duckdb |  128 |       893 |   128 |    46 |     25 |     20 |
| v2.0.0 | chr20:10000000-10100000 |       20 | htslib |   22 |        22 |    22 |    22 |     22 |     22 |
| v2.0.0 | chr20:10000000-10100000 |       20 | duckdb |   20 |        29 |    20 |    19 |     19 |     19 |
| v2.0.0 | chr20:10000000-20000000 |       20 | htslib |   22 |        22 |    22 |    22 |     22 |     22 |
| v2.0.0 | chr20:10000000-20000000 |       20 | duckdb |  128 |       893 |   128 |    46 |     25 |     20 |
| v1.5.6 | chr20:10000000-10100000 |        0 | htslib |   15 |        15 |    15 |    15 |     15 |     15 |
| v1.5.6 | chr20:10000000-10100000 |        0 | duckdb |   10 |        10 |    10 |     9 |      9 |      9 |
| v1.5.6 | chr20:10000000-20000000 |        0 | htslib |   15 |        15 |    15 |    15 |     15 |     15 |
| v1.5.6 | chr20:10000000-20000000 |        0 | duckdb |   16 |        16 |    16 |    36 |     15 |     10 |
| v1.5.6 | chr20:10000000-10100000 |       20 | htslib |   15 |        15 |    15 |    15 |     15 |     15 |
| v1.5.6 | chr20:10000000-10100000 |       20 | duckdb |   10 |        10 |    10 |     9 |      9 |      9 |
| v1.5.6 | chr20:10000000-20000000 |       20 | htslib |   15 |        15 |    15 |    15 |     15 |     15 |
| v1.5.6 | chr20:10000000-20000000 |       20 | duckdb |   16 |        16 |    16 |    36 |     15 |     10 |

Body MiB the server wrote (median):

| DuckDB | region                  | delay ms | route  |  auto | 0.125 MiB | 1 MiB | 4 MiB | 16 MiB | 64 MiB |
|:-------|:------------------------|---------:|:-------|------:|----------:|------:|------:|-------:|-------:|
| v2.0.0 | chr20:10000000-10100000 |        0 | htslib |  26.7 |      24.5 |  27.0 |  29.6 |   57.5 |  153.0 |
| v2.0.0 | chr20:10000000-10100000 |        0 | duckdb |   4.3 |       2.8 |   4.3 |   9.3 |   33.3 |  129.3 |
| v2.0.0 | chr20:10000000-20000000 |        0 | htslib | 137.8 |     133.1 | 136.3 | 143.8 |  155.0 |  222.9 |
| v2.0.0 | chr20:10000000-20000000 |        0 | duckdb | 112.3 |     110.8 | 112.3 | 117.3 |  129.3 |  193.3 |
| v2.0.0 | chr20:10000000-10100000 |       20 | htslib |  43.3 |      40.1 |  43.6 |  42.1 |   68.2 |  162.8 |
| v2.0.0 | chr20:10000000-10100000 |       20 | duckdb |   4.3 |       2.8 |   4.3 |   9.3 |   33.3 |  129.3 |
| v2.0.0 | chr20:10000000-20000000 |       20 | htslib | 151.1 |     149.4 | 151.5 | 158.6 |  162.3 |  237.9 |
| v2.0.0 | chr20:10000000-20000000 |       20 | duckdb | 112.3 |     110.8 | 112.3 | 117.3 |  129.3 |  193.3 |
| v1.5.6 | chr20:10000000-10100000 |        0 | htslib |  21.0 |      18.5 |  21.6 |  25.0 |   52.0 |  147.6 |
| v1.5.6 | chr20:10000000-10100000 |        0 | duckdb |   7.4 |       7.4 |   7.4 |  11.4 |   35.4 |  131.4 |
| v1.5.6 | chr20:10000000-20000000 |        0 | htslib | 131.8 |     128.7 | 131.0 | 138.5 |  150.3 |  218.1 |
| v1.5.6 | chr20:10000000-20000000 |        0 | duckdb | 131.4 |     131.4 | 131.4 | 119.4 |  131.4 |  195.4 |
| v1.5.6 | chr20:10000000-10100000 |       20 | htslib |  31.6 |      27.7 |  29.6 |  35.0 |   57.3 |  155.7 |
| v1.5.6 | chr20:10000000-10100000 |       20 | duckdb |   7.4 |       7.4 |   7.4 |  11.4 |   35.4 |  131.4 |
| v1.5.6 | chr20:10000000-20000000 |       20 | htslib | 141.5 |     138.5 | 140.3 | 150.0 |  156.6 |  226.1 |
| v1.5.6 | chr20:10000000-20000000 |       20 | duckdb | 131.4 |     131.4 | 131.4 | 119.4 |  131.4 |  195.4 |

Peak resident memory across these runs was at most 335 MiB, at a block
size of 64 MiB.

## Scaling

Time exponent for each doubling of the input, `log2(t(2n) / t(n))`, from
the median times above. For region queries the input is the region width
from 2.5 Mb to 5 Mb to 10 Mb; for scans it is the file from 1× to 2× to
4×.

Region width, 1× is the 2.5 Mb region:

| route  | delay ms | threads | records 1x | seconds 1x | exponent 1x to 2x | exponent 2x to 4x | record ratio 1x to 2x | record ratio 2x to 4x |
|:-------|---------:|--------:|-----------:|-----------:|------------------:|------------------:|----------------------:|----------------------:|
| duckdb |        0 |       1 |    582,233 |      0.412 |             0.238 |             0.411 |                 2.009 |                 2.024 |
| htslib |        0 |       1 |    582,233 |      0.108 |             0.729 |             0.820 |                 2.009 |                 2.024 |
| duckdb |       20 |       1 |    582,233 |      1.775 |             0.598 |             0.797 |                 2.009 |                 2.024 |
| htslib |       20 |       1 |    582,233 |      0.593 |             0.157 |             0.252 |                 2.009 |                 2.024 |

File size:

| route  | delay ms | threads | records 1x | seconds 1x | exponent 1x to 2x | exponent 2x to 4x | record ratio 1x to 2x | record ratio 2x to 4x |
|:-------|---------:|--------:|-----------:|-----------:|------------------:|------------------:|----------------------:|----------------------:|
| duckdb |        0 |       1 | 16,816,430 |      2.711 |             0.785 |             0.862 |                     2 |                     2 |
| htslib |        0 |       1 | 16,816,430 |      2.053 |             0.910 |             0.954 |                     2 |                     2 |
| duckdb |        0 |       4 | 16,816,430 |      1.306 |             0.568 |             0.734 |                     2 |                     2 |
| htslib |        0 |       4 | 16,816,430 |      0.851 |             0.820 |             0.973 |                     2 |                     2 |

`STYLE.md` asks for a 1× run of at least 5 s on one thread before a
timing verdict. The longest single-thread 1× run here takes 2.71 s, so
these exponents describe the measured runs and carry no verdict, and
none may be read as a statement about larger inputs than those measured.

## What the measurements show

- **Same rows.** All 122 runs returned the samtools record count and one
  hash per input, on both routes.
- **Request pattern.** On the 10 Mb region the htslib route sent 22 GET
  requests, 5 of them open-ended; the DuckDB route sent 128, all closed
  ranges, the largest 1 MiB. Across the four region widths the DuckDB
  route’s count rises with the bytes read and the htslib route’s count
  stays the same; in the whole-file scans both rise, the DuckDB route’s
  far more. No closed range in any run exceeded 1 MiB. htslib asks its
  backend for at most one hFILE block at a time (`refill_buffer` in
  `third_party/htslib/hfile.c`), and every scan here, the single-thread
  ones included, runs under the 1 MiB profile: a full scan of an indexed
  BAM is cut into pieces (`src/bam_reader.c`, `is_parallel`).
- **Block size.** On the 10 Mb region with 20 ms per request:
  - DuckDB v2.0.0: the DuckDB route took 4.68 s with 128 GETs at the
    automatic size and 1.22 s with 20 GETs at its fastest size, 64 MiB;
    the htslib route took 0.80 s with 22 GETs at the automatic size.
  - DuckDB v1.5.6: the DuckDB route took 1.06 s with 16 GETs at the
    automatic size and 0.89 s with 10 GETs at its fastest size, 64 MiB;
    the htslib route took 0.63 s with 15 GETs at the automatic size.
- **httpfs settings.** `no_external_file_cache`, `grid_request_sizing`,
  `cache_block_16mib`, `grid_and_cache_block_16mib`, `no_keep_alive`,
  `no_connection_caching`, `curl_client` left the GET count of both
  regions equal to the defaults. `force_download` replaced the ranges
  with 8 requests for the whole file: the 10 Mb region with 20 ms per
  request took 1.19 s against 4.64 s at the defaults, at a peak of 1023
  MiB resident and 1731 MiB sent for a 865 MiB file.
- **Per-request cost.** With no added delay the 10 Mb region took 0.32 s
  through htslib and 0.65 s through DuckDB. With 20 ms added per request
  it took 0.79 s and 4.67 s. The extra time per GET is 21 ms on the
  htslib route and 31 ms on the DuckDB route.
- **Connection settings.** With the same GET count,
  `http_keep_alive = false` took 0.07 s on the 100 kb region with no
  added delay against 0.33 s at the defaults, and 0.43 s against 0.65 s
  on the 10 Mb region. With 20 ms per request the 10 Mb region took 4.59
  s against 4.64 s. `httpfs_connection_caching = false` gave times
  within the run-to-run range of those. This report does not establish
  why reusing connections costs time against this server.
- **Memory.** In the route comparison, peak resident memory stayed
  between 147 and 173 MiB, under every ceiling, and no run spilled.

## Not measured

- The block size was varied for region queries on one thread only, not
  for whole-file scans or several threads. The BGZF cache size was not
  varied.
- httpfs settings beyond those in the table, and combinations of them.
  Whether the external file cache serves reads made through the C file
  API was not examined.
- The route comparison, the scans and the httpfs settings ran on the
  DuckDB version of the first table only.
- Real object stores, TLS, bandwidth limits, and connections that fail
  or stall.
- CRAM, SAM, and any reader other than `read_bam`.
- Writing: the DuckDB route is read-only.
