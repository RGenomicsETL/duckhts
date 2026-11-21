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
for. The tuning is fixed in the source; neither a setting nor a function
argument changes it. htslib’s libcurl backend has further controls that
this benchmark leaves at their defaults: the retry and low-speed
environment variables (`HTS_RETRY_MAX`, `HTS_RETRY_DELAY`,
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
| source revision             | 517c7436cac1a7aa5004ea2eb925e5f22ff29de7                                                |
| source tree                 | clean                                                                                   |
| extension md5               | cd70bf3f4669192065e3d76a2f29a6c9                                                        |
| htslib                      | 1.24                                                                                    |
| R                           | R version 4.6.2 (2026-10-29)                                                            |
| DuckDB R package            | 2.0.0                                                                                   |
| httpfs                      | installed build for this DuckDB version, default settings                               |
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
154 MiB. The ceiling for a run is that overhead plus three times the I/O
state of its open handles: one hFILE block, one BGZF block cache and one
32 MiB read-ahead buffer per handle. A region query and a single-thread
scan hold one handle; a scan on four threads holds four. The query
itself keeps no rows: it aggregates as it reads.

## Region queries

One thread, the 1× file, regions that start at chr20:10,000,000.

| region                  | scale | delay ms | threads | route  |   records | median s | range s     | peak RSS MiB | ceiling MiB | GETs | open-ended | closed ranges | largest range MiB | body MiB written |
|:------------------------|:------|---------:|--------:|:-------|----------:|---------:|:------------|-------------:|------------:|-----:|-----------:|--------------:|------------------:|-----------------:|
| chr20:10000000-10100000 | 1x    |        0 |       1 | htslib |    23,544 |    0.051 | 0.043-0.071 |      153.906 |     444.609 |   22 |          5 |             0 |                NA |           23.764 |
| chr20:10000000-10100000 | 1x    |        0 |       1 | duckdb |    23,544 |    0.327 | 0.306-0.334 |      156.840 |     444.609 |   20 |          0 |            20 |                 1 |            4.257 |
| chr20:10000000-12500000 | 1x    |        0 |       1 | htslib |   582,233 |    0.112 | 0.108-0.120 |      153.738 |     444.609 |   22 |          5 |             0 |                NA |           51.764 |
| chr20:10000000-12500000 | 1x    |        0 |       1 | duckdb |   582,233 |    0.412 | 0.404-0.420 |      156.656 |     444.609 |   45 |          0 |            45 |                 1 |           29.257 |
| chr20:10000000-15000000 | 1x    |        0 |       1 | htslib | 1,169,601 |    0.175 | 0.170-0.194 |      154.051 |     444.609 |   22 |          5 |             0 |                NA |           78.201 |
| chr20:10000000-15000000 | 1x    |        0 |       1 | duckdb | 1,169,601 |    0.486 | 0.478-0.505 |      156.973 |     444.609 |   72 |          0 |            72 |                 1 |           56.257 |
| chr20:10000000-20000000 | 1x    |        0 |       1 | htslib | 2,367,216 |    0.322 | 0.311-0.331 |      155.371 |     444.609 |   22 |          5 |             0 |                NA |          134.639 |
| chr20:10000000-20000000 | 1x    |        0 |       1 | duckdb | 2,367,216 |    0.638 | 0.579-0.647 |      158.789 |     444.609 |  128 |          0 |           128 |                 1 |          112.257 |
| chr20:10000000-10100000 | 1x    |       20 |       1 | htslib |    23,544 |    0.530 | 0.528-0.540 |      155.234 |     444.609 |   22 |          5 |             0 |                NA |           38.201 |
| chr20:10000000-10100000 | 1x    |       20 |       1 | duckdb |    23,544 |    0.901 | 0.898-0.902 |      156.812 |     444.609 |   20 |          0 |            20 |                 1 |            4.257 |
| chr20:10000000-12500000 | 1x    |       20 |       1 | htslib |   582,233 |    0.590 | 0.582-0.603 |      155.543 |     444.609 |   22 |          5 |             0 |                NA |           60.576 |
| chr20:10000000-12500000 | 1x    |       20 |       1 | duckdb |   582,233 |    1.789 | 1.768-1.798 |      157.098 |     444.609 |   45 |          0 |            45 |                 1 |           29.257 |
| chr20:10000000-15000000 | 1x    |       20 |       1 | htslib | 1,169,601 |    0.655 | 0.650-0.669 |      156.785 |     444.609 |   22 |          5 |             0 |                NA |           93.014 |
| chr20:10000000-15000000 | 1x    |       20 |       1 | duckdb | 1,169,601 |    2.735 | 2.703-2.739 |      157.066 |     444.609 |   72 |          0 |            72 |                 1 |           56.257 |
| chr20:10000000-20000000 | 1x    |       20 |       1 | htslib | 2,367,216 |    0.795 | 0.786-0.814 |      155.504 |     444.609 |   22 |          5 |             0 |                NA |          147.639 |
| chr20:10000000-20000000 | 1x    |       20 |       1 | duckdb | 2,367,216 |    4.684 | 4.607-4.701 |      158.395 |     444.609 |  128 |          0 |           128 |                 1 |          112.257 |

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
| whole file | 1x    |        0 |       1 | htslib | 16,816,430 |    2.083 | 2.074-2.132 |      155.684 |   1,029.609 |    33 |         16 |             0 |                NA |          994.109 |
| whole file | 1x    |        0 |       1 | duckdb | 16,816,430 |    2.733 | 2.720-2.755 |      155.449 |   1,029.609 |   899 |          0 |           899 |                 1 |          882.228 |
| whole file | 2x    |        0 |       1 | htslib | 33,632,860 |    3.948 | 3.879-3.971 |      154.152 |   1,029.609 |    42 |         25 |             0 |                NA |        1,322.057 |
| whole file | 2x    |        0 |       1 | duckdb | 33,632,860 |    4.775 | 4.727-4.961 |      155.262 |   1,029.609 | 1,152 |          0 |         1,152 |                 1 |        1,134.382 |
| whole file | 4x    |        0 |       1 | htslib | 67,265,720 |    7.604 | 7.571-7.837 |      154.207 |   1,029.609 |    44 |         27 |             0 |                NA |        1,668.958 |
| whole file | 4x    |        0 |       1 | duckdb | 67,265,720 |    8.690 | 8.642-8.715 |      155.492 |   1,029.609 | 1,481 |          0 |         1,481 |                 1 |        1,463.623 |
| whole file | 1x    |        0 |       4 | htslib | 16,816,430 |    0.864 | 0.842-0.871 |      168.688 |   1,317.609 |    96 |         64 |             0 |                NA |        1,402.972 |
| whole file | 1x    |        0 |       4 | duckdb | 16,816,430 |    1.306 | 1.214-1.329 |      168.355 |   1,317.609 |   930 |          0 |           930 |                 1 |          898.801 |
| whole file | 2x    |        0 |       4 | htslib | 33,632,860 |    1.593 | 1.452-1.626 |      169.348 |   1,317.609 |   112 |         80 |             0 |                NA |        1,818.545 |
| whole file | 2x    |        0 |       4 | duckdb | 33,632,860 |    1.908 | 1.853-1.951 |      168.723 |   1,317.609 | 1,204 |          0 |         1,204 |                 1 |        1,173.183 |
| whole file | 4x    |        0 |       4 | htslib | 67,265,720 |    2.931 | 2.927-2.997 |      168.922 |   1,317.609 |   129 |         97 |             0 |                NA |        2,306.208 |
| whole file | 4x    |        0 |       4 | duckdb | 67,265,720 |    3.228 | 3.222-3.311 |      168.285 |   1,317.609 | 1,536 |          0 |         1,536 |                 1 |        1,505.026 |
| whole file | 1x    |       NA |       1 | htslib | 16,816,430 |    1.954 | 1.923-1.956 |      147.633 |   1,029.609 |     0 |          0 |             0 |                NA |            0.000 |
| whole file | 1x    |       NA |       1 | duckdb | 16,816,430 |    1.946 | 1.926-1.954 |      146.934 |   1,029.609 |     0 |          0 |             0 |                NA |            0.000 |

## httpfs settings

The 100 kb and 10 Mb regions of the 1× file on one thread, on the DuckDB
route under httpfs settings that each differ from the defaults in one
respect, 3 runs each. The `htslib` rows are the htslib route at its
defaults, for reference.

| region                  | delay ms | setting                    | SET                                                                                               | median s | range s     | peak RSS MiB | GETs | largest range MiB | body MiB written |
|:------------------------|---------:|:---------------------------|:--------------------------------------------------------------------------------------------------|---------:|:------------|-------------:|-----:|------------------:|-----------------:|
| chr20:10000000-10100000 |        0 | htslib                     |                                                                                                   |    0.065 | 0.053-0.065 |      154.223 |   22 |                NA |           27.451 |
| chr20:10000000-10100000 |        0 | defaults                   |                                                                                                   |    0.321 | 0.272-0.329 |      156.879 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | force_download             | force_download = true                                                                             |    0.713 | 0.712-0.721 |    1,022.516 |    8 |                NA |        1,730.523 |
| chr20:10000000-10100000 |        0 | no_external_file_cache     | enable_external_file_cache = false                                                                |    0.325 | 0.322-0.334 |      156.996 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    0.334 | 0.324-0.334 |      157.258 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    0.336 | 0.325-0.336 |      157.062 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    0.332 | 0.322-0.332 |      157.410 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | no_keep_alive              | http_keep_alive = false                                                                           |    0.071 | 0.057-0.077 |      156.176 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    0.057 | 0.054-0.059 |      155.383 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |        0 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    0.331 | 0.327-0.336 |      157.207 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | htslib                     |                                                                                                   |    0.533 | 0.527-0.540 |      155.203 |   22 |                NA |           42.701 |
| chr20:10000000-10100000 |       20 | defaults                   |                                                                                                   |    0.903 | 0.894-0.903 |      159.512 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | force_download             | force_download = true                                                                             |    0.911 | 0.909-0.920 |    1,023.465 |    8 |                NA |        1,730.523 |
| chr20:10000000-10100000 |       20 | no_external_file_cache     | enable_external_file_cache = false                                                                |    0.899 | 0.897-0.900 |      157.133 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    0.901 | 0.898-0.905 |      157.000 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    0.898 | 0.896-0.905 |      156.941 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    0.898 | 0.894-0.910 |      157.043 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | no_keep_alive              | http_keep_alive = false                                                                           |    0.691 | 0.687-0.691 |      157.332 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    0.688 | 0.688-0.693 |      156.648 |   20 |                 1 |            4.257 |
| chr20:10000000-10100000 |       20 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    0.900 | 0.900-0.900 |      157.008 |   20 |                 1 |            4.257 |
| chr20:10000000-20000000 |        0 | htslib                     |                                                                                                   |    0.319 | 0.317-0.325 |      155.289 |   22 |                NA |          137.701 |
| chr20:10000000-20000000 |        0 | defaults                   |                                                                                                   |    0.645 | 0.582-0.649 |      156.867 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | force_download             | force_download = true                                                                             |    0.967 | 0.961-0.986 |    1,022.660 |    8 |                NA |        1,730.523 |
| chr20:10000000-20000000 |        0 | no_external_file_cache     | enable_external_file_cache = false                                                                |    0.652 | 0.638-0.660 |      156.938 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    0.668 | 0.656-0.674 |      157.086 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    0.656 | 0.645-0.658 |      157.277 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    0.644 | 0.641-0.654 |      157.430 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | no_keep_alive              | http_keep_alive = false                                                                           |    0.443 | 0.438-0.445 |      157.043 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    0.443 | 0.434-0.465 |      158.680 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |        0 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    0.656 | 0.643-0.656 |      156.980 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | htslib                     |                                                                                                   |    0.794 | 0.791-0.797 |      155.805 |   22 |                NA |          148.389 |
| chr20:10000000-20000000 |       20 | defaults                   |                                                                                                   |    4.639 | 4.638-4.661 |      157.109 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | force_download             | force_download = true                                                                             |    1.169 | 1.158-1.186 |    1,022.906 |    8 |                NA |        1,730.523 |
| chr20:10000000-20000000 |       20 | no_external_file_cache     | enable_external_file_cache = false                                                                |    4.689 | 4.683-4.713 |      157.078 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | grid_request_sizing        | external_file_cache_request_sizing = ‘GRID’                                                       |    4.673 | 4.630-4.678 |      157.156 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | cache_block_16mib          | external_file_cache_remote_max_block_size = 16777216                                              |    4.703 | 4.694-4.705 |      157.352 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | grid_and_cache_block_16mib | external_file_cache_request_sizing = ‘GRID’; external_file_cache_remote_max_block_size = 16777216 |    4.697 | 4.675-4.702 |      157.297 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | no_keep_alive              | http_keep_alive = false                                                                           |    4.594 | 4.575-4.601 |      157.242 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | no_connection_caching      | httpfs_connection_caching = false                                                                 |    4.574 | 4.574-4.579 |      157.234 |  128 |                 1 |          112.257 |
| chr20:10000000-20000000 |       20 | curl_client                | httpfs_client_implementation = ‘curl’                                                             |    4.681 | 4.670-4.707 |      160.223 |  128 |                 1 |          112.257 |

## Scaling

Time exponent for each doubling of the input, `log2(t(2n) / t(n))`, from
the median times above. For region queries the input is the region width
from 2.5 Mb to 5 Mb to 10 Mb; for scans it is the file from 1× to 2× to
4×.

Region width, 1× is the 2.5 Mb region:

| route  | delay ms | threads | records 1x | seconds 1x | exponent 1x to 2x | exponent 2x to 4x | record ratio 1x to 2x | record ratio 2x to 4x |
|:-------|---------:|--------:|-----------:|-----------:|------------------:|------------------:|----------------------:|----------------------:|
| duckdb |        0 |       1 |    582,233 |      0.412 |             0.238 |             0.393 |                 2.009 |                 2.024 |
| htslib |        0 |       1 |    582,233 |      0.112 |             0.644 |             0.880 |                 2.009 |                 2.024 |
| duckdb |       20 |       1 |    582,233 |      1.789 |             0.612 |             0.776 |                 2.009 |                 2.024 |
| htslib |       20 |       1 |    582,233 |      0.590 |             0.151 |             0.279 |                 2.009 |                 2.024 |

File size:

| route  | delay ms | threads | records 1x | seconds 1x | exponent 1x to 2x | exponent 2x to 4x | record ratio 1x to 2x | record ratio 2x to 4x |
|:-------|---------:|--------:|-----------:|-----------:|------------------:|------------------:|----------------------:|----------------------:|
| duckdb |        0 |       1 | 16,816,430 |      2.733 |             0.805 |             0.864 |                     2 |                     2 |
| htslib |        0 |       1 | 16,816,430 |      2.083 |             0.922 |             0.946 |                     2 |                     2 |
| duckdb |        0 |       4 | 16,816,430 |      1.306 |             0.547 |             0.759 |                     2 |                     2 |
| htslib |        0 |       4 | 16,816,430 |      0.864 |             0.883 |             0.880 |                     2 |                     2 |

`STYLE.md` asks for a 1× run of at least 5 s on one thread before a
timing verdict. The longest single-thread 1× run here takes 2.73 s, so
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
  BAM is cut into pieces (`src/bam_reader.c`, `is_parallel`). The block
  size was not varied, so its effect on the request count is read from
  the code, not measured.
- **httpfs settings.** `no_external_file_cache`, `grid_request_sizing`,
  `cache_block_16mib`, `grid_and_cache_block_16mib`, `no_keep_alive`,
  `no_connection_caching`, `curl_client` left the GET count of both
  regions equal to the defaults. `force_download` replaced the ranges
  with 8 requests for the whole file: the 10 Mb region with 20 ms per
  request took 1.17 s against 4.64 s at the defaults, at a peak of 1023
  MiB resident and 1731 MiB sent for a 865 MiB file.
- **Per-request cost.** With no added delay the 10 Mb region took 0.32 s
  through htslib and 0.64 s through DuckDB. With 20 ms added per request
  it took 0.79 s and 4.68 s. The extra time per GET is 21 ms on the
  htslib route and 32 ms on the DuckDB route.
- **Connection settings.** With the same GET count,
  `http_keep_alive = false` took 0.07 s on the 100 kb region with no
  added delay against 0.32 s at the defaults, and 0.44 s against 0.64 s
  on the 10 Mb region. With 20 ms per request the 10 Mb region took 4.59
  s against 4.64 s. `httpfs_connection_caching = false` gave times
  within the run-to-run range of those. This report does not establish
  why reusing connections costs time against this server.
- **Memory.** In the route comparison, peak resident memory stayed
  between 146 and 169 MiB, under every ceiling, and no run spilled.

## Not measured

- The hFILE block size and BGZF cache size were not varied: DuckHTS
  fixes them in `src/include/hts_io_tuning.h`. A sweep of the block size
  on both routes, with the scan profile each handle received recorded,
  is the missing measurement that would say how much of the DuckDB
  route’s request count is a tuning choice.
- httpfs settings beyond those in the table, and combinations of them.
  Whether the external file cache serves reads made through the C file
  API was not examined.
- Only the DuckDB version in the table above was run.
- Real object stores, TLS, bandwidth limits, and connections that fail
  or stall.
- CRAM, SAM, and any reader other than `read_bam`.
- Writing: the DuckDB route is read-only.
