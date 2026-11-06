Counting BED fragments by tile and cell barcode in SQL
================

<!-- benchmark_fragment_bin_counts.md is generated from benchmark_fragment_bin_counts.Rmd. -->

@ttriche happy that the package is useful to you!

On the bin count question, for these things i usually start from composing over some sql that does the equivalent work, and if it is not satisfactory add a native kernel. Some tests of this approach are in [duckseq](https://github.com/RGenomicsETL/duckseq/tree/main). Its [AIE demo](https://rgenomicsetl.github.io/duckseq/aie/) just got a [rewrite of its preparation SQL](https://github.com/RGenomicsETL/duckseq/commit/d53d66df415fb045e3473f11afabf267590df417): 4M real PBMC records at 4 threads went from 13.9 s to 4.3 s, against 5.3 s for Gravlax. The two do not store the same thing: the Gravlax archive is 20 times smaller, and it still wins on one thread.

Indeed i thought at one point about adding a fragment count mode to `bam_bin_counts`, if this is useful i can add this functionality and the blacklisting: for the native kernel cgranges should make the blacklisting relatively fast.

For now, the SQL below should work. I am guessing the application from the file names in [\#356](https://github.com/RGenomicsETL/duckhts/issues/356): single-cell CUT&Tag style `chrom, start, end, name` fragments, `name` is the cell barcode, 50 kb hg38 tiles, and a sparse tile by cell matrix at the end. The data is simulated in that shape, it is not a real experiment.

<details>
<summary>
Simulated input
</summary>

``` r
library(data.table)

hg38 <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
  chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
  chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
  chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
  chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415
)
tile_width <- 50000L
n_cells <- 5000L
n_background <- 15000L

set.seed(32026)
random_barcode <- function(n) {
  bases <- matrix(sample(c("A", "C", "G", "T"), n * 16L, TRUE), nrow = n)
  do.call(paste0, as.data.frame(bases))
}
barcodes <- unique(random_barcode(n_cells + n_background))
cell_barcodes <- barcodes[seq_len(n_cells)]
barcode_prob <- c(
  rep(0.9 / n_cells, n_cells),
  rep(0.1 / (length(barcodes) - n_cells), length(barcodes) - n_cells)
)

# Broad marked domains: 70% of fragments fall in about 3,000 domains.
domains <- data.table(chrom = sample(names(hg38), 3000L, TRUE, prob = hg38))
domains[, width := as.integer(runif(.N, 1e5, 5e5))]
domains[, start := as.integer(runif(.N) * (hg38[chrom] - width - 1000))]

simulate_fragments <- function(n, path) {
  in_domain <- runif(n) < 0.7
  domain <- sample(nrow(domains), n, TRUE, prob = domains$width)
  genome_chrom <- sample(names(hg38), n, TRUE, prob = hg38)
  chrom <- ifelse(in_domain, domains$chrom[domain], genome_chrom)
  start <- ifelse(
    in_domain,
    domains$start[domain] + runif(n) * domains$width[domain],
    runif(n) * (hg38[genome_chrom] - 1000)
  )
  nucleosomes <- ifelse(runif(n) < 0.7, 180, 360)
  width <- pmax(30L, as.integer(rnorm(n, nucleosomes, 40)))
  fragments <- data.table(
    chrom = factor(chrom, levels = names(hg38)),
    start = as.integer(start),
    end = as.integer(start) + width,
    name = sample(barcodes, n, TRUE, prob = barcode_prob)
  )
  setorder(fragments, chrom, start, end)
  fwrite(fragments, path, sep = "\t", col.names = FALSE)
}

sizes <- c("4M" = 4e6, "8M" = 8e6, "16M" = 16e6)
fragment_files <- sprintf("fragments_%s.bed.gz", names(sizes))
names(fragment_files) <- names(sizes)
for (size in names(sizes)) {
  simulate_fragments(sizes[[size]], fragment_files[[size]])
}
writeLines(cell_barcodes, "cells.txt")

blacklist <- data.table(chrom = sample(names(hg38), 640L, TRUE, prob = hg38))
blacklist[, start := as.integer(runif(.N) * (hg38[chrom] - 2e5))]
blacklist[, end := start + 100L + as.integer(rexp(.N, 1 / 50000))]
fwrite(blacklist, "blacklist.bed", sep = "\t", col.names = FALSE)

tiles <- rbindlist(lapply(names(hg38), function(chrom) {
  start <- seq(0, hg38[[chrom]] - 1, by = tile_width)
  data.table(chrom, start = start, end = pmin(start + tile_width, hg38[[chrom]]))
}))
fwrite(tiles, "tiles.bed", sep = "\t", col.names = FALSE, scipen = 50)
```

</details>

# Set up

``` r
library(Rduckhts)
library(DBI)

con <- rduckhts_connect()
invisible(dbExecute(con, "SET threads = 4"))
knitr::opts_chunk$set(connection = "con", max.print = 4)
```

`fragments` is a view over the file, so nothing is loaded yet. `end` is reserved and needs quotes. All coordinates stay 0-based half-open.

``` sql
CREATE TEMP VIEW fragments AS
SELECT chrom, start, "end", name
FROM read_bed('fragments_4M.bed.gz')
```

``` sql
CREATE TEMP TABLE cells AS
SELECT barcode, row_number() OVER (ORDER BY barcode) AS cell_id
FROM read_csv('cells.txt', header = false, names = ['barcode'])
```

``` sql
CREATE TEMP TABLE blacklist AS
SELECT chrom, start, "end"
FROM read_bed('blacklist.bed')
```

``` sql
CREATE TEMP TABLE tiles AS
SELECT chrom, start, "end", row_number() OVER (ORDER BY chrom, start) AS tile_id
FROM read_bed('tiles.bed')
```

``` sql
SELECT * FROM duckhts_cgranges_from_table(
  'blacklist_idx', 'blacklist', 'chrom', 'start', 'end'
)
```

<div class="knitsql-table">

| indexed |
|:--------|
| TRUE    |

1 records

</div>

``` sql
SELECT * FROM duckhts_cgranges_from_table_labelled(
  'tiles_idx', 'tiles', 'chrom', 'start', 'end', 'tile_id'
)
```

<div class="knitsql-table">

| indexed |
|:--------|
| TRUE    |

1 records

</div>

# Counts

**One bin per fragment.** Fixed tiles are arithmetic, here on the fragment midpoint.

``` sql
SELECT
  chrom,
  ((start + "end") // 2) // 50000 AS bin,
  name AS barcode,
  count(*) AS n
FROM fragments
WHERE name IN (SELECT barcode FROM cells)
  AND NOT duckhts_cgranges_has_overlap('blacklist_idx', chrom, start, "end")
GROUP BY ALL
```

<div class="knitsql-table">

| chrom |  bin | barcode          |   n |
|:------|-----:|:-----------------|----:|
| chr1  | 4056 | TGGTCGAATGATAACA |   1 |
| chr1  | 4059 | GAGCTGACTTGATCAA |   1 |
| chr1  | 4061 | CCGGGTGTAACCTAGC |   1 |
| chr1  | 4108 | ATGGTGTCAAGCGGCA |   1 |

Displaying records 1 - 4

</div>

**Every overlapped bin.** A fragment that crosses a tile edge counts in both tiles. The `unnest` must sit in a CTE because DuckDB cannot group on it directly.

``` sql
WITH hits AS (
  SELECT
    chrom,
    unnest(range(start // 50000, ("end" - 1) // 50000 + 1)) AS bin,
    name AS barcode
  FROM fragments
  WHERE name IN (SELECT barcode FROM cells)
    AND NOT duckhts_cgranges_has_overlap('blacklist_idx', chrom, start, "end")
)
SELECT chrom, bin, barcode, count(*) AS n
FROM hits
GROUP BY ALL
```

<div class="knitsql-table">

| chrom |  bin | barcode          |   n |
|:------|-----:|:-----------------|----:|
| chr1  | 3919 | TAGGATTATAGCCATG |   1 |
| chr1  | 3941 | TGTCGGTTGGTAAATT |   1 |
| chr1  | 3961 | GATCGAGCCAAGGGAG |   1 |
| chr1  | 3975 | TAACATCCTGTAAGTA |   1 |

Displaying records 1 - 4

</div>

**Any tile or peak BED, with integer ids.** This is the same count through the tile index, so `tiles` can be irregular. Joining `cells` replaces the barcode string with an integer before the group by, and the ids are the matrix coordinates.

``` sql
WITH hits AS (
  SELECT
    c.cell_id,
    unnest(duckhts_cgranges_overlaps_list(
      'tiles_idx', f.chrom, f.start, f."end"
    )) AS tile
  FROM fragments AS f
  JOIN cells AS c ON f.name = c.barcode
  WHERE NOT duckhts_cgranges_has_overlap('blacklist_idx', f.chrom, f.start, f."end")
)
SELECT tile.label::INTEGER AS tile_id, cell_id, count(*) AS n
FROM hits
GROUP BY ALL
```

<div class="knitsql-table">

| tile_id | cell_id |   n |
|--------:|--------:|----:|
|   61738 |      81 |   1 |
|   61756 |      95 |   1 |
|   24475 |     219 |   1 |
|   15740 |     269 |   1 |

Displaying records 1 - 4

</div>

``` r
counts <- dbGetQuery(con, sql_of("tile_cell"))
tile_by_cell <- Matrix::sparseMatrix(
  i = counts$tile_id,
  j = counts$cell_id,
  x = counts$n,
  dims = c(dbGetQuery(con, "SELECT count(*) AS n FROM tiles")$n, n_cells)
)
dim(tile_by_cell)
#> [1] 61775  5000
length(tile_by_cell@x)
#> [1] 3522749
sum(tile_by_cell)
#> [1] 3576016
```

<details>
<summary>
Check against bedtools and data.table
</summary>

`bedtools intersect -v` does the blacklisting, `data.table` does the barcode filter and the counts. All three queries must match on every row.

``` r
kept <- fread(
  cmd = "bedtools intersect -v -a fragments_4M.bed.gz -b blacklist.bed",
  col.names = c("chrom", "start", "end", "barcode")
)
n_blacklisted <- sizes[["4M"]] - nrow(kept)
kept <- kept[barcode %chin% cell_barcodes]

expected_one_bin <- kept[
  , .(n = as.numeric(.N)),
  by = .(chrom, bin = as.numeric(((start + end) %/% 2L) %/% tile_width), barcode)
]

first_bin <- kept$start %/% tile_width
n_bins <- (kept$end - 1L) %/% tile_width - first_bin + 1L
expanded <- kept[rep(seq_len(.N), n_bins)]
expanded[, bin := as.numeric(rep(first_bin, n_bins) + sequence(n_bins) - 1L)]
expected_all_bins <- expanded[, .(n = as.numeric(.N)), by = .(chrom, bin, barcode)]

observed <- function(label) {
  result <- as.data.table(dbGetQuery(con, sql_of(label)))
  result[, n := as.numeric(n)]
  result
}
one_bin <- observed("one_bin")
one_bin[, bin := as.numeric(bin)]
all_bins <- observed("all_bins")
all_bins[, bin := as.numeric(bin)]

tile_key <- as.data.table(dbGetQuery(con, "SELECT tile_id, chrom, start FROM tiles"))
cell_key <- as.data.table(dbGetQuery(con, "SELECT cell_id, barcode FROM cells"))
tile_cell <- observed("tile_cell")
tile_cell <- tile_key[tile_cell, on = "tile_id"]
tile_cell <- cell_key[tile_cell, on = "cell_id"]
tile_cell <- tile_cell[, .(chrom, bin = as.numeric(start %/% tile_width), barcode, n)]

stopifnot(
  fsetequal(one_bin, expected_one_bin),
  fsetequal(all_bins, expected_all_bins),
  fsetequal(tile_cell, expected_all_bins),
  !fsetequal(expected_one_bin, expected_all_bins),
  n_blacklisted > 0
)

data.frame(
  fragments = sizes[["4M"]],
  blacklisted = n_blacklisted,
  kept_in_cells = nrow(kept),
  one_bin_rows = nrow(expected_one_bin),
  all_bins_rows = nrow(expected_all_bins),
  extra_tile_counts = sum(expected_all_bins$n) - nrow(kept)
)
#>   fragments blacklisted kept_in_cells one_bin_rows all_bins_rows
#> 1     4e+06       45262       3559539      3506782       3522749
#>   extra_tile_counts
#> 1             16477
```

</details>

# Bench

`bench::mark()` runs the chunks above unchanged, each one written to a temp table. `reader_only` is the floor: it touches every barcode and does no other work. The Parquet file is written once from the same fragments.

``` sql
SELECT sum(strlen(name)) AS barcode_bytes
FROM fragments
```

<div class="knitsql-table">

| barcode_bytes |
|--------------:|
|       6.4e+07 |

1 records

</div>

The blacklist can also be plain SQL, so this is `one_bin` without cgranges:

``` sql
SELECT
  f.chrom,
  ((f.start + f."end") // 2) // 50000 AS bin,
  f.name AS barcode,
  count(*) AS n
FROM fragments AS f
WHERE f.name IN (SELECT barcode FROM cells)
  AND NOT EXISTS (
    SELECT 1
    FROM blacklist AS b
    WHERE b.chrom = f.chrom AND b.start < f."end" AND f.start < b."end"
  )
GROUP BY ALL
```

<div class="knitsql-table">

| chrom |  bin | barcode          |   n |
|:------|-----:|:-----------------|----:|
| chr1  | 4203 | GATCTGACGCACGACT |   1 |
| chr1  | 4323 | AAAAAAATCGGAGCCA |   1 |
| chr1  | 4643 | TTCGTCCTATTCGCAC |   1 |
| chr1  | 4807 | ACAGAGACTCGCCCTC |   1 |

Displaying records 1 - 4

</div>

And this is `one_bin` with the integer cell id, to price the barcode string:

``` sql
SELECT
  f.chrom,
  ((f.start + f."end") // 2) // 50000 AS bin,
  c.cell_id,
  count(*) AS n
FROM fragments AS f
JOIN cells AS c ON f.name = c.barcode
WHERE NOT duckhts_cgranges_has_overlap('blacklist_idx', f.chrom, f.start, f."end")
GROUP BY ALL
```

<div class="knitsql-table">

| chrom |  bin | cell_id |   n |
|:------|-----:|--------:|----:|
| chr14 | 1347 |      71 |   1 |
| chr14 |  554 |     153 |   1 |
| chr1  | 2510 |     196 |   1 |
| chr18 |  891 |     218 |   1 |

Displaying records 1 - 4

</div>

``` r
use_source <- function(relation) {
  dbExecute(con, paste("CREATE OR REPLACE TEMP VIEW fragments AS", relation))
}
run <- function(label) {
  dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE result AS", sql_of(label)))
}
time_queries <- function() {
  bench::mark(
    reader_only = run("reader_only"),
    one_bin = run("one_bin"),
    one_bin_cell_id = run("one_bin_cell_id"),
    one_bin_sql_blacklist = run("one_bin_sql_blacklist"),
    all_bins = run("all_bins"),
    tile_cell = run("tile_cell"),
    check = FALSE,
    memory = FALSE,
    filter_gc = FALSE,
    min_time = 0,
    min_iterations = 5,
    max_iterations = 5
  )
}

timings <- list()
convert_seconds <- numeric()
parquet_mib <- numeric()
for (size in names(sizes)) {
  parquet_file <- sprintf("fragments_%s.parquet", size)
  convert_seconds[size] <- system.time(dbExecute(con, sprintf(
    "COPY (SELECT chrom, start, \"end\", name FROM read_bed('%s'))
     TO '%s' (FORMAT parquet, COMPRESSION zstd)",
    fragment_files[[size]], parquet_file
  )))[["elapsed"]]
  parquet_mib[size] <- file.size(parquet_file) / 2^20

  use_source(sprintf(
    "SELECT chrom, start, \"end\", name FROM read_bed('%s')", fragment_files[[size]]
  ))
  timings[[paste("read_bed", size)]] <- time_queries()

  use_source(sprintf("SELECT * FROM read_parquet('%s')", parquet_file))
  timings[[paste("parquet", size)]] <- time_queries()
}

median_seconds <- sapply(timings, function(mark) as.numeric(mark$median))
rownames(median_seconds) <- as.character(timings[[1]]$expression)
knitr::kable(median_seconds, digits = 2)
```

|                       | read_bed 4M | parquet 4M | read_bed 8M | parquet 8M | read_bed 16M | parquet 16M |
|:----------------------|------------:|-----------:|------------:|-----------:|-------------:|------------:|
| reader_only           |        0.67 |       0.01 |        1.30 |       0.02 |         2.61 |        0.03 |
| one_bin               |        1.84 |       0.33 |        3.59 |       0.62 |         7.18 |        1.25 |
| one_bin_cell_id       |        1.61 |       0.29 |        3.28 |       0.53 |         6.46 |        1.06 |
| one_bin_sql_blacklist |        2.31 |       0.56 |        4.72 |       1.07 |         9.41 |        2.01 |
| all_bins              |        1.94 |       0.41 |        3.84 |       0.79 |         7.55 |        1.56 |
| tile_cell             |        2.77 |       1.08 |        5.48 |       2.10 |        10.74 |        4.17 |

Median seconds of 5 runs, 4 DuckDB threads. Time grows linearly with fragments. On 16M fragments:

- `read_bed` runs the whole query on one thread: Parquet is 5.7 times faster for `one_bin`. If you count the same file more than once, convert it once.
- The barcode string is not the main cost: integer ids save 15% from Parquet, so i would not chase a custom hash before fixing the reader.
- cgranges blacklisting is 1.6 times faster than the plain SQL anti join from Parquet.
- Fixed tiles do not need the tile index: `tile_cell` takes 2.7 times as long as `all_bins`, but it accepts any BED and gives matrix ids.

| fragments | gz_MiB | parquet_MiB | convert_once_s |
|:----------|-------:|------------:|---------------:|
| 4M        |   51.4 |        24.9 |           1.43 |
| 8M        |  100.7 |        48.1 |           2.79 |
| 16M       |  197.1 |        92.4 |           5.43 |

Peak RSS of a fresh R process running `tile_cell` on 16M fragments, including R and the connection: 1948 MiB from `read_bed`, 594 MiB from Parquet.

# From R

If you would rather write dplyr, these three routes run the `one_bin` count on the same `con`, so the cgranges index is still there. `end` needs no quotes.

**dbplyr** is the simplest. It passes functions it does not know through to DuckDB.

``` r
library(dplyr)

use_source("SELECT * FROM read_parquet('fragments_4M.parquet')")
#> [1] 0

dbplyr_counts <- tbl(con, "fragments") |>
  semi_join(tbl(con, "cells"), by = c(name = "barcode")) |>
  filter(!duckhts_cgranges_has_overlap("blacklist_idx", chrom, start, end)) |>
  mutate(bin = floor((start + end) / 2 / 50000)) |>
  count(chrom, bin, name)
```

**duckplyr** with `prudence = "stingy"` stops instead of quietly falling back to R. It keeps its operators as macros on its own connection, and that connection cannot load the unsigned extension bundled with Rduckhts. So define the four macros this pipeline needs on `con`. `dd$` calls a DuckDB function, and `na_matches = "never"` is safe because a barcode is never NULL.

``` r
library(duckplyr)

duckplyr_macros <- c(
  '"=="(x, y) AS (x == y)',
  '"!"(x) AS (NOT x)',
  '"n"() AS CAST(COUNT(*) AS int32)',
  '"___divide"(x, y) AS CASE WHEN y = 0 THEN CASE WHEN x = 0 THEN CAST(\'NaN\' AS double)
     WHEN x > 0 THEN CAST(\'+Infinity\' AS double) ELSE CAST(\'-Infinity\' AS double) END
     ELSE CAST(x AS double) / y END'
)
for (macro in duckplyr_macros) {
  dbExecute(con, paste("CREATE MACRO", macro))
}

fragments_df <- read_sql_duckdb("SELECT * FROM fragments", con = con, prudence = "stingy")
cells_df <- read_sql_duckdb("SELECT barcode FROM cells", con = con, prudence = "stingy")

duckplyr_counts <- function() {
  fragments_df |>
    semi_join(cells_df, by = c(name = "barcode"), na_matches = "never") |>
    filter(!dd$duckhts_cgranges_has_overlap("blacklist_idx", chrom, start, end)) |>
    mutate(bin = dd$floor((start + end) / 2 / 50000)) |>
    count(chrom, bin, name) |>
    collect()
}
```

All three must match on every row, then all three return the same data frame to R:

``` r
as_counts <- function(counts) {
  counts <- as.data.table(as.data.frame(counts))
  setnames(counts, c("chrom", "bin", "barcode", "n"))
  counts[, `:=`(bin = as.numeric(bin), n = as.numeric(n))]
  counts
}
sql_rows <- as_counts(dbGetQuery(con, sql_of("one_bin")))
stopifnot(
  fsetequal(as_counts(collect(dbplyr_counts)), sql_rows),
  fsetequal(as_counts(duckplyr_counts()), sql_rows)
)

r_routes <- bench::mark(
  sql = dbGetQuery(con, sql_of("one_bin")),
  dbplyr = collect(dbplyr_counts),
  duckplyr = duckplyr_counts(),
  check = FALSE,
  memory = FALSE,
  filter_gc = FALSE,
  min_time = 0,
  min_iterations = 5,
  max_iterations = 5
)
route_seconds <- as.numeric(r_routes$median)
names(route_seconds) <- as.character(r_routes$expression)
knitr::kable(data.frame(route = names(route_seconds), median_s = round(route_seconds, 2)), row.names = FALSE)
```

| route    | median_s |
|:---------|---------:|
| sql      |     0.53 |
| dbplyr   |     0.57 |
| duckplyr |     0.62 |

On 4M fragments, with 3,506,782 rows returned.

**[Rducks](https://github.com/RGenomicsETL/Rducks)** runs an R function inside the query, for a rule that is easier to write in R. R callbacks need DuckDB on one thread, so keep the counting in SQL and use this for the custom rule only. This one labels fragments by nucleosome count, and the check runs the same function outside DuckDB.

``` r
Rducks::rducks_enable(con, threads = "single")

nucleosomes <- function(width) as.integer(round(as.numeric(width) / 180))

widths_df <- read_sql_duckdb(
  "SELECT \"end\" - start AS width FROM fragments",
  con = con,
  prudence = "stingy"
)
nucleosome_counts <- function() {
  Rducks::rducks_with_duckplyr(
    con,
    widths_df |>
      mutate(nucleosomes_spanned = nucleosomes(width)) |>
      count(nucleosomes_spanned) |>
      collect(),
    returns = list(nucleosomes = Rducks::INTEGER),
    mode = "vectorized"
  )
}

rducks_seconds <- system.time(in_query <- nucleosome_counts())[["elapsed"]]
widths <- dbGetQuery(con, "SELECT \"end\" - start AS width FROM fragments")$width
outside <- as.data.frame(table(nucleosomes(widths)), stringsAsFactors = FALSE)
in_query <- as.data.frame(in_query)
in_query <- in_query[order(in_query$nucleosomes_spanned), ]
stopifnot(
  identical(as.integer(outside$Var1), as.integer(in_query$nucleosomes_spanned)),
  identical(as.numeric(outside$Freq), as.numeric(in_query$n))
)
in_query
#>   nucleosomes_spanned       n
#> 1                   0   36451
#> 2                   1 2744651
#> 3                   2 1205335
#> 4                   3   13563
```

The R function took 2.08 s on 4M fragments, one run, on one DuckDB thread.

<details>
<summary>
Run details
</summary>

    #> source revision: 54656aaaaf82f275a67ef8e6deaf202a82a77923 
    #>  Rduckhts: 2.0.2.9000.0.1.1 
    #>  DuckDB: v1.5.5 
    #>  bedtools: bedtools v2.31.1 
    #>  duckplyr: 1.2.1 
    #>  dbplyr: 2.5.2 
    #>  Rducks: 0.1.1.9000 
    #>  load average at end: 3.00 2.12 1.45 2/2109 3405604 
    #>  CPUs available: 8

Cells stay at 5,000 and only the fragment count grows. Render with `taskset -c 8-15 nice -n -15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_fragment_bin_counts.Rmd")'`.

</details>
