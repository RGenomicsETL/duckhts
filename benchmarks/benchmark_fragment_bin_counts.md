Counting BED fragments by tile and cell barcode in SQL
================

<!-- benchmark_fragment_bin_counts.md is generated from benchmark_fragment_bin_counts.Rmd. -->

@ttriche happy that the package is useful to you!

On the bin count question, for these things i usually start by composing some sql that does the equivalent work, and if it is not satisfactory add a native kernel. Some tests of this approach are in [duckseq](https://github.com/RGenomicsETL/duckseq/tree/main).

Indeed i thought of adding a fragment count mode to `bam_bin_counts`. If this is useful i can add it with the blacklisting, cgranges should make that relatively fast.

For now this should work. i am guessing your setup from the file names in [\#356](https://github.com/RGenomicsETL/duckhts/issues/356): `chrom, start, end, name` fragments with the barcode in `name`, 50 kb tiles, and a sparse tile by cell matrix at the end. The data below is simulated.

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

`end` needs quotes in SQL. Coordinates stay 0-based half-open.

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

**One bin per fragment**, by midpoint:

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
| chr1  | 3919 | AGACTCGATTTGCATA |   1 |
| chr1  | 3941 | TAGGGTGATGGGGTGC |   1 |
| chr1  | 3961 | GCCTATGTCGACCTAA |   1 |
| chr1  | 3973 | GGGAGTACAAAGTAGT |   1 |

Displaying records 1 - 4

</div>

**Every overlapped bin**, so a fragment across a tile edge counts in both:

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
| chr1  | 3919 | AGACTCGATTTGCATA |   1 |
| chr1  | 3941 | TAGGGTGATGGGGTGC |   1 |
| chr1  | 3961 | GCCTATGTCGACCTAA |   1 |
| chr1  | 3973 | GGGAGTACAAAGTAGT |   1 |

Displaying records 1 - 4

</div>

**Any tile or peak BED**, with integer ids for the matrix:

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
|   24475 |     219 |   1 |
|   15740 |     269 |   1 |
|    9577 |     281 |   1 |

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

`bench::mark()` runs the queries above from `read_bed` and from a Parquet copy.

<details>
<summary>
Extra queries in the table
</summary>

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

| chrom | bin | barcode          |   n |
|:------|----:|:-----------------|----:|
| chr2  |  94 | GGGTGAGAATCGGAAT |   1 |
| chr2  | 129 | GCAGGGTTGCCCGTAC |   1 |
| chr2  | 309 | GTATTACTACCTCACC |   1 |
| chr2  | 344 | TATGTCGGCACGGTGA |   1 |

Displaying records 1 - 4

</div>

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
| chr18 |  843 |     154 |   1 |
| chrY  |  226 |     246 |   1 |
| chr18 | 1508 |     833 |   1 |
| chr18 |  836 |     922 |   1 |

Displaying records 1 - 4

</div>

</details>

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
for (size in names(sizes)) {
  parquet_file <- sprintf("fragments_%s.parquet", size)
  convert_seconds[size] <- system.time(dbExecute(con, sprintf(
    "COPY (SELECT chrom, start, \"end\", name FROM read_bed('%s'))
     TO '%s' (FORMAT parquet, COMPRESSION zstd)",
    fragment_files[[size]], parquet_file
  )))[["elapsed"]]

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
| reader_only           |        0.67 |       0.01 |        1.32 |       0.02 |         2.60 |        0.03 |
| one_bin               |        1.82 |       0.32 |        3.62 |       0.62 |         7.20 |        1.22 |
| one_bin_cell_id       |        1.61 |       0.28 |        3.25 |       0.53 |         6.44 |        1.05 |
| one_bin_sql_blacklist |        2.33 |       0.57 |        4.75 |       1.06 |         9.38 |        2.01 |
| all_bins              |        1.93 |       0.40 |        3.83 |       0.78 |         7.55 |        1.54 |
| tile_cell             |        2.76 |       1.09 |        5.50 |       2.16 |        10.80 |        4.34 |

Median seconds of 5 runs, 4 threads. Time grows linearly with fragments.

- `read_bed` runs the whole query on one thread. Parquet is 5.9 times faster, so if you count a file more than once, convert it first (5.6 s for 16M).
- Barcode strings cost little: integer ids save 14%.
- cgranges blacklisting is 1.6 times faster than a SQL anti join.

# From R

The same `one_bin` count on the same connection, with no quoting of `end`.

dbplyr:

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

duckplyr, with `prudence = "stingy"` so nothing falls back to R. It needs four of its macros on this connection:

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
| sql      |     0.52 |
| dbplyr   |     0.55 |
| duckplyr |     0.63 |

Median seconds on 4M fragments. All three return the same rows.

[Rducks](https://github.com/RGenomicsETL/Rducks) runs an R function inside the query, on one DuckDB thread, so keep it for rules that are awkward in SQL:

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

That took 2.1 s on 4M fragments.

<details>
<summary>
Run details
</summary>

    #> source revision: b2a0aa055a807201ed4a35770cae75744e45f381 
    #>  Rduckhts: 2.0.2.9000.0.1.1 
    #>  DuckDB: v1.5.5 
    #>  bedtools: bedtools v2.31.1 
    #>  duckplyr: 1.2.1 
    #>  dbplyr: 2.5.2 
    #>  Rducks: 0.1.1.9000 
    #>  load average at end: 2.69 2.16 1.68 1/2124 3417482 
    #>  CPUs available: 8

Peak RSS of a fresh R process running `tile_cell` on 16M fragments: 1946 MiB from `read_bed`, 618 MiB from Parquet. Cells stay at 5,000 and only the fragment count grows. Render with `taskset -c 8-15 nice -n -15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_fragment_bin_counts.Rmd")'`.

</details>
