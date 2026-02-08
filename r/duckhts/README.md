# duckhts (R)

This package is a lightweight harness for building and testing the DuckHTS DuckDB extension from R.

## Bootstrap

```r
source("bootstrap.R")
```

To include conformance datasets:

```r
system("Rscript bootstrap.R --conformance")
```

## Build

```r
duckhts_build()
```

## Load

```r
con <- duckhts_load()
```
