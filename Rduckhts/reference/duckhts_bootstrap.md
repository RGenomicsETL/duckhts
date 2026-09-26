# Bootstrap the duckhts extension sources into the R package

Copies extension source files from the parent duckhts repository into
`inst/duckhts_extension/` so the R package becomes self-contained. Run
this before `R CMD build` to prepare the source tarball.

## Usage

``` r
duckhts_bootstrap(repo_root = NULL)
```

## Arguments

- repo_root:

  Path to the duckhts repository root. Required.

## Value

Invisibly returns the destination directory.
