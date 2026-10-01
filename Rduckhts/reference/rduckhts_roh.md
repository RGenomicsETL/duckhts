# Find Runs of Homozygosity

Detect runs of homozygosity in a VCF or BCF with the two-state hidden
Markov model of \`bcftools roh\`: an autozygous state and a
Hardy-Weinberg state, decoded with the Viterbi algorithm. Segments,
marker counts and start and end positions equal \`bcftools roh\`;
quality is the mean forward-backward phred score, equal to the one
decimal \`bcftools\` prints. The model is a port of \`vcfroh.c\` and
\`HMM.c\` (MIT, Genome Research Ltd), run by the native
\`duckhts_roh_segments()\` kernel over one list per sample and
chromosome.

## Usage

``` r
rduckhts_roh(
  con,
  path,
  af_tag = NULL,
  af_table = NULL,
  genetic_map = NULL,
  hw_to_az = 6.7e-08,
  az_to_hw = 5e-09,
  gt_error = NULL,
  rec_rate = NULL,
  samples = NULL,
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- path:

  One VCF/BCF path or URI.

- af_tag:

  INFO tag holding the ALT allele frequency, for example \`"AF"\`.

- af_table:

  Name of a table or view of allele frequencies, instead of \`af_tag\`.
  Exactly one of the two is required.

- genetic_map:

  Optional name of a genetic-map table or view.

- hw_to_az, az_to_hw:

  Per-base-pair transition probabilities Hardy-Weinberg to autozygous
  (\`-a\`) and back (\`-H\`), each in \`\[0, 1\]\`.

- gt_error:

  Optional phred error for GT-only emissions (\`-G\`), at least 0.

- rec_rate:

  Optional constant recombination rate per base pair (\`-M\`).

- samples:

  Optional HTSlib sample selector: comma-separated inclusion, leading
  \`^\` exclusion, \`"-"\` for all, or \`""\` for none.

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame (or invisible \`TRUE\` when \`table_name\` is given) with
one row per run: \`sample\`, \`chrom\`, \`start\` and \`end\`
(one-based, inclusive, the first and last marker), \`length\` (\`end -
start + 1\`), \`n_markers\` and \`quality\`. Rows are unordered.

## Details

Allele frequencies come from an INFO tag (\`af_tag\`, like \`–AF-tag\`,
which must be declared \`Type=Float,Number=A\`) or from a relation
(\`af_table\`, like \`–AF-file\`) with columns \`chrom\`, \`pos\`,
\`ref\`, \`alt\` and \`af\`, matched to each record by chromosome,
position, REF and its ALT alleles joined by commas. Sites without a
usable frequency (absent, missing or exactly 0) are skipped, as are
records with more than one ALT or none.

Genotype evidence is FORMAT/PL unless \`gt_error\` is given, which uses
the diploid GT calls with that phred error (\`-G\`) and needs no PL.
With no recombination map, transitions are the per-base-pair
probabilities \`hw_to_az\` and \`az_to_hw\` compounded over the physical
distance between sites, exactly as \`bcftools roh\` does without \`-m\`
or \`-M\`. \`rec_rate\` is a constant rate per base pair (\`-M\`).
\`genetic_map\` names a relation with columns \`chrom\`, \`pos\` and
\`cm\` (cumulative centimorgans, the IMPUTE2 format), interpolated as
\`-m\` does; chromosomes it does not cover are skipped.

## Examples

``` r
con <- rduckhts_connect()
path <- system.file("extdata", "roh_fixture.vcf.gz", package = "Rduckhts")
roh <- rduckhts_roh(con, path, af_tag = "AF")
roh[order(roh$sample, roh$chrom, roh$start), ]
#>   sample chrom   start     end  length n_markers  quality
#> 1     S1  chr1 1042655 2396807 1354153       113 31.96439
#> 4     S2  chr2  693881 1722963 1029083       114 32.91658
#> 3     S4  chr1 1876100 3319712 1443613       131 44.83104
#> 2     S4  chr2   23420 2400223 2376804       239 46.64399
DBI::dbDisconnect(con, shutdown = TRUE)
```
