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
  overwrite = FALSE,
  reference_table = NULL,
  proportions_table = NULL,
  af_clamp = 0.001,
  max_sites = 2e+07,
  max_site_bytes = 2^32
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

  Name of a table or view of site allele frequencies, instead of
  \`af_tag\` or the ancestry inputs.

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

- reference_table:

  Long-format reference with chromosome, position, allele_a, allele_b,
  group_id and frequency columns.

- proportions_table:

  Long-format sample ancestry proportions with sample_id, group_id and
  proportion columns. Its group set must exactly match
  \`reference_table\`.

- af_clamp:

  Clamp ancestry-tuned frequencies to \`\[af_clamp, 1-af_clamp\]\`; zero
  disables clamping. Must be in \`\[0, 0.5)\`.

- max_sites:

  Most sites one sample and chromosome may hold, from 1 to 100,000,000.
  A larger group is an error.

- max_site_bytes:

  Most bytes the site buffers of all samples and chromosomes may hold at
  once (16 bytes per site, 24 for read counts), shared by the decodes
  running in the process. Exceeding it is an error: decode fewer samples
  per call or raise it.

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

Alternatively, \`reference_table\` and \`proportions_table\` select the
ancestry-tuned path. The reference is long-format with columns
\`chromosome\`, \`position\`, \`allele_a\`, \`allele_b\`, \`group_id\`
and \`frequency\` (frequency of \`allele_b\`); proportions have
\`sample_id\`, \`group_id\` and \`proportion\`. Their group sets must
match exactly, and each VCF sample must have a proportions row for every
group. Per-site AF is the sum of each group's frequency weighted by the
sample's proportion, with the proportions divided by their sum (which
must be above 0 and at most 1, so \`sum_to_one = FALSE\` results and
rounded proportions are accepted). REF=\`allele_a\`, ALT=\`allele_b\`
uses that AF; reversed alleles use \`1 - AF\`. Reference alleles must be
on the forward strand of the VCF's assembly, as in a FASTA-anchored
panel; no strand flip is attempted, so palindromic (A/T, C/G) sites are
oriented by REF like any other, and a record whose alleles match neither
order is not used and does not claim a repeated position. Contig names
are matched with \`duckhts_contig_key()\` once per distinct name (one
leading \`chr\` removed, M/MT written as MT, X and Y uppercased), so
\`chr1\` and \`1\` match, as do \`chrX\` and \`X\`; accessions, patches
and numeric sex chromosomes are not mapped. \`af_clamp\` limits
nonzero-clamp frequencies to \`\[af_clamp, 1-af_clamp\]\`; the default
keeps zero population frequencies from being treated as impossible, and
zero disables clamping. Only called sites are frequency-weighted. Supply
exactly one frequency source: \`af_tag\`, \`af_table\`, or both ancestry
relations.

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
#> 2     S1  chr1 1042655 2396807 1354153       113 31.96439
#> 1     S2  chr2  693881 1722963 1029083       114 32.91658
#> 4     S4  chr1 1876100 3319712 1443613       131 44.83104
#> 3     S4  chr2   23420 2400223 2376804       239 46.64399
DBI::dbDisconnect(con, shutdown = TRUE)
```
