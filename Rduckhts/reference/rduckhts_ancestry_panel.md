# Create a deterministic ancestry site panel for BAM/CRAM counts

Sites are biallelic unambiguous SNVs present in every reference group
and every PC loading. Positions are positive one-based whole numbers and
all group alleles must be present. With \`candidate_table\`, only
matching candidate loci and alleles are retained. Otherwise one site is
selected per chromosome/spaced genomic window. Results are sorted by
chromosome and position and assigned dense zero-based panel ordinals.
Keep the returned panel SHA-256 and the reference release with the
materialized panel as its versioned identity.

## Usage

``` r
rduckhts_ancestry_panel(
  con,
  reference_table,
  loadings_table,
  table_name,
  assembly,
  candidate_table = NULL,
  spacing_bp = 5000,
  max_sites = 17000,
  overwrite = FALSE
)
```

## Arguments

- con:

  A connection with DuckHTS loaded.

- reference_table, loadings_table:

  Reference products on the connection.

- table_name:

  Destination committed panel relation, visible to the BAM/CRAM panel
  preparation path.

- assembly:

  Reference genome assembly label.

- candidate_table:

  Optional site panel, e.g. a Somalier panel with region, position,
  allele_a and allele_b.

- spacing_bp:

  Genomic window width for deterministic spacing.

- max_sites:

  Maximum number of selected sites. When more sites are eligible, each
  contig keeps at least one site (the largest contigs, when
  \`max_sites\` is below the contig count), the rest are shared by
  largest remainders in proportion to each contig's eligible sites minus
  that one (its remaining capacity, so no contig is offered more than it
  has), and each contig's sites are spread evenly along it.

- overwrite:

  Replace an existing destination.

## Value

A one-row data frame with panel SHA-256 and selected site count.
