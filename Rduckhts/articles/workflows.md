# Worked examples

These examples run on the small test files bundled with Rduckhts. Each
section is self-contained after the setup below; outputs are real, but
the tiny fixtures make the numbers illustrations, not biological
results.

## Setup

``` r
library(DBI)
library(duckdb)
library(Rduckhts)

setup_hts_env()

fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
con <- rduckhts_connect()

rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE)
rduckhts_fastq(con, "reads", fastq_r1, mate_path = fastq_r2, overwrite = TRUE)

dbGetQuery(con, "SELECT COUNT(*) AS n FROM sequences")
```

    ##   n
    ## 1 7

``` r
dbGetQuery(con, "SELECT COUNT(*) AS n FROM reads")
```

    ##    n
    ## 1 10

## Variants

### Region queries

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

rduckhts_bcf(
  con, "variants_idx", bcf_path,
  region = "1:3000150-3000151",
  index_path = bcf_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT CHROM, POS, REF, ALT FROM variants_idx")
```

    ##   CHROM     POS REF ALT
    ## 1     1 3000150   C   T
    ## 2     1 3000151   C   T

``` r
rduckhts_bam(
  con, "bam_idx_reads", bam_path,
  region = "CHROMOSOME_I:1-1000",
  index_path = bam_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM bam_idx_reads")
```

    ##                           QNAME FLAG POS MAPQ
    ## 1 HS18_09653:4:1315:19857:61712  145 914   23
    ## 2 HS18_09653:4:1308:11522:27107  161 934    0

### Parquet conversion and index spans

Region queries can use implicit sidecar indexes or an explicit
`index_path` for custom index names/locations.

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
variants <- dbGetQuery(con, "SELECT * FROM variants LIMIT 5")
variants
```

    ##   CHROM     POS    ID  REF  ALT  QUAL FILTER INFO_TEST   INFO_DP4 INFO_AC
    ## 1     1 3000150  <NA>    C    T  59.2   PASS        NA       NULL       2
    ## 2     1 3000151  <NA>    C    T  59.2   PASS        NA       NULL       2
    ## 3     1 3062915  id3D GTTT    G  12.9    q10        NA 1, 2, 3, 4       2
    ## 4     1 3062915 idSNP    G T, C  12.6   test         5 1, 2, 3, 4    1, 1
    ## 5     1 3106154  <NA> CAAA    C 342.0   PASS        NA       NULL       2
    ##   INFO_AN INFO_INDEL INFO_STR FORMAT_TT_A FORMAT_GT_A FORMAT_GQ_A FORMAT_DP_A
    ## 1       4      FALSE     <NA>        NULL         0/1         245          NA
    ## 2       4      FALSE     <NA>        NULL         0/1         245          32
    ## 3       4       TRUE     test        NULL         0/1         409          35
    ## 4       3      FALSE     <NA>        0, 1         0/1         409          35
    ## 5       4      FALSE     <NA>        NULL         0/1         245          32
    ##                  FORMAT_GL_A FORMAT_TT_B FORMAT_GT_B FORMAT_GQ_B FORMAT_DP_B
    ## 1                       NULL        NULL         0/1         245          NA
    ## 2                       NULL        NULL         0/1         245          32
    ## 3               -20, -5, -20        NULL         0/1         409          35
    ## 4 -20, -5, -20, -20, -5, -20        0, 1           2         409          35
    ## 5                       NULL        NULL         0/1         245          32
    ##    FORMAT_GL_B
    ## 1         NULL
    ## 2         NULL
    ## 3 -20, -5, -20
    ## 4 -20, -5, -20
    ## 5         NULL

``` r
rduckhts_bcf(
  con, "variants_idx", bcf_path,
  region = "1:3000150-3000151",
  index_path = bcf_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT count(*) AS n FROM variants_idx")
```

    ##   n
    ## 1 2

``` r
# Convert to a round-trip-aware Parquet copy with VCF header metadata.
# Extra named metadata is merged into the Parquet key-value metadata.
parquet_path <- tempfile(fileext = ".parquet")
rduckhts_bcf_convert_parquet(
  con,
  bcf_path,
  parquet_path,
  columns = c("CHROM", "POS", "REF", "ALT"),
  metadata = list(project = "demo-cohort", batch = "1"),
  overwrite = TRUE
)

metadata_preview <- dbGetQuery(con, sprintf(
  "SELECT key::VARCHAR AS key, left(value::VARCHAR, 40) AS value_prefix
   FROM parquet_kv_metadata('%s')
   WHERE key::VARCHAR IN ('duckhts_write_format_version', 'duckhts_reader', 'project', 'batch', 'vcf_header')
   ORDER BY key::VARCHAR",
  parquet_path
))
metadata_preview
```

    ##                            key                              value_prefix
    ## 1                        batch                                         1
    ## 2               duckhts_reader                                  read_bcf
    ## 3 duckhts_write_format_version                                         1
    ## 4                      project                               demo-cohort
    ## 5                   vcf_header ##fileformat=VCFv4.1\\x5Cn##FILTER=<ID=PA

``` r
# The same converter helpers support partitioned output for DuckLake-style
# registration of premade Parquet files.
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
gff_parquet_dir <- tempfile("duckhts_gff_parquet_")
rduckhts_gff_convert_parquet(
  con,
  gff_path,
  gff_parquet_dir,
  columns = c("seqname", "source", "feature", "start", "end"),
  partition_by = "feature",
  overwrite = TRUE
)
length(list.files(gff_parquet_dir, pattern = "\\.parquet$", recursive = TRUE))
```

    ## [1] 5

``` r
unlink(c(parquet_path, gff_parquet_dir), recursive = TRUE)

# Span-oriented index view from the same file
index_spans_preview <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
head(index_spans_preview[, c("seqname", "tid", "index_type", "chunk_beg_vo", "chunk_end_vo")], 5)
```

    ##   seqname tid index_type chunk_beg_vo chunk_end_vo
    ## 1       1   0        CSI         1586         1713
    ## 2       1   0        CSI         1713         1973
    ## 3       1   0        CSI         1973         2109
    ## 4       1   0        CSI         2109         2242
    ## 5       1   0        CSI         2242         2372

### Variant normalization

[`rduckhts_bcftools_norm()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_bcftools_norm.md)
wraps the bundled `duckhts_bcftools_norm(...)` table macro and keeps the
original input columns alongside normalized position/reference/ALT
outputs.

``` r
norm_fa <- system.file("extdata", "liftover_repeat_src.fa", package = "Rduckhts")
invisible(DBI::dbExecute(
  con,
  paste(
    "CREATE OR REPLACE TEMP TABLE readme_norm AS SELECT * FROM (VALUES",
    "('chrS', 2, 'T', 'TT,TTT'),",
    "('chrS', 2, 'T', '*,TT')",
    ") AS t(chrom, pos, ref, alt)"
  )
))

norm_out <- rduckhts_bcftools_norm(
  con,
  "readme_norm",
  norm_fa,
  split_multiallelic = TRUE
)
norm_out[order(norm_out$alt, norm_out$alt_index),
         c("chrom", "pos", "ref", "alt", "alt_index", "pos_normed", "ref_normed", "alt_normed", "norm_status")]
```

    ##   chrom pos ref    alt alt_index pos_normed ref_normed alt_normed
    ## 4  chrS   2   T   *,TT         1          2          T          *
    ## 2  chrS   2   T   *,TT         2          1          G         GT
    ## 3  chrS   2   T TT,TTT         1          1          G         GT
    ## 1  chrS   2   T TT,TTT         2          1          G        GTT
    ##        norm_status
    ## 4 SpanningDeletion
    ## 2       Normalized
    ## 3       Normalized
    ## 1       Normalized

### VariantKey and RegionKey

DuckHTS now bundles the official VariantKey / RegionKey C API. The SQL
helper `variantkey(...)` accepts 1-based VCF positions to match bcftools
`%VKX`, while `regionkey(...)` keeps 0-based half-open interval
semantics for span work. Large, ambiguous, and symbolic alleles still
encode through the official hashed nonreversible VariantKey mode, but
those keys do not encode `END`, `SVLEN`, mate breakend coordinates, or
other SV metadata; use RegionKey explicitly for interval/SV span
indexing. See Nicola Asuni (2018) <https://doi.org/10.1101/473744>.

``` r
dbGetQuery(
  con,
  paste(
    "SELECT variantkey_hex(variantkey('1', 324684, 'C', 'G')) AS vkx,",
    "reverse_variantkey(parse_variantkey_hex('08027a2588b00000')) AS reversed"
  )
)
```

    ##                vkx reversed.chrom reversed.chrom_code reversed.pos
    ## 1 08027a2588b00000              1                   1       324684
    ##   reversed.pos0 reversed.ref reversed.alt reversed.refalt_code
    ## 1        324683            C            G            145752064
    ##   reversed.reversible
    ## 1                TRUE

``` r
dbGetQuery(
  con,
  paste(
    "SELECT regionkey_hex(regionkey('X', 1007, 1807, 1)) AS rkx,",
    "are_overlapping_regionkeys(regionkey('X', 1007, 1807, 1), parse_regionkey_hex('b80001f78000387a')) AS overlaps"
  )
)
```

    ##                rkx overlaps
    ## 1 b80001f78000387a     TRUE

### Liftover

``` r
lift_src <- tempfile("duckhts_liftover_src_", fileext = ".fa")
lift_dst <- tempfile("duckhts_liftover_dst_", fileext = ".fa")
lift_chain <- tempfile("duckhts_liftover_", fileext = ".chain")

writeLines(c(
  ">chrF",
  "ACGTACGTAA",
  ">chrR",
  "AACCGGTTAA"
), lift_src)
writeLines(c(
  ">chrLiftF",
  "ACGTACGTAA",
  ">chrLiftR",
  "TTAACCGGTT"
), lift_dst)
writeLines(c(
  "chain 100 chrF 10 + 0 10 chrLiftF 10 + 0 10 1",
  "10",
  "",
  "chain 100 chrR 10 + 0 10 chrLiftR 10 - 0 10 2",
  "10"
), lift_chain)

lift_src_index <- rduckhts_fasta_index(
  con, lift_src, index_path = paste0(lift_src, ".fai")
)
lift_src_index$index_path <- "<tempfile>"
lift_src_index
```

    ##   success index_path
    ## 1    TRUE <tempfile>

``` r
lift_dst_index <- rduckhts_fasta_index(
  con, lift_dst, index_path = paste0(lift_dst, ".fai")
)
lift_dst_index$index_path <- "<tempfile>"
lift_dst_index
```

    ##   success index_path
    ## 1    TRUE <tempfile>

``` r
lifted <- rduckhts_liftover(
  con,
  query = paste(
    "SELECT * FROM (VALUES",
    "('chrF', 2, 'C', 'T'),",
    "('chrR', 2, 'A', 'G'),",
    "('chrF', 11, 'A', 'T')",
    ") AS t(chrom, pos, ref, alt)"
  ),
  chain_path = lift_chain,
  dst_fasta_ref = lift_dst,
  ref_col = "ref",
  alt_col = "alt",
  src_fasta_ref = lift_src
)

lifted[, c(
  "src_chrom", "src_pos", "dest_chrom", "dest_pos",
  "dest_ref", "dest_alt", "mapped", "reverse_complemented",
  "reject_reason", "note"
)]
```

    ##   src_chrom src_pos dest_chrom dest_pos dest_ref dest_alt mapped
    ## 1      chrF       2   chrLiftF        2        C        T   TRUE
    ## 2      chrR       2   chrLiftR        9        T        C   TRUE
    ## 3      chrF      11       <NA>       NA     <NA>     <NA>  FALSE
    ##   reverse_complemented     reject_reason note
    ## 1                FALSE              <NA> <NA>
    ## 2                 TRUE              <NA> <NA>
    ## 3                FALSE SourceRefMismatch <NA>

``` r
unlink(c(lift_src, paste0(lift_src, ".fai"), lift_dst, paste0(lift_dst, ".fai"), lift_chain))
```

### Summary-statistics munging

``` r
munge_fasta <- tempfile("duckhts_munge_", fileext = ".fa")
writeLines(c(
  ">chrF",
  "ACGTACGTAA"
), munge_fasta)
transform(rduckhts_fasta_index(con, munge_fasta, index_path = paste0(munge_fasta, ".fai")),
          index_path = "<tempfile>")
```

    ##   success index_path
    ## 1    TRUE <tempfile>

``` r
munge_out <- rduckhts_munge(
  con,
  query = paste(
    "SELECT * FROM (VALUES",
    "('rs1', 2, 'chrF', 'A', 'C', 0.01, 1.10, 0.20, 0.98, 0.10, 0.01, 1000),",
    "('rs2', 2, 'chrF', 'C', 'A', 0.02, 0.90, -0.20, 0.98, 0.90, 0.01, 1000)",
    ") AS t(SNP, BP, CHR, A1, A2, P, OR_VALUE, BETA, INFO, FRQ, SE, N)"
  ),
  fasta_ref = munge_fasta,
  column_map = c(
    SNP = "SNP", BP = "BP", CHR = "CHR", A1 = "A1", A2 = "A2",
    P = "P", OR = "OR_VALUE", BETA = "BETA", INFO = "INFO", FRQ = "FRQ", SE = "SE", N = "N"
  )
)

munge_out[, c("chrom", "pos", "id", "ref", "alt", "alleles_swapped", "filter", "af", "es", "ns")]
```

    ##   chrom pos  id ref alt alleles_swapped filter  af  es   ns
    ## 1  chrF   2 rs2   C   A            TRUE   <NA> 0.1 0.2 1000
    ## 2  chrF   2 rs1   C   A           FALSE   <NA> 0.1 0.2 1000

``` r
unlink(c(munge_fasta, paste0(munge_fasta, ".fai")))
```

### Polygenic scores

[`rduckhts_score()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_score.md)
computes per-sample polygenic risk scores (PRS) from a genotype VCF/BCF
and one or more GWAS summary statistics files, wrapping
`bcftools_score`.

``` r
vcf_path    <- system.file("extdata", "score_input.vcf",        package = "Rduckhts")
dosage_path <- system.file("extdata", "score_dosage.vcf",       package = "Rduckhts")
sumf_path   <- system.file("extdata", "score_summary.tsv",      package = "Rduckhts")
gwas_path   <- system.file("extdata", "score_gwas_summary.vcf", package = "Rduckhts")

# Hard-call (GT) PRS with PLINK-format summary statistics
# S1: 0×0.5 + 1×(−0.2) + 2×1.0 = 1.8
# S2: 1×0.5 + 2×(−0.2) + 0×1.0 = 0.1
gt_prs <- rduckhts_score(con, vcf_path, sumf_path, use = "GT", columns = "PLINK")
gt_prs[, c("SAMPLE", "score_summary")]
```

    ##   SAMPLE score_summary
    ## 1     S1    1.79999995
    ## 2     S2    0.09999999

``` r
# Multiple TSV/SSF summary files are scored in one genotype scan
sumf_na_path <- system.file("extdata", "score_summary_na.tsv", package = "Rduckhts")
multi_prs <- rduckhts_score(con, vcf_path, c(sumf_path, sumf_na_path),
                            use = "GT", columns = "PLINK")
multi_prs[, c("SAMPLE", "score_summary", "score_summary_na")]
```

    ##   SAMPLE score_summary score_summary_na
    ## 1     S1    1.79999995              2.0
    ## 2     S2    0.09999999              0.5

``` r
# Optional audit log records loaded/matched/allele-mismatch marker counts
sumf_mismatch_path <- system.file("extdata", "score_summary_mismatch.tsv", package = "Rduckhts")
score_log <- tempfile("duckhts_score_", fileext = ".log")
invisible(rduckhts_score(con, vcf_path, c(sumf_path, sumf_mismatch_path),
                         use = "GT", columns = "PLINK", log_path = score_log))
read.delim(score_log, comment.char = "#")[, c("summary_name", "loaded_markers",
                                                "matched_markers", "allele_mismatch_markers")]
```

    ##             summary_name loaded_markers matched_markers allele_mismatch_markers
    ## 1          score_summary              3               3                       0
    ## 2 score_summary_mismatch              3               0                       3

``` r
# Dosage-based PRS (DS field) for imputed genotypes
# S1: 0.1×0.5 + 0.8×(−0.2) + 1.8×1.0 = 1.69
# S2: 1.0×0.5 + 1.9×(−0.2) + 0.2×1.0 = 0.32
ds_prs <- rduckhts_score(con, dosage_path, sumf_path, use = "DS", columns = "PLINK")
ds_prs[, c("SAMPLE", "score_summary")]
```

    ##   SAMPLE score_summary
    ## 1     S1          1.69
    ## 2     S2          0.32

``` r
# GWAS-VCF multi-PRS: each FORMAT/ES sample column becomes a separate PRS track
gwas_prs <- rduckhts_score(con, vcf_path, gwas_path, use = "GT")
gwas_prs[, c("SAMPLE", "PRS_A", "PRS_B")]
```

    ##   SAMPLE      PRS_A PRS_B
    ## 1     S1 1.79999995   1.0
    ## 2     S2 0.09999999   0.3

### Sample relatedness and contamination

The Somalier-derived workflow consumes ordinary typed relations. The
panel defines an ordered A/B allele orientation, and count extraction
preserves one row per sample and panel site. At the second site below,
the bundled VCF is multiallelic: the allele outside the panel’s A/B pair
is counted as `other`.

This two-site fixture keeps the example runnable and shows the complete
API; its estimates are not meaningful population or sample QC results.
Real analyses need a validated genome-wide panel and its matching
population frequencies.

``` r
sites_vcf <- system.file(
  "extdata", "mapping_number_families.vcf", package = "Rduckhts"
)

invisible(dbExecute(con, paste(
  "CREATE TEMP TABLE identity_panel AS SELECT * FROM (VALUES",
  "('GRCh38', 0::UBIGINT, 'chr1', 100::UBIGINT, 'A', 'C', 0.25::DOUBLE),",
  "('GRCh38', 1::UBIGINT, 'chr1', 200::UBIGINT, 'A', 'G', 0.80::DOUBLE)",
  ") p(assembly, site_index, region, position, allele_a, allele_b,",
  "population_b_af)"
)))

rduckhts_somalier_vcf_counts(
  con, sites_vcf, panel_table = "identity_panel", samples = "S1,S2",
  table_name = "identity_counts"
)
dbGetQuery(con, paste(
  "SELECT sample_id, site_index, a, b, other, status",
  "FROM identity_counts ORDER BY sample_id, site_index"
))
```

    ##   sample_id site_index a  b other   status
    ## 1        S1          0 9  3     0 measured
    ## 2        S1          1 5  9     4 measured
    ## 3        S2          0 0 20     0 measured
    ## 4        S2          1 0 10     5 measured

Build one reusable packed sketch per sample, then compare the sketches.
A production result should retain the panel, counts, and sketches so the
reported denominators and identities can be verified later.

``` r
rduckhts_somalier_sketches(
  con, evidence_table = "identity_counts", panel_table = "identity_panel",
  table_name = "identity_sketches", min_depth = 7,
  min_het_balance = 0.20, hom_balance_cutoff = 0.05, max_sites = 2
)

relatedness <- rduckhts_somalier_relatedness(
  con, sketches_table = "identity_sketches", max_sites = 2
)
relatedness[, c(
  "sample_a", "sample_b", "jointly_called", "ibs0", "ibs2", "relatedness"
)]
```

    ##   sample_a sample_b jointly_called ibs0 ibs2 relatedness
    ## 1       S1       S2              1    0    0           0

CHARR uses population-B frequencies from the panel. Matched
contamination is directional, so the requested relation names the
receiver and anchor explicitly.

``` r
charr <- rduckhts_somalier_charr(
  con, evidence_table = "identity_counts", panel_table = "identity_panel",
  frequency_table = "identity_panel", min_depth = 7,
  hom_minor_rate = 0.10, hom_tail_alpha = 0.001, max_sites = 2
)
print(charr[order(charr$sample_id),
            c("sample_id", "status", "usable_sites", "estimate")], row.names = FALSE)
```

    ##  sample_id status usable_sites estimate
    ##         S1     ok            1        1
    ##         S2     ok            1        0

``` r
invisible(dbExecute(con, paste(
  "CREATE TEMP TABLE contamination_pairs AS",
  "SELECT 'S1'::VARCHAR receiver_id, 'S2'::VARCHAR anchor_id"
)))
matched <- rduckhts_somalier_matched_contamination(
  con, evidence_table = "identity_counts", panel_table = "identity_panel",
  frequency_table = "identity_panel", pairs_table = "contamination_pairs",
  min_depth = 7, hom_minor_rate = 0.10, hom_tail_alpha = 0.001,
  max_sites = 2
)
matched[, c(
  "receiver_id", "anchor_id", "status", "usable_sites", "alpha"
)]
```

    ##   receiver_id anchor_id status usable_sites     alpha
    ## 1          S1        S2     ok            1 0.7541658

## Alignments and coverage

### BAM and CRAM

When built with htslib codec, `CRAM` can be opened in addition to `BAM`
files. `index_path` can also be passed for region scans with
non-standard index names.

``` r
cram_path <- system.file("extdata", "range.cram", package = "Rduckhts")
ref_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

rduckhts_bam(con, "cram_reads", cram_path, reference = ref_path, overwrite = TRUE)
cram_reads <- dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM cram_reads LIMIT 5")
cram_reads
```

    ##                           QNAME FLAG  POS MAPQ
    ## 1 HS18_09653:4:1315:19857:61712  145  914   23
    ## 2 HS18_09653:4:1308:11522:27107  161  934    0
    ## 3 HS18_09653:4:2314:14991:85680   83 1020   10
    ## 4 HS18_09653:4:2108:14085:93656  147 1122   60
    ## 5  HS18_09653:4:1303:4347:38100   83 1137   37

``` r
rduckhts_bam(
  con, "bam_idx_reads", bam_path,
  region = "CHROMOSOME_I:1-1000",
  index_path = bam_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT count(*) AS n FROM bam_idx_reads")
```

    ##   n
    ## 1 2

### SAM tags

Standard SAMtags can be exposed as typed columns, and any remaining tags
are available via `AUXILIARY_TAGS`:

``` r
aux_path <- system.file("extdata", "aux_tags.sam.gz", package = "Rduckhts")
rduckhts_bam(con, "aux_reads", aux_path, standard_tags = TRUE, auxiliary_tags = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT RG, NM, map_extract(AUXILIARY_TAGS, 'XZ') AS XZ FROM aux_reads LIMIT 1")
```

    ##   RG NM  XZ
    ## 1 x1  2 foo

### Fixed-bin counts

[`rduckhts_bam_bin_counts()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_bam_bin_counts.md)
exposes the native fixed-width read-start counting kernel. It returns
one row per bin across the selected contig span, including zero-count
bins, which makes it suitable as a dense CNV binning primitive. It can
also add one-pass GC and MAPQ summaries on the same scan.

``` r
mixed_cram <- system.file("extdata", "fixture_mixed.cram", package = "Rduckhts")
fixture_ref <- system.file("extdata", "fixture_ref.fa", package = "Rduckhts")

bin_counts <- rduckhts_bam_bin_counts(
  con,
  mixed_cram,
  5000,
  reference = fixture_ref,
  rmdup = "streaming",
  stats = "gc,mq"
)

bin_counts[, c(
  "bin_id", "count_total", "count_fwd", "count_rev",
  "count_pre", "gc_perc_pre", "gc_perc_post", "mean_mapq_post"
)]
```

    ##    bin_id count_total count_fwd count_rev count_pre gc_perc_pre gc_perc_post
    ## 1       0           2         1         1         4         0.5            0
    ## 2       1           2         1         1         2         0.0            0
    ## 3       2           1         1         0         2         1.0            1
    ## 4       3           0         0         0         0          NA           NA
    ## 5       4           0         0         0         0          NA           NA
    ## 6       5           0         0         0         0          NA           NA
    ## 7       6           0         0         0         0          NA           NA
    ## 8       7           0         0         0         0          NA           NA
    ## 9       8           0         0         0         0          NA           NA
    ## 10      9           0         0         0         0          NA           NA
    ##    mean_mapq_post
    ## 1              60
    ## 2              60
    ## 3              60
    ## 4              NA
    ## 5              NA
    ## 6              NA
    ## 7              NA
    ## 8              NA
    ## 9              NA
    ## 10             NA

### Mosdepth-compatible coverage

[`rduckhts_mosdepth()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_mosdepth.md)
writes mosdepth-style outputs to disk and returns the paths it created.
This example writes windowed fragment coverage from the bundled BAM
fixture and previews the generated regions BED.gz.

``` r
mos_prefix <- tempfile("duckhts_readme_mosdepth_")
mos_out <- rduckhts_mosdepth(
  con,
  prefix = mos_prefix,
  path = bam_path,
  chrom = "CHROMOSOME_II",
  by = "1000",
  no_per_base = TRUE,
  fragment_mode = TRUE,
  use_median = TRUE,
  overwrite = TRUE
)

transform(mos_out[, c("summary_path", "regions_path")],
          summary_path = "<tempfile>", regions_path = "<tempfile>")
```

    ##   summary_path regions_path
    ## 1   <tempfile>   <tempfile>

``` r
utils::read.delim(
  gzfile(mos_out$regions_path[[1]]),
  header = FALSE,
  sep = "\t",
  nrows = 3,
  col.names = c("chrom", "start", "end", "depth")
)
```

    ##           chrom start  end depth
    ## 1 CHROMOSOME_II     0 1000     0
    ## 2 CHROMOSOME_II  1000 2000     5
    ## 3 CHROMOSOME_II  2000 3000     3

``` r
unlink(
  c(
    paste0(mos_prefix, ".mosdepth.summary.txt"),
    paste0(mos_prefix, ".mosdepth.global.dist.txt"),
    paste0(mos_prefix, ".mosdepth.region.dist.txt"),
    paste0(mos_prefix, ".regions.bed.gz"),
    paste0(mos_prefix, ".regions.bed.gz.csi")
  ),
  force = TRUE
)
```

## Sequences

### FASTQ files

Three modes for fastq files, single, paired and interleaved

``` r
r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
interleaved <- system.file("extdata", "interleaved.fq", package = "Rduckhts")
rduckhts_fastq(con, "paired_reads", r1, mate_path = r2, overwrite = TRUE)
rduckhts_fastq(con, "interleaved_reads", interleaved, interleaved = TRUE, overwrite = TRUE)
pairs <- dbGetQuery(con, "SELECT * FROM paired_reads WHERE MATE = 1 LIMIT 5")
pairs
```

    ##                              NAME DESCRIPTION
    ## 1 HS25_09827:2:1201:1505:59795#49        <NA>
    ## 2 HS25_09827:2:1201:1559:70726#49        <NA>
    ## 3 HS25_09827:2:1201:1564:39627#49        <NA>
    ## 4 HS25_09827:2:1201:1565:91731#49        <NA>
    ## 5 HS25_09827:2:1201:1624:69925#49        <NA>
    ##                                                                                               SEQUENCE
    ## 1 CCGTTAGAGCATTTGTTGAAAATGCTTTCCTTGCTCCATGTGATGACTCTGGTGCCCTTGTCAAAAGCCAGCTGGGCCTATTCGTGTGGGTCTGTTTCTG
    ## 2 TTGTTAAAATGACCATACCCAAAGTGATCTACAGACTCAATACAATTTCTATTGAAATACCAATCACACTCTTCACAGAACTAGAAAAACAGTTCTAAAA
    ## 3 ACGCGGCAATCCAATGTGTGAGTTGAGAAGCGGTGAGGAGGGAATCCTAATTTTATGAGCAGGTCAGGACCGTGGGAGATACCTGACACCTGAGATGGTA
    ## 4 GACATGCCATAACATTCATGTTTTATGTGTACAAGTCAATGAATTTTAGTATATTTACAGAGTTGTATGACTGTCTCCACAATCTAATTTTAGGTTTCCA
    ## 5 GCCAGCCTCCTTCTCAATGGTCTTTTTAAACATTATATGAAAACCAGACATTTACATTTGATTTCTTTTTCAATACTATACAGTTCTAAGAGAAAAAACA
    ##                                                                                                QUALITY
    ## 1 CABCFGDEEFFEFHGHGGFFGDIGIJFIFHHGHEIFGHBCGHDIFBE9GIAICGGICFIBFGGHGDGGGHE?GIGDFGGHEGIEJG>;FG<GGHACEFGH
    ## 2 CABEFGFFGFHGGGGJGGFFGKIHHJFIEHHHGIEGGEHJGHDHFGHIGICIJEFIFGIF8GGHKFHGGFEI6GGGFIGHGGIE>EFCFHGGGHEJEAJE
    ## 3 BACCFGBFGFHGGJGHGGFEGHIGIJHFEH:HHEHGHHBGGH9IAGHGFHIFJFFAFGIFDIGHKEIG<C>F,CGD66?7EFI5EEG>EGGGGD5=HH6E
    ## 4 CABFFGFFJFHEGEGJGGDG?FIGHHHBGHHHGIIGHGHGGHDGHFHIDFCIKEGIFHGGII9HFFGGGEEIGGEEHGGEEGDEHFH>FGGGGHAFAHGE
    ## 5 CABEFGFGIFGGGJGHGGFH?FDHGHDHGHEHHJCGHHFHDHDHFGHIGHIFFHGHFGGGI9GHF@IGGH;FICGEFEIHGGIEEFC:DEGGGBDJHHFF
    ##   MATE                         PAIR_ID
    ## 1    1 HS25_09827:2:1201:1505:59795#49
    ## 2    1 HS25_09827:2:1201:1559:70726#49
    ## 3    1 HS25_09827:2:1201:1564:39627#49
    ## 4    1 HS25_09827:2:1201:1565:91731#49
    ## 5    1 HS25_09827:2:1201:1624:69925#49

### FASTQ quality control

FASTQ quality handling has two separate knobs:

- `input_quality_encoding` controls how incoming FASTQ ASCII is decoded.
  The default is modern `phred33`; use `phred64`, `solexa64`, or `auto`
  only for legacy data.
- `quality_representation` controls how qualities are returned to
  DuckDB: canonical Phred+33 text (`"string"`) or numeric Phred arrays
  (`"phred"`).

The flow is:

1.  Decode FASTQ ASCII using `input_quality_encoding`.
2.  Normalize to numeric Phred qualities.
3.  Return either numeric arrays or canonical Phred+33 text.

`BAM`/`CRAM` reads skip the text-decoding step because qualities are
already stored numerically.

For ordinary QC, aggregate the sequence and canonical quality strings
directly. This returns exact global totals plus a small nested per-cycle
relation without generating one SQL row per base.

``` r
fastq_qc <- dbGetQuery(
  con,
  sprintf(
    "WITH q AS (
       SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc
       FROM read_fastq('%s')
     )
     SELECT qc.reads, qc.bases, qc.q30_bases, qc.max_read_length
     FROM q",
    fastq_r1
  )
)
fastq_qc
```

    ##   reads bases q30_bases max_read_length
    ## 1     5   500       475             100

Use numeric quality arrays when the query genuinely needs the full
per-position histogram:

``` r
legacy_fastq <- system.file("extdata", "legacy_phred64.fq", package = "Rduckhts")

rduckhts_detect_quality_encoding(con, legacy_fastq)
```

    ##   format observed_ascii_min observed_ascii_max records_sampled
    ## 1  fastq                104                104               1
    ##       compatible_encodings guessed_encoding is_ambiguous
    ## 1 phred33,phred64,solexa64          phred64         TRUE

``` r
quality_hist <- dbGetQuery(
  con,
  sprintf(
    "WITH q AS (
       SELECT NAME, QUALITY
       FROM read_fastq('%s', quality_representation := 'phred')
     ),
     expanded AS (
       SELECT
         NAME,
         generate_subscripts(QUALITY, 1) AS pos,
         unnest(QUALITY) AS q
       FROM q
     )
     SELECT pos, q AS phred, count(*) AS n_reads
     FROM expanded
     GROUP BY pos, phred
     ORDER BY pos, phred
     LIMIT 12",
    fastq_r1
  )
)
quality_hist
```

    ##    pos phred n_reads
    ## 1    1    33       1
    ## 2    1    34       4
    ## 3    2    32       5
    ## 4    3    33       4
    ## 5    3    34       1
    ## 6    4    34       2
    ## 7    4    36       2
    ## 8    4    37       1
    ## 9    5    37       5
    ## 10   6    38       5
    ## 11   7    33       1
    ## 12   7    35       1

### FASTA region queries

`read_fasta` supports indexed region queries via
`rduckhts_fasta(..., region = ...)`.

``` r
fai_path <- file.path(tempdir(), "duckhts_readme_00000000000000.fai")
fai_info <- rduckhts_fasta_index(con, fasta_path, index_path = fai_path)
fai_info
```

    ##   success                                        index_path
    ## 1    TRUE <tempfile>

``` r
rduckhts_fasta(
  con, "fasta_region", fasta_path,
  region = "CHROMOSOME_I:1-25",
  overwrite = TRUE
)
dbGetQuery(con, "SELECT NAME, length(SEQUENCE) AS n FROM fasta_region")
```

    ##           NAME  n
    ## 1 CHROMOSOME_I 25

``` r
unlink(fai_path)
```

### Sequence functions

The extension also exposes sequence utility UDFs directly in DuckDB SQL,
including 4-bit IUPAC DNA encode/decode helpers. These can be applied to
`SEQUENCE` columns from FASTA and FASTQ scans.

``` r
dbGetQuery(
  con,
  "SELECT
     NAME,
     seq_hash_2bit(substr(SEQUENCE, 1, 12)) AS hash_2bit_prefix,
     seq_encode_4bit(substr(SEQUENCE, 1, 16)) AS codes,
     seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 16))) AS roundtrip
   FROM sequences
   LIMIT 2"
)
```

    ##            NAME hash_2bit_prefix                                          codes
    ## 1  CHROMOSOME_I          9898352 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8
    ## 2 CHROMOSOME_II          6038978 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1
    ##          roundtrip
    ## 1 GCCTAAGCCTAAGCCT
    ## 2 CCTAAGCCTAAGCCTA

``` r
dbGetQuery(
  con,
  "SELECT
     NAME,
     MATE,
     seq_encode_4bit(substr(SEQUENCE, 1, 12)) AS codes,
     seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 12))) AS roundtrip
   FROM reads
   LIMIT 2"
)
```

    ##                              NAME MATE                              codes
    ## 1 HS25_09827:2:1201:1505:59795#49    1 2, 2, 4, 8, 8, 1, 4, 1, 4, 2, 1, 8
    ## 2 HS25_09827:2:1201:1505:59795#49    2 1, 1, 4, 4, 1, 1, 1, 4, 1, 1, 4, 4
    ##      roundtrip
    ## 1 CCGTTAGAGCAT
    ## 2 AAGGAAAGAAGG

## Annotation, intervals and indexes

### GFF and GTF attributes

GFF3 files are read with
[`rduckhts_gff()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_gff.md)
/ SQL `read_gff(...)`; GTF files are read with
[`rduckhts_gtf()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_gtf.md)
/ SQL `read_gtf(...)`. `strict = TRUE` enables GFF3 structural
validation. Attribute decoding can be scalar and raw for legacy
convenience (`attributes_map`), grouped and lossless for multi-values
(`attributes_list`, a DuckDB `MAP(VARCHAR, VARCHAR[])`), or exact
parser-style pairs (`attributes_pairs`, a DuckDB
`LIST<STRUCT(key, value, idx)>`).

The extension-level GFF3 implementation is benchmarked and audited
against [GFFBase](https://github.com/Kuanhao-Chao/gffbase) in the
DuckHTS repo:
<https://github.com/RGenomicsETL/duckhts/blob/develop/benchmarks/benchmark_gffbase_conformance.md>.

``` r
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "genes", gff_path, attributes_map = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT seqname, start, \"end\" FROM genes WHERE feature = 'gene' LIMIT 5")
```

    ##   seqname   start     end
    ## 1       X 2934816 2964270

``` r
gff_attrs_path <- system.file("extdata", "gff_attrs.gff3", package = "Rduckhts")
rduckhts_gff(
  con,
  "gff_attrs",
  gff_attrs_path,
  strict = TRUE,
  attributes_list = TRUE,
  attributes_pairs = TRUE,
  overwrite = TRUE
)
dbGetQuery(con, paste(
  "SELECT seqname, feature,",
  "list_extract(map_extract_value(attributes_list, 'Dbxref'), 1) AS first_dbxref,",
  "list_count(attributes_pairs) AS n_attr_pairs FROM gff_attrs"
))
```

    ##   seqname feature first_dbxref n_attr_pairs
    ## 1    chr1    gene     GeneID:1            6

``` r
gtf_attrs_path <- system.file("extdata", "gtf_attrs.gtf", package = "Rduckhts")
rduckhts_gtf(con, "gtf_attrs", gtf_attrs_path, attributes_list = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT list_extract(map_extract_value(attributes_list, 'note'), 1) AS note FROM gtf_attrs")
```

    ##          note
    ## 1 weird; semi

### BigWig signals

The package bundles libBigWig’s upstream test file so the wrapper,
stored zero-based half-open coordinates, and htslib-style multi-region
semantics are executable offline.
[`rduckhts_bigwig()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_bigwig.md)
materializes a table; pass `table_name = NULL` to create the
`bigwig_data` view instead.

``` r
bigwig_path <- system.file(
  "extdata", "libbigwig_test.bw", package = "Rduckhts", mustWork = TRUE
)
rduckhts_bigwig(
  con, "bigwig_intervals", bigwig_path,
  region = c("1:1-150", "10:201-300"),
  overwrite = TRUE
)
dbGetQuery(
  con,
  paste(
    "SELECT CHROM, START0, END0, round(VALUE::DOUBLE, 1) AS VALUE",
    "FROM bigwig_intervals ORDER BY CHROM, START0"
  )
)
```

    ##   CHROM START0 END0 VALUE
    ## 1     1      0    1   0.1
    ## 2     1      1    2   0.2
    ## 3     1      2    3   0.3
    ## 4     1    100  150   1.4
    ## 5    10    200  300   2.0

### Interval and reference helpers

``` r
bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
fai_path <- file.path(tempdir(), "duckhts_readme_00000000000000.fai")
rduckhts_fasta_index(con, fasta_path, index_path = fai_path)
```

    ##   success                                        index_path
    ## 1    TRUE <tempfile>

``` r
rduckhts_bed(con, "targets", bed_path, overwrite = TRUE)
dbGetQuery(con, "SELECT chrom, start, \"end\", name, block_count FROM targets")
```

    ##            chrom start end    name block_count
    ## 1   CHROMOSOME_I     0  10 target1           2
    ## 2   CHROMOSOME_I    10  20 target2           1
    ## 3  CHROMOSOME_II     0   8 target3          NA
    ## 4 CHROMOSOME_III     0   6 target4           1

``` r
rduckhts_fasta_nuc(con, fasta_path, bed_path = bed_path, index_path = fai_path)
```

    ##            chrom start end pct_at pct_gc num_a num_c num_g num_t num_n
    ## 1   CHROMOSOME_I     0  10  0.400  0.600     2     4     2     2     0
    ## 2   CHROMOSOME_I    10  20  0.500  0.500     4     3     2     1     0
    ## 3  CHROMOSOME_II     0   8  0.375  0.625     2     4     1     1     0
    ## 4 CHROMOSOME_III     0   6  0.500  0.500     2     2     1     1     0
    ##   num_other seq_len
    ## 1         0      10
    ## 2         0      10
    ## 3         0       8
    ## 4         0       6

``` r
rduckhts_fasta_nuc(con, fasta_path, bin_width = 10, region = "CHROMOSOME_I:1-20", index_path = fai_path)
```

    ##          chrom start end pct_at pct_gc num_a num_c num_g num_t num_n num_other
    ## 1 CHROMOSOME_I     0  10    0.4    0.6     2     4     2     2     0         0
    ## 2 CHROMOSOME_I    10  20    0.5    0.5     4     3     2     1     0         0
    ##   seq_len
    ## 1      10
    ## 2      10

``` r
unlink(fai_path)
```

### cgranges interval indexes

The bundled extension also exposes SQL-first `duckhts_cgranges_*` entry
points. These are session-scoped interval indexes that you can populate
either row-wise or in bulk from a table or view with
`duckhts_cgranges_from_table(...)`, which runs on your connection and so
also sees TEMP objects, then query through
`duckhts_cgranges_overlaps(...)`. For row-preserving filters or count
annotations over provider rows, use the vectorized scalar helpers
`duckhts_cgranges_has_overlap(...)` and
`duckhts_cgranges_count_overlaps(...)` directly in queries over
`read_bed(...)`, `read_bam(...)`, `read_bcf(...)`, or regular tables.
For streaming one-row-per-hit expansion while keeping provider columns,
use `duckhts_cgranges_overlaps_list(...)` with `UNNEST(...)` in the
SELECT list, which also covers bulk probing of any relation. There is no
dedicated R wrapper yet, so use them through `DBI`.

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_create('readme_idx') AS ok")
```

    ##     ok
    ## 1 TRUE

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('readme_idx', 'chr1', 10, 20, 'a') AS ok")
```

    ##     ok
    ## 1 TRUE

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('readme_idx', 'chr1', 30, 40, 'b') AS ok")
```

    ##     ok
    ## 1 TRUE

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('readme_idx') AS ok")
```

    ##     ok
    ## 1 TRUE

``` r
DBI::dbGetQuery(
  con,
  paste(
    "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
    "FROM duckhts_cgranges_overlaps('readme_idx', 'chr1', 35, 36, query_row_id := 7)"
  )
)
```

    ##   interval_ordinal label interval_chrom interval_start interval_end
    ## 1                1     b           chr1             30           40

``` r
DBI::dbExecute(
  con,
  paste(
    "CREATE TEMP VIEW readme_targets AS SELECT * FROM (VALUES",
    "('chr2', 100, 110, 'alpha'), ('chr2', 150, 170, 'beta')",
    ") AS t(chrom, start, \"end\", label)"
  )
)
```

    ## [1] 0

``` r
DBI::dbGetQuery(
  con,
  paste(
    "SELECT * FROM duckhts_cgranges_from_table(",
    "'readme_qry_idx', 'readme_targets', 'chrom', 'start', 'end', 'label')"
  )
)
```

    ##   indexed
    ## 1    TRUE

``` r
DBI::dbGetQuery(
  con,
  paste(
    "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
    "FROM duckhts_cgranges_overlaps('readme_qry_idx', 'chr2', 140, 170, mode := 'contain')"
  )
)
```

    ##   interval_ordinal label interval_chrom interval_start interval_end
    ## 1                1  beta           chr2            150          170

``` r
DBI::dbExecute(
  con,
  paste(
    "CREATE TABLE readme_probes AS SELECT * FROM (VALUES",
    "(10, 'chr2', 100, 105),",
    "(20, 'chr2', 160, 161),",
    "(30, 'chr2', 500, 510)",
    ") AS t(probe_id, chrom, start, \"end\")"
  )
)
```

    ## [1] 3

``` r
DBI::dbGetQuery(
  con,
  paste(
    "SELECT probe_id, hit.interval_ordinal, hit.label, hit.label_type,",
    "  hit.interval_chrom, hit.interval_start, hit.interval_end",
    "FROM (",
    "  SELECT p.probe_id,",
    "    unnest(duckhts_cgranges_overlaps_list('readme_qry_idx', p.chrom, p.start, p.\"end\")) AS hit",
    "  FROM readme_probes AS p",
    ")",
    "ORDER BY probe_id, hit.interval_ordinal"
  )
)
```

    ##   probe_id interval_ordinal label label_type interval_chrom interval_start
    ## 1       10                0 alpha    VARCHAR           chr2            100
    ## 2       20                1  beta    VARCHAR           chr2            150
    ##   interval_end
    ## 1          110
    ## 2          170

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_destroy('readme_idx') AS ok")
```

    ##     ok
    ## 1 TRUE

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_destroy('readme_qry_idx') AS ok")
```

    ##     ok
    ## 1 TRUE

### Tabix headers and column types

Use `header = TRUE` to use the first non-meta row as column names, and
`auto_detect = TRUE` / `column_types` to control column typing:

``` r
tabix_header <- system.file("extdata", "header_tabix.tsv.gz", package = "Rduckhts")
tabix_meta <- system.file("extdata", "meta_tabix.tsv.gz", package = "Rduckhts")

rduckhts_tabix(con, "header_tabix", tabix_header, header = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT chrom, pos FROM header_tabix LIMIT 2")
```

    ##   chrom pos
    ## 1  chr1   1
    ## 2  chr1   2

``` r
rduckhts_tabix(con, "typed_tabix", tabix_meta, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT typeof(column1) AS column1_type FROM typed_tabix LIMIT 1")
```

    ##   column1_type
    ## 1       BIGINT

``` r
rduckhts_tabix(con, "typed_tabix_explicit", tabix_header,
               header = TRUE,
               column_types = c("VARCHAR", "BIGINT", "VARCHAR"),
               overwrite = TRUE)
dbGetQuery(con, "SELECT pos + 1 AS pos_plus_one FROM typed_tabix_explicit LIMIT 1")
```

    ##   pos_plus_one
    ## 1            2

### Remote tabix: GTEx eQTL

GTEx eQTL matrices on EBI are tabix-indexed. In browser wasm/webR, this
depends on CORS policy on both the data object and index object. This
remote example is an unevaluated usage snippet.

``` r
gtex_url <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz" 
rduckhts_tabix(con, "gtex_eqtl", gtex_url, region = "1:11868-14409",
                  header = TRUE, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM gtex_eqtl LIMIT 5")
```

### Compression and indexing

``` r
bed_src <- system.file("extdata", "targets.bed", package = "Rduckhts")
bam_src <- system.file("extdata", "range.bam", package = "Rduckhts")
bcf_src <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")

tmp_bed <- tempfile("duckhts_targets_", fileext = ".bed")
tmp_bgz <- paste0(tmp_bed, ".gz")
tmp_tbi <- paste0(tmp_bgz, ".tbi")
tmp_roundtrip <- tempfile("duckhts_targets_roundtrip_", fileext = ".bed")
tmp_bai <- tempfile("duckhts_range_", fileext = ".bam.bai")
tmp_csi <- tempfile("duckhts_variants_", fileext = ".bcf.csi")
file.copy(bed_src, tmp_bed, overwrite = TRUE)
```

    ## [1] TRUE

``` r
bgzip_meta <- rduckhts_bgzip(
  con, tmp_bed,
  output_path = tmp_bgz,
  threads = 1,
  keep = TRUE,
  overwrite = TRUE
)
transform(bgzip_meta[, c("success", "output_path", "bytes_out")],
          output_path = "<tempfile>")
```

    ##   success output_path bytes_out
    ## 1    TRUE  <tempfile>       169

``` r
bgunzip_meta <- rduckhts_bgunzip(
  con, tmp_bgz,
  output_path = tmp_roundtrip,
  threads = 1,
  keep = TRUE,
  overwrite = TRUE
)
bgunzip_meta$output_path <- "<tempfile>"
bgunzip_meta[, c("success", "output_path", "bytes_out")]
```

    ##   success output_path bytes_out
    ## 1    TRUE  <tempfile>       194

``` r
bam_index_meta <- rduckhts_bam_index(
  con, bam_src,
  index_path = tmp_bai,
  threads = 1
)
transform(bam_index_meta, index_path = "<tempfile>")
```

    ##   success index_path index_format
    ## 1    TRUE <tempfile>          BAI

``` r
bcf_index_meta <- rduckhts_bcf_index(
  con, bcf_src,
  index_path = tmp_csi,
  threads = 1
)
transform(bcf_index_meta, index_path = "<tempfile>")
```

    ##   success index_path index_format
    ## 1    TRUE <tempfile>          CSI

``` r
tabix_meta <- rduckhts_tabix_index(
  con, tmp_bgz,
  preset = "bed",
  index_path = tmp_tbi,
  threads = 1
)
transform(tabix_meta, index_path = "<tempfile>")
```

    ##   success index_path index_format
    ## 1    TRUE <tempfile>          TBI

``` r
rduckhts_bed(con, "targets_idx", tmp_bgz, region = "CHROMOSOME_I:1-20", index_path = tmp_tbi, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM targets_idx")
```

    ##          chrom start end    name score strand thick_start thick_end item_rgb
    ## 1 CHROMOSOME_I     0  10 target1   100      +           0        10  255,0,0
    ## 2 CHROMOSOME_I    10  20 target2   200      -          10        20  0,0,255
    ##   block_count block_sizes block_starts extra
    ## 1           2         5,5          0,5  <NA>
    ## 2           1          10            0  <NA>

``` r
unlink(c(tmp_bed, tmp_bgz, tmp_tbi, tmp_roundtrip, tmp_bai, tmp_csi))
```

### HTS header and index metadata

Use metadata helpers to inspect parsed headers, raw header lines, index
summaries, span-oriented index views, and raw index blobs.

``` r
header_meta <- rduckhts_hts_header(con, bcf_path)
head(header_meta[, c("record_type", "id", "number", "value_type")], 5)
```

    ##   record_type   id number value_type
    ## 1  fileformat <NA>   <NA>       <NA>
    ## 2      FILTER PASS   <NA>       <NA>
    ## 3        INFO TEST      1    Integer
    ## 4      FORMAT   TT      A    Integer
    ## 5        INFO  DP4      4    Integer

``` r
header_raw <- rduckhts_hts_header(con, bcf_path, mode = "raw")
head(header_raw[, c("idx", "raw")], 5)
```

    ##   idx
    ## 1   0
    ## 2   1
    ## 3   2
    ## 4   3
    ## 5   4
    ##                                                                                                                                                          raw
    ## 1                                                                                                                                       ##fileformat=VCFv4.1
    ## 2                                                                                                        ##FILTER=<ID=PASS,Description="All filters passed">
    ## 3                                                                                           ##INFO=<ID=TEST,Number=1,Type=Integer,Description="Testing Tag">
    ## 4 ##FORMAT=<ID=TT,Number=A,Type=Integer,Description="Testing Tag, with commas and \\"escapes\\" and escaped escapes combined with \\\\\\"quotes\\\\\\\\\\"">
    ## 5                       ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">

``` r
index_meta <- rduckhts_hts_index(con, bcf_path, index_path = bcf_index_path)
head(index_meta[, c("seqname", "mapped", "unmapped", "index_type")], 5)
```

    ##   seqname mapped unmapped index_type
    ## 1       1     11        0        CSI
    ## 2       2      1        0        CSI
    ## 3       3      1        0        CSI
    ## 4       4      2        0        CSI

``` r
index_spans <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
head(index_spans[, c("seqname", "tid", "index_type", "chunk_beg_vo", "chunk_end_vo")], 5)
```

    ##   seqname tid index_type chunk_beg_vo chunk_end_vo
    ## 1       1   0        CSI         1586         1713
    ## 2       1   0        CSI         1713         1973
    ## 3       1   0        CSI         1973         2109
    ## 4       1   0        CSI         2109         2242
    ## 5       1   0        CSI         2242         2372

``` r
index_raw <- rduckhts_hts_index_raw(con, bcf_path, index_path = bcf_index_path)
head(index_raw, 1)
```

    ## [1] index_type
    ## [2] '<Rduckhts>/extdata/vcf_file.bcf.csi'
    ## [3] raw
    ## <0 rows> (or 0-length row.names)

### Reading many files

The `rduckhts_*_multi` family reads multiple files into a single DuckDB
table with a `filename` column, following the same
`(con, table_name, ...)` convention as the single-file wrappers:

``` r
fq_files <- c(
  system.file("extdata", "r1.fq", package = "Rduckhts"),
  system.file("extdata", "r2.fq", package = "Rduckhts")
)
rduckhts_fastq_multi(con, "fq_multi", fq_files, overwrite = TRUE)
dbGetQuery(con, "SELECT parse_filename(filename) AS file, count(*) AS n FROM fq_multi GROUP BY ALL ORDER BY file")
```

    ##    file n
    ## 1 r1.fq 5
    ## 2 r2.fq 5

Per-file parameters are supported via a `.params` data.frame with a
`file` column and columns matching reader arguments. `NA` values fall
back to the uniform default:

``` r
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_idx  <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

params <- data.frame(
  file       = bam_path,
  region     = "CHROMOSOME_I:1-1000",
  index_path = bam_idx
)
rduckhts_bam_multi(con, "bam_multi", bam_path, .params = params,
                   overwrite = TRUE)
dbGetQuery(con, "SELECT count(*) AS n FROM bam_multi")
```

    ##   n
    ## 1 2

## Diagnostics and linking

### SIMD backends

The bundled extension uses explicit runtime SIMD dispatch for
byte-oriented helper kernels, starting with `seq_gc_content(...)`.
`scalar` is always available and is the portable baseline. Optional
platform backends such as `avx2` or `avx512` should be checked with
[`rduckhts_simd_backend_available()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_simd_backend.md)
before being requested. The `auto` policy resolves each logical kernel
independently from the current compiled-and-CPU-supported capability
mask; use
[`rduckhts_simd_kernel_info()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_simd_backend.md)
for the per-kernel result and `rduckhts_simd_set_backend(con, "auto")`
to return to runtime auto-detection.

``` r
rduckhts_simd_info(con)[, c("backend", "selectable", "compiled", "cpu_supported", "available", "selected")]
```

    ##        backend selectable compiled cpu_supported available selected
    ## 1       scalar       TRUE     TRUE          TRUE      TRUE    FALSE
    ## 2         sse2      FALSE    FALSE          TRUE     FALSE    FALSE
    ## 3        sse41      FALSE    FALSE          TRUE     FALSE    FALSE
    ## 4         avx2       TRUE     TRUE          TRUE      TRUE     TRUE
    ## 5       avx512       TRUE     TRUE         FALSE     FALSE    FALSE
    ## 6         neon       TRUE    FALSE         FALSE     FALSE    FALSE
    ## 7 wasm_simd128       TRUE    FALSE         FALSE     FALSE    FALSE

``` r
rduckhts_simd_kernel_info(con)[, c("kernel", "selected_backend", "scalar_fallback")]
```

    ##            kernel selected_backend scalar_fallback
    ## 1 seq_base_counts             avx2           FALSE
    ## 2 bam_nt16_counts             avx2           FALSE
    ## 3  nt16_gc_counts             avx2           FALSE
    ## 4        fastq_qc             avx2           FALSE

``` r
rduckhts_simd_set_backend(con, "scalar")
```

    ## [1] "scalar"

``` r
DBI::dbGetQuery(
  con,
  paste(
    "SELECT duckhts_simd_requested_backend() AS requested_backend,",
    "duckhts_simd_backend() AS selected_backend,",
    "round(seq_gc_content('ACGTNNacgtnn'), 3) AS gc_content"
  )
)
```

    ##   requested_backend selected_backend gc_content
    ## 1            scalar           scalar        0.5

``` r
data.frame(
  requested_backend = rduckhts_simd_requested_backend(con),
  selected_backend = rduckhts_simd_backend(con)
)
```

    ##   requested_backend selected_backend
    ## 1            scalar           scalar

``` r
restored_backend <- rduckhts_simd_set_backend(con, "auto")
data.frame(
  requested_backend = rduckhts_simd_requested_backend(con),
  selected_backend_known = nzchar(restored_backend)
)
```

    ##   requested_backend selected_backend_known
    ## 1              auto                   TRUE

### Linking against the bundled htslib

[`rduckhts_htslib_config()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_htslib_config.md)
returns the exact installed headers, shared or static library, required
flags, enabled features, and build identity. With no `link` argument it
selects the mode chosen when Rduckhts was configured. Validation loads
DuckHTS and rejects a source/header/runtime version mismatch. A
downstream `configure` script can emit its `PKG_CPPFLAGS` and `PKG_LIBS`
from this one receipt.

``` r
hts_config <- rduckhts_htslib_config()
hts_config[c("contract_version", "htslib_version", "runtime_version", "link")]
```

    ## $contract_version
    ## [1] 1
    ##
    ## $htslib_version
    ## [1] "1.24"
    ##
    ## $runtime_version
    ## [1] "1.24"
    ##
    ## $link
    ## [1] "shared"

``` r
hts_config$features
```

    ## $cram
    ## [1] TRUE
    ##
    ## $zlib
    ## [1] TRUE
    ##
    ## $bzip2
    ## [1] TRUE
    ##
    ## $lzma
    ## [1] TRUE
    ##
    ## $libdeflate
    ## [1] TRUE
    ##
    ## $curl
    ## [1] TRUE
    ##
    ## $openssl
    ## [1] TRUE
    ##
    ## $plugins
    ## [1] TRUE
    ##
    ## $s3
    ## [1] TRUE
    ##
    ## $gcs
    ## [1] TRUE
