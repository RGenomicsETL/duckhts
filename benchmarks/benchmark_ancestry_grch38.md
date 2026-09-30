Ancestry projection: GRCh38 reference product and 30x CRAM acceptance
================

The GRCh38 ancestry reference is derived from the keyed GRCh37 product
by `duckdb_liftover` and checked here against the GRCh37 phased
genotypes of the same 1000 Genomes individuals. Revision:
9eeb223eb30d051c515e72e361667aac0145c2f1.

## GRCh38 reference product

`duckhts_bench_stage_ancestry_grch38_parquet()` lifts each of the
5,816,590 loci of the keyed GRCh37 Parquet with the registered Ensembl
GRCh37-to-GRCh38 chain, the registered `human_g1k_v37` source FASTA and
the registered Broad `Homo_sapiens_assembly38` destination FASTA. The
coded allele of a locus is `allele_b`, the allele that the 21 group
frequencies and 16 PC loadings describe. Liftover reports its
destination spelling through `reverse_complemented`; the derived key is
`allele_a` = the other allele and `allele_b` = the coded allele, and
every frequency and loading is copied unchanged. This is the reversal
convention of `rduckhts_ancestry_proportions()`: a reversed input meets
an unchanged reference row. Rewriting a swapped locus as `1 - f` with
negated loadings is not equivalent once the per-PC correction
coefficients apply (the constant shift of `-sum(U)` in each PC is
multiplied by the coefficient on one side only), so it is not done. The
output is sorted by `(chromosome, position, allele_a, allele_b)` with
the same columns, types and row groups as the GRCh37 product, so the
panel builder and the proportion wrapper consume it unchanged. A
liftover map Parquet retains the source locus of every output locus.

| quantity                              |      loci |
|:--------------------------------------|----------:|
| input_loci                            | 5,816,590 |
| mapped                                | 5,805,811 |
| rejected                              |   803,764 |
| rejected: destination_allele_mismatch |       110 |
| rejected: liftover_ApartAnchors       |         2 |
| rejected: liftover_UnmappedAnchors    |    10,777 |
| rejected: non_snv_source              |   792,875 |
| reverse_complemented                  |     5,863 |
| swapped                               |    14,843 |
| duplicate_destination_dropped         |         0 |
| output_loci                           | 5,012,826 |

`mapped` counts loci for which liftover reported a destination, before
the later classification. Every input locus is in exactly one of
rejected, duplicate-destination dropped or output (803764 + 0 + 5012826
= 5816590). `swapped` and `reverse_complemented` count output loci.
Non-single-base source loci are rejected because ancestry panels select
biallelic SNVs and an indel has no coded allele that survives
left-alignment. The receipt binds these counts to the source Parquet,
chain, both FASTAs, the output and the DuckDB v1.5.5, htslib 1.24 and
Rduckhts 1.5.2.9007.0.1.5 versions:

| field                    | value                                                            |
|:-------------------------|:-----------------------------------------------------------------|
| source_parquet_sha256    | a0bf01ae16db605ff0add59e5ebb747c35454406900c11fd8596f23c167fe02b |
| chain_sha256             | 351de3cd4a01d9fcffd38881981767b697090d2eba876740891b96d5c546b100 |
| source_fasta_sha256      | 2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f |
| destination_fasta_sha256 | 93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564 |
| output_sha256            | 4982810aa69b9c71924e9a02dcd454ccedccc34bd8bcf71506ce9752589847fa |
| liftover_map_sha256      | 3d5975331ee4f49a687ad081ee60ac83da3fe8bbee7721ab6b1de49e1a608c56 |

## Acceptance comparison

The panel is the 1,000-site genome-wide output of
`rduckhts_ancestry_panel()` on the GRCh38 product (5,000 bp windows,
GRCh38 assembly); its SHA-256 is
bd232e576a3889c81550b36c9670ba2d8df614d50e07fc0973c15bc3df8022d9,
verified equal to the staging receipt. Each panel site’s GRCh37 source
locus comes from the liftover map. The GRCh37 side reads the public
phase-3 VCFs (`ALL.chrN.phase3_shapeit2_mvncall_integrated_v5b`) by
indexed `read_bcf(region := ...)` at only those loci, for HG00188,
HG00403, NA18507; keeping REF/ALT that match the source alleles in
either orientation left 3000 sample-locus rows of 3000 possible. The
GRCh38 side is the registered 30x CRAM of each individual, cut by
`samtools view -M -L` to the reads that overlap a panel site and indexed
(`<sample>.panel-sites.cram`), because the whole-genome CRAMs are 14 to
16 GB. Both stagings have a receipt with remote byte counts, ETags and
output SHA-256 values.

Declared tolerance (fixed in this report’s source before any result):
the maximum absolute difference over the 21 groups is at most 0.05 for
each CRAM frequency method (`allele_fraction`, `called_genotype`), both
sides pass the correlation gate (`cor_pred` at least 0.4), and the two
largest groups agree. `min_depth` is 7. Per sample, matched sites are
the panel loci that survive allele matching; gate failures are retained
as `status` and would make the difference undefined.

| sample_id | superpopulation | method          | panel_sites | genotype_matched | cram_matched | genotype_cor | cram_cor | genotype_gate | cram_gate | max_abs_difference | top2_agree | within_tolerance | derivation_matched | derivation_difference |
|:----------|:----------------|:----------------|------------:|-----------------:|-------------:|-------------:|---------:|:--------------|:----------|-------------------:|:-----------|:-----------------|-------------------:|----------------------:|
| HG00188   | EUR             | allele_fraction |        1000 |             1000 |          991 |       0.6860 |   0.6791 | ok            | ok        |             0.0297 | TRUE       | TRUE             |               1000 |                     0 |
| HG00188   | EUR             | called_genotype |        1000 |             1000 |          962 |       0.6860 |   0.6860 | ok            | ok        |             0.0314 | TRUE       | TRUE             |               1000 |                     0 |
| HG00403   | EAS             | allele_fraction |        1000 |             1000 |          992 |       0.7456 |   0.7356 | ok            | ok        |             0.0665 | TRUE       | FALSE            |               1000 |                     0 |
| HG00403   | EAS             | called_genotype |        1000 |             1000 |          968 |       0.7456 |   0.7444 | ok            | ok        |             0.0381 | TRUE       | TRUE             |               1000 |                     0 |
| NA18507   | AFR             | allele_fraction |        1000 |             1000 |          993 |       0.7200 |   0.7192 | ok            | ok        |             0.0000 | TRUE       | TRUE             |               1000 |                     0 |
| NA18507   | AFR             | called_genotype |        1000 |             1000 |          968 |       0.7200 |   0.7213 | ok            | ok        |             0.0025 | FALSE      | FALSE            |               1000 |                     0 |

Gate failures: 0 genotype and 0 CRAM of 6 comparisons. Comparisons
within the declared tolerance: 4 of 6. Largest difference: 0.0665
against the 0.05 tolerance, so the declared acceptance is **not met**
for 2 comparisons: HG00403 allele_fraction; NA18507 called_genotype. The
declared criterion is reported as written and was not adjusted after the
result. The failing rows differ in kind. A maximum difference above the
tolerance (0.0664682 for HG00403) is a real disagreement between the
read allele fractions and the phased genotypes for the East Asian groups
(Asia (East), Japan, Philippines and Sri Lanka absorb each other’s
share). The NA18507 `called_genotype` failure is a property of the
`top2_agree` rule: the GRCh37 result has one non-zero group, so its
second-largest group is a 20-way tie at 0, and the CRAM result’s 0.0025
in Finland changes the second name while the maximum difference is
0.0025.

The derivation control (the same genotypes, keyed by their GRCh38 locus
and destination alleles, against the GRCh38 product) matches the GRCh37
result to 0 at 1000 matched loci for every sample (6 of the panel loci
have swapped allele roles), so the differences in the table come from
comparing 30x read counts with the phased genotypes, not from the
liftover, the swapped or reverse-complemented allele roles, or the
sorted GRCh38 layout.

Proportions, group by group (GRCh37 phased genotypes against the GRCh37
product versus the GRCh38 CRAM against the GRCh38 product):

| sample_id | group_id            | genotype | cram_allele_fraction | cram_called_genotype |
|:----------|:--------------------|---------:|---------------------:|---------------------:|
| HG00188   | Finland             |   0.9437 |               0.9140 |               0.9123 |
| HG00188   | Japan               |   0.0563 |               0.0705 |               0.0877 |
| HG00188   | Africa (East)       |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Africa (North)      |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Africa (South)      |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Africa (West)       |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Ashkenazi           |   0.0000 |               0.0155 |               0.0000 |
| HG00188   | Asia (East)         |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Bangladesh          |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Europe (North East) |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Europe (South East) |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Europe (South West) |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Ireland             |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Italy               |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Middle East         |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Pakistan            |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Philippines         |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Scandinavia         |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | South America       |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | Sri Lanka           |   0.0000 |               0.0000 |               0.0000 |
| HG00188   | United Kingdom      |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Asia (East)         |   0.8421 |               0.7773 |               0.8041 |
| HG00403   | Japan               |   0.1414 |               0.1029 |               0.1309 |
| HG00403   | Philippines         |   0.0157 |               0.0822 |               0.0394 |
| HG00403   | Finland             |   0.0008 |               0.0124 |               0.0000 |
| HG00403   | Africa (East)       |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Africa (North)      |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Africa (South)      |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Africa (West)       |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Ashkenazi           |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Bangladesh          |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Europe (North East) |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Europe (South East) |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Europe (South West) |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Ireland             |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Italy               |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Middle East         |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Pakistan            |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Scandinavia         |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | South America       |   0.0000 |               0.0000 |               0.0000 |
| HG00403   | Sri Lanka           |   0.0000 |               0.0253 |               0.0256 |
| HG00403   | United Kingdom      |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Africa (West)       |   1.0000 |               1.0000 |               0.9975 |
| NA18507   | Africa (East)       |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Africa (North)      |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Africa (South)      |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Ashkenazi           |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Asia (East)         |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Bangladesh          |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Europe (North East) |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Europe (South East) |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Europe (South West) |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Finland             |   0.0000 |               0.0000 |               0.0025 |
| NA18507   | Ireland             |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Italy               |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Japan               |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Middle East         |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Pakistan            |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Philippines         |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Scandinavia         |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | South America       |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | Sri Lanka           |   0.0000 |               0.0000 |               0.0000 |
| NA18507   | United Kingdom      |   0.0000 |               0.0000 |               0.0000 |

## Timing and memory

| step                                                           | seconds |
|:---------------------------------------------------------------|--------:|
| GRCh38 product staging call (cache hit: hashes of the sources) |    51.0 |
| acceptance staging call (cache hit)                            |    70.8 |

| sample_id | method          | genotype_seconds | cram_seconds | process_peak_rss_mib |
|:----------|:----------------|-----------------:|-------------:|---------------------:|
| HG00188   | allele_fraction |             1.55 |         3.73 |              1731.16 |
| HG00188   | called_genotype |             1.55 |         3.82 |              1731.16 |
| HG00403   | allele_fraction |             1.77 |         4.65 |              1731.16 |
| HG00403   | called_genotype |             1.77 |         5.91 |              1731.16 |
| NA18507   | allele_fraction |             2.23 |         5.22 |              1731.16 |
| NA18507   | called_genotype |             2.23 |         5.40 |              1731.16 |

Scale notes. Dimensions: 1,000 panel sites, 3 samples, 21 groups, 16
PCs, 5,012,826 reference loci. Each call reads the panel-site CRAM once,
then joins at most 1,000 counted sites to the 5.0 million-locus Parquet
reference, so the matched relation and the aggregates hold at most 1,000
x 37 numeric values plus 21 x 16 solver terms per sample. The CRAM calls
take 3.7 to 5.9 s and the genotype calls 1.6 to 2.2 s, each a single
measurement in one render process (`VmHWM` is cumulative for that
process, including package loading and the preceding calls), so this
report gives timing and peak RSS without a scaling or speed verdict;
`STYLE.md` asks for at least three fresh-process repetitions and
1x/2x/4x runs, and the calls are near or below its 5 s floor. It is a
correctness benchmark: it does not vary sites, samples or threads, and
it does not replace the 1x/2x/4x measurements of
`benchmark_ancestry.md`. The one-time staging (product derivation,
remote fetches) is measured outside these tables.

Restrictions. The panel has 1,000 genome-wide sites, not the 17,000 of
the default panel, because cutting each 30x CRAM to the panel sites by
indexed remote access costs one region fetch per site. Sites without a
phase-3 record whose alleles match the source alleles, or with a missing
genotype or depth below 7, drop out of the corresponding side; the
matched counts above are the audit of that loss. Phased 1000 Genomes
genotypes and the 30x CRAMs come from different sequencing datasets of
the same individuals, so agreement includes genotype-call and library
differences, not only the GRCh37-to-GRCh38 derivation.
