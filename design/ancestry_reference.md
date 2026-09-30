Status: current implementation guidance.

# Ancestry reference contract

The Figshare products from the bigsnpr ancestry vignette are the reference
frequencies at https://figshare.com/ndownloader/files/31620968 and the
16-component PC loadings at https://figshare.com/ndownloader/files/31620953.
They are paired products: chromosome, position, and reference/alternate allele
orientation must agree. Reference frequencies are one row per variant and group;
loadings are one row per variant and PC. Select biallelic, unambiguous SNVs with
complete group and PC coverage for a BAM/CRAM panel. Record the genome build,
reference release, and `duckhts_somalier_panel_sha256()` for every materialized
panel. `rduckhts_ancestry_panel()` orders selected loci and assigns zero-based
site indices; its default spacing is 5,000 bases with at most 17,000 sites.

`test/data/ancestry_correction_bigsnpr_1.12.21.tsv` is the 16-component
projection-shrinkage correction vector printed in the vignette at
https://privefl.github.io/bigsnpr/articles/ancestry.html. It is paired only
with those loadings. The duckhtsbench registry checks its committed bytes
before staging. The vector must not be applied to different PCA loadings.

The duckhtsbench registry pins the compressed Figshare files by exact byte
count and SHA-256. The typed, sorted Parquet derivation certifies one locus per
row, group frequencies in [0, 1] and finite PC loadings; its receipt binds those
checks to the source hashes, registry derivation, writer version, and output
hash. A table or view passed to the ancestry wrapper has no receipt identity:
wide references require uniqueness, finite model values and group frequencies
in [0, 1] at matched loci. Long references validate their complete keyed rows.
Automated fetches from some hosts receive HTTP 403 or a WAF
browser challenge. Place browser-downloaded files at their registry cache paths
and run `duckhts_bench_fetch(id)` to verify them without contacting Figshare.
A redirected or merged copy is not an acceptable identity. The epilepsy
input and public GRCh37 phase-3 1000 Genomes chr22 genotypes have independent
registry identities. On chr22, their allele agreement with
the reference products must be checked before site selection. The ancestry
comparison must retain its per-sample matching counts and correlation-gate
failures.

## GRCh38 product

`ancestry_reference_grch38_parquet` is the GRCh38 derivation of the keyed
Parquet product, made by `duckhts_bench_stage_ancestry_grch38_parquet()` with
`duckdb_liftover`, the registered GRCh37-to-GRCh38 chain and the registered
source and destination FASTAs. Its layout is that of the GRCh37 product, so the
panel builder and `rduckhts_ancestry_proportions()` consume it unchanged.

The coded allele of a locus is the source `allele_b`: the frequencies and
loadings describe it. Liftover reports its destination spelling through
`reverse_complemented`, and the derived key is `allele_a` = the other allele,
`allele_b` = the coded allele. Every frequency and loading is copied unchanged,
which is the wrapper's reversal convention (a reversed input meets an
unchanged reference row). When the coded allele is the destination reference,
liftover reports a swap; the locus is counted but not rewritten, because
`1 - f` with negated loadings shifts each PC by `-sum(U)`, which the per-PC
correction coefficients scale on one side only. Loci that do not map, are not
single-base substitutions, land off chromosomes 1-22, keep neither destination
allele, or share a destination position with another locus are dropped, and
each is counted.

The receipt (`reference.grch38.parquet.sources.tsv`, columns `field` and
`value`) binds the SHA-256 of the source Parquet, chain, source FASTA and
destination FASTA, the registry derivation, the DuckDB, DuckHTS htslib and
Rduckhts versions and the output SHA-256 to the denominators: input loci,
mapped, rejected by reason, reverse-complemented, swapped,
duplicate-destination dropped and output loci. Every input locus is in exactly
one of rejected, duplicate-destination dropped or output. A companion
`reference.grch38.liftover.parquet` maps each output locus to its GRCh37 source
locus and alleles; its hash is in the receipt. A cached product whose receipt
does not match its sources is an error, not a rebuild.

`benchmarks/benchmark_ancestry_grch38.md` is the acceptance comparison:
ancestry proportions from the registered 30x GRCh38 CRAMs against the GRCh38
product, and from the GRCh37 phase-3 phased genotypes of the same individuals
against the GRCh37 product, at the same genome-wide panel loci. The genotypes
are read by indexed `read_bcf(region := ...)` at the panel loci and the CRAMs are
cut to the panel loci; `duckhts_bench_stage_ancestry_grch38_acceptance()`
stages both with a receipt. The tolerance is declared in the report source; the
per-sample matching counts and correlation-gate results are retained.
