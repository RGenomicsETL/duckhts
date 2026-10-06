# duckhtsbench 0.0.0.9001

- register the `bam-mismatch` workload for
  `benchmark_bam_mismatch_counts.Rmd`: `bam_mismatch_hg00403_chr1_cram`, the
  chr1 alignments of the public 30× CRAM `ancestry_30x_hg00403` as a CRAM
  slice; `bam_mismatch_panel_chr1_source`, the chr1 file of the phased
  1000 Genomes panel; `bam_mismatch_chr1_mask_bcf`, its records without
  genotypes as an indexed BCF; and `bam_mismatch_chr1_mask_half_bcf` and
  `bam_mismatch_chr1_mask_quarter_bcf`, the same records thinned to the
  positions divisible by 2 and by 4 (`modulo` in the identity), for scaling
  the mask records alone. `benchmarks/bam_mismatch_stage.R` stages the derived
  files with the samtools and bcftools of RBCFTools, reading each source from
  its cached copy or, without one, from its registered URL. A slice is pinned
  by its alignment count, position sum and stored-base sum, and a mask by its
  record count and position sum; staging publishes nothing when they differ.
  `make test-bam-mismatch-staging` is the network-free test.

- register `roh_counts_validation_chr20_genotypes`: FORMAT/GT genotypes of
  NA18507, HG00403 and HG00188 at the chr20 ROH evaluation sites, derived from
  `roh_ancestry_chr20_source` and `roh_ancestry_chr20_af_sites` by
  `benchmarks/roh_counts_validation/stage_genotypes.R`, with a record and sample
  identity and a network-free staging test. `benchmark_roh_counts_validation.Rmd`
  reads it with the three `roh_counts_chr20_*` artifacts and checks all four
  identities before decoding.

- register the `genbank-plasmid` workload: the first four
  `plasmid.N.genomic.gbff.gz` parts of NCBI RefSeq release 237, each pinned by
  NCBI's published MD5, and three derived record files joining 1, 2 and 4 of
  them through a new `gunzip_concatenate` transform
  (`duckhts_bench_stage_gunzip_concatenate()`), staged by
  `duckhts_bench_stage_genbank_plasmid()` with a network-free staging test; it
  is the record-count scaling input `benchmark_genbank_named_attributes.Rmd`
  reads. NCBI serves only the current release's sequence files, so the parts
  are registered as `public_current_release` and staging stops before
  downloading, naming the archived release catalog, once NCBI has moved past
  release 237. Gunzip derivations close every connection when a source fails

- register `ancestry_reference_grch38_parquet`, the keyed ancestry reference
  lifted from GRCh37 with `duckdb_liftover` (registered chain and FASTAs), and
  stage it with `duckhts_bench_stage_ancestry_grch38_parquet()`. Its receipt binds
  the source, chain, FASTA and output SHA-256 values and the tool versions to the
  input, mapped, rejected-by-reason, swapped, duplicate-destination and output
  locus counts; a companion map retains each output locus's GRCh37 source.

- stage the GRCh38 acceptance inputs with
  `duckhts_bench_stage_ancestry_grch38_acceptance()`: GRCh37 phase-3 genotypes at
  the panel loci by indexed `read_bcf(region := ...)` and panel-site CRAMs cut
  from the registered 30x CRAMs, with a receipt of remote identities and output
  hashes

- register checksum-verified MANE and GENCODE GFF3 inputs and a pinned GFFBase
  wheel for the feature-database benchmark

- stage the HPRC genotype-reader source and index under genotype-owned registry
  entries, retaining the pinned S3 version IDs and index checksum

- licensed GPL (>= 2), like the other DuckHTS R packages (previously MIT)

- require matching reference/read SHA-256 identities before reusing a staged
  ONT BAM; receipts without identities and changed inputs require derivation, while tool
  version changes alone do not

- record the staged ONT BAM's SHA-256 and byte size in its receipt so in-place
  BAM changes and receipts from the previous format require derivation

# duckhtsbench 0.0.0.9000

- identify both aligned-block and CIGAR-validation benchmarks as consumers of
  the registered ONT inputs

- resolve samtools for ONT staging from PATH or the optional RBCFTools package,
  and test missing executables before checking uncached sources

- register the `ont-ecoli-k12` workload: the NCBI RefSeq E. coli K-12 MG1655
  assembly FASTA pinned by NCBI's published MD5 with its uncompressed form as
  a derived artifact, ENA run `ERR14686255` (25,950 MinION reads, PRJEB86481)
  pinned by ENA's published MD5 and byte size, and the coordinate-sorted BAM
  derived from them with `minimap2 -x map-ont` and `samtools`, staged by
  `duckhts_bench_stage_ont_ecoli()` with a network-free staging test. It is
  the long-read input `benchmark_cigar_aligned_blocks.Rmd` reads.

- document the complete exported registry and staging API and resolve utility
  functions through their owning namespaces, keeping source-package checks clean

- register the network-free `genbank-memory-scaling` fixture workload used to
  separate cumulative record count from largest-record growth in the rendered
  GenBank memory report

- register the `genbank-reader` workload: the NCBI RefSeq E. coli K-12 MG1655
  assembly archive pinned by NCBI's published MD5, SHA-256 and byte size, with
  its uncompressed `.gbff` as a derived artifact, staged by
  `duckhts_bench_stage_genbank()`; rendering resolves the staged path without
  network access

- stage an all-sample HPRC regional genotype cohort as matching VCF.gz/BCF
  representations, with source-version provenance and network-free staging tests

- pin and stage committed small BCF/VCF scan benchmark fixtures and indexes
  through the registry, with byte/hash validation and network-free staging tests

- register and stage the deterministic multi-region VCF, BCF and indexes without
  network access, with a small reconstruction and exact-position test

- introduce the internal benchmark artifact registry and R-native portability check;
  registered VariantKey provider raw sources and derivations use that authority;
  cache reuse validates
  declared publisher byte, MD5, or Ensembl `sum` identities before provenance is written; the
  registered Riker BAM again supports the whole-genome mosdepth benchmark entry point
