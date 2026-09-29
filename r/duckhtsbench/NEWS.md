# duckhtsbench 0.0.0.9001

- register checksum-verified MANE and GENCODE GFF3 inputs and a pinned GFFBase
  wheel for the feature-database benchmark

- stage the HPRC genotype-reader source and index under genotype-owned registry
  entries, retaining the pinned S3 version IDs and index checksum

- licensed GPL (>= 2), like the other DuckHTS R packages (previously MIT)

- require matching reference/read SHA-256 identities before reusing a staged
  ONT BAM; receipts without identities and changed inputs require derivation, while tool
  version changes alone do not

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
