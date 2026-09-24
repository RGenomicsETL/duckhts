# duckhtsbench 0.0.0.9001

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

- repin the Ensembl-116 GRCh38 model after removal of the redundant short-tail
  column and use its logical hash in the registered cache path, allowing fresh
  compilation without reusing or overwriting the previous model artifact

- verify streamed VEP archives against a pinned SHA-256 or MD5, require the
  digest subprocess to succeed, and reject HTTP-metadata-only cache identities;
  shared local-artifact SHA-256/BSD-sum verification also rejects failed checksum
  commands that print matching values

- expose the canonical nested model tables through shared read-only flat SQL
  projections for conformance and benchmark consumers, preserving every exon,
  mature-miRNA segment, peptide edit, transcript label and sequence byte

- stage complete VEP-116 GRCh38 chromosome-21 cache shards and all root metadata
  by streaming the pinned upstream archive, without storing the full compressed
  archive; verify transfer identity and publish a canonical acquisition receipt

- share network-free committed-fixture staging between BCF scan and DuckVEP
  projection workloads, retaining registry byte identities and provenance

- stage an all-sample HPRC regional genotype cohort as matching VCF.gz/BCF
  representations, with source-version provenance and network-free staging tests

- pin and stage committed small BCF/VCF scan benchmark fixtures and indexes
  through the registry, with byte/hash validation and network-free staging tests

- register and stage the deterministic multi-region VCF, BCF and indexes without
  network access, with a small reconstruction and exact-position test

- register and stage the HPRC v2 African-four, Sniffles2 1KGP, and dbVar GRCh38
  chr22 DuckVEP corpora with pinned raw-source identities, source indexes and
  manifests, deterministic derived VCF/index checksums, atomic publication,
  adjacent provenance, and a network-free fake-source CLI test

- add the clean-cache `duckvep_ensembl116_model` producer from checksum-pinned
  public Ensembl 116 core/funcgen MySQL dumps and the matching GRCh38 FASTA;
  validate source metadata and deterministic receipts before provider exports

- introduce the internal benchmark artifact registry and R-native portability check;
  registered VariantKey provider raw sources, derivations, model exports, and
  their network-free staging closure now use that authority; cache reuse validates
  declared publisher byte, MD5, or Ensembl `sum` identities before provenance is written; the
  registered Riker BAM again supports the whole-genome mosdepth benchmark entry point
