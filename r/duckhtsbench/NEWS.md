# duckhtsbench 0.0.0.9000

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
