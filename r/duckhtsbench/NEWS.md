# duckhtsbench 0.0.0.9000

- introduce the internal benchmark artifact registry and R-native portability check;
  registered VariantKey provider raw sources, derivations, model exports, and
  their network-free staging closure now use that authority; cache reuse validates
  declared publisher byte, MD5, or Ensembl `sum` identities before provenance is written; the
  registered Riker BAM again supports the whole-genome mosdepth benchmark entry point
