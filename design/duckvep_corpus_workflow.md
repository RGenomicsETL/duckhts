# DuckVEP corpus and release-conformance policy

Status: current evidence-integrity policy. Scripts, manifests, test fixtures, and rendered
reports are the executable authorities for individual campaigns.

This note controls what a DuckVEP conformance result means. It does not prescribe staging
commands, enumerate transient source locations, or duplicate campaign results.

## Keep four identities separate

Every result names four independent inputs:

1. **Biological source** — the released callset or generated event distribution.
2. **Annotation model** — species, assembly, Ensembl or Ensembl Genomes release,
   transcript-selection policy, reference, and optional core feature model.
3. **Oracle** — the pinned executable VEP, cache, and plugin environment. An Ensembl release
   VCF is a separate release-pipeline product, not a substitute for executable VEP.
4. **Comparison evidence** — the complete pair relation and its aggregation.

A source VCF, DuckVEP model, and VEP cache that share an assembly name are not thereby the
same experiment. Replacing any one of these identities requires a new receipt and result.

## Preserve the full comparison denominator

The evidence path is:

```text
declared source
  -> immutable local source object
  -> deterministic eligibility, normalization, and sampling
  -> one exact VEP input relation
  -> VEP and DuckVEP pair relations
  -> duplicate-key check and full outer comparison
  -> aggregate history and retained counterexamples
```

A comparison key is one input allele or event plus one annotated transcript,
RegulatoryFeature, or MotifFeature. Missing, extra, unresolved, consequence, NMD, HGVSc,
HGVSn, and HGVSp outcomes remain separate metrics. No release mode may discard a discordant
row or treat absence as agreement.

Record each exclusion category before sampling:

```text
source records
  -> decomposed source alleles
  -> eligible analysis alleles
  -> selected analysis alleles
  -> emitted annotation-object pairs
```

Sampling is declared by seed and stratum. A complete eligible corpus uses no sampling; a
rare-state campaign specifies the state minimum rather than depending on an incidental broad
draw. Failed and unresolved rows remain in the pair artifact.

## Retain source identity without ritual rehashing

Verify an upstream checksum, immutable object version, or release manifest when acquiring an
artifact. Preserve that acquisition receipt beside an immutable or read-only local object.
For ordinary reuse, validate cheap path/inventory metadata rather than rehashing a large
unchanged artifact. Verify bytes again only at an explicit transfer, publication, or
release-audit boundary.

Every execution still validates input, object-pair, and result cardinalities. Those
semantic denominators reveal truncation or selective filtering that a repeated byte digest
cannot. Logical model receipts identify prepared biological relations without treating a
physical DuckDB-file rewrite as a model change.

## Separate source records from analysis alleles

A corpus retains both relations:

```text
source_record(source_id, source_record_id, original_contig, POS, REF, ALT[], INFO...)
analysis_allele(source_id, source_record_id, alt_index,
                seq_region, normalized_pos, normalized_ref, normalized_alt)
```

Splitting or minimizing occurs only in `analysis_allele`. The original record, ALT ordinal,
symbolic type, event geometry, confidence coordinates, genotype arrays, and raw ALT remain
traceable. Contig aliases are an explicit relation; failed joins do not license guessed
aliases. Annotation convenience does not license genotype-likelihood remapping.

Network access belongs only in explicit staging scripts. The extension build and ordinary
test suite remain network-free. The
[conformance README](../test/duckvep/conformance/README.md) owns commands, external source
locators, and campaign-specific execution details.

## Use complementary evidence lanes

Fixed witnesses preserve minimized counterexamples. Seeded pure-C properties search cheap
mechanical state spaces with coverage counters. Executable VEP differentials are the
biological authority for the declared scope. Real corpora expose representation and
high-fan-out states that generated cases may miss; generated structural and breakend cases
cross geometry states that an incidental callset may not cover.

A counterexample is minimized, classified, and made a fixed witness or required generated
stratum. Re-running a broad distribution and hoping to rediscover it is not regression
coverage. The checked campaign registry, pair artifacts, and rendered reports own the
current campaign matrix and measured denominators.

## Advance a VEP release as a new authority

Do not overwrite VEP-116 evidence in place. For a later release:

1. pin the VEP, plugin, source, cache, reference, and BioPerl authorities separately;
2. verify the new upstream environment before treating it as an oracle;
3. regenerate derived metadata only from the new pinned source;
4. build a new matching DuckVEP model and compare source/model facts explicitly;
5. run fixed witnesses and property gates before the complete executable matrix;
6. record changed VEP behaviour as a release-labelled erratum with a minimized witness; and
7. measure the same declared workload separately from any changed model or output surface.

A release transition never rescues a result by changing selection, dropping HGVS or
regulatory rows, or comparing unlike corpora.

## Publish a portable release pack deliberately

Cross-machine reproducibility requires a versioned pack, not a historical result with a
digest. A pack needs receipted models, permitted corpus shards, a machine-readable manifest
of source/derivation/assembly/normalization/seed/cardinality identities, and expected
pair-level summaries. Publish only redistributable data; controlled inputs require a
reproducible derivation schema and a portable public or synthetic gate.

Until such a pack is published and verified from a clean download, a campaign may establish
semantic evidence without claiming one-command portability on another machine.
