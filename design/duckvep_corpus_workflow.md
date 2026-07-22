# DuckVEP corpus derivation and release-conformance workflow

Status: current release-conformance and corpus-derivation operating contract.

This note is the mental map for rebuilding DuckVEP evidence on another machine and for
advancing the compatibility target from Ensembl VEP 116 to a later release. It describes
which facts must remain controlled, which repository files execute the work, and where the
current corpus matrix is incomplete. It is not a second command catalogue: the scripts,
manifests, and rendered reports linked below remain the executable authorities.

## The four things a conformance run must not conflate

Every result has four independent inputs:

1. **The biological source** supplies variants or generated event geometry. Examples are a
   dated ClinVar release, GIAB HG002, an HPRC graph callset, or a seed-defined structural
   generator.
2. **The annotation model** supplies transcripts, exons, sequences, translation metadata,
   mature-miRNA segments, and optionally regulatory and motif features. It is identified by
   species, assembly, Ensembl or Ensembl Genomes release, and transcript-selection policy.
3. **The oracle** is the exact Ensembl VEP executable/cache combination, a pinned VEP plugin,
   or an official Ensembl release VCF field. The source corpus is not the oracle.
4. **The comparison evidence** is the complete pair relation and its aggregation. A pair is
   keyed by input allele or event plus the transcript, RegulatoryFeature, or MotifFeature
   that VEP annotates. Missing, extra, unresolved, and discordant pairs stay in the
   denominator.

The data flow is therefore:

```text
declared source
  -> immutable local source object
  -> deterministic eligibility, normalization, and sampling
  -> one exact VEP input relation
  -> VEP oracle pairs + DuckVEP pairs
  -> duplicate-key check and full outer comparison
  -> append-only aggregate history + retained counterexamples
```

A successful row in a report is meaningful only when those four identities are visible.
Conversely, reproducing the same input VCF while silently changing the Ensembl model or VEP
cache is a different experiment.

## Artifact identity without repeated full-file hashing

Cryptographic digests are useful when bytes cross a trust or publication interface. They
are not the biological method and should not consume every run.

A digest proves only that two byte streams are identical. It does not prove that the
assembly, cache release, transcript filter, contig aliases, ALT ordinals, VEP options, or
comparison denominator are correct. Recomputing a digest over a large reference can also
consume the same storage bandwidth needed by the annotation run and evict useful pages from
the filesystem cache. Treating a digest as a universal proxy for provenance is therefore
both incomplete and unnecessarily expensive.

Use this policy:

- **Acquisition:** verify an upstream checksum, immutable object version, or release
  manifest once. Record the resulting local artifact identity in its staging receipt.
- **Ordinary reuse:** use the acquisition receipt from an immutable/content-addressed local
  object, read-only filesystem snapshot, or equivalently enforced store. Path, byte size,
  and modification time are cheap cache-validation metadata, not identity. Do not reread a
  multi-gigabyte FASTA, VCF, cache, or DuckDB model merely to calculate the same digest
  again.
- **Publication or transfer:** verify the bytes once before publishing and once after
  downloading on another machine. An explicit release-audit mode may deliberately pay this
  I/O cost.
- **Execution outputs:** validate row/allele/object cardinalities and comparison totals on
  every run. These denominators detect truncated or selectively filtered analyses more
  directly than repeatedly hashing an unchanged input.

Small checked-in authorities and source files may be byte-verified on every audit because
that cost is negligible. Large-artifact digests belong in one acquisition or publication
manifest and are then referenced, not recomputed as ceremony. The shared evidence helper
caches a digest after its first complete read and reuses it only while the canonical path,
byte size, modification time, and change time are unchanged. Those metadata validate the
cache entry; they do not replace the digest or create artifact identity. Set
`DUCKVEP_ARTIFACT_VERIFY=full` for acquisition, transfer, publication, or an explicit
release audit. `DUCKVEP_EVIDENCE_DIGEST_CACHE` may place the small local cache somewhere
other than its XDG/user-cache default.

## Reproducibility levels

The project deliberately separates four levels instead of claiming that every historical
large-corpus run is already one-command portable.

| Level | What another machine needs | Current state |
| --- | --- | --- |
| In-repository | Git checkout, C/R build dependencies | Fixed SQL/R/C witnesses, randomized property generators, finite-state manifests, and the exact public Ensembl test mirror run offline. |
| Oracle environment | Pinned VEP environment plus the small checked-in fixtures | VEP 116 self-tests and the generated GFF differential are reproducible from the documented environment. |
| External corpus | Declared public source plus its deterministic staging recipe | HPRC, Sniffles2/1000G, and dbVar chromosome-22 recipes are executable; some earlier GIAB, ClinVar, GRCh37, and *P. falciparum* campaigns retain result lineage but still need their acquisition commands consolidated here. |
| Published release pack | Downloadable models, selected source corpora or shards, receipts, and expected pair summaries | Planned for Zenodo or equivalent versioned object storage; not yet published. |

The current strict state-machine gate proves that an executable campaign passed and that
the implementation/harness inputs have not changed since it ran. It does not prove that an
external corpus or model is downloadable on another machine. Cross-machine portability is
a separate release-pack gate that cannot become executable until the external pack is
published. A historical row remains useful evidence, but it is not portable merely because
it contains a digest.

## Corpus families and what each one searches

The checked campaign registry is
[`state_machine_campaigns.tsv`](../test/duckvep/conformance/data/state_machine_campaigns.tsv).
It maps each required state-machine campaign to its evidence ledger and oracle authority.
The larger history and exact denominators are rendered from the CSV data by
[`duckvep_conformance.Rmd`](../benchmarks/duckvep_conformance.Rmd); do not copy those numbers
into this operating contract.

| Corpus family | Derivation and purpose | Executable authority |
| --- | --- | --- |
| Fixed witnesses | Hand-minimized counterexamples for exon/intron, CDS, start/stop, miRNA, exceptional edits, HGVS shifts, structural, and BND states. Every newly found rare state becomes a fixed witness. | `generate_witnesses.R`, SQL/R/C tests, executable VEP differential. |
| Statistical state exploration | A seed-defined distribution over SNV, MNV, insertion, deletion, delins, splice, CDS, terminal-codon, and rare configuration strata. It searches for counterexamples; it is not a population model. | `property_history.R`, `corpus_differential.R`, coverage requirements. |
| GIAB GRCh38 | Public HG002 benchmark variants. Used for ordinary WGS topology, dense transcript expansion, regulation/motif overlap, end-to-end I/O, and comparison with FastVEP. | Indexed-cache VEP 116 and DuckVEP on the identical staged alleles. |
| ClinVar GRCh38 | A dated public release. Coding-enriched and cross-chromosome selections exercise clinical allele shapes; complete dense shards exercise HGVS and NMD rather than rewarding a sparse WGS average. | Indexed-cache VEP 116; NMD additionally uses pinned VEP Plugins release/116. |
| Human GRCh37 | Variants are evaluated against the separate Ensembl 116 GRCh37 model and matching reference/cache. It is not a coordinate flag on the GRCh38 model and carries no invented native MANE state. | VEP 116 GRCh37 cache and the prepared GRCh37 core model. |
| *P. falciparum* | Non-human, compact-genome evidence with Ensembl Genomes release 63, assembly `GCA000002765v3`, and codon tables 1, 4, and 11. | VEP 116 plus its species/cache release 63 data and the matching prepared model. |
| Generated exact SV | Real transcript geometry is crossed with exact DEL, DUP, tandem DUP, INV, INS, and CNV relations around transcript, exon, intron, UTR, CDS, and start/stop sites. | `corpus_differential.R --event-mode structural`. |
| Paired BND | Same- and cross-contig endpoint pairs, four bracket orientations, and transcript/regulatory admission distances. VEP runs with buffer size one to prevent neighboring records changing its coordinate-only candidate tree. | `corpus_differential.R --event-mode breakend`. |
| Real SV callsets | HPRC carried long alleles, Sniffles2 joint 1000G calls, and dbVar regions check producer encodings and large-event distributions that a generator may miss. | `scripts/stage_duckvep_conformance_corpora.sh`; executable differential after staging. |
| Official Ensembl release VCF | Ensembl's `VE` index maps published consequence terms to the original GVF allele before VCF padding. Small release-specific shards provide a fast, Perl-free CI oracle for states present in the release. | `release_vcf_differential.R`; current first fail-closed stratum is literal SNV. |
| HPRC/pangenome long alleles | Carried ALT alleles from four African-ancestry HPRC v2 samples on chromosome 22. This searches long literal and pangenome-derived states; it does not validate graph-coordinate annotation. | Public versioned HPRC object plus the staged sites-only VCF. |

Do not replace the dense ClinVar/HGVS, exact SV/BND, non-human, or GRCh37 lanes with a
single WGS average. They cover different reachable parts of the VEP state machine.

## Deriving a corpus

### 1. Declare the source before downloading it

Record the source organization, release/date, assembly, species, URL or immutable object
version, access terms, and intended role. Prefer sources with their own dated release
manifest. A short-lived signed URL is not a durable release authority unless its underlying
object version is also known.

Network access belongs only in an explicit staging script. The extension build and ordinary
test suite remain network-free.

### 2. Preserve the source object

Keep the downloaded object outside Git and do not rewrite it in place. Store its acquisition
receipt beside it. A derived subset points back to that source receipt and records the
selection command; it does not pretend to be the upstream release.

The existing external structural/pangenome staging entry point is:

```sh
scripts/stage_duckvep_conformance_corpora.sh /data/duckvep-corpora all
```

It currently stages:

- four African-ancestry HPRC v2 samples on chromosome 22, retaining carried ALT alleles,
  splitting multiallelic records, removing genotypes, and mapping `chr22` to `22`
  (`HG02055`, `HG02145`, `HG02723`, and `HG03098`);
- the Sniffles2 2.0.7 joint 1000 Genomes chromosome-22 SV sites; and
- the dated dbVar GRCh38 chromosome-22 variant-region slice.

The script contains the exact public objects and deterministic `bcftools` pipeline. Its
small derived outputs may be byte-verified cheaply; the multi-gigabyte remote sources are
identified by their provider manifest, object version, or ETag rather than downloaded only
to hash them.

Two important small-variant sources have stable public locations even though their download
modes are not yet folded into that staging script:

- GIAB HG002 v4.2.1 GRCh38 is
  [`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz)
  with its adjacent tabix index; and
- the measured ClinVar source is the dated NCBI archive object
  [`clinvar_20260706.vcf.gz`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260706.vcf.gz),
  with the adjacent NCBI MD5 file and tabix index.

Use the dated ClinVar archive, never the moving `clinvar.vcf.gz`, for a retained campaign.
The current benchmark receipts record the exact local identities and denominators for these
two objects; adding `giab-hg002-grch38` and `clinvar-grch38-DATE` modes to the staging script
will consolidate the remaining mechanical download step.

### 3. Normalize representation without losing identity

Keep two relations:

```text
source_record(source_id, source_record_id, original_contig, POS, REF, ALT[], INFO...)
analysis_allele(source_id, source_record_id, alt_index,
                seq_region, normalized_pos, normalized_ref, normalized_alt)
```

Split or minimize alleles only in `analysis_allele`. Preserve the original record, ALT
ordinal, symbolic type, END/mate geometry, confidence coordinates, genotype arrays, and raw
ALT where they exist. A biallelic analysis row remains traceable to its multiallelic source.
Do not remap genotype likelihood arrays merely to make annotation convenient.

Contig aliases are a declared relation. For the human primary assemblies, accept the
deliberate `chr` prefix mapping and the mitochondrial aliases required by the selected
reference. Do not guess aliases from failed joins.

### 4. Derive deterministic eligibility and sampling

Record every exclusion category before sampling. The denominator is:

```text
source records
  -> decomposed source alleles
  -> eligible analysis alleles
  -> selected analysis alleles
  -> emitted annotation-object pairs
```

The corpus runner samples within allele shape and length-change bins using a declared seed.
Use `--sample-per-shape 0` for a complete eligible corpus. A rare-state campaign uses
dedicated strata and explicit minimum counts; it does not rely on a rare draw from a broad
distribution. Failures and unresolved rows remain in the pair artifact even when the run
fails.

For annotation-density performance work, use
[`stage_duckvep_dense_corpus.R`](../scripts/stage_duckvep_dense_corpus.R). It ranks fixed
coordinate tiles from the real model by transcript, coding exon, non-coding RNA,
regulatory/motif, and observed-variant density. `--all-tiles` retains the complete primary
assembly. The resulting database is a workload description, not a conformance oracle.

### 5. Build the matching model

Follow [`duckvep.md`](duckvep.md) for source tables, transcript filtering, FASTA matching,
sequence reconstruction, mature-miRNA, peptide edits, regulation/motif preparation, and
model loading. A corpus may run only when species, assembly, contig coordinate system, and
release match the model and VEP cache.

The current model families are:

- Ensembl 116 human GRCh38, optionally with release-matched funcgen features;
- Ensembl 116 archived human GRCh37; and
- Ensembl Genomes release 63 *P. falciparum* `GCA000002765v3`.

The rendered cardinalities and logical model receipts live in
[`duckvep_model_receipts.csv`](../benchmarks/data/duckvep_model_receipts.csv). A logical
model receipt is preferable to repeatedly hashing the whole DuckDB file: it identifies the
prepared biological relations while permitting harmless physical rewrites of the database.

The VEP oracle cache is a separate upstream artifact from the DuckVEP model. Install it in
an explicit networked staging step from the exact VEP checkout, for example:

```sh
perl .sync/ensembl-vep/INSTALL.pl \
  --AUTO cf --CACHE_VERSION 116 \
  --CACHEDIR /data/vep-cache \
  --SPECIES homo_sapiens --ASSEMBLY GRCh38 \
  --USE_HTTPS_PROTO
```

Repeat with `--ASSEMBLY GRCh37` for the separate human archive cache. The existing
non-human campaign uses the VEP-compatible *P. falciparum* cache directory
`plasmodium_falciparum/63_GCA000002765v3`; retain cache release 63, assembly
`GCA000002765v3`, and the VEP 116 executable as three separate fields. Preserve each
cache's `info.txt`: its source genebuild, assembly, variation release, and regulatory flag
are more informative for biological identity than hashing every chromosome cache shard on
each run.

### 6. Run the oracle and DuckVEP on the same allele relation

The general entry point is:

```sh
make -f Makefile duckvep-corpus-differential DUCKVEP_DIFFERENTIAL_ARGS="\
  --corpus NAME \
  --vcf /data/corpus.vcf.gz \
  --database /data/model.duckdb \
  --model-sql '' \
  --cache-dir /data/vep-cache \
  --fasta /data/reference.fa \
  --assembly ASSEMBLY --species SPECIES \
  --sample-per-shape 0"
```

Use the GFF oracle only for small fixtures and importer diagnosis. Release evidence uses the
matching indexed VEP cache because the core database and GFF importer do not necessarily
admit the same objects. Structural and paired-breakend command examples, NMD plugin usage,
buffer-size policy, spill controls, and output relations are maintained in the
[`conformance README`](../test/duckvep/conformance/README.md).

The comparison projects both engines to unique `(source allele/event, object)` rows and
then full-outer-joins them. Duplicate keys fail before the join. Consequence terms, status,
reason, NMD state, HGVSc, and HGVSp are compared as distinct declared metrics; an absent
value is not silently treated as agreement.

### 7. Publish evidence, not scratch state

Keep the complete pair artifact for failed or novel runs until the mismatch is minimized.
Append aggregate results only after the run has passed its declared gates. The checked
histories are:

- `conformance_history.csv` for consequence, impact, status, and NMD strata;
- `hgvs_history.csv` for HGVSc/HGVSp;
- `property_history.csv` and `property_coverage_history.csv` for randomized C campaigns;
- `state_machine_campaigns.tsv` for the required release matrix.

Rendered reports read these data files. Never hand-edit a rendered `.md` or recollect a
number in prose when the source data already contain it.

## Official Ensembl release VCF/GVF lane

The release VCF is derived from Ensembl Variation GVF. `VE` contains
`Consequence|Index|Feature_type|Feature_id`; `Index` is the stable relation between a
published consequence and the original GVF `Variant_seq`. VCF conversion may pad or
otherwise transform alleles, so indels cannot be joined by comparing rendered CSQ allele
text with VCF ALT text.

The current release-116 chromosome-22 source and Parquet staging measurements are recorded
in [`duckvep_release_parquet.csv`](../benchmarks/data/duckvep_release_parquet.csv). Stage it
with `bench-duckvep-release-parquet`, then run `test-duckvep-release-vcf` as documented in
the conformance README. The current fail-closed implementation covers literal SNVs. Before
claiming the release VCF as a general oracle, implement the producer-faithful GVF
`Variant_seq`/`Index` mapping for indels and other transformed alleles.

Small model-plus-release-VCF shards are the preferred ordinary-CI artifact. They are fast
and do not start Perl. They complement, rather than replace, executable VEP on generated
novel states.

## Moving from VEP 116 to the next release

Treat a new VEP release as a new authority, not an in-place update of old evidence.

1. Pin the new VEP, Ensembl Variation, VEP Plugins, Ensembl/Ensembl Genomes data, cache,
   reference, and BioPerl environment. Preserve the VEP 116 authorities alongside them.
2. Generalize the currently VEP-116-specific upstream manifest/checker so authorities are
   keyed by repository, ref, and commit and one self-test receipt is required per supported
   release. Then mirror the new public test sources in a release-labelled directory and
   record a verbose self-test transcript. Passing upstream tests establishes oracle health,
   not DuckVEP conformance.
3. Re-extract the SO registry/rule order and inspect the semantic diff. Regenerate code only
   from the new pinned authority; never edit generated tables by hand.
4. Build new models for every supported species/assembly. Compare logical relation counts,
   transcript selection, sequence availability, codon tables, peptide edits, mature-miRNA,
   and regulation/motif inputs with the previous release.
5. Run fixed witnesses and C property/state-machine checks first. A changed VEP behavior
   gets a release-labelled erratum and a minimal fixed witness; do not overwrite the VEP
   116 explanation.
6. Reuse a source corpus only when its assembly and normalization relation are unchanged.
   Re-derive official Ensembl release VCF shards because they are release products.
7. Run the complete matrix: generated rare states, GIAB, complete/dense ClinVar HGVS and
   NMD, GRCh37, non-human, regulatory/motif, exact SV, paired BND, and long HPRC alleles.
8. Run performance on the same checked workload and separately measure any changed model or
   output surface. Do not compare unlike corpora or omit HGVS/regulatory rows to rescue a
   headline.
9. Add new campaign rows only after the state transition and outcome relations name the
   newly supported semantics. Preserve every discovered counterexample.
10. Publish the new model/shard pack, then verify it once from a clean second-machine
    download before calling the release portable.

## External release-pack contract

One versioned Zenodo record or equivalent object-store prefix should contain:

- prepared DuckDB models for each release/species/assembly;
- small official release-VCF/GVF oracle shards for ordinary CI;
- public derived corpus shards whose upstream terms permit redistribution;
- a machine-readable manifest with source release/object, derivation script revision and
  arguments, assembly/species, normalization and contig-alias policy, seed, cardinalities,
  logical model identity, byte sizes, and acquisition/publication digests; and
- expected pair-level aggregate summaries and tool/environment versions.

The manifest should contain one digest per published object, computed during publication.
Consumers verify once after download and then use the local receipt. Do not publish licensed
or controlled clinical callsets; publish their derivation schema and use a synthetic/public
equivalent for the portable gate.

## Present gaps before the workflow is fully portable

- Add explicit GIAB and dated-ClinVar modes to the external-corpus staging script, then
  consolidate the historical GRCh37 and *P. falciparum* corpus selection commands into
  equivalent entry points.
- Publish the Ensembl models and small release-VCF oracle shards outside Git.
- Publish a cross-machine release-pack checker; the current state-machine audit establishes
  semantic/current-source evidence, not remote artifact availability.
- Complete producer-faithful non-SNV mapping for the official release VCF/GVF lane.
- Add release-labelled recipes for full ClinVar HGVS, broader species, HPRC long alleles,
  exact SV/BND, and eventually phased Haplosaurus edit sets.

Until those items are closed, the repository is sufficient to understand and rerun the
algorithms and small authorities, while the historical large-corpus results remain
well-described but not yet a self-contained portable release pack.
