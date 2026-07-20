# DuckVEP architecture

Status: current implementation contract. Code and tests are authoritative for implemented
behavior; the remaining gaps are linked by full URL in the relevant sections below.

## Purpose and compatibility

DuckVEP is the deterministic variant-consequence kernel inside DuckHTS. It targets Ensembl
VEP 116 (`57ea5c52340acc1f156267f810ad162e26597082`) with Ensembl core
`c0cf13daa961d80584bad797b2eb0ff3a7500ef3` and variation
`2fb834b987ede3824e200197a838ce11e91aeb4b`. The behavioral target is not
human-only: the kernel accepts every NCBI codon table supported by VEP 116's BioPerl
translator. Full-model acceptance currently covers human GRCh37 and GRCh38; non-human
acceptance covers *Plasmodium falciparum* from Ensembl Genomes 63 paired with VEP 116,
including its codon tables 1, 4, and 11. Evidence for one species or assembly is never
treated as evidence for another.

The kernel emits structured transcript consequences. VEP-compatible CSQ text and HGVS are
edge projections, not the internal representation. Prepared-event and CDS-edit facts are
the shared consequence/HGVS authorities; an HGVS-facing transcript-edit carrier adds
VEP's clipped transcript-slice state without changing those semantic edits. Allocation-free
`c.`/`n.`/`p.` fact/render kernels, a worker-owned indexed-FASTA provider, and cumulative
`duckvep_annotate_hgvs(...)` SQL adapter are implemented for independent literal small
variants. Strict executable-VEP evidence covers the fixed position-one/right-anchor cases
and 56,998 chromosome-21 ClinVar transcript pairs with zero HGVSc/HGVSp differences.
Genomic HGVS, RefSeq RNA-edit shifting, structural/BND HGVS, broader cross-species HGVS
distributions, and compound/phased HGVS remain open.
VEP 116 is the behavioral authority; pure-C properties and bcftools csq supply independent
checks for mechanics and phased edit state.

## Mental model

DuckVEP does two different jobs at two different times.

First, it compiles Ensembl core and funcgen relations plus matching reference sequence into
an immutable execution model. DuckDB keeps the rich source relations, stable identifiers,
attributes, and hashes. The compiler checks them, assembles transcript-oriented sequence,
and publishes compact C arrays containing only facts used repeatedly during annotation.

Second, it runs a sorted stream of variant alleles through that model. Each allele is
interpreted once, matched to candidate transcripts, projected through exon and CDS
coordinates, and reduced to consequence facts. One generated VEP-116 rule program maps
those facts to SO terms. Numeric rows leave the kernel first; identifiers, HGVS,
supplementary annotations, and strings are joined or rendered later.

```text
model construction, once per declared source

Ensembl core + funcgen relations + matching FASTA
                -> checked prepared relations + deterministic receipt
                -> immutable transcript/exon/sequence + regulation/motif arrays

annotation, repeatedly

sorted ALT alleles
    -> raw VCF span + VEP feature span + minimized sequence edit
    -> continuing transcript and regulation/motif sweeps
    -> topology / splice / projection / sequence / interval-feature facts
    -> one VEP-116 fact-to-SO program
    -> compact consequence rows keyed to transcripts or core interval features
    -> SQL joins, identifiers, HGVS, CSQ/JSON, clinical evidence
```

The data placement follows the cost and lifetime of each fact:

| Data | Where it lives and why |
| --- | --- |
| Ensembl identifiers, attributes, source tables, and receipts | DuckDB relations, because they are provenance and joinable data rather than hot-loop state. |
| Transcript spans, exon maps, prepared CDS/flanks, sparse sequence edits, and compact regulatory/motif intervals | One immutable named C model, because every variant reuses them. |
| Variant alleles | Sorted batches, because genome-scale inputs stream and coordinate order makes candidate discovery nearly linear. |
| dbSNP, ClinVar, gnomAD, scores, and interval tracks | DuckDB/Parquet streams or bounded tiles, because these sources can dwarf the transcript model and may be selected independently. |
| Exon cursors, candidate lists, sequence scratch, and output buffers | One worker-owned workspace, because mutable state must not be shared across DuckDB workers or named models. |

Five rules keep the design understandable:

1. There is one biological authority for independent consequences: fact producers followed
   by the generated VEP-116 rule program. A fast path must prove equality with it.
2. The uploaded VCF span, VEP's feature span, and the minimized sequence edit are retained
   separately because VEP predicates genuinely consume different coordinates.
3. Coordinate order is part of the execution contract. cgranges seeds a run; a continuing
   sweep, not a fresh transcript search per allele, handles the common sorted workload.
4. Missing source facts stay explicit. A partial model, unsupported sequence edit, invalid
   projection, or unknown reference is not silently converted into a supported consequence.
5. Compatibility and speed are measured separately. VEP differentials compare the union of
   all emitted variant/object pairs; throughput reports count input alleles, candidate
   work, emitted rows, bytes, model, threads, and materialization.

Phased haplotypes and structural variants reuse projection, sequence editing, and the SO
rules, but they are not disguised as independent small variants. Haplotypes group several
edits before translation; breakends retain two loci. Supplementary annotation and ACMG/AMP
reasoning remain relational work above the deterministic consequence kernel.

## Responsibility split

DuckDB owns:

- staging Ensembl relations and FASTA, filtering transcripts, verifying receipts, and
  constructing canonical model relations;
- sorting and partitioning variant input;
- exact and interval supplementary annotations, distinct from VEP core regulatory/motif consequences;
- stable biological identifier tables, provenance, final joins, and presentation; and
- later clinical evidence and ACMG/AMP reasoning as inspectable relations.

The C kernel owns:

- checked event geometry and allele trimming;
- sorted transcript and regulatory/motif candidate traversal;
- exon, intron, UTR, splice, CDS, and directional topology;
- genomic-to-transcript/CDS projection;
- reference validation, edit application, translation, and consequence facts;
- structural-event geometry and core RegulatoryFeature/MotifFeature predicates; and
- compact consequence rows over numeric model ordinals.

htslib owns VCF/BCF transport and genotype decoding. HGVS consumes the same projected edit
facts as consequence prediction; it must not reconstruct biology from rendered SO strings.

## Code map

- `functions.yaml` is the public SQL contract.
- `src/duckvep/duckvep_ensembl.c` registers the relation-to-model SQL described below.
- `src/duckvep/duckvep_model.c` validates, owns, publishes, pins, and drops named models.
- `src/duckvep/duckvep_annotate.c` and `src/duckvep/duckvep_variant_tile.c` adapt DuckDB
  vectors to the pure C batch interface and materialize results.
- `src/duckvep/kernel/include/duckvep_kernel.h` is the host-neutral kernel ABI; the sibling
  kernel sources own traversal, topology, projection, coding edits, phased edits,
  structural geometry, and fact-to-SO evaluation.
- `test/sql/duckvep_ensembl.test` exercises model preparation and publication;
  `test/sql/duckvep_annotate.test` exercises the SQL adapter; `test/duckvep/property/` and
  `test/duckvep/conformance/` own pure-C properties and VEP differentials.

## Transcript model build

The model builder is a compiler over relations already staged in DuckDB. It is not a
downloader, a MySQL client, or a parser for Ensembl dump archives. The caller first imports
the pinned Ensembl tables and supplies matching reference sequence chunks. Building and
loading are then separate operations:

```text
Ensembl core tables ─────┐
Ensembl funcgen tables ──┼─> validate and prepare relations ─> receipt
reference FASTA ─────────┘                 │
                                           ├─> region projection
                                           ├─> transcript/exon projection
                                           └─> regulation/motif projection
                                                   │
                                                   v
                                      immutable named C model
                                                   │
                                                   v
                                      sorted consequence queries
```

This separation is deliberate. DuckDB performs the large joins, sequence assembly, and
provenance work once. Annotation workers reuse compact immutable arrays and do not carry a
Perl object graph, stable-ID strings, or source-table metadata through the hot loop.
The four preparation/receipt functions are DuckDB table macros registered through the stable C API;
the C registration file does not iterate transcript rows or implement a second importer.

### Inputs

`duckvep_ensembl_regions(...)` and `duckvep_ensembl_transcripts(...)` read these Ensembl
core relations by name from the supplied schema:

- `coord_system` and `seq_region` identify the requested assembly and its regions;
- `gene`, `transcript`, `exon`, and `exon_transcript` define transcript topology;
- `translation` defines coding start and end within ranked exons; and
- `attrib_type`, `seq_region_attrib`, `transcript_attrib`, and `translation_attrib` supply
  codon tables, consequence-relevant flags, and exceptional sequence edits.

`duckvep_ensembl_regulation_features(...)` reads the release-matched funcgen
`regulatory_feature`, `feature_type`, and `motif_feature` relations. Funcgen uses core
sequence-region IDs; the macro joins them to the already prepared region relation, rejects
missing feature types or invalid coordinates, and assigns one dense ordinal space across
RegulatoryFeature and MotifFeature rows. It also reproduces VEP 116's source selection by
discarding `epigenetically_modified_region` (EMAR) RegulatoryFeature rows before ordinals
are assigned; VEP's database annotation source removes those rows before constructing
overlap objects. Stable IDs, feature types, regulatory-build IDs, binding-matrix IDs,
strand, and scores remain cold DuckDB columns. Only region ordinal, start, end, kind, and
feature ordinal enter the C model.

Attribute values are cast to their semantic text form while importing these relations.
This matters for official dumps where a generic CSV/Parquet staging pass may infer a
numeric-only `value` column as an integer; the model compiler does not require callers to
rewrite an otherwise valid staged schema solely to satisfy string predicates.

The reference relation has `(chrom, start, end, seq)`, where `start` is zero-based and
`end` is half-open. It is normally persisted from tiled `fasta_nuc(..., include_seq :=
true)` output. A reduced FASTA deliberately builds a reduced model; the builder never
silently fills absent contigs from another source.

GRCh37 is a separate Ensembl source, not a coordinate option applied to the GRCh38
genebuild. Release 116 publishes `homo_sapiens_core_116_37` on the GRCh37 archive; it uses
the release-116 schema around the frozen release-75/GENCODE-19 annotation. It has GENCODE
attributes but no MANE data, because MANE is defined only on GRCh38. A GRCh37 model must
therefore use matching GRCh37 core relations and reference sequence and must not invent a
MANE mapping. See https://grch37.ensembl.org/index.html and
https://www.ensembl.org/info/genome/genebuild/mane.html.

### Region preparation

`duckvep_ensembl_regions(...)`:

1. verifies that every reference chunk is non-null, has valid coordinates, and contains
   exactly `end - start` bases;
2. verifies that chunks cover each supplied contig continuously from zero;
3. matches each FASTA contig to exactly one same-name, same-length Ensembl region on the
   requested assembly and species; and
4. assigns dense model-local region ordinals from zero while retaining the Ensembl source
   identifiers and coordinate-system fields.

The model contract bounds a region coordinate at `UINT32_MAX` and the number of regions at
65,536. A mismatch is a query error; no partial region model is published.

### Transcript and sequence preparation

`duckvep_ensembl_transcripts(...)` applies the VEP-116 core-source filter before assigning
model ordinals: a transcript must be current and have a non-empty stable ID, and neither
the `artifact` biotype nor a `readthrough_tra` attribute is admitted. This is why a model
built from the full core dump has the same candidate transcript population as VEP's core
cache instead of merely the same coordinate source. For each selected transcript it:

1. joins its gene, optional translation, ranked exons, and selected attributes;
2. assigns dense transcript and gene ordinals while retaining Ensembl numeric IDs, stable
   IDs, versions, and biotypes in the prepared relation;
3. calculates each exon's transcript-oriented cDNA span from its rank;
4. projects the translation start and end from exon-relative coordinates into genomic and
   cDNA coordinates;
5. reconstructs exon sequence from the reference chunks, reverse-complements negative-
   strand exons, joins exons in transcript order, and extracts the CDS;
6. applies Ensembl start-phase preparation: phase `-1` or `0` adds no prefix, while a
   positive phase adds one or two leading `N` bases; and
7. reads the sequence region's `codon_table` attribute, defaulting to table 1 exactly as
   VEP does, and rejects conflicting, malformed, or unsupported table IDs;
8. extracts the complete transcript-oriented spliced sequence before and after the CDS;
9. parses supported `initial_met`, `_selenocysteine`, `amino_acid_sub`, and
   `_stop_codon_rt` Translation SeqEdits into a sparse reference-peptide relation; and
10. parses each mature-miRNA cDNA range from the Ensembl `miRNA` transcript attribute and
    projects it through the ranked exons into one or more genomic segments.

The compact flag word records translation presence, protein-coding/NMD/miRNA biotype,
incomplete CDS ends, selenocysteine, stop readthrough, transcript RNA edits, peptide
edits, MANE, GENCODE, CCDS, and upstream-start state. Translation `_rna_edit` is not one
of `Translation::get_all_SeqEdits()` in VEP 116 and is therefore not treated as a peptide
edit. The standalone model ABI still defines a readthrough-transcript flag for explicitly
prepared alternate models, but the VEP-compatible core importer filters those transcripts
before publication.

For GRCh38 MANE attributes, the prepared DuckDB relation also retains the attribute value
as `mane_select_refseq` or `mane_plus_clinical_refseq`. That value is the paired versioned
RefSeq transcript accession. The builder rejects multiple or empty mappings for a MANE
flag. Only the selection bits enter the resident C model; Ensembl/GENCODE and RefSeq
identifiers stay in the cold relation for late SQL projection and future RefSeq HGVS.

Single-position, single-amino-acid Translation SeqEdits are kept with the reference-derived
CDS and applied only to VEP's reference peptide; the alternate peptide remains the raw
codon translation. Transcript RNA edits and other Translation SeqEdit shapes withhold
sequence with an explicit reason. The transcript, coordinates, and flags remain in the
model. Unsupported reference alphabet has the same fail-closed shape. This is different
from malformed topology: bad coordinates, strand, exon phase or rank, translation bounds,
or incomplete sequence reconstruction abort the build.

The resident model stores CDS bytes and both non-coding transcript flanks in two packed
byte pools with offsets and lengths per transcript. The flank pool has no per-transcript
allocation or padding and is cold on ordinary coding rows. It is read only when a VEP
predicate rebuilds the 5-prime-UTR-plus-CDS or CDS-plus-3-prime-UTR string. The prepared
relation retains the old three-base `post_cds_bases` projection for compatibility, but
production models load the complete flanks.

Peptide edits and mature miRNA ranges are prepared once, not remapped for every variant.
Peptide edits are packed by transcript and protein position and consulted only for an
overlapping reference-peptide window. Ensembl stores mature-miRNA ranges as transcript
cDNA intervals, possibly more than one per transcript. The builder splits a range when it
crosses an exon boundary and returns the resulting inclusive genomic segments in
`mature_mirna_regions`. The resident model packs all segment starts and ends into flat
arrays with one offset per transcript. A non-miRNA transcript pays only the transcript-flag
test; a miRNA candidate scans its usually tiny owned slice.

### Prepared relations and publication

The builder returns one row per transcript and one row per core regulation/motif feature.
The transcript row contains the hot fields, source and stable identifiers, biotypes, the
optional prepared CDS, an ordered nested exon list, and nested mature-miRNA genomic
segments. The feature row contains a dense ordinal, interval, kind, and cold source
metadata. Callers persist the region, transcript, and feature relations and derive the
required sorted loader projections plus optional side projections from them:

- region ordinal, and sequence length when complete-coverage claims are requested;
- transcript ordinal, region, span, strand, gene ordinal, flags, optional CDS fields,
  codon table, and complete pre-CDS and post-CDS sequence;
- transcript ordinal plus each exon's genomic span, cDNA span, phase, and end phase;
- for transcripts with mature-miRNA attributes, transcript ordinal plus each projected
  genomic segment start and end; and
- for supported Translation SeqEdits, transcript ordinal, one-based protein position, and
  replacement amino acid; and
- for RegulatoryFeature and MotifFeature, feature ordinal, region ordinal, inclusive start
  and end, and compact feature kind.

The `core_schema` argument names relations, not a transport. It can point at tables loaded
from Ensembl's tab-separated MySQL dumps, or at a read-only MySQL catalog attached through
DuckDB's `mysql` extension. Downloading, attaching, and staging stay outside the model
builder so extension builds remain offline and the same validation runs for either source.
The builder does not require a MySQL server once those relations and the matching reference
chunks have been persisted.

`duckvep_model_receipt(...)` checks dense ordinals, region/transcript agreement, and every
regulatory/motif interval against its declared region. It
records the declared source, release, assembly, transcript filter, source-manifest hash,
reference hash, model counts including CDS, transcript-flank bases, mature-miRNA
transcripts and projected segments, peptide edits, regulatory regions, and motif features,
and a deterministic hash over every hot model field, including mature-miRNA ranges,
reference-peptide edits, and interval-feature geometry. There is no
timestamp: identical declared inputs must produce the same receipt.

The checked-in acceptance fixtures under `test/data/duckvep/ensembl_core/` are about 120
KiB. The GRCh38 fixture contains complete release-116 MT, `KI270395.1`, and
`HG2047_PATCH` source rows;
it covers mitochondrial codon-table and peptide-edit behavior, an ordinary multi-exon CDS,
three real release-116 MotifFeature rows on `KI270395.1`, and the real
`ENST00000715685` ↔ `NM_032790.4` MANE pair. The GRCh37 fixture contains MT and
`GL000201.1`; it proves sequence-backed coding annotation from the archived GENCODE-19
model and the absence of MANE mappings. The explicit staging script verifies both official
core manifests, the GRCh38 funcgen manifest, assembled reference hashes, deterministic
model receipts, and exact model counts before writing Parquet. The component manifest
hashes are folded into one sorted canonical source-manifest hash. Tests never contact
Ensembl.

`duckvep_model_load(...)` reads committed, non-temporary relations through a private
connection, validates and narrows every value, builds independent transcript and
regulation/motif seed indexes, and
only then publishes the named immutable model. A failed load publishes nothing. Several
models may coexist in one database instance, and their numeric ordinals are meaningful
only within their model. Stable IDs and provenance remain DuckDB columns. This is also the
contract for haplotype-resolved or pangenome paths: each model declares its exact assembly
or path set and sequence hashes; contig aliases and mapping confidence remain explicit
relations rather than being guessed from chromosome spelling.

The transcript query accepts the original 11 columns, the legacy 12-column form ending in
three `post_cds_bases`, or the complete 13-column form ending in `pre_cds_sequence` and
`post_cds_sequence`. Only the complete form may resolve length-changing edits crossing a
CDS boundary. Older models remain loadable and return `missing_transcript_flank` rather
than guessing when such a predicate is reached.

The optional `mature_mirna_query` has three columns: transcript ordinal, inclusive genomic
start, and inclusive genomic end. Rows must be ordered by transcript and start. The loader
proves that each range belongs to a miRNA transcript, stays inside one of its exons, and is
ordered before publishing the model. Omitting the query preserves compatibility for
models that do not carry these Ensembl attributes; it does not invent mature regions.

The optional `peptide_edit_query` has transcript ordinal, one-based protein position, and
one uppercase replacement amino acid. Rows must be unique and ordered by transcript and
position. The loader packs them into the immutable sequence model and proves that every
position lies within its prepared peptide before publication.

The optional `interval_feature_query` has five columns ordered by region, start, and dense
feature ordinal: feature ordinal, region ordinal, inclusive start, inclusive end, and kind
(`1` RegulatoryFeature, `2` MotifFeature). The loader narrows these to a separate immutable
SoA and builds its own cgranges seed index. Cold funcgen metadata is joined later by feature
ordinal; it is not copied into each consequence row.

The loader treats transcript coverage as partial by default. Only a model deliberately
loaded with `transcript_coverage_complete := true` may turn “no loaded transcript here”
into supported `intergenic_variant`; a partial model returns an unresolved result instead.

```sql
CREATE TABLE reference_chunks AS
SELECT chrom, start, "end", seq
FROM fasta_nuc('GRCh38.primary.fa', bin_width := 1048576, include_seq := true);

CREATE TABLE model_regions AS
SELECT * FROM duckvep_ensembl_regions(
  'ensembl_core', 'reference_chunks', 'GRCh38'
);

CREATE TABLE model_transcripts AS
SELECT * FROM duckvep_ensembl_transcripts(
  'ensembl_core', 'reference_chunks', 'GRCh38'
);

CREATE TABLE model_regulation AS
SELECT * FROM duckvep_ensembl_regulation_features(
  'ensembl_funcgen', 'model_regions'
);

CREATE TABLE model_receipt AS
SELECT * FROM duckvep_model_receipt(
  'model_regions', 'model_transcripts',
  'Ensembl', '116', 'GRCh38', source_manifest_sha256,
  reference_sha256, 'VEP 116 core transcript selection',
  regulation_features_table := 'model_regulation'
);
```

### What the builder does not implement

The implemented builder does not yet:

- stage Ensembl `table.sql` and dump files into DuckDB;
- reproduce the remaining non-core VEP 116 transcript sources and selection rules,
  including `estgene` and `otherfeatures`/RefSeq;
- apply transcript RNA edits or arbitrary length-changing/range Translation SeqEdits;
- extend full-model receipts and indexed-cache differentials beyond human GRCh37/38 and
  *P. falciparum* to more Ensembl species, assemblies, and codon-table combinations;
- preserve the complete Ensembl xref, protein-feature, supporting-feature, attribute, and
  alternate-transcript relations; or
- import variation, phenotype, or other supplementary annotations. RegulatoryFeature and
  MotifFeature are implemented core consequence inputs, not supplementary annotations.

Those richer facts should remain typed DuckDB relations joined by numeric source IDs. They
do not belong in every resident C transcript record. Exact VEP-compatible selection and
the richer Ensembl relation set can therefore be separate named products built from the
same staged release. Production-corpus conformance and throughput remain tracked at
https://github.com/RGenomicsETL/duckhts/issues/93 and
https://github.com/RGenomicsETL/duckhts/issues/95. Further species and genetic-code
coverage is tracked at https://github.com/RGenomicsETL/duckhts/issues/119.

## Ownership

- A named model owns transcript, exon, sequence, regulation/motif SoAs, and both interval
  indexes and is immutable after publication.
- The registry is tied to one DuckDB database. Dropping a model fails while a worker pins
  it.
- Each worker owns its workspace, exon cursor state, allele pool, result storage, mutable
  faidx handle, and sequence scratch. The model pins open read descriptors for the validated
  FASTA, `.fai`, and optional `.gzi`. Linux workers reopen those descriptors through
  `/proc/self/fd`; Windows workers use the resolved source while the model retains a
  deny-write handle. Other POSIX systems use independent handles on the resolved source and
  verify their identities around lazy open, avoiding `/dev/fd` handles whose seek offsets
  may be shared. External replacement or in-place mutation of a loaded source is outside the
  immutable-model contract on those systems and fails when the identity change is observed.
  Contained reference requests reuse the worker
  window, and bounded forward read-ahead amortizes faidx fetches over coordinate-sorted
  alleles. Workspaces are pooled only after their call completes.
- No mutable FASTA handle, interval cursor, genotype iterator, or result builder is shared
  between workers.

## Independent-variant execution

`duckvep_annotate(model, seq_region, position, reference, alternate
[, distance] | [, upstream_distance, downstream_distance])` is the current stable-C-API
adapter for independent biallelic small variants. Both direction windows default to VEP's
5,000 bases; one optional distance changes both, while two values configure them separately
and zero disables the corresponding direction (or both when it is the single symmetric
value). These distances extend candidate admission only beyond transcript endpoints; they
do not cap an allele span or clip an event that overlaps a transcript. Literal alleles up
to the compact 65,535-byte slice limit retain uploaded, VEP-feature, and minimized-edit
geometry independently. An ordinary minimized deletion that contains a complete transcript
therefore reaches the same `transcript_ablation` fact as a symbolic deletion, while an
equal-length containing span retains VEP's endpoint-UTR comparison behavior. The adapter
copies one DuckDB vector into compact arrays, splits on
model/contig/window/order changes, seeds the first candidate set through cgranges, and
advances independent sorted transcript and regulation/motif sweeps. The SNV point path
keeps a per-transcript exon rank and advances it monotonically; other transcript spans use
the exhaustive classifier. Regulation/motif rows use exact event overlap with no transcript
flank. Both sweeps share the generic interval-candidate helper, but own separate active
sets because their cardinalities differ. Transcript fast/exhaustive paths and the complete
feature sweep are property-checked against independent or brute-force oracles.

The scalar callback has no expression-local state across DuckDB vectors, so it restarts the
sweep at each vector. It is a real batch interface, not the stateful whole-stream contract.
The latter must expose explicit carry state through a stable host lifecycle rather than
assuming DuckDB happens to preserve scalar callback state.

## One consequence authority

The intended decision chain is:

1. event geometry;
2. transcript topology and projection;
3. one edit/CDS/peptide context;
4. consequence facts; and
5. one static fact-to-SO rule table.

Specialized traversal or SNV codon code is acceptable only when an exhaustive/property
oracle proves identical results. It is not a second biological authority.

Production MNVs and indels use the generalized edit/CDS/peptide context exactly once. A
projection, sequence, capacity, or unsupported-state failure remains explicit and is never
retried through a smaller shape-specific classifier. Narrow direct classifiers remain only
as independent pure-C test references; the annotation path cannot select them after a
context failure. https://github.com/RGenomicsETL/duckhts/issues/93 owns continued
real-VEP differential closure, not a second consequence implementation.

## NMD prediction

`NMD_transcript_variant` remains the VEP core consequence for a variant inside a
transcript already imported with the `nonsense_mediated_decay` biotype. It does not say
that the current variant creates a new NMD substrate.

Variant-induced NMD is a separate compact result derived from VEP Plugins release/116
`NMD.pm` (`0082591268417af618e03850c5ffdc7c09998a5d`). Stop-gained, frameshift,
splice-donor, and splice-acceptor consequences are predicted to escape for an intronless
transcript, an early-CDS event (`cds_end <= 101` in the plugin), an event in the last
exon, or an event in the plugin's inclusive 51-base penultimate-exon-end window. An
eligible projected event matching none of those rules is `triggering`; an eligible event
without coding coordinates is `unresolved`. The SQL result exposes each escape reason as
a boolean.

The plugin does not use the minimized sequence edit for those coordinates. Its
`BaseTranscriptVariation::cds_coords` path projects the full VEP `VariationFeature` and
its exon rules read `VariationFeature::seq_region_end`. For an ordinary span, the first
and last mapped coordinates may enclose an internal mapper gap as long as both endpoints
map. For a pure insertion, the parent `TranscriptVariation` instead retains the empty
feature as a reversed CDS range such as `102,101`; one genomic flank may be sufficient at
an exon edge. This is different from the expanded `101..102` allele range rendered in
VEP JSON, and `NMD.pm` reads the parent's lower `cds_end`. Immediately before the first
coding base, the valid parent range is `1,0`: the plugin tests whether both values are
defined, not whether they are nonzero, and classifies the zero end as an early-CDS escape.

Equal-length alleles retain unchanged uploaded bases in their feature. DuckVEP must
therefore preserve both views: the minimized edit changes the CDS, while the full feature
decides NMD position. The early-CDS fact cached by the coding classifier is reusable only
when both genomic spans are identical; insertions and wider features use the feature
projector. This representation-dependent behavior is pinned by paired executable-plugin
witnesses on both strands.

This is the pinned VEP positional policy, not a direct molecular assay or a claim that
every transcript follows one universal NMD rule.

## Phased edits

The pure C mutation core already rebuilds a CDS from several non-overlapping edits in one
reverse-coordinate pass, translates once, and partitions interactions while the frame is
displaced or the next edit touches the same alternate codon. The missing stream groups edits by
`(model, transcript, sample, phase_set, haplotype)` and retains all contributing variant
IDs.

The stream must preserve the original record/ALT identity, decoded allele indexes,
ploidy, phasing flag, and `PS`/`PID`-like phase-set provenance. The same called local
haplotype may arrive as one MNV, several SNVs, or overlapping SNV/indel records; those
representations must yield the same alternate CDS and peptide when they encode the same
phased edits. Independent annotation of each row is not an acceptable substitute. The
existing phased-SNV-versus-MNV property covers one important subset; mixed replacement and
indel representations still require generated equivalence and executable csq/Haplosaurus
corpora.

Long-read callsets make that contract unavoidable. Clair3 or DeepVariant may emit phased
SNVs/MNPs while Sniffles or cuteSV emits a structural record for the same sample and local
haplotype; the biological edit set can cross those record classes. HBA1/HBA2-like paralogous
loci add a separate ambiguity: an aligner/caller may report one of several near-identical
placements. DuckVEP must retain call, alignment, assembly-path, and phase provenance and
must not manufacture certainty by merging records solely because their nominal coordinates
are close. Haplotype-resolved HPRC assemblies can be loaded as separately receipted models
or explicit paths; comparison to GRCh37/38 remains a mapping relation, not a chromosome-name
alias.

The annotation-model receipt is not enough to reproduce a long-read result. A callset
receipt must also identify the read chemistry, basecaller and model, alignment reference
and aligner, small-variant/SV caller versions and options, phasing method, and any callable
or confidence masks. Agreeing `GT` and phase-set fields do not by themselves make
overlapping records compatible. The phased executor must prove that their reference and
alternate sequences form one consistent local haplotype or return an explicit edit-conflict
result while retaining every source record.

An assembled HPRC haplotype is useful truth evidence, but it is not automatically a
DuckVEP model. Path-coordinate annotation additionally requires transcripts projected or
annotated on that exact path, mapping confidence, and locus-level assembly QC; a nominally
haplotype-resolved assembly can still carry a flagged collapse or misassembly. Incremental
annotation therefore invalidates work by dependency: an independent new allele can be
annotated alone, while a changed phase block, caller interpretation, transcript projection,
or model receipt requires recomputing the affected transcript haplotypes.

Sorted input bounds lifetime: a transcript's phased state can be finalized once the stream
passes its end. Reference paths stay implicit; non-reference paths should share compact
edit prefixes across samples and translate each distinct leaf once. GT/PS decoding,
arbitrary ploidy, transcript-close flushing, and VEP Haplosaurus comparison are tracked at
https://github.com/RGenomicsETL/duckhts/issues/92.

DuckDB may decode genotypes, derive explicit sample/phase-set/haplotype columns, and sort
the input relation by genomic coordinate before streaming it. That SQL grouping is not the
same as biological edit interaction: the C executor still owns transcript lifetime,
same-codon interactions, an open displaced frame, frame restoration, and the complete list
of contributing variants. Benchmarks must include both an already-sorted stream and the
same workload with DuckDB sorting included.

## Structural events

Small edits and structural events share overlap, projection, provenance, and output, but
not one lossy event shape. Deletion, duplication, inversion, CNV, breakend pairs, inserted
sequence, and repeat changes retain their typed geometry. A breakend is two loci, not a
wide interval or symbolic point.

`duckvep_annotate_sv(...)` and `duckvep_annotate_sv_compact(...)` accept exact single-locus
events plus a typed DEL, DUP, tandem-DUP, tandem-repeat (`STR`), INV, INS, CNV, or unknown operation and an
explicit loss/neutral/gain/unknown copy direction. Span operations use one-based inclusive
start/end coordinates. An insertion uses `start = end = P` for the interbase site after
reference base `P`; preparing symbolic VCF therefore removes the left anchor, maps a span
to `start = POS + 1, end = INFO/END`, and maps an insertion to `P = POS`. The adapter
rejects contradictory operation/direction pairs. This single-locus adapter rejects BND
because it is not a single-locus event; the dedicated two-locus adapters below handle it.

VEP expands a bounded `<CNV:TR>` from `RN` plus `RUS`/`RUC` or `RB` into literal alleles
when the result fits its configured structural-size limit. Such an event is then an
ordinary `VariationFeature` and enters DuckVEP's small-variant path. An oversized or
unexpanded repeat remains a structural `tandem_repeat`. DuckVEP's `STR` type preserves
that identity for provenance and later HGVS while reproducing VEP's tandem-duplication
gain/insertion predicates. Raw repeat units and counts remain columns in the surrounding
relation; the consequence kernel consumes the prepared literal allele or exact structural
span, not parser-specific INFO strings.

`duckvep_annotate_breakend(...)` and its compact form accept the local and mate regions and
raw one-based VCF positions in one row. They query the resident cgranges transcript and
regulatory/motif indexes around both loci, merge and deduplicate each object class, and
call shared C evaluators with explicit variant/object pairs. This is intentionally
separate from the sorted single-locus sweep: neither a wide span nor two independent
endpoint calls can reproduce VEP 116.

The exact VEP state is asymmetric. `BaseVCF4::get_start` moves the local BND `POS` to
`POS + 1`; `StructuralVariationFeature::_parse_breakends` retains the mate coordinate.
Ordinary intron, exon, UTR, splice, and coding predicates use only that shifted local
feature. Candidate discovery also creates an overlap allele for the mate, but the mate
changes transcript consequences only through the exceptional `feature_truncation`
predicate. VEP may therefore emit two internal overlap-allele rows for one BND and
transcript. DuckVEP returns their consequence-set union once. A transcript reached only
through the mate has a zero region mask and NULL rich region because no local topology
exists.

VEP's `RegFeat` lane also evaluates both points. A RegulatoryFeature or MotifFeature hit
by the shifted local point, the verbatim mate point, or both produces one DuckVEP object
row. The result is asymmetric: a shifted-local hit retains
`regulatory_region_variant` or `TF_binding_site_variant`; a mate-only hit takes VEP's
generic HIGH-impact `feature_truncation` chromosome-breakpoint branch without requiring
deletion or copy-number loss. Once the mate has discovered an object exactly, VEP also
attaches a shifted local point on the same contig when it is outside but within the fixed
5000-base structural-feature admission distance. That local allele falls back to
`intergenic_variant`, so the object-level result is
`feature_truncation&intergenic_variant`. A close point does not discover an object by
itself. If both points hit one object exactly, the local base term wins. VEP may
materialize duplicate identical rows or distinct allele-level rows for one object;
DuckVEP's public contract is their consequence-set union once per `(event, object)`.

The surrounding DuckDB relation retains event identity, bracket orientation, raw ALT, and
provenance for HGVS, fusion, and round-trip consumers. Orientation does not change the
transcript consequence set, so it is not an ignored kernel argument. The C lane preserves
VEP's fixed 5 kb endpoint admission cap in addition to the configured directional
upstream/downstream distances. These are independent controls despite sharing the same
default number. Raising the caller distance to 10 kb may widen upstream/downstream
transcript terms, but cannot turn an endpoint 5,001 bases from the same transcript or
interval feature into a `StructuralVariationOverlapAllele`. Interval-feature candidate
discovery remains exact; the fixed-cap endpoint is attached only after the mate has
discovered that same object. This does not clamp ordinary transcript predicates: an
overlap allele created by the mate still reads the shifted local feature, so a 10 kb
caller window can emit an upstream/downstream term for a local point beyond the fixed
allele-admission cap. Pure-C, SQL, and R regressions pin 5,000 versus 5,001 bases
under both a 10 kb caller distance and a zero caller distance. In the latter case an
admitted local transcript allele has no directional predicate and falls back to
`intergenic_variant`, which is unioned with a mate-derived `feature_truncation`. Randomized
sweep scenes include zero, 1, 50, 100, 4,999, 5,000, 5,001, 10,000, and 65,535-base
windows rather than treating 5 kb as a maximum allocation size. A seeded executable
differential covering chromosomes 1,
2, 7, 21, and X, all four bracket orientations, same- and cross-chromosome pairs, and
transcript/exon/intron/CDS/flank endpoint states matched all 91,428 transcript pairs from
1,004 generated events.

That differential isolates every BND in VEP with `buffer_size=1`. VEP 116 otherwise puts
mate coordinates from several records into one chromosome-blind interval tree, so a
neighboring BND can change the executable oracle's transcript set. One VEP process is
retained for cache reuse; only the semantic event buffer is isolated. Chromosomes stay
contiguous and positions increase within each chromosome in the generated VCF.

The public SO mask now binds all 41 terms registered by VEP 116. Six regulatory-region and
transcription-factor-binding-site terms are produced by a separate interval-feature
evaluator over the same event geometry. Their five hot columns live in a separate resident
SoA, not in fake transcript rows; cold funcgen metadata remains an ordinary DuckDB relation.
The small-variant and exact single-locus structural adapters advance the transcript and
feature sweeps together and emit typed rows distinguished by `overlap_object`. Regulatory
or motif output therefore does not require a SQL range join, and it preserves the same
resumable output cursor as transcript output. `sequence_variant`, the forty-first registry
term, has no VEP overlap predicate and is retained as registry metadata rather than emitted
to hide an incomplete model.

Raw BND ALT parsing and repeat-unit/count expansion remain outside the consequence kernel.
VEP 116 stores `CIPOS`/`CIEND` and structural inserted-sequence payloads, but its registered
41-term consequence predicates use nominal `POS`/`END` and its structural
`inframe_insertion` branch explicitly has no inserted-sequence implementation. Those
fields therefore remain required provenance and future HGVS/round-trip inputs, not missing
consequence facts. The BND statistical differential is part of the executable-VEP harness;
broader chromosome, species, and real fusion corpora remain continuing evidence rather
than a second implementation. Remaining structural follow-up is tracked at
https://github.com/RGenomicsETL/duckhts/issues/98.

Callers such as Sniffles and cuteSV may provide `CIPOS`/`CIEND`, mate identity, inserted
sequence, copy number, and several records for one event. For strict VEP-116 consequence
parity, an upstream relation may annotate the nominal `POS`/`END` while retaining every
confidence interval and payload beside the result. It must not relabel the nominal span as
experimentally exact, discard those fields, or reuse it as HGVS/fusion geometry. Likewise,
ambiguous placement between paralogous loci is source evidence for downstream SQL; choosing
one locus is not a consequence-kernel inference.

## HGVS

HGVS is a consumer of the projected edit, not a formatter over consequence names. The
current event keeps the uploaded span, VEP feature span, minimized REF/ALT edit, insertion
boundary, and anchor side separately. The transcript model keeps complete spliced
pre-CDS/CDS/post-CDS sequence, and the sequence layer can apply one edit or a grouped edit
set and compare the resulting peptide. Those facts must remain stable while phased work
is added.

The implemented internal layer makes that ownership explicit:

- `duckvep_model_open(...)` owns the canonical prepared transcript, exon, and sequence
  views used by both consequence and HGVS. It derives the first complete reference stop
  once per coding transcript, or validates a supplied immutable cache, so per-row HGVSp
  does not rescan an unchanged CDS;
- consequence evaluation and HGVS both consume the same prepared allele and CDS edit-set
  helpers; the hot consequence path does not construct a wider formatting object;
- `duckvep_transcript_edit_t` is an HGVS-facing carrier built from that prepared allele and
  one transcript. Projection first adds VEP's endpoint-clipped transcript-slice
  coordinates without trimming, reinterpreting, or clamping the semantic allele; CDS
  projection is attached lazily only when reference validation or protein replay needs it;
- typed transcript-coordinate facts distinguish coding `c.` from non-coding `n.` edits,
  preserve insertion and range geometry, and render into caller-owned bounded buffers;
- VEP-compatible transcript 3-prime shifting consumes VEP's exact constrained +/-1000
  genomic-reference view, while complete uploaded-REF validation and adjacent duplication
  lookup use a separately bounded wider view. Neither creates a second allele authority or
  changes the shift bytes; and
- typed protein facts describe equality, substitution, deletion, insertion, delins,
  duplication, frameshift, start loss, and extension before any string is allocated. A
  single CDS edit is exposed as a virtual alternate sequence, and first-stop scans traverse
  unchanged reference spans plus the edit payload rather than materializing the full
  alternate CDS.

`duckvep_annotate_hgvs(...)` exposes those mechanics as a cumulative SQL result: the first
16 fields are the compact consequence row, followed by transcript/protein suffixes,
transcript-direction shift, and separate structured status/reason fields. Stable versioned
transcript and translation identifiers remain ordinary prepared-model columns and are
joined by transcript ordinal rather than copied into every hot result row. A model may bind
an existing indexed FASTA through an exact sequence-region ordinal/name/length relation;
the loader validates the index without creating or modifying it, while each annotation
worker owns its faidx handle and reusable reference-window storage. One bounded fetch
supplies distinct borrowed views for shifting and lookup. The adapter performs
one candidate sweep and renders into DuckDB-owned output vectors. One worker-owned scratch
buffer handles transcript and protein strings in a single render pass; it grows and retries
only when a result does not fit.

Consequence predicate flags are reusable evidence, not a closed-world serialization of the
later HGVS replay. Positive frameshift evidence can complete a length-changing delta, and
absence of frameshift is conclusive for a length-preserving CDS edit. A length-changing
splice-overlapping edit without that positive flag must run the complete delta evaluator;
otherwise VEP-compatible frameshifts can be misrendered as premature stops or delins.

This implemented surface is not yet a full VEP-HGVS compatibility claim. It covers
independent literal small variants and returns explicit unresolved reasons when reference,
projection, transcript flank, tail, or protein facts are unavailable. Fixed executable
VEP witnesses cover position-one right anchoring, endpoint-clipped deletions/delins,
terminal insertion states, both transcript strands, and VEP's short alternate-CDS trimming
bug. A strict chromosome-21 ClinVar run compared 56,998 transcript pairs: 20,782 matched
both HGVSc and HGVSp, 24,089 matched HGVSc with HGVSp absent on both sides, and 12,127 had
both strings absent, with zero discordant, unresolved, missing, or extra HGVS rows. That is
evidence for the exercised GRCh38 independent-event distribution, not for untested species,
assemblies, transcript-source quirks, or other HGVS classes.

The execution split is:

- genomic `g.`/`m.` HGVS is computed once per allele and needs a genomic reference-window
  provider; `VariationFeature::get_all_hgvs_genomic` uses
  `get_3prime_seq_offset`;
- ordinary Ensembl transcript `c.`/`n.` HGVS is per transcript, but VEP first runs
  `TranscriptVariationAllele::_genomic_shift` over a forward-reference window in that
  transcript's strand direction. Its `perform_shift` loop has a hard-coded 1000-base
  search and retains its own allele-length-dependent loop limit. The shifted genomic
  event is then projected back to transcript coordinates. Complete uploaded-REF checking
  and `hgvs_variant_notation` duplication-source comparison are separate consumers of a
  wider bounded reference lookup and must not widen that exact shift slice;
- RefSeq transcripts with RNA-edit attributes may reject reuse of that genomic shift and
  rerun `perform_shift` over the edited transcript sequence. This is a separate model
  capability, not the default Ensembl path; and
- protein `p.` HGVS consumes the alternate peptide difference produced by the same edit
  set used for phased consequence classification.

The hot transcript sweep therefore emits or retains numeric projected-edit facts; it does
not allocate HGVS strings. Rendering is late, after filtering. Before a public phased API
is fixed, the same prepared CDS edit set must become the phased executor's input rather
than a second trimming or projection authority. External VEP-116 differentials must expand
from the current Ensembl/GRCh38 independent-event coverage to RefSeq RNA edits, other
assemblies/species, all shift modes, exact structural events, and compound edits.
Apply-then-diff sequence equivalence remains an independent property oracle. Exact
structural HGVS can later consume typed exact events. BND HGVS
additionally needs the paired relation's mate and orientation facts; imprecise structural
events remain unsupported until their confidence geometry is represented.

The pinned
[ferro-hgvs v0.9.0 source](https://github.com/fulcrumgenomics/ferro-hgvs/tree/278e2c11134e3b49067d0c334f650c7c29db9cbe)
is an independent HGVS-spec oracle and a useful model for structured fuzzing, large ClinVar
corpora, and hermetic reference fixtures. It is not a DuckHTS dependency and does not
replace the VEP 116 executable as the compatibility authority: canonical HGVS,
Mutalyzer/biocommons behavior, and VEP's historical output can legitimately disagree.
The intended differential therefore has three independent observations: exact VEP-116
output, ferro-hgvs parsing/normalization where its contract applies, and DuckVEP's
apply-then-diff sequence replay.

## Supplementary annotations

Exact sources such as dbSNP, ClinVar, and gnomAD should remain sorted numeric-key streams;
range sources should remain interval streams. They share prepared variant keys and final
output, not a universal physical format. DuckDB/Parquet pruning, parallel reads, memory
management, and spill are the first implementation. A custom cache requires a measured
workload that those facilities do not meet.

Variant-level exact and positional joins happen before transcript expansion. The compact
consequence row already carries transcript and gene ordinals, so gene relations are joined
after consequence expansion without repeating a gene lookup inside C. Ensembl
RegulatoryFeature and MotifFeature are core VEP inputs and run in the C consequence kernel.
Protein domains and genuinely supplementary interval sources use range joins or sorted
interval streams, followed by the smallest source-specific predicate needed for a joined
pair. Strings and JSON remain the final projection.

The stateful stable-API integration is tracked at
https://github.com/RGenomicsETL/duckhts/issues/94. Supplementary source plumbing does not
belong as ignored arguments in the consequence-kernel API.

## Validation and performance

- `make test-duckvep-kernel` runs fixed and randomized pure-C properties; ASan and UBSan
  execute the same suite.
- `make test-duckvep-kernel-statistical` raises randomized targets to an explicit seed and
  trial count.
- `DUCKVEP_PROPERTY_ARGS='-t hgvs'` and `'-t haplotype'` are diagnostic filters for local
  iteration. An official history run leaves this argument empty so no favorable family can
  replace the complete state-machine denominator.
- Randomized properties emit named distribution counters in addition to pass/fail status.
  `test/duckvep/conformance/property_history.R` records those counters in a long-form
  append-only ledger, so a passing run cannot hide that a declared edit shape, strand,
  terminal state, shift limit, or haplotype interaction received zero observations.
- The regulation/motif lane compares the complete resumable cursor, at one output row of
  capacity, with a brute-force all-event/all-feature evaluator for randomized small alleles
  and exact structural events. The offline Ensembl fixture independently imports and emits
  the three real release-116 MotifFeature rows on `KI270395.1`.
- `make test-duckvep-differential` compares generated witnesses to pinned VEP 116.
- `make duckvep-corpus-differential` records the union of emitted variant/transcript or
  variant/core-feature pairs, including mismatches, misses, extras, and unresolved rows.
- The corpus runner's small-event mode samples SNVs, MNVs, insertion-like alleles, and
  deletion-like alleles independently by length-change bin. It can split multiallelic ALT
  rows without rewriting genotypes, stratify by raw allele length through greater-than-
  10-kb representations, checksum the complete source, and emit a source-eligibility
  receipt separately from executable-VEP agreement. Structural mode either reads
  a symbolic VCF or generates seeded DEL/DUP/tandem-DUP/INV/CNV/INS events from real model
  geometry at transcript, exon, intron, splice, UTR, CDS, start-codon, and stop-codon
  states. Breakend mode generates both local-anchor-removed and verbatim-mate points at
  RegulatoryFeature and MotifFeature starts, midpoints, and ends when `--regulatory` is
  enabled, in addition to transcript-derived points. Independent seeds can run
  concurrently against one read-only attached model.
- A held-out pure-C run with seed `20260719` executed 100,000 trials per randomized
  property (175 tests; 204,759 assertions) and passed under the ordinary,
  AddressSanitizer, and UndefinedBehaviorSanitizer targets. It found and minimized a rare
  terminal-codon oracle error after 93,064 generated frame-changing cases; the pinned VEP
  source showed that the production kernel was correct, and the corrected oracle retained
  both concrete-local-peptide and endpoint-reconstruction coverage. The independent
  executable-VEP run with seed `20260716` compared 100,268 generated/fixed alleles and all
  100,268 transcript pairs were exact. Its generated set contained SNVs, MNVs, insertions,
  deletions, and delins; duplicate rejection deliberately makes the accepted shape counts
  non-uniform rather than resampling a cosmetically balanced result.
- The current eight-seed GRCh38 structural campaign generated 40,375 events on chromosomes
  1, 2, 6, 11, 17, 21, 22, and X and compared 2,140,911 emitted transcript pairs with the
  executable indexed VEP 116 cache. Every pair was resolved and exact, with no missing or
  extra emission. This is evidence for the sampled exact single-locus states, not a claim
  for imprecise coordinates, STR payloads, or untested species.
- The paired-BND campaign generated 1,004 same- and cross-chromosome events on chromosomes
  1, 2, 7, 21, and X and matched all 91,428 isolated executable-VEP transcript pairs.
  BND oracle buffers contain one event because VEP's chromosome-blind mate-coordinate tree
  otherwise makes a record's output depend on neighboring BNDs.
- Executable-VEP `--regulatory` comparisons cover a 1,196-site chromosome-21 GIAB sample
  (14,955 annotation-object pairs) and 2,700 generated exact structural events
  (120,224 pairs). Both are exact with no unresolved, extra, missing, or discordant row.
  The structural corpus observes every one of the six regulatory/motif SO terms across
  exact, containing, and partial feature geometries. These are declared sampled
  distributions, not evidence for imprecise SV coordinates or untested funcgen releases.
- `make bench-duckvep-release-parquet` reads the official Ensembl variation consequence
  VCF through typed CSQ columns and records complete versus consequence-only Parquet size,
  checksum, cardinality, and elapsed time without committing the large artifacts.
- Throughput reports must state input variants, candidate/object pairs, output rows and
  bytes, threads, model, and machine. Site-level rate and cohort haplotype-update rate are
  different measurements.
- Comparisons with FastVEP report two timings over the same source and output contract:
  annotation from already coordinate-sorted input, and read-plus-sort-plus-annotation.
  Sorting cost is not hidden, but it is not charged only to DuckVEP either. DuckDB may sort
  once and stream ordered chunks directly into the consequence or phased executor.

HGVS and haplotype performance use cumulative lanes over identical prepared input:

| Lane | Required work |
| --- | --- |
| consequence compact | candidate discovery, consequence facts, compact rows |
| consequence + HGVS facts | the compact lane plus transcript-edit and typed DNA/protein facts, without strings |
| consequence + rendered HGVS | the preceding lane plus accession lookup and rendered bytes |
| bcftools local-CSQ projection | independent-event consequences plus its declared local output fields |
| phased edit sets | genotype grouping, active transcript state, unique alternate leaves, combined translation, and attributed rows |

Every rendered record states input records and ALT alleles, candidate pairs, output rows,
rendered bytes, thread count and core pinning, model/corpus identity, and source revision.
Sorting-plus-execution is recorded separately from execution over already ordered input.
A single independent edit is never relabelled as a haplotype benchmark.

The generated reports are the numeric authority. At the latest tested ancestor for each
corpus, the declared GRCh38 dbSNP, GIAB, ClinVar coding, and ClinVar cross-chromosome
samples, GRCh37, and *P. falciparum* are exact against their indexed caches. The two
GRCh38 ClinVar samples contain 604,233 transcript pairs together, with no unresolved,
extra, missing, or discordant row. The separate NMD-plugin differential is exact for all
68,554 eligible transcript pairs, including the 29,416 states both implementations leave
unresolved.

On the final 644,427-transcript GRCh38 model, the current adjacent one-core compact
measurements are about 841,000 input alleles/s for the transcript-only GIAB topology sample
and 801,000/s when the same run also loads all 1,383,580 admitted regulatory/motif features
and emits their overlaps. The latter produces about 9.36 million output rows/s. The nearest
recorded coding-focused baselines are about 332,000/s for repeated coding SNVs, 112,000/s
for repeated coding non-SNVs, and 135,000/s for the repeated mixed coding set. The
paired-breakend lane reaches about 74,000 semantic events/s and 6.92 million emitted
transcript rows/s. These include stable-API list materialization and aggregation but exclude
model load and input staging. The one-million-input-allele target is not met on the final
model; output-row rates are a second denominator, not a substitute for the site-rate target.
Exact revisions, checksums, row counts, resource measurements, and conditions live in the
generated [conformance report](../benchmarks/duckvep_conformance.md) and
[throughput report](../benchmarks/duckvep_throughput.md).

Open conformance and throughput work is tracked at
https://github.com/RGenomicsETL/duckhts/issues/93 and
https://github.com/RGenomicsETL/duckhts/issues/95.
