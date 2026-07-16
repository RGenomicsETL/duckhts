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
edge projections, not the internal representation. VEP 116 is the behavioral authority;
pure-C properties and bcftools csq supply independent checks for mechanics and phased edit
state.

## Mental model

DuckVEP does two different jobs at two different times.

First, it compiles Ensembl relations and matching reference sequence into an immutable
execution model. DuckDB keeps the rich source relations, stable identifiers, attributes,
and hashes. The compiler checks them, assembles transcript-oriented sequence, and publishes
compact C arrays containing only facts used repeatedly during annotation.

Second, it runs a sorted stream of variant alleles through that model. Each allele is
interpreted once, matched to candidate transcripts, projected through exon and CDS
coordinates, and reduced to consequence facts. One generated VEP-116 rule program maps
those facts to SO terms. Numeric rows leave the kernel first; identifiers, HGVS,
supplementary annotations, and strings are joined or rendered later.

```text
model construction, once per declared source

Ensembl core relations + matching FASTA
                -> checked prepared relations + deterministic receipt
                -> immutable transcript/exon/sequence arrays

annotation, repeatedly

sorted ALT alleles
    -> raw VCF span + VEP feature span + minimized sequence edit
    -> continuing transcript sweep
    -> topology / splice / projection / sequence facts
    -> one VEP-116 fact-to-SO program
    -> compact transcript rows
    -> SQL joins, identifiers, HGVS, CSQ/JSON, clinical evidence
```

The data placement follows the cost and lifetime of each fact:

| Data | Where it lives and why |
| --- | --- |
| Ensembl identifiers, attributes, source tables, and receipts | DuckDB relations, because they are provenance and joinable data rather than hot-loop state. |
| Transcript spans, exon maps, prepared CDS/flanks, and sparse sequence edits | One immutable named C model, because every variant reuses them. |
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
   all emitted variant-transcript pairs; throughput reports count input alleles, candidate
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
- exact and interval supplementary annotations;
- stable biological identifier tables, provenance, final joins, and presentation; and
- later clinical evidence and ACMG/AMP reasoning as inspectable relations.

The C kernel owns:

- checked event geometry and allele trimming;
- sorted transcript candidate traversal;
- exon, intron, UTR, splice, CDS, and directional topology;
- genomic-to-transcript/CDS projection;
- reference validation, edit application, translation, and consequence facts;
- structural-event geometry; and
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
Ensembl core tables ──┐
                      ├─> validate and prepare relations ─> receipt
reference FASTA ──────┘                 │
                                        ├─> region projection
                                        ├─> transcript projection
                                        └─> exon projection
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
The three builder functions are DuckDB table macros registered through the stable C API;
the C registration file does not iterate transcript rows or implement a second importer.

### Inputs

`duckvep_ensembl_regions(...)` and `duckvep_ensembl_transcripts(...)` read these Ensembl
core relations by name from the supplied schema:

- `coord_system` and `seq_region` identify the requested assembly and its regions;
- `gene`, `transcript`, `exon`, and `exon_transcript` define transcript topology;
- `translation` defines coding start and end within ranked exons; and
- `attrib_type`, `seq_region_attrib`, `transcript_attrib`, and `translation_attrib` supply
  codon tables, consequence-relevant flags, and exceptional sequence edits.

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

The builder returns one row per transcript. The row contains the hot transcript fields,
source and stable identifiers, biotypes, the optional prepared CDS, an ordered nested
exon list, and the nested mature-miRNA genomic segments. Callers persist two canonical
tables and derive three required sorted loader projections plus optional side projections
from them:

- region ordinal, and sequence length when complete-coverage claims are requested;
- transcript ordinal, region, span, strand, gene ordinal, flags, optional CDS fields,
  codon table, and complete pre-CDS and post-CDS sequence;
- transcript ordinal plus each exon's genomic span, cDNA span, phase, and end phase;
- for transcripts with mature-miRNA attributes, transcript ordinal plus each projected
  genomic segment start and end; and
- for supported Translation SeqEdits, transcript ordinal, one-based protein position, and
  replacement amino acid.

The `core_schema` argument names relations, not a transport. It can point at tables loaded
from Ensembl's tab-separated MySQL dumps, or at a read-only MySQL catalog attached through
DuckDB's `mysql` extension. Downloading, attaching, and staging stay outside the model
builder so extension builds remain offline and the same validation runs for either source.
The builder does not require a MySQL server once those relations and the matching reference
chunks have been persisted.

`duckvep_model_receipt(...)` checks dense ordinals and region/transcript agreement. It
records the declared source, release, assembly, transcript filter, source-manifest hash,
reference hash, model counts including CDS, transcript-flank bases, mature-miRNA
transcripts and projected segments, and peptide edits, and a deterministic hash over every
prepared model field, including mature-miRNA ranges and
reference-peptide edits. There is no
timestamp: identical declared inputs must produce the same receipt.

The checked-in acceptance fixtures under `test/data/duckvep/ensembl_core/` are about 116
KiB. The GRCh38 fixture contains complete release-116 MT and `HG2047_PATCH` source rows;
it covers mitochondrial codon-table and peptide-edit behavior, an ordinary multi-exon CDS,
and the real
`ENST00000715685` ↔ `NM_032790.4` MANE pair. The GRCh37 fixture contains MT and
`GL000201.1`; it proves sequence-backed coding annotation from the archived GENCODE-19
model and the absence of MANE mappings. The explicit staging script verifies both official
dump manifests, assembled reference hashes, deterministic model receipts, and exact model
counts before writing Parquet. Tests never contact Ensembl.

`duckvep_model_load(...)` reads committed, non-temporary relations through a private
connection, validates and narrows every value, builds the transcript interval index, and
only then publishes the named immutable model. A failed load publishes nothing. Several
models may coexist in one database instance, and their numeric ordinals are meaningful
only within their model. Stable IDs and provenance remain DuckDB columns.

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
- import variation, phenotype, regulatory, or other supplementary annotations.

Those richer facts should remain typed DuckDB relations joined by numeric source IDs. They
do not belong in every resident C transcript record. Exact VEP-compatible selection and
the richer Ensembl relation set can therefore be separate named products built from the
same staged release. Production-corpus conformance and throughput remain tracked at
https://github.com/RGenomicsETL/duckhts/issues/93 and
https://github.com/RGenomicsETL/duckhts/issues/95. Further species and genetic-code
coverage is tracked at https://github.com/RGenomicsETL/duckhts/issues/119.

## Ownership

- A named model owns all transcript, exon, sequence, and interval-index storage and is
  immutable after publication.
- The registry is tied to one DuckDB database. Dropping a model fails while a worker pins
  it.
- Each worker owns its workspace, exon cursor state, allele pool, result storage, and
  sequence scratch. Workspaces are pooled only after their call completes.
- No mutable FASTA handle, interval cursor, genotype iterator, or result builder is shared
  between workers.

## Independent-variant execution

`duckvep_annotate(model, seq_region, position, reference, alternate
[, distance] | [, upstream_distance, downstream_distance])` is the current stable-C-API
adapter for independent biallelic small variants. Both direction windows default to VEP's
5,000 bases; one optional distance changes both, while two values configure them separately
and zero disables the corresponding direction (or both when it is the single symmetric
value). The adapter copies one DuckDB vector into compact arrays, splits on
model/contig/window/order changes, seeds the first candidate set through cgranges, and
advances a sorted C sweep. The SNV point path keeps a per-transcript exon rank and advances
it monotonically; other spans use the exhaustive classifier. Both use the same splice
predicates and are property-checked for equality.

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
events plus a typed DEL, DUP, tandem-DUP, INV, INS, CNV, or unknown operation and an
explicit loss/neutral/gain/unknown copy direction. Span operations use one-based inclusive
start/end coordinates. An insertion uses `start = end = P` for the interbase site after
reference base `P`; preparing symbolic VCF therefore removes the left anchor, maps a span
to `start = POS + 1, end = INFO/END`, and maps an insertion to `P = POS`. The adapter
rejects contradictory operation/direction pairs. BND is rejected because it is not a
single-locus event.

`duckvep_annotate_breakend(...)` and its compact form accept the local and mate regions and
raw one-based VCF positions in one row. They query the resident cgranges transcript index
around both loci, merge and deduplicate the transcript candidates, and call the shared C
evaluator with explicit variant/transcript pairs. This is intentionally separate from the
sorted single-locus sweep: neither a wide span nor two independent endpoint calls can
reproduce VEP 116.

The exact VEP state is asymmetric. `BaseVCF4::get_start` moves the local BND `POS` to
`POS + 1`; `StructuralVariationFeature::_parse_breakends` retains the mate coordinate.
Ordinary intron, exon, UTR, splice, and coding predicates use only that shifted local
feature. Candidate discovery also creates an overlap allele for the mate, but the mate
changes transcript consequences only through the exceptional `feature_truncation`
predicate. VEP may therefore emit two internal overlap-allele rows for one BND and
transcript. DuckVEP returns their consequence-set union once. A transcript reached only
through the mate has a zero region mask and NULL rich region because no local topology
exists.

The surrounding DuckDB relation retains event identity, bracket orientation, raw ALT, and
provenance for HGVS, fusion, and round-trip consumers. Orientation does not change the
transcript consequence set, so it is not an ignored kernel argument. The C lane preserves
VEP's fixed 5 kb endpoint admission cap in addition to the configured directional
upstream/downstream distances. A seeded executable differential covering chromosomes 1,
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
evaluator over the same event geometry. Regulatory and motif features remain ordinary
DuckDB relations for range joins; they are not inserted into the resident transcript
arrays. Importing those Ensembl regulation relations and exposing the joined-pair SQL
adapter remain open. `sequence_variant`, the forty-first registry term, has no VEP overlap
predicate and is retained as metadata rather than emitted to hide an incomplete model.

Raw BND ALT parsing, inserted-sequence payloads, imprecise confidence intervals, and STR
repeat-unit/count preparation remain outside the consequence kernel. The BND statistical
differential is now part of the executable-VEP harness; broader chromosome, species, and
real fusion corpora remain continuing evidence rather than a second implementation.
STR records whose repeat unit and count expand to an ordinary small edit should enter the
existing small-variant path once per allele; oversized or underspecified repeats still
need a typed structural input. This work is tracked at
https://github.com/RGenomicsETL/duckhts/issues/98.

## HGVS

HGVS is a consumer of the projected edit, not a formatter over consequence names. The
current event keeps the uploaded span, VEP feature span, minimized REF/ALT edit, insertion
boundary, and anchor side separately. The transcript model keeps complete spliced
pre-CDS/CDS/post-CDS sequence, and the sequence layer can apply one edit or a grouped edit
set and compare the resulting peptide. Those facts must remain stable while phased work
is added.

The execution split is:

- genomic `g.`/`m.` HGVS is computed once per allele and needs a genomic reference-window
  provider; it uses VEP 116's genomic `get_3prime_seq_offset` behavior;
- transcript `c.`/`n.` HGVS is per transcript and uses VEP's separate transcript
  `perform_shift` path, not the genomic shift routine; and
- protein `p.` HGVS consumes the alternate peptide difference produced by the same edit
  set used for phased consequence classification.

The hot transcript sweep therefore emits or retains numeric projected-edit facts; it does
not allocate HGVS strings. Rendering is late, after filtering. Before a public phased API
is fixed, the projected-edit sidecar and external VEP-116 HGVS differential must cover
both shift routines, both strands, exon/intron/UTR positions, repeat runs, position-one
right anchors, and compound edits. Apply-then-diff sequence equivalence is an independent
property oracle. Exact structural HGVS can later consume typed exact events. BND HGVS
additionally needs the paired relation's mate and orientation facts; imprecise structural
events remain unsupported until their confidence geometry is represented.

## Supplementary annotations

Exact sources such as dbSNP, ClinVar, and gnomAD should remain sorted numeric-key streams;
range sources should remain interval streams. They share prepared variant keys and final
output, not a universal physical format. DuckDB/Parquet pruning, parallel reads, memory
management, and spill are the first implementation. A custom cache requires a measured
workload that those facilities do not meet.

Variant-level exact and positional joins happen before transcript expansion. The compact
consequence row already carries transcript and gene ordinals, so gene relations are joined
after consequence expansion without repeating a gene lookup inside C. Regulatory, motif,
domain, and other interval sources use range joins or sorted interval streams, followed by
the smallest source-specific predicate needed for a joined pair. Strings and JSON remain
the final projection.

The stateful stable-API integration is tracked at
https://github.com/RGenomicsETL/duckhts/issues/94. Supplementary source plumbing does not
belong as ignored arguments in the consequence-kernel API.

## Validation and performance

- `make test-duckvep-kernel` runs fixed and randomized pure-C properties; ASan and UBSan
  execute the same suite.
- `make test-duckvep-kernel-statistical` raises randomized targets to an explicit seed and
  trial count.
- `make test-duckvep-differential` compares generated witnesses to pinned VEP 116.
- `make duckvep-corpus-differential` records the union of emitted `(variant, transcript)`
  pairs, including mismatches, misses, extras, and unresolved rows.
- The corpus runner's small-event mode samples SNVs, MNVs, insertion-like alleles, and
  deletion-like alleles independently by length-change bin. Structural mode either reads
  a symbolic VCF or generates seeded DEL/DUP/tandem-DUP/INV/CNV/INS events from real model
  geometry at transcript, exon, intron, splice, UTR, CDS, start-codon, and stop-codon
  states. Independent seeds can run concurrently against one read-only attached model.
- A held-out small-event run with seed `20260716` executed 100,000 trials per randomized C
  property (170 tests; 204,521 assertions) and compared 100,268 generated/fixed alleles
  with executable VEP 116. All 100,268 transcript pairs were exact. Its generated set
  contained SNVs, MNVs, insertions, deletions, and delins; duplicate rejection deliberately
  makes the accepted shape counts non-uniform rather than resampling a cosmetically balanced
  result.
- The current eight-seed GRCh38 structural campaign generated 40,375 events on chromosomes
  1, 2, 6, 11, 17, 21, 22, and X and compared 2,140,911 emitted transcript pairs with the
  executable indexed VEP 116 cache. Every pair was resolved and exact, with no missing or
  extra emission. This is evidence for the sampled exact single-locus states, not a claim
  for imprecise coordinates, STR payloads, or untested species.
- The paired-BND campaign generated 1,004 same- and cross-chromosome events on chromosomes
  1, 2, 7, 21, and X and matched all 91,428 isolated executable-VEP transcript pairs.
  BND oracle buffers contain one event because VEP's chromosome-blind mate-coordinate tree
  otherwise makes a record's output depend on neighboring BNDs.
- `make bench-duckvep-release-parquet` reads the official Ensembl variation consequence
  VCF through typed CSQ columns and records complete versus consequence-only Parquet size,
  checksum, cardinality, and elapsed time without committing the large artifacts.
- Throughput reports must state input variants, variant-transcript pairs, output rows and
  bytes, threads, model, and machine. Site-level rate and cohort haplotype-update rate are
  different measurements.
- Comparisons with FastVEP report two timings over the same source and output contract:
  annotation from already coordinate-sorted input, and read-plus-sort-plus-annotation.
  Sorting cost is not hidden, but it is not charged only to DuckVEP either. DuckDB may sort
  once and stream ordered chunks directly into the consequence or phased executor.

The generated reports are the numeric authority. At the latest tested ancestor for each
corpus, the declared GRCh38 dbSNP, GIAB, ClinVar coding, and ClinVar cross-chromosome
samples, GRCh37, and *P. falciparum* are exact against their indexed caches. The two
GRCh38 ClinVar samples contain 604,233 transcript pairs together, with no unresolved,
extra, missing, or discordant row. The separate NMD-plugin differential is exact for all
68,554 eligible transcript pairs, including the 29,416 states both implementations leave
unresolved.

On the final 644,427-transcript GRCh38 model, the current one-core compact measurements are
about 841,000 input alleles/s for the GIAB topology sample, 332,000/s for repeated coding
SNVs, 111,000/s for repeated coding non-SNVs, and 136,000/s for the repeated mixed coding
set. These include stable-API list materialization and aggregation but exclude model load
and input staging. The one-million-input-allele target is not met on the final model; the
mixed lane's roughly 3.92 million emitted transcript rows/s is a second denominator, not a
substitute for the site-rate target. Exact revisions, checksums, row counts, resource
measurements, and conditions live in the generated
[conformance report](../benchmarks/duckvep_conformance.md) and
[throughput report](../benchmarks/duckvep_throughput.md).

Open conformance and throughput work is tracked at
https://github.com/RGenomicsETL/duckhts/issues/93 and
https://github.com/RGenomicsETL/duckhts/issues/95.
