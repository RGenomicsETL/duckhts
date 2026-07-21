# DuckVEP implementation mental model and frontier

Status: current explanatory companion to [`duckvep.md`](duckvep.md), updated
2026-07-20. This file explains the implementation deeply enough to rebuild the
human mental model; it is not a second specification. The concise current contract
is [`duckvep.md`](duckvep.md); the public SQL contract is
[`functions.yaml`](../functions.yaml); code, tests, generated conformance reports,
and generated throughput reports are authoritative. Durable unfinished work lives
in the linked GitHub issues.

The reason for this document follows the distinction made by antirez: using an
agent without reading every line can still be disciplined if the human controls
the main ideas, interfaces, invariants, failure modes, and tests. The dangerous
case is not simply “code written by an agent”; it is code whose concepts are no
longer held by anyone. This note is intended to make those concepts inspectable.

## 1. What DuckVEP is trying to be

DuckVEP is a library-first C consequence engine embedded in DuckHTS. It is not a
new VCF application, a clone of the VEP command-line interface, or a serialized
copy of VEP's Perl objects.

The central bet is a division of work:

- htslib already knows how to transport and decode VCF/BCF and genotypes;
- DuckDB already knows how to scan columnar data, plan joins, prune Parquet,
  parallelize chunks, bound memory, spill, and expose results through SQL;
- a small native C library should own the repeated transcript geometry,
  projection, sequence editing, translation, and consequence predicates;
- HGVS should consume the same projected edit facts as the consequence engine;
- supplementary annotations should be ordinary exact, positional, interval, or
  gene relations, joined at the correct cardinality;
- clinical evidence and ACMG/AMP reasoning should remain explainable relations,
  not branches hidden inside the consequence loop.

The behavioral target is Ensembl VEP 116. VEP's implementation is expressed as
Perl objects and predicate methods, but the semantics can be represented as a
dependency graph:

```text
uploaded allele
    -> interpreted event geometry
    -> candidate transcripts
    -> transcript topology and coordinate projection
    -> sequence/edit facts when needed
    -> predicate-fact bitset
    -> VEP-116 rule program
    -> compact numeric consequence rows
    -> late SQL/HGVS/string projections
```

This is not an attempt to “improve” VEP biology. Compatibility sometimes means
reproducing states that look odd when viewed outside VEP's call order. Those
observations belong in `design/duckvep_errata.md` so they are not rediscovered or
silently cleaned up later.

## 2. Layer responsibilities

### Ingestion and event preparation

The VCF/BCF layer owns record transport, ALT splitting, genotype decoding,
symbolic-allele parsing, and any requested canonical VCF normalization such as
left alignment. The consequence engine does not rewrite the user's variant
identity.

The kernel does need a lossless semantic interpretation of each ALT. That is not
the same as canonical normalization. For small alleles the current event object
keeps three geometries because VEP uses three geometries:

1. `raw_*`: the uploaded VCF REF span, including required padding;
2. `feature_*`: the span VEP gives its `VariationFeature` after its own
   length-changing-allele minimization rules;
3. `start1/end1` plus `insertion_boundary0`: the fully trimmed edit used for
   transcript/CDS projection and sequence application.

It also retains REF/ALT slice offsets, differing lengths, and whether the
validation anchor is on the left or right. An insertion before contig base 1 is
therefore represented as interbase boundary zero with a right anchor; it does not
invent genomic base zero.

This separation matters. Region predicates, splice mismatch islands, and coding
edits do not always consume the same interval in VEP. Collapsing them into one
“normalized span” produced several earlier disagreements.

HGVS normalization is another operation again. Transcript HGVS 3-prime shifting
and genomic HGVS shifting are separate VEP routines with different behavior. They
must consume projected edits; they must not mutate the consequence engine's event
identity or infer edits back from rendered SO terms.

### DuckDB

DuckDB owns the large, mostly one-time or relational work:

- staging pinned Ensembl tables and matching reference sequence;
- building and checking canonical transcript relations;
- sorting/partitioning input variants;
- exact, positional, interval, and gene supplementary joins;
- stable IDs, display fields, source provenance, filtering, and final rendering;
- keeping large auxiliary sources on disk or spilling when memory is insufficient.

### The C kernel

The kernel owns the repeatedly executed mechanics:

- checked event geometry;
- candidate transcript traversal;
- exon/intron/UTR/CDS/directional and splice topology;
- genomic-to-cDNA/CDS projection;
- reference checking when sequence is available;
- edit application and translation;
- structural span facts;
- VEP predicate facts and fact-to-SO evaluation;
- compact numeric result emission.

The public kernel header contains no DuckDB, htslib, Arrow, or Parquet type. The
kernel accepts borrowed arrays and explicit counts. This keeps the same code usable
by DuckDB, a standalone differential runner, sanitizer/property tests, or a Wasm
wrapper.

## 3. What is relational, what is resident, and what streams

The architecture makes three different data-lifetime choices deliberately.

| Data | Representation and lifetime |
| --- | --- |
| Ensembl source tables, stable IDs, attributes, receipts | typed DuckDB relations; retained for provenance and late joins |
| Transcript geometry and sequence needed for repeated consequences | immutable resident C arrays, named and model-local |
| Variants | sorted tiles/chunks; streamed, not retained as a genome-scale resident object |
| dbSNP, ClinVar, gnomAD, scores, regional tracks | DuckDB/Parquet streams or bounded tiles; never copied into every transcript model |
| Per-worker state | mutable workspace owned by one worker and one model at a time |

This is why the Ensembl SQL relations add value even though the hot kernel is C.
They are the compilation and provenance layer: they preserve rich source facts,
allow alternative transcript products, and produce a deterministic compact model.
The hot C model is an execution image, not the sole biological database.

A parser for an Ensembl indexed VEP cache could be added later, but it should be
an alternative producer of the same canonical prepared relations and receipt. It
should not become another consequence authority or another opaque runtime object
graph.

## 4. Building the transcript model

### Inputs

The current builder is a compiler over relations already staged in DuckDB. It is
not a downloader and does not parse Ensembl dump archives itself. Its required
core inputs are:

- `coord_system` and `seq_region` for assembly regions;
- `gene`, `transcript`, `exon`, and `exon_transcript` for topology;
- `translation` for the CDS endpoints inside ranked exons;
- `attrib_type`, `transcript_attrib`, and `translation_attrib` for relevant flags
  and exceptional sequence states;
- a matching tiled FASTA relation `(chrom, start, end, seq)`.

The tables may have been loaded from Ensembl's tab-separated MySQL dumps or read
from an attached read-only MySQL catalog. The extension build remains offline.

### Region compilation

`duckvep_ensembl_regions(...)` checks each FASTA chunk, requires contiguous
coverage from zero for every supplied contig, and matches each contig to exactly
one same-name, same-length Ensembl region on the requested assembly. It then
assigns dense model-local region ordinals. Missing/mismatched regions abort the
query; the builder does not silently borrow a contig from another assembly.

### Transcript compilation

`duckvep_ensembl_transcripts(...)` applies the VEP-116 core-source selection
before it assigns resident ordinals: a transcript must be current, must have a
non-empty stable ID, must not have biotype `artifact`, and must not carry the
`readthrough_tra` attribute. This is behavioral, not cosmetic: loading every
current core row produced 774 transcript pairs that the indexed VEP cache did
not emit in the chromosome-21 differential. For each selected transcript the
builder:

1. joins the gene, translation, ranked exons, and selected attributes;
2. assigns dense transcript and gene ordinals while retaining Ensembl IDs,
   stable IDs, versions, and biotypes in the relation;
3. computes contiguous transcript-oriented cDNA positions for each exon;
4. maps translation start/end from exon-relative coordinates to genomic and cDNA
   coordinates;
5. reconstructs each exon from FASTA, reverse-complements negative-strand exons,
   and concatenates in transcript order;
6. extracts the CDS and applies Ensembl start-phase padding (`N` for a positive
   start phase);
7. extracts the complete spliced transcript sequence before and after the CDS;
8. reads the sequence region's `codon_table` attribute, defaulting to table 1 as
   VEP does, and selects the corresponding immutable 64-codon translation table;
9. parses supported single-position `initial_met`, `_selenocysteine`,
   `amino_acid_sub`, and `_stop_codon_rt` Translation SeqEdits into a sparse
   reference-peptide overlay;
10. projects each repeated Ensembl `miRNA` cDNA range through the ranked exons
   into one or more inclusive genomic segments; and
11. records compact flags for translation, protein-coding/NMD/miRNA biotype,
   incomplete CDS ends, selenocysteine, stop readthrough, RNA/peptide edits,
   MANE, GENCODE, CCDS, and upstream-start state. The standalone ABI retains a
   readthrough flag for deliberately prepared alternate models, but the
   VEP-compatible core importer filters those transcripts.

Supported single-amino-acid Translation SeqEdits do not rewrite the prepared CDS.
They patch VEP's reference peptide at the affected protein positions; the alternate
peptide remains the raw translation of the edited CDS. This distinction is required
for terminal stop/readthrough and selenocysteine predicates. Transcript RNA edits
and arbitrary range or length-changing Translation SeqEdits still withhold sequence
with an explicit reason. Coordinates and flags remain available. Malformed topology
is different: invalid coordinates, exon order, phase, translation bounds, or
incomplete sequence reconstruction abort the model build.

### GRCh38, GRCh37, GENCODE, and MANE

GRCh38 and GRCh37 are separate source models. GRCh37 is not a coordinate switch
applied to the GRCh38 gene build. Release 116 exposes an archived GRCh37 core
around the frozen release-75/GENCODE-19 annotation. It has no MANE mapping because
MANE is a GRCh38 product; DuckVEP must not fabricate one.

On GRCh38, the prepared cold relation retains versioned RefSeq accessions from
MANE Select and MANE Plus Clinical attributes. Only compact selection flags enter
the resident model. Ensembl/GENCODE/RefSeq identifiers stay relational for late
projection, transcript selection, and future RefSeq HGVS.

### Species and genetic codes

The consequence registry is VEP code, not a human vocabulary: VEP 116 defines 41
terms for every species. Species changes the model facts that feed those rules.
The kernel now accepts every NCBI codon table present in VEP 116's BioPerl
translator. The importer resolves `seq_region_attrib` through
`attrib_type.code = 'codon_table'`, defaults to table 1 when the attribute is
absent exactly as VEP does, records observed tables in the receipt, and fails before
publication for conflicting, malformed, or unsupported values. Translation remains
one indexed 64-codon lookup; there is no species-name or chromosome-name branch.

Complete model builds and indexed-cache differentials currently cover:

| Source | Resident model | Observed codon tables | Latest consequence result |
| --- | ---: | --- | ---: |
| human GRCh38, Ensembl 116 | 644,427 transcripts | 1, 2 | retained dbSNP, GIAB, and both ClinVar corpora exact |
| human GRCh37, Ensembl 116 archive | 195,379 transcripts | 1, 2 | 486,464 / 486,464 exact transcript pairs |
| *P. falciparum*, Ensembl Genomes 63 with VEP 116 libraries | 5,791 transcripts | 1, 4, 11 | 40,732 / 40,732 exact transcript pairs |

Those rows prove only the named assemblies, source releases, tables, and sampled
distributions. Mouse, Drosophila, Arabidopsis, and species using other genetic codes
remain deliberate additions to the acceptance matrix. Tiny source-faithful excerpts
belong in offline package tests; complete caches remain external differential inputs.

Tracked at https://github.com/RGenomicsETL/duckhts/issues/119.

### Receipt and publication

`duckvep_model_receipt(...)` checks dense ordinals and region/transcript agreement,
then records declared source/release/assembly/filter, source and reference hashes,
model counts, and a deterministic hash over every semantic model field, including
projected mature-miRNA ranges, reference-peptide edits, and the optional core
regulatory/motif interval relation. The latter records separate RegulatoryFeature
and MotifFeature counts and hashes the five resident fields after VEP's source
filter has excluded `epigenetically_modified_region` rows. It has no timestamp:
identical inputs must rebuild the same receipt.

`duckvep_model_load(...)` reads three required sorted projections—regions,
transcripts, and exon memberships—through a private DuckDB connection. One optional
query supplies `(transcript_index, mature_mirna_start, mature_mirna_end)`; another
supplies `(transcript_index, protein_position, alternate_amino_acid)` for supported
reference-peptide edits. A third optional query supplies dense regulation-feature
ordinals, model-region ordinals, inclusive start/end coordinates, and a compact
kind distinguishing RegulatoryFeature from MotifFeature. The loader validates and
narrows every value, packs the sparse transcript relations and the independent
interval-feature SoA, builds separate transcript and feature cgranges seed indexes,
opens the pure-C model, and publishes the name only after the complete operation
succeeds. Failed construction publishes nothing.

The transcript projection supports old 11- and 12-column shapes for compatibility,
plus the current 13-column shape containing complete pre-CDS and post-CDS sequence.
Old models remain usable but return `missing_transcript_flank` for predicates that
need pre-CDS or post-CDS transcript sequence absent from those projections.

A model is partial unless explicitly loaded with
`transcript_coverage_complete := TRUE`. Only a complete model may turn “no loaded
transcript candidate” into supported `intergenic_variant`. A partial model returns
an unresolved `no_feature_in_loaded_model` row.

### Fail-closed model invariants

The kernel validates the execution model again at `duckvep_model_open`. Important
invariants include:

- valid nonzero transcript spans and strand exactly +1 or -1;
- transcripts sorted by `(chrom_id, start1)`;
- each transcript's exon slice lies inside the exon arrays;
- exons are non-overlapping, in transcript order, and their genomic envelope is
  exactly the transcript span;
- exon cDNA spans start at 1, preserve exon lengths, and are contiguous;
- phase and end phase are in `-1,0,1,2`;
- CDS boundaries project into exons;
- prepared CDS length agrees with projected cDNA length and start phase;
- codon table is supported;
- CDS/flank offsets are inside their byte pools and bases are A/C/G/T/N;
- complete flank lengths agree with the cDNA projection; and
- mature-miRNA offsets cover the complete side relation monotonically, while
  each segment belongs to a miRNA transcript and lies wholly inside one exon; and
- regulatory/motif ordinals are dense, model regions exist, coordinates are valid
  and sorted, and each compact feature kind is recognized.

These checks happen once. The hot loop can then trust array bounds and topology.

## 5. Resident memory layout and ownership

The resident model is structure-of-arrays (SoA), not one heap object per
transcript. Frequently accessed columns—contig, start/end, strand, flags, exon
slice, CDS endpoints—are independent contiguous arrays. Exons are another set of
parallel arrays. CDS bytes and transcript flanks live in packed byte pools with
per-transcript offset/length arrays. There is no per-transcript sequence allocation
or padding.

The raw array cost can be reasoned about without pretending that array payload is
the same as measured process RSS. Current transcript columns occupy roughly this
many payload bytes:

```text
per transcript = region(2) + span(8) + strand(1) + flags(8)
               + gene(4) + exon slice(6) + CDS genomic bounds(8)
               + CDS offset/length/table(13)
               + pre-flank offset/length(12)
               + post-flank offset/length(12)
               + point-path eligibility(1)

per exon membership = genomic span(8) + cDNA span(8) + phase/end-phase(2)

sequence bytes = prepared CDS bytes + complete spliced pre/post-CDS bytes
index bytes = cgranges index storage
mature-miRNA bytes = 4 * (transcript_count + 1) + 8 * segment_count
peptide-edit bytes = 4 * (transcript_count + 1) + 5 * edit_count
validated-CDS cache = 13 * transcript_count
regulation/motif payload = 11 * interval_feature_count, plus its cgranges index
```

The final GRCh38 model has 2,806 mature-miRNA genomic segments, 389 peptide
edits, 380,818 admitted regulatory regions, and 1,002,762 admitted motif features.
In both sparse transcript relations, the per-transcript offset array costs more
than the sparse payload. That fixed cost buys allocation-free lookup and must be budgeted
for every coexisting named model. The interval-feature payload stays independent
of transcript rows and cold funcgen metadata remains relational. The validated-CDS cache stores three `uint32_t`
coordinates and one validity byte per transcript; its exact element payload on the
644,427-transcript model is 8,377,551 bytes (7.989 MiB).

Because the columns are separate allocations, there is no AoS padding per row;
the honest total still needs measurement on the full model, including allocator,
index, registry, DuckDB relation, and workspace memory.

The registry belongs to one DuckDB database. Named models are immutable and may
coexist in the same process. Ordinals are only meaningful inside one model.
Annotation pins a model; dropping it fails while pinned.

Each concurrent annotation callback borrows one exclusive workspace for its lifetime:

- active transcript indexes and a span-candidate buffer, each sized to the
  largest transcript run on a contig;
- separate active and candidate feature indexes when the model carries core
  regulatory/motif intervals;
- a 16-bit current exon rank per transcript for the monotone point path;
- edit, alternate-CDS, and peptide scratch sized from the model's maximum CDS and
  bounded allele capacity;
- caller-owned result storage.

No mutable FASTA handle, htslib iterator, active set, sequence cache, or result
builder is shared between concurrent callbacks. Returned workspaces may be reused by a
later callback. This is the core thread-safety rule and also allows different models to
execute in the same process.

Independent-event parallelism uses disjoint, internally ordered input branches. Each
concurrent scalar callback borrows one exclusive workspace from the immutable model's
pool; a partition is not permanently bound to a workspace, and the current scalar sweep
restarts at DuckDB vector edges. Stable active-set compaction preserves transcript and
feature order inside each variant list, so the per-input `annotation_index` does not
depend on the vector or ordered branch that processed it. A final global row order is
still an SQL property: consumers that need one use
`ORDER BY input_variant_index, annotation_index` after the parallel union. The scalar
adapter does not manufacture these branches merely because `SET threads` is greater than
one; the caller or a future table-function planner must expose them explicitly.

The measured final GRCh38 artifact is 1.608 GiB. In the pinned benchmark process,
loading it took 3.047 seconds; RSS rose from 0.114 GiB to 4.125 GiB and peaked at
5.041 GiB while DuckDB held its source-table cache and the immutable C image at the
same time. Planned scratch is 3.208 MiB per worker. These are process measurements,
not a claim that the C arrays alone consume four GiB; the generated throughput
report keeps the exact conditions and accounting.

## 6. Independent small-variant execution

The stable SQL surfaces today are:

```text
duckvep_annotate(model_name, seq_region, position, reference, alternate
                 [, distance] | [, upstream_distance, downstream_distance])
duckvep_annotate_compact(model_name, seq_region, position, reference, alternate
                         [, distance] | [, upstream_distance, downstream_distance])
```

It handles independent biallelic A/C/G/T/N small variants. Input should be sorted
by model-local contig and coordinate to earn the sweep behavior. Both directional
windows default to 5,000 bases as in VEP; one optional value changes both, while two
values configure upstream and downstream separately. Zero disables the corresponding
side; a single symmetric zero disables both.

### Step 1: adapt one DuckDB vector

The scalar adapter copies a DuckDB vector into compact arrays and one allele byte
pool. It validates non-null fields, coordinate widths, allele alphabet, and event
shape. Work is split when model, contig, distance configuration, or sort order changes.

The adapter is vectorized in the DuckDB sense, but its expression-local state does
not survive arbitrary future vectors. It therefore seeds/restarts a sweep at each
DuckDB vector. This is already a real batch path, but it is not yet a whole-file
stream with explicit carry state.

Worker workspaces can be returned to the pool and reused. Point and normalized-span
exon ranks therefore have an explicit cross-vector invariant: after a vector whose
prepared coordinates are non-monotone, both rank arrays are reset before the next
non-empty vector. Comparing only the next vector's first coordinate with the prior
vector's final coordinate is insufficient because a transcript skipped after an earlier
forward jump can retain a rank beyond its next admitted event.

### Step 2: prepare event geometry once

Opening the annotation cursor validates the whole batch and constructs one compact
`duckvep_event_t` per variant. REF/ALT prefix/suffix comparison and right-anchor
handling are therefore not repeated for every overlapping transcript.

The kernel does not left-align or rewrite alleles. It only derives the raw,
VEP-feature, and semantic-edit views described earlier.

### Step 3: find candidate transcripts using sortedness

The adapter queries the cgranges index for the first event's point window. This is
the seed. From there the pure-C cursor advances through transcripts and variants
monotonically.

For each event at position `p`, the cursor maintains transcripts whose start is at
or before `p + halo` and whose end plus halo has not fallen before `p`. As `p`
increases, new transcripts are admitted once. Expired transcripts are removed and
cannot become candidates again.

Span ends are not necessarily monotone. A large deletion/CNV could otherwise keep
many future transcripts in the hot active set and poison later SNVs. The current
sweep therefore keeps a persistent point frontier and appends the extra transcript
tail required by the current span only into a separate candidate buffer. That tail
is discarded after the event.

For `V` sorted variants, `T` transcripts admitted by the frontier, and `P` emitted
candidate pairs, the intended traversal cost is approximately `O(V + T + P)` plus
one tile seed, rather than `O(V log T)` independent searches. Biological work is
still proportional to candidate pairs; a dense gene region genuinely creates more
work and output.

The active array uses compact transcript ordinals. Expiry currently removes an
element by replacing it with the last active entry, so active iteration is not a
stable transcript-order walk. The SQL adapter enforces result grouping/order before
materialization. Any candidate-index experiment must compare actual candidate
pairs, order costs, and cache misses, not only asymptotic lookup complexity.

### Step 4: topology and splice classification

Each candidate pair first takes the cheap path. The classifier determines whether
the feature span is upstream, downstream, intronic, exonic, UTR, CDS, complete
feature overlap, partial overlap, or wholly within the feature. Span placement is
not mutually exclusive.

Splice classification is separate because a VEP insertion can have a reversed
interbase interval and a length-changing feature allele can contain distinct
mismatch islands. Donor, acceptor, donor fifth base, donor region,
polypyrimidine tract, and generic splice-region facts are recorded independently.

Sorted SNVs have a specialized traversal: each transcript has a current exon rank
that moves only forward as positions increase. The specialization calls the same
topology/splice predicates as the exhaustive path and is property-tested against
it. It is an optimization, not a second biological implementation.

### Step 5: project and evaluate sequence only when needed

Noncoding/topology-only pairs do not reconstruct or translate sequence. A pair in
the CDS bucket is projected to cDNA/CDS coordinates and escalates to the sequence
delta layer.

The generalized delta path owns MNVs and indels. It can represent several
differing islands, constructs one edit/CDS/peptide context, validates reference
bases, applies the edit in transcript orientation, translates the necessary
sequence, and derives facts such as:

- synonymous or missense;
- start lost/retained;
- stop gained/lost/retained;
- frameshift or restored frame;
- in-frame insertion/deletion;
- general protein alteration or unresolved coding state.

Production no longer retries an unresolved generalized context through an older
shape-specific classifier. Projection failure, missing sequence/flanks, ambiguous
bases, REF mismatch, unsupported compound state, or capacity failure stays
reason-coded and visible. Narrow direct classifiers remain only as pure-C test
oracles.

### Step 6: one fact-to-SO program

Topology, splice, event-shape, structural, and optional sequence facts are packed
into predicate bits. A generated static rule table contains required facts,
forbidden facts, output SO mask, tier, rank, and impact. Tier controls suppression;
rank is severity metadata. This separates “how a fact is computed” from “which VEP
term that fact combination emits”.

An exhausted predicate list for a real transcript overlap receives VEP's
transcript-associated default `intergenic_variant`, even though that looks odd.
No transcript candidate at all is handled separately by model-completeness policy.

### Step 7: NMD as a separate result

`NMD_transcript_variant` is a consequence of the imported transcript already
having the `nonsense_mediated_decay` biotype. It does not claim the current allele
creates NMD.

Variant-induced NMD is an additional compact result based on the pinned VEP
Plugins release/116 `NMD.pm` positional rules. Eligible stop-gained, frameshift,
splice-donor, or splice-acceptor events are predicted escaping for an intronless
transcript, an early coding event, an event in the last exon, or the plugin's
penultimate-exon-end window. Otherwise they are `triggering`; missing coding
projection is `unresolved`. Escape reasons are independent bits.

The plugin uses the complete VEP `VariationFeature`, not the minimized edit used
to change sequence. DuckVEP therefore projects both endpoints of the full feature
for NMD. A pure insertion keeps VEP's reversed parent CDS interval, for example
`102,101`; immediately before the first coding base, `1,0` is valid because the
plugin tests definedness rather than truth. A wider padded allele and its minimal
equivalent can consequently change the same bases but sit on different sides of an
NMD position test. The coding classifier's cached early-CDS fact is reused only
when the full feature and minimized edit occupy exactly the same genomic span.
Executable-plugin witnesses on both strands pin these representation-dependent
rules.

### Step 8: bounded output and late strings

The C kernel appends fixed numeric records to caller-owned storage. A resumable
cursor preserves the current event, candidate slice, and candidate position when
the result buffer fills. The caller drains/resets the buffer and resumes without
recomputing the tile or silently truncating rows.

The rich SQL adapter materializes a list of structs and renders consequence,
impact, region, status, reason, amino-acid, and NMD strings.
`duckvep_annotate_compact(...)` runs the same cursor but exposes transcript/gene
ordinals, SO and region masks, stable numeric status/reason/NMD codes, coordinates,
and amino-acid bytes. It is the filtering/joining floor; strings and stable IDs can
be projected later. The return is still a DuckDB list, so it does not settle the
future appender/direct-vector materialization comparison.

## 7. Cache behavior, vectorization, and the throughput target

The target is not “fast for a toy transcript”. It is at least one million sorted
input ALT alleles per second on one pinned core against the declared full Ensembl
model while emitting every compact transcript consequence. Candidate pairs,
output rows, output bytes, and transcript density must accompany the input rate.

The one-transcript smoke after PR #116 reached more than three million inputs/s.
That only shows the adapter loop is not intrinsically slow. The nearest historical
one-core baselines use the final 644,427-transcript, 5,068,416-exon GRCh38 model on one
pinned i5-13500 core. The generated
[throughput report](../benchmarks/duckvep_throughput.md) includes the stable DuckDB
adapter, compact list materialization, `unnest`, aggregation, and an explicit input
`ORDER BY`; model loading and input staging are outside the timed pass. These exact
workloads predate source `e25c151`, whose directly comparable current-revision evidence
is the annotation-dense matrix below.

| Workload | Input alleles | Output rows | Median seconds | Input alleles/s | Output rows/s |
| --- | ---: | ---: | ---: | ---: | ---: |
| GIAB 1-in-40, transcripts only | 100,957 | 1,174,245 | 0.120 | 841,308 | 9,785,375 |
| same GIAB stream plus 1,383,580 regulatory/motif intervals | 100,957 | 1,179,329 | 0.126 | 801,246 | 9,359,754 |
| repeated coding SNVs | 200,000 | 5,191,000 | 0.603 | 331,675 | 8,608,624 |
| repeated coding non-SNVs | 36,000 | 1,047,524 | 0.322 | 111,801 | 3,253,180 |
| repeated mixed coding alleles | 200,000 | 5,756,720 | 1.481 | 135,044 | 3,887,049 |

The coding sets deliberately repeat a frozen allele sample to stabilize resident
kernel timing, so they favour warm model and sequence caches. They are not file-
ingestion or whole-genome timings. The rows also expose the difference between site
rate and expanded work: the mixed lane writes almost 29 transcript rows per allele
and still exceeds 3.8 million output rows/s, but 135,044 input alleles/s is far below
the declared one-million input floor.

Those historical rows do not establish current-head regression status. They do show
that the floor remained open on the final model and that the coding-heavy gap was large;
the identical workloads still need rerunning after a future hot-path change when a direct
comparison is required. Supplementary joins and final rendering are unmeasured end to
end. Every optimization must name and measure the work it removes—candidate visits,
model gathers, projection, sequence bytes, translation, row writes, or DuckDB
materialization—rather than substituting a smaller fixture or a different denominator.

At realistic density, the likely costs are:

- model gathers for every candidate transcript;
- exon walking and branch-heavy exon/CDS junction cases;
- sequence/CDS bytes for the minority of coding rows;
- writing many consequence rows per input variant;
- building DuckDB nested lists and assigning several strings per row;
- downstream `UNNEST` and joins.

DRAM or SSD bandwidth is the eventual floor, but current code can hit instruction,
branch, cache-miss, and materialization limits first. A performance report must
separate:

1. event preparation (alleles/s and allele bytes/s);
2. candidate discovery (candidate pairs/s and cache misses);
3. topology (pairs/s);
4. coding sequence work (coding pairs/s and edited/translated bases/s);
5. compact C emission (rows/s and bytes/s);
6. stable DuckDB API materialization (rows/s and bytes/s);
7. supplementary scans/merges;
8. final string/CSQ/JSON rendering.

SoA helps cache locality and gives the compiler simple contiguous loops, but this
engine is not uniformly SIMD-shaped. Exon projection and variable edit/translation
work are gather- and branch-heavy. Likely SIMD/autovectorization targets are:

- coarse span comparisons over start/end arrays;
- bulk base validation;
- compact-key comparisons for exact annotations;
- translation of many independent codons/haplotype sequences;
- consequence-mask-to-impact reductions;
- count-only interval queries.

The common sorted SNV path is more likely to win first from smaller hot columns,
predictable branches, inlining, and keeping active state close than from a wide
handwritten SIMD loop.

The current transcript model could be split even more explicitly into a tiny hot
header and cold sequence/flags metadata if counters show repeated cache-line waste.
Likewise, the current result AoS contains fields that can be derived later. Those
changes need `perf`/hardware-counter evidence on the full model.

## 8. Candidate-index research still open

The production path currently uses cgranges once to seed a tile and then a forward
active sweep. That is already different from a binary search for every variant.

An implicit interval tree, SuperIntervals-style backward branch array, or another
flat superset index may still be better for:

- the first seed;
- unsorted interactive queries;
- supplementary interval sources;
- perhaps the whole candidate lane in gene-dense or deeply nested transcript
  regions.

For sorted transcript annotation, an alternative must beat the continuing sweep,
not a straw-man binary search. The benchmark needs candidate pairs/s, false
candidate visits, cycles/pair, L1/L2/LLC misses, branch misses, index bytes, and
point/span workloads across sparse and dense regions. A monotone SuperIntervals
cursor is a legitimate experiment; it is not yet an architectural decision.

## 9. Statistical conformance as anti-reward-hacking

There are several different tests because they catch different failures.

### Fixed witnesses

Minimized examples anchor known VEP states and regressions. They are readable and
fast but can be overfit.

### Randomized pure-C properties

The property suite checks mechanical invariants and equivalence between optimized
and exhaustive paths: projection, reverse-strand behavior, sweep candidates,
splice predicates, edit application, interaction partitioning, and cursor output
splits. ASan and UBSan execute the same library paths.

### Statistical VEP differential

Fresh corpora are generated from explicit seeds and distributions, executed by the
real pinned VEP 116, and compared at the union of emitted `(variant, transcript)`
pairs. Missing rows, extra rows, unresolved rows, and full consequence-set
disagreements all remain in the denominator.

This is an explicit anti-reward-hacking mechanism. The implementation cannot earn
agreement merely by adding branches for a fixed list of known witnesses. The
distribution should itself be versioned and expanded across allele shapes,
lengths, exon/CDS/transcript endpoints, strands, transcript states, and rare
predicate combinations.

PR #116 recorded:

- fixed differential: 268/268 exact;
- four independent 100,000-case corpora: 401,072/401,072 exact in aggregate;
- 3,900,500 randomized property trials with zero failure;
- clean ASan/UBSan, SQL, and R package suites.

Generated exactness is strong evidence for the sampled state space, not proof of
whole-human correctness. The real-data lane now keeps the same union denominator
and publishes full-set, per-SO-term, impact, allele-shape, strand, transcript-state,
and unresolved strata. At the latest tested ancestor for each corpus:

| Corpus | Exact annotation-object consequence sets | Unresolved | Resolved disagreements |
| --- | ---: | ---: | ---: |
| GRCh38 dbSNP | 73,620 / 73,620 | 0 | 0 |
| GRCh38 GIAB | 54,905 / 54,905 | 0 | 0 |
| GRCh38 ClinVar coding | 287,836 / 287,836 | 0 | 0 |
| GRCh38 ClinVar across chromosomes | 316,397 / 316,397 | 0 | 0 |
| GRCh37 | 486,464 / 486,464 | 0 | 0 |
| *P. falciparum* | 40,732 / 40,732 | 0 | 0 |
| GRCh38 paired BND | 91,428 / 91,428 | 0 | 0 |
| GRCh38 GIAB plus core regulation/motif | 14,955 / 14,955 object pairs | 0 | 0 |
| GRCh38 exact SV plus core regulation/motif | 120,224 / 120,224 object pairs | 0 | 0 |
| HPRC v2 GRCh38 chromosome-22 long literal alleles | 4,670 / 4,670 transcript pairs | 0 | 0 |

The HPRC row is a deterministic carried-ALT sample from four African-ancestry
pangenome individuals. It includes uploaded REF/ALT strings longer than VEP's
5,000-base transcript-distance default and closes ordinary complete-transcript
deletion and equal-length-span states. The distance controls candidate admission
outside transcript endpoints; it does not clip an overlapping feature span or
impose an allele-length limit.

The GRCh37 and *P. falciparum* rows closed their declared indexed-cache samples
only after codon-table sourcing, supported Translation SeqEdits, terminal peptide-
edit behavior, and assembly-specific source differences were fixed. They do not
stand in for untested species. The two ClinVar rows closed their retained samples
only after every discrepancy was minimized and classified; fresh independent-
small-variant distributions remain necessary to explore states outside those frozen
corpora. The paired-BND row covers 1,004 generated same- and cross-chromosome events,
all four bracket orientations, and 14 observed SO terms.

At the paired-BND implementation revision, the pure-C property suite contained 43
randomized targets, 4,200,500 trials, and 205,610 assertions with no failure.

Variant-induced NMD has its own executable differential against VEP Plugins
release/116 `NMD.pm`. On the ClinVar chromosome-21 corpus, the core SO result is
exact for 1,331,664 transcript pairs. Among 68,554 NMD-eligible pairs, DuckVEP
matches 6,937 escaping, 32,201 triggering, and 29,416 mutually unresolved states,
with zero missing, extra, or mismatched predictions.

The generated [conformance report](../benchmarks/duckvep_conformance.md) is built from
the append-only CSV ledgers and is the current numeric authority. Exact sampled
corpora remain evidence for their declared distributions, not a license to stop
randomized state exploration.
Phased interactions, more species, raw repeat expansion, producer-specific structural
encodings, and real fusion corpora require their own generated distributions.
Those are not all missing consequence predicates: VEP 116 evaluates structural events
with confidence metadata and inserted-sequence payloads at nominal `POS`/`END`, while
bounded repeat expansion is an input-preparation operation that can turn a
`<CNV:TR>` record into an ordinary literal allele. Core regulatory/motif coverage
now has exact GIAB-small-variant and generated-SV distributions, but other funcgen
releases, species, paired BNDs, and uncertain event representations remain separate
evidence. Independent small-variant exploration also continues; exact frozen
corpora are checkpoints, not a stopping rule.

For the declared VEP-116 consequence contract, the consequence engine is now treated as
proven infrastructure rather than an open predicate-development project. Concretely,
that statement covers the registered 41-term rule program, supported independent literal
small variants, typed exact single-locus structural events, paired breakends, core
RegulatoryFeature/MotifFeature objects, the declared assemblies/species/model receipts,
and the recorded generated and real-data distributions above. It does not claim phased
or compound consequences, raw bounded `<CNV:TR>` expansion, producer-specific
symbolic/BND parsing, genomic HGVS, supplementary databases, or a species/model not named
by a receipt. `CIPOS`/`CIEND` remain preserved uncertainty metadata while strict VEP-116
consequences use nominal `POS`/`END`; no additional confidence-aware consequence policy is
being claimed. Randomized properties and executable-VEP differentials remain permanent
regression gates; new distributions extend the evidence without making day-to-day work
wait on an unbounded claim of exhaustive biology.

## 10. Phased haplotypes and Haplosaurus

Independent variant annotation cannot recover the effect of variants in cis. Two
SNVs can change one codon; two indels can open and later restore a reading frame;
one edit can change the reference seen by the next. Haplotype annotation therefore
uses the same projection/edit core but a different grouping and output unit.

### What exists

The pure-C layer has:

- a projected CDS edit element with original-CDS coordinate, REF/ALT, and variant
  strand;
- interaction partitioning for edits sorted on the original CDS axis;
- grouping rules that keep a block open while cumulative length change displaces
  the reading frame or the next edit touches the same alternate codon;
- a linear reverse-coordinate CDS rebuild, REF validation, strand orientation,
  and one translation of the resulting CDS;
- flags for indel, open frameshift, restored frameshift, and stop truncation.

The normal distinct-buffer path validates the edit set first and then copies each
unchanged CDS segment and alternate base exactly once, so mutation is
`O(CDS bytes + ALT bytes + edits)`. Exact input/output aliasing remains compatible
through the older `memmove` path because cumulative displacement can change sign
within one haplotype; the phased executor owns distinct worker scratch and takes the
linear path.

### What is missing

The streaming executor must decode GT/PS and group projected edits by:

```text
(model, transcript, sample, phase_set, haplotype_lane)
```

It should keep reference paths implicit and allocate state only for non-reference
carriers. With coordinate-sorted input, transcript state can be finalized once the
stream watermark passes the transcript end. Identical sparse edit prefixes across
samples should be hash-consed or represented as a prefix DAG/radix tree so each
distinct final CDS/protein is built and translated once.

The result must retain all contributing variant IDs and carrier counts. It needs
two explicit phase policies:

- `vep116_compat`: reproduce Haplosaurus allele-slot behavior;
- strict phase: require real phase information and respect phase-set boundaries,
  representing unphased ambiguity explicitly.

Haplotypes are not a second consequence engine. At the coding layer they are an
N-edit set passed to the same fact production and rule table. The separate concerns
are grouping, lifetime, deduplication, attribution, and output.

Long-read inputs add contracts outside the edit loop. A reproducible callset needs
its own receipt naming read chemistry, basecaller/model, reference and aligner,
small-variant and SV callers/options, phasing method, and callable/confidence masks.
`GT` plus a shared phase-set value is not proof that overlapping records form a
valid alternate sequence. The executor must reconstruct one consistent local
haplotype or return an edit-conflict result while preserving every source row.
An HPRC assembly is orthogonal truth evidence, not automatically an execution
model: path-coordinate annotation also needs transcript projection, mapping
confidence, and locus-level collapse/misassembly QC. The tracked R3 inputs and QC
are at https://github.com/human-pangenomics/R3_HPRC. Cross-caller reconciliation
must also distinguish two descriptions of one event from two edits that both
belong on the haplotype; sample, phase set, and nearby coordinates are insufficient.
If a phase block or transcript projection changes, incremental annotation
recomputes affected transcript haplotypes rather than treating only the new VCF
row as a safe delta.

Tracked at https://github.com/RGenomicsETL/duckhts/issues/92.

## 11. Structural variants

The engine does not force every event into a single `[start,end]` plus ALT string.
The typed public families now cover exact single-locus DEL, DUP, tandem-DUP,
tandem-repeat (`STR`), INV, INS, CNV, and unknown operations plus paired BND.
Single-locus calls carry an explicit copy-change direction. The adapter rejects
contradictory combinations and the kernel derives feature ablation, amplification,
elongation, truncation, copy-gain/loss, and breakpoint facts from span topology.
Structural `STR` remains a distinct event for provenance and later HGVS but uses
VEP 116's tandem-repeat predicate algebra: copy gain, insertion, and amplification.
When VEP can expand a bounded `<CNV:TR>` from repeat metadata below its structural
size limit, that literal allele instead follows the ordinary small-variant path.
The structural differential explicitly passes and receipts
`VEP --max_sv_size 10000000` so large exact spans in its distribution are not silently
skipped. VEP's own 5,000-base default is also the bounded-repeat expansion limit;
neither value is the transcript-distance control or the separate fixed BND admission
distance.

`duckvep_annotate_sv(...)` and its compact form execute single-locus events through
the sorted transcript sweep. `duckvep_annotate_breakend(...)` and its compact form
accept both region/position pairs in one row, query candidate transcripts and
resident RegulatoryFeature/MotifFeature objects around both loci, merge and
deduplicate each object class, and evaluate the pair as one semantic event. This
preserves VEP 116's asymmetric transcript behavior: ordinary topology comes from
the shifted local locus, while a mate-only transcript can contribute through
`feature_truncation`. RegulatoryFeature and MotifFeature semantics are also
endpoint-asymmetric: a shifted-local hit keeps its ordinary base term, a mate-only
exact hit uses generic `feature_truncation`, and a same-contig shifted local point
outside but within VEP's fixed 5000-base admission distance adds
`intergenic_variant` to that mate-discovered object. The local base term wins if
both points hit one object exactly. VEP may duplicate identical rows or emit
distinct allele-level consequence rows for one object; DuckVEP returns their
union once. Neither one wide span nor
two independent point calls has those semantics.

The paired-BND differential generated 1,004 events on chromosomes 1, 2, 7, 21,
and X, including all four bracket orientations and same- and cross-chromosome
pairs. DuckVEP matched all 91,428 VEP transcript consequence sets with no missing,
extra, unresolved, or discordant row. VEP was run with one BND per semantic buffer
because its otherwise chromosome-blind mate-coordinate interval tree lets a nearby
record change the candidate set. The checked-in harness records that oracle rule.

The surrounding preparation work must retain inserted sequence where known,
imprecise start/end confidence intervals, STR repeat unit and reference/alternate
counts, raw BND ALT/orientation, mate identity, and source provenance. Those facts
belong in typed input relations; the consequence loop must not repeatedly infer
them from strings. They are not unfinished VEP-116 SO predicates: the executable
uses nominal structural coordinates, and its structural insertion predicate does
not inspect the inserted payload. They remain necessary for lossless round trip,
HGVS, uncertainty-aware interpretation, and fusion work. Phased SVs and small edits
must feed the same haplotype grouping model rather than creating a second
compound-event engine.

The executable regression for this rule is
`test/duckvep/conformance/data/structural_confidence_grch38.vcf`. Its 12 records pair
nominal and `IMPRECISE;CIPOS;CIEND` forms of CNV, DEL, DUP, tandem DUP, INV, and INS.
The pinned VEP 116 cache and DuckVEP matched all 466 transcript pairs; within each
engine, all six nominal/imprecise consequence multisets were identical. The
conformance runner retains the declared source ID and confidence fields when it writes
the sampled oracle VCF, so this test cannot accidentally erase the state it claims
to exercise.

The append-only public SO mask now binds all 41 terms in VEP 116's
`%OVERLAP_CONSEQUENCES` registry. Six regulatory/TFBS terms have a separate
interval-feature evaluator over the same event geometry. The release-matched
Ensembl funcgen importer, deterministic receipt, resident feature SoA/index, and
rich/compact emission are implemented for small alleles, exact single-locus
structural events, and paired BNDs. BND feature candidates use both endpoint
points, then enter the same explicit event/object-pair validator as transcript
candidates. The feature rows retain their own ordinal; they are not fake
transcript rows.
The forty-first term, `sequence_variant`, has no VEP overlap predicate and remains
registry metadata rather than a fallback used to conceal missing facts. FastVEP's
49-entry enum is not the VEP contract: it appends eight Nirvana-parity labels after
the 41 VEP entries.

Tracked at https://github.com/RGenomicsETL/duckhts/issues/98.

## 12. Supplementary annotations

The transcript model is small enough to be resident. Dozens of resources such as
dbSNP, ClinVar, gnomAD frequencies, prediction scores, constraint tables, and
regional tracks are not. They should not be copied into every named model or
hidden behind ignored kernel arguments.

The cardinality order should be:

```text
VCF record
  -> one prepared row per ALT
  -> exact and positional annotation while cardinality is still one row/ALT
  -> transcript consequence expansion
  -> transcript ordinal -> gene ordinal
  -> gene annotation once per unique gene
  -> filtering/picking
  -> stable IDs, SO strings, HGVS, CSQ/JSON rendering
```

### Exact sources

Use compact numeric keys and sorted merge when both streams are sorted. Long
alleles can use offsets into a side pool. Source-specific canonical keys (for
example a left-aligned gnomAD key) must not mutate the uploaded event used for
consequence prediction.

For sparse access, partition Parquet by contig/coarse genomic tile. Zone maps and
predicate pushdown should remove untouched row groups; merge within selected
tiles. String chromosome joins are unnecessary: model-local numeric regions or a
declared chromosome dictionary can be used.

### Positional and interval sources

Positional sources are another sorted merge. Interval sources use a continuing
active sweep for sorted streams or a measured flat/tile-local index for sparse
random access. Large interval tracks remain DuckDB/Parquet data, with memory
bounded by active windows/tiles.

### Gene sources

Map transcript ordinals to gene ordinals, collect the unique genes for a batch,
and join gene-level resources once. Repeating the same gene strings for every
transcript consequence before filtering is avoidable work.

### Why DuckDB rather than a new physical format

FastVEP/fastSA and echtvar-like stores demonstrate useful ideas—numeric keys,
binning, compact categorical columns, late strings, batch preload. DuckHTS should
use those ideas through DuckDB first: Parquet pruning, parallel scans, planner
choice, memory accounting, disk spill, and typed SQL composition already exist.

The comparison is not “twenty naive LEFT JOINs” versus a bespoke index. The
DuckDB plan can prefilter bins/tiles, scan only selected row groups, merge sorted
numeric keys, and perform interval/range joins at the right cardinality. A custom
cache or mmap format is justified only by a reproducible workload where this path
cannot meet latency or throughput.

The canonical provider receipt, exact/hashed lane split, assembly gate, and bounded
chromosome/tile scheduling contract are defined in `design/duckvep.md`. Peak RSS and
dbSNP/TOPMed-scale projections are numeric claims and therefore live only in the rendered
`benchmarks/benchmark_variantkey_join_overlap.md` report.

Interactive serving is also compatible with the architecture: keep named
transcript models resident; use prepared DuckDB queries, zone-map-pruned tiles,
and a small bounded hot-tile cache; expose the database through a DuckDB server
such as the existing DuckNNG/duckdb-quack work. Random lookup latency still needs
a direct benchmark against specialized stores.

The integration must use the stable DuckDB extension API. If whole-stream query
state cannot be expressed there, document the exact missing lifecycle operation,
the unstable API candidate, ABI coupling, failure modes, and stable fallback
before using it.

Tracked at https://github.com/RGenomicsETL/duckhts/issues/94.

## 13. Pipeline parallelism

DuckDB data chunks are one source of parallelism, not the only pipeline
opportunity. A future end-to-end plan can overlap:

1. htslib decompression and VCF/BCF decode;
2. ALT preparation and numeric-key construction;
3. supplementary tile prefetch/scan;
4. independent C consequence work on disjoint contig/range partitions;
5. compact result materialization;
6. late identifier/string/HGVS rendering and output compression.

Tiles can be double-buffered so the next input/source tile is prepared while the
current C tile executes. Parallel partitions must preserve a deterministic final
order and may not share mutable model/workspace/source iterators. Haplotype
partitions need transcript-end and phase-set-aware splits; arbitrary range splits
can cut an active haplotype and are not equivalent.

The implemented independent-event adapter has been exercised as one ordered branch and
as four disjoint ordered branches on pinned physical P cores. Concurrent callbacks borrow
exclusive workspaces while sharing one immutable model. Their commutative full-row
fingerprints prove the same `(input_variant_index, annotation_index, public row)` multiset
at transcript distances 0, 5,000, 10,000, and 50,000 bases; they do not prove global
emission order. The widest run expands 517,097 inputs into 88,784,213 rows. Source
inspection establishes immutable model sharing, while the process-level RSS measurement
shows that no model-sized increase was observed for this model, corpus, and four-branch
configuration. It is not allocation attribution or permission to share mutable faidx
handles, iterators, reference windows, or haplotype state.

At source `e25c151`, the default 5,000-base run improves from 2.106 seconds on one
P core (245,535 input alleles/s; 12.59 million output rows/s) to 0.632 seconds on
four P cores and four ordered partitions (818,191 input alleles/s; 41.96 million
output rows/s). Peak RSS changes from 5,446,084 to 5,453,280 KiB. At 50,000 bases,
the same corpus improves from 4.161 to 1.342 seconds while emitting 88,784,213 rows;
the four-worker RSS premium is 8.16 MiB. The rendered throughput report is the numeric
authority for all four distances and their row-multiset fingerprints. Peak RSS is GNU
`/usr/bin/time -v` maximum resident set size for the complete one-process R/DuckDB run;
it is not a per-allocation C-heap measurement.

The current scalar adapter restarts at DuckDB vector edges. A stateful table-function
or another stable-lifecycle surface could preserve the sweep across chunks, but it
must be explicit and ABI-stable. DuckDB Appender versus direct vector production
is a later materialization benchmark, not a consequence-kernel decision.

## 14. Current code map

- `src/duckvep/duckvep_ensembl.c`: registers the SQL relation compiler and receipt.
- `src/duckvep/duckvep_model.c`: loads/validates owned arrays, builds cgranges,
  manages named models and worker workspaces.
- `src/duckvep/duckvep_annotate.c`: stable DuckDB scalar adapter and rich or
  compact list materialization over the same kernel cursor.
- `src/duckvep/duckvep_variant_tile.c`: pure-C sorted tile and allele/ID pools.
- `src/duckvep/kernel/include/duckvep_kernel.h`: host-neutral public C ABI.
- `src/duckvep/kernel/src/duckvep_event.h`: lossless raw/feature/edit geometry.
- `duckvep_sweep.c`: sorted candidate sweep.
- `duckvep_classify.c` and `duckvep_projection.c`: topology, splice, and coordinate
  projection.
- `duckvep_delta.c`, `duckvep_coding.c`, `duckvep_codon.c`: edit/CDS/peptide facts.
- `duckvep_effect.c` plus generated rule/metadata files: fact-to-SO program,
  mature-miRNA predicate, and transcript/variant-induced NMD results.
- `duckvep_haplotype.c`: multi-edit partition, mutation, and translation mechanics.
- `duckvep_sv.c`: typed structural and paired-breakend facts.
- `test/duckvep/property/`: fixed/randomized pure-C mechanics.
- `test/duckvep/conformance/`: pinned VEP extraction, generated witnesses,
  statistical differentials, and ledgers.

## 15. What the merged consequence series established

The implementation arrived in several reviewable verticals, but the resulting
architecture has one authority rather than one authority per pull request:

- https://github.com/RGenomicsETL/duckhts/pull/116 made the generalized
  edit/CDS/peptide context authoritative for production MNVs and indels, added
  complete transcript flanks, explicit reason codes, exact VEP feature geometry,
  and statistical state exploration.
- https://github.com/RGenomicsETL/duckhts/pull/117 added the compact numeric SQL
  projection over the same kernel cursor; it did not add another consequence path.
- https://github.com/RGenomicsETL/duckhts/pull/118 applied the exact VEP core
  transcript filter, added mature-miRNA source attributes and resident segments,
  fixed insertion/incomplete-CDS/UTR states, and recorded the first exact indexed-
  cache ClinVar slice.
- https://github.com/RGenomicsETL/duckhts/pull/120 sourced all VEP/BioPerl codon
  tables and applied the supported single-amino-acid Ensembl Translation SeqEdits.
- https://github.com/RGenomicsETL/duckhts/pull/121 closed the declared GRCh37 and
  *P. falciparum* differentials; https://github.com/RGenomicsETL/duckhts/pull/123
  tied those exact results to the generated receipt and conformance ledgers.
- https://github.com/RGenomicsETL/duckhts/pull/124 cached validated CDS projections,
  https://github.com/RGenomicsETL/duckhts/pull/125 reused resolved intron geometry,
  and https://github.com/RGenomicsETL/duckhts/pull/126 reused an already established
  coding fact for NMD. Each change retained before/after rows in the rendered
  throughput data; controls and regressions stayed visible.
- https://github.com/RGenomicsETL/duckhts/pull/128 reproduced the exact feature-
  coordinate state read by VEP's NMD plugin, including reversed insertion ranges,
  defined zero coordinates, and the narrow safe cache-reuse rule.
- https://github.com/RGenomicsETL/duckhts/pull/129 added typed single-locus
  structural annotation and statistical DEL/DUP/tandem-DUP/INV/CNV/INS campaigns;
  https://github.com/RGenomicsETL/duckhts/pull/130 closed every retained GRCh38
  ClinVar small-variant discrepancy.
- https://github.com/RGenomicsETL/duckhts/pull/131 added the paired-breakend SQL/C
  path, bound the full 41-term VEP registry without inventing a generic fallback,
  and matched all 91,428 transcript pairs in the generated BND campaign.
- https://github.com/RGenomicsETL/duckhts/pull/132 merged release-matched funcgen
  preparation and receipts, a separate resident feature sweep, exact small/SV VEP
  differentials for all six regulatory/TFBS terms, the compact 56-byte tagged result
  row, rendered performance/conformance evidence, and the review hardening that rejects
  a NULL transcript sweep model. The 2026-07-19 consequence-closure salvo then sourced
  VEP's paired-BND feature behavior: RegulatoryFeature and MotifFeature candidates use
  both the shifted local point and verbatim mate point, share the explicit event/object
  pair evaluator, union allele-level rows per object, and reproduce its exact-local,
  exact-mate, both-exact, and mate-exact/local-close states. Structural tandem repeats
  also retain their own typed event while reproducing VEP's gain/insertion predicate
  algebra. The same salvo used HPRC v2 long literal alleles to show that ordinary
  `VariationFeature` deletions can enter the shared complete-feature-ablation fact,
  and that equal-length containing spans retain VEP's endpoint-UTR comparisons even
  when the corresponding UTR interval contains no base. A randomized pure-C property
  varies both strands, transcript lengths, uploaded padding, exact transcript endpoints,
  and representations longer than 5,000 bases.

The mature-miRNA and peptide-edit paths show the intended way to add Ensembl facts:
keep the source attribute visible in the prepared relation and deterministic receipt,
project it once, pack the small execution fact into immutable arrays, and feed a bit
to the existing rule program. There is no per-variant attribute parsing and no
term-specific consequence engine.

## 16. Remaining implementation frontier, in dependency order

### A. Consequence conformance as a permanent regression gate

Issue: https://github.com/RGenomicsETL/duckhts/issues/93

The declared dbSNP, GIAB, both GRCh38 ClinVar, GRCh37, and *P. falciparum*
samples are exact, as are the declared structural, paired-breakend, regulation/motif,
NMD, and long-literal distributions. The registered VEP-116 consequence program is
closed for those contracts. Fresh seeded distributions continue as regression and
scope-extension evidence rather than an open-ended release blocker. The
anti-overfitting claim depends on random placement around exon-intron junctions,
CDS ends, transcript ends, rare transcript attributes, allele lengths, strands,
and call-order states—not only a growing regression file. Any new discrepancy must
be minimized and classified as a wrong source fact, event geometry, projection,
sequence fact, SO rule, VEP source-selection difference, or intentionally
unsupported input before code changes. Do not add a term-shaped branch whose only
evidence is the witness that requested it.

### B. Routine consequence throughput engineering

Issue: https://github.com/RGenomicsETL/duckhts/issues/95

The final GRCh38 receipt, model load/RSS accounting, topology, coding-heavy, and
annotation-dense workloads now exist. Performance remains a normal engineering backlog;
it no longer justifies changing consequence semantics or delaying the next semantic
vertical. The original one-million-input-allele single-core target is not met on every
full-model workload, so reports continue to publish input rate beside candidate and
output-row denominators rather than declaring the engine finished in a performance sense.

Retain stage counters and hardware counters for future work: candidates,
topology-only pairs, CDS projections, splice projections, edited and translated
bases, result bytes, branch misses, cache misses, and DuckDB materialization time.
Profile fresh and warm sequence access separately. The current ranked avenues are:

1. write compact SoA columns directly into DuckDB vectors instead of building a rich
   internal row and then materializing a nested stable-API list;
2. keep the common coding-SNV path lean while centralizing its effect composition so the
   optimized and generalized routes cannot drift;
3. reuse prepared point-classifier facts where profiling shows repeated exon/splice work;
4. separate HGVS fact production from string rendering and reuse allele-level genomic
   shift work across transcript rows where VEP semantics permit it;
5. preserve sweep state across DuckDB vectors or expose a stateful table-function surface;
6. make disjoint ordered partition construction a first-class execution helper and keep
   strings, stable identifiers, and supplementary joins late.

These are profiling-informed avenues, not promises that SIMD is the universal answer.
Exploratory local sampling suggested low branch-miss and cache-miss rates in the candidate
sweep and a small candidate-discovery share, but it is not retained as a checked comparable
benchmark and therefore supports no quantitative performance claim. Exon projection,
variable edits, and nested result writing remain gather/branch/materialization work. Every
optimization retains full-row fingerprints, VEP differentials, held-out properties, and
RSS accounting.

### C. More species and genetic-code combinations

Issue: https://github.com/RGenomicsETL/duckhts/issues/119

The importer and translator are species-general in shape, and the full current
evidence covers GRCh38, GRCh37, and *P. falciparum* with codon tables 1, 2, 4, and
11. Add source-faithful models and indexed-cache samples for mouse, Drosophila,
Arabidopsis, and organisms exercising other BioPerl codon tables. Every new row
must pin the Ensembl/Ensembl-Genomes release, assembly, FASTA, transcript filter,
receipt, observed table IDs, and VEP oracle.

### D. Phased executor and Haplosaurus

Issue: https://github.com/RGenomicsETL/duckhts/issues/92

The linear reverse builder is complete. Add GT/PS grouping, sparse carrier state,
transcript watermarks, prefix sharing, unique-leaf translation, source-call and
phase provenance, explicit conflicting-overlap state, phase policies, public SQL,
and VEP Haplosaurus/bcftools-csq differentials. The input receipt and invalidation
unit must distinguish independent alleles from a changed phase block.

### E. Structural input preparation after consequence closure

Issue: https://github.com/RGenomicsETL/duckhts/issues/98

Exact typed single-locus SVs including structural tandem repeats, paired BND,
all 41 registry bits, and the core regulatory/motif lane for small, structural,
and paired-BND events are implemented. The paired feature lane is pinned by a
shared explicit-pair API plus deterministic and randomized all-feature oracles.
The held-out executable-VEP BND/regulation corpus is the semantic authority for
its chromosome-breakpoint `feature_truncation` behavior and its close-local
`intergenic_variant` fallback, not a second consequence implementation.

What remains here is lossless input preparation and downstream use: expand bounded
repeat records once from repeat unit/count metadata, retain unexpanded repeat
payloads, preserve structural inserted sequence and `CIPOS`/`CIEND`, parse raw BND
ALT/orientation outside the kernel, add real fusion corpora, and compose structural
records with the phased executor. VEP 116's registered consequence predicates use
nominal structural coordinates and do not consume the inserted payload, so these
columns must not be presented as missing SO semantics. Keep generating real-VEP
corpora by event class; do not infer typed facts from ALT strings in the inner loop.

### F. Supplementary streams

Issue: https://github.com/RGenomicsETL/duckhts/issues/94

Build one end-to-end exact source and one interval source with bounded memory and
stable API only. Measure tile pruning, merge/sweep, cache misses, materialization,
and dozens-of-sources behavior before considering a custom physical format.

### G. HGVS as a sibling consumer

The shared lower-level authorities are implemented: one prepared event defines semantic
allele geometry, and one CDS edit-set authority drives sequence replay. Consequence uses
those lower-level facts directly. Independent-event DNA/protein HGVS wraps them in
`duckvep_transcript_edit_t`, an HGVS-facing carrier that additionally records VEP's
endpoint-clipped transcript-slice coordinates. Transcript projection and CDS attachment
are separate stages: the adapter builds the former for HGVSc and attaches the latter lazily
for no-FASTA reference validation or HGVSp. This wider carrier is not constructed in the
consequence hot path, and the phased executor does not exist yet; the future phased
executor must consume the same prepared CDS edit set rather than acquire another
allele-trimming or coordinate authority. Typed `c.`/`n.` facts, VEP's genomic-reference
transcript 3-prime shift, typed `p.` facts, and allocation-free renderers have deterministic
and randomized pure-C coverage.

The kernel-opened model is the canonical prepared model exposed internally to the SQL
adapter. It derives or validates one immutable first-stop coordinate for each coding CDS.
The coding context exposes a single edit as a virtual alternate CDS, scans unchanged spans
and the edited payload codon-by-codon, and borrows the cached reference stop when the edit
cannot affect it. Transcript and protein strings share one worker-owned render buffer and
normally need one render call. These are execution choices beneath the typed edit/fact
contracts; neither HGVS nor the adapter gains a second transcript or sequence authority.

The consequence predicate sidecar is only complete enough to skip a later delta when CDS
length is unchanged or it positively records frameshift. A missing frameshift bit on a
length-changing splice overlap is not negative evidence: the full delta must run. This
condition is centralized in
`duckvep_sequence_delta_consequence_flags_complete_for_hgvs(...)` so future local-CSQ and
phased consumers cannot repeat the closed-world mistake.

`duckvep_annotate_hgvs(...)` is the cumulative public adapter. It returns the compact
consequence row plus DNA/protein suffixes and structured status/reason fields from the same
candidate pass. `duckvep_model_load(...)` may bind an existing indexed FASTA using an exact
sequence-region ordinal/name/length projection; the loader validates it, and worker-local
faidx handles supply reusable genomic reference windows. The published model pins open read
descriptors for the validated FASTA, `.fai`, and optional `.gzi`. Linux workers reopen those
descriptors through `/proc/self/fd`; Windows workers use the resolved source while a
deny-write handle remains open. Other POSIX workers open the resolved source independently
and compare its identity around lazy open rather than sharing `/dev/fd` seek state. External
replacement or in-place mutation is outside the immutable-model contract on those systems
and fails when the identity change is observed. Contained sorted requests reuse a
bounded forward-read-ahead window. Versioned transcript/protein identifiers remain cold
relational columns joined by transcript ordinal.

The adapter derives two borrowed views from one bounded faidx fetch. The first is exactly
VEP 116's constrained +/-1000 `_genomic_shift` slice. The second additionally covers the
complete uploaded REF, shifted endpoint-clipped feature sequence, and both adjacent copied
sources for duplication recognition. Passing the wider lookup view into `perform_shift`
would let retained VCF padding change HGVSc placement; limiting duplication recognition to
the shift view would misrender copied insertions longer than 1000 bases as ordinary `ins`.

This is still not a complete VEP-HGVS promise. Fixed executable witnesses cover
position-one right anchors and endpoint-clipped transcript slices, while a strict GRCh38
chromosome-21 ClinVar run found zero HGVSc/HGVSp differences across 56,998 transcript
pairs. That closes the exercised independent-event distribution, not genomic `g./m.`,
RefSeq RNA-edit-aware transcript shifting, structural/BND HGVS, cross-species/assembly
coverage, or compound/phased HGVS. The renderer must not recover biology by parsing SO
strings or alter the uploaded event used for consequence prediction.

The pinned ferro-hgvs v0.9.0 source at
https://github.com/fulcrumgenomics/ferro-hgvs/tree/278e2c11134e3b49067d0c334f650c7c29db9cbe
is a secondary specification oracle and corpus-design reference. It is deliberately
not linked into the C/CRAN/Wasm implementation and cannot replace executable VEP 116
where VEP preserves historically odd output. The differential matrix is therefore
three-way: VEP output compatibility, ferro-hgvs normalization/parsing where its
contract applies, and independent apply-then-diff sequence replay.

## 17. The next recommended vertical

With the declared VEP-116 consequence predicates and independent-event
`c.`/`n.`/`p.` HGVS distribution closed, the next semantic vertical is the phased
executor:

1. add the bcftools `--local-csq` field projection over the existing independent-event
   transcript-edit facts when VariantStory/Talos compatibility needs it; this is an edge
   projection, not another consequence engine;
2. drive the shared edit IR from GT/PS-grouped transcript state and compare complete
   attributed results with VEP Haplosaurus and haplotype-aware bcftools csq; and
3. render haplotype HGVS from the completed alternate transcript/CDS/protein state.

The phased executor remains tracked in
https://github.com/RGenomicsETL/duckhts/issues/92. Independent small variants,
exact single-locus SVs, and paired BND already share one event/projection/fact
authority; the existing linear multi-edit CDS builder is ready to be driven by a
stream that owns genotype grouping and transcript lifetime.

Start with an explicit input relation containing model, event identity, sample,
GT-derived haplotype lane, phase set/policy, and sorted genomic coordinates. Project
each allele once, keep reference haplotypes implicit, retain sparse non-reference
carrier state until the stream passes the transcript end, share identical edit
prefixes, and translate each unique completed leaf once. Output must retain every
contributing variant and expose VEP-compatible versus strict phase handling as
named policies rather than silently treating `/` as phased.

Validation must compare complete CDS/protein haplotypes and attributed consequence
sets to the pinned Haplosaurus executable across same-codon substitutions, open and
restored frameshifts, overlapping edits, multiple samples, ploidies, phase sets,
transcript strands, and chunk splits. Fresh randomized edit sets remain the anti-
overfitting lane. Randomized families must publish observed state counters as well as
trial totals so a passing seed cannot hide absent edit shapes or interactions.
Performance is cumulative and separately records compact consequences, HGVS numeric
facts, rendered HGVS bytes, bcftools-local projection, and phased execution. The
phased lane reports input records and ALT alleles, projected edits, active carrier
states, unique haplotype leaves, translated bases, output haplotypes, thread/core
pinning, and memory high-water marks for both pre-sorted input and DuckDB sort plus
execution.

The existing GIAB, coding SNV, non-SNV, mixed, structural, and BND workloads remain
controls while that executor is added. The current mixed lane's 135,044 alleles/s
and 3.89 million rows/s describe different constraints, and neither may be hidden.
The rendered benchmark report remains the source for any performance claim.
