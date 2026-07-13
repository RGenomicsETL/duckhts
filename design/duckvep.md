# DuckVEP architecture

Status: current implementation contract. Code and tests are authoritative for implemented
behavior; the remaining gaps are linked by full URL in the relevant sections below.

## Purpose and compatibility

DuckVEP is the deterministic variant-consequence kernel inside DuckHTS. It targets Ensembl
VEP 116 (`57ea5c52340acc1f156267f810ad162e26597082`) with Ensembl core
`c0cf13daa961d80584bad797b2eb0ff3a7500ef3` and variation
`2fb834b987ede3824e200197a838ce11e91aeb4b`, initially for human GRCh37/38 and codon
tables 1 and 2.

The kernel emits structured transcript consequences. VEP-compatible CSQ text and HGVS are
edge projections, not the internal representation. VEP 116 is the behavioral authority;
pure-C properties and bcftools csq supply independent checks for mechanics and phased edit
state.

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
- `attrib_type`, `transcript_attrib`, and `translation_attrib` supply consequence-relevant
  flags and exceptional sequence edits.

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

`duckvep_ensembl_transcripts(...)` currently selects every `is_current = 1` transcript on
the supplied regions. For each transcript it:

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
7. retains up to three transcript-oriented bases after the CDS for VEP's terminal-stop
   rule.

The compact flag word records translation presence, protein-coding/NMD/miRNA biotype,
incomplete CDS ends, selenocysteine, stop readthrough, RNA or peptide edits, MANE,
GENCODE, CCDS, readthrough-transcript, and upstream-start state. An RNA edit on either the
transcript or translation attribute path is the same exceptional state.

For GRCh38 MANE attributes, the prepared DuckDB relation also retains the attribute value
as `mane_select_refseq` or `mane_plus_clinical_refseq`. That value is the paired versioned
RefSeq transcript accession. The builder rejects multiple or empty mappings for a MANE
flag. Only the selection bits enter the resident C model; Ensembl/GENCODE and RefSeq
identifiers stay in the cold relation for late SQL projection and future RefSeq HGVS.

Reference-derived CDS bytes are withheld when Ensembl records an RNA edit,
selenocysteine, stop readthrough, amino-acid substitution, or initial-methionine edit that
the C kernel does not yet implement. The transcript, coordinates, flags, and an explicit
reason remain in the model. Unsupported reference alphabet has the same fail-closed shape.
This is different from malformed topology: bad coordinates, strand, exon phase or rank,
translation bounds, or incomplete sequence reconstruction abort the build.

The three post-CDS bases cost three bytes per sequence-backed transcript. They are enough
to reconstruct the original translation endpoint after a short terminal deletion; they
are not a general UTR cache. An edit needing more sequence remains unresolved.

### Prepared relations and publication

The builder returns one row per transcript. The row contains the hot transcript fields,
source and stable identifiers, biotypes, the optional prepared CDS, and an ordered nested
exon list. Callers persist two canonical tables and derive the three sorted loader
projections from them:

- region ordinal, and sequence length when complete-coverage claims are requested;
- transcript ordinal, region, span, strand, gene ordinal, flags, optional CDS fields,
  codon table, and post-CDS bases; and
- transcript ordinal plus each exon's genomic span, cDNA span, phase, and end phase.

The `core_schema` argument names relations, not a transport. It can point at tables loaded
from Ensembl's tab-separated MySQL dumps, or at a read-only MySQL catalog attached through
DuckDB's `mysql` extension. Downloading, attaching, and staging stay outside the model
builder so extension builds remain offline and the same validation runs for either source.
The builder does not require a MySQL server once those relations and the matching reference
chunks have been persisted.

`duckvep_model_receipt(...)` checks dense ordinals and region/transcript agreement. It
records the declared source, release, assembly, transcript filter, source-manifest hash,
reference hash, model counts, and a deterministic hash over every prepared model field.
There is no timestamp: identical declared inputs must produce the same receipt.

The checked-in acceptance fixtures under `test/data/duckvep/ensembl_core/` are about 116
KiB. The GRCh38 fixture contains complete release-116 MT and `HG2047_PATCH` source rows;
it covers mitochondrial fail-closed edits, an ordinary multi-exon CDS, and the real
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

### Current boundary

The implemented builder does not yet:

- stage Ensembl `table.sql` and dump files into DuckDB;
- reproduce VEP 116's exact transcript selection, including stable-ID, `estgene`,
  `artifact`, `readthrough_tra`, and `otherfeatures`/RefSeq rules;
- apply Ensembl RNA and peptide sequence edits;
- derive codon table from `seq_region_attrib` (mitochondrial names currently select table
  2, all other regions table 1);
- preserve the complete Ensembl xref, protein-feature, supporting-feature, attribute, and
  alternate-transcript relations; or
- import variation, phenotype, regulatory, or other supplementary annotations.

Those richer facts should remain typed DuckDB relations joined by numeric source IDs. They
do not belong in every resident C transcript record. Exact VEP-compatible selection and
the richer Ensembl relation set can therefore be separate named products built from the
same staged release. Production-corpus conformance and throughput remain tracked at
https://github.com/RGenomicsETL/duckhts/issues/95.

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

`duckvep_annotate(model, seq_region, position, reference, alternate [, distance])` is the
current stable-C-API adapter for independent biallelic small variants. It copies one DuckDB
vector into compact arrays, splits on model/contig/distance/order changes, seeds the first
candidate set through cgranges, and advances a sorted C sweep. The SNV point path keeps a
per-transcript exon rank and advances it monotonically; other spans use the exhaustive
classifier. Both use the same splice predicates and are property-checked for equality.

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

The current coding path has debt: MNVs and indels first use the generalized edit/CDS/peptide
context and then call older shape-specific code when that context is unresolved. Header
comments and dead route values also describe a superseded equality gate. Issue #93 owns the
removal of that dual authority: extend the edit interpreter or emit an explicit unresolved
result; do not silently choose a different classifier.

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
a boolean. This is the pinned VEP positional policy, not a direct molecular assay or a
claim that every transcript follows one universal NMD rule.

## Phased edits

The pure C mutation core already applies several non-overlapping CDS edits in reverse CDS
order, translates once, and partitions interactions while the frame is displaced or the
next edit touches the same alternate codon. The missing stream groups edits by
`(model, transcript, sample, phase_set, haplotype)` and retains all contributing variant
IDs.

Sorted input bounds lifetime: a transcript's phased state can be finalized once the stream
passes its end. Reference paths stay implicit; non-reference paths should share compact
edit prefixes across samples and translate each distinct leaf once. GT/PS decoding,
arbitrary ploidy, transcript-close flushing, and VEP Haplosaurus comparison are tracked at
https://github.com/RGenomicsETL/duckhts/issues/92.

## Structural events

Small edits and structural events share overlap, projection, provenance, and output, but
not one lossy event shape. Deletion, duplication, inversion, CNV, breakend pairs, inserted
sequence, and repeat changes retain their typed geometry. A breakend is two loci, not a
wide interval or symbolic point.

The kernel currently classifies single-interval structural geometry; the SQL adapter does
not yet accept the typed event fields. Full VEP-116 structural behavior and the remaining
SO vocabulary are tracked at https://github.com/RGenomicsETL/duckhts/issues/98.

## Supplementary annotations

Exact sources such as dbSNP, ClinVar, and gnomAD should remain sorted numeric-key streams;
range sources should remain interval streams. They share prepared variant keys and final
output, not a universal physical format. DuckDB/Parquet pruning, parallel reads, memory
management, and spill are the first implementation. A custom cache requires a measured
workload that those facilities do not meet.

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
- `make bench-duckvep-release-parquet` reads the official Ensembl variation consequence
  VCF through typed CSQ columns and records complete versus consequence-only Parquet size,
  checksum, cardinality, and elapsed time without committing the large artifacts.
- Throughput reports must state input variants, variant-transcript pairs, output rows and
  bytes, threads, model, and machine. Site-level rate and cohort haplotype-update rate are
  different measurements.

Open conformance and throughput work is tracked at
https://github.com/RGenomicsETL/duckhts/issues/93 and
https://github.com/RGenomicsETL/duckhts/issues/95.
