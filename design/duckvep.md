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

- importing Ensembl relations and FASTA, filtering transcripts, verifying receipts, and
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

## Transcript model

The production model originates from the pinned Ensembl core schema and matching genomic
FASTA. `duckvep_ensembl_regions(...)` and `duckvep_ensembl_transcripts(...)` read the raw
`coord_system`, `seq_region`, `gene`, `transcript`, `exon`, `exon_transcript`, `translation`,
`attrib_type`, `transcript_attrib`, and `translation_attrib` relations directly. GFF plus
FASTA may construct reduced test models, but missing Ensembl facts remain explicit.

DuckDB prepares three dense execution projections:

- sequence-region ordinals;
- transcript spans, strand, gene ordinal, flags, CDS bounds, prepared CDS bytes, codon
  table, and up to three transcript-oriented bases immediately after the CDS; and
- exon membership with genomic and transcript-oriented cDNA spans plus phase.

The Ensembl builder matches every supplied FASTA contig to exactly one same-length region
on the requested assembly, imports every current transcript on those regions, splices both
strands from tiled reference chunks, and follows Ensembl `translateable_seq` phase rules.
Ensembl uses exon phase `-1` when translation begins after 5-prime UTR in the same exon;
that means zero prepended bases, not an invalid model. RNA/peptide edits not yet implemented
by the C kernel keep their flags and CDS span but deliberately withhold CDS bytes with an
auditable reason.

`duckvep_model_receipt(...)` hashes both prepared relations, the reference/source hashes,
and the declared transcript filter without a timestamp. The current `duckvep_model_load`
function then reads three sorted projections from committed, non-temporary relations on a
private connection, validates and narrows them, builds one transcript interval index, and
publishes an immutable named model. Several models may coexist in one database instance;
stable IDs and provenance remain ordinary DuckDB columns.

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

The source database and matching FASTA are staged once; annotation processes only open the
prepared relations. Exact production-corpus conformance and throughput remain tracked at
https://github.com/RGenomicsETL/duckhts/issues/95.

The three post-CDS bases cost three bytes per transcript and let VEP's terminal-stop rule
reconstruct the codon at the original translation endpoint after a short deletion. They
are not a general UTR cache: an edit needing more sequence remains unresolved.

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
- Throughput reports must state input variants, variant-transcript pairs, output rows and
  bytes, threads, model, and machine. Site-level rate and cohort haplotype-update rate are
  different measurements.

Open conformance and throughput work is tracked at
https://github.com/RGenomicsETL/duckhts/issues/93 and
https://github.com/RGenomicsETL/duckhts/issues/95.
