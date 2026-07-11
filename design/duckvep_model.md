# DuckVEP model and execution contract

Status: current implementation guidance. This describes the model DuckDB prepares and the
borrowed arrays consumed by the private C kernel under `src/duckvep/`.

## Compatibility target

- Ensembl VEP 116, commit `57ea5c52340acc1f156267f810ad162e26597082`.
- Ensembl core commit `c0cf13daa961d80584bad797b2eb0ff3a7500ef3`.
- Ensembl variation commit `2fb834b987ede3824e200197a838ce11e91aeb4b`.
- Human GRCh37 and GRCh38; codon tables 1 and 2.

An assembly, an Ensembl annotation release, and a VEP behavior target are separate
identities. Imported data records the source URI, source commit or release, file checksum,
assembly, and FASTA checksum.

## Where transcript models come from

The production importer reads the pinned Ensembl core schema and matching genomic FASTA.
It selects ordinary relations; DuckHTS does not invent a cache file format or mirror the
whole Ensembl database.

The indispensable Ensembl relations are:

- `coord_system`, `seq_region`, and `dna` for assembly identity and reference sequence;
- `gene`, `transcript`, `exon`, and `exon_transcript` for stable identities, genomic spans,
  strand, exon membership, rank, and transcript order;
- `translation` for the selected coding start and end within its boundary exons;
- `attrib_type` and the gene/transcript/translation attribute relations for incomplete CDS,
  RNA edits, selenocysteine, amino-acid substitutions, stop-codon readthrough, MANE, and
  related facts used by VEP predicates; and
- stable-ID/xref relations only for requested output identifiers.

Those relations add facts a generic GFF normally loses: exact translation boundaries,
source exon ranks, multiple translations, incomplete-CDS state, typed sequence edits, and
release-specific selection metadata. GFF plus FASTA remains useful for small tests and
portable reduced models, but it cannot silently claim the same behavior when those facts
are absent.

## DuckDB preparation

DuckDB performs import, filtering, joins, receipt checks, and sequence construction.
Prepared transcript rows are sorted by sequence region and genomic start. Exon membership
is sorted in genomic order with precomputed transcript-oriented cDNA spans.

For each selected translation, preparation builds the transcript-oriented coding sequence
once:

1. join ranked exons to the matching FASTA sequence;
2. reverse-complement minus-strand transcripts;
3. apply declared transcript-level edits in their recorded order;
4. cut the translation boundaries and preserve start-phase padding; and
5. hash the resulting bytes and the inputs that produced them.

The C kernel receives the prepared coding sequence as a borrowed pointer and length. It
does not query FASTA, allocate a transcript string, or rebuild a CDS for each variant.

## Resident model

`duckvep_model_load` takes a name and three SQL queries. They are execution projections,
not a second annotation database:

- sequence regions: one sorted, unique `seq_region UINTEGER` column;
- transcripts: `transcript_index`, sequence region, genomic start/end, strand, gene
  ordinal, flags, optional CDS start/end, optional prepared CDS bytes, and codon table; and
- exons: transcript ordinal, genomic start/end, transcript-oriented cDNA start/end, phase,
  and end phase.

The loader reads the three committed relations in one snapshot on a private connection, so
temporary and uncommitted relations are intentionally out of scope. It checks names, exact
SQL types, ordering, spans, phases, CDS null coupling, and sequence bounds before publishing
the model. Wide DuckDB coordinates are narrowed once to the kernel's human-genome storage:
16-bit sequence-region and exon-count fields, 32-bit one-based inclusive genomic/cDNA
coordinates, 32-bit transcript/gene ordinals, and a shared CDS byte pool. Stable biological
identifiers remain in ordinary DuckDB lookup relations.

Several named models may coexist in one database instance. Model arrays and their cgranges
indexes are immutable after loading. `duckvep_model_drop` refuses to remove a model while a
worker is using it.

## Variant execution

`duckvep_annotate` is a stable-ABI vector function. DuckDB supplies variant columns directly;
the adapter does not execute a nested variant query or materialize the variant relation. For
each DuckDB vector it:

1. validates and copies REF/ALT bytes into one compact pool;
2. splits rows only when model, sequence region, distance, or sort order changes;
3. asks cgranges once for transcripts already crossing the start of each sorted run; and
4. continues with a forward sweep in which transcripts enter once and exon cursors only
   advance.

There is no transcript binary search per variant and no all-exon scan per candidate pair.
Adapter arrays are recycled across vectors, and each model keeps a pool of already-sized
kernel workspaces for concurrent DuckDB workers. Results use bounded, resumable caller-owned
arrays before becoming DuckDB list rows.

The SQL caller should put `unnest(duckvep_annotate(...))` in the `SELECT` list and expand
the returned struct in an outer query. Writing it as a lateral table-function join creates
a more expensive dependent plan and defeats much of the point of the vector adapter.

The stable community-extension ABI does not provide expression-local scalar state across
successive DuckDB vectors. Consequently, the current adapter seeds once per sorted vector,
not once per whole chromosome. DuckDB's prepared-result streaming creation calls are exposed
only by the unstable extension ABI, so DuckHTS does not use them. Carrying the sweep frontier
across vectors must come from a future stable table-in/table-out API or a kernel interface
that accepts an explicit carry state; it is not a reason to make the extension version-locked.

## Supplementary annotations (planned execution rule)

Large exact-match sources such as dbSNP, ClinVar, and gnomAD will remain sorted streams.
One current key per source will participate in a small merge heap. Keys and projected
payload blocks will be discarded behind the variant cursor; the source will never be
resident in full.

Region sources will use the same event sweep. Intervals will enter one global active set
containing their end, source ordinal, and row ordinal, and leave after their end. Memory
will therefore be one projected input block per selected source plus the intervals
overlapping the current position. cgranges may seed a tile or serve a sparse random query,
but dozens of full resident interval trees are not the whole-genome path.

The transcript model will remain shared by all workers. Variant blocks, live transcript
state, supplementary input blocks, and result blocks will be per-worker. Workers will
claim neighboring tiles so their source cursors and decoded blocks are reused instead of
reopened per tile.

## Consequences and phased edits

Alleles are represented as checked immutable byte views and a semantic edit obtained by
maximal common-prefix and common-suffix trimming. The lossless allele/transcript context
owns the projected cDNA/CDS/protein facts, oriented alleles, reference validation, codon
change, and provenance needed by both consequence and future HGVS consumers.

The public vector accepts biallelic A/C/G/T/N substitutions and indels. The ported VEP-116
predicates cover intergenic, upstream/downstream, intronic, noncoding-exon, UTR, splice,
synonymous, missense, start/stop, frameshift, in-frame, and protein-altering cases. When a
coding span exists but its prepared sequence facts cannot settle a predicate, the row says
`status = 'unresolved'` and gives a reason; it is not counted as supported by the statistical
report. A coding REF mismatch is an error. The kernel also contains the single-interval
SV/CNV classifier, but the current SQL vector has no typed structural-event argument yet.

The pure C haplotype path applies all edits to one prepared CDS and translates once. Its
partitioner expects edits already grouped by model, transcript, sample, phase set, and
haplotype; it keeps a block open while the reading frame is displaced or the next edit still
touches the same alternate codon. This covers combined same-codon substitutions and a later
indel that restores a frame. The SQL surface still annotates independent alleles: genotype
decoding, phase grouping, and tile-boundary finalization must be added together rather than
pretending independent calls are haplotype consequences.

## Validation

- `make test-duckvep-kernel` checks sweep sets against brute force, projection against a base
  walk, codons against the genetic code, edits against rebuilt sequences, and fused output
  against the independently checked parts. Sanitizer targets run the same properties.
- `make test-duckvep-kernel-statistical` raises every randomized property to 100,000 trials;
  the seed and trial count are explicit and reproducible.
- `make test-duckvep-differential` runs formal boundary/codon/allele witnesses through VEP
  116 and DuckVEP on the same GFF, FASTA, and transcript, then retains every exact match,
  disagreement, missing/extra emission, and unresolved row.
- `make duckvep-corpus-differential` applies the same union-denominator comparison to a
  deterministic sample from a real VCF and reports exact binomial 95% upper bounds per
  consequence, allele type, and length bin.
