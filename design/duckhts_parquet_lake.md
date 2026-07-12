# DuckHTS Parquet and DuckLake

Status: current contract for DuckHTS Parquet write-format version 1 and open design for native
HTS regeneration and DuckLake-backed cohort workflows. `functions.yaml` is authoritative for
the implemented converter helpers.

DuckLake and Quack evolve independently of DuckHTS. Revalidate their SQL and transactional
behavior against pinned extension versions before making a support claim; this note defines
DuckHTS's storage boundaries, not their APIs.

## Responsibility boundaries

- DuckHTS parses HTS records, preserves or corrects format headers, records source provenance,
  and will own any future semantic HTS writer.
- Parquet provides portable typed columns, row-group statistics, and file key-value metadata.
- DuckLake may catalog immutable Parquet files, snapshots, schemas, and partition values.
- Mutable annotations, reviews, and analysis state belong in ordinary tables, not Parquet
  footers.

DuckHTS Parquet is ordinary Parquet, not a private sparse-array format. A custom cache, anchor
cell layout, or index is justified only by a measured workload that DuckDB projection,
row-group pruning, joins, and explicit interval tables do not meet.

Semantic HTS reconstruction is the target. Byte-identical compressed round trips are not
promised unless a future raw-record mode defines and tests that stronger contract.

## Write-format version 1

Every format-aware conversion can record:

- `duckhts_write_format_version`;
- `duckhts_source_format`, `duckhts_reader`, and `duckhts_source_path`;
- `duckhts_header_corrected`;
- `duckhts_filter`, `duckhts_columns`, and `duckhts_partition_by`; and
- one family header key: `vcf_header`, `sam_header`, `gff_header`, or `tabix_header`.

Raw VCF/BCF header capture includes the final `#CHROM` sample line. Without it, typed columns
and metadata are insufficient to reconstruct sample-bearing VCF/BCF semantics.

Callers may supply corrected header text and additional key-value metadata. The map/list path
is the offline-safe interface; JSON-file convenience may depend on DuckDB's JSON extension.
Caller metadata must not accidentally override reserved `duckhts_*` or family header keys. An
intentional override is a correction and needs provenance.

Published data files are immutable. Updating footer metadata means writing a replacement file
with merged metadata and, when cataloged, registering it in a new snapshot. Footer mutation is
not an append-only analysis log.

An incompatible layout or interpretation change requires a write-format version bump. Added
optional columns do not require a bump when old readers can still interpret the stable core.
The exported column list must be sufficient for a future writer to reject an incomplete
semantic reconstruction rather than guess.

## Layout guidance

Partitioning and row groups solve different problems. Partitions avoid whole files; row-group
statistics avoid ranges within a selected file. Start with sorted or clustered rows and tune
row-group size against actual predicates before adding directories.

For genomic data:

- keep genomic identifiers and positions as real columns;
- prefer coarse chromosome, study, run, or sample-batch partitions when measurements justify
  them;
- avoid default position-bucket partitions and one tiny directory per sample;
- choose per-sample layout only when sample-selective reads dominate and partitions stay large;
- keep primary variant rows simple, with long ranges in explicit interval tables; and
- do not expand gVCF reference/no-call blocks to one row per base.

For mixed cohort interval queries, coarse chromosome plus sample-batch files sorted by position
is the default experiment. Repeated point lookup may justify derived `variant_presence` or
`variant_interval` facts; it does not justify over-partitioning the primary files by variant.

## Cohort ingestion and N+1

New samples or batches should produce new immutable files. Do not rewrite every older sample
merely because a new sample or allele appears. Register completed batches as one logical cohort
table when the deployed catalog supports it, and compact or recluster through a new snapshot
when small files harm planning or locality.

The corresponding read rule is set-oriented SQL. Do not issue one scan per sample, source
file, or review variant when a join against one logical table can express the workload.
Catalog file statistics, partition values, and Parquet statistics are pruning aids; they are
not substitutes for normalized allele keys or interval facts.

A cohort manifest should record source URI, sample identity, reference build, header and source
hashes, index identity, conversion options, write-format version, and catalog snapshot. The
manifest is an ordinary table so it can evolve without rewriting data files.

## Derived facts and reanalysis

Keep observed facts separate from derived evidence:

- variant observations and gVCF blocks preserve source data;
- normalized allele and interval facts support repeated cohort lookup;
- annotation runs record tool, reference, cache/knowledgebase, parameters, input hashes, and
  the snapshot read; and
- clinical interpretations are append-only results keyed by stable variant and case/cohort
  identity plus an analysis identifier.

Reanalysis writes a new evidence result against a recorded snapshot. It does not mutate the
original observed-variant lake. This preserves the ability to explain what was known at an
earlier analysis point.

External annotations belong in typed tables joined by normalized allele identity or genomic
interval. Parquet key-value metadata is suitable for compact file-intrinsic facts, not VUS
state, dashboards, row-level statistics, or a query index.

For gVCF cohorts, retain the fields that distinguish variant observations from reference and
no-call blocks. Derived `gvcf_blocks` and `coverage_evidence` tables may answer absence and
confident-reference questions without repeated per-sample scans or dense expansion.

## Open decisions

- Define the minimum lossless column contract for native VCF, BCF, BAM, GFF, and tabix writers.
- Define grouping and ordering requirements for reconstructing tidy multi-sample VCF data.
- Decide whether a raw-record mode is worth its storage and schema cost.
- Benchmark staged Parquet plus catalog registration against direct DuckLake writes.
- Specify compaction, snapshot retention, and schema-migration procedures for cohort tables.
- Validate a server-gateway deployment separately from any remote metadata-catalog design;
  Quack is optional and is not part of the local or CRAN format contract.
- Define append-friendly cohort statistics and annotation tables before adding fused ingest
  shortcuts.
