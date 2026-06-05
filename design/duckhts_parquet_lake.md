# DuckHTS Parquet lake round-trip format and DuckLake registration notes

Status: current design note for the implemented R/DBI Parquet converter helper surface (`duckhts_*_convert_parquet_sql(...)` plus R wrappers) and future DuckLake/native-writer work. Native VCF/BCF/BAM/GFF/tabix writers remain future work.

## Design stance

DuckHTS Parquet is a typed, ordinary-Parquet export format with enough metadata for future semantic HTS regeneration. DuckLake is the catalog/snapshot layer that can register those Parquet files without rewriting them. Keep those responsibilities separate:

- DuckHTS owns HTS parsing, corrected/preserved headers, source-format metadata, and future HTS writers.
- Parquet owns portable columnar storage and file-level key-value metadata.
- DuckLake owns table catalog state, snapshots/time travel, schema evolution, file statistics, and registration of premade Parquet files.

Non-goals for this design:

- Do not clone TileDB-VCF/Phoebe-style anchor-cell sparse-array storage. Use typed columns, intervals, Parquet statistics, DuckLake metadata, and optional future indexes instead.
- Do not require Quack or any server path for local/package reproducibility. Quack remains an important deployment/gateway option, especially when the DuckLake catalog or an attached DuckLake database is served remotely.
- Do not make memory-backed or mmap-backed hFILE experiments part of the lake format.
- Do not promise byte-identical compressed HTS round trips. The target is semantic reconstruction unless an explicit future `raw_record`/raw-block fidelity mode is selected.

## Goals

DuckHTS Parquet exports should be reusable as ordinary Parquet files and carry enough metadata for future format-aware regeneration:

- the original or corrected spec header;
- the DuckHTS reader family and Parquet write-format version;
- the logical layout (`tidy_format` for VCF/BCF, selected columns, output partitions);
- SQL filter/provenance strings when a subset was exported;
- caller-supplied correction notes or application metadata.

The initial helper surface is:

- `rduckhts_bcf_convert_parquet()` / `duckhts_bcf_convert_parquet_sql()`;
- `rduckhts_bam_convert_parquet()` / `duckhts_bam_convert_parquet_sql()`;
- `rduckhts_gff_convert_parquet()` / `duckhts_gff_convert_parquet_sql()`;
- `rduckhts_tabix_convert_parquet()` / `duckhts_tabix_convert_parquet_sql()`.

The R helpers are thin DBI wrappers. They call extension-owned SQL builders (`duckhts_*_convert_parquet_sql(...)`) that generate executable `COPY ... TO ... (FORMAT PARQUET)` statements and use Parquet key-value metadata, including `duckhts_write_format_version = "1"`.

## Metadata keys in DuckHTS Parquet write-format version 1

Common keys:

- `duckhts_write_format_version`
- `duckhts_source_format`
- `duckhts_reader`
- `duckhts_r_package_version`
- `duckhts_source_path`
- `duckhts_header_corrected`
- `duckhts_filter`
- `duckhts_columns`
- `duckhts_partition_by`

Header keys by family:

- VCF/BCF: `vcf_header`
- SAM/BAM/CRAM: `sam_header`
- GFF: `gff_header`
- generic tabix: `tabix_header`

VCF/BCF raw header extraction includes the final `#CHROM` sample line. This is required for semantic VCF/BCF reconstruction from typed columns plus preserved header metadata.

Callers can pass `header_text` to store a corrected header and set `duckhts_header_corrected = true`. They can pass a named `metadata` list/vector in R, or `metadata := map([...], [...])` in SQL, to add key-value entries such as correction notes, upstream filter descriptions, project identifiers, analysis identifiers, or cohort identifiers. This map/list path is the primary CRAN/offline-safe metadata interface. `metadata_json_file` remains an optional caller-managed convenience when DuckDB's `json` extension is available; map/list metadata remains the CRAN/offline-safe form.

Avoid using caller metadata to casually overwrite reserved `duckhts_*` or header keys. Overriding those keys should mean an intentional correction with provenance.

### Updating Parquet key-value metadata later

Current DuckHTS methods write Parquet key-value metadata at conversion time through `metadata := map(...)`, R named lists, and optional `metadata_json_file`. DuckHTS does not currently expose an in-place "add metadata to existing Parquet" method.

Treat Parquet data files as immutable after publication or DuckLake registration. Changing Parquet key-value metadata means changing the file footer; portable workflows should rewrite the file(s) with merged metadata and register the replacement files in a new DuckLake snapshot. Do not rely on footer-appending tricks as a supported append-only mutation path, especially for object stores, content-addressed artifacts, or broadly interoperable files.

Use Parquet KV metadata for file-intrinsic facts needed to interpret/regenerate that file: source format, DuckHTS write-format version, preserved/corrected header, selected columns, filters, and partition layout. Use DuckLake tables or ordinary side tables for evolving analysis state: VUS review status, annotation runs, clinical classifications, QC decisions, and other metadata that changes after export.

## DuckLake `ducklake_add_data_files(...)` findings

Current DuckLake has a catalog-only add-file path for premade Parquet files:

```sql
CALL ducklake_add_data_files(
  'lake_catalog',
  'variants',
  'path/to/files/**/*.parquet',
  hive_partitioning => true,
  allow_missing => true
);
```

Important behavior observed from `duckdb/ducklake` source and tests:

- the function accepts a single path/glob string or a list of strings;
- it registers existing Parquet files without copying or rewriting data;
- it maps columns by name, can reject extra/missing columns, and supports `ignore_extra_columns` / `allow_missing`;
- it reads Parquet file statistics and stores DuckLake file/column metadata;
- `hive_partitioning` can be forced on/off; automatic mode discovers Hive-style path segments when possible;
- partitioned DuckLake tables are configured with `ALTER TABLE ... SET PARTITIONED BY (...)`;
- for identity partitions such as `CHROM` or `SAMPLE_ID`, files under `CHROM=1/SAMPLE_ID=HG00096/...` can be registered and DuckLake records partition values;
- for transformed partitions (`year(dt)`, `month(dt)`, `bucket(n, col)`), file paths must contain the expected transform key names/values; missing or incompatible partition values are rejected;
- order of Hive path segments does not matter when all required partition keys are present.

Practical DuckHTS pattern:

```sql
ATTACH 'ducklake:catalog.ducklake' AS lake (DATA_PATH 'parquet-data/');

CREATE TABLE lake.variants AS
FROM read_parquet('cohort/**/*.parquet', hive_partitioning = true)
WITH NO DATA;

ALTER TABLE lake.variants SET PARTITIONED BY (CHROM, SAMPLE_ID);

CALL ducklake_add_data_files(
  'lake',
  'variants',
  'cohort/**/*.parquet',
  hive_partitioning => true
);
```

For DuckHTS-generated partitioned Parquet, prefer identity partitions that are already columns in the reader output (`CHROM`, `RNAME`, `seqname`, feature type, sample batch, etc.). `SAMPLE_ID` is useful when files are large enough to justify per-sample pruning, but blindly partitioning tiny samples into many tiny files creates metadata overhead and slow planning.

### Partitioning versus row-group statistics

Directory/Hive partitions and Parquet row-group statistics solve different pruning problems:

- partitions prune whole files/directories before Parquet footers and row groups are inspected;
- Parquet row-group statistics prune inside an already selected file;
- DuckLake records file statistics and partition values in the catalog, while DuckDB's Parquet reader can still use row-group statistics during the scan.

For genomic positions, sorted row groups should be the first tool. If files are sorted or clustered by `CHROM, POS` and `ROW_GROUP_SIZE` is tuned well, row-group min/max statistics usually make a `floor(POS / 1000000)` partition column redundant for the primary variant table. A 1 Mb position-bucket partition can create thousands of partitions per genome, and combining it with `SAMPLE_ID` can explode into many tiny files.

Default guidance:

- do not partition by `chrom_pos_bucket` by default;
- preserve or produce sorted files by `CHROM, POS` where possible;
- tune row-group size against the query workload, remembering that `ROW_GROUP_SIZE` is a row count, not a base-pair width;
- use coarse `CHROM` partitioning only when it measurably reduces file opens/listing on large or object-store datasets;
- use `SAMPLE_ID` partitioning only when sample-selective queries dominate and each partition remains large enough; for large cohorts, prefer a `sample_batch`/`sample_shard` column over one directory per sample;
- for mixed cohort interval queries, `CHROM` plus sample-batch layout and sorted row groups is usually a better starting point than `SAMPLE_ID` plus position buckets;
- if repeated sub-megabase point lookups dominate, build derived `variant_presence`, `variant_interval`, or `gvcf_blocks` tables rather than over-partitioning the primary Parquet files.

For tidy multi-sample VCF output, the tradeoff is workload-dependent. Without `SAMPLE_ID` partitioning, row groups ordered by genomic position may give excellent region pruning across the cohort, but a single-sample query may touch many row groups. With `SAMPLE_ID` partitioning, single-sample reads prune well, but cohort interval queries may touch many sample partitions. This is why the design should benchmark layouts and treat per-sample partitioning as a workload choice, not a universal default.

## Questions this design must answer

### 1. N+1 population-genomics problem

This section uses the broader TileDB-VCF framing of the N+1 problem, not only the SQL anti-pattern. See TileDB's description of the population-genomics N+1 problem: <https://documentation.cloud.tiledb.com/academy/structure/life-sciences/population-genomics/foundation/key-concepts/n-plus-1/>.

The problem has both write-side and read-side forms:

- **CombinedVCF write-side N+1**: adding sample `N+1` to a merged cohort file can require rewriting the whole combined file, effectively injecting the new sample's bytes into many genomic locations. A novel variant in the new sample can also require adding cohort-wide missing/reference information for all existing samples, causing storage growth and poor ingest behavior.
- **Single-sample-file read-side N+1**: keeping every sample in its own VCF/gVCF avoids rewriting old samples, but cohort queries then perform N separate range searches and assemble scattered byte ranges.
- **gVCF semantics**: gVCF reference/no-call blocks help answer whether an existing sample has evidence at a newly observed variant, but the format alone does not solve the storage-layout and cohort-query problem.

DuckHTS should not answer this by cloning TileDB-VCF's sparse-array model, but it must answer the same operational pressure. The intended DuckHTS/DuckLake pattern is:

- append new samples or sample batches as new Parquet files/directories; do not rewrite old sample files just because a new sample arrives;
- register completed files into DuckLake so the catalog, snapshots, file statistics, and partition values describe the cohort as one logical table;
- use coarse partitions and sorted/clustered files to preserve genomic locality where it matters;
- materialize derived cohort facts such as `variant_presence`, `variant_interval`, `gvcf_blocks`, and `coverage_evidence` for repeated point/range clinical questions;
- compact/recluster in the background when many small ingests degrade locality, without blocking readers or invalidating previous snapshots.

The lake must also avoid ordinary N+1 query loops that open or scan one file per variant, per sample, or per VUS. The intended query pattern is one logical DuckLake table per cohort/fact type and set-oriented SQL over all relevant variants/samples:

```sql
WITH vus(chrom, pos, ref, alt) AS (
  VALUES ('7', 140453136, 'A', 'T')
)
SELECT v.SAMPLE_ID, v.CHROM, v.POS, v.REF, v.ALT, v.FILTER
FROM lake.variants v
JOIN vus u
  ON v.CHROM = u.chrom
 AND v.POS = u.pos
 AND v.REF = u.ref
 AND list_contains(v.ALT, u.alt);
```

DuckLake helps by storing file statistics, partition values, and snapshots. DuckHTS should help by encouraging layouts that make pruning possible. Practical rules:

- register premade Parquet files into DuckLake instead of leaving users with ad hoc file globs;
- batch conversion/registration; do not issue one DBI call per source file when a set of files can be converted or registered together;
- partition coarsely by common predicates (`CHROM`, `RNAME`, `seqname`, cohort/run/batch, sometimes `SAMPLE_ID`), not by high-cardinality variant identifiers or positions;
- keep row groups large enough for efficient scans and statistics, but not so large that predicate pruning becomes useless;
- for point-lookup-heavy clinical workloads, plan a derived `variant_presence`/`variant_interval` fact table or secondary index rather than repeatedly probing raw per-sample files.

This design does not solve arbitrary point lookup by file naming tricks. The answer to N+1 is append-friendly files, DuckLake catalog/snapshot management, set-oriented tables, pruning, derived cohort facts, and background compaction/reclustering.

### 2. DuckHTS/DuckLake equivalent story

The equivalent story is SQL-first lakehouse infrastructure rather than a TileDB-VCF sparse-array clone. The current PR supplies the typed Parquet conversion and metadata layer; the cohort-scale layers below are the intended next pieces.

**Ingestion.** DuckHTS readers are the ingest kernels: htslib parses VCF/BCF/gVCF/BAM/GFF/tabix input, and `duckhts_*_convert_parquet_sql(...)` writes typed Parquet with preserved headers, filters, selected columns, partition layout, and source metadata. Cohort ingestion should create or update a manifest table containing source URI, sample IDs, reference build, header hash, index/checksum information, conversion command, and DuckHTS write-format version. For VCF/BCF cohorts, prefer the tidy layout with `SAMPLE_ID` over expanding multi-sample files into a CombinedVCF-style wide layout.

**Append and layout.** New samples or batches become new Parquet files/directories registered into DuckLake. Old files are not rewritten just because a new sample appears. To preserve pruning, batches should be organized by coarse genomic partitions and sample batches rather than by one tiny file per sample. Practical layouts are benchmark-driven, but the starting point is `CHROM` plus a batch/study/run key, sorted or clustered by genomic position and optionally by `SAMPLE_ID`. If many small ingests accumulate, write a new compacted/reclustered snapshot instead of mutating published Parquet files.

**Reads.** DuckLake exposes the cohort as logical tables. Reads use SQL projection/filter pushdown, partition pruning, Parquet row-group statistics, and DuckLake file/column statistics. For repeated point/range clinical questions, use derived facts such as `variant_presence`, `variant_interval`, `gvcf_blocks`, and `coverage_evidence`; do not repeatedly scan every raw per-sample source file.

**Annotations.** Embedded VCF annotations are materialized as typed columns when selected from `read_bcf(...)` output and can participate in SQL filters and joins where the resulting Parquet/DuckDB type supports it. External annotations should live in ordinary DuckDB/DuckLake tables, not Parquet key-value metadata. Use two canonical join shapes: normalized allele keys (`CHROM`, `POS`, `REF`, `ALT`, or a stable variant key) and genomic intervals (`CHROM`, `start`, `end`). VEP/SnpEff-style annotation tables should record one row per allele/transcript/effect plus annotator version, reference build, cache/source versions, and options. A single annotation table can serve many cohorts when those provenance fields match; only novel variants need new annotation.

**Interval annotations and long ranges.** DuckHTS should not inject TileDB-style anchors into the primary variant table. For gVCF blocks, CNVs, genes, regulatory regions, and other long intervals, use explicit interval side tables, coarse bins/region keys, and native interval kernels such as the existing interval/cgranges facilities. This keeps the primary Parquet layout simple while still allowing fast overlap joins.

**Variant statistics.** TileDB-VCF's `variant_stats`, `allele_count`, and `sample_stats` arrays map to DuckLake side tables. Planned tables include normalized allele AC/AN/IAF partials, raw `CHROM`/`POS`/`REF`/`ALT` counts optionally stratified by `FILTER` and genotype class, and per-sample QC summaries similar to bcftools stats/Hail sample_qc. These should be append-friendly partial sums keyed by ingest batch and DuckLake snapshot, with final statistics computed by SQL aggregation. A future fused ingest path can compute them while records are already parsed, but the table semantics must stay explicit: pre/post-filter, genotype model, normalization rules, reference build, and snapshot.

**Compute.** The local compute story is DuckDB SQL, macros, native DuckHTS table functions/UDFs, and R wrappers. The remote/server story is Quack or ducknng-style gateways running DuckDB/DuckLake near the data. DuckHTS should not require a cloud-only UDF system, but server-side SQL/UDF execution is the right analogue for distributed annotation, statistics, and reanalysis workloads.

### 3. VUS reanalysis

VUS reanalysis should not mutate the original observed-variant lake. Treat the lake as immutable observed facts plus append-only analysis/evidence outputs:

- `variants` / `gvcf_blocks` / `alignments`: observed data facts from HTS inputs;
- `annotation_runs`: one row per annotation/interpretation run with tool versions, reference build, VEP/cache/ClinVar/gnomAD versions, parameters, and the DuckLake snapshot read;
- `variant_annotations` or `clinical_interpretations`: append-only results keyed by stable normalized variant identity plus sample/cohort keys and `analysis_id`.

Reanalysis then becomes a new write against a known input snapshot, not a destructive update:

```sql
SELECT *
FROM lake.variants AT (VERSION => 42)
WHERE CHROM = '7' AND POS BETWEEN 140000000 AND 141000000;
```

Use DuckLake snapshots and the run metadata to answer: "what did we know then?", "what changed since then?", and "which VUS need review under the new knowledgebase?" DuckHTS Parquet metadata records provenance for exported files, but clinical VUS lifecycle management belongs in analysis tables layered on top of those facts.

### 4. VCF statistics and summaries

Small immutable source summaries can be stored at conversion time, but Parquet key-value metadata should not become a row-level statistics database. Good file-level metadata candidates are source record count, sample count, contig count, reference build, normalization status, source checksum, and the command/run identifier.

Most VCF statistics should be computed from typed columns or materialized into ordinary tables:

- one-off checks such as row counts, PASS counts, filter distributions, ALT length summaries, or per-contig counts are cheap enough to compute from Parquet columns because DuckDB reads only the needed columns;
- DuckLake already records table/file row counts, file sizes, column min/max/null-count style statistics, and partition values in its catalog when files are written or registered;
- repeated dashboards, QC gates, clinical review counts, and cohort-level summaries should be stored in append-only side tables keyed by `analysis_id` / snapshot, not hidden in Parquet footers;
- statistics needed for pruning or joining must be real columns, partition values, DuckLake metadata, or derived fact/index tables. Arbitrary Parquet KV entries are not a query index.

So the default is: compute simple summaries on demand; materialize expensive/repeated summaries as lake tables; reserve Parquet KV metadata for compact file-intrinsic facts needed to interpret or regenerate the file.

### 5. Performant writes

The performant write path is bulk `COPY` to Parquet, then optional catalog-only registration with DuckLake. Avoid row-at-a-time appends and avoid producing thousands of tiny files.

Recommended write pattern:

1. Convert a meaningful batch to a staging directory using `duckhts_*_convert_parquet_sql(...)` / R wrappers.
2. Use sensible `COMPRESSION`, `ROW_GROUP_SIZE`, and coarse `PARTITION_BY` values.
3. Register the completed files with `ducklake_add_data_files(...)`.
4. Compact or rewrite later if a workload creates too many small files.

Direct `INSERT INTO lake.table SELECT ... FROM read_bcf(...)` may be useful, but the current portable DuckHTS surface is the extension-owned `COPY ... TO Parquet` builder. The design should benchmark direct DuckLake writes before declaring them the default.

### 6. Time travel

Parquet key-value metadata alone does not provide time travel. DuckLake snapshots do.

Once DuckHTS Parquet files are registered or written through DuckLake, DuckLake can expose snapshot history and versioned reads:

```sql
SELECT * FROM ducklake_snapshots('lake');
SELECT * FROM lake.variants AT (VERSION => 42);
```

Retention policy matters. If an analysis or clinical report needs reproducibility, record the input DuckLake snapshot ID and do not expire that snapshot or its files until the retention policy allows it.

### 7. Schema evolution

There are two schema-evolution layers:

1. DuckHTS Parquet write-format evolution, tracked by `duckhts_write_format_version` and reserved metadata keys.
2. DuckLake table evolution, tracked in the DuckLake catalog with snapshots and column mappings.

Rules for DuckHTS:

- keep a stable core schema for each reader family;
- add optional columns rather than changing meanings in place;
- bump `duckhts_write_format_version` for incompatible layout changes;
- preserve the exported column list in `duckhts_columns` so regeneration can detect omitted required fields;
- use DuckLake `ALTER TABLE` for lake-table evolution when the data is cataloged;
- use `allow_missing` / `ignore_extra_columns` only as explicit migration bridges, not as a way to hide semantic incompatibility.

Column renames and nested-field evolution are DuckLake catalog concerns; HTS semantic compatibility remains DuckHTS's responsibility.

### 8. Quack as a DuckLake catalog/gateway path

Quack matters for remote DuckLake deployments, but it should be treated as a deployment option with explicit validation, not as a requirement for the DuckHTS Parquet format.

Details checked from current DuckLake and Quack docs/source:

- DuckLake's documented catalog database choices are DuckDB for single-client local use, SQLite for multiple local clients, and PostgreSQL for multi-user/remote clients.
- Quack is currently described as experimental/pre-release. It exposes a DuckDB server over HTTP(S); clients can `ATTACH 'quack:host' AS remote`, scan remote tables with projection/filter pushdown, execute DDL/writes, and forward transactions.
- DuckLake attaches its metadata database internally with `ATTACH OR REPLACE {METADATA_PATH} AS {METADATA_CATALOG}`. That makes a Quack-backed metadata path plausible, but it is not yet a documented DuckLake catalog backend in the stable DuckLake catalog guidance.

There are two distinct architectures:

1. **DuckLake behind Quack**: run DuckLake on the server, expose the attached DuckLake database through Quack, and let clients query/write through the remote catalog. This is the safer gateway model because data-path credentials, DuckLake metadata changes, and Parquet reads/writes happen in one server-side DuckDB process.
2. **DuckLake metadata over Quack**: clients attach DuckLake whose metadata path is a remote Quack catalog. This needs conformance testing before support claims, especially around transactions, conflict retries, metadata DDL/DML, `ducklake_add_data_files(...)`, and whether every client can access the same `DATA_PATH` object store paths.

For DuckHTS docs, say: Quack is optional for local CRAN/offline workflows, but important for remote catalog/server deployments. Do not promise DuckLake-over-Quack catalog compatibility until there is an integration test that starts a Quack server, attaches DuckLake through it or behind it, registers DuckHTS Parquet, reads snapshots, and exercises a write/commit conflict path.

### 9. gVCF

gVCF is first-class input, but it must not be handled by dense per-base expansion. A gVCF contains both variant records and reference/no-call blocks. The lake design should preserve both without exploding reference blocks into one row per base.

For write-format version 1, DuckHTS conversion preserves gVCF records as VCF/BCF rows. Users who need gVCF semantics must keep the columns that define block meaning, especially `END`/interval fields and relevant sample FORMAT fields. Future derived tables should separate:

- `variants`: non-reference variant observations;
- `gvcf_blocks`: reference/no-call block intervals with sample-level depth/quality fields;
- `variant_presence` / `coverage_evidence`: derived facts for fast cohort questions such as "who has this variant?" and "who has confident reference evidence here?".

This directly supports VUS reanalysis: absence evidence should come from indexed gVCF block/coverage facts, not from N per-sample gVCF scans and not from per-base expansion.

## Future work

- Native `duckhts_vcf_write(...)` / `duckhts_bcf_write(...)` that consumes wide and tidy DuckHTS Parquet layouts.
- Tidy VCF writer grouping by variant key plus `SAMPLE_ID`, initially requiring sorted/contiguous groups.
- Optional `raw_record` columns for exact text preservation when users need stronger fidelity than semantic reconstruction.
- Derived `variant_presence`, `variant_interval`, `gvcf_blocks`, and `coverage_evidence` tables for point lookup, VUS review, and absence-evidence workflows.
- Cohort manifest tables recording source URI, sample IDs, reference build, header hash, index/checksum information, conversion command, and DuckHTS write-format version.
- External annotation table builders for normalized allele joins and genomic-interval joins, including VEP/SnpEff-style one-row-per-allele/transcript/effect outputs with provenance.
- Optional source-summary helpers that compute compact immutable VCF/BAM/GFF statistics at conversion time and/or materialize append-friendly partial-statistics side tables.
- Benchmarks comparing staged Parquet + `ducklake_add_data_files(...)` against direct DuckLake `INSERT/COPY` paths.
- Quack/DuckLake integration tests covering the server-gateway model and, separately, a Quack-backed metadata-catalog experiment before any compatibility claim.
- `duckbcftools_*` operations layered on typed columns plus the native writer.
