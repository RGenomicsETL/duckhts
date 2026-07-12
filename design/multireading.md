# Multi-file reading

Status: current generated-SQL guidance plus open criteria for a possible native multi-file
reader. `functions.yaml` is authoritative for helper signatures and reader parameters.

## Current contract

Use `hts_union_query(...)` in SQL when one reader and one parameter string apply to every
matched file. It generates one branch per file, adds a `filename` column, and combines branches
with `UNION ALL BY NAME`:

```sql
SET VARIABLE q = hts_union_query(
  'read_bcf',
  'cohort/*.vcf.gz',
  'tidy_format := true'
);

SELECT * FROM query(getvariable('q'));
```

The two statements are required because a table macro cannot contain `SET`. An empty match
produces no query string and must be checked before calling `query(...)`.

R users should use the typed multi-file wrappers documented in the generated catalog. Those
wrappers also support explicit file vectors and per-file parameters without asking callers to
construct SQL text.

When files need different SQL arguments, generate the branches from a relation that contains
the file and its parameters. Quote file paths and parameter values as SQL data, not trusted raw
fragments. The built-in helper escapes file paths, but its `reader` and `params` arguments are
SQL-generating inputs and should come from trusted application code.

## Schema rule

Always combine HTS reader branches with `UNION ALL BY NAME`. Fixed-schema inputs also work with
positional unions, but name-based reconciliation is required for VCF/BCF headers that declare
different INFO, FORMAT, or sample columns.

Name matching does not repair semantically incompatible schemas. In particular, positional
tabix columns with the same generated name can mean different things, and columns with one
name but incompatible header-declared types still require an explicit normalization layer.

DuckDB can push projections and filters into individual branches. The cost is one scan node and
one bind per file, so plan construction remains proportional to the number of files.

## Why there is no native glob reader yet

A native reader would improve single-statement ergonomics and might reduce planning overhead,
but it would also own resource and schema behavior that generated SQL currently delegates to
independent reader instances.

The design must resolve:

- a hard bound on simultaneously open local and remote handles;
- index ownership and whether a loaded htslib index is safely shared by worker-local handles;
- one total budget for DuckDB scan workers and htslib decompression workers;
- header probing and `UNION ALL BY NAME`-equivalent schema reconciliation;
- local and object-store glob semantics through APIs available at the stable C boundary;
- failure behavior after some files have emitted rows;
- per-file parameter expressiveness; and
- deterministic source identification and requested-file ordering.

Do not model work as every `(file, contig, thread)` combination with a permanently open handle
and index. That grows file descriptors, index memory, and decompression pools multiplicatively.
A viable design needs a bounded work queue and worker-local handles that are opened, reused, and
closed under an explicit policy.

## Decision gate

Keep generated SQL as the default until a measured workload shows that bind/plan cost or
resource scheduling, rather than I/O and decompression, is the limiting factor. A native reader
must demonstrate a bound and a semantic advantage; removing the two-statement SQL ceremony
alone does not justify duplicating every existing reader's bind and schema logic.
