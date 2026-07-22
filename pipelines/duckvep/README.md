# DuckVEP targets workflow

Status: current optional orchestration for executable corpus campaigns.

This pipeline does not implement another workflow engine or command runner. The ownership
split is:

- `{targets}` owns dependency discovery, file invalidation, dynamic campaign branches,
  resume behavior, saved error workspaces, and the local metadata store;
- `blit` owns shell-script construction and execution of the targets-to-runner and
  micromamba/VEP subprocesses; DuckVEP quotes every dynamic command token before giving
  it to `blit`;
- the DuckVEP conformance scripts own sampling, model/oracle identity, complete pair
  denominators, semantic receipts, and fail-closed acceptance.

Run the checked-in witness campaign from the repository root:

```sh
make duckvep-targets
```

The workflow is optional development infrastructure. Install `{targets}` from CRAN and
the `blit` revision named in
[`test/duckvep/conformance/README.md`](../../test/duckvep/conformance/README.md) before
running it; neither package is linked into the extension or bundled R package.

The pipeline stores paths and receipts, not DuckDB connections, resident DuckVEP models,
or copied corpus objects. It builds the release extension once, binds that build to the
current clean Git revision, and lets all campaign branches reuse the read-only artifact.
Each branch has its own output directory under `.tmp/duckvep_targets/results/`.
The evidence runner hashes the model, FASTA/index, and source VCF it actually consumes
when an executed campaign publishes a new receipt. Those are provenance reads, not a
second invalidation cache; `{targets}` prevents the entire runner from being invoked for
an unchanged branch.
Budget concurrency explicitly: the number of concurrent targets workers multiplies the
per-campaign VEP `fork` and DuckDB `duckdb_threads` values. The checked-in witness uses one
of each, and large-corpus manifests should do the same accounting before enabling a
parallel targets backend.

The external `Rscript`/VEP process is opaque to targets, so its dependency interface is
explicit: a clean source revision, the manifest, the declared local corpus/model/reference
files, the VEP cache metadata marker and acquisition receipt, and the shared extension
receipt. A cache-state target reevaluates the complete installed-file inventory on every
`tar_make()` before downstream evidence is reused.
Adding an input only inside the subprocess without declaring it in the manifest is a
workflow bug.

This is intentionally a thin `{targets}` project, not a DuckVEP workflow engine. DuckVEP
does not implement its own branch scheduler, object store, retry database, timestamp
policy, or digest cache. The evidence runner can still be called directly, but incremental
multi-campaign execution belongs to `{targets}` and external command execution belongs to
`blit`.

Set `DUCKVEP_TARGETS_CAMPAIGNS` to another tab-separated manifest and
`DUCKVEP_TARGETS_OUTPUT` to another result root. Campaign IDs must be unique and path-safe.
Columns use the long option names from `corpus_differential.R` with underscores instead of
hyphens. Boolean fields contain `true` or `false`. Paths may be absolute or relative to the
DuckHTS root. A cache-backed campaign must provide the cache's `cache_info` and
`cache_receipt` files. Create the receipt once after checksum- or
immutable-object-verified acquisition with `make duckvep-cache-receipt`. It binds the exact
species/release/assembly leaf, source identity, `info.txt`, file count, byte total, and a
deterministic relative-path/size/modification/change-time inventory. Each `tar_make()`
restates that inventory, which detects incomplete, replaced, or ordinarily edited shards without rereading tens of
gigabytes of cache contents. Keep the cache in an immutable/read-only store: the metadata
inventory is a practical mutation detector, not an adversarial content hash. The VEP
prefix must likewise be immutable/read-only. The workflow tracks its explicit package
lock and `conda-meta/history`, and the runner verifies the installed package URL set and
uses a clean process environment, but it does not hash every installed package file.

For example, the release-116 GRCh38 archive currently exposes an immutable HTTP ETag and
byte length rather than a cryptographic sidecar:

```sh
make duckvep-cache-receipt DUCKVEP_CACHE_RECEIPT_ARGS="\
  --cache-dir /data/vep-cache \
  --species homo_sapiens --cache-version 116 --assembly GRCh38 \
  --source-url https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz \
  --source-identity http-etag:66fbffe0a-64efc55395049:27644657162 \
  --output /data/receipts/homo_sapiens-116-GRCh38.tsv"
```

If the acquisition source publishes SHA-256, MD5, or BSD `CHECKSUMS`, use
`sha256:HEX`, `md5:HEX`, or `bsd-sum:SUM:BLOCKS` instead. This command inventories an
already installed cache; the download/extraction step remains explicitly outside the
network-free extension build.

`gff_index_policy` defaults to `ignore`, which always restages and indexes the declared
GFF. Set it to `require` only when the manifest intentionally declares the existing `.tbi`
companion. Unlike the standalone runner's `auto` default, targets never changes behavior
merely because an optional index appeared beside an unchanged input.

Branch at a complete campaign or a deliberately coarse corpus shard. Do not create one
target per variant. Carry the campaign ID in every output name rather than relying on
branch order. `{targets}` deliberately chooses the appropriate timestamp-trust behavior
for the filesystem; the script does not override `trust_timestamps`. Large immutable
objects remain where they already live, while `format = "file"` tracks the declared input
and output paths.

`error = "continue"` lets independent campaign branches finish after one branch fails, and
`workspace_on_error = TRUE` makes a failed branch inspectable with `tar_workspace()`. The
checked-in `run.R` wrapper then inspects `tar_errored()` and exits unsuccessfully if any
branch failed. Those features do not relax acceptance: a failed, missing, unresolved, or
discordant campaign cannot satisfy the DuckVEP release-conformance gate.

The campaign-manifest target rejects cross-campaign output collisions, overlap with any
manifest/input/build artifact, output inside the targets store, and arbitrary repo-internal
output roots before producing branches. Each branch carries only its own validated output
set, so changing one campaign does not invalidate unrelated evidence branches. A rerun
writes a complete fresh sibling directory, then replaces the campaign directory with one
atomic rename. An interrupted replacement leaves a discoverable `.previous-*` directory;
the next run restores it when the final directory is absent, or removes it after observing
the complete replacement. The runner reads and writes staging paths during execution but
records the stable published annotation path in the methodology audit.
