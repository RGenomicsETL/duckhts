# `duckhts_mosdepth` compatibility contract

Status: current compatibility contract plus open memory and conformance work for the
implemented `duckhts_mosdepth(...)` table function. `functions.yaml`, source, and executable
tests are authoritative for the public signature and implemented options.

## Compatibility target

DuckHTS targets the supported subset of mosdepth 0.3.13, tag `v0.3.13`, commit
`52813e08fe575c39a76556d1fc5c5399a6141e02`. The repository last recorded a local validation
against that target on 2026-04-15.

This pin limits the claim. Do not describe the function as compatible with arbitrary mosdepth
versions, and re-run differential tests after changing coverage, filtering, formatting, or
output code.

## Output and semantic contract

The function is a side-effecting compatibility surface, distinct from SQL-native coverage
relations. For requested outputs, it preserves mosdepth file naming, header order, formatting,
and interval conventions. Plain-text outputs must compare byte-for-byte; BGZF outputs compare
as decompressed text, and their indexes must be query-equivalent rather than compressed-byte
identical.

The supported contract includes these invariants:

- the default excluded SAM flag mask is `1796`;
- default mode is CIGAR-aware and performs mosdepth-style mate-overlap correction;
- fast and fragment modes retain their upstream counting meanings;
- window and BED regions share the same mean, median, threshold, and distribution intervals;
- summaries and distributions retain upstream ordering and decimal formatting; and
- htslib owns BAM/CRAM transport, index loading, CRAM reference handling, and BGZF/CSI I/O.

`threads` controls htslib decompression workers, matching upstream's thread concept.
`processing_threads` is a separate DuckHTS extension that processes contigs concurrently. It
must not change bytes or row order in the final files and must not be described as an upstream
mosdepth option.

Do not fold this surface into a generic coverage function. A compatibility command may share
private readers and kernels, but its filenames, side effects, and differential oracle remain a
separate public contract.

## Ownership and execution

One invocation owns the complete output set. Configuration is fixed before coverage begins;
each processing worker owns its BAM handle, iterator, coverage scratch, and overlap state.
Shared writers are serialized in header order.

Coverage filtering, CIGAR expansion, overlap correction, and fragment geometry are one
semantic pipeline. Refactors may split helpers but must not create a second implementation of
those rules for a particular output type. Summary, distribution, per-base, quantized, and
threshold writers consume the same finalized coverage facts.

Output creation must remain internal to htslib and first-party code. Do not shell out to
`mosdepth`, `bgzip`, or `tabix` from the extension.

## Validation

Automated SQL and R tests cover the shipped API. Upstream comparison is driven by:

```bash
python3 scripts/mosdepth_conformance.py input.bam \
  --extension build/release/duckhts.duckdb_extension

python3 scripts/mosdepth_benchmark.py input.bam \
  --mode default \
  --extension build/release/duckhts.duckdb_extension \
  --verify
```

A compatibility or performance claim must record the input, exact upstream commit, command,
mode, DuckHTS processing/decompression thread counts, validation date, and unsupported
features. Benchmarking the native function against an old SQL reconstruction is not
conformance evidence.

Docs and scientific output should credit mosdepth and its authors when relying on
mosdepth-compatible behavior.

## Open work

- Bound whole-contig coverage scratch when `processing_threads > 1` without changing output;
  see [`coverage_memory_footprint.md`](coverage_memory_footprint.md).
- Refresh differential coverage on broader real BAM and CRAM inputs after engine changes.
- Add another upstream feature or output format only for a demonstrated downstream need, with
  the pin, supported subset, tests, and attribution updated together.
