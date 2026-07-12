# Packed-state kernels

Status: open design for richer packed-state classification. Existing logical kernels and
backend coverage are defined by `src/include/duckhts_simd_kernels.def` and
[`simd_dispatch_matrix.md`](simd_dispatch_matrix.md); this note does not reserve new APIs.

## Shared shape

Nucleotide text, BAM nt16 sequence, and a future packed genotype representation can share one
implementation pattern:

```text
load packed bytes -> classify states -> form class masks -> reduce counts
```

The reusable unit is a typed logical kernel with a scalar oracle. SQL-facing policy remains
above it: which states are called, whether masking changes a denominator, and whether ambiguous
input is null, an error, or a reported class.

This common shape does not imply one generic untyped function. Text bytes, nt16 nibbles, and
two-bit calls have different validity and padding contracts and should remain separate typed
kernels.

## Current boundary

The shipped text GC kernel intentionally preserves `seq_gc_content(...)` behavior: ASCII case
is folded, `A/C/G/T/N` are accepted, and another byte makes the SQL result null. The shipped
BAM nt16 kernel classifies the packed htslib alphabet for `bam_bin_counts(...)`. Do not describe
the text kernel as a total IUPAC classifier, and do not change its public behavior as an
incidental SIMD refactor.

htslib is the authority for sequence encodings. Any text-to-nt16 work must use the linked
`seq_nt16_table`; BAM work uses the nibble values supplied by `bam_get_seq(...)`. A backend must
not embed a divergent private alphabet.

## Rich text classification investigation

A richer internal classifier may be useful for reference summaries or QC that must retain
ambiguous and soft-masked bases. If a measured consumer requires it, the design should:

- partition the full htslib nt16 space into A/T, C/G, N, ambiguity, and `=` facts;
- retain lowercase as an independent mask overlay rather than destroying it during case fold;
- map every input byte through htslib's table before classification;
- report facts only, leaving called/GC/error policy to the caller; and
- preserve current `seq_gc_content(...)` behavior through its policy adapter and regression
  tests.

One unresolved choice is whether callers need to distinguish a literal `N` from an arbitrary
byte that htslib maps to nt16 `N`. That distinction cannot be recovered after table mapping and
must be justified before adding another byte class.

## Packed genotype investigation

Do not add a genotype reduction kernel until a producer defines the packed row contract. That
contract must specify allele ordering, ploidy, missing and partial calls, phasing, padding bits,
sample masks, and multiallelic handling. Only then can allele count, allele number, genotype
class, and missingness reductions have a stable scalar oracle.

The kernel belongs in the SIMD manifest only when its producer, consumer, bounds, and tests are
present. A hypothetical reader name or C struct in a design note is not an interface.

## Validation gate

For any added packed-state kernel:

- scalar and every runnable backend must be bit-identical for every reported class;
- tests cover tails, odd packing widths, invalid values, padding, lowercase, ambiguity, and
  null policy as applicable;
- dispatch is sampled once per vector, chunk, or worker phase; and
- a standalone kernel benchmark and at least one real consumer workload must show where the
  speedup survives transport and materialization cost.
