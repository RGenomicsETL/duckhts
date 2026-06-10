# Packed-state classification and reduction kernels

Status: open design / implementation contract. This is the L0/L1 pillar under [`north_star.md`](north_star.md). It is the review-before-code gate for extending the SIMD manifest (`src/include/duckhts_simd_kernels.def`). Nothing here is implemented until it lands in the manifest, the backends, `functions.yaml` where public, and tests. The current manifest has exactly one kernel (`seq_base_counts`); this note governs how the next several land coherently instead of as ad-hoc additions.

## The unifying abstraction

Three seemingly different DuckHTS hot loops are the same problem:

| Domain | Alphabet | Missing/masked state | Transport |
| --- | --- | --- | --- |
| nucleotide text | `A C G T` (+ IUPAC) | `N`, soft-mask (lowercase) | one ASCII byte/base |
| BAM sequence | `A C G T` (+ IUPAC) | `N` (`=15`) | one nt16 nibble/base (2 bases/byte) |
| genotype hardcall | `hom-ref het hom-alt` | `missing` | 2 bits/call (PLINK `.bed` layout) |

All are **a packed small alphabet with a missing/masked state**, and all reduce by the same operation: classify each element into a small set of classes, then **masked popcount per class**. The kernel idiom is therefore one shape across the project:

```
load packed bytes -> classify (table lookup) -> per-class bitmasks -> popcount per class -> accumulate
```

This is portable across every backend we ship:

- table lookup: `vpshufb` (AVX2/AVX-512), `vqtbl1q_u8` (NEON), `i8x16.swizzle` (wasm128), scalar LUT array.
- bitmask extraction: `movemask`/compare, or popcount over packed lanes.
- popcount: `__builtin_popcount*` / `vcnt` / `i8x16.popcnt`.

Two non-negotiable properties, both fixing real defects in today's `base_counts`:

1. **Non-destructive.** No `& 0xDF` case-fold that erases soft-masking. Case (mask state) is a class, not discarded information.
2. **Total, never aborting.** Every byte/element maps to some class (including `other`). A single ambiguous base must never void a whole sequence's counts. The current early-return-on-`invalid` behavior is a footgun at reference/population scale.

The kernel is a **byte/element oracle**: it returns per-class counts. **Policy lives in the SQL layer above** (does soft-mask count as called? do `S`/`W` contribute to GC? is `other` an error?). Scalar is the correctness oracle; `auto` must match it bit-for-bit on every class.

## Kernel 1 — `seq_base_counts` rewrite (text nucleotides)

Replaces the current scalar/AVX2 5-compare classifier in `src/simd/duckhts_simd_{scalar,avx2,...}.c`.

### Classification

ASCII letters classify by **high-nibble gate** (`0x40` upper, `0x60` lower → mask state) plus a **low-nibble LUT** that maps `A..Z` to a class code. Classes:

- `gc` — `C`,`G`
- `at` — `A`,`T`
- `n` — `N`
- `iupac_sw` — `S` (∈{C,G}), `W` (∈{A,T}) — the only GC-*determinate* ambiguity codes
- `iupac_other` — `R Y K M B D H V` — GC-ambiguous ambiguity codes
- `masked` — set additionally when the source byte was lowercase (independent of base class)
- `other` — anything else (digits, `*`, `-`, `.`, junk)

`called = gc + at` (ACGT only). `masked` is an overlay count (how many of the classified bases were soft-masked), not a separate partition, so policy can choose to include or exclude it.

### Output struct (richer; supersedes the 3-field struct)

```c
typedef struct {
    uint64_t gc;            /* C,G (upper+lower) */
    uint64_t at;            /* A,T (upper+lower) */
    uint64_t n;             /* N,n */
    uint64_t iupac_sw;      /* S,W,s,w */
    uint64_t iupac_other;   /* other IUPAC ambiguity codes */
    uint64_t masked;        /* overlay: bases that were lowercase */
    uint64_t other;         /* non-nucleotide bytes */
} duckhts_simd_base_class_counts_t;
```

### Backward compatibility with `seq_gc_content` (mandatory)

`seq_gc_content(...)` is public with existing SQL/R tests. The rewrite must keep its observable output **bit-identical by default**. The current semantics, reconstructed from the new counts in the *policy layer* (not the kernel):

- `gc_current = gc` (soft-mask already folded in, matching `& 0xDF`).
- `called_current = gc + at` (soft-masked ACGT already included).
- `invalid_current = (iupac_sw + iupac_other + other) > 0` — current code flags every non-`ACGTN` byte invalid, and IUPAC is non-`ACGTN`.

One behavioral subtlety to verify against existing tests: the current code returns *partial* counts up to the first invalid byte and sets `invalid=1`; the new kernel returns *full* counts. If `seq_gc_content` surfaces `NULL`/error when invalid (likely), partial-vs-full counts are unobservable and we are safe. If any test asserts the partial count, either reproduce it in the policy layer or treat it as a documented, changelogged change. Confirm before landing.

Any *new* behavior (e.g. "don't treat IUPAC as invalid", "expected-GC counting `S`/`W` as ½", "exclude soft-masked") is exposed as an **opt-in parameter or new function**, separately changelogged in both `NEWS.md` files. The kernel change is internal; public behavior changes are explicit.

## Kernel 2 — `bam_nt16_counts` (BAM packed sequence)

Replaces the scalar `count_read_gc_bases` in `src/bam_bin_counts.c`. Input is `bam_get_seq(rec)` — one nt16 nibble per base (`A=1,C=2,G=4,T=8,N=15`, IUPAC `M=3,R=5,S=6,…`). This is the *purest* form of the idiom: a **16-entry nibble LUT** is literally the native encoding — no ASCII gymnastics. Two bases per byte; process high and low nibble per byte.

```c
typedef struct {
    uint64_t a, c, g, t, n;
    uint64_t iupac;   /* any ambiguity nt16 code that is not N */
    uint64_t gc;      /* c + g */
    uint64_t called;  /* a + c + g + t */
} duckhts_simd_bam_nt16_counts_t;
```

Consumers: `bam_bin_counts(..., stats := 'gc')`, markduplicate, future QC. nt16 has no case/soft-mask, so `masked` is absent here; the class LUT is otherwise the same partition.

## Kernel 3 — genotype reductions (packed 2-bit hardcalls)

Operates on a variant-major packed row (`ceil(n_samples/4)` bytes, 2 bits/sample, states `{hom-ref, het, hom-alt, missing}`) produced by the future `read_bcf_gt` reader (see `genotype_matrix_reader.md`, planned). Reductions are masked popcount over bitplanes:

```c
typedef struct {
    uint64_t hom_ref, het, hom_alt, missing;
    uint64_t ac;   /* 2*hom_alt + het (alt allele count, diploid) */
    uint64_t an;   /* 2*(hom_ref+het+hom_alt) (called alleles) */
} duckhts_simd_geno_counts_t;
```

`ac`/`an`/`missing` over a packed row are the per-variant allele-frequency/missingness primitives. Same kernel family; the "alphabet" is 2-bit states and the LUT is a 256-entry byte→4-state-count table (or bitplane popcount).

## Dispatch / manifest integration

Each new kernel follows the existing typed pattern (`src/include/duckhts_simd_internal.h`), keeping the X-macro table/enum in sync and the per-kernel `consider` helpers small and explicit (TigerStyle: typed boilerplate over a `void*` generic). Per kernel:

1. **`.def` line** — e.g. `DUCKHTS_SIMD_KERNEL(BAM_NT16_COUNTS, bam_nt16_counts, BAM_NT16_COUNTS, "bam_nt16_counts")`. This auto-generates the enum entry, the table function-pointer field, and the name-table row.
2. **fn typedef + `DUCKHTS_SIMD_FN_<TAG>` macro** in `duckhts_simd_internal.h`.
3. **`duckhts_simd_builder_consider_<kernel>(builder, cap, backend, priority, fn)`** helper, mirroring `consider_base_counts`.
4. **backend registration** in each TU that implements it (`*_register`). A backend that lacks a kernel simply does not register; the scalar slot fills it. **Scalar must implement every kernel** — it is the fallback and the oracle.
5. **caller snapshots once per chunk** via `duckhts_simd_dispatch_snapshot()` + `_with_table`, exactly as `kmer_udf.c` already does for `base_counts` — never an atomic load per row.

The `seq_base_counts` enum/field is renamed/extended to carry the richer struct; this is internal (no public manifest name churn unless `duckhts_simd_kernel_info()` rows are asserted — update those tests in lockstep).

## Backend gating and portability

- Scalar: always present, all kernels, the oracle.
- AVX2/AVX-512: `__attribute__((target("avx2,popcnt")))` per-function so the TU stays baseline-safe; `vpshufb` for the LUT. AVX-512 can fold the popcount with `vpopcntb` where available (`avx512-vpopcntdq`) — separate cap, separate priority.
- NEON: `vqtbl1q_u8` + `vcntq_u8`.
- wasm128: `i8x16.swizzle` + `i8x16.popcnt`, gated on `DUCKDB_WASM_EXTENSION` with the existing wasm exception/longjmp model untouched.
- Unknown ISA / musl / RISC-V: scalar only until a tested backend with real probes is added. No libc-specific detection paths.

## Testing contract

Backend-agnostic, scalar-as-oracle, per the dispatch-matrix doc:

- **bit-identical scalar vs `auto`** on every class count, over long deterministic inputs, for all three kernels.
- text fixtures must include: pure ACGT, embedded `N`, **soft-masked lowercase runs**, `S`/`W`, other IUPAC, and trailing junk (`*`,`-`,`.`) — asserting non-destructive case handling and no-abort totality.
- `bam_nt16_counts` validated against the text classifier on the same sequence decoded both ways (nt16 round-trip), and against `count_read_gc_bases` on the GC subtotal.
- genotype reductions validated against a naive scalar AC/AN/missing loop over unpacked calls.
- `seq_gc_content` regression: existing SQL/R tests stay green unchanged (the bit-identical default), proving the rewrite is non-breaking.
- chunk-boundary tests (lengths not multiples of 16/32; odd nt16 base counts; sample counts not multiples of 4).

## Out of scope here

- The `read_bcf_gt` reader schema and matrix decode — `genotype_matrix_reader.md` (planned).
- CSQ/INFO delimiter-scan kernel — related idiom (field-split, not class-count); its own future note.
- Alignment kernels (SW/edit distance) — different shape, not packed-state reduction.