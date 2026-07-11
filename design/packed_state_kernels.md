# Packed-state classification and reduction kernels

Status: open design / implementation contract. This is the review-before-code gate for extending the SIMD manifest (`src/include/duckhts_simd_kernels.def`). Nothing here is implemented until it lands in the manifest, the backends, `functions.yaml` where public, and tests. The current manifest has exactly one kernel (`seq_base_counts`); this note governs how the next several land coherently instead of as ad-hoc additions.

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

## Caliber: htslib is the sequence-semantics oracle

We do not invent base semantics; htslib owns them. The authority is three tables, **defined in the linked vendored htslib (`third_party/htslib/hts.c`) and exported as symbols via `htslib/hts.h`**:

- `seq_nt16_table[256]` — ASCII byte → nt16 4-bit code. **Case-folded** (lowercase ≡ uppercase), and **every unrepresentable byte folds to `15` (N)** — there is no "invalid" or "other" in htslib's model; `*`, `-`, digits, punctuation, NUL all become N. (`=`→0, `U`→T, `0-3`→A,C,G,T are htslib quirks we inherit by using the table.)
- `seq_nt16_str[] = "=ACMGRSVTWYHKDBN"` — nt16 code → canonical uppercase char.
- `seq_nt16_int[] = {4,0,1,4,2,4,4,4,3,4,4,4,4,4,4,4}` — nt16 code → 2-bit (`A=0,C=1,G=2,T=3`, everything else `=4`). This is htslib's own "is it ACGT" authority: `called` ⇔ `seq_nt16_int[code] != 4`.

Source-of-truth discipline: C code uses these **exported symbols** (`#include "htslib/hts.h"`); it must not hardcode the 256 values and must not reference `r/Rduckhts/inst/duckhts_extension/htslib`, which is bootstrap-generated output, never a source or caliber. The values quoted above are documentation only.

**Mandate: the text classifier is a vectorized replica of htslib's `seq_nt16_table`, not a reinvention, and is derived from the linked symbol at init** (not a baked-in copy), so it tracks whatever htslib version the build links. Each text byte is mapped to its nt16 code first, then classified by code. This guarantees (a) htslib-faithfulness by construction, (b) that the text path and the BAM path — which already starts from nt16 via `bam_get_seq` — agree exactly, and (c) a trivial conformance oracle: the SIMD LUT must reproduce `seq_nt16_table` for all 256 bytes and the per-class reduction must reconcile to `seq_nt16_int` buckets. The current `& 0xDF` fold and abort-on-non-ACGTN are the two existing divergences this removes (htslib never errors; it folds to N).

## Kernel 1 — `seq_base_counts` rewrite (text nucleotides)

Replaces the current scalar/AVX2 5-compare classifier in `src/simd/duckhts_simd_{scalar,avx2,...}.c`.

### Classification

Map each byte through the `seq_nt16_table` replica, then partition the 16 nt16 codes. Classes (all sub-partitions of htslib's space, so totals reconcile to `seq_nt16_int`):

- `gc` — codes `C(2)`, `G(4)`
- `at` — codes `A(1)`, `T(8)`
- `n` — code `15` (true N **and** every byte htslib folds to N; we do not separate them, because htslib cannot either)
- `iupac` — the ambiguity codes `{M3,R5,S6,V7,W9,Y10,H11,K12,D13,B14}` (optionally expose `S`/`W` separately as the only GC-*determinate* pair; the rest are GC-ambiguous)
- `equals` — code `0` (`=`, the equals-reference symbol)
- `masked` — overlay: source byte was lowercase. **The one intentional superset over htslib** (which case-folds); additive only, changes no base-identity count, and reconciles to htslib exactly when ignored.

`called = gc + at` (⇔ `seq_nt16_int != 4`), bit-identical to what samtools/htslib would call ACGT. There is **no `other` class** — htslib has none.

### Output struct (richer; supersedes the 3-field struct)

```c
typedef struct {
    uint64_t gc;       /* nt16 C,G */
    uint64_t at;       /* nt16 A,T */
    uint64_t n;        /* nt16 15 (N + htslib-unrepresentable) */
    uint64_t iupac;    /* nt16 ambiguity codes (not N, not =) */
    uint64_t equals;   /* nt16 0 (=) */
    uint64_t masked;   /* overlay: bytes that were lowercase */
} duckhts_simd_base_class_counts_t;
/* called = gc + at; total = gc+at+n+iupac+equals (== len) */
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

## Performance harness and methodology

Proving a kernel is faster is a first-class deliverable, not an afterthought. The harness is a DuckDB/R-free C microbenchmark, `test/scripts/bench_simd_kernels.c`, run via `make bench-simd-kernels` (built `-O2` to match the CRAN R package; the extension itself ships `-O3` and the win holds at both). It isolates the kernel from query/transport/decode overhead, gates every compiled+cpu-supported backend against the scalar oracle (a divergence fails the run), and reports scalar-relative throughput.

A speedup number is only honest once it survives three checks; record them with the result:

1. **It is vectorization, not the optimizer.** Disassemble (`objdump -d`) and confirm the SIMD kernel actually emits vector instructions and the scalar does not. For `bam_nt16_counts` the AVX2 body shows `vpcmpeqb`/`vpmovmskb`/`popcnt` on `ymm` registers (5 ACGTN codes × 2 nibbles = 10 compares); the scalar shows only `movzbl`/`add`.
2. **It is not an `-O3` unrolling artifact.** Sweep `-O0/-O2/-O3`; the win must be stable at the shipping levels. `bam_nt16_counts` measures ~4.6× at `-O0` (nothing else optimized) but ~10.4× at `-O2` and ~10.7× at `-O3` — stable where it ships.
3. **The scalar baseline is fair.** The oracle must be a real branchless classifier, not a mispredict-bound `switch`; a branchy baseline inflates the ratio (the same kernel showed a spurious ~125× against a `switch` scalar on random input). `bam_nt16_counts` uses a branchless nt16→class histogram.

Recorded result (this host, `-O2`): `bam_nt16_counts` ≈ 1.75 → 18.2 Gbase/s, **~10.4× over the scalar reference**, AVX2 bit-identical to scalar. Numbers are workload- and host-dependent; the microbench isolates the kernel, while the end-to-end win through `bam_bin_counts` is diluted by BAM decode (`benchmark_simd_bam_gc.Rmd`). Do not advertise the microbench multiplier as an end-to-end speedup.

## Out of scope here

- The `read_bcf_gt` reader schema and matrix decode — `genotype_matrix_reader.md` (planned).
- CSQ/INFO delimiter-scan kernel — related idiom (field-split, not class-count); its own future note.
- Alignment kernels (SW/edit distance) — different shape, not packed-state reduction.
