/*
 * duckvep_kernel.h — the duckvep consequence engine ABI.
 *
 * LIBRARY BOUNDARY. This header exposes NO DuckDB, htslib, Arrow, or Parquet
 * types — only immutable borrowed views (struct-of-arrays + explicit lengths)
 * and plain-old-data result records. The DuckDB extension adapter points the
 * view fields directly at flat DuckDB vectors; a standalone CLI, conformance
 * runner, or WASM wrapper points them at ordinary arrays. The engine never
 * owns, frees, or outlives the borrowed memory, and never transfers pointer
 * ownership across this API. That keeps the engine an independently testable,
 * fuzzable library and protects it from DuckDB/htslib API coupling.
 *
 * C discipline enforced throughout the implementation:
 *   - every borrowed array carries an explicit element count;
 *   - all offsets are validated once, when a model is opened;
 *   - 64-bit arithmetic for pool offsets/sizes, with checked narrowing to 32;
 *   - no variable-length arrays, no per-row malloc (per-worker arenas);
 *   - the model is immutable after preparation; no shared mutable active set;
 *   - status-return APIs only — never silent truncation.
 */
#ifndef DUCKVEP_KERNEL_H
#define DUCKVEP_KERNEL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DUCKVEP_KERNEL_VERSION_MAJOR 0
#define DUCKVEP_KERNEL_VERSION_MINOR 5
#define DUCKVEP_KERNEL_VERSION_PATCH 0

/* --------------------------------------------------------------- status -- */

typedef enum duckvep_status {
    DUCKVEP_OK = 0,
    DUCKVEP_ERR_INVALID_ARG,    /* null/zero-length/contract violation at the boundary */
    DUCKVEP_ERR_MODEL_INVALID,  /* offsets/counts inconsistent at model-open           */
    DUCKVEP_ERR_OUT_OF_RANGE,   /* a 64->32 narrowing or index bound check failed      */
    DUCKVEP_ERR_RESULT_FULL,    /* result builder capacity exhausted (no truncation)   */
    DUCKVEP_ERR_UNSUPPORTED,    /* event cannot be represented by this entry point */
    DUCKVEP_ERR_INTERNAL
} duckvep_status_t;

/* Caller-allocated error detail; no ownership transfer. `message` is a
 * NUL-terminated, fixed-capacity buffer the engine writes into on failure. */
typedef struct duckvep_error {
    duckvep_status_t status;
    uint32_t         where_code;   /* stable engine-internal site id (test anchor) */
    char             message[256];
} duckvep_error_t;

/* ------------------------------------------------------------- vocabulary -
 * Small POD enums referenced by the views/results. The full SO-term bitset
 * lives in a separate duckvep_so.h (added with the predicate kernels). */

typedef enum duckvep_impact {   /* mirrors VEP HIGH/MODERATE/LOW/MODIFIER */
    DUCKVEP_IMPACT_MODIFIER = 0,
    DUCKVEP_IMPACT_LOW,
    DUCKVEP_IMPACT_MODERATE,
    DUCKVEP_IMPACT_HIGH
} duckvep_impact_t;

/* Bits in duckvep_consequence_t.region_mask — which transcript compartments an
 * event span overlaps. Strand-resolved (UPSTREAM/DOWNSTREAM are
 * 5'/3' of the transcript). EXON is a non-coding exon; coding exon positions are
 * CDS or UTR. SPLICE is set iff the engine emits a splice consequence for the
 * position (the authoritative strand-resolved splice fact — NOT mere proximity to
 * any exon boundary), and combines with the others. */
typedef enum duckvep_region_bit {
    DUCKVEP_REGION_UPSTREAM   = 1u << 0,
    DUCKVEP_REGION_DOWNSTREAM = 1u << 1,
    DUCKVEP_REGION_INTRON     = 1u << 2,
    DUCKVEP_REGION_EXON       = 1u << 3,
    DUCKVEP_REGION_CDS        = 1u << 4,
    DUCKVEP_REGION_UTR        = 1u << 5,
    DUCKVEP_REGION_SPLICE     = 1u << 6
} duckvep_region_bit_t;

typedef enum duckvep_variant_kind {
    DUCKVEP_KIND_SNV = 0,
    DUCKVEP_KIND_INS,
    DUCKVEP_KIND_DEL,
    DUCKVEP_KIND_INDEL,
    DUCKVEP_KIND_MNV,
    DUCKVEP_KIND_SV         /* symbolic/structural; dispatched to the SV kernel */
} duckvep_variant_kind_t;

/* Structural-event operation. `variant_kind` selects the event family;
 * this subtype preserves the operation required by VEP's structural predicates.
 * CNV direction is deliberately separate: an undirected copy-number interval is
 * not intrinsically a gain or loss without sample/ploidy context. */
typedef enum duckvep_sv_type {
    DUCKVEP_SV_NONE = 0,          /* non-SV row / no structural subtype */
    DUCKVEP_SV_UNKNOWN,           /* structural event with unknown operation */
    DUCKVEP_SV_INSERTION,
    DUCKVEP_SV_DELETION,
    DUCKVEP_SV_DUPLICATION,
    DUCKVEP_SV_TANDEM_DUPLICATION,
    DUCKVEP_SV_INVERSION,
    DUCKVEP_SV_CNV,
    DUCKVEP_SV_BREAKEND
} duckvep_sv_type_t;

typedef enum duckvep_copy_change {
    DUCKVEP_COPY_CHANGE_UNKNOWN = 0,
    DUCKVEP_COPY_CHANGE_LOSS,
    DUCKVEP_COPY_CHANGE_NEUTRAL,
    DUCKVEP_COPY_CHANGE_GAIN
} duckvep_copy_change_t;

/* Validate one single-locus structural operation/copy-direction pair. This is
 * shared by direct borrowed-view callers and adapter-side transport validation,
 * so contradictory compound semantics cannot enter the predicate kernel. BND is
 * metadata-valid here but remains unsupported by the one-interval executor. */
int duckvep_sv_metadata_valid(duckvep_sv_type_t sv_type,
                              duckvep_copy_change_t copy_change);

/* Distilled transcript facts produced by the Ensembl-cache/prepared-model import
 * layer. The kernel is organism-blind: it never switches on raw biotype strings or
 * attrib codes, only these stable bits. Build-time DuckDB SQL (MySQL extension or
 * FTP-dump read_csv) maps Ensembl's transcript/translation_attrib vocabulary into
 * this bitset for each organism/assembly. */
typedef enum duckvep_tx_flag {
    DUCKVEP_TX_HAS_TRANSLATION        = UINT64_C(1) << 0, /* VEP $feat->translation */
    DUCKVEP_TX_BIOTYPE_PROTEIN_CODING = UINT64_C(1) << 1,
    DUCKVEP_TX_BIOTYPE_NMD            = UINT64_C(1) << 2, /* nonsense_mediated_decay */
    DUCKVEP_TX_BIOTYPE_MIRNA          = UINT64_C(1) << 3,
    DUCKVEP_TX_CDS_START_NF           = UINT64_C(1) << 4,
    DUCKVEP_TX_CDS_END_NF             = UINT64_C(1) << 5,
    DUCKVEP_TX_SELENOCYSTEINE         = UINT64_C(1) << 6,
    DUCKVEP_TX_STOP_CODON_READTHROUGH = UINT64_C(1) << 7,
    DUCKVEP_TX_RNA_EDIT               = UINT64_C(1) << 8,
    DUCKVEP_TX_AMINO_ACID_SUB         = UINT64_C(1) << 9,
    DUCKVEP_TX_MANE_SELECT            = UINT64_C(1) << 10,
    DUCKVEP_TX_MANE_PLUS_CLINICAL     = UINT64_C(1) << 11,
    DUCKVEP_TX_GENCODE_BASIC          = UINT64_C(1) << 12,
    DUCKVEP_TX_GENCODE_PRIMARY        = UINT64_C(1) << 13,
    DUCKVEP_TX_CCDS                   = UINT64_C(1) << 14,
    DUCKVEP_TX_READTHROUGH_TRANSCRIPT = UINT64_C(1) << 15,
    DUCKVEP_TX_UPSTREAM_ATG           = UINT64_C(1) << 16
} duckvep_tx_flag_t;

/* ---------------------------------------------------- borrowed input views
 * Coordinates are 1-based, inclusive (the `1` suffix). `chrom_id` is a
 * model-local contig ordinal — the adapter interns names to ids; the kernel
 * never sees a chromosome string. Every pointer is non-owning, has the stated
 * element count, and must outlive the call. */

typedef struct duckvep_variant_batch {
    const uint16_t *chrom_id;      /* [count] model-local contig ordinal      */
    const uint32_t *pos1;          /* [count] 1-based start                    */
    const uint32_t *end1;          /* [count] 1-based inclusive end            */
    const uint32_t *ref_offset;    /* [count] byte offset into allele_bytes    */
    const uint16_t *ref_length;    /* [count]                                  */
    const uint32_t *alt_offset;    /* [count] byte offset into allele_bytes    */
    const uint16_t *alt_length;    /* [count]                                  */
    const uint8_t  *allele_bytes;  /* shared REF/ALT byte pool (A/C/G/T/N/...)  */
    size_t          allele_bytes_len; /* total bytes in allele_bytes (bounds offsets) */
    const uint8_t  *variant_kind;  /* [count] duckvep_variant_kind_t           */
    const uint8_t  *sv_type;       /* OPTIONAL [count] duckvep_sv_type_t; NULL=UNKNOWN for SV */
    const uint8_t  *copy_change;   /* OPTIONAL [count] duckvep_copy_change_t; CNV direction */
    size_t          count;
} duckvep_variant_batch_t;

typedef struct duckvep_transcript_model {
    const uint16_t *chrom_id;      /* [transcript_count]                       */
    const uint32_t *start1;        /* [transcript_count] transcript span start */
    const uint32_t *end1;          /* [transcript_count] transcript span end   */
    const int8_t   *strand;        /* [transcript_count] +1 / -1               */
    const uint64_t *flags;         /* [transcript_count] duckvep_tx_flag_t bits */
    const uint32_t *exon_offset;   /* [transcript_count] index into exon model */
    const uint16_t *exon_count;    /* [transcript_count]                       */
    const uint32_t *cds_start1;    /* [transcript_count] 0 = non-coding        */
    const uint32_t *cds_end1;      /* [transcript_count]                       */
    size_t          transcript_count;
} duckvep_transcript_model_t;

typedef struct duckvep_exon_model {
    const uint32_t *start1;        /* [exon_count] genomic                     */
    const uint32_t *end1;          /* [exon_count] genomic                     */
    const uint32_t *cdna_start1;   /* [exon_count] cDNA-projected, tx order    */
    const uint32_t *cdna_end1;     /* [exon_count]                             */
    const int8_t   *phase;         /* [exon_count] Ensembl phase (-1/0/1/2)     */
    const int8_t   *end_phase;     /* [exon_count] Ensembl end phase           */
    size_t          exon_count;
} duckvep_exon_model_t;

/* Borrowed per-transcript CDS sequence pool. The adapter pre-extracts each
 * transcript's translateable CDS from FASTA at model-prepare time, transcript-
 * strand oriented, with cds_bytes[cds_offset[t]] == CDS position 1 (including any
 * positive Ensembl start-phase padding, often 'N'). This keeps the kernel
 * string-, FASTA-, and allocation-free: it only ever sees A/C/G/T/N bytes. The
 * pool is OPTIONAL — a structural-only or fully non-coding model passes NULL, and
 * a non-coding transcript carries cds_length[t] == 0. */
typedef struct duckvep_sequence_pool {
    const uint8_t  *cds_bytes;     /* shared CDS byte pool, all transcripts concatenated */
    size_t          cds_bytes_len; /* total bytes in cds_bytes (bounds offset+length)    */
    const uint64_t *cds_offset;    /* [transcript_count] byte offset into cds_bytes      */
    const uint32_t *cds_length;    /* [transcript_count] 0 = non-coding                  */
    const uint8_t  *codon_table;   /* [transcript_count] duckvep_codon_table_t per tx    */
    size_t          transcript_count;
} duckvep_sequence_pool_t;

/* ----------------------------------------------------------- compact output
 * The structured consequence result. Compatibility text and HGVS are separate
 * consumers of projected facts. Position fields are -1 when not applicable. */

typedef enum duckvep_consequence_flag {
    /* The pair overlaps CDS, but the sequence context could not resolve a VEP
     * codon/peptide predicate. The consequence_mask may still contain
     * coding_sequence_variant as the safe
     * fallback; this bit makes that fallback auditable rather than silent. */
    DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED = 1u << 0
} duckvep_consequence_flag_t;

typedef struct duckvep_consequence {
    uint32_t variant_idx;
    uint32_t tx_idx;
    uint32_t gene_idx;

    uint64_t consequence_mask;     /* SO terms as a bitset (duckvep_so.h)      */
    uint32_t region_mask;          /* exon/intron/utr/splice/up/down bits      */
    uint32_t flags;

    uint8_t  impact;               /* duckvep_impact_t                         */
    int32_t  cdna_pos;
    int32_t  cds_pos;
    int32_t  protein_pos;

    uint8_t  aa_ref;
    uint8_t  aa_alt;
} duckvep_consequence_t;

/* ------------------------------------------------------------ opaque handles
 * Prepared once by the adapter, then immutable for the lifetime of a run.
 * Their layouts are private to the implementation translation units. */

typedef struct duckvep_model          duckvep_model_t;          /* immutable transcript/exon/seq model */
typedef struct duckvep_options        duckvep_options_t;        /* distance and splice windows          */
typedef struct duckvep_workspace      duckvep_workspace_t;      /* per-worker scratch arena (no malloc/row) */
typedef struct duckvep_annotate_cursor duckvep_annotate_cursor_t; /* resumable tile-local annotator */

/* The result builder is a PUBLIC, stack-allocatable per-worker sink: it hides no
 * invariant beyond "append into a caller-owned fixed-capacity array, never
 * truncate". The caller owns `rows` (capacity entries); the engine appends and
 * returns DUCKVEP_ERR_RESULT_FULL on overflow. `duckvep_annotate_tile` is the
 * one-shot convenience API and is not resumable; use `duckvep_annotate_cursor_*`
 * when the caller needs to pause on output buffer boundaries without recomputing
 * the tile. A row is never partially written. Use duckvep_result_builder_init;
 * do not poke the fields directly. */
typedef struct duckvep_result_builder {
    duckvep_consequence_t *rows;
    size_t                 capacity;
    size_t                 count;
} duckvep_result_builder_t;

/* ------------------------------------------------------------ default tunables
 * The SINGLE place the engine's distance defaults live (VEP's conventions).
 * duckvep_options_open() fills any zero field from these; nothing else in the
 * engine hardcodes them. Override per-call via duckvep_options_init_t. */
#define DUCKVEP_DEFAULT_UPSTREAM_DIST          5000u /* up-/downstream gene window      */
#define DUCKVEP_DEFAULT_DOWNSTREAM_DIST        5000u
#define DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC   3u    /* bases INTO the exon (1-3)        */
#define DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC 8u    /* bases INTO the intron (1-8)      */

/* ----------------------------------------------------------- options init
 * Plain-old-data init descriptor; the opaque duckvep_options_t is prepared from
 * it once. Zero-initialize and override what you need (0 fields take the
 * DUCKVEP_DEFAULT_* values above at open time). */
typedef struct duckvep_options_init {
    uint32_t upstream_dist;          /* 0 -> DUCKVEP_DEFAULT_UPSTREAM_DIST              */
    uint32_t downstream_dist;        /* 0 -> DUCKVEP_DEFAULT_DOWNSTREAM_DIST            */
    uint32_t splice_region_exonic;   /* 0 -> DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC       */
    uint32_t splice_region_intronic; /* 0 -> DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC     */
    uint32_t halo;                   /* 0 -> max(upstream_dist, downstream_dist)        */
} duckvep_options_init_t;

/* ------------------------------------------------------------- entry points
 * The entire engine surface. No file paths, SQL names, DuckDB handles, or
 * bcf1_t — adapters translate those into the views above before calling.
 * Lifecycle: open model + options + a per-worker workspace once, initialize a
 * result builder over caller memory, then annotate each sorted batch. */

/* Returns a static "major.minor.patch" string. Safe before any model open. */
const char *duckvep_kernel_version(void);

/* Prepare an immutable model from borrowed SoA views. Validates ALL offsets and
 * counts ONCE (exon slices in range, cds within transcript span, transcripts
 * sorted ascending by (chrom_id, start1) for the sweep, sequence-pool lengths
 * consistent). `seq` may be NULL (structural-only / non-coding model). On OK,
 * *out_model is a heap handle that BORROWS the view pointers — the caller must
 * keep the underlying arrays alive until duckvep_model_close. */
duckvep_status_t duckvep_model_open(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    duckvep_model_t                 **out_model,
    duckvep_error_t                  *error);
void duckvep_model_close(duckvep_model_t *model);

/* Prepare options from a POD init (NULL init -> all defaults). */
duckvep_status_t duckvep_options_open(
    const duckvep_options_init_t *init,
    duckvep_options_t           **out_options,
    duckvep_error_t              *error);
void duckvep_options_close(duckvep_options_t *options);

/* Open a per-worker scratch arena sized for `model` (sweep active set plus
 * internal sequence-delta scratch). One workspace per thread; never shared
 * concurrently. */
duckvep_status_t duckvep_workspace_open(
    const duckvep_model_t *model,
    duckvep_workspace_t  **out_workspace,
    duckvep_error_t       *error);
void duckvep_workspace_close(duckvep_workspace_t *workspace);

/* Initialize a result builder over caller-owned storage (`rows` holds `capacity`
 * entries). No allocation. */
void   duckvep_result_builder_init(duckvep_result_builder_t *builder,
                                   duckvep_consequence_t    *rows,
                                   size_t                    capacity);
size_t duckvep_result_builder_count(const duckvep_result_builder_t *builder);
void   duckvep_result_builder_reset(duckvep_result_builder_t *builder);

/* Open a tile-local annotation cursor over borrowed views. The cursor validates
 * the same boundary contract as `duckvep_annotate_tile`, then preserves the sweep
 * and candidate position across calls to `duckvep_annotate_cursor_fill`. It copies
 * the small `variants` view struct but borrows the arrays it points at; it does not
 * own `model`, `options`, or `workspace`; keep all borrowed storage alive
 * until close. The cursor is single-threaded and consumes the supplied workspace
 * for its lifetime.
 *
 * `fill` appends rows into a caller-owned result builder until either the cursor is
 * exhausted (DUCKVEP_OK and `duckvep_annotate_cursor_done(cursor) != 0`) or the
 * builder has no remaining capacity (DUCKVEP_ERR_RESULT_FULL, a resumable pause).
 * The caller drains/resets the builder and calls `fill` again. */
duckvep_status_t duckvep_annotate_cursor_open(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    duckvep_annotate_cursor_t     **out_cursor,
    duckvep_error_t                *error);

duckvep_status_t duckvep_annotate_cursor_fill(
    duckvep_annotate_cursor_t      *cursor,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error);

/* Seed the candidate sweep before the first fill. `transcript_indices` is the
 * unique set overlapping the first event's point window. This lets an adapter
 * pay one interval-index lookup per sorted tile, then use the forward sweep. */
duckvep_status_t duckvep_annotate_cursor_seed(
    duckvep_annotate_cursor_t *cursor,
    const uint32_t            *transcript_indices,
    size_t                     transcript_count,
    duckvep_error_t           *error);

int  duckvep_annotate_cursor_done(const duckvep_annotate_cursor_t *cursor);
void duckvep_annotate_cursor_close(duckvep_annotate_cursor_t *cursor);

/* Annotate one tile-local variant batch against an immutable model,
 * emitting duckvep_consequence_t rows into `results`. `workspace` is this
 * worker's scratch.
 *
 * The boundary validates that `variants` is sorted ascending by
 * (chrom_id, pos1), has valid 1-based intervals, and that REF/ALT slices are present,
 * inside `allele_bytes_len`, and consistent with variant kind/raw REF span whenever
 * sequence predicates or allele-trimmed small-variant topology can read them.
 * Contract violations fail before the sweep; they never yield silently wrong
 * candidates or out-of-bounds allele reads. (`model` transcripts
 * are validated sorted at duckvep_model_open.)
 *
 * Thread-safe across workers as long as each owns its own `workspace` + `results`
 * and the model/options are immutable. On failure, writes `error` (if
 * non-NULL) and returns non-OK; on a full result builder returns
 * DUCKVEP_ERR_RESULT_FULL without truncating a partial row. */
duckvep_status_t duckvep_annotate_tile(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_KERNEL_H */
