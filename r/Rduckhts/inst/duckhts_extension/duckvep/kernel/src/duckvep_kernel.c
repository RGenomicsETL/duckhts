/*
 * duckvep_kernel.c — engine entry points: lifecycle constructors + the fused
 * annotate_tile path.
 *
 * This translation unit (and every other under kernel/) links against NO
 * DuckDB, htslib, Arrow, or Parquet symbol. It sees only borrowed typed arrays
 * (duckvep_kernel.h) and produces compact result rows. annotate_tile fuses the
 * already-oracle-tested predicate kernels: the sorted multi-track sweep
 * (duckvep_sweep) generates candidate (variant, transcript) pairs, and the
 * structural classifier (duckvep_classify) places each pair relative to the
 * transcript; sequence-backed coding and structural/CNV predicate producers feed
 * the generated consequence program. HGVS and general edit-set semantics remain
 * late, separate layers.
 */
#include "duckvep_kernel.h"

#include "duckvep_classify.h"
#include "duckvep_delta.h"
#include "duckvep_effect.h"
#include "duckvep_event.h"
#include "duckvep_so.h"
#include "duckvep_sweep.h"
#include "duckvep_workspace_internal.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DUCKVEP_STRINGIFY_INNER(value) #value
#define DUCKVEP_STRINGIFY(value) DUCKVEP_STRINGIFY_INNER(value)

const char *duckvep_kernel_version(void) {
    return DUCKVEP_STRINGIFY(DUCKVEP_KERNEL_VERSION_MAJOR) "."
           DUCKVEP_STRINGIFY(DUCKVEP_KERNEL_VERSION_MINOR) "."
           DUCKVEP_STRINGIFY(DUCKVEP_KERNEL_VERSION_PATCH);
}

static duckvep_status_t fail(duckvep_error_t *error,
                             duckvep_status_t status,
                             uint32_t where_code,
                             const char *message) {
    if (error != NULL) {
        error->status = status;
        error->where_code = where_code;
        /* Bounded copy; the buffer is fixed-capacity by contract. */
        (void)snprintf(error->message, sizeof error->message, "%s",
                       message != NULL ? message : "");
    }
    return status;
}

/* ------------------------------------------------------------------- model --
 * Immutable, prepared once. Borrows the caller's SoA view pointers (it copies
 * the small view structs by value, but never the underlying arrays). */
struct duckvep_model {
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t       exons;
    duckvep_sequence_pool_t    seq;
    uint8_t                   *point_ordered;
    size_t                     max_transcripts_per_chrom;
    size_t                     max_cds_len;
    int                        has_seq;
};

static int transcript_point_ordered(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx) {

    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    size_t cnt = (size_t)transcripts->exon_count[tx_idx];
    size_t e;

    for (e = 0u; e < cnt; e++) {
        size_t ei = off + e;
        uint32_t es = exons->start1[ei];
        uint32_t ee = exons->end1[ei];

        if (es > ee || es < transcripts->start1[tx_idx] ||
            ee > transcripts->end1[tx_idx])
            return 0;
        if (e == 0u) continue;
        if (transcripts->strand[tx_idx] >= 0) {
            if (exons->end1[ei - 1u] >= es) return 0;
        } else {
            if (ee >= exons->start1[ei - 1u]) return 0;
        }
    }
    return 1;
}

/* where_codes are STABLE engine-internal failure-site ids; tests anchor on them.
 * Never renumber an existing site; append new ones. */
enum {
    DVW_MODEL_NULL_OUT      = 10u,
    DVW_MODEL_NULL_VIEW     = 11u,
    DVW_MODEL_OOM           = 12u,
    DVW_MODEL_EXON_RANGE    = 13u,
    DVW_MODEL_CDS_RANGE     = 14u,
    DVW_MODEL_UNSORTED      = 15u,
    DVW_MODEL_SEQ_COUNT     = 16u,
    DVW_MODEL_SEQ_RANGE     = 17u,
    DVW_MODEL_TX_COUNT      = 18u,
    DVW_OPTIONS_NULL_OUT    = 20u,
    DVW_OPTIONS_OOM         = 21u,
    DVW_WS_NULL_ARG         = 30u,
    DVW_WS_OOM              = 31u,
    DVW_ANN_NULL_MODEL      = 40u,
    DVW_ANN_NULL_VARIANTS   = 41u,
    DVW_ANN_NULL_COLUMNS    = 42u,
    DVW_ANN_NULL_OPTIONS    = 43u,
    DVW_ANN_NULL_WORKSPACE  = 44u,
    DVW_ANN_NULL_RESULTS    = 45u,
    DVW_ANN_RESULT_FULL     = 46u,
    DVW_ANN_SWEEP           = 47u,
    DVW_ANN_UNSORTED        = 48u,
    DVW_ANN_COORD_RANGE     = 49u,
    DVW_ANN_ALLELE_RANGE    = 50u,
    DVW_ANN_VARIANT_COUNT   = 51u,
    DVW_ANN_RESULT_STORAGE  = 52u,
    DVW_ANN_RESULT_COUNT    = 53u,
    DVW_ANN_SV_TYPE         = 54u,
    DVW_CURSOR_NULL_OUT     = 60u,
    DVW_CURSOR_OOM          = 61u,
    DVW_CURSOR_NULL         = 62u,
    DVW_WS_SCRATCH_RANGE    = 63u,
    DVW_WS_SCRATCH_OOM      = 64u,
    DVW_WS_MODEL            = 65u
};

duckvep_status_t duckvep_model_open(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    duckvep_model_t                 **out_model,
    duckvep_error_t                  *error) {

    struct duckvep_model *m;
    size_t t;
    size_t chrom_run = 0u;
    size_t max_chrom_run = 0u;
    size_t max_cds_len = 0u;

    if (out_model == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_OUT,
                    "out_model is NULL");
    }
    *out_model = NULL;
    if (transcripts == NULL || exons == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_VIEW,
                    "transcripts or exons view is NULL");
    }
    if (transcripts->transcript_count > (size_t)UINT32_MAX) {
        return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_MODEL_TX_COUNT,
                    "transcript_count exceeds uint32 tx_idx space");
    }
    if (transcripts->transcript_count > 0u &&
        (transcripts->chrom_id == NULL || transcripts->start1 == NULL ||
         transcripts->end1 == NULL || transcripts->strand == NULL ||
         transcripts->flags == NULL ||
         transcripts->exon_offset == NULL || transcripts->exon_count == NULL ||
         transcripts->cds_start1 == NULL || transcripts->cds_end1 == NULL)) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_VIEW,
                    "transcript view has count>0 but null columns");
    }
    /* If any transcript references exons, the exon arrays the classifier reads
     * must be present. */
    if (exons->exon_count > 0u && (exons->start1 == NULL || exons->end1 == NULL)) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_VIEW,
                    "exon view has count>0 but null start/end columns");
    }

    /* Validate every transcript once: span ordering, exon slice in range, cds
     * within span, and the (chrom_id, start1) sort order the sweep relies on. */
    for (t = 0u; t < transcripts->transcript_count; t++) {
        size_t eoff = (size_t)transcripts->exon_offset[t];
        size_t ecnt = (size_t)transcripts->exon_count[t];
        uint32_t cds_s = transcripts->cds_start1[t];
        uint32_t cds_e = transcripts->cds_end1[t];

        if (transcripts->start1[t] > transcripts->end1[t]) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_CDS_RANGE,
                        "transcript start1 > end1");
        }
        if (eoff > exons->exon_count || ecnt > exons->exon_count - eoff) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_EXON_RANGE,
                        "exon slice out of range for a transcript");
        }
        if (cds_s != 0u) {
            if (cds_s > cds_e || cds_s < transcripts->start1[t] ||
                cds_e > transcripts->end1[t]) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_CDS_RANGE,
                            "cds interval outside the transcript span");
            }
        } else if (cds_e != 0u) {
            /* cds_start1 == 0 is the non-coding sentinel; a stray cds_end1 is malformed. */
            return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_CDS_RANGE,
                        "cds_start1 == 0 (non-coding) but cds_end1 != 0");
        }
        if (t > 0u) {
            uint16_t pc = transcripts->chrom_id[t - 1u];
            uint16_t cc = transcripts->chrom_id[t];
            if (cc < pc ||
                (cc == pc && transcripts->start1[t] < transcripts->start1[t - 1u])) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_UNSORTED,
                            "transcripts not sorted ascending by (chrom_id, start1)");
            }
            if (cc == pc) {
                chrom_run++;
            } else {
                if (chrom_run > max_chrom_run) max_chrom_run = chrom_run;
                chrom_run = 1u;
            }
        } else {
            chrom_run = 1u;
        }
    }
    if (chrom_run > max_chrom_run) max_chrom_run = chrom_run;

    if (seq != NULL) {
        if (seq->transcript_count != transcripts->transcript_count) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_SEQ_COUNT,
                        "sequence pool transcript_count mismatch");
        }
        if (seq->transcript_count > 0u &&
            (seq->cds_offset == NULL || seq->cds_length == NULL ||
             seq->codon_table == NULL || (seq->cds_bytes_len > 0u && seq->cds_bytes == NULL))) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_VIEW,
                        "sequence pool has count>0 but null columns");
        }
        for (t = 0u; t < seq->transcript_count; t++) {
            uint64_t off = seq->cds_offset[t];
            uint64_t len = (uint64_t)seq->cds_length[t];
            if (off > (uint64_t)seq->cds_bytes_len ||
                len > (uint64_t)seq->cds_bytes_len - off) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_SEQ_RANGE,
                            "cds offset/length out of the sequence pool");
            }
            if ((size_t)len > max_cds_len) max_cds_len = (size_t)len;
        }
    }

    m = (struct duckvep_model *)calloc(1u, sizeof *m);
    if (m == NULL) {
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_MODEL_OOM, "model alloc failed");
    }
    m->transcripts = *transcripts;
    m->exons = *exons;
    if (transcripts->transcript_count > 0u) {
        m->point_ordered = (uint8_t *)calloc(transcripts->transcript_count,
                                              sizeof *m->point_ordered);
        if (m->point_ordered == NULL) {
            free(m);
            return fail(error, DUCKVEP_ERR_INTERNAL, DVW_MODEL_OOM,
                        "model point-layout alloc failed");
        }
        for (t = 0u; t < transcripts->transcript_count; t++)
            m->point_ordered[t] = (uint8_t)transcript_point_ordered(
                transcripts, exons, t);
    }
    m->max_transcripts_per_chrom = max_chrom_run;
    m->max_cds_len = max_cds_len;
    if (seq != NULL) {
        m->seq = *seq;
        m->has_seq = 1;
    }
    *out_model = m;
    return DUCKVEP_OK;
}

void duckvep_model_close(duckvep_model_t *model) {
    if (model != NULL) {
        free(model->point_ordered);
        free(model);
    }
}

/* ----------------------------------------------------------------- options --
 * Resolved defaults + a precomputed sweep halo (max of the up/downstream reach,
 * since the sweep window is symmetric and the directional cut happens during
 * classification). */
struct duckvep_options {
    uint32_t upstream_dist;
    uint32_t downstream_dist;
    uint32_t splice_region_exonic;
    uint32_t splice_region_intronic;
    uint32_t halo;
};

duckvep_status_t duckvep_options_open(
    const duckvep_options_init_t *init,
    duckvep_options_t           **out_options,
    duckvep_error_t              *error) {

    struct duckvep_options *o;

    if (out_options == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_OPTIONS_NULL_OUT,
                    "out_options is NULL");
    }
    *out_options = NULL;
    o = (struct duckvep_options *)calloc(1u, sizeof *o);
    if (o == NULL) {
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_OPTIONS_OOM, "options alloc failed");
    }

    o->upstream_dist          = (init != NULL && init->upstream_dist != 0u)          ? init->upstream_dist          : DUCKVEP_DEFAULT_UPSTREAM_DIST;
    o->downstream_dist        = (init != NULL && init->downstream_dist != 0u)        ? init->downstream_dist        : DUCKVEP_DEFAULT_DOWNSTREAM_DIST;
    o->splice_region_exonic   = (init != NULL && init->splice_region_exonic != 0u)   ? init->splice_region_exonic   : DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC;
    o->splice_region_intronic = (init != NULL && init->splice_region_intronic != 0u) ? init->splice_region_intronic : DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC;
    {
        /* The sweep window must cover the directional reach, or upstream/downstream
         * candidates beyond `halo` are dropped before the directional filter can
         * keep them. So halo is never smaller than max(upstream, downstream). */
        uint32_t base = o->upstream_dist > o->downstream_dist ? o->upstream_dist : o->downstream_dist;
        uint32_t want = (init != NULL && init->halo != 0u) ? init->halo : base;
        o->halo = want > base ? want : base;
    }

    *out_options = o;
    return DUCKVEP_OK;
}

void duckvep_options_close(duckvep_options_t *options) {
    free(options);
}

/* --------------------------------------------------------------- workspace --
 * Per-worker scratch. The point active set and current-span candidate slice are
 * each sized to the largest prepared contig transcript run. The delta scratch is
 * sized once from the model's maximum CDS length plus the maximum uint16_t small-
 * variant payload, enough for current single-variant MNV, shortening deletion, and
 * lengthening insertion CodingContext builds. Broader haplotype/grouped edit-set
 * paths will widen this policy explicitly when they are wired into production. */
struct duckvep_workspace {
    const duckvep_model_t    *model;
    uint32_t                *active;
    uint32_t                *candidates;
    uint16_t                *point_exon_rank;
    size_t                   point_exon_count;
    size_t                   active_cap;
    uint16_t                 point_last_chrom;
    uint32_t                 point_last_pos;
    int                      point_run_active;
    duckvep_delta_scratch_t                 delta_scratch;
    duckvep_workspace_delta_route_stats_t   delta_route_stats;
    int                                     delta_route_stats_enabled;
};

static void workspace_point_cursor_reset(duckvep_workspace_t *workspace) {
    if (workspace == NULL || workspace->point_exon_rank == NULL) return;
    memset(workspace->point_exon_rank, 0xff,
           workspace->point_exon_count * sizeof *workspace->point_exon_rank);
}

static void workspace_point_run_begin(
    duckvep_workspace_t           *workspace,
    const duckvep_variant_batch_t *variants) {

    uint16_t first_chrom;
    uint32_t first_pos;
    size_t last;

    if (workspace == NULL || variants == NULL || variants->count == 0u) return;
    first_chrom = variants->chrom_id[0];
    first_pos = variants->pos1[0];
    if (workspace->point_run_active &&
        (first_chrom < workspace->point_last_chrom ||
         (first_chrom == workspace->point_last_chrom &&
          first_pos < workspace->point_last_pos))) {
        workspace_point_cursor_reset(workspace);
    }
    last = variants->count - 1u;
    workspace->point_last_chrom = variants->chrom_id[last];
    workspace->point_last_pos = variants->pos1[last];
    workspace->point_run_active = 1;
}

static int size_add_checked(size_t a, size_t b, size_t *out) {
    if (out == NULL) return 0;
    if (b > SIZE_MAX - a) return 0;
    *out = a + b;
    return 1;
}

static int size_mul_checked(size_t a, size_t b, size_t *out) {
    if (out == NULL) return 0;
    if (a != 0u && b > SIZE_MAX / a) return 0;
    *out = a * b;
    return 1;
}

static void workspace_delta_scratch_free(duckvep_delta_scratch_t *s) {
    if (s == NULL) return;
    free(s->edits);
    free(s->alt_cds);
    free(s->ref_peptide);
    free(s->alt_peptide);
    memset(s, 0, sizeof *s);
}

static duckvep_status_t workspace_delta_scratch_open(
    const duckvep_model_t      *model,
    duckvep_delta_scratch_t    *s,
    duckvep_error_t            *error) {

    size_t max_cds;
    size_t max_alt_cds;
    size_t mnv_span;
    size_t bytes;

    if (model == NULL || s == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_WS_NULL_ARG,
                    "workspace delta scratch args are NULL");
    }
    memset(s, 0, sizeof *s);
    max_cds = model->max_cds_len;
    if (max_cds == 0u) return DUCKVEP_OK;

    if (!size_add_checked(max_cds, (size_t)UINT16_MAX, &max_alt_cds)) {
        memset(s, 0, sizeof *s);
        return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_WS_SCRATCH_RANGE,
                    "workspace delta scratch capacity overflow");
    }
    s->alt_cds_cap = max_alt_cds;
    s->ref_peptide_cap = max_cds / 3u + 1u; /* NUL room for full-CDS peptide. */
    s->alt_peptide_cap = max_alt_cds / 3u + 1u;
    mnv_span = max_cds < (size_t)UINT16_MAX ? max_cds : (size_t)UINT16_MAX;
    s->edits_cap = (mnv_span + 1u) / 2u; /* worst-case alternating MNV islands. */

    if (!size_mul_checked(s->edits_cap, sizeof *s->edits, &bytes) ||
        !size_mul_checked(s->alt_cds_cap, sizeof *s->alt_cds, &bytes) ||
        !size_mul_checked(s->ref_peptide_cap, sizeof *s->ref_peptide, &bytes) ||
        !size_mul_checked(s->alt_peptide_cap, sizeof *s->alt_peptide, &bytes)) {
        memset(s, 0, sizeof *s);
        return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_WS_SCRATCH_RANGE,
                    "workspace delta scratch capacity overflow");
    }

    s->edits = (duckvep_haplotype_edit_t *)calloc(s->edits_cap, sizeof *s->edits);
    s->alt_cds = (uint8_t *)calloc(s->alt_cds_cap, sizeof *s->alt_cds);
    s->ref_peptide = (uint8_t *)calloc(s->ref_peptide_cap, sizeof *s->ref_peptide);
    s->alt_peptide = (uint8_t *)calloc(s->alt_peptide_cap, sizeof *s->alt_peptide);
    if (s->edits == NULL || s->alt_cds == NULL ||
        s->ref_peptide == NULL || s->alt_peptide == NULL) {
        workspace_delta_scratch_free(s);
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_WS_SCRATCH_OOM,
                    "workspace delta scratch alloc failed");
    }
    (void)bytes;
    return DUCKVEP_OK;
}

duckvep_status_t duckvep_workspace_open(
    const duckvep_model_t *model,
    duckvep_workspace_t  **out_workspace,
    duckvep_error_t       *error) {

    struct duckvep_workspace *w;
    duckvep_status_t st;
    size_t cap;
    size_t point_bytes = 0u;

    if (out_workspace == NULL || model == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_WS_NULL_ARG,
                    "out_workspace or model is NULL");
    }
    *out_workspace = NULL;
    /* The sweep active set is chromosome-local. Sizing to the largest contig
     * slice preserves the hard bound without allocating one index per transcript
     * across the whole model. */
    cap = model->max_transcripts_per_chrom;
    if (cap == 0u) cap = 1u; /* always have a non-NULL active buffer */
    if (model->transcripts.transcript_count > 0u &&
        !size_mul_checked(model->transcripts.transcript_count,
                          sizeof *w->point_exon_rank, &point_bytes)) {
        return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_WS_SCRATCH_RANGE,
                    "workspace point cursor capacity overflow");
    }

    w = (struct duckvep_workspace *)calloc(1u, sizeof *w);
    if (w == NULL) {
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_WS_OOM, "workspace alloc failed");
    }
    w->active = (uint32_t *)calloc(cap, sizeof *w->active);
    w->candidates = (uint32_t *)calloc(cap, sizeof *w->candidates);
    if (model->transcripts.transcript_count > 0u) {
        w->point_exon_rank = (uint16_t *)malloc(point_bytes);
    }
    if (w->active == NULL || w->candidates == NULL ||
        (model->transcripts.transcript_count > 0u &&
         w->point_exon_rank == NULL)) {
        free(w->active);
        free(w->candidates);
        free(w->point_exon_rank);
        free(w);
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_WS_OOM,
                    "workspace sweep scratch alloc failed");
    }
    st = workspace_delta_scratch_open(model, &w->delta_scratch, error);
    if (st != DUCKVEP_OK) {
        free(w->active);
        free(w->candidates);
        free(w->point_exon_rank);
        free(w);
        return st;
    }
    w->model = model;
    w->point_exon_count = model->transcripts.transcript_count;
    workspace_point_cursor_reset(w);
    w->active_cap = cap;
    *out_workspace = w;
    return DUCKVEP_OK;
}

void duckvep_workspace_close(duckvep_workspace_t *workspace) {
    if (workspace != NULL) {
        free(workspace->active);
        free(workspace->candidates);
        free(workspace->point_exon_rank);
        workspace_delta_scratch_free(&workspace->delta_scratch);
        free(workspace);
    }
}

DUCKVEP_INTERNAL_API duckvep_delta_scratch_t *duckvep_workspace_delta_scratch(
    duckvep_workspace_t *workspace) {

    return workspace != NULL ? &workspace->delta_scratch : NULL;
}

DUCKVEP_INTERNAL_API void duckvep_workspace_delta_route_stats_reset(duckvep_workspace_t *workspace) {
    if (workspace == NULL) return;
    memset(&workspace->delta_route_stats, 0, sizeof workspace->delta_route_stats);
    workspace->delta_route_stats_enabled = 1;
}

DUCKVEP_INTERNAL_API const duckvep_workspace_delta_route_stats_t *duckvep_workspace_delta_route_stats(
    const duckvep_workspace_t *workspace) {

    return workspace != NULL ? &workspace->delta_route_stats : NULL;
}

/* ----------------------------------------------------------- result builder */
void duckvep_result_builder_init(duckvep_result_builder_t *builder,
                                 duckvep_consequence_t    *rows,
                                 size_t                    capacity) {
    if (builder == NULL) return;
    builder->rows = rows;
    builder->capacity = capacity;
    builder->count = 0u;
}

size_t duckvep_result_builder_count(const duckvep_result_builder_t *builder) {
    return builder != NULL ? builder->count : 0u;
}

void duckvep_result_builder_reset(duckvep_result_builder_t *builder) {
    if (builder != NULL) builder->count = 0u;
}

/* Per-candidate sink context for the sweep. */
struct annotate_ctx {
    const duckvep_model_t         *model;
    const duckvep_variant_batch_t *variants;
    const duckvep_options_t       *options;
    duckvep_workspace_t           *workspace;
    duckvep_delta_scratch_t       *delta_scratch;
    duckvep_workspace_delta_route_stats_t *delta_route_stats;
    duckvep_result_builder_t      *results;
    duckvep_status_t               status; /* first non-OK halts emission */
};

/* The VEP-shaped per-candidate decision: fill the cheap facts (effect ctx),
 * escalate to the sequence delta ONLY for the CDS bucket (lazy), then evaluate the
 * static rule table. No biological special-casing lives here — every SO decision is
 * a rule in duckvep_effect.c. */
static int annotate_pair(uint32_t variant_idx, uint32_t tx_idx, void *vctx) {
    struct annotate_ctx *c = (struct annotate_ctx *)vctx;
    const duckvep_transcript_model_t *tx;
    uint32_t pos;
    int fwd;
    uint32_t dist;
    duckvep_variant_kind_t kind;
    duckvep_event_t event;
    duckvep_effect_ctx_t ectx;
    duckvep_sequence_delta_t delta;
    uint64_t cmask;
    duckvep_consequence_t *row;
    int cds_delta_attempted = 0;

    if (c->status != DUCKVEP_OK) return 0;

    tx = &c->model->transcripts;
    duckvep_event_load(c->variants, (size_t)variant_idx, &event);
    pos = event.raw_start1;
    fwd = tx->strand[tx_idx] >= 0;
    kind = (duckvep_variant_kind_t)event.kind;

    if (kind == DUCKVEP_KIND_SNV && c->model->point_ordered[tx_idx]) {
        duckvep_effect_ctx_fill_point_sorted(
            tx, &c->model->exons, variant_idx, (size_t)tx_idx,
            event.start1, c->options->splice_region_exonic,
            c->options->splice_region_intronic,
            &c->workspace->point_exon_rank[tx_idx], &ectx);
    } else {
        duckvep_effect_ctx_fill(tx, &c->model->exons, variant_idx,
                                (size_t)tx_idx, event.start1, event.end1,
                                event.interbase,
                                c->options->splice_region_exonic,
                                c->options->splice_region_intronic, &ectx);
    }
    duckvep_effect_ctx_apply_event(&ectx, &event);
    if (kind == DUCKVEP_KIND_SV) {
        duckvep_sv_effect_t sv = duckvep_sv_effect_fill(&event, &ectx.region_state);
        duckvep_effect_ctx_apply_sv(&ectx, &sv);
    }

    /* The symmetric sweep halo admits candidates up to max(up,down); enforce the
     * directional up/downstream window here so an asymmetric config is honored. */
    if (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM)) {
        dist = fwd ? tx->start1[tx_idx] - event.end1
                   : event.start1 - tx->end1[tx_idx];
        if (dist > c->options->upstream_dist) return 1;
    } else if (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM)) {
        dist = fwd ? event.start1 - tx->end1[tx_idx]
                   : tx->start1[tx_idx] - event.end1;
        if (dist > c->options->downstream_dist) return 1;
    }

    /* Only the CDS bucket needs a sequence delta. The dispatcher takes borrowed
     * model views directly; the annotation sweep does not classify coding effects. */
    memset(&delta, 0, sizeof delta);
    if (kind != DUCKVEP_KIND_SV &&
        (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_CDS))) {
        const duckvep_sequence_pool_t *seq = c->model->has_seq ? &c->model->seq : NULL;
        cds_delta_attempted = 1;
        duckvep_sequence_delta_route_t route = DUCKVEP_DELTA_ROUTE_DIRECT;
        duckvep_sequence_delta_fill_for_annotation_trace(kind, tx, &c->model->exons,
                                                         seq, c->variants, variant_idx,
                                                         (size_t)tx_idx, pos,
                                                         tx->strand[tx_idx],
                                                         c->delta_scratch,
                                                         c->delta_route_stats != NULL ? &route : NULL,
                                                         &delta);
        if (c->delta_route_stats != NULL) {
            if (route == DUCKVEP_DELTA_ROUTE_MNV_CONTEXT) {
                c->delta_route_stats->mnv_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_MNV_DIRECT_FALLBACK) {
                c->delta_route_stats->mnv_direct_fallback++;
            } else if (route == DUCKVEP_DELTA_ROUTE_DEL_CONTEXT) {
                c->delta_route_stats->del_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_DEL_DIRECT_FALLBACK) {
                c->delta_route_stats->del_direct_fallback++;
            } else if (route == DUCKVEP_DELTA_ROUTE_INS_CONTEXT) {
                c->delta_route_stats->ins_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_INS_DIRECT_FALLBACK) {
                c->delta_route_stats->ins_direct_fallback++;
            } else if (route == DUCKVEP_DELTA_ROUTE_INDEL_CONTEXT) {
                c->delta_route_stats->indel_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_INDEL_DIRECT_FALLBACK) {
                c->delta_route_stats->indel_direct_fallback++;
            }
        }
        duckvep_effect_ctx_apply_delta(&ectx, &delta);
    }

    duckvep_effect_ctx_finalize(&ectx);

    cmask = duckvep_effect_eval(ectx.pre_bits);
    if (cmask == 0u) return 1; /* nothing to emit for this pair */

    if (c->results->count >= c->results->capacity) {
        c->status = DUCKVEP_ERR_RESULT_FULL;
        return 0;
    }
    row = &c->results->rows[c->results->count++];
    memset(row, 0, sizeof *row);
    row->variant_idx = variant_idx;
    row->tx_idx = tx_idx;
    row->gene_idx = 0u; /* gene grouping is an adapter concern for now */
    row->consequence_mask = cmask;
    if (cds_delta_attempted && !delta.valid) {
        row->flags |= (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED;
    }
    row->region_mask = ectx.region;
    row->impact = (uint8_t)duckvep_so_impact(cmask);
    if (delta.valid) {
        row->cdna_pos = delta.cdna_pos;
        row->cds_pos = delta.cds_pos;
        row->protein_pos = delta.protein_pos;
        row->aa_ref = delta.ref_aa;
        row->aa_alt = delta.alt_aa;
    } else {
        row->cdna_pos = -1;
        row->cds_pos = -1;
        row->protein_pos = -1;
    }
    return 1;
}

static int allele_shape_matches_kind(uint8_t        kind,
                                     uint32_t       pos1,
                                     uint32_t       end1,
                                     const uint8_t *ref,
                                     uint16_t       ref_len,
                                     const uint8_t *alt,
                                     uint16_t       alt_len) {
    uint16_t prefix = 0u;
    uint16_t suffix = 0u;
    uint16_t ref_rem;
    uint16_t alt_rem;
    uint16_t ref_diff_len;
    uint16_t alt_diff_len;
    uint64_t span;

    if (ref_len == 0u || alt_len == 0u || ref == NULL || alt == NULL) return 0;
    span = (uint64_t)ref_len - 1u;
    if ((uint64_t)pos1 + span > UINT32_MAX ||
        (uint64_t)end1 != (uint64_t)pos1 + span) {
        return 0;
    }

    while (prefix < ref_len && prefix < alt_len && ref[prefix] == alt[prefix]) {
        prefix++;
    }
    ref_rem = (uint16_t)(ref_len - prefix);
    alt_rem = (uint16_t)(alt_len - prefix);
    while (suffix < ref_rem && suffix < alt_rem &&
           ref[(uint16_t)(ref_len - 1u - suffix)] ==
           alt[(uint16_t)(alt_len - 1u - suffix)]) {
        suffix++;
    }
    ref_diff_len = (uint16_t)(ref_len - prefix - suffix);
    alt_diff_len = (uint16_t)(alt_len - prefix - suffix);

    switch ((duckvep_variant_kind_t)kind) {
    case DUCKVEP_KIND_SNV:
        return ref_len == 1u && alt_len == 1u && ref_diff_len == 1u &&
               alt_diff_len == 1u;
    case DUCKVEP_KIND_INS:
        return prefix > 0u && ref_diff_len == 0u && alt_diff_len > 0u;
    case DUCKVEP_KIND_DEL:
        return ref_diff_len > 0u && alt_diff_len == 0u;
    case DUCKVEP_KIND_MNV:
        return ref_len == alt_len && ref_len >= 2u && ref_diff_len > 0u &&
               ref_diff_len == alt_diff_len;
    case DUCKVEP_KIND_INDEL:
        return ref_diff_len > 0u && alt_diff_len > 0u &&
               ref_diff_len != alt_diff_len;
    case DUCKVEP_KIND_SV:
    default:
        return 0;
    }
}

static duckvep_status_t validate_variant_batch(
    const duckvep_model_t         *model,
    const duckvep_variant_batch_t *variants,
    duckvep_error_t               *error) {

    size_t i;
    int needs_alleles = 0;
    int allele_columns_seen;
    int validate_alleles;

    if (variants->count > (size_t)UINT32_MAX) {
        return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_ANN_VARIANT_COUNT,
                    "variant batch count exceeds uint32 variant_idx space");
    }
    for (i = 0u; i < variants->count; i++) {
        uint8_t kind = variants->variant_kind[i];
        uint8_t sv_type = variants->sv_type != NULL
                            ? variants->sv_type[i]
                            : (kind == (uint8_t)DUCKVEP_KIND_SV
                                ? (uint8_t)DUCKVEP_SV_UNKNOWN
                                : (uint8_t)DUCKVEP_SV_NONE);
        uint8_t copy_change = variants->copy_change != NULL
                                ? variants->copy_change[i]
                                : (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;

        if (variants->pos1[i] == 0u || variants->end1[i] < variants->pos1[i] ||
            kind > (uint8_t)DUCKVEP_KIND_SV) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_COORD_RANGE,
                        "variant has an invalid interval or kind");
        }
        if ((kind != (uint8_t)DUCKVEP_KIND_SV &&
             (sv_type != (uint8_t)DUCKVEP_SV_NONE ||
              copy_change != (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN)) ||
            (kind == (uint8_t)DUCKVEP_KIND_SV &&
             !duckvep_sv_metadata_valid((duckvep_sv_type_t)sv_type,
                                        (duckvep_copy_change_t)copy_change))) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_SV_TYPE,
                        "invalid structural subtype or copy-change value");
        }
        if (kind == (uint8_t)DUCKVEP_KIND_SV &&
            sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
            return fail(error, DUCKVEP_ERR_UNSUPPORTED, DVW_ANN_SV_TYPE,
                        "breakends require paired two-locus coordinates");
        }
        if (i > 0u &&
            (variants->chrom_id[i] < variants->chrom_id[i - 1u] ||
             (variants->chrom_id[i] == variants->chrom_id[i - 1u] &&
              variants->pos1[i] < variants->pos1[i - 1u]))) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_UNSORTED,
                        "variants not sorted ascending by (chrom_id, pos1)");
        }
        if (kind != (uint8_t)DUCKVEP_KIND_SV &&
            (model->has_seq || kind != (uint8_t)DUCKVEP_KIND_SNV)) {
            needs_alleles = 1;
        }
    }

    allele_columns_seen = variants->ref_offset != NULL || variants->ref_length != NULL ||
                          variants->alt_offset != NULL || variants->alt_length != NULL ||
                          variants->allele_bytes != NULL;
    validate_alleles = needs_alleles || allele_columns_seen;

    /* Structural events never enter sequence-delta classification. SNVs also have point
     * topology when the model has no sequence and no allele storage is supplied.
     * Every other small variant needs REF/ALT storage because topology uses the
     * allele-trimmed differing region. If allele storage is supplied for no-seq
     * SNVs, validate it too so kind/allele mismatches cannot silently change
     * topology below the boundary. */
    if (validate_alleles) {
        if (variants->ref_offset == NULL || variants->ref_length == NULL ||
            variants->alt_offset == NULL || variants->alt_length == NULL ||
            variants->allele_bytes == NULL) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_ALLELE_RANGE,
                        "small-variant topology requires REF/ALT storage");
        }
        for (i = 0u; i < variants->count; i++) {
            uint64_t roff;
            uint64_t rlen;
            uint64_t aoff;
            uint64_t alen;
            uint64_t pool;
            if (variants->variant_kind[i] == (uint8_t)DUCKVEP_KIND_SV) continue;
            roff = variants->ref_offset[i];
            rlen = variants->ref_length[i];
            aoff = variants->alt_offset[i];
            alen = variants->alt_length[i];
            pool = (uint64_t)variants->allele_bytes_len;
            if (roff > pool || rlen > pool - roff ||
                aoff > pool || alen > pool - aoff) {
                return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_ANN_ALLELE_RANGE,
                            "REF/ALT slice outside allele_bytes_len");
            }
            if (rlen > UINT16_MAX || alen > UINT16_MAX ||
                !allele_shape_matches_kind(variants->variant_kind[i],
                                           variants->pos1[i],
                                           variants->end1[i],
                                           variants->allele_bytes + (size_t)roff,
                                           (uint16_t)rlen,
                                           variants->allele_bytes + (size_t)aoff,
                                           (uint16_t)alen)) {
                return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_ALLELE_RANGE,
                            "variant kind inconsistent with REF/ALT span");
            }
        }
    }
    return DUCKVEP_OK;
}

struct duckvep_annotate_cursor {
    const duckvep_model_t         *model;
    duckvep_variant_batch_t        variants_view;
    const duckvep_variant_batch_t *variants;
    const duckvep_options_t       *options;
    duckvep_workspace_t           *workspace;
    duckvep_sweep_cursor_t         sweep;
    uint32_t                       variant_idx;
    const uint32_t                *tx_indices;
    size_t                         tx_count;
    size_t                         tx_pos;
    int                            have_slice;
    int                            done;
};

static duckvep_status_t validate_common_annotate_args(
    const duckvep_model_t         *model,
    const duckvep_variant_batch_t *variants,
    const duckvep_options_t       *options,
    duckvep_workspace_t           *workspace,
    duckvep_error_t               *error) {

    duckvep_status_t validation_status;

    if (model == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_MODEL, "model is NULL");
    }
    if (variants == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_VARIANTS,
                    "variants batch is NULL");
    }
    if (variants->count > 0u &&
        (variants->chrom_id == NULL || variants->pos1 == NULL ||
         variants->end1 == NULL || variants->variant_kind == NULL)) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_COLUMNS,
                    "variant batch has count>0 but null required columns");
    }
    validation_status = validate_variant_batch(model, variants, error);
    if (validation_status != DUCKVEP_OK) return validation_status;
    if (options == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_OPTIONS, "options is NULL");
    }
    if (workspace == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_WORKSPACE,
                    "workspace is NULL");
    }
    if (workspace->model != model) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_WS_MODEL,
                    "workspace belongs to a different model");
    }
    return DUCKVEP_OK;
}

static duckvep_status_t validate_result_builder_for_append(
    const duckvep_result_builder_t *results,
    duckvep_error_t                *error) {

    if (results == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_RESULTS, "results is NULL");
    }
    if (results->capacity > 0u && results->rows == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_RESULT_STORAGE,
                    "result builder has capacity>0 but rows is NULL");
    }
    if (results->count > results->capacity) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_RESULT_COUNT,
                    "result builder count exceeds capacity");
    }
    return DUCKVEP_OK;
}

duckvep_status_t duckvep_annotate_cursor_open(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    duckvep_annotate_cursor_t     **out_cursor,
    duckvep_error_t                *error) {

    duckvep_annotate_cursor_t *cursor;
    duckvep_status_t st;

    if (out_cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_CURSOR_NULL_OUT,
                    "out_cursor is NULL");
    }
    *out_cursor = NULL;
    st = validate_common_annotate_args(model, variants, options, workspace, error);
    if (st != DUCKVEP_OK) return st;

    cursor = (duckvep_annotate_cursor_t *)calloc(1u, sizeof *cursor);
    if (cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_CURSOR_OOM,
                    "annotate cursor allocation failed");
    }
    cursor->model = model;
    cursor->variants_view = *variants;
    cursor->variants = &cursor->variants_view;
    cursor->options = options;
    cursor->workspace = workspace;
    workspace_point_run_begin(workspace, variants);
    duckvep_sweep_cursor_init(&cursor->sweep, cursor->variants, &model->transcripts,
                              options->halo, workspace->active, workspace->active_cap,
                              workspace->candidates, workspace->active_cap);
    if (cursor->sweep.status != DUCKVEP_OK) {
        st = fail(error, cursor->sweep.status, DVW_ANN_SWEEP,
                  "candidate sweep initialization failed");
        free(cursor);
        return st;
    }
    *out_cursor = cursor;
    return DUCKVEP_OK;
}

int duckvep_annotate_cursor_done(const duckvep_annotate_cursor_t *cursor) {
    return cursor == NULL || cursor->done;
}

void duckvep_annotate_cursor_close(duckvep_annotate_cursor_t *cursor) {
    free(cursor);
}

duckvep_status_t duckvep_annotate_cursor_seed(
    duckvep_annotate_cursor_t *cursor,
    const uint32_t            *transcript_indices,
    size_t                     transcript_count,
    duckvep_error_t           *error) {

    if (cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_CURSOR_NULL,
                    "annotate cursor is NULL");
    }
    if (cursor->have_slice || cursor->done || cursor->sweep.vi != 0u) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_SWEEP,
                    "annotate cursor can only be seeded before its first fill");
    }
    if (!duckvep_sweep_cursor_seed(&cursor->sweep, transcript_indices,
                                   transcript_count)) {
        return fail(error, cursor->sweep.status, DVW_ANN_SWEEP,
                    "candidate sweep seed is invalid");
    }
    return DUCKVEP_OK;
}

duckvep_status_t duckvep_annotate_cursor_fill(
    duckvep_annotate_cursor_t      *cursor,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error) {

    struct annotate_ctx ctx;
    duckvep_status_t st;

    if (cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_CURSOR_NULL,
                    "annotate cursor is NULL");
    }
    st = validate_result_builder_for_append(results, error);
    if (st != DUCKVEP_OK) return st;
    if (cursor->done) return DUCKVEP_OK;
    if (results->count >= results->capacity) {
        return fail(error, DUCKVEP_ERR_RESULT_FULL, DVW_ANN_RESULT_FULL,
                    "result builder capacity exhausted; cursor paused");
    }

    ctx.model = cursor->model;
    ctx.variants = cursor->variants;
    ctx.options = cursor->options;
    ctx.workspace = cursor->workspace;
    ctx.delta_scratch = duckvep_workspace_delta_scratch(cursor->workspace);
    ctx.delta_route_stats = cursor->workspace->delta_route_stats_enabled
                              ? &cursor->workspace->delta_route_stats
                              : NULL;
    ctx.results = results;
    ctx.status = DUCKVEP_OK;

    for (;;) {
        if (!cursor->have_slice || cursor->tx_pos >= cursor->tx_count) {
            int rc = duckvep_sweep_cursor_next(&cursor->sweep, &cursor->variant_idx,
                                               &cursor->tx_indices, &cursor->tx_count);
            if (rc < 0) {
                return fail(error, cursor->sweep.status, DVW_ANN_SWEEP,
                            "candidate sweep failed");
            }
            if (rc == 0) {
                cursor->done = 1;
                cursor->have_slice = 0;
                return DUCKVEP_OK;
            }
            cursor->have_slice = 1;
            cursor->tx_pos = 0u;
        }

        while (cursor->tx_pos < cursor->tx_count) {
            if (!annotate_pair(cursor->variant_idx, cursor->tx_indices[cursor->tx_pos], &ctx)) {
                if (ctx.status == DUCKVEP_ERR_RESULT_FULL) {
                    return fail(error, DUCKVEP_ERR_RESULT_FULL, DVW_ANN_RESULT_FULL,
                                "result builder capacity exhausted; cursor paused");
                }
                return fail(error, ctx.status, DVW_ANN_SWEEP, "annotate cursor failed");
            }
            cursor->tx_pos++;
            if (results->count >= results->capacity) {
                return fail(error, DUCKVEP_ERR_RESULT_FULL, DVW_ANN_RESULT_FULL,
                            "result builder capacity exhausted; cursor paused");
            }
        }
        cursor->have_slice = 0;
    }
}

duckvep_status_t duckvep_annotate_tile(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error) {

    struct annotate_ctx ctx;
    duckvep_sweep_cursor_t sweep;
    duckvep_status_t sweep_status = DUCKVEP_OK;
    duckvep_status_t validation_status;
    int sweep_rc;
    int stopped = 0;
    uint32_t ctx_variant_idx;
    const uint32_t *ctx_tx_indices;
    size_t ctx_tx_count;

    /* Boundary contract checks first — fail loud, never on a partial result. */
    if (model == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_MODEL, "model is NULL");
    }
    if (variants == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_VARIANTS,
                    "variants batch is NULL");
    }
    if (variants->count > 0u &&
        (variants->chrom_id == NULL || variants->pos1 == NULL ||
         variants->end1 == NULL || variants->variant_kind == NULL)) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_COLUMNS,
                    "variant batch has count>0 but null required columns");
    }
    validation_status = validate_variant_batch(model, variants, error);
    if (validation_status != DUCKVEP_OK) return validation_status;
    if (options == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_OPTIONS, "options is NULL");
    }
    if (workspace == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_WORKSPACE,
                    "workspace is NULL");
    }
    if (workspace->model != model) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_WS_MODEL,
                    "workspace belongs to a different model");
    }
    if (results == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_NULL_RESULTS, "results is NULL");
    }
    if (results->capacity > 0u && results->rows == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_RESULT_STORAGE,
                    "result builder has capacity>0 but rows is NULL");
    }
    if (results->count > results->capacity) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_RESULT_COUNT,
                    "result builder count exceeds capacity");
    }

    ctx.model = model;
    ctx.variants = variants;
    ctx.options = options;
    ctx.workspace = workspace;
    ctx.delta_scratch = duckvep_workspace_delta_scratch(workspace);
    ctx.delta_route_stats = workspace->delta_route_stats_enabled
                              ? &workspace->delta_route_stats
                              : NULL;
    ctx.results = results;
    ctx.status = DUCKVEP_OK;
    workspace_point_run_begin(workspace, variants);

    /* One cursor step per variant, then a direct per-candidate loop in this TU.
     * annotate_pair is therefore inlineable and a full result buffer stops the
     * sweep immediately instead of draining the remaining candidate relation.
     * This entry point is tile-local; the cursor API preserves state only within the
     * supplied batch. */
    duckvep_sweep_cursor_init(&sweep, variants, &model->transcripts, options->halo,
                              workspace->active, workspace->active_cap,
                              workspace->candidates, workspace->active_cap);
    while ((sweep_rc = duckvep_sweep_cursor_next(&sweep, &ctx_variant_idx,
                                                  &ctx_tx_indices, &ctx_tx_count)) > 0) {
        size_t j;
        for (j = 0u; j < ctx_tx_count; j++) {
            if (!annotate_pair(ctx_variant_idx, ctx_tx_indices[j], &ctx)) {
                stopped = 1;
                break;
            }
        }
        if (stopped) break;
    }
    sweep_status = sweep.status;
    (void)sweep_rc;

    if (ctx.status == DUCKVEP_ERR_RESULT_FULL) {
        return fail(error, DUCKVEP_ERR_RESULT_FULL, DVW_ANN_RESULT_FULL,
                    "result builder capacity exhausted");
    }
    if (sweep_status != DUCKVEP_OK) {
        /* The only sweep failure with a non-NULL active buffer is active-set
         * overflow, which cannot happen here (active_cap bounds the largest chromosome slice). */
        return fail(error, sweep_status, DVW_ANN_SWEEP, "candidate sweep failed");
    }
    return DUCKVEP_OK;
}
