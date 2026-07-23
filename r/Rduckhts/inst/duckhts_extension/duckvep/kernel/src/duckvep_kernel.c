/*
 * duckvep_kernel.c — engine entry points: lifecycle constructors + the fused
 * annotate_tile path.
 *
 * This translation unit (and every other under kernel/) links against NO
 * DuckDB, htslib, Arrow, or Parquet symbol. It sees only borrowed typed arrays
 * (duckvep_kernel.h) and produces compact result rows. annotate_tile fuses the
 * already-oracle-tested predicate kernels: two sorted sweeps generate transcript
 * and core-regulation candidates from independent SoAs. The structural classifier
 * places transcript pairs; sequence-backed coding, structural/CNV, regulatory,
 * and motif fact producers feed the generated consequence program. HGVS and
 * general edit-set semantics remain late, separate layers.
 */
#include "duckvep_kernel.h"
#include "duckvep_model_internal.h"

#include "duckvep_classify.h"
#include "duckvep_annotation_internal.h"
#include "duckvep_compat.h"
#include "duckvep_codon.h"
#include "duckvep_delta.h"
#include "duckvep_effect.h"
#include "duckvep_event.h"
#include "duckvep_projection.h"
#include "duckvep_so.h"
#include "duckvep_sweep.h"
#include "duckvep_workspace_internal.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DUCKVEP_STRINGIFY_INNER(value) #value
#define DUCKVEP_STRINGIFY(value) DUCKVEP_STRINGIFY_INNER(value)

#if defined(__GNUC__) || defined(__clang__)
#define DUCKVEP_NOINLINE __attribute__((noinline))
#define DUCKVEP_HOT_ALIGN __attribute__((aligned(64)))
#elif defined(_MSC_VER)
#define DUCKVEP_NOINLINE __declspec(noinline)
#define DUCKVEP_HOT_ALIGN
#else
#define DUCKVEP_NOINLINE
#define DUCKVEP_HOT_ALIGN
#endif

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
    duckvep_interval_feature_model_t interval_features;
    duckvep_exon_model_t       exons;
    duckvep_sequence_pool_t    seq;
    uint32_t                  *cds_cdna_start1;
    uint32_t                  *cds_cdna_end1;
    uint32_t                  *cds_start_exon_index;
    uint8_t                   *cds_phase_offset;
    uint32_t                  *first_stop_position1;
    uint8_t                   *point_ordered;
    uint8_t                   *has_frameshift_intron;
    size_t                     max_transcripts_per_chrom;
    size_t                     max_interval_features_per_chrom;
    size_t                     max_cds_len;
    int                        has_seq;
};

static DUCKVEP_NOINLINE duckvep_status_t model_prepare_derived_layout(
    struct duckvep_model              *model,
    const duckvep_transcript_model_t  *transcripts,
    const duckvep_exon_model_t        *exons,
    duckvep_error_t                   *error);

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

static int transcript_has_frameshift_intron(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx) {

    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    size_t cnt = (size_t)transcripts->exon_count[tx_idx];
    size_t e;

    for (e = 1u; e < cnt; e++) {
        size_t previous = off + e - 1u;
        size_t current = off + e;
        uint32_t gap_start;
        uint32_t gap_end;

        if (transcripts->strand[tx_idx] >= 0) {
            gap_start = exons->end1[previous] + 1u;
            gap_end = exons->start1[current] - 1u;
        } else {
            gap_start = exons->end1[current] + 1u;
            gap_end = exons->start1[previous] - 1u;
        }
        if (gap_end >= gap_start &&
            gap_end - gap_start <=
                DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN)
            return 1;
    }
    return 0;
}

static int event_overlaps_vep_stretched_exon(
    const struct duckvep_model *model,
    size_t                      tx_idx,
    const duckvep_event_t      *event) {

    const duckvep_transcript_model_t *transcripts = &model->transcripts;
    const duckvep_exon_model_t *exons = &model->exons;
    size_t off;
    size_t cnt;
    size_t e;
    uint32_t event_lo;
    uint32_t event_hi;

    if (model->has_frameshift_intron == NULL ||
        model->has_frameshift_intron[tx_idx] == 0u) {
        return 0;
    }

    event_lo = event->feature_start1 < event->feature_end1
        ? event->feature_start1 : event->feature_end1;
    event_hi = event->feature_start1 > event->feature_end1
        ? event->feature_start1 : event->feature_end1;
    off = (size_t)transcripts->exon_offset[tx_idx];
    cnt = (size_t)transcripts->exon_count[tx_idx];
    for (e = 0u; e < cnt; e++) {
        uint32_t exon_start = exons->start1[off + e];
        uint32_t exon_end = exons->end1[off + e];
        uint32_t stretched_start =
            exon_start > DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN
            ? exon_start - DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN : 0u;
        uint32_t stretched_end =
            exon_end > UINT32_MAX - DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN
            ? UINT32_MAX
            : exon_end + DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN;

        if (event_hi >= stretched_start && event_lo <= stretched_end)
            return 1;
    }
    return 0;
}

static int model_exon_for_genomic(
    const duckvep_exon_model_t *exons,
    size_t                      offset,
    size_t                      count,
    uint32_t                    genomic_pos,
    size_t                     *exon_idx) {

    size_t i;

    for (i = offset; i < offset + count; i++) {
        if (genomic_pos >= exons->start1[i] && genomic_pos <= exons->end1[i]) {
            if (exon_idx != NULL) *exon_idx = i;
            return 1;
        }
    }
    return 0;
}

static int model_genomic_to_cdna(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos,
    uint32_t                         *cdna_pos,
    size_t                           *exon_idx) {

    size_t offset = (size_t)transcripts->exon_offset[tx_idx];
    size_t count = (size_t)transcripts->exon_count[tx_idx];
    size_t i;

    if (cdna_pos == NULL || exons->cdna_start1 == NULL ||
        exons->cdna_end1 == NULL) return 0;
    for (i = offset; i < offset + count; i++) {
        uint32_t position;

        if (genomic_pos < exons->start1[i] || genomic_pos > exons->end1[i]) continue;
        position = transcripts->strand[tx_idx] > 0
            ? exons->cdna_start1[i] + genomic_pos - exons->start1[i]
            : exons->cdna_start1[i] + exons->end1[i] - genomic_pos;
        if (position < exons->cdna_start1[i] || position > exons->cdna_end1[i]) {
            return 0;
        }
        *cdna_pos = position;
        if (exon_idx != NULL) *exon_idx = i;
        return 1;
    }
    return 0;
}

static int model_sequence_base_valid(uint8_t base) {
    uint8_t upper = (uint8_t)(base & UINT8_C(0xdf));
    return upper == (uint8_t)'A' || upper == (uint8_t)'C' ||
           upper == (uint8_t)'G' || upper == (uint8_t)'T' ||
           upper == (uint8_t)'N';
}

static int model_peptide_edit_alt_valid(uint8_t amino_acid) {
    return amino_acid == (uint8_t)'*' ||
           (amino_acid >= (uint8_t)'A' && amino_acid <= (uint8_t)'Z');
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
    DVW_OPTIONS_PROFILE     = 22u,
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
    DVW_WS_MODEL            = 65u,
    DVW_MODEL_TX_LAYOUT      = 66u,
    DVW_MODEL_EXON_LAYOUT    = 67u,
    DVW_MODEL_CDNA_LAYOUT    = 68u,
    DVW_MODEL_PHASE          = 69u,
    DVW_MODEL_CDS_PROJECTION = 70u,
    DVW_MODEL_SEQ_CONTRACT   = 71u,
    DVW_MODEL_MIRNA_LAYOUT   = 72u,
    DVW_MODEL_PEPTIDE_EDIT_LAYOUT = 73u,
    DVW_ANN_PAIR_COLUMNS     = 74u,
    DVW_ANN_PAIR_ORDER       = 75u,
    DVW_MODEL_INTERVAL_FEATURE_LAYOUT = 76u,
    DVW_ANN_FEATURE_PAIR_COLUMNS = 77u,
    DVW_ANN_FEATURE_PAIR_ORDER = 78u
};

duckvep_status_t duckvep_model_open_with_interval_features(
    const duckvep_transcript_model_t       *transcripts,
    const duckvep_exon_model_t             *exons,
    const duckvep_sequence_pool_t          *seq,
    const duckvep_interval_feature_model_t *interval_features,
    duckvep_model_t                       **out_model,
    duckvep_error_t                        *error) {

    struct duckvep_model *m;
    duckvep_status_t status;
    size_t t;
    size_t chrom_run = 0u;
    size_t max_chrom_run = 0u;
    size_t feature_chrom_run = 0u;
    size_t max_feature_chrom_run = 0u;
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
    if ((exons->cdna_start1 == NULL) != (exons->cdna_end1 == NULL) ||
        (exons->phase == NULL) != (exons->end_phase == NULL)) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_VIEW,
                    "paired exon cDNA or phase columns are incomplete");
    }
    if (interval_features != NULL) {
        size_t feature;

        if (interval_features->feature_count > (size_t)UINT32_MAX) {
            return fail(error, DUCKVEP_ERR_OUT_OF_RANGE,
                        DVW_MODEL_INTERVAL_FEATURE_LAYOUT,
                        "interval feature count exceeds uint32 index space");
        }
        if (interval_features->feature_count != 0u &&
            (interval_features->chrom_id == NULL ||
             interval_features->start1 == NULL ||
             interval_features->end1 == NULL ||
             interval_features->kind == NULL)) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG,
                        DVW_MODEL_INTERVAL_FEATURE_LAYOUT,
                        "interval feature view has count>0 but null columns");
        }
        for (feature = 0u; feature < interval_features->feature_count;
             feature++) {
            uint8_t kind = interval_features->kind[feature];

            if (interval_features->start1[feature] == 0u ||
                interval_features->start1[feature] >
                    interval_features->end1[feature] ||
                (kind != (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION &&
                 kind != (uint8_t)DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE)) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                            DVW_MODEL_INTERVAL_FEATURE_LAYOUT,
                            "interval feature has invalid coordinates or kind");
            }
            if (feature != 0u) {
                uint16_t previous_chrom =
                    interval_features->chrom_id[feature - 1u];
                uint16_t current_chrom = interval_features->chrom_id[feature];

                if (current_chrom < previous_chrom ||
                    (current_chrom == previous_chrom &&
                     interval_features->start1[feature] <
                         interval_features->start1[feature - 1u])) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_INTERVAL_FEATURE_LAYOUT,
                                "interval features are not sorted by contig and start");
                }
                if (current_chrom == previous_chrom) {
                    feature_chrom_run++;
                } else {
                    if (feature_chrom_run > max_feature_chrom_run)
                        max_feature_chrom_run = feature_chrom_run;
                    feature_chrom_run = 1u;
                }
            } else {
                feature_chrom_run = 1u;
            }
        }
        if (feature_chrom_run > max_feature_chrom_run)
            max_feature_chrom_run = feature_chrom_run;
    }
    {
        int any_mature_column = transcripts->mature_mirna_offset != NULL ||
                                transcripts->mature_mirna_start1 != NULL ||
                                transcripts->mature_mirna_end1 != NULL ||
                                transcripts->mature_mirna_count != 0u;

        if ((transcripts->mature_mirna_start1 == NULL) !=
            (transcripts->mature_mirna_end1 == NULL) ||
            (any_mature_column && transcripts->mature_mirna_offset == NULL) ||
            (transcripts->mature_mirna_count != 0u &&
             transcripts->mature_mirna_start1 == NULL) ||
            transcripts->mature_mirna_count > (size_t)UINT32_MAX) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG,
                        DVW_MODEL_MIRNA_LAYOUT,
                        "mature-miRNA side relation has incomplete columns");
        }
        if (transcripts->mature_mirna_offset != NULL &&
            (transcripts->mature_mirna_offset[0] != 0u ||
             transcripts->mature_mirna_offset[
                 transcripts->transcript_count] !=
                 (uint32_t)transcripts->mature_mirna_count)) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                        DVW_MODEL_MIRNA_LAYOUT,
                        "mature-miRNA row offsets do not cover the side relation");
        }
    }

    /* Validate every transcript once: span ordering, exon slice in range, cds
     * within span, and the (chrom_id, start1) sort order the sweep relies on. */
    for (t = 0u; t < transcripts->transcript_count; t++) {
        size_t eoff = (size_t)transcripts->exon_offset[t];
        size_t ecnt = (size_t)transcripts->exon_count[t];
        uint32_t cds_s = transcripts->cds_start1[t];
        uint32_t cds_e = transcripts->cds_end1[t];

        if (transcripts->start1[t] == 0u ||
            transcripts->start1[t] > transcripts->end1[t] ||
            (transcripts->strand[t] != 1 && transcripts->strand[t] != -1)) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_TX_LAYOUT,
                        "transcript has an invalid span or strand");
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
        if (ecnt > 0u) {
            uint32_t genomic_min = UINT32_MAX;
            uint32_t genomic_max = 0u;
            size_t e;

            for (e = 0u; e < ecnt; e++) {
                size_t ei = eoff + e;
                uint32_t exon_len;

                if (exons->start1[ei] == 0u ||
                    exons->start1[ei] > exons->end1[ei] ||
                    exons->start1[ei] < transcripts->start1[t] ||
                    exons->end1[ei] > transcripts->end1[t]) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_EXON_LAYOUT,
                                "exon lies outside its transcript span");
                }
                if (e > 0u) {
                    size_t previous = ei - 1u;
                    if ((transcripts->strand[t] > 0 &&
                         exons->start1[ei] <= exons->end1[previous]) ||
                        (transcripts->strand[t] < 0 &&
                         exons->end1[ei] >= exons->start1[previous])) {
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_EXON_LAYOUT,
                                    "exons overlap or are not in transcript order");
                    }
                }
                if (exons->start1[ei] < genomic_min) genomic_min = exons->start1[ei];
                if (exons->end1[ei] > genomic_max) genomic_max = exons->end1[ei];
                exon_len = exons->end1[ei] - exons->start1[ei] + 1u;
                if (exons->cdna_start1 != NULL) {
                    uint32_t cdna_len;
                    if (exons->cdna_start1[ei] == 0u ||
                        exons->cdna_start1[ei] > exons->cdna_end1[ei]) {
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_CDNA_LAYOUT,
                                    "exon has an invalid cDNA span");
                    }
                    cdna_len = exons->cdna_end1[ei] -
                               exons->cdna_start1[ei] + 1u;
                    if (cdna_len != exon_len ||
                        (e == 0u && exons->cdna_start1[ei] != 1u) ||
                        (e > 0u &&
                         (exons->cdna_end1[ei - 1u] == UINT32_MAX ||
                          exons->cdna_start1[ei] !=
                          exons->cdna_end1[ei - 1u] + 1u))) {
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_CDNA_LAYOUT,
                                    "exon cDNA spans are not contiguous and length preserving");
                    }
                }
                if (exons->phase != NULL &&
                    (exons->phase[ei] < -1 || exons->phase[ei] > 2 ||
                     exons->end_phase[ei] < -1 || exons->end_phase[ei] > 2)) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_PHASE,
                                "exon phase is outside -1,0,1,2");
                }
            }
            if (genomic_min != transcripts->start1[t] ||
                genomic_max != transcripts->end1[t]) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                            DVW_MODEL_EXON_LAYOUT,
                            "transcript span is not the outer exon envelope");
            }
        }
        if (transcripts->mature_mirna_offset != NULL) {
            size_t begin = transcripts->mature_mirna_offset[t];
            size_t finish = transcripts->mature_mirna_offset[t + 1u];
            size_t segment;

            if (finish < begin || finish > transcripts->mature_mirna_count) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                            DVW_MODEL_MIRNA_LAYOUT,
                            "mature-miRNA row offsets are not monotone");
            }
            if (finish != begin &&
                (transcripts->flags[t] &
                 (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA) == 0u) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                            DVW_MODEL_MIRNA_LAYOUT,
                            "mature-miRNA segment belongs to a non-miRNA transcript");
            }
            for (segment = begin; segment < finish; segment++) {
                uint32_t start1 = transcripts->mature_mirna_start1[segment];
                uint32_t end1 = transcripts->mature_mirna_end1[segment];
                size_t start_exon;
                size_t end_exon;

                if (start1 == 0u || start1 > end1 ||
                    start1 < transcripts->start1[t] ||
                    end1 > transcripts->end1[t] ||
                    (segment != begin &&
                     start1 < transcripts->mature_mirna_start1[segment - 1u]) ||
                    !model_exon_for_genomic(exons, eoff, ecnt, start1,
                                            &start_exon) ||
                    !model_exon_for_genomic(exons, eoff, ecnt, end1,
                                            &end_exon) ||
                    start_exon != end_exon) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_MIRNA_LAYOUT,
                                "mature-miRNA segment is not an ordered exonic interval");
                }
            }
        }
        if (cds_s != 0u &&
            (!model_exon_for_genomic(exons, eoff, ecnt, cds_s, NULL) ||
             !model_exon_for_genomic(exons, eoff, ecnt, cds_e, NULL))) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                        DVW_MODEL_CDS_PROJECTION,
                        "CDS boundary does not project into an exon");
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
        int have_flank_columns;
        int any_flank_column;
        int any_peptide_edit_column;

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
        any_peptide_edit_column = seq->peptide_edit_offset != NULL ||
                                  seq->peptide_edit_position1 != NULL ||
                                  seq->peptide_edit_alt != NULL ||
                                  seq->peptide_edit_count != 0u;
        if ((seq->peptide_edit_position1 == NULL) !=
                (seq->peptide_edit_alt == NULL) ||
            (any_peptide_edit_column && seq->peptide_edit_offset == NULL) ||
            (seq->peptide_edit_count != 0u &&
             seq->peptide_edit_position1 == NULL) ||
            seq->peptide_edit_count > (size_t)UINT32_MAX) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG,
                        DVW_MODEL_PEPTIDE_EDIT_LAYOUT,
                        "peptide-edit side relation has incomplete columns");
        }
        if (seq->peptide_edit_offset != NULL &&
            (seq->peptide_edit_offset[0] != 0u ||
             seq->peptide_edit_offset[seq->transcript_count] !=
                 (uint32_t)seq->peptide_edit_count)) {
            return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                        DVW_MODEL_PEPTIDE_EDIT_LAYOUT,
                        "peptide-edit row offsets do not cover the side relation");
        }
        any_flank_column = seq->pre_cds_offset != NULL ||
                           seq->pre_cds_length != NULL ||
                           seq->post_cds_offset != NULL ||
                           seq->post_cds_length != NULL;
        have_flank_columns = seq->pre_cds_offset != NULL &&
                             seq->pre_cds_length != NULL &&
                             seq->post_cds_offset != NULL &&
                             seq->post_cds_length != NULL;
        if ((any_flank_column && !have_flank_columns) ||
            ((seq->flanks_complete != 0u || seq->flank_bytes_len != 0u ||
              seq->flank_bytes != NULL) && !have_flank_columns) ||
            (seq->flank_bytes_len != 0u && seq->flank_bytes == NULL) ||
            seq->flanks_complete > 1u) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_MODEL_NULL_VIEW,
                        "transcript flank pool has incomplete columns or storage");
        }
        for (t = 0u; t < seq->transcript_count; t++) {
            uint64_t off = seq->cds_offset[t];
            uint64_t len = (uint64_t)seq->cds_length[t];
            uint64_t pre_off = have_flank_columns
                ? seq->pre_cds_offset[t] : 0u;
            uint64_t pre_len = have_flank_columns
                ? (uint64_t)seq->pre_cds_length[t] : 0u;
            uint64_t post_off = have_flank_columns
                ? seq->post_cds_offset[t] : 0u;
            uint64_t post_len = have_flank_columns
                ? (uint64_t)seq->post_cds_length[t] : 0u;
            size_t flank_i;

            if (seq->peptide_edit_offset != NULL) {
                size_t begin = seq->peptide_edit_offset[t];
                size_t finish = seq->peptide_edit_offset[t + 1u];
                size_t edit;
                uint32_t previous_position = 0u;

                if (finish < begin || finish > seq->peptide_edit_count) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_PEPTIDE_EDIT_LAYOUT,
                                "peptide-edit row offsets are not monotone");
                }
                for (edit = begin; edit < finish; edit++) {
                    uint32_t position = seq->peptide_edit_position1[edit];

                    if (position == 0u || position <= previous_position ||
                        position > seq->cds_length[t] / 3u ||
                        !model_peptide_edit_alt_valid(
                            seq->peptide_edit_alt[edit])) {
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_PEPTIDE_EDIT_LAYOUT,
                                    "peptide edits are not ordered, unique, and in range");
                    }
                    previous_position = position;
                }
            }

            if (pre_off > (uint64_t)seq->flank_bytes_len ||
                pre_len > (uint64_t)seq->flank_bytes_len - pre_off ||
                post_off > (uint64_t)seq->flank_bytes_len ||
                post_len > (uint64_t)seq->flank_bytes_len - post_off) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                            DVW_MODEL_SEQ_RANGE,
                            "transcript flank offset/length is out of range");
            }
            for (flank_i = 0u; flank_i < (size_t)pre_len; flank_i++) {
                if (!model_sequence_base_valid(
                    seq->flank_bytes[(size_t)pre_off + flank_i])) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "pre-CDS sequence contains a non-ACGTN base");
                }
            }
            for (flank_i = 0u; flank_i < (size_t)post_len; flank_i++) {
                if (!model_sequence_base_valid(
                    seq->flank_bytes[(size_t)post_off + flank_i])) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "post-CDS sequence contains a non-ACGTN base");
                }
            }
            if (off > (uint64_t)seq->cds_bytes_len ||
                len > (uint64_t)seq->cds_bytes_len - off) {
                return fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_MODEL_SEQ_RANGE,
                            "cds offset/length out of the sequence pool");
            }
            if (len > 0u) {
                uint32_t coding_start_cdna;
                uint32_t coding_end_cdna;
                uint64_t expected_len;
                uint64_t expected_post_len;
                uint64_t expected_pre_len;
                uint32_t coding_start_genomic;
                uint32_t coding_end_genomic;
                size_t coding_start_exon;
                size_t i;
                int8_t phase;
                uint8_t phase_offset;

                if (transcripts->cds_start1[t] == 0u ||
                    exons->cdna_start1 == NULL || exons->phase == NULL) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "sequence-backed transcript lacks CDS projection columns");
                }
                coding_start_genomic = transcripts->strand[t] > 0
                    ? transcripts->cds_start1[t] : transcripts->cds_end1[t];
                coding_end_genomic = transcripts->strand[t] > 0
                    ? transcripts->cds_end1[t] : transcripts->cds_start1[t];
                if (!model_genomic_to_cdna(transcripts, exons, t,
                                           coding_start_genomic,
                                           &coding_start_cdna,
                                           &coding_start_exon) ||
                    !model_genomic_to_cdna(transcripts, exons, t,
                                           coding_end_genomic,
                                           &coding_end_cdna, NULL) ||
                    coding_end_cdna < coding_start_cdna) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_CDS_PROJECTION,
                                "CDS does not project to a contiguous cDNA interval");
                }
                phase = exons->phase[coding_start_exon];
                /* Ensembl uses -1 when translation starts inside an exon
                 * whose genomic start is UTR. Its translateable_seq path
                 * prepends bases only for positive phase, as does projection. */
                phase_offset = phase > 0 ? (uint8_t)phase : 0u;
                expected_len = (uint64_t)coding_end_cdna -
                               (uint64_t)coding_start_cdna + 1u +
                               (uint64_t)phase_offset;
                if (len != expected_len ||
                    !duckvep_codon_table_supported(
                        (duckvep_codon_table_t)seq->codon_table[t])) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "prepared CDS length or codon table is inconsistent");
                }
                expected_pre_len = (uint64_t)coding_start_cdna - 1u;
                expected_post_len =
                    (uint64_t)exons->cdna_end1[
                        (size_t)transcripts->exon_offset[t] +
                        (size_t)transcripts->exon_count[t] - 1u] -
                    (uint64_t)coding_end_cdna;
                if ((seq->flanks_complete != 0u &&
                     (pre_len != expected_pre_len ||
                      post_len != expected_post_len)) ||
                    (seq->flanks_complete == 0u &&
                     (pre_len > expected_pre_len ||
                      post_len > expected_post_len))) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "transcript flank lengths are inconsistent with cDNA projection");
                }
                for (i = 0u; i < (size_t)len; i++) {
                    if (!model_sequence_base_valid(seq->cds_bytes[(size_t)off + i])) {
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_SEQ_CONTRACT,
                                    "prepared CDS contains a non-ACGTN base");
                    }
                }
                if (seq->first_stop_position1 != NULL) {
                    uint32_t first_stop_position1 = 0u;

                    if (!duckvep_cds_first_stop_position1(
                            seq->cds_bytes + (size_t)off, (size_t)len,
                            (duckvep_codon_table_t)seq->codon_table[t],
                            &first_stop_position1) ||
                        seq->first_stop_position1[t] != first_stop_position1) {
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_SEQ_CONTRACT,
                                    "prepared CDS first-stop cache is inconsistent");
                    }
                }
            } else {
                if (pre_len != 0u || post_len != 0u) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "transcript without prepared CDS has sequence flanks");
                }
                if (transcripts->cds_start1[t] == 0u &&
                    seq->codon_table[t] != 0u &&
                    !duckvep_codon_table_supported(
                        (duckvep_codon_table_t)seq->codon_table[t])) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "non-coding transcript has an invalid codon-table value");
                }
                if (seq->first_stop_position1 != NULL &&
                    seq->first_stop_position1[t] != 0u) {
                    return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                DVW_MODEL_SEQ_CONTRACT,
                                "non-coding transcript has a CDS first-stop cache");
                }
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
    status = model_prepare_derived_layout(m, transcripts, exons, error);
    if (status != DUCKVEP_OK) {
        duckvep_model_close(m);
        return status;
    }
    m->max_transcripts_per_chrom = max_chrom_run;
    m->max_interval_features_per_chrom = max_feature_chrom_run;
    m->max_cds_len = max_cds_len;
    if (interval_features != NULL)
        m->interval_features = *interval_features;
    if (seq != NULL) {
        size_t t;

        m->seq = *seq;
        m->has_seq = 1;
        if (seq->transcript_count != 0u) {
            m->first_stop_position1 = (uint32_t *)calloc(
                seq->transcript_count, sizeof *m->first_stop_position1);
            if (m->first_stop_position1 == NULL) {
                duckvep_model_close(m);
                return fail(error, DUCKVEP_ERR_INTERNAL, DVW_MODEL_OOM,
                            "model first-stop cache alloc failed");
            }
            if (seq->first_stop_position1 != NULL) {
                memcpy(m->first_stop_position1, seq->first_stop_position1,
                       seq->transcript_count *
                           sizeof *m->first_stop_position1);
            } else {
                for (t = 0u; t < seq->transcript_count; t++) {
                    uint64_t offset = seq->cds_offset[t];
                    size_t length = (size_t)seq->cds_length[t];

                    if (length != 0u && !duckvep_cds_first_stop_position1(
                            seq->cds_bytes + (size_t)offset, length,
                            (duckvep_codon_table_t)seq->codon_table[t],
                            &m->first_stop_position1[t])) {
                        duckvep_model_close(m);
                        return fail(error, DUCKVEP_ERR_MODEL_INVALID,
                                    DVW_MODEL_SEQ_CONTRACT,
                                    "prepared CDS cannot derive first-stop cache");
                    }
                }
            }
            m->seq.first_stop_position1 = m->first_stop_position1;
        }
    }
    *out_model = m;
    return DUCKVEP_OK;
}

duckvep_status_t duckvep_model_open(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    duckvep_model_t                 **out_model,
    duckvep_error_t                  *error) {

    return duckvep_model_open_with_interval_features(
        transcripts, exons, seq, NULL, out_model, error);
}

void duckvep_model_close(duckvep_model_t *model) {
    if (model != NULL) {
        free(model->first_stop_position1);
        free(model->has_frameshift_intron);
        free(model->point_ordered);
        free(model->cds_phase_offset);
        free(model->cds_start_exon_index);
        free(model->cds_cdna_end1);
        free(model->cds_cdna_start1);
        free(model);
    }
}

DUCKVEP_INTERNAL_API const duckvep_transcript_model_t *
duckvep_model_prepared_transcripts(const duckvep_model_t *model) {
    return model != NULL ? &model->transcripts : NULL;
}

DUCKVEP_INTERNAL_API const duckvep_exon_model_t *
duckvep_model_prepared_exons(const duckvep_model_t *model) {
    return model != NULL ? &model->exons : NULL;
}

DUCKVEP_INTERNAL_API const duckvep_sequence_pool_t *
duckvep_model_prepared_sequences(const duckvep_model_t *model) {
    return model != NULL && model->has_seq ? &model->seq : NULL;
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
    duckvep_compat_profile_t compatibility_profile;
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
    if (init != NULL && !duckvep_compat_profile_valid(
            (duckvep_compat_profile_t)init->compatibility_profile)) {
        free(o);
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_OPTIONS_PROFILE,
                    "compatibility profile is invalid");
    }
    o->compatibility_profile = init != NULL
        ? (duckvep_compat_profile_t)init->compatibility_profile
        : DUCKVEP_COMPAT_VEP_116;

    o->upstream_dist = init != NULL && init->distances_are_explicit
        ? init->upstream_dist
        : (init != NULL && init->upstream_dist != 0u
            ? init->upstream_dist : DUCKVEP_DEFAULT_UPSTREAM_DIST);
    o->downstream_dist = init != NULL && init->distances_are_explicit
        ? init->downstream_dist
        : (init != NULL && init->downstream_dist != 0u
            ? init->downstream_dist : DUCKVEP_DEFAULT_DOWNSTREAM_DIST);
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
 * each sized to the largest prepared contig transcript run. Independent exon
 * cursors keep raw-position SNVs monotone even when a neighboring non-SNV's
 * trimmed differing-region start moves backward within the same DuckDB vector.
 * The delta scratch is sized once from the model's maximum CDS length plus the
 * maximum uint16_t small-variant payload, enough for current single-variant MNV,
 * shortening deletion, and lengthening insertion CodingContext builds. Broader
 * haplotype/grouped edit-set paths will widen this policy explicitly when they
 * are wired into production. */
struct duckvep_workspace {
    const duckvep_model_t    *model;
    uint32_t                *active;
    uint32_t                *candidates;
    uint32_t                *interval_feature_active;
    uint32_t                *interval_feature_candidates;
    uint16_t                *point_exon_rank;
    uint16_t                *span_exon_rank;
    size_t                   point_exon_count;
    size_t                   active_cap;
    size_t                   interval_feature_active_cap;
    uint16_t                 point_last_chrom;
    uint32_t                 point_last_pos;
    int                      point_run_active;
    int                      point_cursor_reset_pending;
    uint16_t                 span_last_chrom;
    uint32_t                 span_last_feature_start;
    int                      span_run_active;
    int                      span_cursor_reset_pending;
    duckvep_delta_scratch_t                 delta_scratch;
    duckvep_workspace_delta_route_stats_t   delta_route_stats;
    int                                     delta_route_stats_enabled;
    int                                     force_generalized_annotation;
};

static void workspace_point_cursor_reset(duckvep_workspace_t *workspace) {
    if (workspace == NULL || workspace->point_exon_rank == NULL) return;
    memset(workspace->point_exon_rank, 0xff,
           workspace->point_exon_count * sizeof *workspace->point_exon_rank);
}

static void workspace_span_cursor_reset(duckvep_workspace_t *workspace) {
    if (workspace == NULL || workspace->span_exon_rank == NULL) return;
    memset(workspace->span_exon_rank, 0xff,
           workspace->point_exon_count * sizeof *workspace->span_exon_rank);
}

DUCKVEP_INTERNAL_API void duckvep_workspace_force_generalized_annotation(
    duckvep_workspace_t *workspace,
    int                   enabled) {

    if (workspace == NULL) return;
    workspace->force_generalized_annotation = enabled != 0;
    workspace_point_cursor_reset(workspace);
    workspace_span_cursor_reset(workspace);
    workspace->point_run_active = 0;
    workspace->span_run_active = 0;
    workspace->point_cursor_reset_pending = 0;
    workspace->span_cursor_reset_pending = 0;
}

static int workspace_sorted_runs_begin(
    duckvep_workspace_t           *workspace,
    const duckvep_variant_batch_t *variants,
    const duckvep_event_t         *events,
    int                           *span_monotonic_out) {

    uint16_t first_chrom;
    uint32_t first_pos;
    uint32_t first_feature_start;
    int point_monotonic = 1;
    int span_monotonic = 1;
    size_t i;
    size_t last;

    if (span_monotonic_out != NULL) *span_monotonic_out = 1;
    if (workspace == NULL || variants == NULL || variants->count == 0u) return 1;
    first_chrom = variants->chrom_id[0];
    first_pos = events != NULL ? events[0].start1 : variants->pos1[0];
    first_feature_start = events != NULL
        ? events[0].feature_start1 : variants->pos1[0];

    /* A non-monotone normalized vector can leave an individual transcript's
     * cursor ahead of the vector's final coordinate when later variants do not
     * admit that transcript.  The next vector can still look globally
     * monotone relative to that final coordinate, so reset before inspecting
     * it rather than carrying transcript-local stale ranks across the DuckDB
     * vector edge.  Rewinds within the current vector continue to use mode 2
     * below and avoid a model-wide reset per event. */
    if (workspace->point_cursor_reset_pending) {
        workspace_point_cursor_reset(workspace);
        workspace->point_cursor_reset_pending = 0;
    }
    if (workspace->span_cursor_reset_pending) {
        workspace_span_cursor_reset(workspace);
        workspace->span_cursor_reset_pending = 0;
    }

    for (i = 1u; i < variants->count; i++) {
        uint32_t previous_pos = events != NULL
            ? events[i - 1u].start1 : variants->pos1[i - 1u];
        uint32_t current_pos = events != NULL
            ? events[i].start1 : variants->pos1[i];
        uint32_t previous_feature_start = events != NULL
            ? events[i - 1u].feature_start1 : variants->pos1[i - 1u];
        uint32_t current_feature_start = events != NULL
            ? events[i].feature_start1 : variants->pos1[i];
        if (variants->chrom_id[i] < variants->chrom_id[i - 1u] ||
            (variants->chrom_id[i] == variants->chrom_id[i - 1u] &&
             current_pos < previous_pos))
            point_monotonic = 0;
        if (variants->chrom_id[i] < variants->chrom_id[i - 1u] ||
            (variants->chrom_id[i] == variants->chrom_id[i - 1u] &&
             current_feature_start < previous_feature_start))
            span_monotonic = 0;
    }
    if (workspace->point_run_active &&
        (first_chrom < workspace->point_last_chrom ||
         (first_chrom == workspace->point_last_chrom &&
          first_pos < workspace->point_last_pos))) {
        workspace_point_cursor_reset(workspace);
    }
    if (!point_monotonic) workspace_point_cursor_reset(workspace);
    if (workspace->span_run_active &&
        (first_chrom < workspace->span_last_chrom ||
         (first_chrom == workspace->span_last_chrom &&
          first_feature_start < workspace->span_last_feature_start))) {
        workspace_span_cursor_reset(workspace);
    }
    if (!span_monotonic) workspace_span_cursor_reset(workspace);
    last = variants->count - 1u;
    workspace->point_last_chrom = variants->chrom_id[last];
    workspace->point_last_pos = events != NULL
        ? events[last].start1 : variants->pos1[last];
    workspace->point_run_active = 1;
    workspace->span_last_chrom = variants->chrom_id[last];
    workspace->span_last_feature_start = events != NULL
        ? events[last].feature_start1 : variants->pos1[last];
    workspace->span_run_active = 1;
    workspace->point_cursor_reset_pending = !point_monotonic;
    workspace->span_cursor_reset_pending = !span_monotonic;
    /* Mode 1 is the minimum-instruction monotone cursor. Mode 2 repairs only
     * the normalized-coordinate rewinds observed in this input vector. A
     * non-monotone vector also schedules one reset before the next vector so a
     * transcript skipped after the forward jump cannot retain an ahead rank.
     * Both modes retain cursor traversal. */
    if (span_monotonic_out != NULL)
        *span_monotonic_out = span_monotonic ? 1 : 2;
    return point_monotonic ? 1 : 2;
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
    size_t interval_feature_cap;
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
    interval_feature_cap = model->max_interval_features_per_chrom;
    if (interval_feature_cap == 0u) interval_feature_cap = 1u;
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
    w->interval_feature_active = (uint32_t *)calloc(
        interval_feature_cap, sizeof *w->interval_feature_active);
    w->interval_feature_candidates = (uint32_t *)calloc(
        interval_feature_cap, sizeof *w->interval_feature_candidates);
    if (model->transcripts.transcript_count > 0u) {
        w->point_exon_rank = (uint16_t *)malloc(point_bytes);
        w->span_exon_rank = (uint16_t *)malloc(point_bytes);
    }
    if (w->active == NULL || w->candidates == NULL ||
        w->interval_feature_active == NULL ||
        w->interval_feature_candidates == NULL ||
        (model->transcripts.transcript_count > 0u &&
         (w->point_exon_rank == NULL || w->span_exon_rank == NULL))) {
        free(w->active);
        free(w->candidates);
        free(w->interval_feature_active);
        free(w->interval_feature_candidates);
        free(w->point_exon_rank);
        free(w->span_exon_rank);
        free(w);
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_WS_OOM,
                    "workspace sweep scratch alloc failed");
    }
    st = workspace_delta_scratch_open(model, &w->delta_scratch, error);
    if (st != DUCKVEP_OK) {
        free(w->active);
        free(w->candidates);
        free(w->interval_feature_active);
        free(w->interval_feature_candidates);
        free(w->point_exon_rank);
        free(w->span_exon_rank);
        free(w);
        return st;
    }
    w->model = model;
    w->point_exon_count = model->transcripts.transcript_count;
    workspace_point_cursor_reset(w);
    workspace_span_cursor_reset(w);
    w->active_cap = cap;
    w->interval_feature_active_cap = interval_feature_cap;
    *out_workspace = w;
    return DUCKVEP_OK;
}

void duckvep_workspace_close(duckvep_workspace_t *workspace) {
    if (workspace != NULL) {
        free(workspace->active);
        free(workspace->candidates);
        free(workspace->interval_feature_active);
        free(workspace->interval_feature_candidates);
        free(workspace->point_exon_rank);
        free(workspace->span_exon_rank);
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
    const duckvep_event_t         *events;
    const duckvep_options_t       *options;
    duckvep_workspace_t           *workspace;
    duckvep_delta_scratch_t       *delta_scratch;
    duckvep_workspace_delta_route_stats_t *delta_route_stats;
    duckvep_result_builder_t      *results;
    duckvep_annotation_observer_fn observer;
    void                          *observer_context;
    duckvep_status_t               status; /* first non-OK halts emission */
    int                            point_sorted_safe;
    int                            span_sorted_safe;
    /* Event normalization and feature-allele slicing are properties of the
     * uploaded allele, not of each overlapping transcript. Dense loci often
     * emit dozens of transcript rows per allele, so prepare these facts once
     * and reuse them until candidate traversal advances to the next variant. */
    uint32_t                       prepared_variant_idx;
    duckvep_event_t                prepared_event;
    const uint8_t                 *prepared_feature_ref;
    const uint8_t                 *prepared_feature_alt;
    uint16_t                       prepared_feature_ref_length;
    uint16_t                       prepared_feature_alt_length;
    int                            prepared_have_feature_alleles;
};

static void annotate_ctx_prepare_variant(
    struct annotate_ctx *ctx,
    uint32_t             variant_idx) {

    if (ctx == NULL || ctx->prepared_variant_idx == variant_idx) return;
    ctx->prepared_variant_idx = variant_idx;
    ctx->prepared_event = ctx->events[variant_idx];
    ctx->prepared_feature_ref = NULL;
    ctx->prepared_feature_alt = NULL;
    ctx->prepared_feature_ref_length = 0u;
    ctx->prepared_feature_alt_length = 0u;
    ctx->prepared_have_feature_alleles = duckvep_event_feature_alleles(
        ctx->variants, variant_idx, &ctx->prepared_event,
        &ctx->prepared_feature_ref, &ctx->prepared_feature_ref_length,
        &ctx->prepared_feature_alt, &ctx->prepared_feature_alt_length);
}

static int insertion_placement_uses_right_flank(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event) {

    uint32_t boundary;
    uint32_t right;
    uint32_t cds_end;
    duckvep_coding_projection_t projection;

    if (event == NULL || !event->interbase ||
        event->anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT) {
        return 0;
    }
    boundary = event->insertion_boundary0;
    right = duckvep_event_right_flank1(event);
    if (boundary == transcripts->end1[tx_idx]) return 1;
    cds_end = transcripts->cds_end1[tx_idx];
    if (cds_end != 0u && boundary == cds_end) return 1;

    /* At an internal coding-exon entrance the VCF padding base is intronic but
     * the inserted sequence maps immediately before the right-hand CDS base.
     * The global CDS entrance is deliberately excluded: VEP's _before_coding
     * special case assigns that boundary to the UTR. */
    if (transcripts->cds_start1[tx_idx] == 0u ||
        right == transcripts->cds_start1[tx_idx] ||
        duckvep_project_coding_base(transcripts, exons, tx_idx,
                                    boundary, &projection)) {
        return 0;
    }
    return duckvep_project_coding_base(transcripts, exons, tx_idx,
                                       right, &projection);
}

static int insertion_flanks_share_exon(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          left1,
    uint32_t                          right1) {

    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    size_t cnt = (size_t)transcripts->exon_count[tx_idx];
    size_t e;

    for (e = 0u; e < cnt; e++) {
        size_t ei = off + e;
        if (left1 >= exons->start1[ei] && right1 <= exons->end1[ei]) {
            return 1;
        }
    }
    return 0;
}

static int sorted_exon_index_at_rank(
    const duckvep_transcript_model_t *transcripts,
    size_t                            tx_idx,
    uint16_t                          exon_rank,
    size_t                           *exon_idx) {

    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    uint32_t cnt = (uint32_t)transcripts->exon_count[tx_idx];

    if (exon_idx == NULL || (uint32_t)exon_rank >= cnt) return 0;
    *exon_idx = transcripts->strand[tx_idx] >= 0
        ? off + (size_t)exon_rank
        : off + (size_t)(cnt - (uint32_t)exon_rank - 1u);
    return 1;
}

static int insertion_flanks_share_sorted_exon(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint16_t                          exon_rank,
    uint32_t                          left1,
    uint32_t                          right1) {

    size_t ei;

    if (!sorted_exon_index_at_rank(
            transcripts, tx_idx, exon_rank, &ei)) return 0;
    return left1 >= exons->start1[ei] && right1 <= exons->end1[ei];
}

/* VEP keeps three insertion questions distinct. Its generic exon overlap uses
 * the reversed feature interval P+1,P, so both flanks must lie in one exon.
 * UTR predicates instead use the transcript mapper and may accept either
 * exonic flank at an internal boundary. Splice predicates retain the reversed
 * interval and are computed elsewhere. A single point placement cannot express
 * all three, so correct only the region facts after the cheap point classifier. */
static void insertion_apply_region_boundaries(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    int                               sorted_span,
    uint16_t                          exon_rank,
    duckvep_region_state_t           *region) {

    uint32_t left1;
    uint32_t right1;
    uint32_t cds_start;
    uint32_t cds_end;
    int flanks_share_exon;

    if (transcripts == NULL || exons == NULL || event == NULL || region == NULL ||
        !event->interbase ||
        event->anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT) {
        return;
    }
    left1 = event->insertion_boundary0;
    right1 = duckvep_event_right_flank1(event);
    cds_start = transcripts->cds_start1[tx_idx];
    cds_end = transcripts->cds_end1[tx_idx];
    flanks_share_exon = sorted_span
        ? insertion_flanks_share_sorted_exon(
            transcripts, exons, tx_idx, exon_rank, left1, right1)
        : insertion_flanks_share_exon(
            transcripts, exons, tx_idx, left1, right1);

    /* VEP 116 clamps a reversed insertion interval to one exonic flank for
     * within_cdna, even though its generic exon predicate still requires both
     * flanks to share an exon. This combination enables feature_elongation plus
     * non_coding_transcript_variant at an internal exon entrance. */
    {
        uint32_t cdna;
        if (duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, left1, &cdna, NULL) ||
            duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, right1, &cdna, NULL)) {
            region->within_cdna = 1u;
        }
    }

    if (cds_start == 0u && region->overlaps_exon &&
        !flanks_share_exon) {
        region->overlaps_exon = 0u;
        region->region_mask &= ~(uint32_t)DUCKVEP_REGION_EXON;
    }

    /* An insertion at an outer transcript edge is upstream/downstream, not UTR.
     * Internal exon boundaries may map through either flank after VEP 116's
     * transcript-coordinate clamping. Only import UTR facts from those flanks;
     * CDS placement remains governed by the established coding projection. */
    if (cds_start != 0u && left1 >= transcripts->start1[tx_idx] &&
        right1 <= transcripts->end1[tx_idx]) {
        int before_coding;
        int after_coding;

        /* Away from an exon or CDS boundary, both flanks have the region facts
         * already carried by the chosen placement point. Only the boundary
         * case needs VEP's two independent mapper queries. */
        if (flanks_share_exon &&
            !(left1 < cds_start && right1 >= cds_start) &&
            !(left1 <= cds_end && right1 > cds_end)) {
            return;
        }
        duckvep_region_state_t left = duckvep_region_classify_span(
            transcripts, exons, tx_idx, left1, left1, 0u, 0u);
        duckvep_region_state_t right = duckvep_region_classify_span(
            transcripts, exons, tx_idx, right1, right1, 0u, 0u);

        if (left.within_cdna || right.within_cdna) region->within_cdna = 1u;

        /* VariationEffect::_before_coding and _after_coding each have an exact
         * insertion exception. The mapper may accept the coding-side flank
         * while the chosen topology point is intronic, so derive these two
         * UTR facts from VEP's reversed feature interval P+1,P rather than from
         * either flank's ordinary base classification. */
        before_coding = event->feature_start1 > event->feature_end1 &&
            event->feature_start1 - event->feature_end1 == 1u &&
            event->feature_start1 == cds_start;
        after_coding = event->feature_start1 > event->feature_end1 &&
            event->feature_start1 - event->feature_end1 == 1u &&
            event->feature_end1 == cds_end;
        if (before_coding || after_coding) {
            /* BaseTranscriptVariation::cds_coords is a mapper Gap for the
             * reversed interval straddling the outer CDS edge. VEP therefore
             * admits the UTR special case but excludes every consequence whose
             * registry entry requires coding=1. */
            region->overlaps_cds = 0u;
            region->region_mask &= ~(uint32_t)DUCKVEP_REGION_CDS;
        }
        if (region->within_cdna) {
            if ((transcripts->strand[tx_idx] >= 0 && before_coding) ||
                (transcripts->strand[tx_idx] < 0 && after_coding)) {
                region->overlaps_utr5 = 1u;
                region->region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
            }
            if ((transcripts->strand[tx_idx] >= 0 && after_coding) ||
                (transcripts->strand[tx_idx] < 0 && before_coding)) {
                region->overlaps_utr3 = 1u;
                region->region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
            }
        }
        if (left.overlaps_utr5 || right.overlaps_utr5) {
            region->overlaps_utr5 = 1u;
            region->region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
        }
        if (left.overlaps_utr3 || right.overlaps_utr3) {
            region->overlaps_utr3 = 1u;
            region->region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
        }
    }
}

static int breakend_endpoint_close_to_transcript(
    const duckvep_transcript_model_t *transcripts,
    size_t                            tx_idx,
    uint16_t                          chrom_id,
    uint32_t                          pos1) {

    uint32_t start1;
    uint32_t end1;

    if (transcripts->chrom_id[tx_idx] != chrom_id || pos1 == 0u) return 0;
    start1 = transcripts->start1[tx_idx];
    end1 = transcripts->end1[tx_idx];
    if (pos1 < start1) return start1 - pos1 <= DUCKVEP_BREAKEND_ALLELE_DISTANCE;
    if (pos1 > end1) return pos1 - end1 <= DUCKVEP_BREAKEND_ALLELE_DISTANCE;
    return 1;
}

static int breakend_endpoint_in_directional_window(
    const duckvep_transcript_model_t *transcripts,
    size_t                            tx_idx,
    uint16_t                          chrom_id,
    uint32_t                          pos1,
    const duckvep_options_t          *options) {

    uint32_t distance;
    uint32_t allowed;

    if (transcripts->chrom_id[tx_idx] != chrom_id || pos1 == 0u) return 0;
    if (pos1 < transcripts->start1[tx_idx]) {
        distance = transcripts->start1[tx_idx] - pos1;
        allowed = transcripts->strand[tx_idx] >= 0
            ? options->upstream_dist : options->downstream_dist;
        return distance <= allowed;
    }
    if (pos1 > transcripts->end1[tx_idx]) {
        distance = pos1 - transcripts->end1[tx_idx];
        allowed = transcripts->strand[tx_idx] >= 0
            ? options->downstream_dist : options->upstream_dist;
        return distance <= allowed;
    }
    return 1;
}

/* A one-base event far outside the outermost splice reach can only emit the
 * strand-relative upstream/downstream term for this transcript. Keep the SO
 * rule program authoritative, but avoid constructing the full region, splice,
 * sequence, and observer state for ordinary rich/compact output. HGVS keeps
 * the generalized path because its observer consumes the complete pair facts. */
static int annotate_far_directional_snv(
    struct annotate_ctx       *c,
    uint32_t                   variant_idx,
    uint32_t                   tx_idx,
    const duckvep_event_t     *event,
    int                       *handled) {

    const duckvep_transcript_model_t *tx = &c->model->transcripts;
    duckvep_consequence_t *row;
    uint32_t pos;
    uint32_t distance;
    uint32_t allowed;
    uint32_t region;
    uint64_t pre;
    uint64_t cmask;
    int upstream;

    *handled = 0;
    if (c->observer != NULL || c->workspace->force_generalized_annotation ||
        event == NULL || event->kind != (uint8_t)DUCKVEP_KIND_SNV ||
        event->feature_start1 != event->feature_end1 ||
        c->model->point_ordered[tx_idx] == 0u) {
        return 1;
    }
    pos = event->feature_start1;
    if (pos >= tx->start1[tx_idx] && pos <= tx->end1[tx_idx]) return 1;
    if (!duckvep_point_outside_splice_reach(
            tx, &c->model->exons, (size_t)tx_idx, pos,
            c->options->splice_region_exonic,
            c->options->splice_region_intronic)) {
        return 1;
    }

    *handled = 1;
    if (c->delta_route_stats != NULL) {
        c->delta_route_stats->far_directional_snv++;
    }
    upstream = tx->strand[tx_idx] >= 0
        ? pos < tx->start1[tx_idx]
        : pos > tx->end1[tx_idx];
    if (pos < tx->start1[tx_idx]) {
        distance = tx->start1[tx_idx] - pos;
    } else {
        distance = pos - tx->end1[tx_idx];
    }
    allowed = upstream ? c->options->upstream_dist
                       : c->options->downstream_dist;
    if (distance > allowed) return 1;

    region = upstream ? (uint32_t)DUCKVEP_REGION_UPSTREAM
                      : (uint32_t)DUCKVEP_REGION_DOWNSTREAM;
    pre = upstream ? DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM)
                   : DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM);
    pre |= tx->cds_start1[tx_idx] != 0u
        ? DUCKVEP_PRE(DUCKVEP_PRE_CODING)
        : DUCKVEP_PRE(DUCKVEP_PRE_NONCODING);
    cmask = duckvep_effect_eval(pre);
    if (cmask == 0u) cmask = DUCKVEP_SO(DUCKVEP_SO_INTERGENIC);

    if (c->results->count >= c->results->capacity) {
        c->status = DUCKVEP_ERR_RESULT_FULL;
        return 0;
    }
    row = &c->results->rows[c->results->count++];
    memset(row, 0, sizeof *row);
    row->variant_idx = variant_idx;
    row->tx_idx = tx_idx;
    row->overlap_object_kind =
        (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
    row->consequence_mask = cmask;
    row->region_mask = region;
    row->impact = (uint8_t)duckvep_so_impact(cmask);
    row->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE;
    row->nmd_prediction = (uint8_t)DUCKVEP_NMD_NOT_APPLICABLE;
    row->cdna_pos = -1;
    row->cds_pos = -1;
    row->protein_pos = -1;
    return 1;
}

/* A sorted SNV strictly inside a transcript but outside CDS and every exact
 * splice predicate has no sequence, structural, insertion/deletion, or NMD
 * prediction work left.  Finalize the already-classified topology through the
 * generated VEP rule program and emit the compact kernel row directly.  Keep
 * miRNA and frameshift-intron transcripts on the generalized path because
 * their VEP predicates consult additional event or mapper state. */
static int annotate_simple_point_snv_classified(
    struct annotate_ctx       *c,
    uint32_t                   variant_idx,
    uint32_t                   tx_idx,
    const duckvep_event_t     *event,
    const duckvep_region_state_t *region_state,
    const duckvep_splice_state_t *splice,
    int                       *handled) {

    const duckvep_transcript_model_t *tx = &c->model->transcripts;
    duckvep_consequence_t *row;
    uint64_t pre;
    uint64_t cmask;
    int coding;
    int noncoding_exon;

    *handled = 0;
    if (c->observer != NULL || event == NULL || region_state == NULL ||
        splice == NULL ||
        event->kind != (uint8_t)DUCKVEP_KIND_SNV ||
        !region_state->within_feature || region_state->overlaps_cds ||
        region_state->within_frameshift_intron || splice->any ||
        (tx->flags[tx_idx] & (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA) != 0u) {
        return 1;
    }
    *handled = 1;
    if (c->delta_route_stats != NULL) {
        c->delta_route_stats->simple_point_snv++;
    }

    /* This is the scalar projection of effect_ctx_set_topology + finalize for
     * an ordinary within-transcript SNV with no CDS, splice, miRNA, structural,
     * or length-change facts. The generated rule evaluator remains the SO
     * authority; the optimized-vs-generalized property compares every emitted
     * row over arbitrary sorted scenes. */
    coding = tx->cds_start1[tx_idx] != 0u;
    pre = coding ? DUCKVEP_PRE(DUCKVEP_PRE_CODING)
                 : DUCKVEP_PRE(DUCKVEP_PRE_NONCODING);
    if (region_state->overlaps_utr5)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_UTR5);
    if (region_state->overlaps_utr3)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_UTR3);
    if (!coding && region_state->overlaps_exon)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_EXON);
    if (region_state->overlaps_intron)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_INTRON);
    if (splice->intronic)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON);
    if ((tx->flags[tx_idx] & (uint64_t)DUCKVEP_TX_BIOTYPE_NMD) != 0u)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NMD_TRANSCRIPT);
    if (!coding) {
        noncoding_exon = region_state->overlaps_exon &&
            !region_state->complete_overlap_feature;
        pre |= noncoding_exon
            ? DUCKVEP_PRE(DUCKVEP_PRE_NONCODING_EXON)
            : DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NONCODING_GENE);
    }
    cmask = duckvep_effect_eval(pre);
    if (cmask == 0u) cmask = DUCKVEP_SO(DUCKVEP_SO_INTERGENIC);

    if (c->results->count >= c->results->capacity) {
        c->status = DUCKVEP_ERR_RESULT_FULL;
        return 0;
    }
    row = &c->results->rows[c->results->count++];
    memset(row, 0, sizeof *row);
    row->variant_idx = variant_idx;
    row->tx_idx = tx_idx;
    row->overlap_object_kind =
        (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
    row->consequence_mask = cmask;
    row->region_mask = region_state->region_mask;
    row->impact = (uint8_t)duckvep_so_impact(cmask);
    row->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE;
    row->nmd_prediction = (uint8_t)DUCKVEP_NMD_NOT_APPLICABLE;
    row->cdna_pos = -1;
    row->cds_pos = -1;
    row->protein_pos = -1;
    return 1;
}

/* The VEP-shaped per-candidate decision: fill the cheap facts (effect ctx),
 * escalate to the sequence delta ONLY for the CDS bucket (lazy), then evaluate the
 * static rule table. No biological special-casing lives here — every SO decision is
 * a rule in duckvep_effect.c. */
/* This is the per-candidate hot loop. Keep its entry on one instruction-cache
 * line so changes to one-time model preparation do not perturb common SNVs. */
static DUCKVEP_HOT_ALIGN int annotate_pair(
    uint32_t variant_idx, uint32_t tx_idx, void *vctx) {
    struct annotate_ctx *c = (struct annotate_ctx *)vctx;
    const duckvep_transcript_model_t *tx;
    uint32_t pos;
    uint32_t topology_start1;
    uint32_t topology_end1;
    const uint8_t *feature_ref = NULL;
    const uint8_t *feature_alt = NULL;
    uint16_t feature_ref_length = 0u;
    uint16_t feature_alt_length = 0u;
    int have_feature_alleles;
    int fwd;
    uint32_t dist;
    duckvep_variant_kind_t kind;
    const duckvep_event_t *event;
    duckvep_effect_ctx_t ectx;
    duckvep_sequence_delta_t delta;
    duckvep_coding_context_t coding_context;
    duckvep_variant_coding_context_status_t coding_context_status;
    duckvep_transcript_edit_t transcript_edit;
    duckvep_transcript_edit_status_t transcript_edit_status;
    duckvep_nmd_result_t nmd;
    uint64_t cmask;
    uint32_t consequence_flags;
    duckvep_consequence_t *row;
    uint16_t projection_exon_rank = UINT16_MAX;
    uint32_t projection_exon_hint = UINT32_MAX;
    int cds_delta_attempted = 0;
    int feature_endpoints_map_to_cdna;
    int feature_has_internal_cdna_gap;
    int feature_mapping_blocks_peptide;
    int has_partial_terminal_codon;
    int breakend;
    int breakend_local_close = 0;
    int breakend_local_in_directional_window = 0;
    int breakend_mate_close = 0;
    int breakend_mate_admits = 0;
    int breakend_local_defaults_intergenic = 0;

    if (c->status != DUCKVEP_OK) return 0;

    tx = &c->model->transcripts;
    annotate_ctx_prepare_variant(c, variant_idx);
    event = &c->prepared_event;
    breakend = event->kind == (uint8_t)DUCKVEP_KIND_SV &&
        event->sv_type == (uint8_t)DUCKVEP_SV_BREAKEND;
    if (breakend) {
        int local_admits;

        breakend_local_close = breakend_endpoint_close_to_transcript(
            tx, (size_t)tx_idx, event->chrom_id, event->feature_start1);
        breakend_mate_close = event->has_mate &&
            breakend_endpoint_close_to_transcript(
                tx, (size_t)tx_idx, event->mate_chrom_id, event->mate_pos1);
        breakend_local_in_directional_window =
            breakend_endpoint_in_directional_window(
                tx, (size_t)tx_idx, event->chrom_id, event->feature_start1,
                c->options);
        local_admits = breakend_local_close &&
            breakend_local_in_directional_window;
        breakend_mate_admits = breakend_mate_close &&
            breakend_endpoint_in_directional_window(
                tx, (size_t)tx_idx, event->mate_chrom_id, event->mate_pos1,
                c->options);
        if (!local_admits && !breakend_mate_admits) return 1;
        breakend_local_defaults_intergenic = breakend_local_close &&
            !breakend_local_in_directional_window &&
            breakend_mate_admits;
    }
    pos = event->start1;
    topology_start1 = event->feature_start1;
    topology_end1 = event->feature_end1;
    have_feature_alleles = c->prepared_have_feature_alleles;
    feature_ref = c->prepared_feature_ref;
    feature_alt = c->prepared_feature_alt;
    feature_ref_length = c->prepared_feature_ref_length;
    feature_alt_length = c->prepared_feature_alt_length;
    if (event->interbase) {
        /* Region placement for an insertion is transcript-relative, while its
         * splice geometry remains VEP's reversed feature interval P+1,P. */
        topology_start1 = event->start1;
        topology_end1 = event->start1;
        if (insertion_placement_uses_right_flank(
                tx, &c->model->exons, (size_t)tx_idx, event)) {
            topology_start1 = duckvep_event_right_flank1(event);
            topology_end1 = topology_start1;
        }
    }
    fwd = tx->strand[tx_idx] >= 0;
    kind = (duckvep_variant_kind_t)event->kind;

    {
        int handled;
        int ok = annotate_far_directional_snv(
            c, variant_idx, tx_idx, event, &handled);

        if (!ok || handled) return ok;
    }

    if (breakend && !breakend_local_in_directional_window) {
        duckvep_region_state_t region;
        duckvep_splice_state_t splice;

        /* A transcript discovered only through the mate still receives a
         * StructuralVariationOverlapAllele. Ordinary predicates keep the
         * local feature, but it has no configured topology here; only the
         * mate-aware truncation predicate can add a non-default term. A local
         * point outside the fixed 5000-base allele-admission range is not
         * sufficient to take this branch when it remains inside the caller's
         * wider directional window: predicates on the mate allele still read
         * that local feature and can emit upstream/downstream. */
        memset(&region, 0, sizeof region);
        memset(&splice, 0, sizeof splice);
        duckvep_effect_ctx_fill_classified(
            tx, variant_idx, (size_t)tx_idx,
            event->feature_start1, event->feature_end1,
            &region, &splice, &ectx);
    } else if (kind == DUCKVEP_KIND_SNV && c->point_sorted_safe &&
        c->model->point_ordered[tx_idx] &&
        event->feature_start1 == event->feature_end1 &&
        (!have_feature_alleles ||
         (feature_ref_length == 1u && feature_alt_length == 1u))) {
        duckvep_region_state_t region;
        duckvep_splice_state_t splice;

        duckvep_classify_point_sorted(
            tx, &c->model->exons, (size_t)tx_idx,
            event->feature_start1,
            c->options->splice_region_exonic,
            c->options->splice_region_intronic,
            (uint8_t)(c->point_sorted_safe == 2),
            &c->workspace->point_exon_rank[tx_idx], &region, &splice);
        projection_exon_rank = c->workspace->point_exon_rank[tx_idx];
        {
            int handled;
            int ok = annotate_simple_point_snv_classified(
                c, variant_idx, tx_idx, event, &region, &splice, &handled);

            if (!ok || handled) return ok;
        }
        duckvep_effect_ctx_fill_classified(
            tx, variant_idx, (size_t)tx_idx,
            event->feature_start1, event->feature_end1,
            &region, &splice, &ectx);
    } else if (kind != DUCKVEP_KIND_SV && have_feature_alleles) {
        int sorted_span = c->span_sorted_safe &&
            c->model->point_ordered[tx_idx];
        uint16_t sorted_rank = c->workspace->span_exon_rank[tx_idx];
        uint16_t *rank_io = event->interbase
            ? &sorted_rank : &c->workspace->span_exon_rank[tx_idx];
        duckvep_region_state_t region = sorted_span
            ? duckvep_region_classify_span_sorted(
                tx, &c->model->exons, (size_t)tx_idx,
                topology_start1, topology_end1, 0u, 0u,
                (uint8_t)(c->span_sorted_safe == 2),
                rank_io)
            : duckvep_region_classify_span(
                tx, &c->model->exons, (size_t)tx_idx,
                topology_start1, topology_end1, 0u, 0u);
        if (sorted_span) {
            sorted_rank = *rank_io;
            projection_exon_rank = sorted_rank;
        }
        insertion_apply_region_boundaries(
            tx, &c->model->exons, (size_t)tx_idx, event,
            sorted_span, sorted_rank, &region);
        duckvep_splice_state_t splice = sorted_span
            ? duckvep_splice_classify_differing_regions_sorted_with_windows(
                tx, &c->model->exons, (size_t)tx_idx,
                event->feature_start1,
                feature_ref, feature_ref_length,
                feature_alt, feature_alt_length,
                c->options->splice_region_exonic,
                c->options->splice_region_intronic, sorted_rank)
            : duckvep_splice_classify_differing_regions_with_windows(
                tx, &c->model->exons, (size_t)tx_idx,
                event->feature_start1,
                feature_ref, feature_ref_length,
                feature_alt, feature_alt_length,
                c->options->splice_region_exonic,
                c->options->splice_region_intronic);
        duckvep_effect_ctx_fill_classified(
            tx, variant_idx, (size_t)tx_idx,
            topology_start1, topology_end1, &region, &splice, &ectx);
    } else if (kind == DUCKVEP_KIND_SV && event->interbase) {
        duckvep_region_state_t region = duckvep_region_classify_span(
            tx, &c->model->exons, (size_t)tx_idx,
            topology_start1, topology_end1, 0u, 0u);
        duckvep_splice_state_t splice;

        insertion_apply_region_boundaries(
            tx, &c->model->exons, (size_t)tx_idx, event,
            0, UINT16_MAX, &region);
        splice = duckvep_splice_classify_span_with_windows(
            tx, &c->model->exons, (size_t)tx_idx,
            event->feature_start1, event->feature_end1, event->interbase,
            c->options->splice_region_exonic,
            c->options->splice_region_intronic);
        duckvep_effect_ctx_fill_classified(
            tx, variant_idx, (size_t)tx_idx,
            topology_start1, topology_end1, &region, &splice, &ectx);
    } else {
        duckvep_effect_ctx_fill_geometry(
            tx, &c->model->exons, variant_idx, (size_t)tx_idx,
            topology_start1, topology_end1,
            event->feature_start1, event->feature_end1,
            event->interbase,
            c->options->splice_region_exonic,
            c->options->splice_region_intronic, &ectx);
    }
    if (c->delta_route_stats != NULL) {
        c->delta_route_stats->generalized_pair++;
    }
    if (ectx.splice.splice_polypyrimidine &&
        event_overlaps_vep_stretched_exon(c->model, (size_t)tx_idx,
                                          event)) {
        duckvep_effect_ctx_apply_exon_gate(&ectx, 1);
    }
    has_partial_terminal_codon = 0;
    if (c->model->has_seq && ectx.region_state.overlaps_cds &&
        ectx.region_state.overlaps_utr3 && have_feature_alleles &&
        feature_ref_length == feature_alt_length &&
        (tx->flags[tx_idx] & (uint64_t)DUCKVEP_TX_CDS_END_NF) == 0u) {
        has_partial_terminal_codon =
            duckvep_transcript_has_partial_terminal_codon(
                &c->model->seq, (size_t)tx_idx);
    }
    duckvep_effect_ctx_apply_event(tx, &ectx, event);
    if (kind == DUCKVEP_KIND_SV) {
        uint32_t translateable_cds_length =
            c->model->has_seq && c->model->seq.cds_length != NULL
            ? c->model->seq.cds_length[tx_idx] : 0u;
        duckvep_sv_effect_t sv = duckvep_sv_effect_fill(
            tx, &c->model->exons, (size_t)tx_idx,
            translateable_cds_length,
            event, &ectx.region_state);
        duckvep_effect_ctx_apply_sv(&ectx, &sv);
    }

    /* VEP maps an equal-length VariationFeature as one uploaded span. An
     * internal intron does not by itself suppress sequence predicates: when
     * both feature endpoints map to CDS, VEP uses those outer CDS coordinates
     * and retains the internal Mapper::Gap. The peptide path is unavailable
     * only when an outer endpoint is a Gap, including transcript/5'-UTR and
     * ordinary CDS-to-3'-UTR cases. Semantic trimming must not turn either
     * unavailable state into a smaller CDS edit. A partial terminal codon is
     * the 3' exception: VEP can still classify its first mapped coding piece. */
    feature_endpoints_map_to_cdna = 0;
    feature_has_internal_cdna_gap = 0;
    if (have_feature_alleles && !event->interbase &&
        event->feature_start1 != 0u &&
        event->feature_end1 >= event->feature_start1) {
        uint32_t feature_start_cdna;
        uint32_t feature_end_cdna;

        feature_endpoints_map_to_cdna =
            duckvep_project_genomic_to_cdna(
                tx, &c->model->exons, (size_t)tx_idx,
                event->feature_start1, &feature_start_cdna, NULL) &&
            duckvep_project_genomic_to_cdna(
                tx, &c->model->exons, (size_t)tx_idx,
                event->feature_end1, &feature_end_cdna, NULL);
        if (feature_endpoints_map_to_cdna) {
            uint64_t genomic_span =
                (uint64_t)event->feature_end1 -
                (uint64_t)event->feature_start1 + 1u;
            uint64_t cdna_span = feature_start_cdna < feature_end_cdna
                ? (uint64_t)feature_end_cdna -
                      (uint64_t)feature_start_cdna + 1u
                : (uint64_t)feature_start_cdna -
                      (uint64_t)feature_end_cdna + 1u;
            feature_has_internal_cdna_gap = genomic_span > cdna_span;
        }
    }
    feature_mapping_blocks_peptide =
        ectx.region_state.overlaps_cds &&
        (ectx.region_state.within_frameshift_intron ||
         (have_feature_alleles &&
          feature_ref_length == feature_alt_length &&
          !feature_has_internal_cdna_gap &&
          ((ectx.region_state.partial_overlap_feature &&
            ectx.region_state.overlaps_utr5) ||
           ectx.region_state.overlaps_intron ||
           (ectx.region_state.overlaps_utr3 &&
            (tx->flags[tx_idx] &
             (uint64_t)DUCKVEP_TX_CDS_END_NF) == 0u &&
            !has_partial_terminal_codon))));

    /* The symmetric sweep halo admits candidates up to max(up,down); enforce the
     * directional up/downstream window here so an asymmetric config is honored. */
    if (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM)) {
        dist = fwd ? tx->start1[tx_idx] - event->feature_end1
                   : event->feature_start1 - tx->end1[tx_idx];
        if (dist > c->options->upstream_dist) {
            if (!breakend || !breakend_mate_admits) return 1;
            breakend_local_defaults_intergenic = breakend_local_close;
            ectx.pre_bits &= ~DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM);
            ectx.region &= ~(uint32_t)DUCKVEP_REGION_UPSTREAM;
            ectx.region_state.region_mask &=
                ~(uint32_t)DUCKVEP_REGION_UPSTREAM;
        }
    } else if (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM)) {
        dist = fwd ? event->feature_start1 - tx->end1[tx_idx]
                   : tx->start1[tx_idx] - event->feature_end1;
        if (dist > c->options->downstream_dist) {
            if (!breakend || !breakend_mate_admits) return 1;
            breakend_local_defaults_intergenic = breakend_local_close;
            ectx.pre_bits &= ~DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM);
            ectx.region &= ~(uint32_t)DUCKVEP_REGION_DOWNSTREAM;
            ectx.region_state.region_mask &=
                ~(uint32_t)DUCKVEP_REGION_DOWNSTREAM;
        }
    }

    /* Reuse the sorted transcript exon cursor for every later projection, not
     * only CDS mutation. HGVS consumes the same absolute hint for contained
     * exon/UTR edits, while its projector independently verifies that the
     * complete uploaded feature span still fits the exon. */
    if (kind != DUCKVEP_KIND_SV && projection_exon_rank != UINT16_MAX &&
        (c->observer != NULL ||
         (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_CDS)) != 0u)) {
        size_t hinted_exon;

        if (sorted_exon_index_at_rank(
                tx, (size_t)tx_idx, projection_exon_rank, &hinted_exon) &&
            hinted_exon <= UINT32_MAX) {
            uint32_t exon_start = c->model->exons.start1[hinted_exon];
            uint32_t exon_end = c->model->exons.end1[hinted_exon];
            int contained;

            if (event->interbase) {
                uint32_t right1 = duckvep_event_right_flank1(event);
                contained = event->insertion_boundary0 >= exon_start &&
                            right1 <= exon_end;
            } else if (event->ref_diff_length == 0u) {
                contained = 0;
            } else {
                uint64_t semantic_end = (uint64_t)event->start1 +
                    (uint64_t)event->ref_diff_length - 1u;
                contained = event->start1 >= exon_start &&
                            semantic_end <= (uint64_t)exon_end;
            }
            if (contained) projection_exon_hint = (uint32_t)hinted_exon;
        }
    }

    /* Only the CDS bucket needs a sequence delta. The dispatcher takes borrowed
     * model views directly; the annotation sweep does not classify coding effects. */
    memset(&delta, 0, sizeof delta);
    if (c->observer != NULL) {
        memset(&coding_context, 0, sizeof coding_context);
        coding_context_status =
            DUCKVEP_VARIANT_CODING_CONTEXT_UNSUPPORTED_KIND;
    }
    if (kind != DUCKVEP_KIND_SV &&
        (ectx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_CDS)) &&
        !feature_mapping_blocks_peptide) {
        const duckvep_sequence_pool_t *seq = c->model->has_seq ? &c->model->seq : NULL;
        cds_delta_attempted = 1;
        duckvep_sequence_delta_route_t route = DUCKVEP_DELTA_ROUTE_DIRECT;
        if (c->observer != NULL) {
            duckvep_sequence_delta_fill_for_annotation_observed(
                kind, tx, &c->model->exons, seq, c->variants, variant_idx,
                (size_t)tx_idx, pos, tx->strand[tx_idx], c->delta_scratch,
                event, ectx.region, projection_exon_hint,
                c->delta_route_stats != NULL ? &route : NULL, &delta,
                &coding_context, &coding_context_status);
            if (coding_context_status == DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
                coding_context.compatibility_profile =
                    (uint8_t)c->options->compatibility_profile;
            }
        } else {
            duckvep_sequence_delta_fill_for_annotation_trace(
                kind, tx, &c->model->exons, seq, c->variants, variant_idx,
                (size_t)tx_idx, pos, tx->strand[tx_idx], c->delta_scratch,
                event, ectx.region, projection_exon_hint,
                c->delta_route_stats != NULL ? &route : NULL, &delta);
        }
        if (c->delta_route_stats != NULL) {
            if (route == DUCKVEP_DELTA_ROUTE_SIMPLE_INDEL) {
                c->delta_route_stats->simple_indel++;
            } else if (route == DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT) {
                c->delta_route_stats->substitution_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_BOUNDARY_CONTEXT) {
                c->delta_route_stats->boundary_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_DEL_CONTEXT) {
                c->delta_route_stats->del_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_INS_CONTEXT) {
                c->delta_route_stats->ins_context++;
            } else if (route == DUCKVEP_DELTA_ROUTE_INDEL_CONTEXT) {
                c->delta_route_stats->indel_context++;
            } else if (route ==
                       DUCKVEP_DELTA_ROUTE_UPLOADED_FEATURE_CONTEXT) {
                c->delta_route_stats->uploaded_feature_context++;
            }
        }
        duckvep_effect_ctx_apply_delta(&ectx, &delta);
    }

    duckvep_effect_ctx_finalize(&ectx);

    cmask = kind == DUCKVEP_KIND_SV
        ? duckvep_effect_eval_structural(ectx.pre_bits)
        : duckvep_effect_eval(ectx.pre_bits);
    /* StructuralVariationOverlap evaluates every endpoint admitted by its
     * fixed 5000-base rule as a separate overlap allele. If the shifted local
     * endpoint is admitted but its directional predicate was disabled by the
     * caller's narrower upstream/downstream window, that local allele falls
     * back to intergenic even when the mate allele independently contributes
     * feature_truncation. The public row is the union of those allele rows. */
    if (breakend_local_defaults_intergenic) {
        cmask |= DUCKVEP_SO(DUCKVEP_SO_INTERGENIC);
    }
    /* BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences assigns
     * VEP 116's default intergenic consequence when a real overlap allele has
     * exhausted its predicate list without a match. Preserve the transcript
     * candidate here. The adapter separately handles an input for which the
     * model produced no transcript candidate at all. */
    if (cmask == 0u) cmask = DUCKVEP_SO(DUCKVEP_SO_INTERGENIC);
    nmd.prediction = (uint8_t)DUCKVEP_NMD_NOT_APPLICABLE;
    nmd.escape_reasons = 0u;
    if ((cmask & DUCKVEP_NMD_ELIGIBLE_SO_MASK) != 0u) {
        duckvep_nmd_predict(tx, &c->model->exons, (size_t)tx_idx, event,
                            cmask, cds_delta_attempted ? &delta : NULL, &nmd);
    }

    if (c->results->count >= c->results->capacity) {
        c->status = DUCKVEP_ERR_RESULT_FULL;
        return 0;
    }
    row = &c->results->rows[c->results->count++];
    memset(row, 0, sizeof *row);
    row->variant_idx = variant_idx;
    row->tx_idx = tx_idx;
    row->gene_idx = 0u; /* gene grouping is an adapter concern for now */
    row->overlap_object_kind =
        (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
    row->consequence_mask = cmask;
    consequence_flags = duckvep_sequence_delta_consequence_flags(
        &delta, cds_delta_attempted);
    /* delta.start_lost describes the reconstructed sequence. VEP's predicate
     * additionally requires its unshifted complete feature interval to
     * overlap the start codon, and caches that result before HGVS-only
     * placement. Keep the sidecar at predicate level rather than leaking the
     * broader sequence fact into HGVSp. */
    if ((consequence_flags &
         (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_START_LOST) != 0u &&
        !duckvep_project_feature_overlaps_start_codon_unshifted(
            tx, &c->model->exons, (size_t)tx_idx, event)) {
        consequence_flags &=
            ~(uint32_t)DUCKVEP_CONSEQUENCE_FLAG_START_LOST;
    }
    row->flags = consequence_flags;
    row->region_mask = ectx.region;
    row->impact = (uint8_t)duckvep_so_impact(cmask);
    row->sequence_status = cds_delta_attempted
        ? delta.sequence_status
        : (uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE;
    row->nmd_prediction = nmd.prediction;
    row->nmd_escape_reasons = nmd.escape_reasons;
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
    if (c->observer != NULL) {
        duckvep_pair_facts_t facts;

        memset(&transcript_edit, 0, sizeof transcript_edit);
        if (row->region_mask == (uint32_t)DUCKVEP_REGION_UPSTREAM ||
            row->region_mask == (uint32_t)DUCKVEP_REGION_DOWNSTREAM) {
            transcript_edit_status = DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT;
        } else {
            transcript_edit_status =
                duckvep_transcript_edit_project_prepared_hint(
                    tx, &c->model->exons, c->variants, variant_idx, tx_idx,
                    event, projection_exon_hint, &transcript_edit);
        }
        memset(&facts, 0, sizeof facts);
        facts.event = event;
        facts.transcript_edit =
            transcript_edit_status == DUCKVEP_TRANSCRIPT_EDIT_OK
                ? &transcript_edit : NULL;
        facts.delta = cds_delta_attempted ? &delta : NULL;
        facts.coding_context =
            coding_context_status == DUCKVEP_VARIANT_CODING_CONTEXT_OK
                ? &coding_context : NULL;
        facts.projection_exon_hint = projection_exon_hint;
        facts.transcript_edit_status = (uint8_t)transcript_edit_status;
        facts.coding_context_status = (uint8_t)coding_context_status;
        facts.cds_delta_attempted = (uint8_t)cds_delta_attempted;
        facts.coding_context_valid = (uint8_t)(
            coding_context_status == DUCKVEP_VARIANT_CODING_CONTEXT_OK);
        if (!c->observer(c->observer_context, c->variants, row, &facts)) {
            c->status = DUCKVEP_ERR_INTERNAL;
            return 0;
        }
    }
    return 1;
}

static int annotate_interval_feature(uint32_t variant_idx,
                                     uint32_t feature_idx,
                                     struct annotate_ctx *c) {
    const duckvep_interval_feature_model_t *features;
    duckvep_consequence_t *row;
    uint64_t consequence_mask;
    uint8_t feature_kind;

    features = &c->model->interval_features;
    if ((size_t)feature_idx >= features->feature_count) {
        c->status = DUCKVEP_ERR_INTERNAL;
        return 0;
    }
    feature_kind = features->kind[feature_idx];
    consequence_mask = duckvep_effect_eval_interval_feature(
        (duckvep_interval_feature_kind_t)feature_kind,
        &c->events[variant_idx], features->chrom_id[feature_idx],
        features->start1[feature_idx],
        features->end1[feature_idx]);
    /* The sweep uses a coarse outer interval. VEP's minimized/reversed
     * insertion geometry can reject a pair at an interval endpoint. */
    if (consequence_mask == 0u) return 1;
    if (c->results->count >= c->results->capacity) {
        c->status = DUCKVEP_ERR_RESULT_FULL;
        return 0;
    }
    row = &c->results->rows[c->results->count++];
    memset(row, 0, sizeof *row);
    row->variant_idx = variant_idx;
    row->gene_idx = UINT32_MAX;
    row->interval_feature_idx = feature_idx;
    row->overlap_object_kind = feature_kind ==
        (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION
        ? (uint8_t)DUCKVEP_OVERLAP_OBJECT_REGULATORY_REGION
        : (uint8_t)DUCKVEP_OVERLAP_OBJECT_TF_BINDING_SITE;
    row->consequence_mask = consequence_mask;
    row->impact = (uint8_t)duckvep_so_impact(consequence_mask);
    row->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE;
    row->nmd_prediction = (uint8_t)DUCKVEP_NMD_NOT_APPLICABLE;
    row->cdna_pos = -1;
    row->cds_pos = -1;
    row->protein_pos = -1;
    if (c->observer != NULL &&
        !c->observer(c->observer_context, c->variants, row, NULL)) {
        c->status = DUCKVEP_ERR_INTERNAL;
        return 0;
    }
    return 1;
}

static int allele_shape_matches_kind(uint8_t        kind,
                                     uint32_t       pos1,
                                     uint32_t       end1,
                                     const uint8_t *ref,
                                     uint16_t       ref_len,
                                     const uint8_t *alt,
                                     uint16_t       alt_len,
                                     duckvep_event_t *event_out) {
    duckvep_event_t event;

    if (!duckvep_event_prepare_small(pos1, ref, ref_len, alt, alt_len, &event) ||
        event.raw_end1 != end1 || event.kind != kind) {
        return 0;
    }
    if (event_out != NULL) *event_out = event;
    return 1;
}

static duckvep_status_t validate_variant_batch(
    const duckvep_model_t         *model,
    const duckvep_variant_batch_t *variants,
    duckvep_event_t               *events,
    int                            allow_breakends,
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
            !duckvep_sv_geometry_valid((duckvep_sv_type_t)sv_type,
                                       variants->pos1[i], variants->end1[i])) {
            return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_COORD_RANGE,
                        "invalid structural event geometry");
        }
        if (kind == (uint8_t)DUCKVEP_KIND_SV &&
            sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
            if (!allow_breakends) {
                return fail(error, DUCKVEP_ERR_UNSUPPORTED, DVW_ANN_SV_TYPE,
                            "breakends require explicit two-locus candidate pairs");
            }
            if (variants->mate_chrom_id == NULL || variants->mate_pos1 == NULL ||
                variants->pos1[i] == UINT32_MAX || variants->mate_pos1[i] == 0u) {
                return fail(error, DUCKVEP_ERR_INVALID_ARG,
                            DVW_ANN_COORD_RANGE,
                            "breakend is missing a valid local or mate locus");
            }
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
        if (events != NULL && kind == (uint8_t)DUCKVEP_KIND_SV) {
            duckvep_event_load(variants, i, &events[i]);
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
                                           (uint16_t)alen,
                                           events != NULL ? &events[i] : NULL)) {
                return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_ALLELE_RANGE,
                            "variant kind inconsistent with REF/ALT span");
            }
            if (events != NULL) {
                events[i].chrom_id = variants->chrom_id[i];
                events[i].sv_type = (uint8_t)DUCKVEP_SV_NONE;
                events[i].copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
            }
        }
    } else if (events != NULL) {
        /* Geometry-only SNV callers may omit allele storage when the model has
         * no sequence. Preserve their point topology without manufacturing REF
         * bytes that no predicate can inspect. */
        for (i = 0u; i < variants->count; i++) {
            if (variants->variant_kind[i] == (uint8_t)DUCKVEP_KIND_SV) continue;
            duckvep_event_load_raw_interval(variants, i, &events[i]);
            events[i].chrom_id = variants->chrom_id[i];
            events[i].kind = (uint8_t)DUCKVEP_KIND_SNV;
            events[i].ref_diff_length = 1u;
            events[i].alt_diff_length = 1u;
            events[i].sv_type = (uint8_t)DUCKVEP_SV_NONE;
            events[i].copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
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
    duckvep_sweep_cursor_t         transcript_sweep;
    duckvep_sweep_cursor_t         interval_feature_sweep;
    uint32_t                       variant_idx;
    const uint32_t                *tx_indices;
    size_t                         tx_count;
    size_t                         tx_pos;
    const uint32_t                *interval_feature_indices;
    size_t                         interval_feature_count;
    size_t                         interval_feature_pos;
    int                            have_slice;
    int                            done;
    int                            point_sorted_safe;
    int                            span_sorted_safe;
    duckvep_annotation_observer_fn observer;
    void                          *observer_context;
    duckvep_event_t                events[];
};

static duckvep_status_t validate_common_annotate_args(
    const duckvep_model_t         *model,
    const duckvep_variant_batch_t *variants,
    const duckvep_options_t       *options,
    duckvep_workspace_t           *workspace,
    duckvep_event_t               *events,
    int                            allow_breakends,
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
    validation_status = validate_variant_batch(
        model, variants, events, allow_breakends, error);
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
    size_t event_bytes;
    size_t cursor_bytes;

    if (out_cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_CURSOR_NULL_OUT,
                    "out_cursor is NULL");
    }
    *out_cursor = NULL;
    if (variants != NULL &&
        (!size_mul_checked(variants->count, sizeof *cursor->events, &event_bytes) ||
         !size_add_checked(sizeof *cursor, event_bytes, &cursor_bytes))) {
        return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_CURSOR_OOM,
                    "annotate cursor event storage overflow");
    }
    if (variants == NULL) cursor_bytes = sizeof *cursor;
    cursor = (duckvep_annotate_cursor_t *)calloc(1u, cursor_bytes);
    if (cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_CURSOR_OOM,
                    "annotate cursor allocation failed");
    }
    st = validate_common_annotate_args(model, variants, options, workspace,
                                       cursor->events, 0, error);
    if (st != DUCKVEP_OK) {
        free(cursor);
        return st;
    }
    cursor->model = model;
    cursor->variants_view = *variants;
    cursor->variants = &cursor->variants_view;
    cursor->options = options;
    cursor->workspace = workspace;
    cursor->point_sorted_safe = workspace_sorted_runs_begin(
        workspace, variants, cursor->events, &cursor->span_sorted_safe);
    if (workspace->force_generalized_annotation) {
        cursor->point_sorted_safe = 0;
        cursor->span_sorted_safe = 0;
    }
    duckvep_sweep_cursor_init(&cursor->transcript_sweep, cursor->variants,
                              &model->transcripts, options->halo,
                              workspace->active, workspace->active_cap,
                              workspace->candidates, workspace->active_cap);
    cursor->transcript_sweep.events = cursor->events;
    duckvep_interval_sweep_cursor_init(
        &cursor->interval_feature_sweep, cursor->variants,
        model->interval_features.chrom_id,
        model->interval_features.start1,
        model->interval_features.end1,
        model->interval_features.feature_count, 0u,
        workspace->interval_feature_active,
        workspace->interval_feature_active_cap,
        workspace->interval_feature_candidates,
        workspace->interval_feature_active_cap);
    cursor->interval_feature_sweep.events = cursor->events;
    if (cursor->transcript_sweep.status != DUCKVEP_OK ||
        cursor->interval_feature_sweep.status != DUCKVEP_OK) {
        duckvep_status_t sweep_status =
            cursor->transcript_sweep.status != DUCKVEP_OK
            ? cursor->transcript_sweep.status
            : cursor->interval_feature_sweep.status;
        st = fail(error, sweep_status, DVW_ANN_SWEEP,
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

void duckvep_annotate_cursor_set_observer(
    duckvep_annotate_cursor_t      *cursor,
    duckvep_annotation_observer_fn  observer,
    void                           *observer_context) {

    if (cursor == NULL) return;
    cursor->observer = observer;
    cursor->observer_context = observer_context;
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
    if (cursor->have_slice || cursor->done ||
        cursor->transcript_sweep.vi != 0u) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_SWEEP,
                    "annotate cursor can only be seeded before its first fill");
    }
    if (!duckvep_sweep_cursor_seed(&cursor->transcript_sweep,
                                   transcript_indices,
                                   transcript_count)) {
        return fail(error, cursor->transcript_sweep.status, DVW_ANN_SWEEP,
                    "candidate sweep seed is invalid");
    }
    return DUCKVEP_OK;
}

duckvep_status_t duckvep_annotate_cursor_seed_interval_features(
    duckvep_annotate_cursor_t *cursor,
    const uint32_t            *feature_indices,
    size_t                     feature_count,
    duckvep_error_t           *error) {

    if (cursor == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_CURSOR_NULL,
                    "annotate cursor is NULL");
    }
    if (cursor->have_slice || cursor->done ||
        cursor->interval_feature_sweep.vi != 0u) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_SWEEP,
                    "interval-feature sweep can only be seeded before its first fill");
    }
    if (!duckvep_sweep_cursor_seed(&cursor->interval_feature_sweep,
                                   feature_indices, feature_count)) {
        return fail(error, cursor->interval_feature_sweep.status,
                    DVW_ANN_SWEEP,
                    "interval-feature sweep seed is invalid");
    }
    return DUCKVEP_OK;
}

static duckvep_status_t annotate_cursor_fill(
    duckvep_annotate_cursor_t      *cursor,
    duckvep_result_builder_t       *results,
    int                             pause_when_full,
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
    ctx.events = cursor->events;
    ctx.options = cursor->options;
    ctx.workspace = cursor->workspace;
    ctx.delta_scratch = duckvep_workspace_delta_scratch(cursor->workspace);
    ctx.delta_route_stats = cursor->workspace->delta_route_stats_enabled
                              ? &cursor->workspace->delta_route_stats
                              : NULL;
    ctx.results = results;
    ctx.observer = cursor->observer;
    ctx.observer_context = cursor->observer_context;
    ctx.status = DUCKVEP_OK;
    ctx.point_sorted_safe = cursor->point_sorted_safe;
    ctx.span_sorted_safe = cursor->span_sorted_safe;
    ctx.prepared_variant_idx = UINT32_MAX;

    for (;;) {
        if (!cursor->have_slice) {
            uint32_t transcript_variant_idx, feature_variant_idx;
            int transcript_rc, feature_rc;

            transcript_rc = duckvep_sweep_cursor_next(
                &cursor->transcript_sweep, &transcript_variant_idx,
                &cursor->tx_indices, &cursor->tx_count);
            if (cursor->model->interval_features.feature_count != 0u) {
                feature_rc = duckvep_sweep_cursor_next(
                    &cursor->interval_feature_sweep, &feature_variant_idx,
                    &cursor->interval_feature_indices,
                    &cursor->interval_feature_count);
            } else {
                feature_rc = transcript_rc;
                feature_variant_idx = transcript_variant_idx;
                cursor->interval_feature_indices = NULL;
                cursor->interval_feature_count = 0u;
            }
            if (transcript_rc < 0 || feature_rc < 0) {
                duckvep_status_t sweep_status = transcript_rc < 0
                    ? cursor->transcript_sweep.status
                    : cursor->interval_feature_sweep.status;
                return fail(error, sweep_status, DVW_ANN_SWEEP,
                            "candidate sweep failed");
            }
            if (transcript_rc != feature_rc ||
                (transcript_rc != 0 &&
                 transcript_variant_idx != feature_variant_idx)) {
                return fail(error, DUCKVEP_ERR_INTERNAL, DVW_ANN_SWEEP,
                            "transcript and interval-feature sweeps diverged");
            }
            if (transcript_rc == 0) {
                cursor->done = 1;
                cursor->have_slice = 0;
                return DUCKVEP_OK;
            }
            cursor->variant_idx = transcript_variant_idx;
            cursor->have_slice = 1;
            cursor->tx_pos = 0u;
            cursor->interval_feature_pos = 0u;
        }

        while (cursor->tx_pos < cursor->tx_count) {
            if (!annotate_pair(
                    cursor->variant_idx,
                    cursor->tx_indices[cursor->tx_pos], &ctx)) {
                if (ctx.status == DUCKVEP_ERR_RESULT_FULL) {
                    return fail(error, DUCKVEP_ERR_RESULT_FULL, DVW_ANN_RESULT_FULL,
                                "result builder capacity exhausted; cursor paused");
                }
                return fail(error, ctx.status, DVW_ANN_SWEEP, "annotate cursor failed");
            }
            cursor->tx_pos++;
            if (pause_when_full && results->count >= results->capacity) {
                return fail(error, DUCKVEP_ERR_RESULT_FULL, DVW_ANN_RESULT_FULL,
                            "result builder capacity exhausted; cursor paused");
            }
        }
        while (cursor->interval_feature_pos <
               cursor->interval_feature_count) {
            if (!annotate_interval_feature(
                    cursor->variant_idx,
                    cursor->interval_feature_indices[
                        cursor->interval_feature_pos], &ctx)) {
                if (ctx.status == DUCKVEP_ERR_RESULT_FULL) {
                    return fail(error, DUCKVEP_ERR_RESULT_FULL,
                                DVW_ANN_RESULT_FULL,
                                "result builder capacity exhausted; cursor paused");
                }
                return fail(error, ctx.status, DVW_ANN_SWEEP,
                            "interval-feature annotation failed");
            }
            cursor->interval_feature_pos++;
            if (pause_when_full && results->count >= results->capacity) {
                return fail(error, DUCKVEP_ERR_RESULT_FULL,
                            DVW_ANN_RESULT_FULL,
                            "result builder capacity exhausted; cursor paused");
            }
        }
        cursor->have_slice = 0;
    }
}

duckvep_status_t duckvep_annotate_cursor_fill(
    duckvep_annotate_cursor_t      *cursor,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error) {

    return annotate_cursor_fill(cursor, results, 1, error);
}

duckvep_status_t duckvep_annotate_tile(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error) {
    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_status_t status;

    status = validate_result_builder_for_append(results, error);
    if (status != DUCKVEP_OK) return status;
    status = duckvep_annotate_cursor_open(model, variants, options, workspace,
                                          &cursor, error);
    if (status != DUCKVEP_OK) return status;
    /* The one-shot API historically returns RESULT_FULL only if another row
     * actually needs storage. A resumable cursor instead pauses as soon as its
     * output chunk fills. Both paths share event preparation and annotation;
     * only their backpressure contract differs. */
    status = annotate_cursor_fill(cursor, results, 0, error);
    duckvep_annotate_cursor_close(cursor);
    return status;
}

static duckvep_status_t annotate_explicit_pairs(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const uint32_t                 *pair_variant_idx,
    const uint32_t                 *pair_object_idx,
    size_t                          pair_count,
    int                             interval_features,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    duckvep_result_builder_t       *results,
    duckvep_error_t                *error) {
    struct annotate_ctx ctx;
    duckvep_event_t *events = NULL;
    duckvep_status_t status;
    size_t event_bytes = 0u;
    size_t object_count;
    size_t pair;
    uint32_t columns_where = interval_features
        ? DVW_ANN_FEATURE_PAIR_COLUMNS : DVW_ANN_PAIR_COLUMNS;
    uint32_t order_where = interval_features
        ? DVW_ANN_FEATURE_PAIR_ORDER : DVW_ANN_PAIR_ORDER;

    status = validate_result_builder_for_append(results, error);
    if (status != DUCKVEP_OK) return status;
    if (pair_count != 0u &&
        (pair_variant_idx == NULL || pair_object_idx == NULL)) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, columns_where,
                    interval_features
                        ? "interval-feature pair view has null required columns"
                        : "candidate pair view has null required columns");
    }
    if (variants != NULL && variants->count != 0u) {
        if (!size_mul_checked(variants->count, sizeof *events, &event_bytes)) {
            return fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_CURSOR_OOM,
                        "candidate-pair event storage overflow");
        }
        events = (duckvep_event_t *)calloc(1u, event_bytes);
        if (events == NULL) {
            return fail(error, DUCKVEP_ERR_INTERNAL, DVW_CURSOR_OOM,
                        "candidate-pair event allocation failed");
        }
    }
    status = validate_common_annotate_args(
        model, variants, options, workspace, events, 1, error);
    if (status != DUCKVEP_OK) {
        free(events);
        return status;
    }
    object_count = interval_features
        ? model->interval_features.feature_count
        : model->transcripts.transcript_count;
    for (pair = 0u; pair < pair_count; pair++) {
        uint32_t variant_idx = pair_variant_idx[pair];
        uint32_t object_idx = pair_object_idx[pair];

        if ((size_t)variant_idx >= variants->count ||
            (size_t)object_idx >= object_count ||
            (pair != 0u &&
             (variant_idx < pair_variant_idx[pair - 1u] ||
              (variant_idx == pair_variant_idx[pair - 1u] &&
               object_idx <= pair_object_idx[pair - 1u])))) {
            free(events);
            return fail(error, DUCKVEP_ERR_INVALID_ARG, order_where,
                        interval_features
                            ? "interval-feature pairs are out of range, unsorted, or duplicated"
                            : "candidate pairs are out of range, unsorted, or duplicated");
        }
    }

    workspace_point_cursor_reset(workspace);
    workspace_span_cursor_reset(workspace);
    ctx.model = model;
    ctx.variants = variants;
    ctx.events = events;
    ctx.options = options;
    ctx.workspace = workspace;
    ctx.delta_scratch = duckvep_workspace_delta_scratch(workspace);
    ctx.delta_route_stats = workspace->delta_route_stats_enabled
                              ? &workspace->delta_route_stats : NULL;
    ctx.results = results;
    ctx.observer = NULL;
    ctx.observer_context = NULL;
    ctx.status = DUCKVEP_OK;
    ctx.point_sorted_safe = 0;
    ctx.span_sorted_safe = 0;
    ctx.prepared_variant_idx = UINT32_MAX;
    for (pair = 0u; pair < pair_count; pair++) {
        int ok = interval_features
            ? annotate_interval_feature(pair_variant_idx[pair],
                                        pair_object_idx[pair], &ctx)
            : annotate_pair(pair_variant_idx[pair], pair_object_idx[pair],
                            &ctx);
        if (!ok) {
            status = ctx.status == DUCKVEP_ERR_RESULT_FULL
                ? fail(error, DUCKVEP_ERR_RESULT_FULL,
                       DVW_ANN_RESULT_FULL,
                       "result builder capacity exhausted")
                : fail(error, ctx.status, order_where,
                       interval_features
                           ? "interval-feature candidate annotation failed"
                           : "candidate-pair annotation failed");
            free(events);
            return status;
        }
    }
    free(events);
    return DUCKVEP_OK;
}

duckvep_status_t duckvep_annotate_pairs(
    const duckvep_model_t            *model,
    const duckvep_variant_batch_t    *variants,
    const duckvep_candidate_pairs_t  *pairs,
    const duckvep_options_t          *options,
    duckvep_workspace_t              *workspace,
    duckvep_result_builder_t         *results,
    duckvep_error_t                  *error) {

    if (pairs == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_ANN_PAIR_COLUMNS,
                    "candidate pair view is null");
    }
    return annotate_explicit_pairs(
        model, variants, pairs->variant_idx, pairs->tx_idx, pairs->count, 0,
        options, workspace, results, error);
}

duckvep_status_t duckvep_annotate_interval_feature_pairs(
    const duckvep_model_t                  *model,
    const duckvep_variant_batch_t          *variants,
    const duckvep_interval_feature_pairs_t *pairs,
    const duckvep_options_t                *options,
    duckvep_workspace_t                    *workspace,
    duckvep_result_builder_t               *results,
    duckvep_error_t                        *error) {

    if (pairs == NULL) {
        return fail(error, DUCKVEP_ERR_INVALID_ARG,
                    DVW_ANN_FEATURE_PAIR_COLUMNS,
                    "interval-feature pair view is null");
    }
    return annotate_explicit_pairs(
        model, variants, pairs->variant_idx, pairs->feature_idx, pairs->count,
        1, options, workspace, results, error);
}

/* Model preparation is deliberately kept after the annotation functions. It
 * runs once per immutable model, while the functions above run for every
 * variant/transcript pair. Keeping cold preparation code out of their layout
 * also makes instruction-cache measurements stable when derived model fields
 * are added. */
static DUCKVEP_NOINLINE duckvep_status_t model_prepare_derived_layout(
    struct duckvep_model              *model,
    const duckvep_transcript_model_t  *transcripts,
    const duckvep_exon_model_t        *exons,
    duckvep_error_t                   *error) {

    size_t t;

    if (transcripts->transcript_count == 0u) return DUCKVEP_OK;
    model->cds_cdna_start1 = (uint32_t *)calloc(
        transcripts->transcript_count, sizeof *model->cds_cdna_start1);
    model->cds_cdna_end1 = (uint32_t *)calloc(
        transcripts->transcript_count, sizeof *model->cds_cdna_end1);
    model->cds_start_exon_index = (uint32_t *)calloc(
        transcripts->transcript_count, sizeof *model->cds_start_exon_index);
    model->cds_phase_offset = (uint8_t *)calloc(
        transcripts->transcript_count, sizeof *model->cds_phase_offset);
    model->point_ordered = (uint8_t *)calloc(
        transcripts->transcript_count, sizeof *model->point_ordered);
    model->has_frameshift_intron = (uint8_t *)calloc(
        transcripts->transcript_count, sizeof *model->has_frameshift_intron);
    if (model->cds_cdna_start1 == NULL || model->cds_cdna_end1 == NULL ||
        model->cds_start_exon_index == NULL ||
        model->cds_phase_offset == NULL || model->point_ordered == NULL ||
        model->has_frameshift_intron == NULL) {
        return fail(error, DUCKVEP_ERR_INTERNAL, DVW_MODEL_OOM,
                    "model derived-layout alloc failed");
    }

    for (t = 0u; t < transcripts->transcript_count; t++) {
        uint32_t coding_start_genomic;
        uint32_t coding_end_genomic;
        size_t coding_start_exon = 0u;

        model->point_ordered[t] = (uint8_t)transcript_point_ordered(
            transcripts, exons, t);
        model->has_frameshift_intron[t] =
            (uint8_t)transcript_has_frameshift_intron(
                transcripts, exons, t);
        if (transcripts->cds_start1[t] == 0u) continue;
        coding_start_genomic = transcripts->strand[t] > 0
            ? transcripts->cds_start1[t] : transcripts->cds_end1[t];
        coding_end_genomic = transcripts->strand[t] > 0
            ? transcripts->cds_end1[t] : transcripts->cds_start1[t];
        if (!model_genomic_to_cdna(
                transcripts, exons, t, coding_start_genomic,
                &model->cds_cdna_start1[t], &coding_start_exon) ||
            !model_genomic_to_cdna(
                transcripts, exons, t, coding_end_genomic,
                &model->cds_cdna_end1[t], NULL) ||
            coding_start_exon > UINT32_MAX ||
            model->cds_cdna_end1[t] < model->cds_cdna_start1[t]) {
            /* Topology-only models may describe a CDS span without the
             * complete exon/cDNA relation needed for sequence projection.
             * Keep that accepted model usable and leave this transcript on
             * the exhaustive path. Sequence-backed coding models proved this
             * relation during validation. */
            model->cds_cdna_start1[t] = 0u;
            model->cds_cdna_end1[t] = 0u;
            continue;
        }
        model->cds_start_exon_index[t] = (uint32_t)coding_start_exon;
        if (exons->phase != NULL && exons->phase[coding_start_exon] > 0) {
            model->cds_phase_offset[t] =
                (uint8_t)exons->phase[coding_start_exon];
        }
    }
    model->transcripts.cds_cdna_start1 = model->cds_cdna_start1;
    model->transcripts.cds_cdna_end1 = model->cds_cdna_end1;
    model->transcripts.cds_start_exon_index =
        model->cds_start_exon_index;
    model->transcripts.cds_phase_offset = model->cds_phase_offset;
    return DUCKVEP_OK;
}
