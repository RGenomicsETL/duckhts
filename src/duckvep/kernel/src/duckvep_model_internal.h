/* duckvep_model_internal.h — prepared immutable model views (INTERNAL). */
#ifndef DUCKVEP_MODEL_INTERNAL_H
#define DUCKVEP_MODEL_INTERNAL_H

#include "duckvep_delta.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Stable model-layout failure-site IDs shared by construction and projection. */
enum {
    DVW_MODEL_NULL_VIEW = 11u,
    DVW_MODEL_EXON_RANGE = 13u,
    DVW_MODEL_CDS_RANGE = 14u,
    DVW_MODEL_TX_LAYOUT = 66u,
    DVW_MODEL_EXON_LAYOUT = 67u,
    DVW_MODEL_CDNA_LAYOUT = 68u,
    DVW_MODEL_PHASE = 69u,
    DVW_MODEL_CDS_PROJECTION = 70u
};

/* Validate one transcript's span, ordered exon/cDNA slice, phases and CDS
 * endpoints. Does not allocate, inspect sequence or trust projection caches.
 * Model construction and borrowed-view replay share these layout checks. */
DUCKVEP_INTERNAL_API duckvep_status_t duckvep_model_validate_transcript_layout(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t *exons, size_t transcript, duckvep_error_t *error);

/* The public constructor borrows the caller's model arrays, then attaches
 * immutable derived projection caches to its private transcript view. Kernel
 * consumers outside duckvep_kernel.c must use this prepared view after open;
 * retaining the pre-open transcript struct would silently bypass those caches
 * and repeat both CDS-endpoint projections for every coding annotation. */
DUCKVEP_INTERNAL_API const duckvep_transcript_model_t *
duckvep_model_prepared_transcripts(const duckvep_model_t *model);

DUCKVEP_INTERNAL_API const duckvep_exon_model_t *
duckvep_model_prepared_exons(const duckvep_model_t *model);

DUCKVEP_INTERNAL_API const duckvep_sequence_pool_t *
duckvep_model_prepared_sequences(const duckvep_model_t *model);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_MODEL_INTERNAL_H */
