/* duckvep_model_internal.h — prepared immutable model views (INTERNAL). */
#ifndef DUCKVEP_MODEL_INTERNAL_H
#define DUCKVEP_MODEL_INTERNAL_H

#include "duckvep_delta.h"

#ifdef __cplusplus
extern "C" {
#endif

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
