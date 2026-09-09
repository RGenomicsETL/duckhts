/*
 * duckvep_annotation_internal.h — synchronous annotation-row observation.
 *
 * This is private adapter plumbing, not part of the stable kernel ABI.  It lets
 * a consumer render a derivative such as HGVS while the exact consequence
 * coding context still borrows worker scratch.  No pointer in the facts may be
 * retained after the callback returns.
 */
#ifndef DUCKVEP_ANNOTATION_INTERNAL_H
#define DUCKVEP_ANNOTATION_INTERNAL_H

#include "duckvep_delta.h"
#include "duckvep_event.h"
#include "duckvep_kernel.h"
#include "duckvep_transcript_edit.h"

/* Canonical callback-scoped facts for one emitted transcript pair.  All
 * pointers borrow stack or worker scratch and expire when the observer
 * returns.  HGVS and future phased collection consume this object instead of
 * independently rediscovering prepared-event or transcript projection state. */
typedef struct duckvep_pair_facts {
    const duckvep_event_t            *event;
    const duckvep_transcript_edit_t  *transcript_edit;
    const duckvep_sequence_delta_t   *delta;
    const duckvep_coding_context_t   *coding_context;
    uint32_t                          projection_exon_hint;
    uint8_t                           transcript_edit_status;
    uint8_t                           coding_context_status;
    uint8_t                           cds_delta_attempted;
    uint8_t                           coding_context_valid;
} duckvep_pair_facts_t;

/* Invoked exactly once for every emitted row, including regulation/motif rows.
 * `facts` is NULL for non-transcript objects. Returning zero aborts the vector;
 * the adapter remains responsible for its own reason-coded error text. */
typedef int (*duckvep_annotation_observer_fn)(
    void                            *observer_context,
    const duckvep_variant_batch_t   *variants,
    const duckvep_consequence_t     *row,
    const duckvep_pair_facts_t       *facts);

DUCKVEP_INTERNAL_API void duckvep_annotate_cursor_set_observer(
    duckvep_annotate_cursor_t          *cursor,
    duckvep_annotation_observer_fn      observer,
    void                               *observer_context);

/* Observe one explicitly selected transcript for a one-row literal-allele
 * batch and its already prepared duckvep_event_prepare_small geometry. The
 * caller owns the matching immutable source bytes and initializes the event's
 * chromosome and non-structural metadata. Full-span source replacements are
 * not this geometry: normalize their HGVS interpretation separately. Uses
 * VEP-116 defaults with zero upstream/downstream reach and the same consequence
 * machinery as the cursor. No allele normalization, sweep state or heap
 * storage is constructed here. The callback may be omitted by transcript
 * admission and must consume its borrowed facts before returning. `projected`
 * may borrow a successful physical CDS edit for this exact prepared allele and
 * immutable model/transcript. Feature and shifted-HGVS coordinates stay distinct. */
DUCKVEP_INTERNAL_API duckvep_status_t duckvep_annotate_pair_observed(
    const duckvep_model_t *model, const duckvep_variant_batch_t *variants,
    const duckvep_event_t *event, uint32_t transcript_index, duckvep_delta_scratch_t *scratch,
    const duckvep_haplotype_edit_t *projected,
    duckvep_annotation_observer_fn observer, void *observer_context,
    duckvep_error_t *error);

#endif /* DUCKVEP_ANNOTATION_INTERNAL_H */
