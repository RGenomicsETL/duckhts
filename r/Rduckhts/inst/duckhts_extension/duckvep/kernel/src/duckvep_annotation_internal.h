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

#endif /* DUCKVEP_ANNOTATION_INTERNAL_H */
