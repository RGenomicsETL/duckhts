/*
 * duckvep_sv.h — structural-variant/CNV predicate facts (INTERNAL).
 *
 * The adapter supplies an explicit structural operation and, for generic CNVs, a
 * sample-aware copy-change direction. This module maps those semantics plus
 * exact span/topology facts to the predicates used by VEP's consequence program.
 */
#ifndef DUCKVEP_SV_H
#define DUCKVEP_SV_H

#include "duckvep_classify.h"
#include "duckvep_event.h"

#include <stdint.h>

typedef struct duckvep_sv_effect {
    uint8_t deletion;
    uint8_t insertion;
    uint8_t copy_number_gain;
    uint8_t copy_number_loss;
    uint8_t chromosome_breakpoint;
    uint8_t feature_ablation;
    uint8_t feature_amplification;
    uint8_t feature_elongation;
    uint8_t feature_truncation;
    uint8_t start_lost;
    uint8_t start_retained;
    uint8_t stop_lost;
    uint8_t frameshift;
    uint8_t inframe_deletion;
} duckvep_sv_effect_t;

duckvep_sv_effect_t duckvep_sv_effect_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          translateable_cds_length,
    const duckvep_event_t            *event,
    const duckvep_region_state_t     *region);

#endif /* DUCKVEP_SV_H */
