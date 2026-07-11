/* duckvep_sv.c — VEP-shaped structural/CNV predicates. */
#include "duckvep_sv.h"

#include <string.h>

int duckvep_sv_metadata_valid(duckvep_sv_type_t sv_type,
                              duckvep_copy_change_t copy_change) {
    if (sv_type <= DUCKVEP_SV_NONE || sv_type > DUCKVEP_SV_BREAKEND ||
        copy_change < DUCKVEP_COPY_CHANGE_UNKNOWN ||
        copy_change > DUCKVEP_COPY_CHANGE_GAIN) {
        return 0;
    }

    switch (sv_type) {
    case DUCKVEP_SV_DELETION:
        return copy_change == DUCKVEP_COPY_CHANGE_UNKNOWN ||
               copy_change == DUCKVEP_COPY_CHANGE_LOSS;
    case DUCKVEP_SV_DUPLICATION:
    case DUCKVEP_SV_TANDEM_DUPLICATION:
        return copy_change == DUCKVEP_COPY_CHANGE_UNKNOWN ||
               copy_change == DUCKVEP_COPY_CHANGE_GAIN;
    case DUCKVEP_SV_INSERTION:
        return copy_change != DUCKVEP_COPY_CHANGE_LOSS;
    case DUCKVEP_SV_INVERSION:
    case DUCKVEP_SV_BREAKEND:
        return copy_change == DUCKVEP_COPY_CHANGE_UNKNOWN ||
               copy_change == DUCKVEP_COPY_CHANGE_NEUTRAL;
    case DUCKVEP_SV_UNKNOWN:
    case DUCKVEP_SV_CNV:
        return 1;
    case DUCKVEP_SV_NONE:
    default:
        return 0;
    }
}

duckvep_sv_effect_t duckvep_sv_effect_fill(
    const duckvep_event_t        *event,
    const duckvep_region_state_t *region) {

    duckvep_sv_effect_t st;
    memset(&st, 0, sizeof st);
    if (event == NULL || region == NULL ||
        event->kind != (uint8_t)DUCKVEP_KIND_SV) {
        return st;
    }

    st.copy_number_gain =
        (uint8_t)(event->copy_change == (uint8_t)DUCKVEP_COPY_CHANGE_GAIN ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_DUPLICATION ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_TANDEM_DUPLICATION);
    st.copy_number_loss =
        (uint8_t)(event->copy_change == (uint8_t)DUCKVEP_COPY_CHANGE_LOSS ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_DELETION);
    st.deletion =
        (uint8_t)(event->sv_type == (uint8_t)DUCKVEP_SV_DELETION ||
                  st.copy_number_loss);
    st.insertion =
        (uint8_t)(event->sv_type == (uint8_t)DUCKVEP_SV_INSERTION ||
                  st.copy_number_gain);
    st.chromosome_breakpoint =
        (uint8_t)(event->sv_type == (uint8_t)DUCKVEP_SV_BREAKEND);

    st.feature_ablation =
        (uint8_t)(region->complete_overlap_feature && st.deletion);
    st.feature_amplification =
        (uint8_t)(region->complete_overlap_feature && st.copy_number_gain);
    st.feature_elongation =
        (uint8_t)(region->within_cdna && region->complete_within_feature &&
                  (st.copy_number_gain || st.insertion));

    /* VEP treats a chromosome breakpoint specially: a local breakend inside the
     * feature truncates it even when the point is intronic. Other losses/deletions
     * must overlap cDNA and be partial-overlap or wholly within the feature. */
    if (st.chromosome_breakpoint) {
        st.feature_truncation = region->within_feature;
    } else {
        st.feature_truncation =
            (uint8_t)(region->within_cdna &&
                      (region->partial_overlap_feature ||
                       region->complete_within_feature) &&
                      (st.copy_number_loss || st.deletion));
    }
    return st;
}
