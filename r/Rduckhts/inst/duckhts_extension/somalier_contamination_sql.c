/* Query-owned CHARR reduction over validated, ordered panel evidence. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckhts_somalier.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define CONTAM_SQL_MAX_SITES UINT64_C(100000000)
#define CONTAM_METHOD_VERSION "somalier-0.3.4-duckhts-1.5.2"

enum charr_input {
    CH_IN_SAMPLE = 0, CH_IN_ASSEMBLY, CH_IN_PANEL, CH_IN_FREQUENCY,
    CH_IN_SITE_INDEX, CH_IN_SITE_COUNT, CH_IN_A, CH_IN_B, CH_IN_OTHER,
    CH_IN_POPULATION_B_AF, CH_IN_MIN_DEPTH, CH_IN_MAX_DEPTH,
    CH_IN_HOM_MINOR_RATE, CH_IN_HOM_TAIL_ALPHA, CH_IN_MAX_SITES,
    CH_IN_COUNT
};

enum charr_output {
    CH_OUT_SAMPLE = 0, CH_OUT_ASSEMBLY, CH_OUT_PANEL, CH_OUT_FREQUENCY,
    CH_OUT_METHOD, CH_OUT_STATUS, CH_OUT_SITE_COUNT, CH_OUT_OBSERVED,
    CH_OUT_UNAVAILABLE, CH_OUT_USABLE, CH_OUT_HOM_A, CH_OUT_HOM_B,
    CH_OUT_ESTIMATE, CH_OUT_MIN_DEPTH, CH_OUT_MAX_DEPTH,
    CH_OUT_HOM_MINOR_RATE, CH_OUT_HOM_TAIL_ALPHA, CH_OUT_MAX_SITES,
    CH_OUT_COUNT
};

static const char *charr_output_names[CH_OUT_COUNT] = {
    "sample_id", "assembly", "panel_sha256", "frequency_sha256",
    "method_version", "status", "site_count", "observed_sites",
    "unavailable_sites", "usable_sites", "usable_hom_a", "usable_hom_b",
    "estimate", "min_depth", "max_depth", "hom_minor_rate",
    "hom_tail_alpha", "max_sites"
};

typedef struct contamination_identity {
    char *sample;
    char *assembly;
    char *panel;
    char *frequency;
    uint32_t sample_len;
    uint32_t assembly_len;
    uint32_t panel_len;
    uint32_t frequency_len;
} contamination_identity_t;

typedef struct contamination_state {
    contamination_identity_t identity;
    duckhts_somalier_contamination_settings_t filter;
    duckhts_somalier_charr_accumulator_t charr;
    uint64_t *seen;
    size_t site_count;
    size_t word_count;
    uint64_t observed_sites;
    uint64_t unavailable_sites;
    uint8_t initialized;
} contamination_state_t;

static bool row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return validity == NULL || duckdb_validity_row_is_valid(validity, row);
}

static bool lowercase_sha256(duckdb_string_t *value) {
    const char *text = duckdb_string_t_data(value);
    if (duckdb_string_t_length(*value) != 64u) return false;
    for (unsigned i = 0u; i < 64u; i++) {
        char c = text[i];
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    }
    return true;
}

static char *copy_string(duckdb_string_t *value, uint32_t *length) {
    char *copy;
    *length = duckdb_string_t_length(*value);
    if (*length == 0u || *length > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES) return NULL;
    copy = duckdb_malloc((size_t)*length + 1u);
    if (copy == NULL) return NULL;
    memcpy(copy, duckdb_string_t_data(value), *length);
    copy[*length] = '\0';
    return copy;
}

static bool same_string(duckdb_string_t *value, const char *copy, uint32_t length) {
    return copy != NULL && duckdb_string_t_length(*value) == length &&
        memcmp(duckdb_string_t_data(value), copy, length) == 0;
}

static bool same_copy(const char *left, uint32_t left_len,
                      const char *right, uint32_t right_len) {
    return left_len == right_len && memcmp(left, right, left_len) == 0;
}

static void identity_release(contamination_identity_t *identity) {
    duckdb_free(identity->sample);
    duckdb_free(identity->assembly);
    duckdb_free(identity->panel);
    duckdb_free(identity->frequency);
    memset(identity, 0, sizeof(*identity));
}

static void state_release(contamination_state_t *state) {
    if (state == NULL) return;
    identity_release(&state->identity);
    duckdb_free(state->seen);
    memset(state, 0, sizeof(*state));
}

static idx_t state_size(duckdb_function_info info) {
    (void)info;
    return sizeof(contamination_state_t);
}

static void state_init(duckdb_function_info info, duckdb_aggregate_state state) {
    (void)info;
    if (state != NULL) memset(state, 0, sizeof(contamination_state_t));
}

static void state_destroy(duckdb_aggregate_state *states, idx_t count) {
    if (states == NULL) return;
    for (idx_t i = 0u; i < count; i++) {
        state_release((contamination_state_t *)states[i]);
    }
}

static bool identity_open(contamination_identity_t *identity,
                          duckdb_string_t *sample,
                          duckdb_string_t *assembly, duckdb_string_t *panel,
                          duckdb_string_t *frequency) {
    if (duckdb_string_t_length(*sample) == 0u ||
        duckdb_string_t_length(*assembly) == 0u ||
        !lowercase_sha256(panel) || !lowercase_sha256(frequency)) return false;
    if (duckdb_string_t_length(*sample) > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES ||
        duckdb_string_t_length(*assembly) > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES) return false;
    identity->sample = copy_string(sample, &identity->sample_len);
    identity->assembly = copy_string(assembly, &identity->assembly_len);
    identity->panel = copy_string(panel, &identity->panel_len);
    identity->frequency = copy_string(frequency, &identity->frequency_len);
    return identity->sample != NULL && identity->assembly != NULL &&
        identity->panel != NULL && identity->frequency != NULL;
}

static bool identity_matches(const contamination_identity_t *identity,
                             duckdb_string_t *sample,
                             duckdb_string_t *assembly, duckdb_string_t *panel,
                             duckdb_string_t *frequency) {
    return same_string(sample, identity->sample, identity->sample_len) &&
        same_string(assembly, identity->assembly, identity->assembly_len) &&
        same_string(panel, identity->panel, identity->panel_len) &&
        same_string(frequency, identity->frequency, identity->frequency_len);
}

static bool identity_same(const contamination_identity_t *left,
                          const contamination_identity_t *right) {
    return same_copy(left->sample, left->sample_len, right->sample, right->sample_len) &&
        same_copy(left->assembly, left->assembly_len,
                  right->assembly, right->assembly_len) &&
        same_copy(left->panel, left->panel_len, right->panel, right->panel_len) &&
        same_copy(left->frequency, left->frequency_len,
                  right->frequency, right->frequency_len);
}

static bool identity_clone(contamination_identity_t *target,
                           const contamination_identity_t *source) {
    const char *strings[4] = {source->sample, source->assembly,
                              source->panel, source->frequency};
    char **targets[4] = {&target->sample, &target->assembly,
                         &target->panel, &target->frequency};
    uint32_t lengths[4] = {source->sample_len, source->assembly_len,
                           source->panel_len, source->frequency_len};
    for (unsigned i = 0u; i < 4u; i++) {
        if (strings[i] == NULL) continue;
        *targets[i] = duckdb_malloc((size_t)lengths[i] + 1u);
        if (*targets[i] == NULL) return false;
        memcpy(*targets[i], strings[i], (size_t)lengths[i] + 1u);
    }
    target->sample_len = source->sample_len;
    target->assembly_len = source->assembly_len;
    target->panel_len = source->panel_len;
    target->frequency_len = source->frequency_len;
    return true;
}

static bool state_open(contamination_state_t *state,
                       duckdb_string_t *sample,
                       duckdb_string_t *assembly, duckdb_string_t *panel,
                       duckdb_string_t *frequency, uint64_t site_count,
                       const duckhts_somalier_contamination_settings_t *filter) {
    size_t words;
    if (site_count == 0u || site_count > CONTAM_SQL_MAX_SITES ||
        site_count > filter->max_sites || site_count > SIZE_MAX) return false;
    words = (size_t)(site_count / 64u + (site_count % 64u != 0u));
    if (words > SIZE_MAX / sizeof(uint64_t)) return false;
    if (!identity_open(&state->identity, sample, assembly, panel, frequency)) {
        state_release(state);
        return false;
    }
    state->seen = duckdb_malloc(words * sizeof(*state->seen));
    if (state->seen == NULL) goto oom;
    memset(state->seen, 0, words * sizeof(*state->seen));
    state->filter = *filter;
    state->site_count = (size_t)site_count;
    state->word_count = words;
    state->initialized = 1u;
    return true;
oom:
    state_release(state);
    return false;
}

static bool filter_same(const duckhts_somalier_contamination_settings_t *left,
                        const duckhts_somalier_contamination_settings_t *right) {
    return left->min_depth == right->min_depth &&
        left->max_depth == right->max_depth && left->max_sites == right->max_sites &&
        left->hom_minor_rate == right->hom_minor_rate &&
        left->hom_tail_alpha == right->hom_tail_alpha;
}

static bool state_matches(const contamination_state_t *state,
                          duckdb_string_t *sample,
                          duckdb_string_t *assembly, duckdb_string_t *panel,
                          duckdb_string_t *frequency, uint64_t site_count,
                          const duckhts_somalier_contamination_settings_t *filter) {
    return state->initialized && state->site_count == (size_t)site_count &&
        identity_matches(&state->identity, sample, assembly, panel, frequency) &&
        filter_same(&state->filter, filter);
}

static bool state_complete(const contamination_state_t *state) {
    if (state->observed_sites != state->site_count) return false;
    for (size_t i = 0u; i < state->word_count; i++) {
        uint64_t required = UINT64_MAX;
        if (i + 1u == state->word_count && state->site_count % 64u != 0u) {
            required >>= 64u - (unsigned)(state->site_count % 64u);
        }
        if (state->seen[i] != required) return false;
    }
    return true;
}

static bool counts_from_vectors(duckdb_vector a_vec, duckdb_vector b_vec,
                                duckdb_vector other_vec, idx_t row,
                                duckhts_somalier_counts_t *counts) {
    bool a = row_valid(a_vec, row);
    bool b = row_valid(b_vec, row);
    bool other = row_valid(other_vec, row);
    if (a != b || a != other) return false;
    memset(counts, 0, sizeof(*counts));
    if (!a) return true;
    uint64_t av = ((uint64_t *)duckdb_vector_get_data(a_vec))[row];
    uint64_t bv = ((uint64_t *)duckdb_vector_get_data(b_vec))[row];
    uint64_t ov = ((uint64_t *)duckdb_vector_get_data(other_vec))[row];
    if (av > UINT32_MAX || bv > UINT32_MAX || ov > UINT32_MAX) return false;
    counts->allele_a = (uint32_t)av;
    counts->allele_b = (uint32_t)bv;
    counts->other = (uint32_t)ov;
    counts->available = 1u;
    return true;
}

static bool mark_site(contamination_state_t *state, uint64_t site_index) {
    size_t word;
    uint64_t bit;
    if (site_index >= state->site_count) return false;
    word = (size_t)(site_index / 64u);
    bit = UINT64_C(1) << (site_index % 64u);
    if ((state->seen[word] & bit) != 0u) return false;
    state->seen[word] |= bit;
    state->observed_sites++;
    return true;
}

static bool site_is_new(const contamination_state_t *state, uint64_t site_index) {
    if (site_index >= state->site_count) return false;
    return (state->seen[site_index / 64u] &
            (UINT64_C(1) << (site_index % 64u))) == 0u;
}

static bool state_clone(contamination_state_t *target,
                        const contamination_state_t *source) {
    size_t words = source->word_count;
    memset(target, 0, sizeof(*target));
    if (!identity_clone(&target->identity, &source->identity)) goto oom;
    target->seen = duckdb_malloc(words * sizeof(*target->seen));
    if (target->seen == NULL) goto oom;
    memcpy(target->seen, source->seen, words * sizeof(*target->seen));
    target->filter = source->filter;
    target->charr = source->charr;
    target->site_count = source->site_count;
    target->word_count = words;
    target->observed_sites = source->observed_sites;
    target->unavailable_sites = source->unavailable_sites;
    target->initialized = 1u;
    return true;
oom:
    state_release(target);
    return false;
}

static bool states_same(const contamination_state_t *left,
                        const contamination_state_t *right) {
    return left->site_count == right->site_count &&
        identity_same(&left->identity, &right->identity) &&
        filter_same(&left->filter, &right->filter);
}

static bool states_combine(contamination_state_t *target,
                           const contamination_state_t *source) {
    duckhts_somalier_status_t status;
    if (!states_same(target, source) ||
        source->observed_sites > target->site_count - target->observed_sites ||
        source->unavailable_sites > target->site_count - target->unavailable_sites) return false;
    for (size_t word = 0u; word < target->word_count; word++) {
        if ((target->seen[word] & source->seen[word]) != 0u) return false;
    }
    if (source->charr.usable_sites >
            target->site_count - target->charr.usable_sites ||
        source->charr.usable_hom_a >
            target->site_count - target->charr.usable_hom_a ||
        source->charr.usable_hom_b >
            target->site_count - target->charr.usable_hom_b) return false;
    status = duckhts_somalier_charr_combine(&target->charr, &source->charr);
    if (status != DUCKHTS_SOMALIER_OK) return false;
    for (size_t word = 0u; word < target->word_count; word++) {
        target->seen[word] |= source->seen[word];
    }
    target->observed_sites += source->observed_sites;
    target->unavailable_sites += source->unavailable_sites;
    return true;
}

static void write_string(duckdb_vector field, idx_t row,
                         const char *value, uint32_t length) {
    duckdb_vector_assign_string_element_len(field, row, value, length);
}

static void write_u64(duckdb_vector field, idx_t row, uint64_t value) {
    ((uint64_t *)duckdb_vector_get_data(field))[row] = value;
}

static void write_double(duckdb_vector field, idx_t row, double value) {
    ((double *)duckdb_vector_get_data(field))[row] = value;
}

static duckdb_logical_type charr_result_type(void) {
    duckdb_logical_type fields[CH_OUT_COUNT];
    duckdb_logical_type result;
    for (unsigned i = 0u; i < CH_OUT_COUNT; i++) {
        duckdb_type kind = i <= CH_OUT_STATUS ? DUCKDB_TYPE_VARCHAR :
            i == CH_OUT_ESTIMATE || i == CH_OUT_HOM_MINOR_RATE ||
            i == CH_OUT_HOM_TAIL_ALPHA ? DUCKDB_TYPE_DOUBLE : DUCKDB_TYPE_UBIGINT;
        fields[i] = duckdb_create_logical_type(kind);
    }
    result = duckdb_create_struct_type(fields, charr_output_names, CH_OUT_COUNT);
    for (unsigned i = 0u; i < CH_OUT_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    return result;
}

static void charr_update(duckdb_function_info info, duckdb_data_chunk input,
                         duckdb_aggregate_state *states) {
    duckdb_vector vector[CH_IN_COUNT];
    duckdb_string_t *sample, *assembly, *panel, *frequency;
    uint64_t *site_index, *site_count, *min_depth, *max_depth, *max_sites;
    double *population_b_af, *hom_minor_rate, *hom_tail_alpha;
    idx_t rows = duckdb_data_chunk_get_size(input);

    if (states == NULL) {
        duckdb_aggregate_function_set_error(info,
            "duckhts_somalier_charr: aggregate state is missing");
        return;
    }
    for (unsigned i = 0u; i < CH_IN_COUNT; i++) {
        vector[i] = duckdb_data_chunk_get_vector(input, i);
    }
    sample = duckdb_vector_get_data(vector[CH_IN_SAMPLE]);
    assembly = duckdb_vector_get_data(vector[CH_IN_ASSEMBLY]);
    panel = duckdb_vector_get_data(vector[CH_IN_PANEL]);
    frequency = duckdb_vector_get_data(vector[CH_IN_FREQUENCY]);
    site_index = duckdb_vector_get_data(vector[CH_IN_SITE_INDEX]);
    site_count = duckdb_vector_get_data(vector[CH_IN_SITE_COUNT]);
    population_b_af = duckdb_vector_get_data(vector[CH_IN_POPULATION_B_AF]);
    min_depth = duckdb_vector_get_data(vector[CH_IN_MIN_DEPTH]);
    max_depth = duckdb_vector_get_data(vector[CH_IN_MAX_DEPTH]);
    hom_minor_rate = duckdb_vector_get_data(vector[CH_IN_HOM_MINOR_RATE]);
    hom_tail_alpha = duckdb_vector_get_data(vector[CH_IN_HOM_TAIL_ALPHA]);
    max_sites = duckdb_vector_get_data(vector[CH_IN_MAX_SITES]);

    for (idx_t row = 0u; row < rows; row++) {
        contamination_state_t *state = (contamination_state_t *)states[row];
        duckhts_somalier_contamination_settings_t settings;
        duckhts_somalier_counts_t counts;
        duckhts_somalier_status_t status;
        bool counts_available;
        uint8_t usable;

        if (state == NULL) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: aggregate state is missing during update");
            return;
        }
        for (unsigned i = 0u; i < CH_IN_COUNT; i++) {
            if (i >= CH_IN_A && i <= CH_IN_OTHER) continue;
            if (!row_valid(vector[i], row)) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_charr: identity, ordinal, frequency, and settings cannot be NULL");
                return;
            }
        }
        if (duckdb_string_t_length(sample[row]) > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES ||
            duckdb_string_t_length(assembly[row]) > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: sample_id and assembly must be at most 1024 bytes");
            return;
        }
        counts_available = row_valid(vector[CH_IN_A], row);
        if (!counts_from_vectors(vector[CH_IN_A], vector[CH_IN_B],
                                 vector[CH_IN_OTHER], row, &counts)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: counts must be all measured or all unavailable and fit UINT32");
            return;
        }
        if (max_sites[row] == 0u || max_sites[row] > CONTAM_SQL_MAX_SITES ||
            max_sites[row] > SIZE_MAX ||
            max_depth[row] > DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: max_sites must be 1..100000000 and max_depth no larger than 1000000");
            return;
        }
        duckhts_somalier_charr_settings_default(&settings);
        settings.min_depth = min_depth[row];
        settings.max_depth = max_depth[row];
        settings.max_sites = (size_t)max_sites[row];
        settings.hom_minor_rate = hom_minor_rate[row];
        settings.hom_tail_alpha = hom_tail_alpha[row];
        if (!isfinite(population_b_af[row]) || population_b_af[row] < 0.0 ||
            population_b_af[row] > 1.0) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: population B frequency must be finite and within [0, 1]");
            return;
        }
        status = duckhts_somalier_contamination_usable(&counts, &settings, &usable);
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_aggregate_function_set_error(info,
                status == DUCKHTS_SOMALIER_LIMIT_EXCEEDED
                    ? "duckhts_somalier_charr: a depth or site setting exceeds its declared limit"
                    : "duckhts_somalier_charr: invalid contamination settings");
            return;
        }
        if (!state->initialized) {
            if (!state_open(state, &sample[row], &assembly[row], &panel[row],
                            &frequency[row], site_count[row], &settings)) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_charr: invalid identity, site limit, or bounded-state allocation failure");
                return;
            }
        } else if (!state_matches(state, &sample[row], &assembly[row],
                                  &panel[row], &frequency[row], site_count[row],
                                  &settings)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: identity, site count, or settings changed within one group");
            return;
        }
        if (!site_is_new(state, site_index[row])) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: duplicate or out-of-range site_index");
            return;
        }
        status = duckhts_somalier_charr_observe(&state->charr, &counts,
            population_b_af[row], &settings);
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_aggregate_function_set_error(info,
                status == DUCKHTS_SOMALIER_LIMIT_EXCEEDED
                    ? "duckhts_somalier_charr: a depth or site setting exceeds its declared limit"
                    : "duckhts_somalier_charr: invalid counts, frequency, or contamination settings");
            return;
        }
        (void)mark_site(state, site_index[row]);
        if (!counts_available) state->unavailable_sites++;
    }
}

static void charr_combine(duckdb_function_info info,
                          duckdb_aggregate_state *source,
                          duckdb_aggregate_state *target, idx_t count) {
    if (source == NULL || target == NULL) {
        duckdb_aggregate_function_set_error(info,
            "duckhts_somalier_charr: aggregate state is missing during combine");
        return;
    }
    for (idx_t row = 0u; row < count; row++) {
        const contamination_state_t *src = (const contamination_state_t *)source[row];
        contamination_state_t *dst = (contamination_state_t *)target[row];
        if (src == NULL || dst == NULL) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: aggregate state is missing during combine");
            return;
        }
        if (!src->initialized) continue;
        if (!dst->initialized) {
            if (!state_clone(dst, src)) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_charr: bounded combined-state allocation failed");
                return;
            }
        } else if (!states_combine(dst, src)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: duplicate ordinal, changed identity/settings, or counter overflow across parallel states");
            return;
        }
    }
}

static void charr_finalize(duckdb_function_info info,
                           duckdb_aggregate_state *source,
                           duckdb_vector result, idx_t count, idx_t offset) {
    duckdb_vector field[CH_OUT_COUNT];
    if (source == NULL) {
        duckdb_aggregate_function_set_error(info,
            "duckhts_somalier_charr: aggregate state is missing during finalize");
        return;
    }
    for (unsigned i = 0u; i < CH_OUT_COUNT; i++) {
        field[i] = duckdb_struct_vector_get_child(result, i);
    }
    duckdb_vector_ensure_validity_writable(result);
    duckdb_vector_ensure_validity_writable(field[CH_OUT_ESTIMATE]);
    for (idx_t i = 0u; i < count; i++) {
        const contamination_state_t *state = (const contamination_state_t *)source[i];
        duckhts_somalier_charr_result_t estimate;
        duckhts_somalier_status_t status;
        idx_t row = offset + i;
        if (state == NULL) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: aggregate state is missing during finalize");
            return;
        }
        if (!state->initialized) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(result), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(result), row);
        if (!state_complete(state)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: missing site_index in completed panel");
            return;
        }
        status = duckhts_somalier_charr_finish(&state->charr, &estimate);
        if (status != DUCKHTS_SOMALIER_OK && status != DUCKHTS_SOMALIER_NO_EVIDENCE) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_charr: invalid accumulated score");
            return;
        }
        write_string(field[CH_OUT_SAMPLE], row, state->identity.sample,
                     state->identity.sample_len);
        write_string(field[CH_OUT_ASSEMBLY], row, state->identity.assembly,
                     state->identity.assembly_len);
        write_string(field[CH_OUT_PANEL], row, state->identity.panel,
                     state->identity.panel_len);
        write_string(field[CH_OUT_FREQUENCY], row, state->identity.frequency,
                     state->identity.frequency_len);
        duckdb_vector_assign_string_element(field[CH_OUT_METHOD], row,
                                             CONTAM_METHOD_VERSION);
        duckdb_vector_assign_string_element(field[CH_OUT_STATUS], row,
            status == DUCKHTS_SOMALIER_NO_EVIDENCE ? "no_evidence" : "ok");
        write_u64(field[CH_OUT_SITE_COUNT], row, state->site_count);
        write_u64(field[CH_OUT_OBSERVED], row, state->observed_sites);
        write_u64(field[CH_OUT_UNAVAILABLE], row, state->unavailable_sites);
        write_u64(field[CH_OUT_USABLE], row, estimate.usable_sites);
        write_u64(field[CH_OUT_HOM_A], row, estimate.usable_hom_a);
        write_u64(field[CH_OUT_HOM_B], row, estimate.usable_hom_b);
        if (status == DUCKHTS_SOMALIER_NO_EVIDENCE) {
            duckdb_validity_set_row_invalid(
                duckdb_vector_get_validity(field[CH_OUT_ESTIMATE]), row);
        } else {
            duckdb_validity_set_row_valid(
                duckdb_vector_get_validity(field[CH_OUT_ESTIMATE]), row);
            write_double(field[CH_OUT_ESTIMATE], row, estimate.estimate);
        }
        write_u64(field[CH_OUT_MIN_DEPTH], row, state->filter.min_depth);
        write_u64(field[CH_OUT_MAX_DEPTH], row, state->filter.max_depth);
        write_double(field[CH_OUT_HOM_MINOR_RATE], row, state->filter.hom_minor_rate);
        write_double(field[CH_OUT_HOM_TAIL_ALPHA], row, state->filter.hom_tail_alpha);
        write_u64(field[CH_OUT_MAX_SITES], row, state->filter.max_sites);
    }
}

static duckdb_aggregate_function create_charr_aggregate(
    duckdb_logical_type varchar, duckdb_logical_type bigint,
    duckdb_logical_type real, duckdb_logical_type result) {
    duckdb_aggregate_function function = duckdb_create_aggregate_function();
    if (function == NULL) return NULL;
    duckdb_aggregate_function_set_name(function, "__duckhts_somalier_charr");
    for (unsigned i = 0u; i < CH_IN_COUNT; i++) {
        duckdb_logical_type type = i <= CH_IN_FREQUENCY ? varchar :
            i == CH_IN_POPULATION_B_AF || i == CH_IN_HOM_MINOR_RATE ||
            i == CH_IN_HOM_TAIL_ALPHA ? real : bigint;
        duckdb_aggregate_function_add_parameter(function, type);
    }
    duckdb_aggregate_function_set_return_type(function, result);
    duckdb_aggregate_function_set_special_handling(function);
    duckdb_aggregate_function_set_functions(function, state_size, state_init,
        charr_update, charr_combine, charr_finalize);
    duckdb_aggregate_function_set_destructor(function, state_destroy);
    return function;
}

void register_duckhts_somalier_contamination_functions(
    duckdb_connection connection) {
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type result = charr_result_type();
    duckdb_aggregate_function charr = create_charr_aggregate(
        varchar, bigint, real, result);

    if (charr != NULL) duckdb_register_aggregate_function(connection, charr);
    duckdb_destroy_aggregate_function(&charr);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&result);
}
