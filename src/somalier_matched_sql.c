/* Pair fitting over DuckDB-owned, panel-aligned sample and frequency profiles. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckhts_somalier.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#define MATCHED_SQL_MAX_SITES UINT64_C(100000000)
#define MATCHED_METHOD_VERSION "somalier-0.3.4-duckhts-1.5.2"

enum filter_input {
    FILTER_IN_A = 0, FILTER_IN_B, FILTER_IN_OTHER, FILTER_IN_MIN_DEPTH,
    FILTER_IN_MAX_DEPTH, FILTER_IN_HOM_MINOR_RATE, FILTER_IN_HOM_TAIL_ALPHA,
    FILTER_IN_COUNT
};

enum filter_output {
    FILTER_OUT_RECEIVER_USABLE = 0, FILTER_OUT_ANCHOR_GENOTYPE,
    FILTER_OUT_COUNT
};

enum profile_field {
    PROFILE_SAMPLE = 0, PROFILE_ASSEMBLY, PROFILE_PANEL, PROFILE_SITE_COUNT,
    PROFILE_OBSERVED, PROFILE_UNAVAILABLE, PROFILE_MIN_DEPTH,
    PROFILE_MAX_DEPTH, PROFILE_HOM_MINOR_RATE, PROFILE_HOM_TAIL_ALPHA,
    PROFILE_RECEIVER_A, PROFILE_RECEIVER_B, PROFILE_RECEIVER_USABLE,
    PROFILE_ANCHOR_GENOTYPE, PROFILE_FIELD_COUNT
};

enum frequency_field {
    FREQUENCY_ASSEMBLY = 0, FREQUENCY_PANEL, FREQUENCY_IDENTITY,
    FREQUENCY_SITE_COUNT, FREQUENCY_VALUES, FREQUENCY_FIELD_COUNT
};

enum matched_input {
    MATCHED_IN_RECEIVER = 0, MATCHED_IN_ANCHOR, MATCHED_IN_FREQUENCY,
    MATCHED_IN_ERROR_RATE, MATCHED_IN_MIN_PROBABILITY,
    MATCHED_IN_MIN_PRIOR_FREQUENCY, MATCHED_IN_ALPHA_MIN,
    MATCHED_IN_ALPHA_MAX, MATCHED_IN_GRID_STEP, MATCHED_IN_REFINE_TOLERANCE,
    MATCHED_IN_MAX_EVALUATIONS, MATCHED_IN_MAX_SITES, MATCHED_IN_COUNT
};

enum matched_settings_input {
    SETTINGS_IN_MIN_DEPTH = 0, SETTINGS_IN_MAX_DEPTH,
    SETTINGS_IN_HOM_MINOR_RATE, SETTINGS_IN_HOM_TAIL_ALPHA,
    SETTINGS_IN_ERROR_RATE, SETTINGS_IN_MIN_PROBABILITY,
    SETTINGS_IN_MIN_PRIOR_FREQUENCY, SETTINGS_IN_ALPHA_MIN,
    SETTINGS_IN_ALPHA_MAX, SETTINGS_IN_GRID_STEP,
    SETTINGS_IN_REFINE_TOLERANCE, SETTINGS_IN_MAX_EVALUATIONS,
    SETTINGS_IN_MAX_SITES, SETTINGS_IN_COUNT
};

enum matched_output {
    MATCHED_OUT_RECEIVER = 0, MATCHED_OUT_ANCHOR, MATCHED_OUT_ASSEMBLY,
    MATCHED_OUT_PANEL, MATCHED_OUT_FREQUENCY, MATCHED_OUT_METHOD,
    MATCHED_OUT_STATUS, MATCHED_OUT_SITE_COUNT, MATCHED_OUT_OBSERVED,
    MATCHED_OUT_RECEIVER_UNAVAILABLE, MATCHED_OUT_ANCHOR_UNAVAILABLE,
    MATCHED_OUT_USABLE, MATCHED_OUT_ALPHA, MATCHED_OUT_RELATIVE_LOG_LIKELIHOOD,
    MATCHED_OUT_EVALUATIONS, MATCHED_OUT_MIN_DEPTH, MATCHED_OUT_MAX_DEPTH,
    MATCHED_OUT_HOM_MINOR_RATE, MATCHED_OUT_HOM_TAIL_ALPHA,
    MATCHED_OUT_ERROR_RATE, MATCHED_OUT_MIN_PROBABILITY,
    MATCHED_OUT_MIN_PRIOR_FREQUENCY, MATCHED_OUT_ALPHA_MIN,
    MATCHED_OUT_ALPHA_MAX, MATCHED_OUT_GRID_STEP,
    MATCHED_OUT_REFINE_TOLERANCE, MATCHED_OUT_MAX_EVALUATIONS,
    MATCHED_OUT_MAX_SITES, MATCHED_OUT_COUNT
};

static const char *filter_output_names[FILTER_OUT_COUNT] = {
    "receiver_usable", "anchor_genotype"
};

static const char *profile_names[PROFILE_FIELD_COUNT] = {
    "sample_id", "assembly", "panel_sha256", "site_count", "observed_sites",
    "unavailable_sites", "min_depth", "max_depth", "hom_minor_rate",
    "hom_tail_alpha", "receiver_a", "receiver_b", "receiver_usable",
    "anchor_genotype"
};

static const char *frequency_names[FREQUENCY_FIELD_COUNT] = {
    "assembly", "panel_sha256", "frequency_sha256", "site_count",
    "population_b_af"
};

static const char *matched_output_names[MATCHED_OUT_COUNT] = {
    "receiver_id", "anchor_id", "assembly", "panel_sha256",
    "frequency_sha256", "method_version", "status", "site_count",
    "observed_sites", "receiver_unavailable_sites",
    "anchor_unavailable_sites", "usable_sites", "alpha",
    "relative_log_likelihood", "evaluations", "min_depth", "max_depth",
    "hom_minor_rate", "hom_tail_alpha", "error_rate", "min_probability",
    "min_prior_frequency", "alpha_min", "alpha_max", "grid_step",
    "refine_tolerance", "max_evaluations", "max_sites"
};

typedef struct matched_profile_view {
    duckdb_string_t *sample;
    duckdb_string_t *assembly;
    duckdb_string_t *panel;
    uint64_t site_count;
    uint64_t observed_sites;
    uint64_t unavailable_sites;
    duckhts_somalier_contamination_settings_t filter;
    const uint32_t *receiver_a;
    const uint32_t *receiver_b;
    const uint8_t *receiver_usable;
    const int8_t *anchor_genotype;
} matched_profile_view_t;

typedef struct frequency_profile_view {
    duckdb_string_t *assembly;
    duckdb_string_t *panel;
    duckdb_string_t *identity;
    uint64_t site_count;
    const double *population_b_af;
} frequency_profile_view_t;

static bool row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return validity == NULL || duckdb_validity_row_is_valid(validity, row);
}

static bool same_string(duckdb_string_t *left, duckdb_string_t *right) {
    uint32_t left_length = duckdb_string_t_length(*left);
    uint32_t right_length = duckdb_string_t_length(*right);
    return left_length == right_length &&
        memcmp(duckdb_string_t_data(left), duckdb_string_t_data(right),
               left_length) == 0;
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

static bool counts_from_vectors(duckdb_vector a_vector, duckdb_vector b_vector,
                                duckdb_vector other_vector, idx_t row,
                                duckhts_somalier_counts_t *counts) {
    bool have_a = row_valid(a_vector, row);
    bool have_b = row_valid(b_vector, row);
    bool have_other = row_valid(other_vector, row);
    uint64_t a;
    uint64_t b;
    uint64_t other;
    if (have_a != have_b || have_a != have_other) return false;
    memset(counts, 0, sizeof(*counts));
    if (!have_a) return true;
    a = ((uint64_t *)duckdb_vector_get_data(a_vector))[row];
    b = ((uint64_t *)duckdb_vector_get_data(b_vector))[row];
    other = ((uint64_t *)duckdb_vector_get_data(other_vector))[row];
    if (a > UINT32_MAX || b > UINT32_MAX || other > UINT32_MAX) return false;
    counts->allele_a = (uint32_t)a;
    counts->allele_b = (uint32_t)b;
    counts->other = (uint32_t)other;
    counts->available = 1u;
    return true;
}

static void filter_scalar(duckdb_function_info info, duckdb_data_chunk input,
                          duckdb_vector output) {
    duckdb_vector in[FILTER_IN_COUNT];
    duckdb_vector usable_vector = duckdb_struct_vector_get_child(
        output, FILTER_OUT_RECEIVER_USABLE);
    duckdb_vector genotype_vector = duckdb_struct_vector_get_child(
        output, FILTER_OUT_ANCHOR_GENOTYPE);
    uint8_t *usable = duckdb_vector_get_data(usable_vector);
    int8_t *genotype = duckdb_vector_get_data(genotype_vector);
    idx_t rows = duckdb_data_chunk_get_size(input);

    for (unsigned i = 0u; i < FILTER_IN_COUNT; i++) {
        in[i] = duckdb_data_chunk_get_vector(input, i);
    }
    for (idx_t row = 0u; row < rows; row++) {
        duckhts_somalier_contamination_settings_t settings;
        duckhts_somalier_counts_t counts;
        duckhts_somalier_status_t status;
        uint64_t max_depth;
        uint64_t max_sites = MATCHED_SQL_MAX_SITES;

        for (unsigned i = FILTER_IN_MIN_DEPTH; i < FILTER_IN_COUNT; i++) {
            if (!row_valid(in[i], row)) {
                duckdb_scalar_function_set_error(info,
                    "duckhts_somalier_matched_contamination: filter settings cannot be NULL");
                return;
            }
        }
        if (!counts_from_vectors(in[FILTER_IN_A], in[FILTER_IN_B],
                                 in[FILTER_IN_OTHER], row, &counts)) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_matched_contamination: count tuple must be all measured or all unavailable and fit UINTEGER");
            return;
        }
        max_depth = ((uint64_t *)duckdb_vector_get_data(in[FILTER_IN_MAX_DEPTH]))[row];
        duckhts_somalier_matched_anchor_settings_default(&settings);
        settings.min_depth =
            ((uint64_t *)duckdb_vector_get_data(in[FILTER_IN_MIN_DEPTH]))[row];
        settings.max_depth = max_depth;
        settings.max_sites = (size_t)max_sites;
        settings.hom_minor_rate =
            ((double *)duckdb_vector_get_data(in[FILTER_IN_HOM_MINOR_RATE]))[row];
        settings.hom_tail_alpha =
            ((double *)duckdb_vector_get_data(in[FILTER_IN_HOM_TAIL_ALPHA]))[row];
        status = duckhts_somalier_contamination_usable(
            &counts, &settings, &usable[row]);
        if (status == DUCKHTS_SOMALIER_OK) {
            duckhts_somalier_genotype_t call;
            status = duckhts_somalier_classify_contamination(
                &counts, &settings, &call);
            genotype[row] = (int8_t)call;
        }
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_scalar_function_set_error(info,
                status == DUCKHTS_SOMALIER_LIMIT_EXCEEDED
                    ? "duckhts_somalier_matched_contamination: a sample depth exceeds max_depth"
                    : "duckhts_somalier_matched_contamination: invalid contamination filter settings");
            return;
        }
    }
}

static void matched_settings_valid_scalar(duckdb_function_info info,
                                          duckdb_data_chunk input,
                                          duckdb_vector output) {
    duckdb_vector in[SETTINGS_IN_COUNT];
    bool *valid = duckdb_vector_get_data(output);
    idx_t rows = duckdb_data_chunk_get_size(input);

    for (unsigned i = 0u; i < SETTINGS_IN_COUNT; i++) {
        in[i] = duckdb_data_chunk_get_vector(input, i);
    }
    for (idx_t row = 0u; row < rows; row++) {
        duckhts_somalier_contamination_settings_t filter;
        duckhts_somalier_matched_settings_t matched;
        duckhts_somalier_counts_t unavailable = {0};
        duckhts_somalier_matched_view_t empty = {0};
        duckhts_somalier_matched_result_t result;
        duckhts_somalier_status_t status;
        uint64_t max_evaluations;
        uint64_t max_sites;
        uint8_t usable;
        bool missing = false;

        valid[row] = false;
        for (unsigned i = 0u; i < SETTINGS_IN_COUNT; i++) {
            if (!row_valid(in[i], row)) missing = true;
        }
        if (missing) continue;
        max_evaluations = ((uint64_t *)duckdb_vector_get_data(
            in[SETTINGS_IN_MAX_EVALUATIONS]))[row];
        max_sites = ((uint64_t *)duckdb_vector_get_data(
            in[SETTINGS_IN_MAX_SITES]))[row];
        if (max_evaluations > SIZE_MAX || max_sites > SIZE_MAX ||
            max_sites > MATCHED_SQL_MAX_SITES) {
            continue;
        }

        duckhts_somalier_matched_anchor_settings_default(&filter);
        filter.min_depth = ((uint64_t *)duckdb_vector_get_data(
            in[SETTINGS_IN_MIN_DEPTH]))[row];
        filter.max_depth = ((uint64_t *)duckdb_vector_get_data(
            in[SETTINGS_IN_MAX_DEPTH]))[row];
        filter.max_sites = (size_t)max_sites;
        filter.hom_minor_rate = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_HOM_MINOR_RATE]))[row];
        filter.hom_tail_alpha = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_HOM_TAIL_ALPHA]))[row];
        if (filter.max_depth > DUCKHTS_SOMALIER_MAX_BINOMIAL_DEPTH ||
            duckhts_somalier_contamination_usable(
                &unavailable, &filter, &usable) != DUCKHTS_SOMALIER_OK) {
            continue;
        }

        duckhts_somalier_matched_settings_default(&matched);
        matched.max_depth = filter.max_depth;
        matched.max_sites = (size_t)max_sites;
        matched.max_evaluations = (size_t)max_evaluations;
        matched.error_rate = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_ERROR_RATE]))[row];
        matched.min_probability = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_MIN_PROBABILITY]))[row];
        matched.min_prior_frequency = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_MIN_PRIOR_FREQUENCY]))[row];
        matched.alpha_min = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_ALPHA_MIN]))[row];
        matched.alpha_max = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_ALPHA_MAX]))[row];
        matched.grid_step = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_GRID_STEP]))[row];
        matched.refine_tolerance = ((double *)duckdb_vector_get_data(
            in[SETTINGS_IN_REFINE_TOLERANCE]))[row];
        status = duckhts_somalier_matched_anchor(&empty, &matched, &result);
        valid[row] = status == DUCKHTS_SOMALIER_NO_EVIDENCE;
    }
    (void)info;
}

static duckdb_logical_type filter_result_type(void) {
    duckdb_logical_type fields[FILTER_OUT_COUNT];
    duckdb_logical_type result;
    fields[FILTER_OUT_RECEIVER_USABLE] =
        duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    fields[FILTER_OUT_ANCHOR_GENOTYPE] =
        duckdb_create_logical_type(DUCKDB_TYPE_TINYINT);
    result = duckdb_create_struct_type(fields, filter_output_names,
                                       FILTER_OUT_COUNT);
    for (unsigned i = 0u; i < FILTER_OUT_COUNT; i++) {
        duckdb_destroy_logical_type(&fields[i]);
    }
    return result;
}

static duckdb_logical_type profile_type(void) {
    duckdb_logical_type fields[PROFILE_FIELD_COUNT];
    duckdb_logical_type uinteger =
        duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type utinyint =
        duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
    duckdb_logical_type tinyint =
        duckdb_create_logical_type(DUCKDB_TYPE_TINYINT);
    duckdb_logical_type result;
    for (unsigned i = 0u; i < PROFILE_FIELD_COUNT; i++) {
        duckdb_type kind = i <= PROFILE_PANEL ? DUCKDB_TYPE_VARCHAR :
            i == PROFILE_HOM_MINOR_RATE || i == PROFILE_HOM_TAIL_ALPHA
                ? DUCKDB_TYPE_DOUBLE : DUCKDB_TYPE_UBIGINT;
        if (i == PROFILE_RECEIVER_A || i == PROFILE_RECEIVER_B) {
            fields[i] = duckdb_create_list_type(uinteger);
        } else if (i == PROFILE_RECEIVER_USABLE) {
            fields[i] = duckdb_create_list_type(utinyint);
        } else if (i == PROFILE_ANCHOR_GENOTYPE) {
            fields[i] = duckdb_create_list_type(tinyint);
        } else {
            fields[i] = duckdb_create_logical_type(kind);
        }
    }
    result = duckdb_create_struct_type(fields, profile_names,
                                       PROFILE_FIELD_COUNT);
    for (unsigned i = 0u; i < PROFILE_FIELD_COUNT; i++) {
        duckdb_destroy_logical_type(&fields[i]);
    }
    duckdb_destroy_logical_type(&uinteger);
    duckdb_destroy_logical_type(&utinyint);
    duckdb_destroy_logical_type(&tinyint);
    return result;
}

static duckdb_logical_type frequency_type(void) {
    duckdb_logical_type fields[FREQUENCY_FIELD_COUNT];
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type result;
    for (unsigned i = 0u; i < FREQUENCY_FIELD_COUNT; i++) {
        fields[i] = i == FREQUENCY_VALUES ? duckdb_create_list_type(real) :
            duckdb_create_logical_type(i <= FREQUENCY_IDENTITY
                ? DUCKDB_TYPE_VARCHAR : DUCKDB_TYPE_UBIGINT);
    }
    result = duckdb_create_struct_type(fields, frequency_names,
                                       FREQUENCY_FIELD_COUNT);
    for (unsigned i = 0u; i < FREQUENCY_FIELD_COUNT; i++) {
        duckdb_destroy_logical_type(&fields[i]);
    }
    duckdb_destroy_logical_type(&real);
    return result;
}

static duckdb_logical_type matched_result_type(void) {
    duckdb_logical_type fields[MATCHED_OUT_COUNT];
    duckdb_logical_type result;
    for (unsigned i = 0u; i < MATCHED_OUT_COUNT; i++) {
        duckdb_type kind = i <= MATCHED_OUT_STATUS ? DUCKDB_TYPE_VARCHAR :
            i == MATCHED_OUT_ALPHA || i == MATCHED_OUT_RELATIVE_LOG_LIKELIHOOD ||
            (i >= MATCHED_OUT_HOM_MINOR_RATE && i <= MATCHED_OUT_REFINE_TOLERANCE)
                ? DUCKDB_TYPE_DOUBLE : DUCKDB_TYPE_UBIGINT;
        fields[i] = duckdb_create_logical_type(kind);
    }
    result = duckdb_create_struct_type(fields, matched_output_names,
                                       MATCHED_OUT_COUNT);
    for (unsigned i = 0u; i < MATCHED_OUT_COUNT; i++) {
        duckdb_destroy_logical_type(&fields[i]);
    }
    return result;
}

static bool list_span(duckdb_vector list, idx_t row, uint64_t expected,
                      duckdb_vector *child, idx_t *offset, const char **error) {
    duckdb_list_entry entry =
        ((duckdb_list_entry *)duckdb_vector_get_data(list))[row];
    idx_t child_size = duckdb_list_vector_get_size(list);
    if (entry.length != expected || entry.offset > child_size ||
        entry.length > child_size - entry.offset) {
        *error = "duckhts_somalier_matched_contamination: profile list length or offset is invalid";
        return false;
    }
    *child = duckdb_list_vector_get_child(list);
    *offset = entry.offset;
    return true;
}

static bool child_has_null(duckdb_vector child, idx_t offset, idx_t count) {
    for (idx_t i = 0u; i < count; i++) {
        if (!row_valid(child, offset + i)) return true;
    }
    return false;
}

static bool read_profile(duckdb_vector input, idx_t row, uint64_t max_sites,
                         matched_profile_view_t *profile, const char **error) {
    duckdb_vector field[PROFILE_FIELD_COUNT];
    duckdb_vector child;
    idx_t offset;
    for (unsigned i = 0u; i < PROFILE_FIELD_COUNT; i++) {
        field[i] = duckdb_struct_vector_get_child(input, i);
        if (!row_valid(field[i], row)) {
            *error = "duckhts_somalier_matched_contamination: profile members cannot be NULL";
            return false;
        }
    }
    profile->sample = &((duckdb_string_t *)duckdb_vector_get_data(
        field[PROFILE_SAMPLE]))[row];
    profile->assembly = &((duckdb_string_t *)duckdb_vector_get_data(
        field[PROFILE_ASSEMBLY]))[row];
    profile->panel = &((duckdb_string_t *)duckdb_vector_get_data(
        field[PROFILE_PANEL]))[row];
    profile->site_count = ((uint64_t *)duckdb_vector_get_data(
        field[PROFILE_SITE_COUNT]))[row];
    profile->observed_sites = ((uint64_t *)duckdb_vector_get_data(
        field[PROFILE_OBSERVED]))[row];
    profile->unavailable_sites = ((uint64_t *)duckdb_vector_get_data(
        field[PROFILE_UNAVAILABLE]))[row];
    if (!duckdb_string_t_length(*profile->sample) ||
        !duckdb_string_t_length(*profile->assembly) ||
        !lowercase_sha256(profile->panel) || profile->site_count == 0u ||
        profile->site_count > max_sites ||
        profile->observed_sites != profile->site_count ||
        profile->unavailable_sites > profile->site_count) {
        *error = "duckhts_somalier_matched_contamination: invalid identities, site limit, or profile counters";
        return false;
    }
    duckhts_somalier_matched_anchor_settings_default(&profile->filter);
    profile->filter.min_depth = ((uint64_t *)duckdb_vector_get_data(
        field[PROFILE_MIN_DEPTH]))[row];
    profile->filter.max_depth = ((uint64_t *)duckdb_vector_get_data(
        field[PROFILE_MAX_DEPTH]))[row];
    profile->filter.max_sites = (size_t)max_sites;
    profile->filter.hom_minor_rate = ((double *)duckdb_vector_get_data(
        field[PROFILE_HOM_MINOR_RATE]))[row];
    profile->filter.hom_tail_alpha = ((double *)duckdb_vector_get_data(
        field[PROFILE_HOM_TAIL_ALPHA]))[row];

    if (!list_span(field[PROFILE_RECEIVER_A], row, profile->site_count,
                   &child, &offset, error) ||
        child_has_null(child, offset, (idx_t)profile->site_count)) return false;
    profile->receiver_a = ((uint32_t *)duckdb_vector_get_data(child)) + offset;
    if (!list_span(field[PROFILE_RECEIVER_B], row, profile->site_count,
                   &child, &offset, error) ||
        child_has_null(child, offset, (idx_t)profile->site_count)) return false;
    profile->receiver_b = ((uint32_t *)duckdb_vector_get_data(child)) + offset;
    if (!list_span(field[PROFILE_RECEIVER_USABLE], row, profile->site_count,
                   &child, &offset, error) ||
        child_has_null(child, offset, (idx_t)profile->site_count)) return false;
    profile->receiver_usable = ((uint8_t *)duckdb_vector_get_data(child)) + offset;
    if (!list_span(field[PROFILE_ANCHOR_GENOTYPE], row, profile->site_count,
                   &child, &offset, error) ||
        child_has_null(child, offset, (idx_t)profile->site_count)) return false;
    profile->anchor_genotype = ((int8_t *)duckdb_vector_get_data(child)) + offset;
    return true;
}

static bool read_frequency(duckdb_vector input, idx_t row, uint64_t max_sites,
                           frequency_profile_view_t *frequency,
                           const char **error) {
    duckdb_vector field[FREQUENCY_FIELD_COUNT];
    duckdb_vector child;
    idx_t offset;
    for (unsigned i = 0u; i < FREQUENCY_FIELD_COUNT; i++) {
        field[i] = duckdb_struct_vector_get_child(input, i);
        if (!row_valid(field[i], row)) {
            *error = "duckhts_somalier_matched_contamination: frequency profile members cannot be NULL";
            return false;
        }
    }
    frequency->assembly = &((duckdb_string_t *)duckdb_vector_get_data(
        field[FREQUENCY_ASSEMBLY]))[row];
    frequency->panel = &((duckdb_string_t *)duckdb_vector_get_data(
        field[FREQUENCY_PANEL]))[row];
    frequency->identity = &((duckdb_string_t *)duckdb_vector_get_data(
        field[FREQUENCY_IDENTITY]))[row];
    frequency->site_count = ((uint64_t *)duckdb_vector_get_data(
        field[FREQUENCY_SITE_COUNT]))[row];
    if (!duckdb_string_t_length(*frequency->assembly) ||
        !lowercase_sha256(frequency->panel) ||
        !lowercase_sha256(frequency->identity) ||
        frequency->site_count == 0u || frequency->site_count > max_sites ||
        !list_span(field[FREQUENCY_VALUES], row, frequency->site_count,
                   &child, &offset, error) ||
        child_has_null(child, offset, (idx_t)frequency->site_count)) {
        *error = "duckhts_somalier_matched_contamination: invalid frequency profile";
        return false;
    }
    frequency->population_b_af =
        ((double *)duckdb_vector_get_data(child)) + offset;
    return true;
}

static bool filter_same(const duckhts_somalier_contamination_settings_t *left,
                        const duckhts_somalier_contamination_settings_t *right) {
    return left->min_depth == right->min_depth &&
        left->max_depth == right->max_depth &&
        left->hom_minor_rate == right->hom_minor_rate &&
        left->hom_tail_alpha == right->hom_tail_alpha;
}

static void write_string(duckdb_vector vector, idx_t row,
                         duckdb_string_t *value) {
    duckdb_vector_assign_string_element_len(vector, row,
        duckdb_string_t_data(value), duckdb_string_t_length(*value));
}

static void matched_scalar(duckdb_function_info info, duckdb_data_chunk input,
                           duckdb_vector output) {
    duckdb_vector in[MATCHED_IN_COUNT];
    duckdb_vector out[MATCHED_OUT_COUNT];
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (unsigned i = 0u; i < MATCHED_IN_COUNT; i++) {
        in[i] = duckdb_data_chunk_get_vector(input, i);
    }
    for (unsigned i = 0u; i < MATCHED_OUT_COUNT; i++) {
        out[i] = duckdb_struct_vector_get_child(output, i);
    }
    duckdb_vector_ensure_validity_writable(out[MATCHED_OUT_ALPHA]);
    duckdb_vector_ensure_validity_writable(
        out[MATCHED_OUT_RELATIVE_LOG_LIKELIHOOD]);

    for (idx_t row = 0u; row < rows; row++) {
        matched_profile_view_t receiver;
        matched_profile_view_t anchor;
        frequency_profile_view_t frequency;
        duckhts_somalier_matched_settings_t settings;
        duckhts_somalier_matched_view_t view;
        duckhts_somalier_matched_result_t result;
        duckhts_somalier_status_t status;
        const char *error = NULL;
        uint64_t max_sites;

        for (unsigned i = MATCHED_IN_RECEIVER; i < MATCHED_IN_COUNT; i++) {
            if (!row_valid(in[i], row)) {
                duckdb_scalar_function_set_error(info,
                    "duckhts_somalier_matched_contamination: profiles and search settings cannot be NULL");
                return;
            }
        }
        max_sites = ((uint64_t *)duckdb_vector_get_data(
            in[MATCHED_IN_MAX_SITES]))[row];
        if (max_sites == 0u || max_sites > MATCHED_SQL_MAX_SITES ||
            !read_profile(in[MATCHED_IN_RECEIVER], row, max_sites,
                          &receiver, &error) ||
            !read_profile(in[MATCHED_IN_ANCHOR], row, max_sites,
                          &anchor, &error) ||
            !read_frequency(in[MATCHED_IN_FREQUENCY], row, max_sites,
                            &frequency, &error)) {
            duckdb_scalar_function_set_error(info, error != NULL ? error :
                "duckhts_somalier_matched_contamination: max_sites must be 1..100000000");
            return;
        }
        if (same_string(receiver.sample, anchor.sample) ||
            !same_string(receiver.assembly, anchor.assembly) ||
            !same_string(receiver.assembly, frequency.assembly) ||
            !same_string(receiver.panel, anchor.panel) ||
            !same_string(receiver.panel, frequency.panel) ||
            receiver.site_count != anchor.site_count ||
            receiver.site_count != frequency.site_count ||
            !filter_same(&receiver.filter, &anchor.filter)) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_matched_contamination: receiver, anchor, panel, frequency, or filter identity mismatch");
            return;
        }
        duckhts_somalier_matched_settings_default(&settings);
        settings.max_depth = receiver.filter.max_depth;
        settings.max_sites = (size_t)max_sites;
        settings.error_rate = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_ERROR_RATE]))[row];
        settings.min_probability = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_MIN_PROBABILITY]))[row];
        settings.min_prior_frequency = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_MIN_PRIOR_FREQUENCY]))[row];
        settings.alpha_min = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_ALPHA_MIN]))[row];
        settings.alpha_max = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_ALPHA_MAX]))[row];
        settings.grid_step = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_GRID_STEP]))[row];
        settings.refine_tolerance = ((double *)duckdb_vector_get_data(
            in[MATCHED_IN_REFINE_TOLERANCE]))[row];
        {
            uint64_t evaluations = ((uint64_t *)duckdb_vector_get_data(
                in[MATCHED_IN_MAX_EVALUATIONS]))[row];
            if (evaluations > SIZE_MAX) {
                duckdb_scalar_function_set_error(info,
                    "duckhts_somalier_matched_contamination: max_evaluations exceeds addressable size");
                return;
            }
            settings.max_evaluations = (size_t)evaluations;
        }
        view.receiver_a = receiver.receiver_a;
        view.receiver_b = receiver.receiver_b;
        view.receiver_usable = receiver.receiver_usable;
        view.anchor_genotype = anchor.anchor_genotype;
        view.population_b_frequency = frequency.population_b_af;
        view.site_count = (size_t)receiver.site_count;
        status = duckhts_somalier_matched_anchor(&view, &settings, &result);
        if (status != DUCKHTS_SOMALIER_OK &&
            status != DUCKHTS_SOMALIER_NO_EVIDENCE) {
            duckdb_scalar_function_set_error(info,
                status == DUCKHTS_SOMALIER_LIMIT_EXCEEDED
                    ? "duckhts_somalier_matched_contamination: workspace or max_evaluations limit exceeded"
                    : "duckhts_somalier_matched_contamination: invalid profile or search settings");
            return;
        }

        write_string(out[MATCHED_OUT_RECEIVER], row, receiver.sample);
        write_string(out[MATCHED_OUT_ANCHOR], row, anchor.sample);
        write_string(out[MATCHED_OUT_ASSEMBLY], row, receiver.assembly);
        write_string(out[MATCHED_OUT_PANEL], row, receiver.panel);
        write_string(out[MATCHED_OUT_FREQUENCY], row, frequency.identity);
        duckdb_vector_assign_string_element(out[MATCHED_OUT_METHOD], row,
                                             MATCHED_METHOD_VERSION);
        duckdb_vector_assign_string_element(out[MATCHED_OUT_STATUS], row,
            status == DUCKHTS_SOMALIER_NO_EVIDENCE ? "no_evidence" : "ok");
#define WRITE_U64(field, value) \
        (((uint64_t *)duckdb_vector_get_data(out[field]))[row] = (uint64_t)(value))
#define WRITE_DOUBLE(field, value) \
        (((double *)duckdb_vector_get_data(out[field]))[row] = (double)(value))
        WRITE_U64(MATCHED_OUT_SITE_COUNT, receiver.site_count);
        WRITE_U64(MATCHED_OUT_OBSERVED, receiver.observed_sites);
        WRITE_U64(MATCHED_OUT_RECEIVER_UNAVAILABLE, receiver.unavailable_sites);
        WRITE_U64(MATCHED_OUT_ANCHOR_UNAVAILABLE, anchor.unavailable_sites);
        WRITE_U64(MATCHED_OUT_USABLE, result.usable_sites);
        if (status == DUCKHTS_SOMALIER_NO_EVIDENCE) {
            duckdb_validity_set_row_invalid(
                duckdb_vector_get_validity(out[MATCHED_OUT_ALPHA]), row);
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(
                out[MATCHED_OUT_RELATIVE_LOG_LIKELIHOOD]), row);
        } else {
            duckdb_validity_set_row_valid(
                duckdb_vector_get_validity(out[MATCHED_OUT_ALPHA]), row);
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(
                out[MATCHED_OUT_RELATIVE_LOG_LIKELIHOOD]), row);
            WRITE_DOUBLE(MATCHED_OUT_ALPHA, result.alpha);
            WRITE_DOUBLE(MATCHED_OUT_RELATIVE_LOG_LIKELIHOOD,
                         result.log_likelihood);
        }
        WRITE_U64(MATCHED_OUT_EVALUATIONS, result.evaluations);
        WRITE_U64(MATCHED_OUT_MIN_DEPTH, receiver.filter.min_depth);
        WRITE_U64(MATCHED_OUT_MAX_DEPTH, receiver.filter.max_depth);
        WRITE_DOUBLE(MATCHED_OUT_HOM_MINOR_RATE,
                     receiver.filter.hom_minor_rate);
        WRITE_DOUBLE(MATCHED_OUT_HOM_TAIL_ALPHA,
                     receiver.filter.hom_tail_alpha);
        WRITE_DOUBLE(MATCHED_OUT_ERROR_RATE, settings.error_rate);
        WRITE_DOUBLE(MATCHED_OUT_MIN_PROBABILITY, settings.min_probability);
        WRITE_DOUBLE(MATCHED_OUT_MIN_PRIOR_FREQUENCY,
                     settings.min_prior_frequency);
        WRITE_DOUBLE(MATCHED_OUT_ALPHA_MIN, settings.alpha_min);
        WRITE_DOUBLE(MATCHED_OUT_ALPHA_MAX, settings.alpha_max);
        WRITE_DOUBLE(MATCHED_OUT_GRID_STEP, settings.grid_step);
        WRITE_DOUBLE(MATCHED_OUT_REFINE_TOLERANCE,
                     settings.refine_tolerance);
        WRITE_U64(MATCHED_OUT_MAX_EVALUATIONS, settings.max_evaluations);
        WRITE_U64(MATCHED_OUT_MAX_SITES, settings.max_sites);
#undef WRITE_U64
#undef WRITE_DOUBLE
    }
}

void register_duckhts_somalier_matched_functions(duckdb_connection connection) {
    duckdb_logical_type ubigint = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type filter_result = filter_result_type();
    duckdb_logical_type profile = profile_type();
    duckdb_logical_type frequency = frequency_type();
    duckdb_logical_type matched_result = matched_result_type();
    duckdb_scalar_function function = duckdb_create_scalar_function();

    duckdb_scalar_function_set_name(function,
        "__duckhts_somalier_contamination_filter");
    for (unsigned i = 0u; i < FILTER_IN_COUNT; i++) {
        duckdb_scalar_function_add_parameter(function,
            i >= FILTER_IN_HOM_MINOR_RATE ? real : ubigint);
    }
    duckdb_scalar_function_set_return_type(function, filter_result);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, filter_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function,
        "__duckhts_somalier_matched_settings_valid");
    for (unsigned i = 0u; i < SETTINGS_IN_COUNT; i++) {
        duckdb_scalar_function_add_parameter(function,
            i < SETTINGS_IN_HOM_MINOR_RATE || i >= SETTINGS_IN_MAX_EVALUATIONS
                ? ubigint : real);
    }
    duckdb_scalar_function_set_return_type(function, boolean);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function,
        matched_settings_valid_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function,
        "__duckhts_somalier_matched_profiles");
    duckdb_scalar_function_add_parameter(function, profile);
    duckdb_scalar_function_add_parameter(function, profile);
    duckdb_scalar_function_add_parameter(function, frequency);
    for (unsigned i = MATCHED_IN_ERROR_RATE;
         i <= MATCHED_IN_REFINE_TOLERANCE; i++) {
        duckdb_scalar_function_add_parameter(function, real);
    }
    duckdb_scalar_function_add_parameter(function, ubigint);
    duckdb_scalar_function_add_parameter(function, ubigint);
    duckdb_scalar_function_set_return_type(function, matched_result);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, matched_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&ubigint);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&filter_result);
    duckdb_destroy_logical_type(&profile);
    duckdb_destroy_logical_type(&frequency);
    duckdb_destroy_logical_type(&matched_result);
}
