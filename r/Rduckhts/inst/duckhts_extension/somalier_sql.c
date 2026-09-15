/* DuckDB owns prepared sketch relations and output vectors. The numerical
 * reductions borrow their list words and use the host-neutral C kernel. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckhts_somalier.h"
#include "duckdb_list.h"

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define SOMALIER_SQL_MAX_SITES UINT64_C(100000000)
#define SOMALIER_SQL_DEFAULT_MAX_SITES UINT64_C(1000000)
#define SOMALIER_METHOD_VERSION "somalier-0.3.4-duckhts-1.5.2"

enum somalier_class_field {
    CLASS_GENOTYPE = 0,
    CLASS_MIDDLING,
    CLASS_UNAVAILABLE,
    CLASS_FIELD_COUNT
};

enum somalier_sketch_field {
    SKETCH_SAMPLE_ID = 0,
    SKETCH_ASSEMBLY,
    SKETCH_PANEL_SHA256,
    SKETCH_SITE_COUNT,
    SKETCH_MIN_DEPTH,
    SKETCH_MIN_HET_BALANCE,
    SKETCH_HOM_BALANCE_CUTOFF,
    SKETCH_MIDDLING_COUNT,
    SKETCH_UNAVAILABLE_COUNT,
    SKETCH_CONTENT_DIGEST_0,
    SKETCH_CONTENT_DIGEST_1,
    SKETCH_COUNT_DIGEST_0,
    SKETCH_COUNT_DIGEST_1,
    SKETCH_HOM_A,
    SKETCH_HET,
    SKETCH_HOM_B,
    SKETCH_FIELD_COUNT
};

enum somalier_pair_field {
    PAIR_SAMPLE_A = 0,
    PAIR_SAMPLE_B,
    PAIR_ASSEMBLY,
    PAIR_PANEL_SHA256,
    PAIR_METHOD_VERSION,
    PAIR_STATUS,
    PAIR_SITE_COUNT,
    PAIR_JOINTLY_CALLED,
    PAIR_IBS0,
    PAIR_IBS2,
    PAIR_SHARED_HETS,
    PAIR_HET_AB,
    PAIR_SHARED_HOM_B,
    PAIR_HET_COUNT_A,
    PAIR_HET_COUNT_B,
    PAIR_HOM_B_COUNT_A,
    PAIR_HOM_B_COUNT_B,
    PAIR_CALLABLE_HOM_COUNT_A,
    PAIR_CALLABLE_HOM_COUNT_B,
    PAIR_MATCHING_HOM_COUNT,
    PAIR_MIDDLING_COUNT_A,
    PAIR_MIDDLING_COUNT_B,
    PAIR_UNAVAILABLE_COUNT_A,
    PAIR_UNAVAILABLE_COUNT_B,
    PAIR_RELATEDNESS,
    PAIR_INFERRED_HOM_CONCORDANCE,
    PAIR_RAW_HOM_B_CONCORDANCE,
    PAIR_P_MIDDLING_A,
    PAIR_P_MIDDLING_B,
    PAIR_ADJUSTED_CONCORDANCE,
    PAIR_FIELD_COUNT
};

enum somalier_sketch_input {
    SKETCH_IN_SAMPLE_ID = 0,
    SKETCH_IN_ASSEMBLY,
    SKETCH_IN_PANEL_SHA256,
    SKETCH_IN_SITE_INDEX,
    SKETCH_IN_SITE_COUNT,
    SKETCH_IN_ALLELE_A,
    SKETCH_IN_ALLELE_B,
    SKETCH_IN_OTHER,
    SKETCH_IN_MIN_DEPTH,
    SKETCH_IN_MIN_HET_BALANCE,
    SKETCH_IN_HOM_BALANCE_CUTOFF,
    SKETCH_IN_MAX_SITES,
    SKETCH_IN_FIELD_COUNT
};

static const char *class_names[CLASS_FIELD_COUNT] = {
    "genotype", "middling", "unavailable"
};

static const char *sketch_names[SKETCH_FIELD_COUNT] = {
    "sample_id", "assembly", "panel_sha256", "site_count", "min_depth",
    "min_het_balance", "hom_balance_cutoff", "middling_balance_count",
    "unavailable_count", "content_digest_0", "content_digest_1",
    "count_digest_0", "count_digest_1",
    "hom_a", "het", "hom_b"
};

static const char *pair_names[PAIR_FIELD_COUNT] = {
    "sample_a", "sample_b", "assembly", "panel_sha256", "method_version",
    "status", "site_count", "jointly_called", "ibs0", "ibs2",
    "shared_hets", "het_ab", "shared_hom_b", "het_count_a", "het_count_b",
    "hom_b_count_a", "hom_b_count_b", "callable_hom_count_a",
    "callable_hom_count_b", "matching_hom_count", "middling_balance_count_a",
    "middling_balance_count_b", "unavailable_count_a", "unavailable_count_b",
    "relatedness", "inferred_hom_concordance", "raw_hom_b_concordance",
    "p_middling_a", "p_middling_b", "adjusted_concordance"
};

static bool sql_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return mask == NULL || duckdb_validity_row_is_valid(mask, row);
}

static int hex_digit(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    return -1;
}

static bool parse_panel_sha256(duckdb_string_t *value, uint8_t digest[32]) {
    const char *text;
    unsigned i;
    if (duckdb_string_t_length(*value) != 64u) return false;
    text = duckdb_string_t_data(value);
    for (i = 0u; i < 32u; i++) {
        int hi = hex_digit(text[2u * i]);
        int lo = hex_digit(text[2u * i + 1u]);
        if (hi < 0 || lo < 0) return false;
        digest[i] = (uint8_t)((hi << 4) | lo);
    }
    return true;
}

static bool string_equal(duckdb_string_t *a, duckdb_string_t *b) {
    uint32_t a_len = duckdb_string_t_length(*a);
    uint32_t b_len = duckdb_string_t_length(*b);
    return a_len == b_len &&
        memcmp(duckdb_string_t_data(a), duckdb_string_t_data(b), a_len) == 0;
}

static bool string_equals_literal(duckdb_string_t *value, const char *literal) {
    size_t length = strlen(literal);
    return duckdb_string_t_length(*value) == length &&
        memcmp(duckdb_string_t_data(value), literal, length) == 0;
}

static void somalier_classify_scalar(duckdb_function_info info,
                                     duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector inputs[6];
    duckdb_vector fields[CLASS_FIELD_COUNT];
    uint64_t *a, *b, *other, *min_depth;
    double *min_het, *hom_cutoff;
    int8_t *genotypes;
    bool *middling, *unavailable;
    idx_t rows = duckdb_data_chunk_get_size(input);

    for (unsigned i = 0u; i < 6u; i++) inputs[i] = duckdb_data_chunk_get_vector(input, i);
    for (unsigned i = 0u; i < CLASS_FIELD_COUNT; i++) {
        fields[i] = duckdb_struct_vector_get_child(output, i);
    }
    a = duckdb_vector_get_data(inputs[0]);
    b = duckdb_vector_get_data(inputs[1]);
    other = duckdb_vector_get_data(inputs[2]);
    min_depth = duckdb_vector_get_data(inputs[3]);
    min_het = duckdb_vector_get_data(inputs[4]);
    hom_cutoff = duckdb_vector_get_data(inputs[5]);
    genotypes = duckdb_vector_get_data(fields[CLASS_GENOTYPE]);
    middling = duckdb_vector_get_data(fields[CLASS_MIDDLING]);
    unavailable = duckdb_vector_get_data(fields[CLASS_UNAVAILABLE]);

    for (idx_t row = 0u; row < rows; row++) {
        duckhts_somalier_relatedness_settings_t settings;
        duckhts_somalier_counts_t counts = {0};
        duckhts_somalier_classification_t classification;
        duckhts_somalier_status_t status;
        bool have_a = sql_valid(inputs[0], row);
        bool have_b = sql_valid(inputs[1], row);
        bool have_other = sql_valid(inputs[2], row);

        if (!sql_valid(inputs[3], row) || !sql_valid(inputs[4], row) ||
            !sql_valid(inputs[5], row) ||
            (have_a != have_b) || (have_a != have_other)) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_classify: counts must be all measured or all unavailable");
            return;
        }
        if (have_a && (a[row] > UINT32_MAX || b[row] > UINT32_MAX ||
                       other[row] > UINT32_MAX)) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_classify: measured counts exceed UINT32_MAX");
            return;
        }
        duckhts_somalier_relatedness_settings_default(&settings);
        settings.min_depth = min_depth[row];
        settings.min_het_balance = min_het[row];
        settings.hom_balance_cutoff = hom_cutoff[row];
        settings.max_sites = (size_t)SOMALIER_SQL_MAX_SITES;
        if (have_a) {
            counts.allele_a = (uint32_t)a[row];
            counts.allele_b = (uint32_t)b[row];
            counts.other = (uint32_t)other[row];
            counts.available = 1u;
        }
        status = duckhts_somalier_classify_relatedness_evidence(
            &counts, &settings, &classification);
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_classify: invalid classification settings");
            return;
        }
        genotypes[row] = (int8_t)classification.genotype;
        middling[row] = classification.middling_balance != 0u;
        unavailable[row] = classification.unavailable != 0u;
    }
}

static duckdb_logical_type class_result_type(void) {
    duckdb_logical_type fields[CLASS_FIELD_COUNT];
    duckdb_logical_type result;
    fields[CLASS_GENOTYPE] = duckdb_create_logical_type(DUCKDB_TYPE_TINYINT);
    fields[CLASS_MIDDLING] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[CLASS_UNAVAILABLE] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    result = duckdb_create_struct_type(fields, class_names, CLASS_FIELD_COUNT);
    for (unsigned i = 0u; i < CLASS_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    return result;
}

typedef struct somalier_sketch_state {
    duckhts_somalier_masks_t masks;
    duckhts_somalier_relatedness_settings_t settings;
    uint64_t *storage;
    uint64_t *seen;
    char *sample_id;
    char *assembly;
    char panel_sha256[65];
    uint32_t sample_id_len;
    uint32_t assembly_len;
    uint64_t observed_sites;
    uint64_t max_sites;
    uint8_t initialized;
} somalier_sketch_state_t;

typedef struct somalier_sketch_config {
    uint8_t explicit_limit;
} somalier_sketch_config_t;

static const somalier_sketch_config_t sketch_default_config = {0u};
static const somalier_sketch_config_t sketch_explicit_config = {1u};

static idx_t somalier_sketch_state_size(duckdb_function_info info) {
    (void)info;
    return (idx_t)sizeof(somalier_sketch_state_t);
}

static void somalier_sketch_state_init(duckdb_function_info info,
                                       duckdb_aggregate_state state) {
    (void)info;
    if (state != NULL) memset(state, 0, sizeof(somalier_sketch_state_t));
}

static void somalier_sketch_state_release(somalier_sketch_state_t *state) {
    if (state == NULL) return;
    duckdb_free(state->sample_id);
    duckdb_free(state->assembly);
    duckdb_free(state->storage);
    memset(state, 0, sizeof(*state));
}

static void somalier_sketch_state_destroy(duckdb_aggregate_state *states,
                                          idx_t count) {
    if (states == NULL) return;
    for (idx_t i = 0u; i < count; i++) {
        somalier_sketch_state_release((somalier_sketch_state_t *)states[i]);
    }
}

static char *copy_duckdb_string(duckdb_string_t *value, uint32_t *length_out) {
    uint32_t length = duckdb_string_t_length(*value);
    char *copy;
    if (length == 0u || length > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES) return NULL;
    copy = duckdb_malloc((size_t)length + 1u);
    if (copy == NULL) return NULL;
    memcpy(copy, duckdb_string_t_data(value), length);
    copy[length] = '\0';
    *length_out = length;
    return copy;
}

static bool string_matches_copy(duckdb_string_t *value, const char *copy,
                                uint32_t copy_len) {
    uint32_t length = duckdb_string_t_length(*value);
    return length == copy_len &&
        memcmp(duckdb_string_t_data(value), copy, length) == 0;
}

static bool somalier_sketch_state_open(
    somalier_sketch_state_t *state, duckdb_string_t *sample,
    duckdb_string_t *assembly, duckdb_string_t *panel, uint64_t site_count,
    uint64_t max_sites, const duckhts_somalier_relatedness_settings_t *settings,
    const char **error) {
    size_t word_count;
    size_t storage_words;
    uint8_t digest[32];

    if (duckdb_string_t_length(*sample) == 0u ||
        duckdb_string_t_length(*assembly) == 0u ||
        !parse_panel_sha256(panel, digest)) {
        *error = "duckhts_somalier_sketch: sample, assembly, and lowercase panel SHA-256 are required";
        return false;
    }
    if (duckdb_string_t_length(*sample) > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES ||
        duckdb_string_t_length(*assembly) > DUCKHTS_SOMALIER_MAX_IDENTITY_BYTES) {
        *error = "duckhts_somalier_sketch: sample_id and assembly must be at most 1024 bytes";
        return false;
    }
    if (site_count == 0u || site_count > max_sites || site_count > SIZE_MAX) {
        *error = "duckhts_somalier_sketch: site_count must be positive and no larger than max_sites";
        return false;
    }
    if (duckhts_somalier_mask_word_count((size_t)site_count, &word_count) !=
        DUCKHTS_SOMALIER_OK || word_count > SIZE_MAX / 4u) {
        *error = "duckhts_somalier_sketch: mask word count overflows addressable storage";
        return false;
    }
    storage_words = word_count * 4u;
    if (storage_words > SIZE_MAX / sizeof(uint64_t)) {
        *error = "duckhts_somalier_sketch: mask storage exceeds addressable memory";
        return false;
    }
    state->sample_id = copy_duckdb_string(sample, &state->sample_id_len);
    state->assembly = copy_duckdb_string(assembly, &state->assembly_len);
    state->storage = duckdb_malloc(storage_words * sizeof(*state->storage));
    if (state->sample_id == NULL || state->assembly == NULL || state->storage == NULL) {
        somalier_sketch_state_release(state);
        *error = "duckhts_somalier_sketch: could not allocate bounded sketch storage";
        return false;
    }
    memset(state->storage, 0, storage_words * sizeof(*state->storage));
    memcpy(state->panel_sha256, duckdb_string_t_data(panel), 64u);
    state->panel_sha256[64] = '\0';
    state->settings = *settings;
    state->settings.max_sites = (size_t)max_sites;
    state->max_sites = max_sites;
    state->masks.hom_a = state->storage;
    state->masks.het = state->storage + word_count;
    state->masks.hom_b = state->storage + 2u * word_count;
    state->seen = state->storage + 3u * word_count;
    state->masks.site_count = (size_t)site_count;
    state->masks.word_count = word_count;
    memcpy(state->masks.identity.panel_sha256, digest, sizeof(digest));
    state->masks.identity.min_depth = settings->min_depth;
    state->masks.identity.min_het_balance = settings->min_het_balance == 0.0
        ? 0.0 : settings->min_het_balance;
    state->masks.identity.hom_balance_cutoff = settings->hom_balance_cutoff == 0.0
        ? 0.0 : settings->hom_balance_cutoff;
    state->initialized = 1u;
    return true;
}

static bool somalier_sketch_metadata_matches(
    const somalier_sketch_state_t *state, duckdb_string_t *sample,
    duckdb_string_t *assembly, duckdb_string_t *panel, uint64_t site_count,
    uint64_t max_sites, const duckhts_somalier_relatedness_settings_t *settings) {
    return state->initialized && state->masks.site_count == (size_t)site_count &&
        state->max_sites == max_sites && settings->min_depth == state->settings.min_depth &&
        settings->min_het_balance == state->settings.min_het_balance &&
        settings->hom_balance_cutoff == state->settings.hom_balance_cutoff &&
        string_matches_copy(sample, state->sample_id, state->sample_id_len) &&
        string_matches_copy(assembly, state->assembly, state->assembly_len) &&
        duckdb_string_t_length(*panel) == 64u &&
        memcmp(duckdb_string_t_data(panel), state->panel_sha256, 64u) == 0;
}

static void somalier_sketch_update(duckdb_function_info info,
                                   duckdb_data_chunk input,
                                   duckdb_aggregate_state *states) {
    const somalier_sketch_config_t *config =
        duckdb_aggregate_function_get_extra_info(info);
    duckdb_vector vectors[SKETCH_IN_FIELD_COUNT];
    duckdb_string_t *samples, *assemblies, *panels;
    uint64_t *site_indices, *site_counts, *allele_a, *allele_b, *other;
    uint64_t *min_depth, *max_sites_data = NULL;
    double *min_het, *hom_cutoff;
    idx_t rows = duckdb_data_chunk_get_size(input);
    unsigned input_count = config != NULL && config->explicit_limit
        ? SKETCH_IN_FIELD_COUNT : SKETCH_IN_MAX_SITES;

    if (states == NULL) {
        duckdb_aggregate_function_set_error(info,
            "duckhts_somalier_sketch: aggregate state is missing");
        return;
    }
    for (unsigned i = 0u; i < input_count; i++) {
        vectors[i] = duckdb_data_chunk_get_vector(input, i);
    }
    samples = duckdb_vector_get_data(vectors[SKETCH_IN_SAMPLE_ID]);
    assemblies = duckdb_vector_get_data(vectors[SKETCH_IN_ASSEMBLY]);
    panels = duckdb_vector_get_data(vectors[SKETCH_IN_PANEL_SHA256]);
    site_indices = duckdb_vector_get_data(vectors[SKETCH_IN_SITE_INDEX]);
    site_counts = duckdb_vector_get_data(vectors[SKETCH_IN_SITE_COUNT]);
    allele_a = duckdb_vector_get_data(vectors[SKETCH_IN_ALLELE_A]);
    allele_b = duckdb_vector_get_data(vectors[SKETCH_IN_ALLELE_B]);
    other = duckdb_vector_get_data(vectors[SKETCH_IN_OTHER]);
    min_depth = duckdb_vector_get_data(vectors[SKETCH_IN_MIN_DEPTH]);
    min_het = duckdb_vector_get_data(vectors[SKETCH_IN_MIN_HET_BALANCE]);
    hom_cutoff = duckdb_vector_get_data(vectors[SKETCH_IN_HOM_BALANCE_CUTOFF]);
    if (config != NULL && config->explicit_limit) {
        max_sites_data = duckdb_vector_get_data(vectors[SKETCH_IN_MAX_SITES]);
    }

    for (idx_t row = 0u; row < rows; row++) {
        somalier_sketch_state_t *state = (somalier_sketch_state_t *)states[row];
        duckhts_somalier_relatedness_settings_t settings;
        duckhts_somalier_counts_t counts = {0};
        duckhts_somalier_classification_t classification;
        duckhts_somalier_status_t status;
        uint64_t site_index;
        uint64_t max_sites = max_sites_data == NULL
            ? SOMALIER_SQL_DEFAULT_MAX_SITES : max_sites_data[row];
        uint64_t bit;
        size_t word;
        const char *error = NULL;
        bool have_a = sql_valid(vectors[SKETCH_IN_ALLELE_A], row);
        bool have_b = sql_valid(vectors[SKETCH_IN_ALLELE_B], row);
        bool have_other = sql_valid(vectors[SKETCH_IN_OTHER], row);

        if (state == NULL) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: aggregate state is missing during update");
            return;
        }
        for (unsigned i = 0u; i < SKETCH_IN_ALLELE_A; i++) {
            if (!sql_valid(vectors[i], row)) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_sketch: identity and site fields cannot be NULL");
                return;
            }
        }
        for (unsigned i = SKETCH_IN_MIN_DEPTH; i < input_count; i++) {
            if (!sql_valid(vectors[i], row)) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_sketch: settings cannot be NULL");
                return;
            }
        }
        if ((have_a != have_b) || (have_a != have_other)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: counts must be all measured or all unavailable");
            return;
        }
        if (max_sites == 0u || max_sites > SOMALIER_SQL_MAX_SITES ||
            (have_a && (allele_a[row] > UINT32_MAX || allele_b[row] > UINT32_MAX ||
                        other[row] > UINT32_MAX))) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: invalid max_sites or count outside UINT32 range");
            return;
        }
        duckhts_somalier_relatedness_settings_default(&settings);
        settings.min_depth = min_depth[row];
        settings.min_het_balance = min_het[row];
        settings.hom_balance_cutoff = hom_cutoff[row];
        settings.max_sites = (size_t)max_sites;
        if (have_a) {
            counts.allele_a = (uint32_t)allele_a[row];
            counts.allele_b = (uint32_t)allele_b[row];
            counts.other = (uint32_t)other[row];
            counts.available = 1u;
        }
        status = duckhts_somalier_classify_relatedness_evidence(
            &counts, &settings, &classification);
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: invalid classification settings");
            return;
        }
        if (!state->initialized) {
            if (!somalier_sketch_state_open(state, &samples[row], &assemblies[row],
                &panels[row], site_counts[row], max_sites, &settings, &error)) {
                duckdb_aggregate_function_set_error(info, error);
                return;
            }
        } else if (!somalier_sketch_metadata_matches(state, &samples[row],
                   &assemblies[row], &panels[row], site_counts[row], max_sites,
                   &settings)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: identity, site count, or settings changed within one sample group");
            return;
        }
        site_index = site_indices[row];
        if (site_index >= site_counts[row]) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: site_index is outside [0, site_count)");
            return;
        }
        word = (size_t)(site_index / 64u);
        bit = UINT64_C(1) << (site_index % 64u);
        if ((state->seen[word] & bit) != 0u) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: duplicate site_index in one sample group");
            return;
        }
        status = duckhts_somalier_count_digest_observe(
            state->masks.count_digest, site_index, &counts);
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: could not update raw-count receipt");
            return;
        }
        state->seen[word] |= bit;
        if (classification.genotype == DUCKHTS_SOMALIER_HOM_A) state->masks.hom_a[word] |= bit;
        if (classification.genotype == DUCKHTS_SOMALIER_HET) state->masks.het[word] |= bit;
        if (classification.genotype == DUCKHTS_SOMALIER_HOM_B) state->masks.hom_b[word] |= bit;
        state->masks.middling_balance_count += classification.middling_balance;
        state->masks.unavailable_count += classification.unavailable;
        state->observed_sites++;
    }
}

static bool somalier_sketch_state_clone(somalier_sketch_state_t *target,
                                        const somalier_sketch_state_t *source) {
    size_t words = source->masks.word_count;
    size_t bytes = words * 4u * sizeof(uint64_t);

    memset(target, 0, sizeof(*target));
    target->sample_id = duckdb_malloc((size_t)source->sample_id_len + 1u);
    target->assembly = duckdb_malloc((size_t)source->assembly_len + 1u);
    target->storage = duckdb_malloc(bytes);
    if (target->sample_id == NULL || target->assembly == NULL || target->storage == NULL) {
        somalier_sketch_state_release(target);
        return false;
    }
    memcpy(target->sample_id, source->sample_id, (size_t)source->sample_id_len + 1u);
    memcpy(target->assembly, source->assembly, (size_t)source->assembly_len + 1u);
    memcpy(target->storage, source->storage, bytes);
    memcpy(target->panel_sha256, source->panel_sha256, sizeof(target->panel_sha256));
    target->sample_id_len = source->sample_id_len;
    target->assembly_len = source->assembly_len;
    target->observed_sites = source->observed_sites;
    target->max_sites = source->max_sites;
    target->settings = source->settings;
    target->masks = source->masks;
    target->masks.hom_a = target->storage;
    target->masks.het = target->storage + words;
    target->masks.hom_b = target->storage + 2u * words;
    target->seen = target->storage + 3u * words;
    target->initialized = 1u;
    return true;
}

static bool somalier_sketch_states_match(const somalier_sketch_state_t *left,
                                         const somalier_sketch_state_t *right) {
    return left->sample_id_len == right->sample_id_len &&
        left->assembly_len == right->assembly_len &&
        memcmp(left->sample_id, right->sample_id, left->sample_id_len) == 0 &&
        memcmp(left->assembly, right->assembly, left->assembly_len) == 0 &&
        memcmp(left->panel_sha256, right->panel_sha256, 64u) == 0 &&
        left->max_sites == right->max_sites &&
        left->masks.site_count == right->masks.site_count &&
        left->settings.min_depth == right->settings.min_depth &&
        left->settings.min_het_balance == right->settings.min_het_balance &&
        left->settings.hom_balance_cutoff == right->settings.hom_balance_cutoff;
}

static void somalier_sketch_combine(duckdb_function_info info,
                                    duckdb_aggregate_state *source,
                                    duckdb_aggregate_state *target,
                                    idx_t count) {
    if (source == NULL || target == NULL) {
        duckdb_aggregate_function_set_error(info,
            "duckhts_somalier_sketch: aggregate state is missing during combine");
        return;
    }
    for (idx_t row = 0u; row < count; row++) {
        const somalier_sketch_state_t *src =
            (const somalier_sketch_state_t *)source[row];
        somalier_sketch_state_t *dst = (somalier_sketch_state_t *)target[row];
        if (src == NULL || dst == NULL) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: aggregate state is missing during combine");
            return;
        }
        if (!src->initialized) continue;
        if (!dst->initialized) {
            if (!somalier_sketch_state_clone(dst, src)) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_sketch: could not allocate bounded combined state");
                return;
            }
            continue;
        }
        if (!somalier_sketch_states_match(dst, src)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: identity or settings changed across parallel states");
            return;
        }
        for (size_t word = 0u; word < dst->masks.word_count; word++) {
            if ((dst->seen[word] & src->seen[word]) != 0u) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_sketch: duplicate site_index across parallel states");
                return;
            }
        }
        if (src->observed_sites > dst->masks.site_count - dst->observed_sites ||
            src->masks.middling_balance_count >
                dst->masks.site_count - dst->masks.middling_balance_count ||
            src->masks.unavailable_count >
                dst->masks.site_count - dst->masks.unavailable_count) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: site or evidence counter exceeds site_count");
            return;
        }
        for (size_t word = 0u; word < dst->masks.word_count; word++) {
            dst->seen[word] |= src->seen[word];
            dst->masks.hom_a[word] |= src->masks.hom_a[word];
            dst->masks.het[word] |= src->masks.het[word];
            dst->masks.hom_b[word] |= src->masks.hom_b[word];
        }
        dst->observed_sites += src->observed_sites;
        dst->masks.middling_balance_count += src->masks.middling_balance_count;
        dst->masks.unavailable_count += src->masks.unavailable_count;
        if (duckhts_somalier_count_digest_combine(
                dst->masks.count_digest, src->masks.count_digest) !=
            DUCKHTS_SOMALIER_OK) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: could not combine raw-count receipts");
            return;
        }
    }
}

static bool somalier_sketch_complete(somalier_sketch_state_t *state) {
    size_t words = state->masks.word_count;
    if (state->observed_sites != state->masks.site_count) return false;
    for (size_t word = 0u; word < words; word++) {
        uint64_t required = UINT64_MAX;
        if (word == words - 1u && state->masks.site_count % 64u != 0u) {
            required >>= 64u - (unsigned)(state->masks.site_count % 64u);
        }
        if (state->seen[word] != required) return false;
    }
    return duckhts_somalier_seal_masks(&state->masks) == DUCKHTS_SOMALIER_OK;
}

static void somalier_sketch_finalize(duckdb_function_info info,
                                     duckdb_aggregate_state *source,
                                     duckdb_vector result,
                                     idx_t count,
                                     idx_t offset) {
    duckdb_vector fields[SKETCH_FIELD_COUNT];
    if (source == NULL) {
        duckdb_aggregate_function_set_error(info,
            "duckhts_somalier_sketch: aggregate state is missing during finalize");
        return;
    }
    for (unsigned i = 0u; i < SKETCH_FIELD_COUNT; i++) {
        fields[i] = duckdb_struct_vector_get_child(result, i);
    }
    duckdb_vector_ensure_validity_writable(result);
    for (idx_t i = 0u; i < count; i++) {
        somalier_sketch_state_t *state =
            (somalier_sketch_state_t *)source[i];
        idx_t row = offset + i;
        duckdb_list_entry entries[3];
        if (state == NULL) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: aggregate state is missing during finalize");
            return;
        }
        if (!state->initialized) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(result), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(result), row);
        if (!somalier_sketch_complete(state)) {
            duckdb_aggregate_function_set_error(info,
                "duckhts_somalier_sketch: missing site_index or corrupt completed mask");
            return;
        }
        for (unsigned mask = 0u; mask < 3u; mask++) {
            if (!duckhts_list_extend(fields[SKETCH_HOM_A + mask],
                (idx_t)state->masks.word_count, &entries[mask])) {
                duckdb_aggregate_function_set_error(info,
                    "duckhts_somalier_sketch: DuckDB could not allocate result mask words");
                return;
            }
        }
        for (unsigned mask = 0u; mask < 3u; mask++) {
            duckdb_vector child = duckdb_list_vector_get_child(fields[SKETCH_HOM_A + mask]);
            uint64_t *data = duckdb_vector_get_data(child);
            const uint64_t *words = mask == 0u ? state->masks.hom_a :
                mask == 1u ? state->masks.het : state->masks.hom_b;
            memcpy(data + entries[mask].offset, words,
                state->masks.word_count * sizeof(*words));
            ((duckdb_list_entry *)duckdb_vector_get_data(fields[SKETCH_HOM_A + mask]))[row] =
                entries[mask];
        }
        duckdb_vector_assign_string_element_len(fields[SKETCH_SAMPLE_ID], row,
            state->sample_id, state->sample_id_len);
        duckdb_vector_assign_string_element_len(fields[SKETCH_ASSEMBLY], row,
            state->assembly, state->assembly_len);
        duckdb_vector_assign_string_element_len(fields[SKETCH_PANEL_SHA256], row,
            state->panel_sha256, 64u);
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_SITE_COUNT]))[row] =
            state->masks.site_count;
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_MIN_DEPTH]))[row] =
            state->masks.identity.min_depth;
        ((double *)duckdb_vector_get_data(fields[SKETCH_MIN_HET_BALANCE]))[row] =
            state->masks.identity.min_het_balance;
        ((double *)duckdb_vector_get_data(fields[SKETCH_HOM_BALANCE_CUTOFF]))[row] =
            state->masks.identity.hom_balance_cutoff;
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_MIDDLING_COUNT]))[row] =
            state->masks.middling_balance_count;
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_UNAVAILABLE_COUNT]))[row] =
            state->masks.unavailable_count;
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_CONTENT_DIGEST_0]))[row] =
            state->masks.content_digest[0];
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_CONTENT_DIGEST_1]))[row] =
            state->masks.content_digest[1];
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_COUNT_DIGEST_0]))[row] =
            state->masks.count_digest[0];
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_COUNT_DIGEST_1]))[row] =
            state->masks.count_digest[1];
    }
}

static duckdb_logical_type sketch_type(void) {
    duckdb_logical_type fields[SKETCH_FIELD_COUNT];
    duckdb_logical_type words = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type result;
    for (unsigned i = 0u; i < SKETCH_FIELD_COUNT; i++) {
        duckdb_type kind = i <= SKETCH_PANEL_SHA256 ? DUCKDB_TYPE_VARCHAR :
            i == SKETCH_MIN_HET_BALANCE || i == SKETCH_HOM_BALANCE_CUTOFF
                ? DUCKDB_TYPE_DOUBLE : DUCKDB_TYPE_UBIGINT;
        fields[i] = i >= SKETCH_HOM_A ? duckdb_create_list_type(words)
            : duckdb_create_logical_type(kind);
    }
    result = duckdb_create_struct_type(fields, sketch_names, SKETCH_FIELD_COUNT);
    for (unsigned i = 0u; i < SKETCH_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&words);
    return result;
}

static duckdb_logical_type pair_result_type(void) {
    duckdb_logical_type fields[PAIR_FIELD_COUNT];
    duckdb_logical_type result;
    for (unsigned i = 0u; i < PAIR_FIELD_COUNT; i++) {
        duckdb_type kind = i <= PAIR_STATUS ? DUCKDB_TYPE_VARCHAR :
            i >= PAIR_RELATEDNESS ? DUCKDB_TYPE_DOUBLE : DUCKDB_TYPE_UBIGINT;
        fields[i] = duckdb_create_logical_type(kind);
    }
    result = duckdb_create_struct_type(fields, pair_names, PAIR_FIELD_COUNT);
    for (unsigned i = 0u; i < PAIR_FIELD_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    return result;
}

typedef struct somalier_sql_sketch {
    duckhts_somalier_masks_t masks;
    duckdb_string_t *sample_id;
    duckdb_string_t *assembly;
    duckdb_string_t *panel_sha256;
} somalier_sql_sketch_t;

static bool read_sketch(duckdb_vector input, idx_t row, uint64_t max_sites,
                        somalier_sql_sketch_t *sketch, const char **error) {
    duckdb_vector fields[SKETCH_FIELD_COUNT];
    uint64_t sites, words;
    for (unsigned i = 0u; i < SKETCH_FIELD_COUNT; i++) {
        fields[i] = duckdb_struct_vector_get_child(input, i);
        if (!sql_valid(fields[i], row)) {
            *error = "duckhts_somalier_relatedness: sketch members cannot be NULL";
            return false;
        }
    }
    sketch->sample_id = &((duckdb_string_t *)duckdb_vector_get_data(fields[SKETCH_SAMPLE_ID]))[row];
    sketch->assembly = &((duckdb_string_t *)duckdb_vector_get_data(fields[SKETCH_ASSEMBLY]))[row];
    sketch->panel_sha256 = &((duckdb_string_t *)duckdb_vector_get_data(fields[SKETCH_PANEL_SHA256]))[row];
    if (!duckdb_string_t_length(*sketch->sample_id) ||
        !duckdb_string_t_length(*sketch->assembly) ||
        !parse_panel_sha256(sketch->panel_sha256, sketch->masks.identity.panel_sha256)) {
        *error = "duckhts_somalier_relatedness: invalid sample, assembly, or panel digest";
        return false;
    }
    sites = ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_SITE_COUNT]))[row];
    if (sites > max_sites || sites > SIZE_MAX) {
        *error = "duckhts_somalier_relatedness: site_count exceeds max_sites";
        return false;
    }
    words = sites / 64u + (sites % 64u != 0u);
    sketch->masks.site_count = (size_t)sites;
    sketch->masks.word_count = (size_t)words;
    sketch->masks.identity.min_depth =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_MIN_DEPTH]))[row];
    sketch->masks.identity.min_het_balance =
        ((double *)duckdb_vector_get_data(fields[SKETCH_MIN_HET_BALANCE]))[row];
    sketch->masks.identity.hom_balance_cutoff =
        ((double *)duckdb_vector_get_data(fields[SKETCH_HOM_BALANCE_CUTOFF]))[row];
    sketch->masks.middling_balance_count =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_MIDDLING_COUNT]))[row];
    sketch->masks.unavailable_count =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_UNAVAILABLE_COUNT]))[row];
    sketch->masks.content_digest[0] =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_CONTENT_DIGEST_0]))[row];
    sketch->masks.content_digest[1] =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_CONTENT_DIGEST_1]))[row];
    sketch->masks.count_digest[0] =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_COUNT_DIGEST_0]))[row];
    sketch->masks.count_digest[1] =
        ((uint64_t *)duckdb_vector_get_data(fields[SKETCH_COUNT_DIGEST_1]))[row];
    for (unsigned i = SKETCH_HOM_A; i <= SKETCH_HOM_B; i++) {
        duckdb_vector child = duckdb_list_vector_get_child(fields[i]);
        duckdb_list_entry entry = ((duckdb_list_entry *)duckdb_vector_get_data(fields[i]))[row];
        idx_t child_size = duckdb_list_vector_get_size(fields[i]);
        uint64_t *data;
        if (entry.length != words || entry.offset > child_size ||
            entry.length > child_size - entry.offset) {
            *error = "duckhts_somalier_relatedness: mask word count or offset is invalid";
            return false;
        }
        for (idx_t j = 0u; j < entry.length; j++) {
            if (!sql_valid(child, entry.offset + j)) {
                *error = "duckhts_somalier_relatedness: mask words cannot be NULL";
                return false;
            }
        }
        data = duckdb_vector_get_data(child);
        if (i == SKETCH_HOM_A) sketch->masks.hom_a = words ? data + entry.offset : NULL;
        if (i == SKETCH_HET) sketch->masks.het = words ? data + entry.offset : NULL;
        if (i == SKETCH_HOM_B) sketch->masks.hom_b = words ? data + entry.offset : NULL;
    }
    if (duckhts_somalier_validate_masks(&sketch->masks) != DUCKHTS_SOMALIER_OK) {
        *error = "duckhts_somalier_relatedness: corrupt sketch masks or settings";
        return false;
    }
    return true;
}

static void pair_write_string(duckdb_vector vector, idx_t row, duckdb_string_t *value) {
    duckdb_vector_assign_string_element_len(vector, row,
        duckdb_string_t_data(value), duckdb_string_t_length(*value));
}

static void pair_write_u64(duckdb_vector vector, idx_t row, uint64_t value) {
    ((uint64_t *)duckdb_vector_get_data(vector))[row] = value;
}

static void pair_write_double(duckdb_vector vector, idx_t row, double value) {
    ((double *)duckdb_vector_get_data(vector))[row] = value;
}

static void somalier_pair_scalar(duckdb_function_info info,
                                 duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector left = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector right = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector limit = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector result_fields[PAIR_FIELD_COUNT];
    uint64_t *max_sites = duckdb_vector_get_data(limit);
    idx_t rows = duckdb_data_chunk_get_size(input);

    for (unsigned i = 0u; i < PAIR_FIELD_COUNT; i++) {
        result_fields[i] = duckdb_struct_vector_get_child(output, i);
    }
    duckdb_vector_ensure_validity_writable(result_fields[PAIR_RELATEDNESS]);
    duckdb_vector_ensure_validity_writable(result_fields[PAIR_INFERRED_HOM_CONCORDANCE]);
    duckdb_vector_ensure_validity_writable(result_fields[PAIR_RAW_HOM_B_CONCORDANCE]);
    duckdb_vector_ensure_validity_writable(result_fields[PAIR_P_MIDDLING_A]);
    duckdb_vector_ensure_validity_writable(result_fields[PAIR_P_MIDDLING_B]);
    duckdb_vector_ensure_validity_writable(result_fields[PAIR_ADJUSTED_CONCORDANCE]);

    for (idx_t row = 0u; row < rows; row++) {
        somalier_sql_sketch_t a = {0}, b = {0};
        duckhts_somalier_pair_stats_t stats;
        duckhts_somalier_status_t status;
        const char *error = NULL;
        bool no_evidence;

        if (!sql_valid(left, row) || !sql_valid(right, row) ||
            !sql_valid(limit, row) || max_sites[row] == 0u ||
            max_sites[row] > SOMALIER_SQL_MAX_SITES) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_relatedness: sketches and max_sites must be non-NULL; max_sites is 1..100000000");
            return;
        }
        if (!read_sketch(left, row, max_sites[row], &a, &error) ||
            !read_sketch(right, row, max_sites[row], &b, &error)) {
            duckdb_scalar_function_set_error(info, error);
            return;
        }
        if (string_equal(a.sample_id, b.sample_id)) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_relatedness: sample identities must differ");
            return;
        }
        if (!string_equal(a.assembly, b.assembly)) {
            duckdb_scalar_function_set_error(info,
                "duckhts_somalier_relatedness: assembly identity mismatch");
            return;
        }
        status = duckhts_somalier_pair_stats(&a.masks, &b.masks, &stats);
        if (status != DUCKHTS_SOMALIER_OK) {
            duckdb_scalar_function_set_error(info,
                status == DUCKHTS_SOMALIER_IDENTITY_MISMATCH
                    ? "duckhts_somalier_relatedness: ordered panel or classification settings mismatch"
                    : "duckhts_somalier_relatedness: invalid paired masks");
            return;
        }
        no_evidence = stats.jointly_called == 0u;
        pair_write_string(result_fields[PAIR_SAMPLE_A], row, a.sample_id);
        pair_write_string(result_fields[PAIR_SAMPLE_B], row, b.sample_id);
        pair_write_string(result_fields[PAIR_ASSEMBLY], row, a.assembly);
        pair_write_string(result_fields[PAIR_PANEL_SHA256], row, a.panel_sha256);
        duckdb_vector_assign_string_element(result_fields[PAIR_METHOD_VERSION], row,
            SOMALIER_METHOD_VERSION);
        duckdb_vector_assign_string_element(result_fields[PAIR_STATUS], row,
            no_evidence ? "no_evidence" : "ok");
        pair_write_u64(result_fields[PAIR_SITE_COUNT], row, a.masks.site_count);
        pair_write_u64(result_fields[PAIR_JOINTLY_CALLED], row, stats.jointly_called);
        pair_write_u64(result_fields[PAIR_IBS0], row, stats.ibs0);
        pair_write_u64(result_fields[PAIR_IBS2], row, stats.ibs2);
        pair_write_u64(result_fields[PAIR_SHARED_HETS], row, stats.shared_hets);
        pair_write_u64(result_fields[PAIR_HET_AB], row, stats.het_ab);
        pair_write_u64(result_fields[PAIR_SHARED_HOM_B], row, stats.shared_hom_b);
        pair_write_u64(result_fields[PAIR_HET_COUNT_A], row, stats.het_count_a);
        pair_write_u64(result_fields[PAIR_HET_COUNT_B], row, stats.het_count_b);
        pair_write_u64(result_fields[PAIR_HOM_B_COUNT_A], row, stats.hom_b_count_a);
        pair_write_u64(result_fields[PAIR_HOM_B_COUNT_B], row, stats.hom_b_count_b);
        pair_write_u64(result_fields[PAIR_CALLABLE_HOM_COUNT_A], row, stats.callable_hom_count_a);
        pair_write_u64(result_fields[PAIR_CALLABLE_HOM_COUNT_B], row, stats.callable_hom_count_b);
        pair_write_u64(result_fields[PAIR_MATCHING_HOM_COUNT], row, stats.matching_hom_count);
        pair_write_u64(result_fields[PAIR_MIDDLING_COUNT_A], row, a.masks.middling_balance_count);
        pair_write_u64(result_fields[PAIR_MIDDLING_COUNT_B], row, b.masks.middling_balance_count);
        pair_write_u64(result_fields[PAIR_UNAVAILABLE_COUNT_A], row, a.masks.unavailable_count);
        pair_write_u64(result_fields[PAIR_UNAVAILABLE_COUNT_B], row, b.masks.unavailable_count);
        if (no_evidence) {
            for (unsigned i = PAIR_RELATEDNESS; i < PAIR_FIELD_COUNT; i++) {
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(result_fields[i]), row);
            }
        } else {
            for (unsigned i = PAIR_RELATEDNESS; i < PAIR_FIELD_COUNT; i++) {
                duckdb_validity_set_row_valid(duckdb_vector_get_validity(result_fields[i]), row);
            }
            pair_write_double(result_fields[PAIR_RELATEDNESS], row, stats.relatedness);
            pair_write_double(result_fields[PAIR_INFERRED_HOM_CONCORDANCE], row,
                stats.inferred_hom_concordance);
            pair_write_double(result_fields[PAIR_RAW_HOM_B_CONCORDANCE], row,
                stats.raw_hom_b_concordance);
            pair_write_double(result_fields[PAIR_P_MIDDLING_A], row, stats.p_middling_a);
            pair_write_double(result_fields[PAIR_P_MIDDLING_B], row, stats.p_middling_b);
            pair_write_double(result_fields[PAIR_ADJUSTED_CONCORDANCE], row,
                stats.adjusted_concordance);
        }
    }
}

static bool read_reported_pair(duckdb_vector field[PAIR_FIELD_COUNT], idx_t row,
                               bool no_evidence,
                               duckhts_somalier_pair_stats_t *stats) {
    for (unsigned i = PAIR_SITE_COUNT; i < PAIR_RELATEDNESS; i++) {
        if (!sql_valid(field[i], row)) return false;
    }
    for (unsigned i = PAIR_RELATEDNESS; i < PAIR_FIELD_COUNT; i++) {
        if (sql_valid(field[i], row) == no_evidence) return false;
    }
    memset(stats, 0, sizeof(*stats));
    stats->jointly_called = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_JOINTLY_CALLED]))[row];
    stats->ibs0 = ((uint64_t *)duckdb_vector_get_data(field[PAIR_IBS0]))[row];
    stats->ibs2 = ((uint64_t *)duckdb_vector_get_data(field[PAIR_IBS2]))[row];
    stats->shared_hets = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_SHARED_HETS]))[row];
    stats->het_ab = ((uint64_t *)duckdb_vector_get_data(field[PAIR_HET_AB]))[row];
    stats->shared_hom_b = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_SHARED_HOM_B]))[row];
    stats->het_count_a = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_HET_COUNT_A]))[row];
    stats->het_count_b = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_HET_COUNT_B]))[row];
    stats->hom_b_count_a = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_HOM_B_COUNT_A]))[row];
    stats->hom_b_count_b = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_HOM_B_COUNT_B]))[row];
    stats->callable_hom_count_a = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_CALLABLE_HOM_COUNT_A]))[row];
    stats->callable_hom_count_b = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_CALLABLE_HOM_COUNT_B]))[row];
    stats->matching_hom_count = ((uint64_t *)duckdb_vector_get_data(
        field[PAIR_MATCHING_HOM_COUNT]))[row];
    if (!no_evidence) {
        stats->relatedness = ((double *)duckdb_vector_get_data(
            field[PAIR_RELATEDNESS]))[row];
        stats->inferred_hom_concordance = ((double *)duckdb_vector_get_data(
            field[PAIR_INFERRED_HOM_CONCORDANCE]))[row];
        stats->raw_hom_b_concordance = ((double *)duckdb_vector_get_data(
            field[PAIR_RAW_HOM_B_CONCORDANCE]))[row];
        stats->p_middling_a = ((double *)duckdb_vector_get_data(
            field[PAIR_P_MIDDLING_A]))[row];
        stats->p_middling_b = ((double *)duckdb_vector_get_data(
            field[PAIR_P_MIDDLING_B]))[row];
        stats->adjusted_concordance = ((double *)duckdb_vector_get_data(
            field[PAIR_ADJUSTED_CONCORDANCE]))[row];
    }
    return true;
}

static void somalier_verify_pair_scalar(duckdb_function_info info,
                                        duckdb_data_chunk input,
                                        duckdb_vector output) {
    duckdb_vector pair = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector left = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector right = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector limit = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector field[PAIR_FIELD_COUNT];
    bool *verified = duckdb_vector_get_data(output);
    uint64_t *max_sites = duckdb_vector_get_data(limit);
    idx_t rows = duckdb_data_chunk_get_size(input);
    (void)info;
    for (unsigned i = 0u; i < PAIR_FIELD_COUNT; i++) {
        field[i] = duckdb_struct_vector_get_child(pair, i);
    }
    for (idx_t row = 0u; row < rows; row++) {
        somalier_sql_sketch_t a = {0}, b = {0};
        duckhts_somalier_pair_stats_t stats;
        duckhts_somalier_pair_stats_t expected;
        duckhts_somalier_status_t reported_status;
        duckdb_string_t *sample_a, *sample_b, *assembly, *panel, *method, *status;
        const char *error = NULL;
        bool no_evidence;
        bool text_valid = true;
        verified[row] = false;
        if (!sql_valid(pair, row) || !sql_valid(left, row) ||
            !sql_valid(right, row) || !sql_valid(limit, row) ||
            max_sites[row] == 0u || max_sites[row] > SOMALIER_SQL_MAX_SITES ||
            !read_sketch(left, row, max_sites[row], &a, &error) ||
            !read_sketch(right, row, max_sites[row], &b, &error)) continue;
        for (unsigned i = PAIR_SAMPLE_A; i <= PAIR_STATUS; i++) {
            if (!sql_valid(field[i], row)) {
                text_valid = false;
                break;
            }
        }
        if (!text_valid) continue;
        sample_a = &((duckdb_string_t *)duckdb_vector_get_data(field[PAIR_SAMPLE_A]))[row];
        sample_b = &((duckdb_string_t *)duckdb_vector_get_data(field[PAIR_SAMPLE_B]))[row];
        assembly = &((duckdb_string_t *)duckdb_vector_get_data(field[PAIR_ASSEMBLY]))[row];
        panel = &((duckdb_string_t *)duckdb_vector_get_data(field[PAIR_PANEL_SHA256]))[row];
        method = &((duckdb_string_t *)duckdb_vector_get_data(field[PAIR_METHOD_VERSION]))[row];
        status = &((duckdb_string_t *)duckdb_vector_get_data(field[PAIR_STATUS]))[row];
        no_evidence = string_equals_literal(status, "no_evidence");
        if (!no_evidence && !string_equals_literal(status, "ok")) continue;
        reported_status = no_evidence
            ? DUCKHTS_SOMALIER_NO_EVIDENCE : DUCKHTS_SOMALIER_OK;
        if (!read_reported_pair(field, row, no_evidence, &stats) ||
            !string_equal(sample_a, a.sample_id) ||
            !string_equal(sample_b, b.sample_id) ||
            string_equal(sample_a, sample_b) ||
            !string_equal(a.assembly, b.assembly) ||
            !string_equal(assembly, a.assembly) ||
            !string_equal(panel, a.panel_sha256) ||
            !string_equal(panel, b.panel_sha256) ||
            !string_equals_literal(method, SOMALIER_METHOD_VERSION) ||
            ((uint64_t *)duckdb_vector_get_data(field[PAIR_SITE_COUNT]))[row] !=
                a.masks.site_count ||
            ((uint64_t *)duckdb_vector_get_data(field[PAIR_SITE_COUNT]))[row] !=
                b.masks.site_count ||
            ((uint64_t *)duckdb_vector_get_data(field[PAIR_MIDDLING_COUNT_A]))[row] !=
                a.masks.middling_balance_count ||
            ((uint64_t *)duckdb_vector_get_data(field[PAIR_MIDDLING_COUNT_B]))[row] !=
                b.masks.middling_balance_count ||
            ((uint64_t *)duckdb_vector_get_data(field[PAIR_UNAVAILABLE_COUNT_A]))[row] !=
                a.masks.unavailable_count ||
            ((uint64_t *)duckdb_vector_get_data(field[PAIR_UNAVAILABLE_COUNT_B]))[row] !=
                b.masks.unavailable_count) continue;
        if (no_evidence) {
            if (duckhts_somalier_pair_stats(&a.masks, &b.masks, &expected) !=
                DUCKHTS_SOMALIER_OK) continue;
            stats.relatedness = expected.relatedness;
            stats.inferred_hom_concordance = expected.inferred_hom_concordance;
            stats.raw_hom_b_concordance = expected.raw_hom_b_concordance;
            stats.p_middling_a = expected.p_middling_a;
            stats.p_middling_b = expected.p_middling_b;
            stats.adjusted_concordance = expected.adjusted_concordance;
        }
        verified[row] = duckhts_somalier_verify_pair_result(
            &a.masks, &b.masks, &stats, reported_status) == DUCKHTS_SOMALIER_OK;
    }
}

static duckdb_aggregate_function create_sketch_overload(
    bool explicit_limit, duckdb_logical_type varchar,
    duckdb_logical_type bigint, duckdb_logical_type real,
    duckdb_logical_type sketch) {
    duckdb_aggregate_function function = duckdb_create_aggregate_function();
    const somalier_sketch_config_t *config = explicit_limit
        ? &sketch_explicit_config : &sketch_default_config;
    if (function == NULL) return NULL;
    duckdb_aggregate_function_set_name(function, "__duckhts_somalier_sketch");
    for (unsigned i = 0u; i < 3u; i++) {
        duckdb_aggregate_function_add_parameter(function, varchar);
    }
    for (unsigned i = 0u; i < 6u; i++) {
        duckdb_aggregate_function_add_parameter(function, bigint);
    }
    duckdb_aggregate_function_add_parameter(function, real);
    duckdb_aggregate_function_add_parameter(function, real);
    if (explicit_limit) duckdb_aggregate_function_add_parameter(function, bigint);
    duckdb_aggregate_function_set_return_type(function, sketch);
    duckdb_aggregate_function_set_special_handling(function);
    duckdb_aggregate_function_set_functions(function, somalier_sketch_state_size,
        somalier_sketch_state_init, somalier_sketch_update,
        somalier_sketch_combine, somalier_sketch_finalize);
    duckdb_aggregate_function_set_destructor(function, somalier_sketch_state_destroy);
    duckdb_aggregate_function_set_extra_info(function, (void *)config, NULL);
    return function;
}

void register_duckhts_somalier_functions(duckdb_connection connection) {
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type classified = class_result_type();
    duckdb_logical_type sketch = sketch_type();
    duckdb_logical_type pair = pair_result_type();
    duckdb_scalar_function function;
    duckdb_aggregate_function_set sketch_set =
        duckdb_create_aggregate_function_set("__duckhts_somalier_sketch");
    duckdb_aggregate_function sketch_default = NULL;
    duckdb_aggregate_function sketch_explicit = NULL;

    if (sketch_set != NULL) {
        sketch_default = create_sketch_overload(false, varchar, bigint, real, sketch);
        sketch_explicit = create_sketch_overload(true, varchar, bigint, real, sketch);
        if (sketch_default != NULL && sketch_explicit != NULL &&
            duckdb_add_aggregate_function_to_set(sketch_set, sketch_default) == DuckDBSuccess &&
            duckdb_add_aggregate_function_to_set(sketch_set, sketch_explicit) == DuckDBSuccess) {
            duckdb_register_aggregate_function_set(connection, sketch_set);
        }
    }
    if (sketch_default != NULL) duckdb_destroy_aggregate_function(&sketch_default);
    if (sketch_explicit != NULL) duckdb_destroy_aggregate_function(&sketch_explicit);
    if (sketch_set != NULL) duckdb_destroy_aggregate_function_set(&sketch_set);

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "duckhts_somalier_classify");
    for (unsigned i = 0u; i < 4u; i++) duckdb_scalar_function_add_parameter(function, bigint);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_set_return_type(function, classified);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, somalier_classify_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "duckhts_somalier_relatedness");
    duckdb_scalar_function_add_parameter(function, sketch);
    duckdb_scalar_function_add_parameter(function, sketch);
    duckdb_scalar_function_add_parameter(function, bigint);
    duckdb_scalar_function_set_return_type(function, pair);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, somalier_pair_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "duckhts_somalier_verify_relatedness");
    duckdb_scalar_function_add_parameter(function, pair);
    duckdb_scalar_function_add_parameter(function, sketch);
    duckdb_scalar_function_add_parameter(function, sketch);
    duckdb_scalar_function_add_parameter(function, bigint);
    {
        duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
        duckdb_scalar_function_set_return_type(function, boolean);
        duckdb_destroy_logical_type(&boolean);
    }
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, somalier_verify_pair_scalar);
    duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&classified);
    duckdb_destroy_logical_type(&sketch);
    duckdb_destroy_logical_type(&pair);
}
