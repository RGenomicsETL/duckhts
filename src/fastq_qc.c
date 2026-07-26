#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "duckhts_simd_internal.h"

#define DUCKHTS_FASTQ_QC_DEFAULT_MAX_CYCLES ((uint32_t)1048576)
#define DUCKHTS_FASTQ_QC_HARD_MAX_CYCLES ((uint32_t)16777216)
#define DUCKHTS_FASTQ_QC_MAX_PHRED ((uint64_t)93)

enum {
    FASTQ_QC_READS = 0,
    FASTQ_QC_BASES,
    FASTQ_QC_Q20,
    FASTQ_QC_Q30,
    FASTQ_QC_Q40,
    FASTQ_QC_QUALITY_SUM,
    FASTQ_QC_A,
    FASTQ_QC_C,
    FASTQ_QC_G,
    FASTQ_QC_T,
    FASTQ_QC_N,
    FASTQ_QC_OTHER,
    FASTQ_QC_MAX_READ_LENGTH,
    FASTQ_QC_CYCLES,
    FASTQ_QC_FIELD_COUNT
};

enum {
    FASTQ_QC_CYCLE_INDEX = 0,
    FASTQ_QC_CYCLE_BASES,
    FASTQ_QC_CYCLE_QUALITY_SUM,
    FASTQ_QC_CYCLE_A,
    FASTQ_QC_CYCLE_C,
    FASTQ_QC_CYCLE_G,
    FASTQ_QC_CYCLE_T,
    FASTQ_QC_CYCLE_N,
    FASTQ_QC_CYCLE_OTHER,
    FASTQ_QC_CYCLE_A_QUALITY_SUM,
    FASTQ_QC_CYCLE_C_QUALITY_SUM,
    FASTQ_QC_CYCLE_G_QUALITY_SUM,
    FASTQ_QC_CYCLE_T_QUALITY_SUM,
    FASTQ_QC_CYCLE_FIELD_COUNT
};

typedef struct {
    int explicit_max_cycles;
} duckhts_fastq_qc_config_t;

typedef struct {
    uint64_t reads;
    uint64_t bases;
    uint64_t q20;
    uint64_t q30;
    uint64_t q40;
    uint64_t quality_sum;
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
    duckhts_simd_fastq_cycle_t *cycles;
    uint32_t cycle_count;
    uint32_t cycle_capacity;
    uint32_t max_cycles;
} duckhts_fastq_qc_state_t;

static const duckhts_fastq_qc_config_t duckhts_fastq_qc_default_config = {0};
static const duckhts_fastq_qc_config_t duckhts_fastq_qc_explicit_config = {1};

static int duckhts_fastq_qc_row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return !validity || duckdb_validity_row_is_valid(validity, row);
}

static int duckhts_fastq_qc_add_overflows(uint64_t left, uint64_t right) {
    return left > UINT64_MAX - right;
}

static int duckhts_fastq_qc_ensure_cycles(duckhts_fastq_qc_state_t *state,
                                          uint32_t needed) {
    duckhts_simd_fastq_cycle_t *grown;
    uint32_t new_capacity;
    size_t old_bytes;
    size_t new_bytes;

    if (needed <= state->cycle_capacity) {
        if (needed > state->cycle_count) state->cycle_count = needed;
        return 1;
    }
    if (needed > state->max_cycles) return 0;

    new_capacity = state->cycle_capacity
        ? state->cycle_capacity
        : (state->max_cycles < 128u ? state->max_cycles : 128u);
    while (new_capacity < needed) {
        if (new_capacity > state->max_cycles / 2u) {
            new_capacity = state->max_cycles;
            break;
        }
        new_capacity *= 2u;
    }
    if (new_capacity < needed) return 0;
    old_bytes = (size_t)state->cycle_capacity * sizeof(*state->cycles);
    new_bytes = (size_t)new_capacity * sizeof(*state->cycles);
    if (new_capacity != 0 &&
        new_bytes / (size_t)new_capacity != sizeof(*state->cycles)) return 0;
    grown = (duckhts_simd_fastq_cycle_t *)realloc(state->cycles, new_bytes);
    if (!grown) return 0;
    memset((unsigned char *)grown + old_bytes, 0, new_bytes - old_bytes);
    state->cycles = grown;
    state->cycle_capacity = new_capacity;
    state->cycle_count = needed;
    return 1;
}

static idx_t duckhts_fastq_qc_state_size(duckdb_function_info info) {
    (void)info;
    return (idx_t)sizeof(duckhts_fastq_qc_state_t);
}

static void duckhts_fastq_qc_state_init(duckdb_function_info info,
                                        duckdb_aggregate_state state) {
    const duckhts_fastq_qc_config_t *config =
        (const duckhts_fastq_qc_config_t *)duckdb_aggregate_function_get_extra_info(info);
    duckhts_fastq_qc_state_t *qc = (duckhts_fastq_qc_state_t *)state;
    if (!qc) return;
    memset(qc, 0, sizeof(*qc));
    if (!config || !config->explicit_max_cycles) {
        qc->max_cycles = DUCKHTS_FASTQ_QC_DEFAULT_MAX_CYCLES;
    }
}

static void duckhts_fastq_qc_state_destroy(duckdb_aggregate_state *states,
                                           idx_t count) {
    if (!states) return;
    for (idx_t i = 0; i < count; i++) {
        duckhts_fastq_qc_state_t *state = (duckhts_fastq_qc_state_t *)states[i];
        if (!state) continue;
        free(state->cycles);
        state->cycles = NULL;
        state->cycle_count = 0;
        state->cycle_capacity = 0;
    }
}

static int duckhts_fastq_qc_set_explicit_limit(duckdb_function_info info,
                                                duckhts_fastq_qc_state_t *state,
                                                int64_t requested_limit) {
    uint32_t limit;
    if (requested_limit <= 0 ||
        requested_limit > (int64_t)DUCKHTS_FASTQ_QC_HARD_MAX_CYCLES) {
        duckdb_aggregate_function_set_error(
            info,
            "duckhts_fastq_qc max_cycles must be between 1 and 16777216");
        return 0;
    }
    limit = (uint32_t)requested_limit;
    if (state->max_cycles == 0) {
        state->max_cycles = limit;
        return 1;
    }
    if (state->max_cycles != limit) {
        duckdb_aggregate_function_set_error(
            info,
            "duckhts_fastq_qc max_cycles must be constant within each aggregate group");
        return 0;
    }
    return 1;
}

static void duckhts_fastq_qc_update(duckdb_function_info info,
                                    duckdb_data_chunk input,
                                    duckdb_aggregate_state *states) {
    const duckhts_fastq_qc_config_t *config =
        (const duckhts_fastq_qc_config_t *)duckdb_aggregate_function_get_extra_info(info);
    const duckhts_simd_dispatch_table_t *simd = duckhts_simd_dispatch_snapshot();
    idx_t count = duckdb_data_chunk_get_size(input);
    duckdb_vector sequence_vector = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector quality_vector = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector max_cycles_vector =
        config && config->explicit_max_cycles ? duckdb_data_chunk_get_vector(input, 2) : NULL;
    duckdb_string_t *sequences =
        (duckdb_string_t *)duckdb_vector_get_data(sequence_vector);
    duckdb_string_t *qualities =
        (duckdb_string_t *)duckdb_vector_get_data(quality_vector);
    int64_t *max_cycles = max_cycles_vector
        ? (int64_t *)duckdb_vector_get_data(max_cycles_vector) : NULL;

    if (!states || !sequences || !qualities) {
        duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc aggregate state is missing");
        return;
    }

    for (idx_t row = 0; row < count; row++) {
        duckhts_fastq_qc_state_t *state = (duckhts_fastq_qc_state_t *)states[row];
        duckhts_simd_fastq_read_t read_stats;
        idx_t sequence_len;
        idx_t quality_len;
        const char *sequence;
        const char *quality;
        uint64_t worst_quality_sum;

        if (!state) continue;
        if (!duckhts_fastq_qc_row_valid(sequence_vector, row) ||
            !duckhts_fastq_qc_row_valid(quality_vector, row) ||
            (max_cycles_vector &&
             !duckhts_fastq_qc_row_valid(max_cycles_vector, row))) {
            continue;
        }
        if (max_cycles_vector &&
            !duckhts_fastq_qc_set_explicit_limit(info, state, max_cycles[row])) {
            return;
        }

        sequence_len = duckdb_string_t_length(sequences[row]);
        quality_len = duckdb_string_t_length(qualities[row]);
        if (sequence_len != quality_len) {
            duckdb_aggregate_function_set_error(
                info,
                "duckhts_fastq_qc requires sequence and quality strings of equal length");
            return;
        }
        if (sequence_len > (idx_t)state->max_cycles) {
            duckdb_aggregate_function_set_error(
                info,
                "duckhts_fastq_qc read length exceeds max_cycles; use the "
                "three-argument overload with an explicit larger limit");
            return;
        }
        if (sequence_len > (idx_t)SIZE_MAX) {
            duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc read is too large for this platform");
            return;
        }
        if (sequence_len > 0 && state->reads >= UINT64_MAX / DUCKHTS_FASTQ_QC_MAX_PHRED) {
            duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc counter overflow");
            return;
        }
        worst_quality_sum = (uint64_t)sequence_len * DUCKHTS_FASTQ_QC_MAX_PHRED;
        if (duckhts_fastq_qc_add_overflows(state->bases, (uint64_t)sequence_len) ||
            duckhts_fastq_qc_add_overflows(state->quality_sum, worst_quality_sum) ||
            state->reads == UINT64_MAX) {
            duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc counter overflow");
            return;
        }
        if (!duckhts_fastq_qc_ensure_cycles(state, (uint32_t)sequence_len)) {
            duckdb_aggregate_function_set_error(
                info, "duckhts_fastq_qc could not allocate bounded per-cycle state");
            return;
        }

        sequence = duckdb_string_t_data(&sequences[row]);
        quality = duckdb_string_t_data(&qualities[row]);
        duckhts_simd_fastq_qc_read_with_table(simd, sequence, quality,
                                               (size_t)sequence_len,
                                               state->cycles, &read_stats);
        if (read_stats.invalid_quality) {
            duckdb_aggregate_function_set_error(
                info,
                "duckhts_fastq_qc quality must be canonical printable Phred+33 text");
            return;
        }

        state->reads++;
        state->bases += (uint64_t)sequence_len;
        state->q20 += read_stats.q20;
        state->q30 += read_stats.q30;
        state->q40 += read_stats.q40;
        state->quality_sum += read_stats.quality_sum;
        state->a += read_stats.a;
        state->c += read_stats.c;
        state->g += read_stats.g;
        state->t += read_stats.t;
        state->n += read_stats.n;
        state->other += read_stats.other;
    }
}

static int duckhts_fastq_qc_combine_pair(duckdb_function_info info,
                                         const duckhts_fastq_qc_state_t *source,
                                         duckhts_fastq_qc_state_t *target) {
    uint64_t combined_reads;

    if (source->max_cycles != 0) {
        if (target->max_cycles == 0) target->max_cycles = source->max_cycles;
        if (target->max_cycles != source->max_cycles) {
            duckdb_aggregate_function_set_error(
                info,
                "duckhts_fastq_qc max_cycles differed while combining aggregate states");
            return 0;
        }
    }
    if (duckhts_fastq_qc_add_overflows(target->reads, source->reads) ||
        duckhts_fastq_qc_add_overflows(target->bases, source->bases) ||
        duckhts_fastq_qc_add_overflows(target->quality_sum, source->quality_sum)) {
        duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc counter overflow while combining states");
        return 0;
    }
    combined_reads = target->reads + source->reads;
    if (source->cycle_count > 0 &&
        combined_reads > UINT64_MAX / DUCKHTS_FASTQ_QC_MAX_PHRED) {
        duckdb_aggregate_function_set_error(
            info, "duckhts_fastq_qc cycle counter overflow while combining states");
        return 0;
    }
    if (!duckhts_fastq_qc_ensure_cycles(target, source->cycle_count)) {
        duckdb_aggregate_function_set_error(
            info, "duckhts_fastq_qc could not allocate bounded per-cycle state while combining");
        return 0;
    }

    target->reads = combined_reads;
    target->bases += source->bases;
    target->q20 += source->q20;
    target->q30 += source->q30;
    target->q40 += source->q40;
    target->quality_sum += source->quality_sum;
    target->a += source->a;
    target->c += source->c;
    target->g += source->g;
    target->t += source->t;
    target->n += source->n;
    target->other += source->other;

    for (uint32_t i = 0; i < source->cycle_count; i++) {
        duckhts_simd_fastq_cycle_t *dst = &target->cycles[i];
        const duckhts_simd_fastq_cycle_t *src = &source->cycles[i];
        dst->quality_sum += src->quality_sum;
        dst->a += src->a;
        dst->c += src->c;
        dst->g += src->g;
        dst->t += src->t;
        dst->n += src->n;
        dst->other += src->other;
        dst->a_quality_sum += src->a_quality_sum;
        dst->c_quality_sum += src->c_quality_sum;
        dst->g_quality_sum += src->g_quality_sum;
        dst->t_quality_sum += src->t_quality_sum;
    }
    return 1;
}

static void duckhts_fastq_qc_combine(duckdb_function_info info,
                                     duckdb_aggregate_state *source,
                                     duckdb_aggregate_state *target,
                                     idx_t count) {
    if (!source || !target) {
        duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc aggregate state is missing");
        return;
    }
    for (idx_t i = 0; i < count; i++) {
        const duckhts_fastq_qc_state_t *src =
            (const duckhts_fastq_qc_state_t *)source[i];
        duckhts_fastq_qc_state_t *dst = (duckhts_fastq_qc_state_t *)target[i];
        if (!src || !dst) continue;
        if (!duckhts_fastq_qc_combine_pair(info, src, dst)) return;
    }
}

static uint64_t duckhts_fastq_qc_cycle_bases(const duckhts_simd_fastq_cycle_t *cycle) {
    return cycle->a + cycle->c + cycle->g + cycle->t + cycle->n + cycle->other;
}

static void duckhts_fastq_qc_finalize(duckdb_function_info info,
                                      duckdb_aggregate_state *source,
                                      duckdb_vector result,
                                      idx_t count,
                                      idx_t offset) {
    duckdb_vector fields[FASTQ_QC_FIELD_COUNT];
    duckdb_vector cycle_fields[FASTQ_QC_CYCLE_FIELD_COUNT];
    uint64_t *result_u64[FASTQ_QC_CYCLES];
    uint64_t *cycle_u64[FASTQ_QC_CYCLE_FIELD_COUNT];
    uint32_t *result_max_read_length;
    uint32_t *result_cycle_index;
    duckdb_list_entry *cycle_entries;
    idx_t child_start;
    idx_t total_cycles = 0;

    if (!source) {
        duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc aggregate state is missing");
        return;
    }
    for (int field = 0; field < FASTQ_QC_FIELD_COUNT; field++) {
        fields[field] = duckdb_struct_vector_get_child(result, (idx_t)field);
    }
    child_start = duckdb_list_vector_get_size(fields[FASTQ_QC_CYCLES]);
    for (idx_t i = 0; i < count; i++) {
        const duckhts_fastq_qc_state_t *state =
            (const duckhts_fastq_qc_state_t *)source[i];
        idx_t n = state ? (idx_t)state->cycle_count : 0;
        if (total_cycles > (idx_t)-1 - n || child_start > (idx_t)-1 - total_cycles - n) {
            duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc result is too large");
            return;
        }
        total_cycles += n;
    }
    if (duckdb_list_vector_reserve(fields[FASTQ_QC_CYCLES], child_start + total_cycles) == DuckDBError ||
        duckdb_list_vector_set_size(fields[FASTQ_QC_CYCLES], child_start + total_cycles) == DuckDBError) {
        duckdb_aggregate_function_set_error(info, "duckhts_fastq_qc could not allocate result cycles");
        return;
    }

    cycle_entries = (duckdb_list_entry *)duckdb_vector_get_data(fields[FASTQ_QC_CYCLES]);
    {
        duckdb_vector cycle_struct = duckdb_list_vector_get_child(fields[FASTQ_QC_CYCLES]);
        for (int field = 0; field < FASTQ_QC_CYCLE_FIELD_COUNT; field++) {
            cycle_fields[field] = duckdb_struct_vector_get_child(cycle_struct, (idx_t)field);
        }
    }
    for (int field = 0; field < FASTQ_QC_CYCLES; field++) {
        result_u64[field] = (uint64_t *)duckdb_vector_get_data(fields[field]);
    }
    result_max_read_length =
        (uint32_t *)duckdb_vector_get_data(fields[FASTQ_QC_MAX_READ_LENGTH]);
    result_cycle_index =
        (uint32_t *)duckdb_vector_get_data(cycle_fields[FASTQ_QC_CYCLE_INDEX]);
    cycle_u64[FASTQ_QC_CYCLE_INDEX] = NULL;
    for (int field = 1; field < FASTQ_QC_CYCLE_FIELD_COUNT; field++) {
        cycle_u64[field] = (uint64_t *)duckdb_vector_get_data(cycle_fields[field]);
    }

    total_cycles = 0;
    for (idx_t i = 0; i < count; i++) {
        const duckhts_fastq_qc_state_t empty = {0};
        const duckhts_fastq_qc_state_t *state = source[i]
            ? (const duckhts_fastq_qc_state_t *)source[i] : &empty;
        idx_t row = offset + i;
        idx_t cycle_offset = child_start + total_cycles;

        result_u64[FASTQ_QC_READS][row] = state->reads;
        result_u64[FASTQ_QC_BASES][row] = state->bases;
        result_u64[FASTQ_QC_Q20][row] = state->q20;
        result_u64[FASTQ_QC_Q30][row] = state->q30;
        result_u64[FASTQ_QC_Q40][row] = state->q40;
        result_u64[FASTQ_QC_QUALITY_SUM][row] = state->quality_sum;
        result_u64[FASTQ_QC_A][row] = state->a;
        result_u64[FASTQ_QC_C][row] = state->c;
        result_u64[FASTQ_QC_G][row] = state->g;
        result_u64[FASTQ_QC_T][row] = state->t;
        result_u64[FASTQ_QC_N][row] = state->n;
        result_u64[FASTQ_QC_OTHER][row] = state->other;
        result_max_read_length[row] = state->cycle_count;
        cycle_entries[row].offset = cycle_offset;
        cycle_entries[row].length = state->cycle_count;

        for (uint32_t cycle_index = 0; cycle_index < state->cycle_count; cycle_index++) {
            const duckhts_simd_fastq_cycle_t *cycle = &state->cycles[cycle_index];
            idx_t out = cycle_offset + (idx_t)cycle_index;
            result_cycle_index[out] = cycle_index + 1u;
            cycle_u64[FASTQ_QC_CYCLE_BASES][out] = duckhts_fastq_qc_cycle_bases(cycle);
            cycle_u64[FASTQ_QC_CYCLE_QUALITY_SUM][out] = cycle->quality_sum;
            cycle_u64[FASTQ_QC_CYCLE_A][out] = cycle->a;
            cycle_u64[FASTQ_QC_CYCLE_C][out] = cycle->c;
            cycle_u64[FASTQ_QC_CYCLE_G][out] = cycle->g;
            cycle_u64[FASTQ_QC_CYCLE_T][out] = cycle->t;
            cycle_u64[FASTQ_QC_CYCLE_N][out] = cycle->n;
            cycle_u64[FASTQ_QC_CYCLE_OTHER][out] = cycle->other;
            cycle_u64[FASTQ_QC_CYCLE_A_QUALITY_SUM][out] = cycle->a_quality_sum;
            cycle_u64[FASTQ_QC_CYCLE_C_QUALITY_SUM][out] = cycle->c_quality_sum;
            cycle_u64[FASTQ_QC_CYCLE_G_QUALITY_SUM][out] = cycle->g_quality_sum;
            cycle_u64[FASTQ_QC_CYCLE_T_QUALITY_SUM][out] = cycle->t_quality_sum;
        }
        total_cycles += state->cycle_count;
    }
}

static duckdb_logical_type duckhts_fastq_qc_result_type(void) {
    static const char *field_names[FASTQ_QC_FIELD_COUNT] = {
        "reads", "bases", "q20_bases", "q30_bases", "q40_bases",
        "quality_sum", "a_bases", "c_bases", "g_bases", "t_bases",
        "n_bases", "other_bases", "max_read_length", "cycles"
    };
    static const char *cycle_names[FASTQ_QC_CYCLE_FIELD_COUNT] = {
        "cycle", "bases", "quality_sum", "a_bases", "c_bases",
        "g_bases", "t_bases", "n_bases", "other_bases",
        "a_quality_sum", "c_quality_sum", "g_quality_sum", "t_quality_sum"
    };
    duckdb_logical_type field_types[FASTQ_QC_FIELD_COUNT];
    duckdb_logical_type cycle_types[FASTQ_QC_CYCLE_FIELD_COUNT];
    duckdb_logical_type cycle_struct;
    duckdb_logical_type cycle_list;
    duckdb_logical_type result;

    cycle_types[FASTQ_QC_CYCLE_INDEX] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    for (int i = 1; i < FASTQ_QC_CYCLE_FIELD_COUNT; i++) {
        cycle_types[i] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    }
    cycle_struct = duckdb_create_struct_type(cycle_types, cycle_names,
                                             FASTQ_QC_CYCLE_FIELD_COUNT);
    cycle_list = duckdb_create_list_type(cycle_struct);
    for (int i = 0; i < FASTQ_QC_CYCLE_FIELD_COUNT; i++) {
        duckdb_destroy_logical_type(&cycle_types[i]);
    }
    duckdb_destroy_logical_type(&cycle_struct);

    for (int i = 0; i < FASTQ_QC_FIELD_COUNT; i++) {
        field_types[i] = i == FASTQ_QC_MAX_READ_LENGTH
            ? duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER)
            : duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    }
    duckdb_destroy_logical_type(&field_types[FASTQ_QC_CYCLES]);
    field_types[FASTQ_QC_CYCLES] = cycle_list;
    result = duckdb_create_struct_type(field_types, field_names, FASTQ_QC_FIELD_COUNT);
    for (int i = 0; i < FASTQ_QC_FIELD_COUNT; i++) {
        duckdb_destroy_logical_type(&field_types[i]);
    }
    return result;
}

static duckdb_aggregate_function duckhts_fastq_qc_create_overload(
    int explicit_max_cycles,
    duckdb_logical_type varchar_type,
    duckdb_logical_type bigint_type,
    duckdb_logical_type result_type) {
    duckdb_aggregate_function function = duckdb_create_aggregate_function();
    const duckhts_fastq_qc_config_t *config = explicit_max_cycles
        ? &duckhts_fastq_qc_explicit_config : &duckhts_fastq_qc_default_config;
    if (!function) return NULL;

    duckdb_aggregate_function_set_name(function, "duckhts_fastq_qc");
    duckdb_aggregate_function_add_parameter(function, varchar_type);
    duckdb_aggregate_function_add_parameter(function, varchar_type);
    if (explicit_max_cycles) {
        duckdb_aggregate_function_add_parameter(function, bigint_type);
    }
    duckdb_aggregate_function_set_return_type(function, result_type);
    duckdb_aggregate_function_set_functions(function,
                                            duckhts_fastq_qc_state_size,
                                            duckhts_fastq_qc_state_init,
                                            duckhts_fastq_qc_update,
                                            duckhts_fastq_qc_combine,
                                            duckhts_fastq_qc_finalize);
    duckdb_aggregate_function_set_destructor(function, duckhts_fastq_qc_state_destroy);
    duckdb_aggregate_function_set_extra_info(function, (void *)config, NULL);
    return function;
}

void register_duckhts_fastq_qc_function(duckdb_connection connection) {
    duckdb_aggregate_function_set set =
        duckdb_create_aggregate_function_set("duckhts_fastq_qc");
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type result_type = duckhts_fastq_qc_result_type();
    duckdb_aggregate_function default_overload = NULL;
    duckdb_aggregate_function explicit_overload = NULL;

    if (!set || !varchar_type || !bigint_type || !result_type) goto cleanup;
    default_overload = duckhts_fastq_qc_create_overload(
        0, varchar_type, bigint_type, result_type);
    explicit_overload = duckhts_fastq_qc_create_overload(
        1, varchar_type, bigint_type, result_type);
    if (!default_overload || !explicit_overload) goto cleanup;
    if (duckdb_add_aggregate_function_to_set(set, default_overload) == DuckDBError ||
        duckdb_add_aggregate_function_to_set(set, explicit_overload) == DuckDBError) {
        goto cleanup;
    }
    duckdb_register_aggregate_function_set(connection, set);

cleanup:
    if (default_overload) duckdb_destroy_aggregate_function(&default_overload);
    if (explicit_overload) duckdb_destroy_aggregate_function(&explicit_overload);
    if (set) duckdb_destroy_aggregate_function_set(&set);
    if (varchar_type) duckdb_destroy_logical_type(&varchar_type);
    if (bigint_type) duckdb_destroy_logical_type(&bigint_type);
    if (result_type) duckdb_destroy_logical_type(&result_type);
}
