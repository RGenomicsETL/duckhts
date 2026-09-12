/* Typed GT phase preparation. Biological lane decisions belong to the native
 * reducer; DuckDB owns list storage. Two input passes require no native scratch. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckvep_phase.h"
#include "duckvep_sql.h"
#include "kernel/src/duckvep_haplotype_stream.h"

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

static bool phase_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return !validity || duckdb_validity_row_is_valid(validity, row);
}

static void phase_scalar(duckdb_function_info info, duckdb_data_chunk input,
                         duckdb_vector output) {
    duckdb_vector alleles = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector phases = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector ps = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector policies = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector allele_values = duckdb_list_vector_get_child(alleles);
    duckdb_vector phase_values = duckdb_list_vector_get_child(phases);
    duckdb_list_entry *allele_lists = duckdb_vector_get_data(alleles);
    duckdb_list_entry *phase_lists = duckdb_vector_get_data(phases);
    int32_t *values = duckdb_vector_get_data(allele_values);
    bool *flags = duckdb_vector_get_data(phase_values);
    int64_t *sets = duckdb_vector_get_data(ps);
    duckdb_string_t *policy_names = duckdb_vector_get_data(policies);
    idx_t rows = duckdb_data_chunk_get_size(input), total = 0u;

    for (idx_t row = 0u; row < rows; row++) {
        bool have_gt = phase_valid(alleles, row), have_phase = phase_valid(phases, row);
        if ((!have_gt && have_phase) ||
            (have_gt && have_phase && allele_lists[row].length != phase_lists[row].length)) {
            duckdb_scalar_function_set_error(info, "duckvep_phase_call: allele and phase lists must have equal length");
            return;
        }
        if (!have_gt) continue;
        idx_t count = allele_lists[row].length;
        if (!count || count > UINT16_MAX || count > UINT64_MAX - total) {
            duckdb_scalar_function_set_error(info, "duckvep_phase_call: ploidy must be between 1 and 65535");
            return;
        }
        total += count;
    }
    if (duckdb_list_vector_reserve(output, total) != DuckDBSuccess ||
        duckdb_list_vector_set_size(output, total) != DuckDBSuccess) {
        duckdb_scalar_function_set_error(info, "duckvep_phase_call: could not reserve output allele slots");
        return;
    }
    duckdb_vector_ensure_validity_writable(output);
    duckdb_list_entry *lists = duckdb_vector_get_data(output);
    duckdb_vector records = duckdb_list_vector_get_child(output), fields[7];
    for (idx_t i = 0u; i < 7u; i++) {
        fields[i] = duckdb_struct_vector_get_child(records, i);
        duckdb_vector_ensure_validity_writable(fields[i]);
    }
    idx_t at = 0u;
    for (idx_t row = 0u; row < rows; row++) {
        duckvep_phase_policy_t policy = DUCKVEP_PHASE_STRICT;
        if (phase_valid(policies, row)) {
            const char *name = duckdb_string_t_data(&policy_names[row]);
            uint32_t length = duckdb_string_t_length(policy_names[row]);
            if (length == 13u && !memcmp(name, "vep116_compat", 13u)) {
                policy = DUCKVEP_PHASE_VEP116_COMPAT;
            } else if (length != 6u || memcmp(name, "strict", 6u)) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: phase_policy must be 'strict' or 'vep116_compat'");
                return;
            }
        }
        lists[row] = (duckdb_list_entry){at, 0u};
        if (!phase_valid(alleles, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        duckvep_phase_summary_t summary = {0};
        idx_t count = allele_lists[row].length, base = allele_lists[row].offset;
        bool have_phase = phase_valid(phases, row);
        idx_t phase_base = have_phase ? phase_lists[row].offset : 0u;
        for (idx_t slot = 0u; slot < count; slot++) {
            bool called = phase_valid(allele_values, base + slot);
            int32_t allele = called ? values[base + slot] : -1;
            uint8_t phased = have_phase && phase_valid(phase_values, phase_base + slot)
                && flags[phase_base + slot];
            if ((called && allele < 0) ||
                duckvep_phase_observe(&summary, allele, phased) != DUCKVEP_PHASE_OK) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: called allele indices must be non-negative INTEGER values");
                return;
            }
        }
        lists[row].length = count;
        uint16_t called_before = 0u;
        for (idx_t slot = 0u; slot < count; slot++, at++) {
            bool called = phase_valid(allele_values, base + slot);
            int32_t allele = called ? values[base + slot] : -1;
            uint8_t phased = have_phase && phase_valid(phase_values, phase_base + slot)
                && flags[phase_base + slot];
            duckvep_phase_assignment_t assignment;
            if (duckvep_phase_assign(&summary, (uint16_t)(slot + 1u), called_before, allele, phased,
                                     policy, &assignment) != DUCKVEP_PHASE_OK) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: invalid decoded phase state");
                return;
            }
            if (called) called_before++;
            for (idx_t i = 0u; i < 7u; i++)
                duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[i]), at);
            ((uint16_t *)duckdb_vector_get_data(fields[0]))[at] = (uint16_t)(slot + 1u);
            ((int32_t *)duckdb_vector_get_data(fields[1]))[at] = allele;
            if (!called) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[1]), at);
            ((uint16_t *)duckdb_vector_get_data(fields[2]))[at] = assignment.lane;
            if (!assignment.lane) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[2]), at);
            ((uint16_t *)duckdb_vector_get_data(fields[3]))[at] = summary.ploidy;
            if (assignment.scope == DUCKVEP_PHASE_SET && phase_valid(ps, row)) {
                ((int64_t *)duckdb_vector_get_data(fields[4]))[at] = sets[row];
            } else {
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[4]), at);
            }
            const char *scope = assignment.scope == DUCKVEP_PHASE_SET ? "phase_set" :
                assignment.scope == DUCKVEP_PHASE_ALL_SETS ? "all_phase_sets" :
                assignment.scope == DUCKVEP_PHASE_ALLELE_SLOT ? "allele_slot" : "unresolved";
            const char *status = assignment.status == DUCKVEP_PHASE_CALLED ? "called" :
                assignment.status == DUCKVEP_PHASE_MISSING ? "missing" : "unphased";
            duckdb_vector_assign_string_element(fields[5], at, scope);
            duckdb_vector_assign_string_element(fields[6], at, status);
        }
    }
}

static void raw_gt_scalar(duckdb_function_info info, duckdb_data_chunk input,
    duckdb_vector output) {
    (void)info;
    duckdb_vector text = duckdb_data_chunk_get_vector(input, 0u);
    duckdb_vector count_vector = duckdb_data_chunk_get_vector(input, 1u);
    duckdb_string_t *strings = duckdb_vector_get_data(text);
    uint32_t *counts = duckdb_vector_get_data(count_vector);
    uint32_t *fields[7];
    for (idx_t i = 0u; i < 7u; i++)
        fields[i] = duckdb_vector_get_data(duckdb_struct_vector_get_child(output, i));
    duckdb_vector_ensure_validity_writable(output);
    for (idx_t row = 0u; row < duckdb_data_chunk_get_size(input); row++) {
        if (!phase_valid(text, row) || !phase_valid(count_vector, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        duckvep_raw_gt_t call = {0};
        duckvep_raw_gt_status_t status = duckvep_phase_parse_vep116_raw(
            (const uint8_t *)duckdb_string_t_data(&strings[row]),
            duckdb_string_t_length(strings[row]), counts[row], &call);
        fields[0][row] = (uint32_t)status;
        fields[1][row] = call.allele_index[0]; fields[2][row] = call.allele_index[1];
        fields[3][row] = call.parsed_slots; fields[4][row] = call.source_ploidy;
        fields[5][row] = call.source_has_missing; fields[6][row] = (uint32_t)call.disposition;
    }
}

static void record_order_scalar(duckdb_function_info info, duckdb_data_chunk input,
    duckdb_vector output) {
    duckdb_vector count_vector = duckdb_data_chunk_get_vector(input, 0u);
    duckdb_vector ordinal_vector = duckdb_data_chunk_get_vector(input, 1u);
    uint64_t *counts = duckdb_vector_get_data(count_vector);
    uint64_t *ordinals = duckdb_vector_get_data(ordinal_vector);
    uint64_t *ranks = duckdb_vector_get_data(output);
    duckdb_vector_ensure_validity_writable(output);
    for (idx_t row = 0u; row < duckdb_data_chunk_get_size(input); row++) {
        if (!phase_valid(count_vector, row) || !phase_valid(ordinal_vector, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        ranks[row] = duckvep_haplotype_record_order(counts[row], ordinals[row]);
        if (!ranks[row]) {
            duckdb_scalar_function_set_error(info, "duckvep_haplotypes: invalid source-buffer ordinal");
            return;
        }
    }
}

static bool register_raw_preparation(duckdb_connection connection) {
    duckdb_logical_type uinteger = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type ubigint = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type fields[7];
    for (idx_t i = 0u; i < 7u; i++) fields[i] = uinteger;
    const char *names[] = {"status", "allele0", "allele1", "parsed_slots", "source_ploidy",
        "source_has_missing", "disposition"};
    duckdb_logical_type result = duckdb_create_struct_type(fields, names, 7u);
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "_duckvep_raw_gt");
    duckdb_scalar_function_add_parameter(function, varchar);
    duckdb_scalar_function_add_parameter(function, uinteger);
    duckdb_scalar_function_set_return_type(function, result);
    duckdb_scalar_function_set_function(function, raw_gt_scalar);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_state state = duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    if (state == DuckDBSuccess) {
        function = duckdb_create_scalar_function();
        duckdb_scalar_function_set_name(function, "_duckvep_record_order");
        duckdb_scalar_function_add_parameter(function, ubigint);
        duckdb_scalar_function_add_parameter(function, ubigint);
        duckdb_scalar_function_set_return_type(function, ubigint);
        duckdb_scalar_function_set_function(function, record_order_scalar);
        duckdb_scalar_function_set_special_handling(function);
        state = duckdb_register_scalar_function(connection, function);
        duckdb_destroy_scalar_function(&function);
    }
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&ubigint);
    duckdb_destroy_logical_type(&uinteger);
    return state == DuckDBSuccess;
}

bool duckvep_register_phase_call(duckdb_connection connection) {
    if (!register_raw_preparation(connection)) return false;
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type ushort = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type alleles = duckdb_create_list_type(integer);
    duckdb_logical_type phases = duckdb_create_list_type(boolean);
    duckdb_logical_type types[] = {ushort, integer, ushort, ushort, bigint, varchar, varchar};
    const char *names[] = {"input_slot", "allele_index", "haplotype_lane", "ploidy",
        "phase_set", "phase_scope", "status"};
    duckdb_logical_type record = duckdb_create_struct_type(types, names, 7u);
    duckdb_logical_type result = duckdb_create_list_type(record);
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "_duckvep_phase_call");
    duckdb_scalar_function_add_parameter(function, alleles);
    duckdb_scalar_function_add_parameter(function, phases);
    duckdb_scalar_function_add_parameter(function, bigint);
    duckdb_scalar_function_add_parameter(function, varchar);
    duckdb_scalar_function_set_return_type(function, result);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, phase_scalar);
    duckdb_state state = duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&record);
    duckdb_destroy_logical_type(&phases);
    duckdb_destroy_logical_type(&alleles);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&ushort);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&integer);
    if (state != DuckDBSuccess) return false;
    /* A list-wise cast of SQLNULL[] can leave a constant NULL child under a
     * flat parent. The stable C callback only flattens that parent in this
     * case, so indexing the child's validity mask treats later NULLs as valid.
     * Element-wise casts materialize typed child slots, including every NULL;
     * the native reducer must never infer phase from uninitialized payload. */
    const char *sql[] = {
        "CREATE OR REPLACE MACRO duckvep_phase_call(alleles, phase_before, ",
        "phase_set := NULL, phase_policy := 'strict') AS ",
        "_duckvep_phase_call(list_transform(alleles, a -> CAST(a AS INTEGER)), ",
        "list_transform(phase_before, p -> CAST(p AS BOOLEAN)), ",
        "CAST(phase_set AS BIGINT), CAST(phase_policy AS VARCHAR))"
    };
    return duckvep_register_sql_parts(connection, sql, sizeof(sql) / sizeof(sql[0]));
}
