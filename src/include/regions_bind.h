#ifndef DUCKHTS_REGIONS_BIND_H
#define DUCKHTS_REGIONS_BIND_H

#include "region_list.h"

/* Shared DuckDB bind adapter for the typed `regions` named parameter of
 * read_bcf and read_geno. Include after duckdb_extension.h /
 * DUCKDB_EXTENSION_EXTERN. Uses only LIST/STRUCT value accessors available
 * at C_STRUCT v1.2.0. The interval representation itself stays in
 * region_list.h and knows nothing about DuckDB. */

/* STRUCT(chrom VARCHAR, start BIGINT, "end" BIGINT)[]; caller destroys. */
static inline duckdb_logical_type duckhts_regions_parameter_type(void) {
    duckdb_logical_type children[3] = {
        duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR),
        duckdb_create_logical_type(DUCKDB_TYPE_BIGINT),
        duckdb_create_logical_type(DUCKDB_TYPE_BIGINT)
    };
    const char *names[3] = {"chrom", "start", "end"};
    duckdb_logical_type item = duckdb_create_struct_type(children, names, 3);
    duckdb_logical_type list = item ? duckdb_create_list_type(item) : NULL;
    if (item) duckdb_destroy_logical_type(&item);
    for (int i = 0; i < 3; i++) duckdb_destroy_logical_type(&children[i]);
    return list;
}

/* Reads the optional `regions` parameter into a finished plan.
 * *present is 0 for an omitted/NULL parameter (plan untouched), else 1, even
 * for an empty list. Returns 0 with a message on any failure; the plan is
 * then destroyed. Zero-initialize the plan before the call. */
static inline int duckhts_regions_bind(duckdb_bind_info info, duckhts_interval_plan_t *plan,
                                       int *present, char *error, size_t error_size) {
    *present = 0;
    duckdb_value list = duckdb_bind_get_named_parameter(info, "regions");
    if (!list) return 1;
    if (duckdb_is_null_value(list)) {
        duckdb_destroy_value(&list);
        return 1;
    }
    *present = 1;
    idx_t count = duckdb_get_list_size(list);
    if (count > DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS) {
        snprintf(error, error_size,
                 "regions: %llu intervals exceed the cap of %u intervals per query",
                 (unsigned long long)count, DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS);
        goto fail;
    }
    for (idx_t i = 0; i < count; i++) {
        duckdb_value item = duckdb_get_list_child(list, i);
        if (!item || duckdb_is_null_value(item)) {
            snprintf(error, error_size, "regions: interval %llu is NULL", (unsigned long long)(i + 1));
            if (item) duckdb_destroy_value(&item);
            goto fail;
        }
        duckdb_value chrom = duckdb_get_struct_child(item, 0);
        duckdb_value start = duckdb_get_struct_child(item, 1);
        duckdb_value end = duckdb_get_struct_child(item, 2);
        int ok = chrom && start && end && !duckdb_is_null_value(chrom) &&
                 !duckdb_is_null_value(start) && !duckdb_is_null_value(end);
        char *name = NULL;
        int alloc_failed = 0;
        if (ok) {
            name = duckdb_get_varchar(chrom);
            alloc_failed = !name;
        }
        if (!ok) {
            snprintf(error, error_size,
                     "regions: interval %llu has a NULL chrom, start or end", (unsigned long long)(i + 1));
        } else if (alloc_failed) {
            snprintf(error, error_size, "regions: out of memory reading interval %llu",
                     (unsigned long long)(i + 1));
            ok = 0;
        } else {
            ok = duckhts_interval_plan_add(plan, name, duckdb_get_int64(start),
                                           duckdb_get_int64(end), (size_t)(i + 1), error, error_size);
        }
        if (name) duckdb_free(name);
        if (chrom) duckdb_destroy_value(&chrom);
        if (start) duckdb_destroy_value(&start);
        if (end) duckdb_destroy_value(&end);
        duckdb_destroy_value(&item);
        if (!ok) goto fail;
    }
    duckdb_destroy_value(&list);
    if (!duckhts_interval_plan_finish(plan, error, error_size)) {
        duckhts_interval_plan_destroy(plan);
        return 0;
    }
    return 1;
fail:
    duckdb_destroy_value(&list);
    duckhts_interval_plan_destroy(plan);
    return 0;
}

#endif
