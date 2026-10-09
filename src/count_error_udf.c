/* __duckhts_count_error_fit(cells, min_block_sites, max_cell_bytes): the
 * count-error fit of one sample from its packed histogram (count_cells.h), as
 * one STRUCT row; a NULL histogram is a sample with no usable site and gives
 * the no_sites row. __duckhts_count_error_valid_args(...): TRUE, or the
 * error the macro's arguments would raise, so that bad arguments fail on
 * empty input too.
 *
 * The pooled fit uses every cell; the relative fit refits the same cells with
 * a parent or child as the second genome. Each block is refitted for
 * contamination with the other parameters held at the pooled values, and the
 * spread of the block estimates gives contamination_sd. Memory: one sorted
 * copy of the cells (16 bytes each) for the whole fit, and the working memory
 * of one problem at a time (count_error_model.h), all charged against
 * max_cell_bytes with the aggregate's buffers (count_cells.c); the BLOB the
 * cells arrive in is DuckDB's.
 */
#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/count_cells.h"
#include "include/count_error_model.h"

#define FIT_NAME "__duckhts_count_error_fit"
#define VALID_ARGS_NAME "__duckhts_count_error_valid_args"
#define FIT_METHOD "duckhts_count_error_fit v1"
#define FIT_FEW_SITES 1000
#define FIT_FEW_BLOCKS 3
#define FIT_ERRLEN 256

enum {
    OUT_SITES = 0,
    OUT_READS,
    OUT_MEAN_DEPTH,
    OUT_BLOCKS,
    OUT_SEQ_ERROR,
    OUT_CONTAMINATION,
    OUT_CONTAMINATION_SD,
    OUT_HOMOZYGOSITY_EXCESS,
    OUT_ALLELE_BALANCE,
    OUT_SPREAD_HOM,
    OUT_SPREAD_HET,
    OUT_ARTEFACT_WEIGHT,
    OUT_LOG_LIKELIHOOD,
    OUT_CONTAMINATION_RELATIVE,
    OUT_LOG_LIKELIHOOD_RELATIVE,
    OUT_STATUS,
    OUT_METHOD,
    OUT_COUNT
};

static const char *output_names[OUT_COUNT] = {
    "sites", "reads", "mean_depth", "blocks", "seq_error", "contamination", "contamination_sd",
    "homozygosity_excess", "allele_balance", "spread_hom", "spread_het", "artefact_weight",
    "log_likelihood", "contamination_relative", "log_likelihood_relative", "status", "method"};

static void fit_error(duckdb_function_info info, const char *message) {
    char full[320];
    snprintf(full, sizeof(full), "duckhts_count_error_fit: %s", message);
    duckdb_scalar_function_set_error(info, full);
}

static bool row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

/* The arguments of the macro, checked as the macro would use them. */
const char *duckhts_count_error_check_args(int has_block_bases, int64_t block_bases, int has_freq_bins,
                                           int64_t freq_bins, int has_max_depth, int64_t max_depth,
                                           int has_min_block_sites, int64_t min_block_sites) {
    if (!has_block_bases || block_bases < 1) return "block_bases must be at least 1";
    if (!has_freq_bins || freq_bins < 1 || freq_bins > 1024) return "freq_bins must be between 1 and 1024";
    if (!has_max_depth || max_depth < 1 || max_depth > DUCKHTS_COUNT_MAX_DEPTH) {
        return "max_depth must be between 1 and 65535";
    }
    if (!has_min_block_sites || min_block_sites < 1) return "min_block_sites must be at least 1";
    return NULL;
}

/* Cells in block, depth, alt and frequency order: the fit then sums them in
 * one order whatever order DuckDB delivered them in, and a block's cells are
 * contiguous. */
static int compare_cells(const void *a, const void *b) {
    const duckhts_count_cell_t *x = a, *y = b;
    if (x->block != y->block) return x->block < y->block ? -1 : 1;
    if (x->depth != y->depth) return x->depth < y->depth ? -1 : 1;
    if (x->alt != y->alt) return x->alt < y->alt ? -1 : 1;
    if (x->af != y->af) return x->af < y->af ? -1 : 1;
    return 0;
}

static void release_cells(duckhts_count_cell_t *cells, size_t count) {
    free(cells);
    duckhts_count_bytes_release((uint64_t)count * sizeof(*cells));
}

/* Reads the cells of one BLOB into an owned, aligned array of `count` cells,
 * charged against max_cell_bytes; release_cells gives both back. The BLOB
 * must hold at least one cell, and every cell must be in range. */
static duckhts_count_cell_t *read_cells(duckdb_string_t *blob, size_t *count, uint64_t max_cell_bytes,
                                        char *error) {
    const char *data = duckdb_string_is_inlined(*blob) ? blob->value.inlined.inlined : blob->value.pointer.ptr;
    uint32_t length = blob->value.inlined.length;
    duckhts_count_cells_header_t header;
    duckhts_count_cell_t *cells;
    uint64_t bytes, would_hold = 0;
    if (length < sizeof(header)) {
        snprintf(error, FIT_ERRLEN, "internal error: histogram is too short");
        return NULL;
    }
    memcpy(&header, data, sizeof(header));
    if (memcmp(header.magic, DUCKHTS_COUNT_CELLS_MAGIC, DUCKHTS_COUNT_CELLS_MAGIC_LENGTH) != 0 ||
        (length - sizeof(header)) % sizeof(duckhts_count_cell_t) != 0) {
        snprintf(error, FIT_ERRLEN, "internal error: histogram has an unknown layout");
        return NULL;
    }
    *count = (length - sizeof(header)) / sizeof(duckhts_count_cell_t);
    if (*count == 0) {
        snprintf(error, FIT_ERRLEN, "internal error: histogram has no cell");
        return NULL;
    }
    bytes = (uint64_t)*count * sizeof(*cells);
    if (!duckhts_count_bytes_charge(bytes, max_cell_bytes, &would_hold)) {
        snprintf(error, FIT_ERRLEN,
                 "a fit would hold %llu bytes of histogram memory, more than max_cell_bytes = %llu; "
                 "fit fewer samples per query or raise max_cell_bytes",
                 (unsigned long long)would_hold, (unsigned long long)max_cell_bytes);
        return NULL;
    }
    cells = malloc((size_t)bytes);
    if (!cells) {
        duckhts_count_bytes_release(bytes);
        snprintf(error, FIT_ERRLEN, "out of memory for %llu bytes of a fit", (unsigned long long)bytes);
        return NULL;
    }
    memcpy(cells, data + sizeof(header), *count * sizeof(*cells));
    /* The function takes any BLOB, so every cell is checked as the aggregate
     * checked it (count_cells.c) before the model indexes tables by it. */
    for (size_t i = 0; i < *count; i++) {
        const duckhts_count_cell_t *cell = &cells[i];
        if (cell->depth < 1 || cell->alt > cell->depth || cell->sites < 1 || !(cell->af > 0 && cell->af < 1)) {
            snprintf(error, FIT_ERRLEN, "internal error: histogram cell %llu is out of range",
                     (unsigned long long)i);
            release_cells(cells, *count);
            return NULL;
        }
    }
    return cells;
}

typedef struct {
    uint64_t sites;
    uint64_t reads;
    int blocks;
    duckhts_count_error_fit_t pooled;
    duckhts_count_error_fit_t relative;
    double contamination_sd; /* NAN without enough blocks */
    const char *status;
} fit_result_t;

/* The whole fit of one sample. Returns 1, or 0 with `error` written. */
static int fit_sample(duckhts_count_cell_t *cells, size_t count, int64_t min_block_sites,
                      uint64_t max_cell_bytes, fit_result_t *result, char *error) {
    size_t at = 0;
    double sum = 0, sum_squares = 0;
    memset(result, 0, sizeof(*result));
    result->contamination_sd = NAN;
    qsort(cells, count, sizeof(*cells), compare_cells);
    for (size_t i = 0; i < count; i++) {
        /* sites * depth is below 2^48 per cell; the sums are reported as
         * BIGINT, so a crafted histogram past that is an error, not a wrap. */
        const uint64_t cell_reads = (uint64_t)cells[i].sites * cells[i].depth;
        if (result->sites > (uint64_t)INT64_MAX - cells[i].sites || result->reads > (uint64_t)INT64_MAX - cell_reads) {
            snprintf(error, FIT_ERRLEN, "internal error: histogram sites or reads exceed BIGINT");
            return 0;
        }
        result->sites += cells[i].sites;
        result->reads += cell_reads;
    }
    if (!duckhts_count_error_fit(cells, count, DUCKHTS_COUNT_ERROR_UNRELATED, max_cell_bytes, &result->pooled,
                                 error, FIT_ERRLEN) ||
        !duckhts_count_error_fit(cells, count, DUCKHTS_COUNT_ERROR_RELATIVE, max_cell_bytes, &result->relative,
                                 error, FIT_ERRLEN)) {
        return 0;
    }
    /* Blocks are contiguous after the sort. */
    while (at < count) {
        size_t end = at;
        uint64_t block_sites = 0;
        duckhts_count_error_fit_t fit;
        while (end < count && cells[end].block == cells[at].block) block_sites += cells[end++].sites;
        if (block_sites >= (uint64_t)min_block_sites) {
            if (!duckhts_count_error_fit_block(cells + at, end - at, &result->pooled.params, max_cell_bytes, &fit,
                                               error, FIT_ERRLEN)) {
                return 0;
            }
            sum += fit.params.contamination;
            sum_squares += fit.params.contamination * fit.params.contamination;
            result->blocks++;
        }
        at = end;
    }
    if (result->blocks >= FIT_FEW_BLOCKS) {
        double mean = sum / result->blocks;
        double variance = (sum_squares - result->blocks * mean * mean) / (result->blocks - 1);
        result->contamination_sd = sqrt((variance > 0 ? variance : 0) / result->blocks);
    }
    if (!result->pooled.converged) {
        result->status = "no_convergence";
    } else if (!result->relative.converged) {
        result->status = "no_convergence_relative";
    } else if (result->pooled.at_bound) {
        result->status = "at_bound";
    } else if (result->sites < FIT_FEW_SITES) {
        result->status = "few_sites";
    } else if (result->blocks < FIT_FEW_BLOCKS) {
        result->status = "few_blocks";
    } else {
        result->status = "ok";
    }
    return 1;
}

/* The result of a sample with no usable site: counts of 0, numbers NULL. */
static void no_sites_result(fit_result_t *result) {
    double *pooled = (double *)&result->pooled.params;
    double *relative = (double *)&result->relative.params;
    memset(result, 0, sizeof(*result));
    for (int i = 0; i < DUCKHTS_COUNT_ERROR_PARAMETERS; i++) pooled[i] = relative[i] = NAN;
    result->pooled.log_likelihood = NAN;
    result->relative.log_likelihood = NAN;
    result->contamination_sd = NAN;
    result->status = "no_sites";
}

static void set_double(duckdb_vector vector, idx_t row, double value) {
    if (isnan(value)) {
        duckdb_vector_ensure_validity_writable(vector);
        duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
    } else {
        ((double *)duckdb_vector_get_data(vector))[row] = value;
    }
}

static void count_error_fit_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector cells_vector = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector min_block_vector = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector max_bytes_vector = duckdb_data_chunk_get_vector(input, 2);
    duckdb_string_t *blobs = duckdb_vector_get_data(cells_vector);
    const int64_t *min_block_sites = duckdb_vector_get_data(min_block_vector);
    const int64_t *max_cell_bytes = duckdb_vector_get_data(max_bytes_vector);
    duckdb_vector out[OUT_COUNT];
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (int i = 0; i < OUT_COUNT; i++) out[i] = duckdb_struct_vector_get_child(output, i);

    for (idx_t row = 0; row < rows; row++) {
        char error[FIT_ERRLEN];
        size_t count = 0;
        duckhts_count_cell_t *cells;
        fit_result_t result;
        int fitted;
        if (!row_valid(min_block_vector, row) || min_block_sites[row] < 1) {
            fit_error(info, "min_block_sites must be at least 1");
            return;
        }
        if (!row_valid(max_bytes_vector, row) || max_cell_bytes[row] < 1) {
            fit_error(info, "max_cell_bytes must be at least 1");
            return;
        }
        if (!row_valid(cells_vector, row)) {
            /* The macro gives a sample with no usable site a NULL histogram. */
            no_sites_result(&result);
        } else {
            cells = read_cells(&blobs[row], &count, (uint64_t)max_cell_bytes[row], error);
            if (!cells) {
                fit_error(info, error);
                return;
            }
            fitted = fit_sample(cells, count, min_block_sites[row], (uint64_t)max_cell_bytes[row], &result, error);
            release_cells(cells, count);
            if (!fitted) {
                fit_error(info, error);
                return;
            }
        }
        ((int64_t *)duckdb_vector_get_data(out[OUT_SITES]))[row] = (int64_t)result.sites;
        ((int64_t *)duckdb_vector_get_data(out[OUT_READS]))[row] = (int64_t)result.reads;
        set_double(out[OUT_MEAN_DEPTH], row, result.sites ? (double)result.reads / (double)result.sites : NAN);
        ((int32_t *)duckdb_vector_get_data(out[OUT_BLOCKS]))[row] = result.blocks;
        set_double(out[OUT_SEQ_ERROR], row, result.pooled.params.seq_error);
        set_double(out[OUT_CONTAMINATION], row, result.pooled.params.contamination);
        set_double(out[OUT_CONTAMINATION_SD], row, result.contamination_sd);
        set_double(out[OUT_HOMOZYGOSITY_EXCESS], row, result.pooled.params.homozygosity_excess);
        set_double(out[OUT_ALLELE_BALANCE], row, result.pooled.params.allele_balance);
        set_double(out[OUT_SPREAD_HOM], row, result.pooled.params.spread_hom);
        set_double(out[OUT_SPREAD_HET], row, result.pooled.params.spread_het);
        set_double(out[OUT_ARTEFACT_WEIGHT], row, result.pooled.params.artefact_weight);
        set_double(out[OUT_LOG_LIKELIHOOD], row, result.pooled.log_likelihood);
        set_double(out[OUT_CONTAMINATION_RELATIVE], row, result.relative.params.contamination);
        set_double(out[OUT_LOG_LIKELIHOOD_RELATIVE], row, result.relative.log_likelihood);
        duckdb_vector_assign_string_element(out[OUT_STATUS], row, result.status);
        duckdb_vector_assign_string_element(out[OUT_METHOD], row, FIT_METHOD);
    }
}

static void valid_args_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector columns[6];
    bool *valid = duckdb_vector_get_data(output);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (int i = 0; i < 6; i++) columns[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < rows; row++) {
        int64_t value[6];
        int has[6];
        const char *error;
        for (int i = 0; i < 6; i++) {
            has[i] = row_valid(columns[i], row);
            value[i] = has[i] ? ((const int64_t *)duckdb_vector_get_data(columns[i]))[row] : 0;
        }
        error = duckhts_count_error_check_args(has[0], value[0], has[1], value[1], has[2], value[2], has[3], value[3]);
        if (!error) error = duckhts_count_check_limits(has[4], value[4], has[5], value[5]);
        if (error) {
            fit_error(info, error);
            return;
        }
        valid[row] = true;
    }
}

static duckdb_logical_type fit_result_type(void) {
    duckdb_logical_type fields[OUT_COUNT];
    duckdb_logical_type result;
    for (int i = 0; i < OUT_COUNT; i++) {
        duckdb_type kind = DUCKDB_TYPE_DOUBLE;
        if (i == OUT_SITES || i == OUT_READS) kind = DUCKDB_TYPE_BIGINT;
        if (i == OUT_BLOCKS) kind = DUCKDB_TYPE_INTEGER;
        if (i == OUT_STATUS || i == OUT_METHOD) kind = DUCKDB_TYPE_VARCHAR;
        fields[i] = duckdb_create_logical_type(kind);
    }
    result = duckdb_create_struct_type(fields, output_names, OUT_COUNT);
    for (int i = 0; i < OUT_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    return result;
}

extern bool register_duckhts_count_cells(duckdb_connection connection);

bool register_duckhts_count_error_functions(duckdb_connection connection) {
    duckdb_logical_type blob = duckdb_create_logical_type(DUCKDB_TYPE_BLOB);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type result = fit_result_type();
    duckdb_scalar_function function;
    bool ok;

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, FIT_NAME);
    duckdb_scalar_function_add_parameter(function, blob);
    duckdb_scalar_function_add_parameter(function, bigint); /* min_block_sites */
    duckdb_scalar_function_add_parameter(function, bigint); /* max_cell_bytes */
    duckdb_scalar_function_set_return_type(function, result);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, count_error_fit_scalar);
    ok = duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, VALID_ARGS_NAME);
    for (int i = 0; i < 6; i++) duckdb_scalar_function_add_parameter(function, bigint);
    duckdb_scalar_function_set_return_type(function, boolean);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, valid_args_scalar);
    ok = ok && duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&blob);
    return ok && register_duckhts_count_cells(connection);
}
