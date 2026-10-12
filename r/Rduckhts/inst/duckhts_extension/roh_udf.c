/* Runs-of-homozygosity kernel: one sample and chromosome is decoded per row.
 *
 *   - duckhts_roh_segments(), public, takes DuckDB lists sorted by position.
 *   - __duckhts_roh_decode(), used by the ROH macros, takes the packed site
 *     list that __duckhts_roh_sites() builds (roh_sites.h).
 *   - __duckhts_roh_decode_stable(), used by duckhts_roh_counts_fitted, also
 *     decodes a read-count list at artefact_weight -/+ artefact_weight_sd
 *     and gives each run the share of its bases that both decodes call, and
 *     the bases that the decode at + artefact_weight_sd adds next to it.
 *   - __duckhts_roh_valid_args() checks the model parameters and memory limits
 *     of a macro call without reading a site, so that invalid arguments fail
 *     on empty input too.
 *
 * One row is decoded at a time, so the kernel workspace follows the longest
 * list, not the chunk.
 */
#if defined(__MINGW32__) && !defined(__USE_MINGW_ANSI_STDIO)
#define __USE_MINGW_ANSI_STDIO 1
#endif
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckdb_list.h"
#include "roh_hmm.h"
#include "roh_sites.h"

#define ROH_NAME "duckhts_roh_segments"
#define ROH_ERRLEN 256
#define ROH_MAX_POS1 INT32_MAX

/* Model parameters: the vectors of a call, and the values of one row. A
 * parameter the evidence kind does not use has no vector and stays 0. The
 * read-count overload of duckhts_roh_segments has no allele_balance,
 * artefact_weight or artefact_weight_sd vector and uses 1/2, 0 and no
 * stability decode. */
typedef struct {
    duckdb_vector gt_error, seq_error, contamination, allele_balance, artefact_weight, artefact_weight_sd;
    duckdb_vector map_pos, map_cm, rec_rate, hw_to_az, az_to_hw;
} roh_model_args_t;

typedef struct {
    double hw_to_az, az_to_hw, gt_error, rec_rate;
    duckhts_roh_count_model_t counts;
    /* Read counts only: with a non-NULL artefact_weight_sd the sites are also
     * decoded at artefact_weight - and + artefact_weight_sd (stability_decode). */
    bool stability;
    double artefact_weight_sd;
} roh_model_t;

/* Accepted ranges of the two fitted numbers that only move the read-count
 * emission (see check_model). */
#define ROH_ALLELE_BALANCE_MIN 0.35
#define ROH_ALLELE_BALANCE_MAX 0.65
#define ROH_ARTEFACT_WEIGHT_MAX 0.05

/* Arguments of duckhts_roh_segments. The evidence kind of an overload is
 * attached to it as extra info: PL lists (INTEGER[][]), dosages with a phred
 * gt_error, or read counts with seq_error and contamination. */
typedef struct {
    duckhts_roh_sites_kind_t kind;
    /* evidence is the PL or dosage list, or for counts the other-allele counts;
     * counted holds the counts of the allele af refers to. */
    duckdb_vector positions, af, evidence, counted;
    roh_model_args_t model;
} roh_args_t;

/* Arguments of __duckhts_roh_decode. */
typedef struct {
    duckdb_vector sites;
    roh_model_args_t model;
} roh_packed_args_t;

typedef struct {
    int32_t *map_pos0;
    double *map_rate;
    size_t map_capacity;
} roh_map_t;

typedef struct {
    bool valid;
    const duckdb_list_entry *entries;
    duckdb_vector child;
    const uint64_t *child_validity;
} roh_list_t;

static bool row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static bool element_valid(const uint64_t *mask, idx_t index) {
    return !mask || duckdb_validity_row_is_valid((uint64_t *)mask, index);
}

static roh_list_t open_list(duckdb_vector vector, idx_t row) {
    roh_list_t list = {0};
    list.valid = row_valid(vector, row);
    list.entries = duckdb_vector_get_data(vector);
    list.child = duckdb_list_vector_get_child(vector);
    list.child_validity = duckdb_vector_get_validity(list.child);
    return list;
}

/* MinGW's GCC treats the printf archetype as Microsoft's format rules, which
 * reject %lld/%llu. Use its C99-conforming stdio and check formats as such.
 * Clang (Windows ARM64) has no gnu_printf archetype and accepts %llu as printf. */
#if defined(__MINGW32__) && !defined(__clang__)
#define ROH_PRINTF_ARCHETYPE gnu_printf
#else
#define ROH_PRINTF_ARCHETYPE printf
#endif
static void roh_error(duckdb_function_info info, const char *format, ...)
    __attribute__((format(ROH_PRINTF_ARCHETYPE, 2, 3)));

static void roh_error(duckdb_function_info info, const char *format, ...) {
    char message[ROH_ERRLEN];
    char full[ROH_ERRLEN + sizeof(ROH_NAME) + 4];
    va_list args;
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    snprintf(full, sizeof(full), "%s: %s", ROH_NAME, message);
    duckdb_scalar_function_set_error(info, full);
}

static void set_row_null(duckdb_vector output, idx_t row) {
    duckdb_vector_ensure_validity_writable(output);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
}

static bool read_double_arg(duckdb_vector vector, idx_t row, double *value) {
    if (!row_valid(vector, row)) return false;
    *value = ((const double *)duckdb_vector_get_data(vector))[row];
    return true;
}

/* Reads the model parameters a kind uses. Returns false when a required one is
 * NULL: the row then has a NULL result. A NULL rec_rate means no rate. */
static bool read_model(duckhts_roh_sites_kind_t kind, const roh_model_args_t *args, idx_t row,
                       roh_model_t *model) {
    memset(model, 0, sizeof(*model));
    if (!read_double_arg(args->hw_to_az, row, &model->hw_to_az) ||
        !read_double_arg(args->az_to_hw, row, &model->az_to_hw)) return false;
    if (kind == DUCKHTS_ROH_SITES_GT && !read_double_arg(args->gt_error, row, &model->gt_error)) {
        return false;
    }
    if (kind == DUCKHTS_ROH_SITES_COUNTS) {
        duckhts_roh_count_model_t *counts = &model->counts;
        if (!read_double_arg(args->seq_error, row, &counts->seq_error) ||
            !read_double_arg(args->contamination, row, &counts->contamination)) return false;
        counts->allele_balance = 0.5;
        counts->artefact_weight = 0.0;
        if (args->allele_balance != NULL &&
            !read_double_arg(args->allele_balance, row, &counts->allele_balance)) return false;
        if (args->artefact_weight != NULL &&
            !read_double_arg(args->artefact_weight, row, &counts->artefact_weight)) return false;
        /* A NULL artefact_weight_sd asks for no stability decode; it is not a
         * NULL result. */
        model->stability = args->artefact_weight_sd != NULL &&
                           read_double_arg(args->artefact_weight_sd, row, &model->artefact_weight_sd);
    }
    if (!read_double_arg(args->rec_rate, row, &model->rec_rate)) model->rec_rate = 0;
    return true;
}

/* Returns NULL, or the text of the first invalid model parameter. */
static const char *check_model(duckhts_roh_sites_kind_t kind, const roh_model_t *model) {
    if (!isfinite(model->hw_to_az) || model->hw_to_az < 0 || model->hw_to_az > 1) {
        return "hw_to_az must be in [0, 1]";
    }
    if (!isfinite(model->az_to_hw) || model->az_to_hw < 0 || model->az_to_hw > 1) {
        return "az_to_hw must be in [0, 1]";
    }
    if (kind == DUCKHTS_ROH_SITES_GT && (!isfinite(model->gt_error) || model->gt_error < 0)) {
        return "gt_error must be a finite phred value of at least 0";
    }
    if (kind == DUCKHTS_ROH_SITES_COUNTS) {
        const duckhts_roh_count_model_t *counts = &model->counts;
        if (!isfinite(counts->seq_error) || counts->seq_error <= 0 || counts->seq_error >= 0.5) {
            return "seq_error must be in (0, 0.5)";
        }
        if (!isfinite(counts->contamination) || counts->contamination < 0 ||
            counts->contamination >= 1) {
            return "contamination must be in [0, 1)";
        }
        /* The fit measures allele_balance 0.4974 to 0.4984 and artefact_weight
         * 0.0018 to 0.0030 on 30x samples, and 0.4960 and 0.0046 with a fifth
         * of the reads from another genome (benchmarks/benchmark_count_error_fit.md).
         * The ranges hold those values with wide margins and exclude what the
         * emission cannot take as a per-sample constant: an allele balance
         * beyond 0.35 or 0.65 is 50 times the mapping bias the fit measures
         * and points to copy number, a mixture or a broken panel, which one
         * number per sample does not describe; and an artefact weight
         * above 0.05 caps one site's likelihood ratio (see
         * duckhts_roh_pdg_from_counts) below 1 + 19 (d + 1) at depth d, so
         * that depth barely informs a site. Out-of-range values are errors,
         * not clipped. */
        if (!isfinite(counts->allele_balance) || counts->allele_balance < ROH_ALLELE_BALANCE_MIN ||
            counts->allele_balance > ROH_ALLELE_BALANCE_MAX) {
            return "allele_balance must be in [0.35, 0.65]";
        }
        if (!isfinite(counts->artefact_weight) || counts->artefact_weight < 0 ||
            counts->artefact_weight > ROH_ARTEFACT_WEIGHT_MAX) {
            return "artefact_weight must be in [0, 0.05]";
        }
        if (model->stability && (!isfinite(model->artefact_weight_sd) || model->artefact_weight_sd < 0)) {
            return "artefact_weight_sd must be finite and at least 0";
        }
    }
    if (!isfinite(model->rec_rate) || model->rec_rate < 0) return "rec_rate must be finite and at least 0";
    return NULL;
}

/* bcftools skips a repeated position, a site without a usable frequency
 * (missing, NaN or exactly 0) and a site without usable genotype evidence. A
 * site without a frequency arrives here with af NaN. */
static bool site_enters_model(bool duplicate, bool usable, double af) {
    return !duplicate && usable && !isnan(af) && af != 0.0;
}

static void free_map(roh_map_t *map) {
    free(map->map_pos0);
    free(map->map_rate);
}

/* Loads the map nodes of one row. Returns 1 with *present set, or 0 after
 * reporting an error. */
static int load_map(duckdb_function_info info, const roh_model_args_t *args, idx_t row,
                    roh_map_t *map, size_t *count) {
    roh_list_t pos_list = open_list(args->map_pos, row);
    roh_list_t cm_list = open_list(args->map_cm, row);
    *count = 0;
    if (!pos_list.valid && !cm_list.valid) return 1;
    if (!pos_list.valid || !cm_list.valid) {
        roh_error(info, "genetic map positions and cM values must both be given or both NULL");
        return 0;
    }
    idx_t n = pos_list.entries[row].length;
    if (n != cm_list.entries[row].length) {
        roh_error(info, "genetic map positions and cM values differ in length (%llu vs %llu)",
                  (unsigned long long)n, (unsigned long long)cm_list.entries[row].length);
        return 0;
    }
    if (n == 0) {
        roh_error(info, "genetic map is empty");
        return 0;
    }
    if (n > map->map_capacity) {
        int32_t *pos0 = realloc(map->map_pos0, n * sizeof(*pos0));
        if (pos0 == NULL) goto oom;
        map->map_pos0 = pos0;
        double *rate = realloc(map->map_rate, n * sizeof(*rate));
        if (rate == NULL) goto oom;
        map->map_rate = rate;
        map->map_capacity = n;
    }
    const int64_t *positions = duckdb_vector_get_data(pos_list.child);
    const double *cm = duckdb_vector_get_data(cm_list.child);
    idx_t pos_offset = pos_list.entries[row].offset;
    idx_t cm_offset = cm_list.entries[row].offset;
    for (idx_t i = 0; i < n; i++) {
        if (!element_valid(pos_list.child_validity, pos_offset + i) ||
            !element_valid(cm_list.child_validity, cm_offset + i)) {
            roh_error(info, "genetic map entry %llu has a NULL position or cM value",
                      (unsigned long long)(i + 1));
            return 0;
        }
        int64_t pos1 = positions[pos_offset + i];
        double value = cm[cm_offset + i];
        if (pos1 < 1 || pos1 > ROH_MAX_POS1) {
            roh_error(info, "genetic map position %lld is outside 1..%d", (long long)pos1,
                      ROH_MAX_POS1);
            return 0;
        }
        if (!isfinite(value)) {
            roh_error(info, "genetic map cM values must be finite");
            return 0;
        }
        if (i > 0 && pos1 - 1 <= map->map_pos0[i - 1]) {
            roh_error(info, "genetic map positions must be strictly ascending");
            return 0;
        }
        map->map_pos0[i] = (int32_t)(pos1 - 1);
        map->map_rate[i] = value * 0.01;
    }
    *count = n;
    return 1;
oom:
    roh_error(info, "out of memory for the genetic map");
    return 0;
}

/* Loads the genetic map of a row into `map` and the transition parameters
 * into `params`, which borrows the map. Returns false after reporting an
 * error. */
static bool load_params(duckdb_function_info info, const roh_model_args_t *args, const roh_model_t *model,
                        idx_t row, roh_map_t *map, duckhts_roh_params_t *params) {
    size_t map_count;
    if (!load_map(info, args, row, map, &map_count)) return false;
    params->hw_to_az = model->hw_to_az;
    params->az_to_hw = model->az_to_hw;
    params->rec_rate = model->rec_rate;
    params->map_pos0 = map_count ? map->map_pos0 : NULL;
    params->map_rate = map_count ? map->map_rate : NULL;
    params->map_n = map_count;
    return true;
}

/* Starts a run of at most `sites` sites. Returns false after reporting an
 * error. */
static bool begin_run(duckdb_function_info info, const duckhts_roh_params_t *params, size_t sites,
                      duckhts_roh_t *roh) {
    if (duckhts_roh_begin(roh, params, sites) != 0) {
        roh_error(info, "out of memory for %llu sites", (unsigned long long)sites);
        return false;
    }
    return true;
}

/* Loads the genetic map of a row and starts a run of at most `sites` sites.
 * Returns false after reporting an error. */
static bool begin_decode(duckdb_function_info info, const roh_model_args_t *args,
                         const roh_model_t *model, idx_t row, size_t sites, duckhts_roh_t *roh,
                         roh_map_t *map) {
    duckhts_roh_params_t params;
    return load_params(info, args, model, row, map, &params) && begin_run(info, &params, sites, roh);
}

/* Decodes one row. Returns 1 when the row holds a result (segments are in roh),
 * 0 after an error, and -1 for a NULL result. */
static int decode_row(duckdb_function_info info, const roh_args_t *args, idx_t row,
                      duckhts_roh_t *roh, roh_map_t *map) {
    roh_model_t model;
    if (!read_model(args->kind, &args->model, row, &model)) return -1;
    const char *model_error = check_model(args->kind, &model);
    if (model_error != NULL) {
        roh_error(info, "%s", model_error);
        return 0;
    }

    bool counts = args->kind == DUCKHTS_ROH_SITES_COUNTS;
    roh_list_t positions = open_list(args->positions, row);
    roh_list_t af = open_list(args->af, row);
    roh_list_t evidence = open_list(args->evidence, row);
    roh_list_t counted = {0};
    if (counts) counted = open_list(args->counted, row);
    if (!positions.valid || !af.valid || !evidence.valid || (counts && !counted.valid)) return -1;
    idx_t n = positions.entries[row].length;
    const char *evidence_name = args->kind == DUCKHTS_ROH_SITES_GT ? "genotype"
                                : counts                       ? "other-allele count"
                                                               : "PL";
    if (af.entries[row].length != n) {
        roh_error(info, "positions and af lists differ in length (%llu vs %llu)",
                  (unsigned long long)n, (unsigned long long)af.entries[row].length);
        return 0;
    }
    if (evidence.entries[row].length != n) {
        roh_error(info, "positions and %s lists differ in length (%llu vs %llu)", evidence_name,
                  (unsigned long long)n, (unsigned long long)evidence.entries[row].length);
        return 0;
    }
    if (counts && counted.entries[row].length != n) {
        roh_error(info, "positions and counted-allele count lists differ in length (%llu vs %llu)",
                  (unsigned long long)n, (unsigned long long)counted.entries[row].length);
        return 0;
    }

    if (!begin_decode(info, &args->model, &model, row, n, roh, map)) return 0;

    const int64_t *pos = duckdb_vector_get_data(positions.child);
    const double *freq = duckdb_vector_get_data(af.child);
    idx_t pos_offset = positions.entries[row].offset;
    idx_t af_offset = af.entries[row].offset;
    idx_t ev_offset = evidence.entries[row].offset;

    /* PL evidence is a list of lists; GT evidence and counts are flat lists. */
    const duckdb_list_entry *pl_entries = NULL;
    const int32_t *pl_values = NULL;
    const uint64_t *pl_validity = NULL;
    const int32_t *flat = NULL;
    const int32_t *counted_values = NULL;
    idx_t counted_offset = 0;
    if (args->kind != DUCKHTS_ROH_SITES_PL) {
        flat = duckdb_vector_get_data(evidence.child);
    }
    if (counts) {
        counted_values = duckdb_vector_get_data(counted.child);
        counted_offset = counted.entries[row].offset;
    }
    if (args->kind == DUCKHTS_ROH_SITES_PL) {
        pl_entries = duckdb_vector_get_data(evidence.child);
        duckdb_vector pl_child = duckdb_list_vector_get_child(evidence.child);
        pl_values = duckdb_vector_get_data(pl_child);
        pl_validity = duckdb_vector_get_validity(pl_child);
    }

    int64_t prev_pos1 = 0;
    for (idx_t i = 0; i < n; i++) {
        if (!element_valid(positions.child_validity, pos_offset + i)) {
            roh_error(info, "positions cannot contain NULL (entry %llu)", (unsigned long long)(i + 1));
            return 0;
        }
        int64_t pos1 = pos[pos_offset + i];
        if (pos1 < 1 || pos1 > ROH_MAX_POS1) {
            roh_error(info, "position %lld is outside 1..%d", (long long)pos1, ROH_MAX_POS1);
            return 0;
        }
        if (pos1 < prev_pos1) {
            roh_error(info, "positions must be sorted ascending (entry %llu is %lld after %lld)",
                      (unsigned long long)(i + 1), (long long)pos1, (long long)prev_pos1);
            return 0;
        }
        bool duplicate = pos1 == prev_pos1;
        prev_pos1 = pos1;

        /* A valid AF element outside [0, 1] is an error even on skipped sites. */
        bool has_af = element_valid(af.child_validity, af_offset + i);
        double alt_freq = has_af ? freq[af_offset + i] : 0;
        if (has_af && !isnan(alt_freq) && (alt_freq < 0 || alt_freq > 1)) {
            roh_error(info, "af must be in [0, 1] (entry %llu is %g)", (unsigned long long)(i + 1),
                      alt_freq);
            return 0;
        }

        double pdg[3];
        bool usable;
        if (args->kind == DUCKHTS_ROH_SITES_GT) {
            idx_t at = ev_offset + i;
            usable = element_valid(evidence.child_validity, at);
            if (usable) {
                int32_t value = flat[at];
                if (value < 0 || value > 2) {
                    roh_error(info, "genotype dosage must be 0, 1 or 2 (entry %llu is %d)",
                              (unsigned long long)(i + 1), value);
                    return 0;
                }
                usable = duckhts_roh_pdg_from_gt(model.gt_error, value, pdg);
            }
        } else if (counts) {
            idx_t at = ev_offset + i;
            idx_t counted_at = counted_offset + i;
            bool other_valid = element_valid(evidence.child_validity, at);
            bool counted_valid = element_valid(counted.child_validity, counted_at);
            /* Each count present is checked on its own, so a negative count is
             * an error even when the other count of its site is NULL. */
            if (other_valid && flat[at] < 0) {
                roh_error(info, "read counts must be at least 0 (entry %llu other-allele count is %d)",
                          (unsigned long long)(i + 1), flat[at]);
                return 0;
            }
            if (counted_valid && counted_values[counted_at] < 0) {
                roh_error(info, "read counts must be at least 0 (entry %llu counted-allele count is %d)",
                          (unsigned long long)(i + 1), counted_values[counted_at]);
                return 0;
            }
            usable = other_valid && counted_valid;
            if (usable) {
                int32_t other_count = flat[at];
                int32_t counted_count = counted_values[counted_at];
                /* The read model needs the site's frequency; sites without a
                 * usable one are skipped below before the emission is used. */
                usable = has_af && !isnan(alt_freq) && alt_freq > 0 &&
                         duckhts_roh_pdg_from_counts(other_count, counted_count, &model.counts,
                                                     alt_freq, pdg);
            }
        } else {
            idx_t at = ev_offset + i;
            usable = element_valid(evidence.child_validity, at) && pl_entries[at].length == 3;
            if (usable) {
                idx_t base = pl_entries[at].offset;
                usable = element_valid(pl_validity, base) && element_valid(pl_validity, base + 1) &&
                         element_valid(pl_validity, base + 2);
                if (usable) {
                    usable = duckhts_roh_pdg_from_pl(roh, pl_values[base], pl_values[base + 1],
                                                     pl_values[base + 2], pdg);
                }
            }
        }

        if (!site_enters_model(duplicate, usable, has_af ? alt_freq : NAN)) continue;
        duckhts_roh_push(roh, (int32_t)(pos1 - 1), pdg, alt_freq);
    }
    if (duckhts_roh_decode(roh) != 0) {
        roh_error(info, "out of memory while decoding");
        return 0;
    }
    return 1;
}

/* A packed site list (roh_sites.h) whose header has been checked; its records
 * are checked as they are decoded. */
typedef struct {
    duckhts_roh_sites_kind_t kind;
    const unsigned char *records;
    size_t size; /* bytes per record */
    size_t n;    /* records */
} roh_packed_list_t;

/* Opens the packed site list of a row. Same return values as decode_row. The
 * list comes from __duckhts_roh_sites(), but the function is callable with any
 * BLOB, so the header is checked here and every record before it is used. */
static int open_packed_list(duckdb_function_info info, duckdb_vector sites, idx_t row, roh_packed_list_t *list) {
    if (!row_valid(sites, row)) return -1;
    duckdb_string_t *blobs = duckdb_vector_get_data(sites);
    const unsigned char *data = (const unsigned char *)duckdb_string_t_data(&blobs[row]);
    size_t length = duckdb_string_t_length(blobs[row]);

    duckhts_roh_sites_header_t header;
    if (length < sizeof(header)) {
        roh_error(info, "sites is not a packed site list");
        return 0;
    }
    memcpy(&header, data, sizeof(header));
    size_t size = (size_t)duckhts_roh_site_size(header.kind);
    if (memcmp(header.magic, DUCKHTS_ROH_SITES_MAGIC, DUCKHTS_ROH_SITES_MAGIC_LENGTH) != 0 ||
        size == 0 || (length - sizeof(header)) % size != 0) {
        roh_error(info, "sites is not a packed site list");
        return 0;
    }
    list->kind = (duckhts_roh_sites_kind_t)header.kind;
    list->size = size;
    list->n = (length - sizeof(header)) / size;
    list->records = data + sizeof(header);
    /* The aggregate never builds a longer list; a hand-made one must not size
     * the decoder workspace either. */
    if (list->n > DUCKHTS_ROH_MAX_SITES) {
        roh_error(info, "sites holds more than %d sites", DUCKHTS_ROH_MAX_SITES);
        return 0;
    }
    return 1;
}

/* Decodes the sites of a packed list with the emission of `model`: starts a
 * run, pushes the usable sites and decodes them, so that the segments are in
 * `roh`. Returns false after reporting an error. */
static bool decode_packed_list(duckdb_function_info info, const roh_packed_list_t *list,
                               const roh_model_t *model, const duckhts_roh_params_t *params,
                               duckhts_roh_t *roh) {
    const duckhts_roh_sites_kind_t kind = list->kind;
    if (!begin_run(info, params, list->n, roh)) return false;

    uint32_t prev_pos1 = 0;
    for (size_t i = 0; i < list->n; i++) {
        /* The BLOB is not aligned for the records, so each one is copied out. */
        duckhts_roh_site_t site;
        duckhts_roh_count_site_t count_site;
        double af;
        uint32_t pos1;
        if (kind == DUCKHTS_ROH_SITES_COUNTS) {
            memcpy(&count_site, list->records + i * list->size, sizeof(count_site));
            af = count_site.af;
            pos1 = count_site.pos1;
        } else {
            memcpy(&site, list->records + i * list->size, sizeof(site));
            af = site.af;
            pos1 = site.pos1;
        }
        if (pos1 < 1 || pos1 > (uint32_t)ROH_MAX_POS1 || pos1 < prev_pos1 ||
            (!isnan(af) && (af < 0 || af > 1))) {
            roh_error(info, "packed site %llu is out of order or out of range",
                      (unsigned long long)(i + 1));
            return 0;
        }
        bool duplicate = pos1 == prev_pos1;
        prev_pos1 = pos1;

        double pdg[3];
        bool usable;
        if (kind == DUCKHTS_ROH_SITES_GT) {
            if (site.evidence[0] > 2) {
                roh_error(info, "packed site %llu has dosage %d", (unsigned long long)(i + 1),
                          (int)site.evidence[0]);
                return false;
            }
            usable = site.usable && duckhts_roh_pdg_from_gt(model->gt_error, site.evidence[0], pdg);
        } else if (kind == DUCKHTS_ROH_SITES_COUNTS) {
            if (count_site.usable && (count_site.other < 0 || count_site.counted < 0)) {
                roh_error(info, "packed site %llu has a negative read count",
                          (unsigned long long)(i + 1));
                return false;
            }
            /* The read model needs the site's frequency. */
            usable = count_site.usable && !isnan(af) && af > 0 &&
                     duckhts_roh_pdg_from_counts(count_site.other, count_site.counted,
                                                 &model->counts, af, pdg);
        } else {
            usable = site.usable != 0;
            if (usable) duckhts_roh_pdg_from_capped_pl(roh, site.evidence, pdg);
        }

        if (!site_enters_model(duplicate, usable, af)) continue;
        duckhts_roh_push(roh, (int32_t)(pos1 - 1), pdg, af);
    }
    if (duckhts_roh_decode(roh) != 0) {
        roh_error(info, "out of memory while decoding");
        return false;
    }
    return true;
}

/* The runs of one decode, copied out of the kernel before the next decode of
 * the same sites replaces them. */
typedef struct {
    duckhts_roh_segment_t *items;
    size_t count;
    size_t capacity;
} roh_run_copy_t;

/* Work space of the stability decode, reused across rows: the runs of the two
 * perturbed decodes, and one fraction and one gained-bases count per run of
 * the central decode. Each array holds at most one entry per site of the
 * longest list, so the site limits bound it. */
typedef struct {
    roh_run_copy_t low;
    roh_run_copy_t high;
    double *fractions;
    int64_t *gained;
    size_t run_capacity; /* entries of fractions and of gained */
    bool active;         /* the row's runs have fractions and gained bases */
} roh_stability_t;

static void free_stability(roh_stability_t *stability) {
    free(stability->low.items);
    free(stability->high.items);
    free(stability->fractions);
    free(stability->gained);
}

/* Makes room for one fraction and one gained count per central run. Returns
 * false after reporting an error. */
static bool reserve_run_values(duckdb_function_info info, roh_stability_t *stability, size_t count) {
    if (count <= stability->run_capacity) return true;
    double *fractions = realloc(stability->fractions, count * sizeof(*fractions));
    if (fractions != NULL) stability->fractions = fractions;
    int64_t *gained = fractions != NULL ? realloc(stability->gained, count * sizeof(*gained)) : NULL;
    if (gained == NULL) {
        roh_error(info, "out of memory for %llu runs", (unsigned long long)count);
        return false;
    }
    stability->gained = gained;
    stability->run_capacity = count;
    return true;
}

/* Copies the runs of the last decode. Returns false after reporting an error. */
static bool copy_runs(duckdb_function_info info, const duckhts_roh_t *roh, roh_run_copy_t *copy) {
    size_t count;
    const duckhts_roh_segment_t *segments = duckhts_roh_segments(roh, &count);
    if (count > copy->capacity) {
        duckhts_roh_segment_t *grown = realloc(copy->items, count * sizeof(*grown));
        if (grown == NULL) {
            roh_error(info, "out of memory for %llu runs", (unsigned long long)count);
            return false;
        }
        copy->items = grown;
        copy->capacity = count;
    }
    if (count > 0) memcpy(copy->items, segments, count * sizeof(*segments));
    copy->count = count;
    return true;
}

static int64_t max3(int64_t a, int64_t b, int64_t c) {
    int64_t m = a > b ? a : b;
    return m > c ? m : c;
}

static int64_t min3(int64_t a, int64_t b, int64_t c) {
    int64_t m = a < b ? a : b;
    return m < c ? m : c;
}

/* For each run of `central`, the share of its bases (one-based, inclusive)
 * that are also in a run of `low` and in a run of `high`. Each list is in
 * position order and its runs do not overlap, so the pieces where a run of
 * low and a run of high overlap are disjoint, and one merge pass over the
 * three lists sums them. */
static void stable_fractions(const duckhts_roh_segment_t *central, size_t count, const roh_run_copy_t *low,
                             const roh_run_copy_t *high, double *fractions) {
    size_t first_low = 0, first_high = 0;
    for (size_t r = 0; r < count; r++) {
        const int64_t start = central[r].start_pos1, end = central[r].end_pos1;
        int64_t shared = 0;
        /* Runs that end before this run cannot meet it or any later one. */
        while (first_low < low->count && low->items[first_low].end_pos1 < start) first_low++;
        while (first_high < high->count && high->items[first_high].end_pos1 < start) first_high++;
        size_t a = first_low, b = first_high;
        while (a < low->count && b < high->count && low->items[a].start_pos1 <= end &&
               high->items[b].start_pos1 <= end) {
            const duckhts_roh_segment_t *l = &low->items[a], *h = &high->items[b];
            const int64_t from = max3(start, l->start_pos1, h->start_pos1);
            const int64_t to = min3(end, l->end_pos1, h->end_pos1);
            if (to >= from) shared += to - from + 1;
            if (l->end_pos1 < h->end_pos1) {
                a++;
            } else {
                b++;
            }
        }
        fractions[r] = (double)shared / (double)(end - start + 1);
    }
}

/* For each run of `central`, the bases that runs of `high` add next to it.
 * A run of high meets a run of central when they share a base or abut (no
 * base between them). Each base of a run of high that meets some run of
 * central, and that lies outside every run of central, is given to one run
 * of central: a base before the first run of central the high run meets goes
 * to that run; any other goes to the last run of central before it. A run of
 * high that meets no run of central gives its bases to no run. Both lists are
 * in position order and their runs do not overlap, so one pass with a
 * pointer into central serves every run of high. */
static void gained_bases(const duckhts_roh_segment_t *central, size_t count, const roh_run_copy_t *high,
                         int64_t *gained) {
    size_t first = 0;
    for (size_t r = 0; r < count; r++) gained[r] = 0;
    for (size_t k = 0; k < high->count; k++) {
        const int64_t start = high->items[k].start_pos1, end = high->items[k].end_pos1;
        /* Central runs that end more than one base before this run meet
         * neither it nor a later one. */
        while (first < count && (int64_t)central[first].end_pos1 + 1 < start) first++;
        if (first == count || (int64_t)central[first].start_pos1 - 1 > end) continue;
        if (central[first].start_pos1 > start) gained[first] += central[first].start_pos1 - start;
        for (size_t r = first; r < count && (int64_t)central[r].start_pos1 - 1 <= end; r++) {
            const bool next_meets = r + 1 < count && (int64_t)central[r + 1].start_pos1 - 1 <= end;
            const int64_t until = next_meets ? (int64_t)central[r + 1].start_pos1 - 1 : end;
            if (until > central[r].end_pos1) gained[r] += until - central[r].end_pos1;
        }
    }
}

/* An artefact weight moved by `shift`, kept in the range check_model accepts. */
static double shifted_artefact_weight(double weight, double shift) {
    const double value = weight + shift;
    if (value < 0) return 0;
    if (value > ROH_ARTEFACT_WEIGHT_MAX) return ROH_ARTEFACT_WEIGHT_MAX;
    return value;
}

/* The stability decode of one list: the sites decoded at artefact_weight
 * - and + artefact_weight_sd (each kept in [0, 0.05]), then at
 * artefact_weight, so that the central runs are in `roh` and each has the
 * share of its bases that both perturbed decodes also call
 * (stable_fractions) and the bases the decode at + artefact_weight_sd adds
 * next to it (gained_bases). The three decodes share the list and the kernel
 * workspace. Returns false after reporting an error. */
static bool stability_decode(duckdb_function_info info, const roh_packed_list_t *list, const roh_model_t *model,
                             const duckhts_roh_params_t *params, duckhts_roh_t *roh,
                             roh_stability_t *stability) {
    roh_model_t perturbed = *model;
    const double weight = model->counts.artefact_weight;
    perturbed.counts.artefact_weight = shifted_artefact_weight(weight, -model->artefact_weight_sd);
    if (!decode_packed_list(info, list, &perturbed, params, roh) || !copy_runs(info, roh, &stability->low)) {
        return false;
    }
    perturbed.counts.artefact_weight = shifted_artefact_weight(weight, model->artefact_weight_sd);
    if (!decode_packed_list(info, list, &perturbed, params, roh) || !copy_runs(info, roh, &stability->high)) {
        return false;
    }
    if (!decode_packed_list(info, list, model, params, roh)) return false;

    size_t count;
    const duckhts_roh_segment_t *central = duckhts_roh_segments(roh, &count);
    if (!reserve_run_values(info, stability, count)) return false;
    stable_fractions(central, count, &stability->low, &stability->high, stability->fractions);
    gained_bases(central, count, &stability->high, stability->gained);
    return true;
}

/* Decodes one packed site list (roh_sites.h), with the stability decode when
 * the row has an artefact_weight_sd. Same return values as decode_row. */
static int decode_packed_row(duckdb_function_info info, const roh_packed_args_t *args, idx_t row,
                             duckhts_roh_t *roh, roh_map_t *map, roh_stability_t *stability) {
    roh_packed_list_t list;
    int opened = open_packed_list(info, args->sites, row, &list);
    if (opened != 1) return opened;

    roh_model_t model;
    if (!read_model(list.kind, &args->model, row, &model)) return -1;
    const char *model_error = check_model(list.kind, &model);
    if (model_error != NULL) {
        roh_error(info, "%s", model_error);
        return 0;
    }
    duckhts_roh_params_t params;
    if (!load_params(info, &args->model, &model, row, map, &params)) return 0;
    stability->active = model.stability;
    if (model.stability) return stability_decode(info, &list, &model, &params, roh, stability) ? 1 : 0;
    return decode_packed_list(info, &list, &model, &params, roh) ? 1 : 0;
}

/* Writes the runs of the last decode as the list of `row`. With
 * `stability_members` the run struct has a fifth and sixth member
 * (__duckhts_roh_decode_stable), stable_fraction and gained_bases, which get
 * the values of `stability`, or NULL when that is NULL. */
static bool write_segments(duckdb_function_info info, duckdb_vector output, idx_t row,
                           const duckhts_roh_t *roh, bool stability_members, const roh_stability_t *stability) {
    size_t count;
    const duckhts_roh_segment_t *segments = duckhts_roh_segments(roh, &count);
    duckdb_list_entry entry;
    if (!duckhts_list_extend(output, count, &entry)) {
        roh_error(info, "output allocation failed");
        return false;
    }
    duckdb_vector child = duckdb_list_vector_get_child(output);
    int64_t *starts = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 0));
    int64_t *ends = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 1));
    int32_t *markers = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 2));
    double *quality = duckdb_vector_get_data(duckdb_struct_vector_get_child(child, 3));
    for (size_t i = 0; i < count; i++) {
        starts[entry.offset + i] = segments[i].start_pos1;
        ends[entry.offset + i] = segments[i].end_pos1;
        markers[entry.offset + i] = segments[i].n_markers;
        quality[entry.offset + i] = segments[i].quality;
    }
    if (stability_members) {
        duckdb_vector stable = duckdb_struct_vector_get_child(child, 4);
        duckdb_vector gained = duckdb_struct_vector_get_child(child, 5);
        double *fractions = duckdb_vector_get_data(stable);
        int64_t *bases = duckdb_vector_get_data(gained);
        if (stability == NULL && count > 0) {
            duckdb_vector_ensure_validity_writable(stable);
            duckdb_vector_ensure_validity_writable(gained);
        }
        for (size_t i = 0; i < count; i++) {
            if (stability != NULL) {
                fractions[entry.offset + i] = stability->fractions[i];
                bases[entry.offset + i] = stability->gained[i];
            } else {
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(stable), entry.offset + i);
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(gained), entry.offset + i);
            }
        }
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(output))[row] = entry;
    return true;
}

static void roh_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    roh_args_t args = {0};
    args.kind = (duckhts_roh_sites_kind_t)(uintptr_t)duckdb_scalar_function_get_extra_info(info);
    if (args.kind != DUCKHTS_ROH_SITES_PL && args.kind != DUCKHTS_ROH_SITES_GT &&
        args.kind != DUCKHTS_ROH_SITES_COUNTS) {
        roh_error(info, "internal error: unknown evidence kind");
        return;
    }
    idx_t at = 0;
    args.positions = duckdb_data_chunk_get_vector(input, at++);
    args.af = duckdb_data_chunk_get_vector(input, at++);
    args.evidence = duckdb_data_chunk_get_vector(input, at++);
    if (args.kind == DUCKHTS_ROH_SITES_GT) {
        args.model.gt_error = duckdb_data_chunk_get_vector(input, at++);
    }
    if (args.kind == DUCKHTS_ROH_SITES_COUNTS) {
        args.counted = duckdb_data_chunk_get_vector(input, at++);
        args.model.seq_error = duckdb_data_chunk_get_vector(input, at++);
        args.model.contamination = duckdb_data_chunk_get_vector(input, at++);
    }
    args.model.map_pos = duckdb_data_chunk_get_vector(input, at++);
    args.model.map_cm = duckdb_data_chunk_get_vector(input, at++);
    args.model.rec_rate = duckdb_data_chunk_get_vector(input, at++);
    args.model.hw_to_az = duckdb_data_chunk_get_vector(input, at++);
    args.model.az_to_hw = duckdb_data_chunk_get_vector(input, at++);

    duckhts_roh_t *roh = duckhts_roh_create();
    roh_map_t map = {0};
    if (roh == NULL) {
        roh_error(info, "out of memory");
        return;
    }
    duckdb_list_vector_set_size(output, 0);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        int status = decode_row(info, &args, row, roh, &map);
        if (status == 0) break;
        if (status < 0) {
            set_row_null(output, row);
            continue;
        }
        if (!write_segments(info, output, row, roh, false, NULL)) break;
    }
    free_map(&map);
    duckhts_roh_destroy(roh);
}

/* __duckhts_roh_decode(sites, gt_error, seq_error, contamination,
 * allele_balance, artefact_weight, map_pos, map_cm, rec_rate, hw_to_az,
 * az_to_hw), and __duckhts_roh_decode_stable with artefact_weight_sd after
 * artefact_weight. The evidence kind is in the site list. The runs of
 * __duckhts_roh_decode_stable have a fifth member, stable_fraction: NULL
 * unless the list holds read counts and artefact_weight_sd is not NULL
 * (stability_decode). The extra info of the function tells the two apart. */
static void roh_packed_scalar(duckdb_function_info info, duckdb_data_chunk input,
                              duckdb_vector output) {
    const bool stable = (uintptr_t)duckdb_scalar_function_get_extra_info(info) != 0;
    roh_packed_args_t args = {0};
    idx_t at = 0;
    args.sites = duckdb_data_chunk_get_vector(input, at++);
    args.model.gt_error = duckdb_data_chunk_get_vector(input, at++);
    args.model.seq_error = duckdb_data_chunk_get_vector(input, at++);
    args.model.contamination = duckdb_data_chunk_get_vector(input, at++);
    args.model.allele_balance = duckdb_data_chunk_get_vector(input, at++);
    args.model.artefact_weight = duckdb_data_chunk_get_vector(input, at++);
    if (stable) args.model.artefact_weight_sd = duckdb_data_chunk_get_vector(input, at++);
    args.model.map_pos = duckdb_data_chunk_get_vector(input, at++);
    args.model.map_cm = duckdb_data_chunk_get_vector(input, at++);
    args.model.rec_rate = duckdb_data_chunk_get_vector(input, at++);
    args.model.hw_to_az = duckdb_data_chunk_get_vector(input, at++);
    args.model.az_to_hw = duckdb_data_chunk_get_vector(input, at++);

    duckhts_roh_t *roh = duckhts_roh_create();
    roh_map_t map = {0};
    roh_stability_t stability = {0};
    if (roh == NULL) {
        roh_error(info, "out of memory");
        return;
    }
    duckdb_list_vector_set_size(output, 0);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        int status = decode_packed_row(info, &args, row, roh, &map, &stability);
        if (status == 0) break;
        if (status < 0) {
            set_row_null(output, row);
            continue;
        }
        if (!write_segments(info, output, row, roh, stable, stability.active ? &stability : NULL)) break;
    }
    free_stability(&stability);
    free_map(&map);
    duckhts_roh_destroy(roh);
}

/* __duckhts_roh_valid_args(kind, gt_error, seq_error, contamination,
 * allele_balance, artefact_weight, artefact_weight_sd, rec_rate, hw_to_az,
 * az_to_hw, max_sites, max_site_bytes): TRUE, or the error the decode would
 * raise. A NULL model parameter is not an error: the decode returns NULL, or
 * for artefact_weight_sd makes no stability decode. */
static void roh_valid_args_scalar(duckdb_function_info info, duckdb_data_chunk input,
                                  duckdb_vector output) {
    roh_model_args_t args = {0};
    idx_t at = 0;
    duckdb_vector kinds = duckdb_data_chunk_get_vector(input, at++);
    args.gt_error = duckdb_data_chunk_get_vector(input, at++);
    args.seq_error = duckdb_data_chunk_get_vector(input, at++);
    args.contamination = duckdb_data_chunk_get_vector(input, at++);
    args.allele_balance = duckdb_data_chunk_get_vector(input, at++);
    args.artefact_weight = duckdb_data_chunk_get_vector(input, at++);
    args.artefact_weight_sd = duckdb_data_chunk_get_vector(input, at++);
    args.rec_rate = duckdb_data_chunk_get_vector(input, at++);
    args.hw_to_az = duckdb_data_chunk_get_vector(input, at++);
    args.az_to_hw = duckdb_data_chunk_get_vector(input, at++);
    duckdb_vector max_sites = duckdb_data_chunk_get_vector(input, at++);
    duckdb_vector max_site_bytes = duckdb_data_chunk_get_vector(input, at++);

    bool *valid = duckdb_vector_get_data(output);
    idx_t rows = duckdb_data_chunk_get_size(input);
    for (idx_t row = 0; row < rows; row++) {
        bool has_max_sites = row_valid(max_sites, row);
        bool has_max_site_bytes = row_valid(max_site_bytes, row);
        const char *limit_error = duckhts_roh_check_limits(
            has_max_sites, has_max_sites ? ((const int64_t *)duckdb_vector_get_data(max_sites))[row] : 0,
            has_max_site_bytes,
            has_max_site_bytes ? ((const int64_t *)duckdb_vector_get_data(max_site_bytes))[row] : 0);
        if (limit_error != NULL) {
            roh_error(info, "%s", limit_error);
            return;
        }
        int32_t kind = row_valid(kinds, row) ? ((const int32_t *)duckdb_vector_get_data(kinds))[row] : 0;
        if (duckhts_roh_site_size((uint32_t)kind) == 0) {
            roh_error(info, "internal error: unknown evidence kind");
            return;
        }
        roh_model_t model;
        if (read_model((duckhts_roh_sites_kind_t)kind, &args, row, &model)) {
            const char *model_error = check_model((duckhts_roh_sites_kind_t)kind, &model);
            if (model_error != NULL) {
                roh_error(info, "%s", model_error);
                return;
            }
        }
        valid[row] = true;
    }
}

/* Parameters, in order: positions, af, then the evidence of the kind (PL lists;
 * dosages and gt_error; other-allele counts, counted-allele counts, seq_error
 * and contamination), then map_pos, map_cm, rec_rate, hw_to_az and az_to_hw. */
static bool register_overload(duckdb_connection connection, duckhts_roh_sites_kind_t kind) {
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type bigint_list = duckdb_create_list_type(bigint);
    duckdb_logical_type real_list = duckdb_create_list_type(real);
    duckdb_logical_type integer_list = duckdb_create_list_type(integer);
    duckdb_logical_type pl_list = duckdb_create_list_type(integer_list);
    duckdb_logical_type member_types[4] = {bigint, bigint, integer, real};
    const char *member_names[4] = {"start", "end", "n_markers", "quality"};
    duckdb_logical_type segment = duckdb_create_struct_type(member_types, member_names, 4);
    duckdb_logical_type segments = duckdb_create_list_type(segment);

    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, ROH_NAME);
    duckdb_scalar_function_add_parameter(function, bigint_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    if (kind == DUCKHTS_ROH_SITES_PL) {
        duckdb_scalar_function_add_parameter(function, pl_list);
    } else if (kind == DUCKHTS_ROH_SITES_GT) {
        duckdb_scalar_function_add_parameter(function, integer_list);
        duckdb_scalar_function_add_parameter(function, real);
    } else {
        duckdb_scalar_function_add_parameter(function, integer_list);
        duckdb_scalar_function_add_parameter(function, integer_list);
        duckdb_scalar_function_add_parameter(function, real);
        duckdb_scalar_function_add_parameter(function, real);
    }
    duckdb_scalar_function_add_parameter(function, bigint_list);
    duckdb_scalar_function_add_parameter(function, real_list);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_set_return_type(function, segments);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_extra_info(function, (void *)(uintptr_t)kind, NULL);
    duckdb_scalar_function_set_function(function, roh_scalar);
    bool ok = duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&segments);
    duckdb_destroy_logical_type(&segment);
    duckdb_destroy_logical_type(&pl_list);
    duckdb_destroy_logical_type(&integer_list);
    duckdb_destroy_logical_type(&real_list);
    duckdb_destroy_logical_type(&bigint_list);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&bigint);
    return ok;
}

/* The runs of __duckhts_roh_decode, those of duckhts_roh_segments; with
 * `stable`, the runs of __duckhts_roh_decode_stable, which add
 * stable_fraction and gained_bases (write_segments). */
static duckdb_logical_type segments_type(bool stable) {
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type member_types[6] = {bigint, bigint, integer, real, real, bigint};
    const char *member_names[6] = {"start", "end", "n_markers", "quality", "stable_fraction", "gained_bases"};
    duckdb_logical_type segment = duckdb_create_struct_type(member_types, member_names, stable ? 6 : 4);
    duckdb_logical_type segments = duckdb_create_list_type(segment);
    duckdb_destroy_logical_type(&segment);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&bigint);
    return segments;
}

/* __duckhts_roh_decode, or with `stable` __duckhts_roh_decode_stable (see
 * roh_packed_scalar). */
static bool register_packed_decode(duckdb_connection connection, bool stable) {
    duckdb_logical_type blob = duckdb_create_logical_type(DUCKDB_TYPE_BLOB);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type bigint_list = duckdb_create_list_type(bigint);
    duckdb_logical_type real_list = duckdb_create_list_type(real);
    duckdb_logical_type segments = segments_type(stable);

    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, stable ? "__duckhts_roh_decode_stable" : "__duckhts_roh_decode");
    duckdb_scalar_function_add_parameter(function, blob);        /* sites */
    duckdb_scalar_function_add_parameter(function, real);        /* gt_error */
    duckdb_scalar_function_add_parameter(function, real);        /* seq_error */
    duckdb_scalar_function_add_parameter(function, real);        /* contamination */
    duckdb_scalar_function_add_parameter(function, real);        /* allele_balance */
    duckdb_scalar_function_add_parameter(function, real);        /* artefact_weight */
    if (stable) duckdb_scalar_function_add_parameter(function, real); /* artefact_weight_sd */
    duckdb_scalar_function_add_parameter(function, bigint_list); /* map_pos */
    duckdb_scalar_function_add_parameter(function, real_list);   /* map_cm */
    duckdb_scalar_function_add_parameter(function, real);        /* rec_rate */
    duckdb_scalar_function_add_parameter(function, real);        /* hw_to_az */
    duckdb_scalar_function_add_parameter(function, real);        /* az_to_hw */
    duckdb_scalar_function_set_return_type(function, segments);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_extra_info(function, (void *)(uintptr_t)stable, NULL);
    duckdb_scalar_function_set_function(function, roh_packed_scalar);
    bool ok = duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&segments);
    duckdb_destroy_logical_type(&real_list);
    duckdb_destroy_logical_type(&bigint_list);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&blob);
    return ok;
}

static bool register_valid_args(duckdb_connection connection) {
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type real = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "__duckhts_roh_valid_args");
    duckdb_scalar_function_add_parameter(function, integer); /* kind */
    /* gt_error, seq_error, contamination, allele_balance, artefact_weight,
     * artefact_weight_sd, rec_rate, hw_to_az, az_to_hw */
    for (int i = 0; i < 9; i++) duckdb_scalar_function_add_parameter(function, real);
    duckdb_scalar_function_add_parameter(function, bigint); /* max_sites */
    duckdb_scalar_function_add_parameter(function, bigint); /* max_site_bytes */
    duckdb_scalar_function_set_return_type(function, boolean);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_function(function, roh_valid_args_scalar);
    bool ok = duckdb_register_scalar_function(connection, function) == DuckDBSuccess;
    duckdb_destroy_scalar_function(&function);

    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&real);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&integer);
    return ok;
}

/* The site-collecting aggregate is in roh_sites.c. */
extern bool register_duckhts_roh_sites(duckdb_connection connection);

bool register_duckhts_roh_functions(duckdb_connection connection) {
    return register_overload(connection, DUCKHTS_ROH_SITES_PL) &&
           register_overload(connection, DUCKHTS_ROH_SITES_GT) &&
           register_overload(connection, DUCKHTS_ROH_SITES_COUNTS) &&
           register_packed_decode(connection, false) && register_packed_decode(connection, true) &&
           register_valid_args(connection) &&
           register_duckhts_roh_sites(connection);
}
