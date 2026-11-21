/*
 * The indexed FASTA of a scalar function. See reference_access.h for the
 * contract.
 */
#include "reference_access.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hts_file_access.h"

DUCKDB_EXTENSION_EXTERN

/* One reference read through the connection. The paths are its key. */
typedef struct {
    char *fasta_path;
    char *fai_path; /* NULL: the FASTA path with .fai appended */
    char *gzi_path; /* NULL: the FASTA path with .gzi appended */
    int indexed_remote_tuned;
    duckhts_reference_entry_t *entry;
} routed_reference_t;

struct duckhts_reference_access {
    duckhts_file_access_t *access;
    /* Used only when access reads through DuckDB; most recently used first. */
    routed_reference_t routed[DUCKHTS_REFERENCE_CACHE_ENTRIES];
    unsigned n_routed;
};

static void routed_reference_close(routed_reference_t *reference) {
    duckhts_reference_entry_destroy(reference->entry);
    free(reference->fasta_path);
    free(reference->fai_path);
    free(reference->gzi_path);
    memset(reference, 0, sizeof(*reference));
}

/* Closes the routed handles before the route they were opened through. */
static void reference_access_destroy(void *pointer) {
    duckhts_reference_access_t *references = pointer;
    if (!references) {
        return;
    }
    for (unsigned i = 0; i < references->n_routed; i++) {
        routed_reference_close(&references->routed[i]);
    }
    duckhts_file_access_destroy(references->access);
    free(references);
}

void duckhts_reference_access_scalar_init(duckdb_init_info info, const char *function_name) {
    duckhts_reference_access_t *references = calloc(1, sizeof(*references));
    if (!references) {
        char message[160];
        snprintf(message, sizeof(message), "%s: out of memory", function_name);
        duckdb_scalar_function_init_set_error(info, message);
        return;
    }
    if (!duckhts_file_access_scalar_init(info, function_name, &references->access)) {
        free(references);
        return;
    }
    duckdb_scalar_function_init_set_state(info, references, reference_access_destroy);
}

static bool same_path(const char *a, const char *b) {
    return a && b ? strcmp(a, b) == 0 : a == b;
}

static char *copy_path(const char *path) {
    const size_t length = strlen(path);
    char *copy = malloc(length + 1);
    if (copy) {
        memcpy(copy, path, length + 1);
    }
    return copy;
}

/*
 * DuckDB's reason for a load of these paths that just failed: the reason for
 * the index, the FASTA or the .gzi, whichever the route's most recent failure
 * concerns.
 */
static const char *routed_load_reason(const duckhts_file_access_t *access, const char *fasta_path,
                                      const char *fai_path, const char *gzi_path) {
    static const char *const suffixes[] = {".fai", ".gzi"};
    const char *const explicit_paths[] = {fai_path, gzi_path};
    const size_t fasta_length = strlen(fasta_path);
    const char *reason = duckhts_file_access_reason(access, fasta_path);
    if (reason[0] || fasta_length > SIZE_MAX - sizeof(".fai")) {
        return reason;
    }
    char *beside = malloc(fasta_length + sizeof(".fai"));
    if (!beside) {
        return reason;
    }
    for (size_t i = 0; i < 2 && !reason[0]; i++) {
        memcpy(beside, fasta_path, fasta_length);
        memcpy(beside + fasta_length, suffixes[i], sizeof(".fai"));
        reason = duckhts_file_access_reason(access, explicit_paths[i] ? explicit_paths[i] : beside);
    }
    free(beside);
    return reason;
}

duckhts_reference_entry_t *duckhts_reference_access_get(
    duckhts_reference_access_t *references, const char *fasta_path, const char *fai_path,
    const char *gzi_path, int tune_indexed_remote_io, const char **error, const char **reason) {
    *error = NULL;
    *reason = "";
    if (!references) {
        *error = "DuckDB provided no function state";
        return NULL;
    }
    if (!duckhts_file_access_uses_duckdb(references->access)) {
        return duckhts_reference_cache_get(fasta_path, fai_path, gzi_path, tune_indexed_remote_io,
                                           error);
    }
    if (!fasta_path || !*fasta_path) {
        *error = "failed to load FASTA index";
        return NULL;
    }
    if (fai_path && !*fai_path) {
        fai_path = NULL;
    }
    if (gzi_path && !*gzi_path) {
        gzi_path = NULL;
    }
    const int tuned = tune_indexed_remote_io != 0;

    routed_reference_t *routed = references->routed;
    for (unsigned i = 0; i < references->n_routed; i++) {
        if (same_path(routed[i].fasta_path, fasta_path) && same_path(routed[i].fai_path, fai_path) &&
            same_path(routed[i].gzi_path, gzi_path) && routed[i].indexed_remote_tuned == tuned) {
            const routed_reference_t found = routed[i];
            memmove(&routed[1], &routed[0], i * sizeof(routed[0]));
            routed[0] = found;
            return found.entry;
        }
    }

    routed_reference_t opened = {0};
    opened.indexed_remote_tuned = tuned;
    opened.fasta_path = copy_path(fasta_path);
    opened.fai_path = fai_path ? copy_path(fai_path) : NULL;
    opened.gzi_path = gzi_path ? copy_path(gzi_path) : NULL;
    if (!opened.fasta_path || (fai_path && !opened.fai_path) || (gzi_path && !opened.gzi_path)) {
        routed_reference_close(&opened);
        *error = "out of memory";
        return NULL;
    }
    /* Close the least recently used handle before opening another one. */
    if (references->n_routed == DUCKHTS_REFERENCE_CACHE_ENTRIES) {
        routed_reference_close(&routed[--references->n_routed]);
    }
    faidx_t *fai = duckhts_file_access_fai_load(
        references->access, fasta_path, fai_path, gzi_path, 0, FAI_FASTA,
        tuned ? DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION : DUCKHTS_HTS_IO_PROFILE_METADATA);
    if (!fai) {
        routed_reference_close(&opened);
        *error = duckhts_reference_load_error(fasta_path, fai_path, gzi_path);
        *reason = routed_load_reason(references->access, fasta_path, fai_path, gzi_path);
        return NULL;
    }
    opened.entry = duckhts_reference_entry_wrap(fai);
    if (!opened.entry) {
        fai_destroy(fai);
        routed_reference_close(&opened);
        *error = "out of memory";
        return NULL;
    }
    memmove(&routed[1], &routed[0], references->n_routed * sizeof(routed[0]));
    routed[0] = opened;
    references->n_routed++;
    return opened.entry;
}
