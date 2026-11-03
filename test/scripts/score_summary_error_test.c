/* Exercise the production summary loaders with deterministic allocation failures.
 * htslib owns its internal allocations; the wrappers cover this translation unit
 * and separately inject its documented FORMAT allocation-error return (-4).
 */
#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <htslib/vcf.h>
#include "duckdb_extension.h"

static void *live[4096];
static size_t nlive;
static int allocation_count;
static int fail_at;
static const char *fail_format;
static int fail_record;

static int allocation_fails(void) {
    return ++allocation_count == fail_at;
}

static void remember(void *ptr) {
    if (ptr) {
        assert(nlive < sizeof(live) / sizeof(live[0]));
        live[nlive++] = ptr;
    }
}

static void forget(void *ptr) {
    for (size_t i = 0; i < nlive; i++) {
        if (live[i] == ptr) {
            live[i] = live[--nlive];
            return;
        }
    }
}

static void *test_malloc(size_t size) {
    if (allocation_fails()) return NULL;
    void *ptr = malloc(size);
    remember(ptr);
    return ptr;
}

static void *test_calloc(size_t count, size_t size) {
    if (allocation_fails()) return NULL;
    void *ptr = calloc(count, size);
    remember(ptr);
    return ptr;
}

static void *test_realloc(void *ptr, size_t size) {
    if (allocation_fails()) return NULL;
    /* Restore tracking if realloc fails without consuming the old pointer. */
    forget(ptr);
    void *next = realloc(ptr, size);
    remember(next ? next : ptr);
    return next;
}

static void test_free(void *ptr) {
    forget(ptr);
    free(ptr);
}

static bcf1_t *test_bcf_init(void) {
    return fail_record ? NULL : bcf_init();
}

static int test_get_format(const bcf_hdr_t *hdr, bcf1_t *rec, const char *tag,
                           void **dst, int *capacity, int type) {
    if (fail_format && strcmp(tag, fail_format) == 0) return -4;
    return bcf_get_format_values(hdr, rec, tag, dst, capacity, type);
}

#define malloc test_malloc
#define calloc test_calloc
#define realloc test_realloc
#define free test_free
#define bcf_init test_bcf_init
#define bcf_get_format_values test_get_format
#include "../../src/score_udf.c"
#undef malloc
#undef calloc
#undef realloc
#undef free
#undef bcf_init
#undef bcf_get_format_values

DUCKDB_EXTENSION_GLOBAL

int main(int argc, char **argv) {
    if (argc < 3) return 2;
    const char *mode = argv[1];
    const char *path = argv[2];
    int inject = argc > 3 ? atoi(argv[3]) : 0;
    int gwas = strcmp(mode, "gwas") == 0;
    bcf_hdr_t *hdr = bcf_hdr_init("r");
    assert(hdr);
    assert(bcf_hdr_append(hdr, "##contig=<ID=1,length=1000000>") == 0);
    assert(bcf_hdr_sync(hdr) == 0);
    duckdb_malloc = test_malloc;
    duckdb_free = test_free;
    fail_at = inject;
    if (argc > 4) {
        if (strcmp(argv[4], "record") == 0) fail_record = 1;
        else fail_format = argv[4];
    }

    char err[512] = {0};
    int rc;
    int markers = 0;
    score_bind_t bind = {0};
    bind.columns_preset = "PLINK";
    bind.use_variant_id = strcmp(mode, "rsid") == 0;
    double threshold = 0.0;
    bind.q_thr_lp = &threshold;
    if (gwas) {
        score_summary_t **summaries = score_load_summary_vcf(&bind, hdr, path, 2, err, sizeof(err));
        rc = summaries ? 0 : -1;
        if (summaries) {
            for (int i = 0; i < 2; i++) {
                markers += summaries[i]->n_markers;
                score_destroy_summary(summaries[i]);
            }
            test_free(summaries);
        }
    } else {
        score_summary_t *summary = score_load_summary(&bind, hdr, path, err, sizeof(err));
        rc = summary ? 0 : -1;
        if (summary) markers = summary->n_markers;
        score_destroy_summary(summary);
    }
    bcf_hdr_destroy(hdr);
    printf("%s allocations=%d rc=%d markers=%d live=%zu error=%s\n",
           mode, allocation_count, rc, markers, nlive, err);
    if (nlive != 0) return 1;
    if (inject || fail_record || fail_format) {
        if (rc == 0 || !strstr(err, "bcftools_score:") ||
            !strstr(err, "out of memory") || !strstr(err, path)) return 1;
    } else {
        if (rc != 0 || markers != (gwas ? 5 : 3)) return 1;
    }
    return 0;
}
