/* Test the shared cache against independently retained FASTA bytes. Include the
 * implementation to inspect its retained-storage bound without a public debug API. */
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <errno.h>
#include <inttypes.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "hts_io_tuning.h"

static atomic_uint live_handles;
static atomic_ulong fetch_calls;
static unsigned tuning_requests;
static void observed_tuning(faidx_t *fai, const char *path, duckhts_hts_io_profile_t profile) {
    tuning_requests++;
    duckhts_apply_remote_faidx_tuning(fai, path, profile);
}
static faidx_t *observed_load(const char *p, const char *f, const char *g,
                             int flags, enum fai_format_options format) {
    faidx_t *fai = fai_load3_format(p, f, g, flags, format);
    if (fai) atomic_fetch_add(&live_handles, 1);
    return fai;
}
static void observed_destroy(faidx_t *fai) {
    fai_destroy(fai);
    atomic_fetch_sub(&live_handles, 1);
}
static char *observed_fetch(const faidx_t *fai, const char *chrom,
                           hts_pos_t beg, hts_pos_t end, hts_pos_t *length) {
    atomic_fetch_add(&fetch_calls, 1);
    return faidx_fetch_seq64(fai, chrom, beg, end, length);
}
#define fai_load3_format observed_load
#define fai_destroy observed_destroy
#define faidx_fetch_seq64 observed_fetch
#define duckhts_apply_remote_faidx_tuning observed_tuning
#include "../../src/reference_cache.c"
#undef fai_load3_format
#undef fai_destroy
#undef faidx_fetch_seq64
#undef duckhts_apply_remote_faidx_tuning

enum { REFERENCES = 12, BASES = 262144, WORKERS = 8 };
static char paths[REFERENCES][1024];
static char expected[REFERENCES][BASES];
static unsigned trials = 2000;
static uint32_t seed = 171;
static pthread_mutex_t start_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t start_cond = PTHREAD_COND_INITIALIZER;
static unsigned ready;
struct worker { unsigned id; int failed; size_t bases; };

#define CHECK(test) do { if (!(test)) { \
    fprintf(stderr, "seed=%u failed line %d: %s\n", seed, __LINE__, #test); \
    return 0; \
} } while (0)

static uint32_t random_u32(uint32_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    return *state;
}

static int write_references(const char *directory) {
    uint32_t state = seed;
    for (unsigned r = 0; r < REFERENCES; r++) {
        int n = snprintf(paths[r], sizeof(paths[r]), "%s/ref-%u.fa%s",
                         directory, r, r % 2 ? ".gz" : "");
        CHECK(n > 0 && (size_t)n < sizeof(paths[r]));
        BGZF *out = bgzf_open(paths[r], r % 2 ? "w" : "wu");
        CHECK(out != NULL);
        CHECK(bgzf_write(out, ">chr1\n", 6) == 6);
        for (size_t p = 0; p < BASES; p += 64) {
            char line[66];
            for (size_t i = 0; i < 64; i++) {
                char base = "ACGTNRYSWKM"[random_u32(&state) % 11];
                expected[r][p + i] = base;
                line[i] = i % 3 ? base : (char)(base + ('a' - 'A'));
            }
            size_t bytes = r % 3 ? 65u : 66u;
            line[64] = r % 3 ? '\n' : '\r'; line[65] = '\n';
            CHECK(bgzf_write(out, line, bytes) == (ssize_t)bytes);
        }
        static const char tail[] = ">1\nACGT\n>MT\nRYACGT\n";
        CHECK(bgzf_write(out, tail, sizeof(tail) - 1) == sizeof(tail) - 1);
        CHECK(bgzf_close(out) == 0);
        CHECK(fai_build(paths[r]) == 0);
    }
    return 1;
}

static int bounded_fetch_edges(faidx_t *fai) {
    char buffer[16], saved[16];
    memset(buffer, 0xa5, sizeof buffer); memcpy(saved, buffer, sizeof buffer);
    hts_pos_t length;
    size_t required;
    CHECK(faidx_fetch_seq64_into(NULL, "1", 0, 3, buffer, sizeof buffer, &length, &required) == -1);
    CHECK(errno == EINVAL && length == -1 && required == 0u);
    CHECK(faidx_fetch_seq64_into(fai, NULL, 0, 3, buffer, sizeof buffer, &length, &required) == -1);
    CHECK(errno == EINVAL && length == -1 && required == 0u);
    CHECK(faidx_fetch_seq64_into(fai, "1", 0, 3, NULL, 1u, &length, &required) == -1);
    CHECK(errno == EINVAL && length == -1 && required == 0u);
    CHECK(faidx_fetch_seq64_into(fai, "1", 0, 3, buffer, sizeof buffer, NULL, &required) == -1);
    CHECK(errno == EINVAL && required == 0u);
    CHECK(faidx_fetch_seq64_into(fai, "1", 0, 3, buffer, sizeof buffer, &length, NULL) == -1);
    CHECK(errno == EINVAL && length == -1);
    CHECK(faidx_fetch_seq64_into(fai, "absent", 0, 3, buffer, sizeof buffer, &length, &required) == -1);
    CHECK(length == -2 && required == 0u && !memcmp(buffer, saved, sizeof buffer));
    static const struct { hts_pos_t begin, end; const char *bases; } cases[] = {
        {-20, 2, "ACG"}, {3, 20, "T"}, {5, 5, ""}, {3, 1, "C"}, {-3, -1, "A"}, {0, 3, "ACGT"},
        {INT64_MAX, INT64_MAX, ""}, {0, INT64_MAX, "ACGT"}
    };
    for (size_t i = 0u; i < sizeof cases / sizeof cases[0]; i++) {
        memset(buffer, 0xa5, sizeof buffer);
        CHECK(faidx_fetch_seq64_into(fai, "1", cases[i].begin, cases[i].end,
            buffer, sizeof buffer, &length, &required) == 0);
        CHECK(length == (hts_pos_t)strlen(cases[i].bases));
        CHECK(required == (size_t)length + 2u);
        CHECK(!strcmp(buffer, cases[i].bases) && (unsigned char)buffer[required] == 0xa5u);
    }
    return 1;
}

static int bounded_fetch_failures(const char *directory) {
    char fasta[1100], index[1100];
    int n = snprintf(fasta, sizeof fasta, "%s/bounded-failures.fa", directory);
    CHECK(n > 0 && (size_t)n < sizeof fasta);
    n = snprintf(index, sizeof index, "%s/bounded-failures.fai", directory);
    CHECK(n > 0 && (size_t)n < sizeof index);
    FILE *file = fopen(fasta, "wb");
    CHECK(file != NULL);
    CHECK(fputs(">ok\nAcGT\n", file) >= 0 && fclose(file) == 0);
    file = fopen(index, "wb");
    CHECK(file != NULL);
    CHECK(fprintf(file,
        "ok\t4\t4\t4\t5\n"
        "truncated\t20\t4\t4\t5\n"
        "zero_width\t4\t4\t0\t1\n"
        "offset_wrap\t20\t%" PRIu64 "\t4\t5\n"
        "line_product\t%" PRIu64 "\t0\t1\t3\n"
        "signed_end\t%" PRIu64 "\t4\t4\t5\n",
        UINT64_MAX - 1u, UINT64_MAX, UINT64_MAX) > 0);
    CHECK(fclose(file) == 0);
    faidx_t *fai = observed_load(fasta, index, NULL, 0, FAI_FASTA);
    CHECK(fai != NULL);
    static const struct {
        const char *name;
        hts_pos_t begin, end;
        int error;
    } failures[] = {
        {"truncated", 0, 19, 0},
        {"zero_width", 0, 3, EINVAL},
        /* Without checked seek arithmetic, the offset wraps to byte 3 and
         * returns the header newline as if it were a reference base. */
        {"offset_wrap", 4, 4, EOVERFLOW},
        {"line_product", INT64_MAX - 1, INT64_MAX - 1, EOVERFLOW},
        {"signed_end", INT64_MAX, INT64_MAX, EOVERFLOW}
    };
    for (size_t i = 0; i < sizeof failures / sizeof failures[0]; i++) {
        unsigned char buffer[66];
        memset(buffer, 0xa5, sizeof buffer);
        hts_pos_t length = 42;
        size_t required = 42;
        errno = 0;
        CHECK(faidx_fetch_seq64_into(fai, failures[i].name, failures[i].begin,
            failures[i].end, (char *)buffer + 1, sizeof buffer - 2, &length, &required) == -1);
        CHECK(length == -1 && (!failures[i].error || errno == failures[i].error));
        CHECK(buffer[0] == 0xa5u && buffer[sizeof buffer - 1] == 0xa5u);
        if (failures[i].error) {
            for (size_t j = 0; j < sizeof buffer; j++) CHECK(buffer[j] == 0xa5u);
        }
        char *allocated = faidx_fetch_seq64(fai, failures[i].name,
            failures[i].begin, failures[i].end, &length);
        CHECK(allocated == NULL && length == -1);
        CHECK(!failures[i].error || errno == failures[i].error);
        CHECK(faidx_fetch_seq64_into(fai, "ok", 0, 3, (char *)buffer + 1,
            sizeof buffer - 2, &length, &required) == 0);
        CHECK(length == 4 && !strcmp((char *)buffer + 1, "AcGT"));
    }
    observed_destroy(fai);
    CHECK(atomic_load(&live_handles) == 0);
    return 1;
}

static int bounds_and_aliases(void) {
    const char *error = NULL;
    for (unsigned r = 0; r < REFERENCES; r++) {
        duckhts_reference_entry_t *entry =
            duckhts_reference_cache_get(paths[r], NULL, NULL, 0, &error);
        CHECK(entry && !error);
        CHECK(bounded_fetch_edges(entry->fai));
        char *bases = duckhts_reference_fetch(entry, "chr1", 10, 20,
                                            DUCKHTS_REFERENCE_EXACT);
        CHECK(bases && !memcmp(bases, expected[r] + 10, 11) && bases[11] == 0);
        free(bases);
        unsigned long calls = atomic_load(&fetch_calls);
        bases = duckhts_reference_fetch(entry, "chr1", 12, 15,
                                       DUCKHTS_REFERENCE_EXACT);
        CHECK(bases && !memcmp(bases, expected[r] + 12, 4));
        CHECK(atomic_load(&fetch_calls) == calls);
        bases[0] = '!'; /* Caller-owned results must not poison the shared window. */
        free(bases);
        bases = duckhts_reference_fetch(entry, "chr1", 12, 15,
                                       DUCKHTS_REFERENCE_EXACT);
        CHECK(bases && !memcmp(bases, expected[r] + 12, 4));
        free(bases);
        for (size_t length = DUCKHTS_REFERENCE_WINDOW_BASES;
             length <= DUCKHTS_REFERENCE_WINDOW_BASES + 1; length++) {
            bases = duckhts_reference_fetch(entry, "chr1", 4001, 4000 + length,
                                           DUCKHTS_REFERENCE_EXACT);
            CHECK(bases && strlen(bases) == length &&
                  !memcmp(bases, expected[r] + 4001, length));
            free(bases);
        }
        bases = duckhts_reference_fetch(entry, "chr1", 123, 70122,
                                       DUCKHTS_REFERENCE_EXACT);
        CHECK(bases && strlen(bases) == 70000 &&
              !memcmp(bases, expected[r] + 123, 70000));
        free(bases);
        reference_cache_t *cache = pthread_getspecific(reference_key);
        CHECK(cache->count <= DUCKHTS_REFERENCE_CACHE_ENTRIES);
        CHECK(atomic_load(&live_handles) <= DUCKHTS_REFERENCE_CACHE_ENTRIES);
        size_t retained = 0;
        for (duckhts_reference_entry_t *e = cache->head; e; e = e->next) {
            if (!e->window_seq) continue;
            size_t length = (size_t)(e->window_end - e->window_beg + 1);
            CHECK(length <= DUCKHTS_REFERENCE_WINDOW_BASES);
            retained += length;
        }
        CHECK(retained <= DUCKHTS_REFERENCE_CACHE_ENTRIES * DUCKHTS_REFERENCE_WINDOW_BASES);
        CHECK(duckhts_reference_fetch(entry, "chr1", BASES - 3, BASES + 5,
                                     DUCKHTS_REFERENCE_EXACT) == NULL);
        bases = duckhts_reference_fetch(entry, "chr1", BASES - 3, BASES + 5,
                                       DUCKHTS_REFERENCE_CLIP_END);
        CHECK(bases && strlen(bases) == 3 &&
              !memcmp(bases, expected[r] + BASES - 3, 3));
        free(bases);
        /* EXACT retries chr1 when the literal contig 1 cannot supply the span;
         * normalization's CLIP_END keeps the first matching contig authoritative. */
        bases = duckhts_reference_fetch(entry, "1", 12, 15, DUCKHTS_REFERENCE_EXACT);
        CHECK(bases && !memcmp(bases, expected[r] + 12, 4));
        free(bases);
        CHECK(duckhts_reference_fetch(entry, "1", 12, 15,
                                     DUCKHTS_REFERENCE_CLIP_END) == NULL);
        bases = duckhts_reference_fetch(entry, "chrM", 0, 5, DUCKHTS_REFERENCE_EXACT);
        CHECK(bases && !strcmp(bases, "RYACGT"));
        free(bases);
        CHECK(duckhts_reference_fetch(entry, "chr1", -1, 4,
                                     DUCKHTS_REFERENCE_EXACT) == NULL);
    }
    /* Reacquire an evicted file and verify file identity, not pointer identity. */
    duckhts_reference_entry_t *entry =
        duckhts_reference_cache_get(paths[0], NULL, NULL, 0, &error);
    CHECK(entry != NULL);
    char *bases = duckhts_reference_fetch(entry, "chr1", 0, 4095, DUCKHTS_REFERENCE_EXACT);
    CHECK(bases && !memcmp(bases, expected[0], 4096));
    free(bases);
    CHECK(duckhts_reference_cache_get(paths[0], "", "", 0, &error) == entry);
    CHECK(tuning_requests == 0);
    duckhts_reference_entry_t *tuned =
        duckhts_reference_cache_get(paths[0], NULL, NULL, 1, &error);
    CHECK(tuned && tuned != entry && tuned->indexed_remote_tuned);
    CHECK(tuning_requests == 1);
    CHECK(duckhts_reference_cache_get(paths[0], NULL, NULL, 0, &error) == entry);
    CHECK(duckhts_reference_cache_get(paths[0], NULL, NULL, 1, &error) == tuned);
    CHECK(!entry->indexed_remote_tuned);
    CHECK(tuning_requests == 1);
    char index[1100];
    snprintf(index, sizeof(index), "%s.fai", paths[0]);
    duckhts_reference_entry_t *explicit_index =
        duckhts_reference_cache_get(paths[0], index, NULL, 0, &error);
    CHECK(explicit_index && explicit_index != entry);
    CHECK(duckhts_reference_cache_get(paths[0], index, NULL, 0, &error) == explicit_index);
    char gzi[1100];
    snprintf(index, sizeof(index), "%s.fai", paths[1]);
    snprintf(gzi, sizeof(gzi), "%s.gzi", paths[1]);
    explicit_index = duckhts_reference_cache_get(paths[1], index, gzi, 0, &error);
    CHECK(explicit_index && !error);
    bases = duckhts_reference_fetch(explicit_index, "chr1", 0, 4095,
                                   DUCKHTS_REFERENCE_EXACT);
    CHECK(bases && !memcmp(bases, expected[1], 4096));
    free(bases);
    reference_cache_t *cache = pthread_getspecific(reference_key);
    CHECK(pthread_setspecific(reference_key, NULL) == 0);
    reference_cache_destroy(cache);
    CHECK(atomic_load(&live_handles) == 0);
    return 1;
}

static int bounded_fetch_matches(faidx_t *fai, unsigned reference,
    size_t position, size_t length) {
    unsigned char buffer[4102], saved[4102];
    memset(buffer, 0xa5, sizeof buffer); memcpy(saved, buffer, sizeof saved);
    hts_pos_t got_length = -1;
    size_t required = 0u;
    errno = 0;
    CHECK(faidx_fetch_seq64_into(fai, "chr1", position, position + length - 1u,
        NULL, 0u, &got_length, &required) == -1);
    CHECK(errno == ENOSPC && got_length == (hts_pos_t)length);
    CHECK(required == length + (reference % 3 ? 2u : 3u));
    CHECK(faidx_fetch_seq64_into(fai, "chr1", position, position + length - 1u,
        (char *)buffer + 1u, required - 1u, &got_length, &required) == -1);
    CHECK(errno == ENOSPC && !memcmp(buffer, saved, sizeof buffer));
    CHECK(faidx_fetch_seq64_into(fai, "chr1", position, position + length - 1u,
        (char *)buffer + 1u, required, &got_length, &required) == 0);
    CHECK(got_length == (hts_pos_t)length && buffer[length + 1u] == 0u);
    CHECK(buffer[0] == 0xa5u && buffer[required + 1u] == 0xa5u);
    for (size_t i = 0u; i < length; i++) {
        unsigned char raw = (unsigned char)expected[reference][position + i];
        if ((position + i) % 64u % 3u == 0u) raw += 'a' - 'A';
        CHECK(buffer[i + 1u] == raw);
    }
    return 1;
}

static void *random_worker(void *pointer) {
    struct worker *w = pointer;
    uint32_t state = seed + w->id;
    pthread_mutex_lock(&start_mutex);
    ready++;
    pthread_cond_broadcast(&start_cond);
    while (ready < WORKERS) pthread_cond_wait(&start_cond, &start_mutex);
    pthread_mutex_unlock(&start_mutex);
    /* All workers acquire the same initially unindexed plain FASTA. */
    const char *initial_error = NULL;
    if (!duckhts_reference_cache_get(paths[0], NULL, NULL, 0, &initial_error)) {
        fprintf(stderr, "index creation: %s\n", initial_error);
        w->failed = 1;
        return NULL;
    }
    for (unsigned i = 0; i < trials; i++) {
        unsigned r = (i / 32 + w->id) % REFERENCES;
        size_t length = 1 + random_u32(&state) % 4096;
        size_t position = random_u32(&state) % (BASES - length);
        const char *error = NULL;
        duckhts_reference_entry_t *entry =
            duckhts_reference_cache_get(paths[r], NULL, NULL, 0, &error);
        char *bases = entry ? duckhts_reference_fetch(entry, "chr1", position,
            position + length - 1, DUCKHTS_REFERENCE_EXACT) : NULL;
        if (!bases || strlen(bases) != length ||
            memcmp(bases, expected[r] + position, length)) {
            fprintf(stderr, "seed=%u worker=%u trial=%u ref=%u pos=%zu len=%zu error=%s\n",
                    seed, w->id, i, r, position, length, error ? error : "");
            w->failed = 1;
            free(bases);
            break;
        }
        free(bases);
        if (!bounded_fetch_matches(entry->fai, r, position, length)) {
            w->failed = 1;
            break;
        }
        w->bases += length;
    }
    return NULL; /* pthread destructor must close every cached handle. */
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) return 2;
    if (argc > 2) seed = (uint32_t)strtoul(argv[2], NULL, 10);
    if (argc > 3) trials = (unsigned)strtoul(argv[3], NULL, 10);
    if (!seed || !trials) return 2;
    if (!write_references(argv[1]) || !bounded_fetch_failures(argv[1]) || !bounds_and_aliases()) return 1;
    char index[1100];
    snprintf(index, sizeof(index), "%s.fai", paths[0]);
    if (remove(index)) return 1;
    pthread_t threads[WORKERS];
    struct worker workers[WORKERS] = {{0}};
    for (unsigned i = 0; i < WORKERS; i++) {
        workers[i].id = i;
        if (pthread_create(&threads[i], NULL, random_worker, &workers[i])) return 2;
    }
    size_t bases = 0;
    for (unsigned i = 0; i < WORKERS; i++) {
        if (pthread_join(threads[i], NULL) || workers[i].failed) return 1;
        bases += workers[i].bases;
    }
    if (atomic_load(&live_handles) != 0) return 1;
    printf("reference cache: seed=%u workers=%u requests=%u bases=%zu fetches=%lu live_handles=0\n",
           seed, WORKERS, WORKERS * trials, bases, atomic_load(&fetch_calls));
    return 0;
}
