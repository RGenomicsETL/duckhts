/* Test the shared cache against independently retained FASTA bytes. Include the
 * implementation to inspect its retained-storage bound without a public debug API. */
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
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
            char line[65];
            for (size_t i = 0; i < 64; i++) {
                char base = "ACGTNRYSWKM"[random_u32(&state) % 11];
                expected[r][p + i] = base;
                line[i] = i % 3 ? base : (char)(base + ('a' - 'A'));
            }
            line[64] = '\n';
            CHECK(bgzf_write(out, line, sizeof(line)) == sizeof(line));
        }
        static const char tail[] = ">1\nACGT\n>MT\nRYACGT\n";
        CHECK(bgzf_write(out, tail, sizeof(tail) - 1) == sizeof(tail) - 1);
        CHECK(bgzf_close(out) == 0);
        CHECK(fai_build(paths[r]) == 0);
    }
    return 1;
}

static int bounds_and_aliases(void) {
    const char *error = NULL;
    for (unsigned r = 0; r < REFERENCES; r++) {
        duckhts_reference_entry_t *entry =
            duckhts_reference_cache_get(paths[r], NULL, NULL, 0, &error);
        CHECK(entry && !error);
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
        w->bases += length;
    }
    return NULL; /* pthread destructor must close every cached handle. */
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) return 2;
    if (argc > 2) seed = (uint32_t)strtoul(argv[2], NULL, 10);
    if (argc > 3) trials = (unsigned)strtoul(argv[3], NULL, 10);
    if (!seed || !trials) return 2;
    if (!write_references(argv[1]) || !bounds_and_aliases()) return 1;
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
