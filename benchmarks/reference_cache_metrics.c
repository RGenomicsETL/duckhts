/* Linux benchmark-only interposition. No production counters or DuckDB API. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <htslib/faidx.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>

static pthread_once_t metrics_once = PTHREAD_ONCE_INIT;
static char *(*real_fetch)(const faidx_t *, const char *, hts_pos_t, hts_pos_t, hts_pos_t *);
static faidx_t *(*real_load)(const char *, const char *, const char *, int, enum fai_format_options);
static faidx_t *(*real_load3)(const char *, const char *, const char *, int);
static void (*real_destroy)(faidx_t *);
static atomic_ulong fetches, fetched_bases, loads, destroys;

static void metrics_init(void) {
    const char *path = getenv("REFERENCE_CACHE_EXTENSION");
    void *handle = path ? dlopen(path, RTLD_NOW | RTLD_NOLOAD) : NULL;
    if (!handle) {
        fprintf(stderr, "reference metrics: extension is not loaded\n");
        exit(2);
    }
    *(void **)(&real_fetch) = dlsym(handle, "faidx_fetch_seq64");
    *(void **)(&real_load) = dlsym(handle, "fai_load3_format");
    *(void **)(&real_load3) = dlsym(handle, "fai_load3");
    *(void **)(&real_destroy) = dlsym(handle, "fai_destroy");
    if (!real_fetch || !real_load || !real_load3 || !real_destroy || real_fetch == faidx_fetch_seq64) {
        fprintf(stderr, "reference metrics: cannot resolve htslib functions\n");
        exit(2);
    }
    dlclose(handle);
}

char *faidx_fetch_seq64(const faidx_t *fai, const char *chrom,
                      hts_pos_t begin, hts_pos_t end, hts_pos_t *length) {
    pthread_once(&metrics_once, metrics_init);
    char *bases = real_fetch(fai, chrom, begin, end, length);
    atomic_fetch_add(&fetches, 1);
    if (bases && *length > 0) atomic_fetch_add(&fetched_bases, (unsigned long)*length);
    return bases;
}

faidx_t *fai_load3_format(const char *path, const char *fai, const char *gzi,
                        int flags, enum fai_format_options format) {
    pthread_once(&metrics_once, metrics_init);
    faidx_t *result = real_load(path, fai, gzi, flags, format);
    if (result) atomic_fetch_add(&loads, 1);
    return result;
}

faidx_t *fai_load3(const char *path, const char *fai, const char *gzi, int flags) {
    pthread_once(&metrics_once, metrics_init);
    faidx_t *result = real_load3(path, fai, gzi, flags);
    if (result) atomic_fetch_add(&loads, 1);
    return result;
}

void fai_destroy(faidx_t *fai) {
    pthread_once(&metrics_once, metrics_init);
    real_destroy(fai);
    atomic_fetch_add(&destroys, 1);
}

/* Called through R's .C only between completed queries. */
void reference_cache_metrics(double *values) {
    values[0] = (double)atomic_load(&fetches);
    values[1] = (double)atomic_load(&fetched_bases);
    values[2] = (double)atomic_load(&loads);
    values[3] = (double)atomic_load(&destroys);
}
