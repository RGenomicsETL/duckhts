#include "reference_cache.h"
#include "hts_io_tuning.h"

#include <ctype.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#define REFERENCE_FETCH_PAD 4096

struct duckhts_reference_entry {
    char *fasta_path;
    char *fai_path;
    char *gzi_path;
    faidx_t *fai;
    int indexed_remote_tuned;
    char *window_contig;
    char *window_seq;
    hts_pos_t window_beg;
    hts_pos_t window_end;
    struct duckhts_reference_entry *next;
};

typedef struct {
    duckhts_reference_entry_t *head;
    unsigned count;
} reference_cache_t;

static pthread_key_t reference_key;
static pthread_once_t reference_once = PTHREAD_ONCE_INIT;
static int reference_key_error;
/* FAI_CREATE may write a common .fai/.gzi file. Only handle acquisition is
 * serialized; every fetch uses a handle owned by the calling thread. */
static pthread_mutex_t reference_index_mutex = PTHREAD_MUTEX_INITIALIZER;

static char *reference_copy(const char *text, size_t length) {
    if (length == SIZE_MAX) return NULL;
    char *copy = malloc(length + 1);
    if (!copy) return NULL;
    memcpy(copy, text, length);
    copy[length] = '\0';
    return copy;
}

static char *reference_dup(const char *text) {
    return text ? reference_copy(text, strlen(text)) : NULL;
}

static int reference_path_equal(const char *a, const char *b) {
    return a && b ? strcmp(a, b) == 0 : a == b;
}

static void reference_entry_destroy(duckhts_reference_entry_t *entry) {
    if (entry->fai) fai_destroy(entry->fai);
    free(entry->fasta_path);
    free(entry->fai_path);
    free(entry->gzi_path);
    free(entry->window_contig);
    free(entry->window_seq);
    free(entry);
}

static void reference_cache_destroy(void *pointer) {
    reference_cache_t *cache = pointer;
    if (!cache) return;
    while (cache->head) {
        duckhts_reference_entry_t *entry = cache->head;
        cache->head = entry->next;
        reference_entry_destroy(entry);
    }
    free(cache);
}

static void reference_key_init(void) {
    reference_key_error = pthread_key_create(&reference_key, reference_cache_destroy);
}

duckhts_reference_entry_t *duckhts_reference_cache_get(
    const char *fasta_path, const char *fai_path, const char *gzi_path,
    int tune_indexed_remote_io, const char **error) {
    reference_cache_t *cache;
    duckhts_reference_entry_t *entry, **link;
    if (!error) return NULL;
    *error = NULL;
    if (!fasta_path || !*fasta_path) {
        *error = "failed to load FASTA index";
        return NULL;
    }
    if (fai_path && !*fai_path) fai_path = NULL;
    if (gzi_path && !*gzi_path) gzi_path = NULL;
    if (pthread_once(&reference_once, reference_key_init) || reference_key_error) {
        *error = "failed to initialize thread-local FASTA cache";
        return NULL;
    }
    cache = pthread_getspecific(reference_key);
    if (!cache) {
        cache = calloc(1, sizeof(*cache));
        if (!cache) goto oom;
        if (pthread_setspecific(reference_key, cache)) {
            free(cache);
            *error = "failed to attach thread-local FASTA cache";
            return NULL;
        }
    }
    for (link = &cache->head; (entry = *link) != NULL; link = &entry->next) {
        if (reference_path_equal(entry->fasta_path, fasta_path) &&
            reference_path_equal(entry->fai_path, fai_path) &&
            reference_path_equal(entry->gzi_path, gzi_path) &&
            entry->indexed_remote_tuned == (tune_indexed_remote_io != 0)) {
            *link = entry->next;
            entry->next = cache->head;
            cache->head = entry;
            return entry;
        }
    }
    entry = calloc(1, sizeof(*entry));
    if (!entry) goto oom;
    entry->fasta_path = reference_dup(fasta_path);
    entry->fai_path = reference_dup(fai_path);
    entry->gzi_path = reference_dup(gzi_path);
    if (!entry->fasta_path || (fai_path && !entry->fai_path) ||
        (gzi_path && !entry->gzi_path)) {
        reference_entry_destroy(entry);
        goto oom;
    }
    /* Evict before opening another handle, including on a failed acquisition. */
    if (cache->count == DUCKHTS_REFERENCE_CACHE_ENTRIES) {
        link = &cache->head;
        while ((*link)->next) link = &(*link)->next;
        reference_entry_destroy(*link);
        *link = NULL;
        cache->count--;
    }
    pthread_mutex_lock(&reference_index_mutex);
    entry->fai = fai_load3_format(fasta_path, fai_path, gzi_path, FAI_CREATE, FAI_FASTA);
    pthread_mutex_unlock(&reference_index_mutex);
    if (!entry->fai) {
        reference_entry_destroy(entry);
        *error = "failed to load FASTA index";
        return NULL;
    }
    entry->indexed_remote_tuned = (tune_indexed_remote_io != 0);
    if (entry->indexed_remote_tuned) {
        duckhts_apply_remote_faidx_tuning(entry->fai, fasta_path,
                                        DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    }
    entry->next = cache->head;
    cache->head = entry;
    cache->count++;
    return entry;
oom:
    *error = "out of memory";
    return NULL;
}

void duckhts_reference_uppercase(char *bases) {
    if (!bases) return;
    for (; *bases; bases++) *bases = (char)toupper((unsigned char)*bases);
}

static char *reference_fetch_name(duckhts_reference_entry_t *entry,
    const char *name, hts_pos_t start0, hts_pos_t end0,
    duckhts_reference_end_policy_t end_policy) {
    hts_pos_t contig_length = faidx_seq_len64(entry->fai, name);
    if (contig_length <= 0 || start0 >= contig_length) return NULL;
    if (end0 >= contig_length) {
        if (end_policy == DUCKHTS_REFERENCE_EXACT) return NULL;
        end0 = contig_length - 1;
    }
    uint64_t length = (uint64_t)(end0 - start0) + 1;
    if (length >= SIZE_MAX) return NULL;
    if (entry->window_seq && strcmp(entry->window_contig, name) == 0 &&
        start0 >= entry->window_beg && end0 <= entry->window_end) {
        return reference_copy(entry->window_seq + (size_t)(start0 - entry->window_beg),
                              (size_t)length);
    }
    hts_pos_t fetch_beg = start0, fetch_end = end0, found = 0;
    if (length <= DUCKHTS_REFERENCE_WINDOW_BASES) {
        hts_pos_t pad = DUCKHTS_REFERENCE_WINDOW_BASES - (hts_pos_t)length;
        if (pad > REFERENCE_FETCH_PAD) pad = REFERENCE_FETCH_PAD;
        fetch_beg = start0 > pad ? start0 - pad : 0;
        hts_pos_t available = contig_length - fetch_beg;
        if (available > DUCKHTS_REFERENCE_WINDOW_BASES)
            available = DUCKHTS_REFERENCE_WINDOW_BASES;
        fetch_end = fetch_beg + available - 1;
    }
    char *bases = faidx_fetch_seq64(entry->fai, name, fetch_beg, fetch_end, &found);
    if (!bases || found != fetch_end - fetch_beg + 1) {
        free(bases);
        return NULL;
    }
    duckhts_reference_uppercase(bases);
    if (length > DUCKHTS_REFERENCE_WINDOW_BASES) return bases;
    char *contig = reference_dup(name);
    if (!contig) {
        free(bases);
        return NULL;
    }
    free(entry->window_contig);
    free(entry->window_seq);
    entry->window_contig = contig;
    entry->window_seq = bases;
    entry->window_beg = fetch_beg;
    entry->window_end = fetch_end;
    return reference_copy(bases + (size_t)(start0 - fetch_beg), (size_t)length);
}

char *duckhts_reference_fetch(duckhts_reference_entry_t *entry,
    const char *chrom, hts_pos_t start0, hts_pos_t end0,
    duckhts_reference_end_policy_t end_policy) {
    static const char *const mitochondrial[] = {"MT", "chrM", "M"};
    const char *aliases[5];
    char *prefixed = NULL, *bases = NULL;
    size_t count = 0;
    if (!entry || !chrom || !*chrom || start0 < 0 || end0 < start0 ||
        (end_policy != DUCKHTS_REFERENCE_EXACT &&
         end_policy != DUCKHTS_REFERENCE_CLIP_END)) return NULL;
    if (entry->window_seq && strcmp(entry->window_contig, chrom) == 0 &&
        start0 >= entry->window_beg && end0 <= entry->window_end) {
        return reference_copy(entry->window_seq + (size_t)(start0 - entry->window_beg),
                              (size_t)(end0 - start0 + 1));
    }
    aliases[count++] = chrom;
    if (strncasecmp(chrom, "chr", 3) == 0) {
        if (chrom[3]) aliases[count++] = chrom + 3;
    } else {
        size_t length = strlen(chrom);
        if (length > SIZE_MAX - 4) return NULL;
        prefixed = malloc(length + 4);
        if (!prefixed) return NULL;
        memcpy(prefixed, "chr", 3);
        memcpy(prefixed + 3, chrom, length + 1);
        aliases[count++] = prefixed;
    }
    if (strcasecmp(chrom, "M") == 0 || strcasecmp(chrom, "MT") == 0 ||
        strcasecmp(chrom, "chrM") == 0) {
        for (size_t i = 0; i < 3; i++) aliases[count++] = mitochondrial[i];
    }
    for (size_t i = 0; i < count; i++) {
        if (faidx_has_seq(entry->fai, aliases[i]) <= 0) continue;
        bases = reference_fetch_name(entry, aliases[i], start0, end0, end_policy);
        if (bases || end_policy == DUCKHTS_REFERENCE_CLIP_END) break;
    }
    free(prefixed);
    return bases;
}
