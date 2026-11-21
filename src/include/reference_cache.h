#ifndef DUCKHTS_REFERENCE_CACHE_H
#define DUCKHTS_REFERENCE_CACHE_H

#include <htslib/faidx.h>
#include <htslib/hts.h>

#define DUCKHTS_REFERENCE_CACHE_ENTRIES 8u
#define DUCKHTS_REFERENCE_WINDOW_BASES 65536u

typedef struct duckhts_reference_entry duckhts_reference_entry_t;

/* Thread-owned handle and window, borrowed until the next cache_get on this
 * thread. Callers must not retain this pointer across scalar invocations.
 * Keys include the FASTA, optional FAI/GZI paths, and transport-tuning choice.
 * At most eight entries and eight 64-KiB sequence windows are retained per
 * thread; thread exit destroys them. htslib index/transport storage and
 * caller-owned fetch results are separate.
 * tune_indexed_remote_io requests the existing indexed remote-I/O tuning;
 * false leaves transport defaults unchanged. Different tuning choices use
 * different entries; a cached handle is never retuned for another caller.
 * On failure, *error receives a static diagnostic. */
duckhts_reference_entry_t *duckhts_reference_cache_get(
    const char *fasta_path, const char *fai_path, const char *gzi_path,
    int tune_indexed_remote_io, const char **error);

/* An entry outside the thread cache, for a handle that must not outlive its
 * opener: a faidx_t read through one connection's file access. Takes ownership
 * of fai on success; on NULL (out of memory) fai stays with the caller. The
 * caller owns the entry, uses it from one thread at a time with
 * duckhts_reference_fetch, and releases it with duckhts_reference_entry_destroy,
 * which closes fai and frees the entry's one sequence window. cache_get never
 * returns such an entry and must not be asked to release one. */
duckhts_reference_entry_t *duckhts_reference_entry_wrap(faidx_t *fai);
void duckhts_reference_entry_destroy(duckhts_reference_entry_t *entry);

/* The diagnostic cache_get reports when the index of these paths does not
 * load. A NULL or empty fai_path or gzi_path names the FASTA path with .fai or
 * .gzi appended. Borrowed, per-thread storage, valid until the next failed
 * cache_get or the next call on this thread. */
const char *duckhts_reference_load_error(
    const char *fasta_path, const char *fai_path, const char *gzi_path);

typedef enum {
    DUCKHTS_REFERENCE_EXACT,       /* require the complete requested interval */
    DUCKHTS_REFERENCE_CLIP_END     /* first matching contig, clip at its end */
} duckhts_reference_end_policy_t;

/* Inclusive, zero-based coordinates. Resolve the existing chr/MT aliases in
 * priority order. Return malloc-owned, NUL-terminated uppercase bases, or NULL.
 * Requests larger than the cache window are fetched without retaining them.
 * IUPAC interpretation belongs to the caller, not this shared reference cache. */
char *duckhts_reference_fetch(duckhts_reference_entry_t *entry,
    const char *chrom, hts_pos_t start0, hts_pos_t end0,
    duckhts_reference_end_policy_t end_policy);

void duckhts_reference_uppercase(char *bases);

#endif
