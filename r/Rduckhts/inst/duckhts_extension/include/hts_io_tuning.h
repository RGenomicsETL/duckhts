#ifndef DUCKHTS_HTS_IO_TUNING_H
#define DUCKHTS_HTS_IO_TUNING_H

#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

/* htslib exports this function but declares it only inside hwrite() in the
 * public header.  Give direct hFILE consumers a normal file-scope prototype. */
extern int hfile_set_blksize(hFILE *fp, size_t bufsiz);

/* Remote-I/O tuning policy.
 *
 * Principle:
 *   - only tune paths that htslib classifies as remote
 *   - pick larger buffers for long streaming scans / index builds
 *   - keep indexed-region / seek-heavy profiles smaller
 *   - use more conservative budgets on Emscripten/browser builds, where our
 *     http/https transport is synchronous XHR inside a worker rather than
 *     libcurl sockets on a native host
 *
 * These settings are advisory. hts_set_opt() may no-op for formats or backends
 * that do not honour a specific option; callers should not treat that as a
 * fatal error.
 */

typedef enum {
    DUCKHTS_HTS_IO_PROFILE_METADATA = 0,
    DUCKHTS_HTS_IO_PROFILE_STREAMING = 1,
    DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION = 2
} duckhts_hts_io_profile_t;

static inline int duckhts_is_remote_path(const char *path) {
    return path && path[0] != '\0' && hisremote(path);
}

static inline int duckhts_remote_block_size_bytes(duckhts_hts_io_profile_t profile) {
#if defined(__EMSCRIPTEN__)
    switch (profile) {
        case DUCKHTS_HTS_IO_PROFILE_STREAMING:
            return 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION:
            return 256 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_METADATA:
        default:
            return 0;
    }
#else
    switch (profile) {
        case DUCKHTS_HTS_IO_PROFILE_STREAMING:
            return 4 * 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION:
            return 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_METADATA:
        default:
            return 0;
    }
#endif
}

static inline int duckhts_remote_cache_size_bytes(duckhts_hts_io_profile_t profile) {
#if defined(__EMSCRIPTEN__)
    switch (profile) {
        case DUCKHTS_HTS_IO_PROFILE_STREAMING:
            return 32 * 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION:
            return 8 * 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_METADATA:
        default:
            return 0;
    }
#else
    switch (profile) {
        case DUCKHTS_HTS_IO_PROFILE_STREAMING:
            return 256 * 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION:
            return 64 * 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_METADATA:
        default:
            return 0;
    }
#endif
}

/* Raw seekable hFILE consumers do not have BGZF's block cache.  Keep their
 * read-ahead smaller than HTS_OPT_BLOCK_SIZE so an index walk does not fetch a
 * megabyte after every non-local seek. */
static inline int duckhts_remote_hfile_block_size_bytes(duckhts_hts_io_profile_t profile) {
#if defined(__EMSCRIPTEN__)
    switch (profile) {
        case DUCKHTS_HTS_IO_PROFILE_STREAMING:
            return 512 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION:
            return 128 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_METADATA:
        default:
            return 0;
    }
#else
    switch (profile) {
        case DUCKHTS_HTS_IO_PROFILE_STREAMING:
            return 1024 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION:
            return 256 * 1024;
        case DUCKHTS_HTS_IO_PROFILE_METADATA:
        default:
            return 0;
    }
#endif
}

static inline void duckhts_apply_remote_hts_tuning(htsFile *fp,
                                                    const char *path,
                                                    duckhts_hts_io_profile_t profile) {
    int block_size;
    int cache_size;

    if (!fp || !duckhts_is_remote_path(path) || profile == DUCKHTS_HTS_IO_PROFILE_METADATA) {
        return;
    }

    block_size = duckhts_remote_block_size_bytes(profile);
    cache_size = duckhts_remote_cache_size_bytes(profile);

    if (block_size > 0) {
        (void)hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, block_size);
    }
    if (cache_size > 0) {
        (void)hts_set_opt(fp, HTS_OPT_CACHE_SIZE, cache_size);
    }
}

static inline void duckhts_apply_remote_hfile_tuning(hFILE *fp,
                                                      const char *path,
                                                      duckhts_hts_io_profile_t profile) {
    int block_size;

    if (!fp || !duckhts_is_remote_path(path) || profile == DUCKHTS_HTS_IO_PROFILE_METADATA) {
        return;
    }

    block_size = duckhts_remote_hfile_block_size_bytes(profile);
    if (block_size > 0) {
        (void)hfile_set_blksize(fp, (size_t)block_size);
    }
}

static inline void duckhts_apply_remote_faidx_tuning(faidx_t *fai,
                                                      const char *path,
                                                      duckhts_hts_io_profile_t profile) {
    int cache_size;

    if (!fai || !duckhts_is_remote_path(path) || profile == DUCKHTS_HTS_IO_PROFILE_METADATA) {
        return;
    }

    cache_size = duckhts_remote_cache_size_bytes(profile);
    if (cache_size > 0) {
        fai_set_cache_size(fai, cache_size);
    }
}

#endif
