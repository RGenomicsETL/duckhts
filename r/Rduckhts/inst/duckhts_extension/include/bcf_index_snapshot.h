#ifndef DUCKHTS_BCF_INDEX_SNAPSHOT_H
#define DUCKHTS_BCF_INDEX_SNAPSHOT_H

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

/* Zero-initialize, load once, then borrow read-only until every iterator has
 * been destroyed. Owns parsed HTSlib index memory, not a data file or header.
 * Query/read operations may share this snapshot; records, headers, iterators,
 * kstrings and file handles must remain worker-local. */
typedef struct {
    hts_idx_t *bcf;
    tbx_t *tabix;
} duckhts_bcf_index_t;

static inline void duckhts_bcf_index_destroy(duckhts_bcf_index_t *index) {
    if (index->bcf) hts_idx_destroy(index->bcf);
    if (index->tabix) tbx_destroy(index->tabix);
    index->bcf = NULL;
    index->tabix = NULL;
}

/* Return 1 for a loaded snapshot, 0 for no usable index, -1 if finalizing its
 * read-only state fails. File/index association is the caller's input contract;
 * retaining an index does not validate an unrelated index supplied initially. */
static inline int duckhts_bcf_index_load(duckhts_bcf_index_t *index,
                                       enum htsExactFormat format,
                                       const char *path, const char *index_path,
                                       int flags) {
    if (format == bcf) {
        index->bcf = bcf_index_load3(path, index_path, flags);
    } else {
        index->tabix = tbx_index_load3(path, index_path, flags);
        if (!index->tabix) index->bcf = bcf_index_load3(path, index_path, flags);
    }
    if (index->tabix && !index->tabix->dict) {
        /* HTSlib 1.24's get_tid(..., is_add=0) lazily creates even an empty
         * dictionary. Do that before publication: query/readrec then only
         * look up names, and all per-query mutations stay in the iterator. */
        (void)tbx_name2id(index->tabix, "");
        if (!index->tabix->dict) {
            duckhts_bcf_index_destroy(index);
            return -1;
        }
    }
    return index->bcf || index->tabix;
}

#endif
