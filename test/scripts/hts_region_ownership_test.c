/* HTSlib multi-region iterators own their region list from the call on,
 * also when creation fails (HTSlib 1.25). Run under ASan with leak detection:
 * a caller free after failure is a double free, a missing HTSlib free a leak. */
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

static int fail_iterator_calloc;
static int failed_iterator_callocs;

void *__real_calloc(size_t count, size_t width);

void *__wrap_calloc(size_t count, size_t width) {
    if (fail_iterator_calloc && count == 1u && width == sizeof(hts_itr_t)) {
        fail_iterator_calloc = 0;
        failed_iterator_callocs++;
        return NULL;
    }
    return __real_calloc(count, width);
}

static hts_reglist_t *make_region_list(const char *name) {
    hts_reglist_t *regions = calloc(1u, sizeof(*regions));
    assert(regions != NULL);
    regions[0].reg = name;
    regions[0].intervals = calloc(1u, sizeof(*regions[0].intervals));
    assert(regions[0].intervals != NULL);
    regions[0].tid = 0;
    regions[0].count = 1u;
    regions[0].min_beg = 0;
    regions[0].max_end = 1;
    regions[0].intervals[0].beg = 0;
    regions[0].intervals[0].end = 1;
    return regions;
}

static int fail_name_lookup(void *header, const char *name) {
    (void)header;
    (void)name;
    return -2;
}

static int fail_query(const hts_idx_t *index, hts_itr_t *iterator) {
    (void)index;
    (void)iterator;
    return -1;
}

static int unused_read(BGZF *file, void *data, void *record, int *tid,
                       hts_pos_t *begin, hts_pos_t *end) {
    (void)file;
    (void)data;
    (void)record;
    (void)tid;
    (void)begin;
    (void)end;
    return -1;
}

static int unused_seek(void *file, int64_t offset, int where) {
    (void)file;
    (void)offset;
    (void)where;
    return -1;
}

static int64_t unused_tell(void *file) {
    (void)file;
    return -1;
}

static void check_failure_takes_list(void) {
    hts_reglist_t *regions = make_region_list("invalid");
    hts_itr_t *iterator = hts_itr_regions(
        NULL, regions, 1, fail_name_lookup, NULL, fail_query, unused_read,
        unused_seek, unused_tell);
    assert(iterator == NULL);

    regions = make_region_list(NULL);
    iterator = hts_itr_regions(
        NULL, regions, 1, fail_name_lookup, NULL, fail_query, unused_read,
        unused_seek, unused_tell);
    assert(iterator == NULL);
}

static void check_sam_allocation_failure(const char *path,
                                         const char *index_path) {
    samFile *file = sam_open(path, "r");
    sam_hdr_t *header;
    hts_idx_t *index;
    hts_reglist_t *regions;
    hts_itr_t *iterator;
    int failures_before = failed_iterator_callocs;

    assert(file != NULL);
    header = sam_hdr_read(file);
    assert(header != NULL);
    index = sam_index_load3(file, path, index_path, HTS_IDX_SILENT_FAIL);
    assert(index != NULL);
    regions = make_region_list(NULL);

    fail_iterator_calloc = 1;
    iterator = sam_itr_regions(index, header, regions, 1u);
    assert(iterator == NULL);
    assert(fail_iterator_calloc == 0);
    assert(failed_iterator_callocs == failures_before + 1);

    hts_idx_destroy(index);
    sam_hdr_destroy(header);
    assert(sam_close(file) == 0);
}

int main(int argc, char **argv) {
    assert(argc == 5);
    check_failure_takes_list();
    check_sam_allocation_failure(argv[1], argv[2]);
    check_sam_allocation_failure(argv[3], argv[4]);
    return 0;
}
