/* Parsed index ownership without DuckDB. Each concurrent scan owns every
 * mutable HTSlib object; the actual production index helper is shared. */
#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../../src/include/bcf_index_snapshot.h"

static const char *rows[] = {
    "chr1\t10\tknown\tA\tC\t60\tPASS\tDP=7\tGT\t0/1\t1/1\n",
    "chr3\t20\tduplicate\tG\tT\t50\tPASS\tDP=8\tGT\t./1\t0/0\n",
    "chr3\t20\tallele\tG\tA\t30\tPASS\tDP=9\tGT\t1/1\t./.\n",
    "chr3\t40\tlast\tT\tG\t.\tPASS\tDP=.\tGT\t1\t0|1\n"
};

typedef struct {
    const duckhts_bcf_index_t *index;
    const char *path;
    int empty;
} scan_input_t;

static void scan(const scan_input_t *input, int selection) {
    htsFile *fp = bcf_open(input->path, "r");
    assert(fp);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *record = bcf_init();
    assert(hdr && record);
    kstring_t line = KS_INITIALIZE, formatted = KS_INITIALIZE;
    char *regions[] = {"chr3:20-40", "chr1:1-10", "chr3:20-20"};
    const duckhts_bcf_index_t *index = input->index;
    const char *single = selection == 1 ? "chr3:20-40" : "absent:1-10";
    hts_itr_t *itr = index->bcf
        ? (selection == 0 ? bcf_itr_regarray(index->bcf, hdr, regions, 3)
                          : bcf_itr_querys(index->bcf, hdr, single))
        : (selection == 0 ? tbx_itr_regarray(index->tabix, regions, 3)
                          : tbx_itr_querys(index->tabix, single));
    unsigned counts[4] = {0};
    int ret = -1;
    if (!input->empty && selection < 2) assert(itr);
    while (itr) {
        if (index->tabix) {
            ret = tbx_itr_next(fp, index->tabix, itr, &line);
            if (ret >= 0) ret = vcf_parse1(&line, hdr, record);
        } else ret = bcf_itr_next(fp, itr, record);
        if (ret < 0) break;
        formatted.l = 0;
        assert(vcf_format1(hdr, record, &formatted) == 0);
        unsigned i;
        for (i = 0; i < 4; i++) if (strcmp(rows[i], formatted.s) == 0) break;
        assert(i < 4);
        counts[i]++;
    }
    assert(ret == -1);
    int present = !input->empty && selection < 2;
    assert(counts[0] == (unsigned)(present && selection == 0));
    assert(counts[1] == (unsigned)(present ? 2 : 0));
    assert(counts[2] == (unsigned)present && counts[3] == (unsigned)present);
    if (itr) hts_itr_destroy(itr);
    free(line.s);
    free(formatted.s);
    bcf_destroy(record);
    bcf_hdr_destroy(hdr);
    assert(bcf_close(fp) == 0);
}

static void *worker(void *arg) {
    for (int iteration = 0; iteration < 100; iteration++) {
        for (int selection = 0; selection < 3; selection++) scan(arg, selection);
    }
    return NULL;
}

/* Counts and dictionary agree, but virtual offsets differ. Demonstrate this
 * directly, not via index size/inode comparisons masquerading as identity. */
static void check_shifted(const duckhts_bcf_index_t *original,
                          const duckhts_bcf_index_t *shifted) {
    const hts_idx_t *a = original->bcf ? original->bcf : original->tabix->idx;
    const hts_idx_t *b = shifted->bcf ? shifted->bcf : shifted->tabix->idx;
    assert(hts_idx_nseq(a) == hts_idx_nseq(b));
    assert(hts_idx_get_n_no_coor(a) == hts_idx_get_n_no_coor(b));
    if (original->tabix) {
        int an = 0, bn = 0;
        const char **names_a = tbx_seqnames(original->tabix, &an);
        const char **names_b = tbx_seqnames(shifted->tabix, &bn);
        assert(names_a && names_b && an == bn);
        for (int i = 0; i < an; i++) assert(strcmp(names_a[i], names_b[i]) == 0);
        free(names_a);
        free(names_b);
    }
    int different = 0;
    for (int tid = 0; tid < hts_idx_nseq(a); tid++) {
        uint64_t am = 0, au = 0, bm = 0, bu = 0;
        int ar = hts_idx_get_stat(a, tid, &am, &au);
        int br = hts_idx_get_stat(b, tid, &bm, &bu);
        assert(ar == br && am == bm && au == bu);
        hts_itr_t *ai = hts_itr_query(a, tid, 0, HTS_POS_MAX, NULL);
        hts_itr_t *bi = hts_itr_query(b, tid, 0, HTS_POS_MAX, NULL);
        if (!ai || !bi) { assert(!ai && !bi); continue; }
        assert(ai->n_off == bi->n_off);
        for (int i = 0; i < ai->n_off; i++) {
            different |= ai->off[i].u != bi->off[i].u || ai->off[i].v != bi->off[i].v;
        }
        hts_itr_destroy(ai);
        hts_itr_destroy(bi);
    }
    assert(different);
}

int main(int argc, char **argv) {
    assert(argc == 2); // Caller-owned temporary directory; no fixture mutation.
    hts_set_log_level(HTS_LOG_ERROR);
    const char *suffixes[] = {"bcf", "full.vcf.gz", "full.vcf.gz"};
    const char *index_suffixes[] = {"csi", "index.tbi", "index.csi"};
    for (int kind = 0; kind < 3; kind++) {
        char path[512], shifted_path[512], empty_path[512], index_path[512], shifted_index[1024];
        snprintf(path, sizeof(path), "test/data/bcf_scan_contigs.%s", suffixes[kind]);
        const char *format_suffix = kind == 0 ? "bcf" : "vcf.gz";
        snprintf(shifted_path, sizeof(shifted_path), "test/data/bcf_scan_contigs.shifted.%s", format_suffix);
        snprintf(empty_path, sizeof(empty_path), "test/data/bcf_scan_contigs.empty.%s", format_suffix);
        assert(snprintf(index_path, sizeof(index_path), "%s/snapshot.index", argv[1]) < (int)sizeof(index_path));
        snprintf(shifted_index, sizeof(shifted_index), "%s.%s", shifted_path, index_suffixes[kind]);
        enum htsExactFormat format = kind == 0 ? bcf : vcf;
        int min_shift = kind == 1 ? 0 : 14;
        assert(bcf_index_build3(path, index_path, min_shift, 0) == 0);
        duckhts_bcf_index_t index = {0}, shifted = {0};
        assert(duckhts_bcf_index_load(&index, format, path, index_path, 0) == 1);
        assert(duckhts_bcf_index_load(&shifted, format, shifted_path, shifted_index, 0) == 1);
        check_shifted(&index, &shifted);
        duckhts_bcf_index_destroy(&shifted);
        assert(unlink(index_path) == 0);
        scan_input_t input = {&index, path, 0};
        scan(&input, 0); // Removal after loading does not invalidate memory.
        assert(bcf_index_build3(empty_path, index_path, min_shift, 0) == 0);
        pthread_t threads[4];
        for (unsigned i = 0; i < 4; i++) assert(pthread_create(&threads[i], NULL, worker, &input) == 0);
        for (unsigned i = 0; i < 4; i++) assert(pthread_join(threads[i], NULL) == 0);
        duckhts_bcf_index_destroy(&index);
        assert(!index.bcf && !index.tabix);
        assert(duckhts_bcf_index_load(&index, format, empty_path, index_path, 0) == 1);
        assert(!index.tabix || index.tabix->dict); // No lazy shared dictionary allocation.
        input.path = empty_path;
        input.empty = 1;
        for (unsigned i = 0; i < 4; i++) assert(pthread_create(&threads[i], NULL, worker, &input) == 0);
        for (unsigned i = 0; i < 4; i++) assert(pthread_join(threads[i], NULL) == 0);
        duckhts_bcf_index_destroy(&index);
        assert(unlink(index_path) == 0);
    }
    puts("BCF index snapshot: CSI/TBI, shifted offsets, removal/replacement, 7200 concurrent exact-row scans: OK");
    return 0;
}
