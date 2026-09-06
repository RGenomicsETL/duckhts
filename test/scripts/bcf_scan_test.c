/* Canonical scanning and decode checks without DuckDB. Each concurrent scan
 * owns its mutable state; records and decoded genotypes have explicit owners. */
#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../../src/include/bcf_scan.h"
#include "../../src/include/bcf_genotypes.h"

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
    duckhts_bcf_scan_t reader = {0};
    char error[512];
    assert(duckhts_bcf_scan_open(&reader, input->path, input->index, 0,
        DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
    bcf_hdr_t *hdr = reader.hdr;
    bcf1_t *record = bcf_init();
    assert(hdr && record);
    kstring_t formatted = KS_INITIALIZE;
    char *regions[] = {"chr3:20-40", "chr1:1-10", "chr3:20-20"};
    char *single = selection == 1 ? "chr3:20-40" : "absent:1-10";
    assert(duckhts_bcf_scan_regions(&reader, selection == 0 ? regions : &single,
                                    selection == 0 ? 3 : 1, error, sizeof(error)));
    unsigned counts[4] = {0};
    int ret = -1;
    if (!input->empty && selection < 2) assert(reader.itr);
    while ((ret = duckhts_bcf_scan_next(&reader, record)) >= 0) {
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
    free(formatted.s);
    bcf_destroy(record);
    duckhts_bcf_scan_close(&reader);
    assert(!reader.fp && !reader.hdr && !reader.itr && !reader.line.s);
    duckhts_bcf_scan_close(&reader); // Partial/repeated cleanup is safe.
}

static void *worker(void *arg) {
    for (int iteration = 0; iteration < 100; iteration++) {
        for (int selection = 0; selection < 3; selection++) scan(arg, selection);
    }
    return NULL;
}

static void genotypes(const char *path, const duckhts_bcf_index_t *index) {
    /* Literal decoded GT encodings, independent of the SQL string formatter.
     * HTSlib 1.24 vcf.c updatephasing/VCF parsing sets the first-slot phase bit
     * for a known haploid or fully phased call, including pre-4.4 inputs.
     * Rows: known, duplicate, allele, last; columns retain original sample IDs. */
    const int32_t expected[4][2][2] = {
        {{2,4},{4,4}}, {{0,4},{2,2}}, {{4,4},{0,0}}, {{5,0},{3,5}}
    };
    const char *ids[] = {"known", "duplicate", "allele", "last"};
    const char *samples[] = {"-", "S1", "S2", "", "^S1,S2", "S2,S1,S2", "^S1", NULL};
    const int selected_counts[] = {2, 1, 1, 0, 0, 2, 1, 2};
    for (int selection = 0; selection < 8; selection++) for (int indexed = 0; indexed < 2; indexed++) {
        duckhts_bcf_scan_t reader = {0};
        char error[512];
        assert(duckhts_bcf_scan_open(&reader, path, index, 0,
            DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
        duckhts_bcf_samples_t selected = {0};
        assert(duckhts_bcf_samples_build(&selected, reader.hdr, samples[selection], error, sizeof(error)));
        assert(selected.original_count == 2 && selected.count == selected_counts[selection]);
        assert(bcf_hdr_nsamples(reader.hdr) == 2); // Building a plan does not mutate its input.
        assert(duckhts_bcf_samples_apply(&selected, reader.hdr, error, sizeof(error)));
        if (indexed) {
            char *regions[] = {"chr1", "chr3", "chr3:20-20"};
            assert(duckhts_bcf_scan_regions(&reader, regions, 3, error, sizeof(error)));
        }
        bcf1_t *record = bcf_init();
        assert(record);
        int32_t *gt = NULL;
        duckhts_bcf_genotypes_t typed = {0};
        int capacity = 0, ret;
        unsigned counts[4] = {0};
        while ((ret = duckhts_bcf_scan_next(&reader, record)) >= 0) {
            assert(bcf_unpack(record, BCF_UN_STR) == 0);
            int row = 0;
            while (row < 4 && strcmp(record->d.id, ids[row]) != 0) row++;
            assert(row < 4);
            counts[row]++;
            int nsamples = bcf_hdr_nsamples(reader.hdr);
            assert(nsamples == selected_counts[selection]);
            assert(record->n_sample == (unsigned)nsamples);
            if (!nsamples) continue; // A site still exists with no selected calls.
            int header_id = bcf_hdr_id2int(reader.hdr, BCF_DT_ID, "GT");
            assert(duckhts_bcf_check_field_type(reader.hdr, record, DUCKHTS_BCF_FIELD_GT,
                header_id, BCF_HT_STR, "test", error, sizeof(error)));
            int values = bcf_get_genotypes(reader.hdr, record, &gt, &capacity);
            assert(values > 0 && values % nsamples == 0);
            assert(duckhts_bcf_genotypes_decode(&typed, reader.hdr, record,
                DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
            assert(typed.samples == nsamples && !typed.ps_present);
            assert(typed.gt_stride == values / nsamples);
            assert(memcmp(typed.gt, gt, (size_t)values * sizeof(*gt)) == 0);
            assert(duckhts_bcf_check_format_width("test", "GT", reader.hdr, record,
                values, nsamples, error, sizeof(error)));
            int stride = values / nsamples;
            for (int sample = 0; sample < nsamples; sample++) {
                int original = selected_counts[selection] == 2 ? sample : selection == 1 ? 0 : 1;
                assert(selected.indices[sample] == (unsigned)original);
                assert(strcmp(selected.names[sample], original ? "S2" : "S1") == 0);
                assert(strcmp(reader.hdr->samples[sample], original ? "S2" : "S1") == 0);
                int slots = row == 3 && original == 0 ? 1 : 2;
                assert(stride >= slots);
                for (int slot = 0; slot < slots; slot++) {
                    if (gt[sample * stride + slot] != expected[row][original][slot]) {
                        fprintf(stderr, "%s selection=%d indexed=%d row=%d sample=%d slot=%d: %d != %d\n",
                            path, selection, indexed, row, original, slot,
                            gt[sample * stride + slot], expected[row][original][slot]);
                    }
                    assert(gt[sample * stride + slot] == expected[row][original][slot]);
                }
                for (int slot = slots; slot < stride; slot++) {
                    assert(gt[sample * stride + slot] == bcf_int32_vector_end);
                }
            }
        }
        assert(ret == -1 && counts[0] == 1 && counts[1] == 2 && counts[2] == 1 && counts[3] == 1);
        free(gt);
        duckhts_bcf_genotypes_destroy(&typed);
        bcf_destroy(record);
        duckhts_bcf_scan_close(&reader);
        duckhts_bcf_samples_destroy(&selected);
        duckhts_bcf_samples_destroy(&selected);
    }
}

static void genotype_values(void) {
    bcf_hdr_t *header = bcf_hdr_init("w");
    bcf1_t *record = bcf_init();
    assert(header && record);
    assert(bcf_hdr_append(header, "##contig=<ID=chrG,length=100>") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">") == 0);
    assert(bcf_hdr_add_sample(header, "S1") == 0 && bcf_hdr_add_sample(header, "S2") == 0);
    assert(bcf_hdr_add_sample(header, NULL) == 0 && bcf_hdr_sync(header) == 0);
    record->rid = 0;
    record->pos = 0;
    kstring_t alleles = KS_INITIALIZE;
    assert(kputc('A', &alleles) >= 0);
    for (int i = 1; i <= 20000; i++) assert(ksprintf(&alleles, ",<ALT%d>", i) >= 0);
    assert(bcf_update_alleles_str(header, record, alleles.s) == 0);
    ks_free(&alleles);
    int32_t gt[] = {bcf_gt_phased(20000), 1, bcf_gt_unphased(1), bcf_gt_unphased(0),
                    bcf_gt_phased(1), bcf_int32_vector_end, bcf_int32_vector_end, bcf_int32_vector_end};
    int32_t ps[] = {INT32_MAX, bcf_int32_missing};
    assert(bcf_update_genotypes(header, record, gt, 8) == 0);
    assert(bcf_update_format_int32(header, record, "PS", ps, 2) == 0);
    assert(bcf_get_fmt(header, record, "GT")->type == BCF_BT_INT32);
    duckhts_bcf_genotypes_t values = {0};
    char error[512];
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(values.samples == 2 && values.gt_stride == 4 && values.ps_present);
    assert(memcmp(values.gt, gt, sizeof(gt)) == 0 && memcmp(values.ps, ps, sizeof(ps)) == 0);
    assert(duckhts_bcf_genotype_ploidy(values.gt, 4) == 4);
    assert(duckhts_bcf_genotype_ploidy(values.gt + 4, 4) == 1);
    assert(duckhts_bcf_genotype_has_alt(values.gt, 4));
    int32_t no_alt[] = {0, 1, bcf_gt_unphased(0)};
    assert(!duckhts_bcf_genotype_has_alt(no_alt, 3));
    int32_t wide_ps[] = {10, 11, 20, bcf_int32_vector_end};
    assert(bcf_update_format_int32(header, record, "PS", wide_ps, 4) == 0);
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(strstr(error, "one value per sample"));
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(values.gt_stride == 4 && !values.ps_present);
    const char *strings[] = {"10", "20"};
    assert(bcf_update_format_string(header, record, "PS", strings, 2) == 0);
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(strstr(error, "encoded BCF type CHAR"));
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(values.gt_stride == 4 && !values.ps_present);
    assert(bcf_update_format_int32(header, record, "PS", ps, 2) == 0);
    // bcf_get_format_values() stops at the first vector-end and fills the
    // remaining decoded slots with vector-end, even if raw padding has values.
    // The typed view follows that HTSlib contract, not a second raw decoder.
    gt[1] = bcf_int32_vector_end;
    assert(bcf_update_genotypes(header, record, gt, 8) == 0);
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(values.gt_stride == 4 && duckhts_bcf_genotype_ploidy(values.gt, 4) == 1);
    assert(values.gt[1] == bcf_int32_vector_end && values.gt[2] == bcf_int32_vector_end &&
           values.gt[3] == bcf_int32_vector_end);
    gt[0] = bcf_gt_phased(20001); // Index outside this record's ALT dictionary.
    assert(bcf_update_genotypes(header, record, gt, 8) == 0);
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(strstr(error, "invalid FORMAT/GT allele"));
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(values.gt_stride == 0 && values.ps_present);
    bcf_fmt_t *format = bcf_get_fmt(header, record, "GT");
    int width = format->n;
    format->n = INT32_MAX; // Refuse signed count multiplication before HTSlib decoding.
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(strstr(error, "exceeds the supported decoded-value capacity"));
    format->n = width;
    duckhts_bcf_samples_t selection = {0};
    assert(duckhts_bcf_samples_build(&selection, header, "S2", error, sizeof(error)));
    char *reversed[] = {"S2", "S1"};
    int mapping[2];
    bcf_hdr_t *changed = bcf_hdr_subset(header, 2, reversed, mapping);
    assert(changed && !duckhts_bcf_samples_apply(&selection, changed, error, sizeof(error)));
    assert(strstr(error, "dictionary changed since bind"));
    bcf_hdr_destroy(changed);
    duckhts_bcf_samples_destroy(&selection);
    assert(!duckhts_bcf_samples_build(&selection, header, "S2,absent", error, sizeof(error)));
    assert(strstr(error, "selector item 2 is not in the header"));
    assert(!selection.names && !selection.indices && !selection.selector);
    duckhts_bcf_genotypes_destroy(&values);
    duckhts_bcf_genotypes_destroy(&values);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
}

static void literal_contigs(void) {
    const char *paths[] = {"test/data/bcf_literal_contigs.bcf", "test/data/bcf_literal_contigs.vcf.gz"};
    const char *names[] = {"chr1", "chr1:100-200", "absent"};
    const char *expected[] = {
        "chr1\t150\tordinary\tA\tC\t60\tPASS\t.\tGT\t0/1\n",
        "chr1:100-200\t10\tliteral\tG\tT\t50\tPASS\t.\tGT\t1|1\n"
    };
    for (int kind = 0; kind < 2; kind++) {
        duckhts_bcf_index_t index = {0};
        assert(duckhts_bcf_index_load(&index, kind == 0 ? bcf : vcf, paths[kind], NULL, 0) == 1);
        duckhts_bcf_scan_t reader = {0};
        char error[512];
        assert(duckhts_bcf_scan_open(&reader, paths[kind], &index, 0,
            DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
        bcf1_t *record = bcf_init();
        assert(record);
        kstring_t text = KS_INITIALIZE;
        for (int contig = 0; contig < 3; contig++) {
            assert(duckhts_bcf_scan_contig(&reader, names[contig], error, sizeof(error)));
            int count = 0, ret;
            while ((ret = duckhts_bcf_scan_next(&reader, record)) >= 0) {
                text.l = 0;
                assert(vcf_format1(reader.hdr, record, &text) == 0);
                assert(contig < 2 && strcmp(text.s, expected[contig]) == 0);
                count++;
            }
            assert(ret == -1 && count == (contig == 0 ? 1 : contig == 1 ? 2 : 0));
        }
        // Selecting an iterator does not consume or invalidate the caller's record.
        char *region = "{chr1}:140-160";
        assert(duckhts_bcf_scan_regions(&reader, &region, 1, error, sizeof(error)));
        text.l = 0;
        assert(vcf_format1(reader.hdr, record, &text) == 0);
        assert(strcmp(text.s, expected[1]) == 0);
        assert(duckhts_bcf_scan_next(&reader, record) >= 0);
        text.l = 0;
        assert(vcf_format1(reader.hdr, record, &text) == 0);
        assert(strcmp(text.s, expected[0]) == 0);
        assert(duckhts_bcf_scan_next(&reader, record) == -1);
        region = "chr1:200-100";
        assert(!duckhts_bcf_scan_regions(&reader, &region, 1, error, sizeof(error)));
        assert(duckhts_bcf_scan_next(&reader, record) == -1); // Never stream on invalid selection.
        free(text.s);
        bcf_destroy(record);
        duckhts_bcf_scan_close(&reader);
        duckhts_bcf_index_destroy(&index);
    }
}

static void decode_errors(void) {
    const char *paths[] = {"bcf_info_type_clash.bcf", "bcf_info_str_clash.bcf", "bcf_format_type_clash.bcf"};
    const char *tags[] = {"DP", "NN", "XX"};
    const char *messages[] = {
        "read_bcf: INFO/DP encoded BCF type CHAR does not match header Type=Integer at chr1:10",
        "read_bcf: INFO/NN encoded BCF type INT8 does not match header Type=String at chr1:10",
        "read_bcf: FORMAT/XX encoded BCF type CHAR does not match header Type=Integer at chr1:10"
    };
    for (int i = 0; i < 3; i++) {
        char path[256], error[512];
        snprintf(path, sizeof(path), "test/data/%s", paths[i]);
        duckhts_bcf_scan_t reader = {0};
        assert(duckhts_bcf_scan_open(&reader, path, NULL, 0,
            DUCKHTS_HTS_IO_PROFILE_METADATA, "read_bcf", error, sizeof(error)));
        bcf1_t *record = bcf_init();
        assert(record && duckhts_bcf_scan_next(&reader, record) >= 0);
        int id = bcf_hdr_id2int(reader.hdr, BCF_DT_ID, tags[i]);
        int type = i == 1 ? BCF_HT_STR : BCF_HT_INT;
        assert(!duckhts_bcf_check_field_type(reader.hdr, record,
            i == 2 ? DUCKHTS_BCF_FIELD_FORMAT : DUCKHTS_BCF_FIELD_INFO,
            id, type, "read_bcf", error, sizeof(error)));
        assert(strcmp(error, messages[i]) == 0);
        for (int ret = -5; ret <= 5; ret++) {
            duckhts_bcf_decode_status_t expected = ret >= 0 || ret == -3 ? DUCKHTS_BCF_DECODE_OK
                : ret == -2 ? DUCKHTS_BCF_DECODE_TYPE_MISMATCH : DUCKHTS_BCF_DECODE_FATAL;
            assert(duckhts_bcf_decode_status("read_bcf", "FORMAT", "XX", reader.hdr, record,
                ret, error, sizeof(error)) == expected);
            if (ret == -4) assert(strstr(error, "out of memory decoding FORMAT/XX at chr1:10"));
        }
        assert(!duckhts_bcf_check_format_width("read_bcf", "XX", reader.hdr, record,
            3, 2, error, sizeof(error)));
        assert(strcmp(error, "read_bcf: FORMAT/XX decoded value count 3 is not divisible by sample count 2 at chr1:10") == 0);
        bcf_destroy(record);
        duckhts_bcf_scan_close(&reader);
    }
    duckhts_bcf_scan_t reader = {0};
    char error[512];
    assert(!duckhts_bcf_scan_open(&reader, "test/data/no-such-bcf-file", NULL, 0,
        DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
    assert(!reader.fp && !reader.hdr);
    assert(duckhts_bcf_scan_open(&reader, "test/data/malformed_bad_pos.vcf", NULL, 0,
        DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
    bcf1_t *record = bcf_init();
    assert(record);
    int ret;
    while ((ret = duckhts_bcf_scan_next(&reader, record)) >= 0) {}
    assert(ret < -1); // A malformed physical record must not become EOF.
    bcf_destroy(record);
    duckhts_bcf_scan_close(&reader);
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
        genotypes(path, &index);
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
    decode_errors();
    genotype_values();
    literal_contigs();
    puts("BCF scanner: CSI/TBI, shifted offsets, 7200 concurrent exact-row scans, selected raw GT and decode errors: OK");
    return 0;
}
