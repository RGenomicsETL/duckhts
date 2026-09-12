/* Canonical scanning and decode checks without DuckDB. Each concurrent scan
 * owns its mutable state; records and decoded genotypes have explicit owners. */
#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <htslib/bgzf.h>
#include <htslib/hts_log.h>
#include "../../src/include/bcf_scan.h"
#include "../../src/include/bcf_genotypes.h"
#include "../../src/include/bcf_format.h"

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
            assert(typed.samples == nsamples && !typed.ps_stride);
            assert(typed.gt_stride == values / nsamples);
            assert(memcmp(typed.gt, gt, (size_t)values * sizeof(*gt)) == 0);
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

static void raw_genotypes(const char *directory) {
    const char *header = "##fileformat=VCFv4.4\n"
        "##contig=<ID=chrR,length=1000>\n"
        "##INFO=<ID=TEXT,Number=1,Type=String,Description=\"Text\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n";
    const char *records[] = {
        "chrR\t10\tprefix\tA\tC,G\t.\t.\t.\tGT\t0|1\t|0|1\t/0|1/2\n",
        "chrR\t10\tprefix\tA\tC,G\t.\t.\t.\tGT\t0|1\t|0|1\t/0|1/2\n",
        "chrR\t20\tmixed\tA\tC,G\t.\t.\t.\tGT:DP\t|0|1/2:8\t.:.\t00|01:9\n",
        "chrR\t30\tabsent\tA\tC,G\t.\t.\t.\tDP\t7\t.\t9\n",
        "chrR\t40\tnone\tA\tC,G\t.\t.\t.\t.\t.\t.\t.\n",
        "chrR\t50\tshort\tA\tC,G\t.\t.\t.\tGT:DP\t1\t./.\t0|1\n",
        "chrR\t60\ttrailing\tA\tC,G\t.\t.\t.\tDP:GT\t5:|1/0\t7\t8:.\n"
    };
    const char *expected[][3] = {
        {"0|1", "|0|1", "/0|1/2"}, {"0|1", "|0|1", "/0|1/2"},
        {"|0|1/2", ".", "00|01"}, {NULL, NULL, NULL}, {NULL, NULL, NULL},
        {"1", "./.", "0|1"}, {"|1/0", NULL, "."}, {"1|0", "|0/1", "."}
    };
    char plain[512], compressed[512], index_path[512];
    assert(snprintf(plain, sizeof(plain), "%s/raw.vcf", directory) < (int)sizeof(plain));
    assert(snprintf(compressed, sizeof(compressed), "%s/raw.vcf.gz", directory) < (int)sizeof(compressed));
    assert(snprintf(index_path, sizeof(index_path), "%s/raw.csi", directory) < (int)sizeof(index_path));
    FILE *file = fopen(plain, "wb");
    assert(file && fputs(header, file) >= 0);
    for (size_t i = 0; i < sizeof(records) / sizeof(*records); i++) assert(fputs(records[i], file) >= 0);
    const char *large_prefix = "chrR\t70\tlarge\tA\tC,G\t.\t.\tTEXT=";
    const char *large_suffix = "\tGT\t1|0\t|0/1\t.\n";
    const size_t line_length = 131070u;
    size_t fill = line_length - strlen(large_prefix) - strlen(large_suffix) + 1u;
    assert(fputs(large_prefix, file) >= 0);
    for (size_t i = 0; i < fill; i++) assert(fputc('A', file) != EOF);
    assert(fputs(large_suffix, file) >= 0 && fclose(file) == 0);
    file = fopen(plain, "rb");
    BGZF *bgzf = bgzf_open(compressed, "w");
    assert(file && bgzf);
    char buffer[4096];
    size_t bytes;
    while ((bytes = fread(buffer, 1, sizeof(buffer), file)) != 0)
        assert(bgzf_write(bgzf, buffer, bytes) == (ssize_t)bytes);
    assert(!ferror(file) && fclose(file) == 0 && bgzf_close(bgzf) == 0);
    assert(bcf_index_build3(compressed, index_path, 14, 0) == 0);
    duckhts_bcf_index_t index = {0};
    assert(duckhts_bcf_index_load(&index, vcf, compressed, index_path, 0) == 1);
    const char *selectors[] = {NULL, "S2", "^S1", "S3,S1", ""};
    for (int format = 0; format < 3; format++) for (int selection = 0; selection < 5; selection++) {
        duckhts_bcf_scan_t reader = {0};
        duckhts_bcf_samples_t samples = {0};
        duckhts_bcf_gt_span_t spans[3];
        char error[512];
        assert(duckhts_bcf_scan_open(&reader, format ? compressed : plain, &index, 0,
            DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
        assert(duckhts_bcf_samples_build(&samples, reader.hdr, selectors[selection], error, sizeof(error)));
        assert(duckhts_bcf_samples_apply(&samples, reader.hdr, error, sizeof(error)));
        reader.raw_samples = &samples;
        reader.raw_gt = samples.count ? spans : NULL;
        /* BGZF getline needs two spare bytes; vcf_parse needs four. The first
         * compressed row therefore resizes after these GT offsets are captured. */
        size_t first_capacity = strlen(records[0]) + 1u;
        reader.line.s = malloc(first_capacity);
        assert(reader.line.s);
        reader.line.m = first_capacity;
        if (format == 2) {
            char *regions[] = {"chrR:10-40", "chrR:20-70", "chrR:10-10"};
            assert(duckhts_bcf_scan_regions(&reader, regions, 3, error, sizeof(error)));
        }
        bcf1_t *record = bcf_init();
        assert(record);
        size_t row = 0;
        int ret;
        while ((ret = duckhts_bcf_scan_next(&reader, record)) >= 0) {
            assert(row < sizeof(expected) / sizeof(*expected) && !reader.raw_gt_error);
            if (!row && format) assert(reader.line.m > first_capacity);
            assert(record->n_sample == (unsigned)samples.count);
            for (int i = 0; i < samples.count; i++) {
                const char *text = expected[row][samples.indices[i]];
                assert((spans[i].offset == SIZE_MAX) == (text == NULL));
                if (text) {
                    assert(spans[i].offset < reader.line.m && spans[i].length < reader.line.m - spans[i].offset);
                    assert(spans[i].length == strlen(text));
                    assert(memcmp(reader.line.s + spans[i].offset, text, spans[i].length) == 0);
                }
            }
            if (row == 7) assert(reader.line.m >= line_length + 4u);
            row++;
        }
        assert(ret == -1 && row == sizeof(expected) / sizeof(*expected));
        bcf_destroy(record);
        duckhts_bcf_scan_close(&reader);
        duckhts_bcf_samples_destroy(&samples);
    }
    duckhts_bcf_index_destroy(&index);
    const char *bad[] = {
        "chrR\t10\tbad\tA\tC\t.\t.\t.\tGT:GT\t0/1:1/1\t0/0:1/1\t1/1:0/0\n",
        "chrR\t10\tbad\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\n",
        "chrR\t10\tbad\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t1/1\n",
        "chrR\t10\tbad\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\tX\n",
        "chrR\t10\tbad\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n"
    };
    const char *reasons[] = {"duplicate FORMAT/GT", "sample column count", "sample columns", NULL,
                            "embedded NUL"};
    for (size_t i = 0; i < sizeof(bad) / sizeof(*bad); i++) {
        file = fopen(plain, "wb");
        assert(file && fputs(header, file) >= 0);
        if (i == 4) {
            size_t prefix = (size_t)(strstr(bad[i], "0/1") - bad[i]) + 1u;
            assert(fwrite(bad[i], 1, prefix, file) == prefix && fputc('\0', file) != EOF);
            assert(fputs(bad[i] + prefix, file) >= 0);
        } else assert(fputs(bad[i], file) >= 0);
        assert(fclose(file) == 0);
        duckhts_bcf_scan_t reader = {0};
        duckhts_bcf_samples_t samples = {0};
        duckhts_bcf_gt_span_t spans[3];
        char error[512];
        assert(duckhts_bcf_scan_open(&reader, plain, NULL, 0,
            DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
        assert(duckhts_bcf_samples_build(&samples, reader.hdr, NULL, error, sizeof(error)));
        reader.raw_samples = &samples;
        reader.raw_gt = spans;
        bcf1_t *record = bcf_init();
        assert(record && duckhts_bcf_scan_next(&reader, record) < -1);
        if (reasons[i]) assert(reader.raw_gt_error && strstr(reader.raw_gt_error, reasons[i]));
        else assert(!reader.raw_gt_error); /* Invalid GT grammar is rejected by HTSlib. */
        bcf_destroy(record);
        duckhts_bcf_scan_close(&reader);
        duckhts_bcf_samples_destroy(&samples);
    }
    const char *text_edges[] = {
        "chrR\t10\tedge\tA\tC,G\t.\t.\t.\tGT\t|0|1\t/0|1/2\t.\r\n",
        "chrR\t10\tedge\tA\tC,G\t.\t.\t.\tGT\t|0|1\t/0|1/2\t.",
        "chrR\t10\tedge\tA\tC,G\t.\t.\t.\tGT\t\t/0|1/2\t.\n"
    };
    const char *edge_gt[] = {"|0|1", "/0|1/2", "."};
    for (int edge = 0; edge < 3; edge++) {
        file = fopen(plain, "wb");
        assert(file);
        for (const char *p = header; *p; p++) {
            if (!edge && *p == '\n') assert(fputc('\r', file) != EOF);
            assert(fputc(*p, file) != EOF);
        }
        assert(fputs(text_edges[edge], file) >= 0 && fclose(file) == 0);
        file = fopen(plain, "rb");
        bgzf = bgzf_open(compressed, "w");
        assert(file && bgzf);
        while ((bytes = fread(buffer, 1, sizeof(buffer), file)) != 0)
            assert(bgzf_write(bgzf, buffer, bytes) == (ssize_t)bytes);
        assert(!ferror(file) && fclose(file) == 0 && bgzf_close(bgzf) == 0);
        assert(bcf_index_build3(compressed, index_path, 14, 0) == 0);
        assert(duckhts_bcf_index_load(&index, vcf, compressed, index_path, 0) == 1);
        for (int format = 0; format < 3; format++) {
            const char *path = format ? compressed : plain;
            duckhts_bcf_scan_t reader = {0};
            duckhts_bcf_samples_t samples = {0};
            duckhts_bcf_gt_span_t spans[3];
            char error[512];
            assert(duckhts_bcf_scan_open(&reader, path, &index, 0,
                DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
            assert(duckhts_bcf_samples_build(&samples, reader.hdr, NULL, error, sizeof(error)));
            reader.raw_samples = &samples;
            reader.raw_gt = spans;
            if (format == 2) {
                char *region = "chrR:10-10";
                assert(duckhts_bcf_scan_regions(&reader, &region, 1, error, sizeof(error)));
            }
            htsFile *oracle = hts_open(path, "r");
            assert(oracle);
            bcf_hdr_t *oracle_header = bcf_hdr_read(oracle);
            bcf1_t *oracle_record = bcf_init(), *record = bcf_init();
            assert(oracle_header && oracle_record && record);
            int expected_status = bcf_read(oracle, oracle_header, oracle_record);
            int status = duckhts_bcf_scan_next(&reader, record);
            assert(status == expected_status && !reader.raw_gt_error);
            if (edge == 2) assert(status < -1); /* HTSlib rejects the explicit empty GT. */
            else {
                assert(status >= 0 && record->n_sample == 3);
                for (int i = 0; i < 3; i++) {
                    assert(spans[i].offset != SIZE_MAX && spans[i].length == strlen(edge_gt[i]));
                    assert(memcmp(reader.line.s + spans[i].offset, edge_gt[i], spans[i].length) == 0);
                }
                assert(duckhts_bcf_scan_next(&reader, record) == -1);
            }
            bcf_destroy(oracle_record);
            bcf_hdr_destroy(oracle_header);
            assert(hts_close(oracle) == 0);
            bcf_destroy(record);
            duckhts_bcf_scan_close(&reader);
            duckhts_bcf_samples_destroy(&samples);
        }
        duckhts_bcf_index_destroy(&index);
    }
    duckhts_bcf_scan_t reader = {0};
    duckhts_bcf_samples_t samples = {0};
    char error[512];
    assert(duckhts_bcf_scan_open(&reader, "test/data/bcf_scan_contigs.bcf", NULL, 0,
        DUCKHTS_HTS_IO_PROFILE_METADATA, "test", error, sizeof(error)));
    reader.raw_samples = &samples;
    bcf1_t *record = bcf_init();
    assert(record && duckhts_bcf_scan_next(&reader, record) < -1);
    assert(reader.raw_gt_error && strstr(reader.raw_gt_error, "BCF does not retain original GT text"));
    bcf_destroy(record);
    duckhts_bcf_scan_close(&reader);
    assert(unlink(plain) == 0 && unlink(compressed) == 0 && unlink(index_path) == 0);
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
    assert(values.samples == 2 && values.gt_stride == 4 && values.ps_stride == 1);
    assert(memcmp(values.gt, gt, sizeof(gt)) == 0 && memcmp(values.ps, ps, sizeof(ps)) == 0);
    assert(duckhts_bcf_genotype_ploidy(values.gt, 4) == 4);
    assert(duckhts_bcf_genotype_ploidy(values.gt + 4, 4) == 1);
    assert(duckhts_bcf_genotype_has_alt(values.gt, 4));
    int32_t no_alt[] = {0, 1, bcf_gt_unphased(0)};
    assert(!duckhts_bcf_genotype_has_alt(no_alt, 3));
    const int32_t scalar_ps[] = {INT32_MAX, bcf_int32_missing, bcf_int32_vector_end};
    for (int stride = 1; stride <= 8; stride++) for (int first = 0; first < 3; first++)
    for (int second = 0; second < 3; second++) {
        int32_t padded[16];
        for (int i = 0; i < 2 * stride; i++) padded[i] = bcf_int32_vector_end;
        padded[0] = scalar_ps[first];
        padded[stride] = scalar_ps[second];
        assert(bcf_update_format_int32(header, record, "PS", padded, 2 * stride) == 0);
        assert(duckhts_bcf_genotypes_decode(&values, header, record,
            DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
        assert(values.gt_stride == 4 && values.ps_stride == stride);
        assert(values.ps[0] == scalar_ps[first] && values.ps[stride] == scalar_ps[second]);
        assert(memcmp(values.ps, padded, 2u * (size_t)stride * sizeof(*padded)) == 0);
        assert(memcmp(values.gt, gt, sizeof(gt)) == 0);
    }
    int32_t wide_ps[] = {10, 11, 20, bcf_int32_vector_end};
    assert(bcf_update_format_int32(header, record, "PS", wide_ps, 4) == 0);
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(strstr(error, "FORMAT/PS has 2 values for header Number=1 at chrG:1 sample S1"));
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(values.gt_stride == 4 && !values.ps_stride);
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_WARN, error, sizeof(error)));
    assert(values.gt_stride == 4 && !values.ps_stride);
    const char *strings[] = {"10", "20"};
    assert(bcf_update_format_string(header, record, "PS", strings, 2) == 0);
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(strstr(error, "encoded BCF type CHAR"));
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(values.gt_stride == 4 && !values.ps_stride);
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
    assert(values.gt_stride == 0 && values.ps_stride == 1);
    bcf_fmt_t *format = bcf_get_fmt(header, record, "GT");
    int width = format->n;
    format->n = INT32_MAX; // Refuse signed count multiplication before HTSlib decoding.
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(strstr(error, "exceeds the supported decoded-value capacity"));
    format->n = width;
    assert(bcf_update_format_string(header, record, "GT", strings, 2) == 0);
    assert(!duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_ERROR, error, sizeof(error)));
    assert(strstr(error, "encoded BCF type CHAR"));
    assert(duckhts_bcf_genotypes_decode(&values, header, record, DUCKHTS_BCF_DECODE_NULL, error, sizeof(error)));
    assert(values.gt_stride == 0 && values.ps_stride == 1);
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

static void format_values(void) {
    bcf_hdr_t *header = bcf_hdr_init("w");
    bcf1_t *record = bcf_init();
    assert(header && record);
    assert(bcf_hdr_append(header, "##contig=<ID=chrF,length=100>") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=VI,Number=.,Type=Integer,Description=\"Integers\">") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=VF,Number=.,Type=Float,Description=\"Floats\">") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=ST,Number=.,Type=String,Description=\"Strings\">") == 0);
    assert(bcf_hdr_add_sample(header, "S1") == 0 && bcf_hdr_add_sample(header, "S2") == 0);
    assert(bcf_hdr_add_sample(header, NULL) == 0 && bcf_hdr_sync(header) == 0);
    record->rid = 0;
    record->pos = 9;
    assert(bcf_update_alleles_str(header, record, "A,C,G") == 0);
    duckhts_bcf_format_t integers = {0}, floats = {0}, strings = {0};
    char error[512];
    int32_t *expected_i = NULL;
    float *expected_f = NULL;
    int capacity_i = 0, capacity_f = 0;
    uint32_t rng = 217;
    for (int trial = 0; trial < 1024; trial++) {
        rng = rng * 1664525u + 1013904223u;
        int width = 1 + (int)(rng % 257);
        int32_t vi[514];
        float vf[514];
        for (int sample = 0; sample < 2; sample++) for (int j = 0; j < width; j++) {
            size_t i = (size_t)sample * width + j;
            rng = rng * 1664525u + 1013904223u;
            vi[i] = (int32_t)(rng % (trial % 3 == 0 ? 120 : trial % 3 == 1 ? 30000 : 1000000));
            vf[i] = (float)vi[i] / 8;
            if (rng % 5 == 0) { vi[i] = bcf_int32_missing; bcf_float_set_missing(vf[i]); }
            if (sample == 1 && j >= width / 2) {
                vi[i] = bcf_int32_vector_end;
                bcf_float_set_vector_end(vf[i]);
            }
        }
        assert(bcf_update_format_int32(header, record, "VI", vi, width * 2) == 0);
        assert(bcf_update_format_float(header, record, "VF", vf, width * 2) == 0);
        int ni = bcf_get_format_int32(header, record, "VI", &expected_i, &capacity_i);
        int nf = bcf_get_format_float(header, record, "VF", &expected_f, &capacity_f);
        integers.loaded = floats.loaded = 0;
        assert(duckhts_bcf_format_decode(&integers, header, record, "VI", BCF_HT_INT,
            DUCKHTS_BCF_DECODE_ERROR, "test", error, sizeof(error)));
        assert(duckhts_bcf_format_decode(&floats, header, record, "VF", BCF_HT_REAL,
            DUCKHTS_BCF_DECODE_ERROR, "test", error, sizeof(error)));
        assert(ni == 2 * width && nf == ni && integers.count == ni && floats.count == nf);
        assert(integers.stride == width && floats.stride == width);
        assert(memcmp(integers.data, expected_i, (size_t)ni * sizeof(*expected_i)) == 0);
        assert(memcmp(floats.data, expected_f, (size_t)nf * sizeof(*expected_f)) == 0);
        assert(memcmp(integers.data, vi, (size_t)ni * sizeof(*vi)) == 0);
        assert(memcmp(floats.data, vf, (size_t)nf * sizeof(*vf)) == 0);
        const char *text[] = {trial % 2 ? "a,.,b" : "longer,.,last", "."};
        assert(bcf_update_format_string(header, record, "ST", text, 2) == 0);
        strings.loaded = 0;
        assert(duckhts_bcf_format_decode(&strings, header, record, "ST", BCF_HT_STR,
            DUCKHTS_BCF_DECODE_ERROR, "test", error, sizeof(error)));
        assert(strcmp(strings.strings[0], text[0]) == 0 && strcmp(strings.strings[1], text[1]) == 0);
        void *retained = integers.data;
        assert(bcf_update_format_int32(header, record, "VI", NULL, 0) == 0);
        assert(duckhts_bcf_format_decode(&integers, header, record, "VI", BCF_HT_INT,
            DUCKHTS_BCF_DECODE_ERROR, "test", error, sizeof(error)));
        assert(integers.count == ni && integers.data == retained); // Cache lasts until explicit reset.
        integers.loaded = 0;
        assert(duckhts_bcf_format_decode(&integers, header, record, "VI", BCF_HT_INT,
            DUCKHTS_BCF_DECODE_ERROR, "test", error, sizeof(error)));
        assert(!integers.count && !integers.stride && integers.data == retained);
    }
    const char *bad[] = {"bad", "payload"};
    assert(bcf_update_format_string(header, record, "VI", bad, 2) == 0);
    integers.loaded = 0;
    assert(!duckhts_bcf_format_decode(&integers, header, record, "VI", BCF_HT_INT,
        DUCKHTS_BCF_DECODE_ERROR, "test", error, sizeof(error)));
    assert(strstr(error, "encoded BCF type CHAR") && !integers.count);
    assert(duckhts_bcf_format_decode(&integers, header, record, "VI", BCF_HT_INT,
        DUCKHTS_BCF_DECODE_NULL, "test", error, sizeof(error)));
    assert(!integers.count);
    bcf_fmt_t *format = bcf_get_fmt(header, record, "ST");
    int width = format->n;
    format->n = INT32_MAX;
    strings.loaded = 0;
    assert(!duckhts_bcf_format_decode(&strings, header, record, "ST", BCF_HT_STR,
        DUCKHTS_BCF_DECODE_NULL, "test", error, sizeof(error)));
    assert(strstr(error, "decoded-value capacity") && !strings.count);
    format->n = width;
    free(expected_i);
    free(expected_f);
    duckhts_bcf_format_destroy(&integers);
    duckhts_bcf_format_destroy(&floats);
    duckhts_bcf_format_destroy(&strings);
    duckhts_bcf_format_destroy(&strings);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
}

static void numeric_scalars(void) {
    bcf_hdr_t *header = bcf_hdr_init("w");
    bcf1_t *record = bcf_init();
    assert(header && record);
    assert(bcf_hdr_append(header, "##contig=<ID=chrS,length=100>") == 0);
    assert(bcf_hdr_append(header, "##INFO=<ID=SI,Number=1,Type=Integer,Description=\"Scalar\">") == 0);
    assert(bcf_hdr_append(header, "##INFO=<ID=SF,Number=1,Type=Float,Description=\"Scalar\">") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=SI,Number=1,Type=Integer,Description=\"Scalar\">") == 0);
    assert(bcf_hdr_append(header, "##FORMAT=<ID=SF,Number=1,Type=Float,Description=\"Scalar\">") == 0);
    assert(bcf_hdr_add_sample(header, "S1") == 0 && bcf_hdr_add_sample(header, "S2") == 0);
    assert(bcf_hdr_add_sample(header, NULL) == 0 && bcf_hdr_sync(header) == 0);
    record->rid = 0;
    record->pos = 9;
    assert(bcf_update_alleles_str(header, record, "A,C") == 0);
    duckhts_bcf_format_t decoded[2] = {{0}, {0}};
    const char *tags[] = {"SI", "SF"};
    const int types[] = {BCF_HT_INT, BCF_HT_REAL};
    const int magnitudes[] = {1, 300, 100000};
    char error[512];
    enum htsLogLevel log_level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_ERROR);
    for (int scale = 0; scale < 3; scale++) for (int first = 0; first <= 3; first++)
    for (int second = 0; second <= 3; second++) for (int missing = 0; missing < 8; missing++) {
        int32_t integers[6];
        float floats[6];
        for (int sample = 0; sample < 2; sample++) for (int slot = 0; slot < 3; slot++) {
            int i = sample * 3 + slot;
            integers[i] = magnitudes[scale] + i;
            floats[i] = (float)integers[i] + 0.5f;
            if (missing & (1 << slot)) {
                integers[i] = bcf_int32_missing;
                bcf_float_set_missing(floats[i]);
            }
            if (slot >= (sample ? second : first)) {
                integers[i] = bcf_int32_vector_end;
                bcf_float_set_vector_end(floats[i]);
            }
        }
        assert(bcf_update_format_int32(header, record, "SI", integers, 6) == 0);
        assert(bcf_update_format_float(header, record, "SF", floats, 6) == 0);
        for (int type = 0; type < 2; type++) {
            const void *data = type ? (const void *)floats : (const void *)integers;
            int id = bcf_hdr_id2int(header, BCF_DT_ID, tags[type]);
            assert(duckhts_bcf_check_scalar_count(header, record, DUCKHTS_BCF_FIELD_INFO,
                id, types[type], data, 3, "test", error, sizeof(error)) == (first <= 1));
            if (first > 1) assert(strstr(error, "INFO/") && strstr(error, "Number=1 at chrS:10"));
            for (int policy = DUCKHTS_BCF_DECODE_NULL; policy <= DUCKHTS_BCF_DECODE_ERROR; policy++) {
                decoded[type].loaded = 0;
                int valid = first <= 1 && second <= 1;
                int ok = duckhts_bcf_format_decode(&decoded[type], header, record, tags[type], types[type],
                    policy, "test", error, sizeof(error));
                assert(ok == (valid || policy != DUCKHTS_BCF_DECODE_ERROR));
                assert(decoded[type].count == (valid ? 6 : 0));
                assert(decoded[type].stride == (valid ? 3 : 0));
                assert(decoded[type].loaded == ok);
                if (valid) assert(memcmp(decoded[type].data, data, 6 * sizeof(int32_t)) == 0);
                else assert(strstr(error, "FORMAT/") && strstr(error, "Number=1 at chrS:10") &&
                            strstr(error, first > 1 ? "sample S1" : "sample S2"));
            }
        }
    }
    hts_set_log_level(log_level);
    duckhts_bcf_format_destroy(&decoded[0]);
    duckhts_bcf_format_destroy(&decoded[1]);
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
            if (ret == -5) assert(strstr(error, "FORMAT/XX exceeds the supported decoded-value capacity at chr1:10"));
        }
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
    raw_genotypes(argv[1]);
    genotype_values();
    format_values();
    numeric_scalars();
    literal_contigs();
    puts("BCF scanner: CSI/TBI, shifted offsets, 7200 concurrent exact-row scans, selected raw GT and decode errors: OK");
    return 0;
}
