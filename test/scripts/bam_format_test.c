/* No DuckDB headers/runtime: test borrowed-record formatters with real HTSlib.
 * The printf substitution is local to this test translation unit. */
#include <htslib/sam.h>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static unsigned calls, fail_at;
static int checked_printf(char *dst, size_t size, const char *format, ...) {
    if (++calls == fail_at) return -1;
    va_list args;
    va_start(args, format);
    int n = vsnprintf(dst, size, format, args);
    va_end(args);
    return n;
}
#define snprintf checked_printf
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC poison malloc calloc realloc free strdup asprintf
#endif
#include "bam_format.h"
#undef snprintf

int main(void) {
    unsigned short_buffers = 0, failed_formats = 0, aux_fields = 0;
    char output[4096];
    size_t capacity, length;
    uint32_t cigar[] = {bam_cigar_gen(UINT32_C(268435455), BAM_CMATCH),
                        bam_cigar_gen(1, BAM_CINS), bam_cigar_gen(2, BAM_CSOFT_CLIP)};
    const char *text = "268435455M1I2S";
    assert(duckhts_bam_cigar_capacity(3, &capacity) && capacity == 31);
    assert(duckhts_bam_cigar_format(cigar, 3, output, strlen(text) + 1, &length));
    assert(length == strlen(text) && !strcmp(output, text));
    for (size_t size = 0; size <= strlen(text); size++) {
        memset(output, 0x5a, sizeof(output));
        assert(!duckhts_bam_cigar_format(cigar, 3, output + 1, size, &length));
        assert(length == 0 && output[0] == 0x5a && output[size + 1] == 0x5a);
        short_buffers++;
    }
    /* Decimal transitions and all HTSlib operation codes, independently
     * formatted by libc. Short capacities above cover the writer's failures. */
    unsigned cigar_cases = 0;
    for (uint32_t power = 1; power <= UINT32_C(100000000); power *= 10) {
        for (uint32_t value = power - 1; value <= power + 1; value++) {
            for (uint32_t op = 0; op < 10; op++) {
                uint32_t encoded = bam_cigar_gen(value, op);
                char expected[16];
                int size = snprintf(expected, sizeof(expected), "%u%c", value, bam_cigar_opchr(encoded));
                assert(size > 0 && (size_t)size < sizeof(expected));
                assert(duckhts_bam_cigar_format(&encoded, 1, output, (size_t)size + 1, &length));
                assert(length == (size_t)size && !strcmp(output, expected));
                cigar_cases++;
            }
        }
    }
    fail_at = 0;
    assert(duckhts_bam_cigar_format(NULL, 0, output, 1, &length));
    assert(length == 0 && output[0] == 0);
    cigar[0] = 15;
    assert(!duckhts_bam_cigar_format(cigar, 1, output, sizeof(output), &length));

    samFile *fp = sam_open("test/data/bam_materialize.sam", "r");
    assert(fp);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *record = bam_init1();
    assert(hdr && record);
    int status;
    while ((status = sam_read1(fp, hdr, record)) >= 0) {
        for (uint8_t *aux = bam_aux_first(record); aux; aux = bam_aux_next(record, aux)) {
            const uint8_t *end = record->data + record->l_data;
            assert(duckhts_bam_aux_capacity(aux, end, &capacity));
            assert(capacity < sizeof(output) - 2);
            calls = 0;
            assert(duckhts_bam_aux_format(aux, output, capacity, &length));
            unsigned format_calls = calls;
            size_t exact_size = length + 1;
            kstring_t oracle = KS_INITIALIZE;
            assert(sam_format_aux1(aux - 2, *aux, aux + 1, end, &oracle));
            assert(oracle.l >= 5 && length == oracle.l - 5);
            assert(!memcmp(output, oracle.s + 5, length));
            ks_free(&oracle);
            for (size_t size = 0; size < exact_size; size++) {
                memset(output, 0x5a, sizeof(output));
                assert(!duckhts_bam_aux_format(aux, output + 1, size, &length));
                assert(length == 0 && output[0] == 0x5a && output[size + 1] == 0x5a);
                short_buffers++;
            }
            assert(duckhts_bam_aux_format(aux, output, exact_size, &length));
            for (unsigned i = 1; i <= format_calls; i++) {
                calls = 0;
                fail_at = i;
                assert(!duckhts_bam_aux_format(aux, output, capacity, &length));
                assert(length == 0 && calls == i);
                failed_formats++;
            }
            fail_at = 0;
            aux_fields++;
        }
    }
    assert(status == -1 && aux_fields == 20);
    bam_destroy1(record);
    sam_hdr_destroy(hdr);
    assert(sam_close(fp) == 0);
    uint8_t unterminated[] = {'Z', 'a'};
    uint8_t oversized[] = {'B', 'I', 255, 255, 255, 255, 0, 0, 0, 0};
    assert(!duckhts_bam_aux_capacity(unterminated, unterminated + sizeof(unterminated), &capacity));
    assert(!duckhts_bam_aux_capacity(oversized, oversized + sizeof(oversized), &capacity));
    printf("bounded BAM formatting: %u CIGAR cases, %u AUX fields, %u short buffers, %u injected printf errors: OK\n",
           cigar_cases, aux_fields, short_buffers, failed_formats);
    return 0;
}
