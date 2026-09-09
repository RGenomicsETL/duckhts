/* Isolated HTSlib INFO allocation failures, independent of DuckDB's allocator.
 * Link the static HTSlib with --wrap=realloc; arm only the decoded buffer. */
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <htslib/vcf.h>

static void *failed_pointer;
static int armed, failures;
void *__real_realloc(void *pointer, size_t bytes);
void *__wrap_realloc(void *pointer, size_t bytes) {
    if (armed && pointer == failed_pointer) {
        armed = 0;
        failures++;
        return NULL;
    }
    return __real_realloc(pointer, bytes);
}

static void check(int type, int growing) {
    const char *names[] = {"Integer", "Float", "String", "Integer"};
    const int types[] = {BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR, BCF_HT_LONG};
    int32_t integers[] = {1, 2, 3, 4, 5, 6, 7, 8};
    int64_t long_integers[] = {1, 2, 3, 4, 5, 6, 7, 8};
    float floats[] = {1.25f, 2.25f, 3.25f, 4.25f, 5.25f, 6.25f, 7.25f, 8.25f};
    const char *strings[] = {"ab", "abcdefgh"};
    bcf_hdr_t *header = bcf_hdr_init("w");
    bcf1_t *record = bcf_init();
    assert(header && record);
    char line[128];
    snprintf(line, sizeof(line), "##INFO=<ID=V,Number=.,Type=%s,Description=\"Values\">", names[type]);
    assert(bcf_hdr_append(header, "##contig=<ID=chrO,length=100>") == 0);
    assert(bcf_hdr_append(header, line) == 0 && bcf_hdr_sync(header) == 0);
    record->rid = 0;
    assert(bcf_update_alleles_str(header, record, "A,C") == 0);
    void *values = NULL;
    int capacity = 0;
    for (int stage = !growing; stage < 2; stage++) {
        int length = stage ? 8 : 2;
        if (type == 0 || type == 3) assert(bcf_update_info_int32(header, record, "V", integers, length) == 0);
        if (type == 1) assert(bcf_update_info_float(header, record, "V", floats, length) == 0);
        if (type == 2) assert(bcf_update_info_string(header, record, "V", strings[stage]) == 0);
        assert(bcf_unpack(record, BCF_UN_INFO) == 0);
        if (stage) {
            void *previous = values;
            int previous_capacity = capacity;
            unsigned char previous_bytes[16];
            size_t bytes = previous ? (type == 2 ? 3 : 2 * (type == 3 ? sizeof(int64_t) : sizeof(int32_t))) : 0;
            if (bytes) memcpy(previous_bytes, previous, bytes);
            failed_pointer = values;
            armed = 1;
            assert(bcf_get_info_values(header, record, "V", &values, &capacity, types[type]) == -4);
            assert(failures == 1 && !armed);
            assert(values == previous && capacity == previous_capacity);
            if (bytes) assert(memcmp(values, previous_bytes, bytes) == 0);
        }
        assert(bcf_get_info_values(header, record, "V", &values, &capacity, types[type]) == length);
        const void *expected = type == 0 ? (const void *)integers : type == 1 ? (const void *)floats :
            type == 2 ? (const void *)strings[stage] : (const void *)long_integers;
        assert(memcmp(values, expected, type == 2 ? (size_t)length + 1 :
                           (size_t)length * (type == 3 ? sizeof(int64_t) : sizeof(int32_t))) == 0);
        if (type == 2 && stage) {
            bcf_info_t *info = bcf_get_info(header, record, "V");
            int old_length = info->len, old_capacity = capacity;
            void *old_values = values;
            info->len = INT_MAX; /* Terminator cannot fit the API's signed capacity. */
            assert(bcf_get_info_values(header, record, "V", &values, &capacity, types[type]) == -4);
            assert(values == old_values && capacity == old_capacity);
            info->len = old_length;
        }
    }
    free(values);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
}

int main(void) {
    struct rlimit core_limit = {0, 0};
    assert(setrlimit(RLIMIT_CORE, &core_limit) == 0);
    int failed = 0;
    for (int type = 0; type < 4; type++) for (int growing = 0; growing < 2; growing++) {
        pid_t child = fork();
        assert(child >= 0);
        if (!child) { check(type, growing); exit(0); }
        int status;
        assert(waitpid(child, &status, 0) == child);
        if (!WIFEXITED(status) || WEXITSTATUS(status)) {
            fprintf(stderr, "INFO allocation failure: type=%d growing=%d status=%d\n", type, growing, status);
            failed++;
        }
    }
    if (!failed) puts("INFO allocation failures: 8 injected errors, retained ownership and exact retry: OK");
    return failed != 0;
}
