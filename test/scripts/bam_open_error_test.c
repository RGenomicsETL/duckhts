/* Exercise worker initialization when a bound input is no longer available.
 * SQL cannot deterministically remove the file between binding and local init.
 * Run with test/scripts/bam_open_error_test.sh after make release. */
#include "../../src/bam_reader.c"

static bam_bind_data_t test_bind;
static char actual_error[1024];
static unsigned error_count;
static unsigned allocations;

static void *test_malloc(size_t bytes) {
    void *data = malloc(bytes);
    if (data) allocations++;
    return data;
}

static void test_free(void *data) {
    if (data) allocations--;
    free(data);
}

static void *get_bind(duckdb_init_info info) {
    (void)info;
    return &test_bind;
}

static idx_t get_column_count(duckdb_init_info info) {
    (void)info;
    return 0;
}

static void set_error(duckdb_init_info info, const char *message) {
    (void)info;
    error_count++;
    snprintf(actual_error, sizeof(actual_error), "%s", message);
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    test_bind.file_path = argv[1];
    test_bind.scan_sequential = 1;
    duckdb_malloc = test_malloc;
    duckdb_free = test_free;
    duckdb_init_get_bind_data = get_bind;
    duckdb_init_get_column_count = get_column_count;
    duckdb_init_set_error = set_error;

    bam_read_local_init(NULL);

    char expected[1024];
    snprintf(expected, sizeof(expected),
             "read_bam: failed to open SAM/BAM/CRAM file: %s", argv[1]);
    if (error_count != 1 || strcmp(actual_error, expected) != 0 || allocations != 0) {
        fprintf(stderr, "expected: %s\nactual: %s\nerrors: %u, live allocations: %u\n",
                expected, actual_error, error_count, allocations);
        return 1;
    }
    return 0;
}
