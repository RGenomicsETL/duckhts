/* Callback-level check: SQL suppresses errored results, so inspect chunk size. */
#include <stdio.h>
#include "../../src/bigwig_reader.c"
DUCKDB_EXTENSION_GLOBAL

static bigwig_global_data_t global;
static bigwig_local_data_t local;
static int fail_next;
static int error_calls;
static char error_message[256];
static idx_t output_size;
static uint32_t output_start[8];

int
bwInit(size_t size)
{
    (void)size;
    abort();
}

bigWigFile_t *
bwOpen(const char *path, CURLcode (*callback)(CURL *), const char *mode)
{
    (void)path;
    (void)callback;
    (void)mode;
    abort();
}

void
bwIteratorDestroy(bwOverlapIterator_t *iterator)
{
    (void)iterator;
}

bwOverlapIterator_t *
bwIteratorNext(bwOverlapIterator_t *iterator)
{
    if (fail_next)
        return NULL;
    iterator->data = NULL;
    return iterator;
}

bwOverlapIterator_t *
bwOverlappingIntervalsIterator(bigWigFile_t *file, const char *chrom,
    uint32_t start, uint32_t end, uint32_t blocks)
{
    (void)file;
    (void)chrom;
    (void)start;
    (void)end;
    (void)blocks;
    return NULL;
}

static void *
mock_global(duckdb_function_info info)
{
    (void)info;
    return &global;
}

static void *
mock_local(duckdb_function_info info)
{
    (void)info;
    return &local;
}

static void
mock_error(duckdb_function_info info, const char *message)
{
    (void)info;
    error_calls++;
    snprintf(error_message, sizeof(error_message), "%s", message);
}

static void
mock_size(duckdb_data_chunk chunk, idx_t size)
{
    (void)chunk;
    output_size = size;
}

static idx_t
mock_vector_size(void)
{
    return 8u;
}

static duckdb_vector
mock_vector(duckdb_data_chunk chunk, idx_t column)
{
    (void)chunk;
    (void)column;
    return (duckdb_vector)output_start;
}

static void *
mock_data(duckdb_vector vector)
{
    return vector;
}

static int
check_case(int next_error, int claim_error, int projected)
{
    char *chrom = "chr1";
    chromList_t chromosomes = {.nKeys = 1, .chrom = &chrom};
    bigWigFile_t file = {.cl = &chromosomes};
    uint32_t start = 10u;
    bwOverlappingIntervals_t intervals = {.l = 1, .start = &start};
    bwOverlapIterator_t iterator = {.data = &intervals, .intervals = &intervals};
    bigwig_job_t job = {.tid = 0};
    bigwig_bind_data_t bind = {.blocks_per_iteration = 1};
    idx_t column = BIGWIG_COL_START;
    global = (bigwig_global_data_t){.jobs = &job, .job_count = 1,
        .next_job = claim_error ? 0 : 1};
    local = (bigwig_local_data_t){.file = &file, .iterator = &iterator,
        .bind = &bind, .column_count = projected ? 1 : 0, .column_ids = &column,
        .owner_start0 = 0, .owner_end0 = 100};
    fail_next = next_error;
    error_calls = 0;
    error_message[0] = '\0';
    output_size = 99;
    output_start[0] = 0;
    read_bigwig_scan(NULL, NULL);
    int expected_error = next_error || claim_error;
    const char *expected_message = next_error
        ? "read_bigwig: failed while reading an indexed range"
        : "read_bigwig: failed to initialize an indexed range";
    if (error_calls != expected_error || output_size != (expected_error ? 0u : 1u) ||
        !local.done || (projected && output_start[0] != start) ||
        (expected_error && strcmp(error_message, expected_message) != 0)) {
        fprintf(stderr, "next_error=%d claim_error=%d projected=%d: "
            "errors=%d rows=%llu done=%d message=%s\n", next_error, claim_error,
            projected, error_calls, (unsigned long long)output_size, local.done,
            error_message);
        return 1;
    }
    output_size = 99;
    read_bigwig_scan(NULL, NULL);
    return output_size != 0u || error_calls != expected_error;
}

int
main(void)
{
    duckdb_function_get_init_data = mock_global;
    duckdb_function_get_local_init_data = mock_local;
    duckdb_function_set_error = mock_error;
    duckdb_data_chunk_set_size = mock_size;
    duckdb_vector_size = mock_vector_size;
    duckdb_data_chunk_get_vector = mock_vector;
    duckdb_vector_get_data = mock_data;
    int failures = 0;
    for (int projected = 0; projected <= 1; projected++) {
        failures += check_case(1, 0, projected);
        failures += check_case(0, 1, projected);
        failures += check_case(0, 0, projected);
    }
    return failures == 0 ? 0 : 1;
}
