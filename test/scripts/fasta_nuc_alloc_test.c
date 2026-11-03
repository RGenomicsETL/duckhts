#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "duckdb_extension.h"
DUCKDB_EXTENSION_GLOBAL

static bool fail_allocation;
static void *test_malloc(size_t bytes) {
    return fail_allocation ? NULL : malloc(bytes);
}
#define malloc test_malloc
#include "../../src/interval_udf.c"
#undef malloc

static fasta_nuc_init_data_t init_data;
static idx_t output_size;
static char scan_message[1024];

static void *get_init(duckdb_function_info info) {
    (void)info;
    return &init_data;
}
static idx_t vector_size(void) {
    return 2048;
}
static idx_t column_count(duckdb_data_chunk chunk) {
    (void)chunk;
    return 0;
}
static void set_size(duckdb_data_chunk chunk, idx_t size) {
    (void)chunk;
    output_size = size;
}
static void set_error(duckdb_function_info info, const char *message) {
    (void)info;
    snprintf(scan_message, sizeof(scan_message), "%s", message);
}

int main(void) {
    duckdb_function_get_init_data = get_init;
    duckdb_vector_size = vector_size;
    duckdb_data_chunk_get_column_count = column_count;
    duckdb_data_chunk_set_size = set_size;
    duckdb_function_set_error = set_error;
    fasta_nuc_bind_data_t bind = {0};
    bind.mode = FASTA_NUC_MODE_BED;
    bind.fasta_path = "test/data/fixture_ref.fa";
    bind.bed_path = "test/data/reader_io_error.bed";
    init_data.bind = &bind;
    init_data.bed_fp = hts_open(bind.bed_path, "r");
    if (!init_data.bed_fp) return 2;
    output_size = 99;
    fail_allocation = true;
    fasta_nuc_scan(NULL, NULL);
    fail_allocation = false;
    int failed = !strstr(scan_message, "fasta_nuc:") ||
        !strstr(scan_message, "out of memory") ||
        !strstr(scan_message, bind.bed_path) || output_size != 0;
    if (failed) {
        fprintf(stderr, "expected allocation error and empty chunk; message='%s', rows=%llu\n",
                scan_message, (unsigned long long)output_size);
    }
    free(init_data.bed_line.s);
    hts_close(init_data.bed_fp);
    return failed;
}
