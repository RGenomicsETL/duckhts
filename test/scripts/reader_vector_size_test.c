/* Exercise the production scan loops with a non-default DuckDB vector size.
 * Build/run with test/scripts/reader_vector_size_test.sh after make release. */
#include <stdio.h>
#include <stdlib.h>

#ifdef TEST_TABIX
#include "../../src/tabix_reader.c"
#else
#include "../../src/interval_udf.c"
#endif

static void *test_bind;
static void *test_init;
static idx_t emitted;
static const idx_t capacity = 17;

static idx_t vector_capacity(void) {
    return capacity;
}

static void *bind_data(duckdb_function_info info) {
    (void)info;
    return test_bind;
}

static void *init_data(duckdb_function_info info) {
    (void)info;
    return test_init;
}

static idx_t column_count(duckdb_data_chunk output) {
    (void)output;
    return 0;
}

static void chunk_size(duckdb_data_chunk output, idx_t size) {
    (void)output;
    emitted = size;
}

static void scan_error(duckdb_function_info info, const char *message) {
    (void)info;
    fprintf(stderr, "scan error: %s\n", message);
    exit(1);
}

static int check_chunk(const char *name) {
    if (emitted == capacity) return 0;
    fprintf(stderr, "%s emitted %llu rows for vector capacity %llu\n", name,
            (unsigned long long)emitted, (unsigned long long)capacity);
    return 1;
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    duckdb_vector_size = vector_capacity;
    duckdb_function_get_bind_data = bind_data;
    duckdb_function_get_init_data = init_data;
    duckdb_data_chunk_get_column_count = column_count;
    duckdb_data_chunk_set_size = chunk_size;
    duckdb_function_set_error = scan_error;
    int failed = 0;
#ifdef TEST_TABIX
    tabix_bind_data_t bind = {0};
    tabix_init_data_t init = {0};
    test_bind = &bind;
    test_init = &init;
    init.fp = hts_open(argv[1], "r");
    if (!init.fp) return 2;
    tabix_scan(NULL, NULL);
    failed = check_chunk("read_tabix");
    free(init.line.s);
    hts_close(init.fp);
#else
    for (int tolerant = 0; tolerant < 2; tolerant++) {
        bed_init_data_t init = {0};
        init.fp = hts_open(argv[1], "r");
        if (!init.fp) return 2;
        if (tolerant) {
            read_bed_rows_tolerant(NULL, NULL, &init, NULL, 0);
        } else {
            read_bed_rows_error(NULL, NULL, &init, NULL, 0);
        }
        failed |= check_chunk(tolerant ? "read_bed tolerant" : "read_bed error");
        free(init.line.s);
        hts_close(init.fp);
    }
    fasta_nuc_bind_data_t bind = {0};
    fasta_nuc_init_data_t init = {0};
    bind.mode = FASTA_NUC_MODE_BINS;
    bind.bin_width = 1;
    init.bind = &bind;
    init.fai = fai_load3("test/data/fixture_ref.fa", NULL, NULL, 0);
    if (!init.fai) return 2;
    init.current_seq_end_exclusive = faidx_seq_len64(init.fai, faidx_iseq(init.fai, 0));
    test_bind = &bind;
    test_init = &init;
    fasta_nuc_scan(NULL, NULL);
    failed |= check_chunk("fasta_nuc");
    fai_destroy(init.fai);
#endif
    return failed;
}
