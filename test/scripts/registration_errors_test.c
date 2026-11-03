/* Fault injection for macro hash diagnostics and DuckDB result ownership. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../src/extension_registration.c"
DUCKDB_EXTENSION_GLOBAL

typedef enum {
    HASH_OK,
    HASH_PREPARE_ERROR,
    HASH_BIND_ERROR,
    HASH_EXECUTE_ERROR,
    HASH_BAD_LENGTH
} hash_failure_t;

static hash_failure_t failure;
static int result_destroys;
static int statement_destroys;
static int chunk_destroys;
static int error_calls;
static int failures;
static char error_message[256];
static char hash[65];
static duckdb_string_t hash_value;

static void
mock_error(duckdb_extension_info info, const char *message)
{
    (void)info;
    error_calls++;
    snprintf(error_message, sizeof(error_message), "%s", message);
}

static duckdb_state
mock_prepare(duckdb_connection connection, const char *sql,
    duckdb_prepared_statement *statement)
{
    (void)connection;
    (void)sql;
    *statement = (duckdb_prepared_statement)&failure;
    return failure == HASH_PREPARE_ERROR ? DuckDBError : DuckDBSuccess;
}

static const char *
mock_prepare_error(duckdb_prepared_statement statement)
{
    (void)statement;
    return "injected hash prepare/bind diagnostic";
}

static void
mock_destroy_prepare(duckdb_prepared_statement *statement)
{
    statement_destroys++;
    *statement = NULL;
}

static duckdb_state
mock_bind(duckdb_prepared_statement statement, idx_t index, const char *value, idx_t length)
{
    (void)statement;
    (void)index;
    (void)value;
    (void)length;
    return failure == HASH_BIND_ERROR ? DuckDBError : DuckDBSuccess;
}

static duckdb_state
mock_execute(duckdb_prepared_statement statement, duckdb_result *result)
{
    (void)statement;
    (void)result;
    return failure == HASH_EXECUTE_ERROR ? DuckDBError : DuckDBSuccess;
}

static const char *
mock_result_error(duckdb_result *result)
{
    (void)result;
    return "injected hash execute diagnostic";
}

static void
mock_destroy_result(duckdb_result *result)
{
    (void)result;
    result_destroys++;
}

static duckdb_data_chunk
mock_fetch(duckdb_result result)
{
    (void)result;
    return (duckdb_data_chunk)&hash_value;
}

static duckdb_vector
mock_vector(duckdb_data_chunk chunk, idx_t column)
{
    (void)column;
    return (duckdb_vector)chunk;
}

static void *
mock_data(duckdb_vector vector)
{
    return vector;
}

static uint32_t
mock_length(duckdb_string_t value)
{
    (void)value;
    return failure == HASH_BAD_LENGTH ? 1u : 64u;
}

static const char *
mock_string(duckdb_string_t *value)
{
    (void)value;
    return hash;
}

static void
mock_destroy_chunk(duckdb_data_chunk *chunk)
{
    chunk_destroys++;
    *chunk = NULL;
}

static void
check_case(hash_failure_t selected, const char *diagnostic,
    int expected_results, int expected_chunks)
{
    struct duckdb_extension_access access = {0};
    access.set_error = mock_error;
    duckhts_registration_t registration = {.access = &access};
    failure = selected;
    result_destroys = 0;
    statement_destroys = 0;
    chunk_destroys = 0;
    error_calls = 0;
    error_message[0] = '\0';
    bool ok = macro_digest(&registration);
    if (ok != (selected == HASH_OK) ||
        (diagnostic != NULL && strstr(error_message, diagnostic) == NULL) ||
        error_calls != (selected == HASH_OK ? 0 : 1) ||
        result_destroys != expected_results || statement_destroys != 1 ||
        chunk_destroys != expected_chunks ||
        (selected == HASH_OK && strcmp(definitions_sha256, hash) != 0)) {
        fprintf(stderr, "case %d: ok=%d error=%s error_calls=%d "
            "result_destroys=%d (expected %d) statement_destroys=%d chunk_destroys=%d\n",
            selected, ok, error_message, error_calls, result_destroys, expected_results,
            statement_destroys, chunk_destroys);
        failures++;
    }
}

int
main(void)
{
    duckdb_prepare = mock_prepare;
    duckdb_prepare_error = mock_prepare_error;
    duckdb_destroy_prepare = mock_destroy_prepare;
    duckdb_bind_varchar_length = mock_bind;
    duckdb_execute_prepared = mock_execute;
    duckdb_result_error = mock_result_error;
    duckdb_destroy_result = mock_destroy_result;
    duckdb_fetch_chunk = mock_fetch;
    duckdb_data_chunk_get_vector = mock_vector;
    duckdb_vector_get_data = mock_data;
    duckdb_string_t_length = mock_length;
    duckdb_string_t_data = mock_string;
    duckdb_destroy_data_chunk = mock_destroy_chunk;
    macros[0].sql = "CREATE OR REPLACE MACRO test_macro() AS 1";
    macro_count = 1;
    memset(hash, 'a', 64u);
    check_case(HASH_PREPARE_ERROR, "injected hash prepare/bind diagnostic", 0, 0);
    check_case(HASH_BIND_ERROR, "injected hash prepare/bind diagnostic", 0, 0);
    check_case(HASH_EXECUTE_ERROR, "injected hash execute diagnostic", 1, 0);
    check_case(HASH_BAD_LENGTH, "Cannot hash macro definitions", 1, 1);
    check_case(HASH_OK, NULL, 1, 1);
    return failures == 0 ? 0 : 1;
}
