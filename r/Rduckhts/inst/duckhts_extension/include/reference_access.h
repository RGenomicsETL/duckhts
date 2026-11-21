#ifndef DUCKHTS_REFERENCE_ACCESS_H
#define DUCKHTS_REFERENCE_ACCESS_H

/*
 * The indexed FASTA of a scalar function, read with the file access of the
 * connection that runs the query (hts_file_access.h).
 *
 * A reference access is the state of one worker for one query. With 'htslib'
 * its references are the entries of the thread cache (reference_cache.h),
 * which builds a missing .fai and keeps handles across queries. With 'duckdb'
 * the FASTA, its .fai and its .gzi are read through the connection, so the
 * reference access owns those handles itself: at most
 * DUCKHTS_REFERENCE_CACHE_ENTRIES of them, the least recently used closed
 * first, and all of them closed when DuckDB destroys the state. That route
 * builds no index. Either way the caller fetches bases with
 * duckhts_reference_fetch.
 */

#include "duckdb_extension.h"
#include "reference_cache.h"

typedef struct duckhts_reference_access duckhts_reference_access_t;

/*
 * The body of a scalar init hook (duckdb_scalar_function_set_init). Captures
 * the file access of the connection and makes a reference access the state of
 * this worker; DuckDB destroys it. An invalid setting is the init error,
 * prefixed with function_name.
 */
void duckhts_reference_access_scalar_init(duckdb_init_info info, const char *function_name);

/*
 * The reference of these paths, borrowed until the next call with this
 * reference access; do not keep it across scalar invocations. references is
 * duckdb_scalar_function_get_state() of a function whose init hook is the one
 * above. The arguments are those of duckhts_reference_cache_get. On failure
 * returns NULL, *error is that function's diagnostic and *reason is DuckDB's
 * reason as " (reason)", or "" when there is none; both are borrowed and valid
 * until the next call on this thread.
 */
duckhts_reference_entry_t *duckhts_reference_access_get(
    duckhts_reference_access_t *references, const char *fasta_path, const char *fai_path,
    const char *gzi_path, int tune_indexed_remote_io, const char **error, const char **reason);

#endif
