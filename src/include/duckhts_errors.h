#ifndef DUCKHTS_ERRORS_H
#define DUCKHTS_ERRORS_H

#include "duckdb_extension.h"
#include "duckhts_attributes.h"

#include <htslib/kstring.h>
#include <stdarg.h>

/* printf-style DuckDB error setters for table-function callbacks. If formatting the
 * message fails, a fixed out-of-memory message is set instead, so the callback still
 * reports an error. Include after DUCKDB_EXTENSION_EXTERN: the setters resolve through
 * the extension's C API table. */
#define DUCKHTS_DEFINE_ERRORF(name, info_type, set_error)                              \
    static inline void name(info_type info, const char *format, ...) DUCKHTS_PRINTF(2, 3); \
    static inline void name(info_type info, const char *format, ...) {                 \
        kstring_t message = KS_INITIALIZE;                                              \
        va_list args;                                                                   \
        va_start(args, format);                                                         \
        int written = kvsprintf(&message, format, args);                                \
        va_end(args);                                                                   \
        set_error(info, written < 0 ? "DuckHTS: out of memory formatting an error" : message.s); \
        ks_free(&message);                                                              \
    }

DUCKHTS_DEFINE_ERRORF(duckhts_bind_errorf, duckdb_bind_info, duckdb_bind_set_error)
DUCKHTS_DEFINE_ERRORF(duckhts_init_errorf, duckdb_init_info, duckdb_init_set_error)
DUCKHTS_DEFINE_ERRORF(duckhts_scan_errorf, duckdb_function_info, duckdb_function_set_error)

#endif
