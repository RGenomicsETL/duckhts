#ifndef DUCKHTS_DUCKDB_FS_HFILE_H
#define DUCKHTS_DUCKDB_FS_HFILE_H

/*
 * htslib reads through the DuckDB file system.
 *
 * htslib opens a file by name and offers a scheme handler no other argument,
 * while a DuckDB file system belongs to one connection: it carries that
 * connection's secrets, settings and registered file systems. A route joins
 * the two. It holds the client context of one table function bind and gives
 * each path a name of the form
 *
 *     duckdbfs:<route id>:<path>
 *
 * that htslib opens through the handler registered here. Every name htslib
 * derives from it, such as an index or a .fai beside the file, keeps the
 * prefix and so reads through the same connection.
 *
 * A route is read-only. Opening a routed name for writing fails with EROFS.
 */

#include <stdbool.h>
#include <stddef.h>

#include "duckdb_extension.h"

#define DUCKHTS_FILE_SYSTEM_OPTION "duckhts_file_system"

typedef struct duckhts_fs_route duckhts_fs_route_t;

/*
 * Registers the duckdbfs scheme with htslib when the extension loads.
 * Returns false when htslib cannot create its table of schemes.
 */
bool register_duckdb_fs_hfile_backend(void);

/*
 * Registers the duckhts_file_system setting: 'htslib' (the default) keeps
 * htslib's own file access, 'duckdb' routes reads through the DuckDB file
 * system. Returns false when DuckDB rejects the registration.
 */
bool duckhts_register_file_system_option(duckdb_connection connection);

/*
 * Reads the duckhts_file_system setting of the binding connection. Sets
 * *route to a new route for 'duckdb' and to NULL for 'htslib'. On failure it
 * sets the bind error, prefixed with function_name, and returns false.
 * The caller releases the route with duckhts_fs_route_destroy.
 */
bool duckhts_fs_route_bind(duckdb_bind_info info, const char *function_name,
                           duckhts_fs_route_t **route);

/*
 * Ends the route. Files it already opened stay readable until htslib closes
 * them; a later open of one of its names fails with ENOENT.
 */
void duckhts_fs_route_destroy(duckhts_fs_route_t *route);

/*
 * The name htslib opens to read path: the routed name, or a copy of path
 * when route is NULL. Returns NULL for a NULL path and when memory runs out.
 * Release it with duckdb_free.
 */
char *duckhts_fs_route_path(const duckhts_fs_route_t *route, const char *path);

/*
 * Copies DuckDB's message for the most recent failed open, read, seek or
 * close of path on this route. Returns false, leaving message empty, when
 * route is NULL or its most recent failure concerns another path.
 */
bool duckhts_fs_route_error(duckhts_fs_route_t *route, const char *path,
                            char *message, size_t capacity);

#endif
