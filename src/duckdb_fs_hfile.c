/*
 * htslib hFILE backend over the DuckDB file system. See duckdb_fs_hfile.h for
 * the contract.
 */
#include "duckdb_fs_hfile.h"

#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include "hfile_internal.h" /* hFILE_backend, hFILE_scheme_handler, hfile_init */
#include "htslib/hfile.h"

DUCKDB_EXTENSION_EXTERN

#define DUCKDBFS_SCHEME "duckdbfs"
#define DUCKDBFS_PREFIX DUCKDBFS_SCHEME ":"

struct duckhts_fs_route {
    uint64_t id;
    duckdb_client_context context;
    /* The fields below are guarded by routes_lock. */
    size_t open_files; /* hFILEs that still refer to the route */
    bool ended;        /* duckhts_fs_route_destroy was called */
    char *error_path;
    char *error_message;
    struct duckhts_fs_route *next;
};

typedef struct {
    hFILE base;
    duckdb_file_handle handle;
    duckhts_fs_route_t *route;
    char *path;
} hFILE_duckdbfs;

static pthread_mutex_t routes_lock = PTHREAD_MUTEX_INITIALIZER;
static duckhts_fs_route_t *routes;
static uint64_t next_route_id = 1;

static char *copy_text(const char *text) {
    const size_t length = strlen(text);
    char *copy = duckdb_malloc(length + 1);
    if (copy) {
        memcpy(copy, text, length + 1);
    }
    return copy;
}

/* Frees a route that no hFILE refers to. The caller holds no lock. */
static void route_free(duckhts_fs_route_t *route) {
    duckdb_destroy_client_context(&route->context);
    duckdb_free(route->error_path);
    duckdb_free(route->error_message);
    duckdb_free(route);
}

/* Returns the live route with this id and counts one more open file on it. */
static duckhts_fs_route_t *route_acquire(uint64_t id) {
    duckhts_fs_route_t *found = NULL;
    pthread_mutex_lock(&routes_lock);
    for (duckhts_fs_route_t *route = routes; route; route = route->next) {
        if (route->id == id) {
            route->open_files++;
            found = route;
            break;
        }
    }
    pthread_mutex_unlock(&routes_lock);
    return found;
}

static void route_release(duckhts_fs_route_t *route) {
    pthread_mutex_lock(&routes_lock);
    route->open_files--;
    const bool unused = route->ended && route->open_files == 0;
    pthread_mutex_unlock(&routes_lock);
    if (unused) {
        route_free(route);
    }
}

/* Takes ownership of error, which may be NULL. */
static void route_record_error(duckhts_fs_route_t *route, const char *path, duckdb_error_data error) {
    const char *text = error ? duckdb_error_data_message(error) : NULL;
    char *message = copy_text(text && text[0] ? text : "the DuckDB file system reported no detail");
    char *failed_path = copy_text(path);
    if (error) {
        duckdb_destroy_error_data(&error);
    }
    if (!message || !failed_path) {
        duckdb_free(message);
        duckdb_free(failed_path);
        return;
    }
    pthread_mutex_lock(&routes_lock);
    char *old_message = route->error_message;
    char *old_path = route->error_path;
    route->error_message = message;
    route->error_path = failed_path;
    pthread_mutex_unlock(&routes_lock);
    duckdb_free(old_message);
    duckdb_free(old_path);
}

static ssize_t duckdbfs_read(hFILE *file, void *buffer, size_t nbytes) {
    hFILE_duckdbfs *fp = (hFILE_duckdbfs *)file;
    const int64_t wanted = nbytes > (size_t)INT64_MAX ? INT64_MAX : (int64_t)nbytes;
    const int64_t got = duckdb_file_handle_read(fp->handle, buffer, wanted);
    if (got < 0) {
        route_record_error(fp->route, fp->path, duckdb_file_handle_error_data(fp->handle));
        errno = EIO;
        return -1;
    }
    return (ssize_t)got;
}

static ssize_t duckdbfs_write(hFILE *file, const void *buffer, size_t nbytes) {
    (void)file;
    (void)buffer;
    (void)nbytes;
    errno = EROFS;
    return -1;
}

static off_t duckdbfs_seek(hFILE *file, off_t offset, int whence) {
    hFILE_duckdbfs *fp = (hFILE_duckdbfs *)file;
    int64_t origin = 0;
    if (whence == SEEK_CUR) {
        origin = duckdb_file_handle_tell(fp->handle);
    } else if (whence == SEEK_END) {
        origin = duckdb_file_handle_size(fp->handle);
    } else if (whence != SEEK_SET) {
        errno = EINVAL;
        return -1;
    }
    if (origin < 0) {
        route_record_error(fp->route, fp->path, duckdb_file_handle_error_data(fp->handle));
        errno = EIO;
        return -1;
    }
    const int64_t delta = (int64_t)offset;
    if ((delta > 0 && origin > INT64_MAX - delta) || (delta < 0 && origin < -delta)) {
        errno = EINVAL;
        return -1;
    }
    const int64_t target = origin + delta;
    if (duckdb_file_handle_seek(fp->handle, target) != DuckDBSuccess) {
        route_record_error(fp->route, fp->path, duckdb_file_handle_error_data(fp->handle));
        errno = EIO;
        return -1;
    }
    return (off_t)target;
}

static int duckdbfs_flush(hFILE *file) {
    (void)file;
    return 0;
}

static int duckdbfs_close(hFILE *file) {
    hFILE_duckdbfs *fp = (hFILE_duckdbfs *)file;
    int status = 0;
    if (duckdb_file_handle_close(fp->handle) != DuckDBSuccess) {
        route_record_error(fp->route, fp->path, duckdb_file_handle_error_data(fp->handle));
        status = -1;
    }
    duckdb_destroy_file_handle(&fp->handle);
    route_release(fp->route);
    duckdb_free(fp->path);
    if (status < 0) {
        errno = EIO;
    }
    return status;
}

static const struct hFILE_backend duckdbfs_backend = {
    duckdbfs_read, duckdbfs_write, duckdbfs_seek, duckdbfs_flush, duckdbfs_close,
};

/*
 * Splits "duckdbfs:<id>:<path>" into its route id and path. Returns NULL for
 * a name that does not have that form.
 */
static const char *parse_routed_name(const char *name, uint64_t *id) {
    const size_t prefix_length = sizeof(DUCKDBFS_PREFIX) - 1;
    if (strncmp(name, DUCKDBFS_PREFIX, prefix_length) != 0) {
        return NULL;
    }
    const char *cursor = name + prefix_length;
    uint64_t value = 0;
    if (*cursor < '0' || *cursor > '9') {
        return NULL;
    }
    for (; *cursor >= '0' && *cursor <= '9'; cursor++) {
        const uint64_t digit = (uint64_t)(*cursor - '0');
        if (value > (UINT64_MAX - digit) / 10) {
            return NULL;
        }
        value = value * 10 + digit;
    }
    if (*cursor != ':' || cursor[1] == '\0') {
        return NULL;
    }
    *id = value;
    return cursor + 1;
}

static hFILE *duckdbfs_open(const char *name, const char *mode) {
    uint64_t id = 0;
    const char *path = parse_routed_name(name, &id);
    if (!path) {
        errno = EINVAL;
        return NULL;
    }
    if (strchr(mode, 'w') || strchr(mode, 'a') || strchr(mode, '+')) {
        errno = EROFS;
        return NULL;
    }
    duckhts_fs_route_t *route = route_acquire(id);
    if (!route) {
        errno = ENOENT;
        return NULL;
    }

    int failure = EIO;
    duckdb_file_handle handle = NULL;
    char *path_copy = copy_text(path);
    duckdb_file_open_options options = duckdb_create_file_open_options();
    duckdb_file_system file_system = duckdb_client_context_get_file_system(route->context);
    if (!path_copy || !options || !file_system) {
        failure = ENOMEM;
    } else if (duckdb_file_open_options_set_flag(options, DUCKDB_FILE_FLAG_READ, true) != DuckDBSuccess) {
        failure = EINVAL;
    } else if (duckdb_file_system_open(file_system, path, options, &handle) != DuckDBSuccess) {
        route_record_error(route, path, duckdb_file_system_error_data(file_system));
        handle = NULL;
    }
    /* An open handle does not depend on the wrapper it was opened from. */
    duckdb_destroy_file_open_options(&options);
    duckdb_destroy_file_system(&file_system);

    hFILE_duckdbfs *fp = NULL;
    if (handle) {
        fp = (hFILE_duckdbfs *)hfile_init(sizeof(*fp), mode, 0);
        if (!fp) {
            failure = ENOMEM;
            duckdb_file_handle_close(handle);
            duckdb_destroy_file_handle(&handle);
        }
    }
    if (!fp) {
        duckdb_free(path_copy);
        route_release(route);
        errno = failure;
        return NULL;
    }
    fp->handle = handle;
    fp->route = route;
    fp->path = path_copy;
    fp->base.backend = &duckdbfs_backend;
    return &fp->base;
}

/*
 * htslib looks for the index of a remote file by opening candidate names and
 * for the index of a local file with stat(). Only the first works for a
 * routed name, whatever storage the path names.
 */
static int duckdbfs_is_remote(const char *name) {
    (void)name;
    return 1;
}

static const struct hFILE_scheme_handler duckdbfs_handler = {
    duckdbfs_open, duckdbfs_is_remote, "DuckHTS DuckDB file system", 50, NULL,
};

bool register_duckdb_fs_hfile_backend(void) {
    /* htslib creates its scheme table on first use and cannot add a handler
     * before that. Listing the schemes creates the table. */
    const char *listed[1];
    int n_listed = 0;
    if (hfile_list_schemes(NULL, listed, &n_listed) < 0) {
        return false;
    }
    hfile_add_scheme_handler(DUCKDBFS_SCHEME, &duckdbfs_handler);
    return true;
}

/*
 * Registers one session setting unless the database already has it, which is
 * the case when another connection loaded the extension first. Takes
 * ownership of default_value.
 */
static bool register_option(duckdb_connection connection, const char *name, duckdb_type type,
                            duckdb_value default_value, const char *description) {
    duckdb_client_context context = NULL;
    duckdb_connection_get_client_context(connection, &context);
    duckdb_value existing = context ? duckdb_client_context_get_config_option(context, name, NULL) : NULL;
    duckdb_config_option option = duckdb_create_config_option();
    duckdb_logical_type logical_type = duckdb_create_logical_type(type);
    bool registered = existing != NULL;
    if (!registered && context && option && logical_type && default_value) {
        duckdb_config_option_set_name(option, name);
        duckdb_config_option_set_type(option, logical_type);
        duckdb_config_option_set_default_value(option, default_value);
        duckdb_config_option_set_default_scope(option, DUCKDB_CONFIG_OPTION_SCOPE_SESSION);
        duckdb_config_option_set_description(option, description);
        registered = duckdb_register_config_option(connection, option) == DuckDBSuccess;
    }
    duckdb_destroy_value(&existing);
    duckdb_destroy_value(&default_value);
    duckdb_destroy_logical_type(&logical_type);
    duckdb_destroy_config_option(&option);
    if (context) {
        duckdb_destroy_client_context(&context);
    }
    return registered;
}

bool duckhts_register_file_access_options(duckdb_connection connection) {
    return register_option(connection, DUCKHTS_FILE_SYSTEM_OPTION, DUCKDB_TYPE_VARCHAR,
                           duckdb_create_varchar("htslib"),
                           "File access for DuckHTS readers: 'htslib' (htslib's own I/O) or 'duckdb' "
                           "(the DuckDB file system of the connection, with its secrets and settings)") &&
           register_option(connection, DUCKHTS_REMOTE_BLOCK_SIZE_OPTION, DUCKDB_TYPE_BIGINT,
                           duckdb_create_int64(0),
                           "Bytes htslib asks for in each read of a remote file; 0 lets DuckHTS choose "
                           "by scan type");
}

static void bind_error(duckdb_bind_info info, const char *function_name, const char *problem) {
    char message[256];
    snprintf(message, sizeof(message), "%s: %s", function_name, problem);
    duckdb_bind_set_error(info, message);
}

bool duckhts_remote_block_size_bind(duckdb_bind_info info, const char *function_name, int *block_size) {
    *block_size = 0;
    duckdb_client_context context = NULL;
    duckdb_table_function_get_client_context(info, &context);
    if (!context) {
        bind_error(info, function_name, "DuckDB provided no client context");
        return false;
    }
    duckdb_value setting = duckdb_client_context_get_config_option(context, DUCKHTS_REMOTE_BLOCK_SIZE_OPTION, NULL);
    duckdb_destroy_client_context(&context);
    const bool is_set = setting && !duckdb_is_null_value(setting);
    const int64_t bytes = is_set ? duckdb_get_int64(setting) : -1;
    duckdb_destroy_value(&setting);
    if (bytes != 0 && (bytes < DUCKHTS_REMOTE_BLOCK_SIZE_MIN || bytes > DUCKHTS_REMOTE_BLOCK_SIZE_MAX)) {
        char problem[160];
        snprintf(problem, sizeof(problem), DUCKHTS_REMOTE_BLOCK_SIZE_OPTION " must be 0 or between %d and %d bytes",
                 DUCKHTS_REMOTE_BLOCK_SIZE_MIN, DUCKHTS_REMOTE_BLOCK_SIZE_MAX);
        bind_error(info, function_name, problem);
        return false;
    }
    *block_size = (int)bytes;
    return true;
}

bool duckhts_fs_route_bind(duckdb_bind_info info, const char *function_name, duckhts_fs_route_t **route) {
    *route = NULL;
    duckdb_client_context context = NULL;
    duckdb_table_function_get_client_context(info, &context);
    if (!context) {
        bind_error(info, function_name, "DuckDB provided no client context");
        return false;
    }

    duckdb_value setting = duckdb_client_context_get_config_option(context, DUCKHTS_FILE_SYSTEM_OPTION, NULL);
    char *choice = setting ? duckdb_get_varchar(setting) : NULL;
    duckdb_destroy_value(&setting);
    const bool use_htslib = choice && strcmp(choice, "htslib") == 0;
    const bool use_duckdb = choice && strcmp(choice, "duckdb") == 0;
    duckdb_free(choice);
    if (!use_duckdb) {
        duckdb_destroy_client_context(&context);
        if (!use_htslib) {
            bind_error(info, function_name, DUCKHTS_FILE_SYSTEM_OPTION " must be 'htslib' or 'duckdb'");
        }
        return use_htslib;
    }

    duckhts_fs_route_t *created = duckdb_malloc(sizeof(*created));
    if (!created) {
        duckdb_destroy_client_context(&context);
        bind_error(info, function_name, "out of memory routing reads through the DuckDB file system");
        return false;
    }
    memset(created, 0, sizeof(*created));
    created->context = context;
    pthread_mutex_lock(&routes_lock);
    created->id = next_route_id++;
    created->next = routes;
    routes = created;
    pthread_mutex_unlock(&routes_lock);
    *route = created;
    return true;
}

void duckhts_fs_route_destroy(duckhts_fs_route_t *route) {
    if (!route) {
        return;
    }
    pthread_mutex_lock(&routes_lock);
    for (duckhts_fs_route_t **link = &routes; *link; link = &(*link)->next) {
        if (*link == route) {
            *link = route->next;
            break;
        }
    }
    route->ended = true;
    const bool unused = route->open_files == 0;
    pthread_mutex_unlock(&routes_lock);
    if (unused) {
        route_free(route);
    }
}

char *duckhts_fs_route_path(const duckhts_fs_route_t *route, const char *path) {
    if (!path) {
        return NULL;
    }
    if (!route) {
        return copy_text(path);
    }
    const size_t path_length = strlen(path);
    /* "duckdbfs:" + up to 20 digits + ":" + path + NUL */
    const size_t prefix_capacity = sizeof(DUCKDBFS_PREFIX) + 21;
    if (path_length > SIZE_MAX - prefix_capacity) {
        return NULL;
    }
    char *name = duckdb_malloc(prefix_capacity + path_length);
    if (name) {
        snprintf(name, prefix_capacity + path_length, DUCKDBFS_PREFIX "%" PRIu64 ":%s", route->id, path);
    }
    return name;
}

bool duckhts_fs_route_error(duckhts_fs_route_t *route, const char *path, char *message, size_t capacity) {
    bool found = false;
    if (capacity > 0) {
        message[0] = '\0';
    }
    if (!route || !path || capacity == 0) {
        return false;
    }
    pthread_mutex_lock(&routes_lock);
    if (route->error_path && strcmp(route->error_path, path) == 0) {
        snprintf(message, capacity, "%s", route->error_message);
        found = true;
    }
    pthread_mutex_unlock(&routes_lock);
    return found;
}
