/*
 * File access for every DuckHTS reader. See hts_file_access.h for the
 * contract.
 */
#include "hts_file_access.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "duckdb_fs_hfile.h"

DUCKDB_EXTENSION_EXTERN

struct duckhts_file_access {
    duckhts_fs_route_t *route; /* NULL: htslib's own file access */
    int remote_block_size;     /* bytes; 0 is the automatic size */
};

/* ================================================================
 * Settings
 * ================================================================ */

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

/*
 * Reads both settings of context into a new file access. Takes ownership of
 * context. On failure it writes the problem, without the function name, to
 * problem and returns NULL.
 */
static duckhts_file_access_t *capture(duckdb_client_context context, char *problem, size_t capacity) {
    if (!context) {
        snprintf(problem, capacity, "DuckDB provided no client context");
        return NULL;
    }

    duckdb_value block_setting =
        duckdb_client_context_get_config_option(context, DUCKHTS_REMOTE_BLOCK_SIZE_OPTION, NULL);
    const bool block_is_set = block_setting && !duckdb_is_null_value(block_setting);
    const int64_t block_bytes = block_is_set ? duckdb_get_int64(block_setting) : -1;
    duckdb_destroy_value(&block_setting);
    if (block_bytes != 0 &&
        (block_bytes < DUCKHTS_REMOTE_BLOCK_SIZE_MIN || block_bytes > DUCKHTS_REMOTE_BLOCK_SIZE_MAX)) {
        snprintf(problem, capacity, DUCKHTS_REMOTE_BLOCK_SIZE_OPTION " must be 0 or between %d and %d bytes",
                 DUCKHTS_REMOTE_BLOCK_SIZE_MIN, DUCKHTS_REMOTE_BLOCK_SIZE_MAX);
        duckdb_destroy_client_context(&context);
        return NULL;
    }

    duckdb_value system_setting =
        duckdb_client_context_get_config_option(context, DUCKHTS_FILE_SYSTEM_OPTION, NULL);
    char *choice = system_setting ? duckdb_get_varchar(system_setting) : NULL;
    duckdb_destroy_value(&system_setting);
    const bool use_htslib = choice && strcmp(choice, "htslib") == 0;
    const bool use_duckdb = choice && strcmp(choice, "duckdb") == 0;
    duckdb_free(choice);
    if (!use_htslib && !use_duckdb) {
        snprintf(problem, capacity, DUCKHTS_FILE_SYSTEM_OPTION " must be 'htslib' or 'duckdb'");
        duckdb_destroy_client_context(&context);
        return NULL;
    }

    duckhts_file_access_t *access = duckdb_malloc(sizeof(*access));
    if (!access) {
        snprintf(problem, capacity, "out of memory capturing the file access settings");
        duckdb_destroy_client_context(&context);
        return NULL;
    }
    access->route = NULL;
    access->remote_block_size = (int)block_bytes;
    if (use_htslib) {
        duckdb_destroy_client_context(&context);
        return access;
    }
    access->route = duckhts_fs_route_create(context);
    if (!access->route) {
        snprintf(problem, capacity, "out of memory routing reads through the DuckDB file system");
        duckdb_free(access);
        return NULL;
    }
    return access;
}

bool duckhts_file_access_bind(duckdb_bind_info info, const char *function_name,
                              duckhts_file_access_t **access) {
    char problem[160];
    char message[256];
    duckdb_client_context context = NULL;
    duckdb_table_function_get_client_context(info, &context);
    *access = capture(context, problem, sizeof(problem));
    if (!*access) {
        snprintf(message, sizeof(message), "%s: %s", function_name, problem);
        duckdb_bind_set_error(info, message);
    }
    return *access != NULL;
}

bool duckhts_file_access_scalar_init(duckdb_init_info info, const char *function_name,
                                     duckhts_file_access_t **access) {
    char problem[160];
    char message[256];
    duckdb_client_context context = NULL;
    duckdb_scalar_function_init_get_client_context(info, &context);
    *access = capture(context, problem, sizeof(problem));
    if (!*access) {
        snprintf(message, sizeof(message), "%s: %s", function_name, problem);
        duckdb_scalar_function_init_set_error(info, message);
    }
    return *access != NULL;
}

void duckhts_file_access_destroy(duckhts_file_access_t *access) {
    if (!access) {
        return;
    }
    duckhts_fs_route_destroy(access->route);
    duckdb_free(access);
}

/* ================================================================
 * Names
 * ================================================================ */

static const duckhts_fs_route_t *route_of(const duckhts_file_access_t *access) {
    return access ? access->route : NULL;
}

bool duckhts_file_access_uses_duckdb(const duckhts_file_access_t *access) {
    return route_of(access) != NULL;
}

static char *copy_span(const char *text, size_t length) {
    char *copy = duckdb_malloc(length + 1);
    if (copy) {
        memcpy(copy, text, length);
        copy[length] = '\0';
    }
    return copy;
}

/*
 * The name htslib opens to read path: a copy of path with 'htslib', the routed
 * name with 'duckdb'. Both halves of "file" HTS_IDX_DELIM "index" are named.
 * Returns NULL for a NULL path and when memory runs out. Release it with
 * duckdb_free.
 */
static char *access_name(const duckhts_file_access_t *access, const char *path) {
    const duckhts_fs_route_t *route = route_of(access);
    const char *delimiter = (route && path) ? strstr(path, HTS_IDX_DELIM) : NULL;
    if (!delimiter) {
        return duckhts_fs_route_path(route, path);
    }

    const size_t delimiter_length = sizeof(HTS_IDX_DELIM) - 1;
    char *file = copy_span(path, (size_t)(delimiter - path));
    char *file_name = file ? duckhts_fs_route_path(route, file) : NULL;
    char *index_name = duckhts_fs_route_path(route, delimiter + delimiter_length);
    char *name = NULL;
    if (file_name && index_name) {
        const size_t file_length = strlen(file_name);
        const size_t index_length = strlen(index_name);
        if (index_length <= SIZE_MAX - file_length - delimiter_length - 1) {
            name = duckdb_malloc(file_length + delimiter_length + index_length + 1);
        }
        if (name) {
            memcpy(name, file_name, file_length);
            memcpy(name + file_length, HTS_IDX_DELIM, delimiter_length);
            memcpy(name + file_length + delimiter_length, index_name, index_length + 1);
        }
    }
    duckdb_free(file);
    duckdb_free(file_name);
    duckdb_free(index_name);
    return name;
}

/* ================================================================
 * Remote tuning
 * ================================================================ */

bool duckhts_file_access_is_remote(const duckhts_file_access_t *access, const char *path) {
    if (!path || path[0] == '\0') {
        return false;
    }
    if (!route_of(access)) {
        return duckhts_is_remote_path(path);
    }
    const char *separator = strstr(path, "://");
    if (!separator || separator == path) {
        return false;
    }
    for (const char *cursor = path; cursor < separator; cursor++) {
        const char c = *cursor;
        const bool scheme_character = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                                      (c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.';
        if (!scheme_character) {
            return false;
        }
    }
    return strncmp(path, "file://", 7) != 0;
}

/* These settings are advisory: htslib ignores an option a format or backend
 * does not honour. */
static void tune_hts(const duckhts_file_access_t *access, htsFile *fp, const char *path,
                     duckhts_hts_io_profile_t profile) {
    if (!duckhts_file_access_is_remote(access, path)) {
        return;
    }
    int block_size = duckhts_remote_block_size_bytes(profile);
    const int cache_size = duckhts_remote_cache_size_bytes(profile);
    if (access && access->remote_block_size > 0) {
        block_size = access->remote_block_size;
    }
    if (block_size > 0) {
        (void)hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, block_size);
    }
    if (cache_size > 0) {
        (void)hts_set_opt(fp, HTS_OPT_CACHE_SIZE, cache_size);
    }
}

static void tune_hfile(const duckhts_file_access_t *access, hFILE *fp, const char *path,
                       duckhts_hts_io_profile_t profile) {
    if (!duckhts_file_access_is_remote(access, path)) {
        return;
    }
    int block_size = duckhts_remote_hfile_block_size_bytes(profile);
    if (access && access->remote_block_size > 0) {
        block_size = access->remote_block_size;
    }
    if (block_size > 0) {
        (void)hfile_set_blksize(fp, (size_t)block_size);
    }
}

/* ================================================================
 * Opens
 * ================================================================ */

htsFile *duckhts_file_access_open(const duckhts_file_access_t *access, const char *path,
                                  duckhts_hts_io_profile_t profile) {
    htsFile *fp = NULL;
    if (!route_of(access)) {
        fp = hts_open(path, "r");
    } else {
        char *name = access_name(access, path);
        if (name) {
            fp = hts_open(name, "r");
            duckdb_free(name);
        }
    }
    if (fp) {
        tune_hts(access, fp, path, profile);
    }
    return fp;
}

int duckhts_file_access_sr_add_reader(const duckhts_file_access_t *access, bcf_srs_t *readers,
                                      const char *path) {
    if (!route_of(access)) {
        return bcf_sr_add_reader(readers, path);
    }
    char *name = access_name(access, path);
    const int status = name ? bcf_sr_add_reader(readers, name) : 0;
    duckdb_free(name);
    return status;
}

hFILE *duckhts_file_access_hopen(const duckhts_file_access_t *access, const char *path,
                                 duckhts_hts_io_profile_t profile) {
    hFILE *fp = NULL;
    if (!route_of(access)) {
        fp = hopen(path, "r");
    } else {
        char *name = access_name(access, path);
        if (name) {
            fp = hopen(name, "r");
            duckdb_free(name);
        }
    }
    if (fp) {
        tune_hfile(access, fp, path, profile);
    }
    return fp;
}

char **duckhts_file_access_readlist(void *access, const char *path, int *n) {
    const duckhts_file_access_t *file_access = access;
    if (!route_of(file_access)) {
        return hts_readlist(path, 1, n);
    }
    hFILE *hfile = duckhts_file_access_hopen(file_access, path, DUCKHTS_HTS_IO_PROFILE_METADATA);
    BGZF *fp = hfile ? bgzf_hopen(hfile, "r") : NULL;
    if (!fp) {
        if (hfile) {
            hclose_abruptly(hfile);
        }
        return NULL;
    }

    char **lines = NULL;
    size_t count = 0;
    size_t capacity = 0;
    kstring_t line = KS_INITIALIZE;
    int status;
    bool failed = false;
    while ((status = bgzf_getline(fp, '\n', &line)) >= 0) {
        if (line.l == 0) {
            continue;
        }
        if (count == capacity) {
            const size_t grown = capacity ? capacity * 2 : 16;
            char **larger = grown <= INT_MAX ? realloc(lines, grown * sizeof(*lines)) : NULL;
            if (!larger) {
                failed = true;
                break;
            }
            lines = larger;
            capacity = grown;
        }
        lines[count] = strdup(line.s);
        if (!lines[count]) {
            failed = true;
            break;
        }
        count++;
    }
    failed = failed || status < -1;
    failed = bgzf_close(fp) != 0 || failed;
    free(line.s);
    if (failed) {
        for (size_t i = 0; i < count; i++) {
            free(lines[i]);
        }
        free(lines);
        return NULL;
    }
    if (!lines) {
        /* An empty list is a valid result; hts_readlist returns a non-NULL array for it. */
        lines = malloc(sizeof(*lines));
        if (!lines) {
            return NULL;
        }
    }
    *n = (int)count;
    return lines;
}

/* ================================================================
 * Indexes
 * ================================================================ */

typedef enum {
    INDEX_OF_BAM, /* BAM and SAM */
    INDEX_OF_CRAM,
    INDEX_OF_BCF,
    INDEX_OF_TABIX
} index_family_t;

#define INDEX_SUFFIX_CAPACITY sizeof(".crai")

/*
 * The index beside a routed file, as a routed name, or NULL when none opens.
 * Candidates follow the order of htslib's hts_idx_check_local: each suffix
 * appended to the path, then in place of the extension of its last component.
 * For a URL, as in htslib's haddextension, the suffix goes before a query or
 * fragment, so a presigned https://host/a.bam?sig=... looks for
 * https://host/a.bam.bai?sig=...; s3 URLs keep a '#' in the key.
 */
static char *find_routed_index(const duckhts_fs_route_t *route, const char *path,
                               index_family_t family) {
    static const char *const bam_suffixes[] = {".csi", ".bai"};
    static const char *const cram_suffixes[] = {".crai"};
    static const char *const bcf_suffixes[] = {".csi"};
    static const char *const tabix_suffixes[] = {".csi", ".tbi"};
    const char *const *suffixes = bam_suffixes;
    size_t n_suffixes = 2;
    if (family == INDEX_OF_CRAM) {
        suffixes = cram_suffixes;
        n_suffixes = 1;
    } else if (family == INDEX_OF_BCF) {
        suffixes = bcf_suffixes;
        n_suffixes = 1;
    } else if (family == INDEX_OF_TABIX) {
        suffixes = tabix_suffixes;
    }

    const size_t path_length = strlen(path);
    /* path_end is where a suffix goes: the end of the path, or for a URL the
     * start of its query or fragment. */
    size_t path_end = path_length;
    if (strstr(path, "://")) {
        const bool s3 = strncmp(path, "s3://", 5) == 0 || strncmp(path, "s3+http://", 10) == 0 ||
                        strncmp(path, "s3+https://", 11) == 0;
        path_end = strcspn(path, s3 ? "?" : "?#");
    }
    size_t stem_length = path_end;
    for (size_t i = path_end; i > 0; i--) {
        if (path[i - 1] == '.') {
            stem_length = i - 1;
            break;
        }
        if (path[i - 1] == '/') {
            break;
        }
    }
    if (path_length > SIZE_MAX - INDEX_SUFFIX_CAPACITY) {
        return NULL;
    }
    char *candidate = duckdb_malloc(path_length + INDEX_SUFFIX_CAPACITY);
    if (!candidate) {
        return NULL;
    }

    char *found = NULL;
    for (size_t s = 0; s < n_suffixes && !found; s++) {
        const size_t bases[] = {path_end, stem_length};
        const size_t n_bases = stem_length == path_end ? 1 : 2;
        const size_t suffix_length = strlen(suffixes[s]);
        for (size_t b = 0; b < n_bases && !found; b++) {
            memcpy(candidate, path, bases[b]);
            memcpy(candidate + bases[b], suffixes[s], suffix_length);
            strcpy(candidate + bases[b] + suffix_length, path + path_end);
            char *name = duckhts_fs_route_path(route, candidate);
            hFILE *probe = name ? hopen(name, "r") : NULL;
            if (probe && hclose(probe) == 0) {
                found = name;
            } else {
                duckdb_free(name);
            }
        }
    }
    duckdb_free(candidate);
    return found;
}

/*
 * The routed names of a file and of its index. A NULL index_path takes the
 * index half of "file" HTS_IDX_DELIM "index" when path has one and otherwise
 * searches beside the file. Returns false, with both names NULL, when there
 * is no index or memory runs out. The caller frees both with duckdb_free.
 */
static bool routed_index_names(const duckhts_fs_route_t *route, const char *path,
                               const char *index_path, index_family_t family, char **file_name,
                               char **index_name) {
    *file_name = NULL;
    *index_name = NULL;
    const char *delimiter = strstr(path, HTS_IDX_DELIM);
    char *file = delimiter ? copy_span(path, (size_t)(delimiter - path)) : NULL;
    const char *file_path = delimiter ? file : path;
    if (file_path) {
        *file_name = duckhts_fs_route_path(route, file_path);
        if (index_path) {
            *index_name = duckhts_fs_route_path(route, index_path);
        } else if (delimiter) {
            *index_name = duckhts_fs_route_path(route, delimiter + sizeof(HTS_IDX_DELIM) - 1);
        } else {
            *index_name = find_routed_index(route, file_path, family);
        }
    }
    duckdb_free(file);
    if (*file_name && *index_name) {
        return true;
    }
    duckdb_free(*file_name);
    duckdb_free(*index_name);
    *file_name = NULL;
    *index_name = NULL;
    return false;
}

hts_idx_t *duckhts_file_access_sam_index_load(const duckhts_file_access_t *access, htsFile *fp,
                                              const char *path, const char *index_path, int flags) {
    const duckhts_fs_route_t *route = route_of(access);
    if (!route) {
        return sam_index_load3(fp, path, index_path, flags);
    }
    const index_family_t family = hts_get_format(fp)->format == cram ? INDEX_OF_CRAM : INDEX_OF_BAM;
    char *file_name = NULL;
    char *index_name = NULL;
    hts_idx_t *index = NULL;
    if (routed_index_names(route, path, index_path, family, &file_name, &index_name)) {
        index = sam_index_load3(fp, file_name, index_name, flags & ~HTS_IDX_SAVE_REMOTE);
    }
    duckdb_free(file_name);
    duckdb_free(index_name);
    return index;
}

hts_idx_t *duckhts_file_access_bcf_index_load(const duckhts_file_access_t *access,
                                              const char *path, const char *index_path, int flags) {
    const duckhts_fs_route_t *route = route_of(access);
    if (!route) {
        return bcf_index_load3(path, index_path, flags);
    }
    char *file_name = NULL;
    char *index_name = NULL;
    hts_idx_t *index = NULL;
    if (routed_index_names(route, path, index_path, INDEX_OF_BCF, &file_name, &index_name)) {
        index = bcf_index_load3(file_name, index_name, flags & ~HTS_IDX_SAVE_REMOTE);
    }
    duckdb_free(file_name);
    duckdb_free(index_name);
    return index;
}

tbx_t *duckhts_file_access_tbx_index_load(const duckhts_file_access_t *access, const char *path,
                                          const char *index_path, int flags) {
    const duckhts_fs_route_t *route = route_of(access);
    if (!route) {
        return tbx_index_load3(path, index_path, flags);
    }
    char *file_name = NULL;
    char *index_name = NULL;
    tbx_t *index = NULL;
    if (routed_index_names(route, path, index_path, INDEX_OF_TABIX, &file_name, &index_name)) {
        index = tbx_index_load3(file_name, index_name, flags & ~HTS_IDX_SAVE_REMOTE);
    }
    duckdb_free(file_name);
    duckdb_free(index_name);
    return index;
}

faidx_t *duckhts_file_access_fai_load(const duckhts_file_access_t *access, const char *path,
                                      const char *fai_path, const char *gzi_path, int flags,
                                      enum fai_format_options format,
                                      duckhts_hts_io_profile_t profile) {
    const duckhts_fs_route_t *route = route_of(access);
    faidx_t *fai = NULL;
    if (!route) {
        fai = fai_load3_format(path, fai_path, gzi_path, flags, format);
    } else {
        char *name = duckhts_fs_route_path(route, path);
        char *fai_name = duckhts_fs_route_path(route, fai_path);
        char *gzi_name = duckhts_fs_route_path(route, gzi_path);
        if (name && (fai_name || !fai_path) && (gzi_name || !gzi_path)) {
            fai = fai_load3_format(name, fai_name, gzi_name, flags & ~FAI_CREATE, format);
        }
        duckdb_free(name);
        duckdb_free(fai_name);
        duckdb_free(gzi_name);
    }
    if (fai && duckhts_file_access_is_remote(access, path)) {
        const int cache_size = duckhts_remote_cache_size_bytes(profile);
        if (cache_size > 0) {
            fai_set_cache_size(fai, cache_size);
        }
    }
    return fai;
}

/* ================================================================
 * Alignment reference
 * ================================================================ */

int duckhts_file_access_set_cram_reference(const duckhts_file_access_t *access, htsFile *fp,
                                           const char *locator) {
    if (!route_of(access)) {
        return hts_set_opt(fp, CRAM_OPT_REFERENCE, locator);
    }
    char *name = access_name(access, locator);
    const int status = name ? hts_set_opt(fp, CRAM_OPT_REFERENCE, name) : -1;
    duckdb_free(name);
    return status;
}

int duckhts_file_access_set_fai_filename(const duckhts_file_access_t *access, htsFile *fp,
                                         const char *locator) {
    if (!route_of(access)) {
        return hts_set_fai_filename(fp, locator);
    }
    char *name = access_name(access, locator);
    const int status = name ? hts_set_fai_filename(fp, name) : -1;
    duckdb_free(name);
    return status;
}

/* ================================================================
 * Errors
 * ================================================================ */

const char *duckhts_file_access_reason(const duckhts_file_access_t *access, const char *path) {
    static _Thread_local char reason[512];
    char detail[480];
    reason[0] = '\0';
    if (access && access->route && duckhts_fs_route_error(access->route, path, detail, sizeof(detail))) {
        snprintf(reason, sizeof(reason), " (%s)", detail);
    }
    return reason;
}
