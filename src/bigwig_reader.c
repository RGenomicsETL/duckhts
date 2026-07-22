/**
 * Indexed BigWig interval reader backed by libBigWig.
 *
 * read_bigwig(path, region := NULL, blocks_per_iteration := 64)
 *   -> CHROM VARCHAR, START0 UINTEGER, END0 UINTEGER, VALUE FLOAT
 *
 * Coordinates are the stored zero-based, half-open BigWig intervals.  A
 * requested range filters overlapping records but does not clip their stored
 * coordinates.  Each scan owns its file handle and iterator.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include "bigWig.h"

#define DUCKHTS_BIGWIG_DEFAULT_BLOCKS 64u
#define DUCKHTS_BIGWIG_MAX_BLOCKS 1048576u
#define DUCKHTS_BIGWIG_MAX_REGIONS 1048576u
#if defined(__EMSCRIPTEN__)
#define DUCKHTS_BIGWIG_REMOTE_BUFFER (128u * 1024u)
#else
#define DUCKHTS_BIGWIG_REMOTE_BUFFER (256u * 1024u)
#endif

enum {
    BIGWIG_COL_CHROM = 0,
    BIGWIG_COL_START,
    BIGWIG_COL_END,
    BIGWIG_COL_VALUE,
    BIGWIG_COL_COUNT
};

typedef struct {
    char *path;
    char *region;
    uint32_t blocks_per_iteration;
} bigwig_bind_data_t;

typedef struct {
    uint32_t tid;
    uint32_t start0;
    uint32_t end0;
    uint32_t owner_start0;
    uint32_t owner_end0;
} bigwig_job_t;

typedef struct {
    bigWigFile_t *seed_file;
    bigwig_job_t *jobs;
    uint32_t job_count;
    uint32_t next_job;
    int seed_claimed;
} bigwig_global_data_t;

typedef struct {
    bigWigFile_t *file;
    bwOverlapIterator_t *iterator;
    bigwig_bind_data_t *bind;
    idx_t *column_ids;
    idx_t column_count;
    uint32_t interval_index;
    uint32_t owner_start0;
    uint32_t owner_end0;
    int done;
} bigwig_local_data_t;

static pthread_once_t bigwig_once = PTHREAD_ONCE_INIT;
static int bigwig_init_status = 1;

static void
duckhts_bigwig_init_once(void)
{
    bigwig_init_status = bwInit(DUCKHTS_BIGWIG_REMOTE_BUFFER);
}

static int
duckhts_bigwig_init(void)
{
    if (pthread_once(&bigwig_once, duckhts_bigwig_init_once) != 0)
        return 0;
    return bigwig_init_status == 0;
}

static void
bigwig_bind_destroy(void *data)
{
    bigwig_bind_data_t *bind = data;

    if (bind == NULL)
        return;
    duckdb_free(bind->path);
    duckdb_free(bind->region);
    duckdb_free(bind);
}

static void
bigwig_global_destroy(void *data)
{
    bigwig_global_data_t *global = data;

    if (global == NULL)
        return;
    bwClose(global->seed_file);
    duckdb_free(global->jobs);
    duckdb_free(global);
}

static void
bigwig_local_destroy(void *data)
{
    bigwig_local_data_t *local = data;

    if (local == NULL)
        return;
    bwIteratorDestroy(local->iterator);
    bwClose(local->file);
    duckdb_free(local->column_ids);
    duckdb_free(local);
}

static int
bigwig_open_checked(const char *path, bigWigFile_t **out)
{
    bigWigFile_t *file;

    if (!duckhts_bigwig_init())
        return 0;
    file = bwOpen(path, NULL, "r");
    if (file == NULL)
        return 0;
    *out = file;
    return 1;
}

static void
bigwig_add_columns(duckdb_bind_info info)
{
    duckdb_logical_type varchar_type;
    duckdb_logical_type uinteger_type;
    duckdb_logical_type float_type;

    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    float_type = duckdb_create_logical_type(DUCKDB_TYPE_FLOAT);
    duckdb_bind_add_result_column(info, "CHROM", varchar_type);
    duckdb_bind_add_result_column(info, "START0", uinteger_type);
    duckdb_bind_add_result_column(info, "END0", uinteger_type);
    duckdb_bind_add_result_column(info, "VALUE", float_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&uinteger_type);
    duckdb_destroy_logical_type(&float_type);
}

static int
bigwig_name_to_id(void *header, const char *name)
{
    uint32_t tid;

    tid = bwGetTid((bigWigFile_t *)header, name);
    return tid == UINT32_MAX || tid > INT32_MAX ? -1 : (int)tid;
}

static void
bigwig_region_tokens_destroy(char **tokens, size_t count)
{
    size_t i;

    if (tokens == NULL)
        return;
    for (i = 0u; i < count; i++)
        free(tokens[i]);
    free(tokens);
}

static int
bigwig_region_compare(const void *left, const void *right)
{
    const hts_reglist_t *a = left;
    const hts_reglist_t *b = right;

    if (a->tid < b->tid)
        return -1;
    return a->tid > b->tid;
}

/* Parse the public comma-separated region syntax with htslib itself.  The
 * second htslib pass groups by reference and merges overlapping or adjacent
 * ranges, which gives this reader the same once-only multi-region contract as
 * the indexed HTS readers. */
static int
bigwig_parse_regions(const char *spec, bigWigFile_t *file,
    hts_reglist_t **out_regions, int *out_count)
{
    const char *cursor, *next, *token_end;
    char **tokens;
    size_t count, capacity;
    hts_reglist_t *regions;
    int region_count;

    *out_regions = NULL;
    *out_count = 0;
    tokens = NULL;
    count = 0u;
    capacity = 0u;
    cursor = spec;
    while (*cursor != '\0') {
        hts_pos_t start0, end0;
        size_t length;
        int tid;

        tid = -1;
        next = hts_parse_region(cursor, &tid, &start0, &end0,
            bigwig_name_to_id, file, HTS_PARSE_LIST);
        if (next == NULL || tid < 0 || next == cursor)
            goto invalid;
        token_end = next;
        if (token_end[-1] == ',')
            token_end--;
        if (token_end == cursor || count >= DUCKHTS_BIGWIG_MAX_REGIONS)
            goto invalid;
        if (count == capacity) {
            char **grown;
            size_t next_capacity;

            next_capacity = capacity == 0u ? 8u : capacity * 2u;
            if (next_capacity > DUCKHTS_BIGWIG_MAX_REGIONS)
                next_capacity = DUCKHTS_BIGWIG_MAX_REGIONS;
            grown = realloc(tokens, next_capacity * sizeof(*tokens));
            if (grown == NULL)
                goto oom;
            tokens = grown;
            capacity = next_capacity;
        }
        length = (size_t)(token_end - cursor);
        tokens[count] = malloc(length + 1u);
        if (tokens[count] == NULL)
            goto oom;
        memcpy(tokens[count], cursor, length);
        tokens[count][length] = '\0';
        count++;
        if (next[-1] == ',' && *next == '\0')
            goto invalid;
        cursor = next;
    }
    if (count == 0u || count > INT_MAX)
        goto invalid;
    region_count = 0;
    regions = hts_reglist_create(tokens, (int)count, &region_count,
        file, bigwig_name_to_id);
    bigwig_region_tokens_destroy(tokens, count);
    if (regions == NULL || region_count < 1)
        return 0;
    qsort(regions, (size_t)region_count, sizeof(*regions),
        bigwig_region_compare);
    *out_regions = regions;
    *out_count = region_count;
    return 1;

invalid:
    bigwig_region_tokens_destroy(tokens, count);
    return 0;
oom:
    bigwig_region_tokens_destroy(tokens, count);
    return -1;
}

static void
read_bigwig_bind(duckdb_bind_info info)
{
    duckdb_value path_value, region_value, blocks_value;
    bigwig_bind_data_t *bind;
    char *path, *region;
    int64_t blocks;

    path_value = duckdb_bind_get_parameter(info, 0);
    path = duckdb_get_varchar(path_value);
    duckdb_destroy_value(&path_value);
    if (path == NULL || path[0] == '\0') {
        duckdb_bind_set_error(info, "read_bigwig requires a file path or URL");
        duckdb_free(path);
        return;
    }

    region = NULL;
    region_value = duckdb_bind_get_named_parameter(info, "region");
    if (region_value != NULL && !duckdb_is_null_value(region_value))
        region = duckdb_get_varchar(region_value);
    if (region_value != NULL)
        duckdb_destroy_value(&region_value);
    if (region != NULL && region[0] == '\0') {
        duckdb_bind_set_error(info, "read_bigwig: region must be non-empty");
        duckdb_free(path);
        duckdb_free(region);
        return;
    }

    bind = duckdb_malloc(sizeof(*bind));
    if (bind == NULL) {
        duckdb_bind_set_error(info, "read_bigwig: out of memory");
        duckdb_free(path);
        duckdb_free(region);
        return;
    }
    memset(bind, 0, sizeof(*bind));
    bind->path = path;
    bind->region = region;
    bind->blocks_per_iteration = DUCKHTS_BIGWIG_DEFAULT_BLOCKS;

    blocks_value = duckdb_bind_get_named_parameter(info,
        "blocks_per_iteration");
    if (blocks_value != NULL && !duckdb_is_null_value(blocks_value)) {
        blocks = duckdb_get_int64(blocks_value);
        if (blocks < 1 || blocks > (int64_t)DUCKHTS_BIGWIG_MAX_BLOCKS) {
            duckdb_destroy_value(&blocks_value);
            duckdb_bind_set_error(info,
                "read_bigwig: blocks_per_iteration must be between 1 and 1048576");
            bigwig_bind_destroy(bind);
            return;
        }
        bind->blocks_per_iteration = (uint32_t)blocks;
    }
    if (blocks_value != NULL)
        duckdb_destroy_value(&blocks_value);

    bigwig_add_columns(info);
    duckdb_bind_set_bind_data(info, bind, bigwig_bind_destroy);
}

static int
bigwig_build_jobs(bigwig_bind_data_t *bind, bigWigFile_t *file,
    bigwig_global_data_t *global)
{
    hts_reglist_t *regions;
    int region_count, parse_status;
    uint64_t capacity;
    uint32_t count;
    int i;

    if (bind->region == NULL) {
        int64_t tid;

        if (file->cl->nKeys < 0 ||
            (uint64_t)file->cl->nKeys > DUCKHTS_BIGWIG_MAX_REGIONS)
            return 0;
        capacity = (uint64_t)file->cl->nKeys;
        global->jobs = capacity == 0u ? NULL :
            duckdb_malloc((size_t)capacity * sizeof(*global->jobs));
        if (capacity != 0u && global->jobs == NULL)
            return -1;
        count = 0u;
        for (tid = 0; tid < file->cl->nKeys; tid++) {
            if (file->cl->len[tid] == 0u)
                continue;
            global->jobs[count].tid = (uint32_t)tid;
            global->jobs[count].start0 = 0u;
            global->jobs[count].end0 = file->cl->len[tid];
            global->jobs[count].owner_start0 = 0u;
            global->jobs[count].owner_end0 = file->cl->len[tid];
            count++;
        }
        global->job_count = count;
        return 1;
    }

    regions = NULL;
    region_count = 0;
    parse_status = bigwig_parse_regions(bind->region, file, &regions,
        &region_count);
    if (parse_status <= 0)
        return parse_status;
    capacity = 0u;
    for (i = 0; i < region_count; i++) {
        if ((uint64_t)regions[i].count >
            DUCKHTS_BIGWIG_MAX_REGIONS - capacity) {
            hts_reglist_free(regions, region_count);
            return 0;
        }
        capacity += regions[i].count;
    }
    if (capacity > DUCKHTS_BIGWIG_MAX_REGIONS) {
        hts_reglist_free(regions, region_count);
        return 0;
    }
    global->jobs = capacity == 0u ? NULL :
        duckdb_malloc((size_t)capacity * sizeof(*global->jobs));
    if (capacity != 0u && global->jobs == NULL) {
        hts_reglist_free(regions, region_count);
        return -1;
    }
    count = 0u;
    for (i = 0; i < region_count; i++) {
        hts_reglist_t *region = &regions[i];
        uint32_t chrom_length = file->cl->len[region->tid];
        uint32_t previous_end0 = 0u;
        uint32_t j;

        for (j = 0u; j < region->count; j++) {
            hts_pair_pos_t *interval = &region->intervals[j];
            uint32_t start0, end0;

            start0 = interval->beg < 0 ? 0u :
                (uint64_t)interval->beg > chrom_length ? chrom_length :
                (uint32_t)interval->beg;
            end0 = interval->end < 0 ? 0u :
                (uint64_t)interval->end > chrom_length ? chrom_length :
                (uint32_t)interval->end;
            if (start0 < end0) {
                global->jobs[count].tid = (uint32_t)region->tid;
                global->jobs[count].start0 = start0;
                global->jobs[count].end0 = end0;
                global->jobs[count].owner_start0 = previous_end0;
                global->jobs[count].owner_end0 = end0;
                count++;
            }
            previous_end0 = end0;
        }
    }
    hts_reglist_free(regions, region_count);
    global->job_count = count;
    return 1;
}

static void
read_bigwig_global_init(duckdb_init_info info)
{
    bigwig_bind_data_t *bind;
    bigwig_global_data_t *global;
    int jobs_status;
    idx_t max_threads;

    bind = duckdb_init_get_bind_data(info);
    global = duckdb_malloc(sizeof(*global));
    if (global == NULL) {
        duckdb_init_set_error(info, "read_bigwig: out of memory");
        return;
    }
    memset(global, 0, sizeof(*global));
    if (!bigwig_open_checked(bind->path, &global->seed_file)) {
        duckdb_init_set_error(info,
            "read_bigwig: failed to open a valid BigWig file");
        bigwig_global_destroy(global);
        return;
    }
    jobs_status = bigwig_build_jobs(bind, global->seed_file, global);
    if (jobs_status <= 0) {
        duckdb_init_set_error(info, jobs_status < 0
            ? "read_bigwig: out of memory while preparing scan ranges"
            : "read_bigwig: invalid or unknown region list");
        bigwig_global_destroy(global);
        return;
    }
    max_threads = global->job_count == 0u ? 1u : global->job_count;
    if (max_threads > 16u)
        max_threads = 16u;
    duckdb_init_set_max_threads(info, max_threads);
    duckdb_init_set_init_data(info, global, bigwig_global_destroy);
}

static void
read_bigwig_local_init(duckdb_init_info info)
{
    bigwig_bind_data_t *bind;
    bigwig_local_data_t *local;
    idx_t column;

    bind = duckdb_init_get_bind_data(info);
    local = duckdb_malloc(sizeof(*local));
    if (local == NULL) {
        duckdb_init_set_error(info, "read_bigwig: out of memory");
        return;
    }
    memset(local, 0, sizeof(*local));
    local->bind = bind;
    local->column_count = duckdb_init_get_column_count(info);
    if (local->column_count != 0u) {
        local->column_ids = duckdb_malloc(sizeof(*local->column_ids) *
            (size_t)local->column_count);
        if (local->column_ids == NULL) {
            duckdb_init_set_error(info, "read_bigwig: out of memory");
            bigwig_local_destroy(local);
            return;
        }
        for (column = 0; column < local->column_count; column++)
            local->column_ids[column] = duckdb_init_get_column_index(info,
                column);
    }
    duckdb_init_set_init_data(info, local, bigwig_local_destroy);
}

static int
bigwig_local_open(bigwig_local_data_t *local, bigwig_global_data_t *global)
{
    if (local->file != NULL)
        return 1;
    if (__sync_bool_compare_and_swap(&global->seed_claimed, 0, 1)) {
        local->file = global->seed_file;
        global->seed_file = NULL;
        return local->file != NULL;
    }
    return bigwig_open_checked(local->bind->path, &local->file);
}

static int
bigwig_claim_job(bigwig_local_data_t *local, bigwig_global_data_t *global)
{
    bigwig_job_t *job;
    uint32_t job_index;

    bwIteratorDestroy(local->iterator);
    local->iterator = NULL;
    local->interval_index = 0u;
    job_index = __sync_fetch_and_add(&global->next_job, 1u);
    if (job_index >= global->job_count)
        return 0;
    job = &global->jobs[job_index];
    if (job->tid >= (uint64_t)local->file->cl->nKeys)
        return -1;
    local->owner_start0 = job->owner_start0;
    local->owner_end0 = job->owner_end0;
    local->iterator = bwOverlappingIntervalsIterator(local->file,
        local->file->cl->chrom[job->tid], job->start0, job->end0,
        local->bind->blocks_per_iteration);
    if (local->iterator == NULL)
        return -1;
    if (local->iterator->blocks == NULL)
        return -1;
    if (((bwOverlapBlock_t *)local->iterator->blocks)->n != 0u &&
        local->iterator->data == NULL)
        return -1;
    return 1;
}

static void
read_bigwig_scan(duckdb_function_info info, duckdb_data_chunk output)
{
    bigwig_global_data_t *global;
    bigwig_local_data_t *local;
    duckdb_vector vectors[BIGWIG_COL_COUNT];
    idx_t row_count, vector_size, column;

    global = duckdb_function_get_init_data(info);
    local = duckdb_function_get_local_init_data(info);
    if (global == NULL || local == NULL || local->done) {
        duckdb_data_chunk_set_size(output, 0u);
        return;
    }
    if (global->job_count == 0u) {
        local->done = 1;
        duckdb_data_chunk_set_size(output, 0u);
        return;
    }
    if (!bigwig_local_open(local, global)) {
        duckdb_function_set_error(info,
            "read_bigwig: failed to open a worker-local BigWig handle");
        local->done = 1;
        duckdb_data_chunk_set_size(output, 0u);
        return;
    }
    if (local->column_count != 0u) {
        for (column = 0; column < local->column_count; column++)
            vectors[column] = duckdb_data_chunk_get_vector(output, column);
    }

    row_count = 0u;
    vector_size = duckdb_vector_size();
    while (row_count < vector_size && !local->done) {
        bwOverlappingIntervals_t *intervals;
        const char *chrom;

        if (local->iterator == NULL || local->iterator->data == NULL) {
            int claim_status = bigwig_claim_job(local, global);

            if (claim_status < 0) {
                duckdb_function_set_error(info,
                    "read_bigwig: failed to initialize an indexed range");
                local->done = 1;
                break;
            }
            if (claim_status == 0) {
                local->done = 1;
                break;
            }
            if (local->iterator->data == NULL)
                continue;
        }
        intervals = local->iterator->intervals;
        if (local->interval_index >= intervals->l) {
            bwOverlapIterator_t *next = bwIteratorNext(local->iterator);

            if (next == NULL) {
                local->iterator = NULL;
                duckdb_function_set_error(info,
                    "read_bigwig: failed while reading an indexed range");
                local->done = 1;
                break;
            }
            local->iterator = next;
            local->interval_index = 0u;
            continue;
        }
        if (intervals->start[local->interval_index] < local->owner_start0 ||
            intervals->start[local->interval_index] >= local->owner_end0) {
            local->interval_index++;
            continue;
        }
        chrom = local->file->cl->chrom[local->iterator->tid];
        for (column = 0; column < local->column_count; column++) {
            idx_t logical = local->column_ids[column];

            switch (logical) {
            case BIGWIG_COL_CHROM:
                duckdb_vector_assign_string_element(vectors[column],
                    row_count, chrom);
                break;
            case BIGWIG_COL_START:
                ((uint32_t *)duckdb_vector_get_data(vectors[column]))[row_count] =
                    intervals->start[local->interval_index];
                break;
            case BIGWIG_COL_END:
                ((uint32_t *)duckdb_vector_get_data(vectors[column]))[row_count] =
                    intervals->end[local->interval_index];
                break;
            case BIGWIG_COL_VALUE:
                ((float *)duckdb_vector_get_data(vectors[column]))[row_count] =
                    intervals->value[local->interval_index];
                break;
            default:
                duckdb_function_set_error(info,
                    "read_bigwig: invalid projected column");
                duckdb_data_chunk_set_size(output, 0u);
                return;
            }
        }
        local->interval_index++;
        row_count++;
    }
    duckdb_data_chunk_set_size(output, row_count);
}

void
register_read_bigwig_function(duckdb_connection connection)
{
    duckdb_table_function function;
    duckdb_logical_type varchar_type, bigint_type;

    function = duckdb_create_table_function();
    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_table_function_set_name(function, "read_bigwig");
    duckdb_table_function_add_parameter(function, varchar_type);
    duckdb_table_function_add_named_parameter(function, "region", varchar_type);
    duckdb_table_function_add_named_parameter(function,
        "blocks_per_iteration", bigint_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_table_function_set_bind(function, read_bigwig_bind);
    duckdb_table_function_set_init(function, read_bigwig_global_init);
    duckdb_table_function_set_local_init(function, read_bigwig_local_init);
    duckdb_table_function_set_function(function, read_bigwig_scan);
    duckdb_table_function_supports_projection_pushdown(function, true);
    duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
}
