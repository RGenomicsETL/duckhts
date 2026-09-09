/* Query-scoped DuckDB adapter for the native phased replay stream. DuckDB owns
 * input materialization, phase-domain aggregation, sorting and output vectors.
 * The C workspace is allocated at init and retains only the active window. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "duckdb_list.h"

#include "duckvep_model.h"
#include "kernel/src/duckvep_haplotype_stream.h"
#include "kernel/src/duckvep_sequence_diff.h"
#include "kernel/src/duckvep_dna.h"
#include "kernel/src/duckvep_effect.h"
#include "kernel/src/duckvep_hgvs.h"
#include "kernel/src/duckvep_annotation_internal.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { LIMIT_EVENTS, LIMIT_TRANSCRIPTS, LIMIT_CARRIERS, LIMIT_PREFIXES, LIMIT_PROJECTIONS,
    LIMIT_ALLELES, LIMIT_LEAF_EVENTS, LIMIT_LEAF_EDITS, LIMIT_SEQUENCE, LIMIT_PLOIDY,
    LIMIT_PHASE_SETS, LIMIT_ALIGNMENT, LIMIT_DIFFERENCES, LIMIT_HGVS_OPERATIONS,
    LIMIT_HGVS_BYTES, LIMIT_HGVS_REFERENCE, LIMIT_WORKSPACE, LIMIT_COUNT };
enum { HAPLOTYPE_LIST_COLUMN = 9, HAPLOTYPE_STOP_COLUMN = 14,
    HAPLOTYPE_HGVSP_COLUMN = 15, HAPLOTYPE_HGVSP_STATUS_COLUMN = 16, HAPLOTYPE_OUTPUT_COLUMNS = 17 };
enum { HAPLOTYPE_BLOCK_EVENT_FIELD = 9, HAPLOTYPE_BLOCK_FIELDS = 10 };
static const char *const limit_names[] = {"max_active_events", "max_active_transcripts",
    "max_active_carriers", "max_active_prefixes", "max_active_projections", "max_allele_bytes",
    "max_leaf_events", "max_leaf_edits", "max_sequence_bases", "max_ploidy", "max_phase_sets",
    "max_alignment_cells", "max_leaf_differences", "max_hgvs_operations", "max_hgvs_bytes",
    "max_hgvs_reference_bytes", "workspace_limit"};
static const uint64_t limit_defaults[] = {16384, 4096, 65536, 262144, 262144, 8388608,
    4096, 65536, 1048576, 64, 1024, 16777216, 65536, 65536, 1048576,
    DUCKVEP_REFERENCE_DEFAULT_BYTES, 268435456};

typedef struct {
    duckvep_registry_t *registry;
    duckvep_model_entry_t *entry;
    char *query;
    duckvep_phase_policy_t policy;
    int source_records, hgvs;
    size_t limits[LIMIT_COUNT];
} haplotype_bind_t;

typedef struct {
    duckdb_result input;
    duckdb_data_chunk chunk;
    idx_t row;
    int have_result, eof, have_call;
    uint32_t last_tx;
    duckvep_haplotype_stream_t stream;
    duckvep_haplotype_stream_buffers_t buffers;
    int32_t *gt;
    uint8_t *phase;
    duckvep_haplotype_phase_set_t *sets;
    duckvep_sequence_diff_scratch_t difference_scratch;
    duckvep_sequence_difference_t *differences;
    uint8_t *difference_reference;
    uint8_t *reference_coding_protein;
    duckvep_translation_t reference_coding_translation;
    uint32_t difference_transcript;
    int have_difference_reference;
    duckvep_hgvs_protein_operation_t *protein_operations;
    duckvep_reference_reader_t reference;
    duckvep_delta_scratch_t hgvs_scratch;
    uint8_t *hgvs_shifted_allele;
    size_t hgvs_shifted_capacity;
    char *hgvsp;
    size_t workspace_bytes;
} haplotype_state_t;

static void haplotype_bind_destroy(void *pointer) {
    haplotype_bind_t *b = pointer;
    if (!b) return;
    duckvep_registry_unpin(b->registry, b->entry);
    duckvep_registry_release(b->registry);
    duckdb_free(b->query);
    free(b);
}

static duckdb_logical_type record_type(const char *const *names, const duckdb_type *ids,
    idx_t count, int last_field_is_list) {
    duckdb_logical_type types[HAPLOTYPE_BLOCK_FIELDS];
    const char *field_names[HAPLOTYPE_BLOCK_FIELDS];
    for (idx_t i = 0u; i < count; i++) {
        types[i] = duckdb_create_logical_type(ids[i]); field_names[i] = names[i];
        if (last_field_is_list && i == count - 1u) {
            duckdb_logical_type element = types[i];
            types[i] = duckdb_create_list_type(element);
            duckdb_destroy_logical_type(&element);
        }
    }
    duckdb_logical_type type = duckdb_create_struct_type(types, field_names, count);
    for (idx_t i = 0u; i < count; i++) duckdb_destroy_logical_type(&types[i]);
    return type;
}

static void bind_record_list(duckdb_bind_info info, const char *name,
    const char *const *fields, const duckdb_type *ids, idx_t count, int last_field_is_list) {
    duckdb_logical_type record = record_type(fields, ids, count, last_field_is_list);
    duckdb_logical_type list = duckdb_create_list_type(record);
    duckdb_bind_add_result_column(info, name, list);
    duckdb_destroy_logical_type(&list);
    duckdb_destroy_logical_type(&record);
}

static void haplotype_bind(duckdb_bind_info info) {
    haplotype_bind_t *b = calloc(1u, sizeof(*b));
    if (!b) { duckdb_bind_set_error(info, "duckvep_haplotypes: bind allocation failed"); return; }
    b->registry = duckdb_bind_get_extra_info(info);
    duckvep_registry_retain(b->registry);
    duckdb_value value = duckdb_bind_get_parameter(info, 0u);
    if (value && !duckdb_is_null_value(value)) b->query = duckdb_get_varchar(value);
    duckdb_destroy_value(&value);
    value = duckdb_bind_get_parameter(info, 1u);
    char *name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    duckdb_destroy_value(&value);
    if (name) b->entry = duckvep_registry_pin(b->registry, name);
    duckdb_free(name);
    if (!b->entry || !b->query || !b->query[0]) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: require a nonempty calls query and loaded model name");
        haplotype_bind_destroy(b); return;
    }
    value = duckdb_bind_get_named_parameter(info, "phase_policy");
    name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    int valid = !value || (name && (!strcmp(name, "strict") || !strcmp(name, "vep116_compat")));
    b->policy = name && !strcmp(name, "vep116_compat") ? DUCKVEP_PHASE_VEP116_COMPAT : DUCKVEP_PHASE_STRICT;
    duckdb_free(name); duckdb_destroy_value(&value);
    if (!valid) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: phase_policy must be 'strict' or 'vep116_compat'");
        haplotype_bind_destroy(b); return;
    }
    value = duckdb_bind_get_named_parameter(info, "input_mode");
    name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    valid = !value || (name && (!strcmp(name, "alt_events") || !strcmp(name, "source_records")));
    b->source_records = name && !strcmp(name, "source_records");
    duckdb_free(name); duckdb_destroy_value(&value);
    if (!valid || (b->source_records && b->policy != DUCKVEP_PHASE_VEP116_COMPAT)) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: input_mode must be 'alt_events' or 'source_records'; source_records requires phase_policy='vep116_compat'");
        haplotype_bind_destroy(b); return;
    }
    value = duckdb_bind_get_named_parameter(info, "hgvs");
    b->hgvs = value && !duckdb_is_null_value(value) && duckdb_get_bool(value);
    duckdb_destroy_value(&value);
    for (unsigned i = 0u; i < LIMIT_COUNT; i++) {
        value = duckdb_bind_get_named_parameter(info, limit_names[i]);
        uint64_t n = value && !duckdb_is_null_value(value) ? duckdb_get_uint64(value) : limit_defaults[i];
        valid = (!value || !duckdb_is_null_value(value)) && n && n <= SIZE_MAX;
        duckdb_destroy_value(&value);
        if (i <= LIMIT_PROJECTIONS && n > (UINT32_C(1) << 29)) valid = 0;
        if (i == LIMIT_PLOIDY && n > UINT16_MAX) valid = 0;
        if (!valid) {
            char error[128]; snprintf(error, sizeof(error), "duckvep_haplotypes: invalid %s", limit_names[i]);
            duckdb_bind_set_error(info, error); haplotype_bind_destroy(b); return;
        }
        b->limits[i] = (size_t)n;
    }
    const char *const names[] = {"transcript_index", "cds", "protein", "sequence_flags",
        "evidence_flags", "projection_status", "sequence_status", "edit_count", "carrier_count"};
    const duckdb_type ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UINTEGER};
    for (unsigned i = 0u; i < 9u; i++) {
        duckdb_logical_type type = duckdb_create_logical_type(ids[i]);
        duckdb_bind_add_result_column(info, names[i], type); duckdb_destroy_logical_type(&type);
    }
    const char *const carrier_names[] = {"sample_index", "phase_set", "haplotype_lane", "ploidy"};
    const duckdb_type carrier_ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_BIGINT,
        DUCKDB_TYPE_USMALLINT, DUCKDB_TYPE_USMALLINT};
    const char *const event_names[] = {"event_index", "seq_region", "position", "reference", "alternate",
        "evidence_flags", "projection_status", "alt_index"};
    const duckdb_type event_ids[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT,
        DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UINTEGER};
    const char *const block_names[] = {"cds_start", "reference", "alternate", "alt_start0",
        "length_change", "sequence_flags", "coding_status", "local_consequence_mask",
        "after_first_stop", "event_indices"};
    const duckdb_type block_ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_BIGINT, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_UBIGINT};
    bind_record_list(info, "carriers", carrier_names, carrier_ids, 4u, 0);
    bind_record_list(info, "contributors", event_names, event_ids, b->source_records ? 8u : 7u, 0);
    bind_record_list(info, "coding_blocks", block_names, block_ids, HAPLOTYPE_BLOCK_FIELDS, 1);
    const char *const difference_names[] = {"ref_start0", "alt_start0", "reference", "alternate",
        "alignment_start0"};
    const duckdb_type difference_ids[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
        DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UBIGINT};
    bind_record_list(info, "cds_differences", difference_names, difference_ids, 5u, 0);
    bind_record_list(info, "protein_differences", difference_names, difference_ids, 5u, 0);
    duckdb_logical_type stop_in_frame_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_bind_add_result_column(info, "stop_in_displaced_frame", stop_in_frame_type);
    duckdb_destroy_logical_type(&stop_in_frame_type);
    duckdb_logical_type string_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "hgvsp", string_type);
    duckdb_bind_add_result_column(info, "hgvsp_status", string_type);
    duckdb_destroy_logical_type(&string_type);
    duckdb_bind_set_bind_data(info, b, haplotype_bind_destroy);
}

static void haplotype_state_destroy(void *pointer) {
    haplotype_state_t *s = pointer;
    if (!s) return;
    if (s->chunk) duckdb_destroy_data_chunk(&s->chunk);
    if (s->have_result) duckdb_destroy_result(&s->input);
    duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    free(b->carriers.transcripts); free(b->carriers.calls); free(b->carriers.prefixes);
    free(b->carriers.active_transcripts); free(b->carriers.transcript_index);
    free(b->carriers.call_index); free(b->carriers.prefix_index);
    free(b->events); free(b->projections); free(b->alleles); free(b->leaf_events);
    free(b->contributors); free(b->edits); free(b->blocks); free(b->cds); free(b->protein);
    free(b->edit_event_ids);
    free(s->difference_scratch.scores); free(s->difference_scratch.trace); free(s->differences);
    free(s->difference_reference); free(b->reference_protein);
    free(s->reference_coding_protein);
    free(s->protein_operations); free(s->hgvsp);
    if (s->reference.fai) fai_destroy(s->reference.fai);
    free(s->reference.bases);
    free(s->hgvs_scratch.edits); free(s->hgvs_scratch.alt_cds);
    free(s->hgvs_scratch.ref_peptide); free(s->hgvs_scratch.alt_peptide);
    free(s->hgvs_shifted_allele);
    free(s->gt); free(s->phase); free(s->sets); free(s);
}

static uint32_t bucket_count(uint32_t capacity) {
    uint32_t n = 2u;
    while (n < capacity * 2u) n *= 2u;
    return n;
}

static int workspace_allocate(haplotype_state_t *s, const haplotype_bind_t *bind) {
    const size_t *n = bind->limits;
    duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    b->carriers.transcript_capacity = (uint32_t)n[LIMIT_TRANSCRIPTS];
    b->carriers.call_capacity = (uint32_t)n[LIMIT_CARRIERS];
    b->carriers.prefix_capacity = (uint32_t)n[LIMIT_PREFIXES];
    b->carriers.transcript_buckets = bucket_count(b->carriers.transcript_capacity);
    b->carriers.call_buckets = bucket_count(b->carriers.call_capacity);
    b->carriers.prefix_buckets = bucket_count(b->carriers.prefix_capacity);
    b->event_capacity = (uint32_t)n[LIMIT_EVENTS]; b->projection_capacity = (uint32_t)n[LIMIT_PROJECTIONS];
    b->allele_capacity = n[LIMIT_ALLELES]; b->leaf_capacity = n[LIMIT_LEAF_EVENTS];
    b->edit_capacity = n[LIMIT_LEAF_EDITS];
    if (n[LIMIT_SEQUENCE] == SIZE_MAX || (bind->hgvs && n[LIMIT_HGVS_BYTES] == SIZE_MAX)) return 0;
    b->cds_capacity = n[LIMIT_SEQUENCE];
    b->protein_capacity = n[LIMIT_SEQUENCE] + 1u;
    if (b->protein_capacity > SIZE_MAX / 2u) return 0;
    s->difference_scratch.score_capacity = b->protein_capacity * 2u;
    s->difference_scratch.trace_capacity = n[LIMIT_ALIGNMENT];
    b->reference_protein_capacity = n[LIMIT_SEQUENCE] / 3u + 2u;
    if (bind->hgvs) {
        s->reference.capacity = bind->entry->model.reference_fasta_path ? n[LIMIT_HGVS_REFERENCE] : 0u;
        /* Alternating differing/retained bases maximize one uint16-length MNV's islands. */
        size_t islands = ((size_t)UINT16_MAX + 1u) / 2u;
        s->hgvs_scratch.edits_cap = n[LIMIT_LEAF_EDITS] < islands ? n[LIMIT_LEAF_EDITS] : islands;
        s->hgvs_scratch.alt_cds_cap = n[LIMIT_SEQUENCE];
        s->hgvs_scratch.ref_peptide_cap = s->hgvs_scratch.alt_peptide_cap = n[LIMIT_SEQUENCE] / 3u + 1u;
        s->hgvs_shifted_capacity = n[LIMIT_ALLELES] < UINT16_MAX + 1u ? n[LIMIT_ALLELES] : UINT16_MAX + 1u;
    }
    /* Count every byte before allocating any of the arrays. Each pointer has
     * one owner and one cleanup site; no allocator is used by the scan loop. */
#define ARRAYS(X) \
    X(b->carriers.transcripts, b->carriers.transcript_capacity) \
    X(b->carriers.calls, b->carriers.call_capacity) \
    X(b->carriers.prefixes, b->carriers.prefix_capacity) \
    X(b->carriers.active_transcripts, b->carriers.transcript_capacity) \
    X(b->carriers.transcript_index, b->carriers.transcript_buckets) \
    X(b->carriers.call_index, b->carriers.call_buckets) \
    X(b->carriers.prefix_index, b->carriers.prefix_buckets) \
    X(b->events, b->event_capacity) X(b->projections, b->projection_capacity) \
    X(b->alleles, b->allele_capacity) X(b->leaf_events, b->leaf_capacity) \
    X(b->contributors, b->leaf_capacity) X(b->edits, b->edit_capacity) \
    X(b->blocks, b->edit_capacity) X(b->edit_event_ids, b->edit_capacity) \
    X(b->cds, b->cds_capacity) X(b->protein, b->protein_capacity) \
    X(s->difference_scratch.scores, s->difference_scratch.score_capacity) \
    X(s->difference_scratch.trace, s->difference_scratch.trace_capacity) \
    X(s->differences, n[LIMIT_DIFFERENCES]) \
    X(s->difference_reference, n[LIMIT_SEQUENCE]) \
    X(b->reference_protein, b->reference_protein_capacity) \
    X(s->reference_coding_protein, b->reference_protein_capacity) \
    X(s->protein_operations, bind->hgvs ? n[LIMIT_HGVS_OPERATIONS] : 0u) \
    X(s->hgvsp, bind->hgvs ? n[LIMIT_HGVS_BYTES] + 1u : 0u) \
    X(s->reference.bases, s->reference.capacity) \
    X(s->hgvs_scratch.edits, s->hgvs_scratch.edits_cap) \
    X(s->hgvs_scratch.alt_cds, s->hgvs_scratch.alt_cds_cap) \
    X(s->hgvs_scratch.ref_peptide, s->hgvs_scratch.ref_peptide_cap) \
    X(s->hgvs_scratch.alt_peptide, s->hgvs_scratch.alt_peptide_cap) \
    X(s->hgvs_shifted_allele, s->hgvs_shifted_capacity) \
    X(s->gt, bind->source_records ? 0u : n[LIMIT_PLOIDY]) \
    X(s->phase, bind->source_records ? 0u : n[LIMIT_PLOIDY]) \
    X(s->sets, bind->source_records ? 0u : n[LIMIT_PHASE_SETS])
#define COUNT(p, count) \
    if ((count) > (n[LIMIT_WORKSPACE] - s->workspace_bytes) / sizeof(*(p))) return 0; \
    s->workspace_bytes += (count) * sizeof(*(p));
    s->workspace_bytes = sizeof(*s);
    if (s->workspace_bytes > n[LIMIT_WORKSPACE]) return 0;
    ARRAYS(COUNT)
#undef COUNT
#define ALLOCATE(p, count) if ((count) && !((p) = malloc((count) * sizeof(*(p))))) return 0;
    ARRAYS(ALLOCATE)
#undef ALLOCATE
#undef ARRAYS
    return 1;
}

static int raw_prepare_query(duckdb_connection connection, const char *sql,
    duckdb_result *result, char *error, size_t error_size) {
    if (duckdb_query(connection, sql, result) == DuckDBSuccess) return 1;
    duckvep_sql_set_error(error, error_size, duckdb_result_error(result));
    return 0;
}

/* Called under the registry preparation lock. The SELECT is evaluated once;
 * DuckDB owns both temporary relations and their spill. No native record array
 * grows with the source stream. Drop only tables created by this invocation. */
static int raw_prepare(haplotype_state_t *s, const haplotype_bind_t *b,
    char *error, size_t error_size) {
    duckdb_connection connection = b->registry->query_connection;
    duckdb_result result = {0};
    duckdb_appender appender = NULL;
    duckdb_data_chunk chunk = NULL;
    int have_raw = 0, have_order = 0, ok = 0;
    if (!raw_prepare_query(connection,
        "CREATE TEMP TABLE __duckvep_haplotype_raw(event_index UBIGINT, seq_region UINTEGER, "
        "position UBIGINT, reference VARCHAR, alternates VARCHAR[], transcript_index UINTEGER, "
        "sample_index UINTEGER, gt VARCHAR)", &result, error, error_size)) goto cleanup;
    have_raw = 1; duckdb_destroy_result(&result);
    if (duckdb_appender_create(connection, NULL, "__duckvep_haplotype_raw", &appender) != DuckDBSuccess)
        goto append_failed;
    while ((chunk = duckdb_fetch_chunk(s->input))) {
        if (duckdb_append_data_chunk(appender, chunk) != DuckDBSuccess) goto append_failed;
        duckdb_destroy_data_chunk(&chunk);
    }
    if (duckdb_result_error(&s->input)) {
        duckvep_sql_set_error(error, error_size, duckdb_result_error(&s->input)); goto cleanup;
    }
    if (duckdb_appender_close(appender) != DuckDBSuccess) goto append_failed;
    duckdb_appender_destroy(&appender);
    duckdb_destroy_result(&s->input); s->have_result = 0;
    if (!raw_prepare_query(connection,
        "CREATE TEMP TABLE __duckvep_haplotype_order(event_index UBIGINT, buffer_id UBIGINT, "
        "ordinal UBIGINT, source_ordinal UBIGINT)", &result, error, error_size)) goto cleanup;
    have_order = 1; duckdb_destroy_result(&result);
    if (!raw_prepare_query(connection,
        "SELECT event_index, first(seq_region), first(position), first(reference), "
        "count(DISTINCT (seq_region,position,reference,alternates)) versions "
        "FROM __duckvep_haplotype_raw GROUP BY event_index ORDER BY 2,3,1",
        &result, error, error_size)) goto cleanup;
    if (duckdb_appender_create(connection, NULL, "__duckvep_haplotype_order", &appender) != DuckDBSuccess)
        goto append_failed;
    duckvep_haplotype_record_plan_t plan;
    if (!duckvep_haplotype_record_plan_init(&plan, &b->entry->model.transcripts)) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid source-planning model");
        goto cleanup;
    }
    uint64_t source_ordinal = 0u;
    while ((chunk = duckdb_fetch_chunk(result))) {
        duckdb_vector v[5];
        for (unsigned i = 0u; i < 5u; i++) v[i] = duckdb_data_chunk_get_vector(chunk, i);
        for (idx_t row = 0u; row < duckdb_data_chunk_get_size(chunk); row++) {
            for (unsigned i = 0u; i < 4u; i++) if (duckvep_row_is_null(v[i], row)) {
                duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: required input source column is NULL");
                goto cleanup;
            }
            if (((int64_t *)duckdb_vector_get_data(v[4]))[row] != 1) {
                duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: inconsistent source record identity");
                goto cleanup;
            }
            uint32_t chrom = ((uint32_t *)duckdb_vector_get_data(v[1]))[row];
            uint64_t pos = ((uint64_t *)duckdb_vector_get_data(v[2]))[row];
            duckdb_string_t ref = ((duckdb_string_t *)duckdb_vector_get_data(v[3]))[row];
            uint32_t length = duckdb_string_t_length(ref);
            uint64_t buffer, ordinal;
            if (chrom > UINT16_MAX || !pos || pos > UINT32_MAX || !length || length > UINT16_MAX ||
                length - 1u > UINT32_MAX - pos || source_ordinal == UINT64_MAX ||
                !duckvep_haplotype_record_plan_next(&plan, (uint16_t)chrom, (uint32_t)pos,
                    (uint32_t)(pos + length - 1u), &buffer, &ordinal)) {
                duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid source span or record count");
                goto cleanup;
            }
            source_ordinal++;
            if (duckdb_append_uint64(appender, ((uint64_t *)duckdb_vector_get_data(v[0]))[row]) != DuckDBSuccess ||
                duckdb_append_uint64(appender, buffer) != DuckDBSuccess ||
                duckdb_append_uint64(appender, ordinal) != DuckDBSuccess ||
                duckdb_append_uint64(appender, source_ordinal) != DuckDBSuccess ||
                duckdb_appender_end_row(appender) != DuckDBSuccess) goto append_failed;
        }
        duckdb_destroy_data_chunk(&chunk);
    }
    if (duckdb_result_error(&result)) {
        duckvep_sql_set_error(error, error_size, duckdb_result_error(&result)); goto cleanup;
    }
    if (duckdb_appender_close(appender) != DuckDBSuccess) goto append_failed;
    duckdb_appender_destroy(&appender); duckdb_destroy_result(&result);
    s->have_result = 1;
    ok = raw_prepare_query(connection,
        "WITH ordered AS (SELECT event_index, CASE WHEN buffer_id=0 THEN source_ordinal ELSE "
        "source_ordinal-ordinal+_duckvep_record_order(count(*) OVER(PARTITION BY buffer_id)::UBIGINT,ordinal) "
        "END replay_order FROM __duckvep_haplotype_order), "
        "genotypes AS (SELECT event_index,sample_index,count(DISTINCT gt) gt_versions "
        "FROM __duckvep_haplotype_raw GROUP BY event_index,sample_index), "
        "calls AS MATERIALIZED (SELECT *, count(*) OVER(PARTITION BY event_index,transcript_index,sample_index) copies, "
        "CASE WHEN alternates IS NULL OR len(alternates)>2147483647 OR "
        "len(list_filter(alternates,a -> a IS NULL OR len(a)=0 OR len(a)>65535))>0 "
        "THEN error('duckvep_haplotypes: invalid source ALT list') ELSE len(alternates) END alt_count "
        "FROM __duckvep_haplotype_raw JOIN ordered USING(event_index)), "
        "parsed AS MATERIALIZED (SELECT *, _duckvep_raw_gt(gt,alt_count::UINTEGER) raw_gt FROM calls), "
        "selected AS (SELECT *, max(replay_order) FILTER(WHERE raw_gt.disposition=3) OVER "
        "(PARTITION BY seq_region,position,reference,alternates,transcript_index) selected_order FROM parsed) "
        "SELECT c.event_index,seq_region,position,reference, "
        "CASE WHEN a.i=0 THEN reference WHEN a.i>alt_count THEN '' ELSE alternates[a.i] END alternate, "
        "(CASE WHEN a.i>alt_count THEN 4294967295 ELSE a.i END)::UINTEGER alt_index, "
        "transcript_index,c.sample_index,raw_gt,NULL::BOOLEAN[] phase_before,NULL::BIGINT phase_set, "
        "alt_count,copies,1::BIGINT versions,gt_versions,replay_order, "
        "(selected_order IS NULL OR replay_order=selected_order) source_selected "
        "FROM selected c LEFT JOIN genotypes g USING(event_index,sample_index),range(0,alt_count+2) a(i) "
        "ORDER BY seq_region,position,event_index,alt_index,transcript_index,sample_index",
        &s->input, error, error_size);
    goto cleanup;
append_failed:
    duckvep_sql_set_error(error, error_size, appender ? duckdb_appender_error(appender)
        : "duckvep_haplotypes: could not stage source records");
cleanup:
    if (chunk) duckdb_destroy_data_chunk(&chunk);
    if (appender) duckdb_appender_destroy(&appender);
    duckdb_destroy_result(&result);
    const char *drops[] = {"DROP TABLE temp.main.__duckvep_haplotype_order", "DROP TABLE temp.main.__duckvep_haplotype_raw"};
    const int created[] = {have_order, have_raw};
    for (unsigned i = 0u; i < 2u; i++) if (created[i]) {
        if (duckdb_query(connection, drops[i], &result) != DuckDBSuccess && ok) {
            duckvep_sql_set_error(error, error_size, duckdb_result_error(&result)); ok = 0;
        }
        duckdb_destroy_result(&result);
    }
    return ok;
}

static int input_open(haplotype_state_t *s, const haplotype_bind_t *b, char *error, size_t error_size) {
    /* One SELECT statement and snapshot on the registry's retained connection.
     * Only query preparation/materialization is serialized. The returned result
     * is owned by this scan and fetched without the registry mutex. TEMP objects
     * and uncommitted writes on the caller's connection are not visible here. */
    const char *prefix =
        "WITH raw AS MATERIALIZED (SELECT event_index::UBIGINT event_index, seq_region::UINTEGER seq_region, "
        "position::UBIGINT AS position, reference::VARCHAR AS reference, alternate::VARCHAR AS alternate, "
        "alt_index::UINTEGER alt_index, transcript_index::UINTEGER transcript_index, sample_index::UINTEGER sample_index, "
        "alleles::INTEGER[] alleles, phase_before::BOOLEAN[] phase_before, phase_set::BIGINT phase_set FROM (";
    const char *middle =
        ") source), calls AS MATERIALIZED (SELECT *, "
        "list_contains(list_transform(duckvep_phase_call(alleles,phase_before,phase_set := phase_set), "
        "a -> a.phase_scope), 'phase_set') scoped FROM raw), domains AS (SELECT transcript_index, sample_index, ";
    const char *domain = b->policy == DUCKVEP_PHASE_STRICT ?
        "coalesce(list(DISTINCT phase_set ORDER BY phase_set NULLS FIRST) FILTER(WHERE scoped), [NULL]::BIGINT[]) AS domain_sets " :
        "[NULL]::BIGINT[] AS domain_sets ";
    /* Validate each identity/domain once, then attach its facts to every call.
     * Duplicate calls remain visible, including calls carrying only REF. */
    const char *suffix =
        ",count(DISTINCT len(alleles)) ploidies FROM calls GROUP BY transcript_index,sample_index), "
        "event_versions AS (SELECT event_index, "
        "count(DISTINCT (seq_region,position,reference,alternate,alt_index)) versions FROM raw GROUP BY event_index) "
        "SELECT c.* EXCLUDE(scoped), d.domain_sets, "
        "count(*) OVER(PARTITION BY c.event_index,transcript_index,sample_index) copies, v.versions, d.ploidies "
        "FROM calls c LEFT JOIN domains d USING(transcript_index,sample_index) "
        "LEFT JOIN event_versions v USING(event_index) "
        "ORDER BY seq_region,position,event_index,transcript_index,sample_index";
    if (b->source_records) {
        prefix = "SELECT event_index::UBIGINT event_index, "
            "seq_region::UINTEGER seq_region, position::UBIGINT AS position, reference::VARCHAR AS reference, "
            "alternates::VARCHAR[] alternates, transcript_index::UINTEGER transcript_index, "
            "sample_index::UINTEGER sample_index, gt::VARCHAR gt FROM (";
        middle = ") source";
        domain = ""; suffix = "";
    }
    size_t qlen = strlen(b->query), overhead = strlen(prefix) + strlen(middle) + strlen(domain) + strlen(suffix) + 1u;
    if (qlen > SIZE_MAX - overhead || qlen + overhead > b->limits[LIMIT_WORKSPACE] - s->workspace_bytes) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: workspace_limit exceeded by calls query text");
        return 0;
    }
    char *sql = malloc(qlen + overhead);
    if (!sql) return 0;
    snprintf(sql, qlen + overhead, "%s%s%s%s%s", prefix, b->query, middle, domain, suffix);
    duckdb_prepared_statement statement = NULL;
    if (!duckvep_registry_query_acquire(b->registry, error, error_size)) { free(sql); return 0; }
    duckdb_extracted_statements extracted = NULL;
    idx_t statements = duckdb_extract_statements(b->registry->query_connection, sql, &extracted);
    int ok = statements == 1u && duckdb_prepare_extracted_statement(
        b->registry->query_connection, extracted, 0u, &statement) == DuckDBSuccess;
    free(sql);
    if (!ok) {
        const char *message = statement ? duckdb_prepare_error(statement) :
            extracted ? duckdb_extract_statements_error(extracted) : NULL;
        duckvep_sql_set_error(error, error_size, message ? message : "calls query must be one SELECT statement");
    } else if (duckdb_prepared_statement_type(statement) != DUCKDB_STATEMENT_TYPE_SELECT) {
        duckvep_sql_set_error(error, error_size, "calls query must be one SELECT statement"); ok = 0;
    } else {
        s->have_result = 1;
        ok = duckdb_execute_prepared(statement, &s->input) == DuckDBSuccess;
        if (!ok) duckvep_sql_set_error(error, error_size, duckdb_result_error(&s->input));
    }
    if (statement) duckdb_destroy_prepare(&statement);
    if (extracted) duckdb_destroy_extracted(&extracted);
    if (ok && b->source_records) ok = raw_prepare(s, b, error, error_size);
    pthread_mutex_unlock(&b->registry->query_mutex);
    return ok;
}

static void haplotype_init(duckdb_init_info info) {
    const haplotype_bind_t *bind = duckdb_init_get_bind_data(info);
    haplotype_state_t *s = calloc(1u, sizeof(*s));
    char error[DUCKVEP_SQL_ERROR_SIZE] = "duckvep_haplotypes: workspace allocation or configured limit exceeded";
    duckdb_init_set_max_threads(info, 1u);
    if (!s || !workspace_allocate(s, bind)) goto failed;
    duckvep_owned_model_t *m = &bind->entry->model;
    if (bind->hgvs && !duckvep_reference_reader_init(&s->reference, m,
            s->reference.bases, s->reference.capacity, error, sizeof(error))) goto failed;
    if (duckvep_haplotype_stream_init(&s->stream, &m->transcripts, &m->exons, &m->sequences,
        &s->buffers) != DUCKVEP_HAPLOTYPE_STREAM_OK) {
        duckvep_sql_set_error(error, sizeof(error), "duckvep_haplotypes: invalid native model/workspace"); goto failed;
    }
    if (!input_open(s, bind, error, sizeof(error))) goto failed;
    duckdb_init_set_init_data(info, s, haplotype_state_destroy);
    return;
failed:
    duckdb_init_set_error(info, error); haplotype_state_destroy(s);
}

static const char *projection_name(duckvep_cds_edit_status_t status) {
    switch (status) {
    case DUCKVEP_CDS_EDIT_OK: return "ok";
    case DUCKVEP_CDS_EDIT_REF_MISMATCH: return "reference_mismatch";
    case DUCKVEP_CDS_EDIT_SOURCE_SHADOWED: return "shadowed_duplicate";
    case DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED: return "source_unmapped";
    case DUCKVEP_CDS_EDIT_INVALID_ARG: return "invalid_argument";
    case DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND: return "unsupported_kind";
    case DUCKVEP_CDS_EDIT_INVALID_EVENT: return "invalid_event";
    case DUCKVEP_CDS_EDIT_OUT_OF_CDS: return "outside_cds";
    case DUCKVEP_CDS_EDIT_NON_CONTIGUOUS: return "non_contiguous";
    case DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL: return "projection_capacity";
    case DUCKVEP_CDS_EDIT_INVALID_ALLELE: return "invalid_allele";
    default: return "unavailable_projection";
    }
}

static const char *sequence_name(duckvep_haplotype_status_t status) {
    switch (status) {
    case DUCKVEP_HAPLOTYPE_OK: return "ok";
    case DUCKVEP_HAPLOTYPE_CONDITIONAL: return "conditional";
    case DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE: return "incomplete_input";
    case DUCKVEP_HAPLOTYPE_EDIT_ORDER: return "edit_conflict";
    case DUCKVEP_HAPLOTYPE_REF_MISMATCH: return "reference_mismatch";
    case DUCKVEP_HAPLOTYPE_INVALID_BASE: return "invalid_base";
    default: return "invalid_sequence";
    }
}

static int exhausted_limit(duckvep_haplotype_stream_status_t status, duckvep_carriers_status_t carrier) {
    switch (status) {
    case DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL: return LIMIT_EVENTS;
    case DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL: return LIMIT_PROJECTIONS;
    case DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL: return LIMIT_ALLELES;
    case DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL: return LIMIT_LEAF_EVENTS;
    case DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL: return LIMIT_LEAF_EDITS;
    case DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL: return LIMIT_SEQUENCE;
    case DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR:
        if (carrier == DUCKVEP_CARRIERS_TRANSCRIPT_FULL) return LIMIT_TRANSCRIPTS;
        if (carrier == DUCKVEP_CARRIERS_CALL_FULL) return LIMIT_CARRIERS;
        if (carrier == DUCKVEP_CARRIERS_PREFIX_FULL) return LIMIT_PREFIXES;
        return -1;
    default: return -1;
    }
}

static void null_cell(duckdb_vector vector, idx_t row) {
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
}

static int prepare_difference_reference(haplotype_state_t *s, const haplotype_bind_t *bind,
    const duckvep_haplotype_leaf_t *leaf, char *error, size_t error_size) {
    uint32_t tx = leaf->carriers.transcript_index;
    if (!leaf->cds || (s->have_difference_reference && tx == s->difference_transcript)) return 1;
    const duckvep_sequence_pool_t *seq = s->stream.sequences;
    size_t length = seq->cds_length[tx];
    if (length > bind->limits[LIMIT_SEQUENCE]) {
        snprintf(error, error_size,
            "duckvep_haplotypes: max_sequence_bases=%zu, reference requires=%zu at transcript %u",
            bind->limits[LIMIT_SEQUENCE], length, tx);
        return 0;
    }
    /* CDS alignment uses replay's canonical spelling; reference protein
     * preparation retains the model bytes for Ensembl's exact stop convention. */
    for (size_t i = 0u; i < length; i++) {
        char base = duckvep_dna_normalize((char)leaf->reference_cds[i], 1);
        if (!base) {
            duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid reference CDS base");
            return 0;
        }
        s->difference_reference[i] = (uint8_t)base;
    }
    duckvep_codon_table_t table = seq->codon_table
        ? (duckvep_codon_table_t)seq->codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
    duckvep_translation_status_t coding = duckvep_translate_cds(leaf->reference_cds,
        length, table, DUCKVEP_TRANSLATION_N_UNKNOWN, s->reference_coding_protein,
        s->buffers.reference_protein_capacity, &s->reference_coding_translation);
    if (coding != DUCKVEP_TRANSLATION_OK) {
        snprintf(error, error_size, "duckvep_haplotypes: reference coding translation status %u at transcript %u",
            (unsigned)coding, tx);
        return 0;
    }
    s->difference_transcript = tx; s->have_difference_reference = 1;
    return 1;
}

/* Both sequence axes reuse the same bounded traceback and descriptor storage.
 * DuckDB copies one list's spans before the next axis resets those descriptors. */
static int append_sequence_differences(duckdb_vector vector, idx_t row, haplotype_state_t *s,
    const haplotype_bind_t *bind, const duckvep_haplotype_leaf_t *leaf, int protein,
    char *error, size_t error_size) {
    int known = leaf->cds && (!protein || leaf->reference_protein);
    const uint8_t *reference = protein ? leaf->reference_protein : s->difference_reference;
    const uint8_t *alternate = protein ? leaf->protein : leaf->cds;
    duckvep_sequence_diff_result_t result = {0};
    if (known) {
        size_t ref_length = protein ? leaf->reference_protein_length
            : s->stream.sequences->cds_length[leaf->carriers.transcript_index];
        size_t alt_length = protein ? leaf->protein_length : leaf->cds_length;
        duckvep_sequence_diff_status_t status = duckvep_sequence_differences(reference, ref_length,
            alternate, alt_length, (leaf->flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u,
            &s->difference_scratch, s->differences, bind->limits[LIMIT_DIFFERENCES], &result);
        if (status != DUCKVEP_SEQUENCE_DIFF_OK) {
            int limit = status == DUCKVEP_SEQUENCE_DIFF_TRACE_FULL ? LIMIT_ALIGNMENT :
                status == DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL ? LIMIT_DIFFERENCES : LIMIT_SEQUENCE;
            size_t required = limit == LIMIT_ALIGNMENT ? result.trace_cells :
                limit == LIMIT_DIFFERENCES ? result.count : alt_length;
            snprintf(error, error_size,
                "duckvep_haplotypes: %s difference status %u, %s=%zu, required=%zu at transcript %u",
                protein ? "protein" : "CDS", (unsigned)status, limit_names[limit],
                bind->limits[limit], required, leaf->carriers.transcript_index);
            return 0;
        }
    }
    duckdb_list_entry entry;
    if (!duckhts_list_extend(vector, result.count, &entry)) return 0;
    idx_t base = entry.offset;
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
    if (!known) null_cell(vector, row);
    duckdb_vector records = duckdb_list_vector_get_child(vector), fields[5];
    duckdb_vector_ensure_validity_writable(records);
    for (unsigned j = 0u; j < 5u; j++) {
        fields[j] = duckdb_struct_vector_get_child(records, j);
        duckdb_vector_ensure_validity_writable(fields[j]);
    }
    for (size_t i = 0u; i < result.count; i++) {
        idx_t at = base + i;
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(records), at);
        for (unsigned j = 0u; j < 5u; j++)
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[j]), at);
        const duckvep_sequence_difference_t *d = &s->differences[i];
        ((uint64_t *)duckdb_vector_get_data(fields[0]))[at] = d->ref_start0;
        ((uint64_t *)duckdb_vector_get_data(fields[1]))[at] = d->alt_start0;
        duckdb_vector_assign_string_element_len(fields[2], at,
            (const char *)reference + d->ref_start0, d->ref_length);
        duckdb_vector_assign_string_element_len(fields[3], at,
            (const char *)alternate + d->alt_start0, d->alt_length);
        ((uint64_t *)duckdb_vector_get_data(fields[4]))[at] = d->alignment_start0;
    }
    return 1;
}

typedef struct {
    haplotype_state_t *state;
    const haplotype_bind_t *bind;
    duckvep_hgvs_status_t status;
    size_t length;
    char *error;
    size_t error_size;
} haplotype_hgvs_observer_t;

static int haplotype_hgvs_observe(void *pointer, const duckvep_variant_batch_t *variants,
    const duckvep_consequence_t *row, const duckvep_pair_facts_t *facts) {
    haplotype_hgvs_observer_t *o = pointer;
    haplotype_state_t *s = o->state;
    const duckvep_owned_model_t *m = &o->bind->entry->model;
    if (!facts || !facts->event) return 0;
    if (facts->transcript_edit_status != DUCKVEP_TRANSCRIPT_EDIT_OK) {
        o->status = facts->transcript_edit_status == DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT
            ? DUCKVEP_HGVS_NOT_APPLICABLE : DUCKVEP_HGVS_INVALID_PROJECTION;
        return 1;
    }
    int available;
    duckvep_hgvs_reference_window_t shift, lookup;
    if (!duckvep_reference_reader_windows(&s->reference, facts->event, &available,
            &shift, &lookup, o->error, o->error_size)) return 0;
    o->status = available ? duckvep_hgvs_uploaded_reference_validate(&lookup, facts->event,
        variants->allele_bytes + variants->ref_offset[0], variants->ref_length[0]) : DUCKVEP_HGVS_OK;
    if (o->status != DUCKVEP_HGVS_OK) return 1;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t dna;
    o->status = duckvep_hgvs_dna_pair_build(&m->transcripts, &m->exons, &m->sequences,
        variants, facts, available ? &shift : NULL, available ? &lookup : NULL,
        &s->hgvs_scratch, &edit, &dna);
    if (o->status != DUCKVEP_HGVS_OK) return 1;
    duckvep_pair_facts_t protein_facts = *facts;
    protein_facts.transcript_edit = &edit;
    duckvep_hgvs_protein_pair_t protein;
    size_t required;
    o->status = duckvep_hgvs_protein_pair_build(&m->transcripts, &m->exons, &m->sequences,
        variants, row, &protein_facts, &dna, available ? &lookup : NULL,
        &s->hgvs_scratch, s->hgvs_shifted_allele, s->hgvs_shifted_capacity, &required, &protein);
    if (o->status != DUCKVEP_HGVS_OK) return 1;
    o->status = duckvep_hgvs_protein_render(&protein.fact, 1, s->hgvsp,
        o->bind->limits[LIMIT_HGVS_BYTES] + 1u, &o->length);
    if (o->status == DUCKVEP_HGVS_BUFFER_TOO_SMALL) {
        snprintf(o->error, o->error_size,
            "duckvep_haplotypes: max_hgvs_bytes=%zu, required=%zu at transcript %u",
            o->bind->limits[LIMIT_HGVS_BYTES], o->length, row->tx_idx);
        return 0;
    }
    return 1;
}

static int append_single_event_hgvsp(duckdb_vector text, duckdb_vector status_vector, idx_t row,
    haplotype_state_t *s, const haplotype_bind_t *bind, const duckvep_haplotype_leaf_t *leaf,
    char *error, size_t error_size) {
    const duckvep_haplotype_contributor_t *contributor = &leaf->contributors[0];
    const duckvep_haplotype_source_t *source = &contributor->source;
    size_t bytes = (size_t)source->ref_len + source->alt_len;
    const duckvep_event_t *event = contributor->prepared;
    duckvep_event_t normalized = {0};
    if (source->source_record) {
        if (!duckvep_event_prepare_small(source->pos1, source->ref, source->ref_len,
                source->alt, source->alt_len, &normalized)) {
            duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid source allele for HGVS");
            return 0;
        }
        normalized.chrom_id = source->chrom_id;
        event = &normalized;
    }
    uint32_t ref_offset = 0u, alt_offset = source->ref_len;
    duckvep_variant_batch_t variant = {.chrom_id = &source->chrom_id, .pos1 = &source->pos1,
        .end1 = &event->raw_end1, .ref_offset = &ref_offset, .alt_offset = &alt_offset,
        .ref_length = &source->ref_len, .alt_length = &source->alt_len,
        .allele_bytes = source->ref, .allele_bytes_len = bytes, .variant_kind = &event->kind,
        .count = 1u};
    haplotype_hgvs_observer_t observer = {.state = s, .bind = bind,
        .status = DUCKVEP_HGVS_NOT_APPLICABLE, .error = error, .error_size = error_size};
    duckvep_error_t native_error = {0};
    if (duckvep_annotate_pair_observed(bind->entry->model.kernel, &variant, event,
            leaf->carriers.transcript_index, &s->hgvs_scratch,
            source->source_record ? NULL : contributor->projected, haplotype_hgvs_observe,
            &observer, &native_error) != DUCKVEP_OK) {
        if (!error[0]) duckvep_sql_set_error(error, error_size, native_error.message);
        return 0;
    }
    const char *name;
    switch (observer.status) {
        case DUCKVEP_HGVS_OK:
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(text), row);
            duckdb_vector_assign_string_element_len(text, row, s->hgvsp, observer.length);
            name = "ok"; break;
        case DUCKVEP_HGVS_NOT_APPLICABLE: name = "not_applicable"; break;
        case DUCKVEP_HGVS_MISSING_REFERENCE: name = "missing_reference"; break;
        case DUCKVEP_HGVS_REFERENCE_MISMATCH: name = "reference_mismatch"; break;
        case DUCKVEP_HGVS_INVALID_ALLELE: name = "invalid_allele"; break;
        case DUCKVEP_HGVS_MISSING_PEPTIDE: name = "missing_peptide"; break;
        case DUCKVEP_HGVS_MISSING_TRANSCRIPT_TAIL: name = "missing_transcript_tail"; break;
        case DUCKVEP_HGVS_MISSING_TRANSCRIPT_FLANK: name = "missing_transcript_flank"; break;
        case DUCKVEP_HGVS_UNSUPPORTED_PROTEIN: name = "unsupported_protein"; break;
        case DUCKVEP_HGVS_UNSUPPORTED_EDIT: name = "unsupported_coding_context"; break;
        case DUCKVEP_HGVS_INVALID_PROJECTION: name = "invalid_projection"; break;
        default:
            snprintf(error, error_size,
                "duckvep_haplotypes: HGVS status %u at event %llu; max_sequence_bases=%zu, max_leaf_edits=%zu, max_allele_bytes=%zu",
                (unsigned)observer.status, (unsigned long long)source->event_id,
                bind->limits[LIMIT_SEQUENCE], bind->limits[LIMIT_LEAF_EDITS], bind->limits[LIMIT_ALLELES]);
            return 0;
    }
    duckdb_vector_assign_string_element(status_vector, row, name);
    return 1;
}

static int append_hgvsp(duckdb_vector text, duckdb_vector status_vector, idx_t row,
    haplotype_state_t *s, const haplotype_bind_t *bind, const duckvep_haplotype_leaf_t *leaf,
    const duckvep_coding_context_t *coding, char *error, size_t error_size) {
    null_cell(text, row);
    const char *name = !bind->hgvs ? "not_requested" :
        !leaf->cds || !leaf->protein ? "unavailable_sequence" :
        leaf->sequence_status != DUCKVEP_HAPLOTYPE_OK ? "incomplete_input" :
        !leaf->reference_protein || !leaf->reference_protein_length ? "missing_reference_protein" :
        leaf->ordered_replacements ? "unsupported_ordered_replacements" : NULL;
    /* A singleton source allele uses the independent-event VEP presentation,
     * including its original anchor and cached predicates. Physical MNV islands
     * do not turn that one source into a compound event. */
    if (!name && leaf->contributor_count == 1u) {
        const duckvep_haplotype_contributor_t *c = &leaf->contributors[0];
        if ((!c->source.source_record || (c->source.allele_index && c->source.allele_index != UINT32_MAX)) &&
            (c->projection_status == DUCKVEP_CDS_EDIT_OK || c->projection_status == DUCKVEP_CDS_EDIT_OUT_OF_CDS))
            return append_single_event_hgvsp(text, status_vector, row, s, bind, leaf, error, error_size);
    }
    size_t count = 0u;
    /* The raw reference-only route already borrows the prepared reference.
     * No codon replay is needed to establish equality of these exact operands. */
    int reference_only = !name && !leaf->edit_count && leaf->protein == leaf->reference_protein;
    if (!name && !reference_only) {
        duckvep_hgvs_protein_reference_t reference = {
            leaf->reference_protein, leaf->reference_protein_length};
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_build(coding,
            &reference, s->buffers.edits, leaf->edit_count, leaf->blocks, leaf->block_count,
            bind->entry->model.transcripts.flags[leaf->carriers.transcript_index],
            s->protein_operations, bind->limits[LIMIT_HGVS_OPERATIONS], &count);
        if (status == DUCKVEP_HGVS_BUFFER_TOO_SMALL) {
            snprintf(error, error_size,
                "duckvep_haplotypes: max_hgvs_operations=%zu exhausted at transcript %u",
                bind->limits[LIMIT_HGVS_OPERATIONS], leaf->carriers.transcript_index);
            return 0;
        }
        switch (status) {
            case DUCKVEP_HGVS_OK: break;
            case DUCKVEP_HGVS_NOT_APPLICABLE: name = "not_applicable"; break;
            case DUCKVEP_HGVS_MISSING_PEPTIDE: name = "missing_peptide"; break;
            case DUCKVEP_HGVS_UNSUPPORTED_PROTEIN: name = "unsupported_protein"; break;
            case DUCKVEP_HGVS_UNSUPPORTED_EDIT: name = "unsupported_coding_context"; break;
            default:
                snprintf(error, error_size, "duckvep_haplotypes: HGVS status %u at transcript %u",
                    (unsigned)status, leaf->carriers.transcript_index);
                return 0;
        }
    }
    if (!name) {
        size_t required;
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_render(s->protein_operations,
            count, 1, s->hgvsp, bind->limits[LIMIT_HGVS_BYTES] + 1u, &required);
        if (status != DUCKVEP_HGVS_OK) {
            snprintf(error, error_size,
                "duckvep_haplotypes: HGVS render status %u, max_hgvs_bytes=%zu, required=%zu at transcript %u",
                (unsigned)status, bind->limits[LIMIT_HGVS_BYTES], required, leaf->carriers.transcript_index);
            return 0;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(text), row);
        duckdb_vector_assign_string_element_len(text, row, s->hgvsp, required);
        name = "ok";
    }
    duckdb_vector_assign_string_element(status_vector, row, name);
    return 1;
}

static int append_leaf(duckdb_data_chunk output, idx_t row, haplotype_state_t *s,
    const haplotype_bind_t *bind, const duckvep_haplotype_leaf_t *leaf,
    char *error, size_t error_size) {
    if (!prepare_difference_reference(s, bind, leaf, error, error_size)) return 0;
    duckvep_coding_context_t coding;
    uint32_t tx = leaf->carriers.transcript_index;
    if (!leaf->ordered_replacements && (leaf->block_count ||
            (bind->hgvs && leaf->cds && leaf->reference_protein &&
             s->reference_coding_translation.length))) {
        const duckvep_owned_model_t *model = &bind->entry->model;
        duckvep_edit_set_t edits = {s->buffers.edits, leaf->edit_count};
        duckvep_haplotype_result_t applied = {leaf->cds_length,
            (int64_t)leaf->cds_length - model->sequences.cds_length[tx],
            leaf->flags & ~(uint32_t)DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED, leaf->edit_count};
        duckvep_codon_table_t table = model->sequences.codon_table
            ? (duckvep_codon_table_t)model->sequences.codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
        const duckvep_event_t *event = NULL;
        if (leaf->edit_count == 1u) {
            for (size_t i = 0u; i < leaf->contributor_count; i++)
                if (leaf->contributors[i].source.event_id == leaf->edit_event_ids[0])
                    event = leaf->contributors[i].prepared;
        }
        if (duckvep_coding_context_open_replay(leaf->reference_cds, model->sequences.cds_length[tx],
                &edits, model->transcripts.strand[tx], table, leaf->cds, &applied,
                s->reference_coding_protein, &s->reference_coding_translation,
                s->buffers.protein, &leaf->translation, &coding) != DUCKVEP_CODING_CONTEXT_OK ||
            duckvep_coding_context_attach_model(&model->transcripts, &model->exons, &model->sequences,
                tx, event, leaf->edit_count == 1u ? edits.edits[0].cds_start : 0u, &coding) !=
                DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
            duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid completed coding context");
            return 0;
        }
    }
    duckdb_vector v[HAPLOTYPE_OUTPUT_COLUMNS];
    for (unsigned i = 0u; i < HAPLOTYPE_OUTPUT_COLUMNS; i++) {
        v[i] = duckdb_data_chunk_get_vector(output, i);
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(v[i]), row);
    }
    ((uint32_t *)duckdb_vector_get_data(v[0]))[row] = leaf->carriers.transcript_index;
    if (leaf->cds) duckdb_vector_assign_string_element_len(v[1], row, (const char *)leaf->cds, leaf->cds_length);
    else null_cell(v[1], row);
    if (leaf->protein) duckdb_vector_assign_string_element_len(v[2], row, (const char *)leaf->protein, leaf->protein_length);
    else null_cell(v[2], row);
    ((uint32_t *)duckdb_vector_get_data(v[3]))[row] = leaf->flags;
    ((uint8_t *)duckdb_vector_get_data(v[4]))[row] = leaf->evidence_flags;
    duckdb_vector_assign_string_element(v[5], row, projection_name(leaf->projection_status));
    duckdb_vector_assign_string_element(v[6], row, leaf->projection_status == DUCKVEP_CDS_EDIT_OK
        ? sequence_name(leaf->sequence_status) : "unavailable_projection");
    ((uint64_t *)duckdb_vector_get_data(v[7]))[row] = leaf->edit_count;
    ((uint32_t *)duckdb_vector_get_data(v[8]))[row] = leaf->carriers.call_count;
    if (leaf->cds && !leaf->ordered_replacements)
        ((bool *)duckdb_vector_get_data(v[HAPLOTYPE_STOP_COLUMN]))[row] = leaf->stop_in_displaced_frame != 0u;
    else null_cell(v[HAPLOTYPE_STOP_COLUMN], row);
    const size_t counts[] = {leaf->carriers.call_count, leaf->contributor_count, leaf->block_count};
    const unsigned field_counts[] = {4u, bind->source_records ? 8u : 7u, HAPLOTYPE_BLOCK_FIELDS};
    for (unsigned list = 0u; list < 3u; list++) {
        duckdb_vector vector = v[HAPLOTYPE_LIST_COLUMN + list];
        size_t count = counts[list];
        duckdb_list_entry entry;
        if (!duckhts_list_extend(vector, count, &entry)) return 0;
        idx_t base = entry.offset;
        ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
        if (list >= 2u && !leaf->cds) null_cell(vector, row);
        duckdb_vector records = duckdb_list_vector_get_child(vector), fields[HAPLOTYPE_BLOCK_FIELDS];
        duckdb_vector_ensure_validity_writable(records);
        for (unsigned j = 0u; j < field_counts[list]; j++) {
            fields[j] = duckdb_struct_vector_get_child(records, j);
            duckdb_vector_ensure_validity_writable(fields[j]);
        }
        idx_t event_base = 0u;
        if (list == 2u && leaf->cds) {
            duckdb_vector event_vector = fields[HAPLOTYPE_BLOCK_EVENT_FIELD];
            duckdb_list_entry events;
            if (!duckhts_list_extend(event_vector, leaf->edit_count, &events)) return 0;
            event_base = events.offset;
            duckdb_vector ids = duckdb_list_vector_get_child(event_vector);
            duckdb_vector_ensure_validity_writable(ids);
            uint64_t *data = duckdb_vector_get_data(ids);
            for (size_t i = 0u; i < leaf->edit_count; i++) {
                data[event_base + i] = leaf->edit_event_ids[i];
                duckdb_validity_set_row_valid(duckdb_vector_get_validity(ids), event_base + i);
            }
        }
        uint32_t call_id = leaf->carriers.first_call;
        for (size_t i = 0u; i < count; i++) {
            idx_t at = base + i;
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(records), at);
            for (unsigned j = 0u; j < field_counts[list]; j++)
                duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[j]), at);
            if (!list) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->stream.carriers, call_id);
                if (!call) return 0;
                ((uint32_t *)duckdb_vector_get_data(fields[0]))[at] = call->key.sample_index;
                ((int64_t *)duckdb_vector_get_data(fields[1]))[at] = call->key.phase_set;
                if (!call->key.phase_set_present) null_cell(fields[1], at);
                ((uint16_t *)duckdb_vector_get_data(fields[2]))[at] = call->key.lane;
                ((uint16_t *)duckdb_vector_get_data(fields[3]))[at] = call->key.ploidy;
                call_id = call->next_leaf;
            } else if (list == 1u) {
                const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
                ((uint64_t *)duckdb_vector_get_data(fields[0]))[at] = c->source.event_id;
                ((uint32_t *)duckdb_vector_get_data(fields[1]))[at] = c->source.chrom_id;
                ((uint64_t *)duckdb_vector_get_data(fields[2]))[at] = c->source.pos1;
                duckdb_vector_assign_string_element_len(fields[3], at, (const char *)c->source.ref, c->source.ref_len);
                duckdb_vector_assign_string_element_len(fields[4], at, (const char *)c->source.alt, c->source.alt_len);
                ((uint8_t *)duckdb_vector_get_data(fields[5]))[at] = c->evidence_flags;
                duckdb_vector_assign_string_element(fields[6], at, projection_name(c->projection_status));
                if (bind->source_records) {
                    if (c->source.allele_index == UINT32_MAX) null_cell(fields[7], at);
                    else ((uint32_t *)duckdb_vector_get_data(fields[7]))[at] = c->source.allele_index;
                }
            } else if (list == 2u) {
                const duckvep_haplotype_block_t *block = &leaf->blocks[i];
                ((uint32_t *)duckdb_vector_get_data(fields[0]))[at] = block->cds_start;
                duckdb_vector_assign_string_element_len(fields[1], at,
                    (const char *)leaf->reference_cds + block->cds_start - 1u, block->ref_len);
                duckdb_vector_assign_string_element_len(fields[2], at,
                    (const char *)leaf->cds + block->alt_start0, block->alt_len);
                ((uint64_t *)duckdb_vector_get_data(fields[3]))[at] = block->alt_start0;
                ((int64_t *)duckdb_vector_get_data(fields[4]))[at] = block->length_diff;
                ((uint32_t *)duckdb_vector_get_data(fields[5]))[at] = block->flags;
                duckvep_sequence_delta_t delta;
                duckvep_context_delta_status_t status = leaf->ordered_replacements
                    ? DUCKVEP_CONTEXT_DELTA_UNSUPPORTED : duckvep_coding_context_block_delta_fill(
                    &coding, s->buffers.edits, leaf->edit_count, block,
                    bind->entry->model.transcripts.flags[tx], &delta);
                const char *name = leaf->ordered_replacements ? "unsupported_ordered_replacements" :
                    status == DUCKVEP_CONTEXT_DELTA_OK ? "ok" :
                    status == DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL ? "missing_transcript_tail" :
                    status == DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK ? "missing_transcript_flank" :
                    status == DUCKVEP_CONTEXT_DELTA_UNSUPPORTED ? "unsupported" : "invalid_argument";
                duckdb_vector_assign_string_element(fields[6], at, name);
                if (status == DUCKVEP_CONTEXT_DELTA_OK)
                    ((uint64_t *)duckdb_vector_get_data(fields[7]))[at] = duckvep_effect_eval_coding_delta(&delta);
                else null_cell(fields[7], at);
                ((bool *)duckdb_vector_get_data(fields[8]))[at] =
                    leaf->translation.first_stop_position1 &&
                    block->alt_start0 / 3u >= leaf->translation.first_stop_position1;
                ((duckdb_list_entry *)duckdb_vector_get_data(fields[HAPLOTYPE_BLOCK_EVENT_FIELD]))[at] =
                    (duckdb_list_entry){event_base + block->edit_begin, block->edit_count};
            }
        }
    }
    return append_sequence_differences(v[12], row, s, bind, leaf, 0, error, error_size) &&
        append_sequence_differences(v[13], row, s, bind, leaf, 1, error, error_size) &&
        append_hgvsp(v[HAPLOTYPE_HGVSP_COLUMN], v[HAPLOTYPE_HGVSP_STATUS_COLUMN], row,
            s, bind, leaf, &coding, error, error_size);
}

static duckvep_haplotype_stream_status_t consume_call(haplotype_state_t *s,
    const haplotype_bind_t *bind, char *error, size_t error_size) {
    duckdb_vector v[17];
    idx_t row = s->row;
    for (unsigned i = 0u; i < (bind->source_records ? 17u : 15u); i++) {
        v[i] = duckdb_data_chunk_get_vector(s->chunk, i);
        if (i != 9u && i != 10u && duckvep_row_is_null(v[i], row)) {
            snprintf(error, error_size, "duckvep_haplotypes: required input column %u is NULL", i + 1u);
            return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        }
    }
    if (((int64_t *)duckdb_vector_get_data(v[12]))[row] != 1 ||
        ((int64_t *)duckdb_vector_get_data(v[13]))[row] != 1 ||
        ((int64_t *)duckdb_vector_get_data(v[14]))[row] != 1) {
        duckvep_sql_set_error(error, error_size, bind->source_records
            ? "duckvep_haplotypes: duplicate call, inconsistent source record identity or source GT"
            : "duckvep_haplotypes: duplicate call, inconsistent event identity or changed sample/transcript ploidy");
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    }
    uint32_t chrom = ((uint32_t *)duckdb_vector_get_data(v[1]))[row];
    uint64_t pos = ((uint64_t *)duckdb_vector_get_data(v[2]))[row];
    duckdb_string_t ref = ((duckdb_string_t *)duckdb_vector_get_data(v[3]))[row];
    duckdb_string_t alt = ((duckdb_string_t *)duckdb_vector_get_data(v[4]))[row];
    uint32_t ref_len = duckdb_string_t_length(ref), alt_len = duckdb_string_t_length(alt);
    uint32_t allele_index = ((uint32_t *)duckdb_vector_get_data(v[5]))[row];
    if (chrom > UINT16_MAX || !pos || pos > UINT32_MAX || !ref_len || ref_len > UINT16_MAX ||
        (!alt_len && !(bind->source_records && allele_index == UINT32_MAX)) || alt_len > UINT16_MAX)
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    duckvep_haplotype_source_t source = {((uint64_t *)duckdb_vector_get_data(v[0]))[row],
        (const uint8_t *)duckdb_string_t_data(&ref), (const uint8_t *)duckdb_string_t_data(&alt),
        (uint32_t)pos, (uint16_t)chrom, (uint16_t)ref_len, (uint16_t)alt_len,
        bind->source_records ? allele_index : 0u, (uint8_t)bind->source_records,
        bind->source_records ? ((uint64_t *)duckdb_vector_get_data(v[15]))[row] : 0u};
    uint32_t tx = ((uint32_t *)duckdb_vector_get_data(v[6]))[row];
    duckvep_haplotype_stream_status_t status;
    int new_event = !s->stream.have_input || source.event_id != s->stream.last_event_id ||
        source.allele_index != s->stream.last_allele_index;
    if (new_event) {
        status = duckvep_haplotype_stream_begin(&s->stream, &source);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        s->have_call = 0;
    }
    if (!s->have_call || tx != s->last_tx) {
        status = duckvep_haplotype_stream_project(&s->stream, tx);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    if (bind->source_records) {
        uint32_t fields[7];
        for (unsigned i = 0u; i < 7u; i++) fields[i] =
            ((uint32_t *)duckdb_vector_get_data(duckdb_struct_vector_get_child(v[8], i)))[row];
        duckvep_raw_gt_status_t parsed = (duckvep_raw_gt_status_t)fields[0];
        duckvep_raw_gt_t call = {{fields[1], fields[2]}, fields[3], (uint16_t)fields[4],
            (uint8_t)fields[5], (duckvep_raw_gt_disposition_t)fields[6]};
        if (parsed != DUCKVEP_RAW_GT_OK) {
            snprintf(error, error_size, "duckvep_haplotypes: raw GT status %u at event %llu, sample %u",
                (unsigned)parsed, (unsigned long long)source.event_id,
                ((uint32_t *)duckdb_vector_get_data(v[7]))[row]);
            return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        }
        if (call.source_ploidy > bind->limits[LIMIT_PLOIDY] || bind->limits[LIMIT_PLOIDY] < 2u) {
            snprintf(error, error_size, "duckvep_haplotypes: max_ploidy=%zu exceeded by raw source/file ploidy",
                bind->limits[LIMIT_PLOIDY]);
            return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        }
        status = duckvep_haplotype_stream_push_raw_call(&s->stream, tx,
            ((uint32_t *)duckdb_vector_get_data(v[7]))[row], &call,
            (uint8_t)((bool *)duckdb_vector_get_data(v[16]))[row]);
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            s->have_call = 1; s->last_tx = tx; s->row++;
        }
        return status;
    }
    duckdb_list_entry gt = ((duckdb_list_entry *)duckdb_vector_get_data(v[8]))[row];
    int have_phase = !duckvep_row_is_null(v[9], row);
    duckdb_list_entry phases = have_phase ? ((duckdb_list_entry *)duckdb_vector_get_data(v[9]))[row] : (duckdb_list_entry){0};
    duckdb_list_entry sets = ((duckdb_list_entry *)duckdb_vector_get_data(v[11]))[row];
    if (!gt.length || gt.length > bind->limits[LIMIT_PLOIDY] ||
        (have_phase && phases.length != gt.length) || sets.length > bind->limits[LIMIT_PHASE_SETS]) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid GT/phase lengths or max_ploidy/max_phase_sets exceeded");
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    }
    duckdb_vector g = duckdb_list_vector_get_child(v[8]), p = duckdb_list_vector_get_child(v[9]);
    for (idx_t i = 0u; i < gt.length; i++) {
        int missing = duckvep_row_is_null(g, gt.offset + i);
        s->gt[i] = missing ? -1 : ((int32_t *)duckdb_vector_get_data(g))[gt.offset + i];
        if (!missing && s->gt[i] < 0) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        s->phase[i] = have_phase && !duckvep_row_is_null(p, phases.offset + i) &&
            ((bool *)duckdb_vector_get_data(p))[phases.offset + i];
    }
    duckdb_vector domains = duckdb_list_vector_get_child(v[11]);
    for (idx_t i = 0u; i < sets.length; i++) {
        s->sets[i].present = !duckvep_row_is_null(domains, sets.offset + i);
        s->sets[i].value = s->sets[i].present ? ((int64_t *)duckdb_vector_get_data(domains))[sets.offset + i] : 0;
    }
    duckvep_haplotype_call_t call = {s->gt, s->phase,
        ((uint32_t *)duckdb_vector_get_data(v[7]))[row], ((uint32_t *)duckdb_vector_get_data(v[5]))[row],
        (uint16_t)gt.length, {0, 0u}, bind->policy};
    call.phase_set.present = !duckvep_row_is_null(v[10], row);
    if (call.phase_set.present) call.phase_set.value = ((int64_t *)duckdb_vector_get_data(v[10]))[row];
    status = duckvep_haplotype_stream_push_call(&s->stream, tx, &call, s->sets, sets.length);
    if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        s->have_call = 1; s->last_tx = tx;
        s->row++;
    }
    return status;
}

static void haplotype_scan(duckdb_function_info info, duckdb_data_chunk output) {
    const haplotype_bind_t *bind = duckdb_function_get_bind_data(info);
    haplotype_state_t *s = duckdb_function_get_init_data(info);
    idx_t rows = 0u, capacity = duckdb_vector_size();
    for (unsigned i = 0u; i < HAPLOTYPE_OUTPUT_COLUMNS; i++)
        duckdb_vector_ensure_validity_writable(duckdb_data_chunk_get_vector(output, i));
    for (unsigned i = HAPLOTYPE_LIST_COLUMN; i < HAPLOTYPE_STOP_COLUMN; i++)
        if (duckdb_list_vector_set_size(duckdb_data_chunk_get_vector(output, i), 0u) != DuckDBSuccess) {
            duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset output list"); return;
        }
    duckdb_vector blocks = duckdb_list_vector_get_child(duckdb_data_chunk_get_vector(output, 11u));
    if (duckdb_list_vector_set_size(duckdb_struct_vector_get_child(blocks, HAPLOTYPE_BLOCK_EVENT_FIELD), 0u) != DuckDBSuccess) {
        duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset block event list"); return;
    }
    char error[DUCKVEP_SQL_ERROR_SIZE] = {0};
    while (rows < capacity) {
        duckvep_haplotype_stream_status_t status;
        if (s->stream.closing) {
            duckvep_haplotype_leaf_t leaf;
            status = duckvep_haplotype_stream_next(&s->stream, &leaf);
            if (status == DUCKVEP_HAPLOTYPE_STREAM_DONE) continue;
            if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                if (!append_leaf(output, rows, s, bind, &leaf, error, sizeof(error))) {
                    duckdb_function_set_error(info, error[0] ? error :
                        "duckvep_haplotypes: output list allocation failed"); return;
                }
                rows++; continue;
            }
        } else {
            if (!s->eof && (!s->chunk || s->row == duckdb_data_chunk_get_size(s->chunk))) {
                if (s->chunk) duckdb_destroy_data_chunk(&s->chunk);
                s->chunk = duckdb_fetch_chunk(s->input); s->row = 0u;
                if (!s->chunk) s->eof = 1;
                else if (!duckdb_data_chunk_get_size(s->chunk)) continue;
            }
            status = s->eof ? duckvep_haplotype_stream_finish(&s->stream) : consume_call(s, bind, error, sizeof(error));
            if (status == DUCKVEP_HAPLOTYPE_STREAM_DONE) break;
            if (status == DUCKVEP_HAPLOTYPE_STREAM_OK || status == DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) continue;
        }
        if (!error[0]) {
            uint64_t event = s->stream.last_event_id, position = s->stream.last_pos1;
            if (!s->eof && s->chunk && s->row < duckdb_data_chunk_get_size(s->chunk)) {
                event = ((uint64_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(s->chunk, 0u)))[s->row];
                position = ((uint64_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(s->chunk, 2u)))[s->row];
            }
            int limit = exhausted_limit(status, s->stream.carrier_error);
            if (limit >= 0) snprintf(error, sizeof(error),
                "duckvep_haplotypes: %s=%zu exhausted at event %llu, position %llu",
                limit_names[limit], bind->limits[limit], (unsigned long long)event, (unsigned long long)position);
            else snprintf(error, sizeof(error),
                "duckvep_haplotypes: native status %u, carrier status %u at event %llu, position %llu; invalid input or candidate",
                (unsigned)status, (unsigned)s->stream.carrier_error,
                (unsigned long long)event, (unsigned long long)position);
        }
        duckdb_function_set_error(info, error); return;
    }
    duckdb_data_chunk_set_size(output, rows);
}

void duckvep_register_haplotypes(duckdb_connection connection, duckvep_registry_t *registry) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type string = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_set_name(function, "duckvep_haplotypes");
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_named_parameter(function, "phase_policy", string);
    duckdb_table_function_add_named_parameter(function, "input_mode", string);
    duckdb_table_function_add_named_parameter(function, "hgvs", boolean);
    for (unsigned i = 0u; i < LIMIT_COUNT; i++)
        duckdb_table_function_add_named_parameter(function, limit_names[i], integer);
    duckvep_registry_retain(registry);
    duckdb_table_function_set_extra_info(function, registry, duckvep_registry_release);
    duckdb_table_function_set_bind(function, haplotype_bind);
    duckdb_table_function_set_init(function, haplotype_init);
    duckdb_table_function_set_function(function, haplotype_scan);
    (void)duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&string); duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&boolean);
}
