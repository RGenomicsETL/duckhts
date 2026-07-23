/* Stable-ABI DuckDB vector adapter for the resident DuckVEP kernel. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckvep_model.h"
#include "kernel/src/duckvep_delta.h"
#include "kernel/src/duckvep_annotation_internal.h"
#include "kernel/src/duckvep_effect.h"
#include "kernel/src/duckvep_event.h"
#include "kernel/src/duckvep_hgvs.h"
#include "kernel/src/duckvep_projection.h"
#include "kernel/src/duckvep_transcript_edit.h"
#include "kernel/src/duckvep_workspace_internal.h"
#include "kernel/include/duckvep_kernel.h"
#include "kernel/include/duckvep_so.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DUCKVEP_REFERENCE_READ_AHEAD 65536u
#define DUCKVEP_HGVS_INITIAL_RENDER_CAPACITY 256u

typedef enum duckvep_hgvs_adapter_reason {
	DUCKVEP_HGVS_ADAPTER_NONE = 0,
	DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE,
	DUCKVEP_HGVS_ADAPTER_MISSING_REFERENCE,
	DUCKVEP_HGVS_ADAPTER_MISSING_SEQUENCE,
	DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION,
	DUCKVEP_HGVS_ADAPTER_REFERENCE_MISMATCH,
	DUCKVEP_HGVS_ADAPTER_INVALID_ALLELE,
	DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT,
	DUCKVEP_HGVS_ADAPTER_INTERNAL_CAPACITY,
	DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_TAIL,
	DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_FLANK,
	DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_PROTEIN
} duckvep_hgvs_adapter_reason_t;

typedef struct duckvep_hgvs_scalar_result {
	uint32_t transcript_offset;
	uint32_t transcript_length;
	uint32_t protein_offset;
	uint32_t protein_length;
	uint32_t shift_offset;
	uint8_t transcript_reason;
	uint8_t protein_reason;
} duckvep_hgvs_scalar_result_t;

typedef struct duckvep_scalar_state {
	duckvep_registry_t *registry;
	duckvep_model_entry_t *entry;
	duckvep_workspace_cache_t *workspace_cache;
	duckvep_workspace_t *workspace;
	duckvep_options_t *options;
	uint32_t options_upstream_distance;
	uint32_t options_downstream_distance;
	int have_options_distances;
	int64_t *interval_hits;
	int64_t interval_hit_capacity;
	uint32_t *seed_transcripts;
	size_t seed_capacity;
	uint16_t *seq_regions;
	uint32_t *positions;
	uint32_t *ends;
	uint16_t *mate_seq_regions;
	uint32_t *mate_positions;
	uint32_t *reference_offsets;
	uint16_t *reference_lengths;
	uint32_t *alternate_offsets;
	uint16_t *alternate_lengths;
	uint8_t *variant_kinds;
	uint8_t *sv_types;
	uint8_t *copy_changes;
	uint8_t *transcript_coverage_complete;
	uint8_t *allele_bytes;
	size_t variant_capacity;
	size_t allele_capacity;
	uint32_t *pair_variant_indices;
	uint32_t *pair_object_indices;
	size_t pair_capacity;
	duckvep_consequence_t *results;
	size_t result_count;
	size_t result_capacity;
	duckvep_consequence_t *result_merge;
	size_t result_merge_capacity;
	duckvep_hgvs_scalar_result_t *hgvs_results;
	size_t hgvs_result_capacity;
	uint8_t *hgvs_allele_scratch;
	size_t hgvs_allele_capacity;
	char *hgvs_render_scratch;
	size_t hgvs_render_capacity;
	char *hgvs_text;
	size_t hgvs_text_size;
	size_t hgvs_text_capacity;
	struct duckvep_scalar_state *next_free;
} duckvep_scalar_state_t;

static int
duckvep_scalar_variant_reserve(duckvep_scalar_state_t *state, size_t needed)
{
	size_t capacity;

	if (needed <= state->variant_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->variant_capacity, needed);
#define DUCKVEP_SCALAR_VARIANT_RESIZE(member) \
	if (!duckvep_sql_resize((void **)&state->member, \
	    sizeof(*state->member), capacity)) \
		return 0
	DUCKVEP_SCALAR_VARIANT_RESIZE(seq_regions);
	DUCKVEP_SCALAR_VARIANT_RESIZE(positions);
	DUCKVEP_SCALAR_VARIANT_RESIZE(ends);
	DUCKVEP_SCALAR_VARIANT_RESIZE(mate_seq_regions);
	DUCKVEP_SCALAR_VARIANT_RESIZE(mate_positions);
	DUCKVEP_SCALAR_VARIANT_RESIZE(reference_offsets);
	DUCKVEP_SCALAR_VARIANT_RESIZE(reference_lengths);
	DUCKVEP_SCALAR_VARIANT_RESIZE(alternate_offsets);
	DUCKVEP_SCALAR_VARIANT_RESIZE(alternate_lengths);
	DUCKVEP_SCALAR_VARIANT_RESIZE(variant_kinds);
	DUCKVEP_SCALAR_VARIANT_RESIZE(sv_types);
	DUCKVEP_SCALAR_VARIANT_RESIZE(copy_changes);
	DUCKVEP_SCALAR_VARIANT_RESIZE(transcript_coverage_complete);
#undef DUCKVEP_SCALAR_VARIANT_RESIZE
	state->variant_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_pair_reserve(duckvep_scalar_state_t *state, size_t needed)
{
	size_t capacity;

	if (needed <= state->pair_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->pair_capacity, needed);
	if (!duckvep_sql_resize((void **)&state->pair_variant_indices,
	    sizeof(*state->pair_variant_indices), capacity) ||
	    !duckvep_sql_resize((void **)&state->pair_object_indices,
	    sizeof(*state->pair_object_indices), capacity))
		return 0;
	state->pair_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_allele_reserve(duckvep_scalar_state_t *state, size_t needed)
{
	size_t capacity;

	if (needed <= state->allele_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->allele_capacity, needed);
	if (!duckvep_sql_resize((void **)&state->allele_bytes,
	    sizeof(*state->allele_bytes), capacity))
		return 0;
	state->allele_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_seed_reserve(duckvep_scalar_state_t *state, size_t needed)
{
	size_t capacity;

	if (needed <= state->seed_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->seed_capacity, needed);
	if (!duckvep_sql_resize((void **)&state->seed_transcripts,
	    sizeof(*state->seed_transcripts), capacity))
		return 0;
	state->seed_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_result_reserve(duckvep_scalar_state_t *state, size_t needed)
{
	size_t capacity;

	if (needed <= state->result_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->result_capacity, needed);
	if (!duckvep_sql_resize((void **)&state->results,
	    sizeof(*state->results), capacity))
		return 0;
	state->result_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_result_merge_reserve(duckvep_scalar_state_t *state,
	size_t needed)
{
	size_t capacity;

	if (needed <= state->result_merge_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->result_merge_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&state->result_merge,
	    sizeof(*state->result_merge), capacity))
		return 0;
	state->result_merge_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_hgvs_result_reserve(duckvep_scalar_state_t *state,
	size_t needed)
{
	size_t capacity;

	if (needed <= state->hgvs_result_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->hgvs_result_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&state->hgvs_results,
	    sizeof(*state->hgvs_results), capacity))
		return 0;
	state->hgvs_result_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_hgvs_allele_reserve(duckvep_scalar_state_t *state,
	size_t needed)
{
	size_t capacity;

	if (needed <= state->hgvs_allele_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->hgvs_allele_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&state->hgvs_allele_scratch,
	    sizeof(*state->hgvs_allele_scratch), capacity))
		return 0;
	state->hgvs_allele_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_hgvs_render_reserve(duckvep_scalar_state_t *state,
	size_t needed)
{
	size_t capacity;

	if (needed <= state->hgvs_render_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(state->hgvs_render_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&state->hgvs_render_scratch,
	    sizeof(*state->hgvs_render_scratch), capacity))
		return 0;
	state->hgvs_render_capacity = capacity;
	return 1;
}

static int
duckvep_scalar_hgvs_text_append(duckvep_scalar_state_t *state,
	const char *text, size_t length, uint32_t *offset_out)
{
	size_t needed, capacity;

	if (state == NULL || text == NULL || offset_out == NULL ||
	    length > UINT32_MAX || state->hgvs_text_size > UINT32_MAX ||
	    length > UINT32_MAX - state->hgvs_text_size)
		return 0;
	needed = state->hgvs_text_size + length;
	if (needed > state->hgvs_text_capacity) {
		capacity = duckvep_sql_next_capacity(state->hgvs_text_capacity,
		    needed);
		if (!duckvep_sql_resize((void **)&state->hgvs_text,
		    sizeof(*state->hgvs_text), capacity))
			return 0;
		state->hgvs_text_capacity = capacity;
	}
	*offset_out = (uint32_t)state->hgvs_text_size;
	memcpy(state->hgvs_text + state->hgvs_text_size, text, length);
	state->hgvs_text_size = needed;
	return 1;
}

static void
duckvep_scalar_state_destroy(void *pointer)
{
	duckvep_scalar_state_t *state;

	state = pointer;
	if (state == NULL)
		return;
	if (state->workspace_cache != NULL)
		duckvep_workspace_cache_destroy(state->workspace_cache);
	else
		duckvep_workspace_close(state->workspace);
	duckvep_options_close(state->options);
	duckvep_registry_unpin(state->registry, state->entry);
	free(state->interval_hits);
	free(state->seed_transcripts);
	free(state->seq_regions);
	free(state->positions);
	free(state->ends);
	free(state->mate_seq_regions);
	free(state->mate_positions);
	free(state->reference_offsets);
	free(state->reference_lengths);
	free(state->alternate_offsets);
	free(state->alternate_lengths);
	free(state->variant_kinds);
	free(state->sv_types);
	free(state->copy_changes);
	free(state->transcript_coverage_complete);
	free(state->allele_bytes);
	free(state->pair_variant_indices);
	free(state->pair_object_indices);
	free(state->results);
	free(state->result_merge);
	free(state->hgvs_results);
	free(state->hgvs_allele_scratch);
	free(state->hgvs_render_scratch);
	free(state->hgvs_text);
	free(state);
}

static void
duckvep_scalar_state_pool_destroy(void *pointer)
{
	duckvep_scalar_state_t *state, *next;

	for (state = pointer; state != NULL; state = next) {
		next = state->next_free;
		state->next_free = NULL;
		duckvep_scalar_state_destroy(state);
	}
}

static duckvep_scalar_state_t *
duckvep_scalar_state_acquire(duckvep_registry_t *registry)
{
	duckvep_scalar_state_t *state;

	pthread_mutex_lock(&registry->mutex);
	state = registry->annotation_state_pool;
	if (state != NULL) {
		registry->annotation_state_pool = state->next_free;
		state->next_free = NULL;
	}
	if (registry->annotation_state_pool_destroy == NULL)
		registry->annotation_state_pool_destroy =
		    duckvep_scalar_state_pool_destroy;
	pthread_mutex_unlock(&registry->mutex);
	if (state == NULL)
		state = calloc(1, sizeof(*state));
	if (state != NULL)
		state->registry = registry;
	return state;
}

static void
duckvep_scalar_release_model(duckvep_scalar_state_t *state)
{
	if (state->workspace_cache != NULL) {
		duckvep_registry_workspace_return(state->registry, state->entry,
		    state->workspace_cache);
		state->workspace_cache = NULL;
		state->workspace = NULL;
	}
	duckvep_registry_unpin(state->registry, state->entry);
	state->entry = NULL;
}

static void
duckvep_scalar_state_release(duckvep_scalar_state_t *state)
{
	duckvep_registry_t *registry;

	if (state == NULL)
		return;
	registry = state->registry;
	duckvep_scalar_release_model(state);
	state->result_count = 0;
	state->hgvs_text_size = 0;
	pthread_mutex_lock(&registry->mutex);
	state->next_free = registry->annotation_state_pool;
	registry->annotation_state_pool = state;
	pthread_mutex_unlock(&registry->mutex);
}

static int
duckvep_validity_is_null(const uint64_t *validity, idx_t row)
{
	return validity != NULL &&
	    ((validity[row / 64] >> (row % 64)) & UINT64_C(1)) == 0;
}

static int
duckvep_vector_string_equals(duckdb_vector vector,
	const uint64_t *validity, idx_t row, const char *text)
{
	duckdb_string_t *strings;
	size_t length;

	if (text == NULL || duckvep_validity_is_null(validity, row))
		return 0;
	strings = duckdb_vector_get_data(vector);
	length = strlen(text);
	return length == (size_t)duckdb_string_t_length(strings[row]) &&
	    memcmp(text, duckdb_string_t_data(&strings[row]), length) == 0;
}

static int
duckvep_scalar_select_model(duckvep_scalar_state_t *state,
	duckdb_vector model_vector, idx_t row, char *error, size_t error_size)
{
	duckvep_model_entry_t *entry;
	char *name;

	name = duckvep_vector_string(model_vector, row);
	if (name == NULL || *name == '\0') {
		free(name);
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: model name must be non-empty");
		return 0;
	}
	if (state->entry != NULL && strcmp(state->entry->name, name) == 0) {
		free(name);
		return 1;
	}
	entry = duckvep_registry_pin(state->registry, name);
	free(name);
	if (entry == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: unknown model name");
		return 0;
	}
	duckvep_scalar_release_model(state);
	state->entry = entry;
	return 1;
}

static int
duckvep_scalar_simple_kind(uint32_t position,
	const uint8_t *reference, uint16_t reference_length,
	const uint8_t *alternate, uint16_t alternate_length, uint8_t *kind)
{
	duckvep_event_t event;

	if (reference == NULL || alternate == NULL || kind == NULL ||
	    reference_length == 0 || alternate_length == 0)
		return 0;
	if (!duckvep_event_prepare_small(position, reference, reference_length,
	    alternate, alternate_length, &event))
		return 0;
	*kind = event.kind;
	return 1;
}

enum duckvep_allele_geometry_field {
	DUCKVEP_GEOMETRY_KIND_CODE = 0,
	DUCKVEP_GEOMETRY_INTERBASE,
	DUCKVEP_GEOMETRY_ANCHOR_SIDE_CODE,
	DUCKVEP_GEOMETRY_RAW_START0,
	DUCKVEP_GEOMETRY_RAW_END0,
	DUCKVEP_GEOMETRY_FEATURE_START0,
	DUCKVEP_GEOMETRY_FEATURE_END0,
	DUCKVEP_GEOMETRY_EDIT_START0,
	DUCKVEP_GEOMETRY_EDIT_END0,
	DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0,
	DUCKVEP_GEOMETRY_REF_DIFF_OFFSET,
	DUCKVEP_GEOMETRY_REF_DIFF_LENGTH,
	DUCKVEP_GEOMETRY_ALT_DIFF_OFFSET,
	DUCKVEP_GEOMETRY_ALT_DIFF_LENGTH,
	DUCKVEP_GEOMETRY_FIELD_COUNT
};

static const char *duckvep_allele_geometry_field_names[] = {
	"kind_code", "interbase", "anchor_side_code", "raw_start0",
	"raw_end0", "feature_start0", "feature_end0", "edit_start0",
	"edit_end0", "insertion_boundary0", "reference_difference_offset",
	"reference_difference_length", "alternate_difference_offset",
	"alternate_difference_length"
};

static int
duckvep_scalar_dna_valid(const char *sequence, size_t length)
{
	size_t index;

	if (sequence == NULL || length == 0 || length > UINT16_MAX)
		return 0;
	for (index = 0; index < length; index++) {
		unsigned char base = (unsigned char)sequence[index];
		if (base >= 'a' && base <= 'z')
			base = (unsigned char)(base - ('a' - 'A'));
		if (base != 'A' && base != 'C' && base != 'G' &&
		    base != 'T' && base != 'N')
			return 0;
	}
	return 1;
}

static void
duckvep_allele_geometry_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckdb_vector position_vector, reference_vector, alternate_vector;
	duckdb_vector fields[DUCKVEP_GEOMETRY_FIELD_COUNT];
	uint64_t *positions;
	duckdb_string_t *references, *alternates;
	idx_t rows, row;
	size_t field;

	position_vector = duckdb_data_chunk_get_vector(input, 0);
	reference_vector = duckdb_data_chunk_get_vector(input, 1);
	alternate_vector = duckdb_data_chunk_get_vector(input, 2);
	positions = duckdb_vector_get_data(position_vector);
	references = duckdb_vector_get_data(reference_vector);
	alternates = duckdb_vector_get_data(alternate_vector);
	for (field = 0; field < DUCKVEP_GEOMETRY_FIELD_COUNT; field++)
		fields[field] = duckdb_struct_vector_get_child(output, (idx_t)field);
	duckdb_vector_ensure_validity_writable(output);
	duckdb_vector_ensure_validity_writable(
	    fields[DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0]);
	rows = duckdb_data_chunk_get_size(input);
	for (row = 0; row < rows; row++) {
		duckvep_event_t event;
		const char *reference, *alternate;
		size_t reference_length, alternate_length;

		if (duckvep_validity_is_null(
		        duckdb_vector_get_validity(position_vector), row) ||
		    duckvep_validity_is_null(
		        duckdb_vector_get_validity(reference_vector), row) ||
		    duckvep_validity_is_null(
		        duckdb_vector_get_validity(alternate_vector), row)) {
			duckdb_validity_set_row_invalid(
			    duckdb_vector_get_validity(output), row);
			continue;
		}
		reference = duckdb_string_t_data(&references[row]);
		alternate = duckdb_string_t_data(&alternates[row]);
		reference_length = (size_t)duckdb_string_t_length(references[row]);
		alternate_length = (size_t)duckdb_string_t_length(alternates[row]);
		if (positions[row] == 0 || positions[row] > UINT32_MAX ||
		    !duckvep_scalar_dna_valid(reference, reference_length) ||
		    !duckvep_scalar_dna_valid(alternate, alternate_length) ||
		    !duckvep_event_prepare_small(
		        (uint32_t)positions[row], (const uint8_t *)reference,
		        (uint16_t)reference_length, (const uint8_t *)alternate,
		        (uint16_t)alternate_length, &event)) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_allele_geometry: position must fit UINTEGER and REF/ALT must be distinct non-empty A/C/G/T/N alleles of at most 65,535 bases");
			return;
		}
		((uint8_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_KIND_CODE]))[row] = event.kind;
		((bool *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_INTERBASE]))[row] = event.interbase != 0u;
		((uint8_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_ANCHOR_SIDE_CODE]))[row] = event.anchor_side;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_RAW_START0]))[row] =
		    (uint64_t)event.raw_start1 - 1u;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_RAW_END0]))[row] = event.raw_end1;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_FEATURE_START0]))[row] =
		    (uint64_t)event.feature_start1 - 1u;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_FEATURE_END0]))[row] = event.feature_end1;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_EDIT_START0]))[row] = event.interbase
		    ? event.insertion_boundary0 : (uint64_t)event.start1 - 1u;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_EDIT_END0]))[row] = event.interbase
		    ? event.insertion_boundary0 : (uint64_t)event.end1;
		if (event.interbase) {
			((uint64_t *)duckdb_vector_get_data(
			    fields[DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0]))[row] =
			    event.insertion_boundary0;
		} else {
			duckdb_validity_set_row_invalid(
			    duckdb_vector_get_validity(
			        fields[DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0]), row);
		}
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_REF_DIFF_OFFSET]))[row] =
		    event.ref_diff_offset;
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_REF_DIFF_LENGTH]))[row] =
		    event.ref_diff_length;
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_ALT_DIFF_OFFSET]))[row] =
		    event.alt_diff_offset;
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_ALT_DIFF_LENGTH]))[row] =
		    event.alt_diff_length;
	}
}

static void
duckvep_register_allele_geometry_scalar(duckdb_connection connection)
{
	duckdb_scalar_function scalar;
	duckdb_logical_type position_type, varchar_type, result_type;
	duckdb_logical_type field_types[DUCKVEP_GEOMETRY_FIELD_COUNT];
	size_t field;

	scalar = duckdb_create_scalar_function();
	position_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	field_types[DUCKVEP_GEOMETRY_KIND_CODE] =
	    duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	field_types[DUCKVEP_GEOMETRY_INTERBASE] =
	    duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	field_types[DUCKVEP_GEOMETRY_ANCHOR_SIDE_CODE] =
	    duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	for (field = DUCKVEP_GEOMETRY_RAW_START0;
	    field <= DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0; field++)
		field_types[field] =
		    duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	for (field = DUCKVEP_GEOMETRY_REF_DIFF_OFFSET;
	    field < DUCKVEP_GEOMETRY_FIELD_COUNT; field++)
		field_types[field] =
		    duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
	result_type = duckdb_create_struct_type(field_types,
	    duckvep_allele_geometry_field_names, DUCKVEP_GEOMETRY_FIELD_COUNT);
	duckdb_scalar_function_set_name(scalar, "duckvep_allele_geometry");
	duckdb_scalar_function_add_parameter(scalar, position_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_set_return_type(scalar, result_type);
	duckdb_scalar_function_set_special_handling(scalar);
	duckdb_scalar_function_set_function(scalar,
	    duckvep_allele_geometry_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);
	duckdb_destroy_logical_type(&position_type);
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&result_type);
	for (field = 0; field < DUCKVEP_GEOMETRY_FIELD_COUNT; field++)
		duckdb_destroy_logical_type(&field_types[field]);
}

static int
duckvep_scalar_dna_copy(uint8_t *destination, const char *source, size_t length)
{
	size_t index;

	for (index = 0; index < length; index++) {
		unsigned char base;

		base = (unsigned char)source[index] & UINT8_C(0xdf);
		if (base != 'A' && base != 'C' && base != 'G' && base != 'T' &&
		    base != 'N')
			return 0;
		destination[index] = base;
	}
	return 1;
}

static int
duckvep_scalar_ascii_equals(const duckdb_string_t *value, const char *text)
{
	const char *bytes;
	size_t index, length, text_length;

	if (value == NULL || text == NULL)
		return 0;
	bytes = duckdb_string_t_data((duckdb_string_t *)value);
	length = (size_t)duckdb_string_t_length(*value);
	text_length = strlen(text);
	if (length != text_length)
		return 0;
	for (index = 0; index < length; index++) {
		unsigned char left, right;

		left = (unsigned char)bytes[index];
		right = (unsigned char)text[index];
		if (left >= 'a' && left <= 'z')
			left = (unsigned char)(left - ('a' - 'A'));
		if (right >= 'a' && right <= 'z')
			right = (unsigned char)(right - ('a' - 'A'));
		if (left != right)
			return 0;
	}
	return 1;
}

static int
duckvep_scalar_sv_type(const duckdb_string_t *value, uint8_t *result)
{
	if (duckvep_scalar_ascii_equals(value, "DEL") ||
	    duckvep_scalar_ascii_equals(value, "DELETION"))
		*result = (uint8_t)DUCKVEP_SV_DELETION;
	else if (duckvep_scalar_ascii_equals(value, "DUP") ||
	    duckvep_scalar_ascii_equals(value, "DUPLICATION"))
		*result = (uint8_t)DUCKVEP_SV_DUPLICATION;
	else if (duckvep_scalar_ascii_equals(value, "TDUP") ||
	    duckvep_scalar_ascii_equals(value, "TANDEM_DUPLICATION"))
		*result = (uint8_t)DUCKVEP_SV_TANDEM_DUPLICATION;
	else if (duckvep_scalar_ascii_equals(value, "STR") ||
	    duckvep_scalar_ascii_equals(value, "TANDEM_REPEAT") ||
	    duckvep_scalar_ascii_equals(value, "CNV:TR"))
		*result = (uint8_t)DUCKVEP_SV_TANDEM_REPEAT;
	else if (duckvep_scalar_ascii_equals(value, "INV") ||
	    duckvep_scalar_ascii_equals(value, "INVERSION"))
		*result = (uint8_t)DUCKVEP_SV_INVERSION;
	else if (duckvep_scalar_ascii_equals(value, "INS") ||
	    duckvep_scalar_ascii_equals(value, "INSERTION"))
		*result = (uint8_t)DUCKVEP_SV_INSERTION;
	else if (duckvep_scalar_ascii_equals(value, "CNV") ||
	    duckvep_scalar_ascii_equals(value, "COPY_NUMBER_VARIATION"))
		*result = (uint8_t)DUCKVEP_SV_CNV;
	else if (duckvep_scalar_ascii_equals(value, "UNKNOWN"))
		*result = (uint8_t)DUCKVEP_SV_UNKNOWN;
	else if (duckvep_scalar_ascii_equals(value, "BND") ||
	    duckvep_scalar_ascii_equals(value, "BREAKEND"))
		*result = (uint8_t)DUCKVEP_SV_BREAKEND;
	else
		return 0;
	return 1;
}

static int
duckvep_scalar_copy_change(const duckdb_string_t *value, uint8_t *result)
{
	if (duckvep_scalar_ascii_equals(value, "UNKNOWN"))
		*result = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
	else if (duckvep_scalar_ascii_equals(value, "LOSS"))
		*result = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
	else if (duckvep_scalar_ascii_equals(value, "NEUTRAL"))
		*result = (uint8_t)DUCKVEP_COPY_CHANGE_NEUTRAL;
	else if (duckvep_scalar_ascii_equals(value, "GAIN"))
		*result = (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
	else
		return 0;
	return 1;
}

static int
duckvep_scalar_prepare_batch(duckvep_scalar_state_t *state,
	duckdb_data_chunk input, duckvep_variant_batch_t *batch,
	char *error, size_t error_size)
{
	duckdb_vector region_vector, position_vector;
	duckdb_vector reference_vector, alternate_vector;
	duckdb_string_t *references, *alternates;
	uint32_t *regions;
	uint64_t *positions;
	uint64_t *region_validity, *position_validity;
	uint64_t *reference_validity, *alternate_validity;
	size_t allele_size;
	idx_t row, rows;

	rows = duckdb_data_chunk_get_size(input);
	region_vector = duckdb_data_chunk_get_vector(input, 1);
	position_vector = duckdb_data_chunk_get_vector(input, 2);
	reference_vector = duckdb_data_chunk_get_vector(input, 3);
	alternate_vector = duckdb_data_chunk_get_vector(input, 4);
	if (!duckvep_scalar_variant_reserve(state, (size_t)rows)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory");
		return 0;
	}
	references = duckdb_vector_get_data(reference_vector);
	alternates = duckdb_vector_get_data(alternate_vector);
	regions = duckdb_vector_get_data(region_vector);
	positions = duckdb_vector_get_data(position_vector);
	region_validity = duckdb_vector_get_validity(region_vector);
	position_validity = duckdb_vector_get_validity(position_vector);
	reference_validity = duckdb_vector_get_validity(reference_vector);
	alternate_validity = duckdb_vector_get_validity(alternate_vector);
	allele_size = 0;
	for (row = 0; row < rows; row++) {
		size_t reference_length, alternate_length;
		uint32_t region;
		uint64_t position;

		if (duckvep_validity_is_null(region_validity, row) ||
		    duckvep_validity_is_null(position_validity, row) ||
		    duckvep_validity_is_null(reference_validity, row) ||
		    duckvep_validity_is_null(alternate_validity, row)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: variant fields cannot be NULL");
			return 0;
		}
		region = regions[row];
		position = positions[row];
		reference_length = (size_t)duckdb_string_t_length(references[row]);
		alternate_length = (size_t)duckdb_string_t_length(alternates[row]);
		if (reference_length == 0 || alternate_length == 0 ||
		    reference_length > UINT16_MAX || alternate_length > UINT16_MAX ||
		    region > UINT16_MAX || position == 0 || position > UINT32_MAX ||
		    reference_length - 1 > UINT32_MAX - position ||
		    reference_length > UINT32_MAX - allele_size ||
		    alternate_length > UINT32_MAX - allele_size - reference_length) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: variant exceeds the compact kernel limits");
			return 0;
		}
		state->seq_regions[row] = (uint16_t)region;
		state->positions[row] = (uint32_t)position;
		state->ends[row] = (uint32_t)(position + reference_length - 1);
		state->reference_offsets[row] = (uint32_t)allele_size;
		state->reference_lengths[row] = (uint16_t)reference_length;
		allele_size += reference_length;
		state->alternate_offsets[row] = (uint32_t)allele_size;
		state->alternate_lengths[row] = (uint16_t)alternate_length;
		allele_size += alternate_length;
	}
	if (!duckvep_scalar_allele_reserve(state, allele_size)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory");
		return 0;
	}
	for (row = 0; row < rows; row++) {
		int unspecified_alt =
		    duckvep_scalar_ascii_equals(&alternates[row], "<*>");

		if (!duckvep_scalar_dna_copy(state->allele_bytes +
		    state->reference_offsets[row],
		    duckdb_string_t_data(&references[row]),
		    state->reference_lengths[row])) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: REF/ALT must be unequal A/C/G/T/N alleles");
			return 0;
		}
		if (unspecified_alt) {
			memcpy(state->allele_bytes + state->alternate_offsets[row],
			    "<*>", 3u);
			state->variant_kinds[row] =
			    (uint8_t)DUCKVEP_KIND_UNSPECIFIED_ALT;
		} else if (!duckvep_scalar_dna_copy(state->allele_bytes +
		    state->alternate_offsets[row],
		    duckdb_string_t_data(&alternates[row]),
		    state->alternate_lengths[row]) ||
		    !duckvep_scalar_simple_kind(state->positions[row],
		    state->allele_bytes +
		    state->reference_offsets[row], state->reference_lengths[row],
		    state->allele_bytes + state->alternate_offsets[row],
		    state->alternate_lengths[row], &state->variant_kinds[row])) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: REF/ALT must be unequal A/C/G/T/N alleles or exact <*>");
			return 0;
		}
		state->sv_types[row] = DUCKVEP_SV_NONE;
		state->copy_changes[row] = DUCKVEP_COPY_CHANGE_UNKNOWN;
	}
	memset(batch, 0, sizeof(*batch));
	batch->chrom_id = state->seq_regions;
	batch->pos1 = state->positions;
	batch->end1 = state->ends;
	batch->ref_offset = state->reference_offsets;
	batch->ref_length = state->reference_lengths;
	batch->alt_offset = state->alternate_offsets;
	batch->alt_length = state->alternate_lengths;
	batch->allele_bytes = state->allele_bytes;
	batch->allele_bytes_len = allele_size;
	batch->variant_kind = state->variant_kinds;
	batch->sv_type = state->sv_types;
	batch->copy_change = state->copy_changes;
	batch->count = (size_t)rows;
	return 1;
}

static int
duckvep_scalar_prepare_sv_batch(duckvep_scalar_state_t *state,
	duckdb_data_chunk input, duckvep_variant_batch_t *batch,
	char *error, size_t error_size)
{
	duckdb_vector region_vector, start_vector, end_vector;
	duckdb_vector type_vector, copy_vector;
	duckdb_string_t *types, *copies;
	uint32_t *regions;
	uint64_t *starts, *ends;
	uint64_t *region_validity, *start_validity, *end_validity;
	uint64_t *type_validity, *copy_validity;
	idx_t row, rows;

	rows = duckdb_data_chunk_get_size(input);
	region_vector = duckdb_data_chunk_get_vector(input, 1);
	start_vector = duckdb_data_chunk_get_vector(input, 2);
	end_vector = duckdb_data_chunk_get_vector(input, 3);
	type_vector = duckdb_data_chunk_get_vector(input, 4);
	copy_vector = duckdb_data_chunk_get_vector(input, 5);
	if (!duckvep_scalar_variant_reserve(state, (size_t)rows)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory");
		return 0;
	}
	regions = duckdb_vector_get_data(region_vector);
	starts = duckdb_vector_get_data(start_vector);
	ends = duckdb_vector_get_data(end_vector);
	types = duckdb_vector_get_data(type_vector);
	copies = duckdb_vector_get_data(copy_vector);
	region_validity = duckdb_vector_get_validity(region_vector);
	start_validity = duckdb_vector_get_validity(start_vector);
	end_validity = duckdb_vector_get_validity(end_vector);
	type_validity = duckdb_vector_get_validity(type_vector);
	copy_validity = duckdb_vector_get_validity(copy_vector);
	for (row = 0; row < rows; row++) {
		uint8_t sv_type, copy_change;

		if (duckvep_validity_is_null(region_validity, row) ||
		    duckvep_validity_is_null(start_validity, row) ||
		    duckvep_validity_is_null(end_validity, row) ||
		    duckvep_validity_is_null(type_validity, row) ||
		    duckvep_validity_is_null(copy_validity, row)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: event fields cannot be NULL");
			return 0;
		}
		if (regions[row] > UINT16_MAX || starts[row] == 0u ||
		    starts[row] > UINT32_MAX || ends[row] < starts[row] ||
		    ends[row] > UINT32_MAX) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: invalid 1-based inclusive interval");
			return 0;
		}
		if (!duckvep_scalar_sv_type(&types[row], &sv_type) ||
		    !duckvep_scalar_copy_change(&copies[row], &copy_change)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: unknown structural type or copy change");
			return 0;
		}
		if (sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: provide both mate coordinates for a breakend");
			return 0;
		}
		if (!duckvep_sv_metadata_valid((duckvep_sv_type_t)sv_type,
		    (duckvep_copy_change_t)copy_change)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: structural type contradicts copy change");
			return 0;
		}
		if (!duckvep_sv_geometry_valid((duckvep_sv_type_t)sv_type,
		    (uint32_t)starts[row], (uint32_t)ends[row])) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: invalid structural event geometry");
			return 0;
		}
		state->seq_regions[row] = (uint16_t)regions[row];
		state->positions[row] = (uint32_t)starts[row];
		state->ends[row] = (uint32_t)ends[row];
		state->variant_kinds[row] = (uint8_t)DUCKVEP_KIND_SV;
		state->sv_types[row] = sv_type;
		state->copy_changes[row] = copy_change;
	}
	memset(batch, 0, sizeof(*batch));
	batch->chrom_id = state->seq_regions;
	batch->pos1 = state->positions;
	batch->end1 = state->ends;
	batch->variant_kind = state->variant_kinds;
	batch->sv_type = state->sv_types;
	batch->copy_change = state->copy_changes;
	batch->count = (size_t)rows;
	return 1;
}

static int
duckvep_scalar_prepare_breakend_batch(duckvep_scalar_state_t *state,
	duckdb_data_chunk input, duckvep_variant_batch_t *batch,
	char *error, size_t error_size)
{
	duckdb_vector region_vector, position_vector;
	duckdb_vector mate_region_vector, mate_position_vector;
	uint32_t *regions, *mate_regions;
	uint64_t *positions, *mate_positions;
	uint64_t *region_validity, *position_validity;
	uint64_t *mate_region_validity, *mate_position_validity;
	idx_t row, rows;

	rows = duckdb_data_chunk_get_size(input);
	region_vector = duckdb_data_chunk_get_vector(input, 1);
	position_vector = duckdb_data_chunk_get_vector(input, 2);
	mate_region_vector = duckdb_data_chunk_get_vector(input, 3);
	mate_position_vector = duckdb_data_chunk_get_vector(input, 4);
	if (!duckvep_scalar_variant_reserve(state, (size_t)rows)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory");
		return 0;
	}
	regions = duckdb_vector_get_data(region_vector);
	positions = duckdb_vector_get_data(position_vector);
	mate_regions = duckdb_vector_get_data(mate_region_vector);
	mate_positions = duckdb_vector_get_data(mate_position_vector);
	region_validity = duckdb_vector_get_validity(region_vector);
	position_validity = duckdb_vector_get_validity(position_vector);
	mate_region_validity = duckdb_vector_get_validity(mate_region_vector);
	mate_position_validity = duckdb_vector_get_validity(mate_position_vector);
	for (row = 0; row < rows; row++) {
		if (duckvep_validity_is_null(region_validity, row) ||
		    duckvep_validity_is_null(position_validity, row) ||
		    duckvep_validity_is_null(mate_region_validity, row) ||
		    duckvep_validity_is_null(mate_position_validity, row)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: both loci are required");
			return 0;
		}
		if (regions[row] > UINT16_MAX || mate_regions[row] > UINT16_MAX ||
		    positions[row] == 0u || positions[row] >= UINT32_MAX ||
		    mate_positions[row] == 0u || mate_positions[row] > UINT32_MAX) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: invalid one-based endpoint coordinate");
			return 0;
		}
		state->seq_regions[row] = (uint16_t)regions[row];
		state->positions[row] = (uint32_t)positions[row];
		state->ends[row] = (uint32_t)positions[row];
		state->mate_seq_regions[row] = (uint16_t)mate_regions[row];
		state->mate_positions[row] = (uint32_t)mate_positions[row];
		state->variant_kinds[row] = (uint8_t)DUCKVEP_KIND_SV;
		state->sv_types[row] = (uint8_t)DUCKVEP_SV_BREAKEND;
		state->copy_changes[row] = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
	}
	memset(batch, 0, sizeof(*batch));
	batch->chrom_id = state->seq_regions;
	batch->pos1 = state->positions;
	batch->end1 = state->ends;
	batch->mate_chrom_id = state->mate_seq_regions;
	batch->mate_pos1 = state->mate_positions;
	batch->variant_kind = state->variant_kinds;
	batch->sv_type = state->sv_types;
	batch->copy_change = state->copy_changes;
	batch->count = (size_t)rows;
	return 1;
}

static int
duckvep_scalar_model_region(const duckvep_owned_model_t *model,
	uint16_t seq_region, uint32_t *sequence_length)
{
	size_t begin, end;

	begin = 0;
	end = model->known_seq_region_count;
	while (begin < end) {
		size_t middle;

		middle = begin + (end - begin) / 2;
		if (model->known_seq_regions[middle] < seq_region)
			begin = middle + 1;
		else
			end = middle;
	}
	if (begin >= model->known_seq_region_count ||
	    model->known_seq_regions[begin] != seq_region)
		return 0;
	if (sequence_length != NULL)
		*sequence_length = model->sequence_lengths[begin];
	return 1;
}

static int
duckvep_scalar_model_region_index(const duckvep_owned_model_t *model,
	uint16_t seq_region, size_t *index_out)
{
	size_t begin, end;

	if (model == NULL || index_out == NULL)
		return 0;
	begin = 0;
	end = model->known_seq_region_count;
	while (begin < end) {
		size_t middle;

		middle = begin + (end - begin) / 2;
		if (model->known_seq_regions[middle] < seq_region)
			begin = middle + 1;
		else
			end = middle;
	}
	if (begin >= model->known_seq_region_count ||
	    model->known_seq_regions[begin] != seq_region)
		return 0;
	*index_out = begin;
	return 1;
}

static int
duckvep_scalar_reference_windows(duckvep_scalar_state_t *state,
	const duckvep_event_t *event, int *available,
	duckvep_hgvs_reference_window_t *shift_window,
	duckvep_hgvs_reference_window_t *lookup_window,
	char *error, size_t error_size)
{
	duckvep_owned_model_t *model;
	duckvep_workspace_cache_t *cache;
	size_t region_index;
	uint32_t shift_start1, shift_end1;
	uint32_t fetch_start1, fetch_end1;
	uint32_t cache_fetch_end1;
	uint32_t sequence_length;
	uint64_t cached_end1;
	duckvep_hgvs_status_t hgvs_status;
	hts_pos_t fetched_length;
	char *bases;

	if (available != NULL)
		*available = 0;
	if (shift_window != NULL)
		memset(shift_window, 0, sizeof(*shift_window));
	if (lookup_window != NULL)
		memset(lookup_window, 0, sizeof(*lookup_window));
	if (state == NULL || state->entry == NULL ||
	    state->workspace_cache == NULL || event == NULL ||
	    available == NULL || shift_window == NULL || lookup_window == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: invalid reference-window state");
		return 0;
	}
	model = &state->entry->model;
	cache = state->workspace_cache;
	if (model->reference_fasta_path == NULL)
		return 1;
	if (!duckvep_scalar_model_region_index(model, event->chrom_id,
	    &region_index) || model->sequence_names == NULL ||
	    model->sequence_names[region_index] == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: sequence-region name is absent from the reference-enabled model");
		return 0;
	}
	sequence_length = model->sequence_lengths[region_index];
	hgvs_status = duckvep_hgvs_genomic_search_interval(event,
	    sequence_length, &shift_start1, &shift_end1);
	if (hgvs_status != DUCKVEP_HGVS_OK) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: semantic edit has no valid VEP genomic shift interval");
		return 0;
	}
	hgvs_status = duckvep_hgvs_reference_fetch_interval(event,
	    sequence_length, &fetch_start1, &fetch_end1);
	if (hgvs_status != DUCKVEP_HGVS_OK ||
	    shift_start1 < fetch_start1 || shift_end1 > fetch_end1) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: semantic edit has no valid bounded reference lookup interval");
		return 0;
	}
	if (cache->reference_fai == NULL) {
		if (!duckvep_model_reference_identity_matches(model)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: model reference FASTA or index changed after model load");
			return 0;
		}
		cache->reference_fai = fai_load3_format(
		    model->reference_fasta_open_path,
		    model->reference_fai_open_path,
		    model->reference_gzi_open_path, 0, FAI_FASTA);
		if (cache->reference_fai == NULL) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: could not open the model reference FASTA/index");
			return 0;
		}
		if (!duckvep_model_reference_identity_matches(model)) {
			fai_destroy(cache->reference_fai);
			cache->reference_fai = NULL;
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: model reference FASTA or index changed while a worker was opening it");
			return 0;
		}
	}
	cached_end1 = cache->reference_length != 0u ?
	    (uint64_t)cache->reference_start1 +
	    (uint64_t)cache->reference_length - 1u : 0u;
	if (cache->reference_bases == NULL ||
	    cache->reference_chrom_id != event->chrom_id ||
	    fetch_start1 < cache->reference_start1 ||
	    (uint64_t)fetch_end1 > cached_end1) {
		if (!duckvep_model_reference_identity_matches(model)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: pinned reference FASTA or index changed after model load");
			return 0;
		}
		cache_fetch_end1 = fetch_end1;
		if ((uint64_t)cache_fetch_end1 + DUCKVEP_REFERENCE_READ_AHEAD <
		    (uint64_t)sequence_length) {
			cache_fetch_end1 += DUCKVEP_REFERENCE_READ_AHEAD;
		} else {
			cache_fetch_end1 = sequence_length;
		}
		fetched_length = -1;
		bases = faidx_fetch_seq64(cache->reference_fai,
		    model->sequence_names[region_index], (hts_pos_t)fetch_start1 - 1,
		    (hts_pos_t)cache_fetch_end1 - 1, &fetched_length);
		if (bases == NULL || fetched_length < 0 ||
		    (uint64_t)fetched_length !=
		    (uint64_t)cache_fetch_end1 - (uint64_t)fetch_start1 + 1u) {
			free(bases);
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: reference FASTA fetch did not return the requested interval");
			return 0;
		}
		if (!duckvep_model_reference_identity_matches(model)) {
			free(bases);
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: pinned reference FASTA or index changed during fetch");
			return 0;
		}
		free(cache->reference_bases);
		cache->reference_bases = bases;
		cache->reference_length = (size_t)fetched_length;
		cache->reference_start1 = fetch_start1;
		cache->reference_chrom_id = event->chrom_id;
	}
	lookup_window->bases = (const uint8_t *)cache->reference_bases;
	lookup_window->length = cache->reference_length;
	lookup_window->start1 = cache->reference_start1;
	lookup_window->chrom_id = cache->reference_chrom_id;
	shift_window->bases = (const uint8_t *)cache->reference_bases +
	    (size_t)(shift_start1 - cache->reference_start1);
	shift_window->length =
	    (size_t)((uint64_t)shift_end1 - (uint64_t)shift_start1 + 1u);
	shift_window->start1 = shift_start1;
	shift_window->chrom_id = cache->reference_chrom_id;
	*available = 1;
	return 1;
}

static uint8_t
duckvep_scalar_hgvs_reason(duckvep_hgvs_status_t status)
{
	switch (status) {
	case DUCKVEP_HGVS_OK:
		return DUCKVEP_HGVS_ADAPTER_NONE;
	case DUCKVEP_HGVS_MISSING_REFERENCE:
		return DUCKVEP_HGVS_ADAPTER_MISSING_REFERENCE;
	case DUCKVEP_HGVS_REFERENCE_MISMATCH:
		return DUCKVEP_HGVS_ADAPTER_REFERENCE_MISMATCH;
	case DUCKVEP_HGVS_INVALID_ALLELE:
		return DUCKVEP_HGVS_ADAPTER_INVALID_ALLELE;
	case DUCKVEP_HGVS_UNSUPPORTED_EDIT:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT;
	case DUCKVEP_HGVS_BUFFER_TOO_SMALL:
		return DUCKVEP_HGVS_ADAPTER_INTERNAL_CAPACITY;
	case DUCKVEP_HGVS_MISSING_PEPTIDE:
		return DUCKVEP_HGVS_ADAPTER_MISSING_SEQUENCE;
	case DUCKVEP_HGVS_NOT_APPLICABLE:
		return DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
	case DUCKVEP_HGVS_UNSUPPORTED_PROTEIN:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_PROTEIN;
	case DUCKVEP_HGVS_INVALID_ARG:
	case DUCKVEP_HGVS_INVALID_PROJECTION:
	case DUCKVEP_HGVS_OUT_OF_RANGE:
	default:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	}
}

static void
duckvep_scalar_hgvs_set_transcript_failure(
	const duckvep_owned_model_t *model, const duckvep_event_t *event,
	uint32_t tx_idx, uint8_t reason, duckvep_hgvs_scalar_result_t *result)
{
	if (result == NULL)
		return;
	result->transcript_reason = reason;
	if (model != NULL && event != NULL &&
	    duckvep_project_feature_has_coding_precondition_unshifted(
	    &model->transcripts, &model->exons, tx_idx, event))
		result->protein_reason = reason;
}

static int
duckvep_scalar_hgvs_event_wholly_outside_transcript(
	const duckvep_event_t *event,
	const duckvep_transcript_model_t *transcripts, uint32_t tx_idx)
{
	uint32_t first_pos1, last_pos1;

	if (event == NULL || transcripts == NULL ||
	    tx_idx >= transcripts->transcript_count ||
	    transcripts->start1 == NULL || transcripts->end1 == NULL)
		return 0;
	if (event->interbase != 0u) {
		first_pos1 = event->insertion_boundary0;
		last_pos1 = duckvep_event_right_flank1(event);
	} else {
		first_pos1 = event->start1;
		last_pos1 = event->end1;
	}
	return last_pos1 < transcripts->start1[tx_idx] ||
	    first_pos1 > transcripts->end1[tx_idx];
}

static uint8_t
duckvep_scalar_hgvs_edit_reason(duckvep_transcript_edit_status_t status,
	const duckvep_event_t *event,
	const duckvep_transcript_model_t *transcripts, uint32_t tx_idx)
{
	switch (status) {
	case DUCKVEP_TRANSCRIPT_EDIT_OK:
		return DUCKVEP_HGVS_ADAPTER_NONE;
	case DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT:
		if (event != NULL && event->interbase != 0u &&
		    transcripts != NULL && tx_idx < transcripts->transcript_count &&
		    transcripts->start1 != NULL && transcripts->end1 != NULL &&
		    (event->insertion_boundary0 < transcripts->start1[tx_idx] ||
		     duckvep_event_right_flank1(event) >
		         transcripts->end1[tx_idx]))
			return DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
		return duckvep_scalar_hgvs_event_wholly_outside_transcript(
		    event, transcripts, tx_idx) ?
		    DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE :
		    DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	case DUCKVEP_TRANSCRIPT_EDIT_UNSUPPORTED_KIND:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT;
	case DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG:
	case DUCKVEP_TRANSCRIPT_EDIT_INVALID_EVENT:
	case DUCKVEP_TRANSCRIPT_EDIT_INVALID_PROJECTION:
	case DUCKVEP_TRANSCRIPT_EDIT_OUT_OF_RANGE:
	default:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	}
}

static uint8_t
duckvep_scalar_hgvs_context_reason(
	duckvep_variant_coding_context_status_t status)
{
	switch (status) {
	case DUCKVEP_VARIANT_CODING_CONTEXT_OK:
		return DUCKVEP_HGVS_ADAPTER_NONE;
	case DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_CDS:
	case DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE:
	case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG:
	case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_EVENT:
	case DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_ORDER:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	case DUCKVEP_VARIANT_CODING_CONTEXT_UNSUPPORTED_KIND:
	case DUCKVEP_VARIANT_CODING_CONTEXT_NON_CONTIGUOUS:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT;
	case DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL:
	case DUCKVEP_VARIANT_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL:
	case DUCKVEP_VARIANT_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL:
	case DUCKVEP_VARIANT_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL:
		return DUCKVEP_HGVS_ADAPTER_INTERNAL_CAPACITY;
	case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ALLELE:
	case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_BASE:
		return DUCKVEP_HGVS_ADAPTER_INVALID_ALLELE;
	case DUCKVEP_VARIANT_CODING_CONTEXT_REF_MISMATCH:
		return DUCKVEP_HGVS_ADAPTER_REFERENCE_MISMATCH;
	default:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	}
}

static const char *
duckvep_scalar_hgvs_reason_name(uint8_t reason)
{
	switch ((duckvep_hgvs_adapter_reason_t)reason) {
	case DUCKVEP_HGVS_ADAPTER_NONE:
		return NULL;
	case DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE:
		return "not_applicable";
	case DUCKVEP_HGVS_ADAPTER_MISSING_REFERENCE:
		return "missing_reference";
	case DUCKVEP_HGVS_ADAPTER_MISSING_SEQUENCE:
		return "missing_sequence";
	case DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION:
		return "invalid_projection";
	case DUCKVEP_HGVS_ADAPTER_REFERENCE_MISMATCH:
		return "reference_mismatch";
	case DUCKVEP_HGVS_ADAPTER_INVALID_ALLELE:
		return "invalid_allele";
	case DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT:
		return "unsupported_edit";
	case DUCKVEP_HGVS_ADAPTER_INTERNAL_CAPACITY:
		return "internal_capacity";
	case DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_TAIL:
		return "missing_transcript_tail";
	case DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_FLANK:
		return "missing_transcript_flank";
	case DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_PROTEIN:
		return "unsupported_protein";
	default:
		return "invalid_status";
	}
}

static int
duckvep_scalar_hgvs_store_dna(duckvep_scalar_state_t *state,
	const duckvep_hgvs_dna_fact_t *fact,
	duckvep_hgvs_scalar_result_t *result, char *error, size_t error_size)
{
	duckvep_hgvs_status_t status;
	size_t required;

	if (state->hgvs_render_capacity == 0u &&
	    !duckvep_scalar_hgvs_render_reserve(state,
	    DUCKVEP_HGVS_INITIAL_RENDER_CAPACITY)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory rendering transcript HGVS");
		return 0;
	}
	required = 0u;
	status = duckvep_hgvs_dna_render_basic(fact,
	    state->hgvs_render_scratch, state->hgvs_render_capacity, &required);
	if (status == DUCKVEP_HGVS_BUFFER_TOO_SMALL) {
		if (required == SIZE_MAX ||
		    !duckvep_scalar_hgvs_render_reserve(state, required + 1u)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: out of memory rendering transcript HGVS");
			return 0;
		}
		status = duckvep_hgvs_dna_render_basic(fact,
		    state->hgvs_render_scratch, state->hgvs_render_capacity,
		    &required);
	}
	if (status != DUCKVEP_HGVS_OK) {
		result->transcript_reason = duckvep_scalar_hgvs_reason(status);
		return 1;
	}
	if (!duckvep_scalar_hgvs_text_append(state,
	    state->hgvs_render_scratch, required, &result->transcript_offset)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: transcript HGVS output exceeds the vector limit");
		return 0;
	}
	result->transcript_length = (uint32_t)required;
	result->shift_offset = (uint32_t)fact->shift_offset;
	result->transcript_reason = DUCKVEP_HGVS_ADAPTER_NONE;
	return 1;
}

static int
duckvep_scalar_hgvs_store_protein(duckvep_scalar_state_t *state,
	const duckvep_hgvs_protein_fact_t *fact,
	duckvep_hgvs_scalar_result_t *result, char *error, size_t error_size)
{
	duckvep_hgvs_status_t status;
	size_t required;

	if (state->hgvs_render_capacity == 0u &&
	    !duckvep_scalar_hgvs_render_reserve(state,
	    DUCKVEP_HGVS_INITIAL_RENDER_CAPACITY)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory rendering protein HGVS");
		return 0;
	}
	required = 0u;
	status = duckvep_hgvs_protein_render(fact, 0,
	    state->hgvs_render_scratch, state->hgvs_render_capacity, &required);
	if (status == DUCKVEP_HGVS_BUFFER_TOO_SMALL) {
		if (required == SIZE_MAX ||
		    !duckvep_scalar_hgvs_render_reserve(state, required + 1u)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: out of memory rendering protein HGVS");
			return 0;
		}
		status = duckvep_hgvs_protein_render(fact, 0,
		    state->hgvs_render_scratch, state->hgvs_render_capacity,
		    &required);
	}
	if (status != DUCKVEP_HGVS_OK) {
		result->protein_reason = duckvep_scalar_hgvs_reason(status);
		return 1;
	}
	if (!duckvep_scalar_hgvs_text_append(state,
	    state->hgvs_render_scratch, required, &result->protein_offset)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: protein HGVS output exceeds the vector limit");
		return 0;
	}
	result->protein_length = (uint32_t)required;
	result->protein_reason = DUCKVEP_HGVS_ADAPTER_NONE;
	return 1;
}

static uint8_t
duckvep_scalar_hgvs_delta_reason(duckvep_context_delta_status_t status)
{
	switch (status) {
	case DUCKVEP_CONTEXT_DELTA_OK:
		return DUCKVEP_HGVS_ADAPTER_NONE;
	case DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL:
		return DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_TAIL;
	case DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK:
		return DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_FLANK;
	case DUCKVEP_CONTEXT_DELTA_UNSUPPORTED:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_PROTEIN;
	case DUCKVEP_CONTEXT_DELTA_INVALID_ARG:
	default:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	}
}

static uint8_t
duckvep_scalar_hgvs_cds_reason(duckvep_cds_edit_status_t status)
{
	switch (status) {
	case DUCKVEP_CDS_EDIT_OK:
		return DUCKVEP_HGVS_ADAPTER_NONE;
	case DUCKVEP_CDS_EDIT_OUT_OF_CDS:
		return DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
	case DUCKVEP_CDS_EDIT_REF_MISMATCH:
		return DUCKVEP_HGVS_ADAPTER_REFERENCE_MISMATCH;
	case DUCKVEP_CDS_EDIT_INVALID_ALLELE:
		return DUCKVEP_HGVS_ADAPTER_INVALID_ALLELE;
	case DUCKVEP_CDS_EDIT_NON_CONTIGUOUS:
	case DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT;
	case DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL:
		return DUCKVEP_HGVS_ADAPTER_INTERNAL_CAPACITY;
	case DUCKVEP_CDS_EDIT_INVALID_ARG:
	case DUCKVEP_CDS_EDIT_INVALID_EVENT:
	default:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	}
}

static uint8_t
duckvep_scalar_hgvs_sequence_reason(uint8_t status)
{
	switch ((duckvep_sequence_status_t)status) {
	case DUCKVEP_SEQUENCE_MISSING:
		return DUCKVEP_HGVS_ADAPTER_MISSING_SEQUENCE;
	case DUCKVEP_SEQUENCE_AMBIGUOUS:
		return DUCKVEP_HGVS_ADAPTER_INVALID_ALLELE;
	case DUCKVEP_SEQUENCE_REFERENCE_MISMATCH:
		return DUCKVEP_HGVS_ADAPTER_REFERENCE_MISMATCH;
	case DUCKVEP_SEQUENCE_NON_CONTIGUOUS_EDIT:
	case DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_EDIT;
	case DUCKVEP_SEQUENCE_INTERNAL_CAPACITY:
		return DUCKVEP_HGVS_ADAPTER_INTERNAL_CAPACITY;
	case DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL:
		return DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_TAIL;
	case DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_FLANK:
		return DUCKVEP_HGVS_ADAPTER_MISSING_TRANSCRIPT_FLANK;
	case DUCKVEP_SEQUENCE_INVALID_PROJECTION:
		return DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	case DUCKVEP_SEQUENCE_NOT_APPLICABLE:
		return DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
	case DUCKVEP_SEQUENCE_RESOLVED:
	default:
		return DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_PROTEIN;
	}
}

/* Build and render one transcript row while the consequence-pass facts still
 * borrows the worker scratch. */
static int
duckvep_scalar_build_hgvs_pair(duckvep_scalar_state_t *state,
	const duckvep_variant_batch_t *batch, size_t variant,
	const duckvep_event_t *event,
	const duckvep_hgvs_reference_window_t *shift_reference,
	const duckvep_hgvs_reference_window_t *lookup_reference,
	duckvep_hgvs_status_t uploaded_reference_status,
	const duckvep_consequence_t *consequence,
	const duckvep_pair_facts_t *facts,
	duckvep_hgvs_scalar_result_t *hgvs_result,
	char *error, size_t error_size)
{
	duckvep_owned_model_t *model;
	duckvep_delta_scratch_t *scratch;
	duckvep_transcript_edit_t transcript_edit;
	duckvep_hgvs_dna_fact_t dna_fact;
	duckvep_hgvs_status_t hgvs_status;
	duckvep_transcript_edit_status_t edit_status;
	duckvep_coding_context_t context;
	const duckvep_coding_context_t *context_ptr;
	duckvep_sequence_delta_t delta;
	duckvep_hgvs_protein_fact_t protein_fact;
	duckvep_variant_coding_context_status_t context_status;
	duckvep_context_delta_status_t delta_status;
	duckvep_feature_substitution_result_t feature_result;
	duckvep_event_t shifted_event;
	duckvep_haplotype_edit_t shifted_edit;
	duckvep_coding_context_t late_context;
	uint32_t tx_idx;
	int protein_coordinates_defined;
	int consequence_predicates_valid;
	int delta_ready;
	int shifted_edit_ready;

	if (state == NULL || state->entry == NULL || state->workspace == NULL ||
	    batch == NULL || event == NULL || consequence == NULL ||
	    hgvs_result == NULL || variant >= batch->count) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: invalid transcript row");
		return 0;
	}
	if (consequence->overlap_object_kind !=
	    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT)
		return 1;
	model = &state->entry->model;
	scratch = duckvep_workspace_delta_scratch(state->workspace);
	if (scratch == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: worker scratch is unavailable");
		return 0;
	}
	tx_idx = consequence->tx_idx;
	consequence_predicates_valid =
	    (consequence->flags & (uint32_t)
	        DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID) != 0u;
	hgvs_result->transcript_reason =
	    DUCKVEP_HGVS_ADAPTER_INVALID_PROJECTION;
	hgvs_result->protein_reason = DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
	if (uploaded_reference_status != DUCKVEP_HGVS_OK) {
		duckvep_scalar_hgvs_set_transcript_failure(
		    model, event, tx_idx,
		    duckvep_scalar_hgvs_reason(uploaded_reference_status),
		    hgvs_result);
		return 1;
	}
	edit_status = facts != NULL
	    ? (duckvep_transcript_edit_status_t)facts->transcript_edit_status
	    : DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
	if (edit_status != DUCKVEP_TRANSCRIPT_EDIT_OK) {
		duckvep_scalar_hgvs_set_transcript_failure(
		    model, event, tx_idx,
		    duckvep_scalar_hgvs_edit_reason(edit_status, event,
		    &model->transcripts, tx_idx), hgvs_result);
		return 1;
	}
	if (facts->transcript_edit == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: successful pair facts omit transcript edit");
		return 0;
	}
	transcript_edit = *facts->transcript_edit;
	if (lookup_reference == NULL)
		(void)duckvep_transcript_edit_cds_fill_prepared(
		    &model->transcripts, &model->exons, &model->sequences,
		    batch, scratch->edits, scratch->edits_cap, &transcript_edit);
	if ((duckvep_cds_edit_status_t)transcript_edit.cds_status ==
	        DUCKVEP_CDS_EDIT_REF_MISMATCH ||
	    (duckvep_cds_edit_status_t)transcript_edit.cds_status ==
	        DUCKVEP_CDS_EDIT_INVALID_ALLELE) {
		hgvs_result->transcript_reason = duckvep_scalar_hgvs_cds_reason(
		    (duckvep_cds_edit_status_t)transcript_edit.cds_status);
		hgvs_result->protein_reason = hgvs_result->transcript_reason;
		return 1;
	}
	if (lookup_reference == NULL &&
	    (event->ref_diff_offset != 0u ||
	     event->ref_diff_length != batch->ref_length[variant])) {
		duckvep_scalar_hgvs_set_transcript_failure(
		    model, event, tx_idx, DUCKVEP_HGVS_ADAPTER_MISSING_REFERENCE,
		    hgvs_result);
		return 1;
	}
	if (lookup_reference == NULL &&
	    (event->kind == (uint8_t)DUCKVEP_KIND_INS ||
	     event->kind == (uint8_t)DUCKVEP_KIND_DEL ||
	     (duckvep_cds_edit_status_t)transcript_edit.cds_status !=
	         DUCKVEP_CDS_EDIT_OK)) {
		duckvep_scalar_hgvs_set_transcript_failure(
		    model, event, tx_idx, DUCKVEP_HGVS_ADAPTER_MISSING_REFERENCE,
		    hgvs_result);
		return 1;
	}
	hgvs_status = duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
	    &model->transcripts, &model->exons, shift_reference,
	    lookup_reference, &transcript_edit, &dna_fact);
	if (hgvs_status != DUCKVEP_HGVS_OK) {
		duckvep_scalar_hgvs_set_transcript_failure(
		    model, event, tx_idx, duckvep_scalar_hgvs_reason(hgvs_status),
		    hgvs_result);
		return 1;
	}
	if (!duckvep_scalar_hgvs_store_dna(state, &dna_fact, hgvs_result,
	    error, error_size))
		return 0;
	if (hgvs_result->transcript_reason != DUCKVEP_HGVS_ADAPTER_NONE)
		return 1;

	hgvs_status = duckvep_hgvs_protein_coordinates_defined(
	    &model->transcripts, &model->exons, &transcript_edit, &dna_fact,
	    &protein_coordinates_defined);
	if (hgvs_status != DUCKVEP_HGVS_OK) {
		hgvs_result->protein_reason = duckvep_scalar_hgvs_reason(hgvs_status);
		return 1;
	}
	if (!protein_coordinates_defined) {
		hgvs_result->protein_reason = DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
		return 1;
	}
	if (model->sequences.cds_length == NULL ||
	    tx_idx >= model->sequences.transcript_count ||
	    model->sequences.cds_length[tx_idx] == 0u) {
		hgvs_result->protein_reason = DUCKVEP_HGVS_ADAPTER_MISSING_SEQUENCE;
		return 1;
	}

	if (consequence_predicates_valid && consequence->protein_pos > 0 &&
	    consequence->aa_ref != 0u && consequence->aa_alt != 0u &&
	    transcript_edit.feature_ref_length ==
	        transcript_edit.feature_alt_length &&
	    (event->kind == (uint8_t)DUCKVEP_KIND_SNV ||
	     event->kind == (uint8_t)DUCKVEP_KIND_MNV)) {
		hgvs_status = duckvep_hgvs_protein_fact_build_single_residue(
		    (uint32_t)consequence->protein_pos, consequence->aa_ref,
		    consequence->aa_alt, consequence->flags,
		    DUCKVEP_COMPAT_VEP_116, &protein_fact);
		if (hgvs_status == DUCKVEP_HGVS_OK) {
			return duckvep_scalar_hgvs_store_protein(
			    state, &protein_fact, hgvs_result, error, error_size);
		}
		if (hgvs_status != DUCKVEP_HGVS_NOT_APPLICABLE) {
			hgvs_result->protein_reason =
			    duckvep_scalar_hgvs_reason(hgvs_status);
			return 1;
		}
	}

	context_status = DUCKVEP_VARIANT_CODING_CONTEXT_UNSUPPORTED_KIND;
	context_ptr = &context;
	delta_ready = 0;
	shifted_edit_ready = 0;
	memset(&shifted_event, 0, sizeof shifted_event);
	memset(&shifted_edit, 0, sizeof shifted_edit);
	/* The consequence callback is synchronous: these pointers are valid until
	 * this function returns.  An unshifted DNA fact describes the same edit, so
	 * reuse both context and delta instead of rebuilding projection, alternate
	 * CDS, and local peptides. */
	if (facts != NULL && facts->coding_context_valid &&
	    facts->coding_context != NULL && dna_fact.shift_offset == 0) {
		context_ptr = facts->coding_context;
		context_status = DUCKVEP_VARIANT_CODING_CONTEXT_OK;
		if (facts->delta != NULL) {
			delta = *facts->delta;
			delta_ready = delta.valid != 0u;
		}
	} else if ((dna_fact.ref_length == 0u) !=
	    (dna_fact.alt_length == 0u)) {
		duckvep_edit_set_t edit_set;
		size_t allele_required;

		allele_required = dna_fact.ref_length != 0u ?
		    (size_t)dna_fact.ref_length : (size_t)dna_fact.alt_length + 1u;
		if (!duckvep_scalar_hgvs_allele_reserve(state, allele_required)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: out of memory rotating a shifted allele");
			return 0;
		}
		hgvs_status = duckvep_hgvs_shifted_cds_edit_build(
		    &model->transcripts, &model->exons, &model->sequences,
		    lookup_reference, &transcript_edit, &dna_fact,
		    state->hgvs_allele_scratch, state->hgvs_allele_capacity,
		    &allele_required, &shifted_event, &shifted_edit);
		if (hgvs_status != DUCKVEP_HGVS_OK) {
			hgvs_result->protein_reason =
			    duckvep_scalar_hgvs_reason(hgvs_status);
			return 1;
		}
		shifted_edit_ready = 1;
		if (shifted_edit.ref_len == 0u && shifted_edit.cds_start >
		    model->sequences.cds_length[tx_idx]) {
			hgvs_result->protein_reason =
			    DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
			return 1;
		}
		edit_set.edits = &shifted_edit;
		edit_set.count = 1u;
		context_status = duckvep_model_coding_context_build(
		    &model->transcripts, &model->exons, &model->sequences, tx_idx,
		    model->transcripts.strand[tx_idx], &shifted_event, &edit_set,
		    scratch->alt_cds, scratch->alt_cds_cap,
		    scratch->ref_peptide, scratch->ref_peptide_cap,
		    scratch->alt_peptide, scratch->alt_peptide_cap, &context);
		context_ptr = &context;
	} else {
		feature_result = duckvep_feature_substitution_context_fill(
		    &model->transcripts, &model->exons, &model->sequences, batch,
		    (uint32_t)variant, tx_idx, model->transcripts.strand[tx_idx],
		    scratch, event,
		    facts != NULL ? facts->projection_exon_hint : UINT32_MAX,
		    &context, &delta);
		if (feature_result == DUCKVEP_FEATURE_SUBSTITUTION_DELTA_ONLY) {
			hgvs_result->protein_reason =
			    duckvep_scalar_hgvs_sequence_reason(delta.sequence_status);
			return 1;
		}
		if (feature_result == DUCKVEP_FEATURE_SUBSTITUTION_CONTEXT_READY) {
			context_status = DUCKVEP_VARIANT_CODING_CONTEXT_OK;
			context_ptr = &context;
			delta_ready = 1;
		} else {
			if (!transcript_edit.cds_built)
				(void)duckvep_transcript_edit_cds_fill_prepared(
				    &model->transcripts, &model->exons,
				    &model->sequences, batch, scratch->edits,
				    scratch->edits_cap, &transcript_edit);
			if ((duckvep_cds_edit_status_t)transcript_edit.cds_status !=
			    DUCKVEP_CDS_EDIT_OK) {
				hgvs_result->protein_reason =
				    duckvep_scalar_hgvs_cds_reason(
				    (duckvep_cds_edit_status_t)transcript_edit.cds_status);
				return 1;
			}
			context_status = duckvep_model_coding_context_build(
			    &model->transcripts, &model->exons, &model->sequences,
			    tx_idx, model->transcripts.strand[tx_idx], event,
			    &transcript_edit.cds_edits, scratch->alt_cds,
			    scratch->alt_cds_cap, scratch->ref_peptide,
			    scratch->ref_peptide_cap, scratch->alt_peptide,
			    scratch->alt_peptide_cap, &context);
			context_ptr = &context;
		}
	}
	if (context_status != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
		hgvs_result->protein_reason =
		    duckvep_scalar_hgvs_context_reason(context_status);
		return 1;
	}
	if (!delta_ready &&
	    duckvep_sequence_delta_consequence_flags_complete_for_hgvs(
	        context_ptr, consequence->flags)) {
		memset(&delta, 0, sizeof delta);
		delta.valid = 1u;
		duckvep_sequence_delta_apply_consequence_flags(
		    consequence->flags, &delta);
		delta_ready = 1;
	}
	if (!delta_ready) {
		delta_status = duckvep_coding_context_delta_fill(
		    context_ptr, model->transcript_flags[tx_idx], &delta);
		if (delta_status != DUCKVEP_CONTEXT_DELTA_OK) {
			hgvs_result->protein_reason =
			    duckvep_scalar_hgvs_delta_reason(delta_status);
			return 1;
		}
	}
	if (delta.valid == 0u) {
		hgvs_result->protein_reason =
		    DUCKVEP_HGVS_ADAPTER_UNSUPPORTED_PROTEIN;
		return 1;
	}
	if (consequence_predicates_valid)
		duckvep_sequence_delta_apply_consequence_flags(
		    consequence->flags, &delta);
	if (!consequence_predicates_valid &&
	    !duckvep_project_feature_overlaps_start_codon_unshifted(
	        &model->transcripts, &model->exons, tx_idx,
	        &transcript_edit.event))
		delta.start_lost = 0u;
	hgvs_status = duckvep_hgvs_protein_fact_build(
	    context_ptr, &delta, &protein_fact);
	if (hgvs_status != DUCKVEP_HGVS_OK) {
		hgvs_result->protein_reason = duckvep_scalar_hgvs_reason(hgvs_status);
		return 1;
	}
	if (protein_fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT &&
	    dna_fact.shift_offset != 0 &&
	    (consequence->region_mask & DUCKVEP_REGION_CDS) == 0u) {
		duckvep_haplotype_edit_t late_edit;
		duckvep_edit_set_t late_set;

		/* hgvs_protein deletes the shift hash before _stop_loss_extra_AA.
		 * Replay the restored original allele at its original CDS coordinate;
		 * most such rows retain Ter?, but the rare positive stop is real. */
		protein_fact.termination_known = 0u;
		protein_fact.termination_distance = 0u;
		if (shifted_edit_ready && dna_fact.shift_offset > 0 &&
		    shifted_edit.cds_start > (uint32_t)dna_fact.shift_offset) {
			late_edit = shifted_edit;
			late_edit.cds_start -= (uint32_t)dna_fact.shift_offset;
			/* Deleting shift_hash also restores the original feature alleles;
			 * the late stop search therefore does not consume the rotated HGVS
			 * display allele used by the first peptide pass. */
			late_edit.ref = transcript_edit.ref;
			late_edit.alt = transcript_edit.alt;
			late_edit.ref_len = transcript_edit.ref_length;
			late_edit.alt_len = transcript_edit.alt_length;
			late_set.edits = &late_edit;
			late_set.count = 1u;
			context_status = duckvep_model_coding_context_build(
			    &model->transcripts, &model->exons, &model->sequences,
			    tx_idx, model->transcripts.strand[tx_idx], event,
			    &late_set, scratch->alt_cds, scratch->alt_cds_cap,
			    scratch->ref_peptide, scratch->ref_peptide_cap,
			    scratch->alt_peptide, scratch->alt_peptide_cap,
			    &late_context);
			if (context_status == DUCKVEP_VARIANT_CODING_CONTEXT_OK)
				(void)duckvep_hgvs_protein_frameshift_termination_replay(
				    &late_context, &protein_fact);
		}
	}
	return duckvep_scalar_hgvs_store_protein(
	    state, &protein_fact, hgvs_result, error, error_size);
}

typedef struct duckvep_scalar_hgvs_observer {
	duckvep_scalar_state_t *state;
	size_t next_result;
	size_t prepared_variant;
	duckvep_hgvs_reference_window_t shift_reference;
	duckvep_hgvs_reference_window_t lookup_reference;
	duckvep_hgvs_status_t uploaded_reference_status;
	int reference_available;
	int variant_prepared;
	char *error;
	size_t error_size;
} duckvep_scalar_hgvs_observer_t;

static int
duckvep_scalar_hgvs_observe(void *observer_context,
	const duckvep_variant_batch_t *batch,
	const duckvep_consequence_t *consequence,
	const duckvep_pair_facts_t *facts)
{
	duckvep_scalar_hgvs_observer_t *observer;
	duckvep_scalar_state_t *state;
	duckvep_hgvs_scalar_result_t *result;
	size_t variant;

	observer = (duckvep_scalar_hgvs_observer_t *)observer_context;
	if (observer == NULL || observer->state == NULL || batch == NULL ||
	    consequence == NULL || consequence->variant_idx >= batch->count) {
		if (observer != NULL)
			duckvep_sql_set_error(observer->error, observer->error_size,
			    "duckvep_annotate: invalid fused annotation row");
		return 0;
	}
	state = observer->state;
	if (!duckvep_scalar_hgvs_result_reserve(
	    state, observer->next_result + 1u)) {
		duckvep_sql_set_error(observer->error, observer->error_size,
		    "duckvep_annotate: out of memory growing HGVS facts");
		return 0;
	}
	result = &state->hgvs_results[observer->next_result++];
	memset(result, 0, sizeof(*result));
	if (consequence->overlap_object_kind !=
	    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT)
		return 1;
	if (facts == NULL || facts->event == NULL) {
		duckvep_sql_set_error(observer->error, observer->error_size,
		    "duckvep_annotate: transcript row is missing its live pair facts");
		return 0;
	}
	/* Directional transcript candidates are deliberately emitted by the
	 * consequence sweep, but VEP has no transcript coordinate or protein
	 * projection for an edit wholly outside the transcript span. Avoid opening
	 * a FASTA window and running the general projector for this common 5 kb
	 * admission case. */
	if (consequence->region_mask == (uint32_t)DUCKVEP_REGION_UPSTREAM ||
	    consequence->region_mask == (uint32_t)DUCKVEP_REGION_DOWNSTREAM) {
		result->transcript_reason = DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
		result->protein_reason = DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE;
		return 1;
	}
	variant = (size_t)consequence->variant_idx;
	if (!observer->variant_prepared || observer->prepared_variant != variant) {
		duckvep_owned_model_t *model = &state->entry->model;

		observer->prepared_variant = variant;
		observer->variant_prepared = 1;
		observer->reference_available = 0;
		observer->uploaded_reference_status = DUCKVEP_HGVS_OK;
		memset(&observer->shift_reference, 0,
		    sizeof(observer->shift_reference));
		memset(&observer->lookup_reference, 0,
		    sizeof(observer->lookup_reference));
		if (model->reference_fasta_path != NULL &&
		    !duckvep_scalar_reference_windows(
		        state, facts->event, &observer->reference_available,
		        &observer->shift_reference, &observer->lookup_reference,
		        observer->error, observer->error_size))
			return 0;
		if (observer->reference_available) {
			observer->uploaded_reference_status =
			    duckvep_hgvs_uploaded_reference_validate(
		        &observer->lookup_reference, facts->event,
			        batch->allele_bytes + batch->ref_offset[variant],
			        (size_t)batch->ref_length[variant]);
		}
	}
	return duckvep_scalar_build_hgvs_pair(
	    state, batch, variant, facts->event,
	    observer->reference_available ? &observer->shift_reference : NULL,
	    observer->reference_available ? &observer->lookup_reference : NULL,
	    observer->uploaded_reference_status, consequence, facts, result,
	    observer->error, observer->error_size);
}

static int
duckvep_scalar_u32_compare(const void *left, const void *right)
{
	uint32_t a, b;

	a = *(const uint32_t *)left;
	b = *(const uint32_t *)right;
	return a < b ? -1 : a > b;
}

/* DuckDB list materialization consumes one contiguous run per input variant.
 * Transcript and regulation-feature pair evaluators each preserve variant
 * order, but appending those two streams does not. Keep the adapter result
 * contract explicit after both evaluators have run. Object kind and ordinal
 * make the within-variant order deterministic without changing semantics. */
static int
duckvep_scalar_bnd_result_compare(const void *left, const void *right)
{
	const duckvep_consequence_t *a, *b;
	uint32_t a_object, b_object;

	a = (const duckvep_consequence_t *)left;
	b = (const duckvep_consequence_t *)right;
	if (a->variant_idx != b->variant_idx)
		return a->variant_idx < b->variant_idx ? -1 : 1;
	if (a->overlap_object_kind != b->overlap_object_kind)
		return a->overlap_object_kind < b->overlap_object_kind ? -1 : 1;
	a_object = a->overlap_object_kind ==
	    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT ?
	    a->tx_idx : a->interval_feature_idx;
	b_object = b->overlap_object_kind ==
	    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT ?
	    b->tx_idx : b->interval_feature_idx;
	if (a_object != b_object)
		return a_object < b_object ? -1 : 1;
	if (a->consequence_mask != b->consequence_mask)
		return a->consequence_mask < b->consequence_mask ? -1 : 1;
	return 0;
}

/* The transcript and interval-feature evaluators each emit a sorted stream.
 * Only the smaller interval-feature stream needs sorting by object kind before
 * the streams are merged. Keeping that stream in worker-owned scratch avoids
 * sorting the already ordered transcript expansion, which is normally several
 * orders of magnitude larger. */
static int
duckvep_scalar_merge_bnd_results(duckvep_scalar_state_t *state,
	size_t begin, size_t middle, size_t end, char *error, size_t error_size)
{
	size_t left, right, write, feature_count;

	feature_count = end - middle;
	if (feature_count == 0u)
		return 1;
	if (feature_count > 1u)
		qsort(state->results + middle, feature_count,
		    sizeof(*state->results), duckvep_scalar_bnd_result_compare);
	if (middle == begin)
		return 1;
	if (!duckvep_scalar_result_merge_reserve(state, feature_count)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory merging result streams");
		return 0;
	}
	memcpy(state->result_merge, state->results + middle,
	    feature_count * sizeof(*state->result_merge));
	left = middle;
	right = feature_count;
	write = end;
	while (left > begin && right > 0u) {
		if (duckvep_scalar_bnd_result_compare(&state->results[left - 1u],
		    &state->result_merge[right - 1u]) > 0)
			state->results[--write] = state->results[--left];
		else
			state->results[--write] = state->result_merge[--right];
	}
	while (right > 0u)
		state->results[--write] = state->result_merge[--right];
	return 1;
}

static int
duckvep_scalar_seed_index(duckvep_scalar_state_t *state,
	cgranges_t *index, int index_complete, uint16_t seq_region,
	uint32_t position, uint32_t halo_distance,
	duckvep_annotate_cursor_t *cursor, int interval_features,
	char *error, size_t error_size)
{
	duckvep_error_t kernel_error;
	int32_t query_start, query_end;
	int64_t hit_count, hit;
	char region_name[16];

	if (!index_complete ||
	    position > (uint32_t)INT32_MAX ||
	    halo_distance >= (uint32_t)INT32_MAX ||
	    position > (uint32_t)INT32_MAX - halo_distance)
		return 1;
	query_start = position > halo_distance + 1 ?
	    (int32_t)(position - halo_distance - 1) : 0;
	query_end = (int32_t)(position + halo_distance);
	(void)snprintf(region_name, sizeof(region_name), "%u",
	    (unsigned)seq_region);
	hit_count = cr_overlap(index, region_name, query_start,
	    query_end, &state->interval_hits, &state->interval_hit_capacity);
	if (hit_count < 0 ||
	    !duckvep_scalar_seed_reserve(state, (size_t)hit_count)) {
		duckvep_sql_set_error(error, error_size, interval_features ?
		    "duckvep_annotate: could not seed the sorted regulation-feature sweep" :
		    "duckvep_annotate: could not seed the sorted transcript sweep");
		return 0;
	}
	for (hit = 0; hit < hit_count; hit++)
		state->seed_transcripts[hit] = (uint32_t)cr_label(
		    index, state->interval_hits[hit]);
	qsort(state->seed_transcripts, (size_t)hit_count,
	    sizeof(*state->seed_transcripts), duckvep_scalar_u32_compare);
	memset(&kernel_error, 0, sizeof(kernel_error));
	if ((interval_features ?
	    duckvep_annotate_cursor_seed_interval_features(cursor,
	    state->seed_transcripts, (size_t)hit_count, &kernel_error) :
	    duckvep_annotate_cursor_seed(cursor, state->seed_transcripts,
	    (size_t)hit_count, &kernel_error)) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "duckvep_annotate: %s",
		    kernel_error.message);
		return 0;
	}
	return 1;
}

static int
duckvep_scalar_seed_cursor(duckvep_scalar_state_t *state,
	uint16_t seq_region, uint32_t position, uint32_t halo_distance,
	duckvep_annotate_cursor_t *cursor, char *error, size_t error_size)
{
	duckvep_owned_model_t *model;

	model = &state->entry->model;
	if (!duckvep_scalar_seed_index(state, model->interval_index,
	    model->interval_index_complete, seq_region, position, halo_distance,
	    cursor, 0, error, error_size))
		return 0;
	if (model->interval_feature_count == 0u)
		return 1;
	return duckvep_scalar_seed_index(state, model->interval_feature_index,
	    model->interval_feature_index_complete, seq_region, position, 0u,
	    cursor, 1, error, error_size);
}

static duckvep_variant_batch_t
duckvep_scalar_batch_slice(const duckvep_variant_batch_t *batch,
	size_t begin, size_t count)
{
	duckvep_variant_batch_t slice;

	slice = *batch;
	slice.chrom_id += begin;
	slice.pos1 += begin;
	slice.end1 += begin;
	if (slice.mate_chrom_id != NULL)
		slice.mate_chrom_id += begin;
	if (slice.mate_pos1 != NULL)
		slice.mate_pos1 += begin;
	if (slice.ref_offset != NULL)
		slice.ref_offset += begin;
	if (slice.ref_length != NULL)
		slice.ref_length += begin;
	if (slice.alt_offset != NULL)
		slice.alt_offset += begin;
	if (slice.alt_length != NULL)
		slice.alt_length += begin;
	slice.variant_kind += begin;
	slice.sv_type += begin;
	slice.copy_change += begin;
	slice.count = count;
	return slice;
}

static int
duckvep_scalar_run(duckvep_scalar_state_t *state,
	const duckvep_variant_batch_t *batch, size_t begin, size_t count,
	uint64_t upstream_distance, uint64_t downstream_distance,
	int with_hgvs, char *error, size_t error_size)
{
	duckvep_variant_batch_t slice;
	duckvep_annotate_cursor_t *cursor;
	duckvep_options_init_t options_init;
	duckvep_error_t kernel_error;
	duckvep_status_t status;
	uint32_t halo_distance;
	uint32_t sequence_length;
	size_t variant;
	duckvep_scalar_hgvs_observer_t hgvs_observer;

	if (upstream_distance > UINT32_MAX || downstream_distance > UINT32_MAX) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: upstream/downstream distance exceeds the uint32 kernel limit");
		return 0;
	}
	halo_distance = upstream_distance > downstream_distance ?
	    (uint32_t)upstream_distance : (uint32_t)downstream_distance;
	if (!duckvep_scalar_model_region(&state->entry->model,
	    batch->chrom_id[begin], &sequence_length)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: seq_region is absent from the loaded model");
		return 0;
	}
	for (variant = begin; variant < begin + count; variant++) {
		state->transcript_coverage_complete[variant] = (uint8_t)
		    state->entry->model.transcript_coverage_complete;
		if (sequence_length != 0 && batch->end1[variant] > sequence_length) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: variant span exceeds sequence-region length");
			return 0;
		}
	}
	if (state->workspace == NULL) {
		state->workspace_cache = duckvep_registry_workspace_take(
		    state->registry, state->entry, error, error_size);
		if (state->workspace_cache == NULL)
			return 0;
		state->workspace = state->workspace_cache->workspace;
	}
	if (!state->have_options_distances ||
	    state->options_upstream_distance != (uint32_t)upstream_distance ||
	    state->options_downstream_distance != (uint32_t)downstream_distance) {
		duckvep_options_close(state->options);
		state->options = NULL;
		state->have_options_distances = 0;
		memset(&options_init, 0, sizeof(options_init));
		options_init.upstream_dist = (uint32_t)upstream_distance;
		options_init.downstream_dist = (uint32_t)downstream_distance;
		options_init.halo = halo_distance;
		options_init.distances_are_explicit = 1u;
		options_init.compatibility_profile =
		    (uint8_t)DUCKVEP_COMPAT_VEP_116;
		memset(&kernel_error, 0, sizeof(kernel_error));
		if (duckvep_options_open(&options_init, &state->options,
		    &kernel_error) != DUCKVEP_OK) {
			(void)snprintf(error, error_size, "duckvep_annotate: %s",
			    kernel_error.message);
			return 0;
		}
		state->options_upstream_distance = (uint32_t)upstream_distance;
		state->options_downstream_distance = (uint32_t)downstream_distance;
		state->have_options_distances = 1;
	}
	slice = duckvep_scalar_batch_slice(batch, begin, count);
	cursor = NULL;
	if (duckvep_annotate_cursor_open(state->entry->model.kernel, &slice,
	    state->options, state->workspace, &cursor,
	    &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "duckvep_annotate: %s",
		    kernel_error.message);
		return 0;
	}
	if (!duckvep_scalar_seed_cursor(state, batch->chrom_id[begin],
	    batch->pos1[begin], halo_distance, cursor, error, error_size)) {
		duckvep_annotate_cursor_close(cursor);
		return 0;
	}
	memset(&hgvs_observer, 0, sizeof(hgvs_observer));
	if (with_hgvs) {
		hgvs_observer.state = state;
		hgvs_observer.next_result = state->result_count;
		hgvs_observer.error = error;
		hgvs_observer.error_size = error_size;
		duckvep_annotate_cursor_set_observer(
		    cursor, duckvep_scalar_hgvs_observe, &hgvs_observer);
	}
	while (!duckvep_annotate_cursor_done(cursor)) {
		duckvep_result_builder_t builder;
		size_t old_count, index;

		old_count = state->result_count;
		if (!duckvep_scalar_result_reserve(state,
		    old_count + (count > 64 ? count : 64))) {
			duckvep_annotate_cursor_close(cursor);
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: out of memory growing results");
			return 0;
		}
		duckvep_result_builder_init(&builder, state->results + old_count,
		    state->result_capacity - old_count);
		memset(&kernel_error, 0, sizeof(kernel_error));
		status = duckvep_annotate_cursor_fill(cursor, &builder,
		    &kernel_error);
		for (index = 0; index < builder.count; index++) {
			duckvep_consequence_t *row;

			/* The builder already aliases the final result-array suffix.
			 * Adjust adapter-owned ordinals in place instead of copying the
			 * complete hot result record onto the same address. */
			row = &builder.rows[index];
			row->variant_idx += (uint32_t)begin;
			if (row->overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT)
				row->gene_idx =
				    state->entry->model.gene_indices[row->tx_idx];
		}
		state->result_count = old_count + builder.count;
		if (with_hgvs &&
		    hgvs_observer.next_result != state->result_count) {
			duckvep_annotate_cursor_close(cursor);
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: fused result streams diverged");
			return 0;
		}
		if (status == DUCKVEP_ERR_RESULT_FULL)
			continue;
		if (status != DUCKVEP_OK) {
			duckvep_annotate_cursor_close(cursor);
			if (error == NULL || error[0] == '\0')
				(void)snprintf(error, error_size, "duckvep_annotate: %s",
				    kernel_error.message);
			return 0;
		}
	}
	duckvep_annotate_cursor_close(cursor);
	return 1;
}

static int
duckvep_scalar_append_endpoint_candidates(duckvep_scalar_state_t *state,
	cgranges_t *index, int index_complete, uint16_t seq_region,
	uint32_t position, uint32_t halo_distance, int interval_features,
	size_t *candidate_count, char *error, size_t error_size)
{
	int32_t query_start, query_end;
	int64_t hit_count, hit;
	char region_name[16];

	if (!index_complete || position > (uint32_t)INT32_MAX ||
	    halo_distance > (uint32_t)INT32_MAX ||
	    position > (uint32_t)INT32_MAX - halo_distance) {
		duckvep_sql_set_error(error, error_size,
		    interval_features ?
		    "duckvep_annotate: resident regulation-feature index cannot represent an endpoint" :
		    "duckvep_annotate: resident transcript index cannot represent an endpoint");
		return 0;
	}
	query_start = position > halo_distance + 1u ?
	    (int32_t)(position - halo_distance - 1u) : 0;
	query_end = (int32_t)(position + halo_distance);
	(void)snprintf(region_name, sizeof(region_name), "%u",
	    (unsigned)seq_region);
	hit_count = cr_overlap(index, region_name, query_start,
	    query_end, &state->interval_hits, &state->interval_hit_capacity);
	if (hit_count < 0 || (uint64_t)hit_count > SIZE_MAX - *candidate_count ||
	    !duckvep_scalar_seed_reserve(state,
	    *candidate_count + (size_t)hit_count)) {
		duckvep_sql_set_error(error, error_size,
		    interval_features ?
		    "duckvep_annotate: could not collect endpoint regulation features" :
		    "duckvep_annotate: could not collect endpoint transcripts");
		return 0;
	}
	for (hit = 0; hit < hit_count; hit++)
		state->seed_transcripts[(*candidate_count)++] =
		    (uint32_t)cr_label(index,
		    state->interval_hits[hit]);
	return 1;
}

/* Candidate discovery remains an adapter concern because cgranges indexes the
 * resident DuckDB-owned model. Exact transcript or interval-feature semantics
 * remain in the kernel. Both object kinds share ordering, deduplication,
 * allocation, and result-index rebasing here so BND support cannot drift into
 * two subtly different join implementations. */
static int
duckvep_scalar_run_breakend_object_pairs(duckvep_scalar_state_t *state,
	const duckvep_variant_batch_t *batch, size_t begin, size_t count,
	cgranges_t *index, int index_complete, uint32_t search_distance,
	int interval_features, char *error, size_t error_size)
{
	duckvep_variant_batch_t slice;
	duckvep_candidate_pairs_t transcript_pairs;
	duckvep_interval_feature_pairs_t feature_pairs;
	duckvep_error_t kernel_error;
	duckvep_result_builder_t builder;
	duckvep_status_t status;
	size_t pair_count, variant, old_count, result;

	pair_count = 0u;
	for (variant = 0u; variant < count; variant++) {
		size_t candidate_count = 0u;
		size_t candidate, unique_count;
		uint32_t local_point = batch->pos1[begin + variant] + 1u;

		if (!duckvep_scalar_append_endpoint_candidates(state, index,
		    index_complete, batch->chrom_id[begin + variant], local_point,
		    search_distance, interval_features, &candidate_count, error,
		    error_size) ||
		    !duckvep_scalar_append_endpoint_candidates(state, index,
		    index_complete, batch->mate_chrom_id[begin + variant],
		    batch->mate_pos1[begin + variant], search_distance,
		    interval_features, &candidate_count, error, error_size))
			return 0;
		if (candidate_count > 1u)
			qsort(state->seed_transcripts, candidate_count,
			    sizeof(*state->seed_transcripts),
			    duckvep_scalar_u32_compare);
		unique_count = 0u;
		for (candidate = 0u; candidate < candidate_count; candidate++) {
			if (candidate == 0u || state->seed_transcripts[candidate] !=
			    state->seed_transcripts[candidate - 1u])
				state->seed_transcripts[unique_count++] =
				    state->seed_transcripts[candidate];
		}
		if (unique_count > SIZE_MAX - pair_count ||
		    !duckvep_scalar_pair_reserve(state, pair_count + unique_count)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: out of memory collecting candidate pairs");
			return 0;
		}
		for (candidate = 0u; candidate < unique_count; candidate++) {
			state->pair_variant_indices[pair_count] = (uint32_t)variant;
			state->pair_object_indices[pair_count] =
			    state->seed_transcripts[candidate];
			pair_count++;
		}
	}
	if (pair_count == 0u)
		return 1;
	old_count = state->result_count;
	if (pair_count > SIZE_MAX - old_count ||
	    !duckvep_scalar_result_reserve(state, old_count + pair_count)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: out of memory growing results");
		return 0;
	}
	slice = duckvep_scalar_batch_slice(batch, begin, count);
	duckvep_result_builder_init(&builder, state->results + old_count,
	    state->result_capacity - old_count);
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (interval_features) {
		feature_pairs.variant_idx = state->pair_variant_indices;
		feature_pairs.feature_idx = state->pair_object_indices;
		feature_pairs.count = pair_count;
		status = duckvep_annotate_interval_feature_pairs(
		    state->entry->model.kernel, &slice, &feature_pairs,
		    state->options, state->workspace, &builder, &kernel_error);
	} else {
		transcript_pairs.variant_idx = state->pair_variant_indices;
		transcript_pairs.tx_idx = state->pair_object_indices;
		transcript_pairs.count = pair_count;
		status = duckvep_annotate_pairs(state->entry->model.kernel, &slice,
		    &transcript_pairs, state->options, state->workspace, &builder,
		    &kernel_error);
	}
	if (status != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "duckvep_annotate: %s",
		    kernel_error.message);
		return 0;
	}
	for (result = 0u; result < builder.count; result++) {
		duckvep_consequence_t *row = &state->results[old_count + result];
		row->variant_idx += (uint32_t)begin;
		if (row->overlap_object_kind ==
		    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT)
			row->gene_idx =
			    state->entry->model.gene_indices[row->tx_idx];
	}
	state->result_count = old_count + builder.count;
	return 1;
}

static int
duckvep_scalar_run_breakends(duckvep_scalar_state_t *state,
	const duckvep_variant_batch_t *batch, size_t begin, size_t count,
	uint64_t upstream_distance, uint64_t downstream_distance,
	char *error, size_t error_size)
{
	duckvep_options_init_t options_init;
	duckvep_error_t kernel_error;
	uint32_t halo_distance;
	uint32_t search_distance;
	size_t result_begin, transcript_end, variant;

	if (upstream_distance > UINT32_MAX || downstream_distance > UINT32_MAX) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: transcript distance exceeds the uint32 kernel limit");
		return 0;
	}
	halo_distance = upstream_distance > downstream_distance ?
	    (uint32_t)upstream_distance : (uint32_t)downstream_distance;
	search_distance = halo_distance < DUCKVEP_BREAKEND_ALLELE_DISTANCE ?
	    halo_distance : DUCKVEP_BREAKEND_ALLELE_DISTANCE;
	for (variant = begin; variant < begin + count; variant++) {
		uint32_t local_length, mate_length;

		state->transcript_coverage_complete[variant] = (uint8_t)
		    state->entry->model.transcript_coverage_complete;
		if (!duckvep_scalar_model_region(&state->entry->model,
		    batch->chrom_id[variant], &local_length) ||
		    !duckvep_scalar_model_region(&state->entry->model,
		    batch->mate_chrom_id[variant], &mate_length)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: endpoint region is absent from the loaded model");
			return 0;
		}
		if ((local_length != 0u && batch->pos1[variant] > local_length) ||
		    (mate_length != 0u && batch->mate_pos1[variant] > mate_length)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: endpoint exceeds sequence-region length");
			return 0;
		}
	}
	if (state->workspace == NULL) {
		state->workspace_cache = duckvep_registry_workspace_take(
		    state->registry, state->entry, error, error_size);
		if (state->workspace_cache == NULL)
			return 0;
		state->workspace = state->workspace_cache->workspace;
	}
	if (!state->have_options_distances ||
	    state->options_upstream_distance != (uint32_t)upstream_distance ||
	    state->options_downstream_distance != (uint32_t)downstream_distance) {
		duckvep_options_close(state->options);
		state->options = NULL;
		state->have_options_distances = 0;
		memset(&options_init, 0, sizeof(options_init));
		options_init.upstream_dist = (uint32_t)upstream_distance;
		options_init.downstream_dist = (uint32_t)downstream_distance;
		options_init.halo = halo_distance;
		options_init.distances_are_explicit = 1u;
		options_init.compatibility_profile =
		    (uint8_t)DUCKVEP_COMPAT_VEP_116;
		memset(&kernel_error, 0, sizeof(kernel_error));
		if (duckvep_options_open(&options_init, &state->options,
		    &kernel_error) != DUCKVEP_OK) {
			(void)snprintf(error, error_size,
			    "duckvep_annotate: %s", kernel_error.message);
			return 0;
		}
		state->options_upstream_distance = (uint32_t)upstream_distance;
		state->options_downstream_distance = (uint32_t)downstream_distance;
		state->have_options_distances = 1;
	}

	result_begin = state->result_count;
	if (!duckvep_scalar_run_breakend_object_pairs(state, batch, begin, count,
	    state->entry->model.interval_index,
	    state->entry->model.interval_index_complete, search_distance, 0,
	    error, error_size))
		return 0;
	transcript_end = state->result_count;
	if (state->entry->model.interval_feature_count != 0u &&
	    !duckvep_scalar_run_breakend_object_pairs(state, batch, begin, count,
	    state->entry->model.interval_feature_index,
	    state->entry->model.interval_feature_index_complete, 0u, 1,
	    error, error_size))
		return 0;
	return duckvep_scalar_merge_bnd_results(state, result_begin,
	    transcript_end, state->result_count, error, error_size);
}

static void
duckvep_scalar_set_null(duckdb_vector vector, uint64_t **validity, idx_t row)
{
	if (*validity == NULL) {
		duckdb_vector_ensure_validity_writable(vector);
		*validity = duckdb_vector_get_validity(vector);
	}
	(*validity)[row / 64] &= ~(UINT64_C(1) << (row % 64));
}

static void
duckvep_scalar_set_null_range(duckdb_vector vector, uint64_t **validity,
	size_t count)
{
	size_t words;

	if (count == 0u)
		return;
	if (*validity == NULL) {
		duckdb_vector_ensure_validity_writable(vector);
		*validity = duckdb_vector_get_validity(vector);
	}
	words = (count + 63u) / 64u;
	memset(*validity, 0, words * sizeof(**validity));
}

static void
duckvep_scalar_set_valid(uint64_t *validity, idx_t row)
{
	validity[row / 64] |= UINT64_C(1) << (row % 64);
}

static void
duckvep_scalar_assign_ascii(duckdb_vector vector, idx_t row,
	const char *text, size_t length)
{
	duckdb_vector_assign_string_element_len(vector, row, text,
	    (idx_t)length);
}

static void
duckvep_scalar_assign_ascii_valid(duckdb_vector vector, uint64_t *validity,
	idx_t row, const char *text, size_t length)
{
	/* The checked DuckDB assignment marks invalid UTF-8 NULL.  Establish the
	 * initially valid state before assigning so a rejection remains NULL;
	 * restoring validity afterwards would expose an uninitialized string_t.
	 * This ordering also works with the DuckDB 1.2 stable extension API, which
	 * predates the explicit unsafe string-assignment entry point. */
	duckvep_scalar_set_valid(validity, row);
	duckvep_scalar_assign_ascii(vector, row, text, length);
}

static int
duckvep_scalar_format_terms(uint64_t terms, char *buffer, size_t size,
	size_t *length)
{
	*length = size == 0 ? 0 : duckvep_so_render(terms, '&', buffer, size);
	return size != 0 && *length < size;
}

static int
duckvep_scalar_format_region(uint32_t region, char *buffer, size_t size,
	size_t *result_length)
{
	static const struct {
		uint32_t bit;
		const char *name;
	} names[] = {
		{DUCKVEP_REGION_UPSTREAM, "upstream"},
		{DUCKVEP_REGION_DOWNSTREAM, "downstream"},
		{DUCKVEP_REGION_INTRON, "intron"},
		{DUCKVEP_REGION_EXON, "non_coding_exon"},
		{DUCKVEP_REGION_CDS, "CDS"},
		{DUCKVEP_REGION_UTR, "UTR"},
		{DUCKVEP_REGION_SPLICE, "splice"}
	};
	size_t index, length;

	if (size == 0)
		return 0;
	length = 0;
	buffer[0] = '\0';
	for (index = 0; index < sizeof(names) / sizeof(names[0]); index++) {
		size_t name_length;

		if ((region & names[index].bit) == 0)
			continue;
		name_length = strlen(names[index].name);
		if ((length != 0 && length + 1 >= size) ||
		    name_length >= size - length)
			return 0;
		if (length != 0)
			buffer[length++] = '&';
		memcpy(buffer + length, names[index].name, name_length);
		length += name_length;
		buffer[length] = '\0';
	}
	*result_length = length;
	return length != 0;
}

static const char *
duckvep_scalar_sequence_reason(uint8_t status)
{
	switch ((duckvep_sequence_status_t)status) {
	case DUCKVEP_SEQUENCE_MISSING:
		return "missing_sequence";
	case DUCKVEP_SEQUENCE_AMBIGUOUS:
		return "ambiguous_sequence";
	case DUCKVEP_SEQUENCE_REFERENCE_MISMATCH:
		return "reference_mismatch";
	case DUCKVEP_SEQUENCE_NON_CONTIGUOUS_EDIT:
		return "non_contiguous_cds_edit";
	case DUCKVEP_SEQUENCE_INVALID_PROJECTION:
		return "invalid_model_projection";
	case DUCKVEP_SEQUENCE_INTERNAL_CAPACITY:
		return "internal_capacity_error";
	case DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL:
		return "missing_transcript_tail";
	case DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_FLANK:
		return "missing_transcript_flank";
	case DUCKVEP_SEQUENCE_NOT_APPLICABLE:
	case DUCKVEP_SEQUENCE_RESOLVED:
	case DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT:
	default:
		return "unsupported_compound_consequence";
	}
}

static const char *
duckvep_scalar_overlap_object_name(uint8_t kind)
{
	switch ((duckvep_overlap_object_kind_t)kind) {
	case DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT:
		return "transcript";
	case DUCKVEP_OVERLAP_OBJECT_REGULATORY_REGION:
		return "regulatory_region";
	case DUCKVEP_OVERLAP_OBJECT_TF_BINDING_SITE:
		return "transcription_factor_binding_site";
	default:
		return NULL;
	}
}

static int
duckvep_result_range_has_transcript(const duckvep_scalar_state_t *state,
	size_t begin, size_t end)
{
	/* Both scalar paths preserve the adapter's within-variant object order:
	 * transcript rows precede regulatory regions and motif features.  The BND
	 * merge comparator makes the same enum order explicit.  Therefore a range
	 * contains a transcript exactly when its first row is a transcript; do not
	 * rescan every expanded row while sizing and writing the DuckDB list. */
	return begin < end &&
	    state->results[begin].overlap_object_kind ==
	    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
}

static int
duckvep_scalar_prepare_output_list(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, size_t *prepared_output_count,
	char *error, size_t error_size)
{
	duckdb_list_entry *lists;
	size_t output_capacity, output_count, source, row;

	/* At most one synthetic no-transcript row is added per input. Reserve that
	 * tight upper bound, then discover exact list offsets and final child size
	 * in one ordered pass instead of scanning every expanded transcript row
	 * once to count and again to assign offsets. */
	if (state->result_count > SIZE_MAX - (size_t)input_rows) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: list result size overflow");
		return 0;
	}
	output_capacity = state->result_count + (size_t)input_rows;
	if (duckdb_list_vector_reserve(output, (idx_t)output_capacity) ==
	    DuckDBError) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: could not allocate list result");
		return 0;
	}
	lists = duckdb_vector_get_data(output);
	output_count = 0;
	source = 0;
	for (row = 0; row < (size_t)input_rows; row++) {
		size_t begin;

		begin = source;
		while (source < state->result_count &&
		    state->results[source].variant_idx == (uint32_t)row)
			source++;
		lists[row].offset = output_count;
		lists[row].length = source - begin;
		if (!duckvep_result_range_has_transcript(state, begin, source))
			lists[row].length++;
		output_count += lists[row].length;
	}
	if (source != state->result_count) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: result rows are out of order");
		return 0;
	}
	if (duckdb_list_vector_set_size(output, (idx_t)output_count) ==
	    DuckDBError) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: could not size list result");
		return 0;
	}
	*prepared_output_count = output_count;
	return 1;
}

/* Compact SQL codes preserve the kernel enum values for sequence failures.
 * Zero means no failure; one is reserved for the adapter-only case where a
 * partial model has no loaded transcript feature for the input. */
enum {
	DUCKVEP_COMPACT_STATUS_SUPPORTED = 0,
	DUCKVEP_COMPACT_STATUS_UNRESOLVED = 1,
	DUCKVEP_COMPACT_REASON_NONE = 0,
	DUCKVEP_COMPACT_REASON_NO_FEATURE = 1
};

static int
duckvep_scalar_write_hgvs_fields(duckvep_scalar_state_t *state,
	duckdb_vector *vectors, uint64_t **validity, size_t first_column,
	size_t output_row, size_t result_index, char *error, size_t error_size)
{
	const duckvep_hgvs_scalar_result_t *hgvs;
	const char *reason;
	uint32_t *hgvs_shifts;

	hgvs = &state->hgvs_results[result_index];
	if ((hgvs->transcript_length != 0u &&
	    (hgvs->transcript_offset > state->hgvs_text_size ||
	    hgvs->transcript_length > state->hgvs_text_size -
	    hgvs->transcript_offset)) ||
	    (hgvs->protein_length != 0u &&
	    (hgvs->protein_offset > state->hgvs_text_size ||
	    hgvs->protein_length > state->hgvs_text_size -
	    hgvs->protein_offset))) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: rendered HGVS slice exceeds the vector text arena");
		return 0;
	}
	hgvs_shifts = duckdb_vector_get_data(vectors[first_column + 2u]);
	if (hgvs->transcript_reason == DUCKVEP_HGVS_ADAPTER_NONE &&
	    hgvs->transcript_length != 0u) {
		duckvep_scalar_assign_ascii_valid(vectors[first_column],
		    validity[first_column], (idx_t)output_row,
		    state->hgvs_text + hgvs->transcript_offset,
		    hgvs->transcript_length);
		hgvs_shifts[output_row] = hgvs->shift_offset;
		duckvep_scalar_set_valid(validity[first_column + 2u],
		    (idx_t)output_row);
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 3u],
		    validity[first_column + 3u], (idx_t)output_row,
		    "supported", sizeof("supported") - 1u);
	} else if (hgvs->transcript_reason ==
	    DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE) {
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 3u],
		    validity[first_column + 3u], (idx_t)output_row,
		    "not_applicable", sizeof("not_applicable") - 1u);
	} else {
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 3u],
		    validity[first_column + 3u], (idx_t)output_row,
		    "unresolved", sizeof("unresolved") - 1u);
		reason = duckvep_scalar_hgvs_reason_name(hgvs->transcript_reason);
		if (reason != NULL)
			duckvep_scalar_assign_ascii_valid(vectors[first_column + 4u],
			    validity[first_column + 4u], (idx_t)output_row,
			    reason, strlen(reason));
	}
	if (hgvs->protein_reason == DUCKVEP_HGVS_ADAPTER_NONE &&
	    hgvs->protein_length != 0u) {
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 1u],
		    validity[first_column + 1u], (idx_t)output_row,
		    state->hgvs_text + hgvs->protein_offset,
		    hgvs->protein_length);
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 5u],
		    validity[first_column + 5u], (idx_t)output_row,
		    "supported", sizeof("supported") - 1u);
	} else if (hgvs->protein_reason ==
	    DUCKVEP_HGVS_ADAPTER_NOT_APPLICABLE) {
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 5u],
		    validity[first_column + 5u], (idx_t)output_row,
		    "not_applicable", sizeof("not_applicable") - 1u);
	} else {
		duckvep_scalar_assign_ascii_valid(vectors[first_column + 5u],
		    validity[first_column + 5u], (idx_t)output_row,
		    "unresolved", sizeof("unresolved") - 1u);
		reason = duckvep_scalar_hgvs_reason_name(hgvs->protein_reason);
		if (reason != NULL)
			duckvep_scalar_assign_ascii_valid(vectors[first_column + 6u],
			    validity[first_column + 6u], (idx_t)output_row,
			    reason, strlen(reason));
	}
	return 1;
}

static int
duckvep_scalar_write_output(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, int with_hgvs,
	char *error, size_t error_size)
{
	duckdb_vector child, vectors[36];
	uint64_t *validity[36];
	uint32_t *transcript_indices, *gene_indices;
	uint32_t *regulation_feature_indices;
	uint32_t *cdna_positions, *cds_positions, *protein_positions;
	bool *nmd_escape_intronless, *nmd_escape_early_cds;
	bool *nmd_escape_last_exon, *nmd_escape_penultimate_exon_end;
	uint64_t *consequence_masks;
	uint32_t *region_masks;
	uint8_t *impact_codes, *status_codes, *reason_codes;
	uint8_t *reference_amino_acids, *alternate_amino_acids;
	uint8_t *nmd_prediction_codes, *nmd_escape_reasons;
	uint8_t *overlap_object_codes;
	size_t output_count, source, row, column, column_count;

	if (!duckvep_scalar_prepare_output_list(state, output, input_rows,
	    &output_count,
	    error, error_size))
		return 0;
	child = duckdb_list_vector_get_child(output);
	column_count = with_hgvs ? 36u : 29u;
	for (column = 0; column < column_count; column++) {
		vectors[column] = duckdb_struct_vector_get_child(child,
		    (idx_t)column);
		validity[column] = duckdb_vector_get_validity(vectors[column]);
	}
	transcript_indices = duckdb_vector_get_data(vectors[0]);
	gene_indices = duckdb_vector_get_data(vectors[1]);
	cdna_positions = duckdb_vector_get_data(vectors[7]);
	cds_positions = duckdb_vector_get_data(vectors[8]);
	protein_positions = duckdb_vector_get_data(vectors[9]);
	nmd_escape_intronless = duckdb_vector_get_data(vectors[13]);
	nmd_escape_early_cds = duckdb_vector_get_data(vectors[14]);
	nmd_escape_last_exon = duckdb_vector_get_data(vectors[15]);
	nmd_escape_penultimate_exon_end = duckdb_vector_get_data(vectors[16]);
	regulation_feature_indices = duckdb_vector_get_data(vectors[17]);
	consequence_masks = duckdb_vector_get_data(vectors[19]);
	region_masks = duckdb_vector_get_data(vectors[20]);
	impact_codes = duckdb_vector_get_data(vectors[21]);
	status_codes = duckdb_vector_get_data(vectors[22]);
	reason_codes = duckdb_vector_get_data(vectors[23]);
	reference_amino_acids = duckdb_vector_get_data(vectors[24]);
	alternate_amino_acids = duckdb_vector_get_data(vectors[25]);
	nmd_prediction_codes = duckdb_vector_get_data(vectors[26]);
	nmd_escape_reasons = duckdb_vector_get_data(vectors[27]);
	overlap_object_codes = duckdb_vector_get_data(vectors[28]);
	duckvep_scalar_set_null_range(vectors[17], &validity[17], output_count);
	for (column = 24; column <= 25; column++)
		duckvep_scalar_set_null_range(vectors[column], &validity[column],
		    output_count);
	if (with_hgvs) {
		for (column = 29; column < 36; column++)
			duckvep_scalar_set_null_range(vectors[column],
			    &validity[column], output_count);
	}
	output_count = 0;
	source = 0;
	for (row = 0; row < (size_t)input_rows; row++) {
		size_t begin, end;

		begin = source;
		while (source < state->result_count &&
		    state->results[source].variant_idx == (uint32_t)row)
			source++;
		end = source;
		if (!duckvep_result_range_has_transcript(state, begin, end)) {
			int complete;

			complete = state->transcript_coverage_complete[row] != 0;
			duckvep_scalar_set_null(vectors[0], &validity[0],
			    (idx_t)output_count);
			duckvep_scalar_set_null(vectors[1], &validity[1],
			    (idx_t)output_count);
			if (complete)
				duckvep_scalar_assign_ascii(vectors[2],
				    (idx_t)output_count, "intergenic_variant",
				    sizeof("intergenic_variant") - 1);
			else
				duckvep_scalar_assign_ascii(vectors[2],
				    (idx_t)output_count, "sequence_variant",
				    sizeof("sequence_variant") - 1);
			duckvep_scalar_assign_ascii(vectors[3], (idx_t)output_count,
			    "MODIFIER", sizeof("MODIFIER") - 1);
			if (complete)
				duckvep_scalar_assign_ascii(vectors[4],
				    (idx_t)output_count, "intergenic",
				    sizeof("intergenic") - 1);
			else
				duckvep_scalar_set_null(vectors[4], &validity[4],
				    (idx_t)output_count);
			duckvep_scalar_assign_ascii(vectors[5], (idx_t)output_count,
			    complete ? "supported" : "unresolved",
			    complete ? sizeof("supported") - 1 :
			    sizeof("unresolved") - 1);
			if (complete)
				duckvep_scalar_set_null(vectors[6], &validity[6],
				    (idx_t)output_count);
			else
				duckvep_scalar_assign_ascii(vectors[6],
				    (idx_t)output_count, "no_feature_in_loaded_model",
				    sizeof("no_feature_in_loaded_model") - 1);
			for (column = 7; column < 17; column++)
				duckvep_scalar_set_null(vectors[column],
				    &validity[column],
				    (idx_t)output_count);
			duckvep_scalar_set_null(vectors[18], &validity[18],
			    (idx_t)output_count);
			consequence_masks[output_count] = complete ?
			    DUCKVEP_SO(DUCKVEP_SO_INTERGENIC) : 0u;
			region_masks[output_count] = 0u;
			impact_codes[output_count] = DUCKVEP_IMPACT_MODIFIER;
			status_codes[output_count] = complete ?
			    DUCKVEP_COMPACT_STATUS_SUPPORTED :
			    DUCKVEP_COMPACT_STATUS_UNRESOLVED;
			reason_codes[output_count] = complete ?
			    DUCKVEP_COMPACT_REASON_NONE :
			    DUCKVEP_COMPACT_REASON_NO_FEATURE;
			nmd_prediction_codes[output_count] =
			    DUCKVEP_NMD_NOT_APPLICABLE;
			nmd_escape_reasons[output_count] = 0u;
			duckvep_scalar_set_null(vectors[28], &validity[28],
			    (idx_t)output_count);
			output_count++;
		}
		for (source = begin; source < end; source++, output_count++) {
			const duckvep_consequence_t *result;
			const char *impact, *status, *reason;
			char consequence[512], region[128], amino_acid[2];
			size_t consequence_length, region_length;

			result = &state->results[source];
			consequence_masks[output_count] = result->consequence_mask;
			region_masks[output_count] = result->region_mask;
			impact_codes[output_count] = result->impact;
			status_codes[output_count] = (result->flags &
			    DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0 ?
			    DUCKVEP_COMPACT_STATUS_UNRESOLVED :
			    DUCKVEP_COMPACT_STATUS_SUPPORTED;
			reason_codes[output_count] = status_codes[output_count] ==
			    DUCKVEP_COMPACT_STATUS_UNRESOLVED ?
			    result->sequence_status : DUCKVEP_COMPACT_REASON_NONE;
			nmd_prediction_codes[output_count] = result->nmd_prediction;
			nmd_escape_reasons[output_count] = result->nmd_escape_reasons;
			overlap_object_codes[output_count] =
			    result->overlap_object_kind;
			if (result->overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT) {
				transcript_indices[output_count] = result->tx_idx;
				gene_indices[output_count] = result->gene_idx;
			} else {
				duckvep_scalar_set_null(vectors[0], &validity[0],
				    (idx_t)output_count);
				duckvep_scalar_set_null(vectors[1], &validity[1],
				    (idx_t)output_count);
				duckvep_scalar_set_valid(validity[17],
				    (idx_t)output_count);
				regulation_feature_indices[output_count] =
				    result->interval_feature_idx;
			}
			{
				const char *overlap_object;

				overlap_object = duckvep_scalar_overlap_object_name(
				    result->overlap_object_kind);
				if (overlap_object == NULL) {
					duckvep_sql_set_error(error, error_size,
					    "duckvep_annotate: invalid overlap object kind");
					return 0;
				}
				duckvep_scalar_assign_ascii(vectors[18],
				    (idx_t)output_count, overlap_object,
				    strlen(overlap_object));
			}
			if (result->consequence_mask == 0u ||
			    !duckvep_scalar_format_terms(result->consequence_mask,
			    consequence, sizeof(consequence), &consequence_length)) {
				(void)snprintf(error, error_size,
				    "duckvep_annotate: could not render consequence mask "
				    "0x%llx for variant %u transcript %u",
				    (unsigned long long)result->consequence_mask,
				    result->variant_idx, result->tx_idx);
				return 0;
			}
			duckvep_scalar_assign_ascii(vectors[2], (idx_t)output_count,
			    consequence, consequence_length);
			impact = duckvep_impact_name(
			    (duckvep_impact_t)result->impact);
			duckvep_scalar_assign_ascii(vectors[3], (idx_t)output_count,
			    impact, strlen(impact));
			if (result->region_mask == 0u) {
				/* A cross-contig breakend can affect a transcript only through
				 * its mate. VEP then has no local topological region to report. */
				duckvep_scalar_set_null(vectors[4], &validity[4],
				    (idx_t)output_count);
			} else {
				if (!duckvep_scalar_format_region(result->region_mask,
				    region, sizeof(region), &region_length)) {
					(void)snprintf(error, error_size,
					    "duckvep_annotate: could not render region mask "
					    "0x%x for variant %u transcript %u",
					    result->region_mask, result->variant_idx,
					    result->tx_idx);
					return 0;
				}
				duckvep_scalar_assign_ascii(vectors[4],
				    (idx_t)output_count, region, region_length);
			}
			status = (result->flags &
			    DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0 ?
			    "unresolved" : "supported";
			duckvep_scalar_assign_ascii(vectors[5], (idx_t)output_count,
			    status, status[0] == 'u' ? sizeof("unresolved") - 1 :
			    sizeof("supported") - 1);
			if (status[0] == 'u') {
				reason = duckvep_scalar_sequence_reason(
				    result->sequence_status);
				duckvep_scalar_assign_ascii(vectors[6],
				    (idx_t)output_count, reason, strlen(reason));
			} else {
				duckvep_scalar_set_null(vectors[6], &validity[6],
				    (idx_t)output_count);
			}
			if (result->cdna_pos >= 0)
				cdna_positions[output_count] =
				    (uint32_t)result->cdna_pos;
			else
				duckvep_scalar_set_null(vectors[7], &validity[7],
				    (idx_t)output_count);
			if (result->cds_pos >= 0)
				cds_positions[output_count] =
				    (uint32_t)result->cds_pos;
			else
				duckvep_scalar_set_null(vectors[8], &validity[8],
				    (idx_t)output_count);
			if (result->protein_pos >= 0)
				protein_positions[output_count] =
				    (uint32_t)result->protein_pos;
			else
				duckvep_scalar_set_null(vectors[9], &validity[9],
				    (idx_t)output_count);
			amino_acid[1] = '\0';
			if (result->aa_ref != 0u) {
				amino_acid[0] = (char)result->aa_ref;
				reference_amino_acids[output_count] = result->aa_ref;
				duckvep_scalar_set_valid(validity[24],
				    (idx_t)output_count);
				duckvep_scalar_assign_ascii(vectors[10],
				    (idx_t)output_count, amino_acid, 1);
			} else {
				duckvep_scalar_set_null(vectors[10], &validity[10],
				    (idx_t)output_count);
			}
			if (result->aa_alt != 0u) {
				amino_acid[0] = (char)result->aa_alt;
				alternate_amino_acids[output_count] = result->aa_alt;
				duckvep_scalar_set_valid(validity[25],
				    (idx_t)output_count);
				duckvep_scalar_assign_ascii(vectors[11],
				    (idx_t)output_count, amino_acid, 1);
			} else {
				duckvep_scalar_set_null(vectors[11], &validity[11],
				    (idx_t)output_count);
			}
			switch ((duckvep_nmd_prediction_t)result->nmd_prediction) {
			case DUCKVEP_NMD_PREDICTED_TRIGGERING:
				duckvep_scalar_assign_ascii(vectors[12],
				    (idx_t)output_count, "triggering",
				    sizeof("triggering") - 1);
				nmd_escape_intronless[output_count] = false;
				nmd_escape_early_cds[output_count] = false;
				nmd_escape_last_exon[output_count] = false;
				nmd_escape_penultimate_exon_end[output_count] = false;
				break;
			case DUCKVEP_NMD_PREDICTED_ESCAPING:
				duckvep_scalar_assign_ascii(vectors[12],
				    (idx_t)output_count, "escaping",
				    sizeof("escaping") - 1);
				nmd_escape_intronless[output_count] =
				    (result->nmd_escape_reasons &
				    DUCKVEP_NMD_ESCAPE_INTRONLESS) != 0;
				nmd_escape_early_cds[output_count] =
				    (result->nmd_escape_reasons &
				    DUCKVEP_NMD_ESCAPE_EARLY_CDS) != 0;
				nmd_escape_last_exon[output_count] =
				    (result->nmd_escape_reasons &
				    DUCKVEP_NMD_ESCAPE_LAST_EXON) != 0;
				nmd_escape_penultimate_exon_end[output_count] =
				    (result->nmd_escape_reasons &
				    DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END) != 0;
				break;
			case DUCKVEP_NMD_UNRESOLVED:
				duckvep_scalar_assign_ascii(vectors[12],
				    (idx_t)output_count, "unresolved",
				    sizeof("unresolved") - 1);
				for (column = 13; column < 17; column++)
					duckvep_scalar_set_null(vectors[column],
					    &validity[column],
					    (idx_t)output_count);
				break;
			case DUCKVEP_NMD_NOT_APPLICABLE:
				for (column = 12; column < 17; column++)
					duckvep_scalar_set_null(vectors[column],
					    &validity[column],
					    (idx_t)output_count);
				break;
			default:
				duckvep_sql_set_error(error, error_size,
				    "duckvep_annotate: invalid kernel NMD prediction");
				return 0;
			}
			if (with_hgvs && result->overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT &&
			    !duckvep_scalar_write_hgvs_fields(state, vectors, validity,
			    29u, output_count, source, error, error_size))
				return 0;
		}
	}
	return 1;
}

static int
duckvep_scalar_write_compact_output(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, int with_hgvs,
	char *error, size_t error_size)
{
	duckdb_vector child, vectors[23];
	uint64_t *validity[23];
	uint32_t *transcript_indices, *gene_indices, *region_masks;
	uint32_t *regulation_feature_indices;
	uint32_t *cdna_positions, *cds_positions, *protein_positions;
	uint64_t *consequence_masks;
	uint8_t *impact_codes, *status_codes, *reason_codes;
	uint8_t *reference_amino_acids, *alternate_amino_acids;
	uint8_t *nmd_prediction_codes, *nmd_escape_reasons;
	uint8_t *overlap_object_codes;
	size_t output_count, source, row, column, column_count;

	if (!duckvep_scalar_prepare_output_list(state, output, input_rows,
	    &output_count,
	    error, error_size))
		return 0;
	child = duckdb_list_vector_get_child(output);
	column_count = with_hgvs ? 23u : 16u;
	for (column = 0; column < column_count; column++) {
		vectors[column] = duckdb_struct_vector_get_child(child,
		    (idx_t)column);
		validity[column] = duckdb_vector_get_validity(vectors[column]);
	}
	transcript_indices = duckdb_vector_get_data(vectors[0]);
	gene_indices = duckdb_vector_get_data(vectors[1]);
	consequence_masks = duckdb_vector_get_data(vectors[2]);
	region_masks = duckdb_vector_get_data(vectors[3]);
	impact_codes = duckdb_vector_get_data(vectors[4]);
	status_codes = duckdb_vector_get_data(vectors[5]);
	reason_codes = duckdb_vector_get_data(vectors[6]);
	cdna_positions = duckdb_vector_get_data(vectors[7]);
	cds_positions = duckdb_vector_get_data(vectors[8]);
	protein_positions = duckdb_vector_get_data(vectors[9]);
	reference_amino_acids = duckdb_vector_get_data(vectors[10]);
	alternate_amino_acids = duckdb_vector_get_data(vectors[11]);
	nmd_prediction_codes = duckdb_vector_get_data(vectors[12]);
	nmd_escape_reasons = duckdb_vector_get_data(vectors[13]);
	regulation_feature_indices = duckdb_vector_get_data(vectors[14]);
	overlap_object_codes = duckdb_vector_get_data(vectors[15]);
	/* cDNA/CDS/protein coordinates and scalar amino acids are absent for the
	 * majority of dense transcript expansion (intronic, UTR, directional,
	 * non-coding, and regulation rows).  Initialize their validity masks once
	 * and mark the smaller resolved subset valid while writing values. */
	for (column = 7; column <= 11; column++)
		duckvep_scalar_set_null_range(vectors[column], &validity[column],
		    output_count);
	duckvep_scalar_set_null_range(vectors[14], &validity[14], output_count);
	if (with_hgvs) {
		for (column = 16; column < 23; column++)
			duckvep_scalar_set_null_range(vectors[column],
			    &validity[column], output_count);
	}
	output_count = 0;
	source = 0;
	for (row = 0; row < (size_t)input_rows; row++) {
		size_t begin, end;

		begin = source;
		while (source < state->result_count &&
		    state->results[source].variant_idx == (uint32_t)row)
			source++;
		end = source;
		if (!duckvep_result_range_has_transcript(state, begin, end)) {
			int complete;

			complete = state->transcript_coverage_complete[row] != 0;
			duckvep_scalar_set_null(vectors[0], &validity[0],
			    (idx_t)output_count);
			duckvep_scalar_set_null(vectors[1], &validity[1],
			    (idx_t)output_count);
			consequence_masks[output_count] = complete ?
			    DUCKVEP_SO(DUCKVEP_SO_INTERGENIC) : 0;
			region_masks[output_count] = 0;
			impact_codes[output_count] = DUCKVEP_IMPACT_MODIFIER;
			status_codes[output_count] = complete ?
			    DUCKVEP_COMPACT_STATUS_SUPPORTED :
			    DUCKVEP_COMPACT_STATUS_UNRESOLVED;
			reason_codes[output_count] = complete ?
			    DUCKVEP_COMPACT_REASON_NONE :
			    DUCKVEP_COMPACT_REASON_NO_FEATURE;
			nmd_prediction_codes[output_count] =
			    DUCKVEP_NMD_NOT_APPLICABLE;
			nmd_escape_reasons[output_count] = 0;
			duckvep_scalar_set_null(vectors[15], &validity[15],
			    (idx_t)output_count);
			output_count++;
		}
		for (source = begin; source < end; source++, output_count++) {
			const duckvep_consequence_t *result;
			int unresolved;

			result = &state->results[source];
			unresolved = (result->flags &
			    DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0;
			if (result->overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT) {
				transcript_indices[output_count] = result->tx_idx;
				gene_indices[output_count] = result->gene_idx;
			} else {
				duckvep_scalar_set_null(vectors[0], &validity[0],
				    (idx_t)output_count);
				duckvep_scalar_set_null(vectors[1], &validity[1],
				    (idx_t)output_count);
				duckvep_scalar_set_valid(validity[14],
				    (idx_t)output_count);
				regulation_feature_indices[output_count] =
				    result->interval_feature_idx;
			}
			overlap_object_codes[output_count] =
			    result->overlap_object_kind;
			consequence_masks[output_count] = result->consequence_mask;
			region_masks[output_count] = result->region_mask;
			impact_codes[output_count] = result->impact;
			status_codes[output_count] = unresolved ?
			    DUCKVEP_COMPACT_STATUS_UNRESOLVED :
			    DUCKVEP_COMPACT_STATUS_SUPPORTED;
			reason_codes[output_count] = unresolved ?
			    result->sequence_status : DUCKVEP_COMPACT_REASON_NONE;
			if (result->cdna_pos >= 0) {
				cdna_positions[output_count] = (uint32_t)result->cdna_pos;
				duckvep_scalar_set_valid(validity[7], (idx_t)output_count);
			}
			if (result->cds_pos >= 0) {
				cds_positions[output_count] = (uint32_t)result->cds_pos;
				duckvep_scalar_set_valid(validity[8], (idx_t)output_count);
			}
			if (result->protein_pos >= 0) {
				protein_positions[output_count] =
				    (uint32_t)result->protein_pos;
				duckvep_scalar_set_valid(validity[9], (idx_t)output_count);
			}
			if (result->aa_ref != 0) {
				reference_amino_acids[output_count] = result->aa_ref;
				duckvep_scalar_set_valid(validity[10], (idx_t)output_count);
			}
			if (result->aa_alt != 0) {
				alternate_amino_acids[output_count] = result->aa_alt;
				duckvep_scalar_set_valid(validity[11], (idx_t)output_count);
			}
			nmd_prediction_codes[output_count] = result->nmd_prediction;
			nmd_escape_reasons[output_count] = result->nmd_escape_reasons;
			if (with_hgvs && result->overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT &&
			    !duckvep_scalar_write_hgvs_fields(state, vectors, validity,
			    16u, output_count, source, error, error_size))
				return 0;
		}
	}
	return 1;
}

typedef enum duckvep_scalar_event_family {
	DUCKVEP_SCALAR_SMALL = 0,
	DUCKVEP_SCALAR_STRUCTURAL,
	DUCKVEP_SCALAR_BREAKEND
} duckvep_scalar_event_family_t;

static void
duckvep_annotate_scalar_execute(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output, int compact, int with_hgvs,
	duckvep_scalar_event_family_t event_family)
{
	duckvep_scalar_state_t *state;
	duckvep_variant_batch_t batch;
	duckdb_vector model_vector, upstream_distance_vector;
	duckdb_vector downstream_distance_vector;
	uint64_t *model_validity, *upstream_distance_validity;
	uint64_t *downstream_distance_validity;
	uint16_t *regions;
	uint32_t *positions;
	uint64_t *upstream_distances;
	uint64_t *downstream_distances;
	idx_t rows;
	idx_t column_count;
	idx_t distance_column;
	size_t begin;
	char error[DUCKVEP_SQL_ERROR_SIZE];

	rows = duckdb_data_chunk_get_size(input);
	if (rows == 0) {
		(void)duckdb_list_vector_set_size(output, 0);
		return;
	}
	memset(error, 0, sizeof(error));
	state = duckvep_scalar_state_acquire(
	    duckdb_scalar_function_get_extra_info(info));
	if (state == NULL) {
		duckdb_scalar_function_set_error(info,
		    "duckvep_annotate: out of memory");
		return;
	}
	if (!duckvep_scalar_result_reserve(state,
	    (size_t)duckdb_vector_size())) {
		duckvep_sql_set_error(error, sizeof(error),
		    "duckvep_annotate: out of memory");
		goto failed;
	}
	model_vector = duckdb_data_chunk_get_vector(input, 0);
	model_validity = duckdb_vector_get_validity(model_vector);
	if (!((event_family == DUCKVEP_SCALAR_STRUCTURAL) ?
	    duckvep_scalar_prepare_sv_batch(state, input, &batch, error,
	    sizeof(error)) :
	    (event_family == DUCKVEP_SCALAR_BREAKEND) ?
	    duckvep_scalar_prepare_breakend_batch(state, input, &batch, error,
	    sizeof(error)) :
	    duckvep_scalar_prepare_batch(state, input, &batch, error,
	    sizeof(error))))
		goto failed;
	regions = state->seq_regions;
	positions = state->positions;
	column_count = duckdb_data_chunk_get_column_count(input);
	distance_column = event_family == DUCKVEP_SCALAR_STRUCTURAL ? 6 : 5;
	upstream_distance_vector = column_count > distance_column ?
	    duckdb_data_chunk_get_vector(input, distance_column) : NULL;
	downstream_distance_vector = column_count > distance_column + 1 ?
	    duckdb_data_chunk_get_vector(input, distance_column + 1) :
	    upstream_distance_vector;
	upstream_distance_validity = upstream_distance_vector != NULL ?
	    duckdb_vector_get_validity(upstream_distance_vector) : NULL;
	downstream_distance_validity = downstream_distance_vector != NULL ?
	    duckdb_vector_get_validity(downstream_distance_vector) : NULL;
	upstream_distances = upstream_distance_vector != NULL ?
	    duckdb_vector_get_data(upstream_distance_vector) : NULL;
	downstream_distances = downstream_distance_vector != NULL ?
	    duckdb_vector_get_data(downstream_distance_vector) : NULL;
	state->result_count = 0;
	state->hgvs_text_size = 0;
	begin = 0;
	while (begin < (size_t)rows) {
		uint64_t upstream_distance;
		uint64_t downstream_distance;
		size_t end;

		if (!duckvep_scalar_select_model(state, model_vector, (idx_t)begin,
		    error, sizeof(error)))
			goto failed;
		if (positions[begin] == 0 ||
		    (upstream_distance_vector != NULL &&
		    duckvep_validity_is_null(upstream_distance_validity,
		    (idx_t)begin)) ||
		    (downstream_distance_vector != NULL &&
		    duckvep_validity_is_null(downstream_distance_validity,
		    (idx_t)begin))) {
			duckvep_sql_set_error(error, sizeof(error),
			    "duckvep_annotate: position and upstream/downstream distances must be non-NULL; position is one-based");
			goto failed;
		}
		upstream_distance = upstream_distances != NULL ?
		    upstream_distances[begin] : DUCKVEP_DEFAULT_UPSTREAM_DIST;
		downstream_distance = downstream_distances != NULL ?
		    downstream_distances[begin] : DUCKVEP_DEFAULT_DOWNSTREAM_DIST;
		end = begin + 1;
		while (end < (size_t)rows) {
			uint64_t next_upstream_distance;
			uint64_t next_downstream_distance;

			if ((upstream_distance_vector != NULL &&
			    duckvep_validity_is_null(upstream_distance_validity,
			    (idx_t)end)) ||
			    (downstream_distance_vector != NULL &&
			    duckvep_validity_is_null(downstream_distance_validity,
			    (idx_t)end)))
				break;
			next_upstream_distance = upstream_distances != NULL ?
			    upstream_distances[end] : DUCKVEP_DEFAULT_UPSTREAM_DIST;
			next_downstream_distance = downstream_distances != NULL ?
			    downstream_distances[end] : DUCKVEP_DEFAULT_DOWNSTREAM_DIST;
			if (!duckvep_vector_string_equals(model_vector, model_validity,
			    (idx_t)end, state->entry->name) ||
			    regions[end] != regions[begin] ||
			    positions[end] < positions[end - 1] ||
			    next_upstream_distance != upstream_distance ||
			    next_downstream_distance != downstream_distance)
				break;
			end++;
		}
		if (event_family == DUCKVEP_SCALAR_BREAKEND) {
			if (!duckvep_scalar_run_breakends(
			    state, &batch, begin, end - begin,
			    upstream_distance, downstream_distance,
			    error, sizeof(error)))
				goto failed;
		} else if (!duckvep_scalar_run(
		    state, &batch, begin, end - begin,
		    upstream_distance, downstream_distance, with_hgvs,
		    error, sizeof(error))) {
			goto failed;
		}
		begin = end;
	}
	if (!(compact ? duckvep_scalar_write_compact_output(state, output, rows,
	    with_hgvs, error, sizeof(error)) : duckvep_scalar_write_output(state,
	    output, rows, with_hgvs, error, sizeof(error))))
		goto failed;
	duckvep_scalar_state_release(state);
	return;

failed:
	duckvep_scalar_state_release(state);
	duckdb_scalar_function_set_error(info,
	    error[0] != '\0' ? error : "duckvep_annotate failed");
}

static void
duckvep_annotate_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 0, 0, 0);
}

static void
duckvep_annotate_compact_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1, 0, 0);
}

static void
duckvep_annotate_hgvs_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1, 1,
	    DUCKVEP_SCALAR_SMALL);
}

static void
duckvep_annotate_rich_hgvs_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 0, 1,
	    DUCKVEP_SCALAR_SMALL);
}

static void
duckvep_annotate_sv_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 0, 0,
	    DUCKVEP_SCALAR_STRUCTURAL);
}

static void
duckvep_annotate_sv_compact_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1, 0,
	    DUCKVEP_SCALAR_STRUCTURAL);
}

static void
duckvep_annotate_breakend_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 0, 0,
	    DUCKVEP_SCALAR_BREAKEND);
}

static void
duckvep_annotate_breakend_compact_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1, 0,
	    DUCKVEP_SCALAR_BREAKEND);
}

static duckdb_logical_type
duckvep_annotation_list_type(int with_hgvs)
{
	const char *names[] = {
		"transcript_index", "gene_index", "consequence", "impact",
		"region", "status", "reason", "cdna_position", "cds_position",
		"protein_position", "reference_amino_acid",
		"alternate_amino_acid", "nmd_prediction",
		"nmd_escape_intronless", "nmd_escape_early_cds",
		"nmd_escape_last_exon", "nmd_escape_penultimate_exon_end",
		"regulation_feature_index", "overlap_object",
		"consequence_mask", "region_mask", "impact_code", "status_code",
		"reason_code", "reference_amino_acid_code",
		"alternate_amino_acid_code", "nmd_prediction_code",
		"nmd_escape_reasons", "overlap_object_code",
		"transcript_hgvs", "protein_hgvs", "hgvs_shift",
		"transcript_hgvs_status", "transcript_hgvs_reason",
		"protein_hgvs_status", "protein_hgvs_reason"
	};
	duckdb_logical_type types[36], structure, list;
	size_t index, column_count;

	types[0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[1] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 2; index <= 6; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	for (index = 7; index <= 9; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[10] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[11] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[12] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	for (index = 13; index < 17; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	types[17] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[18] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[19] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	types[20] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 21; index < 29; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	column_count = with_hgvs ? 36u : 29u;
	if (with_hgvs) {
		types[29] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
		types[30] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
		types[31] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
		for (index = 32; index < 36; index++)
			types[index] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	}
	structure = duckdb_create_struct_type(types, names, column_count);
	list = duckdb_create_list_type(structure);
	for (index = 0; index < column_count; index++)
		duckdb_destroy_logical_type(&types[index]);
	duckdb_destroy_logical_type(&structure);
	return list;
}

static duckdb_logical_type
duckvep_compact_annotation_list_type(void)
{
	const char *names[] = {
		"transcript_index", "gene_index", "consequence_mask",
		"region_mask", "impact_code", "status_code", "reason_code",
		"cdna_position", "cds_position", "protein_position",
		"reference_amino_acid_code", "alternate_amino_acid_code",
		"nmd_prediction_code", "nmd_escape_reasons",
		"regulation_feature_index", "overlap_object_code"
	};
	duckdb_logical_type types[16], structure, list;
	size_t index;

	types[0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[1] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[2] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	types[3] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 4; index <= 6; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	for (index = 7; index <= 9; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 10; index < 14; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	types[14] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[15] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	structure = duckdb_create_struct_type(types, names, 16);
	list = duckdb_create_list_type(structure);
	for (index = 0; index < 16; index++)
		duckdb_destroy_logical_type(&types[index]);
	duckdb_destroy_logical_type(&structure);
	return list;
}

static duckdb_logical_type
duckvep_hgvs_annotation_list_type(void)
{
	const char *names[] = {
		"transcript_index", "gene_index", "consequence_mask",
		"region_mask", "impact_code", "status_code", "reason_code",
		"cdna_position", "cds_position", "protein_position",
		"reference_amino_acid_code", "alternate_amino_acid_code",
		"nmd_prediction_code", "nmd_escape_reasons",
		"regulation_feature_index", "overlap_object_code",
		"transcript_hgvs", "protein_hgvs", "hgvs_shift",
		"transcript_hgvs_status", "transcript_hgvs_reason",
		"protein_hgvs_status", "protein_hgvs_reason"
	};
	duckdb_logical_type types[23], structure, list;
	size_t index;

	types[0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[1] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[2] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	types[3] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 4; index <= 6; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	for (index = 7; index <= 9; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 10; index < 14; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	types[14] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[15] = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	types[16] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[17] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[18] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 19; index < 23; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	structure = duckdb_create_struct_type(types, names, 23);
	list = duckdb_create_list_type(structure);
	for (index = 0; index < 23; index++)
		duckdb_destroy_logical_type(&types[index]);
	duckdb_destroy_logical_type(&structure);
	return list;
}

static void
duckvep_register_annotate_scalar(duckdb_connection connection,
	duckvep_registry_t *registry, duckdb_logical_type varchar_type,
	duckdb_logical_type uinteger_type, duckdb_logical_type ubigint_type,
	int distance_parameters, int compact,
	duckvep_scalar_event_family_t event_family)
{
	duckdb_scalar_function scalar;
	duckdb_logical_type result_type;
	const char *name;

	scalar = duckdb_create_scalar_function();
	result_type = compact ? duckvep_compact_annotation_list_type() :
	    duckvep_annotation_list_type(0);
	if (event_family == DUCKVEP_SCALAR_STRUCTURAL)
		name = compact ? "_duckvep_annotate_structural_compact" :
		    "_duckvep_annotate_structural_rich";
	else if (event_family == DUCKVEP_SCALAR_BREAKEND)
		name = compact ? "_duckvep_annotate_breakend_compact" :
		    "_duckvep_annotate_breakend_rich";
	else
		name = compact ? "_duckvep_annotate_small_compact" :
		    "_duckvep_annotate_small_rich";
	duckdb_scalar_function_set_name(scalar, name);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, uinteger_type);
	duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	if (event_family == DUCKVEP_SCALAR_STRUCTURAL) {
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
	} else if (event_family == DUCKVEP_SCALAR_BREAKEND) {
		duckdb_scalar_function_add_parameter(scalar, uinteger_type);
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	} else {
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
	}
	if (distance_parameters >= 1)
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	if (distance_parameters >= 2)
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	duckdb_scalar_function_set_return_type(scalar, result_type);
	duckdb_scalar_function_set_volatile(scalar);
	duckvep_registry_retain(registry);
	duckdb_scalar_function_set_extra_info(scalar, registry,
	    duckvep_registry_release);
	if (event_family == DUCKVEP_SCALAR_STRUCTURAL)
		duckdb_scalar_function_set_function(scalar, compact ?
		    duckvep_annotate_sv_compact_scalar : duckvep_annotate_sv_scalar);
	else if (event_family == DUCKVEP_SCALAR_BREAKEND)
		duckdb_scalar_function_set_function(scalar, compact ?
		    duckvep_annotate_breakend_compact_scalar :
		    duckvep_annotate_breakend_scalar);
	else
		duckdb_scalar_function_set_function(scalar, compact ?
		    duckvep_annotate_compact_scalar : duckvep_annotate_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);
	duckdb_destroy_logical_type(&result_type);
}

static void
duckvep_register_hgvs_scalar(duckdb_connection connection,
	duckvep_registry_t *registry, duckdb_logical_type varchar_type,
	duckdb_logical_type uinteger_type, duckdb_logical_type ubigint_type,
	int distance_parameters, int rich)
{
	duckdb_scalar_function scalar;
	duckdb_logical_type result_type;

	scalar = duckdb_create_scalar_function();
	result_type = rich ? duckvep_annotation_list_type(1) :
	    duckvep_hgvs_annotation_list_type();
	duckdb_scalar_function_set_name(scalar,
	    rich ? "_duckvep_annotate_small_rich_hgvs" :
	    "_duckvep_annotate_small_hgvs");
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, uinteger_type);
	duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	if (distance_parameters >= 1)
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	if (distance_parameters >= 2)
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	duckdb_scalar_function_set_return_type(scalar, result_type);
	duckdb_scalar_function_set_volatile(scalar);
	duckvep_registry_retain(registry);
	duckdb_scalar_function_set_extra_info(scalar, registry,
	    duckvep_registry_release);
	duckdb_scalar_function_set_function(scalar, rich ?
	    duckvep_annotate_rich_hgvs_scalar : duckvep_annotate_hgvs_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);
	duckdb_destroy_logical_type(&result_type);
}

void
register_duckvep_functions(duckdb_connection connection,
	duckdb_database database)
{
	duckvep_registry_t *registry;
	duckdb_logical_type varchar_type, uinteger_type, ubigint_type;

	registry = duckvep_registry_create(database);
	if (registry == NULL)
		return;
	duckvep_register_model_functions(connection, registry);
	duckvep_register_allele_geometry_scalar(connection);

	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 0, DUCKVEP_SCALAR_SMALL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 0, DUCKVEP_SCALAR_SMALL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 0, DUCKVEP_SCALAR_SMALL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 1, DUCKVEP_SCALAR_SMALL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 1, DUCKVEP_SCALAR_SMALL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 1, DUCKVEP_SCALAR_SMALL);
	duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 0);
	duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 0);
	duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 0);
	duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 1);
	duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 1);
	duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 1);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 0, DUCKVEP_SCALAR_STRUCTURAL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 0, DUCKVEP_SCALAR_STRUCTURAL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 0, DUCKVEP_SCALAR_STRUCTURAL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 1, DUCKVEP_SCALAR_STRUCTURAL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 1, DUCKVEP_SCALAR_STRUCTURAL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 1, DUCKVEP_SCALAR_STRUCTURAL);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 0, DUCKVEP_SCALAR_BREAKEND);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 0, DUCKVEP_SCALAR_BREAKEND);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 0, DUCKVEP_SCALAR_BREAKEND);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, 1, DUCKVEP_SCALAR_BREAKEND);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, 1, DUCKVEP_SCALAR_BREAKEND);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 2, 1, DUCKVEP_SCALAR_BREAKEND);
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&uinteger_type);
	duckdb_destroy_logical_type(&ubigint_type);

	duckvep_registry_release(registry);
}
