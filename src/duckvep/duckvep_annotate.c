/* Stable-ABI DuckDB vector adapter for the resident DuckVEP kernel. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckvep_model.h"
#include "kernel/src/duckvep_effect.h"
#include "kernel/src/duckvep_event.h"
#include "kernel/include/duckvep_kernel.h"
#include "kernel/include/duckvep_so.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
	uint32_t *pair_transcript_indices;
	size_t pair_capacity;
	duckvep_consequence_t *results;
	size_t result_count;
	size_t result_capacity;
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
	    !duckvep_sql_resize((void **)&state->pair_transcript_indices,
	    sizeof(*state->pair_transcript_indices), capacity))
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

static void
duckvep_scalar_state_destroy(void *pointer)
{
	duckvep_scalar_state_t *state;

	state = pointer;
	if (state == NULL)
		return;
	if (state->workspace_cache != NULL) {
		duckvep_workspace_close(state->workspace_cache->workspace);
		free(state->workspace_cache);
	} else
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
	free(state->pair_transcript_indices);
	free(state->results);
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
		if (!duckvep_scalar_dna_copy(state->allele_bytes +
		    state->reference_offsets[row],
		    duckdb_string_t_data(&references[row]),
		    state->reference_lengths[row]) ||
		    !duckvep_scalar_dna_copy(state->allele_bytes +
		    state->alternate_offsets[row],
		    duckdb_string_t_data(&alternates[row]),
		    state->alternate_lengths[row]) ||
		    !duckvep_scalar_simple_kind(state->positions[row],
		    state->allele_bytes +
		    state->reference_offsets[row], state->reference_lengths[row],
		    state->allele_bytes + state->alternate_offsets[row],
		    state->alternate_lengths[row], &state->variant_kinds[row])) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: REF/ALT must be unequal A/C/G/T/N alleles");
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
		    "duckvep_annotate_sv: out of memory");
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
			    "duckvep_annotate_sv: event fields cannot be NULL");
			return 0;
		}
		if (regions[row] > UINT16_MAX || starts[row] == 0u ||
		    starts[row] > UINT32_MAX || ends[row] < starts[row] ||
		    ends[row] > UINT32_MAX) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_sv: invalid 1-based inclusive interval");
			return 0;
		}
		if (!duckvep_scalar_sv_type(&types[row], &sv_type) ||
		    !duckvep_scalar_copy_change(&copies[row], &copy_change)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_sv: unknown structural type or copy change");
			return 0;
		}
		if (sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_sv: use duckvep_annotate_breakend with both loci");
			return 0;
		}
		if (!duckvep_sv_metadata_valid((duckvep_sv_type_t)sv_type,
		    (duckvep_copy_change_t)copy_change)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_sv: structural type contradicts copy change");
			return 0;
		}
		if (!duckvep_sv_geometry_valid((duckvep_sv_type_t)sv_type,
		    (uint32_t)starts[row], (uint32_t)ends[row])) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_sv: invalid structural event geometry");
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
		    "duckvep_annotate_breakend: out of memory");
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
			    "duckvep_annotate_breakend: both loci are required");
			return 0;
		}
		if (regions[row] > UINT16_MAX || mate_regions[row] > UINT16_MAX ||
		    positions[row] == 0u || positions[row] >= UINT32_MAX ||
		    mate_positions[row] == 0u || mate_positions[row] > UINT32_MAX) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_breakend: invalid one-based endpoint coordinate");
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
duckvep_scalar_u32_compare(const void *left, const void *right)
{
	uint32_t a, b;

	a = *(const uint32_t *)left;
	b = *(const uint32_t *)right;
	return a < b ? -1 : a > b;
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
	return duckvep_scalar_seed_index(state, model->interval_index,
	    model->interval_index_complete, seq_region, position, halo_distance,
	    cursor, 0, error, error_size) &&
	    duckvep_scalar_seed_index(state, model->interval_feature_index,
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
	char *error, size_t error_size)
{
	duckvep_variant_batch_t slice;
	duckvep_annotate_cursor_t *cursor;
	duckvep_options_init_t options_init;
	duckvep_error_t kernel_error;
	duckvep_status_t status;
	uint32_t halo_distance;
	uint32_t sequence_length;
	size_t variant;

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
	while (!duckvep_annotate_cursor_done(cursor)) {
		duckvep_result_builder_t builder;
		size_t old_count, index, write;

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
		write = old_count;
		for (index = 0; index < builder.count; index++) {
			duckvep_consequence_t row;

			row = builder.rows[index];
			row.variant_idx += (uint32_t)begin;
			if (row.overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT)
				row.gene_idx =
				    state->entry->model.gene_indices[row.tx_idx];
			state->results[write++] = row;
		}
		state->result_count = write;
		if (status == DUCKVEP_ERR_RESULT_FULL)
			continue;
		if (status != DUCKVEP_OK) {
			duckvep_annotate_cursor_close(cursor);
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
	uint16_t seq_region, uint32_t position, uint32_t halo_distance,
	size_t *candidate_count, char *error, size_t error_size)
{
	duckvep_owned_model_t *model;
	int32_t query_start, query_end;
	int64_t hit_count, hit;
	char region_name[16];

	model = &state->entry->model;
	if (!model->interval_index_complete || position > (uint32_t)INT32_MAX ||
	    halo_distance > (uint32_t)INT32_MAX ||
	    position > (uint32_t)INT32_MAX - halo_distance) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate_breakend: resident transcript index cannot represent an endpoint");
		return 0;
	}
	query_start = position > halo_distance + 1u ?
	    (int32_t)(position - halo_distance - 1u) : 0;
	query_end = (int32_t)(position + halo_distance);
	(void)snprintf(region_name, sizeof(region_name), "%u",
	    (unsigned)seq_region);
	hit_count = cr_overlap(model->interval_index, region_name, query_start,
	    query_end, &state->interval_hits, &state->interval_hit_capacity);
	if (hit_count < 0 || (uint64_t)hit_count > SIZE_MAX - *candidate_count ||
	    !duckvep_scalar_seed_reserve(state,
	    *candidate_count + (size_t)hit_count)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate_breakend: could not collect endpoint transcripts");
		return 0;
	}
	for (hit = 0; hit < hit_count; hit++)
		state->seed_transcripts[(*candidate_count)++] =
		    (uint32_t)cr_label(model->interval_index,
		    state->interval_hits[hit]);
	return 1;
}

static int
duckvep_scalar_run_breakends(duckvep_scalar_state_t *state,
	const duckvep_variant_batch_t *batch, size_t begin, size_t count,
	uint64_t upstream_distance, uint64_t downstream_distance,
	char *error, size_t error_size)
{
	duckvep_variant_batch_t slice;
	duckvep_candidate_pairs_t pairs;
	duckvep_options_init_t options_init;
	duckvep_error_t kernel_error;
	duckvep_result_builder_t builder;
	duckvep_status_t status;
	uint32_t halo_distance;
	uint32_t search_distance;
	size_t pair_count, variant, old_count, result;

	if (upstream_distance > UINT32_MAX || downstream_distance > UINT32_MAX) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate_breakend: transcript distance exceeds the uint32 kernel limit");
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
			    "duckvep_annotate_breakend: endpoint region is absent from the loaded model");
			return 0;
		}
		if ((local_length != 0u && batch->pos1[variant] > local_length) ||
		    (mate_length != 0u && batch->mate_pos1[variant] > mate_length)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate_breakend: endpoint exceeds sequence-region length");
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
		memset(&kernel_error, 0, sizeof(kernel_error));
		if (duckvep_options_open(&options_init, &state->options,
		    &kernel_error) != DUCKVEP_OK) {
			(void)snprintf(error, error_size,
			    "duckvep_annotate_breakend: %s", kernel_error.message);
			return 0;
		}
		state->options_upstream_distance = (uint32_t)upstream_distance;
		state->options_downstream_distance = (uint32_t)downstream_distance;
		state->have_options_distances = 1;
	}

	pair_count = 0u;
	for (variant = 0u; variant < count; variant++) {
		size_t candidate_count = 0u;
		size_t candidate, unique_count;
		uint32_t local_point = batch->pos1[begin + variant] + 1u;

		if (!duckvep_scalar_append_endpoint_candidates(state,
		    batch->chrom_id[begin + variant], local_point,
		    search_distance, &candidate_count, error, error_size) ||
		    !duckvep_scalar_append_endpoint_candidates(state,
		    batch->mate_chrom_id[begin + variant],
		    batch->mate_pos1[begin + variant], search_distance,
		    &candidate_count, error, error_size))
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
			    "duckvep_annotate_breakend: out of memory collecting candidate pairs");
			return 0;
		}
		for (candidate = 0u; candidate < unique_count; candidate++) {
			state->pair_variant_indices[pair_count] = (uint32_t)variant;
			state->pair_transcript_indices[pair_count] =
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
		    "duckvep_annotate_breakend: out of memory growing results");
		return 0;
	}
	slice = duckvep_scalar_batch_slice(batch, begin, count);
	pairs.variant_idx = state->pair_variant_indices;
	pairs.tx_idx = state->pair_transcript_indices;
	pairs.count = pair_count;
	duckvep_result_builder_init(&builder, state->results + old_count,
	    state->result_capacity - old_count);
	memset(&kernel_error, 0, sizeof(kernel_error));
	status = duckvep_annotate_pairs(state->entry->model.kernel, &slice,
	    &pairs, state->options, state->workspace, &builder, &kernel_error);
	if (status != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "duckvep_annotate_breakend: %s",
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
duckvep_scalar_assign_ascii(duckdb_vector vector, idx_t row,
	const char *text, size_t length)
{
	duckdb_vector_assign_string_element_len(vector, row, text,
	    (idx_t)length);
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
	size_t source;

	for (source = begin; source < end; source++) {
		if (state->results[source].overlap_object_kind ==
		    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT)
			return 1;
	}
	return 0;
}

static int
duckvep_scalar_prepare_output_list(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, char *error, size_t error_size)
{
	duckdb_list_entry *lists;
	size_t output_count, source, row;

	output_count = 0;
	source = 0;
	for (row = 0; row < (size_t)input_rows; row++) {
		size_t begin;

		begin = source;
		while (source < state->result_count &&
		    state->results[source].variant_idx == (uint32_t)row)
			source++;
		output_count += source - begin;
		if (!duckvep_result_range_has_transcript(state, begin, source))
			output_count++;
	}
	if (source != state->result_count) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: result rows are out of order");
		return 0;
	}
	if (duckdb_list_vector_reserve(output, (idx_t)output_count) ==
	    DuckDBError || duckdb_list_vector_set_size(output,
	    (idx_t)output_count) == DuckDBError) {
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
	return 1;
}

static int
duckvep_scalar_write_output(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, char *error, size_t error_size)
{
	duckdb_vector child, vectors[19];
	uint64_t *validity[19];
	uint32_t *transcript_indices, *gene_indices;
	uint32_t *regulation_feature_indices;
	uint64_t *cdna_positions, *cds_positions, *protein_positions;
	bool *nmd_escape_intronless, *nmd_escape_early_cds;
	bool *nmd_escape_last_exon, *nmd_escape_penultimate_exon_end;
	size_t output_count, source, row, column;

	if (!duckvep_scalar_prepare_output_list(state, output, input_rows,
	    error, error_size))
		return 0;
	child = duckdb_list_vector_get_child(output);
	for (column = 0; column < 19; column++) {
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
			for (column = 7; column < 19; column++)
				duckvep_scalar_set_null(vectors[column],
				    &validity[column],
				    (idx_t)output_count);
			output_count++;
		}
		for (source = begin; source < end; source++, output_count++) {
			const duckvep_consequence_t *result;
			const char *impact, *status, *reason;
			char consequence[512], region[128], amino_acid[2];
			size_t consequence_length, region_length;

			result = &state->results[source];
			if (result->overlap_object_kind ==
			    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT) {
				transcript_indices[output_count] = result->tx_idx;
				gene_indices[output_count] = result->gene_idx;
				duckvep_scalar_set_null(vectors[17], &validity[17],
				    (idx_t)output_count);
			} else {
				duckvep_scalar_set_null(vectors[0], &validity[0],
				    (idx_t)output_count);
				duckvep_scalar_set_null(vectors[1], &validity[1],
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
				    (uint64_t)result->cdna_pos;
			else
				duckvep_scalar_set_null(vectors[7], &validity[7],
				    (idx_t)output_count);
			if (result->cds_pos >= 0)
				cds_positions[output_count] =
				    (uint64_t)result->cds_pos;
			else
				duckvep_scalar_set_null(vectors[8], &validity[8],
				    (idx_t)output_count);
			if (result->protein_pos >= 0)
				protein_positions[output_count] =
				    (uint64_t)result->protein_pos;
			else
				duckvep_scalar_set_null(vectors[9], &validity[9],
				    (idx_t)output_count);
			amino_acid[1] = '\0';
			if (result->aa_ref != 0u) {
				amino_acid[0] = (char)result->aa_ref;
				duckvep_scalar_assign_ascii(vectors[10],
				    (idx_t)output_count, amino_acid, 1);
			} else {
				duckvep_scalar_set_null(vectors[10], &validity[10],
				    (idx_t)output_count);
			}
			if (result->aa_alt != 0u) {
				amino_acid[0] = (char)result->aa_alt;
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
		}
	}
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
duckvep_scalar_write_compact_output(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, char *error, size_t error_size)
{
	duckdb_vector child, vectors[16];
	uint64_t *validity[16];
	uint32_t *transcript_indices, *gene_indices, *region_masks;
	uint32_t *regulation_feature_indices;
	uint32_t *cdna_positions, *cds_positions, *protein_positions;
	uint64_t *consequence_masks;
	uint8_t *impact_codes, *status_codes, *reason_codes;
	uint8_t *reference_amino_acids, *alternate_amino_acids;
	uint8_t *nmd_prediction_codes, *nmd_escape_reasons;
	uint8_t *overlap_object_codes;
	size_t output_count, source, row, column;

	if (!duckvep_scalar_prepare_output_list(state, output, input_rows,
	    error, error_size))
		return 0;
	child = duckdb_list_vector_get_child(output);
	for (column = 0; column < 16; column++) {
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
			for (column = 7; column <= 11; column++)
				duckvep_scalar_set_null(vectors[column],
				    &validity[column], (idx_t)output_count);
			nmd_prediction_codes[output_count] =
			    DUCKVEP_NMD_NOT_APPLICABLE;
			nmd_escape_reasons[output_count] = 0;
			duckvep_scalar_set_null(vectors[14], &validity[14],
			    (idx_t)output_count);
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
				duckvep_scalar_set_null(vectors[14], &validity[14],
				    (idx_t)output_count);
			} else {
				duckvep_scalar_set_null(vectors[0], &validity[0],
				    (idx_t)output_count);
				duckvep_scalar_set_null(vectors[1], &validity[1],
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
			if (result->cdna_pos >= 0)
				cdna_positions[output_count] = (uint32_t)result->cdna_pos;
			else
				duckvep_scalar_set_null(vectors[7], &validity[7],
				    (idx_t)output_count);
			if (result->cds_pos >= 0)
				cds_positions[output_count] = (uint32_t)result->cds_pos;
			else
				duckvep_scalar_set_null(vectors[8], &validity[8],
				    (idx_t)output_count);
			if (result->protein_pos >= 0)
				protein_positions[output_count] =
				    (uint32_t)result->protein_pos;
			else
				duckvep_scalar_set_null(vectors[9], &validity[9],
				    (idx_t)output_count);
			if (result->aa_ref != 0)
				reference_amino_acids[output_count] = result->aa_ref;
			else
				duckvep_scalar_set_null(vectors[10], &validity[10],
				    (idx_t)output_count);
			if (result->aa_alt != 0)
				alternate_amino_acids[output_count] = result->aa_alt;
			else
				duckvep_scalar_set_null(vectors[11], &validity[11],
				    (idx_t)output_count);
			nmd_prediction_codes[output_count] = result->nmd_prediction;
			nmd_escape_reasons[output_count] = result->nmd_escape_reasons;
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
	duckdb_data_chunk input, duckdb_vector output, int compact,
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
		if (!((event_family == DUCKVEP_SCALAR_BREAKEND ?
		    duckvep_scalar_run_breakends : duckvep_scalar_run)(
		    state, &batch, begin, end - begin,
		    upstream_distance, downstream_distance,
		    error, sizeof(error))))
			goto failed;
		begin = end;
	}
	if (!(compact ? duckvep_scalar_write_compact_output(state, output, rows,
	    error, sizeof(error)) : duckvep_scalar_write_output(state, output,
	    rows, error, sizeof(error))))
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
	duckvep_annotate_scalar_execute(info, input, output, 0, 0);
}

static void
duckvep_annotate_compact_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1, 0);
}

static void
duckvep_annotate_sv_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 0,
	    DUCKVEP_SCALAR_STRUCTURAL);
}

static void
duckvep_annotate_sv_compact_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1,
	    DUCKVEP_SCALAR_STRUCTURAL);
}

static void
duckvep_annotate_breakend_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 0,
	    DUCKVEP_SCALAR_BREAKEND);
}

static void
duckvep_annotate_breakend_compact_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_annotate_scalar_execute(info, input, output, 1,
	    DUCKVEP_SCALAR_BREAKEND);
}

static duckdb_logical_type
duckvep_annotation_list_type(void)
{
	const char *names[] = {
		"transcript_index", "gene_index", "consequence", "impact",
		"region", "status", "reason", "cdna_position", "cds_position",
		"protein_position", "reference_amino_acid",
		"alternate_amino_acid", "nmd_prediction",
		"nmd_escape_intronless", "nmd_escape_early_cds",
		"nmd_escape_last_exon", "nmd_escape_penultimate_exon_end",
		"regulation_feature_index", "overlap_object"
	};
	duckdb_logical_type types[19], structure, list;
	size_t index;

	types[0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[1] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 2; index <= 6; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	for (index = 7; index <= 9; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	types[10] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[11] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[12] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	for (index = 13; index < 17; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	types[17] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[18] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	structure = duckdb_create_struct_type(types, names, 19);
	list = duckdb_create_list_type(structure);
	for (index = 0; index < 19; index++)
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
	    duckvep_annotation_list_type();
	if (event_family == DUCKVEP_SCALAR_STRUCTURAL)
		name = compact ? "duckvep_annotate_sv_compact" :
		    "duckvep_annotate_sv";
	else if (event_family == DUCKVEP_SCALAR_BREAKEND)
		name = compact ? "duckvep_annotate_breakend_compact" :
		    "duckvep_annotate_breakend";
	else
		name = compact ? "duckvep_annotate_compact" :
		    "duckvep_annotate";
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
