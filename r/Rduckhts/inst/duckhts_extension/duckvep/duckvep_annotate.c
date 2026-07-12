/* Stable-ABI DuckDB vector adapter for the resident DuckVEP kernel. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckvep_model.h"
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
	uint32_t options_distance;
	int have_options_distance;
	int64_t *interval_hits;
	int64_t interval_hit_capacity;
	uint32_t *seed_transcripts;
	size_t seed_capacity;
	uint16_t *seq_regions;
	uint32_t *positions;
	uint32_t *ends;
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
	free(state->reference_offsets);
	free(state->reference_lengths);
	free(state->alternate_offsets);
	free(state->alternate_lengths);
	free(state->variant_kinds);
	free(state->sv_types);
	free(state->copy_changes);
	free(state->transcript_coverage_complete);
	free(state->allele_bytes);
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
duckvep_vector_string_equals(duckdb_vector vector, idx_t row,
	const char *text)
{
	duckdb_string_t *strings;
	size_t length;

	if (text == NULL || duckvep_row_is_null(vector, row))
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
duckvep_scalar_prepare_batch(duckvep_scalar_state_t *state,
	duckdb_data_chunk input, duckvep_variant_batch_t *batch,
	char *error, size_t error_size)
{
	duckdb_vector region_vector, position_vector;
	duckdb_vector reference_vector, alternate_vector;
	duckdb_string_t *references, *alternates;
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
	allele_size = 0;
	for (row = 0; row < rows; row++) {
		size_t reference_length, alternate_length;
		uint32_t region;
		uint64_t position;

		if (duckvep_row_is_null(region_vector, row) ||
		    duckvep_row_is_null(position_vector, row) ||
		    duckvep_row_is_null(reference_vector, row) ||
		    duckvep_row_is_null(alternate_vector, row)) {
			duckvep_sql_set_error(error, error_size,
			    "duckvep_annotate: variant fields cannot be NULL");
			return 0;
		}
		region = ((uint32_t *)duckdb_vector_get_data(region_vector))[row];
		position = ((uint64_t *)duckdb_vector_get_data(position_vector))[row];
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
duckvep_scalar_seed_cursor(duckvep_scalar_state_t *state,
	uint16_t seq_region, uint32_t position, uint32_t distance,
	duckvep_annotate_cursor_t *cursor, char *error, size_t error_size)
{
	duckvep_owned_model_t *model;
	duckvep_error_t kernel_error;
	int32_t query_start, query_end;
	int64_t hit_count, hit;
	char region_name[16];

	model = &state->entry->model;
	if (!model->interval_index_complete ||
	    position > (uint32_t)INT32_MAX || distance >= (uint32_t)INT32_MAX ||
	    position > (uint32_t)INT32_MAX - distance)
		return 1;
	query_start = position > distance + 1 ?
	    (int32_t)(position - distance - 1) : 0;
	query_end = (int32_t)(position + distance);
	(void)snprintf(region_name, sizeof(region_name), "%u",
	    (unsigned)seq_region);
	hit_count = cr_overlap(model->interval_index, region_name, query_start,
	    query_end, &state->interval_hits, &state->interval_hit_capacity);
	if (hit_count < 0 ||
	    !duckvep_scalar_seed_reserve(state, (size_t)hit_count)) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: could not seed the sorted transcript sweep");
		return 0;
	}
	for (hit = 0; hit < hit_count; hit++)
		state->seed_transcripts[hit] = (uint32_t)cr_label(
		    model->interval_index, state->interval_hits[hit]);
	qsort(state->seed_transcripts, (size_t)hit_count,
	    sizeof(*state->seed_transcripts), duckvep_scalar_u32_compare);
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (duckvep_annotate_cursor_seed(cursor, state->seed_transcripts,
	    (size_t)hit_count, &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "duckvep_annotate: %s",
		    kernel_error.message);
		return 0;
	}
	return 1;
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
	slice.ref_offset += begin;
	slice.ref_length += begin;
	slice.alt_offset += begin;
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
	uint64_t distance, char *error, size_t error_size)
{
	duckvep_variant_batch_t slice;
	duckvep_annotate_cursor_t *cursor;
	duckvep_options_init_t options_init;
	duckvep_error_t kernel_error;
	duckvep_status_t status;
	uint32_t sweep_distance;
	uint32_t sequence_length;
	size_t variant;

	if (distance > UINT32_MAX) {
		duckvep_sql_set_error(error, error_size,
		    "duckvep_annotate: distance exceeds the uint32 kernel limit");
		return 0;
	}
	/* The kernel uses zero to request its 5 kb default. For an exact public
	 * query, sweep one base on either side and discard the directional rows.
	 * This keeps the cgranges seed and avoids restarting at the first transcript
	 * for every exact-position DuckDB vector. */
	sweep_distance = distance == 0 ? 1u : (uint32_t)distance;
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
	if (!state->have_options_distance ||
	    state->options_distance != sweep_distance) {
		duckvep_options_close(state->options);
		state->options = NULL;
		state->have_options_distance = 0;
		memset(&options_init, 0, sizeof(options_init));
		options_init.upstream_dist = sweep_distance;
		options_init.downstream_dist = sweep_distance;
		options_init.halo = sweep_distance;
		memset(&kernel_error, 0, sizeof(kernel_error));
		if (duckvep_options_open(&options_init, &state->options,
		    &kernel_error) != DUCKVEP_OK) {
			(void)snprintf(error, error_size, "duckvep_annotate: %s",
			    kernel_error.message);
			return 0;
		}
		state->options_distance = sweep_distance;
		state->have_options_distance = 1;
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
	    batch->pos1[begin], sweep_distance, cursor, error, error_size)) {
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
			if (distance == 0 && (row.region_mask &
			    (DUCKVEP_REGION_UPSTREAM |
			    DUCKVEP_REGION_DOWNSTREAM)) != 0)
				continue;
			row.variant_idx += (uint32_t)begin;
			row.gene_idx = state->entry->model.gene_indices[row.tx_idx];
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

static void
duckvep_scalar_set_null(duckdb_vector vector, idx_t row)
{
	uint64_t *validity;

	duckdb_vector_ensure_validity_writable(vector);
	validity = duckdb_vector_get_validity(vector);
	validity[row / 64] &= ~(UINT64_C(1) << (row % 64));
}

static int
duckvep_scalar_format_terms(uint64_t terms, char *buffer, size_t size)
{
	return size != 0 && duckvep_so_render(terms, '&', buffer, size) < size;
}

static int
duckvep_scalar_format_region(uint32_t region, char *buffer, size_t size)
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
	case DUCKVEP_SEQUENCE_NOT_APPLICABLE:
	case DUCKVEP_SEQUENCE_RESOLVED:
	case DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT:
	default:
		return "unsupported_compound_consequence";
	}
}

static int
duckvep_scalar_write_output(duckvep_scalar_state_t *state,
	duckdb_vector output, idx_t input_rows, char *error, size_t error_size)
{
	duckdb_list_entry *lists;
	duckdb_vector child, vectors[12];
	uint32_t *transcript_indices, *gene_indices;
	uint64_t *cdna_positions, *cds_positions, *protein_positions;
	size_t output_count, source, row, column;

	output_count = 0;
	source = 0;
	for (row = 0; row < (size_t)input_rows; row++) {
		size_t begin;

		begin = source;
		while (source < state->result_count &&
		    state->results[source].variant_idx == (uint32_t)row)
			source++;
		output_count += source == begin ? 1 : source - begin;
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
		lists[row].length = source == begin ? 1 : source - begin;
		output_count += lists[row].length;
	}
	child = duckdb_list_vector_get_child(output);
	for (column = 0; column < 12; column++)
		vectors[column] = duckdb_struct_vector_get_child(child,
		    (idx_t)column);
	transcript_indices = duckdb_vector_get_data(vectors[0]);
	gene_indices = duckdb_vector_get_data(vectors[1]);
	cdna_positions = duckdb_vector_get_data(vectors[7]);
	cds_positions = duckdb_vector_get_data(vectors[8]);
	protein_positions = duckdb_vector_get_data(vectors[9]);
	output_count = 0;
	source = 0;
	for (row = 0; row < (size_t)input_rows; row++) {
		size_t begin, end;

		begin = source;
		while (source < state->result_count &&
		    state->results[source].variant_idx == (uint32_t)row)
			source++;
		end = source;
		if (begin == end) {
			duckvep_scalar_set_null(vectors[0], (idx_t)output_count);
			duckvep_scalar_set_null(vectors[1], (idx_t)output_count);
			duckdb_vector_assign_string_element(vectors[2],
			    (idx_t)output_count, state->transcript_coverage_complete[row]
			    ? "intergenic_variant" : "sequence_variant");
			duckdb_vector_assign_string_element(vectors[3],
			    (idx_t)output_count, "MODIFIER");
			if (state->transcript_coverage_complete[row])
				duckdb_vector_assign_string_element(vectors[4],
				    (idx_t)output_count, "intergenic");
			else
				duckvep_scalar_set_null(vectors[4],
				    (idx_t)output_count);
			duckdb_vector_assign_string_element(vectors[5],
			    (idx_t)output_count, state->transcript_coverage_complete[row]
			    ? "supported" : "unresolved");
			if (state->transcript_coverage_complete[row])
				duckvep_scalar_set_null(vectors[6],
				    (idx_t)output_count);
			else
				duckdb_vector_assign_string_element(vectors[6],
				    (idx_t)output_count, "no_feature_in_loaded_model");
			for (column = 7; column < 12; column++)
				duckvep_scalar_set_null(vectors[column],
				    (idx_t)output_count);
			output_count++;
			continue;
		}
		for (source = begin; source < end; source++, output_count++) {
			const duckvep_consequence_t *result;
			const char *status, *reason;
			char consequence[512], region[128], amino_acid[2];

			result = &state->results[source];
			transcript_indices[output_count] = result->tx_idx;
			gene_indices[output_count] = result->gene_idx;
			if (!duckvep_scalar_format_terms(result->consequence_mask,
			    consequence, sizeof(consequence)) ||
			    !duckvep_scalar_format_region(result->region_mask, region,
			    sizeof(region))) {
				duckvep_sql_set_error(error, error_size,
				    "duckvep_annotate: could not render a kernel result");
				return 0;
			}
			duckdb_vector_assign_string_element(vectors[2],
			    (idx_t)output_count, consequence);
			duckdb_vector_assign_string_element(vectors[3],
			    (idx_t)output_count, duckvep_impact_name(
			    (duckvep_impact_t)result->impact));
			duckdb_vector_assign_string_element(vectors[4],
			    (idx_t)output_count, region);
			status = (result->flags &
			    DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0 ?
			    "unresolved" : "supported";
			duckdb_vector_assign_string_element(vectors[5],
			    (idx_t)output_count, status);
			if (status[0] == 'u') {
				reason = duckvep_scalar_sequence_reason(
				    result->sequence_status);
				duckdb_vector_assign_string_element(vectors[6],
				    (idx_t)output_count, reason);
			} else {
				duckvep_scalar_set_null(vectors[6],
				    (idx_t)output_count);
			}
			if (result->cdna_pos >= 0)
				cdna_positions[output_count] =
				    (uint64_t)result->cdna_pos;
			else
				duckvep_scalar_set_null(vectors[7],
				    (idx_t)output_count);
			if (result->cds_pos >= 0)
				cds_positions[output_count] =
				    (uint64_t)result->cds_pos;
			else
				duckvep_scalar_set_null(vectors[8],
				    (idx_t)output_count);
			if (result->protein_pos >= 0)
				protein_positions[output_count] =
				    (uint64_t)result->protein_pos;
			else
				duckvep_scalar_set_null(vectors[9],
				    (idx_t)output_count);
			amino_acid[1] = '\0';
			if (result->aa_ref != 0u) {
				amino_acid[0] = (char)result->aa_ref;
				duckdb_vector_assign_string_element(vectors[10],
				    (idx_t)output_count, amino_acid);
			} else {
				duckvep_scalar_set_null(vectors[10],
				    (idx_t)output_count);
			}
			if (result->aa_alt != 0u) {
				amino_acid[0] = (char)result->aa_alt;
				duckdb_vector_assign_string_element(vectors[11],
				    (idx_t)output_count, amino_acid);
			} else {
				duckvep_scalar_set_null(vectors[11],
				    (idx_t)output_count);
			}
		}
	}
	return 1;
}

static void
duckvep_annotate_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_scalar_state_t *state;
	duckvep_variant_batch_t batch;
	duckdb_vector model_vector, distance_vector;
	uint16_t *regions;
	uint32_t *positions;
	uint64_t *distances;
	idx_t rows;
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
	if (!duckvep_scalar_prepare_batch(state, input, &batch, error,
	    sizeof(error)))
		goto failed;
	regions = state->seq_regions;
	positions = state->positions;
	distance_vector = duckdb_data_chunk_get_column_count(input) == 6 ?
	    duckdb_data_chunk_get_vector(input, 5) : NULL;
	distances = distance_vector != NULL ?
	    duckdb_vector_get_data(distance_vector) : NULL;
	state->result_count = 0;
	begin = 0;
	while (begin < (size_t)rows) {
		uint64_t distance;
		size_t end;

		if (!duckvep_scalar_select_model(state, model_vector, (idx_t)begin,
		    error, sizeof(error)))
			goto failed;
		if (positions[begin] == 0 ||
		    (distance_vector != NULL &&
		    duckvep_row_is_null(distance_vector, (idx_t)begin))) {
			duckvep_sql_set_error(error, sizeof(error),
			    "duckvep_annotate: position and distance must be non-NULL; position is one-based");
			goto failed;
		}
		distance = distances != NULL ? distances[begin] : 5000;
		end = begin + 1;
		while (end < (size_t)rows) {
			uint64_t next_distance;

			if (distance_vector != NULL &&
			    duckvep_row_is_null(distance_vector, (idx_t)end))
				break;
			next_distance = distances != NULL ? distances[end] : 5000;
			if (!duckvep_vector_string_equals(model_vector, (idx_t)end,
			    state->entry->name) ||
			    regions[end] != regions[begin] ||
			    positions[end] < positions[end - 1] ||
			    next_distance != distance)
				break;
			end++;
		}
		if (!duckvep_scalar_run(state, &batch, begin, end - begin,
		    distance, error, sizeof(error)))
			goto failed;
		begin = end;
	}
	if (!duckvep_scalar_write_output(state, output, rows, error,
	    sizeof(error)))
		goto failed;
	duckvep_scalar_state_release(state);
	return;

failed:
	duckvep_scalar_state_release(state);
	duckdb_scalar_function_set_error(info,
	    error[0] != '\0' ? error : "duckvep_annotate failed");
}

static duckdb_logical_type
duckvep_annotation_list_type(void)
{
	const char *names[] = {
		"transcript_index", "gene_index", "consequence", "impact",
		"region", "status", "reason", "cdna_position", "cds_position",
		"protein_position", "reference_amino_acid",
		"alternate_amino_acid"
	};
	duckdb_logical_type types[12], structure, list;
	size_t index;

	types[0] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	types[1] = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	for (index = 2; index <= 6; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	for (index = 7; index <= 9; index++)
		types[index] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	types[10] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	types[11] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	structure = duckdb_create_struct_type(types, names, 12);
	list = duckdb_create_list_type(structure);
	for (index = 0; index < 12; index++)
		duckdb_destroy_logical_type(&types[index]);
	duckdb_destroy_logical_type(&structure);
	return list;
}

static void
duckvep_register_annotate_scalar(duckdb_connection connection,
	duckvep_registry_t *registry, duckdb_logical_type varchar_type,
	duckdb_logical_type uinteger_type, duckdb_logical_type ubigint_type,
	int with_distance)
{
	duckdb_scalar_function scalar;
	duckdb_logical_type result_type;

	scalar = duckdb_create_scalar_function();
	result_type = duckvep_annotation_list_type();
	duckdb_scalar_function_set_name(scalar, "duckvep_annotate");
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, uinteger_type);
	duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	if (with_distance)
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
	duckdb_scalar_function_set_return_type(scalar, result_type);
	duckdb_scalar_function_set_volatile(scalar);
	duckvep_registry_retain(registry);
	duckdb_scalar_function_set_extra_info(scalar, registry,
	    duckvep_registry_release);
	duckdb_scalar_function_set_function(scalar, duckvep_annotate_scalar);
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
	    uinteger_type, ubigint_type, 0);
	duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1);
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&uinteger_type);
	duckdb_destroy_logical_type(&ubigint_type);

	duckvep_registry_release(registry);
}
