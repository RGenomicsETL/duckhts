#include "duckvep_model.h"

DUCKDB_EXTENSION_EXTERN

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct duckvep_query_result {
	duckdb_prepared_statement statement;
	duckdb_result result;
	int have_result;
} duckvep_query_result_t;

typedef struct duckvep_load_bind {
	duckvep_registry_t *registry;
	char *arguments[4];
	int transcript_coverage_complete;
} duckvep_load_bind_t;

typedef struct duckvep_load_state {
	int emitted;
} duckvep_load_state_t;

void
duckvep_sql_set_error(char *error, size_t error_size, const char *message)
{
	if (error != NULL && error_size != 0)
		(void)snprintf(error, error_size, "%s",
		    message != NULL ? message : "unknown error");
}

int
duckvep_sql_resize(void **pointer, size_t width, size_t count)
{
	void *resized;

	if (width != 0 && count > SIZE_MAX / width)
		return 0;
	resized = realloc(*pointer, width * count);
	if (resized == NULL && count != 0)
		return 0;
	*pointer = resized;
	return 1;
}

size_t
duckvep_sql_next_capacity(size_t current, size_t needed)
{
	size_t capacity;

	capacity = current != 0 ? current : 64;
	while (capacity < needed) {
		if (capacity > SIZE_MAX / 2)
			return needed;
		capacity *= 2;
	}
	return capacity;
}

static int
duckvep_model_reserve_regions(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->known_seq_region_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(
	    model->known_seq_region_capacity, needed);
	if (!duckvep_sql_resize((void **)&model->known_seq_regions,
	    sizeof(*model->known_seq_regions), capacity))
		return 0;
	if (!duckvep_sql_resize((void **)&model->sequence_lengths,
	    sizeof(*model->sequence_lengths), capacity))
		return 0;
	model->known_seq_region_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_transcripts(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity, old_capacity;

	if (needed <= model->transcript_capacity)
		return 1;
	old_capacity = model->transcript_capacity;
	capacity = duckvep_sql_next_capacity(model->transcript_capacity,
	    needed);
#define DUCKVEP_RESIZE_TRANSCRIPT(member) \
	if (!duckvep_sql_resize((void **)&model->member, \
	    sizeof(*model->member), capacity)) \
		return 0
	DUCKVEP_RESIZE_TRANSCRIPT(seq_regions);
	DUCKVEP_RESIZE_TRANSCRIPT(transcript_starts);
	DUCKVEP_RESIZE_TRANSCRIPT(transcript_ends);
	DUCKVEP_RESIZE_TRANSCRIPT(strands);
	DUCKVEP_RESIZE_TRANSCRIPT(transcript_flags);
	DUCKVEP_RESIZE_TRANSCRIPT(gene_indices);
	DUCKVEP_RESIZE_TRANSCRIPT(exon_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(exon_counts);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_starts);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_ends);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_sequence_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_sequence_lengths);
	DUCKVEP_RESIZE_TRANSCRIPT(codon_tables);
#undef DUCKVEP_RESIZE_TRANSCRIPT
	if (capacity > SIZE_MAX / DUCKVEP_POST_CDS_BASE_COUNT ||
	    !duckvep_sql_resize((void **)&model->post_cds_bases,
	    sizeof(*model->post_cds_bases),
	    capacity * DUCKVEP_POST_CDS_BASE_COUNT))
		return 0;
	memset(model->post_cds_bases +
	    old_capacity * DUCKVEP_POST_CDS_BASE_COUNT, 0,
	    (capacity - old_capacity) * DUCKVEP_POST_CDS_BASE_COUNT);
	model->transcript_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_exons(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->exon_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->exon_capacity, needed);
#define DUCKVEP_RESIZE_EXON(member) \
	if (!duckvep_sql_resize((void **)&model->member, \
	    sizeof(*model->member), capacity)) \
		return 0
	DUCKVEP_RESIZE_EXON(exon_starts);
	DUCKVEP_RESIZE_EXON(exon_ends);
	DUCKVEP_RESIZE_EXON(exon_cdna_starts);
	DUCKVEP_RESIZE_EXON(exon_cdna_ends);
	DUCKVEP_RESIZE_EXON(exon_phases);
	DUCKVEP_RESIZE_EXON(exon_end_phases);
#undef DUCKVEP_RESIZE_EXON
	model->exon_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_sequence(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->cds_sequence_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->cds_sequence_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->cds_sequence_bytes,
	    sizeof(*model->cds_sequence_bytes), capacity))
		return 0;
	model->cds_sequence_capacity = capacity;
	return 1;
}

static void
duckvep_owned_model_destroy(duckvep_owned_model_t *model)
{
	if (model == NULL)
		return;
	if (model->kernel != NULL)
		duckvep_model_close(model->kernel);
	free(model->known_seq_regions);
	free(model->sequence_lengths);
	free(model->seq_regions);
	free(model->transcript_starts);
	free(model->transcript_ends);
	free(model->strands);
	free(model->transcript_flags);
	free(model->gene_indices);
	free(model->exon_offsets);
	free(model->exon_counts);
	free(model->cds_starts);
	free(model->cds_ends);
	free(model->cds_sequence_offsets);
	free(model->cds_sequence_lengths);
	free(model->codon_tables);
	free(model->post_cds_bases);
	free(model->exon_starts);
	free(model->exon_ends);
	free(model->exon_cdna_starts);
	free(model->exon_cdna_ends);
	free(model->exon_phases);
	free(model->exon_end_phases);
	free(model->cds_sequence_bytes);
	if (model->interval_index != NULL) {
		/* cgranges 0.1.1 does not release its interval array. */
		free(model->interval_index->r);
		model->interval_index->r = NULL;
		cr_destroy(model->interval_index);
	}
	memset(model, 0, sizeof(*model));
}

static void
duckvep_owned_model_publish(duckvep_owned_model_t *model)
{
	model->transcripts.chrom_id = model->seq_regions;
	model->transcripts.start1 = model->transcript_starts;
	model->transcripts.end1 = model->transcript_ends;
	model->transcripts.strand = model->strands;
	model->transcripts.flags = model->transcript_flags;
	model->transcripts.exon_offset = model->exon_offsets;
	model->transcripts.exon_count = model->exon_counts;
	model->transcripts.cds_start1 = model->cds_starts;
	model->transcripts.cds_end1 = model->cds_ends;
	model->exons.start1 = model->exon_starts;
	model->exons.end1 = model->exon_ends;
	model->exons.cdna_start1 = model->exon_cdna_starts;
	model->exons.cdna_end1 = model->exon_cdna_ends;
	model->exons.phase = model->exon_phases;
	model->exons.end_phase = model->exon_end_phases;
	model->sequences.cds_bytes = model->cds_sequence_bytes;
	model->sequences.cds_bytes_len = model->cds_sequence_length;
	model->sequences.cds_offset = model->cds_sequence_offsets;
	model->sequences.cds_length = model->cds_sequence_lengths;
	model->sequences.codon_table = model->codon_tables;
	model->sequences.post_cds_bases = model->post_cds_bases;
	model->sequences.transcript_count = model->transcripts.transcript_count;
}

static int
duckvep_owned_model_index(duckvep_owned_model_t *model, char *error,
	size_t error_size)
{
	size_t index;
	char region_name[16];

	model->interval_index_complete = 0;
	if (model->transcripts.transcript_count > (size_t)INT32_MAX)
		return 1;
	for (index = 0; index < model->transcripts.transcript_count; index++) {
		if (model->transcript_starts[index] > (uint32_t)INT32_MAX ||
		    model->transcript_ends[index] > (uint32_t)INT32_MAX)
			return 1;
	}
	model->interval_index = cr_init();
	if (model->interval_index == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory building transcript interval index");
		return 0;
	}
	for (index = 0; index < model->known_seq_region_count; index++) {
		(void)snprintf(region_name, sizeof(region_name), "%u",
		    model->known_seq_regions[index]);
		(void)cr_add_ctg(model->interval_index, region_name, 0);
	}
	for (index = 0; index < model->transcripts.transcript_count; index++) {
		(void)snprintf(region_name, sizeof(region_name), "%u",
		    model->seq_regions[index]);
		if (cr_add(model->interval_index, region_name,
		    (int32_t)(model->transcript_starts[index] - 1),
		    (int32_t)model->transcript_ends[index],
		    (int32_t)index) == NULL) {
			duckvep_sql_set_error(error, error_size,
			    "could not build transcript interval index");
			return 0;
		}
	}
	cr_index(model->interval_index);
	model->interval_index_complete = 1;
	return 1;
}

static void
duckvep_query_result_close(duckvep_query_result_t *query_result)
{
	if (query_result == NULL)
		return;
	if (query_result->have_result)
		duckdb_destroy_result(&query_result->result);
	if (query_result->statement != NULL)
		duckdb_destroy_prepare(&query_result->statement);
	memset(query_result, 0, sizeof(*query_result));
}

static int
duckvep_query_result_open(duckdb_connection connection, const char *query,
	duckvep_query_result_t *query_result, char *error, size_t error_size)
{
	duckdb_state state;
	const char *message;

	memset(query_result, 0, sizeof(*query_result));
	state = duckdb_prepare(connection, query, &query_result->statement);
	if (state != DuckDBSuccess) {
		message = duckdb_prepare_error(query_result->statement);
		duckvep_sql_set_error(error, error_size, message);
		duckvep_query_result_close(query_result);
		return 0;
	}
	state = duckdb_execute_prepared(query_result->statement,
	    &query_result->result);
	query_result->have_result = 1;
	if (state != DuckDBSuccess) {
		duckvep_sql_set_error(error, error_size,
		    duckdb_result_error(&query_result->result));
		duckvep_query_result_close(query_result);
		return 0;
	}
	return 1;
}

int
duckvep_row_is_null(duckdb_vector vector, idx_t row)
{
	uint64_t *validity;

	validity = duckdb_vector_get_validity(vector);
	return validity != NULL &&
	    ((validity[row / 64] >> (row % 64)) & UINT64_C(1)) == 0;
}

static int
duckvep_result_schema(duckdb_result *result, const char *const *names,
	const duckdb_type *types, size_t count, size_t flexible_string_column,
	char *error, size_t error_size)
{
	size_t column;

	if (duckdb_column_count(result) != (idx_t)count) {
		duckvep_sql_set_error(error, error_size,
		    "query returned the wrong number of columns");
		return 0;
	}
	for (column = 0; column < count; column++) {
		const char *name;
		duckdb_type type;

		name = duckdb_column_name(result, (idx_t)column);
		type = duckdb_column_type(result, (idx_t)column);
		if (name == NULL || strcmp(name, names[column]) != 0) {
			(void)snprintf(error, error_size,
			    "query column %zu must be named %s", column + 1,
			    names[column]);
			return 0;
		}
		if (column == flexible_string_column &&
		    (type == DUCKDB_TYPE_VARCHAR || type == DUCKDB_TYPE_BLOB))
			continue;
		if (type != types[column]) {
			(void)snprintf(error, error_size,
			    "query column %s has the wrong type", names[column]);
			return 0;
		}
	}
	return 1;
}

static int
duckvep_load_regions(duckdb_connection connection, const char *query,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const one_name[] = {"seq_region"};
	static const duckdb_type one_type[] = {DUCKDB_TYPE_UINTEGER};
	static const char *const two_names[] = {"seq_region", "sequence_length"};
	static const duckdb_type two_types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT
	};
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	idx_t column_count;
	int ok;

	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	ok = 0;
	column_count = duckdb_column_count(&query_result.result);
	if (column_count != 1 && column_count != 2) {
		duckvep_sql_set_error(error, error_size,
		    "seq_region query must return seq_region and optional sequence_length");
		goto done;
	}
	if (model->transcript_coverage_complete && column_count != 2) {
		duckvep_sql_set_error(error, error_size,
		    "complete transcript coverage requires sequence_length for every region");
		goto done;
	}
	if (!duckvep_result_schema(&query_result.result,
	    column_count == 1 ? one_name : two_names,
	    column_count == 1 ? one_type : two_types,
	    (size_t)column_count, SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector region_vector, length_vector;
		uint32_t *values;
		uint64_t *lengths;
		idx_t row, rows;

		rows = duckdb_data_chunk_get_size(chunk);
		region_vector = duckdb_data_chunk_get_vector(chunk, 0);
		length_vector = column_count == 2
		    ? duckdb_data_chunk_get_vector(chunk, 1) : NULL;
		values = (uint32_t *)duckdb_vector_get_data(region_vector);
		lengths = length_vector != NULL
		    ? (uint64_t *)duckdb_vector_get_data(length_vector) : NULL;
		for (row = 0; row < rows; row++) {
			size_t index;

			if (duckvep_row_is_null(region_vector, row) ||
			    (length_vector != NULL &&
			    duckvep_row_is_null(length_vector, row))) {
				duckvep_sql_set_error(error, error_size,
				    "seq_region query contains NULL");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			index = model->known_seq_region_count;
			if (values[row] > UINT16_MAX) {
				duckvep_sql_set_error(error, error_size,
				    "seq_region exceeds the compact uint16 model limit");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (lengths != NULL &&
			    (lengths[row] == 0 || lengths[row] > UINT32_MAX)) {
				duckvep_sql_set_error(error, error_size,
				    "sequence_length must fit a positive uint32 coordinate");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (index != 0 && model->known_seq_regions[index - 1] >=
			    values[row]) {
				duckvep_sql_set_error(error, error_size,
				    "seq_region query must be sorted and unique");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (!duckvep_model_reserve_regions(model, index + 1)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading sequence regions");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			model->known_seq_regions[index] = (uint16_t)values[row];
			model->sequence_lengths[index] = lengths != NULL
			    ? (uint32_t)lengths[row] : 0;
			model->known_seq_region_count++;
		}
		duckdb_destroy_data_chunk(&chunk);
	}
	if (model->known_seq_region_count == 0) {
		duckvep_sql_set_error(error, error_size,
		    "seq_region query returned no rows");
		goto done;
	}
	ok = 1;
done:
	duckvep_query_result_close(&query_result);
	return ok;
}

static int
duckvep_load_transcripts(duckdb_connection connection, const char *query,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "seq_region", "transcript_start",
		"transcript_end", "strand", "gene_index", "transcript_flags",
		"cds_start", "cds_end", "cds_sequence", "codon_table",
		"post_cds_bases"
	};
	static const duckdb_type types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
		DUCKDB_TYPE_TINYINT, DUCKDB_TYPE_UINTEGER,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_BLOB,
		DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_BLOB
	};
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	idx_t column_count;
	int ok;

	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	ok = 0;
	column_count = duckdb_column_count(&query_result.result);
	if (column_count != 11 && column_count != 12) {
		duckvep_sql_set_error(error, error_size,
		    "transcript query must return 11 columns and optional post_cds_bases");
		goto done;
	}
	if (!duckvep_result_schema(&query_result.result, names, types,
	    (size_t)column_count, 9,
	    error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector vectors[12];
		uint32_t *transcript_indices, *seq_regions, *gene_indices;
		uint64_t *starts, *ends, *flags, *cds_starts, *cds_ends;
		int8_t *strands;
		duckdb_string_t *sequences, *tails;
		uint8_t *tables;
		idx_t row, rows;
		size_t column;

		rows = duckdb_data_chunk_get_size(chunk);
		for (column = 0; column < (size_t)column_count; column++)
			vectors[column] = duckdb_data_chunk_get_vector(chunk,
			    (idx_t)column);
		transcript_indices = duckdb_vector_get_data(vectors[0]);
		seq_regions = duckdb_vector_get_data(vectors[1]);
		starts = duckdb_vector_get_data(vectors[2]);
		ends = duckdb_vector_get_data(vectors[3]);
		strands = duckdb_vector_get_data(vectors[4]);
		gene_indices = duckdb_vector_get_data(vectors[5]);
		flags = duckdb_vector_get_data(vectors[6]);
		cds_starts = duckdb_vector_get_data(vectors[7]);
		cds_ends = duckdb_vector_get_data(vectors[8]);
		sequences = duckdb_vector_get_data(vectors[9]);
		tables = duckdb_vector_get_data(vectors[10]);
		tails = column_count == 12 ? duckdb_vector_get_data(vectors[11]) : NULL;
		for (row = 0; row < rows; row++) {
			size_t index, sequence_length, sequence_offset, tail_length;
			int cds_nulls, sequence_nulls;

			for (column = 0; column < 7; column++) {
				if (duckvep_row_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "transcript query contains NULL in %s",
					    names[column]);
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			index = model->transcripts.transcript_count;
			if (transcript_indices[row] != (uint32_t)index ||
			    index > UINT32_MAX) {
				duckvep_sql_set_error(error, error_size,
				    "transcript_index must be dense, zero-based, and ordered");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (!duckvep_model_reserve_transcripts(model, index + 1)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading transcripts");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (seq_regions[row] > UINT16_MAX || starts[row] == 0 ||
			    starts[row] > UINT32_MAX || ends[row] > UINT32_MAX ||
			    ends[row] < starts[row] ||
			    (strands[row] != 1 && strands[row] != -1)) {
				duckvep_sql_set_error(error, error_size,
				    "transcript row has an invalid region, span, or strand");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			model->seq_regions[index] = (uint16_t)seq_regions[row];
			model->transcript_starts[index] = (uint32_t)starts[row];
			model->transcript_ends[index] = (uint32_t)ends[row];
			model->strands[index] = strands[row];
			model->transcript_flags[index] = flags[row];
			model->gene_indices[index] = gene_indices[row];
			model->exon_offsets[index] = 0;
			model->exon_counts[index] = 0;
			cds_nulls = duckvep_row_is_null(vectors[7], row) +
			    duckvep_row_is_null(vectors[8], row);
			sequence_nulls = duckvep_row_is_null(vectors[9], row) +
			    duckvep_row_is_null(vectors[10], row);
			if (cds_nulls == 1 || sequence_nulls == 1 ||
			    (cds_nulls == 2 && sequence_nulls == 0)) {
				duckvep_sql_set_error(error, error_size,
				    "CDS span and CDS sequence/table must each be both NULL or both present");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			sequence_offset = model->cds_sequence_length;
			sequence_length = 0;
			if (cds_nulls == 0) {
				if (cds_starts[row] == 0 ||
				    cds_starts[row] > UINT32_MAX ||
				    cds_ends[row] > UINT32_MAX ||
				    cds_ends[row] < cds_starts[row] ||
				    cds_starts[row] < starts[row] ||
				    cds_ends[row] > ends[row]) {
					duckvep_sql_set_error(error, error_size,
					    "coding transcript has an invalid CDS span");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
				model->cds_starts[index] = (uint32_t)cds_starts[row];
				model->cds_ends[index] = (uint32_t)cds_ends[row];
			} else {
				model->cds_starts[index] = 0;
				model->cds_ends[index] = 0;
			}
			model->codon_tables[index] = 1;
			tail_length = 0;
			memset(model->post_cds_bases +
			    index * DUCKVEP_POST_CDS_BASE_COUNT, 0,
			    DUCKVEP_POST_CDS_BASE_COUNT);
			if (sequence_nulls == 0) {
				sequence_length = (size_t)duckdb_string_t_length(
				    sequences[row]);
				if (sequence_length == 0 || sequence_length > UINT32_MAX ||
				    (tables[row] != 1 && tables[row] != 2)) {
					duckvep_sql_set_error(error, error_size,
					    "coding transcript has an invalid CDS sequence or codon table");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
				if (sequence_length > SIZE_MAX - sequence_offset ||
				    !duckvep_model_reserve_sequence(model,
				    sequence_offset + sequence_length)) {
					duckvep_sql_set_error(error, error_size,
					    "out of memory loading coding sequences");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
				memcpy(model->cds_sequence_bytes + sequence_offset,
				    duckdb_string_t_data(&sequences[row]),
				    sequence_length);
				model->codon_tables[index] = tables[row];
			}
			if (tails != NULL && !duckvep_row_is_null(vectors[11], row)) {
				tail_length = (size_t)duckdb_string_t_length(tails[row]);
				if (tail_length > DUCKVEP_POST_CDS_BASE_COUNT ||
				    (tail_length != 0 && cds_nulls != 0)) {
					duckvep_sql_set_error(error, error_size,
					    "post_cds_bases must contain at most three bases for a coding transcript");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
				if (tail_length != 0)
					memcpy(model->post_cds_bases +
					    index * DUCKVEP_POST_CDS_BASE_COUNT,
					    duckdb_string_t_data(&tails[row]),
					    tail_length);
			}
			model->cds_sequence_offsets[index] =
			    (uint64_t)sequence_offset;
			model->cds_sequence_lengths[index] =
			    (uint32_t)sequence_length;
			model->cds_sequence_length += sequence_length;
			model->transcripts.transcript_count++;
		}
		duckdb_destroy_data_chunk(&chunk);
	}
	if (model->transcripts.transcript_count == 0) {
		duckvep_sql_set_error(error, error_size,
		    "transcript query returned no rows");
		goto done;
	}
	ok = 1;
done:
	duckvep_query_result_close(&query_result);
	return ok;
}

static int
duckvep_load_exons(duckdb_connection connection, const char *query,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "exon_start", "exon_end",
		"exon_cdna_start", "exon_cdna_end", "phase", "end_phase"
	};
	static const duckdb_type types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_TINYINT,
		DUCKDB_TYPE_TINYINT
	};
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	uint32_t previous_transcript;
	int have_previous, ok;

	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	have_previous = 0;
	previous_transcript = 0;
	ok = 0;
	if (!duckvep_result_schema(&query_result.result, names, types, 7, SIZE_MAX,
	    error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector vectors[7];
		uint32_t *transcript_indices;
		uint64_t *starts, *ends, *cdna_starts, *cdna_ends;
		int8_t *phases, *end_phases;
		idx_t row, rows;
		size_t column;

		rows = duckdb_data_chunk_get_size(chunk);
		for (column = 0; column < 7; column++)
			vectors[column] = duckdb_data_chunk_get_vector(chunk,
			    (idx_t)column);
		transcript_indices = duckdb_vector_get_data(vectors[0]);
		starts = duckdb_vector_get_data(vectors[1]);
		ends = duckdb_vector_get_data(vectors[2]);
		cdna_starts = duckdb_vector_get_data(vectors[3]);
		cdna_ends = duckdb_vector_get_data(vectors[4]);
		phases = duckdb_vector_get_data(vectors[5]);
		end_phases = duckdb_vector_get_data(vectors[6]);
		for (row = 0; row < rows; row++) {
			uint32_t transcript_index;
			size_t exon_index;

			for (column = 0; column < 7; column++) {
				if (duckvep_row_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "exon query contains NULL in %s",
					    names[column]);
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			transcript_index = transcript_indices[row];
			if (transcript_index >= model->transcripts.transcript_count ||
			    (have_previous && transcript_index < previous_transcript)) {
				duckvep_sql_set_error(error, error_size,
				    "exon query must be grouped by transcript_index");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (model->exon_counts[transcript_index] == UINT16_MAX) {
				duckvep_sql_set_error(error, error_size,
				    "transcript has too many exons");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			exon_index = model->exons.exon_count;
			if (exon_index > UINT32_MAX || starts[row] == 0 ||
			    starts[row] > UINT32_MAX || ends[row] > UINT32_MAX ||
			    ends[row] < starts[row] || cdna_starts[row] == 0 ||
			    cdna_starts[row] > UINT32_MAX ||
			    cdna_ends[row] > UINT32_MAX ||
			    cdna_ends[row] < cdna_starts[row] ||
			    starts[row] < model->transcript_starts[transcript_index] ||
			    ends[row] > model->transcript_ends[transcript_index] ||
			    phases[row] < -1 || phases[row] > 2 ||
			    end_phases[row] < -1 || end_phases[row] > 2) {
				duckvep_sql_set_error(error, error_size,
				    "exon row has invalid coordinates or phase");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (model->exon_counts[transcript_index] != 0) {
				size_t previous_exon;
				int ordered;

				previous_exon = exon_index - 1;
				ordered = cdna_starts[row] >
				    model->exon_cdna_ends[previous_exon];
				if (model->strands[transcript_index] > 0)
					ordered = ordered && starts[row] >
					    model->exon_ends[previous_exon];
				else
					ordered = ordered && ends[row] <
					    model->exon_starts[previous_exon];
				if (!ordered) {
					duckvep_sql_set_error(error, error_size,
					    "exons must be non-overlapping and ordered on the transcript");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			if (!duckvep_model_reserve_exons(model, exon_index + 1)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading exons");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (model->exon_counts[transcript_index] == 0)
				model->exon_offsets[transcript_index] =
				    (uint32_t)exon_index;
			model->exon_starts[exon_index] = (uint32_t)starts[row];
			model->exon_ends[exon_index] = (uint32_t)ends[row];
			model->exon_cdna_starts[exon_index] =
			    (uint32_t)cdna_starts[row];
			model->exon_cdna_ends[exon_index] =
			    (uint32_t)cdna_ends[row];
			model->exon_phases[exon_index] = phases[row];
			model->exon_end_phases[exon_index] = end_phases[row];
			model->exon_counts[transcript_index]++;
			model->exons.exon_count++;
			previous_transcript = transcript_index;
			have_previous = 1;
		}
		duckdb_destroy_data_chunk(&chunk);
	}
	{
		size_t transcript;

		for (transcript = 0; transcript < model->transcripts.transcript_count;
		    transcript++) {
			if (model->exon_counts[transcript] == 0) {
				duckvep_sql_set_error(error, error_size,
				    "every transcript must have at least one exon");
				goto done;
			}
		}
	}
	ok = 1;
done:
	duckvep_query_result_close(&query_result);
	return ok;
}

static int
duckvep_query_command(duckdb_connection connection, const char *sql,
	char *error, size_t error_size)
{
	duckdb_result result;
	duckdb_state state;

	memset(&result, 0, sizeof(result));
	state = duckdb_query(connection, sql, &result);
	if (state != DuckDBSuccess)
		duckvep_sql_set_error(error, error_size,
		    duckdb_result_error(&result));
	duckdb_destroy_result(&result);
	return state == DuckDBSuccess;
}

static int
duckvep_model_load_queries(duckdb_connection connection,
	const char *region_query, const char *transcript_query,
	const char *exon_query, int transcript_coverage_complete,
	duckvep_owned_model_t *model,
	char *error, size_t error_size)
{
	duckvep_error_t kernel_error;
	size_t transcript;
	int ok;

	memset(model, 0, sizeof(*model));
	if (connection == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "model-loading connection is unavailable");
		return 0;
	}
	model->transcript_coverage_complete = transcript_coverage_complete;
	if (!duckvep_query_command(connection, "BEGIN TRANSACTION", error,
	    error_size))
		return 0;
	ok = duckvep_load_regions(connection, region_query, model, error,
	    error_size) &&
	    duckvep_load_transcripts(connection, transcript_query, model, error,
	    error_size) &&
	    duckvep_load_exons(connection, exon_query, model, error, error_size);
	if (ok)
		ok = duckvep_query_command(connection, "COMMIT", error, error_size);
	if (!ok)
		(void)duckvep_query_command(connection, "ROLLBACK", NULL, 0);
	if (!ok) {
		duckvep_owned_model_destroy(model);
		return 0;
	}
	duckvep_owned_model_publish(model);
	for (transcript = 0;
	    transcript < model->transcripts.transcript_count; transcript++) {
		size_t region;
		int found;

		found = 0;
		for (region = 0; region < model->known_seq_region_count; region++) {
			if (model->known_seq_regions[region] ==
			    model->seq_regions[transcript]) {
				found = 1;
				break;
			}
		}
		if (!found) {
			duckvep_sql_set_error(error, error_size,
			    "transcript seq_region is absent from the region query");
			duckvep_owned_model_destroy(model);
			return 0;
		}
	}
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (duckvep_model_open(&model->transcripts, &model->exons,
	    &model->sequences, &model->kernel, &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "invalid transcript model: %s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "kernel validation failed");
		duckvep_owned_model_destroy(model);
		return 0;
	}
	if (!duckvep_owned_model_index(model, error, error_size)) {
		duckvep_owned_model_destroy(model);
		return 0;
	}
	return 1;
}

static void
duckvep_model_entry_destroy(duckvep_model_entry_t *entry)
{
	duckvep_workspace_cache_t *workspace, *next;

	if (entry == NULL)
		return;
	for (workspace = entry->workspaces; workspace != NULL;
	    workspace = next) {
		next = workspace->next;
		duckvep_workspace_close(workspace->workspace);
		free(workspace);
	}
	free(entry->name);
	duckvep_owned_model_destroy(&entry->model);
	free(entry);
}

static duckvep_model_entry_t *
duckvep_registry_find_locked(duckvep_registry_t *registry, const char *name)
{
	duckvep_model_entry_t *entry;

	for (entry = registry->models; entry != NULL; entry = entry->next) {
		if (strcmp(entry->name, name) == 0)
			return entry;
	}
	return NULL;
}

duckvep_model_entry_t *
duckvep_registry_pin(duckvep_registry_t *registry, const char *name)
{
	duckvep_model_entry_t *entry;

	pthread_mutex_lock(&registry->mutex);
	entry = duckvep_registry_find_locked(registry, name);
	if (entry != NULL)
		entry->pins++;
	pthread_mutex_unlock(&registry->mutex);
	return entry;
}

void
duckvep_registry_unpin(duckvep_registry_t *registry,
	duckvep_model_entry_t *entry)
{
	if (registry == NULL || entry == NULL)
		return;
	pthread_mutex_lock(&registry->mutex);
	if (entry->pins != 0)
		entry->pins--;
	pthread_mutex_unlock(&registry->mutex);
}

duckvep_workspace_cache_t *
duckvep_registry_workspace_take(duckvep_registry_t *registry,
	duckvep_model_entry_t *entry, char *error, size_t error_size)
{
	duckvep_workspace_cache_t *cache;
	duckvep_error_t kernel_error;

	if (registry == NULL || entry == NULL)
		return NULL;
	pthread_mutex_lock(&registry->mutex);
	cache = entry->workspaces;
	if (cache != NULL) {
		entry->workspaces = cache->next;
		cache->next = NULL;
	}
	pthread_mutex_unlock(&registry->mutex);
	if (cache != NULL)
		return cache;
	cache = calloc(1, sizeof(*cache));
	if (cache == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory allocating a DuckVEP workspace");
		return NULL;
	}
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (duckvep_workspace_open(entry->model.kernel, &cache->workspace,
	    &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "%s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "could not open a DuckVEP workspace");
		free(cache);
		return NULL;
	}
	return cache;
}

void
duckvep_registry_workspace_return(duckvep_registry_t *registry,
	duckvep_model_entry_t *entry, duckvep_workspace_cache_t *cache)
{
	if (registry == NULL || entry == NULL || cache == NULL)
		return;
	pthread_mutex_lock(&registry->mutex);
	cache->next = entry->workspaces;
	entry->workspaces = cache;
	pthread_mutex_unlock(&registry->mutex);
}

void
duckvep_registry_retain(duckvep_registry_t *registry)
{
	pthread_mutex_lock(&registry->mutex);
	registry->references++;
	pthread_mutex_unlock(&registry->mutex);
}

void
duckvep_registry_release(void *pointer)
{
	duckvep_registry_t *registry;
	duckvep_model_entry_t *entry, *next;
	int destroy;

	registry = pointer;
	if (registry == NULL)
		return;
	pthread_mutex_lock(&registry->mutex);
	if (registry->references != 0)
		registry->references--;
	destroy = registry->references == 0;
	pthread_mutex_unlock(&registry->mutex);
	if (!destroy)
		return;
	for (entry = registry->models; entry != NULL; entry = next) {
		next = entry->next;
		duckvep_model_entry_destroy(entry);
	}
	if (registry->annotation_state_pool_destroy != NULL)
		registry->annotation_state_pool_destroy(
		    registry->annotation_state_pool);
	if (registry->query_connection != NULL)
		duckdb_disconnect(&registry->query_connection);
	pthread_mutex_destroy(&registry->query_mutex);
	pthread_mutex_destroy(&registry->mutex);
	free(registry);
}

char *
duckvep_vector_string(duckdb_vector vector, idx_t row)
{
	duckdb_string_t *strings;
	const char *data;
	uint32_t length;
	char *copy;

	if (duckvep_row_is_null(vector, row))
		return NULL;
	strings = duckdb_vector_get_data(vector);
	length = duckdb_string_t_length(strings[row]);
	data = duckdb_string_t_data(&strings[row]);
	if (memchr(data, '\0', length) != NULL)
		return NULL;
	copy = malloc((size_t)length + 1);
	if (copy == NULL)
		return NULL;
	memcpy(copy, data, length);
	copy[length] = '\0';
	return copy;
}

static char *
duckvep_bind_string(duckdb_bind_info info, idx_t parameter)
{
	duckdb_value value;
	char *string;

	value = duckdb_bind_get_parameter(info, parameter);
	if (value == NULL || duckdb_is_null_value(value)) {
		if (value != NULL)
			duckdb_destroy_value(&value);
		return NULL;
	}
	string = duckdb_get_varchar(value);
	duckdb_destroy_value(&value);
	return string;
}

static char *
duckvep_string_copy(const char *source)
{
	size_t length;
	char *copy;

	if (source == NULL)
		return NULL;
	length = strlen(source);
	copy = malloc(length + 1);
	if (copy != NULL)
		memcpy(copy, source, length + 1);
	return copy;
}

static void
duckvep_model_load_bind_destroy(void *pointer)
{
	duckvep_load_bind_t *bind;
	size_t index;

	bind = pointer;
	if (bind == NULL)
		return;
	for (index = 0; index < 4; index++) {
		if (bind->arguments[index] != NULL)
			duckdb_free(bind->arguments[index]);
	}
	free(bind);
}

static void
duckvep_model_load_bind(duckdb_bind_info info)
{
	duckvep_load_bind_t *bind;
	duckdb_value complete_value;
	duckdb_logical_type bool_type;
	idx_t parameter_count;
	size_t index;

	parameter_count = duckdb_bind_get_parameter_count(info);
	if (parameter_count != 4) {
		duckdb_bind_set_error(info,
		    "duckvep_model_load: expected four positional arguments");
		return;
	}
	bind = calloc(1, sizeof(*bind));
	if (bind == NULL) {
		duckdb_bind_set_error(info, "duckvep_model_load: out of memory");
		return;
	}
	bind->registry = duckdb_bind_get_extra_info(info);
	for (index = 0; index < 4; index++) {
		bind->arguments[index] = duckvep_bind_string(info, (idx_t)index);
		if (bind->arguments[index] == NULL ||
		    bind->arguments[index][0] == '\0') {
			duckdb_bind_set_error(info,
			    "duckvep_model_load: arguments must be non-empty strings");
			duckvep_model_load_bind_destroy(bind);
			return;
		}
	}
	complete_value = duckdb_bind_get_named_parameter(
	    info, "transcript_coverage_complete");
	if (complete_value != NULL) {
		if (duckdb_is_null_value(complete_value)) {
			duckdb_destroy_value(&complete_value);
			duckdb_bind_set_error(info,
			    "duckvep_model_load: transcript_coverage_complete cannot be NULL");
			duckvep_model_load_bind_destroy(bind);
			return;
		}
		bind->transcript_coverage_complete = duckdb_get_bool(complete_value);
		duckdb_destroy_value(&complete_value);
	}
	bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	duckdb_bind_add_result_column(info, "loaded", bool_type);
	duckdb_destroy_logical_type(&bool_type);
	duckdb_bind_set_bind_data(info, bind,
	    duckvep_model_load_bind_destroy);
}

static void
duckvep_model_load_state_destroy(void *pointer)
{
	free(pointer);
}

static void
duckvep_model_load_init(duckdb_init_info info)
{
	duckvep_load_bind_t *bind;
	duckvep_load_state_t *state;
	duckvep_registry_t *registry;
	duckvep_model_entry_t *entry;
	char error[DUCKVEP_SQL_ERROR_SIZE];
	int loaded;

	duckdb_init_set_max_threads(info, 1);
	bind = duckdb_init_get_bind_data(info);
	if (bind == NULL || bind->registry == NULL) {
		duckdb_init_set_error(info,
		    "duckvep_model_load: missing bind state");
		return;
	}
	state = calloc(1, sizeof(*state));
	if (state == NULL) {
		duckdb_init_set_error(info, "duckvep_model_load: out of memory");
		return;
	}
	registry = bind->registry;
	pthread_mutex_lock(&registry->mutex);
	entry = duckvep_registry_find_locked(registry, bind->arguments[0]);
	pthread_mutex_unlock(&registry->mutex);
	if (entry != NULL) {
		free(state);
		duckdb_init_set_error(info,
		    "duckvep_model_load: model name already exists");
		return;
	}
	entry = calloc(1, sizeof(*entry));
	if (entry == NULL ||
	    (entry->name = duckvep_string_copy(bind->arguments[0])) == NULL) {
		duckvep_model_entry_destroy(entry);
		free(state);
		duckdb_init_set_error(info, "duckvep_model_load: out of memory");
		return;
	}
	memset(error, 0, sizeof(error));
	pthread_mutex_lock(&registry->query_mutex);
	loaded = duckvep_model_load_queries(registry->query_connection,
	    bind->arguments[1], bind->arguments[2], bind->arguments[3],
	    bind->transcript_coverage_complete, &entry->model, error,
	    sizeof(error));
	pthread_mutex_unlock(&registry->query_mutex);
	if (!loaded) {
		duckvep_model_entry_destroy(entry);
		free(state);
		duckdb_init_set_error(info, error[0] != '\0' ? error :
		    "duckvep_model_load: model load failed");
		return;
	}
	pthread_mutex_lock(&registry->mutex);
	if (duckvep_registry_find_locked(registry, entry->name) != NULL) {
		pthread_mutex_unlock(&registry->mutex);
		duckvep_model_entry_destroy(entry);
		free(state);
		duckdb_init_set_error(info,
		    "duckvep_model_load: model name was created concurrently");
		return;
	}
	entry->next = registry->models;
	registry->models = entry;
	pthread_mutex_unlock(&registry->mutex);
	duckdb_init_set_init_data(info, state, duckvep_model_load_state_destroy);
}

static void
duckvep_model_load_scan(duckdb_function_info info, duckdb_data_chunk output)
{
	duckvep_load_state_t *state;
	duckdb_vector vector;
	bool *values;

	state = duckdb_function_get_init_data(info);
	if (state->emitted) {
		duckdb_data_chunk_set_size(output, 0);
		return;
	}
	vector = duckdb_data_chunk_get_vector(output, 0);
	values = duckdb_vector_get_data(vector);
	values[0] = true;
	state->emitted = 1;
	duckdb_data_chunk_set_size(output, 1);
}

static void
duckvep_model_drop_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_registry_t *registry;
	duckdb_vector name_vector;
	bool *values;
	idx_t row, rows;

	registry = duckdb_scalar_function_get_extra_info(info);
	rows = duckdb_data_chunk_get_size(input);
	name_vector = duckdb_data_chunk_get_vector(input, 0);
	values = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		duckvep_model_entry_t *entry, *previous;
		char *name;

		name = duckvep_vector_string(name_vector, row);
		if (name == NULL || *name == '\0') {
			free(name);
			duckdb_scalar_function_set_error(info,
			    "duckvep_model_drop: name must be a non-empty string");
			return;
		}
		entry = NULL;
		previous = NULL;
		pthread_mutex_lock(&registry->mutex);
		for (entry = registry->models; entry != NULL;
		    previous = entry, entry = entry->next) {
			if (strcmp(entry->name, name) == 0)
				break;
		}
		if (entry != NULL && entry->pins == 0) {
			if (previous != NULL)
				previous->next = entry->next;
			else
				registry->models = entry->next;
			entry->next = NULL;
		} else {
			entry = NULL;
		}
		pthread_mutex_unlock(&registry->mutex);
		free(name);
		values[row] = entry != NULL;
		duckvep_model_entry_destroy(entry);
	}
}



duckvep_registry_t *
duckvep_registry_create(duckdb_database database)
{
	duckvep_registry_t *registry;

	registry = calloc(1, sizeof(*registry));
	if (registry == NULL)
		return NULL;
	(void)pthread_mutex_init(&registry->mutex, NULL);
	(void)pthread_mutex_init(&registry->query_mutex, NULL);
	registry->database = database;
	registry->references = 1;
	if (duckdb_connect(database, &registry->query_connection) !=
	    DuckDBSuccess || registry->query_connection == NULL) {
		pthread_mutex_destroy(&registry->query_mutex);
		pthread_mutex_destroy(&registry->mutex);
		free(registry);
		return NULL;
	}
	return registry;
}

void
duckvep_register_model_functions(duckdb_connection connection,
	duckvep_registry_t *registry)
{
	duckdb_table_function table;
	duckdb_scalar_function scalar;
	duckdb_logical_type varchar_type, bool_type;

	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	table = duckdb_create_table_function();
	duckdb_table_function_set_name(table, "duckvep_model_load");
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "transcript_coverage_complete", bool_type);
	duckvep_registry_retain(registry);
	duckdb_table_function_set_extra_info(table, registry,
	    duckvep_registry_release);
	duckdb_table_function_set_bind(table, duckvep_model_load_bind);
	duckdb_table_function_set_init(table, duckvep_model_load_init);
	duckdb_table_function_set_function(table, duckvep_model_load_scan);
	(void)duckdb_register_table_function(connection, table);
	duckdb_destroy_table_function(&table);

	scalar = duckdb_create_scalar_function();
	duckdb_scalar_function_set_name(scalar, "duckvep_model_drop");
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_set_return_type(scalar, bool_type);
	duckdb_scalar_function_set_volatile(scalar);
	duckvep_registry_retain(registry);
	duckdb_scalar_function_set_extra_info(scalar, registry,
	    duckvep_registry_release);
	duckdb_scalar_function_set_function(scalar, duckvep_model_drop_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);

	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&bool_type);
}
