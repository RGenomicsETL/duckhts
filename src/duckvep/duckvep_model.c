#include "duckvep_model.h"
#include "kernel/src/duckvep_model_internal.h"

#include <htslib/faidx.h>

DUCKDB_EXTENSION_EXTERN

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#if defined(_WIN32)
#include <io.h>
#include <share.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

typedef struct duckvep_query_result {
	duckdb_prepared_statement statement;
	duckdb_result result;
	int have_result;
} duckvep_query_result_t;

typedef struct duckvep_load_bind {
	duckvep_registry_t *registry;
	char *arguments[4];
	char *mature_mirna_query;
	char *peptide_edit_query;
	char *interval_feature_query;
	char *reference_fasta;
	int transcript_coverage_complete;
} duckvep_load_bind_t;

typedef struct duckvep_load_state {
	int emitted;
} duckvep_load_state_t;

static char *duckvep_string_copy(const char *);

static void
duckvep_reference_descriptor_close(int descriptor)
{
	if (descriptor < 0)
		return;
#if defined(_WIN32)
	(void)_close(descriptor);
#else
	(void)close(descriptor);
#endif
}

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
	if (!duckvep_sql_resize((void **)&model->sequence_names,
	    sizeof(*model->sequence_names), capacity))
		return 0;
	model->known_seq_region_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_transcripts(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->transcript_capacity)
		return 1;
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
	DUCKVEP_RESIZE_TRANSCRIPT(pre_cds_sequence_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(pre_cds_sequence_lengths);
	DUCKVEP_RESIZE_TRANSCRIPT(post_cds_sequence_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(post_cds_sequence_lengths);
#undef DUCKVEP_RESIZE_TRANSCRIPT
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

static int
duckvep_model_reserve_mature_mirna(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->mature_mirna_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->mature_mirna_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->mature_mirna_starts,
	    sizeof(*model->mature_mirna_starts), capacity) ||
	    !duckvep_sql_resize((void **)&model->mature_mirna_ends,
	    sizeof(*model->mature_mirna_ends), capacity))
		return 0;
	model->mature_mirna_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_peptide_edits(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->peptide_edit_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->peptide_edit_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->peptide_edit_positions,
	    sizeof(*model->peptide_edit_positions), capacity) ||
	    !duckvep_sql_resize((void **)&model->peptide_edit_alts,
	    sizeof(*model->peptide_edit_alts), capacity))
		return 0;
	model->peptide_edit_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_flanks(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->flank_sequence_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->flank_sequence_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->flank_sequence_bytes,
	    sizeof(*model->flank_sequence_bytes), capacity))
		return 0;
	model->flank_sequence_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_interval_features(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->interval_feature_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->interval_feature_capacity,
	    needed);
#define DUCKVEP_RESIZE_INTERVAL_FEATURE(member) \
	if (!duckvep_sql_resize((void **)&model->member, \
	    sizeof(*model->member), capacity)) \
		return 0
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_seq_regions);
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_starts);
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_ends);
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_kinds);
#undef DUCKVEP_RESIZE_INTERVAL_FEATURE
	model->interval_feature_capacity = capacity;
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
	if (model->sequence_names != NULL) {
		size_t region;

		for (region = 0; region < model->known_seq_region_count; region++)
			free(model->sequence_names[region]);
	}
	free(model->sequence_names);
	free(model->reference_fasta_path);
	free(model->reference_fai_path);
	free(model->reference_gzi_path);
	free(model->reference_fasta_open_path);
	free(model->reference_fai_open_path);
	free(model->reference_gzi_open_path);
	if (model->reference_descriptors_open) {
		duckvep_reference_descriptor_close(
		    model->reference_fasta_descriptor);
		duckvep_reference_descriptor_close(model->reference_fai_descriptor);
		duckvep_reference_descriptor_close(model->reference_gzi_descriptor);
	}
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
	free(model->pre_cds_sequence_offsets);
	free(model->pre_cds_sequence_lengths);
	free(model->post_cds_sequence_offsets);
	free(model->post_cds_sequence_lengths);
	free(model->exon_starts);
	free(model->exon_ends);
	free(model->exon_cdna_starts);
	free(model->exon_cdna_ends);
	free(model->exon_phases);
	free(model->exon_end_phases);
	free(model->mature_mirna_offsets);
	free(model->mature_mirna_starts);
	free(model->mature_mirna_ends);
	free(model->peptide_edit_offsets);
	free(model->peptide_edit_positions);
	free(model->peptide_edit_alts);
	free(model->cds_sequence_bytes);
	free(model->flank_sequence_bytes);
	free(model->interval_feature_seq_regions);
	free(model->interval_feature_starts);
	free(model->interval_feature_ends);
	free(model->interval_feature_kinds);
	if (model->interval_index != NULL) {
		/* cgranges 0.1.1 does not release its interval array. */
		free(model->interval_index->r);
		model->interval_index->r = NULL;
		cr_destroy(model->interval_index);
	}
	if (model->interval_feature_index != NULL) {
		/* cgranges 0.1.1 does not release its interval array. */
		free(model->interval_feature_index->r);
		model->interval_feature_index->r = NULL;
		cr_destroy(model->interval_feature_index);
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
	model->transcripts.mature_mirna_offset = model->mature_mirna_offsets;
	model->transcripts.mature_mirna_start1 = model->mature_mirna_starts;
	model->transcripts.mature_mirna_end1 = model->mature_mirna_ends;
	model->transcripts.mature_mirna_count = model->mature_mirna_count;
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
	model->sequences.transcript_count = model->transcripts.transcript_count;
	model->sequences.peptide_edit_offset = model->peptide_edit_offsets;
	model->sequences.peptide_edit_position1 = model->peptide_edit_positions;
	model->sequences.peptide_edit_alt = model->peptide_edit_alts;
	model->sequences.peptide_edit_count = model->peptide_edit_count;
	model->sequences.flank_bytes = model->flank_sequence_bytes;
	model->sequences.flank_bytes_len = model->flank_sequence_length;
	model->sequences.pre_cds_offset = model->pre_cds_sequence_offsets;
	model->sequences.pre_cds_length = model->pre_cds_sequence_lengths;
	model->sequences.post_cds_offset = model->post_cds_sequence_offsets;
	model->sequences.post_cds_length = model->post_cds_sequence_lengths;
	model->sequences.flanks_complete =
	    (uint8_t)(model->transcript_flanks_complete != 0);
	model->interval_features.chrom_id =
	    model->interval_feature_seq_regions;
	model->interval_features.start1 = model->interval_feature_starts;
	model->interval_features.end1 = model->interval_feature_ends;
	model->interval_features.kind = model->interval_feature_kinds;
	model->interval_features.feature_count = model->interval_feature_count;
}

static int
duckvep_owned_interval_index(duckvep_owned_model_t *model,
	const uint16_t *seq_regions, const uint32_t *starts,
	const uint32_t *ends, size_t count, const char *object_name,
	cgranges_t **result, int *complete, char *error, size_t error_size)
{
	size_t index;
	char message[128], region_name[16];

	*complete = 0;
	if (count > (size_t)INT32_MAX)
		return 1;
	for (index = 0; index < count; index++) {
		if (starts[index] > (uint32_t)INT32_MAX ||
		    ends[index] > (uint32_t)INT32_MAX)
			return 1;
	}
	*result = cr_init();
	if (*result == NULL) {
		(void)snprintf(message, sizeof(message),
		    "out of memory building %s interval index", object_name);
		duckvep_sql_set_error(error, error_size, message);
		return 0;
	}
	for (index = 0; index < model->known_seq_region_count; index++) {
		(void)snprintf(region_name, sizeof(region_name), "%u",
		    model->known_seq_regions[index]);
		(void)cr_add_ctg(*result, region_name, 0);
	}
	for (index = 0; index < count; index++) {
		(void)snprintf(region_name, sizeof(region_name), "%u",
		    seq_regions[index]);
		if (cr_add(*result, region_name,
		    (int32_t)(starts[index] - 1), (int32_t)ends[index],
		    (int32_t)index) == NULL) {
			(void)snprintf(message, sizeof(message),
			    "could not build %s interval index", object_name);
			duckvep_sql_set_error(error, error_size, message);
			return 0;
		}
	}
	cr_index(*result);
	*complete = 1;
	return 1;
}

static int
duckvep_owned_model_index(duckvep_owned_model_t *model, char *error,
	size_t error_size)
{
	return duckvep_owned_interval_index(model, model->seq_regions,
	    model->transcript_starts, model->transcript_ends,
	    model->transcripts.transcript_count, "transcript",
	    &model->interval_index, &model->interval_index_complete,
	    error, error_size) &&
	    duckvep_owned_interval_index(model,
	    model->interval_feature_seq_regions, model->interval_feature_starts,
	    model->interval_feature_ends, model->interval_feature_count,
	    "regulation feature", &model->interval_feature_index,
	    &model->interval_feature_index_complete, error, error_size);
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
	int require_sequence_names, duckvep_owned_model_t *model,
	char *error, size_t error_size)
{
	static const char *const one_name[] = {"seq_region"};
	static const duckdb_type one_type[] = {DUCKDB_TYPE_UINTEGER};
	static const char *const two_names[] = {"seq_region", "sequence_length"};
	static const duckdb_type two_types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT
	};
	static const char *const three_names[] = {
		"seq_region", "sequence_length", "seq_region_name"
	};
	static const duckdb_type three_types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_VARCHAR
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
	if (column_count != 1 && column_count != 2 && column_count != 3) {
		duckvep_sql_set_error(error, error_size,
		    "seq_region query must return seq_region, optional sequence_length, and optional seq_region_name");
		goto done;
	}
	if (model->transcript_coverage_complete && column_count < 2) {
		duckvep_sql_set_error(error, error_size,
		    "complete transcript coverage requires sequence_length for every region");
		goto done;
	}
	if (require_sequence_names && column_count != 3) {
		duckvep_sql_set_error(error, error_size,
		    "reference_fasta requires seq_region, sequence_length, and seq_region_name in the region query");
		goto done;
	}
	if (!duckvep_result_schema(&query_result.result,
	    column_count == 1 ? one_name :
	    (column_count == 2 ? two_names : three_names),
	    column_count == 1 ? one_type :
	    (column_count == 2 ? two_types : three_types),
	    (size_t)column_count, SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector region_vector, length_vector, name_vector;
		uint32_t *values;
		uint64_t *lengths;
		idx_t row, rows;

		rows = duckdb_data_chunk_get_size(chunk);
		region_vector = duckdb_data_chunk_get_vector(chunk, 0);
		length_vector = column_count == 2
		    ? duckdb_data_chunk_get_vector(chunk, 1) :
		    (column_count == 3 ? duckdb_data_chunk_get_vector(chunk, 1) : NULL);
		name_vector = column_count == 3
		    ? duckdb_data_chunk_get_vector(chunk, 2) : NULL;
		values = (uint32_t *)duckdb_vector_get_data(region_vector);
		lengths = length_vector != NULL
		    ? (uint64_t *)duckdb_vector_get_data(length_vector) : NULL;
		for (row = 0; row < rows; row++) {
			size_t index;

			if (duckvep_row_is_null(region_vector, row) ||
			    (length_vector != NULL &&
			    duckvep_row_is_null(length_vector, row)) ||
			    (name_vector != NULL &&
			    duckvep_row_is_null(name_vector, row))) {
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
			model->sequence_names[index] = NULL;
			if (name_vector != NULL) {
				model->sequence_names[index] =
				    duckvep_vector_string(name_vector, row);
				if (model->sequence_names[index] == NULL ||
				    model->sequence_names[index][0] == '\0') {
					free(model->sequence_names[index]);
					model->sequence_names[index] = NULL;
					duckvep_sql_set_error(error, error_size,
					    "seq_region_name must be a non-empty string without embedded NUL bytes");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
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
duckvep_sequence_name_compare(const void *left, const void *right)
{
	const char *const *a, *const *b;

	a = left;
	b = right;
	return strcmp(*a, *b);
}

static char *
duckvep_reference_path_with_suffix(const char *path, const char *suffix)
{
	size_t path_length, suffix_length;
	char *result;

	if (path == NULL || suffix == NULL)
		return NULL;
	path_length = strlen(path);
	suffix_length = strlen(suffix);
	if (path_length > SIZE_MAX - suffix_length - 1u)
		return NULL;
	result = malloc(path_length + suffix_length + 1u);
	if (result == NULL)
		return NULL;
	memcpy(result, path, path_length);
	memcpy(result + path_length, suffix, suffix_length + 1u);
	return result;
}

static int
duckvep_reference_file_identity_read_descriptor(int descriptor,
	duckvep_reference_file_identity_t *identity)
{
#if defined(_WIN32)
	struct _stat64 status;
#else
	struct stat status;
#endif

	if (identity == NULL)
		return 0;
	memset(identity, 0, sizeof(*identity));
	if (descriptor < 0)
		return 1;
	errno = 0;
#if defined(_WIN32)
	if (_fstat64(descriptor, &status) != 0)
#else
	if (fstat(descriptor, &status) != 0)
#endif
		return 0;
	if (status.st_size < 0)
		return 0;
	identity->device = (uint64_t)status.st_dev;
	identity->inode = (uint64_t)status.st_ino;
	identity->size = (uint64_t)status.st_size;
	identity->mtime_seconds = (int64_t)status.st_mtime;
#if defined(__APPLE__)
	identity->mtime_nanoseconds =
	    (uint32_t)status.st_mtimespec.tv_nsec;
#elif !defined(_WIN32)
	identity->mtime_nanoseconds = (uint32_t)status.st_mtim.tv_nsec;
#endif
	identity->present = 1;
	return 1;
}

static int
duckvep_reference_descriptor_open(const char *path, int required,
	int *descriptor_out, duckvep_reference_file_identity_t *identity)
{
	int descriptor;

	if (path == NULL || descriptor_out == NULL || identity == NULL)
		return 0;
	*descriptor_out = -1;
	memset(identity, 0, sizeof(*identity));
	errno = 0;
#if defined(_WIN32)
	errno_t open_error;

	descriptor = -1;
	open_error = _sopen_s(&descriptor, path, _O_RDONLY | _O_BINARY,
	    _SH_DENYWR, 0);
	if (open_error != 0) {
		if (!required && open_error == ENOENT)
			return 1;
		return 0;
	}
#else
	{
		int flags;

		flags = O_RDONLY;
#if defined(O_CLOEXEC)
		flags |= O_CLOEXEC;
#endif
		descriptor = open(path, flags);
	}
	if (descriptor < 0) {
		if (!required && errno == ENOENT)
			return 1;
		return 0;
	}
#endif
	if (!duckvep_reference_file_identity_read_descriptor(descriptor,
	    identity)) {
		duckvep_reference_descriptor_close(descriptor);
		return 0;
	}
	*descriptor_out = descriptor;
	return 1;
}

static char *
duckvep_reference_descriptor_path(int descriptor, const char *source_path)
{
#if defined(_WIN32)
	HANDLE handle;
	DWORD needed, written;
	char *path;

	(void)source_path;
	if (descriptor < 0)
		return NULL;
	handle = (HANDLE)_get_osfhandle(descriptor);
	if (handle == INVALID_HANDLE_VALUE)
		return NULL;
	needed = GetFinalPathNameByHandleA(handle, NULL, 0,
	    FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);
	if (needed == 0 || needed > SIZE_MAX - 1u)
		return NULL;
	path = malloc((size_t)needed + 1u);
	if (path == NULL)
		return NULL;
	written = GetFinalPathNameByHandleA(handle, path, needed + 1u,
	    FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);
	if (written == 0 || written > needed) {
		free(path);
		return NULL;
	}
	return path;
#else
	if (descriptor < 0)
		return NULL;
#if defined(__linux__)
	char buffer[64];
	int written;

	(void)source_path;
	written = snprintf(buffer, sizeof(buffer), "/proc/self/fd/%d",
	    descriptor);
	if (written < 0 || (size_t)written >= sizeof(buffer))
		return NULL;
	return duckvep_string_copy(buffer);
#elif defined(__APPLE__)
	{
		char path[PATH_MAX];

		(void)source_path;
		if (fcntl(descriptor, F_GETPATH, path) != 0)
			return NULL;
		return duckvep_string_copy(path);
	}
#else
	return duckvep_string_copy(source_path);
#endif
#endif
}

static int
duckvep_reference_file_identity_equal(
	const duckvep_reference_file_identity_t *left,
	const duckvep_reference_file_identity_t *right)
{
	if (left == NULL || right == NULL || left->present != right->present)
		return 0;
	if (!left->present)
		return 1;
	/*
	 * POSIX rename can update ctime without changing the open file or its
	 * contents.  Device/inode detect a different object; size/mtime detect
	 * ordinary in-place mutation without rejecting a safely pinned rename.
	 */
	return left->device == right->device &&
	    left->inode == right->inode && left->size == right->size &&
	    left->mtime_seconds == right->mtime_seconds &&
	    left->mtime_nanoseconds == right->mtime_nanoseconds;
}

#if !defined(__linux__) && !defined(_WIN32)
static int
duckvep_reference_open_path_identity_matches(const char *path,
	const duckvep_reference_file_identity_t *expected)
{
	duckvep_reference_file_identity_t observed;
	int descriptor;
	int matches;

	if (expected == NULL)
		return 0;
	if (!expected->present)
		return path == NULL;
	if (path == NULL || !duckvep_reference_descriptor_open(path, 1,
	    &descriptor, &observed))
		return 0;
	matches = duckvep_reference_file_identity_equal(expected, &observed);
	duckvep_reference_descriptor_close(descriptor);
	return matches;
}
#endif

int
duckvep_model_reference_identity_matches(const duckvep_owned_model_t *model)
{
	duckvep_reference_file_identity_t fasta, fai, gzi;

	if (model == NULL)
		return 0;
	if (model->reference_fasta_path == NULL)
		return 1;
	if (!model->reference_descriptors_open ||
	    !duckvep_reference_file_identity_read_descriptor(
	    model->reference_fasta_descriptor, &fasta) ||
	    !duckvep_reference_file_identity_read_descriptor(
	    model->reference_fai_descriptor, &fai) ||
	    !duckvep_reference_file_identity_read_descriptor(
	    model->reference_gzi_descriptor, &gzi))
		return 0;
#if !defined(__linux__) && !defined(_WIN32)
	if (!duckvep_reference_open_path_identity_matches(
	    model->reference_fasta_open_path,
	    &model->reference_fasta_identity) ||
	    !duckvep_reference_open_path_identity_matches(
	    model->reference_fai_open_path,
	    &model->reference_fai_identity) ||
	    !duckvep_reference_open_path_identity_matches(
	    model->reference_gzi_open_path,
	    &model->reference_gzi_identity))
		return 0;
#endif
	return duckvep_reference_file_identity_equal(
	    &model->reference_fasta_identity, &fasta) &&
	    duckvep_reference_file_identity_equal(
	    &model->reference_fai_identity, &fai) &&
	    duckvep_reference_file_identity_equal(
	    &model->reference_gzi_identity, &gzi);
}

static int
duckvep_validate_reference_fasta(const char *reference_fasta,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	faidx_t *fai;
	char **sorted_names;
	size_t region;
	int ok;

	if (reference_fasta == NULL)
		return 1;
	model->reference_fasta_path = duckvep_string_copy(reference_fasta);
	model->reference_fai_path = duckvep_reference_path_with_suffix(
	    reference_fasta, ".fai");
	model->reference_gzi_path = duckvep_reference_path_with_suffix(
	    reference_fasta, ".gzi");
	if (model->reference_fasta_path == NULL ||
	    model->reference_fai_path == NULL ||
	    model->reference_gzi_path == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory retaining reference FASTA identity paths");
		return 0;
	}
	model->reference_fasta_descriptor = -1;
	model->reference_fai_descriptor = -1;
	model->reference_gzi_descriptor = -1;
	model->reference_descriptors_open = 1;
	if (!duckvep_reference_descriptor_open(model->reference_fasta_path, 1,
	    &model->reference_fasta_descriptor,
	    &model->reference_fasta_identity)) {
		duckvep_sql_set_error(error, error_size,
		    "could not pin the reference FASTA");
		return 0;
	}
	if (!duckvep_reference_descriptor_open(model->reference_fai_path, 1,
	    &model->reference_fai_descriptor,
	    &model->reference_fai_identity)) {
		duckvep_sql_set_error(error, error_size,
		    "could not pin the reference FASTA .fai index");
		return 0;
	}
	if (!duckvep_reference_descriptor_open(model->reference_gzi_path, 0,
	    &model->reference_gzi_descriptor,
	    &model->reference_gzi_identity)) {
		duckvep_sql_set_error(error, error_size,
		    "could not pin the reference FASTA .gzi index");
		return 0;
	}
	model->reference_fasta_open_path = duckvep_reference_descriptor_path(
	    model->reference_fasta_descriptor, model->reference_fasta_path);
	model->reference_fai_open_path = duckvep_reference_descriptor_path(
	    model->reference_fai_descriptor, model->reference_fai_path);
	if (model->reference_gzi_identity.present)
		model->reference_gzi_open_path =
		    duckvep_reference_descriptor_path(
		    model->reference_gzi_descriptor, model->reference_gzi_path);
	if (model->reference_fasta_open_path == NULL ||
	    model->reference_fai_open_path == NULL ||
	    (model->reference_gzi_identity.present &&
	    model->reference_gzi_open_path == NULL)) {
		duckvep_sql_set_error(error, error_size,
		    "could not retain stable reference descriptor paths");
		return 0;
	}
	if (model->known_seq_region_count == 0 ||
	    model->known_seq_region_count > SIZE_MAX / sizeof(*sorted_names)) {
		duckvep_sql_set_error(error, error_size,
		    "reference FASTA region map exceeds addressable memory");
		return 0;
	}
	sorted_names = malloc(model->known_seq_region_count *
	    sizeof(*sorted_names));
	if (sorted_names == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory validating reference FASTA region names");
		return 0;
	}
	for (region = 0; region < model->known_seq_region_count; region++) {
		if (model->sequence_names[region] == NULL ||
		    model->sequence_lengths[region] == 0) {
			free(sorted_names);
			duckvep_sql_set_error(error, error_size,
			    "reference FASTA requires a name and length for every model region");
			return 0;
		}
		sorted_names[region] = model->sequence_names[region];
	}
	qsort(sorted_names, model->known_seq_region_count,
	    sizeof(*sorted_names), duckvep_sequence_name_compare);
	for (region = 1; region < model->known_seq_region_count; region++) {
		if (strcmp(sorted_names[region - 1], sorted_names[region]) == 0) {
			free(sorted_names);
			duckvep_sql_set_error(error, error_size,
			    "seq_region_name values must be unique");
			return 0;
		}
	}
	free(sorted_names);
	if (!duckvep_model_reference_identity_matches(model)) {
		duckvep_sql_set_error(error, error_size,
		    "reference FASTA or index changed while the model was loading");
		return 0;
	}

	fai = fai_load3_format(model->reference_fasta_open_path,
	    model->reference_fai_open_path, model->reference_gzi_open_path,
	    0, FAI_FASTA);
	if (fai == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "could not open the indexed reference FASTA without creating an index");
		return 0;
	}
	ok = 1;
	for (region = 0; region < model->known_seq_region_count; region++) {
		hts_pos_t length;

		length = faidx_seq_len64(fai, model->sequence_names[region]);
		if (length < 0 || (uint64_t)length !=
		    (uint64_t)model->sequence_lengths[region]) {
			(void)snprintf(error, error_size,
			    "reference FASTA contig %s is absent or has the wrong length",
			    model->sequence_names[region]);
			ok = 0;
			break;
		}
	}
	fai_destroy(fai);
	if (!ok)
		return 0;
	if (!duckvep_model_reference_identity_matches(model)) {
		duckvep_sql_set_error(error, error_size,
		    "reference FASTA or index changed while the model was loading");
		return 0;
	}
	return 1;
}

static int
duckvep_load_transcripts(duckdb_connection connection, const char *query,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const base_names[] = {
		"transcript_index", "seq_region", "transcript_start",
		"transcript_end", "strand", "gene_index", "transcript_flags",
		"cds_start", "cds_end", "cds_sequence", "codon_table"
	};
	static const char *const legacy_names[] = {
		"transcript_index", "seq_region", "transcript_start",
		"transcript_end", "strand", "gene_index", "transcript_flags",
		"cds_start", "cds_end", "cds_sequence", "codon_table",
		"post_cds_bases"
	};
	static const char *const complete_names[] = {
		"transcript_index", "seq_region", "transcript_start",
		"transcript_end", "strand", "gene_index", "transcript_flags",
		"cds_start", "cds_end", "cds_sequence", "codon_table",
		"pre_cds_sequence", "post_cds_sequence"
	};
	static const duckdb_type types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
		DUCKDB_TYPE_TINYINT, DUCKDB_TYPE_UINTEGER,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
		DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_BLOB,
		DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_BLOB,
		DUCKDB_TYPE_BLOB
	};
	const char *const *names;
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	idx_t column_count;
	int ok;

	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	ok = 0;
	column_count = duckdb_column_count(&query_result.result);
	if (column_count != 11 && column_count != 12 && column_count != 13) {
		duckvep_sql_set_error(error, error_size,
		    "transcript query must return 11 columns, 12 with legacy post_cds_bases, or 13 with complete pre_cds_sequence and post_cds_sequence");
		goto done;
	}
	names = column_count == 11 ? base_names :
	    (column_count == 12 ? legacy_names : complete_names);
	if (!duckvep_result_schema(&query_result.result, names, types,
	    (size_t)column_count, 9,
	    error, error_size))
		goto done;
	model->transcript_flanks_complete = column_count == 13;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector vectors[13];
		uint32_t *transcript_indices, *seq_regions, *gene_indices;
		uint64_t *starts, *ends, *flags, *cds_starts, *cds_ends;
		int8_t *strands;
		duckdb_string_t *sequences, *legacy_tails, *pre_flanks;
		duckdb_string_t *post_flanks;
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
		legacy_tails = column_count == 12
		    ? duckdb_vector_get_data(vectors[11]) : NULL;
		pre_flanks = column_count == 13
		    ? duckdb_vector_get_data(vectors[11]) : NULL;
		post_flanks = column_count == 13
		    ? duckdb_vector_get_data(vectors[12]) : NULL;
		for (row = 0; row < rows; row++) {
			size_t flank_offset, index, post_length, pre_length;
			size_t sequence_length, sequence_offset;
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
			if (column_count == 13) {
				int pre_null = duckvep_row_is_null(vectors[11], row);
				int post_null = duckvep_row_is_null(vectors[12], row);

				if (pre_null != post_null ||
				    ((sequence_nulls == 0) == (pre_null != 0))) {
					duckvep_sql_set_error(error, error_size,
					    "complete transcript flanks must both be present exactly when CDS sequence is present");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
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
			if (sequence_nulls == 0) {
				sequence_length = (size_t)duckdb_string_t_length(
				    sequences[row]);
				if (sequence_length == 0 || sequence_length > UINT32_MAX ||
				    !duckvep_codon_table_supported(
				    (duckvep_codon_table_t)tables[row])) {
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
			flank_offset = model->flank_sequence_length;
			pre_length = 0;
			post_length = 0;
			if (legacy_tails != NULL &&
			    !duckvep_row_is_null(vectors[11], row)) {
				post_length = (size_t)duckdb_string_t_length(
				    legacy_tails[row]);
				if (post_length > 3u ||
				    (post_length != 0 && sequence_nulls != 0)) {
					duckvep_sql_set_error(error, error_size,
					    "post_cds_bases must contain at most three bases for a sequence-backed coding transcript");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			} else if (pre_flanks != NULL) {
				if (sequence_nulls == 0) {
					pre_length = (size_t)duckdb_string_t_length(
					    pre_flanks[row]);
					post_length = (size_t)duckdb_string_t_length(
					    post_flanks[row]);
				}
				if (pre_length > UINT32_MAX || post_length > UINT32_MAX) {
					duckvep_sql_set_error(error, error_size,
					    "transcript flank exceeds the uint32 model limit");
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			if (pre_length > SIZE_MAX - flank_offset ||
			    post_length > SIZE_MAX - flank_offset - pre_length ||
			    !duckvep_model_reserve_flanks(model,
			    flank_offset + pre_length + post_length)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading transcript flanks");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (pre_length != 0)
				memcpy(model->flank_sequence_bytes + flank_offset,
				    duckdb_string_t_data(&pre_flanks[row]), pre_length);
			if (post_length != 0) {
				duckdb_string_t *source = legacy_tails != NULL
				    ? legacy_tails : post_flanks;

				memcpy(model->flank_sequence_bytes + flank_offset +
				    pre_length, duckdb_string_t_data(&source[row]),
				    post_length);
			}
			model->cds_sequence_offsets[index] =
			    (uint64_t)sequence_offset;
			model->cds_sequence_lengths[index] =
			    (uint32_t)sequence_length;
			model->pre_cds_sequence_offsets[index] =
			    (uint64_t)flank_offset;
			model->pre_cds_sequence_lengths[index] =
			    (uint32_t)pre_length;
			model->post_cds_sequence_offsets[index] =
			    (uint64_t)(flank_offset + pre_length);
			model->post_cds_sequence_lengths[index] =
			    (uint32_t)post_length;
			model->cds_sequence_length += sequence_length;
			model->flank_sequence_length += pre_length + post_length;
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
duckvep_mature_mirna_segment_is_exonic(
	const duckvep_owned_model_t *model, uint32_t transcript_index,
	uint32_t start1, uint32_t end1)
{
	size_t offset, count, exon;

	offset = model->exon_offsets[transcript_index];
	count = model->exon_counts[transcript_index];
	for (exon = offset; exon < offset + count; exon++) {
		if (start1 >= model->exon_starts[exon] &&
		    end1 <= model->exon_ends[exon])
			return 1;
	}
	return 0;
}

static int
duckvep_load_mature_mirna(duckdb_connection connection, const char *query,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "mature_mirna_start", "mature_mirna_end"
	};
	static const duckdb_type types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT
	};
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	uint32_t previous_transcript, previous_start;
	int have_previous, ok;
	size_t transcript;

	if (model->transcripts.transcript_count >
	    SIZE_MAX / sizeof(*model->mature_mirna_offsets) - 1u) {
		duckvep_sql_set_error(error, error_size,
		    "mature-miRNA row-offset array exceeds addressable memory");
		return 0;
	}
	model->mature_mirna_offsets = calloc(
	    model->transcripts.transcript_count + 1u,
	    sizeof(*model->mature_mirna_offsets));
	if (model->mature_mirna_offsets == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading mature-miRNA row offsets");
		return 0;
	}
	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	have_previous = 0;
	previous_transcript = 0;
	previous_start = 0;
	ok = 0;
	if (!duckvep_result_schema(&query_result.result, names, types, 3,
	    SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector vectors[3];
		uint32_t *transcript_indices;
		uint64_t *starts, *ends;
		idx_t row, rows;
		size_t column;

		rows = duckdb_data_chunk_get_size(chunk);
		for (column = 0; column < 3; column++)
			vectors[column] = duckdb_data_chunk_get_vector(chunk,
			    (idx_t)column);
		transcript_indices = duckdb_vector_get_data(vectors[0]);
		starts = duckdb_vector_get_data(vectors[1]);
		ends = duckdb_vector_get_data(vectors[2]);
		for (row = 0; row < rows; row++) {
			uint32_t transcript_index;
			size_t feature_index;

			for (column = 0; column < 3; column++) {
				if (duckvep_row_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "mature-miRNA query contains NULL in %s",
					    names[column]);
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			transcript_index = transcript_indices[row];
			if (transcript_index >= model->transcripts.transcript_count ||
			    (have_previous &&
			    (transcript_index < previous_transcript ||
			    (transcript_index == previous_transcript &&
			    starts[row] < previous_start)))) {
				duckvep_sql_set_error(error, error_size,
				    "mature-miRNA query must be ordered by transcript_index and mature_mirna_start");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (starts[row] == 0 || starts[row] > UINT32_MAX ||
			    ends[row] > UINT32_MAX || ends[row] < starts[row] ||
			    starts[row] < model->transcript_starts[transcript_index] ||
			    ends[row] > model->transcript_ends[transcript_index] ||
			    (model->transcript_flags[transcript_index] &
			    (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA) == 0u ||
			    !duckvep_mature_mirna_segment_is_exonic(model,
			    transcript_index, (uint32_t)starts[row],
			    (uint32_t)ends[row])) {
				duckvep_sql_set_error(error, error_size,
				    "mature-miRNA row is not an exonic interval of a miRNA transcript");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			feature_index = model->mature_mirna_count;
			if (feature_index == (size_t)UINT32_MAX ||
			    !duckvep_model_reserve_mature_mirna(model,
			    feature_index + 1u)) {
				duckvep_sql_set_error(error, error_size,
				    "mature-miRNA side relation exceeds the uint32 model limit");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			model->mature_mirna_starts[feature_index] =
			    (uint32_t)starts[row];
			model->mature_mirna_ends[feature_index] =
			    (uint32_t)ends[row];
			model->mature_mirna_offsets[transcript_index + 1u]++;
			model->mature_mirna_count++;
			previous_transcript = transcript_index;
			previous_start = (uint32_t)starts[row];
			have_previous = 1;
		}
		duckdb_destroy_data_chunk(&chunk);
	}
	for (transcript = 1u;
	    transcript <= model->transcripts.transcript_count; transcript++) {
		model->mature_mirna_offsets[transcript] +=
		    model->mature_mirna_offsets[transcript - 1u];
	}
	ok = 1;
done:
	duckvep_query_result_close(&query_result);
	return ok;
}

static int
duckvep_peptide_edit_alt_valid(uint8_t amino_acid)
{
	return amino_acid == (uint8_t)'*' ||
	    (amino_acid >= (uint8_t)'A' && amino_acid <= (uint8_t)'Z');
}

static int
duckvep_load_peptide_edits(duckdb_connection connection, const char *query,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "protein_position", "alternate_amino_acid"
	};
	static const duckdb_type types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR
	};
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	uint32_t previous_transcript, previous_position;
	int have_previous, ok;
	size_t transcript;

	if (model->transcripts.transcript_count >
	    SIZE_MAX / sizeof(*model->peptide_edit_offsets) - 1u) {
		duckvep_sql_set_error(error, error_size,
		    "peptide-edit row-offset array exceeds addressable memory");
		return 0;
	}
	model->peptide_edit_offsets = calloc(
	    model->transcripts.transcript_count + 1u,
	    sizeof(*model->peptide_edit_offsets));
	if (model->peptide_edit_offsets == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading peptide-edit row offsets");
		return 0;
	}
	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	have_previous = 0;
	previous_transcript = 0;
	previous_position = 0;
	ok = 0;
	if (!duckvep_result_schema(&query_result.result, names, types, 3,
	    SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector vectors[3];
		uint32_t *transcript_indices, *positions;
		duckdb_string_t *alternates;
		idx_t row, rows;
		size_t column;

		rows = duckdb_data_chunk_get_size(chunk);
		for (column = 0; column < 3; column++)
			vectors[column] = duckdb_data_chunk_get_vector(chunk,
			    (idx_t)column);
		transcript_indices = duckdb_vector_get_data(vectors[0]);
		positions = duckdb_vector_get_data(vectors[1]);
		alternates = duckdb_vector_get_data(vectors[2]);
		for (row = 0; row < rows; row++) {
			const char *alternate;
			uint32_t transcript_index, position;
			uint32_t alternate_length;
			size_t edit_index;
			uint8_t amino_acid;

			for (column = 0; column < 3; column++) {
				if (duckvep_row_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "peptide-edit query contains NULL in %s",
					    names[column]);
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			transcript_index = transcript_indices[row];
			position = positions[row];
			if (transcript_index >= model->transcripts.transcript_count ||
			    (have_previous &&
			    (transcript_index < previous_transcript ||
			    (transcript_index == previous_transcript &&
			    position <= previous_position)))) {
				duckvep_sql_set_error(error, error_size,
				    "peptide-edit query must be ordered and unique by transcript_index and protein_position");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			alternate_length = duckdb_string_t_length(alternates[row]);
			alternate = duckdb_string_t_data(&alternates[row]);
			amino_acid = alternate_length == 1u
			    ? (uint8_t)alternate[0] : 0u;
			if (position == 0 ||
			    position > model->cds_sequence_lengths[transcript_index] / 3u ||
			    alternate_length != 1u ||
			    !duckvep_peptide_edit_alt_valid(amino_acid)) {
				duckvep_sql_set_error(error, error_size,
				    "peptide edit must replace one in-range protein position with one uppercase amino-acid code");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			edit_index = model->peptide_edit_count;
			if (edit_index == (size_t)UINT32_MAX ||
			    !duckvep_model_reserve_peptide_edits(model,
			    edit_index + 1u)) {
				duckvep_sql_set_error(error, error_size,
				    "peptide-edit side relation exceeds the uint32 model limit");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			model->peptide_edit_positions[edit_index] = position;
			model->peptide_edit_alts[edit_index] = amino_acid;
			model->peptide_edit_offsets[transcript_index + 1u]++;
			model->peptide_edit_count++;
			previous_transcript = transcript_index;
			previous_position = position;
			have_previous = 1;
		}
		duckdb_destroy_data_chunk(&chunk);
	}
	for (transcript = 1u;
	    transcript <= model->transcripts.transcript_count; transcript++) {
		model->peptide_edit_offsets[transcript] +=
		    model->peptide_edit_offsets[transcript - 1u];
	}
	ok = 1;
done:
	duckvep_query_result_close(&query_result);
	return ok;
}

static int
duckvep_model_region_length(const duckvep_owned_model_t *model,
	uint32_t seq_region, uint32_t *sequence_length)
{
	size_t begin, end;

	begin = 0u;
	end = model->known_seq_region_count;
	while (begin < end) {
		size_t middle;

		middle = begin + (end - begin) / 2u;
		if ((uint32_t)model->known_seq_regions[middle] < seq_region)
			begin = middle + 1u;
		else
			end = middle;
	}
	if (begin >= model->known_seq_region_count ||
	    (uint32_t)model->known_seq_regions[begin] != seq_region)
		return 0;
	if (sequence_length != NULL)
		*sequence_length = model->sequence_lengths[begin];
	return 1;
}

static int
duckvep_load_interval_features(duckdb_connection connection,
	const char *query, duckvep_owned_model_t *model, char *error,
	size_t error_size)
{
	static const char *const names[] = {
		"regulation_feature_index", "seq_region", "feature_start",
		"feature_end", "feature_kind"
	};
	static const duckdb_type types[] = {
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER,
		DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER,
		DUCKDB_TYPE_UTINYINT
	};
	duckvep_query_result_t query_result;
	duckdb_data_chunk chunk;
	uint32_t previous_seq_region, previous_start;
	int have_previous, ok;

	if (!duckvep_query_result_open(connection, query, &query_result, error,
	    error_size))
		return 0;
	ok = 0;
	have_previous = 0;
	previous_seq_region = 0u;
	previous_start = 0u;
	if (!duckvep_result_schema(&query_result.result, names, types, 5,
	    SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckdb_fetch_chunk(query_result.result)) != NULL) {
		duckdb_vector vectors[5];
		uint32_t *indices, *seq_regions, *starts, *ends;
		uint8_t *kinds;
		idx_t row, rows;
		size_t column;

		rows = duckdb_data_chunk_get_size(chunk);
		for (column = 0u; column < 5u; column++)
			vectors[column] = duckdb_data_chunk_get_vector(chunk,
			    (idx_t)column);
		indices = duckdb_vector_get_data(vectors[0]);
		seq_regions = duckdb_vector_get_data(vectors[1]);
		starts = duckdb_vector_get_data(vectors[2]);
		ends = duckdb_vector_get_data(vectors[3]);
		kinds = duckdb_vector_get_data(vectors[4]);
		for (row = 0; row < rows; row++) {
			size_t index;
			uint32_t sequence_length;

			for (column = 0u; column < 5u; column++) {
				if (duckvep_row_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "interval-feature query contains NULL in %s",
					    names[column]);
					duckdb_destroy_data_chunk(&chunk);
					goto done;
				}
			}
			index = model->interval_feature_count;
			sequence_length = 0u;
			if ((uint64_t)index > UINT32_MAX ||
			    indices[row] != (uint32_t)index ||
			    seq_regions[row] > UINT16_MAX || starts[row] == 0u ||
			    ends[row] < starts[row] ||
			    (kinds[row] !=
			    (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION &&
			    kinds[row] !=
			    (uint8_t)DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE) ||
			    !duckvep_model_region_length(model, seq_regions[row],
			    &sequence_length) ||
			    (sequence_length != 0u && ends[row] > sequence_length)) {
				duckvep_sql_set_error(error, error_size,
				    "interval-feature row has an invalid dense index, region, coordinate, or kind");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (have_previous &&
			    (seq_regions[row] < previous_seq_region ||
			    (seq_regions[row] == previous_seq_region &&
			    starts[row] < previous_start))) {
				duckvep_sql_set_error(error, error_size,
				    "interval-feature query must be ordered by seq_region, feature_start, and regulation_feature_index");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			if (!duckvep_model_reserve_interval_features(model,
			    index + 1u)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading interval features");
				duckdb_destroy_data_chunk(&chunk);
				goto done;
			}
			model->interval_feature_seq_regions[index] =
			    (uint16_t)seq_regions[row];
			model->interval_feature_starts[index] = starts[row];
			model->interval_feature_ends[index] = ends[row];
			model->interval_feature_kinds[index] = kinds[row];
			model->interval_feature_count++;
			previous_seq_region = seq_regions[row];
			previous_start = starts[row];
			have_previous = 1;
		}
		duckdb_destroy_data_chunk(&chunk);
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
	const char *exon_query, const char *mature_mirna_query,
	const char *peptide_edit_query, const char *interval_feature_query,
	const char *reference_fasta, int transcript_coverage_complete,
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
	ok = duckvep_load_regions(connection, region_query,
	    reference_fasta != NULL, model, error, error_size) &&
	    duckvep_validate_reference_fasta(reference_fasta, model, error,
	    error_size) &&
	    duckvep_load_transcripts(connection, transcript_query, model, error,
	    error_size) &&
	    duckvep_load_exons(connection, exon_query, model, error, error_size) &&
	    (mature_mirna_query == NULL ||
	    duckvep_load_mature_mirna(connection, mature_mirna_query, model,
	    error, error_size)) &&
	    (peptide_edit_query == NULL ||
	    duckvep_load_peptide_edits(connection, peptide_edit_query, model,
	    error, error_size)) &&
	    (interval_feature_query == NULL ||
	    duckvep_load_interval_features(connection, interval_feature_query,
	    model, error, error_size));
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
	if (duckvep_model_open_with_interval_features(&model->transcripts,
	    &model->exons, &model->sequences, &model->interval_features,
	    &model->kernel, &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "invalid transcript model: %s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "kernel validation failed");
		duckvep_owned_model_destroy(model);
		return 0;
	}
	/* The kernel validates the borrowed model once and owns the canonical
	 * immutable projection caches. Publish those prepared views to the adapter
	 * instead of retaining a second transcript authority that lacks the caches. */
	model->transcripts = *duckvep_model_prepared_transcripts(model->kernel);
	model->exons = *duckvep_model_prepared_exons(model->kernel);
	model->sequences = *duckvep_model_prepared_sequences(model->kernel);
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
		duckvep_workspace_cache_destroy(workspace);
	}
	free(entry->name);
	duckvep_owned_model_destroy(&entry->model);
	free(entry);
}

void
duckvep_workspace_cache_destroy(duckvep_workspace_cache_t *cache)
{
	if (cache == NULL)
		return;
	if (cache->reference_fai != NULL)
		fai_destroy(cache->reference_fai);
	free(cache->reference_bases);
	duckvep_workspace_close(cache->workspace);
	free(cache);
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
	if (bind->mature_mirna_query != NULL)
		duckdb_free(bind->mature_mirna_query);
	if (bind->peptide_edit_query != NULL)
		duckdb_free(bind->peptide_edit_query);
	if (bind->interval_feature_query != NULL)
		duckdb_free(bind->interval_feature_query);
	if (bind->reference_fasta != NULL)
		duckdb_free(bind->reference_fasta);
	free(bind);
}

static void
duckvep_model_load_bind(duckdb_bind_info info)
{
	duckvep_load_bind_t *bind;
	duckdb_value complete_value;
	duckdb_value mature_mirna_value;
	duckdb_value peptide_edit_value;
	duckdb_value interval_feature_value;
	duckdb_value reference_fasta_value;
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
	mature_mirna_value = duckdb_bind_get_named_parameter(
	    info, "mature_mirna_query");
	if (mature_mirna_value != NULL) {
		if (duckdb_is_null_value(mature_mirna_value)) {
			duckdb_destroy_value(&mature_mirna_value);
		} else {
			bind->mature_mirna_query = duckdb_get_varchar(
			    mature_mirna_value);
			duckdb_destroy_value(&mature_mirna_value);
			if (bind->mature_mirna_query == NULL ||
			    bind->mature_mirna_query[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: mature_mirna_query must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	peptide_edit_value = duckdb_bind_get_named_parameter(
	    info, "peptide_edit_query");
	if (peptide_edit_value != NULL) {
		if (duckdb_is_null_value(peptide_edit_value)) {
			duckdb_destroy_value(&peptide_edit_value);
		} else {
			bind->peptide_edit_query = duckdb_get_varchar(
			    peptide_edit_value);
			duckdb_destroy_value(&peptide_edit_value);
			if (bind->peptide_edit_query == NULL ||
			    bind->peptide_edit_query[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: peptide_edit_query must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	interval_feature_value = duckdb_bind_get_named_parameter(
	    info, "interval_feature_query");
	if (interval_feature_value != NULL) {
		if (duckdb_is_null_value(interval_feature_value)) {
			duckdb_destroy_value(&interval_feature_value);
		} else {
			bind->interval_feature_query = duckdb_get_varchar(
			    interval_feature_value);
			duckdb_destroy_value(&interval_feature_value);
			if (bind->interval_feature_query == NULL ||
			    bind->interval_feature_query[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: interval_feature_query must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	reference_fasta_value = duckdb_bind_get_named_parameter(
	    info, "reference_fasta");
	if (reference_fasta_value != NULL) {
		if (duckdb_is_null_value(reference_fasta_value)) {
			duckdb_destroy_value(&reference_fasta_value);
		} else {
			bind->reference_fasta = duckdb_get_varchar(
			    reference_fasta_value);
			duckdb_destroy_value(&reference_fasta_value);
			if (bind->reference_fasta == NULL ||
			    bind->reference_fasta[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: reference_fasta must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
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
	    bind->mature_mirna_query,
	    bind->peptide_edit_query,
	    bind->interval_feature_query,
	    bind->reference_fasta, bind->transcript_coverage_complete,
	    &entry->model, error,
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
	    "mature_mirna_query", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "peptide_edit_query", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "interval_feature_query", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "reference_fasta", varchar_type);
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
