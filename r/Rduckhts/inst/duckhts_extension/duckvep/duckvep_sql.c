/* Relation-oriented DuckVEP SQL surface and static annotation metadata. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckvep_so.h"
#include "duckvep_codon.h"
#include "duckvep_sql.h"

typedef struct {
	idx_t offset;
} duckvep_so_scan_t;

bool
duckvep_register_sql_parts(duckdb_connection connection,
	const char *const *parts, size_t part_count)
{
	duckdb_result result;
	duckdb_state state;
	char *sql;
	size_t index, length, offset;

	length = 0;
	for (index = 0; index < part_count; index++) {
		size_t part_length;

		part_length = strlen(parts[index]);
		if (part_length > SIZE_MAX - length)
			return false;
		length += part_length;
	}
	if (length == SIZE_MAX)
		return false;
	sql = malloc(length + 1);
	if (sql == NULL)
		return false;
	offset = 0;
	for (index = 0; index < part_count; index++) {
		size_t part_length;

		part_length = strlen(parts[index]);
		memcpy(sql + offset, parts[index], part_length);
		offset += part_length;
	}
	sql[offset] = '\0';
	state = duckdb_query(connection, sql, &result);
	free(sql);
	if (state != DuckDBSuccess) {
		const char *error;

		error = duckdb_result_error(&result);
		fprintf(stderr, "[duckhts] failed DuckVEP SQL registration: %s\n",
		    error != NULL ? error : "unknown error");
		duckdb_destroy_result(&result);
		return false;
	}
	duckdb_destroy_result(&result);
	return true;
}

static void
duckvep_so_terms_bind(duckdb_bind_info info)
{
	duckdb_logical_type utinyint_type, ubigint_type, varchar_type;

	utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_bind_add_result_column(info, "bit_index", utinyint_type);
	duckdb_bind_add_result_column(info, "consequence_mask", ubigint_type);
	duckdb_bind_add_result_column(info, "consequence", varchar_type);
	duckdb_bind_add_result_column(info, "impact_code", utinyint_type);
	duckdb_bind_add_result_column(info, "impact", varchar_type);
	duckdb_bind_add_result_column(info, "severity_rank", utinyint_type);
	duckdb_bind_add_result_column(info, "evaluator_tier", utinyint_type);
	duckdb_destroy_logical_type(&utinyint_type);
	duckdb_destroy_logical_type(&ubigint_type);
	duckdb_destroy_logical_type(&varchar_type);
}

static void
duckvep_so_scan_destroy(void *data)
{
	duckdb_free(data);
}

static void
duckvep_so_terms_init(duckdb_init_info info)
{
	duckvep_so_scan_t *scan;

	scan = duckdb_malloc(sizeof(*scan));
	if (scan == NULL) {
		duckdb_init_set_error(info,
		    "duckvep_so_terms: failed to allocate scan state");
		return;
	}
	scan->offset = 0;
	duckdb_init_set_init_data(info, scan, duckvep_so_scan_destroy);
}

static void
duckvep_so_terms_scan(duckdb_function_info info, duckdb_data_chunk output)
{
	duckvep_so_scan_t *scan;
	duckdb_vector bit_vector, mask_vector, consequence_vector;
	duckdb_vector impact_code_vector, impact_vector, rank_vector, tier_vector;
	uint8_t *bits, *impact_codes, *ranks, *tiers;
	uint64_t *masks;
	idx_t count, vector_size;

	scan = duckdb_function_get_init_data(info);
	if (scan == NULL) {
		duckdb_data_chunk_set_size(output, 0);
		return;
	}
	bit_vector = duckdb_data_chunk_get_vector(output, 0);
	mask_vector = duckdb_data_chunk_get_vector(output, 1);
	consequence_vector = duckdb_data_chunk_get_vector(output, 2);
	impact_code_vector = duckdb_data_chunk_get_vector(output, 3);
	impact_vector = duckdb_data_chunk_get_vector(output, 4);
	rank_vector = duckdb_data_chunk_get_vector(output, 5);
	tier_vector = duckdb_data_chunk_get_vector(output, 6);
	bits = duckdb_vector_get_data(bit_vector);
	masks = duckdb_vector_get_data(mask_vector);
	impact_codes = duckdb_vector_get_data(impact_code_vector);
	ranks = duckdb_vector_get_data(rank_vector);
	tiers = duckdb_vector_get_data(tier_vector);
	vector_size = duckdb_vector_size();
	count = 0;
	while (count < vector_size && scan->offset < DUCKVEP_SO_BIT_COUNT) {
		duckvep_so_bit_t bit;
		duckvep_impact_t impact;
		const char *name;

		bit = (duckvep_so_bit_t)scan->offset;
		impact = duckvep_so_bit_impact(bit);
		name = duckvep_so_name(bit);
		bits[count] = (uint8_t)bit;
		masks[count] = DUCKVEP_SO(bit);
		impact_codes[count] = (uint8_t)impact;
		ranks[count] = duckvep_so_rank(bit);
		tiers[count] = duckvep_so_tier(bit);
		duckdb_vector_assign_string_element(consequence_vector, count,
		    name != NULL ? name : "");
		duckdb_vector_assign_string_element(impact_vector, count,
		    duckvep_impact_name(impact));
		count++;
		scan->offset++;
	}
	duckdb_data_chunk_set_size(output, count);
}

static bool
duckvep_register_so_terms(duckdb_connection connection)
{
	duckdb_table_function function;
	duckdb_state state;

	function = duckdb_create_table_function();
	duckdb_table_function_set_name(function, "duckvep_so_terms");
	duckdb_table_function_set_bind(function, duckvep_so_terms_bind);
	duckdb_table_function_set_init(function, duckvep_so_terms_init);
	duckdb_table_function_set_function(function, duckvep_so_terms_scan);
	state = duckdb_register_table_function(connection, function);
	duckdb_destroy_table_function(&function);
	return state == DuckDBSuccess;
}

static bool
duckvep_register_annotate_relation(duckdb_connection connection)
{
	static const char *const sql[] = {
		"CREATE OR REPLACE MACRO duckvep_annotate(events_table, model_name, ",
		"hgvs := false, upstream_distance := 5000, downstream_distance := 5000, ",
		"rich := false) AS TABLE ",
		"WITH parameters AS MATERIALIZED (SELECT CAST(model_name AS VARCHAR) AS model_name, ",
		"coalesce(CAST(hgvs AS BOOLEAN), false) AS include_hgvs, ",
		"coalesce(CAST(rich AS BOOLEAN), false) AS include_rich, ",
		"CAST(upstream_distance AS UBIGINT) AS upstream_distance, ",
		"CAST(downstream_distance AS UBIGINT) AS downstream_distance), ",
		"source AS (SELECT CAST(e.event_index AS UBIGINT) AS __duckvep_event_index, ",
		"CAST(e.seq_region AS UINTEGER) AS __duckvep_seq_region, ",
		"CAST(e.position AS UBIGINT) AS __duckvep_position, ",
		"CAST(e.reference AS VARCHAR) AS __duckvep_reference, ",
		"CAST(e.alternate AS VARCHAR) AS __duckvep_alternate, ",
		"CAST(e.end_position AS UBIGINT) AS __duckvep_end_position, ",
		"upper(CAST(e.structural_type AS VARCHAR)) AS __duckvep_explicit_structural_type, ",
		"upper(CAST(e.copy_change AS VARCHAR)) AS __duckvep_explicit_copy_change, ",
		"CAST(e.mate_seq_region AS UINTEGER) AS __duckvep_mate_seq_region, ",
		"CAST(e.mate_position AS UBIGINT) AS __duckvep_mate_position ",
		"FROM query_table(events_table) e), ",
		"typed_base AS MATERIALIZED (SELECT *, ",
		"CASE __duckvep_alternate ",
		"WHEN '<DEL>' THEN 'DEL' WHEN '<DUP>' THEN 'DUP' ",
		"WHEN '<TDUP>' THEN 'TDUP' WHEN '<STR>' THEN 'STR' ",
		"WHEN '<INV>' THEN 'INV' WHEN '<INS>' THEN 'INS' ",
		"WHEN '<CNV>' THEN 'CNV' WHEN '<UNKNOWN>' THEN 'UNKNOWN' ",
		"ELSE NULL END AS __duckvep_symbolic_structural_type, ",
		"__duckvep_alternate = '<*>' AS __duckvep_unspecified_alt, ",
		"__duckvep_alternate IN ('<NON_REF>', '*', '.') ",
		"AS __duckvep_non_variant FROM source), ",
		"typed AS MATERIALIZED (SELECT *, ",
		"coalesce(__duckvep_explicit_structural_type, ",
		"__duckvep_symbolic_structural_type) AS __duckvep_structural_type ",
		"FROM typed_base), ",
		"classified AS MATERIALIZED (SELECT *, CASE ",
		"WHEN __duckvep_non_variant THEN 'non_variant' ",
		"WHEN __duckvep_unspecified_alt THEN 'small_variant' ",
		"WHEN __duckvep_mate_seq_region IS NOT NULL OR __duckvep_mate_position IS NOT NULL ",
		"THEN 'breakend' WHEN __duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL ",
		"OR (__duckvep_alternate IS NOT NULL AND ",
		"(starts_with(__duckvep_alternate, '<') OR contains(__duckvep_alternate, '[') ",
		"OR contains(__duckvep_alternate, ']'))) ",
		"OR (__duckvep_end_position IS NOT NULL AND ",
		"(__duckvep_reference IS NULL OR __duckvep_alternate IS NULL)) ",
		"THEN 'structural_variant' ",
		"ELSE 'small_variant' END AS __duckvep_event_kind, ",
		"coalesce(__duckvep_explicit_copy_change, CASE __duckvep_structural_type ",
		"WHEN 'DEL' THEN 'LOSS' WHEN 'DUP' THEN 'GAIN' WHEN 'TDUP' THEN 'GAIN' ",
		"ELSE 'UNKNOWN' END) AS __duckvep_copy_change FROM typed), ",
		"validation AS MATERIALIZED (SELECT CASE ",
		"WHEN (SELECT model_name IS NULL OR model_name = '' FROM parameters) ",
		"THEN error('duckvep_annotate: model_name must be non-empty') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_index IS NULL) ",
		"THEN error('duckvep_annotate: event_index is required') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_seq_region IS NULL ",
		"OR __duckvep_position IS NULL OR __duckvep_position = 0) ",
		"THEN error('duckvep_annotate: seq_region and positive one-based position are required') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'non_variant' ",
		"AND (__duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL ",
		"OR __duckvep_mate_seq_region IS NOT NULL ",
		"OR __duckvep_mate_position IS NOT NULL)) ",
		"THEN error('duckvep_annotate: a non-variant gVCF allele cannot carry structural or breakend fields') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_unspecified_alt ",
		"AND (__duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL ",
		"OR __duckvep_mate_seq_region IS NOT NULL ",
		"OR __duckvep_mate_position IS NOT NULL)) ",
		"THEN error('duckvep_annotate: <*> cannot carry structural or breakend fields') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE ",
		"(__duckvep_mate_seq_region IS NULL) != (__duckvep_mate_position IS NULL)) ",
		"THEN error('duckvep_annotate: a breakend requires both mate_seq_region and mate_position') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'breakend' AND ",
		"(__duckvep_end_position IS NOT NULL OR __duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL)) ",
		"THEN error('duckvep_annotate: breakend mate coordinates cannot be mixed with single-locus structural fields') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'breakend' ",
		"AND __duckvep_mate_position = 0) ",
		"THEN error('duckvep_annotate: mate_position must be positive and one-based') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'small_variant' ",
		"AND (__duckvep_reference IS NULL OR __duckvep_reference = '' ",
		"OR __duckvep_alternate IS NULL OR __duckvep_alternate = '')) ",
		"THEN error('duckvep_annotate: literal small variants require reference and alternate') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE ",
		"__duckvep_explicit_structural_type IS NOT NULL ",
		"AND __duckvep_symbolic_structural_type IS NOT NULL ",
		"AND __duckvep_explicit_structural_type <> __duckvep_symbolic_structural_type) ",
		"THEN error('duckvep_annotate: symbolic alternate and structural_type disagree') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'structural_variant' ",
		"AND (__duckvep_structural_type IS NULL OR __duckvep_end_position IS NULL)) ",
		"THEN error('duckvep_annotate: a structural event requires end_position and a supported structural_type') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'structural_variant' ",
		"AND __duckvep_structural_type NOT IN ('DEL','DUP','TDUP','STR','INV','INS','CNV','UNKNOWN')) ",
		"THEN error('duckvep_annotate: unsupported structural_type') ",
		"ELSE true END AS valid), ",
		"validated AS (SELECT classified.* FROM classified CROSS JOIN validation ",
		"WHERE validation.valid), ",
		"small_events AS (SELECT * FROM validated ",
		"WHERE __duckvep_event_kind = 'small_variant'), ",
		"structural_events AS (SELECT * FROM validated ",
		"WHERE __duckvep_event_kind = 'structural_variant'), ",
		"breakend_events AS (SELECT * FROM validated ",
		"WHERE __duckvep_event_kind = 'breakend'), ",
		"small_compact_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_compact(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE (NOT p.include_hgvs OR __duckvep_unspecified_alt) ",
		"AND NOT p.include_rich), ",
		"small_hgvs_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_hgvs(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE p.include_hgvs AND NOT __duckvep_unspecified_alt ",
		"AND NOT p.include_rich), ",
		"small_rich_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_rich(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE (NOT p.include_hgvs OR __duckvep_unspecified_alt) ",
		"AND p.include_rich), ",
		"small_rich_hgvs_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_rich_hgvs(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE p.include_hgvs AND NOT __duckvep_unspecified_alt ",
		"AND p.include_rich), ",
		"structural_compact_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_structural_compact(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_end_position, __duckvep_structural_type, ",
		"__duckvep_copy_change, p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM structural_events e CROSS JOIN parameters p WHERE NOT p.include_rich), ",
		"structural_rich_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_structural_rich(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_end_position, __duckvep_structural_type, ",
		"__duckvep_copy_change, p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM structural_events e CROSS JOIN parameters p WHERE p.include_rich), ",
		"breakend_compact_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_breakend_compact(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_mate_seq_region, __duckvep_mate_position, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM breakend_events e CROSS JOIN parameters p WHERE NOT p.include_rich), ",
		"breakend_rich_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_breakend_rich(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_mate_seq_region, __duckvep_mate_position, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM breakend_events e CROSS JOIN parameters p WHERE p.include_rich), ",
		"small_compact AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* ",
		"FROM small_compact_raw), ",
		"small_hgvs AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* FROM small_hgvs_raw), ",
		"small_rich AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason FROM small_rich_raw), ",
		"small_rich_hgvs AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason ",
		"FROM small_rich_hgvs_raw), ",
		"structural_compact AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* ",
		"FROM structural_compact_raw), ",
		"structural_rich AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason ",
		"FROM structural_rich_raw), ",
		"breakend_compact AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* ",
		"FROM breakend_compact_raw), ",
		"breakend_rich AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason ",
		"FROM breakend_rich_raw) ",
		"SELECT * FROM small_compact UNION ALL BY NAME SELECT * FROM small_hgvs ",
		"UNION ALL BY NAME SELECT * FROM small_rich ",
		"UNION ALL BY NAME SELECT * FROM small_rich_hgvs ",
		"UNION ALL BY NAME SELECT * FROM structural_compact ",
		"UNION ALL BY NAME SELECT * FROM structural_rich ",
		"UNION ALL BY NAME SELECT * FROM breakend_compact ",
		"UNION ALL BY NAME SELECT * FROM breakend_rich"
	};

	return duckvep_register_sql_parts(connection, sql,
	    sizeof(sql) / sizeof(sql[0]));
}

/* SQL presentation uses the kernel's immutable genetic-code authority. The
 * string is indexed by T/C/A/G at each of the three codon positions. */
static void
duckvep_projection_code(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	duckdb_vector vector = duckdb_data_chunk_get_vector(input, 0);
	const uint8_t *codes = duckdb_vector_get_data(vector);
	uint64_t *validity = duckdb_vector_get_validity(vector);
	idx_t count = duckdb_data_chunk_get_size(input);

	duckdb_vector_ensure_validity_writable(output);
	for (idx_t row = 0; row < count; row++) {
		const char *amino_acids;

		if (validity != NULL && !duckdb_validity_row_is_valid(validity, row)) {
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
			continue;
		}
		amino_acids = duckvep_codon_table_amino_acids((duckvep_codon_table_t)codes[row]);
		if (amino_acids == NULL) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_transcript_projection: unsupported genetic code");
			return;
		}
		duckdb_vector_assign_string_element_len(output, row, amino_acids, 64);
	}
}

static bool
duckvep_register_projection_code(duckdb_connection connection)
{
	duckdb_scalar_function function = duckdb_create_scalar_function();
	duckdb_logical_type code_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	duckdb_logical_type text_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_state state;

	duckdb_scalar_function_set_name(function, "__duckvep_projection_code");
	duckdb_scalar_function_add_parameter(function, code_type);
	duckdb_scalar_function_set_return_type(function, text_type);
	duckdb_scalar_function_set_function(function, duckvep_projection_code);
	state = duckdb_register_scalar_function(connection, function);
	duckdb_destroy_scalar_function(&function);
	duckdb_destroy_logical_type(&code_type);
	duckdb_destroy_logical_type(&text_type);
	return state == DuckDBSuccess;
}

static bool
duckvep_register_projection_relation(duckdb_connection connection)
{
	static const char *const sql[] = {
		"CREATE OR REPLACE MACRO __duckvep_projection_base(exons, strand, pos1) AS\n",
		"  list_transform(list_filter(exons, e -> pos1 BETWEEN e.exon_start AND e.exon_end),\n",
		"    e -> e.exon_cdna_start::BIGINT + CASE WHEN strand > 0\n",
		"      THEN pos1::BIGINT - e.exon_start::BIGINT\n",
		"      ELSE e.exon_end::BIGINT - pos1::BIGINT END)[1];\n",
		"\n",
		"CREATE OR REPLACE MACRO __duckvep_projection_peptide(dna, genetic_code) AS (\n",
		"  WITH triplets AS (\n",
		"    SELECT i, substring(upper(dna), 3*i+1, 3) AS triplet\n",
		"    FROM range(length(dna) // 3) r(i)\n",
		"  ), amino_acids AS (\n",
		"    SELECT i, CASE WHEN count(DISTINCT substring(genetic_code, j+1, 1)) = 1\n",
		"      THEN min(substring(genetic_code, j+1, 1)) ELSE 'X' END AS aa\n",
		"    FROM triplets LEFT JOIN range(64) r(j) ON\n",
		"      (triplet[1] = 'N' OR triplet[1] = substring('TCAG', j // 16 + 1, 1)) AND\n",
		"      (triplet[2] = 'N' OR triplet[2] = substring('TCAG', (j // 4) % 4 + 1, 1)) AND\n",
		"      (triplet[3] = 'N' OR triplet[3] = substring('TCAG', j % 4 + 1, 1))\n",
		"    GROUP BY i\n",
		"  ), peptide AS (\n",
		"    SELECT coalesce(string_agg(aa, '' ORDER BY i), '') AS aa FROM amino_acids\n",
		"  ) SELECT CASE WHEN dna IS NULL THEN NULL ELSE\n",
		"    coalesce(nullif(aa || CASE WHEN length(dna) % 3 != 0 AND aa != '*'\n",
		"      THEN 'X' ELSE '' END, ''), '-') END FROM peptide\n",
		");\n",
		"\n",
		"CREATE OR REPLACE MACRO duckvep_transcript_projection(\n",
		"  events_table, annotations_table, transcripts_table\n",
		") AS TABLE\n",
		"WITH events AS MATERIALIZED (\n",
		"  SELECT event_index::UBIGINT AS event_index, seq_region::UINTEGER AS seq_region,\n",
		"    position::UBIGINT AS position, reference::VARCHAR AS reference,\n",
		"    alternate::VARCHAR AS alternate FROM query_table(events_table)\n",
		"), annotations AS MATERIALIZED (\n",
		"  SELECT event_index::UBIGINT AS event_index, transcript_index::UINTEGER AS transcript_index,\n",
		"    consequence_mask::UBIGINT AS consequence_mask FROM query_table(annotations_table)\n",
		"), transcripts AS MATERIALIZED (\n",
		"  SELECT transcript_index, seq_region, transcript_start, transcript_end, cds_start, cds_end,\n",
		"    strand, exons, transcript_flags, peptide_edits, cds_sequence, post_cds_sequence, codon_table\n",
		"  FROM query_table(transcripts_table)\n",
		"), validation AS MATERIALIZED (\n",
		"  SELECT CASE\n",
		"    WHEN EXISTS (SELECT 1 FROM events WHERE event_index IS NULL OR seq_region IS NULL\n",
		"      OR position IS NULL OR position = 0 OR reference IS NULL OR alternate IS NULL)\n",
		"      THEN error('duckvep_transcript_projection: event keys, position and literal alleles are required')\n",
		"    WHEN EXISTS (SELECT event_index FROM events GROUP BY event_index HAVING count(*) != 1)\n",
		"      THEN error('duckvep_transcript_projection: event_index must be unique')\n",
		"    WHEN EXISTS (SELECT transcript_index FROM transcripts GROUP BY transcript_index\n",
		"      HAVING count(*) != 1 OR transcript_index IS NULL)\n",
		"      THEN error('duckvep_transcript_projection: transcript_index must be non-null and unique')\n",
		"    WHEN EXISTS (SELECT 1 FROM annotations a LEFT JOIN events e USING(event_index)\n",
		"      WHERE e.event_index IS NULL OR a.consequence_mask IS NULL)\n",
		"      THEN error('duckvep_transcript_projection: annotation needs a source event and consequence mask')\n",
		"    WHEN EXISTS (SELECT 1 FROM annotations a JOIN events e USING(event_index)\n",
		"      LEFT JOIN transcripts t USING(transcript_index)\n",
		"      WHERE a.transcript_index IS NOT NULL AND\n",
		"        (t.transcript_index IS NULL OR e.seq_region IS DISTINCT FROM t.seq_region))\n",
		"      THEN error('duckvep_transcript_projection: annotation transcript and source region must match the model')\n",
		"    ELSE true END AS valid\n",
		"), joined AS MATERIALIZED (\n",
		"  SELECT a.event_index, a.transcript_index, a.consequence_mask,\n",
		"    e.reference, e.alternate, e.position,\n",
		"    t.transcript_start::BIGINT AS tx_start, t.transcript_end::BIGINT AS tx_end,\n",
		"    t.cds_start::BIGINT AS tx_cds_start, t.cds_end::BIGINT AS tx_cds_end,\n",
		"    t.strand, t.exons, t.transcript_flags, t.peptide_edits,\n",
		"    upper(decode(t.cds_sequence)) AS cds_sequence,\n",
		"    upper(decode(t.post_cds_sequence)) AS post_cds_sequence,\n",
		"    __duckvep_projection_code(t.codon_table) AS genetic_code,\n",
		"    duckvep_allele_geometry(e.position, e.reference, e.alternate) AS geometry\n",
		"  FROM annotations a JOIN events e USING (event_index)\n",
		"  LEFT JOIN transcripts t USING (transcript_index)\n",
		"  CROSS JOIN validation WHERE validation.valid\n",
		"), features AS (\n",
		"  SELECT *, geometry.feature_start0::BIGINT + 1 AS vf_start,\n",
		"    geometry.feature_end0::BIGINT AS vf_end,\n",
		"    upper(CASE WHEN length(reference) = length(alternate) THEN reference\n",
		"      ELSE substring(reference, geometry.reference_difference_offset + 1,\n",
		"        geometry.reference_difference_length) END) AS feature_reference,\n",
		"    upper(CASE WHEN length(reference) = length(alternate) THEN alternate\n",
		"      ELSE substring(alternate, geometry.alternate_difference_offset + 1,\n",
		"        geometry.alternate_difference_length) END) AS feature_alternate,\n",
		"    __duckvep_projection_base(exons, strand,\n",
		"      CASE WHEN strand > 0 THEN tx_cds_start ELSE tx_cds_end END) AS coding_cdna_start,\n",
		"    __duckvep_projection_base(exons, strand,\n",
		"      CASE WHEN strand > 0 THEN tx_cds_end ELSE tx_cds_start END) AS coding_cdna_end,\n",
		/* VEP BaseTranscriptVariation and TranscriptMapper use the first
		 * transcript exon here. CDS sequence padding separately uses the
		 * translation-start exon; substituting that phase changes VEP output. */
		"    greatest(exons[1].phase::BIGINT, 0) AS phase_offset\n",
		"  FROM joined\n",
		"), mapped AS (\n",
		"  SELECT *, vf_end >= tx_start AND vf_start <= tx_end AS within_transcript,\n",
		"    __duckvep_projection_base(exons, strand,\n",
		"      CASE WHEN strand > 0 THEN vf_start ELSE vf_end END) AS first_cdna,\n",
		"    __duckvep_projection_base(exons, strand,\n",
		"      CASE WHEN strand > 0 THEN vf_end ELSE vf_start END) AS last_cdna,\n",
		"    list_filter(range(1, length(exons)+1), i ->\n",
		"      vf_start <= exons[i].exon_end AND vf_end >= exons[i].exon_start) AS exon_hits,\n",
		"    list_filter(range(1, length(exons)), i ->\n",
		"      vf_start <= greatest(exons[i].exon_start, exons[i+1].exon_start)::BIGINT - 1 AND\n",
		"      vf_end >= least(exons[i].exon_end, exons[i+1].exon_end)::BIGINT + 1) AS intron_hits,\n",
		"    CASE WHEN strand > 0 THEN feature_alternate ELSE seq_revcomp(feature_alternate) END\n",
		"      AS transcript_alternate\n",
		"  FROM features\n",
		"), cdna AS (\n",
		"  SELECT *, CASE WHEN within_transcript THEN CASE WHEN vf_start > vf_end\n",
		"      THEN coalesce(first_cdna, last_cdna + 1) ELSE first_cdna END END AS raw_cdna_start,\n",
		"    CASE WHEN within_transcript THEN CASE WHEN vf_start > vf_end\n",
		"      THEN coalesce(last_cdna, first_cdna - 1) ELSE last_cdna END END AS raw_cdna_end\n",
		"  FROM mapped\n",
		"), coding AS (\n",
		"  SELECT *, CASE WHEN raw_cdna_start BETWEEN coding_cdna_start AND coding_cdna_end\n",
		"      AND NOT (vf_start > vf_end AND (raw_cdna_start > coding_cdna_end OR raw_cdna_end < coding_cdna_start))\n",
		"      THEN raw_cdna_start - coding_cdna_start + phase_offset + 1 END AS raw_cds_start,\n",
		"    CASE WHEN raw_cdna_end BETWEEN coding_cdna_start AND coding_cdna_end\n",
		"      AND NOT (vf_start > vf_end AND (raw_cdna_start > coding_cdna_end OR raw_cdna_end < coding_cdna_start))\n",
		"      THEN raw_cdna_end - coding_cdna_start + phase_offset + 1 END AS raw_cds_end\n",
		"  FROM cdna\n",
		"), protein AS (\n",
		"  SELECT *, (raw_cds_start + 2) // 3 AS raw_protein_start,\n",
		"    (raw_cds_end + 2) // 3 AS raw_protein_end,\n",
		"    substring(cds_sequence, 1, raw_cds_start - 1) || transcript_alternate ||\n",
		"      substring(cds_sequence, raw_cds_end + 1) AS alternate_cds\n",
		"  FROM coding\n",
		"), codons AS (\n",
		"  SELECT *, CASE WHEN raw_cds_start IS NOT NULL AND raw_cds_end IS NOT NULL\n",
		"    THEN substring(cds_sequence, raw_protein_start * 3 - 2,\n",
		"      greatest((raw_protein_end - raw_protein_start + 1) * 3, 0)) END AS ref_codons,\n",
		"    CASE WHEN raw_cds_start IS NOT NULL AND raw_cds_end IS NOT NULL\n",
		"    THEN substring((CASE WHEN length(alternate_cds) >= 3 THEN alternate_cds ELSE '' END)\n",
		"      || coalesce(post_cds_sequence, ''), raw_protein_start * 3 - 2,\n",
		"      greatest((raw_protein_end - raw_protein_start + 1) * 3 + length(feature_alternate)\n",
		"        - (raw_cds_end - raw_cds_start + 1), 0)) END AS alt_codons,\n",
		"    (raw_cds_start - 1) % 3 AS changed_offset\n",
		"  FROM protein\n",
		"), translated_raw AS (\n",
		"  SELECT *, CASE WHEN regexp_full_match(feature_reference, '[ACGT]*')\n",
		"      THEN __duckvep_projection_peptide(ref_codons, genetic_code) END AS ref_peptide,\n",
		"    CASE WHEN regexp_full_match(feature_alternate, '[ACGT]*')\n",
		"      THEN __duckvep_projection_peptide(alt_codons, genetic_code) END AS alt_peptide\n",
		"  FROM codons\n",
		"), translated AS (\n",
		"  SELECT * EXCLUDE(ref_peptide), CASE WHEN ref_peptide IS NULL OR ref_peptide = '-'\n",
		"      THEN ref_peptide ELSE array_to_string(list_transform(range(length(ref_peptide)), i ->\n",
		"        coalesce(list_transform(list_filter(peptide_edits, edit ->\n",
		"          edit.protein_position = least(raw_protein_start, raw_protein_end) + i),\n",
		"          edit -> edit.alternate_amino_acid)[1], substring(ref_peptide, i+1, 1))), '')\n",
		"    END AS ref_peptide\n",
		"  FROM translated_raw\n",
		")\n",
		"SELECT event_index, transcript_index,\n",
		"  coalesce(nullif(feature_alternate, ''), '-') AS output_allele,\n",
		"  geometry.interbase AS interbase,\n",
		"  CASE WHEN raw_cdna_start > raw_cdna_end THEN raw_cdna_end ELSE raw_cdna_start END AS cdna_start,\n",
		"  CASE WHEN raw_cdna_start > raw_cdna_end THEN raw_cdna_start ELSE raw_cdna_end END AS cdna_end,\n",
		"  CASE WHEN raw_cds_start > raw_cds_end THEN raw_cds_end ELSE raw_cds_start END AS cds_start,\n",
		"  CASE WHEN raw_cds_start > raw_cds_end THEN raw_cds_start ELSE raw_cds_end END AS cds_end,\n",
		"  CASE WHEN raw_protein_start > raw_protein_end THEN raw_protein_end ELSE raw_protein_start END AS protein_start,\n",
		"  CASE WHEN raw_protein_start > raw_protein_end THEN raw_protein_start ELSE raw_protein_end END AS protein_end,\n",
		"  exon_hits[1]::UINTEGER AS exon_first, exon_hits[-1]::UINTEGER AS exon_last,\n",
		"  length(exons)::UINTEGER AS exon_total,\n",
		"  intron_hits[1]::UINTEGER AS intron_first, intron_hits[-1]::UINTEGER AS intron_last,\n",
		"  CASE WHEN transcript_index IS NOT NULL THEN greatest(length(exons) - 1, 0)::UINTEGER END AS intron_total,\n",
		"  CASE WHEN (consequence_mask & (SELECT bit_or(consequence_mask) FROM duckvep_so_terms()\n",
		"    WHERE consequence IN ('upstream_gene_variant', 'downstream_gene_variant'))) != 0\n",
		"    THEN least(abs(vf_start - tx_start), abs(vf_start - tx_end),\n",
		"      abs(vf_end - tx_start), abs(vf_end - tx_end)) END AS transcript_distance,\n",
		"  (transcript_flags & 16) != 0 AS cds_start_nf,\n",
		"  (transcript_flags & 32) != 0 AS cds_end_nf,\n",
		"  CASE WHEN alt_peptide IS NOT NULL THEN ref_peptide END AS reference_amino_acids,\n",
		"  CASE WHEN ref_peptide IS NOT NULL THEN alt_peptide END AS alternate_amino_acids,\n",
		"  CASE WHEN ref_codons = '' THEN '-' ELSE\n",
		"    lower(substring(ref_codons, 1, changed_offset)) ||\n",
		"    substring(ref_codons, changed_offset + 1, length(feature_reference)) ||\n",
		"    lower(substring(ref_codons, changed_offset + length(feature_reference) + 1)) END\n",
		"    AS reference_codons,\n",
		"  CASE WHEN alt_codons = '' THEN '-' ELSE\n",
		"    lower(substring(alt_codons, 1, changed_offset)) ||\n",
		"    substring(alt_codons, changed_offset + 1, length(feature_alternate)) ||\n",
		"    lower(substring(alt_codons, changed_offset + length(feature_alternate) + 1)) END\n",
		"    AS alternate_codons\n",
		"FROM translated;\n"
	};

	return duckvep_register_sql_parts(connection, sql,
	    sizeof(sql) / sizeof(sql[0]));
}

bool
register_duckvep_sql_functions(duckdb_connection connection)
{
	return duckvep_register_so_terms(connection) &&
	    duckvep_register_projection_code(connection) &&
	    duckvep_register_annotate_relation(connection) &&
	    duckvep_register_projection_relation(connection);
}
