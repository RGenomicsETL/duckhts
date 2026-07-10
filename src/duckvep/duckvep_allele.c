/*
 * Private scalar allele and edit primitives for DuckVEP.
 */
#include "duckvep_allele.h"

#include <limits.h>

static int
duckvep_bytes_view_valid(duckvep_bytes_view_t view)
{
	return view.length == 0 || view.data != NULL;
}

static duckvep_allele_status_t
duckvep_size_to_u64(size_t value, uint64_t *result)
{
	uint64_t narrowed;

	if (result == NULL)
		return DUCKVEP_ALLELE_INVALID;
	narrowed = (uint64_t)value;
	if ((size_t)narrowed != value)
		return DUCKVEP_ALLELE_RANGE;
	*result = narrowed;
	return DUCKVEP_ALLELE_OK;
}

duckvep_allele_status_t
duckvep_span0_make(uint64_t begin, uint64_t end, duckvep_span0_t *result)
{
	if (result == NULL || begin > end)
		return DUCKVEP_ALLELE_INVALID;
	result->begin = begin;
	result->end = end;
	return DUCKVEP_ALLELE_OK;
}

duckvep_allele_status_t
duckvep_span0_from_length(uint64_t begin, uint64_t length,
	duckvep_span0_t *result)
{
	if (result == NULL)
		return DUCKVEP_ALLELE_INVALID;
	if (length > UINT64_MAX - begin)
		return DUCKVEP_ALLELE_RANGE;
	result->begin = begin;
	result->end = begin + length;
	return DUCKVEP_ALLELE_OK;
}

duckvep_allele_status_t
duckvep_span0_length(const duckvep_span0_t *span, uint64_t *result)
{
	if (span == NULL || result == NULL || span->begin > span->end)
		return DUCKVEP_ALLELE_INVALID;
	*result = span->end - span->begin;
	return DUCKVEP_ALLELE_OK;
}

duckvep_allele_status_t
duckvep_bytes_view_make(const void *data, size_t length,
	duckvep_bytes_view_t *result)
{
	if (result == NULL || (length != 0 && data == NULL))
		return DUCKVEP_ALLELE_INVALID;
	result->data = (const uint8_t *)data;
	result->length = length;
	return DUCKVEP_ALLELE_OK;
}

duckvep_allele_status_t
duckvep_edit_view_make(duckvep_span0_t ref_span,
	duckvep_bytes_view_t ref, duckvep_bytes_view_t alt,
	duckvep_edit_view_t *result)
{
	uint64_t span_length, ref_length, alt_length;
	duckvep_allele_status_t status;

	if (result == NULL || !duckvep_bytes_view_valid(ref) ||
	    !duckvep_bytes_view_valid(alt))
		return DUCKVEP_ALLELE_INVALID;
	status = duckvep_span0_length(&ref_span, &span_length);
	if (status != DUCKVEP_ALLELE_OK)
		return status;
	status = duckvep_size_to_u64(ref.length, &ref_length);
	if (status != DUCKVEP_ALLELE_OK)
		return status;
	status = duckvep_size_to_u64(alt.length, &alt_length);
	if (status != DUCKVEP_ALLELE_OK)
		return status;
	(void)alt_length;
	if (span_length != ref_length)
		return DUCKVEP_ALLELE_INVALID;
	result->ref_span = ref_span;
	result->ref = ref;
	result->alt = alt;
	return DUCKVEP_ALLELE_OK;
}

static duckvep_bytes_view_t
duckvep_bytes_view_slice(duckvep_bytes_view_t source, size_t offset,
	size_t length)
{
	duckvep_bytes_view_t result;

	result.data = source.data == NULL ? NULL : source.data + offset;
	result.length = length;
	return result;
}

duckvep_allele_status_t
duckvep_semantic_edit_trim(const duckvep_edit_view_t *raw,
	duckvep_semantic_edit_view_t *result)
{
	duckvep_edit_view_t checked;
	duckvep_span0_t semantic_span;
	size_t prefix, suffix, ref_remaining, alt_remaining;
	duckvep_allele_status_t status;

	if (raw == NULL || result == NULL)
		return DUCKVEP_ALLELE_INVALID;
	status = duckvep_edit_view_make(raw->ref_span, raw->ref, raw->alt,
	    &checked);
	if (status != DUCKVEP_ALLELE_OK)
		return status;

	prefix = 0;
	while (prefix < checked.ref.length && prefix < checked.alt.length &&
	    checked.ref.data[prefix] == checked.alt.data[prefix])
		prefix++;

	ref_remaining = checked.ref.length - prefix;
	alt_remaining = checked.alt.length - prefix;
	suffix = 0;
	while (suffix < ref_remaining && suffix < alt_remaining &&
	    checked.ref.data[checked.ref.length - 1 - suffix] ==
	    checked.alt.data[checked.alt.length - 1 - suffix])
		suffix++;

	status = duckvep_span0_make(
	    checked.ref_span.begin + (uint64_t)prefix,
	    checked.ref_span.end - (uint64_t)suffix, &semantic_span);
	if (status != DUCKVEP_ALLELE_OK)
		return status;
	result->edit.ref_span = semantic_span;
	result->edit.ref = duckvep_bytes_view_slice(checked.ref, prefix,
	    checked.ref.length - prefix - suffix);
	result->edit.alt = duckvep_bytes_view_slice(checked.alt, prefix,
	    checked.alt.length - prefix - suffix);
	result->common_prefix = prefix;
	result->common_suffix = suffix;
	return DUCKVEP_ALLELE_OK;
}

/*
 * Independent C rewrite of
 * Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele::
 * _get_differing_regions() from Ensembl Variation release 116 commit
 * 2fb834b987ede3824e200197a838ce11e91aeb4b.  That module is distributed by
 * VEP release/116.0 commit 57ea5c52340acc1f156267f810ad162e26597082;
 * the pinned installed module's SHA-256 is
 * 493dd598a8f582d31f3df95bdef44a73e98eea8ee183dfaecd36490b6ad2ebea.
 * No upstream implementation text is copied.
 *
 * Perl string XOR extends to the longer operand and treats absent bytes as
 * zero.  VEP bypasses XOR when either allele length is <= 1, and its insertion
 * sentinel {s => 0, e => -1} converts to the empty half-open run [0,0).
 */
static int
duckvep_byte_xor_differs(duckvep_bytes_view_t ref,
	duckvep_bytes_view_t alt, size_t offset)
{
	uint8_t ref_byte, alt_byte;

	ref_byte = offset < ref.length ? ref.data[offset] : 0;
	alt_byte = offset < alt.length ? alt.data[offset] : 0;
	return (uint8_t)(ref_byte ^ alt_byte) != 0;
}

static duckvep_allele_status_t
duckvep_byte_xor_diff_run_count(duckvep_bytes_view_t ref,
	duckvep_bytes_view_t alt, size_t *result)
{
	size_t i, limit, count;
	int in_run;
	uint64_t ignored;

	if (result == NULL || !duckvep_bytes_view_valid(ref) ||
	    !duckvep_bytes_view_valid(alt))
		return DUCKVEP_ALLELE_INVALID;
	limit = ref.length > alt.length ? ref.length : alt.length;
	if (duckvep_size_to_u64(limit, &ignored) != DUCKVEP_ALLELE_OK)
		return DUCKVEP_ALLELE_RANGE;

	if (ref.length <= 1 || alt.length <= 1) {
		*result = 1;
		return DUCKVEP_ALLELE_OK;
	}

	count = 0;
	in_run = 0;
	for (i = 0; i < limit; i++) {
		if (duckvep_byte_xor_differs(ref, alt, i)) {
			if (!in_run) {
				count++;
				in_run = 1;
			}
		} else {
			in_run = 0;
		}
	}
	*result = count;
	return DUCKVEP_ALLELE_OK;
}

static void
duckvep_vep_diff_view_clear(duckvep_vep_diff_view_t *view)
{
	view->algorithm = DUCKVEP_DIFF_ALGORITHM_NONE;
	view->reference.data = NULL;
	view->reference.length = 0;
	view->alternate.data = NULL;
	view->alternate.length = 0;
	view->runs = NULL;
	view->run_count = 0;
}

duckvep_allele_status_t
duckvep_diff_view_build(duckvep_diff_algorithm_t algorithm,
	duckvep_bytes_view_t ref, duckvep_bytes_view_t alt,
	duckvep_span0_t *run_storage, size_t run_capacity,
	duckvep_vep_diff_view_t *result, size_t *required_runs)
{
	duckvep_allele_status_t status;
	size_t required, i, run_start, run_index, limit;
	uint64_t begin, end, ref_length;
	int in_run;

	if (result == NULL || required_runs == NULL)
		return DUCKVEP_ALLELE_INVALID;
	duckvep_vep_diff_view_clear(result);
	*required_runs = 0;
	if (algorithm != DUCKVEP_DIFF_VEP_116_BYTE_XOR)
		return DUCKVEP_ALLELE_INVALID;
	status = duckvep_byte_xor_diff_run_count(ref, alt, &required);
	if (status != DUCKVEP_ALLELE_OK)
		return status;
	*required_runs = required;
	if (required > run_capacity)
		return DUCKVEP_ALLELE_NO_SPACE;
	if (required != 0 && run_storage == NULL)
		return DUCKVEP_ALLELE_INVALID;

	if (ref.length <= 1 || alt.length <= 1) {
		status = duckvep_size_to_u64(ref.length, &ref_length);
		if (status != DUCKVEP_ALLELE_OK)
			return status;
		run_storage[0].begin = 0;
		run_storage[0].end = ref_length;
	} else {
		limit = ref.length > alt.length ? ref.length : alt.length;
		run_start = 0;
		run_index = 0;
		in_run = 0;
		for (i = 0; i < limit; i++) {
			if (duckvep_byte_xor_differs(ref, alt, i)) {
				if (!in_run) {
					run_start = i;
					in_run = 1;
				}
				continue;
			}
			if (!in_run)
				continue;
			status = duckvep_size_to_u64(run_start, &begin);
			if (status != DUCKVEP_ALLELE_OK)
				return status;
			status = duckvep_size_to_u64(i, &end);
			if (status != DUCKVEP_ALLELE_OK)
				return status;
			run_storage[run_index].begin = begin;
			run_storage[run_index].end = end;
			run_index++;
			in_run = 0;
		}
		if (in_run) {
			status = duckvep_size_to_u64(run_start, &begin);
			if (status != DUCKVEP_ALLELE_OK)
				return status;
			status = duckvep_size_to_u64(limit, &end);
			if (status != DUCKVEP_ALLELE_OK)
				return status;
			run_storage[run_index].begin = begin;
			run_storage[run_index].end = end;
			run_index++;
		}
		if (run_index != required)
			return DUCKVEP_ALLELE_INVALID;
	}

	result->algorithm = algorithm;
	result->reference = ref;
	result->alternate = alt;
	result->runs = required == 0 ? NULL : run_storage;
	result->run_count = required;
	return DUCKVEP_ALLELE_OK;
}
