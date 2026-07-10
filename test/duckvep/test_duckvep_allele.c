#include "duckvep_allele.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures;

#define CHECK(condition) check_true((condition), #condition, __LINE__)

static void
check_true(int condition, const char *expression, int line)
{
	if (condition)
		return;
	fprintf(stderr, "line %d: CHECK_FAILED: %s\n", line, expression);
	failures++;
}

static duckvep_bytes_view_t
text_view(const char *text)
{
	duckvep_bytes_view_t result;

	result.data = (const uint8_t *)text;
	result.length = strlen(text);
	return result;
}

#define DUCKVEP_FUZZ_SEED UINT64_C(0x6476637478743136)
#define DUCKVEP_FUZZ_CASES ((size_t)12000)

typedef struct fuzz_rng {
	uint64_t state;
} fuzz_rng_t;

typedef enum fuzz_shape {
	FUZZ_BOTH_EMPTY = 0,
	FUZZ_EMPTY_REF,
	FUZZ_EMPTY_ALT,
	FUZZ_REF_ONE,
	FUZZ_ALT_ONE,
	FUZZ_EQUAL_LENGTH,
	FUZZ_IDENTICAL,
	FUZZ_RETAINED_ISLANDS,
	FUZZ_PADDED_DELINS,
	FUZZ_LONG_INSERTION,
	FUZZ_LONG_DELETION,
	FUZZ_RANDOM_UNEQUAL,
	FUZZ_SHAPE_COUNT
} fuzz_shape_t;

typedef struct fuzz_alleles {
	uint8_t *ref;
	uint8_t *alt;
	size_t ref_length;
	size_t alt_length;
	size_t common_prefix;
	size_t common_suffix;
	fuzz_shape_t shape;
} fuzz_alleles_t;

typedef struct expected_runs {
	duckvep_span0_t *runs;
	size_t count;
} expected_runs_t;

typedef struct owned_bytes {
	uint8_t *data;
	size_t length;
} owned_bytes_t;

static uint64_t
fuzz_next(fuzz_rng_t *rng)
{
	uint64_t value;

	rng->state += UINT64_C(0x9e3779b97f4a7c15);
	value = rng->state;
	value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
	value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
	return value ^ (value >> 31);
}

static size_t
fuzz_bounded(fuzz_rng_t *rng, size_t limit)
{
	if (limit == 0)
		return 0;
	return (size_t)(fuzz_next(rng) % (uint64_t)limit);
}

static uint8_t
fuzz_byte(fuzz_rng_t *rng)
{
	static const uint8_t alphabet[] = {
		0, 'A', 'C', 'G', 'T', 'N', UINT8_C(0x7f), UINT8_C(0xff)
	};

	return alphabet[fuzz_bounded(rng, sizeof(alphabet))];
}

static uint8_t
fuzz_different_byte(fuzz_rng_t *rng, uint8_t value)
{
	uint8_t bit;

	bit = (uint8_t)(UINT8_C(1) << fuzz_bounded(rng, 8));
	return (uint8_t)(value ^ bit);
}

static const char *
fuzz_shape_name(fuzz_shape_t shape)
{
	switch ((unsigned int)shape) {
	case FUZZ_BOTH_EMPTY:
		return "both_empty";
	case FUZZ_EMPTY_REF:
		return "empty_ref";
	case FUZZ_EMPTY_ALT:
		return "empty_alt";
	case FUZZ_REF_ONE:
		return "ref_one";
	case FUZZ_ALT_ONE:
		return "alt_one";
	case FUZZ_EQUAL_LENGTH:
		return "equal_length";
	case FUZZ_IDENTICAL:
		return "identical";
	case FUZZ_RETAINED_ISLANDS:
		return "retained_islands";
	case FUZZ_PADDED_DELINS:
		return "padded_delins";
	case FUZZ_LONG_INSERTION:
		return "long_insertion";
	case FUZZ_LONG_DELETION:
		return "long_deletion";
	case FUZZ_RANDOM_UNEQUAL:
		return "random_unequal";
	default:
		return "invalid";
	}
}

static void
fuzz_alleles_destroy(fuzz_alleles_t *alleles)
{
	free(alleles->ref);
	free(alleles->alt);
	alleles->ref = NULL;
	alleles->alt = NULL;
}

static int
fuzz_alleles_build(fuzz_rng_t *rng, size_t case_index,
	fuzz_alleles_t *result)
{
	size_t i, first, second, ref_core, alt_core;

	memset(result, 0, sizeof(*result));
	result->shape = (fuzz_shape_t)(case_index % (size_t)FUZZ_SHAPE_COUNT);
	switch ((unsigned int)result->shape) {
	case FUZZ_BOTH_EMPTY:
		break;
	case FUZZ_EMPTY_REF:
		result->alt_length = 1 + fuzz_bounded(rng, 128);
		break;
	case FUZZ_EMPTY_ALT:
		result->ref_length = 1 + fuzz_bounded(rng, 128);
		break;
	case FUZZ_REF_ONE:
		result->ref_length = 1;
		result->alt_length = fuzz_bounded(rng, 129);
		break;
	case FUZZ_ALT_ONE:
		result->ref_length = 2 + fuzz_bounded(rng, 127);
		result->alt_length = 1;
		break;
	case FUZZ_EQUAL_LENGTH:
		result->ref_length = 2 + fuzz_bounded(rng, 255);
		result->alt_length = result->ref_length;
		break;
	case FUZZ_IDENTICAL:
		result->ref_length = 2 + fuzz_bounded(rng, 255);
		result->alt_length = result->ref_length;
		break;
	case FUZZ_RETAINED_ISLANDS:
		result->ref_length = 8 + fuzz_bounded(rng, 249);
		result->alt_length = result->ref_length;
		break;
	case FUZZ_PADDED_DELINS:
		result->common_prefix = 1 + fuzz_bounded(rng, 8);
		result->common_suffix = 1 + fuzz_bounded(rng, 8);
		ref_core = 1 + fuzz_bounded(rng, 64);
		alt_core = 1 + fuzz_bounded(rng, 64);
		result->ref_length = result->common_prefix + ref_core +
		    result->common_suffix;
		result->alt_length = result->common_prefix + alt_core +
		    result->common_suffix;
		break;
	case FUZZ_LONG_INSERTION:
		result->ref_length = fuzz_bounded(rng, 2);
		result->alt_length = 512 + fuzz_bounded(rng, 4097);
		break;
	case FUZZ_LONG_DELETION:
		result->ref_length = 512 + fuzz_bounded(rng, 4097);
		result->alt_length = fuzz_bounded(rng, 2);
		break;
	case FUZZ_RANDOM_UNEQUAL:
		result->ref_length = 2 + fuzz_bounded(rng, 4095);
		result->alt_length = 2 + fuzz_bounded(rng, 4095);
		if (result->ref_length == result->alt_length)
			result->alt_length = result->alt_length == 4096 ? 2 :
			    result->alt_length + 1;
		break;
	default:
		return 0;
	}

	if (result->ref_length != 0) {
		result->ref = (uint8_t *)malloc(result->ref_length);
		if (result->ref == NULL)
			goto fail;
	}
	if (result->alt_length != 0) {
		result->alt = (uint8_t *)malloc(result->alt_length);
		if (result->alt == NULL)
			goto fail;
	}
	for (i = 0; i < result->ref_length; i++)
		result->ref[i] = fuzz_byte(rng);
	for (i = 0; i < result->alt_length; i++)
		result->alt[i] = fuzz_byte(rng);

	if (result->shape == FUZZ_IDENTICAL) {
		memcpy(result->alt, result->ref, result->ref_length);
	} else if (result->shape == FUZZ_RETAINED_ISLANDS) {
		memcpy(result->alt, result->ref, result->ref_length);
		first = 1 + fuzz_bounded(rng, result->ref_length / 3 - 1);
		second = 2 * result->ref_length / 3;
		second += fuzz_bounded(rng, result->ref_length - 1 - second);
		result->alt[first] = fuzz_different_byte(rng, result->ref[first]);
		result->alt[second] = fuzz_different_byte(rng, result->ref[second]);
	} else if (result->shape == FUZZ_PADDED_DELINS) {
		memcpy(result->alt, result->ref, result->common_prefix);
		for (i = 0; i < result->common_suffix; i++)
			result->alt[result->alt_length - 1 - i] =
			    result->ref[result->ref_length - 1 - i];
		result->alt[result->common_prefix] = fuzz_different_byte(rng,
		    result->ref[result->common_prefix]);
	}
	return 1;

fail:
	fuzz_alleles_destroy(result);
	return 0;
}

static void
expected_runs_destroy(expected_runs_t *expected)
{
	free(expected->runs);
	expected->runs = NULL;
	expected->count = 0;
}

/*
 * Slow test oracle: materialize the complete Perl-XOR byte mask first, count
 * zero-to-nonzero transitions, then walk that stored mask into runs.  This is
 * deliberately separate from the production streaming state machine.
 */
static int
slow_expected_vep_runs(duckvep_bytes_view_t ref,
	duckvep_bytes_view_t alt, expected_runs_t *result)
{
	uint8_t *mask;
	size_t i, limit, count, run_index, start;
	uint8_t ref_byte, alt_byte;

	result->runs = NULL;
	result->count = 0;
	if ((ref.length != 0 && ref.data == NULL) ||
	    (alt.length != 0 && alt.data == NULL))
		return 0;
	if (ref.length <= 1 || alt.length <= 1) {
		result->runs = (duckvep_span0_t *)malloc(sizeof(*result->runs));
		if (result->runs == NULL)
			return 0;
		result->runs[0].begin = 0;
		result->runs[0].end = (uint64_t)ref.length;
		result->count = 1;
		return 1;
	}

	limit = ref.length > alt.length ? ref.length : alt.length;
	mask = (uint8_t *)malloc(limit);
	if (mask == NULL)
		return 0;
	for (i = 0; i < limit; i++) {
		ref_byte = i < ref.length ? ref.data[i] : 0;
		alt_byte = i < alt.length ? alt.data[i] : 0;
		mask[i] = (uint8_t)(ref_byte ^ alt_byte);
	}

	count = 0;
	for (i = 0; i < limit; i++) {
		if (mask[i] != 0 && (i == 0 || mask[i - 1] == 0))
			count++;
	}
	if (count != 0) {
		result->runs = (duckvep_span0_t *)malloc(
		    count * sizeof(*result->runs));
		if (result->runs == NULL) {
			free(mask);
			return 0;
		}
	}

	i = 0;
	run_index = 0;
	while (i < limit) {
		while (i < limit && mask[i] == 0)
			i++;
		if (i == limit)
			break;
		start = i;
		while (i < limit && mask[i] != 0)
			i++;
		result->runs[run_index].begin = (uint64_t)start;
		result->runs[run_index].end = (uint64_t)i;
		run_index++;
	}
	free(mask);
	if (run_index != count) {
		expected_runs_destroy(result);
		return 0;
	}
	result->count = count;
	return 1;
}

static void
owned_bytes_destroy(owned_bytes_t *bytes)
{
	free(bytes->data);
	bytes->data = NULL;
	bytes->length = 0;
}

/* Independent whole-context splice used only as the trim-equivalence oracle. */
static int
slow_apply_edit(duckvep_bytes_view_t reference,
	const duckvep_edit_view_t *edit, owned_bytes_t *result)
{
	size_t begin, end, suffix_length, output_length;
	uint64_t span_length;

	result->data = NULL;
	result->length = 0;
	if (edit == NULL || edit->ref_span.begin > edit->ref_span.end ||
	    (reference.length != 0 && reference.data == NULL) ||
	    (edit->ref.length != 0 && edit->ref.data == NULL) ||
	    (edit->alt.length != 0 && edit->alt.data == NULL))
		return 0;
	begin = (size_t)edit->ref_span.begin;
	end = (size_t)edit->ref_span.end;
	if ((uint64_t)begin != edit->ref_span.begin ||
	    (uint64_t)end != edit->ref_span.end || end > reference.length)
		return 0;
	span_length = edit->ref_span.end - edit->ref_span.begin;
	if (span_length != (uint64_t)edit->ref.length)
		return 0;
	if (edit->ref.length != 0 &&
	    memcmp(reference.data + begin, edit->ref.data,
	    edit->ref.length) != 0)
		return 0;
	suffix_length = reference.length - end;
	if (edit->alt.length > SIZE_MAX - begin ||
	    suffix_length > SIZE_MAX - begin - edit->alt.length)
		return 0;
	output_length = begin + edit->alt.length + suffix_length;
	if (output_length == 0)
		return 1;
	result->data = (uint8_t *)malloc(output_length);
	if (result->data == NULL)
		return 0;
	if (begin != 0)
		memcpy(result->data, reference.data, begin);
	if (edit->alt.length != 0)
		memcpy(result->data + begin, edit->alt.data, edit->alt.length);
	if (suffix_length != 0)
		memcpy(result->data + begin + edit->alt.length,
		    reference.data + end, suffix_length);
	result->length = output_length;
	return 1;
}

static void
fuzz_failure(uint64_t case_state, size_t case_index,
	const fuzz_alleles_t *alleles, const char *expression, int line)
{
	fprintf(stderr,
	    "FUZZ_FAILED seed=0x%016" PRIx64 " case=%zu state=0x%016" PRIx64
	    " shape=%s ref_len=%zu alt_len=%zu line=%d: %s\n",
	    DUCKVEP_FUZZ_SEED, case_index, case_state,
	    fuzz_shape_name(alleles->shape), alleles->ref_length,
	    alleles->alt_length, line, expression);
	failures++;
}

#define FUZZ_REQUIRE(condition) do {                                      \
	if (!(condition)) {                                                \
		fuzz_failure(case_state, case_index, &alleles, #condition,   \
		    __LINE__);                                                \
		goto cleanup;                                                 \
	}                                                                 \
} while (0)

static int
run_fuzz_case(fuzz_rng_t *rng, uint64_t case_state, size_t case_index)
{
	fuzz_alleles_t alleles;
	expected_runs_t expected = {NULL, 0};
	owned_bytes_t raw_applied = {NULL, 0};
	owned_bytes_t semantic_applied = {NULL, 0};
	duckvep_bytes_view_t ref, alt, context_view;
	duckvep_span0_t raw_span, *actual_runs = NULL;
	duckvep_edit_view_t raw_edit;
	duckvep_semantic_edit_view_t semantic;
	duckvep_vep_diff_view_t actual;
	duckvep_allele_status_t status;
	uint8_t *context = NULL;
	size_t i, required, left, right, context_length;
	int passed = 0;

	memset(&alleles, 0, sizeof(alleles));
	FUZZ_REQUIRE(fuzz_alleles_build(rng, case_index, &alleles));
	ref = (duckvep_bytes_view_t){alleles.ref, alleles.ref_length};
	alt = (duckvep_bytes_view_t){alleles.alt, alleles.alt_length};
	FUZZ_REQUIRE(slow_expected_vep_runs(ref, alt, &expected));

	required = SIZE_MAX;
	status = duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR, ref,
	    alt, NULL, 0, &actual, &required);
	FUZZ_REQUIRE(status == (expected.count == 0 ? DUCKVEP_ALLELE_OK :
	    DUCKVEP_ALLELE_NO_SPACE));
	FUZZ_REQUIRE(required == expected.count);
	if (expected.count != 0) {
		FUZZ_REQUIRE(expected.count <= SIZE_MAX / sizeof(*actual_runs));
		actual_runs = (duckvep_span0_t *)malloc(
		    expected.count * sizeof(*actual_runs));
		FUZZ_REQUIRE(actual_runs != NULL);
	}
	status = duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR, ref,
	    alt, actual_runs, expected.count, &actual, &required);
	FUZZ_REQUIRE(status == DUCKVEP_ALLELE_OK);
	FUZZ_REQUIRE(required == expected.count);
	FUZZ_REQUIRE(actual.algorithm == DUCKVEP_DIFF_VEP_116_BYTE_XOR);
	FUZZ_REQUIRE(actual.reference.data == ref.data);
	FUZZ_REQUIRE(actual.alternate.data == alt.data);
	FUZZ_REQUIRE(actual.run_count == expected.count);
	FUZZ_REQUIRE(actual.runs == (expected.count == 0 ? NULL : actual_runs));
	for (i = 0; i < expected.count; i++) {
		FUZZ_REQUIRE(actual.runs[i].begin == expected.runs[i].begin);
		FUZZ_REQUIRE(actual.runs[i].end == expected.runs[i].end);
	}

	left = fuzz_bounded(rng, 65);
	right = fuzz_bounded(rng, 65);
	FUZZ_REQUIRE(alleles.ref_length <= SIZE_MAX - left);
	FUZZ_REQUIRE(right <= SIZE_MAX - left - alleles.ref_length);
	context_length = left + alleles.ref_length + right;
	if (context_length != 0) {
		context = (uint8_t *)malloc(context_length);
		FUZZ_REQUIRE(context != NULL);
	}
	for (i = 0; i < context_length; i++)
		context[i] = fuzz_byte(rng);
	if (alleles.ref_length != 0)
		memcpy(context + left, alleles.ref, alleles.ref_length);
	context_view = (duckvep_bytes_view_t){context, context_length};
	FUZZ_REQUIRE(duckvep_span0_from_length((uint64_t)left,
	    (uint64_t)alleles.ref_length, &raw_span) == DUCKVEP_ALLELE_OK);
	FUZZ_REQUIRE(duckvep_edit_view_make(raw_span, ref, alt, &raw_edit) ==
	    DUCKVEP_ALLELE_OK);
	FUZZ_REQUIRE(duckvep_semantic_edit_trim(&raw_edit, &semantic) ==
	    DUCKVEP_ALLELE_OK);
	FUZZ_REQUIRE(semantic.common_prefix <= ref.length);
	FUZZ_REQUIRE(semantic.common_prefix <= alt.length);
	FUZZ_REQUIRE(semantic.common_suffix <=
	    ref.length - semantic.common_prefix);
	FUZZ_REQUIRE(semantic.common_suffix <=
	    alt.length - semantic.common_prefix);
	FUZZ_REQUIRE(semantic.edit.ref.length == ref.length -
	    semantic.common_prefix - semantic.common_suffix);
	FUZZ_REQUIRE(semantic.edit.alt.length == alt.length -
	    semantic.common_prefix - semantic.common_suffix);
	FUZZ_REQUIRE(semantic.edit.ref_span.begin == raw_span.begin +
	    (uint64_t)semantic.common_prefix);
	FUZZ_REQUIRE(semantic.edit.ref_span.end == raw_span.end -
	    (uint64_t)semantic.common_suffix);
	FUZZ_REQUIRE(semantic.edit.ref.data == (ref.data == NULL ? NULL :
	    ref.data + semantic.common_prefix));
	FUZZ_REQUIRE(semantic.edit.alt.data == (alt.data == NULL ? NULL :
	    alt.data + semantic.common_prefix));
	FUZZ_REQUIRE(slow_apply_edit(context_view, &raw_edit, &raw_applied));
	FUZZ_REQUIRE(slow_apply_edit(context_view, &semantic.edit,
	    &semantic_applied));
	FUZZ_REQUIRE(raw_applied.length == semantic_applied.length);
	FUZZ_REQUIRE(raw_applied.length == 0 ||
	    memcmp(raw_applied.data, semantic_applied.data,
	    raw_applied.length) == 0);

	passed = 1;
cleanup:
	free(actual_runs);
	free(context);
	expected_runs_destroy(&expected);
	owned_bytes_destroy(&raw_applied);
	owned_bytes_destroy(&semantic_applied);
	fuzz_alleles_destroy(&alleles);
	return passed;
}

#undef FUZZ_REQUIRE

static void
test_deterministic_properties(void)
{
	fuzz_rng_t rng;
	uint64_t case_state;
	size_t case_index;

	rng.state = DUCKVEP_FUZZ_SEED;
	for (case_index = 0; case_index < DUCKVEP_FUZZ_CASES; case_index++) {
		case_state = rng.state;
		if (!run_fuzz_case(&rng, case_state, case_index))
			return;
	}
}

static void
test_checked_spans_and_views(void)
{
	duckvep_span0_t span;
	duckvep_bytes_view_t ref, alt;
	duckvep_edit_view_t edit;
	uint64_t length;

	CHECK(duckvep_span0_make(4, 9, &span) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_span0_length(&span, &length) == DUCKVEP_ALLELE_OK);
	CHECK(length == 5);
	CHECK(duckvep_span0_make(9, 4, &span) == DUCKVEP_ALLELE_INVALID);
	CHECK(duckvep_span0_from_length(UINT64_MAX, 1, &span) ==
	    DUCKVEP_ALLELE_RANGE);
	CHECK(duckvep_span0_from_length(UINT64_MAX, 0, &span) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(span.begin == UINT64_MAX && span.end == UINT64_MAX);

	CHECK(duckvep_bytes_view_make(NULL, 0, &ref) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_bytes_view_make(NULL, 1, &ref) ==
	    DUCKVEP_ALLELE_INVALID);
	CHECK(duckvep_bytes_view_make("AC", 2, &ref) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_bytes_view_make("GT", 2, &alt) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_span0_from_length(10, 2, &span) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_edit_view_make(span, ref, alt, &edit) ==
	    DUCKVEP_ALLELE_OK);
	span.end++;
	CHECK(duckvep_edit_view_make(span, ref, alt, &edit) ==
	    DUCKVEP_ALLELE_INVALID);
}

static void
test_semantic_trim(void)
{
	static const uint8_t insertion_ref[] = "A";
	static const uint8_t insertion_alt[] = "AT";
	static const uint8_t deletion_ref[] = "AT";
	static const uint8_t deletion_alt[] = "A";
	static const uint8_t delins_ref[] = "GACCT";
	static const uint8_t delins_alt[] = "GATCT";
	duckvep_span0_t span;
	duckvep_edit_view_t raw;
	duckvep_semantic_edit_view_t semantic;

	CHECK(duckvep_span0_from_length(100, 1, &span) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_edit_view_make(span,
	    (duckvep_bytes_view_t){insertion_ref, 1},
	    (duckvep_bytes_view_t){insertion_alt, 2}, &raw) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(duckvep_semantic_edit_trim(&raw, &semantic) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(semantic.common_prefix == 1 && semantic.common_suffix == 0);
	CHECK(semantic.edit.ref_span.begin == 101 &&
	    semantic.edit.ref_span.end == 101);
	CHECK(semantic.edit.ref.length == 0 && semantic.edit.alt.length == 1);
	CHECK(semantic.edit.ref.data == insertion_ref + 1);
	CHECK(semantic.edit.alt.data == insertion_alt + 1);
	CHECK(semantic.edit.alt.data[0] == 'T');

	CHECK(duckvep_span0_from_length(100, 2, &span) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_edit_view_make(span,
	    (duckvep_bytes_view_t){deletion_ref, 2},
	    (duckvep_bytes_view_t){deletion_alt, 1}, &raw) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(duckvep_semantic_edit_trim(&raw, &semantic) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(semantic.edit.ref_span.begin == 101 &&
	    semantic.edit.ref_span.end == 102);
	CHECK(semantic.edit.ref.length == 1 && semantic.edit.ref.data[0] == 'T');
	CHECK(semantic.edit.alt.length == 0);

	CHECK(duckvep_span0_from_length(7, 5, &span) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_edit_view_make(span,
	    (duckvep_bytes_view_t){delins_ref, 5},
	    (duckvep_bytes_view_t){delins_alt, 5}, &raw) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(duckvep_semantic_edit_trim(&raw, &semantic) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(semantic.common_prefix == 2 && semantic.common_suffix == 2);
	CHECK(semantic.edit.ref_span.begin == 9 &&
	    semantic.edit.ref_span.end == 10);
	CHECK(semantic.edit.ref.length == 1 && semantic.edit.ref.data[0] == 'C');
	CHECK(semantic.edit.alt.length == 1 && semantic.edit.alt.data[0] == 'T');

	CHECK(duckvep_edit_view_make(span, text_view("GACCT"),
	    text_view("GACCT"), &raw) == DUCKVEP_ALLELE_OK);
	CHECK(duckvep_semantic_edit_trim(&raw, &semantic) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(semantic.common_prefix == 5 && semantic.common_suffix == 0);
	CHECK(semantic.edit.ref_span.begin == 12 &&
	    semantic.edit.ref_span.end == 12);
	CHECK(semantic.edit.ref.length == 0 && semantic.edit.alt.length == 0);
}

static void
check_one_run(duckvep_bytes_view_t ref, duckvep_bytes_view_t alt,
	uint64_t begin, uint64_t end)
{
	duckvep_span0_t storage[2];
	duckvep_vep_diff_view_t view;
	size_t required;

	CHECK(duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR, ref, alt,
	    storage, 2, &view, &required) == DUCKVEP_ALLELE_OK);
	CHECK(required == 1 && view.run_count == 1);
	CHECK(view.algorithm == DUCKVEP_DIFF_VEP_116_BYTE_XOR);
	CHECK(view.reference.data == ref.data && view.alternate.data == alt.data);
	CHECK(view.runs == storage);
	CHECK(view.runs[0].begin == begin && view.runs[0].end == end);
}

static void
test_vep116_diff_runs(void)
{
	static const uint8_t nul_tail_ref[] = {'A', 'B'};
	static const uint8_t nul_tail_alt[] = {'A', 'B', 0, 'C'};
	duckvep_span0_t storage[2];
	duckvep_vep_diff_view_t view;
	size_t required;

	required = SIZE_MAX;
	CHECK(duckvep_diff_view_build(DUCKVEP_DIFF_ALGORITHM_NONE,
	    text_view("A"), text_view("T"), storage, 2, &view, &required) ==
	    DUCKVEP_ALLELE_INVALID);
	CHECK(required == 0 && view.algorithm == DUCKVEP_DIFF_ALGORITHM_NONE);

	CHECK(duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR,
	    text_view("ABCDEF"), text_view("AXCYEF"), storage, 2, &view,
	    &required) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(required == 2 && view.run_count == 2);
	CHECK(view.runs[0].begin == 1 && view.runs[0].end == 2);
	CHECK(view.runs[1].begin == 3 && view.runs[1].end == 4);

	/* Perl's unequal string XOR zero-pads to the longer operand. */
	check_one_run(text_view("AB"), text_view("ABCD"), 2, 4);
	check_one_run(text_view("ABCD"), text_view("AB"), 2, 4);
	check_one_run((duckvep_bytes_view_t){nul_tail_ref, 2},
	    (duckvep_bytes_view_t){nul_tail_alt, 4}, 3, 4);

	/* VEP's inclusive e=-1 insertion sentinel is the empty span [0,0). */
	check_one_run(text_view(""), text_view("T"), 0, 0);
	check_one_run(text_view(""), text_view("TT"), 0, 0);
	check_one_run(text_view("A"), text_view(""), 0, 1);
	check_one_run(text_view("ABC"), text_view(""), 0, 3);
	check_one_run(text_view("ABC"), text_view("T"), 0, 3);
	check_one_run(text_view("A"), text_view("T"), 0, 1);

	CHECK(duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR,
	    text_view("ABCD"), text_view("ABCD"), storage, 2, &view,
	    &required) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(required == 0 && view.run_count == 0 && view.runs == NULL);

	storage[0].begin = 99;
	storage[0].end = 100;
	CHECK(duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR,
	    text_view("ABCDEF"), text_view("AXCYEF"), storage, 1, &view,
	    &required) ==
	    DUCKVEP_ALLELE_NO_SPACE);
	CHECK(required == 2 && view.run_count == 0 && view.runs == NULL);
	CHECK(storage[0].begin == 99 && storage[0].end == 100);
}

static void
test_long_alleles(void)
{
	const size_t length = 70003;
	uint8_t *ref, *alt;
	duckvep_span0_t span, storage[2];
	duckvep_edit_view_t raw;
	duckvep_semantic_edit_view_t semantic;
	duckvep_vep_diff_view_t diff;
	size_t required;

	ref = (uint8_t *)malloc(length);
	alt = (uint8_t *)malloc(length);
	CHECK(ref != NULL && alt != NULL);
	if (ref == NULL || alt == NULL) {
		free(ref);
		free(alt);
		return;
	}
	memset(ref, 'A', length);
	memset(alt, 'A', length);
	ref[65536] = 'C';
	alt[65536] = 'G';
	ref[70001] = 'T';
	alt[70001] = 'C';

	CHECK(duckvep_span0_from_length(17, length, &span) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(duckvep_edit_view_make(span,
	    (duckvep_bytes_view_t){ref, length},
	    (duckvep_bytes_view_t){alt, length}, &raw) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(duckvep_semantic_edit_trim(&raw, &semantic) ==
	    DUCKVEP_ALLELE_OK);
	CHECK(semantic.common_prefix == 65536);
	CHECK(semantic.common_suffix == 1);
	CHECK(semantic.edit.ref.data == ref + 65536);
	CHECK(semantic.edit.ref.length == length - 65537);
	CHECK(semantic.edit.ref_span.begin == 65553);
	CHECK(semantic.edit.ref_span.end == 70019);

	CHECK(duckvep_diff_view_build(DUCKVEP_DIFF_VEP_116_BYTE_XOR, raw.ref,
	    raw.alt, storage, 2, &diff, &required) == DUCKVEP_ALLELE_OK);
	CHECK(required == 2 && diff.run_count == 2);
	CHECK(diff.runs[0].begin == 65536 && diff.runs[0].end == 65537);
	CHECK(diff.runs[1].begin == 70001 && diff.runs[1].end == 70002);

	free(ref);
	free(alt);
}

int
main(void)
{
	test_checked_spans_and_views();
	test_semantic_trim();
	test_vep116_diff_runs();
	test_long_alleles();
	test_deterministic_properties();

	if (failures != 0) {
		fprintf(stderr, "FAILURES=%d\n", failures);
		return 1;
	}
	printf("DUCKVEP_ALLELE_OK=1\n");
	return 0;
}
