/*
 * Private DuckVEP allele and edit primitives.
 *
 * These are builder-time components, not a public extension ABI.  Every byte
 * view is borrowed.  A successful builder publishes only const pointers; the
 * caller owns the input bytes and any run storage for the lifetime of the
 * resulting view.
 */
#ifndef DUCKVEP_ALLELE_H
#define DUCKVEP_ALLELE_H

#include <stddef.h>
#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
#define DUCKVEP_INTERNAL __attribute__((__visibility__("hidden")))
#else
#define DUCKVEP_INTERNAL
#endif

typedef enum duckvep_allele_status {
	DUCKVEP_ALLELE_OK = 0,
	DUCKVEP_ALLELE_INVALID,
	DUCKVEP_ALLELE_RANGE,
	DUCKVEP_ALLELE_NO_SPACE
} duckvep_allele_status_t;

typedef struct duckvep_span0 {
	uint64_t begin;
	uint64_t end;
} duckvep_span0_t;

typedef struct duckvep_bytes_view {
	const uint8_t *data;
	size_t length;
} duckvep_bytes_view_t;

/* One replacement on the forward-reference axis. */
typedef struct duckvep_edit_view {
	duckvep_span0_t ref_span;
	duckvep_bytes_view_t ref;
	duckvep_bytes_view_t alt;
} duckvep_edit_view_t;

/* A borrowed slice of the raw edit after maximal prefix-then-suffix trim. */
typedef struct duckvep_semantic_edit_view {
	duckvep_edit_view_t edit;
	size_t common_prefix;
	size_t common_suffix;
} duckvep_semantic_edit_view_t;

typedef enum duckvep_diff_algorithm {
	DUCKVEP_DIFF_ALGORITHM_NONE = 0,
	DUCKVEP_DIFF_VEP_116_BYTE_XOR = 1
} duckvep_diff_algorithm_t;

/*
 * Runs use offsets into the caller-supplied VEP minimized overlap-allele pair,
 * not genomic coordinates.  The caller decodes VEP's "-" allele sentinel to
 * an empty byte view.  This builder does not trim or substitute the separate
 * semantic edit envelope.  Unequal-length XOR means a run can extend past the
 * shorter allele.  Empty insertion is represented by the one empty run [0,0).
 */
typedef struct duckvep_vep_diff_view {
	duckvep_diff_algorithm_t algorithm;
	duckvep_bytes_view_t reference;
	duckvep_bytes_view_t alternate;
	const duckvep_span0_t *runs;
	size_t run_count;
} duckvep_vep_diff_view_t;

DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_span0_make(uint64_t, uint64_t,
	duckvep_span0_t *);
DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_span0_from_length(uint64_t, uint64_t,
	duckvep_span0_t *);
DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_span0_length(const duckvep_span0_t *,
	uint64_t *);

DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_bytes_view_make(const void *, size_t,
	duckvep_bytes_view_t *);
DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_edit_view_make(duckvep_span0_t,
	duckvep_bytes_view_t, duckvep_bytes_view_t, duckvep_edit_view_t *);
DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_semantic_edit_trim(
	const duckvep_edit_view_t *, duckvep_semantic_edit_view_t *);

/*
 * Build a const borrowed view over caller-owned run_storage for the named
 * compatibility algorithm.  required_runs is set on OK and NO_SPACE;
 * NO_SPACE never publishes or partially fills a view.
 */
DUCKVEP_INTERNAL duckvep_allele_status_t duckvep_diff_view_build(
	duckvep_diff_algorithm_t, duckvep_bytes_view_t, duckvep_bytes_view_t,
	duckvep_span0_t *, size_t, duckvep_vep_diff_view_t *, size_t *);

#endif /* DUCKVEP_ALLELE_H */
