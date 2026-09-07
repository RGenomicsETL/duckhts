/* Test-only R bridge. Production callers own/reuse this storage at query init. */
#include "duckvep_sequence_diff.h"
#include <stdlib.h>
#include <string.h>

void duckhts_test_sequence_differences(char **reference, char **alternate, int *align,
    double *spans, int *capacity, int *status, double *facts) {
    size_t n = strlen(*reference), m = strlen(*alternate);
    *status = DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    memset(facts, 0, 3u * sizeof(*facts));
    if (n > 10000u || m > 10000u || *capacity < 0 || *capacity > 20000) return;
    duckvep_sequence_diff_scratch_t scratch = {0};
    scratch.score_capacity = (m + 1u) * 2u;
    scratch.trace_capacity = (n + 1u) * (m + 1u);
    scratch.scores = malloc(scratch.score_capacity * sizeof(*scratch.scores));
    scratch.trace = malloc(scratch.trace_capacity);
    duckvep_sequence_difference_t *out = calloc((size_t)*capacity + 1u, sizeof(*out));
    if (!scratch.scores || !scratch.trace || !out) goto cleanup;
    duckvep_sequence_diff_result_t result;
    *status = duckvep_sequence_differences((const uint8_t *)*reference, n,
        (const uint8_t *)*alternate, m, *align, &scratch, out, (size_t)*capacity, &result);
    facts[0] = (double)result.count; facts[1] = (double)result.alignment_length;
    facts[2] = (double)result.trace_cells;
    if (*status == DUCKVEP_SEQUENCE_DIFF_OK) for (size_t i = 0u; i < result.count; i++) {
        spans[i * 5u] = (double)out[i].ref_start0;
        spans[i * 5u + 1u] = (double)out[i].alt_start0;
        spans[i * 5u + 2u] = (double)out[i].ref_length;
        spans[i * 5u + 3u] = (double)out[i].alt_length;
        spans[i * 5u + 4u] = (double)out[i].alignment_start0;
    }
cleanup:
    free(out); free(scratch.trace); free(scratch.scores);
}
