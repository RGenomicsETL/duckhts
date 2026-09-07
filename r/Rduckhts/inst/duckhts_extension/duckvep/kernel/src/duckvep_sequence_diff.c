/* Alignment score/tie and difference-run semantics from Ensembl Variation
 * release/116 (Apache-2.0), credited in duckvep_sequence_diff.h. This bounded
 * implementation stores two score rows and an exact band of traceback cells. */
#include "duckvep_sequence_diff.h"

#include <string.h>

enum { DIAGONAL, DELETE_REFERENCE, INSERT_ALTERNATE, MATCH };

static size_t minimum(size_t a, size_t b) { return a < b ? a : b; }
static size_t maximum(size_t a, size_t b) { return a > b ? a : b; }

/* This is a feasible alignment, not a trimmed problem: retain the common
 * prefix/suffix and align the middle positionally. NW still sees every byte,
 * so repeat-associated tie placement is unaffected. */
static uint64_t cost_bound(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m) {
    size_t first = 0u, last = 0u, shared = minimum(n, m);
    while (first < shared && ref[first] == alt[first]) first++;
    while (last < shared - first && ref[n - 1u - last] == alt[m - 1u - last]) last++;
    uint64_t cost = (uint64_t)(maximum(n, m) - shared) * 3u;
    for (size_t i = first; i < shared - last; i++) if (ref[i] != alt[i]) cost += 4u;
    return cost;
}

static void build_trace(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
                        size_t band, size_t stride, const duckvep_sequence_diff_scratch_t *b) {
    const uint64_t infinity = UINT64_MAX / 2u;
    uint64_t *previous = b->scores, *current = b->scores + m + 1u;
    for (size_t j = 0u; j <= m; j++) previous[j] = infinity;
    for (size_t j = 0u; j <= minimum(m, band); j++) {
        previous[j] = (uint64_t)j * 3u; b->trace[j] = INSERT_ALTERNATE;
    }
    for (size_t i = 1u; i <= n; i++) {
        size_t start = i > band ? i - band : 0u;
        size_t end = minimum(m, i + band);
        if (start) current[start - 1u] = infinity;
        if (end > minimum(m, i - 1u + band)) previous[end] = infinity;
        for (size_t j = start; j <= end; j++) {
            uint8_t direction = DELETE_REFERENCE;
            uint64_t score = (uint64_t)i * 3u;
            if (j) {
                uint64_t diagonal = previous[j - 1u] + (ref[i - 1u] == alt[j - 1u] ? 0u : 4u);
                uint64_t insertion = current[j - 1u] + 3u, deletion = previous[j] + 3u;
                /* Minimizing 4*substitutions + 3*gaps is exactly maximizing
                 * matches - substitutions - gaps at fixed endpoint lengths. */
                if (diagonal < insertion && diagonal < deletion) {
                    score = diagonal; direction = DIAGONAL;
                } else if (insertion < deletion) {
                    score = insertion; direction = INSERT_ALTERNATE;
                } else score = deletion;
            }
            current[j] = score;
            b->trace[i * stride + j - start] = direction;
        }
        uint64_t *swap = previous; previous = current; current = swap;
    }
}

/* Walk once to count complete runs and once to fill. Reverse traversal groups
 * exactly the same adjacent gap classes as upstream's forward difference scan. */
static void traceback(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
    const uint8_t *trace, size_t band, size_t stride,
    duckvep_sequence_difference_t *out, duckvep_sequence_diff_result_t *result) {
    size_t i = n, j = m, columns = 0u, count = 0u;
    uint8_t previous = MATCH;
    while (i || j) {
        uint8_t direction;
        if (trace) {
            size_t start = i > band ? i - band : 0u;
            direction = trace[i * stride + j - start];
        } else {
            direction = i == j ? DIAGONAL : i > j ? DELETE_REFERENCE : INSERT_ALTERNATE;
        }
        size_t rn = direction != INSERT_ALTERNATE, an = direction != DELETE_REFERENCE;
        i -= rn; j -= an;
        uint8_t kind = direction == DIAGONAL && ref[i] == alt[j] ? MATCH : direction;
        if (kind != MATCH) {
            if (kind != previous) {
                count++;
                if (out) out[count - 1u] = (duckvep_sequence_difference_t){i, j, 0u, 0u, columns};
            }
            if (out) {
                duckvep_sequence_difference_t *d = &out[count - 1u];
                d->ref_start0 = i; d->alt_start0 = j;
                d->ref_length += rn; d->alt_length += an;
            }
        }
        previous = kind; columns++;
    }
    if (out) {
        for (size_t k = 0u; k < count; k++)
            out[k].alignment_start0 = columns - out[k].alignment_start0 -
                maximum(out[k].ref_length, out[k].alt_length);
        for (size_t k = 0u; k < count / 2u; k++) {
            duckvep_sequence_difference_t swap = out[k];
            out[k] = out[count - 1u - k]; out[count - 1u - k] = swap;
        }
    }
    result->count = count; result->alignment_length = columns;
}

duckvep_sequence_diff_status_t duckvep_sequence_differences(
    const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
    int align_indels, const duckvep_sequence_diff_scratch_t *b,
    duckvep_sequence_difference_t *out, size_t capacity, duckvep_sequence_diff_result_t *result) {
    if (!result) return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    memset(result, 0, sizeof(*result));
    if ((!ref && n) || (!alt && m) || (!out && capacity) ||
        (align_indels != 0 && align_indels != 1) ||
        n > SIZE_MAX / 8u || m > SIZE_MAX / 8u)
        return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    for (size_t i = 0u; i < n; i++)
        if (!ref[i] || ref[i] >= 128u || ref[i] == '-') return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    for (size_t j = 0u; j < m; j++)
        if (!alt[j] || alt[j] >= 128u || alt[j] == '-') return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    size_t band = 0u, stride = 0u;
    const uint8_t *trace = NULL;
    if (align_indels && n && m && (n != m || memcmp(ref, alt, n))) {
        /* Every optimum costs at most this feasible path. Leaving |i-j| <= U/3
         * needs more than U in gap cost alone, so no optimal traceback is lost. */
        band = (size_t)minimum(maximum(n, m), (size_t)(cost_bound(ref, n, alt, m) / 3u));
        stride = band >= m / 2u ? m + 1u : band * 2u + 1u;
        if (n + 1u > SIZE_MAX / stride) return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
        result->trace_cells = (n + 1u) * stride;
        if (!b || !b->trace || result->trace_cells > b->trace_capacity)
            return DUCKVEP_SEQUENCE_DIFF_TRACE_FULL;
        if (!b->scores || m + 1u > b->score_capacity / 2u)
            return DUCKVEP_SEQUENCE_DIFF_SCORE_FULL;
        build_trace(ref, n, alt, m, band, stride, b);
        trace = b->trace;
    }
    traceback(ref, n, alt, m, trace, band, stride, NULL, result);
    if (result->count > capacity) return DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL;
    traceback(ref, n, alt, m, trace, band, stride, out, result);
    return DUCKVEP_SEQUENCE_DIFF_OK;
}
