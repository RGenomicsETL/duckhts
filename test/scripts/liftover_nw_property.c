#include "liftover_nw_limit.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static uint64_t state;
static uint64_t run_seed;

static uint64_t next_u64(void) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static int check_case(size_t s_len, size_t t_len, uint64_t trial) {
    size_t cells = 0, reverse_cells = 0;
    int got = liftover_nw_cell_count(s_len, t_len, &cells);
    int reverse = liftover_nw_cell_count(t_len, s_len, &reverse_cells);
    int expected = s_len < LIFTOVER_NW_MAX_CELLS && t_len < LIFTOVER_NW_MAX_CELLS &&
                   s_len + 1 <= LIFTOVER_NW_MAX_CELLS / (t_len + 1);
    if (got != expected || reverse != expected || (got && (cells != reverse_cells || cells > LIFTOVER_NW_MAX_CELLS))) {
        fprintf(stderr, "seed=%" PRIu64 " trial=%" PRIu64 " s=%zu t=%zu expected=%d got=%d cells=%zu\n",
                run_seed, trial, s_len, t_len, expected, got, cells);
        return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    uint64_t seed = argc > 1 ? strtoull(argv[1], NULL, 10) : 169;
    uint64_t trials = argc > 2 ? strtoull(argv[2], NULL, 10) : 100000;
    static const size_t edges[] = {0, 1, 2, 2047, 2048, 4095, 4096,
                                   LIFTOVER_NW_MAX_CELLS - 1, LIFTOVER_NW_MAX_CELLS, SIZE_MAX};
    run_seed = seed ? seed : 169;
    state = run_seed;
    for (size_t i = 0; i < sizeof(edges) / sizeof(edges[0]); i++)
        for (size_t j = 0; j < sizeof(edges) / sizeof(edges[0]); j++)
            if (!check_case(edges[i], edges[j], i * 100 + j)) return 1;
    for (uint64_t i = 0; i < trials; i++) {
        size_t s_len = (i & 3) ? (size_t)(next_u64() % (LIFTOVER_NW_MAX_CELLS + 4096)) : (size_t)next_u64();
        size_t t_len = (i & 3) ? (size_t)(next_u64() % (LIFTOVER_NW_MAX_CELLS + 4096)) : (size_t)next_u64();
        if (!check_case(s_len, t_len, i)) return 1;
    }
    printf("liftover NW properties: OK (%" PRIu64 " trials, seed=%" PRIu64 ")\n", trials, seed);
    return 0;
}
