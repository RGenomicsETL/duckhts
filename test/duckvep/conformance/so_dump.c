/*
 * so_dump.c — print the C engine's SO vocabulary via the real kernel ABI
 * (name / impact / rank / tier). so_conformance.R diffs this against the
 * VEP-derived class model (data/so_consequences.tsv) so the engine's metadata is
 * checked against Ensembl VEP's authoritative %OVERLAP_CONSEQUENCES spec, not
 * hand-asserted.
 *
 * This file is built by so_conformance.R against the pure C kernel.
 */
#include "duckvep_so.h"

#include <stdio.h>

static const char *impact_name(duckvep_impact_t i) {
    switch (i) {
        case DUCKVEP_IMPACT_HIGH:     return "HIGH";
        case DUCKVEP_IMPACT_MODERATE: return "MODERATE";
        case DUCKVEP_IMPACT_LOW:      return "LOW";
        default:                      return "MODIFIER";
    }
}

int main(void) {
    int b;
    printf("so_term\timpact\trank\ttier\n");
    for (b = 0; b < DUCKVEP_SO_BIT_COUNT; b++) {
        printf("%s\t%s\t%u\t%u\n", duckvep_so_name((duckvep_so_bit_t)b),
               impact_name(duckvep_so_bit_impact((duckvep_so_bit_t)b)),
               (unsigned)duckvep_so_rank((duckvep_so_bit_t)b),
               (unsigned)duckvep_so_tier((duckvep_so_bit_t)b));
    }
    return 0;
}
