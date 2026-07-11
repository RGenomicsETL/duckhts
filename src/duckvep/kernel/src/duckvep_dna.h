/* Small ASCII DNA operations shared by the consequence kernels. */
#ifndef DUCKVEP_DNA_H
#define DUCKVEP_DNA_H

#include <stdint.h>

#define DUCKVEP_DNA_ASCII_MASK \
    ((UINT32_C(1) << ('A' & 31)) | (UINT32_C(1) << ('C' & 31)) | \
     (UINT32_C(1) << ('G' & 31)) | (UINT32_C(1) << ('T' & 31)))

static inline char duckvep_dna_normalize(char base, int allow_n) {
    unsigned char upper = (unsigned char)base & UINT8_C(0xdf);
    uint32_t bit;

    if (upper == 'U') upper = 'T';
    if (upper == 'N') return allow_n ? 'N' : '\0';
    if (upper < 'A' || upper > 'Z') return '\0';
    bit = UINT32_C(1) << (upper & 31);
    return (DUCKVEP_DNA_ASCII_MASK & bit) != 0u ? (char)upper : '\0';
}

static inline char duckvep_dna_complement(char base) {
    /* The caller supplies normalized A/C/G/T: A^21=T and C^4=G. */
    return (char)((unsigned char)base ^
        (((unsigned char)base & 2u) != 0u ? 4u : 21u));
}

static inline int duckvep_dna_codon_code(char base) {
    /* The translation tables use T=0, C=1, A=2, G=3. The low three ASCII
     * bits distinguish the four normalized bases without a branch ladder. */
    static const uint8_t code[8] = {0u, 2u, 0u, 1u, 0u, 0u, 0u, 3u};
    char normalized = duckvep_dna_normalize(base, 0);

    return normalized == '\0' ? -1 : (int)code[(unsigned char)normalized & 7u];
}

#endif /* DUCKVEP_DNA_H */
