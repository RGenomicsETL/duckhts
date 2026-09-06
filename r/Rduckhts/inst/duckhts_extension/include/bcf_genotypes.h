#ifndef DUCKHTS_BCF_GENOTYPES_H
#define DUCKHTS_BCF_GENOTYPES_H

#include "bcf_scan.h"

/* Worker-owned HTSlib decode buffers. Reused across records; destroy with the
 * matching function. Borrowed sample slices expire at the next decode. These
 * are adaptive HTSlib buffers, not an allocation-static execution guarantee. */
typedef struct {
    int32_t *gt;
    int32_t *ps;
    int gt_capacity;
    int ps_capacity;
    int gt_stride; /* zero means absent/invalid GT, not known zero ploidy */
    int ps_present;
    int samples;
} duckhts_bcf_genotypes_t;

void duckhts_bcf_genotypes_destroy(duckhts_bcf_genotypes_t *values);
/* Decode GT once, then scalar PS once. Missing fields are not errors. A bad
 * field becomes absent under null/warn, without discarding a valid other
 * field; allocation and physical decoding failures always return 0. */
int duckhts_bcf_genotypes_decode(duckhts_bcf_genotypes_t *values, bcf_hdr_t *header,
                                 bcf1_t *record, duckhts_bcf_decode_policy_t policy,
                                 char *error, size_t error_size);

/* GT has already been validated. Stops at HTSlib's vector-end sentinel;
 * missing alleles retain slots, and phase flags are the decoded HTSlib flags,
 * including its leading-slot convention for VCF versions before 4.4. */
int duckhts_bcf_genotype_ploidy(const int32_t *gt, int stride);
int duckhts_bcf_genotype_has_alt(const int32_t *gt, int ploidy);

#endif
