/* Test-only R .C bridge. The production reducer owns raw-GT interpretation;
 * this bridge only converts its fixed result into R-owned integer arrays. */
#include "duckvep_phase.h"
#include <string.h>

void duckhts_test_vep116_raw_gt(char **gt, int *alt_count, int *count,
    int *status, int *disposition, int *ploidy, int *missing, int *slots,
    int *first, int *second) {
    for (int i = 0; i < *count; i++) {
        duckvep_raw_gt_t parsed;
        status[i] = duckvep_phase_parse_vep116_raw((const uint8_t *)gt[i], strlen(gt[i]),
            (uint32_t)alt_count[i], &parsed);
        disposition[i] = parsed.disposition;
        ploidy[i] = parsed.source_ploidy;
        missing[i] = parsed.source_has_missing;
        slots[i] = (int)parsed.parsed_slots;
        first[i] = parsed.allele_index[0] == UINT32_MAX ? -1 : (int)parsed.allele_index[0];
        second[i] = parsed.allele_index[1] == UINT32_MAX ? -1 : (int)parsed.allele_index[1];
    }
}
