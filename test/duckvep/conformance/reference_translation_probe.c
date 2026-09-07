/* Test-only bridge for the complete raw translation, before reference rules.
 * No imitation of Ensembl's reference rules belongs in this probe. */
#include "duckvep_codon.h"
#include <string.h>

void duckhts_test_raw_translation(char **cds, int *table, unsigned char *peptide,
    int *capacity, int *status, double *facts) {
    duckvep_translation_t result = {0};
    *status = DUCKVEP_TRANSLATION_INVALID_ARG;
    memset(facts, 0, 3u * sizeof(*facts));
    if (*capacity < 1) return;
    *status = duckvep_translate_cds((const uint8_t *)*cds, strlen(*cds),
        (duckvep_codon_table_t)*table, DUCKVEP_TRANSLATION_N_CONSENSUS, peptide, (size_t)*capacity, &result);
    facts[0] = (double)result.length;
    facts[1] = (double)result.first_stop_position1;
    facts[2] = (double)result.unambiguous;
}
