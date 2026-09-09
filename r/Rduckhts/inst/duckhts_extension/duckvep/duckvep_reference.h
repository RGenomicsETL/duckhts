/* Worker-owned reference transport for borrowed HGVS genomic windows. */
#ifndef DUCKVEP_REFERENCE_H
#define DUCKVEP_REFERENCE_H

#include "kernel/src/duckvep_hgvs.h"
#include <htslib/faidx.h>

struct duckvep_owned_model;

/* Two maximum uint16 literal alleles, VEP's two shift flanks, read-ahead and
 * FASTA line-ending scratch fit the ordinary independent-event workspace. */
#define DUCKVEP_REFERENCE_DEFAULT_BYTES 262144u

typedef struct duckvep_reference_reader {
    const struct duckvep_owned_model *model;
    faidx_t *fai;
    char *bases;
    size_t capacity, length;
    uint32_t start1;
    uint16_t chrom_id;
} duckvep_reference_reader_t;

/* Initialization owns one mutable faidx handle. The pinned model and supplied
 * byte span are borrowed and must outlive the reader. Teardown destroys fai;
 * the caller frees its byte span. A model without FASTA opens no handle. */
int duckvep_reference_reader_init(duckvep_reference_reader_t *reader,
    const struct duckvep_owned_model *model, char *bases, size_t capacity,
    char *error, size_t error_size);

/* No first-party allocation, resizing or freeing. A miss overwrites the fixed
 * byte span. Views expire at the next miss. The exact VEP shift slice is kept
 * separate from wider lookup/read-ahead storage. HTSlib owns transport buffers.
 * Failure clears both views and available; capacity pressure is an error. */
int duckvep_reference_reader_windows(duckvep_reference_reader_t *reader,
    const duckvep_event_t *event, int *available,
    duckvep_hgvs_reference_window_t *shift, duckvep_hgvs_reference_window_t *lookup,
    char *error, size_t error_size);

#endif
