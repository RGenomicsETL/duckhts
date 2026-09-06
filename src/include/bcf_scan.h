#ifndef DUCKHTS_BCF_SCAN_H
#define DUCKHTS_BCF_SCAN_H

#include "bcf_index_snapshot.h"
#include "hts_io_tuning.h"

/* Zero-initialize before open. One worker owns the file, mutable header,
 * iterator and text buffer. The parsed index is borrowed and must outlive the
 * scan. Records are caller-owned, so decoded views remain valid until the
 * caller explicitly advances that record; output chunk edges do not advance it. */
typedef struct {
    htsFile *fp;
    bcf_hdr_t *hdr;
    const duckhts_bcf_index_t *index;
    hts_itr_t *itr;
    kstring_t line;
    int indexed; // An empty indexed selection must never fall back to streaming.
} duckhts_bcf_scan_t;

int duckhts_bcf_scan_open(duckhts_bcf_scan_t *scan, const char *path,
                         const duckhts_bcf_index_t *index, int decompression_threads,
                         duckhts_hts_io_profile_t profile, const char *reader_name,
                         char *error, size_t error_size);
void duckhts_bcf_scan_close(duckhts_bcf_scan_t *scan);

/* Select one native HTSlib region union. Returns 1 on success (including an
 * empty selection), 0 on invalid input or iterator allocation failure.
 * Replaces only this worker's iterator. */
int duckhts_bcf_scan_regions(duckhts_bcf_scan_t *scan, char **regions,
                            unsigned int count, char *error, size_t error_size);

/* Select one literal contig from the scan dictionary, not a region expression. */
int duckhts_bcf_scan_contig(duckhts_bcf_scan_t *scan, const char *name,
                           char *error, size_t error_size);

/* HTSlib status: >=0 record, -1 EOF, <-1 read/parse failure. No unpacking or
 * output materialization beyond that required by the format reader. */
int duckhts_bcf_scan_next(duckhts_bcf_scan_t *scan, bcf1_t *record);

typedef enum {
    DUCKHTS_BCF_FIELD_INFO,
    DUCKHTS_BCF_FIELD_FORMAT,
    DUCKHTS_BCF_FIELD_GT
} duckhts_bcf_field_class_t;

/* Pure validation: no warning callbacks or host state. GT has String header
 * metadata but integer payload; other FORMAT types follow the canonical
 * header-driven decoder. Absent fields are not type failures. */
int duckhts_bcf_check_field_type(bcf_hdr_t *hdr, bcf1_t *record,
                                duckhts_bcf_field_class_t field_class,
                                int header_id, int header_type,
                                const char *reader_name, char *error, size_t error_size);
int duckhts_bcf_check_format_width(const char *reader_name, const char *tag,
                                  bcf_hdr_t *hdr, bcf1_t *record, int values,
                                  int samples, char *error, size_t error_size);

typedef enum {
    DUCKHTS_BCF_DECODE_OK,
    DUCKHTS_BCF_DECODE_TYPE_MISMATCH,
    DUCKHTS_BCF_DECODE_FATAL
} duckhts_bcf_decode_status_t;

/* Classify bcf_get_* return codes. The caller chooses null/warn/error for a
 * type mismatch; undefined header fields and OOM always remain fatal. */
duckhts_bcf_decode_status_t duckhts_bcf_decode_status(
    const char *reader_name, const char *field_class, const char *tag,
    bcf_hdr_t *hdr, bcf1_t *record, int ret, char *error, size_t error_size);

#endif
