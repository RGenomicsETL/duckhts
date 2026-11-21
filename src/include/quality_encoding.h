#ifndef DUCKHTS_QUALITY_ENCODING_H
#define DUCKHTS_QUALITY_ENCODING_H

#include <stddef.h>
#include <stdint.h>

#include "hts_file_access.h"

typedef enum {
    DUCKHTS_QUALITY_ENCODING_AUTO = 0,
    DUCKHTS_QUALITY_ENCODING_PHRED33,
    DUCKHTS_QUALITY_ENCODING_PHRED64,
    DUCKHTS_QUALITY_ENCODING_SOLEXA64
} duckhts_quality_encoding;

typedef enum {
    DUCKHTS_QUALITY_REPR_STRING = 0,
    DUCKHTS_QUALITY_REPR_PHRED
} duckhts_quality_representation;

typedef struct {
    int observed_ascii_min;
    int observed_ascii_max;
    int records_sampled;
    int compatible_phred33;
    int compatible_phred64;
    int compatible_solexa64;
    int is_ambiguous;
    duckhts_quality_encoding guessed_encoding;
} duckhts_quality_detect_result;

const char *duckhts_quality_encoding_name(duckhts_quality_encoding enc);
const char *duckhts_quality_representation_name(duckhts_quality_representation repr);
int duckhts_parse_quality_encoding(const char *text, int allow_auto,
                                   duckhts_quality_encoding *out);
int duckhts_parse_quality_representation(const char *text,
                                         duckhts_quality_representation *out);
uint8_t duckhts_normalize_fastq_quality(uint8_t stored_qual,
                                        duckhts_quality_encoding enc);
/* Reads path through access (hts_file_access.h). */
int duckhts_detect_fastq_quality_encoding(const duckhts_file_access_t *access,
                                          const char *path, int64_t max_records,
                                          duckhts_quality_detect_result *out,
                                          char *errbuf, size_t errbuf_len);

#endif
