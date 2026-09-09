/* Shared, host-neutral FORMAT decoding over HTSlib. */
#include "include/bcf_format.h"
#include <htslib/hts_log.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void duckhts_bcf_format_destroy(duckhts_bcf_format_t *values) {
    free(values->data);
    if (values->strings) free(values->strings[0]);
    free(values->strings);
    memset(values, 0, sizeof(*values));
}

static int accept_mismatch(duckhts_bcf_decode_policy_t policy, const char *error) {
    if (policy == DUCKHTS_BCF_DECODE_ERROR) return 0;
    if (policy == DUCKHTS_BCF_DECODE_WARN) hts_log_warning("%s", error);
    return 1;
}

int duckhts_bcf_format_decode(duckhts_bcf_format_t *values, bcf_hdr_t *header,
                              bcf1_t *record, const char *tag, int header_type,
                              duckhts_bcf_decode_policy_t policy, const char *reader_name,
                              char *error, size_t error_size) {
    if (values->loaded) return 1;
    values->count = values->stride = 0;
    int samples = bcf_hdr_nsamples(header);
    if (!samples) goto absent;
    if (bcf_unpack(record, BCF_UN_FMT) < 0 || record->n_sample != (unsigned)samples) {
        snprintf(error, error_size, "%s: failed to unpack FORMAT or sample count differs from header",
                 reader_name);
        return 0;
    }
    int id = bcf_hdr_id2int(header, BCF_DT_ID, tag);
    bcf_fmt_t *format = id >= 0 ? bcf_get_fmt_id(record, id) : NULL;
    if (!format || !format->p) goto absent;
    int is_gt = strcmp(tag, "GT") == 0 && header_type == BCF_HT_STR;
    if (is_gt && (bcf_hdr_id2length(header, BCF_HL_FMT, id) != BCF_VL_FIXED ||
                  bcf_hdr_id2number(header, BCF_HL_FMT, id) != 1)) {
        snprintf(error, error_size, "%s: FORMAT/GT must have Number=1,Type=String", reader_name);
        goto mismatch;
    }
    if (!duckhts_bcf_check_field_type(header, record,
            is_gt ? DUCKHTS_BCF_FIELD_GT : DUCKHTS_BCF_FIELD_FORMAT,
            id, header_type, reader_name, error, error_size)) goto mismatch;

    int numeric = is_gt || header_type == BCF_HT_INT || header_type == BCF_HT_REAL;
    size_t width = (size_t)format->n + (numeric ? 0 : 1);
    size_t element_bytes = numeric ? sizeof(int32_t) : 1;
    if (format->n < 0 || width > (size_t)INT_MAX / (size_t)samples ||
        width > SIZE_MAX / element_bytes / (size_t)samples ||
        (size_t)samples > SIZE_MAX / sizeof(char *)) {
        snprintf(error, error_size, "%s: FORMAT/%s exceeds the supported decoded-value capacity",
                 reader_name, tag);
        return 0;
    }
    int ret = numeric
        ? bcf_get_format_values(header, record, tag, &values->data, &values->capacity,
                                is_gt ? BCF_HT_INT : header_type)
        : bcf_get_format_string(header, record, tag, &values->strings, &values->capacity);
    duckhts_bcf_decode_status_t status = duckhts_bcf_decode_status(
        reader_name, "FORMAT", tag, header, record, ret, error, error_size);
    if (status == DUCKHTS_BCF_DECODE_FATAL) return 0;
    if (status == DUCKHTS_BCF_DECODE_TYPE_MISMATCH) goto mismatch;
    if (ret <= 0) goto absent;
    if (numeric && !duckhts_bcf_check_format_width(reader_name, tag, header, record,
                                                  ret, samples, error, error_size)) goto mismatch;
    if (numeric && !duckhts_bcf_check_scalar_count(header, record,
            is_gt ? DUCKHTS_BCF_FIELD_GT : DUCKHTS_BCF_FIELD_FORMAT,
            id, header_type, values->data, ret, reader_name, error, error_size)) goto mismatch;
    values->count = ret;
    values->stride = numeric ? ret / samples : 0;
absent:
    values->loaded = 1;
    return 1;
mismatch:
    if (!accept_mismatch(policy, error)) return 0;
    goto absent;
}
