/* Typed genotype mechanics with no DuckDB or R dependency. */
#include "include/bcf_genotypes.h"
#include <htslib/hts_log.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void duckhts_bcf_genotypes_destroy(duckhts_bcf_genotypes_t *values) {
    free(values->gt);
    free(values->ps);
    memset(values, 0, sizeof(*values));
}

int duckhts_bcf_genotype_ploidy(const int32_t *gt, int stride) {
    int ploidy = 0;
    while (ploidy < stride && gt[ploidy] != bcf_int32_vector_end) ploidy++;
    return ploidy;
}

int duckhts_bcf_genotype_has_alt(const int32_t *gt, int ploidy) {
    for (int i = 0; i < ploidy; i++) {
        if (!bcf_gt_is_missing(gt[i]) && bcf_gt_allele(gt[i]) > 0) return 1;
    }
    return 0;
}

static int accept_mismatch(duckhts_bcf_decode_policy_t policy, const char *error) {
    if (policy == DUCKHTS_BCF_DECODE_ERROR) return 0;
    if (policy == DUCKHTS_BCF_DECODE_WARN) hts_log_warning("%s", error);
    return 1;
}

static int decode_field(bcf_hdr_t *header, bcf1_t *record, int is_gt,
                         int32_t **data, int *capacity, int *count,
                         duckhts_bcf_decode_policy_t policy, char *error, size_t error_size) {
    const char *tag = is_gt ? "GT" : "PS";
    *count = 0;
    int id = bcf_hdr_id2int(header, BCF_DT_ID, tag);
    if (id < 0 || !bcf_hdr_idinfo_exists(header, BCF_HL_FMT, id)) return 1;
    if (bcf_hdr_id2type(header, BCF_HL_FMT, id) != (is_gt ? BCF_HT_STR : BCF_HT_INT) ||
        bcf_hdr_id2length(header, BCF_HL_FMT, id) != BCF_VL_FIXED ||
        bcf_hdr_id2number(header, BCF_HL_FMT, id) != 1) {
        snprintf(error, error_size, "read_geno: FORMAT/%s must have Number=1,Type=%s",
                 tag, is_gt ? "String" : "Integer");
        return accept_mismatch(policy, error);
    }
    if (!duckhts_bcf_check_field_type(header, record,
            is_gt ? DUCKHTS_BCF_FIELD_GT : DUCKHTS_BCF_FIELD_FORMAT,
            id, is_gt ? BCF_HT_STR : BCF_HT_INT, "read_geno", error, error_size)) {
        return accept_mismatch(policy, error);
    }
    bcf_fmt_t *format = bcf_get_fmt_id(record, id);
    if (!format || !format->p) return 1;
    int samples = bcf_hdr_nsamples(header);
    if (format->n < 0 || format->n > INT_MAX / samples ||
        (size_t)format->n > SIZE_MAX / sizeof(int32_t) / (size_t)samples) {
        snprintf(error, error_size, "read_geno: FORMAT/%s exceeds the supported decoded-value capacity", tag);
        return 0;
    }
    int32_t *previous = *data;
    int previous_capacity = *capacity;
    int ret = is_gt ? bcf_get_genotypes(header, record, data, capacity)
        : bcf_get_format_int32(header, record, tag, data, capacity);
    /* HTSlib 1.24 overwrites *data with realloc's return. On failure the old
     * allocation is still ours; retain it so worker teardown can release it. */
    if (ret == -4 && !*data) {
        *data = previous;
        *capacity = previous_capacity;
    }
    duckhts_bcf_decode_status_t status = duckhts_bcf_decode_status(
        "read_geno", "FORMAT", tag, header, record, ret, error, error_size);
    if (status == DUCKHTS_BCF_DECODE_FATAL) return 0;
    if (status == DUCKHTS_BCF_DECODE_TYPE_MISMATCH) return accept_mismatch(policy, error);
    if (ret <= 0) return 1;
    if (!duckhts_bcf_check_format_width("read_geno", tag, header, record,
                                         ret, samples, error, error_size)) {
        return accept_mismatch(policy, error);
    }
    if (!is_gt && ret != samples) {
        snprintf(error, error_size,
                 "read_geno: FORMAT/PS requires one value per sample at %s:%lld (got %d for %d samples)",
                 bcf_hdr_id2name(header, record->rid), (long long)record->pos + 1, ret, samples);
        return accept_mismatch(policy, error);
    }
    if (is_gt) {
        int stride = ret / samples;
        for (int sample = 0; sample < samples; sample++) {
            for (int slot = 0; slot < stride; slot++) {
                int32_t allele = (*data)[(size_t)sample * stride + slot];
                if (allele == bcf_int32_vector_end) break;
                if (allele < 0 || (!bcf_gt_is_missing(allele) &&
                                            bcf_gt_allele(allele) >= record->n_allele)) {
                    snprintf(error, error_size,
                        "read_geno: invalid FORMAT/GT allele at %s:%lld sample %s slot %d",
                        bcf_hdr_id2name(header, record->rid), (long long)record->pos + 1,
                        header->samples[sample], slot);
                    return accept_mismatch(policy, error);
                }
            }
        }
    }
    *count = ret;
    return 1;
}

int duckhts_bcf_genotypes_decode(duckhts_bcf_genotypes_t *values, bcf_hdr_t *header,
                                 bcf1_t *record, duckhts_bcf_decode_policy_t policy,
                                 char *error, size_t error_size) {
    values->gt_stride = 0;
    values->ps_present = 0;
    values->samples = bcf_hdr_nsamples(header);
    if (!values->samples) return 1;
    if (bcf_unpack(record, BCF_UN_FMT) < 0 || record->n_sample != (unsigned)values->samples) {
        snprintf(error, error_size, "read_geno: failed to unpack FORMAT or sample count differs from header");
        return 0;
    }
    int gt_count, ps_count;
    if (!decode_field(header, record, 1, &values->gt, &values->gt_capacity, &gt_count,
                       policy, error, error_size) ||
        !decode_field(header, record, 0, &values->ps, &values->ps_capacity, &ps_count,
                       policy, error, error_size)) return 0;
    values->gt_stride = gt_count / values->samples;
    values->ps_present = ps_count > 0;
    return 1;
}
