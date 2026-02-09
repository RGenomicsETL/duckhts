/**
 * VCF Type System for DuckHTS Extension
 * 
 * VCF-spec compliant type definitions matching the nanoarrow
 * vcf_arrow_stream implementation for type consistency.
 * 
 * Based on htslib's bcf_hdr_check_sanity() and VCF specification.
 * Ported from RBCFTools bcf_reader extension.
 */

#ifndef VCF_TYPES_H
#define VCF_TYPES_H

#include <htslib/vcf.h>
#include <string.h>
#include <stdio.h>

/* ================================================================
 * VCF Field Specification (matching nanoarrow's vcf_fmt_spec_t)
 * ================================================================ */

typedef struct {
    const char* name;
    const char* number_str;
    int vl_type;
    int count;
    int type;
} vcf_field_spec_t;

static const char* vcf_type_names[] = {"Flag", "Integer", "Float", "String"};

/* ================================================================
 * Standard FORMAT Field Definitions
 * ================================================================ */

static const vcf_field_spec_t VCF_FORMAT_SPECS[] = {
    {"AD",   "R", BCF_VL_R,     0, BCF_HT_INT},
    {"ADF",  "R", BCF_VL_R,     0, BCF_HT_INT},
    {"ADR",  "R", BCF_VL_R,     0, BCF_HT_INT},
    {"EC",   "A", BCF_VL_A,     0, BCF_HT_INT},
    {"GL",   "G", BCF_VL_G,     0, BCF_HT_REAL},
    {"GP",   "G", BCF_VL_G,     0, BCF_HT_REAL},
    {"PL",   "G", BCF_VL_G,     0, BCF_HT_INT},
    {"PP",   "G", BCF_VL_G,     0, BCF_HT_INT},
    {"DP",   "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"LEN",  "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"FT",   "1", BCF_VL_FIXED, 1, BCF_HT_STR},
    {"GQ",   "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"GT",   "1", BCF_VL_FIXED, 1, BCF_HT_STR},
    {"HQ",   "2", BCF_VL_FIXED, 2, BCF_HT_INT},
    {"MQ",   "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"PQ",   "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"PS",   "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {NULL,   NULL, 0,           0, 0}
};

/* ================================================================
 * Standard INFO Field Definitions
 * ================================================================ */

static const vcf_field_spec_t VCF_INFO_SPECS[] = {
    {"AD",        "R", BCF_VL_R,     0, BCF_HT_INT},
    {"ADF",       "R", BCF_VL_R,     0, BCF_HT_INT},
    {"ADR",       "R", BCF_VL_R,     0, BCF_HT_INT},
    {"AC",        "A", BCF_VL_A,     0, BCF_HT_INT},
    {"AF",        "A", BCF_VL_A,     0, BCF_HT_REAL},
    {"CIGAR",     "A", BCF_VL_A,     0, BCF_HT_STR},
    {"AA",        "1", BCF_VL_FIXED, 1, BCF_HT_STR},
    {"AN",        "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"BQ",        "1", BCF_VL_FIXED, 1, BCF_HT_REAL},
    {"DB",        "0", BCF_VL_FIXED, 0, BCF_HT_FLAG},
    {"DP",        "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"END",       "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"H2",        "0", BCF_VL_FIXED, 0, BCF_HT_FLAG},
    {"H3",        "0", BCF_VL_FIXED, 0, BCF_HT_FLAG},
    {"MQ",        "1", BCF_VL_FIXED, 1, BCF_HT_REAL},
    {"MQ0",       "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"NS",        "1", BCF_VL_FIXED, 1, BCF_HT_INT},
    {"SB",        "4", BCF_VL_FIXED, 4, BCF_HT_INT},
    {"SOMATIC",   "0", BCF_VL_FIXED, 0, BCF_HT_FLAG},
    {"VALIDATED", "0", BCF_VL_FIXED, 0, BCF_HT_FLAG},
    {"1000G",     "0", BCF_VL_FIXED, 0, BCF_HT_FLAG},
    {NULL,        NULL, 0,           0, 0}
};

/* ================================================================
 * Lookup Functions
 * ================================================================ */

static inline const vcf_field_spec_t* vcf_lookup_format_spec(const char* name) {
    for (int i = 0; VCF_FORMAT_SPECS[i].name != NULL; i++) {
        if (strcmp(name, VCF_FORMAT_SPECS[i].name) == 0)
            return &VCF_FORMAT_SPECS[i];
    }
    return NULL;
}

static inline const vcf_field_spec_t* vcf_lookup_info_spec(const char* name) {
    for (int i = 0; VCF_INFO_SPECS[i].name != NULL; i++) {
        if (strcmp(name, VCF_INFO_SPECS[i].name) == 0)
            return &VCF_INFO_SPECS[i];
    }
    return NULL;
}

/* ================================================================
 * Validation
 * ================================================================ */

static inline int vcf_check_number(const vcf_field_spec_t* spec, int header_vl_type) {
    if (!spec) return 0;
    if (spec->vl_type == BCF_VL_FIXED)
        return (header_vl_type != BCF_VL_FIXED);
    return (header_vl_type != spec->vl_type && header_vl_type != BCF_VL_VAR) ? 1 : 0;
}

static inline int vcf_check_type(const vcf_field_spec_t* spec, int header_type) {
    if (!spec) return 0;
    return (header_type != spec->type);
}

/* ================================================================
 * Warning Callback
 * ================================================================ */

typedef void (*vcf_warning_func_t)(const char* msg, void* ctx);

static vcf_warning_func_t g_vcf_warning_func = NULL;
static void* g_vcf_warning_ctx = NULL;

static inline void vcf_set_warning_callback(vcf_warning_func_t func, void* ctx) {
    g_vcf_warning_func = func;
    g_vcf_warning_ctx = ctx;
}

static inline void vcf_emit_warning(const char* msg) {
    if (g_vcf_warning_func)
        g_vcf_warning_func(msg, g_vcf_warning_ctx);
    else
        fprintf(stderr, "Warning: %s\n", msg);
}

static inline int vcf_validate_format_field(const char* field_name,
                                            int header_vl_type,
                                            int header_type,
                                            int* corrected_type) {
    const vcf_field_spec_t* spec = vcf_lookup_format_spec(field_name);
    int corrected_vl_type = header_vl_type;
    *corrected_type = header_type;

    if (spec) {
        if (vcf_check_number(spec, header_vl_type)) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "FORMAT/%s should be Number=%s per VCF spec; correcting schema",
                     field_name, spec->number_str);
            vcf_emit_warning(msg);
            corrected_vl_type = spec->vl_type;
        }
        if (vcf_check_type(spec, header_type)) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "FORMAT/%s should be Type=%s per VCF spec, but header declares Type=%s; using header type",
                     field_name, vcf_type_names[spec->type], vcf_type_names[header_type]);
            vcf_emit_warning(msg);
        }
    }
    return corrected_vl_type;
}

static inline int vcf_validate_info_field(const char* field_name,
                                          int header_vl_type,
                                          int header_type,
                                          int* corrected_type) {
    const vcf_field_spec_t* spec = vcf_lookup_info_spec(field_name);
    int corrected_vl_type = header_vl_type;
    *corrected_type = header_type;

    if (spec) {
        if (vcf_check_number(spec, header_vl_type)) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "INFO/%s should be Number=%s per VCF spec; correcting schema",
                     field_name, spec->number_str);
            vcf_emit_warning(msg);
            corrected_vl_type = spec->vl_type;
        }
        if (vcf_check_type(spec, header_type)) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "INFO/%s should be Type=%s per VCF spec, but header declares Type=%s; using header type",
                     field_name, vcf_type_names[spec->type], vcf_type_names[header_type]);
            vcf_emit_warning(msg);
        }
    }
    return corrected_vl_type;
}

/* ================================================================
 * Type Mapping Utilities
 * ================================================================ */

static inline int vcf_is_list_type(int vl_type) {
    return (vl_type != BCF_VL_FIXED);
}

static inline int vcf_get_expected_count(int vl_type, int n_allele, int ploidy) {
    (void)ploidy;
    switch (vl_type) {
        case BCF_VL_FIXED: return 1;
        case BCF_VL_VAR:   return -1;
        case BCF_VL_A:     return n_allele - 1;
        case BCF_VL_G:     return (n_allele * (n_allele + 1)) / 2;
        case BCF_VL_R:     return n_allele;
        default:           return -1;
    }
}

#endif /* VCF_TYPES_H */
