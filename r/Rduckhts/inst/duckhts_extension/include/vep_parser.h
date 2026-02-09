/**
 * vep_parser.h - VEP/SnpEff/BCSQ Annotation Parser
 *
 * Parses structured annotation fields (CSQ, BCSQ, ANN) from VCF headers
 * and records.  Type inference follows bcftools split-vep conventions.
 * Ported from RBCFTools.
 */

#ifndef VEP_PARSER_H
#define VEP_PARSER_H

#include <htslib/vcf.h>

#ifdef __cplusplus
extern "C" {
#endif

#define VEP_MAX_FIELDS 256
#define VEP_MAX_FIELD_NAME 128
#define VEP_TAG_CSQ   "CSQ"
#define VEP_TAG_BCSQ  "BCSQ"
#define VEP_TAG_ANN   "ANN"

typedef enum {
    VEP_TYPE_STRING  = BCF_HT_STR,
    VEP_TYPE_INTEGER = BCF_HT_INT,
    VEP_TYPE_FLOAT   = BCF_HT_REAL,
    VEP_TYPE_FLAG    = BCF_HT_FLAG
} vep_field_type_t;

typedef struct {
    char* name;
    vep_field_type_t type;
    int index;
    int is_list;
} vep_field_t;

typedef struct vep_schema_t {
    char* tag_name;
    int n_fields;
    vep_field_t* fields;
    int header_id;
} vep_schema_t;

typedef struct {
    char* str_value;
    int32_t int_value;
    float float_value;
    int is_missing;
} vep_value_t;

typedef struct {
    int n_values;
    vep_value_t* values;
} vep_transcript_t;

typedef struct {
    int n_transcripts;
    vep_transcript_t* transcripts;
} vep_record_t;

typedef struct {
    const char* tag;
    const char* columns;
    int transcript_mode;
} vep_options_t;

void vep_options_init(vep_options_t* opts);
vep_schema_t* vep_schema_parse(const bcf_hdr_t* hdr, const char* tag);
void vep_schema_destroy(vep_schema_t* schema);
int vep_schema_get_field_index(const vep_schema_t* schema, const char* name);
const vep_field_t* vep_schema_get_field(const vep_schema_t* schema, int index);
vep_field_type_t vep_infer_type(const char* field_name);
const char* vep_type_name(vep_field_type_t type);
vep_record_t* vep_record_parse(const vep_schema_t* schema, const char* csq_value);
vep_record_t* vep_record_parse_bcf(const vep_schema_t* schema, const bcf_hdr_t* hdr, bcf1_t* rec);
void vep_record_destroy(vep_record_t* record);
const vep_value_t* vep_record_get_value(const vep_record_t* record, int transcript_idx, int field_idx);
const char* vep_detect_tag(const bcf_hdr_t* hdr);
int vep_has_annotation(const bcf_hdr_t* hdr);

#ifdef __cplusplus
}
#endif

#endif /* VEP_PARSER_H */
