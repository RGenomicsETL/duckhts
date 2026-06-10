/**
 * vep_parser.c - VEP/SnpEff/BCSQ Annotation Parser (self-contained)
 *
 * Minimal copy of the RBCFTools parser for use inside the DuckDB extension.
 * Follows bcftools split-vep type inference and parses pipe-delimited CSQ/BCSQ/ANN.
 */

#include "include/vep_parser.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <regex.h>

typedef struct {
    const char* pattern;
    vep_field_type_t type;
} vep_column_type_rule_t;

typedef struct {
    regex_t regex;
    vep_field_type_t type;
    int compiled;
} vep_compiled_rule_t;

static void vep_destroy_custom_type_rules(vep_compiled_rule_t* rules, int count);

static const vep_column_type_rule_t VEP_CSQ_BUILTIN_TYPE_RULES[] = {
    {"DISTANCE",                  VEP_TYPE_INTEGER},
    {"STRAND",                    VEP_TYPE_INTEGER},
    {"TSL",                       VEP_TYPE_INTEGER},
    {"GENE_PHENO",                VEP_TYPE_INTEGER},
    {"HGVS_OFFSET",               VEP_TYPE_INTEGER},
    {".*_POPS",                   VEP_TYPE_STRING},
    {"AF",                        VEP_TYPE_FLOAT},
    {".*_AF",                     VEP_TYPE_FLOAT},
    {"MAX_AF_.*",                 VEP_TYPE_FLOAT},
    {"MOTIF_POS",                 VEP_TYPE_INTEGER},
    {"MOTIF_SCORE_CHANGE",        VEP_TYPE_FLOAT},
    {"existing_InFrame_oORFs",    VEP_TYPE_INTEGER},
    {"existing_OutOfFrame_oORFs", VEP_TYPE_INTEGER},
    {"existing_uORFs",            VEP_TYPE_INTEGER},
    {"SpliceAI_pred_DP_.*",       VEP_TYPE_INTEGER},
    {"SpliceAI_pred_DS_.*",       VEP_TYPE_FLOAT},
    {"SpliceAI_pred_SYMBOL",      VEP_TYPE_STRING},
    {"SpliceAI_pred",             VEP_TYPE_STRING},
    {"SpliceAI_cutoff",           VEP_TYPE_STRING},
    {"CADD_PHRED",                VEP_TYPE_FLOAT},
    {"CADD_RAW",                  VEP_TYPE_FLOAT},
    {"REVEL",                     VEP_TYPE_FLOAT},
    {"am_pathogenicity",          VEP_TYPE_FLOAT},
    {"am_class",                  VEP_TYPE_STRING},
    {"am_protein_variant",        VEP_TYPE_STRING},
    {"am_uniprot_id",             VEP_TYPE_STRING},
    {"am_transcript_id",          VEP_TYPE_STRING},
    {"SIFT",                      VEP_TYPE_STRING},
    {"PolyPhen",                  VEP_TYPE_STRING},
    {"PolyPhen_2",                VEP_TYPE_STRING},
    {NULL,                        VEP_TYPE_STRING}
};

static const vep_column_type_rule_t VEP_ANN_BUILTIN_TYPE_RULES[] = {
    {"Distance", VEP_TYPE_INTEGER},
    {NULL,       VEP_TYPE_STRING}
};

static void* vep_malloc(size_t size) { return malloc(size); }
static void vep_free(void* ptr) { free(ptr); }

static char* vep_strdup(const char* s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char* copy = (char*)vep_malloc(len);
    if (copy) memcpy(copy, s, len);
    return copy;
}

void vep_options_init(vep_options_t* opts) {
    if (!opts) return;
    opts->tag = NULL;
    opts->columns = NULL;
    opts->additional_csq_column_types = NULL;
    opts->transcript_mode = 1; // first
}

static const char* parse_format_string(const char* description, int* n_fields) {
    const char* format = strstr(description, "Format: ");
    if (!format) return NULL;
    format += strlen("Format: ");
    const char* end = strchr(format, '"');
    if (!end) end = format + strlen(format);
    *n_fields = 1;
    for (const char* p = format; p < end; p++) {
        if (*p == '|') (*n_fields)++;
    }
    return format;
}

static char** split_format_fields(const char* format, int n_fields) {
    char** fields = (char**)vep_malloc(n_fields * sizeof(char*));
    if (!fields) return NULL;
    const char* start = format;
    int idx = 0;
    for (const char* p = format; ; p++) {
        if (*p == '|' || *p == '\0' || *p == '"') {
            size_t len = p - start;
            char* name = (char*)vep_malloc(len + 1);
            if (!name) {
                for (int i = 0; i < idx; i++) vep_free(fields[i]);
                vep_free(fields);
                return NULL;
            }
            memcpy(name, start, len);
            name[len] = '\0';
            fields[idx++] = name;
            if (*p == '\0' || *p == '"') break;
            start = p + 1;
        }
    }
    return fields;
}

static int vep_compile_rule(regex_t* regex, const char* pattern) {
    size_t pattern_len = strlen(pattern);
    char* anchored = (char*)vep_malloc(pattern_len + 3);
    if (!anchored) return -1;
    anchored[0] = '^';
    memcpy(anchored + 1, pattern, pattern_len);
    anchored[pattern_len + 1] = '$';
    anchored[pattern_len + 2] = '\0';
    int ret = regcomp(regex, anchored, REG_NOSUB | REG_EXTENDED);
    vep_free(anchored);
    return ret;
}

static const vep_column_type_rule_t* vep_builtin_type_rules(const char* tag_name) {
    if (tag_name && strcmp(tag_name, VEP_TAG_ANN) == 0) return VEP_ANN_BUILTIN_TYPE_RULES;
    return VEP_CSQ_BUILTIN_TYPE_RULES;
}

static vep_field_type_t vep_parse_type_name(const char* type_name, int len, int* ok) {
    *ok = 1;
    if ((len == 3 && strncmp(type_name, "Str", 3) == 0) ||
        (len == 6 && strncmp(type_name, "String", 6) == 0)) {
        return VEP_TYPE_STRING;
    }
    if ((len == 3 && strncmp(type_name, "Int", 3) == 0) ||
        (len == 7 && strncmp(type_name, "Integer", 7) == 0)) {
        return VEP_TYPE_INTEGER;
    }
    if ((len == 4 && strncmp(type_name, "Real", 4) == 0) ||
        (len == 5 && strncmp(type_name, "Float", 5) == 0)) {
        return VEP_TYPE_FLOAT;
    }
    if (len == 4 && strncmp(type_name, "Flag", 4) == 0) {
        return VEP_TYPE_FLAG;
    }
    *ok = 0;
    return VEP_TYPE_STRING;
}

static int vep_parse_custom_type_rules(const char* rules_text,
                                       vep_compiled_rule_t** out_rules,
                                       int* out_count) {
    *out_rules = NULL;
    *out_count = 0;
    if (!rules_text || !*rules_text) return 0;

    char* copy = vep_strdup(rules_text);
    if (!copy) return -1;

    char* cursor = copy;
    while (*cursor) {
        while (*cursor == '\n' || *cursor == '\r' || *cursor == ';') cursor++;
        if (!*cursor) break;

        char* line_end = cursor;
        while (*line_end && *line_end != '\n' && *line_end != '\r' && *line_end != ';') line_end++;
        char saved = *line_end;
        *line_end = '\0';

        while (*cursor && isspace((unsigned char)*cursor)) cursor++;
        if (*cursor && *cursor != '#') {
            char* pattern = cursor;
            char* type_start = pattern;
            while (*type_start && !isspace((unsigned char)*type_start)) type_start++;
            if (!*type_start) {
                vep_destroy_custom_type_rules(*out_rules, *out_count);
                *out_rules = NULL;
                *out_count = 0;
                vep_free(copy);
                return -1;
            }
            *type_start++ = '\0';
            while (*type_start && isspace((unsigned char)*type_start)) type_start++;
            if (!*type_start) {
                vep_destroy_custom_type_rules(*out_rules, *out_count);
                *out_rules = NULL;
                *out_count = 0;
                vep_free(copy);
                return -1;
            }
            char* type_end = type_start;
            while (*type_end && !isspace((unsigned char)*type_end)) type_end++;
            int ok = 0;
            vep_field_type_t type = vep_parse_type_name(type_start, (int)(type_end - type_start), &ok);
            if (!ok) {
                vep_destroy_custom_type_rules(*out_rules, *out_count);
                *out_rules = NULL;
                *out_count = 0;
                vep_free(copy);
                return -1;
            }

            vep_compiled_rule_t* rules = (vep_compiled_rule_t*)realloc(*out_rules, ((*out_count) + 1) * sizeof(**out_rules));
            if (!rules) {
                vep_destroy_custom_type_rules(*out_rules, *out_count);
                *out_rules = NULL;
                *out_count = 0;
                vep_free(copy);
                return -1;
            }
            *out_rules = rules;
            (*out_rules)[*out_count].type = type;
            (*out_rules)[*out_count].compiled = 0;
            if (vep_compile_rule(&(*out_rules)[*out_count].regex, pattern) != 0) {
                vep_destroy_custom_type_rules(*out_rules, *out_count);
                *out_rules = NULL;
                *out_count = 0;
                vep_free(copy);
                return -1;
            }
            (*out_rules)[*out_count].compiled = 1;
            (*out_count)++;
        }

        *line_end = saved;
        cursor = line_end;
    }

    vep_free(copy);
    return 0;
}

static void vep_destroy_custom_type_rules(vep_compiled_rule_t* rules, int count) {
    if (!rules) return;
    for (int i = 0; i < count; i++) {
        if (rules[i].compiled) regfree(&rules[i].regex);
    }
    free(rules);
}

int vep_validate_column_type_rules(const char* rules_text, char* errbuf, size_t errbuf_size) {
    vep_compiled_rule_t* rules = NULL;
    int rule_count = 0;
    int ok = vep_parse_custom_type_rules(rules_text, &rules, &rule_count) == 0;
    vep_destroy_custom_type_rules(rules, rule_count);
    if (!ok && errbuf && errbuf_size > 0) {
        snprintf(errbuf, errbuf_size,
                 "additional_csq_column_types must use bcftools-style 'PATTERN TYPE' entries separated by newlines or ';'");
    }
    return ok ? 1 : 0;
}

vep_field_type_t vep_infer_type(const char* tag_name, const char* field_name, const char* additional_csq_column_types) {
    if (!field_name) return VEP_TYPE_STRING;

    vep_compiled_rule_t* custom_rules = NULL;
    int n_custom_rules = 0;
    if (vep_parse_custom_type_rules(additional_csq_column_types, &custom_rules, &n_custom_rules) == 0) {
        for (int i = 0; i < n_custom_rules; i++) {
            if (regexec(&custom_rules[i].regex, field_name, 0, NULL, 0) == 0) {
                vep_field_type_t type = custom_rules[i].type;
                vep_destroy_custom_type_rules(custom_rules, n_custom_rules);
                return type;
            }
        }
    }
    vep_destroy_custom_type_rules(custom_rules, n_custom_rules);

    const vep_column_type_rule_t* rules = vep_builtin_type_rules(tag_name);
    for (int i = 0; rules[i].pattern != NULL; i++) {
        regex_t regex;
        if (vep_compile_rule(&regex, rules[i].pattern) != 0) continue;
        int matched = regexec(&regex, field_name, 0, NULL, 0) == 0;
        regfree(&regex);
        if (matched) return rules[i].type;
    }

    return VEP_TYPE_STRING;
}

const char* vep_type_name(vep_field_type_t type) {
    switch (type) {
        case VEP_TYPE_INTEGER: return "Integer";
        case VEP_TYPE_FLOAT:   return "Float";
        case VEP_TYPE_FLAG:    return "Flag";
        case VEP_TYPE_STRING:
        default:               return "String";
    }
}

const char* vep_detect_tag(const bcf_hdr_t* hdr) {
    if (!hdr) return NULL;
    const char* tags[] = {
        VEP_TAG_CSQ,
        VEP_TAG_BCSQ,
        VEP_TAG_ANN,
        VEP_TAG_VEP,
        VEP_TAG_vep,
        NULL
    };
    for (int i = 0; tags[i]; i++) {
        int id = bcf_hdr_id2int(hdr, BCF_DT_ID, tags[i]);
        if (id >= 0 && bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id)) {
            return tags[i];
        }
    }
    return NULL;
}

int vep_has_annotation(const bcf_hdr_t* hdr) {
    return vep_detect_tag(hdr) != NULL;
}

vep_schema_t* vep_schema_parse(const bcf_hdr_t* hdr, const char* tag, const char* additional_csq_column_types) {
    if (!hdr) return NULL;
    const char* detected_tag = tag ? tag : vep_detect_tag(hdr);
    if (!detected_tag) return NULL;
    
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, detected_tag);
    if (id < 0 || !bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id)) return NULL;
    
    bcf_hrec_t* hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", detected_tag, NULL);
    if (!hrec) return NULL;
    
    const char* description = NULL;
    for (int i = 0; i < hrec->nkeys; i++) {
        if (strcmp(hrec->keys[i], "Description") == 0) {
            description = hrec->vals[i];
            break;
        }
    }
    if (!description) return NULL;
    
    int n_fields = 0;
    const char* format = parse_format_string(description, &n_fields);
    if (!format || n_fields == 0 || n_fields > VEP_MAX_FIELDS) return NULL;
    
    char** field_names = split_format_fields(format, n_fields);
    if (!field_names) return NULL;
    
    vep_schema_t* schema = (vep_schema_t*)vep_malloc(sizeof(vep_schema_t));
    if (!schema) {
        for (int i = 0; i < n_fields; i++) vep_free(field_names[i]);
        vep_free(field_names);
        return NULL;
    }
    
    schema->tag_name = vep_strdup(detected_tag);
    schema->n_fields = n_fields;
    schema->header_id = id;
    schema->fields = (vep_field_t*)vep_malloc(n_fields * sizeof(vep_field_t));
    
    if (!schema->fields) {
        vep_free(schema->tag_name);
        vep_free(schema);
        for (int i = 0; i < n_fields; i++) vep_free(field_names[i]);
        vep_free(field_names);
        return NULL;
    }
    
    for (int i = 0; i < n_fields; i++) {
        schema->fields[i].name = field_names[i];
        schema->fields[i].type = vep_infer_type(detected_tag, field_names[i], additional_csq_column_types);
        schema->fields[i].index = i;
        schema->fields[i].is_list = (strcmp(field_names[i], "Consequence") == 0 ||
                                     strcmp(field_names[i], "Annotation") == 0 ||
                                     strcmp(field_names[i], "FLAGS") == 0 ||
                                     strcmp(field_names[i], "CLIN_SIG") == 0);
    }
    vep_free(field_names);
    return schema;
}

void vep_schema_destroy(vep_schema_t* schema) {
    if (!schema) return;
    if (schema->fields) {
        for (int i = 0; i < schema->n_fields; i++) vep_free(schema->fields[i].name);
        vep_free(schema->fields);
    }
    vep_free(schema->tag_name);
    vep_free(schema);
}

int vep_schema_get_field_index(const vep_schema_t* schema, const char* name) {
    if (!schema || !name) return -1;
    for (int i = 0; i < schema->n_fields; i++) {
        if (schema->fields[i].name && strcmp(schema->fields[i].name, name) == 0) return i;
    }
    return -1;
}

const vep_field_t* vep_schema_get_field(const vep_schema_t* schema, int index) {
    if (!schema || index < 0 || index >= schema->n_fields) return NULL;
    return &schema->fields[index];
}

int vep_parse_int(const char* str, int32_t* result) {
    if (!str || !*str || strcmp(str, ".") == 0) {
        *result = INT32_MIN;
        return 0;
    }
    char* endptr;
    long val = strtol(str, &endptr, 10);
    if (endptr == str || *endptr != '\0') {
        *result = INT32_MIN;
        return -1;
    }
    *result = (int32_t)val;
    return 1;
}

int vep_parse_float(const char* str, float* result) {
    if (!str || !*str || strcmp(str, ".") == 0) {
        *result = NAN;
        return 0;
    }
    char* endptr;
    double val = strtod(str, &endptr);
    if (endptr == str || *endptr != '\0') {
        *result = NAN;
        return -1;
    }
    *result = (float)val;
    return 1;
}

static int count_char(const char* s, char c) {
    int count = 0;
    while (*s) { if (*s == c) count++; s++; }
    return count;
}

static int vep_field_selected(int field_idx, const int* selected_indices, int n_selected) {
    if (!selected_indices || n_selected <= 0) return 1;
    for (int i = 0; i < n_selected; i++) {
        if (selected_indices[i] == field_idx) return 1;
    }
    return 0;
}

static vep_transcript_t* parse_single_transcript_selected(const vep_schema_t* schema,
                                                          const char* transcript_str,
                                                          const int* selected_indices,
                                                          int n_selected) {
    if (!schema || !transcript_str) return NULL;
    vep_transcript_t* transcript = (vep_transcript_t*)vep_malloc(sizeof(vep_transcript_t));
    if (!transcript) return NULL;
    transcript->n_values = schema->n_fields;
    transcript->values = (vep_value_t*)vep_malloc(schema->n_fields * sizeof(vep_value_t));
    if (!transcript->values) { vep_free(transcript); return NULL; }

    for (int i = 0; i < schema->n_fields; i++) {
        transcript->values[i].str_value = NULL;
        transcript->values[i].int_value = INT32_MIN;
        transcript->values[i].float_value = NAN;
        transcript->values[i].is_missing = 1;
    }

    char* copy = vep_strdup(transcript_str);
    if (!copy) { vep_free(transcript->values); vep_free(transcript); return NULL; }

    int field_idx = 0;
    char* token = copy;
    char* next_pipe;
    while (field_idx < schema->n_fields) {
        next_pipe = strchr(token, '|');
        if (next_pipe) *next_pipe = '\0';
        if (vep_field_selected(field_idx, selected_indices, n_selected)) {
            while (*token && isspace((unsigned char)*token)) token++;
            size_t token_len = strlen(token);
            while (token_len > 0 && isspace((unsigned char)token[token_len - 1])) {
                token[--token_len] = '\0';
            }

            if (token_len > 0 && strcmp(token, ".") != 0) {
                transcript->values[field_idx].str_value = vep_strdup(token);
                transcript->values[field_idx].is_missing = 0;
                const vep_field_t* field = &schema->fields[field_idx];
                if (field->type == VEP_TYPE_INTEGER) vep_parse_int(token, &transcript->values[field_idx].int_value);
                else if (field->type == VEP_TYPE_FLOAT) vep_parse_float(token, &transcript->values[field_idx].float_value);
            }
        }
        field_idx++;
        if (!next_pipe) break;
        token = next_pipe + 1;
    }
    vep_free(copy);
    return transcript;
}

static vep_transcript_t* parse_single_transcript(const vep_schema_t* schema, const char* transcript_str) {
    return parse_single_transcript_selected(schema, transcript_str, NULL, 0);
}

vep_record_t* vep_record_parse_selected(const vep_schema_t* schema, const char* csq_value,
                                         const int* selected_indices, int n_selected) {
    if (!schema || !csq_value || !*csq_value) return NULL;
    int n_transcripts = count_char(csq_value, ',') + 1;
    vep_record_t* record = (vep_record_t*)vep_malloc(sizeof(vep_record_t));
    if (!record) return NULL;
    record->n_transcripts = 0;
    record->transcripts = (vep_transcript_t*)vep_malloc(n_transcripts * sizeof(vep_transcript_t));
    if (!record->transcripts) { vep_free(record); return NULL; }

    char* copy = vep_strdup(csq_value);
    if (!copy) { vep_free(record->transcripts); vep_free(record); return NULL; }
    char* saveptr = NULL;
    char* token = strtok_r(copy, ",", &saveptr);
    while (token && record->n_transcripts < n_transcripts) {
        vep_transcript_t* transcript = parse_single_transcript_selected(schema, token, selected_indices, n_selected);
        if (transcript) {
            record->transcripts[record->n_transcripts] = *transcript;
            vep_free(transcript);
            record->n_transcripts++;
        }
        token = strtok_r(NULL, ",", &saveptr);
    }
    vep_free(copy);
    if (record->n_transcripts == 0) {
        vep_free(record->transcripts);
        vep_free(record);
        return NULL;
    }
    return record;
}

vep_record_t* vep_record_parse(const vep_schema_t* schema, const char* csq_value) {
    return vep_record_parse_selected(schema, csq_value, NULL, 0);
}

vep_record_t* vep_record_parse_bcf(const vep_schema_t* schema, const bcf_hdr_t* hdr, bcf1_t* rec) {
    return vep_record_parse_selected_bcf(schema, hdr, rec, NULL, 0);
}

vep_record_t* vep_record_parse_selected_bcf(const vep_schema_t* schema, const bcf_hdr_t* hdr,
                                             bcf1_t* rec, const int* selected_indices, int n_selected) {
    if (!schema || !hdr || !rec) return NULL;
    char* csq_value = NULL;
    int n_csq = 0;
    int ret = bcf_get_info_string(hdr, rec, schema->tag_name, &csq_value, &n_csq);
    if (ret <= 0 || !csq_value) { free(csq_value); return NULL; }
    vep_record_t* record = vep_record_parse_selected(schema, csq_value, selected_indices, n_selected);
    free(csq_value);
    return record;
}

void vep_record_destroy(vep_record_t* record) {
    if (!record) return;
    if (record->transcripts) {
        for (int i = 0; i < record->n_transcripts; i++) {
            if (record->transcripts[i].values) {
                for (int j = 0; j < record->transcripts[i].n_values; j++) {
                    vep_free(record->transcripts[i].values[j].str_value);
                }
                vep_free(record->transcripts[i].values);
            }
        }
        vep_free(record->transcripts);
    }
    vep_free(record);
}

const vep_value_t* vep_record_get_value(const vep_record_t* record, int transcript_idx, int field_idx) {
    if (!record) return NULL;
    if (transcript_idx < 0 || transcript_idx >= record->n_transcripts) return NULL;
    const vep_transcript_t* transcript = &record->transcripts[transcript_idx];
    if (field_idx < 0 || field_idx >= transcript->n_values) return NULL;
    return &transcript->values[field_idx];
}
