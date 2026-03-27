#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/khash_str2int.h>
#include <htslib/vcf.h>

#ifndef duckdb_malloc
#define duckdb_malloc malloc
#define duckdb_realloc realloc
#define duckdb_free free
#endif

#define SCORE_USE_GT 1
#define SCORE_USE_DS 2
#define SCORE_USE_HDS 3
#define SCORE_USE_AP 4
#define SCORE_USE_GP 5
#define SCORE_USE_AS 6

#define SCORE_HDR_SNP 0
#define SCORE_HDR_BP 1
#define SCORE_HDR_CHR 2
#define SCORE_HDR_A1 3
#define SCORE_HDR_P 5
#define SCORE_HDR_OR 7
#define SCORE_HDR_BETA 8
#define SCORE_HDR_LP 16
#define SCORE_HDR_SIZE 24

KHASH_MAP_INIT_INT64(score64, int)

typedef struct {
    char *hdr_str;
    int hdr_num;
} score_mapping_t;

typedef struct {
    int use_tag;
    int use_variant_id;
    int counts;
    char *bcf_path;
    char *summary_path;
    char *columns_preset;
    char *columns_file;
    double *q_thr_lp;
    int n_q_thr;
    char *prs_name;
} score_bind_t;

typedef struct {
    int n_samples;
    int n_metrics;
    char **sample_names;
    char **metric_names;
    double **metric_values;
    uint8_t *metric_is_count;
    int64_t row_offset;
} score_init_t;

typedef struct {
    char *a1;
    float es;
    float lp;
} score_marker_t;

typedef struct {
    int use_snp;
    int snp_ok;
    int chr_ok;
    int bp_ok;
    int a1_ok;
    int beta_ok;
    int p_ok;
    void *id2idx;
    void *rid_pos2idx;
    score_marker_t *markers;
    int n_markers;
    int m_markers;
} score_summary_t;

static const char *score_col_headers[SCORE_HDR_SIZE] = {
    "SNP", "BP", "CHR", "A1", "A2", "P", "Z", "OR", "BETA", "N", "N_CAS", "N_CON",
    "INFO", "FRQ", "A*", "SE", "LP", "AC", "NEFF", "NEFFDIV2", "HET_I2", "HET_P", "HET_LP", "DIRE"
};

static score_mapping_t score_plink_mapping[] = {
    {"SNP", SCORE_HDR_SNP}, {"BP", SCORE_HDR_BP}, {"CHR", SCORE_HDR_CHR}, {"A1", SCORE_HDR_A1},
    {"P", SCORE_HDR_P}, {"OR", SCORE_HDR_OR}, {"BETA", SCORE_HDR_BETA}
};

static score_mapping_t score_plink2_mapping[] = {
    {"ID", SCORE_HDR_SNP}, {"POS", SCORE_HDR_BP}, {"CHROM", SCORE_HDR_CHR}, {"A1", SCORE_HDR_A1},
    {"P", SCORE_HDR_P}, {"OR", SCORE_HDR_OR}, {"BETA", SCORE_HDR_BETA}, {"LOG10_P", SCORE_HDR_LP}
};

static score_mapping_t score_regenie_mapping[] = {
    {"ID", SCORE_HDR_SNP}, {"GENPOS", SCORE_HDR_BP}, {"CHROM", SCORE_HDR_CHR}, {"ALLELE1", SCORE_HDR_A1},
    {"BETA", SCORE_HDR_BETA}, {"LOG10P", SCORE_HDR_LP}
};

static score_mapping_t score_saige_mapping[] = {
    {"SNPID", SCORE_HDR_SNP}, {"markerID", SCORE_HDR_SNP}, {"POS", SCORE_HDR_BP}, {"CHR", SCORE_HDR_CHR},
    {"Allele2", SCORE_HDR_A1}, {"p.value", SCORE_HDR_P}, {"BETA", SCORE_HDR_BETA}
};

static score_mapping_t score_bolt_mapping[] = {
    {"SNP", SCORE_HDR_SNP}, {"BP", SCORE_HDR_BP}, {"CHR", SCORE_HDR_CHR}, {"ALLELE1", SCORE_HDR_A1},
    {"P_BOLT_LMM", SCORE_HDR_P}, {"P_BOLT_LMM_INF", SCORE_HDR_P}, {"P_LINREG", SCORE_HDR_P}, {"BETA", SCORE_HDR_BETA}
};

static score_mapping_t score_metal_mapping[] = {
    {"MarkerName", SCORE_HDR_SNP}, {"Position", SCORE_HDR_BP}, {"Chromosome", SCORE_HDR_CHR}, {"Allele1", SCORE_HDR_A1},
    {"P-value", SCORE_HDR_P}, {"Effect", SCORE_HDR_BETA}, {"log(P)", SCORE_HDR_LP}
};

static score_mapping_t score_pgs_mapping[] = {
    {"rsID", SCORE_HDR_SNP}, {"chr_position", SCORE_HDR_BP}, {"chr_name", SCORE_HDR_CHR}, {"effect_allele", SCORE_HDR_A1},
    {"OR", SCORE_HDR_OR}, {"HR", SCORE_HDR_OR}, {"effect_weight", SCORE_HDR_BETA}
};

static score_mapping_t score_ssf_mapping[] = {
    {"variant_id", SCORE_HDR_SNP}, {"rsid", SCORE_HDR_SNP}, {"rs_id", SCORE_HDR_SNP},
    {"base_pair_location", SCORE_HDR_BP}, {"chromosome", SCORE_HDR_CHR}, {"effect_allele", SCORE_HDR_A1},
    {"p_value", SCORE_HDR_P}, {"odds_ratio", SCORE_HDR_OR}, {"hazard_ratio", SCORE_HDR_OR}, {"beta", SCORE_HDR_BETA}
};

static int score_is_missing(float f) {
    return isnan(f) || bcf_float_is_missing(f) || bcf_float_is_vector_end(f);
}

/* Port of bcf_hdr_name2id_flexible() from bcftools +score (score.h).
 * Handles chr prefix strip/prepend, numeric sex chromosomes (23->X, 24->Y,
 * 26/MT/chrM), XY/XX/PAR1/PAR2 aliases. */
static int score_hdr_name2id_flexible(const bcf_hdr_t *hdr, const char *chr) {
    int rid;
    char buf[8];
    if (!chr || chr[0] == '\0') return -1;
    rid = bcf_hdr_name2id(hdr, chr);
    if (rid >= 0) return rid;
    /* strip chr prefix */
    if (strncmp(chr, "chr", 3) == 0) {
        rid = bcf_hdr_name2id(hdr, chr + 3);
        if (rid >= 0) return rid;
    }
    /* prepend chr prefix (up to 2-char contig name) */
    if (strlen(chr) <= 2) {
        buf[0] = 'c'; buf[1] = 'h'; buf[2] = 'r';
        buf[3] = chr[0]; buf[4] = chr[1]; buf[5] = '\0';
        rid = bcf_hdr_name2id(hdr, buf);
        if (rid >= 0) return rid;
    }
    /* numeric sex/mito aliases */
    if (strcmp(chr, "23") == 0 || strcmp(chr, "25") == 0 || strcmp(chr, "XY") == 0
        || strcmp(chr, "XX") == 0 || strcmp(chr, "PAR1") == 0 || strcmp(chr, "PAR2") == 0) {
        rid = bcf_hdr_name2id(hdr, "X");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrX");
    } else if (strcmp(chr, "24") == 0) {
        rid = bcf_hdr_name2id(hdr, "Y");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrY");
    } else if (strcmp(chr, "26") == 0 || strcmp(chr, "MT") == 0 || strcmp(chr, "chrM") == 0) {
        rid = bcf_hdr_name2id(hdr, "MT");
        if (rid >= 0) return rid;
        rid = bcf_hdr_name2id(hdr, "chrM");
    }
    return rid;
}

/* Parse a float from a field, treating "." and "NA" as NAN (matches upstream
 * tsv_read_float semantics from score.h). Returns 0 on success, -1 on error. */
static int score_parse_float(const char *s, float *out) {
    char *end = NULL;
    if (!s || *s == '\0') return -1;
    if (*s == '.' && (s[1] == '\0' || s[1] == '\t' || s[1] == '\n')) { *out = NAN; return 0; }
    if (s[0] == 'N' && s[1] == 'A' && (s[2] == '\0' || s[2] == '\t' || s[2] == '\n')) { *out = NAN; return 0; }
    *out = strtof(s, &end);
    if (end == s || *end != '\0') return -1;
    return 0;
}

/* Parse a p-value and compute -log10(p) with improved numerical stability.
 * Port of tsv_read_float_and_minus_log10() from bcftools +score (score.h).
 * For scientific notation, splits at e/E to compute -log10(mantissa) - exponent,
 * avoiding underflow with extremely small p-values (e.g. 5e-324). */
static int score_parse_pvalue_lp(const char *s, float *out) {
    const char *eptr;
    char buf[64];
    char *end = NULL;
    float mantissa, exponent;
    size_t mlen;

    if (!s || *s == '\0') return -1;
    if (*s == '.' && (s[1] == '\0' || s[1] == '\t' || s[1] == '\n')) { *out = NAN; return 0; }
    if (s[0] == 'N' && s[1] == 'A' && (s[2] == '\0' || s[2] == '\t' || s[2] == '\n')) { *out = NAN; return 0; }

    /* find 'e' or 'E' */
    for (eptr = s; *eptr && *eptr != 'e' && *eptr != 'E'; eptr++) ;

    if (*eptr == '\0') {
        /* no exponent — simple parse */
        float pv = strtof(s, &end);
        if (end == s || *end != '\0' || pv <= 0.0f) return -1;
        *out = -log10f(pv);
        return 0;
    }

    /* parse mantissa */
    mlen = (size_t)(eptr - s);
    if (mlen >= sizeof(buf)) return -1;
    memcpy(buf, s, mlen);
    buf[mlen] = '\0';
    mantissa = strtof(buf, &end);
    if (end == buf) return -1;
    /* parse exponent */
    exponent = strtof(eptr + 1, &end);
    if (*end != '\0') return -1;
    /* -log10(mantissa * 10^exponent) = -log10(mantissa) - exponent */
    if (mantissa <= 0.0f) return -1;
    *out = -log10f(mantissa) - exponent;
    return 0;
}

static char *score_dup(const char *s) {
    size_t n;
    char *out;
    if (!s) return NULL;
    n = strlen(s) + 1;
    out = (char *)duckdb_malloc(n);
    if (!out) return NULL;
    memcpy(out, s, n);
    return out;
}

static void score_strtoupper(char *s) {
    while (*s) {
        *s = (char)toupper((unsigned char)*s);
        s++;
    }
}

static score_mapping_t *score_mapping_preset(const char *preset, int *n) {
    if (!preset) return NULL;
    if (strcasecmp(preset, "PLINK") == 0) {
        *n = (int)(sizeof(score_plink_mapping) / sizeof(score_mapping_t));
        return score_plink_mapping;
    }
    if (strcasecmp(preset, "PLINK2") == 0) {
        *n = (int)(sizeof(score_plink2_mapping) / sizeof(score_mapping_t));
        return score_plink2_mapping;
    }
    if (strcasecmp(preset, "REGENIE") == 0) {
        *n = (int)(sizeof(score_regenie_mapping) / sizeof(score_mapping_t));
        return score_regenie_mapping;
    }
    if (strcasecmp(preset, "SAIGE") == 0) {
        *n = (int)(sizeof(score_saige_mapping) / sizeof(score_mapping_t));
        return score_saige_mapping;
    }
    if (strcasecmp(preset, "BOLT") == 0) {
        *n = (int)(sizeof(score_bolt_mapping) / sizeof(score_mapping_t));
        return score_bolt_mapping;
    }
    if (strcasecmp(preset, "METAL") == 0) {
        *n = (int)(sizeof(score_metal_mapping) / sizeof(score_mapping_t));
        return score_metal_mapping;
    }
    if (strcasecmp(preset, "PGS") == 0) {
        *n = (int)(sizeof(score_pgs_mapping) / sizeof(score_mapping_t));
        return score_pgs_mapping;
    }
    if (strcasecmp(preset, "SSF") == 0 || strcasecmp(preset, "GWAS-SSF") == 0) {
        *n = (int)(sizeof(score_ssf_mapping) / sizeof(score_mapping_t));
        return score_ssf_mapping;
    }
    return NULL;
}

static score_mapping_t *score_mapping_file(const char *path, int *n) {
    FILE *fp;
    char line[4096];
    score_mapping_t *mapping = NULL;
    int mapping_n = 0;
    int mapping_m = 0;
    if (!path) return NULL;
    fp = fopen(path, "r");
    if (!fp) return NULL;
    while (fgets(line, sizeof(line), fp)) {
        char *tab;
        int i;
        size_t len = strlen(line);
        if (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[len - 1] = '\0';
        tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = '\0';
        for (i = 0; i < SCORE_HDR_SIZE; i++) {
            if (strcmp(tab + 1, score_col_headers[i]) != 0) continue;
            if (mapping_n == mapping_m) {
                int new_m = mapping_m ? mapping_m * 2 : 16;
                score_mapping_t *new_map = (score_mapping_t *)realloc(mapping, sizeof(score_mapping_t) * (size_t)new_m);
                if (!new_map) {
                    fclose(fp);
                    duckdb_free(mapping);
                    return NULL;
                }
                mapping = new_map;
                mapping_m = new_m;
            }
            mapping[mapping_n].hdr_str = score_dup(line);
            mapping[mapping_n].hdr_num = i;
            mapping_n++;
        }
    }
    fclose(fp);
    *n = mapping_n;
    return mapping;
}

static void score_destroy_mapping(score_mapping_t *mapping, int n, int owns_strings) {
    int i;
    if (!mapping) return;
    if (owns_strings) {
        for (i = 0; i < n; i++) {
            if (mapping[i].hdr_str) duckdb_free(mapping[i].hdr_str);
        }
    }
    duckdb_free(mapping);
}

static int score_parse_q_thresholds(const char *s, double **out_arr, int *out_n) {
    char *tmp;
    char *tok;
    int n = 0;
    int m = 0;
    double *arr = NULL;
    if (!s || !*s) {
        *out_arr = NULL;
        *out_n = 1;
        return 0;
    }
    tmp = score_dup(s);
    if (!tmp) return -1;
    tok = strtok(tmp, ",");
    while (tok) {
        char *end = NULL;
        /* Upstream uses strtof() for threshold parsing (score.c:453).
         * The float→double promotion before -log10() reproduces the exact
         * same boundary-precision behavior as upstream bcftools +score.
         * For very small p-values (< ~1e-38) that underflow to 0 in float,
         * fall back to strtod to avoid silent breakage (upstream would
         * produce +inf threshold, excluding all markers). */
        float pf = strtof(tok, &end);
        double p;
        if (pf == 0.0f && end != tok) {
            /* strtof succeeded but returned 0 — could be underflow.
             * Re-parse with strtod to handle very small p-values. */
            p = strtod(tok, &end);
        } else {
            p = (double)pf;
        }
        while (*end && isspace((unsigned char)*end)) end++;
        if (end == tok || *end != '\0' || p <= 0.0 || p > 1.0) {
            duckdb_free(tmp);
            duckdb_free(arr);
            return -1;
        }
        if (n == m) {
            int new_m = m ? m * 2 : 8;
            double *new_arr = (double *)realloc(arr, sizeof(double) * (size_t)new_m);
            if (!new_arr) {
                duckdb_free(tmp);
                duckdb_free(arr);
                return -1;
            }
            arr = new_arr;
            m = new_m;
        }
        arr[n++] = -log10(p);
        tok = strtok(NULL, ",");
    }
    duckdb_free(tmp);
    if (!n) {
        duckdb_free(arr);
        *out_arr = NULL;
        *out_n = 1;
        return 0;
    }
    *out_arr = arr;
    *out_n = n;
    return 0;
}

static void score_destroy_bind(void *data) {
    score_bind_t *bind = (score_bind_t *)data;
    if (!bind) return;
    if (bind->bcf_path) duckdb_free(bind->bcf_path);
    if (bind->summary_path) duckdb_free(bind->summary_path);
    if (bind->columns_preset) duckdb_free(bind->columns_preset);
    if (bind->columns_file) duckdb_free(bind->columns_file);
    if (bind->q_thr_lp) duckdb_free(bind->q_thr_lp);
    if (bind->prs_name) duckdb_free(bind->prs_name);
    duckdb_free(bind);
}

static void score_destroy_init(void *data) {
    score_init_t *init = (score_init_t *)data;
    int i;
    if (!init) return;
    if (init->sample_names) {
        for (i = 0; i < init->n_samples; i++) {
            if (init->sample_names[i]) duckdb_free(init->sample_names[i]);
        }
        duckdb_free(init->sample_names);
    }
    if (init->metric_names) {
        for (i = 0; i < init->n_metrics; i++) {
            if (init->metric_names[i]) duckdb_free(init->metric_names[i]);
        }
        duckdb_free(init->metric_names);
    }
    if (init->metric_values) {
        for (i = 0; i < init->n_metrics; i++) {
            if (init->metric_values[i]) duckdb_free(init->metric_values[i]);
        }
        duckdb_free(init->metric_values);
    }
    if (init->metric_is_count) duckdb_free(init->metric_is_count);
    duckdb_free(init);
}

static void score_destroy_summary(score_summary_t *summary) {
    int i;
    if (!summary) return;
    for (i = 0; i < summary->n_markers; i++) {
        if (summary->markers && summary->markers[i].a1) duckdb_free(summary->markers[i].a1);
    }
    if (summary->rid_pos2idx) kh_destroy(score64, (khash_t(score64) *)summary->rid_pos2idx);
    if (summary->id2idx) khash_str2int_destroy_free(summary->id2idx);
    if (summary->markers) duckdb_free(summary->markers);
    duckdb_free(summary);
}

static char score_detect_delimiter(const char *line) {
    return strchr(line, '\t') ? '\t' : '\0';
}

static int score_split_fields(char *line, char delimiter, char **fields, int max_fields) {
    int n = 0;
    char *p = line;
    while (*p && n < max_fields) {
        while (*p && ((delimiter == '\0') ? isspace((unsigned char)*p) : (*p == delimiter))) p++;
        if (!*p) break;
        fields[n++] = p;
        if (delimiter == '\0') {
            while (*p && !isspace((unsigned char)*p)) p++;
        } else {
            while (*p && *p != delimiter && *p != '\n' && *p != '\r') p++;
        }
        if (*p) {
            *p = '\0';
            p++;
        }
    }
    return n;
}

static int score_header_index(char **headers, int n_headers, const char *key) {
    int i;
    for (i = 0; i < n_headers; i++) {
        if (headers[i] && strcasecmp(headers[i], key) == 0) return i;
    }
    return -1;
}

static score_summary_t *score_load_summary(const score_bind_t *bind, bcf_hdr_t *bcf_hdr, char *err, size_t err_sz) {
    FILE *fp = NULL;
    char line[32768];
    char *fields[512];
    char delimiter = '\t';
    score_mapping_t *mapping = NULL;
    int mapping_n = 0;
    int mapping_owns = 0;
    int idx_snp = -1, idx_chr = -1, idx_bp = -1, idx_a1 = -1, idx_beta = -1, idx_or = -1, idx_p = -1, idx_lp = -1;
    score_summary_t *summary = NULL;
    void *unused_alleles = NULL;

    if (bind->columns_file && bind->columns_file[0]) {
        mapping = score_mapping_file(bind->columns_file, &mapping_n);
        mapping_owns = 1;
        if (!mapping) {
            snprintf(err, err_sz, "bcftools_score: failed reading columns_file '%s'", bind->columns_file);
            return NULL;
        }
    } else {
        mapping = score_mapping_preset(bind->columns_preset ? bind->columns_preset : "PLINK", &mapping_n);
        if (!mapping) {
            snprintf(err, err_sz, "bcftools_score: unknown preset '%s'", bind->columns_preset ? bind->columns_preset : "");
            return NULL;
        }
    }

    fp = fopen(bind->summary_path, "r");
    if (!fp) {
        snprintf(err, err_sz, "bcftools_score: cannot open summary '%s': %s", bind->summary_path, strerror(errno));
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' && strncmp(line, "#CHR", 4) != 0 && strncmp(line, "#ID", 3) != 0) continue;
        break;
    }
    if (line[0] == '\0') {
        snprintf(err, err_sz, "bcftools_score: empty summary file '%s'", bind->summary_path);
        fclose(fp);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }
    if (line[0] == '#') memmove(line, line + 1, strlen(line));
    delimiter = score_detect_delimiter(line);
    {
        int i;
        int n_fields = score_split_fields(line, delimiter, fields, 512);
        for (i = 0; i < mapping_n; i++) {
            int idx = score_header_index(fields, n_fields, mapping[i].hdr_str);
            if (idx < 0) continue;
            switch (mapping[i].hdr_num) {
            case SCORE_HDR_SNP: idx_snp = idx; break;
            case SCORE_HDR_CHR: idx_chr = idx; break;
            case SCORE_HDR_BP: idx_bp = idx; break;
            case SCORE_HDR_A1: idx_a1 = idx; break;
            case SCORE_HDR_BETA: idx_beta = idx; break;
            case SCORE_HDR_OR: idx_or = idx; break;
            case SCORE_HDR_P: idx_p = idx; break;
            case SCORE_HDR_LP: idx_lp = idx; break;
            default: break;
            }
        }
    }

    summary = (score_summary_t *)duckdb_malloc(sizeof(score_summary_t));
    memset(summary, 0, sizeof(*summary));
    summary->use_snp = bind->use_variant_id || idx_chr < 0 || idx_bp < 0;
    summary->snp_ok = idx_snp >= 0;
    summary->chr_ok = idx_chr >= 0;
    summary->bp_ok = idx_bp >= 0;
    summary->a1_ok = idx_a1 >= 0;
    summary->beta_ok = idx_beta >= 0 || idx_or >= 0;
    summary->p_ok = idx_p >= 0 || idx_lp >= 0;
    if (!summary->a1_ok || !summary->beta_ok) {
        snprintf(err, err_sz, "bcftools_score: summary must include effect allele and effect size/OR columns");
        score_destroy_summary(summary);
        fclose(fp);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }
    if (bind->use_variant_id && !summary->snp_ok) {
        snprintf(err, err_sz, "bcftools_score: summary must include marker name column when use_variant_id=true");
        score_destroy_summary(summary);
        fclose(fp);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }
    if (summary->use_snp && !summary->snp_ok) {
        snprintf(err, err_sz, "bcftools_score: summary must include marker name column when matching by variant ID");
        score_destroy_summary(summary);
        fclose(fp);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }

    if (summary->use_snp) {
        summary->id2idx = khash_str2int_init();
    } else {
        summary->rid_pos2idx = kh_init(score64);
    }

    (void)unused_alleles;

    while (fgets(line, sizeof(line), fp)) {
        int n_fields;
        float es = NAN;
        float lp = NAN;
        score_marker_t marker;
        char *id = NULL;
        int rid = -1;
        int64_t pos = -1;

        marker.a1 = NULL;
        if (line[0] == '#') continue;
        n_fields = score_split_fields(line, delimiter, fields, 512);
        if (n_fields == 0) continue;

        if (idx_a1 < n_fields) {
            marker.a1 = score_dup(fields[idx_a1]);
            if (!marker.a1) continue;
            score_strtoupper(marker.a1);
        } else {
            continue;
        }

        if (idx_beta >= 0 && idx_beta < n_fields) {
            if (score_parse_float(fields[idx_beta], &es) != 0) es = NAN;
        }
        if (isnan(es) && idx_or >= 0 && idx_or < n_fields) {
            float orv = NAN;
            if (score_parse_float(fields[idx_or], &orv) == 0 && !isnan(orv) && orv > 0.0f)
                es = logf(orv);
        }
        if (isnan(es)) {
            free(marker.a1);
            continue;
        }

        if (idx_lp >= 0 && idx_lp < n_fields) {
            if (score_parse_float(fields[idx_lp], &lp) != 0) lp = NAN;
        }
        if (isnan(lp) && idx_p >= 0 && idx_p < n_fields) {
            if (score_parse_pvalue_lp(fields[idx_p], &lp) != 0) lp = NAN;
        }

        marker.es = es;
        marker.lp = lp;

        if (summary->n_markers == summary->m_markers) {
            int new_m = summary->m_markers ? summary->m_markers * 2 : 1024;
            score_marker_t *new_markers = (score_marker_t *)realloc(summary->markers, sizeof(score_marker_t) * (size_t)new_m);
            if (!new_markers) continue;
            summary->markers = new_markers;
            summary->m_markers = new_m;
        }

        if (summary->use_snp) {
            int size;
            int ret;
            if (idx_snp < 0 || idx_snp >= n_fields || fields[idx_snp][0] == '\0') continue;
            id = score_dup(fields[idx_snp]);
            if (!id) continue;
            size = khash_str2int_size(summary->id2idx);
            ret = khash_str2int_inc(summary->id2idx, id);
            if (ret == size) {
                summary->markers[summary->n_markers++] = marker;
            } else {
                duckdb_free(id);
                if (marker.a1) duckdb_free(marker.a1);
            }
        } else {
            char *end = NULL;
            khash_t(score64) *hash;
            khiter_t it;
            int ret;
            uint64_t key;
            if (idx_chr < 0 || idx_chr >= n_fields || idx_bp < 0 || idx_bp >= n_fields) continue;
            rid = score_hdr_name2id_flexible(bcf_hdr, fields[idx_chr]);
            if (rid < 0) continue;
            pos = strtoll(fields[idx_bp], &end, 10) - 1;
            if (end == fields[idx_bp] || *end != '\0' || pos < 0) continue;
            key = (((uint64_t)(uint32_t)rid) << 44) | (uint64_t)pos;
            hash = (khash_t(score64) *)summary->rid_pos2idx;
            it = kh_put(score64, hash, key, &ret);
            if (ret > 0) {
                kh_val(hash, it) = summary->n_markers;
                summary->markers[summary->n_markers++] = marker;
            } else {
                if (marker.a1) duckdb_free(marker.a1);
            }
        }
    }

    fclose(fp);
    if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
    return summary;
}

static int score_resolve_use_tag(score_bind_t *bind, bcf_hdr_t *hdr, char *err, size_t err_sz) {
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    int ds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS");
    int hds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HDS");
    int ap1_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AP1");
    int ap2_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AP2");
    int gp_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GP");
    int as_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AS");

    if (!bind->use_tag) {
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id)) bind->use_tag = SCORE_USE_GT;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, hds_id)) bind->use_tag = SCORE_USE_HDS;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap1_id) && bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap2_id)) bind->use_tag = SCORE_USE_AP;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gp_id)) bind->use_tag = SCORE_USE_GP;
        if (bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ds_id)) bind->use_tag = SCORE_USE_DS;
    }

    if (!bind->use_tag) {
        snprintf(err, err_sz, "bcftools_score: no supported FORMAT tag found (GT/DS/HDS/AP1+AP2/GP/AS)");
        return -1;
    }
    if (bind->use_tag == SCORE_USE_AS && !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, as_id)) {
        snprintf(err, err_sz, "bcftools_score: AS FORMAT not found");
        return -1;
    }
    return 0;
}

static void score_bind(duckdb_bind_info info) {
    duckdb_value bcf_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value summary_val = duckdb_bind_get_parameter(info, 1);
    duckdb_value use_val;
    duckdb_value columns_val;
    duckdb_value columns_file_val;
    duckdb_value variant_id_val;
    duckdb_value q_thr_val;
    duckdb_value counts_val;
    duckdb_logical_type varchar_type;
    duckdb_logical_type bool_type;
    duckdb_logical_type dbl_type;
    duckdb_logical_type bigint_type;
    score_bind_t *bind = (score_bind_t *)duckdb_malloc(sizeof(score_bind_t));

    memset(bind, 0, sizeof(*bind));
    bind->bcf_path = duckdb_get_varchar(bcf_val);
    bind->summary_path = duckdb_get_varchar(summary_val);
    bind->n_q_thr = 1;
    bind->columns_preset = score_dup("PLINK");

    duckdb_destroy_value(&bcf_val);
    duckdb_destroy_value(&summary_val);

    if (!bind->bcf_path || !bind->summary_path || bind->bcf_path[0] == '\0' || bind->summary_path[0] == '\0') {
        duckdb_bind_set_error(info, "bcftools_score: bcf_path and summary_path are required");
        score_destroy_bind(bind);
        return;
    }

    use_val = duckdb_bind_get_named_parameter(info, "use");
    if (use_val && !duckdb_is_null_value(use_val)) {
        char *use_str = duckdb_get_varchar(use_val);
        if (use_str) {
            if (strcasecmp(use_str, "GT") == 0) bind->use_tag = SCORE_USE_GT;
            else if (strcasecmp(use_str, "DS") == 0) bind->use_tag = SCORE_USE_DS;
            else if (strcasecmp(use_str, "HDS") == 0) bind->use_tag = SCORE_USE_HDS;
            else if (strcasecmp(use_str, "AP") == 0) bind->use_tag = SCORE_USE_AP;
            else if (strcasecmp(use_str, "GP") == 0) bind->use_tag = SCORE_USE_GP;
            else if (strcasecmp(use_str, "AS") == 0) bind->use_tag = SCORE_USE_AS;
            else {
                duckdb_bind_set_error(info, "bcftools_score: use must be one of GT,DS,HDS,AP,GP,AS");
                duckdb_free(use_str);
                duckdb_destroy_value(&use_val);
                score_destroy_bind(bind);
                return;
            }
            duckdb_free(use_str);
        }
    }
    if (use_val) duckdb_destroy_value(&use_val);

    columns_val = duckdb_bind_get_named_parameter(info, "columns");
    if (columns_val && !duckdb_is_null_value(columns_val)) {
        char *p = duckdb_get_varchar(columns_val);
        if (p) {
            if (bind->columns_preset) duckdb_free(bind->columns_preset);
            bind->columns_preset = p;
        }
    }
    if (columns_val) duckdb_destroy_value(&columns_val);

    columns_file_val = duckdb_bind_get_named_parameter(info, "columns_file");
    if (columns_file_val && !duckdb_is_null_value(columns_file_val)) {
        bind->columns_file = duckdb_get_varchar(columns_file_val);
    }
    if (columns_file_val) duckdb_destroy_value(&columns_file_val);

    variant_id_val = duckdb_bind_get_named_parameter(info, "use_variant_id");
    if (variant_id_val && !duckdb_is_null_value(variant_id_val)) bind->use_variant_id = duckdb_get_bool(variant_id_val) ? 1 : 0;
    if (variant_id_val) duckdb_destroy_value(&variant_id_val);

    q_thr_val = duckdb_bind_get_named_parameter(info, "q_score_thr");
    if (q_thr_val && !duckdb_is_null_value(q_thr_val)) {
        char *q_str = duckdb_get_varchar(q_thr_val);
        if (score_parse_q_thresholds(q_str, &bind->q_thr_lp, &bind->n_q_thr) != 0) {
            duckdb_bind_set_error(info, "bcftools_score: q_score_thr must be a comma-separated list of p-values in (0,1]");
            duckdb_free(q_str);
            duckdb_destroy_value(&q_thr_val);
            score_destroy_bind(bind);
            return;
        }
        if (q_str) duckdb_free(q_str);
    }
    if (q_thr_val) duckdb_destroy_value(&q_thr_val);

    counts_val = duckdb_bind_get_named_parameter(info, "counts");
    if (counts_val && !duckdb_is_null_value(counts_val)) bind->counts = duckdb_get_bool(counts_val) ? 1 : 0;
    if (counts_val) duckdb_destroy_value(&counts_val);

    {
        const char *slash = strrchr(bind->summary_path, '/');
        const char *base = slash ? slash + 1 : bind->summary_path;
        const char *dot = strrchr(base, '.');
        size_t len = dot && dot > base ? (size_t)(dot - base) : strlen(base);
        bind->prs_name = (char *)duckdb_malloc(len + 1);
        if (!bind->prs_name) {
            duckdb_bind_set_error(info, "bcftools_score: out of memory");
            score_destroy_bind(bind);
            return;
        }
        memcpy(bind->prs_name, base, len);
        bind->prs_name[len] = '\0';
    }

    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    dbl_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "SAMPLE", varchar_type);
    {
        int j;
        for (j = 0; j < bind->n_q_thr; j++) {
            char name[256];
            if (bind->q_thr_lp && bind->n_q_thr > 0) {
                double p = pow(10.0, -bind->q_thr_lp[j]);
                snprintf(name, sizeof(name), "%s_p%.6g", bind->prs_name, p);
            } else {
                snprintf(name, sizeof(name), "%s", bind->prs_name);
            }
            duckdb_bind_add_result_column(info, name, dbl_type);
            if (bind->counts) {
                /* Upstream builds CNT name as "<prs>_CNT" then appends "_p<thr>" */
                char cnt_name[256];
                if (bind->q_thr_lp && bind->n_q_thr > 0) {
                    double p = pow(10.0, -bind->q_thr_lp[j]);
                    snprintf(cnt_name, sizeof(cnt_name), "%s_CNT_p%.6g", bind->prs_name, p);
                } else {
                    snprintf(cnt_name, sizeof(cnt_name), "%s_CNT", bind->prs_name);
                }
                duckdb_bind_add_result_column(info, cnt_name, bigint_type);
            }
        }
    }
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&dbl_type);
    duckdb_destroy_logical_type(&bigint_type);

    duckdb_bind_set_bind_data(info, bind, score_destroy_bind);
}

static void score_init(duckdb_init_info info) {
    score_bind_t *bind = (score_bind_t *)duckdb_init_get_bind_data(info);
    score_init_t *init = (score_init_t *)duckdb_malloc(sizeof(score_init_t));
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    score_summary_t *summary = NULL;
    void *unused_alleles = NULL;
    float *aps = NULL;
    int m_aps = 0;
    int32_t *int32_arr = NULL;
    int m_int32 = 0;
    float *float_arr = NULL;
    int m_float = 0;
    float *missing = NULL;
    int n_samples;
    int n_metrics;
    int i;
    char err[512];

    memset(init, 0, sizeof(*init));
    fp = hts_open(bind->bcf_path, "r");
    if (!fp) {
        duckdb_init_set_error(info, "bcftools_score: failed to open BCF/VCF input");
        duckdb_free(init);
        return;
    }
    hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        duckdb_init_set_error(info, "bcftools_score: failed to read BCF/VCF header");
        duckdb_free(init);
        return;
    }

    if (score_resolve_use_tag(bind, hdr, err, sizeof(err)) != 0) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        duckdb_init_set_error(info, err);
        duckdb_free(init);
        return;
    }

    summary = score_load_summary(bind, hdr, err, sizeof(err));
    if (!summary) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        duckdb_init_set_error(info, err[0] ? err : "bcftools_score: failed to load summary file");
        duckdb_free(init);
        return;
    }

    n_samples = bcf_hdr_nsamples(hdr);
    n_metrics = bind->n_q_thr * (bind->counts ? 2 : 1);
    init->n_samples = n_samples;
    init->n_metrics = n_metrics;
    init->sample_names = (char **)duckdb_malloc(sizeof(char *) * (size_t)n_samples);
    init->metric_names = (char **)duckdb_malloc(sizeof(char *) * (size_t)n_metrics);
    init->metric_values = (double **)duckdb_malloc(sizeof(double *) * (size_t)n_metrics);
    init->metric_is_count = (uint8_t *)duckdb_malloc((size_t)n_metrics);
    if (!init->sample_names || !init->metric_names || !init->metric_values || !init->metric_is_count) {
        score_destroy_summary(summary);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        score_destroy_init(init);
        duckdb_init_set_error(info, "bcftools_score: out of memory");
        return;
    }

    for (i = 0; i < n_samples; i++) {
        init->sample_names[i] = score_dup(hdr->samples[i]);
    }
    {
        int j;
        int metric_idx = 0;
        memset(init->metric_is_count, 0, (size_t)n_metrics);
        for (j = 0; j < bind->n_q_thr; j++) {
            char name[256];
            if (bind->q_thr_lp && bind->n_q_thr > 0) {
                double p = pow(10.0, -bind->q_thr_lp[j]);
                snprintf(name, sizeof(name), "%s_p%.6g", bind->prs_name, p);
            } else {
                snprintf(name, sizeof(name), "%s", bind->prs_name);
            }
            init->metric_names[metric_idx] = score_dup(name);
            init->metric_values[metric_idx] = (double *)duckdb_malloc(sizeof(double) * (size_t)n_samples);
            memset(init->metric_values[metric_idx], 0, sizeof(double) * (size_t)n_samples);
            init->metric_is_count[metric_idx] = 0;
            metric_idx++;
            if (bind->counts) {
                /* Upstream builds CNT name as "<prs>_CNT" then appends "_p<thr>" */
                char cnt_name[256];
                if (bind->q_thr_lp && bind->n_q_thr > 0) {
                    double p = pow(10.0, -bind->q_thr_lp[j]);
                    snprintf(cnt_name, sizeof(cnt_name), "%s_CNT_p%.6g", bind->prs_name, p);
                } else {
                    snprintf(cnt_name, sizeof(cnt_name), "%s_CNT", bind->prs_name);
                }
                init->metric_names[metric_idx] = score_dup(cnt_name);
                init->metric_values[metric_idx] = (double *)duckdb_malloc(sizeof(double) * (size_t)n_samples);
                memset(init->metric_values[metric_idx], 0, sizeof(double) * (size_t)n_samples);
                init->metric_is_count[metric_idx] = 1;
                metric_idx++;
            }
        }
    }

    (void)unused_alleles;
    missing = (float *)duckdb_malloc(sizeof(float) * (size_t)n_samples);
    rec = bcf_init();

    while (bcf_read(fp, hdr, rec) == 0) {
        int marker_idx = -1;
        int idx_allele;
        float es;
        float lp;
        int number;
        int k;
        bcf_unpack(rec, BCF_UN_ALL);
        if (summary->use_snp) {
            if (rec->d.id && rec->d.id[0] != '\0') {
                khash_str2int_get(summary->id2idx, rec->d.id, &marker_idx);
            }
        } else {
            khash_t(score64) *hash = (khash_t(score64) *)summary->rid_pos2idx;
            uint64_t key = (((uint64_t)(uint32_t)rec->rid) << 44) | (uint64_t)rec->pos;
            khiter_t it = kh_get(score64, hash, key);
            if (it != kh_end(hash)) marker_idx = kh_val(hash, it);
        }
        if (marker_idx < 0 || marker_idx >= summary->n_markers) continue;

        memset(missing, 0, sizeof(float) * (size_t)n_samples);
        if (m_aps < rec->n_allele * n_samples) {
            int new_m = rec->n_allele * n_samples;
            float *new_aps = (float *)realloc(aps, sizeof(float) * (size_t)new_m);
            if (!new_aps) continue;
            aps = new_aps;
            m_aps = new_m;
        }
        memset(aps, 0, sizeof(float) * (size_t)(rec->n_allele * n_samples));

        switch (bind->use_tag) {
        case SCORE_USE_GT:
            number = bcf_get_genotypes(hdr, rec, &int32_arr, &m_int32);
            if (number <= 0) continue;
            number /= n_samples;
            if (number < 2) continue;
            for (k = 0; k < n_samples; k++) {
                int32_t *ptr = int32_arr + number * k;
                if (bcf_gt_is_missing(ptr[0])) {
                    missing[k] = 1;
                } else {
                    int a0 = bcf_gt_allele(ptr[0]);
                    if (a0 >= 0 && a0 < rec->n_allele) aps[a0 * n_samples + k] += 1.0f;
                    if (ptr[1] == bcf_int32_vector_end) {
                        /* haploid (e.g. chrX male) — one allele only */
                    } else if (bcf_gt_is_missing(ptr[1])) {
                        missing[k] = 1;
                    } else {
                        int a1 = bcf_gt_allele(ptr[1]);
                        if (a1 >= 0 && a1 < rec->n_allele) aps[a1 * n_samples + k] += 1.0f;
                    }
                }
            }
            break;
        case SCORE_USE_DS:
            number = bcf_get_format_float(hdr, rec, "DS", &float_arr, &m_float);
            if (number <= 0) continue;
            number /= n_samples;
            if (number == 1 && rec->n_allele == 2) {
                for (k = 0; k < n_samples; k++) {
                    if (score_is_missing(float_arr[k])) missing[k] = 1;
                    else {
                        aps[k] += 2.0f - float_arr[k];
                        aps[n_samples + k] += float_arr[k];
                    }
                }
            } else {
                for (k = 0; k < n_samples; k++) {
                    int idx;
                    float *ptr = float_arr + number * k;
                    aps[k] += 2.0f;
                    for (idx = 0; idx < number; idx++) {
                        if (score_is_missing(ptr[idx])) missing[k] = 1;
                        else {
                            aps[k] -= ptr[idx];
                            aps[(idx + 1) * n_samples + k] += ptr[idx];
                        }
                    }
                }
            }
            break;
        case SCORE_USE_HDS:
            number = bcf_get_format_float(hdr, rec, "HDS", &float_arr, &m_float);
            if (number <= 0) continue;
            number /= n_samples;
            if (number != 2 || rec->n_allele != 2) continue;
            for (k = 0; k < n_samples; k++) {
                if (score_is_missing(float_arr[2 * k]) || score_is_missing(float_arr[2 * k + 1])) missing[k] = 1;
                else {
                    aps[k] += 2.0f - float_arr[2 * k] - float_arr[2 * k + 1];
                    aps[n_samples + k] += float_arr[2 * k] + float_arr[2 * k + 1];
                }
            }
            break;
        case SCORE_USE_AP: {
            const char *ap_tags[2] = {"AP1", "AP2"};
            int ap;
            for (ap = 0; ap < 2; ap++) {
                number = bcf_get_format_float(hdr, rec, ap_tags[ap], &float_arr, &m_float);
                if (number <= 0) continue;
                number /= n_samples;
                if (number == 1 && rec->n_allele == 2) {
                    for (k = 0; k < n_samples; k++) {
                        if (score_is_missing(float_arr[k])) missing[k] = 1;
                        else {
                            aps[k] += 1.0f - float_arr[k];
                            aps[n_samples + k] += float_arr[k];
                        }
                    }
                } else {
                    for (k = 0; k < n_samples; k++) {
                        int idx;
                        float *ptr = float_arr + number * k;
                        aps[k] += 1.0f;
                        for (idx = 0; idx < number; idx++) {
                            if (score_is_missing(ptr[idx])) missing[k] = 1;
                            else {
                                aps[k] -= ptr[idx];
                                aps[(idx + 1) * n_samples + k] += ptr[idx];
                            }
                        }
                    }
                }
            }
            break;
        }
        case SCORE_USE_GP:
            number = bcf_get_format_float(hdr, rec, "GP", &float_arr, &m_float);
            if (number <= 0) continue;
            number /= n_samples;
            if (number == 3 && rec->n_allele == 2) {
                for (k = 0; k < n_samples; k++) {
                    float *ptr = float_arr + number * k;
                    if (score_is_missing(ptr[0]) || score_is_missing(ptr[1]) || score_is_missing(ptr[2])) missing[k] = 1;
                    else {
                        aps[k] += 2.0f * ptr[0] + ptr[1];
                        aps[n_samples + k] += ptr[1] + 2.0f * ptr[2];
                    }
                }
            } else {
                for (k = 0; k < n_samples; k++) {
                    int a, b;
                    float *ptr = float_arr + number * k;
                    for (b = 0; b < rec->n_allele; b++) {
                        for (a = 0; a <= b; a++) {
                            int idx = b * (b + 1) / 2 + a;
                            if (score_is_missing(ptr[idx])) missing[k] = 1;
                            else {
                                aps[a * n_samples + k] += ptr[idx];
                                aps[b * n_samples + k] += ptr[idx];
                            }
                        }
                    }
                }
            }
            break;
        case SCORE_USE_AS:
            number = bcf_get_format_int32(hdr, rec, "AS", &int32_arr, &m_int32);
            if (number <= 0) continue;
            number /= n_samples;
            if (number != 1 || rec->n_allele != 2) continue;
            for (k = 0; k < n_samples; k++) {
                if (int32_arr[k] == 0 || int32_arr[k] == bcf_int32_missing) missing[k] = 1;
                else {
                    aps[k] -= (float)int32_arr[k];
                    aps[n_samples + k] += (float)int32_arr[k];
                }
            }
            break;
        default:
            continue;
        }

        {
            const score_marker_t *marker = &summary->markers[marker_idx];
            const char *a1 = marker->a1;
            int j;
            if (!a1) continue;
            for (idx_allele = 0; idx_allele < rec->n_allele; idx_allele++) {
                if (strcmp(a1, rec->d.allele[idx_allele]) == 0) break;
            }
            if (idx_allele == rec->n_allele) continue;

            es = marker->es;
            lp = marker->lp;
            for (j = 0; j < bind->n_q_thr; j++) {
                int out_idx = j * (bind->counts ? 2 : 1);
                if (bind->q_thr_lp && !isnan(lp) && lp < bind->q_thr_lp[j]) continue;
                for (k = 0; k < n_samples; k++) {
                    if (missing[k]) continue;
                    init->metric_values[out_idx][k] += (double)es * (double)aps[idx_allele * n_samples + k];
                    if (bind->counts) init->metric_values[out_idx + 1][k] += 1.0;
                }
            }
        }
    }

    if (int32_arr) free(int32_arr);
    if (float_arr) free(float_arr);
    if (aps) duckdb_free(aps);
    if (missing) duckdb_free(missing);
    if (rec) bcf_destroy(rec);
    if (summary) score_destroy_summary(summary);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);

    duckdb_init_set_init_data(info, init, score_destroy_init);
}

static void score_scan(duckdb_function_info info, duckdb_data_chunk output) {
    score_init_t *init = (score_init_t *)duckdb_function_get_init_data(info);
    idx_t chunk_size = duckdb_vector_size();
    idx_t remaining = (idx_t)(init->n_samples - (int)init->row_offset);
    idx_t n = remaining < chunk_size ? remaining : chunk_size;
    idx_t i;

    if (init->row_offset >= init->n_samples) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    for (i = 0; i < n; i++) {
        idx_t row = init->row_offset + i;
        duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 0), i, init->sample_names[row]);
    }

    {
        int m;
        for (m = 0; m < init->n_metrics; m++) {
            duckdb_vector vec = duckdb_data_chunk_get_vector(output, (idx_t)(m + 1));
            if (init->metric_is_count && init->metric_is_count[m]) {
                int64_t *out_i64 = (int64_t *)duckdb_vector_get_data(vec);
                for (i = 0; i < n; i++) {
                    idx_t row = init->row_offset + i;
                    out_i64[i] = (int64_t)(init->metric_values[m][row]);
                }
            } else {
                double *out = (double *)duckdb_vector_get_data(vec);
                for (i = 0; i < n; i++) {
                    idx_t row = init->row_offset + i;
                    out[i] = init->metric_values[m][row];
                }
            }
        }
    }

    init->row_offset += (int64_t)n;
    duckdb_data_chunk_set_size(output, n);
}

void register_bcftools_score_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_set_name(tf, "bcftools_score");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "use", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "columns", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "columns_file", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "q_score_thr", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "use_variant_id", bool_type);
    duckdb_table_function_add_named_parameter(tf, "counts", bool_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);

    duckdb_table_function_set_bind(tf, score_bind);
    duckdb_table_function_set_init(tf, score_init);
    duckdb_table_function_set_function(tf, score_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
