#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <htslib/hts.h>
#include <htslib/khash_str2int.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include "bcftools_shim.h"
#include "filter.h"
#include "duckhts_column_names.h"
#include "duckhts_bind_params.h"

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
/* Sentinel for columns that upstream recognizes but does not extract data from
 * (NULL setter in upstream score.c:141-164).  These entries ensure the header
 * matcher recognizes the column name and silently consumes it, preventing
 * misidentification of adjacent columns.  Falls through to default:break in
 * the switch at header-match time (S-M5). */
#define SCORE_HDR_UNUSED -1

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
    char **summary_paths;
    int n_summary_paths;
    int m_summary_paths;
    int *summary_prs_counts;
    char *summaries_list_file;
    char *log_path;
    char *columns_preset;
    char *columns_file;
    char *samples;
    char *regions;
    char *regions_file;
    int regions_overlap;
    char *targets;
    char *targets_file;
    int targets_overlap;
    char *apply_filters;
    char *include_expr;
    char *exclude_expr;
    int force_samples;
    double *q_thr_lp;
    int n_q_thr;
    int gwas_vcf_mode;    /* 1 if summaries are GWAS-VCF (FORMAT/ES+LP), 0 for TSV summaries */
    int n_prs;            /* number of PRS output tracks */
    char **prs_names;     /* PRS output names (TSV basenames or GWAS-VCF sample names) */
    char **prs_paths;     /* source summary path for each PRS output track */
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
    int64_t all_markers;
    int64_t duplicate_markers;
} score_summary_t;

typedef struct {
    int rid;
    hts_pos_t beg;
    hts_pos_t end;
} score_interval_t;

typedef struct {
    score_interval_t *intervals;
    size_t count;
    int active;
} score_region_list_t;

typedef struct {
    bcf_srs_t *sr;
    score_region_list_t regions;
    score_region_list_t targets;
    filter_t *include_filt;
    filter_t *exclude_filt;
} score_scan_source_t;

static char *score_to_filter_expr(const char *expr) {
    size_t i;
    size_t j;
    size_t len;
    char *out;
    if (!expr) return NULL;
    len = strlen(expr);
    out = (char *)malloc(len + 1);
    if (!out) return NULL;
    for (i = 0, j = 0; i < len;) {
        if (expr[i] == '>' || expr[i] == '<') {
            char op = expr[i];
            i++;
            if (i < len && expr[i] == op) {
                free(out);
                return strdup("");
            }
            out[j++] = op;
            if (i < len && expr[i] == '=') out[j++] = expr[i++];
            continue;
        }
        out[j++] = expr[i++];
    }
    out[j] = '\0';
    return out;
}

static const char *score_col_headers[SCORE_HDR_SIZE] = {
    "SNP", "BP", "CHR", "A1", "A2", "P", "Z", "OR", "BETA", "N", "N_CAS", "N_CON",
    "INFO", "FRQ", "A*", "SE", "LP", "AC", "NEFF", "NEFFDIV2", "HET_I2", "HET_P", "HET_LP", "DIRE"
};

static score_mapping_t score_plink_mapping[] = {
    {"SNP", SCORE_HDR_SNP}, {"BP", SCORE_HDR_BP}, {"CHR", SCORE_HDR_CHR}, {"A1", SCORE_HDR_A1},
    {"A2", SCORE_HDR_UNUSED}, {"P", SCORE_HDR_P}, {"OR", SCORE_HDR_OR}, {"BETA", SCORE_HDR_BETA},
    {"INFO", SCORE_HDR_UNUSED}, {"F_U", SCORE_HDR_UNUSED}, {"FRQ", SCORE_HDR_UNUSED}, {"SE", SCORE_HDR_UNUSED}
};

static score_mapping_t score_plink2_mapping[] = {
    {"ID", SCORE_HDR_SNP}, {"POS", SCORE_HDR_BP}, {"CHROM", SCORE_HDR_CHR}, {"A1", SCORE_HDR_A1},
    {"AX", SCORE_HDR_UNUSED}, {"P", SCORE_HDR_P}, {"Z_STAT", SCORE_HDR_UNUSED}, {"OR", SCORE_HDR_OR},
    {"BETA", SCORE_HDR_BETA}, {"MACH_R2", SCORE_HDR_UNUSED}, {"A1_FREQ", SCORE_HDR_UNUSED},
    {"SE", SCORE_HDR_UNUSED}, {"LOG(OR)_SE", SCORE_HDR_UNUSED}, {"LOG10_P", SCORE_HDR_LP}
};

static score_mapping_t score_regenie_mapping[] = {
    {"ID", SCORE_HDR_SNP}, {"GENPOS", SCORE_HDR_BP}, {"CHROM", SCORE_HDR_CHR}, {"ALLELE1", SCORE_HDR_A1},
    {"BETA", SCORE_HDR_BETA}, {"N", SCORE_HDR_UNUSED}, {"N_CASES", SCORE_HDR_UNUSED},
    {"N_CONTROLS", SCORE_HDR_UNUSED}, {"INFO", SCORE_HDR_UNUSED}, {"A1FREQ", SCORE_HDR_UNUSED},
    {"ALLELE0", SCORE_HDR_UNUSED}, {"SE", SCORE_HDR_UNUSED}, {"LOG10P", SCORE_HDR_LP},
    {"AC_ALLELE1", SCORE_HDR_UNUSED}
};

static score_mapping_t score_saige_mapping[] = {
    {"SNPID", SCORE_HDR_SNP}, {"markerID", SCORE_HDR_SNP}, {"POS", SCORE_HDR_BP}, {"CHR", SCORE_HDR_CHR},
    {"Allele2", SCORE_HDR_A1}, {"Allele1", SCORE_HDR_UNUSED}, {"p.value", SCORE_HDR_P}, {"BETA", SCORE_HDR_BETA},
    {"N", SCORE_HDR_UNUSED}, {"N_case", SCORE_HDR_UNUSED}, {"N_ctrl", SCORE_HDR_UNUSED},
    {"imputationInfo", SCORE_HDR_UNUSED}, {"AF_Allele2", SCORE_HDR_UNUSED}, {"SE", SCORE_HDR_UNUSED},
    {"AC_Allele2", SCORE_HDR_UNUSED}
};

static score_mapping_t score_bolt_mapping[] = {
    {"SNP", SCORE_HDR_SNP}, {"BP", SCORE_HDR_BP}, {"CHR", SCORE_HDR_CHR}, {"ALLELE1", SCORE_HDR_A1},
    {"P_LINREG", SCORE_HDR_P}, {"P_BOLT_LMM_INF", SCORE_HDR_P}, {"P_BOLT_LMM", SCORE_HDR_P},
    {"BETA", SCORE_HDR_BETA}, {"A1FREQ", SCORE_HDR_UNUSED}, {"ALLELE0", SCORE_HDR_UNUSED},
    {"SE", SCORE_HDR_UNUSED}
};

static score_mapping_t score_metal_mapping[] = {
    {"MarkerName", SCORE_HDR_SNP}, {"Position", SCORE_HDR_BP}, {"Chromosome", SCORE_HDR_CHR},
    {"Allele1", SCORE_HDR_A1}, {"Allele2", SCORE_HDR_UNUSED}, {"P-value", SCORE_HDR_P},
    {"Zscore", SCORE_HDR_UNUSED}, {"Effect", SCORE_HDR_BETA}, {"N", SCORE_HDR_UNUSED},
    {"Freq1", SCORE_HDR_UNUSED}, {"StdErr", SCORE_HDR_UNUSED}, {"log(P)", SCORE_HDR_LP},
    {"Weight", SCORE_HDR_UNUSED}, {"HetISq", SCORE_HDR_UNUSED}, {"HetPVal", SCORE_HDR_UNUSED},
    {"logHetP", SCORE_HDR_UNUSED}, {"Direction", SCORE_HDR_UNUSED}
};

static score_mapping_t score_pgs_mapping[] = {
    {"rsID", SCORE_HDR_SNP}, {"chr_position", SCORE_HDR_BP}, {"chr_name", SCORE_HDR_CHR}, {"effect_allele", SCORE_HDR_A1},
    {"other_allele", SCORE_HDR_UNUSED}, {"OR", SCORE_HDR_OR}, {"HR", SCORE_HDR_OR}, {"effect_weight", SCORE_HDR_BETA},
    {"allelefrequency_effect", SCORE_HDR_UNUSED}
};

static score_mapping_t score_ssf_mapping[] = {
    {"variant_id", SCORE_HDR_SNP}, {"rsid", SCORE_HDR_SNP}, {"rs_id", SCORE_HDR_SNP},
    {"base_pair_location", SCORE_HDR_BP}, {"chromosome", SCORE_HDR_CHR}, {"effect_allele", SCORE_HDR_A1},
    {"other_allele", SCORE_HDR_UNUSED}, {"p_value", SCORE_HDR_P}, {"odds_ratio", SCORE_HDR_OR},
    {"hazard_ratio", SCORE_HDR_OR}, {"beta", SCORE_HDR_BETA}, {"n", SCORE_HDR_UNUSED},
    {"info", SCORE_HDR_UNUSED}, {"effect_allele_frequency", SCORE_HDR_UNUSED},
    {"standard_error", SCORE_HDR_UNUSED}
};

static int score_is_missing(float f) {
    return isnan(f) || bcf_float_is_missing(f) || bcf_float_is_vector_end(f);
}

/* Detect whether a summary file is GWAS-VCF format by checking htslib file type.
 * Returns 1 if VCF/BCF, 0 otherwise. */
static int score_is_vcf_summary(const char *path) {
    htsFile *fp = hts_open(path, "r");
    int is_vcf = 0;
    if (!fp) return 0;
    if (fp->format.category == sequence_data || fp->format.category == variant_data) {
        enum htsExactFormat fmt = fp->format.format;
        if (fmt == vcf || fmt == bcf) is_vcf = 1;
    }
    hts_close(fp);
    return is_vcf;
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

static char *score_prs_name_from_path(const char *path) {
    static const char *ext_str[] = {"gz", "txt", "tsv", "vcf", "bcf"};
    static const int n_ext = (int)(sizeof(ext_str) / sizeof(ext_str[0]));
    char *tmp;
    const char *base;
    char *out;
    int stripped = 1;
    if (!path) return NULL;
    base = strrchr(path, '/');
    base = base ? base + 1 : path;
    tmp = score_dup(base);
    if (!tmp) return NULL;
    while (stripped) {
        char *ptr = strrchr(tmp, '.');
        int j;
        stripped = 0;
        if (!ptr) break;
        for (j = 0; j < n_ext; j++) {
            if (strcmp(ptr + 1, ext_str[j]) == 0) {
                *ptr = '\0';
                stripped = 1;
                break;
            }
        }
    }
    out = score_dup(tmp);
    duckdb_free(tmp);
    return out;
}

static int score_str_endswith(const char *s, const char *suffix) {
    size_t n, m;
    if (!s || !suffix) return 0;
    n = strlen(s);
    m = strlen(suffix);
    return n >= m && strcmp(s + n - m, suffix) == 0;
}

static int score_summary_dir_entry_supported(const char *name) {
    if (!name || !name[0]) return 0;
    if (name[0] == '.') return 0;
    /* Directory mode is a convenience for summary-statistic folders.  Restrict
     * it to regular summary-like filenames so index sidecars (.tbi/.csi) and
     * incidental files do not become PRS tracks. */
    if (score_str_endswith(name, ".bcf")) return 1;
    if (score_str_endswith(name, ".vcf") || score_str_endswith(name, ".vcf.gz") || score_str_endswith(name, ".vcf.bgz")) return 1;
    if (score_str_endswith(name, ".tsv") || score_str_endswith(name, ".tsv.gz") || score_str_endswith(name, ".tsv.bgz")) return 1;
    if (score_str_endswith(name, ".txt") || score_str_endswith(name, ".txt.gz") || score_str_endswith(name, ".txt.bgz")) return 1;
    if (score_str_endswith(name, ".ssf") || score_str_endswith(name, ".ssf.gz") || score_str_endswith(name, ".ssf.bgz")) return 1;
    return 0;
}

static int score_string_ptr_cmp(const void *a, const void *b) {
    const char *sa = *(const char * const *)a;
    const char *sb = *(const char * const *)b;
    if (!sa && !sb) return 0;
    if (!sa) return -1;
    if (!sb) return 1;
    return strcmp(sa, sb);
}

static int score_add_summary_path(score_bind_t *bind, const char *path, char *err, size_t err_sz) {
    char **new_paths;
    char *copy;
    if (!path || !path[0]) return 0;
    if (bind->n_summary_paths == bind->m_summary_paths) {
        int new_m = bind->m_summary_paths ? bind->m_summary_paths * 2 : 4;
        new_paths = (char **)realloc(bind->summary_paths, sizeof(char *) * (size_t)new_m);
        if (!new_paths) {
            snprintf(err, err_sz, "bcftools_score: out of memory");
            return -1;
        }
        bind->summary_paths = new_paths;
        bind->m_summary_paths = new_m;
    }
    copy = score_dup(path);
    if (!copy) {
        snprintf(err, err_sz, "bcftools_score: out of memory");
        return -1;
    }
    bind->summary_paths[bind->n_summary_paths++] = copy;
    return 0;
}

static int score_add_summary_paths_from_list_file(score_bind_t *bind, const char *path, char *err, size_t err_sz) {
    DIR *d;
    htsFile *fp;
    kstring_t str = {0, 0, NULL};
    int added = 0;
    if (!path || !path[0]) return 0;

    /* Directory mode intentionally uses only opendir/readdir/stat/qsort rather
     * than scandir or dirent.d_type: those calls are available in Emscripten's
     * virtual filesystem for webR/wasm, while remote URL list files fall through
     * to hts_open() below and are handled by the wasm HTTP hFILE backend. */
    d = opendir(path);
    if (d) {
        struct dirent *dir;
        char **entries = NULL;
        int n_entries = 0;
        int m_entries = 0;
        while ((dir = readdir(d))) {
            char *joined;
            char **new_entries;
            struct stat st;
            size_t p_len;
            size_t q_len;
            if (strcmp(dir->d_name, ".") == 0 || strcmp(dir->d_name, "..") == 0) continue;
            if (!score_summary_dir_entry_supported(dir->d_name)) continue;
            p_len = strlen(path);
            q_len = strlen(dir->d_name);
            joined = (char *)duckdb_malloc(p_len + 1 + q_len + 1);
            if (!joined) {
                int i;
                closedir(d);
                for (i = 0; i < n_entries; i++) duckdb_free(entries[i]);
                duckdb_free(entries);
                snprintf(err, err_sz, "bcftools_score: out of memory");
                return -1;
            }
            memcpy(joined, path, p_len);
            joined[p_len] = '/';
            memcpy(joined + p_len + 1, dir->d_name, q_len + 1);
            if (stat(joined, &st) != 0 || !S_ISREG(st.st_mode)) {
                duckdb_free(joined);
                continue;
            }
            if (n_entries == m_entries) {
                int new_m = m_entries ? m_entries * 2 : 8;
                new_entries = (char **)realloc(entries, sizeof(char *) * (size_t)new_m);
                if (!new_entries) {
                    int i;
                    duckdb_free(joined);
                    closedir(d);
                    for (i = 0; i < n_entries; i++) duckdb_free(entries[i]);
                    duckdb_free(entries);
                    snprintf(err, err_sz, "bcftools_score: out of memory");
                    return -1;
                }
                entries = new_entries;
                m_entries = new_m;
            }
            entries[n_entries++] = joined;
        }
        closedir(d);
        if (n_entries > 1) qsort(entries, (size_t)n_entries, sizeof(char *), score_string_ptr_cmp);
        for (added = 0; added < n_entries; added++) {
            if (score_add_summary_path(bind, entries[added], err, err_sz) != 0) {
                int i;
                for (i = 0; i < n_entries; i++) duckdb_free(entries[i]);
                duckdb_free(entries);
                return -1;
            }
        }
        {
            int i;
            for (i = 0; i < n_entries; i++) duckdb_free(entries[i]);
        }
        duckdb_free(entries);
        if (!added) {
            snprintf(err, err_sz, "bcftools_score: no supported summary files found in '%s'", path);
            return -1;
        }
        return 0;
    }

    fp = hts_open(path, "r");
    if (!fp) {
        snprintf(err, err_sz, "bcftools_score: failed to open summaries_list_file '%s'", path);
        return -1;
    }
    while (hts_getline(fp, '\n', &str) >= 0) {
        char *s = str.s;
        char *e;
        char *resolved;
        while (*s && isspace((unsigned char)*s)) s++;
        if (!*s || *s == '#') continue;
        e = s + strlen(s);
        while (e > s && isspace((unsigned char)e[-1])) e--;
        *e = '\0';
        resolved = score_dup(s);
        if (!resolved) {
            if (str.s) free(str.s);
            hts_close(fp);
            snprintf(err, err_sz, "bcftools_score: out of memory");
            return -1;
        }
        if (score_add_summary_path(bind, resolved, err, err_sz) != 0) {
            duckdb_free(resolved);
            if (str.s) free(str.s);
            hts_close(fp);
            return -1;
        }
        duckdb_free(resolved);
        added++;
    }
    if (str.s) free(str.s);
    hts_close(fp);
    if (!added) {
        snprintf(err, err_sz, "bcftools_score: no summary paths found in summaries_list_file '%s'", path);
        return -1;
    }
    return 0;
}

static int score_add_summary_paths_from_value(score_bind_t *bind, duckdb_value val, char *err, size_t err_sz) {
    duckdb_logical_type type;
    duckdb_type type_id;
    if (!val || duckdb_is_null_value(val)) return 0;
    type = duckdb_get_value_type(val);
    type_id = duckdb_get_type_id(type);
    if (type_id == DUCKDB_TYPE_LIST || type_id == DUCKDB_TYPE_ARRAY) {
        idx_t i;
        idx_t n = duckdb_get_list_size(val);
        if (n == 0) {
            snprintf(err, err_sz, "bcftools_score: summary_path list must not be empty");
            return -1;
        }
        for (i = 0; i < n; i++) {
            duckdb_value child = duckdb_get_list_child(val, i);
            char *path;
            if (!child || duckdb_is_null_value(child)) {
                if (child) duckdb_destroy_value(&child);
                snprintf(err, err_sz, "bcftools_score: summary_path list must not contain NULL");
                return -1;
            }
            path = duckdb_get_varchar(child);
            if (!path || !path[0]) {
                if (path) duckdb_free(path);
                duckdb_destroy_value(&child);
                snprintf(err, err_sz, "bcftools_score: summary_path list must contain non-empty paths");
                return -1;
            }
            if (score_add_summary_path(bind, path, err, err_sz) != 0) {
                duckdb_free(path);
                duckdb_destroy_value(&child);
                return -1;
            }
            duckdb_free(path);
            duckdb_destroy_value(&child);
        }
        return 0;
    } else {
        char *path = duckdb_get_varchar(val);
        if (!path || !path[0]) {
            if (path) duckdb_free(path);
            return 0;
        }
        if (score_add_summary_path(bind, path, err, err_sz) != 0) {
            duckdb_free(path);
            return -1;
        }
        duckdb_free(path);
        return 0;
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
    int i;
    if (!bind) return;
    if (bind->bcf_path) duckdb_free(bind->bcf_path);
    if (bind->summary_paths) {
        for (i = 0; i < bind->n_summary_paths; i++) {
            if (bind->summary_paths[i]) duckdb_free(bind->summary_paths[i]);
        }
        duckdb_free(bind->summary_paths);
    }
    if (bind->summary_prs_counts) duckdb_free(bind->summary_prs_counts);
    if (bind->summaries_list_file) duckdb_free(bind->summaries_list_file);
    if (bind->log_path) duckdb_free(bind->log_path);
    if (bind->columns_preset) duckdb_free(bind->columns_preset);
    if (bind->columns_file) duckdb_free(bind->columns_file);
    if (bind->samples) duckdb_free(bind->samples);
    if (bind->regions) duckdb_free(bind->regions);
    if (bind->regions_file) duckdb_free(bind->regions_file);
    if (bind->targets) duckdb_free(bind->targets);
    if (bind->targets_file) duckdb_free(bind->targets_file);
    if (bind->apply_filters) duckdb_free(bind->apply_filters);
    if (bind->include_expr) duckdb_free(bind->include_expr);
    if (bind->exclude_expr) duckdb_free(bind->exclude_expr);
    if (bind->q_thr_lp) duckdb_free(bind->q_thr_lp);
    if (bind->prs_names) {
        for (i = 0; i < bind->n_prs; i++) {
            if (bind->prs_names[i]) duckdb_free(bind->prs_names[i]);
        }
        duckdb_free(bind->prs_names);
    }
    if (bind->prs_paths) {
        for (i = 0; i < bind->n_prs; i++) {
            if (bind->prs_paths[i]) duckdb_free(bind->prs_paths[i]);
        }
        duckdb_free(bind->prs_paths);
    }
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

static score_summary_t *score_load_summary(const score_bind_t *bind, bcf_hdr_t *bcf_hdr, const char *summary_path, char *err, size_t err_sz) {
    htsFile *fp = NULL;
    kstring_t str = {0, 0, NULL};
    char *fields[512];
    char delimiter = '\t';
    score_mapping_t *mapping = NULL;
    int mapping_n = 0;
    int mapping_owns = 0;
    int idx_snp = -1, idx_chr = -1, idx_bp = -1, idx_a1 = -1, idx_beta = -1, idx_or = -1, idx_p = -1, idx_lp = -1;
    score_summary_t *summary = NULL;
    char *pending_allele = NULL;
    uint64_t line_no = 0;

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

    /* Use hts_open for transparent .gz/.bgz decompression (S-C3).
     * Upstream: score.c:196. */
    fp = hts_open(summary_path, "r");
    if (!fp) {
        snprintf(err, err_sz, "bcftools_score: cannot open summary '%s': %s", summary_path, strerror(errno));
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }

    while (hts_getline(fp, '\n', &str) >= 0) {
        line_no++;
        if (str.s[0] == '#' && strncmp(str.s, "#CHR", 4) != 0 && strncmp(str.s, "#ID", 3) != 0) continue;
        break;
    }
    if (!str.s || str.l == 0) {
        snprintf(err, err_sz, "bcftools_score: empty summary file '%s'", summary_path);
        hts_close(fp);
        free(str.s);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }
    if (str.s[0] == '#') { memmove(str.s, str.s + 1, str.l - 1); str.l--; str.s[str.l] = '\0'; }
    delimiter = score_detect_delimiter(str.s);
    {
        int i;
        int n_fields = score_split_fields(str.s, delimiter, fields, 512);
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
    if (!summary) goto oom;
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
        hts_close(fp);
        free(str.s);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }
    if (bind->use_variant_id && !summary->snp_ok) {
        snprintf(err, err_sz, "bcftools_score: summary must include marker name column when use_variant_id=true");
        score_destroy_summary(summary);
        hts_close(fp);
        free(str.s);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }
    if (summary->use_snp && !summary->snp_ok) {
        snprintf(err, err_sz, "bcftools_score: summary must include marker name column when matching by variant ID");
        score_destroy_summary(summary);
        hts_close(fp);
        free(str.s);
        if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
        return NULL;
    }

    if (summary->use_snp) {
        summary->id2idx = khash_str2int_init();
        if (!summary->id2idx) goto oom;
    } else {
        summary->rid_pos2idx = kh_init(score64);
        if (!summary->rid_pos2idx) goto oom;
    }

    while (hts_getline(fp, '\n', &str) >= 0) {
        int n_fields;
        float es = NAN;
        float lp = NAN;
        score_marker_t marker;
        char *id = NULL;
        int rid = -1;
        int64_t pos = -1;

        line_no++;
        marker.a1 = NULL;
        if (str.s[0] == '#') continue;
        summary->all_markers++;
        n_fields = score_split_fields(str.s, delimiter, fields, 512);
        if (n_fields == 0) continue;

        if (idx_a1 >= n_fields) continue;

        if (idx_beta >= 0 && idx_beta < n_fields) {
            if (score_parse_float(fields[idx_beta], &es) != 0) es = NAN;
        }
        if (isnan(es) && idx_or >= 0 && idx_or < n_fields) {
            float orv = NAN;
            if (score_parse_float(fields[idx_or], &orv) == 0 && !isnan(orv) && orv > 0.0f)
                es = logf(orv);
        }
        if (isnan(es)) continue;

        if (idx_lp >= 0 && idx_lp < n_fields) {
            if (score_parse_float(fields[idx_lp], &lp) != 0) lp = NAN;
        }
        if (isnan(lp) && idx_p >= 0 && idx_p < n_fields) {
            if (score_parse_pvalue_lp(fields[idx_p], &lp) != 0) lp = NAN;
        }

        if (summary->use_snp) {
            if (idx_snp < 0 || idx_snp >= n_fields || fields[idx_snp][0] == '\0') continue;
        } else {
            char *end;
            if (idx_chr < 0 || idx_chr >= n_fields || idx_bp < 0 || idx_bp >= n_fields) continue;
            rid = score_hdr_name2id_flexible(bcf_hdr, fields[idx_chr]);
            if (rid < 0) continue;
            pos = strtoll(fields[idx_bp], &end, 10) - 1;
            if (end == fields[idx_bp] || *end != '\0' || pos < 0) continue;
        }

        pending_allele = score_dup(fields[idx_a1]);
        if (!pending_allele) goto oom;
        score_strtoupper(pending_allele);
        marker.a1 = pending_allele;
        marker.es = es;
        marker.lp = lp;

        if (summary->n_markers == summary->m_markers) {
            if (summary->m_markers > INT_MAX / 2) goto oom;
            int new_m = summary->m_markers ? summary->m_markers * 2 : 1024;
            if ((size_t)new_m > SIZE_MAX / sizeof(score_marker_t)) goto oom;
            score_marker_t *new_markers = (score_marker_t *)realloc(summary->markers, sizeof(score_marker_t) * (size_t)new_m);
            if (!new_markers) goto oom;
            summary->markers = new_markers;
            summary->m_markers = new_m;
        }

        if (summary->use_snp) {
            int size;
            int ret;
            id = score_dup(fields[idx_snp]);
            if (!id) goto oom;
            size = khash_str2int_size(summary->id2idx);
            ret = khash_str2int_inc(summary->id2idx, id);
            if (ret < 0) {
                duckdb_free(id);
                goto oom;
            }
            if (ret == size) {
                summary->markers[summary->n_markers++] = marker;
            } else {
                summary->duplicate_markers++;
                duckdb_free(id);
                if (marker.a1) duckdb_free(marker.a1);
            }
        } else {
            khash_t(score64) *hash;
            khiter_t it;
            int ret;
            uint64_t key;
            key = (((uint64_t)(uint32_t)rid) << 44) | (uint64_t)pos;
            hash = (khash_t(score64) *)summary->rid_pos2idx;
            it = kh_put(score64, hash, key, &ret);
            if (ret < 0) goto oom;
            if (ret > 0) {
                kh_val(hash, it) = summary->n_markers;
                summary->markers[summary->n_markers++] = marker;
            } else {
                summary->duplicate_markers++;
                if (marker.a1) duckdb_free(marker.a1);
            }
        }
        pending_allele = NULL;
    }

    hts_close(fp);
    free(str.s);
    if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
    return summary;

oom:
    snprintf(err, err_sz, "bcftools_score: out of memory loading summary '%s' at line %llu",
             summary_path, (unsigned long long)line_no);
    duckdb_free(pending_allele);
    score_destroy_summary(summary);
    hts_close(fp);
    free(str.s);
    if (mapping_owns) score_destroy_mapping(mapping, mapping_n, 1);
    return NULL;
}

/* Load GWAS-VCF summary: reads FORMAT/ES and FORMAT/LP per-sample from a VCF/BCF
 * summary file and builds one score_summary_t per PRS (= GWAS-VCF sample).
 * Port of upstream score.c:681-710 (synced reader VCF path), adapted to pre-load
 * into the TSV-mode summary structure for scoring.
 *
 * Returns array of n_prs summaries on success, NULL on error. */
static score_summary_t **score_load_summary_vcf(const score_bind_t *bind,
                                                 bcf_hdr_t *geno_hdr,
                                                 const char *summary_path,
                                                 int n_prs,
                                                 char *err, size_t err_sz) {
    htsFile *sfp = NULL;
    bcf_hdr_t *shdr = NULL;
    bcf1_t *rec = NULL;
    score_summary_t **summaries = NULL;
    float *float_arr = NULL;
    float *lp_arr = NULL;
    char *pending_allele = NULL;
    int m_float = 0;
    int m_lp = 0;
    uint64_t record_no = 0;
    int i;

    sfp = hts_open(summary_path, "r");
    if (!sfp) {
        snprintf(err, err_sz, "bcftools_score: cannot open GWAS-VCF summary '%s'", summary_path);
        return NULL;
    }
    shdr = bcf_hdr_read(sfp);
    if (!shdr) {
        snprintf(err, err_sz, "bcftools_score: cannot read GWAS-VCF header '%s'", summary_path);
        hts_close(sfp);
        return NULL;
    }

    /* Validate FORMAT/ES exists */
    {
        int es_id = bcf_hdr_id2int(shdr, BCF_DT_ID, "ES");
        if (!bcf_hdr_idinfo_exists(shdr, BCF_HL_FMT, es_id)) {
            snprintf(err, err_sz, "bcftools_score: GWAS-VCF '%s' does not include the ES FORMAT field", summary_path);
            bcf_hdr_destroy(shdr);
            hts_close(sfp);
            return NULL;
        }
    }
    /* Validate FORMAT/LP if q_score_thr is set */
    if (bind->q_thr_lp) {
        int lp_id = bcf_hdr_id2int(shdr, BCF_DT_ID, "LP");
        if (!bcf_hdr_idinfo_exists(shdr, BCF_HL_FMT, lp_id)) {
            snprintf(err, err_sz, "bcftools_score: GWAS-VCF '%s' does not include the LP FORMAT field", summary_path);
            bcf_hdr_destroy(shdr);
            hts_close(sfp);
            return NULL;
        }
    }

    if (bcf_hdr_nsamples(shdr) != n_prs) {
        snprintf(err, err_sz, "bcftools_score: GWAS-VCF sample count mismatch (expected %d, got %d)", n_prs, bcf_hdr_nsamples(shdr));
        bcf_hdr_destroy(shdr);
        hts_close(sfp);
        return NULL;
    }

    /* Allocate per-PRS summaries */
    summaries = (score_summary_t **)duckdb_malloc(sizeof(score_summary_t *) * (size_t)n_prs);
    if (!summaries) goto oom;
    memset(summaries, 0, sizeof(score_summary_t *) * (size_t)n_prs);
    for (i = 0; i < n_prs; i++) {
        summaries[i] = (score_summary_t *)duckdb_malloc(sizeof(score_summary_t));
        if (!summaries[i]) goto oom;
        memset(summaries[i], 0, sizeof(score_summary_t));
        summaries[i]->use_snp = 0; /* match by chr:pos */
        summaries[i]->chr_ok = 1;
        summaries[i]->bp_ok = 1;
        summaries[i]->a1_ok = 1;
        summaries[i]->beta_ok = 1;
        summaries[i]->p_ok = bind->q_thr_lp ? 1 : 0;
        summaries[i]->rid_pos2idx = kh_init(score64);
        if (!summaries[i]->rid_pos2idx) goto oom;
    }

    rec = bcf_init1();
    if (!rec) goto oom;
    while (bcf_read(sfp, shdr, rec) == 0) {
        int n_es, n_lp = 0;
        int rid;
        uint64_t key;
        const char *chrom;

        record_no++;
        bcf_unpack(rec, BCF_UN_ALL);
        for (i = 0; i < n_prs; i++) summaries[i]->all_markers++;
        if (rec->n_allele < 2) continue;

        chrom = bcf_hdr_id2name(shdr, rec->rid);
        rid = score_hdr_name2id_flexible(geno_hdr, chrom);
        if (rid < 0) continue;

        key = (((uint64_t)(uint32_t)rid) << 44) | (uint64_t)rec->pos;

        /* htslib 1.24 returns -4 when a FORMAT buffer allocation fails. */
        n_es = bcf_get_format_float(shdr, rec, "ES", &float_arr, &m_float);
        if (n_es == -4) goto oom;
        if (n_es < n_prs) continue;

        /* Read FORMAT/LP if needed */
        if (bind->q_thr_lp) {
            n_lp = bcf_get_format_float(shdr, rec, "LP", &lp_arr, &m_lp);
            if (n_lp == -4) goto oom;
            if (n_lp < n_prs) continue;
        }

        /* Add marker to each PRS summary */
        for (i = 0; i < n_prs; i++) {
            float es = float_arr[i];
            float lp = (lp_arr && i < n_lp) ? lp_arr[i] : NAN;
            score_marker_t marker;
            khash_t(score64) *hash;
            khiter_t it;
            int ret;

            if (score_is_missing(es)) continue;
            if (bind->q_thr_lp && score_is_missing(lp)) continue;

            /* Effect allele = ALT (allele index 1) — upstream convention */
            pending_allele = score_dup(rec->d.allele[1]);
            if (!pending_allele) goto oom;
            score_strtoupper(pending_allele);
            marker.a1 = pending_allele;
            marker.es = es;
            marker.lp = lp;

            /* Grow markers array if needed */
            if (summaries[i]->n_markers == summaries[i]->m_markers) {
                if (summaries[i]->m_markers > INT_MAX / 2) goto oom;
                int new_m = summaries[i]->m_markers ? summaries[i]->m_markers * 2 : 1024;
                if ((size_t)new_m > SIZE_MAX / sizeof(score_marker_t)) goto oom;
                score_marker_t *new_mk = (score_marker_t *)realloc(summaries[i]->markers,
                                                                     sizeof(score_marker_t) * (size_t)new_m);
                if (!new_mk) goto oom;
                summaries[i]->markers = new_mk;
                summaries[i]->m_markers = new_m;
            }

            hash = (khash_t(score64) *)summaries[i]->rid_pos2idx;
            it = kh_put(score64, hash, key, &ret);
            if (ret < 0) goto oom;
            if (ret > 0) {
                kh_val(hash, it) = summaries[i]->n_markers;
                summaries[i]->markers[summaries[i]->n_markers++] = marker;
            } else {
                summaries[i]->duplicate_markers++;
                duckdb_free(marker.a1);
            }
            pending_allele = NULL;
        }
    }

    bcf_destroy1(rec);
    free(float_arr);
    free(lp_arr);
    bcf_hdr_destroy(shdr);
    hts_close(sfp);
    return summaries;

oom:
    snprintf(err, err_sz, "bcftools_score: out of memory loading GWAS-VCF summary '%s' at record %llu",
             summary_path, (unsigned long long)record_no);
    duckdb_free(pending_allele);
    if (summaries) {
        for (i = 0; i < n_prs; i++) score_destroy_summary(summaries[i]);
        duckdb_free(summaries);
    }
    bcf_destroy1(rec);
    free(float_arr);
    free(lp_arr);
    bcf_hdr_destroy(shdr);
    hts_close(sfp);
    return NULL;
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

    /* Validate that the selected tag (whether auto-detected or explicitly
     * requested via use := '...') actually exists in the VCF header.
     * Port of upstream score.c:767-795. */
    switch (bind->use_tag) {
    case SCORE_USE_GT:
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id)) {
            snprintf(err, err_sz, "bcftools_score: VCF does not include the GT FORMAT field");
            return -1;
        }
        break;
    case SCORE_USE_DS:
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ds_id)) {
            snprintf(err, err_sz, "bcftools_score: VCF does not include the DS FORMAT field");
            return -1;
        }
        break;
    case SCORE_USE_HDS:
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, hds_id)) {
            snprintf(err, err_sz, "bcftools_score: VCF does not include the HDS FORMAT field");
            return -1;
        }
        break;
    case SCORE_USE_AP:
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap1_id) || !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, ap2_id)) {
            snprintf(err, err_sz, "bcftools_score: VCF does not include either the AP1 or the AP2 FORMAT fields");
            return -1;
        }
        break;
    case SCORE_USE_GP:
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gp_id)) {
            snprintf(err, err_sz, "bcftools_score: VCF does not include the GP FORMAT field");
            return -1;
        }
        break;
    case SCORE_USE_AS:
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, as_id)) {
            snprintf(err, err_sz, "bcftools_score: VCF does not include the AS FORMAT field");
            return -1;
        }
        break;
    }
    return 0;
}

static int score_record_passes_filter_list(const bcf_hdr_t *hdr, bcf1_t *rec, const char *apply_filters) {
    char *tmp = NULL;
    char *tok = NULL;
    int pass = 0;
    int i;

    if (!apply_filters || apply_filters[0] == '\0') return 1;
    tmp = score_dup(apply_filters);
    if (!tmp) return 0;

    if (rec->d.n_flt == 0) {
        for (tok = strtok(tmp, ","); tok; tok = strtok(NULL, ",")) {
            if (strcmp(tok, ".") == 0 || strcasecmp(tok, "PASS") == 0) {
                pass = 1;
                break;
            }
        }
        duckdb_free(tmp);
        return pass;
    }

    for (i = 0; i < rec->d.n_flt && !pass; i++) {
        const char *flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
        tok = strtok(tmp, ",");
        while (tok) {
            if (strcmp(tok, flt) == 0) {
                pass = 1;
                break;
            }
            tok = strtok(NULL, ",");
        }
        if (!pass) strcpy(tmp, apply_filters);
    }

    duckdb_free(tmp);
    return pass;
}

typedef struct {
    const bcf_hdr_t *hdr;
    char *missing_name; /* owned copy of the most recent name absent from the header */
    int lookup_failed;  /* the most recent lookup found no header contig */
    int out_of_memory;
} score_region_lookup_t;

static int score_region_name2id(void *data, const char *name) {
    score_region_lookup_t *lookup = (score_region_lookup_t *)data;
    int id = bcf_hdr_name2id(lookup->hdr, name);
    lookup->lookup_failed = id < 0;
    if (id < 0) {
        free(lookup->missing_name);
        lookup->missing_name = strdup(name);
        if (!lookup->missing_name) lookup->out_of_memory = 1;
    }
    return id;
}

/* Resolves only the recorded missing name, so a second parse checks the
 * coordinates of an item whose contig the header lacks. */
static int score_region_missing_name2id(void *data, const char *name) {
    return strcmp((const char *)data, name) == 0 ? 0 : -1;
}

/* Start of the next list item after `item`, splitting like hts_parse_region with
 * HTS_PARSE_LIST: a braced name ends at its closing brace, then the item at the next
 * comma. Returns NULL for an unterminated brace. */
static const char *score_region_next_item(const char *item) {
    const char *from = item;
    const char *comma;
    if (*item == '{') {
        from = strchr(item, '}');
        if (!from) return NULL;
    }
    comma = strchr(from, ',');
    return comma ? comma + 1 : item + strlen(item);
}

static int score_interval_compare(const void *a, const void *b) {
    const score_interval_t *left = (const score_interval_t *)a;
    const score_interval_t *right = (const score_interval_t *)b;
    if (left->rid != right->rid) return left->rid < right->rid ? -1 : 1;
    if (left->beg != right->beg) return left->beg < right->beg ? -1 : 1;
    return (left->end > right->end) - (left->end < right->end);
}

/* Store disjoint half-open intervals ordered by header RID and start. Merging
 * also makes their ends increasing within each contig. */
static int score_parse_region_list(const bcf_hdr_t *hdr, const char *regions, score_region_list_t *list) {
    const char *cursor = regions;
    score_region_lookup_t lookup = { hdr, NULL, 0, 0 };
    size_t capacity = 0;
    size_t merged = 0;
    int status = 0;

    if (!regions || !regions[0]) return 0;
    list->active = 1;
    while (*cursor) {
        int tid;
        hts_pos_t beg, end;
        const char *next;
        /* Use HTSlib's bcftools-style coordinates and braced contig quoting,
         * including names containing colons or commas. */
        lookup.lookup_failed = 0;
        next = hts_parse_region(cursor, &tid, &beg, &end, score_region_name2id,
                                &lookup, HTS_PARSE_LIST | HTS_PARSE_ONE_COORD);
        if (lookup.out_of_memory) {
            status = -1;
            break;
        }
        if (!next) {
            /* A contig absent from this VCF matches no record, as in the htslib
             * synced reader. hts_parse_region stops at the failed lookup, so parse
             * the item again with that name resolved: bad coordinates or an
             * unterminated brace remain errors. */
            if (tid != -1 || !lookup.lookup_failed || !lookup.missing_name ||
                !hts_parse_region(cursor, &tid, &beg, &end, score_region_missing_name2id,
                                  lookup.missing_name, HTS_PARSE_LIST | HTS_PARSE_ONE_COORD)) {
                status = -1;
                break;
            }
            next = score_region_next_item(cursor);
            if (!next) {
                status = -1;
                break;
            }
            cursor = next;
            continue;
        }
        if (list->count == capacity) {
            size_t new_capacity = capacity ? capacity * 2 : 16;
            score_interval_t *intervals;
            if (new_capacity < capacity || new_capacity > SIZE_MAX / sizeof(*intervals)) {
                status = -1;
                break;
            }
            intervals = (score_interval_t *)realloc(list->intervals, new_capacity * sizeof(*intervals));
            if (!intervals) {
                status = -1;
                break;
            }
            list->intervals = intervals;
            capacity = new_capacity;
        }
        list->intervals[list->count++] = (score_interval_t){tid, beg, end};
        cursor = next;
    }
    free(lookup.missing_name);
    if (status < 0 || list->count == 0) return status;

    qsort(list->intervals, list->count, sizeof(*list->intervals), score_interval_compare);
    for (size_t i = 0; i < list->count; i++) {
        score_interval_t *current = &list->intervals[i];
        if (merged && list->intervals[merged - 1].rid == current->rid &&
            list->intervals[merged - 1].end >= current->beg) {
            if (list->intervals[merged - 1].end < current->end) {
                list->intervals[merged - 1].end = current->end;
            }
        } else {
            list->intervals[merged++] = *current;
        }
    }
    list->count = merged;
    return 0;
}

static int score_record_in_regions(const bcf1_t *rec, const score_region_list_t *list, int overlap_mode) {
    hts_pos_t end;
    size_t lo = 0, hi = list->count;
    if (!list->active) return 1;
    end = rec->pos + (overlap_mode == 0 ? 1 : (rec->rlen > 0 ? rec->rlen : 1));
    /* The last interval starting before the record's exclusive end is the only
     * candidate needed: earlier disjoint intervals also have earlier ends. */
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        const score_interval_t *interval = &list->intervals[mid];
        if (interval->rid < rec->rid || (interval->rid == rec->rid && interval->beg < end)) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo > 0 && list->intervals[lo - 1].rid == rec->rid && list->intervals[lo - 1].end > rec->pos;
}

static int score_append_token(char **out, size_t *len, size_t *cap, const char *tok) {
    size_t tok_len;
    size_t need;
    char *next;
    if (!tok || !tok[0]) return 0;
    tok_len = strlen(tok);
    need = *len + tok_len + ((*len > 0) ? 1 : 0) + 1;
    if (need > *cap) {
        size_t new_cap = *cap ? *cap : 256;
        while (new_cap < need) new_cap *= 2;
        next = (char *)realloc(*out, new_cap);
        if (!next) return -1;
        *out = next;
        *cap = new_cap;
    }
    if (*len > 0) {
        (*out)[*len] = ',';
        *len += 1;
    }
    memcpy(*out + *len, tok, tok_len);
    *len += tok_len;
    (*out)[*len] = '\0';
    return 0;
}

static int score_parse_region_line(const char *line, int is_bed, char *token, size_t token_sz) {
    char buf[2048];
    char *fields[8];
    int n_fields;
    if (!line || strlen(line) >= sizeof(buf)) return -1;
    snprintf(buf, sizeof(buf), "%s", line);
    n_fields = score_split_fields(buf, '\t', fields, 8);
    if (n_fields <= 0) return -1;
    if (n_fields == 1) {
        if (is_bed || strlen(fields[0]) >= token_sz) return -1;
        snprintf(token, token_sz, "%s", fields[0]);
        return 0;
    }
    {
        const char *chr = fields[0];
        char *tail;
        int written;
        errno = 0;
        int64_t beg = strtoll(fields[1], &tail, 10);
        if (errno == ERANGE || tail == fields[1] || *tail) return -1;
        int64_t end = beg;
        if (is_bed && n_fields < 3) return -1;
        if (n_fields >= 3) {
            errno = 0;
            end = strtoll(fields[2], &tail, 10);
            if (errno == ERANGE || tail == fields[2] || *tail) return -1;
        }
        if (is_bed) {
            if (beg < 0 || beg == INT64_MAX) return -1;
            beg++;
        }
        if (beg < 1 || end < beg) return -1;
        written = snprintf(token, token_sz, "{%s}:%lld-%lld", chr, (long long)beg, (long long)end);
        if (written < 0 || (size_t)written >= token_sz) return -1;
    }
    return 0;
}

static char *score_regions_from_file(const char *path, char *err, size_t err_sz) {
    htsFile *fp;
    kstring_t str = {0, 0, NULL};
    char *out = NULL;
    size_t len = 0;
    size_t cap = 0;
    int have = 0;
    int ret;
    size_t path_len;
    int is_bed;

    if (!path || !path[0]) return NULL;
    path_len = strlen(path);
    is_bed = (path_len >= 4 && strcasecmp(path + path_len - 4, ".bed") == 0) ||
             (path_len >= 7 && strcasecmp(path + path_len - 7, ".bed.gz") == 0);
    fp = hts_open(path, "r");
    if (!fp) {
        snprintf(err, err_sz, "bcftools_score: failed to open region/target file: %s", path);
        return NULL;
    }

    while ((ret = hts_getline(fp, '\n', &str)) >= 0) {
        char token[512];
        char *s = str.s;
        while (*s && isspace((unsigned char)*s)) s++;
        if (!*s || *s == '#') continue;
        int parsed = score_parse_region_line(s, is_bed, token, sizeof(token));
        if (parsed < 0) {
            if (str.s) free(str.s);
            hts_close(fp);
            if (out) duckdb_free(out);
            snprintf(err, err_sz, "bcftools_score: invalid region/target line in file: %s", path);
            return NULL;
        }
        if (score_append_token(&out, &len, &cap, token) != 0) {
            if (str.s) free(str.s);
            hts_close(fp);
            if (out) duckdb_free(out);
            snprintf(err, err_sz, "bcftools_score: out of memory");
            return NULL;
        }
        have = 1;
    }

    if (str.s) free(str.s);
    hts_close(fp);

    if (ret < -1) {
        if (out) duckdb_free(out);
        snprintf(err, err_sz, "bcftools_score: failed to read region/target file: %s", path);
        return NULL;
    }
    if (!have) {
        if (out) duckdb_free(out);
        return score_dup("");
    }
    return out;
}

static int score_record_matches_exprs(bcf1_t *rec, const score_scan_source_t *src, char *err, size_t err_sz) {
    if (src->include_filt) {
        int keep;
        keep = filter_test(src->include_filt, rec, NULL);
        if (keep < 0) {
            const char *detail = filter_last_error(src->include_filt);
            snprintf(err, err_sz, "bcftools_score: failed to evaluate include expression%s%s",
                     detail && detail[0] ? ": " : "", detail && detail[0] ? detail : "");
            return -1;
        }
        if (!keep) return 0;
    }
    if (src->exclude_filt) {
        int drop;
        drop = filter_test(src->exclude_filt, rec, NULL);
        if (drop < 0) {
            const char *detail = filter_last_error(src->exclude_filt);
            snprintf(err, err_sz, "bcftools_score: failed to evaluate exclude expression%s%s",
                     detail && detail[0] ? ": " : "", detail && detail[0] ? detail : "");
            return -1;
        }
        if (drop) return 0;
    }
    return 1;
}

static int score_init_source(const score_bind_t *bind, score_scan_source_t *src, bcf_hdr_t **out_hdr, char *err, size_t err_sz) {
    bcf_hdr_t *hdr0;
    const char **undef = NULL;
    int nundef = 0;
    src->sr = bcf_sr_init();
    if (!src->sr) {
        snprintf(err, err_sz, "bcftools_score: failed to initialize synced reader");
        return -1;
    }
    bcf_sr_set_opt(src->sr, BCF_SR_ALLOW_NO_IDX);

    if (bind->apply_filters && bind->apply_filters[0]) src->sr->apply_filters = bind->apply_filters;
    if (bind->regions_overlap >= 0) bcf_sr_set_opt(src->sr, BCF_SR_REGIONS_OVERLAP, bind->regions_overlap);
    if (bind->targets_overlap >= 0) bcf_sr_set_opt(src->sr, BCF_SR_TARGETS_OVERLAP, bind->targets_overlap);

    if (!bcf_sr_add_reader(src->sr, bind->bcf_path)) {
        snprintf(err, err_sz, "bcftools_score: failed to open BCF/VCF input: %s", bind->bcf_path);
        return -1;
    }

    hdr0 = bcf_sr_get_header(src->sr, 0);
    if (!hdr0) {
        snprintf(err, err_sz, "bcftools_score: failed to read BCF/VCF header");
        return -1;
    }

    if (bind->include_expr && bind->include_expr[0]) {
        char *include_expr = score_to_filter_expr(bind->include_expr);
        if (!include_expr) {
            snprintf(err, err_sz, "bcftools_score: out of memory");
            return -1;
        }
        if (!include_expr[0]) {
            free(include_expr);
            include_expr = NULL;
        }
        if (!include_expr) {
            src->include_filt = NULL;
        } else {
        src->include_filt = filter_parse(hdr0, include_expr);
        free(include_expr);
        if (!src->include_filt) {
            snprintf(err, err_sz, "bcftools_score: failed to evaluate include expression");
            return -1;
        }
        if (filter_status(src->include_filt) != FILTER_OK) {
            undef = filter_list_undef_tags(src->include_filt, &nundef);
            (void)undef;
            (void)nundef;
            if (filter_last_error(src->include_filt) && filter_last_error(src->include_filt)[0]) {
                snprintf(err, err_sz, "bcftools_score: failed to evaluate include expression: %s", filter_last_error(src->include_filt));
            } else {
                snprintf(err, err_sz, "bcftools_score: failed to evaluate include expression");
            }
            return -1;
        }
        }
    }
    if (bind->exclude_expr && bind->exclude_expr[0]) {
        char *exclude_expr = score_to_filter_expr(bind->exclude_expr);
        if (!exclude_expr) {
            snprintf(err, err_sz, "bcftools_score: out of memory");
            return -1;
        }
        if (!exclude_expr[0]) {
            free(exclude_expr);
            exclude_expr = NULL;
        }
        if (!exclude_expr) {
            src->exclude_filt = NULL;
        } else {
        src->exclude_filt = filter_parse(hdr0, exclude_expr);
        free(exclude_expr);
        if (!src->exclude_filt) {
            snprintf(err, err_sz, "bcftools_score: failed to evaluate exclude expression");
            return -1;
        }
        if (filter_status(src->exclude_filt) != FILTER_OK) {
            undef = filter_list_undef_tags(src->exclude_filt, &nundef);
            (void)undef;
            (void)nundef;
            if (filter_last_error(src->exclude_filt) && filter_last_error(src->exclude_filt)[0]) {
                snprintf(err, err_sz, "bcftools_score: failed to evaluate exclude expression: %s", filter_last_error(src->exclude_filt));
            } else {
                snprintf(err, err_sz, "bcftools_score: failed to evaluate exclude expression");
            }
            return -1;
        }
        }
    }

    *out_hdr = bcf_sr_get_header(src->sr, 0);
    if (!*out_hdr) {
        snprintf(err, err_sz, "bcftools_score: failed to read BCF/VCF header");
        return -1;
    }
    return 0;
}

static int score_source_next(score_scan_source_t *src, bcf1_t **out_rec) {
    while (bcf_sr_next_line(src->sr)) {
        bcf1_t *line;
        if (!bcf_sr_has_line(src->sr, 0)) continue;
        line = bcf_sr_get_line(src->sr, 0);
        *out_rec = line;
        return 1;
    }
    return 0;
}

static void score_source_destroy(score_scan_source_t *src) {
    if (!src) return;
    free(src->regions.intervals);
    free(src->targets.intervals);
    if (src->include_filt) filter_destroy(src->include_filt);
    if (src->exclude_filt) filter_destroy(src->exclude_filt);
    if (src->sr) {
        bcf_sr_destroy(src->sr);
        src->sr = NULL;
    }
}

static const char *score_use_tag_label(int use_tag) {
    switch (use_tag) {
    case SCORE_USE_GT: return "genotypes (GT)";
    case SCORE_USE_DS: return "genotype dosages (DS)";
    case SCORE_USE_HDS: return "haploid alternate allele dosage (HDS)";
    case SCORE_USE_AP: return "ALT haplotype probabilities (AP)";
    case SCORE_USE_GP: return "genotype probabilities (GP)";
    case SCORE_USE_AS: return "allelic shifts (AS)";
    default: return "auto-detected";
    }
}

static int score_format_metric_name(const score_bind_t *bind,
                                    const char *pname,
                                    int threshold_idx,
                                    int is_count,
                                    char *out,
                                    size_t out_sz) {
    int n;
    if (!out || out_sz == 0) return -1;
    if (bind->q_thr_lp && bind->n_q_thr > 0) {
        double p = pow(10.0, -bind->q_thr_lp[threshold_idx]);
        n = snprintf(out, out_sz, is_count ? "%s_CNT_p%.6g" : "%s_p%.6g", pname, p);
    } else {
        n = snprintf(out, out_sz, is_count ? "%s_CNT" : "%s", pname);
    }
    return (n < 0 || (size_t)n >= out_sz) ? -1 : 0;
}

static int score_validate_prs_names(const score_bind_t *bind, char *err, size_t err_sz) {
    duckhts_column_names columns;
    char detail[300];
    int rc = 0;
    if (!bind || !bind->prs_names) return 0;

    duckhts_column_names_init(&columns, 0);
    if (duckhts_column_names_claim(&columns, "SAMPLE", NULL, detail, sizeof(detail)) != 0) {
        snprintf(err, err_sz, "bcftools_score: %s", detail);
        rc = -1;
    }
    for (int i = 0; rc == 0 && i < bind->n_prs; i++) {
        if (!bind->prs_names[i] || !bind->prs_names[i][0]) {
            snprintf(err, err_sz, "bcftools_score: summary PRS name must not be empty");
            rc = -1;
        }
        for (int j = 0; rc == 0 && j < bind->n_q_thr; j++) {
            for (int k = 0; rc == 0 && k < (bind->counts ? 2 : 1); k++) {
                char metric_name[256];
                if (score_format_metric_name(bind, bind->prs_names[i], j, k == 1, metric_name, sizeof(metric_name)) != 0) {
                    snprintf(err, err_sz, "bcftools_score: generated output column name for PRS '%s' is too long", bind->prs_names[i]);
                    rc = -1;
                } else if (!metric_name[0]) {
                    snprintf(err, err_sz, "bcftools_score: generated output column name must not be empty");
                    rc = -1;
                } else if (duckhts_column_names_claim(&columns, metric_name, NULL, detail, sizeof(detail)) != 0) {
                    snprintf(err, err_sz, "bcftools_score: duplicate generated output column name: %s; use unique summary filenames or GWAS-VCF samples", detail);
                    rc = -1;
                }
            }
        }
    }
    duckhts_column_names_destroy(&columns);
    return rc;
}

static int score_write_log(const score_bind_t *bind,
                           score_summary_t **summaries,
                           const int64_t *matched_markers,
                           const int64_t *allele_mismatch_markers,
                           int n_samples,
                           char *err, size_t err_sz) {
    FILE *fp;
    int pi;
    if (!bind->log_path || !bind->log_path[0]) return 0;
    fp = fopen(bind->log_path, "w");
    if (!fp) {
        snprintf(err, err_sz, "bcftools_score: cannot write log_path '%s': %s", bind->log_path, strerror(errno));
        return -1;
    }
    fprintf(fp, "# DuckHTS bcftools_score summary\n");
    fprintf(fp, "# mode\t%s\n", bind->gwas_vcf_mode ? "GWAS-VCF" : "TSV");
    fprintf(fp, "# genotype_path\t%s\n", bind->bcf_path ? bind->bcf_path : "");
    fprintf(fp, "# use\t%s\n", score_use_tag_label(bind->use_tag));
    fprintf(fp, "# samples\t%d\n", n_samples);
    fprintf(fp, "# allele_flipping\tnot_performed_by_bcftools_score; use duckhts_munge/bcftools_munge_row outputs to audit REF/ALT swaps before scoring\n");
    fprintf(fp, "summary_name\tsummary_path\tsummary_mode\tmatching_by\tloaded_markers\tall_markers\tmatched_markers\tallele_mismatch_markers\tduplicate_markers\n");
    for (pi = 0; pi < bind->n_prs; pi++) {
        const score_summary_t *summary = summaries[pi];
        fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%lld\t%lld\t%lld\t%lld\n",
                bind->prs_names && bind->prs_names[pi] ? bind->prs_names[pi] : "",
                bind->prs_paths && bind->prs_paths[pi] ? bind->prs_paths[pi] : "",
                bind->gwas_vcf_mode ? "GWAS-VCF" : "TSV",
                (summary && summary->use_snp) ? "marker name" : "chromosome position",
                summary ? summary->n_markers : 0,
                (long long)(summary ? summary->all_markers : 0),
                (long long)(matched_markers ? matched_markers[pi] : 0),
                (long long)(allele_mismatch_markers ? allele_mismatch_markers[pi] : 0),
                (long long)(summary ? summary->duplicate_markers : 0));
    }
    if (fclose(fp) != 0) {
        snprintf(err, err_sz, "bcftools_score: failed closing log_path '%s': %s", bind->log_path, strerror(errno));
        return -1;
    }
    return 0;
}

static void score_bind(duckdb_bind_info info) {
    duckdb_value bcf_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value summary_val = duckdb_bind_get_parameter(info, 1);
    duckdb_logical_type varchar_type;
    duckdb_logical_type bool_type;
    duckdb_logical_type dbl_type;
    duckdb_logical_type bigint_type;
    score_bind_t *bind = (score_bind_t *)duckdb_malloc(sizeof(score_bind_t));

    memset(bind, 0, sizeof(*bind));
    bind->bcf_path = duckdb_get_varchar(bcf_val);
    bind->n_q_thr = 1;
    bind->regions_overlap = -1;
    bind->targets_overlap = -1;
    bind->columns_preset = score_dup("PLINK");

    duckdb_destroy_value(&bcf_val);

    if (!bind->bcf_path || bind->bcf_path[0] == '\0') {
        duckdb_bind_set_error(info, "bcftools_score: bcf_path is required");
        if (summary_val) duckdb_destroy_value(&summary_val);
        score_destroy_bind(bind);
        return;
    }
    {
        char err[512] = {0};
        if (score_add_summary_paths_from_value(bind, summary_val, err, sizeof(err)) != 0) {
            duckdb_bind_set_error(info, err[0] ? err : "bcftools_score: failed to parse summary_path");
            if (summary_val) duckdb_destroy_value(&summary_val);
            score_destroy_bind(bind);
            return;
        }
    }
    if (summary_val) duckdb_destroy_value(&summary_val);

    static const char *const use_names[] = {"GT", "DS", "HDS", "AP", "GP", "AS", NULL};
    static const int use_tags[] = {SCORE_USE_GT, SCORE_USE_DS, SCORE_USE_HDS,
                                   SCORE_USE_AP, SCORE_USE_GP, SCORE_USE_AS};
    char *use_str = NULL;
    char *columns_preset = NULL;
    char *q_str = NULL;
    int read_ok = duckhts_bind_varchar(info, "use", &use_str) &&
                  duckhts_bind_varchar(info, "columns", &columns_preset) &&
                  duckhts_bind_varchar(info, "columns_file", &bind->columns_file) &&
                  duckhts_bind_varchar(info, "q_score_thr", &q_str) &&
                  duckhts_bind_varchar(info, "samples", &bind->samples) &&
                  duckhts_bind_varchar(info, "regions", &bind->regions) &&
                  duckhts_bind_varchar(info, "regions_file", &bind->regions_file) &&
                  duckhts_bind_varchar(info, "targets", &bind->targets) &&
                  duckhts_bind_varchar(info, "targets_file", &bind->targets_file) &&
                  duckhts_bind_varchar(info, "apply_filters", &bind->apply_filters) &&
                  duckhts_bind_varchar(info, "include", &bind->include_expr) &&
                  duckhts_bind_varchar(info, "exclude", &bind->exclude_expr) &&
                  duckhts_bind_varchar(info, "summaries_list_file", &bind->summaries_list_file) &&
                  duckhts_bind_varchar(info, "log_path", &bind->log_path);
    bind->use_variant_id = duckhts_bind_bool(info, "use_variant_id", false);
    bind->counts = duckhts_bind_bool(info, "counts", false);
    bind->force_samples = duckhts_bind_bool(info, "force_samples", false);
    /* INT64_MIN marks an absent overlap mode; the bind keeps -1 for "use the default". */
    int64_t regions_overlap = duckhts_bind_int64(info, "regions_overlap", INT64_MIN);
    int64_t targets_overlap = duckhts_bind_int64(info, "targets_overlap", INT64_MIN);
    const char *error = NULL;
    if (columns_preset) {
        duckdb_free(bind->columns_preset);
        bind->columns_preset = columns_preset;
    }
    if (!read_ok) {
        error = "bcftools_score: out of memory reading parameters";
    } else if (use_str) {
        int i = 0;
        while (use_names[i] && strcasecmp(use_str, use_names[i]) != 0) i++;
        if (use_names[i]) bind->use_tag = use_tags[i];
        else error = "bcftools_score: use must be one of GT,DS,HDS,AP,GP,AS";
    }
    if (!error && q_str && score_parse_q_thresholds(q_str, &bind->q_thr_lp, &bind->n_q_thr) != 0) {
        error = "bcftools_score: q_score_thr must be a comma-separated list of p-values in (0,1]";
    }
    if (!error && regions_overlap != INT64_MIN && (regions_overlap < 0 || regions_overlap > 2)) {
        error = "bcftools_score: regions_overlap must be 0, 1, or 2";
    }
    if (!error && targets_overlap != INT64_MIN && (targets_overlap < 0 || targets_overlap > 2)) {
        error = "bcftools_score: targets_overlap must be 0, 1, or 2";
    }
    if (use_str) duckdb_free(use_str);
    if (q_str) duckdb_free(q_str);
    if (error) {
        duckdb_bind_set_error(info, error);
        score_destroy_bind(bind);
        return;
    }
    if (regions_overlap != INT64_MIN) bind->regions_overlap = (int)regions_overlap;
    if (targets_overlap != INT64_MIN) bind->targets_overlap = (int)targets_overlap;

    if (bind->include_expr && bind->exclude_expr) {
        duckdb_bind_set_error(info, "bcftools_score: only one of include or exclude may be set");
        score_destroy_bind(bind);
        return;
    }
    if (bind->regions && bind->regions_file) {
        duckdb_bind_set_error(info, "bcftools_score: only one of regions or regions_file may be set");
        score_destroy_bind(bind);
        return;
    }
    if (bind->targets && bind->targets_file) {
        duckdb_bind_set_error(info, "bcftools_score: only one of targets or targets_file may be set");
        score_destroy_bind(bind);
        return;
    }
    if (bind->regions_file && bind->regions_file[0]) {
        char err[512] = {0};
        char *parsed = score_regions_from_file(bind->regions_file, err, sizeof(err));
        if (!parsed) {
            duckdb_bind_set_error(info, err[0] ? err : "bcftools_score: failed to parse regions_file");
            score_destroy_bind(bind);
            return;
        }
        bind->regions = parsed;
    }
    if (bind->targets_file && bind->targets_file[0]) {
        char err[512] = {0};
        char *parsed = score_regions_from_file(bind->targets_file, err, sizeof(err));
        if (!parsed) {
            duckdb_bind_set_error(info, err[0] ? err : "bcftools_score: failed to parse targets_file");
            score_destroy_bind(bind);
            return;
        }
        bind->targets = parsed;
    }
    if (bind->summaries_list_file && bind->summaries_list_file[0]) {
        char err[512] = {0};
        if (bind->n_summary_paths > 0) {
            duckdb_bind_set_error(info, "bcftools_score: provide either summary_path(s) or summaries_list_file, not both");
            score_destroy_bind(bind);
            return;
        }
        if (score_add_summary_paths_from_list_file(bind, bind->summaries_list_file, err, sizeof(err)) != 0) {
            duckdb_bind_set_error(info, err[0] ? err : "bcftools_score: failed to parse summaries_list_file");
            score_destroy_bind(bind);
            return;
        }
    }
    if (bind->n_summary_paths == 0) {
        duckdb_bind_set_error(info, "bcftools_score: summary_path(s) or summaries_list_file are required");
        score_destroy_bind(bind);
        return;
    }
    /* Detect GWAS-VCF mode if all summary files are VCF/BCF and no columns_file
     * was provided. Otherwise use TSV mode, where each summary file is one PRS,
     * matching upstream bcftools +score -c/--columns behavior. */
    {
        int si;
        int n_vcf = 0;
        for (si = 0; si < bind->n_summary_paths; si++) {
            if (score_is_vcf_summary(bind->summary_paths[si])) n_vcf++;
        }
        if (n_vcf > 0 && n_vcf != bind->n_summary_paths) {
            duckdb_bind_set_error(info, "bcftools_score: cannot mix GWAS-VCF and TSV summary files in one call");
            score_destroy_bind(bind);
            return;
        }
        bind->gwas_vcf_mode = (n_vcf == bind->n_summary_paths && !bind->columns_file) ? 1 : 0;
    }

    if (bind->gwas_vcf_mode) {
        int si;
        int offset = 0;
        bind->summary_prs_counts = (int *)duckdb_malloc(sizeof(int) * (size_t)bind->n_summary_paths);
        if (!bind->summary_prs_counts) {
            duckdb_bind_set_error(info, "bcftools_score: out of memory");
            score_destroy_bind(bind);
            return;
        }
        for (si = 0; si < bind->n_summary_paths; si++) {
            htsFile *sfp = hts_open(bind->summary_paths[si], "r");
            bcf_hdr_t *shdr;
            int n_file_prs;
            int pi;
            if (!sfp) {
                duckdb_bind_set_error(info, "bcftools_score: cannot open GWAS-VCF summary for header probe");
                score_destroy_bind(bind);
                return;
            }
            shdr = bcf_hdr_read(sfp);
            if (!shdr) {
                hts_close(sfp);
                duckdb_bind_set_error(info, "bcftools_score: cannot read GWAS-VCF header");
                score_destroy_bind(bind);
                return;
            }
            {
                int es_id = bcf_hdr_id2int(shdr, BCF_DT_ID, "ES");
                if (!bcf_hdr_idinfo_exists(shdr, BCF_HL_FMT, es_id)) {
                    bcf_hdr_destroy(shdr);
                    hts_close(sfp);
                    duckdb_bind_set_error(info, "bcftools_score: GWAS-VCF summary does not include the ES FORMAT field");
                    score_destroy_bind(bind);
                    return;
                }
            }
            if (bind->q_thr_lp) {
                int lp_id = bcf_hdr_id2int(shdr, BCF_DT_ID, "LP");
                if (!bcf_hdr_idinfo_exists(shdr, BCF_HL_FMT, lp_id)) {
                    bcf_hdr_destroy(shdr);
                    hts_close(sfp);
                    duckdb_bind_set_error(info, "bcftools_score: GWAS-VCF summary does not include the LP FORMAT field (required by q_score_thr)");
                    score_destroy_bind(bind);
                    return;
                }
            }
            n_file_prs = bcf_hdr_nsamples(shdr);
            if (n_file_prs == 0) {
                bcf_hdr_destroy(shdr);
                hts_close(sfp);
                duckdb_bind_set_error(info, "bcftools_score: GWAS-VCF summary has no samples");
                score_destroy_bind(bind);
                return;
            }
            {
                int old_n_prs = bind->n_prs;
                int new_n_prs = old_n_prs + n_file_prs;
                char **new_names = (char **)realloc(bind->prs_names, sizeof(char *) * (size_t)new_n_prs);
                char **new_paths;
                if (!new_names) {
                    bcf_hdr_destroy(shdr);
                    hts_close(sfp);
                    duckdb_bind_set_error(info, "bcftools_score: out of memory");
                    score_destroy_bind(bind);
                    return;
                }
                bind->prs_names = new_names;
                new_paths = (char **)realloc(bind->prs_paths, sizeof(char *) * (size_t)new_n_prs);
                if (!new_paths) {
                    bcf_hdr_destroy(shdr);
                    hts_close(sfp);
                    duckdb_bind_set_error(info, "bcftools_score: out of memory");
                    score_destroy_bind(bind);
                    return;
                }
                bind->prs_paths = new_paths;
                memset(bind->prs_names + old_n_prs, 0, sizeof(char *) * (size_t)n_file_prs);
                memset(bind->prs_paths + old_n_prs, 0, sizeof(char *) * (size_t)n_file_prs);
                bind->summary_prs_counts[si] = n_file_prs;
                bind->n_prs = new_n_prs;
                offset = old_n_prs;
            }
            for (pi = 0; pi < n_file_prs; pi++) {
                bind->prs_names[offset + pi] = score_dup(shdr->samples[pi]);
                bind->prs_paths[offset + pi] = score_dup(bind->summary_paths[si]);
                if (!bind->prs_names[offset + pi] || !bind->prs_paths[offset + pi]) {
                    bcf_hdr_destroy(shdr);
                    hts_close(sfp);
                    duckdb_bind_set_error(info, "bcftools_score: out of memory");
                    score_destroy_bind(bind);
                    return;
                }
            }
            bcf_hdr_destroy(shdr);
            hts_close(sfp);
        }
    } else {
        int si;
        bind->n_prs = bind->n_summary_paths;
        bind->prs_names = (char **)duckdb_malloc(sizeof(char *) * (size_t)bind->n_prs);
        bind->prs_paths = (char **)duckdb_malloc(sizeof(char *) * (size_t)bind->n_prs);
        bind->summary_prs_counts = (int *)duckdb_malloc(sizeof(int) * (size_t)bind->n_summary_paths);
        if (!bind->prs_names || !bind->prs_paths || !bind->summary_prs_counts) {
            duckdb_bind_set_error(info, "bcftools_score: out of memory");
            score_destroy_bind(bind);
            return;
        }
        memset(bind->prs_names, 0, sizeof(char *) * (size_t)bind->n_prs);
        memset(bind->prs_paths, 0, sizeof(char *) * (size_t)bind->n_prs);
        for (si = 0; si < bind->n_summary_paths; si++) {
            bind->summary_prs_counts[si] = 1;
            bind->prs_names[si] = score_prs_name_from_path(bind->summary_paths[si]);
            bind->prs_paths[si] = score_dup(bind->summary_paths[si]);
            if (!bind->prs_names[si] || !bind->prs_paths[si]) {
                duckdb_bind_set_error(info, "bcftools_score: out of memory");
                score_destroy_bind(bind);
                return;
            }
        }
    }

    {
        char err[512] = {0};
        if (score_validate_prs_names(bind, err, sizeof(err)) != 0) {
            duckdb_bind_set_error(info, err[0] ? err : "bcftools_score: invalid PRS output names");
            score_destroy_bind(bind);
            return;
        }
    }

    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    dbl_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "SAMPLE", varchar_type);
    {
        int pi, j;
        for (pi = 0; pi < bind->n_prs; pi++) {
            const char *pname = bind->prs_names[pi];
            for (j = 0; j < bind->n_q_thr; j++) {
                char name[256];
                score_format_metric_name(bind, pname, j, 0, name, sizeof(name));
                duckdb_bind_add_result_column(info, name, dbl_type);
                if (bind->counts) {
                    char cnt_name[256];
                    score_format_metric_name(bind, pname, j, 1, cnt_name, sizeof(cnt_name));
                    duckdb_bind_add_result_column(info, cnt_name, bigint_type);
                }
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
    score_scan_source_t src = {0};
    int regions_parsed = 0;
    bcf_hdr_t *hdr = NULL;
    score_summary_t **summaries = NULL;    /* one summary index per PRS */
    int64_t *matched_markers = NULL;
    int64_t *allele_mismatch_markers = NULL;
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
    int n_prs;
    int i;
    int scan_failed = 0;
    char err[512];

    memset(init, 0, sizeof(*init));
    if (score_init_source(bind, &src, &hdr, err, sizeof(err)) != 0) {
        score_source_destroy(&src);
        duckdb_init_set_error(info, err[0] ? err : "bcftools_score: failed to initialize input source");
        duckdb_free(init);
        return;
    }

    /* Subset samples if requested (port of upstream score.c:652-667) */
    if (bind->samples && bind->samples[0]) {
        int ret = bcf_hdr_set_samples(hdr, bind->samples, 0);
        if (ret < 0) {
            score_source_destroy(&src);
            duckdb_init_set_error(info, "bcftools_score: error parsing the sample list");
            duckdb_free(init);
            return;
        } else if (ret > 0) {
            if (!bind->force_samples) {
                snprintf(err, sizeof(err), "bcftools_score: sample #%d not found in the header; use force_samples := true to ignore", ret);
                score_source_destroy(&src);
                duckdb_init_set_error(info, err);
                duckdb_free(init);
                return;
            }
        }
        if (bcf_hdr_nsamples(hdr) == 0) {
            score_source_destroy(&src);
            duckdb_init_set_error(info, "bcftools_score: subsetting has removed all samples");
            duckdb_free(init);
            return;
        }
    }

    if (score_resolve_use_tag(bind, hdr, err, sizeof(err)) != 0) {
        score_source_destroy(&src);
        duckdb_init_set_error(info, err);
        duckdb_free(init);
        return;
    }

    /* Load summary data: TSV mode has one PRS per summary file; GWAS-VCF mode
     * has one PRS per FORMAT sample across all summary VCF/BCF files. */
    n_prs = bind->n_prs;
    summaries = (score_summary_t **)duckdb_malloc(sizeof(score_summary_t *) * (size_t)n_prs);
    if (!summaries) {
        score_source_destroy(&src);
        duckdb_init_set_error(info, "bcftools_score: out of memory");
        duckdb_free(init);
        return;
    }
    memset(summaries, 0, sizeof(score_summary_t *) * (size_t)n_prs);
    if (bind->gwas_vcf_mode) {
        int si;
        int offset = 0;
        for (si = 0; si < bind->n_summary_paths; si++) {
            int file_n = bind->summary_prs_counts[si];
            score_summary_t **tmp = score_load_summary_vcf(bind, hdr, bind->summary_paths[si], file_n, err, sizeof(err));
            int pi;
            if (!tmp) {
                score_source_destroy(&src);
                for (i = 0; i < n_prs; i++) if (summaries[i]) score_destroy_summary(summaries[i]);
                duckdb_free(summaries);
                duckdb_init_set_error(info, err[0] ? err : "bcftools_score: failed to load GWAS-VCF summary");
                duckdb_free(init);
                return;
            }
            for (pi = 0; pi < file_n; pi++) summaries[offset + pi] = tmp[pi];
            duckdb_free(tmp);
            offset += file_n;
        }
    } else {
        int pi;
        for (pi = 0; pi < n_prs; pi++) {
            summaries[pi] = score_load_summary(bind, hdr, bind->summary_paths[pi], err, sizeof(err));
            if (!summaries[pi]) {
                score_source_destroy(&src);
                for (i = 0; i < n_prs; i++) if (summaries[i]) score_destroy_summary(summaries[i]);
                duckdb_free(summaries);
                duckdb_init_set_error(info, err[0] ? err : "bcftools_score: failed to load summary file");
                duckdb_free(init);
                return;
            }
        }
    }

    n_samples = bcf_hdr_nsamples(hdr);
    n_metrics = n_prs * bind->n_q_thr * (bind->counts ? 2 : 1);
    init->n_samples = n_samples;
    init->n_metrics = n_metrics;
    init->sample_names = (char **)duckdb_malloc(sizeof(char *) * (size_t)n_samples);
    init->metric_names = (char **)duckdb_malloc(sizeof(char *) * (size_t)n_metrics);
    init->metric_values = (double **)duckdb_malloc(sizeof(double *) * (size_t)n_metrics);
    init->metric_is_count = (uint8_t *)duckdb_malloc((size_t)n_metrics);
    if (!init->sample_names || !init->metric_names || !init->metric_values || !init->metric_is_count) {
        if (summaries) {
            for (i = 0; i < n_prs; i++) if (summaries[i]) score_destroy_summary(summaries[i]);
            duckdb_free(summaries);
        }
        score_source_destroy(&src);
        score_destroy_init(init);
        duckdb_init_set_error(info, "bcftools_score: out of memory");
        return;
    }

    for (i = 0; i < n_samples; i++) {
        init->sample_names[i] = score_dup(hdr->samples[i]);
    }
    {
        int pi, j;
        int metric_idx = 0;
        memset(init->metric_is_count, 0, (size_t)n_metrics);
        for (pi = 0; pi < n_prs; pi++) {
            const char *pname = bind->prs_names[pi];
            for (j = 0; j < bind->n_q_thr; j++) {
                char name[256];
                score_format_metric_name(bind, pname, j, 0, name, sizeof(name));
                init->metric_names[metric_idx] = score_dup(name);
                init->metric_values[metric_idx] = (double *)duckdb_malloc(sizeof(double) * (size_t)n_samples);
                memset(init->metric_values[metric_idx], 0, sizeof(double) * (size_t)n_samples);
                init->metric_is_count[metric_idx] = 0;
                metric_idx++;
                if (bind->counts) {
                    char cnt_name[256];
                    score_format_metric_name(bind, pname, j, 1, cnt_name, sizeof(cnt_name));
                    init->metric_names[metric_idx] = score_dup(cnt_name);
                    init->metric_values[metric_idx] = (double *)duckdb_malloc(sizeof(double) * (size_t)n_samples);
                    memset(init->metric_values[metric_idx], 0, sizeof(double) * (size_t)n_samples);
                    init->metric_is_count[metric_idx] = 1;
                    metric_idx++;
                }
            }
        }
    }

    (void)unused_alleles;
    missing = (float *)duckdb_malloc(sizeof(float) * (size_t)n_samples);
    matched_markers = (int64_t *)duckdb_malloc(sizeof(int64_t) * (size_t)n_prs);
    allele_mismatch_markers = (int64_t *)duckdb_malloc(sizeof(int64_t) * (size_t)n_prs);
    if (!missing || !matched_markers || !allele_mismatch_markers) {
        if (missing) duckdb_free(missing);
        if (matched_markers) duckdb_free(matched_markers);
        if (allele_mismatch_markers) duckdb_free(allele_mismatch_markers);
        score_source_destroy(&src);
        for (i = 0; i < n_prs; i++) if (summaries[i]) score_destroy_summary(summaries[i]);
        duckdb_free(summaries);
        score_destroy_init(init);
        duckdb_init_set_error(info, "bcftools_score: out of memory");
        return;
    }
    memset(matched_markers, 0, sizeof(int64_t) * (size_t)n_prs);
    memset(allele_mismatch_markers, 0, sizeof(int64_t) * (size_t)n_prs);

    /* Per-PRS marker indices: TSV mode has one summary per file; GWAS-VCF mode
     * has one preloaded summary per FORMAT sample across summary VCF files. */
    {
        int *idxs = (int *)duckdb_malloc(sizeof(int) * (size_t)n_prs);
        if (!idxs) {
            if (matched_markers) duckdb_free(matched_markers);
            if (allele_mismatch_markers) duckdb_free(allele_mismatch_markers);
            if (missing) duckdb_free(missing);
            score_source_destroy(&src);
            for (i = 0; i < n_prs; i++) if (summaries[i]) score_destroy_summary(summaries[i]);
            duckdb_free(summaries);
            score_destroy_init(init);
            duckdb_init_set_error(info, "bcftools_score: out of memory");
            return;
        }

        bcf1_t *rec = NULL;
        while (score_source_next(&src, &rec)) {
            int number;
            int k;
            int skip_line = 1;
            bcf_unpack(rec, BCF_UN_ALL);

            if (!score_record_passes_filter_list(hdr, rec, bind->apply_filters)) continue;
            /* Parse on the first FILTER-passing record; empty or fully FILTER-
             * excluded inputs do not evaluate region expressions. */
            if (!regions_parsed) {
                if (score_parse_region_list(hdr, bind->regions, &src.regions) < 0 ||
                    score_parse_region_list(hdr, bind->targets, &src.targets) < 0) {
                    snprintf(err, sizeof(err), "bcftools_score: invalid region/target expression");
                    goto score_record_error;
                }
                regions_parsed = 1;
            }
            int in_regions = score_record_in_regions(rec, &src.regions,
                                                     bind->regions_overlap >= 0 ? bind->regions_overlap : 1);
            int in_targets = score_record_in_regions(rec, &src.targets,
                                                     bind->targets_overlap >= 0 ? bind->targets_overlap : 0);
            if (!in_regions || !in_targets) continue;
            {
                int expr_keep = score_record_matches_exprs(rec, &src, err, sizeof(err));
                if (expr_keep < 0) {
                    goto score_record_error;
                }
                if (expr_keep == 0) continue;
            }
            /* Phase 1: Look up marker index in each PRS summary. */
            {
                int pi;
                for (pi = 0; pi < n_prs; pi++) {
                    idxs[pi] = -1;
                    if (summaries[pi]->use_snp) {
                        if (rec->d.id && rec->d.id[0] != '\0') {
                            khash_str2int_get(summaries[pi]->id2idx, rec->d.id, &idxs[pi]);
                        }
                    } else {
                        khash_t(score64) *hash = (khash_t(score64) *)summaries[pi]->rid_pos2idx;
                        uint64_t key = (((uint64_t)(uint32_t)rec->rid) << 44) | (uint64_t)rec->pos;
                        khiter_t it = kh_get(score64, hash, key);
                        if (it != kh_end(hash)) idxs[pi] = kh_val(hash, it);
                    }
                    if (idxs[pi] >= 0 && idxs[pi] < summaries[pi]->n_markers) skip_line = 0;
                }
            }
            if (skip_line) continue;

            /* Phase 2: Extract genotype/dosage into per-allele-per-sample array */
            memset(missing, 0, sizeof(float) * (size_t)n_samples);
            /* Keep allele/sample offsets representable in the int indexing below. */
            uint64_t n_aps = (uint64_t)rec->n_allele * (uint64_t)n_samples;
            if (n_aps > INT_MAX || n_aps > SIZE_MAX / sizeof(float)) {
                snprintf(err, sizeof(err), "bcftools_score: allele/sample array exceeds allocation or indexing limits");
                goto score_record_error;
            }
            if (m_aps < (int)n_aps) {
                float *new_aps = (float *)realloc(aps, sizeof(float) * (size_t)n_aps);
                if (!new_aps) {
                    snprintf(err, sizeof(err), "bcftools_score: out of memory");
                    goto score_record_error;
                }
                aps = new_aps;
                m_aps = (int)n_aps;
            }
            memset(aps, 0, sizeof(float) * (size_t)n_aps);

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
                if ((int64_t)number != (int64_t)n_samples * (rec->n_allele - 1)) {
                    snprintf(err, sizeof(err), "bcftools_score: invalid DS FORMAT cardinality at %s:%" PRId64,
                             bcf_hdr_id2name(hdr, rec->rid), (int64_t)rec->pos + 1);
                    goto score_record_error;
                }
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
                    if ((int64_t)number != (int64_t)n_samples * (rec->n_allele - 1)) {
                        snprintf(err, sizeof(err), "bcftools_score: invalid %s FORMAT cardinality at %s:%" PRId64,
                                 ap_tags[ap], bcf_hdr_id2name(hdr, rec->rid), (int64_t)rec->pos + 1);
                        goto score_record_error;
                    }
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
                if ((int64_t)number != (int64_t)n_samples * rec->n_allele * (rec->n_allele + 1) / 2) {
                    snprintf(err, sizeof(err), "bcftools_score: invalid GP FORMAT cardinality at %s:%" PRId64,
                             bcf_hdr_id2name(hdr, rec->rid), (int64_t)rec->pos + 1);
                    goto score_record_error;
                }
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
                                int idx = (int)((int64_t)b * (b + 1) / 2 + a);
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

            /* Phase 3: Accumulate scores for each PRS. */
            {
                int pi;
                for (pi = 0; pi < n_prs; pi++) {
                    const score_marker_t *marker;
                    const char *a1;
                    float es, lp;
                    int idx_allele, j;

                    if (idxs[pi] < 0 || idxs[pi] >= summaries[pi]->n_markers) continue;
                    marker = &summaries[pi]->markers[idxs[pi]];
                    a1 = marker->a1;
                    if (!a1) continue;

                    for (idx_allele = 0; idx_allele < rec->n_allele; idx_allele++) {
                        if (strcmp(a1, rec->d.allele[idx_allele]) == 0) break;
                    }
                    if (idx_allele == rec->n_allele) {
                        allele_mismatch_markers[pi]++;
                        continue;
                    }

                    matched_markers[pi]++;
                    es = marker->es;
                    lp = marker->lp;
                    for (j = 0; j < bind->n_q_thr; j++) {
                        int out_idx = (pi * bind->n_q_thr + j) * (bind->counts ? 2 : 1);
                        if (bind->q_thr_lp && !isnan(lp) && lp < bind->q_thr_lp[j]) continue;
                        for (k = 0; k < n_samples; k++) {
                            if (missing[k]) continue;
                            {
                                /* Upstream bcftools +score accumulates scores in float arrays
                                 * and writes %#.6g text output. Preserve the float32 summation
                                 * path before widening to DuckDB DOUBLE so large synthetic PRS
                                 * benchmarks do not drift from the plugin by accumulated double
                                 * precision differences. */
                                float score_value = (float)init->metric_values[out_idx][k];
                                score_value += es * aps[idx_allele * n_samples + k];
                                init->metric_values[out_idx][k] = (double)score_value;
                            }
                            if (bind->counts) init->metric_values[out_idx + 1][k] += 1.0;
                        }
                    }
                }
            }
        }

        goto score_records_done;
score_record_error:
        scan_failed = 1;
score_records_done:
        if (idxs) duckdb_free(idxs);
    }

    if (scan_failed) goto score_cleanup;

    if (bind->log_path && bind->log_path[0]) {
        if (score_write_log(bind, summaries, matched_markers, allele_mismatch_markers, n_samples, err, sizeof(err)) != 0) {
            if (int32_arr) free(int32_arr);
            if (float_arr) free(float_arr);
            if (aps) duckdb_free(aps);
            if (missing) duckdb_free(missing);
            if (matched_markers) duckdb_free(matched_markers);
            if (allele_mismatch_markers) duckdb_free(allele_mismatch_markers);
            score_source_destroy(&src);
            if (summaries) {
                for (i = 0; i < n_prs; i++) {
                    if (summaries[i]) score_destroy_summary(summaries[i]);
                }
                duckdb_free(summaries);
            }
            score_destroy_init(init);
            duckdb_init_set_error(info, err[0] ? err : "bcftools_score: failed to write log_path");
            return;
        }
    }

score_cleanup:
    if (int32_arr) free(int32_arr);
    if (float_arr) free(float_arr);
    if (aps) duckdb_free(aps);
    if (missing) duckdb_free(missing);
    if (matched_markers) duckdb_free(matched_markers);
    if (allele_mismatch_markers) duckdb_free(allele_mismatch_markers);
    score_source_destroy(&src);
    if (summaries) {
        for (i = 0; i < n_prs; i++) {
            if (summaries[i]) score_destroy_summary(summaries[i]);
        }
        duckdb_free(summaries);
    }
    hdr = NULL;
    if (scan_failed) {
        score_destroy_init(init);
        duckdb_init_set_error(info, err);
        return;
    }

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
    duckdb_logical_type any_type = duckdb_create_logical_type(DUCKDB_TYPE_ANY);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);

    duckdb_table_function_set_name(tf, "bcftools_score");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, any_type);
    duckdb_table_function_add_named_parameter(tf, "use", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "columns", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "columns_file", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "q_score_thr", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "summaries_list_file", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "log_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "use_variant_id", bool_type);
    duckdb_table_function_add_named_parameter(tf, "counts", bool_type);
    duckdb_table_function_add_named_parameter(tf, "samples", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "force_samples", bool_type);
    duckdb_table_function_add_named_parameter(tf, "regions", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "regions_file", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "regions_overlap", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "targets", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "targets_file", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "targets_overlap", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "apply_filters", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "include", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "exclude", varchar_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&any_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&bigint_type);

    duckdb_table_function_set_bind(tf, score_bind);
    duckdb_table_function_set_init(tf, score_init);
    duckdb_table_function_set_function(tf, score_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
