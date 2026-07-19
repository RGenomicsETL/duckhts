/**
 * DuckHTS FASTA/FASTQ Reader
 *
 * Table functions for reading FASTA and FASTQ files over htslib transport.
 * htslib's sam_open() auto-detects fasta/fastq format (including .gz).
 * FASTA keeps the htslib SAM-sequence adapter; FASTQ is parsed directly into
 * projected DuckDB columns so text records do not make a temporary bam1_t
 * round trip.  The FASTA adapter gives us access to:
 *   - bam_get_qname(b)       → sequence name
 *   - bam_get_seq(b)         → packed sequence (decode with seq_nt16_str[bam_seqi()])
 *   - bam_get_qual(b)        → quality values (Phred, add 33 for ASCII)
 *   - b->core.l_qseq         → sequence length
 *
 * API reference: htslib-1.23 samples/read_fast.c
 *
 * Schema:
 *   read_fasta(path) → (NAME VARCHAR, DESCRIPTION VARCHAR, SEQUENCE VARCHAR)
 *   read_fastq(path) → (NAME VARCHAR, DESCRIPTION VARCHAR, SEQUENCE VARCHAR, QUALITY VARCHAR)
 *
 * read_fasta named parameters:
 *   - region := 'chr:start-end[,chr2:...]' for indexed faidx fetches
 *   - index_path := explicit .fai path
 *   - gzi_path := explicit BGZF .gzi path for bgzipped FASTA
 *
 * Note: htslib's default FASTA/FASTQ adapter does not expose the header
 * comment as a "CO" tag unless its optional aux parser is enabled.  Keep
 * DESCRIPTION null on the direct FASTQ path to preserve read_fastq output.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>
#include <limits.h>
#include <strings.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>

#include "include/hts_io_tuning.h"
#include "include/seq_encoding.h"
#include "include/quality_encoding.h"

/* ================================================================
 * Column indices
 * ================================================================ */

enum {
    SEQ_COL_NAME = 0,
    SEQ_COL_DESCRIPTION,
    SEQ_COL_SEQUENCE,
    SEQ_COL_QUALITY,  /* FASTQ only */
    SEQ_COL_MATE,     /* FASTQ paired/interleaved only */
    SEQ_COL_PAIR_ID,
    SEQ_COL_MAX
};

/* ================================================================
 * Bind Data
 * ================================================================ */

typedef struct {
    char *file_path;
    char *mate_path;
    char *index_path;
    char *gzi_path;
    char *region;
    char **regions;
    unsigned int n_regions;
    int is_fastq;
    int interleaved;
    int paired;
    int scan_sequential;
    int direct_fastq;    /* input was detected as FASTQ; use the text parser */
    int seq_nt16;       /* 1 if sequence_encoding := 'nt16'; SEQUENCE → LIST(UTINYINT) */
    duckhts_quality_representation qual_repr;
    duckhts_quality_encoding input_quality_encoding;
} seq_bind_data_t;

typedef struct {
    kstring_t line;
    kstring_t name;
    kstring_t sequence;
    kstring_t quality;
    size_t sequence_len;
} fastq_record_t;

/* ================================================================
 * Init Data
 * ================================================================ */

typedef struct {
    samFile *fp;
    sam_hdr_t *hdr;
    bam1_t *rec;
    htsFile *count_fp;
    kstring_t count_line;
    samFile *fp_mate;
    sam_hdr_t *hdr_mate;
    bam1_t *rec_mate;
    htsFile *count_fp_mate;
    fastq_record_t fastq_record;
    fastq_record_t fastq_record_mate;
    int done;
    int is_fastq;
    int interleaved;
    int paired;
    int count_only;
    int direct_fastq;
    int fasta_count_only;
    uint64_t count_remaining;
    int pending_mate;
    int interleaved_mate;
    faidx_t *fai;
    unsigned int n_regions;
    unsigned int next_region_idx;
    char **regions;  /* reference to bind regions */

    idx_t column_count;
    idx_t *column_ids;

    /* Reusable buffer for decoded sequence */
    char *seq_buf;
    size_t seq_buf_cap;
    /* Reusable buffer for decoded quality */
    char *qual_buf;
    size_t qual_buf_cap;
    char *pair_buf;
    size_t pair_buf_cap;
    int seq_nt16;       /* 1 if sequence_encoding := 'nt16' */
    duckhts_quality_representation qual_repr;
    duckhts_quality_encoding input_quality_encoding;
    int fastq_need_name;
    int fastq_need_sequence;
    int fastq_need_quality;
} seq_init_data_t;

/* ================================================================
 * Destructors
 * ================================================================ */

static void destroy_seq_bind(void *data) {
    seq_bind_data_t *b = (seq_bind_data_t *)data;
    if (!b) return;
    if (b->file_path) duckdb_free(b->file_path);
    if (b->mate_path) duckdb_free(b->mate_path);
    if (b->index_path) duckdb_free(b->index_path);
    if (b->gzi_path) duckdb_free(b->gzi_path);
    if (b->region) duckdb_free(b->region);
    if (b->regions) {
        for (unsigned int i = 0; i < b->n_regions; i++) {
            if (b->regions[i]) duckdb_free(b->regions[i]);
        }
        duckdb_free(b->regions);
    }
    duckdb_free(b);
}

static void fastq_record_destroy(fastq_record_t *record) {
    if (!record) return;
    ks_free(&record->line);
    ks_free(&record->name);
    ks_free(&record->sequence);
    ks_free(&record->quality);
}

static void destroy_seq_init(void *data) {
    seq_init_data_t *init = (seq_init_data_t *)data;
    if (!init) return;
    if (init->rec) bam_destroy1(init->rec);
    if (init->hdr) sam_hdr_destroy(init->hdr);
    if (init->fp) sam_close(init->fp);
    if (init->count_fp) hts_close(init->count_fp);
    if (init->fai) fai_destroy(init->fai);
    if (init->rec_mate) bam_destroy1(init->rec_mate);
    if (init->hdr_mate) sam_hdr_destroy(init->hdr_mate);
    if (init->fp_mate) sam_close(init->fp_mate);
    if (init->count_fp_mate) hts_close(init->count_fp_mate);
    ks_free(&init->count_line);
    fastq_record_destroy(&init->fastq_record);
    fastq_record_destroy(&init->fastq_record_mate);
    if (init->column_ids) duckdb_free(init->column_ids);
    if (init->seq_buf) free(init->seq_buf);
    if (init->qual_buf) free(init->qual_buf);
    if (init->pair_buf) free(init->pair_buf);
    duckdb_free(init);
}

/* ================================================================
 * Helpers
 * ================================================================ */

static inline void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *v = duckdb_vector_get_validity(vec);
    v[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static int ensure_buf(char **buf, size_t *cap, size_t need) {
    if (need == SIZE_MAX) return -1;
    size_t n = need + 1;
    if (n > *cap) {
        size_t new_cap = n <= SIZE_MAX / 2 ? n * 2 : n;
        char *new_buf = (char *)realloc(*buf, new_cap);
        if (!new_buf) return -1;
        *buf = new_buf;
        *cap = new_cap;
    }
    return 0;
}

/* Decode packed sequence via seq_decode_to_string() from include/seq_encoding.h */

/* Decode quality to canonical Phred+33 ASCII */
static void decode_qual(const uint8_t *qual, int len, char *buf,
                        duckhts_quality_encoding input_quality_encoding) {
    for (int i = 0; i < len; i++)
        buf[i] = (char)(duckhts_normalize_fastq_quality(qual[i], input_quality_encoding) + 33);
    buf[len] = '\0';
}

static const char *strip_pair_suffix(const char *name, char *buf, size_t buf_cap) {
    if (!name) return "";
    size_t len = strlen(name);
    if (len >= 2 && name[len - 2] == '/' &&
        (name[len - 1] == '1' || name[len - 1] == '2')) {
        len -= 2;
    }
    if (len + 1 > buf_cap) len = buf_cap - 1;
    memcpy(buf, name, len);
    buf[len] = '\0';
    return buf;
}

static char *strdup_duckdb(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char *copy = (char *)duckdb_malloc(len);
    if (copy) memcpy(copy, s, len);
    return copy;
}

static char *default_fasta_index_path(const char *path) {
    size_t path_len;
    char *out;

    if (!path) return NULL;
    path_len = strlen(path);
    out = (char *)duckdb_malloc(path_len + 5);
    if (!out) return NULL;
    memcpy(out, path, path_len);
    memcpy(out + path_len, ".fai", 5);
    return out;
}

static void parse_regions_duckdb(const char *region_str, char ***out_regions, unsigned int *out_count) {
    *out_regions = NULL;
    *out_count = 0;
    if (!region_str || region_str[0] == '\0') return;

    unsigned int count = 1;
    for (const char *p = region_str; *p; p++) if (*p == ',') count++;

    char **arr = (char **)duckdb_malloc(sizeof(char *) * count);
    char *dup = strdup_duckdb(region_str);
    if (!arr || !dup) {
        if (arr) duckdb_free(arr);
        if (dup) duckdb_free(dup);
        return;
    }

    unsigned int idx = 0;
    char *tok = strtok(dup, ",");
    while (tok && idx < count) {
        while (*tok == ' ' || *tok == '\t') tok++;
        size_t len = strlen(tok);
        while (len > 0 && (tok[len - 1] == ' ' || tok[len - 1] == '\t')) tok[--len] = '\0';
        if (len > 0) {
            arr[idx] = strdup_duckdb(tok);
            if (!arr[idx]) {
                for (unsigned int i = 0; i < idx; i++) duckdb_free(arr[i]);
                duckdb_free(arr);
                duckdb_free(dup);
                return;
            }
            idx++;
        }
        tok = strtok(NULL, ",");
    }
    duckdb_free(dup);
    *out_regions = arr;
    *out_count = idx;
}

static int seq_parse_scan_mode(const char *mode, int *scan_sequential) {
    if (!scan_sequential) return 0;
    *scan_sequential = 0;
    if (!mode || mode[0] == '\0' || strcasecmp(mode, "auto") == 0) {
        return 1;
    }
    if (strcasecmp(mode, "sequential") == 0 ||
        strcasecmp(mode, "streaming") == 0 ||
        strcasecmp(mode, "stream") == 0 ||
        strcasecmp(mode, "seq") == 0) {
        *scan_sequential = 1;
        return 1;
    }
    return 0;
}

static int fastq_line_is_header(const kstring_t *line) {
    return line && line->l > 0 && line->s && line->s[0] == '@';
}

static int fastq_line_is_plus(const kstring_t *line) {
    return line && line->l > 0 && line->s && line->s[0] == '+';
}

static int fastq_extract_qname(const kstring_t *line, kstring_t *name,
                               int keep_name) {
    size_t start = 1;
    size_t end = start;
    size_t name_len;

    if (!line) return -1;
    /* The BAM QNAME limit is part of the published NAME/PAIR_ID contract,
     * not FASTQ record framing.  Avoid both the scan and the limit when the
     * query name is neither returned nor needed to compare mate files. */
    if (!keep_name) return 0;
    while (end < line->l && !isspace((unsigned char)line->s[end])) end++;
    name_len = end - start;

    /* Match htslib fastq_parse1(): strip a terminal /<digit> before
     * publishing the query name. */
    if (name_len >= 2 && line->s[end - 2] == '/' &&
        isdigit((unsigned char)line->s[end - 1])) {
        name_len -= 2;
    }
    if (name_len == 0) name_len = 1; /* bam_set1() publishes "*" */
    if (name_len > BAM_MAX_QNAME_LEN) return -2;

    if (ks_resize(name, name_len + 1) < 0) return -1;
    if (end == start) {
        name->s[0] = '*';
    } else {
        memcpy(name->s, line->s + start, name_len);
    }
    name->l = name_len;
    name->s[name->l] = '\0';
    return 0;
}

static int fastq_append_line(kstring_t *dst, const kstring_t *line,
                             int keep, const char *path, const char *field,
                             char *errbuf, size_t errbuf_len) {
    if (!keep) return 0;
    if (kputsn(line->s, line->l, dst) < 0) {
        snprintf(errbuf, errbuf_len,
                 "read_fastq: out of memory parsing %s in %s",
                 field, path ? path : "");
        return -1;
    }
    return 0;
}

/*
 * Read one FASTQ record over htslib transport.  This follows the framing and
 * name rules in htslib fastq_parse1() but retains only projected text fields.
 * Returns 1 for a record, 0 at clean EOF, and -1 with errbuf populated.
 */
static int fastq_read_record(htsFile *fp, fastq_record_t *record,
                             int keep_name, int keep_sequence, int keep_quality,
                             const char *path, char *errbuf, size_t errbuf_len) {
    int ret;
    size_t seq_len = 0;
    size_t remaining;

    record->sequence_len = 0;
    record->sequence.l = 0;
    record->quality.l = 0;
    if (record->sequence.s) record->sequence.s[0] = '\0';
    if (record->quality.s) record->quality.s[0] = '\0';

    ret = hts_getline(fp, KS_SEP_LINE, &record->line);
    if (ret == -1) return 0;
    if (ret < -1) {
        snprintf(errbuf, errbuf_len,
                 "read_fastq: failed while reading %s", path ? path : "");
        return -1;
    }
    if (!fastq_line_is_header(&record->line)) {
        snprintf(errbuf, errbuf_len, "read_fastq: invalid FASTQ header in %s", path ? path : "");
        return -1;
    }
    ret = fastq_extract_qname(&record->line, &record->name, keep_name);
    if (ret == -2) {
        snprintf(errbuf, errbuf_len,
                 "read_fastq: query name exceeds %d bytes in %s",
                 BAM_MAX_QNAME_LEN, path ? path : "");
        return -1;
    }
    if (ret != 0) {
        snprintf(errbuf, errbuf_len, "read_fastq: out of memory parsing FASTQ header in %s", path ? path : "");
        return -1;
    }

    for (;;) {
        ret = hts_getline(fp, KS_SEP_LINE, &record->line);
        if (ret < 0) break;
        if (fastq_line_is_plus(&record->line)) break;
        if (record->line.l > (size_t)INT32_MAX - seq_len) {
            snprintf(errbuf, errbuf_len,
                     "read_fastq: sequence exceeds %d bytes in %s",
                     INT32_MAX, path ? path : "");
            return -1;
        }
        seq_len += record->line.l;
        if (fastq_append_line(&record->sequence, &record->line,
                              keep_sequence, path, "sequence",
                              errbuf, errbuf_len) != 0) {
            return -1;
        }
    }
    if (ret < 0) {
        snprintf(errbuf, errbuf_len,
                 "read_fastq: truncated FASTQ sequence block in %s", path ? path : "");
        return -1;
    }
    record->sequence_len = seq_len;

    /* htslib consumes one quality line even for a zero-length sequence. */
    remaining = seq_len;
    do {
        ret = hts_getline(fp, KS_SEP_LINE, &record->line);
        if (ret < 0) {
            snprintf(errbuf, errbuf_len,
                     "read_fastq: truncated FASTQ quality block in %s",
                     path ? path : "");
            return -1;
        }
        if (record->line.l > remaining) {
            snprintf(errbuf, errbuf_len,
                     "read_fastq: FASTQ quality is longer than sequence in %s",
                     path ? path : "");
            return -1;
        }
        if (fastq_append_line(&record->quality, &record->line,
                              keep_quality, path, "quality",
                              errbuf, errbuf_len) != 0) {
            return -1;
        }
        remaining -= record->line.l;
    } while (remaining > 0);

    return 1;
}

/* ================================================================
 * Bind (shared by fasta_read / fastq_read)
 * ================================================================ */

static void seq_read_bind(duckdb_bind_info info, int is_fastq) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);

    if (!file_path || strlen(file_path) == 0) {
        duckdb_bind_set_error(info,
            is_fastq ? "read_fastq requires a file path"
                     : "read_fasta requires a file path");
        if (file_path) duckdb_free(file_path);
        return;
    }

    /* Verify the file opens and is the expected format */
    samFile *fp = sam_open(file_path, "r");
    if (!fp) {
        char err[512];
        snprintf(err, sizeof(err), "Failed to open file: %s", file_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        return;
    }

    /* Check format — sam_open auto-detects */
    enum htsExactFormat fmt = hts_get_format(fp)->format;
    if (is_fastq && fmt != fastq_format) {
        /* Allow it but warn — htslib will still read it */
    }
    if (!is_fastq && fmt != fasta_format) {
        /* Same — allow but it may not have quality */
    }
    sam_close(fp);

    seq_bind_data_t *bind = (seq_bind_data_t *)duckdb_malloc(sizeof(seq_bind_data_t));
    memset(bind, 0, sizeof(seq_bind_data_t));
    bind->file_path = file_path;
    bind->is_fastq = is_fastq;
    bind->direct_fastq = is_fastq && fmt == fastq_format;
    bind->qual_repr = DUCKHTS_QUALITY_REPR_STRING;
    bind->input_quality_encoding = DUCKHTS_QUALITY_ENCODING_PHRED33;

    duckdb_value scan_mode_val = duckdb_bind_get_named_parameter(info, "scan_mode");
    if (scan_mode_val && !duckdb_is_null_value(scan_mode_val)) {
        char *scan_mode = duckdb_get_varchar(scan_mode_val);
        if (!seq_parse_scan_mode(scan_mode, &bind->scan_sequential)) {
            duckdb_bind_set_error(info,
                is_fastq
                    ? "read_fastq: scan_mode must be 'auto' or 'sequential'"
                    : "read_fasta: scan_mode must be 'auto' or 'sequential'");
            if (scan_mode) duckdb_free(scan_mode);
            duckdb_destroy_value(&scan_mode_val);
            destroy_seq_bind(bind);
            return;
        }
        if (scan_mode) duckdb_free(scan_mode);
    }
    if (scan_mode_val) duckdb_destroy_value(&scan_mode_val);

    if (is_fastq) {
        duckdb_value mate_val = duckdb_bind_get_named_parameter(info, "mate_path");
        if (mate_val && !duckdb_is_null_value(mate_val)) {
            bind->mate_path = duckdb_get_varchar(mate_val);
            bind->paired = 1;
        }
        if (mate_val) duckdb_destroy_value(&mate_val);

        duckdb_value inter_val = duckdb_bind_get_named_parameter(info, "interleaved");
        if (inter_val && !duckdb_is_null_value(inter_val)) {
            bind->interleaved = duckdb_get_bool(inter_val) ? 1 : 0;
        }
        if (inter_val) duckdb_destroy_value(&inter_val);

        if (bind->paired && bind->interleaved) {
            duckdb_bind_set_error(info, "read_fastq: use mate_path or interleaved, not both");
            destroy_seq_bind(bind);
            return;
        }
    } else {
        duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
        if (region_val && !duckdb_is_null_value(region_val)) {
            bind->region = duckdb_get_varchar(region_val);
            parse_regions_duckdb(bind->region, &bind->regions, &bind->n_regions);
        }
        if (region_val) duckdb_destroy_value(&region_val);

        if (bind->scan_sequential && bind->region && bind->region[0] != '\0') {
            duckdb_bind_set_error(info, "read_fasta: scan_mode := 'sequential' is incompatible with region queries");
            destroy_seq_bind(bind);
            return;
        }

        duckdb_value index_val = duckdb_bind_get_named_parameter(info, "index_path");
        if (index_val && !duckdb_is_null_value(index_val)) {
            bind->index_path = duckdb_get_varchar(index_val);
        }
        if (index_val) duckdb_destroy_value(&index_val);

        duckdb_value gzi_val = duckdb_bind_get_named_parameter(info, "gzi_path");
        if (gzi_val && !duckdb_is_null_value(gzi_val)) {
            bind->gzi_path = duckdb_get_varchar(gzi_val);
        }
        if (gzi_val) duckdb_destroy_value(&gzi_val);
    }

    /* sequence_encoding — shared by both FASTA and FASTQ */
    duckdb_value enc_val = duckdb_bind_get_named_parameter(info, "sequence_encoding");
    if (enc_val && !duckdb_is_null_value(enc_val)) {
        char *enc = duckdb_get_varchar(enc_val);
        if (enc) {
            if (strcmp(enc, "nt16") == 0) {
                bind->seq_nt16 = 1;
            } else if (strcmp(enc, "string") != 0) {
                duckdb_bind_set_error(info,
                    is_fastq
                        ? "read_fastq: sequence_encoding must be 'string' or 'nt16'"
                        : "read_fasta: sequence_encoding must be 'string' or 'nt16'");
            }
            duckdb_free(enc);
        }
    }
    if (enc_val) duckdb_destroy_value(&enc_val);

    if (is_fastq) {
        duckdb_value qrepr_val = duckdb_bind_get_named_parameter(info, "quality_representation");
        if (qrepr_val && !duckdb_is_null_value(qrepr_val)) {
            char *repr = duckdb_get_varchar(qrepr_val);
            if (repr) {
                if (duckhts_parse_quality_representation(repr, &bind->qual_repr) != 0) {
                    duckdb_bind_set_error(info,
                        "read_fastq: quality_representation must be 'string' or 'phred'");
                    duckdb_free(repr);
                    if (qrepr_val) duckdb_destroy_value(&qrepr_val);
                    destroy_seq_bind(bind);
                    return;
                }
                duckdb_free(repr);
            }
        }
        if (qrepr_val) duckdb_destroy_value(&qrepr_val);

        duckdb_value qenc_val = duckdb_bind_get_named_parameter(info, "input_quality_encoding");
        if (qenc_val && !duckdb_is_null_value(qenc_val)) {
            char *enc = duckdb_get_varchar(qenc_val);
            if (enc) {
                if (duckhts_parse_quality_encoding(enc, 1, &bind->input_quality_encoding) != 0) {
                    duckdb_bind_set_error(info,
                        "read_fastq: input_quality_encoding must be 'auto', 'phred33', 'phred64', or 'solexa64'");
                    duckdb_free(enc);
                    if (qenc_val) duckdb_destroy_value(&qenc_val);
                    destroy_seq_bind(bind);
                    return;
                }
                duckdb_free(enc);
            }
        }
        if (qenc_val) duckdb_destroy_value(&qenc_val);

        if (bind->input_quality_encoding == DUCKHTS_QUALITY_ENCODING_AUTO) {
            duckhts_quality_detect_result detect;
            char err[512];
            if (duckhts_detect_fastq_quality_encoding(bind->file_path, 10000, &detect, err, sizeof(err)) != 0) {
                duckdb_bind_set_error(info, err);
                destroy_seq_bind(bind);
                return;
            }
            bind->input_quality_encoding = detect.guessed_encoding;
            if (bind->paired) {
                duckhts_quality_detect_result mate_detect;
                if (duckhts_detect_fastq_quality_encoding(bind->mate_path, 10000, &mate_detect, err, sizeof(err)) != 0) {
                    duckdb_bind_set_error(info, err);
                    destroy_seq_bind(bind);
                    return;
                }
                if (mate_detect.guessed_encoding != bind->input_quality_encoding) {
                    duckdb_bind_set_error(info,
                        "read_fastq: auto-detected mate FASTQ quality encoding does not match path");
                    destroy_seq_bind(bind);
                    return;
                }
            }
        }
    }

    /* Define schema */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);

    duckdb_bind_add_result_column(info, "NAME", varchar_type);
    duckdb_bind_add_result_column(info, "DESCRIPTION", varchar_type);
    if (bind->seq_nt16) {
        duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
        duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);
        duckdb_bind_add_result_column(info, "SEQUENCE", list_type);
        duckdb_destroy_logical_type(&utinyint_type);
        duckdb_destroy_logical_type(&list_type);
    } else {
        duckdb_bind_add_result_column(info, "SEQUENCE", varchar_type);
    }
    if (is_fastq) {
        if (bind->qual_repr == DUCKHTS_QUALITY_REPR_PHRED) {
            duckdb_logical_type utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
            duckdb_logical_type list_type = duckdb_create_list_type(utinyint_type);
            duckdb_bind_add_result_column(info, "QUALITY", list_type);
            duckdb_destroy_logical_type(&utinyint_type);
            duckdb_destroy_logical_type(&list_type);
        } else {
            duckdb_bind_add_result_column(info, "QUALITY", varchar_type);
        }
        if (bind->paired || bind->interleaved) {
            duckdb_bind_add_result_column(info, "MATE", usmallint_type);
            duckdb_bind_add_result_column(info, "PAIR_ID", varchar_type);
        }
    }

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&usmallint_type);
    duckdb_bind_set_bind_data(info, bind, destroy_seq_bind);
}

static void fasta_read_bind(duckdb_bind_info info) { seq_read_bind(info, 0); }
static void fastq_read_bind(duckdb_bind_info info) { seq_read_bind(info, 1); }

/* ================================================================
 * Init
 * ================================================================ */

static void seq_read_init(duckdb_init_info info) {
    seq_bind_data_t *bind = (seq_bind_data_t *)duckdb_init_get_bind_data(info);

    seq_init_data_t *init = (seq_init_data_t *)duckdb_malloc(sizeof(seq_init_data_t));
    memset(init, 0, sizeof(seq_init_data_t));

    /* Projection pushdown: record which columns DuckDB actually needs */
    init->column_count = duckdb_init_get_column_count(info);
    if (init->column_count > 0) {
        init->column_ids = (idx_t *)duckdb_malloc(sizeof(idx_t) * init->column_count);
        for (idx_t i = 0; i < init->column_count; i++)
            init->column_ids[i] = duckdb_init_get_column_index(info, i);
    }

    init->is_fastq = bind->is_fastq;
    init->direct_fastq = bind->direct_fastq;
    init->paired = bind->paired;
    init->interleaved = bind->interleaved;
    init->seq_nt16 = bind->seq_nt16;
    init->qual_repr = bind->qual_repr;
    init->input_quality_encoding = bind->input_quality_encoding;
    init->pending_mate = 0;
    init->interleaved_mate = 1;
    init->n_regions = bind->n_regions;
    init->next_region_idx = 0;
    init->regions = bind->regions;

    init->fastq_need_name = bind->paired;
    for (idx_t i = 0; i < init->column_count; i++) {
        switch (init->column_ids[i]) {
        case SEQ_COL_NAME:
        case SEQ_COL_PAIR_ID:
            init->fastq_need_name = 1;
            break;
        case SEQ_COL_SEQUENCE:
            init->fastq_need_sequence = 1;
            break;
        case SEQ_COL_QUALITY:
            init->fastq_need_quality = 1;
            break;
        default:
            break;
        }
    }

    if (init->column_count == 0 && bind->n_regions == 0) {
        if (bind->is_fastq && !bind->paired && !bind->interleaved && !bind->scan_sequential) {
            faidx_t *fqi = fai_load3_format(bind->file_path, NULL, NULL, 0, FAI_FASTQ);
            if (fqi) {
                init->count_only = 1;
                init->count_remaining = (uint64_t)faidx_nseq(fqi);
                init->done = (init->count_remaining == 0);
                fai_destroy(fqi);
                duckdb_init_set_max_threads(info, 1);
                duckdb_init_set_init_data(info, init, destroy_seq_init);
                return;
            }
        } else if (!bind->is_fastq) {
            faidx_t *fai_count = bind->scan_sequential
                ? NULL
                : fai_load3_format(bind->file_path, bind->index_path,
                                   bind->gzi_path, 0, FAI_FASTA);
            if (fai_count) {
                init->count_only = 1;
                init->count_remaining = (uint64_t)faidx_nseq(fai_count);
                init->done = (init->count_remaining == 0);
                fai_destroy(fai_count);
                duckdb_init_set_max_threads(info, 1);
                duckdb_init_set_init_data(info, init, destroy_seq_init);
                return;
            }
            init->fasta_count_only = 1;
            init->count_fp = hts_open(bind->file_path, "r");
            if (!init->count_fp) {
                duckdb_init_set_error(info, "Failed to open FASTA file");
                destroy_seq_init(init);
                return;
            }
            duckhts_apply_remote_hts_tuning(init->count_fp, bind->file_path,
                                            DUCKHTS_HTS_IO_PROFILE_STREAMING);
            init->done = 0;
            duckdb_init_set_max_threads(info, 1);
            duckdb_init_set_init_data(info, init, destroy_seq_init);
            return;
        }
    }

    if (bind->direct_fastq) {
        init->count_fp = hts_open(bind->file_path, "r");
        if (!init->count_fp) {
            duckdb_init_set_error(info, "Failed to open FASTQ file");
            destroy_seq_init(init);
            return;
        }
        duckhts_apply_remote_hts_tuning(init->count_fp, bind->file_path,
                                        DUCKHTS_HTS_IO_PROFILE_STREAMING);
        if (bind->paired) {
            init->count_fp_mate = hts_open(bind->mate_path, "r");
            if (!init->count_fp_mate) {
                duckdb_init_set_error(info, "Failed to open mate FASTQ file");
                destroy_seq_init(init);
                return;
            }
            duckhts_apply_remote_hts_tuning(init->count_fp_mate, bind->mate_path,
                                            DUCKHTS_HTS_IO_PROFILE_STREAMING);
        }
        init->done = 0;
        duckdb_init_set_max_threads(info, 1);
        duckdb_init_set_init_data(info, init, destroy_seq_init);
        return;
    }

    /* sam_open handles .gz / .bgzf transparently */
    init->fp = sam_open(bind->file_path, "r");
    if (!init->fp) {
        duckdb_init_set_error(info, "Failed to open sequence file");
        destroy_seq_init(init);
        return;
    }
    duckhts_apply_remote_hts_tuning(init->fp, bind->file_path,
                                    DUCKHTS_HTS_IO_PROFILE_STREAMING);

    /* sam_hdr_read for FASTA/FASTQ returns a valid (possibly empty) header */
    init->hdr = sam_hdr_read(init->fp);
    if (!init->hdr) {
        sam_close(init->fp); init->fp = NULL;
        duckdb_init_set_error(info, "Failed to read header");
        destroy_seq_init(init);
        return;
    }

    init->rec = bam_init1();

    if (bind->paired) {
        init->fp_mate = sam_open(bind->mate_path, "r");
        if (!init->fp_mate) {
            duckdb_init_set_error(info, "Failed to open mate FASTQ file");
            destroy_seq_init(init);
            return;
        }
        duckhts_apply_remote_hts_tuning(init->fp_mate, bind->mate_path,
                                        DUCKHTS_HTS_IO_PROFILE_STREAMING);
        init->hdr_mate = sam_hdr_read(init->fp_mate);
        if (!init->hdr_mate) {
            duckdb_init_set_error(info, "Failed to read mate FASTQ header");
            destroy_seq_init(init);
            return;
        }
        init->rec_mate = bam_init1();
    }

    if (!bind->is_fastq && bind->n_regions > 0) {
        init->fai = fai_load3_format(bind->file_path, bind->index_path,
                                     bind->gzi_path, 0, FAI_FASTA);
        if (!init->fai) {
            duckdb_init_set_error(info, "read_fasta: region query requires a FASTA index (.fai); run fasta_index(path) first");
            destroy_seq_init(init);
            return;
        }
        duckhts_apply_remote_faidx_tuning(init->fai, bind->file_path,
                                          DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    }

    init->done = 0;

    duckdb_init_set_max_threads(info, 1);
    duckdb_init_set_init_data(info, init, destroy_seq_init);
}

static uint8_t fastq_stored_quality(char ascii_quality) {
    return (uint8_t)((unsigned char)ascii_quality - (unsigned char)'!');
}

static int seq_fastq_materialize_record(duckdb_function_info info,
                                        duckdb_data_chunk output,
                                        seq_init_data_t *init,
                                        const fastq_record_t *record,
                                        idx_t row, int mate) {
    size_t seq_len = record->sequence_len;

    for (idx_t i = 0; i < init->column_count; i++) {
        idx_t col_id = init->column_ids[i];
        duckdb_vector vec = duckdb_data_chunk_get_vector(output, i);

        switch (col_id) {
        case SEQ_COL_NAME:
            duckdb_vector_assign_string_element_len(
                vec, row, record->name.s ? record->name.s : "", record->name.l);
            break;

        case SEQ_COL_DESCRIPTION:
            set_null(vec, row);
            break;

        case SEQ_COL_SEQUENCE:
            if (init->seq_nt16) {
                duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                duckdb_list_entry *list_data =
                    (duckdb_list_entry *)duckdb_vector_get_data(vec);
                idx_t child_offset = duckdb_list_vector_get_size(vec);
                if (seq_len > 0) {
                    duckdb_list_vector_reserve(vec, child_offset + (idx_t)seq_len);
                    duckdb_list_vector_set_size(vec, child_offset + (idx_t)seq_len);
                    uint8_t *child_data =
                        (uint8_t *)duckdb_vector_get_data(child_vec);
                    seq_text_to_nt16(record->sequence.s,
                                     child_data + child_offset,
                                     (idx_t)seq_len);
                }
                list_data[row].offset = child_offset;
                list_data[row].length = (idx_t)seq_len;
            } else if (seq_len > 0) {
                if (ensure_buf(&init->seq_buf, &init->seq_buf_cap, seq_len) != 0) {
                    duckdb_function_set_error(
                        info, "read_fastq: out of memory allocating sequence buffer");
                    return -1;
                }
                for (size_t k = 0; k < seq_len; k++) {
                    uint8_t code = seq_nt16_table[(unsigned char)record->sequence.s[k]];
                    init->seq_buf[k] = seq_nt16_str[code];
                }
                init->seq_buf[seq_len] = '\0';
                duckdb_vector_assign_string_element_len(
                    vec, row, init->seq_buf, (idx_t)seq_len);
            } else {
                duckdb_vector_assign_string_element(vec, row, "");
            }
            break;

        case SEQ_COL_QUALITY:
            if (seq_len > 0 &&
                fastq_stored_quality(record->quality.s[0]) != UINT8_MAX) {
                if (init->qual_repr == DUCKHTS_QUALITY_REPR_PHRED) {
                    duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                    duckdb_list_entry *list_data =
                        (duckdb_list_entry *)duckdb_vector_get_data(vec);
                    idx_t child_offset = duckdb_list_vector_get_size(vec);
                    duckdb_list_vector_reserve(vec,
                                               child_offset + (idx_t)seq_len);
                    duckdb_list_vector_set_size(vec,
                                                child_offset + (idx_t)seq_len);
                    uint8_t *child_data =
                        (uint8_t *)duckdb_vector_get_data(child_vec);
                    for (size_t k = 0; k < seq_len; k++) {
                        child_data[child_offset + k] =
                            duckhts_normalize_fastq_quality(
                                fastq_stored_quality(record->quality.s[k]),
                                init->input_quality_encoding);
                    }
                    list_data[row].offset = child_offset;
                    list_data[row].length = (idx_t)seq_len;
                } else if (init->input_quality_encoding ==
                           DUCKHTS_QUALITY_ENCODING_PHRED33) {
                    duckdb_vector_assign_string_element_len(
                        vec, row, record->quality.s, (idx_t)seq_len);
                } else {
                    if (ensure_buf(&init->qual_buf, &init->qual_buf_cap,
                                   seq_len) != 0) {
                        duckdb_function_set_error(
                            info, "read_fastq: out of memory allocating quality buffer");
                        return -1;
                    }
                    for (size_t k = 0; k < seq_len; k++) {
                        init->qual_buf[k] = (char)(
                            duckhts_normalize_fastq_quality(
                                fastq_stored_quality(record->quality.s[k]),
                                init->input_quality_encoding) + 33);
                    }
                    init->qual_buf[seq_len] = '\0';
                    duckdb_vector_assign_string_element_len(
                        vec, row, init->qual_buf, (idx_t)seq_len);
                }
            } else {
                set_null(vec, row);
            }
            break;

        case SEQ_COL_MATE: {
            uint16_t *data = (uint16_t *)duckdb_vector_get_data(vec);
            data[row] = (uint16_t)mate;
            break;
        }

        case SEQ_COL_PAIR_ID:
            duckdb_vector_assign_string_element_len(
                vec, row, record->name.s ? record->name.s : "", record->name.l);
            break;

        default:
            break;
        }
    }

    return 0;
}

static void seq_fastq_direct_function(duckdb_function_info info,
                                      duckdb_data_chunk output,
                                      seq_init_data_t *init,
                                      const seq_bind_data_t *bind) {
    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;
    char err[512];

    while (row_count < vector_size) {
        fastq_record_t *record;
        int mate = 0;

        if (init->paired) {
            if (init->pending_mate) {
                record = &init->fastq_record_mate;
                mate = 2;
                init->pending_mate = 0;
            } else {
                int r1 = fastq_read_record(
                    init->count_fp, &init->fastq_record,
                    init->fastq_need_name, init->fastq_need_sequence,
                    init->fastq_need_quality, bind->file_path,
                    err, sizeof(err));
                if (r1 < 0) {
                    duckdb_function_set_error(info, err);
                    goto error;
                }
                int r2 = fastq_read_record(
                    init->count_fp_mate, &init->fastq_record_mate,
                    init->fastq_need_name, init->fastq_need_sequence,
                    init->fastq_need_quality, bind->mate_path,
                    err, sizeof(err));
                if (r2 < 0) {
                    duckdb_function_set_error(info, err);
                    goto error;
                }
                if (r1 == 0 || r2 == 0) {
                    if (r1 == 0 && r2 == 0) {
                        init->done = 1;
                        break;
                    }
                    duckdb_function_set_error(
                        info, "read_fastq: mate files have different record counts");
                    goto error;
                }
                if (init->fastq_record.name.l != init->fastq_record_mate.name.l ||
                    memcmp(init->fastq_record.name.s,
                           init->fastq_record_mate.name.s,
                           init->fastq_record.name.l) != 0) {
                    char msg[768];
                    snprintf(msg, sizeof(msg),
                             "read_fastq: mate files out of sync "
                             "(QNAME mismatch: '%s' vs '%s')",
                             init->fastq_record.name.s,
                             init->fastq_record_mate.name.s);
                    duckdb_function_set_error(info, msg);
                    goto error;
                }
                record = &init->fastq_record;
                mate = 1;
                init->pending_mate = 1;
            }
        } else {
            int ret = fastq_read_record(
                init->count_fp, &init->fastq_record,
                init->fastq_need_name, init->fastq_need_sequence,
                init->fastq_need_quality, bind->file_path,
                err, sizeof(err));
            if (ret < 0) {
                duckdb_function_set_error(info, err);
                goto error;
            }
            if (ret == 0) {
                if (init->interleaved && init->interleaved_mate == 2) {
                    duckdb_function_set_error(
                        info, "read_fastq: interleaved file has an unpaired record");
                    goto error;
                }
                init->done = 1;
                break;
            }
            record = &init->fastq_record;
            if (init->interleaved) {
                mate = init->interleaved_mate;
                init->interleaved_mate =
                    init->interleaved_mate == 1 ? 2 : 1;
            }
        }

        if (seq_fastq_materialize_record(info, output, init, record,
                                         row_count, mate) != 0) {
            goto error;
        }
        row_count++;
    }

    duckdb_data_chunk_set_size(output, row_count);
    return;

error:
    init->done = 1;
    duckdb_data_chunk_set_size(output, 0);
}

static void seq_fasta_count_only_function(duckdb_function_info info,
                                          duckdb_data_chunk output,
                                          seq_init_data_t *init,
                                          const seq_bind_data_t *bind) {
    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;
    int ret;

    while (row_count < vector_size) {
        ret = hts_getline(init->count_fp, KS_SEP_LINE, &init->count_line);
        if (ret == -1) {
            init->done = 1;
            break;
        }
        if (ret < -1) {
            char msg[512];
            snprintf(msg, sizeof(msg), "read_fasta: failed while parsing %s",
                     bind->file_path ? bind->file_path : "");
            duckdb_function_set_error(info, msg);
            init->done = 1;
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
        if (init->count_line.l == 0) continue;
        if (init->count_line.s[0] == '>') {
            row_count++;
        }
    }

    duckdb_data_chunk_set_size(output, row_count);
}

/* ================================================================
 * Scan
 *
 * FASTQ takes the direct projected path above.  FASTA uses sam_read1() to
 * read each record into a bam1_t and extracts fields with bam_get_* macros.
 *
 * Reference: htslib-1.23/samples/read_fast.c
 * ================================================================ */

static void seq_read_function(duckdb_function_info info, duckdb_data_chunk output) {
    seq_bind_data_t *bind = (seq_bind_data_t *)duckdb_function_get_bind_data(info);
    seq_init_data_t *init = (seq_init_data_t *)duckdb_function_get_init_data(info);

    if (!init || init->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    if (init->count_only) {
        idx_t row_count = (init->count_remaining > (uint64_t)duckdb_vector_size())
            ? duckdb_vector_size()
            : (idx_t)init->count_remaining;
        init->count_remaining -= row_count;
        if (init->count_remaining == 0) {
            init->done = 1;
        }
        duckdb_data_chunk_set_size(output, row_count);
        return;
    }

    if (init->direct_fastq) {
        seq_fastq_direct_function(info, output, init, bind);
        return;
    }

    if (init->fasta_count_only) {
        seq_fasta_count_only_function(info, output, init, bind);
        return;
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;

    while (row_count < vector_size) {
        if (init->fai && init->n_regions > 0) {
            if (init->next_region_idx >= init->n_regions) {
                init->done = 1;
                break;
            }
            const char *region = init->regions[init->next_region_idx++];
            hts_pos_t len = 0;
            char *seq = fai_fetch64(init->fai, region, &len);
            if (!seq || len < 0) {
                if (seq) free(seq);
                char msg[512];
                snprintf(msg, sizeof(msg), "read_fasta: invalid or missing region '%s'", region ? region : "");
                duckdb_function_set_error(info, msg);
                init->done = 1;
                duckdb_data_chunk_set_size(output, 0);
                return;
            }

            const char *name = region;
            size_t name_len = strlen(region);
            const char *colon = strchr(region, ':');
            if (colon) name_len = (size_t)(colon - region);
            if (ensure_buf(&init->pair_buf, &init->pair_buf_cap,
                           name_len) != 0) {
                free(seq);
                duckdb_function_set_error(info, "read_fasta: out of memory allocating name buffer");
                init->done = 1;
                duckdb_data_chunk_set_size(output, 0);
                return;
            }
            memcpy(init->pair_buf, name, name_len);
            init->pair_buf[name_len] = '\0';

            for (idx_t i = 0; i < init->column_count; i++) {
                idx_t col_id = init->column_ids[i];
                duckdb_vector vec = duckdb_data_chunk_get_vector(output, i);
                if (col_id == SEQ_COL_NAME) {
                    duckdb_vector_assign_string_element(vec, row_count, init->pair_buf);
                } else if (col_id == SEQ_COL_DESCRIPTION) {
                    set_null(vec, row_count);
                } else if (col_id == SEQ_COL_SEQUENCE) {
                    if (init->seq_nt16) {
                        /* Encode text sequence to nt16 LIST(UTINYINT)
                           using shared helper (permissive: unknown → N). */
                        duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                        duckdb_list_entry *list_data =
                            (duckdb_list_entry *)duckdb_vector_get_data(vec);
                        idx_t child_offset = duckdb_list_vector_get_size(vec);
                        if (len > 0) {
                            duckdb_list_vector_reserve(vec, child_offset + len);
                            duckdb_list_vector_set_size(vec, child_offset + len);
                            uint8_t *child_data =
                                (uint8_t *)duckdb_vector_get_data(child_vec);
                            seq_text_to_nt16(seq, child_data + child_offset, (idx_t)len);
                            list_data[row_count].offset = child_offset;
                            list_data[row_count].length = (idx_t)len;
                        } else {
                            list_data[row_count].offset = child_offset;
                            list_data[row_count].length = 0;
                        }
                    } else {
                        duckdb_vector_assign_string_element_len(vec, row_count, seq, len);
                    }
                }
            }
            free(seq);
            row_count++;
            continue;
        }

        bam1_t *b = NULL;
        int mate = 0;
        if (init->paired) {
            if (init->pending_mate) {
                b = init->rec_mate;
                mate = 2;
                init->pending_mate = 0;
            } else {
                int r1 = sam_read1(init->fp, init->hdr, init->rec);
                int r2 = sam_read1(init->fp_mate, init->hdr_mate, init->rec_mate);
                if (r1 < 0 || r2 < 0) {
                    if (r1 < 0 && r2 < 0) {
                        init->done = 1;
                        break;
                    }
                    duckdb_function_set_error(info,
                        "read_fastq: mate files have different record counts");
                    init->done = 1;
                    duckdb_data_chunk_set_size(output, 0);
                    return;
                }

                const char *q1 = bam_get_qname(init->rec);
                const char *q2 = bam_get_qname(init->rec_mate);
                if (!q1 || !q2 || strcmp(q1, q2) != 0) {
                    char msg[256];
                    snprintf(msg, sizeof(msg),
                        "read_fastq: mate files out of sync (QNAME mismatch: '%s' vs '%s')",
                        q1 ? q1 : "", q2 ? q2 : "");
                    duckdb_function_set_error(info, msg);
                    init->done = 1;
                    duckdb_data_chunk_set_size(output, 0);
                    return;
                }
                b = init->rec;
                mate = 1;
                init->pending_mate = 1;
            }
        } else {
            int ret = sam_read1(init->fp, init->hdr, init->rec);
            if (ret < 0) {
                /* -1 = EOF, < -1 = error */
                if (init->interleaved && init->interleaved_mate == 2) {
                    duckdb_function_set_error(info,
                        "read_fastq: interleaved file has an unpaired record");
                    init->done = 1;
                    duckdb_data_chunk_set_size(output, 0);
                    return;
                }
                init->done = 1;
                break;
            }
            b = init->rec;
            if (init->interleaved) {
                mate = init->interleaved_mate;
                init->interleaved_mate = (init->interleaved_mate == 1) ? 2 : 1;
            }
        }
        int seq_len = b->core.l_qseq;

        for (idx_t i = 0; i < init->column_count; i++) {
            idx_t col_id = init->column_ids[i];
            duckdb_vector vec = duckdb_data_chunk_get_vector(output, i);

            switch (col_id) {

            case SEQ_COL_NAME: {
                /* bam_get_qname returns the sequence name */
                const char *name = bam_get_qname(b);
                duckdb_vector_assign_string_element(vec, row_count,
                                                     name ? name : "");
                break;
            }

            case SEQ_COL_DESCRIPTION: {
                /* For FASTA/FASTQ, htslib stores the comment/description
                 * as a "CO" aux tag on the bam record */
                uint8_t *aux = bam_aux_get(b, "CO");
                if (aux) {
                    const char *comment = bam_aux2Z(aux);
                    if (comment) {
                        duckdb_vector_assign_string_element(vec, row_count,
                                                             comment);
                    } else {
                        set_null(vec, row_count);
                    }
                } else {
                    set_null(vec, row_count);
                }
                break;
            }

            case SEQ_COL_SEQUENCE: {
                if (init->seq_nt16) {
                    /* Emit packed nt16 codes as LIST(UTINYINT) */
                    duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                    duckdb_list_entry *list_data =
                        (duckdb_list_entry *)duckdb_vector_get_data(vec);
                    idx_t child_offset = duckdb_list_vector_get_size(vec);
                    if (seq_len > 0) {
                        duckdb_list_vector_reserve(vec, child_offset + seq_len);
                        duckdb_list_vector_set_size(vec, child_offset + seq_len);
                        uint8_t *child_data =
                            (uint8_t *)duckdb_vector_get_data(child_vec);
                        const uint8_t *raw_seq = bam_get_seq(b);
                        for (int k = 0; k < seq_len; k++)
                            child_data[child_offset + k] =
                                (uint8_t)bam_seqi(raw_seq, k);
                        list_data[row_count].offset = child_offset;
                        list_data[row_count].length = (idx_t)seq_len;
                    } else {
                        list_data[row_count].offset = child_offset;
                        list_data[row_count].length = 0;
                    }
                } else {
                    if (seq_len > 0) {
                        if (ensure_buf(&init->seq_buf, &init->seq_buf_cap,
                                       (size_t)seq_len) != 0) {
                            duckdb_function_set_error(info, "read_seq: out of memory allocating sequence buffer");
                            init->done = 1;
                            duckdb_data_chunk_set_size(output, 0);
                            return;
                        }
                        seq_decode_to_string(bam_get_seq(b), seq_len, init->seq_buf);
                        duckdb_vector_assign_string_element(vec, row_count,
                                                             init->seq_buf);
                    } else {
                        duckdb_vector_assign_string_element(vec, row_count, "");
                    }
                }
                break;
            }

            case SEQ_COL_QUALITY: {
                if (seq_len > 0 && bam_get_qual(b)[0] != 255) {
                    if (init->qual_repr == DUCKHTS_QUALITY_REPR_PHRED) {
                        duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                        duckdb_list_entry *list_data =
                            (duckdb_list_entry *)duckdb_vector_get_data(vec);
                        idx_t child_offset = duckdb_list_vector_get_size(vec);
                        const uint8_t *raw_qual = bam_get_qual(b);
                        duckdb_list_vector_reserve(vec, child_offset + seq_len);
                        duckdb_list_vector_set_size(vec, child_offset + seq_len);
                        uint8_t *child_data = (uint8_t *)duckdb_vector_get_data(child_vec);
                        for (int k = 0; k < seq_len; k++) {
                            child_data[child_offset + k] =
                                duckhts_normalize_fastq_quality(raw_qual[k], init->input_quality_encoding);
                        }
                        list_data[row_count].offset = child_offset;
                        list_data[row_count].length = (idx_t)seq_len;
                    } else {
                        if (ensure_buf(&init->qual_buf, &init->qual_buf_cap,
                                       (size_t)seq_len) != 0) {
                            duckdb_function_set_error(info, "read_seq: out of memory allocating quality buffer");
                            init->done = 1;
                            duckdb_data_chunk_set_size(output, 0);
                            return;
                        }
                        decode_qual(bam_get_qual(b), seq_len, init->qual_buf,
                                    init->input_quality_encoding);
                        duckdb_vector_assign_string_element(vec, row_count,
                                                             init->qual_buf);
                    }
                } else {
                    set_null(vec, row_count);
                }
                break;
            }

            case SEQ_COL_MATE: {
                if (init->is_fastq && (init->paired || init->interleaved)) {
                    uint16_t *data = (uint16_t *)duckdb_vector_get_data(vec);
                    data[row_count] = (uint16_t)mate;
                } else {
                    set_null(vec, row_count);
                }
                break;
            }

            case SEQ_COL_PAIR_ID: {
                if (init->is_fastq && (init->paired || init->interleaved)) {
                    const char *name = bam_get_qname(b);
                    if (ensure_buf(&init->pair_buf, &init->pair_buf_cap,
                                   strlen(name)) != 0) {
                        duckdb_function_set_error(info, "read_seq: out of memory allocating pair buffer");
                        init->done = 1;
                        duckdb_data_chunk_set_size(output, 0);
                        return;
                    }
                    const char *pair_id = strip_pair_suffix(name, init->pair_buf, init->pair_buf_cap);
                    duckdb_vector_assign_string_element(vec, row_count, pair_id);
                } else {
                    set_null(vec, row_count);
                }
                break;
            }

            default:
                break;
            } /* switch */
        } /* for columns */

        row_count++;
    } /* while rows */

    duckdb_data_chunk_set_size(output, row_count);
}

/* ================================================================
 * Registration
 * ================================================================ */

void register_read_fasta_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_fasta");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "gzi_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "sequence_encoding", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "scan_mode", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);

    duckdb_table_function_set_bind(tf, fasta_read_bind);
    duckdb_table_function_set_init(tf, seq_read_init);
    duckdb_table_function_set_function(tf, seq_read_function);
    duckdb_table_function_supports_projection_pushdown(tf, true);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

typedef struct {
    char *index_path;
    int emitted;
} fasta_index_bind_t;

static void destroy_fasta_index_bind(void *data) {
    fasta_index_bind_t *b = (fasta_index_bind_t *)data;
    if (!b) return;
    if (b->index_path) duckdb_free(b->index_path);
    duckdb_free(b);
}

static void fasta_index_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);
    if (!file_path || strlen(file_path) == 0) {
        duckdb_bind_set_error(info, "fasta_index requires a file path");
        if (file_path) duckdb_free(file_path);
        return;
    }

    char *index_path = NULL;
    duckdb_value idx_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (idx_val && !duckdb_is_null_value(idx_val)) {
        index_path = duckdb_get_varchar(idx_val);
    }
    if (idx_val) duckdb_destroy_value(&idx_val);

    if (fai_build3(file_path, index_path, NULL) != 0) {
        char err[512];
        snprintf(err, sizeof(err), "fasta_index: failed to build index for %s", file_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        return;
    }

    if (!index_path) {
        index_path = default_fasta_index_path(file_path);
        if (!index_path) {
            duckdb_bind_set_error(info, "fasta_index: out of memory");
            duckdb_free(file_path);
            return;
        }
    }

    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "success", bool_type);
    duckdb_bind_add_result_column(info, "index_path", varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&varchar_type);

    fasta_index_bind_t *bind = (fasta_index_bind_t *)duckdb_malloc(sizeof(fasta_index_bind_t));
    bind->index_path = index_path;
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_fasta_index_bind);
    duckdb_free(file_path);
}

static void fasta_index_init(duckdb_init_info info) {
    fasta_index_bind_t *bind = (fasta_index_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

static void fasta_index_scan(duckdb_function_info info, duckdb_data_chunk output) {
    fasta_index_bind_t *bind = (fasta_index_bind_t *)duckdb_function_get_bind_data(info);
    if (bind->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }
    duckdb_vector success_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector index_vec = duckdb_data_chunk_get_vector(output, 1);
    bool *success_data = (bool *)duckdb_vector_get_data(success_vec);
    success_data[0] = true;
    duckdb_vector_assign_string_element(index_vec, 0, bind->index_path ? bind->index_path : "");
    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

void register_fasta_index_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "fasta_index");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);

    duckdb_table_function_set_bind(tf, fasta_index_bind);
    duckdb_table_function_set_init(tf, fasta_index_init);
    duckdb_table_function_set_function(tf, fasta_index_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_fastq_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_fastq");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "mate_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "sequence_encoding", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "quality_representation", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "input_quality_encoding", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "scan_mode", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);

    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_named_parameter(tf, "interleaved", bool_type);
    duckdb_destroy_logical_type(&bool_type);

    duckdb_table_function_set_bind(tf, fastq_read_bind);
    duckdb_table_function_set_init(tf, seq_read_init);
    duckdb_table_function_set_function(tf, seq_read_function);
    duckdb_table_function_supports_projection_pushdown(tf, true);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
