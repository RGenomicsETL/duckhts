/**
 * DuckHTS FASTA/FASTQ Reader
 *
 * Table functions for reading FASTA and FASTQ files via the htslib SAM API.
 * htslib's sam_open() auto-detects fasta/fastq format (including .gz) and
 * sam_read1() returns records in bam1_t, giving us access to:
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
 * Note: htslib's FASTA/FASTQ reader stores the comment/description in
 * bam_get_l_aux(b) area as a "CO" aux tag when present.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

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
    int is_fastq;
    int interleaved;
    int paired;
} seq_bind_data_t;

/* ================================================================
 * Init Data
 * ================================================================ */

typedef struct {
    samFile *fp;
    sam_hdr_t *hdr;
    bam1_t *rec;
    samFile *fp_mate;
    sam_hdr_t *hdr_mate;
    bam1_t *rec_mate;
    int done;
    int is_fastq;
    int interleaved;
    int paired;
    int pending_mate;
    int interleaved_mate;

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
} seq_init_data_t;

/* ================================================================
 * Destructors
 * ================================================================ */

static void destroy_seq_bind(void *data) {
    seq_bind_data_t *b = (seq_bind_data_t *)data;
    if (!b) return;
    if (b->file_path) duckdb_free(b->file_path);
    if (b->mate_path) duckdb_free(b->mate_path);
    duckdb_free(b);
}

static void destroy_seq_init(void *data) {
    seq_init_data_t *init = (seq_init_data_t *)data;
    if (!init) return;
    if (init->rec) bam_destroy1(init->rec);
    if (init->hdr) sam_hdr_destroy(init->hdr);
    if (init->fp) sam_close(init->fp);
    if (init->rec_mate) bam_destroy1(init->rec_mate);
    if (init->hdr_mate) sam_hdr_destroy(init->hdr_mate);
    if (init->fp_mate) sam_close(init->fp_mate);
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

static void ensure_buf(char **buf, size_t *cap, int need) {
    size_t n = (size_t)(need + 1);
    if (n > *cap) {
        *cap = n * 2;
        *buf = (char *)realloc(*buf, *cap);
    }
}

/* Decode packed sequence via seq_nt16_str[] and bam_seqi() from htslib */
static void decode_seq(const uint8_t *seq_data, int len, char *buf) {
    for (int i = 0; i < len; i++)
        buf[i] = seq_nt16_str[bam_seqi(seq_data, i)];
    buf[len] = '\0';
}

/* Decode quality to Phred+33 ASCII */
static void decode_qual(const uint8_t *qual, int len, char *buf) {
    for (int i = 0; i < len; i++)
        buf[i] = (char)(qual[i] + 33);
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
    }

    /* Define schema */
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);

    duckdb_bind_add_result_column(info, "NAME", varchar_type);
    duckdb_bind_add_result_column(info, "DESCRIPTION", varchar_type);
    duckdb_bind_add_result_column(info, "SEQUENCE", varchar_type);
    if (is_fastq) {
        duckdb_bind_add_result_column(info, "QUALITY", varchar_type);
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

    /* sam_open handles .gz / .bgzf transparently */
    init->fp = sam_open(bind->file_path, "r");
    if (!init->fp) {
        duckdb_init_set_error(info, "Failed to open sequence file");
        duckdb_free(init);
        return;
    }

    /* sam_hdr_read for FASTA/FASTQ returns a valid (possibly empty) header */
    init->hdr = sam_hdr_read(init->fp);
    if (!init->hdr) {
        sam_close(init->fp); init->fp = NULL;
        duckdb_init_set_error(info, "Failed to read header");
        duckdb_free(init);
        return;
    }

    init->rec = bam_init1();
    init->is_fastq = bind->is_fastq;
    init->paired = bind->paired;
    init->interleaved = bind->interleaved;
    init->pending_mate = 0;
    init->interleaved_mate = 1;

    if (bind->paired) {
        init->fp_mate = sam_open(bind->mate_path, "r");
        if (!init->fp_mate) {
            duckdb_init_set_error(info, "Failed to open mate FASTQ file");
            destroy_seq_init(init);
            return;
        }
        init->hdr_mate = sam_hdr_read(init->fp_mate);
        if (!init->hdr_mate) {
            duckdb_init_set_error(info, "Failed to read mate FASTQ header");
            destroy_seq_init(init);
            return;
        }
        init->rec_mate = bam_init1();
    }

    init->done = 0;

    /* Projection pushdown */
    init->column_count = duckdb_init_get_column_count(info);
    init->column_ids = (idx_t *)duckdb_malloc(sizeof(idx_t) * init->column_count);
    for (idx_t i = 0; i < init->column_count; i++)
        init->column_ids[i] = duckdb_init_get_column_index(info, i);

    duckdb_init_set_max_threads(info, 1);
    duckdb_init_set_init_data(info, init, destroy_seq_init);
}

/* ================================================================
 * Scan
 *
 * Uses sam_read1() to read each FASTA/FASTQ record into a bam1_t,
 * then extracts fields using the standard bam_get_* macros.
 *
 * Reference: htslib-1.23/samples/read_fast.c
 * ================================================================ */

static void seq_read_function(duckdb_function_info info, duckdb_data_chunk output) {
    seq_init_data_t *init = (seq_init_data_t *)duckdb_function_get_init_data(info);

    if (!init || init->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;

    while (row_count < vector_size) {
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
                    init->done = 1;
                    break;
                }
                b = init->rec;
                mate = 1;
                init->pending_mate = 1;
            }
        } else {
            int ret = sam_read1(init->fp, init->hdr, init->rec);
        if (ret < 0) {
            /* -1 = EOF, < -1 = error */
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
                if (seq_len > 0) {
                    ensure_buf(&init->seq_buf, &init->seq_buf_cap, seq_len);
                    decode_seq(bam_get_seq(b), seq_len, init->seq_buf);
                    duckdb_vector_assign_string_element(vec, row_count,
                                                         init->seq_buf);
                } else {
                    duckdb_vector_assign_string_element(vec, row_count, "");
                }
                break;
            }

            case SEQ_COL_QUALITY: {
                if (seq_len > 0 && bam_get_qual(b)[0] != 255) {
                    ensure_buf(&init->qual_buf, &init->qual_buf_cap, seq_len);
                    decode_qual(bam_get_qual(b), seq_len, init->qual_buf);
                    duckdb_vector_assign_string_element(vec, row_count,
                                                         init->qual_buf);
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
                    ensure_buf(&init->pair_buf, &init->pair_buf_cap, (int)strlen(name));
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
    duckdb_destroy_logical_type(&varchar_type);

    duckdb_table_function_set_bind(tf, fasta_read_bind);
    duckdb_table_function_set_init(tf, seq_read_init);
    duckdb_table_function_set_function(tf, seq_read_function);
    duckdb_table_function_supports_projection_pushdown(tf, true);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_fastq_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_fastq");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "mate_path", varchar_type);
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
