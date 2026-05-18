/*
 * DuckHTS BAM pileup reader.
 *
 * read_pileup(path, region := 'chr:start-end', index_path := NULL,
 *             min_mapq := 0, flag_mask := 1796)
 *
 * Exposes htslib's pileup engine as a region-scoped DuckDB table function.
 * This is not samtools mpileup text parity. We emit one row per covered
 * position with the observed bases and base qualities after simple BAM-level
 * filtering.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "include/hts_io_tuning.h"

typedef struct {
    char *file_path;
    char *index_path;
    char *region;
    int min_mapq;
    uint16_t flag_mask;
} bam_pileup_bind_data_t;

typedef struct {
    samFile *fp;
    sam_hdr_t *hdr;
    hts_idx_t *idx;
    hts_itr_t *itr;
    bam_plp_t plp;
    int done;
    int min_mapq;
    uint16_t flag_mask;
    int clip_to_region;
    int region_tid;
    hts_pos_t region_beg;
    hts_pos_t region_end;
    kstring_t bases_tmp;
    kstring_t quals_tmp;
} bam_pileup_local_init_data_t;

static int bam_pileup_fetch(void *data, bam1_t *b) {
    bam_pileup_local_init_data_t *local = (bam_pileup_local_init_data_t *)data;

    for (;;) {
        int ret = sam_itr_next(local->fp, local->itr, b);
        if (ret < 0) return ret;
        if (b->core.flag & local->flag_mask) continue;
        if (b->core.qual < local->min_mapq) continue;
        return ret;
    }
}

static void destroy_bam_pileup_bind(void *ptr) {
    bam_pileup_bind_data_t *bind = (bam_pileup_bind_data_t *)ptr;
    if (!bind) return;
    if (bind->file_path) duckdb_free(bind->file_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->region) duckdb_free(bind->region);
    duckdb_free(bind);
}

static void destroy_bam_pileup_local(void *ptr) {
    bam_pileup_local_init_data_t *local = (bam_pileup_local_init_data_t *)ptr;
    if (!local) return;
    if (local->plp) bam_plp_destroy(local->plp);
    if (local->itr) hts_itr_destroy(local->itr);
    if (local->idx) hts_idx_destroy(local->idx);
    if (local->hdr) sam_hdr_destroy(local->hdr);
    if (local->fp) sam_close(local->fp);
    ks_free(&local->bases_tmp);
    ks_free(&local->quals_tmp);
    duckdb_free(local);
}

static void bam_pileup_bind(duckdb_bind_info info) {
    bam_pileup_bind_data_t *bind = (bam_pileup_bind_data_t *)duckdb_malloc(sizeof(*bind));
    memset(bind, 0, sizeof(*bind));
    bind->flag_mask = 1796;

    duckdb_value file_val = duckdb_bind_get_parameter(info, 0);
    if (!file_val) {
        duckdb_bind_set_error(info, "read_pileup: BAM path is required");
        destroy_bam_pileup_bind(bind);
        return;
    }
    bind->file_path = duckdb_get_varchar(file_val);
    duckdb_destroy_value(&file_val);
    if (!bind->file_path || bind->file_path[0] == '\0') {
        duckdb_bind_set_error(info, "read_pileup: BAM path is required");
        destroy_bam_pileup_bind(bind);
        return;
    }

    duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
    if (!region_val || duckdb_is_null_value(region_val)) {
        duckdb_bind_set_error(info, "read_pileup: region := 'chr:start-end' is required");
        if (region_val) duckdb_destroy_value(&region_val);
        destroy_bam_pileup_bind(bind);
        return;
    }
    bind->region = duckdb_get_varchar(region_val);
    duckdb_destroy_value(&region_val);

    duckdb_value index_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (index_val && !duckdb_is_null_value(index_val)) {
        bind->index_path = duckdb_get_varchar(index_val);
    }
    if (index_val) duckdb_destroy_value(&index_val);

    duckdb_value mapq_val = duckdb_bind_get_named_parameter(info, "min_mapq");
    if (mapq_val && !duckdb_is_null_value(mapq_val)) {
        int64_t min_mapq = duckdb_get_int64(mapq_val);
        if (min_mapq < 0 || min_mapq > INT32_MAX) {
            duckdb_bind_set_error(info,
                "read_pileup: min_mapq must be between 0 and INT32_MAX");
        } else {
            bind->min_mapq = (int)min_mapq;
        }
    }
    if (mapq_val) duckdb_destroy_value(&mapq_val);

    duckdb_value flag_val = duckdb_bind_get_named_parameter(info, "flag_mask");
    if (flag_val && !duckdb_is_null_value(flag_val)) {
        int64_t flag_mask = duckdb_get_int64(flag_val);
        if (flag_mask < 0 || flag_mask > UINT16_MAX) {
            duckdb_bind_set_error(info,
                "read_pileup: flag_mask must be between 0 and 65535");
        } else {
            bind->flag_mask = (uint16_t)flag_mask;
        }
    }
    if (flag_val) duckdb_destroy_value(&flag_val);

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_bind_add_result_column(info, "chrom", varchar_type);
    duckdb_bind_add_result_column(info, "pos", bigint_type);
    duckdb_bind_add_result_column(info, "depth", int_type);
    duckdb_bind_add_result_column(info, "bases", varchar_type);
    duckdb_bind_add_result_column(info, "quals", varchar_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&int_type);

    duckdb_bind_set_bind_data(info, bind, destroy_bam_pileup_bind);
}

static void bam_pileup_global_init(duckdb_init_info info) {
    duckdb_init_set_max_threads(info, 1);
}

static void bam_pileup_local_init(duckdb_init_info info) {
    bam_pileup_bind_data_t *bind =
        (bam_pileup_bind_data_t *)duckdb_init_get_bind_data(info);
    bam_pileup_local_init_data_t *local =
        (bam_pileup_local_init_data_t *)duckdb_malloc(sizeof(*local));
    memset(local, 0, sizeof(*local));
    local->min_mapq = bind->min_mapq;
    local->flag_mask = bind->flag_mask;

    local->fp = sam_open(bind->file_path, "r");
    if (!local->fp) {
        duckdb_init_set_error(info, "read_pileup: failed to open BAM");
        destroy_bam_pileup_local(local);
        return;
    }
    duckhts_apply_remote_hts_tuning(local->fp, bind->file_path,
                                    DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    (void)hts_set_threads(local->fp, 2);

    local->hdr = sam_hdr_read(local->fp);
    if (!local->hdr) {
        duckdb_init_set_error(info, "read_pileup: failed to read BAM header");
        destroy_bam_pileup_local(local);
        return;
    }

    local->idx = sam_index_load3(local->fp, bind->file_path, bind->index_path,
                                 HTS_IDX_SILENT_FAIL);
    if (!local->idx) {
        duckdb_init_set_error(info,
            "read_pileup: BAM index (.bai/.csi) is required for region pileups");
        destroy_bam_pileup_local(local);
        return;
    }

    local->itr = sam_itr_querys(local->idx, local->hdr, bind->region);
    if (!local->itr) {
        char err[512];
        snprintf(err, sizeof(err),
                 "read_pileup: failed to parse region '%s'", bind->region);
        duckdb_init_set_error(info, err);
        destroy_bam_pileup_local(local);
        return;
    }

    if (strcmp(bind->region, ".") != 0 && strcmp(bind->region, "*") != 0) {
        const char *remaining = sam_parse_region(local->hdr, bind->region,
                                                 &local->region_tid,
                                                 &local->region_beg,
                                                 &local->region_end,
                                                 HTS_PARSE_THOUSANDS_SEP);
        if (!remaining || *remaining != '\0') {
            char err[512];
            snprintf(err, sizeof(err),
                     "read_pileup: failed to parse region '%s'", bind->region);
            duckdb_init_set_error(info, err);
            destroy_bam_pileup_local(local);
            return;
        }
        local->clip_to_region = 1;
    }

    local->plp = bam_plp_init(bam_pileup_fetch, local);
    if (!local->plp) {
        duckdb_init_set_error(info, "read_pileup: bam_plp_init failed");
        destroy_bam_pileup_local(local);
        return;
    }

    duckdb_init_set_init_data(info, local, destroy_bam_pileup_local);
}

static void bam_pileup_execute(duckdb_function_info info,
                               duckdb_data_chunk output) {
    bam_pileup_local_init_data_t *local =
        (bam_pileup_local_init_data_t *)duckdb_function_get_local_init_data(info);
    if (!local || local->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector depth_vec = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector bases_vec = duckdb_data_chunk_get_vector(output, 3);
    duckdb_vector quals_vec = duckdb_data_chunk_get_vector(output, 4);
    int64_t *pos_data = (int64_t *)duckdb_vector_get_data(pos_vec);
    int32_t *depth_data = (int32_t *)duckdb_vector_get_data(depth_vec);

    idx_t row_count = 0;
    idx_t vector_size = duckdb_vector_size();

    while (row_count < vector_size) {
        int tid = -1;
        int n_plp = 0;
        hts_pos_t pos = 0;
        const bam_pileup1_t *plp = bam_plp64_auto(local->plp, &tid, &pos, &n_plp);
        if (!plp) {
            local->done = 1;
            break;
        }

        if (local->clip_to_region) {
            if (tid != local->region_tid || pos < local->region_beg) {
                continue;
            }
            if (pos >= local->region_end) {
                local->done = 1;
                break;
            }
        }

        local->bases_tmp.l = 0;
        local->quals_tmp.l = 0;

        int actual_depth = 0;
        for (int i = 0; i < n_plp; i++) {
            const bam_pileup1_t *p = &plp[i];
            uint8_t *seq;
            uint8_t qual;
            char base;

            if (p->is_del || p->is_refskip) continue;

            seq = bam_get_seq(p->b);
            base = seq_nt16_str[bam_seqi(seq, p->qpos)];
            if (bam_is_rev(p->b)) {
                base = (char)tolower((unsigned char)base);
            }
            kputc(base, &local->bases_tmp);

            qual = bam_get_qual(p->b)[p->qpos];
            if (qual == 0xff) qual = 0;
            kputc((char)(qual + 33), &local->quals_tmp);
            actual_depth++;
        }

        if (actual_depth == 0) continue;

        duckdb_vector_assign_string_element(chrom_vec, row_count,
            (tid >= 0) ? sam_hdr_tid2name(local->hdr, tid) : "*");
        pos_data[row_count] = (int64_t)(pos + 1);
        depth_data[row_count] = actual_depth;
        duckdb_vector_assign_string_element_len(bases_vec, row_count,
                                                local->bases_tmp.s,
                                                local->bases_tmp.l);
        duckdb_vector_assign_string_element_len(quals_vec, row_count,
                                                local->quals_tmp.s,
                                                local->quals_tmp.l);
        row_count++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

void register_read_pileup_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_pileup");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "min_mapq", int_type);
    duckdb_table_function_add_named_parameter(tf, "flag_mask", int_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int_type);

    duckdb_table_function_set_bind(tf, bam_pileup_bind);
    duckdb_table_function_set_init(tf, bam_pileup_global_init);
    duckdb_table_function_set_local_init(tf, bam_pileup_local_init);
    duckdb_table_function_set_function(tf, bam_pileup_execute);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
