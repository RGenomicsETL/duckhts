/* Canonical HTSlib transport and decode checks, independent of DuckDB/R. */
#include "include/bcf_scan.h"
#include "include/region_list.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>

void duckhts_bcf_scan_close(duckhts_bcf_scan_t *scan) {
    if (!scan) return;
    if (scan->itr) hts_itr_destroy(scan->itr);
    if (scan->hdr) bcf_hdr_destroy(scan->hdr);
    if (scan->fp) hts_close(scan->fp);
    ks_free(&scan->line);
    memset(scan, 0, sizeof(*scan));
}

int duckhts_bcf_scan_open(duckhts_bcf_scan_t *scan, const char *path,
                         const duckhts_bcf_index_t *index, int decompression_threads,
                         duckhts_hts_io_profile_t profile, const char *reader_name,
                         char *error, size_t error_size) {
    const char *name = reader_name ? reader_name : "read_bcf";
    scan->index = index;
    scan->fp = hts_open(path, "r");
    if (!scan->fp) {
        snprintf(error, error_size, "%s: failed to open BCF/VCF file: %s", name, path);
        goto fail;
    }
    duckhts_apply_remote_hts_tuning(scan->fp, path, profile);
    if (decompression_threads > 0 &&
        hts_get_format(scan->fp)->compression != no_compression &&
        hts_set_threads(scan->fp, decompression_threads) < 0) {
        snprintf(error, error_size, "%s: failed to configure BCF/VCF decompression threads", name);
        goto fail;
    }
    scan->hdr = bcf_hdr_read(scan->fp);
    if (!scan->hdr) {
        snprintf(error, error_size, "%s: failed to read BCF/VCF header", name);
        goto fail;
    }
    return 1;
fail:
    duckhts_bcf_scan_close(scan);
    return 0;
}

static int header_name2id(void *hdr, const char *name) {
    return bcf_hdr_name2id(hdr, name);
}

static int tabix_name2id(void *index, const char *name) {
    return tbx_name2id(index, name);
}

int duckhts_bcf_scan_regions(duckhts_bcf_scan_t *scan, char **regions,
                            unsigned int count, char *error, size_t error_size) {
    if (scan->itr) hts_itr_destroy(scan->itr);
    scan->itr = NULL;
    scan->indexed = 1;
    const duckhts_bcf_index_t *index = scan->index;
    if (!index || (!index->bcf && !index->tabix) || !regions || count == 0) {
        snprintf(error, error_size, "BCF scan requires an index and a nonempty region list");
        return 0;
    }
    if (!duckhts_region_list_validate(regions, count,
            index->bcf ? header_name2id : tabix_name2id,
            index->bcf ? (void *)scan->hdr : (void *)index->tabix, error, error_size)) {
        return 0;
    }
    errno = 0;
    if (index->bcf) {
        scan->itr = count == 1
            ? bcf_itr_querys(index->bcf, scan->hdr, regions[0])
            : bcf_itr_regarray(index->bcf, scan->hdr, regions, count);
    } else {
        scan->itr = count == 1
            ? tbx_itr_querys(index->tabix, regions[0])
            : tbx_itr_regarray(index->tabix, regions, count);
    }
    if (!scan->itr && errno == ENOMEM) {
        snprintf(error, error_size, "BCF scan: out of memory creating region iterator");
        return 0;
    }
    return 1;
}

int duckhts_bcf_scan_contig(duckhts_bcf_scan_t *scan, const char *name,
                           char *error, size_t error_size) {
    if (scan->itr) hts_itr_destroy(scan->itr);
    scan->itr = NULL;
    scan->indexed = 1;
    const duckhts_bcf_index_t *index = scan->index;
    if (!index || (!index->bcf && !index->tabix)) {
        snprintf(error, error_size, "BCF scan requires an index for contig queries");
        return 0;
    }
    int tid = index->bcf ? bcf_hdr_name2id(scan->hdr, name) : tbx_name2id(index->tabix, name);
    if (tid < 0) return 1;
    scan->itr = index->bcf
        ? bcf_itr_queryi(index->bcf, tid, 0, HTS_POS_MAX)
        : tbx_itr_queryi(index->tabix, tid, 0, HTS_POS_MAX);
    if (!scan->itr) {
        snprintf(error, error_size, "BCF scan: failed to create iterator for contig %s", name);
        return 0;
    }
    return 1;
}

int duckhts_bcf_scan_next(duckhts_bcf_scan_t *scan, bcf1_t *record) {
    if (!scan->indexed) return bcf_read(scan->fp, scan->hdr, record);
    if (!scan->itr) return -1;
    if (scan->index->tabix) {
        int ret = tbx_itr_next(scan->fp, scan->index->tabix, scan->itr, &scan->line);
        if (ret >= 0) {
            ret = vcf_parse1(&scan->line, scan->hdr, record);
            scan->line.l = 0;
        }
        return ret;
    }
    int ret = bcf_itr_next(scan->fp, scan->itr, record);
    /* bcf_read/vcf_parse already honor header selection. Indexed BCF reads
     * have no header argument, so apply its selection exactly once here. */
    if (ret >= 0 && scan->hdr->keep_samples) ret = bcf_subset_format(scan->hdr, record);
    return ret;
}

static const char *encoded_type_name(int type) {
    switch (type) {
    case BCF_BT_NULL: return "NULL";
    case BCF_BT_INT8: return "INT8";
    case BCF_BT_INT16: return "INT16";
    case BCF_BT_INT32: return "INT32";
    case BCF_BT_INT64: return "INT64";
    case BCF_BT_FLOAT: return "FLOAT";
    case BCF_BT_CHAR: return "CHAR";
    default: return "unknown";
    }
}

static const char *header_type_name(int type) {
    switch (type) {
    case BCF_HT_FLAG: return "Flag";
    case BCF_HT_INT: return "Integer";
    case BCF_HT_REAL: return "Float";
    case BCF_HT_STR: return "String";
    default: return "unknown";
    }
}

static void record_location(bcf_hdr_t *hdr, bcf1_t *record, char *buf, size_t size) {
    const char *chrom = hdr && record && record->rid >= 0 ? bcf_hdr_id2name(hdr, record->rid) : NULL;
    snprintf(buf, size, "%s:%lld", chrom ? chrom : "?", record ? (long long)record->pos + 1 : 0);
}

static int encoded_type_matches(int header_type, int encoded_type) {
    switch (header_type) {
    case BCF_HT_INT:
        return encoded_type == BCF_BT_INT8 || encoded_type == BCF_BT_INT16 || encoded_type == BCF_BT_INT32;
    case BCF_HT_REAL: return encoded_type == BCF_BT_FLOAT;
    case BCF_HT_STR: return encoded_type == BCF_BT_CHAR;
    default: return 1;
    }
}

int duckhts_bcf_check_field_type(bcf_hdr_t *hdr, bcf1_t *record,
                                duckhts_bcf_field_class_t field_class,
                                int header_id, int header_type,
                                const char *reader_name, char *error, size_t error_size) {
    if (!hdr || !record || header_id < 0) return 1;
    int encoded_type, expected_type = header_type;
    if (field_class == DUCKHTS_BCF_FIELD_INFO) {
        if (header_type == BCF_HT_FLAG) return 1;
        bcf_unpack(record, BCF_UN_INFO);
        bcf_info_t *info = bcf_get_info_id(record, header_id);
        if (!info || !info->vptr) return 1;
        encoded_type = info->type;
    } else {
        bcf_unpack(record, BCF_UN_FMT);
        bcf_fmt_t *fmt = bcf_get_fmt_id(record, header_id);
        if (!fmt || !fmt->p) return 1;
        encoded_type = fmt->type;
        expected_type = field_class == DUCKHTS_BCF_FIELD_GT ? BCF_HT_INT
            : header_type == BCF_HT_INT || header_type == BCF_HT_REAL ? header_type : BCF_HT_STR;
    }
    if (encoded_type_matches(expected_type, encoded_type)) return 1;
    char loc[128];
    record_location(hdr, record, loc, sizeof(loc));
    snprintf(error, error_size, "%s: %s/%s encoded BCF type %s does not match header Type=%s at %s",
             reader_name ? reader_name : "read_bcf",
             field_class == DUCKHTS_BCF_FIELD_INFO ? "INFO" : "FORMAT",
             bcf_hdr_int2id(hdr, BCF_DT_ID, header_id), encoded_type_name(encoded_type),
             field_class == DUCKHTS_BCF_FIELD_GT ? "String/GT-encoded-integer" : header_type_name(header_type), loc);
    return 0;
}

int duckhts_bcf_check_format_width(const char *reader_name, const char *tag,
                                  bcf_hdr_t *hdr, bcf1_t *record, int values,
                                  int samples, char *error, size_t error_size) {
    if (values <= 0 || samples <= 0 || values % samples == 0) return 1;
    char loc[128];
    record_location(hdr, record, loc, sizeof(loc));
    snprintf(error, error_size,
             "%s: FORMAT/%s decoded value count %d is not divisible by sample count %d at %s",
             reader_name ? reader_name : "read_bcf", tag ? tag : "?", values, samples, loc);
    return 0;
}

duckhts_bcf_decode_status_t duckhts_bcf_decode_status(
    const char *reader_name, const char *field_class, const char *tag,
    bcf_hdr_t *hdr, bcf1_t *record, int ret, char *error, size_t error_size) {
    if (ret >= 0 || ret == -3) return DUCKHTS_BCF_DECODE_OK;
    char loc[128];
    record_location(hdr, record, loc, sizeof(loc));
    const char *name = reader_name ? reader_name : "read_bcf";
    const char *klass = field_class ? field_class : "field";
    const char *field = tag ? tag : "?";
    if (ret == -1) {
        snprintf(error, error_size, "%s: %s/%s is not defined in the BCF/VCF header at %s", name, klass, field, loc);
    } else if (ret == -2) {
        snprintf(error, error_size, "%s: %s/%s encoded BCF type does not match the header at %s", name, klass, field, loc);
        return DUCKHTS_BCF_DECODE_TYPE_MISMATCH;
    } else if (ret == -4) {
        snprintf(error, error_size, "%s: out of memory decoding %s/%s at %s", name, klass, field, loc);
    } else {
        snprintf(error, error_size, "%s: failed to decode %s/%s at %s (htslib return %d)", name, klass, field, loc, ret);
    }
    return DUCKHTS_BCF_DECODE_FATAL;
}
