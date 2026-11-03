#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

#include "include/quality_encoding.h"

const char *duckhts_quality_encoding_name(duckhts_quality_encoding enc) {
    switch (enc) {
        case DUCKHTS_QUALITY_ENCODING_AUTO: return "auto";
        case DUCKHTS_QUALITY_ENCODING_PHRED33: return "phred33";
        case DUCKHTS_QUALITY_ENCODING_PHRED64: return "phred64";
        case DUCKHTS_QUALITY_ENCODING_SOLEXA64: return "solexa64";
        default: return "unknown";
    }
}

const char *duckhts_quality_representation_name(duckhts_quality_representation repr) {
    switch (repr) {
        case DUCKHTS_QUALITY_REPR_STRING: return "string";
        case DUCKHTS_QUALITY_REPR_PHRED: return "phred";
        default: return "unknown";
    }
}

int duckhts_parse_quality_encoding(const char *text, int allow_auto,
                                   duckhts_quality_encoding *out) {
    if (!text || !out) return -1;
    if (strcmp(text, "auto") == 0 && allow_auto) {
        *out = DUCKHTS_QUALITY_ENCODING_AUTO;
        return 0;
    }
    if (strcmp(text, "phred33") == 0) {
        *out = DUCKHTS_QUALITY_ENCODING_PHRED33;
        return 0;
    }
    if (strcmp(text, "phred64") == 0) {
        *out = DUCKHTS_QUALITY_ENCODING_PHRED64;
        return 0;
    }
    if (strcmp(text, "solexa64") == 0) {
        *out = DUCKHTS_QUALITY_ENCODING_SOLEXA64;
        return 0;
    }
    return -1;
}

int duckhts_parse_quality_representation(const char *text,
                                         duckhts_quality_representation *out) {
    if (!text || !out) return -1;
    if (strcmp(text, "string") == 0) {
        *out = DUCKHTS_QUALITY_REPR_STRING;
        return 0;
    }
    if (strcmp(text, "phred") == 0) {
        *out = DUCKHTS_QUALITY_REPR_PHRED;
        return 0;
    }
    return -1;
}

static uint8_t clamp_u8(int value) {
    if (value < 0) return 0;
    if (value > 255) return 255;
    return (uint8_t)value;
}

static uint8_t solexa_to_phred(int solexa) {
    double phred = 10.0 * log10(1.0 + pow(10.0, solexa / 10.0));
    return clamp_u8((int)lround(phred));
}

uint8_t duckhts_normalize_fastq_quality(uint8_t stored_qual,
                                        duckhts_quality_encoding enc) {
    switch (enc) {
        case DUCKHTS_QUALITY_ENCODING_PHRED64:
            return stored_qual > 31 ? (uint8_t)(stored_qual - 31) : 0;
        case DUCKHTS_QUALITY_ENCODING_SOLEXA64:
            return solexa_to_phred((int)stored_qual - 31);
        case DUCKHTS_QUALITY_ENCODING_PHRED33:
        case DUCKHTS_QUALITY_ENCODING_AUTO:
        default:
            return stored_qual;
    }
}

static int is_fastq_header(const char *line) {
    return line && line[0] == '@';
}

static int is_fastq_plus(const char *line) {
    return line && line[0] == '+';
}

static void update_ascii_range(duckhts_quality_detect_result *out,
                               const char *text, size_t len) {
    size_t i;
    for (i = 0; i < len; i++) {
        unsigned char c = (unsigned char)text[i];
        if ((int)c < out->observed_ascii_min) out->observed_ascii_min = (int)c;
        if ((int)c > out->observed_ascii_max) out->observed_ascii_max = (int)c;
    }
}

static void finalize_quality_detection(duckhts_quality_detect_result *out) {
    if (out->observed_ascii_min == 999 || out->records_sampled == 0) {
        out->observed_ascii_min = 0;
        out->observed_ascii_max = 0;
        out->guessed_encoding = DUCKHTS_QUALITY_ENCODING_PHRED33;
        return;
    }

    out->compatible_phred33 =
        out->observed_ascii_min >= 33 && out->observed_ascii_max <= 126;
    out->compatible_phred64 =
        out->observed_ascii_min >= 64 && out->observed_ascii_max <= 126;
    out->compatible_solexa64 =
        out->observed_ascii_min >= 59 && out->observed_ascii_max <= 126;

    if (out->observed_ascii_min < 59) {
        out->guessed_encoding = DUCKHTS_QUALITY_ENCODING_PHRED33;
        out->is_ambiguous = 0;
    } else if (out->observed_ascii_min < 64) {
        out->guessed_encoding = DUCKHTS_QUALITY_ENCODING_SOLEXA64;
        out->is_ambiguous = out->compatible_phred33 ? 1 : 0;
    } else {
        out->guessed_encoding = DUCKHTS_QUALITY_ENCODING_PHRED64;
        out->is_ambiguous = out->compatible_phred33 ? 1 : 0;
    }
}

int duckhts_detect_fastq_quality_encoding(const char *path, int64_t max_records,
                                          duckhts_quality_detect_result *out,
                                          char *errbuf, size_t errbuf_len) {
    htsFile *fp = NULL;
    kstring_t line = {0, 0, NULL};
    kstring_t seq = {0, 0, NULL};
    int ret = 0;
    int reached_eof = 0;

    if (!path || !out) return -1;
    memset(out, 0, sizeof(*out));
    out->observed_ascii_min = 999;
    if (errbuf && errbuf_len > 0) errbuf[0] = '\0';

    fp = hts_open(path, "r");
    if (!fp) {
        if (errbuf && errbuf_len > 0) {
            snprintf(errbuf, errbuf_len, "detect_quality_encoding: failed to open %s", path);
        }
        return -1;
    }

    while (max_records <= 0 || out->records_sampled < max_records) {
        size_t seq_len = 0;
        size_t qual_len = 0;

        ret = hts_getline(fp, KS_SEP_LINE, &line);
        if (ret == -1) {
            reached_eof = 1;
            break;
        }
        if (ret < -1 || !is_fastq_header(line.s)) {
            if (errbuf && errbuf_len > 0) {
                snprintf(errbuf, errbuf_len,
                         "detect_quality_encoding: invalid FASTQ header in %s", path);
            }
            ret = -1;
            break;
        }

        seq.l = 0;
        for (;;) {
            ret = hts_getline(fp, KS_SEP_LINE, &line);
            if (ret < 0) break;
            if (is_fastq_plus(line.s)) break;
            if (kputsn(line.s, line.l, &seq) < 0) {
                ret = -1;
                break;
            }
        }
        if (ret < 0) {
            if (errbuf && errbuf_len > 0) {
                snprintf(errbuf, errbuf_len,
                         "detect_quality_encoding: truncated FASTQ sequence block in %s", path);
            }
            ret = -1;
            break;
        }

        seq_len = seq.l;
        do {
            ret = hts_getline(fp, KS_SEP_LINE, &line);
            if (ret < 0) break;
            qual_len += line.l;
            update_ascii_range(out, line.s, line.l);
        } while (qual_len < seq_len);
        if (ret < 0 || qual_len != seq_len) {
            if (errbuf && errbuf_len > 0) {
                snprintf(errbuf, errbuf_len,
                         "detect_quality_encoding: truncated FASTQ quality block in %s", path);
            }
            ret = -1;
            break;
        }

        out->records_sampled++;
    }

    ks_free(&line);
    ks_free(&seq);
    hts_close(fp);

    if (ret == -1 && !reached_eof) {
        if (errbuf && errbuf_len > 0 && errbuf[0] == '\0') {
            snprintf(errbuf, errbuf_len,
                     "detect_quality_encoding: failed while parsing %s", path);
        }
        return -1;
    }
    finalize_quality_detection(out);
    return 0;
}
