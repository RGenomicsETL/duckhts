#include "filter.h"
#include <htslib/vcf.h>
#include <stdio.h>
#include <string.h>

static int expect_recovery(filter_t *filter, bcf1_t *rec, const char *message) {
    int result = filter_test(filter, rec, NULL);
    int status = filter_status(filter);
    const char *error = filter_last_error(filter);
    if (result == -1 && status == FILTER_ERR_OTHER && error && strstr(error, message)) return 1;
    fprintf(stderr, "result=%d status=%d expected=%s error=%s\n",
            result, status, message, error ? error : "(null)");
    return 0;
}

int main(void) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf1_t *rec = bcf_init();
    filter_t *filter;
    int gt[] = {bcf_gt_unphased(0), bcf_gt_unphased(0)};
    int dp[] = {10, 11}, fmt_type, info_type, i;
    const char *alleles[] = {"A", "C"};
    if (!hdr || !rec) return 2;
    bcf_hdr_append(hdr, "##fileformat=VCFv4.2");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
    bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=2,Type=Integer,Description=\"Depth\">");
    bcf_hdr_add_sample(hdr, "S1");
    bcf_hdr_add_sample(hdr, "S2");
    bcf_hdr_add_sample(hdr, NULL);
    if (bcf_hdr_sync(hdr) != 0 || bcf_update_alleles(hdr, rec, alleles, 2) != 0 ||
        bcf_update_genotypes(hdr, rec, gt, 2) != 0 ||
        bcf_update_info_int32(hdr, rec, "DP", dp, 2) != 0) return 3;
    bcf_unpack(rec, BCF_UN_ALL);

    filter = filter_parse(hdr, "N_MISSING>0");
    if (!filter || filter_status(filter) != FILTER_OK || rec->n_fmt < 1) return 4;
    fmt_type = rec->d.fmt[0].type;
    rec->d.fmt[0].type = BCF_BT_NULL;
    if (!expect_recovery(filter, rec, "Unsupported FORMAT type")) return 5;
    filter_destroy(filter);
    rec->d.fmt[0].type = fmt_type;

    filter = filter_parse(hdr, "INFO/DP[0]>0");
    if (!filter || filter_status(filter) != FILTER_OK || rec->n_info < 1) return 6;
    info_type = rec->d.info[0].type;
    rec->d.info[0].type = BCF_BT_NULL;
    if (!expect_recovery(filter, rec, "Unsupported INFO type")) return 7;
    filter_destroy(filter);
    rec->d.info[0].type = info_type;

    for (i = 0; i < 1600; i++) {
        const char *expr = i % 4 == 0 ? "INFO/DP>0 & (" :
                           i % 4 == 1 ? "ID~\"[\"" :
                           i % 4 == 2 ? "FORMAT/DP[0-2]" : "FISHER(FORMAT/DP[0-2])";
        const char *message = i % 4 == 0 ? "Could not parse" :
                              i % 4 == 1 ? "Could not compile" : "sample index is too large";
        filter = filter_parse(hdr, expr);
        if (!filter || filter_status(filter) != FILTER_ERR_OTHER ||
            !filter_last_error(filter) || !strstr(filter_last_error(filter), message)) return 8;
        filter_destroy(filter);
    }

    puts("bcftools filter recovery: OK");
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    return 0;
}
