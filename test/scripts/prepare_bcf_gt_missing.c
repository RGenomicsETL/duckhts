/* Reconstruct the GT missing-sentinel BCF from its VCF source.
 * bcf_update_genotypes accepts int32 data and compresses these sentinel values
 * to the smallest BCF integer width. See bcf_gt_missing.test for commands. */
#include <stdio.h>
#include <stdint.h>
#include <htslib/vcf.h>

int main(int argc, char **argv) {
    htsFile *input = NULL;
    htsFile *output = NULL;
    bcf_hdr_t *header = NULL;
    bcf1_t *record = NULL;
    int result = 1;
    int status;
    const int32_t genotypes[][4] = {
        {bcf_int32_missing, bcf_int32_vector_end, bcf_gt_unphased(0), bcf_int32_missing},
        {bcf_int32_missing, bcf_gt_unphased(1), bcf_gt_unphased(0), bcf_gt_phased(1)},
        {bcf_gt_missing, bcf_int32_vector_end, bcf_gt_missing, bcf_gt_phased(-1)}
    };
    size_t row = 0;

    if (argc != 3) {
        fprintf(stderr, "usage: %s input.vcf output.bcf\n", argv[0]);
        return 1;
    }
    input = bcf_open(argv[1], "r");
    if (!input) goto cleanup;
    header = bcf_hdr_read(input);
    record = bcf_init();
    if (!header || !record) goto cleanup;
    output = bcf_open(argv[2], "wb");
    if (!output || bcf_hdr_write(output, header) < 0) goto cleanup;
    while ((status = bcf_read(input, header, record)) == 0) {
        if (row >= sizeof(genotypes) / sizeof(genotypes[0])) goto cleanup;
        if (bcf_update_genotypes(header, record, genotypes[row], 4) < 0) goto cleanup;
        if (bcf_write(output, header, record) < 0) goto cleanup;
        row++;
    }
    if (status != -1 || row != sizeof(genotypes) / sizeof(genotypes[0])) goto cleanup;
    result = 0;

cleanup:
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    if (input && bcf_close(input) < 0) result = 1;
    if (output && bcf_close(output) < 0) result = 1;
    return result;
}
