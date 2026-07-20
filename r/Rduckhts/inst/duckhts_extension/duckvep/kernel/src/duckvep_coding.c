/*
 * duckvep_coding.c — sequence-backed coding SNV codon edit application.
 * See duckvep_coding.h. No allocation; no DuckDB/htslib.
 */
#include "duckvep_coding.h"
#include "duckvep_dna.h"

#include <string.h>

static char norm_dna_base(char c) {
    return duckvep_dna_normalize(c, 0);
}

static char norm_seq_base_or_n(uint8_t b) {
    return duckvep_dna_normalize((char)b, 1);
}

static char complement_base(char b) {
    return duckvep_dna_complement(b);
}

static char orient_genomic_base(char genomic_base, int8_t transcript_strand) {
    char b = norm_dna_base(genomic_base);
    if (b == '\0') return '\0';
    return transcript_strand < 0 ? complement_base(b) : b;
}

static void init_result(duckvep_coding_snv_result_t *out) {
    if (out == NULL) return;
    memset(out, 0, sizeof *out);
    out->aa_ref = 'X';
    out->aa_alt = 'X';
    out->change = (uint32_t)DUCKVEP_CODON_INVALID;
}

static duckvep_coding_snv_status_t fail_result(duckvep_coding_snv_result_t *out,
                                               duckvep_coding_snv_status_t status) {
    init_result(out);
    return status;
}

duckvep_coding_snv_status_t duckvep_coding_snv_from_cds(
    const uint8_t                      *cds_seq,
    size_t                              cds_len,
    const duckvep_coding_projection_t  *projection,
    char                                genomic_ref,
    char                                genomic_alt,
    int8_t                              transcript_strand,
    duckvep_codon_table_t               table,
    duckvep_coding_snv_result_t        *out) {

    size_t codon_idx;
    size_t edit_idx;
    size_t i;
    char ref_tx;
    char alt_tx;
    char seq_ref;
    duckvep_codon_result_t cr;
    uint64_t expected_cds_pos;
    uint64_t expected_codon_start;
    uint32_t expected_protein_pos;

    init_result(out);
    if (cds_seq == NULL || projection == NULL || out == NULL) {
        return DUCKVEP_CODING_SNV_INVALID_ARG;
    }
    if (projection->codon_offset > 2u || projection->codon_start_cds == 0u ||
        projection->cds_pos == 0u || projection->protein_pos == 0u ||
        (transcript_strand != (int8_t)1 && transcript_strand != (int8_t)-1)) {
        return fail_result(out, DUCKVEP_CODING_SNV_INVALID_ARG);
    }
    expected_cds_pos = (uint64_t)projection->codon_start_cds + (uint64_t)projection->codon_offset;
    expected_codon_start = ((uint64_t)projection->protein_pos - 1u) * 3u + 1u;
    expected_protein_pos = (uint32_t)((projection->cds_pos - 1u) / 3u + 1u);
    if (expected_cds_pos > UINT32_MAX || expected_codon_start > UINT32_MAX ||
        projection->cds_pos != (uint32_t)expected_cds_pos ||
        projection->codon_start_cds != (uint32_t)expected_codon_start ||
        projection->protein_pos != expected_protein_pos) {
        return fail_result(out, DUCKVEP_CODING_SNV_INVALID_ARG);
    }

    codon_idx = (size_t)projection->codon_start_cds - 1u;
    edit_idx = (size_t)projection->cds_pos - 1u;
    if (codon_idx > cds_len || cds_len - codon_idx < 3u || edit_idx >= cds_len) {
        return fail_result(out, DUCKVEP_CODING_SNV_CODON_OUT_OF_RANGE);
    }

    ref_tx = orient_genomic_base(genomic_ref, transcript_strand);
    alt_tx = orient_genomic_base(genomic_alt, transcript_strand);
    if (ref_tx == '\0' || alt_tx == '\0') return fail_result(out, DUCKVEP_CODING_SNV_INVALID_BASE);

    for (i = 0u; i < 3u; i++) {
        char b = norm_seq_base_or_n(cds_seq[codon_idx + i]);
        if (b == '\0') return fail_result(out, DUCKVEP_CODING_SNV_INVALID_BASE);
        out->ref_codon[i] = b;
        out->alt_codon[i] = b;
    }
    out->ref_codon[3] = '\0';
    out->alt_codon[3] = '\0';
    out->ref_base_tx = ref_tx;
    out->alt_base_tx = alt_tx;
    out->cds_pos = projection->cds_pos;
    out->codon_start_cds = projection->codon_start_cds;
    out->protein_pos = projection->protein_pos;
    out->codon_offset = projection->codon_offset;

    seq_ref = out->ref_codon[projection->codon_offset];
    if (seq_ref == 'N') return fail_result(out, DUCKVEP_CODING_SNV_INVALID_BASE);
    if (seq_ref != ref_tx) return fail_result(out, DUCKVEP_CODING_SNV_REF_MISMATCH);

    out->alt_codon[projection->codon_offset] = alt_tx;
    cr = duckvep_codon_change_prepared(
        out->ref_codon, out->alt_codon, table);
    out->aa_ref = cr.aa_ref;
    out->aa_alt = cr.aa_alt;
    out->change = cr.change;
    return DUCKVEP_CODING_SNV_OK;
}
