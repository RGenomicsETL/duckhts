/*
 * duckvep_delta.c — sequence-delta producers. See duckvep_delta.h.
 *
 * The codon/peptide FACTS layer of the engine. The annotation sweep remains
 * orchestration-only and all coding edit classification enters here. No
 * DuckDB/htslib/Arrow; no allocation. Builds on the oracle-tested coding kernels
 * (duckvep_project_coding_base + duckvep_coding_snv_from_cds).
 */
#include "duckvep_delta.h"

#include "duckvep_coding.h" /* pulls duckvep_projection.h + duckvep_codon.h */
#include "duckvep_dna.h"
#include "duckvep_event.h"

#include <string.h>

static char delta_norm_base(char c) {
    return duckvep_dna_normalize(c, 1);
}

static int delta_norm_span_unambiguous(const uint8_t *seq, size_t len) {
    size_t i;

    if (len == 0u) return 1;
    if (seq == NULL) return 0;
    for (i = 0u; i < len; i++) {
        char b = delta_norm_base((char)seq[i]);
        if (b != 'A' && b != 'C' && b != 'G' && b != 'T') return 0;
    }
    return 1;
}

static int delta_peptide_span_known_nonstop(const uint8_t *pep, size_t off, size_t len) {
    size_t i;

    if (len == 0u) return 1;
    if (pep == NULL) return 0;
    for (i = 0u; i < len; i++) {
        if (pep[off + i] == (uint8_t)'\0' || pep[off + i] == (uint8_t)'X' ||
            pep[off + i] == (uint8_t)'*') return 0;
    }
    return 1;
}

static int delta_peptide_span_known(const uint8_t *pep, size_t off, size_t len) {
    size_t i;

    if (len == 0u) return 1;
    if (pep == NULL) return 0;
    for (i = 0u; i < len; i++) {
        if (pep[off + i] == (uint8_t)'\0' || pep[off + i] == (uint8_t)'X') {
            return 0;
        }
    }
    return 1;
}

static int delta_peptide_stop_index(
    const uint8_t *pep, size_t off, size_t len, size_t *index) {

    size_t i;

    if (index != NULL) *index = 0u;
    if (pep == NULL || index == NULL) return 0;
    for (i = 0u; i < len; i++) {
        if (pep[off + i] == (uint8_t)'*') {
            *index = i;
            return 1;
        }
    }
    return 0;
}

static int delta_norm_span_equal(const uint8_t *a, const uint8_t *b, size_t len) {
    size_t i;

    if (len == 0u) return 1;
    if (a == NULL || b == NULL) return 0;
    for (i = 0u; i < len; i++) {
        char ab = delta_norm_base((char)a[i]);
        char bb = delta_norm_base((char)b[i]);
        if (ab == '\0' || bb == '\0' || ab != bb) return 0;
    }
    return 1;
}

static char delta_complement_base(char b) {
    return duckvep_dna_complement(b);
}

static char delta_orient_genomic_base(char genomic_base, int8_t transcript_strand) {
    char b = delta_norm_base(genomic_base);
    if (b == 'N') return '\0';
    if (b == '\0') return '\0';
    return transcript_strand < 0 ? delta_complement_base(b) : b;
}

typedef struct delta_edit_view {
    duckvep_event_t event;
    size_t          ref_off;
    size_t          alt_off;
    size_t          anchor_ref_off;
    uint16_t        ref_len;
    uint16_t        alt_len;
} delta_edit_view_t;

static int delta_allele_slice_ok(const duckvep_variant_batch_t *v, uint32_t variant_idx) {
    uint64_t pool;
    uint64_t roff, rlen, aoff, alen;
    if (v == NULL || v->allele_bytes == NULL || v->ref_offset == NULL || v->alt_offset == NULL ||
        v->ref_length == NULL || v->alt_length == NULL) return 0;
    pool = (uint64_t)v->allele_bytes_len;
    roff = (uint64_t)v->ref_offset[variant_idx]; rlen = (uint64_t)v->ref_length[variant_idx];
    aoff = (uint64_t)v->alt_offset[variant_idx]; alen = (uint64_t)v->alt_length[variant_idx];
    return roff <= pool && rlen <= pool - roff && aoff <= pool && alen <= pool - aoff &&
           rlen <= UINT16_MAX && alen <= UINT16_MAX;
}

static int delta_edit_view_load(
    const duckvep_variant_batch_t *v,
    uint32_t                       variant_idx,
    const duckvep_event_t         *prepared_event,
    delta_edit_view_t             *edit) {

    if (!delta_allele_slice_ok(v, variant_idx)) return 0;
    if (prepared_event != NULL) edit->event = *prepared_event;
    else duckvep_event_load(v, (size_t)variant_idx, &edit->event);
    edit->ref_off = (size_t)v->ref_offset[variant_idx] +
                    (size_t)edit->event.ref_diff_offset;
    edit->alt_off = (size_t)v->alt_offset[variant_idx] +
                    (size_t)edit->event.alt_diff_offset;
    edit->ref_len = edit->event.ref_diff_length;
    edit->alt_len = edit->event.alt_diff_length;
    edit->anchor_ref_off = (size_t)v->ref_offset[variant_idx] +
                           (size_t)edit->event.anchor_ref_offset;
    if (edit->event.interbase) {
        if (edit->event.anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT &&
            edit->event.anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_RIGHT) return 0;
    }
    return edit->ref_off <= v->allele_bytes_len &&
           (size_t)edit->ref_len <= v->allele_bytes_len - edit->ref_off &&
           edit->alt_off <= v->allele_bytes_len &&
           (size_t)edit->alt_len <= v->allele_bytes_len - edit->alt_off &&
           edit->anchor_ref_off < v->allele_bytes_len;
}

static int delta_cds_ref_matches(const uint8_t *cds_seq, size_t cds_len,
                                 uint32_t cds_pos, char genomic_ref,
                                 int8_t transcript_strand) {
    char ref_tx;
    char cds_base;
    if (cds_pos == 0u || (size_t)cds_pos > cds_len) return 0;
    ref_tx = delta_orient_genomic_base(genomic_ref, transcript_strand);
    cds_base = delta_norm_base((char)cds_seq[(size_t)cds_pos - 1u]);
    if (ref_tx == '\0' || cds_base == '\0' || cds_base == 'N') return 0;
    return ref_tx == cds_base;
}

static int delta_cds_slice(
    const duckvep_sequence_pool_t *seq,
    size_t                         tx_idx,
    const uint8_t                **cds_seq,
    size_t                        *cds_len) {

    uint64_t off;
    uint64_t len;
    uint64_t pool;

    if (cds_seq != NULL) *cds_seq = NULL;
    if (cds_len != NULL) *cds_len = 0u;
    if (seq == NULL || cds_seq == NULL || cds_len == NULL ||
        seq->cds_bytes == NULL || seq->cds_offset == NULL || seq->cds_length == NULL ||
        tx_idx >= seq->transcript_count) {
        return 0;
    }
    off = seq->cds_offset[tx_idx];
    len = seq->cds_length[tx_idx];
    pool = (uint64_t)seq->cds_bytes_len;
    if (off > pool || len > pool - off) return 0;
    *cds_seq = seq->cds_bytes + (size_t)off;
    *cds_len = (size_t)len;
    return 1;
}

static int delta_allele_bases_valid(
    const duckvep_variant_batch_t *v,
    size_t                         off,
    uint16_t                       len) {

    uint16_t i;
    if (len == 0u) return 1;
    if (v == NULL || v->allele_bytes == NULL || off > v->allele_bytes_len ||
        (size_t)len > v->allele_bytes_len - off) {
        return 0;
    }
    for (i = 0u; i < len; i++) {
        char b = delta_norm_base((char)v->allele_bytes[off + (size_t)i]);
        if (b == '\0' || b == 'N') return 0;
    }
    return 1;
}

static duckvep_cds_edit_status_t duckvep_variant_cds_edit_build_event(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *prepared_event,
    duckvep_haplotype_edit_t         *out) {

    delta_edit_view_t edit;
    const uint8_t *cds_seq;
    size_t cds_len;
    duckvep_coding_projection_t proj;
    uint32_t min_cds;
    uint32_t max_cds;
    uint32_t projected_pos;
    uint16_t j;
    uint8_t kind;

    if (out != NULL) memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || seq == NULL || v == NULL || out == NULL ||
        v->variant_kind == NULL || variant_idx >= v->count ||
        tx_idx >= transcripts->transcript_count ||
        (transcript_strand != (int8_t)1 && transcript_strand != (int8_t)-1)) {
        return DUCKVEP_CDS_EDIT_INVALID_ARG;
    }
    kind = v->variant_kind[variant_idx];
    if (kind != (uint8_t)DUCKVEP_KIND_SNV &&
        kind != (uint8_t)DUCKVEP_KIND_INS &&
        kind != (uint8_t)DUCKVEP_KIND_DEL &&
        kind != (uint8_t)DUCKVEP_KIND_INDEL &&
        kind != (uint8_t)DUCKVEP_KIND_MNV) {
        return DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND;
    }
    if (!delta_cds_slice(seq, tx_idx, &cds_seq, &cds_len) || cds_len > UINT32_MAX) {
        return DUCKVEP_CDS_EDIT_INVALID_ARG;
    }
    if (!delta_edit_view_load(v, variant_idx, prepared_event, &edit)) {
        return DUCKVEP_CDS_EDIT_INVALID_EVENT;
    }
    if (edit.ref_len == 0u && edit.alt_len == 0u) return DUCKVEP_CDS_EDIT_INVALID_EVENT;
    if (!delta_allele_bases_valid(v, edit.alt_off, edit.alt_len)) {
        return DUCKVEP_CDS_EDIT_INVALID_ALLELE;
    }

    out->variant_strand = (int8_t)1;
    out->alt_len = (uint32_t)edit.alt_len;
    out->alt = edit.alt_len > 0u ? v->allele_bytes + edit.alt_off : NULL;

    if (edit.event.interbase) {
        uint32_t boundary;
        uint32_t right_flank;

        if (edit.ref_len != 0u || edit.alt_len == 0u) return DUCKVEP_CDS_EDIT_INVALID_EVENT;
        boundary = edit.event.insertion_boundary0;
        right_flank = duckvep_event_right_flank1(&edit.event);
        projected_pos = edit.event.start1;
        if (duckvep_project_coding_base(transcripts, exons, tx_idx,
                                        projected_pos, &proj)) {
            /* The retained VCF padding base is usable only when it is itself in
             * the CDS. Genomic REF validation outside the CDS belongs upstream. */
            if (!delta_cds_ref_matches(cds_seq, cds_len, proj.cds_pos,
                                       (char)v->allele_bytes[edit.anchor_ref_off],
                                       transcript_strand)) {
                return DUCKVEP_CDS_EDIT_REF_MISMATCH;
            }
        } else {
            projected_pos = edit.event.anchor_side ==
                                (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT
                              ? right_flank
                              : boundary;
            if (projected_pos == 0u ||
                !duckvep_project_coding_base(transcripts, exons, tx_idx,
                                             projected_pos, &proj)) {
                return DUCKVEP_CDS_EDIT_OUT_OF_CDS;
            }
        }

        /* cds_start is the 1-based position before which the inserted bases are
         * written. The same genomic boundary reverses direction on a minus-strand
         * transcript, hence the two flank-specific formulas. */
        if (projected_pos == boundary) {
            out->cds_start = transcript_strand > 0 ? proj.cds_pos + 1u
                                                   : proj.cds_pos;
        } else if (projected_pos == right_flank) {
            out->cds_start = transcript_strand > 0 ? proj.cds_pos
                                                   : proj.cds_pos + 1u;
        } else {
            return DUCKVEP_CDS_EDIT_INVALID_EVENT;
        }
        if (out->cds_start == 0u || (uint64_t)out->cds_start > (uint64_t)cds_len + 1u) {
            return DUCKVEP_CDS_EDIT_OUT_OF_CDS;
        }
        out->ref_len = 0u;
        out->ref = NULL;
        return DUCKVEP_CDS_EDIT_OK;
    }

    if (edit.ref_len == 0u) return DUCKVEP_CDS_EDIT_INVALID_EVENT;
    if (!delta_allele_bases_valid(v, edit.ref_off, edit.ref_len)) {
        return DUCKVEP_CDS_EDIT_INVALID_ALLELE;
    }

    min_cds = UINT32_MAX;
    max_cds = 0u;
    for (j = 0u; j < edit.ref_len; j++) {
        uint32_t gpos = edit.event.start1 + (uint32_t)j;
        if (!duckvep_project_coding_base(transcripts, exons, tx_idx, gpos, &proj)) {
            return DUCKVEP_CDS_EDIT_OUT_OF_CDS;
        }
        if (!delta_cds_ref_matches(cds_seq, cds_len, proj.cds_pos,
                                   (char)v->allele_bytes[edit.ref_off + (size_t)j],
                                   transcript_strand)) {
            return DUCKVEP_CDS_EDIT_REF_MISMATCH;
        }
        if (proj.cds_pos < min_cds) min_cds = proj.cds_pos;
        if (proj.cds_pos > max_cds) max_cds = proj.cds_pos;
    }
    if (min_cds == UINT32_MAX || max_cds < min_cds ||
        max_cds - min_cds + 1u != (uint32_t)edit.ref_len) {
        return DUCKVEP_CDS_EDIT_NON_CONTIGUOUS;
    }
    out->cds_start = min_cds;
    out->ref_len = (uint32_t)edit.ref_len;
    out->ref = v->allele_bytes + edit.ref_off;
    return DUCKVEP_CDS_EDIT_OK;
}

DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t duckvep_variant_cds_edit_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_haplotype_edit_t         *out) {

    return duckvep_variant_cds_edit_build_event(
        transcripts, exons, seq, v, variant_idx, tx_idx, transcript_strand,
        NULL, out);
}

static duckvep_cds_edit_status_t duckvep_variant_cds_edit_set_build_event(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *prepared_event,
    duckvep_haplotype_edit_t         *scratch,
    size_t                            scratch_cap,
    duckvep_edit_set_t               *out) {

    duckvep_haplotype_edit_t edit;
    duckvep_cds_edit_status_t st;
    uint8_t kind;
    uint32_t islands = 0u;
    uint32_t i;
    uint32_t out_i;

    if (out == NULL) return DUCKVEP_CDS_EDIT_INVALID_ARG;
    out->edits = NULL;
    out->count = 0u;
    if (scratch == NULL && scratch_cap > 0u) return DUCKVEP_CDS_EDIT_INVALID_ARG;

    st = duckvep_variant_cds_edit_build_event(
        transcripts, exons, seq, v, variant_idx, tx_idx, transcript_strand,
        prepared_event, &edit);
    if (st != DUCKVEP_CDS_EDIT_OK) return st;
    if (v == NULL || v->variant_kind == NULL || variant_idx >= v->count) {
        return DUCKVEP_CDS_EDIT_INVALID_ARG;
    }
    kind = v->variant_kind[variant_idx];
    if (kind != (uint8_t)DUCKVEP_KIND_MNV || edit.ref_len != edit.alt_len ||
        edit.ref_len == 0u || edit.ref == NULL || edit.alt == NULL) {
        if (scratch_cap == 0u) return DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL;
        scratch[0] = edit;
        out->edits = scratch;
        out->count = 1u;
        return DUCKVEP_CDS_EDIT_OK;
    }

    i = 0u;
    while (i < edit.ref_len) {
        char rb = delta_norm_base((char)edit.ref[i]);
        char ab = delta_norm_base((char)edit.alt[i]);
        if (rb == '\0' || ab == '\0') return DUCKVEP_CDS_EDIT_INVALID_ALLELE;
        if (rb == ab) {
            i++;
            continue;
        }
        islands++;
        i++;
        while (i < edit.ref_len) {
            rb = delta_norm_base((char)edit.ref[i]);
            ab = delta_norm_base((char)edit.alt[i]);
            if (rb == '\0' || ab == '\0') return DUCKVEP_CDS_EDIT_INVALID_ALLELE;
            if (rb == ab) break;
            i++;
        }
    }
    if (islands == 0u) return DUCKVEP_CDS_EDIT_INVALID_EVENT;
    if (scratch_cap < (size_t)islands) return DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL;

    out_i = 0u;
    if (transcript_strand > 0) {
        uint32_t pos = edit.ref_len;
        while (pos > 0u) {
            uint32_t hi;
            uint32_t lo;
            while (pos > 0u &&
                   delta_norm_base((char)edit.ref[(size_t)pos - 1u]) ==
                   delta_norm_base((char)edit.alt[(size_t)pos - 1u])) {
                pos--;
            }
            if (pos == 0u) break;
            hi = pos;
            lo = pos - 1u;
            while (lo > 0u &&
                   delta_norm_base((char)edit.ref[(size_t)lo - 1u]) !=
                   delta_norm_base((char)edit.alt[(size_t)lo - 1u])) {
                lo--;
            }
            scratch[out_i] = edit;
            scratch[out_i].cds_start = edit.cds_start + lo;
            scratch[out_i].ref_len = hi - lo;
            scratch[out_i].alt_len = hi - lo;
            scratch[out_i].ref = edit.ref + lo;
            scratch[out_i].alt = edit.alt + lo;
            out_i++;
            pos = lo;
        }
    } else {
        uint32_t pos = 0u;
        while (pos < edit.ref_len) {
            uint32_t lo;
            uint32_t hi;
            while (pos < edit.ref_len &&
                   delta_norm_base((char)edit.ref[pos]) ==
                   delta_norm_base((char)edit.alt[pos])) {
                pos++;
            }
            if (pos >= edit.ref_len) break;
            lo = pos;
            while (pos < edit.ref_len &&
                   delta_norm_base((char)edit.ref[pos]) !=
                   delta_norm_base((char)edit.alt[pos])) {
                pos++;
            }
            hi = pos;
            scratch[out_i] = edit;
            scratch[out_i].cds_start = edit.cds_start + edit.ref_len - hi;
            scratch[out_i].ref_len = hi - lo;
            scratch[out_i].alt_len = hi - lo;
            scratch[out_i].ref = edit.ref + lo;
            scratch[out_i].alt = edit.alt + lo;
            out_i++;
        }
    }
    if (out_i != islands) return DUCKVEP_CDS_EDIT_INVALID_EVENT;
    out->edits = scratch;
    out->count = (size_t)islands;
    return DUCKVEP_CDS_EDIT_OK;
}

DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t duckvep_variant_cds_edit_set_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_haplotype_edit_t         *scratch,
    size_t                            scratch_cap,
    duckvep_edit_set_t               *out) {

    return duckvep_variant_cds_edit_set_build_event(
        transcripts, exons, seq, v, variant_idx, tx_idx, transcript_strand,
        NULL, scratch, scratch_cap, out);
}

static duckvep_coding_context_status_t delta_context_status_from_haplo(
    duckvep_haplotype_status_t st) {

    switch (st) {
    case DUCKVEP_HAPLOTYPE_OK:
        return DUCKVEP_CODING_CONTEXT_OK;
    case DUCKVEP_HAPLOTYPE_INVALID_ARG:
        return DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    case DUCKVEP_HAPLOTYPE_OUT_OF_RANGE:
        return DUCKVEP_CODING_CONTEXT_OUT_OF_RANGE;
    case DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL:
        return DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL;
    case DUCKVEP_HAPLOTYPE_INVALID_BASE:
        return DUCKVEP_CODING_CONTEXT_INVALID_BASE;
    case DUCKVEP_HAPLOTYPE_REF_MISMATCH:
        return DUCKVEP_CODING_CONTEXT_REF_MISMATCH;
    case DUCKVEP_HAPLOTYPE_EDIT_ORDER:
        return DUCKVEP_CODING_CONTEXT_EDIT_ORDER;
    default:
        return DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    }
}

static duckvep_coding_context_status_t delta_translate_cds_full(
    const uint8_t        *cds,
    size_t                cds_len,
    duckvep_codon_table_t table,
    uint8_t              *pep,
    size_t                pep_cap,
    size_t               *pep_len_out,
    duckvep_coding_context_status_t cap_status) {

    size_t codons;
    size_t i;

    if (pep_len_out != NULL) *pep_len_out = 0u;
    if (cds == NULL || pep == NULL || pep_len_out == NULL) {
        return DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    }
    codons = cds_len / 3u;
    if (pep_cap < codons + 1u) return cap_status;
    for (i = 0u; i < codons; i++) {
        char codon[4];
        uint32_t j;
        int has_n = 0;
        for (j = 0u; j < 3u; j++) {
            char b = delta_norm_base((char)cds[i * 3u + (size_t)j]);
            if (b == '\0') return DUCKVEP_CODING_CONTEXT_INVALID_BASE;
            if (b == 'N') has_n = 1;
            codon[j] = b;
        }
        codon[3] = '\0';
        pep[i] = has_n ? (uint8_t)'X'
                       : (uint8_t)duckvep_translate_codon(codon, table);
    }
    pep[codons] = (uint8_t)'\0';
    *pep_len_out = codons;
    return DUCKVEP_CODING_CONTEXT_OK;
}

static duckvep_coding_context_status_t delta_cds_changed(
    const uint8_t *ref_cds,
    size_t         ref_cds_len,
    const uint8_t *alt_cds,
    size_t         alt_cds_len,
    uint8_t       *changed) {

    size_t i;
    if (changed != NULL) *changed = 0u;
    if (ref_cds == NULL || alt_cds == NULL || changed == NULL) {
        return DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    }
    if (ref_cds_len != alt_cds_len) {
        *changed = 1u;
        return DUCKVEP_CODING_CONTEXT_OK;
    }
    for (i = 0u; i < ref_cds_len; i++) {
        char rb = delta_norm_base((char)ref_cds[i]);
        char ab = delta_norm_base((char)alt_cds[i]);
        if (rb == '\0' || ab == '\0') return DUCKVEP_CODING_CONTEXT_INVALID_BASE;
        if (rb != ab) {
            *changed = 1u;
            return DUCKVEP_CODING_CONTEXT_OK;
        }
    }
    return DUCKVEP_CODING_CONTEXT_OK;
}

static duckvep_coding_context_status_t delta_peptide_diff_window(
    const uint8_t *ref_pep,
    size_t         ref_len,
    const uint8_t *alt_pep,
    size_t         alt_len,
    uint32_t      *ref_first,
    uint32_t      *ref_last,
    uint32_t      *alt_first,
    uint32_t      *alt_last) {

    size_t prefix = 0u;
    size_t suffix = 0u;
    size_t ref_mid;
    size_t alt_mid;

    *ref_first = 0u;
    *ref_last = 0u;
    *alt_first = 0u;
    *alt_last = 0u;
    while (prefix < ref_len && prefix < alt_len && ref_pep[prefix] == alt_pep[prefix]) {
        prefix++;
    }
    if (prefix == ref_len && prefix == alt_len) return DUCKVEP_CODING_CONTEXT_OK;
    while (suffix < ref_len - prefix && suffix < alt_len - prefix &&
           ref_pep[ref_len - 1u - suffix] == alt_pep[alt_len - 1u - suffix]) {
        suffix++;
    }
    ref_mid = ref_len - prefix - suffix;
    alt_mid = alt_len - prefix - suffix;
    if (ref_mid > 0u) {
        if (prefix >= (size_t)UINT32_MAX || ref_len - suffix > (size_t)UINT32_MAX) {
            return DUCKVEP_CODING_CONTEXT_OUT_OF_RANGE;
        }
        *ref_first = (uint32_t)prefix + 1u;
        *ref_last = (uint32_t)(ref_len - suffix);
    }
    if (alt_mid > 0u) {
        if (prefix >= (size_t)UINT32_MAX || alt_len - suffix > (size_t)UINT32_MAX) {
            return DUCKVEP_CODING_CONTEXT_OUT_OF_RANGE;
        }
        *alt_first = (uint32_t)prefix + 1u;
        *alt_last = (uint32_t)(alt_len - suffix);
    }
    return DUCKVEP_CODING_CONTEXT_OK;
}

DUCKVEP_INTERNAL_API duckvep_coding_context_status_t duckvep_coding_context_build(
    const uint8_t               *ref_cds,
    size_t                       ref_cds_len,
    const duckvep_edit_set_t    *edit_set,
    int8_t                       transcript_strand,
    duckvep_codon_table_t        table,
    uint8_t                     *alt_cds_scratch,
    size_t                       alt_cds_cap,
    uint8_t                     *ref_peptide_scratch,
    size_t                       ref_peptide_cap,
    uint8_t                     *alt_peptide_scratch,
    size_t                       alt_peptide_cap,
    duckvep_coding_context_t    *ctx) {

    duckvep_coding_context_t tmp;
    duckvep_haplotype_result_t apply_result;
    duckvep_haplotype_status_t hst;
    duckvep_coding_context_status_t cst;
    size_t alt_cds_len = 0u;
    size_t ref_pep_len = 0u;
    size_t alt_pep_len = 0u;
    uint8_t cds_changed = 0u;

    if (ctx == NULL) return DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    memset(ctx, 0, sizeof *ctx);
    if (ref_cds == NULL || edit_set == NULL || alt_cds_scratch == NULL ||
        ref_peptide_scratch == NULL || alt_peptide_scratch == NULL ||
        (edit_set->count > 0u && edit_set->edits == NULL) ||
        (transcript_strand != (int8_t)1 && transcript_strand != (int8_t)-1)) {
        return DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    }

    hst = duckvep_haplotype_apply_cds_edits(ref_cds, ref_cds_len, edit_set->edits,
                                            edit_set->count, transcript_strand,
                                            alt_cds_scratch, alt_cds_cap,
                                            &alt_cds_len, &apply_result);
    if (hst != DUCKVEP_HAPLOTYPE_OK) return delta_context_status_from_haplo(hst);

    cst = delta_translate_cds_full(ref_cds, ref_cds_len, table, ref_peptide_scratch,
                                   ref_peptide_cap, &ref_pep_len,
                                   DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL);
    if (cst != DUCKVEP_CODING_CONTEXT_OK) return cst;
    cst = delta_translate_cds_full(alt_cds_scratch, alt_cds_len, table,
                                   alt_peptide_scratch, alt_peptide_cap,
                                   &alt_pep_len,
                                   DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL);
    if (cst != DUCKVEP_CODING_CONTEXT_OK) return cst;
    cst = delta_cds_changed(ref_cds, ref_cds_len, alt_cds_scratch, alt_cds_len,
                            &cds_changed);
    if (cst != DUCKVEP_CODING_CONTEXT_OK) return cst;

    memset(&tmp, 0, sizeof tmp);
    tmp.ref_cds = ref_cds;
    tmp.ref_cds_len = ref_cds_len;
    tmp.alt_cds = alt_cds_scratch;
    tmp.alt_cds_len = alt_cds_len;
    tmp.ref_peptide = ref_peptide_scratch;
    tmp.ref_peptide_len = ref_pep_len;
    tmp.alt_peptide = alt_peptide_scratch;
    tmp.alt_peptide_len = alt_pep_len;
    tmp.codon_table = (uint8_t)table;
    tmp.length_diff = apply_result.length_diff;
    tmp.flags = apply_result.flags;
    tmp.applied_edits = apply_result.applied_edits;
    if (edit_set->count == 1u) {
        tmp.has_single_edit = 1u;
        tmp.single_edit_cds_start = edit_set->edits[0].cds_start;
        tmp.single_edit_ref_len = edit_set->edits[0].ref_len;
        tmp.single_edit_alt_len = edit_set->edits[0].alt_len;
    }
    tmp.cds_changed = cds_changed;
    cst = delta_peptide_diff_window(ref_peptide_scratch, ref_pep_len,
                                    alt_peptide_scratch, alt_pep_len,
                                    &tmp.ref_first_changed_codon,
                                    &tmp.ref_last_changed_codon,
                                    &tmp.alt_first_changed_codon,
                                    &tmp.alt_last_changed_codon);
    if (cst != DUCKVEP_CODING_CONTEXT_OK) return cst;
    *ctx = tmp;
    return DUCKVEP_CODING_CONTEXT_OK;
}

static duckvep_variant_coding_context_status_t delta_variant_context_from_edit_status(
    duckvep_cds_edit_status_t st) {

    switch (st) {
    case DUCKVEP_CDS_EDIT_OK:
        return DUCKVEP_VARIANT_CODING_CONTEXT_OK;
    case DUCKVEP_CDS_EDIT_INVALID_ARG:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    case DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND:
        return DUCKVEP_VARIANT_CODING_CONTEXT_UNSUPPORTED_KIND;
    case DUCKVEP_CDS_EDIT_INVALID_EVENT:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_EVENT;
    case DUCKVEP_CDS_EDIT_OUT_OF_CDS:
        return DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_CDS;
    case DUCKVEP_CDS_EDIT_NON_CONTIGUOUS:
        return DUCKVEP_VARIANT_CODING_CONTEXT_NON_CONTIGUOUS;
    case DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL:
        return DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL;
    case DUCKVEP_CDS_EDIT_INVALID_ALLELE:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ALLELE;
    case DUCKVEP_CDS_EDIT_REF_MISMATCH:
        return DUCKVEP_VARIANT_CODING_CONTEXT_REF_MISMATCH;
    default:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    }
}

static duckvep_variant_coding_context_status_t delta_variant_context_from_context_status(
    duckvep_coding_context_status_t st) {

    switch (st) {
    case DUCKVEP_CODING_CONTEXT_OK:
        return DUCKVEP_VARIANT_CODING_CONTEXT_OK;
    case DUCKVEP_CODING_CONTEXT_INVALID_ARG:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    case DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL:
        return DUCKVEP_VARIANT_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL;
    case DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL:
        return DUCKVEP_VARIANT_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL;
    case DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL:
        return DUCKVEP_VARIANT_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL;
    case DUCKVEP_CODING_CONTEXT_OUT_OF_RANGE:
        return DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE;
    case DUCKVEP_CODING_CONTEXT_INVALID_BASE:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_BASE;
    case DUCKVEP_CODING_CONTEXT_REF_MISMATCH:
        return DUCKVEP_VARIANT_CODING_CONTEXT_REF_MISMATCH;
    case DUCKVEP_CODING_CONTEXT_EDIT_ORDER:
        return DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_ORDER;
    default:
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    }
}

static uint8_t delta_sequence_status_from_context(
    duckvep_variant_coding_context_status_t status) {

    switch (status) {
    case DUCKVEP_VARIANT_CODING_CONTEXT_OK:
        return (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ALLELE:
    case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_BASE:
        return (uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS;
    case DUCKVEP_VARIANT_CODING_CONTEXT_REF_MISMATCH:
        return (uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH;
    case DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_CDS:
    case DUCKVEP_VARIANT_CODING_CONTEXT_NON_CONTIGUOUS:
        return (uint8_t)DUCKVEP_SEQUENCE_NON_CONTIGUOUS_EDIT;
    case DUCKVEP_VARIANT_CODING_CONTEXT_UNSUPPORTED_KIND:
    case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_EVENT:
        return (uint8_t)DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT;
    case DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL:
    case DUCKVEP_VARIANT_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL:
    case DUCKVEP_VARIANT_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL:
    case DUCKVEP_VARIANT_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL:
        return (uint8_t)DUCKVEP_SEQUENCE_INTERNAL_CAPACITY;
    case DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG:
    case DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE:
    case DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_ORDER:
    default:
        return (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
    }
}

static uint8_t delta_sequence_status_from_delta(
    duckvep_context_delta_status_t status) {

    switch (status) {
    case DUCKVEP_CONTEXT_DELTA_OK:
        return (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    case DUCKVEP_CONTEXT_DELTA_UNSUPPORTED:
        return (uint8_t)DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT;
    case DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL:
        return (uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL;
    case DUCKVEP_CONTEXT_DELTA_INVALID_ARG:
    default:
        return (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
    }
}

static uint8_t delta_sequence_status_from_snv(
    duckvep_coding_snv_status_t status) {

    switch (status) {
    case DUCKVEP_CODING_SNV_OK:
        return (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    case DUCKVEP_CODING_SNV_INVALID_BASE:
        return (uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS;
    case DUCKVEP_CODING_SNV_REF_MISMATCH:
        return (uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH;
    case DUCKVEP_CODING_SNV_INVALID_ARG:
    case DUCKVEP_CODING_SNV_CODON_OUT_OF_RANGE:
    default:
        return (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
    }
}

static int delta_transcript_has_sequence(
    const duckvep_sequence_pool_t *seq,
    size_t                         tx_idx) {

    return seq != NULL && seq->cds_length != NULL &&
           tx_idx < seq->transcript_count && seq->cds_length[tx_idx] > 0u;
}

static duckvep_variant_coding_context_status_t duckvep_variant_coding_context_build_event(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *prepared_event,
    duckvep_haplotype_edit_t         *edit_scratch,
    size_t                            edit_scratch_cap,
    uint8_t                          *alt_cds_scratch,
    size_t                            alt_cds_cap,
    uint8_t                          *ref_peptide_scratch,
    size_t                            ref_peptide_cap,
    uint8_t                          *alt_peptide_scratch,
    size_t                            alt_peptide_cap,
    duckvep_coding_context_t         *ctx) {

    const uint8_t *ref_cds;
    size_t ref_cds_len;
    duckvep_edit_set_t edit_set;
    duckvep_cds_edit_status_t est;
    duckvep_coding_context_status_t cst;

    if (ctx == NULL) return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    memset(ctx, 0, sizeof *ctx);
    if (seq == NULL || seq->codon_table == NULL || tx_idx >= seq->transcript_count) {
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    }
    if (!delta_cds_slice(seq, tx_idx, &ref_cds, &ref_cds_len)) {
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    }
    est = duckvep_variant_cds_edit_set_build_event(
        transcripts, exons, seq, v, variant_idx, tx_idx, transcript_strand,
        prepared_event, edit_scratch, edit_scratch_cap, &edit_set);
    if (est != DUCKVEP_CDS_EDIT_OK) return delta_variant_context_from_edit_status(est);
    cst = duckvep_coding_context_build(ref_cds, ref_cds_len, &edit_set,
                                       transcript_strand,
                                       (duckvep_codon_table_t)seq->codon_table[tx_idx],
                                       alt_cds_scratch, alt_cds_cap,
                                       ref_peptide_scratch, ref_peptide_cap,
                                       alt_peptide_scratch, alt_peptide_cap, ctx);
    if (cst != DUCKVEP_CODING_CONTEXT_OK) {
        return delta_variant_context_from_context_status(cst);
    }
    if (seq->post_cds_bases != NULL) {
        size_t base;
        uint8_t length = 0u;

        if (tx_idx >
            (SIZE_MAX - (DUCKVEP_POST_CDS_BASE_COUNT - 1u)) /
                DUCKVEP_POST_CDS_BASE_COUNT) {
            return DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE;
        }
        base = tx_idx * DUCKVEP_POST_CDS_BASE_COUNT;
        while (length < DUCKVEP_POST_CDS_BASE_COUNT &&
               seq->post_cds_bases[base + length] != 0u) {
            length++;
        }
        ctx->post_cds_bases = seq->post_cds_bases + base;
        ctx->post_cds_length = length;
    }
    return DUCKVEP_VARIANT_CODING_CONTEXT_OK;
}

DUCKVEP_INTERNAL_API duckvep_variant_coding_context_status_t duckvep_variant_coding_context_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_haplotype_edit_t         *edit_scratch,
    size_t                            edit_scratch_cap,
    uint8_t                          *alt_cds_scratch,
    size_t                            alt_cds_cap,
    uint8_t                          *ref_peptide_scratch,
    size_t                            ref_peptide_cap,
    uint8_t                          *alt_peptide_scratch,
    size_t                            alt_peptide_cap,
    duckvep_coding_context_t         *ctx) {

    duckvep_event_t event;

    if (v == NULL || variant_idx >= v->count) {
        if (ctx != NULL) memset(ctx, 0, sizeof *ctx);
        return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
    }
    duckvep_event_load(v, variant_idx, &event);

    return duckvep_variant_coding_context_build_event(
        transcripts, exons, seq, v, variant_idx, tx_idx, transcript_strand,
        &event, edit_scratch, edit_scratch_cap, alt_cds_scratch, alt_cds_cap,
        ref_peptide_scratch, ref_peptide_cap, alt_peptide_scratch,
        alt_peptide_cap, ctx);
}

static duckvep_context_delta_status_t delta_context_inframe_deletion(
    const duckvep_coding_context_t *ctx,
    duckvep_sequence_delta_t       *delta) {

    uint64_t deleted_len;
    uint64_t deleted_end;
    size_t prefix_len;
    size_t suffix_len;
    uint32_t protein_cds;

    if (ctx == NULL || delta == NULL) return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    if (ctx->ref_cds == NULL || ctx->alt_cds == NULL || ctx->ref_peptide == NULL ||
        ctx->alt_peptide == NULL || ctx->ref_cds_len == 0u ||
        (ctx->ref_cds_len % 3u) != 0u || !ctx->has_single_edit ||
        ctx->applied_edits != 1u || ctx->length_diff >= 0 ||
        ctx->length_diff == INT64_MIN || ctx->single_edit_alt_len != 0u ||
        ctx->single_edit_ref_len == 0u || ctx->cds_changed == 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    deleted_len = (uint64_t)ctx->single_edit_ref_len;
    if (deleted_len != (uint64_t)(-(ctx->length_diff + 1)) + 1u ||
        (deleted_len % 3u) != 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (deleted_len > ctx->ref_cds_len || ctx->alt_cds_len != ctx->ref_cds_len - (size_t)deleted_len ||
        ctx->ref_peptide_len != ctx->ref_cds_len / 3u ||
        ctx->alt_peptide_len != ctx->alt_cds_len / 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    /* A pure in-frame deletion is inframe_deletion in VEP regardless of codon alignment.
     * The start-codon predicate is layered by the caller after this shape is known. */
    deleted_end = (uint64_t)ctx->single_edit_cds_start + deleted_len - 1u;
    if (deleted_end > (uint64_t)ctx->ref_cds_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    prefix_len = (size_t)ctx->single_edit_cds_start - 1u;
    suffix_len = ctx->ref_cds_len - prefix_len - (size_t)deleted_len;
    /* Clean single deletion: deleted bases unambiguous, alt == ref-prefix + ref-suffix. */
    if (!delta_norm_span_unambiguous(ctx->ref_cds + prefix_len, (size_t)deleted_len) ||
        !delta_norm_span_equal(ctx->alt_cds, ctx->ref_cds, prefix_len) ||
        !delta_norm_span_equal(ctx->alt_cds + prefix_len,
                               ctx->ref_cds + prefix_len + (size_t)deleted_len,
                               suffix_len)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    /* Classify from the peptide diff window (already trimmed in the context). Both the ref
     * codons removed/merged and any alt junction codon must be known and NON-stop: a stop in
     * either window means the deletion removed the terminator (stop_lost) or the merged
     * junction is a premature stop (inframe_deletion&stop_gained). Leave those contexts
     * unresolved rather than emit a partial SO set. A codon-aligned deletion has
     * an empty alt window (alt_first_changed_codon == 0); a mid-codon deletion carries a
     * single merged junction codon in the alt window. */
    if (ctx->ref_first_changed_codon == 0u ||
        ctx->ref_last_changed_codon < ctx->ref_first_changed_codon ||
        ctx->ref_last_changed_codon > ctx->ref_peptide_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (!delta_peptide_span_known_nonstop(
            ctx->ref_peptide, (size_t)ctx->ref_first_changed_codon - 1u,
            (size_t)(ctx->ref_last_changed_codon - ctx->ref_first_changed_codon) + 1u)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (ctx->alt_first_changed_codon == 0u && ctx->alt_last_changed_codon != 0u) {
        /* Inconsistent window bounds from a malformed direct context: refuse to classify. */
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (ctx->alt_first_changed_codon != 0u) {
        if (ctx->alt_last_changed_codon < ctx->alt_first_changed_codon ||
            ctx->alt_last_changed_codon > ctx->alt_peptide_len ||
            !delta_peptide_span_known_nonstop(
                ctx->alt_peptide, (size_t)ctx->alt_first_changed_codon - 1u,
                (size_t)(ctx->alt_last_changed_codon - ctx->alt_first_changed_codon) + 1u)) {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
    }
    protein_cds = ctx->single_edit_cds_start;
    if ((uint64_t)((protein_cds - 1u) / 3u) + 1u > (uint64_t)INT32_MAX) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    delta->inframe_deletion = 1u;
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (int32_t)(((protein_cds - 1u) / 3u) + 1u);
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

static int delta_inframe_insertion_site(uint32_t before_cds, size_t cds_len);

typedef struct delta_local_peptide_window {
    size_t ref_nt_off;
    size_t ref_nt_len;
    size_t alt_nt_off;
    size_t alt_nt_len;
    size_t ref_off;
    size_t ref_len;
    size_t alt_off;
    size_t alt_len;
} delta_local_peptide_window_t;

/* TranscriptVariationAllele::codon does not compare whole translated proteins. It rounds the
 * reference edit span out to codon boundaries and adjusts that nucleotide window by the net
 * allele-length change. For an insertion on a codon boundary the reference window is empty;
 * inside a codon it contains that one flanking codon. The corresponding peptide slices are
 * already present in the full CodingContext, so deriving the views costs no translation or
 * allocation. */
static duckvep_context_delta_status_t delta_context_local_peptide_window(
    const duckvep_coding_context_t *ctx,
    delta_local_peptide_window_t   *window) {

    uint64_t edit_start0;
    uint64_t edit_end0;
    uint64_t codon_start0;
    uint64_t codon_end0;
    uint64_t remainder;
    uint64_t ref_nt_len;
    uint64_t alt_nt_len;
    uint64_t length_change;

    if (window != NULL) memset(window, 0, sizeof *window);
    if (ctx == NULL || window == NULL || !ctx->has_single_edit ||
        ctx->single_edit_cds_start == 0u || ctx->length_diff == 0 ||
        ctx->length_diff == INT64_MIN || (ctx->length_diff % 3) != 0 ||
        ctx->ref_cds == NULL || ctx->alt_cds == NULL ||
        ctx->ref_peptide == NULL || ctx->alt_peptide == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }

    edit_start0 = (uint64_t)ctx->single_edit_cds_start - 1u;
    if (edit_start0 > (uint64_t)ctx->ref_cds_len ||
        (uint64_t)ctx->single_edit_ref_len >
            (uint64_t)ctx->ref_cds_len - edit_start0) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    edit_end0 = edit_start0 + (uint64_t)ctx->single_edit_ref_len;
    codon_start0 = edit_start0 - (edit_start0 % 3u);
    codon_end0 = edit_end0;
    remainder = codon_end0 % 3u;
    if (remainder != 0u) {
        if (codon_end0 > UINT64_MAX - (3u - remainder)) {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        codon_end0 += 3u - remainder;
    }
    if (codon_end0 < codon_start0 || codon_end0 > (uint64_t)ctx->ref_cds_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    ref_nt_len = codon_end0 - codon_start0;
    if (ctx->length_diff > 0) {
        length_change = (uint64_t)ctx->length_diff;
        if (ref_nt_len > UINT64_MAX - length_change) {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        alt_nt_len = ref_nt_len + length_change;
    } else {
        length_change = (uint64_t)(-(ctx->length_diff + 1)) + 1u;
        if (length_change > ref_nt_len) {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        alt_nt_len = ref_nt_len - length_change;
    }
    if ((ref_nt_len % 3u) != 0u || (alt_nt_len % 3u) != 0u ||
        codon_start0 / 3u > (uint64_t)SIZE_MAX ||
        ref_nt_len / 3u > (uint64_t)SIZE_MAX ||
        alt_nt_len / 3u > (uint64_t)SIZE_MAX) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    if (codon_start0 > (uint64_t)SIZE_MAX || ref_nt_len > (uint64_t)SIZE_MAX ||
        alt_nt_len > (uint64_t)SIZE_MAX) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    window->ref_nt_off = (size_t)codon_start0;
    window->alt_nt_off = window->ref_nt_off;
    window->ref_nt_len = (size_t)ref_nt_len;
    window->alt_nt_len = (size_t)alt_nt_len;
    if (window->ref_nt_off > ctx->ref_cds_len ||
        window->ref_nt_len > ctx->ref_cds_len - window->ref_nt_off ||
        window->alt_nt_off > ctx->alt_cds_len ||
        window->alt_nt_len > ctx->alt_cds_len - window->alt_nt_off) {
        memset(window, 0, sizeof *window);
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    window->ref_off = (size_t)(codon_start0 / 3u);
    window->alt_off = window->ref_off;
    window->ref_len = (size_t)(ref_nt_len / 3u);
    window->alt_len = (size_t)(alt_nt_len / 3u);
    if (window->ref_off > ctx->ref_peptide_len ||
        window->ref_len > ctx->ref_peptide_len - window->ref_off ||
        window->alt_off > ctx->alt_peptide_len ||
        window->alt_len > ctx->alt_peptide_len - window->alt_off) {
        memset(window, 0, sizeof *window);
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    return DUCKVEP_CONTEXT_DELTA_OK;
}

/* VEP's inframe_insertion test: the reference peptide of its codon-local window is a
 * prefix or suffix of the alternate window. An empty reference window therefore matches a
 * clean codon-boundary insertion. The caller must prove both views are in bounds. */
static int delta_pep_window_prefix_or_suffix(
    const uint8_t *outer, size_t outer_off, size_t outer_len,
    const uint8_t *inner, size_t inner_off, size_t inner_len) {

    size_t i;
    int ok;
    if (inner_len == 0u) return 1;
    if (outer == NULL || inner == NULL || inner_len > outer_len) return 0;
    ok = 1;
    for (i = 0u; i < inner_len; i++) {
        if (outer[outer_off + i] != inner[inner_off + i]) { ok = 0; break; }
    }
    if (ok) return 1;
    ok = 1;
    for (i = 0u; i < inner_len; i++) {
        if (outer[outer_off + outer_len - inner_len + i] != inner[inner_off + i]) { ok = 0; break; }
    }
    return ok;
}

static int delta_norm_window_prefix_or_suffix(
    const uint8_t *outer, size_t outer_len,
    const uint8_t *inner, size_t inner_len) {

    if (outer == NULL || inner == NULL || inner_len > outer_len) return 0;
    if (inner_len == 0u) return 1;
    return delta_norm_span_equal(outer, inner, inner_len) ||
           delta_norm_span_equal(
               outer + outer_len - inner_len, inner, inner_len);
}

static duckvep_context_delta_status_t delta_context_lengthening_predicates(
    const duckvep_coding_context_t *ctx,
    duckvep_sequence_delta_t       *delta) {

    uint64_t ref_edit_len;
    uint64_t alt_edit_len;
    uint64_t length_increase;
    size_t edit_start0;
    size_t ref_suffix_len;
    size_t ref_stop_index = 0u;
    size_t alt_stop_index = 0u;
    size_t alt_shape_len;
    delta_local_peptide_window_t local;
    duckvep_context_delta_status_t status;
    int ref_has_stop;
    int alt_has_stop;
    int stop_retained = 0;
    int raw_preserves_ref;
    int ref_peptide_is_alt_suffix;
    int inframe;

    if (ctx == NULL || delta == NULL) return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    if (ctx->ref_cds == NULL || ctx->alt_cds == NULL || ctx->ref_peptide == NULL ||
        ctx->alt_peptide == NULL || ctx->ref_cds_len == 0u ||
        (ctx->ref_cds_len % 3u) != 0u || !ctx->has_single_edit ||
        ctx->applied_edits != 1u || ctx->length_diff <= 0 ||
        ctx->single_edit_alt_len <= ctx->single_edit_ref_len ||
        ctx->cds_changed == 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    ref_edit_len = (uint64_t)ctx->single_edit_ref_len;
    alt_edit_len = (uint64_t)ctx->single_edit_alt_len;
    length_increase = alt_edit_len - ref_edit_len;
    if (length_increase != (uint64_t)ctx->length_diff ||
        (length_increase % 3u) != 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (length_increase > SIZE_MAX ||
        ctx->ref_cds_len > SIZE_MAX - (size_t)length_increase ||
        ctx->alt_cds_len != ctx->ref_cds_len + (size_t)length_increase ||
        ctx->alt_cds_len < ctx->ref_cds_len || ctx->ref_peptide_len != ctx->ref_cds_len / 3u ||
        ctx->alt_peptide_len != ctx->alt_cds_len / 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (ctx->single_edit_cds_start == 0u) return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    edit_start0 = (size_t)ctx->single_edit_cds_start - 1u;
    if (edit_start0 > ctx->ref_cds_len ||
        (size_t)ctx->single_edit_ref_len > ctx->ref_cds_len - edit_start0 ||
        (size_t)ctx->single_edit_alt_len > ctx->alt_cds_len - edit_start0) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    ref_suffix_len = ctx->ref_cds_len - edit_start0 - (size_t)ctx->single_edit_ref_len;
    if (ctx->alt_cds_len - edit_start0 - (size_t)ctx->single_edit_alt_len !=
        ref_suffix_len ||
        !delta_norm_span_equal(ctx->alt_cds, ctx->ref_cds, edit_start0) ||
        !delta_norm_span_unambiguous(ctx->ref_cds + edit_start0,
                                     (size_t)ctx->single_edit_ref_len) ||
        !delta_norm_span_unambiguous(ctx->alt_cds + edit_start0,
                                     (size_t)ctx->single_edit_alt_len) ||
        !delta_norm_span_equal(
            ctx->alt_cds + edit_start0 + (size_t)ctx->single_edit_alt_len,
            ctx->ref_cds + edit_start0 + (size_t)ctx->single_edit_ref_len,
            ref_suffix_len)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    status = delta_context_local_peptide_window(ctx, &local);
    if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
    if (!delta_peptide_span_known(ctx->ref_peptide, local.ref_off, local.ref_len) ||
        !delta_peptide_span_known(ctx->alt_peptide, local.alt_off, local.alt_len) ||
        (uint64_t)local.ref_off + 1u > (uint64_t)INT32_MAX) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    ref_has_stop = delta_peptide_stop_index(
        ctx->ref_peptide, local.ref_off, local.ref_len, &ref_stop_index);
    alt_has_stop = delta_peptide_stop_index(
        ctx->alt_peptide, local.alt_off, local.alt_len, &alt_stop_index);

    /* VEP's ref_eq_alt_sequence has two retained-stop cases: the changed windows keep a
     * stop at the same peptide position, or the complete reference peptide is preserved
     * and the newly appended peptide begins with a stop. The latter is the common
     * insertion immediately before a terminal stop. */
    if (ref_has_stop && alt_has_stop && ref_stop_index == alt_stop_index) {
        stop_retained = 1;
    } else if (ctx->ref_peptide_len > 0u &&
               ctx->alt_peptide_len > ctx->ref_peptide_len &&
               memcmp(ctx->ref_peptide, ctx->alt_peptide,
                      ctx->ref_peptide_len) == 0 &&
               ctx->alt_peptide[ctx->ref_peptide_len] == (uint8_t)'*') {
        stop_retained = 1;
    }

    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (int32_t)(local.ref_off + 1u);
    if (ref_has_stop && !alt_has_stop) {
        delta->stop_lost = 1u;
    }
    if (stop_retained) delta->stop_retained = 1u;
    else if (alt_has_stop && !ref_has_stop) delta->stop_gained = 1u;

    alt_shape_len = alt_has_stop ? alt_stop_index + 1u : local.alt_len;
    /* VEP trims the alternate at its first stop only for the inframe-insertion predicate.
     * protein_altering_variant independently tests the untrimmed peptide. That asymmetry can
     * intentionally leave stop_gained as the sole coding term when the preserved reference
     * flank occurs only after the new stop. */
    inframe = local.ref_len == 0u ||
        delta_pep_window_prefix_or_suffix(
            ctx->alt_peptide, local.alt_off, alt_shape_len,
            ctx->ref_peptide, local.ref_off, local.ref_len);
    if (local.ref_len == 1u && alt_shape_len == 1u &&
        ctx->ref_peptide[local.ref_off] == (uint8_t)'*' &&
        ctx->alt_peptide[local.alt_off] == (uint8_t)'*') {
        inframe = 0;
    }
    raw_preserves_ref = delta_pep_window_prefix_or_suffix(
        ctx->alt_peptide, local.alt_off, local.alt_len,
        ctx->ref_peptide, local.ref_off, local.ref_len);
    ref_peptide_is_alt_suffix =
        local.ref_len <= local.alt_len &&
        (local.ref_len == 0u ||
         memcmp(ctx->alt_peptide + local.alt_off + local.alt_len - local.ref_len,
                ctx->ref_peptide + local.ref_off, local.ref_len) == 0);
    /* VEP suppresses both insertion shape predicates when the start is lost.
     * It also suppresses inframe_insertion when an insertion before the start
     * retains the complete reference peptide as the alternate suffix: that is
     * an upstream addition, not an insertion into the translated protein. */
    if (delta->start_lost ||
        (delta->start_retained && ref_peptide_is_alt_suffix)) {
        inframe = 0;
    } else if (inframe) {
        delta->inframe_insertion = 1u;
    } else if ((local.ref_len == 0u ||
                ctx->ref_peptide[local.ref_off] != (uint8_t)'*') &&
               (local.alt_len == 0u ||
                ctx->alt_peptide[local.alt_off] != (uint8_t)'*') &&
               !raw_preserves_ref && !delta->start_lost) {
        delta->protein_altering = 1u;
    }
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

static duckvep_context_delta_status_t delta_context_shortening_delins(
    const duckvep_coding_context_t *ctx,
    duckvep_sequence_delta_t       *delta) {

    uint64_t ref_len64;
    uint64_t alt_len64;
    uint64_t deleted_len64;
    uint64_t ref_edit_end;
    size_t ref_len;
    size_t alt_len;
    size_t deleted_len;
    size_t prefix_len;
    size_t ref_suffix_len;
    size_t trim_prefix = 0u;
    size_t trim_suffix = 0u;
    size_t ref_remainder;
    size_t alt_remainder;
    delta_local_peptide_window_t local;
    duckvep_context_delta_status_t status;
    int inframe;
    if (ctx == NULL || delta == NULL) return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    if (ctx->ref_cds == NULL || ctx->alt_cds == NULL || ctx->ref_peptide == NULL ||
        ctx->alt_peptide == NULL || ctx->ref_cds_len == 0u ||
        (ctx->ref_cds_len % 3u) != 0u || !ctx->has_single_edit ||
        ctx->applied_edits != 1u || ctx->length_diff >= 0 ||
        ctx->length_diff == INT64_MIN ||
        ctx->single_edit_ref_len == 0u || ctx->single_edit_alt_len == 0u ||
        ctx->cds_changed == 0u || ctx->single_edit_cds_start == 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    ref_len64 = (uint64_t)ctx->single_edit_ref_len;
    alt_len64 = (uint64_t)ctx->single_edit_alt_len;
    if ((ctx->length_diff % 3) != 0 || ref_len64 > SIZE_MAX ||
        alt_len64 > SIZE_MAX) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    deleted_len64 = (uint64_t)(-(ctx->length_diff + 1)) + 1u;
    if (ref_len64 <= alt_len64 || ref_len64 - alt_len64 != deleted_len64) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (deleted_len64 > SIZE_MAX) return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    ref_len = (size_t)ref_len64;
    alt_len = (size_t)alt_len64;
    deleted_len = (size_t)deleted_len64;
    if (ctx->ref_cds_len < deleted_len ||
        ctx->alt_cds_len != ctx->ref_cds_len - deleted_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if ((ctx->alt_cds_len % 3u) != 0u ||
        ctx->ref_peptide_len != ctx->ref_cds_len / 3u ||
        ctx->alt_peptide_len != ctx->alt_cds_len / 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    ref_edit_end = (uint64_t)ctx->single_edit_cds_start + ref_len64 - 1u;
    if (ctx->single_edit_cds_start <= 3u ||
        ref_edit_end > (uint64_t)ctx->ref_cds_len - 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    prefix_len = (size_t)ctx->single_edit_cds_start - 1u;
    if (prefix_len > ctx->ref_cds_len || ref_len > ctx->ref_cds_len - prefix_len ||
        prefix_len > ctx->alt_cds_len || alt_len > ctx->alt_cds_len - prefix_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    ref_suffix_len = ctx->ref_cds_len - prefix_len - ref_len;
    if (ctx->alt_cds_len - prefix_len - alt_len != ref_suffix_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (!delta_norm_span_equal(ctx->alt_cds, ctx->ref_cds, prefix_len) ||
        !delta_norm_span_unambiguous(ctx->ref_cds + prefix_len, ref_len) ||
        !delta_norm_span_unambiguous(ctx->alt_cds + prefix_len, alt_len) ||
        !delta_norm_span_equal(ctx->alt_cds + prefix_len + alt_len,
                               ctx->ref_cds + prefix_len + ref_len,
                               ref_suffix_len)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    status = delta_context_local_peptide_window(ctx, &local);
    if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
    if (!delta_peptide_span_known_nonstop(
            ctx->ref_peptide, local.ref_off, local.ref_len) ||
        !delta_peptide_span_known_nonstop(
            ctx->alt_peptide, local.alt_off, local.alt_len) ||
        (uint64_t)local.ref_off + 1u > (uint64_t)INT32_MAX) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    /* VariationEffect::inframe_deletion compares the actual codon-local strings. A
     * net -3 delins is not automatically deletion-shaped: the shorter alternate must
     * preserve a prefix/suffix, or disappear completely after common-edge trimming
     * while the remaining reference length is divisible by three. */
    inframe = delta_norm_window_prefix_or_suffix(
        ctx->ref_cds + local.ref_nt_off, local.ref_nt_len,
        ctx->alt_cds + local.alt_nt_off, local.alt_nt_len);
    if (!inframe) {
        while (trim_prefix < local.ref_nt_len &&
               trim_prefix < local.alt_nt_len &&
               delta_norm_span_equal(
                   ctx->ref_cds + local.ref_nt_off + trim_prefix,
                   ctx->alt_cds + local.alt_nt_off + trim_prefix, 1u)) {
            trim_prefix++;
        }
        while (trim_suffix < local.ref_nt_len - trim_prefix &&
               trim_suffix < local.alt_nt_len - trim_prefix &&
               delta_norm_span_equal(
                   ctx->ref_cds + local.ref_nt_off + local.ref_nt_len - 1u - trim_suffix,
                   ctx->alt_cds + local.alt_nt_off + local.alt_nt_len - 1u - trim_suffix,
                   1u)) {
            trim_suffix++;
        }
        ref_remainder = local.ref_nt_len - trim_prefix - trim_suffix;
        alt_remainder = local.alt_nt_len - trim_prefix - trim_suffix;
        inframe = alt_remainder == 0u && (ref_remainder % 3u) == 0u;
    }
    if (inframe) delta->inframe_deletion = 1u;
    else delta->protein_altering = 1u;
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (int32_t)(local.ref_off + 1u);
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

/* VEP local-window stop_gained for a frameshift allele.
 *
 * The critical semantic (read from ensembl-variation
 * TranscriptVariationAllele::codon + Utils/VariationEffect::stop_gained): VEP does NOT
 * retranslate the shifted reading frame to the next downstream stop. It recomputes ONLY
 * the codons overlapping the variant — translation_start..translation_end — extended by
 * the net length delta, taken from the EDITED CDS, and calls stop_gained when a COMPLETE
 * codon in that LOCAL window is a stop the reference window did not already have. That is
 * why only a small minority of frameshifts are stop_gained: a full downstream
 * retranslation would make nearly all of them stop_gained. Reusing the kernel's full alt
 * peptide (delta_translate_cds_full, which keeps internal '*' as ordinary bytes) to scan
 * for a "new" stop would therefore over-fire and regress every clean frameshift.
 *
 * codon_cds_start = translation_start*3 - 2 is codon-aligned in BOTH the ref and alt
 * frames, so the window's whole codons are exactly the already-translated peptide bytes
 * starting at peptide index (translation_start - 1). We read those bytes (which carry the
 * correct std/mito codon table applied at build) rather than re-deriving translation.
 *
 * Returns 1 = stop_gained, 0 = not (caller keeps the bare frameshift fact). */
static int delta_frameshift_local_stop_gained(const duckvep_coding_context_t *ctx) {
    uint32_t first_cds;      /* 1-based CDS position of the first affected base */
    uint64_t last_cds;       /* 1-based; for a pure insertion VEP sets cds_end = cds_start-1 */
    uint64_t tv_tr_start, tv_tr_end;
    uint64_t pep_off;        /* 0-based peptide index of the codon-aligned window start */
    uint64_t off_alt;        /* 0-based, codon-aligned CDS offset of the window start */
    uint64_t ref_codons;     /* codon_len/3; 0 for a codon-boundary insertion */
    uint64_t alt_codons;
    uint64_t avail_alt;
    int64_t  alt_window_len;
    uint64_t i;

    if (ctx == NULL || ctx->alt_peptide == NULL || ctx->ref_peptide == NULL) return 0;
    first_cds = ctx->single_edit_cds_start;
    if (first_cds == 0u) return 0;
    tv_tr_start = ((uint64_t)(first_cds - 1u) / 3u) + 1u;
    pep_off = tv_tr_start - 1u;
    off_alt = pep_off * 3u;

    /* Reference codon span of the window (codon_len/3). For a pure insertion VEP sets
     * cds_end = cds_start - 1; an insertion before the first CDS base (first_cds == 1)
     * spans no reference codon, so handle it explicitly and never form last_cds == 0
     * (whose (last_cds - 1) would underflow). This helper is also called on direct
     * frameshift rows without the start-codon guard, so first_cds == 1 can reach here. */
    if (ctx->single_edit_ref_len > 0u) {
        last_cds = (uint64_t)first_cds + (uint64_t)ctx->single_edit_ref_len - 1u;
        tv_tr_end = ((last_cds - 1u) / 3u) + 1u;
        ref_codons = (tv_tr_end >= tv_tr_start) ? (tv_tr_end - tv_tr_start + 1u) : 0u;
    } else if (first_cds > 1u) {
        last_cds = (uint64_t)first_cds - 1u;
        tv_tr_end = ((last_cds - 1u) / 3u) + 1u;
        ref_codons = (tv_tr_end >= tv_tr_start) ? (tv_tr_end - tv_tr_start + 1u) : 0u;
    } else {
        ref_codons = 0u;
    }

    /* Alt window length = codon_len + net length delta, clamped to the available alt CDS,
     * exactly as VEP's substr(alt_cds, codon_cds_start-1, codon_len + net_delta) would. */
    alt_window_len = (int64_t)(ref_codons * 3u) + ctx->length_diff;
    if (alt_window_len <= 0) return 0;            /* no whole alt codon -> pep 'X', no stop */
    if (off_alt > (uint64_t)ctx->alt_cds_len) return 0;
    avail_alt = (uint64_t)ctx->alt_cds_len - off_alt;
    if ((uint64_t)alt_window_len > avail_alt) alt_window_len = (int64_t)avail_alt;
    alt_codons = (uint64_t)alt_window_len / 3u;

    if (pep_off >= (uint64_t)ctx->alt_peptide_len) return 0;
    if (alt_codons > (uint64_t)ctx->alt_peptide_len - pep_off) {
        alt_codons = (uint64_t)ctx->alt_peptide_len - pep_off;
    }
    {
        int alt_has_stop = 0;
        for (i = 0u; i < alt_codons; i++) {
            if (ctx->alt_peptide[pep_off + i] == (uint8_t)'*') { alt_has_stop = 1; break; }
        }
        if (!alt_has_stop) return 0;
    }
    /* Suppress when the reference window already carried a stop over the same codons. */
    if (pep_off < (uint64_t)ctx->ref_peptide_len) {
        uint64_t rc = ref_codons;
        if (rc > (uint64_t)ctx->ref_peptide_len - pep_off) {
            rc = (uint64_t)ctx->ref_peptide_len - pep_off;
        }
        for (i = 0u; i < rc; i++) {
            if (ctx->ref_peptide[pep_off + i] == (uint8_t)'*') return 0;
        }
    }
    return 1;
}

/* Layer VEP's start-codon predicate on a length-changing edit. The alternate
 * CDS has already been rebuilt, so the resulting start codon is read once here
 * instead of being reconstructed independently by every shape classifier.
 * cds_start is the replaced CDS base, or the base before which an insertion is
 * written; cds_start == 4 is therefore after the complete start codon. */
static duckvep_context_delta_status_t delta_context_start_facts(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    int                            *overlaps_start,
    duckvep_sequence_delta_t       *delta) {

    char b0, b1, b2;

    if (overlaps_start != NULL) *overlaps_start = 0;
    if (ctx == NULL || overlaps_start == NULL || delta == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if ((tx_flags & (uint64_t)DUCKVEP_TX_CDS_START_NF) != 0u ||
        !ctx->has_single_edit || ctx->length_diff == 0 ||
        ctx->single_edit_cds_start == 0u || ctx->single_edit_cds_start > 3u) {
        return DUCKVEP_CONTEXT_DELTA_OK;
    }
    *overlaps_start = 1;
    if (ctx->alt_cds == NULL || ctx->alt_cds_len < 3u) {
        delta->start_lost = 1u;
        return DUCKVEP_CONTEXT_DELTA_OK;
    }
    b0 = delta_norm_base((char)ctx->alt_cds[0]);
    b1 = delta_norm_base((char)ctx->alt_cds[1]);
    b2 = delta_norm_base((char)ctx->alt_cds[2]);
    if (b0 == '\0' || b1 == '\0' || b2 == '\0' ||
        b0 == 'N' || b1 == 'N' || b2 == 'N') {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (b0 == 'A' && b1 == 'T' && b2 == 'G') {
        delta->start_retained = 1u;
    } else {
        delta->start_lost = 1u;
    }
    return DUCKVEP_CONTEXT_DELTA_OK;
}

/* Reproduce VEP's stop predicates for a frame-changing edit that overlaps the
 * terminal codon. VEP rebuilds translateable CDS + 3' UTR, then translates the
 * three bases at the ORIGINAL coding end. A stop retained there suppresses the
 * frameshift term. A lost stop suppresses frameshift only when the first affected
 * reference peptide is itself the stop; an edit beginning upstream emits both
 * frameshift_variant and stop_lost.
 *
 * A pure insertion immediately BEFORE the terminal codon has reversed CDS
 * coordinates [terminal_start, terminal_start-1] and does not overlap it. An
 * insertion after its first or second base does. This distinction is why an
 * insertion at the CDS/stop boundary remains an ordinary frameshift. */
static duckvep_context_delta_status_t delta_context_frameshift_stop_facts(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    int                            *suppress_frameshift,
    duckvep_sequence_delta_t       *delta) {

    uint32_t terminal_start;
    uint64_t edit_end;
    char codon[4];
    size_t i;
    char alt_aa;

    if (suppress_frameshift != NULL) *suppress_frameshift = 0;
    if (ctx == NULL || suppress_frameshift == NULL || delta == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (!ctx->has_single_edit || ctx->ref_cds == NULL || ctx->alt_cds == NULL ||
        ctx->ref_peptide == NULL || ctx->ref_cds_len < 3u ||
        (ctx->ref_cds_len % 3u) != 0u ||
        ctx->ref_peptide_len != ctx->ref_cds_len / 3u ||
        ctx->ref_peptide[ctx->ref_peptide_len - 1u] != (uint8_t)'*' ||
        (tx_flags & (uint64_t)DUCKVEP_TX_CDS_END_NF) != 0u ||
        ctx->single_edit_cds_start == 0u ||
        ctx->ref_cds_len > (size_t)UINT32_MAX) {
        return DUCKVEP_CONTEXT_DELTA_OK;
    }

    terminal_start = (uint32_t)ctx->ref_cds_len - 2u;
    if (ctx->single_edit_ref_len == 0u) {
        if (ctx->single_edit_cds_start <= terminal_start ||
            ctx->single_edit_cds_start > ctx->ref_cds_len) {
            return DUCKVEP_CONTEXT_DELTA_OK;
        }
    } else {
        edit_end = (uint64_t)ctx->single_edit_cds_start +
                   (uint64_t)ctx->single_edit_ref_len - 1u;
        if (ctx->single_edit_cds_start > ctx->ref_cds_len ||
            edit_end < (uint64_t)terminal_start) {
            return DUCKVEP_CONTEXT_DELTA_OK;
        }
    }

    codon[3] = '\0';
    for (i = 0u; i < 3u; i++) {
        size_t position = ctx->ref_cds_len - 3u + i;
        char base;

        if (position < ctx->alt_cds_len) {
            base = delta_norm_base((char)ctx->alt_cds[position]);
        } else {
            size_t tail_position = position - ctx->alt_cds_len;
            if (ctx->post_cds_bases == NULL ||
                tail_position >= (size_t)ctx->post_cds_length) {
                return DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL;
            }
            base = delta_norm_base((char)ctx->post_cds_bases[tail_position]);
        }
        if (base == '\0' || base == 'N') {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        codon[i] = base;
    }
    if (ctx->codon_table != (uint8_t)DUCKVEP_CODON_TABLE_STANDARD &&
        ctx->codon_table != (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    alt_aa = duckvep_translate_codon(
        codon, (duckvep_codon_table_t)ctx->codon_table);
    if (alt_aa == '*') {
        delta->stop_retained = 1u;
        *suppress_frameshift = 1;
    } else {
        delta->stop_lost = 1u;
        *suppress_frameshift =
            ctx->single_edit_cds_start >= terminal_start;
    }
    if (*suppress_frameshift) {
        delta->cdna_pos = -1;
        delta->cds_pos = -1;
        delta->protein_pos = (int32_t)ctx->ref_peptide_len;
        delta->ref_aa = (uint8_t)'*';
        delta->alt_aa = (uint8_t)alt_aa;
        delta->valid = 1u;
    }
    return DUCKVEP_CONTEXT_DELTA_OK;
}

typedef struct delta_vep_peptide_summary {
    size_t  complete_length;
    uint8_t first_aa;
    uint8_t has_stop;
    uint8_t has_x;
    uint8_t starts_with_stop;
    uint8_t exactly_stop;
} delta_vep_peptide_summary_t;

/* Summarize VEP's local peptide without materializing it. TVA::peptide translates
 * complete codons and appends X for a trailing partial codon, except when the
 * complete peptide is exactly "*". That exception is observable at the terminal
 * stop and must precede the length-based frameshift test. */
static duckvep_context_delta_status_t delta_vep_peptide_summarize(
    const uint8_t                  *cds,
    size_t                          cds_len,
    uint8_t                         codon_table,
    delta_vep_peptide_summary_t    *summary) {

    size_t codons;
    size_t i;

    if (summary == NULL || (cds == NULL && cds_len != 0u)) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    memset(summary, 0, sizeof *summary);
    if (codon_table != (uint8_t)DUCKVEP_CODON_TABLE_STANDARD &&
        codon_table != (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }

    codons = cds_len / 3u;
    for (i = 0u; i < codons; i++) {
        char triplet[4];
        size_t j;
        char aa;

        for (j = 0u; j < 3u; j++) {
            char base = delta_norm_base((char)cds[i * 3u + j]);
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
            }
            triplet[j] = base;
        }
        triplet[3] = '\0';
        aa = duckvep_translate_codon(
            triplet, (duckvep_codon_table_t)codon_table);
        if (i == 0u) summary->first_aa = (uint8_t)aa;
        if (aa == '*') summary->has_stop = 1u;
    }
    summary->complete_length = codons;
    summary->starts_with_stop =
        codons > 0u && summary->first_aa == (uint8_t)'*';
    summary->exactly_stop =
        codons == 1u && summary->first_aa == (uint8_t)'*';
    if ((cds_len % 3u) != 0u && !summary->exactly_stop) {
        summary->has_x = 1u;
    }
    return DUCKVEP_CONTEXT_DELTA_OK;
}

/* Terminal insertions expose a VEP predicate-order state that a modulo-three
 * classifier cannot represent. At the boundary before the stop, VEP's reference
 * local peptide is empty. Inside the stop, it is "*", which suppresses frameshift;
 * stop, insertion, and coding_unknown predicates then remain independently true.
 * See design/duckvep_errata.md for the pinned source anchors and witnesses. */
static duckvep_context_delta_status_t delta_context_terminal_stop_insertion(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    int                            *handled,
    duckvep_sequence_delta_t       *delta) {

    uint32_t terminal_start;
    size_t inserted_len;
    size_t prefix_len;
    size_t local_off;
    size_t local_len;
    delta_vep_peptide_summary_t local;
    duckvep_context_delta_status_t status;

    if (handled != NULL) *handled = 0;
    if (ctx == NULL || handled == NULL || delta == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (!ctx->has_single_edit || ctx->applied_edits != 1u ||
        ctx->single_edit_ref_len != 0u || ctx->single_edit_alt_len == 0u ||
        ctx->length_diff <= 0 || ctx->ref_cds == NULL || ctx->alt_cds == NULL ||
        ctx->ref_peptide == NULL || ctx->ref_cds_len < 3u ||
        (ctx->ref_cds_len % 3u) != 0u ||
        ctx->ref_peptide_len != ctx->ref_cds_len / 3u ||
        ctx->ref_peptide[ctx->ref_peptide_len - 1u] != (uint8_t)'*' ||
        (tx_flags & (uint64_t)DUCKVEP_TX_CDS_END_NF) != 0u ||
        ctx->ref_cds_len > (size_t)UINT32_MAX) {
        return DUCKVEP_CONTEXT_DELTA_OK;
    }

    terminal_start = (uint32_t)ctx->ref_cds_len - 2u;
    if (ctx->single_edit_cds_start < terminal_start ||
        ctx->single_edit_cds_start > ctx->ref_cds_len) {
        return DUCKVEP_CONTEXT_DELTA_OK;
    }
    *handled = 1;

    inserted_len = (size_t)ctx->single_edit_alt_len;
    if ((uint64_t)ctx->single_edit_alt_len != (uint64_t)ctx->length_diff ||
        inserted_len > SIZE_MAX - ctx->ref_cds_len ||
        ctx->alt_cds_len != ctx->ref_cds_len + inserted_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    prefix_len = (size_t)ctx->single_edit_cds_start - 1u;
    if (prefix_len > ctx->ref_cds_len ||
        !delta_norm_span_equal(ctx->alt_cds, ctx->ref_cds, prefix_len) ||
        !delta_norm_span_unambiguous(ctx->alt_cds + prefix_len, inserted_len) ||
        !delta_norm_span_equal(ctx->alt_cds + prefix_len + inserted_len,
                               ctx->ref_cds + prefix_len,
                               ctx->ref_cds_len - prefix_len)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    local_off = (size_t)terminal_start - 1u;
    local_len = ctx->single_edit_cds_start == terminal_start
                  ? inserted_len
                  : inserted_len + 3u;
    if (local_off > ctx->alt_cds_len ||
        local_len > ctx->alt_cds_len - local_off) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    status = delta_vep_peptide_summarize(
        ctx->alt_cds + local_off, local_len, ctx->codon_table, &local);
    if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;

    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (int32_t)ctx->ref_peptide_len;
    delta->ref_aa = ctx->single_edit_cds_start == terminal_start
                      ? 0u
                      : (uint8_t)'*';
    delta->alt_aa = local.complete_length > 0u ? local.first_aa : 0u;

    if (ctx->single_edit_cds_start == terminal_start) {
        if (!local.has_x && local.starts_with_stop) {
            delta->stop_retained = 1u;
        } else if ((inserted_len % 3u) != 0u) {
            delta->frameshift = 1u;
        }
        if (!delta->stop_retained && local.has_stop) {
            delta->stop_gained = 1u;
        }
        if (!delta->frameshift) delta->inframe_insertion = 1u;
    } else {
        if (local.has_x) {
            delta_vep_peptide_summary_t original_end;
            status = delta_vep_peptide_summarize(
                ctx->alt_cds + ctx->ref_cds_len - 3u, 3u,
                ctx->codon_table, &original_end);
            if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
            if (!original_end.starts_with_stop) delta->stop_lost = 1u;
        } else if (!local.has_stop) {
            delta->stop_lost = 1u;
        } else if (local.starts_with_stop) {
            delta->stop_retained = 1u;
        }

        if (local.has_stop && !local.exactly_stop) {
            delta->inframe_insertion = 1u;
        }
        if (local.has_x && !delta->stop_lost && !delta->stop_retained) {
            delta->coding_unknown = 1u;
        }
    }
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

/* Coarse frameshift FACT from the CodingContext: a net CDS length change not divisible
 * by three shifts the reading frame, which VEP labels frameshift_variant. Start-codon
 * edits layer start_lost or start_retained_variant on that frame fact. When the
 * frameshift's local recomputed codon is a premature stop, VEP composites
 * frameshift_variant&stop_gained; delta_frameshift_local_stop_gained adds that fact. The
 * region layer adds any splice/UTR composite terms. VEP `--gff`-validated. */
static duckvep_context_delta_status_t delta_context_frameshift(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    duckvep_sequence_delta_t       *delta) {

    if (ctx == NULL || delta == NULL) return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    if (ctx->ref_cds == NULL || ctx->alt_cds == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (ctx->length_diff == 0 || ctx->length_diff == INT64_MIN ||
        (ctx->length_diff % 3) == 0 || ctx->cds_changed == 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (ctx->ref_cds_len == 0u || (ctx->ref_cds_len % 3u) != 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    /* This classifier accepts one edit; grouped haplotypes enter as an edit set elsewhere. */
    if (!ctx->has_single_edit) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    /* A pure insertion after the final CDS base lands at the CDS/3'UTR junction, where VEP
     * composites the frameshift with 3_prime_UTR_variant. The coding frameshift fact alone
     * is an incomplete SO set there, so leave it unresolved. */
    if (ctx->single_edit_ref_len == 0u &&
        (size_t)ctx->single_edit_cds_start > ctx->ref_cds_len) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    delta->frameshift = 1u;
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (ctx->single_edit_cds_start >= 1u)
                           ? (int32_t)(((ctx->single_edit_cds_start - 1u) / 3u) + 1u)
                           : -1;
    delta->ref_aa = 0u;
    delta->alt_aa = 0u;
    if (!delta->stop_lost && !delta->stop_retained &&
        delta_frameshift_local_stop_gained(ctx)) {
        delta->stop_gained = 1u;
    }
    {
        int overlaps_start;
        duckvep_context_delta_status_t status =
            delta_context_start_facts(ctx, tx_flags, &overlaps_start, delta);
        if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
    }
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

DUCKVEP_INTERNAL_API duckvep_context_delta_status_t duckvep_coding_context_delta_fill(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    duckvep_sequence_delta_t       *delta) {

    size_t prefix = 0u;
    size_t suffix = 0u;
    size_t change_end;
    size_t codon_start;
    size_t lo_codon;
    size_t hi_codon;
    size_t c;
    size_t win_len;
    size_t ref_stop = 0u;
    size_t alt_stop = 0u;
    int has_ref_stop = 0;
    int has_alt_stop = 0;
    int win_equal = 1;

    if (delta == NULL) return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    memset(delta, 0, sizeof *delta);
    if (ctx == NULL || ctx->ref_cds == NULL || ctx->alt_cds == NULL ||
        ctx->ref_peptide == NULL || ctx->alt_peptide == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (ctx->length_diff > 0 && ctx->has_single_edit &&
        ctx->single_edit_ref_len == 0u) {
        duckvep_context_delta_status_t terminal_status;
        int handled;

        terminal_status = delta_context_terminal_stop_insertion(
            ctx, tx_flags, &handled, delta);
        if (handled || terminal_status != DUCKVEP_CONTEXT_DELTA_OK) {
            return terminal_status;
        }
    }
    if (ctx->length_diff != 0 && (ctx->length_diff % 3) != 0) {
        duckvep_context_delta_status_t stop_status;
        int suppress_frameshift;

        stop_status = delta_context_frameshift_stop_facts(
            ctx, tx_flags, &suppress_frameshift, delta);
        if (suppress_frameshift || stop_status != DUCKVEP_CONTEXT_DELTA_OK) {
            return stop_status;
        }
        return delta_context_frameshift(ctx, tx_flags, delta);
    }
    if (ctx->length_diff > 0) {
        int overlaps_start;
        duckvep_context_delta_status_t status =
            delta_context_start_facts(ctx, tx_flags, &overlaps_start, delta);
        if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
        return delta_context_lengthening_predicates(ctx, delta);
    }
    if (ctx->length_diff < 0 && ctx->has_single_edit &&
        ctx->single_edit_ref_len != 0u && ctx->single_edit_alt_len != 0u) {
        return delta_context_shortening_delins(ctx, delta);
    }
    if (ctx->length_diff < 0) {
        duckvep_context_delta_status_t status =
            delta_context_inframe_deletion(ctx, delta);
        if (status == DUCKVEP_CONTEXT_DELTA_OK && delta->valid) {
            int overlaps_start;
            status = delta_context_start_facts(
                ctx, tx_flags, &overlaps_start, delta);
        }
        return status;
    }
    if (ctx->ref_cds_len == 0u || ctx->ref_cds_len != ctx->alt_cds_len ||
        ctx->cds_changed == 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    while (prefix < ctx->ref_cds_len) {
        char rb = delta_norm_base((char)ctx->ref_cds[prefix]);
        char ab = delta_norm_base((char)ctx->alt_cds[prefix]);
        if (rb == '\0' || ab == '\0') return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        if (rb != ab) break;
        prefix++;
    }
    if (prefix == ctx->ref_cds_len) return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    while (suffix < ctx->ref_cds_len - prefix) {
        char rb = delta_norm_base((char)ctx->ref_cds[ctx->ref_cds_len - 1u - suffix]);
        char ab = delta_norm_base((char)ctx->alt_cds[ctx->alt_cds_len - 1u - suffix]);
        if (rb == '\0' || ab == '\0') return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        if (rb != ab) break;
        suffix++;
    }

    change_end = ctx->ref_cds_len - suffix;
    codon_start = prefix - (prefix % 3u);
    if (codon_start > SIZE_MAX - 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (change_end == 0u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    /* Pure substitution (length-preserving edit, CDS changed). Locate the affected codon
     * window [lo_codon..hi_codon] from the CDS diff, then classify from the ref/alt peptide
     * windows with VEP's whole-window predicates (VariationEffect.pm _get_peptide_alleles
     * + stop_gained/stop_lost/stop_retained/synonymous/missense/start_lost). Shape-agnostic
     * over the number of changed codons: SNV, same-codon MNV, and cross-codon / N-codon MNV
     * all take this one path. Precedence (from VEP's guards): stop_lost > stop_retained >
     * stop_gained > synonymous/missense; start_lost is overlaid at codon 0 and coexists with
     * a co-occurring stop change (e.g. ATG->TAG => start_lost&stop_gained). */
    lo_codon = codon_start / 3u;
    hi_codon = (change_end - 1u) / 3u;
    if (ctx->ref_peptide_len != ctx->alt_peptide_len) {
        /* Not a pure substitution after full-CDS translation; leave to the length-diff
         * dispatch. Reaching here with unequal peptide lengths would be a contract break. */
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (hi_codon >= ctx->ref_peptide_len || hi_codon >= (size_t)INT32_MAX ||
        (hi_codon + 1u) * 3u > ctx->ref_cds_len ||
        (hi_codon + 1u) * 3u > ctx->alt_cds_len) {
        /* hi_codon >= INT32_MAX keeps the 1-based protein_pos (lo_codon+1) in range; the CDS
         * bound guards the raw-byte codon scan below against a malformed/inconsistent direct
         * context whose peptide length overreaches its CDS. */
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    win_len = hi_codon - lo_codon + 1u;
    /* An ambiguous/unknown residue in the window -> VEP coding_unknown fallback
     * (coding_sequence_variant), matching the incomplete-terminal-codon frontier. Check both
     * the translated peptide (X) and the raw CDS codon bytes (N / non-ACGT): a caller that
     * hands us an inconsistent non-X peptide over an N-containing codon must still not be
     * classified, so we do not trust the peptide alone. */
    for (c = lo_codon; c <= hi_codon; c++) {
        char r = (char)ctx->ref_peptide[c];
        char a = (char)ctx->alt_peptide[c];
        size_t b;
        if (r == '\0' || a == '\0' || r == 'X' || a == 'X') {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        for (b = 0u; b < 3u; b++) {
            char rb = delta_norm_base((char)ctx->ref_cds[c * 3u + b]);
            char ab = delta_norm_base((char)ctx->alt_cds[c * 3u + b]);
            if (rb == '\0' || rb == 'N' || ab == '\0' || ab == 'N') {
                return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
            }
        }
        if (r != a) win_equal = 0;
        if (r == '*' && !has_ref_stop) { ref_stop = c - lo_codon; has_ref_stop = 1; }
        if (a == '*' && !has_alt_stop) { alt_stop = c - lo_codon; has_alt_stop = 1; }
    }

    /* start_lost: the window covers the start codon (codon 1) and the alt no longer encodes
     * Met. VEP treats codon 1 as the start by coordinate (translation_start == 1), so the
     * test is on the alt residue alone — not gated on the reference being Met — matching
     * VariationEffect::start_lost for realistic ATG starts and the genetic-code oracle.
     * CDS_START_NF (annotated-incomplete start) suppresses it, as in VEP. */
    if (lo_codon == 0u && (char)ctx->alt_peptide[0] != 'M' &&
        (tx_flags & (uint64_t)DUCKVEP_TX_CDS_START_NF) == 0u) {
        delta->start_lost = 1u;
    }

    if (has_ref_stop && !has_alt_stop) {
        delta->stop_lost = 1u;
    } else if (has_ref_stop && has_alt_stop && ref_stop == alt_stop) {
        delta->stop_retained = 1u;
    } else if (has_alt_stop && !has_ref_stop) {
        delta->stop_gained = 1u;
    }
    if (!delta->stop_lost && !delta->stop_retained && !delta->stop_gained) {
        if (win_equal) {
            delta->synonymous = 1u;
        } else if (!delta->start_lost) {
            delta->missense = 1u;
        }
    }

    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    /* Single changed codon carries a precise protein position + AA pair; a multi-codon
     * window is reported coarsely (protein_pos = -1, no AA pair) — the SO fact set is what
     * conformance compares, and the range/HGVS detail is rendered later from the context. */
    if (win_len == 1u) {
        delta->protein_pos = (int32_t)(lo_codon + 1u);
        delta->ref_aa = ctx->ref_peptide[lo_codon];
        delta->alt_aa = ctx->alt_peptide[lo_codon];
    } else {
        delta->protein_pos = -1;
    }
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

static int delta_frameshift_body_cds_pos(uint32_t cds_pos, size_t cds_len) {
    return cds_len >= 6u && cds_pos > 3u && (size_t)cds_pos <= cds_len - 3u;
}

static int delta_inframe_insertion_site(uint32_t before_cds, size_t cds_len) {
    return cds_len >= 12u && (cds_len % 3u) == 0u && before_cds > 3u &&
           (size_t)before_cds < cds_len - 3u && (before_cds % 3u) == 0u;
}

static int delta_protein_altering_insertion_site(uint32_t before_cds, size_t cds_len) {
    return cds_len >= 12u && (cds_len % 3u) == 0u && before_cds > 3u &&
           (size_t)before_cds < cds_len - 3u && (before_cds % 3u) != 0u;
}

static char delta_insertion_payload_base(
    const duckvep_variant_batch_t *v,
    const delta_edit_view_t       *edit,
    uint16_t                       payload_i,
    int8_t                         strand) {

    size_t payload_idx;
    payload_idx = strand > 0 ? (size_t)payload_i
                             : (size_t)edit->alt_len - 1u - (size_t)payload_i;
    return delta_orient_genomic_base(
        (char)v->allele_bytes[edit->alt_off + payload_idx], strand);
}

static int delta_insertion_payload_nonstop(
    const duckvep_variant_batch_t *v,
    const delta_edit_view_t       *edit,
    int8_t                         strand,
    duckvep_codon_table_t          table) {

    uint16_t payload_len;
    uint16_t i;
    char codon[4];

    payload_len = edit->alt_len;
    codon[3] = '\0';
    if (payload_len == 0u || (payload_len % 3u) != 0u) return 0;
    for (i = 0u; i < payload_len; i++) {
        char b = delta_insertion_payload_base(v, edit, i, strand);
        if (b == '\0') return 0;
        codon[i % 3u] = b;
        if ((i % 3u) == 2u && duckvep_translate_codon(codon, table) == '*') return 0;
    }
    return 1;
}

static int delta_protein_altering_insertion_nonstop(
    const uint8_t                 *cds_seq,
    uint32_t                       before_cds,
    const duckvep_variant_batch_t *v,
    const delta_edit_view_t       *edit,
    int8_t                         strand,
    duckvep_codon_table_t          table) {

    uint16_t payload_len;
    uint32_t codon_start;
    uint32_t prefix_len;
    uint32_t alt_len;
    uint32_t i;
    char codon[4];

    payload_len = edit->alt_len;
    codon_start = ((before_cds - 1u) / 3u) * 3u + 1u;
    prefix_len = before_cds - codon_start + 1u;
    alt_len = 3u + (uint32_t)payload_len;
    codon[3] = '\0';
    if (prefix_len == 0u || prefix_len >= 3u || (payload_len % 3u) != 0u) return 0;
    for (i = 0u; i < alt_len; i++) {
        char b;
        if (i < prefix_len) {
            b = delta_norm_base((char)cds_seq[(size_t)codon_start - 1u + (size_t)i]);
            if (b == 'N') b = '\0';
        } else if (i < prefix_len + (uint32_t)payload_len) {
            b = delta_insertion_payload_base(
                v, edit, (uint16_t)(i - prefix_len), strand);
        } else {
            uint32_t suffix_i = i - prefix_len - (uint32_t)payload_len;
            b = delta_norm_base(
                (char)cds_seq[(size_t)codon_start - 1u + (size_t)prefix_len + (size_t)suffix_i]);
            if (b == 'N') b = '\0';
        }
        if (b == '\0') return 0;
        codon[i % 3u] = b;
        if ((i % 3u) == 2u && duckvep_translate_codon(codon, table) == '*') return 0;
    }
    return 1;
}

static int delta_body_complete_codon(uint32_t codon_start, size_t cds_len) {
    if (cds_len < 12u || (cds_len % 3u) != 0u) return 0;
    if (codon_start == 0u || (uint64_t)codon_start + 2u > (uint64_t)cds_len) return 0;
    return codon_start > 1u && (uint64_t)codon_start + 2u <= (uint64_t)cds_len - 3u;
}

static int delta_copy_codon_from_cds(
    const uint8_t *cds_seq,
    size_t         cds_len,
    uint32_t       codon_start,
    char           codon[4]) {

    uint32_t j;
    if (codon_start == 0u || (uint64_t)codon_start + 2u > (uint64_t)cds_len) return 0;
    for (j = 0u; j < 3u; j++) {
        char b = delta_norm_base((char)cds_seq[(size_t)codon_start - 1u + (size_t)j]);
        if (b == '\0' || b == 'N') return 0;
        codon[j] = b;
    }
    codon[3] = '\0';
    return 1;
}

static void delta_apply_codon_change(
    const duckvep_transcript_model_t *transcripts,
    size_t                            tx_idx,
    const duckvep_coding_projection_t *proj,
    uint32_t                          change,
    char                              aa_ref,
    char                              aa_alt,
    duckvep_sequence_delta_t         *delta) {

    if (proj->codon_start_cds == 1u && proj->protein_pos == 1u && aa_alt != 'M' &&
        (transcripts->flags[tx_idx] & (uint64_t)DUCKVEP_TX_CDS_START_NF) == 0u) {
        delta->start_lost = 1u;
    }

    if (change & DUCKVEP_CODON_STOP_GAINED)     delta->stop_gained = 1u;
    else if (change & DUCKVEP_CODON_STOP_LOST)  delta->stop_lost = 1u;
    else if (change & DUCKVEP_CODON_MISSENSE) {
        if (!delta->start_lost) delta->missense = 1u;
    }
    else if (change & DUCKVEP_CODON_SYNONYMOUS) {
        if (aa_ref == '*') delta->stop_retained = 1u;
        else               delta->synonymous = 1u;
    } else return;

    delta->cdna_pos = (int32_t)proj->cdna_pos;
    delta->cds_pos = (int32_t)proj->cds_pos;
    delta->protein_pos = (int32_t)proj->protein_pos;
    delta->ref_aa = (uint8_t)aa_ref;
    delta->alt_aa = (uint8_t)aa_alt;
    delta->valid = 1u;
}

/* SNV fast path: ref_len==alt_len==1, the count==1 case of the edit set. Kept
 * private; callers enter through duckvep_sequence_delta_fill. */
static void sequence_delta_fill_snv(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta) {

    duckvep_event_t loaded_event;
    const duckvep_event_t *event;
    duckvep_coding_projection_t proj;
    duckvep_coding_snv_result_t res;
    const uint8_t *cds_seq;
    size_t cds_len;
    char gref, galt;
    duckvep_codon_table_t table;
    duckvep_coding_snv_status_t status;

    delta->valid = 0u;
    delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT;
    if (!delta_transcript_has_sequence(seq, tx_idx)) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_MISSING;
        return;
    }
    if (v->variant_kind == NULL ||
        v->variant_kind[variant_idx] != (uint8_t)DUCKVEP_KIND_SNV) return;
    if (!delta_allele_slice_ok(v, variant_idx)) return;
    if (prepared_event == NULL) {
        duckvep_event_load(v, (size_t)variant_idx, &loaded_event);
        event = &loaded_event;
    } else {
        event = prepared_event;
    }
    if (event->kind != (uint8_t)DUCKVEP_KIND_SNV ||
        event->ref_diff_length != 1u || event->alt_diff_length != 1u) return;
    if (seq->cds_length[tx_idx] == 0u) return;

    pos = event->start1;
    if (!duckvep_project_coding_base(transcripts, exons, tx_idx, pos, &proj)) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
        return;
    }

    cds_seq = seq->cds_bytes + seq->cds_offset[tx_idx];
    cds_len = (size_t)seq->cds_length[tx_idx];
    gref = (char)v->allele_bytes[v->ref_offset[variant_idx] +
                                 event->ref_diff_offset];
    galt = (char)v->allele_bytes[v->alt_offset[variant_idx] +
                                 event->alt_diff_offset];
    table = (duckvep_codon_table_t)seq->codon_table[tx_idx];

    status = duckvep_coding_snv_from_cds(cds_seq, cds_len, &proj, gref, galt,
                                         strand, table, &res);
    delta->sequence_status = delta_sequence_status_from_snv(status);
    if (status != DUCKVEP_CODING_SNV_OK) {
        return;
    }

    delta_apply_codon_change(transcripts, tx_idx, &proj, res.change, res.aa_ref,
                             res.aa_alt, delta);
    if (delta->valid) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    }
}

/* Same-codon MNV path: multiple genomic bases are edited together, but only when the
 * whole REF/ALT span projects to the same coding codon. MNVs spanning codons, introns,
 * UTRs, or incomplete terminal codons deliberately fall back to coding_sequence_variant
 * until the general edit-set/CodingContext path lands. */
static void sequence_delta_fill_mnv(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta) {

    duckvep_coding_projection_t first;
    delta_edit_view_t edit;
    const uint8_t *cds_seq;
    size_t cds_len;
    size_t codon_idx;
    uint16_t len;
    uint16_t j;
    uint8_t seen = 0u;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    duckvep_codon_table_t table;

    (void)pos;
    delta->valid = 0u;
    if (seq == NULL) return;
    if (v->variant_kind == NULL || v->variant_kind[variant_idx] != (uint8_t)DUCKVEP_KIND_MNV) return;
    if (!delta_edit_view_load(v, variant_idx, prepared_event, &edit)) return;
    if (seq->cds_length[tx_idx] == 0u) return;
    if (edit.ref_len == 0u || edit.ref_len != edit.alt_len) return;
    len = edit.ref_len;
    if (len > 3u) return; /* cannot fit inside one codon */

    if (!duckvep_project_coding_base(transcripts, exons, tx_idx,
                                     edit.event.start1, &first)) return;
    cds_seq = seq->cds_bytes + seq->cds_offset[tx_idx];
    cds_len = (size_t)seq->cds_length[tx_idx];
    codon_idx = (size_t)first.codon_start_cds - 1u;
    if (codon_idx > cds_len || cds_len - codon_idx < 3u) return;
    for (j = 0u; j < 3u; j++) {
        char b = delta_norm_base((char)cds_seq[codon_idx + (size_t)j]);
        if (b == '\0' || b == 'N') return;
        ref_codon[j] = b;
        alt_codon[j] = b;
    }
    ref_codon[3] = '\0';
    alt_codon[3] = '\0';

    for (j = 0u; j < len; j++) {
        duckvep_coding_projection_t p;
        char ref_tx;
        char alt_tx;
        uint8_t bit;
        uint32_t gpos = edit.event.start1 + (uint32_t)j;
        if (!duckvep_project_coding_base(transcripts, exons, tx_idx, gpos, &p)) return;
        if (p.codon_start_cds != first.codon_start_cds || p.protein_pos != first.protein_pos) return;
        bit = (uint8_t)(1u << p.codon_offset);
        if (seen & bit) return;
        seen |= bit;
        ref_tx = delta_orient_genomic_base(
            (char)v->allele_bytes[edit.ref_off + (size_t)j], strand);
        alt_tx = delta_orient_genomic_base(
            (char)v->allele_bytes[edit.alt_off + (size_t)j], strand);
        if (ref_tx == '\0' || alt_tx == '\0') return;
        if (ref_codon[p.codon_offset] != ref_tx) return;
        alt_codon[p.codon_offset] = alt_tx;
    }

    table = (duckvep_codon_table_t)seq->codon_table[tx_idx];
    cr = duckvep_codon_change(ref_codon, alt_codon, table);
    if (cr.change & DUCKVEP_CODON_INVALID) return;
    delta_apply_codon_change(transcripts, tx_idx, &first, cr.change, cr.aa_ref,
                             cr.aa_alt, delta);
    /* cDNA/CDS positions are ranges for multi-base alleles in VEP output. The current
     * public result row has scalar fields only, so same-codon MNVs intentionally leave
     * nucleotide coordinate scalars absent until the ABI grows range fields. The
     * protein position is scalar because this path only accepts one-codon edits. */
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
}

/* Narrow cross-codon MNV path: a same-length substitution that touches exactly two
 * adjacent non-terminal CDS body codons can only emit the already-existing missense fact.
 * Stops, first/last codons, intron/CDS-boundary spans, and all-synonymous two-codon edits
 * deliberately fall back until the general edit-set/CodingContext path lands. */
static void sequence_delta_fill_mnv_cross_codon(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta) {

    delta_edit_view_t edit;
    const uint8_t *cds_seq;
    size_t cds_len;
    uint16_t len;
    uint16_t j;
    uint32_t min_cds = UINT32_MAX;
    uint32_t max_cds = 0u;
    uint32_t codon_start[2] = {0u, 0u};
    uint8_t seen[2] = {0u, 0u};
    char ref_codon[2][4];
    char alt_codon[2][4];
    uint32_t codon_count = 0u;
    int aa_changed = 0;
    duckvep_codon_table_t table;

    (void)pos;
    delta->valid = 0u;
    if (seq == NULL) return;
    if (v->variant_kind == NULL || v->variant_kind[variant_idx] != (uint8_t)DUCKVEP_KIND_MNV) return;
    if (!delta_edit_view_load(v, variant_idx, prepared_event, &edit)) return;
    if (seq->cds_length[tx_idx] == 0u) return;
    if (edit.ref_len == 0u || edit.ref_len != edit.alt_len) return;
    len = edit.ref_len;
    if (len < 2u || len > 3u) return;

    cds_seq = seq->cds_bytes + seq->cds_offset[tx_idx];
    cds_len = (size_t)seq->cds_length[tx_idx];
    table = (duckvep_codon_table_t)seq->codon_table[tx_idx];

    for (j = 0u; j < len; j++) {
        duckvep_coding_projection_t p;
        char ref_tx;
        char alt_tx;
        uint32_t k;
        uint32_t slot = UINT32_MAX;
        uint8_t bit;
        uint32_t gpos = edit.event.start1 + (uint32_t)j;

        if (!duckvep_project_coding_base(transcripts, exons, tx_idx, gpos, &p)) return;
        if (!delta_body_complete_codon(p.codon_start_cds, cds_len)) return;
        for (k = 0u; k < codon_count; k++) {
            if (codon_start[k] == p.codon_start_cds) {
                slot = k;
                break;
            }
        }
        if (slot == UINT32_MAX) {
            if (codon_count >= 2u) return;
            slot = codon_count++;
            codon_start[slot] = p.codon_start_cds;
            if (!delta_copy_codon_from_cds(cds_seq, cds_len, p.codon_start_cds,
                                           ref_codon[slot])) return;
            memcpy(alt_codon[slot], ref_codon[slot], sizeof alt_codon[slot]);
        }

        bit = (uint8_t)(1u << p.codon_offset);
        if (seen[slot] & bit) return;
        seen[slot] |= bit;
        ref_tx = delta_orient_genomic_base(
            (char)v->allele_bytes[edit.ref_off + (size_t)j], strand);
        alt_tx = delta_orient_genomic_base(
            (char)v->allele_bytes[edit.alt_off + (size_t)j], strand);
        if (ref_tx == '\0' || alt_tx == '\0') return;
        if (ref_codon[slot][p.codon_offset] != ref_tx) return;
        alt_codon[slot][p.codon_offset] = alt_tx;
        if (p.cds_pos < min_cds) min_cds = p.cds_pos;
        if (p.cds_pos > max_cds) max_cds = p.cds_pos;
    }

    if (codon_count != 2u) return;
    if (max_cds < min_cds || max_cds - min_cds + 1u != (uint32_t)len) return;
    if (codon_start[0] > codon_start[1]) {
        uint32_t tmp = codon_start[0];
        codon_start[0] = codon_start[1];
        codon_start[1] = tmp;
    }
    if (codon_start[1] - codon_start[0] != 3u) return;

    for (j = 0u; j < 2u; j++) {
        duckvep_codon_result_t cr = duckvep_codon_change(
            ref_codon[j], alt_codon[j], table);
        if (cr.change & DUCKVEP_CODON_INVALID) return;
        if (cr.aa_ref == '*' || cr.aa_alt == '*') return;
        if (cr.aa_ref != cr.aa_alt) aa_changed = 1;
    }
    if (!aa_changed) return;

    delta->missense = 1u;
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = -1;
    delta->valid = 1u;
}

/* Narrow normalized-differing-region indel path: non-terminal CDS-body frameshift
 * for INS/DEL/INDEL, codon-boundary in-frame insertion, non-boundary protein-altering
 * insertion, plus codon-aligned in-frame deletion. Raw VCF padding is stripped by
 * delta_edit_view_t, but this still intentionally avoids start/stop-boundary indels,
 * in-frame delins, cross-boundary edit reconstruction, and post-edit peptide comparison
 * until the general edit-set/CodingContext path lands. */
static void sequence_delta_fill_indel_frameshift(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta) {

    delta_edit_view_t edit;
    const uint8_t *cds_seq;
    size_t cds_len;
    int64_t diff;
    int is_frameshift;
    int is_inframe_deletion;
    int is_inframe_insertion;
    uint32_t protein_cds = 0u;
    duckvep_coding_projection_t proj;

    (void)pos;
    delta->valid = 0u;
    if (seq == NULL || seq->cds_length[tx_idx] == 0u) return;
    if (v->variant_kind == NULL) return;
    if (!delta_edit_view_load(v, variant_idx, prepared_event, &edit)) return;
    cds_seq = seq->cds_bytes + seq->cds_offset[tx_idx];
    cds_len = (size_t)seq->cds_length[tx_idx];
    if (edit.ref_len == edit.alt_len) return;
    diff = (int64_t)edit.alt_len - (int64_t)edit.ref_len;
    is_frameshift = (diff % 3) != 0;
    is_inframe_deletion = (diff % 3) == 0 &&
                           v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_DEL;
    is_inframe_insertion = (diff % 3) == 0 &&
                            v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_INS;
    if (!is_frameshift && !is_inframe_deletion && !is_inframe_insertion) return;

    if (v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_INS) {
        uint32_t before_cds;
        uint32_t insertion_cds;
        if (!edit.event.interbase || edit.ref_len != 0u || edit.alt_len == 0u) return;
        if (!duckvep_project_coding_base(transcripts, exons, tx_idx,
                                         edit.event.start1, &proj)) return;
        if (!delta_cds_ref_matches(cds_seq, cds_len, proj.cds_pos,
                                   (char)v->allele_bytes[edit.anchor_ref_off],
                                   strand)) return;
        if (is_frameshift) {
            if (!delta_frameshift_body_cds_pos(proj.cds_pos, cds_len)) return;
            protein_cds = proj.cds_pos;
        } else if (is_inframe_insertion) {
            if ((edit.alt_len % 3u) != 0u) return;
            if (edit.event.anchor_side == (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT) {
                insertion_cds = strand > 0 ? proj.cds_pos + 1u : proj.cds_pos;
            } else {
                insertion_cds = strand > 0 ? proj.cds_pos : proj.cds_pos + 1u;
            }
            if (insertion_cds == 0u) return;
            before_cds = insertion_cds - 1u;
            if (delta_inframe_insertion_site(before_cds, cds_len)) {
                if (!delta_insertion_payload_nonstop(
                        v, &edit, strand,
                        (duckvep_codon_table_t)seq->codon_table[tx_idx])) return;
                delta->inframe_insertion = 1u;
            } else if (delta_protein_altering_insertion_site(before_cds, cds_len)) {
                if (!delta_protein_altering_insertion_nonstop(
                        cds_seq, before_cds, v, &edit, strand,
                        (duckvep_codon_table_t)seq->codon_table[tx_idx])) return;
                delta->protein_altering = 1u;
            } else return;
            protein_cds = before_cds + 1u;
        } else return;
    } else if (v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_DEL ||
               v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_INDEL) {
        uint16_t j;
        uint32_t min_cds = UINT32_MAX;
        uint32_t max_cds = 0u;
        uint32_t deleted_len;
        if (edit.ref_len == 0u) return;
        if (v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_DEL) {
            if (edit.alt_len != 0u) return;
        } else {
            if (!is_frameshift || edit.alt_len == 0u) return;
        }
        for (j = 0u; j < edit.ref_len; j++) {
            uint32_t gpos = edit.event.start1 + (uint32_t)j;
            if (!duckvep_project_coding_base(transcripts, exons, tx_idx, gpos, &proj)) return;
            if (!delta_cds_ref_matches(cds_seq, cds_len, proj.cds_pos,
                                       (char)v->allele_bytes[edit.ref_off + (size_t)j],
                                       strand)) return;
            if (!delta_frameshift_body_cds_pos(proj.cds_pos, cds_len)) return;
            if (proj.cds_pos < min_cds) min_cds = proj.cds_pos;
            if (proj.cds_pos > max_cds) max_cds = proj.cds_pos;
        }
        if (min_cds == UINT32_MAX) return;
        deleted_len = (uint32_t)edit.ref_len;
        if (max_cds < min_cds || max_cds - min_cds + 1u != deleted_len) return;
        if (v->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_INDEL) {
            for (j = 0u; j < edit.alt_len; j++) {
                if (delta_insertion_payload_base(v, &edit, j, strand) == '\0') return;
            }
        } else if (is_inframe_deletion) {
            if (deleted_len % 3u != 0u) return;
            if (((min_cds - 1u) % 3u) != 0u) return;
        }
        protein_cds = min_cds;
    } else return;

    if (is_frameshift) delta->frameshift = 1u;
    else if (is_inframe_deletion) delta->inframe_deletion = 1u;
    else if (!delta->inframe_insertion && !delta->protein_altering) return;
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (int32_t)(((protein_cds - 1u) / 3u) + 1u);
    delta->valid = 1u;
}

static void duckvep_sequence_delta_fill_with_scratch_event(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta) {

    if (delta == NULL) return;
    memset(delta, 0, sizeof *delta);
    delta->sequence_status = delta_transcript_has_sequence(seq, tx_idx)
        ? (uint8_t)DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT
        : (uint8_t)DUCKVEP_SEQUENCE_MISSING;
    switch (kind) {
    case DUCKVEP_KIND_SNV:
        sequence_delta_fill_snv(transcripts, exons, seq, v, variant_idx, tx_idx, pos,
                                strand, prepared_event, delta);
        break;
    case DUCKVEP_KIND_MNV:
        if (scratch != NULL) {
            duckvep_coding_context_t ctx;
            duckvep_variant_coding_context_status_t context_status;
            duckvep_context_delta_status_t delta_status;
            uint64_t tx_flags = 0u;
            if (v == NULL || v->variant_kind == NULL || v->ref_length == NULL ||
                v->alt_length == NULL || variant_idx >= v->count ||
                v->variant_kind[variant_idx] != (uint8_t)DUCKVEP_KIND_MNV ||
                v->ref_length[variant_idx] != v->alt_length[variant_idx]) {
                break;
            }
            if (transcripts != NULL && transcripts->flags != NULL &&
                tx_idx < transcripts->transcript_count) {
                tx_flags = transcripts->flags[tx_idx];
            }
            context_status = duckvep_variant_coding_context_build_event(
                transcripts, exons, seq, v, variant_idx, tx_idx, strand,
                prepared_event, scratch->edits, scratch->edits_cap,
                scratch->alt_cds, scratch->alt_cds_cap,
                scratch->ref_peptide, scratch->ref_peptide_cap,
                scratch->alt_peptide, scratch->alt_peptide_cap, &ctx);
            delta->sequence_status = delta_sequence_status_from_context(context_status);
            if (context_status == DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
                delta_status = duckvep_coding_context_delta_fill(&ctx, tx_flags, delta);
                if (delta_status != DUCKVEP_CONTEXT_DELTA_OK) {
                    delta->sequence_status =
                        delta_sequence_status_from_delta(delta_status);
                }
            }
        } else {
            sequence_delta_fill_mnv(transcripts, exons, seq, v, variant_idx, tx_idx,
                                    pos, strand, prepared_event, delta);
            if (!delta->valid) {
                sequence_delta_fill_mnv_cross_codon(transcripts, exons, seq, v,
                                                    variant_idx, tx_idx, pos, strand,
                                                    prepared_event, delta);
            }
        }
        break;
    case DUCKVEP_KIND_INS:
    case DUCKVEP_KIND_DEL:
    case DUCKVEP_KIND_INDEL:
        if (scratch != NULL) {
            duckvep_coding_context_t ctx;
            duckvep_variant_coding_context_status_t context_status;
            duckvep_context_delta_status_t delta_status;
            uint64_t tx_flags = 0u;
            if (transcripts != NULL && transcripts->flags != NULL &&
                tx_idx < transcripts->transcript_count) {
                tx_flags = transcripts->flags[tx_idx];
            }
            context_status = duckvep_variant_coding_context_build_event(
                transcripts, exons, seq, v, variant_idx, tx_idx, strand,
                prepared_event, scratch->edits, scratch->edits_cap,
                scratch->alt_cds, scratch->alt_cds_cap,
                scratch->ref_peptide, scratch->ref_peptide_cap,
                scratch->alt_peptide, scratch->alt_peptide_cap, &ctx);
            delta->sequence_status = delta_sequence_status_from_context(context_status);
            if (context_status == DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
                delta_status = duckvep_coding_context_delta_fill(&ctx, tx_flags, delta);
                if (delta_status != DUCKVEP_CONTEXT_DELTA_OK) {
                    delta->sequence_status =
                        delta_sequence_status_from_delta(delta_status);
                }
            }
        } else {
            sequence_delta_fill_indel_frameshift(
                transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
                prepared_event, delta);
        }
        break;
    case DUCKVEP_KIND_SV:
    default:
        /* Structural coding edits are not represented by this dispatcher. Leaving
         * the delta invalid keeps the CDS bucket at coding_sequence_variant. */
        break;
    }
    if (delta->valid) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    }
}

DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_with_scratch(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    duckvep_sequence_delta_t         *delta) {

    duckvep_sequence_delta_fill_with_scratch_event(
        kind, transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
        scratch, NULL, delta);
}

DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_sequence_delta_t         *delta) {

    duckvep_sequence_delta_fill_with_scratch(kind, transcripts, exons, seq, v,
                                             variant_idx, tx_idx, pos, strand, NULL,
                                             delta);
}

DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_for_annotation_trace(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_route_t   *route,
    duckvep_sequence_delta_t         *delta) {

    uint8_t context_sequence_status;

    if (route != NULL) *route = DUCKVEP_DELTA_ROUTE_DIRECT;
    if (delta == NULL) return;

    /* SNV (and any non-coding-editable kind, or a call with no edit-set scratch) takes the
     * light projection path: one codon edit, no full-CDS translation. */
    if (scratch == NULL || (kind != DUCKVEP_KIND_MNV && kind != DUCKVEP_KIND_DEL &&
                            kind != DUCKVEP_KIND_INS && kind != DUCKVEP_KIND_INDEL)) {
        duckvep_sequence_delta_fill_with_scratch_event(
            kind, transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
            NULL, prepared_event, delta);
        return;
    }

    /* duckvep_sequence_delta_fill_with_scratch builds the alternate
     * CDS + ref/alt peptides once from the variant's edit set and classifies MNV / DEL / INS
     * / INDEL from the peptide diff (duckvep_coding_context_delta_fill), including the
     * frameshift local-stop composite. Correctness is graded against VEP `--gff`, not
     * against the direct per-shape fallback. */
    duckvep_sequence_delta_fill_with_scratch_event(
        kind, transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
        scratch, prepared_event, delta);
    context_sequence_status = delta->sequence_status;
    if (delta->valid) {
        if (route != NULL) {
            *route = kind == DUCKVEP_KIND_MNV ? DUCKVEP_DELTA_ROUTE_MNV_CONTEXT
                   : kind == DUCKVEP_KIND_DEL ? DUCKVEP_DELTA_ROUTE_DEL_CONTEXT
                   : kind == DUCKVEP_KIND_INS ? DUCKVEP_DELTA_ROUTE_INS_CONTEXT
                                              : DUCKVEP_DELTA_ROUTE_INDEL_CONTEXT;
        }
        return;
    }
    if (context_sequence_status ==
        (uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL) {
        return;
    }

    /* The context classifier deferred this shape or its edit-set build ran out of
     * scratch. Use the direct shape-specific path while that classifier is incomplete.
     * An invalid result renders as coding_sequence_variant (VEP's coding_unknown). */
    duckvep_sequence_delta_fill_with_scratch_event(
        kind, transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
        NULL, prepared_event, delta);
    if (!delta->valid &&
        delta->sequence_status == (uint8_t)DUCKVEP_SEQUENCE_UNSUPPORTED_EDIT &&
        context_sequence_status != (uint8_t)DUCKVEP_SEQUENCE_RESOLVED) {
        delta->sequence_status = context_sequence_status;
    }
    if (route != NULL) {
        *route = kind == DUCKVEP_KIND_MNV ? DUCKVEP_DELTA_ROUTE_MNV_DIRECT_FALLBACK
               : kind == DUCKVEP_KIND_DEL ? DUCKVEP_DELTA_ROUTE_DEL_DIRECT_FALLBACK
               : kind == DUCKVEP_KIND_INS ? DUCKVEP_DELTA_ROUTE_INS_DIRECT_FALLBACK
                                          : DUCKVEP_DELTA_ROUTE_INDEL_DIRECT_FALLBACK;
    }
}

DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_for_annotation(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta) {

    duckvep_sequence_delta_fill_for_annotation_trace(kind, transcripts, exons, seq, v,
                                                     variant_idx, tx_idx, pos, strand,
                                                     scratch, prepared_event, NULL, delta);
}
