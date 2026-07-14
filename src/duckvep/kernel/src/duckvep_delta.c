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

static char delta_orient_genomic_base(char genomic_base,
                                      int8_t transcript_strand);
static int delta_allele_slice_ok(const duckvep_variant_batch_t *v,
                                 uint32_t variant_idx);
static int delta_cds_slice(const duckvep_sequence_pool_t *seq,
                           size_t tx_idx, const uint8_t **cds_seq,
                           size_t *cds_len);
static int delta_transcript_has_sequence(const duckvep_sequence_pool_t *seq,
                                         size_t tx_idx);

typedef struct delta_complete_transcript_sequence {
    const uint8_t *pre_cds;
    const uint8_t *cds;
    const uint8_t *post_cds;
    size_t pre_cds_len;
    size_t cds_len;
    size_t post_cds_len;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint8_t phase_offset;
} delta_complete_transcript_sequence_t;

typedef struct delta_feature_span_edit {
    const uint8_t *ref;
    const uint8_t *alt;
    uint16_t ref_len;
    uint16_t alt_len;
    uint32_t cdna_start;
    uint32_t cdna_end;
    size_t edit_start0;
} delta_feature_span_edit_t;

typedef struct delta_sequence_edit_view {
    const uint8_t *left;
    const uint8_t *right;
    const uint8_t *alt;
    size_t left_len;
    size_t right_len;
    size_t alt_len;
    size_t edit_start;
    size_t edit_ref_len;
    size_t edited_len;
    int8_t strand;
} delta_sequence_edit_view_t;

static int delta_coding_cdna_bounds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                         *coding_start_cdna,
    uint32_t                         *coding_end_cdna,
    uint32_t                         *coding_start_exon) {

    uint32_t start_genomic;
    uint32_t end_genomic;
    uint32_t ignored_exon;

    if (coding_start_cdna != NULL) *coding_start_cdna = 0u;
    if (coding_end_cdna != NULL) *coding_end_cdna = 0u;
    if (coding_start_exon != NULL) *coding_start_exon = 0u;
    if (transcripts == NULL || exons == NULL || coding_start_cdna == NULL ||
        coding_end_cdna == NULL || tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL || transcripts->cds_start1 == NULL ||
        transcripts->cds_end1 == NULL || transcripts->cds_start1[tx_idx] == 0u) {
        return 0;
    }
    start_genomic = transcripts->strand[tx_idx] > 0
        ? transcripts->cds_start1[tx_idx] : transcripts->cds_end1[tx_idx];
    end_genomic = transcripts->strand[tx_idx] > 0
        ? transcripts->cds_end1[tx_idx] : transcripts->cds_start1[tx_idx];
    if (!duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, start_genomic,
            coding_start_cdna, coding_start_exon != NULL
                ? coding_start_exon : &ignored_exon) ||
        !duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, end_genomic,
            coding_end_cdna, &ignored_exon) ||
        *coding_end_cdna < *coding_start_cdna) {
        return 0;
    }
    return 1;
}

/* Returns 1 for a complete view, 0 when the model deliberately lacks complete
 * flanks, and -1 for an internally inconsistent view. Model-open validation
 * proves these invariants for production calls; the checks keep the borrowed
 * pure-C entry points fail-closed as well. */
static int delta_complete_transcript_sequence_load(
    const duckvep_transcript_model_t       *transcripts,
    const duckvep_exon_model_t             *exons,
    const duckvep_sequence_pool_t          *seq,
    size_t                                  tx_idx,
    delta_complete_transcript_sequence_t   *out) {

    const uint8_t *cds;
    size_t cds_len;
    uint64_t pre_off;
    uint64_t post_off;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint32_t coding_start_exon;
    size_t exon_off;
    size_t exon_count;
    uint64_t expected_post;
    int8_t phase;

    if (out != NULL) memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || seq == NULL || out == NULL ||
        tx_idx >= transcripts->transcript_count ||
        !delta_cds_slice(seq, tx_idx, &cds, &cds_len)) {
        return -1;
    }
    if (seq->flanks_complete == 0u) return 0;
    if (seq->pre_cds_offset == NULL || seq->pre_cds_length == NULL ||
        seq->post_cds_offset == NULL || seq->post_cds_length == NULL ||
        !delta_coding_cdna_bounds(transcripts, exons, tx_idx,
                                  &coding_start_cdna, &coding_end_cdna,
                                  &coding_start_exon) ||
        exons->phase == NULL || coding_start_exon >= exons->exon_count) {
        return -1;
    }
    pre_off = seq->pre_cds_offset[tx_idx];
    post_off = seq->post_cds_offset[tx_idx];
    if (pre_off > (uint64_t)seq->flank_bytes_len ||
        (uint64_t)seq->pre_cds_length[tx_idx] >
            (uint64_t)seq->flank_bytes_len - pre_off ||
        post_off > (uint64_t)seq->flank_bytes_len ||
        (uint64_t)seq->post_cds_length[tx_idx] >
            (uint64_t)seq->flank_bytes_len - post_off ||
        ((seq->pre_cds_length[tx_idx] != 0u ||
          seq->post_cds_length[tx_idx] != 0u) && seq->flank_bytes == NULL)) {
        return -1;
    }
    exon_off = (size_t)transcripts->exon_offset[tx_idx];
    exon_count = (size_t)transcripts->exon_count[tx_idx];
    if (exon_count == 0u || exon_off > exons->exon_count ||
        exon_count > exons->exon_count - exon_off ||
        exons->cdna_end1 == NULL) {
        return -1;
    }
    phase = exons->phase[coding_start_exon];
    out->phase_offset = phase > 0 ? (uint8_t)phase : 0u;
    expected_post = (uint64_t)exons->cdna_end1[exon_off + exon_count - 1u] -
                    (uint64_t)coding_end_cdna;
    if ((uint64_t)seq->pre_cds_length[tx_idx] !=
            (uint64_t)coding_start_cdna - 1u ||
        (uint64_t)seq->post_cds_length[tx_idx] != expected_post ||
        (uint64_t)cds_len !=
            (uint64_t)coding_end_cdna - (uint64_t)coding_start_cdna + 1u +
                (uint64_t)out->phase_offset) {
        return -1;
    }
    out->pre_cds = seq->pre_cds_length[tx_idx] != 0u
        ? seq->flank_bytes + (size_t)pre_off : NULL;
    out->cds = cds;
    out->post_cds = seq->post_cds_length[tx_idx] != 0u
        ? seq->flank_bytes + (size_t)post_off : NULL;
    out->pre_cds_len = (size_t)seq->pre_cds_length[tx_idx];
    out->cds_len = cds_len;
    out->post_cds_len = (size_t)seq->post_cds_length[tx_idx];
    out->coding_start_cdna = coding_start_cdna;
    out->coding_end_cdna = coding_end_cdna;
    return 1;
}

static char delta_complete_transcript_cdna_base(
    const delta_complete_transcript_sequence_t *sequence,
    uint32_t                                    cdna_pos) {

    size_t offset;

    if (sequence == NULL || cdna_pos == 0u) return '\0';
    if ((size_t)cdna_pos <= sequence->pre_cds_len) {
        return delta_norm_base((char)sequence->pre_cds[(size_t)cdna_pos - 1u]);
    }
    if (cdna_pos >= sequence->coding_start_cdna &&
        cdna_pos <= sequence->coding_end_cdna) {
        offset = (size_t)sequence->phase_offset +
                 (size_t)(cdna_pos - sequence->coding_start_cdna);
        if (offset >= sequence->cds_len) return '\0';
        return delta_norm_base((char)sequence->cds[offset]);
    }
    if (cdna_pos > sequence->coding_end_cdna) {
        offset = (size_t)(cdna_pos - sequence->coding_end_cdna - 1u);
        if (offset >= sequence->post_cds_len) return '\0';
        return delta_norm_base((char)sequence->post_cds[offset]);
    }
    return '\0';
}

static int delta_feature_span_edit_load(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    delta_feature_span_edit_t        *out) {

    uint32_t first_cdna;
    uint32_t last_cdna;
    uint64_t genomic_span;
    uint64_t cdna_span;

    if (out != NULL) memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || v == NULL || event == NULL ||
        out == NULL || variant_idx >= v->count ||
        tx_idx >= transcripts->transcript_count || event->interbase) {
        return 0;
    }
    if (!duckvep_event_feature_alleles(
            v, variant_idx, event, &out->ref, &out->ref_len,
            &out->alt, &out->alt_len)) {
        return 0;
    }
    if (out->ref_len == 0u || out->ref_len == out->alt_len ||
        event->feature_start1 == 0u ||
        event->feature_end1 < event->feature_start1) {
        return 0;
    }
    genomic_span = (uint64_t)event->feature_end1 -
                   (uint64_t)event->feature_start1 + 1u;
    if (genomic_span != (uint64_t)out->ref_len ||
        !duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, event->feature_start1,
            &first_cdna, NULL) ||
        !duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, event->feature_end1,
            &last_cdna, NULL)) {
        return 0;
    }
    out->cdna_start = first_cdna < last_cdna ? first_cdna : last_cdna;
    out->cdna_end = first_cdna > last_cdna ? first_cdna : last_cdna;
    cdna_span = (uint64_t)out->cdna_end - (uint64_t)out->cdna_start + 1u;
    if (cdna_span != (uint64_t)out->ref_len) return 0;
    out->edit_start0 = (size_t)out->cdna_start - 1u;
    return 1;
}

static char delta_feature_allele_base(
    const uint8_t *allele,
    size_t         length,
    size_t         index,
    int8_t         strand) {

    size_t oriented_index;

    if (allele == NULL || index >= length) return '\0';
    oriented_index = strand > 0 ? index : length - 1u - index;
    return delta_orient_genomic_base((char)allele[oriented_index], strand);
}

static int delta_feature_span_edit_validate(
    const delta_complete_transcript_sequence_t *sequence,
    const delta_feature_span_edit_t             *edit,
    int8_t                                       strand) {

    size_t i;

    if (sequence == NULL || edit == NULL) return 0;
    for (i = 0u; i < (size_t)edit->ref_len; i++) {
        char uploaded = delta_feature_allele_base(
            edit->ref, (size_t)edit->ref_len, i, strand);
        char reference = delta_complete_transcript_cdna_base(
            sequence, edit->cdna_start + (uint32_t)i);
        if (uploaded == '\0' || uploaded == 'N' || reference == '\0' ||
            reference == 'N') return -1;
        if (uploaded != reference) return 0;
    }
    for (i = 0u; i < (size_t)edit->alt_len; i++) {
        char uploaded = delta_feature_allele_base(
            edit->alt, (size_t)edit->alt_len, i, strand);
        if (uploaded == '\0' || uploaded == 'N') return -1;
    }
    return 1;
}

static int delta_sequence_edit_view_open(
    const uint8_t              *left,
    size_t                      left_len,
    const uint8_t              *right,
    size_t                      right_len,
    const uint8_t              *alt,
    size_t                      alt_len,
    size_t                      edit_ref_len,
    size_t                      edit_start,
    int8_t                      strand,
    delta_sequence_edit_view_t *out) {

    size_t reference_len;

    if (out != NULL) memset(out, 0, sizeof *out);
    if (out == NULL ||
        (left_len != 0u && left == NULL) ||
        (right_len != 0u && right == NULL) ||
        (alt_len != 0u && alt == NULL) ||
        (strand != (int8_t)1 && strand != (int8_t)-1) ||
        left_len > SIZE_MAX - right_len) return 0;
    reference_len = left_len + right_len;
    if (edit_start > reference_len ||
        edit_ref_len > reference_len - edit_start ||
        alt_len > SIZE_MAX - (reference_len - edit_ref_len)) {
        return 0;
    }
    out->left = left;
    out->right = right;
    out->alt = alt;
    out->left_len = left_len;
    out->right_len = right_len;
    out->alt_len = alt_len;
    out->edit_start = edit_start;
    out->edit_ref_len = edit_ref_len;
    out->edited_len = reference_len - edit_ref_len + alt_len;
    out->strand = strand;
    return 1;
}

static char delta_sequence_reference_base(
    const delta_sequence_edit_view_t *sequence,
    size_t                            position) {

    if (sequence == NULL) return '\0';
    if (position < sequence->left_len) {
        return delta_norm_base((char)sequence->left[position]);
    }
    position -= sequence->left_len;
    if (position >= sequence->right_len) return '\0';
    return delta_norm_base((char)sequence->right[position]);
}

static char delta_sequence_edited_base(
    const delta_sequence_edit_view_t *sequence,
    size_t                            position) {

    size_t reference_position;

    if (sequence == NULL || position >= sequence->edited_len) return '\0';
    if (position < sequence->edit_start) {
        return delta_sequence_reference_base(sequence, position);
    }
    if (position - sequence->edit_start < sequence->alt_len) {
        return delta_feature_allele_base(
            sequence->alt, sequence->alt_len,
            position - sequence->edit_start, sequence->strand);
    }
    reference_position = position - sequence->alt_len +
                         sequence->edit_ref_len;
    return delta_sequence_reference_base(sequence, reference_position);
}

static int delta_sequence_edited_range_equal(
    const delta_sequence_edit_view_t *sequence,
    size_t                            edited_offset,
    const uint8_t                    *reference,
    size_t                            length) {

    size_t i;

    if (sequence == NULL || (length != 0u && reference == NULL) ||
        edited_offset > sequence->edited_len ||
        length > sequence->edited_len - edited_offset) {
        return 0;
    }
    for (i = 0u; i < length; i++) {
        char edited = delta_sequence_edited_base(sequence, edited_offset + i);
        char expected = delta_norm_base((char)reference[i]);
        if (edited == '\0' || expected == '\0' || edited != expected) return 0;
    }
    return 1;
}

static int delta_sequence_start_codon_is_atg(
    const delta_sequence_edit_view_t *sequence,
    size_t                            start) {

    return delta_sequence_edited_base(sequence, start) == 'A' &&
           delta_sequence_edited_base(sequence, start + 1u) == 'T' &&
           delta_sequence_edited_base(sequence, start + 2u) == 'G';
}

/* Exact string tests from VEP 116::_ins_del_start_altered. The first test is
 * intentionally stronger than “the start codon is still ATG”: VEP also requires
 * the original 5-prime UTR to remain byte-identical. If that test fails, VEP
 * compares the complete original translateable sequence with the same-length
 * suffix of the edited UTR+CDS sequence. */
static int delta_sequence_start_altered(
    const delta_sequence_edit_view_t *sequence) {

    if (sequence == NULL) return 1;
    if (sequence->left_len != 0u &&
        delta_sequence_edited_range_equal(
            sequence, 0u, sequence->left, sequence->left_len) &&
        delta_sequence_start_codon_is_atg(sequence, sequence->left_len)) {
        return 0;
    }
    if (sequence->edited_len < sequence->right_len) return 1;
    return !delta_sequence_edited_range_equal(
        sequence, sequence->edited_len - sequence->right_len,
        sequence->right, sequence->right_len);
}

/* VEP calls its inversion helper for ordinary TranscriptVariationAlleles too.
 * That helper edits the same UTR+CDS string and merely asks whether the codon at
 * the original start offset is still ATG. Keeping this independent from the
 * retained predicate reproduces the observable lost+retained combination. */
static int delta_sequence_start_offset_altered(
    const delta_sequence_edit_view_t *sequence) {

    if (sequence == NULL || sequence->left_len == 0u) return 0;
    return !delta_sequence_start_codon_is_atg(sequence, sequence->left_len);
}

static int delta_sequence_stop_altered(
    const delta_sequence_edit_view_t *sequence,
    duckvep_codon_table_t             codon_table) {

    char codon[4];
    size_t i;

    if (sequence == NULL || sequence->left_len < 3u ||
        sequence->edited_len < sequence->left_len) {
        return 1;
    }
    for (i = 0u; i < 3u; i++) {
        codon[i] = delta_sequence_edited_base(
            sequence, sequence->left_len - 3u + i);
    }
    codon[3] = '\0';
    return duckvep_translate_codon(codon, codon_table) != '*';
}

static int delta_feature_span_overlaps_cds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          feature_start1,
    uint32_t                          feature_end1) {

    uint32_t exon_offset;
    uint16_t exon_count;
    uint16_t i;
    uint32_t cds_start1;
    uint32_t cds_end1;

    if (transcripts == NULL || exons == NULL ||
        transcripts->cds_start1 == NULL || transcripts->cds_end1 == NULL ||
        transcripts->exon_offset == NULL || transcripts->exon_count == NULL ||
        exons->start1 == NULL || exons->end1 == NULL ||
        tx_idx >= transcripts->transcript_count || feature_start1 == 0u ||
        feature_end1 < feature_start1) {
        return 0;
    }
    cds_start1 = transcripts->cds_start1[tx_idx];
    cds_end1 = transcripts->cds_end1[tx_idx];
    exon_offset = transcripts->exon_offset[tx_idx];
    exon_count = transcripts->exon_count[tx_idx];
    if (cds_start1 == 0u || cds_end1 < cds_start1 ||
        exon_offset > exons->exon_count ||
        exon_count > exons->exon_count - exon_offset) {
        return 0;
    }
    for (i = 0u; i < exon_count; i++) {
        size_t exon_idx = (size_t)exon_offset + (size_t)i;
        uint32_t lo = exons->start1[exon_idx];
        uint32_t hi = exons->end1[exon_idx];

        if (lo < cds_start1) lo = cds_start1;
        if (lo < feature_start1) lo = feature_start1;
        if (hi > cds_end1) hi = cds_end1;
        if (hi > feature_end1) hi = feature_end1;
        if (lo <= hi) return 1;
    }
    return 0;
}

/* VEP's length-changing mapper-gap state is not a partial CDS edit. When either
 * transcript-oriented endpoint maps to a Gap, cds_start/cds_end or their cDNA
 * counterparts are undefined and ordinary codon/peptide construction stops.
 * `coding_unknown` then represents any mapped CDS portion. One independent
 * exception runs first: `_overlaps_stop_codon_cil` tests the genomic feature span,
 * while `_ins_del_stop_altered_cil` returns false on the missing endpoint, yielding
 * stop_retained_variant without editing sequence. */
static int delta_feature_length_change_mapper_gap_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            strand,
    const duckvep_event_t            *event,
    duckvep_sequence_delta_t         *delta) {

    uint32_t first_genomic1;
    uint32_t last_genomic1;
    uint32_t first_cdna1;
    uint32_t last_cdna1;
    uint32_t stop_start1;
    uint32_t stop_end1;
    uint64_t tx_flags = 0u;
    int first_is_coordinate;
    int last_is_coordinate;
    int overlaps_cds;
    int overlaps_stop = 0;

    if (transcripts == NULL || exons == NULL || v == NULL || event == NULL ||
        delta == NULL || transcripts->strand == NULL ||
        transcripts->cds_start1 == NULL || transcripts->cds_end1 == NULL ||
        variant_idx >= v->count || tx_idx >= transcripts->transcript_count ||
        (strand != (int8_t)1 && strand != (int8_t)-1) || event->interbase ||
        event->ref_diff_length == event->alt_diff_length ||
        event->feature_start1 == 0u ||
        event->feature_end1 < event->feature_start1) {
        return 0;
    }

    first_genomic1 = strand > 0 ? event->feature_start1 : event->feature_end1;
    last_genomic1 = strand > 0 ? event->feature_end1 : event->feature_start1;
    first_is_coordinate = duckvep_project_genomic_to_cdna(
        transcripts, exons, tx_idx, first_genomic1, &first_cdna1, NULL);
    last_is_coordinate = duckvep_project_genomic_to_cdna(
        transcripts, exons, tx_idx, last_genomic1, &last_cdna1, NULL);
    if (first_is_coordinate && last_is_coordinate) return 0;

    overlaps_cds = delta_feature_span_overlaps_cds(
        transcripts, exons, tx_idx,
        event->feature_start1, event->feature_end1);
    if (!overlaps_cds) return 0;

    if (transcripts->flags != NULL) tx_flags = transcripts->flags[tx_idx];
    if ((tx_flags & (uint64_t)DUCKVEP_TX_CDS_END_NF) == 0u) {
        if (strand > 0) {
            stop_end1 = transcripts->cds_end1[tx_idx];
            stop_start1 = stop_end1 >= 3u ? stop_end1 - 2u : 1u;
        } else {
            stop_start1 = transcripts->cds_start1[tx_idx];
            stop_end1 = stop_start1 <= UINT32_MAX - 2u
                ? stop_start1 + 2u : UINT32_MAX;
        }
        overlaps_stop = event->feature_start1 <= stop_end1 &&
                        event->feature_end1 >= stop_start1;
    }

    memset(delta, 0, sizeof *delta);
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = -1;
    if (overlaps_stop) delta->stop_retained = 1u;
    else delta->coding_unknown = 1u;
    delta->valid = 1u;
    delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    return 1;
}

/* Handle length-changing uploaded features whose REF span crosses one cDNA/CDS
 * boundary. VEP cannot form its ordinary peptide alleles for this mapper-gap
 * shape. Its start and stop predicates instead rebuild UTR+translateable strings,
 * so this path uses the complete cold transcript flanks and never constructs a
 * temporary sequence. Return 1 once the boundary shape is recognized: missing or
 * malformed sequence is then an explicit result, not permission to retry a
 * biologically different CDS-only edit. */
static int delta_feature_length_change_boundary_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            strand,
    const duckvep_event_t            *event,
    duckvep_sequence_delta_t         *delta) {

    delta_complete_transcript_sequence_t transcript_sequence;
    delta_feature_span_edit_t edit;
    delta_sequence_edit_view_t sequence_edit;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint64_t tx_flags = 0u;
    int crosses_start;
    int crosses_stop;
    int sequence_status;
    int allele_status;

    if (delta == NULL ||
        !delta_feature_span_edit_load(
            transcripts, exons, v, variant_idx, tx_idx, event, &edit) ||
        !delta_coding_cdna_bounds(
            transcripts, exons, tx_idx, &coding_start_cdna,
            &coding_end_cdna, NULL)) {
        return 0;
    }
    crosses_start = edit.cdna_start < coding_start_cdna &&
                    edit.cdna_end >= coding_start_cdna &&
                    (uint64_t)edit.cdna_start <=
                        (uint64_t)coding_start_cdna + 2u;
    crosses_stop = edit.cdna_start <= coding_end_cdna &&
                   edit.cdna_end > coding_end_cdna &&
                   (uint64_t)edit.cdna_end + 2u >=
                       (uint64_t)coding_end_cdna;
    if (!crosses_start && !crosses_stop) return 0;
    /* A feature spanning both ends needs the general transcript edit-set path;
     * neither VEP helper alone defines all of its interactions. */
    if (crosses_start && crosses_stop) return 0;

    memset(delta, 0, sizeof *delta);
    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = -1;
    delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
    if (transcripts->flags != NULL) tx_flags = transcripts->flags[tx_idx];
    if ((tx_flags & (uint64_t)DUCKVEP_TX_CDS_START_NF) != 0u) {
        crosses_start = 0;
    }
    if ((tx_flags & (uint64_t)DUCKVEP_TX_CDS_END_NF) != 0u) {
        crosses_stop = 0;
    }
    if (!crosses_start && !crosses_stop) {
        delta->coding_unknown = 1u;
        delta->valid = 1u;
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
        return 1;
    }
    if (!delta_transcript_has_sequence(seq, tx_idx)) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_MISSING;
        return 1;
    }
    sequence_status = delta_complete_transcript_sequence_load(
        transcripts, exons, seq, tx_idx, &transcript_sequence);
    if (sequence_status == 0) {
        delta->sequence_status =
            (uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_FLANK;
        return 1;
    }
    if (sequence_status < 0) return 1;

    allele_status = delta_feature_span_edit_validate(
        &transcript_sequence, &edit, strand);
    if (allele_status < 0) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS;
        return 1;
    }
    if (allele_status == 0) {
        delta->sequence_status =
            (uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH;
        return 1;
    }

    if (crosses_start) {
        int altered;
        int inversion_helper_altered;

        if (!delta_sequence_edit_view_open(
                transcript_sequence.pre_cds,
                transcript_sequence.pre_cds_len,
                transcript_sequence.cds,
                transcript_sequence.cds_len,
                edit.alt, (size_t)edit.alt_len, (size_t)edit.ref_len,
                edit.edit_start0, strand, &sequence_edit)) {
            return 1;
        }
        altered = delta_sequence_start_altered(&sequence_edit);
        inversion_helper_altered =
            delta_sequence_start_offset_altered(&sequence_edit);
        delta->start_retained = (uint8_t)!altered;
        /* A UTR/CDS mapper-gap has no peptide/codon alleles, so VEP's in-frame
         * guards are false. The independent inversion helper still runs. */
        delta->start_lost = (uint8_t)(altered || inversion_helper_altered);
    } else {
        size_t edit_start;
        int altered;
        duckvep_codon_table_t codon_table;

        if (edit.cdna_start < transcript_sequence.coding_start_cdna ||
            seq == NULL || seq->codon_table == NULL) {
            return 1;
        }
        edit_start = (size_t)transcript_sequence.phase_offset +
                     (size_t)(edit.cdna_start -
                              transcript_sequence.coding_start_cdna);
        if (!delta_sequence_edit_view_open(
                transcript_sequence.cds, transcript_sequence.cds_len,
                transcript_sequence.post_cds,
                transcript_sequence.post_cds_len,
                edit.alt, (size_t)edit.alt_len, (size_t)edit.ref_len,
                edit_start, strand, &sequence_edit) ||
            (seq->codon_table[tx_idx] !=
                 (uint8_t)DUCKVEP_CODON_TABLE_STANDARD &&
             seq->codon_table[tx_idx] !=
                 (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO)) {
            return 1;
        }
        codon_table = (duckvep_codon_table_t)seq->codon_table[tx_idx];
        altered = delta_sequence_stop_altered(&sequence_edit, codon_table);
        delta->stop_lost = (uint8_t)altered;
        delta->stop_retained = (uint8_t)!altered;
    }

    delta->valid = 1u;
    delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
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
    /* The host-neutral builder has no separate UTR input. Treat that as a
     * complete empty 5-prime UTR; model-backed callers replace this fact with
     * the imported transcript flank before consequence evaluation. */
    tmp.pre_cds_complete = 1u;
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
    case DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK:
        return (uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_FLANK;
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
    ctx->pre_cds_complete = 0u;
    if (seq->flanks_complete != 0u) {
        uint64_t offset;
        size_t length;

        if (seq->pre_cds_offset == NULL || seq->pre_cds_length == NULL) {
            return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
        }
        offset = seq->pre_cds_offset[tx_idx];
        length = (size_t)seq->pre_cds_length[tx_idx];
        if (offset > (uint64_t)seq->flank_bytes_len ||
            (uint64_t)length > (uint64_t)seq->flank_bytes_len - offset) {
            return DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE;
        }
        if (length != 0u) {
            if (seq->flank_bytes == NULL) {
                return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
            }
            ctx->pre_cds_bases = seq->flank_bytes + (size_t)offset;
            ctx->pre_cds_length = length;
        }
        ctx->pre_cds_complete = 1u;
    } else {
        uint32_t coding_start_cdna;
        uint32_t coding_end_cdna;

        /* A legacy model still proves an empty 5-prime UTR when its CDS starts
         * at cDNA base one. No sequence bytes are missing in that case. */
        if (delta_coding_cdna_bounds(
                transcripts, exons, tx_idx, &coding_start_cdna,
                &coding_end_cdna, NULL) && coding_start_cdna == 1u) {
            ctx->pre_cds_complete = 1u;
        }
    }
    if (seq->post_cds_offset != NULL && seq->post_cds_length != NULL) {
        uint64_t offset = seq->post_cds_offset[tx_idx];
        size_t length = (size_t)seq->post_cds_length[tx_idx];

        if (offset > (uint64_t)seq->flank_bytes_len ||
            (uint64_t)length > (uint64_t)seq->flank_bytes_len - offset) {
            return DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE;
        }
        if (length != 0u) {
            if (seq->flank_bytes == NULL) {
                return DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG;
            }
            ctx->post_cds_bases = seq->flank_bytes + (size_t)offset;
            ctx->post_cds_length = length;
        }
    }
    ctx->post_cds_complete = seq->flanks_complete != 0u;
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

static int delta_inframe_insertion_site(uint32_t before_cds, size_t cds_len);

static int delta_norm_window_prefix_or_suffix(
    const uint8_t *outer, size_t outer_len,
    const uint8_t *inner, size_t inner_len) {

    if (outer == NULL || inner == NULL || inner_len > outer_len) return 0;
    if (inner_len == 0u) return 1;
    return delta_norm_span_equal(outer, inner, inner_len) ||
           delta_norm_span_equal(
               outer + outer_len - inner_len, inner, inner_len);
}

typedef struct delta_vep_local_peptide {
    size_t peptide_offset;
    size_t reference_span_length;
    size_t ref_nt_length;
    size_t alt_nt_length;
    size_t ref_whole_length;
    size_t alt_whole_length;
    size_t ref_length;
    size_t alt_length;
    uint8_t ref_partial_x;
    uint8_t alt_partial_x;
} delta_vep_local_peptide_t;

/* TranscriptVariationAllele::codon rounds the feature to reference codons,
 * changes that nucleotide length by allele_length - feature_length, and takes
 * the resulting slice from the edited CDS. TVA::peptide then appends X for a
 * trailing partial codon, except when the complete peptide is exactly "*".
 * Keep this as the one allocation-free representation of those local strings;
 * start, stop, frame, and insertion predicates all inspect the same bytes. */
static int delta_context_vep_local_peptide_open(
    const duckvep_coding_context_t *ctx,
    delta_vep_local_peptide_t      *view) {

    uint64_t first_cds;
    uint64_t last_cds;
    uint64_t translation_start;
    uint64_t translation_end;
    uint64_t ref_codons;
    uint64_t nt_offset64;
    uint64_t ref_nt_length64;
    uint64_t ref_whole_length64;
    int64_t requested_alt_nt;
    size_t nt_offset;

    if (view != NULL) memset(view, 0, sizeof *view);
    if (ctx == NULL || view == NULL || !ctx->has_single_edit ||
        ctx->single_edit_cds_start == 0u || ctx->length_diff == 0 ||
        ctx->length_diff == INT64_MIN || ctx->ref_cds == NULL ||
        ctx->alt_cds == NULL || ctx->ref_peptide == NULL ||
        ctx->alt_peptide == NULL) {
        return 0;
    }

    first_cds = (uint64_t)ctx->single_edit_cds_start;
    translation_start = ((first_cds - 1u) / 3u) + 1u;
    if (ctx->single_edit_ref_len != 0u) {
        last_cds = first_cds + (uint64_t)ctx->single_edit_ref_len - 1u;
        if (last_cds < first_cds) return 0;
        translation_end = ((last_cds - 1u) / 3u) + 1u;
        ref_codons = translation_end - translation_start + 1u;
    } else if (first_cds > 1u) {
        last_cds = first_cds - 1u;
        translation_end = ((last_cds - 1u) / 3u) + 1u;
        ref_codons = translation_end >= translation_start
            ? translation_end - translation_start + 1u : 0u;
    } else {
        ref_codons = 0u;
    }

    nt_offset64 = (translation_start - 1u) * 3u;
    ref_nt_length64 = ref_codons * 3u;
    if (nt_offset64 > (uint64_t)SIZE_MAX ||
        ref_codons > (uint64_t)SIZE_MAX ||
        ref_nt_length64 > (uint64_t)INT64_MAX ||
        nt_offset64 > (uint64_t)ctx->ref_cds_len) {
        return 0;
    }
    if (ref_nt_length64 > (uint64_t)ctx->ref_cds_len - nt_offset64) {
        /* A CDS_END_NF feature can end one or two nucleotides into the last
         * rounded codon. TVA::codon returns those available bases and
         * TVA::peptide appends X; it does not reject the feature because the
         * requested reference codon extends past the incomplete CDS. */
        if ((ctx->ref_cds_len % 3u) == 0u) return 0;
        ref_nt_length64 = (uint64_t)ctx->ref_cds_len - nt_offset64;
    }
    ref_whole_length64 = ref_nt_length64 / 3u;
    if (nt_offset64 / 3u > (uint64_t)ctx->ref_peptide_len ||
        ref_whole_length64 >
            (uint64_t)ctx->ref_peptide_len - nt_offset64 / 3u) {
        return 0;
    }
    if (ctx->length_diff > 0 &&
        ctx->length_diff > INT64_MAX - (int64_t)ref_nt_length64) {
        return 0;
    }
    requested_alt_nt = (int64_t)ref_nt_length64 + ctx->length_diff;
    if (requested_alt_nt < 0) return 0;

    nt_offset = (size_t)nt_offset64;
    if (nt_offset > ctx->alt_cds_len) return 0;
    view->peptide_offset = nt_offset / 3u;
    view->ref_nt_length = (size_t)ref_nt_length64;
    view->alt_nt_length = (size_t)requested_alt_nt;
    if (view->alt_nt_length > ctx->alt_cds_len - nt_offset) {
        view->alt_nt_length = ctx->alt_cds_len - nt_offset;
    }
    view->ref_whole_length = view->ref_nt_length / 3u;
    view->alt_whole_length = view->alt_nt_length / 3u;
    if (view->peptide_offset > ctx->alt_peptide_len ||
        view->alt_whole_length >
            ctx->alt_peptide_len - view->peptide_offset) {
        return 0;
    }
    view->ref_partial_x = (uint8_t)(
        (view->ref_nt_length % 3u) != 0u &&
        !(view->ref_whole_length == 1u &&
          ctx->ref_peptide[view->peptide_offset] == (uint8_t)'*'));
    view->alt_partial_x = (uint8_t)(
        (view->alt_nt_length % 3u) != 0u &&
        !(view->alt_whole_length == 1u &&
          ctx->alt_peptide[view->peptide_offset] == (uint8_t)'*'));
    if (view->ref_whole_length > SIZE_MAX - (size_t)view->ref_partial_x ||
        view->alt_whole_length > SIZE_MAX - (size_t)view->alt_partial_x) {
        return 0;
    }
    view->ref_length = view->ref_whole_length + (size_t)view->ref_partial_x;
    view->alt_length = view->alt_whole_length + (size_t)view->alt_partial_x;
    view->reference_span_length = view->ref_length;
    return 1;
}

static int delta_context_local_translation_matches(
    const duckvep_coding_context_t  *ctx,
    const delta_vep_local_peptide_t *view,
    int                              alternate) {

    const uint8_t *cds;
    const uint8_t *peptide;
    size_t whole_length;
    size_t nt_offset;
    size_t i;

    if (ctx == NULL || view == NULL ||
        (ctx->codon_table != (uint8_t)DUCKVEP_CODON_TABLE_STANDARD &&
         ctx->codon_table != (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO) ||
        view->peptide_offset > SIZE_MAX / 3u) {
        return 0;
    }
    cds = alternate ? ctx->alt_cds : ctx->ref_cds;
    peptide = alternate ? ctx->alt_peptide : ctx->ref_peptide;
    whole_length = alternate ? view->alt_whole_length
                             : view->ref_whole_length;
    nt_offset = view->peptide_offset * 3u;
    for (i = 0u; i < whole_length; i++) {
        char codon[4];
        size_t b;

        for (b = 0u; b < 3u; b++) {
            codon[b] = delta_norm_base(
                (char)cds[nt_offset + i * 3u + b]);
            if (codon[b] != 'A' && codon[b] != 'C' &&
                codon[b] != 'G' && codon[b] != 'T') {
                return 0;
            }
        }
        codon[3] = '\0';
        if (peptide[view->peptide_offset + i] !=
            (uint8_t)duckvep_translate_codon(
                codon, (duckvep_codon_table_t)ctx->codon_table)) {
            return 0;
        }
    }
    return 1;
}

static uint8_t delta_context_vep_local_peptide_base(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view,
    int                              alternate,
    size_t                           index) {

    size_t whole_length;
    size_t length;
    uint8_t partial_x;
    const uint8_t *peptide;

    if (ctx == NULL || view == NULL) return 0u;
    whole_length = alternate ? view->alt_whole_length
                             : view->ref_whole_length;
    length = alternate ? view->alt_length : view->ref_length;
    partial_x = alternate ? view->alt_partial_x : view->ref_partial_x;
    peptide = alternate ? ctx->alt_peptide : ctx->ref_peptide;
    if (index >= length) return 0u;
    if (index == whole_length && partial_x) return (uint8_t)'X';
    return peptide[view->peptide_offset + index];
}

static int delta_context_vep_local_contains(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view,
    int                              alternate,
    uint8_t                          needle,
    size_t                          *first_index) {

    size_t length = alternate ? view->alt_length : view->ref_length;
    size_t i;
    if (first_index != NULL) *first_index = 0u;
    for (i = 0u; i < length; i++) {
        if (delta_context_vep_local_peptide_base(
                ctx, view, alternate, i) == needle) {
            if (first_index != NULL) *first_index = i;
            return 1;
        }
    }
    return 0;
}

/* The only edge comparison used by VEP's start, insertion, and protein-altering
 * rules: is the complete reference local peptide a prefix or suffix of the
 * selected alternate prefix? An empty reference matches both edges. */
static int delta_context_vep_alt_preserves_ref_edge(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view,
    size_t                           alt_length) {

    size_t i;
    int prefix = 1;
    int suffix = 1;

    if (ctx == NULL || view == NULL || alt_length > view->alt_length ||
        view->ref_length > alt_length) {
        return 0;
    }
    for (i = 0u; i < view->ref_length; i++) {
        uint8_t reference = delta_context_vep_local_peptide_base(
            ctx, view, 0, i);
        if (delta_context_vep_local_peptide_base(
                ctx, view, 1, i) != reference) {
            prefix = 0;
        }
        if (delta_context_vep_local_peptide_base(
                ctx, view, 1, alt_length - view->ref_length + i) != reference) {
            suffix = 0;
        }
    }
    return prefix || suffix;
}

/* VEP's insertion-before-start guard is narrower than its later insertion-shape
 * test: start_retained suppresses inframe_insertion only when the complete
 * reference peptide is the alternate suffix. A reference prefix means the
 * insertion is after the retained start and remains an in-frame insertion. */
static int delta_context_vep_alt_has_ref_suffix(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view,
    size_t                           alt_length) {

    size_t i;

    if (ctx == NULL || view == NULL || alt_length > view->alt_length ||
        view->ref_length > alt_length) {
        return 0;
    }
    for (i = 0u; i < view->ref_length; i++) {
        uint8_t reference = delta_context_vep_local_peptide_base(
            ctx, view, 0, i);
        if (delta_context_vep_local_peptide_base(
                ctx, view, 1, alt_length - view->ref_length + i) != reference) {
            return 0;
        }
    }
    return 1;
}

/* Final start_lost predicate after the two UTR+CDS string tests. */
static int delta_context_start_peptide_altered(
    const duckvep_coding_context_t *ctx) {

    delta_vep_local_peptide_t view;
    if (!delta_context_vep_local_peptide_open(ctx, &view) ||
        view.ref_length == 0u || view.alt_length == 0u ||
        (view.alt_length == 1u &&
         delta_context_vep_local_peptide_base(
             ctx, &view, 1, 0u) == (uint8_t)'X')) {
        return 0;
    }
    return !delta_context_vep_alt_preserves_ref_edge(
        ctx, &view, view.alt_length);
}

/* Layer VEP's independent start predicates on a length-changing edit. VEP
 * evaluates two different views of the same edited UTR+CDS string: retained
 * compares the original CDS with the edited suffix, while its inversion helper
 * checks ATG at the original start offset. They can deliberately disagree and
 * emit both start_lost and start_retained_variant. */
static duckvep_context_delta_status_t delta_context_start_facts(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    int                            *overlaps_start,
    duckvep_sequence_delta_t       *delta) {

    delta_sequence_edit_view_t sequence_edit;
    const uint8_t *alt_edit = NULL;
    size_t edit_start0;
    size_t transcript_edit_start;
    int altered;
    int offset_altered;
    int peptide_altered;

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
    if (ctx->pre_cds_complete == 0u) {
        return DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK;
    }
    if (ctx->ref_cds == NULL || ctx->alt_cds == NULL ||
        (ctx->pre_cds_length != 0u && ctx->pre_cds_bases == NULL)) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    edit_start0 = (size_t)ctx->single_edit_cds_start - 1u;
    if (edit_start0 > ctx->ref_cds_len ||
        (size_t)ctx->single_edit_ref_len > ctx->ref_cds_len - edit_start0 ||
        edit_start0 > ctx->alt_cds_len ||
        (size_t)ctx->single_edit_alt_len > ctx->alt_cds_len - edit_start0 ||
        ctx->pre_cds_length > SIZE_MAX - edit_start0) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (ctx->single_edit_alt_len != 0u) {
        alt_edit = ctx->alt_cds + edit_start0;
    }
    transcript_edit_start = ctx->pre_cds_length + edit_start0;
    if (!delta_sequence_edit_view_open(
            ctx->pre_cds_bases, ctx->pre_cds_length,
            ctx->ref_cds, ctx->ref_cds_len,
            alt_edit, (size_t)ctx->single_edit_alt_len,
            (size_t)ctx->single_edit_ref_len, transcript_edit_start,
            (int8_t)1, &sequence_edit)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    altered = delta_sequence_start_altered(&sequence_edit);
    offset_altered = delta_sequence_start_offset_altered(&sequence_edit);
    delta->start_retained = (uint8_t)!altered;

    peptide_altered = delta_context_start_peptide_altered(ctx);
    if ((altered && !delta->inframe_insertion &&
         !delta->inframe_deletion) || offset_altered || peptide_altered) {
        delta->start_lost = 1u;
    }
    return DUCKVEP_CONTEXT_DELTA_OK;
}

static duckvep_context_delta_status_t delta_context_original_endpoint_stop_altered(
    const duckvep_coding_context_t *ctx,
    int                            *altered,
    uint8_t                        *alternate_aa) {

    char codon[4];
    size_t i;

    if (altered != NULL) *altered = 1;
    if (alternate_aa != NULL) *alternate_aa = 0u;
    if (ctx == NULL || altered == NULL || alternate_aa == NULL ||
        ctx->ref_cds == NULL ||
        ctx->alt_cds == NULL || ctx->ref_cds_len < 3u ||
        (ctx->ref_cds_len % 3u) != 0u ||
        (ctx->codon_table != (uint8_t)DUCKVEP_CODON_TABLE_STANDARD &&
         ctx->codon_table != (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO)) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    /* _ins_del_stop_altered first edits translateable CDS + the complete 3-prime
     * UTR. If that whole edited string is shorter than the original CDS, VEP
     * declares the original stop altered without attempting to translate a
     * codon at an offset that no longer exists. A legacy short-tail model cannot
     * prove that state and must continue to fail closed. */
    if (ctx->alt_cds_len > SIZE_MAX - ctx->post_cds_length) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (ctx->alt_cds_len + ctx->post_cds_length < ctx->ref_cds_len) {
        if (ctx->post_cds_complete == 0u) {
            return DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL;
        }
        *altered = 1;
        return DUCKVEP_CONTEXT_DELTA_OK;
    }

    for (i = 0u; i < 3u; i++) {
        size_t position = ctx->ref_cds_len - 3u + i;
        char base;

        if (position < ctx->alt_cds_len) {
            base = delta_norm_base((char)ctx->alt_cds[position]);
        } else {
            size_t tail_position = position - ctx->alt_cds_len;
            if (ctx->post_cds_bases == NULL ||
                tail_position >= ctx->post_cds_length) {
                return DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL;
            }
            base = delta_norm_base((char)ctx->post_cds_bases[tail_position]);
        }
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        codon[i] = base;
    }
    codon[3] = '\0';
    *alternate_aa = (uint8_t)duckvep_translate_codon(
        codon, (duckvep_codon_table_t)ctx->codon_table);
    *altered = *alternate_aa != (uint8_t)'*';
    return DUCKVEP_CONTEXT_DELTA_OK;
}

static uint8_t delta_context_vep_mutated_peptide_base(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view,
    size_t                           position) {

    size_t alt_end;
    size_t suffix_start;
    size_t reference_position;

    if (ctx == NULL || view == NULL ||
        view->peptide_offset > ctx->ref_peptide_len ||
        view->reference_span_length >
            ctx->ref_peptide_len - view->peptide_offset ||
        view->alt_length > SIZE_MAX - view->peptide_offset) {
        return 0u;
    }
    alt_end = view->peptide_offset + view->alt_length;
    suffix_start = view->peptide_offset + view->reference_span_length;
    if (position < view->peptide_offset) {
        return position < ctx->ref_peptide_len
            ? ctx->ref_peptide[position] : 0u;
    }
    if (position < alt_end) {
        return delta_context_vep_local_peptide_base(
            ctx, view, 1, position - view->peptide_offset);
    }
    reference_position = suffix_start + (position - alt_end);
    return reference_position < ctx->ref_peptide_len
        ? ctx->ref_peptide[reference_position] : 0u;
}

/* VariationEffect::ref_eq_alt_sequence. The local stop-index equality is its
 * common path. The full-peptide splice comparison preserves the odd insertion
 * case where the original peptide remains unchanged and a new trailing stop is
 * appended. */
static int delta_context_vep_ref_eq_alt_sequence(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view) {

    size_t ref_stop = 0u;
    size_t alt_stop = 0u;
    size_t reference_length;
    size_t replaced_length;
    size_t mutated_length;
    size_t suffix_length;
    size_t i;
    int ref_has_stop;
    int alt_has_stop;

    if (ctx == NULL || view == NULL ||
        (view->ref_length == 1u && view->alt_length == 1u &&
         delta_context_vep_local_peptide_base(ctx, view, 0, 0u) ==
             (uint8_t)'X' &&
         delta_context_vep_local_peptide_base(ctx, view, 1, 0u) ==
             (uint8_t)'X')) {
        return 0;
    }
    ref_has_stop = delta_context_vep_local_contains(
        ctx, view, 0, (uint8_t)'*', &ref_stop);
    alt_has_stop = delta_context_vep_local_contains(
        ctx, view, 1, (uint8_t)'*', &alt_stop);
    if (ref_has_stop && alt_has_stop && ref_stop == alt_stop) return 1;

    /* BaseTranscriptVariation::_peptide calls Transcript::translate, whose
     * sequence excludes the ordinary terminal stop. The local TVA peptide may
     * still contain that stop. VEP's full-string branch therefore compares the
     * mutated prefix against ref_peptide[0..terminal_stop), then asks whether
     * the next altered residue is '*'. */
    reference_length = ctx->ref_peptide_len;
    if (reference_length != 0u &&
        ctx->ref_peptide[reference_length - 1u] == (uint8_t)'*') {
        reference_length--;
    }
    if (view->peptide_offset > reference_length) {
        return 0;
    }
    replaced_length = view->reference_span_length;
    if (replaced_length > reference_length - view->peptide_offset) {
        replaced_length = reference_length - view->peptide_offset;
    }
    suffix_length = reference_length - view->peptide_offset - replaced_length;
    if (view->peptide_offset > SIZE_MAX - view->alt_length ||
        view->peptide_offset + view->alt_length > SIZE_MAX - suffix_length) {
        return 0;
    }
    mutated_length = view->peptide_offset + view->alt_length + suffix_length;
    if (mutated_length <= reference_length) return 0;
    for (i = 0u; i < reference_length; i++) {
        if (delta_context_vep_mutated_peptide_base(ctx, view, i) !=
            ctx->ref_peptide[i]) {
            return 0;
        }
    }
    return delta_context_vep_mutated_peptide_base(
        ctx, view, reference_length) == (uint8_t)'*';
}

static int delta_context_vep_local_equal(
    const duckvep_coding_context_t *ctx,
    const delta_vep_local_peptide_t *view) {

    size_t i;
    if (ctx == NULL || view == NULL ||
        view->ref_length != view->alt_length) return 0;
    for (i = 0u; i < view->ref_length; i++) {
        if (delta_context_vep_local_peptide_base(ctx, view, 0, i) !=
            delta_context_vep_local_peptide_base(ctx, view, 1, i)) {
            return 0;
        }
    }
    return 1;
}

/* One authority for every single, contiguous length-changing CDS edit. VEP
 * evaluates start, stop, frame, in-frame, and protein-altering predicates
 * independently over the same local codon/peptide strings. Terminal edits add
 * the original-end reconstruction; ordinary edits use the same local state
 * without that terminal-only predicate. */
static duckvep_context_delta_status_t delta_context_length_change_predicates(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    int                            *handled,
    duckvep_sequence_delta_t       *delta) {

    delta_vep_local_peptide_t view;
    duckvep_context_delta_status_t status;
    uint64_t edit_end = 0u;
    size_t terminal_start;
    size_t ref_stop = 0u;
    size_t alt_stop = 0u;
    size_t alt_shape_length;
    size_t nt_offset;
    size_t trim_prefix = 0u;
    size_t trim_suffix = 0u;
    size_t ref_remainder;
    size_t alt_remainder;
    int endpoint_altered = 1;
    uint8_t endpoint_aa = 0u;
    int ref_has_stop;
    int alt_has_stop;
    int ref_has_x;
    int alt_has_x;
    int start_overlaps;
    int frameshift;
    int inframe_deletion = 0;
    int preserves_ref;
    int terminal_complete;
    int overlaps_terminal = 0;

    if (handled != NULL) *handled = 0;
    if (ctx == NULL || handled == NULL || delta == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    /* The stored peptide contains only complete codons. CDS_END_NF may leave
     * one or two CDS bases beyond it; the local peptide view below represents
     * that legitimate VEP state with a synthetic trailing X. */
    if (!ctx->has_single_edit || ctx->applied_edits != 1u ||
        ctx->length_diff == 0 ||
        ctx->length_diff == INT64_MIN || ctx->ref_cds == NULL ||
        ctx->alt_cds == NULL || ctx->ref_peptide == NULL ||
        ctx->alt_peptide == NULL || ctx->ref_cds_len < 3u ||
        ctx->ref_peptide_len != ctx->ref_cds_len / 3u ||
        ctx->single_edit_cds_start == 0u) {
        return DUCKVEP_CONTEXT_DELTA_OK;
    }
    if (ctx->length_diff > 0) {
        uint64_t increase = (uint64_t)ctx->length_diff;
        if (increase > (uint64_t)SIZE_MAX ||
            ctx->ref_cds_len > SIZE_MAX - (size_t)increase ||
            ctx->alt_cds_len != ctx->ref_cds_len + (size_t)increase ||
            (uint64_t)ctx->single_edit_alt_len !=
                (uint64_t)ctx->single_edit_ref_len + increase) {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
    } else {
        uint64_t decrease = (uint64_t)(-(ctx->length_diff + 1)) + 1u;
        if (decrease > (uint64_t)SIZE_MAX ||
            (size_t)decrease > ctx->ref_cds_len ||
            ctx->alt_cds_len != ctx->ref_cds_len - (size_t)decrease ||
            (uint64_t)ctx->single_edit_ref_len !=
                (uint64_t)ctx->single_edit_alt_len + decrease) {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
    }
    terminal_start = ctx->ref_cds_len - 2u;
    terminal_complete =
        (tx_flags & (uint64_t)DUCKVEP_TX_CDS_END_NF) == 0u &&
        ctx->ref_peptide[ctx->ref_peptide_len - 1u] == (uint8_t)'*';
    if (ctx->single_edit_ref_len == 0u) {
        if ((uint64_t)ctx->single_edit_cds_start >
                (uint64_t)ctx->ref_cds_len) {
            return DUCKVEP_CONTEXT_DELTA_OK;
        }
        overlaps_terminal = terminal_complete &&
            (uint64_t)ctx->single_edit_cds_start > (uint64_t)terminal_start;
    } else {
        edit_end = (uint64_t)ctx->single_edit_cds_start +
                   (uint64_t)ctx->single_edit_ref_len - 1u;
        if ((uint64_t)ctx->single_edit_cds_start >
                (uint64_t)ctx->ref_cds_len ||
            edit_end > (uint64_t)ctx->ref_cds_len) {
            return DUCKVEP_CONTEXT_DELTA_OK;
        }
        overlaps_terminal = terminal_complete &&
            edit_end >= (uint64_t)terminal_start;
    }
    *handled = 1;
    if (!delta_context_vep_local_peptide_open(ctx, &view) ||
        view.peptide_offset > (size_t)INT32_MAX - 1u ||
        view.peptide_offset > SIZE_MAX / 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    nt_offset = view.peptide_offset * 3u;
    /* TVA::peptide also spells a trailing partial codon as X. That state is
     * meaningful and drives several VEP predicates, whereas X translated from
     * an N-bearing codon means the sequence fact is unknowable. Inspect the
     * nucleotide window so those two sources of X cannot be conflated. */
    if (nt_offset > ctx->ref_cds_len ||
        view.ref_nt_length > ctx->ref_cds_len - nt_offset ||
        nt_offset > ctx->alt_cds_len ||
        view.alt_nt_length > ctx->alt_cds_len - nt_offset ||
        !delta_norm_span_unambiguous(
            ctx->ref_cds + nt_offset, view.ref_nt_length) ||
        !delta_norm_span_unambiguous(
            ctx->alt_cds + nt_offset, view.alt_nt_length) ||
        !delta_context_local_translation_matches(ctx, &view, 0) ||
        !delta_context_local_translation_matches(ctx, &view, 1)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }

    ref_has_stop = delta_context_vep_local_contains(
        ctx, &view, 0, (uint8_t)'*', &ref_stop);
    alt_has_stop = delta_context_vep_local_contains(
        ctx, &view, 1, (uint8_t)'*', &alt_stop);
    ref_has_x = delta_context_vep_local_contains(
        ctx, &view, 0, (uint8_t)'X', NULL);
    alt_has_x = delta_context_vep_local_contains(
        ctx, &view, 1, (uint8_t)'X', NULL);

    if (!alt_has_x) {
        delta->stop_lost = (uint8_t)(ref_has_stop && !alt_has_stop);
    } else if (overlaps_terminal && ctx->single_edit_ref_len == 0u) {
        /* An insertion before the terminal codon does not satisfy VEP's ordinary
         * stop-overlap test. Inside that codon, an X-bearing local peptide uses
         * the original-end reconstruction for stop_lost, but stop_retained's
         * CIL predicate remains false when the local value is incomplete. */
        if ((uint64_t)ctx->single_edit_cds_start >
            (uint64_t)terminal_start) {
            status = delta_context_original_endpoint_stop_altered(
                ctx, &endpoint_altered, &endpoint_aa);
            if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
            delta->stop_lost = (uint8_t)endpoint_altered;
        }
    } else if (overlaps_terminal) {
        status = delta_context_original_endpoint_stop_altered(
            ctx, &endpoint_altered, &endpoint_aa);
        if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
        delta->stop_lost = (uint8_t)endpoint_altered;
    }

    if (!delta->stop_lost) {
        if (view.alt_length != 0u && !alt_has_x) {
            delta->stop_retained = (uint8_t)
                delta_context_vep_ref_eq_alt_sequence(ctx, &view);
        } else if (overlaps_terminal && ctx->single_edit_ref_len != 0u) {
            status = delta_context_original_endpoint_stop_altered(
                ctx, &endpoint_altered, &endpoint_aa);
            if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
            delta->stop_retained = (uint8_t)!endpoint_altered;
        }
    }
    if (!delta->stop_lost && !delta->stop_retained &&
        alt_has_stop && !ref_has_stop) {
        delta->stop_gained = 1u;
    }

    frameshift = !delta->stop_retained &&
        !(view.ref_length != 0u &&
          delta_context_vep_local_peptide_base(
              ctx, &view, 0, 0u) == (uint8_t)'*') &&
        (ctx->length_diff % 3) != 0;
    delta->frameshift = (uint8_t)frameshift;

    status = delta_context_start_facts(
        ctx, tx_flags, &start_overlaps, delta);
    if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;

    if (!frameshift && view.alt_nt_length > view.ref_nt_length &&
        !(view.ref_length == 1u && view.alt_length == 1u &&
          delta_context_vep_local_peptide_base(
              ctx, &view, 0, 0u) == (uint8_t)'*' &&
          delta_context_vep_local_peptide_base(
              ctx, &view, 1, 0u) == (uint8_t)'*')) {
        alt_shape_length = view.alt_length;
        if (alt_has_stop && alt_stop + 1u < alt_shape_length) {
            alt_shape_length = alt_stop + 1u;
        }
        preserves_ref = delta_context_vep_alt_preserves_ref_edge(
            ctx, &view, alt_shape_length);
        if (delta->start_lost ||
            (delta->start_retained &&
             delta_context_vep_alt_has_ref_suffix(
                 ctx, &view, view.alt_length))) {
            preserves_ref = 0;
        }
        if (preserves_ref) delta->inframe_insertion = 1u;
    }

    if (!frameshift && view.alt_nt_length < view.ref_nt_length &&
        !(view.ref_length == 1u &&
          delta_context_vep_local_peptide_base(
              ctx, &view, 0, 0u) == (uint8_t)'*')) {
        inframe_deletion = delta_norm_window_prefix_or_suffix(
            ctx->ref_cds + nt_offset, view.ref_nt_length,
            ctx->alt_cds + nt_offset, view.alt_nt_length);
        if (!inframe_deletion) {
            while (trim_prefix < view.ref_nt_length &&
                   trim_prefix < view.alt_nt_length &&
                   delta_norm_span_equal(
                       ctx->ref_cds + nt_offset + trim_prefix,
                       ctx->alt_cds + nt_offset + trim_prefix, 1u)) {
                trim_prefix++;
            }
            while (trim_suffix < view.ref_nt_length - trim_prefix &&
                   trim_suffix < view.alt_nt_length - trim_prefix &&
                   delta_norm_span_equal(
                       ctx->ref_cds + nt_offset + view.ref_nt_length - 1u -
                           trim_suffix,
                       ctx->alt_cds + nt_offset + view.alt_nt_length - 1u -
                           trim_suffix, 1u)) {
                trim_suffix++;
            }
            ref_remainder = view.ref_nt_length - trim_prefix - trim_suffix;
            alt_remainder = view.alt_nt_length - trim_prefix - trim_suffix;
            inframe_deletion = alt_remainder == 0u &&
                               (ref_remainder % 3u) == 0u;
        }
        if (inframe_deletion) delta->inframe_deletion = 1u;
    }

    preserves_ref = delta_context_vep_alt_preserves_ref_edge(
        ctx, &view, view.alt_length);
    if (view.ref_length != view.alt_length &&
        (view.ref_length == 0u ||
         delta_context_vep_local_peptide_base(
             ctx, &view, 0, 0u) != (uint8_t)'*') &&
        (view.alt_length == 0u ||
         delta_context_vep_local_peptide_base(
             ctx, &view, 1, 0u) != (uint8_t)'*') &&
        !preserves_ref && !inframe_deletion && !delta->start_lost &&
        !frameshift) {
        delta->protein_altering = 1u;
    }

    if (!delta->stop_lost && !delta->stop_retained && !delta->stop_gained &&
        !delta->start_lost && !delta->start_retained && !alt_has_x &&
        view.ref_length == view.alt_length) {
        if (delta_context_vep_local_equal(ctx, &view) && !ref_has_x) {
            delta->synonymous = 1u;
        } else if (!ref_has_x) {
            delta->missense = 1u;
        }
    }
    if (alt_has_x && !delta->frameshift && !delta->inframe_deletion &&
        !delta->protein_altering && !delta->start_retained &&
        !delta->start_lost && !delta->stop_retained && !delta->stop_lost) {
        delta->coding_unknown = 1u;
    }

    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    delta->protein_pos = (int32_t)(view.peptide_offset + 1u);
    /* Scalar amino-acid fields describe a one-for-one substitution only. An
     * in-frame insertion/deletion, frameshift, or variable peptide replacement
     * keeps its protein position but leaves both scalar residues invalid. */
    if (!delta->inframe_insertion && !delta->inframe_deletion &&
        !delta->frameshift && !delta->protein_altering) {
        delta->ref_aa = view.ref_length == 1u
            ? delta_context_vep_local_peptide_base(ctx, &view, 0, 0u) : 0u;
        delta->alt_aa = view.alt_length == 1u
            ? delta_context_vep_local_peptide_base(ctx, &view, 1, 0u) : 0u;
    }
    if (alt_has_x && (delta->stop_lost || delta->stop_retained) &&
        !delta->inframe_insertion && !delta->inframe_deletion &&
        !delta->frameshift && !delta->protein_altering) {
        delta->ref_aa = (uint8_t)'*';
        delta->alt_aa = endpoint_aa;
    }
    delta->valid = 1u;
    return DUCKVEP_CONTEXT_DELTA_OK;
}

typedef enum {
    DELTA_SUBSTITUTION_EDIT_WINDOW = 0,
    DELTA_SUBSTITUTION_UPLOADED_CDS_WINDOW,
    DELTA_SUBSTITUTION_UPLOADED_START_WINDOW
} delta_substitution_window_t;

/* One VEP substitution authority over an explicitly selected peptide window.
 * Ordinary edits select the changed CDS codons. A wholly-CDS uploaded feature may
 * select a wider window because VEP does not minimize substitutions by default. A
 * feature crossing from 5' UTR into CDS selects only the start predicate: VEP does
 * not add stop, synonymous, or missense terms from later coding bases in that same
 * feature. The enum makes those three states mutually exclusive. */
static duckvep_context_delta_status_t delta_context_substitution_window(
    const duckvep_coding_context_t *ctx,
    size_t                          lo_codon,
    size_t                          hi_codon,
    uint64_t                        tx_flags,
    delta_substitution_window_t     window,
    duckvep_sequence_delta_t       *delta) {

    size_t window_len;
    size_t c;
    size_t ref_stop = 0u;
    size_t alt_stop = 0u;
    int has_ref_stop = 0;
    int has_alt_stop = 0;
    int peptide_equal = 1;
    int overlaps_start;
    int uploaded_feature_window;
    int start_only;

    if (ctx == NULL || delta == NULL || ctx->ref_cds == NULL ||
        ctx->alt_cds == NULL || ctx->ref_peptide == NULL ||
        ctx->alt_peptide == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (window < DELTA_SUBSTITUTION_EDIT_WINDOW ||
        window > DELTA_SUBSTITUTION_UPLOADED_START_WINDOW) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    uploaded_feature_window = window != DELTA_SUBSTITUTION_EDIT_WINDOW;
    start_only = window == DELTA_SUBSTITUTION_UPLOADED_START_WINDOW;
    memset(delta, 0, sizeof *delta);
    if (ctx->ref_cds_len == 0u || ctx->ref_cds_len != ctx->alt_cds_len ||
        ctx->ref_peptide_len != ctx->alt_peptide_len ||
        hi_codon < lo_codon || hi_codon >= ctx->ref_peptide_len ||
        hi_codon >= (size_t)INT32_MAX ||
        hi_codon >= ctx->ref_cds_len / 3u ||
        hi_codon >= ctx->alt_cds_len / 3u) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    window_len = hi_codon - lo_codon + 1u;
    for (c = lo_codon; c <= hi_codon; c++) {
        uint8_t ref_aa = ctx->ref_peptide[c];
        uint8_t alt_aa = ctx->alt_peptide[c];
        size_t b;

        if (ref_aa == (uint8_t)'\0' || alt_aa == (uint8_t)'\0' ||
            ref_aa == (uint8_t)'X' || alt_aa == (uint8_t)'X') {
            return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
        }
        for (b = 0u; b < 3u; b++) {
            char rb = delta_norm_base((char)ctx->ref_cds[c * 3u + b]);
            char ab = delta_norm_base((char)ctx->alt_cds[c * 3u + b]);
            if (rb == '\0' || rb == 'N' || ab == '\0' || ab == 'N') {
                return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
            }
        }
        if (ref_aa != alt_aa) peptide_equal = 0;
        if (!has_ref_stop && ref_aa == (uint8_t)'*') {
            ref_stop = c - lo_codon;
            has_ref_stop = 1;
        }
        if (!has_alt_stop && alt_aa == (uint8_t)'*') {
            alt_stop = c - lo_codon;
            has_alt_stop = 1;
        }
    }

    overlaps_start = lo_codon == 0u &&
        (tx_flags & (uint64_t)DUCKVEP_TX_CDS_START_NF) == 0u;
    /* A feature may differ only in 5' UTR bases while retaining the uploaded
     * start-codon bases. VEP still evaluates start_retained_variant because the
     * unminimized feature overlaps the start codon. No other peptide window may
     * classify an unchanged CDS. */
    if (ctx->cds_changed == 0u &&
        !(uploaded_feature_window && overlaps_start)) {
        return DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
    }
    if (overlaps_start) {
        if (uploaded_feature_window) {
            char b0 = delta_norm_base((char)ctx->alt_cds[0]);
            char b1 = delta_norm_base((char)ctx->alt_cds[1]);
            char b2 = delta_norm_base((char)ctx->alt_cds[2]);

            /* Equal-length, unequal local peptides are neither a complete prefix nor
             * suffix, so VariationEffect::start_lost is true even when ATG survives. */
            if (!peptide_equal) delta->start_lost = 1u;

            /* start_retained_variant is an independent VEP predicate. */
            if (b0 == 'A' && b1 == 'T' && b2 == 'G') {
                delta->start_retained = 1u;
            }
        } else if (ctx->alt_peptide[0] != (uint8_t)'M') {
            delta->start_lost = 1u;
        }
    }

    if (!start_only) {
        if (has_ref_stop && !has_alt_stop) {
            delta->stop_lost = 1u;
        } else if (has_ref_stop && has_alt_stop && ref_stop == alt_stop) {
            delta->stop_retained = 1u;
        } else if (has_alt_stop && !has_ref_stop) {
            delta->stop_gained = 1u;
        }
        if (!delta->stop_lost && !delta->stop_retained && !delta->stop_gained) {
            if (peptide_equal) {
                if (!delta->start_retained) delta->synonymous = 1u;
            } else if (!delta->start_lost) {
                delta->missense = 1u;
            }
        }
    }

    delta->cdna_pos = -1;
    delta->cds_pos = -1;
    if (window_len == 1u) {
        delta->protein_pos = (int32_t)(lo_codon + 1u);
        delta->ref_aa = ctx->ref_peptide[lo_codon];
        delta->alt_aa = ctx->alt_peptide[lo_codon];
    } else {
        delta->protein_pos = -1;
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

    if (delta == NULL) return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    memset(delta, 0, sizeof *delta);
    if (ctx == NULL || ctx->ref_cds == NULL || ctx->alt_cds == NULL ||
        ctx->ref_peptide == NULL || ctx->alt_peptide == NULL) {
        return DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    }
    if (ctx->length_diff != 0) {
        duckvep_context_delta_status_t status;
        int handled;

        status = delta_context_length_change_predicates(
            ctx, tx_flags, &handled, delta);
        if (status != DUCKVEP_CONTEXT_DELTA_OK) return status;
        return handled ? DUCKVEP_CONTEXT_DELTA_OK
                       : DUCKVEP_CONTEXT_DELTA_UNSUPPORTED;
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
    return delta_context_substitution_window(
        ctx, lo_codon, hi_codon, tx_flags,
        DELTA_SUBSTITUTION_EDIT_WINDOW, delta);
}

/* Try VEP's uploaded-feature peptide view for an equal-length multi-base
 * substitution. A wholly-CDS feature selects every codon touched by its raw span.
 * A contiguous 5'-UTR/CDS feature selects only the start codon: VEP's mapper keeps
 * that translated coordinate while its start predicate separately edits the complete
 * cDNA feature. Once either shape is proven, this path is authoritative so a failure
 * cannot silently retry the smaller semantic edit. */
static int delta_feature_substitution_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *event,
    duckvep_sequence_delta_t         *delta) {

    duckvep_coding_context_t ctx;
    duckvep_coding_context_status_t coding_status;
    duckvep_variant_coding_context_status_t context_status;
    duckvep_context_delta_status_t delta_status;
    duckvep_coding_projection_t start_projection;
    duckvep_coding_projection_t end_projection;
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    const uint8_t *cds_seq;
    size_t cds_len;
    size_t ref_off;
    size_t alt_off;
    uint16_t feature_len;
    uint16_t coding_off = 0u;
    uint16_t coding_len = 0u;
    uint16_t i;
    uint32_t first_cds = UINT32_MAX;
    uint32_t last_cds = 0u;
    uint64_t feature_span;
    uint64_t tx_flags = 0u;
    int utr5_boundary = 0;

    if (transcripts == NULL || exons == NULL || seq == NULL || v == NULL ||
        scratch == NULL || event == NULL || delta == NULL ||
        v->ref_length == NULL || v->alt_length == NULL ||
        variant_idx >= v->count || tx_idx >= transcripts->transcript_count ||
        v->ref_length[variant_idx] != v->alt_length[variant_idx] ||
        v->ref_length[variant_idx] <= 1u || event->interbase ||
        event->feature_start1 == 0u || event->feature_end1 < event->feature_start1) {
        return 0;
    }
    feature_len = v->ref_length[variant_idx];
    feature_span = (uint64_t)event->feature_end1 -
                   (uint64_t)event->feature_start1 + 1u;
    if (feature_span != (uint64_t)feature_len) return 0;

    if (transcripts->flags != NULL) tx_flags = transcripts->flags[tx_idx];

    /* Fast common case: both endpoints are coding and the CDS span is exactly
     * the uploaded genomic span. */
    if (duckvep_project_coding_base(transcripts, exons, tx_idx,
                                    event->feature_start1, &start_projection) &&
        duckvep_project_coding_base(transcripts, exons, tx_idx,
                                    event->feature_end1, &end_projection)) {
        first_cds = start_projection.cds_pos < end_projection.cds_pos
            ? start_projection.cds_pos : end_projection.cds_pos;
        last_cds = start_projection.cds_pos > end_projection.cds_pos
            ? start_projection.cds_pos : end_projection.cds_pos;
        if (first_cds == UINT32_MAX || last_cds < first_cds ||
            (uint64_t)last_cds - (uint64_t)first_cds + 1u !=
                (uint64_t)feature_len) {
            return 0;
        }
        coding_len = feature_len;
    } else {
        uint32_t feature_start_cdna;
        uint32_t feature_end_cdna;
        uint32_t coding_start_cdna;
        uint32_t coding_start_genomic;
        uint64_t cdna_span;
        uint64_t coding_span;

        /* The start codon is not defined for an incomplete CDS. */
        if ((tx_flags & (uint64_t)DUCKVEP_TX_CDS_START_NF) != 0u ||
            transcripts->cds_start1 == NULL || transcripts->cds_end1 == NULL ||
            transcripts->strand == NULL) {
            return 0;
        }
        coding_start_genomic = strand > 0
            ? transcripts->cds_start1[tx_idx]
            : transcripts->cds_end1[tx_idx];
        if (coding_start_genomic == 0u ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, event->feature_start1,
                &feature_start_cdna, NULL) ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, event->feature_end1,
                &feature_end_cdna, NULL) ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, coding_start_genomic,
                &coding_start_cdna, NULL)) {
            return 0;
        }
        cdna_span = feature_start_cdna < feature_end_cdna
            ? (uint64_t)feature_end_cdna - (uint64_t)feature_start_cdna + 1u
            : (uint64_t)feature_start_cdna - (uint64_t)feature_end_cdna + 1u;
        if (cdna_span != (uint64_t)feature_len) return 0;

        if (strand > 0) {
            if (event->feature_start1 >= coding_start_genomic ||
                event->feature_end1 < coding_start_genomic ||
                event->feature_end1 > transcripts->cds_end1[tx_idx] ||
                feature_start_cdna >= coding_start_cdna) {
                return 0;
            }
            coding_span = (uint64_t)event->feature_end1 -
                          (uint64_t)coding_start_genomic + 1u;
            coding_off = (uint16_t)((uint64_t)coding_start_genomic -
                                    (uint64_t)event->feature_start1);
            if (!duckvep_project_coding_base(
                    transcripts, exons, tx_idx, coding_start_genomic,
                    &start_projection) ||
                !duckvep_project_coding_base(
                    transcripts, exons, tx_idx, event->feature_end1,
                    &end_projection)) {
                return 0;
            }
        } else {
            if (event->feature_start1 > coding_start_genomic ||
                event->feature_end1 <= coding_start_genomic ||
                event->feature_start1 < transcripts->cds_start1[tx_idx] ||
                feature_end_cdna >= coding_start_cdna) {
                return 0;
            }
            coding_span = (uint64_t)coding_start_genomic -
                          (uint64_t)event->feature_start1 + 1u;
            if (!duckvep_project_coding_base(
                    transcripts, exons, tx_idx, event->feature_start1,
                    &start_projection) ||
                !duckvep_project_coding_base(
                    transcripts, exons, tx_idx, coding_start_genomic,
                    &end_projection)) {
                return 0;
            }
        }
        if (coding_span == 0u || coding_span > UINT16_MAX ||
            (uint64_t)coding_off + coding_span > (uint64_t)feature_len) {
            return 0;
        }
        coding_len = (uint16_t)coding_span;
        first_cds = start_projection.cds_pos < end_projection.cds_pos
            ? start_projection.cds_pos : end_projection.cds_pos;
        last_cds = start_projection.cds_pos > end_projection.cds_pos
            ? start_projection.cds_pos : end_projection.cds_pos;
        if (first_cds != 1u || last_cds < first_cds ||
            (uint64_t)last_cds - (uint64_t)first_cds + 1u != coding_span) {
            return 0;
        }
        utr5_boundary = 1;
    }

    memset(delta, 0, sizeof *delta);
    if (!delta_allele_slice_ok(v, variant_idx) ||
        !delta_cds_slice(seq, tx_idx, &cds_seq, &cds_len) ||
        seq->codon_table == NULL) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
        return 1;
    }
    ref_off = (size_t)v->ref_offset[variant_idx];
    alt_off = (size_t)v->alt_offset[variant_idx];
    for (i = 0u; i < feature_len; i++) {
        char ref_base = delta_norm_base(
            (char)v->allele_bytes[ref_off + (size_t)i]);
        char alt_base = delta_norm_base(
            (char)v->allele_bytes[alt_off + (size_t)i]);
        if (ref_base == '\0' || ref_base == 'N' ||
            alt_base == '\0' || alt_base == 'N') {
            delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS;
            return 1;
        }
    }
    if (first_cds == 0u || (size_t)first_cds - 1u > cds_len ||
        (size_t)coding_len > cds_len - ((size_t)first_cds - 1u)) {
        delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_INVALID_PROJECTION;
        return 1;
    }
    for (i = 0u; i < coding_len; i++) {
        char cds_base = delta_norm_base(
            (char)cds_seq[(size_t)first_cds - 1u + (size_t)i]);
        if (cds_base == '\0' || cds_base == 'N') {
            delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS;
            return 1;
        }
    }

    memset(&edit, 0, sizeof edit);
    edit.cds_start = first_cds;
    edit.ref_len = (uint32_t)coding_len;
    edit.alt_len = (uint32_t)coding_len;
    edit.ref = v->allele_bytes + ref_off + (size_t)coding_off;
    edit.alt = v->allele_bytes + alt_off + (size_t)coding_off;
    edit.variant_strand = (int8_t)1;
    edit_set.edits = &edit;
    edit_set.count = 1u;
    coding_status = duckvep_coding_context_build(
        cds_seq, cds_len, &edit_set, strand,
        (duckvep_codon_table_t)seq->codon_table[tx_idx],
        scratch->alt_cds, scratch->alt_cds_cap,
        scratch->ref_peptide, scratch->ref_peptide_cap,
        scratch->alt_peptide, scratch->alt_peptide_cap, &ctx);
    context_status = delta_variant_context_from_context_status(coding_status);
    if (context_status != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        memset(delta, 0, sizeof *delta);
        delta->sequence_status = delta_sequence_status_from_context(context_status);
        return 1;
    }
    delta_status = delta_context_substitution_window(
        &ctx, utr5_boundary ? 0u : ((size_t)first_cds - 1u) / 3u,
        utr5_boundary ? 0u : ((size_t)last_cds - 1u) / 3u,
        tx_flags,
        utr5_boundary ? DELTA_SUBSTITUTION_UPLOADED_START_WINDOW
                      : DELTA_SUBSTITUTION_UPLOADED_CDS_WINDOW,
        delta);
    if (delta_status != DUCKVEP_CONTEXT_DELTA_OK || !delta->valid) {
        memset(delta, 0, sizeof *delta);
        delta->sequence_status = delta_sequence_status_from_delta(delta_status);
        return 1;
    }
    delta->sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    return 1;
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

/* Narrow same-codon MNV reference used by the pure-C property suite. Production MNVs use
 * the CodingContext interpreter. Unsupported reference shapes remain invalid. */
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

/* Narrow cross-codon MNV reference used by the pure-C property suite. It accepts only a
 * missense change across two adjacent non-terminal body codons; all other shapes remain
 * invalid. Production MNVs use the CodingContext interpreter. */
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

/* Narrow normalized-differing-region indel reference for the pure-C property suite:
 * non-terminal CDS-body frameshift
 * for INS/DEL/INDEL, codon-boundary in-frame insertion, non-boundary protein-altering
 * insertion, plus codon-aligned in-frame deletion. Raw VCF padding is stripped by
 * delta_edit_view_t. Production indels use the CodingContext interpreter; unsupported
 * reference shapes remain invalid. */
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

    if (route != NULL) *route = DUCKVEP_DELTA_ROUTE_DIRECT;
    if (delta == NULL) return;

    if (delta_feature_length_change_boundary_fill(
            transcripts, exons, seq, v, variant_idx, tx_idx, strand,
            prepared_event, delta)) {
        if (route != NULL) *route = DUCKVEP_DELTA_ROUTE_BOUNDARY_CONTEXT;
        return;
    }

    if (delta_feature_length_change_mapper_gap_fill(
            transcripts, exons, v, variant_idx, tx_idx, strand,
            prepared_event, delta)) {
        if (route != NULL) *route = DUCKVEP_DELTA_ROUTE_BOUNDARY_CONTEXT;
        return;
    }

    if (delta_feature_substitution_fill(
            transcripts, exons, seq, v, variant_idx, tx_idx, strand,
            scratch, prepared_event, delta)) {
        if (route != NULL) *route = DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT;
        return;
    }

    /* SNV (and any non-coding-editable kind, or a call with no edit-set scratch) takes the
     * light projection path: one codon edit, no full-CDS translation. */
    if (scratch == NULL || (kind != DUCKVEP_KIND_MNV && kind != DUCKVEP_KIND_DEL &&
                            kind != DUCKVEP_KIND_INS && kind != DUCKVEP_KIND_INDEL)) {
        duckvep_sequence_delta_fill_with_scratch_event(
            kind, transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
            NULL, prepared_event, delta);
        return;
    }

    /* Build one alternate CDS and one local VEP state. Failed projection,
     * ambiguous sequence, and unsupported multi-edit input remain explicit;
     * retrying a smaller shape-specific classifier would create a second
     * consequence authority. */
    if (route != NULL) {
        *route = kind == DUCKVEP_KIND_MNV ? DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT
               : kind == DUCKVEP_KIND_DEL ? DUCKVEP_DELTA_ROUTE_DEL_CONTEXT
               : kind == DUCKVEP_KIND_INS ? DUCKVEP_DELTA_ROUTE_INS_CONTEXT
                                          : DUCKVEP_DELTA_ROUTE_INDEL_CONTEXT;
    }
    duckvep_sequence_delta_fill_with_scratch_event(
        kind, transcripts, exons, seq, v, variant_idx, tx_idx, pos, strand,
        scratch, prepared_event, delta);
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
