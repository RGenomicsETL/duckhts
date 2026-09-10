#include "duckvep_property.h"

/* ===================================================================== *
 * Codon-bucket FUSION vs the coding-SNV kernel oracle.
 *
 * For a generated single-exon '+' coding transcript (random CDS sequence) and a
 * random single-base SNV at a non-boundary CDS position, annotate_tile's emitted
 * row must EQUAL the result of calling the tested coding-SNV kernel directly:
 * same SO term, same cds_pos / protein_pos / aa. The SNV is kept >=3 bases inside
 * each exon end so no splice_region bit is involved (that interaction is anchored
 * separately). This proves the fused codon path == sweep . classify . coding-SNV
 * composition across all codon-change classes (syn/missense/stop_gained/lost) at
 * scale, not just the deterministic anchor.
 * ===================================================================== */

static void kprop_coding_free(void *instance, void *env) {
    struct kprop_coding *s = (struct kprop_coding *)instance;
    (void)env;
    if (s == NULL) return;
    free(s->cds);
    free(s);
}

static enum theft_alloc_res kprop_coding_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u; /* 3..8 */
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t ci, ri, ai, i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    /* SNV at a non-boundary CDS index: [3 .. cds_len-4] (cds_len>=9). */
    ci = (uint32_t)kprop_bounded(t, cds_len - 6u) + 3u;
    ri = 0u; while (BASES[ri] != (char)s->cds[ci]) ri++;          /* ref base index */
    ai = (ri + 1u + (uint32_t)kprop_bounded(t, 3u)) % 4u;         /* a different base */

    s->chrom = 0u; s->strand = 1;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;

    s->vchrom = 0u; s->vpos = base + ci; s->vend = s->vpos; s->vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s->abytes[0] = s->cds[ci]; s->abytes[1] = (uint8_t)BASES[ai];
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 1u;

    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;

    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_coding_info = {
    .alloc = kprop_coding_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_start_codon_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t ci, ri, ai, i, tries;
    char first[4];
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }

    first[3] = '\0';
    for (tries = 0u; tries < 32u; tries++) {
        first[0] = BASES[kprop_bounded(t, 4u)];
        first[1] = BASES[kprop_bounded(t, 4u)];
        first[2] = BASES[kprop_bounded(t, 4u)];
        if (duckvep_translate_codon(first, DUCKVEP_CODON_TABLE_STANDARD) != '*') break;
    }
    if (tries == 32u) { first[0] = 'A'; first[1] = 'T'; first[2] = 'G'; }
    s->cds[0] = (uint8_t)first[0]; s->cds[1] = (uint8_t)first[1]; s->cds[2] = (uint8_t)first[2];
    for (i = 3u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    ci = (uint32_t)kprop_bounded(t, 3u); /* a base inside the annotated start codon */
    ri = 0u; while (BASES[ri] != (char)s->cds[ci]) ri++;
    ai = (ri + 1u + (uint32_t)kprop_bounded(t, 3u)) % 4u;

    s->chrom = 0u; s->strand = 1; s->flags = 0u;
    s->tstart = base - 20u; s->tend = base + cds_len + 19u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = s->tstart; s->ee = s->tend;
    s->ecds = 1u; s->ecde = cds_len + 40u; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;

    s->vchrom = 0u; s->vpos = base + ci; s->vend = s->vpos; s->vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s->abytes[0] = s->cds[ci]; s->abytes[1] = (uint8_t)BASES[ai];
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 1u;

    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_start_codon_info = {
    .alloc = kprop_start_codon_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_mnv_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon = (uint32_t)kprop_bounded(t, ncodons);
    uint32_t mnv_len = (uint32_t)kprop_bounded(t, 2u) + 2u;
    uint32_t codon_off = mnv_len == 3u ? 0u : (uint32_t)kprop_bounded(t, 2u);
    uint32_t ci = codon * 3u + codon_off;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = 1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;

    s->vchrom = 0u; s->vpos = base + ci; s->vend = s->vpos + mnv_len - 1u;
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < mnv_len; i++) s->abytes[i] = s->cds[ci + i];
    for (i = 0u; i < mnv_len; i++) {
        uint32_t ri = 0u;
        uint32_t ai;
        while (BASES[ri] != (char)s->abytes[i]) ri++;
        ai = (ri + 1u + (uint32_t)kprop_bounded(t, 3u)) % 4u;
        s->abytes[mnv_len + i] = (uint8_t)BASES[ai];
    }
    s->roff = 0u; s->aoff = mnv_len; s->rlen = (uint16_t)mnv_len; s->alen = (uint16_t)mnv_len;

    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = (size_t)mnv_len * 2u;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_mnv_info = {
    .alloc = kprop_mnv_alloc,
    .free  = kprop_coding_free,
};

char kprop_complement_base(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

static uint32_t kprop_cds_pos_for_genomic(const struct kprop_coding *s, uint32_t pos1) {
    uint32_t off0 = pos1 - s->es;
    return s->strand > 0 ? off0 + 1u : s->cds_lenv - off0;
}

static uint32_t kprop_genomic_pos_for_cds(const struct kprop_coding *s, uint32_t cds_pos) {
    uint32_t off0 = s->strand > 0 ? cds_pos - 1u : s->cds_lenv - cds_pos;
    return s->es + off0;
}

static char kprop_genomic_base_at(const struct kprop_coding *s, uint32_t pos1) {
    uint32_t cds_pos = kprop_cds_pos_for_genomic(s, pos1);
    char b = (char)s->cds[cds_pos - 1u];
    return s->strand > 0 ? b : kprop_complement_base(b);
}

static uint8_t kprop_base_not(uint8_t a, uint8_t b) {
    static const uint8_t BASES[4] = {'A', 'C', 'G', 'T'};
    uint32_t i;
    for (i = 0u; i < 4u; i++) {
        if (BASES[i] != a && BASES[i] != b) return BASES[i];
    }
    return (uint8_t)'N';
}

static uint8_t kprop_base_not3(uint8_t a, uint8_t b, uint8_t c) {
    static const uint8_t BASES[4] = {'A', 'C', 'G', 'T'};
    uint32_t i;
    for (i = 0u; i < 4u; i++) {
        if (BASES[i] != a && BASES[i] != b && BASES[i] != c) return BASES[i];
    }
    return (uint8_t)'N';
}

enum {
    KPROP_CDS_EDIT_SNV = 0u,
    KPROP_CDS_EDIT_MNV = 1u,
    KPROP_CDS_EDIT_INS = 2u,
    KPROP_CDS_EDIT_DEL = 3u,
    KPROP_CDS_EDIT_INDEL = 4u
};

enum {
    KPROP_CONTEXT_DELTA_SYNONYMOUS = 0u,
    KPROP_CONTEXT_DELTA_MISSENSE = 1u,
    KPROP_CONTEXT_DELTA_STOP_GAINED = 2u,
    KPROP_CONTEXT_DELTA_STOP_LOST = 3u,
    KPROP_CONTEXT_DELTA_STOP_RETAINED = 4u
};

static uint32_t kprop_pick_cds_start(struct theft *t, uint32_t cds_len,
                                     uint32_t ref_len, uint32_t region) {
    uint32_t max_start = cds_len - ref_len + 1u;
    if (region == KPROP_CDS_EDIT_START) return 1u;
    if (region == KPROP_CDS_EDIT_STOP) return max_start;
    return (uint32_t)kprop_bounded(t, max_start - 6u) + 4u;
}

static void kprop_fill_expected_cds(struct kprop_coding *s, uint32_t cds_start,
                                    uint32_t ref_len, const uint8_t *alt_tx,
                                    uint32_t alt_len) {
    uint32_t prefix = cds_start - 1u;
    uint32_t suffix_start = prefix + ref_len;
    uint32_t suffix_len = s->cds_lenv - suffix_start;
    memcpy(s->expect_cds, s->cds, prefix);
    if (alt_len > 0u) memcpy(s->expect_cds + prefix, alt_tx, alt_len);
    memcpy(s->expect_cds + prefix + alt_len, s->cds + suffix_start, suffix_len);
    s->expect_len = prefix + alt_len + suffix_len;
}

static void kprop_fill_variant_alt_from_tx(struct kprop_coding *s, uint32_t off,
                                           const uint8_t *alt_tx, uint32_t alt_len) {
    uint32_t i;
    for (i = 0u; i < alt_len; i++) {
        uint8_t b = s->strand > 0 ? alt_tx[i]
            : (uint8_t)kprop_complement_base((char)alt_tx[alt_len - 1u - i]);
        s->abytes[off + i] = b;
    }
}

static void kprop_wire_coding_scene(struct kprop_coding *s, uint32_t cds_len) {
    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;
}

TEST nmd_early_cds_fact_matches_exhaustive_projection_known_scene(void) {
    static uint8_t cds[120];
    static const uint32_t deleted_cds[2] = {101u, 102u};
    static const uint8_t expected_fact[2] = {
        (uint8_t)DUCKVEP_NMD_EARLY_CDS_ENDS_THROUGH_101,
        (uint8_t)DUCKVEP_NMD_EARLY_CDS_ENDS_AFTER_101
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[160];
    uint8_t ref_peptide[64];
    uint8_t alt_peptide[64];
    duckvep_delta_scratch_t scratch;
    size_t i;

    memset(&s, 0, sizeof s);
    memset(cds, 'A', sizeof cds);
    cds[0] = (uint8_t)'A'; cds[1] = (uint8_t)'T'; cds[2] = (uint8_t)'G';
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1119u; s.cds_s = 1000u; s.cds_e = 1119u;
    s.es = 1000u; s.ee = 1119u; s.ecds = 1u; s.ecde = 120u;
    s.eph = 0; s.eeph = 0; s.exoff = 0u; s.excnt = 1u;
    s.vchrom = 0u; s.vkind = (uint8_t)DUCKVEP_KIND_DEL;
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 1u;
    s.abytes[0] = (uint8_t)'A'; s.abytes[1] = (uint8_t)'A';
    s.abytes[2] = (uint8_t)'A';
    kprop_wire_coding_scene(&s, sizeof cds);

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_peptide;
    scratch.ref_peptide_cap = sizeof ref_peptide;
    scratch.alt_peptide = alt_peptide;
    scratch.alt_peptide_cap = sizeof alt_peptide;

    for (i = 0u; i < 2u; i++) {
        duckvep_event_t event;
        duckvep_sequence_delta_t delta;
        duckvep_sequence_delta_route_t route;
        duckvep_nmd_result_t cached;
        duckvep_nmd_result_t exhaustive;
        uint64_t frameshift = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);

        s.vpos = s.tstart + deleted_cds[i] - 2u;
        s.vend = s.vpos + 1u;
        duckvep_event_load(&s.v, 0u, &event);
        duckvep_sequence_delta_fill_for_annotation_trace(
            DUCKVEP_KIND_DEL, &s.tx, &s.ex, &s.seq, &s.v,
            0u, 0u, s.vpos, s.strand, &scratch, &event,
            (uint32_t)DUCKVEP_REGION_CDS, 0u, &route, &delta);
        ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SIMPLE_INDEL, route);
        ASSERT(delta.valid && delta.frameshift);
        ASSERT_EQ(expected_fact[i], delta.nmd_early_cds_fact);

        duckvep_nmd_predict(&s.tx, &s.ex, 0u, &event, frameshift,
                            &delta, &cached);
        duckvep_nmd_predict(&s.tx, &s.ex, 0u, &event, frameshift,
                            NULL, &exhaustive);
        ASSERT_EQ(exhaustive.prediction, cached.prediction);
        ASSERT_EQ(exhaustive.escape_reasons, cached.escape_reasons);
    }
    PASS();
}

TEST nmd_later_phase_cds_does_not_reuse_physical_early_fact(void) {
    uint8_t cds[122];
    struct kprop_coding s;
    uint32_t starts[] = {90u, 1000u}, ends[] = {99u, 1119u};
    uint32_t cdna_starts[] = {1u, 11u}, cdna_ends[] = {10u, 130u};
    int8_t phases[] = {-1, 2}, end_phases[] = {-1, 2};
    duckvep_haplotype_edit_t edit;
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t delta;
    duckvep_sequence_delta_route_t route;
    duckvep_event_t event;
    duckvep_nmd_result_t cached, exhaustive;
    int8_t strand, phase;

    for (strand = -1; strand <= 1; strand += 2) {
        for (phase = 0; phase <= 2; phase++) {
            memset(&s, 0, sizeof s);
            memset(cds, 'A', sizeof cds); memset(cds, 'N', (size_t)phase);
            phases[1] = phase; end_phases[1] = phase;
            starts[0] = strand > 0 ? 90u : 2000u; ends[0] = strand > 0 ? 99u : 2009u;
            starts[1] = strand > 0 ? 1000u : 980u; ends[1] = strand > 0 ? 1119u : 1099u;
            s.cds = cds; s.strand = strand;
            s.tstart = strand > 0 ? 90u : 980u; s.tend = strand > 0 ? 1119u : 2009u;
            s.cds_s = starts[1]; s.cds_e = ends[1];
            s.excnt = 2u; s.vpos = strand > 0 ? 1099u : 998u; s.vend = s.vpos + 1u;
            s.vkind = (uint8_t)DUCKVEP_KIND_DEL;
            s.aoff = 2u; s.rlen = 2u; s.alen = 1u;
            memset(s.abytes, strand > 0 ? 'A' : 'T', 3u);
            kprop_wire_coding_scene(&s, 120u + (uint32_t)phase);
            s.ex.start1 = starts; s.ex.end1 = ends;
            s.ex.cdna_start1 = cdna_starts; s.ex.cdna_end1 = cdna_ends;
            s.ex.phase = phases; s.ex.end_phase = end_phases; s.ex.exon_count = 2u;
            memset(&scratch, 0, sizeof scratch);
            scratch.edits = &edit; scratch.edits_cap = 1u;
            duckvep_event_load(&s.v, 0u, &event);
            duckvep_sequence_delta_fill_for_annotation_trace(
                DUCKVEP_KIND_DEL, &s.tx, &s.ex, &s.seq, &s.v,
                0u, 0u, s.vpos, s.strand, &scratch, &event,
                (uint32_t)DUCKVEP_REGION_CDS, 1u, &route, &delta);
            ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SIMPLE_INDEL, route);
            ASSERT(delta.valid && delta.frameshift);
            /* Physical CDS 101+phase is VEP feature CDS 101. The producer now
             * caches the feature-coordinate fact, like the exhaustive mapper. */
            ASSERT_EQ(DUCKVEP_NMD_EARLY_CDS_ENDS_THROUGH_101, delta.nmd_early_cds_fact);
            duckvep_nmd_predict(&s.tx, &s.ex, 0u, &event,
                DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), NULL, &exhaustive);
            ASSERT(exhaustive.escape_reasons & DUCKVEP_NMD_ESCAPE_EARLY_CDS);
            duckvep_nmd_predict(&s.tx, &s.ex, 0u, &event,
                DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), &delta, &cached);
            ASSERT_EQ(exhaustive.prediction, cached.prediction);
            ASSERT_EQ(exhaustive.escape_reasons, cached.escape_reasons);
            /* Retain the negative control: a stale physical-coordinate cache
             * must not override exhaustive NMD projection in a rephased model. */
            if (phase > 0) {
                delta.nmd_early_cds_fact = DUCKVEP_NMD_EARLY_CDS_ENDS_AFTER_101;
                duckvep_nmd_predict(&s.tx, &s.ex, 0u, &event,
                    DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), &delta, &cached);
                ASSERT_EQ(exhaustive.prediction, cached.prediction);
                ASSERT_EQ(exhaustive.escape_reasons, cached.escape_reasons);
            }
        }
    }
    PASS();
}

/* VEP's coding_unknown predicate wins when either local peptide allele contains
 * X. Incomplete-start models carry N padding in the first codon; equal X/X
 * peptide bytes are not evidence for a synonymous consequence. */
TEST sequence_delta_snv_x_peptide_is_coding_unknown_known_scene(void) {
    static uint8_t cds[6] = {'N','C','A', 'T','A','A'};
    struct kprop_coding s;
    duckvep_sequence_delta_t delta;

    memset(&s, 0, sizeof s);
    s.cds = cds;
    s.chrom = 0u;
    s.strand = 1;
    s.flags = (uint64_t)DUCKVEP_TX_CDS_START_NF;
    s.tstart = 1000u; s.tend = 1005u;
    s.cds_s = 1000u; s.cds_e = 1005u;
    s.es = 1000u; s.ee = 1005u;
    s.ecds = 1u; s.ecde = 6u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u;
    s.vchrom = 0u; s.vpos = 1001u; s.vend = 1001u;
    s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s.abytes[0] = (uint8_t)'C'; s.abytes[1] = (uint8_t)'T';
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 1u;
    kprop_wire_coding_scene(&s, sizeof cds);

    duckvep_sequence_delta_fill_for_annotation(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v,
        0u, 0u, s.vpos, s.strand, NULL, NULL, &delta);
    ASSERT(delta.valid && delta.coding_unknown && !delta.synonymous &&
           !delta.missense && !delta.stop_gained && !delta.stop_lost);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);
    ASSERT_EQ((uint8_t)'X', delta.ref_aa);
    ASSERT_EQ((uint8_t)'X', delta.alt_aa);
    PASS();
}

TEST sequence_delta_terminal_substitution_borrows_utr(void) {
    /* Reduced from the four pinned chrDuck terminal substitution witnesses.
     * REF codons end inside CDS; ALT codons can borrow one known UTR base. */
    static uint8_t cds[2][12] = {"ATGCCCTGGTA", "ATGCCCGGTAA"};
    static const uint32_t starts[4] = {9u, 11u, 9u, 10u};
    static const char *refs[4] = {"GT", "A", "TA", "A"};
    static const char *alts[4] = {"TT", "G", "AT", "T"};
    static const char *ref_peptides[4] = {"WX", "X", "GX", "X"};
    static const char *alt_peptides[4] = {"C*", "*", "G*", "*"};
    static const uint8_t post[] = "A";
    uint64_t offset = 0u;
    uint32_t pre_length = 0u, post_length = 1u;
    int8_t strand;
    size_t i;
    for (strand = -1; strand <= 1; strand += 2) for (i = 0u; i < 4u; i++) {
        struct kprop_coding s;
        duckvep_sequence_delta_t delta;
        duckvep_delta_scratch_t scratch;
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edits;
        duckvep_coding_context_t context;
        duckvep_coding_peptide_window_t view;
        uint8_t alt_cds[16], ref_peptide[8], alt_peptide[8];
        size_t b, length = strlen(refs[i]);
        memset(&s, 0, sizeof s);
        memset(&scratch, 0, sizeof scratch);
        scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
        scratch.ref_peptide = ref_peptide; scratch.ref_peptide_cap = sizeof ref_peptide;
        scratch.alt_peptide = alt_peptide; scratch.alt_peptide_cap = sizeof alt_peptide;
        s.cds = cds[i / 2u]; s.strand = strand;
        s.tstart = s.es = 100u; s.tend = s.ee = 111u;
        s.cds_s = strand > 0 ? 100u : 101u;
        s.cds_e = strand > 0 ? 110u : 111u;
        s.ecds = 1u; s.ecde = 12u; s.excnt = 1u;
        s.vpos = strand > 0 ? 99u + starts[i] : 113u - starts[i] - (uint32_t)length;
        s.vend = s.vpos + (uint32_t)length - 1u;
        s.vkind = length == 1u ? DUCKVEP_KIND_SNV : DUCKVEP_KIND_MNV;
        s.rlen = s.alen = (uint16_t)length; s.aoff = (uint32_t)length;
        for (b = 0u; b < length; b++) {
            size_t source = strand > 0 ? b : length - 1u - b;
            s.abytes[b] = (uint8_t)coding_test_genomic_from_tx(refs[i][source], strand);
            s.abytes[length + b] = (uint8_t)coding_test_genomic_from_tx(alts[i][source], strand);
        }
        kprop_wire_coding_scene(&s, 11u);
        s.seq.flank_bytes = post; s.seq.flank_bytes_len = 1u;
        s.seq.pre_cds_offset = s.seq.post_cds_offset = &offset;
        s.seq.pre_cds_length = &pre_length; s.seq.post_cds_length = &post_length;
        s.seq.flanks_complete = 1u;
        duckvep_sequence_delta_fill_for_annotation((duckvep_variant_kind_t)s.vkind,
            &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, &scratch, NULL, &delta);
        ASSERT(delta.valid && !delta.start_lost && !delta.start_retained &&
            !delta.stop_lost && !delta.missense && !delta.synonymous &&
            !delta.frameshift && !delta.inframe_insertion && !delta.inframe_deletion &&
            !delta.protein_altering);
        ASSERT_EQ(i != 2u, delta.stop_gained);
        ASSERT_EQ(i == 2u, delta.stop_retained);
        ASSERT_EQ(i == 1u || i == 3u, delta.partial_codon);
        ASSERT_EQ(i != 2u, delta.coding_unknown);
        ASSERT_EQ(DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);
        /* The materialized mechanics context uses the identical codon view. */
        memset(&edit, 0, sizeof edit);
        edit.cds_start = starts[i]; edit.ref_len = edit.alt_len = (uint32_t)length;
        edit.ref = (const uint8_t *)refs[i]; edit.alt = (const uint8_t *)alts[i];
        edit.variant_strand = strand; edits.edits = &edit; edits.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
            s.cds, 11u, &edits, strand, DUCKVEP_CODON_TABLE_STANDARD,
            alt_cds, sizeof alt_cds, ref_peptide, sizeof ref_peptide,
            alt_peptide, sizeof alt_peptide, &context));
        context.post_cds_bases = post; context.post_cds_length = 1u;
        context.post_cds_complete = 1u;
        ASSERT(duckvep_coding_context_peptide_window_open(&context, &view));
        ASSERT_EQ(strlen(ref_peptides[i]), view.ref_length);
        ASSERT_EQ(strlen(alt_peptides[i]), view.alt_length);
        for (b = 0u; b < view.ref_length; b++) ASSERT_EQ((uint8_t)ref_peptides[i][b],
            duckvep_coding_context_peptide_window_base(&context, &view, 0, b));
        for (b = 0u; b < view.alt_length; b++) ASSERT_EQ((uint8_t)alt_peptides[i][b],
            duckvep_coding_context_peptide_window_base(&context, &view, 1, b));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_delta_fill(&context, 0u, &delta));
        /* Physical edit-set mechanics minimize GT>TT to G>T, whereas the
         * independent VEP feature keeps the retained T and its second codon. */
        ASSERT_EQ(i == 0u, delta.missense);
        ASSERT_EQ(i == 1u || i == 3u, delta.stop_gained);
        ASSERT_EQ(i == 2u, delta.stop_retained);
        context.post_cds_length = 0u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_delta_fill(&context, 0u, &delta));
        ASSERT(!delta.stop_gained && !delta.stop_retained);
        {
            const char *wrong = "ACGT";
            while (*wrong == s.abytes[0] || *wrong == s.abytes[length]) wrong++;
            s.abytes[0] = (uint8_t)*wrong;
        }
        duckvep_sequence_delta_fill_for_annotation((duckvep_variant_kind_t)s.vkind,
            &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, &scratch, NULL, &delta);
        ASSERT(!delta.valid);
        ASSERT_EQ(DUCKVEP_SEQUENCE_REFERENCE_MISMATCH, delta.sequence_status);
    }
    PASS();
}

TEST sequence_delta_later_cds_phase_uses_feature_edits(void) {
    /* Pinned executable witness: projection_fixtures.R's later CDS start,
     * chrDuck:158 C>A. VEP's UTR+CDS test precedes its X-peptide rejection.
     * Mirror the same transcript-oriented model for the reverse-strand control. */
    static const uint8_t unpadded[] = "CGTACGTACGTACGTACGTACGTTACGTACGTACGTACTGGTAA";
    static const uint32_t cdna_start[] = {1u, 27u, 58u};
    static const uint32_t cdna_end[] = {26u, 57u, 88u};
    /* Same retained VEP campaign: body, stop, and final physical CDS base.
     * 0=missense, 1=synonymous, 2=stop gained, 3=stop lost, 4=partial codon. */
    static const uint32_t positions[] = {162u, 162u, 238u, 240u};
    static const char references[] = "CCTA";
    static const char alternates[] = "ATAC";
    static const uint8_t consequences[3][4] = {{0,0,1,4}, {1,0,2,3}, {0,1,1,0}};
    uint8_t cds[sizeof unpadded + 2u];
    uint8_t flanks[44];
    uint64_t pre_offset = 0u, post_offset = 34u;
    uint32_t pre_length = 34u, post_length = 10u;
    int8_t strand;
    int8_t phase;

    memset(flanks, 'A', sizeof flanks);
    memcpy(flanks + 34u, "ACGTACGTAC", 10u);
    for (strand = -1; strand <= 1; strand += 2) {
        for (phase = 0; phase <= 2; phase++) {
            struct kprop_coding s;
            duckvep_sequence_delta_t delta;
            duckvep_coding_projection_t projection;
            uint32_t starts[] = {100u, 150u, 220u};
            uint32_t ends[] = {125u, 180u, 250u};
            int8_t phases[] = {-1, phase, (int8_t)((phase + 23) % 3)};
            int8_t end_phases[] = {-1, phases[2], phases[2]};
            size_t i;

            memset(&s, 0, sizeof s);
            memset(cds, 'N', (size_t)phase);
            memcpy(cds + phase, unpadded, sizeof unpadded - 1u);
            s.cds = cds; s.strand = strand;
            s.tstart = 100u; s.tend = 250u;
            s.cds_s = strand > 0 ? 158u : 110u;
            s.cds_e = strand > 0 ? 240u : 192u;
            s.excnt = 3u;
            s.vpos = strand > 0 ? 158u : 192u; s.vend = s.vpos;
            s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
            s.abytes[0] = strand > 0 ? 'C' : 'G';
            s.abytes[1] = strand > 0 ? 'A' : 'T';
            s.aoff = 1u; s.rlen = 1u; s.alen = 1u;
            kprop_wire_coding_scene(&s, (uint32_t)(sizeof unpadded - 1u + phase));
            if (strand < 0) for (i = 0u; i < 3u; i++) {
                uint32_t start = starts[i];
                starts[i] = 350u - ends[i]; ends[i] = 350u - start;
            }
            s.ex.start1 = starts; s.ex.end1 = ends;
            s.ex.cdna_start1 = cdna_start; s.ex.cdna_end1 = cdna_end;
            s.ex.phase = phases; s.ex.end_phase = end_phases; s.ex.exon_count = 3u;
            s.seq.flank_bytes = flanks; s.seq.flank_bytes_len = sizeof flanks;
            s.seq.pre_cds_offset = &pre_offset; s.seq.pre_cds_length = &pre_length;
            s.seq.post_cds_offset = &post_offset; s.seq.post_cds_length = &post_length;
            s.seq.flanks_complete = 1u;

            ASSERT(duckvep_project_coding_base(&s.tx, &s.ex, 0u, s.vpos, &projection));
            ASSERT_EQ((uint32_t)phase + 1u, projection.cds_pos);
            ASSERT_EQ('C', cds[projection.cds_pos - 1u]);
            duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_SNV,
                &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
            ASSERT(delta.valid && delta.start_lost && !delta.coding_unknown && !delta.missense);
            ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);

            s.flags = (uint64_t)DUCKVEP_TX_CDS_START_NF;
            duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_SNV,
                &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
            ASSERT(delta.valid && !delta.start_lost && !delta.start_retained);
            s.flags = 0u;
            s.abytes[0] = strand > 0 ? 'G' : 'C';
            duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_SNV,
                &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
            ASSERT(!delta.valid && !delta.start_lost);
            ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH, delta.sequence_status);

            for (i = 0u; i < 4u; i++) {
                uint8_t expected = consequences[phase][i];
                const char *wrong_ref = "ACGT";
                s.vpos = strand > 0 ? positions[i] : 350u - positions[i];
                s.vend = s.vpos;
                s.abytes[0] = (uint8_t)coding_test_genomic_from_tx(references[i], strand);
                s.abytes[1] = (uint8_t)coding_test_genomic_from_tx(alternates[i], strand);
                {
                    duckvep_coding_projection_t display;
                    duckvep_event_t event;
                    uint32_t expected_cds = positions[i] <= 180u
                        ? positions[i] - 157u : positions[i] - 196u;
                    uint32_t first, last;
                    ASSERT(duckvep_project_coding_base(&s.tx, &s.ex, 0u, s.vpos, &projection));
                    ASSERT_EQ(expected_cds + (uint32_t)phase, projection.cds_pos);
                    ASSERT(duckvep_project_vep_coding_position(&s.tx, &s.ex, 0u,
                        &projection, &display));
                    ASSERT_EQ(expected_cds, display.cds_pos);
                    ASSERT_EQ((expected_cds - 1u) / 3u + 1u, display.protein_pos);
                    ASSERT_EQ((expected_cds - 1u) % 3u, display.codon_offset);
                    ASSERT_EQ(0u, display.phase_offset);
                    ASSERT_EQ(projection.cdna_pos, display.cdna_pos);
                    duckvep_event_load(&s.v, 0u, &event);
                    ASSERT(duckvep_project_feature_translation_start(&s.tx, &s.ex, 0u,
                        &event, &display));
                    ASSERT_EQ(expected_cds, display.cds_pos);
                    ASSERT(duckvep_project_complete_feature_translation_bounds(
                        &s.tx, &s.ex, 0u, &event, 0u, &first, &last));
                    ASSERT_EQ(expected_cds, first); ASSERT_EQ(expected_cds, last);
                    ASSERT(duckvep_project_feature_to_cds(&s.tx, &s.ex, 0u, &event, &first, &last));
                    ASSERT_EQ(expected_cds, first); ASSERT_EQ(expected_cds, last);
                    ASSERT(duckvep_project_event_to_cds(&s.tx, &s.ex, 0u, &event, &first, &last));
                    ASSERT_EQ(projection.cds_pos, first); ASSERT_EQ(projection.cds_pos, last);
                }
                duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_SNV,
                    &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
                ASSERT(delta.valid && !delta.start_lost && !delta.start_retained &&
                       !delta.stop_retained && !delta.frameshift);
                ASSERT_EQ(expected == 0u, delta.missense);
                ASSERT_EQ(expected == 1u, delta.synonymous);
                ASSERT_EQ(expected == 2u, delta.stop_gained);
                ASSERT_EQ(expected == 3u, delta.stop_lost);
                ASSERT_EQ(expected == 4u, delta.partial_codon);
                ASSERT_EQ(expected == 4u, delta.coding_unknown);
                ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);
                /* The displayed codon's original byte can differ from genomic
                 * REF. It must never become the physical reference oracle. */
                while (*wrong_ref == references[i] || *wrong_ref == alternates[i]) wrong_ref++;
                s.abytes[0] = (uint8_t)coding_test_genomic_from_tx(*wrong_ref, strand);
                duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_SNV,
                    &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
                ASSERT(!delta.valid);
                ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH, delta.sequence_status);
            }
            {
                duckvep_event_t event;
                duckvep_coding_projection_t display;
                uint32_t first, last;
                uint32_t expected_insert = strand > 0 ? 6u : 5u;
                uint8_t feature_phase = 99u;

                s.vpos = strand > 0 ? 162u : 188u; s.vend = s.vpos;
                s.vkind = (uint8_t)DUCKVEP_KIND_INS; s.alen = 2u;
                s.abytes[0] = strand > 0 ? 'C' : 'G';
                s.abytes[1] = s.abytes[0]; s.abytes[2] = 'A';
                duckvep_event_load(&s.v, 0u, &event);
                ASSERT(duckvep_project_event_to_cds(&s.tx, &s.ex, 0u, &event, &first, &last));
                ASSERT_EQ(expected_insert + (uint32_t)phase, first);
                ASSERT_EQ(first, last);
                ASSERT(duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    first, (uint8_t)phase, &first, &feature_phase));
                ASSERT_EQ(expected_insert, first);
                ASSERT_EQ(0u, feature_phase);
                ASSERT(duckvep_project_feature_to_cds(&s.tx, &s.ex, 0u, &event, &first, &last));
                ASSERT_EQ(expected_insert, first); ASSERT_EQ(first - 1u, last);
                ASSERT(duckvep_project_complete_feature_translation_bounds(
                    &s.tx, &s.ex, 0u, &event, 0u, &first, &last));
                ASSERT_EQ(expected_insert, first); ASSERT_EQ(first - 1u, last);

                /* Rephasing rejects malformed phases and checked-add overflow;
                 * failures zero the output even for an in-place conversion. */
                ASSERT(duckvep_project_coding_base(&s.tx, &s.ex, 0u, s.vpos, &projection));
                phases[0] = 3;
                ASSERT(!duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    10u, 0u, &first, &feature_phase));
                ASSERT_EQ(0u, first); ASSERT_EQ(0u, feature_phase);
                ASSERT(!duckvep_project_vep_coding_position(&s.tx, &s.ex, 0u, &projection, &display));
                ASSERT_EQ(0u, display.cds_pos);
                phases[0] = -2;
                ASSERT(!duckvep_project_vep_coding_position(&s.tx, &s.ex, 0u, &projection, &display));
                phases[0] = 2;
                ASSERT(duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    (uint32_t)phase + 5u, (uint8_t)phase, &first, &feature_phase));
                ASSERT_EQ(7u, first); ASSERT_EQ(2u, feature_phase);
                ASSERT(!duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    UINT32_MAX, 0u, &first, &feature_phase));
                ASSERT_EQ(0u, first); ASSERT_EQ(0u, feature_phase);
                projection.cds_pos = UINT32_MAX; projection.phase_offset = 0u;
                ASSERT(!duckvep_project_vep_coding_position(&s.tx, &s.ex, 0u, &projection, &projection));
                ASSERT_EQ(0u, projection.cds_pos);
                phases[0] = -1;
                ASSERT(!duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    2u, 2u, &first, &feature_phase));
                ASSERT(!duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    10u, 3u, &first, &feature_phase));
                first = 99u;
                ASSERT(!duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    10u, 0u, &first, NULL));
                ASSERT_EQ(0u, first);
                feature_phase = 99u;
                ASSERT(!duckvep_project_vep_cds_position(&s.tx, &s.ex, 0u,
                    10u, 0u, NULL, &feature_phase));
                ASSERT_EQ(0u, feature_phase);
                ASSERT(!duckvep_project_vep_coding_position(&s.tx, &s.ex, 0u, &projection, &display));
            }
            {
                /* Same pinned physical campaign, chrDuck:157 ACG>ATT and
                 * ACGT>AATG: edit the full UTR+CDS string at cDNA 34, not a
                 * phase-adjusted/clipped CDS replacement. */
                static const char *refs[] = {"ACG", "ACGT"};
                static const char *alts[] = {"ATT", "AATG"};
                size_t j;
                s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
                for (j = 0u; j < 2u; j++) {
                    size_t n = strlen(refs[j]);
                    size_t b;
                    s.rlen = (uint16_t)n; s.alen = (uint16_t)n;
                    s.aoff = (uint32_t)n;
                    s.vpos = strand > 0 ? 157u : 350u - (156u + (uint32_t)n);
                    s.vend = s.vpos + (uint32_t)n - 1u;
                    for (b = 0u; b < n; b++) {
                        size_t oriented = strand > 0 ? b : n - 1u - b;
                        s.abytes[b] = (uint8_t)coding_test_genomic_from_tx(refs[j][oriented], strand);
                        s.abytes[n + b] = (uint8_t)coding_test_genomic_from_tx(alts[j][oriented], strand);
                    }
                    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
                    ASSERT(delta.valid && !delta.coding_unknown && !delta.missense &&
                           !delta.synonymous && !delta.stop_gained && !delta.stop_lost &&
                           !delta.stop_retained && !delta.frameshift && !delta.protein_altering);
                    ASSERT_EQ(j == 0u, delta.start_lost);
                    ASSERT_EQ(j == 1u, delta.start_retained);
                    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);
                    /* Retained UTR REF must also match the physical transcript.
                     * The smaller semantic CDS edit is unchanged by this mutation. */
                    b = strand > 0 ? 0u : n - 1u;
                    s.abytes[b] = s.abytes[n + b] =
                        (uint8_t)coding_test_genomic_from_tx('C', strand);
                    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, NULL, NULL, &delta);
                    if (delta.valid) fprintf(stderr,
                        "[MNV retained UTR REF] strand=%d phase=%d case=%zu status=%u\n",
                        strand, phase, j, delta.sequence_status);
                    ASSERT(!delta.valid && !delta.start_lost && !delta.start_retained);
                    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH, delta.sequence_status);
                }
            }
            {
                /* Pinned VEP indels: independent start, codon shape, terminal
                 * stop, and CDS-to-UTR edits. The physical edit must survive
                 * opening the differently phased consequence view unchanged. */
                static const uint32_t positions[] = {158u,160u,165u,165u,235u,238u,239u,239u};
                static const char *refs[] = {"C","T","A","ACGT","TGGT","T","AA","AAAC"};
                static const char *alts[] = {"CT","AC","AATG","A","T","TAGGT","A","A"};
                static const uint32_t unpadded_starts[] = {2u,3u,9u,9u,40u,43u,44u};
                static const uint8_t kinds[] = {DUCKVEP_KIND_INS,DUCKVEP_KIND_INDEL,
                    DUCKVEP_KIND_INS,DUCKVEP_KIND_DEL,DUCKVEP_KIND_DEL,DUCKVEP_KIND_INS,
                    DUCKVEP_KIND_DEL,DUCKVEP_KIND_DEL};
                duckvep_haplotype_edit_t physical[2];
                uint8_t alt_cds[64], ref_peptide[32], alt_peptide[32];
                duckvep_delta_scratch_t scratch;
                size_t j;
                memset(&scratch, 0, sizeof scratch);
                scratch.edits = physical; scratch.edits_cap = 2u;
                scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
                scratch.ref_peptide = ref_peptide; scratch.ref_peptide_cap = sizeof ref_peptide;
                scratch.alt_peptide = alt_peptide; scratch.alt_peptide_cap = sizeof alt_peptide;
                for (j = 0u; j < 8u; j++) {
                    duckvep_event_t event;
                    duckvep_coding_context_t context;
                    duckvep_sequence_delta_t standalone;
                    size_t rlen = strlen(refs[j]), alen = strlen(alts[j]);
                    size_t b;
                    s.vkind = kinds[j]; s.rlen = (uint16_t)rlen;
                    s.alen = (uint16_t)alen; s.aoff = (uint32_t)rlen;
                    s.vpos = strand > 0 ? positions[j] : 351u - positions[j] - (uint32_t)rlen;
                    s.vend = s.vpos + (uint32_t)rlen - 1u;
                    for (b = 0u; b < rlen; b++) s.abytes[b] =
                        (uint8_t)coding_test_genomic_from_tx(refs[j][strand > 0 ? b : rlen-1u-b], strand);
                    for (b = 0u; b < alen; b++) s.abytes[rlen+b] =
                        (uint8_t)coding_test_genomic_from_tx(alts[j][strand > 0 ? b : alen-1u-b], strand);
                    if (strand < 0 && kinds[j] != DUCKVEP_KIND_INDEL) {
                        /* A reverse VCF needs its own left genomic anchor.
                         * Naively reversing TGGT>T or T>TAGGT lets prefix-first
                         * trimming choose the other retained T, a different
                         * VEP feature even though the biological edit agrees. */
                        uint32_t after = positions[j] + (uint32_t)rlen;
                        char anchor = after > 240u ? "ACGTACGTAC"[after-241u]
                            : (char)unpadded[after <= 180u ? after-158u : after-197u];
                        s.vpos = 350u - after; s.vend = s.vpos + (uint32_t)rlen - 1u;
                        s.abytes[0] = s.abytes[rlen] = (uint8_t)coding_test_genomic_from_tx(anchor, strand);
                        for (b = 1u; b < rlen; b++) s.abytes[b] =
                            (uint8_t)coding_test_genomic_from_tx(refs[j][rlen-b], strand);
                        for (b = 1u; b < alen; b++) s.abytes[rlen+b] =
                            (uint8_t)coding_test_genomic_from_tx(alts[j][alen-b], strand);
                    }
                    duckvep_event_load(&s.v, 0u, &event);
                    if (j < 7u) {
                        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                            duckvep_variant_feature_coding_context_build_prepared(
                                &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, strand, &event,
                                UINT32_MAX, NULL, physical, 2u, alt_cds, sizeof alt_cds,
                                ref_peptide, sizeof ref_peptide, alt_peptide, sizeof alt_peptide, &context));
                        if (physical[0].cds_start != unpadded_starts[j] + (uint32_t)phase)
                            fprintf(stderr, "[indel coordinate] case=%zu phase=%d strand=%d physical=%u feature=%u\n",
                                j, phase, strand, physical[0].cds_start, context.single_edit_cds_start);
                        ASSERT_EQ(unpadded_starts[j] + (uint32_t)phase, physical[0].cds_start);
                        ASSERT_EQ(unpadded_starts[j], context.single_edit_cds_start);
                        ASSERT_EQ(unpadded_starts[j], context.single_edit_unpadded_start1);
                        ASSERT_EQ((uint8_t)phase, context.cds_phase_padding);
                    }
                    duckvep_sequence_delta_fill_for_annotation((duckvep_variant_kind_t)s.vkind,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, &scratch, &event, &delta);
                    if (!delta.valid) fprintf(stderr, "[later-CDS indel] case=%zu phase=%d strand=%d status=%u\n",
                        j, phase, strand, delta.sequence_status);
                    ASSERT(delta.valid && !delta.synonymous && !delta.missense &&
                        !delta.start_retained && !delta.protein_altering);
                    ASSERT_EQ(j < 2u, delta.start_lost);
                    ASSERT_EQ((j == 2u && phase < 2) || (j == 3u && phase == 0), delta.stop_gained);
                    ASSERT_EQ((j == 2u && phase > 0) || (j == 5u && phase == 0), delta.inframe_insertion);
                    ASSERT_EQ(j == 3u || j == 4u, delta.inframe_deletion);
                    ASSERT_EQ(j < 2u || (j == 5u && phase > 0) || (j == 6u && phase == 2), delta.frameshift);
                    ASSERT_EQ((j == 5u && phase > 0) || (j == 6u && phase == 2) ||
                        (j == 7u && phase > 0), delta.stop_lost);
                    ASSERT_EQ(j == 6u && phase == 1, delta.stop_retained);
                    ASSERT_EQ(j >= 5u && phase == 0, delta.partial_codon);
                    ASSERT_EQ(j >= 5u && phase == 0, delta.coding_unknown);
                    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);
                    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s.vkind,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, &scratch, &standalone);
                    ASSERT_EQ(0, memcmp(&delta, &standalone, sizeof delta));
                    if (j == 1u) {
                        uint8_t saved = cds[(size_t)phase];
                        /* Genuine N immediately after the synthetic prefix is
                         * not known padding, even outside the validated REF. */
                        cds[(size_t)phase] = 'N';
                        duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s.vkind,
                            &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, &scratch, &delta);
                        ASSERT(!delta.valid);
                        cds[(size_t)phase] = saved;
                    }
                    /* Wrong anchor or changed REF must fail before rephasing. */
                    b = kinds[j] == DUCKVEP_KIND_DEL ? 1u : 0u;
                    s.abytes[b] = s.abytes[b] == 'A' ? 'C' : 'A';
                    if (kinds[j] == DUCKVEP_KIND_INS) s.abytes[rlen] = s.abytes[0];
                    if (kinds[j] == DUCKVEP_KIND_INDEL) s.abytes[0] = strand > 0 ? 'G' : 'C';
                    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s.vkind,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand, &scratch, &delta);
                    if (delta.valid) fprintf(stderr, "[indel accepted bad REF] case=%zu phase=%d strand=%d\n",
                        j, phase, strand);
                    ASSERT(!delta.valid);
                    if (delta.sequence_status != (uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH)
                        fprintf(stderr, "[indel bad REF] case=%zu phase=%d strand=%d status=%u\n",
                            j, phase, strand, delta.sequence_status);
                    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH, delta.sequence_status);
                }
            }
            {
                /* Complete uploaded CDS spans from the same pinned campaign:
                 * start, two-codon start, body, and physical CDS endpoint. */
                static const uint32_t mnv_positions[] = {158u, 160u, 165u, 239u, 235u};
                static const char *refs[] = {"CG", "TA", "AC", "AA", "TG"};
                static const char *alts[] = {"GC", "AT", "TG", "CA", "AC"};
                uint8_t alt_cds[64], ref_peptide[32], alt_peptide[32];
                duckvep_delta_scratch_t scratch;
                size_t j;
                memset(&scratch, 0, sizeof scratch);
                scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
                scratch.ref_peptide = ref_peptide; scratch.ref_peptide_cap = sizeof ref_peptide;
                scratch.alt_peptide = alt_peptide; scratch.alt_peptide_cap = sizeof alt_peptide;
                s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
                s.rlen = s.alen = 2u; s.aoff = 2u;
                for (j = 0u; j < 5u; j++) {
                    duckvep_event_t event;
                    duckvep_coding_context_t context;
                    duckvep_sequence_delta_t routed;
                    duckvep_feature_substitution_result_t result;
                    uint32_t unpadded = mnv_positions[j] <= 180u
                        ? mnv_positions[j] - 157u : mnv_positions[j] - 196u;
                    size_t b;
                    s.vpos = strand > 0 ? mnv_positions[j] : 349u - mnv_positions[j];
                    s.vend = s.vpos + 1u;
                    for (b = 0u; b < 2u; b++) {
                        size_t oriented = strand > 0 ? b : 1u - b;
                        s.abytes[b] = (uint8_t)coding_test_genomic_from_tx(refs[j][oriented], strand);
                        s.abytes[2u + b] = (uint8_t)coding_test_genomic_from_tx(alts[j][oriented], strand);
                    }
                    duckvep_event_load(&s.v, 0u, &event);
                    result = duckvep_feature_substitution_context_fill(
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, strand,
                        &scratch, &event, UINT32_MAX, &context, &delta);
                    ASSERT_EQ(DUCKVEP_FEATURE_SUBSTITUTION_CONTEXT_READY, result);
                    ASSERT(delta.valid && !delta.start_retained &&
                           !delta.stop_retained && !delta.frameshift && !delta.protein_altering);
                    ASSERT_EQ(j < 2u, delta.start_lost);
                    ASSERT_EQ((j == 2u && phase == 1) || (j == 4u && phase == 2), delta.synonymous);
                    ASSERT_EQ((j == 2u && phase != 1) || (j == 3u && phase == 2) ||
                              (j == 4u && phase == 0), delta.missense);
                    ASSERT_EQ(j == 4u && phase == 1, delta.stop_gained);
                    ASSERT_EQ(j == 3u && phase == 1, delta.stop_lost);
                    ASSERT_EQ(j == 3u && phase == 0, delta.partial_codon);
                    ASSERT_EQ(j == 3u && phase == 0, delta.coding_unknown);
                    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, delta.sequence_status);
                    if (result == DUCKVEP_FEATURE_SUBSTITUTION_CONTEXT_READY) {
                        ASSERT_EQ(unpadded, context.single_edit_cds_start);
                        ASSERT_EQ((uint8_t)phase, context.cds_phase_padding);
                        ASSERT_EQ(2u, context.single_edit_ref_len);
                        if (j == 4u && phase == 2) ASSERT_EQ(0u, context.cds_changed);
                    }
                    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand,
                        &scratch, NULL, &routed);
                    ASSERT_EQ(0, memcmp(&delta, &routed, sizeof delta));
                    /* Even partial-codon admission must not bypass REF. For
                     * AA>CA mutate the retained A in both alleles, preserving
                     * the smaller semantic edit while making REF wrong. */
                    b = strand > 0 ? 1u : 0u;
                    s.abytes[b] = (uint8_t)coding_test_genomic_from_tx(
                        refs[j][1] == 'A' ? 'G' : 'A', strand);
                    if (j == 3u) s.abytes[b + 2u] = s.abytes[b];
                    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV,
                        &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.vpos, strand,
                        &scratch, NULL, &delta);
                    ASSERT(!delta.valid && !delta.partial_codon);
                    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH, delta.sequence_status);
                }
            }
        }
    }
    PASS();
}

/* VEP's missense and coding_unknown predicates are independent. An equal-length
 * replacement on a CDS_START_NF transcript can preserve the leading X while
 * changing later amino acids, yielding both facts. This is the reduced form of
 * GRCh38 17:75629084 GATGCCAGCAGA>TCTGCCTCTGGG on ENST00000581825. */
TEST coding_context_incomplete_start_mnv_is_unknown_and_missense_known_scene(void) {
    static const uint8_t cds[] = {
        'N','N','T', 'C','T','G', 'C','T','G',
        'G','C','A', 'T','C','A', 'T','A','A'
    };
    static const uint8_t alt[] = {
        'C','C','C', 'A','G','A', 'G','G','C', 'A','G','A'
    };
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[sizeof cds];
    uint8_t ref_peptide[sizeof cds / 3u + 1u];
    uint8_t alt_peptide[sizeof cds / 3u + 1u];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 3u;
    edit.ref = cds + 2u;
    edit.ref_len = sizeof alt;
    edit.alt = alt;
    edit.alt_len = sizeof alt;
    edit.variant_strand = 1;
    edit_set.edits = &edit;
    edit_set.count = 1u;

    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, sizeof cds, &edit_set, 1,
                  DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds,
                  ref_peptide, sizeof ref_peptide,
                  alt_peptide, sizeof alt_peptide, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(
                  &ctx, (uint64_t)DUCKVEP_TX_CDS_START_NF, &delta));
    ASSERT(delta.valid && delta.coding_unknown && delta.missense &&
           !delta.synonymous && !delta.start_lost &&
           !delta.stop_gained && !delta.stop_lost && !delta.stop_retained);
    PASS();
}

static enum theft_alloc_res kprop_cds_edit_builder_alloc(struct theft *t, void *env,
                                                         void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t shape = (uint32_t)kprop_bounded(t, 5u);
    uint32_t region = (uint32_t)kprop_bounded(t, 3u);
    uint32_t ref_len = 1u;
    uint32_t alt_len = 1u;
    uint32_t cds_start = 1u;
    uint8_t alt_tx[8];
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    s->flags = 0u; s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    s->expect_shape = (uint8_t)shape;
    s->expect_region = (uint8_t)region;
    if (shape == KPROP_CDS_EDIT_SNV) {
        ref_len = 1u; alt_len = 1u;
        cds_start = kprop_pick_cds_start(t, cds_len, ref_len, region);
        alt_tx[0] = kprop_base_not(s->cds[cds_start - 1u], (uint8_t)'N');
        s->vkind = (uint8_t)DUCKVEP_KIND_SNV;
        s->vpos = kprop_genomic_pos_for_cds(s, cds_start); s->vend = s->vpos;
        s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
        kprop_fill_variant_alt_from_tx(s, 1u, alt_tx, alt_len);
        s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 1u;
        kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, alt_len);
    } else if (shape == KPROP_CDS_EDIT_MNV || shape == KPROP_CDS_EDIT_INDEL) {
        ref_len = (uint32_t)kprop_bounded(t, 3u) + 1u;
        alt_len = shape == KPROP_CDS_EDIT_MNV
            ? ref_len
            : (uint32_t)kprop_bounded(t, 4u) + 1u;
        if (shape == KPROP_CDS_EDIT_INDEL && alt_len == ref_len) alt_len = alt_len == 4u ? 1u : alt_len + 1u;
        cds_start = kprop_pick_cds_start(t, cds_len, ref_len, region);
        {
            uint8_t b0 = s->cds[cds_start - 1u];
            uint8_t b1 = ref_len > 1u ? s->cds[cds_start] : b0;
            uint8_t b2 = ref_len > 2u ? s->cds[cds_start + 1u] : b1;
            uint8_t alt_base = kprop_base_not3(b0, b1, b2);
            for (i = 0u; i < alt_len; i++) alt_tx[i] = alt_base;
        }
        s->vkind = shape == KPROP_CDS_EDIT_MNV ? (uint8_t)DUCKVEP_KIND_MNV : (uint8_t)DUCKVEP_KIND_INDEL;
        s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                                : kprop_genomic_pos_for_cds(s, cds_start + ref_len - 1u);
        s->vend = s->vpos + ref_len - 1u;
        for (i = 0u; i < ref_len; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
        kprop_fill_variant_alt_from_tx(s, ref_len, alt_tx, alt_len);
        s->roff = 0u; s->aoff = ref_len; s->rlen = (uint16_t)ref_len; s->alen = (uint16_t)alt_len;
        kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, alt_len);
    } else if (shape == KPROP_CDS_EDIT_INS) {
        uint32_t anchor_cds = region == KPROP_CDS_EDIT_START ? 1u
            : (region == KPROP_CDS_EDIT_STOP ? cds_len
                                             : (uint32_t)kprop_bounded(t, cds_len - 6u) + 4u);
        uint32_t insert_cds = s->strand > 0 ? anchor_cds + 1u : anchor_cds;
        alt_len = (uint32_t)kprop_bounded(t, 3u) + 1u;
        for (i = 0u; i < alt_len; i++) alt_tx[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];
        s->vkind = (uint8_t)DUCKVEP_KIND_INS;
        s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds); s->vend = s->vpos;
        s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
        s->abytes[1] = s->abytes[0];
        kprop_fill_variant_alt_from_tx(s, 2u, alt_tx, alt_len);
        s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = (uint16_t)(alt_len + 1u);
        kprop_fill_expected_cds(s, insert_cds, 0u, alt_tx, alt_len);
    } else {
        uint32_t anchor_cds;
        ref_len = (uint32_t)kprop_bounded(t, 3u) + 1u;
        cds_start = (uint32_t)kprop_bounded(t, cds_len - ref_len - 1u) + 2u;
        if (s->strand > 0) {
            anchor_cds = cds_start - 1u;
        } else {
            if (cds_start + ref_len > cds_len) cds_start = cds_len - ref_len;
            anchor_cds = cds_start + ref_len;
        }
        s->expect_region = cds_start <= 3u ? (uint8_t)KPROP_CDS_EDIT_START
            : (cds_start + ref_len - 1u > cds_len - 3u ? (uint8_t)KPROP_CDS_EDIT_STOP
                                                        : (uint8_t)KPROP_CDS_EDIT_BODY);
        s->vkind = (uint8_t)DUCKVEP_KIND_DEL;
        s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
        s->vend = s->vpos + ref_len;
        for (i = 0u; i <= ref_len; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
        s->abytes[ref_len + 1u] = s->abytes[0];
        s->roff = 0u; s->aoff = ref_len + 1u; s->rlen = (uint16_t)(ref_len + 1u); s->alen = 1u;
        kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, 0u);
    }

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_cds_edit_builder_info = {
    .alloc = kprop_cds_edit_builder_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_cds_edit_set_mnv_alloc(struct theft *t, void *env,
                                                         void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t len = (uint32_t)kprop_bounded(t, 5u) + 3u;
    uint32_t region = (uint32_t)kprop_bounded(t, 3u);
    uint32_t cds_start;
    uint8_t alt_tx[8];
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    s->flags = 0u; s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    cds_start = kprop_pick_cds_start(t, cds_len, len, region);
    s->expect_shape = KPROP_CDS_EDIT_MNV;
    s->expect_region = (uint8_t)region;
    for (i = 0u; i < len; i++) {
        uint8_t ref = s->cds[cds_start - 1u + i];
        int diff = (i == 0u || i + 1u == len) ? 1 : (kprop_bounded(t, 2u) == 0u);
        if (i == 1u) diff = 0; /* force at least one retained internal base */
        alt_tx[i] = diff ? kprop_base_not(ref, (uint8_t)'N') : ref;
    }
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                            : kprop_genomic_pos_for_cds(s, cds_start + len - 1u);
    s->vend = s->vpos + len - 1u;
    for (i = 0u; i < len; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    kprop_fill_variant_alt_from_tx(s, len, alt_tx, len);
    s->roff = 0u; s->aoff = len; s->rlen = (uint16_t)len; s->alen = (uint16_t)len;
    kprop_fill_expected_cds(s, cds_start, len, alt_tx, len);
    *instance = s;
    return THEFT_ALLOC_OK;
}

struct theft_type_info kprop_cds_edit_set_mnv_info = {
    .alloc = kprop_cds_edit_set_mnv_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_context_delta_alloc(struct theft *t, void *env,
                                                      void **instance) {
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t category = (uint32_t)kprop_bounded(t, 5u);
    uint32_t codon_idx = (uint32_t)kprop_bounded(t, ncodons - 1u) + 1u;
    uint32_t cds_start = codon_idx * 3u + 1u;
    const uint8_t *ref_tx;
    const uint8_t *alt_tx;
    uint32_t i;
    static const uint8_t syn_ref[3] = {'G','A','A'};
    static const uint8_t syn_alt[3] = {'G','A','G'};
    static const uint8_t mis_ref[3] = {'G','A','A'};
    static const uint8_t mis_alt[3] = {'G','A','C'};
    static const uint8_t sg_ref[3] = {'T','G','G'};
    static const uint8_t sg_alt[3] = {'T','G','A'};
    static const uint8_t sl_ref[3] = {'T','A','A'};
    static const uint8_t sl_alt[3] = {'C','A','A'};
    static const uint8_t sr_ref[3] = {'T','A','A'};
    static const uint8_t sr_alt[3] = {'T','A','G'};
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    s->cds[0] = 'A'; s->cds[1] = 'T'; s->cds[2] = 'G';
    for (i = 3u; i < cds_len; i += 3u) {
        s->cds[i] = 'G'; s->cds[i + 1u] = 'A'; s->cds[i + 2u] = 'A';
    }
    if (category == KPROP_CONTEXT_DELTA_SYNONYMOUS) {
        ref_tx = syn_ref; alt_tx = syn_alt;
    } else if (category == KPROP_CONTEXT_DELTA_MISSENSE) {
        ref_tx = mis_ref; alt_tx = mis_alt;
    } else if (category == KPROP_CONTEXT_DELTA_STOP_GAINED) {
        ref_tx = sg_ref; alt_tx = sg_alt;
    } else if (category == KPROP_CONTEXT_DELTA_STOP_LOST) {
        ref_tx = sl_ref; alt_tx = sl_alt;
    } else {
        ref_tx = sr_ref; alt_tx = sr_alt;
    }
    memcpy(s->cds + (size_t)cds_start - 1u, ref_tx, 3u);

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    s->flags = 0u; s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    s->expect_shape = KPROP_CDS_EDIT_MNV;
    s->expect_region = (uint8_t)category;
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                            : kprop_genomic_pos_for_cds(s, cds_start + 2u);
    s->vend = s->vpos + 2u;
    for (i = 0u; i < 3u; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    kprop_fill_variant_alt_from_tx(s, 3u, alt_tx, 3u);
    s->roff = 0u; s->aoff = 3u; s->rlen = 3u; s->alen = 3u;
    kprop_fill_expected_cds(s, cds_start, 3u, alt_tx, 3u);
    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_context_delta_info = {
    .alloc = kprop_context_delta_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_cross_codon_mnv_alloc(struct theft *t, void *env, void **instance) {
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t codon_start;
    uint32_t first_cds;
    uint32_t mode;
    uint32_t len;
    char alt_tx[3];
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)'A';

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 3u) + 1u; /* body codon, plus next body codon */
    codon_start = codon_idx * 3u + 1u;
    mode = (uint32_t)kprop_bounded(t, 4u);
    if (mode < 2u) {
        /* AAA/TTA, editing last base of codon1 + first of codon2.
         * mode 0: one synonymous codon and one missense codon.
         * mode 1: both codons synonymous, so this narrow slice must fall back. */
        s->cds[(size_t)codon_start - 1u] = (uint8_t)'A';
        s->cds[(size_t)codon_start] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 1u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 2u] = (uint8_t)'T';
        s->cds[(size_t)codon_start + 3u] = (uint8_t)'T';
        s->cds[(size_t)codon_start + 4u] = (uint8_t)'A';
        first_cds = codon_start + 2u;
        len = 2u;
        alt_tx[0] = 'G';
        alt_tx[1] = mode == 0u ? 'G' : 'C';
    } else {
        /* TCA/AAA, editing last two bases of codon1 + first of codon2.
         * mode 2 introduces TAG stop and must fall back; mode 3 is non-stop missense. */
        s->cds[(size_t)codon_start - 1u] = (uint8_t)'T';
        s->cds[(size_t)codon_start] = (uint8_t)'C';
        s->cds[(size_t)codon_start + 1u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 2u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 3u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 4u] = (uint8_t)'A';
        first_cds = codon_start + 1u;
        len = 3u;
        if (mode == 2u) {
            alt_tx[0] = 'A'; alt_tx[1] = 'G'; alt_tx[2] = 'C';
        } else {
            alt_tx[0] = 'G'; alt_tx[1] = 'G'; alt_tx[2] = 'C';
        }
    }

    s->vchrom = 0u;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, first_cds)
                            : kprop_genomic_pos_for_cds(s, first_cds + len - 1u);
    s->vend = s->vpos + len - 1u;
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < len; i++) {
        uint32_t gpos = s->vpos + i;
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, gpos);
        uint32_t tx_i = cds_pos - first_cds;
        char alt = alt_tx[tx_i];
        s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, gpos);
        s->abytes[len + i] = (uint8_t)(s->strand > 0 ? alt : kprop_complement_base(alt));
    }
    s->roff = 0u; s->aoff = len; s->rlen = (uint16_t)len; s->alen = (uint16_t)len;

    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = (size_t)len * 2u;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_cross_codon_mnv_info = {
    .alloc = kprop_cross_codon_mnv_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_frameshift_indel_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    static const uint8_t STOPS[3][3] = {
        {'T','A','A'}, {'T','A','G'}, {'T','G','A'}
    };
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t g_off;
    uint32_t i;
    uint32_t extra;
    uint32_t cds_pos;
    uint32_t mode;
    int force_terminal;
    int force_terminal_nonstop;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];
    {
        uint32_t terminal_scene = (uint32_t)kprop_bounded(t, 8u);
        force_terminal = terminal_scene < 2u;
        force_terminal_nonstop = terminal_scene == 0u;
    }
    s->cds[0] = (uint8_t)'A';
    s->cds[1] = (uint8_t)'T';
    s->cds[2] = (uint8_t)'G';
    for (i = 3u; i + 3u < cds_len; i += 3u) {
        if (duckvep_translate_codon(
                (const char *)(s->cds + i),
                DUCKVEP_CODON_TABLE_STANDARD) == '*') {
            s->cds[i + 1u] = (uint8_t)'C';
        }
    }
    if (force_terminal) {
        if (force_terminal_nonstop) {
            static const uint8_t TERMINAL_NONSTOP[3] = {'C','A','A'};
            memcpy(s->cds + cds_len - 3u, TERMINAL_NONSTOP, 3u);
        } else {
            uint32_t stop = (uint32_t)kprop_bounded(t, 3u);
            memcpy(s->cds + cds_len - 3u, STOPS[stop], 3u);
        }
    }

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    mode = force_terminal_nonstop ? 0u : (uint32_t)kprop_bounded(t, 3u);
    if (mode == 0u) {
        extra = (uint32_t)kprop_bounded(t, 2u) + 1u; /* net +1 or +2, both frameshift */
        if (force_terminal) {
            cds_pos = s->strand > 0 ? cds_len - 2u : cds_len - 1u;
        } else {
            cds_pos = (uint32_t)kprop_bounded(t, cds_len - 6u) + 4u;
        }
        s->vpos = kprop_genomic_pos_for_cds(s, cds_pos);
        s->vend = s->vpos; s->vkind = (uint8_t)DUCKVEP_KIND_INS;
        s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
        s->abytes[1] = s->abytes[0];
        for (i = 0u; i < extra; i++) s->abytes[2u + i] = (uint8_t)BASES[kprop_bounded(t, 4u)];
        s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = (uint16_t)(1u + extra);
    } else if (mode == 1u) {
        extra = (uint32_t)kprop_bounded(t, 2u) + 1u; /* delete 1 or 2 bases after anchor */
        if (force_terminal) {
            /* Exercise the fail-closed VEP endpoint state where a deletion
             * removes the final CDS bases and no complete transcript tail is
             * available for reconstructing the original stop codon. Keep this
             * stratum forward because the minimal one-exon generator has no
             * genomic base beyond the reverse-strand CDS for the VCF anchor. */
            s->strand = 1;
            cds_pos = cds_len - extra + 1u;
        } else {
            cds_pos = (uint32_t)kprop_bounded(t, cds_len - extra - 5u) + 4u;
        }
        if (s->strand > 0) {
            s->vpos = kprop_genomic_pos_for_cds(s, cds_pos - 1u);
        } else {
            s->vpos = kprop_genomic_pos_for_cds(s, cds_pos + extra);
        }
        g_off = s->vpos - base;
        s->vend = s->vpos + extra; s->vkind = (uint8_t)DUCKVEP_KIND_DEL;
        for (i = 0u; i <= extra; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, base + g_off + i);
        s->abytes[extra + 1u] = s->abytes[0];
        s->roff = 0u; s->aoff = extra + 1u; s->rlen = (uint16_t)(extra + 1u); s->alen = 1u;
    } else {
        uint32_t delins_mode = (uint32_t)kprop_bounded(t, 4u);
        uint32_t ref_len = delins_mode < 2u ? 1u : delins_mode;
        uint32_t alt_len = delins_mode < 2u ? delins_mode + 2u : 1u;
        uint32_t first_cds = (uint32_t)kprop_bounded(t, cds_len - ref_len - 5u) + 4u;
        s->vpos = s->strand > 0
            ? kprop_genomic_pos_for_cds(s, first_cds)
            : kprop_genomic_pos_for_cds(s, first_cds + ref_len - 1u);
        s->vend = s->vpos + ref_len - 1u;
        s->vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < ref_len; i++) {
            s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
        }
        for (i = 0u; i < alt_len; i++) {
            s->abytes[ref_len + i] = kprop_base_not(s->abytes[0], s->abytes[ref_len - 1u]);
        }
        s->roff = 0u; s->aoff = ref_len;
        s->rlen = (uint16_t)ref_len; s->alen = (uint16_t)alt_len;
    }

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_frameshift_indel_info = {
    .alloc = kprop_frameshift_indel_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_inframe_deletion_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 3u) + 4u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t del_start_cds;
    uint32_t del_end_cds;
    uint32_t anchor_cds;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 2u) + 1u; /* exclude first/last codon */
    del_start_cds = codon_idx * 3u + 1u;
    del_end_cds = del_start_cds + 2u;
    /* Keep the deleted codon non-stop: real coding transcripts have no internal stop, and
     * removing an internal stop is stop_lost (not a clean in-frame deletion) — the classifier
     * defers those, so exclude them here to keep this generator to clean in-frame deletions.
     * TAA/TAG/TGA -> T[C]A/T[C]G/T[C]A (Ser), a minimal non-stop rewrite of the middle base. */
    {
        uint8_t *dc = &s->cds[(size_t)del_start_cds - 1u];
        if (dc[0] == (uint8_t)'T' &&
            ((dc[1] == (uint8_t)'A' && (dc[2] == (uint8_t)'A' || dc[2] == (uint8_t)'G')) ||
             (dc[1] == (uint8_t)'G' && dc[2] == (uint8_t)'A'))) {
            dc[1] = (uint8_t)'C';
        }
    }
    anchor_cds = s->strand > 0 ? del_start_cds - 1u : del_end_cds + 1u;
    s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
    s->vend = s->vpos + 3u;
    s->vkind = (uint8_t)DUCKVEP_KIND_DEL;
    for (i = 0u; i < 4u; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    s->abytes[4] = s->abytes[0];
    s->roff = 0u; s->aoff = 4u; s->rlen = 4u; s->alen = 1u;

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_inframe_deletion_info = {
    .alloc = kprop_inframe_deletion_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_inframe_insertion_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    static const char INSERTED[3] = {'G', 'C', 'C'}; /* alanine, not a stop codon */
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t before_cds;
    uint32_t anchor_cds;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 4u) + 2u; /* after codon 2..n-2 */
    before_cds = codon_idx * 3u;
    anchor_cds = s->strand > 0 ? before_cds : before_cds + 1u;
    s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
    s->vend = s->vpos;
    s->vkind = (uint8_t)DUCKVEP_KIND_INS;
    s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
    s->abytes[1] = s->abytes[0];
    for (i = 0u; i < 3u; i++) {
        char b = s->strand > 0 ? INSERTED[i] : kprop_complement_base(INSERTED[2u - i]);
        s->abytes[2u + i] = (uint8_t)b;
    }
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 4u;

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_inframe_insertion_info = {
    .alloc = kprop_inframe_insertion_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_delins_shape_alloc(struct theft *t, void *env, void **instance) {
    static const uint8_t ALT3_PRESERVED[3] = { 'A', 'A', 'A' };
    static const uint8_t ALT3_ALTERED[3] = { 'G', 'C', 'C' };
    static const uint8_t ALT6_PRESERVED[6] = { 'A', 'A', 'G', 'G', 'C', 'C' };
    static const uint8_t ALT6_ALTERED[6] = { 'G', 'C', 'C', 'G', 'C', 'T' };
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t lengthening = (uint32_t)kprop_bounded(t, 2u) == 0u;
    uint32_t preserved = (uint32_t)kprop_bounded(t, 2u) == 0u;
    uint32_t codon_idx;
    uint32_t cds_start;
    uint32_t ref_len;
    uint32_t alt_len;
    const uint8_t *alt_tx;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)'A';

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    ref_len = lengthening ? 3u : 6u;
    alt_len = lengthening ? 6u : 3u;
    alt_tx = lengthening
        ? (preserved ? ALT6_PRESERVED : ALT6_ALTERED)
        : (preserved ? ALT3_PRESERVED : ALT3_ALTERED);
    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 3u) + 1u;
    cds_start = codon_idx * 3u + 1u;
    s->expect_shape = lengthening ? (uint8_t)KPROP_CDS_EDIT_INS
                                  : (uint8_t)KPROP_CDS_EDIT_DEL;
    s->expect_region = (uint8_t)KPROP_CDS_EDIT_BODY;
    s->expect_protein_altering = preserved ? 0u : 1u;
    s->vkind = (uint8_t)DUCKVEP_KIND_INDEL;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                            : kprop_genomic_pos_for_cds(s, cds_start + ref_len - 1u);
    s->vend = s->vpos + ref_len - 1u;
    for (i = 0u; i < ref_len; i++) {
        s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    }
    kprop_fill_variant_alt_from_tx(s, ref_len, alt_tx, alt_len);
    s->roff = 0u; s->aoff = ref_len;
    s->rlen = (uint16_t)ref_len; s->alen = (uint16_t)alt_len;
    kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, alt_len);

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_delins_shape_info = {
    .alloc = kprop_delins_shape_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_protein_altering_insertion_alloc(struct theft *t, void *env, void **instance) {
    static const char INSERTED[3] = {'G', 'C', 'C'}; /* AAA + inserted GCC stays non-stop */
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t before_cds;
    uint32_t anchor_cds;
    uint32_t i;
    uint32_t off;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)'A';

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 2u) + 1u; /* codon 2..n-1 */
    off = (uint32_t)kprop_bounded(t, 2u) + 1u;                 /* after base 1 or 2 */
    before_cds = codon_idx * 3u + off;
    anchor_cds = s->strand > 0 ? before_cds : before_cds + 1u;
    s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
    s->vend = s->vpos;
    s->vkind = (uint8_t)DUCKVEP_KIND_INS;
    s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
    s->abytes[1] = s->abytes[0];
    for (i = 0u; i < 3u; i++) {
        char b = s->strand > 0 ? INSERTED[i] : kprop_complement_base(INSERTED[2u - i]);
        s->abytes[2u + i] = (uint8_t)b;
    }
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 4u;

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_protein_altering_insertion_info = {
    .alloc = kprop_protein_altering_insertion_alloc,
    .free  = kprop_coding_free,
};

static uint64_t codon_change_to_so(uint32_t change, char aa_ref) {
    if (change & DUCKVEP_CODON_STOP_GAINED)     return DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    if (change & DUCKVEP_CODON_STOP_LOST)       return DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    if (change & DUCKVEP_CODON_MISSENSE)        return DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
    if (change & DUCKVEP_CODON_SYNONYMOUS) {
        return (aa_ref == '*') ? DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED)
                               : DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    }
    return 0u;
}

static uint64_t start_codon_snv_to_so(
    uint32_t change, char aa_ref, char aa_alt, const char alt_codon[4],
    int has_five_prime_utr) {

    uint64_t mask;
    int retained = alt_codon != NULL && alt_codon[0] == 'A' &&
        alt_codon[1] == 'T' && alt_codon[2] == 'G';

    mask = codon_change_to_so(change, aa_ref);
    if (aa_alt != aa_ref || (has_five_prime_utr && !retained)) {
        mask &= ~DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
        mask |= DUCKVEP_SO(DUCKVEP_SO_START_LOST);
    }
    if (retained) mask |= DUCKVEP_SO(DUCKVEP_SO_START_RETAINED);
    return mask;
}

/* Independent VEP-116 oracle for one codon selected by an equal-length uploaded
 * feature. Unlike a minimized SNV, the start predicates inspect the complete
 * feature: unequal local peptides make start_lost true, while the rebuilt ATG
 * independently makes start_retained_variant true. */
static uint64_t vep_feature_codon_to_so(
    char ref_aa, char alt_aa, const char alt_codon[4], int overlaps_start) {

    uint64_t mask = 0u;
    int start_lost = overlaps_start && ref_aa != alt_aa;
    int start_retained = overlaps_start &&
        alt_codon[0] == 'A' && alt_codon[1] == 'T' && alt_codon[2] == 'G';

    if (start_lost) mask |= DUCKVEP_SO(DUCKVEP_SO_START_LOST);
    if (start_retained) mask |= DUCKVEP_SO(DUCKVEP_SO_START_RETAINED);
    if (ref_aa == '*' && alt_aa != '*') {
        mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    } else if (ref_aa == '*' && alt_aa == '*') {
        mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
    } else if (ref_aa != '*' && alt_aa == '*') {
        mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    } else if (ref_aa == alt_aa) {
        if (!start_retained) mask |= DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    } else if (!start_lost) {
        mask |= DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
    }
    return mask;
}

/* Coverage accounting so "across codon-change classes" is a checked claim, not a
 * hope: the run must actually observe each class at least once (the framework's
 * coverage signal, applied locally). */
static struct { uint32_t syn; uint32_t mis; uint32_t sg; uint32_t sl; uint32_t sr; } g_codon_cov;

static enum theft_trial_res prop_annotate_codon_matches_kernel(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_coding_projection_t proj;
    duckvep_coding_snv_result_t res;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    /* Oracle: the tested kernels invoked directly. */
    if (!duckvep_project_coding_base(&s->tx, &s->ex, 0u, s->vpos, &proj)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_coding_snv_from_cds(s->cds, (size_t)s->cds_lenv, &proj, proj.cds_pos,
                                    (char)s->abytes[0], (char)s->abytes[1], 1,
                                    DUCKVEP_CODON_TABLE_STANDARD, &res) != DUCKVEP_CODING_SNV_OK) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

    if (res.change & DUCKVEP_CODON_STOP_GAINED)     g_codon_cov.sg++;
    else if (res.change & DUCKVEP_CODON_STOP_LOST)  g_codon_cov.sl++;
    else if (res.change & DUCKVEP_CODON_MISSENSE)   g_codon_cov.mis++;
    else if (res.change & DUCKVEP_CODON_SYNONYMOUS) {
        if (res.aa_ref == '*') g_codon_cov.sr++; else g_codon_cov.syn++;
    }

    if (rows[0].consequence_mask != codon_change_to_so(res.change, res.aa_ref)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cds_pos != (int32_t)res.cds_pos) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)res.protein_pos) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_ref != (uint8_t)res.aa_ref) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_alt != (uint8_t)res.aa_alt) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

static struct {
    uint32_t sl;
    uint32_t sl_sg;
    uint32_t sl_syn;
    uint32_t syn;
    uint32_t retained;
    uint32_t lost_and_retained;
} g_start_cov;

static enum theft_trial_res prop_annotate_start_lost_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_coding_projection_t proj;
    duckvep_coding_snv_result_t res;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    if (!duckvep_project_coding_base(&s->tx, &s->ex, 0u, s->vpos, &proj)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (proj.codon_start_cds != 1u || proj.protein_pos != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_coding_snv_from_cds(s->cds, (size_t)s->cds_lenv, &proj, proj.cds_pos,
                                    (char)s->abytes[0], (char)s->abytes[1], 1,
                                    DUCKVEP_CODON_TABLE_STANDARD, &res) != DUCKVEP_CODING_SNV_OK) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    want = start_codon_snv_to_so(
        res.change, res.aa_ref, res.aa_alt, res.alt_codon,
        proj.cdna_pos > proj.cds_pos);
    if (want & DUCKVEP_SO(DUCKVEP_SO_START_RETAINED)) {
        g_start_cov.retained++;
    }
    if (want & DUCKVEP_SO(DUCKVEP_SO_START_LOST)) {
        g_start_cov.sl++;
        if (want & DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)) g_start_cov.sl_sg++;
        if (want & DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)) g_start_cov.sl_syn++;
        if (want & DUCKVEP_SO(DUCKVEP_SO_START_RETAINED)) {
            g_start_cov.lost_and_retained++;
        }
    }
    if (want & DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)) g_start_cov.syn++;

    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != 1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_ref != (uint8_t)res.aa_ref) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_alt != (uint8_t)res.aa_alt) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_codon_matches_kernel_for_any_cds_snv(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile codon refinement == coding-SNV kernel oracle";
    cfg.prop1 = prop_annotate_codon_matches_kernel;
    cfg.type_info[0] = &kprop_coding_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_codon_cov, 0, sizeof g_codon_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    /* The "across codon-change classes" claim is only honest if the run actually
     * hit each class. (Deterministic under the fixed seed.) */
    ASSERT(g_codon_cov.syn > 0u);
    ASSERT(g_codon_cov.mis > 0u);
    ASSERT(g_codon_cov.sg > 0u);
    ASSERT(g_codon_cov.sl > 0u);
    ASSERT(g_codon_cov.sr > 0u); /* stop->stop synonymous reaches stop_retained */
    fprintf(stderr, "[codon coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u\n",
            g_codon_cov.syn, g_codon_cov.mis, g_codon_cov.sg, g_codon_cov.sl, g_codon_cov.sr);
    PASS();
}

static struct { uint32_t len2; uint32_t len3; } g_mnv_cov;

static enum theft_trial_res prop_annotate_mnv_matches_codon_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    uint32_t cds0;
    uint32_t codon_start;
    uint32_t protein_pos;
    uint8_t off0;
    uint32_t i;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    /* Independent single-exon '+' oracle: compute the CDS/codon coordinates directly,
     * edit the codon slice, then classify via the canonical genetic-code helper. */
    cds0 = s->vpos - s->tstart + 1u;
    codon_start = ((cds0 - 1u) / 3u) * 3u + 1u;
    protein_pos = ((codon_start - 1u) / 3u) + 1u;
    off0 = (uint8_t)(cds0 - codon_start);
    ref_codon[0] = (char)s->cds[codon_start - 1u];
    ref_codon[1] = (char)s->cds[codon_start];
    ref_codon[2] = (char)s->cds[codon_start + 1u];
    ref_codon[3] = '\0';
    alt_codon[0] = ref_codon[0]; alt_codon[1] = ref_codon[1]; alt_codon[2] = ref_codon[2]; alt_codon[3] = '\0';
    for (i = 0u; i < (uint32_t)s->alen; i++) {
        alt_codon[off0 + i] = (char)s->abytes[s->aoff + i];
    }
    cr = duckvep_codon_change(ref_codon, alt_codon, DUCKVEP_CODON_TABLE_STANDARD);
    want = (codon_start == 1u && protein_pos == 1u)
        ? vep_feature_codon_to_so(cr.aa_ref, cr.aa_alt, alt_codon, 1)
        : codon_change_to_so(cr.change, cr.aa_ref);

    if (s->alen == 2u) g_mnv_cov.len2++; else if (s->alen == 3u) g_mnv_cov.len3++;
    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)protein_pos) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_ref != (uint8_t)cr.aa_ref) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_alt != (uint8_t)cr.aa_alt) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_mnv_same_codon_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile same-codon MNV == codon oracle";
    cfg.prop1 = prop_annotate_mnv_matches_codon_oracle;
    cfg.type_info[0] = &kprop_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_mnv_cov, 0, sizeof g_mnv_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_mnv_cov.len2 > 0u);
    ASSERT(g_mnv_cov.len3 > 0u);
    fprintf(stderr, "[mnv coverage] len2=%u len3=%u\n", g_mnv_cov.len2, g_mnv_cov.len3);
    PASS();
}

static int kprop_cross_codon_mnv_oracle_missense(const struct kprop_coding *s) {
    uint32_t codon_start[2] = {0u, 0u};
    uint8_t seen[2] = {0u, 0u};
    char ref_codon[2][4];
    char alt_codon[2][4];
    uint32_t codon_count = 0u;
    uint32_t j;
    int changed = 0;

    for (j = 0u; j < (uint32_t)s->rlen; j++) {
        uint32_t gpos = s->vpos + j;
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, gpos);
        uint32_t start = ((cds_pos - 1u) / 3u) * 3u + 1u;
        uint32_t off = cds_pos - start;
        uint32_t k;
        uint32_t slot = UINT32_MAX;
        uint8_t bit;
        char ref_tx;
        char alt_tx;

        if (start <= 1u || start + 2u > s->cds_lenv - 3u) return 0;
        for (k = 0u; k < codon_count; k++) {
            if (codon_start[k] == start) { slot = k; break; }
        }
        if (slot == UINT32_MAX) {
            if (codon_count >= 2u) return 0;
            slot = codon_count++;
            codon_start[slot] = start;
            ref_codon[slot][0] = (char)s->cds[start - 1u];
            ref_codon[slot][1] = (char)s->cds[start];
            ref_codon[slot][2] = (char)s->cds[start + 1u];
            ref_codon[slot][3] = '\0';
            memcpy(alt_codon[slot], ref_codon[slot], sizeof alt_codon[slot]);
        }
        bit = (uint8_t)(1u << off);
        if (seen[slot] & bit) return 0;
        seen[slot] |= bit;
        ref_tx = s->strand > 0 ? (char)s->abytes[s->roff + j]
                               : kprop_complement_base((char)s->abytes[s->roff + j]);
        alt_tx = s->strand > 0 ? (char)s->abytes[s->aoff + j]
                               : kprop_complement_base((char)s->abytes[s->aoff + j]);
        if (ref_codon[slot][off] != ref_tx) return 0;
        alt_codon[slot][off] = alt_tx;
    }
    if (codon_count != 2u) return 0;
    if (codon_start[0] > codon_start[1]) {
        uint32_t tmp = codon_start[0];
        codon_start[0] = codon_start[1];
        codon_start[1] = tmp;
    }
    if (codon_start[1] - codon_start[0] != 3u) return 0;

    for (j = 0u; j < 2u; j++) {
        duckvep_codon_result_t cr = duckvep_codon_change(
            ref_codon[j], alt_codon[j], DUCKVEP_CODON_TABLE_STANDARD);
        if (cr.change & DUCKVEP_CODON_INVALID) return 0;
        if (cr.aa_ref == '*' || cr.aa_alt == '*') return 0;
        if (cr.aa_ref != cr.aa_alt) changed = 1;
    }
    return changed;
}

/* Full independent two-codon-window oracle for a body MNV: re-collects the (<=2) affected
 * body codons exactly like kprop_cross_codon_mnv_oracle_missense, then classifies the whole
 * window with VEP's precedence (stop_lost > stop_retained > stop_gained > synonymous/
 * missense) from the genetic code. `.valid == 0` means genuinely unsupported (N base, invalid
 * codon, or not two adjacent body codons) — where the kernel falls back to
 * coding_sequence_variant. Since the generator never touches codon 0, start_lost never
 * applies here. This is the oracle the generalized window classifier is graded against. */
struct kprop_cc_facts {
    int valid;
    int synonymous, missense, stop_gained, stop_lost, stop_retained;
};

static struct kprop_cc_facts kprop_cross_codon_mnv_oracle_facts(const struct kprop_coding *s) {
    struct kprop_cc_facts f;
    uint32_t codon_start[2] = {0u, 0u};
    uint8_t seen[2] = {0u, 0u};
    char ref_codon[2][4];
    char alt_codon[2][4];
    uint32_t codon_count = 0u;
    uint32_t j;
    int has_ref_stop = 0;
    int has_alt_stop = 0;
    int ref_stop_idx = -1;
    int alt_stop_idx = -1;
    int win_equal = 1;

    memset(&f, 0, sizeof f);
    for (j = 0u; j < (uint32_t)s->rlen; j++) {
        uint32_t gpos = s->vpos + j;
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, gpos);
        uint32_t start = ((cds_pos - 1u) / 3u) * 3u + 1u;
        uint32_t off = cds_pos - start;
        uint32_t k;
        uint32_t slot = UINT32_MAX;
        uint8_t bit;
        char ref_tx;
        char alt_tx;

        if (start <= 1u || start + 2u > s->cds_lenv - 3u) return f;
        for (k = 0u; k < codon_count; k++) {
            if (codon_start[k] == start) { slot = k; break; }
        }
        if (slot == UINT32_MAX) {
            if (codon_count >= 2u) return f;
            slot = codon_count++;
            codon_start[slot] = start;
            ref_codon[slot][0] = (char)s->cds[start - 1u];
            ref_codon[slot][1] = (char)s->cds[start];
            ref_codon[slot][2] = (char)s->cds[start + 1u];
            ref_codon[slot][3] = '\0';
            memcpy(alt_codon[slot], ref_codon[slot], sizeof alt_codon[slot]);
        }
        bit = (uint8_t)(1u << off);
        if (seen[slot] & bit) return f;
        seen[slot] |= bit;
        ref_tx = s->strand > 0 ? (char)s->abytes[s->roff + j]
                               : kprop_complement_base((char)s->abytes[s->roff + j]);
        alt_tx = s->strand > 0 ? (char)s->abytes[s->aoff + j]
                               : kprop_complement_base((char)s->abytes[s->aoff + j]);
        if (ref_codon[slot][off] != ref_tx) return f;
        alt_codon[slot][off] = alt_tx;
    }
    if (codon_count != 2u) return f;
    if (codon_start[0] > codon_start[1]) {
        uint32_t tmp = codon_start[0];
        codon_start[0] = codon_start[1];
        codon_start[1] = tmp;
    }
    if (codon_start[1] - codon_start[0] != 3u) return f;

    for (j = 0u; j < 2u; j++) {
        duckvep_codon_result_t cr = duckvep_codon_change(
            ref_codon[j], alt_codon[j], DUCKVEP_CODON_TABLE_STANDARD);
        if (cr.change & DUCKVEP_CODON_INVALID) return f;
        if (cr.aa_ref == 'X' || cr.aa_alt == 'X') return f;
        if (cr.aa_ref != cr.aa_alt) win_equal = 0;
        if (cr.aa_ref == '*' && !has_ref_stop) { ref_stop_idx = (int)j; has_ref_stop = 1; }
        if (cr.aa_alt == '*' && !has_alt_stop) { alt_stop_idx = (int)j; has_alt_stop = 1; }
    }
    f.valid = 1;
    if (has_ref_stop && !has_alt_stop) f.stop_lost = 1;
    else if (has_ref_stop && has_alt_stop && ref_stop_idx == alt_stop_idx) f.stop_retained = 1;
    else if (has_alt_stop && !has_ref_stop) f.stop_gained = 1;
    else if (win_equal) f.synonymous = 1;
    else f.missense = 1;
    return f;
}

/* SO consequence mask that annotate_tile emits for a resolved two-codon window (coarse:
 * protein_pos -1, no AA pair), or coding_sequence_variant when the oracle is unsupported. */
static uint64_t kprop_cross_codon_facts_to_so(const struct kprop_cc_facts *f) {
    if (!f->valid) return DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE);
    if (f->stop_gained) return DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    if (f->stop_lost) return DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    if (f->stop_retained) return DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
    if (f->synonymous) return DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    return DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
}

static struct { uint32_t missense; uint32_t synonymous; uint32_t stop_gained; uint32_t fwd; uint32_t rev; uint32_t len2; uint32_t len3; } g_cross_mnv_cov;

static enum theft_trial_res prop_annotate_cross_codon_mnv_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    struct kprop_cc_facts f;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    f = kprop_cross_codon_mnv_oracle_facts(s);
    want = kprop_cross_codon_facts_to_so(&f);
    if (f.missense) g_cross_mnv_cov.missense++;
    else if (f.synonymous) g_cross_mnv_cov.synonymous++;
    else if (f.stop_gained) g_cross_mnv_cov.stop_gained++;
    if (s->strand > 0) g_cross_mnv_cov.fwd++; else g_cross_mnv_cov.rev++;
    if (s->rlen == 2u) g_cross_mnv_cov.len2++; else if (s->rlen == 3u) g_cross_mnv_cov.len3++;

    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (stats->substitution_context != 1u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1 || rows[0].protein_pos != -1) {
        tr = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (rows[0].aa_ref != (uint8_t)0u || rows[0].aa_alt != (uint8_t)0u) {
        tr = THEFT_TRIAL_FAIL;
        goto done;
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_cross_codon_mnv_missense_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile two-codon body MNV missense == codon-window oracle";
    cfg.prop1 = prop_annotate_cross_codon_mnv_matches_oracle;
    cfg.type_info[0] = &kprop_cross_codon_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cross_mnv_cov, 0, sizeof g_cross_mnv_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cross_mnv_cov.missense > 0u);
    ASSERT(g_cross_mnv_cov.synonymous > 0u);
    ASSERT(g_cross_mnv_cov.stop_gained > 0u);
    ASSERT(g_cross_mnv_cov.fwd > 0u);
    ASSERT(g_cross_mnv_cov.rev > 0u);
    ASSERT(g_cross_mnv_cov.len2 > 0u);
    ASSERT(g_cross_mnv_cov.len3 > 0u);
    fprintf(stderr, "[cross-mnv coverage] missense=%u synonymous=%u stop_gained=%u fwd=%u rev=%u len2=%u len3=%u\n",
            g_cross_mnv_cov.missense, g_cross_mnv_cov.synonymous, g_cross_mnv_cov.stop_gained,
            g_cross_mnv_cov.fwd, g_cross_mnv_cov.rev,
            g_cross_mnv_cov.len2, g_cross_mnv_cov.len3);
    PASS();
}

static int kprop_cds_edits_equal(const duckvep_haplotype_edit_t *a,
                                  const duckvep_haplotype_edit_t *b) {
    return a != NULL && b != NULL &&
           a->cds_start == b->cds_start &&
           a->ref_len == b->ref_len &&
           a->ref == b->ref &&
           a->alt_len == b->alt_len &&
           a->alt == b->alt &&
           a->variant_strand == b->variant_strand;
}

static int kprop_context_zero(const duckvep_coding_context_t *ctx) {
    duckvep_coding_context_t z;
    memset(&z, 0, sizeof z);
    return ctx != NULL && memcmp(ctx, &z, sizeof z) == 0;
}

static char kprop_context_norm_base(uint8_t c) {
    switch ((char)c) {
        case 'A': case 'a': return 'A';
        case 'C': case 'c': return 'C';
        case 'G': case 'g': return 'G';
        case 'T': case 't': case 'U': case 'u': return 'T';
        case 'N': case 'n': return 'N';
        default: return '\0';
    }
}

int kprop_translate_full_oracle(const uint8_t *cds, size_t cds_len,
                                       duckvep_codon_table_t table,
                                       uint8_t *pep, size_t *pep_len) {
    size_t codons = cds_len / 3u;
    size_t i;
    if (cds == NULL || pep == NULL || pep_len == NULL) return 0;
    for (i = 0u; i < codons; i++) {
        char codon[4];
        uint32_t j;
        int has_n = 0;
        for (j = 0u; j < 3u; j++) {
            char b = kprop_context_norm_base(cds[i * 3u + (size_t)j]);
            if (b == '\0') return 0;
            if (b == 'N') has_n = 1;
            codon[j] = b;
        }
        codon[3] = '\0';
        pep[i] = has_n ? (uint8_t)'X'
                       : (uint8_t)duckvep_translate_codon(codon, table);
    }
    pep[codons] = (uint8_t)'\0';
    *pep_len = codons;
    return 1;
}

static int kprop_cds_changed_oracle(const uint8_t *ref_cds, size_t ref_len,
                                    const uint8_t *alt_cds, size_t alt_len) {
    size_t i;
    if (ref_len != alt_len) return 1;
    for (i = 0u; i < ref_len; i++) {
        char rb = kprop_context_norm_base(ref_cds[i]);
        char ab = kprop_context_norm_base(alt_cds[i]);
        if (rb == '\0' || ab == '\0') return -1;
        if (rb != ab) return 1;
    }
    return 0;
}

static void kprop_peptide_window_oracle(const uint8_t *ref_pep, size_t ref_len,
                                        const uint8_t *alt_pep, size_t alt_len,
                                        uint32_t *rf, uint32_t *rl,
                                        uint32_t *af, uint32_t *al) {
    size_t prefix = 0u;
    size_t suffix = 0u;
    size_t ref_mid;
    size_t alt_mid;
    *rf = 0u; *rl = 0u; *af = 0u; *al = 0u;
    while (prefix < ref_len && prefix < alt_len && ref_pep[prefix] == alt_pep[prefix]) {
        prefix++;
    }
    if (prefix == ref_len && prefix == alt_len) return;
    while (suffix < ref_len - prefix && suffix < alt_len - prefix &&
           ref_pep[ref_len - 1u - suffix] == alt_pep[alt_len - 1u - suffix]) {
        suffix++;
    }
    ref_mid = ref_len - prefix - suffix;
    alt_mid = alt_len - prefix - suffix;
    if (ref_mid > 0u) { *rf = (uint32_t)prefix + 1u; *rl = (uint32_t)(ref_len - suffix); }
    if (alt_mid > 0u) { *af = (uint32_t)prefix + 1u; *al = (uint32_t)(alt_len - suffix); }
}

static uint32_t kprop_context_flags_oracle(const duckvep_edit_set_t *edit_set) {
    uint32_t flags = 0u;
    int saw_frameshift = 0;
    int64_t total_diff = 0;
    size_t i;
    for (i = 0u; i < edit_set->count; i++) {
        int64_t d = (int64_t)edit_set->edits[i].alt_len -
                    (int64_t)edit_set->edits[i].ref_len;
        total_diff += d;
        if (d != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if ((d % 3) != 0) saw_frameshift = 1;
    }
    if (saw_frameshift) {
        if ((total_diff % 3) == 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
        else flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    }
    return flags;
}

static uint32_t kprop_single_variant_flags_oracle(uint32_t ref_len, uint32_t alt_len) {
    int64_t d = (int64_t)alt_len - (int64_t)ref_len;
    uint32_t flags = 0u;
    if (d != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
    if ((d % 3) != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    return flags;
}

TEST cds_edit_builder_projects_exon_boundary_insertions_on_both_strands(void) {
    static const uint16_t chrom[2] = {0u, 1u};
    static const uint32_t tx_start[2] = {100u, 100u};
    static const uint32_t tx_end[2] = {250u, 250u};
    static const int8_t strand[2] = {1, -1};
    static const uint32_t exon_offset[2] = {0u, 2u};
    static const uint16_t exon_count[2] = {2u, 2u};
    static const uint32_t cds_start[2] = {120u, 110u};
    static const uint32_t cds_end[2] = {240u, 230u};
    static const uint32_t exon_start[4] = {100u, 200u, 200u, 100u};
    static const uint32_t exon_end[4] = {150u, 250u, 250u, 150u};
    static const uint32_t cdna_start[4] = {1u, 52u, 1u, 52u};
    static const uint32_t cdna_end[4] = {51u, 102u, 51u, 102u};
    static const int8_t phase[4] = {0, 0, 0, 0};
    static const uint64_t cds_offset[2] = {0u, 72u};
    static const uint32_t cds_length[2] = {72u, 72u};
    static const uint8_t codon_table[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint32_t pos[2] = {199u, 199u};
    static const uint8_t kind[2] = {
        (uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t alleles[6] = {'G','G','T', 'G','G','T'};
    static const uint32_t ref_offset[2] = {0u, 3u};
    static const uint32_t alt_offset[2] = {1u, 4u};
    static const uint16_t ref_length[2] = {1u, 1u};
    static const uint16_t alt_length[2] = {2u, 2u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t variants;
    duckvep_haplotype_edit_t edit;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t results;
    duckvep_error_t error;
    uint8_t cds[144];
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    memset(cds, 'A', sizeof cds);
    cds[0] = cds[72] = 'A'; cds[1] = cds[73] = 'T'; cds[2] = cds[74] = 'G';
    cds[69] = cds[141] = 'T'; cds[70] = cds[142] = 'A';
    cds[71] = cds[143] = 'A';
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end; tx.strand = strand;
    tx.flags = k_zero_flags; tx.exon_offset = exon_offset; tx.exon_count = exon_count;
    tx.cds_start1 = cds_start; tx.cds_end1 = cds_end; tx.transcript_count = 2u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 4u;
    seq.cds_bytes = cds; seq.cds_bytes_len = sizeof cds;
    seq.cds_offset = cds_offset; seq.cds_length = cds_length;
    seq.codon_table = codon_table; seq.transcript_count = 2u;
    variants.chrom_id = chrom; variants.pos1 = pos; variants.end1 = pos;
    variants.variant_kind = kind; variants.allele_bytes = alleles;
    variants.allele_bytes_len = sizeof alleles; variants.ref_offset = ref_offset;
    variants.alt_offset = alt_offset; variants.ref_length = ref_length;
    variants.alt_length = alt_length; variants.count = 2u;

    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_build(&tx, &exons, &seq, &variants,
                                                 (uint32_t)i, i, strand[i], &edit));
        ASSERT_EQ(32u, edit.cds_start);
        ASSERT_EQ(0u, edit.ref_len);
        ASSERT_EQ(1u, edit.alt_len);
    }

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&results, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &results, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&results));
    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                  DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION),
                  rows[i].consequence_mask);
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[i].sequence_status);
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

TEST cds_edit_builder_checks_borrowed_sequence_extent(void) {
    uint8_t full_cds[12];
    memset(full_cds, 'A', sizeof full_cds);
    for (int strand = 1; strand >= -1; strand -= 2) {
        for (uint32_t length = 1u; length <= sizeof full_cds; length++) {
            /* A valid pool slice need not cover the borrowed model's mapped
             * CDS. Check both physical allocation ends and neighbouring bytes
             * that belong to another transcript in a larger sequence pool. */
            for (int padded = 0; padded < 2; padded++) {
                size_t bytes = padded ? sizeof full_cds : length;
                uint8_t *cds = malloc(bytes);
                ASSERT(cds != NULL);
                memset(cds, 'A', bytes);
                struct kprop_coding s = {0};
                s.cds = full_cds; s.strand = (int8_t)strand;
                s.tstart = s.cds_s = s.es = 1000u;
                s.tend = s.cds_e = s.ee = 1011u;
                s.ecds = 1u; s.ecde = 12u; s.excnt = 1u;
                kprop_wire_coding_scene(&s, length);
                s.seq.cds_bytes = cds; s.seq.cds_bytes_len = bytes;
                for (int insertion = 0; insertion < 2; insertion++) {
                    for (uint32_t width = 1u; width <= (insertion ? 1u : 3u); width++) {
                        for (uint32_t first = 1u; first + width - 1u <= 12u; first++) {
                            uint32_t coordinate = strand > 0 || insertion ? first : first + width - 1u;
                            s.vpos = strand > 0 ? s.es + coordinate - 1u : s.ee - coordinate + 1u;
                            s.vend = s.vpos + width - 1u;
                            s.vkind = insertion ? DUCKVEP_KIND_INS
                                : width == 1u ? DUCKVEP_KIND_SNV : DUCKVEP_KIND_MNV;
                            s.roff = 0u; s.rlen = (uint16_t)width;
                            s.aoff = width; s.alen = insertion ? 2u : (uint16_t)width;
                            memset(s.abytes, strand > 0 ? 'A' : 'T', width);
                            memset(s.abytes + width, 'C', s.alen);
                            if (insertion) s.abytes[width] = s.abytes[0];
                            duckvep_haplotype_edit_t edit;
                            duckvep_cds_edit_status_t status = duckvep_variant_cds_edit_build(
                                &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u, s.strand, &edit);
                            ASSERT_EQ(first + width - 1u <= length
                                ? DUCKVEP_CDS_EDIT_OK : DUCKVEP_CDS_EDIT_OUT_OF_CDS, status);
                            duckvep_event_t event = {0};
                            ASSERT(duckvep_event_prepare_small(s.vpos, s.abytes, s.rlen,
                                s.abytes + s.aoff, s.alen, &event));
                            duckvep_prepared_cds_allele_t allele = {&event,
                                s.abytes + event.ref_diff_offset,
                                s.abytes + s.aoff + event.alt_diff_offset,
                                s.abytes + event.anchor_ref_offset,
                                event.ref_diff_length, event.alt_diff_length, 1};
                            duckvep_haplotype_edit_t hinted;
                            ASSERT_EQ(status, duckvep_cds_edit_build_prepared_allele(
                                &s.tx, &s.ex, &s.seq, 0u, s.strand, &allele, 0u, &hinted));
                            if (status == DUCKVEP_CDS_EDIT_OK) {
                                ASSERT_EQ(first + (insertion && strand > 0), edit.cds_start);
                                ASSERT_EQ(insertion ? 0u : width, edit.ref_len);
                                ASSERT(kprop_cds_edits_equal(&edit, &hinted));
                            }
                        }
                    }
                }
                free(cds);
            }
        }
    }
    PASS();
}

TEST cds_edit_builder_known_scene(void) {
    static uint8_t plus_cds[12] = {'A','T','G','A','A','A','C','C','C','T','A','A'};
    static uint8_t minus_cds[12] = {'A','T','G','A','A','A','C','C','C','T','A','A'};
    struct kprop_coding s;
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_edit_t scratch[3];
    duckvep_edit_set_t edit_set;
    duckvep_haplotype_result_t result;
    uint8_t mutated[32];
    size_t mutated_len;
    uint8_t payload[2] = {'A','G'};

    memset(&s, 0, sizeof s);
    s.cds = plus_cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1011u; s.cds_s = 1000u; s.cds_e = 1011u;
    s.es = 1000u; s.ee = 1011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 12u);
    s.vpos = 1000u; s.vend = 1000u; s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s.abytes[0] = 'A'; s.abytes[1] = 'C';
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 1u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
                                             s.strand, &edit));
    ASSERT_EQ(1u, edit.cds_start);
    ASSERT_EQ(1u, edit.ref_len);
    ASSERT_EQ(1u, edit.alt_len);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, scratch, 1u,
                                                 &edit_set));
    ASSERT_EQ(1u, edit_set.count);
    ASSERT(edit_set.edits == scratch);
    ASSERT(kprop_cds_edits_equal(&edit, &scratch[0]));
    edit_set.edits = scratch;
    edit_set.count = 99u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, NULL, 0u,
                                                 &edit_set));
    ASSERT(edit_set.edits == NULL);
    ASSERT_EQ(0u, edit_set.count);

    memset(&s, 0, sizeof s);
    s.cds = plus_cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1011u; s.cds_s = 1000u; s.cds_e = 1011u;
    s.es = 1000u; s.ee = 1011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 12u);
    {
        uint8_t alt_mnv[5] = {'T','A','T','C','G'}; /* diff islands at CDS 4,6,8 */
        uint32_t i;
        s.vpos = kprop_genomic_pos_for_cds(&s, 4u); s.vend = s.vpos + 4u;
        s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
        for (i = 0u; i < 5u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 5u, alt_mnv, 5u);
        s.roff = 0u; s.aoff = 5u; s.rlen = 5u; s.alen = 5u;
        kprop_fill_expected_cds(&s, 4u, 5u, alt_mnv, 5u);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, &edit));
        ASSERT_EQ(5u, edit.ref_len);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 2u,
                                                     &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 3u,
                                                     &edit_set));
        ASSERT_EQ(3u, edit_set.count);
        ASSERT_EQ(8u, edit_set.edits[0].cds_start);
        ASSERT_EQ(6u, edit_set.edits[1].cds_start);
        ASSERT_EQ(4u, edit_set.edits[2].cds_start);
        ASSERT_EQ(1u, edit_set.edits[0].ref_len);
        ASSERT_EQ(1u, edit_set.edits[1].ref_len);
        ASSERT_EQ(1u, edit_set.edits[2].ref_len);
        ASSERT(edit_set.edits[0].ref == s.abytes + 4u);
        ASSERT(edit_set.edits[1].ref == s.abytes + 2u);
        ASSERT(edit_set.edits[2].ref == s.abytes);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_apply_cds_edits(s.cds, s.cds_lenv,
                                                    edit_set.edits, edit_set.count,
                                                    s.strand, mutated, sizeof mutated,
                                                    &mutated_len, &result));
        ASSERT_EQ((size_t)s.expect_len, mutated_len);
        ASSERT(memcmp(mutated, s.expect_cds, mutated_len) == 0);
        s.abytes[1] = kprop_base_not(s.abytes[1], (uint8_t)'N');
        s.abytes[6] = s.abytes[1]; /* retained internal mismatch must still validate */
        ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 3u,
                                                     &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
    }

    memset(&s, 0, sizeof s);
    s.cds = minus_cds; s.chrom = 1u; s.strand = -1; s.flags = 0u;
    s.tstart = 2000u; s.tend = 2011u; s.cds_s = 2000u; s.cds_e = 2011u;
    s.es = 2000u; s.ee = 2011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 1u;
    kprop_wire_coding_scene(&s, 12u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 5u); s.vend = s.vpos;
    s.vkind = (uint8_t)DUCKVEP_KIND_INS;
    s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
    s.abytes[1] = s.abytes[0];
    kprop_fill_variant_alt_from_tx(&s, 2u, payload, 2u);
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 3u;
    kprop_fill_expected_cds(&s, 5u, 0u, payload, 2u);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
                                             s.strand, &edit));
    ASSERT_EQ(5u, edit.cds_start);
    ASSERT_EQ(0u, edit.ref_len);
    ASSERT_EQ(2u, edit.alt_len);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, scratch, 1u,
                                                 &edit_set));
    ASSERT_EQ(1u, edit_set.count);
    ASSERT(kprop_cds_edits_equal(&edit, &scratch[0]));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits(s.cds, s.cds_lenv,
                                                edit_set.edits, edit_set.count,
                                                s.strand, mutated, sizeof mutated,
                                                &mutated_len, &result));
    ASSERT_EQ((size_t)s.expect_len, mutated_len);
    ASSERT(memcmp(mutated, s.expect_cds, mutated_len) == 0);
    s.abytes[0] = kprop_base_not(s.abytes[0], (uint8_t)'N');
    s.abytes[1] = s.abytes[0];
    ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH,
              duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
                                             s.strand, &edit));
    edit_set.edits = scratch;
    edit_set.count = 99u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, scratch, 1u,
                                                 &edit_set));
    ASSERT(edit_set.edits == NULL);
    ASSERT_EQ(0u, edit_set.count);

    memset(&s, 0, sizeof s);
    s.cds = minus_cds; s.chrom = 1u; s.strand = -1; s.flags = 0u;
    s.tstart = 2000u; s.tend = 2011u; s.cds_s = 2000u; s.cds_e = 2011u;
    s.es = 2000u; s.ee = 2011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 1u;
    kprop_wire_coding_scene(&s, 12u);
    {
        uint8_t alt_mnv[5] = {'T','A','T','C','G'}; /* transcript-order diff CDS 4,6,8 */
        uint32_t i;
        s.vpos = kprop_genomic_pos_for_cds(&s, 8u); s.vend = s.vpos + 4u;
        s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
        for (i = 0u; i < 5u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 5u, alt_mnv, 5u);
        s.roff = 0u; s.aoff = 5u; s.rlen = 5u; s.alen = 5u;
        kprop_fill_expected_cds(&s, 4u, 5u, alt_mnv, 5u);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 3u,
                                                     &edit_set));
        ASSERT_EQ(3u, edit_set.count);
        ASSERT_EQ(8u, edit_set.edits[0].cds_start);
        ASSERT_EQ(6u, edit_set.edits[1].cds_start);
        ASSERT_EQ(4u, edit_set.edits[2].cds_start);
        ASSERT(edit_set.edits[0].ref == s.abytes);
        ASSERT(edit_set.edits[1].ref == s.abytes + 2u);
        ASSERT(edit_set.edits[2].ref == s.abytes + 4u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_apply_cds_edits(s.cds, s.cds_lenv,
                                                    edit_set.edits, edit_set.count,
                                                    s.strand, mutated, sizeof mutated,
                                                    &mutated_len, &result));
        ASSERT_EQ((size_t)s.expect_len, mutated_len);
        ASSERT(memcmp(mutated, s.expect_cds, mutated_len) == 0);
    }

    {
        static const uint16_t tchrom[1] = {2u};
        static const uint32_t tstart[1] = {3000u};
        static const uint32_t tend[1] = {3014u};
        static const int8_t tstrand[1] = {1};
        static const uint64_t tflags[1] = {0u};
        static const uint32_t texoff[1] = {0u};
        static const uint16_t texcnt[1] = {2u};
        static const uint32_t tcds_s[1] = {3000u};
        static const uint32_t tcds_e[1] = {3014u};
        static const uint32_t estart[2] = {3000u, 3010u};
        static const uint32_t eend[2] = {3002u, 3014u};
        static const uint32_t ecdna_s[2] = {1u, 4u};
        static const uint32_t ecdna_e[2] = {3u, 8u};
        static const int8_t ephase[2] = {0, 0};
        static const uint8_t cds[8] = {'A','T','G','A','A','A','C','C'};
        static const uint64_t cds_off[1] = {0u};
        static const uint32_t cds_len[1] = {8u};
        static const uint8_t abytes[18] = {
            'G','A','A','A','A','A','A','A','T',
            'C','C','C','C','C','C','C','C','C'
        };
        static const uint16_t vchrom[1] = {2u};
        static const uint32_t vpos[1] = {3002u};
        static const uint32_t vend[1] = {3010u};
        static const uint32_t roff[1] = {0u};
        static const uint32_t aoff[1] = {9u};
        static const uint16_t rlen[1] = {9u};
        static const uint16_t alen[1] = {9u};
        static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_MNV};
        duckvep_transcript_model_t tx;
        duckvep_exon_model_t ex;
        duckvep_sequence_pool_t seq;
        duckvep_variant_batch_t v;
        memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
        memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v);
        tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
        tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
        tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
        ex.start1 = estart; ex.end1 = eend; ex.cdna_start1 = ecdna_s; ex.cdna_end1 = ecdna_e;
        ex.phase = ephase; ex.end_phase = ephase; ex.exon_count = 2u;
        seq.cds_bytes = cds; seq.cds_bytes_len = sizeof cds;
        seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.transcript_count = 1u;
        v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.ref_offset = roff; v.alt_offset = aoff;
        v.ref_length = rlen; v.alt_length = alen; v.allele_bytes = abytes;
        v.allele_bytes_len = sizeof abytes; v.variant_kind = vkind; v.count = 1u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                  duckvep_variant_cds_edit_build(&tx, &ex, &seq, &v, 0u, 0u, 1, &edit));
        edit_set.edits = scratch;
        edit_set.count = 99u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                  duckvep_variant_cds_edit_set_build(&tx, &ex, &seq, &v, 0u, 0u,
                                                     1, scratch, 1u, &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
    }

    {
        duckvep_haplotype_edit_t cap_scratch[1];
        edit_set.edits = cap_scratch;
        edit_set.count = 99u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_INVALID_ARG,
                  duckvep_variant_cds_edit_set_build(NULL, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand,
                                                     cap_scratch, 1u, &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
    }

    PASS();
}

TEST cds_edit_noncoding_without_sequence_pool(void) {
    uint32_t cds_start = 0u, cds_end = 0u, cds_length = 0u;
    uint64_t cds_offset = 0u;
    duckvep_transcript_model_t tx = {0};
    duckvep_exon_model_t exons = {0};
    duckvep_sequence_pool_t seq = {0};
    duckvep_event_t event;
    duckvep_haplotype_edit_t edit;
    tx.transcript_count = seq.transcript_count = 1u;
    tx.cds_start1 = &cds_start; tx.cds_end1 = &cds_end;
    seq.cds_offset = &cds_offset; seq.cds_length = &cds_length;
    ASSERT(duckvep_event_prepare_small(100u, (const uint8_t *)"A", 1u,
        (const uint8_t *)"C", 1u, &event));
    duckvep_prepared_cds_allele_t allele = {&event, (const uint8_t *)"A",
        (const uint8_t *)"C", (const uint8_t *)"A", 1u, 1u, 1};
    for (int strand = -1; strand <= 1; strand += 2) {
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS,
            duckvep_cds_edit_build_prepared_allele(&tx, &exons, &seq, 0u,
                (int8_t)strand, &allele, UINT32_MAX, &edit));
        ASSERT_EQ(0u, edit.cds_start);
        ASSERT(edit.ref == NULL && edit.alt == NULL);
    }
    cds_start = 100u; cds_end = 111u; cds_length = 12u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_INVALID_ARG,
        duckvep_cds_edit_build_prepared_allele(&tx, &exons, &seq, 0u,
            1, &allele, UINT32_MAX, &edit));
    PASS();
}

TEST coding_context_open_replay_borrows_complete_sequences(void) {
    static const uint8_t ref[] = "ATGAAATAAGCC", alt[] = "ATGAAGTAAGTT";
    static const uint8_t rp[] = "MK*A", ap[] = "MK*V";
    const duckvep_haplotype_edit_t edits[] = {
        {10u, 3u, ref + 9u, 3u, alt + 9u, 1},
        {4u, 3u, ref + 3u, 3u, alt + 3u, 1}
    };
    duckvep_edit_set_t set = {edits, 2u};
    duckvep_haplotype_result_t applied = {12u, 0, 0u, 2u};
    duckvep_translation_t rt = {4u, 3u, 1u}, at = rt;
    duckvep_coding_context_t ctx;
#define OPEN_REPLAY() duckvep_coding_context_open_replay(ref, 12u, &set, 1, \
    DUCKVEP_CODON_TABLE_STANDARD, alt, &applied, rp, &rt, ap, &at, &ctx)
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, OPEN_REPLAY());
    ASSERT(ctx.ref_cds == ref && ctx.alt_cds == alt);
    ASSERT(ctx.ref_peptide == rp && ctx.alt_peptide == ap);
    ASSERT_EQ(4u, ctx.alt_peptide_len);
    ASSERT_EQ(3u, ctx.alt_first_stop_position1);
    ASSERT_EQ(4u, ctx.ref_first_changed_codon);
    ASSERT_EQ(4u, ctx.alt_last_changed_codon);
    ASSERT_EQ(2u, ctx.applied_edits);
    ASSERT(ctx.cds_changed && !ctx.has_single_edit);
    for (unsigned mutation = 0u; mutation < 8u; mutation++) {
        applied = (duckvep_haplotype_result_t){12u, 0, 0u, 2u};
        rt = at = (duckvep_translation_t){4u, 3u, 1u};
        switch (mutation) {
        case 0: applied.applied_edits = 1u; break;
        case 1: applied.length_diff = 1; break;
        case 2: applied.cds_len = 13u; break;
        case 3: rt.length = 3u; break;
        case 4: at.length = 3u; break; /* A displayed first-stop prefix is not complete. */
        case 5: rt.first_stop_position1 = 5u; break;
        case 6: at.first_stop_position1 = 5u; break;
        case 7: at.first_stop_position1 = 2u; break;
        }
        memset(&ctx, 0xA5, sizeof ctx);
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_INVALID_ARG, OPEN_REPLAY());
        ASSERT(kprop_context_zero(&ctx));
    }
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_INVALID_ARG,
        duckvep_coding_context_open_replay(ref, 12u, &set, 1,
            DUCKVEP_CODON_TABLE_STANDARD, alt, NULL, rp, &rt, ap, &at, &ctx));
    ASSERT(kprop_context_zero(&ctx));
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_INVALID_ARG,
        duckvep_coding_context_open_replay(ref, 12u, &set, 1,
            DUCKVEP_CODON_TABLE_STANDARD, alt, &applied, rp, NULL, ap, &at, &ctx));
    ASSERT(kprop_context_zero(&ctx));
#undef OPEN_REPLAY
    PASS();
}

TEST coding_context_known_scene(void) {
    {
        static const uint8_t ref_cds[9] = {'A','T','G','T','A','A','G','A','A'};
        static const uint8_t ref_base[1] = {'A'};
        static const uint8_t alt_base[1] = {'C'};
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        duckvep_translation_t hres;
        uint8_t alt_cds[16];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];
        size_t trunc_len = 0u;

        edit.cds_start = 8u; edit.ref_len = 1u; edit.ref = ref_base;
        edit.alt_len = 1u; edit.alt = alt_base; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        memset(ref_pep, 0xA5, sizeof ref_pep);
        memset(alt_pep, 0x5A, sizeof alt_pep);
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ((size_t)9u, ctx.ref_cds_len);
        ASSERT_EQ((size_t)9u, ctx.alt_cds_len);
        ASSERT_EQ((size_t)3u, ctx.ref_peptide_len);
        ASSERT_EQ((size_t)3u, ctx.alt_peptide_len);
        ASSERT(memcmp(ctx.ref_peptide, "M*E", 3u) == 0);
        ASSERT(memcmp(ctx.alt_peptide, "M*A", 3u) == 0);
        ASSERT_EQ((uint8_t)'\0', ctx.ref_peptide[ctx.ref_peptide_len]);
        ASSERT_EQ((uint8_t)'\0', ctx.alt_peptide[ctx.alt_peptide_len]);
        ASSERT_EQ(1u, ctx.cds_changed);
        ASSERT_EQ((size_t)1u, ctx.applied_edits);
        ASSERT_EQ(0, ctx.length_diff);
        ASSERT_EQ(3u, ctx.ref_first_changed_codon);
        ASSERT_EQ(3u, ctx.ref_last_changed_codon);
        ASSERT_EQ(3u, ctx.alt_first_changed_codon);
        ASSERT_EQ(3u, ctx.alt_last_changed_codon);
        ASSERT_EQ(DUCKVEP_TRANSLATION_OK,
                  duckvep_translate_cds(ctx.alt_cds, ctx.alt_cds_len,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_pep, sizeof alt_pep, &hres));
        trunc_len = hres.first_stop_position1 ? hres.first_stop_position1 : hres.length;
        ASSERT_EQ((size_t)2u, trunc_len);
        ASSERT(memcmp(alt_pep, "M*", 2u) == 0);
        ASSERT(trunc_len < hres.length);

        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, 8u, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT(kprop_context_zero(&ctx));
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, 3u,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT(kprop_context_zero(&ctx));
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep, alt_pep, 3u,
                                               &ctx));
        ASSERT(kprop_context_zero(&ctx));
    }

    {
        static const uint8_t ref_cds[6] = {'G','C','T','G','A','A'};
        static const uint8_t ref_base[1] = {'T'};
        static const uint8_t alt_base[1] = {'C'};
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        uint8_t alt_cds[8];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];

        edit.cds_start = 3u; edit.ref_len = 1u; edit.ref = ref_base;
        edit.alt_len = 1u; edit.alt = alt_base; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(1u, ctx.cds_changed);
        ASSERT_EQ((size_t)2u, ctx.ref_peptide_len);
        ASSERT_EQ((size_t)2u, ctx.alt_peptide_len);
        ASSERT(memcmp(ctx.ref_peptide, ctx.alt_peptide, 2u) == 0);
        ASSERT_EQ(0u, ctx.ref_first_changed_codon);
        ASSERT_EQ(0u, ctx.ref_last_changed_codon);
        ASSERT_EQ(0u, ctx.alt_first_changed_codon);
        ASSERT_EQ(0u, ctx.alt_last_changed_codon);
    }

    {
        static const uint8_t ref_cds[8] = {'A','T','G','G','A','A','C','T'};
        static const uint8_t ref_base[1] = {'T'};
        static const uint8_t alt_base[1] = {'A'};
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        uint8_t alt_cds[8];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];

        edit.cds_start = 8u; edit.ref_len = 1u; edit.ref = ref_base;
        edit.alt_len = 1u; edit.alt = alt_base; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ((size_t)8u, ctx.alt_cds_len);
        ASSERT_EQ((size_t)2u, ctx.ref_peptide_len);
        ASSERT_EQ((size_t)2u, ctx.alt_peptide_len);
        ASSERT(memcmp(ctx.ref_peptide, "ME", 2u) == 0);
        ASSERT(memcmp(ctx.alt_peptide, "ME", 2u) == 0);
        ASSERT_EQ(1u, ctx.cds_changed);
        ASSERT_EQ(0u, ctx.ref_first_changed_codon);
        ASSERT_EQ(0u, ctx.alt_first_changed_codon);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','N','G','A','A'};
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        uint8_t alt_cds[8];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];
        edit_set.edits = NULL; edit_set.count = 0u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(0u, ctx.cds_changed);
        ASSERT(memcmp(ctx.ref_peptide, "XE", 2u) == 0);
        ASSERT(memcmp(ctx.alt_peptide, "XE", 2u) == 0);
    }

    PASS();
}

TEST variant_coding_context_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'T','A','A',  'G','A','A',  'C','C','C',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t edit_scratch[3];
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[24];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    uint8_t alt_mnv[5] = {'A','A','C','C','G'};
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 7u); s.vend = s.vpos + 4u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 5u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 5u, alt_mnv, 5u);
    s.roff = 0u; s.aoff = 5u; s.rlen = 5u; s.alen = 5u;
    kprop_fill_expected_cds(&s, 7u, 5u, alt_mnv, 5u);

    ctx.ref_cds = cds;
    ctx.ref_cds_len = 99u;
    ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL,
              duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                   0u, 0u, s.strand,
                                                   edit_scratch, 2u,
                                                   alt_cds, sizeof alt_cds,
                                                   ref_pep, sizeof ref_pep,
                                                   alt_pep, sizeof alt_pep,
                                                   &ctx));
    ASSERT(kprop_context_zero(&ctx));

    memset(ref_pep, 0xA5, sizeof ref_pep);
    memset(alt_pep, 0x5A, sizeof alt_pep);
    ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
              duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                   0u, 0u, s.strand,
                                                   edit_scratch, 3u,
                                                   alt_cds, sizeof alt_cds,
                                                   ref_pep, sizeof ref_pep,
                                                   alt_pep, sizeof alt_pep,
                                                   &ctx));
    ASSERT(ctx.ref_cds == s.cds);
    ASSERT(ctx.alt_cds == alt_cds);
    ASSERT_EQ((size_t)15u, ctx.ref_cds_len);
    ASSERT_EQ((size_t)15u, ctx.alt_cds_len);
    ASSERT(memcmp(ctx.alt_cds, s.expect_cds, (size_t)s.expect_len) == 0);
    ASSERT_EQ((size_t)5u, ctx.ref_peptide_len);
    ASSERT_EQ((size_t)5u, ctx.alt_peptide_len);
    ASSERT(memcmp(ctx.ref_peptide, "M*EPF", 5u) == 0);
    ASSERT(memcmp(ctx.alt_peptide, "M*NRF", 5u) == 0);
    ASSERT_EQ((uint8_t)'\0', ctx.ref_peptide[ctx.ref_peptide_len]);
    ASSERT_EQ((uint8_t)'\0', ctx.alt_peptide[ctx.alt_peptide_len]);
    ASSERT_EQ((size_t)3u, ctx.applied_edits);
    ASSERT_EQ(1u, ctx.cds_changed);
    ASSERT_EQ(3u, ctx.ref_first_changed_codon);
    ASSERT_EQ(4u, ctx.ref_last_changed_codon);
    ASSERT_EQ(3u, ctx.alt_first_changed_codon);
    ASSERT_EQ(4u, ctx.alt_last_changed_codon);

    {
        static uint8_t mito_cds[6] = {'A','T','A','G','A','A'};
        uint8_t ref_base[1] = {'A'};
        uint8_t alt_base[1] = {'G'};
        memset(&s, 0, sizeof s);
        s.cds = mito_cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
        s.tstart = 2000u; s.tend = 2005u; s.cds_s = 2000u; s.cds_e = 2005u;
        s.es = 2000u; s.ee = 2005u; s.ecds = 1u; s.ecde = 6u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 6u);
        s.ctab = (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO;
        s.vpos = kprop_genomic_pos_for_cds(&s, 6u); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
        s.abytes[0] = ref_base[0]; s.abytes[1] = alt_base[0];
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 1u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edit_scratch, 1u,
                                                       alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(memcmp(ctx.ref_peptide, "ME", 2u) == 0);
    }

    PASS();
}

TEST coding_context_delta_known_scene(void) {
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    {
        static const uint8_t ref_cds[6] = {'A','T','G','G','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','G','A','G'};
        static const uint8_t pep[3] = {'M','E','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.synonymous && !delta.start_lost);
        ASSERT_EQ(2, delta.protein_pos);
        ASSERT_EQ((uint8_t)'E', delta.ref_aa);
        ASSERT_EQ((uint8_t)'E', delta.alt_aa);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','G','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','G','A','C'};
        static const uint8_t ref_pep[3] = {'M','E','\0'};
        static const uint8_t alt_pep[3] = {'M','D','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.missense && !delta.start_lost);
        ASSERT_EQ((uint8_t)'E', delta.ref_aa);
        ASSERT_EQ((uint8_t)'D', delta.alt_aa);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','T','G','G'};
        static const uint8_t alt_cds[6] = {'A','T','G','T','G','A'};
        static const uint8_t ref_pep[3] = {'M','W','\0'};
        static const uint8_t alt_pep[3] = {'M','*','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && !delta.start_lost);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','T','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','C','A','A'};
        static const uint8_t ref_pep[3] = {'M','*','\0'};
        static const uint8_t alt_pep[3] = {'M','Q','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_lost && !delta.start_lost);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','T','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','T','A','G'};
        static const uint8_t pep[3] = {'M','*','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained && !delta.start_lost);
    }

    {
        static const uint8_t ref_cds[3] = {'A','T','G'};
        static const uint8_t alt_cds[3] = {'G','T','G'};
        static const uint8_t ref_pep[2] = {'M','\0'};
        static const uint8_t alt_pep[2] = {'V','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_lost && !delta.missense);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx,
                                                    (uint64_t)DUCKVEP_TX_CDS_START_NF,
                                                    &delta));
        ASSERT(delta.valid && !delta.start_lost && delta.missense);
    }

    {
        /* Ensembl annotates supported non-ATG starts. VEP compares the
         * alternate peptide with that annotated reference peptide; it does
         * not assume that every complete reference start translates to M. */
        static const uint8_t ref_cds[3] = {'C','T','G'};
        static const uint8_t alt_cds[3] = {'A','T','G'};
        static const uint8_t ref_pep[2] = {'L','\0'};
        static const uint8_t alt_pep[2] = {'M','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_lost && !delta.missense);
    }

    {
        /* Two-codon body substitution (codons 2-3, past the start): the generalized window
         * classifier resolves it as missense_variant — EE -> DD — where the old two-codon
         * slice bailed to coding_sequence_variant. Coarse window: protein_pos -1, no AA pair. */
        static const uint8_t ref_cds[9] = {'A','T','G','G','A','A','G','A','A'};
        static const uint8_t alt_cds[9] = {'A','T','G','G','A','C','G','A','C'};
        static const uint8_t ref_pep[4] = {'M','E','E','\0'};
        static const uint8_t alt_pep[4] = {'M','D','D','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 9u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 9u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 3u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 3u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.missense && !delta.synonymous && !delta.stop_gained &&
               !delta.stop_lost && !delta.stop_retained && !delta.start_lost);
        ASSERT(delta.protein_pos == -1 && delta.ref_aa == (uint8_t)0u &&
               delta.alt_aa == (uint8_t)0u);
    }

    {
        static const uint8_t ref_cds[3] = {'G','A','N'};
        static const uint8_t alt_cds[3] = {'G','A','C'};
        static const uint8_t ref_pep[2] = {'X','\0'};
        static const uint8_t alt_pep[2] = {'D','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
        /* Consensus operands require an explicit table and matching raw
         * peptides. Source-allele validation belongs to the projection layer. */
        ctx.codon_table = DUCKVEP_CODON_TABLE_STANDARD;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_lost && !delta.missense);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK, duckvep_coding_context_delta_fill(
            &ctx, (uint64_t)DUCKVEP_TX_CDS_START_NF, &delta));
        ASSERT(delta.valid && !delta.start_lost && delta.missense && delta.coding_unknown);
        /* An invented reference peptide must still fail with a valid table. */
        {
            static const uint8_t fake_ref_pep[2] = {'E','\0'};
            memset(&ctx, 0, sizeof ctx);
            ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
            ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
            ctx.ref_peptide = fake_ref_pep; ctx.ref_peptide_len = 1u;
            ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
            ctx.codon_table = DUCKVEP_CODON_TABLE_STANDARD;
            ctx.cds_changed = 1u;
        }
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    {
        static const uint8_t ref_cds[3] = {'G','A','A'};
        static const uint8_t alt_cds[4] = {'G','A','A','A'};
        static const uint8_t pep[2] = {'E','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 4u;
        ctx.ref_peptide = pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ctx.length_diff = 1;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
} g_cds_edit_builder_cov;

static enum theft_trial_res prop_cds_edit_builder_matches_splice_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_result_t result;
    duckvep_event_t event;
    uint32_t projected_cds_start = 0u;
    uint32_t projected_cds_end = 0u;
    uint8_t mutated[64];
    size_t mutated_len = 0u;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;

    {
        duckvep_cds_edit_status_t st = duckvep_variant_cds_edit_build(
            &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand, &edit);
        if (st != DUCKVEP_CDS_EDIT_OK) {
            fprintf(stderr, "\n[cds-edit-builder fail] build status=%d shape=%u strand=%d pos=%u end=%u ref_len=%u alt_len=%u\n",
                    (int)st, (unsigned)s->expect_shape, (int)s->strand,
                    s->vpos, s->vend, (unsigned)s->rlen, (unsigned)s->alen);
            return THEFT_TRIAL_FAIL;
        }
    }
    if (!duckvep_event_prepare_small(
            s->vpos, s->v.allele_bytes + s->roff, s->rlen,
            s->v.allele_bytes + s->aoff, s->alen, &event) ||
        !duckvep_project_event_to_cds(&s->tx, &s->ex, 0u, &event,
                                      &projected_cds_start,
                                      &projected_cds_end) ||
        projected_cds_start != edit.cds_start ||
        projected_cds_end != edit.cds_start +
            (edit.ref_len == 0u ? 0u : edit.ref_len - 1u)) {
        fprintf(stderr,
                "\n[cds-edit-builder fail] shared event projection disagrees shape=%u strand=%d edit=%u+%u projected=%u-%u\n",
                (unsigned)s->expect_shape, (int)s->strand,
                edit.cds_start, edit.ref_len, projected_cds_start,
                projected_cds_end);
        return THEFT_TRIAL_FAIL;
    }
    {
        duckvep_haplotype_status_t hst = duckvep_haplotype_apply_cds_edits(
            s->cds, s->cds_lenv, &edit, 1u, s->strand, mutated, sizeof mutated,
            &mutated_len, &result);
        if (hst != DUCKVEP_HAPLOTYPE_OK) {
            fprintf(stderr, "\n[cds-edit-builder fail] haplo status=%d shape=%u strand=%d cds_start=%u ref_len=%u alt_len=%u\n",
                    (int)hst, (unsigned)s->expect_shape, (int)s->strand,
                    edit.cds_start, edit.ref_len, edit.alt_len);
            return THEFT_TRIAL_FAIL;
        }
    }
    if (mutated_len != (size_t)s->expect_len ||
        memcmp(mutated, s->expect_cds, (size_t)s->expect_len) != 0) {
        fprintf(stderr,
                "\n[cds-edit-builder fail] mismatch shape=%u strand=%d vpos=%u end=%u edit_start=%u ref_len=%u alt_len=%u got_len=%zu want_len=%u\n",
                (unsigned)s->expect_shape, (int)s->strand, s->vpos, s->vend,
                edit.cds_start, edit.ref_len, edit.alt_len, mutated_len,
                (unsigned)s->expect_len);
        tr = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_cds_edit_builder_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_cds_edit_builder_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_cds_edit_builder_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_cds_edit_builder_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_cds_edit_builder_cov.indel++;
    else tr = THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_cds_edit_builder_cov.fwd++; else g_cds_edit_builder_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START) g_cds_edit_builder_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY) g_cds_edit_builder_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP) g_cds_edit_builder_cov.stop++;
    else tr = THEFT_TRIAL_FAIL;
done:
    return tr;
}

TEST cds_edit_builder_matches_direct_splice_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant CDS edit builder == direct CDS splice oracle";
    cfg.prop1 = prop_cds_edit_builder_matches_splice_oracle;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cds_edit_builder_cov, 0, sizeof g_cds_edit_builder_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cds_edit_builder_cov.snv > 0u);
    ASSERT(g_cds_edit_builder_cov.mnv > 0u);
    ASSERT(g_cds_edit_builder_cov.ins > 0u);
    ASSERT(g_cds_edit_builder_cov.del > 0u);
    ASSERT(g_cds_edit_builder_cov.indel > 0u);
    ASSERT(g_cds_edit_builder_cov.fwd > 0u);
    ASSERT(g_cds_edit_builder_cov.rev > 0u);
    ASSERT(g_cds_edit_builder_cov.start > 0u);
    ASSERT(g_cds_edit_builder_cov.body > 0u);
    ASSERT(g_cds_edit_builder_cov.stop > 0u);
    fprintf(stderr,
            "[cds-edit-builder coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u start=%u body=%u stop=%u\n",
            g_cds_edit_builder_cov.snv, g_cds_edit_builder_cov.mnv,
            g_cds_edit_builder_cov.ins, g_cds_edit_builder_cov.del,
            g_cds_edit_builder_cov.indel, g_cds_edit_builder_cov.fwd,
            g_cds_edit_builder_cov.rev, g_cds_edit_builder_cov.start,
            g_cds_edit_builder_cov.body, g_cds_edit_builder_cov.stop);
    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
    uint32_t cap0;
} g_cds_edit_set_cov;

static enum theft_trial_res prop_cds_edit_set_matches_single_builder(struct theft *t,
                                                                     void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_edit_t scratch[1];
    duckvep_edit_set_t edit_set;
    duckvep_haplotype_result_t result;
    uint8_t mutated[64];
    size_t mutated_len = 0u;
    duckvep_cds_edit_status_t single_st;
    duckvep_cds_edit_status_t set_st;
    (void)t;

    single_st = duckvep_variant_cds_edit_build(&s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->strand, &edit);
    edit_set.edits = scratch;
    edit_set.count = 99u;
    set_st = duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                0u, 0u, s->strand, scratch, 1u,
                                                &edit_set);
    if (single_st != set_st) return THEFT_TRIAL_FAIL;
    if (set_st != DUCKVEP_CDS_EDIT_OK) {
        return edit_set.edits == NULL && edit_set.count == 0u
            ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
    }
    if (edit_set.edits != scratch || edit_set.count != 1u ||
        !kprop_cds_edits_equal(&edit, &scratch[0])) {
        return THEFT_TRIAL_FAIL;
    }
    {
        duckvep_edit_set_t too_small;
        too_small.edits = scratch;
        too_small.count = 99u;
        if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->strand, NULL, 0u,
                                               &too_small) != DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (too_small.edits != NULL || too_small.count != 0u) return THEFT_TRIAL_FAIL;
        g_cds_edit_set_cov.cap0++;
    }
    if (duckvep_haplotype_apply_cds_edits(s->cds, s->cds_lenv, edit_set.edits,
                                          edit_set.count, s->strand, mutated,
                                          sizeof mutated, &mutated_len,
                                          &result) != DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (mutated_len != (size_t)s->expect_len ||
        memcmp(mutated, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_cds_edit_set_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_cds_edit_set_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_cds_edit_set_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_cds_edit_set_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_cds_edit_set_cov.indel++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_cds_edit_set_cov.fwd++; else g_cds_edit_set_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START) g_cds_edit_set_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY) g_cds_edit_set_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP) g_cds_edit_set_cov.stop++;
    else return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST cds_edit_set_builder_matches_single_edit_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant CDS edit-set builder == single-edit splice oracle";
    cfg.prop1 = prop_cds_edit_set_matches_single_builder;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cds_edit_set_cov, 0, sizeof g_cds_edit_set_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cds_edit_set_cov.snv > 0u);
    ASSERT(g_cds_edit_set_cov.mnv > 0u);
    ASSERT(g_cds_edit_set_cov.ins > 0u);
    ASSERT(g_cds_edit_set_cov.del > 0u);
    ASSERT(g_cds_edit_set_cov.indel > 0u);
    ASSERT(g_cds_edit_set_cov.fwd > 0u);
    ASSERT(g_cds_edit_set_cov.rev > 0u);
    ASSERT(g_cds_edit_set_cov.start > 0u);
    ASSERT(g_cds_edit_set_cov.body > 0u);
    ASSERT(g_cds_edit_set_cov.stop > 0u);
    ASSERT(g_cds_edit_set_cov.cap0 > 0u);
    fprintf(stderr,
            "[cds-edit-set coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u start=%u body=%u stop=%u cap0=%u\n",
            g_cds_edit_set_cov.snv, g_cds_edit_set_cov.mnv,
            g_cds_edit_set_cov.ins, g_cds_edit_set_cov.del,
            g_cds_edit_set_cov.indel, g_cds_edit_set_cov.fwd,
            g_cds_edit_set_cov.rev, g_cds_edit_set_cov.start,
            g_cds_edit_set_cov.body, g_cds_edit_set_cov.stop,
            g_cds_edit_set_cov.cap0);
    PASS();
}

static struct {
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
    uint32_t multi;
    uint32_t capfail;
} g_cds_edit_set_mnv_cov;

static enum theft_trial_res prop_cds_edit_set_splits_mnv_islands(struct theft *t,
                                                                 void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t scratch[8];
    duckvep_edit_set_t edit_set;
    duckvep_haplotype_result_t result;
    uint8_t mutated[64];
    uint8_t covered[64];
    uint32_t exp_start[64];
    uint32_t exp_len[64];
    size_t exp_n = 0u;
    size_t mutated_len = 0u;
    size_t i;
    (void)t;

    memset(covered, 0, sizeof covered);
    if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                           0u, 0u, s->strand, scratch, 8u,
                                           &edit_set) != DUCKVEP_CDS_EDIT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (edit_set.count < 2u) return THEFT_TRIAL_FAIL;
    i = 0u;
    while (i < s->cds_lenv) {
        uint32_t start;
        uint32_t len;
        while (i < s->cds_lenv && s->cds[i] == s->expect_cds[i]) i++;
        if (i >= s->cds_lenv) break;
        start = (uint32_t)i;
        len = 0u;
        while (i < s->cds_lenv && s->cds[i] != s->expect_cds[i]) {
            i++;
            len++;
        }
        exp_start[exp_n] = start;
        exp_len[exp_n] = len;
        exp_n++;
    }
    if (edit_set.count != exp_n) return THEFT_TRIAL_FAIL;
    for (i = 0u; i < edit_set.count; i++) {
        uint32_t j;
        size_t exp_i = exp_n - 1u - i;
        if (edit_set.edits[i].cds_start != exp_start[exp_i] + 1u ||
            edit_set.edits[i].ref_len != exp_len[exp_i] ||
            edit_set.edits[i].alt_len != exp_len[exp_i]) {
            return THEFT_TRIAL_FAIL;
        }
        if (i + 1u < edit_set.count &&
            edit_set.edits[i].cds_start <= edit_set.edits[i + 1u].cds_start) {
            return THEFT_TRIAL_FAIL;
        }
        for (j = 0u; j < edit_set.edits[i].ref_len; j++) {
            uint32_t cds_pos = edit_set.edits[i].cds_start + j;
            if (cds_pos == 0u || cds_pos > s->cds_lenv) return THEFT_TRIAL_FAIL;
            covered[cds_pos - 1u] = 1u;
        }
    }
    for (i = 0u; i < s->cds_lenv; i++) {
        int diff = s->cds[i] != s->expect_cds[i];
        if ((covered[i] != 0u) != diff) return THEFT_TRIAL_FAIL;
    }
    if (edit_set.count > 1u) {
        duckvep_edit_set_t fail_set;
        fail_set.edits = scratch;
        fail_set.count = 99u;
        if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->strand, scratch,
                                               edit_set.count - 1u,
                                               &fail_set) != DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (fail_set.edits != NULL || fail_set.count != 0u) return THEFT_TRIAL_FAIL;
        g_cds_edit_set_mnv_cov.capfail++;
    }
    if (duckvep_haplotype_apply_cds_edits(s->cds, s->cds_lenv, edit_set.edits,
                                          edit_set.count, s->strand, mutated,
                                          sizeof mutated, &mutated_len,
                                          &result) != DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (mutated_len != (size_t)s->expect_len ||
        memcmp(mutated, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_cds_edit_set_mnv_cov.fwd++; else g_cds_edit_set_mnv_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START) g_cds_edit_set_mnv_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY) g_cds_edit_set_mnv_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP) g_cds_edit_set_mnv_cov.stop++;
    else return THEFT_TRIAL_FAIL;
    g_cds_edit_set_mnv_cov.multi++;
    return THEFT_TRIAL_PASS;
}

TEST cds_edit_set_builder_splits_mnv_diff_islands(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant CDS edit-set builder splits MNV diff islands";
    cfg.prop1 = prop_cds_edit_set_splits_mnv_islands;
    cfg.type_info[0] = &kprop_cds_edit_set_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cds_edit_set_mnv_cov, 0, sizeof g_cds_edit_set_mnv_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cds_edit_set_mnv_cov.fwd > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.rev > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.start > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.body > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.stop > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.multi > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.capfail > 0u);
    fprintf(stderr,
            "[cds-edit-set-mnv coverage] fwd=%u rev=%u start=%u body=%u stop=%u multi=%u capfail=%u\n",
            g_cds_edit_set_mnv_cov.fwd, g_cds_edit_set_mnv_cov.rev,
            g_cds_edit_set_mnv_cov.start, g_cds_edit_set_mnv_cov.body,
            g_cds_edit_set_mnv_cov.stop, g_cds_edit_set_mnv_cov.multi,
            g_cds_edit_set_mnv_cov.capfail);
    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t pep_same;
    uint32_t pep_diff;
    uint32_t capfail;
} g_coding_context_cov;

static enum theft_trial_res prop_coding_context_matches_direct_oracles(struct theft *t,
                                                                       void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t scratch[8];
    duckvep_edit_set_t edit_set;
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint8_t ref_want[32];
    uint8_t alt_want[32];
    size_t ref_want_len = 0u;
    size_t alt_want_len = 0u;
    int cds_changed;
    uint32_t rf, rl, af, al;
    (void)t;

    if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                           0u, 0u, s->strand, scratch, 8u,
                                           &edit_set) != DUCKVEP_CDS_EDIT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    memset(ref_pep, 0xA5, sizeof ref_pep);
    memset(alt_pep, 0x5A, sizeof alt_pep);
    if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                     (duckvep_codon_table_t)s->ctab,
                                     alt_cds, sizeof alt_cds,
                                     ref_pep, sizeof ref_pep,
                                     alt_pep, sizeof alt_pep,
                                     &ctx) != DUCKVEP_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.alt_cds != alt_cds || ctx.ref_cds != s->cds) return THEFT_TRIAL_FAIL;
    if (ctx.alt_cds_len != (size_t)s->expect_len ||
        memcmp(ctx.alt_cds, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->cds, s->cds_lenv,
                                     (duckvep_codon_table_t)s->ctab,
                                     ref_want, &ref_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->expect_cds, (size_t)s->expect_len,
                                     (duckvep_codon_table_t)s->ctab,
                                     alt_want, &alt_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide_len != ref_want_len || ctx.alt_peptide_len != alt_want_len) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide[ctx.ref_peptide_len] != (uint8_t)'\0' ||
        ctx.alt_peptide[ctx.alt_peptide_len] != (uint8_t)'\0') {
        return THEFT_TRIAL_FAIL;
    }
    if (memcmp(ctx.ref_peptide, ref_want, ref_want_len) != 0 ||
        memcmp(ctx.alt_peptide, alt_want, alt_want_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    cds_changed = kprop_cds_changed_oracle(s->cds, s->cds_lenv,
                                           s->expect_cds, (size_t)s->expect_len);
    if (cds_changed < 0 || ctx.cds_changed != (uint8_t)cds_changed) return THEFT_TRIAL_FAIL;
    if (ctx.applied_edits != edit_set.count) return THEFT_TRIAL_FAIL;
    if (ctx.length_diff != (int64_t)s->expect_len - (int64_t)s->cds_lenv) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.flags != kprop_context_flags_oracle(&edit_set)) return THEFT_TRIAL_FAIL;
    kprop_peptide_window_oracle(ref_want, ref_want_len, alt_want, alt_want_len,
                                &rf, &rl, &af, &al);
    if (ctx.ref_first_changed_codon != rf || ctx.ref_last_changed_codon != rl ||
        ctx.alt_first_changed_codon != af || ctx.alt_last_changed_codon != al) {
        return THEFT_TRIAL_FAIL;
    }
    /* The borrowed path consumes the independent oracle's already completed
     * sequences, not the builder's buffers or another replay. Keep every
     * original biological/capacity check and its generator unchanged. */
    duckvep_haplotype_result_t applied = {(size_t)s->expect_len,
        (int64_t)s->expect_len - (int64_t)s->cds_lenv,
        kprop_context_flags_oracle(&edit_set), edit_set.count};
    duckvep_translation_t rt = {ref_want_len, 0u, 1u}, at = {alt_want_len, 0u, 1u};
    for (size_t i = 0u; i < ref_want_len; i++)
        if (ref_want[i] == '*') { rt.first_stop_position1 = i + 1u; break; }
    for (size_t i = 0u; i < alt_want_len; i++)
        if (alt_want[i] == '*') { at.first_stop_position1 = i + 1u; break; }
    duckvep_coding_context_t borrowed, expected = ctx;
    expected.alt_cds = s->expect_cds;
    expected.ref_peptide = ref_want;
    expected.alt_peptide = alt_want;
    if (duckvep_coding_context_open_replay(s->cds, s->cds_lenv, &edit_set, s->strand,
            (duckvep_codon_table_t)s->ctab, s->expect_cds, &applied,
            ref_want, &rt, alt_want, &at, &borrowed) != DUCKVEP_CODING_CONTEXT_OK ||
        memcmp(&borrowed, &expected, sizeof borrowed) != 0) return THEFT_TRIAL_FAIL;
    if (ctx.alt_cds_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                         (duckvep_codon_table_t)s->ctab,
                                         alt_cds, ctx.alt_cds_len - 1u,
                                         ref_pep, sizeof ref_pep,
                                         alt_pep, sizeof alt_pep,
                                         &fail_ctx) !=
            DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_coding_context_cov.capfail++;
    }
    if (ctx.ref_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                         (duckvep_codon_table_t)s->ctab,
                                         alt_cds, sizeof alt_cds,
                                         ref_pep, ctx.ref_peptide_len,
                                         alt_pep, sizeof alt_pep,
                                         &fail_ctx) !=
            DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_coding_context_cov.capfail++;
    }
    if (ctx.alt_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                         (duckvep_codon_table_t)s->ctab,
                                         alt_cds, sizeof alt_cds,
                                         ref_pep, sizeof ref_pep,
                                         alt_pep, ctx.alt_peptide_len,
                                         &fail_ctx) !=
            DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_coding_context_cov.capfail++;
    }

    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_coding_context_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_coding_context_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_coding_context_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_coding_context_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_coding_context_cov.indel++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_coding_context_cov.fwd++; else g_coding_context_cov.rev++;
    if (rf == 0u && af == 0u) g_coding_context_cov.pep_same++;
    else g_coding_context_cov.pep_diff++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_matches_direct_oracles(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context == direct CDS splice + full peptide oracles";
    cfg.prop1 = prop_coding_context_matches_direct_oracles;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_coding_context_cov, 0, sizeof g_coding_context_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_coding_context_cov.snv > 0u);
    ASSERT(g_coding_context_cov.mnv > 0u);
    ASSERT(g_coding_context_cov.ins > 0u);
    ASSERT(g_coding_context_cov.del > 0u);
    ASSERT(g_coding_context_cov.indel > 0u);
    ASSERT(g_coding_context_cov.fwd > 0u);
    ASSERT(g_coding_context_cov.rev > 0u);
    ASSERT(g_coding_context_cov.pep_diff > 0u);
    ASSERT(g_coding_context_cov.capfail > 0u);
    fprintf(stderr,
            "[coding-context coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u pep_same=%u pep_diff=%u capfail=%u\n",
            g_coding_context_cov.snv, g_coding_context_cov.mnv,
            g_coding_context_cov.ins, g_coding_context_cov.del,
            g_coding_context_cov.indel, g_coding_context_cov.fwd,
            g_coding_context_cov.rev, g_coding_context_cov.pep_same,
            g_coding_context_cov.pep_diff, g_coding_context_cov.capfail);
    PASS();
}

int kprop_hgvs_protein_fact_replay(
    const duckvep_hgvs_protein_fact_t *fact,
    const uint8_t                     *reference,
    size_t                             reference_length,
    uint8_t                           *alternate,
    size_t                             alternate_capacity,
    size_t                            *alternate_length) {

    size_t edit_start0;
    size_t removed_length;
    size_t inserted_length;
    size_t suffix_start0;
    size_t output_length;
    size_t i;

    if (alternate_length != NULL) *alternate_length = 0u;
    if (fact == NULL || reference == NULL || alternate == NULL ||
        alternate_length == NULL || fact->first_position1 == 0u ||
        fact->last_position1 == 0u) {
        return 0;
    }
    if (fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL) {
        if (reference_length > alternate_capacity) return 0;
        memcpy(alternate, reference, reference_length);
        *alternate_length = reference_length;
        return 1;
    }

    inserted_length = 0u;
    removed_length = 0u;
    if (fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION ||
        fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION ||
        fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELINS) {
        if (fact->last_position1 < fact->first_position1) return 0;
        edit_start0 = (size_t)fact->first_position1 - 1u;
        removed_length = (size_t)fact->last_position1 -
                         (size_t)fact->first_position1 + 1u;
        inserted_length = fact->shape ==
                (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION
            ? 0u : fact->alt_length;
        if (removed_length != fact->ref_length) return 0;
        for (i = 0u; i < removed_length; i++) {
            uint8_t residue;
            if (edit_start0 + i >= reference_length ||
                duckvep_hgvs_protein_base(fact, 0, i, &residue) !=
                    DUCKVEP_HGVS_OK ||
                residue != reference[edit_start0 + i]) {
                return 0;
            }
        }
    } else if (fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION) {
        if (fact->last_position1 != fact->first_position1 + 1u) return 0;
        edit_start0 = (size_t)fact->first_position1;
        inserted_length = fact->alt_length;
    } else if (fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DUPLICATION) {
        size_t source_start0;
        if (fact->last_position1 < fact->first_position1) return 0;
        edit_start0 = (size_t)fact->last_position1;
        source_start0 = (size_t)fact->first_position1 - 1u;
        inserted_length = fact->alt_length;
        if ((size_t)fact->last_position1 -
                (size_t)fact->first_position1 + 1u != inserted_length) {
            return 0;
        }
        for (i = 0u; i < inserted_length; i++) {
            uint8_t residue;
            if (source_start0 + i >= reference_length ||
                duckvep_hgvs_protein_base(fact, 1, i, &residue) !=
                    DUCKVEP_HGVS_OK ||
                residue != reference[source_start0 + i]) {
                return 0;
            }
        }
    } else {
        return 0;
    }

    if (edit_start0 > reference_length ||
        removed_length > reference_length - edit_start0 ||
        inserted_length > SIZE_MAX -
            (reference_length - removed_length)) {
        return 0;
    }
    output_length = reference_length - removed_length + inserted_length;
    if (output_length > alternate_capacity) return 0;
    suffix_start0 = edit_start0 + removed_length;
    memcpy(alternate, reference, edit_start0);
    for (i = 0u; i < inserted_length; i++) {
        if (duckvep_hgvs_protein_base(fact, 1, i,
                                      alternate + edit_start0 + i) !=
            DUCKVEP_HGVS_OK) {
            return 0;
        }
    }
    memcpy(alternate + edit_start0 + inserted_length,
           reference + suffix_start0,
           reference_length - suffix_start0);
    *alternate_length = output_length;
    return 1;
}

static int kprop_hgvs_is_vep_absent_terminal_insertion(
    const duckvep_coding_context_t *context) {

    duckvep_coding_peptide_window_t window;
    size_t prefix = 0u;
    size_t suffix = 0u;
    size_t ref_remaining;
    size_t alt_remaining;
    size_t first_position1;
    size_t last_position1;
    size_t low;
    size_t reference_length;

    if (context == NULL ||
        !duckvep_coding_context_peptide_window_open(context, &window)) {
        return 0;
    }
    while (prefix < window.ref_length && prefix < window.alt_length) {
        uint8_t reference = duckvep_coding_context_peptide_window_base(
            context, &window, 0, prefix);
        uint8_t alternate = duckvep_coding_context_peptide_window_base(
            context, &window, 1, prefix);
        if (reference == 0u || alternate == 0u ||
            (reference == (uint8_t)'*' && alternate == (uint8_t)'*')) {
            return 0;
        }
        if (reference != alternate) break;
        prefix++;
    }
    while (suffix < window.ref_length - prefix &&
           suffix < window.alt_length - prefix) {
        uint8_t reference = duckvep_coding_context_peptide_window_base(
            context, &window, 0, window.ref_length - 1u - suffix);
        uint8_t alternate = duckvep_coding_context_peptide_window_base(
            context, &window, 1, window.alt_length - 1u - suffix);
        if (reference == 0u || alternate == 0u || reference != alternate) {
            break;
        }
        suffix++;
    }
    ref_remaining = window.ref_length - prefix - suffix;
    alt_remaining = window.alt_length - prefix - suffix;
    if (ref_remaining != 0u || alt_remaining == 0u) return 0;
    first_position1 = window.ref_peptide_offset + 1u + prefix;
    last_position1 = window.ref_peptide_offset + window.ref_length - suffix;
    low = first_position1 < last_position1
        ? first_position1 : last_position1;
    reference_length = context->ref_peptide_len;
    if (reference_length != 0u &&
        duckvep_coding_context_peptide_base(
            context, 0, reference_length - 1u) == (uint8_t)'*') {
        reference_length--;
    }
    if (low == 0u || low > reference_length) return 1;
    if (low == reference_length) {
        uint8_t original_reference_first =
            window.ref_length == 0u ? 0u :
            duckvep_coding_context_peptide_window_base(
                context, &window, 0, 0u);
        return duckvep_coding_context_peptide_base(
                   context, 0, reference_length) == 0u ||
               original_reference_first != (uint8_t)'*';
    }
    return 0;
}

static struct {
    uint32_t replayed;
    uint32_t equal;
    uint32_t substitution;
    uint32_t deletion;
    uint32_t insertion;
    uint32_t delins;
    uint32_t duplication;
    uint32_t special;
    uint32_t terminal_not_applicable;
    uint32_t vep_stop_equal;
    uint32_t vep_position_zero;
    uint32_t fwd;
    uint32_t rev;
} g_hgvs_protein_replay_cov;

static enum theft_trial_res prop_hgvs_protein_fact_replays_translated_edit(
    struct theft *t,
    void         *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t scratch[8];
    duckvep_edit_set_t edit_set;
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t fact;
    uint8_t alt_cds[80];
    uint8_t ref_peptide[32];
    uint8_t alt_peptide[32];
    uint8_t ref_want[32];
    uint8_t alt_want[32];
    uint8_t replayed[32];
    size_t ref_want_length = 0u;
    size_t alt_want_length = 0u;
    size_t replayed_length = 0u;
    duckvep_context_delta_status_t delta_status;
    duckvep_hgvs_status_t hgvs_status;
    (void)t;

    if (duckvep_variant_cds_edit_set_build(
            &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand,
            scratch, sizeof scratch / sizeof scratch[0], &edit_set) !=
            DUCKVEP_CDS_EDIT_OK ||
        duckvep_coding_context_build(
            s->cds, s->cds_lenv, &edit_set, s->strand,
            (duckvep_codon_table_t)s->ctab,
            alt_cds, sizeof alt_cds,
            ref_peptide, sizeof ref_peptide,
            alt_peptide, sizeof alt_peptide, &context) !=
            DUCKVEP_CODING_CONTEXT_OK ||
        !kprop_translate_full_oracle(
            s->cds, s->cds_lenv, (duckvep_codon_table_t)s->ctab,
            ref_want, &ref_want_length) ||
        !kprop_translate_full_oracle(
            s->expect_cds, (size_t)s->expect_len,
            (duckvep_codon_table_t)s->ctab,
            alt_want, &alt_want_length)) {
        return THEFT_TRIAL_FAIL;
    }
    delta_status = duckvep_coding_context_delta_fill(
        &context, 0u, &delta);
    if (delta_status != DUCKVEP_CONTEXT_DELTA_OK) {
        g_hgvs_protein_replay_cov.special++;
        return THEFT_TRIAL_PASS;
    }
    if (delta.valid == 0u) return THEFT_TRIAL_FAIL;
    if ((context.length_diff % 3) != 0 || delta.frameshift ||
        delta.start_lost || delta.stop_lost) {
        g_hgvs_protein_replay_cov.special++;
        return THEFT_TRIAL_PASS;
    }
    hgvs_status = duckvep_hgvs_protein_fact_build(
        &context, &delta, &fact);
    if (hgvs_status != DUCKVEP_HGVS_OK) {
        if (hgvs_status == DUCKVEP_HGVS_NOT_APPLICABLE &&
            kprop_hgvs_is_vep_absent_terminal_insertion(&context)) {
            g_hgvs_protein_replay_cov.terminal_not_applicable++;
            return THEFT_TRIAL_PASS;
        }
        fprintf(stderr,
            "[HGVSp unexpected status] status=%u kind=%u strand=%d "
            "diff=%" PRId64 " edit=%u/%u/%u ref=%.*s alt=%.*s\n",
            (unsigned)hgvs_status, (unsigned)s->vkind, (int)s->strand,
            context.length_diff,
            context.single_edit_cds_start,
            context.single_edit_ref_len,
            context.single_edit_alt_len,
            (int)ref_want_length, (const char *)ref_want,
            (int)alt_want_length, (const char *)alt_want);
        return THEFT_TRIAL_FAIL;
    }
    if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT ||
        fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST ||
        fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_EXTENSION) {
        g_hgvs_protein_replay_cov.special++;
        return THEFT_TRIAL_PASS;
    }
    if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL &&
        (ref_want_length != alt_want_length ||
         memcmp(ref_want, alt_want, ref_want_length) != 0)) {
        /* VEP's _clip_alleles returns equality as soon as both local peptide
         * strings meet a stop, even when later translated codons differ. */
        g_hgvs_protein_replay_cov.vep_stop_equal++;
        return THEFT_TRIAL_PASS;
    }
    if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION &&
        fact.first_position1 == 0u && fact.last_position1 == 1u) {
        uint8_t terminal_reference;
        size_t reference_length = context.ref_peptide_len;
        if (reference_length != 0u &&
            ref_want[reference_length - 1u] == (uint8_t)'*') {
            reference_length--;
        }
        if (reference_length == 0u) return THEFT_TRIAL_FAIL;
        terminal_reference = ref_want[reference_length - 1u];
        if (fact.reference_first != terminal_reference ||
            fact.reference_last != terminal_reference) {
            return THEFT_TRIAL_FAIL;
        }
        /* Perl substr(_peptide, -1, 2) supplies presentation residues from
         * the final reference amino acid. This VEP compatibility string is
         * not a semantic edit over HGVS position zero, so the independent
         * replay oracle records it separately. */
        g_hgvs_protein_replay_cov.vep_position_zero++;
        g_hgvs_protein_replay_cov.special++;
        return THEFT_TRIAL_PASS;
    }
    if (!kprop_hgvs_protein_fact_replay(
            &fact, ref_want, ref_want_length,
            replayed, sizeof replayed, &replayed_length) ||
        replayed_length != alt_want_length ||
        memcmp(replayed, alt_want, alt_want_length) != 0) {
        fprintf(stderr,
            "[HGVSp mismatch] kind=%u strand=%d shape=%u pos=%u..%u "
            "lens=%zu/%zu ref=%.*s alt=%.*s replay=%.*s\n",
            (unsigned)s->vkind, (int)s->strand, (unsigned)fact.shape,
            fact.first_position1, fact.last_position1,
            fact.ref_length, fact.alt_length,
            (int)ref_want_length, (const char *)ref_want,
            (int)alt_want_length, (const char *)alt_want,
            (int)replayed_length, (const char *)replayed);
        return THEFT_TRIAL_FAIL;
    }

    if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL) {
        g_hgvs_protein_replay_cov.equal++;
    } else if (fact.shape ==
               (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION) {
        g_hgvs_protein_replay_cov.substitution++;
    } else if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION) {
        g_hgvs_protein_replay_cov.deletion++;
    } else if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION) {
        g_hgvs_protein_replay_cov.insertion++;
    } else if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELINS) {
        g_hgvs_protein_replay_cov.delins++;
    } else if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DUPLICATION) {
        g_hgvs_protein_replay_cov.duplication++;
    } else {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_hgvs_protein_replay_cov.fwd++;
    else g_hgvs_protein_replay_cov.rev++;
    g_hgvs_protein_replay_cov.replayed++;
    return THEFT_TRIAL_PASS;
}

TEST hgvs_protein_facts_replay_independent_translation(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    memset(&g_hgvs_protein_replay_cov, 0,
           sizeof g_hgvs_protein_replay_cov);
    cfg.name = "HGVSp fact replay == independently translated edited CDS";
    cfg.prop1 = prop_hgvs_protein_fact_replays_translated_edit;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_hgvs_protein_replay_cov.replayed > 0u);
    ASSERT(g_hgvs_protein_replay_cov.equal > 0u);
    ASSERT(g_hgvs_protein_replay_cov.substitution > 0u);
    ASSERT(g_hgvs_protein_replay_cov.deletion > 0u);
    ASSERT(g_hgvs_protein_replay_cov.insertion > 0u);
    ASSERT(g_hgvs_protein_replay_cov.delins > 0u);
    ASSERT(g_hgvs_protein_replay_cov.duplication > 0u);
    ASSERT(g_hgvs_protein_replay_cov.terminal_not_applicable > 0u);
    ASSERT(g_hgvs_protein_replay_cov.vep_position_zero > 0u);
    ASSERT(g_hgvs_protein_replay_cov.fwd > 0u);
    ASSERT(g_hgvs_protein_replay_cov.rev > 0u);
    fprintf(stderr,
        "[HGVSp replay coverage] replayed=%u equal=%u sub=%u del=%u "
        "ins=%u delins=%u dup=%u special=%u terminal_not_applicable=%u "
        "vep_stop_equal=%u vep_position_zero=%u "
        "fwd=%u rev=%u\n",
        g_hgvs_protein_replay_cov.replayed,
        g_hgvs_protein_replay_cov.equal,
        g_hgvs_protein_replay_cov.substitution,
        g_hgvs_protein_replay_cov.deletion,
        g_hgvs_protein_replay_cov.insertion,
        g_hgvs_protein_replay_cov.delins,
        g_hgvs_protein_replay_cov.duplication,
        g_hgvs_protein_replay_cov.special,
        g_hgvs_protein_replay_cov.terminal_not_applicable,
        g_hgvs_protein_replay_cov.vep_stop_equal,
        g_hgvs_protein_replay_cov.vep_position_zero,
        g_hgvs_protein_replay_cov.fwd,
        g_hgvs_protein_replay_cov.rev);
    PASS();
}

static struct {
    uint32_t eligible;
    uint32_t insertion;
    uint32_t deletion;
    uint32_t delins;
    uint32_t fwd;
    uint32_t rev;
    uint32_t frameshift;
    uint32_t immediate_stop;
    uint32_t equal_stop;
    uint32_t shortened;
    uint32_t termination_known;
    uint32_t termination_unknown;
    uint32_t non_frameshift;
} g_hgvs_frameshift_cov;

static enum theft_trial_res prop_hgvs_frameshift_matches_extended_translation(
    struct theft *t,
    void         *arg1) {

    static const uint8_t SAFE_CODONS[4][3] = {
        {'G','C','C'}, {'C','A','A'}, {'G','A','T'}, {'T','T','T'}
    };
    static const uint8_t STOP_CODON[3] = {'T','A','A'};
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[8];
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t fact;
    uint8_t alt_cds[80];
    uint8_t ref_peptide[32];
    uint8_t alt_peptide[32];
    uint8_t tail[24];
    uint8_t extended_cds[128];
    uint8_t ref_want[32];
    uint8_t alt_want[64];
    size_t tail_codons = (size_t)kprop_bounded(t, 6u) + 1u;
    size_t tail_length = (tail_codons + 1u) * 3u;
    size_t extended_alt_length;
    size_t extended_cds_length;
    size_t ref_want_length = 0u;
    size_t alt_want_length = 0u;
    size_t reference_length;
    size_t position0;
    size_t i;
    uint8_t expected_reference;
    uint8_t expected_alternate;
    uint8_t expected_shape;
    uint8_t termination_known = 0u;
    uint32_t termination_distance = 0u;

    for (i = 0u; i < tail_codons; i++) {
        memcpy(tail + i * 3u,
               SAFE_CODONS[kprop_bounded(t, 4u)], 3u);
    }
    memcpy(tail + tail_codons * 3u, STOP_CODON, 3u);
    if (duckvep_variant_physical_coding_context_build(
            &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand,
            edits, sizeof edits / sizeof edits[0],
            alt_cds, sizeof alt_cds,
            ref_peptide, sizeof ref_peptide,
            alt_peptide, sizeof alt_peptide, &context) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    context.post_cds_bases = tail;
    context.post_cds_length = tail_length;
    context.post_cds_complete = 1u;
    if (duckvep_coding_context_delta_fill(&context, 0u, &delta) !=
            DUCKVEP_CONTEXT_DELTA_OK || delta.valid == 0u) {
        return THEFT_TRIAL_FAIL;
    }
    if (!delta.frameshift) {
        g_hgvs_frameshift_cov.non_frameshift++;
        return THEFT_TRIAL_PASS;
    }
    if (duckvep_hgvs_protein_fact_build(&context, &delta, &fact) !=
            DUCKVEP_HGVS_OK) {
        return THEFT_TRIAL_FAIL;
    }

    /* VEP 116's _trim_incomplete_codon assigns, rather than compares, its
     * keep length. Any alternate CDS of at least one codon therefore reaches
     * _get_fs_peptides untrimmed, including its terminal remainder. A one- or
     * two-base alternate CDS is assigned keep_length zero and is dropped
     * before the downstream transcript sequence is appended. */
    extended_alt_length = context.alt_cds_len < 3u
        ? 0u : context.alt_cds_len;
    if (extended_alt_length > sizeof extended_cds ||
        tail_length > sizeof extended_cds - extended_alt_length) {
        return THEFT_TRIAL_ERROR;
    }
    memcpy(extended_cds, context.alt_cds, extended_alt_length);
    memcpy(extended_cds + extended_alt_length, tail, tail_length);
    extended_cds_length = extended_alt_length + tail_length;
    if (!kprop_translate_full_oracle(
            context.ref_cds, context.ref_cds_len,
            (duckvep_codon_table_t)context.codon_table,
            ref_want, &ref_want_length) ||
        !kprop_translate_full_oracle(
            extended_cds, extended_cds_length,
            (duckvep_codon_table_t)context.codon_table,
            alt_want, &alt_want_length)) {
        return THEFT_TRIAL_FAIL;
    }
    reference_length = ref_want_length;
    if (reference_length != 0u &&
        ref_want[reference_length - 1u] == (uint8_t)'*') {
        reference_length--;
    }
    position0 = (size_t)(
        (context.single_edit_cds_start - 1u) / 3u);
    if (position0 >= alt_want_length) {
        if (position0 > reference_length ||
            fact.shape != (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION) {
            return THEFT_TRIAL_FAIL;
        }
        g_hgvs_frameshift_cov.shortened++;
        expected_shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION;
        expected_reference = position0 < reference_length
            ? ref_want[position0] : (uint8_t)'*';
        expected_alternate = 0u;
    } else {
        for (;;) {
            expected_reference = position0 < reference_length
                ? ref_want[position0]
                : position0 == reference_length ? (uint8_t)'*' : 0u;
            expected_alternate = alt_want[position0];
            if (expected_reference == (uint8_t)'*' &&
                expected_alternate == (uint8_t)'*') {
                expected_shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL;
                g_hgvs_frameshift_cov.equal_stop++;
                break;
            }
            if (expected_reference == 0u ||
                expected_reference != expected_alternate) {
                expected_shape = expected_alternate == (uint8_t)'*'
                    ? (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION
                    : (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
                if (expected_shape ==
                        (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION) {
                    g_hgvs_frameshift_cov.immediate_stop++;
                } else {
                    g_hgvs_frameshift_cov.frameshift++;
                }
                break;
            }
            position0++;
            if (position0 >= alt_want_length) return THEFT_TRIAL_FAIL;
        }
    }
    if (fact.shape != expected_shape ||
        fact.first_position1 != (uint32_t)position0 + 1u ||
        fact.last_position1 != fact.first_position1 ||
        fact.reference_first != expected_reference ||
        (expected_shape != (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION &&
         fact.alternate_first != expected_alternate)) {
        return THEFT_TRIAL_FAIL;
    }
    if (expected_shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT) {
        for (i = 0u; i < alt_want_length; i++) {
            int64_t distance;
            if (alt_want[i] != (uint8_t)'*') continue;
            distance = (int64_t)i + 1 - (int64_t)position0;
            if (distance > 0 && (uint64_t)distance <= UINT32_MAX) {
                termination_known = 1u;
                termination_distance = (uint32_t)distance;
            }
            break;
        }
        if (fact.termination_known != termination_known ||
            (termination_known &&
             fact.termination_distance != termination_distance)) {
            return THEFT_TRIAL_FAIL;
        }
        if (termination_known) g_hgvs_frameshift_cov.termination_known++;
        else g_hgvs_frameshift_cov.termination_unknown++;
    }
    if (s->vkind == (uint8_t)DUCKVEP_KIND_INS) {
        g_hgvs_frameshift_cov.insertion++;
    } else if (s->vkind == (uint8_t)DUCKVEP_KIND_DEL) {
        g_hgvs_frameshift_cov.deletion++;
    } else if (s->vkind == (uint8_t)DUCKVEP_KIND_INDEL) {
        g_hgvs_frameshift_cov.delins++;
    } else {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_hgvs_frameshift_cov.fwd++;
    else g_hgvs_frameshift_cov.rev++;
    g_hgvs_frameshift_cov.eligible++;
    return THEFT_TRIAL_PASS;
}

TEST hgvs_frameshift_facts_match_extended_translation(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    memset(&g_hgvs_frameshift_cov, 0, sizeof g_hgvs_frameshift_cov);
    cfg.name = "HGVSp frameshift fact == independently extended translation";
    cfg.prop1 = prop_hgvs_frameshift_matches_extended_translation;
    cfg.type_info[0] = &kprop_frameshift_indel_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_hgvs_frameshift_cov.eligible > 0u);
    ASSERT(g_hgvs_frameshift_cov.insertion > 0u);
    ASSERT(g_hgvs_frameshift_cov.deletion > 0u);
    ASSERT(g_hgvs_frameshift_cov.delins > 0u);
    ASSERT(g_hgvs_frameshift_cov.fwd > 0u);
    ASSERT(g_hgvs_frameshift_cov.rev > 0u);
    ASSERT(g_hgvs_frameshift_cov.frameshift > 0u);
    ASSERT(g_hgvs_frameshift_cov.immediate_stop > 0u);
    ASSERT(g_hgvs_frameshift_cov.termination_known > 0u);
    fprintf(stderr,
        "[HGVSp frameshift coverage] eligible=%u ins=%u del=%u delins=%u "
        "fwd=%u rev=%u fs=%u immediate_stop=%u equal_stop=%u shortened=%u "
        "ter_known=%u ter_unknown=%u non_fs=%u\n",
        g_hgvs_frameshift_cov.eligible,
        g_hgvs_frameshift_cov.insertion,
        g_hgvs_frameshift_cov.deletion,
        g_hgvs_frameshift_cov.delins,
        g_hgvs_frameshift_cov.fwd,
        g_hgvs_frameshift_cov.rev,
        g_hgvs_frameshift_cov.frameshift,
        g_hgvs_frameshift_cov.immediate_stop,
        g_hgvs_frameshift_cov.equal_stop,
        g_hgvs_frameshift_cov.shortened,
        g_hgvs_frameshift_cov.termination_known,
        g_hgvs_frameshift_cov.termination_unknown,
        g_hgvs_frameshift_cov.non_frameshift);
    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t pep_same;
    uint32_t pep_diff;
    uint32_t capfail;
} g_variant_coding_context_cov;

static enum theft_trial_res prop_variant_coding_context_matches_oracles(struct theft *t,
                                                                        void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit_scratch[8];
    duckvep_haplotype_edit_t expected_edits[8];
    duckvep_edit_set_t expected_set;
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint8_t ref_want[32];
    uint8_t alt_want[32];
    size_t ref_want_len = 0u;
    size_t alt_want_len = 0u;
    int cds_changed;
    uint32_t rf, rl, af, al;
    (void)t;

    if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                           0u, 0u, s->strand, expected_edits, 8u,
                                           &expected_set) != DUCKVEP_CDS_EDIT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand, edit_scratch, 8u,
                                             alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_cds != s->cds || ctx.alt_cds != alt_cds) return THEFT_TRIAL_FAIL;
    if (ctx.alt_cds_len != (size_t)s->expect_len ||
        memcmp(ctx.alt_cds, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->cds, s->cds_lenv,
                                     (duckvep_codon_table_t)s->ctab,
                                     ref_want, &ref_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->expect_cds, (size_t)s->expect_len,
                                     (duckvep_codon_table_t)s->ctab,
                                     alt_want, &alt_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide_len != ref_want_len || ctx.alt_peptide_len != alt_want_len) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide[ctx.ref_peptide_len] != (uint8_t)'\0' ||
        ctx.alt_peptide[ctx.alt_peptide_len] != (uint8_t)'\0') {
        return THEFT_TRIAL_FAIL;
    }
    if (memcmp(ctx.ref_peptide, ref_want, ref_want_len) != 0 ||
        memcmp(ctx.alt_peptide, alt_want, alt_want_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    cds_changed = kprop_cds_changed_oracle(s->cds, s->cds_lenv,
                                           s->expect_cds, (size_t)s->expect_len);
    if (cds_changed < 0 || ctx.cds_changed != (uint8_t)cds_changed) return THEFT_TRIAL_FAIL;
    if (ctx.applied_edits != expected_set.count) return THEFT_TRIAL_FAIL;
    if (ctx.length_diff != (int64_t)s->expect_len - (int64_t)s->cds_lenv) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.flags != kprop_context_flags_oracle(&expected_set)) return THEFT_TRIAL_FAIL;
    if (ctx.flags != kprop_single_variant_flags_oracle(s->cds_lenv, s->expect_len)) {
        return THEFT_TRIAL_FAIL;
    }
    kprop_peptide_window_oracle(ref_want, ref_want_len, alt_want, alt_want_len,
                                &rf, &rl, &af, &al);
    if (ctx.ref_first_changed_codon != rf || ctx.ref_last_changed_codon != rl ||
        ctx.alt_first_changed_codon != af || ctx.alt_last_changed_codon != al) {
        return THEFT_TRIAL_FAIL;
    }

    {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, NULL, 0u,
                                                 alt_cds, sizeof alt_cds,
                                                 ref_pep, sizeof ref_pep,
                                                 alt_pep, sizeof alt_pep,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }
    if (ctx.alt_cds_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, edit_scratch, 8u,
                                                 alt_cds, ctx.alt_cds_len - 1u,
                                                 ref_pep, sizeof ref_pep,
                                                 alt_pep, sizeof alt_pep,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }
    if (ctx.ref_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, edit_scratch, 8u,
                                                 alt_cds, sizeof alt_cds,
                                                 ref_pep, ctx.ref_peptide_len,
                                                 alt_pep, sizeof alt_pep,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }
    if (ctx.alt_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, edit_scratch, 8u,
                                                 alt_cds, sizeof alt_cds,
                                                 ref_pep, sizeof ref_pep,
                                                 alt_pep, ctx.alt_peptide_len,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }

    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_variant_coding_context_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_variant_coding_context_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_variant_coding_context_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_variant_coding_context_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_variant_coding_context_cov.indel++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_variant_coding_context_cov.fwd++;
    else g_variant_coding_context_cov.rev++;
    if (rf == 0u && af == 0u) g_variant_coding_context_cov.pep_same++;
    else g_variant_coding_context_cov.pep_diff++;
    return THEFT_TRIAL_PASS;
}

TEST variant_coding_context_matches_oracles(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant coding context == direct CDS splice + full peptide oracles";
    cfg.prop1 = prop_variant_coding_context_matches_oracles;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_variant_coding_context_cov, 0, sizeof g_variant_coding_context_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_variant_coding_context_cov.snv > 0u);
    ASSERT(g_variant_coding_context_cov.mnv > 0u);
    ASSERT(g_variant_coding_context_cov.ins > 0u);
    ASSERT(g_variant_coding_context_cov.del > 0u);
    ASSERT(g_variant_coding_context_cov.indel > 0u);
    ASSERT(g_variant_coding_context_cov.fwd > 0u);
    ASSERT(g_variant_coding_context_cov.rev > 0u);
    ASSERT(g_variant_coding_context_cov.pep_diff > 0u);
    ASSERT(g_variant_coding_context_cov.capfail > 0u);
    fprintf(stderr,
            "[variant-coding-context coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u pep_same=%u pep_diff=%u capfail=%u\n",
            g_variant_coding_context_cov.snv, g_variant_coding_context_cov.mnv,
            g_variant_coding_context_cov.ins, g_variant_coding_context_cov.del,
            g_variant_coding_context_cov.indel, g_variant_coding_context_cov.fwd,
            g_variant_coding_context_cov.rev, g_variant_coding_context_cov.pep_same,
            g_variant_coding_context_cov.pep_diff,
            g_variant_coding_context_cov.capfail);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
} g_context_delta_cov;

static enum theft_trial_res prop_context_delta_matches_codon_oracle(struct theft *t,
                                                                    void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit_scratch[4];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    size_t first_diff = 0u;
    size_t codon_start;
    size_t codon_idx;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    uint32_t j;
    (void)t;

    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edit_scratch, 4u,
                                             alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_coding_context_delta_fill(&ctx, s->flags, &delta) !=
        DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    while (first_diff < s->cds_lenv && s->cds[first_diff] == s->expect_cds[first_diff]) {
        first_diff++;
    }
    if (first_diff >= s->cds_lenv) return THEFT_TRIAL_FAIL;
    codon_start = first_diff - (first_diff % 3u);
    codon_idx = codon_start / 3u;
    for (j = 0u; j < 3u; j++) {
        ref_codon[j] = (char)s->cds[codon_start + (size_t)j];
        alt_codon[j] = (char)s->expect_cds[codon_start + (size_t)j];
    }
    ref_codon[3] = '\0';
    alt_codon[3] = '\0';
    cr = duckvep_codon_change(ref_codon, alt_codon, DUCKVEP_CODON_TABLE_STANDARD);
    if (cr.change & DUCKVEP_CODON_INVALID) return THEFT_TRIAL_FAIL;
    if (!delta.valid || delta.cdna_pos != -1 || delta.cds_pos != -1 ||
        delta.protein_pos != (int32_t)codon_idx + 1 ||
        delta.ref_aa != (uint8_t)cr.aa_ref || delta.alt_aa != (uint8_t)cr.aa_alt) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) {
        if (!delta.synonymous || delta.missense || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.syn++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) {
        if (!delta.missense || delta.synonymous || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.mis++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) {
        if (!delta.stop_gained || delta.synonymous || delta.missense ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.stop_gained++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) {
        if (!delta.stop_lost || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.stop_lost++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) {
        if (!delta.stop_retained || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_lost) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.stop_retained++;
    } else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_context_delta_cov.fwd++;
    else g_context_delta_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_matches_codon_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delta == single-codon oracle";
    cfg.prop1 = prop_context_delta_matches_codon_oracle;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_delta_cov, 0, sizeof g_context_delta_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_delta_cov.syn > 0u);
    ASSERT(g_context_delta_cov.mis > 0u);
    ASSERT(g_context_delta_cov.stop_gained > 0u);
    ASSERT(g_context_delta_cov.stop_lost > 0u);
    ASSERT(g_context_delta_cov.stop_retained > 0u);
    ASSERT(g_context_delta_cov.fwd > 0u);
    ASSERT(g_context_delta_cov.rev > 0u);
    fprintf(stderr,
            "[context-delta coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u\n",
            g_context_delta_cov.syn, g_context_delta_cov.mis,
            g_context_delta_cov.stop_gained, g_context_delta_cov.stop_lost,
            g_context_delta_cov.stop_retained, g_context_delta_cov.fwd,
            g_context_delta_cov.rev);
    PASS();
}

static int kprop_delta_is_frameshift_at(const duckvep_sequence_delta_t *d,
                                        int32_t protein_pos) {
    return d != NULL && d->valid && d->frameshift &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->inframe_deletion && !d->inframe_insertion && !d->protein_altering &&
           !d->coding_unknown && !d->partial_codon &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

/* Independent oracle for VEP's LOCAL-window frameshift stop_gained. Re-derives the
 * codon window from the raw CDS bytes (substr(alt_cds, codon_cds_start-1,
 * codon_len + net_delta)) and tests whether a WHOLE codon there is a standard stop the
 * reference window lacked. Deliberately re-translates from bytes (stops enumerated by
 * hand for the standard table) rather than reading the kernel's pre-translated peptides,
 * so it does not restate the implementation. Only valid for the standard codon table. */
static int kprop_is_standard_stop(const uint8_t *b) {
    return b[0] == (uint8_t)'T' &&
           ((b[1] == (uint8_t)'A' && (b[2] == (uint8_t)'A' || b[2] == (uint8_t)'G')) ||
            (b[1] == (uint8_t)'G' && b[2] == (uint8_t)'A'));
}

static int kprop_frameshift_local_stop_oracle(
    const uint8_t *ref_cds, size_t ref_cds_len,
    const uint8_t *alt_cds, size_t alt_cds_len,
    uint32_t cds_start, uint32_t ref_len, int64_t length_diff) {

    uint64_t first = cds_start;
    uint64_t last, tv_s, tv_e, codon_cds_start, codon_len, off, i, whole;
    int64_t win;
    int alt_stop = 0;
    int ref_stop = 0;

    if (first == 0u) return 0;
    tv_s = ((first - 1u) / 3u) + 1u;
    codon_cds_start = tv_s * 3u - 2u;
    off = codon_cds_start - 1u;
    /* Pure insertion before the first CDS base spans no reference codon; never form
     * last == 0 (whose (last - 1) would underflow), matching the kernel helper. */
    if (ref_len > 0u) {
        last = first + ref_len - 1u;
        tv_e = ((last - 1u) / 3u) + 1u;
        codon_len = (tv_e >= tv_s) ? (tv_e - tv_s + 1u) * 3u : 0u;
    } else if (first > 1u) {
        last = first - 1u;
        tv_e = ((last - 1u) / 3u) + 1u;
        codon_len = (tv_e >= tv_s) ? (tv_e - tv_s + 1u) * 3u : 0u;
    } else {
        codon_len = 0u;
    }

    win = (int64_t)codon_len + length_diff;
    if (win > 0 && off <= alt_cds_len) {
        uint64_t avail = (uint64_t)alt_cds_len - off;
        if ((uint64_t)win > avail) win = (int64_t)avail;
        whole = (uint64_t)win / 3u;
        for (i = 0u; i < whole; i++) {
            if (kprop_is_standard_stop(alt_cds + off + i * 3u)) { alt_stop = 1; break; }
        }
    }
    if (codon_len > 0u && off <= ref_cds_len) {
        uint64_t avail = (uint64_t)ref_cds_len - off;
        uint64_t rw = codon_len;
        if (rw > avail) rw = avail;
        whole = rw / 3u;
        for (i = 0u; i < whole; i++) {
            if (kprop_is_standard_stop(ref_cds + off + i * 3u)) { ref_stop = 1; break; }
        }
    }
    return alt_stop && !ref_stop;
}

/* Deterministic anchor for VEP's local-window frameshift stop_gained: the recomputed
 * codon AT THE JUNCTION being a stop composites frameshift_variant&stop_gained, while a
 * frameshift whose junction codon is not a stop stays a bare frameshift EVEN WHEN the
 * shifted downstream frame contains a stop — the discriminator proving the engine uses
 * VEP's local window, not a downstream retranslation. */
TEST coding_context_delta_frameshift_stop_gained_scene(void) {
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[40];
    uint8_t ref_pep[20];
    uint8_t alt_pep[20];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;

    /* (A) +1 frameshift whose recomputed codon 3 becomes TAG (stop). Replace the 'C' at
     *     CDS pos 7 with "TA": ATG AAA CGT GAA -> ATG AAA TAG TGA A. */
    {
        static uint8_t cdsA[12] = {
            'A','T','G',  'A','A','A',  'C','G','T',  'G','A','A'
        };
        static const uint8_t insTA[2] = { 'T','A' };
        edit.cds_start = 7u; edit.ref = cdsA + 6u; edit.ref_len = 1u;
        edit.alt = insTA; edit.alt_len = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsA, sizeof cdsA, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
        ASSERT_EQ(3, delta.protein_pos);
        /* oracle agreement */
        ASSERT_EQ(1, kprop_frameshift_local_stop_oracle(cdsA, sizeof cdsA, ctx.alt_cds,
                                                         ctx.alt_cds_len, 7u, 1u,
                                                         ctx.length_diff));
    }

    /* (B) +1 frameshift whose recomputed codon 3 is CAA (not a stop) but the shifted
     *     frame carries TGA one codon downstream. stop_gained MUST be 0. Replace the 'C'
     *     at CDS pos 7 with "CA": ATG AAA CAT GAA -> ATG AAA CAA TGA A. */
    {
        static uint8_t cdsB[12] = {
            'A','T','G',  'A','A','A',  'C','A','T',  'G','A','A'
        };
        static const uint8_t insCA[2] = { 'C','A' };
        int downstream_stop = 0;
        size_t i;
        edit.cds_start = 7u; edit.ref = cdsB + 6u; edit.ref_len = 1u;
        edit.alt = insCA; edit.alt_len = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsB, sizeof cdsB, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.stop_gained);
        ASSERT_EQ(3, delta.protein_pos);
        /* non-vacuous: a downstream stop really is present in the shifted alt peptide,
         * proving the local window (not full retranslation) drove the 0 verdict. */
        for (i = 3u; i < ctx.alt_peptide_len; i++) {
            if (ctx.alt_peptide[i] == (uint8_t)'*') { downstream_stop = 1; break; }
        }
        ASSERT(downstream_stop);
        ASSERT_EQ(0, kprop_frameshift_local_stop_oracle(cdsB, sizeof cdsB, ctx.alt_cds,
                                                         ctx.alt_cds_len, 7u, 1u,
                                                         ctx.length_diff));
    }

    /* (C) PURE insertion (ref_len 0) — the dominant real ClinVar frameshift&stop_gained
     *     shape (e.g. C>CTTTAA). Insert "TAAG" before CDS pos 7 (codon boundary,
     *     codon_len 0): ATG AAA GGG CCC -> ATG AAA TAA GGG GCC C, recomputed codon 3 = TAA
     *     (stop). Exercises the ref_len==0 window branch and the first_cds>1 path. */
    {
        static uint8_t cdsC[12] = {
            'A','T','G',  'A','A','A',  'G','G','G',  'C','C','C'
        };
        static const uint8_t insTAAG[4] = { 'T','A','A','G' };
        edit.cds_start = 7u; edit.ref = NULL; edit.ref_len = 0u;
        edit.alt = insTAAG; edit.alt_len = 4u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsC, sizeof cdsC, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
        ASSERT_EQ(3, delta.protein_pos);
        ASSERT_EQ(1, kprop_frameshift_local_stop_oracle(cdsC, sizeof cdsC, ctx.alt_cds,
                                                        ctx.alt_cds_len, 7u, 0u,
                                                        ctx.length_diff));
    }

    /* (D) Codon-table awareness: the recomputed junction codon is TGA, a stop under the
     *     standard table but tryptophan under the vertebrate-mitochondrial table. The
     *     helper reads the peptide translated with the build's table, so stop_gained must
     *     flip. Replace the 'C' at CDS pos 7 with "TG": ATG AAA CAG GAA -> ATG AAA TGA GGA A. */
    {
        static uint8_t cdsD[12] = {
            'A','T','G',  'A','A','A',  'C','A','G',  'G','A','A'
        };
        static const uint8_t insTG[2] = { 'T','G' };
        edit.cds_start = 7u; edit.ref = cdsD + 6u; edit.ref_len = 1u;
        edit.alt = insTG; edit.alt_len = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsD, sizeof cdsD, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);

        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsD, sizeof cdsD, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_VERT_MITO,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.stop_gained);
    }

    /* (E) Regression guard for the pure-insertion first_cds==1 window (formerly a
     *     (last_cds-1) underflow). Insert "TAAC" before CDS pos 1 on a CDS_start_NF
     *     transcript (so the start-codon guard is bypassed): recomputed codon 1 = TAA
     *     (stop). Must resolve cleanly with no out-of-bounds read. */
    {
        static uint8_t cdsE[12] = {
            'A','T','G',  'A','A','A',  'G','G','G',  'C','C','C'
        };
        static const uint8_t insTAAC[4] = { 'T','A','A','C' };
        edit.cds_start = 1u; edit.ref = NULL; edit.ref_len = 0u;
        edit.alt = insTAAC; edit.alt_len = 4u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsE, sizeof cdsE, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx,
                                                    (uint64_t)DUCKVEP_TX_CDS_START_NF,
                                                    &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
        ASSERT_EQ(1, kprop_frameshift_local_stop_oracle(cdsE, sizeof cdsE, ctx.alt_cds,
                                                        ctx.alt_cds_len, 1u, 0u,
                                                        ctx.length_diff));
    }

    PASS();
}

/* Deterministic anchor for the general CodingContext frameshift fact. A net CDS length
 * change not divisible by three resolves to frameshift; a start-codon edit also carries
 * VEP's start_lost fact unless CDS_START_NF suppresses that predicate. */
TEST coding_context_delta_frameshift_known_scene(void) {
    /* M K P G F *  (start codon intact unless the edit lands in codon 1). */
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'T','A','A'
    };
    static const uint8_t ins1[2]  = { 'C','A' };  /* replace 1 base with 2 -> +1 */
    static const uint8_t del2a[1] = { 'C' };      /* replace 2 bases with 1 -> -1 */
    static const uint8_t incomplete_start_cds[18] = {
        'N','T','A',  'A','A','A',  'C','C','C',
        'G','G','G',  'T','T','T',  'T','A','A'
    };
    static const uint8_t insert_c[1] = { 'C' };
    static const uint8_t pre_cds[3] = { 'C','C','C' };
    static uint8_t repeated_start_cds[18] = {
        'A','T','G',  'G','T','A',  'C','G','T',
        'A','C','G',  'T','A','C',  'G','T','A'
    };
    static const uint8_t non_stop_terminal_cds[6] = {
        'A','T','G', 'C','C','C'
    };
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[40];
    uint8_t ref_pep[20];
    uint8_t alt_pep[20];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;

    /* (1) +1 frameshift at codon 3 (body), ATG intact -> frameshift at protein pos 3. */
    edit.cds_start = 7u; edit.ref = cds + 6u; edit.ref_len = 1u;
    edit.alt = ins1; edit.alt_len = 2u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_frameshift_at(&delta, 3));

    /* (2) -1 frameshift at codon 4 (body) -> frameshift at protein pos 4. */
    edit.cds_start = 10u; edit.ref = cds + 9u; edit.ref_len = 2u;
    edit.alt = del2a; edit.alt_len = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_frameshift_at(&delta, 4));

    /* (3) +1 frameshift disrupting the start codon (cds 1). */
    edit.cds_start = 1u; edit.ref = cds; edit.ref_len = 1u;
    edit.alt = ins1; edit.alt_len = 2u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ctx.pre_cds_bases = pre_cds;
    ctx.pre_cds_length = sizeof pre_cds;
    ctx.pre_cds_complete = 1u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && delta.start_lost &&
           delta.start_retained);

    /* (4) same start-disrupting edit but CDS_start_NF set: no start_lost to add, so the
     *     frameshift fact is emitted. */
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx,
                                                (uint64_t)DUCKVEP_TX_CDS_START_NF, &delta));
    ASSERT(delta.valid && delta.frameshift && !delta.start_lost &&
           !delta.start_retained);

    /* (5) Ensembl prefixes a CDS_START_NF sequence with synthetic N bases to
     * preserve its phase. VEP still calls a +1 insertion or -1 deletion in
     * that leading X codon a frameshift from CDS-coordinate lengths. Peptide-
     * dependent facts remain absent. Exercise both transcript strands. */
    {
        static const int8_t strands[2] = { 1, -1 };
        size_t strand_idx;

        for (strand_idx = 0u; strand_idx < 2u; strand_idx++) {
            edit.variant_strand = strands[strand_idx];
            edit.cds_start = 3u;
            edit.ref = NULL;
            edit.ref_len = 0u;
            edit.alt = insert_c;
            edit.alt_len = sizeof insert_c;
            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                      duckvep_coding_context_build(
                          incomplete_start_cds, sizeof incomplete_start_cds,
                          &edit_set, strands[strand_idx],
                          DUCKVEP_CODON_TABLE_STANDARD,
                          alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                          alt_pep, sizeof alt_pep, &ctx));
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                      duckvep_coding_context_delta_fill(
                          &ctx, (uint64_t)DUCKVEP_TX_CDS_START_NF, &delta));
            ASSERT(kprop_delta_is_frameshift_at(&delta, 1));

            edit.ref = incomplete_start_cds + 2u;
            edit.ref_len = 1u;
            edit.alt = NULL;
            edit.alt_len = 0u;
            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                      duckvep_coding_context_build(
                          incomplete_start_cds, sizeof incomplete_start_cds,
                          &edit_set, strands[strand_idx],
                          DUCKVEP_CODON_TABLE_STANDARD,
                          alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                          alt_pep, sizeof alt_pep, &ctx));
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                      duckvep_coding_context_delta_fill(
                          &ctx, (uint64_t)DUCKVEP_TX_CDS_START_NF, &delta));
            ASSERT(kprop_delta_is_frameshift_at(&delta, 1));
        }
    }

    /* (6) VEP-116 repeat-deletion witness, normalized from
     * chrDuck:121 TGGTAC>T. Deleting CDS bases 3..7 exposes the next G, so ATG
     * remains at the original start offset and start_retained_variant is true.
     * VEP's independent local codon/peptide comparison also calls start_lost.
     * Preserve both facts; neither is an inference from the other. */
    edit.variant_strand = 1;
    edit.cds_start = 3u;
    edit.ref = repeated_start_cds + 2u;
    edit.ref_len = 5u;
    edit.alt = NULL;
    edit.alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  repeated_start_cds, sizeof repeated_start_cds,
                  &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ctx.pre_cds_bases = pre_cds;
    ctx.pre_cds_length = sizeof pre_cds;
    ctx.pre_cds_complete = 1u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && delta.start_lost &&
           delta.start_retained);

    /* (7) VEP's _inv_start_altered helper is conditional on a 5-prime UTR.
     * Without pre-CDS bases, deleting the complete first codon is an ordinary
     * in-frame deletion and must not acquire start_lost from the now-empty
     * original start offset. */
    edit.variant_strand = 1;
    edit.cds_start = 1u;
    edit.ref = cds;
    edit.ref_len = 3u;
    edit.alt = NULL;
    edit.alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, sizeof cds, &edit_set, 1,
                  DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ctx.pre_cds_complete = 1u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.inframe_deletion && !delta.start_lost);

    /* (8) VEP's terminal-overlap gate is coordinate-based. A complete
     * annotated CDS without CDS_END_NF may end in a non-stop codon; when the
     * local alternate peptide is X, _ins_del_stop_altered still reconstructs
     * the original endpoint and reports stop_lost. Exercise deletion and
     * insertion because both coexist with the independent frameshift fact. */
    edit.cds_start = 6u;
    edit.ref = non_stop_terminal_cds + 5u;
    edit.ref_len = 1u;
    edit.alt = NULL;
    edit.alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  non_stop_terminal_cds, sizeof non_stop_terminal_cds,
                  &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ctx.post_cds_complete = 1u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && delta.stop_lost &&
           !delta.stop_retained);

    edit.ref = NULL;
    edit.ref_len = 0u;
    edit.alt = insert_c;
    edit.alt_len = sizeof insert_c;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  non_stop_terminal_cds, sizeof non_stop_terminal_cds,
                  &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ctx.post_cds_complete = 1u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && delta.stop_lost &&
           !delta.stop_retained);

    PASS();
}

/* Ensembl CDS_END_NF transcripts legitimately end with one or two untranslated
 * CDS bases. VEP still evaluates ordinary length-changing edits against the
 * complete-codon peptide plus its synthetic trailing X; the incomplete tail is
 * not a reason to collapse an otherwise simple edit into an unresolved result. */
TEST coding_context_delta_cds_end_nf_known_scene(void) {
    /* COVERAGE_WITNESS: HGVSp replay coverage/vep_stop_equal */
    static uint8_t cds_mod1[13] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T'
    };
    static uint8_t cds_mod2[14] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T'
    };
    static const uint8_t replace_c_ca[2] = { 'C','A' };
    static const uint8_t insert_gcc[3] = { 'G','C','C' };
    static const uint8_t replace_t_a[1] = { 'A' };
    static const uint8_t replace_t_ta[2] = { 'T','A' };
    static const uint8_t insert_gcca[4] = { 'G','C','C','A' };
    static const uint8_t insert_tag[3] = { 'T','A','G' };
    static const uint8_t insert_taga[4] = { 'T','A','G','A' };
    static const uint8_t insert_agt[3] = { 'A','G','T' };
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[32];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    duckvep_coding_context_t ctx;
    duckvep_coding_peptide_window_t peptide_window;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t protein_fact;
    char protein_hgvs[32];
    size_t protein_hgvs_required = 0u;
    size_t strand_idx;
    static const int8_t strands[2] = { 1, -1 };

    memset(&edit, 0, sizeof edit);
    edit_set.edits = &edit;
    edit_set.count = 1u;

    for (strand_idx = 0u; strand_idx < 2u; strand_idx++) {
        int8_t strand = strands[strand_idx];

        /* A +1 body edit remains a frameshift for both possible incomplete-tail
         * lengths and on both genomic strands. */
        edit.cds_start = 7u;
        edit.ref = cds_mod1 + 6u;
        edit.ref_len = 1u;
        edit.alt = replace_c_ca;
        edit.alt_len = sizeof replace_c_ca;
        edit.variant_strand = strand;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod1, sizeof cds_mod1, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.inframe_insertion &&
               !delta.inframe_deletion);

        edit.ref = cds_mod2 + 6u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.inframe_insertion &&
               !delta.inframe_deletion);

        /* Complete-codon insertion/deletion facts are likewise independent of
         * the trailing partial codon. */
        edit.cds_start = 7u;
        edit.ref = NULL;
        edit.ref_len = 0u;
        edit.alt = insert_gcc;
        edit.alt_len = sizeof insert_gcc;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod1, sizeof cds_mod1, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && !delta.frameshift);

        edit.ref = cds_mod2 + 6u;
        edit.ref_len = 3u;
        edit.alt = NULL;
        edit.alt_len = 0u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.inframe_deletion && !delta.frameshift);

        /* The affected rounded codon itself reaches the partial tail. VEP's
         * local reference peptide is one complete residue plus synthetic X;
         * deleting the three interior bases remains an in-frame deletion. */
        edit.cds_start = 11u;
        edit.ref = cds_mod2 + 10u;
        edit.ref_len = 3u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.inframe_deletion && !delta.frameshift);
        ASSERT(!delta.partial_codon);

        /* VEP's predicate is tied to the first affected peptide coordinate.
         * An equal-length edit beginning in the final two-base codon yields
         * partial-codon plus coding-unknown, never missense. */
        edit.cds_start = 13u;
        edit.ref = cds_mod2 + 12u;
        edit.ref_len = 1u;
        edit.alt = replace_t_a;
        edit.alt_len = sizeof replace_t_a;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.coding_unknown &&
               !delta.synonymous && !delta.missense && !delta.frameshift);

        /* A +1 replacement still is not a frameshift when its VEP peptide
         * coordinate is the incomplete terminal codon. */
        edit.alt = replace_t_ta;
        edit.alt_len = sizeof replace_t_ta;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.coding_unknown &&
               !delta.frameshift && !delta.inframe_deletion);

        /* VEP does not guard inframe_insertion with partial_codon. A four-base
         * insertion leaves a synthetic X and therefore coexists with
         * coding_unknown. */
        edit.ref = NULL;
        edit.ref_len = 0u;
        edit.alt = insert_gcca;
        edit.alt_len = sizeof insert_gcca;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.inframe_insertion &&
               delta.coding_unknown && !delta.frameshift);

        /* ClinVar 1:45013701:C:CTAG on ENST00000650713 inserts before CDS
         * base 280, the first base of its partial codon (not between bases
         * 280 and 281). VEP's consequence codons are "-/TAG" and peptides
         * empty reference peptide versus stop. The smaller CDS preserves
         * that codon-start placement. */
        edit.cds_start = 13u;
        edit.alt = insert_tag;
        edit.alt_len = sizeof insert_tag;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT(duckvep_coding_context_peptide_window_open(
            &ctx, &peptide_window));
        ASSERT_EQ(0u, peptide_window.ref_length);
        ASSERT_EQ(1u, peptide_window.alt_length);
        ASSERT_EQ((uint8_t)'*',
                  duckvep_coding_context_peptide_window_base(
                      &ctx, &peptide_window, 1, 0u));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid);
        ASSERT(delta.partial_codon);
        ASSERT(delta.inframe_insertion);
        ASSERT(delta.stop_gained);
        ASSERT(!delta.coding_unknown);
        ASSERT(!delta.frameshift);

        /* Retain the distinct internal site: pinned VEP instead spells
         * "tt/tTAGt" and "X/LX", with no stop_gained. */
        edit.cds_start = 14u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT(duckvep_coding_context_peptide_window_open(&ctx, &peptide_window));
        ASSERT_EQ(1u, peptide_window.ref_length);
        ASSERT_EQ(2u, peptide_window.alt_length);
        ASSERT_EQ((uint8_t)'L', duckvep_coding_context_peptide_window_base(
            &ctx, &peptide_window, 1, 0u));
        ASSERT_EQ((uint8_t)'X', duckvep_coding_context_peptide_window_base(
            &ctx, &peptide_window, 1, 1u));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.inframe_insertion &&
               delta.coding_unknown && !delta.stop_gained && !delta.frameshift);

        /* An internal TAGA payload gives "tt/tTAGAt" and "X/LD":
         * protein_altering replaces the insertion-only stop prediction. */
        edit.alt = insert_taga;
        edit.alt_len = sizeof insert_taga;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT(duckvep_coding_context_peptide_window_open(
            &ctx, &peptide_window));
        ASSERT_EQ(1u, peptide_window.ref_length);
        ASSERT_EQ(2u, peptide_window.alt_whole_length);
        ASSERT_EQ(0u, peptide_window.alt_partial_x);
        ASSERT_EQ(2u, peptide_window.alt_length);
        ASSERT_EQ((uint8_t)'L',
                  duckvep_coding_context_peptide_window_base(
                      &ctx, &peptide_window, 1, 0u));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon &&
               !delta.stop_gained && !delta.inframe_insertion &&
               !delta.coding_unknown && !delta.frameshift &&
               !delta.inframe_deletion && delta.protein_altering);

        /* A 3-prime UTR participates in TVA's ALT codon slice, independently
         * of the shorter reference slice. tt/tAt becomes tt/tAta when the
         * first UTR base is A: X/YX, not X/Y. Physical CDS APIs still stop
         * at the edited CDS endpoint. */
        edit.alt = replace_t_a;
        edit.alt_len = sizeof replace_t_a;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ctx.post_cds_bases = replace_t_a;
        ctx.post_cds_length = sizeof replace_t_a;
        ctx.post_cds_complete = 1u;
        ASSERT(duckvep_coding_context_peptide_window_open(&ctx, &peptide_window));
        ASSERT_EQ(2u, peptide_window.ref_nt_length);
        ASSERT_EQ(4u, peptide_window.alt_nt_length);
        ASSERT_EQ(2u, peptide_window.alt_length);
        ASSERT_EQ((uint8_t)'Y', duckvep_coding_context_peptide_window_base(
            &ctx, &peptide_window, 1, 0u));
        ASSERT_EQ((uint8_t)'X', duckvep_coding_context_peptide_window_base(
            &ctx, &peptide_window, 1, 1u));
        ASSERT_EQ('\0', duckvep_coding_context_cds_base(&ctx, 1, ctx.alt_cds_len));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK, duckvep_coding_context_delta_fill(
            &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.inframe_insertion &&
               delta.coding_unknown && !delta.stop_gained && !delta.protein_altering);
        ctx.post_cds_length = SIZE_MAX;
        ASSERT(!duckvep_coding_context_peptide_window_open(&ctx, &peptide_window));
        ctx.post_cds_length = 1u;
        ctx.post_cds_bases = NULL;
        ASSERT(!duckvep_coding_context_peptide_window_open(&ctx, &peptide_window));

        /* HGVS 3-prime placement rotates the same insertion to AGT. VEP's
         * protein path independently rounds the edited CDS from the preceding
         * partial base, sees TAG, normalizes reference X against alternate
         * stop, and renders equality even though the unshifted consequence
         * path above emitted stop_gained. */
        edit.alt = insert_agt;
        edit.alt_len = sizeof insert_agt;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.stop_gained &&
               delta.coding_unknown && !delta.inframe_insertion);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_fact_build(
                      &ctx, &delta, &protein_fact));
        ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_EQUAL, protein_fact.shape);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_render(
                      &protein_fact, 0, protein_hgvs, sizeof protein_hgvs,
                      &protein_hgvs_required));
        ASSERT_EQ(0, strcmp("p.Ter5=", protein_hgvs));

        /* Strict mode keeps the internal site's X reference residue instead
         * of VEP's Xaa-to-Ter normalization. It is a substitution, not the
         * empty-reference insertion previously inferred at this site. */
        ctx.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_fact_build(
                      &ctx, &delta, &protein_fact));
        ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_SUBSTITUTION, protein_fact.shape);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_render(
                      &protein_fact, 0, protein_hgvs, sizeof protein_hgvs,
                      &protein_hgvs_required));
        ASSERT_EQ(0, strcmp("p.Xaa5Ter", protein_hgvs));

        /* Removing the two-base partial codon is protein-altering in VEP, not
         * a frameshift or an in-frame deletion. */
        edit.cds_start = 13u;
        edit.ref = cds_mod2 + 12u;
        edit.ref_len = 2u;
        edit.alt = NULL;
        edit.alt_len = 0u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      cds_mod2, sizeof cds_mod2, &edit_set, strand,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(
                      &ctx, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta));
        ASSERT(delta.valid && delta.partial_codon && delta.protein_altering &&
               !delta.frameshift && !delta.inframe_deletion &&
               !delta.coding_unknown);
    }

    PASS();
}

#define KPROP_PARTIAL_INSERTION_MAX 9u

struct kprop_partial_terminal_insertion {
    uint8_t tail_length;
    uint8_t site_offset;
    uint8_t insertion_length;
    uint8_t codon_table;
    int8_t transcript_strand;
    uint8_t genomic_alt[KPROP_PARTIAL_INSERTION_MAX];
};

static enum theft_alloc_res kprop_partial_terminal_insertion_alloc(
    struct theft *t, void *env, void **instance) {

    static const uint8_t bases[4] = {'A', 'C', 'G', 'T'};
    struct kprop_partial_terminal_insertion *s;
    size_t i;
    (void)env;

    s = (struct kprop_partial_terminal_insertion *)calloc(1u, sizeof *s);
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->tail_length = (uint8_t)(1u + kprop_bounded(t, 2u));
    s->site_offset = (uint8_t)kprop_bounded(
        t, (uint64_t)s->tail_length + 1u);
    s->insertion_length = (uint8_t)(
        1u + kprop_bounded(t, KPROP_PARTIAL_INSERTION_MAX));
    s->codon_table = (uint8_t)(kprop_bounded(t, 2u) == 0u
        ? DUCKVEP_CODON_TABLE_STANDARD : DUCKVEP_CODON_TABLE_VERT_MITO);
    s->transcript_strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    for (i = 0u; i < (size_t)s->insertion_length; i++) {
        s->genomic_alt[i] = bases[kprop_bounded(t, 4u)];
    }
    *instance = s;
    return THEFT_ALLOC_OK;
}

static void kprop_partial_terminal_insertion_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static struct theft_type_info kprop_partial_terminal_insertion_info = {
    .alloc = kprop_partial_terminal_insertion_alloc,
    .free = kprop_partial_terminal_insertion_free,
};

static struct {
    uint32_t tail1;
    uint32_t tail2;
    uint32_t site_first;
    uint32_t site_internal;
    uint32_t after_tail_rejected;
    uint32_t length_mod0;
    uint32_t length_mod1;
    uint32_t length_mod2;
    uint32_t standard;
    uint32_t mitochondrial;
    uint32_t same_orientation;
    uint32_t reverse_orientation;
    uint32_t stop;
    uint32_t nonstop;
} g_partial_terminal_insertion_cov;

static uint8_t kprop_partial_terminal_oriented_base(
    const struct kprop_partial_terminal_insertion *s, size_t index) {

    if (s->transcript_strand > 0) return s->genomic_alt[index];
    return (uint8_t)kprop_complement_base(
        (char)s->genomic_alt[(size_t)s->insertion_length - 1u - index]);
}

static enum theft_trial_res prop_partial_terminal_insertion_strata(
    struct theft *t, void *arg1) {

    static const uint8_t cds_mod1[13] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T'
    };
    static const uint8_t cds_mod2[14] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T'
    };
    const struct kprop_partial_terminal_insertion *s =
        (const struct kprop_partial_terminal_insertion *)arg1;
    const uint8_t *cds = s->tail_length == 1u ? cds_mod1 : cds_mod2;
    size_t cds_length = s->tail_length == 1u ? sizeof cds_mod1 : sizeof cds_mod2;
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    duckvep_coding_context_t context;
    duckvep_coding_peptide_window_t window;
    duckvep_sequence_delta_t delta;
    duckvep_context_delta_status_t delta_status;
    uint8_t alt_cds[32];
    uint8_t ref_peptide[16];
    uint8_t alt_peptide[16];
    uint8_t expected_peptide[8];
    size_t expected_ref_nt = s->site_offset == 0u ? 0u : s->tail_length;
    size_t expected_alt_nt = expected_ref_nt + s->insertion_length;
    size_t expected_whole = expected_alt_nt / 3u;
    size_t expected_ref_length = expected_ref_nt != 0u ? 1u : 0u;
    size_t expected_alt_length;
    size_t i;
    int saw_stop = 0;
    int expected_insertion;
    int expected_protein_altering;
    int expected_unknown;
    uint8_t expected_partial_x;
    (void)t;

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 13u + (uint32_t)s->site_offset;
    edit.alt = s->genomic_alt;
    edit.alt_len = s->insertion_length;
    edit.variant_strand = 1;
    edit_set.edits = &edit;
    edit_set.count = 1u;

    if (duckvep_coding_context_build(
            cds, cds_length, &edit_set, s->transcript_strand,
            (duckvep_codon_table_t)s->codon_table,
            alt_cds, sizeof alt_cds,
            ref_peptide, sizeof ref_peptide,
            alt_peptide, sizeof alt_peptide, &context) !=
        DUCKVEP_CODING_CONTEXT_OK) {
        fprintf(stderr,
                "\n[terminal-partial failure] build tail=%u site=%u len=%u "
                "table=%u strand=%d\n",
                s->tail_length, s->site_offset, s->insertion_length,
                s->codon_table, s->transcript_strand);
        return THEFT_TRIAL_FAIL;
    }
    if (!duckvep_coding_context_peptide_window_open(&context, &window) ||
        window.ref_nt_length != expected_ref_nt ||
        window.ref_length != expected_ref_length ||
        window.alt_nt_length != expected_alt_nt ||
        window.alt_whole_length != expected_whole) {
        fprintf(stderr,
                "\n[terminal-partial failure] window tail=%u site=%u len=%u "
                "table=%u strand=%d ref_nt=%zu ref=%zu alt_nt=%zu whole=%zu\n",
                s->tail_length, s->site_offset, s->insertion_length,
                s->codon_table, s->transcript_strand,
                window.ref_nt_length, window.ref_length,
                window.alt_nt_length, window.alt_whole_length);
        return THEFT_TRIAL_FAIL;
    }
    /* TranscriptMapper::genomic2pep rounds the reversed insertion endpoints
     * independently. Only a codon-start insertion has zero reference codons.
     * For internal sites, splice the oriented payload between the retained T
     * bases; this oracle does not read the implementation's alternate CDS. */
    for (i = 0u; i < expected_whole; i++) {
        char codon[4];
        uint8_t expected;
        size_t b;
        for (b = 0u; b < 3u; b++) {
            size_t position = i * 3u + b;
            if (s->site_offset == 0u) {
                codon[b] = (char)kprop_partial_terminal_oriented_base(s, position);
            } else if (position < s->site_offset ||
                       position >= s->site_offset + s->insertion_length) {
                codon[b] = 'T';
            } else {
                codon[b] = (char)kprop_partial_terminal_oriented_base(
                    s, position - s->site_offset);
            }
        }
        codon[3] = '\0';
        expected = (uint8_t)duckvep_translate_codon(
            codon, (duckvep_codon_table_t)s->codon_table);
        if (duckvep_coding_context_peptide_window_base(
                &context, &window, 1, i) != expected) {
            fprintf(stderr,
                    "\n[terminal-partial failure] translation tail=%u site=%u "
                    "len=%u table=%u strand=%d codon=%s got=%c want=%c\n",
                    s->tail_length, s->site_offset, s->insertion_length,
                    s->codon_table, s->transcript_strand, codon,
                    duckvep_coding_context_peptide_window_base(
                        &context, &window, 1, i), expected);
            return THEFT_TRIAL_FAIL;
        }
        expected_peptide[i] = expected;
        if (expected == (uint8_t)'*') saw_stop = 1;
    }
    expected_partial_x = (uint8_t)(
        (expected_alt_nt % 3u) != 0u &&
        !(expected_whole == 1u && expected_peptide[0] == (uint8_t)'*'));
    expected_alt_length = expected_whole + expected_partial_x;
    if (expected_partial_x) expected_peptide[expected_whole] = (uint8_t)'X';
    if (window.alt_partial_x != expected_partial_x ||
        window.alt_length != expected_alt_length) {
        fprintf(stderr,
                "\n[terminal-partial failure] peptide-shape tail=%u site=%u "
                "len=%u table=%u strand=%d whole=%zu partial_x=%u/%u "
                "alt_length=%zu\n",
                s->tail_length, s->site_offset, s->insertion_length,
                s->codon_table, s->transcript_strand,
                window.alt_whole_length, window.alt_partial_x,
                expected_partial_x, window.alt_length);
        return THEFT_TRIAL_FAIL;
    }
    expected_insertion = expected_ref_length == 0u ||
        (!saw_stop && expected_partial_x);
    expected_protein_altering = expected_ref_length != expected_alt_length &&
        expected_ref_length != 0u && expected_peptide[0] != (uint8_t)'*' &&
        !expected_partial_x;
    expected_unknown = !expected_protein_altering &&
        (expected_ref_length != 0u || expected_partial_x);
    delta_status = duckvep_coding_context_delta_fill(
        &context, (uint64_t)DUCKVEP_TX_CDS_END_NF, &delta);
    if (s->site_offset == s->tail_length) {
        if (delta_status != DUCKVEP_CONTEXT_DELTA_UNSUPPORTED || delta.valid) {
            fprintf(stderr,
                    "\n[terminal-partial failure] after-tail tail=%u len=%u "
                    "table=%u strand=%d status=%d valid=%u\n",
                    s->tail_length, s->insertion_length, s->codon_table,
                    s->transcript_strand, (int)delta_status, delta.valid);
            return THEFT_TRIAL_FAIL;
        }
    } else if (delta_status != DUCKVEP_CONTEXT_DELTA_OK ||
               !delta.valid || !delta.partial_codon || delta.frameshift ||
               delta.inframe_deletion ||
               delta.inframe_insertion != (uint8_t)expected_insertion ||
               delta.stop_gained != (uint8_t)saw_stop ||
               delta.coding_unknown != (uint8_t)expected_unknown ||
               delta.stop_lost || delta.stop_retained || delta.start_lost ||
               delta.start_retained ||
               delta.protein_altering != (uint8_t)expected_protein_altering ||
               delta.synonymous || delta.missense ||
               delta.ref_aa != (uint8_t)(
                   !expected_insertion && !expected_protein_altering &&
                   expected_ref_length == 1u ? 'X' : 0u) ||
               delta.alt_aa != (uint8_t)(
                   !expected_insertion && !expected_protein_altering &&
                   expected_alt_length == 1u ? expected_peptide[0] : 0u) ||
               delta.protein_pos != (int32_t)(window.ref_peptide_offset + 1u)) {
        fprintf(stderr,
                "\n[terminal-partial failure] delta tail=%u site=%u len=%u "
                "table=%u strand=%d status=%d valid=%u partial=%u frameshift=%u "
                "inframe_del=%u inframe_ins=%u stop_gained=%u/%d "
                "coding_unknown=%u/%u protein_altering=%u protein_pos=%d/%zu\n",
                s->tail_length, s->site_offset, s->insertion_length,
                s->codon_table, s->transcript_strand, (int)delta_status,
                delta.valid,
                delta.partial_codon, delta.frameshift,
                delta.inframe_deletion, delta.inframe_insertion,
                delta.stop_gained, saw_stop, delta.coding_unknown,
                (unsigned)expected_unknown, delta.protein_altering,
                delta.protein_pos, window.ref_peptide_offset + 1u);
        return THEFT_TRIAL_FAIL;
    }

    if (s->tail_length == 1u) g_partial_terminal_insertion_cov.tail1++;
    else g_partial_terminal_insertion_cov.tail2++;
    if (s->site_offset == 0u) g_partial_terminal_insertion_cov.site_first++;
    else if (s->site_offset == s->tail_length)
        g_partial_terminal_insertion_cov.after_tail_rejected++;
    else g_partial_terminal_insertion_cov.site_internal++;
    if (s->insertion_length % 3u == 0u)
        g_partial_terminal_insertion_cov.length_mod0++;
    else if (s->insertion_length % 3u == 1u)
        g_partial_terminal_insertion_cov.length_mod1++;
    else g_partial_terminal_insertion_cov.length_mod2++;
    if (s->codon_table == (uint8_t)DUCKVEP_CODON_TABLE_STANDARD)
        g_partial_terminal_insertion_cov.standard++;
    else g_partial_terminal_insertion_cov.mitochondrial++;
    if (s->transcript_strand > 0)
        g_partial_terminal_insertion_cov.same_orientation++;
    else g_partial_terminal_insertion_cov.reverse_orientation++;
    if (window.alt_whole_length != 0u) {
        if (saw_stop) g_partial_terminal_insertion_cov.stop++;
        else g_partial_terminal_insertion_cov.nonstop++;
    }
    return THEFT_TRIAL_PASS;
}

TEST partial_terminal_insertion_covers_generated_vep_strata(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    memset(&g_partial_terminal_insertion_cov, 0,
           sizeof g_partial_terminal_insertion_cov);
    cfg.name =
        "terminal partial-codon insertion == codon-rounded VEP translation oracle";
    cfg.prop1 = prop_partial_terminal_insertion_strata;
    cfg.type_info[0] = &kprop_partial_terminal_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_partial_terminal_insertion_cov.tail1 > 0u);
    ASSERT(g_partial_terminal_insertion_cov.tail2 > 0u);
    ASSERT(g_partial_terminal_insertion_cov.site_first > 0u);
    ASSERT(g_partial_terminal_insertion_cov.site_internal > 0u);
    ASSERT(g_partial_terminal_insertion_cov.after_tail_rejected > 0u);
    ASSERT(g_partial_terminal_insertion_cov.length_mod0 > 0u);
    ASSERT(g_partial_terminal_insertion_cov.length_mod1 > 0u);
    ASSERT(g_partial_terminal_insertion_cov.length_mod2 > 0u);
    ASSERT(g_partial_terminal_insertion_cov.standard > 0u);
    ASSERT(g_partial_terminal_insertion_cov.mitochondrial > 0u);
    ASSERT(g_partial_terminal_insertion_cov.same_orientation > 0u);
    ASSERT(g_partial_terminal_insertion_cov.reverse_orientation > 0u);
    ASSERT(g_partial_terminal_insertion_cov.stop > 0u);
    ASSERT(g_partial_terminal_insertion_cov.nonstop > 0u);
    fprintf(stderr,
            "[terminal-partial-insertion coverage] tail1=%u tail2=%u "
            "site_first=%u site_internal=%u after_tail_rejected=%u "
            "length_mod0=%u length_mod1=%u length_mod2=%u "
            "standard=%u mitochondrial=%u same_orientation=%u "
            "reverse_orientation=%u stop=%u nonstop=%u\n",
            g_partial_terminal_insertion_cov.tail1,
            g_partial_terminal_insertion_cov.tail2,
            g_partial_terminal_insertion_cov.site_first,
            g_partial_terminal_insertion_cov.site_internal,
            g_partial_terminal_insertion_cov.after_tail_rejected,
            g_partial_terminal_insertion_cov.length_mod0,
            g_partial_terminal_insertion_cov.length_mod1,
            g_partial_terminal_insertion_cov.length_mod2,
            g_partial_terminal_insertion_cov.standard,
            g_partial_terminal_insertion_cov.mitochondrial,
            g_partial_terminal_insertion_cov.same_orientation,
            g_partial_terminal_insertion_cov.reverse_orientation,
            g_partial_terminal_insertion_cov.stop,
            g_partial_terminal_insertion_cov.nonstop);
    PASS();
}

/* VEP-116 states discovered independently by the statistical differentials.
 * They share one cause: stop, frame, insertion, and protein-shape predicates
 * inspect the same local peptide but remain independently true. */
TEST coding_context_delta_terminal_combinations_known_scene(void) {
    static const uint8_t cds[] =
        "ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA";
    static const uint8_t protein_alt[] =
        "GCGTTATACCGTATACCCACGAGTTACAGTGAACCTAGGC";
    static const uint8_t lost_alt[] =
        "GCCACCACTGGGCCTAAATCGACG";
    static const uint8_t retained_alt[] =
        "GCGAGCGGTACTTGTTCTTGTCTCGCTTTGGGGTGCCACTTGAACAT";
    static const uint8_t retained_before_w[] =
        "TGGTAGCCATTATTGGGCTGTGCTCACCCTAGTATTGCG";
    static const uint8_t retained_inside_w[] =
        "GGTAGCCGCCAGCGGCCGGTCAGAATC";
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[160];
    uint8_t ref_pep[80];
    uint8_t alt_pep[80];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    size_t cds_len = sizeof cds - 1u;

    ASSERT_EQ(72u, cds_len);
    memset(&edit, 0, sizeof edit);
    edit.variant_strand = 1;
    edit_set.edits = &edit;
    edit_set.count = 1u;

    /* chrDuck:233 ACTGGT>A. The alternate local peptide is exactly "*"
     * even though the allele changes frame, so VEP emits frameshift alone. */
    edit.cds_start = 66u;
    edit.ref = cds + 65u;
    edit.ref_len = 5u;
    edit.alt = NULL;
    edit.alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, cds_len, &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && !delta.stop_lost &&
           !delta.stop_retained && !delta.inframe_insertion &&
           !delta.protein_altering);

    /* chrDuck:206 long delins. The terminal stop remains at the same local
     * index, suppressing frameshift, while the altered peptide shape remains a
     * separate protein_altering_variant. */
    edit.cds_start = 41u;
    edit.ref = cds + 40u;
    edit.ref_len = 30u;
    edit.alt = protein_alt;
    edit.alt_len = (uint32_t)(sizeof protein_alt - 1u);
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, cds_len, &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.stop_retained && delta.protein_altering &&
           !delta.frameshift && !delta.stop_lost &&
           !delta.inframe_insertion);

    /* chrDuck:236 GGTAA>GGT... The reference local peptide starts with the
     * stop, suppressing frameshift; the trimmed alternate still ends in that
     * stop, so inframe_insertion and stop_lost coexist. */
    edit.cds_start = 71u;
    edit.ref = cds + 70u;
    edit.ref_len = 2u;
    edit.alt = lost_alt;
    edit.alt_len = (uint32_t)(sizeof lost_alt - 1u);
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, cds_len, &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.inframe_insertion && delta.stop_lost &&
           !delta.frameshift && !delta.stop_retained &&
           !delta.protein_altering);

    /* chrDuck:239 AA>A... Both local peptides start with the stop at index
     * zero. The longer alternate retains it and remains insertion-shaped. */
    edit.cds_start = 72u;
    edit.ref = cds + 71u;
    edit.ref_len = 1u;
    edit.alt = retained_alt;
    edit.alt_len = (uint32_t)(sizeof retained_alt - 1u);
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, cds_len, &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
           !delta.frameshift && !delta.stop_lost &&
           !delta.protein_altering);

    /* Held-out seed 197: a pure insertion immediately before the final TGG
     * replaces the end of Transcript::translate with an identical W prefix and
     * places '*' at the first position beyond that reference peptide. VEP's
     * ref_eq_alt_sequence therefore calls the original stop retained. */
    edit.cds_start = 67u;
    edit.ref = NULL;
    edit.ref_len = 0u;
    edit.alt = retained_before_w;
    edit.alt_len = (uint32_t)(sizeof retained_before_w - 1u);
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, cds_len, &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
           !delta.stop_gained && !delta.stop_lost && !delta.frameshift);

    /* The insertion one base later widens the local codon but reaches the same
     * full-peptide state: the reference W prefix survives and the next residue
     * is the retained stop. */
    edit.cds_start = 68u;
    edit.alt = retained_inside_w;
    edit.alt_len = (uint32_t)(sizeof retained_inside_w - 1u);
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(
                  cds, cds_len, &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                  alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                  alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
           !delta.stop_gained && !delta.stop_lost && !delta.frameshift);

    {
        static const uint8_t short_terminal_cds[] = {
            'A','T','G', 'T','A','A'
        };
        uint32_t terminal_edit_position1 = 2u;
        uint8_t terminal_edit_alt = (uint8_t)'X';

        /* VEP's direct reference-peptide predicates see Translation SeqEdits,
         * but its X/undefined length-change fallback does not: the CIL pair
         * checks genomic terminal coordinates, edits raw translateable DNA,
         * and translates that raw endpoint. Deleting one base from TAA is
         * therefore stop_lost even when a SeqEdit renders the local reference
         * residue X. Do not collapse these two upstream authorities. */
        edit.cds_start = 6u;
        edit.ref = short_terminal_cds + 5u;
        edit.ref_len = 1u;
        edit.alt = NULL;
        edit.alt_len = 0u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      short_terminal_cds, sizeof short_terminal_cds,
                      &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                      alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                      alt_pep, sizeof alt_pep, &ctx));
        ctx.ref_peptide_edit_position1 = &terminal_edit_position1;
        ctx.ref_peptide_edit_alt = &terminal_edit_alt;
        ctx.ref_peptide_edit_count = 1u;
        ctx.post_cds_complete = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_lost && !delta.stop_retained);
    }

    {
        static const uint8_t reverse_cil_cds[] = {
            'A','T','G', 'A','A','A', 'A','T','G', 'C','C','C', 'T','G','A'
        };
        static const uint8_t genomic_insert[] = {
            'A','T','C','A','G','C','C','T'
        };

        /* VEP 116's X-peptide fallback uses _overlaps_stop_codon_cil: on a
         * reverse-strand transcript it extends the minimized insertion by the
         * inserted genomic length before testing the terminal codon. The live
         * _ins_del_stop_altered_cil path applies the insertion to the CDS and
         * retranslates the original endpoint. These bases preserve that stop,
         * so M/IG*X is stop-retained plus protein-altering rather than
         * frameshift plus stop-gained. This is the reduced witness for GRCh38
         * 7:148807690 C>CATCAGCCT / ENST00001066230. */
        edit.cds_start = 9u;
        edit.ref = NULL;
        edit.ref_len = 0u;
        edit.alt = genomic_insert;
        edit.alt_len = sizeof genomic_insert;
        edit.variant_strand = 1;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      reverse_cil_cds, sizeof reverse_cil_cds,
                      &edit_set, -1, DUCKVEP_CODON_TABLE_STANDARD,
                      alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                      alt_pep, sizeof alt_pep, &ctx));
        ctx.insertion_length_reaches_terminal_stop = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained &&
               delta.protein_altering && !delta.frameshift &&
               !delta.stop_gained && !delta.stop_lost &&
               !delta.inframe_insertion);
    }

    {
        static const uint8_t concrete_stop_cds[] = {
            'A','T','G', 'G','C','A', 'T','G','G', 'A','G','T'
        };
        static const uint8_t inserted_at[] = { 'A', 'T' };

        /* Held-out seed 20260719: the insertion length reaches the original
         * CDS endpoint, and retranslating that endpoint still yields '*'.
         * VEP does not consult its CIL fallback because the codon-local
         * alternate peptide is the concrete value "*" with no X.
         * ref_eq_alt_sequence is false for local W versus local *, so
         * stop_gained remains true
         * and stop_retained does not suppress the frameshift. */
        edit.cds_start = 9u;
        edit.ref = NULL;
        edit.ref_len = 0u;
        edit.alt = inserted_at;
        edit.alt_len = sizeof inserted_at;
        edit.variant_strand = 1;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      concrete_stop_cds, sizeof concrete_stop_cds,
                      &edit_set, -1, DUCKVEP_CODON_TABLE_STANDARD,
                      alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                      alt_pep, sizeof alt_pep, &ctx));
        ctx.insertion_length_reaches_terminal_stop = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && delta.frameshift &&
               !delta.stop_retained && !delta.stop_lost &&
               !delta.inframe_insertion && !delta.protein_altering);
    }

    PASS();
}

/* Randomized frame-arithmetic oracle: over random valid CDS + single body edits, the
 * general CodingContext emits the frameshift fact exactly when the net length change is
 * not divisible by three, and never flags an in-frame length change as a frameshift. The
 * oracle is the frame arithmetic itself, independent of the classifier implementation. */
TEST coding_context_delta_frameshift_matches_length_oracle(void) {
    static const char bases[4] = { 'A', 'C', 'G', 'T' };
    uint64_t rng = kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED) ^ UINT64_C(0xF00DFACE);
    uint64_t trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    uint64_t t;
    unsigned fs_seen = 0u;
    unsigned inframe_seen = 0u;
    unsigned fs_stop_seen = 0u;

    for (t = 0u; t < trials; t++) {
        uint8_t cds[36];
        uint8_t alt_bytes[4];
        uint8_t alt_cds[64];
        uint8_t ref_pep[32];
        uint8_t alt_pep[32];
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;
        size_t ncodon, cds_len, i;
        uint32_t codon_idx, cds_start;
        unsigned ref_len, alt_len;
        int net;

        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        ncodon = 7u + (size_t)((rng >> 33) % 6u); /* 7..12 codons */
        cds_len = ncodon * 3u;
        cds[0] = 'A'; cds[1] = 'T'; cds[2] = 'G';
        for (i = 3u; i < cds_len; i++) {
            rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
            cds[i] = (uint8_t)bases[(rng >> 33) % 4u];
        }
        /* Ordinary coding transcripts do not carry an unannotated internal
         * stop. VEP's frameshift predicate deliberately returns false when the
         * first affected reference peptide begins with '*', so that state is a
         * separate predicate test rather than part of this frame-arithmetic
         * property. */
        for (i = 3u; i + 2u < cds_len; i += 3u) {
            if (duckvep_translate_codon(
                    (const char *)(cds + i),
                    DUCKVEP_CODON_TABLE_STANDARD) == '*') {
                cds[i + 1u] = (uint8_t)'C';
            }
        }
        /* Edit in a body codon: not codon 1, and leaving >= 2 codons of tail. */
        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        codon_idx = 1u + (uint32_t)((rng >> 33) % (ncodon - 3u)); /* 1..ncodon-3 */
        cds_start = codon_idx * 3u + 1u;
        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        ref_len = 1u + (unsigned)((rng >> 33) % 4u); /* 1..4 */
        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        alt_len = 1u + (unsigned)((rng >> 33) % 4u); /* 1..4 */
        for (i = 0u; i < alt_len; i++) {
            rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
            alt_bytes[i] = (uint8_t)bases[(rng >> 33) % 4u];
        }

        memset(&edit, 0, sizeof edit);
        edit.variant_strand = 1; edit.cds_start = cds_start;
        edit.ref = cds + (cds_start - 1u); edit.ref_len = (uint16_t)ref_len;
        edit.alt = alt_bytes; edit.alt_len = (uint16_t)alt_len;
        edit_set.edits = &edit; edit_set.count = 1u;

        if (duckvep_coding_context_build(cds, cds_len, &edit_set, 1,
                                         DUCKVEP_CODON_TABLE_STANDARD,
                                         alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                         alt_pep, sizeof alt_pep, &ctx) !=
            DUCKVEP_CODING_CONTEXT_OK) {
            continue;
        }
        (void)duckvep_coding_context_delta_fill(&ctx, 0u, &delta);
        net = (int)alt_len - (int)ref_len;
        if (net != 0 && (net % 3) != 0) {
            int stop_oracle;
            ASSERT(delta.valid && delta.frameshift);
            ASSERT(!delta.inframe_insertion && !delta.inframe_deletion &&
                   !delta.missense && !delta.synonymous && !delta.start_lost);
            /* stop_gained must equal VEP's independent local-window oracle, never the
             * downstream retranslation the full alt peptide would suggest. */
            stop_oracle = kprop_frameshift_local_stop_oracle(cds, cds_len, ctx.alt_cds,
                                                             ctx.alt_cds_len, cds_start,
                                                             ref_len, ctx.length_diff);
            ASSERT_EQ(stop_oracle, delta.stop_gained ? 1 : 0);
            if (delta.stop_gained) fs_stop_seen++;
            fs_seen++;
        } else if (net != 0) {
            if (delta.valid) ASSERT(!delta.frameshift);
            inframe_seen++;
        }
    }
    ASSERT(fs_seen > 0u);
    ASSERT(inframe_seen > 0u);
    fprintf(stderr,
            "[frameshift length-oracle coverage] frameshift=%u inframe_len=%u stop_gained=%u\n",
            fs_seen, inframe_seen, fs_stop_seen);
    PASS();
}

static int kprop_sequence_delta_equal(const duckvep_sequence_delta_t *a,
                                      const duckvep_sequence_delta_t *b) {
    return a != NULL && b != NULL &&
           a->cdna_pos == b->cdna_pos &&
           a->cds_pos == b->cds_pos &&
           a->protein_pos == b->protein_pos &&
           a->ref_aa == b->ref_aa &&
           a->alt_aa == b->alt_aa &&
           a->synonymous == b->synonymous &&
           a->missense == b->missense &&
           a->stop_gained == b->stop_gained &&
           a->stop_lost == b->stop_lost &&
           a->stop_retained == b->stop_retained &&
           a->start_lost == b->start_lost &&
           a->start_retained == b->start_retained &&
           a->frameshift == b->frameshift &&
           a->inframe_deletion == b->inframe_deletion &&
           a->inframe_insertion == b->inframe_insertion &&
           a->protein_altering == b->protein_altering &&
           a->coding_unknown == b->coding_unknown &&
           a->partial_codon == b->partial_codon &&
           a->valid == b->valid;
}

static int kprop_delta_is_coarse_cross_codon_missense(const duckvep_sequence_delta_t *d) {
    return d != NULL && d->valid && d->missense && !d->synonymous &&
           !d->stop_gained && !d->stop_lost && !d->stop_retained &&
           !d->start_lost && !d->start_retained && !d->frameshift &&
           !d->inframe_deletion &&
           !d->inframe_insertion && !d->protein_altering && !d->coding_unknown &&
           !d->partial_codon &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == -1 &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

static int kprop_delta_is_inframe_deletion_at(const duckvep_sequence_delta_t *d,
                                               int32_t protein_pos) {
    return d != NULL && d->valid && d->inframe_deletion &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->frameshift &&
           !d->inframe_insertion && !d->protein_altering && !d->coding_unknown &&
           !d->partial_codon &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

static int kprop_delta_is_inframe_insertion_at(const duckvep_sequence_delta_t *d,
                                                int32_t protein_pos) {
    return d != NULL && d->valid && d->inframe_insertion &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->frameshift &&
           !d->inframe_deletion && !d->protein_altering && !d->coding_unknown &&
           !d->partial_codon &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

static int kprop_delta_is_protein_altering_at(const duckvep_sequence_delta_t *d,
                                               int32_t protein_pos) {
    return d != NULL && d->valid && d->protein_altering &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->frameshift && !d->inframe_deletion && !d->inframe_insertion &&
           !d->coding_unknown && !d->partial_codon &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

TEST coding_context_delta_inframe_deletion_known_scene(void) {
    static const uint8_t short_cds[9] = {
        'A','T','G',  'A','A','A',  'T','T','T'
    };
    static const uint8_t mixed_short_cds[9] = {
        'a','U','g',  'a','a','a',  'u','U','U'
    };
    static const uint8_t repeat_cds[15] = {
        'A','T','G',  'A','A','A',  'A','A','A',  'C','C','C',  'T','T','T'
    };
    static const uint8_t ambiguous_cds[9] = {
        'A','T','G',  'N','N','N',  'T','T','T'
    };
    static const uint8_t short_alt_cds[6] = {
        'A','T','G',  'T','T','T'
    };
    static const uint8_t fake_ref_pep[4] = { 'M', 'K', 'F', '\0' };
    static const uint8_t fake_alt_pep[3] = { 'M', 'F', '\0' };
    static const uint8_t x_ref_pep[4] = { 'M', 'X', 'F', '\0' };
    static const uint8_t multi_cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','A',  'T','T','T'
    };
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_edit_t edits[2];
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[32];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 3u; edit.ref = short_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(short_cds, sizeof short_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 3u; edit.ref = mixed_short_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(mixed_short_cds, sizeof mixed_short_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 3u; edit.ref = repeat_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(repeat_cds, sizeof repeat_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT(ctx.has_single_edit);
    ASSERT_EQ(4u, ctx.single_edit_cds_start);
    ASSERT_EQ(3u, ctx.single_edit_ref_len);
    ASSERT_EQ(0u, ctx.single_edit_alt_len);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 6u; edit.ref = multi_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(multi_cds, sizeof multi_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    /* Non-codon-aligned pure in-frame deletion (cds_start 5, mid-codon): removes CDS 5-7 from
     * M K P G L F, merging codons 2-3 into M [T] G L F — a clean one-residue in-frame deletion
     * with a changed junction codon. VEP's codon-allele trim (AAACCC vs ACC -> alt empties,
     * ref len 3 % 3 == 0) calls this inframe_deletion regardless of alignment; the generalized
     * classifier now resolves it (the old codon-aligned-only slice bailed to unsupported).
     * protein_pos is the first affected codon, ((5-1)/3)+1 = 2. */
    edit.cds_start = 5u; edit.ref_len = 3u; edit.ref = multi_cds + 4u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(multi_cds, sizeof multi_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    /* Non-aligned deletion whose merged junction codon is a premature stop:
     * VEP evaluates inframe_deletion and stop_gained independently over the
     * same codon window, so both facts apply. */
    {
        static const uint8_t junction_stop_cds[18] = {
            'A','T','G',  'T','C','A',  'G','A','A',  'G','G','G',  'T','T','T',  'A','A','A'
        };
        edit.cds_start = 5u; edit.ref_len = 3u; edit.ref = junction_stop_cds + 4u;
        edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(junction_stop_cds, sizeof junction_stop_cds,
                                               &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_deletion && delta.stop_gained &&
               !delta.frameshift && !delta.protein_altering);
    }

    /* Removing the terminal stop codon: the ref local peptide is "*" and the
     * alternate is empty, so VEP resolves stop_lost rather than an in-frame
     * deletion. The shared terminal evaluator owns that distinction. */
    {
        static const uint8_t stop_del_cds[12] = {
            'A','T','G',  'A','A','A',  'C','C','C',  'T','A','A'
        };
        edit.cds_start = 10u; edit.ref_len = 3u; edit.ref = stop_del_cds + 9u;
        edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(stop_del_cds, sizeof stop_del_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_lost && !delta.inframe_deletion &&
               !delta.frameshift);
    }

    /* Deleting the complete start codon is both inframe_deletion and start_lost
     * in VEP 116 when evaluated in its transcript UTR+CDS sequence. */
    {
        static const uint8_t start_del_cds[12] = {
            'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G'
        };
        static const uint8_t pre_cds[3] = { 'C','C','C' };
        edit.cds_start = 1u; edit.ref_len = 3u; edit.ref = start_del_cds;
        edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(start_del_cds, sizeof start_del_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ctx.pre_cds_bases = pre_cds;
        ctx.pre_cds_length = sizeof pre_cds;
        ctx.pre_cds_complete = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_deletion && delta.start_lost &&
               !delta.start_retained);
    }

    edits[0].cds_start = 10u; edits[0].ref_len = 3u; edits[0].ref = multi_cds + 9u;
    edits[0].alt_len = 0u; edits[0].alt = NULL; edits[0].variant_strand = 1;
    edits[1].cds_start = 4u; edits[1].ref_len = 3u; edits[1].ref = multi_cds + 3u;
    edits[1].alt_len = 0u; edits[1].alt = NULL; edits[1].variant_strand = 1;
    edit_set.edits = edits; edit_set.count = 2u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(multi_cds, sizeof multi_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = multi_cds; ctx.ref_cds_len = sizeof multi_cds;
    ctx.alt_cds = multi_cds; ctx.alt_cds_len = sizeof multi_cds;
    ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 6u;
    ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 6u;
    ctx.length_diff = -3;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = ambiguous_cds; ctx.ref_cds_len = sizeof ambiguous_cds;
    ctx.alt_cds = short_alt_cds; ctx.alt_cds_len = sizeof short_alt_cds;
    ctx.ref_peptide = fake_ref_pep; ctx.ref_peptide_len = 3u;
    ctx.alt_peptide = fake_alt_pep; ctx.alt_peptide_len = 2u;
    ctx.length_diff = -3;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = short_cds; ctx.ref_cds_len = sizeof short_cds;
    ctx.alt_cds = short_alt_cds; ctx.alt_cds_len = sizeof short_alt_cds;
    ctx.ref_peptide = x_ref_pep; ctx.ref_peptide_len = 3u;
    ctx.alt_peptide = fake_alt_pep; ctx.alt_peptide_len = 2u;
    ctx.length_diff = -3;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = multi_cds; ctx.ref_cds_len = sizeof multi_cds;
    ctx.alt_cds = multi_cds; ctx.alt_cds_len = sizeof multi_cds;
    ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 6u;
    ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 6u;
    ctx.length_diff = INT64_MIN;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);
    PASS();
}

TEST coding_context_delta_inframe_insertion_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const uint8_t direct_cds[12] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'T','T','T'
    };
    static const uint8_t terminal_cds[12] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'T','A','A'
    };
    static const uint8_t internal_stop_cds[15] = {
        'A','T','G',  'T','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const uint8_t stop_window_cds[12] = {
        'A','T','G',  'A','A','A',  'T','A','C',  'G','G','G'
    };
    static const uint8_t insert_gcc[3] = { 'G','C','C' };
    static const uint8_t insert_gcc2[6] = { 'G','C','C', 'G','C','C' };
    static const uint8_t insert_stop[3] = { 'T','A','A' };
    static const uint8_t insert_stop_before_ref_flank[3] = { 'A','T','A' };
    static const uint8_t insert_ref_flank_before_stop[6] = {
        'A','C','T', 'A','A','G'
    };
    static const uint8_t insert_atg[3] = { 'A','T','G' };
    static const uint8_t insert_after_retained_start[3] = { 'G','T','T' };
    static const uint8_t insert_start_retained[9] = {
        'G','C','G', 'T','T','G', 'G','C','A'
    };
    static const uint8_t insert_start_retained_stop[9] = {
        'G','T','C', 'A','T','C', 'C','T','A'
    };
    static const uint8_t terminal_insert_agc[3] = { 'A','G','C' };
    static const uint8_t terminal_insert_aggt[4] = { 'A','G','G','T' };
    static const uint8_t terminal_insert_stop_lost[11] = {
        'C','G','A','T','G','T','T','A','T','G','A'
    };
    static const uint8_t terminal_before_taaa[4] = { 'T','A','A','A' };
    static const uint8_t terminal_before_c[1] = { 'C' };
    static const uint8_t terminal_before_stop_gained[6] = {
        'A','A','A','T','A','A'
    };
    static const uint8_t pre_cds[3] = { 'C','C','C' };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_haplotype_edit_t edits[4];
        uint8_t alt_cds[32];
        uint8_t ref_pep[16];
        uint8_t alt_pep[16];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;
        uint32_t anchor_cds;
        uint32_t i;

        memset(&s, 0, sizeof s);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        anchor_cds = s.strand > 0 ? 6u : 7u;
        s.vpos = kprop_genomic_pos_for_cds(&s, anchor_cds); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_INS;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        s.abytes[1] = s.abytes[0];
        for (i = 0u; i < 3u; i++) {
            char b = s.strand > 0 ? (char)insert_gcc[i]
                                  : kprop_complement_base((char)insert_gcc[2u - i]);
            s.abytes[2u + i] = (uint8_t)b;
        }
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edits, 4u, alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(ctx.has_single_edit);
        ASSERT_EQ(7u, ctx.single_edit_cds_start);
        ASSERT_EQ(0u, ctx.single_edit_ref_len);
        ASSERT_EQ(3u, ctx.single_edit_alt_len);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 3));
    }

    {
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        uint8_t alt_cds[32];
        uint8_t ref_pep[16];
        uint8_t alt_pep[16];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;

        memset(&edit, 0, sizeof edit);
        edit.cds_start = 7u; edit.ref_len = 0u; edit.ref = NULL;
        edit.alt_len = 6u; edit.alt = insert_gcc2; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 3));

        edit.alt_len = 3u; edit.alt = insert_stop;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_gained &&
               !delta.stop_lost && !delta.stop_retained);

        /* The two VEP predicates deliberately see different peptide values. Inserting
         * ATA before the last base of TAC makes the local peptide Y -> *Y. The in-frame
         * predicate truncates the alternate to '*', while protein_altering sees raw '*Y';
         * neither shape term applies, leaving stop_gained alone. */
        edit.cds_start = 9u;
        edit.alt_len = sizeof insert_stop_before_ref_flank;
        edit.alt = insert_stop_before_ref_flank;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      stop_window_cds, sizeof stop_window_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && !delta.inframe_insertion &&
               !delta.protein_altering && !delta.frameshift);

        /* Nearby positive control: inserting ACTAAG before the second base of TAC makes
         * Y -> Y*D. The preserved Y occurs before the new stop, so VEP emits both terms. */
        edit.cds_start = 8u;
        edit.alt_len = sizeof insert_ref_flank_before_stop;
        edit.alt = insert_ref_flank_before_stop;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      stop_window_cds, sizeof stop_window_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && delta.inframe_insertion &&
               !delta.protein_altering && !delta.frameshift);

        /* VEP first establishes stop_retained, then inframe_insertion trims the
         * alternate *P peptide to its first stop. The resulting * prefix still
         * matches the reference, so both predicates apply. */
        edit.cds_start = 6u;
        edit.alt_len = sizeof insert_gcc;
        edit.alt = insert_gcc;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      internal_stop_cds, sizeof internal_stop_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained && delta.inframe_insertion &&
               !delta.protein_altering && !delta.stop_gained && !delta.stop_lost);

        /* Moving the same ATG insertion across the terminal TAA distinguishes
         * stop_lost from inframe_insertion&stop_retained_variant. */
        edit.alt = insert_atg;
        edit.alt_len = sizeof insert_atg;
        edit.cds_start = 11u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_lost && !delta.stop_gained &&
               !delta.stop_retained && !delta.inframe_insertion);

        edit.cds_start = 12u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained && delta.inframe_insertion &&
               !delta.stop_gained && !delta.stop_lost);

        /* VEP's terminal-stop insertion predicates are deliberately not reducible
         * to length modulo three. These are executable witnesses from the pinned
         * VEP-116 state machine; see ERRATA.md. */
        edit.cds_start = 11u;
        edit.alt_len = sizeof terminal_insert_agc;
        edit.alt = terminal_insert_agc;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
               !delta.frameshift && !delta.stop_lost && !delta.coding_unknown);

        edit.alt_len = sizeof terminal_insert_aggt;
        edit.alt = terminal_insert_aggt;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.coding_unknown &&
               !delta.frameshift && !delta.stop_lost && !delta.stop_retained);

        edit.alt_len = sizeof terminal_insert_stop_lost;
        edit.alt = terminal_insert_stop_lost;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_lost &&
               !delta.frameshift && !delta.stop_retained && !delta.coding_unknown);

        edit.cds_start = 10u;
        edit.alt_len = sizeof terminal_before_taaa;
        edit.alt = terminal_before_taaa;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
               !delta.frameshift && !delta.stop_gained);

        edit.alt_len = sizeof terminal_before_c;
        edit.alt = terminal_before_c;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.inframe_insertion &&
               !delta.stop_retained && !delta.stop_gained);

        edit.alt_len = sizeof terminal_before_stop_gained;
        edit.alt = terminal_before_stop_gained;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_gained &&
               !delta.frameshift && !delta.stop_retained);

        /* The coding context is transcript-oriented; the same predicate state
         * must be independent of genomic strand. */
        edit.cds_start = 11u;
        edit.alt_len = sizeof terminal_insert_aggt;
        edit.alt = terminal_insert_aggt;
        edit.variant_strand = -1;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, -1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.coding_unknown &&
               !delta.frameshift && !delta.stop_lost && !delta.stop_retained);
        edit.variant_strand = 1;

        /* Mid-codon insertion that PRESERVES the flanking residue is inframe_insertion, not
         * protein_altering and not a bail: GCC after CDS 5 makes codon 2 AAA->AAG (still Lys)
         * and inserts Pro, so the ref window is empty. VEP calls this inframe_insertion; the old
         * codon-boundary-only classifier bailed. protein_pos = (5/3)+1 = 2. */
        edit.alt_len = sizeof insert_gcc;
        edit.alt = insert_gcc;
        edit.cds_start = 6u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 2));

        /* Insertion immediately after the complete start codon (before_cds == 3): the start Met
         * is preserved, Ala inserted after it -> inframe_insertion at protein_pos 2. The old
         * classifier required before_cds > 3 and bailed. */
        edit.cds_start = 4u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 2));

        /* Held-out seed 71, chrDuck:121 T>TGTT after anchor removal. Inserting
         * GTT before CDS base 3 keeps ATG as the alternate prefix and adds a
         * residue after it. start_retained_variant therefore does not suppress
         * inframe_insertion; VEP suppresses only when the reference peptide is
         * the alternate suffix (an insertion before the translated start). */
        edit.cds_start = 3u;
        edit.alt_len = sizeof insert_after_retained_start;
        edit.alt = insert_after_retained_start;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ctx.pre_cds_bases = pre_cds;
        ctx.pre_cds_length = sizeof pre_cds;
        ctx.pre_cds_complete = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_retained && delta.inframe_insertion &&
               !delta.start_lost && !delta.stop_gained);

        /* VEP 116 evaluates the start and insertion predicates independently.
         * Inserting nine bases before CDS base 3 can preserve ATG while adding
         * residues: both start_retained_variant and inframe_insertion apply. */
        edit.alt_len = sizeof insert_start_retained;
        edit.alt = insert_start_retained;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ctx.pre_cds_bases = pre_cds;
        ctx.pre_cds_length = sizeof pre_cds;
        ctx.pre_cds_complete = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_retained && delta.inframe_insertion &&
               !delta.start_lost && !delta.stop_gained);

        /* The same state can gain a stop: AT + GTCATCCTA + the original G
         * translates MSS*. VEP emits all three facts. */
        edit.alt = insert_start_retained_stop;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ctx.pre_cds_bases = pre_cds;
        ctx.pre_cds_length = sizeof pre_cds;
        ctx.pre_cds_complete = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_retained && delta.inframe_insertion &&
               delta.stop_gained && !delta.start_lost);

        /* Insertion between the penultimate and last codon (before_cds == 9, one short of the
         * old before_cds < cds_len-3 bound): flank-preserving -> inframe_insertion at pp 4. */
        edit.cds_start = 10u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 4));
    }
    PASS();
}

TEST coding_context_delta_delins_known_scene(void) {
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t alt3[3] = { 'G','C','C' };
    static const uint8_t alt6[6] = { 'G','C','C',  'G','C','T' };
    static const uint8_t stop_alt6[6] = { 'G','C','C',  'T','A','A' };
    static const uint8_t ambiguous_cds[18] = {
        'A','T','G',  'A','A','A',  'N','N','N',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t ambiguous_alt_cds[21] = {
        'A','T','G',  'A','A','A',  'G','C','C',  'G','C','T',
        'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t fake_ref_pep[6] = { 'M', 'K', 'X', 'G', 'F', 'E' };
    static const uint8_t fake_alt_pep[7] = { 'M', 'K', 'A', 'A', 'G', 'F', 'E' };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_haplotype_edit_t edits[4];
        uint8_t alt_cds[40];
        uint8_t ref_pep[20];
        uint8_t alt_pep[20];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;
        uint32_t cds_start;
        uint32_t i;

        memset(&s, 0, sizeof s);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
        s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 18u);

        cds_start = 7u;
        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 2u);
        s.vend = s.vpos + 2u; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < 3u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 3u, alt6, 6u);
        s.roff = 0u; s.aoff = 3u; s.rlen = 3u; s.alen = 6u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edits, 4u, alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(ctx.has_single_edit);
        ASSERT_EQ(7u, ctx.single_edit_cds_start);
        ASSERT_EQ(3u, ctx.single_edit_ref_len);
        ASSERT_EQ(6u, ctx.single_edit_alt_len);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 5u);
        s.vend = s.vpos + 5u;
        for (i = 0u; i < 6u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 6u, alt3, 3u);
        s.roff = 0u; s.aoff = 6u; s.rlen = 6u; s.alen = 3u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edits, 4u, alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(ctx.has_single_edit);
        ASSERT_EQ(7u, ctx.single_edit_cds_start);
        ASSERT_EQ(6u, ctx.single_edit_ref_len);
        ASSERT_EQ(3u, ctx.single_edit_alt_len);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));
    }

    {
        duckvep_haplotype_edit_t edit;
        duckvep_haplotype_edit_t edits[2];
        duckvep_edit_set_t edit_set;
        uint8_t alt_cds[40];
        uint8_t ref_pep[20];
        uint8_t alt_pep[20];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;

        memset(&edit, 0, sizeof edit);
        edit.ref_len = 3u; edit.ref = cds; edit.alt_len = 6u; edit.alt = alt6;
        edit.variant_strand = 1; edit_set.edits = &edit; edit_set.count = 1u;

        edit.cds_start = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_lost && !delta.inframe_insertion &&
               !delta.protein_altering);

        edit.cds_start = 16u; edit.ref = cds + 15u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 6));

        edit.cds_start = 8u; edit.ref = cds + 7u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        edit.cds_start = 7u; edit.ref = cds + 6u; edit.alt = stop_alt6;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.protein_altering && delta.stop_gained &&
               !delta.inframe_insertion && !delta.inframe_deletion);

        edits[0].cds_start = 13u; edits[0].ref_len = 3u; edits[0].ref = cds + 12u;
        edits[0].alt_len = 3u; edits[0].alt = alt3; edits[0].variant_strand = 1;
        edits[1].cds_start = 7u; edits[1].ref_len = 3u; edits[1].ref = cds + 6u;
        edits[1].alt_len = 6u; edits[1].alt = alt6; edits[1].variant_strand = 1;
        edit_set.edits = edits; edit_set.count = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);

        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ambiguous_cds; ctx.ref_cds_len = sizeof ambiguous_cds;
        ctx.alt_cds = ambiguous_alt_cds; ctx.alt_cds_len = sizeof ambiguous_alt_cds;
        ctx.ref_peptide = fake_ref_pep; ctx.ref_peptide_len = sizeof fake_ref_pep;
        ctx.alt_peptide = fake_alt_pep; ctx.alt_peptide_len = sizeof fake_alt_pep;
        ctx.length_diff = 3; ctx.cds_changed = 1u;
        ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
        ctx.single_edit_cds_start = 7u; ctx.single_edit_ref_len = 3u;
        ctx.single_edit_alt_len = 6u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
        s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 18u);
        s.vpos = kprop_genomic_pos_for_cds(&s, 7u); s.vend = s.vpos + 2u;
        s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < 3u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 3u, alt6, 6u);
        s.roff = 0u; s.aoff = 3u; s.rlen = 3u; s.alen = 6u;
        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        /* Replacing P with AA increases the CDS by one codon but preserves neither
         * peptide edge. VEP calls the local shape protein_altering_variant, not an
         * in-frame insertion inferred from net length alone. */
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING), rows[0].consequence_mask);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(0u, stats->substitution_context);
        ASSERT_EQ(0u, stats->del_context);
        ASSERT_EQ(0u, stats->ins_context);
        ASSERT_EQ(1u, stats->indel_context);
        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }
    PASS();
}

TEST sequence_delta_with_scratch_indel_known_scene(void) {
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t alt3[3] = { 'G','C','C' };
    static const uint8_t alt6[6] = { 'G','C','C',  'G','C','T' };
    static const uint8_t alt2[2] = { 'G','T' };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_haplotype_edit_t edits[4];
        uint8_t alt_cds[40];
        uint8_t ref_pep[20];
        uint8_t alt_pep[20];
        duckvep_delta_scratch_t scratch;
        duckvep_sequence_delta_t delta;
        uint32_t cds_start = 7u;
        uint32_t i;

        memset(&s, 0, sizeof s);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
        s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 18u);
        memset(&scratch, 0, sizeof scratch);
        scratch.edits = edits; scratch.edits_cap = 4u;
        scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
        scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
        scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 2u);
        s.vend = s.vpos + 2u; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < 3u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 3u, alt6, 6u);
        s.roff = 0u; s.aoff = 3u; s.rlen = 3u; s.alen = 6u;
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, &scratch, &delta);
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 5u);
        s.vend = s.vpos + 5u;
        for (i = 0u; i < 6u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 6u, alt3, 3u);
        s.roff = 0u; s.aoff = 6u; s.rlen = 6u; s.alen = 3u;
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, &scratch, &delta);
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start);
        s.vend = s.vpos; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        kprop_fill_variant_alt_from_tx(&s, 1u, alt2, 2u);
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 2u;
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, &scratch, &delta);
        /* +1 delins at codon 3 (body) with the ATG start intact: the general
         * CodingContext now resolves the frameshift the direct body-only path skipped. */
        ASSERT(kprop_delta_is_frameshift_at(&delta, 3));
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, NULL, &delta);
        ASSERT(delta.valid);
        ASSERT(delta.frameshift);
        ASSERT(!delta.inframe_insertion);
        ASSERT(!delta.inframe_deletion);
    }
    PASS();
}

TEST annotate_delins_boundary_no_route_known_scene(void) {
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t alt3[3] = { 'G','C','C' };
    static const uint8_t alt6[6] = { 'G','C','C',  'G','C','T' };
    static const uint32_t starts[4] = { 1u, 16u, 13u, 8u };
    static const uint32_t ref_lens[4] = { 3u, 3u, 6u, 3u };
    static const uint32_t alt_lens[4] = { 6u, 6u, 3u, 6u };
    static const uint64_t expected_masks[4] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING)
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;
    size_t strand_idx;

    for (strand_idx = 0u; strand_idx < 2u; strand_idx++) {
        for (case_idx = 0u; case_idx < 4u; case_idx++) {
            struct kprop_coding s;
            duckvep_model_t *model = NULL;
            duckvep_options_t *opts = NULL;
            duckvep_workspace_t *ws = NULL;
            const duckvep_workspace_delta_route_stats_t *stats;
            duckvep_error_t err;
            duckvep_consequence_t rows[2];
            duckvep_result_builder_t rb;
            const uint8_t *alt_tx = alt_lens[case_idx] == 3u ? alt3 : alt6;
            uint32_t cds_start = starts[case_idx];
            uint32_t ref_len = ref_lens[case_idx];
            uint32_t alt_len = alt_lens[case_idx];
            uint32_t i;

            memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
            s.cds = cds; s.chrom = 0u; s.strand = strands[strand_idx]; s.flags = 0u;
            s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
            s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
            s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
            kprop_wire_coding_scene(&s, 18u);
            s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                                  : kprop_genomic_pos_for_cds(&s, cds_start + ref_len - 1u);
            s.vend = s.vpos + ref_len - 1u; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
            for (i = 0u; i < ref_len; i++) {
                s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
            }
            kprop_fill_variant_alt_from_tx(&s, ref_len, alt_tx, alt_len);
            s.roff = 0u; s.aoff = ref_len;
            s.rlen = (uint16_t)ref_len; s.alen = (uint16_t)alt_len;

            ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
            ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
            ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
            duckvep_workspace_delta_route_stats_reset(ws);
            duckvep_result_builder_init(&rb, rows, 2u);
            ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
            ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
            ASSERT_EQ(expected_masks[case_idx], rows[0].consequence_mask);
            stats = duckvep_workspace_delta_route_stats(ws);
            ASSERT(stats != NULL);
            ASSERT_EQ(0u, stats->substitution_context);
            ASSERT_EQ(0u, stats->del_context);
            ASSERT_EQ(0u, stats->ins_context);

            duckvep_workspace_close(ws);
            duckvep_options_close(opts);
            duckvep_model_close(model);
        }
    }
    PASS();
}

TEST annotate_inframe_insertion_route_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint8_t alt_tx[3] = { 'G','C','C' };
        uint32_t anchor_cds;
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        anchor_cds = s.strand > 0 ? 6u : 7u;
        s.vpos = kprop_genomic_pos_for_cds(&s, anchor_cds); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_INS;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        s.abytes[1] = s.abytes[0];
        for (i = 0u; i < 3u; i++) {
            char b = s.strand > 0 ? (char)alt_tx[i]
                                  : kprop_complement_base((char)alt_tx[2u - i]);
            s.abytes[2u + i] = (uint8_t)b;
        }
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;

        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION), rows[0].consequence_mask);
        ASSERT_EQ(-1, rows[0].cdna_pos);
        ASSERT_EQ(-1, rows[0].cds_pos);
        ASSERT_EQ(3, rows[0].protein_pos);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(1u, stats->simple_indel);
        ASSERT_EQ(0u, stats->substitution_context);
        ASSERT_EQ(0u, stats->del_context);
        ASSERT_EQ(0u, stats->ins_context);

        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }

    {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint8_t alt_tx[3] = { 'T','A','A' };
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        s.vpos = kprop_genomic_pos_for_cds(&s, 6u); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_INS;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        s.abytes[1] = s.abytes[0];
        for (i = 0u; i < 3u; i++) s.abytes[2u + i] = alt_tx[i];
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;

        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) |
                      DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
                  rows[0].consequence_mask);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(0u, stats->simple_indel);
        ASSERT_EQ(1u, stats->ins_context);

        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }
    PASS();
}

TEST annotate_inframe_deletion_route_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint32_t anchor_cds;
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        anchor_cds = s.strand > 0 ? 3u : 7u;
        s.vpos = kprop_genomic_pos_for_cds(&s, anchor_cds); s.vend = s.vpos + 3u;
        s.vkind = (uint8_t)DUCKVEP_KIND_DEL;
        for (i = 0u; i < 4u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        s.abytes[4] = s.abytes[0];
        s.roff = 0u; s.aoff = 4u; s.rlen = 4u; s.alen = 1u;

        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION), rows[0].consequence_mask);
        ASSERT_EQ(-1, rows[0].cdna_pos);
        ASSERT_EQ(-1, rows[0].cds_pos);
        ASSERT_EQ(2, rows[0].protein_pos);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(1u, stats->simple_indel);
        ASSERT_EQ(0u, stats->substitution_context);
        ASSERT_EQ(0u, stats->del_context);

        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }
    PASS();
}

TEST workspace_delta_scratch_caps_known(void) {
    static const uint16_t chrom[3] = {0u, 0u, 0u};
    static const uint32_t start1[3] = {100u, 200u, 300u};
    static const uint32_t end1[3] = {100u, 208u, 313u};
    static const int8_t strand[3] = {1, 1, -1};
    static const uint64_t flags[3] = {0u, 0u, 0u};
    static const uint32_t exon_off[3] = {0u, 0u, 1u};
    static const uint16_t exon_cnt[3] = {0u, 1u, 1u};
    static const uint32_t cds_s[3] = {0u, 200u, 300u};
    static const uint32_t cds_e[3] = {0u, 208u, 313u};
    static const uint64_t cds_off[3] = {0u, 0u, 9u};
    static const uint32_t cds_len[3] = {0u, 9u, 14u};
    static const uint8_t table[3] = {0u, 1u, 1u};
    static const uint32_t exon_start[2] = {200u, 300u};
    static const uint32_t exon_end[2] = {208u, 313u};
    static const uint32_t exon_cdna_start[2] = {1u, 1u};
    static const uint32_t exon_cdna_end[2] = {9u, 14u};
    static const int8_t exon_phase[2] = {0, 0};
    static const uint8_t cds_bytes[23] = {
        'A','T','G','G','A','A','T','A','A',
        'A','T','G','C','C','C','G','G','G','T','T','T','A','A'
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_sequence_pool_t seq;
    duckvep_model_t *model = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_delta_scratch_t *scratch;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex); memset(&seq, 0, sizeof seq);
    memset(&err, 0, sizeof err);
    tx.chrom_id = chrom; tx.start1 = start1; tx.end1 = end1; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exon_off; tx.exon_count = exon_cnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 3u;
    ex.start1 = exon_start; ex.end1 = exon_end;
    ex.cdna_start1 = exon_cdna_start; ex.cdna_end1 = exon_cdna_end;
    ex.phase = exon_phase; ex.end_phase = exon_phase; ex.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = table;
    seq.transcript_count = 3u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(scratch != NULL);
    ASSERT_EQ((size_t)14u + (size_t)UINT16_MAX, scratch->alt_cds_cap);
    ASSERT_EQ(5u, scratch->ref_peptide_cap);
    ASSERT_EQ(((size_t)14u + (size_t)UINT16_MAX) / 3u + 1u, scratch->alt_peptide_cap);
    ASSERT_EQ(7u, scratch->edits_cap);
    ASSERT(scratch->alt_cds != NULL);
    ASSERT(scratch->ref_peptide != NULL);
    ASSERT(scratch->alt_peptide != NULL);
    ASSERT(scratch->edits != NULL);
    duckvep_workspace_close(ws); ws = NULL;
    duckvep_model_close(model); model = NULL;

    memset(&seq, 0, sizeof seq); memset(&err, 0, sizeof err);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(scratch != NULL);
    ASSERT_EQ(0u, scratch->alt_cds_cap);
    ASSERT_EQ(0u, scratch->ref_peptide_cap);
    ASSERT_EQ(0u, scratch->alt_peptide_cap);
    ASSERT_EQ(0u, scratch->edits_cap);
    ASSERT(scratch->alt_cds == NULL);
    ASSERT(scratch->ref_peptide == NULL);
    ASSERT(scratch->alt_peptide == NULL);
    ASSERT(scratch->edits == NULL);
    duckvep_workspace_close(ws); ws = NULL;
    duckvep_model_close(model); model = NULL;

    memset(&seq, 0, sizeof seq); memset(&err, 0, sizeof err);
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = table;
    seq.cds_bytes = NULL; seq.cds_bytes_len = 0u; seq.transcript_count = 3u;
    {
        static const uint64_t zero_off[3] = {0u, 0u, 0u};
        static const uint32_t zero_len[3] = {0u, 0u, 0u};
        seq.cds_offset = zero_off;
        seq.cds_length = zero_len;
        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, &seq, NULL, &model, &err));
    }
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(scratch != NULL);
    ASSERT_EQ(0u, scratch->alt_cds_cap);
    ASSERT_EQ(0u, scratch->ref_peptide_cap);
    ASSERT_EQ(0u, scratch->alt_peptide_cap);
    ASSERT_EQ(0u, scratch->edits_cap);
    ASSERT(scratch->alt_cds == NULL);
    ASSERT(scratch->ref_peptide == NULL);
    ASSERT(scratch->alt_peptide == NULL);
    ASSERT(scratch->edits == NULL);
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    PASS();
}

TEST workspace_delta_scratch_builds_lengthening_context(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_model_t *model = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_delta_scratch_t *ws_scratch;
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    duckvep_error_t err;
    uint8_t alt_tx[3] = {'G','C','C'};
    uint32_t i;

    memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 6u); s.vend = s.vpos;
    s.vkind = (uint8_t)DUCKVEP_KIND_INS;
    s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
    s.abytes[1] = s.abytes[0];
    for (i = 0u; i < 3u; i++) s.abytes[2u + i] = alt_tx[i];
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ws_scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(ws_scratch != NULL);
    ASSERT(ws_scratch->alt_cds_cap >= (size_t)s.cds_lenv + (size_t)UINT16_MAX);
    ASSERT(ws_scratch->alt_peptide_cap >= ws_scratch->alt_cds_cap / 3u + 1u);
    ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
              duckvep_variant_physical_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                   0u, 0u, s.strand,
                                                   ws_scratch->edits,
                                                   ws_scratch->edits_cap,
                                                   ws_scratch->alt_cds,
                                                   ws_scratch->alt_cds_cap,
                                                   ws_scratch->ref_peptide,
                                                   ws_scratch->ref_peptide_cap,
                                                   ws_scratch->alt_peptide,
                                                   ws_scratch->alt_peptide_cap,
                                                   &ctx));
    ASSERT_EQ(3, ctx.length_diff);
    ASSERT_EQ((size_t)18u, ctx.alt_cds_len);
    ASSERT_EQ((size_t)5u, ctx.ref_peptide_len);
    ASSERT_EQ((size_t)6u, ctx.alt_peptide_len);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 3));
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    PASS();
}

TEST workspace_delta_scratch_usable_for_mnv(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_model_t *model = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_delta_scratch_t *ws_scratch;
    duckvep_haplotype_edit_t stack_edits[4];
    uint8_t stack_alt_cds[32];
    uint8_t stack_ref_pep[16];
    uint8_t stack_alt_pep[16];
    duckvep_delta_scratch_t stack_scratch;
    duckvep_sequence_delta_t got;
    duckvep_sequence_delta_t want;
    duckvep_error_t err;
    uint8_t alt_tx[2] = {'C','C'};
    uint32_t i;

    memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 4u); s.vend = s.vpos + 1u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 2u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 2u, alt_tx, 2u);
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ws_scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(ws_scratch != NULL);
    ASSERT(ws_scratch->edits_cap >= 1u);
    ASSERT(ws_scratch->alt_cds_cap >= s.cds_lenv);
    ASSERT(ws_scratch->ref_peptide_cap >= s.cds_lenv / 3u + 1u);
    ASSERT(ws_scratch->alt_peptide_cap >= s.cds_lenv / 3u + 1u);

    memset(&stack_scratch, 0, sizeof stack_scratch);
    stack_scratch.edits = stack_edits; stack_scratch.edits_cap = 4u;
    stack_scratch.alt_cds = stack_alt_cds; stack_scratch.alt_cds_cap = sizeof stack_alt_cds;
    stack_scratch.ref_peptide = stack_ref_pep; stack_scratch.ref_peptide_cap = sizeof stack_ref_pep;
    stack_scratch.alt_peptide = stack_alt_pep; stack_scratch.alt_peptide_cap = sizeof stack_alt_pep;

    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &stack_scratch, &want);
    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             ws_scratch, &got);
    ASSERT(want.valid);
    ASSERT(got.valid);
    ASSERT_EQ(want.synonymous, got.synonymous);
    ASSERT_EQ(want.missense, got.missense);
    ASSERT_EQ(want.stop_gained, got.stop_gained);
    ASSERT_EQ(want.stop_lost, got.stop_lost);
    ASSERT_EQ(want.stop_retained, got.stop_retained);
    ASSERT_EQ(want.start_lost, got.start_lost);
    ASSERT_EQ(want.start_retained, got.start_retained);
    ASSERT_EQ(want.frameshift, got.frameshift);
    ASSERT_EQ(want.inframe_deletion, got.inframe_deletion);
    ASSERT_EQ(want.inframe_insertion, got.inframe_insertion);
    ASSERT_EQ(want.protein_altering, got.protein_altering);
    ASSERT_EQ(want.coding_unknown, got.coding_unknown);
    ASSERT_EQ(want.partial_codon, got.partial_codon);
    ASSERT_EQ(want.protein_pos, got.protein_pos);
    ASSERT_EQ(want.cdna_pos, got.cdna_pos);
    ASSERT_EQ(want.cds_pos, got.cds_pos);
    ASSERT_EQ(want.ref_aa, got.ref_aa);
    ASSERT_EQ(want.alt_aa, got.alt_aa);
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    PASS();
}

TEST sequence_delta_scratch_rejects_unequal_mnv_kind(void) {
    static uint8_t cds[9] = {
        'A','T','G',  'A','A','A',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t stack_edits[4];
    uint8_t stack_alt_cds[32];
    uint8_t stack_ref_pep[16];
    uint8_t stack_alt_pep[16];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t delta;
    uint8_t alt_tx[1] = { 'G' };
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1008u; s.cds_s = 1000u; s.cds_e = 1008u;
    s.es = 1000u; s.ee = 1008u; s.ecds = 1u; s.ecde = 9u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 9u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 3u); s.vend = s.vpos + 3u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 4u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 4u, alt_tx, 1u);
    s.roff = 0u; s.aoff = 4u; s.rlen = 4u; s.alen = 1u;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = stack_edits; scratch.edits_cap = 4u;
    scratch.alt_cds = stack_alt_cds; scratch.alt_cds_cap = sizeof stack_alt_cds;
    scratch.ref_peptide = stack_ref_pep; scratch.ref_peptide_cap = sizeof stack_ref_pep;
    scratch.alt_peptide = stack_alt_pep; scratch.alt_peptide_cap = sizeof stack_alt_pep;

    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &scratch, &delta);
    ASSERT(!delta.valid);
    PASS();
}

TEST sequence_delta_annotation_wrapper_del_insufficient_scratch_known(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t stack_edits[1];
    uint8_t stack_alt_cds[32];
    uint8_t stack_ref_pep[16];
    uint8_t stack_alt_pep[16];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 3u); s.vend = s.vpos + 3u;
    s.vkind = (uint8_t)DUCKVEP_KIND_DEL;
    for (i = 0u; i < 4u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    s.abytes[4] = s.abytes[0];
    s.roff = 0u; s.aoff = 4u; s.rlen = 4u; s.alen = 1u;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = stack_edits; scratch.edits_cap = 0u;
    scratch.alt_cds = stack_alt_cds; scratch.alt_cds_cap = sizeof stack_alt_cds;
    scratch.ref_peptide = stack_ref_pep; scratch.ref_peptide_cap = sizeof stack_ref_pep;
    scratch.alt_peptide = stack_alt_pep; scratch.alt_peptide_cap = sizeof stack_alt_pep;

    duckvep_sequence_delta_fill(DUCKVEP_KIND_DEL, &s.tx, &s.ex, &s.seq,
                                &s.v, 0u, 0u, s.vpos, s.strand, &shape);
    duckvep_sequence_delta_fill_for_annotation_trace(DUCKVEP_KIND_DEL, &s.tx, &s.ex,
                                                     &s.seq, &s.v, 0u, 0u, s.vpos,
                                                     s.strand, &scratch, NULL,
                                                     UINT32_MAX, UINT32_MAX,
                                                     &route, &routed);
    ASSERT(kprop_delta_is_inframe_deletion_at(&shape, 2));
    ASSERT(!routed.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_INTERNAL_CAPACITY,
              routed.sequence_status);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_DEL_CONTEXT, route);
    PASS();
}

TEST sequence_delta_annotation_wrapper_start_lost_mnv(void) {
    static uint8_t cds[9] = {'A','T','G', 'G','A','A', 'T','A','A'};
    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[32];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t routed;
    uint8_t alt_tx[2] = {'C','C'};
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1008u; s.cds_s = 1000u; s.cds_e = 1008u;
    s.es = 1000u; s.ee = 1008u; s.ecds = 1u; s.ecde = 9u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 9u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 1u); s.vend = s.vpos + 1u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 2u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 2u, alt_tx, 2u);
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 2u;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq, &s.v,
                                0u, 0u, s.vpos, s.strand, &shape);
    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &scratch, &direct);
    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV, &s.tx, &s.ex,
                                               &s.seq, &s.v, 0u, 0u, s.vpos,
                                               s.strand, &scratch, NULL, &routed);
    ASSERT(shape.valid);
    ASSERT(direct.valid);
    ASSERT(routed.valid);
    ASSERT(kprop_sequence_delta_equal(&shape, &direct));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT(routed.start_lost);
    ASSERT(!routed.missense);
    ASSERT_EQ((uint8_t)'M', routed.ref_aa);
    ASSERT_EQ(1, routed.protein_pos);

    s.flags = (uint64_t)DUCKVEP_TX_CDS_START_NF;
    kprop_wire_coding_scene(&s, 9u);
    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV, &s.tx, &s.ex,
                                               &s.seq, &s.v, 0u, 0u, s.vpos,
                                               s.strand, &scratch, NULL, &routed);
    ASSERT(routed.valid);
    ASSERT(!routed.start_lost);
    ASSERT(routed.missense);
    PASS();
}

static int kprop_cross_codon_scene_deltas(
    int8_t                         strand,
    uint8_t                       *cds,
    uint32_t                       cds_len,
    uint32_t                       first_cds,
    uint32_t                       ref_len,
    const uint8_t                 *alt_tx,
    uint32_t                       alt_len,
    duckvep_sequence_delta_t      *direct,
    duckvep_sequence_delta_t      *shape,
    duckvep_sequence_delta_t      *routed,
    duckvep_sequence_delta_route_t *route) {

    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    uint32_t i;

    if (cds == NULL || alt_tx == NULL || direct == NULL || shape == NULL || routed == NULL ||
        route == NULL || ref_len == 0u || ref_len + alt_len > sizeof s.abytes) {
        return 0;
    }
    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = strand; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1000u + cds_len - 1u;
    s.cds_s = s.tstart; s.cds_e = s.tend;
    s.es = s.tstart; s.ee = s.tend; s.ecds = 1u; s.ecde = cds_len;
    s.eph = 0; s.eeph = 0; s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, cds_len);
    s.vpos = strand > 0 ? kprop_genomic_pos_for_cds(&s, first_cds)
                         : kprop_genomic_pos_for_cds(&s, first_cds + ref_len - 1u);
    s.vend = s.vpos + ref_len - 1u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < ref_len; i++) {
        s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    }
    kprop_fill_variant_alt_from_tx(&s, ref_len, alt_tx, alt_len);
    s.roff = 0u; s.aoff = ref_len; s.rlen = (uint16_t)ref_len; s.alen = (uint16_t)alt_len;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &scratch, direct);
    duckvep_sequence_delta_fill(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq, &s.v,
                                0u, 0u, s.vpos, s.strand, shape);
    duckvep_sequence_delta_fill_for_annotation_trace(DUCKVEP_KIND_MNV, &s.tx, &s.ex,
                                                     &s.seq, &s.v, 0u, 0u, s.vpos,
                                                     s.strand, &scratch, NULL,
                                                     UINT32_MAX, UINT32_MAX,
                                                     route, routed);
    return 1;
}

TEST sequence_delta_with_scratch_cross_codon_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint8_t alt_tx[2] = {'C','C'};

    ASSERT(kprop_cross_codon_scene_deltas(1, cds, 15u, 6u, 2u, alt_tx, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&direct));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&shape));
    ASSERT(kprop_sequence_delta_equal(&shape, &direct));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    PASS();
}

TEST sequence_delta_with_scratch_cross_codon_reverse_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint8_t alt_tx[2] = {'C','C'};

    ASSERT(kprop_cross_codon_scene_deltas(-1, cds, 15u, 6u, 2u, alt_tx, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&direct));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&shape));
    ASSERT(kprop_sequence_delta_equal(&shape, &direct));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    PASS();
}

TEST sequence_delta_with_scratch_cross_codon_negative_scenes(void) {
    static uint8_t syn_cds[15] = {
        'A','T','G',  'A','A','A',  'T','T','A',  'G','A','A',  'T','T','T'
    };
    static uint8_t stop_cds[15] = {
        'A','T','G',  'T','C','A',  'A','A','A',  'G','A','A',  'T','T','T'
    };
    static uint8_t terminal_cds[9] = {
        'A','T','G',  'G','A','A',  'G','A','A'
    };
    static uint8_t wide_cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    static uint8_t length_cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint8_t syn_alt[2] = {'G','C'};
    uint8_t stop_alt[3] = {'A','G','C'};
    uint8_t terminal_alt[2] = {'C','C'};
    uint8_t wide_alt[4] = {'C','C','C','C'};
    uint8_t length_alt[1] = {'C'};

    /* These two-codon windows were unsupported by the old two-codon-missense-only slice and
     * fell back to coding_sequence_variant. The generalized window classifier now resolves
     * each one authoritatively (synonymous / stop_gained / missense, incl. the terminal-codon
     * and >3-base-wide windows). The narrow direct reference still only knows the missense
     * cross-codon window, so it stays invalid on the synonymous and stop-gained cases — which
     * is exactly the capability the interpreter adds. The router echoes the interpreter. */

    /* both codons synonymous -> synonymous_variant */
    ASSERT(kprop_cross_codon_scene_deltas(1, syn_cds, 15u, 6u, 2u, syn_alt, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.synonymous && direct.protein_pos == -1 &&
           !direct.missense && !direct.stop_gained && !direct.stop_lost &&
           !direct.stop_retained && !direct.start_lost && !direct.start_retained);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* junction codon becomes a stop -> stop_gained */
    ASSERT(kprop_cross_codon_scene_deltas(1, stop_cds, 15u, 5u, 3u, stop_alt, 3u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.stop_gained && direct.protein_pos == -1 &&
           !direct.missense && !direct.synonymous && !direct.stop_lost &&
           !direct.stop_retained && !direct.start_lost && !direct.start_retained);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* window includes the terminal (last) codon -> missense */
    ASSERT(kprop_cross_codon_scene_deltas(1, terminal_cds, 9u, 6u, 2u, terminal_alt, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.missense && direct.protein_pos == -1);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* four-base window spanning two codons -> missense */
    ASSERT(kprop_cross_codon_scene_deltas(1, wide_cds, 15u, 5u, 4u, wide_alt, 4u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.missense && direct.protein_pos == -1);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* unequal ref/alt length on a KIND_MNV is a disguised delins: genuinely unsupported. */
    ASSERT(kprop_cross_codon_scene_deltas(1, length_cds, 15u, 6u, 2u, length_alt, 1u,
                                          &direct, &shape, &routed, &route));
    ASSERT(!direct.valid);
    ASSERT(!shape.valid);
    ASSERT(!routed.valid);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    PASS();
}

TEST feature_substitution_window_fails_closed_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'N','A','A',  'C','C','C',  'T','G','G',  'T','A','A'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_event_t event;
    duckvep_sequence_delta_t delta;
    duckvep_sequence_delta_route_t route;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u;
    s.eph = 0; s.eeph = 0; s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 2u;
    kprop_wire_coding_scene(&s, sizeof cds);

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    /* The semantic edit is G>A at CDS position 3. The complete uploaded feature
     * also retains an N in the next codon, so VEP cannot use the widened peptide
     * window. Do not retry the apparently valid one-base edit. */
    s.vpos = 1002u; s.vend = 1003u;
    s.abytes[0] = 'G'; s.abytes[1] = 'N';
    s.abytes[2] = 'A'; s.abytes[3] = 'N';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event,
        UINT32_MAX, UINT32_MAX, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS, delta.sequence_status);

    /* The retained A is not the transcript's G at CDS position 12. The changed
     * T>C base itself is valid, so validating only the trimmed edit would miss
     * this reference mismatch and emit a supported terminal-stop consequence. */
    s.vpos = 1011u; s.vend = 1012u;
    s.abytes[0] = 'A'; s.abytes[1] = 'T';
    s.abytes[2] = 'A'; s.abytes[3] = 'C';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event,
        UINT32_MAX, UINT32_MAX, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH,
              delta.sequence_status);

    /* On the reverse strand genomic 1002-1003 maps to CDS positions 13-12.
     * The same checks must follow genomic byte order while comparing against
     * the reverse-complemented transcript bases. */
    s.strand = -1;
    kprop_wire_coding_scene(&s, sizeof cds);
    s.vpos = 1002u; s.vend = 1003u;
    s.abytes[0] = 'G'; s.abytes[1] = 'C';
    s.abytes[2] = 'G'; s.abytes[3] = 'T';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event,
        UINT32_MAX, UINT32_MAX, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH,
              delta.sequence_status);

    s.abytes[0] = 'N'; s.abytes[1] = 'C';
    s.abytes[2] = 'N'; s.abytes[3] = 'T';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event,
        UINT32_MAX, UINT32_MAX, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS, delta.sequence_status);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
} g_delta_wrapper_cov;

static enum theft_trial_res prop_delta_annotation_wrapper_matches_direct_shape(struct theft *t,
                                                                         void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    (void)t;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill((duckvep_variant_kind_t)s->vkind, &s->tx, &s->ex,
                                &s->seq, &s->v, 0u, 0u, s->vpos, s->strand,
                                &shape);
    duckvep_sequence_delta_fill_for_annotation((duckvep_variant_kind_t)s->vkind,
                                               &s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->vpos, s->strand,
                                               &scratch, NULL, &routed);
    if (!shape.valid || !routed.valid) return THEFT_TRIAL_FAIL;
    if (!kprop_sequence_delta_equal(&shape, &routed)) return THEFT_TRIAL_FAIL;

    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) g_delta_wrapper_cov.syn++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) g_delta_wrapper_cov.mis++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) g_delta_wrapper_cov.stop_gained++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) g_delta_wrapper_cov.stop_lost++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) g_delta_wrapper_cov.stop_retained++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_delta_wrapper_cov.fwd++;
    else g_delta_wrapper_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_annotation_wrapper_matches_direct_shape(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta annotation wrapper MNV == direct shape";
    cfg.prop1 = prop_delta_annotation_wrapper_matches_direct_shape;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_wrapper_cov, 0, sizeof g_delta_wrapper_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_wrapper_cov.syn > 0u);
    ASSERT(g_delta_wrapper_cov.mis > 0u);
    ASSERT(g_delta_wrapper_cov.stop_gained > 0u);
    ASSERT(g_delta_wrapper_cov.stop_lost > 0u);
    ASSERT(g_delta_wrapper_cov.stop_retained > 0u);
    ASSERT(g_delta_wrapper_cov.fwd > 0u);
    ASSERT(g_delta_wrapper_cov.rev > 0u);
    fprintf(stderr,
            "[delta-wrapper coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u\n",
            g_delta_wrapper_cov.syn, g_delta_wrapper_cov.mis,
            g_delta_wrapper_cov.stop_gained, g_delta_wrapper_cov.stop_lost,
            g_delta_wrapper_cov.stop_retained, g_delta_wrapper_cov.fwd,
            g_delta_wrapper_cov.rev);
    PASS();
}

static struct {
    uint32_t shape[5];
    uint32_t fwd;
    uint32_t rev;
} g_delta_exon_hint_cov;

static enum theft_trial_res prop_delta_exon_hint_matches_unhinted(struct theft *t,
                                                                  void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[96];
    uint8_t ref_pep[40];
    uint8_t alt_pep[40];
    duckvep_delta_scratch_t scratch;
    duckvep_event_t event;
    duckvep_sequence_delta_t unhinted;
    duckvep_sequence_delta_t hinted;
    duckvep_sequence_delta_route_t unhinted_route;
    duckvep_sequence_delta_route_t hinted_route;
    (void)t;

    if (s->expect_shape > KPROP_CDS_EDIT_INDEL) return THEFT_TRIAL_FAIL;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;
    duckvep_event_load(&s->v, 0u, &event);

    duckvep_sequence_delta_fill_for_annotation_trace(
        (duckvep_variant_kind_t)s->vkind, &s->tx, &s->ex, &s->seq, &s->v,
        0u, 0u, s->vpos, s->strand, &scratch, &event,
        (uint32_t)DUCKVEP_REGION_CDS, UINT32_MAX,
        &unhinted_route, &unhinted);
    duckvep_sequence_delta_fill_for_annotation_trace(
        (duckvep_variant_kind_t)s->vkind, &s->tx, &s->ex, &s->seq, &s->v,
        0u, 0u, s->vpos, s->strand, &scratch, &event,
        (uint32_t)DUCKVEP_REGION_CDS, 0u,
        &hinted_route, &hinted);

    if (unhinted_route != hinted_route ||
        !kprop_sequence_delta_equal(&unhinted, &hinted)) {
        return THEFT_TRIAL_FAIL;
    }
    duckvep_haplotype_edit_t projected;
    if (duckvep_variant_cds_edit_build(&s->tx, &s->ex, &s->seq, &s->v,
            0u, 0u, s->strand, &projected) != DUCKVEP_CDS_EDIT_OK)
        return THEFT_TRIAL_FAIL;
    duckvep_haplotype_edit_t unchanged = projected;
    struct kprop_context_snapshot computed, retained;
    duckvep_sequence_delta_t observed[2];
    duckvep_sequence_delta_route_t routes[2];
    duckvep_variant_coding_context_status_t statuses[2];
    for (unsigned cached = 0u; cached < 2u; cached++) {
        duckvep_coding_context_t context;
        duckvep_sequence_delta_fill_for_annotation_observed(
            (duckvep_variant_kind_t)s->vkind, &s->tx, &s->ex, &s->seq, &s->v,
            0u, 0u, s->vpos, s->strand, &scratch, &event,
            (uint32_t)DUCKVEP_REGION_CDS, UINT32_MAX, cached ? &projected : NULL,
            &routes[cached], &observed[cached], &context, &statuses[cached]);
        if (statuses[cached] == DUCKVEP_VARIANT_CODING_CONTEXT_OK &&
            !kprop_context_snapshot(&context, cached ? &retained : &computed))
            return THEFT_TRIAL_FAIL;
    }
    if (routes[0] != routes[1] || statuses[0] != statuses[1] ||
        memcmp(&observed[0], &observed[1], sizeof observed[0]) ||
        memcmp(&unchanged, &projected, sizeof projected) ||
        (statuses[0] == DUCKVEP_VARIANT_CODING_CONTEXT_OK &&
         memcmp(&computed, &retained, sizeof computed))) return THEFT_TRIAL_FAIL;
    g_delta_exon_hint_cov.shape[s->expect_shape]++;
    if (s->strand > 0) g_delta_exon_hint_cov.fwd++;
    else g_delta_exon_hint_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_annotation_exon_hint_matches_unhinted(void) {
    struct theft_run_config cfg;
    uint32_t shape;

    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta exon hint == unhinted projection";
    cfg.prop1 = prop_delta_exon_hint_matches_unhinted;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_exon_hint_cov, 0, sizeof g_delta_exon_hint_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    for (shape = 0u; shape <= KPROP_CDS_EDIT_INDEL; shape++) {
        ASSERT(g_delta_exon_hint_cov.shape[shape] > 0u);
    }
    ASSERT(g_delta_exon_hint_cov.fwd > 0u);
    ASSERT(g_delta_exon_hint_cov.rev > 0u);
    fprintf(stderr,
            "[delta-exon-hint coverage] snv=%u mnv=%u ins=%u del=%u "
            "indel=%u fwd=%u rev=%u\n",
            g_delta_exon_hint_cov.shape[KPROP_CDS_EDIT_SNV],
            g_delta_exon_hint_cov.shape[KPROP_CDS_EDIT_MNV],
            g_delta_exon_hint_cov.shape[KPROP_CDS_EDIT_INS],
            g_delta_exon_hint_cov.shape[KPROP_CDS_EDIT_DEL],
            g_delta_exon_hint_cov.shape[KPROP_CDS_EDIT_INDEL],
            g_delta_exon_hint_cov.fwd, g_delta_exon_hint_cov.rev);
    PASS();
}

static struct {
    uint32_t fast;
    uint32_t fallback;
    uint32_t frameshift;
    uint32_t inframe_insertion;
    uint32_t inframe_deletion;
    uint32_t insertion;
    uint32_t deletion;
    uint32_t delins;
    uint32_t forward;
    uint32_t reverse;
} g_simple_indel_equivalence_cov;

/* The optimized route is not a second consequence authority. Force the same
 * event through the generalized CodingContext interpreter by requesting its
 * borrowed context, then compare every public delta field plus the sequence
 * resolution and cached NMD projection facts. The generator includes start,
 * body, and terminal edits on both strands, so excluded states also prove that
 * the dispatcher falls back without changing results. */
static enum theft_trial_res prop_simple_indel_route_matches_coding_context(
    struct theft *t, void *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t authoritative_edits[4];
    duckvep_haplotype_edit_t routed_edits[4];
    uint8_t authoritative_alt_cds[96];
    uint8_t routed_alt_cds[96];
    uint8_t authoritative_ref_peptide[40];
    uint8_t authoritative_alt_peptide[40];
    uint8_t routed_ref_peptide[40];
    uint8_t routed_alt_peptide[40];
    duckvep_delta_scratch_t authoritative_scratch;
    duckvep_delta_scratch_t routed_scratch;
    duckvep_event_t event;
    duckvep_sequence_delta_t authoritative;
    duckvep_sequence_delta_t routed;
    duckvep_coding_context_t context;
    duckvep_variant_coding_context_status_t context_status;
    duckvep_sequence_delta_route_t authoritative_route;
    duckvep_sequence_delta_route_t routed_route;
    (void)t;

    if (s->expect_shape != KPROP_CDS_EDIT_INS &&
        s->expect_shape != KPROP_CDS_EDIT_DEL &&
        s->expect_shape != KPROP_CDS_EDIT_INDEL) {
        return THEFT_TRIAL_PASS;
    }

    memset(&authoritative_scratch, 0, sizeof authoritative_scratch);
    authoritative_scratch.edits = authoritative_edits;
    authoritative_scratch.edits_cap = 4u;
    authoritative_scratch.alt_cds = authoritative_alt_cds;
    authoritative_scratch.alt_cds_cap = sizeof authoritative_alt_cds;
    authoritative_scratch.ref_peptide = authoritative_ref_peptide;
    authoritative_scratch.ref_peptide_cap = sizeof authoritative_ref_peptide;
    authoritative_scratch.alt_peptide = authoritative_alt_peptide;
    authoritative_scratch.alt_peptide_cap = sizeof authoritative_alt_peptide;

    memset(&routed_scratch, 0, sizeof routed_scratch);
    routed_scratch.edits = routed_edits;
    routed_scratch.edits_cap = 4u;
    routed_scratch.alt_cds = routed_alt_cds;
    routed_scratch.alt_cds_cap = sizeof routed_alt_cds;
    routed_scratch.ref_peptide = routed_ref_peptide;
    routed_scratch.ref_peptide_cap = sizeof routed_ref_peptide;
    routed_scratch.alt_peptide = routed_alt_peptide;
    routed_scratch.alt_peptide_cap = sizeof routed_alt_peptide;

    duckvep_event_load(&s->v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_observed(
        (duckvep_variant_kind_t)s->vkind,
        &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->vpos, s->strand,
        &authoritative_scratch, &event, (uint32_t)DUCKVEP_REGION_CDS, 0u, NULL,
        &authoritative_route, &authoritative, &context, &context_status);
    duckvep_sequence_delta_fill_for_annotation_trace(
        (duckvep_variant_kind_t)s->vkind,
        &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->vpos, s->strand,
        &routed_scratch, &event, (uint32_t)DUCKVEP_REGION_CDS, 0u,
        &routed_route, &routed);

    if (authoritative_route == DUCKVEP_DELTA_ROUTE_SIMPLE_INDEL ||
        !kprop_sequence_delta_equal(&authoritative, &routed) ||
        authoritative.sequence_status != routed.sequence_status ||
        authoritative.nmd_early_cds_fact != routed.nmd_early_cds_fact) {
        return THEFT_TRIAL_FAIL;
    }

    if (routed_route != DUCKVEP_DELTA_ROUTE_SIMPLE_INDEL) {
        g_simple_indel_equivalence_cov.fallback++;
        return THEFT_TRIAL_PASS;
    }
    if (context_status != DUCKVEP_VARIANT_CODING_CONTEXT_OK ||
        !routed.valid) {
        return THEFT_TRIAL_FAIL;
    }

    g_simple_indel_equivalence_cov.fast++;
    if (routed.frameshift) g_simple_indel_equivalence_cov.frameshift++;
    if (routed.inframe_insertion) {
        g_simple_indel_equivalence_cov.inframe_insertion++;
    }
    if (routed.inframe_deletion) {
        g_simple_indel_equivalence_cov.inframe_deletion++;
    }
    if (s->expect_shape == KPROP_CDS_EDIT_INS) {
        g_simple_indel_equivalence_cov.insertion++;
    } else if (s->expect_shape == KPROP_CDS_EDIT_DEL) {
        g_simple_indel_equivalence_cov.deletion++;
    } else {
        g_simple_indel_equivalence_cov.delins++;
    }
    if (s->strand > 0) g_simple_indel_equivalence_cov.forward++;
    else g_simple_indel_equivalence_cov.reverse++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_simple_indel_route_matches_generalized_context(void) {
    struct theft_run_config cfg;

    memset(&cfg, 0, sizeof cfg);
    cfg.name = "simple indel route == generalized CodingContext";
    cfg.prop1 = prop_simple_indel_route_matches_coding_context;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_simple_indel_equivalence_cov, 0,
           sizeof g_simple_indel_equivalence_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_simple_indel_equivalence_cov.fast > 0u);
    ASSERT(g_simple_indel_equivalence_cov.fallback > 0u);
    ASSERT(g_simple_indel_equivalence_cov.frameshift > 0u);
    ASSERT(g_simple_indel_equivalence_cov.inframe_insertion > 0u);
    ASSERT(g_simple_indel_equivalence_cov.inframe_deletion > 0u);
    ASSERT(g_simple_indel_equivalence_cov.insertion > 0u);
    ASSERT(g_simple_indel_equivalence_cov.deletion > 0u);
    ASSERT(g_simple_indel_equivalence_cov.delins > 0u);
    ASSERT(g_simple_indel_equivalence_cov.forward > 0u);
    ASSERT(g_simple_indel_equivalence_cov.reverse > 0u);
    fprintf(stderr,
            "[simple-indel equivalence coverage] fast=%u fallback=%u "
            "frameshift=%u inframe_ins=%u inframe_del=%u "
            "ins=%u del=%u delins=%u fwd=%u rev=%u\n",
            g_simple_indel_equivalence_cov.fast,
            g_simple_indel_equivalence_cov.fallback,
            g_simple_indel_equivalence_cov.frameshift,
            g_simple_indel_equivalence_cov.inframe_insertion,
            g_simple_indel_equivalence_cov.inframe_deletion,
            g_simple_indel_equivalence_cov.insertion,
            g_simple_indel_equivalence_cov.deletion,
            g_simple_indel_equivalence_cov.delins,
            g_simple_indel_equivalence_cov.forward,
            g_simple_indel_equivalence_cov.reverse);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
    uint32_t full;
} g_cursor_route_cov;

static enum theft_trial_res prop_cursor_padded_snv_matches_tile(struct theft *t,
                                                                void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_variant_batch_t variants = s->v;
    const uint8_t semantic_kind = (uint8_t)DUCKVEP_KIND_SNV;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    int saw_full = 0;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    /* The corpus stores whole codons whose differing region is one base. At the
     * public boundary these are SNVs, regardless of their padded representation. */
    variants.variant_kind = &semantic_kind;
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &variants, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &variants, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (!saw_full || tile_n != cursor_n || tile_n == 0u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (tile_stats.substitution_context != 1u ||
        tile_stats.substitution_context != cursor_stats.substitution_context) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    {
        size_t i;
        for (i = 0u; i < tile_n; i++) {
            if (!consequence_rows_equal(&tile_rows[i], &cursor_rows[i])) {
                res = THEFT_TRIAL_FAIL; goto done;
            }
        }
    }

    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) g_cursor_route_cov.syn++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) g_cursor_route_cov.mis++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) g_cursor_route_cov.stop_gained++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) g_cursor_route_cov.stop_lost++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) g_cursor_route_cov.stop_retained++;
    else { res = THEFT_TRIAL_FAIL; goto done; }
    if (s->strand > 0) g_cursor_route_cov.fwd++;
    else g_cursor_route_cov.rev++;
    if (saw_full) g_cursor_route_cov.full++;

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_padded_snv_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor padded SNV == tile under output splits";
    cfg.prop1 = prop_cursor_padded_snv_matches_tile;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_route_cov, 0, sizeof g_cursor_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_route_cov.syn > 0u);
    ASSERT(g_cursor_route_cov.mis > 0u);
    ASSERT(g_cursor_route_cov.stop_gained > 0u);
    ASSERT(g_cursor_route_cov.stop_lost > 0u);
    ASSERT(g_cursor_route_cov.stop_retained > 0u);
    ASSERT(g_cursor_route_cov.fwd > 0u);
    ASSERT(g_cursor_route_cov.rev > 0u);
    ASSERT(g_cursor_route_cov.full > 0u);
    fprintf(stderr,
            "[cursor-route coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u full=%u\n",
            g_cursor_route_cov.syn, g_cursor_route_cov.mis,
            g_cursor_route_cov.stop_gained, g_cursor_route_cov.stop_lost,
            g_cursor_route_cov.stop_retained, g_cursor_route_cov.fwd,
            g_cursor_route_cov.rev, g_cursor_route_cov.full);
    PASS();
}

static struct {
    uint32_t context;
    uint32_t fwd;
    uint32_t rev;
    uint32_t len2;
    uint32_t len3;
} g_cursor_cross_route_cov;

static enum theft_trial_res prop_cursor_cross_codon_mnv_route_matches_tile(
    struct theft *t, void *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (tile_n != cursor_n || tile_n != 1u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (!consequence_rows_equal(&tile_rows[0], &cursor_rows[0])) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (tile_stats.substitution_context != cursor_stats.substitution_context) {
        res = THEFT_TRIAL_FAIL; goto done;
    }

    /* Every generator mode (missense / synonymous / stop-gained cross-codon window)
     * resolves through the context interpreter. The
     * property under test is that a chunked cursor and a single tile agree on both the row
     * and the route stats for any output split — verified above. */
    if (tile_stats.substitution_context != 1u) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    g_cursor_cross_route_cov.context++;
    if (s->strand > 0) g_cursor_cross_route_cov.fwd++; else g_cursor_cross_route_cov.rev++;
    if (s->rlen == 2u) g_cursor_cross_route_cov.len2++;
    else if (s->rlen == 3u) g_cursor_cross_route_cov.len3++;
    else { res = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_cross_codon_mnv_route_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor cross-codon MNV route == tile";
    cfg.prop1 = prop_cursor_cross_codon_mnv_route_matches_tile;
    cfg.type_info[0] = &kprop_cross_codon_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_cross_route_cov, 0, sizeof g_cursor_cross_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_cross_route_cov.context > 0u);
    ASSERT(g_cursor_cross_route_cov.fwd > 0u);
    ASSERT(g_cursor_cross_route_cov.rev > 0u);
    ASSERT(g_cursor_cross_route_cov.len2 > 0u);
    ASSERT(g_cursor_cross_route_cov.len3 > 0u);
    fprintf(stderr,
            "[cursor-cross-route coverage] context=%u fwd=%u rev=%u len2=%u len3=%u\n",
            g_cursor_cross_route_cov.context,
            g_cursor_cross_route_cov.fwd, g_cursor_cross_route_cov.rev,
            g_cursor_cross_route_cov.len2, g_cursor_cross_route_cov.len3);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
    uint32_t capfail;
} g_delta_scratch_cov;

static enum theft_trial_res prop_delta_scratch_mnv_matches_codon_oracle(struct theft *t,
                                                                        void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_delta_scratch_t too_small;
    duckvep_sequence_delta_t delta;
    size_t first_diff = 0u;
    size_t codon_start;
    size_t codon_idx;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    uint32_t j;
    (void)t;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    too_small = scratch;
    /* Full-feature MNVs borrow the uploaded allele, not a materialized edit
     * array. Exhaust the local ALT storage that this producer actually uses. */
    too_small.alt_cds_cap = 0u;
    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s->vkind,
                                             &s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->vpos, s->strand,
                                             &too_small, &delta);
    if (delta.valid) return THEFT_TRIAL_FAIL;
    g_delta_scratch_cov.capfail++;

    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s->vkind,
                                             &s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->vpos, s->strand,
                                             &scratch, &delta);
    if (!delta.valid) return THEFT_TRIAL_FAIL;

    while (first_diff < s->cds_lenv && s->cds[first_diff] == s->expect_cds[first_diff]) {
        first_diff++;
    }
    if (first_diff >= s->cds_lenv) return THEFT_TRIAL_FAIL;
    codon_start = first_diff - (first_diff % 3u);
    codon_idx = codon_start / 3u;
    for (j = 0u; j < 3u; j++) {
        ref_codon[j] = (char)s->cds[codon_start + (size_t)j];
        alt_codon[j] = (char)s->expect_cds[codon_start + (size_t)j];
    }
    ref_codon[3] = '\0';
    alt_codon[3] = '\0';
    cr = duckvep_codon_change(ref_codon, alt_codon, DUCKVEP_CODON_TABLE_STANDARD);
    if (cr.change & DUCKVEP_CODON_INVALID) return THEFT_TRIAL_FAIL;
    if (delta.protein_pos != (int32_t)(codon_idx + 1u) ||
        delta.ref_aa != (uint8_t)cr.aa_ref || delta.alt_aa != (uint8_t)cr.aa_alt ||
        delta.cdna_pos != -1 || delta.cds_pos != -1) {
        return THEFT_TRIAL_FAIL;
    }
    if (delta.start_lost || delta.start_retained || delta.frameshift ||
        delta.inframe_deletion ||
        delta.inframe_insertion || delta.protein_altering || delta.coding_unknown ||
        delta.partial_codon) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) {
        if (!delta.synonymous || delta.missense || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.syn++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) {
        if (!delta.missense || delta.synonymous || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.mis++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) {
        if (!delta.stop_gained || delta.synonymous || delta.missense ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.stop_gained++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) {
        if (!delta.stop_lost || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.stop_lost++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) {
        if (!delta.stop_retained || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_lost) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.stop_retained++;
    } else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_delta_scratch_cov.fwd++;
    else g_delta_scratch_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_with_scratch_mnv_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta scratch MNV == single-codon oracle";
    cfg.prop1 = prop_delta_scratch_mnv_matches_codon_oracle;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_scratch_cov, 0, sizeof g_delta_scratch_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_scratch_cov.syn > 0u);
    ASSERT(g_delta_scratch_cov.mis > 0u);
    ASSERT(g_delta_scratch_cov.stop_gained > 0u);
    ASSERT(g_delta_scratch_cov.stop_lost > 0u);
    ASSERT(g_delta_scratch_cov.stop_retained > 0u);
    ASSERT(g_delta_scratch_cov.fwd > 0u);
    ASSERT(g_delta_scratch_cov.rev > 0u);
    ASSERT(g_delta_scratch_cov.capfail > 0u);
    fprintf(stderr,
            "[delta-scratch coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u capfail=%u\n",
            g_delta_scratch_cov.syn, g_delta_scratch_cov.mis,
            g_delta_scratch_cov.stop_gained, g_delta_scratch_cov.stop_lost,
            g_delta_scratch_cov.stop_retained, g_delta_scratch_cov.fwd,
            g_delta_scratch_cov.rev, g_delta_scratch_cov.capfail);
    PASS();
}

static struct { uint32_t missense; uint32_t synonymous; uint32_t stop_gained; uint32_t fwd; uint32_t rev; uint32_t len2; uint32_t len3; } g_delta_cross_scratch_cov;

static enum theft_trial_res prop_delta_scratch_cross_codon_mnv_matches_oracle(
    struct theft *t, void *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    struct kprop_cc_facts f;
    (void)t;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s->vkind,
                                             &s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->vpos, s->strand,
                                             &scratch, &direct);
    duckvep_sequence_delta_fill((duckvep_variant_kind_t)s->vkind,
                                &s->tx, &s->ex, &s->seq, &s->v,
                                0u, 0u, s->vpos, s->strand, &shape);
    duckvep_sequence_delta_fill_for_annotation_trace((duckvep_variant_kind_t)s->vkind,
                                                     &s->tx, &s->ex, &s->seq, &s->v,
                                                     0u, 0u, s->vpos, s->strand,
                                                     &scratch, NULL, UINT32_MAX,
                                                     UINT32_MAX,
                                                     &route, &routed);

    /* The generator only builds two adjacent body codons, so the window oracle always
     * resolves; the authoritative interpreter emits it directly and the router echoes it. */
    f = kprop_cross_codon_mnv_oracle_facts(s);
    if (!f.valid) return THEFT_TRIAL_FAIL;
    if (route != DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT) return THEFT_TRIAL_FAIL;
    if (!direct.valid || !routed.valid) return THEFT_TRIAL_FAIL;
    if (!kprop_sequence_delta_equal(&direct, &routed)) return THEFT_TRIAL_FAIL;
    /* Coarse multi-codon window: protein_pos -1, no AA pair. */
    if (direct.protein_pos != -1 || direct.cdna_pos != -1 || direct.cds_pos != -1 ||
        direct.ref_aa != (uint8_t)0u || direct.alt_aa != (uint8_t)0u) {
        return THEFT_TRIAL_FAIL;
    }
    /* SO facts equal the independent window oracle. */
    if (direct.synonymous != (uint8_t)f.synonymous ||
        direct.missense != (uint8_t)f.missense ||
        direct.stop_gained != (uint8_t)f.stop_gained ||
        direct.stop_lost != (uint8_t)f.stop_lost ||
        direct.stop_retained != (uint8_t)f.stop_retained ||
        direct.start_lost || direct.start_retained || direct.frameshift ||
        direct.inframe_deletion ||
        direct.inframe_insertion || direct.protein_altering || direct.coding_unknown ||
        direct.partial_codon) {
        return THEFT_TRIAL_FAIL;
    }
    /* The narrow direct reference only resolves the missense cross-codon window; the
     * synonymous and stop-gained windows are exactly what the authoritative interpreter adds. */
    if (f.missense) {
        if (!shape.valid || !kprop_delta_is_coarse_cross_codon_missense(&shape)) {
            return THEFT_TRIAL_FAIL;
        }
        g_delta_cross_scratch_cov.missense++;
    } else if (f.synonymous) {
        if (shape.valid) return THEFT_TRIAL_FAIL;
        g_delta_cross_scratch_cov.synonymous++;
    } else if (f.stop_gained) {
        if (shape.valid) return THEFT_TRIAL_FAIL;
        g_delta_cross_scratch_cov.stop_gained++;
    } else {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_delta_cross_scratch_cov.fwd++; else g_delta_cross_scratch_cov.rev++;
    if (s->rlen == 2u) g_delta_cross_scratch_cov.len2++;
    else if (s->rlen == 3u) g_delta_cross_scratch_cov.len3++;
    else return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_with_scratch_cross_codon_mnv_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta scratch two-codon MNV window == codon-window oracle";
    cfg.prop1 = prop_delta_scratch_cross_codon_mnv_matches_oracle;
    cfg.type_info[0] = &kprop_cross_codon_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_cross_scratch_cov, 0, sizeof g_delta_cross_scratch_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_cross_scratch_cov.missense > 0u);
    ASSERT(g_delta_cross_scratch_cov.synonymous > 0u);
    ASSERT(g_delta_cross_scratch_cov.stop_gained > 0u);
    ASSERT(g_delta_cross_scratch_cov.fwd > 0u);
    ASSERT(g_delta_cross_scratch_cov.rev > 0u);
    ASSERT(g_delta_cross_scratch_cov.len2 > 0u);
    ASSERT(g_delta_cross_scratch_cov.len3 > 0u);
    fprintf(stderr,
            "[delta-cross-scratch coverage] missense=%u synonymous=%u stop_gained=%u fwd=%u rev=%u len2=%u len3=%u\n",
            g_delta_cross_scratch_cov.missense, g_delta_cross_scratch_cov.synonymous,
            g_delta_cross_scratch_cov.stop_gained, g_delta_cross_scratch_cov.fwd,
            g_delta_cross_scratch_cov.rev,
            g_delta_cross_scratch_cov.len2, g_delta_cross_scratch_cov.len3);
    PASS();
}

static struct {
    uint32_t ins;
    uint32_t del;
    uint32_t delins;
    uint32_t delins_plus1;
    uint32_t delins_plus2;
    uint32_t delins_minus1;
    uint32_t delins_minus2;
    uint32_t rev;
    uint32_t stop_gained;
    uint32_t terminal_stop;
    uint32_t terminal_nonstop;
    uint32_t terminal_reverse;
    uint32_t terminal_missing_tail;
    uint32_t terminal_cil_retained;
    uint32_t terminal_cil_protein_altering;
} g_frameshift_cov;

static int kprop_pep_window_prefix_or_suffix(
    const uint8_t *outer, uint32_t outer_off, uint32_t outer_len,
    const uint8_t *inner, uint32_t inner_off, uint32_t inner_len);

/* Rebuild VEP's local peptide alleles from the public coding-context fields.
 * This is deliberately a test-side coordinate derivation: it rounds the
 * uploaded CDS span to codons, applies the net length change, and spells a
 * trailing partial codon X. */
static int kprop_vep_local_peptides(
    const duckvep_coding_context_t *ctx,
    uint8_t ref_local[32], size_t *ref_local_len,
    uint8_t alt_local[32], size_t *alt_local_len) {

    uint64_t first;
    uint64_t translation_start;
    uint64_t translation_end;
    uint64_t ref_codons;
    uint64_t nt_offset;
    uint64_t ref_nt_length;
    int64_t requested_alt_nt;
    size_t alt_nt_length;
    size_t ref_whole;
    size_t alt_whole;
    size_t ref_peptide_offset;
    size_t i;

    if (ctx == NULL || ref_local_len == NULL || alt_local_len == NULL ||
        ctx->single_edit_cds_start == 0u || ctx->length_diff == INT64_MIN) {
        return 0;
    }
    first = (uint64_t)ctx->single_edit_cds_start;
    translation_start = ((first - 1u) / 3u) + 1u;
    if (ctx->single_edit_ref_len != 0u) {
        uint64_t last = first + (uint64_t)ctx->single_edit_ref_len - 1u;
        if (last < first) return 0;
        translation_end = ((last - 1u) / 3u) + 1u;
        ref_codons = translation_end - translation_start + 1u;
    } else if (first > 1u) {
        uint64_t last = first - 1u;
        translation_end = ((last - 1u) / 3u) + 1u;
        ref_codons = translation_end >= translation_start
            ? translation_end - translation_start + 1u : 0u;
    } else {
        ref_codons = 0u;
    }
    nt_offset = (translation_start - 1u) * 3u;
    ref_nt_length = ref_codons * 3u;
    if (nt_offset > (uint64_t)ctx->ref_cds_len) return 0;
    if (ref_nt_length > (uint64_t)ctx->ref_cds_len - nt_offset) {
        ref_nt_length = (uint64_t)ctx->ref_cds_len - nt_offset;
    }
    if (ref_nt_length > (uint64_t)INT64_MAX ||
        (ctx->length_diff > 0 &&
         ctx->length_diff > INT64_MAX - (int64_t)ref_nt_length)) {
        return 0;
    }
    requested_alt_nt = (int64_t)ref_nt_length + ctx->length_diff;
    if (requested_alt_nt < 0 || nt_offset > (uint64_t)ctx->alt_cds_len) {
        return 0;
    }
    alt_nt_length = (size_t)requested_alt_nt;
    if (alt_nt_length > ctx->alt_cds_len - (size_t)nt_offset) {
        alt_nt_length = ctx->alt_cds_len - (size_t)nt_offset;
    }
    ref_whole = (size_t)ref_nt_length / 3u;
    alt_whole = alt_nt_length / 3u;
    ref_peptide_offset = (size_t)nt_offset / 3u;
    if (ref_peptide_offset > ctx->ref_peptide_len ||
        ref_whole > ctx->ref_peptide_len - ref_peptide_offset ||
        ref_peptide_offset > ctx->alt_peptide_len ||
        alt_whole > ctx->alt_peptide_len - ref_peptide_offset ||
        ref_whole + 1u > 32u || alt_whole + 1u > 32u) {
        return 0;
    }
    for (i = 0u; i < ref_whole; i++) {
        ref_local[i] = ctx->ref_peptide[ref_peptide_offset + i];
    }
    for (i = 0u; i < alt_whole; i++) {
        alt_local[i] = ctx->alt_peptide[ref_peptide_offset + i];
    }
    *ref_local_len = ref_whole;
    *alt_local_len = alt_whole;
    if ((ref_nt_length % 3u) != 0u &&
        !(ref_whole == 1u && ref_local[0] == (uint8_t)'*')) {
        ref_local[(*ref_local_len)++] = (uint8_t)'X';
    }
    if ((alt_nt_length % 3u) != 0u &&
        !(alt_whole == 1u && alt_local[0] == (uint8_t)'*')) {
        alt_local[(*alt_local_len)++] = (uint8_t)'X';
    }
    return 1;
}

/* VariationEffect.pm::stop_retained only reaches
 * _ins_del_stop_altered_cil when the codon-local alternate peptide is absent,
 * empty, or contains X. A concrete peptide, including the literal value "*",
 * stays on ref_eq_alt_sequence even when insertion length reaches the original
 * terminal codon. Return -1 only when the local view cannot be reconstructed. */
static int kprop_vep_stop_retained_uses_cil(
    const duckvep_coding_context_t *ctx) {
    uint8_t ref_local[32];
    uint8_t alt_local[32];
    size_t ref_local_len;
    size_t alt_local_len;
    size_t i;

    if (!kprop_vep_local_peptides(
            ctx, ref_local, &ref_local_len, alt_local, &alt_local_len)) {
        return -1;
    }
    (void)ref_local_len;
    if (alt_local_len == 0u) return 1;
    for (i = 0u; i < alt_local_len; i++) {
        if (alt_local[i] == (uint8_t)'X') return 1;
    }
    return 0;
}

static enum theft_trial_res prop_annotate_frameshift_indel_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_haplotype_edit_t edit_scratch[8];
    duckvep_coding_context_t fsctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint32_t protein_cds;
    uint32_t j;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    /* Production may reuse the sequence classifier's first-101-CDS-base fact.
     * Compare its public NMD result with the independent exhaustive projector
     * for every randomized insertion, deletion, delins, and strand. */
    {
        duckvep_event_t event;
        duckvep_nmd_result_t exhaustive;

        duckvep_event_load(&s->v, 0u, &event);
        duckvep_nmd_predict(&s->tx, &s->ex, 0u, &event,
                            rows[0].consequence_mask, NULL, &exhaustive);
        if (rows[0].nmd_prediction != exhaustive.prediction ||
            rows[0].nmd_escape_reasons != exhaustive.escape_reasons) {
            tr = THEFT_TRIAL_FAIL;
            goto done;
        }
    }

    if (s->vkind == (uint8_t)DUCKVEP_KIND_INS) {
        protein_cds = kprop_cds_pos_for_genomic(s, s->vpos);
        g_frameshift_cov.ins++;
    } else if (s->vkind == (uint8_t)DUCKVEP_KIND_DEL) {
        protein_cds = UINT32_MAX;
        for (j = 1u; j < (uint32_t)s->rlen; j++) {
            uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
            if (cds_pos < protein_cds) protein_cds = cds_pos;
        }
        if (protein_cds == UINT32_MAX) { tr = THEFT_TRIAL_FAIL; goto done; }
        g_frameshift_cov.del++;
    } else if (s->vkind == (uint8_t)DUCKVEP_KIND_INDEL) {
        protein_cds = UINT32_MAX;
        for (j = 0u; j < (uint32_t)s->rlen; j++) {
            uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
            if (cds_pos < protein_cds) protein_cds = cds_pos;
        }
        if (protein_cds == UINT32_MAX) { tr = THEFT_TRIAL_FAIL; goto done; }
        g_frameshift_cov.delins++;
        if (s->alen == s->rlen + 1u) g_frameshift_cov.delins_plus1++;
        else if (s->alen == s->rlen + 2u) g_frameshift_cov.delins_plus2++;
        else if (s->rlen == s->alen + 1u) g_frameshift_cov.delins_minus1++;
        else if (s->rlen == s->alen + 2u) g_frameshift_cov.delins_minus2++;
    } else { tr = THEFT_TRIAL_FAIL; goto done; }
    if (s->strand < 0) g_frameshift_cov.rev++;

    /* Frameshift fact plus VEP's local-window stop composite. Rebuild the same coding
     * context the engine uses, then independently re-derive stop_gained from the CDS bytes
     * (kprop_frameshift_local_stop_oracle) for the standard table. For the mitochondrial
     * table the standard-stop enumeration does not apply, so require the frameshift bit and
     * allow at most the stop composite. */
    {
        int have_ctx = duckvep_variant_physical_coding_context_build(
                           &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand, edit_scratch, 8u,
                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                           alt_pep, sizeof alt_pep, &fsctx) ==
                       DUCKVEP_VARIANT_CODING_CONTEXT_OK;
        int overlaps_terminal_stop = 0;
        int reaches_terminal_cil = 0;
        uint32_t terminal_start = 0u;
        if (have_ctx && fsctx.ref_cds_len >= 3u) {
            terminal_start = (uint32_t)fsctx.ref_cds_len - 2u;
        }
        if (have_ctx && terminal_start > 0u) {
            if (fsctx.single_edit_ref_len == 0u) {
                overlaps_terminal_stop =
                    fsctx.single_edit_cds_start > terminal_start &&
                    fsctx.single_edit_cds_start <= fsctx.ref_cds_len;
            } else {
                uint64_t edit_end = (uint64_t)fsctx.single_edit_cds_start +
                    (uint64_t)fsctx.single_edit_ref_len - 1u;
                overlaps_terminal_stop =
                    fsctx.single_edit_cds_start <= fsctx.ref_cds_len &&
                    edit_end >= (uint64_t)terminal_start;
            }
            reaches_terminal_cil =
                fsctx.insertion_length_reaches_terminal_stop != 0u;
        }
        if (overlaps_terminal_stop) {
            if (fsctx.alt_cds_len >= fsctx.ref_cds_len) {
                char codon[3];
                char alt_aa;
                int ref_local_begins_stop;
                uint64_t want;
                memcpy(codon, fsctx.alt_cds + fsctx.ref_cds_len - 3u,
                       sizeof codon);
                alt_aa = duckvep_translate_codon(
                    codon, DUCKVEP_CODON_TABLE_STANDARD);
                memcpy(codon, fsctx.ref_cds + fsctx.ref_cds_len - 3u,
                       sizeof codon);
                ref_local_begins_stop =
                    fsctx.single_edit_cds_start >= terminal_start &&
                    duckvep_translate_codon(
                        codon, DUCKVEP_CODON_TABLE_STANDARD) == '*';
                if (alt_aa == '*') {
                    want = DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
                } else if (ref_local_begins_stop) {
                    want = DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
                } else {
                    want = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                           DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
                }
                if (rows[0].consequence_mask != want ||
                    rows[0].sequence_status !=
                        (uint8_t)DUCKVEP_SEQUENCE_RESOLVED) {
                    fprintf(stderr,
                            "[terminal endpoint mismatch] kind=%u strand=%d "
                            "cds_start=%u ref_len=%u alt_len=%u ref_cds=%zu "
                            "alt_cds=%zu ref_stop=%d alt_aa=%c want=%" PRIu64
                            " got=%" PRIu64 " status=%u\n",
                            (unsigned)s->vkind, (int)s->strand,
                            fsctx.single_edit_cds_start,
                            fsctx.single_edit_ref_len,
                            fsctx.single_edit_alt_len,
                            fsctx.ref_cds_len, fsctx.alt_cds_len,
                            ref_local_begins_stop, alt_aa,
                            want, rows[0].consequence_mask,
                            (unsigned)rows[0].sequence_status);
                    tr = THEFT_TRIAL_FAIL;
                } else if ((want & DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)) == 0u &&
                           (rows[0].protein_pos != (int32_t)fsctx.ref_peptide_len ||
                            rows[0].aa_ref != (uint8_t)'*' ||
                            rows[0].aa_alt != (uint8_t)alt_aa)) {
                    fprintf(stderr,
                            "[terminal endpoint payload mismatch] kind=%u "
                            "strand=%d cds_start=%u ref_len=%u alt_len=%u "
                            "ref_cds=%zu alt_cds=%zu want=%" PRIu64
                            " protein=%d/%zu aa=%c/%c expected_alt=%c\n",
                            (unsigned)s->vkind, (int)s->strand,
                            fsctx.single_edit_cds_start,
                            fsctx.single_edit_ref_len,
                            fsctx.single_edit_alt_len,
                            fsctx.ref_cds_len, fsctx.alt_cds_len, want,
                            rows[0].protein_pos, fsctx.ref_peptide_len,
                            (char)rows[0].aa_ref, (char)rows[0].aa_alt,
                            alt_aa);
                    tr = THEFT_TRIAL_FAIL;
                } else {
                    g_frameshift_cov.terminal_stop++;
                    if (!ref_local_begins_stop) {
                        g_frameshift_cov.terminal_nonstop++;
                    }
                    if (s->strand < 0) g_frameshift_cov.terminal_reverse++;
                }
            } else {
                if (rows[0].consequence_mask !=
                        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) ||
                    rows[0].sequence_status !=
                        (uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL ||
                    rows[0].protein_pos != -1) {
                    tr = THEFT_TRIAL_FAIL;
                } else {
                    g_frameshift_cov.terminal_missing_tail++;
                }
            }
            goto done;
        }
        if (reaches_terminal_cil &&
            fsctx.alt_cds_len >= fsctx.ref_cds_len) {
            char codon[3];
            char endpoint_aa;
            int stop_retained_uses_cil;

            memcpy(codon, fsctx.alt_cds + fsctx.ref_cds_len - 3u,
                   sizeof codon);
            endpoint_aa = duckvep_translate_codon(
                codon, DUCKVEP_CODON_TABLE_STANDARD);
            stop_retained_uses_cil =
                kprop_vep_stop_retained_uses_cil(&fsctx);
            if (stop_retained_uses_cil < 0) {
                tr = THEFT_TRIAL_FAIL;
                goto done;
            }
            if (endpoint_aa == '*' && stop_retained_uses_cil != 0) {
                uint8_t ref_local[32];
                uint8_t alt_local[32];
                size_t ref_local_len;
                size_t alt_local_len;
                size_t alt_shape_len;
                size_t i;
                int inframe;
                uint64_t want = DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);

                if (!kprop_vep_local_peptides(
                        &fsctx, ref_local, &ref_local_len,
                        alt_local, &alt_local_len)) {
                    tr = THEFT_TRIAL_FAIL;
                    goto done;
                }
                alt_shape_len = alt_local_len;
                for (i = 0u; i < alt_local_len; i++) {
                    if (alt_local[i] == (uint8_t)'*') {
                        alt_shape_len = i + 1u;
                        break;
                    }
                }
                inframe = ref_local_len == 0u ||
                    kprop_pep_window_prefix_or_suffix(
                        alt_local, 0u, (uint32_t)alt_shape_len,
                        ref_local, 0u, (uint32_t)ref_local_len);
                if (inframe) {
                    want |= DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION);
                }
                /* VEP evaluates protein_altering_variant independently after
                 * stop_retained has suppressed frameshift.  For an insertion
                 * away from the start codon, unequal local peptides that
                 * preserve neither edge therefore retain the stop and are
                 * also protein-altering.  Do not clean up that odd but real
                 * predicate combination in the oracle. */
                if (fsctx.single_edit_cds_start > 3u &&
                    ref_local_len != 0u && alt_local_len != 0u &&
                    ref_local_len != alt_local_len &&
                    ref_local[0] != (uint8_t)'*' &&
                    alt_local[0] != (uint8_t)'*' &&
                    !kprop_pep_window_prefix_or_suffix(
                        alt_local, 0u, (uint32_t)alt_local_len,
                        ref_local, 0u, (uint32_t)ref_local_len)) {
                    want |= DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
                }
                if (rows[0].consequence_mask != want ||
                    rows[0].sequence_status !=
                        (uint8_t)DUCKVEP_SEQUENCE_RESOLVED) {
                    fprintf(stderr,
                            "[terminal CIL mismatch] kind=%u strand=%d "
                            "genomic=%u-%u raw_ref=%.*s raw_alt=%.*s "
                            "cds_start=%u ref_len=%u alt_len=%u ref_cds=%zu "
                            "alt_cds=%zu local_ref=%.*s local_alt=%.*s "
                            "inframe=%d "
                            "cds=%.*s ref_pep=%.*s alt_pep=%.*s want=%" PRIu64
                            " got=%" PRIu64 " status=%u\n",
                            (unsigned)s->vkind, (int)s->strand,
                            s->vpos, s->vend,
                            (int)s->rlen, s->abytes + s->roff,
                            (int)s->alen, s->abytes + s->aoff,
                            fsctx.single_edit_cds_start,
                            fsctx.single_edit_ref_len,
                            fsctx.single_edit_alt_len,
                            fsctx.ref_cds_len, fsctx.alt_cds_len,
                            (int)ref_local_len, ref_local,
                            (int)alt_local_len, alt_local, inframe,
                            (int)s->cds_lenv, s->cds,
                            (int)fsctx.ref_peptide_len, fsctx.ref_peptide,
                            (int)fsctx.alt_peptide_len, fsctx.alt_peptide, want,
                            rows[0].consequence_mask,
                            (unsigned)rows[0].sequence_status);
                    tr = THEFT_TRIAL_FAIL;
                } else {
                    memcpy(codon,
                           fsctx.ref_cds + fsctx.ref_cds_len - 3u,
                           sizeof codon);
                    if (duckvep_translate_codon(
                            codon, DUCKVEP_CODON_TABLE_STANDARD) != '*') {
                        g_frameshift_cov.terminal_nonstop++;
                    }
                    g_frameshift_cov.terminal_cil_retained++;
                    if ((want & DUCKVEP_SO(
                            DUCKVEP_SO_PROTEIN_ALTERING)) != 0u) {
                        g_frameshift_cov.terminal_cil_protein_altering++;
                    }
                }
                goto done;
            }
        }
        if (have_ctx && (duckvep_codon_table_t)s->ctab == DUCKVEP_CODON_TABLE_STANDARD) {
            int expected_stop = kprop_frameshift_local_stop_oracle(
                fsctx.ref_cds, fsctx.ref_cds_len, fsctx.alt_cds, fsctx.alt_cds_len,
                fsctx.single_edit_cds_start, fsctx.single_edit_ref_len, fsctx.length_diff);
            uint64_t want = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                            (expected_stop ? DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED) : 0u);
            if (rows[0].consequence_mask != want) {
                fprintf(stderr,
                        "[frameshift body mismatch] kind=%u strand=%d "
                        "cds_start=%u ref_len=%u alt_len=%u ref_cds=%zu "
                        "alt_cds=%zu expected_stop=%d want=%" PRIu64
                        " got=%" PRIu64 " status=%u\n",
                        (unsigned)s->vkind, (int)s->strand,
                        fsctx.single_edit_cds_start,
                        fsctx.single_edit_ref_len,
                        fsctx.single_edit_alt_len,
                        fsctx.ref_cds_len, fsctx.alt_cds_len,
                        expected_stop, want, rows[0].consequence_mask,
                        (unsigned)rows[0].sequence_status);
                tr = THEFT_TRIAL_FAIL;
                goto done;
            }
            if (expected_stop) g_frameshift_cov.stop_gained++;
        } else {
            uint64_t allowed = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                               DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
            if ((rows[0].consequence_mask & DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)) == 0u ||
                (rows[0].consequence_mask & ~allowed) != 0u) { tr = THEFT_TRIAL_FAIL; goto done; }
        }
        if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
        /* Frameshift protein_position = the first affected codon. The authoritative path
         * anchors it to the edit's CDS start (single_edit_cds_start) — the first shifted base,
         * which is VEP's frameshift protein position and, for an insertion at a codon boundary,
         * differs from the anchor base's codon. The edit's cds_start is itself validated
         * against an independent oracle by the cds-edit-builder property tests; when the
         * standalone context build is unavailable we fall back to the genomic-derived codon. */
        {
            int32_t want_pp = have_ctx
                ? (int32_t)(((fsctx.single_edit_cds_start - 1u) / 3u) + 1u)
                : (int32_t)(((protein_cds - 1u) / 3u) + 1u);
            if (rows[0].protein_pos != want_pp) { tr = THEFT_TRIAL_FAIL; goto done; }
        }
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_frameshift_indel_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile simple frameshift indel == CDS-position oracle";
    cfg.prop1 = prop_annotate_frameshift_indel_matches_oracle;
    cfg.type_info[0] = &kprop_frameshift_indel_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_frameshift_cov, 0, sizeof g_frameshift_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_frameshift_cov.ins > 0u);
    ASSERT(g_frameshift_cov.del > 0u);
    ASSERT(g_frameshift_cov.delins > 0u);
    ASSERT(g_frameshift_cov.delins_plus1 > 0u);
    ASSERT(g_frameshift_cov.delins_plus2 > 0u);
    ASSERT(g_frameshift_cov.delins_minus1 > 0u);
    ASSERT(g_frameshift_cov.delins_minus2 > 0u);
    ASSERT(g_frameshift_cov.rev > 0u);
    ASSERT(g_frameshift_cov.terminal_stop > 0u);
    ASSERT(g_frameshift_cov.terminal_nonstop > 0u);
    ASSERT(g_frameshift_cov.terminal_reverse > 0u);
    ASSERT(g_frameshift_cov.terminal_missing_tail > 0u);
    ASSERT(g_frameshift_cov.terminal_cil_retained > 0u);
    fprintf(stderr,
            "[frameshift coverage] ins=%u del=%u delins=%u(+1=%u +2=%u -1=%u -2=%u) reverse=%u stop_gained=%u terminal_endpoint=%u terminal_nonstop=%u terminal_reverse=%u terminal_missing_tail=%u terminal_cil_retained=%u terminal_cil_protein_altering=%u\n",
            g_frameshift_cov.ins, g_frameshift_cov.del, g_frameshift_cov.delins,
            g_frameshift_cov.delins_plus1, g_frameshift_cov.delins_plus2,
            g_frameshift_cov.delins_minus1, g_frameshift_cov.delins_minus2,
            g_frameshift_cov.rev, g_frameshift_cov.stop_gained,
            g_frameshift_cov.terminal_stop,
            g_frameshift_cov.terminal_nonstop,
            g_frameshift_cov.terminal_reverse,
            g_frameshift_cov.terminal_missing_tail,
            g_frameshift_cov.terminal_cil_retained,
            g_frameshift_cov.terminal_cil_protein_altering);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_inframe_deletion_cov;

static enum theft_trial_res prop_annotate_inframe_deletion_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    uint32_t min_cds = UINT32_MAX;
    uint32_t j;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    duckvep_workspace_delta_route_stats_reset(ws);

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    for (j = 1u; j < (uint32_t)s->rlen; j++) {
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
        if (cds_pos < min_cds) min_cds = cds_pos;
    }
    if (min_cds == UINT32_MAX || ((min_cds - 1u) % 3u) != 0u) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (s->strand > 0) g_inframe_deletion_cov.fwd++; else g_inframe_deletion_cov.rev++;

    if (rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)(((min_cds - 1u) / 3u) + 1u)) { tr = THEFT_TRIAL_FAIL; goto done; }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL || stats->simple_indel != 1u ||
        stats->del_context != 0u ||
        stats->substitution_context != 0u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_inframe_deletion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile codon-aligned in-frame deletion == CDS-position oracle";
    cfg.prop1 = prop_annotate_inframe_deletion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_deletion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_inframe_deletion_cov, 0, sizeof g_inframe_deletion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_inframe_deletion_cov.fwd > 0u);
    ASSERT(g_inframe_deletion_cov.rev > 0u);
    fprintf(stderr, "[inframe_deletion coverage] forward=%u reverse=%u\n",
            g_inframe_deletion_cov.fwd, g_inframe_deletion_cov.rev);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; uint32_t full; } g_cursor_del_route_cov;

static enum theft_trial_res prop_cursor_del_route_matches_tile(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    uint32_t min_cds = UINT32_MAX;
    uint32_t j;
    int saw_full = 0;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (!saw_full || tile_n != cursor_n || tile_n != 1u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (tile_stats.simple_indel != 1u ||
        tile_stats.del_context != 0u ||
        tile_stats.substitution_context != 0u ||
        tile_stats.simple_indel != cursor_stats.simple_indel ||
        tile_stats.del_context != cursor_stats.del_context) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (!consequence_rows_equal(&tile_rows[0], &cursor_rows[0])) { res = THEFT_TRIAL_FAIL; goto done; }
    for (j = 1u; j < (uint32_t)s->rlen; j++) {
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
        if (cds_pos < min_cds) min_cds = cds_pos;
    }
    if (min_cds == UINT32_MAX || tile_rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION) ||
        tile_rows[0].protein_pos != (int32_t)(((min_cds - 1u) / 3u) + 1u)) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (s->strand > 0) g_cursor_del_route_cov.fwd++;
    else g_cursor_del_route_cov.rev++;
    if (saw_full) g_cursor_del_route_cov.full++;

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_del_route_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor DEL route == tile under output splits";
    cfg.prop1 = prop_cursor_del_route_matches_tile;
    cfg.type_info[0] = &kprop_inframe_deletion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_del_route_cov, 0, sizeof g_cursor_del_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_del_route_cov.fwd > 0u);
    ASSERT(g_cursor_del_route_cov.rev > 0u);
    ASSERT(g_cursor_del_route_cov.full > 0u);
    fprintf(stderr, "[cursor-del-route coverage] forward=%u reverse=%u full=%u\n",
            g_cursor_del_route_cov.fwd, g_cursor_del_route_cov.rev,
            g_cursor_del_route_cov.full);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_context_inframe_deletion_cov;

static enum theft_trial_res prop_context_inframe_deletion_matches_oracle(struct theft *t,
                                                                        void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    uint32_t min_cds = UINT32_MAX;
    uint32_t j;
    (void)t;

    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!ctx.has_single_edit || ctx.single_edit_alt_len != 0u ||
        ctx.single_edit_ref_len != 3u || ctx.applied_edits != 1u) {
        return THEFT_TRIAL_FAIL;
    }
    for (j = 1u; j < (uint32_t)s->rlen; j++) {
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
        if (cds_pos < min_cds) min_cds = cds_pos;
    }
    if (min_cds == UINT32_MAX) return THEFT_TRIAL_FAIL;
    if (ctx.single_edit_cds_start != min_cds) return THEFT_TRIAL_FAIL;
    if (duckvep_coding_context_delta_fill(&ctx, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_delta_is_inframe_deletion_at(
            &delta, (int32_t)(((min_cds - 1u) / 3u) + 1u))) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_context_inframe_deletion_cov.fwd++;
    else g_context_inframe_deletion_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_inframe_deletion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delta in-frame deletion == edit-origin oracle";
    cfg.prop1 = prop_context_inframe_deletion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_deletion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_inframe_deletion_cov, 0, sizeof g_context_inframe_deletion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_inframe_deletion_cov.fwd > 0u);
    ASSERT(g_context_inframe_deletion_cov.rev > 0u);
    fprintf(stderr, "[context-inframe-deletion coverage] forward=%u reverse=%u\n",
            g_context_inframe_deletion_cov.fwd, g_context_inframe_deletion_cov.rev);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_context_inframe_insertion_cov;

static enum theft_trial_res prop_context_inframe_insertion_matches_oracle(struct theft *t,
                                                                         void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    uint32_t anchor_cds;
    uint32_t before_cds;
    (void)t;

    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) != 0u || before_cds <= 3u || before_cds >= s->cds_lenv - 3u) {
        return THEFT_TRIAL_FAIL;
    }
    if (!ctx.has_single_edit || ctx.single_edit_ref_len != 0u ||
        ctx.single_edit_alt_len != 3u || ctx.applied_edits != 1u ||
        ctx.single_edit_cds_start != before_cds + 1u) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_coding_context_delta_fill(&ctx, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_delta_is_inframe_insertion_at(
            &delta, (int32_t)(before_cds / 3u + 1u))) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_context_inframe_insertion_cov.fwd++;
    else g_context_inframe_insertion_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_inframe_insertion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delta in-frame insertion == edit-origin oracle";
    cfg.prop1 = prop_context_inframe_insertion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_inframe_insertion_cov, 0, sizeof g_context_inframe_insertion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_inframe_insertion_cov.fwd > 0u);
    ASSERT(g_context_inframe_insertion_cov.rev > 0u);
    fprintf(stderr, "[context-inframe-insertion coverage] forward=%u reverse=%u\n",
            g_context_inframe_insertion_cov.fwd, g_context_inframe_insertion_cov.rev);
    PASS();
}

static struct {
    uint32_t fwd;
    uint32_t rev;
    uint32_t lengthen;
    uint32_t shorten;
    uint32_t inframe;
    uint32_t protein_altering;
} g_context_delins_shape_cov;

static enum theft_trial_res prop_context_delins_shape_matches_oracle(struct theft *t,
                                                                    void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[96];
    uint8_t ref_pep[40];
    uint8_t alt_pep[40];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    int32_t protein_pos;
    (void)t;

    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!ctx.has_single_edit || ctx.applied_edits != 1u ||
        (ctx.single_edit_ref_len == 0u && ctx.single_edit_alt_len == 0u) ||
        ctx.length_diff == 0 || (ctx.length_diff % 3) != 0 ||
        ctx.single_edit_cds_start <= 3u) {
        return THEFT_TRIAL_FAIL;
    }
    if ((uint64_t)ctx.single_edit_cds_start +
            (uint64_t)ctx.single_edit_ref_len - 1u >
        (uint64_t)ctx.ref_cds_len - 3u) {
        return THEFT_TRIAL_FAIL;
    }
    protein_pos = (int32_t)(((ctx.single_edit_cds_start - 1u) / 3u) + 1u);
    if (duckvep_coding_context_delta_fill(&ctx, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_protein_altering) {
        if (!kprop_delta_is_protein_altering_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
        g_context_delins_shape_cov.protein_altering++;
    } else if (ctx.length_diff > 0) {
        if (!kprop_delta_is_inframe_insertion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
        g_context_delins_shape_cov.inframe++;
    } else {
        if (!kprop_delta_is_inframe_deletion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
        g_context_delins_shape_cov.inframe++;
    }
    if (ctx.length_diff > 0) g_context_delins_shape_cov.lengthen++;
    else g_context_delins_shape_cov.shorten++;
    if (s->strand > 0) g_context_delins_shape_cov.fwd++;
    else g_context_delins_shape_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_delins_shape_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delins shape == local-edge oracle";
    cfg.prop1 = prop_context_delins_shape_matches_oracle;
    cfg.type_info[0] = &kprop_delins_shape_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_delins_shape_cov, 0, sizeof g_context_delins_shape_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_delins_shape_cov.fwd > 0u);
    ASSERT(g_context_delins_shape_cov.rev > 0u);
    ASSERT(g_context_delins_shape_cov.lengthen > 0u);
    ASSERT(g_context_delins_shape_cov.shorten > 0u);
    ASSERT(g_context_delins_shape_cov.inframe > 0u);
    ASSERT(g_context_delins_shape_cov.protein_altering > 0u);
    fprintf(stderr,
            "[context-delins-shape coverage] forward=%u reverse=%u lengthen=%u "
            "shorten=%u inframe=%u protein_altering=%u\n",
            g_context_delins_shape_cov.fwd, g_context_delins_shape_cov.rev,
            g_context_delins_shape_cov.lengthen, g_context_delins_shape_cov.shorten,
            g_context_delins_shape_cov.inframe,
            g_context_delins_shape_cov.protein_altering);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; uint32_t lengthen; uint32_t shorten; }
    g_delta_scratch_indel_cov;

static enum theft_trial_res prop_delta_scratch_indel_matches_oracle(struct theft *t,
                                                                   void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[96];
    uint8_t ref_pep[40];
    uint8_t alt_pep[40];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t delta;
    duckvep_coding_context_t ctx;
    int32_t protein_pos;
    (void)t;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    protein_pos = (int32_t)(((ctx.single_edit_cds_start - 1u) / 3u) + 1u);
    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s->tx, &s->ex,
                                             &s->seq, &s->v, 0u, 0u, s->vpos,
                                             s->strand, &scratch, &delta);
    if (s->expect_protein_altering) {
        if (!kprop_delta_is_protein_altering_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
    } else if (ctx.length_diff > 0) {
        if (!kprop_delta_is_inframe_insertion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
    } else if (ctx.length_diff < 0) {
        if (!kprop_delta_is_inframe_deletion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
    } else return THEFT_TRIAL_FAIL;
    if (ctx.length_diff > 0) g_delta_scratch_indel_cov.lengthen++;
    else g_delta_scratch_indel_cov.shorten++;
    if (s->strand > 0) g_delta_scratch_indel_cov.fwd++;
    else g_delta_scratch_indel_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_with_scratch_indel_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta scratch INDEL == local delins-shape oracle";
    cfg.prop1 = prop_delta_scratch_indel_matches_oracle;
    cfg.type_info[0] = &kprop_delins_shape_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_scratch_indel_cov, 0, sizeof g_delta_scratch_indel_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_scratch_indel_cov.fwd > 0u);
    ASSERT(g_delta_scratch_indel_cov.rev > 0u);
    ASSERT(g_delta_scratch_indel_cov.lengthen > 0u);
    ASSERT(g_delta_scratch_indel_cov.shorten > 0u);
    fprintf(stderr,
            "[delta-scratch-indel coverage] forward=%u reverse=%u lengthen=%u shorten=%u\n",
            g_delta_scratch_indel_cov.fwd, g_delta_scratch_indel_cov.rev,
            g_delta_scratch_indel_cov.lengthen, g_delta_scratch_indel_cov.shorten);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_inframe_insertion_cov;

static enum theft_trial_res prop_annotate_inframe_insertion_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    uint32_t anchor_cds;
    uint32_t before_cds;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    duckvep_workspace_delta_route_stats_reset(ws);

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) != 0u || before_cds <= 3u || before_cds >= s->cds_lenv - 3u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    if (s->strand > 0) g_inframe_insertion_cov.fwd++; else g_inframe_insertion_cov.rev++;

    if (rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)((before_cds / 3u) + 1u)) { tr = THEFT_TRIAL_FAIL; goto done; }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL || stats->simple_indel != 1u ||
        stats->substitution_context != 0u ||
        stats->del_context != 0u ||
        stats->ins_context != 0u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_inframe_insertion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile codon-boundary in-frame insertion == CDS-position oracle";
    cfg.prop1 = prop_annotate_inframe_insertion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_inframe_insertion_cov, 0, sizeof g_inframe_insertion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_inframe_insertion_cov.fwd > 0u);
    ASSERT(g_inframe_insertion_cov.rev > 0u);
    fprintf(stderr, "[inframe_insertion coverage] forward=%u reverse=%u\n",
            g_inframe_insertion_cov.fwd, g_inframe_insertion_cov.rev);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; uint32_t full; } g_cursor_ins_route_cov;

static enum theft_trial_res prop_cursor_ins_route_matches_tile(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    uint32_t anchor_cds;
    uint32_t before_cds;
    int saw_full = 0;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (!saw_full || tile_n != cursor_n || tile_n != 1u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (tile_stats.simple_indel != 1u ||
        tile_stats.ins_context != 0u ||
        tile_stats.substitution_context != 0u ||
        tile_stats.del_context != 0u ||
        tile_stats.simple_indel != cursor_stats.simple_indel ||
        tile_stats.ins_context != cursor_stats.ins_context) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (!consequence_rows_equal(&tile_rows[0], &cursor_rows[0])) { res = THEFT_TRIAL_FAIL; goto done; }
    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) != 0u ||
        tile_rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) ||
        tile_rows[0].protein_pos != (int32_t)((before_cds / 3u) + 1u)) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (s->strand > 0) g_cursor_ins_route_cov.fwd++;
    else g_cursor_ins_route_cov.rev++;
    if (saw_full) g_cursor_ins_route_cov.full++;

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_ins_route_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor INS route == tile under output splits";
    cfg.prop1 = prop_cursor_ins_route_matches_tile;
    cfg.type_info[0] = &kprop_inframe_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_ins_route_cov, 0, sizeof g_cursor_ins_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_ins_route_cov.fwd > 0u);
    ASSERT(g_cursor_ins_route_cov.rev > 0u);
    ASSERT(g_cursor_ins_route_cov.full > 0u);
    fprintf(stderr, "[cursor-ins-route coverage] forward=%u reverse=%u full=%u\n",
            g_cursor_ins_route_cov.fwd, g_cursor_ins_route_cov.rev,
            g_cursor_ins_route_cov.full);
    PASS();
}

/* Independent VEP inframe_insertion test: is the ref peptide window a prefix OR suffix of the
 * alt window (an empty ref window trivially yes)? Mirrors VariationEffect::inframe_insertion,
 * re-derived here so the test does not restate the kernel helper. */
static int kprop_pep_window_prefix_or_suffix(
    const uint8_t *outer, uint32_t outer_off, uint32_t outer_len,
    const uint8_t *inner, uint32_t inner_off, uint32_t inner_len) {
    uint32_t i;
    int ok;
    if (inner_len == 0u) return 1;
    if (inner_len > outer_len) return 0;
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

static struct { uint32_t fwd; uint32_t rev; uint32_t inframe; uint32_t altering; } g_protein_altering_cov;

static enum theft_trial_res prop_annotate_protein_altering_insertion_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_haplotype_edit_t edit_scratch[4];
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint32_t anchor_cds;
    uint32_t before_cds;
    uint32_t rf, rl, af, al;
    int inframe;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) == 0u || before_cds <= 3u || before_cds >= s->cds_lenv - 3u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

    /* Independent expected term: build the coding context (the alt CDS construction is itself
     * validated by the cds-edit-builder oracle), diff the peptides, and apply VEP's prefix/
     * suffix rule. A non-boundary insertion is inframe_insertion when the flanking residues are
     * preserved (empty ref window, or ref window prefix/suffix of alt) and protein_altering when
     * the junction residue also changes. */
    if (duckvep_variant_physical_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand,
                                             edit_scratch, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    kprop_peptide_window_oracle(ctx.ref_peptide, ctx.ref_peptide_len,
                                ctx.alt_peptide, ctx.alt_peptide_len, &rf, &rl, &af, &al);
    if (af == 0u) { tr = THEFT_TRIAL_FAIL; goto done; }
    inframe = (rf == 0u) ||
              kprop_pep_window_prefix_or_suffix(ctx.alt_peptide, af - 1u, al - af + 1u,
                                                ctx.ref_peptide, rf - 1u, rl - rf + 1u);
    want = inframe ? DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)
                   : DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
    if (s->strand > 0) g_protein_altering_cov.fwd++; else g_protein_altering_cov.rev++;
    if (inframe) g_protein_altering_cov.inframe++; else g_protein_altering_cov.altering++;

    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)((before_cds / 3u) + 1u)) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_protein_altering_insertion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile non-boundary in-frame insertion == peptide-window oracle";
    cfg.prop1 = prop_annotate_protein_altering_insertion_matches_oracle;
    cfg.type_info[0] = &kprop_protein_altering_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_protein_altering_cov, 0, sizeof g_protein_altering_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_protein_altering_cov.fwd > 0u);
    ASSERT(g_protein_altering_cov.rev > 0u);
    ASSERT(g_protein_altering_cov.inframe > 0u);
    ASSERT(g_protein_altering_cov.altering > 0u);
    fprintf(stderr, "[non-boundary insertion coverage] forward=%u reverse=%u inframe_insertion=%u protein_altering=%u\n",
            g_protein_altering_cov.fwd, g_protein_altering_cov.rev,
            g_protein_altering_cov.inframe, g_protein_altering_cov.altering);
    PASS();
}

TEST annotate_start_lost_matches_oracle_for_any_start_codon_snv(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile start_lost SNV == start-codon oracle";
    cfg.prop1 = prop_annotate_start_lost_matches_oracle;
    cfg.type_info[0] = &kprop_start_codon_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_start_cov, 0, sizeof g_start_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_start_cov.sl > 0u);
    ASSERT(g_start_cov.sl_sg > 0u);
    ASSERT(g_start_cov.sl_syn > 0u);
    ASSERT(g_start_cov.syn > 0u);
    ASSERT(g_start_cov.retained > 0u);
    ASSERT(g_start_cov.lost_and_retained > 0u);
    fprintf(stderr,
            "[start-codon coverage] start_lost=%u co_stop_gained=%u "
            "co_synonymous=%u synonymous=%u start_retained=%u "
            "lost_and_retained=%u\n",
            g_start_cov.sl, g_start_cov.sl_sg, g_start_cov.sl_syn,
            g_start_cov.syn,
            g_start_cov.retained, g_start_cov.lost_and_retained);
    PASS();
}
