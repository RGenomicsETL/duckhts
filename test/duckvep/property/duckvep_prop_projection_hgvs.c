#include "duckvep_property.h"

/* ===================================================================== *
 * Coding coordinate projection (genomic -> cDNA -> CDS -> peptide).
 *
 * Oracle: a brute-force transcript-order base walk over small random exon
 * models. This deliberately does NOT use the model's cdna_start1/cdna_end1
 * arrays, while the implementation does, so it catches adapter/order/strand
 * mistakes instead of restating the implementation. The phase anchor covers the
 * positive Ensembl start-exon phase convention audited from VEP/Ensembl source:
 * a positive phase shifts CDS numbering by that many bases.
 * ===================================================================== */


static uint32_t proj_brute_cdna_len(const struct kprop_proj_scene *s) {
    uint32_t n = 0u;
    uint32_t i;
    for (i = 0u; i < (uint32_t)s->excnt; i++) n += (uint32_t)(s->ee[i] - s->es[i] + 1u);
    return n;
}

static int proj_brute_genomic_to_cdna(const struct kprop_proj_scene *s, uint32_t genomic,
                                      uint32_t *cdna, uint32_t *exon_idx) {
    uint32_t c = 1u;
    uint32_t i;
    int fwd = s->strand >= 0;
    for (i = 0u; i < (uint32_t)s->excnt; i++) {
        uint32_t len = (uint32_t)(s->ee[i] - s->es[i] + 1u);
        uint32_t j;
        for (j = 0u; j < len; j++, c++) {
            uint32_t g = fwd ? (uint32_t)(s->es[i] + j) : (uint32_t)(s->ee[i] - j);
            if (g == genomic) {
                *cdna = c;
                if (exon_idx != NULL) *exon_idx = i;
                return 1;
            }
        }
    }
    return 0;
}

int proj_brute_cdna_to_genomic(const struct kprop_proj_scene *s, uint32_t cdna,
                                      uint32_t *genomic, uint32_t *exon_idx) {
    uint32_t c = 1u;
    uint32_t i;
    int fwd = s->strand >= 0;
    for (i = 0u; i < (uint32_t)s->excnt; i++) {
        uint32_t len = (uint32_t)(s->ee[i] - s->es[i] + 1u);
        uint32_t j;
        for (j = 0u; j < len; j++, c++) {
            if (c == cdna) {
                *genomic = fwd ? (uint32_t)(s->es[i] + j) : (uint32_t)(s->ee[i] - j);
                if (exon_idx != NULL) *exon_idx = i;
                return 1;
            }
        }
    }
    return 0;
}

static int proj_brute_coding(const struct kprop_proj_scene *s, uint32_t genomic,
                             duckvep_coding_projection_t *out) {
    duckvep_coding_projection_t r;
    uint32_t cdna = 0u, exon_idx = 0u;
    uint32_t start_cdna = 0u, end_cdna = 0u, start_exon = 0u, dummy = 0u;
    uint32_t cds;
    uint32_t start_g;
    uint32_t end_g;
    uint8_t phase;
    memset(&r, 0, sizeof r);
    if (s->cds_s == 0u || s->cds_e == 0u) return 0;
    if (!proj_brute_genomic_to_cdna(s, genomic, &cdna, &exon_idx)) return 0;
    start_g = s->strand >= 0 ? s->cds_s : s->cds_e;
    end_g   = s->strand >= 0 ? s->cds_e : s->cds_s;
    if (!proj_brute_genomic_to_cdna(s, start_g, &start_cdna, &start_exon)) return 0;
    if (!proj_brute_genomic_to_cdna(s, end_g, &end_cdna, &dummy)) return 0;
    if (end_cdna < start_cdna || cdna < start_cdna || cdna > end_cdna) return 0;
    phase = s->phase[start_exon] > 0 ? (uint8_t)s->phase[start_exon] : 0u;
    cds = (uint32_t)(cdna - start_cdna + 1u + (uint32_t)phase);
    r.cdna_pos = cdna;
    r.cds_pos = cds;
    r.protein_pos = (uint32_t)((cds - 1u) / 3u + 1u);
    r.codon_offset = (uint8_t)((cds - 1u) % 3u);
    r.codon_start_cds = (uint32_t)(cds - (uint32_t)r.codon_offset);
    r.exon_idx = exon_idx;
    r.phase_offset = phase;
    *out = r;
    return 1;
}

void kprop_proj_scene_finish(struct kprop_proj_scene *s) {
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e;
    s->tx.transcript_count = 1u;
    s->ex.start1 = s->es; s->ex.end1 = s->ee; s->ex.cdna_start1 = s->cs; s->ex.cdna_end1 = s->ce;
    s->ex.phase = s->phase; s->ex.end_phase = s->end_phase; s->ex.exon_count = (size_t)s->excnt;
}

static enum theft_alloc_res kprop_proj_alloc(struct theft *t, void *env, void **instance) {
    struct kprop_proj_scene *s = (struct kprop_proj_scene *)calloc(1u, sizeof *s);
    uint32_t nex = (uint32_t)kprop_bounded(t, KPROP_MAX_PROJ_EXONS) + 1u;
    uint32_t asc_s[KPROP_MAX_PROJ_EXONS] = {0u};
    uint32_t asc_e[KPROP_MAX_PROJ_EXONS] = {0u};
    uint32_t base = (uint32_t)kprop_bounded(t, 10000u) + 100u;
    uint32_t cursor = base;
    uint32_t cdna = 1u;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;

    for (i = 0u; i < nex; i++) {
        uint32_t len = (uint32_t)kprop_bounded(t, 12u) + 3u;
        uint32_t gap = (uint32_t)kprop_bounded(t, 12u) + 1u;
        asc_s[i] = cursor;
        asc_e[i] = cursor + len - 1u;
        cursor = asc_e[i] + gap + 1u;
    }

    s->chrom = 0u;
    s->tstart = asc_s[0];
    s->tend = asc_e[nex - 1u];
    s->strand = (kprop_bounded(t, 2u) == 0u) ? (int8_t)1 : (int8_t)-1;
    s->exoff = 0u;
    s->excnt = (uint16_t)nex;
    for (i = 0u; i < nex; i++) {
        uint32_t src = s->strand >= 0 ? i : (uint32_t)(nex - 1u - i);
        uint32_t len = (uint32_t)(asc_e[src] - asc_s[src] + 1u);
        s->es[i] = asc_s[src];
        s->ee[i] = asc_e[src];
        s->cs[i] = cdna;
        s->ce[i] = cdna + len - 1u;
        s->phase[i] = (int8_t)-1;
        s->end_phase[i] = (int8_t)-1;
        cdna += len;
    }
    kprop_proj_scene_finish(s);

    if (kprop_bounded(t, 4u) != 0u) { /* mostly coding, but keep non-coding cases */
        uint32_t total = proj_brute_cdna_len(s);
        uint32_t c_start = (uint32_t)kprop_bounded(t, total) + 1u;
        uint32_t c_len = (uint32_t)kprop_bounded(t, (uint64_t)(total - c_start + 1u)) + 1u;
        uint32_t c_end = c_start + c_len - 1u;
        uint32_t g_start = 0u, g_end = 0u, start_exon = 0u;
        (void)proj_brute_cdna_to_genomic(s, c_start, &g_start, &start_exon);
        (void)proj_brute_cdna_to_genomic(s, c_end, &g_end, NULL);
        s->cds_s = g_start < g_end ? g_start : g_end;
        s->cds_e = g_start < g_end ? g_end : g_start;
        s->phase[start_exon] = (int8_t)kprop_bounded(t, 3u); /* Ensembl phase 0/1/2 */
    }

    *instance = s;
    return THEFT_ALLOC_OK;
}

static void kprop_proj_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_proj_info = {
    .alloc = kprop_proj_alloc,
    .free  = kprop_proj_free,
};

static enum theft_trial_res prop_projection_matches_bruteforce(struct theft *t, void *arg1) {
    const struct kprop_proj_scene *s = (const struct kprop_proj_scene *)arg1;
    uint32_t total = proj_brute_cdna_len(s);
    uint32_t cdna;
    uint32_t outside = s->tstart > 1u ? (uint32_t)(s->tstart - 1u) : (uint32_t)(s->tend + 1u);
    (void)t;

    for (cdna = 1u; cdna <= total; cdna++) {
        uint32_t want_g = 0u, want_exon = 0u;
        uint32_t got_g = 0u, got_exon = 0u;
        uint32_t got_cdna = 0u;
        duckvep_coding_projection_t want_cp, got_cp;
        int want_coding;
        int got_coding;
        if (!proj_brute_cdna_to_genomic(s, cdna, &want_g, &want_exon)) return THEFT_TRIAL_FAIL;
        if (!duckvep_project_cdna_to_genomic(&s->tx, &s->ex, 0u, cdna, &got_g, &got_exon)) return THEFT_TRIAL_FAIL;
        if (got_g != want_g || got_exon != want_exon) return THEFT_TRIAL_FAIL;
        if (!duckvep_project_genomic_to_cdna(&s->tx, &s->ex, 0u, want_g, &got_cdna, &got_exon)) return THEFT_TRIAL_FAIL;
        if (got_cdna != cdna || got_exon != want_exon) return THEFT_TRIAL_FAIL;

        want_coding = proj_brute_coding(s, want_g, &want_cp);
        got_coding = duckvep_project_coding_base(&s->tx, &s->ex, 0u, want_g, &got_cp);
        if (want_coding != got_coding) return THEFT_TRIAL_FAIL;
        if (want_coding) {
            if (got_cp.cdna_pos != want_cp.cdna_pos || got_cp.cds_pos != want_cp.cds_pos ||
                got_cp.protein_pos != want_cp.protein_pos ||
                got_cp.codon_start_cds != want_cp.codon_start_cds ||
                got_cp.exon_idx != want_cp.exon_idx ||
                got_cp.codon_offset != want_cp.codon_offset ||
                got_cp.phase_offset != want_cp.phase_offset) {
                return THEFT_TRIAL_FAIL;
            }
        }
    }

    if (duckvep_project_genomic_to_cdna(&s->tx, &s->ex, 0u, outside, &cdna, NULL)) return THEFT_TRIAL_FAIL;
    if (duckvep_project_coding_base(&s->tx, &s->ex, 0u, outside, &(duckvep_coding_projection_t){0})) return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

static int proj_brute_transcript_coordinate(
    const struct kprop_proj_scene      *s,
    uint32_t                            genomic,
    duckvep_transcript_coordinate_t    *out) {

    duckvep_transcript_coordinate_t r;
    uint32_t cdna = 0u;
    uint32_t exon_idx = 0u;
    uint32_t lower_idx = UINT32_MAX;
    uint32_t upper_idx = UINT32_MAX;
    uint32_t lower_end = 0u;
    uint32_t upper_start = UINT32_MAX;
    uint32_t i;

    memset(&r, 0, sizeof r);
    if (genomic < s->tstart || genomic > s->tend) return 0;
    r.genomic_pos1 = genomic;
    if (proj_brute_genomic_to_cdna(s, genomic, &cdna, &exon_idx)) {
        r.cdna_anchor1 = cdna;
        r.exon_idx = exon_idx;
        r.exonic = 1u;
        *out = r;
        return 1;
    }
    for (i = 0u; i < (uint32_t)s->excnt; i++) {
        if (s->ee[i] < genomic &&
            (lower_idx == UINT32_MAX || s->ee[i] > lower_end)) {
            lower_idx = i;
            lower_end = s->ee[i];
        }
        if (s->es[i] > genomic &&
            (upper_idx == UINT32_MAX || s->es[i] < upper_start)) {
            upper_idx = i;
            upper_start = s->es[i];
        }
    }
    if (lower_idx == UINT32_MAX || upper_idx == UINT32_MAX) return 0;
    {
        uint32_t lower_distance = genomic - lower_end;
        uint32_t upper_distance = upper_start - genomic;
        if (lower_distance < upper_distance ||
            (lower_distance == upper_distance && s->strand > 0)) {
            r.exon_idx = lower_idx;
            if (s->strand > 0) {
                r.cdna_anchor1 = s->ce[lower_idx];
                r.intron_offset = (int32_t)lower_distance;
            } else {
                r.cdna_anchor1 = s->cs[lower_idx];
                r.intron_offset = -(int32_t)lower_distance;
            }
        } else {
            r.exon_idx = upper_idx;
            if (s->strand > 0) {
                r.cdna_anchor1 = s->cs[upper_idx];
                r.intron_offset = -(int32_t)upper_distance;
            } else {
                r.cdna_anchor1 = s->ce[upper_idx];
                r.intron_offset = (int32_t)upper_distance;
            }
        }
    }
    *out = r;
    return 1;
}

static enum theft_trial_res prop_transcript_coordinate_matches_bruteforce(
    struct theft *t, void *arg1) {

    const struct kprop_proj_scene *s = (const struct kprop_proj_scene *)arg1;
    uint32_t genomic;
    (void)t;

    for (genomic = s->tstart; genomic <= s->tend; genomic++) {
        duckvep_transcript_coordinate_t want;
        duckvep_transcript_coordinate_t got;
        duckvep_transcript_edit_status_t status;
        if (!proj_brute_transcript_coordinate(s, genomic, &want)) {
            return THEFT_TRIAL_FAIL;
        }
        status = duckvep_project_transcript_coordinate(
            &s->tx, &s->ex, 0u, genomic, &got);
        if (status != DUCKVEP_TRANSCRIPT_EDIT_OK ||
            got.genomic_pos1 != want.genomic_pos1 ||
            got.cdna_anchor1 != want.cdna_anchor1 ||
            got.exon_idx != want.exon_idx ||
            got.intron_offset != want.intron_offset ||
            got.exonic != want.exonic) {
            return THEFT_TRIAL_FAIL;
        }
        if (genomic == UINT32_MAX) break;
    }
    return THEFT_TRIAL_PASS;
}

TEST transcript_coordinate_known_intronic_ties(void) {
    duckvep_transcript_coordinate_t coordinate;

    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 208u; s.strand = (int8_t)1;
        s.excnt = 2u;
        s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
        s.es[1] = 199u; s.ee[1] = 208u; s.cs[1] = 11u; s.ce[1] = 20u;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
                  duckvep_project_transcript_coordinate(
                      &s.tx, &s.ex, 0u, 110u, &coordinate));
        ASSERT_EQ(10u, coordinate.cdna_anchor1);
        ASSERT_EQ(1, coordinate.intron_offset);
        ASSERT_EQ(0u, (uint32_t)coordinate.exonic);
        /* 154 is equidistant from exon bases 109 and 199. VEP anchors the
         * tie to the lower-genomic exon on a forward transcript. */
        ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
                  duckvep_project_transcript_coordinate(
                      &s.tx, &s.ex, 0u, 154u, &coordinate));
        ASSERT_EQ(10u, coordinate.cdna_anchor1);
        ASSERT_EQ(45, coordinate.intron_offset);
    }
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 208u; s.strand = (int8_t)-1;
        s.excnt = 2u;
        s.es[0] = 199u; s.ee[0] = 208u; s.cs[0] = 1u; s.ce[0] = 10u;
        s.es[1] = 100u; s.ee[1] = 109u; s.cs[1] = 11u; s.ce[1] = 20u;
        kprop_proj_scene_finish(&s);
        /* The same genomic tie anchors to the higher-genomic exon on a
         * reverse transcript, which is +45 from its transcript-order end. */
        ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
                  duckvep_project_transcript_coordinate(
                      &s.tx, &s.ex, 0u, 154u, &coordinate));
        ASSERT_EQ(10u, coordinate.cdna_anchor1);
        ASSERT_EQ(45, coordinate.intron_offset);
    }
    PASS();
}

TEST transcript_edit_orders_endpoints_and_reuses_cds_edits(void) {
    static const uint8_t alleles[] = {'A', 'G'};
    static const uint16_t chrom[] = {0u};
    static const uint32_t pos[] = {203u};
    static const uint32_t end[] = {203u};
    static const uint32_t roff[] = {0u};
    static const uint16_t rlen[] = {1u};
    static const uint32_t aoff[] = {1u};
    static const uint16_t alen[] = {1u};
    static const uint8_t kind[] = {(uint8_t)DUCKVEP_KIND_SNV};
    static const uint8_t cds[] = {'T','T','T','T','T','C','C','C','C','C','C','C'};
    static const uint64_t cds_off[] = {0u};
    static const uint32_t cds_len[] = {12u};
    static const uint8_t table[] = {1u};
    struct kprop_proj_scene s;
    duckvep_variant_batch_t variants;
    duckvep_sequence_pool_t sequences;
    duckvep_haplotype_edit_t scratch[1];
    duckvep_event_t event;
    duckvep_transcript_edit_t edit;
    duckvep_transcript_edit_t no_scratch;
    duckvep_transcript_edit_t projected;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)-1;
    s.excnt = 2u; s.cds_s = 104u; s.cds_e = 205u;
    s.es[0] = 200u; s.ee[0] = 209u; s.cs[0] = 1u; s.ce[0] = 10u; s.phase[0] = 0;
    s.es[1] = 100u; s.ee[1] = 109u; s.cs[1] = 11u; s.ce[1] = 20u; s.phase[1] = 0;
    kprop_proj_scene_finish(&s);
    memset(&variants, 0, sizeof variants);
    variants.chrom_id = chrom; variants.pos1 = pos; variants.end1 = end;
    variants.ref_offset = roff; variants.ref_length = rlen;
    variants.alt_offset = aoff; variants.alt_length = alen;
    variants.allele_bytes = alleles; variants.allele_bytes_len = sizeof alleles;
    variants.variant_kind = kind; variants.count = 1u;
    duckvep_event_load(&variants, 0u, &event);
    memset(&sequences, 0, sizeof sequences);
    sequences.cds_bytes = cds; sequences.cds_bytes_len = sizeof cds;
    sequences.cds_offset = cds_off; sequences.cds_length = cds_len;
    sequences.codon_table = table; sequences.transcript_count = 1u;

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              duckvep_transcript_edit_build(
                  &s.tx, &s.ex, &sequences, &variants, 0u, 0u,
                  scratch, 1u, &edit));
    ASSERT_EQ(7u, edit.first.cdna_anchor1);
    ASSERT_EQ(7u, edit.last.cdna_anchor1);
    ASSERT(edit.ref == alleles);
    ASSERT(edit.alt == alleles + 1u);
    ASSERT_EQ(1u, (uint32_t)edit.ref_length);
    ASSERT_EQ(1u, (uint32_t)edit.alt_length);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, edit.cds_status);
    ASSERT_EQ(1u, (uint32_t)edit.cds_built);
    ASSERT_EQ(1u, edit.cds_edits.count);
    ASSERT_EQ(3u, edit.cds_edits.edits[0].cds_start);
    ASSERT_EQ(1u, edit.cds_edits.edits[0].ref_len);
    ASSERT_EQ(1u, edit.cds_edits.edits[0].alt_len);
    ASSERT(edit.cds_edits.edits[0].ref == alleles);
    ASSERT(edit.cds_edits.edits[0].alt == alleles + 1u);

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              duckvep_transcript_edit_build(
                  &s.tx, &s.ex, &sequences, &variants, 0u, 0u,
                  NULL, 0u, &no_scratch));
    ASSERT_EQ(edit.first.cdna_anchor1, no_scratch.first.cdna_anchor1);
    ASSERT_EQ(edit.last.cdna_anchor1, no_scratch.last.cdna_anchor1);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL, no_scratch.cds_status);
    ASSERT_EQ(1u, (uint32_t)no_scratch.cds_built);
    ASSERT(no_scratch.cds_edits.edits == NULL);
    ASSERT_EQ(0u, no_scratch.cds_edits.count);

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              duckvep_transcript_edit_project_prepared(
                  &s.tx, &s.ex, &variants, 0u, 0u, &event, &projected));
    ASSERT_EQ(0u, (uint32_t)projected.cds_built);
    ASSERT_EQ(edit.first.cdna_anchor1, projected.first.cdna_anchor1);
    ASSERT_EQ(edit.last.cdna_anchor1, projected.last.cdna_anchor1);
    ASSERT(projected.ref == edit.ref);
    ASSERT(projected.alt == edit.alt);
    ASSERT_EQ(edit.ref_length, projected.ref_length);
    ASSERT_EQ(edit.alt_length, projected.alt_length);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_transcript_edit_cds_fill_prepared(
                  &s.tx, &s.ex, &sequences, &variants,
                  scratch, 1u, &projected));
    ASSERT_EQ(1u, (uint32_t)projected.cds_built);
    ASSERT_EQ(edit.cds_status, projected.cds_status);
    ASSERT_EQ(edit.cds_edits.count, projected.cds_edits.count);
    ASSERT_EQ(edit.cds_edits.edits[0].cds_start,
              projected.cds_edits.edits[0].cds_start);
    ASSERT_EQ(edit.cds_edits.edits[0].ref_len,
              projected.cds_edits.edits[0].ref_len);
    ASSERT_EQ(edit.cds_edits.edits[0].alt_len,
              projected.cds_edits.edits[0].alt_len);
    ASSERT(edit.cds_edits.edits[0].ref == projected.cds_edits.edits[0].ref);
    ASSERT(edit.cds_edits.edits[0].alt == projected.cds_edits.edits[0].alt);
    PASS();
}

TEST transcript_edit_rejects_unrepresentable_transcript_index(void) {
#if SIZE_MAX > UINT32_MAX
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_transcript_edit_t edit;
    static const uint8_t kind[] = {(uint8_t)DUCKVEP_KIND_SNV};

    memset(&transcripts, 0, sizeof transcripts);
    memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants);
    transcripts.transcript_count = (size_t)UINT32_MAX + 2u;
    variants.variant_kind = kind;
    variants.count = 1u;
    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG,
              duckvep_transcript_edit_build(
                  &transcripts, &exons, NULL, &variants, 0u,
                  (size_t)UINT32_MAX + 1u, NULL, 0u, &edit));
#endif
    PASS();
}

TEST transcript_edit_rejects_missing_transcript_spans(void) {
    static const uint8_t alleles[] = {'A', 'G'};
    static const uint16_t chrom[] = {0u};
    static const uint32_t pos[] = {105u};
    static const uint32_t end[] = {105u};
    static const uint32_t roff[] = {0u};
    static const uint16_t rlen[] = {1u};
    static const uint32_t aoff[] = {1u};
    static const uint16_t alen[] = {1u};
    static const uint8_t kind[] = {(uint8_t)DUCKVEP_KIND_SNV};
    struct kprop_proj_scene s;
    duckvep_variant_batch_t variants;
    duckvep_event_t event;
    duckvep_transcript_edit_t edit;
    const uint32_t *saved_start1;
    const uint32_t *saved_end1;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 109u; s.strand = (int8_t)1;
    s.excnt = 1u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.phase[0] = -1;
    kprop_proj_scene_finish(&s);
    memset(&variants, 0, sizeof variants);
    variants.chrom_id = chrom; variants.pos1 = pos; variants.end1 = end;
    variants.ref_offset = roff; variants.ref_length = rlen;
    variants.alt_offset = aoff; variants.alt_length = alen;
    variants.allele_bytes = alleles; variants.allele_bytes_len = sizeof alleles;
    variants.variant_kind = kind; variants.count = 1u;
    duckvep_event_load(&variants, 0u, &event);

    saved_start1 = s.tx.start1;
    saved_end1 = s.tx.end1;
    s.tx.start1 = NULL;
    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG,
              duckvep_transcript_edit_build_prepared(
                  &s.tx, &s.ex, NULL, &variants, 0u, 0u, &event,
                  NULL, 0u, &edit));
    s.tx.start1 = saved_start1;
    s.tx.end1 = NULL;
    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG,
              duckvep_transcript_edit_build_prepared(
                  &s.tx, &s.ex, NULL, &variants, 0u, 0u, &event,
                  NULL, 0u, &edit));
    s.tx.end1 = saved_end1;
    PASS();
}

TEST transcript_edit_prepared_event_is_the_geometry_authority(void) {
    static const uint8_t alleles[] = {'A', 'G'};
    static const uint16_t chrom[] = {0u};
    static const uint32_t pos[] = {105u};
    static const uint32_t end[] = {105u};
    static const uint32_t roff[] = {0u};
    static const uint16_t rlen[] = {1u};
    static const uint32_t aoff[] = {1u};
    static const uint16_t alen[] = {1u};
    static const uint8_t kind[] = {(uint8_t)DUCKVEP_KIND_SNV};
    struct kprop_proj_scene s;
    duckvep_variant_batch_t variants;
    duckvep_event_t prepared;
    duckvep_transcript_edit_t raw_edit;
    duckvep_transcript_edit_t prepared_edit;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 109u; s.strand = (int8_t)1;
    s.excnt = 1u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.phase[0] = -1;
    kprop_proj_scene_finish(&s);
    memset(&variants, 0, sizeof variants);
    variants.chrom_id = chrom; variants.pos1 = pos; variants.end1 = end;
    variants.ref_offset = roff; variants.ref_length = rlen;
    variants.alt_offset = aoff; variants.alt_length = alen;
    variants.allele_bytes = alleles; variants.allele_bytes_len = sizeof alleles;
    variants.variant_kind = kind; variants.count = 1u;

    duckvep_event_load(&variants, 0u, &prepared);
    prepared.start1 = 106u;
    prepared.end1 = 106u;
    prepared.feature_start1 = 106u;
    prepared.feature_end1 = 106u;
    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              duckvep_transcript_edit_build(
                  &s.tx, &s.ex, NULL, &variants, 0u, 0u,
                  NULL, 0u, &raw_edit));
    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              duckvep_transcript_edit_build_prepared(
                  &s.tx, &s.ex, NULL, &variants, 0u, 0u, &prepared,
                  NULL, 0u, &prepared_edit));
    ASSERT_EQ(6u, raw_edit.first.cdna_anchor1);
    ASSERT_EQ(7u, prepared_edit.first.cdna_anchor1);
    ASSERT_EQ(106u, prepared_edit.event.start1);
    PASS();
}

TEST transcript_edit_retains_raw_feature_and_semantic_alleles(void) {
    static const uint8_t alleles[] = {'G','A','C', 'G','G','T'};
    static const uint16_t chrom[] = {0u};
    static const uint32_t pos[] = {100u};
    static const uint32_t end[] = {102u};
    static const uint32_t roff[] = {0u};
    static const uint16_t rlen[] = {3u};
    static const uint32_t aoff[] = {3u};
    static const uint16_t alen[] = {3u};
    static const uint8_t kind[] = {(uint8_t)DUCKVEP_KIND_MNV};
    struct kprop_proj_scene s;
    duckvep_variant_batch_t variants;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    char rendered[64];
    size_t required;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 109u; s.strand = (int8_t)1;
    s.excnt = 1u; s.cds_s = 100u; s.cds_e = 109u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);
    memset(&variants, 0, sizeof variants);
    variants.chrom_id = chrom; variants.pos1 = pos; variants.end1 = end;
    variants.ref_offset = roff; variants.ref_length = rlen;
    variants.alt_offset = aoff; variants.alt_length = alen;
    variants.allele_bytes = alleles; variants.allele_bytes_len = sizeof alleles;
    variants.variant_kind = kind; variants.count = 1u;

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              duckvep_transcript_edit_build(
                  &s.tx, &s.ex, NULL, &variants, 0u, 0u,
                  NULL, 0u, &edit));
    ASSERT(edit.raw_ref == alleles);
    ASSERT(edit.raw_alt == alleles + 3u);
    ASSERT_EQ(3u, (uint32_t)edit.raw_ref_length);
    ASSERT_EQ(3u, (uint32_t)edit.raw_alt_length);
    ASSERT(edit.feature_ref == edit.raw_ref);
    ASSERT(edit.feature_alt == edit.raw_alt);
    ASSERT_EQ(3u, (uint32_t)edit.feature_ref_length);
    ASSERT_EQ(3u, (uint32_t)edit.feature_alt_length);
    ASSERT_EQ(1u, edit.feature_first.cdna_anchor1);
    ASSERT_EQ(3u, edit.feature_last.cdna_anchor1);
    ASSERT(edit.ref == alleles + 1u);
    ASSERT(edit.alt == alleles + 4u);
    ASSERT_EQ(2u, (uint32_t)edit.ref_length);
    ASSERT_EQ(2u, (uint32_t)edit.alt_length);
    ASSERT_EQ(2u, edit.first.cdna_anchor1);
    ASSERT_EQ(3u, edit.last.cdna_anchor1);

    /* AC>GT is an inversion after clipping, but VEP classifies the complete
     * equal-length feature GAC>GGT first; that feature is not an inversion. */
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted(
                  &s.tx, &s.ex, NULL, &edit, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_DNA_REPLACEMENT, fact.shape);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.2_3delinsGT", rendered));
    PASS();
}

TEST hgvs_transcript_coordinate_numbering_known(void) {
    struct kprop_proj_scene s;
    duckvep_transcript_coordinate_t coordinate;
    duckvep_hgvs_coordinate_t hgvs;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)1;
    s.excnt = 2u; s.cds_s = 104u; s.cds_e = 205u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.es[1] = 200u; s.ee[1] = 209u; s.cs[1] = 11u; s.ce[1] = 20u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);
    memset(&coordinate, 0, sizeof coordinate);
    coordinate.exonic = 1u;

    coordinate.cdna_anchor1 = 1u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_C, hgvs.kind);
    ASSERT_EQ(-4, hgvs.base);
    coordinate.cdna_anchor1 = 4u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(-1, hgvs.base);
    coordinate.cdna_anchor1 = 5u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(1, hgvs.base);
    coordinate.cdna_anchor1 = 16u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_C, hgvs.kind);
    ASSERT_EQ(12, hgvs.base);
    coordinate.cdna_anchor1 = 17u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_C_STAR, hgvs.kind);
    ASSERT_EQ(1, hgvs.base);

    coordinate.exonic = 0u;
    coordinate.cdna_anchor1 = 10u;
    coordinate.intron_offset = 45;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_C, hgvs.kind);
    ASSERT_EQ(6, hgvs.base);
    ASSERT_EQ(45, hgvs.intron_offset);
    coordinate.cdna_anchor1 = 16u;
    coordinate.intron_offset = 1;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_C_STAR, hgvs.kind);
    ASSERT_EQ(0, hgvs.base);
    ASSERT_EQ(1, hgvs.intron_offset);

    /* VEP's generic _get_cDNA_position path does not add the positive
     * first-CDS-exon phase. hgvs_transcript has a separate phase-aware fast
     * path for literal exonic SNPs, tested below. */
    s.phase[0] = 2;
    s.flags = (uint64_t)DUCKVEP_TX_CDS_START_NF;
    coordinate.exonic = 1u;
    coordinate.intron_offset = 0;
    coordinate.cdna_anchor1 = 1u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_C, hgvs.kind);
    ASSERT_EQ(-4, hgvs.base);
    coordinate.cdna_anchor1 = 5u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(1, hgvs.base);

    s.cds_s = 0u;
    s.cds_e = 0u;
    coordinate.exonic = 0u;
    coordinate.cdna_anchor1 = 10u;
    coordinate.intron_offset = 45;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_coordinate_from_transcript(
        &s.tx, &s.ex, 0u, &coordinate, &hgvs));
    ASSERT_EQ(DUCKVEP_HGVS_COORDINATE_N, hgvs.kind);
    ASSERT_EQ(10, hgvs.base);
    ASSERT_EQ(45, hgvs.intron_offset);
    PASS();
}

TEST hgvs_exonic_snp_phase_fast_path_is_representation_specific(void) {
    static const uint8_t ref_a[] = {'A'};
    static const uint8_t alt_g[] = {'G'};
    static const uint8_t feature_ref[] = {'C', 'A'};
    static const uint8_t feature_alt[] = {'C', 'G'};
    struct kprop_proj_scene s;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    char rendered[64];
    size_t required;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 109u; s.strand = (int8_t)1;
    s.excnt = 1u; s.cds_s = 100u; s.cds_e = 109u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.phase[0] = 2;
    s.flags = (uint64_t)DUCKVEP_TX_CDS_START_NF;
    kprop_proj_scene_finish(&s);

    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u;
    edit.transcript_strand = (int8_t)1;
    edit.first.genomic_pos1 = 100u;
    edit.first.cdna_anchor1 = 1u;
    edit.first.exonic = 1u;
    edit.last = edit.first;
    edit.feature_first = edit.first;
    edit.feature_last = edit.first;
    edit.ref = ref_a; edit.ref_length = 1u;
    edit.alt = alt_g; edit.alt_length = 1u;
    edit.feature_ref = ref_a; edit.feature_ref_length = 1u;
    edit.feature_alt = alt_g; edit.feature_alt_length = 1u;

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build(&s.tx, &s.ex, &edit, &fact));
    ASSERT_EQ(3, fact.first.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.3A>G", rendered));

    /* The same one-base semantic edit inside a two-base uploaded feature is
     * not VEP var_class SNP. It follows _get_cDNA_position and stays c.2. */
    edit.first.genomic_pos1 = 101u;
    edit.first.cdna_anchor1 = 2u;
    edit.last = edit.first;
    edit.feature_first.genomic_pos1 = 100u;
    edit.feature_first.cdna_anchor1 = 1u;
    edit.feature_last.genomic_pos1 = 101u;
    edit.feature_last.cdna_anchor1 = 2u;
    edit.feature_ref = feature_ref; edit.feature_ref_length = 2u;
    edit.feature_alt = feature_alt; edit.feature_alt_length = 2u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build(&s.tx, &s.ex, &edit, &fact));
    ASSERT_EQ(2, fact.first.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.2A>G", rendered));
    PASS();
}

TEST hgvs_protein_mapper_endpoints_define_applicability(void) {
    static const uint8_t ref[] = {'A', 'A', 'A'};
    static const uint8_t alt[] = {'C', 'C', 'C'};
    struct kprop_proj_scene s;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    int defined;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)1;
    s.excnt = 2u; s.cds_s = 100u; s.cds_e = 209u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.es[1] = 200u; s.ee[1] = 209u; s.cs[1] = 11u; s.ce[1] = 20u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);

    memset(&edit, 0, sizeof edit);
    memset(&fact, 0, sizeof fact);
    edit.tx_idx = 0u;
    edit.transcript_strand = (int8_t)1;
    edit.event.kind = (uint8_t)DUCKVEP_KIND_MNV;
    edit.event.feature_start1 = 108u;
    edit.event.feature_end1 = 110u;
    fact.ref = ref; fact.ref_length = 3u;
    fact.alt = alt; fact.alt_length = 3u;
    fact.transcript_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_coordinates_defined(
                  &s.tx, &s.ex, &edit, &fact, &defined));
    ASSERT_EQ(0, defined);

    edit.event.feature_end1 = 109u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_coordinates_defined(
                  &s.tx, &s.ex, &edit, &fact, &defined));
    ASSERT_EQ(1, defined);

    /* The same check consumes the shifted semantic placement for an indel. */
    edit.event.kind = (uint8_t)DUCKVEP_KIND_DEL;
    edit.event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_NONE;
    fact.shape = (uint8_t)DUCKVEP_HGVS_DNA_DELETION;
    fact.ref_length = 2u;
    fact.alt_length = 0u;
    fact.placed_start1 = 108u;
    fact.placed_end1 = 109u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_coordinates_defined(
                  &s.tx, &s.ex, &edit, &fact, &defined));
    ASSERT_EQ(1, defined);
    fact.ref_length = 3u;
    fact.placed_end1 = 110u;
    edit.event.feature_end1 = 110u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_coordinates_defined(
                  &s.tx, &s.ex, &edit, &fact, &defined));
    ASSERT_EQ(0, defined);

    /* HGVSc clamps a partially overlapping VariationFeature, but VEP's
     * genomic2pep() applicability check still sees the complete feature. */
    edit.event.feature_start1 = 108u;
    edit.event.feature_end1 = 210u;
    fact.ref_length = 2u;
    fact.placed_start1 = 108u;
    fact.placed_end1 = 109u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_coordinates_defined(
                  &s.tx, &s.ex, &edit, &fact, &defined));
    ASSERT_EQ(0, defined);
    PASS();
}

TEST hgvs_protein_pair_reuses_fused_facts_and_bounds_shift_scratch(void) {
    static const uint8_t genome[] =
        "AAAAAAAAAAATGGGTCCTTAAAAAGAACAATAATAACTAGCTGAAAAAAAAAAA";
    static const uint32_t positions[] = {18u, 19u, 20u, 21u, 20u, 17u};
    static const char *alleles[] = {"CCA", "TTT", "TTT", "ATA", "TTG", "CCC"};
    static const char *expected[] = {NULL, "p.Ter4LeufsTer9", "p.Ter4delinsLeuTer",
        "p.Ter4delinsLeuTer", "p.Ter4=", "p.Ter4LeufsTer9"};
    struct kprop_proj_scene s = {0};
    s.chrom = 0u; s.tstart = 11u; s.tend = 45u; s.strand = 1;
    s.excnt = 1u; s.cds_s = 11u; s.cds_e = 22u;
    s.es[0] = 11u; s.ee[0] = 45u; s.cs[0] = 1u; s.ce[0] = 35u;
    s.flags = DUCKVEP_TX_HAS_TRANSLATION | DUCKVEP_TX_BIOTYPE_PROTEIN_CODING;
    kprop_proj_scene_finish(&s);
    uint64_t offset = 0u, post_offset = 12u;
    uint32_t length = 12u, empty = 0u, post_length = 23u;
    uint8_t table = DUCKVEP_CODON_TABLE_STANDARD;
    duckvep_sequence_pool_t sequences = {0};
    sequences.cds_bytes = genome + 10u; sequences.cds_bytes_len = length;
    sequences.cds_offset = &offset; sequences.cds_length = &length;
    sequences.codon_table = &table; sequences.transcript_count = 1u;
    sequences.flank_bytes = genome + 10u; sequences.flank_bytes_len = 35u;
    sequences.pre_cds_offset = &offset; sequences.pre_cds_length = &empty;
    sequences.post_cds_offset = &post_offset; sequences.post_cds_length = &post_length;
    sequences.flanks_complete = 1u;
    duckvep_hgvs_reference_window_t reference = {genome, sizeof genome - 1u, 1u, 0u};
    uint32_t ro = 0u, ao = 1u;
    uint16_t rl = 1u, al = 2u;
    uint8_t kind = DUCKVEP_KIND_INS;
    duckvep_variant_batch_t variants = {0};
    variants.chrom_id = &s.chrom; variants.ref_offset = &ro; variants.alt_offset = &ao;
    variants.ref_length = &rl; variants.alt_length = &al; variants.variant_kind = &kind;
    variants.count = 1u; variants.allele_bytes_len = 3u;
    for (size_t i = 0u; i < sizeof positions / sizeof positions[0]; i++) {
        variants.pos1 = variants.end1 = positions + i;
        variants.allele_bytes = (const uint8_t *)alleles[i];
        duckvep_event_t event;
        ASSERT(duckvep_event_prepare_small(positions[i], variants.allele_bytes, rl,
            variants.allele_bytes + ao, al, &event));
        event.chrom_id = 0u;
        duckvep_haplotype_edit_t edits[4];
        duckvep_transcript_edit_t edit;
        ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK, duckvep_transcript_edit_build_prepared(
            &s.tx, &s.ex, &sequences, &variants, 0u, 0u, &event, edits, 4u, &edit));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, edit.cds_status);
        uint8_t cds[64], rp[32], ap[32];
        duckvep_delta_scratch_t scratch = {edits, 4u, cds, sizeof cds,
            rp, sizeof rp, ap, sizeof ap};
        duckvep_coding_context_t context;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK, duckvep_model_coding_context_build(
            &s.tx, &s.ex, &sequences, 0u, 1, &event, &edit.cds_edits,
            cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &context));
        duckvep_sequence_delta_t delta;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK, duckvep_coding_context_delta_fill(
            &context, s.flags, &delta));
        duckvep_pair_facts_t facts = {0};
        facts.event = &event; facts.transcript_edit = &edit;
        facts.transcript_edit_status = DUCKVEP_TRANSCRIPT_EDIT_OK;
        facts.coding_context = &context; facts.coding_context_valid = 1u;
        facts.delta = &delta; facts.projection_exon_hint = 0u;
        duckvep_consequence_t row = {0};
        row.overlap_object_kind = DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
        row.region_mask = DUCKVEP_REGION_CDS;
        row.flags = duckvep_sequence_delta_consequence_flags(&delta, 1);
        duckvep_hgvs_dna_fact_t dna;
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
            &s.tx, &s.ex, &reference, &reference, &edit, &dna));
        duckvep_transcript_edit_t placed_edit;
        duckvep_hgvs_dna_fact_t placed_dna;
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_dna_pair_build(&s.tx, &s.ex, &sequences,
            &variants, &facts, &reference, &reference, &scratch, &placed_edit, &placed_dna));
        ASSERT_MEM_EQ(&dna, &placed_dna, sizeof dna);
        ASSERT_EQ(DUCKVEP_HGVS_MISSING_REFERENCE, duckvep_hgvs_dna_pair_build(&s.tx, &s.ex,
            &sequences, &variants, &facts, NULL, NULL, &scratch, &placed_edit, &placed_dna));
        duckvep_coding_context_t saved_context = context;
        duckvep_sequence_delta_t saved_delta = delta;
        duckvep_transcript_edit_t saved_edit = edit;
        duckvep_hgvs_protein_pair_t pair;
        duckvep_hgvs_protein_fact_t zero = {0};
        size_t required;
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_pair_build(
            &s.tx, &s.ex, &sequences, &variants, &row, &facts, &dna, &reference,
            &scratch, NULL, 0u, &required, &pair);
        if (!expected[i]) {
            ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE, status);
            ASSERT_MEM_EQ(&zero, &pair.fact, sizeof zero);
            ASSERT_EQ(0u, required);
        } else {
            uint8_t rotated[4] = {0xa5u, 0xa5u, 0xa5u, 0xa5u};
            if (dna.shift_offset) {
                ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, status);
                ASSERT_EQ(2u, required);
                ASSERT_MEM_EQ(&zero, &pair.fact, sizeof zero);
                ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_pair_build(
                    &s.tx, &s.ex, &sequences, &variants, &row, &facts, &dna, &reference,
                    &scratch, rotated + 1u, 1u, &required, &pair));
                ASSERT_EQ(0xa5u, rotated[1]);
                status = duckvep_hgvs_protein_pair_build(&s.tx, &s.ex, &sequences, &variants,
                    &row, &facts, &dna, &reference, &scratch, rotated + 1u, 2u, &required, &pair);
                ASSERT(pair.fact.context == &pair.context);
            } else {
                ASSERT_EQ(0u, required);
                ASSERT(pair.fact.context == &context);
            }
            ASSERT_EQ(DUCKVEP_HGVS_OK, status);
            ASSERT_EQ(0xa5u, rotated[0]); ASSERT_EQ(0xa5u, rotated[3]);
            char rendered[128];
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                &pair.fact, 0, rendered, sizeof rendered, &required));
            ASSERT_STR_EQ(expected[i], rendered);
        }
        ASSERT_MEM_EQ(&saved_context, &context, sizeof context);
        ASSERT_MEM_EQ(&saved_delta, &delta, sizeof delta);
        ASSERT_MEM_EQ(&saved_edit, &edit, sizeof edit);
        facts.transcript_edit = NULL;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_protein_pair_build(
            &s.tx, &s.ex, &sequences, &variants, &row, &facts, &dna, &reference,
            &scratch, NULL, 0u, &required, &pair));
        ASSERT_MEM_EQ(&zero, &pair.fact, sizeof zero);
        ASSERT_EQ(0u, required);
    }
    PASS();
}

TEST hgvs_retained_n_insertion_anchor_matches_vep(void) {
    /* Pinned Parser/VCF.pm removes the shared anchor before allele eligibility.
     * Actual VEP CLI observations retain these original N>NGCC records, with
     * GCN->A and NCN->X consensus translation. The canonical anchor is a control. */
    static const struct {
        const char *genome;
        uint32_t position1;
        const char *alleles;
        const char *protein;
    } cases[] = {
        {"AAAAAAAAAA" "ATGGCNGCCTAA" "AAAAAAAAAA", 16u, "NNGCC", "p.Ala2dup"},
        {"AAAAAAAAAA" "ATGNCNGCCTAA" "AAAAAAAAAA", 14u, "NNGCC", "p.Ter2_Ala3insPro"},
        {"AAAAAAAAAA" "ATGGCTGCCTAA" "AAAAAAAAAA", 16u, "TTGCC", "p.Ala2dup"}
    };
    for (size_t i = 0u; i < sizeof cases / sizeof cases[0]; i++) {
        struct kprop_proj_scene s = {0};
        s.chrom = 0u;
        s.tstart = s.cds_s = s.es[0] = 11u;
        s.tend = s.cds_e = s.ee[0] = 22u;
        s.cs[0] = 1u;
        s.ce[0] = 12u;
        s.strand = 1;
        s.excnt = 1u;
        s.flags = DUCKVEP_TX_HAS_TRANSLATION | DUCKVEP_TX_BIOTYPE_PROTEIN_CODING;
        kprop_proj_scene_finish(&s);
        uint64_t offset = 0u;
        uint32_t length = 12u, empty = 0u;
        uint8_t table = DUCKVEP_CODON_TABLE_STANDARD;
        duckvep_sequence_pool_t sequences = {0};
        sequences.cds_bytes = (const uint8_t *)cases[i].genome + 10u;
        sequences.cds_bytes_len = length;
        sequences.cds_offset = &offset;
        sequences.cds_length = &length;
        sequences.codon_table = &table;
        sequences.transcript_count = 1u;
        sequences.flank_bytes = sequences.cds_bytes;
        sequences.flank_bytes_len = length;
        sequences.pre_cds_offset = sequences.post_cds_offset = &offset;
        sequences.pre_cds_length = sequences.post_cds_length = &empty;
        sequences.flanks_complete = 1u;
        duckvep_hgvs_reference_window_t reference = {
            (const uint8_t *)cases[i].genome, strlen(cases[i].genome), 1u, 0u};
        uint32_t ro = 0u, ao = 1u;
        uint16_t rl = 1u, al = 4u;
        uint8_t kind = DUCKVEP_KIND_INS;
        duckvep_variant_batch_t variants = {0};
        variants.chrom_id = &s.chrom;
        variants.pos1 = variants.end1 = &cases[i].position1;
        variants.ref_offset = &ro;
        variants.alt_offset = &ao;
        variants.ref_length = &rl;
        variants.alt_length = &al;
        variants.variant_kind = &kind;
        variants.count = 1u;
        variants.allele_bytes = (const uint8_t *)cases[i].alleles;
        variants.allele_bytes_len = 5u;
        duckvep_event_t event;
        ASSERT(duckvep_event_prepare_small(cases[i].position1, variants.allele_bytes,
            rl, variants.allele_bytes + ao, al, &event));
        event.chrom_id = 0u;
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, variants.allele_bytes, rl));
        duckvep_haplotype_edit_t edits[4];
        duckvep_transcript_edit_t edit;
        ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK, duckvep_transcript_edit_build_prepared(
            &s.tx, &s.ex, &sequences, &variants, 0u, 0u, &event, edits, 4u, &edit));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, edit.cds_status);
        uint8_t cds[64], rp[32], ap[32];
        duckvep_delta_scratch_t scratch = {edits, 4u, cds, sizeof cds,
            rp, sizeof rp, ap, sizeof ap};
        duckvep_coding_context_t context;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK, duckvep_model_coding_context_build(
            &s.tx, &s.ex, &sequences, 0u, 1, &event, &edit.cds_edits,
            cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &context));
        duckvep_sequence_delta_t delta;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_delta_fill(&context, s.flags, &delta));
        duckvep_pair_facts_t facts = {0};
        facts.event = &event;
        facts.transcript_edit = &edit;
        facts.transcript_edit_status = DUCKVEP_TRANSCRIPT_EDIT_OK;
        facts.coding_context = &context;
        facts.delta = &delta;
        facts.projection_exon_hint = 0u;
        duckvep_consequence_t row = {0};
        row.overlap_object_kind = DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
        row.region_mask = DUCKVEP_REGION_CDS;
        row.flags = duckvep_sequence_delta_consequence_flags(&delta, 1);
        duckvep_hgvs_dna_fact_t dna;
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
            &s.tx, &s.ex, &reference, &reference, &edit, &dna));
        ASSERT_EQ(0, dna.shift_offset);
        for (uint8_t reuse_context = 0u; reuse_context < 2u; reuse_context++) {
            facts.coding_context_valid = reuse_context;
            uint8_t alleles[6];
            memset(alleles, 0xa5, sizeof alleles);
            duckvep_hgvs_protein_pair_t pair;
            size_t required;
            if (!reuse_context) {
                ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_pair_build(
                    &s.tx, &s.ex, &sequences, &variants, &row, &facts, &dna, &reference,
                    &scratch, alleles + 1u, 3u, &required, &pair));
                ASSERT_EQ(4u, required);
                for (size_t j = 0u; j < sizeof alleles; j++) ASSERT_EQ(0xa5u, alleles[j]);
            }
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_pair_build(
                &s.tx, &s.ex, &sequences, &variants, &row, &facts, &dna, &reference,
                &scratch, alleles + 1u, 4u, &required, &pair));
            ASSERT_EQ(reuse_context ? 0u : 4u, required);
            ASSERT(pair.fact.context == (reuse_context ? &context : &pair.context));
            ASSERT_EQ(0xa5u, alleles[0]);
            ASSERT_EQ(0xa5u, alleles[5]);
            char rendered[64];
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                &pair.fact, 0, rendered, sizeof rendered, &required));
            ASSERT_STR_EQ(cases[i].protein, rendered);
        }
        /* Permission to read the retained anchor does not admit an N in the
         * inserted payload or the genomic sequence used for shifting. */
        duckvep_hgvs_dna_fact_t unknown_payload = dna;
        unknown_payload.alt = (const uint8_t *)"GNC";
        facts.coding_context_valid = 0u;
        uint8_t allele_scratch[4];
        duckvep_hgvs_protein_pair_t pair;
        size_t required;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE, duckvep_hgvs_protein_pair_build(
            &s.tx, &s.ex, &sequences, &variants, &row, &facts, &unknown_payload, &reference,
            &scratch, allele_scratch, sizeof allele_scratch, &required, &pair));
        uint8_t unknown_shift_bases[32];
        ASSERT_EQ(sizeof unknown_shift_bases, reference.length);
        memcpy(unknown_shift_bases, reference.bases, sizeof unknown_shift_bases);
        unknown_shift_bases[0] = (uint8_t)'N';
        duckvep_hgvs_reference_window_t unknown_shift = reference;
        unknown_shift.bases = unknown_shift_bases;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE,
            duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                &s.tx, &s.ex, &unknown_shift, &reference, &edit, &dna));
    }
    PASS();
}

TEST hgvs_dna_facts_orient_reverse_alleles_once(void) {
    static const uint8_t ref[] = {'A', 'C'};
    static const uint8_t alt[] = {'G'};
    struct kprop_proj_scene s;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    uint8_t base;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)-1;
    s.excnt = 2u; s.cds_s = 104u; s.cds_e = 205u;
    s.es[0] = 200u; s.ee[0] = 209u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.es[1] = 100u; s.ee[1] = 109u; s.cs[1] = 11u; s.ce[1] = 20u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);
    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u;
    edit.first.cdna_anchor1 = 6u; edit.first.exonic = 1u;
    edit.last.cdna_anchor1 = 7u; edit.last.exonic = 1u;
    edit.ref = ref; edit.ref_length = 2u;
    edit.alt = alt; edit.alt_length = 1u;
    edit.transcript_strand = (int8_t)-1;

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build(&s.tx, &s.ex, &edit, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_NUMBERING_C, fact.numbering);
    ASSERT_EQ(DUCKVEP_HGVS_DNA_REPLACEMENT, fact.shape);
    ASSERT_EQ(2, fact.first.base);
    ASSERT_EQ(3, fact.last.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_base(&fact, 0, 0u, &base));
    ASSERT_EQ('G', base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_base(&fact, 0, 1u, &base));
    ASSERT_EQ('T', base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_base(&fact, 1, 0u, &base));
    ASSERT_EQ('C', base);
    {
        char rendered[64];
        char small[6];
        size_t required;
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_dna_render_basic(
                      &fact, rendered, sizeof rendered, &required));
        ASSERT_EQ(0, strcmp("c.2_3delinsC", rendered));
        ASSERT_EQ(strlen(rendered), required);
        ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL,
                  duckvep_hgvs_dna_render_basic(
                      &fact, small, sizeof small, &required));
        ASSERT_EQ(0, strcmp("c.2_3", small));
        ASSERT_EQ(strlen("c.2_3delinsC"), required);
    }
    PASS();
}

TEST hgvs_basic_renderer_known_coordinate_forms(void) {
    static const uint8_t ref_a[] = {'A'};
    static const uint8_t alt_g[] = {'G'};
    static const uint8_t alt_ac[] = {'A', 'C'};
    duckvep_hgvs_dna_fact_t fact;
    char rendered[64];
    size_t required;

    memset(&fact, 0, sizeof fact);
    fact.first.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C;
    fact.first.base = 1;
    fact.last = fact.first;
    fact.ref = ref_a; fact.ref_length = 1u;
    fact.alt = alt_g; fact.alt_length = 1u;
    fact.transcript_strand = (int8_t)1;
    fact.numbering = (uint8_t)DUCKVEP_HGVS_NUMBERING_C;
    fact.shape = (uint8_t)DUCKVEP_HGVS_DNA_SUBSTITUTION;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.1A>G", rendered));

    fact.first.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C_STAR;
    fact.first.base = 0;
    fact.first.intron_offset = 1;
    fact.last = fact.first;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.*1A>G", rendered));

    memset(&fact, 0, sizeof fact);
    fact.first.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_N;
    fact.first.base = 10;
    fact.last.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_N;
    fact.last.base = 11;
    fact.alt = alt_ac; fact.alt_length = 2u;
    fact.transcript_strand = (int8_t)1;
    fact.numbering = (uint8_t)DUCKVEP_HGVS_NUMBERING_N;
    fact.shape = (uint8_t)DUCKVEP_HGVS_DNA_INSERTION;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("n.10_11insAC", rendered));

    fact.shift_offset = -1;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    PASS();
}

TEST hgvs_dna_renderer_reports_long_output_before_retry(void) {
    enum { ALTERNATE_LENGTH = 1405, EXPECTED_PREFIX_LENGTH = 10 };
    uint8_t alternate[ALTERNATE_LENGTH];
    duckvep_hgvs_dna_fact_t fact;
    char small[256];
    char rendered[ALTERNATE_LENGTH + EXPECTED_PREFIX_LENGTH + 1u];
    size_t required;
    size_t i;

    memset(alternate, 'A', sizeof alternate);
    memset(&fact, 0, sizeof fact);
    fact.first.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C;
    fact.first.base = 10;
    fact.last.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C;
    fact.last.base = 11;
    fact.alt = alternate;
    fact.alt_length = ALTERNATE_LENGTH;
    fact.transcript_strand = (int8_t)1;
    fact.numbering = (uint8_t)DUCKVEP_HGVS_NUMBERING_C;
    fact.shape = (uint8_t)DUCKVEP_HGVS_DNA_INSERTION;

    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL,
              duckvep_hgvs_dna_render_basic(
                  &fact, small, sizeof small, &required));
    ASSERT_EQ((size_t)(EXPECTED_PREFIX_LENGTH + ALTERNATE_LENGTH), required);
    ASSERT_EQ(sizeof small - 1u, strlen(small));
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL,
              duckvep_hgvs_dna_render_basic(
                  &fact, NULL, 0u, &required));
    ASSERT_EQ((size_t)(EXPECTED_PREFIX_LENGTH + ALTERNATE_LENGTH), required);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, memcmp(rendered, "c.10_11ins", EXPECTED_PREFIX_LENGTH));
    for (i = EXPECTED_PREFIX_LENGTH; i < required; i++) {
        ASSERT_EQ('A', rendered[i]);
    }
    ASSERT_EQ('\0', rendered[required]);
    PASS();
}

TEST hgvs_shift_reproduces_vep_short_region_overlap(void) {
    static const uint8_t reference_bytes[] = {
        'C','A','A','A','A','A','G','G','G','G','G','G'
    };
    static const uint8_t deleted[] = {'A', 'A'};
    static const uint8_t inserted[] = {'A'};
    static const uint8_t inserted_g[] = {'G'};
    struct kprop_proj_scene s;
    duckvep_hgvs_reference_window_t reference;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    char rendered[64];
    size_t required;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 111u; s.strand = (int8_t)1;
    s.excnt = 1u; s.cds_s = 100u; s.cds_e = 111u;
    s.es[0] = 100u; s.ee[0] = 111u; s.cs[0] = 1u; s.ce[0] = 12u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);
    memset(&reference, 0, sizeof reference);
    reference.bases = reference_bytes;
    reference.length = sizeof reference_bytes;
    reference.start1 = 100u;
    reference.chrom_id = 0u;

    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u; edit.transcript_strand = (int8_t)1;
    edit.first.cdna_anchor1 = 2u; edit.first.exonic = 1u;
    edit.last.cdna_anchor1 = 3u; edit.last.exonic = 1u;
    edit.ref = deleted; edit.ref_length = 2u;
    edit.event.chrom_id = 0u;
    edit.event.start1 = 101u;
    edit.event.end1 = 102u;
    edit.event.feature_start1 = 101u;
    edit.event.feature_end1 = 102u;
    edit.event.ref_diff_length = 2u;
    edit.event.kind = (uint8_t)DUCKVEP_KIND_DEL;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted(
                  &s.tx, &s.ex, &reference, &edit, &fact));
    /* VEP takes both the first and last 1,000 bases of the constrained
     * sequence-region slice.  On this 12-base region those strings are the
     * same whole sequence, so the leading C stops a deletion that a clean
     * event-adjacent walk would have shifted across the following A run. */
    ASSERT_EQ(0, fact.shift_offset);
    ASSERT_EQ(2, fact.first.base);
    ASSERT_EQ(3, fact.last.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.2_3del", rendered));

    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u; edit.transcript_strand = (int8_t)1;
    edit.first.cdna_anchor1 = 1u; edit.first.exonic = 1u;
    edit.last.cdna_anchor1 = 2u; edit.last.exonic = 1u;
    edit.alt = inserted; edit.alt_length = 1u;
    edit.event.chrom_id = 0u;
    edit.event.start1 = 100u;
    edit.event.end1 = 100u;
    edit.event.insertion_boundary0 = 100u;
    edit.event.interbase = 1u;
    edit.event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    edit.event.feature_start1 = 101u;
    edit.event.feature_end1 = 100u;
    edit.event.alt_diff_length = 1u;
    edit.event.kind = (uint8_t)DUCKVEP_KIND_INS;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted(
                  &s.tx, &s.ex, &reference, &edit, &fact));
    ASSERT_EQ(0, fact.shift_offset);
    ASSERT_EQ(DUCKVEP_HGVS_DNA_DUPLICATION, fact.shape);
    ASSERT_EQ(2, fact.first.base);
    ASSERT_EQ(2, fact.last.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.2dup", rendered));

    /* At the transcript endpoint an ordinary insertion has no right mapper
     * flank. VEP can nevertheless name a duplication from the final copied
     * transcript base because dup syntax projects that source span. */
    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u; edit.transcript_strand = (int8_t)1;
    edit.first.genomic_pos1 = 111u;
    edit.first.cdna_anchor1 = 12u;
    edit.first.exonic = 1u;
    edit.last = edit.first;
    edit.alt = inserted_g; edit.alt_length = 1u;
    edit.event.chrom_id = 0u;
    edit.event.start1 = 111u;
    edit.event.end1 = 111u;
    edit.event.insertion_boundary0 = 111u;
    edit.event.interbase = 1u;
    edit.event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    edit.event.feature_start1 = 112u;
    edit.event.feature_end1 = 111u;
    edit.event.alt_diff_length = 1u;
    edit.event.kind = (uint8_t)DUCKVEP_KIND_INS;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted(
                  &s.tx, &s.ex, &reference, &edit, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_DNA_DUPLICATION, fact.shape);
    ASSERT_EQ(12, fact.first.base);
    ASSERT_EQ(12, fact.last.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.12dup", rendered));
    PASS();
}

TEST hgvs_large_duplication_uses_lookup_beyond_shift_slice(void) {
    enum { TRANSCRIPT_LENGTH = 4000, DUPLICATION_LENGTH = 1001 };
    uint8_t reference_bytes[TRANSCRIPT_LENGTH];
    uint8_t inserted[DUPLICATION_LENGTH];
    struct kprop_proj_scene s;
    duckvep_hgvs_reference_window_t shift_reference;
    duckvep_hgvs_reference_window_t lookup_reference;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    uint32_t shift_start1;
    uint32_t shift_end1;
    uint32_t lookup_start1;
    uint32_t lookup_end1;
    char rendered[64];
    size_t required;

    memset(reference_bytes, 'G', sizeof reference_bytes);
    memset(inserted, 'A', sizeof inserted);
    /* The copied 5-prime source starts one base before VEP's exact shift
     * slice. The first 3-prime comparison base is C, so the insertion does
     * not shift before hgvs_variant_notation checks both adjacent sources. */
    memset(reference_bytes + 1499u, 'A', DUPLICATION_LENGTH);
    reference_bytes[2500u] = 'C';

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 1000u; s.tend = 4999u;
    s.strand = (int8_t)1; s.excnt = 1u;
    s.cds_s = 1000u; s.cds_e = 4999u;
    s.es[0] = 1000u; s.ee[0] = 4999u;
    s.cs[0] = 1u; s.ce[0] = TRANSCRIPT_LENGTH; s.phase[0] = 0;
    kprop_proj_scene_finish(&s);

    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u; edit.transcript_strand = (int8_t)1;
    edit.first.genomic_pos1 = 3499u;
    edit.first.cdna_anchor1 = 2500u; edit.first.exonic = 1u;
    edit.last.genomic_pos1 = 3500u;
    edit.last.cdna_anchor1 = 2501u; edit.last.exonic = 1u;
    edit.alt = inserted; edit.alt_length = DUPLICATION_LENGTH;
    edit.event.chrom_id = 0u;
    edit.event.start1 = 3499u; edit.event.end1 = 3499u;
    edit.event.insertion_boundary0 = 3499u;
    edit.event.interbase = 1u;
    edit.event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    edit.event.feature_start1 = 3500u;
    edit.event.feature_end1 = 3499u;
    edit.event.alt_diff_length = DUPLICATION_LENGTH;
    edit.event.kind = (uint8_t)DUCKVEP_KIND_INS;

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_genomic_search_interval(
                  &edit.event, TRANSCRIPT_LENGTH + 999u,
                  &shift_start1, &shift_end1));
    ASSERT_EQ(2500u, shift_start1);
    ASSERT_EQ(4499u, shift_end1);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_reference_fetch_interval(
                  &edit.event, TRANSCRIPT_LENGTH + 999u,
                  &lookup_start1, &lookup_end1));
    ASSERT_EQ(1499u, lookup_start1);
    ASSERT_EQ(4999u, lookup_end1);

    memset(&lookup_reference, 0, sizeof lookup_reference);
    lookup_reference.bases = reference_bytes;
    lookup_reference.length = sizeof reference_bytes;
    lookup_reference.start1 = 1000u;
    lookup_reference.chrom_id = 0u;
    shift_reference.bases = reference_bytes +
        (size_t)(shift_start1 - lookup_reference.start1);
    shift_reference.length =
        (size_t)((uint64_t)shift_end1 - shift_start1 + 1u);
    shift_reference.start1 = shift_start1;
    shift_reference.chrom_id = 0u;

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                  &s.tx, &s.ex, &shift_reference, &lookup_reference,
                  &edit, &fact));
    ASSERT_EQ(0, fact.shift_offset);
    ASSERT_EQ(DUCKVEP_HGVS_DNA_DUPLICATION, fact.shape);
    ASSERT_EQ(1500, fact.first.base);
    ASSERT_EQ(2500, fact.last.base);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.1500_2500dup", rendered));
    PASS();
}

/* Randomized HGVS 3-prime shift oracle. The generator owns one forward genomic
 * reference window and its transcript-oriented view, then stores the uploaded
 * allele in genomic orientation. The oracle walks the byte arrays directly;
 * it does not call projection, allele-orientation, or window helpers. */
#define KPROP_HGVS_SHIFT_MAX_TRANSCRIPT 72u
#define KPROP_HGVS_SHIFT_MAX_PATTERN 6u

struct kprop_hgvs_shift_scene {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_sequence_pool_t seq;
    duckvep_hgvs_reference_window_t reference;
    duckvep_transcript_edit_t edit;

    uint16_t chrom;
    uint32_t tx_start;
    uint32_t tx_end;
    int8_t strand;
    uint64_t flags;
    uint32_t exon_offset;
    uint16_t exon_count;
    uint32_t cds_start;
    uint32_t cds_end;
    uint32_t exon_start;
    uint32_t exon_end;
    uint32_t cdna_start;
    uint32_t cdna_end;
    int8_t phase;
    int8_t end_phase;
    uint64_t cds_offset;
    uint32_t cds_length;
    uint8_t codon_table;
    uint64_t pre_cds_offset;
    uint64_t post_cds_offset;
    uint32_t pre_cds_length;
    uint32_t post_cds_length;

    uint8_t transcript[KPROP_HGVS_SHIFT_MAX_TRANSCRIPT];
    uint8_t genomic_reference[KPROP_HGVS_SHIFT_MAX_TRANSCRIPT];
    uint8_t genomic_allele[KPROP_HGVS_SHIFT_MAX_PATTERN];
    uint8_t oriented_allele[KPROP_HGVS_SHIFT_MAX_PATTERN];

    uint32_t transcript_length;
    uint32_t first_cdna;
    uint32_t last_cdna;
    uint16_t allele_length;
    uint8_t insertion;
};

uint8_t kprop_hgvs_complement(uint8_t base) {
    switch (base) {
        case (uint8_t)'A': return (uint8_t)'T';
        case (uint8_t)'C': return (uint8_t)'G';
        case (uint8_t)'G': return (uint8_t)'C';
        case (uint8_t)'T': return (uint8_t)'A';
        default: return 0u;
    }
}

static void kprop_hgvs_shift_finish(struct kprop_hgvs_shift_scene *s) {
    s->tx.chrom_id = &s->chrom;
    s->tx.start1 = &s->tx_start;
    s->tx.end1 = &s->tx_end;
    s->tx.strand = &s->strand;
    s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exon_offset;
    s->tx.exon_count = &s->exon_count;
    s->tx.cds_start1 = &s->cds_start;
    s->tx.cds_end1 = &s->cds_end;
    s->tx.transcript_count = 1u;

    s->ex.start1 = &s->exon_start;
    s->ex.end1 = &s->exon_end;
    s->ex.cdna_start1 = &s->cdna_start;
    s->ex.cdna_end1 = &s->cdna_end;
    s->ex.phase = &s->phase;
    s->ex.end_phase = &s->end_phase;
    s->ex.exon_count = 1u;

    s->reference.bases = s->genomic_reference;
    s->reference.length = s->transcript_length;
    s->reference.start1 = s->tx_start;
    s->reference.chrom_id = s->chrom;

    s->cds_offset = 0u;
    s->cds_length = s->transcript_length;
    s->codon_table = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->seq.cds_bytes = s->transcript;
    s->seq.cds_bytes_len = s->transcript_length;
    s->seq.cds_offset = &s->cds_offset;
    s->seq.cds_length = &s->cds_length;
    s->seq.codon_table = &s->codon_table;
    s->seq.transcript_count = 1u;
    s->seq.pre_cds_offset = &s->pre_cds_offset;
    s->seq.pre_cds_length = &s->pre_cds_length;
    s->seq.post_cds_offset = &s->post_cds_offset;
    s->seq.post_cds_length = &s->post_cds_length;
    s->seq.flanks_complete = 1u;
}

static enum theft_alloc_res kprop_hgvs_shift_alloc(
    struct theft *t, void *env, void **instance) {

    static const uint8_t BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_hgvs_shift_scene *s =
        (struct kprop_hgvs_shift_scene *)calloc(1u, sizeof *s);
    uint32_t i;
    uint32_t repeat;
    uint32_t repeat_cap;
    uint32_t repeat_start;
    uint32_t mismatch_position;
    uint32_t genomic_first;
    uint32_t genomic_last;
    uint32_t genomic_low;
    uint32_t genomic_high;
    int force_terminal_duplication;
    (void)env;

    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->transcript_length =
        ((uint32_t)kprop_bounded(t, 19u) + 6u) * 3u;
    s->allele_length =
        (uint16_t)((uint32_t)kprop_bounded(t,
            KPROP_HGVS_SHIFT_MAX_PATTERN) + 1u);
    force_terminal_duplication = kprop_bounded(t, 32u) == 0u;
    s->insertion = force_terminal_duplication
        ? (uint8_t)1u : (uint8_t)kprop_bounded(t, 2u);
    s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;

    for (i = 0u; i < s->transcript_length; i++) {
        s->transcript[i] = BASES[kprop_bounded(t, 4u)];
    }
    for (i = 0u; i < (uint32_t)s->allele_length; i++) {
        s->oriented_allele[i] = BASES[kprop_bounded(t, 4u)];
    }

    if (s->insertion != 0u) {
        s->first_cdna = force_terminal_duplication
            ? s->transcript_length / 2u
            : (uint32_t)kprop_bounded(
                t, s->transcript_length - 3u) + 2u;
        s->last_cdna = s->first_cdna + 1u;
        repeat_start = s->last_cdna;
        repeat_cap = s->transcript_length - s->first_cdna;
    } else {
        s->first_cdna = (uint32_t)kprop_bounded(
            t, s->transcript_length - (uint32_t)s->allele_length - 2u) + 2u;
        s->last_cdna = s->first_cdna + (uint32_t)s->allele_length - 1u;
        for (i = 0u; i < (uint32_t)s->allele_length; i++) {
            s->transcript[s->first_cdna - 1u + i] = s->oriented_allele[i];
        }
        repeat_start = s->last_cdna + 1u;
        repeat_cap = s->transcript_length - s->last_cdna;
    }

    /* Keep a mismatching reference base inside the transcript so every
     * generated shifted event remains projectable. The VEP loop-limit cases
     * still stop before that mismatch when allele_length > 1. */
    repeat = (uint32_t)kprop_bounded(t, repeat_cap);
    for (i = 0u; i < repeat; i++) {
        s->transcript[repeat_start - 1u + i] =
            s->oriented_allele[i % (uint32_t)s->allele_length];
    }
    mismatch_position = repeat_start + repeat;
    if (mismatch_position <= s->transcript_length) {
        uint8_t expected =
            s->oriented_allele[repeat % (uint32_t)s->allele_length];
        uint8_t replacement = BASES[kprop_bounded(t, 3u)];
        if (replacement == expected) replacement = (uint8_t)'T';
        if (replacement == expected) replacement = (uint8_t)'G';
        s->transcript[mismatch_position - 1u] = replacement;
    }
    if (force_terminal_duplication) {
        uint32_t target_shift =
            s->transcript_length - s->first_cdna;
        uint8_t expected = s->oriented_allele[
            target_shift % (uint32_t)s->allele_length];
        uint8_t replacement = BASES[kprop_bounded(t, 3u)];

        /* VEP's short-region pre/post slices overlap. Construct a dedicated
         * stratum whose byte walk reaches the one shift placing the insertion
         * just beyond the transcript, while the copied source ending at the
         * transcript endpoint remains a printable duplication. This state was
         * previously reached only a few times per 100,000 random scenes. */
        for (i = 0u; i < target_shift; i++) {
            s->transcript[i] = s->oriented_allele[
                i % (uint32_t)s->allele_length];
        }
        if (replacement == expected) replacement = (uint8_t)'T';
        if (replacement == expected) replacement = (uint8_t)'G';
        s->transcript[target_shift] = replacement;
        for (i = 0u; i < (uint32_t)s->allele_length; i++) {
            s->transcript[
                s->transcript_length - (uint32_t)s->allele_length + i
            ] = s->oriented_allele[
                (i + target_shift) % (uint32_t)s->allele_length
            ];
        }
    }

    for (i = 0u; i < (uint32_t)s->allele_length; i++) {
        s->genomic_allele[i] = s->strand > 0
            ? s->oriented_allele[i]
            : kprop_hgvs_complement(
                s->oriented_allele[(uint32_t)s->allele_length - 1u - i]);
    }
    for (i = 0u; i < s->transcript_length; i++) {
        s->genomic_reference[i] = s->strand > 0
            ? s->transcript[i]
            : kprop_hgvs_complement(
                s->transcript[s->transcript_length - 1u - i]);
    }

    s->chrom = 0u;
    s->tx_start = 1000u;
    s->tx_end = s->tx_start + s->transcript_length - 1u;
    s->flags = 0u;
    s->exon_offset = 0u;
    s->exon_count = 1u;
    s->cds_start = s->tx_start;
    s->cds_end = s->tx_end;
    s->exon_start = s->tx_start;
    s->exon_end = s->tx_end;
    s->cdna_start = 1u;
    s->cdna_end = s->transcript_length;
    s->phase = 0;
    s->end_phase = 0;
    kprop_hgvs_shift_finish(s);

    memset(&s->edit, 0, sizeof s->edit);
    s->edit.tx_idx = 0u;
    s->edit.transcript_strand = s->strand;
    s->edit.first.cdna_anchor1 = s->first_cdna;
    s->edit.first.exonic = 1u;
    s->edit.last.cdna_anchor1 = s->last_cdna;
    s->edit.last.exonic = 1u;
    genomic_first = s->strand > 0
        ? s->tx_start + s->first_cdna - 1u
        : s->tx_end - s->first_cdna + 1u;
    genomic_last = s->strand > 0
        ? s->tx_start + s->last_cdna - 1u
        : s->tx_end - s->last_cdna + 1u;
    genomic_low = genomic_first < genomic_last ? genomic_first : genomic_last;
    genomic_high = genomic_first > genomic_last ? genomic_first : genomic_last;
    s->edit.event.chrom_id = s->chrom;
    if (s->insertion != 0u) {
        s->edit.alt = s->genomic_allele;
        s->edit.alt_length = s->allele_length;
        s->edit.event.start1 = genomic_low;
        s->edit.event.end1 = genomic_low;
        s->edit.event.insertion_boundary0 = genomic_low;
        s->edit.event.interbase = 1u;
        s->edit.event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
        s->edit.event.feature_start1 = genomic_low + 1u;
        s->edit.event.feature_end1 = genomic_low;
        s->edit.event.alt_diff_length = s->allele_length;
        s->edit.event.kind = (uint8_t)DUCKVEP_KIND_INS;
    } else {
        s->edit.ref = s->genomic_allele;
        s->edit.ref_length = s->allele_length;
        s->edit.event.start1 = genomic_low;
        s->edit.event.end1 = genomic_high;
        s->edit.event.feature_start1 = genomic_low;
        s->edit.event.feature_end1 = genomic_high;
        s->edit.event.ref_diff_length = s->allele_length;
        s->edit.event.kind = (uint8_t)DUCKVEP_KIND_DEL;
    }
    *instance = s;
    return THEFT_ALLOC_OK;
}

static void kprop_hgvs_shift_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static struct theft_type_info kprop_hgvs_shift_info = {
    .alloc = kprop_hgvs_shift_alloc,
    .free = kprop_hgvs_shift_free,
};

static uint64_t kprop_hgvs_shift_forward_count;
static uint64_t kprop_hgvs_shift_reverse_count;
static uint64_t kprop_hgvs_shift_insertion_count;
static uint64_t kprop_hgvs_shift_deletion_count;
static uint64_t kprop_hgvs_shift_duplication_count;
static uint64_t kprop_hgvs_shift_at_vep_limit_count;
static uint64_t kprop_hgvs_shift_rotated_count;
static uint64_t kprop_hgvs_shift_composed_count;
static uint64_t kprop_hgvs_shift_protein_count;
static uint64_t kprop_hgvs_shift_nonlocal_ref_replay_count;
static uint64_t kprop_hgvs_shift_terminal_duplication_count;

static enum theft_trial_res prop_hgvs_shift_matches_genomic_oracle(
    struct theft *t, void *arg1) {

    const struct kprop_hgvs_shift_scene *s =
        (const struct kprop_hgvs_shift_scene *)arg1;
    duckvep_hgvs_dna_fact_t fact;
    uint32_t shift = 0u;
    uint32_t available;
    uint32_t limit;
    uint32_t expected_first;
    uint32_t expected_last;
    uint8_t expected_shape;
    uint32_t i;
    uint8_t allele_scratch[KPROP_HGVS_SHIFT_MAX_PATTERN + 1u];
    uint8_t applied_cds[KPROP_HGVS_SHIFT_MAX_TRANSCRIPT +
                        KPROP_HGVS_SHIFT_MAX_PATTERN];
    uint8_t expected_cds[KPROP_HGVS_SHIFT_MAX_TRANSCRIPT +
                         KPROP_HGVS_SHIFT_MAX_PATTERN];
    uint8_t context_alt_cds[KPROP_HGVS_SHIFT_MAX_TRANSCRIPT +
                            KPROP_HGVS_SHIFT_MAX_PATTERN];
    uint8_t context_ref_peptide[KPROP_HGVS_SHIFT_MAX_TRANSCRIPT / 3u + 2u];
    uint8_t context_alt_peptide[(KPROP_HGVS_SHIFT_MAX_TRANSCRIPT +
                                 KPROP_HGVS_SHIFT_MAX_PATTERN) / 3u + 2u];
    duckvep_event_t shifted_event;
    duckvep_haplotype_edit_t shifted_edit;
    duckvep_haplotype_result_t apply_result;
    duckvep_edit_set_t edit_set;
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t protein_fact;
    size_t required = 0u;
    size_t applied_length = 0u;
    size_t expected_length;
    size_t prefix_length;
    size_t suffix_start;
    size_t expected_peptide_length;
    int shifted_reference_matches = 1;
    (void)t;
    memset(&fact, 0, sizeof fact);

#define KPROP_HGVS_SHIFT_FAIL() do { \
    fprintf(stderr, \
        "[HGVS shift oracle failure] line=%d ins=%u strand=%d" \
        " transcript_len=%u allele_len=%u start=%u_%u shift=%u limit=%u" \
        " want_shape=%u got_shape=%u got_shift=%d" \
        " transcript=%.*s genomic=%.*s allele=%.*s oriented=%.*s\n", \
        __LINE__, (unsigned)s->insertion, (int)s->strand, \
        s->transcript_length, (unsigned)s->allele_length, \
        s->first_cdna, s->last_cdna, shift, limit, \
        (unsigned)expected_shape, (unsigned)fact.shape, fact.shift_offset, \
        (int)s->transcript_length, (const char *)s->transcript, \
        (int)s->transcript_length, (const char *)s->genomic_reference, \
        (int)s->allele_length, (const char *)s->genomic_allele, \
        (int)s->allele_length, (const char *)s->oriented_allele); \
    return THEFT_TRIAL_FAIL; \
} while (0)

    /* Independent reproduction of TranscriptVariationAllele::_genomic_shift:
     * after constraining the event +/- 1,000 slice to the sequence region,
     * VEP names its first 1,000 bases pre_seq and its last 1,000 bases
     * post_seq.  These overlap completely in this deliberately short random
     * sequence region.  Do not replace this with a transcript-relative walk;
     * the overlapping pre/post-slice behavior is observable VEP 116 output. */
    available = s->transcript_length;
    limit = available < (uint32_t)s->allele_length
        ? available
        : available - (uint32_t)s->allele_length + 1u;
    while (shift < limit) {
        uint32_t reference_index = s->strand > 0
            ? shift : s->transcript_length - 1u - shift;
        uint32_t allele_index = s->strand > 0
            ? shift % (uint32_t)s->allele_length
            : (uint32_t)s->allele_length - 1u -
                shift % (uint32_t)s->allele_length;
        if (s->genomic_reference[reference_index] !=
            s->genomic_allele[allele_index]) {
            break;
        }
        shift++;
    }
    expected_first = s->first_cdna + shift;
    expected_last = s->last_cdna + shift;
    expected_shape = s->insertion != 0u
        ? (uint8_t)DUCKVEP_HGVS_DNA_INSERTION
        : (uint8_t)DUCKVEP_HGVS_DNA_DELETION;

    if (s->insertion != 0u) {
        int direction;
        for (direction = 0; direction < 2; direction++) {
            uint32_t source_start;
            uint32_t source_end;
            int duplicated = 1;
            if (direction == 0) {
                source_start = s->last_cdna + shift;
                if (source_start > s->transcript_length ||
                    (uint32_t)s->allele_length >
                        s->transcript_length - source_start + 1u) {
                    continue;
                }
                source_end = source_start +
                    (uint32_t)s->allele_length - 1u;
            } else {
                source_end = s->first_cdna + shift;
                if (source_end > s->transcript_length ||
                    (uint32_t)s->allele_length > source_end) {
                    continue;
                }
                source_start = source_end -
                    (uint32_t)s->allele_length + 1u;
            }
            for (i = 0u; i < (uint32_t)s->allele_length; i++) {
                uint8_t rotated = s->oriented_allele[
                    (i + shift) % (uint32_t)s->allele_length];
                if (s->transcript[source_start - 1u + i] != rotated) {
                    duplicated = 0;
                    break;
                }
            }
            if (duplicated) {
                expected_shape = (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION;
                expected_first = source_start;
                expected_last = source_end;
                break;
            }
        }
    }

    {
        duckvep_hgvs_status_t got_status =
            duckvep_hgvs_dna_fact_build_genomic_shifted(
                &s->tx, &s->ex, &s->reference, &s->edit, &fact);
        duckvep_hgvs_status_t expected_status =
            expected_last > s->transcript_length
                ? DUCKVEP_HGVS_NOT_APPLICABLE : DUCKVEP_HGVS_OK;
        if (got_status != expected_status) {
            fprintf(stderr,
                "[HGVS shift status mismatch] status=%u want=%u shift=%u"
                " ins=%u strand=%d len=%u allele=%u start=%u_%u\n",
                (unsigned)got_status, (unsigned)expected_status, shift,
                (unsigned)s->insertion, (int)s->strand,
                s->transcript_length, (unsigned)s->allele_length,
                s->first_cdna, s->last_cdna);
            KPROP_HGVS_SHIFT_FAIL();
        }
        if (expected_status == DUCKVEP_HGVS_NOT_APPLICABLE) {
            return THEFT_TRIAL_PASS;
        }
        if (got_status != DUCKVEP_HGVS_OK ||
        fact.shift_offset != (int32_t)shift ||
        fact.shape != expected_shape ||
        fact.first.kind != (uint8_t)DUCKVEP_HGVS_COORDINATE_C ||
        fact.last.kind != (uint8_t)DUCKVEP_HGVS_COORDINATE_C ||
        fact.first.base != (int64_t)expected_first ||
        fact.last.base != (int64_t)expected_last) {
            fprintf(stderr,
                "[HGVS shift mismatch] status=%u want_shift=%u got_shift=%d"
                " ins=%u strand=%d len=%u allele=%u want_shape=%u got_shape=%u"
                " start=%u_%u want=%u_%u got=%" PRId64 "_%" PRId64 "\n",
                (unsigned)got_status, shift, fact.shift_offset,
                (unsigned)s->insertion, (int)s->strand,
                s->transcript_length, (unsigned)s->allele_length,
                (unsigned)expected_shape, (unsigned)fact.shape,
                s->first_cdna, s->last_cdna, expected_first, expected_last,
                fact.first.base, fact.last.base);
            KPROP_HGVS_SHIFT_FAIL();
        }
    }
    for (i = 0u; i < (uint32_t)s->allele_length; i++) {
        uint8_t got = 0u;
        uint8_t expected = s->oriented_allele[
            (i + shift) % (uint32_t)s->allele_length];
        if (duckvep_hgvs_dna_base(
                &fact, s->insertion != 0u, i, &got) != DUCKVEP_HGVS_OK ||
            got != expected) {
            KPROP_HGVS_SHIFT_FAIL();
        }
    }

    /* A too-small allele buffer must not publish a partial genomic event or
     * CDS edit. This remains part of the same randomized distribution instead
     * of a separate hand-picked capacity example. */
    memset(&shifted_event, 0xA5, sizeof shifted_event);
    memset(&shifted_edit, 0xA5, sizeof shifted_edit);
    if (duckvep_hgvs_shifted_cds_edit_build(
            &s->tx, &s->ex, &s->seq, &s->reference, &s->edit, &fact,
            allele_scratch, (size_t)s->allele_length -
                (s->insertion == 0u ? 1u : 0u),
            &required, &shifted_event, &shifted_edit) !=
            DUCKVEP_HGVS_BUFFER_TOO_SMALL ||
        required != (size_t)s->allele_length +
            (s->insertion != 0u ? 1u : 0u) ||
        shifted_event.start1 != 0u || shifted_edit.cds_start != 0u) {
        KPROP_HGVS_SHIFT_FAIL();
    }
    if (s->insertion == 0u) {
        uint32_t expected_cds_start = s->first_cdna + shift;
        for (i = 0u; i < (uint32_t)s->allele_length; i++) {
            uint8_t expected = s->oriented_allele[
                (i + shift) % (uint32_t)s->allele_length];
            if (s->transcript[expected_cds_start - 1u + i] != expected) {
                shifted_reference_matches = 0;
                break;
            }
        }
    }
    if (s->insertion != 0u &&
        (uint64_t)s->last_cdna + (uint64_t)shift >
            (uint64_t)s->transcript_length) {
        /* A copied source can remain printable as c.*dup even when the
         * shifted insertion point itself is beyond this synthetic transcript.
         * Production checks the shifted genomic2pep mapper endpoints before
         * attempting CDS composition, so the low-level composition status is
         * intentionally outside this property's contract. The reference view
         * here is only transcript-sized and would also make that status depend
         * on which genomic VCF anchor happens to be retained. */
        if (fact.shape != (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION) {
            KPROP_HGVS_SHIFT_FAIL();
        }
        kprop_hgvs_shift_terminal_duplication_count++;
        if (s->strand > 0) kprop_hgvs_shift_forward_count++;
        else kprop_hgvs_shift_reverse_count++;
        kprop_hgvs_shift_insertion_count++;
        kprop_hgvs_shift_duplication_count++;
        if (shift == limit) kprop_hgvs_shift_at_vep_limit_count++;
        if (shift % (uint32_t)s->allele_length != 0u) {
            kprop_hgvs_shift_rotated_count++;
        }
        return THEFT_TRIAL_PASS;
    }
    {
        duckvep_hgvs_status_t shifted_status =
            duckvep_hgvs_shifted_cds_edit_build(
                &s->tx, &s->ex, &s->seq, &s->reference, &s->edit, &fact,
                allele_scratch, sizeof allele_scratch, &required,
                &shifted_event, &shifted_edit);
        if (shifted_status != DUCKVEP_HGVS_OK) {
            KPROP_HGVS_SHIFT_FAIL();
        }
        if (!shifted_reference_matches) {
            kprop_hgvs_shift_nonlocal_ref_replay_count++;
        }
    }
    if (s->insertion != 0u) {
        uint32_t expected_cds_start = s->last_cdna + shift;
        if (shifted_event.insertion_boundary0 !=
                fact.placed_insertion_boundary0 ||
            shifted_edit.cds_start != expected_cds_start ||
            shifted_edit.ref_len != 0u ||
            shifted_edit.alt_len != (uint32_t)s->allele_length) {
            KPROP_HGVS_SHIFT_FAIL();
        }
        prefix_length = (size_t)expected_cds_start - 1u;
        suffix_start = prefix_length;
    } else {
        uint32_t expected_cds_start = s->first_cdna + shift;
        if (shifted_event.start1 != fact.placed_start1 ||
            shifted_event.end1 != fact.placed_end1 ||
            shifted_edit.cds_start != expected_cds_start ||
            shifted_edit.ref_len != (uint32_t)s->allele_length ||
            shifted_edit.alt_len != 0u) {
            KPROP_HGVS_SHIFT_FAIL();
        }
        prefix_length = (size_t)expected_cds_start - 1u;
        suffix_start = prefix_length + (size_t)s->allele_length;
    }
    expected_length = prefix_length +
        (s->insertion != 0u ? (size_t)s->allele_length : 0u) +
        ((size_t)s->transcript_length - suffix_start);
    memcpy(expected_cds, s->transcript, prefix_length);
    if (s->insertion != 0u) {
        for (i = 0u; i < (uint32_t)s->allele_length; i++) {
            expected_cds[prefix_length + i] = s->oriented_allele[
                (i + shift) % (uint32_t)s->allele_length];
        }
    }
    memcpy(expected_cds + prefix_length +
               (s->insertion != 0u ? (size_t)s->allele_length : 0u),
           s->transcript + suffix_start,
           (size_t)s->transcript_length - suffix_start);
    if (duckvep_haplotype_apply_cds_edits(
            s->transcript, s->transcript_length, &shifted_edit, 1u,
            s->strand, applied_cds, sizeof applied_cds, &applied_length,
            &apply_result) != DUCKVEP_HAPLOTYPE_OK ||
        applied_length != expected_length ||
        memcmp(applied_cds, expected_cds, expected_length) != 0) {
        KPROP_HGVS_SHIFT_FAIL();
    }

    edit_set.edits = &shifted_edit;
    edit_set.count = 1u;
    if (duckvep_model_coding_context_build(
            &s->tx, &s->ex, &s->seq, 0u, s->strand, &shifted_event,
            &edit_set, context_alt_cds, sizeof context_alt_cds,
            context_ref_peptide, sizeof context_ref_peptide,
            context_alt_peptide, sizeof context_alt_peptide, &context) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        KPROP_HGVS_SHIFT_FAIL();
    }
    expected_peptide_length = expected_length / 3u;
    if (context.alt_peptide_len != expected_peptide_length) {
        KPROP_HGVS_SHIFT_FAIL();
    }
    for (i = 0u; i < expected_peptide_length; i++) {
        char codon[3];
        char expected_amino_acid;
        codon[0] = (char)expected_cds[(size_t)i * 3u];
        codon[1] = (char)expected_cds[(size_t)i * 3u + 1u];
        codon[2] = (char)expected_cds[(size_t)i * 3u + 2u];
        expected_amino_acid = duckvep_translate_codon(
            codon, DUCKVEP_CODON_TABLE_STANDARD);
        if (duckvep_coding_context_peptide_base(&context, 1, i) !=
            (uint8_t)expected_amino_acid) {
            KPROP_HGVS_SHIFT_FAIL();
        }
    }
    kprop_hgvs_shift_composed_count++;
    if (duckvep_coding_context_delta_fill(&context, 0u, &delta) ==
            DUCKVEP_CONTEXT_DELTA_OK && delta.valid != 0u &&
        duckvep_hgvs_protein_fact_build(&context, &delta, &protein_fact) ==
            DUCKVEP_HGVS_OK) {
        kprop_hgvs_shift_protein_count++;
    }

    if (s->strand > 0) kprop_hgvs_shift_forward_count++;
    else kprop_hgvs_shift_reverse_count++;
    if (s->insertion != 0u) kprop_hgvs_shift_insertion_count++;
    else kprop_hgvs_shift_deletion_count++;
    if (expected_shape == (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION) {
        kprop_hgvs_shift_duplication_count++;
    }
    if (shift == limit) kprop_hgvs_shift_at_vep_limit_count++;
    if (shift % (uint32_t)s->allele_length != 0u) {
        kprop_hgvs_shift_rotated_count++;
    }
#undef KPROP_HGVS_SHIFT_FAIL
    return THEFT_TRIAL_PASS;
}

TEST hgvs_genomic_shift_matches_reference_oracle_for_any_indel(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    kprop_hgvs_shift_forward_count = 0u;
    kprop_hgvs_shift_reverse_count = 0u;
    kprop_hgvs_shift_insertion_count = 0u;
    kprop_hgvs_shift_deletion_count = 0u;
    kprop_hgvs_shift_duplication_count = 0u;
    kprop_hgvs_shift_at_vep_limit_count = 0u;
    kprop_hgvs_shift_rotated_count = 0u;
    kprop_hgvs_shift_composed_count = 0u;
    kprop_hgvs_shift_protein_count = 0u;
    kprop_hgvs_shift_nonlocal_ref_replay_count = 0u;
    kprop_hgvs_shift_terminal_duplication_count = 0u;
    cfg.name = "HGVS genomic 3-prime shift == independent reference byte-walk";
    cfg.prop1 = prop_hgvs_shift_matches_genomic_oracle;
    cfg.type_info[0] = &kprop_hgvs_shift_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    fprintf(stderr,
        "[HGVS shift coverage] fwd=%" PRIu64 " rev=%" PRIu64
        " ins=%" PRIu64 " del=%" PRIu64 " dup=%" PRIu64
        " at_vep_limit=%" PRIu64 " rotated=%" PRIu64
        " composed=%" PRIu64 " protein=%" PRIu64
        " nonlocal_ref_replay=%" PRIu64
        " terminal_duplication=%" PRIu64 "\n",
        kprop_hgvs_shift_forward_count, kprop_hgvs_shift_reverse_count,
        kprop_hgvs_shift_insertion_count, kprop_hgvs_shift_deletion_count,
        kprop_hgvs_shift_duplication_count,
        kprop_hgvs_shift_at_vep_limit_count,
        kprop_hgvs_shift_rotated_count,
        kprop_hgvs_shift_composed_count,
        kprop_hgvs_shift_protein_count,
        kprop_hgvs_shift_nonlocal_ref_replay_count,
        kprop_hgvs_shift_terminal_duplication_count);
    ASSERT(kprop_hgvs_shift_nonlocal_ref_replay_count > 0u);
    ASSERT(kprop_hgvs_shift_terminal_duplication_count > 0u);
    PASS();
}

TEST hgvs_clamped_feature_preserves_vep_preclip_multiplication_order(void) {
    static const uint8_t genomic_reference[] = {
        'A','T','G', 'A','A','A', 'T','A','A', 'C','G','T','A'
    };
    static const uint8_t alt_cc[] = {'C', 'C'};
    static const uint8_t alt_ccc[] = {'C', 'C', 'C'};
    static const uint8_t alt_acc[] = {'A', 'C', 'C'};
    static const uint8_t alt_cac[] = {'C', 'A', 'C'};
    struct kprop_proj_scene s;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_reference_window_t reference;
    duckvep_hgvs_dna_fact_t fact;
    char rendered[64];
    size_t required = 0u;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 109u; s.strand = (int8_t)1;
    s.excnt = 1u; s.cds_s = 100u; s.cds_e = 106u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);
    memset(&reference, 0, sizeof reference);
    reference.bases = genomic_reference;
    reference.length = sizeof genomic_reference;
    reference.start1 = 100u;
    reference.chrom_id = 0u;

    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u;
    edit.transcript_strand = (int8_t)1;
    edit.first.cdna_anchor1 = 10u; edit.first.exonic = 1u;
    edit.last = edit.first;
    edit.feature_first = edit.first;
    edit.feature_last = edit.first;
    edit.feature_first.genomic_pos1 = 109u;
    edit.feature_last.genomic_pos1 = 109u;
    edit.event.chrom_id = 0u;
    edit.event.kind = (uint8_t)DUCKVEP_KIND_INDEL;

    /* VEP calls hgvs_variant_notation() before _clip_alleles(). The clamped
     * reference is C, so terminal CG>CC is directly C>CC and remains a
     * printable duplication even though the differing genomic base is past
     * the transcript endpoint. */
    edit.ref = genomic_reference + 10u;
    edit.ref_length = 1u;
    edit.alt = alt_cc + 1u;
    edit.alt_length = 1u;
    edit.feature_ref = genomic_reference + 9u;
    edit.feature_ref_length = 2u;
    edit.feature_alt = alt_cc;
    edit.feature_alt_length = 2u;
    edit.event.start1 = 110u;
    edit.event.end1 = 110u;
    edit.event.feature_start1 = 109u;
    edit.event.feature_end1 = 110u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                  &s.tx, &s.ex, NULL, &reference, &edit, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_DNA_DUPLICATION, fact.shape);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.*3dup", rendered));

    /* Only type 'dup' bypasses TranscriptVariationAllele::_clip_alleles.
     * Three copies are clipped to an insertion past the transcript end;
     * executable VEP 116 emits no HGVSc (seed-27182818 CGT>CCC witness). */
    edit.alt = alt_ccc + 1u;
    edit.alt_length = 2u;
    edit.feature_alt = alt_ccc;
    edit.feature_alt_length = 3u;
    ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE,
              duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                  &s.tx, &s.ex, NULL, &reference, &edit, &fact));
    ASSERT_EQ(0u, fact.shape);

    /* ACG>ACC and CGT>CAC become insertions only after clipping. Their
     * insertion coordinates fall beyond this transcript, so VEP emits no
     * HGVSc. Copying sequence at the endpoint must not promote either one to
     * a duplication after that failed projection. */
    edit.ref = genomic_reference + 10u;
    edit.ref_length = 1u;
    edit.alt = alt_acc + 2u;
    edit.alt_length = 1u;
    edit.feature_ref = genomic_reference + 8u;
    edit.feature_ref_length = 3u;
    edit.feature_alt = alt_acc;
    edit.feature_alt_length = 3u;
    edit.feature_first.cdna_anchor1 = 9u;
    edit.feature_first.genomic_pos1 = 108u;
    edit.feature_last.cdna_anchor1 = 10u;
    edit.feature_last.genomic_pos1 = 109u;
    edit.event.start1 = 110u;
    edit.event.end1 = 110u;
    edit.event.feature_start1 = 108u;
    edit.event.feature_end1 = 110u;
    ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE,
              duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                  &s.tx, &s.ex, NULL, &reference, &edit, &fact));

    edit.ref = genomic_reference + 10u;
    edit.ref_length = 2u;
    edit.alt = alt_cac + 1u;
    edit.alt_length = 2u;
    edit.feature_ref = genomic_reference + 9u;
    edit.feature_ref_length = 3u;
    edit.feature_alt = alt_cac;
    edit.feature_alt_length = 3u;
    edit.feature_first.cdna_anchor1 = 10u;
    edit.feature_first.genomic_pos1 = 109u;
    edit.feature_last = edit.feature_first;
    edit.event.start1 = 110u;
    edit.event.end1 = 111u;
    edit.event.feature_start1 = 109u;
    edit.event.feature_end1 = 111u;
    ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE,
              duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                  &s.tx, &s.ex, NULL, &reference, &edit, &fact));
    PASS();
}

TEST hgvs_clamped_multiplication_projects_clipped_ends_on_both_strands(void) {
    /* Source-derived VEP-116 rule: type 'dup' alone skips allele clipping.
     * Cross every nucleotide, strand, and copy count without altering the
     * random generator that discovered the terminal CGT>CCC counterexample. */
    for (unsigned base = 0u; base < 4u; base++) {
        for (int strand = -1; strand <= 1; strand += 2) {
            for (uint16_t copies = 2u; copies <= 16u; copies++) {
                struct kprop_proj_scene s = {0};
                uint8_t genome[40], alleles[32];
                uint16_t chrom = 0u;
                uint32_t pos = strand > 0 ? 103u : 101u - copies;
                uint32_t end = pos + copies - 1u, ref_offset = 0u, alt_offset = 16u;
                uint8_t kind = copies == 2u ? DUCKVEP_KIND_SNV : DUCKVEP_KIND_MNV;
                duckvep_variant_batch_t batch = {0};
                duckvep_transcript_edit_t edit;
                duckvep_hgvs_dna_fact_t fact;
                duckvep_hgvs_reference_window_t reference = {0};
                memset(genome, 'A', sizeof genome);
                memset(alleles, "ACGT"[base], sizeof alleles);
                size_t terminal = strand > 0 ? 20u : 17u;
                genome[terminal] = "ACGT"[base];
                for (size_t j = 1u; j < copies; j++)
                    genome[strand > 0 ? terminal + j : terminal - j] = "ACGT"[(base + 1u) % 4u];
                s.chrom = chrom; s.tstart = 100u; s.tend = 103u; s.strand = (int8_t)strand;
                s.excnt = 1u; s.es[0] = 100u; s.ee[0] = 103u;
                s.cs[0] = 1u; s.ce[0] = 4u; s.phase[0] = -1;
                kprop_proj_scene_finish(&s);
                batch.count = 1u; batch.chrom_id = &chrom; batch.pos1 = &pos;
                memcpy(alleles, genome + pos - 83u, copies);
                batch.variant_kind = &kind; batch.end1 = &end;
                batch.ref_offset = &ref_offset; batch.alt_offset = &alt_offset;
                batch.ref_length = &copies; batch.alt_length = &copies;
                batch.allele_bytes = alleles; batch.allele_bytes_len = 16u + copies;
                reference.bases = genome; reference.length = sizeof genome;
                reference.start1 = 83u; reference.chrom_id = chrom;
                duckvep_transcript_edit_status_t projected = duckvep_transcript_edit_build(
                    &s.tx, &s.ex, NULL, &batch, 0u, 0u, NULL, 0u, &edit);
                if (projected != DUCKVEP_TRANSCRIPT_EDIT_OK)
                    fprintf(stderr, "multiplication base=%u strand=%d copies=%u projection=%d\n",
                        base, strand, copies, projected);
                ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK, projected);
                memset(&fact, 0xa5, sizeof fact);
                duckvep_hgvs_status_t status = duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
                    &s.tx, &s.ex, NULL, &reference, &edit, &fact);
                if (copies == 2u) {
                    ASSERT_EQ(DUCKVEP_HGVS_OK, status);
                    ASSERT_EQ(DUCKVEP_HGVS_DNA_DUPLICATION, fact.shape);
                    ASSERT_EQ(4, fact.first.base);
                    ASSERT_EQ(4, fact.last.base);
                } else {
                    duckvep_hgvs_dna_fact_t empty = {0};
                    ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE, status);
                    ASSERT_EQ(0, memcmp(&empty, &fact, sizeof fact));
                }
            }
        }
    }
    PASS();
}

TEST hgvs_true_feature_inversion_is_not_delins(void) {
    static const uint8_t ref[] = {'A', 'C'};
    static const uint8_t alt[] = {'G', 'T'};
    struct kprop_proj_scene s;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t fact;
    char rendered[64];
    size_t required;

    memset(&s, 0, sizeof s);
    s.chrom = 0u; s.tstart = 100u; s.tend = 109u; s.strand = (int8_t)1;
    s.excnt = 1u; s.cds_s = 100u; s.cds_e = 109u;
    s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u; s.ce[0] = 10u;
    s.phase[0] = 0;
    kprop_proj_scene_finish(&s);
    memset(&edit, 0, sizeof edit);
    edit.tx_idx = 0u; edit.transcript_strand = (int8_t)1;
    edit.first.cdna_anchor1 = 1u; edit.first.exonic = 1u;
    edit.last.cdna_anchor1 = 2u; edit.last.exonic = 1u;
    edit.feature_first = edit.first; edit.feature_last = edit.last;
    edit.ref = ref; edit.ref_length = 2u;
    edit.alt = alt; edit.alt_length = 2u;
    edit.feature_ref = ref; edit.feature_ref_length = 2u;
    edit.feature_alt = alt; edit.feature_alt_length = 2u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_fact_build_genomic_shifted(
                  &s.tx, &s.ex, NULL, &edit, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_DNA_INVERSION, fact.shape);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_dna_render_basic(
                  &fact, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("c.1_2inv", rendered));
    PASS();
}

static int kprop_hgvs_protein_render_scene(
    const uint8_t *ref_cds,
    size_t         ref_cds_length,
    uint32_t       cds_start,
    const uint8_t *ref,
    uint32_t       ref_length,
    const uint8_t *alt,
    uint32_t       alt_length,
    const uint8_t *post_cds,
    size_t         post_cds_length,
    int            predicted,
    char          *rendered,
    size_t         rendered_capacity,
    uint8_t       *shape_out,
    size_t        *required_out) {

    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t fact;
    uint8_t alt_cds[128];
    uint8_t ref_peptide[64];
    uint8_t alt_peptide[64];

    memset(&edit, 0, sizeof edit);
    edit.cds_start = cds_start;
    edit.ref = ref;
    edit.ref_len = ref_length;
    edit.alt = alt;
    edit.alt_len = alt_length;
    edit.variant_strand = (int8_t)1;
    edit_set.edits = &edit;
    edit_set.count = 1u;
    if (duckvep_coding_context_build(
            ref_cds, ref_cds_length, &edit_set, (int8_t)1,
            DUCKVEP_CODON_TABLE_STANDARD,
            alt_cds, sizeof alt_cds,
            ref_peptide, sizeof ref_peptide,
            alt_peptide, sizeof alt_peptide, &context) !=
            DUCKVEP_CODING_CONTEXT_OK) {
        return 0;
    }
    context.post_cds_bases = post_cds;
    context.post_cds_length = post_cds_length;
    context.post_cds_complete = 1u;
    if (duckvep_coding_context_delta_fill(&context, 0u, &delta) !=
            DUCKVEP_CONTEXT_DELTA_OK || delta.valid == 0u ||
        duckvep_hgvs_protein_fact_build(&context, &delta, &fact) !=
            DUCKVEP_HGVS_OK ||
        duckvep_hgvs_protein_render(
            &fact, predicted, rendered, rendered_capacity, required_out) !=
            DUCKVEP_HGVS_OK) {
        return 0;
    }
    if (shape_out != NULL) *shape_out = fact.shape;
    return 1;
}

TEST hgvs_single_residue_sidecar_matches_core_shapes(void) {
    duckvep_hgvs_protein_fact_t fact;
    char rendered[64];
    char small[8];
    size_t required = 0u;
    uint32_t valid =
        (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID;

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build_single_residue(
                  2u, (uint8_t)'E', (uint8_t)'D', valid,
                  DUCKVEP_COMPAT_VEP_116, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL,
              duckvep_hgvs_protein_render(
                  &fact, 0, small, sizeof small, &required));
    ASSERT_EQ(strlen("p.Glu2Asp"), required);
    ASSERT_EQ(0, strcmp("p.Glu2A", small));
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Glu2Asp", rendered));

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build_single_residue(
                  2u, (uint8_t)'E', (uint8_t)'E', valid,
                  DUCKVEP_COMPAT_VEP_116, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 1, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.(Glu2=)", rendered));

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build_single_residue(
                  2u, (uint8_t)'E', (uint8_t)'*', valid,
                  DUCKVEP_COMPAT_VEP_116, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Glu2Ter", rendered));

    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build_single_residue(
                  1u, (uint8_t)'M', (uint8_t)'V',
                  valid | (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_START_LOST,
                  DUCKVEP_COMPAT_VEP_116, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Met1?", rendered));

    ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE,
              duckvep_hgvs_protein_fact_build_single_residue(
                  4u, (uint8_t)'*', (uint8_t)'Q',
                  valid | (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_STOP_LOST,
                  DUCKVEP_COMPAT_VEP_116, &fact));
    PASS();
}

TEST hgvs_sidecar_requires_frameshift_proof_for_length_change(void) {
    duckvep_coding_context_t context;
    uint32_t valid =
        (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID;

    memset(&context, 0, sizeof context);
    ASSERT(duckvep_sequence_delta_consequence_flags_complete_for_hgvs(
        &context, valid));

    context.length_diff = -2;
    ASSERT_FALSE(duckvep_sequence_delta_consequence_flags_complete_for_hgvs(
        &context, valid));
    ASSERT(duckvep_sequence_delta_consequence_flags_complete_for_hgvs(
        &context,
        valid | (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT));
    ASSERT_FALSE(duckvep_sequence_delta_consequence_flags_complete_for_hgvs(
        &context, 0u));
    PASS();
}

TEST compatibility_policy_inventory_is_versioned(void) {
    const uint32_t all_vep116_language_leaks =
        (uint32_t)(DUCKVEP_COMPAT_HGVS_INCOMPLETE_CODON_ASSIGNMENT |
                   DUCKVEP_COMPAT_HGVS_ALTERNATE_CDS_STANDARD_TABLE |
                   DUCKVEP_COMPAT_HGVS_TERMINAL_PARTIAL_INSERTION |
                   DUCKVEP_COMPAT_HGVS_NEGATIVE_SUBSTR |
                   DUCKVEP_COMPAT_HGVS_XAA_AS_TER);
    duckvep_compat_policy_t vep116 =
        duckvep_compat_policy(DUCKVEP_COMPAT_VEP_116);
    duckvep_compat_policy_t strict =
        duckvep_compat_policy(DUCKVEP_COMPAT_STRICT);

    ASSERT_EQ(all_vep116_language_leaks, vep116.flags);
    ASSERT_EQ(0u, strict.flags);
    ASSERT_EQ(DUCKVEP_CODON_TABLE_STANDARD,
              duckvep_compat_hgvs_alternate_codon_table(
                  DUCKVEP_COMPAT_VEP_116,
                  DUCKVEP_CODON_TABLE_VERT_MITO));
    ASSERT_EQ(DUCKVEP_CODON_TABLE_VERT_MITO,
              duckvep_compat_hgvs_alternate_codon_table(
                  DUCKVEP_COMPAT_STRICT,
                  DUCKVEP_CODON_TABLE_VERT_MITO));
    PASS();
}

TEST compatibility_policy_rejects_unknown_profiles_everywhere(void) {
    static const uint8_t ref_cds[] = {'A', 'T', 'G'};
    static const uint8_t ref_peptide[] = {'M'};
    static const uint8_t alt_peptide[] = {'V'};
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t fact;
    uint32_t valid =
        (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID;
    char rendered[32];
    size_t required = 0u;

    memset(&context, 0, sizeof context);
    memset(&delta, 0, sizeof delta);
    memset(&fact, 0, sizeof fact);
    context.ref_cds = ref_cds;
    context.ref_cds_len = sizeof ref_cds;
    context.ref_peptide = ref_peptide;
    context.ref_peptide_len = sizeof ref_peptide;
    context.alt_peptide = alt_peptide;
    context.alt_peptide_len = sizeof alt_peptide;
    delta.valid = 1u;

    /* The same complete context succeeds with a valid profile, proving that
     * INVALID_ARG below is profile rejection rather than incidental malformed
     * peptide state. */
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build(
                  &context, &delta, &fact));
    context.compatibility_profile = UINT8_MAX;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_protein_fact_build(
                  &context, &delta, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_protein_fact_build_single_residue(
                  1u, (uint8_t)'M', (uint8_t)'V', valid,
                  (duckvep_compat_profile_t)UINT8_MAX, &fact));

    memset(&context, 0, sizeof context);
    memset(&fact, 0, sizeof fact);
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
    fact.first_position1 = 2u;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    context.compatibility_profile = UINT8_MAX;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_protein_frameshift_termination_replay(
                  &context, &fact));
    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_protein_frameshift_termination_replay(
                  &context, &fact));

    memset(&fact, 0, sizeof fact);
    fact.compatibility_profile = UINT8_MAX;
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION;
    fact.first_position1 = 1u;
    fact.last_position1 = 1u;
    fact.reference_first = (uint8_t)'M';
    fact.alternate_first = (uint8_t)'V';
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));

    memset(&context, 0, sizeof context);
    memset(&fact, 0, sizeof fact);
    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    fact.context = &context;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION;
    fact.first_position1 = 1u;
    fact.last_position1 = 1u;
    fact.reference_first = (uint8_t)'M';
    fact.alternate_first = (uint8_t)'V';
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    PASS();
}

TEST strict_compatibility_keeps_xaa_distinct_from_termination(void) {
    static const uint8_t ref_cds[] = {'T', 'A', 'A'};
    static const uint8_t alt_cds[] = {'N', 'N', 'N'};
    static const uint8_t ref_peptide[] = {'*'};
    static const uint8_t alt_peptide[] = {'X', 'A'};
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t fact;
    char rendered[64];
    size_t required = 0u;

    memset(&context, 0, sizeof context);
    memset(&delta, 0, sizeof delta);
    context.ref_cds = ref_cds;
    context.ref_cds_len = sizeof ref_cds;
    context.alt_cds = alt_cds;
    context.alt_cds_len = sizeof alt_cds;
    context.ref_peptide = ref_peptide;
    context.ref_peptide_len = sizeof ref_peptide;
    context.alt_peptide = alt_peptide;
    context.alt_peptide_len = 1u;
    context.has_single_edit = 1u;
    context.single_edit_cds_start = 1u;
    context.single_edit_ref_len = 3u;
    context.single_edit_alt_len = 3u;
    delta.valid = 1u;

    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build(&context, &delta, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_EQUAL, fact.shape);

    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build(&context, &delta, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_SUBSTITUTION, fact.shape);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Ter1Xaa", rendered));

    memset(&fact, 0, sizeof fact);
    context.ref_peptide_len = sizeof ref_peptide;
    context.alt_peptide_len = sizeof alt_peptide;
    fact.context = &context;
    fact.window.ref_whole_length = sizeof ref_peptide;
    fact.window.ref_length = sizeof ref_peptide;
    fact.window.alt_whole_length = sizeof alt_peptide;
    fact.window.alt_length = sizeof alt_peptide;
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_DELINS;
    fact.first_position1 = 1u;
    fact.last_position1 = 1u;
    fact.ref_length = sizeof ref_peptide;
    fact.alt_length = sizeof alt_peptide;
    fact.reference_first = (uint8_t)'*';

    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Ter1delinsTer", rendered));

    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Ter1delinsXaaAla", rendered));
    PASS();
}

TEST hgvs_protein_facts_render_core_vep_shapes(void) {
    static const uint8_t cds[] = {
        'A','T','G', 'G','A','A', 'T','T','T', 'T','A','A'
    };
    static const uint8_t short_cds[] = {
        'A','T','G', 'G','A','A', 'T','A','A'
    };
    static const uint8_t terminal_cds[] = {
        'T','G','G', 'T','A','A'
    };
    static const uint8_t terminal_delins[] = {
        'T','C','A', 'C','G','T', 'C','G','T', 'T','A','A'
    };
    static const uint8_t terminal_trp_cds[] = {
        'T','A','T', 'T','G','G', 'T','A','A'
    };
    static const uint8_t a[] = {'A'};
    static const uint8_t c[] = {'C'};
    static const uint8_t g[] = {'G'};
    static const uint8_t t[] = {'T'};
    static const uint8_t gaa[] = {'G','A','A'};
    static const uint8_t gcc[] = {'G','C','C'};
    static const uint8_t gaa_ttt[] = {'G','A','A','T','T','T'};
    static const uint8_t post_stop[] = {'T','A','A'};
    static const uint8_t post_gly_stop[] = {
        'G','G','G', 'T','A','A'
    };
    char rendered[96];
    char small[8];
    size_t required = 0u;
    uint8_t shape = 0u;

    /* These complete fact-builder scenes are minimized from executable VEP
     * 116 differential rows in
     * conformance/data/hgvs_compatibility_witnesses.tsv, extracted from the
     * complete state_exploration_seed_31415927 HGVS differential:
     * chrDuck:119:G:GATA -> p.Trp0_Trp1insIle and
     * chrDuck:120:A:ACCA -> p.Trp1_?0 on DUCK1-201.  They exercise the
     * negative-substr policy where the fact is owned, rather than fabricating
     * position-zero facts only for the renderer. */
    {
        static const uint8_t ref_cds[] = {'T', 'G', 'G', 'T', 'A', 'A'};
        static const uint8_t ref_peptide[] = {'W', '*'};
        static const uint8_t inserted_peptide[] = {'I', 'W', '*'};
        static const uint8_t start_lost_peptide[] = {'T', 'W', '*'};
        duckvep_coding_context_t context;
        duckvep_sequence_delta_t delta;
        duckvep_hgvs_protein_fact_t fact;

        memset(&context, 0, sizeof context);
        memset(&delta, 0, sizeof delta);
        context.ref_cds = ref_cds;
        context.ref_cds_len = sizeof ref_cds;
        context.ref_peptide = ref_peptide;
        context.ref_peptide_len = sizeof ref_peptide;
        context.alt_peptide = inserted_peptide;
        context.alt_peptide_len = sizeof inserted_peptide;
        delta.valid = 1u;
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_fact_build(
                      &context, &delta, &fact));
        ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_INSERTION, fact.shape);
        ASSERT_EQ(0u, fact.first_position1);
        ASSERT_EQ(1u, fact.last_position1);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_render(
                      &fact, 0, rendered, sizeof rendered, &required));
        ASSERT_EQ(0, strcmp("p.Trp0_Trp1insIle", rendered));

        context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
        ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE,
                  duckvep_hgvs_protein_fact_build(
                      &context, &delta, &fact));

        context.alt_peptide = start_lost_peptide;
        context.alt_peptide_len = sizeof start_lost_peptide;
        context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
        delta.start_lost = 1u;
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_fact_build(
                      &context, &delta, &fact));
        ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_START_LOST, fact.shape);
        ASSERT_EQ(1u, fact.first_position1);
        ASSERT_EQ(0u, fact.last_position1);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_render(
                      &fact, 0, rendered, sizeof rendered, &required));
        ASSERT_EQ(0, strcmp("p.Trp1_?0", rendered));
        context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
        ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE,
                  duckvep_hgvs_protein_fact_build(
                      &context, &delta, &fact));
    }

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 6u, a, 1u, c, 1u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_SUBSTITUTION, shape);
    ASSERT_EQ(0, strcmp("p.Glu2Asp", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 6u, a, 1u, g, 1u, NULL, 0u, 1,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_EQUAL, shape);
    ASSERT_EQ(0, strcmp("p.(Glu2=)", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 4u, g, 1u, t, 1u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_SUBSTITUTION, shape);
    ASSERT_EQ(0, strcmp("p.Glu2Ter", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 1u, a, 1u, g, 1u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_START_LOST, shape);
    ASSERT_EQ(0, strcmp("p.Met1?", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 4u, gaa, 3u, NULL, 0u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_DELETION, shape);
    ASSERT_EQ(0, strcmp("p.Glu2del", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 7u, NULL, 0u, gcc, 3u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_INSERTION, shape);
    ASSERT_EQ(0, strcmp("p.Glu2_Phe3insAla", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 7u, NULL, 0u, gaa, 3u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_DUPLICATION, shape);
    ASSERT_EQ(0, strcmp("p.Glu2dup", rendered));

    /* Executable VEP 116 does not shift an inserted Trp across the final
     * translated Trp: _get_surrounding_peptides() returns undef when its
     * post-variant position equals the peptide length.  The nucleotide event
     * can still be c.4_6dup while HGVSp remains an insertion rather than a
     * peptide duplication.  This is the minimized state discovered by the
     * held-out seed-161803399 differential (C>CTGG at chrDuck:234). */
    ASSERT(kprop_hgvs_protein_render_scene(
        terminal_trp_cds, sizeof terminal_trp_cds, 4u,
        NULL, 0u, terminal_trp_cds + 3u, 3u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_INSERTION, shape);
    ASSERT_EQ(0, strcmp("p.Tyr1_Trp2insTrp", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 4u, gaa_ttt, 6u, gcc, 3u, NULL, 0u, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_DELINS, shape);
    ASSERT_EQ(0, strcmp("p.Glu2_Phe3delinsAla", rendered));

    /* VEP 116 caches the complete local reference peptide before clipping.
     * Here a nucleotide delins changes the terminal stop into a longer
     * peptide ending in the same stop. Peptide clipping leaves an insertion,
     * and the cached leading '*' supplies its second presentation flank. */
    ASSERT(kprop_hgvs_protein_render_scene(
        terminal_cds, sizeof terminal_cds, 4u,
        terminal_cds + 3u, 3u,
        terminal_delins, sizeof terminal_delins,
        NULL, 0u, 0, rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_INSERTION, shape);
    ASSERT_EQ(0, strcmp("p.Trp1_Ter2insSerArgArg", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 4u, NULL, 0u, c, 1u,
        post_stop, sizeof post_stop, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_FRAMESHIFT, shape);
    ASSERT_EQ(0, strcmp("p.Glu2ArgfsTer?", rendered));

    ASSERT(kprop_hgvs_protein_render_scene(
        short_cds, sizeof short_cds, 7u, t, 1u, c, 1u,
        post_gly_stop, sizeof post_gly_stop, 0,
        rendered, sizeof rendered, &shape, &required));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_EXTENSION, shape);
    ASSERT_EQ(0, strcmp("p.Ter3GlnextTer2", rendered));

    ASSERT(!kprop_hgvs_protein_render_scene(
        cds, sizeof cds, 6u, a, 1u, c, 1u, NULL, 0u, 0,
        small, sizeof small, &shape, &required));
    ASSERT_EQ(0, strcmp("p.Glu2A", small));
    ASSERT_EQ(strlen("p.Glu2Asp"), required);
    PASS();
}

TEST hgvs_short_alternate_cds_reproduces_vep_trim_assignment(void) {
    /* COVERAGE_WITNESS: HGVSp frameshift coverage/shortened */
    static const uint8_t ref_cds[3] = {'A','T','G'};
    static const uint8_t alternate_cds[2] = {'A','T'};
    static const uint8_t ref_peptide[1] = {'M'};
    static const uint8_t empty_alt_peptide[1] = {0u};
    static const uint8_t post_cds[3] = {'T','A','A'};
    size_t alternate_length;

    for (alternate_length = 1u; alternate_length <= 2u;
         alternate_length++) {
        duckvep_coding_context_t context;
        duckvep_sequence_delta_t delta;
        duckvep_hgvs_protein_fact_t fact;
        char rendered[32];
        size_t required = 0u;

        memset(&context, 0, sizeof context);
        memset(&delta, 0, sizeof delta);
        context.ref_cds = ref_cds;
        context.ref_cds_len = sizeof ref_cds;
        context.alt_cds = alternate_cds;
        context.alt_cds_len = alternate_length;
        context.ref_peptide = ref_peptide;
        context.ref_peptide_len = sizeof ref_peptide;
        context.alt_peptide = empty_alt_peptide;
        context.alt_peptide_len = 0u;
        context.length_diff = (int64_t)alternate_length - 3;
        context.applied_edits = 1u;
        context.has_single_edit = 1u;
        context.cds_changed = 1u;
        context.single_edit_cds_start = (uint32_t)alternate_length + 1u;
        context.single_edit_ref_len = (uint32_t)(3u - alternate_length);
        context.single_edit_alt_len = 0u;
        context.codon_table = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
        context.post_cds_bases = post_cds;
        context.post_cds_length = sizeof post_cds;
        context.post_cds_complete = 1u;
        delta.valid = 1u;
        delta.frameshift = 1u;

        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_fact_build(
                      &context, &delta, &fact));
        ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_SUBSTITUTION, fact.shape);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_render(
                      &fact, 0, rendered, sizeof rendered, &required));
        ASSERT_EQ(0, strcmp("p.Met1Ter", rendered));
    }

    /* For length >= 3, VEP's assignment keeps an incomplete alternate CDS
     * instead of trimming it to complete codons. This changes where the
     * appended transcript tail is read. Strict mode trims the fourth base and
     * therefore sees the TAA immediately; VEP 116 retains it and sees ATA. */
    {
        static const uint8_t alternate_cds_four[4] = {'A','T','G','A'};
        static const uint8_t post_stop[3] = {'T','A','A'};
        duckvep_coding_context_t context;
        duckvep_hgvs_protein_fact_t fact;

        memset(&context, 0, sizeof context);
        context.alt_cds = alternate_cds_four;
        context.alt_cds_len = sizeof alternate_cds_four;
        context.codon_table = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
        context.post_cds_bases = post_stop;
        context.post_cds_length = sizeof post_stop;
        context.post_cds_complete = 1u;

        memset(&fact, 0, sizeof fact);
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
        fact.first_position1 = 2u;
        context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
        fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_frameshift_termination_replay(
                      &context, &fact));
        ASSERT_FALSE(fact.termination_known);

        context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
        fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
        ASSERT_EQ(DUCKVEP_HGVS_OK,
                  duckvep_hgvs_protein_frameshift_termination_replay(
                      &context, &fact));
        ASSERT(fact.termination_known);
        ASSERT_EQ(1u, fact.termination_distance);
    }
    PASS();
}

TEST hgvs_equal_local_peptides_do_not_hide_frameshift(void) {
    /* Pinned original padding witnesses 4/10/14: ACN>AN, NCN>NN and ACC>AC
     * at genomic 14 in complete CDS 11..22 models. In witness 10 the local
     * peptides are X/X, but _get_hgvs_protein_type still selects frameshift;
     * _get_fs_peptides finds the first full-translation difference at Ala3. */
    static const struct {
        const char *cds;
        const char *reference;
        const char *alternate;
        uint8_t local_reference;
        const char *hgvs;
    } cases[] = {
        {"ATGACNGCCTAA", "ACN", "AN", 'T', "p.Thr2Ter"},
        {"ATGNCNGCCTAA", "NCN", "NN", 'X', "p.Ala3ProfsTer?"},
        {"ATGACCGCCTAA", "ACC", "AC", 'T', "p.Ala3ProfsTer?"}
    };
    static const duckvep_compat_profile_t profiles[] = {DUCKVEP_COMPAT_STRICT, DUCKVEP_COMPAT_VEP_116};
    for (size_t i = 0u; i < sizeof cases / sizeof cases[0]; i++) {
        duckvep_event_t event;
        ASSERT(duckvep_event_prepare_small(14u, (const uint8_t *)cases[i].reference, 3u,
            (const uint8_t *)cases[i].alternate, 2u, &event));
        duckvep_haplotype_edit_t edit = {0};
        edit.cds_start = event.feature_start1 - 10u;
        edit.ref = (const uint8_t *)cases[i].reference + event.feature_allele_offset;
        edit.ref_len = event.ref_diff_length;
        edit.alt = (const uint8_t *)cases[i].alternate + event.feature_allele_offset;
        edit.alt_len = event.alt_diff_length;
        edit.variant_strand = 1;
        duckvep_edit_set_t edits = {&edit, 1u};
        uint8_t alt_cds[32], ref_peptide[16], alt_peptide[16];
        duckvep_coding_context_t context;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
            (const uint8_t *)cases[i].cds, 12u, &edits, 1, DUCKVEP_CODON_TABLE_STANDARD,
            alt_cds, sizeof alt_cds, ref_peptide, sizeof ref_peptide,
            alt_peptide, sizeof alt_peptide, &context));
        context.post_cds_complete = 1u;
        duckvep_coding_peptide_window_t window;
        ASSERT(duckvep_coding_context_peptide_window_open(&context, &window));
        ASSERT_EQ(1u, window.ref_length);
        ASSERT_EQ(1u, window.alt_length);
        ASSERT_EQ(cases[i].local_reference,
            duckvep_coding_context_peptide_window_base(&context, &window, 0, 0u));
        ASSERT_EQ('X', duckvep_coding_context_peptide_window_base(&context, &window, 1, 0u));
        duckvep_sequence_delta_t delta;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK, duckvep_coding_context_delta_fill(&context, 0u, &delta));
        ASSERT(delta.valid);
        ASSERT(delta.frameshift);
        ASSERT_FALSE(delta.stop_retained);
        for (size_t p = 0u; p < sizeof profiles / sizeof profiles[0]; p++) {
            context.compatibility_profile = (uint8_t)profiles[p];
            duckvep_hgvs_protein_fact_t fact;
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_fact_build(&context, &delta, &fact));
            ASSERT_EQ(window.ref_length, fact.window.ref_length);
            ASSERT_EQ(window.alt_length, fact.window.alt_length);
            if (i == 1u) {
                ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_FRAMESHIFT, fact.shape);
                ASSERT_EQ(3u, fact.first_position1);
                ASSERT_EQ('A', fact.reference_first);
                ASSERT_EQ('P', fact.alternate_first);
            }
            char rendered[64];
            size_t required;
            ASSERT_EQ(DUCKVEP_HGVS_OK,
                duckvep_hgvs_protein_render(&fact, 0, rendered, sizeof rendered, &required));
            ASSERT_STR_EQ(i == 0u && profiles[p] == DUCKVEP_COMPAT_STRICT
                ? "p.Thr2XaafsTer?" : cases[i].hgvs, rendered);
        }
    }
    PASS();
}

TEST hgvs_frameshift_xaa_immediate_stop_is_compatibility_policy(void) {
    /* Original VEP event 32261: ATGNNAGCCTAA, genomic 13 G>GAC.
     * Its actual HGVS shift produces ATGNCANAGCCTAA; _get_fs_peptides
     * reports reference A, alternate X, position 3. VEP converts Xaa to Ter
     * before choosing immediate-stop formatting. Strict mode keeps Xaa. */
    static const uint8_t cds[] = "ATGNNAGCCTAA";
    static const uint8_t shifted_insertion[] = "CA";
    duckvep_haplotype_edit_t edit = {0};
    duckvep_edit_set_t edits = {&edit, 1u};
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    uint8_t alt_cds[32], ref_peptide[16], alt_peptide[16];

    edit.cds_start = 5u;
    edit.alt = shifted_insertion;
    edit.alt_len = 2u;
    edit.variant_strand = 1;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        cds, sizeof cds - 1u, &edits, 1, DUCKVEP_CODON_TABLE_STANDARD,
        alt_cds, sizeof alt_cds, ref_peptide, sizeof ref_peptide,
        alt_peptide, sizeof alt_peptide, &context));
    ASSERT_EQ(14u, context.alt_cds_len);
    ASSERT_MEM_EQ("ATGNCANAGCCTAA", context.alt_cds, 14u);
    ASSERT_EQ('A', duckvep_coding_context_peptide_base(&context, 0, 2u));
    ASSERT_EQ('X', duckvep_coding_context_peptide_base(&context, 1, 2u));
    context.post_cds_complete = 1u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
        duckvep_coding_context_delta_fill(&context, 0u, &delta));
    ASSERT(delta.valid);
    ASSERT(delta.frameshift);
    ASSERT_FALSE(delta.stop_gained);
    for (size_t profile = 0u; profile < 2u; profile++) {
        duckvep_hgvs_protein_fact_t fact;
        char rendered[64];
        size_t required = 0u;
        context.compatibility_profile = (uint8_t)(profile == 0u
            ? DUCKVEP_COMPAT_VEP_116 : DUCKVEP_COMPAT_STRICT);
        ASSERT_EQ(DUCKVEP_HGVS_OK,
            duckvep_hgvs_protein_fact_build(&context, &delta, &fact));
        ASSERT_EQ(profile == 0u ? DUCKVEP_HGVS_PROTEIN_SUBSTITUTION
            : DUCKVEP_HGVS_PROTEIN_FRAMESHIFT, fact.shape);
        ASSERT_EQ('X', fact.alternate_first);
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
            &fact, 0, rendered, sizeof rendered, &required));
        ASSERT_STR_EQ(profile == 0u ? "p.Ala3Ter" : "p.Ala3XaafsTer?", rendered);
    }
    PASS();
}

TEST hgvs_alternate_cds_preserves_reference_codon_table(void) {
    /* Original records from ambiguous_indel_consensus, events
     * 4/17/18/73/74/76/284/353/354/356, retained in
     * conformance/data/indel_translation_witnesses.jsonl.gz.
     * VEP 116 _get_fs_peptides()
     * translates alternate CDS with BioPerl's default table, while _peptide
     * and the local consequence predicates retain the transcript table.
     * Strict expectations below test our declared-table policy, not VEP. */
    static const uint8_t cds[] = "ATGAAAGCCTAA";
    static const struct {
        uint8_t table;
        uint32_t position;
        const char *reference;
        const char *alternate;
        const char *vep_hgvsp;
        const char *strict_hgvsp;
        uint8_t reference_second;
        uint8_t alternate_second;
        uint8_t frameshift;
        uint8_t stop_gained;
    } cases[] = {
        {1u, 13u, "G", "GT", "p.Lys2Ter", "p.Lys2Ter", 'K', '*', 1u, 0u},
        {1u, 14u, "A", "AG", "p.Lys2ArgfsTer?", "p.Lys2ArgfsTer?", 'K', 'R', 1u, 0u},
        {1u, 14u, "A", "AT", "p.Lys2IlefsTer?", "p.Lys2IlefsTer?", 'K', 'I', 1u, 0u},
        {2u, 14u, "A", "AG", "p.Lys2ArgfsTer?", "p.Lys2Ter", 'K', '*', 1u, 1u},
        {2u, 14u, "A", "AT", "p.Lys2IlefsTer?", "p.Lys2MetfsTer?", 'K', 'M', 1u, 0u},
        {2u, 14u, "A", "AGCC", "p.Lys2delinsSerGln", "p.Lys2delinsSerGln", 'K', 'S', 0u, 0u},
        {6u, 13u, "G", "GT", "p.Lys2Ter", "p.Lys2GlnfsTer?", 'K', 'Q', 1u, 0u},
        {9u, 14u, "A", "AG", "p.Asn2ArgfsTer?", "p.Asn2SerfsTer?", 'N', 'S', 1u, 0u},
        {9u, 14u, "A", "AT", "p.Asn2IlefsTer?", "p.Asn2IlefsTer?", 'N', 'I', 1u, 0u},
        {9u, 14u, "A", "AGCC", "p.Asn2delinsSerGln", "p.Asn2delinsSerGln", 'N', 'S', 0u, 0u}
    };
    struct kprop_proj_scene scene = {0};
    duckvep_sequence_pool_t sequences = {0};
    uint64_t offset = 0u;
    uint32_t length = sizeof cds - 1u;
    uint32_t empty = 0u;
    uint8_t table = 1u;

    scene.tstart = scene.cds_s = scene.es[0] = 11u;
    scene.tend = scene.cds_e = scene.ee[0] = 22u;
    scene.strand = 1;
    scene.excnt = 1u;
    scene.cs[0] = 1u;
    scene.ce[0] = length;
    scene.flags = DUCKVEP_TX_HAS_TRANSLATION | DUCKVEP_TX_BIOTYPE_PROTEIN_CODING;
    kprop_proj_scene_finish(&scene);
    sequences.cds_bytes = cds;
    sequences.cds_bytes_len = length;
    sequences.cds_offset = &offset;
    sequences.cds_length = &length;
    sequences.codon_table = &table;
    sequences.transcript_count = 1u;
    sequences.pre_cds_offset = sequences.post_cds_offset = &offset;
    sequences.pre_cds_length = sequences.post_cds_length = &empty;
    sequences.flanks_complete = 1u;

    for (size_t i = 0u; i < sizeof cases / sizeof cases[0]; i++) {
        duckvep_event_t event;
        duckvep_transcript_edit_t edit;
        duckvep_haplotype_edit_t edits[4];
        duckvep_coding_context_t context;
        duckvep_sequence_delta_t delta;
        uint8_t alt_cds[32], ref_peptide[16], alt_peptide[16];
        duckvep_variant_batch_t variants = {0};
        uint8_t alleles[8];
        uint32_t ref_offset = 0u, alt_offset = 1u;
        uint16_t ref_length = 1u;
        uint16_t alt_length = (uint16_t)strlen(cases[i].alternate);
        uint8_t kind = DUCKVEP_KIND_INS;

        table = cases[i].table;
        alleles[0] = (uint8_t)cases[i].reference[0];
        memcpy(alleles + alt_offset, cases[i].alternate, alt_length);
        variants.chrom_id = &scene.chrom;
        variants.pos1 = variants.end1 = &cases[i].position;
        variants.ref_offset = &ref_offset;
        variants.alt_offset = &alt_offset;
        variants.ref_length = &ref_length;
        variants.alt_length = &alt_length;
        variants.variant_kind = &kind;
        variants.allele_bytes = alleles;
        variants.allele_bytes_len = (size_t)ref_length + alt_length;
        variants.count = 1u;
        ASSERT(duckvep_event_prepare_small(cases[i].position,
            (const uint8_t *)cases[i].reference, strlen(cases[i].reference),
            (const uint8_t *)cases[i].alternate, strlen(cases[i].alternate), &event));
        event.chrom_id = 0u;
        ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
            duckvep_transcript_edit_build_prepared(&scene.tx, &scene.ex,
                &sequences, &variants, 0u, 0u, &event, edits, 4u, &edit));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, edit.cds_status);
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
            duckvep_model_coding_context_build(&scene.tx, &scene.ex, &sequences,
                0u, 1, &event, &edit.cds_edits, alt_cds, sizeof alt_cds,
                ref_peptide, sizeof ref_peptide, alt_peptide, sizeof alt_peptide, &context));
        ASSERT_EQ(cases[i].reference_second,
            duckvep_coding_context_peptide_base(&context, 0, 1u));
        ASSERT_EQ(cases[i].alternate_second,
            duckvep_coding_context_peptide_base(&context, 1, 1u));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_delta_fill(&context, scene.flags, &delta));
        ASSERT(delta.valid);
        ASSERT_EQ(cases[i].frameshift, delta.frameshift);
        ASSERT_EQ(cases[i].stop_gained, delta.stop_gained);
        for (size_t profile = 0u; profile < 2u; profile++) {
            duckvep_hgvs_protein_fact_t fact;
            char rendered[64];
            size_t required = 0u;
            context.compatibility_profile = (uint8_t)(profile == 0u
                ? DUCKVEP_COMPAT_VEP_116 : DUCKVEP_COMPAT_STRICT);
            ASSERT_EQ(DUCKVEP_HGVS_OK,
                duckvep_hgvs_protein_fact_build(&context, &delta, &fact));
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                &fact, 0, rendered, sizeof rendered, &required));
            ASSERT_STR_EQ(profile == 0u ? cases[i].vep_hgvsp : cases[i].strict_hgvsp,
                rendered);
        }
    }
    PASS();
}

TEST hgvs_late_stop_search_reproduces_vep_standard_table_and_precedence(void) {
    static const uint8_t mt_alt_cds[12] = {
        'A','T','G', 'T','G','A', 'G','C','C', 'T','A','A'
    };
    static const uint8_t mt_virtual_ref_cds[12] = {
        'T','G','A', 'G','C','C', 'G','C','C', 'T','A','A'
    };
    static const uint8_t mt_virtual_ref_peptide[4] = {'W', 'A', 'A', '*'};
    static const uint8_t virtual_alt_base[1] = {'A'};
    static const uint8_t ref_peptide[2] = {'M', '*'};
    static const uint8_t terminal_ref_cds[3] = {'T','A','A'};
    static const uint8_t terminal_alt_cds[2] = {'T','A'};
    static const uint8_t terminal_ref_peptide[1] = {'*'};
    static const uint8_t empty_alt_peptide[1] = {0u};
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    duckvep_hgvs_protein_fact_t fact;
    char rendered[48];
    size_t required = 0u;

    /* Ordinary mitochondrial translation treats TGA as Trp. VEP's late HGVS
     * termination search nevertheless uses BioPerl's default table 1, finds
     * that TGA immediately, and reports one changed residue rather than three. */
    memset(&context, 0, sizeof context);
    memset(&fact, 0, sizeof fact);
    context.alt_cds = mt_alt_cds;
    context.alt_cds_len = sizeof mt_alt_cds;
    context.ref_peptide = ref_peptide;
    context.ref_peptide_len = sizeof ref_peptide;
    context.codon_table =
        (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO;
    context.post_cds_complete = 1u;
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
    fact.first_position1 = 2u;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_frameshift_termination_replay(
                  &context, &fact));
    ASSERT(fact.termination_known);
    ASSERT_EQ(1u, fact.termination_distance);

    /* Strict internal oracle mode keeps the transcript's table 2 and therefore
     * reaches the terminal TAA instead of treating mitochondrial TGA as stop. */
    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    memset(&fact, 0, sizeof fact);
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
    fact.first_position1 = 2u;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_frameshift_termination_replay(
                  &context, &fact));
    ASSERT(fact.termination_known);
    ASSERT_EQ(3u, fact.termination_distance);

    /* The virtual single-edit scanner normally starts at the edited codon
     * after using a reference first-stop fact. That fact was established with
     * table 2 and cannot hide the upstream TGA seen by VEP's table-1 replay. */
    memset(&context, 0, sizeof context);
    memset(&fact, 0, sizeof fact);
    context.ref_cds = mt_virtual_ref_cds;
    context.ref_cds_len = sizeof mt_virtual_ref_cds;
    context.alt_cds_len = sizeof mt_virtual_ref_cds;
    context.ref_peptide = mt_virtual_ref_peptide;
    context.ref_peptide_len = sizeof mt_virtual_ref_peptide;
    context.single_edit_alt = virtual_alt_base;
    context.virtual_single_edit = 1u;
    context.has_single_edit = 1u;
    context.single_edit_cds_start = 7u;
    context.single_edit_ref_len = 1u;
    context.single_edit_alt_len = 1u;
    context.transcript_strand = 1;
    context.single_edit_variant_strand = 1;
    context.codon_table = (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO;
    context.ref_first_stop_known = 1u;
    context.ref_first_stop_position1 = 4u;
    context.post_cds_complete = 1u;
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
    fact.first_position1 = 3u;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_VEP_116;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_frameshift_termination_replay(
                  &context, &fact));
    ASSERT_FALSE(fact.termination_known);
    ASSERT_EQ(0u, fact.termination_distance);

    context.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    memset(&fact, 0, sizeof fact);
    fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
    fact.first_position1 = 3u;
    fact.compatibility_profile = (uint8_t)DUCKVEP_COMPAT_STRICT;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_frameshift_termination_replay(
                  &context, &fact));
    ASSERT(fact.termination_known);
    ASSERT_EQ(2u, fact.termination_distance);

    /* The formatter checks stop_lost+deletion before frameshift. A terminal
     * peptide removed by a frame-changing edit is therefore a del-extension. */
    memset(&context, 0, sizeof context);
    memset(&delta, 0, sizeof delta);
    context.ref_cds = terminal_ref_cds;
    context.ref_cds_len = sizeof terminal_ref_cds;
    context.alt_cds = terminal_alt_cds;
    context.alt_cds_len = sizeof terminal_alt_cds;
    context.ref_peptide = terminal_ref_peptide;
    context.ref_peptide_len = sizeof terminal_ref_peptide;
    context.alt_peptide = empty_alt_peptide;
    context.alt_peptide_len = 0u;
    context.codon_table = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    context.post_cds_complete = 1u;
    delta.valid = 1u;
    delta.frameshift = 1u;
    delta.stop_lost = 1u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_fact_build(&context, &delta, &fact));
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_EXTENSION, fact.shape);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_protein_render(
                  &fact, 0, rendered, sizeof rendered, &required));
    ASSERT_EQ(0, strcmp("p.Ter1delextTer?", rendered));
    PASS();
}

TEST hgvs_genomic_search_interval_preserves_vep_insertion_geometry(void) {
    duckvep_event_t event;
    uint32_t start1 = 0u;
    uint32_t end1 = 0u;

    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_DEL;
    event.start1 = 5000u;
    event.end1 = 5002u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_genomic_search_interval(
                  &event, 10000u, &start1, &end1));
    ASSERT_EQ(4000u, start1);
    ASSERT_EQ(6002u, end1);

    /* Retained uploaded padding is a separate REF assertion. It widens the
     * lookup fetch but must not alter VEP's exact genomic-shift slice. */
    event.raw_start1 = 2500u;
    event.raw_end1 = 8000u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_genomic_search_interval(
                  &event, 10000u, &start1, &end1));
    ASSERT_EQ(4000u, start1);
    ASSERT_EQ(6002u, end1);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_reference_fetch_interval(
                  &event, 10000u, &start1, &end1));
    ASSERT_EQ(2500u, start1);
    ASSERT_EQ(8000u, end1);

    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_INS;
    event.interbase = 1u;
    event.insertion_boundary0 = 5000u;
    event.alt_diff_length = 1u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_genomic_search_interval(
                  &event, 10000u, &start1, &end1));
    ASSERT_EQ(4001u, start1);
    ASSERT_EQ(6000u, end1);

    /* A copied source may be longer than VEP's 1,000-base shift flank. The
     * lookup interval covers both adjacent candidates without changing the
     * exact slice above. */
    event.alt_diff_length = 1001u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_reference_fetch_interval(
                  &event, 10000u, &start1, &end1));
    ASSERT_EQ(3000u, start1);
    ASSERT_EQ(7001u, end1);
    event.alt_diff_length = 1u;

    /* Once constrained, VEP's first/last 1000-base slices overlap on a short
     * sequence region.  The helper must request that complete region. */
    event.insertion_boundary0 = 124u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_genomic_search_interval(
                  &event, 260u, &start1, &end1));
    ASSERT_EQ(1u, start1);
    ASSERT_EQ(260u, end1);

    event.insertion_boundary0 = 10000u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_genomic_search_interval(
                  &event, 10000u, &start1, &end1));
    ASSERT_EQ(9001u, start1);
    ASSERT_EQ(10000u, end1);

    event.insertion_boundary0 = 10001u;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_PROJECTION,
              duckvep_hgvs_genomic_search_interval(
                  &event, 10000u, &start1, &end1));
    PASS();
}

TEST hgvs_uploaded_reference_validates_vcf_anchor_and_padding(void) {
    static const uint8_t reference_bases[] = {'T', 'A', 'C', 'G'};
    static const uint8_t correct_anchor[] = {'T'};
    static const uint8_t wrong_anchor[] = {'A'};
    static const uint8_t correct_padded_ref[] = {'T', 'A'};
    static const uint8_t wrong_padded_ref[] = {'G', 'A'};
    duckvep_hgvs_reference_window_t reference;
    duckvep_event_t event;

    memset(&reference, 0, sizeof reference);
    reference.bases = reference_bases;
    reference.length = sizeof reference_bases;
    reference.start1 = 124u;
    reference.chrom_id = 1u;

    memset(&event, 0, sizeof event);
    event.chrom_id = 1u;
    event.raw_start1 = 124u;
    event.raw_end1 = 124u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_uploaded_reference_validate(
                  &reference, &event, correct_anchor,
                  sizeof correct_anchor));
    ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH,
              duckvep_hgvs_uploaded_reference_validate(
                  &reference, &event, wrong_anchor, sizeof wrong_anchor));

    event.raw_end1 = 125u;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
              duckvep_hgvs_uploaded_reference_validate(
                  &reference, &event, correct_padded_ref,
                  sizeof correct_padded_ref));
    ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH,
              duckvep_hgvs_uploaded_reference_validate(
                  &reference, &event, wrong_padded_ref,
                  sizeof wrong_padded_ref));

    /* An erased N anchor must match literally. An untrimmed N substitution,
     * other ambiguous symbols, and a missing anchor remain invalid. */
    static const uint8_t n_anchor[] = "N";
    static const uint8_t n_insertion[] = "NGCC";
    ASSERT(duckvep_event_prepare_small(124u, n_anchor, 1u, n_insertion, 4u, &event));
    event.chrom_id = reference.chrom_id;
    ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, n_anchor, 1u));
    reference.bases = n_anchor;
    reference.length = 1u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, n_anchor, 1u));
    ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, correct_anchor, 1u));
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, (const uint8_t *)"R", 1u));
    reference.bases = (const uint8_t *)"R";
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, n_anchor, 1u));
    reference.bases = n_anchor;
    reference.length = 0u;
    ASSERT_EQ(DUCKVEP_HGVS_MISSING_REFERENCE, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, n_anchor, 1u));
    reference.length = 1u;
    ASSERT(duckvep_event_prepare_small(124u, n_anchor, 1u,
        (const uint8_t *)"A", 1u, &event));
    event.chrom_id = reference.chrom_id;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE, duckvep_hgvs_uploaded_reference_validate(
        &reference, &event, n_anchor, 1u));

    /* Original VEP deletion 5400, NAAG>N, and delins 5404, NAAG>NACGT,
     * erase the same matching N before parsed-allele eligibility. The whole
     * original REF still has to match, not merely the removed interval. */
    static const uint8_t deletion_ref[] = "NAAG";
    static const char *retained_alternates[] = {"N", "NA", "NACGT"};
    reference.bases = deletion_ref;
    reference.length = sizeof deletion_ref - 1u;
    for (size_t i = 0u; i < sizeof retained_alternates / sizeof retained_alternates[0]; i++) {
        ASSERT(duckvep_event_prepare_small(124u, deletion_ref, 4u,
            (const uint8_t *)retained_alternates[i],
            (uint16_t)strlen(retained_alternates[i]), &event));
        event.chrom_id = reference.chrom_id;
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, deletion_ref, 4u));
        reference.bases = (const uint8_t *)"AAAG";
        ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, deletion_ref, 4u));
        reference.bases = (const uint8_t *)"NACG";
        ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, deletion_ref, 4u));
        reference.bases = deletion_ref;
    }
    static const struct {
        const char *reference;
        const char *alternate;
    } invalid[] = {{"RAAG", "R"}, {"CN", "C"}, {"NAAG", "NCCG"}};
    for (size_t i = 0u; i < sizeof invalid / sizeof invalid[0]; i++) {
        reference.bases = (const uint8_t *)invalid[i].reference;
        reference.length = strlen(invalid[i].reference);
        ASSERT(duckvep_event_prepare_small(124u, reference.bases,
            (uint16_t)reference.length, (const uint8_t *)invalid[i].alternate,
            (uint16_t)strlen(invalid[i].alternate), &event));
        event.chrom_id = reference.chrom_id;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, reference.bases, reference.length));
    }
    /* Actual VEP 116 Parser.pm::minimise_alleles observations: ACN>ATCN
     * erases CN, not just N. ACN>AGTN keeps C/GT. Both ends may contain N. */
    static const struct {
        const char *reference;
        const char *alternate;
        uint16_t offset;
        uint16_t ref_length;
        uint16_t alt_length;
    } padding[] = {
        {"ACN", "ATCN", 1u, 0u, 1u},
        {"CN", "TCN", 0u, 0u, 1u},
        {"ACN", "AGTN", 1u, 1u, 2u},
        {"ACN", "AN", 1u, 1u, 0u},
        {"N", "NT", 1u, 0u, 1u},
        {"NCN", "NTCN", 1u, 0u, 1u},
        {"NCN", "NGTN", 1u, 1u, 2u},
        {"NCN", "NN", 1u, 1u, 0u}
    };
    for (size_t i = 0u; i < sizeof padding / sizeof padding[0]; i++) {
        reference.bases = (const uint8_t *)padding[i].reference;
        reference.length = strlen(padding[i].reference);
        ASSERT(duckvep_event_prepare_small(124u, reference.bases,
            (uint16_t)reference.length, (const uint8_t *)padding[i].alternate,
            (uint16_t)strlen(padding[i].alternate), &event));
        event.chrom_id = reference.chrom_id;
        ASSERT_EQ(padding[i].offset, event.feature_allele_offset);
        ASSERT_EQ(padding[i].ref_length, event.ref_diff_length);
        ASSERT_EQ(padding[i].alt_length, event.alt_diff_length);
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, reference.bases, reference.length));
        event.feature_allele_offset = (uint16_t)(reference.length + 1u);
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, reference.bases, reference.length));
        event.feature_allele_offset = (uint16_t)reference.length;
        event.ref_diff_length = 1u;
        event.alt_diff_length = 2u;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, reference.bases, reference.length));
    }
    PASS();
}

TEST hgvs_uploaded_reference_padding_enumeration(void) {
    /* Independent construction: shared A/G/N padding surrounds canonical
     * differing payloads. Their disjoint alphabets fix the expected interval
     * without asking native preparation to supply the expectation. Reverse
     * complements exchange prefix/suffix roles, including right anchors. */
    static const char *prefixes[] = {"", "A", "N", "AN", "NA", "NN", "AAN", "NAA"};
    static const char *suffixes[] = {"", "G", "N", "GN", "NG", "NN", "GGN", "NGG"};
    static const struct {
        const char *reference;
        const char *alternate;
    } payloads[] = {{"", "C"}, {"C", ""}, {"C", "TT"}, {"CC", "T"}};
    size_t cases = 0u;
    for (size_t p = 0u; p < sizeof prefixes / sizeof prefixes[0]; p++) {
        for (size_t s = 0u; s < sizeof suffixes / sizeof suffixes[0]; s++) {
            if (!p && !s) continue; /* Original VCF REF and ALT must be nonempty. */
            for (size_t shape = 0u; shape < sizeof payloads / sizeof payloads[0]; shape++) {
                char ref[16], alt[16];
                snprintf(ref, sizeof ref, "%s%s%s", prefixes[p], payloads[shape].reference, suffixes[s]);
                snprintf(alt, sizeof alt, "%s%s%s", prefixes[p], payloads[shape].alternate, suffixes[s]);
                size_t ref_length = strlen(ref), alt_length = strlen(alt);
                for (int reverse = 0; reverse < 2; reverse++) {
                    uint8_t uploaded[16], alternate[16], genome[16];
                    for (size_t i = 0u; i < ref_length; i++) {
                        uploaded[i] = reverse ? (uint8_t)kprop_complement_base(ref[ref_length - 1u - i])
                            : (uint8_t)ref[i];
                    }
                    for (size_t i = 0u; i < alt_length; i++) {
                        alternate[i] = reverse ? (uint8_t)kprop_complement_base(alt[alt_length - 1u - i])
                            : (uint8_t)alt[i];
                    }
                    memcpy(genome, uploaded, ref_length);
                    duckvep_hgvs_reference_window_t reference = {genome, ref_length, 124u, 1u};
                    duckvep_event_t event;
                    ASSERT(duckvep_event_prepare_small(124u, uploaded, (uint16_t)ref_length,
                        alternate, (uint16_t)alt_length, &event));
                    event.chrom_id = reference.chrom_id;
                    ASSERT_EQ(strlen(reverse ? suffixes[s] : prefixes[p]), event.feature_allele_offset);
                    ASSERT_EQ(strlen(payloads[shape].reference), event.ref_diff_length);
                    ASSERT_EQ(strlen(payloads[shape].alternate), event.alt_diff_length);
                    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_uploaded_reference_validate(
                        &reference, &event, uploaded, ref_length));
                    for (size_t i = 0u; i < ref_length; i++) {
                        genome[i] = uploaded[i] == (uint8_t)'A' ? (uint8_t)'C' : (uint8_t)'A';
                        ASSERT_EQ(DUCKVEP_HGVS_REFERENCE_MISMATCH, duckvep_hgvs_uploaded_reference_validate(
                            &reference, &event, uploaded, ref_length));
                        genome[i] = uploaded[i];
                    }
                    /* The same unminimized replacement has no erased padding;
                     * it cannot borrow the independent-event N permission. */
                    ASSERT(duckvep_event_prepare_replacement(124u, uploaded, (uint16_t)ref_length,
                        alternate, (uint16_t)alt_length, &event));
                    event.chrom_id = reference.chrom_id;
                    ASSERT_EQ(memchr(uploaded, 'N', ref_length) ? DUCKVEP_HGVS_INVALID_ALLELE : DUCKVEP_HGVS_OK,
                        duckvep_hgvs_uploaded_reference_validate(&reference, &event, uploaded, ref_length));
                    cases++;
                }
            }
        }
    }
    ASSERT_EQ(504u, cases);
    /* N inside the differing interval, equal-length shared N padding, and
     * other ambiguous symbols are not licensed by prefix/suffix minimization. */
    static const struct {
        const char *reference;
        const char *alternate;
    } invalid[] = {{"ANC", "AC"}, {"ACN", "ATN"}, {"NCN", "NTN"},
        {"ACR", "ATCR"}, {"RCA", "RTA"}};
    for (size_t i = 0u; i < sizeof invalid / sizeof invalid[0]; i++) {
        const uint8_t *ref = (const uint8_t *)invalid[i].reference;
        duckvep_hgvs_reference_window_t reference = {ref, strlen(invalid[i].reference), 124u, 1u};
        duckvep_event_t event;
        ASSERT(duckvep_event_prepare_small(124u, ref, (uint16_t)reference.length,
            (const uint8_t *)invalid[i].alternate, (uint16_t)strlen(invalid[i].alternate), &event));
        event.chrom_id = reference.chrom_id;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ALLELE, duckvep_hgvs_uploaded_reference_validate(
            &reference, &event, ref, reference.length));
    }
    PASS();
}

TEST hgvs_shift_stops_at_vep_1000_base_cap(void) {
    /* COVERAGE_WITNESS: HGVS shift coverage/at_vep_limit */
    uint8_t reference_bytes[2001];
    static const uint8_t allele = (uint8_t)'A';
    uint16_t chrom = 0u;
    uint32_t tx_start = 1000u;
    uint32_t tx_end = 4000u;
    int8_t strand = 1;
    uint64_t tx_flags = 0u;
    uint32_t exon_offset = 0u;
    uint16_t exon_count = 1u;
    uint32_t cds_start = 1000u;
    uint32_t cds_end = 4000u;
    uint32_t exon_start = 1000u;
    uint32_t exon_end = 4000u;
    uint32_t cdna_start = 1u;
    uint32_t cdna_end = 3001u;
    int8_t phase = 0;
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_hgvs_reference_window_t reference;
    duckvep_transcript_edit_t edit;
    duckvep_transcript_coordinate_t low;
    duckvep_transcript_coordinate_t high;
    duckvep_hgvs_dna_fact_t fact;
    uint32_t fetch_start;
    uint32_t fetch_end;
    size_t shape;
    size_t direction;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    tx.chrom_id = &chrom;
    tx.start1 = &tx_start;
    tx.end1 = &tx_end;
    tx.strand = &strand;
    tx.flags = &tx_flags;
    tx.exon_offset = &exon_offset;
    tx.exon_count = &exon_count;
    tx.cds_start1 = &cds_start;
    tx.cds_end1 = &cds_end;
    tx.transcript_count = 1u;
    exons.start1 = &exon_start;
    exons.end1 = &exon_end;
    exons.cdna_start1 = &cdna_start;
    exons.cdna_end1 = &cdna_end;
    exons.phase = &phase;
    exons.end_phase = &phase;
    exons.exon_count = 1u;

    for (shape = 0u; shape < 2u; shape++) {
        for (direction = 0u; direction < 2u; direction++) {
            strand = direction == 0u ? (int8_t)1 : (int8_t)-1;
            memset(&edit, 0, sizeof edit);
            edit.tx_idx = 0u;
            edit.transcript_strand = strand;
            edit.event.chrom_id = chrom;
            if (shape == 0u) {
                edit.ref = &allele;
                edit.ref_length = 1u;
                edit.event.start1 = 2500u;
                edit.event.end1 = 2500u;
                edit.event.feature_start1 = 2500u;
                edit.event.feature_end1 = 2500u;
                edit.event.ref_diff_length = 1u;
                edit.event.kind = (uint8_t)DUCKVEP_KIND_DEL;
                ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
                          duckvep_project_transcript_coordinate(
                              &tx, &exons, 0u, 2500u, &edit.first));
                edit.last = edit.first;
            } else {
                edit.alt = &allele;
                edit.alt_length = 1u;
                edit.event.start1 = 2500u;
                edit.event.end1 = 2500u;
                edit.event.insertion_boundary0 = 2500u;
                edit.event.interbase = 1u;
                edit.event.anchor_side =
                    (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
                edit.event.feature_start1 = 2501u;
                edit.event.feature_end1 = 2500u;
                edit.event.alt_diff_length = 1u;
                edit.event.kind = (uint8_t)DUCKVEP_KIND_INS;
                ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
                          duckvep_project_transcript_coordinate(
                              &tx, &exons, 0u, 2500u, &low));
                ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
                          duckvep_project_transcript_coordinate(
                              &tx, &exons, 0u, 2501u, &high));
                edit.first = strand > 0 ? low : high;
                edit.last = strand > 0 ? high : low;
            }
            ASSERT_EQ(DUCKVEP_HGVS_OK,
                      duckvep_hgvs_genomic_search_interval(
                          &edit.event, 5000u, &fetch_start, &fetch_end));
            ASSERT(fetch_end >= fetch_start);
            reference.start1 = fetch_start;
            reference.length =
                (size_t)((uint64_t)fetch_end - fetch_start + 1u);
            reference.chrom_id = chrom;
            reference.bases = reference_bytes;
            ASSERT(reference.length <= sizeof reference_bytes);
            memset(reference_bytes, 'A', reference.length);

            ASSERT_EQ(DUCKVEP_HGVS_OK,
                      duckvep_hgvs_dna_fact_build_genomic_shifted(
                          &tx, &exons, &reference, &edit, &fact));
            ASSERT_EQ(1000, fact.shift_offset);

            reference_bytes[strand > 0 ? reference.length - 1u : 0u] = 'C';
            ASSERT_EQ(DUCKVEP_HGVS_OK,
                      duckvep_hgvs_dna_fact_build_genomic_shifted(
                          &tx, &exons, &reference, &edit, &fact));
            ASSERT_EQ(999, fact.shift_offset);
        }
    }
    PASS();
}

TEST projection_known_forward_reverse_and_phase(void) {
    duckvep_coding_projection_t p;
    duckvep_event_t event;
    uint32_t cdna = 0u, genomic = 0u, exon_idx = 0u;
    uint32_t cds_start = 0u, cds_end = 0u;

    /* Forward: exons [100,109]=>cDNA 1..10, [200,209]=>11..20; CDS 103..205. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)1;
        s.excnt = 2u; s.cds_s = 103u; s.cds_e = 205u;
        s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u;  s.ce[0] = 10u; s.phase[0] = 0;
        s.es[1] = 200u; s.ee[1] = 209u; s.cs[1] = 11u; s.ce[1] = 20u; s.phase[1] = 0;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(1, duckvep_project_genomic_to_cdna(&s.tx, &s.ex, 0u, 100u, &cdna, &exon_idx));
        ASSERT_EQ(1u, cdna); ASSERT_EQ(0u, exon_idx);
        ASSERT_EQ(1, duckvep_project_cdna_to_genomic(&s.tx, &s.ex, 0u, 11u, &genomic, &exon_idx));
        ASSERT_EQ(200u, genomic); ASSERT_EQ(1u, exon_idx);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 105u, &p));
        ASSERT_EQ(6u, p.cdna_pos); ASSERT_EQ(3u, p.cds_pos); ASSERT_EQ(1u, p.protein_pos);
        ASSERT_EQ(1u, p.codon_start_cds); ASSERT_EQ(2u, (uint32_t)p.codon_offset);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 200u, &p));
        ASSERT_EQ(11u, p.cdna_pos); ASSERT_EQ(8u, p.cds_pos); ASSERT_EQ(3u, p.protein_pos);
        ASSERT_EQ(0, duckvep_project_genomic_to_cdna(&s.tx, &s.ex, 0u, 150u, &cdna, NULL));

        /* BaseTranscriptVariation projects the full feature endpoints and
         * permits an internal mapper Gap. The semantic edit projector remains
         * contiguous and sees only the changed base. */
        memset(&event, 0, sizeof event);
        event.start1 = 105u;
        event.end1 = 105u;
        event.ref_diff_length = 1u;
        event.feature_start1 = 105u;
        event.feature_end1 = 203u;
        ASSERT_EQ(1, duckvep_project_event_to_cds(
                         &s.tx, &s.ex, 0u, &event, &cds_start, &cds_end));
        ASSERT_EQ(3u, cds_start); ASSERT_EQ(3u, cds_end);
        ASSERT_EQ(1, duckvep_project_feature_to_cds(
                         &s.tx, &s.ex, 0u, &event, &cds_start, &cds_end));
        ASSERT_EQ(3u, cds_start); ASSERT_EQ(11u, cds_end);

        event.interbase = 1u;
        event.kind = (uint8_t)DUCKVEP_KIND_INS;
        event.feature_start1 = 106u;
        event.feature_end1 = 105u;
        event.insertion_boundary0 = 105u;
        event.start1 = 105u;
        event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
        ASSERT_EQ(1, duckvep_project_feature_to_cds(
                         &s.tx, &s.ex, 0u, &event, &cds_start, &cds_end));
        ASSERT_EQ(4u, cds_start); ASSERT_EQ(3u, cds_end);

        /* Insertion immediately before the first CDS base has reversed cDNA
         * 4,3 and does not satisfy VEP's overlap predicate. Moving it one base
         * right yields 5,4 and overlaps coding start 4. */
        event.feature_start1 = 103u;
        event.feature_end1 = 102u;
        event.insertion_boundary0 = 102u;
        event.start1 = 102u;
        ASSERT_EQ(1,
                  duckvep_project_feature_has_coding_precondition_unshifted(
                      &s.tx, &s.ex, 0u, &event));
        ASSERT_EQ(0,
                  duckvep_project_feature_overlaps_start_codon_unshifted(
                      &s.tx, &s.ex, 0u, &event));
        ASSERT_EQ(0,
                  duckvep_project_complete_feature_translation_bounds(
                      &s.tx, &s.ex, 0u, &event, 0u,
                      &cds_start, &cds_end));
        ASSERT_EQ(1,
                  duckvep_project_complete_feature_translation_bounds(
                      &s.tx, &s.ex, 0u, &event, 1u,
                      &cds_start, &cds_end));
        ASSERT_EQ(2u, cds_start); ASSERT_EQ(1u, cds_end);
        event.feature_start1 = 104u;
        event.feature_end1 = 103u;
        event.insertion_boundary0 = 103u;
        event.start1 = 103u;
        ASSERT_EQ(1,
                  duckvep_project_feature_overlaps_start_codon_unshifted(
                      &s.tx, &s.ex, 0u, &event));

        event.feature_start1 = 110u;
        event.feature_end1 = 109u;
        event.insertion_boundary0 = 109u;
        event.start1 = 109u;
        ASSERT_EQ(1,
                  duckvep_project_complete_feature_translation_bounds(
                      &s.tx, &s.ex, 0u, &event, 0u,
                      &cds_start, &cds_end));
        ASSERT_EQ(1, duckvep_project_feature_to_cds(
                         &s.tx, &s.ex, 0u, &event, &cds_start, &cds_end));
        ASSERT_EQ(8u, cds_start); ASSERT_EQ(7u, cds_end);

        event.interbase = 0u;
        event.feature_start1 = 110u;
        event.feature_end1 = 200u;
        ASSERT_EQ(0, duckvep_project_feature_to_cds(
                         &s.tx, &s.ex, 0u, &event, &cds_start, &cds_end));
    }

    /* Reverse: transcript order is high-to-low. CDS genomic 104..205 starts at 205. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)-1;
        s.excnt = 2u; s.cds_s = 104u; s.cds_e = 205u;
        s.es[0] = 200u; s.ee[0] = 209u; s.cs[0] = 1u;  s.ce[0] = 10u; s.phase[0] = 0;
        s.es[1] = 100u; s.ee[1] = 109u; s.cs[1] = 11u; s.ce[1] = 20u; s.phase[1] = 0;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(1, duckvep_project_genomic_to_cdna(&s.tx, &s.ex, 0u, 205u, &cdna, &exon_idx));
        ASSERT_EQ(5u, cdna); ASSERT_EQ(0u, exon_idx);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 203u, &p));
        ASSERT_EQ(7u, p.cdna_pos); ASSERT_EQ(3u, p.cds_pos); ASSERT_EQ(1u, p.protein_pos);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 104u, &p));
        ASSERT_EQ(16u, p.cdna_pos); ASSERT_EQ(12u, p.cds_pos); ASSERT_EQ(4u, p.protein_pos);
    }

    /* Positive Ensembl phase convention: phase 2 means first coding base is CDS position 3. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 110u; s.strand = (int8_t)1;
        s.excnt = 1u; s.cds_s = 100u; s.cds_e = 110u;
        s.es[0] = 100u; s.ee[0] = 110u; s.cs[0] = 1u; s.ce[0] = 11u; s.phase[0] = 2;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 100u, &p));
        ASSERT_EQ(1u, p.cdna_pos); ASSERT_EQ(3u, p.cds_pos); ASSERT_EQ(1u, p.protein_pos);
        ASSERT_EQ(1u, p.codon_start_cds); ASSERT_EQ(2u, (uint32_t)p.codon_offset);
        ASSERT_EQ(2u, (uint32_t)p.phase_offset);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 101u, &p));
        ASSERT_EQ(4u, p.cds_pos); ASSERT_EQ(2u, p.protein_pos); ASSERT_EQ(4u, p.codon_start_cds);
    }

    /* Non-zero exon_offset / tx_idx>0: returned exon_idx is absolute in exon model. */
    {
        static uint16_t chrom[2] = {0u, 0u};
        static uint32_t tstart[2] = {1u, 100u};
        static uint32_t tend[2] = {10u, 209u};
        static int8_t strand[2] = {(int8_t)1, (int8_t)1};
        static uint64_t flags[2] = {0u, 0u};
        static uint32_t exoff[2] = {0u, 1u};
        static uint16_t excnt[2] = {1u, 2u};
        static uint32_t cds_s[2] = {0u, 103u};
        static uint32_t cds_e[2] = {0u, 205u};
        static uint32_t es[3] = {1u, 100u, 200u};
        static uint32_t ee[3] = {10u, 109u, 209u};
        static uint32_t cs[3] = {1u, 1u, 11u};
        static uint32_t ce[3] = {10u, 10u, 20u};
        static int8_t phase[3] = {0, 0, 0};
        duckvep_transcript_model_t tx;
        duckvep_exon_model_t ex;
        memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
        tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
        tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
        tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
        ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = cs; ex.cdna_end1 = ce;
        ex.phase = phase; ex.end_phase = phase; ex.exon_count = 3u;
        ASSERT_EQ(1, duckvep_project_genomic_to_cdna(&tx, &ex, 1u, 200u, &cdna, &exon_idx));
        ASSERT_EQ(11u, cdna); ASSERT_EQ(2u, exon_idx);
        ASSERT_EQ(1, duckvep_project_coding_base(&tx, &ex, 1u, 200u, &p));
        ASSERT_EQ(8u, p.cds_pos); ASSERT_EQ(2u, p.exon_idx);
    }
    PASS();
}

/* VEP 116's independent-event path collapses an uploaded literal feature whose
 * two endpoints map to CDS across an intron into one outer-CDS replacement.
 * This compatibility edit is deliberately separate from the phased edit-set
 * projector, which must preserve exon/intron structure. */
TEST vep116_outer_cds_edit_spans_internal_intron_both_strands(void) {
    static const uint8_t forward_cds[13] = {
        'A','T','G','A','A','A','C','C','C','G','G','T','A'
    };
    static const uint8_t reverse_cds[12] = {
        'A','T','G','A','A','A','C','C','C','G','T','A'
    };
    static const uint8_t alternate[6] = {'G','G','G','A','A','A'};
    duckvep_haplotype_edit_t edit;
    duckvep_event_t event;

    {
        struct kprop_proj_scene s;
        duckvep_sequence_pool_t seq;
        uint64_t cds_offset = 0u;
        uint32_t cds_length = (uint32_t)sizeof forward_cds;

        memset(&s, 0, sizeof s);
        memset(&seq, 0, sizeof seq);
        memset(&event, 0, sizeof event);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u;
        s.strand = (int8_t)1; s.excnt = 2u;
        s.cds_s = 103u; s.cds_e = 205u;
        s.es[0] = 100u; s.ee[0] = 109u;
        s.cs[0] = 1u; s.ce[0] = 10u;
        s.es[1] = 200u; s.ee[1] = 209u;
        s.cs[1] = 11u; s.ce[1] = 20u;
        kprop_proj_scene_finish(&s);
        seq.cds_bytes = forward_cds;
        seq.cds_bytes_len = sizeof forward_cds;
        seq.cds_offset = &cds_offset;
        seq.cds_length = &cds_length;
        seq.transcript_count = 1u;

        event.feature_start1 = 105u;
        event.feature_end1 = 203u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_compat_vep116_outer_cds_edit_build(
                      &s.tx, &s.ex, &seq, 0u, (int8_t)1, &event,
                      alternate, (uint32_t)sizeof alternate, (int8_t)1,
                      &edit));
        ASSERT_EQ(3u, edit.cds_start);
        ASSERT_EQ(9u, edit.ref_len);
        ASSERT_EQ((uint32_t)sizeof alternate, edit.alt_len);
        ASSERT(edit.ref == forward_cds + 2u);
        ASSERT(edit.alt == alternate);
        ASSERT_EQ((int8_t)1, edit.variant_strand);

        event.feature_start1 = 105u;
        event.feature_end1 = 108u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_NON_CONTIGUOUS,
                  duckvep_compat_vep116_outer_cds_edit_build(
                      &s.tx, &s.ex, &seq, 0u, (int8_t)1, &event,
                      alternate, (uint32_t)sizeof alternate, (int8_t)1,
                      &edit));
        event.feature_start1 = 150u;
        event.feature_end1 = 203u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                  duckvep_compat_vep116_outer_cds_edit_build(
                      &s.tx, &s.ex, &seq, 0u, (int8_t)1, &event,
                      alternate, (uint32_t)sizeof alternate, (int8_t)1,
                      &edit));
    }

    {
        struct kprop_proj_scene s;
        duckvep_sequence_pool_t seq;
        uint64_t cds_offset = 0u;
        uint32_t cds_length = (uint32_t)sizeof reverse_cds;

        memset(&s, 0, sizeof s);
        memset(&seq, 0, sizeof seq);
        memset(&event, 0, sizeof event);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u;
        s.strand = (int8_t)-1; s.excnt = 2u;
        s.cds_s = 104u; s.cds_e = 205u;
        s.es[0] = 200u; s.ee[0] = 209u;
        s.cs[0] = 1u; s.ce[0] = 10u;
        s.es[1] = 100u; s.ee[1] = 109u;
        s.cs[1] = 11u; s.ce[1] = 20u;
        kprop_proj_scene_finish(&s);
        seq.cds_bytes = reverse_cds;
        seq.cds_bytes_len = sizeof reverse_cds;
        seq.cds_offset = &cds_offset;
        seq.cds_length = &cds_length;
        seq.transcript_count = 1u;

        event.feature_start1 = 106u;
        event.feature_end1 = 203u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_compat_vep116_outer_cds_edit_build(
                      &s.tx, &s.ex, &seq, 0u, (int8_t)-1, &event,
                      alternate, (uint32_t)sizeof alternate, (int8_t)-1,
                      &edit));
        ASSERT_EQ(3u, edit.cds_start);
        ASSERT_EQ(8u, edit.ref_len);
        ASSERT_EQ((uint32_t)sizeof alternate, edit.alt_len);
        ASSERT(edit.ref == reverse_cds + 2u);
        ASSERT(edit.alt == alternate);
        ASSERT_EQ((int8_t)-1, edit.variant_strand);
    }

    PASS();
}

TEST projection_start_codon_insert_drops_one_mapper_gap(void) {
    duckvep_event_t event;

    /* The start codon is split after cDNA base two. On either strand, an
     * insertion next to either exon edge maps to cDNA 3,2 after map_insert
     * discards the intronic Gap, and therefore overlaps cDNA 1..3. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u;
        s.strand = (int8_t)1; s.excnt = 2u;
        s.cds_s = 100u; s.cds_e = 209u;
        s.es[0] = 100u; s.ee[0] = 101u;
        s.cs[0] = 1u; s.ce[0] = 2u;
        s.es[1] = 200u; s.ee[1] = 209u;
        s.cs[1] = 3u; s.ce[1] = 12u;
        kprop_proj_scene_finish(&s);

        memset(&event, 0, sizeof event);
        event.kind = (uint8_t)DUCKVEP_KIND_INS;
        event.interbase = 1u;
        event.feature_start1 = 102u;
        event.feature_end1 = 101u;
        ASSERT_EQ(1,
                  duckvep_project_feature_overlaps_start_codon_unshifted(
                      &s.tx, &s.ex, 0u, &event));
        event.feature_start1 = 200u;
        event.feature_end1 = 199u;
        ASSERT_EQ(1,
                  duckvep_project_feature_overlaps_start_codon_unshifted(
                      &s.tx, &s.ex, 0u, &event));
    }

    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 201u;
        s.strand = (int8_t)-1; s.excnt = 2u;
        s.cds_s = 100u; s.cds_e = 201u;
        s.es[0] = 200u; s.ee[0] = 201u;
        s.cs[0] = 1u; s.ce[0] = 2u;
        s.es[1] = 100u; s.ee[1] = 109u;
        s.cs[1] = 3u; s.ce[1] = 12u;
        kprop_proj_scene_finish(&s);

        memset(&event, 0, sizeof event);
        event.kind = (uint8_t)DUCKVEP_KIND_INS;
        event.interbase = 1u;
        event.feature_start1 = 200u;
        event.feature_end1 = 199u;
        ASSERT_EQ(1,
                  duckvep_project_feature_overlaps_start_codon_unshifted(
                      &s.tx, &s.ex, 0u, &event));
        event.feature_start1 = 110u;
        event.feature_end1 = 109u;
        ASSERT_EQ(1,
                  duckvep_project_feature_overlaps_start_codon_unshifted(
                      &s.tx, &s.ex, 0u, &event));
    }
    PASS();
}

TEST projection_matches_bruteforce_for_any_small_transcript(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coordinate projection == brute-force transcript-order base walk";
    cfg.prop1 = prop_projection_matches_bruteforce;
    cfg.type_info[0] = &kprop_proj_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

TEST transcript_coordinate_matches_bruteforce_for_any_small_transcript(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "transcript coordinate == brute-force exon/intron walk";
    cfg.prop1 = prop_transcript_coordinate_matches_bruteforce;
    cfg.type_info[0] = &kprop_proj_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}
