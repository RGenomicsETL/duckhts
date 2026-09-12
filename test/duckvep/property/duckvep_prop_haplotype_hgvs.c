#include "duckvep_property.h"

static struct {
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
    uint32_t one_codon;
    uint32_t several_codons;
} g_haplotype_mnv_equivalence_cov;

static enum theft_trial_res prop_haplotype_snv_set_matches_equivalent_mnv(
    struct theft *t, void *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t island_edits[8];
    duckvep_haplotype_edit_t snv_edits[8];
    duckvep_haplotype_edit_t whole_edit;
    duckvep_edit_set_t island_set;
    duckvep_edit_set_t snv_set;
    duckvep_edit_set_t whole_set;
    duckvep_coding_context_t island_ctx;
    duckvep_coding_context_t snv_ctx;
    duckvep_coding_context_t whole_ctx;
    duckvep_sequence_delta_t island_delta;
    duckvep_sequence_delta_t snv_delta;
    duckvep_sequence_delta_t whole_delta;
    uint8_t island_alt_cds[80], snv_alt_cds[80], whole_alt_cds[80];
    uint8_t island_ref_pep[32], snv_ref_pep[32], whole_ref_pep[32];
    uint8_t island_alt_pep[32], snv_alt_pep[32], whole_alt_pep[32];
    size_t first = SIZE_MAX;
    size_t last = 0u;
    size_t snv_count = 0u;
    size_t i;
    duckvep_context_delta_status_t island_status;
    duckvep_context_delta_status_t snv_status;
    duckvep_context_delta_status_t whole_status;
    (void)t;

    if (duckvep_variant_cds_edit_set_build(
            &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand,
            island_edits, 8u, &island_set) != DUCKVEP_CDS_EDIT_OK ||
        island_set.count < 2u) {
        return THEFT_TRIAL_FAIL;
    }

    for (i = 0u; i < s->cds_lenv; i++) {
        if (s->cds[i] == s->expect_cds[i]) continue;
        if (first == SIZE_MAX) first = i;
        last = i;
    }
    if (first == SIZE_MAX || last - first + 1u > 8u) {
        return THEFT_TRIAL_FAIL;
    }

    for (i = last + 1u; i > first; i--) {
        size_t pos = i - 1u;
        duckvep_haplotype_edit_t *edit;
        if (s->cds[pos] == s->expect_cds[pos]) continue;
        if (snv_count >= sizeof snv_edits / sizeof snv_edits[0]) {
            return THEFT_TRIAL_FAIL;
        }
        edit = &snv_edits[snv_count++];
        edit->cds_start = (uint32_t)pos + 1u;
        edit->ref_len = 1u;
        edit->ref = s->cds + pos;
        edit->alt_len = 1u;
        edit->alt = s->expect_cds + pos;
        edit->variant_strand = s->strand;
    }
    if (snv_count < 2u) return THEFT_TRIAL_FAIL;

    whole_edit.cds_start = (uint32_t)first + 1u;
    whole_edit.ref_len = (uint32_t)(last - first + 1u);
    whole_edit.ref = s->cds + first;
    whole_edit.alt_len = whole_edit.ref_len;
    whole_edit.alt = s->expect_cds + first;
    whole_edit.variant_strand = s->strand;
    snv_set.edits = snv_edits;
    snv_set.count = snv_count;
    whole_set.edits = &whole_edit;
    whole_set.count = 1u;

#define BUILD_EQUIVALENT_CONTEXT(prefix, set)                                      \
    do {                                                                            \
        if (duckvep_coding_context_build(                                           \
                s->cds, s->cds_lenv, &(set), s->strand,                             \
                (duckvep_codon_table_t)s->ctab,                                     \
                prefix##_alt_cds, sizeof prefix##_alt_cds,                          \
                prefix##_ref_pep, sizeof prefix##_ref_pep,                          \
                prefix##_alt_pep, sizeof prefix##_alt_pep,                          \
                &prefix##_ctx) != DUCKVEP_CODING_CONTEXT_OK) {                      \
            return THEFT_TRIAL_FAIL;                                                \
        }                                                                           \
    } while (0)

    BUILD_EQUIVALENT_CONTEXT(island, island_set);
    BUILD_EQUIVALENT_CONTEXT(snv, snv_set);
    BUILD_EQUIVALENT_CONTEXT(whole, whole_set);
#undef BUILD_EQUIVALENT_CONTEXT

    if (island_ctx.alt_cds_len != snv_ctx.alt_cds_len ||
        island_ctx.alt_cds_len != whole_ctx.alt_cds_len ||
        memcmp(island_ctx.alt_cds, snv_ctx.alt_cds,
               island_ctx.alt_cds_len) != 0 ||
        memcmp(island_ctx.alt_cds, whole_ctx.alt_cds,
               island_ctx.alt_cds_len) != 0 ||
        island_ctx.ref_peptide_len != snv_ctx.ref_peptide_len ||
        island_ctx.ref_peptide_len != whole_ctx.ref_peptide_len ||
        island_ctx.alt_peptide_len != snv_ctx.alt_peptide_len ||
        island_ctx.alt_peptide_len != whole_ctx.alt_peptide_len ||
        memcmp(island_ctx.ref_peptide, snv_ctx.ref_peptide,
               island_ctx.ref_peptide_len) != 0 ||
        memcmp(island_ctx.ref_peptide, whole_ctx.ref_peptide,
               island_ctx.ref_peptide_len) != 0 ||
        memcmp(island_ctx.alt_peptide, snv_ctx.alt_peptide,
               island_ctx.alt_peptide_len) != 0 ||
        memcmp(island_ctx.alt_peptide, whole_ctx.alt_peptide,
               island_ctx.alt_peptide_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }

    island_status = duckvep_coding_context_delta_fill(
        &island_ctx, s->flags, &island_delta);
    snv_status = duckvep_coding_context_delta_fill(
        &snv_ctx, s->flags, &snv_delta);
    whole_status = duckvep_coding_context_delta_fill(
        &whole_ctx, s->flags, &whole_delta);
    if (island_status != DUCKVEP_CONTEXT_DELTA_OK ||
        snv_status != island_status || whole_status != island_status ||
        memcmp(&island_delta, &snv_delta, sizeof island_delta) != 0 ||
        memcmp(&island_delta, &whole_delta, sizeof island_delta) != 0) {
        return THEFT_TRIAL_FAIL;
    }

    if (s->strand > 0) g_haplotype_mnv_equivalence_cov.fwd++;
    else g_haplotype_mnv_equivalence_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START)
        g_haplotype_mnv_equivalence_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY)
        g_haplotype_mnv_equivalence_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP)
        g_haplotype_mnv_equivalence_cov.stop++;
    else return THEFT_TRIAL_FAIL;
    if (first / 3u == last / 3u)
        g_haplotype_mnv_equivalence_cov.one_codon++;
    else
        g_haplotype_mnv_equivalence_cov.several_codons++;
    return THEFT_TRIAL_PASS;
}

TEST haplotype_block_windows_keep_both_peptide_axes(void) {
    static const uint8_t reference[] = "ATGAAACCCGGGTTTTAAGCC";
    static const uint8_t inserted[] = "AGTC";
    size_t cases = 0u;
    for (size_t length = 18u; length <= 20u; length++) {
        for (int shift = -3; shift <= 3; shift += 3) {
            for (uint32_t position = 10u; position <= length; position++) {
                for (uint32_t alt_length = 0u; alt_length <= 4u; alt_length++) {
                    for (int strand = -1; strand <= 1; strand += 2) {
                        uint8_t ref[2][4], alt[2][4], cds[32], rp[16], ap[16];
                        duckvep_coding_context_t ctx;
                        duckvep_haplotype_edit_t edits[2] = {
                            {position, 1u, ref[0], alt_length, alt[0], 1},
                            {4u, shift > 0 ? 0u : shift < 0 ? 3u : 1u, ref[1],
                             shift < 0 ? 0u : shift > 0 ? 3u : 1u, alt[1], 1}
                        };
                        const uint8_t *tx_alt[2] = {inserted,
                            (const uint8_t *)(shift > 0 ? "CCC" : "T")};
                        for (size_t e = 0u; e < 2u; e++) {
                            for (uint32_t i = 0u; i < edits[e].ref_len; i++) {
                                size_t at = edits[e].cds_start - 1u +
                                    (strand > 0 ? i : edits[e].ref_len - 1u - i);
                                ref[e][i] = strand > 0 ? reference[at]
                                    : (uint8_t)kprop_complement_base((char)reference[at]);
                            }
                            for (uint32_t i = 0u; i < edits[e].alt_len; i++) {
                                uint8_t b = tx_alt[e][strand > 0 ? i : edits[e].alt_len - 1u - i];
                                alt[e][i] = strand > 0 ? b : (uint8_t)kprop_complement_base((char)b);
                            }
                        }
                        duckvep_edit_set_t set = {edits, 2u};
                        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                            reference, length, &set, (int8_t)strand, DUCKVEP_CODON_TABLE_STANDARD,
                            cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
                        duckvep_haplotype_edit_t ascending[2] = {edits[1], edits[0]};
                        duckvep_haplotype_block_t blocks[2];
                        size_t count;
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                            duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &count));
                        ASSERT_EQ(2u, count);
                        duckvep_coding_context_t saved = ctx;
                        duckvep_coding_peptide_window_t view;
                        ASSERT(duckvep_coding_context_block_window_open(&ctx, &blocks[1], &view));
                        ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
                        size_t ref_begin = ((size_t)position - 1u) / 3u * 3u;
                        size_t alt_begin = (size_t)((int64_t)ref_begin + shift);
                        size_t expected_ref_nt = length - ref_begin;
                        if (expected_ref_nt > 3u) expected_ref_nt = 3u;
                        size_t expected_alt_nt = ctx.alt_cds_len - alt_begin;
                        if (expected_alt_nt > alt_length + 2u) expected_alt_nt = alt_length + 2u;
                        ASSERT_EQ(ref_begin / 3u, view.ref_peptide_offset);
                        ASSERT_EQ(alt_begin / 3u, view.alt_peptide_offset);
                        ASSERT_EQ(expected_ref_nt, view.ref_nt_length);
                        ASSERT_EQ(expected_alt_nt, view.alt_nt_length);
                        for (int side = 0; side < 2; side++) {
                            size_t offset = side ? alt_begin / 3u : ref_begin / 3u;
                            size_t nts = side ? expected_alt_nt : expected_ref_nt;
                            const uint8_t *peptide = side ? ap : rp;
                            size_t whole = nts / 3u;
                            int partial = nts % 3u && !(whole == 1u && peptide[offset] == '*');
                            size_t residues = whole + (size_t)partial;
                            ASSERT_EQ(residues, side ? view.alt_length : view.ref_length);
                            for (size_t i = 0u; i < residues; i++) {
                                uint8_t expected = i == whole ? 'X' : peptide[offset + i];
                                ASSERT_EQ(expected, duckvep_coding_context_peptide_window_base(
                                    &ctx, &view, side, i));
                            }
                            ASSERT_EQ(0u, duckvep_coding_context_peptide_window_base(
                                &ctx, &view, side, residues));
                        }
                        if (view.ref_whole_length) {
                            uint32_t curated_position = (uint32_t)view.ref_peptide_offset + 1u;
                            const uint8_t curated = 'U';
                            uint8_t alternate = duckvep_coding_context_peptide_window_base(
                                &ctx, &view, 1, 0u);
                            ctx.ref_peptide_edit_position1 = &curated_position;
                            ctx.ref_peptide_edit_alt = &curated;
                            ctx.ref_peptide_edit_count = 1u;
                            ASSERT_EQ('U', duckvep_coding_context_peptide_window_base(&ctx, &view, 0, 0u));
                            ASSERT_EQ(alternate, duckvep_coding_context_peptide_window_base(&ctx, &view, 1, 0u));
                            ctx.ref_peptide_edit_count = 0u;
                            ctx.ref_peptide_edit_position1 = NULL;
                            ctx.ref_peptide_edit_alt = NULL;
                        }
                        for (int bad = 0; bad < 5; bad++) {
                            duckvep_haplotype_block_t invalid = blocks[1];
                            duckvep_coding_peptide_window_t cleared = {0};
                            if (bad == 0) invalid.length_diff++;
                            if (bad == 1) invalid.alt_start0++;
                            if (bad == 2) invalid.edit_count = ctx.applied_edits + 1u;
                            if (bad == 3) invalid.alt_start0 = SIZE_MAX;
                            if (bad == 4) invalid.cds_start = UINT32_MAX;
                            memset(&view, 0xff, sizeof view);
                            ASSERT(!duckvep_coding_context_block_window_open(&ctx, &invalid, &view));
                            ASSERT_MEM_EQ(&cleared, &view, sizeof view);
                        }
                        /* A real singleton must expose identical strings/coordinates
                         * through the single-event and interaction-block entry points. */
                        set.count = 1u;
                        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                            reference, length, &set, (int8_t)strand, DUCKVEP_CODON_TABLE_STANDARD,
                            cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                            duckvep_haplotype_partition(edits, 1u, blocks, 2u, &count));
                        duckvep_coding_peptide_window_t single;
                        ASSERT(duckvep_coding_context_peptide_window_open(&ctx, &single));
                        ASSERT(duckvep_coding_context_block_window_open(&ctx, blocks, &view));
                        ASSERT_MEM_EQ(&single, &view, sizeof view);
                        cases++;
                    }
                }
            }
        }
    }
    ASSERT_EQ(900u, cases);

    duckvep_coding_context_t ctx = {0};
    duckvep_haplotype_block_t block = {0};
    duckvep_coding_peptide_window_t view, empty = {0};
    memset(&view, 0xff, sizeof view);
    ASSERT(!duckvep_coding_context_block_window_open(&ctx, &block, &view));
    ASSERT_MEM_EQ(&empty, &view, sizeof view);
    ASSERT(!duckvep_coding_context_block_window_open(NULL, &block, &view));
    ASSERT(!duckvep_coding_context_block_window_open(&ctx, NULL, &view));
    ASSERT(!duckvep_coding_context_block_window_open(&ctx, &block, NULL));
    ctx.ref_peptide = (const uint8_t *)"MK"; ctx.ref_peptide_len = 2u;
    view.ref_peptide_offset = SIZE_MAX; view.ref_length = view.ref_whole_length = 2u;
    ASSERT_EQ(0u, duckvep_coding_context_peptide_window_base(&ctx, &view, 0, 1u));
    PASS();
}

TEST haplotype_restoring_edits_keep_changed_substitution_block(void) {
    static const uint8_t reference[] = "ATGGGTGGTGCTGATGATGCTGATGCTGATGGTTAA";
    for (int strand = -1; strand <= 1; strand += 2) {
        duckvep_haplotype_edit_t ascending[3] = {
            {4u, 3u, reference + 3u, 0u, NULL, (int8_t)strand},
            {11u, 1u, reference + 10u, 1u, (const uint8_t *)"G", (int8_t)strand},
            {13u, 0u, NULL, 3u, (const uint8_t *)"GCT", (int8_t)strand}
        };
        duckvep_haplotype_edit_t descending[3] = {ascending[2], ascending[1], ascending[0]};
        duckvep_edit_set_t set = {descending, 3u};
        uint8_t cds[48], rp[24], ap[24];
        duckvep_coding_context_t ctx;
        duckvep_haplotype_block_t blocks[3];
        duckvep_sequence_delta_t delta;
        size_t count;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
            reference, sizeof reference - 1u, &set, (int8_t)strand, DUCKVEP_CODON_TABLE_STANDARD,
            cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
        ASSERT_EQ(0u, ctx.cds_changed);
        ASSERT_EQ(sizeof reference - 1u, ctx.alt_cds_len);
        ASSERT_MEM_EQ(reference, cds, ctx.alt_cds_len);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
            duckvep_haplotype_partition(ascending, 3u, blocks, 3u, &count));
        ASSERT_EQ(3u, count);
        duckvep_coding_context_t saved = ctx;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_block_delta_fill(&ctx, ascending, 3u, blocks + 1u, 0u, &delta));
        ASSERT(delta.valid && delta.missense);
        ASSERT_EQ('A', delta.ref_aa);
        ASSERT_EQ('G', delta.alt_aa);
        ASSERT_EQ(4, delta.protein_pos);
        ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
            duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    }
    PASS();
}

TEST hgvs_haplotype_frameshift_uses_shifted_alternate_axis(void) {
    static const uint8_t reference[] = "ATGGGTCCTGCTGAACAATAA";
    duckvep_haplotype_edit_t ascending[2] = {
        {4u, 3u, reference + 3u, 0u, NULL, 1},
        {10u, 0u, NULL, 1u, (const uint8_t *)"T", 1}
    };
    duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
    duckvep_edit_set_t set = {descending, 2u};
    uint8_t cds[32], rp[16], ap[16];
    duckvep_coding_context_t ctx;
    duckvep_haplotype_block_t blocks[2];
    duckvep_hgvs_protein_operation_t facts[2];
    size_t blocks_used = 0u, facts_used = 0u, required = 0u;
    char rendered[64];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ctx.pre_cds_complete = ctx.post_cds_complete = 1u;
    ASSERT_EQ(19u, ctx.alt_cds_len);
    ASSERT_MEM_EQ("ATGCCTTGCTGAACAATAA", cds, ctx.alt_cds_len);
    ASSERT_MEM_EQ("MGPAEQ*", rp, ctx.ref_peptide_len);
    ASSERT_MEM_EQ("MPC*TI", ap, ctx.alt_peptide_len);
    ASSERT_EQ(4u, ctx.alt_first_stop_position1);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &blocks_used));
    ASSERT_EQ(2u, blocks_used);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, blocks, blocks_used, 0u, facts, 2u, &facts_used));
    ASSERT_EQ(2u, facts_used);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
        duckvep_hgvs_protein_render(&facts[0].fact, 0, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.Gly2del", rendered);
    ASSERT_EQ(3u, facts[1].fact.window.ref_peptide_offset);
    ASSERT_EQ(2u, facts[1].fact.window.alt_peptide_offset);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
        duckvep_hgvs_protein_render(&facts[1].fact, 0, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.Ala4CysfsTer2", rendered);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        facts, facts_used, 1, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.[(Gly2del;Ala4CysfsTer2)]", rendered);
    /* Full source-record replacements retain VCF anchor bases. The anchor
     * occupies CDS position 3 but does not change the start codon. */
    ascending[0] = (duckvep_haplotype_edit_t){3u, 4u, reference + 2u,
        1u, (const uint8_t *)"G", 1};
    ascending[1] = (duckvep_haplotype_edit_t){9u, 1u, reference + 8u,
        2u, (const uint8_t *)"TT", 1};
    descending[0] = ascending[1]; descending[1] = ascending[0];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ctx.pre_cds_complete = ctx.post_cds_complete = 1u;
    ASSERT_MEM_EQ("ATGCCTTGCTGAACAATAA", cds, ctx.alt_cds_len);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &blocks_used));
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, blocks, blocks_used, 0u, facts, 2u, &facts_used));
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        facts, facts_used, 1, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.[(Gly2del;Ala4CysfsTer2)]", rendered);
    size_t exact = required;
    for (size_t capacity = 0u; capacity <= exact + 1u; capacity++) {
        char guarded[80];
        memset(guarded, '!', sizeof guarded);
        ASSERT_EQ(capacity <= exact ? DUCKVEP_HGVS_BUFFER_TOO_SMALL : DUCKVEP_HGVS_OK,
            duckvep_hgvs_protein_haplotype_render(facts, facts_used, 1,
                guarded + 1u, capacity, &required));
        ASSERT_EQ(exact, required);
        ASSERT_EQ('!', guarded[0]);
        ASSERT_EQ('!', guarded[capacity + 1u]);
        if (capacity) ASSERT_EQ('\0', guarded[capacity]);
    }
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL,
        duckvep_hgvs_protein_haplotype_render(facts, facts_used, 1, NULL, 0u, &required));
    ASSERT_EQ(exact, required);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
        duckvep_hgvs_protein_haplotype_render(facts, facts_used, 0, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.[Gly2del;Ala4CysfsTer2]", rendered);
    duckvep_hgvs_protein_operation_t wrong[2] = {facts[1], facts[0]};
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
        duckvep_hgvs_protein_haplotype_render(wrong, 2u, 1, rendered, sizeof rendered, &required));
    ASSERT_EQ(0u, required);
    ASSERT_STR_EQ("", rendered);
    /* Inserting before either T at the CCT/T junction produces the same CDS.
     * The first position shares a coding block with the earlier deletion;
     * this must not absorb that deletion into the protein frameshift. */
    ascending[0] = (duckvep_haplotype_edit_t){4u, 3u, reference + 3u, 0u, NULL, 1};
    ascending[1] = (duckvep_haplotype_edit_t){9u, 0u, NULL, 1u, (const uint8_t *)"T", 1};
    descending[0] = ascending[1]; descending[1] = ascending[0];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ctx.pre_cds_complete = ctx.post_cds_complete = 1u;
    ASSERT_MEM_EQ("ATGCCTTGCTGAACAATAA", cds, ctx.alt_cds_len);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &blocks_used));
    ASSERT_EQ(1u, blocks_used);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, blocks, blocks_used, 0u, facts, 2u, &facts_used));
    ASSERT_EQ(2u, facts_used);
    ASSERT_EQ(1u, facts[0].span.edit_count);
    ASSERT_EQ(1u, facts[1].span.edit_count);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        facts, facts_used, 1, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.[(Gly2del;Ala4CysfsTer2)]", rendered);
    PASS();
}

TEST hgvs_haplotype_start_contrast_uses_unpadded_cds(void) {
    static const uint8_t dna[] = "ATGGGTCCTGCTGAACAATAA";
    for (unsigned phase = 0u; phase < 3u; phase++) {
        for (unsigned base = 0u; base < 4u; base++) for (int strand = -1; strand <= 1; strand += 2) {
            uint8_t reference[24], cds[40], rp[16], ap[16];
            memset(reference, 'N', phase);
            memcpy(reference + phase, dna, sizeof dna - 1u);
            size_t length = phase + sizeof dna - 1u;
            uint8_t alternate[4] = {"ACGT"[base], 'A', 'A', 'A'};
            uint8_t first_ref = 'G', second_ref = 'G', second_alt = 'T';
            if (strand < 0) {
                uint8_t saved[4]; memcpy(saved, alternate, sizeof saved);
                for (size_t i = 0u; i < 4u; i++)
                    alternate[i] = (uint8_t)kprop_complement_base((char)saved[3u - i]);
                first_ref = second_ref = 'C'; second_alt = 'A';
            }
            duckvep_haplotype_edit_t edits[2] = {
                {phase + 3u, 1u, &first_ref, 4u, alternate, 1},
                {phase + 10u, 1u, &second_ref, 1u, &second_alt, 1}
            };
            duckvep_haplotype_edit_t descending[2] = {edits[1], edits[0]};
            duckvep_edit_set_t set = {descending, 2u};
            duckvep_coding_context_t ctx;
            duckvep_haplotype_block_t blocks[2];
            duckvep_hgvs_protein_operation_t operations[4];
            size_t block_count, operation_count;
            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                reference, length, &set, (int8_t)strand, DUCKVEP_CODON_TABLE_STANDARD,
                cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
            ctx.cds_phase_padding = (uint8_t)phase;
            ctx.pre_cds_complete = ctx.post_cds_complete = 1u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                duckvep_haplotype_partition(edits, 2u, blocks, 2u, &block_count));
            duckvep_sequence_delta_t delta;
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                duckvep_coding_context_block_delta_fill(&ctx, edits, 2u, blocks, 0u, &delta));
            ASSERT(delta.start_lost);
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
                &ctx, NULL, edits, 2u, blocks, block_count, 0u,
                operations, 4u, &operation_count));
            ASSERT(operation_count > 0u);
            ASSERT_EQ("ACGT"[base] != 'G',
                operations[0].fact.shape == DUCKVEP_HGVS_PROTEIN_START_LOST);
        }
    }
    PASS();
}

TEST hgvs_haplotype_renderer_validates_complete_allele(void) {
    duckvep_hgvs_protein_operation_t operations[2] = {0};
    uint32_t valid = DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_fact_build_single_residue(
        2u, 'R', 'G', valid, DUCKVEP_COMPAT_VEP_116,
        &operations[0].fact));
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_fact_build_single_residue(
        4u, 'Q', '*', valid, DUCKVEP_COMPAT_VEP_116,
        &operations[1].fact));
    char output[80]; size_t required;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        operations, 2u, 1, output, sizeof output, &required));
    /* HGVS 21.1.4 protein/alleles: prediction parentheses enclose all cis edits. */
    ASSERT_STR_EQ("p.[(Arg2Gly;Gln4Ter)]", output);
    for (int predicted = 0; predicted < 2; predicted++) {
        const char *expected = predicted ? "p.(=)" : "p.=";
        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
            NULL, 0u, predicted, output, sizeof output, &required));
        ASSERT_STR_EQ(expected, output);
        ASSERT_EQ(strlen(expected), required);
        ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_haplotype_render(
            NULL, 0u, predicted, NULL, 0u, &required));
        ASSERT_EQ(strlen(expected), required);
    }
    for (unsigned invalid = 0u; invalid < 10u; invalid++) {
        duckvep_hgvs_protein_operation_t changed[2] = {operations[0], operations[1]};
        const duckvep_hgvs_protein_operation_t *input = changed;
        int predicted = 1;
        switch (invalid) {
            case 0: input = NULL; break;
            case 1: predicted = -1; break;
            case 2: predicted = 2; break;
            case 3: changed[1].fact.shape = DUCKVEP_HGVS_PROTEIN_EQUAL; break;
            case 4: changed[1].fact.compatibility_profile = UINT8_MAX; break;
            case 5: changed[1].fact.first_position1 = 3u; break;
            case 6: changed[1].fact.first_position1 = 0u; break;
            case 7: changed[1] = changed[0]; break;
            case 8: changed[1].fact.reference.bases = (const uint8_t *)"RQ"; break;
            case 9: changed[1].fact.reference.length = 2u; break;
        }
        required = SIZE_MAX; memset(output, '!', sizeof output);
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_protein_haplotype_render(
            input, 2u, predicted, output, sizeof output, &required));
        ASSERT_EQ(0u, required);
        ASSERT_STR_EQ("", output);
    }
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_protein_haplotype_render(
        operations, 2u, 1, NULL, 1u, &required));
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_protein_haplotype_render(
        operations, 2u, 1, output, sizeof output, NULL));
    PASS();
}

/* Independently constructed blocks have no contract for struct padding. */
static int kprop_hgvs_block_equal(const duckvep_haplotype_block_t *a,
    const duckvep_haplotype_block_t *b) {
    return a->edit_begin == b->edit_begin && a->edit_count == b->edit_count &&
        a->cds_start == b->cds_start && a->ref_len == b->ref_len &&
        a->alt_start0 == b->alt_start0 && a->alt_len == b->alt_len &&
        a->length_diff == b->length_diff && a->flags == b->flags;
}

TEST hgvs_haplotype_prepared_reference_is_borrowed(void) {
    static const uint8_t reference[] = "CTGGCCTAA";
    uint8_t prepared[] = "MA*", cds[16], rp[8], ap[8];
    duckvep_edit_set_t empty = {NULL, 0u};
    duckvep_coding_context_t ctx;
    duckvep_hgvs_protein_operation_t operations[2];
    size_t count, required;
    char text[64];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(reference,
        sizeof reference - 1u, &empty, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_MEM_EQ("LA*", rp, ctx.ref_peptide_len);
    unsigned char saved_context[sizeof ctx];
    memcpy(saved_context, &ctx, sizeof ctx);
    duckvep_hgvs_protein_reference_t view = {prepared, sizeof prepared - 1u};
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, operations, 1u, &count));
    ASSERT_EQ(1u, count);
    ASSERT_EQ(0u, operations[0].span.edit_count);
    ASSERT(operations[0].fact.reference.bases == prepared);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        operations, count, 1, text, sizeof text, &required));
    ASSERT_STR_EQ("p.(Met1Leu)", text);
    ASSERT_MEM_EQ(saved_context, &ctx, sizeof ctx);
    ASSERT_MEM_EQ("MA*", prepared, view.length);
    ASSERT_MEM_EQ("LA*", rp, ctx.ref_peptide_len);
    ASSERT_MEM_EQ("LA*", ap, ctx.alt_peptide_len);
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, NULL, 0u, &count));
    ASSERT_EQ(0u, count);
    for (unsigned invalid = 0u; invalid < 5u; invalid++) {
        duckvep_hgvs_protein_reference_t bad = view;
        duckvep_hgvs_status_t expected = DUCKVEP_HGVS_MISSING_PEPTIDE;
        uint8_t incomplete[] = {'M', 0u, '*'};
        if (invalid == 0u) { bad.bases = NULL; expected = DUCKVEP_HGVS_INVALID_ARG; }
        if (invalid == 1u) bad.length = 0u;
        if (invalid == 2u) bad.length++;
        if (invalid == 3u) { bad.length = SIZE_MAX; expected = DUCKVEP_HGVS_OUT_OF_RANGE; }
        if (invalid == 4u) { bad.bases = incomplete; expected = DUCKVEP_HGVS_MISSING_PEPTIDE; }
        count = SIZE_MAX;
        ASSERT_EQ(expected, duckvep_hgvs_protein_haplotype_build(&ctx, &bad,
            NULL, 0u, NULL, 0u, 0u, operations, 2u, &count));
        ASSERT_EQ(0u, count);
    }
    view.bases = NULL; view.length = 0u;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, operations, 2u, &count));
    ASSERT_EQ(0u, count);
    ASSERT_MEM_EQ(saved_context, &ctx, sizeof ctx);
    rp[1] = ap[1] = 0u;
    ASSERT_EQ(DUCKVEP_HGVS_MISSING_PEPTIDE, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, operations, 2u, &count));
    ASSERT_EQ(0u, count);
    rp[1] = ap[1] = 'A';
    view.bases = (const uint8_t *)"MAW*"; view.length = 4u;
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, operations, 1u, &count));
    ASSERT_EQ(0u, count);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, operations, 2u, &count));
    ASSERT_EQ(2u, count);
    ASSERT_EQ(0u, operations[0].span.edit_count);
    ASSERT_EQ(0u, operations[1].span.edit_count);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        operations, count, 1, text, sizeof text, &required));
    ASSERT_STR_EQ("p.[(Met1Leu;Trp3Ter)]", text);
    view.bases = (const uint8_t *)"MA"; view.length = 2u;
    ASSERT_EQ(DUCKVEP_HGVS_NOT_APPLICABLE, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        NULL, 0u, NULL, 0u, 0u, operations, 2u, &count));
    ASSERT_EQ(0u, count);
    ASSERT_MEM_EQ(saved_context, &ctx, sizeof ctx);
    PASS();
}

TEST hgvs_haplotype_restored_frame_splits_unchanged_residue(void) {
    static const uint8_t reference[] = "ATGCGGCATTTCTATGAATAA";
    duckvep_haplotype_edit_t ascending[2] = {
        {4u, 0u, NULL, 1u, (const uint8_t *)"C", 1},
        {13u, 1u, reference + 12u, 0u, NULL, 1}
    };
    duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
    duckvep_edit_set_t edits = {descending, 2u};
    uint8_t cds[32], rp[16], ap[16], replayed[2][16];
    duckvep_coding_context_t ctx;
    duckvep_haplotype_block_t block;
    duckvep_hgvs_protein_operation_t operations[2];
    size_t blocks_used, count, required, replayed_length;
    char text[96];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(reference,
        sizeof reference - 1u, &edits, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_MEM_EQ("MRHFYE*", rp, ctx.ref_peptide_len);
    ASSERT_MEM_EQ("MPAFHE*", ap, ctx.alt_peptide_len);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, &block, 1u, &blocks_used));
    ASSERT_EQ(1u, blocks_used);
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, &block, 1u, 0u, operations, 1u, &count));
    ASSERT_EQ(0u, count);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, &block, 1u, 0u, operations, 2u, &count));
    ASSERT_EQ(2u, count);
    ASSERT(kprop_hgvs_block_equal(&block, &operations[0].span));
    ASSERT(kprop_hgvs_block_equal(&block, &operations[1].span));
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        operations, count, 1, text, sizeof text, &required));
    ASSERT_STR_EQ("p.[(Arg2_His3delinsProAla;Tyr5His)]", text);
    ASSERT(kprop_hgvs_protein_fact_replay(&operations[1].fact, rp, ctx.ref_peptide_len,
        replayed[0], sizeof replayed[0], &replayed_length));
    ASSERT(kprop_hgvs_protein_fact_replay(&operations[0].fact, replayed[0], replayed_length,
        replayed[1], sizeof replayed[1], &replayed_length));
    ASSERT_EQ(ctx.alt_peptide_len, replayed_length);
    ASSERT_MEM_EQ(ap, replayed[1], replayed_length);
    PASS();
}

TEST hgvs_haplotype_equal_suffix_requires_aligned_protein_axes(void) {
    static const uint8_t reference[] = "ATGGATGATTAA";
    for (int strand = -1; strand <= 1; strand += 2) {
        for (unsigned deletion = 0u; deletion < 2u; deletion++) {
            duckvep_haplotype_edit_t edit = {4u, deletion ? 3u : 0u,
                deletion ? reference + 3u : NULL, deletion ? 0u : 3u,
                deletion ? NULL : (const uint8_t *)"GAT", (int8_t)strand};
            duckvep_edit_set_t set = {&edit, 1u};
            uint8_t cds[32], rp[16], ap[16], replayed[16];
            duckvep_coding_context_t ctx;
            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(reference,
                sizeof reference - 1u, &set, (int8_t)strand, DUCKVEP_CODON_TABLE_STANDARD,
                cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
            ctx.post_cds_complete = 1u;
            duckvep_haplotype_block_t block;
            size_t block_count, count, required, replayed_length;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                duckvep_haplotype_partition(&edit, 1u, &block, 1u, &block_count));
            ASSERT_EQ(1u, block_count);
            duckvep_haplotype_block_t saved_block = block;
            unsigned char saved_context[sizeof ctx];
            memcpy(saved_context, &ctx, sizeof ctx);
            duckvep_hgvs_protein_reference_t view = {
                (const uint8_t *)(deletion ? "MDD*" : "MDDD*"), deletion ? 4u : 5u};
            duckvep_hgvs_protein_operation_t operation;
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
                &edit, 1u, &block, 1u, 0u, &operation, 1u, &count));
            ASSERT_EQ(deletion, count);
            ASSERT_MEM_EQ(saved_context, &ctx, sizeof ctx);
            ASSERT(kprop_hgvs_block_equal(&saved_block, &block));
            ASSERT_EQ(1u, ctx.applied_edits);
            ASSERT(ctx.cds_changed);
            if (deletion) {
                ASSERT_EQ(1u, operation.span.edit_count);
                ASSERT(kprop_hgvs_protein_fact_replay(&operation.fact, view.bases, view.length,
                    replayed, sizeof replayed, &replayed_length));
                ASSERT_EQ(ctx.alt_peptide_len, replayed_length);
                ASSERT_MEM_EQ(ap, replayed, replayed_length);
            } else {
                ASSERT_EQ(view.length, ctx.alt_peptide_len);
                ASSERT_MEM_EQ(view.bases, ap, view.length);
            }
            char text[64];
            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
                &operation, count, 1, text, sizeof text, &required));
            ASSERT_STR_EQ(deletion ? "p.(Asp2del)" : "p.(=)", text);
            /* Equal displayed protein does not bypass physical edit validation. */
            block.edit_count = 0u;
            count = SIZE_MAX;
            ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
                &edit, 1u, &block, 1u, 0u, &operation, 1u, &count));
            ASSERT_EQ(0u, count);
        }
    }
    PASS();
}

TEST hgvs_haplotype_duplication_uses_supplied_reference(void) {
    /* The compound API replays its declared or explicitly supplied reference.
     * The independent TVA duplication query's default-table CDS translation
     * must not replace that authority, including sparse Translation SeqEdits. */
    enum { SPARSE_SEQEDIT = 1u, SUPPLIED_REFERENCE = 2u };
    static const struct {
        const char *cds;
        uint8_t table;
        uint32_t insertion1;
        const char *reference;
        const char *curated;
        const char *alternate;
        duckvep_hgvs_protein_shape_t shape;
    } cases[] = {
        {"ATGAGAGCCTAA", 5u, 5u, "MSA*", "MRA*", "MSRA*", DUCKVEP_HGVS_PROTEIN_INSERTION},
        {"ATGCTGGCCTAA", 26u, 7u, "MAA*", "MLA*", "MAAA*", DUCKVEP_HGVS_PROTEIN_DUPLICATION}
    };
    static const duckvep_compat_profile_t profiles[] = {DUCKVEP_COMPAT_VEP_116, DUCKVEP_COMPAT_STRICT};
    size_t checked = 0u;
    for (size_t i = 0u; i < sizeof cases / sizeof cases[0]; i++) {
        for (unsigned mode = 0u; mode < 4u; mode++) {
            for (size_t p = 0u; p < sizeof profiles / sizeof profiles[0]; p++) {
                duckvep_haplotype_edit_t edit = {cases[i].insertion1, 0u, NULL,
                    3u, (const uint8_t *)"GCC", 1};
                duckvep_edit_set_t set = {&edit, 1u};
                uint8_t cds[32], rp[16], ap[16], replayed[2][16];
                duckvep_coding_context_t context;
                ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                    (const uint8_t *)cases[i].cds, 12u, &set, 1,
                    (duckvep_codon_table_t)cases[i].table,
                    cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &context));
                context.post_cds_complete = 1u;
                context.compatibility_profile = (uint8_t)profiles[p];
                ASSERT_EQ(4u, context.ref_peptide_len);
                ASSERT_MEM_EQ(cases[i].reference, rp, 4u);
                ASSERT_EQ(5u, context.alt_peptide_len);
                ASSERT_MEM_EQ(cases[i].alternate, ap, 5u);
                uint32_t curated_position = 2u;
                uint8_t curated_residue = (uint8_t)cases[i].curated[1];
                if (mode & SPARSE_SEQEDIT) {
                    context.ref_peptide_edit_position1 = &curated_position;
                    context.ref_peptide_edit_alt = &curated_residue;
                    context.ref_peptide_edit_count = 1u;
                }
                const uint8_t *reference = (const uint8_t *)(mode ? cases[i].curated : cases[i].reference);
                duckvep_hgvs_protein_reference_t view = {reference, 4u};
                duckvep_haplotype_block_t block;
                size_t blocks, count;
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                    duckvep_haplotype_partition(&edit, 1u, &block, 1u, &blocks));
                ASSERT_EQ(1u, blocks);
                duckvep_coding_context_t saved_context = context;
                duckvep_haplotype_edit_t saved_edit = edit;
                duckvep_haplotype_block_t saved_block = block;
                duckvep_hgvs_protein_operation_t operations[3];
                ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
                    &context, (mode & SUPPLIED_REFERENCE) ? &view : NULL, &edit, 1u, &block, 1u,
                    0u, operations, 3u, &count));
                ASSERT(count > 0u && count <= 3u);
                if (!mode) {
                    ASSERT_EQ(1u, count);
                    ASSERT_EQ(cases[i].shape, operations[0].fact.shape);
                }
                memcpy(replayed[0], reference, 4u);
                size_t length = 4u;
                unsigned side = 0u;
                for (size_t j = count; j > 0u; j--) {
                    ASSERT(kprop_hgvs_protein_fact_replay(&operations[j - 1u].fact,
                        replayed[side], length, replayed[side ^ 1u], sizeof replayed[0], &length));
                    side ^= 1u;
                }
                ASSERT_EQ(5u, length);
                ASSERT_MEM_EQ(cases[i].alternate, replayed[side], length);
                ASSERT_MEM_EQ(&saved_context, &context, sizeof context);
                ASSERT_MEM_EQ(&saved_edit, &edit, sizeof edit);
                ASSERT_MEM_EQ(&saved_block, &block, sizeof block);
                ASSERT_MEM_EQ(cases[i].reference, rp, 4u);
                ASSERT_MEM_EQ(cases[i].alternate, ap, 5u);
                checked++;
            }
        }
    }
    ASSERT_EQ(16u, checked);
    PASS();
}

TEST hgvs_haplotype_terminal_insertion_replays_complete_sequence(void) {
    static const uint8_t reference[] = "ATGGATGCTTAA", prepared[] = "MDAD*";
    duckvep_haplotype_edit_t edit = {4u, 0u, NULL, 6u, (const uint8_t *)"GATGCT", 1};
    duckvep_edit_set_t set = {&edit, 1u};
    uint8_t cds[32], rp[16], ap[16], replayed[16];
    duckvep_coding_context_t ctx;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(reference,
        sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ctx.post_cds_complete = 1u;
    duckvep_haplotype_block_t block;
    size_t blocks, count, required, replayed_length;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(&edit, 1u, &block, 1u, &blocks));
    duckvep_hgvs_protein_reference_t view = {prepared, sizeof prepared - 1u};
    duckvep_hgvs_protein_operation_t operation;
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        &edit, 1u, &block, blocks, 0u, &operation, 1u, &count));
    ASSERT_EQ(1u, count);
    ASSERT_EQ(1u, operation.span.edit_count);
    ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_INSERTION, operation.fact.shape);
    ASSERT_EQ(1u, operation.fact.alt_length);
    ASSERT(kprop_hgvs_protein_fact_replay(&operation.fact, prepared, sizeof prepared - 1u,
        replayed, sizeof replayed, &replayed_length));
    ASSERT_EQ(6u, replayed_length);
    ASSERT_MEM_EQ("MDADA*", replayed, replayed_length);
    ASSERT_EQ(ctx.alt_peptide_len, replayed_length);
    ASSERT_MEM_EQ(ap, replayed, replayed_length);
    char text[64];
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        &operation, count, 1, text, sizeof text, &required));
    ASSERT_STR_EQ("p.(Asp4_Ter5insAla)", text);
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        &edit, 1u, &block, blocks, 0u, NULL, 0u, &count));
    ASSERT_EQ(0u, count);
    /* A protein-neutral physical edit does not turn an isolated terminal
     * substitution into a deletion of the curated reference residue. */
    edit = (duckvep_haplotype_edit_t){6u, 1u, (const uint8_t *)"C",
        1u, (const uint8_t *)"T", 1};
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        (const uint8_t *)"ATGGCCTAA", 9u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ctx.post_cds_complete = 1u;
    view.bases = (const uint8_t *)"MAW*"; view.length = 4u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(&edit, 1u, &block, 1u, &blocks));
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(&ctx, &view,
        &edit, 1u, &block, blocks, 0u, &operation, 1u, &count));
    ASSERT_EQ(1u, count);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        &operation, count, 1, text, sizeof text, &required));
    ASSERT_STR_EQ("p.(Trp3Ter)", text);
    PASS();
}

TEST hgvs_haplotype_normalization_combines_interacting_spans(void) {
    static const uint8_t reference[] = "ATGGCTGCTGCTGCTGCTGAATAA";
    duckvep_haplotype_edit_t ascending[2] = {
        {4u, 0u, NULL, 3u, (const uint8_t *)"GCT", 1},
        {10u, 3u, reference + 9u, 3u, (const uint8_t *)"GAT", 1}
    };
    duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
    duckvep_edit_set_t set = {descending, 2u};
    uint8_t cds[40], rp[16], ap[16], replayed[16];
    duckvep_coding_context_t ctx;
    duckvep_haplotype_block_t blocks[2];
    duckvep_hgvs_protein_operation_t operations[2];
    size_t blocks_used, count, required, replayed_length;
    char rendered[64];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_EQ(8u, ctx.ref_peptide_len);
    ASSERT_EQ(9u, ctx.alt_peptide_len);
    ASSERT_MEM_EQ("MAAAAAE*", rp, 8u);
    ASSERT_MEM_EQ("MAAADAAE*", ap, 9u);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &blocks_used));
    ASSERT_EQ(2u, blocks_used);
    duckvep_haplotype_block_t saved[2] = {blocks[0], blocks[1]};
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, blocks, 2u, 0u, operations, 2u, &count));
    ASSERT_EQ(1u, count);
    ASSERT_EQ(0u, operations[0].span.edit_begin);
    ASSERT_EQ(2u, operations[0].span.edit_count);
    ASSERT_MEM_EQ(saved, blocks, sizeof blocks);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
        &operations[0].fact, 0, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.Ala4_Ala5insAsp", rendered);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        operations, count, 1, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.(Ala4_Ala5insAsp)", rendered);
    ASSERT(kprop_hgvs_protein_fact_replay(&operations[0].fact, rp, ctx.ref_peptide_len,
        replayed, sizeof replayed, &replayed_length));
    ASSERT_EQ(ctx.alt_peptide_len, replayed_length);
    ASSERT_MEM_EQ(ap, replayed, replayed_length);
    /* Normalizing the two source edits independently moves the insertion
     * across the substitution. Those individually valid descriptions no
     * longer reproduce the complete carried protein. */
    duckvep_coding_context_t independent[2];
    duckvep_hgvs_protein_fact_t separate[2];
    uint8_t separate_cds[2][40], separate_ref[2][16], separate_alt[2][16];
    for (size_t i = 0u; i < 2u; i++) {
        duckvep_edit_set_t one = {ascending + i, 1u};
        duckvep_sequence_delta_t delta;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
            reference, sizeof reference - 1u, &one, 1, DUCKVEP_CODON_TABLE_STANDARD,
            separate_cds[i], sizeof separate_cds[i], separate_ref[i], sizeof separate_ref[i],
            separate_alt[i], sizeof separate_alt[i], independent + i));
        independent[i].pre_cds_complete = independent[i].post_cds_complete = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_delta_fill(independent + i, 0u, &delta));
        ASSERT_EQ(DUCKVEP_HGVS_OK,
            duckvep_hgvs_protein_fact_build(independent + i, &delta, separate + i));
        ASSERT_EQ(DUCKVEP_HGVS_OK,
            duckvep_hgvs_protein_render(separate + i, 0, rendered, sizeof rendered, &required));
        ASSERT_STR_EQ(i ? "p.Ala4Asp" : "p.Ala6dup", rendered);
    }
    uint8_t wrong[2][16]; size_t wrong_length;
    ASSERT(kprop_hgvs_protein_fact_replay(separate, rp, ctx.ref_peptide_len,
        wrong[0], sizeof wrong[0], &wrong_length));
    ASSERT(kprop_hgvs_protein_fact_replay(separate + 1u, wrong[0], wrong_length,
        wrong[1], sizeof wrong[1], &wrong_length));
    ASSERT_EQ(ctx.alt_peptide_len, wrong_length);
    ASSERT_MEM_EQ("MAADAAAE*", wrong[1], wrong_length);
    ASSERT(memcmp(ap, wrong[1], wrong_length) != 0);
    PASS();
}

static enum theft_alloc_res kprop_hgvs_inframe_haplotype_alloc_mode(
    struct theft *t, void **instance, size_t restoration_cell) {
    static const uint8_t codons[4][3] = {{'G','C','T'}, {'G','G','T'}, {'G','A','T'}, {'C','C','T'}};
    int restoring = restoration_cell != SIZE_MAX;
    struct kprop_haplo_case *c = calloc(1u, sizeof *c);
    if (!c) return THEFT_ALLOC_ERROR;
    memcpy(c->ref, "ATG", 3u);
    for (size_t at = 3u; at < KPROP_HAPLO_CDS_LEN - 3u; at += 3u)
        memcpy(c->ref + at, codons[kprop_bounded(t, 4u)], 3u);
    memcpy(c->ref + KPROP_HAPLO_CDS_LEN - 3u, "TAA", 3u);
    c->transcript_strand = kprop_bounded(t, 2u) ? 1 : -1;
    size_t target = 2u + (size_t)kprop_bounded(t, KPROP_HAPLO_MAX_EDITS - 1u);
    size_t codon = 1u;
    size_t repeated = 0u, replaced = 0u;
    size_t orientations = 0u;
    if (restoring) {
        /* REF A-A-B becomes A-A-B through deletion of the first A,
         * B-to-A replacement and insertion of B. The middle edit compares
         * different codons despite identical complete CDS sequences. */
        target = 3u;
        codon += restoration_cell % 7u; restoration_cell /= 7u;
        repeated = restoration_cell % 4u; restoration_cell /= 4u;
        replaced = (repeated + 1u + restoration_cell % 3u) % 4u; restoration_cell /= 3u;
        c->transcript_strand = restoration_cell % 2u ? 1 : -1;
        orientations = restoration_cell / 2u;
        memcpy(c->ref + 3u * codon, codons[repeated], 3u);
        memcpy(c->ref + 3u * (codon + 1u), codons[repeated], 3u);
        memcpy(c->ref + 3u * (codon + 2u), codons[replaced], 3u);
    }
    size_t codon_end = KPROP_HAPLO_CDS_LEN / 3u - (restoring ? 1u : 2u);
    while (c->edit_count < target && codon < codon_end) {
        size_t index = c->edit_count++;
        duckvep_haplotype_edit_t *e = c->edits + index;
        unsigned shape = restoring ? (unsigned)((index + 1u) % 3u) :
            (unsigned)kprop_bounded(t, 3u);
        e->cds_start = (uint32_t)(3u * codon + 1u);
        e->ref_len = shape == 0u ? 0u : 3u;
        e->alt_len = shape == 1u ? 0u : 3u;
        e->variant_strand = restoring ? ((orientations >> index) & 1u ? 1 : -1) :
            (kprop_bounded(t, 2u) ? 1 : -1);
        int reverse = e->variant_strand != c->transcript_strand;
        size_t alt_index = restoring ? (index == 1u ? repeated : replaced) :
            (size_t)kprop_bounded(t, 4u);
        if (shape == 2u && !memcmp(codons[alt_index], c->ref + 3u * codon, 3u))
            alt_index = (alt_index + 1u + (size_t)kprop_bounded(t, 3u)) % 4u;
        for (size_t b = 0u; b < 3u; b++) {
            size_t at = reverse ? 2u - b : b;
            c->ref_alleles[index][at] = haplo_test_variant_from_tx_base(
                (char)c->ref[3u * codon + b], reverse);
            c->alt_alleles[index][at] = haplo_test_variant_from_tx_base(
                (char)codons[alt_index][b], reverse);
        }
        e->ref = e->ref_len ? c->ref_alleles[index] : NULL;
        e->alt = e->alt_len ? c->alt_alleles[index] : NULL;
        codon += restoring ? (index == 0u ? 2u : 1u) : 1u + (size_t)kprop_bounded(t, 2u);
    }
    for (size_t i = 0u; i < c->edit_count / 2u; i++) {
        duckvep_haplotype_edit_t tmp = c->edits[i];
        c->edits[i] = c->edits[c->edit_count - 1u - i];
        c->edits[c->edit_count - 1u - i] = tmp;
    }
    *instance = c;
    return THEFT_ALLOC_OK;
}

static enum theft_alloc_res kprop_hgvs_inframe_haplotype_alloc(
    struct theft *t, void *env, void **instance) {
    (void)env;
    return kprop_hgvs_inframe_haplotype_alloc_mode(t, instance, SIZE_MAX);
}

enum { HGVS_RESTORATION_CELLS = 7 * 4 * 3 * 2 * 8 };
struct hgvs_restoration_quota { size_t next, cells[HGVS_RESTORATION_CELLS]; };

static enum theft_alloc_res kprop_hgvs_restoring_haplotype_alloc(
    struct theft *t, void *env, void **instance) {
    struct hgvs_restoration_quota *quota = env;
    size_t cell = quota->next++ % HGVS_RESTORATION_CELLS;
    enum theft_alloc_res status = kprop_hgvs_inframe_haplotype_alloc_mode(t, instance, cell);
    if (status == THEFT_ALLOC_OK) quota->cells[cell]++;
    return status;
}

static struct {
    size_t cases, merged, split, forward, reverse, restored, restored_changed_block, shapes[10];
} hgvs_haplotype_cov;

struct hgvs_separated_sampling { size_t attempts, rejected; };

static enum theft_alloc_res kprop_hgvs_separated_frame_alloc(
    struct theft *t, void *env, void **instance) {
    static const uint8_t bases[] = "ACGT";
    struct hgvs_separated_sampling *sampling = env;
    struct kprop_haplo_case *c = calloc(1u, sizeof *c);
    if (!c) return THEFT_ALLOC_ERROR;
    memcpy(c->ref, "ATG", 3u);
    for (size_t i = 3u; i < KPROP_HAPLO_CDS_LEN - 3u; i += 3u)
        memcpy(c->ref + i, kprop_bounded(t, 2u) ? "GCT" : "GGT", 3u);
    memcpy(c->ref + KPROP_HAPLO_CDS_LEN - 3u, "TAA", 3u);
    uint8_t ref[12], alt[12], rp[5], ap[5], inserted = 0u;
    int found = 0;
    /* Conditional domain: four translated codons, no stop, and at least two
     * differing runs separated by an unchanged residue after frame restoration.
     * Count every rejected proposal; generation never consults the HGVS builder. */
    for (size_t attempt = 0u; attempt < 4096u; attempt++) {
        sampling->attempts++;
        for (size_t i = 0u; i < sizeof ref; i++) ref[i] = bases[kprop_bounded(t, 4u)];
        inserted = bases[kprop_bounded(t, 4u)];
        alt[0] = inserted; memcpy(alt + 1u, ref, 9u); memcpy(alt + 10u, ref + 10u, 2u);
        size_t rn, an;
        if (!kprop_translate_full_oracle(ref, sizeof ref, DUCKVEP_CODON_TABLE_STANDARD, rp, &rn) ||
            !kprop_translate_full_oracle(alt, sizeof alt, DUCKVEP_CODON_TABLE_STANDARD, ap, &an)) {
            free(c); return THEFT_ALLOC_ERROR;
        }
        size_t runs = 0u; int changed = 0, stop = 0;
        for (size_t i = 0u; i < 4u; i++) {
            if (rp[i] == '*' || ap[i] == '*') stop = 1;
            if (rp[i] != ap[i] && !changed) runs++;
            changed = rp[i] != ap[i];
        }
        if (!stop && runs > 1u) { found = 1; break; }
        sampling->rejected++;
    }
    if (!found) { free(c); return THEFT_ALLOC_ERROR; }
    size_t start0 = 3u * (1u + (size_t)kprop_bounded(t, 7u));
    memcpy(c->ref + start0, ref, sizeof ref);
    c->transcript_strand = kprop_bounded(t, 2u) ? 1 : -1;
    c->edit_count = 2u;
    for (size_t i = 0u; i < 2u; i++) {
        duckvep_haplotype_edit_t *edit = c->edits + i;
        edit->variant_strand = kprop_bounded(t, 2u) ? 1 : -1;
        int reverse = edit->variant_strand != c->transcript_strand;
        edit->cds_start = (uint32_t)(start0 + (i ? 1u : 10u));
        if (i) {
            c->alt_alleles[i][0] = haplo_test_variant_from_tx_base((char)inserted, reverse);
            edit->alt_len = 1u; edit->alt = c->alt_alleles[i];
        } else {
            c->ref_alleles[i][0] = haplo_test_variant_from_tx_base((char)ref[9], reverse);
            edit->ref_len = 1u; edit->ref = c->ref_alleles[i];
        }
    }
    *instance = c;
    return THEFT_ALLOC_OK;
}

static enum theft_trial_res hgvs_haplotype_replays_complete_protein(
    const struct kprop_haplo_case *c, uint32_t curated_position, uint8_t curated_residue) {
    uint8_t cds[80], rp[32], ap[32], wanted_cds[80], wanted[32], replayed[2][80];
    size_t wanted_cds_length, wanted_length, count = 0u, block_count = 0u;
    int64_t difference; uint32_t flags;
    duckvep_edit_set_t set = {c->edits, c->edit_count};
    duckvep_coding_context_t ctx;
    duckvep_haplotype_edit_t ascending[KPROP_HAPLO_MAX_EDITS];
    duckvep_haplotype_block_t blocks[KPROP_HAPLO_MAX_EDITS];
    duckvep_hgvs_protein_operation_t operations[KPROP_HAPLO_MAX_EDITS + 2u];
    if (!haplo_oracle_rebuild(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
            c->transcript_strand, wanted_cds, sizeof wanted_cds, &wanted_cds_length,
            &difference, &flags) ||
        !kprop_translate_full_oracle(wanted_cds, wanted_cds_length,
            DUCKVEP_CODON_TABLE_STANDARD, wanted, &wanted_length)) return THEFT_TRIAL_ERROR;
    if (duckvep_coding_context_build(c->ref, KPROP_HAPLO_CDS_LEN, &set, c->transcript_strand,
            DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx) !=
            DUCKVEP_CODING_CONTEXT_OK || ctx.alt_cds_len != wanted_cds_length ||
        memcmp(cds, wanted_cds, wanted_cds_length) || ctx.alt_peptide_len != wanted_length ||
        memcmp(ap, wanted, wanted_length)) return THEFT_TRIAL_FAIL;
    ctx.post_cds_complete = 1u;
    if (curated_position) {
        if (curated_position > ctx.ref_peptide_len) return THEFT_TRIAL_ERROR;
        ctx.ref_peptide_edit_count = 1u;
        ctx.ref_peptide_edit_position1 = &curated_position;
        ctx.ref_peptide_edit_alt = &curated_residue;
    }
    for (size_t i = 0u; i < c->edit_count; i++) ascending[i] = c->edits[c->edit_count - 1u - i];
    if (duckvep_haplotype_partition(ascending, c->edit_count, blocks, KPROP_HAPLO_MAX_EDITS,
            &block_count) != DUCKVEP_HAPLOTYPE_OK) return THEFT_TRIAL_FAIL;
    duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_build(&ctx, NULL, ascending,
        c->edit_count, blocks, block_count, 0u, operations, KPROP_HAPLO_MAX_EDITS + 2u, &count);
    if (status != DUCKVEP_HGVS_OK) goto fail;
    int restored_changed_block = 0;
    if (!ctx.cds_changed) {
        for (size_t i = 0u; i < block_count; i++) {
            const duckvep_haplotype_block_t *block = blocks + i;
            if (block->ref_len == 0u || block->ref_len != block->alt_len ||
                !memcmp(c->ref + block->cds_start - 1u, cds + block->alt_start0,
                    block->ref_len)) continue;
            duckvep_sequence_delta_t delta;
            if (duckvep_coding_context_block_delta_fill(&ctx, ascending, c->edit_count,
                    block, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK || !delta.valid ||
                !delta.missense) goto fail;
            restored_changed_block = 1;
        }
    }
    memcpy(replayed[0], rp, ctx.ref_peptide_len);
    if (curated_position) replayed[0][curated_position - 1u] = curated_residue;
    size_t replayed_length = ctx.ref_peptide_len;
    unsigned side = 0u;
    for (size_t i = count; i > 0u; i--) {
        const duckvep_hgvs_protein_operation_t *operation = operations + i - 1u;
        int intersects;
        if (!operation->span.edit_count && (!curated_position ||
                operation->fact.shape != DUCKVEP_HGVS_PROTEIN_SUBSTITUTION ||
                operation->fact.first_position1 != curated_position)) goto fail;
        if ((operation->span.edit_count && duckvep_haplotype_block_frame_intersects(
                ascending, c->edit_count, &operation->span, 0u, 0u, &intersects) !=
                DUCKVEP_HAPLOTYPE_OK) ||
            (i < count && operation->span.edit_count && operations[i].span.edit_count &&
             operation->span.edit_begin + operation->span.edit_count >
                operations[i].span.edit_begin &&
                (operation->span.edit_begin != operations[i].span.edit_begin ||
                 operation->span.edit_count != operations[i].span.edit_count)) ||
            !kprop_hgvs_protein_fact_replay(&operation->fact, replayed[side], replayed_length,
                replayed[side ^ 1u], sizeof replayed[0], &replayed_length)) goto fail;
        side ^= 1u;
    }
    if (replayed_length != wanted_length || memcmp(replayed[side], wanted, wanted_length)) goto fail;
    if (!ctx.cds_changed && !curated_position && count) goto fail;
    char rendered[512], expected[512];
    size_t required, at = 0u;
    if (!count) {
        memcpy(expected, "p.(=)", 6u);
    } else {
        at = (size_t)snprintf(expected, sizeof expected, count > 1u ? "p.[(" : "p.");
        for (size_t i = 0u; i < count; i++) {
            char part[128]; size_t part_length;
            if (duckvep_hgvs_protein_render(&operations[i].fact, count == 1u, part, sizeof part,
                    &part_length) != DUCKVEP_HGVS_OK || part_length < 3u) goto fail;
            int written = snprintf(expected + at, sizeof expected - at, "%s%s",
                i ? ";" : "", part + 2u);
            if (written < 0 || (size_t)written >= sizeof expected - at) goto fail;
            at += (size_t)written;
        }
        if (count > 1u) {
            if (at + 2u >= sizeof expected) goto fail;
            expected[at++] = ')'; expected[at++] = ']'; expected[at] = '\0';
        }
    }
    if (duckvep_hgvs_protein_haplotype_render(operations, count, 1, rendered,
            sizeof rendered, &required) != DUCKVEP_HGVS_OK ||
        required != strlen(expected) || strcmp(rendered, expected)) goto fail;
    if (duckvep_hgvs_protein_haplotype_render(operations, count, 1, NULL, 0u,
            &required) != DUCKVEP_HGVS_BUFFER_TOO_SMALL || required != strlen(expected)) goto fail;
    if (curated_position) {
        uint8_t prepared[32];
        memcpy(prepared, rp, ctx.ref_peptide_len);
        prepared[curated_position - 1u] = curated_residue;
        duckvep_hgvs_protein_reference_t view = {prepared, ctx.ref_peptide_len};
        duckvep_hgvs_protein_operation_t borrowed[KPROP_HAPLO_MAX_EDITS + 2u];
        size_t borrowed_count;
        unsigned char saved[sizeof ctx];
        memcpy(saved, &ctx, sizeof ctx);
        if (duckvep_hgvs_protein_haplotype_build(&ctx, &view, ascending, c->edit_count,
                blocks, block_count, 0u, borrowed, KPROP_HAPLO_MAX_EDITS + 2u,
                &borrowed_count) != DUCKVEP_HGVS_OK || borrowed_count != count ||
            memcmp(saved, &ctx, sizeof ctx) ||
            duckvep_hgvs_protein_haplotype_render(borrowed, borrowed_count, 1,
                rendered, sizeof rendered, &required) != DUCKVEP_HGVS_OK ||
            strcmp(rendered, expected)) goto fail;
        memcpy(replayed[0], prepared, ctx.ref_peptide_len);
        replayed_length = ctx.ref_peptide_len; side = 0u;
        for (size_t i = borrowed_count; i > 0u; i--) {
            if (borrowed[i - 1u].fact.reference.bases != prepared ||
                !kprop_hgvs_block_equal(&borrowed[i - 1u].span, &operations[i - 1u].span) ||
                !kprop_hgvs_protein_fact_replay(&borrowed[i - 1u].fact,
                    replayed[side], replayed_length, replayed[side ^ 1u],
                    sizeof replayed[0], &replayed_length)) goto fail;
            side ^= 1u;
        }
        if (replayed_length != wanted_length || memcmp(replayed[side], wanted, wanted_length)) goto fail;
    }
    for (size_t i = 0u; i < count; i++) {
        if (operations[i].fact.shape >= 10u) return THEFT_TRIAL_FAIL;
        hgvs_haplotype_cov.shapes[operations[i].fact.shape]++;
        if (operations[i].span.edit_count > 1u) hgvs_haplotype_cov.merged++;
    }
    for (size_t i = 1u; i < count; i++) {
        if (operations[i - 1u].span.edit_begin == operations[i].span.edit_begin &&
            operations[i - 1u].span.edit_count == operations[i].span.edit_count) {
            hgvs_haplotype_cov.split++;
            break;
        }
    }
    hgvs_haplotype_cov.cases++;
    if (!ctx.cds_changed) hgvs_haplotype_cov.restored++;
    if (restored_changed_block) hgvs_haplotype_cov.restored_changed_block++;
    if (c->transcript_strand > 0) hgvs_haplotype_cov.forward++;
    else hgvs_haplotype_cov.reverse++;
    return THEFT_TRIAL_PASS;
fail:
    fprintf(stderr, "[compound HGVSp replay] status=%u strand=%d ref=%.*s alt=%.*s curated=%u:%c\n",
        (unsigned)status, (int)c->transcript_strand, (int)KPROP_HAPLO_CDS_LEN, c->ref,
        (int)wanted_cds_length, wanted_cds, curated_position,
        curated_position ? (char)curated_residue : '-');
    for (size_t i = 0u; i < c->edit_count; i++)
        fprintf(stderr, "  edit cds_start=%u ref=%.*s alt=%.*s variant_strand=%d\n",
            ascending[i].cds_start, (int)ascending[i].ref_len,
            ascending[i].ref ? (const char *)ascending[i].ref : "",
            (int)ascending[i].alt_len, ascending[i].alt ? (const char *)ascending[i].alt : "",
            (int)ascending[i].variant_strand);
    for (size_t i = 0u; i < count; i++) {
        char rendered[256]; size_t required;
        duckvep_hgvs_protein_render(&operations[i].fact, 0, rendered, sizeof rendered, &required);
        fprintf(stderr, "  operation %s edits=%zu+%zu\n", rendered,
            operations[i].span.edit_begin, operations[i].span.edit_count);
    }
    return THEFT_TRIAL_FAIL;
}

static enum theft_trial_res prop_hgvs_inframe_haplotype_replays_complete_protein(
    struct theft *t, void *arg) {
    (void)t;
    return hgvs_haplotype_replays_complete_protein(arg, 0u, 0u);
}

#define HGVS_CURATION_CELLS 80u
struct hgvs_curation_quota {
    size_t next, cells[HGVS_CURATION_CELLS];
    struct hgvs_separated_sampling sampling;
};
struct hgvs_curated_case {
    struct kprop_haplo_case *haplotype;
    uint32_t position1;
};
static enum theft_alloc_res kprop_hgvs_curated_alloc(struct theft *t, void *env, void **instance) {
    struct hgvs_curation_quota *quota = env;
    size_t cell = quota->next++ % HGVS_CURATION_CELLS;
    size_t route = cell / 20u;
    struct hgvs_curated_case *c = calloc(1u, sizeof *c);
    if (!c) return THEFT_ALLOC_ERROR;
    void *generated = NULL;
    enum theft_alloc_res status = route == 3u
        ? kprop_hgvs_separated_frame_alloc(t, &quota->sampling, &generated)
        : kprop_hgvs_inframe_haplotype_alloc_mode(t, &generated,
            route == 2u ? (size_t)kprop_bounded(t, HGVS_RESTORATION_CELLS) : SIZE_MAX);
    if (status != THEFT_ALLOC_OK) { free(c); return status; }
    c->haplotype = generated;
    c->position1 = (uint32_t)(2u + cell % 10u);
    int strand = (cell / 10u) % 2u ? 1 : -1;
    if (strand != c->haplotype->transcript_strand) {
        c->haplotype->transcript_strand *= -1;
        for (size_t i = 0u; i < c->haplotype->edit_count; i++)
            c->haplotype->edits[i].variant_strand *= -1;
    }
    if (!route) c->haplotype->edit_count = 0u;
    quota->cells[cell]++;
    *instance = c;
    return THEFT_ALLOC_OK;
}
static void kprop_hgvs_curated_free(void *instance, void *env) {
    (void)env;
    struct hgvs_curated_case *c = instance;
    free(c->haplotype);
    free(c);
}
static enum theft_trial_res prop_hgvs_curated_replay(struct theft *t, void *arg) {
    (void)t;
    const struct hgvs_curated_case *c = arg;
    return hgvs_haplotype_replays_complete_protein(c->haplotype, c->position1, 'U');
}

TEST hgvs_haplotype_curated_reference_replay(void) {
    struct hgvs_curation_quota quota = {0};
    struct theft_type_info info = {.alloc = kprop_hgvs_curated_alloc,
        .free = kprop_hgvs_curated_free, .env = &quota};
    struct theft_run_config cfg = {0};
    cfg.name = "compound HGVSp curated reference == complete protein replay";
    cfg.prop1 = prop_hgvs_curated_replay;
    cfg.type_info[0] = &info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    quota.next = cfg.seed % HGVS_CURATION_CELLS;
    memset(&hgvs_haplotype_cov, 0, sizeof hgvs_haplotype_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.cases);
    size_t minimum = SIZE_MAX, seen = 0u, total = 0u;
    for (size_t i = 0u; i < HGVS_CURATION_CELLS; i++) {
        if (quota.cells[i]) seen++;
        if (quota.cells[i] < minimum) minimum = quota.cells[i];
        total += quota.cells[i];
    }
    ASSERT_EQ(cfg.trials, total);
    ASSERT_EQ(cfg.trials / HGVS_CURATION_CELLS, minimum);
    ASSERT_EQ(cfg.trials < HGVS_CURATION_CELLS ? cfg.trials : HGVS_CURATION_CELLS, seen);
    fprintf(stderr, "[compound HGVSp curated-reference coverage] cases=%zu cells=%zu min_per_cell=%zu "
        "forward=%zu reverse=%zu attempts=%zu rejected=%zu prepared_views=%zu\n", hgvs_haplotype_cov.cases,
        seen, minimum, hgvs_haplotype_cov.forward, hgvs_haplotype_cov.reverse,
        quota.sampling.attempts, quota.sampling.rejected, hgvs_haplotype_cov.cases);
    PASS();
}

#define HGVS_TERMINAL_REFERENCE_CELLS (3u * 21u * 2u * 4u)
#define HGVS_TERMINAL_REPEAT_CELLS (20u * 20u * 3u * 2u * 2u)
struct hgvs_terminal_reference_case {
    struct kprop_haplo_case *haplotype;
    uint8_t residue;
    uint8_t inserted[6];
    int complete_replay;
};
struct hgvs_terminal_reference_quota {
    size_t next, cases, passed, forward, reverse, cells[HGVS_TERMINAL_REPEAT_CELLS];
    struct hgvs_separated_sampling sampling;
};

static enum theft_alloc_res kprop_hgvs_terminal_reference_alloc(
    struct theft *t, void *env, void **instance) {
    static const uint8_t residues[] = "ACDEFGHIKLMNPQRSTVWYU";
    static const uint8_t stops[3][3] = {{'T','A','A'}, {'T','A','G'}, {'T','G','A'}};
    struct hgvs_terminal_reference_quota *quota = env;
    size_t cell = quota->next++ % HGVS_TERMINAL_REFERENCE_CELLS;
    size_t route = cell / (3u * 21u * 2u);
    struct hgvs_terminal_reference_case *c = calloc(1u, sizeof *c);
    if (!c) return THEFT_ALLOC_ERROR;
    void *generated = NULL;
    enum theft_alloc_res status = route == 3u
        ? kprop_hgvs_separated_frame_alloc(t, &quota->sampling, &generated)
        : kprop_hgvs_inframe_haplotype_alloc_mode(t, &generated,
            route == 2u ? (size_t)kprop_bounded(t, HGVS_RESTORATION_CELLS) : SIZE_MAX);
    if (status != THEFT_ALLOC_OK) { free(c); return status; }
    c->haplotype = generated;
    c->residue = residues[(cell / 3u) % 21u];
    memcpy(c->haplotype->ref + KPROP_HAPLO_CDS_LEN - 3u, stops[cell % 3u], 3u);
    int strand = (cell / (3u * 21u)) % 2u ? 1 : -1;
    if (strand != c->haplotype->transcript_strand) {
        c->haplotype->transcript_strand *= -1;
        for (size_t i = 0u; i < c->haplotype->edit_count; i++)
            c->haplotype->edits[i].variant_strand *= -1;
    }
    if (!route) c->haplotype->edit_count = 0u;
    quota->cells[cell]++;
    *instance = c;
    return THEFT_ALLOC_OK;
}

static enum theft_alloc_res kprop_hgvs_terminal_repeat_alloc(
    struct theft *t, void *env, void **instance) {
    static const uint8_t residues[] = "ACDEFGHIKLMNPQRSTVWY";
    static const char *codons[] = {"GCT","TGT","GAT","GAA","TTT","GGT","CAT","ATT",
        "AAA","CTG","ATG","AAT","CCT","CAA","CGT","TCT","ACT","GTT","TGG","TAT"};
    static const char *stops[] = {"TAA","TAG","TGA"};
    struct hgvs_terminal_reference_quota *quota = env;
    size_t cell = quota->next++ % HGVS_TERMINAL_REPEAT_CELLS;
    size_t x = cell % 20u, y = (cell / 20u) % 20u;
    struct hgvs_terminal_reference_case *input = calloc(1u, sizeof *input);
    if (!input) return THEFT_ALLOC_ERROR;
    struct kprop_haplo_case *c = calloc(1u, sizeof *c);
    if (!c) { free(input); return THEFT_ALLOC_ERROR; }
    input->haplotype = c;
    input->residue = residues[x];
    input->complete_replay = 1;
    size_t copies = 1u + kprop_bounded(t, 5u), prefix = 10u - 2u * copies;
    memcpy(c->ref, "ATG", 3u);
    for (size_t i = 0u; i < prefix; i++)
        memcpy(c->ref + 3u + 3u * i, codons[kprop_bounded(t, 20u)], 3u);
    uint8_t motif[6];
    memcpy(motif, codons[x], 3u); memcpy(motif + 3u, codons[y], 3u);
    for (size_t i = 0u; i < copies; i++) memcpy(c->ref + 3u + 3u * prefix + 6u * i, motif, 6u);
    memcpy(c->ref + KPROP_HAPLO_CDS_LEN - 3u, stops[(cell / 400u) % 3u], 3u);
    c->transcript_strand = (cell / 1200u) % 2u ? 1 : -1;
    int8_t variant_strand = (cell / 2400u) % 2u ? 1 : -1;
    for (size_t i = 0u; i < 6u; i++) input->inserted[i] = (uint8_t)haplo_test_oriented_base(
        motif, 6u, (uint32_t)i, variant_strand != c->transcript_strand);
    c->edit_count = 1u;
    c->edits[0] = (duckvep_haplotype_edit_t){
        (uint32_t)(4u + 3u * prefix + 6u * kprop_bounded(t, copies + 1u)),
        0u, NULL, 6u, input->inserted, variant_strand};
    quota->cells[cell]++;
    *instance = input;
    return THEFT_ALLOC_OK;
}

static void kprop_hgvs_terminal_reference_free(void *instance, void *env) {
    (void)env;
    struct hgvs_terminal_reference_case *c = instance;
    free(c->haplotype);
    free(c);
}

static enum theft_trial_res prop_hgvs_terminal_reference_replay(struct theft *t, void *arg) {
    struct hgvs_terminal_reference_quota *quota = theft_hook_get_env(t);
    const struct hgvs_terminal_reference_case *input = arg;
    const struct kprop_haplo_case *c = input->haplotype;
    quota->cases++;
    if (c->transcript_strand > 0) quota->forward++;
    else quota->reverse++;
    uint8_t cds[80], rp[32], ap[32], prepared[32], wanted_reference[32];
    uint8_t coding_reference[32];
    duckvep_translation_t coding_translation;
    uint8_t wanted_cds[80], wanted[32], replayed[2][80];
    size_t wanted_cds_length, wanted_length, reference_length, prepared_length;
    int64_t difference; uint32_t flags;
    if (!haplo_oracle_rebuild(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
            c->transcript_strand, wanted_cds, sizeof wanted_cds, &wanted_cds_length,
            &difference, &flags) ||
        !kprop_translate_full_oracle(wanted_cds, wanted_cds_length,
            DUCKVEP_CODON_TABLE_STANDARD, wanted, &wanted_length) ||
        !kprop_translate_full_oracle(c->ref, KPROP_HAPLO_CDS_LEN,
            DUCKVEP_CODON_TABLE_STANDARD, wanted_reference, &reference_length)) return THEFT_TRIAL_ERROR;
    uint32_t position = (uint32_t)reference_length;
    if (!reference_length || wanted_reference[reference_length - 1u] != '*') return THEFT_TRIAL_ERROR;
    wanted_reference[reference_length - 1u] = input->residue;
    wanted_reference[reference_length++] = '*';
    if (duckvep_haplotype_reference_proteins(c->ref, KPROP_HAPLO_CDS_LEN,
            DUCKVEP_CODON_TABLE_STANDARD, &position, &input->residue, 1u,
            prepared, coding_reference, sizeof prepared, &prepared_length, &coding_translation) !=
            DUCKVEP_HAPLOTYPE_OK ||
        prepared_length != reference_length || memcmp(prepared, wanted_reference, reference_length))
        return THEFT_TRIAL_FAIL;
    duckvep_edit_set_t set = {c->edits, c->edit_count};
    duckvep_coding_context_t ctx;
    if (duckvep_coding_context_build(c->ref, KPROP_HAPLO_CDS_LEN, &set, c->transcript_strand,
            DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx) !=
            DUCKVEP_CODING_CONTEXT_OK || ctx.alt_cds_len != wanted_cds_length ||
        memcmp(cds, wanted_cds, wanted_cds_length) || ctx.alt_peptide_len != wanted_length ||
        memcmp(ap, wanted, wanted_length)) return THEFT_TRIAL_FAIL;
    ctx.post_cds_complete = 1u;
    ctx.ref_peptide_edit_count = 1u;
    ctx.ref_peptide_edit_position1 = &position;
    ctx.ref_peptide_edit_alt = &input->residue;
    duckvep_haplotype_edit_t ascending[KPROP_HAPLO_MAX_EDITS];
    for (size_t i = 0u; i < c->edit_count; i++) ascending[i] = c->edits[c->edit_count - 1u - i];
    duckvep_haplotype_block_t blocks[KPROP_HAPLO_MAX_EDITS];
    size_t block_count, count = 0u;
    if (duckvep_haplotype_partition(ascending, c->edit_count, blocks, KPROP_HAPLO_MAX_EDITS,
            &block_count) != DUCKVEP_HAPLOTYPE_OK) return THEFT_TRIAL_FAIL;
    unsigned char saved[sizeof ctx];
    memcpy(saved, &ctx, sizeof ctx);
    duckvep_hgvs_protein_reference_t view = {prepared, prepared_length};
    duckvep_hgvs_protein_operation_t operations[KPROP_HAPLO_MAX_EDITS + 4u];
    duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_build(&ctx, &view, ascending,
        c->edit_count, blocks, block_count, 0u, operations, KPROP_HAPLO_MAX_EDITS + 4u, &count);
    if (status != DUCKVEP_HGVS_OK || memcmp(saved, &ctx, sizeof ctx) ||
        memcmp(prepared, wanted_reference, prepared_length)) goto fail;
    memcpy(replayed[0], prepared, prepared_length);
    size_t replayed_length = prepared_length;
    unsigned side = 0u;
    for (size_t i = count; i > 0u; i--) {
        const duckvep_hgvs_protein_operation_t *operation = operations + i - 1u;
        int intersects;
        if (!operation->span.edit_count && operation->fact.first_position1 != position) goto fail;
        if (operation->span.edit_count && duckvep_haplotype_block_frame_intersects(ascending,
                c->edit_count, &operation->span, 0u, 0u, &intersects) != DUCKVEP_HAPLOTYPE_OK) goto fail;
        if (!kprop_hgvs_protein_fact_replay(&operation->fact, replayed[side], replayed_length,
                replayed[side ^ 1u], sizeof replayed[0], &replayed_length)) goto fail;
        side ^= 1u;
    }
    /* These protein operands use the displayed first-stop prefix. Keep the
     * complete replay above, then apply that presentation rule on both sides. */
    if (input->complete_replay && (wanted_length != replayed_length ||
            memcmp(wanted, replayed[side], wanted_length))) goto fail;
    for (size_t i = 0u; i < wanted_length; i++) if (wanted[i] == '*') { wanted_length = i + 1u; break; }
    for (size_t i = 0u; i < replayed_length; i++)
        if (replayed[side][i] == '*') { replayed_length = i + 1u; break; }
    if (wanted_length != replayed_length || memcmp(wanted, replayed[side], wanted_length)) goto fail;
    char text[512]; size_t required;
    if (duckvep_hgvs_protein_haplotype_render(operations, count, 1, text, sizeof text, &required) !=
            DUCKVEP_HGVS_OK || !required) goto fail;
    quota->passed++;
    return THEFT_TRIAL_PASS;
fail:
    fprintf(stderr, "[compound HGVSp terminal reference] status=%u ref=%.*s alt=%.*s curated=%c\n",
        (unsigned)status, (int)KPROP_HAPLO_CDS_LEN, c->ref, (int)wanted_cds_length,
        wanted_cds, input->residue);
    for (size_t i = 0u; i < c->edit_count; i++)
        fprintf(stderr, "  edit cds_start=%u ref=%.*s alt=%.*s strand=%d transcript_strand=%d\n",
            ascending[i].cds_start, (int)ascending[i].ref_len,
            ascending[i].ref ? (const char *)ascending[i].ref : "",
            (int)ascending[i].alt_len, ascending[i].alt ? (const char *)ascending[i].alt : "",
            (int)ascending[i].variant_strand, (int)c->transcript_strand);
    return THEFT_TRIAL_FAIL;
}

static int kprop_hgvs_terminal_reference_trials(int repeats) {
    struct hgvs_terminal_reference_quota quota = {0};
    size_t cell_count = repeats ? HGVS_TERMINAL_REPEAT_CELLS : HGVS_TERMINAL_REFERENCE_CELLS;
    struct theft_type_info info = {.alloc = repeats ? kprop_hgvs_terminal_repeat_alloc
            : kprop_hgvs_terminal_reference_alloc,
        .free = kprop_hgvs_terminal_reference_free, .env = &quota};
    struct theft_run_config cfg = {0};
    cfg.name = repeats ? "compound HGVSp terminal repeat == complete protein replay"
        : "compound HGVSp terminal reference == displayed protein replay";
    cfg.prop1 = prop_hgvs_terminal_reference_replay;
    cfg.type_info[0] = &info;
    cfg.hooks.env = &quota;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    quota.next = cfg.seed % cell_count;
    enum theft_run_res result = theft_run(&cfg);
    size_t minimum = SIZE_MAX, seen = 0u, total = 0u;
    for (size_t i = 0u; i < cell_count; i++) {
        if (quota.cells[i]) seen++;
        if (quota.cells[i] < minimum) minimum = quota.cells[i];
        total += quota.cells[i];
    }
    fprintf(stderr, "[compound HGVSp %s coverage] cases=%zu passed=%zu generated=%zu "
        "cells=%zu min_per_cell=%zu "
        "forward=%zu reverse=%zu",
        repeats ? "terminal-repeat" : "terminal-reference", quota.cases, quota.passed, total, seen, minimum,
        quota.forward, quota.reverse);
    if (!repeats) fprintf(stderr, " attempts=%zu rejected=%zu", quota.sampling.attempts, quota.sampling.rejected);
    fputc('\n', stderr);
    ASSERT_EQ(THEFT_RUN_PASS, result);
    ASSERT_EQ(cfg.trials, quota.cases);
    ASSERT_EQ(cfg.trials, quota.passed);
    ASSERT_EQ(cfg.trials, total);
    ASSERT_EQ(cfg.trials / cell_count, minimum);
    ASSERT_EQ(cfg.trials < cell_count ? cfg.trials : cell_count, seen);
    PASS();
}

TEST hgvs_haplotype_terminal_reference_replay(void) {
    return kprop_hgvs_terminal_reference_trials(0);
}

TEST hgvs_haplotype_terminal_repeat_insertions_replay(void) {
    return kprop_hgvs_terminal_reference_trials(1);
}

TEST hgvs_haplotype_inframe_operations_replay_complete_protein(void) {
    struct theft_type_info info = {.alloc = kprop_hgvs_inframe_haplotype_alloc, .free = kprop_haplo_free};
    struct theft_run_config cfg = {0};
    cfg.name = "compound HGVSp operations == literal in-frame CDS replay and independent translation";
    cfg.prop1 = prop_hgvs_inframe_haplotype_replays_complete_protein;
    cfg.type_info[0] = &info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&hgvs_haplotype_cov, 0, sizeof hgvs_haplotype_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.cases);
    ASSERT(hgvs_haplotype_cov.merged && hgvs_haplotype_cov.forward && hgvs_haplotype_cov.reverse);
    for (unsigned shape = DUCKVEP_HGVS_PROTEIN_SUBSTITUTION;
         shape <= DUCKVEP_HGVS_PROTEIN_DUPLICATION; shape++) ASSERT(hgvs_haplotype_cov.shapes[shape]);
    fprintf(stderr, "[compound HGVSp replay coverage] cases=%zu merged=%zu forward=%zu reverse=%zu "
        "sub=%zu del=%zu ins=%zu delins=%zu dup=%zu\n", hgvs_haplotype_cov.cases,
        hgvs_haplotype_cov.merged, hgvs_haplotype_cov.forward, hgvs_haplotype_cov.reverse,
        hgvs_haplotype_cov.shapes[DUCKVEP_HGVS_PROTEIN_SUBSTITUTION],
        hgvs_haplotype_cov.shapes[DUCKVEP_HGVS_PROTEIN_DELETION],
        hgvs_haplotype_cov.shapes[DUCKVEP_HGVS_PROTEIN_INSERTION],
        hgvs_haplotype_cov.shapes[DUCKVEP_HGVS_PROTEIN_DELINS],
        hgvs_haplotype_cov.shapes[DUCKVEP_HGVS_PROTEIN_DUPLICATION]);
    PASS();
}

TEST hgvs_haplotype_restored_cds_keeps_local_changes(void) {
    /* Seven coding locations x twelve distinct codon pairs x two transcript
     * strands x eight physical-allele orientations. Randomize flanking CDS
     * within every cell; no cell depends on a chance occurrence of restoration. */
    struct hgvs_restoration_quota quota = {0};
    struct theft_type_info info = {.alloc = kprop_hgvs_restoring_haplotype_alloc,
        .free = kprop_haplo_free, .env = &quota};
    struct theft_run_config cfg = {0};
    cfg.name = "compound HGVSp restored CDS == complete replay with changed local blocks";
    cfg.prop1 = prop_hgvs_inframe_haplotype_replays_complete_protein;
    cfg.type_info[0] = &info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    quota.next = cfg.seed % HGVS_RESTORATION_CELLS;
    memset(&hgvs_haplotype_cov, 0, sizeof hgvs_haplotype_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.cases);
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.restored);
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.restored_changed_block);
    ASSERT(hgvs_haplotype_cov.forward && hgvs_haplotype_cov.reverse);
    size_t minimum = SIZE_MAX, cells_seen = 0u, total = 0u;
    for (size_t i = 0u; i < HGVS_RESTORATION_CELLS; i++) {
        if (quota.cells[i]) cells_seen++;
        if (quota.cells[i] < minimum) minimum = quota.cells[i];
        total += quota.cells[i];
    }
    ASSERT_EQ(cfg.trials, total);
    ASSERT_EQ(cfg.trials / HGVS_RESTORATION_CELLS, minimum);
    ASSERT_EQ(cfg.trials < HGVS_RESTORATION_CELLS ? cfg.trials : HGVS_RESTORATION_CELLS, cells_seen);
    fprintf(stderr, "[compound HGVSp restored-CDS coverage] cases=%zu restored=%zu "
        "changed_block=%zu forward=%zu reverse=%zu cells=%zu min_per_cell=%zu\n",
        hgvs_haplotype_cov.cases,
        hgvs_haplotype_cov.restored, hgvs_haplotype_cov.restored_changed_block,
        hgvs_haplotype_cov.forward, hgvs_haplotype_cov.reverse, cells_seen, minimum);
    PASS();
}

TEST hgvs_haplotype_separated_restored_frame_changes_replay(void) {
    struct hgvs_separated_sampling sampling = {0};
    struct theft_type_info info = {.alloc = kprop_hgvs_separated_frame_alloc,
        .free = kprop_haplo_free, .env = &sampling};
    struct theft_run_config cfg = {0};
    cfg.name = "compound HGVSp separated restored-frame changes == complete protein replay";
    cfg.prop1 = prop_hgvs_inframe_haplotype_replays_complete_protein;
    cfg.type_info[0] = &info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&hgvs_haplotype_cov, 0, sizeof hgvs_haplotype_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.cases);
    ASSERT_EQ(cfg.trials, sampling.attempts - sampling.rejected);
    ASSERT_EQ(cfg.trials, hgvs_haplotype_cov.split);
    ASSERT(hgvs_haplotype_cov.forward && hgvs_haplotype_cov.reverse);
    fprintf(stderr, "[compound HGVSp separated-frame coverage] cases=%zu attempts=%zu "
        "rejected=%zu split=%zu forward=%zu reverse=%zu\n", hgvs_haplotype_cov.cases,
        sampling.attempts, sampling.rejected, hgvs_haplotype_cov.split,
        hgvs_haplotype_cov.forward, hgvs_haplotype_cov.reverse);
    PASS();
}

TEST hgvs_haplotype_checks_complete_inputs_before_output(void) {
    static const uint8_t reference[] = "ATGGGTCCTGCTGAACAATAA";
    duckvep_haplotype_edit_t ascending[2] = {
        {10u, 3u, reference + 9u, 3u, (const uint8_t *)"TAA", 1},
        {16u, 3u, reference + 15u, 3u, (const uint8_t *)"GAA", 1}
    };
    duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
    duckvep_edit_set_t set = {descending, 2u};
    uint8_t cds[32], rp[16], ap[16];
    duckvep_coding_context_t ctx;
    duckvep_haplotype_block_t blocks[2];
    duckvep_hgvs_protein_operation_t facts[2];
    size_t blocks_used, facts_used;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_EQ(4u, ctx.alt_first_stop_position1);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &blocks_used));
    ASSERT_EQ(2u, blocks_used);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
        &ctx, NULL, ascending, 2u, blocks, 2u, 0u, facts, 2u, &facts_used));
    ASSERT_EQ(1u, facts_used);
    char rendered[32]; size_t required;
    ASSERT_EQ(DUCKVEP_HGVS_OK,
        duckvep_hgvs_protein_render(&facts[0].fact, 0, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.Ala4Ter", rendered);
    for (unsigned bad = 0u; bad < 21u; bad++) {
        duckvep_coding_context_t changed = ctx;
        duckvep_haplotype_block_t corrupt[2] = {blocks[0], blocks[1]};
        const duckvep_haplotype_edit_t *edits = ascending;
        const duckvep_haplotype_block_t *input_blocks = corrupt;
        size_t input_count = 2u, capacity = 2u;
        duckvep_hgvs_status_t expected = DUCKVEP_HGVS_INVALID_ARG;
        switch (bad) {
            case 0: changed.ref_cds = NULL; break;
            case 1: changed.alt_cds = NULL; break;
            case 2: changed.ref_peptide = NULL; break;
            case 3: changed.alt_peptide = NULL; break;
            case 4: changed.virtual_single_edit = 1u; break;
            case 5: changed.compatibility_profile = UINT8_MAX; break;
            case 6: changed.codon_table = UINT8_MAX; break;
            case 7: changed.ref_peptide_len++; break;
            case 8: changed.alt_peptide_len++; break;
            case 9: changed.alt_first_stop_position1 = 2u; break;
            case 10: changed.length_diff++; break;
            case 11: changed.applied_edits++; break;
            case 12: edits = NULL; break;
            case 13: input_blocks = NULL; break;
            case 14: capacity = 1u; expected = DUCKVEP_HGVS_BUFFER_TOO_SMALL; break;
            case 15: corrupt[1].flags = DUCKVEP_HAPLOTYPE_FLAG_INDEL;
                expected = DUCKVEP_HGVS_INVALID_PROJECTION; break;
            case 16: corrupt[1].alt_start0++;
                expected = DUCKVEP_HGVS_INVALID_PROJECTION; break;
            case 17: corrupt[1].edit_begin = 0u; break;
            case 18: input_count = 1u; expected = DUCKVEP_HGVS_INVALID_PROJECTION; break;
            case 19: corrupt[0].ref_len = 9u;
                expected = DUCKVEP_HGVS_INVALID_PROJECTION; break;
            case 20: corrupt[0].alt_len = 9u;
                expected = DUCKVEP_HGVS_INVALID_PROJECTION; break;
        }
        memset(facts, 0xa5, sizeof facts);
        duckvep_hgvs_protein_operation_t saved[2]; memcpy(saved, facts, sizeof saved);
        facts_used = SIZE_MAX;
        ASSERT_EQ(expected, duckvep_hgvs_protein_haplotype_build(&changed, NULL, edits, 2u,
            input_blocks, input_count, 0u, facts, capacity, &facts_used));
        ASSERT_EQ(0u, facts_used);
        ASSERT_MEM_EQ(saved, facts, sizeof facts);
    }
    set.count = 0u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_EQ(DUCKVEP_HGVS_OK,
        duckvep_hgvs_protein_haplotype_build(&ctx, NULL, NULL, 0u, NULL, 0u, 0u,
            NULL, 0u, &facts_used));
    ASSERT_EQ(0u, facts_used);
    uint32_t position1 = 2u; uint8_t curated = 'U';
    ctx.ref_peptide_edit_position1 = &position1;
    ctx.ref_peptide_edit_alt = &curated;
    ctx.ref_peptide_edit_count = 1u;
    ASSERT_EQ(DUCKVEP_HGVS_BUFFER_TOO_SMALL,
        duckvep_hgvs_protein_haplotype_build(&ctx, NULL, NULL, 0u, NULL, 0u, 0u,
            NULL, 0u, &facts_used));
    ASSERT_EQ(0u, facts_used);
    ASSERT_EQ(DUCKVEP_HGVS_OK,
        duckvep_hgvs_protein_haplotype_build(&ctx, NULL, NULL, 0u, NULL, 0u, 0u,
            facts, 2u, &facts_used));
    ASSERT_EQ(1u, facts_used);
    ASSERT_EQ(0u, facts[0].span.edit_count);
    ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
        facts, facts_used, 0, rendered, sizeof rendered, &required));
    ASSERT_STR_EQ("p.Sec2Gly", rendered);
    for (unsigned bad = 0u; bad < 6u; bad++) {
        duckvep_coding_context_t corrupt = ctx;
        uint32_t positions[2] = {2u, 2u};
        uint8_t residues[2] = {'U', 'W'};
        corrupt.ref_peptide_edit_position1 = positions;
        corrupt.ref_peptide_edit_alt = residues;
        switch (bad) {
            case 0: corrupt.ref_peptide_edit_position1 = NULL; break;
            case 1: corrupt.ref_peptide_edit_alt = NULL; break;
            case 2: corrupt.ref_peptide_edit_count = ctx.ref_peptide_len + 1u; break;
            case 3: positions[0] = 0u; break;
            case 4: positions[0] = (uint32_t)ctx.ref_peptide_len + 1u; break;
            case 5: corrupt.ref_peptide_edit_count = 2u; break;
        }
        facts_used = SIZE_MAX;
        ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
            duckvep_hgvs_protein_haplotype_build(&corrupt, NULL, NULL, 0u, NULL, 0u, 0u,
                facts, 2u, &facts_used));
        ASSERT_EQ(0u, facts_used);
    }
    ctx.compatibility_profile = UINT8_MAX;
    ASSERT_EQ(DUCKVEP_HGVS_INVALID_ARG,
        duckvep_hgvs_protein_haplotype_build(&ctx, NULL, NULL, 0u, NULL, 0u, 0u,
            NULL, 0u, &facts_used));
    ASSERT_EQ(0u, facts_used);
    PASS();
}

TEST hgvs_haplotype_shifted_frameshift_codon_landscape(void) {
    static const uint8_t dna[] = "ACGT";
    static const uint8_t tail[] = "TAACTAGCTGA";
    static const uint32_t positions1[] = {10u, 11u, 12u, 9u};
    size_t cases = 0u, equal = 0u, immediate_stop = 0u, frameshift = 0u, stop_window = 0u;
    size_t anchored = 0u, anchor_presentation_differences = 0u;
    for (unsigned pair = 0u; pair < 4096u; pair++) {
        uint8_t reference[] = "ATGGGTCCTAAAAAAGAACAATAA";
        for (unsigned b = 0u; b < 6u; b++)
            reference[9u + b] = dna[(pair >> (2u * b)) & 3u];
        for (unsigned phase = 0u; phase < 4u; phase++) {
            for (unsigned base = 0u; base < 4u; base++) {
                for (int shift = -1; shift <= 1; shift += 2) {
                    for (int strand = -1; strand <= 1; strand += 2) {
                        uint8_t ref[3], ins[3], inserted = dna[base];
                        for (unsigned b = 0u; b < 3u; b++) {
                            ref[b] = strand > 0 ? reference[3u + b] :
                                (uint8_t)kprop_complement_base((char)reference[5u - b]);
                            ins[b] = strand > 0 ? 'A' : 'T';
                        }
                        if (strand < 0) inserted = (uint8_t)kprop_complement_base((char)inserted);
                        duckvep_haplotype_edit_t ascending[2] = {
                            {4u, shift < 0 ? 3u : 0u, ref,
                             shift > 0 ? 3u : 0u, ins, 1},
                            {positions1[phase], 0u, NULL, 1u, &inserted, 1}
                        };
                        duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
                        duckvep_edit_set_t set = {descending, 2u};
                        uint8_t cds[40], rp[16], ap[16], extended[64], want_ref[16], want_alt[24];
                        size_t ref_length = 0u, alt_length = 0u, blocks_used = 0u, facts_used = 0u;
                        duckvep_coding_context_t ctx;
                        duckvep_haplotype_block_t blocks[2];
                        duckvep_hgvs_protein_operation_t facts[2];
                        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                            reference, sizeof reference - 1u, &set, (int8_t)strand,
                            DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds,
                            rp, sizeof rp, ap, sizeof ap, &ctx));
                        ctx.pre_cds_complete = ctx.post_cds_complete = 1u;
                        ctx.post_cds_bases = tail;
                        ctx.post_cds_length = sizeof tail - 1u;
                        size_t written = 0u;
                        for (size_t b = 0u; b < sizeof reference - 1u; b++) {
                            if (shift < 0 && b >= 3u && b < 6u) continue;
                            if (shift > 0 && b == 3u) {
                                memcpy(extended + written, "AAA", 3u); written += 3u;
                            }
                            if (b + 1u == positions1[phase]) extended[written++] = dna[base];
                            extended[written++] = reference[b];
                        }
                        ASSERT_EQ(written, ctx.alt_cds_len);
                        ASSERT_MEM_EQ(extended, cds, written);
                        memcpy(extended + written, tail, sizeof tail - 1u);
                        ASSERT(kprop_translate_full_oracle(reference, sizeof reference - 1u,
                            DUCKVEP_CODON_TABLE_STANDARD, want_ref, &ref_length));
                        ASSERT(kprop_translate_full_oracle(extended,
                            ctx.alt_cds_len + sizeof tail - 1u, DUCKVEP_CODON_TABLE_STANDARD,
                            want_alt, &alt_length));
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                            duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &blocks_used));
                        ASSERT_EQ(phase == 3u && shift < 0 ? 1u : 2u, blocks_used);
                        uint8_t single_cds[40], single_rp[16], single_ap[16];
                        duckvep_coding_context_t single;
                        duckvep_sequence_delta_t delta;
                        duckvep_hgvs_protein_fact_t expected;
                        duckvep_edit_set_t one_edit = {ascending + 1u, 1u};
                        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                            reference, sizeof reference - 1u, &one_edit, (int8_t)strand,
                            DUCKVEP_CODON_TABLE_STANDARD, single_cds, sizeof single_cds,
                            single_rp, sizeof single_rp, single_ap, sizeof single_ap, &single));
                        single.pre_cds_complete = single.post_cds_complete = 1u;
                        single.post_cds_bases = tail;
                        single.post_cds_length = sizeof tail - 1u;
                        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                            duckvep_coding_context_delta_fill(&single, 0u, &delta));
                        ASSERT_EQ(DUCKVEP_HGVS_OK,
                            duckvep_hgvs_protein_fact_build(&single, &delta, &expected));
                        duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_build(
                            &ctx, NULL, ascending, 2u, blocks, blocks_used, 0u, facts, 2u, &facts_used);
                        size_t r = 3u, a = (size_t)(3 + shift);
                        while (r < ref_length && a < alt_length &&
                               want_ref[r] == want_alt[a] && want_ref[r] != '*') { r++; a++; }
                        if (status != DUCKVEP_HGVS_OK) {
                            fprintf(stderr, "[compound HGVSp status] pair=%u phase=%u base=%u "
                                "shift=%d strand=%d status=%u ref=%.*s alt=%.*s\n", pair,
                                phase, base, shift, strand, (unsigned)status,
                                (int)ref_length, want_ref, (int)alt_length, want_alt);
                        }
                        ASSERT_EQ(DUCKVEP_HGVS_OK, status);
                        /* Retained anchors preserve the completed CDS. VEP-116
                         * protein presentation instead follows each source edit's
                         * local peptide window; see ERRATA.md. Check every raw
                         * operation against its own isolated edit, including all
                         * representations whose HGVS differs from the minimal edit. */
                        char minimal_text[128]; size_t minimal_required;
                        ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
                            facts, facts_used, 1, minimal_text, sizeof minimal_text, &minimal_required));
                        for (unsigned right = 0u; right < 2u; right++) {
                            uint8_t anchored_ref[2][4], anchored_alt[2][4];
                            duckvep_haplotype_edit_t raw[2] = {
                                {3u + right, shift < 0 ? 4u : 1u, anchored_ref[0],
                                 shift > 0 ? 4u : 1u, anchored_alt[0], 1},
                                {positions1[phase] - 1u + right, 1u, anchored_ref[1],
                                 2u, anchored_alt[1], 1}
                            };
                            for (size_t e = 0u; e < 2u; e++) {
                                memcpy(anchored_ref[e], reference + raw[e].cds_start - 1u,
                                    raw[e].ref_len);
                                size_t anchor = right ? raw[e].ref_len - 1u : 0u;
                                anchored_alt[e][right ? raw[e].alt_len - 1u : 0u] =
                                    anchored_ref[e][anchor];
                                if (!e && shift > 0) memcpy(anchored_alt[e] + !right, "AAA", 3u);
                                if (e) anchored_alt[e][!right] = dna[base];
                                if (strand < 0) {
                                    uint8_t saved_ref[4], saved_alt[4];
                                    memcpy(saved_ref, anchored_ref[e], raw[e].ref_len);
                                    memcpy(saved_alt, anchored_alt[e], raw[e].alt_len);
                                    for (size_t b = 0u; b < raw[e].ref_len; b++)
                                        anchored_ref[e][b] = (uint8_t)kprop_complement_base(
                                            (char)saved_ref[raw[e].ref_len - 1u - b]);
                                    for (size_t b = 0u; b < raw[e].alt_len; b++)
                                        anchored_alt[e][b] = (uint8_t)kprop_complement_base(
                                            (char)saved_alt[raw[e].alt_len - 1u - b]);
                                }
                            }
                            duckvep_haplotype_edit_t raw_descending[2] = {raw[1], raw[0]};
                            duckvep_edit_set_t raw_set = {raw_descending, 2u};
                            uint8_t raw_cds[40], raw_rp[16], raw_ap[16];
                            duckvep_coding_context_t raw_context;
                            duckvep_haplotype_block_t raw_blocks[2];
                            duckvep_hgvs_protein_operation_t raw_facts[2];
                            size_t raw_blocks_used, raw_facts_used, raw_required;
                            char raw_text[128];
                            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                                reference, sizeof reference - 1u, &raw_set, (int8_t)strand,
                                DUCKVEP_CODON_TABLE_STANDARD, raw_cds, sizeof raw_cds,
                                raw_rp, sizeof raw_rp, raw_ap, sizeof raw_ap, &raw_context));
                            raw_context.pre_cds_complete = raw_context.post_cds_complete = 1u;
                            raw_context.post_cds_bases = tail;
                            raw_context.post_cds_length = sizeof tail - 1u;
                            ASSERT_EQ(ctx.alt_cds_len, raw_context.alt_cds_len);
                            ASSERT_MEM_EQ(cds, raw_cds, ctx.alt_cds_len);
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                                duckvep_haplotype_partition(raw, 2u, raw_blocks, 2u, &raw_blocks_used));
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_build(
                                &raw_context, NULL, raw, 2u, raw_blocks, raw_blocks_used, 0u,
                                raw_facts, 2u, &raw_facts_used));
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_haplotype_render(
                                raw_facts, raw_facts_used, 1, raw_text, sizeof raw_text, &raw_required));
                            duckvep_edit_set_t raw_one = {raw + 1u, 1u};
                            duckvep_coding_context_t raw_single;
                            uint8_t raw_single_cds[40], raw_single_rp[16], raw_single_ap[16];
                            duckvep_sequence_delta_t raw_delta;
                            duckvep_hgvs_protein_fact_t raw_expected;
                            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                                reference, sizeof reference - 1u, &raw_one, (int8_t)strand,
                                DUCKVEP_CODON_TABLE_STANDARD,
                                raw_single_cds, sizeof raw_single_cds,
                                raw_single_rp, sizeof raw_single_rp,
                                raw_single_ap, sizeof raw_single_ap, &raw_single));
                            raw_single.pre_cds_complete = raw_single.post_cds_complete = 1u;
                            raw_single.post_cds_bases = tail;
                            raw_single.post_cds_length = sizeof tail - 1u;
                            ASSERT_EQ(single.alt_cds_len, raw_single.alt_cds_len);
                            ASSERT_MEM_EQ(single_cds, raw_single_cds, single.alt_cds_len);
                            ASSERT_MEM_EQ(single_ap, raw_single_ap, single.alt_peptide_len);
                            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                                duckvep_coding_context_delta_fill(&raw_single, 0u, &raw_delta));
                            ASSERT_EQ(DUCKVEP_HGVS_OK,
                                duckvep_hgvs_protein_fact_build(&raw_single, &raw_delta, &raw_expected));
                            ASSERT_EQ(raw_expected.shape == DUCKVEP_HGVS_PROTEIN_EQUAL ? 1u : 2u,
                                raw_facts_used);
                            char first_text[96], raw_first_text[96], one_text[96], raw_one_text[96];
                            size_t first_required, raw_first_required, one_required, raw_one_required;
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                &facts[0].fact, 0, first_text, sizeof first_text, &first_required));
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                &raw_facts[0].fact, 0, raw_first_text, sizeof raw_first_text,
                                &raw_first_required));
                            ASSERT_STR_EQ(first_text, raw_first_text);
                            ASSERT_EQ(first_required, raw_first_required);
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                &expected, 0, one_text, sizeof one_text, &one_required));
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                &raw_expected, 0, raw_one_text, sizeof raw_one_text, &raw_one_required));
                            if (raw_facts_used == 2u) {
                                char actual_text[96]; size_t actual_required;
                                const duckvep_hgvs_protein_fact_t *actual = &raw_facts[1].fact;
                                ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                    actual, 0, actual_text, sizeof actual_text, &actual_required));
                                ASSERT_STR_EQ(raw_one_text, actual_text);
                                ASSERT_EQ(raw_one_required, actual_required);
                                ASSERT_EQ(raw_expected.shape, actual->shape);
                                ASSERT_EQ(raw_expected.first_position1, actual->first_position1);
                                ASSERT_EQ(raw_expected.last_position1, actual->last_position1);
                                ASSERT_EQ(raw_expected.termination_distance, actual->termination_distance);
                                ASSERT_EQ(raw_expected.termination_known, actual->termination_known);
                            }
                            if (strcmp(minimal_text, raw_text) || minimal_required != raw_required) {
                                ASSERT(strcmp(one_text, raw_one_text) != 0);
                                anchor_presentation_differences++;
                            } else ASSERT_STR_EQ(one_text, raw_one_text);
                            anchored++;
                        }
                        if (expected.shape == DUCKVEP_HGVS_PROTEIN_EQUAL) {
                            ASSERT_EQ(1u, facts_used);
                            equal++;
                        } else {
                            ASSERT_EQ(2u, facts_used);
                            const duckvep_hgvs_protein_fact_t *fact = &facts[1].fact;
                            char rendered[96], isolated[96]; size_t required, isolated_required;
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                fact, 0, rendered, sizeof rendered, &required));
                            ASSERT_EQ(DUCKVEP_HGVS_OK, duckvep_hgvs_protein_render(
                                &expected, 0, isolated, sizeof isolated, &isolated_required));
                            ASSERT_STR_EQ(isolated, rendered);
                            ASSERT_EQ(isolated_required, required);
                            ASSERT_EQ(expected.shape, fact->shape);
                            ASSERT_EQ(expected.first_position1, fact->first_position1);
                            ASSERT_EQ(expected.last_position1, fact->last_position1);
                            ASSERT_EQ(expected.termination_distance, fact->termination_distance);
                            ASSERT_EQ(expected.termination_known, fact->termination_known);
                            /* An edit starting in a stop-containing peptide has VEP's
                             * stop-window presentation, including its synthetic X/Ter.
                             * It is checked against the complete singleton operation;
                             * raw translation alone does not define that HGVS spelling. */
                            if (!delta.frameshift) { stop_window++; cases++; continue; }
                            ASSERT(r < ref_length && a < alt_length);
                            ASSERT_EQ(r + 1u, fact->first_position1);
                            ASSERT_EQ(want_ref[r], fact->reference_first);
                            ASSERT_EQ(want_alt[a], fact->alternate_first);
                            if (want_alt[a] == '*') {
                                ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_SUBSTITUTION, fact->shape);
                                immediate_stop++;
                            } else {
                                ASSERT_EQ(DUCKVEP_HGVS_PROTEIN_FRAMESHIFT, fact->shape);
                                size_t stop = 0u;
                                while (stop < alt_length && want_alt[stop] != '*') stop++;
                                ASSERT_EQ(stop < alt_length && stop >= a, fact->termination_known);
                                ASSERT_EQ(stop < alt_length && stop >= a ? stop - a + 1u : 0u,
                                    fact->termination_distance);
                                frameshift++;
                            }
                        }
                        cases++;
                    }
                }
            }
        }
    }
    ASSERT_EQ(262144u, cases);
    ASSERT_EQ(524288u, anchored);
    ASSERT(equal && immediate_stop && frameshift && stop_window);
    ASSERT_EQ(cases, equal + immediate_stop + frameshift + stop_window);
    fprintf(stderr, "[compound HGVSp codon landscape] cases=%zu equal_stop=%zu "
        "immediate_stop=%zu frameshift=%zu stop_window=%zu anchored=%zu "
        "anchor_presentation_differences=%zu\n",
        cases, equal, immediate_stop, frameshift, stop_window, anchored,
        anchor_presentation_differences);
    ASSERT_EQ(3072u, anchor_presentation_differences);
    PASS();
}

TEST haplotype_substitution_blocks_reuse_local_coding_predicates(void) {
    static const uint8_t bases[] = "ACGT";
    size_t cases = 0u, with_stop_before = 0u;
    size_t synonymous = 0u, missense = 0u, gained = 0u, lost = 0u, retained = 0u;
    /* Every distinct reference/alternate codon, with a preceding -3/0/+3 edit
     * and both allele orientations. The independent comparator is one actual
     * substitution of the same codon; production receives the physical SNVs,
     * not that test-only equivalent MNV. Stops in the prefix do not erase local
     * predicates: they are a separate path fact for the eventual SO consumer. */
    for (unsigned r = 0u; r < 64u; r++) {
        for (unsigned a = 0u; a < 64u; a++) {
            if (r == a) continue;
            uint8_t reference[] = "ATGGGGCCCAAACCCTAA";
            uint8_t alternate[3];
            for (unsigned i = 0u; i < 3u; i++) {
                reference[9u + i] = bases[(r >> (4u - 2u * i)) & 3u];
                alternate[i] = bases[(a >> (4u - 2u * i)) & 3u];
            }
            for (int shift = -3; shift <= 3; shift += 3) {
                for (int strand = -1; strand <= 1; strand += 2) {
                    uint8_t ref[4][3] = {{0}}, alt[4][3] = {{0}}, mref[3], malt[3];
                    uint8_t cds[32], rp[12], ap[12], mcds[32], mrp[12], map[12];
                    duckvep_haplotype_edit_t ascending[4], descending[4];
                    size_t n = 1u;
                    ascending[0] = (duckvep_haplotype_edit_t){4u,
                        shift > 0 ? 0u : 3u, ref[0], shift < 0 ? 0u : 3u, alt[0], 1};
                    for (size_t i = 0u; i < 3u; i++) {
                        if (reference[9u + i] == alternate[i]) continue;
                        ascending[n] = (duckvep_haplotype_edit_t){(uint32_t)(10u + i),
                            1u, ref[n], 1u, alt[n], 1};
                        ref[n][0] = reference[9u + i];
                        alt[n][0] = alternate[i];
                        n++;
                    }
                    memcpy(ref[0], "GGG", 3u); memcpy(alt[0], "TGA", 3u);
                    memcpy(mref, reference + 9u, 3u); memcpy(malt, alternate, 3u);
                    if (strand < 0) {
                        for (size_t e = 0u; e < n; e++) {
                            uint8_t saved_ref[3], saved_alt[3];
                            memcpy(saved_ref, ref[e], 3u); memcpy(saved_alt, alt[e], 3u);
                            for (size_t i = 0u; i < ascending[e].ref_len; i++)
                                ref[e][i] = (uint8_t)kprop_complement_base(
                                    (char)saved_ref[ascending[e].ref_len - 1u - i]);
                            for (size_t i = 0u; i < ascending[e].alt_len; i++)
                                alt[e][i] = (uint8_t)kprop_complement_base(
                                    (char)saved_alt[ascending[e].alt_len - 1u - i]);
                        }
                        for (size_t i = 0u; i < 3u; i++) {
                            mref[i] = (uint8_t)kprop_complement_base((char)reference[11u - i]);
                            malt[i] = (uint8_t)kprop_complement_base((char)alternate[2u - i]);
                        }
                    }
                    for (size_t i = 0u; i < n; i++) descending[i] = ascending[n - 1u - i];
                    duckvep_edit_set_t set = {descending, n};
                    duckvep_coding_context_t ctx, independent;
                    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                        reference, sizeof reference - 1u, &set, (int8_t)strand,
                        DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds,
                        rp, sizeof rp, ap, sizeof ap, &ctx));
                    duckvep_haplotype_edit_t mnv = {10u, 3u, mref, 3u, malt, 1};
                    set = (duckvep_edit_set_t){&mnv, 1u};
                    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                        reference, sizeof reference - 1u, &set, (int8_t)strand,
                        DUCKVEP_CODON_TABLE_STANDARD, mcds, sizeof mcds,
                        mrp, sizeof mrp, map, sizeof map, &independent));
                    duckvep_haplotype_block_t blocks[4];
                    size_t count;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                        duckvep_haplotype_partition(ascending, n, blocks, 4u, &count));
                    ASSERT_EQ(2u, count);
                    ASSERT_EQ(0u, blocks[1].flags);
                    ASSERT_EQ(n - 1u, blocks[1].edit_count);
                    ASSERT_EQ((int64_t)blocks[1].cds_start - 1 + shift,
                              (int64_t)blocks[1].alt_start0);
                    duckvep_sequence_delta_t actual, expected;
                    duckvep_coding_context_t saved = ctx;
                    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                        duckvep_coding_context_block_delta_fill(&ctx, ascending, n, blocks + 1u, 0u, &actual));
                    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                        duckvep_coding_context_delta_fill(&independent, 0u, &expected));
                    ASSERT_MEM_EQ(&expected, &actual, sizeof actual);
                    ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
                    ASSERT_EQ(4, actual.protein_pos);
                    synonymous += actual.synonymous; missense += actual.missense;
                    gained += actual.stop_gained; lost += actual.stop_lost;
                    retained += actual.stop_retained;
                    size_t first_stop = 0u;
                    while (first_stop < ctx.alt_peptide_len && ap[first_stop] != '*') first_stop++;
                    with_stop_before += first_stop < blocks[1].alt_start0 / 3u;
                    for (unsigned bad = 0u; bad < 4u; bad++) {
                        duckvep_haplotype_block_t invalid = blocks[1];
                        duckvep_sequence_delta_t empty = {0};
                        if (bad == 0u) invalid.flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
                        if (bad == 1u) invalid.length_diff = 3;
                        if (bad == 2u) invalid.alt_start0 = SIZE_MAX;
                        if (bad == 3u) invalid.edit_count = n + 1u;
                        memset(&actual, 0xff, sizeof actual);
                        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                            duckvep_coding_context_block_delta_fill(&ctx, ascending, n, &invalid, 0u, &actual));
                        ASSERT_MEM_EQ(&empty, &actual, sizeof actual);
                    }
                    cases++;
                }
            }
        }
    }
    ASSERT_EQ(24192u, cases);
    ASSERT_EQ(16128u, with_stop_before);
    ASSERT_EQ(1044u, synonymous); ASSERT_EQ(20916u, missense);
    ASSERT_EQ(1098u, gained); ASSERT_EQ(1098u, lost); ASSERT_EQ(36u, retained);
    duckvep_sequence_delta_t delta, empty = {0};
    duckvep_coding_context_t ctx = {0};
    duckvep_haplotype_block_t block = {0};
    memset(&delta, 0xff, sizeof delta);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_INVALID_ARG,
        duckvep_coding_context_block_delta_fill(NULL, NULL, 0u, &block, 0u, &delta));
    ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
    memset(&delta, 0xff, sizeof delta);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_INVALID_ARG,
        duckvep_coding_context_block_delta_fill(&ctx, NULL, 0u, NULL, 0u, &delta));
    ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_INVALID_ARG,
        duckvep_coding_context_block_delta_fill(&ctx, NULL, 0u, &block, 0u, NULL));
    /* A changed prefix does not turn a later identity edit into a substitution. */
    const uint8_t reference[] = "ATGGGGCCCAAACCCTAA";
    const uint8_t inserted[] = "TGA";
    uint8_t cds[24], rp[8], ap[8];
    duckvep_haplotype_edit_t edits[2] = {
        {10u, 3u, reference + 9u, 3u, reference + 9u, 1},
        {4u, 0u, NULL, 3u, inserted, 1}
    };
    duckvep_edit_set_t set = {edits, 2u};
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    duckvep_haplotype_edit_t ascending[2] = {edits[1], edits[0]};
    duckvep_haplotype_block_t blocks[2];
    size_t count;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &count));
    ASSERT_EQ(2u, count);
    ASSERT(ctx.cds_changed);
    memset(&delta, 0xff, sizeof delta);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
        duckvep_coding_context_block_delta_fill(&ctx, ascending, 2u, blocks + 1u, 0u, &delta));
    ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
    PASS();
}

TEST haplotype_restoring_indels_can_recreate_reference_sequence(void) {
    static const uint8_t bases[] = "ACGT";
    static const uint32_t widths[] = {1u, 2u, 4u, 5u};
    size_t cases = 0u;
    for (size_t base = 0u; base < 4u; base++) {
        uint8_t reference[28] = "ATG";
        memset(reference + 3u, bases[base], 18u);
        memcpy(reference + 21u, "CCCTAA", 6u);
        for (size_t width = 0u; width < 4u; width++) {
            uint32_t n = widths[width];
            for (uint32_t deleted = 4u; deleted <= 12u; deleted++) {
                for (uint32_t restored = deleted + n; restored <= 21u; restored++) {
                    for (int strand = -1; strand <= 1; strand += 2) {
                        uint8_t allele[5];
                        memset(allele, strand > 0 ? bases[base]
                            : (uint8_t)kprop_complement_base((char)bases[base]), n);
                        duckvep_haplotype_edit_t ascending[2] = {
                            {deleted, n, allele, 0u, NULL, 1},
                            {restored, 0u, NULL, n, allele, 1}};
                        duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
                        duckvep_edit_set_t set = {descending, 2u};
                        duckvep_coding_context_t ctx;
                        uint8_t cds[40], rp[16], ap[16];
                        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                            reference, sizeof reference - 1u, &set, (int8_t)strand,
                            DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds,
                            rp, sizeof rp, ap, sizeof ap, &ctx));
                        ASSERT_EQ(sizeof reference - 1u, ctx.alt_cds_len);
                        ASSERT_MEM_EQ(reference, cds, sizeof reference - 1u);
                        ASSERT_MEM_EQ(rp, ap, ctx.ref_peptide_len);
                        ASSERT(!ctx.cds_changed);
                        ASSERT_EQ(2u, ctx.applied_edits);
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_FLAG_INDEL |
                            DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT, ctx.flags);
                        duckvep_haplotype_block_t blocks[2];
                        size_t count;
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                            duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &count));
                        ASSERT_EQ(1u, count);
                        ASSERT_EQ(2u, blocks[0].edit_count);
                        ASSERT_EQ(ctx.flags, blocks[0].flags);
                        ASSERT_EQ(0, blocks[0].length_diff);
                        duckvep_coding_context_t saved = ctx;
                        duckvep_sequence_delta_t delta;
                        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                            duckvep_coding_context_block_delta_fill(
                                &ctx, ascending, 2u, blocks, 0u, &delta));
                        ASSERT(delta.valid && delta.synonymous && !delta.missense &&
                            !delta.frameshift && !delta.inframe_deletion &&
                            !delta.inframe_insertion && !delta.protein_altering &&
                            !delta.stop_gained && !delta.stop_lost && !delta.stop_retained);
                        ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
                        /* Sequence equality does not license the forbidden
                         * whole-context substitution approximation. */
                        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                            duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
                        ASSERT(!delta.valid);
                        cases++;
                    }
                }
            }
        }
    }
    ASSERT_EQ(3168u, cases);
    PASS();
}

TEST haplotype_indel_blocks_reuse_local_coding_predicates(void) {
    static const uint8_t bases[] = "ACGT";
    size_t cases = 0u;
    for (unsigned codon = 0u; codon < 64u; codon++) {
        uint8_t reference[] = "ATGGGGCCCAAACCCTTTGGGCCCTAA";
        for (unsigned i = 0u; i < 3u; i++)
            reference[9u + i] = bases[(codon >> (4u - 2u * i)) & 3u];
        for (uint32_t position = 10u; position <= 14u; position++) {
            for (uint32_t r = 0u; r <= 3u; r++) {
                for (uint32_t a = 0u; a <= 4u; a++) {
                    if (r == a) continue;
                    for (int shift = -3; shift <= 3; shift += 3) {
                        for (int strand = -1; strand <= 1; strand += 2) {
                            uint8_t ref[2][3] = {{0}}, alt[2][4] = {{0}};
                            duckvep_haplotype_edit_t ascending[2] = {
                                {4u, shift > 0 ? 0u : 3u, ref[0], shift < 0 ? 0u : 3u, alt[0], 1},
                                {position, r, ref[1], a, alt[1], 1}};
                            const uint8_t *replacements[2] = {(const uint8_t *)"TGA", (const uint8_t *)"TACG"};
                            for (size_t e = 0u; e < 2u; e++) {
                                for (size_t i = 0u; i < ascending[e].ref_len; i++) {
                                    size_t at = ascending[e].cds_start - 1u +
                                        (strand > 0 ? i : ascending[e].ref_len - 1u - i);
                                    ref[e][i] = strand > 0 ? reference[at]
                                        : (uint8_t)kprop_complement_base((char)reference[at]);
                                }
                                for (size_t i = 0u; i < ascending[e].alt_len; i++) {
                                    uint8_t b = replacements[e][strand > 0 ? i : ascending[e].alt_len - 1u - i];
                                    alt[e][i] = strand > 0 ? b : (uint8_t)kprop_complement_base((char)b);
                                }
                            }
                            duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
                            uint8_t cds[40], rp[16], ap[16], only_cds[40], only_rp[16], only_ap[16];
                            duckvep_coding_context_t ctx, independent;
                            duckvep_edit_set_t set = {descending, 2u};
                            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                                reference, sizeof reference - 1u, &set, (int8_t)strand,
                                DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
                            set = (duckvep_edit_set_t){ascending + 1u, 1u};
                            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                                reference, sizeof reference - 1u, &set, (int8_t)strand,
                                DUCKVEP_CODON_TABLE_STANDARD, only_cds, sizeof only_cds,
                                only_rp, sizeof only_rp, only_ap, sizeof only_ap, &independent));
                            duckvep_haplotype_block_t blocks[2];
                            size_t count;
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                                duckvep_haplotype_partition(ascending, 2u, blocks, 2u, &count));
                            ASSERT_EQ(2u, count);
                            duckvep_coding_context_t saved = ctx;
                            duckvep_sequence_delta_t actual, expected;
                            duckvep_context_delta_status_t expected_status =
                                duckvep_coding_context_delta_fill(&independent, 0u, &expected);
                            ASSERT_EQ(expected_status, duckvep_coding_context_block_delta_fill(
                                &ctx, ascending, 2u, blocks + 1u, 0u, &actual));
                            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK, expected_status);
                            ASSERT_MEM_EQ(&expected, &actual, sizeof actual);
                            ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
                            cases++;
                        }
                    }
                }
            }
        }
    }
    ASSERT_EQ(30720u, cases);
    /* Unchanged seed-173 DHT000002. A DNA-restoring pair has already reached
     * a stop while displaced. It must not acquire in-frame or missense facts. */
    static const uint8_t reference[] =
        "ATGGGCTTGCCTAGCTTAAAACACTTGGTAAACTTTCTCGTGGTTACAACCTCGTTAGGGCTGCAT"
        "CTTCTCTACATGGAGAAGTGCGAAGTCCGGGTAAGGACCGGGTGTTATAACCCAGCAGACAAGCT"
        "GAAGGGGCGCGTGTACTTAGCTGCCCGCTCGAGGCATTGGTCCCCGTAA";
    duckvep_haplotype_edit_t edits[2] = {
        {10u, 0u, NULL, 1u, (const uint8_t *)"T", 1},
        {40u, 1u, reference + 39u, 0u, NULL, 1}};
    duckvep_haplotype_edit_t descending[2] = {edits[1], edits[0]};
    duckvep_edit_set_t set = {descending, 2u};
    duckvep_coding_context_t ctx;
    uint8_t cds[192], rp[64], ap[64];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        reference, sizeof reference - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_EQ(5u, ctx.alt_first_stop_position1);
    duckvep_haplotype_block_t block;
    size_t count;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(edits, 2u, &block, 1u, &count));
    ASSERT_EQ(1u, count);
    duckvep_sequence_delta_t delta;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
        duckvep_coding_context_block_delta_fill(&ctx, edits, 2u, &block, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
    ASSERT(!delta.inframe_insertion && !delta.inframe_deletion && !delta.missense && !delta.synonymous);
    /* The whole-context substitution approximation stays forbidden. */
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED, duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);
    static const uint8_t later_stop[] = "ATGAAATAACCC";
    duckvep_haplotype_edit_t before_stop[2] = {
        {5u, 0u, NULL, 1u, (const uint8_t *)"C", 1},
        {6u, 1u, later_stop + 5u, 0u, NULL, 1}};
    descending[0] = before_stop[1]; descending[1] = before_stop[0];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        later_stop, sizeof later_stop - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_EQ(3u, ctx.alt_first_stop_position1);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(before_stop, 2u, &block, 1u, &count));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
        duckvep_coding_context_block_delta_fill(&ctx, before_stop, 2u, &block, 0u, &delta));
    ASSERT(delta.valid && delta.missense && !delta.frameshift && !delta.stop_gained);
    /* Adjacent deletions leave no alternate base inside the displaced span. */
    before_stop[0] = (duckvep_haplotype_edit_t){4u, 1u, later_stop + 3u, 0u, NULL, 1};
    before_stop[1] = (duckvep_haplotype_edit_t){5u, 2u, later_stop + 4u, 0u, NULL, 1};
    descending[0] = before_stop[1]; descending[1] = before_stop[0];
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
        later_stop, sizeof later_stop - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
        cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
    ASSERT_EQ(2u, ctx.alt_first_stop_position1);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(before_stop, 2u, &block, 1u, &count));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
        duckvep_coding_context_block_delta_fill(&ctx, before_stop, 2u, &block, 0u, &delta));
    ASSERT(delta.valid && delta.inframe_deletion && !delta.frameshift && !delta.stop_gained);
    /* Actual compound endpoint edits retain physical frame geometry while
     * borrowing their rebuilt sequence spans for start/stop predicates. */
    for (uint32_t start1 = 2u; start1 <= 11u; start1 += 9u) {
        before_stop[0] = (duckvep_haplotype_edit_t){start1, 1u, later_stop + start1 - 1u, 0u, NULL, 1};
        before_stop[1] = (duckvep_haplotype_edit_t){start1 + 1u, 0u, NULL, 1u, (const uint8_t *)"G", 1};
        descending[0] = before_stop[1]; descending[1] = before_stop[0];
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
            later_stop, sizeof later_stop - 1u, &set, 1, DUCKVEP_CODON_TABLE_STANDARD,
            cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, duckvep_haplotype_partition(before_stop, 2u, &block, 1u, &count));
        ASSERT_EQ(1u, count);
        duckvep_coding_context_t saved = ctx;
        duckvep_sequence_delta_t empty = {0};
        memset(&delta, 0xff, sizeof delta);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
            duckvep_coding_context_block_delta_fill(&ctx, before_stop, 2u, &block, 0u, &delta));
        ASSERT(delta.valid && !delta.frameshift && !delta.stop_gained && !delta.stop_lost);
        if (start1 == 2u) {
            ASSERT(delta.start_lost && !delta.start_retained && !delta.missense);
            ctx.pre_cds_complete = 0u;
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK,
                duckvep_coding_context_block_delta_fill(&ctx, before_stop, 2u, &block, 0u, &delta));
            ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
            ctx.pre_cds_complete = saved.pre_cds_complete;
        } else {
            ASSERT(delta.missense && !delta.start_lost && !delta.start_retained);
        }
        ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
    }
    PASS();
}

TEST haplotype_endpoint_deletions_match_independent_span(void) {
    static const uint8_t bases[] = "ACGT";
    static const uint8_t flank[] = "TAATGA";
    const uint32_t positions[] = {1u, 2u, 3u, 23u, 24u, 25u};
    size_t cases = 0u, supported = 0u, missing_tail = 0u;
    for (unsigned codon = 0u; codon < 64u; codon++) {
        for (size_t p = 0u; p < sizeof positions / sizeof *positions; p++) {
            uint8_t reference[] = "ATGGGGCCCAAACCCTTTGGGCCCTAA";
            uint32_t position = positions[p];
            for (unsigned i = 0u; i < 3u; i++)
                reference[(p < 3u ? 0u : 24u) + i] = bases[(codon >> (4u - 2u * i)) & 3u];
            for (int shift = -3; shift <= 3; shift += 3) {
                for (int strand = -1; strand <= 1; strand += 2) {
                    uint8_t ref[3][3] = {{0}}, alt[3] = {0}, combined_ref[3];
                    duckvep_haplotype_edit_t edited[3] = {
                        {position, 1u, ref[0], 0u, NULL, 1},
                        {position + 1u, 2u, ref[1], 0u, NULL, 1},
                        {p < 3u ? 19u : 4u, shift > 0 ? 0u : 3u, ref[2],
                            shift < 0 ? 0u : 3u, alt, 1}};
                    for (size_t e = 0u; e < 3u; e++) {
                        for (size_t i = 0u; i < edited[e].ref_len; i++) {
                            size_t at = edited[e].cds_start - 1u +
                                (strand > 0 ? i : edited[e].ref_len - 1u - i);
                            ref[e][i] = strand > 0 ? reference[at]
                                : (uint8_t)kprop_complement_base((char)reference[at]);
                        }
                    }
                    for (size_t i = 0u; i < 3u; i++) {
                        uint8_t b = ((const uint8_t *)"TGC")[strand > 0 ? i : 2u - i];
                        alt[i] = strand > 0 ? b : (uint8_t)kprop_complement_base((char)b);
                        b = reference[position - 1u + (strand > 0 ? i : 2u - i)];
                        combined_ref[i] = strand > 0 ? b : (uint8_t)kprop_complement_base((char)b);
                    }
                    duckvep_haplotype_edit_t ascending[3];
                    for (size_t e = 0u; e < 3u; e++)
                        ascending[e] = edited[p < 3u ? e : (e + 2u) % 3u];
                    duckvep_haplotype_edit_t descending[3] = {ascending[2], ascending[1], ascending[0]};
                    duckvep_haplotype_edit_t combined = {position, 3u, combined_ref, 0u, NULL, 1};
                    duckvep_edit_set_t set = {descending, 3u};
                    duckvep_coding_context_t ctx, independent;
                    uint8_t cds[40], rp[16], ap[16], only_cds[40], only_rp[16], only_ap[16];
                    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                        reference, sizeof reference - 1u, &set, (int8_t)strand,
                        DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
                    set = (duckvep_edit_set_t){&combined, 1u};
                    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                        reference, sizeof reference - 1u, &set, (int8_t)strand,
                        DUCKVEP_CODON_TABLE_STANDARD, only_cds, sizeof only_cds,
                        only_rp, sizeof only_rp, only_ap, sizeof only_ap, &independent));
                    duckvep_haplotype_block_t blocks[3];
                    size_t count;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                        duckvep_haplotype_partition(ascending, 3u, blocks, 3u, &count));
                    ASSERT_EQ(2u, count);
                    const duckvep_haplotype_block_t *block = blocks + (p < 3u ? 0u : 1u);
                    ASSERT_EQ(2u, block->edit_count);
                    for (size_t tail = 0u; tail <= sizeof flank - 1u; tail++) {
                        for (uint8_t complete = 0u; complete <= 1u; complete++) {
                            for (size_t utr = 0u; utr <= 3u; utr += 3u) {
                                ctx.pre_cds_bases = independent.pre_cds_bases = flank;
                                ctx.pre_cds_length = independent.pre_cds_length = utr;
                                ctx.post_cds_bases = independent.post_cds_bases = flank;
                                ctx.post_cds_length = independent.post_cds_length = tail;
                                ctx.post_cds_complete = independent.post_cds_complete = complete;
                                duckvep_coding_context_t saved = ctx;
                                duckvep_sequence_delta_t actual, expected, empty = {0};
                                duckvep_context_delta_status_t status =
                                    duckvep_coding_context_delta_fill(&independent, 0u, &expected);
                                ASSERT_EQ(status, duckvep_coding_context_block_delta_fill(
                                    &ctx, ascending, 3u, block, 0u, &actual));
                                if (status == DUCKVEP_CONTEXT_DELTA_OK) {
                                    ASSERT_MEM_EQ(&expected, &actual, sizeof actual);
                                    supported++;
                                } else {
                                    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL, status);
                                    ASSERT_MEM_EQ(&empty, &actual, sizeof actual);
                                    missing_tail++;
                                }
                                ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
                                cases++;
                            }
                        }
                    }
                }
            }
        }
    }
    ASSERT_EQ(64512u, cases);
    ASSERT(supported != 0u && missing_tail != 0u);
    PASS();
}

TEST haplotype_endpoint_frames_follow_physical_edits(void) {
    static const uint8_t reference[] = "ATGAAACCCTAA";
    const uint32_t starts[] = {1u, 11u, 11u};
    const uint32_t restores[] = {4u, 12u, 12u};
    const uint8_t inserted[] = {'A', 'G', 'C'};
    for (size_t scene = 0u; scene < 3u; scene++) {
        for (int strand = -1; strand <= 1; strand += 2) {
            uint8_t ref = reference[starts[scene] - 1u], alt = inserted[scene];
            if (strand < 0) {
                ref = (uint8_t)kprop_complement_base((char)ref);
                alt = (uint8_t)kprop_complement_base((char)alt);
            }
            duckvep_haplotype_edit_t ascending[2] = {
                {starts[scene], 1u, &ref, 0u, NULL, 1},
                {restores[scene], 0u, NULL, 1u, &alt, 1}};
            duckvep_haplotype_edit_t descending[2] = {ascending[1], ascending[0]};
            duckvep_edit_set_t set = {descending, 2u};
            duckvep_coding_context_t ctx;
            uint8_t cds[16], rp[8], ap[8];
            ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                reference, sizeof reference - 1u, &set, (int8_t)strand,
                DUCKVEP_CODON_TABLE_STANDARD, cds, sizeof cds, rp, sizeof rp, ap, sizeof ap, &ctx));
            duckvep_haplotype_block_t block;
            size_t count;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                duckvep_haplotype_partition(ascending, 2u, &block, 1u, &count));
            ASSERT_EQ(1u, count);
            ASSERT_EQ(0, block.length_diff);
            ASSERT(block.flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL);
            duckvep_coding_context_t saved = ctx;
            duckvep_sequence_delta_t delta, empty = {0};
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                duckvep_coding_context_block_delta_fill(&ctx, ascending, 2u, &block, 0u, &delta));
            ASSERT(delta.valid);
            ASSERT_EQ(scene == 0u, delta.start_lost);
            ASSERT_EQ(scene == 0u, delta.frameshift);
            ASSERT_EQ(scene == 0u, delta.stop_gained);
            ASSERT_EQ(scene == 1u, delta.stop_retained);
            ASSERT_EQ(scene == 2u, delta.stop_lost);
            ASSERT(!delta.start_retained && !delta.inframe_insertion && !delta.inframe_deletion);
            ASSERT(!delta.protein_altering && !delta.synonymous && !delta.missense);
            ASSERT_MEM_EQ(&saved, &ctx, sizeof ctx);
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
            ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
            ctx.alt_cds = NULL;
            ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                duckvep_coding_context_block_delta_fill(&ctx, ascending, 2u, &block, 0u, &delta));
            ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
        }
    }
    PASS();
}

TEST haplotype_frame_spans_match_rebuilt_base_markers(void) {
    const uint8_t bases[] = "AAAA";
    size_t cases = 0u, queries = 0u, displaced = 0u;
    for (int prefix = -3; prefix <= 3; prefix += 3) {
        for (uint32_t second = 10u; second <= 18u; second++) {
            for (uint32_t r1 = 0u; r1 <= 3u; r1++) {
                for (uint32_t a1 = 0u; a1 <= 4u; a1++) {
                    if (r1 == 0u && a1 == 0u) continue;
                    for (uint32_t r2 = 0u; r2 <= 3u; r2++) {
                        for (uint32_t a2 = 0u; a2 <= 4u; a2++) {
                            if (r2 == 0u && a2 == 0u) continue;
                            duckvep_haplotype_edit_t edits[3];
                            size_t n = 0u;
                            if (prefix) edits[n++] = (duckvep_haplotype_edit_t){1u,
                                prefix < 0 ? 3u : 0u, bases, prefix > 0 ? 3u : 0u, bases, 1};
                            edits[n++] = (duckvep_haplotype_edit_t){7u, r1, bases, a1, bases, 1};
                            edits[n++] = (duckvep_haplotype_edit_t){second, r2, bases, a2, bases, -1};
                            /* Independent per-base sidecar reconstruction: an unchanged
                             * base's two coordinates give its frame. Replacement bases
                             * retain either the entering or leaving displaced frame. */
                            uint8_t marker[64] = {0};
                            size_t read = 0u, out = 0u;
                            for (size_t e = 0u; e < n; e++) {
                                size_t start0 = edits[e].cds_start - 1u;
                                while (read < start0) {
                                    marker[out] = ((int64_t)out - (int64_t)read) % 3 != 0;
                                    read++; out++;
                                }
                                int before = ((int64_t)out - (int64_t)read) % 3 != 0;
                                int after = ((int64_t)out + edits[e].alt_len -
                                    (int64_t)read - edits[e].ref_len) % 3 != 0;
                                for (size_t a = 0u; a < edits[e].alt_len; a++)
                                    marker[out++] = (uint8_t)(before || after);
                                read += edits[e].ref_len;
                            }
                            while (read < 42u) {
                                marker[out] = ((int64_t)out - (int64_t)read) % 3 != 0;
                                read++; out++;
                            }
                            duckvep_haplotype_block_t blocks[3];
                            size_t count;
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                                duckvep_haplotype_partition(edits, n, blocks, 3u, &count));
                            for (size_t b = 0u; b < count; b++) {
                                size_t end = b + 1u < count ? blocks[b + 1u].alt_start0 : out;
                                for (size_t p = blocks[b].alt_start0; p <= end; p++) {
                                    for (size_t length = 0u; length <= 3u && length <= end - p; length++) {
                                        int expected = 0, observed = -1;
                                        for (size_t i = p; i < p + length; i++) expected |= marker[i];
                                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                                            duckvep_haplotype_block_frame_intersects(
                                                edits, n, blocks + b, p, length, &observed));
                                        ASSERT_EQ(expected, observed);
                                        queries++; displaced += (size_t)observed;
                                    }
                                }
                                int observed = 1;
                                ASSERT_EQ(DUCKVEP_HAPLOTYPE_OUT_OF_RANGE,
                                    duckvep_haplotype_block_frame_intersects(
                                        edits, n, blocks + b, SIZE_MAX, 1u, &observed));
                                ASSERT_EQ(0, observed);
                                duckvep_haplotype_block_t invalid = blocks[b];
                                invalid.flags ^= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
                                observed = 1;
                                ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG,
                                    duckvep_haplotype_block_frame_intersects(
                                        edits, n, &invalid, 0u, out, &observed));
                                ASSERT_EQ(0, observed);
                            }
                            cases++;
                        }
                    }
                }
            }
        }
    }
    ASSERT_EQ(9747u, cases); ASSERT(queries > 1000000u); ASSERT(displaced > 0u);
    /* Adjacent 1-base and 2-base deletions restore the frame without leaving
     * any alternate bases inside the displacement. A spanning codon query
     * must not mistake the empty interval for a translated excursion. */
    duckvep_haplotype_edit_t deletions[2] = {
        {4u, 1u, bases, 0u, NULL, 1}, {5u, 2u, bases, 0u, NULL, 1}};
    duckvep_haplotype_block_t block;
    size_t count;
    int observed;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_partition(deletions, 2u, &block, 1u, &count));
    ASSERT_EQ(1u, count);
    ASSERT(block.flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
        duckvep_haplotype_block_frame_intersects(deletions, 2u, &block, 2u, 3u, &observed));
    ASSERT_EQ(0, observed);
    PASS();
}

TEST haplotype_compound_indels_are_not_substitution_facts(void) {
    /* Unchanged Haplosaurus generator seed 173, DHT000002: genomic 211 G>GA
     * and 180 AC>A on the reverse strand. The frame returns to zero at CDS 40,
     * but the combined protein is already MGLS*. Pinned bcftools 1.23 csq
     * reports stop_gained&frameshift, not a length-preserving substitution.
     * Extend this witness across both allele orientations and net length signs. */
    static const uint8_t cds[] =
        "ATGGGCTTGCCTAGCTTAAAACACTTGGTAAACTTTCTCGTGGTTACAACCTCGTTAGGGCTGCAT"
        "CTTCTCTACATGGAGAAGTGCGAAGTCCGGGTAAGGACCGGGTGTTATAACCCAGCAGACAAGCT"
        "GAAGGGGCGCGTGTACTTAGCTGCCCGCTCGAGGCATTGGTCCCCGTAA";
    static const uint8_t inserted[] = "TGG";
    for (int strand = -1; strand <= 1; strand += 2) {
        for (uint32_t inserted_len = 1u; inserted_len <= 3u; inserted_len++) {
            for (uint32_t deleted_len = 1u; deleted_len <= 3u; deleted_len++) {
                uint8_t ref[3], alt[3], alt_cds[192], ref_peptide[64], alt_peptide[64];
                duckvep_coding_context_t ctx;
                duckvep_sequence_delta_t delta;
                for (uint32_t i = 0u; i < deleted_len; i++) {
                    ref[i] = strand > 0 ? cds[39u + i]
                        : (uint8_t)kprop_complement_base((char)cds[39u + deleted_len - 1u - i]);
                }
                for (uint32_t i = 0u; i < inserted_len; i++) {
                    alt[i] = strand > 0 ? inserted[i]
                        : (uint8_t)kprop_complement_base((char)inserted[inserted_len - 1u - i]);
                }
                duckvep_haplotype_edit_t edits[2] = {
                    {40u, deleted_len, ref, 0u, NULL, 1},
                    {10u, 0u, NULL, inserted_len, alt, 1}
                };
                duckvep_edit_set_t set = {edits, 2u};
                ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK, duckvep_coding_context_build(
                    cds, sizeof cds - 1u, &set, (int8_t)strand, DUCKVEP_CODON_TABLE_STANDARD,
                    alt_cds, sizeof alt_cds, ref_peptide, sizeof ref_peptide,
                    alt_peptide, sizeof alt_peptide, &ctx));
                ASSERT_EQ(2u, ctx.applied_edits);
                ASSERT(!ctx.has_single_edit);
                ASSERT_EQ((int64_t)inserted_len - deleted_len, ctx.length_diff);
                ASSERT(ctx.flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL);
                if (inserted_len == 1u && deleted_len == 1u) {
                    ASSERT(ctx.flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT);
                    ASSERT_EQ(sizeof cds - 1u, ctx.alt_cds_len);
                    ASSERT_EQ(60u, ctx.alt_peptide_len);
                    ASSERT_MEM_EQ("MGLS*", ctx.alt_peptide, 5u);
                    ASSERT(ctx.alt_peptide[5u] != 0u);
                }
                memset(&delta, 0xff, sizeof delta);
                ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                    duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
                duckvep_sequence_delta_t empty = {0};
                ASSERT_MEM_EQ(&empty, &delta, sizeof delta);
            }
        }
    }
    PASS();
}

TEST haplotype_snv_set_matches_equivalent_mnv_coding_facts(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "phased SNV set == equivalent MNV coding facts";
    cfg.prop1 = prop_haplotype_snv_set_matches_equivalent_mnv;
    cfg.type_info[0] = &kprop_cds_edit_set_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_haplotype_mnv_equivalence_cov, 0,
           sizeof g_haplotype_mnv_equivalence_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_haplotype_mnv_equivalence_cov.fwd > 0u);
    ASSERT(g_haplotype_mnv_equivalence_cov.rev > 0u);
    ASSERT(g_haplotype_mnv_equivalence_cov.start > 0u);
    ASSERT(g_haplotype_mnv_equivalence_cov.body > 0u);
    ASSERT(g_haplotype_mnv_equivalence_cov.stop > 0u);
    ASSERT(g_haplotype_mnv_equivalence_cov.one_codon > 0u);
    ASSERT(g_haplotype_mnv_equivalence_cov.several_codons > 0u);
    fprintf(stderr,
            "[haplotype-MNV equivalence coverage] fwd=%u rev=%u start=%u "
            "body=%u stop=%u one_codon=%u several_codons=%u\n",
            g_haplotype_mnv_equivalence_cov.fwd,
            g_haplotype_mnv_equivalence_cov.rev,
            g_haplotype_mnv_equivalence_cov.start,
            g_haplotype_mnv_equivalence_cov.body,
            g_haplotype_mnv_equivalence_cov.stop,
            g_haplotype_mnv_equivalence_cov.one_codon,
            g_haplotype_mnv_equivalence_cov.several_codons);
    PASS();
}
