#include "duckvep_property.h"

/* ===================================================================== *
 * Codon translation + coding-change classification, incl. the mitochondrial
 * table. The canonical genetic code is the independent oracle.
 * ===================================================================== */

#define MITO DUCKVEP_CODON_TABLE_VERT_MITO

TEST codon_translate_known(void) {
    /* standard (transl_table=1) */
    ASSERT_EQ('M', duckvep_translate_codon("ATG", STD));
    ASSERT_EQ('W', duckvep_translate_codon("TGG", STD));
    ASSERT_EQ('*', duckvep_translate_codon("TAA", STD));
    ASSERT_EQ('*', duckvep_translate_codon("TAG", STD));
    ASSERT_EQ('*', duckvep_translate_codon("TGA", STD));
    ASSERT_EQ('I', duckvep_translate_codon("ATA", STD));
    ASSERT_EQ('R', duckvep_translate_codon("AGA", STD));
    ASSERT_EQ('R', duckvep_translate_codon("AGG", STD));
    ASSERT_EQ('G', duckvep_translate_codon("ggg", STD)); /* case-insensitive */
    ASSERT_EQ('X', duckvep_translate_codon("ANG", STD)); /* unresolved consensus */
    ASSERT_EQ('A', duckvep_translate_codon("GCN", STD));
    ASSERT_EQ('T', duckvep_translate_codon("acn", STD));
    ASSERT_EQ('X', duckvep_translate_codon("AAN", STD));
    ASSERT_EQ('X', duckvep_translate_codon("A?G", STD));
    ASSERT_EQ('X', duckvep_translate_codon(NULL, STD));
    /* vertebrate mitochondrial (transl_table=2) — the four edits */
    ASSERT_EQ('W', duckvep_translate_codon("TGA", MITO));
    ASSERT_EQ('M', duckvep_translate_codon("ATA", MITO));
    ASSERT_EQ('*', duckvep_translate_codon("AGA", MITO));
    ASSERT_EQ('*', duckvep_translate_codon("AGG", MITO));
    ASSERT_EQ('M', duckvep_translate_codon("ATG", MITO)); /* start unchanged */
    /* One discriminating codon for every additional BioPerl/VEP table. */
    ASSERT_EQ('T', duckvep_translate_codon("CTA", (duckvep_codon_table_t)3));
    ASSERT_EQ('W', duckvep_translate_codon("TGA", (duckvep_codon_table_t)4));
    ASSERT_EQ('S', duckvep_translate_codon("AGA", (duckvep_codon_table_t)5));
    ASSERT_EQ('Q', duckvep_translate_codon("TAA", (duckvep_codon_table_t)6));
    ASSERT_EQ('N', duckvep_translate_codon("AAA", (duckvep_codon_table_t)9));
    ASSERT_EQ('C', duckvep_translate_codon("TGA", (duckvep_codon_table_t)10));
    ASSERT_EQ('*', duckvep_translate_codon("TGA", (duckvep_codon_table_t)11));
    ASSERT_EQ('S', duckvep_translate_codon("CTG", (duckvep_codon_table_t)12));
    ASSERT_EQ('G', duckvep_translate_codon("AGA", (duckvep_codon_table_t)13));
    ASSERT_EQ('Y', duckvep_translate_codon("TAA", (duckvep_codon_table_t)14));
    ASSERT_EQ('L', duckvep_translate_codon("TAG", (duckvep_codon_table_t)16));
    ASSERT_EQ('N', duckvep_translate_codon("AAA", (duckvep_codon_table_t)21));
    ASSERT_EQ('*', duckvep_translate_codon("TCA", (duckvep_codon_table_t)22));
    ASSERT_EQ('*', duckvep_translate_codon("TTA", (duckvep_codon_table_t)23));
    ASSERT_EQ('K', duckvep_translate_codon("AGG", (duckvep_codon_table_t)24));
    ASSERT_EQ('G', duckvep_translate_codon("TGA", (duckvep_codon_table_t)25));
    ASSERT_EQ('A', duckvep_translate_codon("CTG", (duckvep_codon_table_t)26));
    ASSERT_EQ('Q', duckvep_translate_codon("TAG", (duckvep_codon_table_t)27));
    ASSERT_EQ('Q', duckvep_translate_codon("TAA", (duckvep_codon_table_t)28));
    ASSERT_EQ('Y', duckvep_translate_codon("TAG", (duckvep_codon_table_t)29));
    ASSERT_EQ('E', duckvep_translate_codon("TAA", (duckvep_codon_table_t)30));
    ASSERT_EQ('W', duckvep_translate_codon("TGA", (duckvep_codon_table_t)31));
    ASSERT_EQ('X', duckvep_translate_codon("ATG", (duckvep_codon_table_t)8));
    ASSERT_EQ('S', duckvep_translate_codon("TCN", STD));
    PASS();
}

TEST codon_supported_set_matches_vep116_bioperl(void) {
    static const uint8_t supported[32] = {
        [1] = 1u, [2] = 1u, [3] = 1u, [4] = 1u, [5] = 1u, [6] = 1u,
        [9] = 1u, [10] = 1u, [11] = 1u, [12] = 1u, [13] = 1u,
        [14] = 1u, [16] = 1u, [21] = 1u, [22] = 1u, [23] = 1u,
        [24] = 1u, [25] = 1u, [26] = 1u, [27] = 1u, [28] = 1u,
        [29] = 1u, [30] = 1u, [31] = 1u
    };
    size_t id;
    for (id = 0u; id < sizeof supported; id++) {
        ASSERT_EQ(supported[id] != 0u,
                  duckvep_codon_table_supported((duckvep_codon_table_t)id));
    }
    ASSERT_EQ(0, duckvep_codon_table_supported((duckvep_codon_table_t)32));
    PASS();
}

TEST codon_change_known(void) {
    duckvep_codon_result_t r;
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS, duckvep_codon_change("AAA", "AAA", STD).change);
    r = duckvep_codon_change("AAA", "AAC", STD); /* K -> N */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, r.change);
    ASSERT_EQ('K', r.aa_ref); ASSERT_EQ('N', r.aa_alt);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, duckvep_codon_change("TAC", "TAA", STD).change); /* Y->* */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_LOST,   duckvep_codon_change("TAA", "TAC", STD).change); /* *->Y */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID,     duckvep_codon_change("ANG", "AAA", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, duckvep_codon_change("GCN", "ACN", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS, duckvep_codon_change("GCN", "GCT", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID, duckvep_codon_change("AAN", "ACN", STD).change);
    /* MT matters: TGG->TGA is stop_gained nuclear but synonymous in mito (TGA=W). */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, duckvep_codon_change("TGG", "TGA", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS,  duckvep_codon_change("TGG", "TGA", MITO).change);
    /* CGA->AGA: synonymous nuclear (R->R) but stop_gained in mito (AGA=*). */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS,  duckvep_codon_change("CGA", "AGA", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, duckvep_codon_change("CGA", "AGA", MITO).change);
    PASS();
}

struct kprop_codon { char ref[4]; char alt[4]; duckvep_codon_table_t table; };

static enum theft_alloc_res kprop_codon_alloc(struct theft *t, void *env, void **instance) {
    static const char bases5[5] = {'A', 'C', 'G', 'T', 'N'};
    struct kprop_codon *c = (struct kprop_codon *)calloc(1u, sizeof *c);
    int i;
    (void)env;
    if (c == NULL) return THEFT_ALLOC_ERROR;
    for (i = 0; i < 3; i++) {
        c->ref[i] = bases5[kprop_bounded(t, 5u)];
        c->alt[i] = bases5[kprop_bounded(t, 5u)];
    }
    c->table = (kprop_bounded(t, 2u) == 0u) ? STD : MITO;
    *instance = c;
    return THEFT_ALLOC_OK;
}
static void kprop_codon_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_codon_info = {
    .alloc = kprop_codon_alloc,
    .free  = kprop_codon_free,
};

/* The classification must be consistent with translation: exactly one change
 * bit, and each bit's defining relationship between aa_ref/aa_alt holds. */
static enum theft_trial_res prop_codon_change_consistent(struct theft *t, void *arg1) {
    const struct kprop_codon *c = (const struct kprop_codon *)arg1;
    char ar = duckvep_translate_codon(c->ref, c->table);
    char aa = duckvep_translate_codon(c->alt, c->table);
    duckvep_codon_result_t r = duckvep_codon_change(c->ref, c->alt, c->table);
    duckvep_codon_result_t prepared =
        duckvep_codon_change_prepared(c->ref, c->alt, c->table);
    (void)t;
    if (prepared.aa_ref != r.aa_ref || prepared.aa_alt != r.aa_alt ||
        prepared.change != r.change) return THEFT_TRIAL_FAIL;
    if (r.aa_ref != ar || r.aa_alt != aa) return THEFT_TRIAL_FAIL;
    if (ar == 'X' || aa == 'X') {
        return (r.change == (uint32_t)DUCKVEP_CODON_INVALID) ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
    }
    if (popcount_u32(r.change) != 1) return THEFT_TRIAL_FAIL;
    switch (r.change) {
        case DUCKVEP_CODON_SYNONYMOUS:  if (ar != aa) return THEFT_TRIAL_FAIL; break;
        case DUCKVEP_CODON_MISSENSE:    if (ar == aa || ar == '*' || aa == '*') return THEFT_TRIAL_FAIL; break;
        case DUCKVEP_CODON_STOP_GAINED: if (!(aa == '*' && ar != '*')) return THEFT_TRIAL_FAIL; break;
        case DUCKVEP_CODON_STOP_LOST:   if (!(ar == '*' && aa != '*')) return THEFT_TRIAL_FAIL; break;
        default: return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST codon_change_consistent_with_translation(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "codon change classification consistent with translation";
    cfg.prop1 = prop_codon_change_consistent;
    cfg.type_info[0] = &kprop_codon_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ===================================================================== *
 * Sequence-backed SNV codon edit application.
 *
 * Oracle: a codon-slice edit derived from the generated CDS position — copy the
 * reference codon from the generated CDS and replace exactly the projected base
 * after strand orientation — then classify via the already independently-tested
 * genetic-code table. This pins the hot-path contract between projection and the
 * codon classifier without claiming indel/haplotype support yet.
 * ===================================================================== */

char coding_test_comp(char b) {
    switch (b) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        default: return 'N';
    }
}
char coding_test_genomic_from_tx(char tx_base, int8_t strand) {
    return strand < 0 ? coding_test_comp(tx_base) : tx_base;
}

TEST coding_snv_from_cds_known_cases(void) {
    duckvep_coding_projection_t p;
    duckvep_coding_snv_result_t r;

    memset(&p, 0, sizeof p);
    p.cds_pos = 4u; p.codon_start_cds = 4u; p.protein_pos = 2u; p.codon_offset = 0u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"ATGGCTTGG", 9u, &p, p.cds_pos,
                                          'G', 'A', (int8_t)1, STD, &r));
    ASSERT_STR_EQ("GCT", r.ref_codon);
    ASSERT_STR_EQ("ACT", r.alt_codon);
    ASSERT_EQ('A', r.aa_ref); ASSERT_EQ('T', r.aa_alt);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, r.change);
    ASSERT_EQ(4u, r.cds_pos); ASSERT_EQ(2u, r.protein_pos);

    /* Negative strand: genomic C->A orients to transcript G->T, TGG->TGT. */
    memset(&p, 0, sizeof p);
    p.cds_pos = 3u; p.codon_start_cds = 1u; p.protein_pos = 1u; p.codon_offset = 2u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"TGG", 3u, &p, p.cds_pos,
                                          'C', 'A', (int8_t)-1, STD, &r));
    ASSERT_STR_EQ("TGG", r.ref_codon);
    ASSERT_STR_EQ("TGT", r.alt_codon);
    ASSERT_EQ('G', r.ref_base_tx); ASSERT_EQ('T', r.alt_base_tx);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, r.change);

    /* Codon table still matters after sequence-backed edit: TGG->TGA. */
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"TGG", 3u, &p, p.cds_pos,
                                          'G', 'A', (int8_t)1, STD, &r));
    ASSERT_STR_EQ("TGA", r.alt_codon);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, r.change);
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"TGG", 3u, &p, p.cds_pos,
                                          'G', 'A', (int8_t)1, MITO, &r));
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS, r.change);

    ASSERT_EQ(DUCKVEP_CODING_SNV_REF_MISMATCH,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p, p.cds_pos,
                                          'C', 'G', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID, r.change);
    ASSERT_EQ('X', r.aa_ref);

    p.codon_start_cds = 2u; p.cds_pos = 2u; p.codon_offset = 0u; p.protein_pos = 1u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAAAAA", 6u, &p, p.cds_pos,
                                          'A', 'G', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);

    p.codon_start_cds = 4u; p.cds_pos = 4u; p.codon_offset = 0u; p.protein_pos = 2u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_CODON_OUT_OF_RANGE,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAAAA", 5u, &p, p.cds_pos,
                                          'A', 'G', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);

    p.codon_start_cds = 1u; p.cds_pos = 3u; p.codon_offset = 2u; p.protein_pos = 1u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_BASE,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p, p.cds_pos,
                                          'A', 'N', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);

    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p, p.cds_pos,
                                          'A', 'G', (int8_t)0, STD, &r));
    p.protein_pos = 99u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p, p.cds_pos,
                                          'A', 'G', (int8_t)1, STD, &r));
    p.protein_pos = 1u;

    /* N in non-edited phase-padding positions is allowed but classifies invalid. */
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNA", 3u, &p, p.cds_pos,
                                          'A', 'G', (int8_t)1, STD, &r));
    ASSERT_STR_EQ("NNA", r.ref_codon);
    ASSERT_STR_EQ("NNG", r.alt_codon);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID, r.change);

    /* VEP edits CDS position one on a later phase-2 CDS, but genomic REF is
     * still the physical C at position three. Padding is not the REF oracle. */
    p.cds_pos = 1u; p.codon_offset = 0u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNC", 3u, &p, 3u,
                                          'C', 'A', 1, STD, &r));
    ASSERT_STR_EQ("NNC", r.ref_codon);
    ASSERT_STR_EQ("ANC", r.alt_codon);
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNC", 3u, &p, 3u,
                                          'G', 'T', -1, STD, &r));
    ASSERT_STR_EQ("ANC", r.alt_codon);
    ASSERT_EQ(DUCKVEP_CODING_SNV_REF_MISMATCH,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNC", 3u, &p, 3u,
                                          'T', 'A', 1, STD, &r));
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNC", 3u, &p, 0u,
                                          'C', 'A', 1, STD, &r));
    ASSERT_EQ(DUCKVEP_CODING_SNV_CODON_OUT_OF_RANGE,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNC", 3u, &p, 4u,
                                          'C', 'A', 1, STD, &r));
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_BASE,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNC", 3u, &p, 1u,
                                          'C', 'A', 1, STD, &r));
    PASS();
}

#define KPROP_MAX_CODING_CDS 90u

struct kprop_coding_snv {
    uint8_t seq[KPROP_MAX_CODING_CDS];
    size_t len;
    duckvep_coding_projection_t p;
    char genomic_ref, genomic_alt;
    char ref_tx, alt_tx;
    char ref_codon[4], alt_codon[4];
    duckvep_codon_table_t table;
    int8_t strand;
};

static enum theft_alloc_res kprop_coding_snv_alloc(struct theft *t, void *env, void **instance) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding_snv *s = (struct kprop_coding_snv *)calloc(1u, sizeof *s);
    uint32_t codons = (uint32_t)kprop_bounded(t, 20u) + 1u;
    uint32_t codon_idx = (uint32_t)kprop_bounded(t, codons);
    uint32_t offset = (uint32_t)kprop_bounded(t, 3u);
    size_t i;
    uint32_t alt_pick;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->len = (size_t)codons * 3u;
    if (s->len > KPROP_MAX_CODING_CDS) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < s->len; i++) s->seq[i] = (uint8_t)bases[kprop_bounded(t, 4u)];
    s->strand = (kprop_bounded(t, 2u) == 0u) ? (int8_t)1 : (int8_t)-1;
    s->table = (kprop_bounded(t, 2u) == 0u) ? STD : MITO;
    s->p.codon_start_cds = codon_idx * 3u + 1u;
    s->p.codon_offset = (uint8_t)offset;
    s->p.cds_pos = s->p.codon_start_cds + offset;
    s->p.protein_pos = codon_idx + 1u;
    s->ref_tx = (char)s->seq[s->p.cds_pos - 1u];
    alt_pick = (uint32_t)kprop_bounded(t, 3u) + 1u;
    s->alt_tx = bases[((s->ref_tx == 'A' ? 0u : s->ref_tx == 'C' ? 1u : s->ref_tx == 'G' ? 2u : 3u) + alt_pick) % 4u];
    s->genomic_ref = coding_test_genomic_from_tx(s->ref_tx, s->strand);
    s->genomic_alt = coding_test_genomic_from_tx(s->alt_tx, s->strand);
    for (i = 0u; i < 3u; i++) {
        s->ref_codon[i] = (char)s->seq[s->p.codon_start_cds - 1u + (uint32_t)i];
        s->alt_codon[i] = s->ref_codon[i];
    }
    s->ref_codon[3] = '\0'; s->alt_codon[3] = '\0';
    s->alt_codon[offset] = s->alt_tx;
    *instance = s;
    return THEFT_ALLOC_OK;
}
static void kprop_coding_snv_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_coding_snv_info = {
    .alloc = kprop_coding_snv_alloc,
    .free  = kprop_coding_snv_free,
};

static enum theft_trial_res prop_coding_snv_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding_snv *s = (const struct kprop_coding_snv *)arg1;
    duckvep_coding_snv_result_t r;
    duckvep_codon_result_t cr;
    (void)t;
    if (duckvep_coding_snv_from_cds(s->seq, s->len, &s->p, s->p.cds_pos, s->genomic_ref, s->genomic_alt,
                                    s->strand, s->table, &r) != DUCKVEP_CODING_SNV_OK) {
        return THEFT_TRIAL_FAIL;
    }
    cr = duckvep_codon_change(s->ref_codon, s->alt_codon, s->table);
    if (strcmp(r.ref_codon, s->ref_codon) != 0) return THEFT_TRIAL_FAIL;
    if (strcmp(r.alt_codon, s->alt_codon) != 0) return THEFT_TRIAL_FAIL;
    if (r.ref_base_tx != s->ref_tx || r.alt_base_tx != s->alt_tx) return THEFT_TRIAL_FAIL;
    if (r.aa_ref != cr.aa_ref || r.aa_alt != cr.aa_alt || r.change != cr.change) return THEFT_TRIAL_FAIL;
    if (r.cds_pos != s->p.cds_pos || r.codon_start_cds != s->p.codon_start_cds ||
        r.protein_pos != s->p.protein_pos || r.codon_offset != s->p.codon_offset) {
        return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST coding_snv_from_cds_matches_oracle_for_any_valid_snv(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence-backed SNV codon edit == codon-slice edit oracle";
    cfg.prop1 = prop_coding_snv_matches_oracle;
    cfg.type_info[0] = &kprop_coding_snv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}
