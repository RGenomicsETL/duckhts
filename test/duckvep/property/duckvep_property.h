/* Shared fixtures for the host-neutral Greatest/theft property families. */
#ifndef DUCKVEP_PROPERTY_H
#define DUCKVEP_PROPERTY_H

#include "duckvep_kernel.h"
#include "duckvep_phase.h"
#include "duckvep_sweep.h"
#include "duckvep_classify.h"
#include "duckvep_effect.h"
#include "duckvep_event.h"
#include "duckvep_so.h"
#include "duckvep_projection.h"
#include "duckvep_delta.h"
#include "duckvep_transcript_edit.h"
#include "duckvep_hgvs.h"
#include "duckvep_codon.h"
#include "duckvep_coding.h"
#include "duckvep_haplotype.h"
#include "duckvep_sequence_diff.h"
#include "duckvep_carriers.h"
#include "duckvep_haplotype_stream.h"
#include "duckvep_sv.h"
#include "duckvep_annotation_internal.h"
#include "duckvep_workspace_internal.h"

#include "greatest.h"
#include "theft.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define KPROP_DEFAULT_TRIALS 3000u
#define KPROP_DEFAULT_SEED   UINT64_C(0xd0c0ffee12345678)
#define STD  DUCKVEP_CODON_TABLE_STANDARD

/* Immutable zero transcript flags for small inline model fixtures. */
static const uint64_t k_zero_flags[8] = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};
#define KPROP_MAX_VARIANTS   64u

static inline uint64_t kprop_env_u64(const char *name, uint64_t dflt) {
    const char *v = getenv(name);
    char *end = NULL;
    unsigned long long parsed;
    if (v == NULL || v[0] == '\0') return dflt;
    parsed = strtoull(v, &end, 0);
    if (end == NULL || *end != '\0') return dflt;
    return (uint64_t)parsed;
}

static inline uint64_t kprop_bounded(struct theft *t, uint64_t bound) {
    if (bound <= 1u) return 0u;
    return theft_random_bits(t, 64) % bound; /* modulo bias is fine for tests */
}

#define KPROP_MAX_TX     32u
#define KPROP_SWEEP_HALO 100u
#define KPROP_COORD_MAX  2000u
#define KPROP_NCHROM     3u
#define KPROP_MAX_PAIRS  (KPROP_MAX_VARIANTS * KPROP_MAX_TX)

struct kprop_scene {
    duckvep_variant_batch_t    v;
    duckvep_transcript_model_t tx;
    uint32_t                   halo;  /* per-scene, so the window size is fuzzed */
    /* variant SoA */
    uint16_t *vchrom; uint32_t *vpos; uint32_t *vend; uint8_t *vkind;
    uint32_t *vroff;  uint16_t *vrlen; uint32_t *vaoff; uint16_t *valen; uint8_t *vbytes;
    /* transcript SoA */
    uint16_t *tchrom; uint32_t *tstart; uint32_t *tend; int8_t *tstrand; uint64_t *tflags;
    uint32_t *texoff; uint16_t *texcnt; uint32_t *tcds_s; uint32_t *tcds_e;
};

struct pair_collector { uint64_t *buf; size_t n; size_t cap; };

struct kprop_allele_sweep_scene {
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    uint32_t halo;
    uint16_t *vchrom;
    uint32_t *vpos;
    uint32_t *vend;
    uint8_t  *vkind;
    uint32_t *vroff;
    uint16_t *vrlen;
    uint32_t *vaoff;
    uint16_t *valen;
    uint8_t  *vbytes;
    size_t    vbytes_len;
    uint16_t *tchrom;
    uint32_t *tstart;
    uint32_t *tend;
    int8_t   *tstrand;
    uint64_t *tflags;
    uint32_t *texoff;
    uint16_t *texcnt;
    uint32_t *tcds_s;
    uint32_t *tcds_e;
};

#define KPROP_MAX_PROJ_EXONS 5u

struct kprop_proj_scene {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t       ex;
    uint16_t chrom; uint32_t tstart; uint32_t tend; int8_t strand; uint64_t flags;
    uint32_t exoff; uint16_t excnt; uint32_t cds_s; uint32_t cds_e;
    uint32_t es[KPROP_MAX_PROJ_EXONS];
    uint32_t ee[KPROP_MAX_PROJ_EXONS];
    uint32_t cs[KPROP_MAX_PROJ_EXONS];
    uint32_t ce[KPROP_MAX_PROJ_EXONS];
    int8_t phase[KPROP_MAX_PROJ_EXONS];
    int8_t end_phase[KPROP_MAX_PROJ_EXONS];
};

#define KPROP_HAPLO_MAX_EDITS 5u
#define KPROP_HAPLO_CDS_LEN   36u
#define KPROP_HAPLO_CAP       80u
#define KPROP_HAPLO_ALLELE_MAX 4u

struct kprop_haplo_case {
    uint8_t ref[KPROP_HAPLO_CDS_LEN + 1u];
    duckvep_haplotype_edit_t edits[KPROP_HAPLO_MAX_EDITS];
    uint8_t ref_alleles[KPROP_HAPLO_MAX_EDITS][3u];
    uint8_t alt_alleles[KPROP_HAPLO_MAX_EDITS][KPROP_HAPLO_ALLELE_MAX];
    size_t edit_count;
    int8_t transcript_strand;
};

/* Own the complete small-fixture views before an observer's scratch expires.
 * Oversized views fail the comparison instead of comparing a prefix. */
struct kprop_context_snapshot {
    duckvep_coding_context_t metadata;
    duckvep_coding_peptide_window_t window;
    int window_present;
    uint8_t cds[2][128], protein[2][128], local_peptide[2][128];
};

#define KPROP_MAX_CODONS 8u

struct kprop_coding {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t       ex;
    duckvep_sequence_pool_t    seq;
    duckvep_variant_batch_t    v;
    uint16_t chrom; uint32_t tstart; uint32_t tend; int8_t strand; uint64_t flags;
    uint32_t exoff; uint16_t excnt; uint32_t cds_s; uint32_t cds_e;
    uint32_t es; uint32_t ee; uint32_t ecds; uint32_t ecde; int8_t eph; int8_t eeph;
    uint8_t *cds; uint64_t cds_off0; uint32_t cds_lenv; uint8_t ctab;
    uint16_t vchrom; uint32_t vpos; uint32_t vend; uint8_t vkind;
    uint8_t abytes[16]; uint32_t roff; uint32_t aoff; uint16_t rlen; uint16_t alen;
    uint8_t expect_cds[64]; uint32_t expect_len; uint8_t expect_shape; uint8_t expect_region;
    uint8_t expect_protein_altering;
};

enum {
    KPROP_CDS_EDIT_START = 0u,
    KPROP_CDS_EDIT_BODY = 1u,
    KPROP_CDS_EDIT_STOP = 2u
};

extern struct theft_type_info kprop_allele_sweep_scene_info;

int popcount_u32(uint32_t x);

char coding_test_comp(char b);

uint8_t kprop_hgvs_complement(uint8_t base);

int pair_sink(uint32_t vi, uint32_t ti, void *ctx);

extern struct theft_type_info kprop_scene_info;

void kprop_proj_scene_finish(struct kprop_proj_scene *s);

int proj_brute_cdna_to_genomic(const struct kprop_proj_scene *s, uint32_t cdna,
                               uint32_t *genomic, uint32_t *exon_idx);

char coding_test_genomic_from_tx(char tx_base, int8_t strand);

char kprop_complement_base(char b);

int kprop_hgvs_protein_fact_replay(
    const duckvep_hgvs_protein_fact_t *fact,
    const uint8_t                     *reference,
    size_t                             reference_length,
    uint8_t                           *alternate,
    size_t                             alternate_capacity,
    size_t                            *alternate_length);

uint8_t haplo_test_variant_from_tx_base(char b, int reverse_complement);

int kprop_translate_full_oracle(const uint8_t *cds, size_t cds_len,
                                duckvep_codon_table_t table,
                                uint8_t *pep, size_t *pep_len);

int haplo_oracle_rebuild(const uint8_t *ref, size_t ref_len,
                         const duckvep_haplotype_edit_t *edits, size_t edit_count,
                         int8_t transcript_strand,
                         uint8_t *out, size_t out_cap, size_t *out_len,
                         int64_t *length_diff, uint32_t *flags);

char haplo_test_oriented_base(const uint8_t *seq, uint32_t len, uint32_t idx,
                              int reverse_complement);

void kprop_haplo_free(void *instance, void *env);

extern struct theft_type_info kprop_cds_edit_set_mnv_info;

int kprop_context_snapshot(const duckvep_coding_context_t *context,
                            struct kprop_context_snapshot *out);

int consequence_rows_equal(const duckvep_consequence_t *a,
                            const duckvep_consequence_t *b);

/* Greatest tests have external linkage across test-family translation units. */
#undef TEST
#define TEST enum greatest_test_res
#define DUCKVEP_PROPERTY_TEST(name) TEST name(void);
#include "duckvep_property_tests.def"
#undef DUCKVEP_PROPERTY_TEST

#endif
