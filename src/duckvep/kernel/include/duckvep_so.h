/*
 * duckvep_so.h — Sequence Ontology consequence vocabulary (PUBLIC ABI).
 *
 * The kernel emits consequences as a bitset (duckvep_consequence_t.consequence_mask):
 * one bit per SO term. This header is the source of truth for the stable public
 * bit indices. VEP term/impact/rank/tier metadata is generated from the pinned
 * class-model extract and joined through conformance/data/so_bit_bindings.tsv.
 * Tests anchor on the indices, so they are ASSIGNED ONCE and never reordered;
 * new terms append at the next free index. The metadata ledger pins provenance,
 * not behavioral conformance; only the VEP differential can prove that.
 *
 * Impacts mirror Ensembl VEP's HIGH/MODERATE/LOW/MODIFIER ranking (see
 * https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).
 * SET by duckvep_annotate_tile today: the structural subset
 * (up/down/intron/utr/non-coding) plus the VEP-source-grounded SNV-point splice
 * subset (splice_region / splice_donor / splice_acceptor / splice_donor_5th_base /
 * splice_donor_region / splice_polypyrimidine_tract — from duckvep_splice_classify;
 * not yet measured against a VEP --gff differential) plus the codon-SNV subset
 * (synonymous / missense / stop_gained / stop_lost / stop_retained / start_lost),
 * narrow two-codon MNV missense, and sequence-backed frameshift / in-frame-deletion
 * indel slices, plus
 * full-transcript structural deletion/CNV-loss and duplication/CNV-gain terms.
 * DECLARED but NOT yet set: incomplete_terminal_codon and the transcript-biotype
 * terms NMD_transcript_variant / mature_miRNA_variant (need the distilled tx_flags
 * facts + side relations). See conformance/README for the current measured scope.
 */
#ifndef DUCKVEP_SO_H
#define DUCKVEP_SO_H

#include "duckvep_kernel.h" /* duckvep_impact_t */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Stable SO term bit indices. NEVER reorder; append new terms only. The set fits
 * in a uint64_t consequence_mask (room for 64 terms; 32 used). */
typedef enum duckvep_so_bit {
    /* --- structural (set today by the structural fusion) --- */
    DUCKVEP_SO_UPSTREAM_GENE              = 0,
    DUCKVEP_SO_DOWNSTREAM_GENE            = 1,
    DUCKVEP_SO_INTRON                     = 2,
    DUCKVEP_SO_5_PRIME_UTR                = 3,
    DUCKVEP_SO_3_PRIME_UTR                = 4,
    DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON = 5,
    DUCKVEP_SO_NON_CODING_TRANSCRIPT      = 6,
    DUCKVEP_SO_SPLICE_REGION              = 7,
    DUCKVEP_SO_SPLICE_DONOR               = 8,
    DUCKVEP_SO_SPLICE_ACCEPTOR            = 9,
    DUCKVEP_SO_SPLICE_DONOR_5TH_BASE      = 10,
    DUCKVEP_SO_CODING_SEQUENCE            = 11,
    /* --- coding substitutions (sequence-backed subset implemented) --- */
    DUCKVEP_SO_SYNONYMOUS                 = 12,
    DUCKVEP_SO_MISSENSE                   = 13,
    DUCKVEP_SO_STOP_GAINED                = 14,
    DUCKVEP_SO_STOP_LOST                  = 15,
    DUCKVEP_SO_START_LOST                 = 16,
    DUCKVEP_SO_STOP_RETAINED              = 17,
    DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON  = 18,
    /* --- additional splice terms (set by the splice-state classifier) --- */
    DUCKVEP_SO_SPLICE_DONOR_REGION        = 19,
    DUCKVEP_SO_SPLICE_POLYPYRIMIDINE      = 20,
    /* --- transcript-biotype terms (declared now; set when the import carries the
     * distilled biotype facts — DUCKVEP_TX_BIOTYPE_NMD / _MIRNA (duckvep_kernel.h) +
     * the mature-miRNA side relation, plus rule-table tier suppression) --- */
    DUCKVEP_SO_NMD_TRANSCRIPT             = 21,
    DUCKVEP_SO_MATURE_MIRNA               = 22,
    DUCKVEP_SO_FRAMESHIFT                 = 23,
    DUCKVEP_SO_INFRAME_DELETION           = 24,
    /* --- structural/CNV terms (structural full-span event geometry) --- */
    DUCKVEP_SO_TRANSCRIPT_ABLATION        = 25,
    DUCKVEP_SO_TRANSCRIPT_AMPLIFICATION   = 26,
    DUCKVEP_SO_FEATURE_ELONGATION         = 27,
    DUCKVEP_SO_FEATURE_TRUNCATION         = 28,
    DUCKVEP_SO_CODING_TRANSCRIPT          = 29,
    DUCKVEP_SO_INFRAME_INSERTION          = 30,
    DUCKVEP_SO_PROTEIN_ALTERING           = 31,
    DUCKVEP_SO_BIT_COUNT                  = 32
} duckvep_so_bit_t;

/* The single-bit mask for a term index. */
#define DUCKVEP_SO(bit) (UINT64_C(1) << (bit))

/* Stable lowercase SO term name (e.g. "intron_variant", "missense_variant"), or
 * NULL if `bit` is out of range. The adapter renders the VCF/Tab CSQ string from
 * these; the kernel never builds strings in the hot path. */
const char *duckvep_so_name(duckvep_so_bit_t bit);

/* Metadata copied from VEP's %OVERLAP_CONSEQUENCES class model. Rank is
 * severity order (smaller is more severe); tier controls evaluator suppression.
 * Out-of-range rank/tier -> UINT8_MAX. */
uint8_t duckvep_so_rank(duckvep_so_bit_t bit);
uint8_t duckvep_so_tier(duckvep_so_bit_t bit);

/* Impact of a single SO term. Out-of-range -> DUCKVEP_IMPACT_MODIFIER. */
duckvep_impact_t duckvep_so_bit_impact(duckvep_so_bit_t bit);

/* The maximum impact across all terms set in `mask` (VEP's "most severe" ranking
 * applied to one consequence row). An empty mask -> DUCKVEP_IMPACT_MODIFIER. */
duckvep_impact_t duckvep_so_impact(uint64_t mask);

/* The VEP impact label ("HIGH"/"MODERATE"/"LOW"/"MODIFIER") for an impact value.
 * Out-of-range -> "MODIFIER". A stable string in static storage (never freed). */
const char *duckvep_impact_name(duckvep_impact_t impact);

/* Render the SO terms set in `mask` as their lowercase names joined by `sep` (VEP
 * joins co-occurring consequences with '&'), in VEP severity-rank order, into `buf` as a
 * NUL-terminated string that never overruns `buflen`. Returns the number of chars the
 * full string WOULD occupy excluding the NUL (snprintf semantics), so the caller can
 * detect + grow on truncation. An empty mask yields "" and returns 0. This is the
 * OUTPUT-TIME render the adapter builds CSQ strings from (called once per emitted row,
 * NOT in the annotate hot loop) — the kernel still spells no SO strings while sweeping. */
size_t duckvep_so_render(uint64_t mask, char sep, char *buf, size_t buflen);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_SO_H */
