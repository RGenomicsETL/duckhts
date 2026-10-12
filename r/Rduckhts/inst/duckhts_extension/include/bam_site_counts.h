#ifndef DUCKHTS_BAM_SITE_COUNTS_H
#define DUCKHTS_BAM_SITE_COUNTS_H

#include <stddef.h>
#include <stdint.h>

#include <htslib/sam.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Somalier v0.3.4 obtains BAM evidence through hileup v0.1.0 at commit
 * 3133320d9660620a4a7f4da5647c85bf2b2738a6 (MIT, Brent Pedersen).
 * The hileup policy below preserves its encounter-order mate suppression.
 * Records and qnames are borrowed for the duration of one reduction call. */

typedef enum duckhts_bam_site_overlap_policy {
    DUCKHTS_BAM_SITE_OVERLAP_NONE = 0,
    DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0 = 1
} duckhts_bam_site_overlap_policy_t;

typedef enum duckhts_bam_site_status {
    DUCKHTS_BAM_SITE_OK = 0,
    DUCKHTS_BAM_SITE_INVALID_ARGUMENT,
    DUCKHTS_BAM_SITE_INVALID_PILEUP,
    DUCKHTS_BAM_SITE_COUNT_OVERFLOW,
    DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED
} duckhts_bam_site_status_t;

/* pos0 is zero-based. A and B are distinct uppercase A/C/G/T alleles in
 * lexical order, matching the canonical Somalier panel orientation. */
typedef struct duckhts_bam_site {
    hts_pos_t pos0;
    char allele_a;
    char allele_b;
} duckhts_bam_site_t;

typedef struct duckhts_bam_site_count_config {
    uint8_t min_baseq;
    duckhts_bam_site_overlap_policy_t overlap_policy;
} duckhts_bam_site_count_config_t;

typedef struct duckhts_bam_site_counts {
    uint32_t allele_a;
    uint32_t allele_b;
    uint32_t other;
} duckhts_bam_site_counts_t;

/* Slots borrow qname storage from the caller: from pileup records during one
 * duckhts_bam_site_count_pileup call, or from the caller's own entries when it
 * uses the table directly. The caller owns the hash table. */
typedef struct duckhts_bam_site_overlap_slot {
    const char *qname;
    uint64_t qname_hash;
    uint32_t qname_len;
    uint32_t generation;
} duckhts_bam_site_overlap_slot_t;

typedef struct duckhts_bam_site_overlap_scratch {
    duckhts_bam_site_overlap_slot_t *slots;
    size_t max_entries;
    size_t slot_capacity;
    uint32_t generation;
} duckhts_bam_site_overlap_scratch_t;

/* Returns the open-addressed table size required for max_entries active
 * qnames, or zero when 2 * max_entries + 1 overflows size_t. */
int duckhts_bam_site_overlap_slot_capacity(size_t max_entries,
                                           size_t *slot_capacity);

/* The qname table on its own, for callers that track mates across records.
 * The caller counts live entries and inserts no more than max_entries.
 *
 * begin empties the table in O(1) by moving to the next generation; it fails
 * when the slot table is smaller than slot_capacity(max_entries). find stores
 * the slot of a live qname in found_slot and returns 1, or returns 0. insert
 * returns 0 only when the table is full or qname_len exceeds UINT32_MAX.
 * remove frees a slot that find returned; it moves other live slots, so a
 * slot index is valid only until the next insert or remove. Every inserted
 * qname stays borrowed until it is removed or the table is begun again. */
uint64_t duckhts_bam_site_qname_hash(const char *qname, size_t qname_len);
int duckhts_bam_site_overlap_begin(duckhts_bam_site_overlap_scratch_t *scratch);
int duckhts_bam_site_overlap_find(const duckhts_bam_site_overlap_scratch_t *scratch,
                                  const char *qname,
                                  size_t qname_len,
                                  uint64_t hash,
                                  size_t *found_slot);
int duckhts_bam_site_overlap_insert(duckhts_bam_site_overlap_scratch_t *scratch,
                                    const char *qname,
                                    size_t qname_len,
                                    uint64_t hash);
void duckhts_bam_site_overlap_remove(duckhts_bam_site_overlap_scratch_t *scratch,
                                     size_t removed_slot);

/* MAPQ and general SAM-flag filtering belong to the reader adapter. This
 * kernel skips deletion/refskip entries, applies base quality, uppercases the
 * observed BAM base, and counts every non-A/B observation as other. Missing
 * QUAL (0xff) is accepted only when min_baseq is zero.
 *
 * overlap_scratch may be NULL for OVERLAP_NONE. Hileup overlap mode requires
 * a zero-initialized scratch descriptor and a zero-initialized slot table of
 * at least duckhts_bam_site_overlap_slot_capacity(max_entries) slots. Zero
 * max_entries is valid when no record needs tracking. counts is assigned only
 * on success. */
duckhts_bam_site_status_t duckhts_bam_site_count_pileup(
    const duckhts_bam_site_t *site,
    const bam_pileup1_t *pileup,
    size_t pileup_count,
    const duckhts_bam_site_count_config_t *config,
    duckhts_bam_site_overlap_scratch_t *overlap_scratch,
    duckhts_bam_site_counts_t *counts);

const char *duckhts_bam_site_status_string(duckhts_bam_site_status_t status);

#ifdef __cplusplus
}
#endif

#endif
