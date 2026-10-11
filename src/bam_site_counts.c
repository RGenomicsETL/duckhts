#include "include/bam_site_counts.h"

#include <ctype.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

_Static_assert(sizeof(uint32_t) * CHAR_BIT == 32,
               "BAM site counts require 32-bit uint32_t");

static int canonical_base(char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

static int checked_add_size(size_t left, size_t right, size_t *result) {
    if (left > SIZE_MAX - right) return 0;
    *result = left + right;
    return 1;
}

static int checked_multiply_size(size_t left, size_t right, size_t *result) {
    if (left != 0 && right > SIZE_MAX / left) return 0;
    *result = left * right;
    return 1;
}

static int record_has_observed_base(const bam1_t *record, int32_t qpos) {
    size_t data_bytes;
    size_t cigar_bytes;
    size_t sequence_bytes;
    size_t required;

    if (!record || !record->data || record->l_data < 0 ||
        record->core.l_qseq <= 0 || qpos < 0 || qpos >= record->core.l_qseq) {
        return 0;
    }

    data_bytes = (size_t)record->l_data;
    if (!checked_multiply_size((size_t)record->core.n_cigar,
                               sizeof(uint32_t), &cigar_bytes)) {
        return 0;
    }
    sequence_bytes = ((size_t)record->core.l_qseq + 1u) / 2u;

    required = (size_t)record->core.l_qname;
    if (!checked_add_size(required, cigar_bytes, &required) ||
        !checked_add_size(required, sequence_bytes, &required) ||
        !checked_add_size(required, (size_t)record->core.l_qseq, &required)) {
        return 0;
    }
    return required <= data_bytes;
}

static int record_qname(const bam1_t *record,
                        const char **qname,
                        size_t *qname_len) {
    size_t i;
    size_t storage_len;

    if (!record || !record->data || record->l_data < 0 ||
        record->core.l_qname == 0 ||
        (size_t)record->core.l_qname > (size_t)record->l_data) {
        return 0;
    }

    storage_len = (size_t)record->core.l_qname;
    *qname = bam_get_qname(record);
    for (i = 0; i < storage_len; i++) {
        if ((*qname)[i] == '\0') {
            *qname_len = i;
            return 1;
        }
    }
    return 0;
}

uint64_t duckhts_bam_site_qname_hash(const char *qname, size_t qname_len) {
    uint64_t hash = UINT64_C(14695981039346656037);
    size_t i;

    for (i = 0; i < qname_len; i++) {
        hash ^= (uint8_t)qname[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

int duckhts_bam_site_overlap_slot_capacity(size_t max_entries,
                                           size_t *slot_capacity) {
    if (!slot_capacity || max_entries > (SIZE_MAX - 1u) / 2u) return 0;
    *slot_capacity = max_entries == 0u ? 0u : max_entries * 2u + 1u;
    return 1;
}

static size_t overlap_probe_distance(size_t home,
                                     size_t position,
                                     size_t capacity) {
    return position >= home ? position - home : capacity - home + position;
}

int duckhts_bam_site_overlap_begin(duckhts_bam_site_overlap_scratch_t *scratch) {
    size_t required;

    if (!scratch ||
        !duckhts_bam_site_overlap_slot_capacity(scratch->max_entries,
                                                &required) ||
        scratch->slot_capacity < required ||
        scratch->slot_capacity > SIZE_MAX / sizeof(*scratch->slots) ||
        (required != 0u && !scratch->slots)) {
        return 0;
    }
    scratch->generation++;
    if (scratch->generation == 0u) {
        memset(scratch->slots, 0,
               scratch->slot_capacity * sizeof(*scratch->slots));
        scratch->generation = 1u;
    }
    return 1;
}

int duckhts_bam_site_overlap_find(const duckhts_bam_site_overlap_scratch_t *scratch,
                                  const char *qname,
                                  size_t qname_len,
                                  uint64_t hash,
                                  size_t *found_slot) {
    size_t index;
    size_t scanned;

    if (scratch->slot_capacity == 0u) return 0;
    index = (size_t)(hash % scratch->slot_capacity);
    for (scanned = 0u; scanned < scratch->slot_capacity; scanned++) {
        const duckhts_bam_site_overlap_slot_t *slot = &scratch->slots[index];
        if (slot->generation != scratch->generation) return 0;
        if (slot->qname_hash == hash && slot->qname_len == qname_len &&
            memcmp(slot->qname, qname, qname_len) == 0) {
            *found_slot = index;
            return 1;
        }
        index++;
        if (index == scratch->slot_capacity) index = 0u;
    }
    return 0;
}

int duckhts_bam_site_overlap_insert(duckhts_bam_site_overlap_scratch_t *scratch,
                                    const char *qname,
                                    size_t qname_len,
                                    uint64_t hash) {
    size_t index;
    size_t scanned;

    if (qname_len > UINT32_MAX || scratch->slot_capacity == 0u) return 0;
    index = (size_t)(hash % scratch->slot_capacity);
    for (scanned = 0u; scanned < scratch->slot_capacity; scanned++) {
        duckhts_bam_site_overlap_slot_t *slot = &scratch->slots[index];
        if (slot->generation != scratch->generation) {
            slot->qname = qname;
            slot->qname_hash = hash;
            slot->qname_len = (uint32_t)qname_len;
            slot->generation = scratch->generation;
            return 1;
        }
        index++;
        if (index == scratch->slot_capacity) index = 0u;
    }
    return 0;
}

void duckhts_bam_site_overlap_remove(duckhts_bam_site_overlap_scratch_t *scratch,
                                     size_t removed_slot) {
    size_t hole = removed_slot;
    size_t next = removed_slot;

    for (;;) {
        size_t home;

        next++;
        if (next == scratch->slot_capacity) next = 0u;
        if (scratch->slots[next].generation != scratch->generation) break;
        home = (size_t)(scratch->slots[next].qname_hash %
                        scratch->slot_capacity);
        if (overlap_probe_distance(home, hole, scratch->slot_capacity) <
            overlap_probe_distance(home, next, scratch->slot_capacity)) {
            scratch->slots[hole] = scratch->slots[next];
            hole = next;
        }
    }
    memset(&scratch->slots[hole], 0, sizeof(scratch->slots[hole]));
}

static int hileup_tracks_mate(const bam1_t *record, hts_pos_t site_pos0) {
    uint16_t non_primary = BAM_FSECONDARY | BAM_FSUPPLEMENTARY;

    return (record->core.flag & non_primary) == 0 &&
           bam_endpos(record) > record->core.mpos &&
           record->core.mpos <= site_pos0 &&
           record->core.tid == record->core.mtid &&
           record->core.pos <= record->core.mpos;
}

static duckhts_bam_site_status_t increment_count(uint32_t *count) {
    if (*count == UINT32_MAX) return DUCKHTS_BAM_SITE_COUNT_OVERFLOW;
    (*count)++;
    return DUCKHTS_BAM_SITE_OK;
}

duckhts_bam_site_status_t duckhts_bam_site_count_pileup(
    const duckhts_bam_site_t *site,
    const bam_pileup1_t *pileup,
    size_t pileup_count,
    const duckhts_bam_site_count_config_t *config,
    duckhts_bam_site_overlap_scratch_t *overlap_scratch,
    duckhts_bam_site_counts_t *counts) {
    duckhts_bam_site_counts_t local_counts = {0, 0, 0};
    size_t overlap_used = 0;
    size_t i;

    if (!site || !config || !counts || (pileup_count != 0 && !pileup)) {
        return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
    }
    if (site->pos0 < 0 || !canonical_base(site->allele_a) ||
        !canonical_base(site->allele_b) || site->allele_a >= site->allele_b) {
        return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
    }
    if (config->overlap_policy != DUCKHTS_BAM_SITE_OVERLAP_NONE &&
        config->overlap_policy != DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0) {
        return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
    }
#if SIZE_MAX > UINT32_MAX
    if (pileup_count > (size_t)UINT32_MAX) {
        return DUCKHTS_BAM_SITE_COUNT_OVERFLOW;
    }
#endif
    if (config->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0) {
        if (!duckhts_bam_site_overlap_begin(overlap_scratch)) {
            return DUCKHTS_BAM_SITE_INVALID_ARGUMENT;
        }
    }

    for (i = 0; i < pileup_count; i++) {
        const bam_pileup1_t *observation = &pileup[i];
        const bam1_t *record = observation->b;
        duckhts_bam_site_status_t status;
        uint8_t quality;
        char base;

        if (!record) return DUCKHTS_BAM_SITE_INVALID_PILEUP;
        if (observation->is_del || observation->is_refskip) continue;
        if (!record_has_observed_base(record, observation->qpos)) {
            return DUCKHTS_BAM_SITE_INVALID_PILEUP;
        }

        /* Pinned hileup order: gaps never consume a name; a same-contig name
         * match is removed before base-quality filtering. */
        if (config->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0) {
            const char *qname;
            size_t qname_len;
            size_t found_slot;
            uint64_t hash;

            if (!record_qname(record, &qname, &qname_len)) {
                return DUCKHTS_BAM_SITE_INVALID_PILEUP;
            }
            hash = duckhts_bam_site_qname_hash(qname, qname_len);
            if (record->core.tid == record->core.mtid &&
                duckhts_bam_site_overlap_find(overlap_scratch, qname, qname_len,
                                              hash, &found_slot)) {
                overlap_used--;
                duckhts_bam_site_overlap_remove(overlap_scratch, found_slot);
                continue;
            }
        }

        quality = bam_get_qual(record)[observation->qpos];
        if (quality == 0xff) {
            if (config->min_baseq != 0) continue;
        } else if (quality < config->min_baseq) {
            continue;
        }

        base = seq_nt16_str[bam_seqi(bam_get_seq(record), observation->qpos)];
        base = (char)toupper((unsigned char)base);
        if (base == site->allele_a) {
            status = increment_count(&local_counts.allele_a);
        } else if (base == site->allele_b) {
            status = increment_count(&local_counts.allele_b);
        } else {
            status = increment_count(&local_counts.other);
        }
        if (status != DUCKHTS_BAM_SITE_OK) return status;

        /* Only an observation that was counted can seed later suppression. */
        if (config->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0 &&
            hileup_tracks_mate(record, site->pos0)) {
            const char *qname;
            size_t qname_len;

            uint64_t hash;

            if (overlap_used == overlap_scratch->max_entries) {
                return DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED;
            }
            if (!record_qname(record, &qname, &qname_len)) {
                return DUCKHTS_BAM_SITE_INVALID_PILEUP;
            }
            hash = duckhts_bam_site_qname_hash(qname, qname_len);
            if (!duckhts_bam_site_overlap_insert(overlap_scratch, qname,
                                                 qname_len, hash)) {
                return DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED;
            }
            overlap_used++;
        }
    }

    *counts = local_counts;
    return DUCKHTS_BAM_SITE_OK;
}

const char *duckhts_bam_site_status_string(duckhts_bam_site_status_t status) {
    switch (status) {
        case DUCKHTS_BAM_SITE_OK:
            return "ok";
        case DUCKHTS_BAM_SITE_INVALID_ARGUMENT:
            return "invalid argument";
        case DUCKHTS_BAM_SITE_INVALID_PILEUP:
            return "invalid pileup";
        case DUCKHTS_BAM_SITE_COUNT_OVERFLOW:
            return "count overflow";
        case DUCKHTS_BAM_SITE_SCRATCH_EXHAUSTED:
            return "overlap scratch exhausted";
    }
    return "unknown BAM site count status";
}
