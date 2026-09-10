#include "duckvep_property.h"

/* Exhaust actual slot permutations, independently of the reducer's equality
 * shortcut. A missing allele permits every allele in this finite alphabet. */
static int phase_next_permutation(uint8_t *order, uint16_t count) {
    int left = (int)count - 2;
    while (left >= 0 && order[left] >= order[left + 1]) left--;
    if (left < 0) return 0;
    int right = (int)count - 1;
    while (order[right] <= order[left]) right--;
    uint8_t tmp = order[left]; order[left] = order[right]; order[right] = tmp;
    for (int a = left + 1, b = (int)count - 1; a < b; a++, b--) {
        tmp = order[a]; order[a] = order[b]; order[b] = tmp;
    }
    return 1;
}

TEST phase_assignments_match_all_permitted_slot_permutations(void) {
    uint32_t cases = 0u;
    for (uint16_t ploidy = 1u; ploidy <= 5u; ploidy++) {
        uint32_t genotypes = UINT32_C(1) << (2u * ploidy);
        for (uint32_t genotype = 0u; genotype < genotypes; genotype++) {
            for (uint32_t bits = 0u; bits < (UINT32_C(1) << ploidy); bits++) {
                cases++;
                int32_t alleles[5];
                uint8_t phase[5], order[5], possible[5] = {0};
                uint8_t global = 0u;
                uint16_t unphased = 0u;
                duckvep_phase_summary_t summary = {0};
                for (uint16_t i = 0u; i < ploidy; i++) {
                    alleles[i] = (int32_t)((genotype >> (2u * i)) & 3u) - 1;
                    phase[i] = (uint8_t)((bits >> i) & 1u);
                    order[i] = (uint8_t)i;
                    if (!phase[i]) unphased++;
                    global |= alleles[i] < 0 ? 7u : (uint8_t)(1u << alleles[i]);
                    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_observe(&summary, alleles[i], phase[i]));
                }
                do {
                    int allowed = 1;
                    for (uint16_t i = 0u; i < ploidy; i++)
                        if (phase[i] && order[i] != i) allowed = 0;
                    if (!allowed) continue;
                    for (uint16_t i = 0u; i < ploidy; i++) {
                        int32_t value = alleles[order[i]];
                        possible[i] |= value < 0 ? 7u : (uint8_t)(1u << value);
                    }
                } while (phase_next_permutation(order, ploidy));
                int invariant = ploidy == 1u || (global & (global - 1u)) == 0u;
                uint16_t called_before = 0u;
                for (uint16_t i = 0u; i < ploidy; i++) {
                    duckvep_phase_assignment_t out;
                    int resolved = invariant || phase[i] || unphased == 1u ||
                        (possible[i] & (possible[i] - 1u)) == 0u;
                    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_assign(&summary, i + 1u, called_before,
                        alleles[i], phase[i], DUCKVEP_PHASE_STRICT, &out));
                    ASSERT_EQ(resolved ? i + 1u : 0u, out.lane);
                    ASSERT_EQ(invariant ? DUCKVEP_PHASE_ALL_SETS :
                        resolved ? DUCKVEP_PHASE_SET : DUCKVEP_PHASE_UNRESOLVED, out.scope);
                    ASSERT_EQ(alleles[i] < 0 ? DUCKVEP_PHASE_MISSING :
                        resolved ? DUCKVEP_PHASE_CALLED : DUCKVEP_PHASE_UNPHASED, out.status);
                    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_assign(&summary, i + 1u, called_before,
                        alleles[i], phase[i], DUCKVEP_PHASE_VEP116_COMPAT, &out));
                    ASSERT_EQ(alleles[i] < 0 ? 0u : called_before + 1u, out.lane);
                    ASSERT_EQ(DUCKVEP_PHASE_ALLELE_SLOT, out.scope);
                    duckvep_phase_call_status_t expected_status = alleles[i] < 0
                        ? DUCKVEP_PHASE_MISSING : DUCKVEP_PHASE_CALLED;
                    ASSERT_EQ(expected_status, out.status);
                    if (alleles[i] >= 0) called_before++;
                }
            }
        }
    }
    ASSERT_EQ(37448u, cases);
    PASS();
}

TEST raw_gt_source_spelling_and_file_slots_are_distinct(void) {
    static const struct {
        const char *gt;
        uint16_t ploidy;
        uint8_t missing;
        uint32_t slots, first, second;
        duckvep_raw_gt_disposition_t disposition;
    } cases[] = {
        {"0|1", 2, 0, 2, 0, 1, DUCKVEP_RAW_GT_RETAINED},
        {"|0|1", 2, 0, 3, 0, 0, DUCKVEP_RAW_GT_RETAINED},
        {"0/1|2", 3, 0, 2, 0, 2, DUCKVEP_RAW_GT_RETAINED},
        {"0|1/2", 3, 0, 2, 0, 1, DUCKVEP_RAW_GT_RETAINED},
        {".|1", 2, 1, 1, 1, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"1|.", 2, 1, 1, 1, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"./.", 2, 1, 0, UINT32_MAX, UINT32_MAX, DUCKVEP_RAW_GT_OMITTED_EMPTY},
        {"|./.", 2, 1, 2, 0, 0, DUCKVEP_RAW_GT_RETAINED},
        {"./.|1", 3, 1, 2, 0, 1, DUCKVEP_RAW_GT_RETAINED},
        {"|.", 1, 1, 1, 0, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"1", 1, 0, 1, 1, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"0", 1, 0, 0, UINT32_MAX, UINT32_MAX, DUCKVEP_RAW_GT_OMITTED_REFERENCE},
        {"0|0/0", 3, 0, 0, UINT32_MAX, UINT32_MAX, DUCKVEP_RAW_GT_OMITTED_REFERENCE},
        {"|0|0", 2, 0, 3, 0, 0, DUCKVEP_RAW_GT_RETAINED},
        {"01/2", 2, 0, 2, 1, 2, DUCKVEP_RAW_GT_RETAINED},
        {"/1|0", 2, 0, 2, 0, 0, DUCKVEP_RAW_GT_RETAINED}
    };
    for (size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); i++) {
        duckvep_raw_gt_t parsed;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
            (const uint8_t *)cases[i].gt, strlen(cases[i].gt), 2u, &parsed));
        ASSERT_EQ(cases[i].ploidy, parsed.source_ploidy);
        ASSERT_EQ(cases[i].missing, parsed.source_has_missing);
        ASSERT_EQ(cases[i].slots, parsed.parsed_slots);
        ASSERT_EQ(cases[i].first, parsed.allele_index[0]);
        ASSERT_EQ(cases[i].second, parsed.allele_index[1]);
        ASSERT_EQ(cases[i].disposition, parsed.disposition);
    }
    static const char *invalid[] = {"|", "/", "||1", "1|", "1.", "-1", "0//1", "x", "1\\2"};
    duckvep_raw_gt_t parsed, zero = {0};
    memset(&zero, 0, sizeof(zero));
    for (size_t i = 0u; i < sizeof(invalid) / sizeof(invalid[0]); i++) {
        memset(&parsed, 0xff, sizeof(parsed));
        ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_SYNTAX, duckvep_phase_parse_vep116_raw(
            (const uint8_t *)invalid[i], strlen(invalid[i]), 2u, &parsed));
        ASSERT_MEM_EQ(&zero, &parsed, sizeof(zero));
    }
    ASSERT_EQ(DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE, duckvep_phase_parse_vep116_raw(
        (const uint8_t *)"0|3", 3u, 2u, &parsed));
    ASSERT_MEM_EQ(&zero, &parsed, sizeof(zero));
    ASSERT_EQ(DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE, duckvep_phase_parse_vep116_raw(
        (const uint8_t *)"2147483648", 10u, INT32_MAX, &parsed));
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_ARG, duckvep_phase_parse_vep116_raw(NULL, 1u, 2u, &parsed));
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_ARG, duckvep_phase_parse_vep116_raw(
        (const uint8_t *)"0", 1u, UINT32_MAX, &parsed));
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_ARG, duckvep_phase_parse_vep116_raw(
        (const uint8_t *)"0", 1u, 2u, NULL));
    const uint8_t embedded_nul[] = {'0', '|', 0};
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_SYNTAX, duckvep_phase_parse_vep116_raw(
        embedded_nul, sizeof(embedded_nul), 2u, &parsed));
    uint8_t *many = malloc(2u * (size_t)UINT16_MAX + 1u);
    ASSERT(many != NULL);
    for (size_t i = 0u; i < 2u * (size_t)UINT16_MAX + 1u; i++) many[i] = i % 2u ? '|' : '0';
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
        many, 2u * (size_t)UINT16_MAX - 1u, 0u, &parsed));
    ASSERT_EQ(UINT16_MAX, parsed.source_ploidy);
    ASSERT_EQ(DUCKVEP_RAW_GT_PLOIDY_LIMIT, duckvep_phase_parse_vep116_raw(
        many, 2u * (size_t)UINT16_MAX + 1u, 0u, &parsed));
    ASSERT_MEM_EQ(&zero, &parsed, sizeof(zero));
    free(many);
    PASS();
}

TEST phase_ploidy_and_invalid_observations_are_checked(void) {
    duckvep_phase_summary_t summary = {0}, before;
    duckvep_phase_assignment_t out;
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG,
        duckvep_phase_assign(&summary, 1u, 0u, 0, 0u, DUCKVEP_PHASE_STRICT, &out));
    for (uint32_t i = 0u; i < UINT16_MAX; i++)
        ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_observe(&summary, INT32_MAX, 0u));
    before = summary;
    ASSERT_EQ(DUCKVEP_PHASE_PLOIDY_LIMIT, duckvep_phase_observe(&summary, 1, 1u));
    ASSERT_EQ(0, memcmp(&before, &summary, sizeof summary));
    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_assign(&summary, UINT16_MAX, UINT16_MAX - 1u,
        INT32_MAX, 0u, DUCKVEP_PHASE_STRICT, &out));
    ASSERT_EQ(UINT16_MAX, out.lane);
    ASSERT_EQ(DUCKVEP_PHASE_ALL_SETS, out.scope);
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_observe(&summary, -2, 0u));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_observe(&summary, 0, 2u));
    ASSERT_EQ(0, memcmp(&before, &summary, sizeof summary));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_observe(NULL, 0, 0u));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_assign(&summary, 0u, 0u,
        INT32_MAX, 0u, DUCKVEP_PHASE_STRICT, &out));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_assign(&summary, 1u, 0u,
        INT32_MAX, 0u, (duckvep_phase_policy_t)99, &out));
    PASS();
}

struct carrier_test_pool {
    duckvep_carrier_transcript_t transcripts[3];
    duckvep_carrier_call_t calls[16];
    duckvep_carrier_prefix_t prefixes[64];
    uint32_t active[3];
    duckvep_carrier_bucket_t tx_index[8], call_index[32], prefix_index[128];
};

static duckvep_carrier_buffers_t carrier_test_buffers(struct carrier_test_pool *pool) {
    duckvep_carrier_buffers_t b = {pool->transcripts, pool->calls, pool->prefixes,
        pool->active, pool->tx_index, pool->call_index, pool->prefix_index,
        3u, 16u, 64u, 8u, 32u, 128u};
    return b;
}

static duckvep_carrier_key_t carrier_test_key(uint32_t sample) {
    duckvep_carrier_key_t key = {sample, sample % 3u == 2u ? -5 : 0,
        (uint16_t)(sample % 4u + 1u), 4u, (uint8_t)(sample % 3u != 0u)};
    return key;
}

struct haplotype_stream_scene {
    duckvep_transcript_model_t model;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t sequences;
    uint16_t chrom[64], exon_count[64];
    uint32_t starts[64], ends[64], offsets[64], cdna_starts[64], cdna_ends[64], lengths[64];
    uint64_t sequence_offsets[64];
    int8_t strands[64];
    uint8_t reference[12];
    struct carrier_test_pool carrier_pool;
    duckvep_haplotype_stored_event_t events[8];
    duckvep_haplotype_projection_t projections[16];
    duckvep_haplotype_contributor_t contributors[8];
    duckvep_haplotype_edit_t edits[8];
    uint64_t edit_event_ids[9]; /* Last slot is a provenance-write canary. */
    duckvep_haplotype_block_t blocks[8];
    duckvep_carrier_event_t leaf_events[8];
    uint8_t alleles[128], cds[128], protein[128], reference_protein[128];
    uint8_t reference_coding_protein[128];
    duckvep_haplotype_stream_buffers_t buffers;
    duckvep_haplotype_stream_t stream;
};

static void haplotype_stream_scene_prepare(struct haplotype_stream_scene *f, uint32_t count) {
    memset(f, 0, sizeof(*f));
    memset(f->reference, 'A', sizeof(f->reference));
    for (uint32_t i = 0u; i < count; i++) {
        f->starts[i] = 100u + (i / 2u) * 12u + (i % 2u) * 4u;
        f->ends[i] = f->starts[i] + 11u;
        f->strands[i] = 1;
        f->offsets[i] = i;
        f->exon_count[i] = 1u;
        f->cdna_starts[i] = 1u;
        f->cdna_ends[i] = f->lengths[i] = 12u;
    }
    f->model.transcript_count = count; f->model.chrom_id = f->chrom;
    f->model.start1 = f->model.cds_start1 = f->starts;
    f->model.end1 = f->model.cds_end1 = f->ends;
    f->model.strand = f->strands; f->model.exon_offset = f->offsets;
    f->model.exon_count = f->exon_count;
    f->exons.exon_count = count; f->exons.start1 = f->starts; f->exons.end1 = f->ends;
    f->exons.cdna_start1 = f->cdna_starts; f->exons.cdna_end1 = f->cdna_ends;
    f->sequences.transcript_count = count; f->sequences.cds_bytes = f->reference;
    f->sequences.cds_bytes_len = sizeof(f->reference);
    f->sequences.cds_offset = f->sequence_offsets; f->sequences.cds_length = f->lengths;
    f->buffers = (duckvep_haplotype_stream_buffers_t){
        .carriers = carrier_test_buffers(&f->carrier_pool),
        .events = f->events, .projections = f->projections, .alleles = f->alleles,
        .event_capacity = 8u, .projection_capacity = 16u, .allele_capacity = sizeof(f->alleles),
        .leaf_events = f->leaf_events, .contributors = f->contributors, .edits = f->edits,
        .blocks = f->blocks, .edit_event_ids = f->edit_event_ids,
        .leaf_capacity = 8u, .edit_capacity = 8u, .cds = f->cds, .protein = f->protein,
        .cds_capacity = sizeof(f->cds), .protein_capacity = sizeof(f->protein),
        .reference_protein = f->reference_protein,
        .reference_coding_protein = f->reference_coding_protein,
        .reference_protein_capacity = sizeof(f->reference_protein)};
}

static duckvep_haplotype_stream_status_t haplotype_test_begin_candidates(
    duckvep_haplotype_stream_t *s, const duckvep_haplotype_source_t *source,
    const uint32_t *candidates, size_t count) {
    duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_begin(s, source);
    for (size_t i = 0u; status == DUCKVEP_HAPLOTYPE_STREAM_OK && i < count; i++)
        status = duckvep_haplotype_stream_project(s, candidates[i]);
    return status;
}

/* A complete, fully phased test GT has one ALT lane and REF elsewhere.
 * Exercise the same call interpretation used by the SQL adapter. */
static duckvep_haplotype_stream_status_t haplotype_test_push_called(
    duckvep_haplotype_stream_t *s, const duckvep_carrier_key_t *key) {
    if (s->error) return s->error;
    if (!s->have_current) return DUCKVEP_HAPLOTYPE_STREAM_OK;
    int32_t gt[8] = {0};
    uint8_t phase[8] = {1u,1u,1u,1u,1u,1u,1u,1u};
    if (!key->lane || key->lane > key->ploidy || key->ploidy > 8u)
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    gt[key->lane - 1u] = 1;
    duckvep_haplotype_phase_set_t set = {key->phase_set, key->phase_set_present};
    duckvep_haplotype_call_t call = {gt, phase, key->sample_index, 1u,
        key->ploidy, set, DUCKVEP_PHASE_STRICT};
    const duckvep_haplotype_stored_event_t *event = &s->buffers.events[s->current_event];
    for (uint32_t i = 0u; i < event->projection_count; i++) {
        size_t at = (event->projection_begin + i) % s->buffers.projection_capacity;
        duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_push_call(s,
            s->buffers.projections[at].transcript_index, &call, &set, 1u);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

TEST haplotype_record_order_matches_red_black_insertions(void) {
    /* Independent pointer-tree insertion and preorder traversal. Sorted keys
     * append on the right spine; fix-up uses rotations and uncle recolouring. */
    struct node { unsigned parent, left, right, red; } nodes[2049] = {{0}};
    unsigned root = 0u, stack[2048];
    ASSERT_EQ(0u, duckvep_haplotype_record_order(0u, 1u));
    ASSERT_EQ(0u, duckvep_haplotype_record_order(1u, 0u));
    ASSERT_EQ(0u, duckvep_haplotype_record_order(1u, 2u));
    for (unsigned n = 1u; n <= 2048u; n++) {
        nodes[n].parent = n - 1u; nodes[n].red = 1u;
        if (n == 1u) root = n;
        else nodes[n - 1u].right = n;
        unsigned current = n;
        while (nodes[nodes[current].parent].red) {
            unsigned parent = nodes[current].parent, grand = nodes[parent].parent;
            ASSERT_EQ(parent, nodes[grand].right);
            unsigned uncle = nodes[grand].left;
            if (nodes[uncle].red) {
                nodes[parent].red = nodes[uncle].red = 0u;
                nodes[grand].red = 1u; current = grand;
            } else {
                nodes[parent].red = 0u; nodes[grand].red = 1u;
                nodes[grand].right = nodes[parent].left;
                if (nodes[parent].left) nodes[nodes[parent].left].parent = grand;
                unsigned above = nodes[grand].parent;
                nodes[parent].parent = above;
                if (!above) root = parent;
                else if (nodes[above].left == grand) nodes[above].left = parent;
                else nodes[above].right = parent;
                nodes[parent].left = grand; nodes[grand].parent = parent;
            }
        }
        nodes[root].red = 0u;
        size_t pending = 1u; stack[0] = root;
        uint64_t rank = 0u;
        while (pending) {
            unsigned at = stack[--pending];
            ASSERT_EQ(++rank, duckvep_haplotype_record_order(n, at));
            if (nodes[at].right) stack[pending++] = nodes[at].right;
            if (nodes[at].left) stack[pending++] = nodes[at].left;
        }
        ASSERT_EQ(n, rank);
    }
    for (unsigned bit = 1u; bit < 64u; bit++) {
        uint64_t count = UINT64_MAX >> (63u - bit);
        uint64_t ranks[64];
        ASSERT_EQ(count, duckvep_haplotype_record_order(count, count));
        for (unsigned i = 0u; i <= bit; i++) {
            ranks[i] = duckvep_haplotype_record_order(count, UINT64_C(1) << i);
            ASSERT(ranks[i] && ranks[i] <= count);
            for (unsigned j = 0u; j < i; j++) ASSERT(ranks[j] != ranks[i]);
        }
    }
    PASS();
}

TEST haplotype_record_plan_closes_transcripts_not_later_record_ends(void) {
    uint16_t chrom[] = {0u, 0u, 0u, 1u};
    uint32_t starts[] = {10u, 18u, 40u, 5u}, ends[] = {20u, 30u, 50u, 9u};
    duckvep_transcript_model_t model = {.chrom_id = chrom, .start1 = starts,
        .end1 = ends, .transcript_count = 4u};
    duckvep_haplotype_record_plan_t plan;
    uint64_t buffer, ordinal;
    ASSERT(duckvep_haplotype_record_plan_init(&plan, &model));
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 1u, 9u, &buffer, &ordinal));
    ASSERT_EQ(0u, buffer); ASSERT_EQ(0u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 9u, 10u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(1u, ordinal); ASSERT_EQ(30u, plan.end1);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 29u, 100u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(2u, ordinal); ASSERT_EQ(30u, plan.end1);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 31u, 39u, &buffer, &ordinal));
    ASSERT_EQ(0u, buffer); ASSERT_EQ(0u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 40u, 40u, &buffer, &ordinal));
    ASSERT_EQ(2u, buffer); ASSERT_EQ(1u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 1u, 5u, 5u, &buffer, &ordinal));
    ASSERT_EQ(3u, buffer); ASSERT_EQ(1u, ordinal);
    ASSERT_FALSE(duckvep_haplotype_record_plan_next(&plan, 1u, 4u, 5u, &buffer, &ordinal));
    ASSERT_EQ(0u, buffer); ASSERT_EQ(0u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_init(&plan, &model));
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 9u, 41u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(50u, plan.end1);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 35u, 35u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(2u, ordinal);
    starts[1] = 1u;
    ASSERT_FALSE(duckvep_haplotype_record_plan_init(&plan, &model));
    PASS();
}

TEST haplotype_record_plan_matches_dense_closure(void) {
    uint64_t seeds[] = {173u, 20260906u};
    for (unsigned seed = 0u; seed < 2u; seed++) {
        uint64_t state = seeds[seed];
        for (unsigned trial = 0u; trial < 4096u; trial++) {
            uint16_t chrom[32];
            uint32_t starts[32], ends[32];
            for (unsigned i = 0u; i < 32u; i++) {
                state = state * UINT64_C(6364136223846793005) + 1u;
                chrom[i] = (uint16_t)(i / 16u);
                starts[i] = 1u + (i % 16u) * 16u;
                ends[i] = starts[i] + (uint32_t)((state >> 32u) % 32u);
            }
            duckvep_transcript_model_t model = {.chrom_id = chrom, .start1 = starts,
                .end1 = ends, .transcript_count = 32u};
            duckvep_haplotype_record_plan_t plan;
            ASSERT(duckvep_haplotype_record_plan_init(&plan, &model));
            uint64_t expected_buffer = 0u, expected_ordinal = 0u;
            for (uint16_t region = 0u; region < 2u; region++) {
                uint32_t start = 1u, expected_end = 0u;
                for (unsigned record = 0u; record < 64u; record++) {
                    state = state * UINT64_C(6364136223846793005) + 1u;
                    start += (uint32_t)((state >> 32u) % 8u);
                    uint32_t end = start + (uint32_t)((state >> 48u) % 64u);
                    if (record % 13u == 0u) end += 200u;
                    if (start > expected_end) {
                        uint32_t reachable = end, previous;
                        int overlaps = 0;
                        do {
                            previous = reachable;
                            for (unsigned i = 0u; i < 32u; i++) {
                                if (chrom[i] != region || starts[i] > reachable || ends[i] < start) continue;
                                overlaps = 1;
                                if (ends[i] > reachable) reachable = ends[i];
                            }
                        } while (reachable != previous);
                        expected_end = overlaps ? reachable : 0u;
                        if (overlaps) { expected_buffer++; expected_ordinal = 0u; }
                    }
                    uint64_t buffer, ordinal;
                    ASSERT(duckvep_haplotype_record_plan_next(&plan, region, start, end, &buffer, &ordinal));
                    if (expected_end) expected_ordinal++;
                    ASSERT_EQ(expected_end ? expected_buffer : 0u, buffer);
                    ASSERT_EQ(expected_end ? expected_ordinal : 0u, ordinal);
                    ASSERT_EQ(expected_end, plan.end1);
                }
            }
        }
    }
    PASS();
}

TEST haplotype_stream_raw_record_selection_and_tie_order(void) {
    for (unsigned planned = 0u; planned < 2u; planned++) {
        for (unsigned shadowed = 0u; shadowed < 2u; shadowed++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            duckvep_raw_gt_t call;
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                (const uint8_t *)"1|1", 3u, 1u, &call));
            for (unsigned i = 0u; i < 2u; i++) {
                duckvep_haplotype_source_t source = {i + 1u, (const uint8_t *)"A",
                    (const uint8_t *)(i ? "G" : "C"), 100u, 0u, 1u, 1u, 1u, 1u,
                    planned ? 2u - i : 0u};
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                    &f.stream, 0u, 0u, &call, (uint8_t)!(i && shadowed)));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(planned || shadowed ? 'C' : 'G', leaf.cds[0]);
            ASSERT_EQ(shadowed ? 1u : 2u, leaf.edit_count);
            ASSERT_EQ(2u, leaf.contributor_count);
            ASSERT_EQ(shadowed ? DUCKVEP_CDS_EDIT_SOURCE_SHADOWED : DUCKVEP_CDS_EDIT_OK,
                leaf.contributors[1].projection_status);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
        }
    }
    PASS();
}

TEST haplotype_stream_raw_record_slots_keep_conditional_deletion_provenance(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    duckvep_raw_gt_t calls[2];
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
        (const uint8_t *)".|1", 3u, 1u, &calls[0]));
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
        (const uint8_t *)"1|1", 3u, 1u, &calls[1]));
    duckvep_haplotype_source_t sources[] = {
        {11u, (const uint8_t *)"A", (const uint8_t *)"C", 101u, 0u, 1u, 1u, 1u, 1u, 0u},
        {11u, (const uint8_t *)"A", (const uint8_t *)"", 101u, 0u, 1u, 0u, UINT32_MAX, 1u, 0u},
        {12u, (const uint8_t *)"A", (const uint8_t *)"G", 104u, 0u, 1u, 1u, 1u, 1u, 0u}
    };
    for (size_t i = 0; i < 3u; i++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &sources[i]));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
            s, 0u, 7u, &calls[i == 2u], 1u));
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    unsigned lanes = 0u;
    while (duckvep_haplotype_stream_next(s, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(2u, leaf.contributor_count);
        ASSERT_EQ(2u, leaf.edit_count);
        ASSERT_EQ(11u, leaf.edit_event_ids[0]); ASSERT_EQ(12u, leaf.edit_event_ids[1]);
        const duckvep_carrier_call_t *carrier = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
        ASSERT(carrier != NULL); ASSERT_EQ(7u, carrier->key.sample_index);
        ASSERT_EQ(2u, carrier->key.ploidy); ASSERT_EQ(0u, carrier->key.phase_set_present);
        ASSERT_EQ(1u, leaf.carriers.call_count);
        unsigned lane = carrier->key.lane;
        ASSERT(lane == 1u || lane == 2u); lanes |= 1u << lane;
        ASSERT_EQ(lane == 1u ? 0 : -1, leaf.nominal_length_diff);
        const char *expected = lane == 1u ? "ACAAGAAAAAAA" : "AAAGAAAAAAA";
        ASSERT_EQ(strlen(expected), leaf.cds_length);
        ASSERT_MEM_EQ(expected, leaf.cds, leaf.cds_length);
        ASSERT_EQ(lane == 1u ? 1u : UINT32_MAX, leaf.contributors[0].source.allele_index);
        ASSERT_EQ(lane == 1u ? 1u : 0u, leaf.contributors[0].source.alt_len);
        ASSERT_EQ(lane == 1u ? DUCKVEP_CARRIER_CALLED : 0u,
            leaf.contributors[0].evidence_flags & DUCKVEP_CARRIER_CALLED);
        ASSERT(leaf.evidence_flags & DUCKVEP_CARRIER_CONDITIONAL);
        ASSERT(leaf.evidence_flags & DUCKVEP_CARRIER_MISSING);
    }
    ASSERT_EQ(6u, lanes);
    ASSERT_EQ(2u, s->completed_leaves); ASSERT_EQ(3u, s->projected_events);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    PASS();
}

TEST haplotype_stream_raw_record_reference_observations_require_matching_cds(void) {
    for (unsigned mismatch = 0u; mismatch < 2u; mismatch++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        const uint8_t *ref = (const uint8_t *)(mismatch ? "CAA" : "AAA");
        duckvep_haplotype_source_t source = {1u, ref, ref, 100u, 0u, 3u, 3u, 0u, 1u, 0u};
        duckvep_raw_gt_t call;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw((const uint8_t *)".", 1u, 1u, &call));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(2u, leaf.carriers.call_count); ASSERT_EQ(1u, leaf.contributor_count);
        ASSERT_EQ(0u, leaf.edit_count); ASSERT_EQ(0u, leaf.block_count);
        ASSERT_EQ(0, leaf.nominal_length_diff);
        ASSERT_EQ(DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_CONDITIONAL, leaf.evidence_flags);
        ASSERT_EQ(0u, leaf.contributors[0].source.allele_index);
        if (mismatch) {
            ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH, leaf.projection_status);
            ASSERT(leaf.cds == NULL); ASSERT(leaf.protein == NULL);
        } else {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
            ASSERT_EQ(12u, leaf.cds_length); ASSERT_MEM_EQ(f.reference, leaf.cds, 12u);
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    }
    PASS();
}

TEST haplotype_stream_retained_reference_replaces_but_omitted_missing_does_not(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        for (unsigned missing = 0u; missing < 2u; missing++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.strands[0] = reverse ? -1 : 1;
            duckvep_haplotype_stream_t *s = &f.stream;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                s, &f.model, &f.exons, &f.sequences, &f.buffers));
            const uint8_t *ref = (const uint8_t *)(reverse ? "TTT" : "AAA");
            const uint8_t *alt = (const uint8_t *)(reverse ? "G" : "C");
            duckvep_haplotype_source_t sources[] = {
                {1u, ref, ref, 100u, 0u, 3u, 3u, 0u, 1u, 0u},
                {2u, ref, alt, 101u, 0u, 1u, 1u, 1u, 1u, 0u}
            };
            for (size_t i = 0u; i < 2u; i++) {
                duckvep_raw_gt_t call;
                const char *gt = i ? "1|1" : missing ? ".|." : "0|1";
                ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                    (const uint8_t *)gt, strlen(gt), 1u, &call));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &sources[i]));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
            duckvep_haplotype_leaf_t leaf;
            unsigned checked = 0u;
            while (duckvep_haplotype_stream_next(s, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                const duckvep_carrier_call_t *carrier = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
                ASSERT(carrier != NULL);
                if (!missing && carrier->key.lane == 2u) continue;
                ASSERT_EQ(12u, leaf.cds_length);
                ASSERT_EQ(0, leaf.nominal_length_diff);
                ASSERT_EQ(2u, leaf.contributor_count);
                ASSERT_EQ(missing ? 1u : 2u, leaf.edit_count);
                ASSERT_EQ(missing ? 0u : 1u, leaf.ordered_replacements);
                ASSERT_EQ(missing ? DUCKVEP_HAPLOTYPE_CONDITIONAL : DUCKVEP_HAPLOTYPE_OK,
                    leaf.sequence_status);
                uint8_t expected[12]; memcpy(expected, f.reference, 12u);
                if (missing) expected[reverse ? 10u : 1u] = 'C';
                ASSERT_MEM_EQ(expected, leaf.cds, 12u);
                ASSERT_EQ(0u, leaf.evidence_flags & DUCKVEP_CARRIER_REFERENCE_REPLAY);
                if (!missing) {
                    ASSERT_EQ(1u, leaf.block_count);
                    ASSERT_EQ(2u, leaf.blocks[0].edit_count);
                    ASSERT_EQ(0, leaf.blocks[0].length_diff);
                    ASSERT_EQ(1u, leaf.contributors[0].source_replaced);
                }
                checked++;
            }
            ASSERT_EQ(1u, checked);
        }
    }
    PASS();
}

TEST haplotype_stream_nominal_lengths_survive_clipped_and_disjoint_replay(void) {
    static const struct {
        uint32_t start[2];
        uint16_t ref_len[2];
        const char *alt[2];
        int64_t nominal;
        size_t length, edits, sources;
        uint32_t flags;
    } cases[] = {
        /* The undefined slot deletes the last base; the following substitution
         * inserts at that emptied position while retaining nominal REF length. */
        {{12u, 12u}, {1u, 1u}, {"", "C"}, -1, 12u, 2u, 2u, 3u},
        /* The six-base source sees only three bases after the other deletion.
         * Its ALT already matches, but its nominal -3 remains in the sum. */
        {{7u, 9u}, {6u, 4u}, {"AAA", "A"}, -6, 9u, 1u, 2u, 1u},
        {{9u, 0u}, {4u, 0u}, {"A", ""}, -3, 9u, 1u, 1u, 1u},
        {{4u, 9u}, {1u, 2u}, {"AAA", "A"}, 1, 13u, 2u, 2u, 3u}
    };
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        for (size_t scenario = 0u; scenario < sizeof(cases) / sizeof(cases[0]); scenario++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.strands[0] = reverse ? -1 : 1;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            uint8_t refs[2][6], alts[2][3];
            duckvep_haplotype_source_t sources[2];
            for (size_t i = 0u; i < cases[scenario].sources; i++) {
                size_t alt_length = strlen(cases[scenario].alt[i]);
                memset(refs[i], reverse ? 'T' : 'A', cases[scenario].ref_len[i]);
                for (size_t j = 0u; j < alt_length; j++) alts[i][j] = (uint8_t)
                    haplo_test_variant_from_tx_base(cases[scenario].alt[i][reverse ? alt_length - 1u - j : j],
                        (int)reverse);
                sources[i] = (duckvep_haplotype_source_t){0u, refs[i], alts[i],
                    reverse ? 113u - cases[scenario].start[i] - cases[scenario].ref_len[i]
                            : 99u + cases[scenario].start[i],
                    0u, cases[scenario].ref_len[i], (uint16_t)alt_length,
                    scenario == 0u && i == 0u ? UINT32_MAX : 1u, 1u, 0u};
            }
            if (cases[scenario].sources == 2u && sources[0].pos1 > sources[1].pos1) {
                duckvep_haplotype_source_t tmp = sources[0]; sources[0] = sources[1]; sources[1] = tmp;
            }
            for (size_t i = 0u; i < cases[scenario].sources; i++) {
                sources[i].event_id = i + 1u;
                const char *gt = scenario ? "1|1" : sources[i].allele_index == UINT32_MAX ? ".|1" : "0|1";
                duckvep_raw_gt_t call;
                ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                    (const uint8_t *)gt, strlen(gt), 1u, &call));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &sources[i]));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                    &f.stream, 0u, 0u, &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
            ASSERT_EQ(scenario ? DUCKVEP_HAPLOTYPE_OK : DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
            ASSERT_EQ(cases[scenario].nominal, leaf.nominal_length_diff);
            ASSERT_EQ(cases[scenario].length, leaf.cds_length);
            ASSERT_EQ(cases[scenario].edits, leaf.edit_count);
            ASSERT_EQ(cases[scenario].sources, leaf.contributor_count);
            ASSERT_EQ(cases[scenario].flags, leaf.flags);
            ASSERT_EQ(scenario < 2u, leaf.ordered_replacements);
            ASSERT_EQ(scenario ? 2u : 1u, leaf.carriers.call_count);
            for (size_t i = 0u; i < leaf.cds_length; i++)
                ASSERT_EQ(!scenario && i == 11u ? 'C' : 'A', leaf.cds[i]);
            if (scenario < 2u) ASSERT(leaf.nominal_length_diff != (int64_t)leaf.cds_length - 12);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
        }
    }
    PASS();
}

TEST haplotype_stream_raw_record_identity_and_evidence_are_checked(void) {
    for (unsigned scenario = 0u; scenario < 6u; scenario++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)"A", (const uint8_t *)"C",
            101u, 0u, 1u, 1u, 1u, 1u, 0u};
        duckvep_raw_gt_t call;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw((const uint8_t *)"1", 1u, 1u, &call));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
        duckvep_haplotype_stream_status_t expected = DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        if (scenario == 0u) {
            source.allele_index = 2u; source.ref = (const uint8_t *)"G";
            ASSERT_EQ(expected, duckvep_haplotype_stream_begin(s, &source));
        } else if (scenario == 1u) {
            expected = DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER;
            ASSERT_EQ(expected, duckvep_haplotype_stream_begin(s, &source));
        } else if (scenario == 2u) {
            source.allele_index = UINT32_MAX; source.alt = (const uint8_t *)""; source.alt_len = 0u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
            duckvep_carrier_key_t key = {0u, 0, 1u, 2u, 0u};
            ASSERT_EQ(expected, haplotype_test_push_called(s, &key));
        } else if (scenario == 3u) {
            call.allele_index[1] = 0u; /* A retained single slot has no second allele. */
            ASSERT_EQ(expected, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
        } else if (scenario == 4u) {
            int32_t gt = 1;
            duckvep_haplotype_call_t decoded = {&gt, NULL, 0u, 1u, 1u, {0, 0u}, DUCKVEP_PHASE_STRICT};
            ASSERT_EQ(expected, duckvep_haplotype_stream_push_call(s, 0u, &decoded, NULL, 0u));
        } else {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
            expected = DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR;
            ASSERT_EQ(expected, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
            ASSERT_EQ(DUCKVEP_CARRIERS_DUPLICATE_CALL, s->carrier_error);
        }
        ASSERT_EQ(expected, duckvep_haplotype_stream_finish(s));
    }
    PASS();
}

TEST haplotype_stream_recycles_owned_alleles_and_projection_slots(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 64u);
    /* Two live transcripts suffice for all 64 inputs. Unequal allele sizes
     * force wrap padding while older events and projections are still live. */
    f.buffers.event_capacity = f.buffers.projection_capacity = 2u;
    f.buffers.allele_capacity = 8u;
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    uint32_t drained = 0u;
    for (uint32_t i = 0u; i <= 64u; i++) {
        uint8_t ref = 'A', alt[3] = {'A', 'A', 'A'};
        if (i % 2u) alt[0] = 'C';
        duckvep_haplotype_source_t source = {UINT64_MAX - i, &ref, alt,
            i < 64u ? f.starts[i] : 0u, 0u, 1u, (uint16_t)(i % 2u ? 1u : 3u), 0u, 0u, 0u};
        duckvep_haplotype_stream_status_t status;
        while ((status = i == 64u ? duckvep_haplotype_stream_finish(s)
            : haplotype_test_begin_candidates(s, &source, &i, 1u)) ==
                DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) {
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
            ASSERT_EQ(drained, leaf.carriers.transcript_index);
            ASSERT_EQ(1u, leaf.contributor_count);
            ASSERT_EQ(UINT64_MAX - drained, leaf.contributors[0].source.event_id);
            ASSERT_EQ('A', leaf.contributors[0].source.ref[0]);
            ASSERT_EQ(drained % 2u ? 'C' : 'A', leaf.contributors[0].source.alt[0]);
            const duckvep_haplotype_contributor_t *contributor = &leaf.contributors[0];
            const duckvep_haplotype_source_t *retained = &contributor->source;
            ASSERT(retained->ref >= f.alleles);
            ASSERT(retained->alt == retained->ref + retained->ref_len);
            ASSERT(retained->alt + retained->alt_len <= f.alleles + f.buffers.allele_capacity);
            duckvep_event_t prepared = {0};
            ASSERT(duckvep_event_prepare_small(retained->pos1, retained->ref, retained->ref_len,
                retained->alt, retained->alt_len, &prepared));
            prepared.chrom_id = retained->chrom_id;
            ASSERT(contributor->prepared != NULL);
            ASSERT_MEM_EQ(&prepared, contributor->prepared, sizeof prepared);
            ASSERT(contributor->projected != NULL);
            ASSERT_EQ(drained % 2u ? 1u : 2u, contributor->projected->cds_start);
            ASSERT_EQ(drained % 2u ? 1u : 0u, contributor->projected->ref_len);
            ASSERT_EQ(drained % 2u ? 1u : 2u, contributor->projected->alt_len);
            ASSERT(contributor->projected->ref == (drained % 2u ? retained->ref : NULL));
            ASSERT(contributor->projected->alt == retained->alt + prepared.alt_diff_offset);
            ASSERT_EQ(1, contributor->projected->variant_strand);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
            ASSERT_EQ(drained % 2u ? 12u : 14u, leaf.cds_length);
            uint8_t expected[14]; memset(expected, 'A', sizeof(expected));
            if (drained % 2u) expected[0] = 'C';
            ASSERT_EQ(0, memcmp(expected, leaf.cds, leaf.cds_length));
            ASSERT_EQ(4u, leaf.protein_length);
            ASSERT_EQ(0, memcmp(drained % 2u ? "QKKK" : "KKKK", leaf.protein, 4u));
            ASSERT_EQ(3u, leaf.carriers.call_count);
            uint32_t seen = 0u;
            for (uint32_t id = leaf.carriers.first_call; id;) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, id);
                ASSERT(call && call->key.sample_index < 3u);
                uint32_t sample = call->key.sample_index;
                ASSERT_EQ(0u, seen & (1u << sample)); seen |= 1u << sample;
                ASSERT_EQ(sample + 1u, call->key.lane);
                ASSERT_EQ(2u * sample + 1u, call->key.ploidy);
                ASSERT_EQ(sample % 2u, call->key.phase_set_present);
                ASSERT_EQ(sample % 2u ? -17 : 0, call->key.phase_set);
                id = call->next_leaf;
            }
            ASSERT_EQ(7u, seen);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
            drained++;
        }
        ASSERT_EQ(i == 64u ? DUCKVEP_HAPLOTYPE_STREAM_DONE : DUCKVEP_HAPLOTYPE_STREAM_OK, status);
        if (i == 64u) break;
        /* Destroy the input immediately, before any carrier is appended. */
        ref = 'N'; memset(alt, 'N', sizeof(alt));
        for (uint32_t sample = 0u; sample < 3u; sample++) {
            duckvep_carrier_key_t key = {sample, -17, (uint16_t)(sample + 1u),
                (uint16_t)(2u * sample + 1u), (uint8_t)(sample % 2u)};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
    }
    ASSERT_EQ(64u, drained); ASSERT_EQ(64u, s->input_events);
    ASSERT_EQ(64u, s->projected_events); ASSERT_EQ(64u, s->completed_leaves);
    ASSERT_EQ(64u * 13u, s->translated_bases);
    ASSERT_EQ(2u, s->peak_events); ASSERT_EQ(2u, s->peak_projections);
    ASSERT_EQ(8u, s->peak_alleles);
    ASSERT_EQ(0u, s->event_count); ASSERT_EQ(0u, s->projection_count);
    ASSERT_EQ(0u, s->allele_count); ASSERT_EQ(0u, s->carriers.call_count);
    PASS();
}

TEST haplotype_stream_projects_once_and_preserves_occupied_ancestors(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 2u);
    f.starts[1] = f.starts[0]; f.ends[1] = f.ends[0];
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    uint32_t candidates[] = {0u, 1u};
    /* An earlier upload can project downstream of the next uploaded allele
     * after unchanged REF/ALT prefixes are trimmed. Preserve the raw alleles. */
    duckvep_haplotype_source_t source = {200u, (const uint8_t *)"AAAAAA", (const uint8_t *)"AAAAAC",
        100u, 0u, 6u, 6u, 0u, 0u, 0u};
    duckvep_carrier_key_t keys[] = {{0u, 0, 1u, 2u, 0u}, {42u, 0, 2u, 4u, 1u}};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, candidates, 2u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &keys[0]));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &keys[1]));
    source.event_id++; source.pos1 = 102u; source.alt = (const uint8_t *)"G";
    source.ref = (const uint8_t *)"A"; source.ref_len = source.alt_len = 1u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, candidates, 2u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &keys[0]));
    uint32_t paths = 0u;
    duckvep_haplotype_stream_status_t status;
    while ((status = duckvep_haplotype_stream_finish(s)) == DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) {
        duckvep_haplotype_leaf_t leaf;
        while ((status = duckvep_haplotype_stream_next(s, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
            ASSERT_EQ(1u, leaf.carriers.call_count);
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
            ASSERT(call); ASSERT(leaf.carriers.transcript_index < 2u);
            uint32_t bit = leaf.carriers.transcript_index * 2u + (call->key.sample_index == 42u);
            ASSERT_EQ(0u, paths & (1u << bit)); paths |= 1u << bit;
            size_t count = call->key.sample_index == 42u ? 1u : 2u;
            ASSERT_EQ(count, leaf.contributor_count);
            ASSERT_EQ(200u, leaf.contributors[0].source.event_id);
            ASSERT_EQ(6u, leaf.contributors[0].source.alt_len);
            ASSERT_EQ(0, memcmp("AAAAAC", leaf.contributors[0].source.alt, 6u));
            if (count == 2u) ASSERT_EQ(201u, leaf.contributors[1].source.event_id);
            uint8_t expected[12]; memset(expected, 'A', sizeof(expected));
            expected[5] = 'C'; if (count == 2u) expected[2] = 'G';
            ASSERT_EQ(12u, leaf.cds_length);
            ASSERT_EQ(0, memcmp(expected, leaf.cds, sizeof(expected)));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
    ASSERT_EQ(15u, paths); ASSERT_EQ(4u, s->projected_events);
    ASSERT_EQ(4u, s->completed_leaves); ASSERT_EQ(48u, s->translated_bases);
    PASS();
}

TEST haplotype_stream_keeps_noncoding_contributors_without_poisoning_cds(void) {
    /* Two exons with three UTR bases at each transcript end. Every small
     * genomic span is classified independently by its literal CDS-base set;
     * strict replay withholds crossing spans, while raw compatibility retains
     * their omission as conditional evidence and replays the mapped source. */
    for (unsigned raw = 0u; raw < 2u; raw++) {
        for (int strand = -1; strand <= 1; strand += 2) {
            for (uint32_t pos = 100u; pos <= 151u; pos++) {
                for (uint16_t n = 1u; n <= 52u && pos + n - 1u <= 151u; n++) {
                    if (pos <= 105u && pos + n - 1u >= 105u) continue;
                    struct haplotype_stream_scene f;
                    haplotype_stream_scene_prepare(&f, 1u);
                    uint32_t cds_start = 103u, cds_end = 148u;
                    uint32_t exon_start[] = {100u, 143u}, exon_end[] = {108u, 151u};
                    uint32_t cdna_start[] = {1u, 10u}, cdna_end[] = {9u, 18u};
                    if (strand < 0) {
                        exon_start[0] = 143u; exon_start[1] = 100u;
                        exon_end[0] = 151u; exon_end[1] = 108u;
                    }
                    f.ends[0] = 151u; f.strands[0] = (int8_t)strand;
                    f.model.cds_start1 = &cds_start; f.model.cds_end1 = &cds_end;
                    f.exon_count[0] = 2u;
                    f.exons = (duckvep_exon_model_t){.exon_count = 2u,
                        .start1 = exon_start, .end1 = exon_end,
                        .cdna_start1 = cdna_start, .cdna_end1 = cdna_end};
                    memset(f.reference, strand > 0 ? 'A' : 'T', sizeof(f.reference));
                    duckvep_haplotype_stream_t *s = &f.stream;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                        s, &f.model, &f.exons, &f.sequences, &f.buffers));
                    uint8_t reference[52], alternate[52];
                    memset(reference, 'A', sizeof(reference));
                    memset(alternate, 'C', sizeof(alternate));
                    duckvep_haplotype_source_t events[] = {
                        {1u, reference, alternate, pos, 0u, n, n, raw, (uint8_t)raw, 0u},
                        {2u, (const uint8_t *)"A", (const uint8_t *)"G", 105u, 0u, 1u, 1u,
                            raw, (uint8_t)raw, 0u}
                    };
                    if (pos > 105u) {
                        duckvep_haplotype_source_t swap = events[0];
                        events[0] = events[1]; events[1] = swap;
                    }
                    duckvep_carrier_key_t key = {0u, 0, 1u, 1u, 0u};
                    uint32_t tx = 0u;
                    for (size_t i = 0u; i < 2u; i++) {
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                            haplotype_test_begin_candidates(s, events + i, &tx, 1u));
                        if (raw) {
                            duckvep_raw_gt_t call;
                            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                                (const uint8_t *)"1|1", 3u, 1u, &call));
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
                        } else {
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                haplotype_test_push_called(s, &key));
                        }
                    }
                    uint32_t coding_bases = 0u;
                    uint8_t expected[12];
                    memcpy(expected, f.reference, sizeof(expected));
                    expected[strand > 0 ? 2u : 9u] = strand > 0 ? 'G' : 'C';
                    for (uint32_t p = pos; p < pos + n; p++) {
                        if ((p >= 103u && p <= 108u) || (p >= 143u && p <= 148u)) {
                            uint32_t offset = p <= 108u ? p - 103u : p - 143u + 6u;
                            expected[strand > 0 ? offset : 11u - offset] = strand > 0 ? 'C' : 'G';
                            coding_bases++;
                        }
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
                        duckvep_haplotype_stream_finish(s));
                    duckvep_haplotype_leaf_t leaf;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
                    ASSERT_EQ(2u, leaf.contributor_count);
                    ASSERT_EQ(events[0].event_id, leaf.contributors[0].source.event_id);
                    ASSERT_EQ(events[1].event_id, leaf.contributors[1].source.event_id);
                    int omitted = raw && coding_bases && coding_bases != n;
                    ASSERT_EQ(coding_bases == n ? DUCKVEP_CDS_EDIT_OK : omitted
                        ? DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED : DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                        leaf.contributors[pos < 105u ? 0u : 1u].projection_status);
                    if (coding_bases == 0u || coding_bases == n || omitted) {
                        if (omitted) {
                            memcpy(expected, f.reference, sizeof(expected));
                            expected[strand > 0 ? 2u : 9u] = strand > 0 ? 'G' : 'C';
                            ASSERT(leaf.contributors[pos < 105u ? 0u : 1u].evidence_flags &
                                DUCKVEP_CARRIER_CONDITIONAL);
                        }
                        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                        ASSERT_EQ(omitted ? DUCKVEP_HAPLOTYPE_CONDITIONAL : DUCKVEP_HAPLOTYPE_OK,
                            leaf.sequence_status);
                        ASSERT_EQ(12u, leaf.cds_length);
                        ASSERT_EQ(0, memcmp(expected, leaf.cds, sizeof(expected)));
                        ASSERT_EQ(coding_bases && !omitted ? 2u : 1u, leaf.edit_count);
                    } else {
                        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS, leaf.projection_status);
                        ASSERT(leaf.cds == NULL && leaf.protein == NULL && leaf.blocks == NULL);
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
                }
            }
        }
    }
    PASS();
}

TEST haplotype_source_omission_checks_layout_sequence_and_ref(void) {
    for (int strand = -1; strand <= 1; strand += 2) {
        for (unsigned scenario = 0u; scenario < 13u; scenario++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            uint32_t starts[] = {100u, 146u}, ends[] = {105u, 151u};
            uint32_t cs[] = {1u, 7u}, ce[] = {6u, 12u};
            if (strand < 0) { starts[0] = 146u; starts[1] = 100u; ends[0] = 151u; ends[1] = 105u; }
            f.ends[0] = 151u; f.strands[0] = (int8_t)strand; f.exon_count[0] = 2u;
            f.exons = (duckvep_exon_model_t){.exon_count = 2u,
                .start1 = starts, .end1 = ends, .cdna_start1 = cs, .cdna_end1 = ce};
            memset(f.reference, strand > 0 ? 'A' : 'T', sizeof(f.reference));
            uint8_t ref[44]; memset(ref, 'A', sizeof(ref));
            uint8_t alt = 'C';
            duckvep_cds_edit_status_t expected = DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED;
            if (scenario == 1u || scenario == 2u) {
                ref[scenario == 1u ? 0u : 43u] = 'C'; expected = DUCKVEP_CDS_EDIT_REF_MISMATCH;
            }
            if (scenario == 3u) { alt = 'Z'; expected = DUCKVEP_CDS_EDIT_INVALID_ALLELE; }
            if (scenario == 4u) { f.sequence_offsets[0] = 13u; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 5u) { f.lengths[0] = 11u; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 6u) { ce[0]--; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 7u) { ends[0]++; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 8u) { f.exons.cdna_start1 = NULL; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            uint32_t cached_start = 2u, cached_end = 13u, cached_exon = 0u;
            uint8_t cached_phase = 0u;
            if (scenario == 9u) {
                f.model.cds_cdna_start1 = &cached_start; f.model.cds_cdna_end1 = &cached_end;
                f.model.cds_start_exon_index = &cached_exon; f.model.cds_phase_offset = &cached_phase;
                expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
            }
            if (scenario == 10u) { f.lengths[0] = 0u; expected = DUCKVEP_CDS_EDIT_OUT_OF_CDS; }
            if (scenario == 11u) ref[20] = 'C'; /* Intronic REF is not represented in this CDS pool. */
            if (scenario == 12u) { ref[20] = 'N'; expected = DUCKVEP_CDS_EDIT_INVALID_ALLELE; }
            duckvep_event_t event;
            ASSERT(duckvep_event_prepare_replacement(104u, ref, sizeof(ref), &alt, 1u, &event));
            duckvep_prepared_cds_allele_t allele = {&event, ref, &alt, ref, sizeof(ref), 1u, 1};
            duckvep_haplotype_edit_t edit;
            duckvep_cds_edit_status_t status = duckvep_compat_vep116_source_cds_edit_build(
                &f.model, &f.exons, &f.sequences, 0u, (int8_t)strand, &allele, &edit);
            if (status != expected) fprintf(stderr, "source mapping strand=%d scenario=%u\n", strand, scenario);
            ASSERT_EQ_FMT(expected, status, "%d");
            if (expected == DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED) {
                duckvep_haplotype_edit_t zero = {0};
                ASSERT_MEM_EQ(&zero, &edit, sizeof(edit));
            }
            if (scenario == 5u) {
                /* Intronic context cannot conceal an inconsistent CDS extent. */
                ASSERT(duckvep_event_prepare_replacement(120u, ref, 1u, &alt, 1u, &event));
                allele.ref_length = allele.alt_length = 1u;
                ASSERT_EQ(DUCKVEP_CDS_EDIT_INVALID_ARG, duckvep_compat_vep116_source_cds_edit_build(
                    &f.model, &f.exons, &f.sequences, 0u, (int8_t)strand, &allele, &edit));
            }
        }
    }
    PASS();
}

TEST haplotype_stream_unmapped_sources_retain_raw_slot_evidence(void) {
    const char *spellings[] = {"1|1", "0|1", ".|1", ".", "0|0"};
    for (int strand = -1; strand <= 1; strand += 2) {
        for (unsigned mode = 0u; mode < 5u; mode++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            uint32_t starts[] = {100u, 146u}, ends[] = {105u, 151u};
            uint32_t cs[] = {1u, 7u}, ce[] = {6u, 12u};
            if (strand < 0) { starts[0] = 146u; starts[1] = 100u; ends[0] = 151u; ends[1] = 105u; }
            f.ends[0] = 151u; f.strands[0] = (int8_t)strand; f.exon_count[0] = 2u;
            f.exons = (duckvep_exon_model_t){.exon_count = 2u,
                .start1 = starts, .end1 = ends, .cdna_start1 = cs, .cdna_end1 = ce};
            memset(f.reference, strand > 0 ? 'A' : 'T', sizeof(f.reference));
            /* Three interpretations of one record occupy 44+44+44+1+44 bytes. */
            uint8_t allele_storage[256];
            f.buffers.alleles = allele_storage; f.buffers.allele_capacity = sizeof(allele_storage);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            uint8_t ref[44]; memset(ref, 'A', sizeof(ref));
            duckvep_raw_gt_t call, anchor;
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                (const uint8_t *)spellings[mode], strlen(spellings[mode]), 1u, &call));
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                (const uint8_t *)"1|1", 3u, 1u, &anchor));
            for (unsigned interpretation = 0u; interpretation < 4u; interpretation++) {
                duckvep_haplotype_source_t source = {1u, ref,
                    interpretation == 0u ? ref : (const uint8_t *)(interpretation == 1u ? "C" : ""),
                    104u, 0u, sizeof(ref), interpretation == 0u ? sizeof(ref) : interpretation == 1u ? 1u : 0u,
                    interpretation < 2u ? interpretation : UINT32_MAX, 1u, 0u};
                if (interpretation == 3u) source = (duckvep_haplotype_source_t){
                    2u, (const uint8_t *)"A", (const uint8_t *)"G", 150u, 0u, 1u, 1u, 1u, 1u, 0u};
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                    &f.stream, 0u, 0u, interpretation == 3u ? &anchor : &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            unsigned carriers = 0u;
            duckvep_haplotype_stream_status_t status;
            while ((status = duckvep_haplotype_stream_next(&f.stream, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                ASSERT_EQ(mode == 4u ? DUCKVEP_HAPLOTYPE_OK : DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
                ASSERT_EQ(12u, leaf.cds_length);
                const char *expected_cds = strand > 0 ? "AAAAAAAAAAGA" : "TCTTTTTTTTTT";
                ASSERT_MEM_EQ(expected_cds, leaf.cds, 12u);
                ASSERT_EQ(1u, leaf.edit_count); ASSERT_EQ(2u, leaf.edit_event_ids[0]);
                ASSERT_EQ(mode == 4u ? 1u : 2u, leaf.contributor_count);
                if (mode != 4u) {
                    ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
                    ASSERT_EQ(DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED, leaf.contributors[0].projection_status);
                    ASSERT(leaf.contributors[0].evidence_flags & DUCKVEP_CARRIER_CONDITIONAL);
                }
                carriers += leaf.carriers.call_count;
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status); ASSERT_EQ(2u, carriers);
        }
    }
    PASS();
}

TEST haplotype_stream_reference_only_protein_uses_call_retention(void) {
    const char *references[] = {"ATGAAACCCTAA", "ATGAAATAACCC", "CTGAAACCCTAA", "ATGGCNCCCTAA"};
    const char *curated[] = {"MKP*", "MK*P", "MKP*", "MAP*"};
    const char *raw[] = {"MKP*", "MK*", "LKP*", "MAP*"};
    const char *alternate[] = {"MQP*", "MQ*", "LQP*", "MPP*"};
    const char *coding[] = {"MKP*", "MK*P", "LKP*", "MXP*"};
    for (unsigned reference = 0u; reference < 4u; reference++) {
        for (int strand = -1; strand <= 1; strand += 2) {
            for (unsigned peptide_edit = 0u; peptide_edit < 2u; peptide_edit++) {
                /* Missing-only, coding, shadowed, unmapped, UTR and intronic calls. */
                for (unsigned mode = 0u; mode < 7u; mode++) {
                    struct haplotype_stream_scene f;
                    haplotype_stream_scene_prepare(&f, 1u);
                    memcpy(f.reference, references[reference], 12u);
                    f.ends[0] = 117u; f.cdna_ends[0] = 18u; f.strands[0] = (int8_t)strand;
                    uint32_t cds_start = 103u, cds_end = 114u;
                    f.model.cds_start1 = &cds_start; f.model.cds_end1 = &cds_end;
                    uint32_t exon_starts[] = {100u, mode == 5u ? 121u : 133u};
                    uint32_t exon_ends[] = {108u, exon_starts[1] + 8u};
                    uint32_t cdna_starts[] = {1u, 10u}, cdna_ends[] = {9u, 18u};
                    if (mode >= 5u) {
                        f.ends[0] = exon_ends[1]; cds_end = exon_starts[1] + 5u;
                        if (strand < 0) {
                            exon_starts[0] = exon_starts[1]; exon_starts[1] = 100u;
                            exon_ends[0] = exon_ends[1]; exon_ends[1] = 108u;
                        }
                        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
                        f.exons.start1 = exon_starts; f.exons.end1 = exon_ends;
                        f.exons.cdna_start1 = cdna_starts; f.exons.cdna_end1 = cdna_ends;
                    }
                    uint32_t peptide_offsets[] = {0u, peptide_edit}, peptide_position = 2u;
                    uint8_t peptide_alt = 'W';
                    f.sequences.peptide_edit_offset = peptide_offsets;
                    f.sequences.peptide_edit_position1 = &peptide_position;
                    f.sequences.peptide_edit_alt = &peptide_alt;
                    f.sequences.peptide_edit_count = peptide_edit;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
                    duckvep_raw_gt_t missing, call;
                    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                        (const uint8_t *)".", 1u, 1u, &missing));
                    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep116_raw(
                        (const uint8_t *)(mode ? "0|1" : "."), mode ? 3u : 1u, 1u, &call));
                    duckvep_haplotype_source_t source = {1u, (const uint8_t *)"C",
                        (const uint8_t *)"C", 100u, 0u, 1u, 1u, 0u, 1u, 0u};
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                        duckvep_haplotype_stream_push_raw_call(&f.stream, 0u, 0u, &missing, 1u));
                    uint8_t ref[3], alt[3];
                    uint32_t pos = mode >= 5u ? 115u : mode == 3u ? 102u
                        : mode == 4u ? 101u : strand > 0 ? 106u : 111u;
                    uint16_t length = mode == 3u ? 3u : 1u;
                    for (uint16_t i = 0u; i < length; i++) {
                        uint32_t p = pos + i;
                        ref[i] = mode >= 5u || p < cds_start ? 'C' : strand > 0 ? f.reference[p - cds_start]
                            : kprop_hgvs_complement(f.reference[cds_end - p]);
                        alt[i] = mode == 4u ? 'G' : strand > 0 ? 'C' : 'G';
                    }
                    for (uint32_t allele = 0u; allele < 2u; allele++) {
                        source = (duckvep_haplotype_source_t){2u, ref, allele ? alt : ref,
                            pos, 0u, length, length, allele, 1u, 0u};
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                            &f.stream, 0u, 0u, &call, mode != 2u));
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
                        duckvep_haplotype_stream_finish(&f.stream));
                    duckvep_haplotype_leaf_t leaf;
                    unsigned seen = 0u;
                    duckvep_haplotype_stream_status_t status;
                    char expected_reference[5]; memcpy(expected_reference, curated[reference], 5u);
                    if (peptide_edit) expected_reference[1] = 'W';
                    while ((status = duckvep_haplotype_stream_next(&f.stream, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
                        ASSERT_EQ(4u, leaf.reference_protein_length);
                        ASSERT_MEM_EQ(expected_reference, leaf.reference_protein, 4u);
                        ASSERT_MEM_EQ(coding[reference], leaf.reference_coding_protein, 5u);
                        ASSERT_EQ(4u, leaf.reference_coding_translation.length);
                        ASSERT_EQ(reference == 1u ? 3u : 4u,
                            leaf.reference_coding_translation.first_stop_position1);
                        ASSERT_EQ(reference != 3u, leaf.reference_coding_translation.unambiguous);
                        for (uint32_t id = leaf.carriers.first_call; id;) {
                            const duckvep_carrier_call_t *carrier = duckvep_carriers_call(&f.stream.carriers, id);
                            ASSERT(carrier); ASSERT(carrier->key.lane == 1u || carrier->key.lane == 2u);
                            unsigned bit = 1u << (carrier->key.lane - 1u);
                            ASSERT_EQ(0u, seen & bit); seen |= bit;
                            int edited = mode == 1u && carrier->key.lane == 2u;
                            int reference_route = !mode || mode >= 5u;
                            const char *expected = reference_route ? expected_reference
                                : edited ? alternate[reference] : raw[reference];
                            ASSERT_EQ(strlen(expected), leaf.protein_length);
                            ASSERT_MEM_EQ(expected, leaf.protein, leaf.protein_length);
                            ASSERT_EQ((size_t)edited, leaf.edit_count);
                            ASSERT_EQ(reference_route, leaf.protein == leaf.reference_protein);
                            if (reference_route) ASSERT_EQ(0u, leaf.flags & DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED);
                            id = carrier->next_leaf;
                        }
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status); ASSERT_EQ(3u, seen);
                }
            }
        }
    }
    PASS();
}

TEST haplotype_stream_retains_conflicts_and_latches_resource_errors(void) {
    struct haplotype_stream_scene f;
    duckvep_carrier_key_t key = {0u, 1, 1u, 2u, 1u};
    uint32_t tx = 0u;
    /* One bad REF, two conflicting edits, and each named capacity limit. */
    for (unsigned scenario = 0u; scenario < 13u; scenario++) {
        haplotype_stream_scene_prepare(&f, 1u);
        uint32_t peptide_offsets[] = {1u, 0u};
        if (scenario == 2u) f.buffers.event_capacity = 1u;
        if (scenario == 3u) f.buffers.projection_capacity = 1u;
        if (scenario == 4u) f.buffers.allele_capacity = 3u;
        if (scenario == 5u) f.buffers.leaf_capacity = 1u;
        if (scenario == 6u) f.buffers.cds_capacity = 11u;
        if (scenario == 7u) f.buffers.protein_capacity = 4u;
        if (scenario == 9u) f.buffers.edit_capacity = 1u;
        if (scenario == 10u) f.buffers.reference_protein_capacity = 5u;
        if (scenario == 11u) f.sequences.peptide_edit_offset = peptide_offsets;
        if (scenario == 12u) f.sequences.peptide_edit_count = 1u;
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)(scenario ? "A" : "C"),
            (const uint8_t *)"G", 100u, 0u, 1u, 1u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        duckvep_haplotype_stream_status_t status = DUCKVEP_HAPLOTYPE_STREAM_OK;
        if ((scenario >= 1u && scenario <= 5u) || scenario == 9u) {
            source.event_id = 2u; source.alt = (const uint8_t *)"C";
            if (scenario != 1u) source.pos1++;
            status = haplotype_test_begin_candidates(s, &source, &tx, 1u);
            if (scenario == 2u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL, status);
            else if (scenario == 3u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL, status);
            else if (scenario == 4u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL, status);
            else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, status);
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
            }
        } else if (scenario == 8u) {
            status = haplotype_test_begin_candidates(s, &source, &tx, 1u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER, status);
        }
        duckvep_haplotype_leaf_t leaf;
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
            status = duckvep_haplotype_stream_next(s, &leaf);
            if (scenario == 5u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL, status);
            else if (scenario == 9u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL, status);
            else if (scenario == 6u || scenario == 7u || scenario == 10u)
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, status);
            else if (scenario >= 11u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG, status);
            else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, status);
                ASSERT_EQ(scenario == 0u ? DUCKVEP_CDS_EDIT_REF_MISMATCH : DUCKVEP_CDS_EDIT_OK,
                          leaf.projection_status);
                if (scenario == 1u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER, leaf.sequence_status);
                ASSERT_EQ(scenario == 0u ? 1u : 2u, leaf.contributor_count);
                ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
                if (scenario == 1u) ASSERT_EQ(2u, leaf.contributors[1].source.event_id);
                ASSERT(!leaf.cds && !leaf.protein);
                ASSERT_EQ(0u, leaf.cds_length); ASSERT_EQ(0u, leaf.protein_length);
                ASSERT_EQ(0, leaf.nominal_length_diff);
                ASSERT_EQ(1u, leaf.carriers.call_count);
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
            }
        }
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) {
            ASSERT_EQ(status, haplotype_test_push_called(s, &key));
            ASSERT_EQ(status, duckvep_haplotype_stream_finish(s));
            memset(&leaf, 0xa5, sizeof(leaf));
            ASSERT_EQ(status, duckvep_haplotype_stream_next(s, &leaf));
            duckvep_haplotype_leaf_t zero_leaf = {0};
            ASSERT_EQ(0, memcmp(&leaf, &zero_leaf, sizeof(leaf)));
        }
    }
    PASS();
}

TEST haplotype_stream_mnv_context_does_not_conflict_with_another_edit(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        if (reverse) {
            memset(f.reference, 'T', sizeof(f.reference));
            f.strands[0] = -1;
        }
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_carrier_key_t key = {0u, 10, 1u, 2u, 1u};
        uint32_t tx = 0u;
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)"AAA", (const uint8_t *)"CAC",
            100u, 0u, 3u, 3u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        source.event_id = 2u; source.pos1 = 101u; source.ref_len = source.alt_len = 1u;
        source.ref = (const uint8_t *)"A"; source.alt = (const uint8_t *)"G";
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
        ASSERT_EQ(2u, leaf.contributor_count); ASSERT_EQ(3u, leaf.edit_count);
        ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
        ASSERT_EQ(0, memcmp("CAC", leaf.contributors[0].source.alt, 3u));
        ASSERT_EQ(2u, leaf.contributors[1].source.event_id);
        ASSERT_EQ(12u, leaf.cds_length); ASSERT_EQ(4u, leaf.protein_length);
        ASSERT_EQ(0, memcmp(reverse ? "TTTTTTTTTGCG" : "CGCAAAAAAAAA", leaf.cds, 12u));
        ASSERT_EQ(0, memcmp(reverse ? "FFFA" : "RKKK", leaf.protein, 4u));
        ASSERT_EQ(2u, s->projected_events); ASSERT_EQ(1u, s->completed_leaves);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    }
    PASS();
}

TEST haplotype_stream_edit_provenance_survives_islands_sorting_and_blocks(void) {
    /* Exhaust eight genomic positions independently: unchanged, part of one
     * uploaded MNV, or a separate SNV. The direct digit/run oracle never uses
     * the production splitter or sort. Repeat on both transcript strands. */
    for (unsigned reverse = 0u; reverse < 2u; reverse++) for (unsigned code = 1u; code < 6561u; code++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        if (reverse) { memset(f.reference, 'T', sizeof(f.reference)); f.strands[0] = -1; }
        uint8_t digits[8], alternate[8], expected_cds[12];
        unsigned value = code, has_mnv = 0u;
        size_t source_count = 0u;
        memset(expected_cds, reverse ? 'T' : 'A', sizeof(expected_cds));
        for (unsigned i = 0u; i < 8u; i++) {
            digits[i] = (uint8_t)(value % 3u); value /= 3u;
            alternate[i] = digits[i] == 1u ? 'C' : 'A';
            has_mnv |= digits[i] == 1u;
            source_count += digits[i] == 2u;
            if (digits[i]) expected_cds[reverse ? 11u - i : i] =
                digits[i] == 1u ? (reverse ? 'G' : 'C') : (reverse ? 'C' : 'G');
        }
        source_count += has_mnv;
        f.edit_event_ids[8] = UINT64_C(0xface0123456789ab);
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_carrier_key_t key = {0u, 0, 1u, 1u, 0u};
        uint32_t tx = 0u;
        if (has_mnv) {
            duckvep_haplotype_source_t source = {17u, (const uint8_t *)"AAAAAAAA", alternate,
                100u, 0u, 8u, 8u, 0u, 0u, 0u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
        for (unsigned i = 0u; i < 8u; i++) if (digits[i] == 2u) {
            duckvep_haplotype_source_t source = {UINT64_MAX - i, (const uint8_t *)"A",
                (const uint8_t *)"G", 100u + i, 0u, 1u, 1u, 0u, 0u, 0u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
        ASSERT_EQ(source_count, leaf.contributor_count);
        ASSERT_EQ(12u, leaf.cds_length);
        ASSERT_EQ(0, memcmp(expected_cds, leaf.cds, sizeof(expected_cds)));
        size_t edit = 0u;
        for (unsigned i = 0u; i < 8u;) {
            if (!digits[i]) { i++; continue; }
            unsigned begin = i++;
            if (digits[begin] == 1u) while (i < 8u && digits[i] == 1u) i++;
            ASSERT(edit < leaf.edit_count);
            size_t at = reverse ? leaf.edit_count - 1u - edit : edit;
            ASSERT_EQ(digits[begin] == 1u ? 17u : UINT64_MAX - begin, leaf.edit_event_ids[at]);
            ASSERT_EQ(reverse ? 13u - i : begin + 1u, f.edits[at].cds_start);
            ASSERT_EQ(i - begin, f.edits[at].ref_len);
            ASSERT_EQ(i - begin, f.edits[at].alt_len);
            edit++;
        }
        ASSERT_EQ(edit, leaf.edit_count);
        size_t covered = 0u;
        for (size_t i = 0u; i < leaf.block_count; i++) {
            ASSERT_EQ(covered, leaf.blocks[i].edit_begin);
            ASSERT(leaf.blocks[i].edit_count > 0u && leaf.blocks[i].edit_count <= edit - covered);
            covered += leaf.blocks[i].edit_count;
        }
        ASSERT_EQ(edit, covered);
        ASSERT_EQ(UINT64_C(0xface0123456789ab), f.edit_event_ids[8]);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    }
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    f.buffers.edit_event_ids = NULL;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG, duckvep_haplotype_stream_init(
        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
    PASS();
}

TEST haplotype_stream_calls_match_permitted_genotypes(void) {
    uint32_t cases = 0u;
    /* This is an additional finite lane, independent of the existing generated
     * carrier/edit properties. Enumerate genotype permutations rather than
     * calling the production phase reducer to construct expected evidence. */
    for (uint16_t ploidy = 1u; ploidy <= 4u; ploidy++) {
        for (uint32_t genotype = 0u; genotype < (1u << (2u * ploidy)); genotype++) {
            for (uint32_t bits = 0u; bits < (1u << ploidy); bits++) {
                int32_t alleles[4];
                uint8_t phase[4], order[4], possible[4] = {0}, missing_at[4] = {0};
                unsigned missing = 0u, unphased = 0u, homogeneous = 1u;
                for (uint16_t i = 0u; i < ploidy; i++) {
                    alleles[i] = (int32_t)((genotype >> (2u * i)) & 3u) - 1;
                    phase[i] = (uint8_t)((bits >> i) & 1u);
                    order[i] = (uint8_t)i;
                    if (alleles[i] < 0) missing++;
                    if (alleles[i] < 0 || alleles[i] != alleles[0]) homogeneous = 0u;
                    if (!phase[i]) unphased++;
                }
                do {
                    int allowed = 1;
                    for (uint16_t i = 0u; i < ploidy; i++)
                        if (phase[i] && order[i] != i) allowed = 0;
                    if (!allowed) continue;
                    for (uint16_t i = 0u; i < ploidy; i++) {
                        int32_t allele = alleles[order[i]];
                        possible[i] |= allele < 0 ? 7u : (uint8_t)(1u << allele);
                        if (allele < 0) missing_at[i] = 1u;
                    }
                } while (phase_next_permutation(order, ploidy));
                for (unsigned compat = 0u; compat < 2u; compat++) {
                    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
                        for (uint32_t alt = 1u; alt <= 2u; alt++) {
                            cases++;
                            struct haplotype_stream_scene f;
                            haplotype_stream_scene_prepare(&f, 1u);
                            if (reverse) {
                                memset(f.reference, 'T', sizeof(f.reference));
                                f.strands[0] = -1;
                            }
                            duckvep_haplotype_stream_t *s = &f.stream;
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                                s, &f.model, &f.exons, &f.sequences, &f.buffers));
                            duckvep_haplotype_phase_set_t sets[] = {{0, 0u}, {-3, 1u}, {0, 1u}, {17, 1u}};
                            size_t set_count = compat ? 1u : 4u;
                            duckvep_haplotype_call_t call = {alleles, bits ? phase : NULL, 0u,
                                alt, ploidy, {0, 1u}, compat ? DUCKVEP_PHASE_VEP116_COMPAT : DUCKVEP_PHASE_STRICT};
                            duckvep_haplotype_source_t source = {71u, (const uint8_t *)"A",
                                (const uint8_t *)(alt == 1u ? "C" : "G"), 103u, 0u, 1u, 1u, 0u, 0u, 0u};
                            uint32_t tx = 0u;
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                haplotype_test_begin_candidates(s, &source, &tx, 1u));
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                duckvep_haplotype_stream_push_call(s, tx, &call, sets, set_count));
                            uint8_t expected[4][4] = {{0}}, seen[4][4] = {{0}};
                            unsigned expected_calls = 0u;
                            if (compat) {
                                unsigned lane = 0u;
                                if (missing) for (unsigned i = 0u; i < ploidy; i++)
                                    expected[0][i] = DUCKVEP_CARRIER_MISSING;
                                for (unsigned i = 0u; i < ploidy; i++) {
                                    if (alleles[i] < 0) continue;
                                    if (alleles[i] == (int32_t)alt) expected[0][lane] |= DUCKVEP_CARRIER_CALLED;
                                    lane++;
                                }
                            } else {
                                for (unsigned i = 0u; i < ploidy; i++) {
                                    int resolved = ploidy == 1u || homogeneous || phase[i] || unphased == 1u ||
                                        (possible[i] & (possible[i] - 1u)) == 0u;
                                    uint8_t evidence = 0u;
                                    if (resolved) {
                                        if (alleles[i] < 0) evidence = DUCKVEP_CARRIER_MISSING;
                                        else if (alleles[i] == (int32_t)alt) evidence = DUCKVEP_CARRIER_CALLED;
                                    } else {
                                        if (possible[i] & (1u << alt)) evidence = DUCKVEP_CARRIER_UNPHASED;
                                        if (missing_at[i]) evidence |= DUCKVEP_CARRIER_MISSING;
                                    }
                                    for (unsigned set = 0u; set < 4u; set++) {
                                        if (set == 2u || ploidy == 1u || homogeneous || !bits)
                                            expected[set][i] = evidence;
                                    }
                                }
                            }
                            for (size_t set = 0u; set < set_count; set++)
                                for (unsigned i = 0u; i < ploidy; i++)
                                    if (expected[set][i]) expected_calls++;
                            ASSERT_EQ(expected_calls, s->carriers.call_count);
                            duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_finish(s);
                            unsigned seen_flags = 0u;
                            if (expected_calls) {
                                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, status);
                                duckvep_haplotype_leaf_t leaf;
                                while ((status = duckvep_haplotype_stream_next(s, &leaf)) ==
                                        DUCKVEP_HAPLOTYPE_STREAM_OK) {
                                    ASSERT(leaf.evidence_flags > 0u && leaf.evidence_flags < 8u);
                                    ASSERT_EQ(0u, seen_flags & (1u << leaf.evidence_flags));
                                    seen_flags |= 1u << leaf.evidence_flags;
                                    ASSERT_EQ(1u, leaf.contributor_count);
                                    ASSERT_EQ(71u, leaf.contributors[0].source.event_id);
                                    ASSERT_EQ(leaf.evidence_flags, leaf.contributors[0].evidence_flags);
                                    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                                    ASSERT_EQ(!!(leaf.evidence_flags & DUCKVEP_CARRIER_CALLED), leaf.edit_count);
                                    if (leaf.evidence_flags != DUCKVEP_CARRIER_CALLED) {
                                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE, leaf.sequence_status);
                                        ASSERT(!leaf.cds && !leaf.protein);
                                        ASSERT_EQ(0u, leaf.cds_length); ASSERT_EQ(0u, leaf.protein_length);
                                    } else {
                                        uint8_t cds[12]; memset(cds, reverse ? 'T' : 'A', sizeof(cds));
                                        cds[reverse ? 8u : 3u] = (uint8_t)(reverse ? (alt == 1u ? 'G' : 'C') :
                                                                                                 (alt == 1u ? 'C' : 'G'));
                                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
                                        ASSERT_EQ(12u, leaf.cds_length); ASSERT_EQ(4u, leaf.protein_length);
                                        ASSERT_EQ(0, memcmp(cds, leaf.cds, sizeof(cds)));
                                        const char *protein = reverse ? (alt == 1u ? "FFLF" : "FFFF") :
                                                                                        (alt == 1u ? "KQKK" : "KEKK");
                                        ASSERT_EQ(0, memcmp(protein, leaf.protein, 4u));
                                    }
                                    unsigned carriers = 0u;
                                    for (uint32_t id = leaf.carriers.first_call; id;) {
                                        const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, id);
                                        ASSERT(c); ASSERT_EQ(0u, c->key.sample_index);
                                        ASSERT_EQ(ploidy, c->key.ploidy);
                                        ASSERT(c->key.lane > 0u && c->key.lane <= ploidy);
                                        size_t set = 0u;
                                        while (set < set_count && (sets[set].present != c->key.phase_set_present ||
                                            (sets[set].present && sets[set].value != c->key.phase_set))) set++;
                                        ASSERT(set < set_count);
                                        ASSERT_EQ(0u, seen[set][c->key.lane - 1u]);
                                        seen[set][c->key.lane - 1u] = leaf.evidence_flags;
                                        carriers++; id = c->next_leaf;
                                    }
                                    ASSERT_EQ(carriers, leaf.carriers.call_count);
                                }
                                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
                                status = duckvep_haplotype_stream_finish(s);
                            }
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
                            ASSERT_EQ(0, memcmp(expected, seen, sizeof(expected)));
                            ASSERT_EQ((seen_flags & (1u << DUCKVEP_CARRIER_CALLED)) ? 12u : 0u,
                                      s->translated_bases);
                            ASSERT_EQ(1u, s->projected_events);
                        }
                    }
                }
            }
        }
    }
    ASSERT_EQ(37440u, cases);
    PASS();
}

TEST haplotype_stream_preserves_later_phase_sets_and_uncertain_prefixes(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    const duckvep_haplotype_phase_set_t sets[] = {{10, 1u}, {20, 1u}};
    uint8_t expected[3][2][2][3] = {{{{0}}}}, seen[3][2][2] = {{{0}}};
    uint32_t tx = 0u;
    for (uint32_t event = 0u; event < 3u; event++) {
        duckvep_haplotype_source_t source = {event + 1u, (const uint8_t *)"A",
            (const uint8_t *)(event == 1u ? "G" : "C"), 100u + event, 0u, 1u, 1u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        for (uint32_t sample = 0u; sample < 3u; sample++) {
            int32_t gt[2] = {event == 2u ? 0 : 1, event == 1u ? 0 : 1};
            uint8_t phase[2] = {1u, 1u};
            if (event == 1u && sample == 2u) gt[1] = -1;
            if (event == 2u && sample == 1u) phase[0] = phase[1] = 0u;
            duckvep_haplotype_call_t call = {gt, phase, sample, 1u, 2u,
                event ? sets[event - 1u] : (duckvep_haplotype_phase_set_t){0, 0u}, DUCKVEP_PHASE_STRICT};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_call(s, tx, &call, sets, 2u));
            for (unsigned set = 0u; set < 2u; set++) {
                for (unsigned lane = 0u; lane < 2u; lane++) {
                    if (!event) expected[sample][set][lane][event] = DUCKVEP_CARRIER_CALLED;
                    else if (event == 1u && set == 0u) {
                        if (lane == 0u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_CALLED;
                        else if (sample == 2u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_MISSING;
                    } else if (event == 2u) {
                        if (sample == 1u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_UNPHASED;
                        else if (set == 1u && lane == 1u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_CALLED;
                    }
                }
            }
        }
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    duckvep_haplotype_stream_status_t status;
    unsigned leaves = 0u, calls = 0u;
    while ((status = duckvep_haplotype_stream_next(s, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        leaves++;
        for (uint32_t id = leaf.carriers.first_call; id;) {
            const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, id);
            ASSERT(c && c->key.sample_index < 3u && c->key.phase_set_present);
            ASSERT(c->key.phase_set == 10 || c->key.phase_set == 20);
            ASSERT(c->key.lane == 1u || c->key.lane == 2u);
            ASSERT_EQ(2u, c->key.ploidy);
            unsigned sample = c->key.sample_index, set = c->key.phase_set == 20, lane = c->key.lane - 1u;
            ASSERT_EQ(0u, seen[sample][set][lane]); seen[sample][set][lane] = 1u;
            unsigned count = 0u;
            uint8_t evidence = 0u, cds[12]; memset(cds, 'A', sizeof(cds));
            for (unsigned event = 0u; event < 3u; event++) {
                uint8_t flags = expected[sample][set][lane][event];
                if (!flags) continue;
                ASSERT(count < leaf.contributor_count);
                ASSERT_EQ(event + 1u, leaf.contributors[count].source.event_id);
                ASSERT_EQ(flags, leaf.contributors[count].evidence_flags);
                evidence |= flags; count++;
                if (flags & DUCKVEP_CARRIER_CALLED) cds[event] = event == 1u ? 'G' : 'C';
            }
            ASSERT_EQ(count, leaf.contributor_count); ASSERT_EQ(evidence, leaf.evidence_flags);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            if (evidence == DUCKVEP_CARRIER_CALLED) {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
                ASSERT_EQ(12u, leaf.cds_length);
                ASSERT_EQ(0, memcmp(cds, leaf.cds, sizeof(cds)));
            } else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE, leaf.sequence_status);
                ASSERT(!leaf.cds && !leaf.protein);
            }
            calls++; id = c->next_leaf;
        }
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
    ASSERT_EQ(12u, calls); ASSERT_EQ(6u, leaves);
    ASSERT_EQ(3u, s->projected_events); ASSERT_EQ(36u, s->translated_bases);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    PASS();
}

TEST haplotype_stream_rejects_invalid_calls_and_latches_partial_broadcast(void) {
    for (unsigned scenario = 0u; scenario < 15u; scenario++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        if (scenario == 14u) f.buffers.carriers.call_capacity = 1u;
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        uint32_t tx = 0u;
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)"A", (const uint8_t *)"C",
            100u, 0u, 1u, 1u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        int32_t gt[2] = {1, 0};
        uint8_t phase[2] = {1u, 1u};
        duckvep_haplotype_phase_set_t sets[] = {{0, 0u}, {10, 1u}};
        duckvep_haplotype_call_t call = {gt, phase, 0u, 1u, 2u, {10, 1u}, DUCKVEP_PHASE_STRICT};
        if (scenario == 0u) call.alleles = NULL;
        if (scenario == 1u) call.ploidy = 0u;
        if (scenario == 2u) call.alt_index = 0u;
        if (scenario == 3u) call.alt_index = UINT32_MAX;
        if (scenario == 4u) call.phase_set.present = 2u;
        if (scenario == 5u) call.policy = (duckvep_phase_policy_t)99;
        if (scenario == 6u) sets[1].present = 0u; /* Repeated absent set. */
        if (scenario == 7u) { sets[0] = sets[1]; sets[1].value = 5; }
        if (scenario == 8u) call.phase_set.value = 11;
        if (scenario == 9u) gt[1] = -2;
        if (scenario == 10u) phase[1] = 2u;
        if (scenario == 11u) call.policy = DUCKVEP_PHASE_VEP116_COMPAT;
        if (scenario == 12u) tx = 1u; /* Not projected by begin. */
        if (scenario == 13u) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_call(s, tx, &call, sets, 2u));
            call.policy = DUCKVEP_PHASE_VEP116_COMPAT;
        }
        if (scenario == 14u) gt[1] = 1; /* First lane lands, second exceeds the call pool. */
        duckvep_haplotype_stream_status_t expected = scenario == 14u
            ? DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR : DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        ASSERT_EQ(expected, duckvep_haplotype_stream_push_call(s, tx, &call, sets,
            scenario == 13u ? 1u : 2u));
        ASSERT_EQ(scenario >= 13u ? 1u : 0u, s->carriers.call_count);
        if (scenario == 14u) ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, s->carrier_error);
        ASSERT_EQ(expected, duckvep_haplotype_stream_push_call(s, tx, &call, NULL, 0u));
        ASSERT_EQ(expected, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf, zero = {0}; memset(&leaf, 0xa5, sizeof(leaf));
        ASSERT_EQ(expected, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(0, memcmp(&leaf, &zero, sizeof(leaf)));
    }
    PASS();
}

/* An independent dense bit matrix is the oracle, not the prefix index. Each
 * bit names a carried event; the empty mask is an implicit reference path. */
static int carrier_test_drain(duckvep_carriers_t *s, uint32_t tx,
                              const uint8_t masks[8], uint64_t event_base) {
    uint64_t seen_paths = 0u;
    unsigned seen_samples = 0u;
    duckvep_carrier_leaf_t leaf;
    duckvep_carriers_status_t status;
    while ((status = duckvep_carriers_next_leaf(s, &leaf)) == DUCKVEP_CARRIERS_OK) {
        duckvep_carrier_event_t events[6], before[6];
        size_t required = 0u;
        memset(events, 0xa5, sizeof events);
        memcpy(before, events, sizeof before);
        if (leaf.transcript_index != tx || !leaf.event_count || !leaf.call_count ||
            duckvep_carriers_leaf_events(s, leaf.id, events, 0u, &required) !=
                DUCKVEP_CARRIERS_OUTPUT_FULL || required != leaf.event_count ||
            memcmp(events, before, sizeof before) != 0 ||
            duckvep_carriers_leaf_events(s, leaf.id, events, 6u, &required) !=
                DUCKVEP_CARRIERS_OK) return 0;
        unsigned mask = 0u;
        for (size_t i = 0u; i < required; i++) {
            if (events[i].evidence_flags != DUCKVEP_CARRIER_CALLED ||
                events[i].event_id < event_base || events[i].event_id >= event_base + 6u ||
                (i && events[i].event_id <= events[i - 1u].event_id)) return 0;
            mask |= 1u << (unsigned)(events[i].event_id - event_base);
        }
        if (!mask || (seen_paths & (UINT64_C(1) << mask))) return 0;
        seen_paths |= UINT64_C(1) << mask;
        uint32_t count = 0u;
        for (uint32_t id = leaf.first_call; id;) {
            const duckvep_carrier_call_t *call = duckvep_carriers_call(s, id);
            if (!call || ++count > 8u || call->key.sample_index >= 8u) return 0;
            uint32_t sample = call->key.sample_index;
            duckvep_carrier_key_t key = carrier_test_key(sample);
            if ((seen_samples & (1u << sample)) || masks[sample] != mask ||
                call->key.lane != key.lane || call->key.ploidy != key.ploidy ||
                call->key.phase_set != key.phase_set ||
                call->key.phase_set_present != key.phase_set_present) return 0;
            seen_samples |= 1u << sample;
            id = call->next_leaf;
        }
        if (count != leaf.call_count) return 0;
    }
    if (status != DUCKVEP_CARRIERS_DONE) return 0;
    for (unsigned sample = 0u; sample < 8u; sample++) {
        if (!!(seen_samples & (1u << sample)) != !!masks[sample]) return 0;
    }
    return duckvep_carriers_release(s) == DUCKVEP_CARRIERS_OK;
}

TEST carrier_stream_lifetime_keys_and_capacity(void) {
    uint16_t chrom[] = {1u, 1u, 2u};
    uint32_t end[] = {3u, 8u, 4u};
    duckvep_transcript_model_t model = {0};
    model.transcript_count = 3u; model.chrom_id = chrom; model.end1 = end;
    struct carrier_test_pool pool, before_pool;
    duckvep_carrier_buffers_t b = carrier_test_buffers(&pool);
    duckvep_carriers_t s, before;
    uint32_t completed;
    duckvep_carrier_key_t key = carrier_test_key(0u);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 1u, 10u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 1u, &key, DUCKVEP_CARRIER_CALLED));
    /* NULL PS and an explicitly present zero are different carrier keys. */
    key.phase_set_present = 1u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(3u, s.call_count);
    ASSERT_EQ(2u, s.prefix_count);
    before = s; memcpy(&before_pool, &pool, sizeof pool);
    ASSERT_EQ(DUCKVEP_CARRIERS_DUPLICATE_CALL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    key.ploidy = 2u;
    ASSERT_EQ(DUCKVEP_CARRIERS_PLOIDY_CONFLICT, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0, memcmp(&before, &s, sizeof s));
    ASSERT_EQ(0, memcmp(&before_pool, &pool, sizeof pool));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 3u, 11u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY,
              duckvep_carriers_advance(&s, 1u, 4u, 12u, &completed));
    ASSERT_EQ(0u, completed);
    ASSERT_EQ(DUCKVEP_CARRIERS_INPUT_ORDER,
              duckvep_carriers_advance(&s, 1u, 5u, 12u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_release(&s));
    duckvep_carrier_leaf_t leaf;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(2u, leaf.call_count);
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    ASSERT_EQ(1u, s.transcript_count); ASSERT_EQ(1u, s.call_count); ASSERT_EQ(1u, s.prefix_count);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 4u, 12u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_INPUT_ORDER,
              duckvep_carriers_advance(&s, 1u, 4u, 11u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY,
              duckvep_carriers_advance(&s, 2u, 1u, 13u, &completed));
    ASSERT_EQ(1u, completed);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 2u, 1u, 13u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 2u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, duckvep_carriers_finish(&s, &completed));
    ASSERT_EQ(2u, completed);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_finish(&s, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_finish(&s, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_INPUT_ORDER,
              duckvep_carriers_advance(&s, 2u, 2u, 14u, &completed));
    ASSERT_EQ(0u, s.transcript_count); ASSERT_EQ(0u, s.call_count); ASSERT_EQ(0u, s.prefix_count);

    b.transcript_capacity = 1u; b.transcript_buckets = 2u;
    b.call_capacity = 2u; b.call_buckets = 4u;
    b.prefix_capacity = 2u; b.prefix_buckets = 4u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 1u, 1u, &completed));
    for (uint32_t sample = 0u; sample < 2u; sample++) {
        key = carrier_test_key(sample);
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    }
    before = s; memcpy(&before_pool, &pool, sizeof pool);
    key = carrier_test_key(2u);
    ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_FULL, duckvep_carriers_push(&s, 1u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0, memcmp(&before, &s, sizeof s));
    ASSERT_EQ(0, memcmp(&before_pool, &pool, sizeof pool));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 2u, 2u, &completed));
    for (uint32_t sample = 0u; sample < 2u; sample++) {
        key = carrier_test_key(sample);
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    }
    ASSERT_EQ(2u, s.prefix_count); /* Reuse succeeds even when no slot is free. */
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 3u, 3u, &completed));
    before = s; memcpy(&before_pool, &pool, sizeof pool);
    ASSERT_EQ(DUCKVEP_CARRIERS_PREFIX_FULL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0, memcmp(&before, &s, sizeof s));
    ASSERT_EQ(0, memcmp(&before_pool, &pool, sizeof pool));
    b.call_capacity = 0u; b.call_buckets = 0u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 1u, 1u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0u, s.transcript_count); ASSERT_EQ(0u, s.prefix_count);
    b.transcript_buckets = 3u;
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_init(&s, &model, &b));
    PASS();
}

struct carrier_matrix_case { uint8_t masks[2][8]; };

TEST carrier_stream_reuses_slots_after_many_transcripts(void) {
    enum { TRANSCRIPTS = 1000 };
    uint16_t chrom[TRANSCRIPTS] = {0};
    uint32_t end[TRANSCRIPTS];
    for (unsigned i = 0u; i < TRANSCRIPTS; i++) end[i] = i + 1u;
    duckvep_transcript_model_t model = {0};
    model.transcript_count = TRANSCRIPTS; model.chrom_id = chrom; model.end1 = end;
    struct carrier_test_pool pool;
    duckvep_carrier_buffers_t b = carrier_test_buffers(&pool);
    b.transcript_capacity = 1u; b.transcript_buckets = 2u;
    b.call_capacity = 8u; b.call_buckets = 16u;
    b.prefix_capacity = 1u; b.prefix_buckets = 2u;
    duckvep_carriers_t s;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    for (uint32_t tx = 0u; tx <= TRANSCRIPTS; tx++) {
        uint32_t completed;
        duckvep_carriers_status_t status = tx == TRANSCRIPTS
            ? duckvep_carriers_finish(&s, &completed)
            : duckvep_carriers_advance(&s, 0u, tx + 1u, UINT64_MAX - tx, &completed);
        if (tx) {
            ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, status);
            ASSERT_EQ(tx - 1u, completed);
            duckvep_carrier_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(8u, leaf.call_count);
            duckvep_carrier_event_t event;
            size_t required;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK,
                      duckvep_carriers_leaf_events(&s, leaf.id, &event, 1u, &required));
            ASSERT_EQ(1u, required);
            ASSERT_EQ(UINT64_MAX - (tx - 1u), event.event_id);
            ASSERT_EQ(DUCKVEP_CARRIER_CALLED, event.evidence_flags);
            unsigned seen = 0u;
            for (uint32_t id = leaf.first_call; id;) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s, id);
                ASSERT(call != NULL);
                uint32_t sample = call->key.sample_index - (tx - 1u) * 8u;
                ASSERT(sample < 8u);
                ASSERT(!(seen & (1u << sample)));
                seen |= 1u << sample;
                ASSERT_EQ(sample == 7u ? UINT16_MAX : sample + 1u, call->key.ploidy);
                ASSERT_EQ(call->key.ploidy, call->key.lane);
                id = call->next_leaf;
            }
            ASSERT_EQ(255u, seen);
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
            ASSERT_EQ(0u, s.transcript_count); ASSERT_EQ(0u, s.call_count); ASSERT_EQ(0u, s.prefix_count);
            status = tx == TRANSCRIPTS ? duckvep_carriers_finish(&s, &completed)
                : duckvep_carriers_advance(&s, 0u, tx + 1u, UINT64_MAX - tx, &completed);
        }
        if (tx == TRANSCRIPTS) {
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, status);
            break;
        }
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, status);
        for (uint32_t sample = 0u; sample < 8u; sample++) {
            duckvep_carrier_key_t key = carrier_test_key(tx * 8u + sample);
            key.ploidy = sample == 7u ? UINT16_MAX : (uint16_t)(sample + 1u);
            key.lane = key.ploidy;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, tx, &key, DUCKVEP_CARRIER_CALLED));
        }
    }
    ASSERT_EQ(1u, s.peak_transcripts); ASSERT_EQ(8u, s.peak_calls); ASSERT_EQ(1u, s.peak_prefixes);
    PASS();
}

TEST carrier_stream_expiry_heap_matches_dense_active_set(void) {
    enum { TRANSCRIPTS = 256, SLOTS = 64 };
    uint16_t chrom[TRANSCRIPTS] = {0};
    uint32_t end[TRANSCRIPTS], start[TRANSCRIPTS], active[SLOTS];
    uint8_t live[TRANSCRIPTS] = {0};
    duckvep_carrier_transcript_t transcripts[SLOTS];
    duckvep_carrier_call_t calls[SLOTS];
    duckvep_carrier_prefix_t prefixes[SLOTS];
    duckvep_carrier_bucket_t tx_index[2 * SLOTS], call_index[2 * SLOTS], prefix_index[2 * SLOTS];
    duckvep_carrier_buffers_t b = {transcripts, calls, prefixes, active,
        tx_index, call_index, prefix_index, SLOTS, SLOTS, SLOTS,
        2 * SLOTS, 2 * SLOTS, 2 * SLOTS};
    /* Model ordinals are not opening order. First fill every slot with tied
     * ends; then interleave openings and expiry with up to 64 live positions. */
    for (uint32_t pos = 1u; pos <= TRANSCRIPTS; pos++) {
        uint32_t tx = (pos * 73u) % TRANSCRIPTS;
        start[tx] = pos;
        end[tx] = pos <= SLOTS ? SLOTS : pos + (tx * 29u) % SLOTS;
    }
    duckvep_transcript_model_t model = {0};
    model.transcript_count = TRANSCRIPTS; model.chrom_id = chrom; model.end1 = end;
    duckvep_carriers_t s;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    uint32_t drained = 0u;
    for (uint32_t pos = 1u; pos <= TRANSCRIPTS + 1u; pos++) {
        int eof = pos > TRANSCRIPTS;
        for (;;) {
            uint32_t expected = UINT32_MAX, completed;
            for (uint32_t tx = 0u; tx < TRANSCRIPTS; tx++) {
                if (live[tx] && (expected == UINT32_MAX || end[tx] < end[expected])) expected = tx;
            }
            duckvep_carriers_status_t status = eof ? duckvep_carriers_finish(&s, &completed)
                : duckvep_carriers_advance(&s, 0u, pos, pos, &completed);
            if (expected == UINT32_MAX || (!eof && end[expected] >= pos)) {
                ASSERT_EQ(eof ? DUCKVEP_CARRIERS_DONE : DUCKVEP_CARRIERS_OK, status);
                break;
            }
            ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, status);
            ASSERT_EQ(expected, completed);
            duckvep_carrier_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(expected, leaf.transcript_index);
            ASSERT_EQ(1u, leaf.call_count);
            duckvep_carrier_event_t event;
            size_t required;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK,
                      duckvep_carriers_leaf_events(&s, leaf.id, &event, 1u, &required));
            ASSERT_EQ(1u, required); ASSERT_EQ(start[expected], event.event_id);
            ASSERT_EQ(DUCKVEP_CARRIER_CALLED, event.evidence_flags);
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&s, leaf.first_call);
            ASSERT(call != NULL); ASSERT_EQ(expected, call->key.sample_index);
            ASSERT_EQ(0u, call->next_leaf);
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
            live[expected] = 0u;
            drained++;
        }
        if (!eof) {
            uint32_t tx = (pos * 73u) % TRANSCRIPTS, completed;
            duckvep_carrier_key_t key = carrier_test_key(tx);
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, tx, &key, DUCKVEP_CARRIER_CALLED));
            live[tx] = 1u;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 0u, pos, pos, &completed));
        }
        uint32_t live_count = 0u, earliest = UINT32_MAX;
        for (uint32_t tx = 0u; tx < TRANSCRIPTS; tx++) if (live[tx]) {
            live_count++;
            if (earliest == UINT32_MAX || end[tx] < end[earliest]) earliest = tx;
        }
        ASSERT_EQ(live_count, s.transcript_count);
        ASSERT_EQ(live_count, s.call_count); ASSERT_EQ(live_count, s.prefix_count);
        if (live_count) {
            ASSERT_EQ(earliest, transcripts[active[0] - 1u].transcript_index);
            for (uint32_t child = 1u; child < live_count; child++) {
                uint32_t a = transcripts[active[(child - 1u) / 2u] - 1u].transcript_index;
                uint32_t z = transcripts[active[child] - 1u].transcript_index;
                ASSERT(end[a] < end[z] || (end[a] == end[z] && a < z));
            }
        }
    }
    ASSERT_EQ(TRANSCRIPTS, drained);
    ASSERT_EQ(SLOTS, s.peak_transcripts); ASSERT_EQ(SLOTS, s.peak_calls); ASSERT_EQ(SLOTS, s.peak_prefixes);
    PASS();
}

static enum theft_alloc_res carrier_matrix_alloc(struct theft *t, void *env, void **instance) {
    (void)env;
    struct carrier_matrix_case *c = malloc(sizeof(*c));
    if (!c) return THEFT_ALLOC_ERROR;
    for (unsigned tx = 0u; tx < 2u; tx++) for (unsigned sample = 0u; sample < 8u; sample++) {
        c->masks[tx][sample] = (uint8_t)kprop_bounded(t, tx ? 64u : 16u);
    }
    /* Include cohort prefix sharing on every trial, not just by chance. */
    for (unsigned tx = 0u; tx < 2u; tx++) c->masks[tx][7] = c->masks[tx][0];
    *instance = c;
    return THEFT_ALLOC_OK;
}

static void carrier_matrix_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static enum theft_trial_res prop_carrier_stream_matches_dense_matrix(struct theft *t, void *arg) {
    (void)t;
    const struct carrier_matrix_case *c = arg;
    uint16_t chrom[2][2] = {{1u, 1u}, {2u, 2u}};
    uint32_t end[2][2] = {{3u, 6u}, {4u, 6u}};
    duckvep_transcript_model_t models[2] = {{0}, {0}};
    uint8_t expected[2][2][8];
    struct carrier_test_pool pools[2];
    duckvep_carriers_t streams[2];
    const uint64_t base = UINT64_MAX - 6u;
    for (unsigned run = 0u; run < 2u; run++) {
        models[run].transcript_count = 2u;
        models[run].chrom_id = chrom[run]; models[run].end1 = end[run];
        memcpy(expected[run], c->masks, sizeof c->masks);
        for (unsigned sample = 0u; sample < 8u; sample++) expected[run][0][sample] &= run ? 15u : 7u;
        duckvep_carrier_buffers_t b = carrier_test_buffers(&pools[run]);
        if (duckvep_carriers_init(&streams[run], &models[run], &b) != DUCKVEP_CARRIERS_OK)
            return THEFT_TRIAL_ERROR;
    }
    /* Interleave different models with identical local ordinals but different
     * contigs/transcript endpoints. One also resumes between every carrier row. */
    for (uint32_t pos = 1u; pos <= 6u; pos++) for (unsigned run = 0u; run < 2u; run++) {
        duckvep_carriers_t *s = &streams[run];
        uint32_t tx;
        duckvep_carriers_status_t status;
        while ((status = duckvep_carriers_advance(s, (uint16_t)(run + 1u), pos, base + pos - 1u, &tx)) ==
                   DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
            if (!carrier_test_drain(s, tx, expected[run][tx], base)) return THEFT_TRIAL_FAIL;
        }
        if (status != DUCKVEP_CARRIERS_OK) return THEFT_TRIAL_FAIL;
        for (tx = 0u; tx < 2u; tx++) for (uint32_t sample = 0u; sample < 8u; sample++) {
            if (!(expected[run][tx][sample] & (1u << (pos - 1u)))) continue;
            uint32_t completed;
            if (run && duckvep_carriers_advance(s, (uint16_t)(run + 1u), pos, base + pos - 1u, &completed) !=
                           DUCKVEP_CARRIERS_OK) return THEFT_TRIAL_FAIL;
            duckvep_carrier_key_t key = carrier_test_key(sample);
            if (duckvep_carriers_push(s, tx, &key, DUCKVEP_CARRIER_CALLED) != DUCKVEP_CARRIERS_OK) return THEFT_TRIAL_FAIL;
        }
    }
    for (unsigned run = 0u; run < 2u; run++) {
        duckvep_carriers_t *s = &streams[run];
        uint32_t tx;
        duckvep_carriers_status_t status;
        while ((status = duckvep_carriers_finish(s, &tx)) == DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
            if (!carrier_test_drain(s, tx, expected[run][tx], base)) return THEFT_TRIAL_FAIL;
        }
        if (status != DUCKVEP_CARRIERS_DONE || s->transcript_count || s->call_count || s->prefix_count ||
            s->peak_transcripts > 2u || s->peak_calls > 16u || s->peak_prefixes > 64u)
            return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST carrier_stream_matches_dense_matrix_across_batches(void) {
    struct theft_type_info type = {.alloc = carrier_matrix_alloc, .free = carrier_matrix_free};
    struct theft_run_config cfg = {0};
    cfg.name = "sparse carrier paths == dense event matrix across input batches";
    cfg.prop1 = prop_carrier_stream_matches_dense_matrix;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* The dense genotype masks and direct genomic-to-CDS index formula are
 * independent of the replay store, projector and prefix traversal. The same
 * model-local ordinal and event IDs coexist on opposite strands/contigs. */
static enum theft_trial_res prop_haplotype_stream_matches_dense_models(struct theft *t, void *arg) {
    (void)t;
    const struct carrier_matrix_case *c = arg;
    struct haplotype_stream_scene f[2];
    for (unsigned run = 0u; run < 2u; run++) {
        haplotype_stream_scene_prepare(&f[run], 1u);
        f[run].chrom[0] = (uint16_t)(run + 1u);
        if (run) {
            memset(f[run].reference, 'T', sizeof(f[run].reference));
            f[run].strands[0] = -1;
        }
        if (duckvep_haplotype_stream_init(&f[run].stream, &f[run].model,
            &f[run].exons, &f[run].sequences, &f[run].buffers) != DUCKVEP_HAPLOTYPE_STREAM_OK)
            return THEFT_TRIAL_ERROR;
    }
    uint32_t tx = 0u;
    for (uint32_t event = 0u; event < 6u; event++) for (unsigned run = 0u; run < 2u; run++) {
        uint8_t ref = 'A', alt = 'C';
        duckvep_haplotype_source_t source = {UINT64_MAX - 5u + event, &ref, &alt,
            100u + event, (uint16_t)(run + 1u), 1u, 1u, 0u, 0u, 0u};
        duckvep_haplotype_stream_t *s = &f[run].stream;
        if (haplotype_test_begin_candidates(s, &source, &tx, 1u) != DUCKVEP_HAPLOTYPE_STREAM_OK)
            return THEFT_TRIAL_FAIL;
        ref = alt = 'N';
        for (uint32_t sample = 0u; sample < 8u; sample++) {
            if (!(c->masks[run][sample] & (1u << event))) continue;
            duckvep_carrier_key_t key = carrier_test_key(sample);
            if (haplotype_test_push_called(s, &key) != DUCKVEP_HAPLOTYPE_STREAM_OK)
                return THEFT_TRIAL_FAIL;
        }
    }
    for (unsigned run = 0u; run < 2u; run++) {
        duckvep_haplotype_stream_t *s = &f[run].stream;
        unsigned seen_samples = 0u, path_count = 0u;
        uint64_t seen_paths = 0u;
        uint8_t reference_peptide[5];
        duckvep_translation_t reference_translation;
        if (duckvep_translate_cds(f[run].reference, sizeof(f[run].reference),
                DUCKVEP_CODON_TABLE_STANDARD, DUCKVEP_TRANSLATION_N_CONSENSUS,
                reference_peptide, sizeof(reference_peptide), &reference_translation) !=
            DUCKVEP_TRANSLATION_OK) return THEFT_TRIAL_ERROR;
        duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_finish(s);
        if (s->carriers.call_count) {
            if (status != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) return THEFT_TRIAL_FAIL;
            duckvep_haplotype_leaf_t leaf;
            while ((status = duckvep_haplotype_stream_next(s, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                if (leaf.carriers.transcript_index || leaf.projection_status != DUCKVEP_CDS_EDIT_OK ||
                    leaf.sequence_status != DUCKVEP_HAPLOTYPE_OK || leaf.cds_length != 12u ||
                    leaf.protein_length != 4u || leaf.flags || !leaf.carriers.call_count)
                    return THEFT_TRIAL_FAIL;
                unsigned mask = 0u;
                uint8_t expected[12]; memset(expected, run ? 'T' : 'A', sizeof(expected));
                for (size_t i = 0u; i < leaf.contributor_count; i++) {
                    const duckvep_haplotype_source_t *source = &leaf.contributors[i].source;
                    if (source->event_id < UINT64_MAX - 5u || source->chrom_id != run + 1u ||
                        (i && source->event_id <= leaf.contributors[i - 1u].source.event_id) ||
                        source->ref[0] != 'A' || source->alt[0] != 'C') return THEFT_TRIAL_FAIL;
                    unsigned event = (unsigned)(source->event_id - (UINT64_MAX - 5u));
                    if (source->pos1 != 100u + event) return THEFT_TRIAL_FAIL;
                    mask |= 1u << event;
                    expected[run ? 11u - event : event] = run ? 'G' : 'C';
                }
                if (!mask || memcmp(expected, leaf.cds, sizeof(expected)) ||
                    (seen_paths & (UINT64_C(1) << mask))) return THEFT_TRIAL_FAIL;
                /* Consume the stream's actual completed leaf, including its
                 * ascending edits, without a second apply or translation. */
                duckvep_edit_set_t set = {s->buffers.edits, leaf.edit_count};
                duckvep_haplotype_result_t applied = {leaf.cds_length, 0, leaf.flags, leaf.edit_count};
                duckvep_coding_context_t context;
                if (duckvep_coding_context_open_replay(leaf.reference_cds, 12u, &set,
                        f[run].strands[0], DUCKVEP_CODON_TABLE_STANDARD, leaf.cds, &applied,
                        reference_peptide, &reference_translation, leaf.protein, &leaf.translation,
                        &context) != DUCKVEP_CODING_CONTEXT_OK ||
                    context.alt_cds != leaf.cds || context.alt_peptide != leaf.protein ||
                    context.ref_peptide != reference_peptide || context.applied_edits != leaf.edit_count)
                    return THEFT_TRIAL_FAIL;
                for (size_t block = 0u; block < leaf.block_count; block++) {
                    duckvep_sequence_delta_t facts;
                    if (duckvep_coding_context_block_delta_fill(&context, set.edits, set.count,
                            leaf.blocks + block, 0u, &facts) != DUCKVEP_CONTEXT_DELTA_OK || !facts.valid)
                        return THEFT_TRIAL_FAIL;
                }
                seen_paths |= UINT64_C(1) << mask;
                path_count++;
                unsigned count = 0u;
                for (uint32_t id = leaf.carriers.first_call; id;) {
                    const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, id);
                    if (!call || call->key.sample_index >= 8u || ++count > 8u) return THEFT_TRIAL_FAIL;
                    unsigned sample = call->key.sample_index;
                    duckvep_carrier_key_t key = carrier_test_key(sample);
                    if ((seen_samples & (1u << sample)) || c->masks[run][sample] != mask ||
                        call->key.lane != key.lane || call->key.ploidy != key.ploidy ||
                        call->key.phase_set_present != key.phase_set_present ||
                        call->key.phase_set != (key.phase_set_present ? key.phase_set : 0))
                        return THEFT_TRIAL_FAIL;
                    seen_samples |= 1u << sample;
                    id = call->next_leaf;
                }
                if (count != leaf.carriers.call_count) return THEFT_TRIAL_FAIL;
            }
            if (status != DUCKVEP_HAPLOTYPE_STREAM_DONE) return THEFT_TRIAL_FAIL;
            status = duckvep_haplotype_stream_finish(s);
        }
        if (status != DUCKVEP_HAPLOTYPE_STREAM_DONE || s->event_count || s->projection_count ||
            s->allele_count || s->completed_leaves != path_count || s->projected_events != 6u ||
            s->translated_bases != 12u * path_count) return THEFT_TRIAL_FAIL;
        for (unsigned sample = 0u; sample < 8u; sample++)
            if (!!(seen_samples & (1u << sample)) != !!c->masks[run][sample]) return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST haplotype_stream_matches_dense_models_across_batches(void) {
    struct theft_type_info type = {.alloc = carrier_matrix_alloc, .free = carrier_matrix_free};
    struct theft_run_config cfg = {0};
    cfg.name = "owned haplotype replay == dense genomic edits in coexisting models";
    cfg.prop1 = prop_haplotype_stream_matches_dense_models;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}
