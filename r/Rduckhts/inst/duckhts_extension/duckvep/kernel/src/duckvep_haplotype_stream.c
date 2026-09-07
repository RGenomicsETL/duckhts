#include "duckvep_haplotype_stream.h"

#include <string.h>

static duckvep_haplotype_stream_status_t fail(
    duckvep_haplotype_stream_t *s, duckvep_haplotype_stream_status_t status) {
    if (s) s->error = status;
    return status;
}

static duckvep_haplotype_stream_status_t carrier_fail(
    duckvep_haplotype_stream_t *s, duckvep_carriers_status_t status) {
    s->carrier_error = status;
    return fail(s, status == DUCKVEP_CARRIERS_INPUT_ORDER
        ? DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER : DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR);
}

/* Both arguments are at most capacity; avoid overflowing an address-size sum. */
static size_t ring_add(size_t begin, size_t count, size_t capacity) {
    return count >= capacity - begin ? count - (capacity - begin) : begin + count;
}

static int valid_array(const void *p, size_t n, size_t width) {
    return p && n && n <= SIZE_MAX / width;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_init(
    duckvep_haplotype_stream_t *s, const duckvep_transcript_model_t *tx,
    const duckvep_exon_model_t *exons, const duckvep_sequence_pool_t *seq,
    const duckvep_haplotype_stream_buffers_t *b) {
    if (!s) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    memset(s, 0, sizeof(*s));
    if (!tx || !exons || !seq || !b || !tx->start1 || !tx->strand ||
        seq->transcript_count != tx->transcript_count ||
        !seq->cds_offset || !seq->cds_length || !seq->cds_bytes ||
        !valid_array(b->events, b->event_capacity, sizeof(*b->events)) ||
        !valid_array(b->projections, b->projection_capacity, sizeof(*b->projections)) ||
        !valid_array(b->alleles, b->allele_capacity, 1u) ||
        !valid_array(b->leaf_events, b->leaf_capacity, sizeof(*b->leaf_events)) ||
        !valid_array(b->contributors, b->leaf_capacity, sizeof(*b->contributors)) ||
        !valid_array(b->edits, b->edit_capacity, sizeof(*b->edits)) ||
        !valid_array(b->cds, b->cds_capacity, 1u) ||
        !valid_array(b->protein, b->protein_capacity, 1u))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    duckvep_carriers_status_t status = duckvep_carriers_init(&s->carriers, tx, &b->carriers);
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    s->exons = exons;
    s->sequences = seq;
    s->buffers = *b;
    s->initialized = 1u;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

/* Called only after the carrier watermark has drained every earlier transcript.
 * An old event can pin younger events behind it, but nothing beyond the active
 * genomic window is retained and all three rings are reclaimed together. */
static void reclaim(duckvep_haplotype_stream_t *s, uint16_t chrom, uint32_t pos1, int eof) {
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    while (s->event_count) {
        const duckvep_haplotype_stored_event_t *event = &b->events[s->event_begin];
        if (!eof && (event->source.chrom_id > chrom ||
            (event->source.chrom_id == chrom && event->last_end1 >= pos1))) break;
        s->allele_begin = ring_add(s->allele_begin, event->allele_consumed, b->allele_capacity);
        s->allele_count -= event->allele_consumed;
        s->projection_begin = (uint32_t)ring_add(s->projection_begin,
            event->projection_count, b->projection_capacity);
        s->projection_count -= event->projection_count;
        s->event_begin = (uint32_t)ring_add(s->event_begin, 1u, b->event_capacity);
        s->event_count--;
    }
    if (!s->event_count) {
        s->event_begin = s->projection_begin = 0u;
        s->allele_begin = 0u;
    }
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_begin(
    duckvep_haplotype_stream_t *s, const duckvep_haplotype_source_t *source,
    const uint32_t *candidates, size_t count) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!source || !source->pos1 || !source->ref || !source->alt ||
        !source->ref_len || !source->alt_len || (count && !candidates) ||
        s->serial == UINT64_MAX || count > UINT32_MAX ||
        (uint32_t)source->ref_len - 1u > UINT32_MAX - source->pos1 ||
        count > UINT64_MAX - s->projected_events)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    if (s->have_input && (source->chrom_id < s->last_chrom ||
        (source->chrom_id == s->last_chrom && (source->pos1 < s->last_pos1 ||
         (source->pos1 == s->last_pos1 && source->event_id <= s->last_event_id)))))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER);
    duckvep_event_t prepared;
    if (!duckvep_event_prepare_small(source->pos1, source->ref, source->ref_len,
                                     source->alt, source->alt_len, &prepared))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    prepared.chrom_id = source->chrom_id;
    const duckvep_transcript_model_t *model = s->carriers.model;
    uint32_t last_end1 = 0u;
    for (size_t i = 0u; i < count; i++) {
        uint32_t tx = candidates[i];
        if (tx >= model->transcript_count || (i && tx <= candidates[i - 1u]) ||
            model->chrom_id[tx] != source->chrom_id || model->end1[tx] < source->pos1 ||
            (model->start1[tx] > prepared.raw_end1 &&
             model->start1[tx] > duckvep_event_feature_max1(&prepared)))
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
        if (model->end1[tx] > last_end1) last_end1 = model->end1[tx];
    }
    uint32_t completed;
    duckvep_carriers_status_t status = duckvep_carriers_advance(&s->carriers,
        source->chrom_id, source->pos1, s->serial + 1u, &completed);
    if (status == DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
        s->closing = 1u;
        return DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY;
    }
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    reclaim(s, source->chrom_id, source->pos1, 0);
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    if (count) {
        if (s->event_count == b->event_capacity)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL);
        if (count > b->projection_capacity - s->projection_count)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL);
        size_t bytes = (size_t)source->ref_len + source->alt_len;
        size_t at = ring_add(s->allele_begin, s->allele_count, b->allele_capacity);
        size_t padding = bytes > b->allele_capacity - at ? b->allele_capacity - at : 0u;
        size_t available = b->allele_capacity - s->allele_count;
        if (padding > available || bytes > available - padding)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL);
        if (padding) at = 0u;
        uint32_t event_at = (uint32_t)ring_add(s->event_begin, s->event_count, b->event_capacity);
        duckvep_haplotype_stored_event_t *stored = &b->events[event_at];
        *stored = (duckvep_haplotype_stored_event_t){0};
        stored->source = *source;
        stored->source.ref = b->alleles + at;
        stored->source.alt = b->alleles + at + source->ref_len;
        memcpy(b->alleles + at, source->ref, source->ref_len);
        memcpy(b->alleles + at + source->ref_len, source->alt, source->alt_len);
        stored->serial = s->serial + 1u;
        stored->allele_consumed = padding + bytes;
        stored->projection_begin = (uint32_t)ring_add(s->projection_begin,
            s->projection_count, b->projection_capacity);
        stored->projection_count = (uint32_t)count;
        stored->last_end1 = last_end1;
        duckvep_prepared_cds_allele_t allele = {&prepared,
            stored->source.ref + prepared.ref_diff_offset,
            stored->source.alt + prepared.alt_diff_offset,
            stored->source.ref + prepared.anchor_ref_offset,
            prepared.ref_diff_length, prepared.alt_diff_length, 1};
        for (size_t i = 0u; i < count; i++) {
            size_t projection_at = ring_add(stored->projection_begin, i, b->projection_capacity);
            duckvep_haplotype_projection_t *p = &b->projections[projection_at];
            p->transcript_index = candidates[i];
            memset(&p->edit, 0, sizeof(p->edit));
            p->status = duckvep_cds_edit_build_prepared_allele(model, s->exons, s->sequences,
                candidates[i], model->strand[candidates[i]], &allele, UINT32_MAX, &p->edit);
        }
        s->current_event = event_at;
        s->event_count++;
        s->projection_count += (uint32_t)count;
        s->allele_count += stored->allele_consumed;
        if (s->event_count > s->peak_events) s->peak_events = s->event_count;
        if (s->projection_count > s->peak_projections) s->peak_projections = s->projection_count;
        if (s->allele_count > s->peak_alleles) s->peak_alleles = s->allele_count;
    }
    s->have_current = count != 0u;
    s->have_input = 1u;
    s->serial++;
    s->last_event_id = source->event_id;
    s->last_pos1 = source->pos1;
    s->last_chrom = source->chrom_id;
    s->input_events++;
    s->projected_events += count;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push(
    duckvep_haplotype_stream_t *s, const duckvep_carrier_key_t *key, uint8_t evidence) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!s->have_input || s->closing || s->carriers.pending || s->carriers.finished || !key ||
        !key->lane || key->lane > key->ploidy || key->phase_set_present > 1u ||
        !evidence || (evidence & ~(DUCKVEP_CARRIER_CALLED | DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_UNPHASED)))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    if (!s->have_current) return DUCKVEP_HAPLOTYPE_STREAM_OK;
    const duckvep_haplotype_stored_event_t *event = &s->buffers.events[s->current_event];
    for (uint32_t i = 0u; i < event->projection_count; i++) {
        size_t at = ring_add(event->projection_begin, i, s->buffers.projection_capacity);
        duckvep_carriers_status_t status = duckvep_carriers_push(&s->carriers,
            s->buffers.projections[at].transcript_index, key, evidence);
        if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_finish(duckvep_haplotype_stream_t *s) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    uint32_t completed;
    duckvep_carriers_status_t status = duckvep_carriers_finish(&s->carriers, &completed);
    if (status == DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
        s->closing = 1u;
        return DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY;
    }
    if (status != DUCKVEP_CARRIERS_DONE) return carrier_fail(s, status);
    reclaim(s, 0u, 0u, 1);
    s->have_current = 0u;
    return DUCKVEP_HAPLOTYPE_STREAM_DONE;
}

static const duckvep_haplotype_stored_event_t *find_event(
    const duckvep_haplotype_stream_t *s, uint64_t serial) {
    size_t lo = 0u, hi = s->event_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        size_t at = ring_add(s->event_begin, mid, s->buffers.event_capacity);
        const duckvep_haplotype_stored_event_t *e = &s->buffers.events[at];
        if (e->serial == serial) return e;
        if (e->serial < serial) lo = mid + 1u;
        else hi = mid;
    }
    return NULL;
}

static const duckvep_haplotype_projection_t *find_projection(
    const duckvep_haplotype_stream_t *s, const duckvep_haplotype_stored_event_t *e, uint32_t tx) {
    size_t lo = 0u, hi = e->projection_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        size_t at = ring_add(e->projection_begin, mid, s->buffers.projection_capacity);
        const duckvep_haplotype_projection_t *p = &s->buffers.projections[at];
        if (p->transcript_index == tx) return p;
        if (p->transcript_index < tx) lo = mid + 1u;
        else hi = mid;
    }
    return NULL;
}

static int same_phase_set(duckvep_haplotype_phase_set_t a, duckvep_haplotype_phase_set_t b) {
    return a.present == b.present && (!a.present || a.value == b.value);
}

static duckvep_haplotype_stream_status_t push_call_lane(
    duckvep_haplotype_stream_t *s, uint32_t tx, const duckvep_haplotype_call_t *call,
    duckvep_haplotype_phase_set_t set, uint16_t lane, uint8_t evidence) {
    duckvep_carrier_key_t key = {call->sample_index, set.value, lane, call->ploidy, set.present};
    duckvep_carriers_status_t status = duckvep_carriers_push(&s->carriers, tx, &key, evidence);
    return status == DUCKVEP_CARRIERS_OK ? DUCKVEP_HAPLOTYPE_STREAM_OK : carrier_fail(s, status);
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push_call(
    duckvep_haplotype_stream_t *s, uint32_t tx, const duckvep_haplotype_call_t *call,
    const duckvep_haplotype_phase_set_t *sets, size_t set_count) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!s->have_current || s->closing || s->carriers.pending || s->carriers.finished ||
        !call || !call->alleles || !call->ploidy || !call->alt_index ||
        call->alt_index > INT32_MAX || call->phase_set.present > 1u ||
        (call->policy != DUCKVEP_PHASE_STRICT && call->policy != DUCKVEP_PHASE_VEP116_COMPAT) ||
        (s->have_phase_policy && s->phase_policy != call->policy) ||
        (set_count && !sets) || set_count > SIZE_MAX / sizeof(*sets) ||
        !find_projection(s, &s->buffers.events[s->current_event], tx))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    const duckvep_haplotype_phase_set_t absent = {0, 0u};
    if (!set_count) { sets = &absent; set_count = 1u; }
    size_t declared_set = set_count;
    for (size_t i = 0u; i < set_count; i++) {
        if (sets[i].present > 1u || (i && (!sets[i].present ||
            (sets[i - 1u].present && sets[i].value <= sets[i - 1u].value))))
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
        if (same_phase_set(sets[i], call->phase_set)) declared_set = i;
    }
    if (call->policy == DUCKVEP_PHASE_VEP116_COMPAT &&
        (set_count != 1u || sets[0].present))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);

    duckvep_phase_summary_t summary = {0};
    uint8_t missing = 0u, pool_missing = 0u, pool_alt = 0u;
    for (uint32_t slot = 0u; slot < call->ploidy; slot++) {
        int32_t allele = call->alleles[slot];
        uint8_t phase = call->phase_before ? call->phase_before[slot] : 0u;
        if (duckvep_phase_observe(&summary, allele, phase) != DUCKVEP_PHASE_OK)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
        if (allele < 0) missing = 1u;
        if (!phase) {
            if (allele < 0) pool_missing = 1u;
            if (allele == (int32_t)call->alt_index) pool_alt = 1u;
        }
    }
    int broadcast = summary.ploidy == 1u || summary.homozygous ||
        summary.unphased_count == summary.ploidy;
    if (call->policy == DUCKVEP_PHASE_STRICT && !broadcast && declared_set == set_count)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    s->phase_policy = call->policy;
    s->have_phase_policy = 1u;

    uint16_t called_before = 0u;
    for (uint32_t slot = 0u; slot < call->ploidy; slot++) {
        int32_t allele = call->alleles[slot];
        uint8_t phase = call->phase_before ? call->phase_before[slot] : 0u;
        duckvep_phase_assignment_t assignment;
        if (duckvep_phase_assign(&summary, (uint16_t)(slot + 1u), called_before,
            allele, phase, call->policy, &assignment) != DUCKVEP_PHASE_OK)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
        if (allele >= 0) called_before++;
        uint8_t evidence = 0u;
        if (call->policy == DUCKVEP_PHASE_VEP116_COMPAT) {
            if (allele < 0) continue; /* Missing slots do not consume compacted lanes. */
            if (missing) evidence = DUCKVEP_CARRIER_MISSING;
            if (allele == (int32_t)call->alt_index) evidence |= DUCKVEP_CARRIER_CALLED;
        } else if (assignment.scope == DUCKVEP_PHASE_UNRESOLVED) {
            if (pool_alt || pool_missing) evidence = DUCKVEP_CARRIER_UNPHASED;
            if (pool_missing) evidence |= DUCKVEP_CARRIER_MISSING;
            assignment.lane = (uint16_t)(slot + 1u);
        } else if (allele < 0) {
            evidence = DUCKVEP_CARRIER_MISSING;
        } else if (allele == (int32_t)call->alt_index) {
            evidence = DUCKVEP_CARRIER_CALLED;
        }
        if (!evidence) continue;
        size_t first = broadcast || call->policy == DUCKVEP_PHASE_VEP116_COMPAT ? 0u : declared_set;
        size_t end = broadcast ? set_count : first + 1u;
        for (size_t i = first; i < end; i++) {
            duckvep_haplotype_stream_status_t status = push_call_lane(s, tx, call, sets[i],
                assignment.lane, evidence);
            if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        }
    }
    /* Compaction loses the original missing slot positions, not their evidence.
     * Include the remaining lanes so an all-missing GT is never implicit REF. */
    if (call->policy == DUCKVEP_PHASE_VEP116_COMPAT && missing) {
        for (uint32_t lane = (uint32_t)called_before + 1u; lane <= call->ploidy; lane++) {
            duckvep_haplotype_stream_status_t status = push_call_lane(s, tx, call, absent,
                (uint16_t)lane, DUCKVEP_CARRIER_MISSING);
            if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        }
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

static void sift_edit_min(duckvep_haplotype_edit_t *edits, size_t root, size_t count) {
    duckvep_haplotype_edit_t value = edits[root];
    while (root < count / 2u) {
        size_t child = root * 2u + 1u;
        if (child + 1u < count && edits[child + 1u].cds_start < edits[child].cds_start) child++;
        if (value.cds_start <= edits[child].cds_start) break;
        edits[root] = edits[child];
        root = child;
    }
    edits[root] = value;
}

/* Genomic upload order is not CDS edit order: trimming retained REF can move
 * an earlier upload past the next edit, and reverse transcripts invert it.
 * A typed in-place heap keeps scratch constant; libc qsort may allocate. */
static void sort_edits_descending(duckvep_haplotype_edit_t *edits, size_t count) {
    for (size_t root = count / 2u; root > 0u; root--) sift_edit_min(edits, root - 1u, count);
    for (size_t remaining = count; remaining > 1u; remaining--) {
        duckvep_haplotype_edit_t first = edits[0];
        edits[0] = edits[remaining - 1u];
        edits[remaining - 1u] = first;
        sift_edit_min(edits, 0u, remaining - 1u);
    }
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_next(
    duckvep_haplotype_stream_t *s, duckvep_haplotype_leaf_t *out) {
    if (out) memset(out, 0, sizeof(*out));
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!out || !s->closing) return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    duckvep_haplotype_leaf_t leaf = {0};
    duckvep_carriers_status_t status = duckvep_carriers_next_leaf(&s->carriers, &leaf.carriers);
    if (status == DUCKVEP_CARRIERS_DONE) {
        status = duckvep_carriers_release(&s->carriers);
        if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
        s->closing = 0u;
        return DUCKVEP_HAPLOTYPE_STREAM_DONE;
    }
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    size_t count;
    status = duckvep_carriers_leaf_events(&s->carriers, leaf.carriers.id,
        b->leaf_events, b->leaf_capacity, &count);
    if (status == DUCKVEP_CARRIERS_OUTPUT_FULL)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL);
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    uint32_t tx = leaf.carriers.transcript_index;
    for (size_t i = 0u; i < count; i++) {
        const duckvep_haplotype_stored_event_t *e = find_event(s, b->leaf_events[i].event_id);
        const duckvep_haplotype_projection_t *p = e ? find_projection(s, e, tx) : NULL;
        if (!p) return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
        uint8_t evidence = b->leaf_events[i].evidence_flags;
        b->contributors[i] = (duckvep_haplotype_contributor_t){e->source, p->status, evidence};
        leaf.evidence_flags |= evidence;
        if (leaf.projection_status == DUCKVEP_CDS_EDIT_OK) leaf.projection_status = p->status;
        if (p->status == DUCKVEP_CDS_EDIT_OK && (evidence & DUCKVEP_CARRIER_CALLED)) {
            duckvep_edit_set_t edits;
            duckvep_cds_edit_status_t split = duckvep_projected_cds_edit_set_build(&p->edit,
                s->carriers.model->strand[tx], b->edits + leaf.edit_count,
                b->edit_capacity - leaf.edit_count, &edits);
            if (split == DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL)
                return fail(s, DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL);
            if (split == DUCKVEP_CDS_EDIT_OK) leaf.edit_count += edits.count;
            else if (leaf.projection_status == DUCKVEP_CDS_EDIT_OK) leaf.projection_status = split;
        }
    }
    leaf.contributors = b->contributors;
    leaf.contributor_count = count;
    if (leaf.evidence_flags & (DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_UNPHASED))
        leaf.sequence_status = DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE;
    if (leaf.projection_status == DUCKVEP_CDS_EDIT_OK && leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
        const duckvep_sequence_pool_t *seq = s->sequences;
        uint64_t offset = seq->cds_offset[tx];
        size_t length = seq->cds_length[tx];
        if (offset > seq->cds_bytes_len || length > seq->cds_bytes_len - offset)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
        sort_edits_descending(b->edits, leaf.edit_count);
        duckvep_haplotype_result_t applied, translated;
        leaf.sequence_status = duckvep_haplotype_apply_cds_edits(seq->cds_bytes + (size_t)offset,
            length, b->edits, leaf.edit_count, s->carriers.model->strand[tx], b->cds, b->cds_capacity,
            &leaf.cds_length, &applied);
        if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
            duckvep_codon_table_t table = seq->codon_table
                ? (duckvep_codon_table_t)seq->codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
            leaf.sequence_status = duckvep_haplotype_translate_cds(b->cds, leaf.cds_length,
                table, b->protein, b->protein_capacity, &leaf.protein_length, &translated);
            if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
                leaf.cds = b->cds;
                leaf.protein = b->protein;
                leaf.flags = applied.flags | translated.flags;
                if (leaf.cds_length > UINT64_MAX - s->translated_bases)
                    return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                s->translated_bases += leaf.cds_length;
            }
        }
        if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL);
        if (leaf.sequence_status != DUCKVEP_HAPLOTYPE_OK)
            leaf.cds_length = leaf.protein_length = 0u;
    }
    s->completed_leaves++;
    *out = leaf;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}
