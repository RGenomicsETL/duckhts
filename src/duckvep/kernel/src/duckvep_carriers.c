#include "duckvep_carriers.h"

#include <string.h>

static uint64_t mix(uint64_t x) {
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static uint64_t call_hash(uint32_t tx, const duckvep_carrier_key_t *key) {
    return mix(((uint64_t)tx << 32) | key->sample_index) ^
        mix((uint64_t)(key->phase_set_present ? key->phase_set : 0)) ^
        mix(((uint64_t)key->lane << 1) | key->phase_set_present);
}

static uint64_t prefix_hash(uint32_t tx, uint32_t parent, uint64_t event) {
    return mix(((uint64_t)tx << 32) | parent) ^ mix(event);
}

/* Backshift deletion keeps lookups bounded by live entries rather than a
 * growing history of tombstones after thousands of transcript closures. */
static void remove_bucket(duckvep_carrier_bucket_t *buckets, uint32_t count, uint32_t at) {
    uint32_t mask = count - 1u;
    uint32_t next = (at + 1u) & mask;
    while (buckets[next].id) {
        uint32_t home = (uint32_t)buckets[next].hash & mask;
        if (((next - home) & mask) >= ((next - at) & mask)) {
            buckets[at] = buckets[next];
            at = next;
        }
        next = (next + 1u) & mask;
    }
    buckets[at].id = 0u;
}

static uint32_t find_transcript(const duckvep_carriers_t *s, uint32_t tx, uint32_t *at) {
    uint32_t count = s->buffers.transcript_buckets;
    uint64_t hash = mix(tx);
    *at = count ? (uint32_t)hash & (count - 1u) : 0u;
    if (!count) return 0u;
    while (s->buffers.transcript_index[*at].id) {
        uint32_t id = s->buffers.transcript_index[*at].id;
        if (s->buffers.transcripts[id - 1u].transcript_index == tx) return id;
        *at = (*at + 1u) & (count - 1u);
    }
    return 0u;
}

static uint32_t find_call(const duckvep_carriers_t *s, uint32_t tx,
                         const duckvep_carrier_key_t *key, uint32_t *at) {
    uint32_t count = s->buffers.call_buckets;
    uint64_t hash = call_hash(tx, key);
    *at = count ? (uint32_t)hash & (count - 1u) : 0u;
    if (!count) return 0u;
    while (s->buffers.call_index[*at].id) {
        uint32_t id = s->buffers.call_index[*at].id;
        const duckvep_carrier_call_t *call = &s->buffers.calls[id - 1u];
        if (call->transcript == tx && call->key.sample_index == key->sample_index &&
            call->key.lane == key->lane && call->key.phase_set_present == key->phase_set_present &&
            (!key->phase_set_present || call->key.phase_set == key->phase_set)) return id;
        *at = (*at + 1u) & (count - 1u);
    }
    return 0u;
}

static uint32_t find_prefix(const duckvep_carriers_t *s, uint32_t tx, uint32_t parent,
                           uint64_t event, uint32_t *at) {
    uint32_t count = s->buffers.prefix_buckets;
    uint64_t hash = prefix_hash(tx, parent, event);
    *at = count ? (uint32_t)hash & (count - 1u) : 0u;
    if (!count) return 0u;
    while (s->buffers.prefix_index[*at].id) {
        uint32_t id = s->buffers.prefix_index[*at].id;
        const duckvep_carrier_prefix_t *prefix = &s->buffers.prefixes[id - 1u];
        if (prefix->transcript == tx && prefix->parent == parent && prefix->event_id == event) return id;
        *at = (*at + 1u) & (count - 1u);
    }
    return 0u;
}

static int valid_pool(const void *slots, uint32_t capacity, size_t size,
                      const duckvep_carrier_bucket_t *buckets, uint32_t bucket_count) {
    if (!capacity) return bucket_count == 0u;
    if (bucket_count && sizeof(*buckets) > SIZE_MAX / (size_t)bucket_count) return 0;
    return slots && buckets && capacity <= UINT32_MAX / 2u &&
        capacity <= SIZE_MAX / size &&
        bucket_count >= 2u * capacity && (bucket_count & (bucket_count - 1u)) == 0u;
}

/* The active array is an expiry min-heap. Input advances inspect only its
 * earliest transcript; opening/closing one path updates O(log active) slots. */
static int expires_before(const duckvep_carriers_t *s, uint32_t left, uint32_t right) {
    uint32_t a = s->buffers.transcripts[left - 1u].transcript_index;
    uint32_t b = s->buffers.transcripts[right - 1u].transcript_index;
    if (s->model->chrom_id[a] != s->model->chrom_id[b])
        return s->model->chrom_id[a] < s->model->chrom_id[b];
    if (s->model->end1[a] != s->model->end1[b]) return s->model->end1[a] < s->model->end1[b];
    return a < b;
}

static void open_transcript(duckvep_carriers_t *s, uint32_t id) {
    uint32_t at = s->transcript_count++;
    while (at) {
        uint32_t parent = (at - 1u) / 2u;
        uint32_t parent_id = s->buffers.active_transcripts[parent];
        if (!expires_before(s, id, parent_id)) break;
        s->buffers.active_transcripts[at] = parent_id;
        at = parent;
    }
    s->buffers.active_transcripts[at] = id;
}

static void close_transcript(duckvep_carriers_t *s) {
    uint32_t count = --s->transcript_count;
    uint32_t last = s->buffers.active_transcripts[count];
    s->buffers.active_transcripts[count] = 0u;
    if (!count) return;
    uint32_t at = 0u;
    while (2u * at + 1u < count) {
        uint32_t child = 2u * at + 1u;
        if (child + 1u < count && expires_before(s, s->buffers.active_transcripts[child + 1u],
                                                s->buffers.active_transcripts[child])) child++;
        uint32_t child_id = s->buffers.active_transcripts[child];
        if (!expires_before(s, child_id, last)) break;
        s->buffers.active_transcripts[at] = child_id;
        at = child;
    }
    s->buffers.active_transcripts[at] = last;
}

duckvep_carriers_status_t duckvep_carriers_init(
    duckvep_carriers_t *s, const duckvep_transcript_model_t *model,
    const duckvep_carrier_buffers_t *b) {
    if (!s) return DUCKVEP_CARRIERS_INVALID_ARG;
    memset(s, 0, sizeof(*s));
    if (!model || !b || model->transcript_count > UINT32_MAX ||
        (model->transcript_count && (!model->chrom_id || !model->end1)) ||
        (b->transcript_capacity && !b->active_transcripts) ||
        !valid_pool(b->transcripts, b->transcript_capacity, sizeof(*b->transcripts),
                    b->transcript_index, b->transcript_buckets) ||
        !valid_pool(b->calls, b->call_capacity, sizeof(*b->calls), b->call_index, b->call_buckets) ||
        !valid_pool(b->prefixes, b->prefix_capacity, sizeof(*b->prefixes),
                    b->prefix_index, b->prefix_buckets)) return DUCKVEP_CARRIERS_INVALID_ARG;
    s->model = model;
    s->buffers = *b;
    if (b->transcript_capacity) memset(b->active_transcripts, 0,
                                       (size_t)b->transcript_capacity * sizeof(*b->active_transcripts));
    for (uint32_t i = 0u; i < b->transcript_capacity; i++) {
        memset(&b->transcripts[i], 0, sizeof(b->transcripts[i]));
        b->transcripts[i].next_free = i + 1u < b->transcript_capacity ? i + 2u : 0u;
    }
    for (uint32_t i = 0u; i < b->call_capacity; i++) {
        memset(&b->calls[i], 0, sizeof(b->calls[i]));
        b->calls[i].next_free = i + 1u < b->call_capacity ? i + 2u : 0u;
    }
    for (uint32_t i = 0u; i < b->prefix_capacity; i++) {
        memset(&b->prefixes[i], 0, sizeof(b->prefixes[i]));
        b->prefixes[i].next_free = i + 1u < b->prefix_capacity ? i + 2u : 0u;
    }
    if (b->transcript_buckets) memset(b->transcript_index, 0,
                                      (size_t)b->transcript_buckets * sizeof(*b->transcript_index));
    if (b->call_buckets) memset(b->call_index, 0, (size_t)b->call_buckets * sizeof(*b->call_index));
    if (b->prefix_buckets) memset(b->prefix_index, 0,
                                 (size_t)b->prefix_buckets * sizeof(*b->prefix_index));
    s->free_transcript = b->transcript_capacity ? 1u : 0u;
    s->free_call = b->call_capacity ? 1u : 0u;
    s->free_prefix = b->prefix_capacity ? 1u : 0u;
    s->initialized = 1u;
    return DUCKVEP_CARRIERS_OK;
}

static duckvep_carriers_status_t advance(
    duckvep_carriers_t *s, uint16_t chrom, uint32_t pos1, uint64_t event,
    int eof, uint32_t *completed) {
    if (completed) *completed = UINT32_MAX;
    if (!s || !s->initialized || !completed || (!eof && !pos1)) return DUCKVEP_CARRIERS_INVALID_ARG;
    if (s->finished) return eof ? DUCKVEP_CARRIERS_DONE : DUCKVEP_CARRIERS_INPUT_ORDER;
    if (s->pending) {
        if (s->pending_eof != eof || (!eof && (chrom != s->pending_chrom ||
            pos1 != s->pending_pos1 || event != s->pending_event_id))) return DUCKVEP_CARRIERS_INPUT_ORDER;
    } else {
        if (!eof && s->have_event && (chrom < s->chrom ||
            (chrom == s->chrom && (pos1 < s->pos1 ||
             (pos1 == s->pos1 && event < s->event_id))))) return DUCKVEP_CARRIERS_INPUT_ORDER;
        s->pending = 1u;
        s->pending_eof = (uint8_t)eof;
        s->pending_chrom = chrom;
        s->pending_pos1 = pos1;
        s->pending_event_id = event;
    }
    if (s->closing) {
        *completed = s->buffers.transcripts[s->closing - 1u].transcript_index;
        return DUCKVEP_CARRIERS_TRANSCRIPT_READY;
    }
    if (s->transcript_count) {
        uint32_t id = s->buffers.active_transcripts[0];
        const duckvep_carrier_transcript_t *tx = &s->buffers.transcripts[id - 1u];
        uint32_t ordinal = tx->transcript_index;
        if (eof || s->model->chrom_id[ordinal] < chrom ||
            (s->model->chrom_id[ordinal] == chrom && s->model->end1[ordinal] < pos1)) {
            s->closing = id;
            s->next_leaf = tx->first_prefix;
            s->drained = 0u;
            *completed = ordinal;
            return DUCKVEP_CARRIERS_TRANSCRIPT_READY;
        }
    }
    s->pending = 0u;
    if (eof) {
        s->finished = 1u;
        return DUCKVEP_CARRIERS_DONE;
    }
    s->chrom = chrom;
    s->pos1 = pos1;
    s->event_id = event;
    s->have_event = 1u;
    return DUCKVEP_CARRIERS_OK;
}

duckvep_carriers_status_t duckvep_carriers_advance(
    duckvep_carriers_t *s, uint16_t chrom, uint32_t pos1, uint64_t event,
    uint32_t *completed) {
    return advance(s, chrom, pos1, event, 0, completed);
}

duckvep_carriers_status_t duckvep_carriers_finish(duckvep_carriers_t *s, uint32_t *completed) {
    return advance(s, 0u, 0u, 0u, 1, completed);
}

duckvep_carriers_status_t duckvep_carriers_push(
    duckvep_carriers_t *s, uint32_t tx_index, const duckvep_carrier_key_t *key) {
    uint32_t tx_at, call_at, prefix_at;
    if (!s || !s->initialized || !key || !s->have_event || s->finished || s->pending ||
        tx_index >= s->model->transcript_count || !key->lane || key->lane > key->ploidy ||
        key->phase_set_present > 1u || s->model->chrom_id[tx_index] != s->chrom ||
        s->model->end1[tx_index] < s->pos1) return DUCKVEP_CARRIERS_INVALID_ARG;
    uint32_t tx_id = find_transcript(s, tx_index, &tx_at);
    int new_tx = !tx_id;
    if (new_tx) {
        if (!s->free_transcript) return DUCKVEP_CARRIERS_TRANSCRIPT_FULL;
        tx_id = s->free_transcript;
    }
    uint32_t call_id = find_call(s, tx_id, key, &call_at);
    uint32_t parent = call_id ? s->buffers.calls[call_id - 1u].leaf : 0u;
    if (call_id && s->buffers.calls[call_id - 1u].key.ploidy != key->ploidy)
        return DUCKVEP_CARRIERS_PLOIDY_CONFLICT;
    if (parent && s->buffers.prefixes[parent - 1u].event_id == s->event_id)
        return DUCKVEP_CARRIERS_DUPLICATE_CALL;
    if (!call_id && !s->free_call) return DUCKVEP_CARRIERS_CALL_FULL;
    uint32_t prefix_id = find_prefix(s, tx_id, parent, s->event_id, &prefix_at);
    if (!prefix_id && !s->free_prefix) return DUCKVEP_CARRIERS_PREFIX_FULL;

    /* Every capacity and key check precedes publication of any new slot. */
    duckvep_carrier_transcript_t *tx = &s->buffers.transcripts[tx_id - 1u];
    if (new_tx) {
        s->free_transcript = tx->next_free;
        memset(tx, 0, sizeof(*tx));
        tx->transcript_index = tx_index;
        open_transcript(s, tx_id);
        s->buffers.transcript_index[tx_at] = (duckvep_carrier_bucket_t){mix(tx_index), tx_id};
    }
    if (!prefix_id) {
        prefix_id = s->free_prefix;
        duckvep_carrier_prefix_t *prefix = &s->buffers.prefixes[prefix_id - 1u];
        s->free_prefix = prefix->next_free;
        memset(prefix, 0, sizeof(*prefix));
        prefix->event_id = s->event_id;
        prefix->transcript = tx_id;
        prefix->parent = parent;
        prefix->depth = parent ? s->buffers.prefixes[parent - 1u].depth + 1u : 1u;
        prefix->next_transcript = tx->first_prefix;
        tx->first_prefix = prefix_id;
        s->buffers.prefix_index[prefix_at] = (duckvep_carrier_bucket_t){
            prefix_hash(tx_id, parent, s->event_id), prefix_id};
        s->prefix_count++;
    }
    if (!call_id) {
        call_id = s->free_call;
        duckvep_carrier_call_t *call = &s->buffers.calls[call_id - 1u];
        s->free_call = call->next_free;
        memset(call, 0, sizeof(*call));
        call->key = *key;
        if (!key->phase_set_present) call->key.phase_set = 0;
        call->transcript = tx_id;
        call->next_transcript = tx->first_call;
        tx->first_call = call_id;
        s->buffers.call_index[call_at] = (duckvep_carrier_bucket_t){call_hash(tx_id, key), call_id};
        s->call_count++;
    }
    duckvep_carrier_call_t *call = &s->buffers.calls[call_id - 1u];
    if (parent) {
        duckvep_carrier_prefix_t *old = &s->buffers.prefixes[parent - 1u];
        if (call->previous_leaf) s->buffers.calls[call->previous_leaf - 1u].next_leaf = call->next_leaf;
        else old->first_call = call->next_leaf;
        if (call->next_leaf) s->buffers.calls[call->next_leaf - 1u].previous_leaf = call->previous_leaf;
        old->call_count--;
    }
    duckvep_carrier_prefix_t *leaf = &s->buffers.prefixes[prefix_id - 1u];
    call->leaf = prefix_id;
    call->previous_leaf = 0u;
    call->next_leaf = leaf->first_call;
    if (leaf->first_call) s->buffers.calls[leaf->first_call - 1u].previous_leaf = call_id;
    leaf->first_call = call_id;
    leaf->call_count++;
    if (s->transcript_count > s->peak_transcripts) s->peak_transcripts = s->transcript_count;
    if (s->call_count > s->peak_calls) s->peak_calls = s->call_count;
    if (s->prefix_count > s->peak_prefixes) s->peak_prefixes = s->prefix_count;
    return DUCKVEP_CARRIERS_OK;
}

duckvep_carriers_status_t duckvep_carriers_next_leaf(
    duckvep_carriers_t *s, duckvep_carrier_leaf_t *leaf) {
    if (leaf) memset(leaf, 0, sizeof(*leaf));
    if (!s || !s->initialized || !s->closing || !leaf) return DUCKVEP_CARRIERS_INVALID_ARG;
    while (s->next_leaf) {
        uint32_t id = s->next_leaf;
        const duckvep_carrier_prefix_t *prefix = &s->buffers.prefixes[id - 1u];
        s->next_leaf = prefix->next_transcript;
        if (!prefix->call_count) continue;
        *leaf = (duckvep_carrier_leaf_t){id,
            s->buffers.transcripts[s->closing - 1u].transcript_index,
            prefix->depth, prefix->first_call, prefix->call_count};
        return DUCKVEP_CARRIERS_OK;
    }
    s->drained = 1u;
    return DUCKVEP_CARRIERS_DONE;
}

duckvep_carriers_status_t duckvep_carriers_leaf_events(
    const duckvep_carriers_t *s, uint32_t id, uint64_t *events,
    size_t capacity, size_t *required) {
    if (required) *required = 0u;
    if (!s || !s->initialized || !s->closing || !required || !id ||
        id > s->buffers.prefix_capacity || (capacity && !events)) return DUCKVEP_CARRIERS_INVALID_ARG;
    const duckvep_carrier_prefix_t *prefix = &s->buffers.prefixes[id - 1u];
    if (prefix->transcript != s->closing || !prefix->call_count) return DUCKVEP_CARRIERS_INVALID_ARG;
    *required = prefix->depth;
    if (capacity < *required) return DUCKVEP_CARRIERS_OUTPUT_FULL;
    size_t at = *required;
    while (id) {
        prefix = &s->buffers.prefixes[id - 1u];
        events[--at] = prefix->event_id;
        id = prefix->parent;
    }
    return DUCKVEP_CARRIERS_OK;
}

const duckvep_carrier_call_t *duckvep_carriers_call(const duckvep_carriers_t *s, uint32_t id) {
    if (!s || !s->initialized || !s->closing || !id || id > s->buffers.call_capacity) return NULL;
    const duckvep_carrier_call_t *call = &s->buffers.calls[id - 1u];
    return call->transcript == s->closing ? call : NULL;
}

duckvep_carriers_status_t duckvep_carriers_release(duckvep_carriers_t *s) {
    if (!s || !s->initialized || !s->closing || !s->drained) return DUCKVEP_CARRIERS_INVALID_ARG;
    uint32_t tx_id = s->closing;
    duckvep_carrier_transcript_t *tx = &s->buffers.transcripts[tx_id - 1u];
    uint32_t id = tx->first_call, at;
    while (id) {
        duckvep_carrier_call_t *call = &s->buffers.calls[id - 1u];
        uint32_t next = call->next_transcript;
        (void)find_call(s, tx_id, &call->key, &at);
        remove_bucket(s->buffers.call_index, s->buffers.call_buckets, at);
        call->transcript = 0u;
        call->next_free = s->free_call;
        s->free_call = id;
        s->call_count--;
        id = next;
    }
    id = tx->first_prefix;
    while (id) {
        duckvep_carrier_prefix_t *prefix = &s->buffers.prefixes[id - 1u];
        uint32_t next = prefix->next_transcript;
        (void)find_prefix(s, tx_id, prefix->parent, prefix->event_id, &at);
        remove_bucket(s->buffers.prefix_index, s->buffers.prefix_buckets, at);
        prefix->transcript = 0u;
        prefix->next_free = s->free_prefix;
        s->free_prefix = id;
        s->prefix_count--;
        id = next;
    }
    (void)find_transcript(s, tx->transcript_index, &at);
    remove_bucket(s->buffers.transcript_index, s->buffers.transcript_buckets, at);
    close_transcript(s);
    tx->next_free = s->free_transcript;
    s->free_transcript = tx_id;
    s->closing = 0u;
    s->drained = 0u;
    return DUCKVEP_CARRIERS_OK;
}
