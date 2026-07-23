/*
 * duckvep_variant_tile.c — see the header. Pure C: storage + run-discipline
 * validation, no DuckDB and no interner. Variant arrays are fixed-capacity
 * (FULL_CAPACITY when count == cap); the allele and variant_id byte pools grow so a
 * single row always fits. Run-validation state persists across reset_rows.
 */
#include "duckvep_variant_tile.h"
#include "kernel/src/duckvep_event.h"

#include <stdlib.h>
#include <string.h>

struct duckvep_variant_tile {
    /* kernel batch SoA, parallel, length == capacity (count rows used) */
    uint16_t *chrom_id;
    uint32_t *pos1;
    uint32_t *end1;
    uint32_t *ref_offset;
    uint16_t *ref_length;
    uint32_t *alt_offset;
    uint16_t *alt_length;
    uint8_t  *variant_kind;
    uint8_t  *sv_type;
    uint8_t  *copy_change;
    size_t    capacity;
    size_t    count;

    /* shared REF/ALT byte pool (growable) */
    uint8_t *allele_bytes;
    size_t   allele_len;
    size_t   allele_cap;

    /* variant_id metadata pool (growable), parallel offset/length arrays sized to capacity */
    uint8_t *id_bytes;
    size_t  id_len;
    size_t  id_cap;
    size_t *id_offset; /* [capacity] */
    size_t *id_length; /* [capacity] */

    /* run-validation state (persists across reset_rows) */
    int      run_active;  /* 1 once a run is open */
    uint16_t run_chrom;   /* current open contig id */
    uint32_t last_pos1;   /* last pos1 seen in the run */
    uint8_t *closed;      /* [max_chroms] bitset-as-bytes: 1 = run for this id has closed */
    size_t   max_chroms;

    char *error; /* owned, last hard validation error */
};

static void tile_set_error(duckvep_variant_tile_t *t, const char *msg) {
    size_t n;
    char  *p;
    if (t == NULL) return;
    free(t->error);
    t->error = NULL;
    if (msg == NULL) return;
    n = strlen(msg);
    p = (char *)malloc(n + 1u);
    if (p != NULL) {
        memcpy(p, msg, n + 1u);
        t->error = p;
    }
}

/* Ensure a growable byte pool has room for `need` more bytes. Returns 1 on success. */
static int pool_reserve(uint8_t **buf, size_t *cap, size_t used, size_t need) {
    size_t want;
    size_t newcap;
    uint8_t *nb;
    if (need > SIZE_MAX - used) return 0; /* size_t overflow */
    want = used + need;
    if (want <= *cap) return 1;
    newcap = (*cap != 0u) ? *cap : 64u;
    while (newcap < want) {
        if (newcap > (SIZE_MAX / 2u)) {
            newcap = want;
            break;
        }
        newcap *= 2u;
    }
    nb = (uint8_t *)realloc(*buf, newcap);
    if (nb == NULL) return 0;
    *buf = nb;
    *cap = newcap;
    return 1;
}

duckvep_variant_tile_t *duckvep_variant_tile_create(size_t variant_capacity,
                                                    size_t max_chroms) {
    duckvep_variant_tile_t *t;
    if (variant_capacity == 0u || max_chroms == 0u) return NULL;
    if (variant_capacity > (size_t)UINT32_MAX) return NULL; /* kernel variant_idx is uint32_t */
    t = (duckvep_variant_tile_t *)calloc(1u, sizeof *t);
    if (t == NULL) return NULL;
    t->capacity = variant_capacity;
    t->max_chroms = max_chroms;

    t->chrom_id = (uint16_t *)calloc(variant_capacity, sizeof *t->chrom_id);
    t->pos1 = (uint32_t *)calloc(variant_capacity, sizeof *t->pos1);
    t->end1 = (uint32_t *)calloc(variant_capacity, sizeof *t->end1);
    t->ref_offset = (uint32_t *)calloc(variant_capacity, sizeof *t->ref_offset);
    t->ref_length = (uint16_t *)calloc(variant_capacity, sizeof *t->ref_length);
    t->alt_offset = (uint32_t *)calloc(variant_capacity, sizeof *t->alt_offset);
    t->alt_length = (uint16_t *)calloc(variant_capacity, sizeof *t->alt_length);
    t->variant_kind = (uint8_t *)calloc(variant_capacity, sizeof *t->variant_kind);
    t->sv_type = (uint8_t *)calloc(variant_capacity, sizeof *t->sv_type);
    t->copy_change = (uint8_t *)calloc(variant_capacity, sizeof *t->copy_change);
    t->id_offset = (size_t *)calloc(variant_capacity, sizeof *t->id_offset);
    t->id_length = (size_t *)calloc(variant_capacity, sizeof *t->id_length);
    t->closed = (uint8_t *)calloc(max_chroms, sizeof *t->closed);

    if (t->chrom_id == NULL || t->pos1 == NULL || t->end1 == NULL ||
        t->ref_offset == NULL || t->ref_length == NULL || t->alt_offset == NULL ||
        t->alt_length == NULL || t->variant_kind == NULL || t->sv_type == NULL ||
        t->copy_change == NULL || t->id_offset == NULL || t->id_length == NULL ||
        t->closed == NULL) {
        duckvep_variant_tile_destroy(t);
        return NULL;
    }
    return t;
}

void duckvep_variant_tile_destroy(duckvep_variant_tile_t *t) {
    if (t == NULL) return;
    free(t->chrom_id);
    free(t->pos1);
    free(t->end1);
    free(t->ref_offset);
    free(t->ref_length);
    free(t->alt_offset);
    free(t->alt_length);
    free(t->variant_kind);
    free(t->sv_type);
    free(t->copy_change);
    free(t->allele_bytes);
    free(t->id_bytes);
    free(t->id_offset);
    free(t->id_length);
    free(t->closed);
    free(t->error);
    free(t);
}

duckvep_tile_status_t duckvep_variant_tile_append_event(
    duckvep_variant_tile_t *t, uint16_t chrom_id, uint32_t pos1, uint32_t end1,
    uint8_t kind, uint8_t sv_type, uint8_t copy_change,
    const char *ref, size_t ref_len, const char *alt, size_t alt_len,
    const char *variant_id, size_t variant_id_len) {
    duckvep_event_t event;
    size_t i;

    if (t == NULL) return DUCKVEP_TILE_INVALID;
    tile_set_error(t, NULL); /* error reflects only THIS append's outcome */

    /* 1. hard field validation: small VCF alleles or a single-interval structural event. */
    if (ref == NULL || alt == NULL || variant_id == NULL) {
        tile_set_error(t, "null ref/alt/variant_id");
        return DUCKVEP_TILE_INVALID;
    }
    if (ref_len > (size_t)UINT16_MAX || alt_len > (size_t)UINT16_MAX) {
        tile_set_error(t, "REF/ALT length exceeds uint16 storage");
        return DUCKVEP_TILE_INVALID;
    }
    if (kind <= (uint8_t)DUCKVEP_KIND_MNV) {
        if (ref_len == 0u || alt_len == 0u ||
            !duckvep_event_prepare_small(pos1, (const uint8_t *)ref,
                                         (uint16_t)ref_len,
                                         (const uint8_t *)alt,
                                         (uint16_t)alt_len, &event) ||
            event.raw_end1 != end1 || event.kind != kind) {
            tile_set_error(t,
                "small-variant kind or raw span is inconsistent with REF/ALT");
            return DUCKVEP_TILE_INVALID;
        }
    } else if (kind == (uint8_t)DUCKVEP_KIND_UNSPECIFIED_ALT) {
        if (!duckvep_event_unspecified_alt_shape_matches(
                pos1, end1, (const uint8_t *)ref, (uint16_t)ref_len,
                (const uint8_t *)alt, (uint16_t)alt_len)) {
            tile_set_error(t,
                "unspecified alternate must be an exact <*> allele with a valid REF span");
            return DUCKVEP_TILE_INVALID;
        }
    } else if (kind == (uint8_t)DUCKVEP_KIND_SV) {
        if (!duckvep_sv_geometry_valid((duckvep_sv_type_t)sv_type,
                                       pos1, end1)) {
            tile_set_error(t, "SV requires a valid 1-based inclusive interval");
            return DUCKVEP_TILE_INVALID;
        }
        if (!duckvep_sv_metadata_valid((duckvep_sv_type_t)sv_type,
                                        (duckvep_copy_change_t)copy_change)) {
            tile_set_error(t, "invalid or contradictory SV subtype/copy-change values");
            return DUCKVEP_TILE_INVALID;
        }
        if (sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
            tile_set_error(t, "BND requires paired two-locus coordinates");
            return DUCKVEP_TILE_INVALID;
        }
    } else {
        tile_set_error(t, "unknown variant kind");
        return DUCKVEP_TILE_INVALID;
    }
    if (kind != (uint8_t)DUCKVEP_KIND_SV &&
        (sv_type != (uint8_t)DUCKVEP_SV_NONE ||
         copy_change != (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN)) {
        tile_set_error(t, "structural metadata supplied for a non-SV row");
        return DUCKVEP_TILE_INVALID;
    }
    if (chrom_id >= t->max_chroms) {
        tile_set_error(t, "chrom_id out of interned id space");
        return DUCKVEP_TILE_INVALID;
    }

    /* 2. run discipline (contiguous chrom runs; nondecreasing pos1; no reappearance) */
    if (!t->run_active) {
        if (t->closed[chrom_id]) {
            tile_set_error(t, "contig reappeared after its run closed");
            return DUCKVEP_TILE_INVALID;
        }
        /* open the run with this chrom; fall through to capacity + store */
    } else if (chrom_id == t->run_chrom) {
        if (pos1 < t->last_pos1) {
            tile_set_error(t, "pos1 decreased within a contig run (input not sorted)");
            return DUCKVEP_TILE_INVALID;
        }
    } else { /* chrom changed while a run is open */
        if (t->closed[chrom_id]) {
            tile_set_error(t, "contig reappeared after its run closed");
            return DUCKVEP_TILE_INVALID;
        }
        return DUCKVEP_TILE_FULL_CHROM; /* flush + close_run, then retry this row */
    }

    /* 3. capacity (variant count only; byte pools grow) */
    if (t->count == t->capacity) {
        return DUCKVEP_TILE_FULL_CAPACITY; /* flush + reset_rows, then retry this row */
    }

    /* 4. store. Allele offsets are uint32_t into allele_bytes — guard the offset
     * space before reserving (ref_offset/alt_offset would otherwise wrap). */
    if (t->allele_len > (size_t)UINT32_MAX - (ref_len + alt_len)) {
        tile_set_error(t, "allele pool exceeds uint32 offset space");
        return DUCKVEP_TILE_INVALID;
    }
    /* Allele pool: ref then alt, contiguous. */
    if (!pool_reserve(&t->allele_bytes, &t->allele_cap, t->allele_len, ref_len + alt_len)) {
        tile_set_error(t, "allele pool OOM");
        return DUCKVEP_TILE_INVALID;
    }
    if (!pool_reserve(&t->id_bytes, &t->id_cap, t->id_len, variant_id_len)) {
        tile_set_error(t, "variant_id pool OOM");
        return DUCKVEP_TILE_INVALID;
    }

    i = t->count;
    t->chrom_id[i] = chrom_id;
    t->pos1[i] = pos1;
    t->end1[i] = end1;
    t->variant_kind[i] = kind;
    t->sv_type[i] = sv_type;
    t->copy_change[i] = copy_change;

    t->ref_offset[i] = (uint32_t)t->allele_len;
    t->ref_length[i] = (uint16_t)ref_len;
    if (ref_len != 0u) {
        memcpy(t->allele_bytes + t->allele_len, ref, ref_len);
        t->allele_len += ref_len;
    }

    t->alt_offset[i] = (uint32_t)t->allele_len;
    t->alt_length[i] = (uint16_t)alt_len;
    if (alt_len != 0u) {
        memcpy(t->allele_bytes + t->allele_len, alt, alt_len);
        t->allele_len += alt_len;
    }

    t->id_offset[i] = t->id_len;
    t->id_length[i] = variant_id_len;
    if (variant_id_len != 0u) {
        memcpy(t->id_bytes + t->id_len, variant_id, variant_id_len);
        t->id_len += variant_id_len;
    }

    t->run_active = 1;
    t->run_chrom = chrom_id;
    t->last_pos1 = pos1;
    t->count = i + 1u;
    return DUCKVEP_TILE_APPENDED;
}

duckvep_tile_status_t duckvep_variant_tile_append(
    duckvep_variant_tile_t *t, uint16_t chrom_id, uint32_t pos1, uint32_t end1,
    uint8_t kind, const char *ref, size_t ref_len, const char *alt, size_t alt_len,
    const char *variant_id, size_t variant_id_len) {
    return duckvep_variant_tile_append_event(
        t, chrom_id, pos1, end1, kind,
        kind == (uint8_t)DUCKVEP_KIND_SV ? (uint8_t)DUCKVEP_SV_UNKNOWN
                                         : (uint8_t)DUCKVEP_SV_NONE,
        (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN,
        ref, ref_len, alt, alt_len, variant_id, variant_id_len);
}

void duckvep_variant_tile_reset_rows(duckvep_variant_tile_t *t) {
    if (t == NULL) return;
    t->count = 0u;
    t->allele_len = 0u;
    t->id_len = 0u;
    /* run_active / run_chrom / last_pos1 / closed preserved (SAME chrom continues) */
}

void duckvep_variant_tile_close_run(duckvep_variant_tile_t *t) {
    if (t == NULL) return;
    if (t->run_active && t->run_chrom < t->max_chroms) {
        t->closed[t->run_chrom] = 1u;
    }
    t->run_active = 0;
    duckvep_variant_tile_reset_rows(t);
}

void duckvep_variant_tile_batch(const duckvep_variant_tile_t *t,
                                duckvep_variant_batch_t *out) {
    if (out == NULL) return;
    memset(out, 0, sizeof *out);
    if (t == NULL) return;
    out->chrom_id = t->chrom_id;
    out->pos1 = t->pos1;
    out->end1 = t->end1;
    out->ref_offset = t->ref_offset;
    out->ref_length = t->ref_length;
    out->alt_offset = t->alt_offset;
    out->alt_length = t->alt_length;
    out->allele_bytes = t->allele_bytes;
    out->allele_bytes_len = t->allele_len;
    out->variant_kind = t->variant_kind;
    out->sv_type = t->sv_type;
    out->copy_change = t->copy_change;
    out->count = t->count;
}

size_t duckvep_variant_tile_count(const duckvep_variant_tile_t *t) {
    return t != NULL ? t->count : 0u;
}

const char *duckvep_variant_tile_variant_id(const duckvep_variant_tile_t *t,
                                            size_t idx, size_t *len) {
    if (t == NULL || idx >= t->count) {
        if (len != NULL) *len = 0u;
        return NULL;
    }
    if (len != NULL) *len = t->id_length[idx];
    if (t->id_bytes == NULL) return ""; /* all ids zero-length: avoid NULL+offset UB */
    return (const char *)(t->id_bytes + t->id_offset[idx]);
}

const char *duckvep_variant_tile_error(const duckvep_variant_tile_t *t) {
    return t != NULL ? t->error : NULL;
}
