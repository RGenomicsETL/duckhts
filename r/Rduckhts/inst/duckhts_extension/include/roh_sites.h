/* Packed per-site evidence for one runs-of-homozygosity decode.
 *
 * __duckhts_roh_sites() collects the sites of one sample and chromosome into
 * these records and returns them, sorted by position, as one BLOB: a header
 * followed by records of the header's kind. __duckhts_roh_decode() reads that
 * BLOB. The layout is private to those two functions and is never persisted.
 *
 * This header has no DuckDB dependency.
 */
#ifndef DUCKHTS_ROH_SITES_H
#define DUCKHTS_ROH_SITES_H

#include <stdint.h>

/* Genotype evidence of a site list. Zero is not a kind. */
typedef enum {
    DUCKHTS_ROH_SITES_PL = 1,    /* phred-scaled genotype likelihoods */
    DUCKHTS_ROH_SITES_GT = 2,    /* a called dosage 0, 1 or 2 */
    DUCKHTS_ROH_SITES_COUNTS = 3 /* allele read counts */
} duckhts_roh_sites_kind_t;

#define DUCKHTS_ROH_SITES_MAGIC "DHRS"
#define DUCKHTS_ROH_SITES_MAGIC_LENGTH 4

typedef struct {
    char magic[DUCKHTS_ROH_SITES_MAGIC_LENGTH];
    uint8_t kind; /* duckhts_roh_sites_kind_t */
    uint8_t reserved[3];
} duckhts_roh_sites_header_t;

/* A site with PL or GT evidence. */
typedef struct {
    double af;           /* frequency of the alternate allele; NaN when the site has none */
    uint32_t pos1;       /* one-based position */
    uint8_t evidence[3]; /* PL: RR, RA and AA, each capped at 255. GT: the dosage in [0]. */
    uint8_t usable;      /* 0 when the genotype evidence is missing or unusable */
} duckhts_roh_site_t;

/* A site with read counts. */
typedef struct {
    double af;       /* frequency of the counted allele; NaN when the site has none */
    uint32_t pos1;   /* one-based position */
    int32_t other;   /* reads of the other allele */
    int32_t counted; /* reads of the allele af refers to */
    uint32_t usable; /* 0 when either count is missing */
} duckhts_roh_count_site_t;

_Static_assert(sizeof(duckhts_roh_sites_header_t) == 8, "ROH site header is 8 bytes");
_Static_assert(sizeof(duckhts_roh_site_t) == 16, "a PL or GT site is 16 bytes");
_Static_assert(sizeof(duckhts_roh_count_site_t) == 24, "a read-count site is 24 bytes");

/* A BLOB length is 32 bits, so one list holds fewer than 2^32 bytes. */
#define DUCKHTS_ROH_MAX_SITES 100000000

/* Record size of a kind, or 0 for an unknown kind. */
static inline uint64_t duckhts_roh_site_size(uint32_t kind) {
    if (kind == DUCKHTS_ROH_SITES_PL || kind == DUCKHTS_ROH_SITES_GT) return sizeof(duckhts_roh_site_t);
    if (kind == DUCKHTS_ROH_SITES_COUNTS) return sizeof(duckhts_roh_count_site_t);
    return 0;
}

/* Checks the two memory limits of a decode. Returns NULL when they are valid,
 * otherwise the error text. has_* is 0 for a SQL NULL. */
static inline const char *duckhts_roh_check_limits(int has_max_sites, int64_t max_sites,
                                                   int has_max_site_bytes, int64_t max_site_bytes) {
    if (!has_max_sites || max_sites < 1 || max_sites > DUCKHTS_ROH_MAX_SITES) {
        return "max_sites must be between 1 and 100000000";
    }
    if (!has_max_site_bytes || max_site_bytes < 1) return "max_site_bytes must be at least 1";
    return (const char *)0;
}

#endif
