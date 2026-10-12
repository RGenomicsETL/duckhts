#ifndef DUCKHTS_BED_LINE_H
#define DUCKHTS_BED_LINE_H

#include <stdint.h>

/* Line-level BED helpers shared by read_bed, fasta_nuc(bed_path := ...),
 * duckhts_mosdepth(by := <BED>) and duckhts_bam_bed_coverage.
 *
 * Input is one line from hts_getline() or tbx_itr_next(). htslib opens plain,
 * gzip and BGZF text alike and strips the trailing "\n" and one "\r" before
 * these helpers see the line (third_party/htslib/kstring.c kgetline2, line 336;
 * third_party/htslib/bgzf.c bgzf_getline, line 2351). The helpers read the line
 * as a NUL-terminated C string, so an embedded NUL ends the line.
 *
 * Callers keep their own policy: which fields they require, coordinate checks,
 * error messages, and what to do with a malformed line. */

/* Returns 1 for a line that carries no interval: NULL, empty, a "#" comment,
 * or a line that starts with "track" or "browser" (prefix match, as UCSC BED
 * headers are written). Returns 0 otherwise. */
int duckhts_bed_line_is_meta(const char *line);

/* Returns a borrowed view of the 0-based tab-delimited field `index` of `line`
 * and stores its byte length in *len. A field may be empty (len 0). Returns
 * NULL with *len = 0 when the line has fewer than index + 1 fields. */
const char *duckhts_bed_line_field(const char *line, int index, int *len);

/* Parses the whole span [s, s + len) as a base-10 integer with strtoll() rules:
 * leading whitespace and a sign are accepted, trailing bytes are not. Overflow
 * saturates at INT64_MIN or INT64_MAX as strtoll() does; ERANGE is not checked.
 * The span has no length limit. Returns 1 and sets *out on success; returns 0
 * for a NULL or empty span, trailing bytes, or allocation failure. */
int duckhts_bed_parse_int64_span(const char *s, int len, int64_t *out);

#endif
