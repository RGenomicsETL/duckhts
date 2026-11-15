/* No DuckDB or HTSlib: check the BED line helpers against libc.
 *
 * duckhts_bed_parse_int64_span must agree with strtoll() in the C locale on a
 * NUL-terminated copy of the span, where a span is valid when strtoll()
 * consumes at least one digit and stops at the terminator.
 * duckhts_bed_line_split must agree with duckhts_bed_line_field and with a
 * strchr() walk of the line. */
#include "bed_line.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned failures;
static unsigned long checks;

#define CHECK(condition) do { \
    checks++; \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
        failures++; \
    } \
} while (0)

static uint64_t random_state = UINT64_C(0x9e3779b97f4a7c15);
static uint64_t next_random(void) {
    random_state ^= random_state << 13;
    random_state ^= random_state >> 7;
    random_state ^= random_state << 17;
    return random_state;
}

/* The libc oracle. Returns 1 and sets *out when strtoll() reads the whole
 * span as one integer. */
static int strtoll_span(const char *s, int len, int64_t *out) {
    char *copy = (char *)malloc((size_t)len + 1);
    char *end = NULL;
    const char *digits;
    long long value;
    int ok;

    if (!copy) abort();
    memcpy(copy, s, (size_t)len);
    copy[len] = '\0';
    value = strtoll(copy, &end, 10);
    /* strtoll() reports "no conversion" by leaving end at the start; an
     * embedded NUL ends its input before the span does. */
    digits = end;
    ok = end != copy && end == copy + len && digits[-1] >= '0' && digits[-1] <= '9';
    if (ok) *out = (int64_t)value;
    free(copy);
    return ok;
}

static void check_parse(const char *s, int len) {
    int64_t expected = 0, actual = 0;
    int expected_ok = strtoll_span(s, len, &expected);
    int actual_ok = duckhts_bed_parse_int64_span(s, len, &actual);
    checks++;
    if (expected_ok != actual_ok || (expected_ok && expected != actual)) {
        fprintf(stderr, "parse mismatch on span of %d bytes \"%.*s\": "
                "strtoll %d %" PRId64 ", helper %d %" PRId64 "\n",
                len, len, s, expected_ok, expected, actual_ok, actual);
        failures++;
    }
}

static void parse_examples(void) {
    static const char *const examples[] = {
        "0", "1", "9", "10", "007", "-0", "+0", "-1", "+1", "123456789",
        "9223372036854775806", "9223372036854775807", "9223372036854775808",
        "9223372036854775809", "9223372036854775817", "18446744073709551615",
        "18446744073709551616", "92233720368547758070",
        "-9223372036854775807", "-9223372036854775808", "-9223372036854775809",
        "-18446744073709551616", "99999999999999999999999999999999999999",
        "-99999999999999999999999999999999999999",
        " 1", "\t1", "\n1", "\v1", "\f1", "\r1", "  \t -5", " +7",
        "", " ", "+", "-", "+-1", "-+1", "--1", "++1", " ", "\t",
        "1 ", "1\t", "1x", "x1", "1.0", "1e3", "0x10", "1,000", "1 2",
        "- 1", "+ 1", "1-", "1+", "\xa0" "1", "\x1f" "1", "\x80", "12\xff",
        "/", ":", "0/", "0:",
    };
    int64_t out = 77;
    for (size_t i = 0; i < sizeof(examples) / sizeof(*examples); i++) {
        check_parse(examples[i], (int)strlen(examples[i]));
    }
    /* The span, not a terminator, ends the number. */
    CHECK(duckhts_bed_parse_int64_span("12345", 3, &out) && out == 123);
    CHECK(duckhts_bed_parse_int64_span("42\t99", 2, &out) && out == 42);
    out = 77;
    CHECK(!duckhts_bed_parse_int64_span(NULL, 1, &out) && out == 77);
    CHECK(!duckhts_bed_parse_int64_span("1", 0, &out) && out == 77);
    CHECK(!duckhts_bed_parse_int64_span("1", -1, &out) && out == 77);
    CHECK(!duckhts_bed_parse_int64_span("1", 1, NULL));
    CHECK(!duckhts_bed_parse_int64_span("1x", 2, &out) && out == 77);
    /* A long span is read in place: 4096 digits saturate, a bad last byte fails. */
    {
        enum { LONG_DIGITS = 4096 };
        char digits[LONG_DIGITS + 1];
        memset(digits, '7', LONG_DIGITS);
        digits[LONG_DIGITS] = '\0';
        check_parse(digits, LONG_DIGITS);
        CHECK(duckhts_bed_parse_int64_span(digits, LONG_DIGITS, &out) && out == INT64_MAX);
        digits[0] = '-';
        check_parse(digits, LONG_DIGITS);
        CHECK(duckhts_bed_parse_int64_span(digits, LONG_DIGITS, &out) && out == INT64_MIN);
        memset(digits, '0', LONG_DIGITS);
        digits[LONG_DIGITS - 1] = '5';
        check_parse(digits, LONG_DIGITS);
        CHECK(duckhts_bed_parse_int64_span(digits, LONG_DIGITS, &out) && out == 5);
        digits[LONG_DIGITS - 1] = 'x';
        check_parse(digits, LONG_DIGITS);
    }
}

/* Every span of up to three bytes over every byte value except NUL, which
 * the line helpers never pass inside a field. */
static void parse_short_spans(void) {
    char span[3];
    for (int a = 1; a < 256; a++) {
        span[0] = (char)a;
        check_parse(span, 1);
        for (int b = 1; b < 256; b++) {
            span[1] = (char)b;
            check_parse(span, 2);
        }
    }
    static const char alphabet[] = "0159+- \t\n\v\f\rx.";
    const int n = (int)sizeof(alphabet) - 1;
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            for (int c = 0; c < n; c++) {
                span[0] = alphabet[a];
                span[1] = alphabet[b];
                span[2] = alphabet[c];
                check_parse(span, 3);
            }
        }
    }
}

static void parse_random_spans(void) {
    static const char alphabet[] = "0123456789012345678901234567890123456789+- \tx";
    const int n = (int)sizeof(alphabet) - 1;
    char span[48];
    for (int trial = 0; trial < 400000; trial++) {
        int len = 1 + (int)(next_random() % 40);
        for (int i = 0; i < len; i++) span[i] = alphabet[next_random() % (uint64_t)n];
        check_parse(span, len);
    }
    /* Values around the int64 limits, with and without a sign. */
    for (int trial = 0; trial < 200000; trial++) {
        uint64_t magnitude = (uint64_t)INT64_MAX - 1000 + next_random() % 2001;
        int len = snprintf(span, sizeof(span), "%s%" PRIu64,
                           next_random() % 2 ? "-" : "", magnitude);
        check_parse(span, len);
    }
    for (int trial = 0; trial < 200000; trial++) {
        int len = snprintf(span, sizeof(span), "%" PRId64, (int64_t)next_random());
        check_parse(span, len);
    }
}

/* Splits `line` with strchr() and compares with the helper for every
 * max_fields from 1 to 16, and with duckhts_bed_line_field per field. */
static void check_split(const char *line) {
    enum { MAX_FIELDS = 16 };
    const char *expected_starts[MAX_FIELDS];
    int expected_lens[MAX_FIELDS];
    int total = 0;
    const char *start = line;
    for (;;) {
        const char *tab = strchr(start, '\t');
        if (total < MAX_FIELDS) {
            expected_starts[total] = start;
            expected_lens[total] = tab ? (int)(tab - start) : (int)strlen(start);
        }
        total++;
        if (!tab) break;
        start = tab + 1;
    }
    for (int max_fields = 1; max_fields <= MAX_FIELDS; max_fields++) {
        const char *starts[MAX_FIELDS];
        int lens[MAX_FIELDS];
        int expected_count = total < max_fields ? total : max_fields;
        int count = duckhts_bed_line_split(line, max_fields, starts, lens);
        CHECK(count == expected_count);
        for (int i = 0; i < count && i < expected_count; i++) {
            int field_len = -1;
            const char *field = duckhts_bed_line_field(line, i, &field_len);
            CHECK(starts[i] == expected_starts[i]);
            CHECK(lens[i] == expected_lens[i]);
            CHECK(field == starts[i] && field_len == lens[i]);
        }
    }
}

static void split_examples(void) {
    static const char *const lines[] = {
        "", "chr1", "\t", "\t\t", "chr1\t", "\tchr1", "chr1\t10\t20",
        "chr1\t10\t20\t", "chr1\t\t20", "chr1 10 20", "chr1\t10\t20\tname\t0\t+",
        "chr1\t10\t20\tn\t0\t+\t10\t20\t0,0,0\t1\t10,\t0,\textra1\textra2\textra3",
        "a\tb\tc\td\te\tf\tg\th\ti\tj\tk\tl\tm\tn\to\tp\tq\tr",
    };
    const char *starts[2] = {NULL, NULL};
    int lens[2] = {-1, -1};
    for (size_t i = 0; i < sizeof(lines) / sizeof(*lines); i++) check_split(lines[i]);
    CHECK(duckhts_bed_line_split(NULL, 2, starts, lens) == 0);
    CHECK(duckhts_bed_line_split("a\tb", 0, starts, lens) == 0);
    CHECK(duckhts_bed_line_split("a\tb", -1, starts, lens) == 0);
    CHECK(starts[0] == NULL && lens[0] == -1);
    /* The last stored field ends at its own tab even when more fields follow. */
    CHECK(duckhts_bed_line_split("ab\tcde\tf", 2, starts, lens) == 2);
    CHECK(lens[0] == 2 && lens[1] == 3 && starts[1][0] == 'c');
}

static void split_random_lines(void) {
    static const char alphabet[] = "\t\t\tchr0123456789 ,+-";
    const int n = (int)sizeof(alphabet) - 1;
    char line[96];
    for (int trial = 0; trial < 20000; trial++) {
        int len = (int)(next_random() % 90);
        for (int i = 0; i < len; i++) line[i] = alphabet[next_random() % (uint64_t)n];
        line[len] = '\0';
        check_split(line);
    }
}

int main(void) {
    parse_examples();
    parse_short_spans();
    parse_random_spans();
    split_examples();
    split_random_lines();
    if (failures) {
        fprintf(stderr, "bed line: %u failures in %lu checks\n", failures, checks);
        return 1;
    }
    printf("bed line: %lu checks passed\n", checks);
    return 0;
}
