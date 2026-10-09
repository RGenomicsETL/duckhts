/* SQL fragments shared by the macros that read allele counts from a caller
 * relation (duckhts_roh_counts, duckhts_count_error_fit). */
#ifndef DUCKHTS_COUNT_SQL_H
#define DUCKHTS_COUNT_SQL_H

/* A value is cast to `type` only when it is a whole number: the cast alone
 * would round 2.5 to 3. The test is on the value's own type, so a DECIMAL
 * with more digits than a DOUBLE holds is tested exactly; the value must be
 * numeric. */
#define DUCKHTS_WHOLE_AS(column, type, message) \
    "CASE WHEN " column " != trunc(" column ") THEN error('" message "') " \
    "ELSE CAST(" column " AS " type ") END"
#define DUCKHTS_WHOLE_COUNT_AS(column, type) \
    DUCKHTS_WHOLE_AS(column, type, "read counts must be whole numbers")
#define DUCKHTS_WHOLE_COUNT(column) DUCKHTS_WHOLE_COUNT_AS(column, "INTEGER")
#define DUCKHTS_WHOLE_OPTION(name) DUCKHTS_WHOLE_AS(name, "BIGINT", name " must be a whole number")

#endif
