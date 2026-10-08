/* SQL fragments shared by the macros that read allele counts from a caller
 * relation (duckhts_roh_counts, duckhts_count_error_fit). */
#ifndef DUCKHTS_COUNT_SQL_H
#define DUCKHTS_COUNT_SQL_H

/* A count is cast to `type` only when it is a whole number: the cast alone
 * would round 2.5 to 3. The test is on the value's own type, so a DECIMAL
 * with more digits than a DOUBLE holds is tested exactly; a count must be
 * numeric. */
#define DUCKHTS_WHOLE_COUNT_AS(column, type) \
    "CASE WHEN " column " != trunc(" column ") " \
    "THEN error('read counts must be whole numbers') ELSE CAST(" column " AS " type ") END"
#define DUCKHTS_WHOLE_COUNT(column) DUCKHTS_WHOLE_COUNT_AS(column, "INTEGER")

/* An option is used only when it is a whole number, for the same reason. */
#define DUCKHTS_WHOLE_OPTION(name) \
    "CASE WHEN " name " != trunc(" name ") " \
    "THEN error('" name " must be a whole number') ELSE CAST(" name " AS BIGINT) END"

#endif
