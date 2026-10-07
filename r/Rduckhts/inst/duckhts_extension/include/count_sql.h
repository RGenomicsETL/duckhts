/* SQL fragments shared by the macros that read allele counts from a caller
 * relation (duckhts_roh_counts, duckhts_count_error_fit). */
#ifndef DUCKHTS_COUNT_SQL_H
#define DUCKHTS_COUNT_SQL_H

/* A count is cast to INTEGER only when it is a whole number: the cast alone
 * would round 2.5 to 3. */
#define DUCKHTS_WHOLE_COUNT(column) \
    "CASE WHEN CAST(" column " AS DOUBLE) != trunc(CAST(" column " AS DOUBLE)) " \
    "THEN error('read counts must be whole numbers') ELSE CAST(" column " AS INTEGER) END"

#endif
