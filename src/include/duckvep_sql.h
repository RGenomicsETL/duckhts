#ifndef DUCKVEP_SQL_H
#define DUCKVEP_SQL_H

#include "duckdb_extension.h"

#include <stdbool.h>
#include <stddef.h>

bool duckvep_register_sql_parts(duckdb_connection connection,
	const char *const *parts, size_t part_count);

#endif
