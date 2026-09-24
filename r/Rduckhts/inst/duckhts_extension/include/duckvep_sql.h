#ifndef DUCKVEP_SQL_H
#define DUCKVEP_SQL_H

#include "duckhts_registration.h"

bool duckvep_register_phase_kernels(duckdb_connection connection);
bool duckvep_register_phase_call(duckhts_registration_t *registration);

#endif
