#ifndef DUCKHTS_SOMALIER_PANEL_SQL_H
#define DUCKHTS_SOMALIER_PANEL_SQL_H

/* Returns the CREATE OR REPLACE MACRO text for duckhts_somalier_panel_sha256
   in a malloc'd string the caller frees, or NULL when out of memory. */
char *duckhts_somalier_panel_sha256_macro_sql(void);

#endif
