#include "duckdb_extension.h"
#include "named_attribute_columns.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

/* Only projection mapping is exercised here; DuckDB-facing paths are covered by SQL tests. */
duckdb_ext_api_v1 duckdb_ext_api;

int main(void) {
    duckhts_attribute_key keys[] = {
        {"Parent", 6}, {"ID", 2}, {"gene_type", 9}
    };
    duckhts_attribute_columns columns = {keys, 3, 12};
    idx_t ids[] = {14, 3, 12, UINT64_MAX, 13};
    idx_t count = 0;
    duckhts_projected_attribute *projected = duckhts_attribute_columns_project(
        &columns, ids, sizeof(ids) / sizeof(ids[0]), &count);
    assert(projected && count == 3);
    assert(projected[0].key == &keys[2] && projected[0].output_column == 0);
    assert(projected[1].key == &keys[0] && projected[1].output_column == 2);
    assert(projected[2].key == &keys[1] && projected[2].output_column == 4);
    free(projected);

    idx_t other[] = {3, 11, UINT64_MAX};
    projected = duckhts_attribute_columns_project(
        &columns, other, sizeof(other) / sizeof(other[0]), &count);
    assert(projected == NULL && count == 0);
    return 0;
}
