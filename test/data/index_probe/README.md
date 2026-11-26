Copies of `region_union.vcf.gz`, `fixture_mixed.bam` and their indexes, each
beside a one-byte `.csi` that is not an index. The index search of
`SET duckhts_file_system = 'duckdb'` tries `.csi` first, so these files test that
it skips a candidate whose content is not the index format its suffix names.
duckdb-wasm opens a missing name as such a one-byte file
(https://github.com/duckdb/duckdb-wasm/issues/2332).
