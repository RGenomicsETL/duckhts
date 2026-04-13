
# Rduckhts 1.1.6.9000-0.0.6 (Development version)

- Add multi-file reading wrappers: `rduckhts_bam_multi`, `rduckhts_bcf_multi`, `rduckhts_fastq_multi`, `rduckhts_fasta_multi`, `rduckhts_bed_multi`, `rduckhts_tabix_multi`, `rduckhts_gff_multi`, `rduckhts_gtf_multi`. Each follows the standard `(con, table_name, files, ..., overwrite)` convention, creates a DuckDB table with a `filename` column, and accepts an optional `.params` data.frame for per-file parameter overrides (e.g. per-sample regions or index paths). File expansion uses DuckDB's `glob()` so S3 URLs work transparently.
- Add bundled `hts_union_query(reader, pattern, params)` SQL scalar macro for pure-SQL multi-file reading via `SELECT * FROM query(hts_union_query('read_bam', '*.bam'))`.
- Clarify the package README's browser/webR documentation: `README.Rmd` now covers the full `Module.duckhtsWasmHttpConfig` parameter set (`headers`, `allowHosts`, `enforceHostAllowlist`, `withCredentials`, `allowInsecureAuth`), explicitly notes that webR consumers can set that config from R via `webr::eval_js()` without editing the host page, and covers practical wasm/browser behaviors such as same-origin setup, CORS requirements, `.csi` to `.tbi` fallback, and non-fatal `Range` warnings under the local `http.server` harness.
- Use one extension-owned Emscripten compatibility header in the package wasm/webR build: `configure` now includes the shared header from `src/include/` via the bootstrapped `inst/duckhts_extension/include/wasm_socket_compat.h` copy, keeping the bundled browser build aligned with the extension sources without changing native package builds.
- Make the bundled wasm extension self-contained with respect to `htslib`: the Emscripten/webR `configure` path now builds only `libhts.a`, links `duckhts.duckdb_extension` directly against that static archive, and no longer relies on runtime loading of bundled `libhts.so*` files in webR/browser environments.
- Add a browser-native wasm `http` / `https` backend in the bundled extension: `src/wasm_http_hfile.c` now registers a synchronous XHR-backed `htslib` scheme handler from the DuckDB extension entry point, so browser wasm builds can read same-origin and CORS-enabled remote HTS URLs without going through libcurl sockets.
- Keep wasm `libcurl` disabled in `configure`: `r-wasm/webr` ships `/opt/webr/wasm/lib/libcurl.a` and the emcc link test against it passes, but libcurl's `connect()` calls from a SIDE_MODULE still trigger a webR Emscripten message-bus error (`resolved is not a function`) on first network use, so the package-owned XHR backend is the supported wasm HTTP path.
- Harden wasm browser HTTP range behavior in the bundled extension: `wasm_http_hfile.c` now caches object sizes from `Content-Range`/`Content-Length`, clamps range requests when size is known, short-circuits reads at/after EOF, and uses a `GET Range: bytes=0-0` fallback for `SEEK_END` size discovery when `HEAD` metadata is unavailable; this avoids cross-origin 416 failures on `.tbi` index EOF probes (including GTEx tabix in webR/browser).
- Harden non-Range wasm/browser HTTP fallback in the bundled extension: when ranged reads receive `200 OK`, `wasm_http_hfile.c` now caches the full object per open handle and serves later reads from that in-memory cache to avoid repeated full downloads, while still emitting one-time warnings when Range is ignored and when large fallback payloads (>=64 MiB) are used.
- Add optional wasm/browser request-header configuration in the bundled extension via `Module.duckhtsWasmHttpConfig`: supports custom headers (including bearer auth), host allowlisting, optional `withCredentials`, and a default HTTPS-only guard that blocks `Authorization` on non-HTTPS URLs unless `allowInsecureAuth` is explicitly enabled.
- Extend `Module.duckhtsWasmHttpConfig` with `enforceHostAllowlist` in the bundled wasm backend: when enabled, requests to hosts outside `allowHosts` are blocked instead of merely omitting configured headers.
- Fix the bundled wasm side-module final link during `configure`: preserve webR/Emscripten `${LDFLAGS}` on the final `duckhts.duckdb_extension` link so the `SIDE_MODULE` settings reach the extension itself, and export `duckhts_init_c_api` explicitly for DuckDB's loader. This fixes webR/browser `rduckhts_load()` failures where DuckDB could not find a usable init export in `duckhts.duckdb_extension`.
- Set the bundled extension metadata platform to `linux_i686_musl` for the Emscripten/webR path in `configure`, matching the platform value you are using for browser-side loading tests.
- Fix Wasm package builds under `rwasm` / r-universe: the package `configure` script now preserves injected `NAME=VALUE` cache overrides, forwards explicit `--build` / `--host` triplets into the vendored `htslib` `./configure`, forwards webR's Emscripten port flags for `zlib`/`bzip2`, seeds wasm-safe Autoconf cache results for `zlib`/`bzip2`/socket probes, injects a tiny Emscripten-only socket compatibility shim for `recv`/`send`/`closesocket`, and disables the optional `htslib` features that are not available in the stock webR/r-universe wasm toolchain (`libcurl`, `S3`, `GCS`, `lzma`, `plugins`); this fixes the original `ac_cv_func_getrandom=no: command not found` failure and the subsequent nested `htslib` cross-compile probe failures without changing native configure behavior.
- Fix bundled wasm extension artifacts: the package/browser wasm build now includes vendored `htslib` in the linked archive, avoiding unresolved symbols such as `bcf_readrec` at `LOAD`.

# Rduckhts 1.1.6-0.0.2 (2026-04-09)

- Fix `test_bam_file_offset`: cast `COUNT(*)` results to `INTEGER` in SQL so the DuckDB driver returns R `integer` rather than `numeric` (BIGINT maps to double in the duckdb R driver), restoring `expect_identical` assertions.

# Rduckhts 1.1.6-0.0.1 (2026-04-09)

- Fix bundled `read_hts_index_spans(...)` / `rduckhts_hts_index_spans()`: the span view now returns real chunk rows from CSI/TBI/BAI indexes, including populated `bin`, `chunk_beg_vo`, `chunk_end_vo`, `chunk_bytes`, `seq_start`, and `seq_end` values instead of placeholder `NA`s; BCF-backed calls also avoid the previous noisy `tbx` probe warning on `.csi` indexes.
- Add `FILE_OFFSET` column to `rduckhts_bam()` / `read_bam(...)`: exposes the BGZF virtual file offset after each record. Zero runtime overhead (macro over already-open struct fields). Enables `ORDER BY FILE_OFFSET` in SQL `LAG()` / `LAST_VALUE()` window functions to reproduce exact BAM file order for streaming deduplication algorithms. Together with the `//` integer-division operator and `LAST_VALUE(... IGNORE NULLS)`, this permits exact replication of WisecondorX's larp/larp2 state machine in pure SQL, confirmed at 0 mismatches across 25,115 non-zero bins on a real NIPT BAM.

# Rduckhts 1.1.5-0.0.1 (2026-04-08)

- Fix bundled `bcftools_liftover(...)` / `rduckhts_liftover()` cache and realignment hardening: per-thread chain/FASTA contexts are now bounded instead of accumulating for the lifetime of worker threads, and scalar left-alignment no longer reuses stale traceback state after failed/empty alignments.
- Fix bundled `read_bam(...)` / `rduckhts_bam()` and `read_bcf(...)` / `rduckhts_bcf()` indexed parallel full scans when headers contain leading empty contigs: contig claiming now retries iteratively instead of recursively, and the BAM reader no longer returns an empty chunk after successfully handing off to the next contig.
- Fix bundled Windows builds under MinGW and Rtools: vendored `htslib` configuration now distinguishes `windows_amd64_mingw` from `windows_amd64_rtools`, keeping the smaller `configure.win`-style library set on MinGW while restoring the fuller static `libcurl` dependency closure needed on Rtools. `CURL_STATICLIB` remains on built objects rather than `./configure` probes.
- Fix bundled Windows `windows_amd64_rtools` builds: the package build now pins `CC`/`AR`/`RANLIB` from `R CMD config`, avoiding mixed compiler/library selection when vendored `htslib` is configured, and keeps the MinGW static-libcurl configuration aligned with Rtools `libcurl.a`.
- Fix bundled `read_bcf(...)` / `rduckhts_bcf()` mapping of fixed-count INFO/FORMAT arrays: exact-cardinality fields such as `Number=2` and `Number=4` now materialize as DuckDB array/list columns instead of silently dropping all but the first value.
- Fix bundled `read_bcf(...)` / `rduckhts_bcf()` handling of string FORMAT lists such as DRAGEN `FORMAT/LAA`: `Number != 1` string FORMAT fields now materialize as `VARCHAR[]` instead of triggering DuckDB internal assertion failures.
- Fix bundled `duckdb_munge(...)` / `rduckhts_munge()` multithreaded FASTA lookups: FASTA index handles are now thread-local and FASTA fetches are synchronized in `munge`, avoiding intermittent `fai_retrieve` failures and aborts when `fasta_ref` is used with `PRAGMA threads > 1`.
- Add `rduckhts_score()`: polygenic risk score computation backed by the `bcftools +score` plugin, supporting GT/DS/HDS/AP/GP/AS dosage modes, all major GWAS summary presets (PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS, SSF/GWAS-SSF), GWAS-VCF multi-PRS scoring, p-value thresholding, sample subsetting, and region/filter controls.
- Add `rduckhts_munge()`: GWAS summary statistics normalization backed by `bcftools +munge`, with FASTA reference allele resolution, swap-aware effect/frequency transforms, and METAL meta-analysis column support.
- Add `rduckhts_liftover()`: variant coordinate liftover backed by `bcftools +liftover` using UCSC chain files, with full indel normalization, INFO/END lifting, and MT passthrough.
- Add `rduckhts_bed()` for BED3–BED12 interval files and `rduckhts_fasta_nuc()` for nucleotide composition over BED intervals or fixed-width bins.
- Add compression and index helpers: `rduckhts_bgzip()`, `rduckhts_bgunzip()`, `rduckhts_bam_index()`, `rduckhts_bcf_index()`, and `rduckhts_tabix_index()`.
- Add HTS metadata readers: `rduckhts_hts_header()`, `rduckhts_hts_index()`, `rduckhts_hts_index_spans()`, and `rduckhts_hts_index_raw()`.
- Add quality encoding controls to `rduckhts_bam()` and `rduckhts_fastq()` (`quality_representation`, `input_quality_encoding`) and `rduckhts_detect_quality_encoding()` for heuristic FASTQ encoding detection.
- Add `sequence_encoding := 'nt16'` parameter to `rduckhts_bam()`, `rduckhts_fasta()`, and `rduckhts_fastq()` for raw htslib nt16 sequence output as `UTINYINT[]`.
- Add SAM flag helpers `sam_flag_bits()` and `sam_flag_has()`, CIGAR utility functions, and `is_forward_aligned()`.
- Bundle duckhts 1.1.5 extension.

# Rduckhts 0.1.3-0.0.2.9000

# Rduckhts 0.1.3-0.0.2

- Conditionaly enable plugins in windows

- Updates the configure script to avoid check faillure on CRAN MacOS 

- Update the extension version to 0.1.3

# Rduckhts  0.1.2-0.1.5

- Fixed inadvertant removal of libexec
- Updated the plugin to add header table functions

# Rduckhts 0.1.2-0.1.4

- CRAN Submission

# Rduckhts 0.1.2-0.0.9000

- Different fixes for CRAN submission
    - Updated DESCRIPTION Title/Description formatting and added HTSlib reference.
    - Removed default write paths in bootstrap/build helpers; now require explicit paths.
    - setup_hts_env now accepts an explicit plugins_dir parameter.
    - duckhts_build now accepts a make argument (GNU make required).

- modified configure to attemp to support wasm
- Update bootstrapped extension code to match `duckhts` 0.1.2.
- Add SAMtags + auxiliary tag support (standard_tags, auxiliary_tags).
- Add tabix header/typing options (header, header_names, auto_detect, column_types).


# Rduckhts 0.1.1-0.0.3

- make the build single threaded

# Rduckhts 0.1.1-0.0.3

- misspeling correction

# Rduckhts 0.1.1-0.0.2

- CRAN resubmission: apply DuckDB C API header patch to avoid strict-prototypes warnings.

# Rduckhts 0.1.1-0.0.1

- CRAN Submission

- Bump bundled duckhts extension version to 0.1.1.

- Initial development release.
- Bundles the DuckHTS DuckDB extension and htslib for HTS file readers.
- Adds table-creation helpers for VCF/BCF, BAM/CRAM, FASTA/FASTQ, GFF/GTF, and tabix.
