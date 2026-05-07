
# Rduckhts 1.2.1-0.1.0 (2026-05-07)
- compile bundled DuckHTS extension sources with `-Wpedantic` during Unix and Windows package builds while leaving vendored `htslib` on its upstream warning flags
- fix the bundled non-Emscripten `wasm_http_hfile.c` translation unit so native package builds do not warn about an empty source file under pedantic C diagnostics
- harden Windows `configure.win` libcurl detection: the package now requires a successful `curl_easy_init` link using the detected `pkg-config` libcurl dependency closure before enabling htslib remote URL support, and otherwise disables libcurl/S3/GCS cleanly

# Rduckhts 1.2.0-0.1.0 (2026-05-07)
- expose richer bundled GFF/GTF parsed attribute outputs through `rduckhts_gff()` / `rduckhts_gtf()` and multi-file wrappers: `attributes_list = TRUE` returns `MAP(VARCHAR, VARCHAR[])` with grouped multi-values and GFF3 percent-decoding, while `attributes_pairs = TRUE` returns `LIST<STRUCT(key VARCHAR, value VARCHAR, idx INTEGER)>` for exact key/value/index records; `attributes_map = TRUE` remains the backward-compatible raw scalar map
- expose bundled `read_gff(..., strict := true)` through `rduckhts_gff(strict = TRUE)` and `rduckhts_gff_multi(strict = TRUE)`, enabling GFF3 structural validation from R/DBI workflows, including wrong field counts and malformed attribute segments, while keeping the default GFF reader permissive for existing ingestion pipelines
- extend bundled `rduckhts_score()` / `bcftools_score(...)` so `summary_path` can be a character vector or callers can use `summaries_list_file`; multiple TSV/SSF summaries are scored in one genotype scan, `log_path` can write per-PRS matching/audit counts for loaded, matched, allele-mismatch, and duplicate markers, `summaries_list_file` directory scans are deterministic and ignore index sidecars, generated score/count column names are validated for uniqueness, and score accumulation now follows upstream `bcftools +score` float32 summation more closely
- collapse the generated package README function catalog behind a disclosure widget so package users can jump to quick-start and workflow examples more easily
- refresh the package README release docs: clarify the bundled `htslib` 1.23.1/system-requirements wording and redact transient temp-file paths in rendered example output so regenerated README diffs stay deterministic
- add bundled `duckhts_cgranges_overlaps_list(...)`, a vectorized scalar overlap expander that returns LIST-of-STRUCT hit records so DBI queries can expand provider rows with `UNNEST(...)` without generated bulk-probe SQL; package tests cover one-row-per-hit expansion over regular tables and bundled BED data, and the existing `duckhts_cgranges_overlaps_bulk(...)` probe path now also handles DuckDB string vector lengths safely
- fix bundled `duckhts_cgranges_from_query(...)` ingestion of DuckDB string vectors by respecting string lengths instead of assuming NUL-terminated buffers; this fixes cgranges construction from providers such as `read_bed(...)` with longer chromosome names and adds package regression coverage
- add bundled vectorized scalar cgranges probe helpers `duckhts_cgranges_has_overlap(...)` and `duckhts_cgranges_count_overlaps(...)`, enabling DBI queries to stream provider rows through an already-finalized session cgranges index for filtering/count annotations without the materializing `overlaps_bulk` query-string path; add package-level coverage for overlap, contain, and NULL probe semantics
- add bundled `duckhts_cgranges_overlaps_bulk(...)` for SQL-first bulk cgranges probing from R/DBI sessions: one table-function call now streams a query of probe intervals through a finalized cgranges index, supports `mode = 'overlap'|'contain'`, accepts an optional `query_row_id_col`, and otherwise emits 1-based probe ordinals as `query_row_id`; add package-level regression coverage for the new bulk path
- document bundled `duckhts_cgranges_*` entry points in the generated function catalog and package README, add bundled DBI smoke coverage for the session-scoped cgranges registry API, and include a packaged overlap-conformance script reference for `bedtk`-style parity checks
- fix bundled `rduckhts_fasta_nuc()` / `fasta_nuc(...)` GC and AT percentages for intervals containing `N`: `pct_gc` and `pct_at` now use only informative `A/C/G/T` bases in the denominator, so ambiguous bases no longer depress reported bin/interval composition percentages; add bundled regression coverage
- add bundled C-built cgranges bulk-ingest support via `duckhts_cgranges_from_query(...)`, which runs the source query on an extension-owned DuckDB connection and builds the cgranges index in C before publishing it to the session registry; `duckhts_cgranges_from_table(...)` remains deferred for now
- bundle htslib 1.23.1 in the package for the upstream CRAM decoder and GZI validation security fixes, including the wasm/browser-exposed parsing path shipped through `Rduckhts`
- add `rduckhts_bam_bed_coverage()`, bundling native `duckhts_bam_bed_coverage(...)` for samtools coverage-like regional summaries over BED targets with DuckHTS-specific pre/post-filter columns and read-mode strand-specific post summaries; bundled SQL/tinytest coverage now checks expected outputs on the packaged mixed BAM fixture, and `fragment_mode` / `processing_threads` are exposed but currently reserved for later phases
- reduce bundled `rduckhts_bam_bed_coverage()` / `duckhts_bam_bed_coverage(...)` peak memory by allocating and freeing per-region working depth buffers during scan processing instead of retaining them for the whole BED, tile large target intervals internally when computing covered-base breadth, keep the tiled implementation single-pass, align `min_depth > 1` mean-depth behavior with `samtools coverage`, and expose `decompression_threads` so package callers can set htslib BAM/CRAM decode worker counts explicitly
- add `rduckhts_samtools_idxstats()`, bundling native `duckhts_samtools_idxstats(...)` for samtools idxstats-compatible BAM/CRAM/SAM summaries with indexed BAM fast-paths and scan fallback; package SQL/tinytest coverage now checks BAM fast-path output, CRAM fallback output, explicit `index_path`, and overwrite errors
- improve package-source hygiene for local development: ignore generated `README.html`, `.Rcheck`, staged `duckhts_extension/htslib` build outputs, wasm/webR harness byproducts, and stray root-level index files under `r/Rduckhts/`; add top-level `make clean_local` to purge the reproducible package-side artifacts
- add `processing_threads` parameter to `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` for parallel contig processing: workers claim contigs atomically and write output in header order; on the NA12878 WGS benchmark with 2 processing threads, fast mode is 1.38x faster, default mode 1.40x faster, and fragment mode 1.61x faster than mosdepth v0.3.13, all byte-identical; new default is `processing_threads = 2`
- change `rduckhts_mosdepth()` defaults to `threads = 2` (decompression) and `processing_threads = 2` (parallel contigs) for better out-of-the-box WGS performance
- ship htslib public headers and static library in the installed package under `duckhts_extension/htslib/{include,lib}/`; add `inst/htslib_config.R` (generated from `htslib_config.R.in` at configure time) providing `htslib_cflags()`, `htslib_libs()`, `htslib_rpath()`, and `htslib_version()` for downstream R packages that link against the bundled htslib
- fix `configure.win` to stage htslib headers into `include/htslib/` alongside `lib/`, matching Unix configure
- change bundled `bam_bin_counts(...)` / `rduckhts_bam_bin_counts()` to return a dense fixed-bin layout across each selected contig span, including zero-count bins up to the contig end instead of only observed bins; this gives downstream CNV/sample serializers stable per-contig bin shapes, and the package docs/tests now describe and validate the dense contract
- add `rduckhts_bam_bin_counts()` and bundle native `bam_bin_counts(...)` fixed-width BAM/CRAM binning in the package. The new wrapper exposes `mapq`, `require_flags`, `exclude_flags`, and `rmdup = "none"|"flag"|"streaming"` duplicate handling, always returns per-bin forward/reverse totals, and can add per-bin GC/MAPQ summaries via `stats = "gc"`, `"mq"`, or `"gc,mq"`; bundled extdata now includes the tiny WisecondorX BAM/CRAM fixtures used by the new SQL/R tests, and the package README now includes a native bin-count example
- add `rduckhts_mosdepth()` examples to the package README, including windowed fragment coverage output and preview of the generated BED.gz regions file, and refresh the generated function-catalog text so the packaged mosdepth description matches the current v0.3.13 parity surface
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` to cover the pinned local `mosdepth 0.3.13` option surface for indexed BAM/CRAM input: `fragment_mode = TRUE` now matches upstream `--fragment-mode` full-fragment insert coverage for proper pairs, default mode is supported with CIGAR-aware coverage plus mate-overlap correction, `read_groups = "..."` filters RG tags, `min_frag_len` / `max_frag_len` filter absolute template length, and `use_median = TRUE` switches `by = "<window|bed>"` outputs from mean to median; add bundled SQL/R/conformance coverage for BAM and CRAM fast/fragment/default/median cases.
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` fast-mode parity with `quantize = "..."`, writing mosdepth-style `.quantized.bed.gz` + CSI output, and add bundled tests for quantized output plus explicit `by = "<bed>"` validation.
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` fast-mode parity with `thresholds = "..."` for `by = "<window|bed>"`, writing mosdepth-style `.thresholds.bed.gz` + CSI outputs; also align window/BED mean accumulation and window-region distribution bucketing with upstream mosdepth's current implementation behavior, and add bundled SQL/R/native-conformance coverage for the new outputs.
- Bundle upstream mosdepth edge-case fixtures (`big`, `empty-tids`, `overlapping-pairs`, `ovl`, `nanopore`, and related BED files) in `inst/extdata/` for stronger mosdepth parity testing, and record Brent Pedersen as the original mosdepth author in the package metadata/copyright bundle.
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)`: the native mosdepth-compatible fast-mode rewrite now accepts indexed CRAM input via `fasta = ...` when required by htslib, and exposes `precision_digits = 2` as an explicit wrapper argument instead of relying on the `MOSDEPTH_PRECISION` environment variable; add bundled BAM/CRAM tests plus explicit precision validation.
- Expand `README.Rmd` with runnable compression/indexing examples covering `rduckhts_bgzip()`, `rduckhts_bgunzip()`, `rduckhts_bam_index()`, `rduckhts_bcf_index()`, and `rduckhts_tabix_index()`, then regenerate the rendered package README outputs.
- Add `decompression_threads` to `rduckhts_bam()` and `rduckhts_bam_multi()`, matching the bundled `read_bam(..., decompression_threads := 2)` SQL parameter. The previous hardcoded htslib worker-thread count is now the documented default, and `0` disables per-file worker threads.
- Speed up bundled zero-column `COUNT(*)` queries across the HTS readers: `read_bam(...)`, `read_bcf(...)`, `read_tabix(...)`, `read_gff(...)`, `read_gtf(...)`, and indexed `read_bed(...)` now use index metadata for full-file count-only scans when DuckDB projects no output columns; `read_fasta(...)` uses `faidx` sequence counts when an index is available and otherwise counts FASTA headers directly; `read_fastq(...)` continues to count raw FASTQ records directly when no projected columns are needed, while preserving paired/interleaved validation errors.
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
