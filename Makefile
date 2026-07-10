.PHONY: clean clean_all clean_local function_catalog \
	test-duckvep-allele test-duckvep-model

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Main extension configuration
EXTENSION_NAME=duckhts

# Set to 1 to enable Unstable API (binaries will only work on TARGET_DUCKDB_VERSION, forwards compatibility will be broken)
# WARNING: When set to 1, the duckdb_extension.h from the TARGET_DUCKDB_VERSION must be used, using any other version of
#          the header is unsafe.
USE_UNSTABLE_C_API=0

# The DuckDB C API version for extension metadata (stable API = v1.2.0)
TARGET_DUCKDB_VERSION=v1.2.0

# The DuckDB release to fetch headers from
DUCKDB_HEADER_VERSION=v1.5.3

# For MinGW/Rtools we build vendored htslib ourselves.
# Do not inherit the generic DuckDB CI vcpkg + Ninja path here.
ifeq ($(DUCKDB_PLATFORM),windows_amd64_mingw)
override GEN=
override VCPKG_TOOLCHAIN_PATH=
override VCPKG_TARGET_TRIPLET=
override VCPKG_HOST_TRIPLET=
endif
ifeq ($(DUCKDB_PLATFORM),windows_amd64_rtools)
override GEN=
override VCPKG_TOOLCHAIN_PATH=
override VCPKG_TARGET_TRIPLET=
override VCPKG_HOST_TRIPLET=
endif

all: configure release

# Include makefiles from DuckDB
include extension-ci-tools/makefiles/c_api_extensions/base.Makefile
include extension-ci-tools/makefiles/c_api_extensions/c_cpp.Makefile

configure: venv platform extension_version

debug: build_extension_library_debug build_extension_with_metadata_debug
release: build_extension_library_release build_extension_with_metadata_release

# duckdb-wasm is compiled with native WASM EH and WASM_BIGINT.
# Without -fwasm-exceptions at link time, emcc generates dynCall_* wrappers locally
# in the SIDE_MODULE for any call_indirect with i64 types (e.g. htslib's off_t seek
# callback). duckdb-wasm already defines dynCall_* with its own signatures, so a
# locally-defined dynCall_jiji etc. causes "indirect call signature mismatch" at
# extension load time.
#
# Without -sWASM_BIGINT, emcc legalizes i64 imports/exports into split i32 ABIs and
# injects setTempRet0/getTempRet0 shims. duckdb-wasm uses native i64/BigInt imports,
# so the legalized ABI also mismatches at extension load time. Override the CI-tools
# link step to preserve the same EH + BigInt ABI as duckdb-wasm itself.
ifneq ($(DUCKDB_WASM_PLATFORM),)
link_wasm_debug:
	@WASM_LINK_RSP=""; \
	if [ -f "$(EXTENSION_BUILD_PATH)/debug/wasm_link_inputs.rsp" ]; then WASM_LINK_RSP="@$(EXTENSION_BUILD_PATH)/debug/wasm_link_inputs.rsp"; fi; \
	emcc $(EXTENSION_BUILD_PATH)/debug/$(EXTENSION_LIB_FILENAME) $$WASM_LINK_RSP -o $(EXTENSION_BUILD_PATH)/debug/$(EXTENSION_FILENAME_NO_METADATA) -O3 -g -fwasm-exceptions -sWASM_BIGINT -sSIDE_MODULE=2 -sEXPORTED_FUNCTIONS="_$(EXTENSION_NAME)_init_c_api"

link_wasm_release:
	@WASM_LINK_RSP=""; \
	if [ -f "$(EXTENSION_BUILD_PATH)/release/wasm_link_inputs.rsp" ]; then WASM_LINK_RSP="@$(EXTENSION_BUILD_PATH)/release/wasm_link_inputs.rsp"; fi; \
	emcc $(EXTENSION_BUILD_PATH)/release/$(EXTENSION_LIB_FILENAME) $$WASM_LINK_RSP -o $(EXTENSION_BUILD_PATH)/release/$(EXTENSION_FILENAME_NO_METADATA) -O3 -fwasm-exceptions -sWASM_BIGINT -sSIDE_MODULE=2 -sEXPORTED_FUNCTIONS="_$(EXTENSION_NAME)_init_c_api"
endif

test: test_debug
test_debug: test-duckvep-allele test_extension_debug
test_release: test-duckvep-allele test_extension_release

# Override header fetch to use the actual DuckDB release version, not the C API version
update_duckdb_headers_custom:
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb.h', 'duckdb_capi/duckdb.h')"
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb_extension.h', 'duckdb_capi/duckdb_extension.h')"

clean: clean_build clean_cmake
clean_all: clean clean_configure
clean_local:
	rm -rf .duckdb-wasm-local-artifacts .duckdb_wasm_docker_work .webr-local-artifacts duckdb-wasm-local-site .tmp_top_wasm
	rm -rf Rduckhts.Rcheck r/Rduckhts.Rcheck scripts/mosdepth_conformance
	rm -f README.html r/Rduckhts/README.html benchmarks/*.html benchmark_results.csv
	rm -f *.tar.gz *.tgz *.tar.bz2 r/Rduckhts/*.tar.gz r/Rduckhts/*.tgz
	rm -f *.idxstats.txt
	rm -f mosdepth_*_sqltest*
	rm -f test_bgzip_input.bed test_bgzip_input.bed.gz test_bgzip_input.bed.gz.tbi
	rm -f test_bgzip_roundtrip.bed test_formatcols.vcf.gz.tbi test_range.bam.bai
	rm -f test_targets.bed.gz test_targets.bed.gz.tbi
	rm -f test/data/nuc_with_n.fa.fai

# Render README.md from README.Rmd (GitHub-flavored markdown)
function_catalog:
	python3 scripts/render_function_catalog.py
rdm: function_catalog
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
bench:
	Rscript -e "rmarkdown::render('benchmarks/Benchmark.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"
bench-lift:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_liftover.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-score:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_score.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-norm:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_norm.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

stage-norm-1000g-dragen-gvcf:
	scripts/stage_norm_1000g_dragen_gvcf.sh

bench-variantkey:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_variantkey_conformance.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-variantkey-join:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_variantkey_join_overlap.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-munge:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_munge.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-mosdepth:
	Rscript -e "rmarkdown::render('benchmarks/Benchmarks_mosdepth.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-gffbase:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_gffbase_conformance.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-simd-seq-gc:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_simd_seq_gc.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-simd-bam-gc:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_simd_bam_gc.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

# Standalone C microbenchmark for the byte-oriented SIMD kernels: isolates each
# kernel from DuckDB/R, gates correctness against the scalar oracle, and reports
# scalar-relative throughput.  Built at -O2 to match the CRAN R package's
# optimization level (the extension itself ships -O3; the SIMD win holds at
# both).  The per-backend kernels carry their own ISA target attributes, so the
# baseline -O2 here does not gate the AVX2 path.  Optional: BENCH_ARGS="bases iters".
bench-simd-kernels:
	cc -O2 -I src/include -o build/bench_simd_kernels test/scripts/bench_simd_kernels.c
	./build/bench_simd_kernels $(BENCH_ARGS)

# Private scalar contract test for semantic edit trimming and the exact
# VEP-116 differing-region compatibility view.  This deliberately has no
# DuckDB dependency so sanitizers and alternate C compilers can exercise it.
test-duckvep-allele:
	mkdir -p build
	$(CC) -std=c99 -O2 -Wall -Wextra -Wpedantic -Werror \
		-Wconversion -Wsign-conversion -Wshadow -Wstrict-prototypes \
		-I src/duckvep \
		src/duckvep/duckvep_allele.c \
		test/duckvep/test_duckvep_allele.c \
		-o build/test_duckvep_allele
	./build/test_duckvep_allele

test-duckvep-model:
	Rscript test/duckvep_model.R

# Build the duckdb-wasm extension (Docker) and run the headless Playwright smoke
# test that loads it in a real browser and asserts the SIMD kernels resolve on
# wasm (scalar fallback).  Same build path as CI's wasm-playwright workflow.
wasm-playwright-test:
	SERVE=0 bash scripts/start_duckdb_wasm_local_test.sh
	cd test/wasm && npm ci && npx playwright install --with-deps chromium && npx playwright test
