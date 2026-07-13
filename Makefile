.PHONY: clean clean_all clean_local function_catalog \
	test-duckvep-kernel test-duckvep-kernel-asan \
	test-duckvep-kernel-ubsan test-duckvep-kernel-statistical \
	duckvep-generated-check test-duckvep-so-conformance \
	duckvep-so-spec duckvep-so-spec-check \
	test-duckvep-witnesses test-duckvep-differential \
	test-duckvep-state-exploration \
	duckvep-corpus-differential duckvep-statistical-report \
	duckvep-record-conformance duckvep-record-properties \
	bench-duckvep-throughput bench-duckvep-release-parquet \
	duckvep-render-reports \
	test-simd-kernels bench-simd-kernels

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Main extension configuration
EXTENSION_NAME=duckhts

# Set to 1 to enable Unstable API (binaries will only work on TARGET_DUCKDB_VERSION, forwards compatibility will be broken)
# WARNING: When set to 1, the duckdb_extension.h from the TARGET_DUCKDB_VERSION must be used, using any other version of
#          the header is unsafe.
USE_UNSTABLE_C_API=0

# The root descriptor is the extension-version authority.  Supplying this value
# prevents extension-ci-tools from falling back to the current Git commit.
EXTENSION_VERSION:=$(shell sed -n 's/^version:[[:space:]]*//p' description.yml)

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

extension_version:
	@$(VERSION_COMMAND)

# Metadata packaging must observe the current root descriptor even when callers
# invoke `make release` or `make debug` without a preceding configure step.
build_extension_with_metadata_debug build_extension_with_metadata_release: extension_version

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
test_debug: test-duckvep-kernel test-simd-kernels test_extension_debug
test_release: test-duckvep-kernel test-simd-kernels test_extension_release

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

# Standalone SIMD contracts link the backend translation units normally and
# obtain private functions only through their registrar callbacks.  The normal
# test is deterministic and quick.  The benchmark repeats that gate before
# timing both nt16 kernels.  Optional: BENCH_ARGS="bases iterations".
SIMD_KERNEL_TEST_SOURCES = \
	test/scripts/bench_simd_kernels.c \
	src/simd/duckhts_simd_scalar.c \
	src/simd/duckhts_simd_avx2.c \
	src/simd/duckhts_simd_neon.c \
	src/simd/duckhts_simd_wasm_simd128.c
SIMD_KERNEL_TEST_CFLAGS = -std=c99 -O2 -Wall -Wextra -Wpedantic -Werror \
	-Wconversion -Wsign-conversion -Wshadow -Wstrict-prototypes

test-simd-kernels:
	@set -e; \
	tmp=$$(mktemp -d "$${TMPDIR:-/tmp}/duckhts-simd-kernels.XXXXXX"); \
	trap 'rm -rf "$$tmp"' EXIT HUP INT TERM; \
	$(CC) $(SIMD_KERNEL_TEST_CFLAGS) -I src/include \
		$(SIMD_KERNEL_TEST_SOURCES) -o "$$tmp/simd_kernels"; \
	"$$tmp/simd_kernels" --check

bench-simd-kernels:
	@set -e; \
	tmp=$$(mktemp -d "$${TMPDIR:-/tmp}/duckhts-simd-kernels.XXXXXX"); \
	trap 'rm -rf "$$tmp"' EXIT HUP INT TERM; \
	$(CC) $(SIMD_KERNEL_TEST_CFLAGS) -I src/include \
		$(SIMD_KERNEL_TEST_SOURCES) -o "$$tmp/simd_kernels"; \
	"$$tmp/simd_kernels" --check; \
	"$$tmp/simd_kernels" --bench $(BENCH_ARGS)

# Pure-C consequence-engine tests ported with the kernel from
# /root/duckvep-c@9f922c8.  The engine is tested through its borrowed-array ABI,
# without DuckDB or htslib.  Each randomized property has an independent
# brute-force, base-walk, genetic-code, rebuild, or composition oracle.
DUCKVEP_KERNEL_SOURCES = \
	src/duckvep/kernel/src/duckvep_kernel.c \
	src/duckvep/kernel/src/duckvep_so.c \
	src/duckvep/kernel/src/duckvep_sweep.c \
	src/duckvep/kernel/src/duckvep_classify.c \
	src/duckvep/kernel/src/duckvep_effect.c \
	src/duckvep/kernel/src/duckvep_sv.c \
	src/duckvep/kernel/src/duckvep_delta.c \
	src/duckvep/kernel/src/duckvep_projection.c \
	src/duckvep/kernel/src/duckvep_codon.c \
	src/duckvep/kernel/src/duckvep_coding.c \
	src/duckvep/kernel/src/duckvep_haplotype.c \
	src/duckvep/duckvep_variant_tile.c
DUCKVEP_PROPERTY_DRIVER = test/duckvep/property/duckvep_kernel_prop.c
DUCKVEP_THEFT_PATCH = test/duckvep/vendor/patches/theft-mingw-no-fork.patch
DUCKVEP_PROPERTY_CPPFLAGS ?=
DUCKVEP_PROPERTY_CFLAGS = -std=c99 -g -O1 -Wall -Wextra \
	-Wno-unused-function -D_DEFAULT_SOURCE -DTHEFT_USE_FLOATING_POINT=0 \
	-I test/duckvep/vendor/greatest \
	-I src/duckvep/kernel/include \
	-I src/duckvep/kernel/src \
	-I src/duckvep

VEP_PREFIX ?= /root/miniconda3/envs/vep
VEP_RUN ?= micromamba run -p $(VEP_PREFIX)
VEP_INSTALL ?= $(firstword $(wildcard $(VEP_PREFIX)/share/ensembl-vep-116*))
VEP_CONSTANTS ?= $(VEP_INSTALL)/Bio/EnsEMBL/Variation/Utils/Constants.pm
VEP_CONSTANTS_SHA256 ?= 98021460cfada22118c6d6b7865bbed3b25c501ca484ba2399380662f1012051

duckvep-so-spec:
	$(VEP_RUN) perl test/duckvep/conformance/extract_so_spec.pl \
		--sha256 $(VEP_CONSTANTS_SHA256) \
		--output test/duckvep/conformance/data/so_consequences.tsv \
		$(VEP_CONSTANTS)

duckvep-so-spec-check:
	$(VEP_RUN) perl test/duckvep/conformance/extract_so_spec.pl \
		--sha256 $(VEP_CONSTANTS_SHA256) \
		--check test/duckvep/conformance/data/so_consequences.tsv \
		$(VEP_CONSTANTS)

duckvep-generated-check:
	perl test/duckvep/conformance/generate_effect_rules.pl --check \
		src/duckvep/kernel/src/duckvep_effect_rules.inc
	perl test/duckvep/conformance/generate_so_metadata.pl --check \
		src/duckvep/kernel/src/duckvep_so_metadata.inc

test-duckvep-kernel: duckvep-generated-check
	@set -e; \
	tmp=$$(mktemp -d "$${TMPDIR:-/tmp}/duckhts-duckvep-kernel.XXXXXX"); \
	trap 'rm -rf "$$tmp"' EXIT HUP INT TERM; \
	cp -R test/duckvep/vendor/theft "$$tmp/theft"; \
	patch --silent --fuzz=0 -d "$$tmp/theft" -p1 < $(DUCKVEP_THEFT_PATCH); \
	$(CC) $(DUCKVEP_PROPERTY_CPPFLAGS) $(DUCKVEP_PROPERTY_CFLAGS) \
		-I "$$tmp/theft/inc" -I "$$tmp/theft/src" \
		$(DUCKVEP_KERNEL_SOURCES) $(DUCKVEP_PROPERTY_DRIVER) \
		"$$tmp"/theft/src/*.c \
		-pthread -o "$$tmp/duckvep_kernel_property"; \
	"$$tmp/duckvep_kernel_property"

test-duckvep-kernel-asan: duckvep-generated-check
	@set -e; \
	tmp=$$(mktemp -d "$${TMPDIR:-/tmp}/duckhts-duckvep-kernel-asan.XXXXXX"); \
	trap 'rm -rf "$$tmp"' EXIT HUP INT TERM; \
	cp -R test/duckvep/vendor/theft "$$tmp/theft"; \
	patch --silent --fuzz=0 -d "$$tmp/theft" -p1 < $(DUCKVEP_THEFT_PATCH); \
	$(CC) $(DUCKVEP_PROPERTY_CPPFLAGS) $(DUCKVEP_PROPERTY_CFLAGS) -fsanitize=address \
		-fno-omit-frame-pointer \
		-I "$$tmp/theft/inc" -I "$$tmp/theft/src" \
		$(DUCKVEP_KERNEL_SOURCES) $(DUCKVEP_PROPERTY_DRIVER) \
		"$$tmp"/theft/src/*.c \
		-pthread -fsanitize=address -o "$$tmp/duckvep_kernel_property"; \
	ASAN_OPTIONS=detect_leaks=1:abort_on_error=1 \
		"$$tmp/duckvep_kernel_property"

test-duckvep-kernel-ubsan: duckvep-generated-check
	@set -e; \
	tmp=$$(mktemp -d "$${TMPDIR:-/tmp}/duckhts-duckvep-kernel-ubsan.XXXXXX"); \
	trap 'rm -rf "$$tmp"' EXIT HUP INT TERM; \
	cp -R test/duckvep/vendor/theft "$$tmp/theft"; \
	patch --silent --fuzz=0 -d "$$tmp/theft" -p1 < $(DUCKVEP_THEFT_PATCH); \
	$(CC) $(DUCKVEP_PROPERTY_CPPFLAGS) $(DUCKVEP_PROPERTY_CFLAGS) -fsanitize=undefined \
		-fno-omit-frame-pointer \
		-I "$$tmp/theft/inc" -I "$$tmp/theft/src" \
		$(DUCKVEP_KERNEL_SOURCES) $(DUCKVEP_PROPERTY_DRIVER) \
		"$$tmp"/theft/src/*.c \
		-pthread -fsanitize=undefined -o "$$tmp/duckvep_kernel_property"; \
	UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
		"$$tmp/duckvep_kernel_property"

# Longer deterministic run for rare-state exploration.  Override either value
# to reproduce or extend a run, for example:
#   make test-duckvep-kernel-statistical DUCKVEP_PROP_TRIALS=1000000 DUCKVEP_PROP_SEED=17
test-duckvep-kernel-statistical:
	DUCKVEP_PROP_TRIALS=$${DUCKVEP_PROP_TRIALS:-100000} \
	DUCKVEP_PROP_SEED=$${DUCKVEP_PROP_SEED:-0xd0c0ffee12345678} \
		$(MAKE) test-duckvep-kernel

test-duckvep-so-conformance:
	Rscript test/duckvep/conformance/so_conformance.R .

test-duckvep-witnesses: release
	Rscript test/duckvep/conformance/generate_witnesses.R \
		--ext build/release/duckhts.duckdb_extension --check

# VEP 116 is the sole behavioral oracle.  The short target runs all formal
# witnesses against the checked-in transcript fixture.  The corpus target uses
# the same runner for a real VCF and a prepared model database; pass paths and
# sampling size through DUCKVEP_DIFFERENTIAL_ARGS.
test-duckvep-differential: release
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--extension build/release/duckhts.duckdb_extension \
		--sample-per-shape 0

# Reproducible rare-state exploration. The generated VCF, pair-level Parquet, and
# statistical summaries remain under the ignored results directory so every failure
# retains its seed and exact allele while large artifacts stay out of git.
test-duckvep-state-exploration: release
	@set -e; \
	cases=$${DUCKVEP_STATE_CASES:-20000}; \
	seed=$${DUCKVEP_STATE_SEED:-17}; \
	max_len=$${DUCKVEP_STATE_MAX_LENGTH:-10}; \
	trials=$${DUCKVEP_PROP_TRIALS:-100000}; \
	DUCKVEP_PROP_TRIALS=$$trials DUCKVEP_PROP_SEED=$$seed \
		$(MAKE) test-duckvep-kernel; \
	vcf=test/duckvep/conformance/results/state_exploration_seed_$${seed}.vcf; \
	Rscript test/duckvep/conformance/generate_witnesses.R \
		--ext build/release/duckhts.duckdb_extension \
		--random-cases $$cases --max-random-length $$max_len \
		--seed $$seed --out $$vcf; \
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--corpus state_exploration_seed_$${seed} \
		--vcf $$vcf --extension build/release/duckhts.duckdb_extension \
		--sample-per-shape 0 --seed $$seed

duckvep-corpus-differential: release
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--extension build/release/duckhts.duckdb_extension \
		$(DUCKVEP_DIFFERENTIAL_ARGS)

# Reads a Parquet annotation dump produced by the corpus differential. This is
# deliberately not in the ordinary test target because it needs external data.
duckvep-statistical-report:
	Rscript test/duckvep/conformance/statistical_conformance.R \
		$(DUCKVEP_STATISTICAL_ARGS)

# Regenerate the real VEP witness output, then replace this revision's rows in
# the append-only audit ledger. No counts are entered by hand.
duckvep-record-conformance: test-duckvep-differential
	Rscript test/duckvep/conformance/statistical_conformance.R \
		--annotations test/duckvep/conformance/results/witnesses_annotations.parquet \
		--history test/duckvep/conformance/data/conformance_history.csv

duckvep-record-properties:
	Rscript test/duckvep/conformance/property_history.R \
		$(DUCKVEP_PROPERTY_HISTORY_ARGS)

# `configure` refreshes extension metadata from description.yml before timing.
bench-duckvep-throughput: configure release
	Rscript benchmarks/duckvep_throughput.R $(DUCKVEP_THROUGHPUT_ARGS)

# Materialize an official Ensembl consequence VCF in both complete typed and
# narrow oracle projections. The large input and Parquet files remain external.
bench-duckvep-release-parquet: release
	Rscript benchmarks/duckvep_release_parquet.R \
		--extension build/release/duckhts.duckdb_extension \
		$(DUCKVEP_RELEASE_PARQUET_ARGS)

duckvep-render-reports:
	DUCKHTS_REPO_ROOT=$(PROJ_DIR) Rscript -e \
		"rmarkdown::render('benchmarks/duckvep_conformance.Rmd', quiet = TRUE)"
	DUCKHTS_REPO_ROOT=$(PROJ_DIR) Rscript -e \
		"rmarkdown::render('benchmarks/duckvep_throughput.Rmd', quiet = TRUE)"

# Build the duckdb-wasm extension (Docker) and run the headless Playwright smoke
# test that loads it in a real browser and asserts the SIMD kernels resolve on
# wasm (scalar fallback).  Same build path as CI's wasm-playwright workflow.
wasm-playwright-test:
	SERVE=0 bash scripts/start_duckdb_wasm_local_test.sh
	cd test/wasm && npm ci && npx playwright install --with-deps chromium && npx playwright test
