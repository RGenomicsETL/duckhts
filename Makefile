# =============================================================================
# Project configuration
# =============================================================================

.PHONY: help docs clean clean_all clean_local function_catalog \
	test-duckvep-kernel test-duckvep-kernel-asan \
	test-duckvep-kernel-ubsan test-duckvep-kernel-statistical \
	duckvep-generated-check duckvep-upstream-git-check \
	duckvep-state-current-check duckvep-release-conformance-audit \
	test-duckvep-targets-contract duckvep-targets duckvep-cache-receipt \
	test-duckvep-so-conformance \
	duckvep-so-spec duckvep-so-spec-check \
	test-duckvep-witnesses test-duckvep-differential \
	test-duckvep-gvcf-differential \
	test-duckvep-state-exploration \
	test-duckvep-release-vcf \
	duckvep-corpus-differential duckvep-statistical-report \
	duckvep-record-conformance duckvep-record-properties \
	bench-duckvep-throughput bench-duckvep-release-parquet bench-mosdepth \
	duckvep-render-reports \
	test-simd-kernels bench-simd-kernels \
	test-sqllogictest-debug test-sqllogictest-release \
	check-benchmark-portability \
	stage-norm-1000g-dragen-gvcf stage-liftover-references \
	stage-giab-v4.2.1 stage-riker-wgs stage-duckvep-conformance-corpora \
	stage-gffbase stage-duckbedqc-data stage-variantkey-providers \
	test-cache-paths test-benchmark-registry test-variantkey-provider-staging \
	test-duckvep-corpus-staging test-cgranges-benchmark-r \
	test-liftover-property test-liftover-property-asan \
	test-liftover-property-ubsan test-liftover-fuzz test-liftover-fuzz-debug \
	test-bcftools-filter-recovery test-sanitized-extension test-sanitizers

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Extension configuration
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

# -----------------------------------------------------------------------------
# Platform overrides
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Core extension build
# -----------------------------------------------------------------------------

all: configure release

help:
	@printf '%s\n' \
		'Build: make [debug|release|test|test_release]' \
		'Docs: make [docs|function_catalog]' \
		'SIMD: make [test-simd-kernels|bench-simd-kernels]' \
		'DuckVEP: make [test-duckvep-kernel|test-duckvep-differential]' \
		'Reports: make [bench-duckvep-throughput|duckvep-render-reports]' \
		'Data: make stage-[giab-v4.2.1|liftover-references|norm-1000g-dragen-gvcf|riker-wgs|gffbase|duckbedqc-data|variantkey-providers]' \
		'Wasm: make wasm-playwright-test' \
		'Cleanup: make [clean|clean_all|clean_local]'

# DuckDB supplies the generic C API extension targets.
include extension-ci-tools/makefiles/c_api_extensions/base.Makefile
include extension-ci-tools/makefiles/c_api_extensions/c_cpp.Makefile

# Route both the project and inherited extension-ci-tools test targets through
# the cleanup-aware executor. `make test_extension_{debug,release}` therefore
# cannot leave SQL-produced indexes, compressed files, or reports behind.
TEST_RUNNER=$(PYTHON_VENV_BIN) scripts/run_sqllogictest.py

configure: venv platform extension_version

extension_version:
	@$(VERSION_COMMAND)

# Metadata packaging must observe the current root descriptor even when callers
# invoke `make release` or `make debug` without a preceding configure step.
build_extension_with_metadata_debug build_extension_with_metadata_release: extension_version

debug: build_extension_library_debug build_extension_with_metadata_debug
release: build_extension_library_release build_extension_with_metadata_release

# -----------------------------------------------------------------------------
# WebAssembly link ABI
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# General tests and maintenance
# -----------------------------------------------------------------------------

test: test_debug
test_debug: test-cache-paths test-duckvep-kernel test-simd-kernels test-liftover-property test-liftover-fuzz-debug test-sqllogictest-debug
test_release: test-cache-paths test-duckvep-kernel test-simd-kernels test-liftover-property test-liftover-fuzz test-bcftools-filter-recovery test-sqllogictest-release
test_release: test-reference-cache
ifneq ($(filter linux_%,$(or $(DUCKDB_PLATFORM),$(shell sed -n '1p' configure/platform.txt 2>/dev/null))),)
test_release: test-reader-alloc
endif

# Distribution CI tests the same binary inside and outside Docker. Its CMake
# cache contains container-absolute paths, so build this host-only shim afresh.
define run_reader_alloc_test
	@set -e; tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; \
		$(CC) -std=c11 -O1 -g -Wall -Wextra -Werror -fPIC -shared \
			-Iduckdb_capi -Isrc/include test/scripts/reader_alloc_probe.c -ldl \
			-o "$$tmp/reader_alloc_probe.so"; \
		$(1) "$$tmp/reader_alloc_probe.so"
endef

.PHONY: test-reader-alloc test-reader-alloc-r
test-reader-alloc:
	$(call run_reader_alloc_test,./configure/venv/bin/python3 test/scripts/reader_alloc_test.py --extension build/release/duckhts.duckdb_extension --probe)

test-reader-alloc-r:
	$(call run_reader_alloc_test,Rscript test/scripts/reader_alloc_test.R)

.PHONY: test-region-list
test-region-list:
	cmake --build cmake_build/release --target duckhts_region_list_test
	./cmake_build/release/duckhts_region_list_test

.PHONY: test-reference-cache test-reference-cache-asan test-reference-cache-ubsan test-reference-cache-tsan
define run_reference_cache_test
	@if [ "$(DUCKDB_PLATFORM)" = "windows_amd64_mingw" ]; then \
		echo "Skipping POSIX reference probe on MinGW; R package coverage remains enabled"; \
	else \
		set -e; tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; \
		$(CC) -std=c11 -D_DEFAULT_SOURCE -O1 -g -Wall -Wextra -Werror $(1) \
			-Isrc/include -Ithird_party/htslib test/scripts/reference_cache_test.c \
			-Lbuild/release -Wl,-rpath,$(PROJ_DIR)build/release -lduckhts -pthread \
			-o "$$tmp/reference_cache_test"; \
		$(2) "$$tmp/reference_cache_test" "$$tmp" \
			$${REFERENCE_CACHE_SEED:-171} $${REFERENCE_CACHE_TRIALS:-2000}; \
	fi
endef

test-reference-cache:
	$(call run_reference_cache_test,,)

test-reference-cache-asan:
	$(call run_reference_cache_test,-fsanitize=address -fno-omit-frame-pointer,ASAN_OPTIONS=detect_leaks=1:halt_on_error=1)

test-reference-cache-ubsan:
	$(call run_reference_cache_test,-fsanitize=undefined -fno-sanitize-recover=undefined,UBSAN_OPTIONS=halt_on_error=1)

test-reference-cache-tsan:
	$(call run_reference_cache_test,-fsanitize=thread,TSAN_OPTIONS=halt_on_error=1)

define run_liftover_property
	@set -e; tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; \
		$${CC:-cc} -std=c11 -Wall -Wextra -Werror $(1) -Isrc/include \
			test/scripts/liftover_nw_property.c -o "$$tmp/liftover_nw_property"; \
		$(2) "$$tmp/liftover_nw_property" $${LIFTOVER_PROP_SEED:-169} $${LIFTOVER_PROP_TRIALS:-100000}
endef

test-liftover-property:
	$(call run_liftover_property,,)

test-liftover-property-asan:
	$(call run_liftover_property,-fsanitize=address -fno-omit-frame-pointer,ASAN_OPTIONS=detect_leaks=1)

test-liftover-property-ubsan:
	$(call run_liftover_property,-fsanitize=undefined -fno-sanitize-recover=undefined,UBSAN_OPTIONS=halt_on_error=1)

test-liftover-fuzz:
	@if [ "$(DUCKDB_PLATFORM)" = "windows_amd64_mingw" ]; then \
		echo "Skipping Python DuckDB fuzz gate for the MinGW extension platform"; \
	else \
		$(PYTHON_VENV_BIN) test/scripts/fuzz_liftover_sql.py \
			--extension build/release/$(EXTENSION_NAME).duckdb_extension \
			--seed $${LIFTOVER_FUZZ_SEED:-169} --trials $${LIFTOVER_FUZZ_TRIALS:-250}; \
	fi

test-liftover-fuzz-debug:
	@if [ "$(DUCKDB_PLATFORM)" = "windows_amd64_mingw" ]; then \
		echo "Skipping Python DuckDB fuzz gate for the MinGW extension platform"; \
	else \
		$(PYTHON_VENV_BIN) test/scripts/fuzz_liftover_sql.py \
			--extension build/debug/$(EXTENSION_NAME).duckdb_extension \
			--seed $${LIFTOVER_FUZZ_SEED:-169} --trials $${LIFTOVER_FUZZ_TRIALS:-250}; \
	fi

test-sanitized-extension:
	bash scripts/test_sanitized_extension.sh $${SANITIZER:-asan}

test-sanitizers:
	$(MAKE) test-sanitized-extension SANITIZER=asan
	$(MAKE) test-sanitized-extension SANITIZER=ubsan

test-bcftools-filter-recovery:
	@if [ "$(DUCKDB_PLATFORM)" = "windows_amd64_mingw" ]; then \
		echo "Skipping filter recovery probe on MinGW"; \
	else \
		set -e; tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; \
		$${CC:-cc} -std=c11 -Wall -Wextra -Werror -Isrc/include -Ithird_party/htslib \
			test/scripts/bcftools_filter_recovery_probe.c -Lbuild/release \
			-Wl,-rpath,$(PROJ_DIR)build/release -lduckhts -o "$$tmp/filter_recovery"; \
		"$$tmp/filter_recovery"; \
	fi

test-cache-paths:
	bash test/scripts/test_duckhts_cache.sh
	bash test/scripts/test_staging_cache.sh
	bash test/scripts/test_liftover_registry_batch.sh
	bash test/scripts/test_conformance_plugin_cache.sh

test-benchmark-registry: test-variantkey-provider-staging test-duckvep-corpus-staging
	@set -e; tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; \
	(cd "$$tmp" && R CMD build --no-build-vignettes --no-manual "$(PROJ_DIR)r/duckhtsbench"); \
	R CMD INSTALL -l "$$tmp" "$$tmp"/duckhtsbench_*.tar.gz; \
	Rscript -e '.libPaths(c("'"$$tmp"'", .libPaths())); tinytest::test_package("duckhtsbench", testdir = "tinytest")'

test-variantkey-provider-staging:
	bash test/scripts/test_variantkey_provider_staging.sh

test-duckvep-corpus-staging:
	bash test/scripts/test_duckvep_corpus_staging.sh

test-cgranges-benchmark-r:
	bash test/scripts/test_cgranges_benchmark_r.sh

test-sqllogictest-debug: check_configure
	@if [ "$(DUCKDB_PLATFORM)" = "windows_amd64_mingw" ]; then \
		echo "Skipping SQLLogicTest: the Python DuckDB wheel is windows_amd64, not windows_amd64_mingw"; \
	else \
		$(PYTHON_VENV_BIN) scripts/run_sqllogictest.py \
			--test-dir test/sql \
			--external-extension build/debug/$(EXTENSION_NAME).duckdb_extension; \
	fi

test-sqllogictest-release: check_configure
	@if [ "$(DUCKDB_PLATFORM)" = "windows_amd64_mingw" ]; then \
		echo "Skipping SQLLogicTest: the Python DuckDB wheel is windows_amd64, not windows_amd64_mingw"; \
	else \
		$(PYTHON_VENV_BIN) scripts/run_sqllogictest.py \
			--test-dir test/sql \
			--external-extension build/release/$(EXTENSION_NAME).duckdb_extension; \
	fi

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
	rm -f test_bgzip_roundtrip.bed test_formatcols.parquet test_formatcols.vcf.gz.tbi
	rm -f test_range.bam.bai test_targets.bed.gz test_targets.bed.gz.tbi
	rm -f test_auto.fa test_auto.fa.fai test/data/nuc_with_n.fa.fai
	rm -f test/tmp_ce.fa.gz test/tmp_ce.fa.gz.fai test/tmp_ce.fa.gz.gzi

# -----------------------------------------------------------------------------
# Documentation and benchmark renders
# -----------------------------------------------------------------------------

# Render README.md from README.Rmd (GitHub-flavored markdown).
docs: check-benchmark-portability rdm

check-benchmark-portability:
	Rscript r/duckhtsbench/scripts/check_benchmark_portability.R

function_catalog:
	python3 scripts/render_function_catalog.py
rdm: function_catalog
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
bench-lift:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_liftover.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-score:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_score.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-norm:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_norm.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

stage-norm-1000g-dragen-gvcf:
	bash scripts/stage_norm_1000g_dragen_gvcf.sh

stage-liftover-references:
	bash scripts/stage_liftover_references.sh

stage-giab-v4.2.1:
	bash scripts/stage_giab_benchmark_vcfs.sh

stage-riker-wgs:
	bash scripts/stage_riker_wgs_bam.sh

stage-duckvep-conformance-corpora:
	bash scripts/stage_duckvep_conformance_corpora.sh

stage-gffbase:
	bash scripts/stage_gffbase.sh

stage-duckbedqc-data:
	bash scripts/stage_duckbedqc_data.sh

stage-variantkey-providers:
	Rscript r/duckhtsbench/scripts/stage_variantkey_providers.R

bench-mosdepth:
	Rscript -e "rmarkdown::render('benchmarks/Benchmarks_mosdepth.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-variantkey:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_variantkey_conformance.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-variantkey-join:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_variantkey_join_overlap.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-munge:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_munge.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-gffbase:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_gffbase_conformance.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-simd-seq-gc:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_simd_seq_gc.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-simd-bam-gc:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_simd_bam_gc.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

bench-fastq-reader:
	Rscript -e "rmarkdown::render('benchmarks/benchmark_fastq_reader.Rmd', output_format = 'github_document', knit_root_dir = normalizePath('.'))"

# =============================================================================
# SIMD kernel contracts
# =============================================================================

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

# =============================================================================
# DuckVEP kernel and conformance
# =============================================================================

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
	src/duckvep/kernel/src/duckvep_transcript_edit.c \
	src/duckvep/kernel/src/duckvep_hgvs.c \
	src/duckvep/kernel/src/duckvep_codon.c \
	src/duckvep/kernel/src/duckvep_coding.c \
	src/duckvep/kernel/src/duckvep_haplotype.c
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
	perl test/duckvep/conformance/check_state_machine_contract.pl
	perl test/duckvep/upstream/check_sources.pl

duckvep-upstream-git-check:
	perl test/duckvep/upstream/check_sources.pl --verify-upstream-git

duckvep-state-current-check:
	perl test/duckvep/conformance/check_state_machine_contract.pl --require-current

duckvep-release-conformance-audit: duckvep-generated-check \
		duckvep-upstream-git-check duckvep-state-current-check

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
	"$$tmp/duckvep_kernel_property" $(DUCKVEP_PROPERTY_ARGS)

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
		"$$tmp/duckvep_kernel_property" $(DUCKVEP_PROPERTY_ARGS)

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
		"$$tmp/duckvep_kernel_property" $(DUCKVEP_PROPERTY_ARGS)

# Longer deterministic run for rare-state exploration.  Override either value
# to reproduce or extend a run, for example:
#   make test-duckvep-kernel-statistical DUCKVEP_PROP_TRIALS=1000000 DUCKVEP_PROP_SEED=17
# Greatest test-name filtering remains available through
# DUCKVEP_PROPERTY_ARGS, for example `-t hgvs` or `-t haplotype`.  An official
# property-history run leaves that argument empty so the ledger covers the
# complete state machine rather than a favorable subset.
test-duckvep-kernel-statistical:
	DUCKVEP_PROP_TRIALS=$${DUCKVEP_PROP_TRIALS:-100000} \
	DUCKVEP_PROP_SEED=$${DUCKVEP_PROP_SEED:-0xd0c0ffee12345678} \
		$(MAKE) test-duckvep-kernel

test-duckvep-so-conformance:
	Rscript test/duckvep/conformance/so_conformance.R .

test-duckvep-witnesses: release
	Rscript test/duckvep/conformance/generate_witnesses.R \
		--ext build/release/duckhts.duckdb_extension --check

# -----------------------------------------------------------------------------
# Executable VEP differentials
# -----------------------------------------------------------------------------

# VEP 116 is the sole behavioral oracle.  The short target runs all formal
# witnesses against the checked-in transcript fixture.  The corpus target uses
# the same runner for a real VCF and a prepared model database; pass paths and
# sampling size through DUCKVEP_DIFFERENTIAL_ARGS.
test-duckvep-differential: release
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--extension build/release/duckhts.duckdb_extension \
		--sample-per-shape 0 --hgvs

# VEP receives one single-ALT record for each expanded input allele. This pins
# both ALT orders in mixed gVCF records without asking either engine to infer
# genotype remapping or record-level block semantics.
test-duckvep-gvcf-differential: release
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--corpus gvcf_semantics \
		--vcf test/duckvep/conformance/data/gvcf_semantics.vcf \
		--extension build/release/duckhts.duckdb_extension \
		--split-multiallelic --sample-per-shape 0 --distance 0 \
		--max-allele-length 171

# Reproducible rare-state exploration. The generated VCF, pair-level Parquet, and
# statistical summaries remain under the ignored results directory so every failure
# retains its seed and exact allele while large artifacts stay out of git.
test-duckvep-state-exploration: release
	@set -e; \
	cases=$${DUCKVEP_STATE_CASES:-20000}; \
	seed=$${DUCKVEP_STATE_SEED:-17}; \
	max_len=$${DUCKVEP_STATE_MAX_LENGTH:-10}; \
	trials=$${DUCKVEP_PROP_TRIALS:-100000}; \
	Rscript test/duckvep/conformance/property_history.R \
		--trials $$trials --seed $$seed \
		--history test/duckvep/conformance/results/state_exploration_seed_$${seed}_properties.csv \
		--coverage-history test/duckvep/conformance/results/state_exploration_seed_$${seed}_coverage.csv \
		--failure-log-dir test/duckvep/conformance/results; \
	vcf=test/duckvep/conformance/results/state_exploration_seed_$${seed}.vcf; \
	Rscript test/duckvep/conformance/generate_witnesses.R \
		--ext build/release/duckhts.duckdb_extension \
		--random-cases $$cases --max-random-length $$max_len \
		--seed $$seed --out $$vcf; \
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--corpus state_exploration_seed_$${seed} \
		--vcf $$vcf --extension build/release/duckhts.duckdb_extension \
		--sample-per-shape 0 --seed $$seed \
		--max-allele-length $$((max_len + 1)) --hgvs

# Compare a receipt-matched DuckVEP model with the lossless VE relation in an
# official Ensembl variation release VCF. Pass the source/model/receipt paths
# through DUCKVEP_RELEASE_DIFFERENTIAL_ARGS; the runner rebuilds and binds the
# extension to the clean source revision unless explicitly put in diagnostic
# mode.
test-duckvep-release-vcf:
	Rscript test/duckvep/conformance/release_vcf_differential.R \
		--extension build/release/duckhts.duckdb_extension \
		$(DUCKVEP_RELEASE_DIFFERENTIAL_ARGS)

duckvep-corpus-differential: release
	VEP_PREFIX=$(VEP_PREFIX) Rscript test/duckvep/conformance/corpus_differential.R \
		--extension build/release/duckhts.duckdb_extension \
		$(DUCKVEP_DIFFERENTIAL_ARGS)

# -----------------------------------------------------------------------------
# Optional corpus campaigns and evidence receipts
# -----------------------------------------------------------------------------

# Optional coarse-grained campaign orchestration. {targets} owns invalidation and
# resume behavior; corpus_differential.R and blit retain semantic and process ownership.
duckvep-targets:
	Rscript pipelines/duckvep/run.R

# Inventory a VEP cache once after checksum-verified acquisition. Runtime
# campaigns recheck the compact path/size/mtime inventory, not every cache byte.
duckvep-cache-receipt:
	Rscript scripts/duckvep_cache_receipt.R $(DUCKVEP_CACHE_RECEIPT_ARGS)

test-duckvep-targets-contract:
	Rscript test/duckvep/conformance/targets_contract.R

# Reads a Parquet annotation dump produced by the corpus differential. This is
# deliberately not in the ordinary test target because it needs external data.
duckvep-statistical-report:
	Rscript test/duckvep/conformance/statistical_conformance.R \
		$(DUCKVEP_STATISTICAL_ARGS)

# Regenerate the real VEP witness output, then replace this revision's rows in
# the append-only audit ledger. No counts are entered by hand.
duckvep-record-conformance: test-duckvep-differential
	Rscript test/duckvep/conformance/hgvs_history.R \
		--pairs test/duckvep/conformance/results/witnesses_hgvs_pairs.parquet \
		--history test/duckvep/conformance/data/hgvs_history.csv
	Rscript test/duckvep/conformance/statistical_conformance.R \
		--annotations test/duckvep/conformance/results/witnesses_annotations.parquet \
		--history test/duckvep/conformance/data/conformance_history.csv

duckvep-record-properties:
	Rscript test/duckvep/conformance/property_history.R \
		$(DUCKVEP_PROPERTY_HISTORY_ARGS)

# -----------------------------------------------------------------------------
# DuckVEP benchmark and report rendering
# -----------------------------------------------------------------------------

# `configure` refreshes extension metadata from description.yml before timing.
bench-duckvep-throughput: configure release
	Rscript benchmarks/duckvep_throughput.R $(DUCKVEP_THROUGHPUT_ARGS)

# Materialize an official Ensembl consequence VCF in both complete typed and
# narrow release-product projections. The large input and Parquet files remain external.
bench-duckvep-release-parquet: release
	Rscript benchmarks/duckvep_release_parquet.R \
		--extension build/release/duckhts.duckdb_extension \
		$(DUCKVEP_RELEASE_PARQUET_ARGS)

duckvep-render-reports:
	DUCKHTS_REPO_ROOT=$(PROJ_DIR) Rscript -e \
		"rmarkdown::render('benchmarks/duckvep_conformance.Rmd', quiet = TRUE)"
	DUCKHTS_REPO_ROOT=$(PROJ_DIR) Rscript -e \
		"rmarkdown::render('benchmarks/duckvep_throughput.Rmd', quiet = TRUE)"

# =============================================================================
# Browser WebAssembly smoke test
# =============================================================================

# Build the duckdb-wasm extension (Docker) and run the headless Playwright smoke
# test that loads it in a real browser and asserts the SIMD kernels resolve on
# wasm (scalar fallback).  Same build path as CI's wasm-playwright workflow.
wasm-playwright-test:
	@bash -c 'set -euo pipefail; \
		source scripts/duckhts_cache.sh; \
		export NPM_CONFIG_CACHE="$$DUCKHTS_CACHE_DIR/npm"; \
		export PLAYWRIGHT_BROWSERS_PATH="$$DUCKHTS_CACHE_DIR/playwright"; \
		export DUCKHTS_WASM_SITE_ROOT="$$DUCKHTS_CACHE_DIR/wasm/local-artifacts/site"; \
		SERVE=0 bash scripts/start_duckdb_wasm_local_test.sh; \
		cd test/wasm && npm ci && npx playwright install --with-deps chromium && npx playwright test'
