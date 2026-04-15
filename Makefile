.PHONY: clean clean_all function_catalog

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
DUCKDB_HEADER_VERSION=v1.4.3

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
test_debug: test_extension_debug
test_release: test_extension_release

# Override header fetch to use the actual DuckDB release version, not the C API version
update_duckdb_headers_custom:
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb.h', 'duckdb_capi/duckdb.h')"
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb_extension.h', 'duckdb_capi/duckdb_extension.h')"

clean: clean_build clean_cmake
clean_all: clean clean_configure

# Render README.md from README.Rmd (GitHub-flavored markdown)
function_catalog:
	python3 scripts/render_function_catalog.py
rdm: function_catalog
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
bench:
	Rscript -e "rmarkdown::render('Benchmark.Rmd', output_format = 'github_document')"
bench-lift:
	Rscript -e "rmarkdown::render('benchmark_liftover.Rmd', output_format = 'github_document')"

bench-score:
	Rscript -e "rmarkdown::render('benchmark_score.Rmd', output_format = 'github_document')"

bench-munge:
	Rscript -e "rmarkdown::render('benchmark_munge.Rmd', output_format = 'github_document')"

bench-mosdepth:
	Rscript -e "rmarkdown::render('Benchmarks_mosdepth.Rmd', output_format = 'github_document')"
