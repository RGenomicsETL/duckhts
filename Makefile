.PHONY: clean clean_all

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

all: configure release

# Include makefiles from DuckDB
include extension-ci-tools/makefiles/c_api_extensions/base.Makefile
include extension-ci-tools/makefiles/c_api_extensions/c_cpp.Makefile

configure: venv platform extension_version

# ---------------------------------------------------------------------------
# Install htslib build dependencies (system dev headers/libs).
# This runs once before cmake; it is a no-op when deps are already present.
# We detect the package manager to stay portable across CI images.
# ---------------------------------------------------------------------------
.PHONY: install_htslib_deps
install_htslib_deps:
	@if [ ! -f /usr/include/zlib.h ] && [ ! -f /usr/local/include/zlib.h ]; then \
		echo "==> Installing htslib build dependencies..."; \
		if command -v apk  >/dev/null 2>&1; then apk add --no-cache zlib-dev bzip2-dev xz-dev curl-dev openssl-dev; \
		elif command -v apt-get >/dev/null 2>&1; then apt-get update -qq && apt-get install -yqq zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libdeflate-dev 2>/dev/null || true; \
		elif command -v yum >/dev/null 2>&1; then yum install -y zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel 2>/dev/null || true; \
		fi; \
	else \
		echo "==> htslib build dependencies already present"; \
	fi

debug: install_htslib_deps build_extension_library_debug build_extension_with_metadata_debug
release: install_htslib_deps build_extension_library_release build_extension_with_metadata_release

test: test_debug
test_debug: test_extension_debug
test_release: test_extension_release

# Override header fetch to use the actual DuckDB release version, not the C API version
update_duckdb_headers:
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb.h', 'duckdb_capi/duckdb.h')"
	$(PYTHON_VENV_BIN) -c "import urllib.request;urllib.request.urlretrieve('https://raw.githubusercontent.com/duckdb/duckdb/$(DUCKDB_HEADER_VERSION)/src/include/duckdb_extension.h', 'duckdb_capi/duckdb_extension.h')"

clean: clean_build clean_cmake
clean_all: clean clean_configure