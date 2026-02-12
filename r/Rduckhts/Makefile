# Makefile for Rduckhts package
# Based on R package conventions and adapted for CRAN compatibility

# Package name
PKG_NAME := Rduckhts
PKG_VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

# Default target
.PHONY: all
all: build

# Install package dependencies (if needed)
.PHONY: deps
deps:
	R -e "if (!require('tinytest', quietly = TRUE)) install.packages('tinytest', repos = 'https://cran.r-project.org')"

# Build the package
.PHONY: build
build: deps
	R CMD build .

.PHONY: install check build rd test test-integration clean cran-ready help

install: build
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

.PHONY: check

check: build
	R CMD check --as-cran  $(PKG_NAME)_$(PKG_VERSION).tar.gz
readme:
	R -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
rd: 
	R -e 'roxygen2::roxygenize(".",load_code=NULL)'

test: install
	R -e "tinytest::test_package('$(PKG_NAME)')"

test2: 
	R -e "library(Rduckhts); tinytest::run_test_dir('inst/tinytest')"
test-integration: install
	R -e "source('inst/tinytest/test_integration.R')"

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PKG_NAME)_*.tar.gz
	rm -rf $(PKG_NAME).Rcheck/
	rm -f src/*.o src/*.so
	cd inst/duckhts_extension/htslib && make clean 2>/dev/null || true
	rm -f inst/duckhts_extension/config.log inst/duckhts_extension/config.status
	rm -f inst/duckhts_extension/htslib/libhts.a inst/duckhts_extension/htslib/libhts.so* inst/duckhts_extension/htslib/htslib.pc
	rm -f inst/extdata/*.duckdb_extension*
	rm -f inst/extdata/*.tmp
	rm -f inst/extdata/*.raw

# Prepare for CRAN submission
.PHONY: cran-ready
cran-ready: check-cran
	@echo "Package ready for CRAN submission"
	@echo "Check the results above for any NOTEs, WARNINGs, or ERRORs"

# Help target
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all        - Build the package (default)"
	@echo "  deps       - Install package dependencies"
	@echo "  build      - Build source package"
	@echo "  install    - Build and install package locally"
	@echo "  check      - Run basic R CMD check"
	@echo "  check-cran - Run CRAN-style check"
	@echo "  test       - Run package tests"
	@echo "  clean      - Clean build artifacts"
	@echo "  cran-ready - Full CRAN check and readiness report"
	@echo "  help       - Show this help message"