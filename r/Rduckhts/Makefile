# Makefile for Rduckhts package
# Based on R package conventions and adapted for CRAN compatibility

# Package name
PKG_NAME := Rduckhts

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

# Install the package locally
.PHONY: install
install: build
	R CMD INSTALL $(PKG_NAME)_*.tar.gz

# Check the package (basic check)
.PHONY: check
check:
	R CMD check --no-manual .
rd: 
	R -e 'roxygen2::roxygenize(".",load_code=NULL)'
# Full CRAN check
.PHONY: check-cran
check-cran:
	R CMD check --as-cran .

# Run tests
.PHONY: test test-integration
test: deps
	R -e "tinytest::test_package('$(PKG_NAME)')"

test-integration: deps
	R -e "source('inst/tinytest/test_integration.R')"

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PKG_NAME)_*.tar.gz
	rm -rf $(PKG_NAME).Rcheck/
	rm -f src/*.o src/*.so

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