# Makefile for fmrireg R package
# ===================================

.PHONY: help install check test document coverage coverage-report clean vignettes deps

# Default target
help:
	@echo "fmrireg Makefile"
	@echo "================"
	@echo ""
	@echo "Available targets:"
	@echo "  make install         - Install package and dependencies"
	@echo "  make deps            - Install all dependencies"
	@echo "  make document        - Generate documentation with roxygen2"
	@echo "  make check           - Run R CMD check"
	@echo "  make check-cran      - Run R CMD check with CRAN settings"
	@echo "  make test            - Run all tests"
	@echo "  make test-file FILE=<file> - Run specific test file"
	@echo "  make coverage        - Generate test coverage report (interactive)"
	@echo "  make coverage-report - Generate and open HTML coverage report"
	@echo "  make coverage-ci     - Generate coverage for CI (codecov format)"
	@echo "  make vignettes       - Build vignettes"
	@echo "  make clean           - Clean built files"
	@echo "  make clean-all       - Deep clean (including installed package)"
	@echo "  make spell           - Check spelling"
	@echo "  make url-check       - Check URLs in documentation"
	@echo "  make lint            - Run lintr checks"
	@echo "  make format          - Format code with styler"
	@echo ""

# Installation targets
# ---------------------

install: document
	@echo "Installing fmrireg..."
	R CMD INSTALL --no-multiarch --with-keep.source .

deps:
	@echo "Installing dependencies..."
	Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')"
	Rscript -e "remotes::install_deps(dependencies = TRUE, upgrade = 'never')"

# Documentation targets
# ----------------------

document:
	@echo "Generating documentation..."
	Rscript -e "devtools::document()"

vignettes:
	@echo "Building vignettes..."
	Rscript -e "devtools::build_vignettes()"

# Testing targets
# ----------------

test:
	@echo "Running all tests..."
	Rscript -e "devtools::test()"

test-file:
ifndef FILE
	@echo "Error: FILE not specified. Usage: make test-file FILE=test_ar_integration.R"
	@exit 1
endif
	@echo "Running tests in $(FILE)..."
	Rscript -e "testthat::test_file('tests/testthat/$(FILE)')"

# Coverage targets
# -----------------

coverage: document
	@echo "Generating test coverage report (interactive)..."
	@echo ""
	Rscript -e '\
		if (!requireNamespace("covr", quietly = TRUE)) { \
			stop("covr package not installed. Install with: install.packages(\"covr\")"); \
		}; \
		cov <- covr::package_coverage(type = c("tests", "examples"), quiet = FALSE); \
		print(cov); \
		cat("\n"); \
		covr::report(cov, browse = FALSE); \
		invisible(cov)'

coverage-report: document
	@echo "Generating HTML coverage report..."
	@mkdir -p coverage
	Rscript -e '\
		if (!requireNamespace("covr", quietly = TRUE)) { \
			stop("covr package not installed. Install with: install.packages(\"covr\")"); \
		}; \
		cov <- covr::package_coverage(type = c("tests", "examples"), quiet = FALSE); \
		print(cov); \
		cat("\n=================================\n"); \
		cat("Coverage Summary:\n"); \
		cat("=================================\n"); \
		print(covr::percent_coverage(cov)); \
		cat("\n"); \
		cat("Generating HTML report...\n"); \
		report_file <- covr::report(cov, file = "coverage/coverage.html", browse = FALSE); \
		cat(sprintf("Coverage report saved to: %s\n", normalizePath(report_file))); \
		cat("\nOpening in browser...\n"); \
		browseURL(report_file)'

coverage-ci:
	@echo "Generating coverage for CI (codecov format)..."
	Rscript -e '\
		if (!requireNamespace("covr", quietly = TRUE)) { \
			stop("covr package not installed. Install with: install.packages(\"covr\")"); \
		}; \
		cov <- covr::package_coverage(type = c("tests", "examples"), quiet = FALSE); \
		print(cov); \
		cat("\n"); \
		cat(sprintf("Total coverage: %.2f%%\n", covr::percent_coverage(cov))); \
		covr::codecov(coverage = cov)'

coverage-summary: document
	@echo "Quick coverage summary..."
	Rscript -e '\
		cov <- covr::package_coverage(type = "tests", quiet = TRUE); \
		cat(sprintf("\nTotal test coverage: %.2f%%\n\n", covr::percent_coverage(cov))); \
		df <- covr::tally_coverage(cov, by = "line"); \
		by_file <- aggregate(value ~ filename, data = df, FUN = function(x) sum(x > 0) / length(x) * 100); \
		by_file <- by_file[order(-by_file$$value), ]; \
		names(by_file) <- c("File", "Coverage %"); \
		print(by_file, row.names = FALSE)'

# Check targets
# --------------

check: document
	@echo "Running R CMD check..."
	Rscript -e "devtools::check()"

check-cran: document
	@echo "Running R CMD check with CRAN settings..."
	Rscript -e "devtools::check(cran = TRUE)"

check-fast: document
	@echo "Running fast R CMD check (no vignettes, no examples)..."
	Rscript -e "devtools::check(vignettes = FALSE, run_dont_test = FALSE)"

# Code quality targets
# ---------------------

spell:
	@echo "Checking spelling..."
	Rscript -e "devtools::spell_check()"

url-check:
	@echo "Checking URLs..."
	Rscript -e '\
		if (!requireNamespace("urlchecker", quietly = TRUE)) { \
			install.packages("urlchecker"); \
		}; \
		urlchecker::url_check()'

lint:
	@echo "Running lintr checks..."
	Rscript -e '\
		if (!requireNamespace("lintr", quietly = TRUE)) { \
			stop("lintr not installed. Install with: install.packages(\"lintr\")"); \
		}; \
		lintr::lint_package()'

format:
	@echo "Formatting R code with styler..."
	Rscript -e '\
		if (!requireNamespace("styler", quietly = TRUE)) { \
			stop("styler not installed. Install with: install.packages(\"styler\")"); \
		}; \
		styler::style_pkg()'

# Build targets
# --------------

build: document
	@echo "Building source package..."
	R CMD build .

build-binary: document
	@echo "Building binary package..."
	R CMD INSTALL --build .

# Cleanup targets
# ----------------

clean:
	@echo "Cleaning built files..."
	@rm -rf src/*.o src/*.so src/*.dll
	@rm -rf man/*.Rd
	@rm -rf *.tar.gz
	@rm -rf ..Rcheck
	@rm -rf coverage
	@rm -rf doc
	@rm -rf Meta

clean-all: clean
	@echo "Deep cleaning..."
	@Rscript -e "remove.packages('fmrireg')" 2>/dev/null || true

# CI/CD targets
# --------------

ci-deps:
	@echo "Installing CI dependencies..."
	Rscript -e "remotes::install_deps(dependencies = TRUE, upgrade = 'always')"
	Rscript -e "install.packages(c('covr', 'devtools', 'testthat'))"

ci-check: document
	@echo "Running CI checks..."
	Rscript -e "devtools::check(error_on = 'warning', check_dir = '.', cran = TRUE)"

ci-test:
	@echo "Running CI tests..."
	Rscript -e "devtools::test(reporter = 'summary')"

# Development targets
# --------------------

dev-deps:
	@echo "Installing development dependencies..."
	Rscript -e '\
		pkgs <- c("devtools", "testthat", "covr", "lintr", "styler", \
		          "urlchecker", "roxygen2", "knitr", "rmarkdown"); \
		install.packages(pkgs)'

quick: document test
	@echo "Quick development cycle complete!"

watch:
	@echo "Watching for changes and running tests..."
	@while true; do \
		inotifywait -r -e modify,create,delete R/ tests/ src/ 2>/dev/null && \
		make quick; \
	done

# Debugging targets
# ------------------

debug-test:
ifndef FILE
	@echo "Error: FILE not specified. Usage: make debug-test FILE=test_ar_integration.R"
	@exit 1
endif
	@echo "Running test with browser() support..."
	Rscript -e "testthat::test_file('tests/testthat/$(FILE)', reporter = 'progress')"

show-coverage-files:
	@echo "Files with coverage data:"
	Rscript -e '\
		cov <- covr::package_coverage(type = "tests", quiet = TRUE); \
		df <- covr::tally_coverage(cov, by = "line"); \
		files <- unique(df$$filename); \
		cat(paste(files, collapse = "\n")); \
		cat("\n")'

# Info targets
# -------------

info:
	@echo "Package: fmrireg"
	@echo "Version: $$(grep '^Version:' DESCRIPTION | sed 's/Version: //')"
	@echo "R version: $$(R --version | head -1)"
	@echo "Installed: $$(Rscript -e 'cat(ifelse(requireNamespace(\"fmrireg\", quietly=TRUE), \"Yes\", \"No\"))' 2>/dev/null || echo 'No')"
	@echo ""
	@echo "Test files: $$(ls tests/testthat/test*.R 2>/dev/null | wc -l | tr -d ' ')"
	@echo "R files: $$(ls R/*.R 2>/dev/null | wc -l | tr -d ' ')"
	@echo "C++ files: $$(ls src/*.cpp 2>/dev/null | wc -l | tr -d ' ')"

# Package release targets
# ------------------------

release-check: clean-all deps document check-cran spell url-check
	@echo "Pre-release checks complete!"
	@echo ""
	@echo "Next steps:"
	@echo "  1. Review check results"
	@echo "  2. Update NEWS.md"
	@echo "  3. Update version in DESCRIPTION"
	@echo "  4. Run: make coverage-report"
	@echo "  5. Commit and tag release"

wincheck:
	@echo "Submitting to win-builder..."
	Rscript -e "devtools::check_win_devel()"

rhub-check:
	@echo "Checking on rhub..."
	Rscript -e "rhub::check_for_cran()"
