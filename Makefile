# --- Makefile parameters

# Default rule
.DEFAULT_GOAL := fast-test

# Package directories
JULIA_PKG_DIR=src
PYTHON_PKG_DIR=pysrc/xtallography

# Documentation directory
DOCS_DIR=docs

# Julia environment
export JULIA_PROJECT = @.

# Python testing parameters
NPROCS=auto
PYTEST_OPTIONS=-n ${NPROCS}

# --- Testing rules

.PHONY: test fast-test \
	julia-tests julia-coverage \
	python-tests python-coverage

## Run all tests.
test:
	@echo Removing old coverage files
	find . -name "*.jl.*.cov" -exec rm -f {} \;  # Julia coverage
	find . -name "*.coverage.*" -exec rm -f {} \;  # Python coverage
	rm -rf coverage htmlcov coverage.xml  # Python coverage
	@echo Running Julia unit tests
	@$(MAKE) julia-tests
	@echo Running Python unit tests
	@$(MAKE) python-tests
	@echo Generating code coverage reports
	@$(MAKE) julia-coverage
	@$(MAKE) python-coverage


## Run tests in fail-fast mode (i.e., stop at first failure)
fast-test: export JLTEST_FAIL_FAST=true
fast-test: test

## Run only Julia unit tests.
julia-tests:
	jltest --code-coverage test/runtests.jl
	@echo

## Generate basic Julia coverage report
julia-coverage:
	@jlcoverage

## Run only Python unit tests.
python-tests:
	pytest ${PYTEST_OPTIONS}
	@make python-lint
	@echo

.coverage:
	-make python-tests

## Generate basic Python coverage report
python-coverage: .coverage
	coverage report -m

# --- Code quality rules

.PHONY: jlcodestyle python-lint

## Check codestyle for Julia code.
jlcodestyle:
	@echo Checking code style
	@jlcodestyle -v $(JULIA_PKG_DIR)

## Run codestyle and lint checks for Python code.
python-lint:
	@echo
	@SHELL=/bin/bash; \
	MESSAGE=" flake8 start "; \
	BOOKENDS="$$(( (`tput cols` - $${#MESSAGE}) / 2))"; \
	yes "" | head -n $$BOOKENDS | tr \\n "="; \
	printf "$$MESSAGE"; \
	yes "" | head -n $$BOOKENDS | tr \\n "=";
	@echo
	-flake8 ${PYTHON_PKG_DIR}
	@SHELL=/bin/bash; \
	MESSAGE=" flake8 end "; \
	BOOKENDS="$$(( (`tput cols` - $${#MESSAGE}) / 2))"; \
	yes "" | head -n $$BOOKENDS | tr \\n "="; \
	printf "$$MESSAGE"; \
	yes "" | head -n $$BOOKENDS | tr \\n "=";
	@echo
	@echo

# --- Documentation rules

.PHONY: docs

## Generate package documentation.
docs:
	julia --project=${DOCS_DIR} --color=yes -e 'using Pkg; Pkg.update()'
	julia --project=${DOCS_DIR} --color=yes --compile=min -O0 ${DOCS_DIR}/make.jl
	pdoc --math ${PYTHON_PKG_DIR} -o ${DOCS_DIR}/build/python

# --- Utility rules

.PHONY: clean spotless

## Remove files and directories automatically generated during development and testing
## (e.g., coverage files).
clean:
	find . -name "*.jl.*.cov" -exec rm -f {} \;  # Julia coverage
	find . -type d -name "__pycache__" -delete  # compiled python
	find . -type f -name "*.py[co]" -delete  # compiled python
	rm -rf .cache .pytest_cache  # pytest
	find . -name "*.coverage.*" -exec rm -f {} \;  # Python coverage
	rm -rf coverage htmlcov coverage.xml  # Python coverage


## Remove all automatically generated files and directories (e.g., coverage files, package
## documentation, and `Manifest.toml` files).
spotless: clean
	@echo Removing auto-generated package documentatoin
	rm -rf docs/build/
	@echo Removing Manifest.toml files
	find . -name "Manifest.toml" -exec rm -rf {} \;

# --- Makefile Self-Documentation

# Inspired by
# <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
#
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>

.PHONY: help

## Display this list of available rules
help:
	@echo "$$(tput bold)Default rule:$$(tput sgr0) ${.DEFAULT_GOAL}"
	@echo
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
