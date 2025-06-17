

# Extracting current package version
VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

.PHONY: document
document:
	Rscript -e "devtools::document()"

.PHONY: clean
clean:
	-rm src/*.so
	-rm src/*.o
	-rm vignettes/*.html
	-rm -r vignettes/*_files/

.PHONY: install build
build: clean document
	@echo Building current version: $(VERSION)
	(cd ../ && R CMD build transitreg)
install: build
	@echo Installing current version: $(VERSION)
	(cd ../ && R CMD INSTALL transitreg_$(VERSION).tar.gz)
check: build
	@echo Checking current version: $(VERSION)
	(cd ../ && R CMD check --as-cran transitreg_$(VERSION).tar.gz)

.PHONY: test
test: clean install
	Rscript -e "library('transitreg'); tinytest::test_all()"

.PHONY: coverage
coverage: install
	Rscript -e "covr::report(covr::package_coverage(line_exclusions = list('src/init.c'), function_exclusions = list('message\\\\s*\\\\(', 'plot.transitreg')), file = \"_coverage.html\")"

.PHONY: readme docs
readme:
	quarto render README.qmd --to gfm

docs: readme
	Rscript -e "altdoc::render_docs(verbose = FALSE, parallel = TRUE)"
	##Rscript -e "altdoc::render_docs(verbose = TRUE)"
