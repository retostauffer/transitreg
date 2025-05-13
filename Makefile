

# 
VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

.PHONY: document
document:
	Rscript -e "devtools::document()"

.PHONY: test
test: clean install
	Rscript -e "library('transitreg'); tinytest::test_all()"

.PHONY: install
install: clean document
	@echo Installing current version: $(VERSION)
	(cd ../ && \
		R CMD build --no-build-vignettes transitreg && \
		R CMD INSTALL transitreg_$(VERSION).tar.gz)

.PHONY: coverage
coverage: install
	Rscript -e "covr::report(covr::package_coverage(line_exclusions = list('src/init.c'), function_exclusions = list('message\\\\s*\\\\(', 'plot.transitreg')), file = \"_coverage.html\")"

.PHONY: check clean
clean:
	-rm src/*.so
	-rm src/*.o
	-rm vignettes/*.html
	-rm -r vignettes/*_files/
check: clean document
	@echo Checking current version: $(VERSION)
	(cd ../ && \
		R CMD build --resave-data --no-build-vignettes transitreg && \
		R CMD check --as-cran transitreg_$(VERSION).tar.gz)

.PHONY: readme docs
readme:
	quarto render README.qmd --to gfm

docs: readme
	Rscript -e "altdoc::render_docs(verbose = FALSE, parallel = TRUE)"
	##Rscript -e "altdoc::render_docs(verbose = TRUE)"
