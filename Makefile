

# 
VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

.PHONY: document
document:
	Rscript -e "devtools::document()"

.PHONY: test
test: install
	Rscript -e "library('transitreg'); tinytest::test_all()"

.PHONY: install
install: document
	@echo Installing current version: $(VERSION)
	(cd ../ && \
		R CMD build --no-build-vignettes transitreg && \
		R CMD INSTALL transitreg_$(VERSION).tar.gz)

.PHONY: coverage
coverage: install
	Rscript -e "covr::report(file = \"_coverage.html\")"

.PHONY: check clean
clean:
	-rm src/*.so
	-rm src/*.o
check: clean document
	@echo Checking current version: $(VERSION)
	(cd ../ && \
		R CMD build --no-build-vignettes transitreg && \
		R CMD check --as-cran transitreg_$(VERSION).tar.gz)
