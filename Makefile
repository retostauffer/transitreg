

# 
VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

.PHONY: install
install:
	@echo Current version $(VERSION)
	(cd ../ && \
		R CMD build --no-build-vignettes TransitionModels && \
		R CMD INSTALL TransitionModels_$(VERSION).tar.gz)
