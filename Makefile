# countmut top-level Makefile.
#
# The heavy computation lives in the bundled C core (backend/).  Build it with
# `make backend`.  The `countmut` console script (Python wrapper) auto-builds it
# on first use if it is missing.

.PHONY: backend test help
backend:
	$(MAKE) -C backend

test:
	python -m pytest tests/

help:
	@echo "make backend   # compile the C countmut_core binary"
	@echo "make test      # run the pytest suite"
