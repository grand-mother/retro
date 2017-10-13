export DEPS_DIR := $(PWD)/deps
export LIB_DIR := $(PWD)/lib
export PDF_DIR := $(abspath deps/ent/data/pdf)
export TURTLE_USE_PNG := 0

.PHONY:  all bin clean lib python

all: bin lib python

bin: bin/danton

clean:
	@rm -rf bin lib

lib: lib/libalouette.so lib/libdanton.so lib/libent.so lib/libjsmn.a           \
	lib/libpumas.so lib/libturtle.so

bin/danton: deps/danton/src deps/danton/include lib
	@mkdir -p bin
	@$(MAKE) -C "deps/danton" -e bin/danton
	@mv deps/danton/bin/danton bin

define build_library
	echo "o Building $(1) ..."
	@$(MAKE) --directory="$(DEPS_DIR)/$(1)" -e $(2)
	@mkdir -p lib && mv $(DEPS_DIR)/$(1)/$(3)/*.$(4) lib
	@echo "--> Done"
endef

lib/lib%.so: deps/%/src deps/%/include
	@$(call build_library,$*,,lib,so)

lib/libdanton.so:
	@$(call build_library,danton,lib/libdanton.so,lib,so)

lib/libjsmn.a: deps/jsmn/jsmn.h
	@$(call build_library,jsmn,libjsmn.a,.,a)
	@$(MAKE) --directory="$(DEPS_DIR)/jsmn" clean

python: lib/python/grand_tour.py lib/python/danton.py

lib/python/grand_tour.py: deps/grand-tour/lib/python/grand_tour.py
	@mkdir -p lib/python && cp $< $@

lib/python/danton.py: deps/danton/lib/python/danton.py
	@mkdir -p lib/python && cp $< $@
