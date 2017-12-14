export DEPS_DIR := $(PWD)/deps
export LIB_DIR := $(PWD)/lib
export PDF_DIR := $(abspath deps/ent/data/pdf)
export TURTLE_USE_PNG := 0
unexport CFLAGS

PLUGINS := $(wildcard plugins/*)

.PHONY:  all bin clean lib python plugins $(PLUGINS)

all: bin lib python

bin: bin/danton bin/retro-run

clean:
	@rm -rf bin lib/*.so lib/*.a lib/python/*.so
	@find lib/python -type l -delete

DANTON_LIBS := lib/libalouette.so lib/libdanton.so lib/libent.so lib/libjsmn.a \
	lib/libpumas.so lib/libturtle.so
lib: $(DANTON_LIBS)

bin/danton: deps/danton/src deps/danton/include $(DANTON_LIBS)
	echo "o Building danton executable ..."
	@mkdir -p bin
	@$(MAKE) -C "deps/danton" -e bin/danton
	@mv deps/danton/bin/danton bin
	@echo "--> Done"

bin/retro-run:
	@mkdir -p bin
	@cd bin && ln -s ../scripts/retro-run.py retro-run

define build_library
	echo "o Building $(1) ..."
	@$(MAKE) --directory="$(DEPS_DIR)/$(1)" -e $(2)
	mv $(DEPS_DIR)/$(1)/$(3)/*.$(4) lib
	@echo "--> Done"
endef

lib/lib%.so: deps/%/src deps/%/include
	@$(call build_library,$*,,lib,so)

lib/libdanton.so: deps/danton/src deps/danton/include
	@$(call build_library,danton,lib/libdanton.so,lib,so)

lib/libjsmn.a: deps/jsmn/jsmn.h
	@$(call build_library,jsmn,libjsmn.a,.,a)
	@$(MAKE) --directory="$(DEPS_DIR)/jsmn" clean

python: lib/python/grand_tour.so lib/python/danton.py

lib/python/grand_tour.so: deps/grand-tour/src/grand-tour.c
	echo "o Building grand-tour ..."
	@$(MAKE) --directory="$(DEPS_DIR)/grand-tour" LIB_DIR=$(PWD)/lib/python
	@echo "--> Done"

lib/python/danton.py:
	@cd lib/python && ln -s ../../deps/danton/lib/python/danton.py danton.py

plugins: $(PLUGINS)

$(PLUGINS):
	@echo "o Building the $(notdir $@) plugin ..."
	@make --directory=$@
	@echo "--> Done"
