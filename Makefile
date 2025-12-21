OUT_BIN_DIR := bin

CC  ?= clang
CXX ?= clang++
UNAME_S := $(shell uname -s)

$(shell mkdir -p $(OUT_BIN_DIR))

BCALM_SRC_DIR     := src/external/bcalm
BCALM_BUILD_DIR   := build_bcalm
BCALM_PATCH_FILE  := patches/bcalm/CMakeLists.txt.gatb

BCALM_CFLAGS := -O3 -DH5_USE_110_API \
                -Wno-implicit-function-declaration \
                -Wno-incompatible-pointer-types

ifeq ($(UNAME_S),Darwin)
  BCALM_CMAKE_OSX := -DCMAKE_OSX_ARCHITECTURES=arm64
else
  BCALM_CMAKE_OSX :=
endif

KMC_SRC_DIR  := src/external/kmc
KMC_MAKEFILE := patches/kmc/Makefile
M2S_DIR      := src/merge2stats
MM_DIR       := src/MiniMatcher
POST_MM_DIR  := src/postprocess_MiniMatcher

.PHONY: all clean distclean help
.PHONY: bcalm clean_bcalm kmc clean_kmc
.PHONY: merge2stats minimatcher postprocess_mm

all: bcalm kmc merge2stats minimatcher postprocess_mm
	@echo "=> [BiTUGA] All components built successfully."

help:
	@echo "Available commands:"
	@echo "  make             - Builds everything (externals & own tools)"
	@echo "  make clean       - Cleans all build files"
	@echo "  make bcalm       - Builds only BCALM"
	@echo "  make kmc         - Builds only KMC"
	@echo "  make merge2stats - Builds merge2stats"
	@echo "  make minimatcher - Builds MiniMatcher"
	@echo "  make postprocess_mm - Builds postprocess_MiniMatcher"


bcalm: $(OUT_BIN_DIR)/bcalm

$(OUT_BIN_DIR)/bcalm: $(BCALM_BUILD_DIR)/stamp-build
	@mkdir -p $(OUT_BIN_DIR)
	@cp -f $(BCALM_BUILD_DIR)/bcalm $(OUT_BIN_DIR)/bcalm
	@echo "=> [BCALM] Installed to: $(OUT_BIN_DIR)/bcalm"

$(BCALM_BUILD_DIR)/stamp-build: $(BCALM_BUILD_DIR)/stamp-config
	@echo "=> [BCALM] Compiling..."
	$$(command -v nproc >/dev/null 2>&1 && echo "make -j$$(nproc)" || echo "make -j4") -C $(BCALM_BUILD_DIR)
	@touch $@

$(BCALM_BUILD_DIR)/stamp-config:
	@mkdir -p $(BCALM_BUILD_DIR)
	@git submodule update --init --recursive
	@cp -f $(BCALM_PATCH_FILE) $(BCALM_SRC_DIR)/gatb-core/gatb-core/CMakeLists.txt
	@cd $(BCALM_BUILD_DIR) && \
	  CC=$(CC) CXX=$(CXX) CFLAGS="$(BCALM_CFLAGS)" \
	  cmake $(CURDIR)/$(BCALM_SRC_DIR) -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release $(BCALM_CMAKE_OSX) -Wno-dev
	@touch $@

clean_bcalm:
	-$(RM) -r $(BCALM_BUILD_DIR)
	-$(RM) -f $(OUT_BIN_DIR)/bcalm

kmc:
	@echo "=> [KMC] Building KMC..."
	@git submodule update --init --recursive
	@cd $(KMC_SRC_DIR) && git checkout . > /dev/null 2>&1 || true
	@echo "=> [KMC] Patching headers..."
	@perl -pi -e 's/#include "defs.h"/#include "defs.h"\n#include <vector>\n#include <regex>/' $(KMC_SRC_DIR)/kmc_tools/tokenizer.h
	@perl -pi -e 's/#include <ext\/algorithm>/#include <algorithm>/' $(KMC_SRC_DIR)/kmc_core/defs.h
	@perl -pi -e 's/using __gnu_cxx::copy_n;/using std::copy_n;/' $(KMC_SRC_DIR)/kmc_core/defs.h
	@perl -pi -e 's/#include <ext\/algorithm>/#include <algorithm>/' $(KMC_SRC_DIR)/kmc_api/kmer_defs.h
	@perl -pi -e 's/using __gnu_cxx::copy_n;/using std::copy_n;/' $(KMC_SRC_DIR)/kmc_api/kmer_defs.h
	@cp $(KMC_MAKEFILE) $(KMC_SRC_DIR)/Makefile
	@mkdir -p $(KMC_SRC_DIR)/include
	@mkdir -p $(KMC_SRC_DIR)/lib
	@cp $(KMC_SRC_DIR)/kmc_api/*.h $(KMC_SRC_DIR)/include/ 2>/dev/null || true
	@cd $(KMC_SRC_DIR) && $(MAKE)
	@echo "=> [KMC] Build complete."

clean_kmc:
	-cd $(KMC_SRC_DIR) && git checkout .
	-cd $(KMC_SRC_DIR) && git clean -fd
	-$(RM) -f $(OUT_BIN_DIR)/kmc $(OUT_BIN_DIR)/kmc_dump $(OUT_BIN_DIR)/kmc_tools


merge2stats:
	@echo "=> [BiTUGA] Building merge2stats..."
	$(MAKE) -C $(M2S_DIR)

minimatcher:
	@echo "=> [BiTUGA] Building MiniMatcher..."
	$(MAKE) -C $(MM_DIR)

postprocess_mm:
	@echo "=> [BiTUGA] Building postprocess_MiniMatcher..."
	$(MAKE) -C $(POST_MM_DIR)


clean: clean_bcalm clean_kmc
	@echo "=> [BiTUGA] Cleaning own tools..."
	-$(MAKE) -C $(M2S_DIR) clean
	-$(MAKE) -C $(MM_DIR) clean
	-$(MAKE) -C $(POST_MM_DIR) clean
	-$(RM) -rf $(OUT_BIN_DIR)/*.dSYM
	@echo "=> [BiTUGA] All clean."

distclean: clean