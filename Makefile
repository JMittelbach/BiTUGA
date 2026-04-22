OUT_BIN_DIR := bin

export LC_ALL ?= C
export LANG ?= C
export LANGUAGE ?= C

CC  ?= cc
CXX ?= c++
UNAME_S := $(shell uname -s)

$(shell mkdir -p $(OUT_BIN_DIR))

BCALM_SRC_DIR     := src/external/bcalm
BCALM_BUILD_DIR   := build_bcalm
BCALM_PATCH_FILE  := patches/bcalm/CMakeLists.txt.gatb

BCALM_CFLAGS := -O3 -DH5_USE_110_API \
                -Wno-implicit-function-declaration \
                -Wno-incompatible-pointer-types

ifeq ($(UNAME_S),Darwin)
  BCALM_CMAKE_OSX := -DCMAKE_OSX_ARCHITECTURES=$(shell uname -m)
  BCALM_SDKROOT   := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
  ifneq ($(BCALM_SDKROOT),)
    BCALM_CMAKE_OSX += -DCMAKE_OSX_SYSROOT=$(BCALM_SDKROOT)
  endif
else
  BCALM_CMAKE_OSX :=
endif

KMC_SRC_DIR  := src/external/kmc
KMC_MAKEFILE := patches/kmc/Makefile
M2S_DIR      := src/merge2stats
MM_DIR       := src/MiniMatcher
POST_MM_DIR  := src/postprocess_MiniMatcher

.PHONY: all build clean distclean help test deps
.PHONY: bcalm clean_bcalm kmc clean_kmc
.PHONY: merge2stats minimatcher postprocess_mm

all: deps build test
	@echo "=> [BiTUGA] Full workflow finished (deps + build + test)."

build: bcalm kmc merge2stats minimatcher postprocess_mm
	@echo "=> [BiTUGA] All components built successfully."

help:
	@echo "Available commands:"
	@echo "  make / make all  - Run deps + build + smoke tests"
	@echo "  make build       - Builds everything (externals & own tools)"
	@echo "  make deps        - Initializes only required submodules"
	@echo "  make test        - Builds everything and runs smoke checks"
	@echo "  make clean       - Cleans all build files"
	@echo "  make bcalm       - Builds only BCALM"
	@echo "  make kmc         - Builds only KMC"
	@echo "  make merge2stats - Builds merge2stats"
	@echo "  make minimatcher - Builds MiniMatcher"
	@echo "  make postprocess_mm - Builds postprocess_MiniMatcher"

test: build
	@echo "=> [BiTUGA] Running smoke tests..."
	@set -e; \
	for exe in \
	  $(OUT_BIN_DIR)/bcalm \
	  $(OUT_BIN_DIR)/kmc \
	  $(OUT_BIN_DIR)/kmc_tools \
	  $(OUT_BIN_DIR)/merge2stats \
	  $(OUT_BIN_DIR)/nt_mini_matcher.x \
	  $(OUT_BIN_DIR)/postprocess_minimatcher \
	  $(OUT_BIN_DIR)/merge_matches; do \
	  if [ ! -x "$$exe" ]; then \
	    echo "ERROR: Missing executable: $$exe"; \
	    exit 1; \
	  fi; \
	done; \
	help_out=$$(mktemp); \
	bash ./BiTUGA.sh --help >$$help_out 2>&1 || true; \
	grep -q "Usage:" $$help_out || { echo "ERROR: BiTUGA help output missing"; cat $$help_out; rm -f $$help_out; exit 1; }; \
	rm -f $$help_out; \
	echo "=> [BiTUGA] Smoke tests passed."

deps:
	@if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
	  git -c submodule.recurse=false submodule update --init $(BCALM_SRC_DIR) $(KMC_SRC_DIR); \
	  if [ -d "$(BCALM_SRC_DIR)" ]; then \
	    git -C "$(BCALM_SRC_DIR)" -c submodule.recurse=false submodule update --init gatb-core; \
	  fi; \
	  if [ -d "$(KMC_SRC_DIR)" ] && git -C "$(KMC_SRC_DIR)" rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
	    git -C "$(KMC_SRC_DIR)" submodule deinit -f -- 3rd_party/cloudflare > /dev/null 2>&1 || true; \
	    rm -rf "$(KMC_SRC_DIR)/3rd_party/cloudflare"; \
	    rm -rf ".git/modules/$(KMC_SRC_DIR)/modules/3rd_party/cloudflare"; \
	  fi; \
	else \
	  echo "=> [BiTUGA] No git worktree detected, skipping dependency initialization."; \
	fi


bcalm: deps $(OUT_BIN_DIR)/bcalm

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
	@cp -f $(BCALM_PATCH_FILE) $(BCALM_SRC_DIR)/gatb-core/gatb-core/CMakeLists.txt
	@cd $(BCALM_BUILD_DIR) && \
	  bcalm_cc="$(CC)"; \
	  bcalm_cxx="$(CXX)"; \
	  tdir=$$(mktemp -d); \
	  printf '%s\n' 'int main(void){return 0;}' > "$$tdir/probe.cpp"; \
	  probe_ok=0; \
	  if command -v "$$bcalm_cxx" >/dev/null 2>&1; then \
	    if [ "$(UNAME_S)" = "Darwin" ]; then \
	      "$$bcalm_cxx" -arch $(shell uname -m) "$$tdir/probe.cpp" -o "$$tdir/probe" >/dev/null 2>&1 && probe_ok=1; \
	    else \
	      "$$bcalm_cxx" "$$tdir/probe.cpp" -o "$$tdir/probe" >/dev/null 2>&1 && probe_ok=1; \
	    fi; \
	  fi; \
	  if [ "$$probe_ok" -ne 1 ]; then \
	    for cand in c++ clang++ g++ g++-15 g++-14 g++-13; do \
	      command -v "$$cand" >/dev/null 2>&1 || continue; \
	      if [ "$(UNAME_S)" = "Darwin" ]; then \
	        "$$cand" -arch $(shell uname -m) "$$tdir/probe.cpp" -o "$$tdir/probe" >/dev/null 2>&1 || continue; \
	      else \
	        "$$cand" "$$tdir/probe.cpp" -o "$$tdir/probe" >/dev/null 2>&1 || continue; \
	      fi; \
	      bcalm_cxx="$$cand"; \
	      case "$$cand" in \
	        clang++) bcalm_cc=clang ;; \
	        g++|g++-*) bcalm_cc=gcc ;; \
	        c++) bcalm_cc=cc ;; \
	      esac; \
	      echo "=> [BCALM] Falling back to '$$bcalm_cxx' for BCALM configuration."; \
	      probe_ok=1; \
	      break; \
	    done; \
	  fi; \
	  if [ "$$probe_ok" -ne 1 ]; then \
	    echo "ERROR: No working C++ compiler found for BCALM."; \
	    rm -rf "$$tdir"; \
	    exit 1; \
	  fi; \
	  rm -rf "$$tdir"; \
	  CC="$$bcalm_cc" CXX="$$bcalm_cxx" CFLAGS="$(BCALM_CFLAGS)" \
	  cmake $(CURDIR)/$(BCALM_SRC_DIR) -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release $(BCALM_CMAKE_OSX) -Wno-dev
	@touch $@

clean_bcalm:
	-$(RM) -r $(BCALM_BUILD_DIR)
	-$(RM) -f $(OUT_BIN_DIR)/bcalm
	@if [ -d "$(BCALM_SRC_DIR)/gatb-core" ] && git -C "$(BCALM_SRC_DIR)/gatb-core" rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
	  git -C "$(BCALM_SRC_DIR)/gatb-core" checkout -- gatb-core/CMakeLists.txt >/dev/null 2>&1 || true; \
	fi

kmc: deps
	@echo "=> [KMC] Building KMC..."
	@if [ -d "$(KMC_SRC_DIR)" ] && git -C "$(KMC_SRC_DIR)" rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
	  git -C "$(KMC_SRC_DIR)" checkout . > /dev/null 2>&1 || true; \
	  git -C "$(KMC_SRC_DIR)" submodule deinit -f -- 3rd_party/cloudflare > /dev/null 2>&1 || true; \
	  rm -rf "$(KMC_SRC_DIR)/3rd_party/cloudflare"; \
	  rm -rf ".git/modules/$(KMC_SRC_DIR)/modules/3rd_party/cloudflare"; \
	fi
	@echo "=> [KMC] Applying portability patches..."
	@grep -q '^#include <regex>$$' $(KMC_SRC_DIR)/kmc_tools/tokenizer.h || \
	  perl -0pi -e 's/#include "defs.h"/#include "defs.h"\n#include <vector>\n#include <regex>/' $(KMC_SRC_DIR)/kmc_tools/tokenizer.h
	@perl -pi -e 's/#include <ext\/algorithm>/#include <algorithm>/' $(KMC_SRC_DIR)/kmc_core/defs.h
	@perl -pi -e 's/using __gnu_cxx::copy_n;/using std::copy_n;/' $(KMC_SRC_DIR)/kmc_core/defs.h
	@perl -pi -e 's/#include <ext\/algorithm>/#include <algorithm>/' $(KMC_SRC_DIR)/kmc_api/kmer_defs.h
	@perl -pi -e 's/using __gnu_cxx::copy_n;/using std::copy_n;/' $(KMC_SRC_DIR)/kmc_api/kmer_defs.h
	@perl -pi -e 's/inline void to_string_impl\(RandomAccessIterator iter\)(?: const)*/inline void to_string_impl(RandomAccessIterator iter) const/' $(KMC_SRC_DIR)/kmc_api/kmer_api.h
	@perl -pi -e 's/inline std::string to_string\(\)(?: const)*/inline std::string to_string() const/' $(KMC_SRC_DIR)/kmc_api/kmer_api.h
	@perl -pi -e 's/inline void to_string\(char \*str\)(?: const)*/inline void to_string(char *str) const/' $(KMC_SRC_DIR)/kmc_api/kmer_api.h
	@perl -pi -e 's/inline void to_long\(std::vector<uint64>& kmer\)(?: const)*/inline void to_long(std::vector<uint64>& kmer) const/' $(KMC_SRC_DIR)/kmc_api/kmer_api.h
	@perl -pi -e 's/inline void to_string\(std::string &str\)(?: const)*/inline void to_string(std::string &str) const/' $(KMC_SRC_DIR)/kmc_api/kmer_api.h
	@perl -pi -e 's@#include \"\.\./3rd_party/cloudflare/zlib\.h\"@#include <zlib.h>@' $(KMC_SRC_DIR)/kmc_core/fastq_reader.h
	@perl -pi -e 's@#include \"\.\./3rd_party/cloudflare/zlib\.h\"@#include <zlib.h>@' $(KMC_SRC_DIR)/kmc_tools/fastq_reader.h
	@mkdir -p $(KMC_SRC_DIR)/include
	@mkdir -p $(KMC_SRC_DIR)/lib
	@cp $(KMC_SRC_DIR)/kmc_api/*.h $(KMC_SRC_DIR)/include/ 2>/dev/null || true
	@$(MAKE) -C "$(KMC_SRC_DIR)" -f "$(CURDIR)/$(KMC_MAKEFILE)" kmc kmc_tools
	@echo "=> [KMC] Build complete."

clean_kmc:
	@if [ -d "$(KMC_SRC_DIR)" ]; then \
	  $(MAKE) -C "$(KMC_SRC_DIR)" -f "$(CURDIR)/$(KMC_MAKEFILE)" clean; \
	  if git -C "$(KMC_SRC_DIR)" rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
	    git -C "$(KMC_SRC_DIR)" checkout .; \
	    git -C "$(KMC_SRC_DIR)" clean -fd; \
	  fi; \
	else \
	  echo "=> [KMC] Source tree not found under $(KMC_SRC_DIR); skipping."; \
	fi
	@$(RM) -f $(OUT_BIN_DIR)/kmc $(OUT_BIN_DIR)/kmc_dump $(OUT_BIN_DIR)/kmc_tools


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
