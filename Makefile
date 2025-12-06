# Top-level Makefile to build BCALM from the BiTUGA root

OUT_BIN_DIR := bin
BCALM_DIR   := src/external/bcalm
BUILD_DIR   := build/bcalm

UNAME_S := $(shell uname -s)

CC  ?= clang
CXX ?= clang++

CFLAGS_EXTRA := -O3 -DH5_USE_110_API \
                -Wno-implicit-function-declaration \
                -Wno-incompatible-pointer-types

ifeq ($(UNAME_S),Darwin)
  CMAKE_OSX := -DCMAKE_OSX_ARCHITECTURES=arm64
else
  CMAKE_OSX :=
endif

.PHONY: all bcalm clean-bcalm distclean-bcalm

# Default target: build BCALM
all: bcalm

bcalm: $(OUT_BIN_DIR)/bcalm

$(OUT_BIN_DIR)/bcalm: $(BUILD_DIR)/stamp-build
	@mkdir -p $(OUT_BIN_DIR)
	@cp -f $(BUILD_DIR)/bcalm $(OUT_BIN_DIR)/bcalm
	@echo "=> installed: $(OUT_BIN_DIR)/bcalm"

$(BUILD_DIR)/stamp-build: $(BUILD_DIR)/stamp-config
	$$(command -v nproc >/dev/null 2>&1 && echo "make -j$$(nproc)" || echo "make -j4") -C $(BUILD_DIR)
	@touch $@

$(BUILD_DIR)/stamp-config:
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && \
	  CC=$(CC) CXX=$(CXX) CFLAGS="$(CFLAGS_EXTRA)" \
	  cmake $(abspath $(BCALM_DIR)) \
	        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
	        -DCMAKE_BUILD_TYPE=Release $(CMAKE_OSX) -Wno-dev
	@touch $@

clean-bcalm:
	-$(RM) -r $(BUILD_DIR)

distclean-bcalm: clean-bcalm
	-$(RM) -f $(OUT_BIN_DIR)/bcalm
