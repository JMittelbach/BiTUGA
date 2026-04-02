CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
LDFLAGS ?=

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
  ifneq ($(SDKROOT),)
    CXXFLAGS += -isysroot $(SDKROOT)
    LDFLAGS  += -Wl,-syslibroot,$(SDKROOT)
  endif
endif

BIN_DIR = ../../bin

all: $(BIN_DIR)/postprocess_minimatcher $(BIN_DIR)/merge_matches

$(BIN_DIR)/postprocess_minimatcher: postprocess_minimatcher.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

$(BIN_DIR)/merge_matches: merge_matches.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f $(BIN_DIR)/postprocess_minimatcher $(BIN_DIR)/merge_matches
