ROOT      := $(abspath ../..)
OUTDIR    := $(ROOT)/bin
KMC_DIR   := ../external/kmc

CXX       ?= c++
override LD := $(CXX)
CXXFLAGS  ?= -O3 -std=gnu++17 -Wall -Wextra -Wno-unused-parameter
CPPFLAGS  += -I $(KMC_DIR)
LDFLAGS   ?=

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SDKROOT    := $(shell xcrun --sdk macosx --show-sdk-path)
  LDFLAGS    += -Wl,-syslibroot,$(SDKROOT)
  LDFLAGS    += -L/opt/homebrew/opt/zlib/lib -L/opt/homebrew/opt/bzip2/lib
  PTHREAD_LIB :=
else
  PTHREAD_LIB := -pthread
endif

LDLIBS    += $(KMC_DIR)/lib/libkmc_core.a -lz -lbz2 $(PTHREAD_LIB)

BIN       := $(OUTDIR)/merge2stats

SRC       := main.cpp
HDRS      := kmc_index.hpp merge_engine.hpp filters.hpp complexity.hpp output.hpp
OBJ       := $(SRC:.cpp=.o)

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	@mkdir -p $(OUTDIR)
	$(LD) $(LDFLAGS) $(OBJ) -o $@ $(LDLIBS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	@$(RM) -f $(OBJ) $(BIN)
