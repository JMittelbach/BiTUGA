GTTL=.
include ${GTTL}/compile_flags.mk
CFLAGS := $(filter-out -march=%,$(CFLAGS)) -mtune=generic
ifeq ($(MACHINE),x86_64)
  CFLAGS += -mno-avx -mno-avx2 -march=x86-64
endif
NT_MINI_MATCHER_OBJ=option_parser.o nt_mini_matcher.o
OUT_BIN_DIR=../../bin
BIN_TARGET=$(OUT_BIN_DIR)/nt_mini_matcher.x


all: $(BIN_TARGET)

$(BIN_TARGET): ${NT_MINI_MATCHER_OBJ}
	@mkdir -p $(OUT_BIN_DIR)
	@echo "${CXX} ${NT_MINI_MATCHER_OBJ} => $@"
	@$(LD) ${LDFLAGS} ${NT_MINI_MATCHER_OBJ} -o $@ ${LDLIBS}

include ${GTTL}/compile_targets.mk

.PHONY:code_check
code_check:
	@code_check.py -wt `find . -name '*.[hc]pp'`

.PHONY:clean
clean:
	@${RM} *.[odx] $(BIN_TARGET)
	@${RM} -r __pycache__

.PHONY:cxx_options
cxx_options:
	@echo $(CXX) $(TIME_OPTION) -c "<file.cpp>" -o "<file.o>" $(CXXFLAGS) $(CFLAGS) $(CPPFLAGS) -MT "<file.o>" -MMD -MP -MF "<file.d>"
	@echo $(LD) $(LDFLAGS) "<file.o>" -o "<file.x>" $(LDLIBS)
