CXX = c++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra

BIN_DIR = ../../bin

all: $(BIN_DIR)/postprocess_minimatcher $(BIN_DIR)/merge_matches

$(BIN_DIR)/postprocess_minimatcher: postprocess_minimatcher.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN_DIR)/merge_matches: merge_matches.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(BIN_DIR)/postprocess_minimatcher $(BIN_DIR)/merge_matches
