CXX := D:/MingGW/ucrt64/bin/g++.exe
CXXFLAGS := -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra -Wpedantic -MMD -MP -Isrc
LDFLAGS := -fopenmp -static-libgcc -static-libstdc++

BUILD := build
CORE_SRC := src/bookshelf.cpp src/wirelength.cpp src/spectral.cpp src/electric.cpp
CORE_OBJ := $(CORE_SRC:src/%.cpp=$(BUILD)/%.o)

.PHONY: all clean test

all: dreamplace_cpp.exe

$(BUILD):
	if not exist "$(BUILD)" mkdir "$(BUILD)"

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/main.o: src/main.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

dreamplace_cpp.exe: $(CORE_OBJ) $(BUILD)/optimizer.o $(BUILD)/nonsmooth.o $(BUILD)/legalizer.o $(BUILD)/main.o
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BUILD)/core_test.exe: tests/core_test.cpp $(CORE_OBJ) $(BUILD)/legalizer.o | $(BUILD)
	$(CXX) $(CXXFLAGS) tests/core_test.cpp $(CORE_OBJ) $(BUILD)/legalizer.o -o $@ $(LDFLAGS)

test: $(BUILD)/core_test.exe
	./$(BUILD)/core_test.exe

$(BUILD)/legalize_checkpoint.exe: tests/legalize_checkpoint.cpp $(CORE_OBJ) $(BUILD)/legalizer.o | $(BUILD)
	$(CXX) $(CXXFLAGS) tests/legalize_checkpoint.cpp $(CORE_OBJ) $(BUILD)/legalizer.o -o $@ $(LDFLAGS)

clean:
	rm -f $(BUILD)/*.o $(BUILD)/*.exe dreamplace_cpp.exe

-include $(wildcard $(BUILD)/*.d)
