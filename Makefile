CXX := D:/MingGW/ucrt64/bin/g++.exe
CXXFLAGS := -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra -Wpedantic -MMD -MP -Iinclude
LDFLAGS := -fopenmp -static-libgcc -static-libstdc++

BUILD := build
CORE_SRC := src/batch_acceptance.cpp src/bisection.cpp src/bookshelf.cpp src/coarse_flow.cpp src/compact_recovery.cpp src/density.cpp src/density_coordinate.cpp src/hpwl.cpp src/lambda_controller.cpp src/optimizer.cpp src/placer.cpp src/recovery.cpp src/swap_recovery.cpp src/transport.cpp
CORE_OBJ := $(CORE_SRC:src/%.cpp=$(BUILD)/%.o)

.PHONY: all clean test

all: $(BUILD)/epsilon_active.exe

$(BUILD):
	if not exist "$(BUILD)" mkdir "$(BUILD)"

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/main.o: src/main.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/epsilon_active.exe: $(CORE_OBJ) $(BUILD)/main.o
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BUILD)/core_tests.exe: tests/core_tests.cpp $(CORE_OBJ) | $(BUILD)
	$(CXX) $(CXXFLAGS) tests/core_tests.cpp $(CORE_OBJ) -o $@ $(LDFLAGS)

test: $(BUILD)/core_tests.exe
	./$(BUILD)/core_tests.exe

clean:
	del /q $(BUILD)\*.o $(BUILD)\*.d $(BUILD)\*.exe 2>NUL

-include $(wildcard $(BUILD)/*.d)
