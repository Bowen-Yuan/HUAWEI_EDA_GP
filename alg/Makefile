ifeq ($(OS),Windows_NT)
MINGW_HOME ?= D:/MingGW/ucrt64
CXX := $(MINGW_HOME)/bin/g++.exe
export PATH := $(MINGW_HOME)/bin;$(PATH)

EXEEXT := .exe
CXXFLAGS := -std=c++17 -O3 -march=native -ffast-math -funroll-loops \
            -fopenmp -Wall -Wextra -DNDEBUG
LDFLAGS := -fopenmp -lm -static-libgcc -static-libstdc++
else
CXX := clang++
OMP_PREFIX := /opt/homebrew/opt/libomp

EXEEXT :=
CXXFLAGS := -std=c++17 -O3 -march=native -ffast-math -funroll-loops \
            -Xclang -fopenmp -I$(OMP_PREFIX)/include \
            -Wall -Wextra -DNDEBUG
LDFLAGS := -L$(OMP_PREFIX)/lib -lomp -lm
endif

SRCDIR := src
OBJDIR := build
TARGET := nsp_placer$(EXEEXT)

SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ifeq ($(OS),Windows_NT)
$(OBJDIR):
	cmd /C if not exist "$(OBJDIR)" mkdir "$(OBJDIR)"

clean:
	cmd /C if exist "$(OBJDIR)" rmdir /S /Q "$(OBJDIR)"
	cmd /C if exist "$(TARGET)" del /Q "$(TARGET)"
else
$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET)
endif

test: $(TARGET)
	./$(TARGET) ispd2005/adaptec1/adaptec1
