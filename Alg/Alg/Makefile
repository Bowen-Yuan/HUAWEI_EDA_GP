CXX := clang++
OMP_PREFIX := /opt/homebrew/opt/libomp
CXXFLAGS := -std=c++17 -O3 -march=native -ffast-math -funroll-loops \
            -Xclang -fopenmp -I$(OMP_PREFIX)/include \
            -Wall -Wextra -DNDEBUG
LDFLAGS := -L$(OMP_PREFIX)/lib -lomp -lm

SRCDIR := src
OBJDIR := build
TARGET := nsp_placer

SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

test: $(TARGET)
	./$(TARGET) ispd2005/adaptec1/adaptec1

clean:
	rm -rf $(OBJDIR) $(TARGET)
