.DEFAULT_GOAL := all

CXX ?= g++
CPPFLAGS ?=
BASE_CXXFLAGS := -O3 -std=c++11 -Wall -Wextra -Wpedantic
CXXFLAGS ?=
LDFLAGS ?=
LDLIBS ?=

TARGET := src/libpad_python.so
SOURCE := src/cpp/CC_PAD_python_lib.cc
DEPENDENCIES := src/cpp/CU_PAD_code.cc src/cpp/CU_kdtree_with_index.cc

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(SOURCE) $(DEPENDENCIES) Makefile
	$(CXX) $(CPPFLAGS) $(BASE_CXXFLAGS) $(CXXFLAGS) -fPIC -fvisibility=hidden -pthread -shared $(SOURCE) $(LDFLAGS) -Wl,--no-undefined -o $@ $(LDLIBS)

clean:
	$(RM) $(TARGET)
