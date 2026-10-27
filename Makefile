CXX = g++
CPPFLAGS ?=
CXXFLAGS ?= -O2 -Wall -Wextra
LDFLAGS ?=
LDLIBS ?=

SOURCE_DIR := source_for_Cxx_shared_library
TARGET := PAD_Cxx_shared_library.so
SOURCES := $(SOURCE_DIR)/CC_PAD_python_lib.cc \
           $(SOURCE_DIR)/CU_PAD_code.cc \
           $(SOURCE_DIR)/CU_kdtree_with_index.cc

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	@printf 'Compiling %s...\n' '$@'
	@printf '%s\n' '$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++11 -fPIC -pthread -shared $(LDFLAGS) -o $@ $(SOURCE_DIR)/CC_PAD_python_lib.cc $(LDLIBS)'
	@if $(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++11 -fPIC -pthread -shared $(LDFLAGS) -o $@ $(SOURCE_DIR)/CC_PAD_python_lib.cc $(LDLIBS); then \
		printf 'Success: %s is ready to use.\n' '$@'; \
	else \
		status=$$?; \
		printf 'Error: compilation of %s failed. See the compiler messages above.\n' '$@' >&2; \
		exit $$status; \
	fi

clean:
	$(RM) $(TARGET)
