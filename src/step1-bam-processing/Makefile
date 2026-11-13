VERSION := 0.0.2

# Installation directory
prefix ?= /usr/local/bin

# Compiler
CXX := g++

# HTSlib directory
# Can be overridden using:
# make HTSLIB=/path/to/htslib-1.18
HTSLIB ?= dependencies/htslib-1.18

# Convert to absolute path for runtime linking
HTSLIB_ABS := $(abspath $(HTSLIB))

# Compilation flags
CPPFLAGS := -I$(HTSLIB_ABS)
CXXFLAGS := -Wall -Wextra -g -O3 -pthread -std=c++14

# Linker flags
# Embed the absolute HTSlib path into the executable
LDFLAGS := -L$(HTSLIB_ABS) -Wl,-rpath,$(HTSLIB_ABS)

# Libraries
LDLIBS := -lhts -lm -lz -lpthread

# Source files
SOURCES := $(wildcard *.cpp msi/*.cpp methylation/*.cpp)

# Header dependencies
HEADERS := $(wildcard *.h msi/*.h methylation/*.h)

# Target
TARGET := BAMConvertor


.PHONY: all clean install check

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCES) -o $@ $(LDFLAGS) $(LDLIBS)

check:
	@test -f "$(HTSLIB_ABS)/htslib/sam.h" || \
		(echo "ERROR: HTSlib headers not found in $(HTSLIB_ABS)" && exit 1)

	@test -e "$(HTSLIB_ABS)/libhts.so" -o -e "$(HTSLIB_ABS)/libhts.a" || \
		(echo "ERROR: HTSlib library not found in $(HTSLIB_ABS)" && exit 1)

	@echo "HTSlib found: $(HTSLIB_ABS)"

clean:
	rm -f *.o $(TARGET)

install: $(TARGET)
	install -d $(DESTDIR)$(prefix)
	install -m 755 $(TARGET) $(DESTDIR)$(prefix)/
