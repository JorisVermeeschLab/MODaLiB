prefix ?= /usr/local/bin
CC ?= g++
LIBS ?=  -Ldependencies/htslib-1.18 -lhts -lm -lz -lstdc++ -lpthread
INCLUDES ?= -Idependencies/htslib-1.18
CPPFLAGS ?= -Wall -g -O3 -pthread --std=c++14
VERSION = 0.0.2

CPPFLAGS += $(INCLUDES)

.PHONY: all clean install version.h

.SUFFIXES:.c .o .cpp

all: BAMConvertor

BAMConvertor: read_info.h serialization.h *.cpp methylation/*.cpp msi/*.cpp
	$(CC) $(CPPFLAGS) $(LIBS) -o BAMConvertor  *.cpp msi/*.cpp methylation/*.cpp -lm -lz -lpthread -lhts

clean:
	rm -f *.o BAMConvertor

install: BAMConvertor
	install BAMConvertor $(prefix)
