prefix ?= /usr/local/bin
CC ?= g++
LIBS ?=  -Ldependencies/htslib-1.13 -lhts -lm -lz -lstdc++ -lpthread
INCLUDES ?= -Idependencies/htslib-1.13
CPPFLAGS ?= -Wall -g -O3 -pthread
OBJS = perRead.o
VERSION = 0.6.1

CPPFLAGS += $(INCLUDES)

.PHONY: all clean install version.h

.SUFFIXES:.c .o .cpp

all: BAMConvertor

BAMConvertor:
	$(CC) $(CPPFLAGS) $(LIBS) -o BAMConvertor  main.cpp -lm -lz -lpthread -lhts

clean:
	rm -f *.o BAMConvertor

install: BAMConvertor
	install BAMConvertor $(prefix)
