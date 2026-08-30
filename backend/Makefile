# Makefile for the countmut C core (htslib-subset backend, from lh3/minipileup).
#
# Builds `countmut_core`, a self-contained binary that does all the counting.
CC ?= cc
CFLAGS ?= -O3 -Wall -Wno-unused-result -pthread
OBJS	= razf.o bgzf.o hts.o sam.o faidx.o bedidx.o countmut_core.o countmut_core_main.o
LIBS	= -lpthread -lz -lm
PROG	= countmut_core

.PHONY: all clean asan

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# AddressSanitizer build for leak/UB auditing
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O1 -g -Wall -fsanitize=address -fno-omit-frame-pointer" LIBS="-lpthread -lz -lm -fsanitize=address" $(PROG)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROG) *.o
