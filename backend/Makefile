# Makefile for the countmut C core (htslib-subset backend, from lh3/minipileup).
#
# Builds `countmut_core`, a self-contained binary that does all the counting.
CC ?= cc
LUA_CFLAGS := $(shell pkg-config --cflags lua5.4 2>/dev/null)
ifeq ($(LUA_CFLAGS),)
LUA_CFLAGS := -I/usr/include/lua5.4
endif
LUA_LIBS := $(shell pkg-config --libs lua5.4 2>/dev/null)
ifeq ($(LUA_LIBS),)
LUA_LIBS := -llua5.4
endif
CFLAGS ?= -O3 -Wall -Wno-unused-result -pthread $(LUA_CFLAGS)
OBJS	= razf.o bgzf.o hts.o sam.o faidx.o bedidx.o countmut_expr.o countmut_core.o countmut_core_main.o
LIBS	= -lpthread -lz -lm $(LUA_LIBS)
PROG	= countmut_core

.PHONY: all clean asan

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# AddressSanitizer build for leak/UB auditing
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O1 -g -Wall -fsanitize=address -fno-omit-frame-pointer $(LUA_CFLAGS)" LIBS="-lpthread -lz -lm -fsanitize=address $(LUA_LIBS)" $(PROG)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROG) *.o
