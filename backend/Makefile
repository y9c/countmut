# Makefile for the countmut C core (links the full htslib vendored as a git
# submodule, for BAM/CRAM support).
#
# Builds `countmut_core`, a binary that does all the counting.  Requires Lua 5.4.
# The full htslib is built from the vendored submodule (../htslib) and linked.
CC ?= cc
LUA_CFLAGS := $(shell pkg-config --cflags lua5.4 2>/dev/null)
ifeq ($(LUA_CFLAGS),)
LUA_CFLAGS := $(shell pkg-config --cflags lua5.3 2>/dev/null)
endif
ifeq ($(LUA_CFLAGS),)
LUA_CFLAGS := -I/usr/include/lua5.4
endif
LUA_LIBS := $(shell pkg-config --libs lua5.4 2>/dev/null)
ifeq ($(LUA_LIBS),)
LUA_LIBS := $(shell pkg-config --libs lua5.3 2>/dev/null)
endif
ifeq ($(LUA_LIBS),)
LUA_LIBS := -llua5.4
endif
# Full htslib vendored as a git submodule (../htslib).  Override with HTSLIB_DIR.
HTSLIB_DIR ?= ../htslib
CFLAGS ?= -O3 -Wall -Wno-unused-result -pthread $(LUA_CFLAGS) -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib
OBJS	= bedidx.o countmut_expr.o countmut_core.o countmut_core_main.o
# Statically link htslib (libhts.a) so the binary is self-contained (no runtime
# libhts.so dependency) -- the PyPI wheel bundles the binary and it must run
# without a separate htslib install.  libhts.a pulls in its own deps.
LIBS	= $(HTSLIB_DIR)/libhts.a -lpthread -lz -lm -lbz2 -llzma -lcurl -lcrypto -ldeflate $(LUA_LIBS)
PROG	= countmut_core

.PHONY: all clean asan htslib

all: htslib $(PROG)

htslib:
	$(MAKE) -C $(HTSLIB_DIR)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# AddressSanitizer build for leak/UB auditing
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O1 -g -Wall -fsanitize=address -fno-omit-frame-pointer $(LUA_CFLAGS) -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib" LIBS="-L$(HTSLIB_DIR) -lhts -lpthread -lz -lm -fsanitize=address $(LUA_LIBS)" $(PROG)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROG) *.o