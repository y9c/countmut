# Makefile for the countmut C core (links the full htslib vendored as a git
# submodule, for BAM/CRAM support).
#
# Builds `countmut_core`, a binary that does all the counting.  Requires Lua 5.4.
# The full htslib is built from the vendored submodule (../htslib) and linked.
CC ?= cc
# Vendored Lua 5.4: the source is bundled at backend/lua/ and compiled into
# liblua.a, statically linked into countmut_core so the wheel is fully
# self-contained (no system liblua dependency on any host).
LUA_DIR ?= lua
LUA_CFLAGS := -I$(LUA_DIR)/src
LUA_LIBS := $(LUA_DIR)/src/liblua.a -lm -ldl
# Full htslib vendored as a git submodule (../htslib).  Override with HTSLIB_DIR.
HTSLIB_DIR ?= ../htslib
CFLAGS ?= -O3 -Wall -Wno-unused-result -pthread $(LUA_CFLAGS) -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib
OBJS	= bedidx.o countmut_expr.o countmut_core.o countmut_core_main.o
# Statically link htslib (libhts.a) so the binary is self-contained (no runtime
# libhts.so dependency) -- the PyPI wheel bundles the binary and it must run
# without a separate htslib install.  libhts.a pulls in its own deps.
LIBS	= $(HTSLIB_DIR)/libhts.a -lpthread -lz -lm -lbz2 -llzma -lcurl -lcrypto -ldeflate $(LUA_LIBS)
# rpath $ORIGIN: the wheel bundles liblua5.4.so.0 into countmut/_core/, so the
# binary finds it next to itself (no system liblua5.4 needed).
RPATH	= -Wl,-rpath,$$ORIGIN
PROG	= countmut_core

.PHONY: all clean asan htslib lua

all: htslib lua $(PROG)

htslib:
	$(MAKE) -C $(HTSLIB_DIR)

lua:
	$(MAKE) -C $(LUA_DIR) linux MYCFLAGS='-fPIC'

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(RPATH) -o $@ $^ $(LIBS)

# AddressSanitizer build for leak/UB auditing
asan:
	$(MAKE) clean
	$(MAKE) CFLAGS="-O1 -g -Wall -fsanitize=address -fno-omit-frame-pointer $(LUA_CFLAGS) -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib" LIBS="-L$(HTSLIB_DIR) -lhts -lpthread -lz -lm -fsanitize=address $(LUA_LIBS)" $(PROG)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROG) *.o