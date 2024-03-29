# LVB
#
# (c) Copyright 2003-2012 by Daniel Barker.
# (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
# (c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
# (c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
# and Chris Wood.
# (c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
# Fernando Guntoro, Maximilian Strobl and Chris Wood.
# (c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
# Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
# and Martyn Winn.
# (c) Copyright 2023 by Joseph Guscott and Daniel Barker.
# (c) Copyright 2022 by Joseph Guscott and Daniel Barker.
# (c) Copyright 2023 by Joseph Guscott and Daniel Barker.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# Makefile - makefile for LVB

# Requires GNU make
# Requires pod2html (this will usually be present if Perl is installed)
# Requires LibreOffice (and soffice must be on PATH)
#
# It is always safer to make from scratch after changing the source code
# or Makefile, e.g.
#
# make clean
#
# first.
#
# For 64-bit Linux on AMD or Intel-based hardware, compile with
#
# make
#
# AFTER BUILDING LVB, RUN TESTS IMMEDIATELY: DO
#
# make test
#
#
# Directories
TEST_DIR = ./test/src
DOCS_PROG_DIR = ./docs_programmer
LVB_SRC_DIR = ./src
#LVB_READ_FILE_DIR = ./LVB_READ_FILES/src

### define compilers and compiler options
G++ = g++
CC = g++
CXX = g++

# No warnings

CFLAGS += -O3 -fopenmp -std=c++11
CXXFLAGS += -O3 -fopenmp -std=c++11

# Warnings

# CFLAGS += -Wall -g -O3 -fopenmp -std=c++11
# CXXFLAGS += -Wall -g -O3 -fopenmp -std=c++11

# Profiling

# CFLAGS += -g -O3 -fopenmp -std=c++11 -ftest-coverage -fprofile-arcs -pg
# CXXFLAGS += -g -O3 -fopenmp -std=c++11 -ftest-coverage -fprofile-arcs -pg

# System-dependent macros - OK for Linux and UNIX-like systems, for others will
# require change
LM = -lm		# UNIX
RANLIB = ranlib	# UNIX
EXE =			# UNIX
OBJ = o			# UNIX
LIB_EXT = a		# UNIX

%.OBJ : %.c
	$(CC) $(CFLAGS) $<

# LVB library
LVB_LIB = libLVB.$(LIB_EXT)
LIBS += $(LVB_LIB)
LVB_PROG = lvb$(EXE)

# Object files that will go into the LVB library

# C files

LVB_LIB_OBJS = $(LVB_SRC_DIR)/Admin.$(OBJ) \
               $(LVB_SRC_DIR)/Cleanup.$(OBJ) \
			   $(LVB_SRC_DIR)/Clock.$(OBJ) \
               $(LVB_SRC_DIR)/DataOperations.$(OBJ) \
               $(LVB_SRC_DIR)/Error.$(OBJ) \
               $(LVB_SRC_DIR)/FileOperations.$(OBJ) \
			   $(LVB_SRC_DIR)/LogFile.$(OBJ) \
               $(LVB_SRC_DIR)/MemoryOperations.$(OBJ) \
               $(LVB_SRC_DIR)/MyMaths.$(OBJ) \
               $(LVB_SRC_DIR)/RandomNumberGenerator.$(OBJ) \
               $(LVB_SRC_DIR)/TreeEvaluation.$(OBJ) \
			   $(LVB_SRC_DIR)/SearchParameters.$(OBJ) \
               $(LVB_SRC_DIR)/Solve.$(OBJ) \
               $(LVB_SRC_DIR)/Sops.$(OBJ) \
			   $(LVB_SRC_DIR)/StartingTemperature.$(OBJ) \
			   $(LVB_SRC_DIR)/Treestack.$(OBJ) \
               $(LVB_SRC_DIR)/TreeOperations.$(OBJ) \
               $(LVB_SRC_DIR)/Wrapper.$(OBJ)

# C++ files

LVB_READ_FILE_OBJS = 	$(LVB_SRC_DIR)/CommandLineParser.$(OBJ) \
						$(LVB_SRC_DIR)/Hash.$(OBJ) \
						$(LVB_SRC_DIR)/MSAInput.$(OBJ) \
						$(LVB_SRC_DIR)/Print.$(OBJ) \
						$(LVB_SRC_DIR)/Verbose.$(OBJ)
						

LVB_LIB_OBJS_OUTPUT = $(LVB_LIB_OBJS)

# Object files that are used directly and will not go into the library

LVB_PROG_OBJS = $(LVB_SRC_DIR)/Main.$(OBJ)

# Documentation files

# TEST_MANUAL = $(DOCS_PROG_DIR)/go.html

# All documentation files

DOCS = $(DOCS_PROGRAMMER) \
       $(TEST_MANUAL) \
       $(LVB_MANUAL)

# Default rule for documentation files (not suitable for user's
# documentation files, which are derived from *.pod not *.c, or for
# go.html, which is derived from a different directory)

%.html : $(LVB_LIB)
	pod2html $(notdir $*).c > $@

all : lvb	# allow 'make all' as synonym for 'make lvb'

lvb : LVB_PROG $(DOCS)

#$(LVB_MANUAL) : LVBManual.odt
#	soffice --headless --convert-to pdf:writer_pdf_Export LVBManual.odt

LVB_PROG : $(LVB_LIB) $(LVB_PROG_OBJS)
	$(G++) $(CFLAGS) $(LDFLAGS) -o $(LVB_PROG) $(LVB_PROG_OBJS) $(LVB_READ_FILE_OBJS) $(LIBS) $(LM)

$(LVB_LIB) : $(LVB_LIB_OBJS) $(LVB_READ_FILE_OBJS)
	ar rv $@ $(LVB_LIB_OBJS_OUTPUT)
	$(RANLIB) $(LVB_LIB)

# If the main test script has changed, we should run the tests
$(TEST_MANUAL) : $(TEST_DIR)/COMMON/go
	pod2html $< >$@

test : FORCE
	cd test/src/COMMON ; env LVB_EXECUTABLE="`pwd`/../../../$(LVB_PROG)" LVB_LIBRARY="`pwd`/../../../$(LVB_LIB)" LVB_OTHERLIBS="$(LM)" LVB_HEADER_PATH=".." GPLUSPLUS="$(G++)" LINKERCPLUSPLUS="$(G++)" CFLAGS="$(CFLAGS)" ./go; cd ..;
	
tests : test	# allow 'make tests' as synonym for 'make test'

# make clean to remove files generated by a previous make, leaving
# source files
clean : FORCE
	rm -f $(LVB_PROG) \
	rm -f test/test*/*.o \
	rm -f test/test*/testprog.exe \
	rm pod2htmd.tmp \
	rm pod2htmi.tmp \
	$(LVB_LIB) \
	$(LVB_LIB_OBJS) \
	$(LVB_READ_FILE_OBJS) \
	$(LVB_PROG_OBJS) \
	$(DOCS)

FORCE:
