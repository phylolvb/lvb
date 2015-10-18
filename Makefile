# LVB
# 
# (c) Copyright 2003-2012 by Daniel Barker
# (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
# (c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
# (c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
#  and Chris Wood.
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
# For 32-bit AMD or Intel-based hardware or entirely different 32-bit CPUs
# (e.g. the Raspberry Pi), uncomment the line that adds -DCOMPILE_32_BITS to
# CFLAGS; also, comment-out the line that adds -msse4.2.
#
# For OS X with the Intel C++ compiler, uncomment the definitions of CC,
# G++ and CFLAGS as instructed below, then compile as for Linux.

# Directories
TEST_DIR = ./tests
DOCS_PROG_DIR = ./docs_programmer
LVB_READ_FILE_DIR = ./LVB_READ_FILES/src

### define a c++ compiler to your platform 
G++ = g++

# UNCOMMENT THIS FOR 32-BIT LINUX (E.G. OLD INTEL CPUS OR RASPBERRY PI)
#CFLAGS += -DCOMPILE_32_BITS
# COMMENT THIS FOR 32-BIT (E.G. OLD INTEL CPUS OR RASPBERRY PI):
CFLAGS += -msse4.2

# UNCOMMENT THESE FOR OS X WITH THE INTEL C++ COMPILER:
G++ = icpc
CC = icc
CFLAGS += -openmp-link static

# General options
CFLAGS += -DLVB	 	# Must be present
CFLAGS += -O2 -Wall # -ansi	# Assumes GNU C compiler
#CFLAGS += -fprofile-arcs -ftest-coverage -ansi
#CFLAGS += -g -std=c99
#CFLAGS += -O3 -std=c99 -fopenmp -ftree-loop-distribution -fvariable-expansion-in-unroller -ftree-vectorizer-verbose=6 -msse4.2
CFLAGS += -O3 -std=c99 -fopenmp

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
LVB_LIB = liblvb.$(LIB_EXT)
LIBS += $(LVB_LIB)

LVB_PROG = lvb$(EXE)

# Object files that will go into the LVB library

LVB_LIB_OBJS = admin.$(OBJ) \
               treestack.$(OBJ) \
               cleanup.$(OBJ) \
               datops.$(OBJ) \
               err.$(OBJ) \
               fops.$(OBJ) \
               getparam.$(OBJ) \
               getstartt.$(OBJ) \
               mops.$(OBJ) \
               mymaths.$(OBJ) \
               myuni.$(OBJ) \
               parsim.$(OBJ) \
               randpint.$(OBJ) \
               solve.$(OBJ) \
               sops.$(OBJ) \
               trops.$(OBJ) \
               wrapper.$(OBJ)

LVB_READ_FILE_OBJS = 	$(LVB_READ_FILE_DIR)/CReadFiles.$(OBJ) \
			$(LVB_READ_FILE_DIR)/ReadFile.$(OBJ)

LVB_LIB_OBJS_OUTPUT = $(LVB_LIB_OBJS)

# Object files that are used directly and will not go into the library

LVB_PROG_OBJS = main.$(OBJ)

# Documentation files

LVB_MANUAL = lvb_manual.pdf

DOCS_PROGRAMMER = $(TEST_MANUAL) \
                  $(DOCS_PROG_DIR)/main.html \
                  $(DOCS_PROG_DIR)/treestack.html \
		  $(DOCS_PROG_DIR)/cleanup.html \
		  $(DOCS_PROG_DIR)/datops.html \
		  $(DOCS_PROG_DIR)/err.html \
		  $(DOCS_PROG_DIR)/fops.html \
		  $(DOCS_PROG_DIR)/getparam.html \
		  $(DOCS_PROG_DIR)/admin.html \
		  $(DOCS_PROG_DIR)/mops.html \
		  $(DOCS_PROG_DIR)/mymaths.html \
		  $(DOCS_PROG_DIR)/myuni.html \
		  $(DOCS_PROG_DIR)/parsim.html \
		  $(DOCS_PROG_DIR)/randpint.html \
		  $(DOCS_PROG_DIR)/solve.html \
		  $(DOCS_PROG_DIR)/sops.html \
		  $(DOCS_PROG_DIR)/trops.html \
		  $(DOCS_PROG_DIR)/wrapper.html

TEST_MANUAL = $(DOCS_PROG_DIR)/go.html

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

$(LVB_MANUAL) : lvb_manual.odt
	soffice --headless --convert-to pdf:writer_pdf_Export lvb_manual.odt

LVB_PROG : $(LVB_LIB) $(LVB_PROG_OBJS)
	$(G++) $(CFLAGS) $(LDFLAGS) -o $(LVB_PROG) $(LVB_PROG_OBJS) $(LVB_READ_FILE_OBJS) $(LIBS) $(LM)

$(LVB_LIB) : $(LVB_LIB_OBJS) $(LVB_READ_FILE_OBJS) 
	ar rv $@ $(LVB_LIB_OBJS_OUTPUT) 
	$(RANLIB) $(LVB_LIB)

# If the main test script has changed, we should run the tests
$(TEST_MANUAL) : $(TEST_DIR)/go
	pod2html $< >$@

test : FORCE
	cd tests ; env LVB_EXECUTABLE="`pwd`/../$(LVB_PROG)" LVB_LIBRARY="`pwd`/../$(LVB_LIB)" LVB_OTHERLIBS="$(LM)" LVB_HEADER_PATH=".." GPLUSPLUS="$(G++)" CC="$(CC)" CFLAGS="$(CFLAGS)" ./go; cd ..;

tests : test	# allow 'make tests' as synonym for 'make test'

# make clean to remove files generated by a previous make, leaving
# source files
clean : FORCE
	rm -f $(LVB_PROG) \
	rm -f tests/test*/*.o \
	rm -f tests/test*/testprog.exe \
	rm pod2htmd.tmp \
	rm pod2htmi.tmp \
	$(LVB_LIB) \
	$(LVB_LIB_OBJS) \
	$(LVB_READ_FILE_OBJS) \
	$(LVB_PROG_OBJS) \
	$(DOCS)

FORCE:
