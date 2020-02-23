/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2019 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/* ----------------------------------------------------------------------
   MR-MPI = MapReduce-MPI library
   http://www.cs.sandia.gov/~sjplimp/mapreduce.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the modified Berkeley Software Distribution (BSD) License.

   See the README file in the top-level MapReduce directory.
------------------------------------------------------------------------- */

#ifndef KEYVALUE_H_
#define KEYVALUE_H_

//#ifdef MAP_REDUCE_SINGLE

  #include <mpi.h>
  #include "stdio.h"
  #include "stdint.h"

  namespace MAPREDUCE_NS {

class KeyValue {
  friend class MapReduce;

 public:
    uint64_t nkv;                   // # of KV pairs in entire KV on this proc
    uint64_t ksize;                 // exact size of all key data
    uint64_t vsize;                 // exact size of all value data
    uint64_t esize;                 // exact size of all data in KV
    uint64_t fsize;                 // size of KV file
    int msize;                      // size of largest KV pair across all procs

    char *page;                     // in-memory page
    int memtag;                     // memory page ID
    int npage;                      // # of pages in entire KV

    KeyValue(class MapReduce *, int, int,
      class Memory *, class Error *, MPI_Comm);
    ~KeyValue();

    void allocate();
    void set_page(uint64_t, char *, int);
    void deallocate(int);
    void truncate(int, int, uint64_t);
    void copy(KeyValue *);
    void append();
    void complete();
    void complete_dummy();
    int request_info(char **);
    int request_page(int, uint64_t &, uint64_t &, uint64_t &);
    void overwrite_page(int);
    void close_file();

    void add(char *, int, char *, int);
    void add(int, char *, int, char *, int);
    void add(int, char *, int *, char *, int *);

    void print(FILE *, int, int, int);

 private:
    MapReduce *mr;
    MPI_Comm comm;
    class Memory *memory;
    class Error *error;
    int me;

    uint64_t pagesize;                // size of page
    int kalign, valign;                // alignment for keys & values
    int talign;                       // alignment of entire KV pair
    int kalignm1, valignm1, talignm1;   // alignments-1 for masking
    int twolenbytes;                  // size of single key,value lengths

    // in-memory page

    int nkey;                         // # of KV pairs in page
    uint64_t keysize;                 // exact size of key data in page
    uint64_t valuesize;               // exact size of value data in page
    uint64_t alignsize;               // current size of page with alignment

    // virtual pages

    struct Page {
      uint64_t keysize;               // exact size of keys
      uint64_t valuesize;             // exact size of values
      uint64_t exactsize;             // exact size of all data in page
      uint64_t alignsize;             // aligned size of all data in page
      uint64_t filesize;              // rounded-up alignsize for file I/O
      uint64_t fileoffset;            // summed filesize of all previous pages
      int nkey;                       // # of KV pairs
    };

    Page *pages;                      // list of pages
    int maxpage;                      // max # of pages currently allocated

    // file info

    char *filename;                   // filename to store KV if needed
    FILE *fp;                         // file ptr
    int fileflag;                     // 1 if file exists, 0 if not

    // private methods

    void add(KeyValue *);
    void add(int, char *);
    void add(char *);
    void add(int, char *, uint64_t, uint64_t, uint64_t);
    void add(int, char *, int, int);

    void init_page();
    void create_page();
    void write_page();
    void read_page(int, int);
    uint64_t roundup(uint64_t, int);
};

  }  // namespace MAPREDUCE_NS

//#endif

#endif  // KEYVALUE_H_
