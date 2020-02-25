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

/* ========== fops.c - file operations ========== */

#include "lvb.h"

Lvb_bool file_exists(const char *const nam) {
  Lvb_bool val;  // return value
  FILE *fp = fopen(nam, "r");  // store file ptr so can close it

  if (fp == NULL) {
    val = LVB_FALSE;
  } else {
    fclose(fp);
    val = LVB_TRUE;
  }

  return val;
}  // end

FILE *clnopen(const char *const nam, const char *const mod) {
  FILE *fp;  // file

  fp = fopen(nam, mod);
  if (fp == NULL) {
    if (strcmp(mod, "w") == 0)
      crash("cannot create file '%s'", nam);
    else if (strcmp(mod, "r") == 0)
      crash("cannot open file '%s' for reading", nam);
    else if (strcmp(mod, "a") == 0)
      crash("cannot open file '%s' for appending to", nam);
    else  // rare mode
      crash("cannot open file '%s' with mode '%s'", nam, mod);
    }

    return fp;
}  // end clnopen()

void clnclose(FILE *const fp, const char *const fnam) {
  if (fp != NULL) {
    if (ferror(fp) != 0)
      crash("file error on file '%s'", fnam);
    if (fclose(fp) != 0)
      crash("cannot close file '%s'", fnam);
    }
}  // end clnclose()

void clnremove(const char *const fnam) {
  if (remove(fnam) != 0)
  scream("cannot delete file '%s'", fnam);
}  // end clnremove

char *f2str(FILE *const stream) {
  char *input;                                // input string
  unsigned long inbytes;                      // bytes for string
  unsigned long off;                          // position in string
  unsigned long offmax = 0UL;                 // length of string
  const unsigned long maxom = ULONG_MAX - 3;  // maximum initial offmax

  /* calculate file size and allocate appropriately */
  while (getc(stream) != EOF) {
    offmax++;
    if (offmax >= maxom)  // crash while value has meaning
      crash("input is too long");
  }
  if (ferror(stream) != 0)
    crash("file error on reading file");
  inbytes = offmax + 2UL;  // '\0', possible '\n'
  input = (char *) alloc(inbytes, "input");

  // get string
  rewind(stream);
  for (off = 0; off < offmax; off++)
    input[off] = (char) getc(stream);
  if (ferror(stream) != 0)
    crash("file error on reading file");

  // terminate string, also adding newline at end if not present
  if (input[off-1] != '\n') {
    input[off] = '\n';
    off++;
  }
  input[off] = '\0';

  return input;
}  // end f2str()
