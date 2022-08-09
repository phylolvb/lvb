/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.

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

/* ========== LogFile.c - logfile parameters ========== */

#include "LogFile.h"

bool LogFileExists(const char *filename) {
  struct stat buffer;
  
  return (stat (filename, &buffer) == 0);
}

void PrintLogFile(long iter, long trees_output_total, long final_length, double overall_time_taken) {
  Lvb_bool logfile_exists = LVB_FALSE;

  if (LogFileExists ("logfile.tsv")) {
    logfile_exists = LVB_TRUE;
  }
  
  FILE * logfile;
  logfile = fopen ("logfile.tsv","a+");
  
  if (logfile_exists == LVB_FALSE)
	{
    fprintf (logfile, "Implementation\tRearrangements\tTopologies\tScore\tRuntime\n");
  }
		fprintf (logfile, "%s\t%ld\t%ld\t%ld\t%.2lf\n", LVB_IMPLEMENTATION, iter, trees_output_total, final_length, overall_time_taken);
		fclose(logfile);
}

#ifdef LVB_MPI
  void PrintMPILogFile(long iter, long trees_output_total, long final_length, double overall_time_taken, int seedMPI){
  Lvb_bool logfile_exists = LVB_FALSE;

  if (LogFileExists ("logfile.tsv")) {
    logfile_exists = LVB_TRUE;
  }
  
  FILE * logfile;
  logfile = fopen ("logfile.tsv","a+");
  
  if (logfile_exists == LVB_FALSE){
    fprintf (logfile, "Implementation\tSeed\tRearrangements\tTopologies\tScore\tRuntime\n");
  }
		fprintf (logfile, "%s\t%d\t%ld\t%ld\t%ld\t%.2lf\n", LVB_IMPLEMENTATION, seedMPI, iter, trees_output_total, final_length, overall_time_taken);
		fclose(logfile);
  }
#endif