/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.

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


#define FORMAT_PHYLIP 		0
#define FORMAT_FASTA 		1
#define FORMAT_NEXUS 		2
#define FORMAT_MSF 			3
#define FORMAT_CLUSTAL 		4

typedef enum { LVB_FALSE, LVB_TRUE } Lvb_bool;	/* boolean type */
#define LVB_FNAMSIZE 2000		/* maximum bytes for file names */

#define MAX_BOOTSTRAPS 1000000	/* max. bootstrap replicates */


/* matrix and associated information */
typedef struct data
{
    int n_threads_getplen;  /* number of possible threads in getplen function */
    int n_slice_size_getplen;  /* slice size in getplen function, usually m/n_threads_getplen  */
    long m;				/* number of columns */
    long original_m;	/* number of columns read from matrix*/
    long n;				/* number of rows */
    long nbranches; 	/* number of possible braches */
    long bytes;
    long tree_bytes;	/* length the tree in bytes */
    long tree_bytes_whitout_sset;	/* length the tree in bytes whitout sset */
    long nwords;
    long nsets;			/* sets per tree */
    long mssz;			/* maximum objects per set */
    char **row;			/* array of row strings */
    char **rowtitle;	/* array of row title strings */
} *Dataptr, DataStructure;

/* unchangeable types */

/* user- or programmer-configurable parameters */
typedef struct
{
    int seed;							/* seed for random number generator */
    int cooling_schedule;   			/* cooling schedule: 0 is geometric, 1 is linear */
 /* int algorithm_selection */       /* algorithm selection: 0 is TBR, 1 is no TBR */
    int n_file_format;					/* number of file format, must be FORMAT_PHYLIP, FORMAT_FASTA, FORMAT_NEXUS, FORMAT_MSF, FORMAT_CLUSTAL*/
    int n_processors_available;			/* number of processors available */
    long verbose;						/* verboseness level */
    long bootstraps;					/* number of bootstrap replicates */
    int n_number_max_trees;				/* number of bootstrap replicates */
    char file_name_in[LVB_FNAMSIZE];	/* input file name */
    char file_name_out[LVB_FNAMSIZE];	/* output file name */
} Params;
