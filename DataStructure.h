/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
All rights reserved.

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

/* ========== DataStructure.h - definition of data structures ========== */

#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#define FORMAT_PHYLIP 		0
#define FORMAT_FASTA 		1
#define FORMAT_NEXUS 		2

#ifndef NP_Implementation
#define FORMAT_MSF 		3
#define FORMAT_CLUSTAL 		4
#else
#define FORMAT_CLUSTAL          3

#define MAX_BOOTSTRAPS 1000000 // max bootstrap replicates
#endif

#define LVB_FNAMSIZE 2000		/* maximum bytes for file names */

/* these flags is to read and save states in specfic time points */
#define DONT_SAVE_READ_STATES		        0		/* dont read and save states, default parameter */
#define DO_SAVE_READ_STATES			1		/* try to read and save states */

#define CHECK_POINT_PROCESS_FINISHED		1		/* process is finished, don't need to run again */
#define CHECK_POINT_PROCESS_NOT_FINISHED	0		/* process not finished yet, need to load states and run */

#define CHECK_POINT_READ_STATE_FILES		1		/* the state files exist and are OK */
#define CHECK_POINT_NOT_READ_STATE_FILES	0		/* the state files don't exist and are corrupted */
/* END save read states flags */

typedef enum { LVB_FALSE, LVB_TRUE } Lvb_bool;	/* boolean type */

typedef struct data
{
     long m;				/* number of columns */
     long original_m;	                /* number of columns read from matrix*/
     long n;				/* number of rows */
     
     long nbranches; 	                /* number of possible braches */
     long bytes;
     long tree_bytes;	                /* length the tree in bytes */
     
     long nwords;
     long nsets;	                /* sets per tree */
     long mssz;	                        /* maximum objects per set */
     int n_threads_getplen;  	        /* number of possible threads in getplen function */
     int n_slice_size_getplen;          /* slice size in getplen function, usually m/n_threads_getplen  */

     #ifndef NP_Implementation
     long max_length_seq_name; 	        /* length of the sequence names */
     long tree_bytes_without_sset;	/* length the tree in bytes without sset */
     long min_len_tree;                 // minimum length of tree given matrix
     #else
     long tree_bytes_without_sset;	/* length the tree in bytes without sset */
     char **row;                        // array of row strings
     char **rowtitle;                   // array of row title strings
     #endif
} *Dataptr, DataStructure;

#ifndef NP_Implementation
typedef struct seq_data
{
     char **row;
     char **rowtitle;
}*DataSeqPtr, DataSeqStructure;
#endif

/* user- or programmer-configurable parameters */
typedef struct
{
     long verbose;		                /* verboseness level */
     int seed;			                /* seed for random number generator */
     int cooling_schedule;                      /* cooling schedule: 0 is geometric, 1 is linear */
     int n_file_format;		                /* number of file format, must be FORMAT_PHYLIP, FORMAT_FASTA, FORMAT_NEXUS, FORMAT_MSF, FORMAT_CLUSTAL*/
     int n_processors_available;	        /* number of processors available */

      #ifndef NP_Implementation
      int n_seeds_need_to_try;	                /* number of seeds that go to try, minimum is the number of mpi process */
     int n_flag_save_read_states;		/* flag to save/read the states, if 1 when starts try to read the last states, if not find */
							/* start from the begin and is going to save the states each pre-defined time schedule */
							/* it can go to 0 when file is zero or corrupted */
      int n_flag_is_finished_process;		/* is set to one when the process is finished  */
							/* 	set to zero is necessary to start from a specific state */
      int n_flag_is_possible_read_state_files;	/* if is possible to read the states or not. */
      int n_checkpoint_interval;		/* value in seconds when a checkpoint file is going to be saved, default(CHECKPOINT_INTERVAL)*/
      int n_make_test;				/* it is only used for tests */      
      #else
      int algorithm_selection;                  // algorithm selection: 0 = SPR+NNI, 2 = SPR+NNI+TBR, 3 = SPR+NNI+TBR
      long bootstraps;                          // number of bootstrap replicates
      long n_number_max_trees;                  //maximum number of trees saved?
      #endif
      char file_name_in[LVB_FNAMSIZE];	        /* input file name */
      char file_name_out[LVB_FNAMSIZE];	        /* output file name */
} Params;

#ifndef NP_Implementation
/* structure to use sending temperature and number of iterations to master process */
typedef struct
{
        int n_iterations;		/* number of iterations */
        int n_seed;			/* seed for this temperature and iteration */
        long l_length;			/* length of the tree */
        double temperature;		/* temperature */
} SendInfoToMaster;

/* structure to use sending if is to continue and new seed if it is */
typedef struct
{
        int n_seed;				/* new seed to start the process again */
        int n_is_to_continue;	                /* if it is to start the process */
        int n_process_tried;	                /* id of the seed tried */
} RecvInfoFromMaster;

/* structures for calculation of averages temperatures and std */
typedef struct IndividualTemperature
{
        double d_temperature;	                                /* temperature */
        int n_try_process;		                        /* number of process tried, is sequential...*/
        int n_seed;				                /* seed for this temperature and iteration */
							                /* the seed is here because it's easier to perform the algorithm */
        long l_length;			                        /* length of the tree */
        struct IndividualTemperature *p_next_temperature;       /* next memory structure */
}IndividualTemperature;

typedef struct IterationTemperature
{
        int n_iteration;
        IndividualTemperature *p_temperature;
        struct IterationTemperature *p_next_iteration;
}IterationTemperature;

#endif

#endif /* DATASTRUCTURE_H */
