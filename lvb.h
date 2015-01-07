/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
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

/* ********** lvb.h - main header for lvb ********** */

#ifndef LVB_LVB_H
#define LVB_LVB_H

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "myuni.h"
#include "mymaths.h"
#include "DataStructure.h"

/* the program */
#define PROGNAM "lvb"			/* program file name */
#define LVB_VERSION "IN DEVELOPMENT"	/* version of program */
#define LVB_SUBVERSION "(2014)"		/* version details e.g. date */

/* set if is to compile with 64 or 32 */
#define COMPILE_64_BITS				/* the default is 32 bits */
//#define	MPI_SEND_ONLY_MATRIX_NAMES	/* if defined only send the names of the matrix
									   /* sometimes the data matrix are huge and it's only necessary to pass
    								   the names of the other process */

/* DNA bases: bits to set in statesets */
#define A_BIT 0b0001		/* (1U << 0) */
#define C_BIT 0b0010		/* (1U << 1) */
#define G_BIT 0b0100		/* (1U << 2) */
#define T_BIT 0b1000		/* (1U << 3) */

#define NIBBLE_MASK 		017			/* space for one stateset, all bits set to 1 */
#define NIBBLE_WIDTH 		4			/* width of nibble in bits */
#define NIBBLE_WIDTH_BITS	2			/* bitwise multiply the NIBBLE_WIDTH */

#ifdef COMPILE_64_BITS
	typedef uint64_t Lvb_bit_lentgh;								/* define 64 bits */
	#define NUMBER_OF_BITS										64
	#define LENGTH_WORD											16	/* length of number packed bases */
	#define LENGTH_WORD_BITS_MULTIPLY							4	/* multiply of number packed bases */
	#define MINIMUM_WORDS_PER_SLICE_GETPLEN						30  /* minimum words per slice that run gplen threading */
	#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		60 /* need to have this size to activate the threading */
	#define MASK_SEVEN											0x7777777777777777U
	#define MASK_EIGHT											0x8888888888888888U
#else		/* default 32 bits */
	typedef uint32_t Lvb_bit_lentgh;								/* define 32 bits */
	#define NUMBER_OF_BITS										32
	#define LENGTH_WORD											8	/* length of number packed bases */
	#define LENGTH_WORD_BITS_MULTIPLY							3	/* multiply of number packed bases */
	#define MINIMUM_WORDS_PER_SLICE_GETPLEN						30	/* minimum words per slice that run gplen threading */
	#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		60 /* need to have this size to activate the threading */
	#define MASK_SEVEN											0x77777777U
	#define MASK_EIGHT											0x88888888U
#endif

/* values some people may feel the dangerous urge to change */
#define LVB_INPUTSTRING_SIZE 2000	/* max. bytes for interactive input */
#define UNSET (-1)			/* value of integral vars when unset */
#define STAT_LOG_INTERVAL 50000	/* min. interval for progress log */
#define REROOT_INTERVAL 1000		/* change root every ... updates */


/* limits that could be changed but, if increased enormously, might lead to
 * some trouble at some point */
#define MAX_N 1000000	/* max. rows */
#define MAX_M 5000000	/* max. cols */

/* implementation-independent limits */
#define LVB_EPS 1E-8		/* 0.0 < DBL_EPSILON < LVB_EPS */
#define MIN_M 1L		/* min. no. of characters for any analysis */
#define MAX_BRANCHES (2 * MAX_N - 3)	/* max. branches per tree */
#define MIN_BRANCHES (2 * MIN_N - 3)	/* max. branches per tree */
#define MIN_N 5L		/* min. no. of objs, for rearrangeable tree */
#define MAX_ALLOC ((size_t) (INT_MAX - 2))	/* max. bytes per dyn. alloc. */
#define MAXSTATES 5		/* max. "true" states in data matrix */

/* limits that could be changed but are likely to be OK */
#define FROZEN_T 0.0001		/* consider system frozen if temp < FROZEN_T */


/* MPI definitions... */
#define MPI_MAIN_PROCESS	0		/* main process */

#define	MPI_TAG_MATRIX					1
#define	MPI_TAG_NAME_AND_SEQ_DATA		2
#define	MPI_TAG_BINARY_DATA				3
#define MPI_TAG_PARAMS					4
/* END MPI definitions... */



/* branch of tree */
typedef struct
{
    long parent;		/* parent branch number, UNSET in root */
    long left;			/* index of first child in tree array */
    long right;			/* index of second child in tree array */
    long changes;		/* changes associated with this branch */
    Lvb_bit_lentgh *sset;	/* statesets for all sites */
} Branch;

/* tree stacks */
typedef struct
{
    Branch *tree;	/* pointer to first branch in tree array */
    long root;		/* root of tree */
} Treestack_element;

typedef struct
{
    long size;			/* number of trees currently allocated for */
    long next;			/* next unused element of stack */
    Treestack_element *stack;	/* pointer to first element in stack */
} Treestack;


/* simulated annealing parameters */
#define MAXACCEPT_SLOW 5L	/* maxaccept for "slow" searches */
#define MAXPROPOSE_SLOW 2000L	/* maxpropose for "slow" searches */
#define MAXFAIL_SLOW 40L	/* maxfail for "slow" searches */

/* fixed file names */
#define MATFNAM "infile"	/* matrix file name */
#define OUTTREEFNAM "outtree"	/* overall best trees */

/* verbose-mode file name bases (run-specific suffixes will be used) */
#define LENFNAM "stat"		/* current tree and length file name prefix */
#define RESFNAM "res"		/* cycle's results file name prefix */
#define SUMFNAM "sum"		/* summary of trees per run file name */
#define TREE1FNAM "ini"		/* cycle's initial tree file name prefix */


/* assert-like macro, differing in that it writes to standard output,
 * calls crash() not abort(), and works whether or not NDEBUG is defined */
#define lvb_assert(test) ((void) ((test) || (lvb_assertion_fail(#test, __FILE__, __LINE__), 0)))

/* PHYLIP global data */
//extern long chars;	/* defined in dnapars.c */

/* LVB global functions */
void *alloc(const size_t, const char *const);
long anneal(Dataptr restrict, Treestack *, const Branch *const, Params rcstruct, long, const double,
 const long, const long, const long, FILE *const, const long *, long *, int, Lvb_bool);
long arbreroot(Dataptr, Branch *const, const long);
long bytes_per_row(const long);
long childadd(Branch *const, const long, const long);
long cistrcmp(const char *const, const char *const);
Lvb_bool cleanup(void);
void clnclose(FILE *const, const char *const);
FILE *clnopen(const char *const, const char *const);
void clnremove(const char *const);
void crash(const char *const, ...);
void defaults_params(Params *const prms);
long deterministic_hillclimb(Dataptr, Treestack *, const Branch *const, Params rcstruct,
	long, FILE * const, const long *, long *, int myMPIid, Lvb_bool);
void dna_makebin(Dataptr restrict, DataSeqPtr matrix_seq, Lvb_bit_lentgh **);
void dnapars_wrapper(void);
char *f2str(FILE *const);
Lvb_bool file_exists(const char *const);
void get_bootstrap_weights(long *, long, long);
double get_initial_t(Dataptr, const Branch *const, Params rcstruct, long, const long *, Lvb_bool);
int getparam(Params *, int argc, char **argv);

long getplen(Dataptr restrict, Branch *, Params rcstruct, const long, const long *restrict, long *restrict p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs);

void alloc_memory_to_getplen(Dataptr matrix, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
double get_predicted_length(double, double, long, long, long, long);
double get_predicted_trees(double, double, long, long, long, long);
long getroot(const Branch *const);
void lvb_assertion_fail(const char *, const char *, int);
void lvb_initialize(void);
Dataptr lvb_matrin(const char *);
long lvb_reroot(Dataptr restrict, Branch *const barray, const long oldroot, const long newroot, Lvb_bool b_with_sset);
void lvb_treeprint (Dataptr, DataSeqPtr restrict matrix_seq_data, FILE *const, const Branch *const, const long);
void matchange(Dataptr, DataSeqPtr, const Params);
Dataptr matrin(const char *const);
void mutate_deterministic(Dataptr restrict, Branch *const, const Branch *const, long, long, Lvb_bool);
void mutate_spr(Dataptr restrict, Branch *const, const Branch *const, long);
void mutate_nni(Dataptr restrict, Branch *const, const Branch *const, long);
char *nextnonwspc(const char *);
void nodeclear(Branch *const, const long);
long objreroot(Branch *const, const long, const long);
void params_change(Params *);
int phylip_dna_matrin(char *, int, Dataptr, DataSeqPtr);
void phylip_mat_dims_in(char *, int, long *, long *, int *);
void randtree(Dataptr, Branch *const);
long randpint(const long);
void rowfree(DataSeqPtr, int n_lines);
void scream(const char *const, ...);
void ss_init(Dataptr, Branch *, Lvb_bit_lentgh **);
char *supper(char *const s);
Branch *treealloc(Dataptr restrict, Lvb_bool b_with_sset);
long tree_bytes(Dataptr restrict matrix);
long tree_bytes_whitout_sset(Dataptr restrict matrix);
void treeclear(Dataptr, Branch *const);
void treecopy(Dataptr restrict, Branch *const, const Branch *const, Lvb_bool b_with_sset);
long treecmp(Dataptr, const Branch *const, const long, const Branch *const, long);
void treedump(Dataptr, FILE *const, const Branch *const, Lvb_bool b_with_sset);
void treedump_screen(Dataptr matrix, const Branch *const tree);
void treestack_clear(Treestack *);
long treestack_cnt(Treestack);
long treestack_dump(Dataptr, Treestack *, FILE *const);
void treestack_free(Treestack *);
Treestack treestack_new(void);
long treestack_transfer(Dataptr, Treestack *, Treestack *);
long treestack_pop(Dataptr, Branch *, long *, Treestack *);
long treestack_print(Dataptr, DataSeqPtr restrict matrix_seq_data, Treestack *, FILE *const, Lvb_bool);
long treestack_push(Dataptr, Treestack *, const Branch *const, const long);
void treeswap(Branch **const, long *const, Branch **const, long *const);
void uint32_dump(FILE *, Lvb_bit_lentgh);
long words_per_row(const long);

#endif /* LVB_LVB_H */
