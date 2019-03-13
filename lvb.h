/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== lvb.h - main header for lvb ========== */

//#NP_Implementation

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
#include "myuni.h"
#include "mymaths.h"
#include "DataStructure.h"
#include <time.h>

/* the program */
#define PROGNAM "lvb"			/* program file name */
#define LVB_VERSION "3.5"		/* version of program */
#define LVB_SUBVERSION "(February 2019)"	/* version details e.g. date */

/* set if is to compile with 64 or 32 */
#ifndef COMPILE_32_BITS
	#define COMPILE_64_BITS				/* the default is 32 bits */
#endif

/* use TBR branch-Swapping Algorithm? */
/* #define TBR */ 

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
	#define MINIMUM_WORDS_PER_SLICE_GETPLEN						20  /* minimum words per slice that run gplen threading */
	#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		40 /* need to have this size to activate the threading */
	#define MASK_SEVEN											0x7777777777777777U
	#define MASK_EIGHT											0x8888888888888888U

	/* limits that could be changed but, if increased enormously, might lead to
	 * some trouble at some point */
	#define MAX_N 1000000		/* max. rows */
	#define MAX_M 5000000		/* max. cols */
#else		/* default 32 bits */
	typedef uint32_t Lvb_bit_lentgh;								/* define 32 bits */
	#define NUMBER_OF_BITS										32
	#define LENGTH_WORD											8	/* length of number packed bases */
	#define LENGTH_WORD_BITS_MULTIPLY							3	/* multiply of number packed bases */
	#define MINIMUM_WORDS_PER_SLICE_GETPLEN						20	/* minimum words per slice that run gplen threading */
	#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		40 /* need to have this size to activate the threading */
	#define MASK_SEVEN											0x77777777U
	#define MASK_EIGHT											0x88888888U

	/* limits that could be changed but, if increased enormously, might lead to
		 * some trouble at some point */
	#define MAX_N 15000		/* max. rows */
	#define MAX_M 30000		/* max. cols */

#endif

/* values some people may feel the dangerous urge to change */
#define LVB_INPUTSTRING_SIZE 2000	/* max. bytes for interactive input */
#define STAT_LOG_INTERVAL 50000		/* min. interval for progress log */
#define REROOT_INTERVAL 1000		/* change root every ... updates */


/* implementation-independent limits */
#define LVB_EPS 1E-8					/* 0.0 < DBL_EPSILON < LVB_EPS */
#define MIN_M 1L						/* min. no. of characters for any analysis */
#define UNSET (-1)						/* value of integral vars when unset */
#define MAX_BRANCHES (2 * MAX_N - 3)	/* max. branches per tree */
#define MIN_BRANCHES (2 * MIN_N - 3)	/* max. branches per tree */
#define MIN_N 5L						/* min. no. of objs, for rearrangeable tree */
#define MAX_ALLOC ((size_t) (INT_MAX - 2))	/* max. bytes per dyn. alloc. */
#define MAXSTATES 5						/* max. "true" states in data matrix */

/* limits that could be changed but are likely to be OK */
#define FROZEN_T 0.0001		/* consider system frozen if temp < FROZEN_T */


typedef	struct	/* object set derived from a cladogram */
{
	int *set;	/* arrays of object sets */
	int cnt;	/* sizes of object sets */
}	Objset;

/* branch of tree */
typedef struct
{
    int parent;		/* parent branch number, UNSET in root */
    int left;			/* index of first child in tree array */
    int right;			/* index of second child in tree array */
    int changes;		/* changes associated with this branch */
    Lvb_bit_lentgh *sset;	/* statesets for all sites */
} Branch;

/* tree stacks */
typedef struct
{
	int root;		/* root of tree */
    Branch *tree;	/* pointer to first branch in tree array */
    Objset *p_sset;	/* array with sset with the root always on zero */
} Treestack_element;

typedef struct
{
	int size;			/* number of trees currently allocated for */
	int next;			/* next unused element of stack */
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
long anneal(Dataptr restrict, Treestack *, Treestack *, const Branch *const, Params rcstruct, long, const double,
 const long, const long, const long, FILE *const, const long *, long *, Lvb_bool);
int arbreroot(Dataptr, Branch *const, const int);
long bytes_per_row(const long);
int childadd(Branch *const, const int, const int);
long cistrcmp(const char *const, const char *const);
Lvb_bool cleanup(void);
void clnclose(FILE *const, const char *const);
FILE *clnopen(const char *const, const char *const);
void clnremove(const char *const);
void copy_sset(Dataptr restrict matrix, Objset *p_sset_1);
void crash(const char *const, ...);
void defaults_params(Params *const prms);
long deterministic_hillclimb(Dataptr, Treestack *, const Branch *const, Params rcstruct,
	long, FILE * const, const long *, long *, Lvb_bool);
void dna_makebin(Dataptr restrict, Lvb_bit_lentgh **);
void dnapars_wrapper(void);
void dump_stack_to_screen(Dataptr matrix, Treestack *sp);
void dump_objset_to_screen(Dataptr matrix, Objset *oset_1);
void dump_objset_to_screen_sset_2(Dataptr matrix);
char *f2str(FILE *const);
Lvb_bool file_exists(const char *const);
void get_bootstrap_weights(long *, long, long);
double get_initial_t(Dataptr, const Branch *const, Params rcstruct, long, const long *, Lvb_bool);
long getminlen(const Dataptr);
void getparam(Params *, int argc, char **argv);
long getplen(Dataptr restrict, Branch *, Params rcstruct, const long, const long *restrict, long *restrict p_todo_arr,
			long *p_todo_arr_sum_changes, int *p_runs);
void alloc_memory_to_getplen(Dataptr matrix, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
double get_predicted_length(double, double, long, long, long, long);
double get_predicted_trees(double, double, long, long, long, long);
long getroot(const Branch *const);
void lvb_assertion_fail(const char *, const char *, int);
void lvb_initialize(void);
Dataptr lvb_matrin(const char *);
long lvb_reroot(Dataptr restrict, Branch *const barray, const int oldroot, const int newroot, Lvb_bool b_with_sset);
void lvb_treeprint (Dataptr, FILE *const, const Branch *const, const long);
void makesets(Dataptr restrict, const Branch *const tree_2, const int root);
void matchange(Dataptr, const Params);
Dataptr matrin(const char *const);
void mutate_deterministic(Dataptr restrict, Branch *const, const Branch *const, long, long, Lvb_bool);
void mutate_spr(Dataptr restrict, Branch *const, const Branch *const, int);
void mutate_nni(Dataptr restrict, Branch *const, const Branch *const, int);
void mutate_tbr(Dataptr restrict, Branch *const, const Branch *const, int);
char *nextnonwspc(const char *);
void nodeclear(Branch *const, const long);
long objreroot(Branch *const, const long, const long);
void params_change(Params *);
void phylip_dna_matrin(char *, int, Dataptr);
void phylip_mat_dims_in(char *, int, long *, long *, int *);
void randtree(Dataptr, Branch *const);
int randpint(const long);
void rowfree(Dataptr);
void scream(const char *const, ...);
long setstcmp_with_sset2(Dataptr matrix, Objset *const oset_1);
void ss_init(Dataptr, Branch *, Lvb_bit_lentgh **);
char *supper(char *const s);
Branch *treealloc(Dataptr restrict, Lvb_bool b_with_sset);
long tree_bytes(Dataptr restrict matrix);
long tree_bytes_without_sset(Dataptr restrict matrix);
void treeclear(Dataptr, Branch *const);
void treecopy(Dataptr restrict, Branch *const, const Branch *const, Lvb_bool b_with_sset);
long treecmp(Dataptr restrict, Objset *, const Branch *const, Lvb_bool b_first);
void treedump(Dataptr, FILE *const, const Branch *const, Lvb_bool b_with_sset);
void treedump_screen(Dataptr matrix, const Branch *const tree);
void treestack_clear(Treestack *);
long treestack_cnt(Treestack);
long treestack_dump(Dataptr, Treestack *, FILE *const);
void treestack_free(Dataptr restrict matrix, Treestack *);
Treestack treestack_new(void);
long treestack_transfer(Dataptr, Treestack *, Treestack *, Lvb_bool b_with_sset);
long treestack_pop(Dataptr, Branch *, long *, Treestack *, Lvb_bool b_with_sset);
int treestack_print(Dataptr, Treestack *, FILE *const, Lvb_bool onerandom);
long treestack_push(Dataptr, Treestack *, const Branch *const, const int, Lvb_bool b_with_sset);
void treeswap(Branch **const, long *const, Branch **const, long *const);
void uint32_dump(FILE *, Lvb_bit_lentgh);
long words_per_row(const long);
int count(Branch *const, int);
int addtoarray(Branch *const, int, int *, int);

#endif /* LVB_LVB_H */

//#endif

#ifdef MPI_Implementation

#ifndef LVB_LVB_H
#define LVB_LVB_H

#include "DataStructure.h"
#include <mpi.h>
#include "mapreduce.h"
#include "blockmacros.h"
#include "keyvalue.h"

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
#include "myuni.h"
#include "mymaths.h"

#ifdef MAP_Reduce
	#include <omp.h>
	#include "sys/stat.h"

	using namespace MAPREDUCE_NS;
	using namespace std;
	#define __STDC_LIMIT_MACROS
#endif

/* the program */
#define PROGNAM "lvb"			/* program file name */
#define LVB_VERSION "IN DEVELOPMENT"	/* version of program */
#define LVB_SUBVERSION "$Id: 2f1a0126358acf0a108b58851f7a368fcfe57f7f $"		/* version details e.g. date */


/* set if is to compile with 64 or 32 bits */
#ifndef COMPILE_32_BITS
	#define COMPILE_64_BITS
#endif
#define	MPI_SEND_ONLY_MATRIX_NAMES	/* if defined only send the names of the matrix */
									/* sometimes the data matrix are huge and it's only necessary to pass */
    								/* the names of the other process */


/* DNA bases: bits to set in statesets */
#define A_BIT 0b0001		/* (1U << 0) */
#define C_BIT 0b0010		/* (1U << 1) */
#define G_BIT 0b0100		/* (1U << 2) */
#define T_BIT 0b1000		/* (1U << 3) */

#define NIBBLE_MASK 		017			/* space for one stateset, all bits set to 1 */
#define NIBBLE_WIDTH 		4			/* width of nibble in bits */
#define NIBBLE_WIDTH_BITS	2			/* bitwise multiply the NIBBLE_WIDTH */

#ifdef COMPILE_64_BITS
	typedef uint64_t Lvb_bit_lentgh;							/* define 64 bits */
	#define NUMBER_OF_BITS										64
	#define LENGTH_WORD											16			/* length of number packed bases */
	#define LENGTH_WORD_BITS_MULTIPLY							4			/* multiply of number packed bases */
	#define MINIMUM_WORDS_PER_SLICE_GETPLEN						30  			/* minimum words per slice that run gplen threading */
	#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		60 			/* need to have this size to activate the threading */
	#define MASK_SEVEN											0x7777777777777777U
	#define MASK_EIGHT											0x8888888888888888U
#else		/* default 32 bits */
	typedef uint32_t Lvb_bit_lentgh;							/* define 32 bits */
	#define NUMBER_OF_BITS										32
	#define LENGTH_WORD											8			/* length of number packed bases */
	#define LENGTH_WORD_BITS_MULTIPLY							3			/* multiply of number packed bases */
	#define MINIMUM_WORDS_PER_SLICE_GETPLEN						30			/* minimum words per slice that run gplen threading */
	#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		60 			/* need to have this size to activate the threading */
	#define MASK_SEVEN											0x77777777U
	#define MASK_EIGHT											0x88888888U
#endif

/* values some people may feel the dangerous urge to change */
#define LVB_INPUTSTRING_SIZE 2000	/* max. bytes for interactive input */
#define UNSET (-1)			/* value of integral vars when unset */
#define CHECKPOINT_INTERVAL 1800	/* checkpoint ~every ... seconds */
#define CHECKPOINT_FNAM_BASE "lvb_checkpoint"
#define STAT_LOG_INTERVAL 50000	/* min. interval for progress log */
#define REROOT_INTERVAL 1000		/* change root every ... updates */


/* limits that could be changed but, if increased enormously, might lead to
 * some trouble at some point */
#define MAX_N 1000000	/* max. rows */
#define MAX_M 5000000	/* max. cols */

/* implementation-independent limits */
#define LVB_EPS 1E-11		/* 0.0 < DBL_EPSILON < LVB_EPS */
#define MIN_M 1L		/* min. no. of characters for any analysis */
#define MAX_BRANCHES (2 * MAX_N - 3)	/* max. branches per tree */
#define MIN_BRANCHES (2 * MIN_N - 3)	/* max. branches per tree */
#define MIN_N 5L		/* min. no. of objs, for rearrangeable tree */
#define MAX_ALLOC ((size_t) (INT_MAX - 2))	/* max. bytes per dyn. alloc. */
#define MAXSTATES 5		/* max. "true" states in data matrix */

/* limits that could be changed but are likely to be OK */
#define INITIAL_INCREMENT 0.00001	/* step size to get initial temp. */
#define FROZEN_T 0.0001		/* helps to decide whether system is frozen */


/* MPI definitions... */
#define MPI_MAIN_PROCESS	0		/* main process */

#define	MPI_TAG_MATRIX					1
#define	MPI_TAG_NAME_AND_SEQ_DATA		2
#define	MPI_TAG_BINARY_DATA				3
#define MPI_TAG_PARAMS					4
#define MPI_TAG_SEND_TEMP_MASTER		5
#define MPI_TAG_SEND_FINISHED			6
#define MPI_TAG_SEND_RESTART			7
#define MPI_TAG_MANAGEMENT_MASTER		8


#define MPI_FINISHED							0x00
#define MPI_IS_TO_RESTART_ANNEAL				0x01
#define MPI_IS_TO_CONTINUE_ANNEAL				0x02
#define MPI_IS_NOT_TO_RESTART					0x03
#define MPI_IS_TO_CONTINUE						0x04

/* END MPI definitions... */

/* anneal state */
#define MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN			0x00
#define MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT				0x01
#define MESSAGE_ANNEAL_FINISHED_AND_REPEAT					0x02
#define MESSAGE_ANNEAL_KILLED								0x03
#define MESSAGE_ANNEAL_KILLED_AND_REPEAT					0x04
#define MESSAGE_ANNEAL_KILLED_AND_NOT_REPEAT				0x05
#define MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE		0x06
#define MESSAGE_ANNEAL_STOP_PROCESS							0x07
#define MESSAGE_BEGIN_CONTROL								0x08
/* anneal state */


/* Define calc iterations algorithm */
#define CALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS		5	/* the number of chunk of iterations that need to be */
																	/* to start release */
																	/* Ex: if 3 the algorithm doesn't kill the process in the first three chunk of temperatures */
#define CALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS			1	/* if for a specific process the length exceeds this many SD from the mean */
																	/* then the process need to restart with other seed */
/* END  Define calc iterations algorithm */


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
#define MAXACCEPT_MIN 5L		/* minimum value for maxaccept */
#ifdef MAP_Reduce
#define MAXACCEPT_MAX 5L		/* maximum value for maxaccept */
#else
#define MAXACCEPT_MAX 500L		/* maximum value for maxaccept */
#endif	/* MAP_REDUCE_SINGLE */
#define MAXPROPOSE_SLOW 2000L	/* maxpropose for "slow" searches */
#define MAXFAIL_SLOW 40L		/* maxfail for "slow" searches */
#define GRAD_GEOM 0.99	/* for relationship between t(n) and t(n+1) */
#define GRAD_LINEAR (3.64 * 1E-8)	/* temperature gradient */

/* fixed file names */
#define MATFNAM 				"infile"	/* matrix file name */
#define OUTTREEFNAM 			"outtree"	/* overall best trees */

/* verbose-mode file name bases (run-specific suffixes will be used) */
#define LENFNAM 				"stat"		/* current tree and length file name prefix */
#define RESFNAM 				"res"		/* cycle's results file name prefix */
#define SUMFNAM 				"sum"		/* summary of trees per run file name */
#define TREE1FNAM 				"ini"		/* cycle's initial tree file name prefix */


/* assert-like macro, differing in that it writes to standard output,
 * calls crash() not abort(), and works whether or not NDEBUG is defined */
#define lvb_assert(test) ((void) ((test) || (lvb_assertion_fail(#test, __FILE__, __LINE__), 0)))

/* PHYLIP global data */
//extern long chars;	/* defined in dnapars.c */

#ifdef MAP_Reduce

	struct MISC {
		int rank,nprocs;

		int ID;
		long num;
		bool SB;

		int ntrees;
		long nsets;
		long mssz;

		int *count;
	};


	#ifdef __cplusplus
		#define restrict    /* nothing */
	#endif

#endif

/* LVB global functions */
void *alloc(const size_t, const char *const);

#ifdef MAP_Reduce
	long anneal(Dataptr restrict, Treestack *, const Branch *const, Params *p_rcstruct, long, const double,
		const long, const long, const long, FILE *const, long *, Lvb_bool, MISC *misc, MapReduce *mrStackTree, MapReduce *mrBuffer);
#else
	long anneal(Dataptr restrict, Treestack *, const Branch *const, Params *p_rcstruct, long, const double,
		const long, const long, const long, FILE *const, long *, int, int *p_n_state_progress, int *p_n_number_tried_seed, Lvb_bool);
#endif

long arbreroot(Dataptr, Branch *const, const long);
long bytes_per_row(const long);
long childadd(Branch *const, const long, const long);
long cistrcmp(const char *const, const char *const);
Lvb_bool cleanup(void);
void clnclose(FILE *const, const char *const);
FILE *clnopen(const char *const, const char *const);
void clnremove(const char *const);
void crash(const char *const, ...);
void defaults_params(Params *const);
unsigned long checkpoint_uni(FILE *);
unsigned long restore_uni(FILE *);
void checkpoint_treestack(FILE *, Treestack *, Dataptr, Lvb_bool b_with_sset);
void restore_treestack(FILE *, Treestack *, Dataptr, Lvb_bool b_with_sset);
void dna_makebin(Dataptr restrict, DataSeqPtr matrix_seq, Lvb_bit_lentgh **);
#ifdef MAP_Reduce
	long deterministic_hillclimb(Dataptr, Treestack *, const Branch *const, Params rcstruct,
			long, FILE * const, long *, Lvb_bool, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer);
#else
	long deterministic_hillclimb(Dataptr, Treestack *, const Branch *const, Params rcstruct,
			long, FILE * const, long *, int myMPIid, Lvb_bool);
#endif

void dnapars_wrapper(void);
char *f2str(FILE *const);
Lvb_bool file_exists(const char *const);
void get_bootstrap_weights(long *, long, long);

double get_initial_t(Dataptr, const Branch *const, Params rcstruct, long, int myMPIid, Lvb_bool);
int getparam(Params *, int argc, char **argv);
long getplen(Dataptr restrict, Branch *, Params rcstruct, const long, long *restrict p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs);

void alloc_memory_to_getplen(Dataptr matrix, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
double get_predicted_length(double, double, long, long, long, long);
double get_predicted_trees(double, double, long, long, long, long);
long get_random_maxaccept(void);
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
void rowfree(DataSeqPtr, int n_lines);
int phylip_dna_matrin(char *, int, Dataptr, DataSeqPtr);
void phylip_mat_dims_in(char *, int, long *, long *, int *);

void randtree(Dataptr, Branch *const);
long randpint(const long);

void scream(const char *const, ...);
void ss_init(Dataptr, Branch *, Lvb_bit_lentgh **);
char *supper(char *const s);
Branch *treealloc(Dataptr restrict, Lvb_bool b_with_sset);
long tree_bytes(Dataptr restrict matrix);
long tree_bytes_whitout_sset(Dataptr restrict matrix);
void treeclear(Dataptr, Branch *const);
void treecopy(Dataptr restrict, Branch *const, const Branch *const, Lvb_bool b_with_sset);
long treecmp(Dataptr matrix, const Branch *const tree_1, const Branch *const tree_2, long root, Lvb_bool b_First);
void treedump(Dataptr, FILE *const, const Branch *const, Lvb_bool b_with_sset);
void treedump_b(Dataptr, FILE *const, const Branch *const, Lvb_bool);
void treedump_screen(Dataptr matrix, const Branch *const tree);
void treestack_clear(Treestack *);
long treestack_cnt(Treestack);
long treestack_dump(Dataptr, Treestack *, FILE *const);
void treestack_free(Treestack *);
long treestack_transfer(Dataptr, Treestack *, Treestack *, Lvb_bool b_with_sset);
long treestack_pop(Dataptr, Branch *, long *, Treestack *, Lvb_bool b_with_sset);
Treestack *treestack_new(void);
long treestack_print(Dataptr, DataSeqPtr restrict matrix_seq_data, Treestack *, FILE *const, Lvb_bool);

#ifdef MAP_Reduce
	uint64_t tree_setpush(Dataptr matrix, const Branch *const tree, const long root, MapReduce *mrObj, MISC *misc);
	void map_clean(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
	void reduce_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
	void reduce_sets(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
	void reduce_filter(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
	void print_sets(Dataptr matrix, Treestack *sp, MISC *misc);
#endif
IterationTemperature *get_alloc_main_calc_iterations(void);
void add_temperature_cal_iterations(IterationTemperature *p_data, SendInfoToMaster *p_info_temp, int n_process);
Lvb_bool is_possible_to_continue(IterationTemperature *p_data, double d_temperature, int n_iteration, int l_tree_length, int n_max_number_process, int n_count_call);
void release_main_calc_iterations(IterationTemperature *p_data);

long treestack_push_only(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset);
long treestack_push(Dataptr, Treestack *, const Branch *const, long, Lvb_bool b_with_sset);
void treeswap(Branch **const, long *const, Branch **const, long *const);
void uint32_dump(FILE *, Lvb_bit_lentgh);
long words_per_row(const long);


#endif /* LVB_LVB_H */

#endif
