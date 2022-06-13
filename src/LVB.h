/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Chang Sik Kim,
Maximilian Strobl and Martyn Winn
All rights reserved.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2020 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.

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

/* ========== LVB.h - main header for lvb ========== */

#ifndef LVB_LVB_H
#define LVB_LVB_H

#include <vector>

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
#include <sys/utsname.h>
#include <time.h>

#include "DataStructure.h"
#include "LogFile.h"
#include "MyMaths.h"
#include "RandomNumberGenerator.h"
#include "sys/stat.h"

#ifdef LVB_MAPREDUCE
	#include <iostream>
	#include <mpi.h>
	#include <omp.h>
	#include "MapReduce.h"

	using namespace MAPREDUCE_NS;
	#define __STDC_LIMIT_MACROS
#endif

/* DNA bases: bits to set in statesets */
#define A_BIT 0b0001		/* (1U << 0) */
#define C_BIT 0b0010		/* (1U << 1) */
#define G_BIT 0b0100		/* (1U << 2) */
#define T_BIT 0b1000		/* (1U << 3) */

#define NIBBLE_MASK 		017			/* space for one stateset, all bits set to 1 */
#define NIBBLE_WIDTH 		4			/* width of nibble in bits */
#define NIBBLE_WIDTH_BITS	2			/* bitwise multiply the NIBBLE_WIDTH */

typedef uint64_t Lvb_bit_length;								/* define 64 bits */
#define NUMBER_OF_BITS										64
#define LENGTH_WORD											16	/* length of number packed bases */
#define LENGTH_WORD_BITS_MULTIPLY							4	/* multiply of number packed bases */
#define MINIMUM_WORDS_PER_SLICE_GETPLEN						30  /* minimum words per slice that run gplen threading */
#define MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING		60 /* need to have this size to activate the threading */
#define MASK_SEVEN											0x7777777777777777U
#define MASK_EIGHT											0x8888888888888888U

/* limits that could be changed but, if increased enormously, might lead to
 * some trouble at some point */
#define MAX_N 1000000		/* max. rows */
#define MAX_M 5000000		/* max. cols */

/* values some people may feel the dangerous urge to change */
#define LVB_INPUTSTRING_SIZE 2000	/* max. bytes for interactive input */
#define STAT_LOG_INTERVAL 50000		/* min. interval for progress log */
#define REROOT_INTERVAL 1000		/* change root every ... updates */


/* implementation-independent limits */
#define LVB_EPS 1E-11					/* 0.0 < DBL_EPSILON < LVB_EPS */
#define MIN_M 1L						/* min. no. of characters for any analysis */
#define UNSET (-1)						/* value of integral vars when unset */
#define MAX_BRANCHES (2 * MAX_N - 3)	/* max. branches per tree */
#define MIN_BRANCHES (2 * MIN_N - 3)	/* max. branches per tree */
#define MIN_N 5L						/* min. no. of objs, for rearrangeable tree */
#define MAX_ALLOC ((size_t) (INT_MAX - 2))	/* max. bytes per dyn. alloc. */
#define MAXSTATES 5						/* max. "true" states in data MSA */

/* limits that could be changed but are likely to be OK */
#define FROZEN_T 0.0001		/* consider system frozen if temp < FROZEN_T */

typedef	struct	/* object set derived from a cladogram */
{
	long *set;	/* arrays of object sets */
	long cnt;	/* sizes of object sets */
}	Objset;

/* branch of tree */
typedef struct
{
    long parent;		/* parent branch number, UNSET in root */
    long left;			/* index of first child in tree array */
    long right;			/* index of second child in tree array */
    long changes;		/* changes associated with this branch */
    Lvb_bit_length *sitestate;	/* statesets for all sites */
} TREESTACK_TREE_NODES; // node

/* tree stacks */
typedef struct
{
	TREESTACK_TREE_NODES *tree;	/* pointer to first branch in tree array */
	long root;		/* root of tree */
    Objset *p_sitestate;	/* array with sitestate with the root always on zero */
} TREESTACK_TREES;

typedef struct
{
	long size;			/* number of trees currently allocated for */
	long next;			/* next unused element of stack */
    TREESTACK_TREES *stack;	/* pointer to first element in stack */
} TREESTACK;

/* simulated annealing parameters */
#define MAXACCEPT_MIN 5L		/* minimum value for maxaccept */
#define MAXACCEPT_MAX 5L		/* maximum value for maxaccept */
#define MAXACCEPT_SLOW 5L	/* maxaccept for "slow" searches */
#define MAXPROPOSE_SLOW 2000L	/* maxpropose for "slow" searches */
#define MAXFAIL_SLOW 40L	/* maxfail for "slow" searches */

/* fixed file names */
#define MATFNAM "infile"	/* MSA file name */
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

#ifdef LVB_MAPREDUCE
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

#endif

#ifdef __cplusplus
	#define restrict    /* nothing */
#endif

void *alloc(const size_t, const char *const);
long arbreroot(Dataptr, TREESTACK_TREE_NODES *const, const long);
long bytes_per_row(const long);
long childadd(TREESTACK_TREE_NODES *const, const long, const long);
long cistrcmp(const char *const, const char *const);
Lvb_bool cleanup(void);
void clnclose(FILE *const, const char *const);
FILE *clnopen(const char *const, const char *const);
void clnremove(const char *const);
void crash(const char *const, ...);
void dnapars_wrapper(void);
char *f2str(FILE *const);
Lvb_bool file_exists(const char *const);
void getparam(Parameters *, int argc, char **argv);
long getplen(Dataptr restrict, TREESTACK_TREE_NODES *, Parameters rcstruct, const long, long *restrict p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs);
void alloc_memory_to_getplen(Dataptr MSA, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
double get_predicted_length(double, double, long, long, long, long);
double get_predicted_trees(double, double, long, long, long, long);
long getroot(const TREESTACK_TREE_NODES *const);
void lvb_assertion_fail(const char *, const char *, int);
void lvb_initialize(void);
Dataptr lvb_matrin(const char *);
long lvb_reroot(Dataptr restrict, TREESTACK_TREE_NODES *const BranchArray, const long oldroot, const long newroot, Lvb_bool b_with_sitestate);
void lvb_treeprint (Dataptr, FILE *const, const TREESTACK_TREE_NODES *const, const long);

void matchange(Dataptr, const Parameters);
Dataptr matrin(const char *const);
void mutate_deterministic(Dataptr restrict, TREESTACK_TREE_NODES *const, const TREESTACK_TREE_NODES *const, long, long, Lvb_bool);
void mutate_spr(Dataptr restrict, TREESTACK_TREE_NODES *const, const TREESTACK_TREE_NODES *const, long);
void mutate_nni(Dataptr restrict, TREESTACK_TREE_NODES *const, const TREESTACK_TREE_NODES *const, long);
void mutate_tbr(Dataptr restrict, TREESTACK_TREE_NODES *const, const TREESTACK_TREE_NODES *const, long);
char *nextnonwspc(const char *);
void nodeclear(TREESTACK_TREE_NODES *const, const long);
long objreroot(TREESTACK_TREE_NODES *const, const long, const long);
void params_change(Parameters *);
void phylip_dna_matrin(char *, int, Dataptr);
void phylip_mat_dims_in(char *, int, long *, long *, int *);
void PullRandomTree(Dataptr, TREESTACK_TREE_NODES *const);
long randpint(const long);
void rowfree(Dataptr);
void scream(const char *const, ...);
void ss_init(Dataptr, TREESTACK_TREE_NODES *, Lvb_bit_length **);
char *supper(char *const s);
TREESTACK_TREE_NODES *treealloc(Dataptr restrict, Lvb_bool b_with_sitestate);
long tree_bytes(Dataptr restrict MSA);
long tree_bytes_without_sitestate(Dataptr restrict MSA);
void treeclear(Dataptr, TREESTACK_TREE_NODES *const);
void treecopy(Dataptr restrict, TREESTACK_TREE_NODES *const, const TREESTACK_TREE_NODES *const, Lvb_bool b_with_sitestate);
void treedump(Dataptr, FILE *const, const TREESTACK_TREE_NODES *const, Lvb_bool b_with_sitestate);
void treedump_screen(Dataptr MSA, const TREESTACK_TREE_NODES *const tree);
void ClearTreestack(TREESTACK *);
long CountTreestack(TREESTACK);
void FreeTreestackMemory(Dataptr restrict MSA, TREESTACK *);
TREESTACK CreateNewTreestack(void);
long PullTreefromTreestack(Dataptr, TREESTACK_TREE_NODES *, long *, TREESTACK *, Lvb_bool b_with_sitestate);
long CompareTreeToTreestack(Dataptr, TREESTACK *, const TREESTACK_TREE_NODES *const, const long, Lvb_bool b_with_sitestate);
int PrintTreestack(Dataptr, TREESTACK *, FILE *const, Lvb_bool onerandom);
void SwapTrees(TREESTACK_TREE_NODES **const, long *const, TREESTACK_TREE_NODES **const, long *const);
long words_per_row(const long);
int count(TREESTACK_TREE_NODES *const, int);
int addtoarray(TREESTACK_TREE_NODES *const, int, int *, int);
void copy_sitestate(Dataptr restrict MSA, Objset *p_sitestate_1);
void DNAToBinary(Dataptr restrict, Lvb_bit_length **);
void makesets(Dataptr restrict, const TREESTACK_TREE_NODES *const tree_2, const long root);
long setstcmp_with_sitestate2(Dataptr MSA, Objset *const oset_1);
long TopologyComparison(Dataptr restrict, Objset *, const TREESTACK_TREE_NODES *const, Lvb_bool b_first);
double StartingTemperature(Dataptr, const TREESTACK_TREE_NODES *const, Parameters rcstruct, long, Lvb_bool);
long PushCurrentTreeToStack(Dataptr, TREESTACK *, const TREESTACK_TREE_NODES *const, const long, Lvb_bool b_with_sitestate);

#ifdef LVB_MAPREDUCE
long Anneal(Dataptr restrict, TREESTACK *, TREESTACK *, const TREESTACK_TREE_NODES *const, Parameters rcstruct, long, const double,
	const long, const long, const long, FILE *const, long *, Lvb_bool, MISC *misc, MapReduce *mrStackTree, MapReduce *mrBuffer);

void defaults_params(Parameters *const);
long deterministic_hillclimb(Dataptr, TREESTACK *, const TREESTACK_TREE_NODES *const, Parameters rcstruct,
	long, FILE * const, long *, Lvb_bool, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer);
uint64_t tree_setpush(Dataptr MSA, const TREESTACK_TREE_NODES *const tree, const long root, MapReduce *mrObj, MISC *misc);
void map_clean(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void reduce_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_sets(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void reduce_filter(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void print_sets(Dataptr MSA, TREESTACK *sp, MISC *misc);

#else
long Anneal(Dataptr restrict, TREESTACK *, TREESTACK *, const TREESTACK_TREE_NODES *const, Parameters rcstruct, long, const double,
 const long, const long, const long, FILE *const, long *, Lvb_bool);

void defaults_params(Parameters *const prms);
long deterministic_hillclimb(Dataptr, TREESTACK *, const TREESTACK_TREE_NODES *const, Parameters rcstruct,
	long, FILE * const, long *, Lvb_bool);
void dump_stack_to_screen(Dataptr MSA, TREESTACK *sp);

#endif

#endif /* LVB_LVB_H */
