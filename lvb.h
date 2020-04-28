/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maxi
milian Strobl and Chris Wood.
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

/* ********** lvb.h - main header for lvb ********** */

#ifndef LVB_LVB_H
#define LVB_LVB_H

#include "DataStructure.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include "sys/stat.h"
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "myuni.h"
#include "mymaths.h"
#include <sys/utsname.h>
#include "clock.h"
#include "log.h"

#if defined (LVB_MAPREDUCE) || defined (LVB_PARALLEL_SEARCH)
#include <mpi.h>
#endif

#ifdef LVB_MAPREDUCE
#include "mapreduce.h"
#include <omp.h>
#include <iostream>

using namespace MAPREDUCE_NS;
using namespace std;
#define __STDC_LIMIT_MACROS

#endif


/* the program */
#define PROGNAM "lvb"			/* program file name */
#define LVB_VERSION "3.5"		/* version of program */
#define LVB_SUBVERSION "(February 2019)"	/* version details e.g. date */

#ifdef LVB_PARALLEL_SEARCH
#define	MPI_SEND_ONLY_MATRIX_NAMES	/* if defined only send the names of the matrix */
									/* sometimes the data matrix are huge and it's only necessary to pass */
    								/* the names of the other process */
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
#define STAT_LOG_INTERVAL 100000		/* min. interval for progress log */
#ifdef LVB_PARALLEL_SEARCH
#define UNSET (-1)			/* value of integral vars when unset */
#define CHECKPOINT_INTERVAL 1800	/* checkpoint ~every ... seconds */
#define CHECKPOINT_FNAM_BASE "lvb_checkpoint"
#endif
#define REROOT_INTERVAL 1000		/* change root every ... updates */


/* implementation-independent limits */
#define LVB_EPS 1E-11					/* 0.0 < DBL_EPSILON < LVB_EPS */
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
	long *set;	/* arrays of object sets */
	long cnt;	/* sizes of object sets */
}	Objset;

#ifdef LVB_PARALLEL_SEARCH
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
#endif

/* branch of tree */
typedef struct
{
    long parent;		/* parent branch number, UNSET in root */
    long left;			/* index of first child in tree array */
    long right;			/* index of second child in tree array */
    long changes;		/* changes associated with this branch */
    Lvb_bit_length *sset;	/* statesets for all sites */
} Branch;

/* tree stacks */
typedef struct
{
	Branch *tree;	/* pointer to first branch in tree array */
	long root;		/* root of tree */
    Objset *p_sset;	/* array with sset with the root always on zero */
} Treestack_element;

typedef struct
{
	long size;			/* number of trees currently allocated for */
	long next;			/* next unused element of stack */
    Treestack_element *stack;	/* pointer to first element in stack */
} Treestack;

/* simulated annealing parameters */
#define MAXACCEPT_MIN 5L		/* minimum value for maxaccept */
#ifdef LVB_PARALLEL_SEARCH
#define MAXACCEPT_MAX 500L		/* maximum value for maxaccept */
#else
#define MAXACCEPT_MAX 5L		/* maximum value for maxaccept */
#endif
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
 * calls CrashVerbosely() not abort(), and works whether or not NDEBUG is defined */
#define lvb_assert(test) ((void) ((test) || (LVBAssertFail(#test, __FILE__, __LINE__), 0)))

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

void *Alloc(const size_t, const char *const);
long RandomTreeRoot(Dataptr, Branch *const, const long);
long BytesPerRow(const long);
long CompareStrings(const char *const, const char *const);
Lvb_bool CleanExit(void);
void CheckFileClosure(FILE *const, const char *const);
FILE *CheckFileOpening(const char *const, const char *const);
void CrashVerbosely(const char *const, ...);
double AnnealStartingTemperature(Dataptr, const Branch *const, Params rcstruct, long, Lvb_bool);
long MinimumTreeLength(const Dataptr);
void PassSearchParameters(Params *, int argc, char **argv);
long CurrentTreeLength(Dataptr restrict, Branch *, Params rcstruct, const long, long *restrict p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs);
void AllocCurrentTreeLength(Dataptr matrix, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
void FreeCurrentTreeLength(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs);
void LVBAssertFail(const char *, const char *, int);
long RerootTree(Dataptr restrict, Branch *const barray, const long oldroot, const long newroot, Lvb_bool b_with_sset);
void PrintCurrentTree (Dataptr, FILE *const, const Branch *const, const long);
void CutMatrixColumns(Dataptr, const Params);
void SingleTopologyChange(Dataptr restrict, Branch *const, const Branch *const, long, long, Lvb_bool);
void HeuristicSPR(Dataptr restrict, Branch *const, const Branch *const, long);
void HeuristicNNI(Dataptr restrict, Branch *const, const Branch *const, long);
void HeuristicTBR(Dataptr restrict, Branch *const, const Branch *const, long);
char *TestNextNonWhiteSpaceCharacter(const char *);
void NodeClear(Branch *const, const long);
void CheckDNAMatrixInput(char *, int, Dataptr);
void ReadDNAMatrix(char *, int, long *, long *, int *);
void GenerateRandomTree(Dataptr, Branch *const);
long RandomNumberGenerator(const long);
void FreeRowStrings(Dataptr);
void PrintError(const char *const, ...);
void CopyCurrentStates(Dataptr, Branch *, Lvb_bit_length **);
char *ConvertToUpperCase(char *const s);
Branch *AllocBlankTreeArray(Dataptr restrict, Lvb_bool b_with_sset);
long TreeBytes(Dataptr restrict matrix);
long TreeBytesWithoutStates(Dataptr restrict matrix);
void ClearTreeArray(Dataptr, Branch *const);
void CopyTreeArray(Dataptr restrict, Branch *const, const Branch *const, Lvb_bool b_with_sset);
void PushTreeToStream(Dataptr, FILE *const, const Branch *const, Lvb_bool b_with_sset);
void PushTreeToScreen(Dataptr matrix, const Branch *const tree);
void ClearTreestack(Treestack *);
void FreeTreestack(Dataptr restrict matrix, Treestack *);
Treestack NewTreestack(void);
long PullTreeOffTreestack(Dataptr, Branch *, long *, Treestack *, Lvb_bool b_with_sset);
long PushTreeOntoTreestack(Dataptr, Treestack *, const Branch *const, const long, Lvb_bool b_with_sset);
int PrintTreestack(Dataptr, Treestack *, FILE *const, Lvb_bool onerandom);
void SwapTrees(Branch **const, long *const, Branch **const, long *const);
long WordsPerRow(const long);
int CountSubtreeNodes(Branch*const, int);
int AddSubtreeNodesToArray(Branch *const, int, int *, int);
void PushObjsectStatestoScreen(Dataptr matrix, Objset *oset_1);
void CopySingleStates(Dataptr restrict matrix, Objset *p_sset_1);
void ConvertDNAToStates(Dataptr restrict, Lvb_bit_length **);
void MakeArrayOfStateSets(Dataptr restrict, const Branch *const tree_2, const long root);
long CompareObjsectStatesWithSecondStates(Dataptr matrix, Objset *const oset_1);
long CompareTrees(Dataptr restrict, Objset *, const Branch *const, Lvb_bool b_first);
void PrintCopyright();
void PrintBanner();
void GetSystemTime();
bool LogFileExists(const char *filename);

#ifdef LVB_MAPREDUCE
long Anneal(Dataptr restrict, Treestack *, Treestack *, const Branch *const, Params rcstruct, long, const double,
	const long, const long, const long, FILE *const, long *, Lvb_bool, MISC *misc, MapReduce *mrStackTree, MapReduce *mrBuffer);

void GetDefaultParameters(Params *const);
long HillClimbingOptimization(Dataptr, Treestack *, const Branch *const, Params rcstruct,
	long, FILE * const, long *, Lvb_bool, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer);
uint64_t TreeSetPush(Dataptr matrix, const Branch *const tree, const long root, MapReduce *mrObj, MISC *misc);
void MapClean(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);
void ReduceCount(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void ReduceFilter(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);
void PrintSets(Dataptr matrix, Treestack *sp, MISC *misc);
long PushTreeOntoTreestackOnly(Dataptr, Treestack *, const Branch *const, const long, Lvb_bool b_with_sset);

#elif LVB_PARALLEL_SEARCH
long Anneal(Dataptr restrict, Treestack *, const Branch *const, Params *p_rcstruct, long, const double,
	const long, const long, const long, FILE *const, long *, int, int *p_n_state_progress, int *p_n_number_tried_seed, Lvb_bool);
long HillClimbingOptimization(Dataptr, Treestack *, const Branch *const, Params rcstruct,
	long, FILE * const, long *, int myMPIid, Lvb_bool);
unsigned long checkpoint_uni(FILE *);
unsigned long restore_uni(FILE *);
void checkpoint_treestack(FILE *, Treestack *, Dataptr, Lvb_bool b_with_sset);
void restore_treestack(FILE *, Treestack *, Dataptr, Lvb_bool b_with_sset);
IterationTemperature *get_alloc_main_calc_iterations(void);
void add_temperature_cal_iterations(IterationTemperature *p_data, SendInfoToMaster *p_info_temp, int n_process);
Lvb_bool is_possible_to_continue(IterationTemperature *p_data, double d_temperature, int n_iteration, int l_tree_length, int n_max_number_process, int n_count_call);
void release_main_calc_iterations(IterationTemperature *p_data);

#else
long Anneal(Dataptr restrict, Treestack *, Treestack *, const Branch *const, Params rcstruct, long, const double,
 const long, const long, const long, FILE *const, long *, Lvb_bool);

void GetDefaultParameters(Params *const prms);
long HillClimbingOptimization(Dataptr, Treestack *, const Branch *const, Params rcstruct,
	long, FILE * const, long *, Lvb_bool);

#endif

#endif /* LVB_LVB_H */