/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** lvb.h - main header for lvb ********** */

#ifndef LVB_LVB_H
#define LVB_LVB_H

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
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

/* the program */
#define PROGNAM "lvb"			/* program file name */
#define LVB_VERSION "3.1"		/* version of program */
#define LVB_SUBVERSION "(2 June 2014)"	/* version details e.g. date */

/* DNA bases: bits to set in statesets */
#define A_BIT (1U << 0)
#define C_BIT (1U << 1)
#define G_BIT (1U << 2)
#define T_BIT (1U << 3)
#define O_BIT (1U << 4)

/* values some people may feel the dangerous urge to change */
#define LVB_FNAMSIZE 2000		/* maximum bytes for file names */
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
#define MAX_BOOTSTRAPS 1000000	/* max. bootstrap replicates */
#define FROZEN_T 0.0001		/* consider system frozen if temp < FROZEN_T */

/* unchangeable types */
typedef enum { LVB_FALSE, LVB_TRUE } Lvb_bool;	/* boolean type */

/* verboseness level (0 = nonverbose, 1 = verbose */
#define VERBOSE_OUTPUT LVB_FALSE

/* branch of tree */
typedef struct
{
    long parent;		/* parent branch number, UNSET in root */
    long left;			/* index of first child in tree array */
    long right;			/* index of second child in tree array */
    long changes;		/* changes associated with this branch */
    unsigned char *sset;	/* statesets for all sites */

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

/* user- or programmer-configurable parameters */
typedef struct
{
    int seed;			/* seed for random number generator */
    long verbose;		/* verboseness level */
    long bootstraps;		/* number of bootstrap replicates */
    Lvb_bool fifthstate;	/* if LVB_TRUE, '-' is 'O'; otherwise is '?' */
    int cooling_schedule;   /* cooling schedule: 0 is geometric, 1 is linear */
    char *p_file_name;		/* input file name */
    int n_file_format;		/* number of file format, must be FORMAT_PHYLIP, FORMAT_FASTA, FORMAT_NEXUS, FORMAT_MSF, FORMAT_CLUSTAL*/
    char *p_file_name_out;	/* output file name */
} Params;

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


#define FORMAT_PHYLIP 		0
#define FORMAT_FASTA 		1
#define FORMAT_NEXUS 		2
#define FORMAT_MSF 			3
#define FORMAT_CLUSTAL 		4

/* assert-like macro, differing in that it writes to standard output,
 * calls crash() not abort(), and works whether or not NDEBUG is defined */
#define lvb_assert(test) ((void) ((test) || (lvb_assertion_fail(#test, __FILE__, __LINE__), 0)))

/* PHYLIP global data */
extern long chars;	/* defined in dnapars.c */

/* LVB global functions */
void *alloc(const size_t, const char *const);
long anneal(Dataptr, Treestack *, const Branch *const, long, const double,
 const long, const long, const long, FILE *const, const long *, long *, const int, Lvb_bool);
long arbreroot(Dataptr, Branch *const, const long);
long childadd(Branch *const, const long, const long);
long cistrcmp(const char *const, const char *const);
Lvb_bool cleanup(void);
void clnclose(FILE *const, const char *const);
FILE *clnopen(const char *const, const char *const);
void clnremove(const char *const);
void crash(const char *const, ...);
void defaults_params(Params *const prms);
long deterministic_hillclimb(Dataptr, Treestack *, const Branch *const, long,
    FILE * const, const long *, long *, Lvb_bool);
void dna_makebin(const Dataptr, Lvb_bool, unsigned char **);
void dnapars_wrapper(void);
char *f2str(FILE *const);
Lvb_bool file_exists(const char *const);
void get_bootstrap_weights(long *, long, long);
double get_initial_t(Dataptr, const Branch *const, long, const long *, Lvb_bool);
long getminlen(const Dataptr);
void getparam(Params *);
long getplen(Dataptr, Branch *, const long, const long *);
double get_predicted_length(double, double, long, long, long, long);
double get_predicted_trees(double, double, long, long, long, long);
long getroot(const Branch *const);
void lvb_assertion_fail(const char *, const char *, int);
void lvb_initialize(void);
Dataptr lvb_matrin(const char *);
long lvb_reroot(Dataptr, Branch *const barray, const long oldroot, const long newroot);
void lvb_treeprint (Dataptr, FILE *const, const Branch *const, const long);
void matchange(Dataptr, const Params);
Dataptr matrin(const char *const);
void mutate_deterministic(Dataptr, Branch *const, const Branch *const, long, long, Lvb_bool);
void mutate_spr(Dataptr, Branch *const, const Branch *const, long);
void mutate_nni(Dataptr, Branch *const, const Branch *const, long);
char *nextnonwspc(const char *);
void nodeclear(Branch *const, const long);
long objreroot(Branch *const, const long, const long);
void params_change(Params *);
void phylip_dna_matrin(char *, int, Dataptr);
void phylip_mat_dims_in(char *, int, long *, long *, int *);
void randtree(Dataptr, Branch *const);
long randpint(const long);
void rowfree(Dataptr);
char *salloc(const long, const char *const);
void scream(const char *const, ...);
void ss_init(Dataptr, Branch *, unsigned char **);
char *supper(char *const s);
Branch *treealloc(Dataptr);
void treeclear(Dataptr, Branch *const);
void treecopy(Dataptr, Branch *const, const Branch *const);
long treecmp(Dataptr, const Branch *const, const long, const Branch *const, long);
void treedump(Dataptr, FILE *const, const Branch *const);
void treestack_clear(Treestack *);
long treestack_cnt(Treestack);
long treestack_dump(Dataptr, Treestack *, FILE *const);
void treestack_free(Treestack *);
Treestack treestack_new(void);
long treestack_transfer(Dataptr, Treestack *, Treestack *);
long treestack_pop(Dataptr, Branch *, long *, Treestack *);
long treestack_print(Dataptr, Treestack *, FILE *const, Lvb_bool);
long treestack_push(Dataptr, Treestack *, const Branch *const, const long);
void treeswap(Branch **const, long *const, Branch **const, long *const);


#endif /* LVB_LVB_H */
