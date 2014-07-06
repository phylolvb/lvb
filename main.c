/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** main.c - LVB ********** */

#include "lvb.h"

static Treestack bstack_overall;	/* overall best tree stack */

static void check_stdout(void)
/* Flush standard output, and crash verbosely on error. */
{
    if (fflush(stdout) == EOF)
        crash("write error on standard output");	/* may not work! */
    if (ferror(stdout))
	crash("file error on standard output");		/* may not work! */
}	/* end check_stdout() */

static void smessg(long start, long cycle)
/* print cycle start message */
{
    printf("Beginning start %ld cycle %ld\n", start, cycle);
    check_stdout();

} /* end smessg() */

static void writeinf(Params prms)
/* write initial details to standard output */
{
    printf("\n");

    /* configurable details */
    printf("format               = ");
    if (prms.interleaved == LVB_TRUE)
	printf("INTERLEAVED\n");
    else if (prms.interleaved == LVB_FALSE)
	printf("SEQUENTIAL\n");
    else
	lvb_assert(0);
    printf("treatment of '-'     = ");
    if (prms.fifthstate == LVB_TRUE)
	printf("FIFTH STATE\n");
    else
	printf("UNKNOWN\n");
    printf("cooling schedule     = ");
    if(prms.cooling_schedule == 0)
	printf("GEOMETRIC\n");
    else
	printf("LINEAR\n");
    printf("seed                 = %d\n", prms.seed);
    printf("bootstrap replicates = %ld\n", prms.bootstraps);

} /* end writeinf() */

static void logtree1(const Branch *const barray, const long start, const long cycle, long root)
/* log initial tree for cycle cycle of start start (in barray) to outfp */
{
    static char outfnam[LVB_FNAMSIZE]; 	/* current file name */
    int fnamlen;			/* length of current file name */
    FILE *outfp;			/* output file */

    fnamlen = sprintf(outfnam, "%s_start%ld_cycle%ld", TREE1FNAM, start, cycle);
    lvb_assert(fnamlen < LVB_FNAMSIZE);	/* shut door if horse bolted */

    /* create tree file */
    outfp = clnopen(outfnam, "w");
    lvb_treeprint(outfp, barray, root);
    clnclose(outfp, outfnam);

} /* end logtree1() */

static long getsoln(Params rcstruct, const long *weight_arr, long *iter_p,
 Lvb_bool log_progress)
/* get and output solution(s) according to parameters in rcstruct;
 * return length of shortest tree(s) found, using weights in weight_arr */
{
    int cooling_schedule = rcstruct.cooling_schedule; /* cooling schedule */
    extern Dataptr matrix;		/* data matrix */
    static char fnam[LVB_FNAMSIZE];	/* current file name */
    long fnamlen;			/* length of current file name */
    long i;				/* loop counter */
    double t0;		/* SA cooling cycle initial temp */
    long maxaccept = MAXACCEPT_SLOW;	/* SA cooling cycle maxaccept */
    long maxpropose = MAXPROPOSE_SLOW;	/* SA cooling cycle maxpropose */
    long maxfail = MAXFAIL_SLOW;	/* SA cooling cycly maxfail */
    long treec;				/* number of trees found */
    long treelength = LONG_MAX;		/* length of each tree found */
    long initroot;			/* initial tree's root */
    FILE *sumfp;			/* best length file */
    FILE *resfp;			/* results file */
    Branch *tree;			/* initial tree */
    Branch *user_tree_ptr = NULL;	/* user-specified initial tree */
    static unsigned char *enc_mat[MAX_N] = { NULL };	/* encoded data mat. */

	/* NOTE: These variables and their values are "dummies" and are no longer used in the
	 * current version of LVB. However, in order to keep the formatting of the output
	 * compatible with that of previous versions of LVB these variables will continue to be
	 * used and written to the summary files.
	 */
    long cyc = 0;				/* current cycle number */
    long start = 0;			/* current random (re)start number */
 
    /* dynamic "local" heap memory */
    tree = treealloc(brcnt(matrix->n), matrix->m);

    /* Allocation of the initial encoded matrix is non-contiguous because
     * this matrix isn't used much, so any performance penalty won't matter. */
    for (i = 0; i < matrix->n; i++)
        enc_mat[i] = alloc(matrix->m * sizeof(unsigned char), "state sets");
    
    dna_makebin(matrix, rcstruct.fifthstate, enc_mat);

    /* open and entitle statistics file shared by all cycles
     * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */

    if (rcstruct.verbose == LVB_TRUE) {
		sumfp = clnopen(SUMFNAM, "w");
		fprintf(sumfp, "StartNo\tCycleNo\tCycInit\tCycBest\tCycTrees\n");
    }
    else
    {
        sumfp = NULL;
    }
	
    /* determine starting temperature */
    randtree(tree);	/* initialise required variables */
    ss_init(tree, enc_mat, brcnt(matrix->n), matrix->m);
    initroot = 0;
    t0 = get_initial_t(tree, initroot, enc_mat, matrix->m, matrix->n,
	weight_arr, log_progress);
    

    randtree(tree);	/* begin from scratch */
    ss_init(tree, enc_mat, brcnt(matrix->n), matrix->m);
    initroot = 0;

    if (rcstruct.verbose)
	smessg(start, cyc);
    check_stdout();

    /* start cycles's entry in sum file
     * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions.  */
    if(rcstruct.verbose == LVB_TRUE)
    {
		fprintf(sumfp, "%ld\t%ld\t%ld\t", start, cyc, getplen(tree,
			initroot,matrix->m, matrix->n, weight_arr));
		logtree1(tree, start, cyc, initroot);
    }

    /* find solution(s) */
    treelength = anneal(&bstack_overall, tree, initroot, t0, maxaccept,
	maxpropose, maxfail, stdout, enc_mat, matrix->m, matrix->n, weight_arr,
	iter_p, cooling_schedule, log_progress);
    treestack_pop(tree, &initroot, &bstack_overall);
    treestack_push(&bstack_overall, tree, initroot);
    treelength = deterministic_hillclimb(&bstack_overall, tree,
	initroot, stdout, matrix->m, matrix->n, weight_arr, iter_p,
	log_progress);

	/* log this cycle's solution and its details 
	 * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */
    if (rcstruct.verbose == LVB_TRUE)
    {
	fnamlen = sprintf(fnam, "%s_start%ld_cycle%ld", RESFNAM, start, cyc);
	lvb_assert(fnamlen < LVB_FNAMSIZE);	/* really too late */
	resfp = clnopen(fnam, "w");
	treec = treestack_print(&bstack_overall, resfp, LVB_FALSE);
	clnclose(resfp, fnam);
	fprintf(sumfp, "%ld\t%ld\n", treelength, treec);

	/* won't use length summary file until end of next cycle */
	fflush(sumfp);
	if (ferror(sumfp))
	{
	    crash("write error on file %s", SUMFNAM);
	}
    }


    if (rcstruct.verbose == LVB_TRUE)
    	printf("Ending start %ld cycle %ld\n", start, cyc);
    check_stdout();

    if (rcstruct.verbose == LVB_TRUE)
    	clnclose(sumfp, SUMFNAM);

    /* "local" dynamic heap memory */
    free(tree);
    if (user_tree_ptr != NULL)
    	free(user_tree_ptr);

    return treelength;

} /* end getsoln() */

static void logstim(void)
/* log start time with message */
{
    time_t tim;	/* time */

    tim = time(NULL);
    printf("Starting at: %s", ctime(&tim));

} /* end logstim() */

int main(void)
{
    extern Dataptr matrix;	/* data matrix */
    int val;			/* return value */
    Params rcstruct;		/* configurable parameters */
    long i;			/* loop counter */
    long m;			/* sites per sequence */
    long n;			/* sequences in the data matrix */
    long iter;			/* iterations of annealing algorithm */
    long replicate_no = 0L;	/* current bootstrap replicate number */
    long trees_output_total = 0L;	/* number of trees output, overall */
    long trees_output;		/* number of trees output for current rep. */
    double total_iter = 0.0;	/* total iterations across all replicates */
    long final_length;		/* length of shortest tree(s) found */
    long m_including_constcols;	/* site count before constant sites removed */
    FILE *outtreefp;		/* best trees found overall */
    static long weight_arr[MAX_M];	/* weights for sites */
    Lvb_bool log_progress;	/* whether or not to log anneal search */

    /* global files */

    /* entitle standard output */
    printf("\n* This is %s version %s %s *\n"
	"* written by Daniel Barker and Maximilian Strobl *\n", PROGNAM,
	LVB_VERSION, LVB_SUBVERSION);
    printf("\n");
    printf("Download and support:\n");
    printf("http://eggg.st-andrews.ac.uk/lvb\n");
    printf("\n");
    printf("Literature reference:\n");
    printf("Barker, D. 2004. LVB: Parsimony and simulated annealing in the "
	"search for\n"
	"phylogenetic trees. Bioinformatics, 20, 274-275.\n");
    printf("\n");

    lvb_initialize();

    getparam(&rcstruct);
    phylip_mat_dims_in(&n, &m);
    logstim();

    matrix = phylip_dna_matrin(rcstruct.interleaved);
    lvb_assert((matrix->m == m) && (matrix->n == n));

    /* "file-local" dynamic heap memory: set up best tree stacks */
    bstack_overall = treestack_new();

    writeinf(rcstruct);
    m_including_constcols = matrix->m;
    matchange(matrix, rcstruct, rcstruct.verbose);	/* cut columns */

    if (rcstruct.verbose == LVB_TRUE)
    {
	printf("getminlen: %ld\n\n", getminlen(matrix));
    }
    rinit(rcstruct.seed);
    if (rcstruct.bootstraps > 0)
    {
	log_progress = LVB_FALSE;
	printf("\nReplicate:      Rearrangements: Trees output:   Length:\n");
    }
    else
	log_progress = LVB_TRUE;
    outtreefp = clnopen(OUTTREEFNAM, "w");
    do
    {
	iter = 0;
	if (rcstruct.bootstraps > 0)
	{
	    get_bootstrap_weights(weight_arr, matrix->m,
		m_including_constcols - matrix->m);
	}
	else
	{
	    for (i = 0; i < matrix->m; i++)
		weight_arr[i] = 1;
	}

	final_length = getsoln(rcstruct, weight_arr, &iter, log_progress);
	if (rcstruct.bootstraps > 0)
	    trees_output = treestack_print(&bstack_overall, outtreefp,
		LVB_TRUE);
	else
	    trees_output = treestack_print(&bstack_overall, outtreefp,
		LVB_FALSE);
	trees_output_total += trees_output;
	treestack_clear(&bstack_overall);
	replicate_no++;
	if (rcstruct.bootstraps > 0)
	{
	    printf("%-16ld%-16ld%-16ld%ld\n", replicate_no, iter, trees_output,
		final_length);
	    total_iter += (double) iter;
	} else
	    printf("\nRearrangements tried: %ld\n", iter);
    } while (replicate_no < rcstruct.bootstraps);
    clnclose(outtreefp, OUTTREEFNAM);

    printf("\n");

    if (rcstruct.bootstraps > 0)
	printf("Total rearrangements tried across all replicates: %g\n\n",
	 total_iter);

    if ((trees_output_total == 1L) && (rcstruct.bootstraps == 0))
    {
	printf("1 most parsimonious tree of length %ld written to file "
	    "'%s'\n", final_length, OUTTREEFNAM);
    }
    else
    {
	if (rcstruct.bootstraps > 0)
	{
	    lvb_assert(trees_output_total == rcstruct.bootstraps);
	    printf("%ld trees written to file '%s'\n", trees_output_total,
		OUTTREEFNAM);
	}
	else
	{
	    printf("%ld equally parsimonious trees of length %ld written to "
         "file '%s'\n", trees_output_total, final_length, OUTTREEFNAM);
	}
    }

    rowfree(matrix);

    /* "file-local" dynamic heap memory */
    treestack_free(&bstack_overall);

    if (cleanup() == LVB_TRUE)
	val = EXIT_FAILURE;
    else
	val = EXIT_SUCCESS;

    return val;

} /* end main() */
