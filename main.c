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

/* ========== main.c - LVB ========== */

#include "admin.h"
#include "lvb.h"

static Treestack bstack_overall;	/* overall best tree stack */
static Treestack stack_treevo;

static void check_stdout(void)
/* Flush standard output, and crash verbosely on error. */
{
    if (fflush(stdout) == EOF)
        CrashVerbosely("write error on standard output");	/* may not work! */
    if (ferror(stdout))
    	CrashVerbosely("file error on standard output");		/* may not work! */
}	/* end check_stdout() */

static void smessg(long start, long cycle)
/* print cycle start message */
{
    // printf("Beginning cycle \n\n");
    check_stdout();

} /* end smessg() */

static void writeinf(Params prms, Dataptr matrix, int argc, char **argv, int n_process)

/* write initial details to standard output */
{
	struct utsname buffer;
	errno = 0;
	if (uname(&buffer) !=0) 
	{
		perror("uname");
		exit(EXIT_FAILURE);
	}

	printf("Executing: ");
	printf("' ");
	for (int i = 0; i < argc; ++i)
	printf("%s ", argv[i]);
	printf("' at: ");
	GetSystemTime();
	printf("\n");

	printf("Analysis Properties: \n");
	printf("  Alignment:          '%s'\n", prms.file_name_in);
	printf("  MSA format:          ");
	if (prms.n_file_format == FORMAT_PHYLIP) printf("PHYLIP\n");
    else if (prms.n_file_format == FORMAT_FASTA) printf("FASTA\n");
    else if (prms.n_file_format == FORMAT_NEXUS) printf("NEXUS\n");
    else if (prms.n_file_format == FORMAT_CLUSTAL) printf("CLUSTAL\n");
    else{
    	fprintf (stderr, "Error, input format file not recognized\n");
    	abort();
    }

	printf("  MSA size:            %ld x %ld\n", matrix->n, matrix->original_m);
	printf("  Seed:                %d\n", prms.seed);
	printf("  Cooling schedule:    ");
    if(prms.cooling_schedule == 0) printf("GEOMETRIC\n");
    else printf("LINEAR\n");
	printf("  Algorithm: ");
    if(prms.algorithm_selection == 0) printf("          0 (SN)\n");
    else if(prms.algorithm_selection == 1) printf("          1 (SEQ-TNS)\n");
    else if(prms.algorithm_selection == 2) printf("          2 (PBS)\n");

	printf("\nParallelisation Properties: \n");
	printf("PThreads:              %d\n", prms.n_processors_available);
	printf("\n================================================================================\n");	
	printf("\nInitialising search: \n");

}


static void logtree1(Dataptr matrix, const Branch *const barray, const long start, const long cycle, long root)
/* log initial tree for cycle cycle of start start (in barray) to outfp */
{
    static char outfnam[LVB_FNAMSIZE]; 	/* current file name */
    int fnamlen;			/* length of current file name */
    FILE *outfp;			/* output file */

    fnamlen = sprintf(outfnam, "%s_start%ld_cycle%ld", TREE1FNAM, start, cycle);
    lvb_assert(fnamlen < LVB_FNAMSIZE);	/* shut door if horse bolted */

    /* create tree file */
    outfp = CheckFileOpening(outfnam, "w");
    PrintCurrentTree(matrix, outfp, barray, root);
    CheckFileClosure(outfp, outfnam);

} /* end logtree1() */

#ifdef LVB_MAPREDUCE
static long getsoln(Dataptr restrict matrix, Params rcstruct, long *iter_p, Lvb_bool log_progress,
				MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer)
#else
static long getsoln(Dataptr restrict matrix, Params rcstruct, long *iter_p, Lvb_bool log_progress)

#endif
/* get and output solution(s) according to parameters in rcstruct;
 * return length of shortest tree(s) found, using weights in weight_arr */
{
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
    Lvb_bit_length **enc_mat;	/* encoded data mat. */
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */
	#ifdef LVB_MAPREDUCE
	int *total_count;
	#endif

    /* NOTE: These variables and their values are "dummies" and are no longer
     * used in the current version of LVB. However, in order to keep the
     * formatting of the output compatible with that of previous versions of
     * LVB these variables will continue to be used and written to the summary
     * files.  */
    long cyc = 0;	/* current cycle number */
    long start = 0;	/* current random (re)start number */
 
    /* dynamic "local" heap memory */
    tree = AllocBlankTreeArray(matrix, LVB_TRUE);

    /* Allocation of the initial encoded matrix is non-contiguous because
     * this matrix isn't used much, so any performance penalty won't matter. */
    enc_mat = (Lvb_bit_length **) malloc((matrix->n) * sizeof(Lvb_bit_length *));
    for (i = 0; i < matrix->n; i++) 
		enc_mat[i] = (Lvb_bit_length *) Alloc(matrix->bytes, "state sets");
    ConvertDNAToStates(matrix, enc_mat);

    /* open and entitle statistics file shared by all cycles
     * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */

	#ifdef LVB_MAPREDUCE
	if (misc->rank == 0) {
	#endif

    if (rcstruct.verbose == LVB_TRUE) {
		sumfp = CheckFileOpening(SUMFNAM, "w");
		fprintf(sumfp, "StartNo\tCycleNo\tCycInit\tCycBest\tCycTrees\n");
    }
    else{
        sumfp = NULL;
    }
	#ifdef LVB_MAPREDUCE
	}
	#endif
	
    /* determine starting temperature */
    GenerateRandomTree(matrix, tree);	/* initialise required variables */
    CopyCurrentStates(matrix, tree, enc_mat);
    initroot = 0;
    t0 = AnnealStartingTemperature(matrix, tree, rcstruct, initroot, log_progress);

    GenerateRandomTree(matrix, tree);	/* begin from scratch */
    CopyCurrentStates(matrix, tree, enc_mat);
    initroot = 0;

    if (rcstruct.verbose) smessg(start, cyc);
    	check_stdout();

    /* start cycles's entry in sum file
     * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions.  */
	#ifdef LVB_MAPREDUCE
	if (misc->rank == 0) {
	#endif
    if(rcstruct.verbose == LVB_TRUE) {
        AllocCurrentTreeLength(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
		fprintf(sumfp, "%ld\t%ld\t%ld\t", start, cyc, CurrentTreeLength(matrix, tree, rcstruct, initroot, p_todo_arr, p_todo_arr_sum_changes, p_runs));
		FreeCurrentTreeLength(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
		logtree1(matrix, tree, start, cyc, initroot);
    }
	#ifdef LVB_MAPREDUCE
	}

		MPI_Barrier(MPI_COMM_WORLD);
		/* find solution(s) */
		maxaccept = RandomNumberGenerator(MAXACCEPT_MAX - MAXACCEPT_MIN) + MAXACCEPT_MIN;
		// printf("\nmaxaccept:%ld\n", maxaccept);
		treelength = Anneal(matrix, &bstack_overall, &stack_treevo, tree, rcstruct, initroot, t0, maxaccept,
				maxpropose, maxfail, stdout, iter_p, log_progress, misc, mrTreeStack, mrBuffer );

		long val = PullTreeOffTreestack(matrix, tree, &initroot, &bstack_overall, LVB_FALSE);
		PushTreeOntoTreestack(matrix, &bstack_overall, tree, initroot, LVB_FALSE);

		if(val ==  1) {
			misc->SB = 0;
			TreeSetPush(matrix, tree, initroot, mrBuffer, misc);
			mrTreeStack->add(mrBuffer);
			mrTreeStack->collate(NULL);
			mrTreeStack->reduce(ReduceFilter, NULL);

			mrBuffer->add(mrTreeStack);
			mrBuffer->collate(NULL);

			misc->count = (int *) Alloc( (bstack_overall.next+1) * sizeof(int), "integer array for tree compare using MapReduce");
			total_count = (int *) Alloc( (bstack_overall.next+1) * sizeof(int), "integer array for tree compare using MapReduce");
			for(int i=0; i<=bstack_overall.next; i++) misc->count[i] = 0;
			mrBuffer->reduce(ReduceCount, misc);

			for(int i=0; i<=bstack_overall.next; i++) total_count[i] = 0;
			MPI_Reduce( misc->count, total_count, bstack_overall.next+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

			int check_cmp = 1;
			if (misc->rank == 0) {
				for(int i=1; i<=bstack_overall.next; i++) {
					if (misc->nsets == total_count[i]) {
						check_cmp = 0; /* different */
						break;
					}
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&check_cmp, 1, MPI_INT, 0,    MPI_COMM_WORLD);
			if (check_cmp == 1) {
		//	  PushTreeOntoTreestackOnly(matrix, &bstack_overall, tree, initroot);
			  misc->ID = bstack_overall.next;
				  misc->SB = 1;
				  TreeSetPush(matrix, tree, initroot, mrBuffer, misc);
				  mrTreeStack->add(mrBuffer);
			}

			free(misc->count);
			free(total_count);
		//} else {
		//	PushTreeOntoTreestackOnly(matrix, &bstack_overall, tree, initroot);
		//	misc->ID = bstack_overall.next;
		//	misc->SB = 1;
		//        TreeSetPush(matrix, tree, initroot, mrBuffer, misc);
		//        mrTreeStack->add(mrBuffer);
		}

		treelength = HillClimbingOptimization(matrix, &bstack_overall, tree, rcstruct, initroot, stdout,
				iter_p, log_progress, misc, mrTreeStack, mrBuffer);
	#else
	    /* find solution(s) */
    treelength = Anneal(matrix, &bstack_overall, &stack_treevo, tree, rcstruct, initroot, t0, maxaccept, 
    maxpropose, maxfail, stdout, iter_p, log_progress);
    PullTreeOffTreestack(matrix, tree, &initroot, &bstack_overall, LVB_FALSE);
    PushTreeOntoTreestack(matrix, &bstack_overall, tree, initroot, LVB_FALSE);
	treelength = HillClimbingOptimization(matrix, &bstack_overall, tree, rcstruct, initroot, stdout,
		iter_p, log_progress);

	#endif

	/* log this cycle's solution and its details 
	 * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */

	#ifdef LVB_MAPREDUCE
	if (misc->rank == 0 ) {
	#endif
    if (rcstruct.verbose == LVB_TRUE){
		fnamlen = sprintf(fnam, "%s_start%ld_cycle%ld", RESFNAM, start, cyc);
		lvb_assert(fnamlen < LVB_FNAMSIZE);	/* really too late */
		resfp = CheckFileOpening(fnam, "w");
		treec = PrintTreestack(matrix, &bstack_overall, resfp, LVB_FALSE);
		CheckFileClosure(resfp, fnam);
		fprintf(sumfp, "%ld\t%ld\n", treelength, treec);

		/* won't use length summary file until end of next cycle */
		fflush(sumfp);
		if (ferror(sumfp)){
			CrashVerbosely("write error on file %s", SUMFNAM);
		}
    }


    if (rcstruct.verbose == LVB_TRUE) // printf("Ending start %ld cycle %ld\n", start, cyc);
    check_stdout();

    if (rcstruct.verbose == LVB_TRUE) CheckFileClosure(sumfp, SUMFNAM);
	#ifdef LVB_MAPREDUCE
	}
	#endif
    /* "local" dynamic heap memory */
    free(tree);
	for (i = 0; i < matrix->n; i++) free(enc_mat[i]);
    free(enc_mat);

    return treelength;

} /* end getsoln() */

/* set the number of processors to use */
void calc_distribution_processors(Dataptr matrix, Params rcstruct){
	int n_threads_temp = 0;
	if (matrix->nwords > MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING){
		do{
			n_threads_temp ++;
			matrix->n_slice_size_getplen = matrix->nwords / n_threads_temp;
		}while (matrix->n_slice_size_getplen > MINIMUM_WORDS_PER_SLICE_GETPLEN && n_threads_temp != rcstruct.n_processors_available);

		if (matrix->n_slice_size_getplen > MINIMUM_WORDS_PER_SLICE_GETPLEN){
			matrix->n_slice_size_getplen = matrix->nwords / n_threads_temp;
			matrix->n_threads_getplen = n_threads_temp;
		}
		else{
			matrix->n_threads_getplen = n_threads_temp - 1;
			matrix->n_slice_size_getplen = matrix->nwords / matrix->n_threads_getplen;
		}
	}
	else{
		matrix->n_threads_getplen = 1; /* need to pass for 1 thread because the number of words is to low */
	}
	// only to protect
	if (matrix->n_threads_getplen < 1) matrix->n_threads_getplen = 1;
}

int main(int argc, char **argv)
{
	Dataptr matrix;	/* data matrix */
	int val;			/* return value */
	Params rcstruct;		/* configurable parameters */
	long iter;			/* iterations of annealing algorithm */
	long trees_output_total = 0L;	/* number of trees output, overall */
	long trees_output;		/* number of trees output for current rep. */
	long final_length;		/* length of shortest tree(s) found */
	FILE *outtreefp;		/* best trees found overall */
	outtreefp = (FILE *) Alloc(sizeof(FILE), "alloc FILE");
	Lvb_bool log_progress;	/* whether or not to log anneal search */
	 

	#ifdef LVB_MAPREDUCE

		MPI_Init(&argc,&argv);

		MISC misc;

		MPI_Comm_rank(MPI_COMM_WORLD,&misc.rank);
		MPI_Comm_size(MPI_COMM_WORLD,&misc.nprocs);
		
		MapReduce *mrTreeStack = new MapReduce(MPI_COMM_WORLD);
		mrTreeStack->memsize = 1024;
		mrTreeStack->verbosity = 0;
		mrTreeStack->timer = 0;

		MapReduce *mrBuffer = new MapReduce(MPI_COMM_WORLD);
		mrBuffer->memsize = 1024;
		mrBuffer->verbosity = 0;
		mrBuffer->timer = 0;

		if(misc.rank == 0) {
	#else

		int nprocs = 1; /* Processor count*/

	#endif

    /* entitle standard output */
    PrintCopyright();
	PrintBanner();

	#ifdef LVB_MAPREDUCE
		}
	#endif

    /* start timer */ 
    clock_t Start = 0;
	clock_t End = 0;
    double Overall_Time_taken;

    Start = clock();
    LVBPreChecks();

    PassSearchParameters(&rcstruct, argc, argv);

    /* read and alloc space to the data structure */
	#ifdef LVB_MAPREDUCE
	matrix = (data *) Alloc(sizeof(DataStructure), "alloc data structure");
	matrix->row = NULL;
	matrix->rowtitle = NULL;
	#else
	matrix = Alloc(sizeof(DataStructure), "alloc data structure");
	#endif
    CheckDNAMatrixInput(rcstruct.file_name_in, rcstruct.n_file_format, matrix);

    /* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
	bstack_overall = NewTreestack();
	if(rcstruct.algorithm_selection ==2) 
    stack_treevo = NewTreestack();
        
    CutMatrixColumns(matrix, rcstruct);	/* cut columns */
	#ifdef LVB_MAPREDUCE
    writeinf(rcstruct, matrix, argc, argv, misc.nprocs);
	#else
	writeinf(rcstruct, matrix, argc, argv, nprocs);
	#endif
    calc_distribution_processors(matrix, rcstruct);

    if (rcstruct.verbose == LVB_TRUE) {
		#ifdef LVB_MAPREDUCE
    	if(misc.rank == 0) printf("Based on matrix provided, maximum parsimony tree length: %ld\n\n", MinimumTreeLength(matrix));
		#else
		printf("Based on matrix provided, maximum parsimony tree length: %ld\n\n", MinimumTreeLength(matrix));
		#endif
    }
    rinit(rcstruct.seed);
	log_progress = LVB_TRUE;

	#ifdef LVB_MAPREDUCE
	if (misc.rank == 0)
	#endif
    outtreefp = CheckFileOpening(rcstruct.file_name_out, "w");
    FILE * treEvo;
	treEvo = (FILE *) Alloc(sizeof(FILE), "alloc FILE");
    if(rcstruct.algorithm_selection ==2)
    treEvo = fopen ("treEvo.tre","w");
		iter = 0;     
		#ifdef LVB_MAPREDUCE
		final_length = getsoln(matrix, rcstruct, &iter, log_progress, &misc, mrTreeStack, mrBuffer);
	    if (misc.rank == 0) {
	       trees_output = PrintTreestack(matrix, &bstack_overall, outtreefp, LVB_FALSE);
	    }

		#else
		final_length = getsoln(matrix, rcstruct, &iter, log_progress);
		trees_output = PrintTreestack(matrix, &bstack_overall, outtreefp, LVB_FALSE);		
		#endif
		
		trees_output_total += trees_output;
        if(rcstruct.algorithm_selection ==2)
		PrintTreestack(matrix, &stack_treevo, treEvo, LVB_FALSE);
        ClearTreestack(&bstack_overall);
		printf("--------------------------------------------------------\n");
		#ifdef LVB_MAPREDUCE
		/* clean the TreeStack and buffer */
	    mrTreeStack->map( mrTreeStack, MapClean, NULL );
	    mrBuffer->map( mrBuffer, MapClean, NULL );
	    /* END clean the TreeStack and buffer */
		#endif

   if(rcstruct.algorithm_selection ==2)
    fclose(treEvo);
	
	#ifdef LVB_MAPREDUCE
	if (misc.rank == 0) {
	#endif
	CheckFileClosure(outtreefp, rcstruct.file_name_out);

    End = clock();
	#ifdef LVB_MAPREDUCE
	}
	#endif
	Overall_Time_taken = ((double) (End - Start)) /CLOCKS_PER_SEC;

	if (LogFileExists ("logfile.tsv"))
	{
		FILE * logfile;
    	logfile = fopen ("logfile.tsv","a+");
		fprintf (logfile, "%s\t%ld\t%ld\t%ld\t%.2lf\n", LVB_IMPLEMENTATION, iter, trees_output_total, final_length, Overall_Time_taken);
		fclose(logfile);
	}
	else {
		FILE * logfile;
	    logfile = fopen ("logfile.tsv","a+");
		fprintf (logfile, "Implementation\tRearrangements\tTopologies\tScore\tRuntime\n");
		fprintf (logfile, "%s\t%ld\t%ld\t%ld\t%.2lf\n", LVB_IMPLEMENTATION, iter, trees_output_total, final_length, Overall_Time_taken);
		fclose(logfile);
	}

	

	printf("\nSearch Complete\n");
	printf("\n================================================================================\n");
	printf("\nSearch Results:\n");
	printf("  Rearrangements evaluated: %ld\n", iter);
	printf("  Topologies recovered:     %ld\n", trees_output_total);
	printf("  Tree score:               %ld\n", final_length);
	printf("  Total runtime (seconds):  %.2lf\n", Overall_Time_taken);
	printf("\nAll topologies written to '%s'\n", rcstruct.file_name_out);


	/* "file-local" dynamic heap memory */
    if (rcstruct.algorithm_selection ==2)
    FreeTreestack(matrix, &stack_treevo);
	FreeTreestack(matrix, &bstack_overall);
    FreeRowStrings(matrix);
    free(matrix);

    if (CleanExit() == LVB_TRUE) val = EXIT_FAILURE;
    else val = EXIT_SUCCESS;

	#ifdef LVB_MAPREDUCE
	FreeTreestack(matrix, &bstack_overall);
	    MPI_Barrier(MPI_COMM_WORLD);

	    delete mrTreeStack;
	    delete mrBuffer;

	    MPI_Finalize();
	#endif

    return val;

} /* end main() */