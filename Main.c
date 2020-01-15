/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== main.c - LVB ========== */

#include "Lvb.h"

#ifndef NP_Implementation
#include <inttypes.h>
#include "Store_states.h"
#endif

	static Treestack bstack_overall;	/* overall best tree stack */
	static Treestack stack_treevo;

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


static void writeinf(Params prms, Dataptr matrix, int argc, char **argv
#ifndef NP_Implementation
, int myMPIid, int n_process
#endif
)
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
	log_Time();
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
	#ifndef NP_Implementation
	printf("MPI processes:         %d\n", n_process);
	printf("Current process ID:    %d\n", myMPIid);
	#endif
	printf("PThreads:              %d\n", prms.n_processors_available);

#ifdef MPI_Implementation
#ifndef MAP_REDUCE_SINGLE
	printf("#seeds to try           = %d\n", prms.n_seeds_need_to_try);
	printf("checkpoint interval (s) = %d\n", prms.n_checkpoint_interval);
#endif
#endif
	printf("\n================================================================================\n");	
	printf("\nInitialising search: \n");

	// if (prms.n_flag_save_read_states == DONT_SAVE_READ_STATES) printf("Don't read and save states at a specific time points\n");
	// else printf("It is going to read and save states at a specific time points\n\n");
	// printf("Output File          = %s\n", prms.file_name_out);
    // printf("Bootstrap Replicates = %ld\n", prms.bootstraps);
	//printf("    After cut =%ld\n", matrix->m);

	// Causes linker error 
	// printf("Host Machine: %s %s, %s (%d processors configured and %d processors available)\n\n", buffer.nodename, buffer.machine, buffer.sysname, get_nprocs_conf(), get_nprocs());
}

static void logtree1(Dataptr matrix, const Branch *const barray, const long start, const long cycle, long root
#ifndef NP_Implementation
, DataSeqPtr restrict matrix_seq_data
#endif
)
/* log initial tree for cycle cycle of start start (in barray) to outfp */
{
    static char outfnam[LVB_FNAMSIZE]; 	/* current file name */
    int fnamlen;			/* length of current file name */
    FILE *outfp;			/* output file */

    fnamlen = sprintf(outfnam, "%s_start%ld_cycle%ld", TREE1FNAM, start, cycle);
    lvb_assert(fnamlen < LVB_FNAMSIZE);	/* shut door if horse bolted */

    /* create tree file */
    outfp = clnopen(outfnam, "w");
	
    lvb_treeprint(matrix, outfp, barray, root
	#ifndef NP_Implementation
	, matrix_seq_data
	#endif
	);
    clnclose(outfp, outfnam);

} /* end logtree1() */

static long getsoln(Dataptr restrict matrix, Params rcstruct, int myMPIid, Lvb_bool log_progress
#ifndef NP_Implementation
, DataSeqPtr restrict matrix_seq_data
#ifdef MAP_REDUCE_SINGLE
, long *iter_p, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer
#else
, Treestack* bstack_overall, Lvb_bit_length **enc_mat
#endif
#else
, const long *weight_arr, long *iter_p
#endif
)

/* get and output solution(s) according to parameters in rcstruct;
	 * return length of shortest tree(s) found */
	{
		static char fnam[LVB_FNAMSIZE]; /* current file name */
		double t0;	/* SA cooling cycle initial temp */

		long maxpropose = MAXPROPOSE_SLOW;	/* SA cooling cycle maxpropose */
		long maxaccept = MAXACCEPT_SLOW;		/* SA cooling cycle maxaccept */
		long maxfail = MAXFAIL_SLOW;	/* SA cooling cycly maxfail */
		long treelength = LONG_MAX;		/* length of each tree found */
		long initroot;			/* initial tree's root */
		FILE *sumfp;			/* best length file */
		Branch *tree;			/* initial tree */
		long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nodes */
	    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
	    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */




		#ifndef NP_Implementation
		#ifndef MAP_REDUCE_SINGLE // MPI
	    long trees_output;		/* number of trees output for current rep. */
	    FILE *outtreefp;		/* best trees found overall */
	    long l_iterations = 0;			/* iterations of annealing algorithm */
	    int n_state_progress;		/* has a state to see if is necessary to repeat or finish */
	    int n_number_tried_seed_next = myMPIid; 	/* it has the number of tried seed */
	    int n_number_tried_seed = myMPIid;	 	/* it has the number of tried seed */
	    long cyc = 0;	/* current cycle number */
	    long start = 0;	/* current random (re)start number */

		#else // MR
		long fnamlen;			/* length of current file name */	
		long i;				/* loop counter */
		long treec;				/* number of trees found */
		FILE *resfp;			/* results file */
		Lvb_bit_length **enc_mat;		/* encoded data mat. */
		int *total_count;
		#endif

		#else // NP
	    long fnamlen;			/* length of current file name */
    	long i;				/* loop counter */
    	long treec;				/* number of trees found */
    	FILE *resfp;			/* results file */
    	Lvb_bit_length **enc_mat;	/* encoded data mat. */
		#endif

	#if defined (MAP_REDUCE_SINGLE) || defined (NP_Implementation)

	/* NOTE: These variables and their values are "dummies" and are no longer
		 * used in the current version of LVB. However, in order to keep the
		 * formatting of the output compatible with that of previous versions of
		 * LVB these variables will continue to be used and written to the summary
		 * files.  */
		long cyc = 0;	/* current cycle number */
		long start = 0;	/* current random (re)start number */

		/* dynamic "local" heap memory */
		tree = treealloc(matrix, LVB_TRUE);

		/* Allocation of the initial encoded matrix is non-contiguous because
		 * this matrix isn't used much, so any performance penalty won't matter. */
		enc_mat = (Lvb_bit_length **) malloc((matrix->n) * sizeof(Lvb_bit_length *));
		for (i = 0; i < matrix->n; i++)	enc_mat[i] = 
		#ifdef MAP_REDUCE_SINGLE
		(Lvb_bit_length *) 
		#endif
		alloc(matrix->bytes, "state sets");
		dna_makebin(matrix, 
		#ifdef MAP_REDUCE_SINGLE
		matrix_seq_data,
		#endif
		enc_mat);

		/* open and entitle statistics file shared by all cycles
		 * NOTE: There are no cycles anymore in the current version
		 * of LVB. The code bellow is purely to keep the output consistent
		 * with that of previous versions. */
		
		#ifdef MAP_REDUCE_SINGLE
		if (misc->rank == 0) {
		#endif
			if (rcstruct.verbose == LVB_TRUE) {
				sumfp = clnopen(SUMFNAM, "w");
				fprintf(sumfp, "StartNo\tCycleNo\tCycInit\tCycBest\tCycTrees\n");
			}
			else{
				sumfp = NULL;
			}
		#ifdef MAP_REDUCE_SINGLE
		}
		#endif

		/* determine starting temperature */
		randtree(matrix, tree);	/* initialise required variables */
		ss_init(matrix, tree, enc_mat);
		initroot = 0;
		#ifdef MAP_REDUCE_SINGLE
		t0 = get_initial_t(matrix, tree, rcstruct, initroot, misc->rank, log_progress);
		#endif
		#ifdef NP_Implementation
		t0 = get_initial_t(matrix, tree, rcstruct, initroot, myMPIid, log_progress, weight_arr);
		#endif
		// to = 0.01

		randtree(matrix, tree);	/* begin from scratch */
		ss_init(matrix, tree, enc_mat);
		initroot = 0;

		#ifdef MAP_REDUCE_SINGLE
		if (misc->rank == 0) {
		#endif
			if (rcstruct.verbose) smessg(start, cyc);
			check_stdout();
		#ifdef MAP_REDUCE_SINGLE
		}
		#endif

		/* start cycles's entry in sum file
		 * NOTE: There are no cycles anymore in the current version
		 * of LVB. The code bellow is purely to keep the output consistent
		 * with that of previous versions.  */
		#ifdef MAP_REDUCE_SINGLE
		if (misc->rank == 0) {
		#endif
			if(rcstruct.verbose == LVB_TRUE) {
				alloc_memory_to_getplen(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
				#ifdef MAP_REDUCE_SINGLE
				fprintf(sumfp, "%ld\t%ld\t%ld\t", start, cyc, getplen(matrix, tree, rcstruct, initroot, p_todo_arr, p_todo_arr_sum_changes, p_runs));
				#endif
				#ifdef NP_Implementation
				fprintf(sumfp, "%ld\t%ld\t%ld\t", start, cyc, getplen(matrix, tree, rcstruct, initroot, p_todo_arr, p_todo_arr_sum_changes, p_runs, weight_arr));
				#endif
				free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
				logtree1(matrix, tree, start, cyc, initroot
				#ifdef MAP_REDUCE_SINGLE
				, matrix_seq_data
				#endif
				);
			}
		#ifdef MAP_REDUCE_SINGLE
		}

		MPI_Barrier(MPI_COMM_WORLD);
		/* find solution(s) */
		maxaccept = get_random_maxaccept();
		// printf("\nmaxaccept:%ld\n", maxaccept);
		treelength = anneal(matrix, &bstack_overall, &stack_treevo, tree, rcstruct, &rcstruct, initroot, t0, maxaccept,
				maxpropose, maxfail, stdout, iter_p, myMPIid, log_progress, misc, mrTreeStack, mrBuffer );
		#endif

		#ifdef NP_Implementation
		treelength = anneal(matrix, &bstack_overall, &stack_treevo, tree, rcstruct, &rcstruct, initroot, t0, maxaccept,
    	maxpropose, maxfail, stdout, iter_p, myMPIid, log_progress, weight_arr);
		treestack_pop(matrix, tree, &initroot, &bstack_overall, LVB_FALSE);
		treestack_push(matrix, &bstack_overall, tree, initroot, LVB_FALSE); //identical
		#endif

		#ifdef MAP_REDUCE_SINGLE
		long val = treestack_pop(matrix, tree, &initroot, &bstack_overall, LVB_FALSE);
		treestack_push(matrix, &bstack_overall, tree, initroot, LVB_FALSE);

		if(val ==  1) {
			misc->SB = 0;
			tree_setpush(matrix, tree, initroot, mrBuffer, misc);
			mrTreeStack->add(mrBuffer);
			mrTreeStack->collate(NULL);
			mrTreeStack->reduce(reduce_filter, NULL);

			mrBuffer->add(mrTreeStack);
			mrBuffer->collate(NULL);

			misc->count = (int *) alloc( (bstack_overall.next+1) * sizeof(int), "integer array for tree compare using MapReduce");
			total_count = (int *) alloc( (bstack_overall.next+1) * sizeof(int), "integer array for tree compare using MapReduce");
			for(int i=0; i<=bstack_overall.next; i++) misc->count[i] = 0;
			mrBuffer->reduce(reduce_count, misc);

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
		//	  treestack_push_only(matrix, &bstack_overall, tree, initroot);
			  misc->ID = bstack_overall.next;
				  misc->SB = 1;
				  tree_setpush(matrix, tree, initroot, mrBuffer, misc);
				  mrTreeStack->add(mrBuffer);
			}

			free(misc->count);
			free(total_count);
		//} else {
		//	treestack_push_only(matrix, &bstack_overall, tree, initroot);
		//	misc->ID = bstack_overall.next;
		//	misc->SB = 1;
		//        tree_setpush(matrix, tree, initroot, mrBuffer, misc);
		//        mrTreeStack->add(mrBuffer);
		}

		treelength = deterministic_hillclimb(matrix, &bstack_overall, tree, rcstruct, initroot, stdout,
				iter_p, myMPIid, log_progress, misc, mrTreeStack, mrBuffer);

		#endif

		/* log this cycle's solution and its details
		 * NOTE: There are no cycles anymore in the current version
		 * of LVB. The code bellow is purely to keep the output consistent
		 * with that of previous versions. */

		#ifdef MAP_REDUCE_SINGLE
	if (misc->rank == 0) {
		#endif
		if (rcstruct.verbose == LVB_TRUE){
			fnamlen = sprintf(fnam, "%s_start%ld_cycle%ld", RESFNAM, start, cyc);
			lvb_assert(fnamlen < LVB_FNAMSIZE);	/* really too late */
			resfp = clnopen(fnam, "w");
		#ifdef MAP_REDUCE_SINGLE
			treec = treestack_print(matrix, matrix_seq_data, &bstack_overall, resfp, LVB_FALSE);
		#endif

		#ifdef NP_Implementation
		treec = treestack_print(matrix, &bstack_overall, resfp, LVB_FALSE);
		#endif
			clnclose(resfp, fnam);
			fprintf(sumfp, "%ld\t%ld\n", treelength, treec);

			/* won't use length summary file until end of next cycle */
			fflush(sumfp);
			if (ferror(sumfp)){
				crash("write error on file %s", SUMFNAM);
			}
		}

		if (rcstruct.verbose == LVB_TRUE) printf("Ending start %ld cycle %ld\n", start, cyc);
		check_stdout();

		if (rcstruct.verbose == LVB_TRUE) clnclose(sumfp, SUMFNAM);
		#ifdef MAP_REDUCE_SINGLE
	}
		#endif
		/* "local" dynamic heap memory */
		free(tree);
		for (i = 0; i < matrix->n; i++) free(enc_mat[i]);
		free(enc_mat);

		return treelength;

	} /* end getsoln() */

	#else
	do{
			/* dynamic "local" heap memory */
			tree = treealloc(matrix, LVB_TRUE);

			/* open and entitle statistics file shared by all cycles
			 * NOTE: There are no cycles anymore in the current version
			 * of LVB. The code bellow is purely to keep the output consistent
			 * with that of previous versions. */

			if (rcstruct.verbose == LVB_TRUE) {
				sumfp = clnopen(SUMFNAM, "w");
				fprintf(sumfp, "Process\tStartNo\tCycleNo\tCycInit\tCycBest\tCycTrees\n");
			}
			else{
				sumfp = NULL;
			}

			/* IMPORTANT: always necessary to initialize this tree, even if it is going to read a checkpoint file */
			randtree(matrix, tree);	/* initialise required variables */
			ss_init(matrix, tree, enc_mat);

			/* is going to start from beginning */
			if (rcstruct.n_flag_is_possible_read_state_files != CHECK_POINT_READ_STATE_FILES){
				/* determine starting temperature because it will be start from begin*/
				initroot = 0;
				t0 = get_initial_t(matrix, tree, rcstruct, initroot, myMPIid, log_progress);
			//    t0 = 0.18540001000004463;

				randtree(matrix, tree);	/* begin from scratch */
				ss_init(matrix, tree, enc_mat);
				initroot = 0;

				if (rcstruct.verbose) smessg(start, cyc);
					check_stdout();

				/* start cycles's entry in sum file
				 * NOTE: There are no cycles anymore in the current version
				 * of LVB. The code bellow is purely to keep the output consistent
				 * with that of previous versions.  */
				if(rcstruct.verbose == LVB_TRUE) {
					alloc_memory_to_getplen(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
					fprintf(sumfp, "%d\t%ld\t%ld\t%ld\t", myMPIid, start, cyc, getplen(matrix, tree, rcstruct, initroot, p_todo_arr, p_todo_arr_sum_changes, p_runs));
					free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
					logtree1(matrix, tree, start, cyc, initroot, matrix_seq_data);
				}
			}
			else{
				t0 = 0; /* it's only a problem of initialization
				 	 	 	 it's going to be read from anneal checkpoint file*/
			}


			/* find solution(s) */
			n_state_progress = 0;	/* there's no state in beginning */
			n_number_tried_seed = n_number_tried_seed_next;
			maxaccept = get_random_maxaccept();
			printf("\nProcess:%d    maxaccept:%ld\n", myMPIid, maxaccept);
			treelength = anneal(matrix, bstack_overall, &stack_treevo, tree, rcstruct, &rcstruct, initroot, t0, maxaccept,
					maxpropose, maxfail, stdout, &l_iterations, myMPIid, log_progress, &n_state_progress, &n_number_tried_seed_next);

	/* 		Several possible outputs */
	/*		ANNEAL_FINISHED_AND_NOT_REPEAT		0x01
			ANNEAL_FINISHED_AND_REPEAT			0x02
			ANNEAL_KILLED_AND_REPEAT			0x03
			ANNEAL_KILLED_AND_NOT_REPEAT		0x04 */

			/* is is killed is not necesary any data */
			if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT || n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT){

				/* work on the trees */
				int l_pop = treestack_pop(matrix, tree, &initroot, bstack_overall, LVB_FALSE);
				if (l_pop == 0){
					printf("\nProcess:%d    Error: can't pop any tree from treestack.   Rearrangements tried: %ld\n", myMPIid, l_iterations);
				}
				else{
					treestack_push(matrix, bstack_overall, tree, initroot, LVB_FALSE);
					treelength = deterministic_hillclimb(matrix, bstack_overall, tree, rcstruct, initroot, stdout, &l_iterations, myMPIid, log_progress);
					/* save it */
					sprintf(fnam, "%s_%d", rcstruct.file_name_out, n_number_tried_seed); /* name of output file for this process */
					outtreefp = clnopen(fnam, "w");
					trees_output = treestack_print(matrix, matrix_seq_data, bstack_overall, outtreefp, LVB_FALSE);
					clnclose(outtreefp, fnam);
					printf("\nProcess:%d   Rearrangements tried: %ld\n", myMPIid, l_iterations);
					if (trees_output == 1L) { printf("1 most parsimonious tree of length %ld written to file '%s'\n", treelength, fnam); }
					else { printf("%ld equally parsimonious1 trees of length %ld written to file '%s'\n", trees_output, treelength, fnam); }
				}
			}
			if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT || n_state_progress == MESSAGE_ANNEAL_KILLED_AND_REPEAT){
				l_iterations = 0;		/* start iterations from zero */
				free(tree);
				treestack_free(bstack_overall);
				printf("Process:%d   try seed number process:%d   new seed:%d", myMPIid, n_number_tried_seed_next, rcstruct.seed);
				rinit(rcstruct.seed); /* at this point the structure has a need see passed by master process */
			}
			else{
				/* Save the finish state file */
				if (rcstruct.n_flag_save_read_states == DO_SAVE_READ_STATES){
					save_finish_state_file(&rcstruct, myMPIid);
				}
				break; /* it is not necessary to repeat again */
			}

			check_stdout();
			if (rcstruct.verbose == LVB_TRUE) printf("Ending start %ld cycle %ld\n", start, cyc);
			if (rcstruct.verbose == LVB_TRUE) clnclose(sumfp, SUMFNAM);
	    } while (1);

	    /* "local" dynamic heap memory */
	    free(tree);
	    return treelength;

	} /* end getsoln() */

	#endif

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
		#ifndef NP_Implementation
		matrix->n_slice_size_getplen = 0;	/* it doens't matter this value for 1 thread */
		#endif
		matrix->n_threads_getplen = 1; /* need to pass for 1 thread because the number of words is to low */
	}
	#ifndef NP_Implementation
}

/* TODO, need to control the seeds that where processed*/
int get_other_seed_to_run_a_process(){
	return (int) (rand() % (unsigned long) MAX_SEED);
	#else
	if (matrix->n_threads_getplen < 1) matrix->n_threads_getplen = 1;
	#endif
}

#ifndef NP_Implementation
#ifndef MAP_REDUCE_SINGLE
	void print_data(Dataptr p_lvbmat, int n_thread){
		printf("###############################\n thread: %d\n", n_thread);
		printf("n_threads_getplen: %d\n", p_lvbmat->n_threads_getplen);
		printf("n_slice_size_getplen: %d\n", p_lvbmat->n_slice_size_getplen);
		printf("m: %ld\n", p_lvbmat->m);
		printf("n: %ld\n", p_lvbmat->n);

		printf("max_length_seq_name: %ld\n", p_lvbmat->max_length_seq_name);
		printf("min_len_tree: %ld\n", p_lvbmat->min_len_tree);
	}

	void print_data_seq(DataSeqPtr p_lvbmat, int n_size, int n_thread){

		printf("########### SEQ ####################\n thread: %d\n", n_thread);
		int i;
		for(i = 0; i < n_size; ++i){
			printf("%s\n", p_lvbmat->rowtitle[i]);
			printf("%s  %s\n", p_lvbmat->rowtitle[i], p_lvbmat->row[i]);
		}
	}

	void print_binary_data(Lvb_bit_length **p_enc_data, int n_size, int nwords, int n_thread){

		printf("########### SEQ ####################\n thread: %d\n", n_thread);
		int i, x;
		for(i = 0; i < n_size; ++i){
			printf("	rank:%d  u: ", i);
			for(x = 0; x < nwords; x++){ printf("%" PRIu64, p_enc_data[i][x]); }
			printf("\n");
		}
	}

	/*  get temperatures from processes */
	void get_temperature_and_control_process_from_other_process(int num_procs, int n_seeds_to_try){

		int nProcessFinished = 0;
		MPI_Request *pHandleTemperatureRecv;
		MPI_Request *pHandleFinishProcess;
		MPI_Request *pHandleManagementProcess;
		MPI_Status mpi_status;
		IterationTemperature *p_calc_iterations;
		int *pIntFinishedProcessChecked, *pIntFinishedProcess, *pIntProcessControlNumberRunning;
		int nFlag, i, n_seeds_tried = 0;

		pHandleTemperatureRecv = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for temperature...");
		pHandleFinishProcess = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for finish process...");
		pHandleManagementProcess = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for management process...");
		pIntProcessControlNumberRunning = (int *) alloc(num_procs * sizeof(int), "Buffer with index between process and tried process, but the ones that are running...");
		pIntFinishedProcess = (int *) alloc(num_procs * sizeof(int), "Buffer for finished processes...");
		pIntFinishedProcessChecked = (int *) alloc(num_procs * sizeof(int), "Buffer for processes states...");
		memset(pIntFinishedProcessChecked, 0, num_procs * sizeof(int)); 	/* control the state of a process, but for all seeds tried */
		memset(pIntProcessControlNumberRunning, 0, num_procs * sizeof(int));		/* has the index number of the process that is running */
		memset(pIntFinishedProcess, 0, num_procs * sizeof(int));		/* control if a specific process is finished */
		for (i = 0; i < num_procs; i++) { *(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN; } /* all the process start with this state */

		/* structure to use sending temperature and number of interactions to master process */
		int				nItems = 3;
		int          	blocklengths[3] = {2, 1, 1};
		MPI_Datatype 	types[3] = {MPI_INT, MPI_LONG, MPI_DOUBLE};
		MPI_Datatype 	mpi_recv_data;
		MPI_Aint     	displacements[3];
		displacements[0] = offsetof(SendInfoToMaster, n_iterations);
		displacements[1] = offsetof(SendInfoToMaster, l_length);
		displacements[2] = offsetof(SendInfoToMaster, temperature);
		MPI_Type_create_struct(nItems, blocklengths, displacements, types, &mpi_recv_data);
		MPI_Type_commit(&mpi_recv_data);

		nItems = 1;
		int          	blocklengths_2[1] = {3};
		MPI_Datatype 	types_2[1] = {MPI_INT};
		MPI_Datatype 	mpi_send_data;
		MPI_Aint     	displacements_2[1];
		displacements_2[0] = offsetof(RecvInfoFromMaster, n_seed);
		MPI_Type_create_struct(nItems, blocklengths_2, displacements_2, types_2, &mpi_send_data);
		MPI_Type_commit(&mpi_send_data);

		/* structure with temperature and interaction data */
		SendInfoToMaster **p_info_temperature;
		p_info_temperature = (SendInfoToMaster **) malloc(sizeof(SendInfoToMaster *) * (num_procs - 1));
		for (i = 0; i < num_procs - 1; i++){
			*(p_info_temperature + i) = malloc(sizeof(SendInfoToMaster));
		}

		/* structure with management data */
		RecvInfoFromMaster **p_info_manage;
		p_info_manage = (RecvInfoFromMaster **) malloc(sizeof(RecvInfoFromMaster *) * (num_procs - 1));
		for (i = 0; i < num_procs - 1; i++){
			*(p_info_manage + i) = malloc(sizeof(RecvInfoFromMaster));
		}

		/* first get handles for all receivers */
		for (i = 1; i < num_procs; i++) {
			MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
			MPI_Irecv(pIntFinishedProcess + i, 1, MPI_INT, i, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, pHandleFinishProcess + i);
			*(pHandleManagementProcess + i) = 0;
			*(pIntProcessControlNumberRunning + i) = i - 1;   /* the index number of ones that are running */
		}

		/* alloc memory main calc iterations*/
		p_calc_iterations = get_alloc_main_calc_iterations();
		n_seeds_tried = num_procs - 1;	/* we try always these number of seeds in beginning */
		while (1){

			for (i = 1; i < num_procs; i++) {

				/* if is equal to ANNEAL_STOP_PROCESS it is not necessary to do anything else to this process */
				/* because the limit of tried process is reached */
				if (*(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE && *(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS){

					MPI_Test(pHandleTemperatureRecv + i, &nFlag, &mpi_status);
					if (nFlag == 1){	/* message received */
						printf("Process:%d    main process getting temperature and length from process:%d   temperature:%-15.8g   length:%ld   iteration:%d\n", mpi_status.MPI_SOURCE,
							i, (*(p_info_temperature + (i - 1)))->temperature, (*(p_info_temperature + (i - 1)))->l_length,
							    (*(p_info_temperature + (i - 1)))->n_iterations);

						/* all the processes are synchronized by the number of iterations... */
						/* upload the data to the structure */
						/* make calculations to how many standard deviations has */
						/* if it is less than X standard deviations can continue */
						/* add data from this process to structure calc temperature iterations */
						add_temperature_cal_iterations(p_calc_iterations, *(p_info_temperature + (i - 1)), *(pIntProcessControlNumberRunning + i));
						if (is_possible_to_continue(p_calc_iterations, (*(p_info_temperature + (i - 1)))->temperature,
								(*(p_info_temperature + (i - 1)))->n_iterations, (*(p_info_temperature + (i - 1)))->l_length, num_procs, 0)){

							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_CONTINUE_ANNEAL;
							/* need to check this one, to see if the previous message was delivered... */
							if (*(pHandleManagementProcess + i) != 0){
								/* printf("Process:%d   wait management from process:%d\n", 0, i); */
								MPI_Wait(pHandleManagementProcess + i, MPI_STATUS_IGNORE);
							}
						}
						else{
							/* if the previous one was not delivered it can be canceled... */
							if (*(pHandleManagementProcess + i) != 0){
								MPI_Test(pHandleManagementProcess + i, &nFlag, &mpi_status);
								if (nFlag == 0) MPI_Cancel(pHandleManagementProcess + i);
							}

							*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE;
							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;
						}

						/* send management message */
						/*printf("Process:%d   send management to process:%d\n", 0, i); */
						MPI_Isend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);

						/* launch other receive temperature */
						/* printf("Process:%d   receive temperature from process:%d\n", 0, i); */
						if ((*(p_info_manage + (i - 1)))->n_is_to_continue != MPI_IS_TO_RESTART_ANNEAL)
							MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
					}
				}

				/* test the ones that finish, and all of them need to send this message...  */
				if (*(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS){
					MPI_Test(pHandleFinishProcess + i, &nFlag, &mpi_status);
					if (nFlag == 1 && *(pIntFinishedProcess + i) == MPI_FINISHED){
						MPI_Test(pHandleTemperatureRecv + i, &nFlag, &mpi_status);
						if (nFlag == 0) MPI_Cancel(pHandleTemperatureRecv + i);
						if (*(pHandleManagementProcess + i) != 0) MPI_Cancel(pHandleManagementProcess + i);
						printf("Process:%d    finish\n", i);
						if (n_seeds_tried < n_seeds_to_try){		/* try another seed */
							n_seeds_tried += 1;									/* try another seed */
							*(pIntProcessControlNumberRunning + i) = n_seeds_tried;
							*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN;

							(*(p_info_manage + (i - 1)))->n_seed = get_other_seed_to_run_a_process();
							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;
							(*(p_info_manage + (i - 1)))->n_process_tried = n_seeds_tried;
							printf("Process:%d    send to restart\n", i);
							MPI_Send(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_SEND_RESTART, MPI_COMM_WORLD);

							/* launch other wait messages */
							MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
							MPI_Irecv(pIntFinishedProcess + i, 1, MPI_INT, i, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, pHandleFinishProcess + i);
						}
						else{
							*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS;
							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_NOT_TO_RESTART;
							printf("Process:%d    isn't to restart\n", i);
							MPI_Send(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_SEND_RESTART, MPI_COMM_WORLD);

							/* there's no need other MPI messages */
							*(pHandleTemperatureRecv + i) = 0;
							*(pHandleFinishProcess + i) = 0;
						}
						nProcessFinished += 1;
					}
				}
			}

			/* is it finished, all of then?  */
			if (nProcessFinished == n_seeds_to_try){
				/* All processes are finished */
				/* cancel the requests of temperature */
				for (i = 1; i < num_procs; i++){
					if (*(pHandleTemperatureRecv + i) != 0) MPI_Cancel(pHandleTemperatureRecv + i);
					if (*(pHandleManagementProcess + i) != 0) MPI_Cancel(pHandleManagementProcess + i);
				}
				break;
			}
		}

		/* release memory main calc iterations */
		release_main_calc_iterations(p_calc_iterations);

		/* release other memory */
		for (i = 0; i < num_procs - 1; i++){ free(*(p_info_temperature + i)); }
		free(p_info_temperature);
		for (i = 0; i < num_procs - 1; i++){ free(*(p_info_manage + i)); }
		free(p_info_manage);

		free(pHandleTemperatureRecv);		/* all asynchronous messages */
		free(pHandleManagementProcess);		/* all asynchronous messages */
		free(pHandleFinishProcess);			/* all asynchronous messages */
		free(pIntFinishedProcess);
		free(pIntProcessControlNumberRunning);
		free(pIntFinishedProcessChecked);
	}

#endif
#endif

	int main(int argc, char **argv)
	{
		Dataptr matrix;	/* data matrix */
		Params rcstruct;		/* configurable parameters */
		Lvb_bool log_progress;	/* whether or not to log anneal search */
		int myMPIid;

	#if defined (MAP_REDUCE_SINGLE) || defined (NP_Implementation)

	int val;			/* return value */
	long iter;			/* iterations of annealing algorithm */
	long trees_output_total = 0L;	/* number of trees output, overall */
	long trees_output;		/* number of trees output for current rep. */
	long final_length;		/* length of shortest tree(s) found */
	FILE *outtreefp;		/* best trees found overall */

	#endif

#ifndef NP_Implementation
#ifndef MAP_REDUCE_SINGLE
	    
	    DataSeqPtr matrix_seq_data;
	    		/* configurable parameters */
	    long i, n_buffer_size_matrix, n_buffer_size_binary;			/* loop counter */
	    int position, n_error_code = EXIT_SUCCESS; /* return value */
	    Lvb_bit_length **enc_mat = NULL;/* encoded data mat. */
	    
	    Treestack *p_bstack_overall = NULL;	/* overall best tree stack */
	    char *pack_data = NULL;
	    Lvb_bit_length *p_pack_data_binary = NULL;
	    /* global files */

	    /* define mpi */
	    MPI_Status	status;		/* return status for receive */
	    int num_procs, mpi_err;
	    mpi_err = MPI_Init(&argc, &argv);
	    mpi_err |= MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
		mpi_err |= MPI_Comm_rank(MPI_COMM_WORLD, &myMPIid);
		if (mpi_err != 0) {
			printf("Error: MPI initialization failed\n");
			exit(1);
		}
		if (num_procs < 2) {
			printf("Number of process must be greater than 1 - try again.\n");
			MPI_Finalize();
			exit(1);
		}

		matrix = alloc(sizeof(DataStructure), "alloc data structure");
		matrix_seq_data = alloc(sizeof(DataSeqStructure), "alloc data structure");
		matrix_seq_data->row = NULL;
		matrix_seq_data->rowtitle = NULL;

		/* define the matrix structure MPI */
		int				nItems = 2;
		int          	blocklengths[2] = {12, 2};
		MPI_Datatype 	types[2] = {MPI_LONG, MPI_INT};
		MPI_Datatype 	mpi_matrix;
		MPI_Aint     	displacements[2];
		displacements[0] = offsetof(DataStructure, m);
		displacements[1] = offsetof(DataStructure, n_threads_getplen);
		MPI_Type_create_struct(nItems, blocklengths, displacements, types, &mpi_matrix);
		MPI_Type_commit(&mpi_matrix);

		nItems = 3;
		int          	blocklengthsPar[3] = {1, 10, LVB_FNAMSIZE + LVB_FNAMSIZE};
		MPI_Datatype 	typesPar[3] = {MPI_LONG, MPI_INT, MPI_CHAR};
		MPI_Datatype 	mpi_params;
		MPI_Aint     	displacementsPar[3];
		displacementsPar[0] = offsetof(Params, verbose);
		displacementsPar[1] = offsetof(Params, seed);
		displacementsPar[2] = offsetof(Params, file_name_in);
		MPI_Type_create_struct(nItems, blocklengthsPar, displacementsPar, typesPar, &mpi_params);
		MPI_Type_commit(&mpi_params);
		/* end define mpi */

		if (myMPIid == MPI_MAIN_PROCESS){
#else
		/* MapReduce version */
		MPI_Init(&argc,&argv);
		MISC misc;
		MPI_Comm_rank(MPI_COMM_WORLD,&misc.rank);
		MPI_Comm_size(MPI_COMM_WORLD,&misc.nprocs);
		DataSeqPtr matrix_seq_data;
	    int n_error_code = EXIT_SUCCESS; /* return value */
		MapReduce *mrTreeStack = new MapReduce(MPI_COMM_WORLD);
		mrTreeStack->memsize = 1024;
		mrTreeStack->verbosity = 0;
		mrTreeStack->timer = 0;
		MapReduce *mrBuffer = new MapReduce(MPI_COMM_WORLD);
		mrBuffer->memsize = 1024;
		mrBuffer->verbosity = 0;
		mrBuffer->timer = 0;

		if(misc.rank == 0) {
#endif
#else
	
    long i;			/* loop counter */
    long replicate_no = 0L;	/* current bootstrap replicate number */
    double total_iter = 0.0;	/* total iterations across all replicates */
    long *weight_arr;  		/* weights for sites */
#endif

	// entitle standard output
	#if defined (MAP_REDUCE_SINGLE) || defined (NP_Implementation)
		print_LVB_COPYRIGHT();
		print_LVB_INFO();
	#endif

#ifdef MAP_REDUCE_SINGLE
		}
#endif
	/* start timer */
    clock_t Start, End;
   	double Overall_Time_taken;
    Start = clock();
	lvb_initialize();

#if defined (MAP_REDUCE_SINGLE) || defined (NP_Implementation)

myMPIid = 0;

#ifdef MAP_REDUCE_SINGLE

		n_error_code = getparam(&rcstruct, argc, argv);
	    /* read and alloc space to the data structure */
	    matrix = (data *) alloc(sizeof(DataStructure), "alloc data structure");
	    matrix_seq_data = (seq_data *) alloc(sizeof(DataSeqStructure), "alloc data structure");
	    matrix_seq_data->row = NULL;
	    matrix_seq_data->rowtitle = NULL;
	    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, matrix, matrix_seq_data);
	    /* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
	    bstack_overall = *(treestack_new());
		if(rcstruct.algorithm_selection ==2)
 	 	stack_treevo = *(treestack_new());
	    matchange(matrix, matrix_seq_data, rcstruct);	/* cut columns */
            writeinf(rcstruct, matrix, argc, argv, misc.rank, misc.nprocs);
	    calc_distribution_processors(matrix, rcstruct);

	    if (rcstruct.verbose == LVB_TRUE) {
	    	if (misc.rank == 0) printf("getminlen: %ld\n\n", matrix->min_len_tree);
	    }
	    rinit(rcstruct.seed);
            log_progress = LVB_TRUE;
				
	    if (misc.rank == 0) outtreefp = clnopen(rcstruct.file_name_out, "w");
		FILE * treEvo;
		if(rcstruct.algorithm_selection ==2)
		treEvo = fopen ("treEvo.tre","w");
	    iter = 0;
		final_length = getsoln(matrix, rcstruct,  myMPIid, log_progress, matrix_seq_data, &iter, &misc, mrTreeStack, mrBuffer);
	    if (misc.rank == 0) {
	       trees_output = treestack_print(matrix, matrix_seq_data, &bstack_overall, outtreefp, LVB_FALSE);
	    }
	    trees_output_total += trees_output;
	    treestack_clear(&bstack_overall);

	    /* clean the TreeStack and buffer */
	    mrTreeStack->map( mrTreeStack, map_clean, NULL );
	    mrBuffer->map( mrBuffer, map_clean, NULL );
	    /* END clean the TreeStack and buffer */

	    if (misc.rank == 0) {
#else

    getparam(&rcstruct, argc, argv);

    /* read and alloc space to the data structure */
    matrix = alloc(sizeof(DataStructure), "alloc data structure");
    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, matrix);

    /* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
    bstack_overall = treestack_new();
    if(rcstruct.algorithm_selection ==2)
    stack_treevo = treestack_new();

    matchange(matrix, rcstruct);	/* cut columns */
    writeinf(rcstruct, matrix, argc, argv);
    calc_distribution_processors(matrix, rcstruct);

    if (rcstruct.verbose == LVB_TRUE) {
    	printf("Based on matrix provided, maximum parsimony tree length: %ld\n\n", getminlen(matrix));
    }
    rinit(rcstruct.seed);
    if (rcstruct.bootstraps > 0) {
    	log_progress = LVB_TRUE;
    	printf("Temperature:   Rearrangement: TreeStack size: Length:\n");
    }
    else log_progress = LVB_TRUE;

    weight_arr = (long*) alloc(sizeof(long) * matrix->m, "alloc data structure");
    outtreefp = clnopen(rcstruct.file_name_out, "w");
    FILE * treEvo;
    if(rcstruct.algorithm_selection ==2)
    treEvo = fopen ("treEvo.tre","w");
    do{
		iter = 0;
		if (rcstruct.bootstraps > 0){
			get_bootstrap_weights(weight_arr, matrix->m, matrix->original_m - matrix->m);
		}
		else{
			for (i = 0; i < matrix->m; i++) weight_arr[i] = 1;    // logstim();    // logstim();ss, weight_arr, &iter);
		}

		final_length = getsoln(matrix, rcstruct,  myMPIid, log_progress, weight_arr, &iter);

		if (rcstruct.bootstraps > 0) trees_output = treestack_print(matrix, &bstack_overall, outtreefp, LVB_TRUE);
		else trees_output = treestack_print(matrix, &bstack_overall, outtreefp, LVB_FALSE);

		trees_output_total += trees_output;
        if(rcstruct.algorithm_selection ==2)
		treestack_print(matrix, &stack_treevo, treEvo, LVB_FALSE);
        treestack_clear(&bstack_overall);
		replicate_no++;
		if (rcstruct.bootstraps > 0) {
			printf("\n\nReplicate %ld complete:\n\nRearrangements tried: %-16ld\nTrees saved:          %-16ld\nLength:               %ld\n\n", replicate_no, iter, trees_output, final_length);
            if (replicate_no < rcstruct.bootstraps)
            printf("Temperature:   Rearrangement: TreeStack size: Length:\n");
			total_iter += (double) iter;
		}
} while (replicate_no < rcstruct.bootstraps);

#endif


	if(rcstruct.algorithm_selection ==2)
		fclose(treEvo);
	    clnclose(outtreefp, rcstruct.file_name_out);

	End = clock();
	Overall_Time_taken = ((double) (End - Start)) /CLOCKS_PER_SEC;

	if (logfile_exists ("logfile.tsv"))
	{
		FILE * logfile;
    	logfile = fopen ("logfile.tsv","a+");
		fprintf (logfile, "%ld\t%ld\t%ld\t%.2lf\n", iter, trees_output_total, final_length, Overall_Time_taken);
		fclose(logfile);
	}
	else {
		FILE * logfile;
	    logfile = fopen ("logfile.tsv","a+");
		fprintf (logfile, "Rearrangments\tTopologies\tScore\tRuntime\n");
		fprintf (logfile, "%ld\t%ld\t%ld\t%.2lf\n", iter, trees_output_total, final_length, Overall_Time_taken);
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

#ifdef MAP_REDUCE_SINGLE
	}
#endif

/* "file-local" dynamic heap memory */

#ifdef MAP_REDUCE_SINGLE

if(rcstruct.algorithm_selection ==2)
treestack_free(&stack_treevo);
treestack_free(&bstack_overall);
rowfree(matrix_seq_data, matrix->n);
free(matrix);

MPI_Barrier(MPI_COMM_WORLD);

delete mrTreeStack;
delete mrBuffer;

MPI_Finalize();
#else
if (rcstruct.algorithm_selection ==2)
treestack_free(matrix, &stack_treevo);
treestack_free(matrix, &bstack_overall);
rowfree(matrix);
free(matrix);
free(weight_arr); 
#endif

if (cleanup() == LVB_TRUE) val = EXIT_FAILURE;
    else val = EXIT_SUCCESS;

return val;

} /* end main() */

#else

			n_error_code = getparam(&rcstruct, argc, argv);
			if (n_error_code == EXIT_SUCCESS){
				/* read and alloc space to the data structure */
				n_error_code = phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, matrix, matrix_seq_data);
				if (n_error_code == EXIT_SUCCESS){

					/* send message to continue */
					MPI_Bcast(&n_error_code, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_COMM_WORLD);

					/* start the main process */
					matchange(matrix, matrix_seq_data, rcstruct);	/* cut columns */
					if (rcstruct.n_seeds_need_to_try < (num_procs - 1)) { rcstruct.n_seeds_need_to_try = num_procs - 1; }
					writeinf(rcstruct, matrix, myMPIid, num_procs, argc, argv);
					calc_distribution_processors(matrix, rcstruct);

					/* Allocation of the initial encoded matrix is non-contiguous because
					 * this matrix isn't used much, so any performance penalty won't matter. */
					enc_mat = (Lvb_bit_length **) alloc((matrix->n) * sizeof(Lvb_bit_length *), "state sets");
					for (i = 0; i < matrix->n; i++) enc_mat[i] = alloc(matrix->bytes, "state sets");
					dna_makebin(matrix, matrix_seq_data, enc_mat);

					/* pack data for seq matrix */
			#ifdef MPI_SEND_ONLY_MATRIX_NAMES
					n_buffer_size_matrix = sizeof(char) * matrix->n * (matrix->max_length_seq_name + 1);
			#else
					n_buffer_size_matrix = sizeof(char) * matrix->n * (matrix->m + 1) + sizeof(char) * matrix->n * (matrix->max_length_seq_name + 1);
			#endif
					pack_data = (char *) alloc(n_buffer_size_matrix, "sequences packing");
					/* packing data */
					position = 0;
					for (i = 0; i < matrix->n; i++){
			#ifndef MPI_SEND_ONLY_MATRIX_NAMES
						MPI_Pack(matrix_seq_data->row[i], matrix->m + 1, MPI_CHAR, pack_data, n_buffer_size_matrix, &position, MPI_COMM_WORLD);
			#endif
						MPI_Pack(matrix_seq_data->rowtitle[i], matrix->max_length_seq_name + 1, MPI_CHAR, pack_data, n_buffer_size_matrix, &position, MPI_COMM_WORLD);
					}
					/* END pack data for seq matrix */

					/* pack binary for seq data */
					n_buffer_size_binary = matrix->bytes * matrix->n;
					p_pack_data_binary = (Lvb_bit_length *) alloc(n_buffer_size_binary, "binary packing");
					position = 0;
					for (i = 0; i < matrix->n; i++){
						MPI_Pack(enc_mat[i], matrix->bytes, MPI_CHAR, p_pack_data_binary, n_buffer_size_binary, &position, MPI_COMM_WORLD);
					}
					/* END pack binary for seq data */

					for (i = 1; i < num_procs; i++) {

						/* send matrix to each process */
						MPI_Send(matrix, 1, mpi_matrix, i, MPI_TAG_MATRIX, MPI_COMM_WORLD); //, &request);

						/* send binary data of each sample */
						rcstruct.seed = get_other_seed_to_run_a_process();
						MPI_Send(&rcstruct, 1, mpi_params, i, MPI_TAG_PARAMS, MPI_COMM_WORLD); //, &request);

						/* send names and/or seq and names */
						MPI_Send(pack_data, n_buffer_size_matrix, MPI_PACKED, i, MPI_TAG_NAME_AND_SEQ_DATA, MPI_COMM_WORLD); //, &request);

						/* send binary data */
						MPI_Send(p_pack_data_binary, n_buffer_size_binary, MPI_PACKED, i, MPI_TAG_BINARY_DATA, MPI_COMM_WORLD); //, &request);
					}

					/*  get temperatures from processes */
					/* and wait until all processes are finished */
					get_temperature_and_control_process_from_other_process(num_procs, rcstruct.n_seeds_need_to_try);

					/* clean memory */
					for (i = 0; i < matrix->n; i++) free(enc_mat[i]);
					free(enc_mat);
					rowfree(matrix_seq_data, matrix->n);
					free(pack_data);
					free(p_pack_data_binary);
					free(p_bstack_overall);
				}
			}
			if (n_error_code != EXIT_SUCCESS){
				/* kill the others */
				MPI_Bcast(&n_error_code, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_COMM_WORLD);
			}
		}
		else{

			/* Is to proceed or to finish? */
			MPI_Bcast(&n_error_code, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_COMM_WORLD);
			if (n_error_code == EXIT_SUCCESS){

				MPI_Recv(matrix, 1, mpi_matrix, MPI_MAIN_PROCESS, MPI_TAG_MATRIX, MPI_COMM_WORLD, &status);
				//print_data(matrix, myMPIid); /* print data.... */

				MPI_Recv(&rcstruct, 1, mpi_params, MPI_MAIN_PROCESS, MPI_TAG_PARAMS, MPI_COMM_WORLD, &status);
				//writeinf(rcstruct, matrix, myMPIid, num_procs); /* print data.... */

				matrix_seq_data->rowtitle = (char **) malloc((matrix->n) * sizeof(char *));
		#ifdef MPI_SEND_ONLY_MATRIX_NAMES
				n_buffer_size_matrix = sizeof(char) * matrix->n * (matrix->max_length_seq_name + 1);
				matrix_seq_data->row = 0L;
		#else
				n_buffer_size_matrix = sizeof(char) * matrix->n * (matrix->m + 1) + sizeof(char) * matrix->n * (matrix->max_length_seq_name + 1);
				/* array for row strings */
				matrix_seq_data->row = (char **) malloc((matrix->n) * sizeof(char *));
		#endif

				/* send names and/or seq and names */
				pack_data = (char *) alloc(n_buffer_size_matrix, "sequences packing");
				MPI_Recv(pack_data, n_buffer_size_matrix, MPI_CHAR, MPI_MAIN_PROCESS, MPI_TAG_NAME_AND_SEQ_DATA, MPI_COMM_WORLD, &status);

				/* array for row title strings */
				position = 0;
				for ( i = 0; i < matrix->n; i++){
		#ifndef MPI_SEND_ONLY_MATRIX_NAMES
					matrix_seq_data->row[i] = (char*) malloc(sizeof(char) * (matrix->m + 1));
					memcpy((char*) (matrix_seq_data->row[i]), pack_data + position, matrix->m + 1);
					position += matrix->m + 1;
		#endif
					matrix_seq_data->rowtitle[i] = (char*) malloc(sizeof(char) * (matrix->max_length_seq_name + 1));
					memcpy((char*) (matrix_seq_data->rowtitle[i]), pack_data + position, matrix->max_length_seq_name + 1);
					position += matrix->max_length_seq_name + 1;
				}
				/* print data.... */
		//		print_data_seq(matrix_seq_data, matrix->n, myMPIid);
				/* END send names and/or seq and names */

				/* send binary data */
				n_buffer_size_binary = matrix->bytes * matrix->n;
				p_pack_data_binary = (Lvb_bit_length *) alloc(n_buffer_size_binary, "binary packing");
				MPI_Recv(p_pack_data_binary, n_buffer_size_binary, MPI_CHAR, MPI_MAIN_PROCESS, MPI_TAG_BINARY_DATA, MPI_COMM_WORLD, &status);

				enc_mat = (Lvb_bit_length **) alloc((matrix->n) * sizeof(Lvb_bit_length *), "state sets");
				for (i = 0; i < matrix->n; i++){
					enc_mat[i] = (Lvb_bit_length *) alloc(sizeof(char) * matrix->bytes, "binary packing");
					memcpy((Lvb_bit_length *) (enc_mat[i]), p_pack_data_binary + i * matrix->nwords, matrix->bytes);
				}
		//		print_binary_data(enc_mat, matrix->n, matrix->nwords, myMPIid);
				/* END send binary data */

				/* try to read the state file */
				if (rcstruct.n_flag_save_read_states == DO_SAVE_READ_STATES){
					/// Try to test the consistency of state file
					printf("Testing consistency checkpoint file MPIid %d.", myMPIid);
					if (is_state_file_exist_and_consistency(myMPIid) == LVB_TRUE){
						printf(" Checkpoint file OK.\n");
						if (is_process_ended_by_MPIid(myMPIid)){
							rcstruct.n_flag_is_finished_process = CHECK_POINT_PROCESS_FINISHED;
							printf("State file has a flag that indicates the MPIid %d is already end\n", myMPIid);
						}
						/* IMPORTANT:
						 * If the seed is defined randomly in the first call, it will be different in this calling.
						 * So, if the state file exist, the seed will be replaced by the previous one.
						 */
						else if (! is_parameters_are_same_from_state(&rcstruct, myMPIid, LVB_TRUE)){
							printf("State file has different parameters than this ones, Start again MPIid %d\n", myMPIid);
							rcstruct.n_flag_is_possible_read_state_files = CHECK_POINT_NOT_READ_STATE_FILES;
						}
						else{
							rcstruct.n_flag_is_possible_read_state_files = CHECK_POINT_READ_STATE_FILES;
						}
						/*
						 * SO, if the flag rcstruct.n_flag_is_possible_read_state_files = CHECK_POINT_READ_STATE_FILES  AND
						 * 		flag rcstruct.n_flag_save_read_states == DO_SAVE_READ_STATES
						 * 		IS GOING TO READ THE STATE FILES.
						 */
					}
					else{
						printf("\nState file not found or corrupted for MPIid: %d\n", myMPIid);
						rcstruct.n_flag_is_possible_read_state_files = CHECK_POINT_NOT_READ_STATE_FILES; /* start again from beginning, some problem with initial file */
					}
				}
				else{
					/* this value is default but ti is better to stay here */
					rcstruct.n_flag_is_possible_read_state_files = CHECK_POINT_NOT_READ_STATE_FILES;
				}
				/* test if the process is finished */
				if (rcstruct.n_flag_is_finished_process == CHECK_POINT_PROCESS_NOT_FINISHED){

					rinit(rcstruct.seed);
					p_bstack_overall = treestack_new();
					log_progress = LVB_TRUE;
					/* get the length and the tree */
					getsoln(matrix, rcstruct,  myMPIid, log_progress, matrix_seq_data, p_bstack_overall, enc_mat);
					/* "file-local" dynamic heap memory */
					treestack_free(p_bstack_overall);

					/* only print the time end the process finish */
					cleanup();
				}
				else{	/* send message to root node to say that this one is finished */
					/* only print the time end the process finish */
					int n_finish_message = MPI_FINISHED;
					MPI_Request request_handle_send = 0;
					MPI_Isend(&n_finish_message, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, &request_handle_send);
					MPI_Wait(&request_handle_send, MPI_STATUS_IGNORE); /* need to do this because the receiver is asynchronous */
					cleanup();
				}

				/* clean memory */
				for (i = 0; i < matrix->n; i++) free(enc_mat[i]);
				free(enc_mat);
				rowfree(matrix_seq_data, matrix->n);
				free(pack_data);
				free(p_pack_data_binary);
				free(p_bstack_overall);
			}
		}
		MPI_Type_free(&mpi_matrix);
		MPI_Type_free(&mpi_params);
		free(matrix_seq_data);
	    free(matrix);
	    if (rcstruct.verbose == LVB_TRUE) {
	        printf("Process finish: %d\n", myMPIid);
	    }

	    /* shut down MPI */
	    MPI_Finalize();
	    return n_error_code;

	} /* end main() */

#endif