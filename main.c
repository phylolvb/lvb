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


/* ********** main.c - LVB ********** */

#include "lvb.h"
#include <inttypes.h>

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

static void writeinf(Params prms, int myMPIid, int n_process)
/* write initial details to standard output */
{
    printf("\n#########\nProcess ID: %d\n", myMPIid);

    printf("cooling schedule     = ");
    if(prms.cooling_schedule == 0) printf("GEOMETRIC\n");
    else printf("LINEAR\n");

    printf("main seed            = %d\n", prms.seed);
    printf("bootstrap replicates = %ld\n", prms.bootstraps);
    printf("input file name      = %s\n", prms.file_name_in);
    printf("output file name     = %s\n", prms.file_name_out);
    if (prms.n_file_format == FORMAT_PHYLIP) printf("format input file    = phylip\n");
    else if (prms.n_file_format == FORMAT_FASTA) printf("format input file    = fasta\n");
    else if (prms.n_file_format == FORMAT_NEXUS) printf("format input file    = nexus\n");
    else if (prms.n_file_format == FORMAT_MSF) printf("format input file    = msf\n");
    else if (prms.n_file_format == FORMAT_CLUSTAL) printf("format input file    = clustal\n");
    else{
    	fprintf (stderr, "Error, input format file not recognized\n");
    	abort();
    }
    printf("bootstrap replicates = %ld\n", prms.bootstraps);
    printf("threads              = %d\n", prms.n_processors_available);
    printf("mpi process          = %d\n", n_process);
} /* end writeinf() */


static void logtree1(Dataptr matrix, DataSeqPtr restrict matrix_seq_data, const Branch *const barray, const long start, const long cycle, long root)
/* log initial tree for cycle cycle of start start (in barray) to outfp */
{
    static char outfnam[LVB_FNAMSIZE]; 	/* current file name */
    int fnamlen;			/* length of current file name */
    FILE *outfp;			/* output file */

    fnamlen = sprintf(outfnam, "%s_start%ld_cycle%ld", TREE1FNAM, start, cycle);
    lvb_assert(fnamlen < LVB_FNAMSIZE);	/* shut door if horse bolted */

    /* create tree file */
    outfp = clnopen(outfnam, "w");
    lvb_treeprint(matrix, matrix_seq_data, outfp, barray, root);
    clnclose(outfp, outfnam);

} /* end logtree1() */

static long getsoln(Dataptr restrict matrix, DataSeqPtr restrict matrix_seq_data, Params rcstruct, const long *weight_arr, long *iter_p,
		Treestack* bstack_overall, Lvb_bit_lentgh **enc_mat, int myMPIid, Lvb_bool log_progress)
/* get and output solution(s) according to parameters in rcstruct;
 * return length of shortest tree(s) found, using weights in weight_arr */
{
    static char fnam[LVB_FNAMSIZE];	/* current file name */
    long fnamlen;			/* length of current file name */
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
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */

    /* NOTE: These variables and their values are "dummies" and are no longer
     * used in the current version of LVB. However, in order to keep the
     * formatting of the output compatible with that of previous versions of
     * LVB these variables will continue to be used and written to the summary
     * files.  */
    long cyc = 0;	/* current cycle number */
    long start = 0;	/* current random (re)start number */
 
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
	
    /* determine starting temperature */
    randtree(matrix, tree);	/* initialise required variables */
    ss_init(matrix, tree, enc_mat);
    initroot = 0;
    t0 = get_initial_t(matrix, tree, rcstruct, initroot, weight_arr, myMPIid, log_progress);
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
		fprintf(sumfp, "%d\t%ld\t%ld\t%ld\t", myMPIid, start, cyc, getplen(matrix, tree, rcstruct, initroot, weight_arr, p_todo_arr, p_todo_arr_sum_changes, p_runs));
		free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
		logtree1(matrix, matrix_seq_data, tree, start, cyc, initroot);
    }

    /* find solution(s) */
    treelength = anneal(matrix, bstack_overall, tree, rcstruct, initroot, t0, maxaccept,
    		maxpropose, maxfail, stdout, weight_arr, iter_p, myMPIid, log_progress);
    treestack_pop(matrix, tree, &initroot, bstack_overall);
    treestack_push(matrix, bstack_overall, tree, initroot);
    treelength = deterministic_hillclimb(matrix, bstack_overall, tree, rcstruct, initroot, stdout,
    		weight_arr, iter_p, myMPIid, log_progress);

	/* log this cycle's solution and its details 
	 * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */
    if (rcstruct.verbose == LVB_TRUE){
		fnamlen = sprintf(fnam, "%s_start%ld_cycle%ld", RESFNAM, start, cyc);
		lvb_assert(fnamlen < LVB_FNAMSIZE);	/* really too late */
		resfp = clnopen(fnam, "w");
		treec = treestack_print(matrix, matrix_seq_data, bstack_overall, resfp, LVB_FALSE);
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

    /* "local" dynamic heap memory */
    free(tree);
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
		matrix->n_slice_size_getplen = 0;	/* it doens't matter this value for 1 thread */
		matrix->n_threads_getplen = 1; /* need to pass for 1 thread because the number of words is to low */
	}
	printf("\nthreads that will be used  = %d\n", matrix->n_threads_getplen);
	printf("(because is related with the size of the data)\n");
}

int get_other_seed_to_run_a_process(){
	return (int) (rand() % (unsigned long) MAX_SEED);
}


static void logstim(void)
/* log start time with message */
{
    time_t tim;	/* time */

    tim = time(NULL);
    printf("Starting at: %s", ctime(&tim));

} /* end logstim() */

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
#ifdef MPI_SEND_ONLY_MATRIX_NAMES
		printf("%s\n", p_lvbmat->rowtitle[i]);
#else
		printf("%s  %s\n", p_lvbmat->rowtitle[i], p_lvbmat->row[i]);
#endif
	}
}

void print_binary_data(Lvb_bit_lentgh **p_enc_data, int n_size, int nwords, int n_thread){

	printf("########### SEQ ####################\n thread: %d\n", n_thread);
	int i, x;
	for(i = 0; i < n_size; ++i){
		printf("	rank:%d  u: ", i);
		for(x = 0; x < nwords; x++){ printf("%" PRIu64, p_enc_data[i][x]); }
		printf("\n");
	}
}


int main(int argc, char **argv)
{
    Dataptr matrix;	/* data matrix */
    DataSeqPtr matrix_seq_data;
    int val = EXIT_SUCCESS;			/* return value */
    Params rcstruct;		/* configurable parameters */
    long i, n_buffer_size_matrix, n_buffer_size_binary;			/* loop counter */
    int position, n_error_code;
    long iter;			/* iterations of annealing algorithm */
    long replicate_no = 0L;	/* current bootstrap replicate number */
    long trees_output_total = 0L;	/* number of trees output, overall */
    long trees_output;		/* number of trees output for current rep. */
    double total_iter = 0.0;	/* total iterations across all replicates */
    long final_length;		/* length of shortest tree(s) found */
    FILE *outtreefp;		/* best trees found overall */
    long *weight_arr;
    Lvb_bit_lentgh **enc_mat;/* encoded data mat. */
    Lvb_bool log_progress;	/* whether or not to log anneal search */
    Treestack bstack_overall;	/* overall best tree stack */
    char *pack_data;
    Lvb_bit_lentgh *p_pack_data_binary;
    /* global files */

    /* define mpi */
    MPI_Status	status ;		/* return status for receive */
    int num_procs, myMPIid, mpi_err;
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

	/* define the matrix structure MPI */
	int				nItems = 2;
	int          	blocklengths[2] = {10, 2};
	MPI_Datatype 	types[2] = {MPI_LONG, MPI_INT};
	MPI_Datatype 	mpi_matrix;
	MPI_Aint     	displacements[2];
	displacements[0] = offsetof(DataStructure, m);
	displacements[1] = offsetof(DataStructure, n_threads_getplen);
	MPI_Type_create_struct(nItems, blocklengths, displacements, types, &mpi_matrix);
	MPI_Type_commit(&mpi_matrix);


	nItems = 3;
	int          	blocklengthsPar [3] = {9, 4, LVB_FNAMSIZE + LVB_FNAMSIZE};
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
		/* entitle standard output */
		printf("\nLVB\n\n"
		"(c) Copyright 2003-2012 by Daniel Barker\n"
		"(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl\n"
		"(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian\n"
		"Strobl.\n"
		"All rights reserved.\n"
		"\n"
		"Redistribution and use in source and binary forms, with or without\n"
		"modification, are permitted provided that the following conditions\n"
		"are met:\n"
		"\n"
		"1. Redistributions of source code must retain the above copyright\n"
		"notice, this list of conditions and the following disclaimer.\n"
		"\n"
		"2. Redistributions in binary form must reproduce the above\n"
		"copyright notice, this list of conditions and the following\n"
		"disclaimer in the documentation and/or other materials provided\n"
		"with the distribution.\n"
		"\n"
		"3. Neither the name of the copyright holder nor the names of its\n"
		"contributors may be used to endorse or promote products derived\n"
		"from this software without specific prior written permission.\n"
		"\n"
		"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
		"\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
		"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
		"FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
		"COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,\n"
		"INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"
		"(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n"
		"SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)\n"
		"HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,\n"
		"STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)\n"
		"ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF\n"
		"ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n");
		printf("* This is %s version %s %s *\n\n", PROGNAM, LVB_VERSION, LVB_SUBVERSION);
		printf("Literature reference:\n"
		"Barker, D. 2004. LVB: Parsimony and simulated annealing in the\n"
		"search for phylogenetic trees. Bioinformatics, 20, 274-275.\n\n");
		printf("Download and support:\n"
		"http://eggg.st-andrews.ac.uk/lvb\n\n");

		lvb_initialize();
		n_error_code = getparam(&rcstruct, argc, argv);
		if (n_error_code != 0){
			MPI_Abort(MPI_COMM_WORLD, n_error_code);
			exit(n_error_code);
		}
		logstim();

		/* read and alloc space to the data structure */
		n_error_code = phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, matrix, matrix_seq_data);
		if (n_error_code != 0){
			MPI_Abort(MPI_COMM_WORLD, n_error_code);
			exit(n_error_code);
		}
		writeinf(rcstruct, myMPIid, num_procs);
		matchange(matrix, matrix_seq_data, rcstruct);	/* cut columns */
		calc_distribution_processors(matrix, rcstruct);

	    /* Allocation of the initial encoded matrix is non-contiguous because
	     * this matrix isn't used much, so any performance penalty won't matter. */
		enc_mat = (Lvb_bit_lentgh **) alloc((matrix->n) * sizeof(Lvb_bit_lentgh *), "state sets");
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
	    p_pack_data_binary = (Lvb_bit_lentgh *) alloc(n_buffer_size_binary, "binary packing");
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

	}
	else{

		MPI_Recv(matrix, 1, mpi_matrix, MPI_MAIN_PROCESS, MPI_TAG_MATRIX, MPI_COMM_WORLD, &status);
//		print_data(matrix, myMPIid); /* print data.... */

		MPI_Recv(&rcstruct, 1, mpi_params, MPI_MAIN_PROCESS, MPI_TAG_PARAMS, MPI_COMM_WORLD, &status);
		sprintf(rcstruct.file_name_out, "%s_%d", rcstruct.file_name_out, myMPIid);
//		writeinf(rcstruct, myMPIid, num_procs); /* print data.... */

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
		p_pack_data_binary = (Lvb_bit_lentgh *) alloc(n_buffer_size_binary, "binary packing");
		MPI_Recv(p_pack_data_binary, n_buffer_size_binary, MPI_CHAR, MPI_MAIN_PROCESS, MPI_TAG_BINARY_DATA, MPI_COMM_WORLD, &status);

		enc_mat = (Lvb_bit_lentgh **) alloc((matrix->n) * sizeof(Lvb_bit_lentgh *), "state sets");
		for (i = 0; i < matrix->n; i++){
			enc_mat[i] = (Lvb_bit_lentgh *) alloc(sizeof(char) * matrix->bytes, "binary packing");
			memcpy((Lvb_bit_lentgh *) (enc_mat[i]), p_pack_data_binary + i * matrix->nwords, matrix->bytes);
		}
//		print_binary_data(enc_mat, matrix->n, matrix->nwords, myMPIid);
		/* END send binary data */

		rinit(rcstruct.seed);

		/* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
		bstack_overall = treestack_new();

		if (rcstruct.bootstraps > 0) {
			log_progress = LVB_FALSE;
			printf("\nProcess: %d    Replicate:      Rearrangements: Trees output:   Length:\n", myMPIid);
		}
		else log_progress = LVB_TRUE;

		weight_arr = (long*) alloc(sizeof(long) * matrix->m, "alloc data structure");
		outtreefp = clnopen(rcstruct.file_name_out, "w");
		do{
			iter = 0;
			if (rcstruct.bootstraps > 0){
				get_bootstrap_weights(weight_arr, matrix->m, matrix->original_m - matrix->m);
			}
			else{
				for (i = 0; i < matrix->m; i++) weight_arr[i] = 1;
			}
			final_length = getsoln(matrix, matrix_seq_data, rcstruct, weight_arr, &iter, &bstack_overall, enc_mat, myMPIid, log_progress);
			if (rcstruct.bootstraps > 0) trees_output = treestack_print(matrix, matrix_seq_data, &bstack_overall, outtreefp, LVB_TRUE);
			else trees_output = treestack_print(matrix, matrix_seq_data, &bstack_overall, outtreefp, LVB_FALSE);

			trees_output_total += trees_output;
			treestack_clear(&bstack_overall);
			replicate_no++;
			if (rcstruct.bootstraps > 0) {
				printf("%-16ld%-16ld%-16ld%ld\n", replicate_no, iter, trees_output, final_length);
				total_iter += (double) iter;
			}
			else  printf("\nRearrangements tried: %ld\n", iter);
		} while (replicate_no < rcstruct.bootstraps);

		clnclose(outtreefp, rcstruct.file_name_out);

		printf("\n");
		if (rcstruct.bootstraps > 0)
			printf("Total rearrangements tried across all replicates: %g\n\n", total_iter);

		if ((trees_output_total == 1L) && (rcstruct.bootstraps == 0)) {
			printf("1 most parsimonious tree of length %ld written to file '%s'\n", final_length, rcstruct.file_name_out);
		}
		else {
			if (rcstruct.bootstraps > 0){
				lvb_assert(trees_output_total == rcstruct.bootstraps);
				printf("%ld trees written to file '%s'\n", trees_output_total,
				rcstruct.file_name_out);
			}
			else{
				printf("%ld equally parsimonious trees of length %ld written to "
						"file '%s'\n", trees_output_total, final_length, rcstruct.file_name_out);
			}
		}

		/* "file-local" dynamic heap memory */
		treestack_free(&bstack_overall);
		free(weight_arr);

		/* only print the time */
		cleanup();
	}
	MPI_Type_free(&mpi_matrix);
	MPI_Type_free(&mpi_params);

	for (i = 0; i < matrix->n; i++) free(enc_mat[i]);
	free(enc_mat);
	rowfree(matrix_seq_data, matrix->n);
	free(matrix_seq_data);
	free(pack_data);
	free(p_pack_data_binary);
    free(matrix);

    if (rcstruct.verbose == LVB_TRUE) {
        printf("Process finish: %d\n", myMPIid);
    }
 //   printf("Process finish: %d\n", myMPIid);

    /* shut down MPI */
    MPI_Finalize();
    return val;

} /* end main() */
