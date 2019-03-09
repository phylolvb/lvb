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

/* ========== ReadFile.cpp - command-line inputs ========== */

#include "ReadFile.h"

//#ifdef NP_Implementation
void read_file(char *file_name, int n_file_type, Dataptr p_lvbmat){
//#endif

#ifdef MPI_Implementation
int read_file(char *file_name, int n_file_type, Dataptr p_lvbmat, DataSeqPtr p_lvbmat_seq){
#endif

	CReadFiles readFiles = CReadFiles();

	#ifdef MPI_Implemenation
	int n_error_code
	#endif

	/// read file
	std::string sz_file_name = std::string(file_name);
	
	//ifdef NP_Implementation
	readFiles.read_file(file_name, n_file_type);
	//#endif

	#ifdef MPI_Implementation
	n_error_code = readFiles.read_file(file_name, n_file_type);
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
	#endif

    /* scalars */
	//cout << (long) readFiles.get_length_sequences() << endl;
    p_lvbmat->m = (long) readFiles.get_length_sequences();
    p_lvbmat->original_m = p_lvbmat->m;
    p_lvbmat->n = (long) readFiles.get_number_seqs();
    p_lvbmat->nbranches = brcnt(p_lvbmat->n); 		/* possible number of braches */

	#ifdef MPI_Implementation
	p_lvbmat->max_length_seq_name = readFiles.get_max_length_seq_name();
	#endif

    /* it is used in tree compare */
    p_lvbmat->nsets = p_lvbmat->n - 3;				/* sets per tree */
    p_lvbmat->mssz = p_lvbmat->n - 2;				/* maximum objects per set */

	//#ifdef NP_Implementation
    /* array for row title strings */
    p_lvbmat->rowtitle = (char **) malloc((p_lvbmat->n) * sizeof(char *));
    if (p_lvbmat->rowtitle == NULL) readFiles.exit_error(1 , "Fail to allocate memory...");

    /* array for row strings */
    p_lvbmat->row = (char **) malloc((p_lvbmat->n) * sizeof(char *));
    if (p_lvbmat->row == NULL) readFiles.exit_error(1 , "Fail to allocate memory...");

    /* we want null-terminated strings, so we cannot simply point to
     * the same, non-null-terminated arrays as are found in PHYLIP's
     * data structures */
    for (int i = 0; i < p_lvbmat->n; i++){
    	p_lvbmat->rowtitle[i] = (char*) malloc(sizeof(char) * (readFiles.get_max_length_seq_name() + 1));
    	p_lvbmat->row[i] = (char*) malloc(sizeof(char) * (p_lvbmat->m + 1));
    }
    for (int i = 0; i < p_lvbmat->n; i++) {
        for (int j = 0; j < p_lvbmat->m; j++) p_lvbmat->row[i][j] = readFiles.get_char_sequences(i, j);
        p_lvbmat->row[i][p_lvbmat->m] = '\0';
    }
    for (int i = 0; i < p_lvbmat->n; i++) {
        for (int j = 0; j < readFiles.get_length_seq_name(i); j++) p_lvbmat->rowtitle[i][j] = readFiles.get_char_seq_name(i, j);
        p_lvbmat->rowtitle[i][readFiles.get_length_seq_name(i)] = '\0';
    }
/*	std::string file_name_out = "/home/mmp/Downloads/file_nexus_nex_out.fas";
	readFiles.save_file(file_name_out);*/
}

	void free_lvbmat_structure(DataStructure *p_lvbmat){

		if (p_lvbmat->row != NULL){
			for(int i = 0; i < p_lvbmat->n; ++i) free(p_lvbmat->row[i]);
			free(p_lvbmat->row);
			p_lvbmat->row = NULL;
		}
		if (p_lvbmat->rowtitle != NULL){
			for(int i = 0; i < p_lvbmat->n; ++i) free(p_lvbmat->rowtitle[i]);
			free(p_lvbmat->rowtitle);
			p_lvbmat->rowtitle = NULL;
		}
	//	free(p_lvbmat);
	//	p_lvbmat = NULL;
}
	//#endif

	#ifdef MPI_Implementation

	/* array for row title strings */
		p_lvbmat_seq->rowtitle = (char **) malloc((p_lvbmat->n) * sizeof(char *));
		if (p_lvbmat_seq->rowtitle == NULL) return readFiles.exit_error(1 , "Fail to allocate memory...");

		/* array for row strings */
		p_lvbmat_seq->row = (char **) malloc((p_lvbmat->n) * sizeof(char *));
		if (p_lvbmat_seq->row == NULL) return readFiles.exit_error(1 , "Fail to allocate memory...");

		/* we want null-terminated strings, so we cannot simply point to
		 * the same, non-null-terminated arrays as are found in PHYLIP's
		 * data structures */
		for (int i = 0; i < p_lvbmat->n; i++){
			p_lvbmat_seq->rowtitle[i] = (char*) malloc(sizeof(char) * (readFiles.get_max_length_seq_name() + 1));
			p_lvbmat_seq->row[i] = (char*) malloc(sizeof(char) * (p_lvbmat->m + 1));
		}
		for (int i = 0; i < p_lvbmat->n; i++) {
			for (int j = 0; j < p_lvbmat->m; j++) p_lvbmat_seq->row[i][j] = readFiles.get_char_sequences(i, j);
			p_lvbmat_seq->row[i][p_lvbmat->m] = '\0';
		}
		for (int i = 0; i < p_lvbmat->n; i++) {
			for (int j = 0; j < readFiles.get_length_seq_name(i); j++) p_lvbmat_seq->rowtitle[i][j] = readFiles.get_char_seq_name(i, j);
			p_lvbmat_seq->rowtitle[i][readFiles.get_length_seq_name(i)] = '\0';
		}
	/*	std::string file_name_out = "/home/mmp/Downloads/file_nexus_nex_out.fas";
		readFiles.save_file(file_name_out);*/
		return EXIT_SUCCESS;
	}

	void free_lvbmat_structure(DataSeqStructure *p_lvbmat_seq, int n_size){

		if (p_lvbmat_seq->row != NULL){
			for(int i = 0; i < n_size; ++i) free(p_lvbmat_seq->row[i]);
			free(p_lvbmat_seq->row);
			p_lvbmat_seq->row = NULL;
		}
		if (p_lvbmat_seq->rowtitle != NULL){
			for(int i = 0; i < n_size; ++i) free(p_lvbmat_seq->rowtitle[i]);
			free(p_lvbmat_seq->rowtitle);
			p_lvbmat_seq->rowtitle = NULL;
		}
	//	free(p_lvbmat);
	//	p_lvbmat = NULL;
	}

long brcnt(long n){
	return (n << 1) - 3;
}; /* return number of branches in unrooted binary tree structure containing n tips */

	#endif

void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name){

	CReadFiles readFiles = CReadFiles();
	/// read file
	std::string sz_file_name = std::string(file_name);
	readFiles.read_file(file_name, n_file_type);

	*species_ptr = (long) readFiles.get_number_seqs();
	*sites_ptr = (long) readFiles.get_length_sequences();
	*max_length_name = readFiles.get_max_length_seq_name();
//	std::string file_name_out = "/home/mmp/workspace/LVB_READ_FILES/aln_files/file_nexus_nex_out.fas";
//	readFiles.save_file(file_name_out);

}

void print_formats_available(){

	printf("\nFormats available to read: phylip, fasta, nexus, and clustal");

}

/* "dcbvsi:o:f:p" */
//#ifdef NP_Implementation 
void usage(char *p_file_name){
	printf("lvb seeks parsimonious trees from an aligned nucleotide data matrix.\n"
			"It uses heuristic searches consisting of simulated annealing followed\n"
			"by hill-climbing.\n\n");

	printf("Usage: lvb [dcbvsioftp]\n");
	printf("\n    -b (0) bootstrap replicates, as an integer in the range 1 to %ld"
	"\n       inclusive.", (long) MAX_BOOTSTRAPS);
	printf("\n       Default (0).");
	printf("\n    -a algorithm zero (0) NNI + SPR,\n"
	        "       algorithm one (1) NNI + SPR + TBR.\n"
			"		algorithm two (2) Point based calculations.\n");
	printf("\n    -c [g|l] (g) cooling schedule. The schedule chosen will\n"
			"       affect the quality and speed of the simulated annealing search.\n"
			"       The GEOMETRIC (g) schedule will take significantly less time,\n"
			"       but may produce lower quality results. The LINEAR (l) schedule may\n"
			"       produce higher quality results, at the cost of increased runtime.\n"
			"       Default (g), the GEOMETRIC schedule.\n");
	printf("    -i input file name.\n"
			"       Default: 'infile'.\n");
	printf("    -o output file name.\n"
			"       Default: 'outfile'.\n");
	printf("    -s specify a random number seed, or use default.\n"
			"       Default: it is taken from the system clock.\n");
	printf("    -v turn verbose mode on.\n");
	printf("    -f [clustal|fasta|nexus|phylip].\n"
			"       Default is 'PHYLIP' format\n");
	printf("    -p (1) Number of threads to be used.\n"
			"       Default: one thread.\n");
	printf("    -t (0) Max number of trees to be saved. Once reached, LVB will stop."
				"          Default: all of them will be saved.\n");
	printf("    -h print help.\n");
	printf("    -? print help.\n");
	exit(0);
}
//#endif

#ifdef MPI_Implementation

void usage(char *p_file_name){
#ifdef MAP_Reduce
	printf("Usage: lvb [cvsSCiofph]\n");
#else
	printf("Usage: lvb [cvsSCNiofph]\n");
#endif

	printf("lvb seeks parsimonious trees from an aligned nucleotide data matrix.\n"
			"It uses heuristic searches consisting of simulated annealing followed by hill-climbing.\n\n");

	printf("\n    -c [g|l] (g) cooling schedule. The schedule chosen\n"
			"       will affect the quality and speed of the simulated annealing\n"
			"       search. The GEOMETRIC (g) schedule will take significantly less time,\n"
			"       but may produce lower quality results. The LINEAR (l) schedule may\n"
			"       produce higher quality results, at the cost of increased runtime.\n"
			"       default (g), is the GEOMETRIC schedule.\n");
	printf("    -i input file name.\n"
			"       default: 'infile'\n");
	printf("    -o output file name.\n"
			"       default: 'outfile'\n");
	printf("    -s specify a random number seed, or use default.\n"
			"       default: it is taken from the system clock.\n");
#ifndef MAP_REDUCE_SINGLE
	printf("    -N specify the number of seeds to try, need to be greater than number of mpi process.\n"
			"       default: it is the number of mpi process.\n");
#endif
	printf("    -C value in seconds when a checkpoint file is going to be saved.\n"
					"       default: 1800\n");
	printf("    -S if defined it is going to read and save states at a specific time points.\n"
			"       If you want to read a state the command line need to be exactly the previous\n        on save state point plus -S flag.\n"
			"       The state files need to be in the same directory where the lvb_mpi is called.\n"
			"       The states files has the names \"state_<number of mpi process>.dat>\"\n"
			"       default: it is not going to save/read states.\n");
	printf("    -v [t|f] (f) verbose.\n");
	printf("    -f [phylip|fasta|nexus|msf|clustal] (phylip) file format of input file.\n"
			"       default: phylip format\n");
	printf("    -p (1) Threads available."
			"       default: only one thread available\n");
	printf("    -h print this help.\n");
	printf("    -? print this help.\n");
}

#endif

//#ifdef NP_Implementation
void read_parameters(Params *prms, int argc, char **argv){

	int c;
	opterr = 0;

	while ((c = getopt (argc, argv, "t:c:b:vs:i:o:f:p:a:")) != -1){
		switch (c)
		{
			case 'c':	/* cooling schedule */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -c [g|l]\n", optopt);
					usage(argv[0]);
					exit(1);
				}
				if (strcmp(optarg, "g") == 0 || strcmp(optarg, "G") == 0) prms->cooling_schedule = 0;
				else if (strcmp(optarg, "l") == 0 || strcmp(optarg, "L") == 0) prms->cooling_schedule = 1;
				else{
					fprintf (stderr, "Unknown cooling schedule option\nPlease, choose between Geometric (g) or Linear (l).");
					exit(1);
				}
				break;
			case 'a':	/* algorithm selection */
				if (optarg == NULL){
					fprintf (stderr, "Option -%a requires an argument -a [0|1]\n", optopt);
					usage(argv[0]);
					exit(1);
				}
				if (strcmp(optarg, "0") == 0 || strcmp(optarg, "0") == 0) prms->algorithm_selection = 0;
				else if (strcmp(optarg, "1") == 0 || strcmp(optarg, "1") == 0) prms->algorithm_selection = 1;
				else if (strcmp(optarg, "2") == 0 || strcmp(optarg, "2") == 0) prms->algorithm_selection = 2;
				else{
					fprintf (stderr, "Unknown algorithm option\nPlease, choose between SN (0) or SEQ-TNS (1).");
					exit(1);
				}
				break;
			case 'b':	/* bootstrap */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -b <int>\n", optopt);
					usage(argv[0]);
				}
				prms->bootstraps = atoi(optarg);
				break;
			case 'v':	/* verbose */
				prms->verbose = LVB_TRUE;
				break;
			case 's':	/* seed */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -s <int>\n", optopt);
					usage(argv[0]);
				}
				prms->seed = atoi(optarg);
				break;
			case 'i':	/* file name in */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -i <file name>\n", optopt);
					usage(argv[0]);
				}
				if (strlen(optarg) > LVB_FNAMSIZE){
					fprintf (stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
					exit(1);
				}
				strcpy(prms->file_name_in, optarg);
				break;
			case 'o':	/* file name out */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -o <file name>\n", optopt);
					usage(argv[0]);
				}
				if (strlen(optarg) > LVB_FNAMSIZE){
					fprintf (stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
					exit(1);
				}
				strcpy(prms->file_name_out, optarg);
				break;
			case 'f':	/* format */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -f [phylip|fasta|nexus|clustal]\n", optopt);
					usage(argv[0]);
				}
				if (strcmp(optarg, "phylip") == 0){
					prms->n_file_format = FORMAT_PHYLIP;
				}
				else if (strcmp(optarg, "fasta") == 0){
					prms->n_file_format = FORMAT_FASTA;
				}
				else if (strcmp(optarg, "nexus") == 0){
					prms->n_file_format = FORMAT_NEXUS;
				}
				else if (strcmp(optarg, "clustal") == 0){
					prms->n_file_format = FORMAT_CLUSTAL;
				}
				else{
					fprintf (stderr, "Unknown file format.");
					print_formats_available();
					exit(1);
				}
				break;
			case 'p':
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -p <file name>\n", optopt);
					usage(argv[0]);
					exit(1);
				}
				prms->n_processors_available = atoi(optarg);
				if (prms->n_processors_available < 1) prms->n_processors_available = 1;
				break;
			case 't':
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -p <file name>\n", optopt);
					usage(argv[0]);
					exit(1);
				}
				prms->n_number_max_trees = atoi(optarg);
				if (prms->n_number_max_trees < 1) prms->n_number_max_trees = 0;
				break;
			case '?':
			case 'h':
			default:
	            usage(argv[0]);
		}
	}
}

//#endif

#ifdef MPI_Implementation

int read_parameters(Params *prms, int argc, char **argv){

	int c;
	opterr = 0;

	while ((c = getopt (argc, argv, "t:c:vs:i:o:f:p:N:SC:h?")) != -1){
		switch (c)
		{
			case 'c':	/* cooling schedule */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -c [g|l]\n", optopt);
					usage(argv[0]);
					return 1;
				}
				if (strcmp(optarg, "g") == 0 || strcmp(optarg, "G") == 0) prms->cooling_schedule = 0;
				else if (strcmp(optarg, "l") == 0 || strcmp(optarg, "L") == 0) prms->cooling_schedule = 1;
				else{
					fprintf (stderr, "Unknown cooling schedule option\nPlease, choose between Geometric (g) or Linear (l).");
					return 1;
				}
				break;
			case 'v':	/* verbose */
				prms->verbose = LVB_TRUE;
				break;
			case 's':	/* seed */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -s <int>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				prms->seed = atoi(optarg);
				break;
			case 'S':	/* number of seeds to try */
				prms->n_flag_save_read_states = DO_SAVE_READ_STATES;
				break;
#ifndef MAP_REDUCE_SINGLE
			case 'N':	/* number of seeds to try */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -N <int>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				prms->n_seeds_need_to_try = atoi(optarg);
				break;
			case 'C':	/* checkpoint in seconds */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -C <int>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				prms->n_checkpoint_interval = atoi(optarg);
				break;
			/* this options is not visible to the regular users, only for making tests */
			case 't':	/* only used for tests */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -t <int>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				prms->n_make_test = atoi(optarg);
				break;
#endif
			case 'i':	/* file name in */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -i <file name>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				if (strlen(optarg) > LVB_FNAMSIZE){
					fprintf (stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
					return 1;
				}
				strcpy(prms->file_name_in, optarg);
				break;
			case 'o':	/* file name out */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -o <file name>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				if (strlen(optarg) > LVB_FNAMSIZE){
					fprintf (stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
					return 1;
				}
				strcpy(prms->file_name_out, optarg);
				break;
			case 'f':	/* format */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -f [phylip|fasta|nexus|msf|clustal]\n", optopt);
					usage(argv[0]);
					return 1;
				}
				if (strcmp(optarg, "phylip") == 0){
					prms->n_file_format = FORMAT_PHYLIP;
				}
				else if (strcmp(optarg, "fasta") == 0){
					prms->n_file_format = FORMAT_FASTA;
				}
				else if (strcmp(optarg, "nexus") == 0){
					prms->n_file_format = FORMAT_NEXUS;
				}
				else if (strcmp(optarg, "msf") == 0){
					prms->n_file_format = FORMAT_MSF;
				}
				else if (strcmp(optarg, "clustal") == 0){
					prms->n_file_format = FORMAT_CLUSTAL;
				}
				else{
					fprintf (stderr, "Unknown file format.");
					print_formats_available();
					return 1;
				}
				break;
			case 'p':
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -p <file name>\n", optopt);
					usage(argv[0]);
					return 1;
				}
				prms->n_processors_available = atoi(optarg);
				if (prms->n_processors_available < 1) prms->n_processors_available = 1;
				break;
			case '?':
			case 'h':
			default:
	            usage(argv[0]);
	            return 1;
		}
	}
	return EXIT_SUCCESS;
}

#endif
