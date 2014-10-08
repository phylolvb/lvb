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

#include "ReadFile.h"

void read_file(char *file_name, int n_file_type, Dataptr p_lvbmat){

	CReadFiles readFiles = CReadFiles();
	/// read file
	std::string sz_file_name = std::string(file_name);
	readFiles.read_file(file_name, n_file_type);

    /* scalars */
	//cout << (long) readFiles.get_length_sequences() << endl;
    p_lvbmat->m = (long) readFiles.get_length_sequences();
    p_lvbmat->original_m = p_lvbmat->m;
    p_lvbmat->n = (long) readFiles.get_number_seqs();
    p_lvbmat->nbranches = brcnt(p_lvbmat->n); 		/* possible number of braches */

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



void print_formats_available(){

	printf("\nFormats available to read: phylip, fasta, nexus, msf and clustal");

}

/* "dcbvsi:o:f:p" */
void usage(char *p_file_name){
	printf("Usage: lvb [dcbvsiofp]\n");
	printf("lvb seeks parsimonious trees from an aligned nucleotide data matrix.\n"
			"It uses heuristic searches consisting of simulated annealing followed by hill-climbing.\n\n");

	printf("\n    -b (0) bootstrap replicates, as an integer in the range 1 to %ld inclusive.", (long) MAX_BOOTSTRAPS);
	printf("\n       default (0)");
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
	printf("    -v [t|f] (f) verbose.\n");
	printf("    -f [phylip|fasta|nexus|msf|clustal] (phylip) file format of input file.\n"
			"       default: phylip format\n");
	printf("    -p (1) Threads available."
			"       default: only one thread available\n");
	printf("    -h print this help.");
	printf("    -? print this help.");
	abort();
}


void read_parameters(Params *prms, int argc, char **argv){

	int c;
	opterr = 0;

	while ((c = getopt (argc, argv, "c:b:vs:i:o:f:p:")) != -1){
		switch (c)
		{
			case 'c':	/* cooling schedule */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -c [g|l]\n", optopt);
					usage(argv[0]);
					abort();
				}
				if (strcmp(optarg, "g") == 0 || strcmp(optarg, "G") == 0) prms->cooling_schedule = 0;
				else if (strcmp(optarg, "l") == 0 || strcmp(optarg, "L") == 0) prms->cooling_schedule = 1;
				else{
					fprintf (stderr, "Unknown cooling schedule option\nPlease, choose between Geometric (g) or Linear (l).");
					abort();
				}
				break;
			case 'b':	/* bootstrap */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -b <int>\n", optopt);
					usage(argv[0]);
					abort();
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
					abort();
				}
				prms->seed = atoi(optarg);
				break;
			case 'i':	/* file name in */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -i <file name>\n", optopt);
					usage(argv[0]);
					abort();
				}
				if (strlen(optarg) > LVB_FNAMSIZE){
					fprintf (stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
					abort();
				}
				strcpy(prms->file_name_in, optarg);
				break;
			case 'o':	/* file name out */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -o <file name>\n", optopt);
					usage(argv[0]);
					abort();
				}
				if (strlen(optarg) > LVB_FNAMSIZE){
					fprintf (stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
					abort();
				}
				strcpy(prms->file_name_out, optarg);
				break;
			case 'f':	/* format */
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -f [phylip|fasta|nexus|msf|clustal]\n", optopt);
					usage(argv[0]);
					abort();
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
					abort();
				}
				break;
			case 'p':
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -p <file name>\n", optopt);
					usage(argv[0]);
					abort();
				}
				prms->n_processors_available = atoi(optarg);
				break;
			case '?':
			case 'h':
			default:
	            usage(argv[0]);
		}
	}
}



