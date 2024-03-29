/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== CommandLineParser.cpp - command-line inputs ========== */

#include "CommandLineParser.h"

int read_file(char *file_name, int n_file_type, Dataptr p_lvbmat)
{

	CReadFiles readFiles = CReadFiles();
	int n_error_code;

	/// read file
	std::string sz_file_name = std::string(file_name);
	n_error_code = readFiles.read_file(file_name, n_file_type);
	if (n_error_code != EXIT_SUCCESS)
		return n_error_code;

	/* scalars */
	// cout << (long) readFiles.get_length_sequences() << endl;
	p_lvbmat->m = (long)readFiles.get_length_sequences();
	p_lvbmat->original_m = p_lvbmat->m;
	p_lvbmat->n = (long)readFiles.get_number_seqs();
	p_lvbmat->numberofpossiblebranches = brcnt(p_lvbmat->n); /* possible number of braches */
	p_lvbmat->max_length_seq_name = readFiles.get_max_length_seq_name();

	/* it is used in tree compare */
	p_lvbmat->nsets = p_lvbmat->n - 3; /* sets per tree */
	p_lvbmat->mssz = p_lvbmat->n - 2;  /* maximum objects per set */

	/* array for row title strings */
	p_lvbmat->rowtitle = (char **)malloc((p_lvbmat->n) * sizeof(char *));
	if (p_lvbmat->rowtitle == NULL)
		return readFiles.exit_error(1, "Fail to allocate memory...");

	/* array for row strings */
	p_lvbmat->row = (char **)malloc((p_lvbmat->n) * sizeof(char *));
	if (p_lvbmat->row == NULL)
		return readFiles.exit_error(1, "Fail to allocate memory...");

	/* we want null-terminated strings, so we cannot simply point to
	 * the same, non-null-terminated arrays as are found in PHYLIP's
	 * data structures */
	for (int i = 0; i < p_lvbmat->n; i++)
	{
		p_lvbmat->rowtitle[i] = (char *)malloc(sizeof(char) * (readFiles.get_max_length_seq_name() + 1));
		p_lvbmat->row[i] = (char *)malloc(sizeof(char) * (p_lvbmat->m + 1));
	}
	for (int i = 0; i < p_lvbmat->n; i++)
	{
		for (int j = 0; j < p_lvbmat->m; j++)
			p_lvbmat->row[i][j] = readFiles.get_char_sequences(i, j);
		p_lvbmat->row[i][p_lvbmat->m] = '\0';
	}
	for (int i = 0; i < p_lvbmat->n; i++)
	{
		for (int j = 0; j < readFiles.get_length_seq_name(i); j++)
			p_lvbmat->rowtitle[i][j] = readFiles.get_char_seq_name(i, j);
		p_lvbmat->rowtitle[i][readFiles.get_length_seq_name(i)] = '\0';
	}
	/*	std::string file_name_out = "/home/mmp/Downloads/file_nexus_nex_out.fas";
		readFiles.save_file(file_name_out);*/
	return EXIT_SUCCESS;
}

void free_lvbmat_structure(DataStructure *p_lvbmat, int n_size)
{

	if (p_lvbmat->row != NULL)
	{
		for (int i = 0; i < n_size; ++i)
			free(p_lvbmat->row[i]);
		free(p_lvbmat->row);
		p_lvbmat->row = NULL;
	}
	if (p_lvbmat->rowtitle != NULL)
	{
		for (int i = 0; i < n_size; ++i)
			free(p_lvbmat->rowtitle[i]);
		free(p_lvbmat->rowtitle);
		p_lvbmat->rowtitle = NULL;
	}
	//	free(p_lvbmat);
	//	p_lvbmat = NULL;
}

long brcnt(long n)
{
	return (n << 1) - 3;
}; /* return number of branches in unrooted binary tree structure containing n tips */

void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name)
{

	CReadFiles readFiles = CReadFiles();
	/// read file
	std::string sz_file_name = std::string(file_name);
	readFiles.read_file(file_name, n_file_type);

	*species_ptr = (long)readFiles.get_number_seqs();
	*sites_ptr = (long)readFiles.get_length_sequences();
	*max_length_name = readFiles.get_max_length_seq_name();
	//	std::string file_name_out = "/home/mmp/workspace/LVB_READ_FILES/aln_files/file_nexus_nex_out.fas";
	//	readFiles.save_file(file_name_out);
}

void print_formats_available()
{

	printf("\nFormats available to read: phylip, fasta, nexus, msf and clustal");
}

/* "dcbvsi:o:f:p" */
void usage(char *p_file_name)
{
	printf(" Usage: lvb -i <alignment> [options]\n");

	printf("\n");

	printf(" Help options: \n");
	printf("    Display program instructions      -h \n");
	printf("    Display program instructions      -? \n");
	// printf("  Display program information     -# \n");

	printf("\n");

	printf(" General Preferences: \n");
	printf("   Verbose mode                       -v [ON|OFF]         Engage verbose mode; default: OFF\n");

	printf("\n");

	printf(" Input/Output Preferences: \n");
	printf("    Input file                        -i [FILE]           Input file name \n");
	printf("    Input file format                 -f [STRING]         'phylip'|'fasta'|'nexus'|'clustal' \n");
	printf("    Output file                       -o [FILE]           Output file name; default: 'outfile'\n");

	printf("\n");

	printf(" Search Origin Preferences: \n");
	printf("    Starting seed                     -s [VALUE]          Specify a starting seed; default: generated from system clock\n");

	printf("\n");

	printf(" Topological Rearrangement Preferences: \n");
	printf("     Rearrangement algorithm          -a [VALUE]          Iterations of NNI + SPR (0)\n"
		   "                                                          TBR followed by iterations of NNI + SPR (1)\n"
		   "                                                          Pseudo counts determining probability of NNI, SPR, or TBR iterations (2)\n");

	printf("\n");

	printf(" Simulated Annealing Preferences: \n");
	printf("    Cooling schedule                  -c [G|L]            SA cooling schedule; GEOMETRIC (G) or LINEAR (L); default: G \n");

	printf("\n");

	printf(" Parallelisation Preferences: \n");
	printf("    Fine-grain multithreading         -p [VALUE]          Number of threads requested; default: one thread \n");

	printf("\n");

	exit(0);
}

void read_parameters(Parameters *prms, int argc, char **argv)
{

	int c;
	opterr = 0;

	while ((c = getopt(argc, argv, "t:c:b:vs:i:o:f:a:p:N:SC:h?")) != -1)
	{
		switch (c)
		{
		case 'c': /* cooling schedule */
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -c [g|l]\n", optopt);
				usage(argv[0]);
			}
			if (strcmp(optarg, "g") == 0 || strcmp(optarg, "G") == 0)
				prms->cooling_schedule = 0;
			else if (strcmp(optarg, "l") == 0 || strcmp(optarg, "L") == 0)
				prms->cooling_schedule = 1;
			else
			{
				fprintf(stderr, "Unknown cooling schedule option\nPlease, choose between Geometric (g) or Linear (l).");
			}
			break;
		case 'a': /* algorithm selection */
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%d requires an argument -a [0|1|2]\n", optopt);
				usage(argv[0]);
			}
			if (strcmp(optarg, "0") == 0 || strcmp(optarg, "0") == 0)
				prms->algorithm_selection = 0;
			else if (strcmp(optarg, "1") == 0 || strcmp(optarg, "1") == 0)
				prms->algorithm_selection = 1;
			else if (strcmp(optarg, "2") == 0 || strcmp(optarg, "2") == 0)
				prms->algorithm_selection = 2;
			else
			{
				fprintf(stderr, "Unknown algorithm option\nPlease, choose between SN (0), SEQ-TNS (1), or PBS (2).");
			}
			break;
		case 'v': /* verbose */
			prms->verbose = LVB_TRUE;
			break;
		case 's': /* seed */
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -s <int>\n", optopt);
				usage(argv[0]);
			}
			prms->seed = atoi(optarg);
			break;
		case 'i': /* file name in */
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -i <file name>\n", optopt);
				usage(argv[0]);
			}
			if (strlen(optarg) > LVB_FNAMSIZE)
			{
				fprintf(stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
			}
			strcpy(prms->file_name_in, optarg);
			break;
		case 'o': /* file name out */
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -o <file name>\n", optopt);
				usage(argv[0]);
			}
			if (strlen(optarg) > LVB_FNAMSIZE)
			{
				fprintf(stderr, "Error, the length file name greater than %d\n", LVB_FNAMSIZE);
			}
			strcpy(prms->file_name_out, optarg);
			break;
		case 'f': /* format */
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -f [phylip|fasta|nexus|msf|clustal]\n", optopt);
				usage(argv[0]);
			}
			if (strcmp(optarg, "phylip") == 0)
			{
				prms->n_file_format = FORMAT_PHYLIP;
			}
			else if (strcmp(optarg, "fasta") == 0)
			{
				prms->n_file_format = FORMAT_FASTA;
			}
			else if (strcmp(optarg, "nexus") == 0)
			{
				prms->n_file_format = FORMAT_NEXUS;
			}
			else if (strcmp(optarg, "clustal") == 0)
			{
				prms->n_file_format = FORMAT_CLUSTAL;
			}
			else
			{
				fprintf(stderr, "Unknown file format.");
				print_formats_available();
			}
			break;
		case 'p':
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -p <int>\n", optopt);
				usage(argv[0]);
			}
			prms->n_processors_available = atoi(optarg);
			if (prms->n_processors_available < 1)
				prms->n_processors_available = 1;
			break;
		case 't':
			if (optarg == NULL)
			{
				fprintf(stderr, "Option -%c requires an argument -p <file name>\n", optopt);
				usage(argv[0]);
				exit(1);
			}
			prms->n_number_max_trees = atoi(optarg);
			if (prms->n_number_max_trees < 1)
				prms->n_number_max_trees = 0;
			break;
		case '?':
		case 'h':
		default:
			usage(argv[0]);
		}
	}
}
