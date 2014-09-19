/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** getparam.c - get and set configurable parameters ********** */

#include "lvb.h"
#include <unistd.h>

static int get_default_seed(void)
/* return a default integer in the interval [0..MAX_SEED], obtained from the
 * system clock, or exit with an error message if the system time is
 * unavailable */
{
    time_t tim;			/* system time */
    unsigned long ul_seed;	/* seed value obtained from system time */

    tim = time(NULL);
    lvb_assert(tim != -1);
    ul_seed = (unsigned long) tim;
    ul_seed = ul_seed % (1UL + (unsigned long) MAX_SEED);
    lvb_assert(ul_seed <= MAX_SEED);
    return (int) ul_seed;

} /* end get_default_seed() */


void defaults_params(Params *const prms)
/* set seed in *prms to unacceptable value, and other parameters to their
 * defaults_params from lvb.h */
{
    prms->bootstraps = 0;	/* sensible default */

    /* meaningful value that is not user-configurable */
    prms->verbose = LVB_FALSE;

    /* dash Unkown */
    prms->fifthstate = LVB_FALSE;

    /* cooling schecdule Generic */
    prms->cooling_schedule = 0;
    /* default value that will usually be used */
    prms->seed = get_default_seed();

    strcpy(prms->file_name_in, "infile");
    strcpy(prms->file_name_out, "outfile");
    prms->n_file_format = FORMAT_PHYLIP;
    prms->n_processors_available = 1;

} /* end defaults_params() */


void print_formats_available(){

	printf("\nFormats available to read: phylip, fasta, nexus, msf and clustal");

}

// "dcbvsi:o:f:p"
void usage(char *p_file_name){
	printf("Usage: lvb [dcbvsiofp]\n");
	printf("lvb seeks parsimonious trees from an aligned nucleotide data matrix.\n"
			"It uses heuristic searches consisting of simulated annealing followed by hill-climbing.\n\n");
	printf("\n    -d [u|f] (u) treatment of gaps represented by '-' in the data matrix.\n"
			"        Unknown (u): ‘-’ is treated as equivalent to ‘?’, i.e., as an ambiguous site that may contain ‘A’ or ‘C’ or ‘G’ or ‘T’ or ‘O’.\n"
			"        Fifth state (f): ‘-’ is treated as equivalent to a gap.\n"
			"        default (u)");

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

	while ((c = getopt (argc, argv, "d:c:b:vs:i:o:f:p:")) != -1){
		switch (c)
		{
			case 'd':
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
					usage(argv[0]);
					abort();
				}
				if (strcmp(optarg, "u") == 0 || strcmp(optarg, "U") == 0) prms->fifthstate = LVB_FALSE; // dash in alignments is UNKNOWN
				else if (strcmp(optarg, "f") == 0 || strcmp(optarg, "F") == 0) prms->fifthstate = LVB_TRUE; // dash in alignments is UNKNOWN
				else{
					fprintf (stderr, "Option -%c requires an argument -d [u|f].\n", optopt);
					print_formats_available();
					abort();
				}
				break;
			case 'c':	// cooling schedule
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
			case 'b':	// bootstrap
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -b <int>\n", optopt);
					usage(argv[0]);
					abort();
				}
				prms->bootstraps = atoi(optarg);
				break;
			case 'v':	// verbose
				prms->verbose = LVB_TRUE;
				break;
			case 's':	// seed
				if (optarg == NULL){
					fprintf (stderr, "Option -%c requires an argument -s <int>\n", optopt);
					usage(argv[0]);
					abort();
				}
				prms->seed = atoi(optarg);
				break;
			case 'i':	// file name in
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
			case 'o':	// file name out
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
			case 'f':	// format
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


void getparam(Params *prms, int argc, char **argv)
/* Get configuration parameters. This function fills *prms with
 * run-time configuration parameters */
{
	defaults_params(prms);
 /*   user_adjust(prms);*/
	read_parameters(prms, argc, argv);

} /* end getparam() */
