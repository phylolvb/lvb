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

/* ========== SearchParameters.c - get and set configurable parameters ========== */

#include "SearchParameters.h"

/* it is in CommandLineParser.cpp library */
void read_parameters(Parameters *prms, int argc, char **argv);

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


void defaults_params(Parameters *const prms)
/* set seed in *prms to unacceptable value, and other parameters to their
 * defaults_params from LVB.h */
{
    /* meaningful value that is not user-configurable */
    prms->verbose = LVB_FALSE;

    /* cooling schecdule Generic */
    prms->cooling_schedule = 0;
    /* default value that will usually be used */
    prms->seed = get_default_seed();
    /* original branch-swapping algorithm */
    prms->algorithm_selection = 1;

    strcpy(prms->file_name_in, "infile");
    strcpy(prms->file_name_out, OUTTREEFNAM);
    prms->n_file_format = FORMAT_PHYLIP;
    prms->n_processors_available = omp_get_max_threads();
	prms->n_number_max_trees = 0;			/* default, keep all EPT */

} /* end defaults_params() */

void getparam(Parameters *prms, int argc, char **argv)
/* Get configuration parameters. This function fills *prms with
 * run-time configuration parameters */
{
	defaults_params(prms);
 /*   user_adjust(prms);*/
	read_parameters(prms, argc, argv);

} /* end getparam() */


#ifdef LVB_MPI
	void writeinf(Parameters prms, Dataptr MSA, int argc, char **argv, int clusterSize)
#else
	void writeinf(Parameters prms, Dataptr MSA, int argc, char **argv)
#endif

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
	LogTime();
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

	printf("  MSA size:            %ld x %ld\n", MSA->n, MSA->original_m);
	printf("  Seed:                %d\n", prms.seed);
	printf("  Cooling schedule:    ");
    if(prms.cooling_schedule == 0) printf("GEOMETRIC\n");
    else printf("LINEAR\n");
	printf("  Algorithm: ");
    if(prms.algorithm_selection == 0) printf("          0 (SN)\n");
    else if(prms.algorithm_selection == 1) printf("          1 (SEQ-TNS)\n");
    else if(prms.algorithm_selection == 2) printf("          2 (PBS)\n");

	printf("\nParallelisation Properties: \n");

	printf("  Additional MPI Seeds: %d: ", clusterSize - 1);

	for(int i = 1; i < clusterSize; i++) {
		printf("%d ", prms.seed + i);
	}

	printf("\n================================================================================\n");
	printf("\nInitialising search: \n");
}
