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

/* ========== search_parameters.c - get and set configurable parameters ========== */

#include "lvb.h"

#include <unistd.h>

#ifdef MAP_REDUCE_SINGLE
	#include "input_options.h"
#else
	int read_parameters(Params *prms, int argc, char **argv);
#endif

#ifndef NP_Implementation
long get_random_maxaccept(void)
/* return a random integer value, for use as the maxaccept parameter in
 * the simulated annealing search */
{
    return randpint(MAXACCEPT_MAX - MAXACCEPT_MIN) + MAXACCEPT_MIN;
}
#endif

int get_default_seed(void)
/* return a default integer in the interval [0..MAX_SEED], obtained from the
 * system clock, or exit with an error message if the system time is
 * unavailable */
{
    time_t tim;		        /* system time */
    unsigned long ul_seed;	/* seed */

    tim = time(NULL);
    lvb_assert(tim != -1);
    #ifndef NP_Implementation
    ul_seed = (unsigned long) tim;
    ul_seed = ul_seed % (1UL + (unsigned long) MAX_SEED);
    #else
    srand(tim);
    ul_seed = (unsigned long) rand() % (1UL + (unsigned long) MAX_SEED);
    #endif
    lvb_assert(ul_seed <= MAX_SEED);
    return (int) ul_seed;

} /* end get_default_seed() */

void defaults_params(Params *const prms)
/* set seed in *prms to unacceptable value, and other parameters to their
 * defaults_params from lvb.h */
{
#ifdef MPI_Implementation 
    prms->n_seeds_need_to_try = 1;
    prms->n_checkpoint_interval = CHECKPOINT_INTERVAL;
#endif

#ifndef NP_Implementation
    /* by default dont read and save states */
    prms->n_flag_save_read_states = DONT_SAVE_READ_STATES;
    prms->n_flag_is_finished_process = CHECK_POINT_PROCESS_NOT_FINISHED;
    prms->n_flag_is_possible_read_state_files = CHECK_POINT_NOT_READ_STATE_FILES;
#endif

    /* meaningful value that is not user-configurable */
    prms->verbose = LVB_FALSE;

    /* cooling schecdule Generic */
    prms->cooling_schedule = 0;
    /* default value that will usually be used */
    prms->seed = get_default_seed();
    // original branch-swapping algorithm
    prms->algorithm_selection = 0;
    // default log interval
    prms->STAT_LOG_INTERVAL = 50000;

    strcpy(prms->file_name_in, "infile");
    strcpy(prms->file_name_out, OUTTREEFNAM);
    prms->n_file_format = FORMAT_PHYLIP;
    prms->n_processors_available = 1;
} /* end defaults_params() */

int Search_Parameters(Params *prms, int argc, char **argv)
/* Get configuration parameters. This function fills *prms with
   run-time configuration parameters */
{
        defaults_params(prms);

#ifndef NP_Implementation

#ifdef MPI_Implementation
        int n_default_seed = prms->seed;
#endif
        int n_error_code = read_parameters(prms, argc, argv);

#ifdef MPI_Implementation
        /* change initial seed because the user defined one */
        if (prms->seed != n_default_seed){
	        srand(prms->seed);
        }
#endif
        return n_error_code;
#else
    read_parameters(prms, argc, argv);
#endif
} /* end Search_Parameters() */