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

/* ========== getparam.c - get and set configurable parameters ========== */

#include "lvb.h"
#include <unistd.h>

/* it is in ReadFile.cpp library */
//#ifdef NP_Implementation
void read_parameters(Params *prms, int argc, char **argv);
//#endif

#ifdef MPI_Implementation
#ifdef MAP_Reduce
	#include "LVB_READ_FILES/src/ReadFile.h"
#else
	int read_parameters(Params *prms, int argc, char **argv);
#endif

long get_random_maxaccept(void)
/* return a random integer value, for use as the maxaccept parameter in
 * the simulated annealing search */
{
    return randpint(MAXACCEPT_MAX - MAXACCEPT_MIN) + MAXACCEPT_MIN;
}

int get_default_seed(void)
#endif

//#ifdef NP_Implementation
static int get_default_seed(void)
//#endif

/* return a default integer in the interval [0..MAX_SEED], obtained from the
 * system clock, or exit with an error message if the system time is
 * unavailable */
{
    time_t tim;			/* system time */
    unsigned long ul_seed;	/* seed value obtained from system time */

    tim = time(NULL);
    lvb_assert(tim != -1);
    //#ifdef NP_Implementation
    ul_seed = (unsigned long) tim;
    ul_seed = ul_seed % (1UL + (unsigned long) MAX_SEED);
    //#endif
    #ifdef MPI_Implementation
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
    //#ifdef NP_Implementation
    prms->bootstraps = 0;	/* sensible default */
    //#endif

    #ifdef MPI_Implementation
    #ifndef MAP_Reduce
    prms->n_seeds_need_to_try = 1;
    prms->n_checkpoint_interval = CHECKPOINT_INTERVAL;
    #endif

    /* by default dont read and save states */
    prms->n_flag_save_read_states = DONT_SAVE_READ_STATES;
    prms->n_flag_is_finished_process = CHECK_POINT_PROCESS_NOT_FINISHED;
    prms->n_flag_is_possible_read_state_files = CHECK_POINT_NOT_READ_STATE_FILES;
    /* meaningful value that is not user-configurable */
    #endif

    prms->verbose = LVB_FALSE;

    /* cooling schecdule Generic */
    prms->cooling_schedule = 0;
    /* default value that will usually be used */
    prms->seed = get_default_seed();
    //#ifdef NP_Implementation
    /* original branch-swapping algorithm */
    prms->algorithm_selection = 0;
    //#endif

    strcpy(prms->file_name_in, "infile");
    strcpy(prms->file_name_out, OUTTREEFNAM);
    prms->n_file_format = FORMAT_PHYLIP;
    prms->n_processors_available = 1;
    //#ifdef NP_Implementation
    prms->n_number_max_trees = 0;			/* by default keepall of them */
    //#endif

    #ifdef MPI_Implementation
    /* it is olnly used for tests */
    prms->n_make_test = 0;	/* by default is not making any test */
    #endif

} /* end defaults_params() */

//#ifdef NP_Implementation
void getparam(Params *prms, int argc, char **argv)
//#endif

#ifdef MPI_Implementation
void getparam(Params *prms, int argc, char **argv)
#endif

/* Get configuration parameters. This function fills *prms with
 * run-time configuration parameters */
{
	defaults_params(prms);
 /*   user_adjust(prms);*/
    //#ifdef NP_Implementation
	read_parameters(prms, argc, argv);
    //#endif

    #ifdef MPI_Implementation
    #ifndef MAP_Reduce
        int n_default_seed = prms->seed;
    #endif
        int n_error_code = read_parameters(prms, argc, argv);

    #ifndef MAP_Reduce
        /* change initial seed because the user defined one */
        if (prms->seed != n_default_seed){
	        srand(prms->seed);
        }
    #endif
        return n_error_code;
    #endif

} /* end getparam() */
