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

/* ********** getparam.c - get and set configurable parameters ********** */

#include "lvb.h"
#include <unistd.h>

/* it is in ReadFile.cpp library */
int read_parameters(Params *prms, int argc, char **argv);

int get_default_seed(void)
/* return a default integer in the interval [0..MAX_SEED], obtained from the
 * system clock, or exit with an error message if the system time is
 * unavailable */
{
	srand(time(NULL));
    return (int) (rand() % (unsigned long) MAX_SEED);
} /* end get_default_seed() */


void defaults_params(Params *const prms)
/* set seed in *prms to unacceptable value, and other parameters to their
 * defaults_params from lvb.h */
{
    prms->bootstraps = 0;	/* sensible default */

    /* meaningful value that is not user-configurable */
    prms->verbose = LVB_FALSE;

    /* cooling schecdule Generic */
    prms->cooling_schedule = 0;
    /* default value that will usually be used */
    prms->seed = get_default_seed();

    strcpy(prms->file_name_in, "infile");
    strcpy(prms->file_name_out, OUTTREEFNAM);
    prms->n_file_format = FORMAT_PHYLIP;
    prms->n_processors_available = 1;

} /* end defaults_params() */

int getparam(Params *prms, int argc, char **argv)
/* Get configuration parameters. This function fills *prms with
 * run-time configuration parameters */
{
	defaults_params(prms);

	int n_default_seed = prms->seed;
 /*   user_adjust(prms);*/
	int n_error_code = read_parameters(prms, argc, argv);

	/* change initial seed because the user defined one */
	if (prms->seed != n_default_seed){
		srand(prms->seed);
	}
	return n_error_code;
} /* end getparam() */
