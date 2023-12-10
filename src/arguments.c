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

/* ========== arguments.c - get and set configurable arguments ========== */

#include "arguments.h"

/* it is in CommandLineParser.cpp library */
void ParseArguments(Arguments *args, int argc, char **argv);

static int SetDefaultSeed(void)
/* return a default integer in the interval [0..MAX_SEED], obtained from the
 * system clock, or exit with an error message if the system time is
 * unavailable */
{
	time_t tim;			   /* system time */
	unsigned long ul_seed; /* seed value obtained from system time */

	tim = time(NULL);
	lvb_assert(tim != -1);
	ul_seed = (unsigned long)tim;
	ul_seed = ul_seed % (1UL + (unsigned long)MAX_SEED);
	lvb_assert(ul_seed <= MAX_SEED);
	return (int)ul_seed;

} /* end SetDefaultSeed() */

void SetDefaultArgumentsValues(Arguments *const args)
/* set seed in *args to unacceptable value, and other arguments to their
 * SetDefaultArgumentsValues from LVB.h */
{
	/* meaningful value that is not user-configurable */
	args->verbose = LVB_FALSE;

	/* cooling schecdule Generic */
	args->cooling_schedule = 0;
	/* default value that will usually be used */
	args->seed = SetDefaultSeed();
	/* original branch-swapping algorithm */
	args->algorithm_selection = 1;

	strcpy(args->file_name_in, "infile");
	strcpy(args->file_name_out, OUTTREEFNAM);
	args->n_file_format = FORMAT_PHYLIP;
	args->num_threads = omp_get_max_threads();
	args->n_number_max_trees = 0; /* default, keep all EPT */

} /* end SetDefaultArgumentsValues() */

void GetArguments(Arguments *args, int argc, char **argv)
/* Get configuration arguments. This function fills *args with
 * run-time configuration arguments */
{
	SetDefaultArgumentsValues(args);
	ParseArguments(args, argc, argv);

} /* end GetArguments() */
