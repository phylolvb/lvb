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

#include <Lvb.h>

/* Test for randpint(). Simply checks it works and doesn't always give
 * the same answer. (There's perhaps a miniscule chance it would give
 * the same answer, but if so, it really is miniscule.) The random
 * number generator is seeded from the clock so it's usually going to
 * be a different test every run. */

#define LOOP_CNT 5000000	/* iterations of main test loop */
#define UPPER_LIM 20000		/* arbitrary moderately large number */

int main(void)
{
    long i;				/* loop counter */
    long rand_val;			/* random number */
    long first_rand_val;		/* first random number */
    time_t tim;				/* system time */
    unsigned long ul_seed;		/* seed, from system time */
    Lvb_bool all_same = LVB_TRUE;	/* all 'random' values same */

    lvb_initialize();

    /* seed random number generator from system clock */
    tim = time(NULL);
    lvb_assert(tim != -1);
    ul_seed = (unsigned long) tim;
    ul_seed = ul_seed % (1UL + (unsigned long) MAX_SEED);
    lvb_assert(ul_seed <= MAX_SEED);
    rinit((int) ul_seed);
 
    first_rand_val = randpint(UPPER_LIM);
    for (i = 0; i < LOOP_CNT; i++)
    {
        rand_val = randpint(UPPER_LIM);
	lvb_assert(rand_val <= UPPER_LIM);
	lvb_assert(rand_val >= 0);
	if (rand_val != first_rand_val)
	    all_same = LVB_FALSE;
    }

    if (all_same == LVB_FALSE)
    {
        printf("test passed\n");
	return EXIT_SUCCESS;
    }
    else
    {
        printf("test failed\n");
	return EXIT_FAILURE;
    }
}
