/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

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
