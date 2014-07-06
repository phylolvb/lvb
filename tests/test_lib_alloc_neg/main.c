/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* Allocate the largest array allowed by alloc(). If allocation fails,
 * print the string "FATAL ERROR: out of memory: cannot allocate for
 * massive array\n" to standard output. If allocation succeeds, check
 * we can write to the whole array, then print "test passed\n". Any
 * output other than one of these two strings indicates failure. */

#include <lvb.h>

#define ALLOCATION MAX_ALLOC
#define SMALL_ALLOCATION 1000000	/* doubles */

char *p = NULL;			/* pointer to large array */
double *small_p = NULL;		/* pointer to small array */

int main(void)
{
    size_t i;	/* loop counter */

    lvb_initialize();
    p = alloc(ALLOCATION, "massive array");

    /* It is likely we do not reach this stage. If we have, check that
     * the allocation really did succeed, by writing to the whole
     * array, freeing it, and testing allocation of a smaller array. */
    for (i = 0; i < ALLOCATION; i++)
    	p[i] = ' ';
    free(p);
    small_p = alloc(SMALL_ALLOCATION * sizeof(double), "second array");
    for (i = 0; i < SMALL_ALLOCATION; i++)
        small_p[i] = DBL_MAX;
    free(small_p);

    printf("test passed\n");
    return EXIT_SUCCESS;
}
