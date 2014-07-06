/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include "lvb.h"

/* An array of TEST_ALLOC_LONGS long ints is allocated and the
 * test expects this allocation to succeed. */
#define TEST_ALLOC_LONGS 1000000

/* Subsequently, that array is freed and an array of TEST_ALLOC_CHARS
 * characters is allocated, and the test expects this to succeed as
 * well. */
#define TEST_ALLOC_CHARS 70000

long *lp = NULL;	/* pointer to start of an allocated array */
char *cp = NULL;	/* pointer to start of an allocated array */

int main(void)
{
    long i;	/* loop counter */

    lvb_initialize();

    /* Check that a zero byte allocation returns a NULL pointer. */
    lp = alloc(0, "zero-byte array");
    lvb_assert(lp == NULL);

    /* Check that a fairly large allocation succeeds and that the
     * memory allocated may be written to. */
    lp = alloc(TEST_ALLOC_LONGS * sizeof(long), "test array 1");
    lvb_assert(lp != NULL);
    for (i = 0; i < TEST_ALLOC_LONGS; i++)
	lp[i] = 1;

    free(lp);

    /* Basic check that heap is OK after free: allocate something else
     * and check we may write to it */
    cp = alloc(TEST_ALLOC_CHARS, "test array 2");
    for (i = 0; i < TEST_ALLOC_CHARS; i++)
    	cp[i] = 'X';

    free(cp);
    
    printf("test passed\n");
    return EXIT_SUCCESS;
}
