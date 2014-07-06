/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

/* check pow_wrapper crashes verbosely on range error if result is too
 * large to represent */

int main(void)
{
    lvb_initialize();

    pow_wrapper(DBL_MAX, 2);	/* should cause range error */

    /* if we get this far, pow_wrapper() has failed to crash so the
     * test has failed */
    printf("test failed\n");

    exit(EXIT_SUCCESS);	/* we want failure: so program success
                         * indicates test failed */
}
