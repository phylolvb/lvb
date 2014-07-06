/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

/* check pow_wrapper crashes verbosely on domain error if first
 * parameter is negative and second parameter is non-integral */

int main(void)
{
    lvb_initialize();

    pow_wrapper(-3, 1.5);	/* should cause domain error */

    /* if we get this far, pow_wrapper() has failed to crash so the
     * test has failed */
    printf("test failed\n");

    exit(EXIT_SUCCESS);	/* we want failure: so program success
                         * indicates test failed */
}
