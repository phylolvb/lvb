/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

/* test for crash() */

int main(void)
{
    long i = 32447881;

    lvb_initialize();

    crash("%ld is causing trouble", i);

    /* if we reach here, crash() isn't working */
    printf("test failed\n");
    return 0;	/* program success indicates crash()'s failure */
}
