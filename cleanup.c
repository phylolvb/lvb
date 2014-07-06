/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** cleanup.c - prepare for exit ********** */

#include "lvb.h"

Lvb_bool cleanup(void)
/* prevent apparent memory leaks to help debugging, log end time; return
 * LVB_TRUE on write error to stdout, LVB_FALSE otherwise */
{
    time_t endtim;	/* time at end of run */
    Lvb_bool val;	/* return value */

    endtim = time(NULL);
    printf("\n");
    printf("Ending at: %s", ctime(&endtim));
    printf("\n");

    /* log file won't be used again */
    fflush(stdout);
    if (ferror(stdout) != 0)
	val = LVB_TRUE;
    else
        val = LVB_FALSE;

    return val;
} /* end cleanup() */
