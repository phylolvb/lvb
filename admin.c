/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

admin.c - LVB library data and administration

=cut

**********/

#include "lvb.h"

/**********

=head1 matrix - DATA MATRIX

=head2 SYNOPSIS

Dataptr matrix;

=head2 DESCRIPTION

Pointer to the data matrix to use with the library during the current
program's execution. It is not possible to change the matrix used once
it has been read from file.

To discourage direct access, C<matrix> is not declared in C<lvb.h>. It
should be declared in functions as required.

=cut

**********/

Dataptr matrix = NULL;

static void functionality_check(void)
/* To the extent possible, check that standard functions and data types match
 * LVB's expectations. Crash verbosely if they are found not to. */
{
    /* time() is expected to work without error for logging the start
     * and end time, and for generating the default random number seed */
    if (time(NULL) == -1)
	crash("cannot get system time");

    /* if the system is not 32-bit, 64-bit, or more, some limits will be
     * less than documented and there may be memory allocation constraints
     * that LVB does not allow for */
    if ((((long) INT_MAX) < 2147483647L)
     || ((sizeof(void *) * CHAR_BIT) < 32)
     || ((sizeof(size_t) * CHAR_BIT) < 32))
    {
        crash("program requires at least a 32-bit system");
    }

    /* LVB_EPS is assumed to be bigger than DBL_EPSILON in code that guards
     * against floating-point arithmetic problems */
    if (DBL_EPSILON >= LVB_EPS)
        crash("program requires greater floating point precision");

    /* DBL_MANT_DIG is checked in rinit() so check not necessary here */

} /* end functionality_check() */

/**********

=head1 lvb_initialize - INITIALIZE lvb LIBRARY

=head2 SYNOPSIS

void lvb_initialize(void);

=head2 DESCRIPTION

Initializes the LVB library. Must be called once, before any other LVB
functions.

Currently, this function just checks that some features of the system
are suitable for use with LVB. If not, it crashes verbosely.

=cut

**********/

void lvb_initialize(void)
{
    functionality_check();

} /* end lvb_initialize() */
