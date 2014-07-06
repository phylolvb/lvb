/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

randpint.c - random positive integer operations

=cut

**********/

#include "lvb.h"

/**********

=head1 randpint - GET RANDOM POSITIVE INTEGER

=head2 SYNOPSIS

    long randpint(long upper);

=head2 DESCRIPTION

Returns a positive integer with a user-specified upper limit.

Before calling C<randpint()>, ensure the random number generator has
been initialized by calling C<rinit()>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item upper

The random number will be in the range [0..C<upper>] inclusive.

=back

=head3 RETURN

A random integer in the interval [0..C<upper>].

=cut

**********/

long randpint(const long upper)
{
    double frand;	/* random real */
    double fupper;	/* upper limit */
    long rand;		/* return value */

    lvb_assert(upper >= 0);

    fupper = (double) upper;
    frand = uni();
    frand = frand * fupper;		/* scale to right range */
    rand = (long) (frand + 0.5);	/* round to nearest integer */

    /* guard against arithmetic inaccuracy */
    if (rand < 0)
	rand = 0;
    else if (rand > upper)
	rand = upper;

    return rand;

} /* end randpint() */
