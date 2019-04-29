/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== randpint.c - random positive integer operations ========== */

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
#ifdef NP_Implementation
int randpint(const long upper)
#endif

#ifdef MPI_Implementation
long randpint(const long upper)
#endif
{
    double frand;	/* random real */
    double fupper;	/* upper limit */
    #ifdef NP_Implementation
    int rand;		/* return value */
    #endif

    #ifdef MPI_Implementation
    long rand;		/* return value */
    #endif
    lvb_assert(upper >= 0);

    fupper = (double) upper;
    frand = uni();
    frand = frand * fupper;		/* scale to right range */
    #ifdef NP_Implementation
    rand = (int) (frand + 0.5);	/* round to nearest integer */
    #endif

    #ifdef MPI_Implementation
    rand = (long) (frand + 0.5);	/* round to nearest integer */
    #endif
    /* guard against arithmetic inaccuracy */
    if (rand < 0)
    rand = 0;
    else if (rand > upper)
    rand = upper;
    
    return rand;

} /* end randpint() */
