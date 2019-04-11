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

/* ========== admin.c - LVB library data and administration ========== */

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
    if ((((long) INT_MAX) < 2147483647L) || ((sizeof(void *) * CHAR_BIT) < 32) || ((sizeof(size_t) * CHAR_BIT) < 32)) {
        crash("program requires at least a 32-bit system");
    }

    /* LVB_EPS is assumed to be bigger than DBL_EPSILON in code that guards
     * against floating-point arithmetic problems */
    if (DBL_EPSILON >= LVB_EPS)
        crash("program requires greater floating point precision");

    #ifdef MPI_Implementation
    if (!((LVB_EPS + INITIAL_INCREMENT) != INITIAL_INCREMENT))
        crash("LVB_EPS and INITIAL_INCREMENT are incompatible with floating point precision");
    #endif

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
