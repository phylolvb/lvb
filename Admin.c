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

#include "Lvb.h"



static void functionality_check(void)
/* To the extent possible, check that standard functions and data types match
 * LVB's expectations. Crash verbosely if they are found not to. */
{
    if (time(NULL) == -1) crash("cannot get system time");

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
    #endif // MPI_Implementation

} // end functionality_check()

// checks that some features of the system are suitable for use with LVB

void lvb_initialize(void)
{
    functionality_check();

} // end lvb_initialize()
