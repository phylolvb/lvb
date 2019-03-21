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

/* ========== cleanup.c - prepare for exit ========== */

#include "lvb.h"

#ifdef NP_Implementation

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
    if (ferror(stdout) != 0) val = LVB_TRUE;
    else val = LVB_FALSE;
    return val;
} /* end cleanup() */

#endif // #ifdef NP_Implementation //

#ifdef MPI_Implementation

Lvb_bool cleanup(void)
/* prevent apparent memory leaks to help debugging, log end time; return
 * LVB_TRUE on write error to stdout, LVB_FALSE otherwise */
{
    time_t endtim;	/* time at end of run */
    Lvb_bool val = LVB_TRUE;	/* return value */

    endtim = time(NULL);
    printf("\n");
    printf("Ending at: %s", ctime(&endtim));
    printf("\n");

    /* log file won't be used again */
#ifdef MAP_REDUCE_SINGLE
    fflush(stdout);
    if (ferror(stdout) != 0) val = LVB_TRUE;
    else val = LVB_FALSE;
#endif

    return val;
} /* end cleanup() */

#endif

