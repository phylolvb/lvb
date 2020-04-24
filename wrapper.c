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

/* ========== wrapper.c - LVB to PHYLIP interface ========== */

#include "lvb.h"

#ifdef LVB_MAPREDUCE

#include "ReadFile.h"

#else

void ReadFile(char *file_name, int n_file_format, Dataptr p_lvbmat);
void ReadDNAMatrixExternal(char *file_name, int n_file_format, long *species_ptr, long *sites_ptr, int *max_length_name);

#endif

void CheckDNAMatrixInput(char *p_file_name, int n_file_format, Dataptr lvbmat)
{
	ReadFile(p_file_name, n_file_format, lvbmat);

    /* check number of sequences is in range for LVB */
    if (lvbmat->n < MIN_N) CrashVerbosely("The data matrix must have at least %ld sequences.", MIN_N);
    else if (lvbmat->n > MAX_N) CrashVerbosely("The data matrix must have no more than %ld sequences.", MAX_N);
    /* check number of sites is in range for LVB */
    else if (lvbmat->m < MIN_M) CrashVerbosely("The data matrix must have at least %ld sites.", MIN_M);
    else if (lvbmat->m > MAX_M) CrashVerbosely("The data matrix must have no more than %ld sites.", MAX_M);

    /* maximum number of object sets per tree */
    lvb_assert (lvbmat->nsets <= (MAX_N - 3));

} /* end CheckDNAMatrixInput() */

void ReadDNAMatrix(char *p_file_name, int n_file_format, long *species_ptr, long *sites_ptr, int *max_length_name){

	ReadDNAMatrixExternal(p_file_name, n_file_format, species_ptr, sites_ptr, max_length_name);
}
