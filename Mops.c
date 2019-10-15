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

/* ========== mops.c - basic memory operation ========== */

#include "Lvb.h"

void *alloc(const size_t bytes, const char *const msg)
{
    void *p;	/* pointer to first byte of new memory, if any */

    if ((bytes == 0) || (bytes > MAX_ALLOC))
        p = NULL;
    else
    {
        p = malloc(bytes);
        if (p == NULL)
        crash("out of memory: cannot allocate for %s", msg);
    }
    return p;

}	/* end alloc() */

void alloc_memory_to_getplen(Dataptr matrix, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs) {
	*p_todo_arr 		    = (long *) alloc((matrix->nbranches - matrix->n) * sizeof(long), "alloc to count runs");
    *p_todo_arr_sum_changes = (long *) alloc(matrix->n_threads_getplen * (1 + matrix->nbranches - matrix->n) * sizeof(long), "alloc to count runs");
	*p_runs 	        	= (int *) alloc(matrix->n_threads_getplen * (matrix->nbranches - matrix->n) * sizeof(int), "alloc to count runs");
}

void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs) {
	free(*p_todo_arr);
	free(*p_todo_arr_sum_changes);
	free(*p_runs);
}