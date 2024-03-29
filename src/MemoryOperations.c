/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== MemoryOperations.c - basic memory operation ========== */

#include "MemoryOperations.h"

/**********

=head1 alloc - ALLOCATE DYNAMIC HEAP MEMORY

=head2 SYNOPSIS

void *alloc(const size_t bytes, const char *const msg);

=head2 DESCRIPTION

B<alloc> is a wrapper for the standard library function C<malloc>. It
allocates new memory, or, on failure, crashes with an error message. It
also assures a system-independent return value, C<NULL>, if a zero-byte
allocation is requested.

B<alloc> will always fail without calling C<malloc> if the request
exceeds C<MAX_ALLOC> bytes, where C<MAX_ALLOC> is defined in C<LVB.h>.
This avoids straining C<malloc> beyond typical use. C<MAX_ALLOC> is a
large number that will comfortably fit in a 32-bit signed integer.

The allocated memory may be freed with the standard library function
C<free>.

The new memory is not initialized.

=head2 PARAMETERS

=head3 INPUT

=over4

=item bytes

The number of bytes to allocate. If C<bytes> is zero, no memory is
allocated.

=item msg

Pointer to the first text character in a string that describes the
object being allocated. On allocation failure, B<alloc> will crash with
an error message of 'FATAL ERROR: out of memory: cannot allocate for ',
followed by this string.

=back

=head2 RETURN

Returns pointer to the first byte of the newly allocated memory, which
is suitably aligned for storage of any object, or C<NULL> if a
zero-byte allocation was requested.

=cut

**********/

void *alloc(const size_t bytes, const char *const msg)
{
	void *p; /* pointer to first byte of new memory, if any */

	if ((bytes == 0) || (bytes > MAX_ALLOC))
		p = NULL;
	else
	{
		p = malloc(bytes);
		if (p == NULL)
			crash("out of memory: cannot allocate for %s", msg);
	}
	return p;

} /* end alloc() */

void alloc_memory_to_getplen(Dataptr MSA, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs)
{
	*p_todo_arr = (long *)alloc((MSA->numberofpossiblebranches - MSA->n) * sizeof(long), "alloc to count runs");
	*p_todo_arr_sum_changes = (long *)alloc(MSA->n_threads_getplen * (1 + MSA->numberofpossiblebranches - MSA->n) * sizeof(long), "alloc to count runs");
	*p_runs = (int *)alloc(MSA->n_threads_getplen * (MSA->numberofpossiblebranches - MSA->n) * sizeof(int), "alloc to count runs");
}

void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs)
{
	free(*p_todo_arr);
	free(*p_todo_arr_sum_changes);
	free(*p_runs);
}

/* set the number of processors to use */
void calc_distribution_processors(Dataptr MSA, Parameters rcstruct)
{
	int n_threads_temp = 0;
	if (MSA->nwords > MINIMUM_SIZE_NUMBER_WORDS_TO_ACTIVATE_THREADING)
	{
		do
		{
			n_threads_temp++;
			MSA->n_slice_size_getplen = MSA->nwords / n_threads_temp;
		} while (MSA->n_slice_size_getplen > MINIMUM_WORDS_PER_SLICE_GETPLEN && n_threads_temp != rcstruct.n_processors_available);

		if (MSA->n_slice_size_getplen > MINIMUM_WORDS_PER_SLICE_GETPLEN)
		{
			MSA->n_slice_size_getplen = MSA->nwords / n_threads_temp;
			MSA->n_threads_getplen = n_threads_temp;
		}
		else
		{
			MSA->n_threads_getplen = n_threads_temp - 1;
			MSA->n_slice_size_getplen = MSA->nwords / MSA->n_threads_getplen;
		}
	}
	else
	{
		MSA->n_threads_getplen = 1; /* need to pass for 1 thread because the number of words is to low */
	}
	// only to protect
	if (MSA->n_threads_getplen < 1)
		MSA->n_threads_getplen = 1;
}
