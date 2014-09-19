/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

mops.c - basic memory operation

=cut

*********/

#include "lvb.h"

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
exceeds C<MAX_ALLOC> bytes, where C<MAX_ALLOC> is defined in C<lvb.h>.
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

void alloc_memory_to_getplen(Dataptr matrix, long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs){
	*p_todo_arr = (long *) alloc(matrix->n_threads_try_combination * (matrix->nbranches - matrix->n) * sizeof(long), "alloc to count runs");
    *p_todo_arr_sum_changes = (long *) alloc(matrix->n_threads_try_combination * matrix->n_threads_getplen * (1 + matrix->nbranches - matrix->n) * sizeof(long), "alloc to count runs");
	*p_runs = (int *) alloc(matrix->n_threads_try_combination * matrix->n_threads_getplen * (matrix->nbranches - matrix->n) * sizeof(int), "alloc to count runs");
}

void free_memory_to_getplen(long **p_todo_arr, long **p_todo_arr_sum_changes, int **p_runs){
	free(*p_todo_arr);
	free(*p_todo_arr_sum_changes);
	free(*p_runs);
}
